# ============================================================================
#  Bayesian Precision Factor Analysis for Time Series
#  Simulation C: Data Generated from Kastner & Huber (2020) DGP
#                VAR(1) with factor stochastic volatility on errors
#                (their model exactly — our model is misspecified here)
#
#  DGP follows Kastner & Huber (2020) Section 5:
#    - VAR(1), p=20
#    - Three sparsity scenarios for B: sparse (1%), intermediate (10%), dense (80%)
#    - Single-factor SV on errors: Sigma_t = lambda*lambda' * exp(h_t) + diag(exp(g_t))
#    - Factor loadings: lambda ~ N(0.001, 0.001^2)
#    - Idiosyncratic log-variances: mu=-12, persistence ~ U(0.85, 0.98), sigma=0.1-0.3
#    - Factor log-variance: highly persistent, rho=0.99, sigma=0.1
#
#  Evaluation: held-out prediction MSPE (no static true Omega exists)
#    - Train on first 80% of observations
#    - At each test time point, hold out 30% of variables randomly
#    - Predict held-out variables given observed using estimated Omega
#
#  Methods:
#    1. PrecFact (no VAR)
#    2. VAR + PrecFact (sequential)
#    3. Joint sampler (proposed)
#    4. Kastner (VAR + FSV) — oracle method for this DGP
#    5. Sample covariance baseline
#
#  5 replicates per scenario
# ============================================================================

library(Rcpp)
library(RcppArmadillo)
library(PrecFact)
library(bayesianVARs)
library(MASS)

Rcpp::sourceCpp("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/var_prec.cpp")

# ============================================================================
#  SECTION 1: DGP — Kastner & Huber (2020)
# ============================================================================

simulate_kastner_dgp <- function(T, p, sparsity=0.10, seed=NULL) {
   if (!is.null(seed)) set.seed(seed)
   
   # ── VAR coefficient matrix B ──────────────────────────────────────────────
   # intercept ~ N(0.01, 0.01^2), diagonal always nonzero
   # off-diagonal nonzero with probability = sparsity
   # sparse:       mu_D=0.15, sd_D=0.15, mu_O=0.01, sd_O=0.01 (typo-corrected)
   # intermediate: mu_D=0.15, sd_D=0.15, mu_O=0.1,  sd_O=0.1
   # dense:        mu_D=0.15, sd_D=0.15, mu_O=0.1,  sd_O=0.1
   if (sparsity <= 0.01) {
      mu_O <- 0.01; sd_O <- 0.01
   } else {
      mu_O <- 0.1;  sd_O <- 0.1
   }
   mu_D <- 0.15; sd_D <- 0.15
   
   B <- matrix(0, p, p)
   # diagonal always nonzero
   diag(B) <- rnorm(p, mu_D, sd_D)
   # off-diagonal with sparsity probability
   for (i in seq_len(p)) {
      for (j in seq_len(p)) {
         if (i != j && runif(1) < sparsity) {
            B[i,j] <- rnorm(1, mu_O, sd_O)
         }
      }
   }
   # scale to ensure stationarity
   sr <- max(abs(eigen(B, only.values=TRUE)$values))
   if (sr >= 1) B <- B * 0.95 / sr
   
   # ── Single-factor SV error structure ──────────────────────────────────────
   # factor loadings: lambda ~ N(0.001, 0.001^2)
   lambda_fac <- rnorm(p, 0.001, 0.001)
   
   # idiosyncratic log-variance parameters
   mu_sv    <- rep(-12, p)
   rho_sv   <- runif(p, 0.85, 0.98)
   sigma_sv <- runif(p, 0.1, 0.3)
   
   # factor log-variance parameters
   rho_fac   <- 0.99
   sigma_fac <- 0.1
   
   # simulate log-variances
   h_fac <- numeric(T)      # factor log-variance
   h_idi <- matrix(0, T, p) # idiosyncratic log-variances
   
   h_fac[1] <- rnorm(1, 0, sigma_fac / sqrt(1 - rho_fac^2))
   for (t in 2:T) h_fac[t] <- rho_fac * h_fac[t-1] + rnorm(1, 0, sigma_fac)
   
   for (j in seq_len(p)) {
      h_idi[1,j] <- rnorm(1, mu_sv[j], sigma_sv[j] / sqrt(1 - rho_sv[j]^2))
      for (t in 2:T)
         h_idi[t,j] <- mu_sv[j] + rho_sv[j]*(h_idi[t-1,j]-mu_sv[j]) +
            rnorm(1, 0, sigma_sv[j])
   }
   
   # simulate errors: epsilon_t ~ N(0, Sigma_t)
   # Sigma_t = lambda*lambda' * exp(h_fac_t) + diag(exp(h_idi_t))
   epsilon <- matrix(0, T, p)
   for (t in seq_len(T)) {
      Sigma_t <- tcrossprod(lambda_fac) * exp(h_fac[t]) +
         diag(exp(h_idi[t,]))
      epsilon[t,] <- MASS::mvrnorm(1, rep(0,p), Sigma_t)
   }
   
   # ── Simulate VAR(1) ───────────────────────────────────────────────────────
   Y      <- matrix(0, T, p)
   Y[1,]  <- epsilon[1,]
   for (t in 2:T) Y[t,] <- B %*% Y[t-1,] + epsilon[t,]
   
   # time-average Sigma for reference (not used in evaluation)
   Sigma_avg <- Reduce("+", lapply(seq_len(T), function(t)
      tcrossprod(lambda_fac) * exp(h_fac[t]) + diag(exp(h_idi[t,])))) / T
   Omega_avg <- solve(Sigma_avg)
   
   list(Y=Y, B=B, Omega_avg=Omega_avg, Sigma_avg=Sigma_avg,
        lambda_fac=lambda_fac, h_fac=h_fac, h_idi=h_idi)
}

# ============================================================================
#  SECTION 2: Helper functions
# ============================================================================

# conditional mean prediction using precision matrix
predict_held_out <- function(Omega, y_obs, obs_idx, hid_idx) {
   as.vector(-solve(Omega[hid_idx, hid_idx]) %*%
                Omega[hid_idx, obs_idx] %*% y_obs)
}

# PrecFact wrapper
fit_precfact <- function(Y_input, nrun=5000, thin=5) {
   res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=5,
                   nrun=nrun, thin=thin, Y=scale(Y_input),
                   alldiff=0, a_dir=0.1, b_dir=0.1)
   fact_prec_processing(res, burnout=(nrun/2)/thin, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# VAR residuals
get_var_residuals <- function(Y) {
   set.seed(99)
   p     <- ncol(Y)
   Y_jit <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-6), nrow(Y), p)
   mod   <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
   B_hat <- apply(mod$PHI, c(1,2), mean)
   list(resid=mod$Y - mod$X %*% B_hat)
}

# joint sampler wrapper
fit_joint <- function(Y, nrun=5000, thin=5) {
   set.seed(99)
   p     <- ncol(Y)
   Y_jit <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-6), nrow(Y), p)
   res   <- var_precision_sampler(Y=Y_jit, lags=1, k=5,
                                  nrun=nrun, burn=nrun/2, thin=thin)
   fact_prec_processing(res, burnout=(nrun/2)/thin, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# Kastner: VAR + factor SV
fit_kastner <- function(Y) {
   set.seed(99)
   p     <- ncol(Y)
   Y_jit <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-6), nrow(Y), p)
   sigma_prior <- specify_prior_sigma(
      M                       = p,
      type                    = "factor",
      factor_factors          = 1L,
      factor_priorfacloadtype = "colwiseng",
      factor_heteroskedastic  = TRUE,
      quiet                   = TRUE
   )
   mod       <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000,
                     prior_sigma=sigma_prior, quiet=TRUE)
   vc        <- vcov(mod)
   Sigma_hat <- apply(vc[1,,,], c(1,2), mean)
   list(Omega=solve(Sigma_hat))
}

# ── Run predictive evaluation on one dataset ──────────────────────────────────
run_prediction <- function(sim, T_train_frac=0.8, seed_pred=42) {
   Y       <- scale(sim$Y)
   T_obs   <- nrow(Y)
   p       <- ncol(Y)
   T_train <- floor(T_train_frac * T_obs)
   T_test  <- T_obs - T_train
   Y_train <- Y[1:T_train, ]
   Y_test  <- Y[(T_train+1):T_obs, ]
   n_hid   <- floor(0.3 * p)
   
   # fit all methods on training data
   cat("      PrecFact (no VAR)...\n")
   pf_raw  <- fit_precfact(Y_train)
   
   cat("      VAR + PrecFact...\n")
   var_out <- get_var_residuals(Y_train)
   pf_var  <- fit_precfact(var_out$resid)
   
   cat("      Joint sampler...\n")
   joint   <- fit_joint(Y_train)
   
   cat("      Kastner...\n")
   kast    <- fit_kastner(Y_train)
   
   Omega_samp <- solve(cov(Y_train))
   
   # predictive evaluation
   set.seed(seed_pred)
   mspe <- matrix(0, T_test, 5)
   colnames(mspe) <- c("PrecFact","VAR+PrecFact","Joint","Kastner","SampCov")
   
   for (t in seq_len(T_test)) {
      y_t     <- Y_test[t,]
      hid_idx <- sort(sample(p, n_hid))
      obs_idx <- setdiff(seq_len(p), hid_idx)
      y_obs   <- y_t[obs_idx]
      y_true  <- y_t[hid_idx]
      
      mspe[t,] <- c(
         mean((y_true - predict_held_out(pf_raw$Omega,  y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(pf_var$Omega,  y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(joint$Omega,   y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(kast$Omega,    y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(Omega_samp,    y_obs, obs_idx, hid_idx))^2)
      )
   }
   
   colMeans(mspe)
}

# summarize replicates
summarize_reps <- function(rep_list, scenario_label) {
   mat <- do.call(rbind, rep_list)
   data.frame(
      Scenario     = scenario_label,
      Method       = colnames(mat),
      MSPE         = sprintf("%.4f (%.4f)", colMeans(mat), apply(mat, 2, sd))
   )
}

# ============================================================================
#  SECTION 3: Run simulation C
#  Fixed: p=20, T=200. Varying: sparsity in {0.01, 0.10, 0.80}. 5 replicates.
# ============================================================================

cat("\n============================================================\n")
cat("SIMULATION C: Kastner & Huber DGP (factor SV errors)\n")
cat("Evaluation: held-out prediction MSPE\n")
cat("p=20, T=200, 5 replicates per sparsity scenario\n")
cat("============================================================\n")

sparsity_levels <- c(0.01, 0.10, 0.80)
sparsity_labels <- c("Sparse (1%)", "Intermediate (10%)", "Dense (80%)")
n_reps          <- 5
results_C       <- list()

for (s in seq_along(sparsity_levels)) {
   sp    <- sparsity_levels[s]
   label <- sparsity_labels[s]
   cat("\n  Scenario:", label, "\n")
   reps  <- vector("list", n_reps)
   
   for (r in seq_len(n_reps)) {
      cat("    replicate", r, "/", n_reps, "\n")
      sim      <- simulate_kastner_dgp(T=200, p=20, sparsity=sp, seed=r)
      reps[[r]] <- run_prediction(sim)
   }
   results_C[[label]] <- reps
}

# ============================================================================
#  SECTION 4: Results
# ============================================================================

summary_C <- do.call(rbind, lapply(names(results_C), function(label)
   summarize_reps(results_C[[label]], label)))

cat("\n============================================================\n")
cat("SIMULATION C RESULTS\n")
cat("Mean Squared Prediction Error (lower is better)\n")
cat("Hold-out 30% of variables per test time point\n")
cat("============================================================\n")
print(summary_C, row.names=FALSE)