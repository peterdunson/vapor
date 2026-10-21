# ============================================================================
#  Simulation C: Kastner & Huber (2020) DGP, VAR(1) with factor SV errors
#  Held-out prediction MSPE; 25 replicates; paired comparison against VAPOR
# ============================================================================

library(VAPOR)
library(PrecFact)
library(bayesianVARs)
library(MASS)

OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec"

N_REPS <- 25
NRUN   <- 5000
THIN   <- 5

add_jitter <- function(Y, seed, sd=1e-6) {
   set.seed(seed)
   Y + matrix(rnorm(prod(dim(Y)), 0, sd), nrow(Y), ncol(Y))
}

# ── DGP: Kastner & Huber (2020) ──────────────────────────────────────────────
simulate_kastner_dgp <- function(T, p, sparsity=0.10, seed=NULL) {
   if (!is.null(seed)) set.seed(seed)
   if (sparsity <= 0.01) { mu_O <- 0.01; sd_O <- 0.01 } else { mu_O <- 0.1; sd_O <- 0.1 }
   mu_D <- 0.15; sd_D <- 0.15
   
   B <- matrix(0, p, p)
   diag(B) <- rnorm(p, mu_D, sd_D)
   for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i != j && runif(1) < sparsity) B[i,j] <- rnorm(1, mu_O, sd_O)
   }
   sr <- max(abs(eigen(B, only.values=TRUE)$values))
   if (sr >= 1) B <- B * 0.95 / sr
   
   lambda_fac <- rnorm(p, 0.001, 0.001)
   mu_sv    <- rep(-12, p)
   rho_sv   <- runif(p, 0.85, 0.98)
   sigma_sv <- runif(p, 0.1, 0.3)
   rho_fac  <- 0.99; sigma_fac <- 0.1
   
   h_fac <- numeric(T); h_idi <- matrix(0, T, p)
   h_fac[1] <- rnorm(1, 0, sigma_fac / sqrt(1 - rho_fac^2))
   for (t in 2:T) h_fac[t] <- rho_fac * h_fac[t-1] + rnorm(1, 0, sigma_fac)
   for (j in seq_len(p)) {
      h_idi[1,j] <- rnorm(1, mu_sv[j], sigma_sv[j] / sqrt(1 - rho_sv[j]^2))
      for (t in 2:T)
         h_idi[t,j] <- mu_sv[j] + rho_sv[j]*(h_idi[t-1,j]-mu_sv[j]) +
            rnorm(1, 0, sigma_sv[j])
   }
   
   epsilon <- matrix(0, T, p)
   for (t in seq_len(T)) {
      Sigma_t <- tcrossprod(lambda_fac) * exp(h_fac[t]) + diag(exp(h_idi[t,]))
      epsilon[t,] <- MASS::mvrnorm(1, rep(0,p), Sigma_t)
   }
   
   Y <- matrix(0, T, p); Y[1,] <- epsilon[1,]
   for (t in 2:T) Y[t,] <- B %*% Y[t-1,] + epsilon[t,]
   
   Sigma_avg <- Reduce("+", lapply(seq_len(T), function(t)
      tcrossprod(lambda_fac) * exp(h_fac[t]) + diag(exp(h_idi[t,])))) / T
   list(Y=Y, B=B, Omega_avg=solve(Sigma_avg), Sigma_avg=Sigma_avg,
        lambda_fac=lambda_fac, h_fac=h_fac, h_idi=h_idi)
}

# ── Conditional mean prediction ──────────────────────────────────────────────
predict_held_out <- function(Omega, y_obs, obs_idx, hid_idx) {
   as.vector(-solve(Omega[hid_idx, hid_idx]) %*% Omega[hid_idx, obs_idx] %*% y_obs)
}

# ── Fitters ──────────────────────────────────────────────────────────────────
fit_precfact <- function(Y_input, k=5) {
   res <- prec_est(as=1.0, bs=1.0, a=0.5, b=0.5, k=k,
                   nrun=NRUN, thin=THIN, Y=scale(Y_input),
                   alldiff=0, a_dir=0.1, b_dir=0.1)
   fact_prec_processing(res, burnout=(NRUN/2)/THIN, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

get_var_residuals <- function(Y, seed) {
   Y_jit <- add_jitter(Y, seed)
   mod   <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
   B_hat <- apply(mod$PHI, c(1,2), mean)
   list(resid=mod$Y - mod$X %*% B_hat)
}

fit_vapor <- function(Y, k=5, seed=99) {
   has_init <- "init" %in% names(formals(vapor))
   fit <- NULL
   for (attempt in 0:2) {
      Y_use <- if (attempt == 0) Y else add_jitter(Y, seed + 1000*attempt, 1e-4)
      args  <- list(Y_use, lags=1, k=k, rank="full", prior_B="dl",
                    nrun=NRUN, burn=NRUN/2, thin=THIN, seed=seed + attempt)
      if (has_init) args$init <- "random"
      fit <- tryCatch(do.call(vapor, args), error=function(e) NULL)
      if (!is.null(fit)) break
   }
   if (is.null(fit)) stop("vapor failed after 3 attempts, seed ", seed)
   res_like <- list(Lambda=fit$Lambda, residuals=fit$residual_prec)
   out <- fact_prec_processing(res_like, burnout=0, alpha=0.1, beta=0.85,
                               epsilon_grid=seq(0.05,0.1,length.out=20))
   diag(out$P) <- 0; out
}

# one factor here, matching the DGP; this is the oracle specification
fit_kastner <- function(Y, seed) {
   Y_jit <- add_jitter(Y, seed)
   p <- ncol(Y)
   sigma_prior <- specify_prior_sigma(M=p, type="factor", factor_factors=1L,
                                      factor_priorfacloadtype="colwiseng",
                                      factor_heteroskedastic=TRUE, quiet=TRUE)
   mod       <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000,
                     prior_sigma=sigma_prior, quiet=TRUE)
   vc        <- vcov(mod)
   Sigma_hat <- apply(vc[1,,,], c(1,2), mean)
   list(Omega=solve(Sigma_hat))
}

# ── One replicate ────────────────────────────────────────────────────────────
METHODS <- c("PrecFact (no VAR)", "VAR + PrecFact", "VAPOR",
             "Kastner (VAR+FSV)", "Sample covariance")

run_prediction <- function(sim, seed, T_train_frac=0.8) {
   Y       <- scale(sim$Y)
   T_obs   <- nrow(Y); p <- ncol(Y)
   T_train <- floor(T_train_frac * T_obs)
   Y_train <- Y[1:T_train, ]; Y_test <- Y[(T_train+1):T_obs, ]
   n_hid   <- floor(0.3 * p)
   
   pf_raw  <- fit_precfact(Y_train)
   var_out <- get_var_residuals(Y_train, seed)
   pf_var  <- fit_precfact(var_out$resid)
   vp      <- fit_vapor(Y_train, seed=seed)
   kast    <- fit_kastner(Y_train, seed)
   Omega_samp <- solve(cov(Y_train))
   
   # hold-out design varies across replicates, shared across methods within one
   set.seed(1000 + seed)
   T_test <- nrow(Y_test)
   mspe <- matrix(0, T_test, length(METHODS), dimnames=list(NULL, METHODS))
   for (t in seq_len(T_test)) {
      y_t     <- Y_test[t,]
      hid_idx <- sort(sample(p, n_hid))
      obs_idx <- setdiff(seq_len(p), hid_idx)
      y_obs   <- y_t[obs_idx]; y_true <- y_t[hid_idx]
      mspe[t,] <- c(
         mean((y_true - predict_held_out(pf_raw$Omega, y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(pf_var$Omega, y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(vp$Omega,     y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(kast$Omega,   y_obs, obs_idx, hid_idx))^2),
         mean((y_true - predict_held_out(Omega_samp,   y_obs, obs_idx, hid_idx))^2))
   }
   colMeans(mspe)
}

# ── Summaries ────────────────────────────────────────────────────────────────
summarize_reps <- function(rep_list, label) {
   mat <- do.call(rbind, rep_list)
   data.frame(Scenario=label, Method=colnames(mat),
              MSPE=sprintf("%.4f (%.4f)", colMeans(mat), apply(mat, 2, sd)))
}

paired_diffs <- function(rep_list, label) {
   mat <- do.call(rbind, rep_list)
   v   <- mat[, "VAPOR"]
   do.call(rbind, lapply(setdiff(colnames(mat), "VAPOR"), function(m) {
      d <- mat[, m] - v                          # positive favours VAPOR
      n <- length(d); se <- sd(d)/sqrt(n); tt <- mean(d)/se
      data.frame(Scenario=label, Method=m,
                 Diff=round(mean(d), 4), SE=round(se, 4),
                 t=round(tt, 2),
                 p=format.pval(2*pt(-abs(tt), df=n-1), digits=2, eps=1e-4),
                 Wins=sprintf("%d / %d", sum(d > 0), n))
   }))
}

# ============================================================================
#  Run: p=20, T=200, sparsity in {0.01, 0.10, 0.80}
# ============================================================================
sparsity_levels <- c(0.01, 0.10, 0.80)
sparsity_labels <- c("Sparse (1%)", "Intermediate (10%)", "Dense (80%)")
results_C <- list()

for (s in seq_along(sparsity_levels)) {
   sp <- sparsity_levels[s]; label <- sparsity_labels[s]
   cat(sprintf("\n=== %s ===\n", label))
   reps <- vector("list", N_REPS)
   for (r in seq_len(N_REPS)) {
      cat(sprintf("\r  rep %d/%d", r, N_REPS)); flush.console()
      sim <- simulate_kastner_dgp(T=200, p=20, sparsity=sp, seed=r)
      reps[[r]] <- run_prediction(sim, seed=r)
   }
   cat("\n")
   results_C[[label]] <- reps
   saveRDS(results_C, file.path(path.expand(OUT), "sim_C_partial.rds"))
}

summary_C <- do.call(rbind, lapply(names(results_C), function(l) summarize_reps(results_C[[l]], l)))
paired_C  <- do.call(rbind, lapply(names(results_C), function(l) paired_diffs(results_C[[l]], l)))

cat("\n--- Simulation C: held-out MSPE, mean (sd) over", N_REPS, "replicates ---\n")
print(summary_C, row.names=FALSE)
cat("\n--- Paired MSPE differences vs VAPOR (positive favours VAPOR) ---\n")
print(paired_C, row.names=FALSE)

saveRDS(list(raw=results_C, summary=summary_C, paired=paired_C, n_reps=N_REPS),
        file.path(path.expand(OUT), "sim_C_results.rds"))
cat("\nSaved to sim_C_results.rds\n")