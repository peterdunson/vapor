# ============================================================================
#  Bayesian Precision Factor Analysis for Time Series
#  Simulation E: Varying Dimension p
#                Data generated from proposed model (Simulation A DGP)
#
#  DGP: VAR(1) with Omega = Lambda_0*Lambda_0^T + I_p (RSM structure)
#
#  Methods compared:
#    1. PrecFact (no VAR)
#    2. Glasso (no VAR)
#    3. VAR + Glasso
#    4. VAR + PrecFact (sequential)
#    5. Joint sampler (proposed)
#    6. Kastner: VAR + factor stochastic volatility (inverted Sigma)
#
#  Settings:
#    p=10,  k=2,  pr=0.25, T=200, rho=0.95
#    p=20,  k=5,  pr=0.20, T=200, rho=0.95
#    p=50,  k=8,  pr=0.08, T=200, rho=0.95
#    p=100, k=12, pr=0.04, T=500, rho=0.95  (larger T for identifiability)
#
#  val=1.5 (baseline signal strength), 5 replicates per scenario
# ============================================================================

library(Rcpp)
library(RcppArmadillo)
library(PrecFact)
library(bayesianVARs)
library(glasso)
library(MASS)

Rcpp::sourceCpp("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/var_prec.cpp")
source("~/Desktop/precision_factor/PrecFact/vignettes/simulate_data_fxns.R")

# ============================================================================
#  SECTION 1: Helper functions
# ============================================================================

# ── Data generator ────────────────────────────────────────────────────────────
simulate_var1 <- function(T, p, target_radius=0.95, k=5,
                          pr=0.2, val=1.5, seed=NULL) {
   if (!is.null(seed)) set.seed(seed)
   lambda <- simulate_lambda_sparse(k, p, pr=pr, val=val)
   Omega  <- tcrossprod(lambda) + diag(p)
   Sigma  <- chol2inv(chol(Omega))
   B      <- matrix(0, p, p)
   B[sample(p^2, floor(0.1*p^2))] <- rnorm(floor(0.1*p^2), 0, 0.3)
   sr     <- max(abs(eigen(B, only.values=TRUE)$values))
   B      <- B * (target_radius / sr)
   Y      <- matrix(0, T, p)
   eps    <- MASS::mvrnorm(T, rep(0,p), Sigma)
   Y[1,]  <- eps[1,]
   for (t in 2:T) Y[t,] <- B %*% Y[t-1,] + eps[t,]
   list(Y=Y, B=B, Omega_true=Omega, Sigma_true=Sigma,
        spectral_radius=max(abs(eigen(B, only.values=TRUE)$values)),
        n_edges=sum((Omega - diag(diag(Omega))) != 0)/2)
}

# ── Adaptive jitter ───────────────────────────────────────────────────────────
make_jitter <- function(Y) {
   T_obs <- nrow(Y); p <- ncol(Y)
   jitter <- ifelse(T_obs / p < 5, 1e-3, 1e-6)
   Y + matrix(rnorm(prod(dim(Y)), 0, jitter), T_obs, p)
}

# ── Metrics ───────────────────────────────────────────────────────────────────
compute_metrics <- function(Omega_est, P_est, Omega_true, method) {
   G_true <- (Omega_true != 0)*1; diag(G_true) <- 0
   G_est  <- P_est*1;             diag(G_est)  <- 0
   TP <- sum(G_est==1 & G_true==1)
   FP <- sum(G_est==1 & G_true==0)
   TN <- sum(G_est==0 & G_true==0)
   FN <- sum(G_est==0 & G_true==1)
   data.frame(
      Method      = method,
      Sensitivity = round(TP/(TP+FN), 3),
      Specificity = round(TN/(TN+FP), 3),
      Frobenius   = round(norm(cov2cor(Omega_est) - cov2cor(Omega_true), "F"), 3)
   )
}

# ── Oracle Glasso ─────────────────────────────────────────────────────────────
fit_glasso_oracle <- function(Y_input, Omega_true,
                              rho_grid=10^seq(-2, 0, length=20)) {
   S    <- cov(scale(Y_input))
   best <- list(frob=Inf)
   for (rho in rho_grid) {
      gl <- tryCatch(glasso(S, rho=rho, penalize.diagonal=FALSE),
                     error=function(e) NULL)
      if (is.null(gl)) next
      frob <- norm(cov2cor(gl$wi) - cov2cor(Omega_true), "F")
      if (frob < best$frob)
         best <- list(Omega=gl$wi, P=(gl$wi!=0)*1, frob=frob)
   }
   best
}

# ── PrecFact wrapper ──────────────────────────────────────────────────────────
fit_precfact <- function(Y_input, k=5, nrun=5000, thin=5) {
   res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=k,
                   nrun=nrun, thin=thin, Y=scale(Y_input),
                   alldiff=0, a_dir=0.1, b_dir=0.1)
   fact_prec_processing(res, burnout=(nrun/2)/thin, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# ── VAR residuals ─────────────────────────────────────────────────────────────
get_var_residuals <- function(Y) {
   set.seed(99)
   Y_jit <- make_jitter(Y)
   mod   <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
   B_hat <- apply(mod$PHI, c(1,2), mean)
   list(resid=mod$Y - mod$X %*% B_hat)
}

# ── Joint sampler wrapper ─────────────────────────────────────────────────────
fit_joint <- function(Y, k=5, nrun=5000, thin=5) {
   set.seed(99)
   Y_jit <- make_jitter(Y)
   res   <- var_precision_sampler(Y=Y_jit, lags=1, k=k,
                                  nrun=nrun, burn=nrun/2, thin=thin)
   fact_prec_processing(res, burnout=(nrun/2)/thin, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# ── Kastner: VAR + factor SV ──────────────────────────────────────────────────
fit_kastner <- function(Y, Omega_true, k=5) {
   set.seed(99)
   p     <- ncol(Y)
   Y_jit <- make_jitter(Y)
   sigma_prior <- specify_prior_sigma(
      M                       = p,
      type                    = "factor",
      factor_factors          = k,
      factor_priorfacloadtype = "colwiseng",
      factor_heteroskedastic  = TRUE,
      quiet                   = TRUE
   )
   mod       <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000,
                     prior_sigma=sigma_prior, quiet=TRUE)
   vc        <- vcov(mod)
   Sigma_hat <- apply(vc[1,,,], c(1,2), mean)
   Omega_hat <- solve(Sigma_hat)
   G_true    <- (Omega_true != 0)*1; diag(G_true) <- 0
   n_true    <- sum(G_true)/2
   thresh    <- quantile(abs(Omega_hat[upper.tri(Omega_hat)]),
                         1 - n_true/choose(p,2))
   P_hat     <- (abs(Omega_hat) > thresh)*1
   list(Omega=Omega_hat, P=P_hat)
}

# ── Run all methods on one dataset ────────────────────────────────────────────
run_one <- function(sim, k=5) {
   results <- list()
   
   pf_raw <- fit_precfact(sim$Y, k=k)
   results[[1]] <- compute_metrics(pf_raw$Omega, pf_raw$P,
                                   sim$Omega_true, "PrecFact (no VAR)")
   
   gl_raw <- fit_glasso_oracle(sim$Y, sim$Omega_true)
   results[[2]] <- compute_metrics(gl_raw$Omega, gl_raw$P,
                                   sim$Omega_true, "Glasso (no VAR)")
   
   var_out <- get_var_residuals(sim$Y)
   
   gl_var <- fit_glasso_oracle(var_out$resid, sim$Omega_true)
   results[[3]] <- compute_metrics(gl_var$Omega, gl_var$P,
                                   sim$Omega_true, "VAR + Glasso")
   
   pf_var <- fit_precfact(var_out$resid, k=k)
   results[[4]] <- compute_metrics(pf_var$Omega, pf_var$P,
                                   sim$Omega_true, "VAR + PrecFact")
   
   joint <- fit_joint(sim$Y, k=k)
   results[[5]] <- compute_metrics(joint$Omega, joint$P,
                                   sim$Omega_true, "Joint sampler")
   
   kast <- fit_kastner(sim$Y, sim$Omega_true, k=k)
   results[[6]] <- compute_metrics(kast$Omega, kast$P,
                                   sim$Omega_true, "Kastner (VAR+FSV)")
   
   do.call(rbind, results)
}

# ── Summarize replicates ──────────────────────────────────────────────────────
summarize_reps <- function(rep_list, scenario_label) {
   all_m   <- do.call(rbind, rep_list)
   methods <- unique(all_m$Method)
   rows <- lapply(methods, function(m) {
      sub <- all_m[all_m$Method == m, ]
      data.frame(
         Scenario = scenario_label,
         Method   = m,
         Sens     = sprintf("%.3f (%.3f)", mean(sub$Sensitivity), sd(sub$Sensitivity)),
         Spec     = sprintf("%.3f (%.3f)", mean(sub$Specificity), sd(sub$Specificity)),
         Frob     = sprintf("%.3f (%.3f)", mean(sub$Frobenius),   sd(sub$Frobenius))
      )
   })
   do.call(rbind, rows)
}

# ============================================================================
#  SECTION 2: Simulation E — varying dimension
# ============================================================================

cat("\n============================================================\n")
cat("SIMULATION E: Varying Dimension p\n")
cat("DGP: VAR(1) with sparse precision factor Omega\n")
cat("Fixed: rho=0.95, val=1.5, 5 replicates\n")
cat("============================================================\n")

# scenario settings: p, k, pr, T
settings <- list(
   list(p=10,  k=2,  pr=0.25, T=200, label="p=10,  k=2"),
   list(p=20,  k=5,  pr=0.20, T=200, label="p=20,  k=5"),
   list(p=50,  k=8,  pr=0.08, T=200, label="p=50,  k=8"),
   list(p=100, k=12, pr=0.04, T=500, label="p=100, k=12")
)

n_reps <- 5
res_E  <- list()

for (s in settings) {
   cat("\n  Scenario:", s$label, "| T=", s$T, "\n")
   reps <- vector("list", n_reps)
   
   for (r in seq_len(n_reps)) {
      cat("    replicate", r, "/", n_reps, "\n")
      sim <- simulate_var1(T=s$T, p=s$p, target_radius=0.95,
                           k=s$k, pr=s$pr, val=1.5, seed=r)
      cat("      true edges:", sim$n_edges, "\n")
      reps[[r]] <- run_one(sim, k=s$k)
   }
   res_E[[s$label]] <- reps
}

# ============================================================================
#  SECTION 3: Results
# ============================================================================

summary_E <- do.call(rbind, lapply(names(res_E), function(label)
   summarize_reps(res_E[[label]], label)))

cat("\n============================================================\n")
cat("SIMULATION E RESULTS\n")
cat("Mean (sd) over 5 replicates\n")
cat("============================================================\n")
print(summary_E, row.names=FALSE)