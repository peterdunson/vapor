# ============================================================================
#  Bayesian Precision Factor Analysis for Time Series
#  Simulation B: Data Generated from VAR(1) with Sparse Precision Matrix
#                (no factor structure — model misspecification setting)
#
#  DGP: Omega is a directly sparse positive definite matrix,
#       NOT generated from a factor model. Tests robustness of
#       precision factor methods under misspecification.
#
#  Methods compared:
#    1. PrecFact (no VAR)
#    2. Glasso (no VAR)
#    3. VAR + Glasso
#    4. VAR + PrecFact (sequential)
#    5. Joint sampler (proposed)
#    6. Kastner: VAR + factor stochastic volatility (inverted Sigma)
#
#  Experiments:
#    A. Varying autocorrelation: rho in {0.50, 0.89, 0.95}, T=200, p=20
#    B. Varying sample size:     T in {100, 200, 500}, rho=0.95, p=20
#
#  5 replicates per scenario
# ============================================================================

library(Rcpp)
library(RcppArmadillo)
library(PrecFact)
library(bayesianVARs)
library(glasso)
library(MASS)

Rcpp::sourceCpp("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/var_prec.cpp")

# ============================================================================
#  SECTION 1: Helper functions
# ============================================================================

# ── Sparse PD precision matrix generator ─────────────────────────────────────
# Generates a directly sparse Omega with no factor structure
# Uses random graph with edge probability edge_prob
simulate_sparse_omega <- function(p, edge_prob=0.1, seed=NULL) {
   if (!is.null(seed)) set.seed(seed)
   
   # generate sparse symmetric matrix via random graph
   Omega <- matrix(0, p, p)
   for (i in seq_len(p-1)) {
      for (j in (i+1):p) {
         if (runif(1) < edge_prob) {
            val <- runif(1, 0.3, 0.7) * sample(c(-1,1), 1)
            Omega[i,j] <- val
            Omega[j,i] <- val
         }
      }
   }
   
   # make positive definite by adding scaled identity
   eigs    <- eigen(Omega, only.values=TRUE)$values
   min_eig <- min(eigs)
   if (min_eig <= 0) {
      Omega <- Omega + (abs(min_eig) + 0.5) * diag(p)
   } else {
      Omega <- Omega + 0.1 * diag(p)
   }
   
   Omega
}

# ── Data generator (simulation B DGP) ────────────────────────────────────────
simulate_var1_sparse <- function(T, p, target_radius=0.89,
                                 edge_prob=0.1, seed=NULL) {
   if (!is.null(seed)) set.seed(seed)
   
   Omega  <- simulate_sparse_omega(p, edge_prob=edge_prob)
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
        n_edges=sum(Omega[upper.tri(Omega)] != 0))
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
fit_precfact <- function(Y_input, nrun=5000, thin=5) {
   res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=5,
                   nrun=nrun, thin=thin, Y=scale(Y_input),
                   alldiff=0, a_dir=0.1, b_dir=0.1)
   fact_prec_processing(res, burnout=(nrun/2)/thin, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# ── VAR residuals ─────────────────────────────────────────────────────────────
get_var_residuals <- function(Y) {
   set.seed(99)
   p      <- ncol(Y)
   T_obs  <- nrow(Y)
   jitter <- ifelse(T_obs < 100, 1e-3, 1e-6)
   Y_jit  <- Y + matrix(rnorm(prod(dim(Y)), 0, jitter), T_obs, p)
   mod    <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
   B_hat  <- apply(mod$PHI, c(1,2), mean)
   list(resid=mod$Y - mod$X %*% B_hat, mod=mod)
}

# ── Joint sampler wrapper ─────────────────────────────────────────────────────
fit_joint <- function(Y, nrun=5000, thin=5) {
   set.seed(99)
   p      <- ncol(Y)
   T_obs  <- nrow(Y)
   # scale jitter with T to handle near-singular cases at small T
   jitter <- ifelse(T_obs < 100, 1e-3, 1e-6)
   Y_jit  <- Y + matrix(rnorm(prod(dim(Y)), 0, jitter), T_obs, p)
   res    <- var_precision_sampler(Y=Y_jit, lags=1, k=5,
                                   nrun=nrun, burn=nrun/2, thin=thin)
   fact_prec_processing(res, burnout=(nrun/2)/thin, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# ── Kastner: VAR + factor SV, invert Sigma ───────────────────────────────────
fit_kastner <- function(Y, Omega_true, k=5) {
   set.seed(99)
   p      <- ncol(Y)
   T_obs  <- nrow(Y)
   jitter <- ifelse(T_obs < 100, 1e-3, 1e-6)
   Y_jit  <- Y + matrix(rnorm(prod(dim(Y)), 0, jitter), T_obs, p)
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

# ── Run all 6 methods on one dataset ─────────────────────────────────────────
run_one <- function(sim) {
   results <- list()
   
   pf_raw <- fit_precfact(sim$Y)
   results[[1]] <- compute_metrics(pf_raw$Omega, pf_raw$P,
                                   sim$Omega_true, "PrecFact (no VAR)")
   
   gl_raw <- fit_glasso_oracle(sim$Y, sim$Omega_true)
   results[[2]] <- compute_metrics(gl_raw$Omega, gl_raw$P,
                                   sim$Omega_true, "Glasso (no VAR)")
   
   var_out <- get_var_residuals(sim$Y)
   
   gl_var <- fit_glasso_oracle(var_out$resid, sim$Omega_true)
   results[[3]] <- compute_metrics(gl_var$Omega, gl_var$P,
                                   sim$Omega_true, "VAR + Glasso")
   
   pf_var <- fit_precfact(var_out$resid)
   results[[4]] <- compute_metrics(pf_var$Omega, pf_var$P,
                                   sim$Omega_true, "VAR + PrecFact")
   
   joint <- fit_joint(sim$Y)
   results[[5]] <- compute_metrics(joint$Omega, joint$P,
                                   sim$Omega_true, "Joint sampler")
   
   kast <- fit_kastner(sim$Y, sim$Omega_true, k=5)
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
#  SECTION 2: Experiment A — varying autocorrelation
#  Fixed: T=200, p=20. Varying: rho in {0.50, 0.89, 0.95}. 5 replicates.
# ============================================================================

cat("\n============================================================\n")
cat("SIMULATION B, EXPERIMENT A: Varying autocorrelation\n")
cat("DGP: sparse Omega (no factor structure), T=200, p=20, 5 reps\n")
cat("============================================================\n")

radii <- c(0.10, 0.50, 0.75, 0.89, 0.95, 0.99)
n_reps <- 5
res_A  <- list()

for (rho in radii) {
   label <- paste0("rho=", rho)
   cat("\n  Scenario:", label, "\n")
   reps <- vector("list", n_reps)
   for (r in seq_len(n_reps)) {
      cat("    replicate", r, "/", n_reps, "\n")
      sim       <- simulate_var1_sparse(T=200, p=20, target_radius=rho,
                                        edge_prob=0.1, seed=r)
      cat("      true edges:", sim$n_edges, "\n")
      reps[[r]] <- run_one(sim)
   }
   res_A[[label]] <- reps
}

summary_A <- do.call(rbind, lapply(names(res_A), function(label)
   summarize_reps(res_A[[label]], label)))

cat("\n--- Simulation B, Experiment A: mean (sd) over 5 replicates ---\n")
print(summary_A, row.names=FALSE)

# ============================================================================
#  SECTION 3: Experiment B — varying sample size
#  Fixed: rho=0.95, p=20. Varying: T in {100, 200, 500}. 5 replicates.
# ============================================================================

cat("\n============================================================\n")
cat("SIMULATION B, EXPERIMENT B: Varying sample size\n")
cat("DGP: sparse Omega (no factor structure), rho=0.95, p=20, 5 reps\n")
cat("============================================================\n")

Ts    <- c(50, 100, 200, 500, 1000)
res_B <- list()

for (Ti in Ts) {
   label <- paste0("T=", Ti)
   cat("\n  Scenario:", label, "\n")
   reps <- vector("list", n_reps)
   for (r in seq_len(n_reps)) {
      cat("    replicate", r, "/", n_reps, "\n")
      sim       <- simulate_var1_sparse(T=Ti, p=20, target_radius=0.95,
                                        edge_prob=0.1, seed=r)
      reps[[r]] <- run_one(sim)
   }
   res_B[[label]] <- reps
}

summary_B <- do.call(rbind, lapply(names(res_B), function(label)
   summarize_reps(res_B[[label]], label)))

cat("\n--- Simulation B, Experiment B: mean (sd) over 5 replicates ---\n")
print(summary_B, row.names=FALSE)

# ============================================================================
#  SECTION 4: Final summary
# ============================================================================

cat("\n============================================================\n")
cat("FINAL SUMMARY — SIMULATION B\n")
cat("============================================================\n")
cat("\nExperiment A: varying autocorrelation (T=200, p=20)\n")
print(summary_A, row.names=FALSE)
cat("\nExperiment B: varying sample size (rho=0.95, p=20)\n")
print(summary_B, row.names=FALSE)