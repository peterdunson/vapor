# ============================================================================
#  Simulation D: varying signal strength in Lambda_0
#  VAR(1) with Omega = Lambda_0 Lambda_0^T + I_p, Simulation A DGP
#  25 replicates; paired per-replicate comparison against VAPOR
# ============================================================================

library(VAPOR)
library(PrecFact)
library(bayesianVARs)
library(glasso)
library(MASS)

#  in-repo copy; previously sourced from
#  ~/Desktop/precision_factor/PrecFact/vignettes/simulate_data_fxns.R
source("/Users/peterdunson/Desktop/vapor/simulation/simulate_data_fxns.R")
OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec"

N_REPS <- 25
NRUN   <- 5000
THIN   <- 5

add_jitter <- function(Y, seed, sd=1e-6) {
   set.seed(seed)
   Y + matrix(rnorm(prod(dim(Y)), 0, sd), nrow(Y), ncol(Y))
}

# ── Data generator ───────────────────────────────────────────────────────────
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

# ── Metrics (no rounding) ────────────────────────────────────────────────────
compute_metrics <- function(Omega_est, P_est, Omega_true, method) {
   G_true <- (Omega_true != 0)*1; diag(G_true) <- 0
   G_est  <- P_est*1;             diag(G_est)  <- 0
   TP <- sum(G_est==1 & G_true==1); FP <- sum(G_est==1 & G_true==0)
   TN <- sum(G_est==0 & G_true==0); FN <- sum(G_est==0 & G_true==1)
   data.frame(Method      = method,
              Sensitivity = TP/(TP+FN),
              Specificity = TN/(TN+FP),
              Frobenius   = norm(cov2cor(Omega_est) - cov2cor(Omega_true), "F"))
}

# ── Oracle Glasso ────────────────────────────────────────────────────────────
fit_glasso_oracle <- function(Y_input, Omega_true,
                              rho_grid=10^seq(-2, 0, length=20)) {
   S <- cov(scale(Y_input)); best <- list(frob=Inf)
   for (rho in rho_grid) {
      gl <- tryCatch(glasso(S, rho=rho, penalize.diagonal=FALSE),
                     error=function(e) NULL)
      if (is.null(gl)) next
      frob <- norm(cov2cor(gl$wi) - cov2cor(Omega_true), "F")
      if (frob < best$frob) best <- list(Omega=gl$wi, P=(gl$wi!=0)*1, frob=frob)
   }
   best
}

# ── PrecFact ─────────────────────────────────────────────────────────────────
fit_precfact <- function(Y_input, k=5) {
   res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=k,
                   nrun=NRUN, thin=THIN, Y=scale(Y_input),
                   alldiff=0, a_dir=0.1, b_dir=0.1)
   fact_prec_processing(res, burnout=(NRUN/2)/THIN, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

# ── VAR residuals ────────────────────────────────────────────────────────────
get_var_residuals <- function(Y, seed) {
   Y_jit <- add_jitter(Y, seed)
   mod   <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
   B_hat <- apply(mod$PHI, c(1,2), mean)
   list(resid=mod$Y - mod$X %*% B_hat)
}

# ── VAPOR ────────────────────────────────────────────────────────────────────
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

# ── Kastner ──────────────────────────────────────────────────────────────────
fit_kastner <- function(Y, Omega_true, k=5, seed=99) {
   Y_jit <- add_jitter(Y, seed)
   p <- ncol(Y)
   sigma_prior <- specify_prior_sigma(M=p, type="factor", factor_factors=k,
                                      factor_priorfacloadtype="colwiseng",
                                      factor_heteroskedastic=TRUE, quiet=TRUE)
   mod <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000,
               prior_sigma=sigma_prior, quiet=TRUE)
   vc        <- vcov(mod)
   Sigma_hat <- apply(vc[1,,,], c(1,2), mean)
   Omega_hat <- solve(Sigma_hat)
   G_true <- (Omega_true != 0)*1; diag(G_true) <- 0
   n_true <- sum(G_true)/2
   thresh <- quantile(abs(Omega_hat[upper.tri(Omega_hat)]),
                      1 - n_true/choose(p,2))
   list(Omega=Omega_hat, P=(abs(Omega_hat) > thresh)*1)
}

# ── One replicate, all six methods ───────────────────────────────────────────
run_one <- function(sim, seed) {
   r <- list()
   pf_raw <- fit_precfact(sim$Y)
   r[[1]] <- compute_metrics(pf_raw$Omega, pf_raw$P, sim$Omega_true, "PrecFact (no VAR)")
   gl_raw <- fit_glasso_oracle(sim$Y, sim$Omega_true)
   r[[2]] <- compute_metrics(gl_raw$Omega, gl_raw$P, sim$Omega_true, "Glasso (no VAR)")
   var_out <- get_var_residuals(sim$Y, seed)
   gl_var  <- fit_glasso_oracle(var_out$resid, sim$Omega_true)
   r[[3]] <- compute_metrics(gl_var$Omega, gl_var$P, sim$Omega_true, "VAR + Glasso")
   pf_var  <- fit_precfact(var_out$resid)
   r[[4]] <- compute_metrics(pf_var$Omega, pf_var$P, sim$Omega_true, "VAR + PrecFact")
   vp <- fit_vapor(sim$Y, seed=seed)
   r[[5]] <- compute_metrics(vp$Omega, vp$P, sim$Omega_true, "VAPOR")
   ks <- fit_kastner(sim$Y, sim$Omega_true, seed=seed)
   r[[6]] <- compute_metrics(ks$Omega, ks$P, sim$Omega_true, "Kastner (VAR+FSV)")
   do.call(rbind, r)
}

# ── Summaries ────────────────────────────────────────────────────────────────
summarize_reps <- function(rep_list, label) {
   all_m <- do.call(rbind, rep_list)
   do.call(rbind, lapply(unique(all_m$Method), function(m) {
      s <- all_m[all_m$Method == m, ]
      data.frame(Scenario=label, Method=m,
                 Sens=sprintf("%.3f (%.3f)", mean(s$Sensitivity), sd(s$Sensitivity)),
                 Spec=sprintf("%.3f (%.3f)", mean(s$Specificity), sd(s$Specificity)),
                 Frob=sprintf("%.3f (%.3f)", mean(s$Frobenius),   sd(s$Frobenius)))
   }))
}

paired_diffs <- function(rep_list, label) {
   fr <- t(sapply(rep_list, function(d) setNames(d$Frobenius, d$Method)))
   v  <- fr[, "VAPOR"]
   do.call(rbind, lapply(setdiff(colnames(fr), "VAPOR"), function(m) {
      d <- fr[, m] - v                           # positive favours VAPOR
      n <- length(d); se <- sd(d)/sqrt(n); tt <- mean(d)/se
      data.frame(Scenario=label, Method=m,
                 Diff=round(mean(d), 3), SE=round(se, 3),
                 t=round(tt, 2),
                 p=format.pval(2*pt(-abs(tt), df=n-1), digits=2, eps=1e-4),
                 Wins=sprintf("%d / %d", sum(d > 0), n))
   }))
}

# ============================================================================
#  Run: p=20, k=5, T=200, rho(B)=0.95, val in {0.5, 1.0, 1.5, 2.5, 4.0}
# ============================================================================
vals  <- c(0.5, 1.0, 1.5, 2.5, 4.0)
res_D <- list()
snr   <- list()

for (v in vals) {
   label <- paste0("val=", v)
   cat(sprintf("\n=== %s ===\n", label))
   reps <- vector("list", N_REPS)
   pcm  <- numeric(N_REPS)
   for (r in seq_len(N_REPS)) {
      cat(sprintf("\r  rep %d/%d", r, N_REPS)); flush.console()
      sim <- simulate_var1(T=200, p=20, target_radius=0.95,
                           k=5, pr=0.2, val=v, seed=r)
      pc <- -cov2cor(sim$Omega_true); diag(pc) <- 0
      pcm[r] <- mean(abs(pc[pc != 0]))
      reps[[r]] <- run_one(sim, seed=r)
   }
   cat("\n")
   res_D[[label]] <- reps
   snr[[label]]   <- pcm
   cat(sprintf("  mean |partial corr| over replicates: %.3f\n", mean(pcm)))
   saveRDS(res_D, file.path(path.expand(OUT), "sim_D_partial.rds"))
}

summary_D <- do.call(rbind, lapply(names(res_D), function(l) summarize_reps(res_D[[l]], l)))
paired_D  <- do.call(rbind, lapply(names(res_D), function(l) paired_diffs(res_D[[l]], l)))

snr_tab <- data.frame(Scenario = names(snr),
                      MeanAbsPartialCorr = round(sapply(snr, mean), 3),
                      SD = round(sapply(snr, sd), 3))

cat("\n--- Simulation D: mean (sd) over", N_REPS, "replicates ---\n")
print(summary_D, row.names=FALSE)
cat("\n--- Paired Frobenius differences vs VAPOR (positive favours VAPOR) ---\n")
print(paired_D, row.names=FALSE)
cat("\n--- Signal strength by scenario ---\n")
print(snr_tab, row.names=FALSE)

saveRDS(list(raw=res_D, summary=summary_D, paired=paired_D,
             snr=snr, snr_summary=snr_tab, n_reps=N_REPS),
        file.path(path.expand(OUT), "sim_D_results.rds"))
cat("\nSaved to sim_D_results.rds\n")