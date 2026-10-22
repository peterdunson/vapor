# ============================================================================
#  Simulation G (T<p): reduced-rank VAPOR under approximately low-rank dynamics
#  25 replicates; paired per-replicate comparison against VAPOR (dl, rank 3)
# ============================================================================

library(VAPOR)
library(bayesianVARs)
library(glasso)
library(MASS)
library(pheatmap)
library(gridExtra)
library(grid)

#  in-repo copy; previously sourced from
#  ~/Desktop/precision_factor/PrecFact/vignettes/simulate_data_fxns.R
source("/Users/peterdunson/Desktop/vapor/simulation/simulate_data_fxns.R")

R0           <- 3
PERTURB_FRAC <- 0.05
SPEC_RADIUS  <- 0.7
FIT_RANKS    <- c(3)
K_FACTORS    <- 5
N_SEEDS      <- 25
NRUN         <- 5000
THIN         <- 5

REF_METHOD <- sprintf("VAPOR dl r=%d", FIT_RANKS[1])

OUTDIR      <- "/Users/peterdunson/Desktop/vapor/simulation/results/fits"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
HEATMAP_DIR <- file.path("/Users/peterdunson/Desktop/vapor/simulation/results/figures", "sim_g_heatmaps")
dir.create(HEATMAP_DIR, recursive = TRUE, showWarnings = FALSE)

#  p = 200, T = 140 was dropped: it ran 17 h to replicate 19 of 25 and was
#  cancelled.  Simulation F is reported on three dimensions.
SETTINGS <- list(
   list(p=50,  T=35,  block="T<p sweep (T/p=0.7)"),
   list(p=100, T=70,  block="T<p sweep (T/p=0.7)"),
   list(p=150, T=105, block="T<p sweep (T/p=0.7)")
)

# ============================================================================
#  DGP
# ============================================================================

simulate_lowrank_calibrated <- function(T, p, target_radius = SPEC_RADIUS,
                                        k = K_FACTORS, pr = 0.15, val = 1.5,
                                        r0 = R0, perturb_frac = PERTURB_FRAC, seed = NULL) {
   if (!is.null(seed)) set.seed(seed)
   lambda <- simulate_lambda_sparse(k, p, pr = pr, val = val)
   Omega  <- tcrossprod(lambda) + diag(p)
   Sigma  <- chol2inv(chol(Omega))
   U0 <- matrix(rnorm(p * r0), p, r0); V0 <- matrix(rnorm(p * r0), p, r0)
   B_lowrank <- U0 %*% t(V0)
   N <- matrix(rnorm(p * p), p, p); N <- N * (perturb_frac * norm(B_lowrank,"F")/norm(N,"F"))
   B <- B_lowrank + N
   sr <- max(abs(eigen(B, only.values=TRUE)$values)); B <- B * (target_radius/sr)
   Y <- matrix(0,T,p); e <- MASS::mvrnorm(T, rep(0,p), Sigma); Y[1,] <- e[1,]
   for (t in 2:T) Y[t,] <- B %*% Y[t-1,] + e[t,]
   sv <- svd(B)$d
   list(Y=Y, B=B, Omega_true=Omega,
        eff_rank_90 = which(cumsum(sv^2)/sum(sv^2) >= 0.90)[1],
        n_edges = sum((Omega - diag(diag(Omega))) != 0)/2)
}

make_jitter <- function(Y, seed) {
   set.seed(seed)
   T_obs <- nrow(Y); p <- ncol(Y)
   sd_j  <- if (T_obs / p < 5) 1e-3 else 1e-6
   Y + matrix(rnorm(prod(dim(Y)), 0, sd_j), T_obs, p)
}

to_pc <- function(Om) { pc <- -cov2cor(Om); diag(pc) <- NA; pc }

b_rel_frob <- function(B_est, B_true) {
   if (is.null(B_est)) return(NA)
   p <- nrow(B_true)
   cand <- list(B_est, t(B_est))
   for (M in cand) {
      if (nrow(M) == p && ncol(M) == p)      return(norm(M - B_true, "F") / norm(B_true, "F"))
      if (nrow(M) == p + 1 && ncol(M) == p)  return(norm(M[1:p, ] - B_true, "F") / norm(B_true, "F"))
      if (nrow(M) == p && ncol(M) == p + 1)  return(norm(M[, 1:p] - B_true, "F") / norm(B_true, "F"))
   }
   NA
}

omega_metrics <- function(Omega_est, Omega_true) {
   if (is.null(Omega_est)) return(c(Frob=NA, Sens=NA, Spec=NA))
   frob <- norm(cov2cor(Omega_est) - cov2cor(Omega_true), "F")
   ut <- upper.tri(Omega_true)
   g_true <- (Omega_true[ut] != 0)
   n_true <- sum(g_true)
   pc_est <- abs(to_pc(Omega_est)); diag(pc_est) <- 0
   vals   <- pc_est[ut]
   thresh <- sort(vals, decreasing=TRUE)[n_true]
   g_est  <- (vals >= thresh)
   TP <- sum(g_est & g_true); FP <- sum(g_est & !g_true)
   TN <- sum(!g_est & !g_true); FN <- sum(!g_est & g_true)
   c(Frob=frob, Sens=TP/(TP+FN), Spec=TN/(TN+FP))
}

# ============================================================================
#  fit all methods on one dataset
# ============================================================================

fit_all <- function(sim, p, seed) {
   B_true <- sim$B
   rows <- list()
   pc_panels <- list(true = to_pc(sim$Omega_true))
   
   for (rk in FIT_RANKS) {
      fit_dl <- vapor(sim$Y, lags=1, k=K_FACTORS, rank=rk, prior_B="dl",
                      nrun=NRUN, burn=NRUN/2, thin=THIN, seed=seed)
      fit_rg <- vapor(sim$Y, lags=1, k=K_FACTORS, rank=rk, prior_B="ridge",
                      nrun=NRUN, burn=NRUN/2, thin=THIN, seed=seed)
      B_dl <- apply(fit_dl$B, c(1,2), mean)
      B_rg <- apply(fit_rg$B, c(1,2), mean)
      
      rows[[sprintf("VAPOR dl r=%d", rk)]] <-
         c(omega_metrics(fit_dl$Omega, sim$Omega_true), B_Frob=b_rel_frob(B_dl, B_true))
      rows[[sprintf("VAPOR ridge r=%d", rk)]] <-
         c(omega_metrics(fit_rg$Omega, sim$Omega_true), B_Frob=b_rel_frob(B_rg, B_true))
      
      pc_panels[[sprintf("dl_r%d", rk)]] <- to_pc(fit_dl$Omega)
      pc_panels[[sprintf("rg_r%d", rk)]] <- to_pc(fit_rg$Omega)
   }
   
   # Kastner
   sigma_prior <- specify_prior_sigma(M=p, type="factor", factor_factors=K_FACTORS,
                                      factor_priorfacloadtype="colwiseng",
                                      factor_heteroskedastic=TRUE, quiet=TRUE)
   mod_k <- tryCatch(bvar(make_jitter(sim$Y, seed), lags=1L, draws=2000, burnin=1000,
                          prior_sigma=sigma_prior, quiet=TRUE),
                     error=function(e){cat("    Kastner failed\n"); NULL})
   if (!is.null(mod_k)) {
      Omega_kast <- solve(apply(vcov(mod_k)[1,,,], c(1,2), mean))
      B_kast     <- apply(mod_k$PHI, c(1,2), mean)
   } else { Omega_kast <- NULL; B_kast <- NULL }
   rows[["Kastner"]] <- c(omega_metrics(Omega_kast, sim$Omega_true), B_Frob=NA)
   
   # Glasso
   S <- cov(scale(sim$Y)); best <- list(frob=Inf, Omega=diag(p))
   for (rho in 10^seq(-2.5,0,length=25)) {
      gl <- tryCatch(glasso(S, rho=rho, penalize.diagonal=FALSE), error=function(e) NULL)
      if (is.null(gl)) next
      fr <- norm(cov2cor(gl$wi) - cov2cor(sim$Omega_true), "F")
      if (fr < best$frob) best <- list(Omega=gl$wi, frob=fr)
   }
   Omega_glasso <- best$Omega
   rows[["Glasso"]] <- c(omega_metrics(Omega_glasso, sim$Omega_true), B_Frob=NA)
   
   pc_panels[["kast"]]   <- if (is.null(Omega_kast)) matrix(NA,p,p) else to_pc(Omega_kast)
   pc_panels[["glasso"]] <- to_pc(Omega_glasso)
   
   list(metrics = do.call(rbind, rows), pc = pc_panels)
}

# ============================================================================
#  paired differences
# ============================================================================

paired_diffs <- function(seed_metrics, label, ref = REF_METHOD) {
   fr <- t(sapply(seed_metrics, function(M) M[, "Frob"]))
   v  <- fr[, ref]
   do.call(rbind, lapply(setdiff(colnames(fr), ref), function(m) {
      d <- fr[, m] - v                                      # positive favours VAPOR
      d <- d[is.finite(d)]
      n <- length(d)
      if (n < 2) return(NULL)
      se <- sd(d)/sqrt(n); tt <- mean(d)/se
      data.frame(Scenario=label, Method=m,
                 Diff=round(mean(d), 3), SE=round(se, 3), t=round(tt, 2),
                 p=format.pval(2*pt(-abs(tt), df=n-1), digits=2, eps=1e-4),
                 Wins=sprintf("%d / %d", sum(d > 0), n))
   }))
}

# ============================================================================
#  heatmap panel builder
# ============================================================================

pal  <- colorRampPalette(c("blue","white","red"))(100)
brks <- seq(-1, 1, length.out = 101)
mk <- function(pc, title, legend=FALSE)
   pheatmap(pc, main=title, cluster_rows=FALSE, cluster_cols=FALSE,
            color=pal, breaks=brks, na_col="grey90",
            show_rownames=FALSE, show_colnames=FALSE, legend=legend, silent=TRUE)$gtable

save_heatmap <- function(pcs, tag, outfile) {
   panels <- list(mk(pcs$true, "TRUE"))
   for (rk in FIT_RANKS) {
      panels[[length(panels)+1]] <- mk(pcs[[sprintf("dl_r%d", rk)]], sprintf("VAPOR dl r=%d", rk))
      panels[[length(panels)+1]] <- mk(pcs[[sprintf("rg_r%d", rk)]], sprintf("VAPOR ridge r=%d", rk))
   }
   panels[[length(panels)+1]] <- mk(pcs$kast, "Kastner")
   panels[[length(panels)+1]] <- mk(pcs$glasso, "Glasso", legend=TRUE)
   ncol <- length(panels)
   png(outfile, width=380*ncol, height=450)
   grid.arrange(grobs=panels, ncol=ncol, top=textGrob(tag, gp=gpar(fontsize=15, fontface="bold")))
   dev.off()
}

# ============================================================================
#  run
# ============================================================================

all_rows   <- list()
all_paired <- list()
raw_G      <- list()

for (s in SETTINGS) {
   p <- s$p; Tval <- s$T
   label <- sprintf("p=%d, T=%d", p, Tval)
   cat(sprintf("\n=== %s  (T/p=%.1f)  [%s] ===\n", label, Tval/p, s$block))
   
   seed_metrics <- vector("list", N_SEEDS)
   for (sd_i in seq_len(N_SEEDS)) {
      cat(sprintf("\r  seed %d/%d", sd_i, N_SEEDS)); flush.console()
      sim <- simulate_lowrank_calibrated(T=Tval, p=p, seed=sd_i)
      if (sd_i == 1)
         cat(sprintf("\n    (Omega edges=%d | B eff rank@90%%=%d | true rank=%d)\n",
                     sim$n_edges, sim$eff_rank_90, qr(sim$B)$rank))
      f <- fit_all(sim, p, seed=sd_i)
      seed_metrics[[sd_i]] <- f$metrics
      
      if (sd_i == 1) {
         tag <- sprintf("p=%d  T=%d  (T/p=%.1f)", p, Tval, Tval/p)
         outfile <- file.path(HEATMAP_DIR, sprintf("sim_g_heatmap_p%d_T%d.png", p, Tval))
         save_heatmap(f$pc, tag, outfile)
         cat(sprintf("    heatmap -> %s\n", outfile))
      }
   }
   cat("\n")
   raw_G[[label]] <- seed_metrics
   
   methods <- rownames(seed_metrics[[1]])
   cols    <- colnames(seed_metrics[[1]])
   agg <- sapply(cols, function(cc)
      sapply(methods, function(mm)
         mean(sapply(seed_metrics, function(M) M[mm, cc]), na.rm=TRUE)))
   agg_sd <- sapply(cols, function(cc)
      sapply(methods, function(mm)
         sd(sapply(seed_metrics, function(M) M[mm, cc]), na.rm=TRUE)))
   
   ms <- function(cc) ifelse(is.na(agg[,cc]), "--",
                             sprintf("%.3f (%.3f)", agg[,cc], agg_sd[,cc]))
   disp <- data.frame(Method=methods, Omega_Frob=ms("Frob"), B_Frob=ms("B_Frob"),
                      Sens=ms("Sens"), Spec=ms("Spec"))
   cat(sprintf("  RESULTS %s, mean (sd) over %d replicates:\n", label, N_SEEDS))
   print(disp, row.names=FALSE)
   
   pd <- paired_diffs(seed_metrics, label)
   cat(sprintf("\n  Paired Omega-Frobenius differences vs %s (positive favours VAPOR):\n", REF_METHOD))
   print(pd, row.names=FALSE)
   all_paired[[length(all_paired)+1]] <- pd
   
   for (mm in methods)
      all_rows[[length(all_rows)+1]] <- data.frame(
         p=p, T=Tval, TP_ratio=round(Tval/p,2), block=s$block, Method=mm,
         Omega_Frob=agg[mm,"Frob"],     Omega_Frob_sd=agg_sd[mm,"Frob"],
         B_Frob=agg[mm,"B_Frob"],       B_Frob_sd=agg_sd[mm,"B_Frob"],
         Sens=agg[mm,"Sens"],           Sens_sd=agg_sd[mm,"Sens"],
         Spec=agg[mm,"Spec"],           Spec_sd=agg_sd[mm,"Spec"])
   
   saveRDS(list(results=do.call(rbind, all_rows),
                paired=do.call(rbind, all_paired), raw=raw_G),
           file.path(OUTDIR, "sim_g_partial.rds"))
}

results <- do.call(rbind, all_rows)
paired  <- do.call(rbind, all_paired)

cat("\n============================================================\n")
cat("FULL RESULTS TABLE\n")
cat("============================================================\n")
print(results, row.names=FALSE, digits=3)
cat("\n--- Paired Omega-Frobenius differences vs", REF_METHOD, "---\n")
print(paired, row.names=FALSE)

saveRDS(list(results=results, paired=paired, raw=raw_G, settings=SETTINGS,
             n_seeds=N_SEEDS, FIT_RANKS=FIT_RANKS, r0=R0, k=K_FACTORS),
        file.path(OUTDIR, "sim_g_results.rds"))
cat(sprintf("\nSaved results + %d heatmaps to %s\n", length(SETTINGS), OUTDIR))