# ============================================================================
#  p = 29 rank and k table (Table S5).  Rscript applications/final_ranks_p29.R
#
#  Cells: k in {5, 8, 10, 15} x rank in {full, 2, 3, 4}, 3 seeds each, 100,000
#  iterations.  Columns are the ones Table S5 reports, computed exactly as the
#  superseded ranks_p29.R computed them: Pr(spectral radius of B > 1), the
#  median spectral radius, the posterior predictive replicated maximum
#  eigenvalue of the correlation matrix and its p-value, held-out MSPE, and the
#  number of edges under the 95% credible-interval rule.  Reported as the mean
#  over the 3 seeds, with the SD alongside.
#
#  Cells run concurrently, up to physical cores minus one.  Every cell's Omega
#  is saved.  Writes ranks_p29_final.rds.  Every setting is in CFG below.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)

CFG <- list(
   LABEL = "p = 29 rank and k table",
   P = 29L,
   CSV = "/Users/peterdunson/Desktop/vapor/applications/data/2026-05-qd.csv",
   OUT = FITS_DIR,
   OUT_RDS = "ranks_p29_final.rds",

   K_GRID    = c(5L, 8L, 10L, 15L),
   RANK_GRID = list("full", 2L, 3L, 4L),
   SEEDS = 101:103,                      # 3 seeds per cell
   NRUN = 100000L, BURN = 20000L, THIN = 80L,
   JIT  = 5e-3,                          # the config D jitter
   VAR_JIT = 1e-3, WORKERS = NULL,

   PR_EXPL_DRAWS = 400L,                 # draws used for Pr(rho(B) > 1)
   PPC_DRAWS = 200L, PPC_SEED = 1L,      # draws used for the eigenvalue PPC

   N_HID = 100L, HID_SEED_BASE = 90000L, HID_FRAC = 0.3,
   #  corrected panel, matching final_p29.R
   WIN   = c(1996 + 0/12, 2025 + 11/12),
   TRAIN = c(1996 + 0/12, 2020 + 11/12),
   TEST  = c(2021 + 2/12, 2025 + 11/12),
   VARS  = c("GDPC1","INDPRO","IPFINAL","IPCONGD","IPBUSEQ","PAYEMS","UEMPMEAN",
      "CLAIMSx","CUMFNS","PCECC96","CMRMTSPLx","BUSINVx","CPIAUCSL","PPICMM","GDPCTPI",
      "OILPRICEx","FEDFUNDS","GS1","GS10","BAA","TB3MS","M2REAL","BUSLOANSx","NONREVSLx",
      "S.P.500","VIXCLSx","GPDIC1","UNRATE","WPSFD49207"))

D <- build_panel(CFG)
Y_train <- D$train; Y_test <- D$test; vars <- D$vars; p <- D$p
hid_sets <- make_hid_sets(CFG$HID_SEED_BASE, CFG$N_HID, p, nrow(Y_test), CFG$HID_FRAC)

one_fit <- function(k, rank, seed){
   f <- NULL; jj <- CFG$JIT; a <- 0L
   while (is.null(f) && a < 4L){ a <- a + 1L
      f <- tryCatch(vapor(Y_train, lags=1, k=k, rank=rank, prior_B="dl",
                          nrun=CFG$NRUN, burn=CFG$BURN, thin=CFG$THIN,
                          seed=seed, jitter=jj, init="random"),
                    error = function(e) NULL)
      jj <- jj*2 }
   if (!is.null(f)) attr(f, "attempts") <- a
   f }

CELLS <- do.call(rbind, lapply(CFG$RANK_GRID, function(r)
   do.call(rbind, lapply(CFG$K_GRID, function(k)
      do.call(rbind, lapply(CFG$SEEDS, function(s)
         data.frame(rank_label = if (identical(r,"full")) "Full" else as.character(r),
                    k = k, seed = s, stringsAsFactors = FALSE)))))))
rank_of <- function(lab) if (identical(lab, "Full")) "full" else as.integer(lab)
W <- n_workers(CFG$WORKERS)
cat(sprintf("%d cells (%d k x %d rank x %d seeds) at %d iterations, %d worker(s)\n",
            nrow(CELLS), length(CFG$K_GRID), length(CFG$RANK_GRID),
            length(CFG$SEEDS), CFG$NRUN, W))

RES <- parallel::mclapply(seq_len(nrow(CELLS)), function(i){
   c_ <- CELLS[i,]; t0 <- Sys.time()
   f <- one_fit(c_$k, rank_of(c_$rank_label), c_$seed)
   el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
   if (is.null(f)) return(list(rank_label=c_$rank_label, k=c_$k, seed=c_$seed,
                               status="failed", minutes=el))
   Om <- omega_hat(f); M <- pair_draws(f, p); sel <- ci_indices(M)
   re <- pr_explosive(f, p, CFG$PR_EXPL_DRAWS)
   pe <- ppc_maxeig(f, Y_train, p, CFG$PPC_DRAWS, CFG$PPC_SEED)
   list(rank_label=c_$rank_label, k=c_$k, seed=c_$seed, status="completed",
        minutes=el, attempts=attr(f,"attempts"), Omega=Om,
        pr_expl=unname(re["pr"]), med_rho=unname(re["med"]),
        rep_eig=unname(pe["rep"]), ppp=unname(pe["ppp"]), obs_eig=unname(pe["obs"]),
        mspe=mean(loss_avg(Om, hid_sets, Y_test, p)),
        n_edges=length(sel), edges=pair_labels(vars)[sel],
        floor_hits=f$floor_hits)
}, mc.cores = W, mc.preschedule = FALSE)

ok <- Filter(function(z) !inherits(z,"try-error") && !is.null(z) &&
                         identical(z$status,"completed"), RES)
long <- do.call(rbind, lapply(ok, function(z) data.frame(
   rank=z$rank_label, k=z$k, seed=z$seed, pr_expl=z$pr_expl, med_rho=z$med_rho,
   rep_eig=z$rep_eig, ppp=z$ppp, mspe=z$mspe, edges=z$n_edges,
   minutes=z$minutes, series_floored=sum(z$floor_hits > 0))))
cat(sprintf("\n%d of %d cells completed | observed max eigenvalue = %.2f\n",
            length(ok), nrow(CELLS), mean(vapply(ok, function(z) z$obs_eig, numeric(1)))))

ms <- function(v, d) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(v), sd(v))
tab <- do.call(rbind, lapply(split(long, list(long$rank, long$k), drop=TRUE),
   function(d) data.frame(rank=d$rank[1], k=d$k[1], n_seeds=nrow(d),
      `Pr_rho_gt_1` = ms(d$pr_expl, 3),
      `median_rho`  = ms(d$med_rho, 3),
      `PPC_max_eig` = ms(d$rep_eig, 2),
      `PPC_p`       = ms(d$ppp, 3),
      MSPE          = ms(d$mspe, 4),
      edges         = ms(d$edges, 1))))
tab$rank <- factor(tab$rank, levels = c("Full","2","3","4"))
tab <- tab[order(tab$rank, tab$k), ]
cat("\n=== Table S5: p = 29 rank and k, mean (SD) over 3 seeds ===\n")
print(tab, row.names=FALSE)

cat("\n--- LaTeX rows ---\n")
agg <- aggregate(cbind(pr_expl, med_rho, rep_eig, ppp, mspe, edges) ~ rank + k,
                 data = long, FUN = mean)
agg$rank <- factor(agg$rank, levels = c("Full","2","3","4"))
agg <- agg[order(agg$rank, agg$k), ]
for (i in seq_len(nrow(agg))) with(agg[i,],
   cat(sprintf("%s & %d & %.3f & %.3f & %.2f & %.3f & %.4f & %.1f \\\\\n",
               rank, k, pr_expl, med_rho, rep_eig, ppp, mspe, edges)))

saveRDS(list(config=CFG, cells=CELLS, results=RES, long=long, table=tab, agg=agg,
             run_at=format(Sys.time())), file.path(CFG$OUT, CFG$OUT_RDS))
cat(sprintf("\nsaved %s\n", CFG$OUT_RDS))
