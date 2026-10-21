# ============================================================================
#  Items 2-4: corrected posterior predictive check, extended edge recovery and
#  group purity, and the within-group audit of VAPOR's edges.
#     Rscript applications/recheck_ppc_recovery.R
#
#  Reads saved objects only (p29_final.rds, p100_final.rds, and
#  kastner_p100_retry.rds when it exists).  Fits nothing.  Appends its tables
#  to results_summary.md and saves recheck_tables.rds.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 205)
APP <- FITS_DIR
MD  <- file.path(TAB_DIR, "results_summary.md")
CFG <- list(PPC_N_REP = 400L, BIG = 1e3, GROUPS = file.path(dirname(RES_DIR), "fredqd_groups.csv"),
            TOPS = list(`29` = c(7L, 28L, 70L), `100` = c(25L, 100L, 250L)),
            PPC_SEED = c(`29` = 7029L, `100` = 7001L))
OUT <- list()

# ============================================================================
#  2a. the nine summaries: definition now in force, against the paper
# ============================================================================
defn <- data.frame(
   statistic = c("max_eig","eig2","mean_abs_cor","q90_abs_cor","max_abs_pcor",
                 "prop_pcor_02","mean_sd","kurtosis","resid_ac1"),
   paper = c("largest eigenvalue of cor(y_rep)",
             "second eigenvalue of cor(y_rep)",
             "mean |off-diagonal cor(y_rep)|",
             "90th percentile of |off-diagonal cor(y_rep)|",
             "max |partial cor| from solve(cov(y_rep))",
             "share of |partial cor| below 0.2",
             "mean over series of sd(y_rep)",
             "mean over series of the 4th standardized moment of y_rep",
             "mean over series of lag-1 autocorrelation of y_rep - B_s y_{t-1}"),
   now = c("largest eigenvalue of cor(y_rep)",
           "second eigenvalue of cor(y_rep)",
           "mean |off-diagonal cor(y_rep)|",
           "90th percentile of |off-diagonal cor(y_rep)|",
           "max |partial cor| from solve(cov(y_rep))",
           "share of |partial cor| below 0.2",
           "mean over series of sd(y_rep)",
           "mean over series of the 4th standardized moment of y_rep",
           "mean over series of lag-1 autocorrelation of y_rep - B_s y_{t-1}"),
   changed = c(rep("panel definition changed", 8), "yes, and panel definition changed"),
   stringsAsFactors = FALSE)
scheme <- data.frame(
   element = c("replication", "panel rows", "resid_ac1 input"),
   paper = c("one step ahead: y_rep_t ~ N(B_s y_{t-1}, Omega_s^-1), observed lag",
             "t = 2..T, so T-1 rows",
             "y_rep - B_s y_{t-1}"),
   before = c("recursive: y_rep_t = B_s y_rep_{t-1} + e_t, simulated lag",
              "t = 1..T with row 1 set to the observed row",
              "y_rep itself"),
   now = c("one step ahead, observed lag (matches the paper)",
           "t = 2..T (matches the paper)",
           "y_rep - B_s y_{t-1} (matches the paper)"),
   stringsAsFactors = FALSE)
cat("=== 2a. replication scheme ===\n"); print(scheme, row.names = FALSE)
cat("\n=== 2a. the nine summaries ===\n"); print(defn, row.names = FALSE)
OUT$ppc_scheme <- scheme; OUT$ppc_defn <- defn

for (d in c(29L, 100L)){
   f <- file.path(APP, sprintf("p%d_final.rds", d)); if (!file.exists(f)) next
   R <- readRDS(f); ch <- R$draws$vapor; Y <- R$Y_train; p <- R$p
   key <- as.character(d)
   z <- ppc_table(ch, Y, p, CFG$PPC_N_REP, CFG$PPC_SEED[[key]],
                  BIG = CFG$BIG, return_detail = TRUE)
   # -- 2b -----------------------------------------------------------------
   det <- z$detail; bad <- z$degenerate
   cat(sprintf("\n=== 2b. p = %d: draws whose replicated panel has any |value| > %.0e ===\n",
               d, CFG$BIG))
   b2 <- data.frame(n_draws_sampled = length(det), n_degenerate = sum(bad),
                    share = round(mean(bad), 3))
   print(b2, row.names = FALSE)
   if (any(bad)){
      dd <- do.call(rbind, lapply(det[bad], function(q) data.frame(
         chain = q$chain, draw = q$draw, max_abs_yrep = signif(q$max_abs, 4),
         n_resid_prec_at_floor = q$at_floor,
         condition_number_Omega = signif(q$kappa, 4))))
      dd <- dd[order(-dd$max_abs_yrep), ]
      print(head(dd, 25), row.names = FALSE)
      cat(sprintf("of the %d degenerate draws, %d have at least one residual precision at the 1e-8 floor\n",
                  sum(bad), sum(dd$n_resid_prec_at_floor > 0)))
      cat(sprintf("condition number of Omega_s: median %.4g over degenerate draws, %.4g over the rest\n",
                  median(dd$condition_number_Omega, na.rm=TRUE),
                  median(vapply(det[!bad], function(q) q$kappa, numeric(1)), na.rm=TRUE)))
      OUT[[sprintf("ppc_degenerate_p%d", d)]] <- dd
   } else cat("none\n")
   OUT[[sprintf("ppc_degen_count_p%d", d)]] <- b2
   # -- 2c -----------------------------------------------------------------
   t_all <- z$table
   cat(sprintf("\n=== 2c. p = %d: posterior predictive check, all %d draws ===\n",
               d, attr(t_all, "n_draws")))
   print(t_all, row.names = FALSE)
   OUT[[sprintf("ppc_all_p%d", d)]] <- t_all
   if (any(bad)){
      t_ex <- ppc_table(ch, Y, p, CFG$PPC_N_REP, CFG$PPC_SEED[[key]],
                        BIG = CFG$BIG, exclude_degenerate = TRUE)
      cat(sprintf("\n=== 2c. p = %d: excluding the %d degenerate draws (%d kept) ===\n",
                  d, sum(bad), attr(t_ex, "n_draws")))
      print(t_ex, row.names = FALSE)
      OUT[[sprintf("ppc_clean_p%d", d)]] <- t_ex }

   # ========================================================================
   #  3. edge recovery and group purity, with every available method
   # ========================================================================
   vars <- R$vars
   mean_pc <- function(oms) if (!length(oms)) NULL else
      Reduce(`+`, lapply(oms, pcor_from_omega))/length(oms)
   #  item 4: ONE convention.  Every Bayesian method is the mean of its
   #  per-draw partial correlations pooled over all chains; the two
   #  deterministic methods keep their point estimates.  At p = 100 Kastner's
   #  only available chains are the 20,000-iteration run, which stored no
   #  per-draw precisions, so that one column stays on the mean-precision
   #  estimate and the provenance table below says so.
   p <- length(vars)
   PCs <- list(VAPOR = pcor_mean_draws(R$draws$vapor, p),
               `PrecFact (no VAR)` = pcor_mean_draws(R$draws$precfact, p),
               `VAR + PrecFact`    = pcor_mean_draws(R$draws$var_precfact, p),
               `VAR + Glasso`      = mean_pc(R$omegas[["VAR + Glasso"]]),
               `Sample covariance` = mean_pc(R$omegas[["Sample covariance"]]),
               Kastner             = pcor_mean_prec(R$draws$kastner, p))
   if (d == 100L){
      kf <- file.path(APP, "kastner_p100_retry.rds")
      if (file.exists(kf)){
         KR <- readRDS(kf)
         if (length(KR$chains)){
            PCs$Kastner <- mean_pc(lapply(KR$chains, function(q) q$Omega))   # no per-draw precisions stored
            cat(sprintf("\nKastner p = 100: %d retry chain(s) (seeds %s) used for the comparisons\n",
                        length(KR$chains),
                        paste(vapply(KR$chains, function(q) q$seed, numeric(1)), collapse=", "))) } } }
   if (is.null(PCs$Kastner)){
      kf <- file.path(APP, "kastner_p100_20k.rds")
      if (file.exists(kf)){ KK <- readRDS(kf)
         if (length(KK$chains)) PCs$Kastner <- mean_pc(lapply(KK$chains, function(q) q$Omega)) } }
   PCs <- PCs[!vapply(PCs, is.null, logical(1))]
   cat("\n  partial-correlation convention:\n")
   print(do.call(rbind, lapply(names(PCs), function(nm){ z <- PCs[[nm]]
      data.frame(method = nm, draws = if (is.null(attr(z,"n_draws"))) NA_integer_ else attr(z,"n_draws"),
                 estimate = if (!is.null(attr(z,"n_draws"))) "mean of per-draw pcor" else
                            if (nm == "Kastner") "pcor of each chain's mean precision" else "point estimate") })),
      row.names = FALSE)
   ed <- R$long$ci_edges
   rk <- edge_recovery_ranks(ed, PCs, vars)
   rownames(rk) <- paste(ed$v1, ed$v2, sep="--")
   tops <- CFG$TOPS[[key]]
   summ <- do.call(rbind, lapply(colnames(rk), function(m){
      v <- rk[,m]
      setNames(data.frame(m, round(median(v),1),
                          sum(v <= tops[1]), sum(v <= tops[2]), sum(v <= tops[3])),
               c("method","median_rank",
                 sprintf("in_top_%d", tops[1]), sprintf("in_top_%d", tops[2]),
                 sprintf("in_top_%d", tops[3]))) }))
   cat(sprintf("\n=== 3. p = %d: rank of each of VAPOR's %d edges in each method's ordering ===\n",
               d, nrow(rk)))
   print(rk)
   cat(sprintf("\n=== 3. p = %d: summary over the %d edges ===\n", d, nrow(rk)))
   print(summ, row.names = FALSE)
   OUT[[sprintf("recovery_p%d", d)]] <- rk; OUT[[sprintf("recovery_summary_p%d", d)]] <- summ

   if (d == 100L){
      pur <- group_purity(PCs, vars, CFG$GROUPS, 100L)
      cat("\n=== 3. p = 100: share of the top 100 pairs within one official FRED-QD group ===\n")
      print(pur, row.names = FALSE)
      OUT$purity_p100 <- pur

      # ====================================================================
      #  4. are all of VAPOR's edges within a single group?
      # ====================================================================
      g <- read.csv(CFG$GROUPS, stringsAsFactors = FALSE)
      gc <- setNames(g$group_code, g$series); gn <- setNames(g$group_name, g$series)
      aud <- data.frame(edge = paste(ed$v1, ed$v2, sep="--"), pcor = ed$pcor,
                        group_1 = gn[ed$v1], group_2 = gn[ed$v2],
                        same = gc[ed$v1] == gc[ed$v2], row.names = NULL)
      cat(sprintf("\n=== 4. p = 100: %d VAPOR edges, %d within one group, %d crossing ===\n",
                  nrow(aud), sum(aud$same), sum(!aud$same)))
      print(aud[!aud$same, c("edge","pcor","group_1","group_2")], row.names = FALSE)
      OUT$edge_group_audit <- aud } }

# ============================================================================
#  append to results_summary.md
# ============================================================================
cat("\n\n", file = MD, append = TRUE)
md_section(MD, "PPC replication scheme, paper vs before vs now", OUT$ppc_scheme)
md_section(MD, "PPC nine summaries, paper vs now", OUT$ppc_defn)
for (d in c(29,100)){
   md_section(MD, sprintf("p = %d PPC degenerate-draw count", d), OUT[[sprintf("ppc_degen_count_p%d", d)]])
   md_section(MD, sprintf("p = %d PPC degenerate draws", d), OUT[[sprintf("ppc_degenerate_p%d", d)]])
   md_section(MD, sprintf("p = %d PPC (corrected, all draws)", d), OUT[[sprintf("ppc_all_p%d", d)]])
   md_section(MD, sprintf("p = %d PPC (corrected, degenerate draws excluded)", d), OUT[[sprintf("ppc_clean_p%d", d)]])
   md_section(MD, sprintf("p = %d edge recovery ranks (all methods)", d), OUT[[sprintf("recovery_p%d", d)]])
   md_section(MD, sprintf("p = %d edge recovery summary", d), OUT[[sprintf("recovery_summary_p%d", d)]]) }
md_section(MD, "p = 100 group purity (all methods)", OUT$purity_p100)
md_section(MD, "p = 100 VAPOR edges, official FRED-QD group of each endpoint", OUT$edge_group_audit)
saveRDS(OUT, file.path(APP, "recheck_tables.rds"))
cat(sprintf("\nappended to %s and saved recheck_tables.rds\n", MD))
