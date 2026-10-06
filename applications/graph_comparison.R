# ============================================================================
#  Graph comparison across methods, p = 29 and p = 100.
#  Reads only the saved results files; refits nothing.
#
#  "VAR + Glasso" and "VAR + PrecFact" are the RANK-3 fits. The full-rank
#  versions are saved in the results file but excluded from R$order, so they
#  appear in no table here.
#
#  MATCHED MCMC EFFORT. Every stochastic method -- VAPOR and all three
#  stochastic competitors -- is represented by partial correlations POOLED
#  over 4 long chains (20,000 iterations, 10,000 burn-in, thin 10), read from
#  p<tag>_long_competitors.rds. This replaces the earlier 3-production-fit
#  averaging, so no method gets a sampling-effort advantage. Deterministic
#  methods (glasso, sample covariance) have a single fit. If the long-chain
#  file is missing the script falls back to the 3-fit average and says so.
#
#  SECTOR PURITY is gated behind RUN_SECTOR. It needs official FRED-QD group
#  codes, which are NOT in the CSV (that file carries only `factors` and
#  `transform` rows), so they come from applications/fredqd_groups.csv, built
#  separately from the McCracken-Ng appendix. Leave RUN_SECTOR FALSE until
#  that file has been checked.
# ============================================================================
APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")

RUN_SECTOR <- TRUE
GROUPS_CSV <- "applications/fredqd_groups.csv"
THRESH     <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

SPEC <- list(
   p29  = list(res = "p29_results.rds",  out = "graph_comparison_p29.rds",
               topK = c(7, 25, 50),   chains = "p29_chains_final.rds",
               long = "p29_long_competitors.rds",
               rates = c("FEDFUNDS","TB3MS","GS1","GS5","GS10")),
   p100 = list(res = "p100_results.rds", out = "graph_comparison_p100.rds",
               topK = c(25, 50, 100), chains = "p100_chains_final.rds",
               long = "p100_long_competitors.rds",
               rates = c("FEDFUNDS","TB3MS","GS1","GS5","GS10")))

# ---- helpers ---------------------------------------------------------------
ut_index <- function(p) which(upper.tri(matrix(0,p,p)), arr.ind = TRUE)
pair_lab <- function(i, j, vars) paste(vars[pmin(i,j)], vars[pmax(i,j)], sep="--")

#  rank every pair by |partial correlation|, 1 = strongest
rank_vector <- function(PC, ut){
   v <- abs(PC[cbind(ut[,1], ut[,2])])
   o <- order(v, decreasing = TRUE); rk <- integer(length(v)); rk[o] <- seq_along(v)
   rk }

#  posterior 95% interval for a handful of pairs, from the saved long chains
ci_for_pairs <- function(chains_path, ij){
   if (is.null(chains_path) || !file.exists(chains_path) || !nrow(ij)) return(NULL)
   ch <- readRDS(chains_path)
   if (!is.null(ch$chains)) ch <- ch$chains      # {chains, config} shape
   acc <- matrix(NA_real_, 0, nrow(ij))
   for (f in ch){
      nd <- dim(f$Lambda)[3]
      blk <- matrix(NA_real_, nd, nrow(ij))
      for (s in seq_len(nd)){
         Om <- tcrossprod(f$Lambda[,,s]); diag(Om) <- diag(Om) + f$residual_prec[s,]
         d <- sqrt(diag(Om)); pc <- -Om/outer(d,d)
         blk[s,] <- pc[cbind(ij[,1], ij[,2])] }
      acc <- rbind(acc, blk) }
   rm(ch); gc(verbose = FALSE)
   data.frame(mean = colMeans(acc),
              lo = apply(acc, 2, quantile, 0.025),
              hi = apply(acc, 2, quantile, 0.975),
              n_draws = nrow(acc)) }

analyse <- function(tag, spec){
   f <- file.path(APP, spec$res)
   if (!file.exists(f)){ cat(sprintf("\n[%s] %s not found, skipped\n", tag, spec$res)); return(NULL) }
   R <- readRDS(f)
   vars <- R$vars; p <- R$p; ut <- ut_index(p)
   labs <- pair_lab(ut[,1], ut[,2], vars)
   cat(sprintf("\n\n################  %s  (p = %d, %d pairs)  ################\n",
               tag, p, nrow(ut)))

   # ---- one partial-correlation matrix per method --------------------------
   LC <- NULL
   lcp <- file.path(APP, spec$long)
   if (file.exists(lcp)) LC <- readRDS(lcp) else
      cat(sprintf("\n[%s is missing: falling back to 3-fit averaging]\n", spec$long))

   PC <- list(); SRC <- character(0)
   for (k in R$order){
      m <- R$methods[[k]]
      if (!length(m$pcor)) next               # random walk has no Omega
      if (!is.null(LC) && !is.null(LC$pooled_pcor[[k]])){
         PC[[m$label]] <- LC$pooled_pcor[[k]]
         SRC[m$label] <- sprintf("pooled 4 long chains (%d draws)",
                                 LC$diagnostics$pooled_draws[LC$diagnostics$Method == k])
      } else if (length(m$pcor) > 1){
         PC[[m$label]] <- Reduce(`+`, m$pcor)/length(m$pcor)
         SRC[m$label] <- sprintf("mean of %d production fits", length(m$pcor))
      } else {
         PC[[m$label]] <- m$pcor[[1]]
         SRC[m$label] <- "single fit (deterministic)" } }

   #  VAPOR is the pooled long-chain mean, which is what pc_vapor already is
   vnm <- "VAPOR (pooled long chains)"
   PC[[vnm]] <- R$long$pc_vapor
   SRC[vnm] <- "pooled 4 long chains (pc_vapor)"
   if (!is.null(LC$pooled_pcor[["VAPOR (rank 3)"]])){
      d <- max(abs(LC$pooled_pcor[["VAPOR (rank 3)"]] - R$long$pc_vapor))
      cat(sprintf("\ncheck: long-chain VAPOR pooling vs stored pc_vapor, max abs diff %.2e\n", d)) }
   PC <- PC[!grepl("^VAPOR \\(rank", names(PC))]
   SRC <- SRC[names(SRC) %in% names(PC)]
   cat("\nmatrix used per method:\n")
   print(data.frame(Method = names(SRC), matrix = unname(SRC), row.names = NULL),
         row.names = FALSE)
   RK <- lapply(PC, rank_vector, ut)

   ci_edges <- R$long$ci_edges
   vi <- cbind(match(ci_edges$v1, vars), match(ci_edges$v2, vars))
   vkeys <- pair_lab(vi[,1], vi[,2], vars)
   nE <- length(vkeys)
   cat(sprintf("\nVAPOR reports %d edges under the 95%% credible-interval rule\n", nE))

   # ========================================================================
   #  EDGE COUNTS BY THRESHOLD
   # ========================================================================
   edge_tab <- data.frame(threshold = THRESH)
   for (nm in names(PC))
      edge_tab[[nm]] <- sapply(THRESH, function(t)
         sum(abs(PC[[nm]][upper.tri(PC[[nm]])]) > t, na.rm = TRUE))
   cat("\n=== |partial correlation| > threshold ===\n")
   print(edge_tab, row.names = FALSE)

   # ========================================================================
   #  EDGE RECOVERY: where do VAPOR's CI edges rank under each method?
   # ========================================================================
   K1 <- nE; K4 <- 4*nE; K10 <- 10*nE
   rec <- do.call(rbind, lapply(names(RK), function(nm){
      r <- RK[[nm]][match(vkeys, labs)]
      data.frame(Method = nm, median_rank = median(r),
                 top_K = sum(r <= K1), top_4K = sum(r <= K4), top_10K = sum(r <= K10)) }))
   cat(sprintf("\n=== rank of VAPOR's %d CI edges (K=%d, 4K=%d, 10K=%d) ===\n",
               nE, K1, K4, K10))
   print(rec, row.names = FALSE)

   # ========================================================================
   #  OVERLAP WITH RANK-3 GLASSO
   # ========================================================================
   gl3 <- grep("Glasso", names(RK), value = TRUE)   # the rank-3 fit
   overlap <- NULL; missed <- NULL
   if (length(gl3)){
      r3 <- RK[[gl3[1]]]
      inK  <- sum(r3[match(vkeys, labs)] <= nE)
      in2K <- sum(r3[match(vkeys, labs)] <= 2*nE)
      overlap <- data.frame(vapor_edges = nE, in_glasso_topK = inK, in_glasso_top2K = in2K)
      cat(sprintf("\n=== VAPOR CI edges inside rank-3 glasso's top K / top 2K ===\n"))
      print(overlap, row.names = FALSE)

      #  glasso's top K that VAPOR does not select
      topK_lab <- labs[order(r3)][seq_len(nE)]
      extra <- setdiff(topK_lab, vkeys)
      if (length(extra)){
         ij <- do.call(rbind, lapply(extra, function(L){
            s <- strsplit(L, "--", fixed = TRUE)[[1]]
            c(match(s[1], vars), match(s[2], vars)) }))
         vp <- R$long$pc_vapor[cbind(ij[,1], ij[,2])]
         ci <- ci_for_pairs(if (is.null(spec$chains)) NULL else file.path(APP, spec$chains), ij)
         missed <- data.frame(pair = extra,
                              glasso_pcor = round(PC[[gl3[1]]][cbind(ij[,1], ij[,2])], 3),
                              vapor_mean  = round(vp, 3),
                              vapor_lo = if (is.null(ci)) NA else round(ci$lo, 3),
                              vapor_hi = if (is.null(ci)) NA else round(ci$hi, 3))
         cat(sprintf("\n=== rank-3 glasso's top %d that VAPOR does not select ===\n", nE))
         if (is.null(ci))
            cat("   (VAPOR 95% interval needs the saved long chains; none for this dimension)\n")
         print(missed, row.names = FALSE) } }

   # ========================================================================
   #  SIGNS OF VAPOR'S NEGATIVE EDGES UNDER EVERY METHOD
   # ========================================================================
   negi <- which(ci_edges$pcor < 0)
   signs <- NULL
   if (length(negi)){
      ijn <- vi[negi, , drop = FALSE]
      signs <- data.frame(pair = vkeys[negi], vapor_pcor = ci_edges$pcor[negi])
      for (nm in names(PC)){
         signs[[paste0(nm, " sign")]] <- sign(PC[[nm]][cbind(ijn[,1], ijn[,2])])
         signs[[paste0(nm, " rank")]] <- RK[[nm]][match(vkeys[negi], labs)] }
      cat("\n=== VAPOR's negative edges: sign and rank under each method ===\n")
      print(signs, row.names = FALSE)
   } else cat("\n(VAPOR reports no negative edges at this dimension)\n")

   # ========================================================================
   #  TERM STRUCTURE
   # ========================================================================
   rates <- spec$rates[spec$rates %in% vars]
   ts_tab <- NULL
   if (length(rates) > 1){
      ri <- match(rates, vars)
      cmb <- t(combn(seq_along(rates), 2))
      ts_tab <- data.frame(
         pair = paste(rates[cmb[,1]], rates[cmb[,2]], sep="--"),
         gap  = cmb[,2] - cmb[,1],
         adjacent = (cmb[,2] - cmb[,1]) == 1)
      for (nm in names(PC))
         ts_tab[[nm]] <- round(PC[[nm]][cbind(ri[cmb[,1]], ri[cmb[,2]])], 3)
      cat(sprintf("\n=== term structure: %s ===\n", paste(rates, collapse=", ")))
      print(ts_tab, row.names = FALSE)
      for (nm in grep("Glasso", names(PC), value = TRUE)){
         z <- ts_tab[[nm]][!ts_tab$adjacent]
         cat(sprintf("  %-28s non-adjacent entries exactly zero: %d of %d\n",
                     nm, sum(z == 0), length(z))) } }

   # ========================================================================
   #  GLASSO'S OWN GRAPH
   # ========================================================================
   gl_tab <- do.call(rbind, lapply(grep("Glasso", R$order, value = TRUE), function(k){
      m <- R$methods[[k]]
      data.frame(Method = m$label, rho = round(m$extra$rho, 4),
                 nonzero_offdiag = m$extra$n_nonzero) }))
   cat("\n=== glasso's own graph at the BIC penalty ===\n")
   print(gl_tab, row.names = FALSE)

   # ========================================================================
   #  SECTOR PURITY  (gated: needs the official FRED-QD groups)
   # ========================================================================
   purity <- cross_pairs <- NULL
   if (RUN_SECTOR){
      if (!file.exists(GROUPS_CSV)) stop("RUN_SECTOR is TRUE but ", GROUPS_CSV, " is missing")
      g <- read.csv(GROUPS_CSV, stringsAsFactors = FALSE)
      grp <- setNames(g$group_code, g$series)
      miss <- setdiff(vars, names(grp))
      if (length(miss)) stop("no group for: ", paste(miss, collapse=", "))
      gi <- grp[vars]
      same <- gi[ut[,1]] == gi[ut[,2]]
      purity <- NULL; cross_pairs <- list()
      for (nm in names(RK)) for (K in spec$topK){
         sel <- order(RK[[nm]])[seq_len(K)]
         purity <- rbind(purity, data.frame(Method = nm, topK = K,
                                            within_group = sum(same[sel]),
                                            fraction = round(mean(same[sel]), 3)))
         cross_pairs[[paste(nm, K)]] <- labs[sel][!same[sel]] }
      #  purity of VAPOR's own CI edge set, not just its top-K by magnitude
      vloc <- match(vkeys, labs)
      ci_purity <- data.frame(n_edges = nE, within_group = sum(same[vloc]),
                              fraction = round(mean(same[vloc]), 3))
      cat(sprintf("\n=== group purity of VAPOR's %d CI edges themselves ===\n", nE))
      print(ci_purity, row.names = FALSE)
      xc <- which(!same[vloc])
      if (length(xc)){
         gn <- setNames(g$group_name, g$series)
         cross_ci <- data.frame(
            pair = vkeys[xc],
            pcor = ci_edges$pcor[xc],
            group_1 = gn[vars[vi[xc,1]]], group_2 = gn[vars[vi[xc,2]]])
         cat("cross-group CI edges:\n"); print(cross_ci, row.names = FALSE)
      } else { cross_ci <- NULL; cat("cross-group CI edges: none\n") }

      cat("\n=== sector purity: fraction of top-K pairs within one FRED-QD group ===\n")
      print(purity, row.names = FALSE)
      cat("\ncross-group pairs:\n")
      for (nm in names(cross_pairs))
         cat(sprintf("  %-36s %s\n", nm, paste(cross_pairs[[nm]], collapse=", ")))
      #  within-group fraction of glasso's own selected support
      for (k in grep("Glasso", R$order, value = TRUE)){
         E <- R$methods[[k]]$extra$nonzero
         sel <- E[cbind(ut[,1], ut[,2])]
         cat(sprintf("  %-28s support within group: %.3f (%d of %d)\n",
                     R$methods[[k]]$label, mean(same[sel]), sum(same[sel]), sum(sel))) }
   } else {
      ci_purity <- cross_ci <- NULL
      cat("\n[sector purity skipped: RUN_SECTOR is FALSE]\n") }

   # ========================================================================
   #  NAMED PAIRS: value and rank under every method
   # ========================================================================
   negs <- vkeys[ci_edges$pcor < 0]
   want <- intersect(c("GS5--GS10", "TB3MS--GS1", "GS5--MORTGAGE30US", negs), labs)
   named <- NULL
   if (length(want)){
      wi <- match(want, labs)
      named <- do.call(rbind, lapply(names(PC), function(nm)
         data.frame(Method = nm, pair = want,
                    pcor = round(PC[[nm]][cbind(ut[wi,1], ut[wi,2])], 3),
                    rank = RK[[nm]][wi])))
      named <- named[order(match(named$pair, want), named$Method), ]
      cat(sprintf("\n=== named pairs: value and rank of %d ===\n", nrow(ut)))
      print(named, row.names = FALSE) }

   # ========================================================================
   #  R-HAT / ESS on the matched long chains
   # ========================================================================
   if (!is.null(LC)){
      cat("\n=== R-hat / ESS on the monitored log-likelihood (4 long chains) ===\n")
      cat("  ", LC$loglik_formula, "\n\n", sep = "")
      print(LC$diagnostics, row.names = FALSE) }

   res <- list(tag = tag, p = p, vars = vars, n_vapor_edges = nE,
               named_pairs = named,
               matrix_source = SRC,
               long_diagnostics = if (is.null(LC)) NULL else LC$diagnostics,
               edge_counts = edge_tab, edge_recovery = rec,
               overlap = overlap, glasso_topK_missed_by_vapor = missed,
               signs = signs, term_structure = ts_tab, glasso_graph = gl_tab,
               sector_purity = purity, cross_group_pairs = cross_pairs,
               ci_edge_purity = ci_purity, cross_group_ci_edges = cross_ci,
               ran_sector = RUN_SECTOR, source = spec$res,
               run_at = format(Sys.time()))
   saveRDS(res, file.path(APP, spec$out))
   cat(sprintf("\nSaved %s\n", spec$out))
   invisible(res) }

out <- lapply(names(SPEC), function(tag) analyse(tag, SPEC[[tag]]))
invisible(out)
