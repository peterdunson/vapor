# ============================================================================
#  Checks for numbers stated in the text of the paper, reproduced from the
#  in-repo data and the installed package alone.
#     Rscript applications/paper_checks.R
#  Writes results/tables/paper_checks.txt.
#
#  WHY THIS SCRIPT EXISTS.  Three statements in the text had no producer in the
#  repository; they were recorded only in loose notes under results/checks/.
#  Each is reproduced here so the claim has a script behind it:
#
#    (a) "a third of its mass below 1e-4", of the Gamma(0.1, 0.1) base measure
#    (b) "R^2 = 0.997" and "smallest eigenvalue 0.0013", the NIPA identity in
#        the p = 100 training panel.  These ALSO appear in results_summary.md,
#        written by final_nipa_p100.R, which reads p100_results.rds; this
#        script instead rebuilds the panel with build_panel(), so the check
#        holds even before the main fits have been run.
#    (c) "stayed near zero for long stretches", the idiosyncratic precisions
#        under the OLD Gamma(0.1, 0.1) base measure at p = 29.  The shipped
#        default is Gamma(1, 1), under which the bound is never reached, so
#        this is the one check that has to refit: four chains at as = bs = 0.1
#        with the main p = 29 settings otherwise.  About 5 minutes.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)

OUT <- file.path(TAB_DIR, "paper_checks.txt")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
sink(OUT, split = TRUE)
on.exit(sink(), add = TRUE)
cat(sprintf("paper_checks.R | %s | %s\n", R.version.string, format(Sys.time())))
cat(sprintf("VAPOR %s\n", as.character(utils::packageVersion("VAPOR"))))

#  the CFG of each dimension, read from its own script so there is one source
read_cfg <- function(dim){
   src <- sprintf("/Users/peterdunson/Desktop/vapor/applications/final_p%d.R", dim)
   for (q in parse(src))
      if (is.call(q) && identical(as.character(q[[1]])[1], "<-") &&
          identical(as.character(q[[2]]), "CFG")) return(eval(q[[3]]))
   stop(sprintf("no CFG <- list(...) in %s", basename(src))) }

# ============================================================================
#  (a) the Gamma(0.1, 0.1) base measure: prior mass below each level
# ============================================================================
cat("\n\n=== (a) Gamma(shape = 0.1, rate = 0.1) prior mass below each level ===\n")
cat("    the base measure of the DP on the idiosyncratic precisions as it was\n")
cat("    before the default moved to Gamma(1, 1)\n\n")
thr <- c(1e-2, 1e-4, 1e-6, 1e-8)
A <- data.frame(threshold = sprintf("%.0e", thr),
                prior_prob = round(stats::pgamma(thr, shape = 0.1, rate = 0.1), 4),
                one_in = sprintf("1 in %.1f", 1 / stats::pgamma(thr, 0.1, 0.1)))
print(A, row.names = FALSE)
mass_1e4 <- stats::pgamma(1e-4, shape = 0.1, rate = 0.1)
cat(sprintf("\n  pgamma(1e-4, 0.1, 0.1) = %.10f\n", mass_1e4))
cat(sprintf("  median = %.4e | mean = %.1f\n",
            stats::qgamma(0.5, 0.1, 0.1), 0.1 / 0.1))
cat(sprintf("  for comparison, Gamma(1, 1) puts %.3e below 1e-4\n",
            stats::pgamma(1e-4, shape = 1, rate = 1)))

# ============================================================================
#  (b) the NIPA accounting identity in the p = 100 training panel
# ============================================================================
cat("\n\n=== (b) NIPA accounting identity, p = 100 training panel ===\n")
CFG100 <- read_cfg(100)
D <- build_panel(CFG100)
S <- c("GDPC1", "PCECC96", "GPDIC1", "GCEC1", "EXPGSC1", "IMPGSC1")
stopifnot(all(S %in% D$vars))
Y <- D$train[, S, drop = FALSE]
cat(sprintf("\n  %d rows (%s), transformed and standardized on the training window\n",
            nrow(Y), paste(range(rownames(Y)), collapse = " .. ")))
cat(sprintf("  series: %s\n", paste(S, collapse = ", ")))

d  <- as.data.frame(Y)
m  <- stats::lm(GDPC1 ~ PCECC96 + GPDIC1 + GCEC1 + EXPGSC1 + IMPGSC1, data = d)
sm <- summary(m)
cat("\n  GDPC1 on the five expenditure components:\n")
print(data.frame(term = names(stats::coef(m)),
                 estimate = round(stats::coef(m), 4),
                 std_error = round(sm$coefficients[, 2], 4),
                 t = round(sm$coefficients[, 3], 2),
                 p = signif(sm$coefficients[, 4], 3), row.names = NULL),
      row.names = FALSE)
cat("\n")
print(data.frame(n = stats::nobs(m), df_residual = stats::df.residual(m),
                 R2 = round(sm$r.squared, 6), adj_R2 = round(sm$adj.r.squared, 6),
                 residual_SD = round(sm$sigma, 6),
                 F_stat = round(sm$fstatistic[1], 1), row.names = NULL),
      row.names = FALSE)

C  <- stats::cor(Y)
ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
cat("\n  the 6 x 6 sample correlation matrix:\n"); print(round(C, 3))
cat("\n  all six eigenvalues:\n")
print(data.frame(index = 1:6, eigenvalue = signif(ev, 6)), row.names = FALSE)
cat(sprintf("\n  smallest eigenvalue: %.8f\n", min(ev)))
cat(sprintf("  condition number (largest / smallest): %.2f | determinant: %.4e\n",
            max(ev) / min(ev), det(C)))
R2_val <- sm$r.squared; eig_val <- min(ev)

# ============================================================================
#  (c) p = 29 under BOTH base measures: how long the precisions stay near
#      zero, and what that does to the convergence of the partial correlations
# ============================================================================
cat("\n\n=== (c) p = 29, Gamma(0.1, 0.1) against the shipped Gamma(1, 1) ===\n")
CFG29 <- read_cfg(29)
cat(sprintf("\n  settings, from final_p29.R except as and bs: k = %d, rank = %d,\n",
            CFG29$K, CFG29$RANK))
cat(sprintf("  nrun = %d, burn = %d, thin = %d, jitter = %g, seeds %s,\n",
            CFG29$NRUN, CFG29$BURN, CFG29$THIN, CFG29$JIT,
            paste(range(CFG29$SEEDS), collapse = " to ")))
cat("  a_dir = b_dir = 0.1 in both; the two runs differ only in as and bs\n")
D29 <- build_panel(CFG29)
Y29 <- D29$train; vars29 <- D29$vars; p29 <- D29$p
LAB <- pair_labels(vars29)
UT  <- which(upper.tri(matrix(0, p29, p29)), arr.ind = TRUE)
cat(sprintf("  %d series, %d partial correlations\n", p29, length(LAB)))

LEVEL <- 1e-4
#  longest run of consecutive retained draws strictly below LEVEL
longest_run <- function(x){
   r <- rle(x < LEVEL)
   if (!any(r$values)) return(0L)
   max(r$lengths[r$values]) }

#  Fit the four chains at one base measure and reduce each immediately to the
#  two things used below: the idiosyncratic precision draws and the per-draw
#  partial correlations of every pair.
fit_setting <- function(ab){
   PS <- list(); DR <- list(); MIN <- list()
   for (sd_i in CFG29$SEEDS){
      t0 <- Sys.time()
      f <- vapor(Y29, lags = 1, k = CFG29$K, rank = CFG29$RANK, prior_B = "dl",
                 nrun = CFG29$NRUN, burn = CFG29$BURN, thin = CFG29$THIN,
                 seed = sd_i, jitter = CFG29$JIT, init = "random",
                 as = ab, bs = ab)
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      PS[[length(PS) + 1L]] <- f$residual_prec
      DR[[length(DR) + 1L]] <- pair_draws(f, p29)
      MIN[[length(MIN) + 1L]] <- data.frame(
         seed = sd_i, draws = nrow(f$residual_prec),
         share_below_1e4 = round(mean(f$residual_prec < LEVEL), 4),
         series_with_any = sum(colMeans(f$residual_prec < LEVEL) > 0),
         max_longest_run = max(apply(f$residual_prec, 2, longest_run)),
         floor_iters = f$floor_iters,
         min_prec = signif(min(f$residual_prec), 3), minutes = round(el, 2))
      cat(sprintf("  as = bs = %-3g seed %d: %.2f min\n", ab, sd_i, el))
      flush.console(); rm(f); invisible(gc(FALSE)) }
   list(ps = PS, dr = DR, per_chain = do.call(rbind, MIN)) }

#  per series: fraction of retained draws below LEVEL and the longest run,
#  per chain
series_table <- function(z){
   do.call(rbind, lapply(seq_along(z$ps), function(i){
      ps <- z$ps[[i]]
      data.frame(seed = CFG29$SEEDS[i], series = vars29,
                 frac_below_1e4 = round(colMeans(ps < LEVEL), 4),
                 longest_run = apply(ps, 2, longest_run),
                 min_prec = signif(apply(ps, 2, min), 3), row.names = NULL) })) }

#  R-hat and bulk ESS of every pair's per-draw partial correlation, with the
#  four chains as the chain dimension
conv_table <- function(z){
   M <- lapply(z$dr, as.matrix)
   do.call(rbind, lapply(seq_along(LAB), function(j){
      X <- do.call(cbind, lapply(M, function(a) a[, j]))   # draws x chains
      data.frame(pair = LAB[j], i = UT[j, 1], j2 = UT[j, 2],
                 Rhat = posterior::rhat(X),
                 bulk_ESS = posterior::ess_bulk(X), row.names = NULL) })) }

#  the split: pairs with at least one endpoint that ever goes below LEVEL
split_summary <- function(CV, flagged){
   touches <- vars29[CV$i] %in% flagged | vars29[CV$j2] %in% flagged
   grp <- ifelse(touches, ">=1 endpoint below 1e-4", "neither endpoint below 1e-4")
   do.call(rbind, lapply(c(">=1 endpoint below 1e-4", "neither endpoint below 1e-4"),
      function(g){
         k <- CV[grp == g, , drop = FALSE]
         if (!nrow(k)) return(data.frame(group = g, n_pairs = 0L,
            median_Rhat = NA_real_, max_Rhat = NA_real_, median_bulk_ESS = NA_real_,
            min_bulk_ESS = NA_real_, n_Rhat_gt_1.05 = 0L))
         data.frame(group = g, n_pairs = nrow(k),
            median_Rhat = round(median(k$Rhat), 3), max_Rhat = round(max(k$Rhat), 3),
            median_bulk_ESS = round(median(k$bulk_ESS), 0),
            min_bulk_ESS = round(min(k$bulk_ESS), 0),
            n_Rhat_gt_1.05 = sum(k$Rhat > 1.05)) })) }

RES <- list()
for (ab in c(0.1, 1.0)){
   tag <- sprintf("as = bs = %g", ab)
   cat(sprintf("\n  --- %s ---\n", tag))
   z  <- fit_setting(ab)
   ST <- series_table(z)
   CV <- conv_table(z)
   #  flagged from this run's own draws, pooled over its four chains
   flagged <- sort(unique(ST$series[ST$frac_below_1e4 > 0]))
   RES[[tag]] <- list(per_chain = z$per_chain, series = ST, conv = CV,
                      flagged = flagged)
   rm(z); invisible(gc(FALSE)) }

for (tag in names(RES)){
   R <- RES[[tag]]
   cat(sprintf("\n\n  ===== %s =====\n", tag))
   cat("\n  per chain:\n"); print(R$per_chain, row.names = FALSE)
   cat(sprintf("\n  series with any draw below 1e-4, pooled over the four chains: %d of %d\n",
               length(R$flagged), p29))
   if (length(R$flagged)) cat(sprintf("    %s\n", paste(R$flagged, collapse = ", ")))
   cat("\n  the 10 longest runs over all chains and series:\n")
   print(head(R$series[order(-R$series$longest_run), ], 10), row.names = FALSE)
   cat(sprintf("\n  pooled share of the 4 x %d x %d retained draws below 1e-4: %.4f\n",
               R$per_chain$draws[1], p29, mean(R$series$frac_below_1e4)))
   cat(sprintf("  longest run of consecutive draws below 1e-4, any series: %d\n",
               max(R$series$longest_run)))
   cat("\n  partial-correlation convergence, split on the flagged series\n")
   cat("  (R-hat and bulk ESS from the per-draw partial correlations, four chains):\n")
   print(split_summary(R$conv, R$flagged), row.names = FALSE)
   cat(sprintf("\n  all %d pairs together: median R-hat %.3f, max %.3f | median bulk ESS %.0f, min %.0f | R-hat > 1.05: %d\n",
               nrow(R$conv), median(R$conv$Rhat), max(R$conv$Rhat),
               median(R$conv$bulk_ESS), min(R$conv$bulk_ESS),
               sum(R$conv$Rhat > 1.05))) }

#  The split is only comparable across the two runs if the grouping is held
#  fixed, so repeat it with the Gamma(0.1, 0.1) flagged set applied to both.
FL01 <- RES[["as = bs = 0.1"]]$flagged
cat("\n\n  ===== both runs, split on the SAME series set =====\n")
cat(sprintf("  grouping fixed to the %d series flagged under Gamma(0.1, 0.1)\n",
            length(FL01)))
CMP <- do.call(rbind, lapply(names(RES), function(tag){
   x <- split_summary(RES[[tag]]$conv, FL01); cbind(setting = tag, x) }))
cat("\n"); print(CMP, row.names = FALSE)

TB <- RES[["as = bs = 0.1"]]$series

# ============================================================================
#  the three quoted statements, next to what this script computed
# ============================================================================
cat("\n\n=== the quoted statements and what this script reproduces ===\n\n")
Q <- data.frame(
   statement = c("a third of its mass below 1e-4",
                 "R^2 = 0.997",
                 "smallest eigenvalue 0.0013",
                 "stayed near zero for long stretches"),
   quoted = c("1/3 = 0.3333", "0.997", "0.0013", "(qualitative)"),
   computed = c(sprintf("%.4f", mass_1e4), sprintf("%.6f", R2_val),
                sprintf("%.6f", eig_val),
                sprintf("longest run %d draws, %.1f%% of draws below 1e-4",
                        max(TB$longest_run), 100 * mean(TB$frac_below_1e4))),
   row.names = NULL)
print(Q, row.names = FALSE)
cat(sprintf("\nwrote %s\n", OUT))
