# ============================================================================
#  Item 2k: the chain-length table at either dimension, refitted on the
#  current panel with the current base measure.
#     Rscript applications/final_chain_length.R <29|100>
#
#  chain_length_source.rds is a summary of the old vapor_debug runs on the
#  2000Q2 panel and carries no fits, so it cannot be re-derived: configs B, A
#  and C are refitted here, four chains each, and config D reuses the main
#  p29_results.rds chains.  The p = 100 block of that file is left alone, its
#  panel being unchanged.
#
#  Quantities, matching the definitions in the file this replaces:
#    pairs_selected  pairs selected by AT LEAST ONE chain's own 95% CI
#    pooled_edges    pairs selected by the pooled 95% CI over all four chains
#    max/median_Rhat, min_bulk_ESS   over the pairs_selected set
#    MSPE            per-chain mean held-out loss, mean (SD) over chains
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)
DIM <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (!DIM %in% c(29L, 100L)) stop("usage: Rscript final_chain_length.R <29|100>")
CFG <- list(
   CFGS = list(B = c(nrun = 7500L,   burn = 2500L,  thin = 5L),
               A = c(nrun = 20000L,  burn = 10000L, thin = 10L),
               C = c(nrun = 50000L,  burn = 10000L, thin = 40L),
               D = c(nrun = 100000L, burn = 20000L, thin = 80L)),
   SEEDS = 101:104, K = if (DIM == 29L) 10L else 25L, RANK = 3L, A_L = 0.5, JIT = 5e-3,
   N_HID = 100L, HID_SEED_BASE = 90000L, HID_FRAC = 0.3)
LEN <- vapply(CFG$CFGS, function(z) sprintf("%d/%d/%d", z[1], z[2], z[3]), character(1))

R <- readRDS(file.path(FITS_DIR, sprintf("p%d_results.rds", DIM)))
vars <- R$vars; p <- length(vars); lab <- pair_labels(vars)
Y_train <- R$Y_train; Y_test <- R$Y_test
hid_sets <- make_hid_sets(CFG$HID_SEED_BASE, CFG$N_HID, p, nrow(Y_test), CFG$HID_FRAC)
cat(sprintf("p = %d: train %d x %d, test %d rows\n", DIM,
            nrow(Y_train), p, nrow(Y_test)))

one <- function(z, seed){
   f <- NULL; jj <- CFG$JIT; a <- 0L
   while (is.null(f) && a < 4L){ a <- a + 1L
      f <- tryCatch(vapor(Y_train, lags = 1, k = CFG$K, rank = CFG$RANK, prior_B = "dl",
                          a_L = CFG$A_L, b_L = CFG$A_L, nrun = z[["nrun"]],
                          burn = z[["burn"]], thin = z[["thin"]], seed = seed,
                          jitter = jj, init = "random"),
                    error = function(e) NULL)
      jj <- jj*2 }
   f }
W <- n_workers(NULL)
TODO <- c("B","A","C")
cat(sprintf("fitting configs %s, 4 chains each, on %d worker(s)\n",
            paste(TODO, collapse=", "), W))
jobs <- expand.grid(cfg = TODO, seed = CFG$SEEDS, stringsAsFactors = FALSE)
FIT <- parallel::mclapply(seq_len(nrow(jobs)), function(i){
   g <- jobs[i,]; t0 <- Sys.time()
   f <- one(CFG$CFGS[[g$cfg]], g$seed)
   list(cfg = g$cfg, seed = g$seed, fit = f,
        minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}, mc.cores = W, mc.preschedule = FALSE)

summarise <- function(chains, minutes, nm){
   DR <- lapply(chains, pair_draws, p = p)
   per <- lapply(DR, ci_indices)
   sel <- sort(unique(unlist(per)))                       # at least one chain
   pool <- do.call(rbind, DR); pooled <- ci_indices(pool) # pooled rule
   if (!length(sel)) return(NULL)
   Rh <- vapply(sel, function(j){
      Mx <- do.call(cbind, lapply(DR, function(M) M[,j])); posterior::rhat(Mx) }, numeric(1))
   Es <- vapply(sel, function(j){
      Mx <- do.call(cbind, lapply(DR, function(M) M[,j])); posterior::ess_bulk(Mx) }, numeric(1))
   ms <- vapply(chains, function(f) mean(loss_avg(omega_hat(f), hid_sets, Y_test, p)), numeric(1))
   list(nm = nm, setting = LEN[[nm]], sel = sel, edges = lab[pooled], lab = lab,
        max_rhat = max(Rh), at_max = lab[sel][which.max(Rh)],
        med_rhat = median(Rh), min_ess = min(Es), at_min = lab[sel][which.min(Es)],
        mspe = ms, minutes = minutes) }

rows <- list()
for (nm in TODO){
   ii <- which(vapply(FIT, function(z) identical(z$cfg, nm), logical(1)))
   ch <- lapply(FIT[ii], function(z) z$fit)
   mn <- vapply(FIT[ii], function(z) z$minutes, numeric(1))
   ok <- !vapply(ch, is.null, logical(1))
   cat(sprintf("  config %s: %d of %d chains completed\n", nm, sum(ok), length(ok)))
   if (!any(ok)) next
   rows[[nm]] <- summarise(ch[ok], mn[ok], nm) }
rows[["D"]] <- summarise(R$draws$vapor,
   vapply(R$draws$vapor, function(f) attr(f, "minutes"), numeric(1)), "D")

tab <- do.call(rbind, lapply(c("B","A","C","D"), function(nm){
   r <- rows[[nm]]; if (is.null(r)) return(NULL)
   data.frame(iterations = CFG$CFGS[[nm]][["nrun"]], nrun_burn_thin = r$setting,
              pairs_selected = length(r$sel), max_Rhat = round(r$max_rhat, 3),
              median_Rhat = round(r$med_rhat, 3), min_bulk_ESS = round(r$min_ess, 0),
              pooled_edges = length(r$edges),
              MSPE = sprintf("%.4f (%.4f)", mean(r$mspe), sd(r$mspe))) }))
tab <- tab[order(tab$iterations), ]
cat(sprintf("\n=== 2k. p = %d chain length ===\n", DIM)); print(tab, row.names = FALSE)
cm <- outer(names(rows), names(rows), Vectorize(function(a,b)
   length(intersect(rows[[a]]$edges, rows[[b]]$edges))))
dimnames(cm) <- list(names(rows), names(rows))
cat("\n  pooled edges shared between configurations:\n"); print(cm)
saveRDS(list(config = CFG, rows = rows, table = tab, common = cm,
             run_at = format(Sys.time())),
        file.path(FITS_DIR, sprintf("chain_length_p%d.rds", DIM)))
cat(sprintf("\nsaved chain_length_p%d.rds\n", DIM))
