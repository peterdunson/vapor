# ============================================================================
#  Rebuild results/fits/reused_p<dim>.rds from the data alone.
#     Rscript applications/rebuild_reused.R 29
#     Rscript applications/rebuild_reused.R 100
#
#  WHY THIS SCRIPT EXISTS.  reused_p29.rds and reused_p100.rds are the only
#  inputs to run_all.sh that nothing in the repository produced: they were
#  EXTRACTED (their own `source` field says "p<d>_results.rds |
#  fullrank_p<d>.rds") by a script that was never committed, from two fit
#  objects that no longer exist.  This script is a RECONSTRUCTION, not the
#  recovered original.  It refits the four reused competitors directly, using
#  the same panel, the same CFG and the same functions as the main fits, so the
#  object it writes is produced by the same code path the main pipeline uses:
#
#     Kastner             kastner_one(), from final_shared.R
#     VAR + Glasso        glasso_bic() on the full-rank residuals and on the
#                         rank-3 truncated residuals, the two first stages the
#                         pipeline needs
#     Sample covariance   solve(cov(Y_train) + 1e-8 I); the ridge is required,
#                         see SAMPLE_RIDGE below
#     Random walk         loss_rw()
#
#  WHAT IT CANNOT REPRODUCE EXACTLY.  Two things, both reported at the end:
#    1. p = 29's stored Kastner seeds are 101, 102, 104, 105, not 101:104 --
#       seed 103 aborted in the original run and the retry ladder replaced it.
#       This script uses the ladder too, so which seeds it lands on depends on
#       which chains abort in THIS run.  Pass SEEDS= to force a set.
#    2. p = 100's stored Kastner chains did NOT come from kastner_one() at
#       100,000 iterations.  Its `kastner_source` field names
#       kastner_p100_20k_draws.rds: 20,000 iterations, 10,000 burn-in, thin 10,
#       1000 retained draws per chain, because the 100,000-iteration chains
#       failed.  At p = 100 this script therefore defers to
#       rebuild_kastner_p100_20k.R and reuses its chains; run that first.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)
DIM  <- if (length(args)) as.integer(args[1]) else 29L
stopifnot(DIM %in% c(29L, 100L))

#  ---- CFG, taken from the dimension's own script so there is ONE definition -
#  final_p29.R / final_p100.R source final_pipeline.R and would run the whole
#  fit, so only their `CFG <- list(...)` assignment is evaluated here.
SRC <- sprintf("/Users/peterdunson/Desktop/vapor/applications/final_p%d.R", DIM)
CFG <- NULL
for (q in parse(SRC))
   if (is.call(q) && identical(as.character(q[[1]])[1], "<-") &&
       identical(as.character(q[[2]]), "CFG")) { CFG <- eval(q[[3]]); break }
if (is.null(CFG)) stop(sprintf("no CFG <- list(...) found in %s", basename(SRC)))
cat(sprintf("CFG read from %s: p = %d, K = %d, rank = %d, NRUN = %d, BURN = %d, THIN = %d\n",
            basename(SRC), CFG$P, CFG$K, CFG$RANK, CFG$NRUN, CFG$BURN, CFG$THIN))

OUT_RDS <- file.path(FITS_DIR, sprintf("reused_p%d.rds", DIM))

#  ---- panel, built by final_shared.R exactly as the main fits build it ------
D <- build_panel(CFG)
Y_train <- D$train; Y_test <- D$test; vars <- D$vars; p <- D$p
T_train <- nrow(Y_train); T_test <- nrow(Y_test)
Y_last_train <- Y_train[T_train, ]
stopifnot(p == CFG$P)
#  build_panel() returns no label, so form the stored one from the row names
qlab <- function(dstr){ z <- as.numeric(strsplit(dstr, "/")[[1]])
   sprintf("%.0fQ%.0f", z[3], ceiling(z[1]/3)) }
span <- function(M) sprintf("%s..%s, %d rows", qlab(rownames(M)[1]),
                            qlab(rownames(M)[nrow(M)]), nrow(M))
train_label <- span(Y_train); test_label <- span(Y_test)
cat(sprintf("panel: train %s | test %s\n", train_label, test_label))

#  ---- held-out sets, the same 100 as everywhere else ------------------------
hid_sets <- make_hid_sets(CFG$HID_SEED_BASE, CFG$N_HID, p, T_test, CFG$HID_FRAC)
cat(sprintf("hid_sets: %d sets x %d quarters, %d of %d series hidden per quarter\n",
            length(hid_sets), length(hid_sets[[1]]), length(hid_sets[[1]][[1]]), p))

#  ---- first-stage VAR: the two residual matrices the glasso rows need -------
vr <- fit_var_mean(Y_train, seed = CFG$SEEDS[1], draws = CFG$VAR_DRAWS,
                   burnin = CFG$VAR_BURNIN, jit = CFG$VAR_JIT)
rt <- rank_truncate_B(vr$B_hat, p, CFG$RANK)
resid_full <- vr$Y - vr$X %*% vr$B_hat     # MAIN rows
resid_rr   <- vr$Y - vr$X %*% rt$B         # supplement rows

#  ---- VAR + Glasso, both first stages ---------------------------------------
gl_full <- glasso_bic(resid_full, CFG$RHO_GRID, T_train)
gl_rr3  <- glasso_bic(resid_rr,   CFG$RHO_GRID, T_train)
dimnames(gl_full$Omega) <- dimnames(gl_rr3$Omega) <- list(vars, vars)
cat(sprintf("VAR + Glasso: full-rank rho %.6f (%d edges) | rank-3 rho %.6f (%d edges)\n",
            gl_full$rho, gl_full$n_nonzero, gl_rr3$rho, gl_rr3$n_nonzero))

#  ---- sample covariance and random walk -------------------------------------
#  A 1e-8 ridge, not a bare solve().  At p = 100, T_train = 100 the sample
#  covariance has rank 99 and solve() fails outright; the ridge is what the
#  stored objects were built with.  Verified: solve(cov(Y_train) + 1e-8 I)
#  reproduces the stored sample_cov_omega with max|diff| = 0 at BOTH p = 29
#  (where cov is full rank, and a bare solve differs by 1.3e-3) and p = 100.
SAMPLE_RIDGE <- 1e-8
Om_sample <- solve(cov(Y_train) + diag(SAMPLE_RIDGE, p))
dimnames(Om_sample) <- list(vars, vars)
loss_sample <- loss_avg(Om_sample, hid_sets, Y_test, p)
loss_rwalk  <- loss_rw(hid_sets, Y_test, Y_last_train, p)
loss_gl_rr3 <- loss_avg(gl_rr3$Omega, hid_sets, Y_test, p)

#  ---- Kastner ---------------------------------------------------------------
#  p = 100 reuses the 20,000-iteration chains; see the header.
K20 <- file.path(FITS_DIR, "kastner_p100_20k_draws.rds")
if (DIM == 100L) {
   if (!file.exists(K20))
      stop("p = 100 Kastner comes from kastner_p100_20k_draws.rds; run ",
           "applications/rebuild_kastner_p100_20k.R first")
   KD <- readRDS(K20)
   CH_kast <- KD$chains
   kast_source <- sprintf("%s (%d iterations, %d retained draws/chain)",
                          basename(K20), KD$config$ITER, CH_kast[[1]]$n_draws)
   cat(sprintf("Kastner: %d chain(s) reused from %s\n", length(CH_kast), basename(K20)))
} else {
   seeds <- if (nzchar(Sys.getenv("SEEDS")))
      as.integer(strsplit(Sys.getenv("SEEDS"), ",")[[1]]) else CFG$SEEDS
   kl <- run_with_ladder(function(sd) kastner_one(Y_train, CFG, sd),
                         seeds, CFG$MAX_REPL, "Kastner")
   CH_kast <- kl$chains
   kast_source <- sprintf("kastner_one() at %d iterations, %d burn-in, thin %d",
                          CFG$NRUN, CFG$BURN, CFG$THIN)
   cat(sprintf("Kastner: %d chain(s), seeds %s\n", length(CH_kast),
               paste(vapply(CH_kast, function(q){
                  s <- attr(q, "seed"); if (is.null(s)) q$seed else s }, numeric(1)),
                  collapse = ", ")))
}
kast_omegas <- lapply(CH_kast, function(q){
   O <- q$Omega; dimnames(O) <- list(vars, vars); O })
kast_losses <- lapply(kast_omegas, loss_avg, hid_sets = hid_sets,
                      Y_test = Y_test, p = p)

#  ---- write, in the field layout the pipeline reads -------------------------
out <- list(
   dim   = DIM,
   panel = list(train = train_label, test = test_label,
                Y_train = Y_train, Y_test = Y_test,
                digest = digest_mat(Y_train), vars = vars),
   kastner            = CH_kast,
   kastner_omegas     = kast_omegas,
   kastner_losses     = kast_losses,
   glasso_rank3       = gl_rr3,
   glasso_rank3_omega = list(gl_rr3$Omega),
   glasso_rank3_loss  = list(loss_gl_rr3),
   glasso_fullrank    = gl_full,
   sample_cov_omega   = list(Om_sample),
   sample_cov_loss    = list(loss_sample),
   random_walk_loss   = list(loss_rwalk),
   hid_sets           = hid_sets,
   source             = sprintf("rebuild_reused.R %d (reconstruction)", DIM),
   extracted_at       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
   kastner_source     = kast_source)
saveRDS(out, OUT_RDS)
cat(sprintf("\nwrote %s\n", OUT_RDS))

#  ---- comparison against whatever was there before, if a copy was kept ------
PREV <- paste0(OUT_RDS, ".prev")
if (file.exists(PREV)) {
   z <- readRDS(PREV)
   cmp <- function(nm, a, b) data.frame(quantity = nm,
      stored = signif(a, 10), rerun = signif(b, 10),
      abs_diff = signif(abs(a - b), 4), row.names = NULL)
   TAB <- rbind(
      cmp("glasso full-rank rho",  z$glasso_fullrank$rho, gl_full$rho),
      cmp("glasso rank-3 rho",     z$glasso_rank3$rho,    gl_rr3$rho),
      cmp("glasso full-rank edges", z$glasso_fullrank$n_nonzero, gl_full$n_nonzero),
      cmp("glasso rank-3 edges",   z$glasso_rank3$n_nonzero,     gl_rr3$n_nonzero),
      cmp("sample cov MSPE",       mean(z$sample_cov_loss[[1]]), mean(loss_sample)),
      cmp("random walk MSPE",      mean(z$random_walk_loss[[1]]), mean(loss_rwalk)),
      cmp("glasso rank-3 MSPE",    mean(z$glasso_rank3_loss[[1]]), mean(loss_gl_rr3)),
      cmp("Kastner MSPE (mean over chains)",
          mean(vapply(z$kastner_losses, mean, numeric(1))),
          mean(vapply(kast_losses, mean, numeric(1)))),
      cmp("Y_train digest sum",    z$panel$digest[["sum"]], digest_mat(Y_train)[["sum"]]))
   cat("\n=== stored vs rerun ===\n"); print(TAB, row.names = FALSE)
}
