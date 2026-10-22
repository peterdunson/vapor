# ============================================================================
#  Rebuild results/fits/kastner_p100_20k.rds and kastner_p100_20k_draws.rds.
#     Rscript applications/rebuild_kastner_p100_20k.R
#
#  WHY THIS SCRIPT EXISTS.  Both objects are inputs to run_all.sh that nothing
#  in the repository produced; the script that wrote them was never committed
#  and is not recoverable.  This is a RECONSTRUCTION from the config the two
#  .rds files themselves carry, reproduced verbatim below:
#
#     P 100 | K 25 | ITER 20000 | BURN 10000 | THIN 10 | NRUN 10000
#     VAR_JIT 0.001 | SEEDS_ALL 101..120 | BATCH 4 | NEED 4
#     WIN 1996Q1..2025Q4 | TRAIN 1996Q1..2020Q4 | TEST 2021Q1..2025Q4
#
#  THE SHORT CHAIN IS DELIBERATE.  Kastner at p = 100 is the one method that is
#  NOT run at the 100,000 iterations the other methods use: its chains abort
#  inside bayesianVARs' factor-SV sampler at that length.  20,000 iterations
#  with 10,000 burn-in and thin 10 retains 1,000 draws per chain, matching the
#  retained count of the other methods.
#
#  TWO OBJECTS, because the per-draw precisions are too large to keep inline:
#     kastner_p100_20k.rds        config, attempts, chains (Omega + settings)
#     kastner_p100_20k_draws.rds  the same chains plus pcor_mean_draws and a
#                                 prec_draws_file per chain, written under
#                                 results/fits/kastner_p100_prec_draws/
#
#  The windows and the series list come from final_p100.R's own CFG, so the
#  panel is the one build_panel() gives the main fits.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)

#  ---- CFG, from final_p100.R, overridden only where the stored config differs
SRC <- "/Users/peterdunson/Desktop/vapor/applications/final_p100.R"
CFG <- NULL
for (q in parse(SRC))
   if (is.call(q) && identical(as.character(q[[1]])[1], "<-") &&
       identical(as.character(q[[2]]), "CFG")) { CFG <- eval(q[[3]]); break }
if (is.null(CFG)) stop("no CFG <- list(...) found in final_p100.R")

#  the stored config, verbatim
CFG$ITER <- 20000L; CFG$BURN <- 10000L; CFG$THIN <- 10L; CFG$NRUN <- 10000L
CFG$SEEDS_ALL <- 101:120; CFG$BATCH <- 4L; CFG$NEED <- 4L
SMOKE <- nzchar(Sys.getenv("SMOKE"))
if (SMOKE) { CFG$ITER <- 400L; CFG$BURN <- 200L; CFG$THIN <- 2L; CFG$NRUN <- 200L
   CFG$NEED <- 1L; CFG$SEEDS_ALL <- 101:104
   cat("SMOKE: 400 iterations, 200 burn-in, thin 2, 1 chain. Code-path check\n",
       "       only; it writes *_SMOKE.rds and its numbers are meaningless.\n", sep = "") }

#  SMOKE writes to its own names so a code-path check cannot overwrite the
#  shipped objects the paper reports.
TAG     <- if (SMOKE) "_SMOKE" else ""
OUT_1   <- file.path(FITS_DIR, sprintf("kastner_p100_20k%s.rds", TAG))
OUT_2   <- file.path(FITS_DIR, sprintf("kastner_p100_20k_draws%s.rds", TAG))
DRAWDIR <- file.path(FITS_DIR, sprintf("kastner_p100_prec_draws%s", TAG))
dir.create(DRAWDIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("K = %d | ITER = %d | BURN = %d | THIN = %d | retained/chain = %d\n",
            CFG$K, CFG$ITER, CFG$BURN, CFG$THIN, (CFG$ITER - CFG$BURN) / CFG$THIN))

#  ---- panel -----------------------------------------------------------------
D <- build_panel(CFG)
Y_train <- D$train; Y_test <- D$test; vars <- D$vars; p <- D$p
stopifnot(p == 100L)

#  ---- one chain.  bvar() counts `draws` POST burn-in, so draws = ITER - BURN
one <- function(seed){
   set.seed(seed)
   Yj <- Y_train + matrix(rnorm(prod(dim(Y_train)), 0, CFG$VAR_JIT),
                          nrow(Y_train), ncol(Y_train))
   sp <- specify_prior_sigma(M = p, type = "factor",
                             factor_factors = as.integer(CFG$K),
                             factor_priorfacloadtype = "colwiseng",
                             factor_heteroskedastic = TRUE, quiet = TRUE)
   m <- bvar(Yj, lags = 1L, draws = as.integer(CFG$ITER - CFG$BURN),
             burnin = as.integer(CFG$BURN), thin = as.integer(CFG$THIN),
             prior_sigma = sp, quiet = TRUE)
   V <- vcov(m); nd <- dim(V)[4]
   acc <- matrix(0, p, p); pc <- matrix(0, p, p)
   prec <- array(NA_real_, c(p, p, nd))
   for (i in seq_len(nd)){
      Pi <- solve(V[1, , , i]); prec[, , i] <- Pi
      acc <- acc + Pi; pc <- pc + pcor_from_omega(Pi) }
   dimnames(acc) <- dimnames(pc) <- list(vars, vars)
   list(Omega = acc / nd, pcor_mean_draws = pc / nd, prec_draws = prec,
        n_draws = nd,
        settings = list(total_iterations = CFG$ITER, burnin = CFG$BURN,
                        thin = CFG$THIN,
                        draws_arg = as.integer(CFG$ITER - CFG$BURN),
                        retained = nd, factor_factors = as.integer(CFG$K))) }

#  ---- run until NEED chains have completed, over SEEDS_ALL ------------------
chains <- list(); attempts <- list()
for (s in CFG$SEEDS_ALL){
   if (length(chains) >= CFG$NEED) break
   t0 <- Sys.time()
   z <- tryCatch(one(s), error = function(e) structure(list(msg = conditionMessage(e)),
                                                      class = "chain_err"))
   el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
   if (inherits(z, "chain_err")){
      cat(sprintf("  seed %d FAILED after %.1f min: %s\n", s, el, z$msg))
      attempts[[length(attempts) + 1L]] <-
         data.frame(seed = s, status = "error", minutes = el, error = z$msg)
      next }
   #  stored as a BASENAME, not a path: final_shared.R's pcor_mean_prec()
   #  resolves it under results/fits/kastner_p100_prec_draws/
   bn <- sprintf("prec_draws_seed%d.rds", s)
   saveRDS(z$prec_draws, file.path(DRAWDIR, bn))
   z$prec_draws <- NULL; z$seed <- s; z$status <- "completed"; z$minutes <- el
   z$prec_draws_file <- bn
   chains[[length(chains) + 1L]] <- z
   attempts[[length(attempts) + 1L]] <-
      data.frame(seed = s, status = "completed", minutes = el, error = NA_character_)
   cat(sprintf("  seed %d: %.1f min, %d draws\n", s, el, z$n_draws)); flush.console() }
AT <- do.call(rbind, attempts)
if (length(chains) < CFG$NEED)
   stop(sprintf("only %d of %d chains completed", length(chains), CFG$NEED))

#  ---- write ------------------------------------------------------------------
cfg_store <- CFG[c("P","K","ITER","BURN","THIN","NRUN","VAR_JIT","SEEDS_ALL",
                   "BATCH","NEED","CSV","OUT_RDS","WIN","TRAIN","TEST")]
cfg_store$OUT_RDS <- basename(OUT_1)
cfg_store$VARS    <- vars
now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

saveRDS(list(config = cfg_store, attempts = AT,
             #  the stored layout: seed, status, minutes, Omega, n_draws
             chains = lapply(chains, function(q)
                q[c("seed","status","minutes","Omega","n_draws")]),
             run_at = now), OUT_1)
saveRDS(list(config = cfg_store,
             #  the stored layout, in order
             chains = lapply(chains, function(q)
                q[c("seed","status","minutes","Omega","pcor_mean_draws",
                    "n_draws","prec_draws_file")]),
             mspe_compare = NULL, draw_dir = DRAWDIR, run_at = now), OUT_2)
cat(sprintf("\nwrote %s\nwrote %s\ndraws in %s\n", OUT_1, OUT_2, DRAWDIR))
print(AT, row.names = FALSE)
