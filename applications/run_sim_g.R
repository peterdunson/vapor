# ============================================================================
#  Simulation G with the floored release build, at its existing settings.
#     Rscript applications/run_sim_g.R
#
#  sim_g.R is the script behind supplement table tab:sim-f and is the only
#  simulation the driver's CFG$SIMS does not cover.  It is run UNCHANGED, in
#  its own Rscript process, exactly as final_simulations.R isolates the other
#  six: vapor() is shadowed in the child's global environment by a wrapper
#  that calls VAPOR::vapor with the arguments it was given, records
#  floor_hits/floor_first, and returns the fit untouched, so the settings,
#  the seeds and the RNG stream are identical to running the script directly.
#
#  Writes sim_g_floor_report.rds.  The pre-floor sim_g_results.rds is copied
#  to results/fits/sim_g_results_prefloor.rds first, because the script
#  overwrites it in place.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)
SIMD <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/simulations")
SCRIPT <- "/Users/peterdunson/Desktop/vapor/simulation/sim_g.R"
stopifnot(file.exists(SCRIPT))

#  single-run lock.  A queued launcher may still be waiting to start this
#  script; the lock makes a second start a no-op rather than a second fit.
LOCK <- file.path(FITS_DIR, "sim_g.lock")
if (file.exists(LOCK)){
   cat(sprintf("sim_g already started at %s (%s exists); nothing to do\n",
               readLines(LOCK, warn = FALSE)[1], basename(LOCK)))
   quit(save = "no", status = 0) }
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), LOCK)

#  keep the pre-floor output: the script writes sim_g_results.rds in place
old <- file.path(SIMD, "sim_g_results.rds")
if (file.exists(old)){
   bk <- file.path(FITS_DIR, "sim_g_results_prefloor.rds")
   if (!file.exists(bk)){ file.copy(old, bk)
      cat(sprintf("pre-floor output copied to %s (written %s)\n", basename(bk),
                  format(file.mtime(old), "%Y-%m-%d %H:%M"))) }
}

lf  <- tempfile(fileext = ".rds")
drv <- file.path(tempdir(), "sim_g_driver.R")
writeLines(c(
   'suppressMessages({library(VAPOR); library(PrecFact)})',
   'graphics.off()',
   sprintf('LOGF <- "%s"; NM <- "sim_g"', lf),
   'ROWS <- list()',
   'vapor <- function(...){',
   '   f <- VAPOR::vapor(...)',
   '   ROWS[[length(ROWS)+1L]] <<- list(simulation = NM, p = nrow(f$Lambda),',
   '      series_floored = sum(f$floor_hits > 0),',
   '      iterations_floored = sum(f$floor_hits),',
   '      first_iteration = if (any(f$floor_first > 0)) min(f$floor_first[f$floor_first > 0]) else NA_integer_)',
   '   f }',
   sprintf('st <- tryCatch({ source("%s"); "completed" },', SCRIPT),
   '              error = function(e) paste("error:", conditionMessage(e)))',
   'graphics.off()',
   'saveRDS(list(status = st, rows = ROWS), LOGF)'), drv)

cat(sprintf("sim_g start %s\n", format(Sys.time(), "%Y-%m-%d %H:%M")))
t0 <- Sys.time()
system2("Rscript", drv, stdout = "", stderr = "")
el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
st <- "error: child produced no log"; rows <- list()
if (file.exists(lf)){ z <- readRDS(lf); st <- z$status; rows <- z$rows }
cat(sprintf("\nsim_g %s in %.1f minutes (%d vapor() calls logged)\n", st, el, length(rows)))

per_call <- if (length(rows)) do.call(rbind, lapply(rows, as.data.frame)) else NULL
if (is.null(per_call)){
   per_p <- data.frame(note = "no vapor() call was logged")
} else {
   per_p <- do.call(rbind, lapply(split(per_call, per_call$p), function(d)
      data.frame(p = d$p[1], vapor_fits = nrow(d),
                 fits_with_any_floor = sum(d$series_floored > 0),
                 mean_series_floored = round(mean(d$series_floored), 2),
                 max_series_floored  = max(d$series_floored),
                 total_floored_draws = sum(d$iterations_floored),
                 earliest_iteration  = if (all(is.na(d$first_iteration))) NA_integer_
                                       else min(d$first_iteration, na.rm = TRUE))))
   per_p <- per_p[order(per_p$p), ]
   cat("\n=== sim_g residual-precision floor hits, by dimension ===\n")
   print(per_p, row.names = FALSE)
   cat(sprintf("\ntotal: %d of %d fits floored at least one series\n",
               sum(per_call$series_floored > 0), nrow(per_call)))
}
saveRDS(list(status = st, minutes = el, per_call = per_call, per_p = per_p,
             run_at = format(Sys.time())),
        file.path(FITS_DIR, "sim_g_floor_report.rds"))
cat("saved sim_g_floor_report.rds\n")
