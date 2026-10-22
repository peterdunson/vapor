# ============================================================================
#  Simulations A-F with the floored build, at their existing settings.
#     Rscript simulation/final_simulations.R
#
#  The simulation scripts in simulation/ are run UNCHANGED, so every setting
#  (N_REPS, NRUN, THIN, burn-in, seeds, DGPs) is exactly what it is now; the
#  only difference is that the installed VAPOR and PrecFact samplers now floor
#  the idiosyncratic precision at 1e-8.
#
#  The existing sim_*.R call vapor() and keep only the processed output, so the
#  floor counters would be discarded.  Rather than edit simulation/, this script
#  shadows vapor() with a logging wrapper in the environment each script is
#  sourced into: the wrapper calls VAPOR::vapor with the arguments it was given,
#  records floor_hits, and returns the fit untouched.  Settings, seeds and the
#  RNG stream are therefore identical to running the scripts directly.
#
#  Writes simulations_floor_report.rds and leaves each sim_*.R to write its own
#  output exactly where it already does.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
suppressMessages({library(VAPOR); library(PrecFact)})
options(width = 200)

CFG <- list(
   SIM_DIR = SIM_DIR,
   #  every simulation that feeds a paper table.  The script letters match the
   #  paper's: sim_a.R is Simulation A ... sim_f.R is Simulation F.
   #  sim_f is NOT in this queue: its three (p, T) settings are run as three
   #  independent processes by run_sim_f_setting.R and combined by
   #  combine_sim_f.R.  sim_f.R re-seeds from the replicate index alone, so
   #  per-setting runs reproduce the sequential result exactly.
   SIMS    = c("sim_a", "sim_b", "sim_c", "sim_d", "sim_e"),
   OUT     = SIM_FITS,
   OUT_RDS = "simulations_floor_report.rds")

FLOOR_LOG <- new.env(parent = emptyenv()); FLOOR_LOG$rows <- list()

#  Each simulation runs in its OWN R process.  In the first full run all six
#  were sourced into a single session and sim_f -- the sixth -- died with
#  "cannot open the connection" at its first png() after the earlier scripts
#  had run; sim_f completes normally in a fresh session.  Isolating the
#  processes removes the accumulated session state, and graphics.off() in the
#  child is a second guard.  The floor counters are written to a file by the
#  child and read back here.
run_one <- function(nm){
   path <- file.path(CFG$SIM_DIR, paste0(nm, ".R"))
   if (!file.exists(path)) return(data.frame(simulation = nm, status = "script missing",
                                             minutes = NA_real_))
   cat(sprintf("\n================ %s ================\n", nm)); flush.console()
   lf <- tempfile(fileext = ".rds")
   drv <- tempfile(fileext = ".R")
   writeLines(c(
      'suppressMessages({library(VAPOR); library(PrecFact)})',
      'graphics.off()',
      sprintf('LOGF <- "%s"; NM <- "%s"', lf, nm),
      'ROWS <- list()',
      'vapor <- function(...){',
      '   f <- VAPOR::vapor(...)',
      '   ROWS[[length(ROWS)+1L]] <<- list(simulation = NM, p = nrow(f$Lambda),',
      '      series_floored = sum(f$floor_hits > 0),',
      '      iterations_floored = sum(f$floor_hits),',
      '      first_iteration = if (any(f$floor_first > 0)) min(f$floor_first[f$floor_first > 0]) else NA_integer_)',
      '   f }',
      sprintf('st <- tryCatch({ source("%s"); "completed" },', path),
      '              error = function(e) paste("error:", conditionMessage(e)))',
      'graphics.off()',
      'saveRDS(list(status = st, rows = ROWS), LOGF)'), drv)
   t0 <- Sys.time()
   system2("Rscript", drv, stdout = "", stderr = "")
   el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
   st <- "error: child produced no log"
   if (file.exists(lf)){ z <- readRDS(lf); st <- z$status
      for (r in z$rows) FLOOR_LOG$rows[[length(FLOOR_LOG$rows)+1L]] <- r }
   unlink(c(lf, drv))
   data.frame(simulation = nm, status = st, minutes = el) }

#  Each sim_*.R now writes its own .rds straight into simulation/results/fits
#  and its per-setting heatmaps into simulation/results/figures, so there is
#  nothing left to relocate after a run; the files each simulation produced are
#  reported from the directory listing instead.
seen <- function() list.files(SIM_FITS, pattern = "[.]rds$")
t_begin <- Sys.time()
before <- seen()
status <- do.call(rbind, lapply(CFG$SIMS, function(nm){
   st <- run_one(nm)
   now <- seen(); st$wrote <- paste(setdiff(now, before), collapse = ", ")
   before <<- now; st }))
cat("\n=== simulation status ===\n"); print(status, row.names = FALSE)

rows <- FLOOR_LOG$rows
if (!length(rows)){
   per_call <- NULL; per_sim <- data.frame(note = "no vapor() call was logged")
} else {
   per_call <- do.call(rbind, lapply(rows, as.data.frame))
   per_sim <- do.call(rbind, lapply(split(per_call, per_call$simulation), function(d)
      data.frame(simulation = d$simulation[1], vapor_fits = nrow(d),
                 fits_with_any_floor = sum(d$series_floored > 0),
                 mean_series_floored = round(mean(d$series_floored), 2),
                 max_series_floored = max(d$series_floored),
                 total_iterations_floored = sum(d$iterations_floored),
                 earliest_iteration = suppressWarnings(min(d$first_iteration, na.rm = TRUE)))))
}
cat("\n=== residual-precision floor hits per simulation ===\n")
print(per_sim, row.names = FALSE)
saveRDS(list(config = CFG, status = status, per_call = per_call, per_sim = per_sim,
             run_at = format(Sys.time())), file.path(CFG$OUT, CFG$OUT_RDS))
cat(sprintf("\nsaved %s\n", CFG$OUT_RDS))
