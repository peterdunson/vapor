# ============================================================================
#  Combine the four per-setting sim_g parts into the single object the table
#  code expects.
#     Rscript applications/combine_sim_g.R
#  Writes sim_g_results.rds with the same fields sim_g.R writes:
#  results, paired, raw, settings, n_seeds, FIT_RANKS, r0, k.
# ============================================================================
PARTS <- "/Users/peterdunson/Desktop/vapor/simulation/sim_g_parts"
OUT   <- "/Users/peterdunson/Desktop/vapor/simulation/sim_g_parts/sim_g_results.rds"
WANT  <- list(c(50,35), c(100,70), c(150,105))   # p = 200 dropped
fs <- vapply(WANT, function(z)
   file.path(PARTS, sprintf("sim_g_part_p%d_T%d.rds", z[1], z[2])), character(1))
have <- file.exists(fs)
if (!all(have))
   stop(sprintf("missing %d of 3 parts: %s", sum(!have),
                paste(basename(fs[!have]), collapse = ", ")))
Z <- lapply(fs, readRDS)
tm <- do.call(rbind, lapply(Z, function(z) data.frame(
   setting = z$label, start = z$start, finish = z$finish,
   minutes = round(z$minutes, 1))))
cat("=== per-setting wall clock ===\n"); print(tm, row.names = FALSE)
cat(sprintf("\ntotal compute if run sequentially: %.1f h | longest setting: %.1f h\n",
            sum(tm$minutes)/60, max(tm$minutes)/60))
out <- list(
   results  = do.call(rbind, lapply(Z, function(z) z$results)),
   paired   = do.call(rbind, lapply(Z, function(z) z$paired)),
   raw      = setNames(lapply(Z, function(z) z$raw), vapply(Z, function(z) z$label, character(1))),
   settings = lapply(WANT, function(z) list(p = z[1], T = z[2],
                                            block = "T<p sweep (T/p=0.7)")),
   n_seeds  = 25L, FIT_RANKS = 3L, r0 = 3L, k = 5L,
   per_setting_timing = tm, run_at = format(Sys.time()),
   note = "four settings run as independent processes; sim_g.R re-seeds from the replicate index alone")
saveRDS(out, OUT)
cat(sprintf("\nwrote %s | %d result rows, %d paired rows\n", OUT,
            nrow(out$results), nrow(out$paired)))
print(unique(out$results[, c("p","T")]), row.names = FALSE)
