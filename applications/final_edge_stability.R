# ============================================================================
#  Item 2j: edge stability counted over DISTINCT fits only, plus how many of
#  seeds 101-106 select each edge at the reported configuration.
#     Rscript applications/final_edge_stability.R
#  Reads sweep_refit_p*.rds (which carries every setting's edge set).
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 205)
#  the reported configuration appears three times in the sweep grid (once in
#  each of the k, rank and a_L blocks); it is one fit and is counted once
DISTINCT <- list(
   `29`  = c("k=5", "k=8", "k=10", "k=15", "rank=2", "rank=4", "a_L=0.3", "a_L=1.0"),
   `100` = c("k=15", "k=25", "k=35", "rank=2", "rank=4", "a_L=0.3", "a_L=1.0"))
REPORTED <- list(`29` = "k=10", `100` = "k=25")
SEEDS <- sprintf("seed=%d", 101:106)
OUT <- list()

for (d in c(29L, 100L)){
   k <- as.character(d)
   #  prefer the sweep's own edge sets (item 18 now stores them); fall back to
   #  the separate refit object where that is all there is
   f <- file.path(FITS_DIR, sprintf("sweeps_p%d_final.rds", d))
   if (!file.exists(f)) f <- file.path(FITS_DIR, sprintf("sweep_refit_p%d.rds", d))
   if (!file.exists(f)){ cat(sprintf("p%d: no sweep object present\n", d)); next }
   ES <- readRDS(f)$edge_sets
   if (is.null(ES)){ cat(sprintf("p%d: %s carries no edge_sets\n", d, basename(f))); next }
   cat(sprintf("\np%d: edge sets read from %s\n", d, basename(f)))
   RPT <- { ff <- file.path(FITS_DIR, sprintf("p%d_results.rds", d))
            if (!file.exists(ff)) character(0) else {
               z <- readRDS(ff); e <- z$long$ci_edges; v <- z$vars
               i <- match(e$v1, v); j <- match(e$v2, v)
               paste(v[pmin(i,j)], v[pmax(i,j)], sep = "--") } }
   use <- DISTINCT[[k]]
   miss <- setdiff(use, names(ES))
   if (length(miss)) cat(sprintf("p%d: settings absent from the sweep: %s\n", d,
                                 paste(miss, collapse=", ")))
   use <- intersect(use, names(ES))
   tb <- table(unlist(ES[use]))
   st <- data.frame(edge = names(tb), n_settings = as.integer(tb),
                    out_of = length(use), reported = names(tb) %in% RPT,
                    row.names = NULL)
   st <- st[order(-st$n_settings, st$edge), ]
   cat(sprintf("\n=== 19. p = %d: selected edges, by fit ===\n", d))
   for (k2 in use){ e <- sort(ES[[k2]])
      cat(sprintf("  %-9s %d edges: %s\n", k2, length(e),
                  if (length(e)) paste(e, collapse=", ") else "(none)")) }
   cat(sprintf("\n=== 19. p = %d: edge stability over %d DISTINCT fits ===\n", d, length(use)))
   cat(sprintf("    %s\n", paste(use, collapse = ", ")))
   cat(sprintf("    the reported configuration (%s) is counted once, not three times\n",
               REPORTED[[k]]))
   print(st, row.names = FALSE)
   OUT[[sprintf("stability_p%d", d)]] <- st
   OUT[[sprintf("settings_p%d", d)]] <-
      data.frame(out_of = length(use), settings = paste(use, collapse = ", "))

   #  seeds 101-106 at the reported configuration
   sd_use <- intersect(SEEDS, names(ES))
   tb2 <- table(unlist(ES[sd_use]))
   sel <- data.frame(edge = names(tb2), n_seeds = as.integer(tb2),
                     out_of = length(sd_use), reported = names(tb2) %in% RPT,
                     row.names = NULL)
   sel <- sel[order(-sel$n_seeds, sel$edge), ]
   cat(sprintf("\n=== 19. p = %d: selected edges, by seed ===\n", d))
   for (k2 in sd_use){ e <- sort(ES[[k2]])
      cat(sprintf("  %-9s %d edges: %s\n", k2, length(e),
                  if (length(e)) paste(e, collapse=", ") else "(none)")) }
   cat(sprintf("\n=== 19. p = %d: seeds selecting each edge at the reported configuration ===\n", d))
   cat(sprintf("    %s\n", paste(sd_use, collapse = ", ")))
   print(sel, row.names = FALSE)
   cat(sprintf("    edges selected by all %d seeds: %d | by exactly one: %d\n",
               length(sd_use), sum(sel$n_seeds == length(sd_use)), sum(sel$n_seeds == 1)))
   OUT[[sprintf("seeds_p%d", d)]] <- sel
}
saveRDS(OUT, file.path(FITS_DIR, "edge_stability.rds"))
cat("\nsaved edge_stability.rds\n")
