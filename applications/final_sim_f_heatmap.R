# ============================================================================
#  sim_f_heatmap_p150_T150.png, restyled to match sim_a_heatmap_rho95_T100.png.
#     Rscript applications/final_sim_f_heatmap.R
#
#  Same replicate and same data as sim_f.R draws (setting p = 150, T = 150,
#  replicate seed 1): sim_f.R re-seeds from the replicate index alone, so
#  refitting that one replicate here reproduces its panels exactly.  Only the
#  styling changes: sim_a's ggplot theme, diverging scale and horizontal
#  colour bar, with geom_raster so there are no cell borders -- at p = 150 the
#  per-cell borders pheatmap draws fill the panel with a grey mesh.
# ============================================================================
suppressMessages({library(VAPOR); library(PrecFact); library(bayesianVARs)
                  library(glasso); library(MASS); library(ggplot2); library(scales)})
SRC <- "/Users/peterdunson/Desktop/vapor/simulation/sim_f.R"
OUT <- "/Users/peterdunson/Desktop/vapor/applications/results/figures/sim_f_heatmap_p150_T150.png"
P <- 150L; TT <- 150L; REP <- 1L
cat(sprintf("VAPOR %s | as = %s, bs = %s | p = %d, T = %d, replicate %d\n",
            as.character(utils::packageVersion("VAPOR")),
            formals(vapor)$as, formals(vapor)$bs, P, TT, REP)); flush.console()

#  sim_f.R's definitions only (everything before its first top-level loop)
e <- new.env(parent = globalenv())
for (q in parse(SRC)){
   if (is.call(q) && identical(as.character(q[[1]])[1], "for")) break
   try(eval(q, envir = e), silent = TRUE) }
for (nm in ls(e)) assign(nm, get(nm, envir = e), envir = globalenv())

CACHE <- "/Users/peterdunson/Desktop/vapor/simulation/sim_g_parts/sim_f_heatmap_p150_panels.rds"
if (file.exists(CACHE)){
   pcs <- readRDS(CACHE)
   cat(sprintf("panels read from %s (no refit)\n", basename(CACHE)))
} else {
   t0 <- Sys.time()
   sim <- simulate_lowrank_calibrated(T = TT, p = P, seed = REP)
   fit <- fit_all(sim, P, seed = REP)
   pcs <- fit$pc
   saveRDS(pcs, CACHE)
   cat(sprintf("refit done in %.1f min; panels cached\n",
               as.numeric(difftime(Sys.time(), t0, units="mins")))) }

mats <- list(`True`          = pcs$true,
             `VAPOR (DL)`    = pcs[[sprintf("dl_r%d", FIT_RANKS[1])]],
             `VAPOR (ridge)` = pcs[[sprintf("rg_r%d", FIT_RANKS[1])]],
             `Kastner`       = pcs$kast,
             `Glasso`        = pcs$glasso)

#  ORIGINAL ordering: the figure this replaces was drawn by pheatmap with
#  cluster_rows = FALSE and cluster_cols = FALSE, i.e. the natural variable
#  order.  Keep it, so only the styling differs.
ord <- seq_len(P)
df <- do.call(rbind, lapply(names(mats), function(nm){
   M <- mats[[nm]][ord, ord]
   data.frame(row = rep(seq_len(P), times = P), col = rep(seq_len(P), each = P),
              value = as.vector(M), method = factor(nm, levels = names(mats))) }))
lim <- 1

g <- ggplot(df, aes(col, row, fill = value)) +
   geom_raster() +                                  # no cell borders
   facet_wrap(~ method, nrow = 1) +
   scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
                        midpoint = 0, limits = c(-lim, lim),
                        oob = scales::squish, na.value = "grey80", name = NULL,
                        breaks = c(-lim, 0, lim),
                        labels = sprintf("%.1f", c(-lim, 0, lim)),
                        guide = guide_colourbar(title = NULL, ticks = FALSE)) +
   scale_x_continuous(expand = c(0,0)) +
   scale_y_reverse(expand = c(0,0)) +
   coord_fixed() +
   theme_void(base_size = 10) +
   theme(strip.text        = element_text(size = 10, margin = margin(b = 4)),
         panel.background  = element_rect(fill = "white", colour = NA),
         plot.background   = element_rect(fill = "white", colour = NA),
         panel.border      = element_rect(fill = NA, colour = "grey45", linewidth = 0.4),
         panel.spacing     = unit(6, "pt"),
         legend.position   = "bottom",
         legend.key.width  = unit(40, "pt"),
         legend.key.height = unit(5, "pt"),
         legend.text       = element_text(size = 8),
         legend.margin     = margin(t = 2),
         plot.margin       = margin(3, 3, 3, 3))
ggsave(OUT, g, width = 6.5, height = 1.85, units = "in", dpi = 600, bg = "white")
cat(sprintf("wrote %s\n", OUT))
