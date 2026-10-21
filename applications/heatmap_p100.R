# ============================================================================
#  p=100 partial-correlation heatmaps.
#  Reads p100_results.rds, and p100_convention.rds for the two Bayesian panels;
#  refits nothing.
#  Boxes mark VAPOR's reported edges (95% credible interval excludes zero).
#  Box size scales with the dimension, so a single cell stays visible.
#  Labels are off by default: 35 edges over 100 variables would be
#  unreadable, and the network figure names them all.
# ============================================================================
RES_DIR  <- "/Users/peterdunson/Desktop/vapor/applications/results"
FIG_DIR  <- file.path(RES_DIR, "figures")
FITS_DIR <- file.path(RES_DIR, "fits")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
APP <- FITS_DIR            # fits are read from here
fit <- readRDS(file.path(APP, "p100_results.rds"))
#  item 4: every Bayesian panel is the mean of its per-draw partial
#  correlations pooled over all four chains.  p100_convention.rds carries those
#  matrices, computed from the SAME saved chains with no refit; only the data
#  behind the panels changes, no style element does.
CV <- file.path(APP, "p100_convention.rds")
if (file.exists(CV)){
   PCV <- readRDS(CV)$PCs
   fit$long$pc_vapor <- PCV$VAPOR
   cat(sprintf("VAPOR panel: mean of per-draw partial correlations, %d chains, %d draws\n",
               attr(PCV$VAPOR, "n_chains"), attr(PCV$VAPOR, "n_draws")))
   #  Kastner now has its per-draw precisions stored, so its panel is on the
   #  same convention; PC_KAST short-circuits the mean-precision path below.
   PC_KAST <- PCV$Kastner
} else cat("panels use p100_results.rds, whose pc_vapor and methods\ are already the one-convention matrices\n")
OUT <- file.path(FIG_DIR, "heatmap_p100_r3k25.png")

LABEL_TOP <- 0       # label the N strongest edges; 0 for none
LAB_CEX   <- 0.75

vars <- fit$vars; n <- length(vars)
ed   <- fit$long$ci_edges
ei   <- cbind(match(ed$v1, vars), match(ed$v2, vars))
up   <- t(apply(ei, 1, function(z) c(min(z), max(z))))   # upper triangle only
BOX_HW <- 0.5                             # ~1.67 at p=100
cat(sprintf("p=%d, %d reported edges, box half-width %.2f cells\n",
            n, nrow(ei), BOX_HW))

#  Kastner, when present, comes from THIS fit's own chains; see the note at
#  the kast_pc block below.
#  The "VAR + Glasso" panel is the RANK-3 fit, matching what that name now
#  means everywhere else. The title stays "VAR + Glasso".
avg_pc <- function(m) Reduce(`+`, m$pcor)/length(m$pcor)
#  Kastner comes from THIS study's own chains only.  At the pipeline's 100,000
#  iterations none completed (0 of 8 seeds), so the panel is taken from the
#  dedicated 20,000-iteration run in kastner_p100_20k.rds and its title says so
#  -- the panel is NOT at the same chain length as the VAPOR panel beside it.
#  (A pre-existing p100_long_competitors.rds from an earlier study was never
#  read -- it was run on a different standardization -- and is not in the repo.)
KF <- file.path(APP, "kastner_p100_20k.rds")
kast_pc <- NULL; kast_title <- "Kastner"
if (exists("PC_KAST") && !is.null(PC_KAST)){
   kast_pc <- PC_KAST
   cat("Kastner panel: mean of per-draw partial correlations, 4 chains, 4000 draws\n") }
if (is.null(kast_pc) && length(fit$methods[["Kastner (VAR+FSV)"]]$pcor) > 0){
   kast_pc <- avg_pc(fit$methods[["Kastner (VAR+FSV)"]])
   cat(sprintf("Kastner panel: %d chain(s) at the pipeline length\n",
               length(fit$methods[["Kastner (VAR+FSV)"]]$pcor)))
} else if (is.null(kast_pc) && file.exists(KF)){
   KK <- readRDS(KF)
   if (length(KK$chains)){
      kast_pc <- Reduce(`+`, lapply(KK$chains, function(z){
         Om <- z$Omega; d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc
      }))/length(KK$chains)
      #  the caption states the 20,000 iterations; the panel is just "Kastner"
      kast_title <- "Kastner"
      cat(sprintf("Kastner panel: %d chain(s) from %s at %d iterations\n",
                  length(KK$chains), basename(KF), KK$config$ITER)) } }
if (is.null(kast_pc)) cat("Kastner panel OMITTED: no Kastner chain available\n")
nodiag <- function(M){ diag(M) <- NA; M }
panels <- c(list("VAPOR" = nodiag(fit$long$pc_vapor)),
            if (!is.null(kast_pc)) setNames(list(nodiag(kast_pc)), kast_title),
            list("VAR + Glasso"      = nodiag(fit$methods[["VAR + Glasso"]]$pcor[[1]]),
                 "Sample covariance" = nodiag(fit$methods[["Sample covariance"]]$pcor[[1]])))
NP <- length(panels)

# ----------------------------------------------------------------------------
#  LAYOUT.  Drawn at the PRINTED size, 6.5 in wide, so every cex below is a
#  true point size (cex = pt/12 at the 12 pt device pointsize); nothing is
#  scaled by the figure width.  Laid out in inches with par(fig=), so all four
#  matrices are exactly square and exactly the same size.  Matches
#  heatmap_p29.R point for point.
# ----------------------------------------------------------------------------
W_IN      <- 6.5
PS        <- 12
pt        <- function(x) x/PS
TITLE_PT  <- 9; TICK_PT <- 6; BARLAB_PT <- 6
FAM       <- "serif"

LM <- 0.03; RM <- 0.03
GAP <- 0.06; GAP_BAR <- 0.22; BAR_W <- 0.12; TICK_W <- 0.24; PAD_R <- 0.05
TITLE_H <- 0.20; PAD_T <- 0.13; PAD_B <- 0.05

SQ <- (W_IN - (LM + RM + (NP-1)*GAP + GAP_BAR + BAR_W + TICK_W + PAD_R))/NP
H_IN <- PAD_T + TITLE_H + SQ + PAD_B
cat(sprintf("printed size %.2f x %.2f in; square %.3f in\n", W_IN, H_IN, SQ))

px <- numeric(NP); px[1] <- LM
if (NP > 1) for (q in 2:NP) px[q] <- px[q-1] + SQ + GAP
bar_x <- px[NP] + SQ + GAP_BAR
y0 <- PAD_B; y1 <- y0 + SQ

pal <- colorRampPalette(c("blue","white","red"))(101)
BOX_HW <- 0.5

png(OUT, width = W_IN, height = H_IN, units = "in", res = 600, bg = "white",
    pointsize = PS)
par(family = FAM)
nms <- names(panels)
for (q in seq_along(panels)){
   par(fig = c(px[q], px[q] + SQ, y0, y1)/c(W_IN, W_IN, H_IN, H_IN),
       mar = c(0,0,0,0), new = (q > 1), xpd = NA)
   M <- panels[[q]]; M[M > 1] <- 1; M[M < -1] <- -1
   image(x = 1:n, y = 1:n, z = t(M[n:1, ]), col = pal, zlim = c(-1,1),
         axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
   box(col = "grey60", lwd = 0.5)                      # thin grey panel border
   text(mean(par("usr")[1:2]), n + 0.5 + (TITLE_H*0.52)/(SQ/n),
        nms[q], cex = pt(TITLE_PT), font = 2, family = FAM, adj = c(0.5, 0))
   if (nms[q] != "VAPOR") next
   for (e in seq_len(nrow(up))){
      i <- up[e,1]; j <- up[e,2]
      rect(j-BOX_HW, n-i+1-BOX_HW, j+BOX_HW, n-i+1+BOX_HW,
           border = "black", lwd = 0.4) } }

par(fig = c(bar_x, bar_x + BAR_W, y0, y1)/c(W_IN, W_IN, H_IN, H_IN),
    mar = c(0,0,0,0), new = TRUE, xpd = NA)
image(x = 1, y = seq(-1, 1, length.out = 101),
      z = matrix(seq(-1, 1, length.out = 101), nrow = 1), col = pal,
      axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
box(col = "grey60", lwd = 0.5)
axis(4, at = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5)),
     las = 1, cex.axis = pt(TICK_PT), tick = FALSE, line = -0.62,
     family = FAM, col.axis = "black")
lh <- (BARLAB_PT*1.18/72)/(SQ/2)
text(1, 1 + 0.5*lh, "Correlation", cex = pt(BARLAB_PT), family = FAM, adj = c(0.5, 0))
text(1, 1 + 1.5*lh, "Partial",     cex = pt(BARLAB_PT), family = FAM, adj = c(0.5, 0))
dev.off()
cat(sprintf("wrote %s (%d panels)\n", OUT, NP))
