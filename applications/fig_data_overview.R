# ============================================================================
#  Three descriptive figures for the FRED-QD panels, drawn from the SAME data
#  the fits use: build_panel() in final_shared.R, standardized with the
#  training moments.
#     Rscript applications/fig_data_overview.R
#  Writes results/figures/data_panel_heatmap.png, data_eigen_spectrum.png and
#  data_overview.png.  Reads no fit and refits nothing.
#
#  Style follows termstructure_p100.png: base graphics, the default sans
#  family, bold panel headings via mtext(font = 2).  All three are 6.5 in wide
#  at 600 dpi; pointsize 9 keeps the apparent text size of the 13.5 in figure
#  at this width.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 205)

WHICH <- commandArgs(trailingOnly = TRUE)
if (!length(WHICH)) WHICH <- c("heatmap", "spectrum", "overview")
want <- function(x) x %in% WHICH

W_IN <- 6.5; DPI <- 600; PS <- 9
CLIP <- 3
PAL  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
RECESS <- list(c(2001.00, 2001.75), c(2007.75, 2009.25), c(2020.00, 2020.25))
SPLIT  <- 2020.875            # between 2020Q4 and 2021Q1

read_cfg <- function(d){
   for (q in parse(sprintf("/Users/peterdunson/Desktop/vapor/applications/final_p%d.R", d)))
      if (is.call(q) && identical(as.character(q[[1]])[1], "<-") &&
          identical(as.character(q[[2]]), "CFG")) return(eval(q[[3]]))
   stop("no CFG") }

#  the full 1996Q1..2025Q4 panel: train and test are the two halves of the one
#  matrix build_panel() standardizes with the training mean and SD
full_panel <- function(d){
   D <- build_panel(read_cfg(d))
   Z <- rbind(D$train, D$test)
   qd <- function(s){ z <- as.numeric(strsplit(s, "/")[[1]]); z[3] + (z[1]-1)/12 }
   tt <- vapply(rownames(Z), qd, numeric(1))
   stopifnot(!is.unsorted(tt), all(abs(diff(tt) - 0.25) < 1e-8))
   list(Z = Z, t = tt, vars = D$vars, p = D$p,
        n_train = nrow(D$train), n_test = nrow(D$test), train = D$train) }

G   <- read.csv("/Users/peterdunson/Desktop/vapor/applications/fredqd_groups.csv",
                stringsAsFactors = FALSE)
GC  <- setNames(G$group_code, G$series); GN <- setNames(G$group_name, G$series)

P100 <- full_panel(100L); P29 <- full_panel(29L)
cat(sprintf("\np = 100 panel: %d quarters (%d train + %d test), %d series\n",
            nrow(P100$Z), P100$n_train, P100$n_test, P100$p))

# ---------------------------------------------------------------------------
#  row order: group, then PC1 loading of the TRAINING correlation matrix
# ---------------------------------------------------------------------------
Ctr <- cor(P100$train)
e1  <- eigen(Ctr, symmetric = TRUE)$vectors[, 1]
if (sum(e1) < 0) e1 <- -e1                      # sign convention: mostly positive
names(e1) <- P100$vars
ORD <- do.call(c, lapply(sort(unique(GC[P100$vars])), function(g){
   s <- P100$vars[GC[P100$vars] == g]
   s[order(-e1[s])] }))
gcode <- GC[ORD]; gname <- GN[ORD]
runs  <- rle(as.integer(gcode))
GTAB  <- data.frame(order = seq_along(runs$values), group_code = runs$values,
                    group_name = GN[ORD][cumsum(runs$lengths)],
                    n_series = runs$lengths, row.names = NULL)
#  compact names for the narrow panel of data_overview.png, where the full
#  names do not fit beside a 4.2 in strip
SHORTMAP <- c("NIPA"="NIPA", "Industrial Production"="Ind. production",
   "Employment and Unemployment"="Employment", "Housing"="Housing",
   "Inventories, Orders, and Sales"="Inventories", "Prices"="Prices",
   "Earnings and Productivity"="Earnings", "Interest Rates"="Interest rates",
   "Money and Credit"="Money, credit", "Exchange Rates"="Exchange rates",
   "Other"="Other", "Stock Markets"="Stock markets")
GTAB$short <- unname(SHORTMAP[GTAB$group_name])
cat("\n=== group order and size, top to bottom ===\n"); print(GTAB, row.names = FALSE)

M <- t(P100$Z[, ORD, drop = FALSE])             # series x quarter
cat(sprintf("\nheat strip: %d series x %d quarters\n", nrow(M), ncol(M)))

# ---------------------------------------------------------------------------
#  clipping report
# ---------------------------------------------------------------------------
nclip <- sum(abs(M) > CLIP)
cat(sprintf("\n=== cells clipped at +/- %g ===\n  %d of %d (%.3f%%)\n",
            CLIP, nclip, length(M), 100*nclip/length(M)))
idx <- which(abs(M) > CLIP, arr.ind = TRUE)
qlab <- function(tv) sprintf("%.0fQ%.0f", floor(tv + 1e-8),
                             floor((tv - floor(tv + 1e-8))*4 + 1e-8) + 1)
TOP <- data.frame(series = rownames(M)[idx[,1]], quarter = qlab(P100$t[idx[,2]]),
                  value = M[idx], row.names = NULL)
TOP <- head(TOP[order(-abs(TOP$value)), ], 10); TOP$value <- round(TOP$value, 3)
cat("\n  10 largest absolute values:\n"); print(TOP, row.names = FALSE)

# ---------------------------------------------------------------------------
#  eigen spectra
# ---------------------------------------------------------------------------
ev100 <- eigen(Ctr, symmetric = TRUE, only.values = TRUE)$values
ev29  <- eigen(cor(P29$train), symmetric = TRUE, only.values = TRUE)$values
Tn    <- nrow(P100$train)
ref   <- function(p) (1 + sqrt(p/Tn))^2
cat(sprintf("\n=== eigenvalues of the training correlation matrix (T = %d) ===\n", Tn))
print(data.frame(index = 1:10, p29 = round(ev29[1:10], 4),
                 p100 = round(ev100[1:10], 4)), row.names = FALSE)
cat(sprintf("\n  reference (1 + sqrt(p/T))^2: p = 29 -> %.4f | p = 100 -> %.4f\n",
            ref(29), ref(100)))
cat(sprintf("  eigenvalues above it: p = 29 -> %d of 29 | p = 100 -> %d of 100\n",
            sum(ev29 > ref(29)), sum(ev100 > ref(100))))
cat(sprintf("  smallest eigenvalue at p = 100: %.4e (rank %d of 100)\n",
            min(ev100), qr(cor(P100$train))$rank))

#  the singular one is drawn as a marker on the axis floor, not plotted
pos100 <- ev100[ev100 > 1e-12]; FLOOR <- min(c(ev29, pos100)) / 2.2
COL29 <- "#2E5E9E"; COL100 <- "#B2182B"

# ---------------------------------------------------------------------------
#  drawing helpers
# ---------------------------------------------------------------------------
#  decimal years are not an exactly regular grid in floating point, which
#  useRaster = TRUE requires, so the strip is drawn on column indices and the
#  time axis is mapped onto them.
col_of <- function(y) approx(P100$t, seq_along(P100$t), xout = y, rule = 2)$y

draw_strip <- function(cex_lab = 0.62, tag = NULL, mar = c(2.0, 8.8, 1.5, 0.6),
                       short = FALSE){
   tt <- P100$t; nq <- length(tt); nr <- nrow(M)
   par(mar = mar, mgp = c(1.5, 0.35, 0), tcl = -0.2)
   Mc <- pmax(pmin(M, CLIP), -CLIP)
   image(x = seq_len(nq), y = seq_len(nr), z = t(Mc[nr:1, , drop = FALSE]),
         col = PAL, zlim = c(-CLIP, CLIP), axes = FALSE, xlab = "", ylab = "",
         useRaster = TRUE)
   #  NBER recessions
   for (r in RECESS)
      rect(col_of(r[1]) - 0.5, 0.5, col_of(r[2]) + 0.5, nr + 0.5,
           col = adjustcolor("grey35", alpha.f = 0.16), border = NA)
   #  group separators and names
   cuts <- nr - cumsum(runs$lengths) + 0.5
   abline(h = cuts[-length(cuts)], col = "grey25", lwd = 0.5)
   #  group names: a label for a one- or two-series group sits almost on top
   #  of its neighbour, so push them apart to a minimum spacing first
   mids <- nr - (cumsum(runs$lengths) - runs$lengths/2) + 0.5
   #  two passes: push down to the minimum gap, then clamp the bottom and push
   #  back up, so the last group keeps its label instead of falling off the axis
   ng   <- length(mids)
   gap  <- min(nr * 0.034, (nr - 1) / max(ng - 1, 1))
   lo   <- 0.5 + gap/2; hi <- nr + 0.5 - gap/2
   ypos <- mids
   for (i in seq_len(ng)[-1])
      if (ypos[i-1] - ypos[i] < gap) ypos[i] <- ypos[i-1] - gap
   if (ypos[ng] < lo){
      ypos[ng] <- lo
      for (i in rev(seq_len(ng))[-1])
         if (ypos[i] - ypos[i+1] < gap) ypos[i] <- ypos[i+1] + gap }
   ypos <- pmin(pmax(ypos, lo), hi)
   labs <- if (short) GTAB$short else GTAB$group_name
   axis(2, at = ypos, labels = labs, tick = FALSE, las = 1,
        cex.axis = cex_lab, line = -0.4)
   #  a leader where a label had to move
   moved <- abs(ypos - mids) > 1e-8
   if (any(moved)) segments(par("usr")[1], ypos[moved], 0.4, mids[moved],
                            col = "grey55", lwd = 0.4, xpd = NA)
   #  train | test, between 2020Q4 and 2021Q1
   xs <- col_of(SPLIT)
   abline(v = xs, lty = 2, lwd = 0.9, col = "grey15")
   text(xs, nr + 0.5, "train | test", pos = 2, offset = 0.25,
        cex = cex_lab * 0.95, col = "grey15", xpd = NA)
   yrs <- seq(2000, 2025, by = 5)
   axis(1, at = col_of(yrs), labels = yrs, cex.axis = cex_lab,
        lwd = 0, lwd.ticks = 0.6)
   box(col = "grey35", lwd = 0.6)
   if (!is.null(tag)) mtext(tag, side = 3, adj = 0, line = 0.25, font = 2,
                            cex = cex_lab * 1.15) }

draw_cbar <- function(cex_lab = 0.62, mar = c(1.9, 8.8, 0.3, 0.6)){
   par(mar = mar)
   z <- seq(-CLIP, CLIP, length.out = 201)
   image(x = z, y = 1, z = matrix(z, ncol = 1), col = PAL, zlim = c(-CLIP, CLIP),
         axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
   box(col = "grey35", lwd = 0.5)
   axis(1, at = c(-3, -1.5, 0, 1.5, 3), cex.axis = cex_lab, lwd = 0,
        lwd.ticks = 0.6, mgp = c(1.5, 0.25, 0), tcl = -0.18)
   mtext(sprintf("standardized value (clipped at +/- %g)", CLIP), side = 1,
         line = 1.05, cex = cex_lab) }

draw_eigen <- function(cex_lab = 0.62, tag = NULL,
                       mar = c(2.9, 3.4, 1.5, 0.8), lab_line = c(1.55, 2.05)){
   par(mar = mar, mgp = c(1.9, 0.4, 0), tcl = -0.2)
   ylim <- c(FLOOR, max(ev100) * 1.35)
   plot(NA, xlim = c(1, 100), ylim = ylim, log = "y", axes = FALSE,
        xlab = "", ylab = "")
   abline(h = ref(29),  lty = 2, lwd = 0.8, col = COL29)
   abline(h = ref(100), lty = 2, lwd = 0.8, col = COL100)
   lines(seq_along(ev29),  ev29,  col = COL29,  lwd = 1.3)
   points(seq_along(ev29), ev29,  col = COL29,  pch = 16, cex = 0.32)
   lines(seq_along(pos100), pos100, col = COL100, lwd = 1.3)
   points(seq_along(pos100), pos100, col = COL100, pch = 16, cex = 0.32)
   #  the singular eigenvalue: a distinct marker on the floor
   points(100, FLOOR, pch = 4, col = COL100, cex = 0.8, lwd = 1.5)
   text(100, FLOOR, "rank 99: singular", pos = 2, offset = 0.35,
        cex = cex_lab * 0.92, col = COL100)
   axis(1, at = c(1, 20, 40, 60, 80, 100), cex.axis = cex_lab,
        lwd = 0, lwd.ticks = 0.6)
   ex <- c(1e-3, 1e-2, 1e-1, 1, 10, 100)
   ex <- ex[ex >= ylim[1] & ex <= ylim[2]]
   axis(2, at = ex, labels = formatC(ex, format = "g"), las = 1,
        cex.axis = cex_lab, lwd = 0, lwd.ticks = 0.6)
   mtext("eigenvalue index", side = 1, line = lab_line[1], cex = cex_lab)
   mtext("eigenvalue", side = 2, line = lab_line[2], cex = cex_lab)
   #  one annotation for both dashed lines, set above the upper one and away
   #  from the legend
   text(46, ref(100) * 1.9, "largest eigenvalue under independence",
        cex = cex_lab * 0.92, col = "grey25")
   legend("bottomleft", legend = c("p = 29", "p = 100"), col = c(COL29, COL100),
          lwd = 1.3, pch = 16, pt.cex = 0.32, bty = "n", cex = cex_lab,
          inset = c(0.02, 0.04), seg.len = 1.4)
   box(col = "grey35", lwd = 0.6)
   if (!is.null(tag)) mtext(tag, side = 3, adj = 0, line = 0.25, font = 2,
                            cex = cex_lab * 1.15) }

# ---------------------------------------------------------------------------
#  1. the heat strip
# ---------------------------------------------------------------------------
O1 <- file.path(FIG_DIR, "data_panel_heatmap.png")
if (want("heatmap")) {
png(O1, width = W_IN, height = 4.0, units = "in", res = DPI, bg = "white",
    pointsize = PS)
layout(matrix(1:2, ncol = 1), heights = c(1, 0.17))
draw_strip(cex_lab = 0.58); draw_cbar(cex_lab = 0.58)
dev.off(); cat(sprintf("\nwrote %s\n", O1)) }

# ---------------------------------------------------------------------------
#  2. the eigen spectrum
# ---------------------------------------------------------------------------
O2 <- file.path(FIG_DIR, "data_eigen_spectrum.png")
if (want("spectrum")) {
png(O2, width = W_IN, height = 3.0, units = "in", res = DPI, bg = "white",
    pointsize = PS)
draw_eigen()
dev.off(); cat(sprintf("wrote %s\n", O2)) }

# ---------------------------------------------------------------------------
#  3. the two side by side.  Drawn at the printed size, 6.5 in wide, with the
#     same point-size convention as the other paper figures (heatmap_p100.R):
#     device pointsize 12, so cex = pt/12 is a true point size.
# ---------------------------------------------------------------------------
O3 <- file.path(FIG_DIR, "data_overview.png")
if (want("overview")) {
H_IN <- 3.4; PS3 <- 12
FAM  <- "sans"
pt   <- function(x) x / 12               # point size -> cex
BODY <- pt(8); TAG <- pt(10)
SEP  <- "grey55"
#  device-fraction geometry
HX  <- c(0.205, 0.600)                   # heat strip, left and right
HY  <- c(0.228, 0.852)                   # heat strip, bottom and top
RY  <- c(0.864, 0.896)                   # recession strip
KY  <- c(0.906, 0.952)                   # train / test brackets
BX  <- c(0.737, 0.990); BY <- c(0.228, 0.852)    # (b): same top and bottom
CBW <- 0.40 * diff(HX)
CBX <- c(HX[1] + (diff(HX) - CBW)/2, HX[1] + (diff(HX) + CBW)/2)
#  Colour bar: TOP fixed so the gap to the year axis is unchanged, thickness
#  set in inches, and the bottom margin sized so the plot region is exactly
#  that thick; the tick labels and title then sit in the margin below it.
CB_TOP   <- 0.1240                      # device fraction, bar top
CB_THICK <- 0.12 / H_IN                 # 0.12 in of a 3.4 in canvas
CB_MARLN <- 1.397                        # bottom margin, in mar lines
CBF <- c(0.012, CB_TOP)                 # colour-bar figure region, bar + text

png(O3, width = W_IN, height = H_IN, units = "in", res = DPI, bg = "white",
    pointsize = PS3)
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", family = FAM)
#  This figure reorders the GROUP BLOCKS so that no two small groups are
#  adjacent; within each group the row order is unchanged.  Two labels of
#  adjacent blocks collide whenever (n_i + n_j)/2 is below one line height, so
#  every small group is placed between two large ones and each label can then
#  sit exactly on its own block midpoint.  data_panel_heatmap.png keeps the
#  original order by group code and is not affected.
GORDER <- c("NIPA", "Other", "Industrial production", "Earnings", "Employment",
            "Exchange rates", "Money and credit", "Stock markets", "Housing",
            "Inventories and orders", "Prices", "Interest rates")
FULLNAME <- c(NIPA = "NIPA", `Industrial production` = "Industrial Production",
   Employment = "Employment and Unemployment", Housing = "Housing",
   `Inventories and orders` = "Inventories, Orders, and Sales", Prices = "Prices",
   Earnings = "Earnings and Productivity", `Interest rates` = "Interest Rates",
   `Money and credit` = "Money and Credit", `Exchange rates` = "Exchange Rates",
   Other = "Other", `Stock markets` = "Stock Markets")
ORD2 <- do.call(c, lapply(GORDER, function(g){
   s <- P100$vars[GN[P100$vars] == FULLNAME[[g]]]
   s[order(-e1[s])] }))
stopifnot(setequal(ORD2, P100$vars), length(ORD2) == P100$p)
M2    <- t(P100$Z[, ORD2, drop = FALSE])
runs2 <- rle(as.character(GN[ORD2]))
GLAB  <- GORDER
stopifnot(identical(unname(FULLNAME[GLAB]), runs2$values))

nq <- length(P100$t); nr <- nrow(M2)
Mc <- pmax(pmin(M2, CLIP), -CLIP)
xs <- col_of(SPLIT)

## ---- (a) the heat strip ---------------------------------------------------
par(fig = c(HX, HY), new = FALSE)
image(x = seq_len(nq), y = seq_len(nr), z = t(Mc[nr:1, , drop = FALSE]),
      col = PAL, zlim = c(-CLIP, CLIP), axes = FALSE, xlab = "", ylab = "",
      useRaster = TRUE)
usr <- par("usr")
cuts <- nr - cumsum(runs2$lengths) + 0.5
abline(h = cuts[-length(cuts)], col = SEP, lwd = 0.3)
abline(v = xs, lty = 2, lwd = 0.8, col = "grey20")
yrs <- seq(2000, 2025, by = 5)
axis(1, at = col_of(yrs), labels = yrs, cex.axis = BODY, lwd = 0,
     lwd.ticks = 0.4, mgp = c(1.4, 0.22, 0), tcl = -0.14, family = FAM)
box(col = "grey35", lwd = 0.4)

#  Group labels: all at 8 pt, right-aligned at one fixed gap from the strip,
#  each centred exactly on its own block midpoint.  The reordering above makes
#  that feasible, so there is no displacement and no leader tick.
blk_hi <- nr - c(0, cumsum(runs2$lengths)[-length(runs2$lengths)]) + 0.5  # top row
blk_lo <- nr - cumsum(runs2$lengths) + 0.5                                # bottom
mids   <- (blk_lo + blk_hi) / 2
ng     <- length(mids)
row_pt <- diff(HY) * H_IN * 72 / nr          # printed points per strip row
hgt    <- 8 / row_pt * 1.03                  # one 8 pt line, in strip rows
ypos   <- mids                               # exactly the midpoints
GAPX   <- 0.016 * nq
text(usr[1] - GAPX, ypos, GLAB, adj = 1, cex = BODY, family = FAM, xpd = NA)

#  ---- checks ---------------------------------------------------------------
half <- hgt/2
chk <- data.frame(group = GLAB, n = runs2$lengths,
   block_lo = round(blk_lo, 1), block_hi = round(blk_hi, 1),
   mid = round(mids, 1), y = round(ypos, 1),
   centre_in_block = ypos >= blk_lo & ypos <= blk_hi, row.names = NULL)
ovl <- c(ypos[-ng] - half < ypos[-1] + half, FALSE)   # bbox of i meets i+1
chk$overlaps_next <- ovl
chk$in_strip <- ypos - half >= 0.5 - 1e-6 & ypos + half <= nr + 0.5 + 1e-6
cat("\n=== (a) group label checks ===\n"); print(chk, row.names = FALSE)
cat(sprintf("  every centre inside its own block: %s\n", all(chk$centre_in_block)))
cat(sprintf("  no two label boxes overlap:        %s\n", !any(chk$overlaps_next)))
cat(sprintf("  all inside the strip:              %s\n", all(chk$in_strip)))
cat(sprintf("  smallest gap between centres %.2f rows, one line is %.2f rows\n",
            min(ypos[-ng] - ypos[-1]), hgt))

## ---- (a) recession strip --------------------------------------------------
par(fig = c(HX, RY), new = TRUE)
plot(NA, xlim = c(0.5, nq + 0.5), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
for (r in RECESS)
   rect(col_of(r[1]) - 0.5, 0, col_of(r[2]) + 0.5, 1, col = "grey62", border = NA)
abline(v = xs, lty = 2, lwd = 0.8, col = "grey20")
box(col = "grey35", lwd = 0.4)
text(0.5 - GAPX, 0.5, "recession", adj = 1, cex = BODY, family = FAM, xpd = NA)

## ---- (a) training / test brackets ----------------------------------------
par(fig = c(HX, KY), new = TRUE)
plot(NA, xlim = c(0.5, nq + 0.5), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
GAPB <- 0.012 * nq                            # the gap at the boundary
yb <- 0.26; tick <- 0.17
brk <- function(x0, x1, lab){
   segments(x0, yb, x1, yb, lwd = 0.5, col = "grey25", xpd = NA)
   segments(c(x0, x1), yb, c(x0, x1), yb - tick, lwd = 0.5, col = "grey25", xpd = NA)
   text((x0 + x1)/2, yb + 0.34, lab, cex = BODY, col = "grey20",
        family = FAM, xpd = NA) }
brk(0.5, xs - GAPB, "training")
brk(xs + GAPB, nq + 0.5, "test")
segments(xs, 0, xs, yb, lty = 2, lwd = 0.8, col = "grey20", xpd = NA)

## ---- (a) colour bar -------------------------------------------------------
#  The bar needs its tick labels and title BELOW it, so the figure region runs
#  down to the foot of the page and a real bottom margin holds the text; with
#  mar = 0 the title was pushed off the canvas.
par(fig = c(CBX, CBF), mar = c(CB_MARLN, 0, 0, 0), new = TRUE)
zz <- seq(-CLIP, CLIP, length.out = 201)
image(x = zz, y = 1, z = matrix(zz, ncol = 1), col = PAL, zlim = c(-CLIP, CLIP),
      axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
box(col = "grey35", lwd = 0.4)
#  tick labels close under the bar, title below them, both inside the margin
axis(1, at = c(-3, 0, 3), cex.axis = BODY, lwd = 0, lwd.ticks = 0.4,
     mgp = c(1.4, -0.04, 0), tcl = -0.07, family = FAM)
mtext("standardized value (clipped at ±3)", side = 1, line = 0.52,
      cex = BODY, family = FAM)
par(mar = c(0, 0, 0, 0))

## ---- (b) the eigen spectrum ----------------------------------------------
C29 <- "#2E5E9E"; C100 <- "#B2182B"
#  pos100 already holds only the 99 nonzero eigenvalues, so the curve stops at
#  index 99; the zero one is not plotted and needs no marker or axis break.
#  The range spans the plotted values with a small margin either side.
posmin <- min(pos100); posmax <- max(ev100)
ybot   <- posmin / 2.0; ymax <- posmax * 2.0
TICKS  <- c(10, 1, 1e-2, 1e-4, 1e-6)
TLAB   <- expression(10^1, 10^0, 10^-2, 10^-4, 10^-6)
stopifnot(min(TICKS) >= ybot, max(TICKS) <= ymax)   # ticks cover the full range
par(fig = c(BX, BY), new = TRUE)
plot(NA, xlim = c(-2, 104), ylim = c(ybot, ymax), log = "y", axes = FALSE,
     xlab = "", ylab = "")
abline(h = ref(29),  lty = 2, lwd = 0.8, col = C29)
abline(h = ref(100), lty = 2, lwd = 0.8, col = C100)
lines(seq_along(ev29),   ev29,   col = C29,  lwd = 1.2)
lines(seq_along(pos100), pos100, col = C100, lwd = 1.2)
axis(1, at = c(1, 25, 50, 75, 100), cex.axis = BODY, lwd = 0, lwd.ticks = 0.4,
     mgp = c(1.4, 0.22, 0), tcl = -0.14, family = FAM)
axis(2, at = TICKS, labels = TLAB, las = 1, cex.axis = BODY, lwd = 0,
     lwd.ticks = 0.4, mgp = c(1.4, 0.32, 0), tcl = -0.14, family = FAM)
mtext("eigenvalue index", side = 1, line = 1.12, cex = BODY, family = FAM)
mtext("eigenvalue (log scale)", side = 2, line = 1.55, cex = BODY, family = FAM)
bu <- par("usr")
legend("bottomleft", inset = c(0.02, 0.17), bty = "n", cex = BODY,
       legend = c("p = 29", "p = 100", "noise bound, p = 29", "noise bound, p = 100"),
       col = c(C29, C100, C29, C100), lty = c(1, 1, 2, 2),
       lwd = c(1.2, 1.2, 0.8, 0.8), seg.len = 1.4, y.intersp = 0.98,
       x.intersp = 0.6, text.font = 1)
box(col = "grey35", lwd = 0.4)

#  ---- check: does the legend box cover either curve? ----------------------
lg <- legend("bottomleft", inset = c(0.02, 0.17), bty = "n", cex = BODY,
       legend = c("p = 29", "p = 100", "noise bound, p = 29", "noise bound, p = 100"),
       col = c(C29, C100, C29, C100), lty = c(1, 1, 2, 2),
       lwd = c(1.2, 1.2, 0.8, 0.8), seg.len = 1.4, y.intersp = 0.98,
       x.intersp = 0.6, plot = FALSE)
inbox <- function(x, y) x >= lg$rect$left & x <= lg$rect$left + lg$rect$w &
                        y <= lg$rect$top  & y >= lg$rect$top - lg$rect$h
hit29  <- sum(inbox(seq_along(ev29), ev29))
hit100 <- sum(inbox(seq_along(pos100), pos100))
cat(sprintf("\n=== (b) legend checks ===\n  legend box %.1f wide of %.1f panel units\n",
            lg$rect$w, bu[2] - bu[1]))
cat(sprintf("  curve points inside the legend box: p=29 %d | p=100 %d\n", hit29, hit100))
cat(sprintf("  legend right edge %.1f, panel right edge %.1f -> fits: %s\n",
            lg$rect$left + lg$rect$w, bu[2], lg$rect$left + lg$rect$w <= bu[2]))

## ---- panel tags -----------------------------------------------------------
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
TAGY <- 0.968
text(HX[1] - 0.038, TAGY, "(a)", font = 2, cex = TAG, adj = 0, family = FAM, xpd = NA)
text(BX[1] - 0.072, TAGY, "(b)", font = 2, cex = TAG, adj = 0, family = FAM, xpd = NA)
dev.off(); cat(sprintf("wrote %s\n", O3))

#  ---- text size, expected against rendered -------------------------------
#  Every text call in this block uses cex = BODY or cex = TAG against a
#  pointsize-12 device, so cex is a true point size.  A capital has a cap
#  height of about 0.717 em in this sans face.
cap_px <- function(pt) pt * 0.717 / 72 * DPI
TS <- data.frame(
   element = c("all body text (cex = BODY)", "panel tags (cex = TAG)"),
   cex = c(BODY, TAG), pointsize = PS3,
   pt = c(BODY, TAG) * PS3,
   expected_cap_px = round(cap_px(c(BODY, TAG) * PS3), 1),
   expected_cap_mm = round(c(BODY, TAG) * PS3 * 0.717 / 72 * 25.4, 2),
   row.names = NULL)
cat("\n=== (5) text size, expected at 6.5 in printed width and", DPI, "dpi ===\n")
print(TS, row.names = FALSE) }

dims <- function(f){ d <- readBin(f, "raw", 33)
   c(w = sum(as.integer(d[17:20]) * 256^(3:0)),
     h = sum(as.integer(d[21:24]) * 256^(3:0))) }
ff <- c(O1, O2, O3)[file.exists(c(O1, O2, O3))]
cat("\n=== files on disk ===\n")
print(data.frame(file = basename(ff), bytes = file.size(ff),
                 width_px  = vapply(ff, function(f) dims(f)[["w"]], numeric(1)),
                 height_px = vapply(ff, function(f) dims(f)[["h"]], numeric(1)),
                 row.names = NULL))
