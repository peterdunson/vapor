# ============================================================================
#  p = 29 partial-correlation heatmaps.  Reads p29_results.rds; refits nothing.
#
#  Boxes appear on the VAPOR panel only: the 95% credible-interval rule needs a
#  posterior and only VAPOR supplies one.  Each box is named in the panel's own
#  side margin, on the side its box sits nearest, with a short leader.
#
#  The figure is drawn at its PRINTED size, 6.5 in wide, so every cex below is
#  a true point size (cex = pt/12 at the default 12 pt device pointsize) and
#  nothing is scaled by the figure width.  The layout is laid out in inches and
#  each panel is placed with par(fig=), so all four matrices are exactly square
#  and exactly the same size.
# ============================================================================
RES_DIR  <- "/Users/peterdunson/Desktop/vapor/applications/results"
FIG_DIR  <- file.path(RES_DIR, "figures")
FITS_DIR <- file.path(RES_DIR, "fits")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
APP <- FITS_DIR            # fits are read from here
fit <- readRDS(file.path(APP, "p29_results.rds"))
OUT <- file.path(FIG_DIR, "heatmap_p29_r3k10.png")
vars <- fit$vars; n <- length(vars)
ed   <- fit$long$ci_edges
ei   <- cbind(match(ed$v1, vars), match(ed$v2, vars))
up   <- t(apply(ei, 1, function(z) c(min(z), max(z))))      # upper triangle
ord  <- order(-abs(ed$pcor))                                 # number by strength
up   <- up[ord, , drop = FALSE]; ed <- ed[ord, , drop = FALSE]
cat(sprintf("%d reported edges\n", nrow(up)))

avg_pc <- function(m) Reduce(`+`, m$pcor)/length(m$pcor)
have_kast <- length(fit$methods[["Kastner (VAR+FSV)"]]$pcor) > 0
kast_pc <- if (have_kast){
   cat(sprintf("Kastner panel: mean of %d chain(s)\n",
               length(fit$methods[["Kastner (VAR+FSV)"]]$pcor)))
   avg_pc(fit$methods[["Kastner (VAR+FSV)"]])
} else { cat("Kastner panel OMITTED: no Kastner chain completed\n"); NULL }

nodiag <- function(M){ diag(M) <- NA; M }
panels <- c(list("VAPOR" = nodiag(fit$long$pc_vapor)),
            if (!is.null(kast_pc)) list("Kastner" = nodiag(kast_pc)),
            list("VAR + Glasso"      = nodiag(fit$methods[["VAR + Glasso"]]$pcor[[1]]),
                 "Sample covariance" = nodiag(fit$methods[["Sample covariance"]]$pcor[[1]])))
NP  <- length(panels)

# ----------------------------------------------------------------------------
#  LAYOUT, in inches at the printed size
# ----------------------------------------------------------------------------
W_IN      <- 6.5
PS        <- 12                        # device pointsize
pt        <- function(x) x/PS          # cex for a true point size
TITLE_PT  <- 9; TICK_PT <- 6; BARLAB_PT <- 6; EDGE_PT <- 5.5
FAM       <- "serif"
CSI       <- PS/72                     # one line of text, inches

LM <- 0.03; RM <- 0.03                 # equal panels, evenly spaced
GAP       <- 0.06                      # between panels
GAP_BAR   <- 0.22                      # before the colour bar
BAR_W     <- 0.12                      # colour bar width
TICK_W    <- 0.24                      # room for the bar's tick labels
PAD_R     <- 0.05                      # clearance at the right edge
TITLE_H   <- 0.20                      # band above each panel for its title
PAD_T     <- 0.13                      # clearance above the titles / bar label
PAD_B     <- 0.05

SQ <- (W_IN - (LM + RM + (NP-1)*GAP + GAP_BAR + BAR_W + TICK_W + PAD_R))/NP
H_IN <- PAD_T + TITLE_H + SQ + PAD_B
cat(sprintf("printed size %.2f x %.2f in; square %.3f in; cell %.4f in\n",
            W_IN, H_IN, SQ, SQ/n))

#  left edge of each panel
px <- numeric(NP); px[1] <- LM
if (NP > 1) for (q in 2:NP) px[q] <- px[q-1] + SQ + GAP
bar_x <- px[NP] + SQ + GAP_BAR
y0 <- PAD_B; y1 <- y0 + SQ            # panels occupy this vertical band

pal    <- colorRampPalette(c("blue","white","red"))(101)
BOX_HW <- 0.5
LAB    <- paste(ed$v1, ed$v2, sep = "-")

# ----------------------------------------------------------------------------
#  LABELS.  Each box is named INSIDE the VAPOR panel, next to its own box, by
#  a search over offsets around the box: shortest leader first, and among
#  those, the position whose cells are blankest.  A candidate is rejected if
#  its rectangle would leave the panel, cover any box, or touch another label,
#  so no label can overlap a box or another label and none leaves the panel.
# ----------------------------------------------------------------------------
LAB_PT   <- 5.0                        # label text, points at the printed size
LEAD_LWD <- 0.8/0.75                   # 0.8 pt: lwd 1 is 1/96 in = 0.75 pt
OFFS <- local({
   o <- list()
   for (r in c(2.2, 2.8, 3.4, 4.2, 5.0, 6.0, 7.2, 8.6, 10.2, 12.0))
      for (a in seq(0, 345, by = 15))
         o[[length(o)+1L]] <- c(r*cos(a*pi/180), r*sin(a*pi/180))
   o })

seg_hits_rect <- function(p1, p2, r){
   if (max(p1[1],p2[1]) < r[1] || min(p1[1],p2[1]) > r[2] ||
       max(p1[2],p2[2]) < r[3] || min(p1[2],p2[2]) > r[4]) return(FALSE)
   ins <- function(p) p[1] >= r[1] && p[1] <= r[2] && p[2] >= r[3] && p[2] <= r[4]
   if (ins(p1) || ins(p2)) return(TRUE)
   ccw <- function(a,b,c) (c[2]-a[2])*(b[1]-a[1]) > (b[2]-a[2])*(c[1]-a[1])
   si <- function(a,b,c,d) (ccw(a,c,d) != ccw(b,c,d)) && (ccw(a,b,c) != ccw(a,b,d))
   cn <- list(c(r[1],r[3]), c(r[2],r[3]), c(r[2],r[4]), c(r[1],r[4]))
   for (k in 1:4) if (si(p1, p2, cn[[k]], cn[[k %% 4 + 1]])) return(TRUE)
   FALSE }

#  do two open segments cross?  Used to keep leaders from crossing each other.
seg_hits_seg <- function(a1, a2, b1, b2){
   ccw <- function(a,b,c) (c[2]-a[2])*(b[1]-a[1]) - (b[2]-a[2])*(c[1]-a[1])
   d1 <- ccw(b1,b2,a1); d2 <- ccw(b1,b2,a2)
   d3 <- ccw(a1,a2,b1); d4 <- ccw(a1,a2,b2)
   ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) }

#  NUMBER KEY.  The 13 labels are far too wide to sit on a 29-unit panel
#  without collisions, so each boxed cell carries a small bold number instead,
#  in the order of the edge table (descending |pooled mean|).  The key is
#  printed to the console for the caption.  Placement is a short offset search:
#  nearest position whose number clears every box, every number already placed
#  and the panel border.
NUM_PT <- 4.2                          # number text, points at the printed size
NOFF <- local({
   o <- list()
   for (r in c(1.15, 1.5, 1.9, 2.4, 3.0, 3.8))
      for (ang in seq(0, 345, by = 15))
         o[[length(o)+1L]] <- c(r*cos(ang*pi/180), r*sin(ang*pi/180))
   o })

place_numbers <- function(up, n, nums){
   boxes <- cbind(up[,2], n - up[,1] + 1)
   hw <- strwidth(nums, cex = pt(NUM_PT), font = 2, units = "user")/2
   hh <- strheight("8", cex = pt(NUM_PT), font = 2, units = "user")/2
   out <- matrix(NA_real_, nrow(up), 2); nfail <- 0L
   for (e in seq_len(nrow(up))){
      cx <- boxes[e,1]; cy <- boxes[e,2]; chosen <- NULL; best <- Inf
      for (d in NOFF){
         px <- cx + d[1]; py <- cy + d[2]
         if (px - hw[e] < 0.7 || px + hw[e] > n + 0.3) next
         if (py - hh < 0.7 || py + hh > n + 0.3) next
         r <- c(px-hw[e]-0.28, px+hw[e]+0.28, py-hh-0.28, py+hh+0.28)
         if (any(boxes[,1] + BOX_HW > r[1] & boxes[,1] - BOX_HW < r[2] &
                 boxes[,2] + BOX_HW > r[3] & boxes[,2] - BOX_HW < r[4])) next
         bad <- FALSE
         for (k in which(!is.na(out[,1])))
            if (abs(out[k,1]-px) < (hw[k]+hw[e]+0.45) &&
                abs(out[k,2]-py) < (2*hh+0.45)){ bad <- TRUE; break }
         if (bad) next
         sc <- sqrt(d[1]^2 + d[2]^2)
         if (sc < best){ best <- sc; chosen <- c(px, py) } }
      if (is.null(chosen)){ nfail <- nfail + 1L
         chosen <- c(min(n - hw[e], cx + 1.3), cy) }
      out[e,] <- chosen }
   #  report the realised geometry
   nov <- 0L
   for (e in seq_len(nrow(up))) for (k in seq_len(nrow(up))) if (k > e)
      if (abs(out[k,1]-out[e,1]) < (hw[k]+hw[e]) &&
          abs(out[k,2]-out[e,2]) < 2*hh) nov <- nov + 1L
   nbox <- 0L
   for (e in seq_len(nrow(up))){
      r <- c(out[e,1]-hw[e], out[e,1]+hw[e], out[e,2]-hh, out[e,2]+hh)
      if (any(boxes[,1] + BOX_HW > r[1] & boxes[,1] - BOX_HW < r[2] &
              boxes[,2] + BOX_HW > r[3] & boxes[,2] - BOX_HW < r[4])) nbox <- nbox + 1L }
   cat(sprintf("number key: %d numbers | number-number overlaps: %d | numbers touching a box: %d | unplaced by the search: %d\n",
               nrow(up), nov, nbox, nfail))
   list(pos = out, hw = hw, hh = hh) }

png(OUT, width = W_IN, height = H_IN, units = "in", res = 600, bg = "white",
    pointsize = PS)
par(family = FAM)
nms <- names(panels)
P <- NULL

for (q in seq_along(panels)){
   par(fig = c(px[q], px[q] + SQ, y0, y1)/c(W_IN, W_IN, H_IN, H_IN),
       mar = c(0,0,0,0), new = (q > 1), xpd = NA)
   M <- panels[[q]]; M[M > 1] <- 1; M[M < -1] <- -1
   image(x = 1:n, y = 1:n, z = t(M[n:1, ]), col = pal, zlim = c(-1,1),
         axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
   box(col = "grey60", lwd = 0.5)                       # thin grey panel border
   #  titles all sit at the same height, measured from the top of the panels
   text(mean(par("usr")[1:2]), n + 0.5 + (TITLE_H*0.52)/(SQ/n),
        nms[q], cex = pt(TITLE_PT), font = 2, family = FAM, adj = c(0.5, 0))
   if (nms[q] != "VAPOR") next
   NUMS <- as.character(seq_len(nrow(up)))
   P <- place_numbers(up, n, NUMS)
   boxes <- cbind(up[,2], n - up[,1] + 1)
   for (e in seq_len(nrow(up))){
      rect(boxes[e,1]-BOX_HW, boxes[e,2]-BOX_HW,
           boxes[e,1]+BOX_HW, boxes[e,2]+BOX_HW, border = "black", lwd = 0.6)
      text(P$pos[e,1], P$pos[e,2], NUMS[e], cex = pt(NUM_PT), font = 2,
           family = FAM, col = "black") } }

# ---- colour bar: same height as the panels, narrow ------------------------
par(fig = c(bar_x, bar_x + BAR_W, y0, y1)/c(W_IN, W_IN, H_IN, H_IN),
    mar = c(0,0,0,0), new = TRUE, xpd = NA)
image(x = 1, y = seq(-1, 1, length.out = 101),
      z = matrix(seq(-1, 1, length.out = 101), nrow = 1), col = pal,
      axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
box(col = "grey60", lwd = 0.5)
axis(4, at = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5)),
     las = 1, cex.axis = pt(TICK_PT), tick = FALSE, line = -0.62,
     family = FAM, col.axis = "black")
#  label above the bar, normal weight, normal line spacing
lh <- (BARLAB_PT*1.18/72)/(SQ/2/1)          # one line, in user units of the bar
text(1, 1 + 0.5*lh, "Correlation", cex = pt(BARLAB_PT), family = FAM, adj = c(0.5, 0))
text(1, 1 + 1.5*lh, "Partial",     cex = pt(BARLAB_PT), family = FAM, adj = c(0.5, 0))
dev.off()
cat(sprintf("wrote %s (%d panels, %d labelled edges)\n", OUT, NP, nrow(up)))

#  caption key, LaTeX-ready
ed_ord <- fit$long$ci_edges
#  inside $...$ a leading "-" already typesets as a minus sign
#  the key names series by their short readable label, from the one shared
#  file, so the caption carries no FRED-QD mnemonic either
KEY_SHORT <- local({
   d <- read.csv("/Users/peterdunson/Desktop/vapor/applications/series_labels.csv",
                 stringsAsFactors = FALSE)
   setNames(d$short_label, d$mnemonic) })
stopifnot(all(c(ed_ord$v1, ed_ord$v2) %in% names(KEY_SHORT)))
tex <- function(x) gsub("&", "\\\\&", x)      # the only LaTeX special in the labels
KEYTXT <- paste(sprintf("%d %s--%s ($%.2f$)", seq_len(nrow(ed_ord)),
                        tex(KEY_SHORT[ed_ord$v1]), tex(KEY_SHORT[ed_ord$v2]),
                        ed_ord$pcor), collapse = ", ")
cat("\n--- caption key, LaTeX-ready ---\n")
cat(KEYTXT, "\n")
