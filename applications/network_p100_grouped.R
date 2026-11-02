# ============================================================================
#  p = 100 conditional dependence graph, laid out BY CONNECTED COMPONENT.
#  Reads p100_results.rds; refits nothing.
#
#  The old version placed every node by hand, which broke as soon as the fit
#  reported a series the layout did not know.  This version computes the layout:
#
#    * each connected component is laid out on its own, compactly, with a
#      circle for components of 3+ nodes and a short bar for pairs;
#    * components are packed into a grid ordered by the official FRED-QD group
#      of their members, each inside its own cell, so no edge can reach across
#      another component or pass through a node that is not its endpoint;
#    * edge labels sit at the midpoint of the edge on a white background, and
#      are suppressed only if the midpoint would land on a node;
#    * node colours and legend keys are generated from ONE palette indexed by
#      group name, so the legend cannot disagree with the nodes.
# ============================================================================
suppressMessages(library(igraph))
RES_DIR  <- "/Users/peterdunson/Desktop/vapor/applications/results"
FIG_DIR  <- file.path(RES_DIR, "figures")
FITS_DIR <- file.path(RES_DIR, "fits")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
APP <- FITS_DIR            # fits are read from here
fit <- readRDS(file.path(APP, "p100_results.rds"))
OUT <- file.path(FIG_DIR, "network_p100_grouped.png")
GROUPS <- "/Users/peterdunson/Desktop/vapor/applications/fredqd_groups.csv"
e <- fit$long$ci_edges
cat(sprintf("%d reported edges\n", nrow(e)))

#  official FRED-QD groups
g <- read.csv(GROUPS, stringsAsFactors = FALSE)
gname <- setNames(g$group_name, g$series); gcode <- setNames(g$group_code, g$series)
NODES <- sort(unique(c(e$v1, e$v2)))
miss <- setdiff(NODES, names(gname))
if (length(miss)) stop("no official group for: ", paste(miss, collapse = ", "))
sect <- gname[NODES]; scode <- gcode[NODES]

gr <- graph_from_data_frame(data.frame(from = e$v1, to = e$v2, pcor = e$pcor),
                            directed = FALSE,
                            vertices = data.frame(name = NODES, sector = unname(sect)))

#  ---- one palette, used for both the nodes and the legend ------------------
lv <- unique(data.frame(code = unname(scode), name = unname(sect),
                        stringsAsFactors = FALSE))
lv <- lv[order(lv$code), ]
PALETTE <- c("#8FB8DE","#9FD3A8","#F2C14E","#E8927C","#C3A6D8",
             "#7FC9C4","#F2A6C3","#BFD56B","#D9B38C","#A8A8D8")
pal <- setNames(PALETTE[seq_len(nrow(lv))], lv$name)
V(gr)$sector <- sect[V(gr)$name]
V(gr)$color  <- unname(pal[V(gr)$sector])

#  ---- per-component layout -------------------------------------------------
em <- ends(gr, E(gr))
comp <- components(gr)
cl <- lapply(seq_len(comp$no), function(k) names(which(comp$membership == k)))
#  order components by the group of their members, then by size
ckey <- vapply(cl, function(v) min(scode[v]), numeric(1))
csz  <- vapply(cl, length, integer(1))
ord  <- order(ckey, -csz)
cl <- cl[ord]
cat(sprintf("%d connected components, sizes %s\n", length(cl),
            paste(vapply(cl, length, integer(1)), collapse = ", ")))

#  ---- GEOMETRY -------------------------------------------------------------
#  Everything below is computed in INCHES and drawn with base graphics, with
#  the plotting window set so that one user unit is one inch.  That makes
#  strwidth() directly comparable to node radii and spacings, so label
#  collisions can be resolved exactly rather than by trial and error.
#
#  Each node's label is placed on the far side of the node from its component's
#  centre, which keeps every label out of the interior where the edges and the
#  edge values sit.  The NIPA component is given an explicit planar embedding
#  (GDPC1 in the middle, EXPGSC1 and IMPGSC1 opposite, GPDIC1 beyond them) so
#  none of its seven edges cross, plus an explicit label direction for the
#  interior node.
W_IN <- 6.6; H_LEGEND <- 0.62
R_NODE <- 0.055                 # node radius, inches
PAD    <- 0.035                 # gap between node edge and its label
LCEX   <- 0.68                  # ~8.2 pt at the default 12 pt pointsize
ECEX   <- 0.62
ROW_MAX <- 6.45                 # usable width for one row of components
GAP_X  <- 0.13; GAP_Y <- 0.13

PLANAR <- list(nodes = c("GDPC1","GPDIC1","EXPGSC1","IMPGSC1","GCEC1"),
   #  this component carries seven edges in a small area, so it is drawn
   #  larger than the others: otherwise its edge values have nowhere to sit
   xy = rbind(GDPC1 = c( 0.00,  0.00), GPDIC1 = c( 0.00, -0.64),
              EXPGSC1 = c(-0.64, -0.32), IMPGSC1 = c( 0.64, -0.32),
              GCEC1 = c( 0.00,  0.50)),
   dir = rbind(GDPC1 = c(-1, 0.35), GPDIC1 = c(0, -1), EXPGSC1 = c(-1, 0),
               IMPGSC1 = c(1, 0), GCEC1 = c(0, 1)))

#  does the segment p1-p2 touch the rectangle r = c(x0, x1, y0, y1)?
seg_hits_rect <- function(p1, p2, r){
   if (max(p1[1],p2[1]) < r[1] || min(p1[1],p2[1]) > r[2] ||
       max(p1[2],p2[2]) < r[3] || min(p1[2],p2[2]) > r[4]) return(FALSE)
   inside <- function(p) p[1] >= r[1] && p[1] <= r[2] && p[2] >= r[3] && p[2] <= r[4]
   if (inside(p1) || inside(p2)) return(TRUE)
   ccw <- function(a,b,c) (c[2]-a[2])*(b[1]-a[1]) > (b[2]-a[2])*(c[1]-a[1])
   segint <- function(a,b,c,d) (ccw(a,c,d) != ccw(b,c,d)) && (ccw(a,b,c) != ccw(a,b,d))
   cn <- list(c(r[1],r[3]), c(r[2],r[3]), c(r[2],r[4]), c(r[1],r[4]))
   for (k in 1:4) if (segint(p1, p2, cn[[k]], cn[[k %% 4 + 1]])) return(TRUE)
   FALSE }

unit_layout <- function(v){
   if (setequal(v, PLANAR$nodes))
      return(list(xy = PLANAR$xy[v, , drop = FALSE], dir = PLANAR$dir[v, , drop = FALSE]))
   m <- length(v)
   #  horizontal separation is driven by the LABEL widths, so two labels that
   #  sit under neighbouring nodes cannot run into each other
   sep2 <- function(a, b) max(0.50, (lw[a] + lw[b])/2 + 0.09)/2
   #  A two-node component is drawn VERTICALLY with each name to the right of
   #  its own node.  Side by side, two ~0.6 in labels force the component to
   #  1.2 in; stacked this way it is about 0.8 in, which is what lets the
   #  figure be landscape.
   xy <- if (m == 2) cbind(c(0, 0), c(0.25, -0.25))
         else if (m == 3){ h <- sep2(v[1], v[2]); cbind(c(-h, h, 0.00), c(-0.17, -0.17, 0.30)) }
         else { th <- seq(0, 2*pi, length.out = m+1)[-(m+1)] + pi/2
                #  neighbours on the circle are 2 r sin(pi/m) apart, so the
                #  radius has to grow with the label width or two labels on
                #  the same side run into each other
                r <- max(0.42, (max(lw[v]) + 0.10) / (2.6 * sin(pi/m)))
                cbind(r*cos(th), r*sin(th)) }
   rownames(xy) <- v
   #  a label goes BELOW its node unless the node has a neighbour below it, in
   #  which case it goes above.  Vertical placement keeps components narrow,
   #  which is what lets the figure fit the text width.
   #  On a ring (m >= 4) the label goes RADIALLY outward, away from the
   #  component's centre, which is what keeps labels off the interior and off
   #  each other.  Pushing every label straight up or down instead stacked
   #  the two left-hand rate labels on top of one another.
   dr <- if (m == 2) rbind(c(1, 0), c(1, 0))
      else if (m >= 4) t(apply(xy, 1, function(q){
         n <- sqrt(sum(q^2)); if (n < 1e-9) c(0, -1) else q/n }))
      else t(vapply(v, function(u){
         nb <- intersect(unique(c(em[em[,1]==u,2], em[em[,2]==u,1])), v)
         below <- length(nb) && any(xy[nb,2] < xy[u,2] - 1e-9)
         if (below) c(0, 1) else c(0, -1) }, numeric(2)))
   rownames(dr) <- v
   list(xy = xy, dir = dr) }

#  ---------------------------------------------------------------------------
#  Node labels: a short readable name instead of the FRED-QD mnemonic.  Each is
#  a shortening of that series' description in varlist_descriptions.csv, kept
#  to about 20 characters so the layout, which is driven by label widths, still
#  resolves.  The graph's vertex names stay the mnemonics.
#  ---------------------------------------------------------------------------
SHORT <- c(
  AAA = "Aaa corp. bond",            AMDMNOx = "Orders, durables",
  ANDENOx = "Orders, capital gds",   BAA = "Baa corp. bond",
  CPIAUCSL = "CPI, all items",       CPITRNSL = "CPI, transport",
  CUMFNS = "Capacity util., mfg",    DMANEMP = "Emp., durable mfg",
  FEDFUNDS = "Fed funds",            GS1 = "1-yr Treasury",
  GS10 = "10-yr Treasury",           GS5 = "5-yr Treasury",
  HOUST = "Housing starts",          HOUSTS = "Housing starts, S.",
  HOUSTW = "Housing starts, W.",     INDPRO = "Industrial prod.",
  IPBUSEQ = "IP, business eq.",      IPCONGD = "IP, consumer gds",
  IPFINAL = "IP, final products",    IPMAT = "IP, materials",
  M1REAL = "Real M1",                M2REAL = "Real M2",
  MANEMP = "Emp., manufacturing",    MORTGAGE30US = "30-yr mortgage",
  PCDGx = "PCE, durables",           PCECC96 = "PCE, total",
  PCECTPI = "PCE price index",       PCESVx = "PCE, services",
  PPIACO = "PPI, all commod.",       TB3MS = "3-mo T-bill",
  TCU = "Capacity util., total",     UEMP27OV = "Unemp. 27+ weeks",
  UNRATE = "Unemployment rate",      UNRATELTx = "Unemp. rate, 15+ wk",
  UNRATESTx = "Short-term unemp.",   USCONS = "Emp., construction",
  USGOOD = "Emp., goods-prod.",      USLAH = "Emp., leisure/hosp.",
  USPRIV = "Emp., total private",    WPSID61 = "PPI, intermediate")
miss2 <- setdiff(NODES, names(SHORT))
if (length(miss2)) stop(sprintf("no short label for: %s", paste(miss2, collapse = ", ")))
DISP <- SHORT[NODES]
cat(sprintf("node labels: %d, longest %d characters (%s)\n", length(DISP),
            max(nchar(DISP)), names(which.max(nchar(DISP)))))
print(data.frame(mnemonic = NODES, label = unname(DISP),
                 nchar = nchar(DISP), row.names = NULL))

png(OUT, width = W_IN, height = 7.4, units = "in", res = 600, bg = "white")
plot.new()                                   # so strwidth() is available
lw <- setNames(strwidth(DISP, cex = LCEX, units = "inches"), NODES)
lh <- strheight("X", cex = LCEX, units = "inches")

#  per-component geometry, in inches, centred on (0,0)
GEO <- lapply(cl, function(v){
   u <- unit_layout(v)
   xy <- u$xy; dr <- u$dir
   lab <- matrix(NA_real_, length(v), 2, dimnames = list(v, NULL))
   adjx <- adjy <- numeric(length(v))
   for (i in seq_along(v)){
      d <- dr[i,]; if (sqrt(sum(d^2)) < 1e-9) d <- c(0,-1)
      d <- d/sqrt(sum(d^2))
      lab[i,] <- xy[i,] + d*(R_NODE + PAD + 0.5*lh*abs(d[2]))
      adjx[i] <- if (abs(d[1]) < 0.35) 0.5 else if (d[1] > 0) 0 else 1
      adjy[i] <- if (abs(d[2]) < 0.35) 0.5 else if (d[2] > 0) 0 else 1 }
   #  push a label out until it clears every edge of its own component: this is
   #  what was putting PCECC96 across its edge and USGOOD on the DMANEMP edge
   ein0 <- which(em[,1] %in% v & em[,2] %in% v)
   for (i in seq_along(v)){
      guard <- 0
      repeat {
         r <- c(lab[i,1] - lw[v[i]]*adjx[i] - 0.012, lab[i,1] + lw[v[i]]*(1-adjx[i]) + 0.012,
                lab[i,2] - lh*adjy[i] - 0.008,       lab[i,2] + lh*(1-adjy[i]) + 0.008)
         hit <- FALSE
         for (ii in ein0)
            if (seg_hits_rect(xy[em[ii,1],], xy[em[ii,2],], r)){ hit <- TRUE; break }
         guard <- guard + 1
         if (!hit || guard > 12) break
         lab[i,] <- lab[i,] + dr[i,]*(lh*0.55) } }

   #  if two labels in this component would touch, push the later one further
   #  out along its own direction until it clears
   ord <- order(lab[,2], decreasing = TRUE)
   for (a in seq_along(ord)) for (b in seq_len(a-1)){
      i <- ord[a]; j <- ord[b]
      repeat {
         xi <- c(lab[i,1] - lw[v[i]]*adjx[i], lab[i,1] + lw[v[i]]*(1-adjx[i]))
         xj <- c(lab[j,1] - lw[v[j]]*adjx[j], lab[j,1] + lw[v[j]]*(1-adjx[j]))
         yi <- c(lab[i,2] - lh*adjy[i], lab[i,2] + lh*(1-adjy[i]))
         yj <- c(lab[j,2] - lh*adjy[j], lab[j,2] + lh*(1-adjy[j]))
         if (xi[2] < xj[1]-0.02 || xj[2] < xi[1]-0.02 ||
             yi[2] < yj[1]-0.01 || yj[2] < yi[1]-0.01) break
         d <- dr[i,]; lab[i,] <- lab[i,] + d*(lh*1.12) } }
   #  bounding box over nodes and their labels
   x0 <- min(xy[,1]-R_NODE, lab[,1] - lw[v]*adjx)
   x1 <- max(xy[,1]+R_NODE, lab[,1] + lw[v]*(1-adjx))
   y0 <- min(xy[,2]-R_NODE, lab[,2] - lh*adjy)
   y1 <- max(xy[,2]+R_NODE, lab[,2] + lh*(1-adjy))
   list(v = v, xy = xy, lab = lab, adjx = adjx, adjy = adjy,
        w = x1-x0, h = y1-y0, cx = (x0+x1)/2, cy = (y0+y1)/2) })

#  Pack components into rows, balanced by width but NEVER wider than ROW_MAX.
#  The cap matters: the plot is drawn with asp = 1, so a row wider than the
#  panel shrinks the whole layout while the labels keep their point size, and
#  labels that were computed to clear each other end up overlapping.
NROWS <- 3L
wv <- vapply(GEO, function(g) g$w, numeric(1))
target <- min(ROW_MAX, (sum(wv) + GAP_X*(length(wv) - NROWS))/NROWS)
rows <- list(); cur <- integer(0); w <- 0
for (k in seq_along(GEO)){
   add <- wv[k] + if (length(cur)) GAP_X else 0
   over_cap <- (w + add) > ROW_MAX
   #  only split a row once it is genuinely full: splitting on the share
   #  estimate alone leaves a nearly empty last row
   over_share <- length(rows) < NROWS - 1L && w + add*0.92 > target
   if (length(cur) && (over_cap || over_share)){
      rows[[length(rows)+1L]] <- cur; cur <- integer(0); w <- 0; add <- wv[k] }
   cur <- c(cur, k); w <- w + add }
if (length(cur)) rows[[length(rows)+1L]] <- cur
#  Greedy packing cannot look back, so a component that arrives last can be
#  stranded on a row of its own even when an earlier row has room.  Move any
#  such stragglers forward into the first row that can take them.
row_w <- function(ks) if (!length(ks)) 0 else sum(wv[ks]) + GAP_X*(length(ks)-1)
repeat {
   if (length(rows) < 2) break
   last <- rows[[length(rows)]]
   if (row_w(last) > ROW_MAX/2) break          # a genuinely full row, leave it
   moved <- FALSE
   for (k in last){
      for (r in seq_len(length(rows)-1L))
         if (row_w(c(rows[[r]], k)) <= ROW_MAX){
            rows[[r]] <- c(rows[[r]], k)
            rows[[length(rows)]] <- setdiff(rows[[length(rows)]], k)
            moved <- TRUE; break }
      if (moved) break }
   if (!moved) break
   if (!length(rows[[length(rows)]])) rows[[length(rows)]] <- NULL }
rw <- vapply(rows, function(ks) sum(wv[ks]) + GAP_X*(length(ks)-1), numeric(1))
cat(sprintf("component widths: %s\n", paste(sprintf("%s=%.2f", vapply(cl, function(v) v[1], character(1)), wv), collapse=", ")))
cat(sprintf("row widths (in): %s | cap %.2f\n",
            paste(sprintf("%.2f", rw), collapse = ", "), ROW_MAX))

OFF <- matrix(NA_real_, length(GEO), 2)
ytop <- 0
for (r in seq_along(rows)){
   ks <- rows[[r]]; ws <- vapply(GEO[ks], function(g) g$w, numeric(1))
   hs <- vapply(GEO[ks], function(g) g$h, numeric(1))
   total <- sum(ws) + GAP_X*(length(ks)-1)
   x <- -total/2
   for (i in seq_along(ks)){
      k <- ks[i]
      OFF[k,1] <- x + ws[i]/2 - GEO[[k]]$cx
      OFF[k,2] <- ytop - max(hs)/2 - GEO[[k]]$cy
      x <- x + ws[i] + GAP_X }
   ytop <- ytop - max(hs) - GAP_Y }
ytop <- ytop + GAP_Y          # no trailing gap under the last row
H_PLOT <- abs(ytop) + 0.30

NP <- do.call(rbind, lapply(seq_along(GEO), function(k){
   g <- GEO[[k]]; m <- sweep(g$xy, 2, OFF[k,], "+"); rownames(m) <- g$v; m }))
LP <- do.call(rbind, lapply(seq_along(GEO), function(k){
   g <- GEO[[k]]; m <- sweep(g$lab, 2, OFF[k,], "+"); rownames(m) <- g$v; m }))
ADJ <- do.call(rbind, lapply(GEO, function(g) cbind(g$adjx, g$adjy)))
rownames(ADJ) <- unlist(lapply(GEO, function(g) g$v))

#  ---- DRAW -----------------------------------------------------------------
dev.off()
H_IN <- H_PLOT + H_LEGEND + 0.12
png(OUT, width = W_IN, height = H_IN, units = "in", res = 600, bg = "white")
layout(matrix(1:2, ncol = 1), heights = c(H_PLOT + 0.12, H_LEGEND))
par(mar = c(0.05, 0.05, 0.05, 0.05))
#  The drawing range must use the ALIGNMENT-aware extent of each label: a
#  right-aligned name occupies [x, x+width], not [x-width, x+width].  Taking
#  the symmetric bound inflates the x range, and since the plot is drawn with
#  asp = 1 that silently compresses the whole layout, which is what made
#  labels computed to clear each other overlap on the page.
lwv <- lw[rownames(LP)]
xr <- range(c(NP[,1] - R_NODE, NP[,1] + R_NODE,
              LP[,1] - lwv*ADJ[,1], LP[,1] + lwv*(1 - ADJ[,1])))
yr <- range(c(NP[,2] - R_NODE, NP[,2] + R_NODE,
              LP[,2] - lh*ADJ[,2], LP[,2] + lh*(1 - ADJ[,2])))
plot(NA, xlim = c(xr[1]-0.05, xr[2]+0.05), ylim = c(yr[1]-0.04, yr[2]+0.04),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
cat(sprintf("plot span %.2f x %.2f in; scale %.3f in per unit\n",
            diff(xr)+0.10, diff(yr)+0.08,
            min((W_IN-0.10)/(diff(xr)+0.10), (H_PLOT+0.02)/(diff(yr)+0.08))))

circle <- function(x, y, r, col){
   th <- seq(0, 2*pi, length.out = 60)
   polygon(x + r*cos(th), y + r*sin(th), col = col, border = "grey30", lwd = 0.7) }

for (i in seq_len(nrow(em))){
   a <- NP[em[i,1], ]; b <- NP[em[i,2], ]
   segments(a[1], a[2], b[1], b[2],
            lwd = 0.8 + 2.4*abs(E(gr)$pcor[i]),
            col = if (E(gr)$pcor[i] > 0) "#2E5E9E" else "#B2182B") }
#  Edge values.  Each is tried at the midpoint first and then slid along its
#  own edge; a position is rejected if it would cover a node, a node label, or
#  a value already placed.  ECEX is about 7.4 pt at this width.
LABR <- do.call(rbind, lapply(rownames(LP), function(v)
   c(LP[v,1] - lw[v]*ADJ[v,1] - 0.010, LP[v,1] + lw[v]*(1-ADJ[v,1]) + 0.010,
     LP[v,2] - lh*ADJ[v,2] - 0.006,    LP[v,2] + lh*(1-ADJ[v,2]) + 0.006)))
TS <- c(0.50, 0.42, 0.58, 0.34, 0.66, 0.28, 0.72, 0.22, 0.78, 0.16, 0.84)
vrect <- NULL; VPOS <- matrix(NA_real_, nrow(em), 2); VLAB <- character(nrow(em))
for (i in seq_len(nrow(em))){
   a <- NP[em[i,1], ]; b <- NP[em[i,2], ]
   lab2 <- sprintf("%+.2f", E(gr)$pcor[i]); VLAB[i] <- lab2
   hw2 <- strwidth(lab2, cex = ECEX)/2 + 0.011
   hh2 <- strheight(lab2, cex = ECEX)/2 + 0.008
   d <- b - a; L <- sqrt(sum(d^2)); perp <- if (L > 0) c(-d[2], d[1])/L else c(0,1)
   cand <- list()
   #  slide ALONG the edge first and only then step off it: a value nudged
   #  perpendicular reads as belonging to whichever edge it lands near
   for (o in c(0, 0.055, -0.055, 0.10, -0.10, 0.145, -0.145, 0.19, -0.19)) for (t in TS)
      cand[[length(cand)+1L]] <- list(t = t, o = o)
   for (cc in cand){
      m <- a + cc$t*(b - a) + perp*cc$o
      r <- c(m[1]-hw2, m[1]+hw2, m[2]-hh2, m[2]+hh2)
      #  clear of every node disc
      if (any(NP[,1] > r[1]-R_NODE & NP[,1] < r[2]+R_NODE &
              NP[,2] > r[3]-R_NODE & NP[,2] < r[4]+R_NODE)){
         next }
      #  clear of every node label
      if (any(LABR[,1] < r[2] & LABR[,2] > r[1] & LABR[,3] < r[4] & LABR[,4] > r[3])){
         next }
      #  clear of every value already placed
      if (!is.null(vrect) &&
          any(vrect[,1] < r[2] & vrect[,2] > r[1] & vrect[,3] < r[4] & vrect[,4] > r[3])){
         next }
      #  a value sitting where another edge crosses reads as belonging to that
      #  edge, so reject any position a different edge passes through
      cross <- FALSE
      for (j in seq_len(nrow(em))) if (j != i)
         if (seg_hits_rect(NP[em[j,1],], NP[em[j,2],], r)){ cross <- TRUE; break }
      if (cross) next
      VPOS[i,] <- m; vrect <- rbind(vrect, r); break } }
nplaced <- sum(!is.na(VPOS[,1]))
cat(sprintf("edge values placed: %d of %d\n", nplaced, nrow(em)))
if (nplaced < nrow(em)) cat(sprintf("  unplaced: %s\n", paste(apply(em[is.na(VPOS[,1]), , drop=FALSE], 1, paste, collapse="--"), collapse=", ")))
for (i in seq_len(nrow(em))) if (!is.na(VPOS[i,1])){
   hw2 <- strwidth(VLAB[i], cex = ECEX)/2 + 0.011
   hh2 <- strheight(VLAB[i], cex = ECEX)/2 + 0.008
   rect(VPOS[i,1]-hw2, VPOS[i,2]-hh2, VPOS[i,1]+hw2, VPOS[i,2]+hh2,
        col = "white", border = NA)
   text(VPOS[i,1], VPOS[i,2], VLAB[i], cex = ECEX, col = "grey15") }

for (v in rownames(NP))
   circle(NP[v,1], NP[v,2], R_NODE, unname(pal[sect[v]]))
for (v in rownames(LP))
   text(LP[v,1], LP[v,2], DISP[v], cex = LCEX, col = "grey10",
        adj = c(ADJ[v,1], ADJ[v,2]))

par(mar = c(0.05, 0.2, 0.05, 0.2))
plot(NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
legend("top", ncol = 4, bty = "n", cex = 0.56, pt.cex = 1.15,
       legend = names(pal), pt.bg = unname(pal), pch = 21, col = "grey30",
       x.intersp = 0.9, y.intersp = 1.1)
legend("bottom", ncol = 2, bty = "n", cex = 0.56, seg.len = 1.8,
       legend = c("positive partial correlation", "negative partial correlation"),
       col = c("#2E5E9E", "#B2182B"), lwd = 2.4, x.intersp = 0.9)
dev.off()

same <- sect[em[,1]] == sect[em[,2]]
cat(sprintf("%d edges | %d within group | %d crossing | %d rows (no edge labels)\n",
            nrow(em), sum(same), sum(!same), length(rows)))
cat(sprintf("wrote %s (%.2f x %.2f in, label cex %.2f)\n", OUT, W_IN, H_IN, LCEX))
