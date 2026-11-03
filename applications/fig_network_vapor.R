# ============================================================================
#  The same 26-edge p = 100 dependence graph as network_p100_grouped.png, with
#  the node positions chosen so the drawing spells VAPOR.
#     Rscript applications/fig_network_vapor.R
#  Writes results/figures/network_p100_vapor.png.
#
#  Nothing about the graph changes: the same 40 nodes, the same 26 edges, the
#  same edge colours (blue positive, red negative) and the same FRED-QD group
#  colours and legend as the grouped figure.  Only the coordinates differ, and
#  every edge is still a straight segment between its own two endpoints.
#
#  The graph has SIXTEEN connected components, not thirteen: two triangles
#  (IP, PCE), two bent paths of three (housing, prices), one six-node tree
#  (the interest rates) and eleven pairs.  All sixteen are used:
#
#    V  housing path, bent at HOUST
#    A  the interest-rate tree: GS1 is the apex, GS1-TB3MS-FEDFUNDS the left
#       leg, GS1-GS5-GS10 the right leg, GS5-MORTGAGE30US the crossbar, with
#       MORTGAGE30US placed on the left leg so the bar spans the letter
#    P  PCE triangle as the bowl; AAA-BAA and M1REAL-M2REAL stacked as the stem
#    O  an octagon: the prices path supplies two adjacent sides meeting at
#       CPITRNSL, six pairs supply the other six
#    R  IP triangle as the bowl; TCU-CUMFNS and INDPRO-IPMAT as the stem;
#       AMDMNOx-ANDENOx as the diagonal leg
# ============================================================================
suppressMessages(library(igraph))
RES_DIR  <- "/Users/peterdunson/Desktop/vapor/applications/results"
FIG_DIR  <- file.path(RES_DIR, "figures")
FITS_DIR <- file.path(RES_DIR, "fits")
GROUPS   <- "/Users/peterdunson/Desktop/vapor/applications/fredqd_groups.csv"
OUT      <- file.path(FIG_DIR, "network_p100_vapor.png")

fit <- readRDS(file.path(FITS_DIR, "p100_results.rds"))
e   <- fit$long$ci_edges
g   <- read.csv(GROUPS, stringsAsFactors = FALSE)
gname <- setNames(g$group_name, g$series); gcode <- setNames(g$group_code, g$series)
NODES <- sort(unique(c(e$v1, e$v2)))
sect  <- gname[NODES]; scode <- gcode[NODES]
cat(sprintf("%d nodes, %d edges\n", length(NODES), nrow(e)))

#  the same palette and legend order as network_p100_grouped.R
lv <- unique(data.frame(code = unname(scode), name = unname(sect),
                        stringsAsFactors = FALSE))
lv <- lv[order(lv$code), ]
PALETTE <- c("#8FB8DE","#9FD3A8","#F2C14E","#E8927C","#C3A6D8",
             "#7FC9C4","#F2A6C3","#BFD56B","#D9B38C","#A8A8D8")
pal <- setNames(PALETTE[seq_len(nrow(lv))], lv$name)

# ---------------------------------------------------------------------------
#  letter geometry.  Each letter is laid out in its own box, x in [0,1] and
#  y in [0,2], then shifted right by its offset.
# ---------------------------------------------------------------------------
P <- matrix(NA_real_, length(NODES), 2, dimnames = list(NODES, c("x","y")))
put <- function(nm, x, y) P[nm, ] <<- c(x, y)

## ---- V: the housing path, bent at HOUST ----------------------------------
put("HOUSTW", 0.00, 2.00); put("HOUST", 0.50, 0.00); put("HOUSTS", 1.00, 2.00)

## ---- A: the interest-rate tree -------------------------------------------
#  both legs are straight: TB3MS sits on GS1-FEDFUNDS, GS5 on GS1-GS10, and
#  MORTGAGE30US sits on the left leg at GS5's height so GS5-MORTGAGE30US is a
#  horizontal crossbar between the two legs
AX <- 1.18
put("GS1",          AX + 0.500, 2.00)
put("TB3MS",        AX + 0.375, 1.50)
put("FEDFUNDS",     AX + 0.000, 0.00)
put("GS5",          AX + 0.775, 0.90)
put("GS10",         AX + 1.000, 0.00)
put("MORTGAGE30US", AX + 0.225, 0.90)

## ---- P: PCE triangle as the bowl, two pairs as the stem ------------------
PX <- 2.42
put("AAA",     PX + 0.15, 2.00); put("BAA",    PX + 0.15, 1.06)
put("M1REAL",  PX + 0.15, 0.94); put("M2REAL", PX + 0.15, 0.00)
put("PCDGx",   PX + 0.22, 2.00)          # bowl, top-left
put("PCECC96", PX + 0.95, 1.56)          # bowl, right
put("PCESVx",  PX + 0.22, 1.12)          # bowl, bottom-left

## ---- O: an octagon -------------------------------------------------------
OX <- 3.56; ocx <- OX + 0.50; ocy <- 1.00; orx <- 0.42; ory <- 0.95
VT <- t(vapply(0:7, function(k){ th <- pi/2 - k*2*pi/8
   c(ocx + orx*cos(th), ocy + ory*sin(th)) }, numeric(2)))
shrink <- function(a, b, f) list(a = a + f*(b - a), b = b - f*(b - a))
SH <- 0.11
side <- function(k) list(a = VT[k + 1, ], b = VT[(k + 1) %% 8 + 1, ])
#  the prices path supplies two adjacent sides, meeting exactly at CPITRNSL
s0 <- side(0); s1 <- side(1)
put("PCECTPI",  (s0$a + SH*(s0$b - s0$a))[1], (s0$a + SH*(s0$b - s0$a))[2])
put("CPITRNSL", VT[2,1], VT[2,2])
put("CPIAUCSL", (s1$b - SH*(s1$b - s1$a))[1], (s1$b - SH*(s1$b - s1$a))[2])
#  six pairs on the remaining six sides
RING <- list(c("UNRATE","UNRATESTx"), c("MANEMP","DMANEMP"),
             c("USGOOD","USCONS"),    c("UNRATELTx","UEMP27OV"),
             c("USPRIV","USLAH"),     c("PPIACO","WPSID61"))
for (i in seq_along(RING)){
   s <- side(i + 1); q <- shrink(s$a, s$b, SH)
   put(RING[[i]][1], q$a[1], q$a[2]); put(RING[[i]][2], q$b[1], q$b[2]) }

## ---- R: IP triangle as the bowl, two pairs as the stem, one as the leg ---
RX <- 4.72
put("TCU",     RX + 0.15, 2.00); put("CUMFNS", RX + 0.15, 1.06)
put("INDPRO",  RX + 0.15, 0.94); put("IPMAT",  RX + 0.15, 0.00)
put("IPBUSEQ", RX + 0.22, 2.00)          # bowl, top-left
put("IPFINAL", RX + 0.95, 1.56)          # bowl, right
put("IPCONGD", RX + 0.22, 1.12)          # bowl, bottom-left
put("AMDMNOx", RX + 0.30, 1.06)          # the leg
put("ANDENOx", RX + 0.98, 0.00)

stopifnot(!anyNA(P))
cat(sprintf("positions set for all %d nodes\n", nrow(P)))

# ---------------------------------------------------------------------------
#  draw
# ---------------------------------------------------------------------------
W_IN <- 6.5; H_IN <- 2.5; DPI <- 600
H_LEG <- 0.46                                  # legend strip at the foot
R_NODE <- 0.030                                # node radius, in letter units
png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI, bg = "white",
    pointsize = 10)
par(fig = c(0, 1, H_LEG/H_IN, 1), mar = c(0.02, 0.05, 0.10, 0.05))
xr <- range(P[,1]); yr <- range(P[,2])
plot(NA, xlim = c(xr[1] - 0.10, xr[2] + 0.10), ylim = c(yr[1] - 0.14, yr[2] + 0.14),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
for (i in seq_len(nrow(e))){
   a <- P[e$v1[i], ]; b <- P[e$v2[i], ]
   segments(a[1], a[2], b[1], b[2], lwd = 1.0 + 3.0*abs(e$pcor[i]),
            col = if (e$pcor[i] > 0) "#2E5E9E" else "#B2182B") }
th <- seq(0, 2*pi, length.out = 48)
for (v in NODES)
   polygon(P[v,1] + R_NODE*cos(th), P[v,2] + R_NODE*sin(th),
           col = unname(pal[sect[v]]), border = "grey30", lwd = 0.55)

par(fig = c(0, 1, 0, H_LEG/H_IN), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
legend("top", ncol = 4, bty = "n", cex = 0.52, pt.cex = 1.0,
       legend = lv$name, pch = 21, pt.bg = unname(pal[lv$name]), col = "grey30")
legend("bottom", ncol = 2, bty = "n", cex = 0.52, seg.len = 1.8, lwd = 2,
       legend = c("positive partial correlation", "negative partial correlation"),
       col = c("#2E5E9E", "#B2182B"), lty = 1)
dev.off()
cat(sprintf("wrote %s (%.2f x %.2f in, %d dpi)\n", OUT, W_IN, H_IN, DPI))
