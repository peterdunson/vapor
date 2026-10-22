# ============================================================================
#  Figure 1: marginal correlation among the rate series against the
#  conditional dependence chain recovered by VAPOR, p = 100.
#  Left panel needs the training data, so the FRED-QD prep is repeated here.
#  Right panel reads ci_edges, so the chain and its labels follow the fit.
# ============================================================================
RES_DIR  <- "/Users/peterdunson/Desktop/vapor/applications/results"
FIG_DIR  <- file.path(RES_DIR, "figures")
FITS_DIR <- file.path(RES_DIR, "fits")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
APP <- FITS_DIR            # fits are read from here
fit <- readRDS(file.path(APP, "p100_results.rds"))
OUT <- file.path(FIG_DIR, "termstructure_p100.png")

RATES <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTGAGE30US")
SHOW  <- RATES                                  # full names on both panels

# -- rebuild the training panel ---------------------------------------------
raw <- read.csv("/Users/peterdunson/Desktop/vapor/applications/data/2026-05-qd.csv")
tf <- function(x, tcode){ x <- as.numeric(x)
switch(as.character(tcode),"1"=x,"2"=c(NA,diff(x)),"3"=c(NA,NA,diff(diff(x))),
       "4"=log(x),"5"=c(NA,diff(log(x))),"6"=c(NA,NA,diff(diff(log(x)))), x) }
tc <- as.numeric(raw[2,-1]); names(tc) <- colnames(raw)[-1]
dat <- raw[-c(1,2),]; dts <- dat[,1]
dat <- apply(dat[,-1], 2, as.numeric)
tr  <- suppressWarnings(mapply(tf, as.data.frame(dat), tc))
tr  <- as.matrix(tr); rownames(tr) <- dts

pd <- function(s){ z <- as.numeric(strsplit(s,"/")[[1]]); z[3] + (z[1]-1)/12 }
dv <- vapply(rownames(tr), pd, numeric(1))
Y  <- tr[, fit$vars]
k  <- which(dv >= 1996 & dv <= 2025 + 11/12); Y <- Y[k,,drop=FALSE]
ok <- apply(Y, 1, function(r) !any(is.na(r))); Y <- Y[ok,,drop=FALSE]
dk <- dv[k][ok]
Ytr <- scale(Y)[dk >= 1996 & dk <= 2020 + 11/12, ]

C <- cor(Ytr[, RATES]); dimnames(C) <- list(SHOW, SHOW)
cat("marginal correlations among the rate series:\n"); print(round(C, 2))

# -- the conditional chain, read from the fit -------------------------------
ed  <- fit$ci_edges
get <- function(a, b){
   i <- which((ed$v1==a & ed$v2==b) | (ed$v1==b & ed$v2==a))
   if (!length(i)) NA else ed$pcor[i[1]] }
chain <- c(get("FEDFUNDS","TB3MS"), get("TB3MS","GS1"),
           get("GS1","GS5"),        get("GS5","GS10"))
mort  <- get("GS5","MORTGAGE30US")
cat(sprintf("\nchain: %s | mortgage link %.3f\n",
            paste(sprintf("%.3f", chain), collapse=", "), mort))
if (any(is.na(chain))) warning("an adjacent-maturity pair is not reported; the chain is broken")

nonadj <- list(c("FEDFUNDS","GS1"), c("FEDFUNDS","GS5"), c("FEDFUNDS","GS10"),
               c("TB3MS","GS5"), c("TB3MS","GS10"), c("GS1","GS10"))
n_abs <- sum(vapply(nonadj, function(z) is.na(get(z[1], z[2])), logical(1)))
cat(sprintf("non-adjacent pairs absent: %d of 6\n", n_abs))

# -- draw --------------------------------------------------------------------
# ============================================================================
#  PLOT.  Both panels come from the current pooled fit: the right panel reads
#  ci_edges and pc_vapor, so the chain, its labels and the one missing link
#  follow the fit rather than a hand-written story.  Selected edges are solid;
#  an adjacent-maturity pair that the 95% rule does NOT select is drawn dashed
#  grey with its posterior mean, so the break in the chain is visible.
#  The subtitle is gone: the caption carries it.
# ============================================================================
ed   <- fit$long$ci_edges
vars <- fit$vars
#  Labels are the POOLED POSTERIOR MEAN of the partial correlation, computed
#  from the retained draws of all four chains -- the same quantity ci_edges
#  reports.  fit$long$pc_vapor is the partial correlation OF the mean
#  precision, which is a different number (0.658 vs 0.636 for FEDFUNDS-TB3MS),
#  and using it would put values on the figure that the edge table contradicts.
NEED <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTGAGE30US")
ni <- match(NEED, vars)
pool <- do.call(rbind, lapply(fit$draws$vapor, function(f){
   nd <- dim(f$Lambda)[3]
   t(vapply(seq_len(nd), function(s){
      Om <- tcrossprod(f$Lambda[,,s]); diag(Om) <- diag(Om) + f$residual_prec[s,]
      dd <- sqrt(diag(Om)); pc <- -Om/outer(dd,dd)
      pc[ni, ni][upper.tri(diag(length(ni)))] }, numeric(choose(length(NEED),2)))) }))
pcm <- matrix(0, length(NEED), length(NEED), dimnames = list(NEED, NEED))
pcm[upper.tri(pcm)] <- colMeans(pool)
pcm <- pcm + t(pcm)
ekey <- paste(pmin(ed$v1, ed$v2), pmax(ed$v1, ed$v2), sep = "--")
selected <- function(a, b) paste(min(a,b), max(a,b), sep = "--") %in% ekey
pcor_of  <- function(a, b) pcm[a, b]

CHAIN <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
LABS  <- CHAIN

png(OUT, width = 13.5, height = 6.0, units = "in", res = 300, bg = "white")
layout(matrix(1:2, nrow = 1), widths = c(1, 1.15))

## left: marginal correlation heatmap
par(mar = c(1, 10.5, 11.0, 1))
pal <- colorRampPalette(c("#FDECE4", "#B2182B"))(101)
n <- nrow(C)
image(1:n, 1:n, t(C[n:1, ]), col = pal, zlim = c(0, 1), axes = FALSE,
      xlab = "", ylab = "", useRaster = TRUE)
for (i in 1:n) for (j in 1:n){
   v <- C[i, j]
   text(j, n - i + 1, sprintf("%.2f", v), cex = 1.15,
        col = if (v > 0.72) "white" else "grey15") }
axis(3, at = 1:n, labels = SHOW, tick = FALSE, las = 2, cex.axis = 1.02, line = -0.6)
axis(2, at = n:1, labels = SHOW, tick = FALSE, las = 1, cex.axis = 1.02, line = -0.6)
mtext("Marginal correlation", side = 3, line = 8.6, font = 2, cex = 1.69)
abline(h = 0:n + 0.5, v = 0:n + 0.5, col = "white", lwd = 2)

## right: the conditional chain, from the current graph
par(mar = c(1, 1, 11.0, 1))
plot(NA, xlim = c(0.3, 5.7), ylim = c(-1.9, 1.5), axes = FALSE, xlab = "", ylab = "")
mtext("Conditional dependence (VAPOR)", side = 3, line = 8.6, font = 2, cex = 1.69)

xs <- seq_along(CHAIN); ys <- rep(0.55, length(CHAIN))
SOLID <- "#2E5E9E"; DASH <- "grey55"; NODE_BG <- "#D9E6F5"
NODE_CEX <- 4.4; LAB_DROP <- 0.40          # labels sit below their node
n_dashed <- 0L
for (k in seq_len(length(CHAIN) - 1L)){
   a <- CHAIN[k]; b <- CHAIN[k+1]; sel <- selected(a, b)
   if (!sel) n_dashed <- n_dashed + 1L
   segments(xs[k], ys[k], xs[k+1], ys[k+1],
            lwd = if (sel) 7 else 3.4, col = if (sel) SOLID else DASH,
            lty = if (sel) 1 else 2)
   lab <- sprintf("%.2f", pcor_of(a, b))
   xm <- mean(xs[k:(k+1)]); ym <- 0.55
   rect(xm - 0.26, ym + 0.17, xm + 0.26, ym + 0.45, col = "white", border = NA)
   text(xm, ym + 0.31, lab, cex = 1.25, font = if (sel) 2 else 1,
        col = if (sel) "grey10" else DASH) }

if (!selected("GS5", "MORTGAGE30US")) n_dashed <- n_dashed + 1L
any_dashed <- n_dashed > 0L
if (selected("GS5", "MORTGAGE30US")){
   segments(4, 0.55, 4, -0.95, lwd = 7, col = SOLID)
   rect(4.04, -0.34, 4.66, -0.06, col = "white", border = NA)
   text(4.35, -0.20, sprintf("%.2f", pcor_of("GS5","MORTGAGE30US")),
        cex = 1.25, font = 2, col = "grey10")
   points(4, -0.95, pch = 21, bg = NODE_BG, col = "grey25", cex = NODE_CEX, lwd = 1.6)
   text(4, -0.95 - LAB_DROP, "MORTGAGE30US", cex = 1.02) }

points(xs, ys, pch = 21, bg = NODE_BG, col = "grey25", cex = NODE_CEX, lwd = 1.6)
#  GS5 carries both the dashed GS1-GS5 link and the branch down to
#  MORTGAGE30US, so its name goes ABOVE-LEFT of the node, high enough to clear
#  the edge values on either side; every other label sits below its node.
lab_dx <- rep(0, length(CHAIN)); lab_dy <- rep(-LAB_DROP, length(CHAIN))
lab_dx[CHAIN == "GS5"] <- -0.22
lab_dy[CHAIN == "GS5"] <-  0.62
text(xs + lab_dx, ys + lab_dy, LABS, cex = 1.02)
#  the dashed entry is shown only when at least one dashed edge was drawn
LG <- c("selected by the 95% rule", "95% interval includes zero")
LC <- c(SOLID, DASH); LW <- c(7, 3.4); LT <- c(1, 2)
if (!any_dashed){ LG <- LG[1]; LC <- LC[1]; LW <- LW[1]; LT <- LT[1] }
legend("bottomleft", inset = c(0.01, 0.00), bty = "n", cex = 1.05,
       legend = LG, col = LC, lwd = LW, lty = LT, seg.len = 2.2)
cat(sprintf("dashed edges drawn: %d; legend entries: %d\n", n_dashed, length(LG)))
dev.off()
cat(sprintf("wrote %s\n", OUT))
cat("chain links:\n")
for (k in seq_len(length(CHAIN)-1L))
   cat(sprintf("  %-12s %-12s pcor %+.3f  %s\n", CHAIN[k], CHAIN[k+1],
               pcor_of(CHAIN[k], CHAIN[k+1]),
               if (selected(CHAIN[k], CHAIN[k+1])) "selected" else "NOT selected (dashed)"))
