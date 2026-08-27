# ============================================================================
#  Figure 1: marginal correlation among the rate series against the
#  conditional dependence chain recovered by VAPOR, p = 100.
#  Left panel needs the training data, so the FRED-QD prep is repeated here.
#  Right panel reads ci_edges, so the chain and its labels follow the fit.
# ============================================================================
APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
fit <- readRDS(file.path(APP, "p100_full_results.rds"))
OUT <- file.path(APP, "termstructure_p100.png")

RATES <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTGAGE30US")
SHOW  <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTG30")   # display labels

# -- rebuild the training panel ---------------------------------------------
raw <- read.csv("~/Desktop/2026-05-qd.csv")
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
png(OUT, width = 13.5, height = 6.0, units = "in", res = 300, bg = "white")
layout(matrix(1:2, nrow = 1), widths = c(1, 1.15))

## left: correlation heatmap
par(mar = c(1, 7.5, 8.5, 1))
pal <- colorRampPalette(c("#FDECE4", "#B2182B"))(101)
n <- nrow(C)
image(1:n, 1:n, t(C[n:1, ]), col = pal, zlim = c(0, 1), axes = FALSE,
      xlab = "", ylab = "", useRaster = TRUE)
for (i in 1:n) for (j in 1:n){
   v <- C[i, j]
   text(j, n - i + 1, sprintf("%.2f", v), cex = 1.15,
        col = if (v > 0.72) "white" else "grey15") }
axis(3, at = 1:n, labels = SHOW, tick = FALSE, las = 2, cex.axis = 1.15, line = -0.6)
axis(2, at = n:1, labels = SHOW, tick = FALSE, las = 1, cex.axis = 1.15, line = -0.6)
mtext("Marginal correlation", side = 3, line = 6.2, font = 2, cex = 1.35)
abline(h = 0:n + 0.5, v = 0:n + 0.5, col = "white", lwd = 2)

## right: conditional chain
par(mar = c(1, 1, 8.5, 1))
plot(NA, xlim = c(0.3, 5.7), ylim = c(-1.5, 1.4), axes = FALSE, xlab = "", ylab = "")
mtext("Conditional dependence (VAPOR)", side = 3, line = 6.2, font = 2, cex = 1.35)
mtext(sprintf("adjacent maturities only; all %s non-adjacent pairs absent",
              c("no","one","two","three","four","five","six")[n_abs + 1]),
      side = 3, line = 4.2, cex = 1.0, col = "grey35")

xs <- 1:5; ys <- rep(0.55, 5)
nm <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
segments(xs[-5], ys[-5], xs[-1], ys[-1], lwd = 7, col = "#2E5E9E")
segments(4, 0.55, 4, -0.95, lwd = 7, col = "#2E5E9E")
points(4, -0.95, pch = 21, bg = "#2E8B3E", col = "grey20", cex = 5.5, lwd = 1.5)
points(xs, ys, pch = 21, bg = "#F2C230", col = "grey20", cex = 5.5, lwd = 1.5)
text(xs, ys - 0.33, nm, cex = 1.1)
text(4, -1.28, "MORTG30", cex = 1.1)
text((xs[-5] + xs[-1]) / 2, 0.78, sprintf("%.2f", chain), cex = 1.1, col = "#1F3F6E")
text(4.3, -0.2, sprintf("%.2f", mort), cex = 1.1, col = "#1F3F6E", adj = 0)
dev.off()
cat("wrote", OUT, "\n")