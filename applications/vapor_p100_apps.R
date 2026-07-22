# ============================================================================
#  VAPOR p=100: edge-overlap comparison, network figure, motivation figure
#  Requires p100_full_results.rds. Section C additionally needs Y_train.
# ============================================================================
library(igraph)

OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application"
r   <- readRDS(file.path(path.expand(OUT), "p100_full_results.rds"))
ce  <- r$ci_edges
vars_qd <- r$vars
p   <- r$p

pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }

PC <- list(
   VAPOR   = pcor_from_omega(r$raw_omega$vapor),
   Kastner = pcor_from_omega(r$raw_omega$kast),
   Glasso  = pcor_from_omega(r$raw_omega$glasso),
   Sample  = pcor_from_omega(r$raw_omega$samp)
)

# ============================================================================
#  SECTION A: do competitors find VAPOR's 23 edges?
# ============================================================================

ut  <- which(upper.tri(matrix(0,p,p)), arr.ind = TRUE)          # all pairs
key <- function(i,j) paste(pmin(i,j), pmax(i,j), sep = "-")

# VAPOR's reported edge set
ei <- cbind(match(ce$v1, vars_qd), match(ce$v2, vars_qd))
vapor_keys <- key(ei[,1], ei[,2])
nE <- nrow(ei)

# rank every pair by |partial correlation| under each estimator
rank_tab <- lapply(PC, function(M){
   v <- abs(M[cbind(ut[,1], ut[,2])])
   o <- order(v, decreasing = TRUE)
   rk <- integer(length(v)); rk[o] <- seq_along(v)
   setNames(rk, key(ut[,1], ut[,2]))
})

# where does each VAPOR edge rank under each competitor?
ranks <- sapply(rank_tab, function(rk) rk[vapor_keys])
rownames(ranks) <- paste(ce$v1, ce$v2, sep = "-")

cat("=== rank of each VAPOR-reported edge under each estimator ===\n")
cat("    (1 = strongest partial correlation of all 4,950 pairs)\n\n")
print(cbind(pcor = ce$pcor, as.data.frame(ranks)))

# overlap at several cutoffs
cutoffs <- c(nE, 50, 100, 250)
ov <- sapply(cutoffs, function(K)
   sapply(rank_tab, function(rk) sum(rk[vapor_keys] <= K)))
dimnames(ov) <- list(names(PC), paste0("top", cutoffs))
cat(sprintf("\n=== how many of VAPOR's %d edges fall in each estimator's top-K ===\n", nE))
print(ov)

cat("\n=== median rank of VAPOR's edges under each estimator ===\n")
print(round(apply(ranks, 2, median), 1))

# and the reverse: what does each competitor rank top-23 that VAPOR does not report?
cat("\n=== each estimator's top-5 pairs (for contrast) ===\n")
for (nm in names(PC)) {
   rk <- rank_tab[[nm]]
   top <- names(sort(rk))[1:5]
   ij  <- do.call(rbind, lapply(strsplit(top, "-"), as.integer))
   cat(sprintf("%-8s %s\n", nm,
               paste(sprintf("%s-%s", vars_qd[ij[,1]], vars_qd[ij[,2]]), collapse = ", ")))
}

# ============================================================================
#  SECTION B: network figure of the 23 reported edges
# ============================================================================

sector_map <- c(
   GDPC1="Output", GPDIC1="Output",
   INDPRO="Industrial", IPMAT="Industrial", IPFINAL="Industrial",
   IPCONGD="Industrial", IPBUSEQ="Industrial", IPNMAT="Industrial",
   IPMANSICS="Industrial",
   AMDMNOx="Orders", ANDENOx="Orders",
   UNRATE="Labor", UNRATESTx="Labor", UNRATELTx="Labor", UEMP27OV="Labor",
   MANEMP="Labor", DMANEMP="Labor", USGOOD="Labor", USCONS="Labor",
   PCECC96="Consumption", PCESVx="Consumption", PCDGx="Consumption",
   CPIAUCSL="Prices", CPITRNSL="Prices", PPIACO="Prices", WPSID61="Prices",
   FEDFUNDS="Rates", TB3MS="Rates", GS1="Rates", GS5="Rates", GS10="Rates",
   AAA="Credit", BAA="Credit", MORTGAGE30US="Credit",
   HOUST="Housing", HOUSTW="Housing"
)

el <- data.frame(from = ce$v1, to = ce$v2, pcor = ce$pcor, stringsAsFactors = FALSE)
g  <- graph_from_data_frame(el, directed = FALSE)

sec <- sector_map[V(g)$name]; sec[is.na(sec)] <- "Other"
sectors <- sort(unique(sec))
pal <- setNames(
   c("#4E79A7","#59A14F","#EDC948","#E15759","#B07AA1","#76B7B2","#FF9DA7","#9C755F")[seq_along(sectors)],
   sectors)

V(g)$color <- pal[sec]
V(g)$size  <- 14
V(g)$label.cex <- 0.62
V(g)$label.color <- "black"
V(g)$frame.color <- "grey30"
E(g)$width <- 1 + 5 * abs(E(g)$pcor)
E(g)$color <- ifelse(E(g)$pcor > 0, "#2C6FBB", "#C0392B")
E(g)$lty   <- ifelse(E(g)$pcor > 0, 1, 2)

set.seed(11)
lay <- layout_with_fr(g, weights = abs(E(g)$pcor) * 3)

png(file.path(path.expand(OUT), "network_p100.png"), width = 2000, height = 1500, res = 200)
par(mar = c(0,0,1,0))
plot(g, layout = lay, vertex.label.family = "sans")
legend("bottomleft", legend = sectors, pt.bg = pal[sectors], pch = 21,
       pt.cex = 1.6, bty = "n", cex = 0.85, title = "Sector")
legend("bottomright",
       legend = c("positive partial correlation", "negative partial correlation"),
       col = c("#2C6FBB", "#C0392B"), lty = c(1,2), lwd = 2, bty = "n", cex = 0.85)
dev.off()
cat("\nSaved network_p100.png\n")

# structural facts to quote in the text
cross <- el[sector_map[el$from] != sector_map[el$to], ]
cat(sprintf("\nedges: %d | within-sector: %d | cross-sector: %d\n",
            nrow(el), nrow(el) - nrow(cross), nrow(cross)))
cat("cross-sector edges:\n"); print(cross, row.names = FALSE)
cat(sprintf("negative edges: %d\n", sum(el$pcor < 0)))
print(el[el$pcor < 0, ], row.names = FALSE)

# the maturity ladder: which rate pairs are connected vs absent
rates <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
present <- apply(el[,1:2], 1, function(z) paste(sort(z), collapse="|"))
cat("\n=== term-structure adjacency ===\n")
for (i in 1:(length(rates)-1)) for (j in (i+1):length(rates)) {
   k <- paste(sort(c(rates[i], rates[j])), collapse="|")
   cat(sprintf("  %-10s %-10s gap %d  %s\n", rates[i], rates[j], j-i,
               ifelse(k %in% present, "EDGE", "-")))
}

# ============================================================================
#  SECTION C: three-panel motivation figure
#  Panel 1: lag-1 autocorrelation across series
#  Panel 2: density of |innovation correlations| after removing dynamics
#  Panel 3: eigenvalue decay of the sample covariance
# ============================================================================

if (!exists("Y_train")) {
   message("Y_train not in session — skipping Section C.")
} else {
   
   Ytr <- as.matrix(Y_train); T_tr <- nrow(Ytr); p_tr <- ncol(Ytr)
   
   # panel 1
   ac1 <- apply(Ytr, 2, function(x) acf(x, lag.max = 1, plot = FALSE)$acf[2])
   
   # panel 2: residuals from a ridge VAR(1)
   X <- Ytr[-T_tr, , drop = FALSE]; Ym <- Ytr[-1, , drop = FALSE]
   Bhat <- t(solve(crossprod(X) + 1e-2 * nrow(X) * diag(p_tr), crossprod(X, Ym)))
   E <- Ym - X %*% t(Bhat)
   Craw <- cor(Ytr);  Cinn <- cor(E)
   raw_off <- abs(Craw[upper.tri(Craw)]); inn_off <- abs(Cinn[upper.tri(Cinn)])
   
   # panel 3
   ev <- eigen(cov(Ytr), symmetric = TRUE, only.values = TRUE)$values
   ev <- pmax(ev, .Machine$double.eps)
   
   png(file.path(path.expand(OUT), "motivation_p100.png"), width = 1800, height = 560, res = 170)
   par(mfrow = c(1,3), mar = c(4.2,4.2,3,1), mgp = c(2.5,0.7,0), bty = "l", las = 1)
   
   hist(ac1, breaks = 20, col = "grey85", border = "white",
        main = "Series are autocorrelated", xlab = "Lag-1 autocorrelation", ylab = "Number of series")
   abline(v = median(ac1), col = "#C0392B", lwd = 2)
   
   plot(density(raw_off), lwd = 2, col = "grey45",
        xlim = c(0,1), main = "Innovations stay dependent",
        xlab = "|correlation|", ylab = "Density")
   lines(density(inn_off), lwd = 2, col = "#2C6FBB")
   legend("topright", c("raw series", "VAR innovations"),
          col = c("grey45", "#2C6FBB"), lwd = 2, bty = "n", cex = 0.9)
   
   plot(ev, log = "y", type = "b", pch = 16, cex = 0.5, col = "#2C6FBB",
        main = sprintf("Covariance is ill-conditioned (T/p = %.2f)", T_tr/p_tr),
        xlab = "Eigenvalue index", ylab = "Eigenvalue (log scale)")
   abline(h = ev[1]/1e3, col = "grey70", lty = 3)
   
   dev.off(); cat("\nSaved motivation_p100.png\n")
   
   cat(sprintf("\nmedian lag-1 autocorrelation: %.2f | share > 0.5: %.0f%%\n",
               median(ac1), 100*mean(ac1 > 0.5)))
   cat(sprintf("mean |corr|: raw %.3f -> innovations %.3f\n",
               mean(raw_off), mean(inn_off)))
   cat(sprintf("condition number of sample covariance: %.3g\n", ev[1]/ev[length(ev)]))
}












# ============================================================================
#  Term-structure figure: marginal correlation vs conditional graph
# ============================================================================
OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application"
stopifnot(exists("Y_train"))
r  <- readRDS(file.path(path.expand(OUT), "p100_full_results.rds"))
ce <- r$ci_edges; vars_qd <- r$vars

rates <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTGAGE30US")
lab   <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTG30")
stopifnot(all(rates %in% vars_qd))

Cm <- cor(Y_train[, rates])

pc_of <- function(a, b){
   i <- which((ce$v1==a & ce$v2==b) | (ce$v1==b & ce$v2==a))
   if (length(i)) ce$pcor[i] else NA_real_
}
edges <- data.frame(
   from = c("FEDFUNDS","TB3MS","GS1","GS5","GS5"),
   to   = c("TB3MS","GS1","GS5","GS10","MORTGAGE30US"),
   stringsAsFactors = FALSE)
edges$pcor <- mapply(pc_of, edges$from, edges$to)
stopifnot(!any(is.na(edges$pcor)))

png(file.path(path.expand(OUT), "termstructure_p100.png"),
    width = 1900, height = 860, res = 190)
layout(matrix(1:2, nrow = 1), widths = c(1, 1.15))

# ---- Panel A: marginal correlations ----------------------------------------
par(mar = c(1.5, 6.2, 6.2, 1.5))          # generous top margin for labels + title
n <- length(rates)
pal <- colorRampPalette(c("white", "#F4A582", "#B2182B"))(101)
plot(NA, xlim = c(0.5, n+0.5), ylim = c(0.5, n+0.5), axes = FALSE,
     xlab = "", ylab = "", asp = 1, main = "")
for (i in 1:n) for (j in 1:n) {
   v <- Cm[i, j]
   rect(j-0.5, n-i+0.5, j+0.5, n-i+1.5,
        col = pal[1 + round(100*max(0, min(1, abs(v))))], border = "white")
   text(j, n-i+1, sprintf("%.2f", v), cex = 0.62,
        col = if (abs(v) > 0.7) "white" else "grey20")
}
axis(2, at = n:1, labels = lab, tick = FALSE, las = 1, cex.axis = 0.72, line = -0.6)
axis(3, at = 1:n, labels = lab, tick = FALSE, las = 2, cex.axis = 0.72, line = -0.9)
mtext("Marginal correlation", side = 3, line = 4.6, font = 2, cex = 1.02)

# ---- Panel B: conditional graph --------------------------------------------
par(mar = c(1.5, 1, 6.2, 1))
pos <- rbind(
   FEDFUNDS     = c(1.0, 1.0), TB3MS = c(2.0, 1.0), GS1 = c(3.0, 1.0),
   GS5          = c(4.0, 1.0), GS10  = c(5.0, 1.0),
   MORTGAGE30US = c(4.0, 0.25))
plot(NA, xlim = c(0.45, 5.55), ylim = c(0.05, 1.30), axes = FALSE,
     xlab = "", ylab = "", main = "")
for (k in seq_len(nrow(edges))) {
   a <- pos[edges$from[k], ]; b <- pos[edges$to[k], ]
   segments(a[1], a[2], b[1], b[2], lwd = 1 + 9*abs(edges$pcor[k]), col = "#2C6FBB")
   if (edges$to[k] == "MORTGAGE30US")
      text(b[1] + 0.07, (a[2]+b[2])/2, sprintf("%.2f", edges$pcor[k]),
           cex = 0.68, col = "#1B4F72", pos = 4)
   else
      text((a[1]+b[1])/2, a[2] + 0.11, sprintf("%.2f", edges$pcor[k]),
           cex = 0.68, col = "#1B4F72")
}
cols <- ifelse(rownames(pos) == "MORTGAGE30US", "#59A14F", "#EDC948")
points(pos[,1], pos[,2], pch = 21, bg = cols, col = "grey25", cex = 4.4, lwd = 1.2)
text(pos[,1], pos[,2] - 0.155, lab, cex = 0.68)
mtext("Conditional dependence (VAPOR)", side = 3, line = 4.6, font = 2, cex = 1.02)
mtext("adjacent maturities only; all six non-adjacent pairs absent",
      side = 3, line = 3.1, cex = 0.72, col = "grey35")
dev.off()
cat("Saved termstructure_p100.png\n")

off <- Cm[1:5, 1:5][upper.tri(Cm[1:5, 1:5])]
cat(sprintf("marginal among five yields: range %.2f to %.2f, median %.2f\n",
            min(off), max(off), median(off)))