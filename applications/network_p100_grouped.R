# ============================================================================
#  p=100 conditional dependence graph, hand-placed layout
#  Updated for the 25-edge graph: IPNMAT and IPMANSICS drop out; TCU, CUMFNS,
#  and PCNDx enter. Everything reads ci_edges, so magnitudes and widths follow
#  the fit automatically.
# ============================================================================
library(igraph)

APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
fit <- readRDS(file.path(APP, "p100_full_results.rds"))
e   <- fit$ci_edges
OUT <- file.path(APP, "network_p100_grouped.png")

sect <- c(
   GDPC1="Real activity", GPDIC1="Real activity", INDPRO="Real activity",
   IPFINAL="Real activity", IPCONGD="Real activity", IPBUSEQ="Real activity",
   IPMAT="Real activity", AMDMNOx="Real activity", ANDENOx="Real activity",
   TCU="Real activity", CUMFNS="Real activity",
   HOUST="Housing", HOUSTW="Housing", MORTGAGE30US="Housing",
   PCECC96="Consumption", PCDGx="Consumption", PCESVx="Consumption",
   PCNDx="Consumption",
   UNRATE="Labor", UNRATESTx="Labor", UNRATELTx="Labor", UEMP27OV="Labor",
   MANEMP="Labor", DMANEMP="Labor", USGOOD="Labor", USCONS="Labor",
   CPIAUCSL="Prices", CPITRNSL="Prices", PPIACO="Prices", WPSID61="Prices",
   FEDFUNDS="Rates & credit", TB3MS="Rates & credit", GS1="Rates & credit",
   GS5="Rates & credit", GS10="Rates & credit", AAA="Rates & credit",
   BAA="Rates & credit"
)

# fail loudly if the fit contains a series the sector table does not know
missing <- setdiff(unique(c(e$v1, e$v2)), names(sect))
if (length(missing)) stop("add to sect: ", paste(missing, collapse=", "))

sector_levels <- c("Real activity","Housing","Consumption","Labor","Prices","Rates & credit")
pal <- setNames(c("#AEC7E8","#9ED6D0","#A8D5A2","#D6B3D6","#F2D98C","#F2A6A6"), sector_levels)

g <- graph_from_data_frame(e[, c("v1","v2","pcor")], directed = FALSE,
                           vertices = data.frame(name = names(sect), sector = unname(sect)))
V(g)$sector <- sect[V(g)$name]
V(g)$color  <- pal[V(g)$sector]

# -- hand-placed coordinates ------------------------------------------------
pos <- rbind(
   # row 4 (top): three two-node pairs
   INDPRO=c(1,7),    IPMAT=c(2.6,7),
   TCU=c(4.4,7),     CUMFNS=c(6,7),
   AMDMNOx=c(7.8,7), ANDENOx=c(9.4,7),
   # row 3
   GDPC1=c(1,5.4),   GPDIC1=c(2.6,5.4),
   HOUST=c(4.4,5.4), HOUSTW=c(6,5.4),
   USGOOD=c(7.8,5.4),USCONS=c(9.4,5.4),
   # row 2
   UNRATE=c(1,3.8),  UNRATESTx=c(2.6,3.8),
   UNRATELTx=c(4.4,3.8), UEMP27OV=c(6,3.8),
   MANEMP=c(7.8,3.8),DMANEMP=c(9.4,3.8),
   # row 1: the IP triangle, the consumption group, price and credit dyads
   IPCONGD=c(1.6,2.1), IPBUSEQ=c(0.8,1.1), IPFINAL=c(2.4,1.1),
   PCDGx=c(4.4,2.1),   PCESVx=c(3.6,1.1),  PCECC96=c(5.2,1.1),
   PCNDx=c(4.5,0.15),
   CPIAUCSL=c(7.3,1.6),CPITRNSL=c(8.7,1.6),
   PPIACO=c(9.4,2.1),  WPSID61=c(9.4,1.0),
   AAA=c(6.2,2.4),     BAA=c(6.2,0.7),
   # bottom row: the interest-rate chain, spread wide
   GS10=c(1.2,-1.1), GS5=c(2.9,-1.1), GS1=c(4.6,-1.1),
   TB3MS=c(6.3,-1.1), FEDFUNDS=c(8.0,-1.1), MORTGAGE30US=c(2.9,0.3)
)
lay <- pos[V(g)$name, ]

V(g)$size         <- 13
V(g)$label.cex    <- 0.62
V(g)$label.color  <- "black"
V(g)$label.family <- "sans"
V(g)$frame.color  <- "grey25"

ends_mat  <- ends(g, E(g))
same_sect <- V(g)$sector[match(ends_mat[,1], V(g)$name)] ==
   V(g)$sector[match(ends_mat[,2], V(g)$name)]
E(g)$crossing <- !same_sect
E(g)$width  <- 1.5 + 5 * abs(E(g)$pcor)
E(g)$lty    <- ifelse(E(g)$pcor < 0, 2, 1)
E(g)$color  <- ifelse(E(g)$crossing, "grey10", "grey60")

cat(sprintf("%d edges | %d within sector | %d crossing\n",
            ecount(g), sum(same_sect), sum(!same_sect)))
cat("crossing: ", paste(apply(ends_mat[!same_sect, , drop=FALSE], 1,
                              paste, collapse="--"), collapse="; "), "\n")

png(OUT, width = 12, height = 9, units = "in", res = 300, bg = "white")
par(mar = c(2.5, 0.5, 0.5, 0.5))
plot(g, layout = lay, rescale = TRUE, asp = 0, edge.label = NA)

xy <- layout.norm(lay, -1, 1, -1, 1); rownames(xy) <- V(g)$name
el <- ends(g, E(g))
for (i in seq_len(ecount(g))) {
   mx <- mean(xy[el[i,],1]); my <- mean(xy[el[i,],2])
   lab <- sprintf("%.2f", E(g)$pcor[i])
   w <- strwidth(lab, cex=0.6)*1.35; h <- strheight(lab, cex=0.6)*1.7
   rect(mx-w/2, my-h/2, mx+w/2, my+h/2, col="white", border=NA)
   text(mx, my, lab, cex=0.6, col="grey20")
}
legend("bottom", legend=sector_levels, pt.bg=pal[sector_levels], pch=21,
       pt.cex=1.7, cex=0.8, bty="n", horiz=TRUE, inset=-0.01, x.intersp=0.6)
dev.off()
cat("wrote", OUT, "\n")