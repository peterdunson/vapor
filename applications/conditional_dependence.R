# ============================================================================
#  p=100 conditional dependence graph, large and readable
# ============================================================================
library(igraph)

fit <- readRDS("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application/p100_full_results.rds")
e   <- fit$ci_edges

OUT <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application/network_p100_grouped.png")

sect <- c(
   GDPC1="Real activity", GPDIC1="Real activity", INDPRO="Real activity",
   IPFINAL="Real activity", IPCONGD="Real activity", IPBUSEQ="Real activity",
   IPMAT="Real activity", IPNMAT="Real activity", IPMANSICS="Real activity",
   AMDMNOx="Real activity", ANDENOx="Real activity",
   HOUST="Housing", HOUSTW="Housing", MORTGAGE30US="Housing",
   PCECC96="Consumption", PCDGx="Consumption", PCESVx="Consumption",
   UNRATE="Labor", UNRATESTx="Labor", UNRATELTx="Labor", UEMP27OV="Labor",
   MANEMP="Labor", DMANEMP="Labor", USGOOD="Labor", USCONS="Labor",
   CPIAUCSL="Prices", CPITRNSL="Prices", PPIACO="Prices", WPSID61="Prices",
   FEDFUNDS="Rates & credit", TB3MS="Rates & credit", GS1="Rates & credit",
   GS5="Rates & credit", GS10="Rates & credit", AAA="Rates & credit", BAA="Rates & credit"
)

sector_levels <- c("Real activity","Housing","Consumption","Labor","Prices","Rates & credit")
pal <- setNames(c("#AEC7E8","#9ED6D0","#A8D5A2","#D6B3D6","#F2D98C","#F2A6A6"), sector_levels)

g <- graph_from_data_frame(e[, c("v1","v2","pcor")], directed = FALSE,
                           vertices = data.frame(name = names(sect), sector = unname(sect)))
V(g)$sector <- sect[V(g)$name]
V(g)$color  <- pal[V(g)$sector]

# wider horizontal spacing (x * 1.25) and taller vertical gaps
pos <- rbind(
   INDPRO=c(1.0,9.0),   IPMAT=c(3.2,9.0),
   IPNMAT=c(5.6,9.0),   IPMANSICS=c(7.8,9.0),
   AMDMNOx=c(10.2,9.0), ANDENOx=c(12.4,9.0),
   GDPC1=c(1.0,7.0),    GPDIC1=c(3.2,7.0),
   HOUST=c(5.6,7.0),    HOUSTW=c(7.8,7.0),
   USGOOD=c(10.2,7.0),  USCONS=c(12.4,7.0),
   UNRATE=c(1.0,5.0),   UNRATESTx=c(3.2,5.0),
   UNRATELTx=c(5.6,5.0),UEMP27OV=c(7.8,5.0),
   MANEMP=c(10.2,5.0),  DMANEMP=c(12.4,5.0),
   IPCONGD=c(2.0,3.1),  IPBUSEQ=c(1.0,1.7), IPFINAL=c(3.0,1.7),
   PCDGx=c(5.6,2.9),    PCESVx=c(4.6,1.5),  PCECC96=c(6.6,1.5),
   CPIAUCSL=c(8.2,2.2), CPITRNSL=c(10.5,2.2),
   PPIACO=c(13.0,3.0),  WPSID61=c(13.0,1.3),
   AAA=c(11.5,3.6),     BAA=c(11.5,0.8),
   GS10=c(1.4,-1.2),  GS5=c(3.6,-1.2), GS1=c(5.9,-1.2),
   TB3MS=c(8.2,-1.2), FEDFUNDS=c(10.5,-1.2), MORTGAGE30US=c(3.6,0.4)
)
lay <- pos[V(g)$name, ]

V(g)$size        <- 20
V(g)$label.cex   <- 1.05
V(g)$label.color <- "black"
V(g)$label.family <- "sans"
V(g)$frame.color <- "grey25"

ends_mat  <- ends(g, E(g))
same_sect <- V(g)$sector[match(ends_mat[,1], V(g)$name)] ==
   V(g)$sector[match(ends_mat[,2], V(g)$name)]
E(g)$crossing <- !same_sect
E(g)$width  <- 2.5 + 7 * abs(E(g)$pcor)
E(g)$lty    <- ifelse(E(g)$pcor < 0, 2, 1)
E(g)$color  <- ifelse(E(g)$crossing, "grey10", "grey55")

png(OUT, width = 16, height = 13, units = "in", res = 300, bg = "white")
par(mar = c(2, 1, 1, 1))
plot(g, layout = lay, rescale = TRUE, asp = 0, edge.label = NA)

# edge numbers, offset perpendicular to each edge, no white box
xy <- layout.norm(lay, -1, 1, -1, 1); rownames(xy) <- V(g)$name
el <- ends(g, E(g))
for (i in seq_len(ecount(g))) {
   x1 <- xy[el[i,1],1]; y1 <- xy[el[i,1],2]
   x2 <- xy[el[i,2],1]; y2 <- xy[el[i,2],2]
   mx <- (x1+x2)/2; my <- (y1+y2)/2
   dx <- x2-x1; dy <- y2-y1; L <- sqrt(dx^2+dy^2)
   ox <- -dy/L * 0.035; oy <- dx/L * 0.035     # small perpendicular offset
   text(mx+ox, my+oy, sprintf("%.2f", E(g)$pcor[i]),
        cex=1.15, col="grey15", font=2)
}
legend(x = mean(par("usr")[1:2]), y = par("usr")[3] + 0.02*diff(par("usr")[3:4]),
       xjust = 0.5, yjust = 1,
       legend = sector_levels, pt.bg = pal[sector_levels], pch = 21,
       pt.cex = 3, cex = 1.6, bty = "n", horiz = TRUE, xpd = NA, x.intersp = 0.7)
dev.off()
cat("wrote", OUT, "\n")