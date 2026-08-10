# ============================================================================
#  p=100 partial-correlation heatmaps.
#  Reads p100_full_results.rds; refits nothing.
#  Boxes mark VAPOR's reported edges (95% credible interval excludes zero).
#  Box size scales with the dimension, so a single cell stays visible.
#  Labels are off by default: 23 edges over 100 variables would be
#  unreadable, and the network figure names them all.
# ============================================================================
APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
fit <- readRDS(file.path(APP, "p100_full_results.rds"))
OUT <- file.path(APP, "heatmap_p100_r3k25.png")

LABEL_TOP <- 0       # label the N strongest edges; 0 for none
LAB_CEX   <- 0.75

vars <- fit$vars; n <- length(vars)
ed   <- fit$ci_edges
ei   <- cbind(match(ed$v1, vars), match(ed$v2, vars))
up   <- t(apply(ei, 1, function(z) c(min(z), max(z))))   # upper triangle only
BOX_HW <- 0.5                             # ~1.67 at p=100
cat(sprintf("p=%d, %d reported edges, box half-width %.2f cells\n",
            n, nrow(ei), BOX_HW))

to_pc <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- NA; pc }
panels <- list("VAPOR"             = { M <- fit$pc_vapor; diag(M) <- NA; M },
               "Kastner"           = to_pc(fit$raw_omega$kast),
               "VAR + Glasso"      = to_pc(fit$raw_omega$glasso),
               "Sample covariance" = to_pc(fit$raw_omega$samp))
pal <- colorRampPalette(c("blue","white","red"))(101)

#  cell (i,j) maps to x = j, y = n - i + 1
draw_panel <- function(M, title, annotate = FALSE){
   M[M >  1] <-  1; M[M < -1] <- -1
   image(x = 1:n, y = 1:n, z = t(M[n:1, ]), col = pal, zlim = c(-1,1),
         axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
   title(main = title, family = "serif", font.main = 1, cex.main = 2.0, line = 0.6)
   if (!annotate) return(invisible())
   for (e in seq_len(nrow(up))){
      i <- up[e,1]; j <- up[e,2]; cx <- j; cy <- n - i + 1
      rect(cx-BOX_HW, cy-BOX_HW, cx+BOX_HW, cy+BOX_HW, border="black", lwd=1.2)
      if (e <= LABEL_TOP){
         text(cx, cy + BOX_HW + 0.9, vars[j], adj=c(0.5,0),
              cex=LAB_CEX, family="serif")
         text(cx - BOX_HW - 0.9, cy, vars[i], adj=c(1,0.5),
              cex=LAB_CEX, family="serif") } } }

png(OUT, width = 17.0, height = 4.9, units = "in", res = 300, bg = "white")
layout(matrix(1:5, nrow = 1), widths = c(1, 1, 1, 1, 0.26))
nms <- names(panels)
for (q in seq_along(panels)){
   par(mar = if (LABEL_TOP > 0 && nms[q] == "VAPOR") c(0.8, 6.0, 3.2, 2.0)
       else c(0.8, 0.8, 3.2, 0.8), xpd = NA)
   draw_panel(panels[[q]], nms[q], annotate = (nms[q] == "VAPOR")) }

par(mar = c(4.2, 1.2, 5.4, 3.4), xpd = FALSE)
image(x = 1, y = seq(-1, 1, length.out = 101),
      z = matrix(seq(-1, 1, length.out = 101), 1),
      col = pal, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
axis(4, at = c(-1, 0, 1), labels = c("-1","0","1"), las = 1,
     cex.axis = 1.4, tick = FALSE, line = -0.5, family = "serif")
box(lwd = 0.6)
dev.off()
cat("wrote", OUT, "\n")