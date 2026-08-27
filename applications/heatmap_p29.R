# ============================================================================
#  p=29 partial-correlation heatmaps.
#  Boxes appear on the VAPOR panel only: the 95% credible-interval rule needs
#  a posterior, and only VAPOR supplies one. Each boxed cell carries the
#  column variable above it and the row variable to its left.
# ============================================================================
APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
fit <- readRDS(file.path(APP, "p29_full_results.rds"))
OUT <- file.path(APP, "heatmap_p29_r3k10.png")

LAB_CEX <- 1.15          # raise or lower to taste

vars <- fit$vars; n <- length(vars)
ed   <- fit$ci_edges
ei   <- cbind(match(ed$v1, vars), match(ed$v2, vars))
cat(sprintf("%d reported edges\n", nrow(ei)))

to_pc <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- NA; pc }
panels <- list("VAPOR"             = { M <- fit$pc_vapor; diag(M) <- NA; M },
               "Kastner"           = to_pc(fit$raw_omega$kast),
               "VAR + Glasso"      = to_pc(fit$raw_omega$glasso),
               "Sample covariance" = to_pc(fit$raw_omega$samp))
pal <- colorRampPalette(c("blue","white","red"))(101)

#  Draw each edge once, in the upper triangle, so labels do not collide with
#  the mirrored cell. Cell (i,j) maps to x = j, y = n - i + 1.
up <- t(apply(ei, 1, function(z) c(min(z), max(z))))   # row = smaller index

draw_panel <- function(M, title, annotate = FALSE){
   M[M >  1] <-  1; M[M < -1] <- -1
   image(x = 1:n, y = 1:n, z = t(M[n:1, ]), col = pal, zlim = c(-1,1),
         axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
   title(main = title, family = "serif", font.main = 1, cex.main = 2.0, line = 0.6)
   if (!annotate) return(invisible())
   for (e in seq_len(nrow(up))){
      i <- up[e,1]; j <- up[e,2]                       # i < j: upper triangle
      cx <- j; cy <- n - i + 1
      rect(cx-0.5, cy-0.5, cx+0.5, cy+0.5, border="black", lwd=2.4)
      #  column variable, above the box
      text(cx, cy + 0.8, vars[j], adj = c(0.5, 0), cex = LAB_CEX,
           family = "serif", col = "black")
      #  row variable, to the left of the box
      text(cx - 0.8, cy, vars[i], adj = c(1, 0.5), cex = LAB_CEX,
           family = "serif", col = "black") } }

png(OUT, width = 19.0, height = 5.4, units = "in", res = 300, bg = "white")
layout(matrix(1:5, nrow = 1), widths = c(1.30, 1, 1, 1, 0.26))
nms <- names(panels)
for (q in seq_along(panels)){
   par(mar = if (nms[q] == "VAPOR") c(0.8, 7.0, 3.2, 2.2) else c(0.8, 0.8, 3.2, 0.8),
       xpd = NA)
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