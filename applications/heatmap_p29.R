library(pheatmap)
library(grid)
library(gridExtra)

fit <- readRDS("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application/p29_full_results.rds")
OUT <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application/heatmap_p29_fullrank.png")

to_pc <- function(Om) { d <- sqrt(diag(Om)); pc <- -Om / outer(d, d); diag(pc) <- NA; pc }

panels <- list(
   "VAPOR"             = to_pc(fit$raw_omega$vapor),
   "Kastner"           = to_pc(fit$raw_omega$kast),
   "VAR + Glasso"      = to_pc(fit$raw_omega$glasso),
   "Sample covariance" = to_pc(fit$raw_omega$samp)
)

pal <- colorRampPalette(c("blue", "white", "red"))(100)

lim_all <- 1
brks    <- seq(-1, 1, length.out = 101)

mk <- function(M, title, show_legend = FALSE) {
   M[M >  lim_all] <-  lim_all; M[M < -lim_all] <- -lim_all
   gt <- pheatmap(M, main = title,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  color = pal, breaks = brks, na_col = "white",
                  border_color = NA,
                  show_rownames = FALSE, show_colnames = FALSE,
                  legend = show_legend,
                  legend_breaks = c(-1, 0, 1),
                  legend_labels = c("-1", "0", "1"),
                  silent = TRUE,
                  fontsize = 22)$gtable
   
   # find the title grob and set it to plain weight
   ti <- which(sapply(gt$grobs, function(g) "text" %in% class(g) && !is.null(g$label) &&
                         any(g$label == title)))
   if (length(ti)) {
      gt$grobs[[ti[1]]]$gp$fontface <- NULL
      gt$grobs[[ti[1]]]$gp$font <- 1L
      gt$grobs[[ti[1]]]$gp$fontfamily <- "serif"
   }
   gt
}

show  <- c(FALSE, FALSE, FALSE, TRUE)
grobs <- Map(mk, panels, names(panels), show)

png(OUT, width = 16.6, height = 4.6, units = "in", res = 300, bg = "white")
grid.arrange(grobs = grobs, ncol = 4, widths = c(1, 1, 1, 1.12))
dev.off()
cat("wrote", OUT, "\n")