# ============================================================================
#  Item 4: the p = 100 specific tables.  NO REFIT: every number comes from
#  p100_results.rds, which already carries the one-convention PCs.
#     Rscript applications/final_p100_tables.R
#  Writes p100_convention.rds and the regenerated heatmap.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 210)
R <- readRDS(file.path(FITS_DIR, "p100_results.rds"))
vars <- R$vars; p <- length(vars); lab <- pair_labels(vars)
ut <- which(upper.tri(matrix(0,p,p)), arr.ind = TRUE)
OUT <- list()

PCs <- R$long$pc_all
stopifnot(!is.null(PCs))
K <- readRDS(file.path(FITS_DIR, "kastner_p100_20k_draws.rds"))
kast_src <- "mean of per-draw partial correlations"; kast_n <- 4000L
prov <- do.call(rbind, lapply(names(PCs), function(nm){ z <- PCs[[nm]]
   data.frame(method = nm,
      estimate = if (nm == "Kastner") kast_src
                 else if (is.null(attr(z,"n_draws"))) "point estimate"
                 else "mean of per-draw partial correlations",
      chains = if (nm == "Kastner") length(K$chains) else
               if (is.null(attr(z,"n_chains"))) 0L else attr(z,"n_chains"),
      draws = if (nm == "Kastner") kast_n else
              if (is.null(attr(z,"n_draws"))) NA_integer_ else attr(z,"n_draws")) }))
cat("=== 24. p = 100 partial-correlation convention, per method ===\n")
print(prov, row.names = FALSE); OUT$provenance <- prov

thr <- threshold_counts(PCs, R$config$THRESHOLDS)
cat("\n=== 24. p = 100 threshold edge counts, one convention ===\n")
print(thr, row.names = FALSE); OUT$thresholds <- thr

ed <- R$long$ci_edges
rk <- edge_recovery_ranks(ed, PCs, vars)
rownames(rk) <- paste(ed$v1, ed$v2, sep = "--")
cat(sprintf("\n=== 24. p = 100 rank of each of VAPOR's %d edges, one convention ===\n", nrow(rk)))
print(rk); OUT$ranks <- rk
TOPS <- c(25L, 100L, 250L)
summ <- do.call(rbind, lapply(colnames(rk), function(m){ v <- rk[,m]
   setNames(data.frame(m, round(median(v),1), sum(v <= TOPS[1]), sum(v <= TOPS[2]),
                       sum(v <= TOPS[3])),
            c("method","median_rank", sprintf("in_top_%d", TOPS))) }))
cat("\n=== 24. p = 100 edge recovery summary, one convention ===\n")
print(summ, row.names = FALSE); OUT$summary <- summ

pur <- group_purity(PCs, vars, R$config$GROUPS, 100L)
cat("\n=== 24. p = 100 group purity of the top 100 pairs, one convention ===\n")
print(pur, row.names = FALSE); OUT$purity <- pur

#  the ranks previously quoted
RATE <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTGAGE30US")
is_rate <- ed$v1 %in% RATE & ed$v2 %in% RATE
neg <- which(ed$pcor < 0)
q1 <- rk[is_rate, , drop = FALSE]
cat("\n=== 24. ranks of the rate edges, every method ===\n"); print(q1); OUT$rate_edges <- q1
q2 <- rk[neg, , drop = FALSE]
cat(sprintf("\n=== 24. ranks of the %d negative edges, every method ===\n", length(neg)))
print(q2); OUT$negative_edges <- q2
g5 <- grep("^GS5--MORTGAGE30US$|^MORTGAGE30US--GS5$", rownames(rk))
if (length(g5)){
   cat("\n=== 24. GS5--MORTGAGE30US, every method ===\n")
   print(rk[g5, , drop = FALSE]); OUT$gs5_mortgage <- rk[g5, , drop = FALSE] }
#  Kastner's exact chain length (item 5, no refit)
kcfg <- data.frame(fit = "kastner_p100_20k_draws.rds",
   total_iterations = K$config$ITER,
   draws_arg = K$config$NRUN, burnin = K$config$BURN, thin = K$config$THIN,
   retained_per_chain = paste(unique(vapply(K$chains, function(z) z$n_draws, integer(1))), collapse=","),
   chains = length(K$chains), factor_factors = K$config$K)
cat("\n=== 5. the existing p = 100 Kastner 20k fit, exact settings ===\n")
print(kcfg, row.names = FALSE); OUT$kastner_cfg <- kcfg
OUT$PCs <- PCs            # the matrices themselves, for the heatmap panel
saveRDS(OUT, file.path(FITS_DIR, "p100_tables_extra.rds"))
cat("\nsaved p100_tables_extra.rds\n")
