# ============================================================================
#  Items 2c and 3: all pairs ranked by |pooled mean per-draw partial
#  correlation|, and the rate-series block.
#     Rscript applications/final_pairs.R <29|100>
#  Reads p29_results.rds; fits nothing.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)
DIM <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (!DIM %in% c(29L, 100L)) stop("usage: Rscript final_pairs.R <29|100>")
f <- readRDS(file.path(FITS_DIR, sprintf("p%d_results.rds", DIM)))
vars <- f$vars; p <- length(vars); lab <- pair_labels(vars)
DR <- lapply(f$draws$vapor, pair_draws, p = p)
pool <- do.call(rbind, DR)
sel <- ci_indices(pool)
A <- data.frame(pair = lab, mean = round(colMeans(pool), 4),
                lo = round(apply(pool, 2, quantile, .025), 4),
                hi = round(apply(pool, 2, quantile, .975), 4),
                selected = seq_along(lab) %in% sel)
A$rank <- rank(-abs(A$mean), ties.method = "first")
A <- A[order(A$rank), ]
cat(sprintf("p = %d | %d chains, %d pooled draws | %d of %d pairs selected\n", DIM,
            length(DR), nrow(pool), length(sel), nrow(A)))
cat(sprintf("\n=== 2c. p = %d ranks 1-10 of all %d pairs by |pooled mean per-draw pcor| ===\n", DIM, nrow(A)))
print(A[1:10, c("rank","pair","mean","lo","hi","selected")], row.names = FALSE)
G <- if (DIM == 29L) c("FEDFUNDS","TB3MS","GS1","GS10") else
     c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
want <- apply(t(combn(G,2)), 1, function(z){ i <- match(z, vars)
   paste(vars[min(i)], vars[max(i)], sep="--") })
cat(sprintf("\n=== p = %d: all %d pairs among %s ===\n", DIM, choose(length(G),2), paste(G, collapse=", ")))
B <- A[match(want, A$pair), c("pair","mean","lo","hi","selected","rank")]
print(B, row.names = FALSE)
saveRDS(list(all_pairs = A, top10 = A[1:10,], rates = B),
        file.path(FITS_DIR, sprintf("pairs_p%d.rds", DIM)))
cat(sprintf("\nsaved pairs_p%d.rds\n", DIM))
