# ============================================================================
#  Item 4, NIPA block: the per-chain behaviour of the NIPA partial
#  correlations in the new p = 100 fit, and the accounting identity in the
#  training panel.
#     Rscript applications/final_nipa_p100.R
#  Read-only; refits nothing.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 205)
sink(file.path(RES_DIR, "logs", "nipa_p100.txt"), split = TRUE)
R <- readRDS(file.path(FITS_DIR, "p100_results.rds"))
S <- c("GDPC1","PCECC96","GPDIC1","GCEC1","EXPGSC1","IMPGSC1")
stopifnot(all(S %in% R$vars))
Y <- R$Y_train[, S, drop = FALSE]
cat(sprintf("p = 100 training panel, %d rows (%s), transformed and standardized on the training window\n",
            nrow(Y), paste(range(rownames(Y)), collapse = " .. ")))
cat(sprintf("series: %s\n", paste(S, collapse = ", ")))
cat(sprintf("column means %.2e to %.2e | column SDs %.4f to %.4f\n",
            min(colMeans(Y)), max(colMeans(Y)), min(apply(Y,2,sd)), max(apply(Y,2,sd))))


# ---------------------------------------------------------------------------
#  per-chain NIPA partial correlations in the new fit
# ---------------------------------------------------------------------------
vars <- R$vars; pp <- length(vars); lab <- pair_labels(vars)
ch <- R$draws$vapor; seeds <- vapply(ch, function(f) attr(f,"seed"), integer(1))
DR <- lapply(ch, pair_draws, p = pp)
PAIRS <- c("GDPC1--GPDIC1","GDPC1--IMPGSC1","GDPC1--EXPGSC1","GPDIC1--IMPGSC1")
jj <- match(PAIRS, lab); stopifnot(!anyNA(jj))
cat("\n=== 4. p = 100 NIPA block, per chain ===\n")
T0 <- do.call(rbind, lapply(seq_along(ch), function(c)
   setNames(data.frame(seeds[c], t(round(colMeans(DR[[c]][, jj, drop = FALSE]), 4)),
                       round(mean(DR[[c]][, jj[1]] > 0.5), 3)),
            c("seed", PAIRS, "frac_GDPC1_GPDIC1_gt_0.5"))))
print(T0, row.names = FALSE)
cat("\n  pooled R-hat and bulk ESS:\n")
T0b <- do.call(rbind, lapply(seq_along(jj), function(i){
   M <- do.call(cbind, lapply(DR, function(X) X[, jj[i]]))
   data.frame(pair = PAIRS[i], pooled_mean = round(mean(M), 4),
              Rhat = round(posterior::rhat(M), 3),
              bulk_ESS = round(posterior::ess_bulk(M), 0),
              tail_ESS = round(posterior::ess_tail(M), 0)) }))
print(T0b, row.names = FALSE)
cat("\n=== 4. p = 100 sample correlations among the five rate series ===\n")
R5 <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
CR <- cor(R$Y_train[, R5]); pw <- CR[upper.tri(CR)]
print(data.frame(n_pairs = length(pw), min = round(min(pw),3),
                 median = round(median(pw),3), max = round(max(pw),3)), row.names = FALSE)
OUTX <- list(nipa_per_chain = T0, nipa_pooled = T0b,
             rate_cor = CR, rate_summary = c(min(pw), median(pw), max(pw)))

cat("\n=== 4. GDPC1 on the other five (accounting identity) ===\n")
d <- as.data.frame(Y)
m <- lm(GDPC1 ~ PCECC96 + GPDIC1 + GCEC1 + EXPGSC1 + IMPGSC1, data = d)
sm <- summary(m)
T1 <- data.frame(term = names(coef(m)), estimate = round(coef(m), 4),
                 std_error = round(sm$coefficients[, 2], 4),
                 t = round(sm$coefficients[, 3], 2),
                 p = signif(sm$coefficients[, 4], 3), row.names = NULL)
print(T1, row.names = FALSE)
T2 <- data.frame(n = nobs(m), df_residual = df.residual(m),
                 R2 = round(sm$r.squared, 6), adj_R2 = round(sm$adj.r.squared, 6),
                 residual_SD = round(sm$sigma, 6),
                 residual_SD_as_share_of_GDPC1_SD = round(sm$sigma/sd(d$GDPC1), 6),
                 F_stat = round(sm$fstatistic[1], 1))
cat("\n"); print(T2, row.names = FALSE)

cat("\n=== 2b. the 6 x 6 sample correlation matrix ===\n")
C <- cor(Y); print(round(C, 3))
ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
cat("\n  all six eigenvalues:\n")
print(data.frame(index = 1:6, eigenvalue = signif(ev, 6)), row.names = FALSE)
T3 <- data.frame(rank = c("smallest","2nd smallest","3rd smallest"),
                 eigenvalue = signif(rev(ev)[1:3], 6))
cat("\n  smallest three:\n"); print(T3, row.names = FALSE)
cat(sprintf("\n  condition number (largest / smallest eigenvalue): %.2f\n", ev[1]/ev[6]))
cat(sprintf("  kappa(C) by the 2-norm: %.2f | determinant: %.4e\n",
            kappa(C, exact = TRUE), det(C)))
saveRDS(c(OUTX, list(coefficients = T1, fit = T2, correlation = C, eigenvalues = ev,
          smallest_three = T3, condition = ev[1]/ev[6])),
        file.path(FITS_DIR, "nipa_p100.rds"))
sink()
