# ============================================================================
#  Bayesian Precision Factor Analysis for Time Series
#  Application: FRED-MD Macroeconomic Data — Full Variable Set
#
#  Data: all available monthly macro variables, March 1992 to March 2026
#  Method: Joint sampler (proposed)
# ============================================================================

library(Rcpp)
library(RcppArmadillo)
library(PrecFact)
library(pheatmap)
library(MASS)

Rcpp::sourceCpp("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/var_prec.cpp")


# ============================================================================
#  SECTION 1: Load and preprocess FRED-MD
# ============================================================================

fred_raw <- read.csv("~/Desktop/2026-05-md.csv")

tcodes        <- as.numeric(fred_raw[1, -1])
names(tcodes) <- colnames(fred_raw)[-1]

fred_data        <- fred_raw[-1, ]
dates            <- fred_data[, 1]
fred_data        <- fred_data[, -1]
fred_data        <- apply(fred_data, 2, as.numeric)

transform_series <- function(x, tcode) {
   x <- as.numeric(x)
   switch(as.character(tcode),
          "1" = x,
          "2" = c(NA, diff(x)),
          "3" = c(NA, NA, diff(diff(x))),
          "4" = log(x),
          "5" = c(NA, diff(log(x))),
          "6" = c(NA, NA, diff(diff(log(x)))),
          "7" = c(NA, diff(log(x) / lag(log(x)) - 1)),
          x
   )
}

fred_trans           <- mapply(transform_series, as.data.frame(fred_data), tcodes)
fred_trans           <- as.matrix(fred_trans)
rownames(fred_trans) <- dates

fred_clean <- fred_trans[apply(fred_trans, 1, function(r) !any(is.na(r))), ]
fred_clean <- fred_clean[, apply(fred_clean, 2, function(c) !any(is.na(c)))]

cat("Dimensions after cleaning:", dim(fred_clean), "\n")
cat("Date range:", rownames(fred_clean)[1], "to",
    rownames(fred_clean)[nrow(fred_clean)], "\n")

# ============================================================================
#  SECTION 2: Prepare full variable set
# ============================================================================

var_sds   <- apply(fred_clean, 2, sd)
vars_keep <- names(var_sds[var_sds > 0.001])
Y_full    <- fred_clean[, vars_keep]
Y_scaled  <- scale(Y_full)
p         <- ncol(Y_scaled)
T_obs     <- nrow(Y_scaled)

cat("Variables after variance filter:", p, "\n")
cat("Final dimensions:", dim(Y_scaled), "\n")

set.seed(1)
Y_jit <- Y_scaled + matrix(rnorm(prod(dim(Y_scaled)), 0, 1e-6), T_obs, p)

# ============================================================================
#  SECTION 3: Joint sampler
# ============================================================================

cat("\nRunning joint sampler (T=", T_obs, ", p=", p, ")...\n")
res_joint <- var_precision_sampler(Y=Y_jit, lags=1, k=10,
                                   nrun=10000, burn=5000, thin=5)
out_joint <- fact_prec_processing(res_joint, burnout=1000, alpha=0.1,
                                  beta=0.85,
                                  epsilon_grid=seq(0.05, 0.1, length=20))
cat("pFDR:", out_joint$pFDR, "\n")
cat("Edges:", sum(out_joint$P)/2, "\n")

# ============================================================================
#  SECTION 4: Heatmap
# ============================================================================

pc_joint <- -cov2cor(out_joint$Omega * out_joint$P)
diag(pc_joint) <- NA
pc_joint <- round(pc_joint, 2)
rownames(pc_joint) <- colnames(pc_joint) <- vars_keep

# cluster by joint sampler ordering
hc        <- hclust(dist(replace(pc_joint, is.na(pc_joint), 0)))
var_order <- hc$order
pc_joint_r <- pc_joint[var_order, var_order]

pheatmap(pc_joint_r,
         main         = "Joint Sampler: FRED-MD Partial Correlation Matrix",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color        = colorRampPalette(c("blue","white","red"))(100),
         breaks       = seq(-1, 1, length.out=101),
         na_col       = "grey90",
         display_numbers = FALSE,
         fontsize     = 4)

# ============================================================================
#  SECTION 5: Edge summary and top edges
# ============================================================================

cat("\n============================================================\n")
cat("EDGE SUMMARY\n")
cat("============================================================\n")
cat("Joint sampler:", sum(out_joint$P)/2, "edges | pFDR:", out_joint$pFDR, "\n")

cat("\nTop 20 strongest edges:\n")
edges     <- which(out_joint$P == 1 & upper.tri(out_joint$P), arr.ind=TRUE)
edge_vals <- sapply(seq_len(nrow(edges)), function(i)
   abs(pc_joint[edges[i,1], edges[i,2]]))
top20 <- order(edge_vals, decreasing=TRUE)[1:min(20, length(edge_vals))]
for (i in top20) {
   r <- edges[i,1]; cc <- edges[i,2]
   cat(sprintf("  %-25s -- %-25s (%.2f)\n",
               vars_keep[r], vars_keep[cc], pc_joint[r, cc]))
}







