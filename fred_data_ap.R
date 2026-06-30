# ============================================================================
#  Bayesian Precision Factor Analysis for Time Series
#  Application: FRED-MD Macroeconomic Data
#
#  Data: 29 monthly macro variables, March 1992 to March 2026 (T=391)
#  Methods compared:
#    1. PrecFact (no VAR)
#    2. Glasso (no VAR)
#    3. VAR + Glasso
#    4. VAR + PrecFact
#    5. Joint sampler (proposed)
# ============================================================================

library(Rcpp)
library(RcppArmadillo)
library(PrecFact)
library(bayesianVARs)
library(glasso)
library(pheatmap)
library(MASS)

Rcpp::sourceCpp("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/var_prec.cpp")

# ============================================================================
#  SECTION 1: Load and preprocess FRED-MD
# ============================================================================

fred_raw <- read.csv("~/Desktop/2026-05-md.csv")

# extract transformation codes (row 1)
tcodes        <- as.numeric(fred_raw[1, -1])
names(tcodes) <- colnames(fred_raw)[-1]

# remove transformation row
fred_data        <- fred_raw[-1, ]
dates            <- fred_data[, 1]
fred_data        <- fred_data[, -1]
fred_data        <- apply(fred_data, 2, as.numeric)

# apply McCracken-Ng stationarity transformations
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

# remove rows and columns with any NAs
fred_clean <- fred_trans[apply(fred_trans, 1, function(r) !any(is.na(r))), ]
fred_clean <- fred_clean[, apply(fred_clean, 2, function(c) !any(is.na(c)))]

cat("Dimensions after cleaning:", dim(fred_clean), "\n")
cat("Date range:", rownames(fred_clean)[1], "to",
    rownames(fred_clean)[nrow(fred_clean)], "\n")

# ============================================================================
#  SECTION 2: Select 29 representative variables
# ============================================================================

vars_30 <- c(
   # output and production
   "INDPRO", "IPFINAL", "IPCONGD", "IPBUSEQ", "CUMFNS",
   # employment
   "PAYEMS", "UNEMRATE", "UEMPMEAN", "CLAIMSx", "HWI", "AWOTMAN",
   # consumption and orders
   "DPCERA3M086SBEA", "CMRMTSPLx", "AMDMNOx", "BUSINVx",
   # prices
   "CPIAUCSL", "PPICMM", "PCEPI", "OILPRICEx", "WPSFD49207",
   # interest rates
   "FEDFUNDS", "TB3MS", "GS1", "GS10", "BAA",
   # money and credit
   "M2SL", "BUSLOANS", "NONREVSL",
   # financial markets
   "S.P.500", "VIXCLSx"
)

vars_30    <- vars_30[vars_30 %in% colnames(fred_clean)]
Y_fred     <- fred_clean[, vars_30]
Y_scaled   <- scale(Y_fred)

cat("Variables selected:", length(vars_30), "\n")
cat("Final dimensions:", dim(Y_scaled), "\n")

# readable variable labels
var_labels <- c(
   IPFINAL           = "IP: Final Products",
   TB3MS             = "3-Month T-Bill",
   GS10              = "10-Year Treasury",
   IPCONGD           = "IP: Consumer Goods",
   IPBUSEQ           = "IP: Business Equip",
   FEDFUNDS          = "Fed Funds Rate",
   GS1               = "1-Year Treasury",
   BAA               = "BAA Corp Bond",
   S.P.500           = "S&P 500",
   DPCERA3M086SBEA   = "Real Consumption",
   VIXCLSx           = "VIX",
   NONREVSL          = "Nonrevolving Credit",
   BUSLOANS          = "Business Loans",
   M2SL              = "M2 Money Supply",
   OILPRICEx         = "Oil Price",
   PPICMM            = "PPI: Metals",
   BUSINVx           = "Business Inventories",
   AMDMNOx           = "Durable Goods Orders",
   CMRMTSPLx         = "Retail Sales",
   AWOTMAN           = "Overtime: Mfg",
   HWI               = "Help Wanted Index",
   CLAIMSx           = "Unemployment Claims",
   PAYEMS            = "Nonfarm Payrolls",
   UEMPMEAN          = "Avg Unemployment Duration",
   INDPRO            = "Industrial Production",
   CUMFNS            = "Capacity Utilization",
   CPIAUCSL          = "CPI",
   PCEPI             = "PCE Inflation",
   WPSFD49207        = "PPI: Finished Goods"
)

# ============================================================================
#  SECTION 3: Fit all methods
# ============================================================================

# ── 1. PrecFact (no VAR) ─────────────────────────────────────────────────────
cat("\n1. PrecFact (no VAR)...\n")
res_pf_raw <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=5,
                       nrun=10000, thin=5, Y=Y_scaled,
                       alldiff=0, a_dir=0.1, b_dir=0.1)
out_pf_raw <- fact_prec_processing(res_pf_raw, burnout=1000, alpha=0.1,
                                   beta=0.85,
                                   epsilon_grid=seq(0.05, 0.1, length=20))
cat("  pFDR:", out_pf_raw$pFDR, "| Edges:", sum(out_pf_raw$P)/2, "\n")

# ── 2. Glasso (no VAR) ───────────────────────────────────────────────────────
cat("\n2. Glasso (no VAR)...\n")
S_raw    <- cov(Y_scaled)
rho_grid <- 10^seq(-2, 0, length=30)
gl_scores <- sapply(rho_grid, function(rho) {
   gl <- tryCatch(glasso(S_raw, rho=rho, penalize.diagonal=FALSE),
                  error=function(e) NULL)
   if (is.null(gl)) return(Inf)
   sum(gl$wi != 0) / 2   # number of edges as proxy
})
# pick rho giving ~35 edges (matching joint sampler) for fair visual comparison
rho_glasso_raw <- rho_grid[which.min(abs(gl_scores - 35))]
gl_raw   <- glasso(S_raw, rho=rho_glasso_raw, penalize.diagonal=FALSE)
out_gl_raw <- list(Omega=gl_raw$wi, P=(gl_raw$wi != 0)*1)
cat("  rho:", round(rho_glasso_raw, 4), "| Edges:", sum(out_gl_raw$P)/2, "\n")

# ── 3. VAR + Glasso ──────────────────────────────────────────────────────────
cat("\n3. VAR + Glasso...\n")
# add tiny jitter to handle near-collinearity in FRED-MD
set.seed(1)
Y_scaled_jit <- Y_scaled + matrix(rnorm(prod(dim(Y_scaled)), 0, 1e-6),
                                  nrow(Y_scaled), ncol(Y_scaled))
mod_var   <- bvar(Y_scaled_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
B_hat     <- apply(mod_var$PHI, c(1,2), mean)
resid_var <- mod_var$Y - mod_var$X %*% B_hat
S_resid  <- cov(resid_var)
gl_scores_var <- sapply(rho_grid, function(rho) {
   gl <- tryCatch(glasso(S_resid, rho=rho, penalize.diagonal=FALSE),
                  error=function(e) NULL)
   if (is.null(gl)) return(Inf)
   sum(gl$wi != 0) / 2
})
rho_glasso_var <- rho_grid[which.min(abs(gl_scores_var - 35))]
gl_var   <- glasso(S_resid, rho=rho_glasso_var, penalize.diagonal=FALSE)
out_gl_var <- list(Omega=gl_var$wi, P=(gl_var$wi != 0)*1)
cat("  rho:", round(rho_glasso_var, 4), "| Edges:", sum(out_gl_var$P)/2, "\n")

# ── 4. VAR + PrecFact ────────────────────────────────────────────────────────
cat("\n4. VAR + PrecFact...\n")
res_pf_var <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=5,
                       nrun=10000, thin=5, Y=scale(resid_var),
                       alldiff=0, a_dir=0.1, b_dir=0.1)
out_pf_var <- fact_prec_processing(res_pf_var, burnout=1000, alpha=0.1,
                                   beta=0.85,
                                   epsilon_grid=seq(0.05, 0.1, length=20))
cat("  pFDR:", out_pf_var$pFDR, "| Edges:", sum(out_pf_var$P)/2, "\n")

# ── 5. Joint sampler ─────────────────────────────────────────────────────────
cat("\n5. Joint sampler...\n")
res_joint <- var_precision_sampler(Y=Y_scaled_jit, lags=1, k=5,
                                   nrun=10000, burn=5000, thin=5)
out_joint <- fact_prec_processing(res_joint, burnout=1000, alpha=0.1,
                                  beta=0.85,
                                  epsilon_grid=seq(0.05, 0.1, length=20))
cat("  pFDR:", out_joint$pFDR, "| Edges:", sum(out_joint$P)/2, "\n")

# ============================================================================
#  SECTION 4: Heatmap comparison
# ============================================================================

make_partcor <- function(Omega, P=NULL) {
   if (!is.null(P)) Omega <- Omega * P
   m <- -cov2cor(Omega)
   diag(m) <- NA
   round(m, 2)
}

apply_labels <- function(mat) {
   rownames(mat) <- var_labels[rownames(mat)]
   colnames(mat) <- var_labels[colnames(mat)]
   mat
}

col    <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-1, 1, length.out=101)
hmap   <- function(mat, title) {
   pheatmap(apply_labels(mat),
            main            = title,
            cluster_rows    = TRUE,
            cluster_cols    = TRUE,
            color           = col,
            breaks          = breaks,
            na_col          = "grey90",
            display_numbers = mat,
            number_color    = "black",
            fontsize_number = 6,
            fontsize        = 8)
}

pc_pf_raw  <- make_partcor(out_pf_raw$Omega, out_pf_raw$P)
pc_gl_raw  <- make_partcor(out_gl_raw$Omega, out_gl_raw$P)
pc_gl_var  <- make_partcor(out_gl_var$Omega, out_gl_var$P)
pc_pf_var  <- make_partcor(out_pf_var$Omega, out_pf_var$P)
pc_joint   <- make_partcor(out_joint$Omega,  out_joint$P)

rownames(pc_pf_raw) <- colnames(pc_pf_raw) <- vars_30
rownames(pc_gl_raw) <- colnames(pc_gl_raw) <- vars_30
rownames(pc_gl_var) <- colnames(pc_gl_var) <- vars_30
rownames(pc_pf_var) <- colnames(pc_pf_var) <- vars_30
rownames(pc_joint)  <- colnames(pc_joint)  <- vars_30

hmap(pc_pf_raw, "PrecFact (no VAR)")
hmap(pc_gl_raw, "Glasso (no VAR)")
hmap(pc_gl_var, "VAR + Glasso")
hmap(pc_pf_var, "VAR + PrecFact")
hmap(pc_joint,  "Joint Sampler (proposed)")

# ============================================================================
#  SECTION 5: Edge summary
# ============================================================================

cat("\n============================================================\n")
cat("EDGE COUNT SUMMARY\n")
cat("============================================================\n")
cat("PrecFact (no VAR): ", sum(out_pf_raw$P)/2, "edges | pFDR:", out_pf_raw$pFDR, "\n")
cat("Glasso (no VAR):   ", sum(out_gl_raw$P)/2, "edges\n")
cat("VAR + Glasso:      ", sum(out_gl_var$P)/2, "edges\n")
cat("VAR + PrecFact:    ", sum(out_pf_var$P)/2, "edges | pFDR:", out_pf_var$pFDR, "\n")
cat("Joint sampler:     ", sum(out_joint$P)/2,  "edges | pFDR:", out_joint$pFDR,  "\n")

# edges selected by joint sampler
cat("\nEdges selected by joint sampler:\n")
edges <- which(out_joint$P == 1 & upper.tri(out_joint$P), arr.ind=TRUE)
for (i in seq_len(nrow(edges))) {
   r <- edges[i, 1]; cc <- edges[i, 2]
   v1 <- var_labels[vars_30[r]]
   v2 <- var_labels[vars_30[cc]]
   cat(sprintf("  %s -- %s (%.2f)\n", v1, v2,
               pc_joint[r, cc]))
}




# get consistent ordering from joint sampler clustering
hc <- hclust(dist(pc_joint))
var_order <- hc$order
vars_ordered <- vars_30[var_order]

# reorder all matrices to match
reorder_mat <- function(mat) mat[var_order, var_order]

pc_pf_raw_r <- reorder_mat(pc_pf_raw)
pc_gl_raw_r <- reorder_mat(pc_gl_raw)
pc_gl_var_r <- reorder_mat(pc_gl_var)
pc_pf_var_r <- reorder_mat(pc_pf_var)
pc_joint_r  <- reorder_mat(pc_joint)

# heatmap with fixed ordering
hmap_fixed <- function(mat, title) {
   m <- mat
   rownames(m) <- var_labels[rownames(m)]
   colnames(m) <- var_labels[colnames(m)]
   pheatmap(m,
            main            = title,
            cluster_rows    = FALSE,
            cluster_cols    = FALSE,
            color           = colorRampPalette(c("blue","white","red"))(100),
            breaks          = seq(-1, 1, length.out=101),
            na_col          = "grey90",
            display_numbers = m,
            number_color    = "black",
            fontsize_number = 6,
            fontsize        = 8)
}

hmap_fixed(pc_pf_raw_r, "PrecFact (no VAR)")
hmap_fixed(pc_gl_raw_r, "Glasso (no VAR)")
hmap_fixed(pc_gl_var_r, "VAR + Glasso")
hmap_fixed(pc_pf_var_r, "VAR + PrecFact")
hmap_fixed(pc_joint_r,  "Joint Sampler (proposed)")






library(gridExtra)

p1 <- pheatmap(apply_labels(pc_pf_var_r),
               main            = "VAR + PrecFact",
               cluster_rows    = FALSE,
               cluster_cols    = FALSE,
               color           = colorRampPalette(c("blue","white","red"))(100),
               breaks          = seq(-1, 1, length.out=101),
               na_col          = "grey90",
               display_numbers = round(apply_labels(pc_pf_var_r), 2),
               number_color    = "black",
               fontsize_number = 5,
               fontsize        = 7,
               silent          = TRUE)

p2 <- pheatmap(apply_labels(pc_joint_r),
               main            = "Joint Sampler (proposed)",
               cluster_rows    = FALSE,
               cluster_cols    = FALSE,
               color           = colorRampPalette(c("blue","white","red"))(100),
               breaks          = seq(-1, 1, length.out=101),
               na_col          = "grey90",
               display_numbers = round(apply_labels(pc_joint_r), 2),
               number_color    = "black",
               fontsize_number = 5,
               fontsize        = 7,
               silent          = TRUE)

grid.arrange(p1$gtable, p2$gtable, ncol=2)
