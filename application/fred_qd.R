# ============================================================================
#  VAPOR: Bayesian Precision Factor Analysis for Vector Autoregressions
#  Application: FRED-QD Quarterly Macroeconomic Data
#
#  Data: 29 quarterly macro variables, Q2 2000 to Q4 2025 (T=97)
#  Train: Q2 2000 to Q4 2020 (T=77)
#  Test:  Q1 2021 to Q4 2025 (T=20)
#
#  Evaluation: held-out prediction MSPE
#    At each test quarter, hold out 30% of variables randomly
#    Predict held-out variables from observed using estimated precision matrix
#
#  Methods:
#    1. PrecFact (no VAR)
#    2. VAR + Glasso
#    3. VAR + PrecFact (sequential)
#    4. Joint sampler / VAPOR (proposed)
#    5. Kastner (VAR + factor SV)
#    6. Sample covariance baseline
#    7. Random walk baseline
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
#  SECTION 1: Load and preprocess FRED-QD
# ============================================================================

fred_qd_raw <- read.csv("~/Desktop/2026-05-qd.csv")

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

# row 1 is "factors", row 2 is "transform"
tcodes_qd        <- as.numeric(fred_qd_raw[2, -1])
names(tcodes_qd) <- colnames(fred_qd_raw)[-1]
fred_qd_data     <- fred_qd_raw[-c(1,2), ]
dates_qd         <- fred_qd_data[, 1]
fred_qd_data     <- fred_qd_data[, -1]
fred_qd_data     <- apply(fred_qd_data, 2, as.numeric)

fred_qd_trans           <- mapply(transform_series,
                                  as.data.frame(fred_qd_data), tcodes_qd)
fred_qd_trans           <- as.matrix(fred_qd_trans)
rownames(fred_qd_trans) <- dates_qd

fred_qd_clean <- fred_qd_trans[apply(fred_qd_trans, 1, function(r) !any(is.na(r))), ]
fred_qd_clean <- fred_qd_clean[, apply(fred_qd_clean, 2, function(c) !any(is.na(c)))]

cat("Dimensions after cleaning:", dim(fred_qd_clean), "\n")
cat("Date range:", rownames(fred_qd_clean)[1], "to",
    rownames(fred_qd_clean)[nrow(fred_qd_clean)], "\n")

# ============================================================================
#  SECTION 2: Select 29 representative variables
# ============================================================================

vars_qd <- c(
   "GDPC1", "INDPRO", "IPFINAL", "IPCONGD", "IPBUSEQ",
   "PAYEMS", "UEMPMEAN", "CLAIMSx", "CUMFNS",
   "PCECC96", "CMRMTSPLx", "BUSINVx",
   "CPIAUCSL", "PPICMM", "GDPCTPI", "OILPRICEx",
   "FEDFUNDS", "GS1", "GS10", "BAA", "TB3MS",
   "M2REAL", "BUSLOANSx", "NONREVSLx",
   "S.P.500", "VIXCLSx",
   "GPDIC1", "UNRATE", "WPSFD49207"
)

vars_qd     <- vars_qd[vars_qd %in% colnames(fred_qd_clean)]
Y_qd        <- fred_qd_clean[, vars_qd]
Y_qd_scaled <- scale(Y_qd)
p           <- ncol(Y_qd_scaled)
T_obs       <- nrow(Y_qd_scaled)

cat("Variables selected:", p, "\n")
cat("Final dimensions:", dim(Y_qd_scaled), "\n")
cat("T/p ratio:", round(T_obs/p, 2), "\n")

var_labels_qd <- c(
   GDPC1      = "Real GDP",
   INDPRO     = "Industrial Production",
   IPFINAL    = "IP: Final Products",
   IPCONGD    = "IP: Consumer Goods",
   IPBUSEQ    = "IP: Business Equip",
   PAYEMS     = "Nonfarm Payrolls",
   UEMPMEAN   = "Avg Unemployment Duration",
   CLAIMSx    = "Unemployment Claims",
   CUMFNS     = "Capacity Utilization",
   PCECC96    = "Real PCE",
   CMRMTSPLx  = "Retail Sales",
   BUSINVx    = "Business Inventories",
   CPIAUCSL   = "CPI",
   PPICMM     = "PPI: Metals",
   GDPCTPI    = "GDP Deflator",
   OILPRICEx  = "Oil Price",
   FEDFUNDS   = "Fed Funds Rate",
   GS1        = "1-Year Treasury",
   GS10       = "10-Year Treasury",
   BAA        = "BAA Corp Bond",
   TB3MS      = "3-Month T-Bill",
   M2REAL     = "Real M2",
   BUSLOANSx  = "Business Loans",
   NONREVSLx  = "Nonrevolving Credit",
   S.P.500    = "S&P 500",
   VIXCLSx    = "VIX",
   GPDIC1     = "Real Gross Investment",
   UNRATE     = "Unemployment Rate",
   WPSFD49207 = "PPI: Finished Goods"
)

# ============================================================================
#  SECTION 3: Train/test split
# ============================================================================

T_train     <- floor(0.8 * T_obs)
T_test      <- T_obs - T_train
Y_train     <- Y_qd_scaled[1:T_train, ]
Y_test      <- Y_qd_scaled[(T_train+1):T_obs, ]

cat("\nTrain:", T_train, "quarters |",
    rownames(Y_qd_scaled)[1], "to",
    rownames(Y_qd_scaled)[T_train], "\n")
cat("Test: ", T_test, "quarters |",
    rownames(Y_qd_scaled)[T_train+1], "to",
    rownames(Y_qd_scaled)[T_obs], "\n")

set.seed(1)
Y_train_jit <- Y_train + matrix(rnorm(prod(dim(Y_train)), 0, 1e-3),
                                T_train, p)

# last training observation for random walk baseline
Y_last_train <- Y_train[T_train, ]

# ============================================================================
#  SECTION 4: Fit all methods on training data
# ============================================================================

# ── 1. PrecFact (no VAR) ─────────────────────────────────────────────────────
cat("\n1. PrecFact (no VAR)...\n")
res_pf <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=5,
                   nrun=7500, thin=5, Y=Y_train,
                   alldiff=0, a_dir=0.1, b_dir=0.1)
out_pf <- fact_prec_processing(res_pf, burnout=750, alpha=0.1,
                               beta=0.85,
                               epsilon_grid=seq(0.05,0.1,length=20))
cat("  pFDR:", out_pf$pFDR, "| Edges:", sum(out_pf$P)/2, "\n")

# ── 2. VAR + Glasso ──────────────────────────────────────────────────────────
cat("\n2. VAR + Glasso...\n")
mod_var   <- bvar(Y_train_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
B_hat     <- apply(mod_var$PHI, c(1,2), mean)
resid_var <- mod_var$Y - mod_var$X %*% B_hat

# cross-validated glasso: use BIC-like criterion over rho grid
S_resid   <- cov(scale(resid_var))
rho_grid  <- 10^seq(-2, 0, length=30)
bic_vals  <- sapply(rho_grid, function(rho) {
   gl <- tryCatch(glasso(S_resid, rho=rho, penalize.diagonal=FALSE),
                  error=function(e) NULL)
   if (is.null(gl)) return(Inf)
   # BIC: -2*loglik + log(T)*edges
   Omega_gl <- gl$wi
   n_edges  <- sum(Omega_gl[upper.tri(Omega_gl)] != 0)
   loglik   <- (T_train/2) * (log(det(Omega_gl)) -
                                 sum(diag(S_resid %*% Omega_gl)))
   -2*loglik + log(T_train) * n_edges
})
rho_bic    <- rho_grid[which.min(bic_vals)]
gl_var     <- glasso(S_resid, rho=rho_bic, penalize.diagonal=FALSE)
Omega_glasso <- gl_var$wi
cat("  rho (BIC):", round(rho_bic, 4),
    "| Edges:", sum(gl_var$wi[upper.tri(gl_var$wi)] != 0), "\n")

# ── 3. VAR + PrecFact (sequential) ───────────────────────────────────────────
cat("\n3. VAR + PrecFact (sequential)...\n")
res_pf_var <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=5,
                       nrun=7500, thin=5, Y=scale(resid_var),
                       alldiff=0, a_dir=0.1, b_dir=0.1)
out_pf_var <- fact_prec_processing(res_pf_var, burnout=750, alpha=0.1,
                                   beta=0.85,
                                   epsilon_grid=seq(0.05,0.1,length=20))
cat("  pFDR:", out_pf_var$pFDR, "| Edges:", sum(out_pf_var$P)/2, "\n")

# ── 4. Joint sampler (VAPOR) ─────────────────────────────────────────────────
cat("\n4. Joint sampler (VAPOR)...\n")
res_joint <- var_precision_sampler(Y=Y_train_jit, lags=1, k=5,
                                   nrun=7500, burn=2500, thin=5)
out_joint <- fact_prec_processing(res_joint, burnout=750, alpha=0.1,
                                  beta=0.85,
                                  epsilon_grid=seq(0.05,0.1,length=20))
cat("  pFDR:", out_joint$pFDR, "| Edges:", sum(out_joint$P)/2, "\n")

# ── 5. Kastner (VAR + factor SV) ─────────────────────────────────────────────
cat("\n5. Kastner (VAR+FSV)...\n")
sigma_prior <- specify_prior_sigma(
   M                       = p,
   type                    = "factor",
   factor_factors          = 5L,
   factor_priorfacloadtype = "colwiseng",
   factor_heteroskedastic  = TRUE,
   quiet                   = TRUE
)
mod_kast   <- bvar(Y_train_jit, lags=1L, draws=2000, burnin=1000,
                   prior_sigma=sigma_prior, quiet=TRUE)
vc         <- vcov(mod_kast)
Sigma_kast <- apply(vc[1,,,], c(1,2), mean)
Omega_kast <- solve(Sigma_kast)
cat("  Done\n")

# ── 6. Sample covariance baseline ────────────────────────────────────────────
Omega_samp <- solve(cov(Y_train))

# ============================================================================
#  SECTION 5: Predictive evaluation
# ============================================================================

predict_held_out <- function(Omega, y_obs, obs_idx, hid_idx) {
   as.vector(-solve(Omega[hid_idx, hid_idx]) %*%
                Omega[hid_idx, obs_idx] %*% y_obs)
}

set.seed(42)
n_hid <- floor(0.3 * p)
mspe  <- matrix(0, T_test, 7)
colnames(mspe) <- c("PrecFact", "VAR+Glasso", "VAR+PrecFact",
                    "VAPOR", "Kastner", "SampCov", "RandomWalk")

# previous quarter for random walk (starts with last training obs)
y_prev <- Y_last_train

for (t in seq_len(T_test)) {
   y_t     <- Y_test[t, ]
   hid_idx <- sort(sample(p, n_hid))
   obs_idx <- setdiff(seq_len(p), hid_idx)
   y_obs   <- y_t[obs_idx]
   y_true  <- y_t[hid_idx]
   
   # random walk: predict hidden variables as their previous quarter value
   pred_rw <- y_prev[hid_idx]
   
   mspe[t, ] <- c(
      mean((y_true - predict_held_out(out_pf$Omega,    y_obs, obs_idx, hid_idx))^2),
      mean((y_true - predict_held_out(Omega_glasso,    y_obs, obs_idx, hid_idx))^2),
      mean((y_true - predict_held_out(out_pf_var$Omega,y_obs, obs_idx, hid_idx))^2),
      mean((y_true - predict_held_out(out_joint$Omega, y_obs, obs_idx, hid_idx))^2),
      mean((y_true - predict_held_out(Omega_kast,      y_obs, obs_idx, hid_idx))^2),
      mean((y_true - predict_held_out(Omega_samp,      y_obs, obs_idx, hid_idx))^2),
      mean((y_true - pred_rw)^2)
   )
   
   # update previous quarter for random walk
   y_prev <- y_t
}

# ============================================================================
#  SECTION 6: Results
# ============================================================================

cat("\n============================================================\n")
cat("FRED-QD PREDICTIVE EVALUATION RESULTS\n")
cat("Train: Q2 2000 - Q4 2020 (T=77) | Test: Q1 2021 - Q4 2025 (T=20)\n")
cat("Hold out", n_hid, "of", p, "variables per quarter\n")
cat("============================================================\n")

results <- data.frame(
   Method = c("PrecFact (no VAR)", "VAR + Glasso", "VAR + PrecFact",
              "VAPOR (proposed)", "Kastner (VAR+FSV)",
              "Sample Cov", "Random Walk"),
   MSPE = round(colMeans(mspe), 4),
   SD   = round(apply(mspe, 2, sd), 4)
)
results$Result <- sprintf("%.4f (%.4f)", results$MSPE, results$SD)
print(results[, c("Method", "Result")], row.names=FALSE)

cat("\nEdge counts:\n")
cat("PrecFact (no VAR): ", sum(out_pf$P)/2,
    "| pFDR:", out_pf$pFDR, "\n")
cat("VAR + Glasso:      ", sum(gl_var$wi[upper.tri(gl_var$wi)] != 0),
    "| rho:", round(rho_bic, 4), "\n")
cat("VAR + PrecFact:    ", sum(out_pf_var$P)/2,
    "| pFDR:", out_pf_var$pFDR, "\n")
cat("VAPOR:             ", sum(out_joint$P)/2,
    "| pFDR:", out_joint$pFDR, "\n")

# ============================================================================
#  SECTION 7: Heatmap of VAPOR precision graph
# ============================================================================

pc_joint <- -cov2cor(out_joint$Omega * out_joint$P)
diag(pc_joint) <- NA
pc_joint <- round(pc_joint, 2)
rownames(pc_joint) <- colnames(pc_joint) <- var_labels_qd[vars_qd]

pheatmap(pc_joint,
         main            = "FRED-QD Partial Correlation Matrix\n(VAPOR, Q2 2000 - Q4 2020)",
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         color           = colorRampPalette(c("blue","white","red"))(100),
         breaks          = seq(-1, 1, length.out=101),
         na_col          = "grey90",
         display_numbers = TRUE,
         number_color    = "black",
         fontsize_number = 6,
         fontsize        = 8)
