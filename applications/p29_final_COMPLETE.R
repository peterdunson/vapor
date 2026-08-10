# ============================================================================
#  VAPOR p=29 FRED-QD -- full rerun, REDUCED RANK (r = 3, k = 10)
#
#  EDGE RULE: an edge is reported when the central 95% posterior interval of
#  the partial correlation excludes zero. This rule is used EVERYWHERE in this
#  script -- reported edges, sensitivity, stability, the justification sweep,
#  and the rank comparison -- so the paper applies one criterion at both
#  dimensions. The pFDR procedure is not used for selection here.
#
#  Other changes from the full-rank version:
#   (0) Rank/k justification sweep reporting Pr(rho(B)>1), the posterior
#       predictive check, and threshold counts against Kastner.
#   (1) Reduced rank replaces full rank. Full rank puts 55-77% of posterior
#       mass on explosive coefficient matrices and badly under-reproduces the
#       leading eigenvalue of the correlation matrix.
#   (2) KASTNER IS REFIT AT THE MATCHING FACTOR COUNT. Comparing VAPOR at
#       k = 10 against Kastner fixed at 5 factors is not like-for-like; with
#       both at k = 10, VAPOR is sparsest at every threshold. The paper must
#       state that the factor counts are matched.
#   (3) The posterior predictive check uses summaries of y_t, not of residuals
#       formed with a point estimate of B, and reports nine statistics.
#   (4) Held-out MSPE is averaged over several draws of the held-out sets; the
#       single-draw number is also reported for comparability.
#   (5) Diebold-Mariano comparisons are computed here rather than separately.
#
#  Output object keeps the same names as before, so the existing heatmap
#  script works unchanged (only the output filename should be updated).
# ============================================================================
library(VAPOR); library(PrecFact); library(bayesianVARs); library(glasso); library(coda)

# ---- configuration ---------------------------------------------------------
SEED   <- 32
RANK   <- 3
K      <- 10         # VAPOR factor count; Kastner is fit at the same value
R_REPS <- 3
NRUN   <- 7500; BURN <- 2500; THIN <- 5
NRUN_L <- 20000; BURN_L <- 10000
JIT    <- 1e-3
N_HID  <- 20         # held-out set replications
RUN_JUSTIFICATION <- TRUE
OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application"

# ---- data ------------------------------------------------------------------
fred_qd_raw <- read.csv("~/Desktop/2026-05-qd.csv")
transform_series <- function(x, tcode){ x <- as.numeric(x)
switch(as.character(tcode),"1"=x,"2"=c(NA,diff(x)),"3"=c(NA,NA,diff(diff(x))),
       "4"=log(x),"5"=c(NA,diff(log(x))),"6"=c(NA,NA,diff(diff(log(x)))),
       "7"=c(NA,diff(log(x)/lag(log(x))-1)),x) }
tcodes_qd <- as.numeric(fred_qd_raw[2,-1]); names(tcodes_qd) <- colnames(fred_qd_raw)[-1]
fred_qd_data <- fred_qd_raw[-c(1,2),]; dates_qd <- fred_qd_data[,1]
fred_qd_data <- apply(fred_qd_data[,-1],2,as.numeric)
fred_qd_trans <- suppressWarnings(mapply(transform_series, as.data.frame(fred_qd_data), tcodes_qd))
fred_qd_trans <- as.matrix(fred_qd_trans); rownames(fred_qd_trans) <- dates_qd
fred_qd_clean <- fred_qd_trans[apply(fred_qd_trans,1,function(r) !any(is.na(r))),]
fred_qd_clean <- fred_qd_clean[,apply(fred_qd_clean,2,function(c) !any(is.na(c)))]

vars_qd <- c("GDPC1","INDPRO","IPFINAL","IPCONGD","IPBUSEQ","PAYEMS","UEMPMEAN",
             "CLAIMSx","CUMFNS","PCECC96","CMRMTSPLx","BUSINVx","CPIAUCSL","PPICMM","GDPCTPI",
             "OILPRICEx","FEDFUNDS","GS1","GS10","BAA","TB3MS","M2REAL","BUSLOANSx","NONREVSLx",
             "S.P.500","VIXCLSx","GPDIC1","UNRATE","WPSFD49207")
vars_qd <- vars_qd[vars_qd %in% colnames(fred_qd_clean)]
Y_qd <- fred_qd_clean[, vars_qd]; p <- ncol(Y_qd)

parse_qd_date <- function(s){ z <- as.numeric(strsplit(s,"/")[[1]]); z[3]+(z[1]-1)/12 }
date_vals <- vapply(rownames(Y_qd), parse_qd_date, numeric(1))
train_idx <- which(date_vals >= 2000+(6-1)/12 & date_vals <= 2020+(12-1)/12)
test_idx  <- which(date_vals >= 2021+(3-1)/12 & date_vals <= 2025+(12-1)/12)
Y_all <- scale(Y_qd); Y_train <- Y_all[train_idx,]; Y_test <- Y_all[test_idx,]
T_train <- nrow(Y_train); T_test <- nrow(Y_test); Y_last_train <- Y_train[T_train,]
cat(sprintf("p=%d  T_train=%d  T/p=%.2f  T_test=%d | rank=%d k=%d\n",
            p, T_train, T_train/p, T_test, RANK, K))

# ---- shared helpers --------------------------------------------------------
omega_one <- function(fit, s){ Om <- tcrossprod(fit$Lambda[,,s])
diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om }
omega_hat <- function(fit) posterior_mean_omega(fit$Lambda, fit$residual_prec)
pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }
fixB <- function(B){ if (nrow(B)!=p) B <- t(B)
if (ncol(B)!=p) B <- B[,seq_len(p),drop=FALSE]; B }
spec_rad <- function(B) max(abs(eigen(B, only.values=TRUE)$values))
chol_sig <- function(Om){
   Sig <- tryCatch(chol2inv(chol(Om)), error=function(e) MASS::ginv(Om))
   tryCatch(chol(Sig), error=function(e) chol(Sig + diag(1e-8,p))) }

# ---- THE EDGE RULE: 95% credible interval excluding zero -------------------
#  Accepts one fit or a list of fits; draws are pooled across fits.
ci_set <- function(fits){
   if (!is.null(fits$Lambda)) fits <- list(fits)
   stk <- unlist(lapply(fits, function(f)
      lapply(seq_len(dim(f$Lambda)[3]), function(s) pcor_from_omega(omega_one(f,s)))),
      recursive = FALSE)
   PCa <- simplify2array(stk)
   lo <- apply(PCa,c(1,2),quantile,0.025); hi <- apply(PCa,c(1,2),quantile,0.975)
   E <- (lo > 0 | hi < 0); diag(E) <- FALSE
   list(E = E, pm = apply(PCa,c(1,2),mean), n = sum(E)/2) }

ci_table <- function(fits){
   cs <- ci_set(fits)
   ei <- which(cs$E & upper.tri(cs$E), arr.ind = TRUE)
   if (!nrow(ei)) return(data.frame(v1=character(0), v2=character(0), pcor=numeric(0)))
   d <- data.frame(v1 = vars_qd[ei[,1]], v2 = vars_qd[ei[,2]],
                   pcor = round(cs$pm[ei], 3), stringsAsFactors = FALSE)
   d[order(-abs(d$pcor)), ] }

ci_names <- function(fits){
   cs <- ci_set(fits); ei <- which(cs$E & upper.tri(cs$E), arr.ind = TRUE)
   if (!nrow(ei)) return(character(0))
   sort(paste(vars_qd[pmin(ei[,1],ei[,2])], vars_qd[pmax(ei[,1],ei[,2])], sep="--")) }

fit_vapor <- function(kk = K, rk = RANK, sd_i = SEED, nrun = NRUN, burn = BURN,
                      aL = NULL, Y = Y_train){
   f <- NULL; jj <- JIT; a <- 0
   while (is.null(f) && a < 4){ a <- a + 1
   args <- list(Y, lags=1, k=kk, rank=rk, prior_B="dl",
                nrun=nrun, burn=burn, thin=THIN, seed=sd_i)
   if (!is.null(aL)) { args$a_L <- aL; args$b_L <- aL }
   if (a > 1) args$jitter <- jj
   f <- tryCatch(do.call(vapor, args),
                 error=function(e){ message("  retry with jitter"); NULL})
   jj <- jj*4 }
   f }

fit_kastner <- function(kk, Yj){
   sp <- specify_prior_sigma(M=p, type="factor", factor_factors=as.integer(kk),
                             factor_priorfacloadtype="colwiseng",
                             factor_heteroskedastic=TRUE, quiet=TRUE)
   m <- tryCatch(bvar(Yj, lags=1L, draws=2000, burnin=1000, prior_sigma=sp, quiet=TRUE),
                 error=function(e) NULL)
   if (is.null(m)) return(NULL)
   solve(apply(vcov(m)[1,,,], c(1,2), mean)) }

pr_explosive <- function(fit, n = 400){
   nd <- dim(fit$B)[3]; idx <- unique(round(seq(1,nd,length.out=min(n,nd))))
   rs <- vapply(idx, function(s) spec_rad(fixB(fit$B[,,s])), numeric(1))
   c(pr_expl = mean(rs > 1), median_rho = median(rs),
     q05 = quantile(rs,.05,names=FALSE), q95 = quantile(rs,.95,names=FALSE)) }

count_edges <- function(PC, thr) sum(abs(PC[upper.tri(PC)]) > thr, na.rm=TRUE)
THRESH <- c(0.05,0.1,0.2,0.3,0.4,0.5)

# ---- held-out prediction, averaged over draws of the hidden sets -----------
predict_held_out <- function(Omega, y_obs, obs_idx, hid_idx)
   as.vector(-solve(Omega[hid_idx,hid_idx]) %*% Omega[hid_idx,obs_idx] %*% y_obs)

set.seed(SEED); n_hid <- floor(0.3*p)
hid_sets <- lapply(seq_len(N_HID), function(b)
   lapply(seq_len(T_test), function(t) sort(sample(p, n_hid))))
hid_list <- hid_sets[[1]]        # the single draw used in the previous version

loss_one <- function(Omega, hl) vapply(seq_len(T_test), function(t){
   hid <- hl[[t]]; obs <- setdiff(seq_len(p), hid); yt <- Y_test[t,]
   mean((yt[hid] - predict_held_out(Omega, yt[obs], obs, hid))^2)}, numeric(1))
loss_avg <- function(Omega) rowMeans(vapply(hid_sets, function(hl) loss_one(Omega,hl),
                                            numeric(T_test)))
mspe_single <- function(Omega) mean(loss_one(Omega, hid_list))
mspe_avg    <- function(Omega) mean(loss_avg(Omega))

# ---- summary statistics for the posterior predictive check ----------------
summ_stats <- function(M){
   C  <- cor(M); ev <- sort(eigen(C, symmetric=TRUE, only.values=TRUE)$values, decreasing=TRUE)
   off <- C[upper.tri(C)]
   PCm <- tryCatch({ Omh <- solve(cov(M)); pcor_from_omega(Omh) }, error=function(e) NULL)
   pv  <- if (is.null(PCm)) c(NA,NA) else {
      o <- abs(PCm[upper.tri(PCm)]); c(max(o), mean(o > 0.2)) }
   c(max_eig      = ev[1],
     eig2         = ev[2],
     mean_abs_cor = mean(abs(off)),
     q90_abs_cor  = quantile(abs(off), .90, names=FALSE),
     max_abs_pcor = pv[1],
     prop_pcor_02 = pv[2],
     mean_sd      = mean(apply(M, 2, sd)),
     kurtosis     = mean(apply(M, 2, function(z) mean((z-mean(z))^4)/var(z)^2))) }

ac1_of <- function(M) mean(vapply(seq_len(ncol(M)),
                                  function(j) cor(M[-1,j], M[-nrow(M),j]), numeric(1)))

#  Nine summaries. The first eight are computed on y directly. The ninth is a
#  Bayesian residual check: under draw s the residuals y_t - B_s y_{t-1} should
#  be white, and observed and replicated residuals are formed identically.
ppc_data <- function(fit, Y = Y_train, n_draws = 300, seed = 1){
   Tn <- nrow(Y); nd <- dim(fit$Lambda)[3]
   idx <- unique(round(seq(1, nd, length.out=min(n_draws,nd))))
   Ylag <- Y[-Tn,,drop=FALSE]; Ycur <- Y[-1,,drop=FALSE]
   set.seed(seed)
   out <- lapply(idx, function(s){
      B <- fixB(fit$B[,,s]); L <- chol_sig(omega_one(fit,s))
      M <- Ylag %*% t(B)
      Yrep <- M + matrix(rnorm((Tn-1)*p), Tn-1, p) %*% L
      list(rep = c(summ_stats(Yrep), resid_ac1 = ac1_of(Yrep - M)),
           obs = c(summ_stats(Ycur), resid_ac1 = ac1_of(Ycur - M))) })
   RP <- t(vapply(out, `[[`, numeric(9), "rep"))
   OB <- t(vapply(out, `[[`, numeric(9), "obs"))
   res <- data.frame(statistic = colnames(RP),
                     observed  = round(colMeans(OB), 3),
                     rep_mean  = round(colMeans(RP), 3),
                     lo        = round(apply(RP, 2, quantile, .025), 3),
                     hi        = round(apply(RP, 2, quantile, .975), 3),
                     ppp       = round(colMeans(RP >= OB), 3),
                     row.names = NULL)
   attr(res, "replicates") <- RP
   res }

# ============================================================================
#  1. BASELINES  (Kastner at the matching factor count)
# ============================================================================
cat("\n=== 1. baselines ===\n")
set.seed(SEED)
Y_train_jit <- Y_train + matrix(rnorm(prod(dim(Y_train)),0,1e-3), T_train, p)
mod_var <- bvar(Y_train_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
B_hat <- apply(mod_var$PHI,c(1,2),mean); resid_var <- mod_var$Y - mod_var$X %*% B_hat
S_resid <- cov(scale(resid_var)); rho_grid <- 10^seq(-2,0,length=30)
bic_vals <- sapply(rho_grid, function(rr){
   gl <- tryCatch(glasso(S_resid, rho=rr, penalize.diagonal=FALSE), error=function(e) NULL)
   if (is.null(gl)) return(Inf); Om <- gl$wi; ne <- sum(Om[upper.tri(Om)] != 0)
   ll <- (T_train/2)*(log(det(Om)) - sum(diag(S_resid %*% Om))); -2*ll + log(T_train)*ne })
gl_var <- glasso(S_resid, rho=rho_grid[which.min(bic_vals)], penalize.diagonal=FALSE)
Omega_glasso <- gl_var$wi; glasso_edges <- sum(gl_var$wi[upper.tri(gl_var$wi)] != 0)

Omega_kast <- fit_kastner(K, Y_train_jit)
if (is.null(Omega_kast)) stop("Kastner fit failed")
Omega_samp <- solve(cov(Y_train))
cat(sprintf("Kastner fit with %d factors, matching VAPOR\n", K))

loss_glasso <- loss_avg(Omega_glasso); loss_kast <- loss_avg(Omega_kast)
loss_samp   <- loss_avg(Omega_samp)
loss_rw <- rowMeans(vapply(hid_sets, function(hl)
   vapply(seq_len(T_test), function(t){ hid <- hl[[t]]
   y_prev <- if (t==1) Y_last_train else Y_test[t-1,]
   mean((Y_test[t,hid] - y_prev[hid])^2)}, numeric(1)), numeric(T_test)))
cat(sprintf("condition numbers -- Kastner %.3g | Glasso %.3g | sample %.3g\n",
            kappa(Omega_kast, exact=TRUE), kappa(Omega_glasso, exact=TRUE),
            kappa(Omega_samp, exact=TRUE)))

# ============================================================================
#  0. RANK / K JUSTIFICATION  (supplement table)
# ============================================================================
if (RUN_JUSTIFICATION){
   cat("\n=== 0. rank / k justification ===\n")
   K_GRID <- c(5, 8, 10, 15)
   kast_by_k <- setNames(lapply(K_GRID, function(kk){
      cat(sprintf("\r  Kastner at k=%d ", kk)); flush.console()
      fit_kastner(kk, Y_train_jit) }), as.character(K_GRID))
   cat("\n")
   
   just <- NULL
   for (rk in list("full", 2, 3)) for (kk in K_GRID) for (sd_i in c(33,34,35)){
      lab <- if (identical(rk,"full")) "full" else as.character(rk)
      cat(sprintf("\r  rank=%-4s k=%-2d seed=%d   ", lab, kk, sd_i)); flush.console()
      f <- fit_vapor(kk, rk, sd_i); if (is.null(f)) next
      re <- pr_explosive(f); pc <- ppc_data(f, n_draws=150)
      PCv <- pcor_from_omega(omega_hat(f)); diag(PCv) <- NA
      Omk <- kast_by_k[[as.character(kk)]]
      wins <- if (is.null(Omk)) NA else {
         PCk <- pcor_from_omega(Omk); diag(PCk) <- NA
         sum(vapply(THRESH, function(t) count_edges(PCv,t) <= count_edges(PCk,t), logical(1))) }
      just <- rbind(just, data.frame(
         rank = lab, k = kk, seed = sd_i,
         pr_expl = round(re["pr_expl"],3), median_rho = round(re["median_rho"],3),
         n_fail  = sum(pc$ppp < 0.05 | pc$ppp > 0.95),
         ppp_eig = pc$ppp[pc$statistic=="max_eig"],
         rep_eig = pc$rep_mean[pc$statistic=="max_eig"],
         obs_eig = pc$observed[pc$statistic=="max_eig"],
         mspe1   = round(mspe_single(omega_hat(f)),4),
         ci_edges = ci_set(f)$n,
         sparser_than_kastner = wins)) }
   cat("\n")
   just_agg <- aggregate(cbind(pr_expl, median_rho, n_fail, ppp_eig, rep_eig,
                               mspe1, ci_edges, sparser_than_kastner) ~ rank + k,
                         data = just, FUN = mean)
   just_agg <- just_agg[order(just_agg$rank, just_agg$k), ]
   cat(sprintf("\n--- averaged over 3 seeds (observed max eigenvalue = %.2f) ---\n",
               mean(just$obs_eig)))
   print(transform(just_agg, pr_expl=round(pr_expl,3), median_rho=round(median_rho,3),
                   n_fail=round(n_fail,1), ppp_eig=round(ppp_eig,3),
                   rep_eig=round(rep_eig,2), mspe1=round(mspe1,4),
                   ci_edges=round(ci_edges,1),
                   sparser_than_kastner=round(sparser_than_kastner,1)), row.names=FALSE)
   cat("\nRead: full rank places most posterior mass on explosive B and badly\n",
       "under-reproduces the leading eigenvalue. Reduced rank fixes both.\n",
       "sparser_than_kastner counts thresholds (of 6) where VAPOR reports no\n",
       "more nonzero partial correlations than Kastner at the same k.\n",
       "ci_edges is the count under the 95% credible-interval rule, per fit.\n")
} else { just <- just_agg <- kast_by_k <- NULL }

# ============================================================================
#  2. MAIN FITS
# ============================================================================
cat("\n=== 2. main fits (rank =", RANK, ", k =", K, ") ===\n")
pf_mspe <- pfvar_mspe <- v_mspe1 <- v_mspe <- v_edges <- numeric(R_REPS)
vapor_fits <- vector("list", R_REPS); vapor_times <- numeric(R_REPS)
v_loss <- pf_loss <- pfv_loss <- vector("list", R_REPS)

run_pf <- function(Ymat, seed0){
   out <- NULL; a <- 0
   while (is.null(out) && a < 6){ a <- a + 1; set.seed(seed0 + a*7919)
   Yj <- if (a==1) Ymat else Ymat + matrix(rnorm(length(Ymat),0,1e-6), nrow(Ymat))
   out <- tryCatch({ res <- prec_est(as=0.1,bs=0.1,a=0.5,b=0.5,k=K,nrun=NRUN,thin=THIN,
                                     Y=Yj, alldiff=0, a_dir=0.1, b_dir=0.1)
   o <- fact_prec_processing(res, burnout=750, alpha=0.1, beta=0.85,
                             epsilon_grid=seq(0.05,0.1,length=20))
   diag(o$P) <- 0; o }, error=function(e) NULL) }
   if (is.null(out)) stop("prec_est failed"); out }

for (r in seq_len(R_REPS)){
   cat(sprintf("\r  rep %d/%d", r, R_REPS)); flush.console()
   o1 <- run_pf(Y_train, SEED + r*100)
   pf_loss[[r]] <- loss_avg(o1$Omega); pf_mspe[r] <- mean(pf_loss[[r]])
   o2 <- run_pf(scale(resid_var), SEED + r*100 + 1)
   pfv_loss[[r]] <- loss_avg(o2$Omega); pfvar_mspe[r] <- mean(pfv_loss[[r]])
   tt <- system.time( f <- fit_vapor(sd_i = SEED + r) )
   if (is.null(f)) stop(sprintf("VAPOR rep %d failed", r))
   vapor_times[r] <- tt["elapsed"]
   Omv <- omega_hat(f)
   v_loss[[r]] <- loss_avg(Omv)
   v_mspe[r] <- mean(v_loss[[r]]); v_mspe1[r] <- mspe_single(Omv)
   v_edges[r] <- ci_set(f)$n; vapor_fits[[r]] <- f }
cat("\n")

fmt <- function(v) sprintf("%.4f (%.4f)", mean(v), sd(v))
results <- data.frame(
   Method = c("PrecFact (no VAR)","VAR + Glasso","VAR + PrecFact","Kastner (VAR+FSV)",
              "Sample covariance","Random walk","VAPOR (proposed)"),
   MSPE = c(fmt(pf_mspe), sprintf("%.4f", mean(loss_glasso)), fmt(pfvar_mspe),
            sprintf("%.4f", mean(loss_kast)), sprintf("%.4f", mean(loss_samp)),
            sprintf("%.4f", mean(loss_rw)), fmt(v_mspe)))
cat("\n=== p=29 MSPE, averaged over", N_HID, "held-out draws ===\n")
print(results, row.names=FALSE)
cat(sprintf("\nVAPOR on the single held-out draw only: %.4f (sd %.4f)\n",
            mean(v_mspe1), sd(v_mspe1)))
cat(sprintf("VAPOR edges (95%% CI, per fit): %.1f (sd %.1f) | mean fit %.1fs\n",
            mean(v_edges), sd(v_edges), mean(vapor_times)))
re_main <- pr_explosive(vapor_fits[[1]])
cat(sprintf("Pr(rho(B) > 1): %.3f | median rho %.3f [%.2f, %.2f]\n",
            re_main["pr_expl"], re_main["median_rho"], re_main["q05"], re_main["q95"]))

# ---- Diebold-Mariano with the Harvey-Leybourne-Newbold correction ---------
dm_test <- function(loss_comp, loss_vapor){
   d <- loss_comp - loss_vapor; n <- length(d)
   se <- sd(d)/sqrt(n); stat <- (mean(d)/se) * sqrt((n - 1)/n)   # h = 1
   data.frame(Difference = round(mean(d),3), SE = round(se,3), DM = round(stat,2),
              p = round(2*pt(-abs(stat), df=n-1),3),
              Won = sprintf("%d / %d", sum(d > 0), n)) }
lv <- Reduce(`+`, v_loss)/R_REPS
dm <- rbind(
   cbind(Method="VAR + Glasso",      MSPE=round(mean(loss_glasso),3), dm_test(loss_glasso, lv)),
   cbind(Method="Kastner (VAR+FSV)", MSPE=round(mean(loss_kast),3),   dm_test(loss_kast, lv)),
   cbind(Method="Sample covariance", MSPE=round(mean(loss_samp),3),   dm_test(loss_samp, lv)),
   cbind(Method="Random walk",       MSPE=round(mean(loss_rw),3),     dm_test(loss_rw, lv)),
   cbind(Method="VAR + PrecFact",    MSPE=round(mean(Reduce(`+`,pfv_loss)/R_REPS),3),
         dm_test(Reduce(`+`,pfv_loss)/R_REPS, lv)),
   cbind(Method="PrecFact (no VAR)", MSPE=round(mean(Reduce(`+`,pf_loss)/R_REPS),3),
         dm_test(Reduce(`+`,pf_loss)/R_REPS, lv)))
cat("\n=== paired predictive comparison (positive favours VAPOR) ===\n")
print(dm, row.names=FALSE)

# ============================================================================
#  3. MCMC DIAGNOSTICS
# ============================================================================
cat("\n=== 3. convergence (4 chains) ===\n")
chains <- Filter(Negate(is.null),
                 lapply(1:4, function(s) fit_vapor(sd_i = 100 + s, nrun = NRUN_L, burn = BURN_L)))

functionals <- function(fit){
   nd <- dim(fit$Lambda)[3]
   pairs <- list(PAYEMS_UNRATE  = c("PAYEMS","UNRATE"),
                 FEDFUNDS_TB3MS = c("FEDFUNDS","TB3MS"),
                 GDPC1_GPDIC1   = c("GDPC1","GPDIC1"),
                 GS10_BAA       = c("GS10","BAA"),
                 GDPC1_PCECC96  = c("GDPC1","PCECC96"),
                 CPI_PPI        = c("CPIAUCSL","WPSFD49207"),
                 PCECC96_GPDIC1 = c("PCECC96","GPDIC1"))
   pairs <- pairs[vapply(pairs, function(pr) all(pr %in% vars_qd), logical(1))]
   idx <- lapply(pairs, function(pr) c(which(vars_qd==pr[1]), which(vars_qd==pr[2])))
   mat <- matrix(0, nd, length(pairs)); colnames(mat) <- names(pairs)
   rho <- numeric(nd)
   for (s in seq_len(nd)){
      PC <- pcor_from_omega(omega_one(fit,s))
      for (j in seq_along(idx)) mat[s,j] <- PC[idx[[j]][1], idx[[j]][2]]
      rho[s] <- spec_rad(fixB(fit$B[,,s])) }
   data.frame(loglik = colSums(fit$loglik), mat, rho_B = rho) }

chain_fns <- lapply(chains, functionals)
mc_list <- mcmc.list(lapply(chain_fns, function(d) mcmc(as.matrix(d))))
ess_vals <- effectiveSize(mc_list)
rhat_vals <- tryCatch(gelman.diag(mc_list, autoburnin=FALSE, multivariate=FALSE)$psrf[,1],
                      error=function(e) rep(NA, ncol(chain_fns[[1]])))
diag_tab <- data.frame(quantity=names(ess_vals), ESS=round(ess_vals,0), Rhat=round(rhat_vals,3))
print(diag_tab, row.names=FALSE)

edge_rows <- setdiff(seq_len(nrow(diag_tab)), c(1, nrow(diag_tab)))   # drop loglik and rho_B
rh <- diag_tab$Rhat[edge_rows]
cat(sprintf("\nchains: 4 x %d iterations (burn %d, thin %d)\n", NRUN_L, BURN_L, THIN))
cat(sprintf("of %d reported edges: Rhat<1.1 %d | 1.1-1.5 %d | >1.5 %d | median %.2f\n",
            length(rh), sum(rh<1.1), sum(rh>=1.1 & rh<=1.5), sum(rh>1.5), median(rh)))

plot_names <- c("loglik", "rho_B", "FEDFUNDS_TB3MS",
                "PAYEMS_UNRATE", "GDPC1_GPDIC1", "CPI_PPI")
plot_labs  <- c(loglik         = "Log-likelihood",
                rho_B          = expression(rho(B)),
                FEDFUNDS_TB3MS = "FEDFUNDS - TB3MS",
                PAYEMS_UNRATE  = "PAYEMS - UNRATE",
                GDPC1_GPDIC1   = "GDPC1 - GPDIC1",
                CPI_PPI        = "CPIAUCSL - WPSFD49207")
cols <- c("#1b1b1b","#c0392b","#2874a6","#1e8449")

png(file.path(path.expand(OUT),"diagnostics_p29_trace.png"),
    width=1500, height=760, res=140)
par(mfrow=c(2,3), mar=c(3.6,3.8,2.8,0.8), mgp=c(2.3,0.7,0),
    cex.main=1.0, cex.lab=0.95, cex.axis=0.85, bty="l", las=1)
for (nm in plot_names){
   i  <- which(diag_tab$quantity == nm)
   rh <- diag_tab$Rhat[i]; es <- diag_tab$ESS[i]
   yr <- range(sapply(chain_fns, function(d) d[[nm]]))
   ylb <- if (nm == "loglik") "Log-likelihood"
   else if (nm == "rho_B") "Spectral radius"
   else "Partial correlation"
   plot(chain_fns[[1]][[nm]], type="l", col=cols[1], lwd=0.5, ylim=yr,
        xlab="Draw (thinned)", ylab=ylb,
        main=bquote(.(plot_labs[[nm]]) ~ ~ (hat(R) == .(sprintf("%.2f", rh)) *
                                               ", ESS" ~ .(es))))
   for (cc in seq_along(chain_fns)[-1])
      lines(chain_fns[[cc]][[nm]], col=cols[cc], lwd=0.5)
   if (nm == "rho_B") abline(h=1, col="#c0392b", lty=2, lwd=1.2) }
dev.off(); cat("saved diagnostics_p29_trace.png\n")

# ============================================================================
#  4. REPORTED EDGES  (95% credible interval, pooled over all fits)
# ============================================================================
all_fits <- chains
cs_all   <- ci_set(all_fits)
edge_ci  <- cs_all$E
ci_edges <- ci_table(all_fits)
cat(sprintf("\n=== 4. reported edges: 95%% credible interval excludes zero ===\n"))
cat(sprintf("    %d of %d pairs, pooled over %d fits\n",
            cs_all$n, p*(p-1)/2, length(all_fits)))
print(ci_edges, row.names=FALSE)

# ============================================================================
#  5. POSTERIOR PREDICTIVE CHECK
# ============================================================================
cat("\n=== 5. posterior predictive check (summaries of y, nine statistics) ===\n")
ppc29 <- ppc_data(vapor_fits[[1]])
print(ppc29, row.names=FALSE)
cat(sprintf("\nstatistics outside [0.05, 0.95]: %d of %d\n",
            sum(ppc29$ppp < 0.05 | ppc29$ppp > 0.95), nrow(ppc29)))

RPm <- attr(ppc29, "replicates")
png(file.path(path.expand(OUT),"diagnostics_p29_ppc.png"), width=1500, height=900, res=130)
par(mfrow=c(3,3), mar=c(4,3.4,3,1), mgp=c(2.2,0.7,0), bty="l", las=1)
for (j in seq_len(9)){
   obs_j <- ppc29$observed[j]
   hist(RPm[,j], breaks=25, col="grey85", border="grey60",
        xlim=range(c(RPm[,j], obs_j)),
        main=sprintf("%s  (ppp = %.2f)", ppc29$statistic[j], ppc29$ppp[j]),
        xlab=ppc29$statistic[j], ylab="", cex.main=0.95)
   abline(v=obs_j, col="#c0392b", lwd=2) }
dev.off(); cat("saved diagnostics_p29_ppc.png\n")

# ============================================================================
#  6. SENSITIVITY  (edge counts under the credible-interval rule)
# ============================================================================
cat("\n=== 6. sensitivity ===\n")
sens_row <- function(tag, f) data.frame(
   setting = tag, MSPE = round(mspe_avg(omega_hat(f)),4),
   edges = ci_set(f)$n, pr_expl = round(pr_explosive(f)["pr_expl"],3))
sens <- rbind(
   do.call(rbind, lapply(c(5,8,10,15), function(kk) sens_row(sprintf("k=%d",kk), fit_vapor(kk=kk)))),
   do.call(rbind, lapply(c(2,3,4),     function(rr) sens_row(sprintf("rank=%d",rr), fit_vapor(rk=rr)))),
   do.call(rbind, lapply(list(c(5000,2500),c(7500,2500),c(15000,5000)),
                         function(nb) sens_row(sprintf("nrun=%d",nb[1]),
                                               fit_vapor(nrun=nb[1], burn=nb[2])))),
   do.call(rbind, lapply(c(0.3,0.5,1.0), function(aL) sens_row(sprintf("a_L=%.1f",aL),
                                                               fit_vapor(aL=aL)))))
print(sens, row.names=FALSE)
cat(sprintf("MSPE range [%.4f, %.4f], sd %.4f | edges %d to %d\n",
            min(sens$MSPE), max(sens$MSPE), sd(sens$MSPE),
            min(sens$edges), max(sens$edges)))

seed_res <- t(vapply(1:6, function(s){
   f <- fit_vapor(sd_i = s); c(mspe = mspe_avg(omega_hat(f)), edges = ci_set(f)$n) },
   numeric(2)))
cat(sprintf("Seeds: MSPE mean %.4f sd %.4f range [%.4f, %.4f] | edges %.0f to %.0f\n",
            mean(seed_res[,1]), sd(seed_res[,1]), min(seed_res[,1]), max(seed_res[,1]),
            min(seed_res[,2]), max(seed_res[,2])))

stab_sets <- c(
   lapply(c(5,8,10,15),   function(kk) ci_names(fit_vapor(kk = kk))),
   lapply(c(2,3,4),       function(rr) ci_names(fit_vapor(rk = rr))),
   lapply(c(0.3,0.5,1.0), function(aL) ci_names(fit_vapor(aL = aL))))
allE <- sort(unique(unlist(stab_sets)))
memb <- vapply(stab_sets, function(s) allE %in% s, logical(length(allE)))
if (is.null(dim(memb))) memb <- matrix(memb, nrow=length(allE))
stability <- data.frame(edge=allE, n_settings=rowSums(memb), row.names=NULL)
stability <- stability[order(-stability$n_settings), ]
cat("\n--- edge presence across", length(stab_sets), "settings (95% CI rule) ---\n")
print(stability, row.names=FALSE)

# ============================================================================
#  7. AGREEMENT WITH COMPETING ESTIMATORS
# ============================================================================
Om_vapor <- omega_hat(vapor_fits[[1]])
PC29 <- list(VAPOR = cs_all$pm, Kastner = pcor_from_omega(Omega_kast),
             Glasso = pcor_from_omega(Omega_glasso), Sample = pcor_from_omega(Omega_samp))

ed29 <- ci_edges          # the reported set
cat("\n=== 7. rank of each reported edge under each estimator ===\n")
ut <- which(upper.tri(matrix(0,p,p)), arr.ind=TRUE)
key <- function(i,j) paste(pmin(i,j), pmax(i,j), sep="-")
rank_tab <- lapply(PC29, function(M){
   v <- abs(M[cbind(ut[,1],ut[,2])]); o <- order(v, decreasing=TRUE)
   rk <- integer(length(v)); rk[o] <- seq_along(v); setNames(rk, key(ut[,1],ut[,2])) })
ekeys <- key(match(ed29$v1,vars_qd), match(ed29$v2,vars_qd))
edge_ranks <- sapply(rank_tab, function(rk) rk[ekeys])
rownames(edge_ranks) <- paste(ed29$v1, ed29$v2, sep="-")
nE <- nrow(ed29)
cat(sprintf("(of %d pairs)\n", nrow(ut)))
print(cbind(pcor=ed29$pcor, as.data.frame(edge_ranks)))
cat("\nmedian rank:\n"); print(round(apply(edge_ranks,2,median),1))
cat(sprintf("\nin top-%d / top-25 / top-50:\n", nE))
print(sapply(rank_tab, function(rk) c(top_nE=sum(rk[ekeys]<=nE),
                                      top25=sum(rk[ekeys]<=25), top50=sum(rk[ekeys]<=50))))

rates29 <- intersect(c("FEDFUNDS","TB3MS","GS1","GS5","GS10"), vars_qd)
pres <- apply(ed29[,1:2], 1, function(z) paste(sort(z), collapse="|"))
cat("\n--- rate block (reported edges) ---\n")
for (i in seq_len(length(rates29)-1)) for (j in (i+1):length(rates29))
   cat(sprintf("  %-9s %-9s  %s\n", rates29[i], rates29[j],
               ifelse(paste(sort(c(rates29[i],rates29[j])),collapse="|") %in% pres, "EDGE","-")))
cat("\nmarginal correlations among rate series:\n")
print(round(cor(Y_train[, rates29]), 2))

# ============================================================================
#  8. EDGE-COUNT TABLE (Table 2) AND SAVE
# ============================================================================
edge_tab <- data.frame(threshold=THRESH,
                       VAPOR     = sapply(THRESH, function(t) count_edges(PC29$VAPOR,t)),
                       Kastner   = sapply(THRESH, function(t) count_edges(PC29$Kastner,t)),
                       Glasso    = sapply(THRESH, function(t) count_edges(PC29$Glasso,t)),
                       SampleCov = sapply(THRESH, function(t) count_edges(PC29$Sample,t)))
cat(sprintf("\n=== 8. edge counts (|pcor| > threshold), out of %d pairs ===\n", choose(p,2)))
cat(sprintf("    VAPOR and Kastner both fit with %d factors\n", K))
print(edge_tab, row.names=FALSE)
cat(sprintf("\nVAPOR at or below Kastner in %d of %d rows; below Glasso in %d; below sample in %d\n",
            sum(edge_tab$VAPOR <= edge_tab$Kastner), nrow(edge_tab),
            sum(edge_tab$VAPOR <= edge_tab$Glasso),
            sum(edge_tab$VAPOR <= edge_tab$SampleCov)))

saveRDS(list(config = list(rank=RANK, k=K, n_hid_draws=N_HID, kastner_factors=K,
                           edge_rule="95% credible interval excludes zero"),
             results=results, dm=dm, vapor_times=vapor_times, diag_tab=diag_tab,
             edge_ci=edge_ci, ci_edges=ci_edges, reported_edges=ed29,
             edge_ranks=edge_ranks, edge_tab=edge_tab,
             sens=sens, seed_res=seed_res, stability=stability,
             justification=just, justification_agg=just_agg,
             ppc=ppc29, pr_explosive=re_main, pc_vapor=cs_all$pm,
             raw_omega=list(vapor=Om_vapor, glasso=Omega_glasso,
                            kast=Omega_kast, samp=Omega_samp),
             loss=list(vapor=v_loss, glasso=loss_glasso, kast=loss_kast,
                       samp=loss_samp, rw=loss_rw, pf=pf_loss, pfvar=pfv_loss),
             hid_sets=hid_sets, Y_test=Y_test,
             vars=vars_qd, p=p, T_train=T_train),
        file.path(path.expand(OUT), "p29_full_results.rds"))
cat("\nSaved p29_full_results.rds\n")
cat("Heatmap script works unchanged; rename its output to heatmap_p29_r3k10.png\n")


