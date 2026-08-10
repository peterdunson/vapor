# ============================================================================
#  VAPOR p=100 FRED-QD -- full rerun (reduced-rank r=3, k=25) + diagnostics
# ============================================================================
library(VAPOR); library(PrecFact); library(bayesianVARs); library(glasso); library(coda)

SEED <- 2024; R_REPS <- 3; K <- 25; RANK <- 3
NRUN  <- 5000;  BURN  <- 2500;  THIN  <- 5
NRUN_L<- 20000; BURN_L<- 10000; THIN_L<- 10
JIT <- 5e-3
N_HID <- 20
OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application"

# ---- data prep -------------------------------------------------------------
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

vars_qd <- c("GDPC1","PCECC96","PCDGx","PCESVx","PCNDx","GPDIC1","FPIx","PNFIx","GCEC1","EXPGSC1",
             "IMPGSC1","DPIC96","INDPRO","IPFINAL","IPCONGD","IPMAT","IPDMAT","IPNMAT","IPBUSEQ","IPMANSICS",
             "IPFUELS","TCU","CUMFNS","PAYEMS","USPRIV","MANEMP","USGOOD","DMANEMP","USCONS","USFIRE","USINFO",
             "USPBS","USLAH","USMINE","USTPU","USGOVT","USTRADE","CE16OV","CIVPART","UNRATE","UNRATESTx",
             "UNRATELTx","UEMPMEAN","UEMPLT5","UEMP27OV","CLAIMSx","HWIURATIOx","AWHMAN","AWHNONAG","AHETPIx",
             "CES2000000008x","HOUST","HOUST5F","PERMIT","HOUSTS","HOUSTW","USSTHPI","SPCS10RSA","CMRMTSPLx",
             "RSAFSx","AMDMNOx","ANDENOx","BUSINVx","ISRATIOx","PCECTPI","PCEPILFE","CPIAUCSL","CPILFESL",
             "PPIACO","WPSFD49207","WPSID61","OILPRICEx","CPITRNSL","CPIMEDSL","FEDFUNDS","TB3MS","GS1","GS5",
             "GS10","AAA","BAA","MORTGAGE30US","GS10TB3Mx","BAA10YM","M1REAL","M2REAL","BOGMBASEREALx","TOTRESNS",
             "BUSLOANSx","CONSUMERx","REALLNx","NONREVSLx","S.P.500","NASDAQCOM","VIXCLSx","NIKKEI225","EXJPUSx",
             "EXUSUKx","EXCAUSx","UMCSENTx")
vars_qd <- vars_qd[vars_qd %in% colnames(fred_qd_trans)]

pd <- function(s){ z <- as.numeric(strsplit(s,"/")[[1]]); z[3] + (z[1]-1)/12 }
dv <- vapply(rownames(fred_qd_trans), pd, numeric(1)); Ysub <- fred_qd_trans[, vars_qd]
keep <- which(dv >= 1996+(1-1)/12 & dv <= 2025+(12-1)/12); Ysub <- Ysub[keep,,drop=FALSE]
ok <- apply(Ysub,1,function(r) !any(is.na(r))); Ysub <- Ysub[ok,,drop=FALSE]; dvk <- dv[keep][ok]
train_idx <- which(dvk >= 1996+(1-1)/12 & dvk <= 2020+(12-1)/12)
test_idx  <- which(dvk >= 2021+(3-1)/12 & dvk <= 2025+(12-1)/12)
Yall <- scale(Ysub); Y_train <- Yall[train_idx,]; Y_test <- Yall[test_idx,]
p <- ncol(Y_train); T_train <- nrow(Y_train); T_test <- nrow(Y_test)
Y_last_train <- Y_train[T_train,]
cat(sprintf("p=%d  T_train=%d  T/p=%.2f  test=%d\n", p, T_train, T_train/p, T_test))

# ---- helpers ---------------------------------------------------------------
omega_one <- function(fit,s){ Om <- tcrossprod(fit$Lambda[,,s])
diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om }
omega_hat <- function(fit) posterior_mean_omega(fit$Lambda, fit$residual_prec)
pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }
fixB <- function(B){ if (nrow(B)!=p) B <- t(B)
if (ncol(B)!=p) B <- B[,seq_len(p),drop=FALSE]; B }
spec_rad <- function(B) max(abs(eigen(B, only.values=TRUE)$values))
chol_sig <- function(Om){
   Sig <- tryCatch(chol2inv(chol(Om)), error=function(e) MASS::ginv(Om))
   tryCatch(chol(Sig), error=function(e) chol(Sig + diag(1e-8,p))) }

#  the reporting rule: 95% credible interval excluding zero
ci_set <- function(fits){
   if (!is.null(fits$Lambda)) fits <- list(fits)
   stk <- unlist(lapply(fits, function(f)
      lapply(seq_len(dim(f$Lambda)[3]), function(s) pcor_from_omega(omega_one(f,s)))),
      recursive = FALSE)
   PCa <- simplify2array(stk)
   lo <- apply(PCa,c(1,2),quantile,0.025); hi <- apply(PCa,c(1,2),quantile,0.975)
   E <- (lo > 0 | hi < 0); diag(E) <- FALSE
   list(E = E, pm = apply(PCa,c(1,2),mean), n = sum(E)/2) }

pr_explosive <- function(fit, n = 400){
   nd <- dim(fit$B)[3]; idx <- unique(round(seq(1,nd,length.out=min(n,nd))))
   rs <- vapply(idx, function(s) spec_rad(fixB(fit$B[,,s])), numeric(1))
   c(pr_expl = mean(rs > 1), median_rho = median(rs)) }

# ---- MSPE harness, averaged over draws of the held-out sets ----------------
predict_held_out <- function(Om,y_obs,obs,hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)
set.seed(SEED); n_hid <- floor(0.3*p)
hid_sets <- lapply(seq_len(N_HID), function(b)
   lapply(seq_len(T_test), function(t) sort(sample(p, n_hid))))
hid_list <- hid_sets[[1]]

loss_one <- function(Om, hl) vapply(seq_len(T_test), function(t){
   hid <- hl[[t]]; obs <- setdiff(seq_len(p), hid); yt <- Y_test[t,]
   mean((yt[hid] - predict_held_out(Om, yt[obs], obs, hid))^2)}, numeric(1))
loss_avg <- function(Om) rowMeans(vapply(hid_sets, function(hl) loss_one(Om,hl),
                                         numeric(T_test)))
mspe_avg <- function(Om) mean(loss_avg(Om))

# ---- posterior predictive check (seven summaries; T = p here) -------------
summ_stats <- function(M){
   C <- cor(M); ev <- sort(eigen(C, symmetric=TRUE, only.values=TRUE)$values, decreasing=TRUE)
   off <- C[upper.tri(C)]
   c(max_eig      = ev[1],
     eig2         = ev[2],
     mean_abs_cor = mean(abs(off)),
     q90_abs_cor  = quantile(abs(off), .90, names=FALSE),
     mean_sd      = mean(apply(M, 2, sd)),
     kurtosis     = mean(apply(M, 2, function(z) mean((z-mean(z))^4)/var(z)^2))) }
ac1_of <- function(M) mean(vapply(seq_len(ncol(M)),
                                  function(j) cor(M[-1,j], M[-nrow(M),j]), numeric(1)))

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
   nq <- length(out[[1]]$rep)
   RP <- t(vapply(out, `[[`, numeric(nq), "rep"))
   OB <- t(vapply(out, `[[`, numeric(nq), "obs"))
   res <- data.frame(statistic = colnames(RP),
                     observed  = round(colMeans(OB), 3),
                     rep_mean  = round(colMeans(RP), 3),
                     lo        = round(apply(RP, 2, quantile, .025), 3),
                     hi        = round(apply(RP, 2, quantile, .975), 3),
                     ppp       = round(colMeans(RP >= OB), 3),
                     row.names = NULL)
   attr(res, "replicates") <- RP
   res }

# ---- fitters ---------------------------------------------------------------
fit_vapor <- function(kk = K, rk = RANK, sd_i = SEED, nrun = NRUN, burn = BURN,
                      thin = THIN){
   f <- NULL; jj <- JIT; a <- 0
   while (is.null(f) && a < 4){ a <- a + 1
   f <- tryCatch(vapor(Y_train, lags=1, k=kk, rank=rk, prior_B="dl",
                       nrun=nrun, burn=burn, thin=thin, seed=sd_i,
                       jitter=jj, init="random"),
                 error=function(e){ message("  retry with larger jitter"); NULL})
   jj <- jj*2 }
   f }

fit_kastner <- function(sd_i){
   set.seed(sd_i)
   Yj <- Y_train + matrix(rnorm(prod(dim(Y_train)),0,1e-3), T_train, p)
   sp <- specify_prior_sigma(M=p, type="factor", factor_factors=as.integer(K),
                             factor_priorfacloadtype="colwiseng",
                             factor_heteroskedastic=TRUE, quiet=TRUE)
   m <- tryCatch(bvar(Yj, lags=1L, draws=2000, burnin=1000, prior_sigma=sp, quiet=TRUE),
                 error=function(e) NULL)
   if (is.null(m)) return(NULL)
   solve(apply(vcov(m)[1,,,], c(1,2), mean)) }

run_pf <- function(Ymat, seed0){
   out <- NULL; a <- 0
   while (is.null(out) && a < 6){ a <- a + 1; set.seed(seed0 + a*7919)
   Yj <- if (a==1) Ymat else Ymat + matrix(rnorm(length(Ymat),0,1e-6), nrow(Ymat))
   out <- tryCatch({ res <- prec_est(as=0.1,bs=0.1,a=0.5,b=0.5,k=K,nrun=NRUN,thin=THIN,
                                     Y=Yj, alldiff=0, a_dir=0.1, b_dir=0.1)
   o <- fact_prec_processing(res, burnout=250, alpha=0.1, beta=0.85,
                             epsilon_grid=seq(0.05,0.1,length=20))
   diag(o$P) <- 0; o }, error=function(e) NULL) }
   if (is.null(out)) stop("prec_est failed after retries"); out }

# ============================================================================
#  1. DETERMINISTIC BASELINES
# ============================================================================
cat("\n=== 1. baselines ===\n")
set.seed(SEED)
Y_train_jit <- Y_train + matrix(rnorm(prod(dim(Y_train)),0,1e-3), T_train, p)
mod_var <- bvar(Y_train_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
B_hat <- apply(mod_var$PHI,c(1,2),mean); resid_var <- mod_var$Y - mod_var$X %*% B_hat
S_resid <- cov(scale(resid_var)); rho_grid <- 10^seq(-2,0,length=30)
bic_vals <- sapply(rho_grid, function(rho){
   gl <- tryCatch(glasso(S_resid,rho=rho,penalize.diagonal=FALSE), error=function(e) NULL)
   if (is.null(gl)) return(Inf); Om <- gl$wi; ne <- sum(Om[upper.tri(Om)] != 0)
   ll <- (T_train/2)*(log(det(Om)) - sum(diag(S_resid %*% Om))); -2*ll + log(T_train)*ne })
gl_fit <- glasso(S_resid, rho=rho_grid[which.min(bic_vals)], penalize.diagonal=FALSE)
Omega_glasso <- gl_fit$wi
Omega_samp <- tryCatch(solve(cov(Y_train)), error=function(e) MASS::ginv(cov(Y_train)))

loss_glasso <- loss_avg(Omega_glasso); loss_samp <- loss_avg(Omega_samp)
loss_rw <- rowMeans(vapply(hid_sets, function(hl)
   vapply(seq_len(T_test), function(t){ hid <- hl[[t]]
   y_prev <- if (t==1) Y_last_train else Y_test[t-1,]
   mean((Y_test[t,hid] - y_prev[hid])^2)}, numeric(1)), numeric(T_test)))

# ---- Kastner, averaged over repeated fits ---------------------------------
cat("Kastner: fitting", R_REPS, "times\n")
kast_fits <- Filter(Negate(is.null), lapply(seq_len(R_REPS), function(r) fit_kastner(SEED + r)))
if (!length(kast_fits)) stop("all Kastner fits failed")
kast_loss <- lapply(kast_fits, loss_avg)
kast_mspe <- vapply(kast_loss, mean, numeric(1))
Omega_kast <- kast_fits[[1]]
loss_kast  <- Reduce(`+`, kast_loss)/length(kast_loss)
cat(sprintf("Kastner MSPE across %d fits: %s\n", length(kast_fits),
            paste(sprintf("%.3f", kast_mspe), collapse=", ")))
cat(sprintf("condition numbers -- Kastner %.3g | Glasso %.3g | sample %.3g\n",
            kappa(Omega_kast, exact=TRUE), kappa(Omega_glasso, exact=TRUE),
            kappa(Omega_samp, exact=TRUE)))

# ============================================================================
#  2. MAIN FITS
# ============================================================================
cat("\n=== 2. main fits (rank =", RANK, ", k =", K, ") ===\n")
pf_mspe <- pfvar_mspe <- v_mspe <- v_edges <- numeric(R_REPS)
vapor_fits <- vector("list", R_REPS); vapor_times <- numeric(R_REPS)
v_loss <- pf_loss <- pfv_loss <- vector("list", R_REPS)

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
   v_loss[[r]] <- loss_avg(Omv); v_mspe[r] <- mean(v_loss[[r]])
   v_edges[r] <- ci_set(f)$n; vapor_fits[[r]] <- f }
cat("\n")

fmt <- function(v) sprintf("%.4f (%.4f)", mean(v), sd(v))
results <- data.frame(
   Method = c("PrecFact (no VAR)","VAR + Glasso","VAR + PrecFact","Kastner (VAR+FSV)",
              "Sample covariance","Random walk","VAPOR (proposed)"),
   MSPE = c(fmt(pf_mspe), sprintf("%.4f", mean(loss_glasso)), fmt(pfvar_mspe),
            fmt(kast_mspe), sprintf("%.4f", mean(loss_samp)),
            sprintf("%.4f", mean(loss_rw)), fmt(v_mspe)))
cat("\n=== p=100 MSPE, averaged over", N_HID, "held-out draws ===\n")
print(results, row.names=FALSE)
cat(sprintf("VAPOR edges (95%% CI, per fit): %.1f (sd %.1f) | mean fit %.1fs\n",
            mean(v_edges), sd(v_edges), mean(vapor_times)))
re_main <- pr_explosive(vapor_fits[[1]])
cat(sprintf("Pr(rho(B) > 1): %.3f | median rho %.3f\n",
            re_main["pr_expl"], re_main["median_rho"]))

# ---- Diebold-Mariano -------------------------------------------------------
dm_test <- function(loss_comp, loss_vapor){
   d <- loss_comp - loss_vapor; n <- length(d)
   se <- sd(d)/sqrt(n); stat <- (mean(d)/se) * sqrt((n - 1)/n)
   data.frame(Difference = round(mean(d),3), SE = round(se,3), DM = round(stat,2),
              p = round(2*pt(-abs(stat), df=n-1),3),
              Won = sprintf("%d / %d", sum(d > 0), n)) }
lv <- Reduce(`+`, v_loss)/R_REPS
dm <- rbind(
   cbind(Method="VAR + Glasso",      MSPE=round(mean(loss_glasso),3), dm_test(loss_glasso, lv)),
   cbind(Method="Kastner (VAR+FSV)", MSPE=round(mean(loss_kast),3),   dm_test(loss_kast, lv)),
   cbind(Method="Sample covariance", MSPE=round(mean(loss_samp),3),   dm_test(loss_samp, lv)),
   cbind(Method="VAR + PrecFact",    MSPE=round(mean(Reduce(`+`,pfv_loss)/R_REPS),3),
         dm_test(Reduce(`+`,pfv_loss)/R_REPS, lv)),
   cbind(Method="PrecFact (no VAR)", MSPE=round(mean(Reduce(`+`,pf_loss)/R_REPS),3),
         dm_test(Reduce(`+`,pf_loss)/R_REPS, lv)),
   cbind(Method="Random walk",       MSPE=round(mean(loss_rw),3),     dm_test(loss_rw, lv)))
cat("\n=== paired predictive comparison (positive favours VAPOR) ===\n")
print(dm, row.names=FALSE)

# ============================================================================
#  3. CONVERGENCE CHAINS
# ============================================================================
cat("\n=== 3. convergence (4 chains of", NRUN_L, ") ===\n")
chains_long <- Filter(Negate(is.null),
                      lapply(1:4, function(s){ cat(sprintf("\r  chain %d/4", s)); flush.console()
                         fit_vapor(sd_i = 100 + s, nrun = NRUN_L, burn = BURN_L, thin = THIN_L) }))
cat("\n")
saveRDS(chains_long, file.path(path.expand(OUT), "p100_chains_final.rds"))

# ============================================================================
#  4. REPORTED EDGES
# ============================================================================
all_fits <- chains_long
cs_all <- ci_set(all_fits); edge_ci <- cs_all$E
idx <- which(edge_ci & upper.tri(edge_ci), arr.ind=TRUE)
ci_edges <- data.frame(v1=vars_qd[idx[,1]], v2=vars_qd[idx[,2]],
                       pcor=round(cs_all$pm[idx],3))
ci_edges <- ci_edges[order(-abs(ci_edges$pcor)), ]
cat(sprintf("\n=== 4. reported edges: %d of %d pairs, pooled over %d fits ===\n",
            cs_all$n, p*(p-1)/2, length(all_fits)))
print(ci_edges, row.names=FALSE)

# ---- per-edge diagnostics --------------------------------------------------
pc_arrays <- lapply(chains_long, function(f){
   nd <- dim(f$Lambda)[3]
   simplify2array(lapply(seq_len(nd), function(s) pcor_from_omega(omega_one(f,s)))) })
diag_of <- function(x){ mc <- mcmc.list(lapply(x, mcmc))
c(ESS = sum(effectiveSize(mc)),
  Rhat = tryCatch(gelman.diag(mc, autoburnin=FALSE)$psrf[1], error=function(e) NA)) }

rho_series <- lapply(chains_long, function(f)
   vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s])), numeric(1)))
series <- list(loglik = lapply(chains_long, function(f) colSums(f$loglik)))
for (e in seq_len(nrow(ci_edges))){
   i <- match(ci_edges$v1[e], vars_qd); j <- match(ci_edges$v2[e], vars_qd)
   series[[paste(ci_edges$v1[e], ci_edges$v2[e], sep="_")]] <-
      lapply(pc_arrays, function(A) A[i,j,]) }
series[["rho_B"]] <- rho_series

diag_all <- do.call(rbind, lapply(series, diag_of))
diag_tab <- data.frame(quantity=rownames(diag_all), ESS=round(diag_all[,"ESS"],0),
                       Rhat=round(diag_all[,"Rhat"],3), row.names=NULL)
cat("\n=== per-quantity diagnostics (20k, 4 chains) ===\n"); print(diag_tab, row.names=FALSE)
er <- diag_tab$Rhat[-c(1, nrow(diag_tab))]
cat(sprintf("\nof %d reported edges: Rhat<1.1 %d | 1.1-1.5 %d | >1.5 %d | median %.2f\n",
            length(er), sum(er<1.1,na.rm=TRUE), sum(er>=1.1 & er<=1.5,na.rm=TRUE),
            sum(er>1.5,na.rm=TRUE), median(er,na.rm=TRUE)))
edge_rhat <- data.frame(v1=ci_edges$v1, v2=ci_edges$v2, pcor=ci_edges$pcor,
                        ESS=diag_tab$ESS[2:(nrow(ci_edges)+1)],
                        Rhat=diag_tab$Rhat[2:(nrow(ci_edges)+1)])

# ---- trace figure: six panels spanning the range ---------------------------
qn <- diag_tab$quantity
worst <- qn[-c(1,length(qn))][order(-diag_tab$Rhat[-c(1,nrow(diag_tab))])][1:3]
plot_names <- c("loglik", "rho_B", worst,
                qn[-c(1,length(qn))][order(diag_tab$Rhat[-c(1,nrow(diag_tab))])][1])
cols <- c("#1b1b1b","#c0392b","#2874a6","#1e8449")
png(file.path(path.expand(OUT),"diagnostics_p100_trace.png"),
    width=1500, height=760, res=140)
par(mfrow=c(2,3), mar=c(3.6,3.8,2.8,0.8), mgp=c(2.3,0.7,0),
    cex.main=1.0, cex.lab=0.95, cex.axis=0.85, bty="l", las=1)
for (nm in plot_names){
   i <- which(diag_tab$quantity == nm)
   rh <- diag_tab$Rhat[i]; es <- diag_tab$ESS[i]
   x <- series[[nm]]; yr <- range(unlist(x))
   ylb <- if (nm=="loglik") "Log-likelihood" else
      if (nm=="rho_B") "Spectral radius" else "Partial correlation"
   plot(x[[1]], type="l", col=cols[1], lwd=0.5, ylim=yr,
        xlab="Draw (thinned)", ylab=ylb,
        main=bquote(.(gsub("_"," - ",nm)) ~ ~ (hat(R) == .(sprintf("%.2f",rh)) *
                                                  ", ESS" ~ .(es))))
   for (cc in 2:length(x)) lines(x[[cc]], col=cols[cc], lwd=0.5)
   if (nm=="rho_B") abline(h=1, col="#c0392b", lty=2, lwd=1.2) }
dev.off(); cat("Saved diagnostics_p100_trace.png\n")

# ============================================================================
#  5. POSTERIOR PREDICTIVE CHECK
# ============================================================================
cat("\n=== 5. posterior predictive check ===\n")
ppc100 <- ppc_data(vapor_fits[[1]])
print(ppc100, row.names=FALSE)
cat(sprintf("\nstatistics outside [0.05, 0.95]: %d of %d\n",
            sum(ppc100$ppp < 0.05 | ppc100$ppp > 0.95), nrow(ppc100)))
RPm <- attr(ppc100, "replicates"); nq <- nrow(ppc100)
png(file.path(path.expand(OUT),"diagnostics_p100_ppc.png"),
    width=1500, height=700, res=130)
par(mfrow=c(ceiling(nq/4),4), mar=c(4,3.4,3,1), mgp=c(2.2,0.7,0), bty="l", las=1)
for (j in seq_len(nq)){
   obs_j <- ppc100$observed[j]
   hist(RPm[,j], breaks=25, col="grey85", border="grey60",
        xlim=range(c(RPm[,j], obs_j)),
        main=sprintf("%s  (ppp = %.2f)", ppc100$statistic[j], ppc100$ppp[j]),
        xlab=ppc100$statistic[j], ylab="", cex.main=0.95)
   abline(v=obs_j, col="#c0392b", lwd=2) }
dev.off(); cat("Saved diagnostics_p100_ppc.png\n")

# ============================================================================
#  6. SENSITIVITY
# ============================================================================
cat("\n=== 6. sensitivity ===\n")
sens_row <- function(tag, f) data.frame(setting=tag,
                                        MSPE=round(mspe_avg(omega_hat(f)),4), edges=ci_set(f)$n,
                                        pr_expl=round(pr_explosive(f)["pr_expl"],3))
fit_vapor_aL <- function(aL){
   f <- NULL; jj <- JIT; a <- 0
   while (is.null(f) && a < 4){ a <- a + 1
   f <- tryCatch(vapor(Y_train, lags=1, k=K, rank=RANK, prior_B="dl",
                       a_L=aL, b_L=aL, nrun=NRUN, burn=BURN, thin=THIN,
                       seed=SEED, jitter=jj, init="random"),
                 error=function(e) NULL); jj <- jj*2 }
   f }
sens <- rbind(
   do.call(rbind, lapply(c(15,25,35), function(kk) sens_row(sprintf("k=%d",kk), fit_vapor(kk=kk)))),
   do.call(rbind, lapply(c(2,3,4),    function(rr) sens_row(sprintf("rank=%d",rr), fit_vapor(rk=rr)))),
   do.call(rbind, lapply(list(c(2500,1250),c(5000,2500),c(10000,5000)),
                         function(nb) sens_row(sprintf("nrun=%d",nb[1]),
                                               fit_vapor(nrun=nb[1], burn=nb[2])))),
   do.call(rbind, lapply(c(0.3,0.5,1.0), function(aL) sens_row(sprintf("a_L=%.1f",aL),
                                                               fit_vapor_aL(aL)))))
print(sens, row.names=FALSE)
cat(sprintf("MSPE range [%.4f, %.4f] sd %.4f | edges %d to %d\n",
            min(sens$MSPE), max(sens$MSPE), sd(sens$MSPE),
            min(sens$edges), max(sens$edges)))
seed_res <- t(vapply(1:6, function(s){
   f <- fit_vapor(sd_i=s); c(mspe=mspe_avg(omega_hat(f)), edges=ci_set(f)$n) }, numeric(2)))
cat(sprintf("Seeds: MSPE mean %.4f sd %.4f range [%.4f, %.4f] | edges %.0f to %.0f\n",
            mean(seed_res[,1]), sd(seed_res[,1]), min(seed_res[,1]), max(seed_res[,1]),
            min(seed_res[,2]), max(seed_res[,2])))

# ============================================================================
#  7. AGREEMENT WITH COMPETING ESTIMATORS
# ============================================================================
Om_vapor <- omega_hat(vapor_fits[[1]])
PC100 <- list(VAPOR = cs_all$pm, Kastner = pcor_from_omega(Omega_kast),
              Glasso = pcor_from_omega(Omega_glasso), Sample = pcor_from_omega(Omega_samp))
ut  <- which(upper.tri(matrix(0,p,p)), arr.ind=TRUE)
key <- function(i,j) paste(pmin(i,j), pmax(i,j), sep="-")
rank_tab <- lapply(PC100, function(M){
   v <- abs(M[cbind(ut[,1],ut[,2])]); o <- order(v, decreasing=TRUE)
   rk <- integer(length(v)); rk[o] <- seq_along(v); setNames(rk, key(ut[,1],ut[,2])) })
ekeys <- key(match(ci_edges$v1,vars_qd), match(ci_edges$v2,vars_qd))
edge_ranks <- sapply(rank_tab, function(rk) rk[ekeys])
rownames(edge_ranks) <- paste(ci_edges$v1, ci_edges$v2, sep="-")
nE <- nrow(ci_edges)
cat(sprintf("\n=== 7. rank of each reported edge (of %d pairs) ===\n", nrow(ut)))
print(cbind(pcor=ci_edges$pcor, as.data.frame(edge_ranks)))
cat("\nmedian rank:\n"); print(round(apply(edge_ranks,2,median),1))
cat(sprintf("\nin top-%d / top-100 / top-250:\n", nE))
print(sapply(rank_tab, function(rk) c(top_nE=sum(rk[ekeys]<=nE),
                                      top100=sum(rk[ekeys]<=100), top250=sum(rk[ekeys]<=250))))

rates <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
pres <- apply(ci_edges[,1:2], 1, function(z) paste(sort(z), collapse="|"))
cat("\n=== term-structure adjacency ===\n")
for (i in 1:(length(rates)-1)) for (j in (i+1):length(rates))
   cat(sprintf("  %-9s %-9s gap %d  %s\n", rates[i], rates[j], j-i,
               ifelse(paste(sort(c(rates[i],rates[j])),collapse="|") %in% pres, "EDGE","-")))

# ============================================================================
#  8. EDGE-COUNT TABLE AND SAVE
# ============================================================================
THRESH <- c(0.05,0.1,0.2,0.3,0.4,0.5)
count_edges <- function(PC,thr) sum(abs(PC[upper.tri(PC)]) > thr, na.rm=TRUE)
edge_tab <- data.frame(threshold=THRESH,
                       VAPOR     = sapply(THRESH, function(t) count_edges(PC100$VAPOR,t)),
                       Kastner   = sapply(THRESH, function(t) count_edges(PC100$Kastner,t)),
                       Glasso    = sapply(THRESH, function(t) count_edges(PC100$Glasso,t)),
                       SampleCov = sapply(THRESH, function(t) count_edges(PC100$Sample,t)))
cat(sprintf("\n=== 8. edge counts, out of %d pairs (VAPOR and Kastner both at k=%d) ===\n",
            choose(p,2), K))
print(edge_tab, row.names=FALSE)

saveRDS(list(config=list(rank=RANK,k=K,n_hid_draws=N_HID,kastner_fits=length(kast_fits),
                         edge_rule="95% credible interval excludes zero"),
             results=results, dm=dm, vapor_times=vapor_times, diag_tab=diag_tab,
             edge_ci=edge_ci, ci_edges=ci_edges, edge_rhat=edge_rhat,
             edge_ranks=edge_ranks, edge_tab=edge_tab,
             sens=sens, seed_res=seed_res, ppc=ppc100, pr_explosive=re_main,
             pc_vapor=cs_all$pm,
             kast_mspe=kast_mspe,
             raw_omega=list(vapor=Om_vapor, glasso=Omega_glasso,
                            kast=Omega_kast, samp=Omega_samp),
             loss=list(vapor=v_loss, glasso=loss_glasso, kast=loss_kast,
                       samp=loss_samp, rw=loss_rw, pf=pf_loss, pfvar=pfv_loss),
             hid_sets=hid_sets, Y_test=Y_test,
             vars=vars_qd, p=p, T_train=T_train),
        file.path(path.expand(OUT), "p100_full_results.rds"))
cat("\nSaved p100_full_results.rds\n")
cat("Rerun the heatmap and term-structure figure scripts against the new object.\n")