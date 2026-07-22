# ============================================================================
#  VAPOR p=100 FRED-QD — full rerun (reduced-rank r=3, k=25) + diagnostics A-E
#  Headline/timing at 5k; convergence chains at 20k.
# ============================================================================
library(VAPOR); library(PrecFact); library(bayesianVARs); library(glasso); library(coda)

SEED <- 2024; R_REPS <- 3; K <- 25; RANK <- 3
NRUN  <- 5000;  BURN  <- 2500;  THIN  <- 5      # headline + sweeps
NRUN_L<- 20000; BURN_L<- 10000; THIN_L<- 10     # convergence chains
JIT <- 5e-3
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
fred_qd_trans <- mapply(transform_series, as.data.frame(fred_qd_data), tcodes_qd)
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

# ---- MSPE harness ----------------------------------------------------------
predict_held_out <- function(Om,y_obs,obs,hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)
set.seed(SEED); n_hid <- floor(0.3*p)
hid_list <- lapply(seq_len(T_test), function(t) sort(sample(p, n_hid)))
mspe_for_omega <- function(Om) sapply(seq_len(T_test), function(t){
   hid <- hid_list[[t]]; obs <- setdiff(seq_len(p), hid); yt <- Y_test[t,]
   mean((yt[hid] - predict_held_out(Om, yt[obs], obs, hid))^2)})
select_graph_vapor <- function(fit, alpha=0.1, beta=0.85,
                               epsilon_grid=seq(0.05,0.1,length.out=20), burnout=0){
   res_like <- list(Lambda=fit$Lambda, residuals=fit$residual_prec)
   out <- fact_prec_processing(res_like, burnout=burnout, alpha=alpha, beta=beta,
                               epsilon_grid=epsilon_grid)
   diag(out$P) <- 0; out }
omega_one <- function(fit,s){ Om <- tcrossprod(fit$Lambda[,,s])
diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om }
pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }

# ---- deterministic baselines ----------------------------------------------
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
Omega_glasso <- gl_fit$wi; glasso_edges <- sum(gl_fit$wi[upper.tri(gl_fit$wi)] != 0)
sigma_prior <- specify_prior_sigma(M=p, type="factor", factor_factors=K,
                                   factor_priorfacloadtype="colwiseng",
                                   factor_heteroskedastic=TRUE, quiet=TRUE)
mod_kast <- tryCatch(bvar(Y_train_jit, lags=1L, draws=2000, burnin=1000,
                          prior_sigma=sigma_prior, quiet=TRUE), error=function(e) NULL)
Omega_kast <- if (!is.null(mod_kast)) { vc <- vcov(mod_kast); solve(apply(vc[1,,,],c(1,2),mean)) } else matrix(NA,p,p)
Omega_samp <- tryCatch(solve(cov(Y_train)), error=function(e) MASS::ginv(cov(Y_train)))
mspe_glasso <- mspe_for_omega(Omega_glasso)
mspe_kast   <- if (!any(is.na(Omega_kast))) mspe_for_omega(Omega_kast) else rep(NA,T_test)
mspe_samp   <- mspe_for_omega(Omega_samp)
mspe_rw <- sapply(seq_len(T_test), function(t){ hid <- hid_list[[t]]
y_prev <- if (t==1) Y_last_train else Y_test[t-1,]
mean((Y_test[t,hid] - y_prev[hid])^2)})

# ---- PrecFact retry wrapper (chol failures are intermittent) ---------------
run_pf <- function(Ymat, seed0){
   out <- NULL; attempt <- 0
   while (is.null(out) && attempt < 6){
      attempt <- attempt + 1
      set.seed(seed0 + attempt*7919)
      Yj <- if (attempt == 1) Ymat else Ymat + matrix(rnorm(length(Ymat),0,1e-6), nrow(Ymat))
      out <- tryCatch({
         res <- prec_est(as=0.1,bs=0.1,a=0.5,b=0.5,k=K,nrun=NRUN,thin=THIN,
                         Y=Yj, alldiff=0, a_dir=0.1, b_dir=0.1)
         op <- fact_prec_processing(res, burnout=250, alpha=0.1, beta=0.85,
                                    epsilon_grid=seq(0.05,0.1,length=20))
         diag(op$P) <- 0; op
      }, error=function(e){
         message(sprintf("  prec_est attempt %d failed (%s), retrying",
                         attempt, conditionMessage(e))); NULL})
   }
   if (is.null(out)) stop("prec_est failed after retries")
   out
}

# ---- SECTION 1: headline MSPE ---------------------------------------------
pf_mspe <- pf_edges <- pfvar_mspe <- pfvar_edges <- v_mspe <- v_edges <- numeric(R_REPS)
vapor_fits <- vector("list", R_REPS); vapor_times <- numeric(R_REPS)
v_loss <- vector("list", R_REPS)
for (r in seq_len(R_REPS)){
   cat(sprintf("\rHeadline rep %d/%d", r, R_REPS)); flush.console()
   
   out_pf <- run_pf(Y_train, SEED + r*100)
   pf_mspe[r] <- mean(mspe_for_omega(out_pf$Omega)); pf_edges[r] <- sum(out_pf$P)/2
   
   out_pfv <- run_pf(scale(resid_var), SEED + r*100 + 1)
   pfvar_mspe[r] <- mean(mspe_for_omega(out_pfv$Omega)); pfvar_edges[r] <- sum(out_pfv$P)/2
   
   fit <- NULL; jj <- JIT; attempt <- 0
   while (is.null(fit) && attempt < 4){
      attempt <- attempt + 1
      tt <- system.time( fit <- tryCatch(
         vapor(Y_train, lags=1, k=K, rank=RANK, prior_B="dl",
               nrun=NRUN, burn=BURN, thin=THIN, seed=SEED+r, jitter=jj, init="random"),
         error=function(e){ message(sprintf("\n rep %d jitter %.0e failed, retrying", r, jj)); NULL}) )
      jj <- jj*2
   }
   if (is.null(fit)) stop(sprintf("rep %d failed after retries", r))
   vapor_times[r] <- tt["elapsed"]
   sel <- select_graph_vapor(fit)
   v_loss[[r]] <- mspe_for_omega(sel$Omega)
   v_mspe[r] <- mean(v_loss[[r]]); v_edges[r] <- sum(sel$P)/2
   vapor_fits[[r]] <- fit
}
cat("\n")
fmt <- function(v) sprintf("%.4f (%.4f)", mean(v), sd(v))
results <- data.frame(
   Method = c("PrecFact (no VAR)","VAR + Glasso","VAR + PrecFact","Kastner (VAR+FSV)",
              "Sample covariance","Random walk","VAPOR (proposed)"),
   MSPE = c(fmt(pf_mspe), sprintf("%.4f",mean(mspe_glasso)), fmt(pfvar_mspe),
            sprintf("%.4f",mean(mspe_kast)), sprintf("%.4f",mean(mspe_samp)),
            sprintf("%.4f",mean(mspe_rw)), fmt(v_mspe)))
cat("\n=== p=100 MSPE (VAPOR at 5k) ===\n"); print(results, row.names=FALSE)
cat(sprintf("VAPOR edges (pFDR): %.1f (sd %.1f) | mean fit %.1fs\n",
            mean(v_edges), sd(v_edges), mean(vapor_times)))

# ---- A: 20k convergence chains --------------------------------------------
run_chain <- function(s){
   jj <- JIT; fit <- NULL; a <- 0
   while (is.null(fit) && a < 5){ a <- a + 1
   fit <- tryCatch(vapor(Y_train, lags=1, k=K, rank=RANK, prior_B="dl",
                         nrun=NRUN_L, burn=BURN_L, thin=THIN_L,
                         seed=100+s, jitter=jj, init="random"),
                   error=function(e){ message(sprintf("\n chain %d retry", s)); NULL}); jj <- jj*2 }
   if (is.null(fit)) stop(sprintf("chain %d failed", s)); fit
}
chains_long <- lapply(1:4, function(s){
   cat(sprintf("\r 20k chain %d/4", s)); flush.console(); run_chain(s)})
cat("\n")
saveRDS(chains_long, file.path(path.expand(OUT), "p100_chains_final.rds"))

# ---- B: 95% credible-interval edges (pool 5k fits + 20k chains) ------------
all_fits <- c(vapor_fits, chains_long); pc_stack <- list()
for (fit in all_fits){ nd <- dim(fit$Lambda)[3]
for (s in seq_len(nd)) pc_stack[[length(pc_stack)+1]] <- pcor_from_omega(omega_one(fit,s)) }
PCarr <- simplify2array(pc_stack)
lo <- apply(PCarr,c(1,2),quantile,0.025); hi <- apply(PCarr,c(1,2),quantile,0.975)
edge_ci <- (lo > 0 | hi < 0); diag(edge_ci) <- FALSE
cat(sprintf("\n=== B: 95%% interval edges ===\nedges: %d (of %d)\n",
            sum(edge_ci)/2, p*(p-1)/2))
idx <- which(edge_ci & upper.tri(edge_ci), arr.ind=TRUE)
pm <- apply(PCarr,c(1,2),mean)
ci_edges <- data.frame(v1=vars_qd[idx[,1]], v2=vars_qd[idx[,2]], pcor=round(pm[idx],3))
ci_edges <- ci_edges[order(-abs(ci_edges$pcor)), ]
print(ci_edges, row.names=FALSE)

# ---- A (cont.): per-edge and panel diagnostics -----------------------------
pc_arrays <- lapply(chains_long, function(f){
   nd <- dim(f$Lambda)[3]
   simplify2array(lapply(seq_len(nd), function(s) pcor_from_omega(omega_one(f,s)))) })

diag_of <- function(x){ mc <- mcmc.list(lapply(x, mcmc))
c(ESS = sum(effectiveSize(mc)),
  Rhat = tryCatch(gelman.diag(mc, autoburnin=FALSE)$psrf[1], error=function(e) NA)) }

vix <- function(nm){ i <- match(nm, vars_qd)
if (is.na(i)) stop(sprintf("'%s' missing", nm)); i }
pair_defs <- list(GDPC1_GPDIC1      = c("GDPC1","GPDIC1"),
                  FEDFUNDS_TB3MS    = c("FEDFUNDS","TB3MS"),
                  AAA_BAA           = c("AAA","BAA"),
                  CPIAUCSL_CPITRNSL = c("CPIAUCSL","CPITRNSL"),
                  GS5_GS10          = c("GS5","GS10"))
series <- list(loglik = lapply(chains_long, function(f) colSums(f$loglik)))
for (nm in names(pair_defs)){
   i <- vix(pair_defs[[nm]][1]); j <- vix(pair_defs[[nm]][2])
   series[[nm]] <- lapply(pc_arrays, function(A) A[i,j,]) }

diag_tab <- do.call(rbind, lapply(series, diag_of))
diag_tab <- data.frame(quantity=rownames(diag_tab), ESS=round(diag_tab[,"ESS"],0),
                       Rhat=round(diag_tab[,"Rhat"],3), row.names=NULL)
cat("\n=== A: MCMC diagnostics (20k, 4 chains) ===\n"); print(diag_tab, row.names=FALSE)

# per-edge R-hat across all reported CI edges
allE <- t(apply(cbind(match(ci_edges$v1,vars_qd), match(ci_edges$v2,vars_qd)), 1,
                function(z) diag_of(lapply(pc_arrays, function(A) A[z[1],z[2],]))))
edge_rhat <- data.frame(v1=ci_edges$v1, v2=ci_edges$v2, pcor=ci_edges$pcor,
                        ESS=round(allE[,"ESS"],0), Rhat=round(allE[,"Rhat"],3))
cat("\n=== per-edge diagnostics (ordered by magnitude) ===\n")
print(edge_rhat, row.names=FALSE)
cat(sprintf("\nof %d reported edges: Rhat<1.1 %d | 1.1-1.5 %d | >1.5 %d | median %.2f\n",
            nrow(edge_rhat), sum(allE[,"Rhat"]<1.1,na.rm=TRUE),
            sum(allE[,"Rhat"]>=1.1 & allE[,"Rhat"]<=1.5,na.rm=TRUE),
            sum(allE[,"Rhat"]>1.5,na.rm=TRUE), median(allE[,"Rhat"],na.rm=TRUE)))

# publication trace figure
cols <- c("#1b1b1b","#c0392b","#2874a6","#1e8449")
lab  <- c(loglik="Log-likelihood", GDPC1_GPDIC1="GDPC1 - GPDIC1",
          FEDFUNDS_TB3MS="FEDFUNDS - TB3MS", AAA_BAA="AAA - BAA",
          CPIAUCSL_CPITRNSL="CPIAUCSL - CPITRNSL", GS5_GS10="GS5 - GS10")
png(file.path(path.expand(OUT),"diagnostics_p100_trace.png"), width=1500, height=850, res=140)
par(mfrow=c(2,3), mar=c(3.6,3.8,2.8,0.8), mgp=c(2.3,0.7,0), cex.main=1.0,
    cex.lab=0.95, cex.axis=0.85, bty="l", las=1)
for (nm in names(series)){
   x <- series[[nm]]; yr <- range(unlist(x))
   rh <- diag_tab$Rhat[diag_tab$quantity==nm]; es <- diag_tab$ESS[diag_tab$quantity==nm]
   plot(x[[1]], type="l", col=cols[1], lwd=0.5, ylim=yr, xlab="Draw (thinned)",
        ylab=if (nm=="loglik") "Log-likelihood" else "Partial correlation",
        main=bquote(.(lab[[nm]]) ~ ~ (hat(R) == .(sprintf("%.2f",rh)) * ", ESS" ~ .(es))))
   for (c in 2:4) lines(x[[c]], col=cols[c], lwd=0.5) }
dev.off(); cat("Saved diagnostics_p100_trace.png\n")

# ---- C: posterior predictive check ----------------------------------------
fit <- vapor_fits[[1]]; Bm <- apply(fit$B,c(1,2),mean)
E_obs <- Y_train[-1,] - Y_train[-T_train,] %*% t(Bm)
obs_maxeig <- max(eigen(cor(E_obs), symmetric=TRUE, only.values=TRUE)$values)
obs_pcbar <- { PC <- pcor_from_omega(posterior_mean_omega(fit$Lambda, fit$residual_prec))
mean(abs(PC[upper.tri(PC)])) }
nd <- dim(fit$Lambda)[3]; Teff <- nrow(E_obs)
rep_maxeig <- rep_pcbar <- numeric(nd)
set.seed(SEED)
for (s in seq_len(nd)){
   Om <- omega_one(fit,s)
   Sig <- tryCatch(solve(Om), error=function(e) MASS::ginv(Om))
   L <- tryCatch(chol(Sig), error=function(e) chol(Sig + diag(1e-6,p)))
   Erep <- matrix(rnorm(Teff*p),Teff,p) %*% L
   rep_maxeig[s] <- max(eigen(cor(Erep), symmetric=TRUE, only.values=TRUE)$values)
   PC <- pcor_from_omega(Om); rep_pcbar[s] <- mean(abs(PC[upper.tri(PC)])) }
cat("\n=== C: posterior predictive check ===\n")
cat(sprintf("max eig    obs %.2f | rep %.2f [%.2f, %.2f]\n", obs_maxeig,
            mean(rep_maxeig), quantile(rep_maxeig,.025), quantile(rep_maxeig,.975)))
cat(sprintf("mean|pcor| obs %.3f | rep %.3f [%.3f, %.3f]\n", obs_pcbar,
            mean(rep_pcbar), quantile(rep_pcbar,.025), quantile(rep_pcbar,.975)))
png(file.path(path.expand(OUT),"diagnostics_p100_ppc.png"), width=1000, height=450)
par(mfrow=c(1,2))
hist(rep_maxeig, breaks=30, main="PPC: max eigenvalue", xlab="max eig(cor)", col="grey85",
     xlim=range(c(rep_maxeig,obs_maxeig))); abline(v=obs_maxeig, col=2, lwd=2)
hist(rep_pcbar, breaks=30, main="PPC: mean|partial cor|", xlab="mean|pcor|", col="grey85",
     xlim=range(c(rep_pcbar,obs_pcbar))); abline(v=obs_pcbar, col=2, lwd=2)
dev.off(); cat("Saved diagnostics_p100_ppc.png\n")

# ---- D: sensitivity --------------------------------------------------------
sens_k <- do.call(rbind, lapply(c(15,25,35), function(kk){
   f <- vapor(Y_train, lags=1, k=kk, rank=RANK, prior_B="dl", nrun=NRUN, burn=BURN,
              thin=THIN, seed=SEED, jitter=JIT, init="random")
   sel <- select_graph_vapor(f)
   data.frame(setting=sprintf("k=%d",kk),
              MSPE=round(mean(mspe_for_omega(sel$Omega)),4), edges=sum(sel$P)/2)}))
sens_r <- do.call(rbind, lapply(c(2,3,4), function(rr){
   f <- vapor(Y_train, lags=1, k=K, rank=rr, prior_B="dl", nrun=NRUN, burn=BURN,
              thin=THIN, seed=SEED, jitter=JIT, init="random")
   sel <- select_graph_vapor(f)
   data.frame(setting=sprintf("rank=%d",rr),
              MSPE=round(mean(mspe_for_omega(sel$Omega)),4), edges=sum(sel$P)/2)}))
sens <- rbind(sens_k, sens_r)
cat("\n=== D: sensitivity ===\n"); print(sens, row.names=FALSE)
cat(sprintf("MSPE range [%.4f, %.4f] sd %.4f\n", min(sens$MSPE), max(sens$MSPE), sd(sens$MSPE)))
seed_mspe <- sapply(1:6, function(s){
   f <- vapor(Y_train, lags=1, k=K, rank=RANK, prior_B="dl", nrun=NRUN, burn=BURN,
              thin=THIN, seed=s, jitter=JIT, init="random")
   mean(mspe_for_omega(select_graph_vapor(f)$Omega))})
cat(sprintf("Seed: mean %.4f sd %.4f range [%.4f, %.4f]\n",
            mean(seed_mspe), sd(seed_mspe), min(seed_mspe), max(seed_mspe)))

# ---- E: do competing estimators recover the same edges? -------------------
Om_vapor <- posterior_mean_omega(vapor_fits[[1]]$Lambda, vapor_fits[[1]]$residual_prec)
PC100 <- list(VAPOR   = pcor_from_omega(Om_vapor),
              Kastner = pcor_from_omega(Omega_kast),
              Glasso  = pcor_from_omega(Omega_glasso),
              Sample  = pcor_from_omega(Omega_samp))

ut  <- which(upper.tri(matrix(0,p,p)), arr.ind=TRUE)
key <- function(i,j) paste(pmin(i,j), pmax(i,j), sep="-")
rank_tab <- lapply(PC100, function(M){
   v <- abs(M[cbind(ut[,1],ut[,2])]); o <- order(v, decreasing=TRUE)
   rk <- integer(length(v)); rk[o] <- seq_along(v); setNames(rk, key(ut[,1],ut[,2])) })
ekeys <- key(match(ci_edges$v1,vars_qd), match(ci_edges$v2,vars_qd))
edge_ranks <- sapply(rank_tab, function(rk) rk[ekeys])
rownames(edge_ranks) <- paste(ci_edges$v1, ci_edges$v2, sep="-")
nE <- nrow(ci_edges)
cat(sprintf("\n=== E: rank of each reported edge (of %d pairs) ===\n", nrow(ut)))
print(cbind(pcor=ci_edges$pcor, as.data.frame(edge_ranks)))
cat("\nmedian rank:\n"); print(round(apply(edge_ranks,2,median),1))
cat(sprintf("\nin each method's top-%d / top-100 / top-250:\n", nE))
print(sapply(rank_tab, function(rk)
   c(top_nE=sum(rk[ekeys]<=nE), top100=sum(rk[ekeys]<=100), top250=sum(rk[ekeys]<=250))))

# term-structure adjacency
rates <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10")
pres <- apply(ci_edges[,1:2], 1, function(z) paste(sort(z), collapse="|"))
cat("\n=== term-structure adjacency ===\n")
for (i in 1:(length(rates)-1)) for (j in (i+1):length(rates))
   cat(sprintf("  %-9s %-9s gap %d  %s\n", rates[i], rates[j], j-i,
               ifelse(paste(sort(c(rates[i],rates[j])),collapse="|") %in% pres, "EDGE","-")))
cat("\nmarginal correlations among rate series:\n")
print(round(cor(Y_train[, c(rates,"MORTGAGE30US")]), 2))

# ---- save ------------------------------------------------------------------
saveRDS(list(results=results, vapor_times=vapor_times, diag_tab=diag_tab,
             edge_ci=edge_ci, ci_edges=ci_edges, edge_rhat=edge_rhat,
             sens=sens, seed_mspe=seed_mspe,
             ppc=list(obs_maxeig=obs_maxeig, rep_maxeig=rep_maxeig,
                      obs_pcbar=obs_pcbar, rep_pcbar=rep_pcbar),
             raw_omega=list(vapor=Om_vapor, glasso=Omega_glasso,
                            kast=Omega_kast, samp=Omega_samp),
             edge_ranks=edge_ranks,
             loss=list(vapor=v_loss, glasso=mspe_glasso, kast=mspe_kast,
                       samp=mspe_samp, rw=mspe_rw),
             hid_list=hid_list, Y_test=Y_test,
             vars=vars_qd, p=p, T_train=T_train),
        file.path(path.expand(OUT), "p100_full_results.rds"))
cat("\nSaved p100_full_results.rds\n")

# ---- four-panel heatmap ----------------------------------------------------
png(file.path(path.expand(OUT),"heatmap_p100_fullrank.png"), width=1600, height=440)
par(mfrow=c(1,4), mar=c(2,2,3,1))
offd <- upper.tri(PC100$VAPOR)
rng <- max(abs(c(PC100$VAPOR[offd], PC100$Kastner[offd], PC100$Glasso[offd], PC100$Sample[offd])))
col_pc <- colorRampPalette(c("blue","white","red"))(101)
imgpc <- function(M,ttl){ diag(M) <- 0
image(t(M[nrow(M):1,]), col=col_pc, zlim=c(-rng,rng), axes=FALSE, main=ttl) }
imgpc(PC100$VAPOR,"VAPOR"); imgpc(PC100$Kastner,"Kastner")
imgpc(PC100$Glasso,"VAR + Glasso"); imgpc(PC100$Sample,"Sample covariance")
dev.off(); cat("Saved heatmap_p100_fullrank.png\n")

# ---- edge-count table ------------------------------------------------------
thresholds <- c(0.05,0.1,0.2,0.3,0.4,0.5)
count_edges <- function(PC,thr){ off <- PC[upper.tri(PC)]; sum(abs(off) > thr) }
edge_tab <- data.frame(threshold=thresholds,
                       VAPOR     = sapply(thresholds, function(t) count_edges(PC100$VAPOR,t)),
                       Kastner   = sapply(thresholds, function(t) count_edges(PC100$Kastner,t)),
                       Glasso    = sapply(thresholds, function(t) count_edges(PC100$Glasso,t)),
                       SampleCov = sapply(thresholds, function(t) count_edges(PC100$Sample,t)))
cat("\nEdge counts (|pcor| > threshold), out of", choose(p,2), "pairs:\n")
print(edge_tab, row.names=FALSE)
r <- readRDS(file.path(path.expand(OUT),"p100_full_results.rds"))
r$edge_tab <- edge_tab; saveRDS(r, file.path(path.expand(OUT),"p100_full_results.rds"))

# ---- term-structure figure -------------------------------------------------
rates6 <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTGAGE30US")
lab6   <- c("FEDFUNDS","TB3MS","GS1","GS5","GS10","MORTG30")
Cm <- cor(Y_train[, rates6])
pc_of <- function(a,b){ i <- which((ci_edges$v1==a & ci_edges$v2==b) |
                                      (ci_edges$v1==b & ci_edges$v2==a))
if (length(i)) ci_edges$pcor[i] else NA_real_ }
edges <- data.frame(from=c("FEDFUNDS","TB3MS","GS1","GS5","GS5"),
                    to  =c("TB3MS","GS1","GS5","GS10","MORTGAGE30US"),
                    stringsAsFactors=FALSE)
edges$pcor <- mapply(pc_of, edges$from, edges$to)

png(file.path(path.expand(OUT),"termstructure_p100.png"), width=1900, height=860, res=190)
layout(matrix(1:2, nrow=1), widths=c(1,1.15))
par(mar=c(1.5,6.2,6.2,1.5))
n <- length(rates6); pal <- colorRampPalette(c("white","#F4A582","#B2182B"))(101)
plot(NA, xlim=c(0.5,n+0.5), ylim=c(0.5,n+0.5), axes=FALSE, xlab="", ylab="", asp=1, main="")
for (i in 1:n) for (j in 1:n){ v <- Cm[i,j]
rect(j-0.5, n-i+0.5, j+0.5, n-i+1.5,
     col=pal[1+round(100*max(0,min(1,abs(v))))], border="white")
text(j, n-i+1, sprintf("%.2f",v), cex=0.62, col=if (abs(v)>0.7) "white" else "grey20") }
axis(2, at=n:1, labels=lab6, tick=FALSE, las=1, cex.axis=0.72, line=-0.6)
axis(3, at=1:n, labels=lab6, tick=FALSE, las=2, cex.axis=0.72, line=-0.9)
mtext("Marginal correlation", side=3, line=4.6, font=2, cex=1.02)

par(mar=c(1.5,1,6.2,1))
pos <- rbind(FEDFUNDS=c(1,1), TB3MS=c(2,1), GS1=c(3,1), GS5=c(4,1), GS10=c(5,1),
             MORTGAGE30US=c(4,0.25))
plot(NA, xlim=c(0.45,5.55), ylim=c(0.05,1.30), axes=FALSE, xlab="", ylab="", main="")
for (k2 in seq_len(nrow(edges))){
   a <- pos[edges$from[k2],]; b <- pos[edges$to[k2],]
   segments(a[1],a[2],b[1],b[2], lwd=1+9*abs(edges$pcor[k2]), col="#2C6FBB")
   if (edges$to[k2]=="MORTGAGE30US")
      text(b[1]+0.07, (a[2]+b[2])/2, sprintf("%.2f",edges$pcor[k2]), cex=0.68, col="#1B4F72", pos=4)
   else text((a[1]+b[1])/2, a[2]+0.11, sprintf("%.2f",edges$pcor[k2]), cex=0.68, col="#1B4F72") }
cols6 <- ifelse(rownames(pos)=="MORTGAGE30US", "#59A14F", "#EDC948")
points(pos[,1], pos[,2], pch=21, bg=cols6, col="grey25", cex=4.4, lwd=1.2)
text(pos[,1], pos[,2]-0.155, lab6, cex=0.68)
mtext("Conditional dependence (VAPOR)", side=3, line=4.6, font=2, cex=1.02)
mtext("adjacent maturities only; all six non-adjacent pairs absent",
      side=3, line=3.1, cex=0.72, col="grey35")
dev.off(); cat("Saved termstructure_p100.png\n")