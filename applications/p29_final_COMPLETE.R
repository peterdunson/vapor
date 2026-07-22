# ============================================================================
#  VAPOR p=29 FRED-QD — full rerun (full-rank DL proposed) + diagnostics A-E
# ============================================================================
library(VAPOR); library(PrecFact); library(bayesianVARs); library(glasso); library(loo)

SEED <- 32; R_REPS <- 3; NRUN <- 7500; BURN <- 2500; THIN <- 5
OUT  <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application"

fred_qd_raw <- read.csv("~/Desktop/2026-05-qd.csv")
transform_series <- function(x, tcode){ x<-as.numeric(x)
switch(as.character(tcode),"1"=x,"2"=c(NA,diff(x)),"3"=c(NA,NA,diff(diff(x))),
       "4"=log(x),"5"=c(NA,diff(log(x))),"6"=c(NA,NA,diff(diff(log(x)))),
       "7"=c(NA,diff(log(x)/lag(log(x))-1)),x) }
tcodes_qd<-as.numeric(fred_qd_raw[2,-1]); names(tcodes_qd)<-colnames(fred_qd_raw)[-1]
fred_qd_data<-fred_qd_raw[-c(1,2),]; dates_qd<-fred_qd_data[,1]
fred_qd_data<-apply(fred_qd_data[,-1],2,as.numeric)
fred_qd_trans<-mapply(transform_series,as.data.frame(fred_qd_data),tcodes_qd)
fred_qd_trans<-as.matrix(fred_qd_trans); rownames(fred_qd_trans)<-dates_qd
fred_qd_clean<-fred_qd_trans[apply(fred_qd_trans,1,function(r)!any(is.na(r))),]
fred_qd_clean<-fred_qd_clean[,apply(fred_qd_clean,2,function(c)!any(is.na(c)))]

vars_qd<-c("GDPC1","INDPRO","IPFINAL","IPCONGD","IPBUSEQ","PAYEMS","UEMPMEAN",
           "CLAIMSx","CUMFNS","PCECC96","CMRMTSPLx","BUSINVx","CPIAUCSL","PPICMM","GDPCTPI",
           "OILPRICEx","FEDFUNDS","GS1","GS10","BAA","TB3MS","M2REAL","BUSLOANSx","NONREVSLx",
           "S.P.500","VIXCLSx","GPDIC1","UNRATE","WPSFD49207")
vars_qd<-vars_qd[vars_qd %in% colnames(fred_qd_clean)]
Y_qd<-fred_qd_clean[,vars_qd]; p<-ncol(Y_qd)

parse_qd_date<-function(s){z<-as.numeric(strsplit(s,"/")[[1]]);z[3]+(z[1]-1)/12}
date_vals<-vapply(rownames(Y_qd),parse_qd_date,numeric(1))
train_idx<-which(date_vals>=2000+(6-1)/12 & date_vals<=2020+(12-1)/12)
test_idx <-which(date_vals>=2021+(3-1)/12 & date_vals<=2025+(12-1)/12)
Y_all<-scale(Y_qd); Y_train<-Y_all[train_idx,]; Y_test<-Y_all[test_idx,]
T_train<-nrow(Y_train); T_test<-nrow(Y_test); Y_last_train<-Y_train[T_train,]
cat(sprintf("p=%d  T_train=%d  T/p=%.2f\n",p,T_train,T_train/p))

predict_held_out<-function(Omega,y_obs,obs_idx,hid_idx)
   as.vector(-solve(Omega[hid_idx,hid_idx])%*%Omega[hid_idx,obs_idx]%*%y_obs)
set.seed(SEED); n_hid<-floor(0.3*p)
hid_list<-lapply(seq_len(T_test),function(t)sort(sample(p,n_hid)))
mspe_for_omega<-function(Omega) sapply(seq_len(T_test),function(t){
   hid<-hid_list[[t]]; obs<-setdiff(seq_len(p),hid); yt<-Y_test[t,]
   mean((yt[hid]-predict_held_out(Omega,yt[obs],obs,hid))^2)})

# deterministic baselines
set.seed(SEED)
Y_train_jit<-Y_train+matrix(rnorm(prod(dim(Y_train)),0,1e-3),T_train,p)
mod_var<-bvar(Y_train_jit,lags=1L,draws=2000,burnin=1000,quiet=TRUE)
B_hat<-apply(mod_var$PHI,c(1,2),mean); resid_var<-mod_var$Y-mod_var$X%*%B_hat
S_resid<-cov(scale(resid_var)); rho_grid<-10^seq(-2,0,length=30)
bic_vals<-sapply(rho_grid,function(rho){gl<-tryCatch(glasso(S_resid,rho=rho,penalize.diagonal=FALSE),error=function(e)NULL)
if(is.null(gl))return(Inf); Om<-gl$wi; ne<-sum(Om[upper.tri(Om)]!=0)
ll<-(T_train/2)*(log(det(Om))-sum(diag(S_resid%*%Om))); -2*ll+log(T_train)*ne})
gl_var<-glasso(S_resid,rho=rho_grid[which.min(bic_vals)],penalize.diagonal=FALSE)
Omega_glasso<-gl_var$wi; glasso_edges<-sum(gl_var$wi[upper.tri(gl_var$wi)]!=0)
sigma_prior<-specify_prior_sigma(M=p,type="factor",factor_factors=5L,
                                 factor_priorfacloadtype="colwiseng",factor_heteroskedastic=TRUE,quiet=TRUE)
mod_kast<-bvar(Y_train_jit,lags=1L,draws=2000,burnin=1000,prior_sigma=sigma_prior,quiet=TRUE)
vc<-vcov(mod_kast); Omega_kast<-solve(apply(vc[1,,,],c(1,2),mean))
Omega_samp<-solve(cov(Y_train))
mspe_glasso<-mspe_for_omega(Omega_glasso); mspe_kast<-mspe_for_omega(Omega_kast)
mspe_samp<-mspe_for_omega(Omega_samp)
mspe_rw<-sapply(seq_len(T_test),function(t){hid<-hid_list[[t]]
y_prev<-if(t==1)Y_last_train else Y_test[t-1,]; mean((Y_test[t,hid]-y_prev[hid])^2)})

select_graph_vapor<-function(fit,alpha=0.1,beta=0.85,
                             epsilon_grid=seq(0.05,0.1,length.out=20),burnout=0){
   res_like<-list(Lambda=fit$Lambda,residuals=fit$residual_prec)
   out<-fact_prec_processing(res_like,burnout=burnout,alpha=alpha,beta=beta,epsilon_grid=epsilon_grid)
   diag(out$P)<-0; out }

# stochastic methods averaged over R_REPS; keep VAPOR fits for diagnostics
pf_mspe<-pf_edges<-pfvar_mspe<-pfvar_edges<-v_mspe<-v_edges<-numeric(R_REPS)
vapor_fits<-vector("list",R_REPS); vapor_times<-numeric(R_REPS)
for(r in seq_len(R_REPS)){
   cat(sprintf("\rRep %d/%d",r,R_REPS)); flush.console()
   res_pf<-prec_est(as=0.1,bs=0.1,a=0.5,b=0.5,k=5,nrun=NRUN,thin=THIN,Y=Y_train,alldiff=0,a_dir=0.1,b_dir=0.1)
   out_pf<-fact_prec_processing(res_pf,burnout=750,alpha=0.1,beta=0.85,epsilon_grid=seq(0.05,0.1,length=20)); diag(out_pf$P)<-0
   pf_mspe[r]<-mean(mspe_for_omega(out_pf$Omega)); pf_edges[r]<-sum(out_pf$P)/2
   res_pfv<-prec_est(as=0.1,bs=0.1,a=0.5,b=0.5,k=5,nrun=NRUN,thin=THIN,Y=scale(resid_var),alldiff=0,a_dir=0.1,b_dir=0.1)
   out_pfv<-fact_prec_processing(res_pfv,burnout=750,alpha=0.1,beta=0.85,epsilon_grid=seq(0.05,0.1,length=20)); diag(out_pfv$P)<-0
   pfvar_mspe[r]<-mean(mspe_for_omega(out_pfv$Omega)); pfvar_edges[r]<-sum(out_pfv$P)/2
   tt<-system.time( fit<-vapor(Y_train,lags=1,k=5,rank="full",prior_B="dl",nrun=NRUN,burn=BURN,thin=THIN,seed=SEED+r) )
   vapor_times[r]<-tt["elapsed"]
   sel<-select_graph_vapor(fit); v_mspe[r]<-mean(mspe_for_omega(sel$Omega)); v_edges[r]<-sum(sel$P)/2
   vapor_fits[[r]]<-fit
}
cat("\n")
fmt<-function(v)sprintf("%.4f (%.4f)",mean(v),sd(v))
results<-data.frame(
   Method=c("PrecFact (no VAR)","VAR + Glasso","VAR + PrecFact","Kastner (VAR+FSV)",
            "Sample covariance","Random walk","VAPOR (proposed)"),
   MSPE=c(fmt(pf_mspe),sprintf("%.4f",mean(mspe_glasso)),fmt(pfvar_mspe),
          sprintf("%.4f",mean(mspe_kast)),sprintf("%.4f",mean(mspe_samp)),
          sprintf("%.4f",mean(mspe_rw)),fmt(v_mspe)))
cat("\n=== p=29 MSPE ===\n"); print(results,row.names=FALSE)
cat(sprintf("VAPOR edges (pFDR): %.1f (sd %.1f)\n",mean(v_edges),sd(v_edges)))
cat(sprintf("VAPOR mean fit time: %.1f s\n",mean(vapor_times)))


# ---- A: MCMC diagnostics ---------------------------------------------------
library(coda)
omega_one <- function(fit, s){
   Om <- tcrossprod(fit$Lambda[, , s])
   diag(Om) <- diag(Om) + fit$residual_prec[s, ]
   Om
}
pcor_from_omega<-function(Om){ d<-sqrt(diag(Om)); pc<- -Om/outer(d,d); diag(pc)<-1; pc }

functionals<-function(fit){
   nd<-dim(fit$Lambda)[3]
   pairs<-list(TB3MS_GS1  =c("TB3MS","GS1"),
               GS10_BAA   =c("GS10","BAA"),
               FEDFUNDS_TB3MS=c("FEDFUNDS","TB3MS"),
               CPI_PPI    =c("CPIAUCSL","WPSFD49207"),
               CPI_OIL    =c("CPIAUCSL","OILPRICEx"))
   idx<-lapply(pairs,function(pr) c(which(vars_qd==pr[1]),which(vars_qd==pr[2])))
   mat<-matrix(0,nd,length(pairs)); colnames(mat)<-names(pairs)
   for(s in seq_len(nd)){
      Om<-omega_one(fit,s); PC<-pcor_from_omega(Om)
      for(j in seq_along(idx)) mat[s,j]<-PC[idx[[j]][1],idx[[j]][2]]
   }
   data.frame(loglik=colSums(fit$loglik), mat)
}
set.seed(1)
chains<-lapply(1:4,function(s) vapor(Y_train,lags=1,k=5,rank="full",prior_B="dl",
                                     nrun=NRUN,burn=BURN,thin=THIN,seed=100+s))
chain_fns<-lapply(chains,functionals)

mc_list<-mcmc.list(lapply(chain_fns,function(d) mcmc(as.matrix(d))))
ess_vals<-effectiveSize(mc_list)
rhat_vals<-tryCatch(gelman.diag(mc_list,autoburnin=FALSE,multivariate=FALSE)$psrf[,1],
                    error=function(e) rep(NA,ncol(chain_fns[[1]])))
cat("\n=== A: MCMC diagnostics (4 chains) ===\n")
diag_tab<-data.frame(quantity=names(ess_vals),
                     ESS=round(ess_vals,0), Rhat=round(rhat_vals,3))
print(diag_tab,row.names=FALSE)
cat(sprintf("\nComp time: mean %.1f s/fit (n=%d), total %.1f s\n",
            mean(vapor_times),length(vapor_times),sum(vapor_times)))

plot_names <- c("loglik","TB3MS_GS1","GS10_BAA","FEDFUNDS_TB3MS","CPI_PPI","CPI_OIL")
png(file.path(path.expand(OUT),"diagnostics_p29_trace.png"),width=1400,height=900)
par(mfrow=c(2,3))
for(nm in plot_names){
   yr <- range(sapply(chain_fns, function(d) d[[nm]]))
   plot(chain_fns[[1]][[nm]],type="l",col=1,main=nm,xlab="draw",ylab=nm,ylim=yr)
   for(cc in 2:4) lines(chain_fns[[cc]][[nm]],col=cc)
}
dev.off()
cat("Saved diagnostics_p29_trace.png\n")


# ---- B: 95% credible-interval edge selection -------------------------------
all_fits <- c(vapor_fits, chains)
pc_stack <- list()
for(fit in all_fits){
   nd<-dim(fit$Lambda)[3]
   for(s in seq_len(nd)){
      Om<-omega_one(fit,s)
      pc_stack[[length(pc_stack)+1]]<-pcor_from_omega(Om)
   }
}
PCarr<-simplify2array(pc_stack)
lo<-apply(PCarr,c(1,2),quantile,0.025)
hi<-apply(PCarr,c(1,2),quantile,0.975)
edge_ci<-(lo>0 | hi<0); diag(edge_ci)<-FALSE
cat(sprintf("\n=== B: 95%% interval edges ===\nedges: %d (of %d possible)\n",
            sum(edge_ci)/2, p*(p-1)/2))
idx<-which(edge_ci & upper.tri(edge_ci),arr.ind=TRUE)
if(nrow(idx)>0){
   pm<-apply(PCarr,c(1,2),mean)
   ci_edges<-data.frame(v1=vars_qd[idx[,1]],v2=vars_qd[idx[,2]],
                        pcor=round(pm[idx],3))
   ci_edges<-ci_edges[order(-abs(ci_edges$pcor)),]
   print(head(ci_edges,20),row.names=FALSE)
}

# ---- C: posterior predictive check -----------------------------------------
fit<-vapor_fits[[1]]
Bm<-apply(fit$B,c(1,2),mean)
E_obs<-Y_train[-1,]-Y_train[-T_train,]%*%t(Bm)
obs_maxeig<-max(eigen(cor(E_obs),symmetric=TRUE,only.values=TRUE)$values)
obs_pcbar<-{PC<-pcor_from_omega(posterior_mean_omega(fit$Lambda, fit$residual_prec))
mean(abs(PC[upper.tri(PC)]))}

nd<-dim(fit$Lambda)[3]; Teff<-nrow(E_obs)
rep_maxeig<-numeric(nd); rep_pcbar<-numeric(nd)
set.seed(SEED)
for(s in seq_len(nd)){
   Om<-omega_one(fit,s)
   Sig<-solve(Om); L<-chol(Sig)
   Erep<-matrix(rnorm(Teff*p),Teff,p)%*%L
   rep_maxeig[s]<-max(eigen(cor(Erep),symmetric=TRUE,only.values=TRUE)$values)
   PC<-pcor_from_omega(Om); rep_pcbar[s]<-mean(abs(PC[upper.tri(PC)]))
}
cat("\n=== C: posterior predictive check ===\n")
cat(sprintf("max eigenvalue  observed %.2f | replicated %.2f [%.2f, %.2f]\n",
            obs_maxeig,mean(rep_maxeig),quantile(rep_maxeig,.025),quantile(rep_maxeig,.975)))
cat(sprintf("mean|pcor|      observed %.3f | replicated %.3f [%.3f, %.3f]\n",
            obs_pcbar,mean(rep_pcbar),quantile(rep_pcbar,.025),quantile(rep_pcbar,.975)))
png(file.path(path.expand(OUT),"diagnostics_p29_ppc.png"),width=1000,height=450)
par(mfrow=c(1,2))
hist(rep_maxeig,breaks=30,main="PPC: max eigenvalue",xlab="max eig(cor)",col="grey85",
     xlim=range(c(rep_maxeig,obs_maxeig))); abline(v=obs_maxeig,col=2,lwd=2)
hist(rep_pcbar,breaks=30,main="PPC: mean|partial cor|",xlab="mean|pcor|",col="grey85",
     xlim=range(c(rep_pcbar,obs_pcbar))); abline(v=obs_pcbar,col=2,lwd=2)
dev.off()
cat("Saved diagnostics_p29_ppc.png (red line = observed)\n")


# ---- D: sensitivity to tunable settings ------------------------------------
sens_k <- do.call(rbind, lapply(c(3,5,8), function(kk){
   fit <- vapor(Y_train, lags=1, k=kk, rank="full", prior_B="dl",
                nrun=NRUN, burn=BURN, thin=THIN, seed=SEED)
   sel <- select_graph_vapor(fit)
   data.frame(setting=sprintf("k=%d", kk),
              MSPE=round(mean(mspe_for_omega(sel$Omega)),4), edges=sum(sel$P)/2)
}))

sens_mcmc <- do.call(rbind, lapply(list(c(5000,2500), c(7500,2500), c(15000,5000)),
                                   function(nb){
                                      fit <- vapor(Y_train, lags=1, k=5, rank="full", prior_B="dl",
                                                   nrun=nb[1], burn=nb[2], thin=THIN, seed=SEED)
                                      sel <- select_graph_vapor(fit)
                                      data.frame(setting=sprintf("nrun=%d", nb[1]),
                                                 MSPE=round(mean(mspe_for_omega(sel$Omega)),4), edges=sum(sel$P)/2)
                                   }))

aL_vals <- c(0.3, 0.5, 1.0)
sens_aL <- do.call(rbind, lapply(aL_vals, function(aL){
   fit <- vapor(Y_train, lags=1, k=5, rank="full", prior_B="dl",
                a_L=aL, b_L=aL, nrun=NRUN, burn=BURN, thin=THIN, seed=SEED)
   sel <- select_graph_vapor(fit)
   data.frame(setting=sprintf("a_L=%.1f", aL),
              MSPE=round(mean(mspe_for_omega(sel$Omega)),4), edges=sum(sel$P)/2)
}))
aL_tunable <- diff(range(sens_aL$MSPE)) > 1e-6

sens <- rbind(sens_k, sens_mcmc, if (aL_tunable) sens_aL else NULL)
cat("\n=== D: sensitivity (tunable settings) ===\n")
print(sens, row.names=FALSE)
if (!aL_tunable) cat("(a_L had no effect on MSPE; loading-prior shape is effectively fixed)\n")
cat(sprintf("MSPE range [%.4f, %.4f], sd %.4f across settings\n",
            min(sens$MSPE), max(sens$MSPE), sd(sens$MSPE)))

seed_mspe <- sapply(1:6, function(s){
   fit <- vapor(Y_train, lags=1, k=5, rank="full", prior_B="dl",
                nrun=NRUN, burn=BURN, thin=THIN, seed=s)
   mean(mspe_for_omega(select_graph_vapor(fit)$Omega))
})
cat(sprintf("Seed sensitivity: mean %.4f sd %.4f range [%.4f, %.4f]\n",
            mean(seed_mspe), sd(seed_mspe), min(seed_mspe), max(seed_mspe)))


# ---- E: do competing estimators recover the same edges? --------------------
Om_vapor <- posterior_mean_omega(vapor_fits[[1]]$Lambda, vapor_fits[[1]]$residual_prec)
PC29 <- list(VAPOR   = pcor_from_omega(Om_vapor),
             Kastner = pcor_from_omega(Omega_kast),
             Glasso  = pcor_from_omega(Omega_glasso),
             Sample  = pcor_from_omega(Omega_samp))

# reported FDR edge set, regenerated from the fit
sel29 <- select_graph_vapor(vapor_fits[[1]])
ei    <- which(sel29$P == 1 & upper.tri(sel29$P), arr.ind = TRUE)
ed29  <- data.frame(v1 = vars_qd[ei[,1]], v2 = vars_qd[ei[,2]],
                    pcor = round(PC29$VAPOR[ei], 3), stringsAsFactors = FALSE)
ed29  <- ed29[order(-abs(ed29$pcor)), ]
cat("\n=== E: FDR-selected edges (regenerated from fit) ===\n")
print(ed29, row.names = FALSE)

ut  <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
key <- function(i,j) paste(pmin(i,j), pmax(i,j), sep = "-")
rank_tab <- lapply(PC29, function(M){
   v <- abs(M[cbind(ut[,1], ut[,2])]); o <- order(v, decreasing = TRUE)
   rk <- integer(length(v)); rk[o] <- seq_along(v); setNames(rk, key(ut[,1], ut[,2]))
})
ekeys <- key(match(ed29$v1, vars_qd), match(ed29$v2, vars_qd))
edge_ranks <- sapply(rank_tab, function(rk) rk[ekeys])
rownames(edge_ranks) <- paste(ed29$v1, ed29$v2, sep = "-")
nE <- nrow(ed29)

cat(sprintf("\n=== rank of each reported edge (of %d pairs) ===\n", nrow(ut)))
print(cbind(pcor = ed29$pcor, as.data.frame(edge_ranks)))
cat("\nmedian rank:\n"); print(round(apply(edge_ranks, 2, median), 1))
cat(sprintf("\nin each method's top-%d / top-25 / top-50:\n", nE))
print(sapply(rank_tab, function(rk)
   c(top_nE = sum(rk[ekeys] <= nE), top25 = sum(rk[ekeys] <= 25),
     top50  = sum(rk[ekeys] <= 50))))

# same comparison for the 95% CI edge set
if (exists("ci_edges") && nrow(ci_edges) > 0) {
   ckeys <- key(match(ci_edges$v1, vars_qd), match(ci_edges$v2, vars_qd))
   ci_ranks <- sapply(rank_tab, function(rk) rk[ckeys])
   rownames(ci_ranks) <- paste(ci_edges$v1, ci_edges$v2, sep = "-")
   cat("\n=== rank of each 95% CI edge ===\n")
   print(cbind(pcor = ci_edges$pcor, as.data.frame(ci_ranks)))
   cat("\nmedian rank (CI edges):\n"); print(round(apply(ci_ranks, 2, median), 1))
} else ci_ranks <- NULL

# rate block: chain or near-complete?
rates29 <- intersect(c("FEDFUNDS","TB3MS","GS1","GS5","GS10"), vars_qd)
cat("\n=== rate series present at p=29:", paste(rates29, collapse=", "), "===\n")
pres <- apply(ed29[,1:2], 1, function(z) paste(sort(z), collapse="|"))
for (i in seq_len(length(rates29)-1)) for (j in (i+1):length(rates29))
   cat(sprintf("  %-9s %-9s  %s\n", rates29[i], rates29[j],
               ifelse(paste(sort(c(rates29[i], rates29[j])), collapse="|") %in% pres, "EDGE", "-")))
cat(sprintf("rate-block density: %d of %d possible pairs\n",
            sum(pres %in% apply(t(combn(rates29,2)),1,function(z) paste(sort(z),collapse="|"))),
            choose(length(rates29),2)))

# marginal correlations among the rate series, for contrast
cat("\nmarginal correlations among rate series:\n")
print(round(cor(Y_train[, rates29]), 2))


# ---- save ------------------------------------------------------------------
saveRDS(list(results=results, vapor_times=vapor_times, diag_tab=diag_tab,
             edge_ci=edge_ci, ci_edges=if(exists("ci_edges")) ci_edges else NULL,
             sens=sens, seed_mspe=seed_mspe, aL_tunable=aL_tunable,
             ppc=list(obs_maxeig=obs_maxeig, rep_maxeig=rep_maxeig,
                      obs_pcbar=obs_pcbar, rep_pcbar=rep_pcbar),
             raw_omega=list(vapor=Om_vapor, glasso=Omega_glasso,
                            kast=Omega_kast, samp=Omega_samp),
             fdr_edges=ed29, edge_ranks=edge_ranks, ci_ranks=ci_ranks,
             vars=vars_qd, p=p, T_train=T_train),
        file.path(path.expand(OUT), "p29_full_results.rds"))
cat("\nSaved p29_full_results.rds\n")


# ============================================================================
#  p=29 edge-count table + four-panel heatmap
# ============================================================================
library(pheatmap); library(gridExtra); library(grid)

PC_vapor <- PC29$VAPOR
PC_kast  <- PC29$Kastner
PC_glas  <- PC29$Glasso
PC_samp  <- PC29$Sample

thresholds <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
count_edges <- function(PC, thr){ off <- PC[upper.tri(PC)]; sum(abs(off) > thr) }
edge_tab <- data.frame(
   threshold = thresholds,
   VAPOR     = sapply(thresholds, function(t) count_edges(PC_vapor, t)),
   Kastner   = sapply(thresholds, function(t) count_edges(PC_kast, t)),
   Glasso    = sapply(thresholds, function(t) count_edges(PC_glas, t)),
   SampleCov = sapply(thresholds, function(t) count_edges(PC_samp, t))
)
cat("\nEdge counts (|pcor| > threshold), out of", choose(p,2), "pairs:\n")
print(edge_tab, row.names=FALSE)

mk_pc <- function(PC){ M<-PC; diag(M)<-0; M }
pal  <- colorRampPalette(c("blue","white","red"))(100); brks <- seq(-1,1,length=101)
mk <- function(M,ttl,leg=FALSE) pheatmap(mk_pc(M), cluster_rows=FALSE, cluster_cols=FALSE,
                                         silent=TRUE, color=pal, breaks=brks, show_rownames=FALSE, show_colnames=FALSE,
                                         legend=leg, main=ttl)$gtable

png(file.path(path.expand(OUT), "heatmap_p29_fullrank.png"), width=1600, height=440)
grid.arrange(mk(PC_vapor,"VAPOR"), mk(PC_kast,"Kastner"), mk(PC_glas,"VAR + Glasso"),
             mk(PC_samp,"Sample covariance", TRUE), ncol=4)
dev.off()
cat("Saved heatmap_p29_fullrank.png\n")

# append the edge-count table to the rds
r29 <- readRDS(file.path(path.expand(OUT), "p29_full_results.rds"))
r29$edge_tab <- edge_tab
saveRDS(r29, file.path(path.expand(OUT), "p29_full_results.rds"))
cat("Appended edge_tab to p29_full_results.rds\n")