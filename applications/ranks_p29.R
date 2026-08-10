# ============================================================================
#  Rank/k table for Supplement S1, with MSPE computed the same way as Table 1
#  (averaged over 20 draws of the held-out sets).
# ============================================================================
library(VAPOR); library(PrecFact); library(MASS)

SEED <- 32; NRUN <- 7500; BURN <- 2500; THIN <- 5; JIT <- 1e-3
SEEDS <- c(33,34,35); K_GRID <- c(5,8,10,15); N_HID <- 20
OUT <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")

# ---- data ------------------------------------------------------------------
raw <- read.csv("~/Desktop/2026-05-qd.csv")
tf <- function(x,tc){ x<-as.numeric(x)
switch(as.character(tc),"1"=x,"2"=c(NA,diff(x)),"3"=c(NA,NA,diff(diff(x))),
       "4"=log(x),"5"=c(NA,diff(log(x))),"6"=c(NA,NA,diff(diff(log(x)))),
       "7"=c(NA,diff(log(x)/lag(log(x))-1)),x) }
tc <- as.numeric(raw[2,-1]); names(tc) <- colnames(raw)[-1]
d <- raw[-c(1,2),]; dts <- d[,1]; d <- apply(d[,-1],2,as.numeric)
TRm <- suppressWarnings(as.matrix(mapply(tf, as.data.frame(d), tc))); rownames(TRm) <- dts
cl <- TRm[apply(TRm,1,function(r)!any(is.na(r))),]
cl <- cl[,apply(cl,2,function(c)!any(is.na(c)))]
vars <- c("GDPC1","INDPRO","IPFINAL","IPCONGD","IPBUSEQ","PAYEMS","UEMPMEAN","CLAIMSx",
          "CUMFNS","PCECC96","CMRMTSPLx","BUSINVx","CPIAUCSL","PPICMM","GDPCTPI","OILPRICEx",
          "FEDFUNDS","GS1","GS10","BAA","TB3MS","M2REAL","BUSLOANSx","NONREVSLx","S.P.500",
          "VIXCLSx","GPDIC1","UNRATE","WPSFD49207")
vars <- vars[vars %in% colnames(cl)]
Y <- cl[,vars]; p <- ncol(Y)
pd <- function(s){ z<-as.numeric(strsplit(s,"/")[[1]]); z[3]+(z[1]-1)/12 }
dv <- vapply(rownames(Y), pd, numeric(1)); Ya <- scale(Y)
Ytr <- Ya[which(dv>=2000+5/12 & dv<=2020+11/12),]
Yte <- Ya[which(dv>=2021+2/12 & dv<=2025+11/12),]
Tn <- nrow(Ytr); Tte <- nrow(Yte)

# ---- helpers ---------------------------------------------------------------
om1 <- function(f,s){ O <- tcrossprod(f$Lambda[,,s]); diag(O) <- diag(O)+f$residual_prec[s,]; O }
oh  <- function(f) posterior_mean_omega(f$Lambda, f$residual_prec)
pc  <- function(O){ dd<-sqrt(diag(O)); P<- -O/outer(dd,dd); diag(P)<-1; P }
fixB<- function(B){ if(nrow(B)!=p) B<-t(B); if(ncol(B)!=p) B<-B[,seq_len(p),drop=FALSE]; B }
sr  <- function(B) max(abs(eigen(B, only.values=TRUE)$values))

set.seed(SEED); nh <- floor(0.3*p)
hs <- lapply(seq_len(N_HID), function(b) lapply(seq_len(Tte), function(t) sort(sample(p,nh))))
mspe1 <- function(Om,hl) mean(vapply(seq_len(Tte), function(t){
   h<-hl[[t]]; o<-setdiff(seq_len(p),h); yt<-Yte[t,]
   mean((yt[h] + as.vector(solve(Om[h,h]) %*% Om[h,o] %*% yt[o]))^2)}, numeric(1)))
mspe_avg <- function(Om) mean(vapply(hs, function(hl) mspe1(Om,hl), numeric(1)))

ci_n <- function(f){ nd<-dim(f$Lambda)[3]
A <- simplify2array(lapply(seq_len(nd), function(s) pc(om1(f,s))))
lo<-apply(A,c(1,2),quantile,.025); hi<-apply(A,c(1,2),quantile,.975)
E<-(lo>0|hi<0); diag(E)<-FALSE; sum(E)/2 }

pr_ex <- function(f){ nd<-dim(f$B)[3]; ix<-unique(round(seq(1,nd,length.out=min(400,nd))))
rs <- vapply(ix, function(s) sr(fixB(f$B[,,s])), numeric(1))
c(pr=mean(rs>1), med=median(rs)) }

maxeig <- function(M) max(eigen(cor(M), symmetric=TRUE, only.values=TRUE)$values)
ppc_eig <- function(f, n_draws=200, seed=1){
   nd<-dim(f$Lambda)[3]; ix<-unique(round(seq(1,nd,length.out=min(n_draws,nd))))
   Yl<-Ytr[-Tn,,drop=FALSE]; obs<-maxeig(Ytr[-1,,drop=FALSE]); set.seed(seed)
   r <- vapply(ix, function(s){
      Sg<-tryCatch(chol2inv(chol(om1(f,s))),error=function(e) ginv(om1(f,s)))
      L <-tryCatch(chol(Sg),error=function(e) chol(Sg+diag(1e-8,p)))
      maxeig(Yl %*% t(fixB(f$B[,,s])) + matrix(rnorm((Tn-1)*p),Tn-1,p) %*% L)}, numeric(1))
   c(rep=mean(r), ppp=mean(r>=obs), obs=obs) }

fitv <- function(kk,rk,sd){ f<-NULL; jj<-JIT; a<-0
while(is.null(f) && a<4){ a<-a+1
args<-list(Ytr,lags=1,k=kk,rank=rk,prior_B="dl",nrun=NRUN,burn=BURN,thin=THIN,seed=sd)
if(a>1) args$jitter<-jj
f<-tryCatch(do.call(vapor,args), error=function(e) NULL); jj<-jj*4 }; f }

# ---- sweep -----------------------------------------------------------------
res <- NULL
for (rk in list("full",2,3)) for (kk in K_GRID) for (sd in SEEDS){
   lab <- if(identical(rk,"full")) "Full" else as.character(rk)
   cat(sprintf("\r rank=%-4s k=%-2d seed=%d   ", lab, kk, sd)); flush.console()
   f <- fitv(kk, rk, sd); if (is.null(f)) next
   re <- pr_ex(f); pe <- ppc_eig(f)
   res <- rbind(res, data.frame(rank=lab, k=kk, seed=sd,
                                pr_expl=re["pr"], med_rho=re["med"], rep_eig=pe["rep"], ppp=pe["ppp"],
                                obs_eig=pe["obs"], mspe=mspe_avg(oh(f)), edges=ci_n(f)))
}
cat("\n")
agg <- aggregate(cbind(pr_expl, med_rho, rep_eig, ppp, mspe, edges) ~ rank + k,
                 data=res, FUN=mean)
agg$rank <- factor(agg$rank, levels=c("Full","2","3"))
agg <- agg[order(agg$rank, agg$k), ]
cat(sprintf("\nobserved max eigenvalue = %.2f\n", mean(res$obs_eig)))
cat("\n--- LaTeX rows ---\n")
for (i in seq_len(nrow(agg))) with(agg[i,],
                                   cat(sprintf("%s & %d & %.3f & %.3f & %.2f & %.3f & %.3f & %.1f \\\\\n",
                                               rank, k, pr_expl, med_rho, rep_eig, ppp, mspe, edges)))
saveRDS(list(grid=res, agg=agg), file.path(OUT,"rank_table_p29.rds"))
cat("\nsaved rank_table_p29.rds\n")