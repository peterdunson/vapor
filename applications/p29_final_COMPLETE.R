# ============================================================================
#  p = 29 FRED-QD production rerun
#
#  Production rerun. Every setting lives in the CFG block below and is passed
#  explicitly to the fitting functions; nothing reads settings out of the
#  global environment, which is what made the earlier PrecFact rows
#  unreproducible.
#
#  RUN_LONG = FALSE re-runs only the production fits and loads pc_vapor, the
#  CI edges, R-hat/ESS, the posterior predictive check and the sensitivity
#  sweep from p29_full_results.rds. Set it to TRUE to refit the four long
#  chains as well (much slower).
#
#  Rows: VAPOR (r=3), VAR + Glasso at full rank and r=3, VAR + PrecFact at
#  full rank and r=3, PrecFact (no VAR), Kastner (VAR+FSV, full rank only --
#  bayesianVARs has no reduced-rank option), sample covariance, random walk.
#  The r=3 rows share one first-stage VAR: B_rr is the rank-3 SVD truncation
#  of the same B_hat the full-rank glasso row uses.
#
#  SEEDS: every stochastic row -- VAPOR, all three PrecFact rows and Kastner
#  -- uses the published VAPOR seed set SEED+1, SEED+2, SEED+3. The glasso,
#  sample-covariance and random-walk rows are deterministic; the first-stage
#  VAR they rest on is fit with SEED.
# ============================================================================
library(VAPOR); library(PrecFact); library(bayesianVARs); library(glasso)
library(coda); library(MASS)

# ============================================================================
#  CONFIG -- the only place settings are defined
# ============================================================================
CFG <- list(
   LABEL      = "p = 29 FRED-QD production rerun",
   RUN_LONG   = FALSE,          # TRUE also refits the four long chains

   SEED       = 32,
   K          = 10,            # precision factors; Kastner uses the same count
   RANK       = 3,
   R_REPS     = 3,

   NRUN       = 7500,  BURN   = 2500,  THIN   = 5,     # production MCMC
   NRUN_L     = 20000, BURN_L = 10000, THIN_L = 10,    # long chains, RUN_LONG only

   JIT        = 1e-3,          # vapor jitter
   VAR_JIT    = 1e-3,           # jitter for the bvar / Kastner fits
   #  100 held-out draws, with the per-draw seed recorded, at BOTH dimensions.
   #  Each draw hides a fresh 30% of the series in every test quarter. NOTE:
   #  the draws index random maskings of the same 20 test quarters, so they
   #  sharpen the loss estimate but do NOT add independent observations --
   #  the Diebold-Mariano test still has 20 units.
   N_HID      = 100, HID_SEED_BASE = 90000, HID_FRAC = 0.3,
   RHO_GRID   = 10^seq(-2, 0, length = 30),
   #  Production MCMC is matched to VAPOR at BOTH dimensions: 7500 iterations,
   #  2500 burn-in, thin 5, 1000 retained draws for every stochastic method.
   #  PrecFact: prec_est stores NRUN/THIN = 1500 draws, so burnout 500 leaves
   #  1000 and discards exactly BURN = 2500 iterations.
   #  Kastner: bvar takes burnin and draws separately, so 2500 + 5000 with
   #  thin 5 is 7500 iterations and 1000 retained.
   PF_BURNOUT = 500,
   KAST_DRAWS = 5000, KAST_BURNIN = 2500, KAST_THIN = 5,
   #  The first-stage VAR is NOT a competitor and is deliberately left at its
   #  original length: changing it would move B_hat, and with it the glasso
   #  and VAR + PrecFact residuals.
   VAR_DRAWS = 2000, VAR_BURNIN = 1000,
   SAMP_INV   = "solve",

   OUT        = "~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application",
   PREV_RDS   = "p29_full_results.rds",
   OUT_RDS    = "p29_results.rds",
   VERIFY     = NULL)

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

cat(sprintf("\n%s | RUN_LONG=%s | p=%d  T_train=%d  T/p=%.2f  T_test=%d\n",
            CFG$LABEL, CFG$RUN_LONG, p, T_train, T_train/p, T_test))
cat(sprintf("seed %d | k %d | rank %d | nrun %d burn %d thin %d | %d production fits\n",
            CFG$SEED, CFG$K, CFG$RANK, CFG$NRUN, CFG$BURN, CFG$THIN, CFG$R_REPS))

# ============================================================================
#  HELPERS.  Every function takes what it needs as an argument; nothing reads
#  the config or the data out of the global environment.
# ============================================================================
omega_one <- function(fit, s){ Om <- tcrossprod(fit$Lambda[,,s])
diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om }
omega_hat <- function(fit) posterior_mean_omega(fit$Lambda, fit$residual_prec)
pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }
fixB <- function(B, p){ if (nrow(B)!=p) B <- t(B)
if (ncol(B)!=p) B <- B[,seq_len(p),drop=FALSE]; B }
spec_rad <- function(B) max(abs(eigen(B, only.values=TRUE)$values))
chol_sig <- function(Om, p){
   Sig <- tryCatch(chol2inv(chol(Om)), error=function(e) MASS::ginv(Om))
   tryCatch(chol(Sig), error=function(e) chol(Sig + diag(1e-8,p))) }

#  the reporting rule: 95% credible interval of the partial correlation
#  excluding zero. Accepts one fit or a list of fits; draws are pooled.
ci_set <- function(fits){
   if (!is.null(fits$Lambda)) fits <- list(fits)
   stk <- unlist(lapply(fits, function(f)
      lapply(seq_len(dim(f$Lambda)[3]), function(s) pcor_from_omega(omega_one(f,s)))),
      recursive = FALSE)
   PCa <- simplify2array(stk)
   lo <- apply(PCa,c(1,2),quantile,0.025); hi <- apply(PCa,c(1,2),quantile,0.975)
   E <- (lo > 0 | hi < 0); diag(E) <- FALSE
   list(E = E, pm = apply(PCa,c(1,2),mean), n = sum(E)/2) }

pr_explosive <- function(fit, p, n = 400){
   nd <- dim(fit$B)[3]; idx <- unique(round(seq(1,nd,length.out=min(n,nd))))
   rs <- vapply(idx, function(s) spec_rad(fixB(fit$B[,,s], p)), numeric(1))
   c(pr_expl = mean(rs > 1), median_rho = median(rs)) }

count_edges <- function(PC, thr) sum(abs(PC[upper.tri(PC)]) > thr, na.rm=TRUE)

# ============================================================================
#  EVALUATION HARNESS.  One set of held-out draws, shared by every method.
# ============================================================================
#  one recorded seed per draw, rather than one stream for all draws, so any
#  single held-out set can be reproduced on its own
make_hid_sets <- function(base_seed, n_draws, p, T_test, frac){
   n_hid <- floor(frac*p)
   out <- lapply(seq_len(n_draws), function(b){
      set.seed(base_seed + b)
      lapply(seq_len(T_test), function(t) sort(sample(p, n_hid))) })
   attr(out, "seeds") <- base_seed + seq_len(n_draws)
   out }

predict_held_out <- function(Om, y_obs, obs, hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)

loss_one <- function(Om, hl, Y_test, p) vapply(seq_along(hl), function(t){
   hid <- hl[[t]]; obs <- setdiff(seq_len(p), hid); yt <- Y_test[t,]
   mean((yt[hid] - predict_held_out(Om, yt[obs], obs, hid))^2)}, numeric(1))

#  per-quarter loss vector, averaged over the draws of the held-out sets
loss_avg <- function(Om, hid_sets, Y_test, p)
   rowMeans(vapply(hid_sets, function(hl) loss_one(Om, hl, Y_test, p),
                   numeric(nrow(Y_test))))

loss_rw <- function(hid_sets, Y_test, Y_last_train, p)
   rowMeans(vapply(hid_sets, function(hl)
      vapply(seq_len(nrow(Y_test)), function(t){ hid <- hl[[t]]
      y_prev <- if (t==1) Y_last_train else Y_test[t-1,]
      mean((Y_test[t,hid] - y_prev[hid])^2)}, numeric(1)), numeric(nrow(Y_test))))

#  Diebold-Mariano with the Harvey-Leybourne-Newbold correction, h = 1
dm_test <- function(loss_comp, loss_vapor){
   d <- loss_comp - loss_vapor; n <- length(d)
   se <- sd(d)/sqrt(n); stat <- (mean(d)/se) * sqrt((n - 1)/n)
   data.frame(Difference = round(mean(d),3), SE = round(se,3), DM = round(stat,2),
              p = round(2*pt(-abs(stat), df=n-1),3),
              Won = sprintf("%d / %d", sum(d > 0), n)) }

# ============================================================================
#  FITTERS.  All settings explicit.
# ============================================================================
fit_vapor_one <- function(Y, k, rank, seed, nrun, burn, thin, jitter, max_try = 4){
   f <- NULL; jj <- jitter; a <- 0
   while (is.null(f) && a < max_try){ a <- a + 1
   f <- tryCatch(vapor(Y, lags=1, k=k, rank=rank, prior_B="dl",
                       nrun=nrun, burn=burn, thin=thin, seed=seed,
                       jitter=jj, init="random"),
                 error=function(e){ message("  vapor retry with larger jitter"); NULL })
   jj <- jj*2 }
   if (is.null(f)) stop(sprintf("vapor failed at seed %d", seed))
   attr(f, "attempts") <- a
   f }

#  VAR(1) posterior mean from bayesianVARs, plus the design matrices.
#  PHI is (p*lags + 1) x p with the intercept in the last row.
fit_var_mean <- function(Y, seed, draws, burnin, jit = 1e-3){
   set.seed(seed)
   Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, jit), nrow(Y), ncol(Y))
   m <- bvar(Yj, lags=1L, draws=draws, burnin=burnin, quiet=TRUE)
   B_hat <- apply(m$PHI, c(1,2), mean)
   list(B_hat = B_hat, X = m$X, Y = m$Y, Y_jit = Yj,
        resid = m$Y - m$X %*% B_hat) }

#  rank-r truncation of the lag block; the intercept row is left alone
rank_truncate_B <- function(B_hat, p, r){
   lag_rows <- seq_len(p)
   sv <- svd(B_hat[lag_rows,,drop=FALSE])
   B <- B_hat
   B[lag_rows,] <- sv$u[,seq_len(r),drop=FALSE] %*% diag(sv$d[seq_len(r)], r) %*%
      t(sv$v[,seq_len(r),drop=FALSE])
   list(B = B, sv = sv$d,
        frob_kept = sum(sv$d[seq_len(r)]^2)/sum(sv$d^2)) }

glasso_bic <- function(R, rho_grid, T_train){
   S <- cov(scale(R))
   bic_vals <- sapply(rho_grid, function(rho){
      gl <- tryCatch(glasso(S, rho=rho, penalize.diagonal=FALSE), error=function(e) NULL)
      if (is.null(gl)) return(Inf)
      Om <- gl$wi; ne <- sum(Om[upper.tri(Om)] != 0)
      ll <- (T_train/2)*(log(det(Om)) - sum(diag(S %*% Om))); -2*ll + log(T_train)*ne })
   rho <- rho_grid[which.min(bic_vals)]
   gl <- glasso(S, rho=rho, penalize.diagonal=FALSE)
   E <- gl$wi != 0; diag(E) <- FALSE
   list(Omega = gl$wi, rho = rho, nonzero = E, n_nonzero = as.integer(sum(E)/2)) }

fit_kastner_one <- function(Y, k, seed, draws, burnin, thin = 1L, jit = 1e-3){
   set.seed(seed)
   Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, jit), nrow(Y), ncol(Y))
   sp <- specify_prior_sigma(M=ncol(Y), type="factor", factor_factors=as.integer(k),
                             factor_priorfacloadtype="colwiseng",
                             factor_heteroskedastic=TRUE, quiet=TRUE)
   m <- tryCatch(bvar(Yj, lags=1L, draws=draws, burnin=burnin, thin=thin,
                      prior_sigma=sp, quiet=TRUE), error=function(e) NULL)
   if (is.null(m)) return(NULL)
   #  ONE RULE FOR EVERY BAYESIAN METHOD: the point estimate is the posterior
   #  mean of the per-draw PRECISION matrix. This previously returned
   #  solve(apply(vcov(m)[1,,,], c(1,2), mean)) -- the inverse of the mean
   #  COVARIANCE -- which made Kastner the only method on a different rule.
   #  The two differ by Jensen's inequality, materially so at p = 100.
   V <- vcov(m); nd <- dim(V)[4]; q <- ncol(Y)
   acc <- matrix(0, q, q)
   for (s in seq_len(nd)) acc <- acc + solve(V[1,,,s])
   acc/nd }

#  Settings are explicit arguments and are recorded on the result, rather than
#  read from whatever K / NRUN / THIN hold in the global environment when the
#  call runs. Reading them implicitly is how an interactive session can produce
#  numbers the script cannot reproduce later. The attempt count is recorded
#  too: the retry path reseeds AND jitters, so a silent retry changes the
#  answer for that rep.
run_pf <- function(Ymat, seed0, k, nrun, thin, burnout, alpha = 0.1, beta = 0.85,
                   epsilon_grid = seq(0.05,0.1,length=20), max_try = 6){
   out <- NULL; a <- 0
   while (is.null(out) && a < max_try){ a <- a + 1; set.seed(seed0 + a*7919)
   Yj <- if (a==1) Ymat else Ymat + matrix(rnorm(length(Ymat),0,1e-6), nrow(Ymat))
   out <- tryCatch({ res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=k,
                                     nrun=nrun, thin=thin, Y=Yj,
                                     alldiff=0, a_dir=0.1, b_dir=0.1)
   o <- fact_prec_processing(res, burnout=burnout, alpha=alpha, beta=beta,
                             epsilon_grid=epsilon_grid)
   diag(o$P) <- 0; o }, error=function(e) NULL) }
   if (is.null(out)) stop(sprintf("prec_est failed after %d attempts at seed %d", max_try, seed0))
   if (a > 1) warning(sprintf("run_pf seed %d: %d attempts; the retry path reseeds and jitters",
                              seed0, a))
   out$attempts <- a
   out$settings <- list(seed0=seed0, k=k, nrun=nrun, thin=thin, burnout=burnout,
                        alpha=alpha, beta=beta)
   out }

# ============================================================================
#  METHOD RECORDS.  One entry per row of the MSPE table.
# ============================================================================
make_record <- function(label, omegas, losses, extra = list()){
   list(label = label,
        omega = omegas,
        pcor  = if (length(omegas)) lapply(omegas, pcor_from_omega) else list(),
        loss  = losses,
        mspe  = vapply(losses, mean, numeric(1)),
        n_fits = length(losses),
        extra = extra) }

hid_sets <- make_hid_sets(CFG$HID_SEED_BASE, CFG$N_HID, p, T_test, CFG$HID_FRAC)
cat(sprintf("held-out evaluation: %d draws, seeds %d..%d, %d%% hidden per quarter\n",
            CFG$N_HID, min(attr(hid_sets,"seeds")), max(attr(hid_sets,"seeds")),
            round(100*CFG$HID_FRAC)))
M <- list()

# ---- first-stage VAR, shared by the glasso and PrecFact rows ---------------
cat("\n=== first-stage VAR(1) and its rank-", CFG$RANK, " truncation ===\n", sep="")
vr <- fit_var_mean(Y_train, seed = CFG$SEED, draws = CFG$VAR_DRAWS,
                   burnin = CFG$VAR_BURNIN, jit = CFG$VAR_JIT)
rt <- rank_truncate_B(vr$B_hat, p, CFG$RANK)
resid_full <- vr$resid
resid_rr   <- vr$Y - vr$X %*% rt$B
cat(sprintf("singular values of the lag block: %s ...\n",
            paste(sprintf("%.2f", head(rt$sv, 6)), collapse=", ")))
cat(sprintf("rank-%d keeps %.1f%% of the squared Frobenius norm\n",
            CFG$RANK, 100*rt$frob_kept))
cat(sprintf("spectral radius: full %.3f | rank %d %.3f\n",
            spec_rad(vr$B_hat[seq_len(p),]), CFG$RANK, spec_rad(rt$B[seq_len(p),])))

# ---- VAPOR, seeds SEED, SEED+1, SEED+2 -------------------------------------
cat(sprintf("\n=== VAPOR (rank %d, k %d, %d iterations) x %d ===\n",
            CFG$RANK, CFG$K, CFG$NRUN, CFG$R_REPS))
vapor_seeds <- CFG$SEED + seq_len(CFG$R_REPS)
vapor_fits <- vector("list", CFG$R_REPS); vapor_times <- numeric(CFG$R_REPS)
for (r in seq_len(CFG$R_REPS)){
   tt <- system.time(f <- fit_vapor_one(Y_train, k=CFG$K, rank=CFG$RANK,
                                        seed=vapor_seeds[r], nrun=CFG$NRUN,
                                        burn=CFG$BURN, thin=CFG$THIN, jitter=CFG$JIT))
   vapor_times[r] <- tt["elapsed"]; vapor_fits[[r]] <- f
   cat(sprintf("  seed %d: %.0fs\n", vapor_seeds[r], vapor_times[r])) }
vapor_om   <- lapply(vapor_fits, omega_hat)
vapor_loss <- lapply(vapor_om, loss_avg, hid_sets, Y_test, p)
vapor_ci   <- lapply(vapor_fits, function(f) ci_set(f)$n)
M[["VAPOR (rank 3)"]] <- make_record(
   sprintf("VAPOR (rank %d, k %d)", CFG$RANK, CFG$K), vapor_om, vapor_loss,
   extra = list(seeds = vapor_seeds, times = vapor_times,
                ci_edges_per_fit = unlist(vapor_ci),
                pr_explosive = pr_explosive(vapor_fits[[1]], p)))

# ---- VAR + Glasso, full rank and rank 3 ------------------------------------
cat("\n=== VAR + Glasso ===\n")
gl_full <- glasso_bic(resid_full, CFG$RHO_GRID, T_train)
gl_rr   <- glasso_bic(resid_rr,   CFG$RHO_GRID, T_train)
cat(sprintf("full rank: rho %.4f | %d nonzero off-diagonals\n", gl_full$rho, gl_full$n_nonzero))
cat(sprintf("rank %d   : rho %.4f | %d nonzero off-diagonals\n",
            CFG$RANK, gl_rr$rho, gl_rr$n_nonzero))
M[["VAR + Glasso (full rank)"]] <- make_record(
   "VAR + Glasso (full rank)", list(gl_full$Omega),
   list(loss_avg(gl_full$Omega, hid_sets, Y_test, p)),
   extra = list(rho = gl_full$rho, nonzero = gl_full$nonzero,
                n_nonzero = gl_full$n_nonzero))
#  "VAR + Glasso" now MEANS the rank-3 version. The full-rank fit is still
#  computed and saved (the paper reports it in one sentence) but is kept out
#  of every printed table via ORDER.
cond_gl_rr <- c(precision   = kappa(gl_rr$Omega, exact = TRUE),
                implied_cov = kappa(solve(gl_rr$Omega), exact = TRUE))
cat(sprintf("rank-%d glasso condition number: precision %.4g | implied covariance %.4g\n",
            CFG$RANK, cond_gl_rr["precision"], cond_gl_rr["implied_cov"]))
M[["VAR + Glasso"]] <- make_record(
   "VAR + Glasso", list(gl_rr$Omega),
   list(loss_avg(gl_rr$Omega, hid_sets, Y_test, p)),
   extra = list(rank = CFG$RANK, rho = gl_rr$rho, nonzero = gl_rr$nonzero,
                n_nonzero = gl_rr$n_nonzero, condition = cond_gl_rr))

# ---- PrecFact: no VAR, full-rank residuals, rank-3 residuals ---------------
cat(sprintf("\n=== PrecFact (k %d, burnout %d) ===\n", CFG$K, CFG$PF_BURNOUT))
pf_run <- function(Ymat, seeds){
   lapply(seeds, function(s) run_pf(Ymat, seed0 = s, k = CFG$K, nrun = CFG$NRUN,
                                    thin = CFG$THIN, burnout = CFG$PF_BURNOUT)) }
#  one seed set for every stochastic method, matching the published VAPOR runs
seeds_pf <- seeds_pfv <- seeds_pfv3 <- CFG$SEED + seq_len(CFG$R_REPS)
pf_raw  <- pf_run(Y_train,           seeds_pf)
pf_full <- pf_run(scale(resid_full), seeds_pfv)
pf_rr   <- pf_run(scale(resid_rr),   seeds_pfv3)
cat(sprintf("attempts -- no VAR: %s | full rank: %s | rank %d: %s\n",
            paste(vapply(pf_raw, function(z) z$attempts, numeric(1)), collapse=","),
            paste(vapply(pf_full, function(z) z$attempts, numeric(1)), collapse=","),
            CFG$RANK,
            paste(vapply(pf_rr, function(z) z$attempts, numeric(1)), collapse=",")))
pf_rec <- function(lab, fits, seeds) make_record(
   lab, lapply(fits, function(z) z$Omega),
   lapply(fits, function(z) loss_avg(z$Omega, hid_sets, Y_test, p)),
   extra = list(seeds = seeds,
                pfdr_graph = lapply(fits, function(z) z$P),
                pfdr_edges = vapply(fits, function(z) sum(z$P)/2, numeric(1)),
                attempts = vapply(fits, function(z) z$attempts, numeric(1)),
                settings = lapply(fits, function(z) z$settings)))
M[["PrecFact (no VAR)"]]        <- pf_rec("PrecFact (no VAR)", pf_raw, seeds_pf)
M[["VAR + PrecFact (full rank)"]] <- pf_rec("VAR + PrecFact (full rank)", pf_full, seeds_pfv)
#  likewise, "VAR + PrecFact" now MEANS the rank-3 version
M[["VAR + PrecFact"]] <- pf_rec("VAR + PrecFact", pf_rr, seeds_pfv3)

# ---- Kastner (VAR + FSV), full rank; the package has no reduced-rank option -
cat(sprintf("\n=== Kastner (VAR+FSV), %d factors ===\n", CFG$K))
kast_seeds <- CFG$SEED + seq_len(CFG$R_REPS)
kast_fits <- Filter(Negate(is.null),
                    lapply(kast_seeds, function(s)
                       fit_kastner_one(Y_train, k=CFG$K, seed=s,
                                       draws=CFG$KAST_DRAWS, burnin=CFG$KAST_BURNIN,
                                       thin=CFG$KAST_THIN, jit=CFG$VAR_JIT)))
if (!length(kast_fits)) stop("all Kastner fits failed")
M[["Kastner (VAR+FSV)"]] <- make_record(
   "Kastner (VAR+FSV)", kast_fits,
   lapply(kast_fits, loss_avg, hid_sets, Y_test, p),
   extra = list(seeds = kast_seeds, n_fits = length(kast_fits),
                condition = vapply(kast_fits, kappa, numeric(1), exact=TRUE)))

cat(sprintf("Kastner per run -- |pcor| > 0.05: %s\n",
            paste(vapply(kast_fits, function(O){
               pc <- pcor_from_omega(O); sum(abs(pc[upper.tri(pc)]) > 0.05) }, numeric(1)),
                  collapse = ", ")))
cat(sprintf("Kastner condition numbers: %s\n",
            paste(sprintf("%.4g", vapply(kast_fits, kappa, numeric(1), exact = TRUE)),
                  collapse = ", ")))

# ---- sample covariance and random walk -------------------------------------
Omega_samp <- if (identical(CFG$SAMP_INV, "ginv")) MASS::ginv(cov(Y_train)) else
   solve(cov(Y_train))
M[["Sample covariance"]] <- make_record(
   "Sample covariance", list(Omega_samp),
   list(loss_avg(Omega_samp, hid_sets, Y_test, p)),
   extra = list(inverse = CFG$SAMP_INV, condition = kappa(Omega_samp, exact=TRUE)))
M[["Random walk"]] <- make_record(
   "Random walk", list(), list(loss_rw(hid_sets, Y_test, Y_last_train, p)))

# ---- cross-check the full-rank rows against previously verified values -----
#  These two rows are not new work: they were verified to reproduce earlier.
#  They are refit here so the script stays self-contained, and checked against
#  the stored loss vectors so any drift is caught rather than assumed away.
if (length(CFG$VERIFY)){
   cat("\n=== cross-check against previously verified values ===\n")
   for (v in CFG$VERIFY){
      f <- file.path(path.expand(CFG$OUT), v$file)
      if (!file.exists(f)){
         cat(sprintf("  %-28s SKIPPED, %s not found\n", v$method, v$file)); next }
      ref  <- Reduce(function(a,b) a[[b]], v$path, readRDS(f))
      ours <- M[[v$method]]$loss
      refl <- if (is.list(ref)) ref else list(ref)
      ok <- length(refl) == length(ours) && all(mapply(function(a,b)
         isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-10)), refl, ours))
      cat(sprintf("  %-28s reproduces: %-3s (%s)\n", v$method,
                  if (ok) "YES" else "NO", v$file)) } }

# ============================================================================
#  (a) MSPE TABLE, paper order
# ============================================================================
#  Paper order. "VAR + Glasso" and "VAR + PrecFact" are the RANK-3 fits. The
#  full-rank fits live in M and in the saved object under glasso_fullrank and
#  precfact_fullrank, and appear in no printed table.
ORDER <- c("PrecFact (no VAR)", "VAR + Glasso", "VAR + PrecFact",
           "Kastner (VAR+FSV)", "Sample covariance", "Random walk", "VAPOR (rank 3)")
fmt_mspe <- function(v) if (length(v) > 1) sprintf("%.4f (%.4f)", mean(v), sd(v)) else
   sprintf("%.4f", mean(v))
mspe_tab <- data.frame(
   Method = vapply(ORDER, function(k) M[[k]]$label, character(1)),
   MSPE   = vapply(ORDER, function(k) fmt_mspe(M[[k]]$mspe), character(1)),
   Fits   = vapply(ORDER, function(k) M[[k]]$n_fits, integer(1)),
   row.names = NULL)
cat(sprintf("\n=== (a) MSPE over %d test quarters, averaged over %d held-out draws ===\n",
            T_test, CFG$N_HID))
print(mspe_tab, row.names = FALSE)

seed_tab <- do.call(rbind, lapply(ORDER, function(k){
   sd_k <- M[[k]]$extra$seeds
   data.frame(Method = M[[k]]$label,
              seeds  = if (is.null(sd_k)) "deterministic" else paste(sd_k, collapse=", ")) }))
cat("\n=== seeds per row (first-stage VAR fit with seed ", CFG$SEED, ") ===\n", sep="")
print(seed_tab, row.names = FALSE)

# ============================================================================
#  (b) DIEBOLD-MARIANO vs VAPOR, HLN corrected
# ============================================================================
mean_loss <- function(k) Reduce(`+`, M[[k]]$loss)/length(M[[k]]$loss)
lv <- mean_loss("VAPOR (rank 3)")
dm_tab <- do.call(rbind, lapply(setdiff(ORDER, "VAPOR (rank 3)"), function(k)
   cbind(Method = M[[k]]$label, MSPE = round(mean(mean_loss(k)),3),
         dm_test(mean_loss(k), lv))))
cat("\n=== (b) Diebold-Mariano vs VAPOR, HLN corrected (positive favours VAPOR) ===\n")
print(dm_tab, row.names = FALSE)

# ============================================================================
#  LONG-RUN SECTIONS: run them, or load them from the previous results file
# ============================================================================
prev_path <- file.path(path.expand(CFG$OUT), CFG$PREV_RDS)
if (!CFG$RUN_LONG){
   if (!file.exists(prev_path))
      stop(sprintf("RUN_LONG is FALSE but %s is missing", prev_path))
   prev <- readRDS(prev_path)
   long <- list(source = CFG$PREV_RDS, rerun = FALSE,
                pc_vapor = prev$pc_vapor, ci_edges = prev$ci_edges,
                edge_ci = prev$edge_ci, diag_tab = prev$diag_tab,
                ppc = prev$ppc, sens = prev$sens, stability = prev$stability,
                edge_rhat = prev$edge_rhat, edge_ranks = prev$edge_ranks,
                justification = prev$justification)
   cat(sprintf("\n[RUN_LONG FALSE] pc_vapor, CI edges (%d), R-hat/ESS, PPC and sensitivity\n",
               nrow(prev$ci_edges)))
   cat(sprintf("                 loaded from %s\n", CFG$PREV_RDS))
} else {
   cat(sprintf("\n=== long chains: 4 x %d, burn %d, thin %d ===\n",
               CFG$NRUN_L, CFG$BURN_L, CFG$THIN_L))
   chains <- Filter(Negate(is.null), lapply(1:4, function(s){
      cat(sprintf("\r  chain %d/4", s)); flush.console()
      fit_vapor_one(Y_train, k=CFG$K, rank=CFG$RANK, seed=100+s,
                    nrun=CFG$NRUN_L, burn=CFG$BURN_L, thin=CFG$THIN_L,
                    jitter=CFG$JIT) }))
   cat("\n")
   cs <- ci_set(chains)
   idx <- which(cs$E & upper.tri(cs$E), arr.ind=TRUE)
   ci_edges <- data.frame(v1=vars_qd[idx[,1]], v2=vars_qd[idx[,2]],
                          pcor=round(cs$pm[idx],3))
   ci_edges <- ci_edges[order(-abs(ci_edges$pcor)), ]
   pc_arrays <- lapply(chains, function(f){
      nd <- dim(f$Lambda)[3]
      simplify2array(lapply(seq_len(nd), function(s) pcor_from_omega(omega_one(f,s)))) })
   series <- list(loglik = lapply(chains, function(f) colSums(f$loglik)))
   for (e in seq_len(nrow(ci_edges))){
      i <- match(ci_edges$v1[e], vars_qd); j <- match(ci_edges$v2[e], vars_qd)
      series[[paste(ci_edges$v1[e], ci_edges$v2[e], sep="_")]] <-
         lapply(pc_arrays, function(A) A[i,j,]) }
   series[["rho_B"]] <- lapply(chains, function(f)
      vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s], p)), numeric(1)))
   diag_of <- function(x){ mc <- mcmc.list(lapply(x, mcmc))
   c(ESS = sum(effectiveSize(mc)),
     Rhat = tryCatch(gelman.diag(mc, autoburnin=FALSE)$psrf[1], error=function(e) NA)) }
   da <- do.call(rbind, lapply(series, diag_of))
   diag_tab <- data.frame(quantity=rownames(da), ESS=round(da[,"ESS"],0),
                          Rhat=round(da[,"Rhat"],3), row.names=NULL)
   print(diag_tab, row.names=FALSE)
   long <- list(source = "this run", rerun = TRUE, pc_vapor = cs$pm,
                ci_edges = ci_edges, edge_ci = cs$E, diag_tab = diag_tab,
                chains_kept = length(chains))
}

# ============================================================================
#  SAVE
# ============================================================================
out <- list(
   config     = CFG,
   provenance = list(R = R.version.string,
                     packages = vapply(c("VAPOR","PrecFact","bayesianVARs","glasso","coda"),
                                       function(pk) as.character(packageVersion(pk)),
                                       character(1)),
                     run_at = format(Sys.time()),
                     sessionInfo = utils::capture.output(utils::sessionInfo())),
   methods    = M,
   order      = ORDER,
   #  reported in one sentence of the paper; deliberately absent from ORDER
   glasso_fullrank   = M[["VAR + Glasso (full rank)"]],
   precfact_fullrank = M[["VAR + PrecFact (full rank)"]],
   mspe_table = mspe_tab, seed_table = seed_tab,
   dm_table   = dm_tab,
   long       = long,
   first_stage = list(B_hat = vr$B_hat, B_rr = rt$B, singular_values = rt$sv,
                      frob_kept = rt$frob_kept,
                      resid_full = resid_full, resid_rr = resid_rr),
   hid_sets = hid_sets, hid_seeds = attr(hid_sets, "seeds"),
   Y_test = Y_test, Y_train = Y_train,
   vars = vars_qd, p = p, T_train = T_train, T_test = T_test)
saveRDS(out, file.path(path.expand(CFG$OUT), CFG$OUT_RDS))
cat(sprintf("\nSaved %s\n", CFG$OUT_RDS))
