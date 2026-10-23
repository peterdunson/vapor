# ============================================================================
#  SHARED HELPERS for the final FRED-QD application.
#
#  Sourced by final_p29.R, final_p100.R, final_sweeps.R, final_ranks_p29.R.
#  Nothing here reads the configuration or the data out of the global
#  environment: every function takes what it needs as an argument.
#
#  Two things differ from the superseded p*_final_COMPLETE.R scripts:
#    1. Centering and scaling use the TRAINING-window mean and sd only (all
#       training quarters), applied to the training and the test rows.  The old
#       scripts called scale() on the stacked train+test matrix.
#    2. Every stochastic method runs 4 chains at the same length, and a chain
#       that errors is replaced by the next seed, at most 3 replacements per
#       method and dimension.  The idiosyncratic-precision floor is now in the
#       installed VAPOR and PrecFact samplers, so errors should be rare; the
#       ladder records them when they happen rather than aborting the run.
#
#  The point estimate for EVERY Bayesian method is the posterior mean of the
#  per-draw precision matrix.
# ============================================================================
suppressMessages({
   library(VAPOR); library(PrecFact); library(bayesianVARs)
   library(glasso); library(posterior); library(MASS)
})

# ============================================================================
#  OUTPUT LOCATIONS.  Everything this study writes lives inside the repository,
#  under applications/results/.  Nothing is written to the precision_factor
#  tree any more; the raw FRED-QD csv is the only thing read from outside.
# ============================================================================
RES_DIR  <- "/Users/peterdunson/Desktop/vapor/applications/results"
FIG_DIR  <- file.path(RES_DIR, "figures")
TAB_DIR  <- file.path(RES_DIR, "tables")
FITS_DIR <- file.path(RES_DIR, "fits")
LOG_DIR  <- file.path(RES_DIR, "logs")
#  Simulation code and simulation output live under simulation/; the tree above
#  holds only the two empirical applications.  The one file both halves write is
#  results_summary.md, which stays in TAB_DIR because it is a single combined
#  summary of the paper.
SIM_DIR   <- "/Users/peterdunson/Desktop/vapor/simulation"
SIM_RES   <- file.path(SIM_DIR, "results")
SIM_FIG   <- file.path(SIM_RES, "figures")
SIM_TAB   <- file.path(SIM_RES, "tables")
SIM_FITS  <- file.path(SIM_RES, "fits")
SIM_LOG   <- file.path(SIM_RES, "logs")
SIM_PARTS <- file.path(SIM_RES, "sim_f_parts")
for (d in c(FIG_DIR, TAB_DIR, FITS_DIR, LOG_DIR,
            SIM_FIG, SIM_TAB, SIM_FITS, SIM_LOG, SIM_PARTS))
   if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ============================================================================
#  DATA
# ============================================================================
#  Same series, transformations and windows as the production scripts.  The one
#  change is the standardization window, see the header.
build_panel <- function(CFG){
   raw <- read.csv(path.expand(CFG$CSV))
   tf <- function(x, tcode){ x <- as.numeric(x)
      switch(as.character(tcode), "1"=x, "2"=c(NA,diff(x)), "3"=c(NA,NA,diff(diff(x))),
             "4"=log(x), "5"=c(NA,diff(log(x))), "6"=c(NA,NA,diff(diff(log(x)))),
             "7"=c(NA,diff(log(x)/lag(log(x))-1)), x) }
   tc  <- as.numeric(raw[2,-1]); names(tc) <- colnames(raw)[-1]
   dat <- raw[-c(1,2),]; dts <- dat[,1]; dat <- apply(dat[,-1], 2, as.numeric)
   TRm <- as.matrix(suppressWarnings(mapply(tf, as.data.frame(dat), tc)))
   rownames(TRm) <- dts
   pd <- function(s){ z <- as.numeric(strsplit(s,"/")[[1]]); z[3] + (z[1]-1)/12 }

   #  ONE branch for both dimensions: the series are selected first and the NA
   #  row filter is applied to the SELECTED series only.  Previously p = 29 ran
   #  the filter over the whole 245-column matrix, which dropped every quarter
   #  in which any FRED-QD series was missing even when all 29 were present:
   #  that cost 2008Q1-2009Q1 (NONBORRES) and every quarter before 2000Q2
   #  (SPCS20RSA), and left 2007Q4 paired with 2009Q2 in the lagged design.
   vars <- CFG$VARS[CFG$VARS %in% colnames(TRm)]
   dv   <- vapply(rownames(TRm), pd, numeric(1))
   Ysub <- TRm[, vars]
   keep <- which(dv >= CFG$WIN[1] & dv <= CFG$WIN[2]); Ysub <- Ysub[keep,,drop=FALSE]
   ok   <- apply(Ysub,1,function(r) !any(is.na(r))); Ysub <- Ysub[ok,,drop=FALSE]
   dvk  <- dv[keep][ok]
   tr <- which(dvk >= CFG$TRAIN[1] & dvk <= CFG$TRAIN[2])
   te <- which(dvk >= CFG$TEST[1]  & dvk <= CFG$TEST[2])

   #  training-window moments, applied to every row
   mu <- colMeans(Ysub[tr,,drop=FALSE]); sd <- apply(Ysub[tr,,drop=FALSE], 2, sd)
   Z  <- sweep(sweep(Ysub, 2, mu, "-"), 2, sd, "/")
   qlab <- function(s){ z <- as.numeric(strsplit(s,"/")[[1]])
      sprintf("%.0fQ%.0f", z[3], ceiling(z[1]/3)) }
   out <- list(train = Z[tr,,drop=FALSE], test = Z[te,,drop=FALSE],
               vars = vars, p = length(vars), mu = mu, sd = sd,
               raw_train = Ysub[tr,,drop=FALSE])
   cat(sprintf("p = %d | standardization: training quarters only (%s..%s, %d rows)\n",
               out$p, qlab(rownames(out$train)[1]),
               qlab(rownames(out$train)[nrow(out$train)]), nrow(out$train)))
   cat(sprintf("train %d x %d | test %d rows (%s..%s)\n", nrow(out$train), out$p,
               nrow(out$test), qlab(rownames(out$test)[1]),
               qlab(rownames(out$test)[nrow(out$test)])))
   gaps <- function(M) sort(unique(round(diff(vapply(rownames(M), pd, numeric(1)))*12)))
   out$gaps <- list(train = gaps(out$train), test = gaps(out$test))
   cat(sprintf("gaps in months: train %s | test %s | all consecutive quarters: %s\n",
               paste(out$gaps$train, collapse=","), paste(out$gaps$test, collapse=","),
               all(c(out$gaps$train, out$gaps$test) == 3)))
   out
}

#  a cheap content fingerprint, used to tell whether a cached fit belongs to
#  the matrix in hand; not a cryptographic hash, just enough to catch a
#  different panel
digest_mat <- function(M) c(dim = dim(M), sum = sum(M), ssq = sum(M^2),
                            first = unname(M[1,1]), last = unname(M[nrow(M), ncol(M)]))

# ============================================================================
#  PRECISION, PARTIAL CORRELATION, EDGES
# ============================================================================
omega_one <- function(fit, s){ Om <- tcrossprod(fit$Lambda[,,s])
   diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om }
#  convention (a): posterior mean of the per-draw precision
omega_hat <- function(fit){ nd <- dim(fit$Lambda)[3]; p <- nrow(fit$Lambda)
   acc <- matrix(0,p,p); for (s in seq_len(nd)) acc <- acc + omega_one(fit,s); acc/nd }
pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }
fixB <- function(B, p){ if (nrow(B)!=p) B <- t(B)
   if (ncol(B)!=p) B <- B[,seq_len(p),drop=FALSE]; B }
spec_rad <- function(B) max(abs(eigen(B, only.values=TRUE)$values))

#  every pair's partial-correlation draws for one fit, as a draws x pairs matrix
pair_draws <- function(fit, p){
   ut <- which(upper.tri(matrix(0,p,p)), arr.ind=TRUE)
   nd <- dim(fit$Lambda)[3]; M <- matrix(NA_real_, nd, nrow(ut))
   for (s in seq_len(nd)) M[s,] <- pcor_from_omega(omega_one(fit,s))[cbind(ut[,1],ut[,2])]
   attr(M, "ut") <- ut; M }
pair_labels <- function(vars){
   p <- length(vars); ut <- which(upper.tri(matrix(0,p,p)), arr.ind=TRUE)
   paste(vars[pmin(ut[,1],ut[,2])], vars[pmax(ut[,1],ut[,2])], sep="--") }
#  the reporting rule: pooled 95% credible interval of the partial correlation
#  excluding zero
ci_indices <- function(M){ lo <- apply(M,2,quantile,0.025); hi <- apply(M,2,quantile,0.975)
   which(lo > 0 | hi < 0) }

# ============================================================================
#  HELD-OUT EVALUATION.  One set of held-out draws, shared by every method.
#  One recorded seed per set, so any single set can be reproduced on its own.
# ============================================================================
make_hid_sets <- function(base_seed, n_draws, p, T_test, frac){
   n_hid <- floor(frac*p)
   out <- lapply(seq_len(n_draws), function(b){ set.seed(base_seed + b)
      lapply(seq_len(T_test), function(t) sort(sample(p, n_hid))) })
   attr(out, "seeds") <- base_seed + seq_len(n_draws); out }
predict_held_out <- function(Om, y_obs, obs, hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)
loss_one <- function(Om, hl, Y_test, p) vapply(seq_along(hl), function(t){
   hid <- hl[[t]]; obs <- setdiff(seq_len(p), hid); yt <- Y_test[t,]
   mean((yt[hid] - predict_held_out(Om, yt[obs], obs, hid))^2)}, numeric(1))
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
#  CHAIN RUNNERS.  Each returns list(chains, log, replacements), where the log
#  records every seed that errored and its message.  Seeds advance past
#  SEEDS[4] on failure, at most MAX_REPL times per method and dimension.
# ============================================================================
run_with_ladder <- function(fit_fun, seeds, max_repl, label){
   queue <- seeds; done <- list(); log <- list(); nrepl <- 0L
   nxt <- max(seeds) + 1L
   while (length(queue) && length(done) < length(seeds)){
      s <- queue[1]; queue <- queue[-1]
      t0 <- Sys.time()
      f <- tryCatch(fit_fun(s),
                    error = function(e) structure(list(msg=conditionMessage(e)),
                                                  class="ladder_err"))
      el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
      if (inherits(f, "ladder_err")){
         cat(sprintf("  %s seed %d FAILED after %.1f min: %s\n", label, s, el, f$msg))
         log[[length(log)+1L]] <- list(method=label, seed=s, status="error",
                                       error=f$msg, minutes=el)
         if (nrepl < max_repl){ nrepl <- nrepl + 1L; queue <- c(queue, nxt)
            cat(sprintf("  %s replacement %d of %d: queueing seed %d\n",
                        label, nrepl, max_repl, nxt)); nxt <- nxt + 1L
         } else cat(sprintf("  %s replacement budget exhausted\n", label))
         next }
      attr(f, "seed") <- s; attr(f, "minutes") <- el
      done[[length(done)+1L]] <- f
      log[[length(log)+1L]] <- list(method=label, seed=s, status="completed",
                                    minutes=el)
      cat(sprintf("  %s seed %d: %.1f min\n", label, s, el)); flush.console() }
   list(chains = done, log = log, replacements = nrepl) }

# ----------------------------------------------------------------------------
#  PARALLEL CHAIN RUNNER.  Every (method, seed) pair is one job; jobs run
#  concurrently on min(physical cores - 1, WORKERS) workers, across methods as
#  well as seeds.  Each job sets its own seed inside the sampler, so a forked
#  worker produces exactly the chain the sequential runner would.
#
#  Failures are handled in waves: after a wave, every method that lost a chain
#  and still has replacement budget queues its next seed, and the next wave
#  runs those concurrently.  The ladder's accounting is therefore per method,
#  as before, but the waiting is shared.
#
#  MEMORY.  At p = 100 and 1,250 retained draws a VAPOR chain is about 130 MB
#  and a Kastner chain about 100 MB, so 7 concurrent chains plus the collected
#  results run to roughly 3 GB.  Lower CFG$WORKERS if the machine is smaller.
# ----------------------------------------------------------------------------
n_workers <- function(cap = NULL){
   pc <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
   if (is.na(pc)) pc <- parallel::detectCores()
   w <- max(1L, as.integer(pc) - 1L)
   if (!is.null(cap)) w <- min(w, as.integer(cap))
   w }

#  max_repl may be a single number or a named vector, one entry per method
run_methods_parallel <- function(method_funs, seeds, max_repl, workers){
   mn <- names(method_funs)
   budget <- if (length(max_repl) == 1L) setNames(rep(as.integer(max_repl), length(mn)), mn)
             else { b <- setNames(rep(0L, length(mn)), mn)
                    b[names(max_repl)] <- as.integer(max_repl); b }
   chains <- setNames(vector("list", length(mn)), mn)
   logs   <- list()
   nrepl  <- setNames(integer(length(mn)), mn)
   nxt    <- setNames(rep(max(seeds) + 1L, length(mn)), mn)
   jobs <- unlist(lapply(mn, function(m)
      lapply(seeds, function(s) list(method = m, seed = s))), recursive = FALSE)
   wave <- 0L
   while (length(jobs)){
      wave <- wave + 1L
      cat(sprintf("  wave %d: %d chain(s) on %d worker(s)\n", wave, length(jobs), workers))
      flush.console()
      res <- parallel::mclapply(jobs, function(j){
         t0 <- Sys.time()
         out <- tryCatch(method_funs[[j$method]](j$seed),
                         error = function(e) structure(list(msg = conditionMessage(e)),
                                                       class = "ladder_err"))
         list(job = j, out = out,
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
      }, mc.cores = workers, mc.preschedule = FALSE)
      jobs <- list()
      for (i in seq_along(res)){
         r <- res[[i]]
         #  a worker that died takes its result with it; treat it as an error
         if (inherits(r, "try-error") || is.null(r$out) || inherits(r$out, "ladder_err")){
            m <- if (is.null(r$job)) jobs_method_fallback(res, i) else r$job$method
            sd <- if (is.null(r$job)) NA_integer_ else r$job$seed
            msg <- if (inherits(r, "try-error")) as.character(r)
                   else if (is.null(r$out)) "worker returned NULL"
                   else r$out$msg
            cat(sprintf("  %s seed %s FAILED: %s\n", m, sd, msg))
            logs[[length(logs)+1L]] <- list(method = m, seed = sd, status = "error",
               error = msg, minutes = if (is.null(r$minutes)) NA_real_ else r$minutes)
            if (nrepl[[m]] < budget[[m]]){
               nrepl[[m]] <- nrepl[[m]] + 1L
               cat(sprintf("  %s replacement %d of %d: queueing seed %d\n",
                           m, nrepl[[m]], budget[[m]], nxt[[m]]))
               jobs[[length(jobs)+1L]] <- list(method = m, seed = nxt[[m]])
               nxt[[m]] <- nxt[[m]] + 1L
            } else cat(sprintf("  %s replacement budget exhausted\n", m))
            next }
         f <- r$out
         attr(f, "seed") <- r$job$seed; attr(f, "minutes") <- r$minutes
         chains[[r$job$method]] <- c(chains[[r$job$method]], list(f))
         logs[[length(logs)+1L]] <- list(method = r$job$method, seed = r$job$seed,
                                         status = "completed", minutes = r$minutes)
         cat(sprintf("  %s seed %d: %.1f min\n", r$job$method, r$job$seed, r$minutes)) }
      flush.console() }
   #  keep each method's chains in seed order
   chains <- lapply(chains, function(l)
      if (length(l)) l[order(vapply(l, function(f) attr(f,"seed"), numeric(1)))] else l)
   list(chains = chains, log = logs, replacements = nrepl) }
jobs_method_fallback <- function(res, i) "unknown"

# ----------------------------------------------------------------------------
#  Per-method closures, for the parallel runner.  Each takes a seed and returns
#  the fit; settings come from CFG, the data from the arguments.
# ----------------------------------------------------------------------------
#  resid_full is the full-rank bayesianVARs posterior-mean first stage, which
#  the MAIN VAR + PrecFact and VAR + Glasso rows now use; resid_rr is the
#  rank-3 truncation, kept for the supplement rows.  Kastner is reused, not
#  refitted, so it is not in this list any more.
method_closures <- function(Y_train, resid_full, resid_rr, CFG){
   list(
      VAPOR = function(s) vapor(Y_train, lags=1, k=CFG$K, rank=CFG$RANK,
                                prior_B="dl", nrun=CFG$NRUN, burn=CFG$BURN,
                                thin=CFG$THIN, seed=s, jitter=CFG$JIT, init="random"),
      PrecFact = function(s) pf_one(Y_train, CFG, s),
      `VAR+PrecFact` = function(s) pf_one(scale(resid_full), CFG, s),
      `VAR+PrecFact (rank 3)` = function(s) pf_one(scale(resid_rr), CFG, s)) }

#  PrecFact: prec_est takes no burn argument, so the burn-in is applied as
#  burnout in units of RETAINED draws, burnout = BURN / THIN.
pf_one <- function(Ymat, CFG, seed){
   set.seed(seed)
   #  DP base measure Gamma(1, 1), matching VAPOR's new package default.
   res <- prec_est(as=1.0, bs=1.0, a=0.5, b=0.5, k=CFG$K,
                   nrun=CFG$NRUN, thin=CFG$THIN, Y=Ymat,
                   alldiff=0, a_dir=0.1, b_dir=0.1)
   o <- fact_prec_processing(res, burnout=CFG$BURN/CFG$THIN, alpha=0.1,
                             beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
   diag(o$P) <- 0
   keep <- seq.int(CFG$BURN/CFG$THIN + 1L, dim(res$Lambda)[3])
   list(Omega = o$Omega, P = o$P,
        Lambda = res$Lambda[,,keep,drop=FALSE],
        residual_prec = res$residuals[keep,,drop=FALSE]) }

#  Kastner: matching draws / burn-in / thin, point estimate = mean of the
#  per-draw precision (NOT the inverse of the mean covariance)
kastner_one <- function(Y, CFG, seed){
   set.seed(seed)
   Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, CFG$VAR_JIT), nrow(Y), ncol(Y))
   sp <- specify_prior_sigma(M=ncol(Y), type="factor",
                             factor_factors=as.integer(CFG$K),
                             factor_priorfacloadtype="colwiseng",
                             factor_heteroskedastic=TRUE, quiet=TRUE)
   #  item 5: bayesianVARs::bvar counts `draws` POST burn-in, so draws = NRUN
   #  used to run NRUN + BURN = 120,000 total iterations and retain NRUN/THIN =
   #  1250 draws.  Passing draws = NRUN - BURN makes the TOTAL NRUN = 100,000
   #  with BURN = 20,000 discarded and (NRUN - BURN)/THIN = 1000 retained,
   #  matching VAPOR, PrecFact and VAR + PrecFact exactly.
   k_draws <- as.integer(CFG$NRUN - CFG$BURN)
   m <- bvar(Yj, lags=1L, draws=k_draws, burnin=as.integer(CFG$BURN),
             thin=as.integer(CFG$THIN), prior_sigma=sp, quiet=TRUE)
   V <- vcov(m); nd <- dim(V)[4]; q <- ncol(Y)
   acc <- matrix(0, q, q); prec <- array(NA_real_, c(q,q,nd))
   for (i in seq_len(nd)){ Pi <- solve(V[1,,,i]); prec[,,i] <- Pi; acc <- acc + Pi }
   list(Omega = acc/nd, prec_draws = prec,
        settings = list(total_iterations = CFG$NRUN, burnin = CFG$BURN,
                        thin = CFG$THIN, draws_arg = k_draws, retained = nd,
                        factor_factors = as.integer(CFG$K))) }

# ----------------------------------------------------------------------------
#  TABLE S5 QUANTITIES, verbatim from the superseded ranks_p29.R
# ----------------------------------------------------------------------------
#  Pr(spectral radius of B > 1) and its median, over up to 400 draws
pr_explosive <- function(f, p, n_draws = 400L){
   nd <- dim(f$B)[3]; ix <- unique(round(seq(1, nd, length.out = min(n_draws, nd))))
   rs <- vapply(ix, function(s) spec_rad(fixB(f$B[,,s], p)), numeric(1))
   c(pr = mean(rs > 1), med = median(rs)) }

#  posterior predictive check on the largest eigenvalue of the correlation
#  matrix of a replicated panel
ppc_maxeig <- function(f, Y, p, n_draws = 200L, seed = 1L){
   maxeig <- function(M) max(eigen(cor(M), symmetric = TRUE, only.values = TRUE)$values)
   nd <- dim(f$Lambda)[3]; ix <- unique(round(seq(1, nd, length.out = min(n_draws, nd))))
   Tn <- nrow(Y); Yl <- Y[-Tn, , drop = FALSE]; obs <- maxeig(Y[-1, , drop = FALSE])
   set.seed(seed)
   r <- vapply(ix, function(s){
      Om <- omega_one(f, s)
      Sg <- tryCatch(chol2inv(chol(Om)), error = function(e) MASS::ginv(Om))
      L  <- tryCatch(chol(Sg), error = function(e) chol(Sg + diag(1e-8, p)))
      maxeig(Yl %*% t(fixB(f$B[,,s], p)) +
             matrix(rnorm((Tn-1)*p), Tn-1, p) %*% L) }, numeric(1))
   c(rep = mean(r), ppp = mean(r >= obs), obs = obs) }

#  Sequential fallbacks, kept so a single method can be re-run on its own.
#  They call exactly the one-chain functions the parallel runner uses, so the
#  two paths cannot drift apart.
vapor_chains <- function(Y, CFG)
   run_with_ladder(function(s) vapor(Y, lags=1, k=CFG$K, rank=CFG$RANK,
                                     prior_B="dl", nrun=CFG$NRUN, burn=CFG$BURN,
                                     thin=CFG$THIN, seed=s, jitter=CFG$JIT,
                                     init="random"),
                   CFG$SEEDS, CFG$MAX_REPL, "VAPOR")
pf_chains <- function(Ymat, CFG, label)
   run_with_ladder(function(s) pf_one(Ymat, CFG, s), CFG$SEEDS, CFG$MAX_REPL, label)
kastner_chains <- function(Y, CFG)
   run_with_ladder(function(s) kastner_one(Y, CFG, s), CFG$SEEDS, CFG$MAX_REPL, "Kastner")

# ============================================================================
#  NON-STOCHASTIC AND FIRST-STAGE PIECES, unchanged from the production scripts
# ============================================================================
fit_var_mean <- function(Y, seed, draws, burnin, jit){
   set.seed(seed); Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, jit), nrow(Y), ncol(Y))
   m <- bvar(Yj, lags=1L, draws=draws, burnin=burnin, quiet=TRUE)
   B_hat <- apply(m$PHI, c(1,2), mean)
   list(B_hat=B_hat, X=m$X, Y=m$Y, resid=m$Y - m$X %*% B_hat) }
#  rank-r truncation of the lag block; the intercept row is left alone
rank_truncate_B <- function(B_hat, p, r){
   lag_rows <- seq_len(p); sv <- svd(B_hat[lag_rows,,drop=FALSE]); B <- B_hat
   B[lag_rows,] <- sv$u[,seq_len(r),drop=FALSE] %*% diag(sv$d[seq_len(r)], r) %*%
      t(sv$v[,seq_len(r),drop=FALSE])
   list(B=B, sv=sv$d, frob_kept=sum(sv$d[seq_len(r)]^2)/sum(sv$d^2)) }
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
   list(Omega=gl$wi, rho=rho, nonzero=E, n_nonzero=as.integer(sum(E)/2)) }

# ============================================================================
#  TABLES
# ============================================================================
#  per-chain and mean (SD) MSPE, plus DM against VAPOR
mspe_and_dm <- function(losses, order, vapor_key = "VAPOR (rank 3)"){
   #  a method whose every chain errored has no loss vector; it stays in the
   #  table with an explicit "no chain completed" rather than silently
   #  becoming NaN, and it is left out of the Diebold-Mariano comparison
   have <- vapply(order, function(nm) length(losses[[nm]]) > 0, logical(1))
   lv <- Reduce(`+`, losses[[vapor_key]])/length(losses[[vapor_key]])
   per <- do.call(rbind, lapply(order, function(nm){
      if (!have[[nm]])
         return(data.frame(Method = nm, n_chains = 0L,
                           per_chain = "no chain completed", MSPE = "-"))
      m <- vapply(losses[[nm]], mean, numeric(1))
      data.frame(Method = nm, n_chains = length(m),
                 per_chain = paste(sprintf("%.4f", m), collapse=", "),
                 MSPE = if (length(m) > 1) sprintf("%.4f (%.4f)", mean(m), sd(m))
                        else sprintf("%.4f", mean(m))) }))
   dm <- do.call(rbind, lapply(order, function(nm){
      if (identical(nm, vapor_key) || !have[[nm]])
         cbind(data.frame(Method=nm),
               data.frame(Difference=NA, SE=NA, DM=NA, p=NA,
                          Won = if (have[[nm]]) "-" else "no chain completed"))
      else cbind(data.frame(Method=nm),
                 dm_test(Reduce(`+`, losses[[nm]])/length(losses[[nm]]), lv)) }))
   list(per_chain = per, dm = dm) }

#  VAPOR edge list from the pooled draws, with between-chain R-hat, bulk ESS
#  and the number of chains that select each edge on its own
edge_table <- function(chains, vars){
   p <- length(vars); lab <- pair_labels(vars)
   DR <- lapply(chains, pair_draws, p = p)
   sets <- lapply(DR, ci_indices)
   pool <- do.call(rbind, DR); pooled <- ci_indices(pool)
   tab <- do.call(rbind, lapply(pooled, function(j){
      Mx <- do.call(cbind, lapply(DR, function(M) M[,j]))
      data.frame(edge = lab[j],
                 pcor = round(mean(pool[,j]), 3),
                 lo = round(quantile(pool[,j], .025), 3),
                 hi = round(quantile(pool[,j], .975), 3),
                 Rhat = round(posterior::rhat(Mx), 3),
                 bulk_ESS = round(posterior::ess_bulk(Mx), 0),
                 tail_ESS = round(posterior::ess_tail(Mx), 0),
                 chains_selecting = sum(vapply(sets, function(S) j %in% S, logical(1)))) }))
   tab <- tab[order(-abs(tab$pcor)), ]
   ut <- attr(DR[[1]], "ut")
   E <- matrix(FALSE, p, p)
   E[cbind(ut[pooled,1], ut[pooled,2])] <- TRUE
   E[cbind(ut[pooled,2], ut[pooled,1])] <- TRUE
   list(table = tab, E = E, pooled_set = pooled, per_chain_sets = sets,
        pc_mean = { acc <- matrix(0,p,p)
           for (f in chains) acc <- acc + pcor_from_omega(omega_hat(f))
           acc/length(chains) },
        draws = DR) }

#  ============================================================================
#  ONE partial-correlation convention for every reported summary (item 4):
#  for a Bayesian method it is the MEAN OF THE PER-DRAW PARTIAL CORRELATIONS,
#  pooled over all chains.  Not the partial correlations of the mean precision,
#  and never a single chain.  VAR + Glasso and the sample covariance have no
#  draws and keep their point estimates.  The MSPE path is untouched: it still
#  uses the per-chain precision point estimate in OM.
#  ============================================================================
#  chains that carry Lambda + residual_prec (VAPOR, PrecFact, VAR + PrecFact)
pcor_mean_draws <- function(chains, p){
   if (!length(chains)) return(NULL)
   ut <- which(upper.tri(matrix(0,p,p)), arr.ind = TRUE)
   M <- do.call(rbind, lapply(chains, pair_draws, p = p))
   v <- colMeans(M)
   PC <- matrix(NA_real_, p, p)
   PC[cbind(ut[,1], ut[,2])] <- v; PC[cbind(ut[,2], ut[,1])] <- v
   diag(PC) <- 1
   attr(PC, "n_chains") <- length(chains); attr(PC, "n_draws") <- nrow(M); PC }
#  Kastner, whose draws are stored as per-draw precisions
#  Kastner chains carry their per-draw precisions in one of three shapes:
#  the raw array (prec_draws), a precomputed within-chain mean of the per-draw
#  partial correlations (pcor_mean_draws), or a file holding the array.  All
#  three give the same pooled mean when the chains have equal draw counts.
pcor_mean_prec <- function(chains, p){
   if (!length(chains)) return(NULL)
   acc <- matrix(0, p, p); n <- 0L; nch <- 0L
   for (z in chains){
      if (!is.null(z$pcor_mean_draws)){
         nd <- if (!is.null(z$n_draws)) z$n_draws else 1L
         acc <- acc + z$pcor_mean_draws * nd; n <- n + nd; nch <- nch + 1L; next }
      #  the p = 29 chains carry their precision draws inline.  A third shape,
      #  a `prec_draws_file` basename read from results/fits/kastner_p100_prec_draws/,
      #  was removed: that directory has never existed in this repository, and
      #  every chain that names one also carries pcor_mean_draws inline, so the
      #  branch above always won.  Verified: the Kastner p = 100 edge counts,
      #  all 26 per-edge ranks and the recovery summary in results_summary.md
      #  are reproduced exactly without it.
      A <- z$prec_draws
      if (is.null(A)) next
      nd <- dim(A)[3]
      for (s in seq_len(nd)) acc <- acc + pcor_from_omega(A[,,s])
      n <- n + nd; nch <- nch + 1L }
   if (!n) return(NULL)
   PC <- acc/n
   attr(PC, "n_chains") <- nch; attr(PC, "n_draws") <- n; PC }
#  the full PCs list under the one convention
pcs_one_convention <- function(CH, gl_Omega, Om_sample, p, gl_rr3_Omega = NULL){
   L <- list(
      VAPOR                     = pcor_mean_draws(CH$vapor,  p),
      Kastner                   = pcor_mean_prec (CH$kastner, p),
      `VAR + Glasso`            = pcor_from_omega(gl_Omega),
      `Sample covariance`       = pcor_from_omega(Om_sample),
      `PrecFact (no VAR)`       = pcor_mean_draws(CH$pf_raw, p),
      `VAR + PrecFact`          = pcor_mean_draws(CH$pf_rr,  p),
      `VAR + Glasso (rank 3)`   = if (!is.null(gl_rr3_Omega)) pcor_from_omega(gl_rr3_Omega),
      `VAR + PrecFact (rank 3)` = pcor_mean_draws(CH$pf_rr3, p))
   L[!vapply(L, is.null, logical(1))] }

#  edge counts at thresholds on |mean per-draw partial correlation|
threshold_counts <- function(PCs, thresholds){
   do.call(cbind, c(list(data.frame(threshold = thresholds)),
      lapply(names(PCs), function(nm){
         v <- vapply(thresholds, function(t)
            sum(abs(PCs[[nm]][upper.tri(PCs[[nm]])]) > t), numeric(1))
         setNames(data.frame(v), nm) }))) }

#  rank of each reported edge within each method's own ordering of |pcor|
edge_recovery_ranks <- function(edges, PCs, vars){
   idx <- cbind(match(edges$v1, vars), match(edges$v2, vars))
   do.call(cbind, lapply(PCs, function(PC){
      v <- abs(PC[upper.tri(PC)])
      ut <- which(upper.tri(PC), arr.ind=TRUE)
      key <- paste(pmin(ut[,1],ut[,2]), pmax(ut[,1],ut[,2]), sep="-")
      want <- paste(pmin(idx[,1],idx[,2]), pmax(idx[,1],idx[,2]), sep="-")
      r <- rank(-v, ties.method="min")
      r[match(want, key)] })) }

#  share of the top-N pairs that join two series from the same official
#  FRED-QD group
group_purity <- function(PCs, vars, groups_csv, top_n){
   g <- read.csv(groups_csv, stringsAsFactors = FALSE)
   #  fredqd_groups.csv carries series, group_code, group_name, source
   stopifnot(all(c("series","group_code") %in% names(g)))
   grp <- setNames(g$group_code, g$series)[vars]
   if (anyNA(grp)) warning(sprintf("%d series absent from %s", sum(is.na(grp)),
                                   basename(groups_csv)))
   do.call(rbind, lapply(names(PCs), function(nm){
      PC <- PCs[[nm]]; ut <- which(upper.tri(PC), arr.ind=TRUE)
      v <- abs(PC[upper.tri(PC)])
      o <- order(-v)[seq_len(min(top_n, length(v)))]
      same <- grp[ut[o,1]] == grp[ut[o,2]]
      data.frame(method = nm, top_n = length(o),
                 same_group = sum(same, na.rm=TRUE),
                 purity = round(mean(same, na.rm=TRUE), 3)) })) }

#  5 x 5 pooled partial-correlation block with 95% intervals
rate_block <- function(DR, vars, series){
   p <- length(vars); ut <- attr(DR[[1]], "ut"); pool <- do.call(rbind, DR)
   M <- matrix("", length(series), length(series), dimnames=list(series, series))
   for (a in seq_along(series)) for (b in seq_along(series)){
      if (a == b){ M[a,b] <- "1"; next }
      ia <- match(series[a], vars); ib <- match(series[b], vars)
      j <- which(ut[,1] == min(ia,ib) & ut[,2] == max(ia,ib))
      M[a,b] <- sprintf("%.3f [%.3f,%.3f]", mean(pool[,j]),
                        quantile(pool[,j],.025), quantile(pool[,j],.975)) }
   M }

#  convergence table: log-likelihood, spectral radius of B, and every reported
#  edge, all on one common set of chains
convergence_table <- function(chains, edges, vars, p){
   q <- list(loglik = lapply(chains, function(f) colSums(f$loglik)),
             rho_B  = lapply(chains, function(f)
                vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s], p)),
                       numeric(1))))
   for (e in seq_len(nrow(edges))){
      i <- match(edges$v1[e], vars); j <- match(edges$v2[e], vars)
      q[[paste(edges$v1[e], edges$v2[e], sep="_")]] <-
         lapply(chains, function(f) vapply(seq_len(dim(f$Lambda)[3]),
            function(s) pcor_from_omega(omega_one(f,s))[i,j], numeric(1))) }
   do.call(rbind, lapply(names(q), function(nm){
      M <- do.call(cbind, q[[nm]])
      data.frame(quantity = nm, Rhat = round(posterior::rhat(M), 3),
                 bulk_ESS = round(posterior::ess_bulk(M), 0),
                 tail_ESS = round(posterior::ess_tail(M), 0)) })) }

# ============================================================================
#  POSTERIOR PREDICTIVE CHECK.
#
#  NOTE ON PROVENANCE.  The nine summaries below are the nine rows of the
#  previously saved ppc table (max_eig, eig2, mean_abs_cor, q90_abs_cor,
#  max_abs_pcor, prop_pcor_02, mean_sd, kurtosis, resid_ac1), but the code that
#  produced that table is not in the repository -- it was run interactively and
#  only its output was saved.  The definitions here are therefore a restatement,
#  written out explicitly so they are reproducible from now on.  They are not
#  guaranteed to match the old numbers to the digit.
# ============================================================================
#  The nine summaries, computed on a panel Y.  `resid` supplies the residual
#  matrix y - B y_{t-1} that resid_ac1 is defined on; when it is absent that
#  row falls back to the autocorrelation of Y itself and is flagged by the
#  caller.  Every other summary is a function of the panel alone.
ppc_stats <- function(Y, resid = NULL){
   C  <- cor(Y); oc <- abs(C[upper.tri(C)])
   ev <- eigen(C, symmetric=TRUE, only.values=TRUE)$values
   S  <- cov(Y); Om <- tryCatch(solve(S), error=function(e) MASS::ginv(S))
   PC <- pcor_from_omega(Om); op <- abs(PC[upper.tri(PC)])
   kurt <- mean(apply(Y, 2, function(x){ z <- (x-mean(x))/sd(x); mean(z^4) }))
   R <- if (is.null(resid)) Y else resid
   ac1 <- mean(apply(R, 2, function(x) cor(x[-1], x[-length(x)])))
   c(max_eig      = ev[1],
     eig2         = ev[2],
     mean_abs_cor = mean(oc),
     q90_abs_cor  = unname(quantile(oc, 0.90)),
     max_abs_pcor = max(op),
     prop_pcor_02 = mean(op < 0.2),
     mean_sd      = mean(apply(Y, 2, sd)),
     kurtosis     = kurt,
     resid_ac1    = ac1) }

#  ONE-STEP-AHEAD replication, as the paper describes it: for each draw s,
#  y^rep_t ~ N(B_s y_{t-1}, Omega_s^{-1}) with y_{t-1} the OBSERVED lag, not the
#  previously simulated value.  The replicated panel is rows 2..T.  Because
#  y^rep - B_s y_{t-1} is exactly the simulated innovation, resid_ac1 on the
#  replicate is the autocorrelation of that innovation.
#
#  This replaces a recursive simulation (y^rep_t = B_s y^rep_{t-1} + e_t) and a
#  resid_ac1 taken on the panel rather than on the residual.  Both differed
#  from the paper; the recursive version compounds B_s over 99 periods, which
#  for a draw near the unit circle produces panels orders of magnitude larger
#  than the data.
ppc_one <- function(fit, s, Y, p){
   B  <- fixB(fit$B[,,s], p); Tn <- nrow(Y)
   M  <- Y[-Tn, , drop=FALSE] %*% t(B)            # one-step-ahead means, t = 2..T
   Om <- omega_one(fit, s)
   Sig <- tryCatch(chol2inv(chol(Om)), error=function(e) MASS::ginv(Om))
   L   <- tryCatch(chol(Sig), error=function(e) chol(Sig + diag(1e-8, p)))
   E   <- matrix(rnorm((Tn-1)*p), Tn-1, p) %*% L   # innovations, rows 2..T
   Yrep <- M + E
   colnames(Yrep) <- colnames(E) <- colnames(Y)
   Yobs <- Y[-1, , drop=FALSE]
   list(Yrep = Yrep, resid_rep = E, Yobs = Yobs, resid_obs = Yobs - M,
        Omega = Om, residual_prec = fit$residual_prec[s, ]) }

#  The check.  For every sampled draw both the observed and the replicated
#  statistics are computed with that draw's B_s, since resid_ac1 depends on it;
#  the other eight observed values are constant across draws.  A draw whose
#  replicated panel contains any |value| > BIG is reported separately and can
#  be excluded, which is what `exclude_degenerate` does.
ppc_table <- function(chains, Y, p, n_rep, seed, BIG = 1e3,
                      exclude_degenerate = FALSE, return_detail = FALSE){
   set.seed(seed)
   picks <- do.call(rbind, lapply(seq_along(chains), function(ci)
      cbind(ci, sample(dim(chains[[ci]]$Lambda)[3], ceiling(n_rep/length(chains)),
                       replace = TRUE))))
   picks <- picks[seq_len(min(n_rep, nrow(picks))), , drop = FALSE]
   D <- lapply(seq_len(nrow(picks)), function(i){
      r <- picks[i,]; z <- ppc_one(chains[[r[1]]], r[2], Y, p)
      list(rep = ppc_stats(z$Yrep, z$resid_rep),
           obs = ppc_stats(z$Yobs, z$resid_obs),
           max_abs = max(abs(z$Yrep)),
           at_floor = sum(z$residual_prec <= 1e-8),
           kappa = tryCatch(kappa(z$Omega, exact = TRUE), error = function(e) NA_real_),
           chain = r[1], draw = r[2]) })
   bad <- vapply(D, function(z) z$max_abs > BIG, logical(1))
   use <- if (exclude_degenerate) !bad else rep(TRUE, length(D))
   REP <- do.call(rbind, lapply(D[use], function(z) z$rep))
   OBS <- do.call(rbind, lapply(D[use], function(z) z$obs))
   tab <- data.frame(statistic = colnames(REP),
      observed = round(colMeans(OBS), 3),
      rep_mean = round(colMeans(REP), 3),
      lo = round(apply(REP, 2, quantile, 0.025), 3),
      hi = round(apply(REP, 2, quantile, 0.975), 3),
      ppp = round(colMeans(REP >= OBS), 3), row.names = NULL)
   attr(tab, "n_draws") <- sum(use); attr(tab, "n_degenerate") <- sum(bad)
   if (return_detail) list(table = tab, detail = D, degenerate = bad) else tab }

# ============================================================================
#  FIGURES.  Trace and posterior-predictive panels; the four application
#  figures are produced by the existing heatmap / network / term-structure
#  scripts, which the driver calls after the pipeline has written its .rds.
# ============================================================================
#  Trace panels with one shared legend BELOW the panels, and each panel titled
#  with its own R-hat and bulk ESS so the diagnostic sits next to the picture.
fig_trace <- function(chains, edges, vars, p, out, n_edge = 6){
   top <- head(edges, n_edge)
   series <- c(list(`Log-likelihood` = lapply(chains, function(f) colSums(f$loglik)),
                    `Spectral radius of B` = lapply(chains, function(f)
                       vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s],p)),
                              numeric(1)))),
               setNames(lapply(seq_len(nrow(top)), function(e){
                  i <- match(top$v1[e], vars); j <- match(top$v2[e], vars)
                  lapply(chains, function(f) vapply(seq_len(dim(f$Lambda)[3]),
                     function(s) pcor_from_omega(omega_one(f,s))[i,j], numeric(1))) }),
                  paste(top$v1, top$v2, sep = " - ")))
   np <- length(series)
   png(out, width = 12, height = 2.05*np + 1.1, units = "in", res = 200, bg = "white")
   on.exit(dev.off())
   #  one row per quantity, plus a short final row that holds the legend only
   layout(matrix(seq_len(np + 1L), ncol = 1), heights = c(rep(1, np), 0.42))
   cols <- c("#1b6ca8", "#c0392b", "#27773c", "#8e44ad")
   for (k in seq_along(series)){
      nm <- names(series)[k]; last <- (k == length(series))
      X <- series[[nm]]; M <- do.call(cbind, X)
      rh <- round(posterior::rhat(M), 3); es <- round(posterior::ess_bulk(M), 0)
      #  only the bottom panel carries the x-axis title, so the label has room
      #  instead of being squeezed between panels
      par(mar = c(if (last) 3.8 else 2.2, 4.8, 2.6, 1.0), mgp = c(2.4, 0.7, 0))
      yl <- range(unlist(X), finite = TRUE)
      plot(NA, xlim = c(1, length(X[[1]])), ylim = yl,
           xlab = if (last) "retained draw" else "",
           ylab = "", main = sprintf("%s   (R-hat %.3f, ESS %d)", nm, rh, es),
           cex.main = 1.2, font.main = 1)
      for (c in seq_along(X)) lines(X[[c]], col = cols[(c-1) %% length(cols) + 1], lwd = 0.8) }
   par(mar = c(0, 0, 0, 0))
   plot(NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
   legend("center", horiz = TRUE, bty = "n", cex = 1.15, seg.len = 2.4,
          legend = sprintf("chain %d (seed %s)", seq_along(chains),
                           vapply(chains, function(f){
                              sd <- attr(f, "seed")
                              if (is.null(sd)) "?" else as.character(sd) }, character(1))),
          col = cols[seq_along(chains)], lwd = 2.2)
   invisible(out) }

#  Readable names for the nine summaries, used in the figure titles.
PPC_NAMES <- c(
   max_eig      = "Largest eigenvalue",
   eig2         = "Second eigenvalue",
   mean_abs_cor = "Mean |correlation|",
   q90_abs_cor  = "90th pct |correlation|",
   max_abs_pcor = "Max |partial correlation|",
   prop_pcor_02 = "Share |partial cor| < 0.2",
   mean_sd      = "Mean standard deviation",
   kurtosis     = "Kurtosis",
   resid_ac1    = "Residual autocorrelation")

fig_ppc <- function(ppc, chains, Y, p, out, n_rep, seed, BIG = 1e3){
   z <- ppc_table(chains, Y, p, n_rep, seed, BIG = BIG, return_detail = TRUE)
   REP <- do.call(rbind, lapply(z$detail, function(d) d$rep))
   OBS <- colMeans(do.call(rbind, lapply(z$detail, function(d) d$obs)))
   png(out, width = 12, height = 9.6, units = "in", res = 200, bg = "white")
   on.exit(dev.off())
   par(mfrow = c(3,3), mar = c(3.6, 3.8, 3.0, 1.0), mgp = c(2.3, 0.7, 0),
       oma = c(0, 0, 2.2, 0))
   for (k in seq_along(OBS)){
      st <- names(OBS)[k]
      nm <- if (st %in% names(PPC_NAMES)) PPC_NAMES[[st]] else st
      h <- hist(REP[,k], breaks = 24, plot = FALSE)
      plot(h, col = "grey88", border = "grey65",
           main = sprintf("%s (p = %.2f)", nm, z$table$ppp[k]),
           cex.main = 1.25, font.main = 1,
           xlab = "replicated value", ylab = "frequency",
           xlim = range(c(h$breaks, OBS[k])))
      abline(v = OBS[k], col = "#c0392b", lwd = 2.6) }
   mtext(sprintf("Posterior predictive check, one step ahead, p = %d (%d draws; red rule = observed)",
                 p, attr(z$table, "n_draws")),
         outer = TRUE, line = 0.4, cex = 1.15)
   invisible(out) }

# ============================================================================
#  MARKDOWN
# ============================================================================
md_section <- function(file, title, x, note = NULL){
   cat(sprintf("\n## %s\n\n", title), file = file, append = TRUE)
   if (!is.null(note)) cat(sprintf("%s\n\n", note), file = file, append = TRUE)
   cat("```\n", file = file, append = TRUE)
   cat(paste(utils::capture.output(
      if (is.data.frame(x)) print(x, row.names = FALSE) else print(x)),
      collapse = "\n"), file = file, append = TRUE)
   cat("\n```\n", file = file, append = TRUE) }

#  ---------------------------------------------------------------------------
#  Provenance banner, printed once by every script that sources this file.
#  ---------------------------------------------------------------------------
.vapor_provenance <- function(){
   gh <- tryCatch(system2("git", c("-C", "/Users/peterdunson/Desktop/vapor",
                                   "rev-parse", "--short", "HEAD"),
                          stdout = TRUE, stderr = FALSE)[1],
                  error = function(e) NA_character_)
   dirty <- tryCatch(length(system2("git", c("-C", "/Users/peterdunson/Desktop/vapor",
                                             "status", "--porcelain"),
                                    stdout = TRUE, stderr = FALSE)) > 0,
                     error = function(e) NA)
   fo <- formals(VAPOR::vapor)
   cat(sprintf("VAPOR %s (%s) | git %s%s | DP base measure as = %s, bs = %s, a_dir = %s, b_dir = %s | %s\n",
               as.character(utils::packageVersion("VAPOR")),
               find.package("VAPOR"), gh, if (isTRUE(dirty)) " +uncommitted" else "",
               fo$as, fo$bs, fo$a_dir, fo$b_dir, format(Sys.time(), "%Y-%m-%d %H:%M"))) }

suppressMessages(.vapor_provenance())
