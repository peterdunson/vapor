# ============================================================================
#  Matched-effort dependence comparison: 4 long chains per stochastic method
#
#  Every stochastic competitor gets the same MCMC effort as VAPOR's long
#  chains: 20,000 iterations, 10,000 burn-in, thin 10, so 1,000 retained
#  draws per chain and 4,000 pooled. Competitor chain seeds are
#  SEED+101 .. SEED+104, distinct from the production fits.
#
#  MSPE and Diebold-Mariano are NOT touched here: those stay with the three
#  production fits per method, in p29_results.rds / p100_results.rds.
#
#  WHAT THIS PRODUCES, per dimension, in p<tag>_long_competitors.rds:
#    * pooled posterior-mean partial correlations per method, formed exactly
#      as pc_vapor is: the mean of pcor(Omega_s) over all draws of all chains
#    * the monitored log-likelihood series per chain
#    * R-hat and ESS per method
#
#  THE MONITORED QUANTITY (identical function for all four methods):
#     loglik(Omega) = (n/2) * [ log|Omega| - tr(S Omega) ] - (n p / 2) log(2 pi)
#  with S = E'E / n and E the residuals of the SHARED rank-3 first-stage
#  VAR(1) (first_stage$resid_rr). The same E, the same S and the same function
#  are used for Kastner, both PrecFact rows and VAPOR. R-hat compares chains
#  WITHIN a method, so holding S fixed across methods makes the diagnostic
#  comparable; it is a convergence monitor, not a cross-method likelihood
#  comparison (PrecFact (no VAR) models the raw series, not innovations).
#
#  CACHING: each chain is summarised to a ~100 KB object and cached under
#  CACHE_DIR, keyed by dimension/method/seed. Re-running skips finished chains,
#  so an interrupted run resumes instead of starting over.
#
#  COST: one Kastner long chain at p = 100 was measured at 9.1 minutes, so
#  the Kastner block at p = 100 is ~36 minutes on its own. Chains run
#  sequentially; the cache makes a resumed run skip whatever finished.
# ============================================================================
suppressMessages({library(VAPOR); library(PrecFact); library(bayesianVARs)
                  library(coda)})

APP       <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
CACHE_DIR <- file.path(APP, "long_chain_cache")
dir.create(CACHE_DIR, showWarnings = FALSE)

#  Chains run SEQUENTIALLY. Each chain is reduced to summaries and its model
#  object dropped before the next starts, so peak memory stays at one chain
#  (~0.5 GB for p=100 Kastner) rather than four at once on a 8 GB machine.

#  SMOKE = TRUE runs the whole pipeline at toy length to validate the code
#  paths in seconds. It writes to a separate cache so it cannot contaminate
#  real results, and its numbers are meaningless.
SMOKE <- as.logical(Sys.getenv("LC_SMOKE", unset = "FALSE"))
if (SMOKE){ CACHE_DIR <- file.path(CACHE_DIR, "smoke"); dir.create(CACHE_DIR, showWarnings = FALSE) }

#  nrun = total iterations, burn discarded, thin applied -> retained draws.
#  `draws` is bvar's post-burn-in iteration count (bvar takes burnin and draws
#  separately, unlike vapor's single nrun), so draws = nrun - burn.
LONG <- if (SMOKE) list(nrun = 400L,   burn = 200L,   thin = 2L,  draws = 200L) else
                   list(nrun = 20000L, burn = 10000L, thin = 10L, draws = 10000L)

SPEC <- list(
   p29  = list(res = "p29_results.rds",  out = "p29_long_competitors.rds",
               SEED = 32,   K = 10, RANK = 3, VAR_JIT = 1e-3, JIT = 1e-3,
               vapor_chains_file = NULL,       # never saved; refit (cheap at p=29)
               vapor_save_file = "p29_chains_final.rds",
               vapor_chain_seeds = 100 + 1:4),
   p100 = list(res = "p100_results.rds", out = "p100_long_competitors.rds",
               SEED = 2024, K = 25, RANK = 3, VAR_JIT = 1e-3, JIT = 5e-3,
               vapor_chains_file = "p100_chains_final.rds",
               vapor_save_file = NULL,         # already saved
               vapor_chain_seeds = 100 + 1:4))

# ---- the monitored quantity -------------------------------------------------
LOGLIK_FORMULA <-
   "loglik(Omega) = (n/2)[log|Omega| - tr(S Omega)] - (n p/2) log(2 pi), S = E'E/n, E = rank-3 first-stage VAR(1) residuals"

make_loglik_fn <- function(E){
   n <- nrow(E); p <- ncol(E); S <- crossprod(E)/n
   function(Om){
      ch <- tryCatch(chol(Om), error = function(e) NULL)
      ld <- if (is.null(ch)) as.numeric(determinant(Om, logarithm = TRUE)$modulus)
            else 2*sum(log(diag(ch)))
      (n/2)*(ld - sum(S*Om)) - (n*p/2)*log(2*pi) } }

pcor_of <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }

#  accumulate pcor and the monitored series over a chain's draws, so the
#  caller never has to hold the draws themselves
summarise_draws <- function(omega_at, n_draws, p, loglik_fn){
   acc <- matrix(0, p, p); ll <- numeric(n_draws)
   for (s in seq_len(n_draws)){
      Om <- omega_at(s)
      acc <- acc + pcor_of(Om)
      ll[s] <- loglik_fn(Om) }
   list(pcor_sum = acc, n_draws = n_draws, loglik = ll) }

# ============================================================================
#  one chain per method
# ============================================================================
#  RETRY POLICY, shared by both samplers: attempt 1 is exactly the plain
#  set.seed(seed) fit, so chains cached before this policy existed stay valid.
#  Only a failure triggers a reseed-and-jitter retry, mirroring run_pf. The
#  attempt count is recorded and warned about, because a retry changes that
#  chain's draws.
chain_kastner <- function(seed, Y, k, p, loglik_fn, long, max_try = 6){
   sp <- specify_prior_sigma(M = p, type = "factor", factor_factors = as.integer(k),
                             factor_priorfacloadtype = "colwiseng",
                             factor_heteroskedastic = TRUE, quiet = TRUE)
   m <- NULL; a <- 0
   while (is.null(m) && a < max_try){
      a <- a + 1
      if (a == 1){ set.seed(seed)
         Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-3), nrow(Y), ncol(Y))
      } else { set.seed(seed + a*7919)
         Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-3 * 2^(a-1)), nrow(Y), ncol(Y)) }
      m <- tryCatch(
   #  bvar counts `draws` post-burn-in and thins them, so burn + draws
   #  iterations total with draws/thin retained: matched to VAPOR's long chain
         bvar(Yj, lags = 1L, draws = long$draws, burnin = long$burn,
              thin = long$thin, prior_sigma = sp, quiet = TRUE),
         error = function(e){ message(sprintf("      bvar attempt %d failed: %s",
                                              a, conditionMessage(e))); NULL }) }
   if (is.null(m)) stop(sprintf("bvar failed after %d attempts at seed %d", max_try, seed))
   if (a > 1) warning(sprintf("Kastner seed %d needed %d attempts", seed, a))
   V <- vcov(m); nd <- dim(V)[4]
   out <- summarise_draws(function(s) solve(V[1,,,s]), nd, p, loglik_fn)
   out$attempts <- a
   rm(m, V); gc(verbose = FALSE)
   out }

chain_precfact <- function(seed, Ymat, k, p, loglik_fn, long, max_try = 6){
   res <- NULL; a <- 0
   while (is.null(res) && a < max_try){
      a <- a + 1
      if (a == 1){ set.seed(seed); Yj <- Ymat }
      else { set.seed(seed + a*7919)
         Yj <- Ymat + matrix(rnorm(length(Ymat), 0, 1e-6), nrow(Ymat)) }
      res <- tryCatch(
         prec_est(as = 0.1, bs = 0.1, a = 0.5, b = 0.5, k = k,
                  nrun = long$nrun, thin = long$thin, Y = Yj,
                  alldiff = 0, a_dir = 0.1, b_dir = 0.1),
         error = function(e){ message(sprintf("      prec_est attempt %d failed: %s",
                                              a, conditionMessage(e))); NULL }) }
   if (is.null(res)) stop(sprintf("prec_est failed after %d attempts at seed %d", max_try, seed))
   if (a > 1) warning(sprintf("PrecFact seed %d needed %d attempts", seed, a))
   nd_all <- dim(res$Lambda)[3]
   keep <- (floor(long$burn/long$thin) + 1L):nd_all      # drop the burn-in
   out <- summarise_draws(function(i){
      s <- keep[i]
      Om <- tcrossprod(res$Lambda[,,s]); diag(Om) <- diag(Om) + res$residuals[s,]; Om },
      length(keep), p, loglik_fn)
   out$attempts <- a
   rm(res); gc(verbose = FALSE)
   out }

chain_vapor_fit <- function(seed, Y, k, rank, p, loglik_fn, long, jitter,
                            keep_fit = FALSE){
   f <- vapor(Y, lags = 1, k = k, rank = rank, prior_B = "dl",
              nrun = long$nrun, burn = long$burn, thin = long$thin,
              seed = seed, jitter = jitter, init = "random")
   nd <- dim(f$Lambda)[3]
   out <- summarise_draws(function(s){
      Om <- tcrossprod(f$Lambda[,,s]); diag(Om) <- diag(Om) + f$residual_prec[s,]; Om },
      nd, p, loglik_fn)
   if (keep_fit) out$fit <- f else { rm(f); gc(verbose = FALSE) }
   out }

#  the 95% credible-interval edge rule, over a list of chains
ci_set_chains <- function(fits){
   stk <- unlist(lapply(fits, function(f)
      lapply(seq_len(dim(f$Lambda)[3]), function(s){
         Om <- tcrossprod(f$Lambda[,,s]); diag(Om) <- diag(Om) + f$residual_prec[s,]
         pcor_of(Om) })), recursive = FALSE)
   A <- simplify2array(stk)
   lo <- apply(A, c(1,2), quantile, 0.025); hi <- apply(A, c(1,2), quantile, 0.975)
   E <- (lo > 0 | hi < 0); diag(E) <- FALSE
   list(E = E, pm = apply(A, c(1,2), mean)) }

summarise_existing_vapor <- function(fit, p, loglik_fn){
   nd <- dim(fit$Lambda)[3]
   summarise_draws(function(s){
      Om <- tcrossprod(fit$Lambda[,,s]); diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om },
      nd, p, loglik_fn) }

# ============================================================================
#  driver
# ============================================================================
cached <- function(tag, method, seed, fn){
   f <- file.path(CACHE_DIR, sprintf("%s_%s_seed%d.rds", tag, method, seed))
   if (file.exists(f)){ cat(sprintf("    seed %d: cached\n", seed)); return(readRDS(f)) }
   t0 <- Sys.time()
   z <- fn()
   z$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
   z$seed <- seed
   saveRDS(z, f)
   cat(sprintf("    seed %d: %.1f s (%d draws%s)\n", seed, z$seconds, z$n_draws,
               if (!is.null(z$attempts) && z$attempts > 1)
                  sprintf(", %d attempts", z$attempts) else ""))
   z }

run_dimension <- function(tag, sp){
   R <- readRDS(file.path(APP, sp$res))
   Y_train <- R$Y_train; p <- R$p; vars <- R$vars
   E <- R$first_stage$resid_rr
   loglik_fn <- make_loglik_fn(E)
   resid_rr_scaled <- scale(E)
   cs <- sp$SEED + 101:104

   cat(sprintf("\n################  %s  (p = %d)  ################\n", tag, p))
   cat(sprintf("long chains: %d iterations, %d burn-in, thin %d -> %d draws/chain\n",
               LONG$nrun, LONG$burn, LONG$thin, LONG$burn/LONG$thin))
   cat(sprintf("competitor seeds: %s | sequential%s\n",
               paste(cs, collapse = ", "), if (SMOKE) "  [SMOKE]" else ""))

   jobs <- list(
      "Kastner (VAR+FSV)" = lapply(cs, function(s) function()
         chain_kastner(s, Y_train, sp$K, p, loglik_fn, LONG)),
      "PrecFact (no VAR)" = lapply(cs, function(s) function()
         chain_precfact(s, Y_train, sp$K, p, loglik_fn, LONG)),
      "VAR + PrecFact"    = lapply(cs, function(s) function()
         chain_precfact(s, resid_rr_scaled, sp$K, p, loglik_fn, LONG)))

   res <- list()
   for (meth in names(jobs)){
      cat(sprintf("\n  -- %s --\n", meth))
      fns <- jobs[[meth]]
      res[[meth]] <- lapply(seq_along(cs), function(i)
         cached(tag, gsub("[^A-Za-z0-9]", "", meth), cs[i], fns[[i]])) }

   # ---- VAPOR: reuse its saved long chains where they exist ----------------
   cat("\n  -- VAPOR (rank 3) --\n")
   if (!is.null(sp$vapor_chains_file) &&
       file.exists(file.path(APP, sp$vapor_chains_file)) && !SMOKE){
      cat(sprintf("    reading %s\n", sp$vapor_chains_file))
      ch <- readRDS(file.path(APP, sp$vapor_chains_file))
      res[["VAPOR (rank 3)"]] <- lapply(seq_along(ch), function(i){
         z <- summarise_existing_vapor(ch[[i]], p, loglik_fn)
         z$seed <- sp$vapor_chain_seeds[i]; z$seconds <- NA; z })
      rm(ch); gc(verbose = FALSE)
   } else {
      vs <- sp$vapor_chain_seeds
      zs <- lapply(vs, function(s){
         t0 <- Sys.time()
         z <- chain_vapor_fit(s, Y_train, sp$K, sp$RANK, p, loglik_fn, LONG, sp$JIT,
                              keep_fit = !is.null(sp$vapor_save_file))
         z$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs")); z$seed <- s
         cat(sprintf("    seed %d: %.1f s (%d draws)\n", s, z$seconds, z$n_draws))
         z })
      if (!is.null(sp$vapor_save_file) && !SMOKE){
         fits <- lapply(zs, function(z) z$fit)
         saveRDS(list(chains = fits,
                      config = list(nrun = LONG$nrun, burn = LONG$burn, thin = LONG$thin,
                                    draws_per_chain = zs[[1]]$n_draws,
                                    seeds = vs, k = sp$K, rank = sp$RANK,
                                    jitter = sp$JIT, saved_at = format(Sys.time()))),
                 file.path(APP, sp$vapor_save_file))
         cat(sprintf("    saved %s (thin %d recorded in config)\n",
                     sp$vapor_save_file, LONG$thin))

         # ---- does the refit move the reported graph? (reporting only) -----
         newcs <- ci_set_chains(fits)
         oldpm <- R$long$pc_vapor
         cat(sprintf("    max |new pooled pcor - stored pc_vapor| = %.4f\n",
                     max(abs(newcs$pm - oldpm))))
         #  BOTH sides alphabetical: an earlier version compared variable-index
         #  ordering on one side with alphabetical on the other and reported
         #  spurious enter/leave pairs that were the same edges reversed.
         canon <- function(a, b) apply(cbind(a, b), 1, function(z) paste(sort(z), collapse = "--"))
         ij <- which(newcs$E & upper.tri(newcs$E), arr.ind = TRUE)
         new_e <- sort(if (!nrow(ij)) character(0) else canon(vars[ij[,1]], vars[ij[,2]]))
         old_e <- sort(canon(R$long$ci_edges$v1, R$long$ci_edges$v2))
         cat(sprintf("    CI edges: reported %d | under the new chains %d | identical: %s\n",
                     length(old_e), length(new_e), setequal(old_e, new_e)))
         if (!setequal(old_e, new_e)){
            cat("      would ENTER:", paste(setdiff(new_e, old_e), collapse = ", "), "\n")
            cat("      would LEAVE:", paste(setdiff(old_e, new_e), collapse = ", "), "\n") }
         cat("    (the reported graph is NOT changed; this is a check only)\n")
         for (i in seq_along(zs)) zs[[i]]$fit <- NULL
         rm(fits); gc(verbose = FALSE) }
      res[["VAPOR (rank 3)"]] <- zs }

   # ---- pool, and diagnose -------------------------------------------------
   pooled <- list(); diag_tab <- NULL; series <- list()
   for (meth in names(res)){
      zs <- res[[meth]]
      bad <- vapply(zs, function(z) inherits(z, "try-error") || is.null(z$pcor_sum), logical(1))
      if (any(bad)) stop(sprintf("%s: %d chain(s) failed", meth, sum(bad)))
      tot <- sum(vapply(zs, function(z) z$n_draws, numeric(1)))
      pooled[[meth]] <- Reduce(`+`, lapply(zs, function(z) z$pcor_sum))/tot
      ll <- lapply(zs, function(z) mcmc(z$loglik))
      series[[meth]] <- lapply(zs, function(z) z$loglik)
      mcl <- mcmc.list(ll)
      rh <- tryCatch(gelman.diag(mcl, autoburnin = FALSE)$psrf[1], error = function(e) NA_real_)
      diag_tab <- rbind(diag_tab, data.frame(
         Method = meth, chains = length(zs), draws_per_chain = zs[[1]]$n_draws,
         pooled_draws = tot,
         Rhat_loglik = round(rh, 3),
         ESS_loglik = round(sum(effectiveSize(mcl)), 0),
         mean_loglik = round(mean(unlist(series[[meth]])), 1),
         minutes = round(sum(vapply(zs, function(z)
            if (is.na(z$seconds)) 0 else z$seconds, numeric(1)))/60, 1))) }

   cat("\n=== R-hat / ESS on the monitored log-likelihood ===\n")
   cat("  ", LOGLIK_FORMULA, "\n\n", sep = "")
   print(diag_tab, row.names = FALSE)

   out <- list(tag = tag, p = p, vars = vars,
               pooled_pcor = pooled, diagnostics = diag_tab, loglik_series = series,
               loglik_formula = LOGLIK_FORMULA,
               long = LONG, competitor_seeds = cs,
               vapor_chain_seeds = sp$vapor_chain_seeds,
               vapor_source = if (is.null(sp$vapor_chains_file)) "refit here" else sp$vapor_chains_file,
               smoke = SMOKE, run_at = format(Sys.time()))
   saveRDS(out, file.path(APP, if (SMOKE) paste0("SMOKE_", sp$out) else sp$out))
   cat(sprintf("\nSaved %s\n", if (SMOKE) paste0("SMOKE_", sp$out) else sp$out))
   invisible(out) }

invisible(lapply(names(SPEC), function(tag) run_dimension(tag, SPEC[[tag]])))
