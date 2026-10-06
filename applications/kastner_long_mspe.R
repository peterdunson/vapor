# ============================================================================
#  Held-out MSPE for Kastner using the precision matrix implied by its POOLED
#  four long chains, next to the 3-fit production MSPE.
#
#  WHY A REFIT. p<tag>_long_competitors.rds stores the pooled partial
#  correlations and the monitored log-likelihood, but not Omega. A partial
#  correlation matrix cannot be turned back into Omega: it drops the diagonal
#  scales, and the conditional mean -Omega_hh^{-1} Omega_ho y_o depends on
#  them. So the four Kastner long chains are refit, with the same seeds and
#  settings, and this time the pooled Omega IS cached -- subsequent runs load
#  it instead of refitting.
#
#  TWO CONVENTIONS, both reported, because they are not the same matrix:
#    mean_precision  = mean over draws of Omega_s           (how VAPOR's
#                      posterior_mean_omega is formed)
#    inv_mean_cov    = solve(mean over draws of Sigma_s)    (how the Kastner
#                      production fits form Omega)
#  Jensen's inequality makes these differ; the gap is reported.
#
#  Evaluation: the same 100 held-out sets as everywhere else, seeds
#  90001..90100, 30% of series hidden per quarter.
# ============================================================================
suppressMessages({library(bayesianVARs)})
APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
N_DRAWS <- 100L; HID_BASE <- 90000L; HID_FRAC <- 0.3
LONG <- list(burn = 10000L, draws = 10000L, thin = 10L)

SPEC <- list(
   p29  = list(res = "p29_results.rds",  cache = "p29_kastner_long_omega.rds",  SEED = 32),
   p100 = list(res = "p100_results.rds", cache = "p100_kastner_long_omega.rds", SEED = 2024))

predict_held_out <- function(Om, y_obs, obs, hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)
make_hid_sets <- function(base_seed, n_draws, p, T_test, frac){
   n_hid <- floor(frac*p)
   lapply(seq_len(n_draws), function(b){ set.seed(base_seed + b)
      lapply(seq_len(T_test), function(t) sort(sample(p, n_hid))) }) }
mspe_avg <- function(Om, hid, Y_test, p)
   mean(vapply(hid, function(hl)
      mean(vapply(seq_along(hl), function(t){
         h <- hl[[t]]; o <- setdiff(seq_len(p), h); yt <- Y_test[t,]
         mean((yt[h] - predict_held_out(Om, yt[o], o, h))^2) }, numeric(1))), numeric(1)))

#  one Kastner long chain, accumulating BOTH pooled quantities; same retry
#  ladder as the long-chain script so the draws match
kastner_chain_sums <- function(seed, Y, k, p, long, max_try = 6){
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
      m <- tryCatch(bvar(Yj, lags = 1L, draws = long$draws, burnin = long$burn,
                         thin = long$thin, prior_sigma = sp, quiet = TRUE),
                    error = function(e) NULL) }
   if (is.null(m)) stop(sprintf("bvar failed at seed %d", seed))
   V <- vcov(m); nd <- dim(V)[4]
   os <- matrix(0, p, p); cs <- matrix(0, p, p)
   for (s in seq_len(nd)){ Sg <- V[1,,,s]; cs <- cs + Sg; os <- os + solve(Sg) }
   rm(m, V); gc(verbose = FALSE)
   list(omega_sum = os, cov_sum = cs, n = nd, attempts = a, seed = seed) }

for (tg in names(SPEC)){
   sp <- SPEC[[tg]]
   R <- readRDS(file.path(APP, sp$res)); C <- R$config
   Y_train <- R$Y_train; Y_test <- R$Y_test; p <- R$p
   hid <- make_hid_sets(HID_BASE, N_DRAWS, p, nrow(Y_test), HID_FRAC)
   cpath <- file.path(APP, sp$cache)
   seeds <- sp$SEED + 101:104

   cat(sprintf("\n################  %s  (p = %d)  ################\n", tg, p))
   if (file.exists(cpath)){
      z <- readRDS(cpath); cat(sprintf("loaded pooled Omega from %s\n", sp$cache))
   } else {
      cat(sprintf("refitting %d Kastner long chains (seeds %s), k = %d\n",
                  length(seeds), paste(seeds, collapse = ", "), C$K))
      parts <- lapply(seeds, function(s){
         t0 <- Sys.time()
         q <- kastner_chain_sums(s, Y_train, C$K, p, LONG)
         cat(sprintf("   seed %d: %.1f min (%d draws%s)\n", s,
                     as.numeric(difftime(Sys.time(), t0, units = "mins")), q$n,
                     if (q$attempts > 1) sprintf(", %d attempts", q$attempts) else ""))
         flush.console(); q })
      tot <- sum(vapply(parts, function(q) q$n, numeric(1)))
      z <- list(mean_precision = Reduce(`+`, lapply(parts, function(q) q$omega_sum))/tot,
                mean_cov       = Reduce(`+`, lapply(parts, function(q) q$cov_sum))/tot,
                pooled_draws = tot, seeds = seeds,
                attempts = vapply(parts, function(q) q$attempts, numeric(1)),
                long = LONG, k = C$K, run_at = format(Sys.time()))
      saveRDS(z, cpath); cat(sprintf("saved %s\n", sp$cache)) }

   Om_mean <- z$mean_precision
   Om_inv  <- solve(z$mean_cov)
   prod_fits <- R$methods[["Kastner (VAR+FSV)"]]$omega
   prod_mspe <- vapply(prod_fits, mspe_avg, numeric(1), hid, Y_test, p)

   out <- data.frame(
      Estimate = c("pooled long chains: mean of per-draw Omega",
                   "pooled long chains: inverse of mean covariance",
                   "production: 3 fits, mean (sd)"),
      MSPE = c(sprintf("%.4f", mspe_avg(Om_mean, hid, Y_test, p)),
               sprintf("%.4f", mspe_avg(Om_inv,  hid, Y_test, p)),
               sprintf("%.4f (%.4f)", mean(prod_mspe), sd(prod_mspe))),
      draws = c(z$pooled_draws, z$pooled_draws, 3*1000),
      condition = c(sprintf("%.4g", kappa(Om_mean, exact = TRUE)),
                    sprintf("%.4g", kappa(Om_inv,  exact = TRUE)),
                    sprintf("%.4g", mean(vapply(prod_fits, kappa, numeric(1), exact = TRUE)))))
   cat(sprintf("\n=== Kastner held-out MSPE, %d held-out sets (seeds %d..%d) ===\n",
               N_DRAWS, HID_BASE+1, HID_BASE+N_DRAWS))
   print(out, row.names = FALSE)
   cat(sprintf("per production fit: %s\n", paste(sprintf("%.4f", prod_mspe), collapse = ", ")))
   cat(sprintf("max |mean_precision - inv_mean_cov| = %.4g\n", max(abs(Om_mean - Om_inv))))
}
