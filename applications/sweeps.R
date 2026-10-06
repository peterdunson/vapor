# ============================================================================
#  Sensitivity sweep, seed sweep, and (p = 29) the rank/k justification table,
#  re-evaluated with the 100 held-out sets.
#
#  SCOPE. This script refits ONLY the sweep models. It does not touch the
#  production fits, the long chains, the posterior predictive check, the
#  reported edges or the graph comparison, and it writes to its own file
#  (p<tag>_sweeps.rds). Nothing else is read from it.
#
#  WHY REFIT AT ALL. The sweep fits were never saved -- only their summary
#  numbers were -- so the MSPE column cannot be re-evaluated at 100 draws
#  without them. The fits are deterministic given their seeds and settings,
#  so refitting reproduces the original draws exactly; the script verifies
#  this by checking the draw-INDEPENDENT columns (edge counts, Pr(rho>1))
#  against the stored values. This time every Omega is saved.
#
#  THE PPC COLUMNS of the rank/k table (n_fail, ppp_eig, rep_eig, obs_eig)
#  and sparser_than_kastner do not depend on the held-out sets and are NOT
#  recomputed: they are carried over from the stored table, which the
#  determinism check above justifies.
# ============================================================================
suppressMessages({library(VAPOR); library(MASS)})

APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
N_DRAWS <- 100L; HID_BASE <- 90000L; HID_FRAC <- 0.3

SPEC <- list(
   p29  = list(res = "p29_results.rds",  prev = "p29_full_results.rds",
               out = "p29_sweeps.rds",  K_GRID = c(5, 8, 10, 15),
               justification = TRUE),
   p100 = list(res = "p100_results.rds", prev = "p100_full_results.rds",
               out = "p100_sweeps.rds", K_GRID = c(15, 25, 35),
               justification = FALSE))

RANKS  <- c(2, 3, 4)
NRUNS  <- list(c(5000, 2500), c(7500, 2500), c(15000, 5000))
A_LS   <- c(0.3, 0.5, 1.0)
SEEDS6 <- 1:6

# ---- helpers (same definitions as the production scripts) -------------------
omega_one <- function(fit, s){ Om <- tcrossprod(fit$Lambda[,,s])
diag(Om) <- diag(Om) + fit$residual_prec[s,]; Om }
omega_hat <- function(fit) posterior_mean_omega(fit$Lambda, fit$residual_prec)
pcor_from_omega <- function(Om){ d <- sqrt(diag(Om)); pc <- -Om/outer(d,d); diag(pc) <- 1; pc }
fixB <- function(B, p){ if (nrow(B)!=p) B <- t(B)
if (ncol(B)!=p) B <- B[,seq_len(p),drop=FALSE]; B }
spec_rad <- function(B) max(abs(eigen(B, only.values=TRUE)$values))

ci_n <- function(fit){
   nd <- dim(fit$Lambda)[3]
   A <- simplify2array(lapply(seq_len(nd), function(s) pcor_from_omega(omega_one(fit,s))))
   lo <- apply(A,c(1,2),quantile,0.025); hi <- apply(A,c(1,2),quantile,0.975)
   E <- (lo > 0 | hi < 0); diag(E) <- FALSE; sum(E)/2 }

pr_expl <- function(fit, p, n = 400){
   nd <- dim(fit$B)[3]; idx <- unique(round(seq(1,nd,length.out=min(n,nd))))
   rs <- vapply(idx, function(s) spec_rad(fixB(fit$B[,,s], p)), numeric(1))
   mean(rs > 1) }

make_hid_sets <- function(base_seed, n_draws, p, T_test, frac){
   n_hid <- floor(frac*p)
   lapply(seq_len(n_draws), function(b){ set.seed(base_seed + b)
      lapply(seq_len(T_test), function(t) sort(sample(p, n_hid))) }) }
predict_held_out <- function(Om, y_obs, obs, hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)
mspe_avg <- function(Om, hid_sets, Y_test, p)
   mean(vapply(hid_sets, function(hl)
      mean(vapply(seq_along(hl), function(t){
         h <- hl[[t]]; o <- setdiff(seq_len(p), h); yt <- Y_test[t,]
         mean((yt[h] - predict_held_out(Om, yt[o], o, h))^2) }, numeric(1))), numeric(1)))

#  the sweep fitter: same retry ladder as the production one, plus a_L
fit_v <- function(Y, k, rank, seed, nrun, burn, thin, jitter, aL = NULL, max_try = 4){
   f <- NULL; jj <- jitter; a <- 0
   while (is.null(f) && a < max_try){ a <- a + 1
      args <- list(Y, lags=1, k=k, rank=rank, prior_B="dl", nrun=nrun, burn=burn,
                   thin=thin, seed=seed, jitter=jj, init="random")
      if (!is.null(aL)) { args$a_L <- aL; args$b_L <- aL }
      f <- tryCatch(do.call(vapor, args), error=function(e) NULL)
      jj <- jj*2 }
   if (is.null(f)) stop(sprintf("vapor failed: k=%s rank=%s seed=%s", k, rank, seed))
   f }

for (tg in names(SPEC)){
   sp <- SPEC[[tg]]
   R  <- readRDS(file.path(APP, sp$res)); C <- R$config
   P  <- readRDS(file.path(APP, sp$prev))
   Y_train <- R$Y_train; Y_test <- R$Y_test; p <- R$p; T_test <- nrow(Y_test)
   hid <- make_hid_sets(HID_BASE, N_DRAWS, p, T_test, HID_FRAC)
   ev  <- function(Om) mspe_avg(Om, hid, Y_test, p)

   cat(sprintf("\n\n################  %s  (p = %d)  ################\n", tg, p))
   cat(sprintf("evaluation: %d held-out draws, seeds %d..%d\n",
               N_DRAWS, HID_BASE+1, HID_BASE+N_DRAWS))

   store <- list()
   row <- function(tag, f){
      Om <- omega_hat(f)
      store[[tag]] <<- Om
      data.frame(setting = tag, MSPE = round(ev(Om),4),
                 edges = ci_n(f), pr_expl = round(pr_expl(f, p),3)) }

   # ---- sensitivity --------------------------------------------------------
   cat("\n-- sensitivity --\n"); flush.console()
   sens <- NULL
   for (kk in sp$K_GRID){ cat(sprintf("\r   k=%d    ", kk)); flush.console()
      sens <- rbind(sens, row(sprintf("k=%d", kk),
         fit_v(Y_train, kk, C$RANK, C$SEED, C$NRUN, C$BURN, C$THIN, C$JIT))) }
   for (rr in RANKS){ cat(sprintf("\r   rank=%d ", rr)); flush.console()
      sens <- rbind(sens, row(sprintf("rank=%d", rr),
         fit_v(Y_train, C$K, rr, C$SEED, C$NRUN, C$BURN, C$THIN, C$JIT))) }
   for (nb in NRUNS){ cat(sprintf("\r   nrun=%d ", nb[1])); flush.console()
      sens <- rbind(sens, row(sprintf("nrun=%d", nb[1]),
         fit_v(Y_train, C$K, C$RANK, C$SEED, nb[1], nb[2], C$THIN, C$JIT))) }
   for (aL in A_LS){ cat(sprintf("\r   a_L=%.1f ", aL)); flush.console()
      sens <- rbind(sens, row(sprintf("a_L=%.1f", aL),
         fit_v(Y_train, C$K, C$RANK, C$SEED, C$NRUN, C$BURN, C$THIN, C$JIT, aL = aL))) }
   cat("\n")

   old_sens <- P$sens
   cmp <- data.frame(setting = sens$setting,
                     MSPE_new = sens$MSPE,  MSPE_old = old_sens$MSPE[match(sens$setting, old_sens$setting)],
                     edges_new = sens$edges, edges_old = old_sens$edges[match(sens$setting, old_sens$setting)],
                     prexpl_new = sens$pr_expl, prexpl_old = old_sens$pr_expl[match(sens$setting, old_sens$setting)])
   cat("\n=== sensitivity: new (100 draws) vs old (20 draws) ===\n"); print(cmp, row.names = FALSE)
   cat(sprintf("MSPE range new [%.4f, %.4f] | old [%.4f, %.4f]\n",
               min(sens$MSPE), max(sens$MSPE), min(old_sens$MSPE), max(old_sens$MSPE)))
   cat(sprintf("edges new %d to %d | old %d to %d\n", min(sens$edges), max(sens$edges),
               min(old_sens$edges), max(old_sens$edges)))
   cat(sprintf("edge counts reproduce the stored sweep exactly: %s\n",
               identical(as.integer(cmp$edges_new), as.integer(cmp$edges_old))))

   # ---- seed sweep ---------------------------------------------------------
   cat("\n-- seeds --\n"); flush.console()
   sd_tab <- NULL
   for (s in SEEDS6){ cat(sprintf("\r   seed=%d ", s)); flush.console()
      f <- fit_v(Y_train, C$K, C$RANK, s, C$NRUN, C$BURN, C$THIN, C$JIT)
      Om <- omega_hat(f); store[[sprintf("seed%d", s)]] <- Om
      sd_tab <- rbind(sd_tab, data.frame(seed = s, mspe = round(ev(Om),4), edges = ci_n(f))) }
   cat("\n")
   os <- P$seed_res
   sc <- data.frame(seed = sd_tab$seed, mspe_new = sd_tab$mspe, mspe_old = round(os[,1],4),
                    edges_new = sd_tab$edges, edges_old = os[,2])
   cat("\n=== seed sweep: new vs old ===\n"); print(sc, row.names = FALSE)
   cat(sprintf("new: mean %.4f sd %.4f range [%.4f, %.4f] | old: mean %.4f sd %.4f range [%.4f, %.4f]\n",
               mean(sd_tab$mspe), sd(sd_tab$mspe), min(sd_tab$mspe), max(sd_tab$mspe),
               mean(os[,1]), sd(os[,1]), min(os[,1]), max(os[,1])))
   cat(sprintf("edge counts reproduce the stored sweep exactly: %s\n",
               identical(as.integer(sc$edges_new), as.integer(sc$edges_old))))

   # ---- p29 rank / k justification ----------------------------------------
   just <- NULL
   if (sp$justification){
      cat("\n-- rank / k justification --\n"); flush.console()
      for (rk in list("full", 2, 3)) for (kk in sp$K_GRID) for (sd_i in c(33,34,35)){
         lab <- if (identical(rk,"full")) "full" else as.character(rk)
         cat(sprintf("\r   rank=%-4s k=%-2d seed=%d ", lab, kk, sd_i)); flush.console()
         f <- fit_v(Y_train, kk, rk, sd_i, C$NRUN, C$BURN, C$THIN, C$JIT)
         Om <- omega_hat(f); store[[sprintf("just_%s_k%d_s%d", lab, kk, sd_i)]] <- Om
         just <- rbind(just, data.frame(rank = lab, k = kk, seed = sd_i,
                                        mspe1 = round(ev(Om),4), ci_edges = ci_n(f),
                                        pr_expl = round(pr_expl(f, p),3))) }
      cat("\n")
      oj <- P$justification
      key <- function(d) paste(d$rank, d$k, d$seed, sep="|")
      m <- match(key(just), key(oj))
      just$mspe1_old <- oj$mspe1[m]; just$ci_edges_old <- oj$ci_edges[m]
      just$pr_expl_old <- oj$pr_expl[m]
      cat(sprintf("\nedge counts reproduce the stored table exactly: %s\n",
                  identical(as.integer(just$ci_edges), as.integer(just$ci_edges_old))))
      agg_new <- aggregate(cbind(mspe1, ci_edges) ~ rank + k, data = just, FUN = mean)
      agg_old <- aggregate(cbind(mspe1_old, ci_edges_old) ~ rank + k, data = just, FUN = mean)
      ag <- merge(agg_new, agg_old, by = c("rank","k"))
      ag <- ag[order(match(ag$rank, c("full","2","3")), ag$k), ]
      ag$mspe1 <- round(ag$mspe1, 4); ag$mspe1_old <- round(ag$mspe1_old, 4)
      cat("\n=== rank/k table, averaged over seeds 33-35: new vs old ===\n")
      cat("    (PPC columns are carried over from the stored table, not recomputed)\n")
      print(ag, row.names = FALSE) }

   saveRDS(list(tag = tg, p = p, n_draws = N_DRAWS,
                hid_seeds = HID_BASE + seq_len(N_DRAWS),
                sensitivity = sens, sensitivity_vs_old = cmp,
                seeds = sd_tab, seeds_vs_old = sc,
                justification = just,
                omegas = store,            # every sweep fit's Omega, saved this time
                config = C, run_at = format(Sys.time())),
           file.path(APP, sp$out))
   cat(sprintf("\nSaved %s (%d Omegas stored)\n", sp$out, length(store)))
}
