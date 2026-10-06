# ============================================================================
#  Which single precision matrix does each Bayesian method hand to the MSPE
#  harness, and how much does that choice matter?
#
#  (a) mean of per-draw precision      Omega_a = mean_s Omega_s
#  (b) inverse of mean covariance      Omega_b = ( mean_s Omega_s^{-1} )^{-1}
#
#  As used by the production scripts:
#    VAPOR              (a)  posterior_mean_omega(fit$Lambda, fit$residual_prec)
#    PrecFact (no VAR)  (a)  fact_prec_processing: Reduce("+", prec_mats)/length(.)
#    VAR + PrecFact     (a)  same
#    Kastner (VAR+FSV)  (b)  solve(apply(vcov(m)[1,,,], c(1,2), mean))
#  so the published MSPE table mixes conventions: Kastner alone uses (b).
#
#  NO PRODUCTION FIT STORES ITS DRAWS -- only the single Omega it produced --
#  so neither convention can be recovered for the other side without refitting.
#  Every production fit is therefore refit here with its exact seed and
#  settings, accumulating BOTH sums. The script checks that each method's own
#  convention reproduces the stored Omega, which is what licenses the rest.
#  Nothing saved is overwritten; results go to p<tag>_convention_audit.rds.
#
#  Evaluation: the same 100 held-out sets, seeds 90001..90100.
# ============================================================================
suppressMessages({library(VAPOR); library(PrecFact); library(bayesianVARs); library(MASS)})
APP <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec/application")
N_DRAWS <- 100L; HID_BASE <- 90000L; HID_FRAC <- 0.3

predict_held_out <- function(Om, y_obs, obs, hid)
   as.vector(-solve(Om[hid,hid]) %*% Om[hid,obs] %*% y_obs)
make_hid_sets <- function(base, n, p, Tt, frac){ nh <- floor(frac*p)
   lapply(seq_len(n), function(b){ set.seed(base+b)
      lapply(seq_len(Tt), function(t) sort(sample(p, nh))) }) }
#  per-quarter loss vector, averaged over the held-out draws
qloss <- function(Om, hid, Y_test, p)
   rowMeans(vapply(hid, function(hl) vapply(seq_along(hl), function(t){
      h <- hl[[t]]; o <- setdiff(seq_len(p), h); yt <- Y_test[t,]
      mean((yt[h] - predict_held_out(Om, yt[o], o, h))^2) }, numeric(1)),
      numeric(nrow(Y_test))))

#  accumulate both sums from a draw accessor
both_sums <- function(omega_at, nd, p){
   os <- matrix(0,p,p); cs <- matrix(0,p,p)
   for (s in seq_len(nd)){ Om <- omega_at(s); os <- os + Om; cs <- cs + solve(Om) }
   list(a = os/nd, b = solve(cs/nd), n = nd) }

for (tg in c("p29","p100")){
   R <- readRDS(file.path(APP, sprintf("%s_results.rds", tg))); C <- R$config
   Y_train <- R$Y_train; Y_test <- R$Y_test; p <- R$p
   residrr <- scale(R$first_stage$resid_rr)
   hid <- make_hid_sets(HID_BASE, N_DRAWS, p, nrow(Y_test), HID_FRAC)
   seeds <- C$SEED + 1:3
   cat(sprintf("\n\n################  %s  (p = %d)  ################\n", tg, p))
   cat(sprintf("refitting the 12 production fits (seeds %s) to recover the draws\n",
               paste(seeds, collapse = ", ")))

   fits <- list()
   # ---- VAPOR -------------------------------------------------------------
   fits[["VAPOR"]] <- lapply(seeds, function(s){
      f <- vapor(Y_train, lags=1, k=C$K, rank=C$RANK, prior_B="dl", nrun=C$NRUN,
                 burn=C$BURN, thin=C$THIN, seed=s, jitter=C$JIT, init="random")
      z <- both_sums(function(i){ Om <- tcrossprod(f$Lambda[,,i])
         diag(Om) <- diag(Om) + f$residual_prec[i,]; Om }, dim(f$Lambda)[3], p)
      cat(sprintf("   VAPOR seed %d done\n", s)); flush.console(); z })

   # ---- PrecFact rows: run_pf seeds the first attempt as seed0 + 7919 -----
   pf_one <- function(Ymat, seed0){
      set.seed(seed0 + 7919)
      res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=C$K, nrun=C$NRUN, thin=C$THIN,
                      Y=Ymat, alldiff=0, a_dir=0.1, b_dir=0.1)
      keep <- (C$PF_BURNOUT + 1L):(C$NRUN/C$THIN)
      both_sums(function(i){ s <- keep[i]
         Om <- tcrossprod(res$Lambda[,,s]); diag(Om) <- diag(Om) + res$residuals[s,]; Om },
         length(keep), p) }
   fits[["PrecFact (no VAR)"]] <- lapply(seeds, function(s){
      z <- pf_one(Y_train, s); cat(sprintf("   PrecFact seed %d done\n", s)); flush.console(); z })
   fits[["VAR + PrecFact"]] <- lapply(seeds, function(s){
      z <- pf_one(residrr, s); cat(sprintf("   VAR+PrecFact seed %d done\n", s)); flush.console(); z })

   # ---- Kastner -----------------------------------------------------------
   fits[["Kastner (VAR+FSV)"]] <- lapply(seeds, function(s){
      set.seed(s)
      Yj <- Y_train + matrix(rnorm(prod(dim(Y_train)), 0, C$VAR_JIT), nrow(Y_train), p)
      sp <- specify_prior_sigma(M=p, type="factor", factor_factors=as.integer(C$K),
                                factor_priorfacloadtype="colwiseng",
                                factor_heteroskedastic=TRUE, quiet=TRUE)
      m <- bvar(Yj, lags=1L, draws=C$KAST_DRAWS, burnin=C$KAST_BURNIN,
                thin=C$KAST_THIN, prior_sigma=sp, quiet=TRUE)
      V <- vcov(m); nd <- dim(V)[4]
      os <- matrix(0,p,p); cs <- matrix(0,p,p)
      for (i in seq_len(nd)){ Sg <- V[1,,,i]; cs <- cs + Sg; os <- os + solve(Sg) }
      rm(m, V); gc(verbose=FALSE)
      cat(sprintf("   Kastner seed %d done\n", s)); flush.console()
      list(a = os/nd, b = solve(cs/nd), n = nd) })

   # ---- does each method's own convention reproduce the stored Omega? -----
   own <- c(VAPOR="a", "PrecFact (no VAR)"="a", "VAR + PrecFact"="a", "Kastner (VAR+FSV)"="b")
   stored_key <- c(VAPOR="VAPOR (rank 3)", "PrecFact (no VAR)"="PrecFact (no VAR)",
                   "VAR + PrecFact"="VAR + PrecFact", "Kastner (VAR+FSV)"="Kastner (VAR+FSV)")
   chk <- do.call(rbind, lapply(names(fits), function(nm){
      st <- R$methods[[stored_key[[nm]]]]$omega
      d <- vapply(seq_along(st), function(i)
         max(abs(fits[[nm]][[i]][[own[[nm]]]] - st[[i]])), numeric(1))
      data.frame(Method = nm, own_convention = own[[nm]],
                 max_abs_diff_vs_stored = signif(max(d), 3)) }))
   cat("\n--- refit check: each method's own convention vs the stored Omega ---\n")
   print(chk, row.names = FALSE)

   # ---- tables under each convention --------------------------------------
   for (cv in c("a","b")){
      ql <- lapply(fits, function(fs) lapply(fs, function(z) qloss(z[[cv]], hid, Y_test, p)))
      mspe <- lapply(ql, function(L) vapply(L, mean, numeric(1)))
      lv <- Reduce(`+`, ql[["VAPOR"]])/3
      dm <- function(lc){ d <- lc - lv; n <- length(d); se <- sd(d)/sqrt(n)
         s <- (mean(d)/se)*sqrt((n-1)/n)
         c(Diff = round(mean(d),3), SE = round(se,3), DM = round(s,2),
           p = round(2*pt(-abs(s), df=n-1),3), Won = sum(d > 0)) }
      tab <- do.call(rbind, lapply(names(fits), function(nm){
         z <- if (nm == "VAPOR") c(Diff=NA,SE=NA,DM=NA,p=NA,Won=NA)
              else dm(Reduce(`+`, ql[[nm]])/3)
         data.frame(Method = nm,
                    MSPE = sprintf("%.4f (%.4f)", mean(mspe[[nm]]), sd(mspe[[nm]])),
                    Diff = z["Diff"], SE = z["SE"], DM = z["DM"], p = z["p"],
                    Won = if (is.na(z["Won"])) "--" else sprintf("%d / %d", z["Won"], nrow(Y_test))) }))
      cat(sprintf("\n=== %s, convention (%s): %s ===\n", tg, cv,
                  if (cv=="a") "mean of per-draw precision" else "inverse of mean covariance"))
      print(tab, row.names = FALSE)
      assign(paste0("tab_", cv), tab) }

   saveRDS(list(tag=tg, p=p, conventions=fits, check=chk,
                table_a=tab_a, table_b=tab_b, seeds=seeds,
                hid_seeds=HID_BASE+seq_len(N_DRAWS), run_at=format(Sys.time())),
           file.path(APP, sprintf("%s_convention_audit.rds", tg)))
   cat(sprintf("\nSaved %s_convention_audit.rds\n", tg))
}
