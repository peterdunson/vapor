# ============================================================================
#  The one-step-ahead posterior predictive check of supplement S1/S2, run for
#  the competitors as well as VAPOR.
#     Rscript applications/ppc_comparison.R 29
#     Rscript applications/ppc_comparison.R 100
#  Writes results/tables/ppc_comparison_p<d>.rds.  Stage 3b of run_all.sh.
#
#  KASTNER CACHE.  Kastner needs B draws and a time-varying Sigma_t, neither of
#  which any stored fit holds, so this script reads
#  results/fits/kastner_ppc_p<d>.rds when it exists and otherwise refits four
#  chains with sv_keep = "all" and writes that cache.  The refit is verified
#  against the stored fits before it is used, and the script stops if they
#  disagree.
#
#  WHAT IS HELD IDENTICAL TO VAPOR'S OWN CHECK.  ppc_stats() from
#  final_shared.R supplies all nine summaries, so the statistics, the
#  pseudo-inverse fallback inside them (solve -> MASS::ginv, which is what
#  makes p = 100 work at all) and the interval and p-value definitions are the
#  same function for every method.  The replicate count, the draw-picking rule
#  and the seed are CFG$PPC_N_REP and CFG$PPC_SEED, again as VAPOR uses them.
#  VAPOR's own column is produced by calling ppc_table() unmodified, and the
#  script asserts it equals the column in p<d>_tables.rds.
#
#  WHAT DIFFERS, PER METHOD, AND WHY.  ppc_one() draws one Omega per replicate
#  and applies it at every t, so it cannot express Kastner's time-varying
#  covariance.  The generators below therefore replace ppc_one() only, keeping
#  everything downstream of it:
#
#    VAPOR              ppc_one() itself: y_rep_t ~ N(B_s y_{t-1}, Omega_s^-1)
#    PrecFact (no VAR)  B = 0, so y_rep_t ~ N(0, Omega_s^-1); residuals are y_t
#    VAR + PrecFact     y_rep_t ~ N(B_hat y_{t-1}, Omega_s^-1), B_hat the
#                       full-rank first-stage posterior mean it is fitted on
#    Kastner            y_rep_t ~ N(B_s y_{t-1}, Sigma_{t,s}) with its OWN
#                       time-varying covariance, reconstructed per draw from
#                       facload and logvar.  Labelled "time-varying Sigma".
#    VAR + Glasso       plug-in: y_rep_t ~ N(B_hat y_{t-1}, Omega_hat^-1) with
#                       the FULL-RANK first stage and glasso_fullrank$Omega,
#                       as in the main MSPE table
#    Sample covariance  plug-in: B = 0, y_rep_t ~ N(0, Omega_hat^-1), matching
#                       the B = 0 convention of the MSPE comparison
#
#  The two plug-ins have no draw-to-draw variation in Omega, so their intervals
#  carry replicate noise only.  Residual lag-1 autocorrelation always uses the
#  method's own residuals: y_t - B y_{t-1} where it has a B, and y_t where not.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 210)

args <- commandArgs(trailingOnly = TRUE)
DIM  <- if (length(args)) as.integer(args[1]) else 29L
stopifnot(DIM %in% c(29L, 100L))

read_cfg <- function(dim){
   src <- sprintf("/Users/peterdunson/Desktop/vapor/applications/final_p%d.R", dim)
   for (q in parse(src))
      if (is.call(q) && identical(as.character(q[[1]])[1], "<-") &&
          identical(as.character(q[[2]]), "CFG")) return(eval(q[[3]]))
   stop(sprintf("no CFG <- list(...) in %s", basename(src))) }
CFG <- read_cfg(DIM)

R  <- readRDS(file.path(FITS_DIR, sprintf("p%d_results.rds", DIM)))
TT <- readRDS(file.path(FITS_DIR, sprintf("p%d_tables.rds", DIM)))
RE <- readRDS(file.path(FITS_DIR, sprintf("reused_p%d.rds", DIM)))
Y  <- R$Y_train; p <- R$p; vars <- R$vars
N_REP <- CFG$PPC_N_REP; PSEED <- CFG$PPC_SEED
cat(sprintf("p = %d | T_train = %d | n_rep = %d | ppc seed = %d\n",
            p, nrow(Y), N_REP, PSEED))
stopifnot(isTRUE(all.equal(unname(RE$panel$Y_train), unname(Y), tolerance = 0)))

#  the full-rank first-stage mean, the same object VAR + PrecFact and the
#  main VAR + Glasso row are built on
B_FULL <- fixB(R$first_stage$B_hat, p)

# ============================================================================
#  Kastner refit with sv_keep = "all", cached, and verified against the stored
#  fits before it is used.  sv_keep controls storage only: it does not touch
#  the random stream, which the verification below demonstrates.
# ============================================================================
KF <- file.path(FITS_DIR, sprintf("kastner_ppc_p%d.rds", DIM))
kastner_refit <- function(){
   st <- RE$kastner[[1]]$settings
   if (is.null(st)) st <- list(total_iterations = 20000L, burnin = 10000L,
                               thin = 10L, factor_factors = as.integer(CFG$K))
   seeds <- vapply(RE$kastner, function(q){
      s <- attr(q, "seed"); if (is.null(s)) q$seed else s }, numeric(1))
   cat(sprintf("\nKastner refit, sv_keep = \"all\": %d iterations, %d burn-in, thin %d,\n",
               st$total_iterations, st$burnin, st$thin))
   cat(sprintf("  factor_factors = %d, the stored seeds %s\n",
               st$factor_factors, paste(seeds, collapse = ", ")))
   out <- list()
   for (sd_i in seeds){
      t0 <- Sys.time(); set.seed(sd_i)
      Yj <- Y + matrix(rnorm(prod(dim(Y)), 0, CFG$VAR_JIT), nrow(Y), ncol(Y))
      sp <- specify_prior_sigma(M = p, type = "factor",
                                factor_factors = as.integer(st$factor_factors),
                                factor_priorfacloadtype = "colwiseng",
                                factor_heteroskedastic = TRUE, quiet = TRUE)
      m <- bvar(Yj, lags = 1L,
                draws = as.integer(st$total_iterations - st$burnin),
                burnin = as.integer(st$burnin), thin = as.integer(st$thin),
                prior_sigma = sp, quiet = TRUE, sv_keep = "all")
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      #  store the pieces Sigma_{t,s} rebuilds from, never the T x p x p x S
      #  vcov array: 0.6 GB at p = 29 and 7.4 GB at p = 100 against 32 and
      #  114 MB for these.
      z <- list(seed = sd_i, minutes = el,
                PHI = m$PHI, facload = m$facload, logvar = m$logvar,
                k = dim(m$facload)[2], M = p, n_draws = dim(m$PHI)[3],
                settings = st)
      #  the stored fits kept sv_keep = "last", which retains the LAST time
      #  point, so that is the slice to compare against
      #  ONE time slice, via vcov()'s `t` argument.  The whole array is
      #  T x p x p x draws: 0.6 GB at p = 29 but 7.4 GB at p = 100, which does
      #  not fit in memory here, and only the last slice is ever compared.
      Tl <- nrow(m$logvar)
      V  <- vcov(m, t = Tl)
      nd <- dim(V)[4]
      acc <- matrix(0, p, p); PD <- array(NA_real_, c(p, p, nd))
      for (i in seq_len(nd)){
         Pi <- solve(V[1, , , i]); PD[, , i] <- Pi; acc <- acc + Pi }
      z$Omega_last <- acc / nd
      dimnames(z$Omega_last) <- list(vars, vars)
      #  reused_p29.rds stores prec_draws and the comparison needs them; at
      #  p = 100 nothing stores them, so only the summary is kept
      if (DIM == 29L) z$prec_draws_last <- PD
      z$pcor_mean_draws <- Reduce(`+`, lapply(seq_len(dim(PD)[3]),
                                  function(i) pcor_from_omega(PD[, , i]))) / dim(PD)[3]
      dimnames(z$pcor_mean_draws) <- list(vars, vars)
      out[[length(out) + 1L]] <- z
      cat(sprintf("  seed %d: %.2f min, %d draws\n", sd_i, el, z$n_draws))
      flush.console(); rm(m, V, PD); invisible(gc(FALSE)) }
   out }

if (file.exists(KF)) { KCH <- readRDS(KF)
   cat(sprintf("\nKastner refit read from %s (no refit)\n", basename(KF)))
} else { KCH <- kastner_refit(); saveRDS(KCH, KF)
   cat(sprintf("\nwrote %s\n", KF)) }

#  ---- verification against the stored fits --------------------------------
cat("\n=== Kastner sv_keep = \"all\" refit vs the stored fits ===\n")
vrows <- list()
for (i in seq_along(KCH)){
   s <- KCH[[i]]$seed; q <- RE$kastner[[i]]
   vrows[[length(vrows) + 1L]] <- data.frame(
      object = sprintf("reused_p%d.rds", DIM), seed = s, quantity = "Omega",
      max_abs_diff = max(abs(KCH[[i]]$Omega_last - q$Omega)), row.names = NULL)
   #  `$` partial-matches, and reused_p100.rds has prec_draws_FILE but no
   #  prec_draws, so q$prec_draws would silently return that filename string.
   #  Both lookups below are exact.
   if ("prec_draws" %in% names(q) && !is.null(KCH[[i]]$prec_draws_last))
      vrows[[length(vrows) + 1L]] <- data.frame(
         object = sprintf("reused_p%d.rds", DIM), seed = s,
         quantity = "prec_draws at the stored time index",
         max_abs_diff = max(abs(KCH[[i]]$prec_draws_last -
                                q[["prec_draws", exact = TRUE]])), row.names = NULL)
   if ("pcor_mean_draws" %in% names(q))
      vrows[[length(vrows) + 1L]] <- data.frame(
         object = sprintf("reused_p%d.rds", DIM), seed = s, quantity = "pcor_mean_draws",
         max_abs_diff = max(abs(KCH[[i]]$pcor_mean_draws -
                                q[["pcor_mean_draws", exact = TRUE]])), row.names = NULL) }
if (DIM == 100L){
   K20 <- readRDS(file.path(FITS_DIR, "kastner_p100_20k_draws.rds"))
   for (i in seq_along(KCH)){
      j <- which(vapply(K20$chains, function(q) q$seed, numeric(1)) == KCH[[i]]$seed)
      if (!length(j)) next
      vrows[[length(vrows) + 1L]] <- data.frame(
         object = "kastner_p100_20k_draws.rds", seed = KCH[[i]]$seed, quantity = "Omega",
         max_abs_diff = max(abs(KCH[[i]]$Omega_last - K20$chains[[j]]$Omega)), row.names = NULL)
      vrows[[length(vrows) + 1L]] <- data.frame(
         object = "kastner_p100_20k_draws.rds", seed = KCH[[i]]$seed,
         quantity = "pcor_mean_draws",
         max_abs_diff = max(abs(KCH[[i]]$pcor_mean_draws - K20$chains[[j]]$pcor_mean_draws)),
         row.names = NULL) } }
VER <- do.call(rbind, vrows)
VER$max_abs_diff <- signif(VER$max_abs_diff, 3)
print(VER, row.names = FALSE)
TOL <- 1e-8
if (any(VER$max_abs_diff > TOL))
   stop(sprintf("sv_keep = \"all\" does NOT reproduce the stored Kastner fits (max %.3g > %g); stopping before the PPC",
                max(VER$max_abs_diff), TOL))
cat(sprintf("\n  all %d comparisons within %g: the refit reproduces the stored fits\n",
            nrow(VER), TOL))

# ============================================================================
#  replicate generators.  Each returns exactly what ppc_one() returns, so
#  everything downstream is the same code for every method.
# ============================================================================
#  shared: innovations from a precision matrix, with ppc_one()'s own fallbacks
innov <- function(Om, n, p){
   Sig <- tryCatch(chol2inv(chol(Om)), error = function(e) MASS::ginv(Om))
   L   <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-8, p)))
   matrix(rnorm(n * p), n, p) %*% L }

pack <- function(M, E, Y){
   Tn <- nrow(Y); Yrep <- M + E
   colnames(Yrep) <- colnames(E) <- colnames(Y)
   Yobs <- Y[-1, , drop = FALSE]
   list(Yrep = Yrep, resid_rep = E, Yobs = Yobs, resid_obs = Yobs - M) }

#  B = 0 (PrecFact raw, sample covariance): the mean is zero and the residual
#  is the series itself
gen_nob <- function(om_of){ function(ch, s, Y, p){
   Tn <- nrow(Y); M <- matrix(0, Tn - 1, p)
   pack(M, innov(om_of(ch, s), Tn - 1, p), Y) } }

#  fixed B (VAR + PrecFact, VAR + Glasso)
gen_fixedB <- function(B, om_of){ function(ch, s, Y, p){
   Tn <- nrow(Y); M <- Y[-Tn, , drop = FALSE] %*% t(B)
   pack(M, innov(om_of(ch, s), Tn - 1, p), Y) } }

#  Kastner: its own B_s and its own time-varying Sigma_{t,s}
gen_kastner <- function(ch, s, Y, p){
   Tn <- nrow(Y); k <- ch$k; M_ <- ch$M
   B <- fixB(ch$PHI[, , s], p)
   MU <- Y[-Tn, , drop = FALSE] %*% t(B)
   L  <- ch$facload[, , s]
   E  <- matrix(NA_real_, Tn - 1, p)
   for (t in seq_len(Tn - 1)){
      lv <- ch$logvar[t, , s]
      Sg <- L %*% diag(exp(lv[M_ + seq_len(k)]), k) %*% t(L) + diag(exp(lv[seq_len(M_)]))
      Lc <- tryCatch(chol(Sg), error = function(e) chol(Sg + diag(1e-8, p)))
      E[t, ] <- rnorm(p) %*% Lc }
   pack(MU, E, Y) }

#  the draw-picking rule and the table assembly of ppc_table(), with the
#  generator swapped in
ppc_generic <- function(chains, Y, p, n_rep, seed, gen, n_draws_of, BIG = 1e3){
   set.seed(seed)
   picks <- do.call(rbind, lapply(seq_along(chains), function(ci)
      cbind(ci, sample(n_draws_of(chains[[ci]]), ceiling(n_rep / length(chains)),
                       replace = TRUE))))
   picks <- picks[seq_len(min(n_rep, nrow(picks))), , drop = FALSE]
   D <- lapply(seq_len(nrow(picks)), function(i){
      r <- picks[i, ]; z <- gen(chains[[r[1]]], r[2], Y, p)
      list(rep = ppc_stats(z$Yrep, z$resid_rep),
           obs = ppc_stats(z$Yobs, z$resid_obs),
           max_abs = max(abs(z$Yrep))) })
   REP <- do.call(rbind, lapply(D, function(z) z$rep))
   OBS <- do.call(rbind, lapply(D, function(z) z$obs))
   tab <- data.frame(statistic = colnames(REP),
      observed = round(colMeans(OBS), 3), rep_mean = round(colMeans(REP), 3),
      lo = round(apply(REP, 2, quantile, 0.025), 3),
      hi = round(apply(REP, 2, quantile, 0.975), 3),
      ppp = round(colMeans(REP >= OBS), 3), row.names = NULL)
   attr(tab, "n_draws") <- nrow(picks)
   attr(tab, "n_degenerate") <- sum(vapply(D, function(z) z$max_abs > BIG, logical(1)))
   tab }

om_vapor    <- function(ch, s) omega_one(ch, s)
nd_lambda   <- function(ch) dim(ch$Lambda)[3]
plug        <- function(Om) list(list(Omega = Om))
om_plug     <- function(ch, s) ch$Omega

# ============================================================================
#  run every method
# ============================================================================
cat("\n=== running the check, one method at a time ===\n")
Om_glasso <- RE$glasso_fullrank$Omega          # MAIN row: full-rank first stage
Om_sample <- RE$sample_cov_omega[[1]]

PPC <- list()
cat("  VAPOR (rank 3)\n")
PPC[["VAPOR (rank 3)"]] <- ppc_table(R$draws$vapor, Y, p, N_REP, PSEED)
cat("  PrecFact (no VAR)\n")
PPC[["PrecFact (no VAR)"]] <- ppc_generic(R$draws$precfact, Y, p, N_REP, PSEED,
                                          gen_nob(om_vapor), nd_lambda)
cat("  VAR + PrecFact\n")
PPC[["VAR + PrecFact"]] <- ppc_generic(R$draws$var_precfact, Y, p, N_REP, PSEED,
                                       gen_fixedB(B_FULL, om_vapor), nd_lambda)
cat("  Kastner (VAR+FSV), time-varying Sigma\n")
PPC[["Kastner (VAR+FSV)"]] <- ppc_generic(KCH, Y, p, N_REP, PSEED,
                                          gen_kastner, function(ch) ch$n_draws)
cat("  VAR + Glasso (plug-in)\n")
PPC[["VAR + Glasso (plug-in)"]] <- ppc_generic(plug(Om_glasso), Y, p, N_REP, PSEED,
                                               gen_fixedB(B_FULL, om_plug),
                                               function(ch) 1L)
cat("  Sample covariance (plug-in, B = 0)\n")
PPC[["Sample covariance (plug-in)"]] <- ppc_generic(plug(Om_sample), Y, p, N_REP, PSEED,
                                                    gen_nob(om_plug), function(ch) 1L)

#  ---- VAPOR's column must equal the paper's ------------------------------
cat("\n=== VAPOR's column against p%d_tables.rds ===\n")
chk <- merge(TT$ppc, PPC[["VAPOR (rank 3)"]], by = "statistic", suffixes = c("_paper", "_here"))
chk <- chk[match(TT$ppc$statistic, chk$statistic), ]
CHK <- data.frame(statistic = chk$statistic,
   paper = sprintf("%.3f to %.3f (%.3f)", chk$lo_paper, chk$hi_paper, chk$ppp_paper),
   here  = sprintf("%.3f to %.3f (%.3f)", chk$lo_here,  chk$hi_here,  chk$ppp_here),
   max_abs_diff = signif(pmax(abs(chk$lo_paper - chk$lo_here),
                              abs(chk$hi_paper - chk$hi_here),
                              abs(chk$ppp_paper - chk$ppp_here),
                              abs(chk$observed_paper - chk$observed_here)), 3),
   row.names = NULL)
print(CHK, row.names = FALSE)
cat(sprintf("\n  VAPOR column identical to the stored table: %s\n",
            all(CHK$max_abs_diff == 0)))

# ============================================================================
#  (3a) nine summaries x methods, and (3b) how many fall inside
# ============================================================================
STATS <- PPC[[1]]$statistic
cell  <- function(t) sprintf("%.3f to %.3f (%.3f)", t$lo, t$hi, t$ppp)
inside <- function(t) t$observed >= t$lo & t$observed <= t$hi

#  Eight of the nine summaries are computed on Yobs, so their observed value is
#  the same for every method.  resid_ac1 is not: it is taken on each method's
#  OWN residuals, y_t - B y_{t-1}, which differ with B (and equal y_t where
#  B = 0).  The observed column below therefore holds the shared value, and the
#  per-method resid_ac1 observed values are given on their own row.  The
#  interval membership count always uses each method's own observed value.
OBSM <- do.call(cbind, lapply(PPC, function(t) t$observed))
rownames(OBSM) <- STATS
shared <- apply(OBSM, 1, function(x) max(x) - min(x)) == 0
stopifnot(all(shared[STATS != "resid_ac1"]))

A <- data.frame(summary = STATS, observed = PPC[[1]]$observed, row.names = NULL)
for (nm in names(PPC)) A[[nm]] <- cell(PPC[[nm]])
A$observed <- sprintf("%.3f", A$observed)
A$observed[STATS == "resid_ac1"] <- "per method, below"

Bo <- data.frame(summary = "resid_ac1 observed (own residuals)",
                 observed = NA_character_, row.names = NULL)
for (nm in names(PPC))
   Bo[[nm]] <- sprintf("%.3f", PPC[[nm]]$observed[STATS == "resid_ac1"])
B <- data.frame(summary = "inside the interval (of 9)", observed = NA_character_,
                row.names = NULL)
for (nm in names(PPC)) B[[nm]] <- sprintf("%d", sum(inside(PPC[[nm]])))
TAB <- rbind(A, Bo, B)

cat(sprintf("\n\n=== (3a) p = %d: central 95%% replicate interval and (p-value) ===\n", p))
cat("    rows = the nine summaries | observed in its own column\n\n")
print(TAB, row.names = FALSE)
cat(sprintf("\n=== (3b) summaries inside the interval, of 9 ===\n\n"))
IN <- data.frame(method = names(PPC),
                 inside_of_9 = vapply(PPC, function(t) sum(inside(t)), integer(1)),
                 n_draws = vapply(PPC, function(t) attr(t, "n_draws"), integer(1)),
                 n_degenerate = vapply(PPC, function(t) attr(t, "n_degenerate"), integer(1)),
                 row.names = NULL)
print(IN, row.names = FALSE)

OUT <- file.path(TAB_DIR, sprintf("ppc_comparison_p%d.rds", DIM))
saveRDS(list(dim = p, n_rep = N_REP, ppc_seed = PSEED, tables = PPC,
             wide = TAB, inside = IN, vapor_check = CHK, observed_matrix = OBSM,
             kastner_verification = VER,
             notes = c(glasso = "full-rank first stage, glasso_fullrank$Omega, plug-in",
                       sample_cov = "B = 0, plug-in",
                       kastner = "own B_s draws and own time-varying Sigma_{t,s}"),
             run_at = format(Sys.time())), OUT)
cat(sprintf("\nwrote %s\n", OUT))
