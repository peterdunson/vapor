# ============================================================================
#  One (p, T) setting of sim_f, in its own process.
#     Rscript simulation/run_sim_f_setting.R <p> <T>
#
#  sim_f.R re-seeds from the replicate index alone -- set.seed(sd_i) before the
#  DGP, before each vapor() call and before the Kastner jitter -- so a setting
#  run on its own reproduces the sequential result exactly.  This script
#  evaluates sim_f.R's definitions (everything before its first top-level for
#  loop), then runs the replicate loop for one setting only and writes
#  simulation/results/sim_f_parts/sim_f_part_p<p>_T<T>.rds.
# ============================================================================
suppressMessages({library(VAPOR); library(PrecFact); library(bayesianVARs)})
A <- commandArgs(trailingOnly = TRUE)
P <- as.integer(A[1]); TT <- as.integer(A[2])
SRC  <- "/Users/peterdunson/Desktop/vapor/simulation/sim_f.R"
PART <- sprintf("/Users/peterdunson/Desktop/vapor/simulation/results/sim_f_parts/sim_f_part_p%d_T%d.rds", P, TT)
t0 <- Sys.time()
cat(sprintf("VAPOR %s | as = %s, bs = %s | setting p = %d, T = %d | start %s\n",
            as.character(utils::packageVersion("VAPOR")),
            formals(vapor)$as, formals(vapor)$bs, P, TT,
            format(t0, "%Y-%m-%d %H:%M:%S"))); flush.console()
e <- new.env(parent = globalenv())
for (q in parse(SRC)){
   if (is.call(q) && identical(as.character(q[[1]])[1], "for")) break
   try(eval(q, envir = e), silent = TRUE) }
for (nm in ls(e)) assign(nm, get(nm, envir = e), envir = globalenv())
label <- sprintf("p=%d, T=%d", P, TT)
#  Resume support.  Each replicate is seeded from sd_i alone, so a replicate
#  computed in an earlier process is bit-identical to the same replicate
#  computed now; the per-replicate cache below therefore changes nothing about
#  the results, only how much has to be recomputed after an interruption.
CACHE <- sub("[.]rds$", "_seeds.rds", PART)
sm <- vector("list", N_SEEDS)
if (file.exists(CACHE)){
   old <- readRDS(CACHE)
   for (k in seq_len(N_SEEDS)) if (!is.null(old[[k]])) sm[[k]] <- old[[k]]
   cat(sprintf("resuming: %d of %d replicates already saved in %s\n",
               sum(!vapply(sm, is.null, logical(1))), N_SEEDS, basename(CACHE))) }
for (sd_i in seq_len(N_SEEDS)){
   if (!is.null(sm[[sd_i]])) next
   ts <- Sys.time()
   sim <- simulate_lowrank_calibrated(T = TT, p = P, seed = sd_i)
   sm[[sd_i]] <- fit_all(sim, P, seed = sd_i)$metrics
   saveRDS(sm, CACHE)                      # save after EVERY replicate
   cat(sprintf("  seed %d/%d done in %.1f min (%s)\n", sd_i, N_SEEDS,
               as.numeric(difftime(Sys.time(), ts, units = "mins")),
               format(Sys.time(), "%H:%M:%S"))); flush.console() }
methods <- rownames(sm[[1]]); cols <- colnames(sm[[1]])
mn <- outer(methods, cols, Vectorize(function(mm, cc)
   mean(sapply(sm, function(M) M[mm, cc]), na.rm = TRUE)))
sdv <- outer(methods, cols, Vectorize(function(mm, cc)
   sd(sapply(sm, function(M) M[mm, cc]), na.rm = TRUE)))
dimnames(mn) <- dimnames(sdv) <- list(methods, cols)
res <- data.frame(p = P, T = TT, TP_ratio = round(TT/P, 2),
                  block = "T<p sweep (T/p=0.7)", Method = methods,
                  Omega_Frob = mn[,"Frob"], Omega_Frob_sd = sdv[,"Frob"],
                  B_Frob = mn[,"B_Frob"], B_Frob_sd = sdv[,"B_Frob"],
                  Sens = mn[,"Sens"], Sens_sd = sdv[,"Sens"],
                  Spec = mn[,"Spec"], Spec_sd = sdv[,"Spec"],
                  row.names = NULL, stringsAsFactors = FALSE)
t1 <- Sys.time()
saveRDS(list(p = P, T = TT, label = label, results = res,
             paired = paired_diffs(sm, label), raw = sm,
             start = format(t0), finish = format(t1),
             minutes = as.numeric(difftime(t1, t0, units = "mins"))), PART)
cat(sprintf("finish %s | %.1f minutes | wrote %s\n",
            format(t1, "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(t1, t0, units = "mins")), basename(PART)))
