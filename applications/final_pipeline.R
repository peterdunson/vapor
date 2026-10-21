# ============================================================================
#  SHARED PIPELINE.  Sourced by final_p29.R and final_p100.R AFTER they define
#  CFG and ORDER; it sets no parameter of its own.
#
#  Order of work: data, held-out sets, first-stage VAR and its rank-3
#  truncation, the four stochastic methods (4 chains each, seeds 101-104, with
#  the replacement ladder), the three deterministic rows, then every table and
#  the two diagnostic figures.  Saves raw draws per chain.
# ============================================================================
stopifnot(exists("CFG"), exists("ORDER"))
t_start <- Sys.time()
cat(sprintf("\n%s\n", CFG$LABEL))
cat(sprintf("seeds %s | nrun %d burn %d thin %d | k %d rank %d | %d replacements allowed\n",
            paste(range(CFG$SEEDS), collapse="-"), CFG$NRUN, CFG$BURN, CFG$THIN,
            CFG$K, CFG$RANK, CFG$MAX_REPL))

# ---- data ------------------------------------------------------------------
D <- build_panel(CFG)
Y_train <- D$train; Y_test <- D$test; vars <- D$vars; p <- D$p
T_train <- nrow(Y_train); T_test <- nrow(Y_test)
Y_last_train <- Y_train[T_train, ]

#  ---- reused competitor fits, never refitted -------------------------------
RE_F <- file.path(FITS_DIR, sprintf("reused_p%d.rds", CFG$P))
stopifnot(file.exists(RE_F))
REUSED <- readRDS(RE_F)
stopifnot(isTRUE(all.equal(unname(REUSED$panel$Y_train), unname(Y_train), tolerance = 0)),
          isTRUE(all.equal(unname(REUSED$panel$Y_test),  unname(Y_test),  tolerance = 0)))
cat(sprintf("\nreused fits from %s, verified on this panel (%s):\n", basename(RE_F),
            REUSED$panel$train))
cat(sprintf("  Kastner         %d chain(s), seeds %s\n", length(REUSED$kastner),
            paste(vapply(REUSED$kastner, function(q){
               s <- attr(q, "seed"); if (is.null(s)) as.integer(q$seed) else as.integer(s) },
               integer(1)), collapse = ", ")))
cat(sprintf("  VAR + Glasso    full-rank rho %.4f | rank-3 rho %.4f\n",
            REUSED$glasso_fullrank$rho, REUSED$glasso_rank3$rho))
cat("  Sample covariance, Random walk: from the same object\n")

hid_sets <- make_hid_sets(CFG$HID_SEED_BASE, CFG$N_HID, p, T_test, CFG$HID_FRAC)
stopifnot(isTRUE(all.equal(hid_sets, REUSED$hid_sets)))
cat(sprintf("held-out evaluation: %d sets, seeds %d..%d, %d of %d hidden per quarter\n",
            CFG$N_HID, min(attr(hid_sets,"seeds")), max(attr(hid_sets,"seeds")),
            floor(CFG$HID_FRAC*p), p))

# ---- first-stage VAR and its rank-r truncation -----------------------------
cat(sprintf("\n=== first-stage VAR(1), %d draws, and its rank-%d truncation ===\n",
            CFG$VAR_DRAWS, CFG$RANK))
vr <- fit_var_mean(Y_train, seed = CFG$SEEDS[1], draws = CFG$VAR_DRAWS,
                   burnin = CFG$VAR_BURNIN, jit = CFG$VAR_JIT)
rt <- rank_truncate_B(vr$B_hat, p, CFG$RANK)
resid_rr   <- vr$Y - vr$X %*% rt$B          # rank-3 truncation, supplement rows
resid_full <- vr$Y - vr$X %*% vr$B_hat      # full-rank posterior mean, MAIN rows
cat(sprintf("singular values: %s ...\n", paste(sprintf("%.2f", head(rt$sv,6)), collapse=", ")))
cat(sprintf("rank %d keeps %.1f%% of the squared Frobenius norm | spectral radius %.3f\n",
            CFG$RANK, 100*rt$frob_kept, spec_rad(rt$B[seq_len(p),])))

# ---- the four stochastic methods, all chains concurrently -----------------
#  The chains are cached as soon as the last one returns.  They are by far the
#  most expensive part of the run, and an error further down -- in a table, a
#  figure, or a column name -- must never cost them again.  The cache is reused
#  only if it was written by the same settings and the same data.
cache_path <- file.path(CFG$OUT, sub("_final", "_chains_cache", CFG$OUT_RDS))
cache_key <- list(K=CFG$K, RANK=CFG$RANK, NRUN=CFG$NRUN, BURN=CFG$BURN,
                  THIN=CFG$THIN, JIT=CFG$JIT, SEEDS=CFG$SEEDS,
                  MAX_REPL=CFG$MAX_REPL, Y=digest_mat(Y_train))
#  USE_CACHE now DEFAULTS TO FALSE: the chain caches written by the earlier
#  runs are in applications/_archive/pre_base11/fits/ and are not read back,
#  so a run refits from scratch unless a config explicitly asks for the cache.
USE_CACHE <- if (is.null(CFG$USE_CACHE)) FALSE else isTRUE(CFG$USE_CACHE)
CA <- if (USE_CACHE && file.exists(cache_path)) readRDS(cache_path) else NULL
if (!USE_CACHE) cat("chain cache disabled (USE_CACHE = FALSE); chains will be refitted\n")
if (!is.null(CA) && !isTRUE(all.equal(CA$key, cache_key))){
   cat("chain cache present but its settings or data differ; refitting\n"); CA <- NULL }
if (!is.null(CA)){
   cat(sprintf("\n=== chains read from %s (fitted %s) ===\n",
               basename(cache_path), CA$run_at))
   CH <- CA$CH; LOG <- CA$LOG
} else {
   W <- n_workers(CFG$WORKERS)
   cat(sprintf("\n=== 4 refitted methods x %d seeds = %d chains, %d worker(s) ===\n",
               length(CFG$SEEDS), 4L*length(CFG$SEEDS), W))
   MC <- run_methods_parallel(method_closures(Y_train, resid_full, resid_rr, CFG),
                              CFG$SEEDS, CFG$MAX_REPL, W)
   CH <- list(vapor   = MC$chains$VAPOR,       pf_raw  = MC$chains$PrecFact,
              pf_rr   = MC$chains$`VAR+PrecFact`,
              pf_rr3  = MC$chains$`VAR+PrecFact (rank 3)`,
              kastner = REUSED$kastner)
   LOG <- MC$log
   saveRDS(list(key = cache_key, CH = CH, LOG = LOG, run_at = format(Sys.time())),
           cache_path)
   cat(sprintf("cached the chains in %s\n", basename(cache_path))) }
cat(sprintf("chains available: %s\n",
            paste(sprintf("%s %d", names(CH), vapply(CH, length, integer(1))),
                  collapse=" | ")))

chain_log <- do.call(rbind, lapply(LOG, function(r) data.frame(
   method = r$method, seed = r$seed, status = r$status,
   minutes = round(r$minutes, 2),
   error = if (is.null(r$error)) NA_character_ else r$error)))
cat("\n=== chain log ===\n"); print(chain_log, row.names = FALSE)
if (!length(CH$vapor)) stop("no VAPOR chain completed; nothing to report")

#  floor hits, per VAPOR chain
floor_tab <- do.call(rbind, lapply(CH$vapor, function(f){
   ix <- which(f$floor_hits > 0)
   if (!length(ix)) return(data.frame(seed=attr(f,"seed"), series="(none)",
                                      iterations_hit=0L, first_iteration=NA_integer_))
   data.frame(seed=attr(f,"seed"), series=names(f$floor_hits)[ix],
              iterations_hit=as.integer(f$floor_hits[ix]),
              first_iteration=as.integer(f$floor_first[ix])) }))
#  item 17: per chain, the number of ITERATIONS in which at least one series was
#  floored.  floor_iters is a new scalar returned by the sampler; it is not
#  recoverable from floor_hits, whose per-series counts cannot be intersected.
floor_iter_tab <- do.call(rbind, lapply(CH$vapor, function(f) data.frame(
   seed = attr(f, "seed"),
   series_floored = sum(f$floor_hits > 0),
   series_iteration_hits = sum(f$floor_hits),
   iterations_with_any_floor = if (is.null(f$floor_iters)) NA_integer_ else f$floor_iters,
   earliest_iteration = if (any(f$floor_first > 0)) min(f$floor_first[f$floor_first > 0]) else NA_integer_,
   nrun = f$settings$nrun, burn = f$settings$burn)))
cat("\n=== residual-precision floor, iterations with any series floored (item 17) ===\n")
print(floor_iter_tab, row.names = FALSE)

# ---- deterministic rows ----------------------------------------------------
cat("\n=== VAR + Glasso and sample covariance ===\n")
gl      <- REUSED$glasso_fullrank      # MAIN: full-rank first stage
gl_rr3  <- REUSED$glasso_rank3         # supplement: rank-3 first stage
cat(sprintf("glasso full rank: rho %.4f | %d nonzero off-diagonals\n",
            gl$rho, gl$n_nonzero))
cat(sprintf("glasso rank %d:    rho %.4f | %d nonzero off-diagonals\n",
            CFG$RANK, gl_rr3$rho, gl_rr3$n_nonzero))
Om_sample <- REUSED$sample_cov_omega[[1]]

# ---- one Omega list and one loss list per method --------------------------
OM <- list(`VAPOR (rank 3)`           = lapply(CH$vapor,  omega_hat),
           `PrecFact (no VAR)`        = lapply(CH$pf_raw, function(z) z$Omega),
           `VAR + PrecFact`           = lapply(CH$pf_rr,  function(z) z$Omega),
           `VAR + PrecFact (rank 3)`  = lapply(CH$pf_rr3, function(z) z$Omega),
           `Kastner (VAR+FSV)`        = lapply(CH$kastner,function(z) z$Omega),
           `VAR + Glasso`             = list(gl$Omega),
           `VAR + Glasso (rank 3)`    = list(gl_rr3$Omega),
           `Sample covariance`        = list(Om_sample))
LOSS <- lapply(OM, function(l) lapply(l, loss_avg, hid_sets, Y_test, p))
LOSS[["Random walk"]] <- list(loss_rw(hid_sets, Y_test, Y_last_train, p))

# ============================================================================
#  TABLES
# ============================================================================
#  VAPOR's MSPE against the config D value measured in the vapor_debug runs,
#  which used the same data, the same 100,000/20,000/80 settings, the same
#  seeds and the same jitter.  Anything above TOL means this pipeline is not
#  reproducing that run and the difference has to be explained before the
#  numbers are used.
vapor_mspe <- vapply(LOSS[["VAPOR (rank 3)"]], mean, numeric(1))
mspe_check <- data.frame(
   quantity = c("config D reference", "this run, mean over chains", "difference"),
   value = sprintf("%.6f", c(CFG$MSPE_REF, mean(vapor_mspe),
                             mean(vapor_mspe) - CFG$MSPE_REF)))
cat(sprintf("\n=== VAPOR MSPE against the config D reference (tolerance %.0e) ===\n",
            CFG$MSPE_TOL))
print(mspe_check, row.names = FALSE)
if (is.na(CFG$MSPE_REF))
   cat("  no reference set for this panel; the check is skipped\n") else
if (abs(mean(vapor_mspe) - CFG$MSPE_REF) > CFG$MSPE_TOL)
   cat(sprintf("  *** FLAG: |difference| = %.3e exceeds %.0e ***\n",
               abs(mean(vapor_mspe) - CFG$MSPE_REF), CFG$MSPE_TOL)) else
   cat("  matches the config D reference within tolerance\n")

md <- mspe_and_dm(LOSS, ORDER)
cat(sprintf("\n=== (1) MSPE over %d test quarters, %d held-out sets ===\n", T_test, CFG$N_HID))
print(md$per_chain, row.names = FALSE)
cat("\n=== (2) Diebold-Mariano vs VAPOR, HLN corrected (positive favours VAPOR) ===\n")
print(md$dm, row.names = FALSE)

ET <- edge_table(CH$vapor, vars)
ci_edges <- data.frame(
   v1 = sub("--.*", "", ET$table$edge), v2 = sub(".*--", "", ET$table$edge),
   pcor = ET$table$pcor, row.names = NULL)
cat(sprintf("\n=== (3) VAPOR edges, pooled 95%% CI rule: %d edges ===\n", nrow(ET$table)))
print(ET$table, row.names = FALSE)

#  item 4: ONE convention.  Every Bayesian method is the mean of its per-draw
#  partial correlations pooled over all four chains; VAR + Glasso and the
#  sample covariance keep their point estimates.  No entry is a single chain.
PCs <- pcs_one_convention(CH, gl$Omega, Om_sample, p, gl_rr3$Omega)
pc_prov <- do.call(rbind, lapply(names(PCs), function(nm){
   z <- PCs[[nm]]
   data.frame(method = nm,
      estimate = if (is.null(attr(z,"n_draws"))) "point estimate"
                 else "mean of per-draw partial correlations",
      chains = if (is.null(attr(z,"n_chains"))) 0L else attr(z,"n_chains"),
      draws  = if (is.null(attr(z,"n_draws")))  0L else attr(z,"n_draws")) }))
cat("\n=== partial-correlation convention used by every table and figure below ===\n")
print(pc_prov, row.names = FALSE)
thr <- threshold_counts(PCs, CFG$THRESHOLDS)
cat("\n=== (4) edge counts at thresholds on |posterior mean partial correlation| ===\n")
print(thr, row.names = FALSE)

rk <- edge_recovery_ranks(ci_edges, PCs, vars)
rownames(rk) <- ET$table$edge
cat("\n=== (5) rank of each VAPOR edge within each method's own ordering ===\n")
print(rk)

#  ORDER uses the long method labels; PCs uses the short ones
PC_FOR <- setNames(lapply(ORDER, function(nm) PCs[[switch(nm,
   `VAPOR (rank 3)` = "VAPOR", `Kastner (VAR+FSV)` = "Kastner", nm)]]), ORDER)
stopifnot(all(c("VAR + Glasso (rank 3)","VAR + PrecFact (rank 3)") %in% names(PCs)))

purity <- if (!is.null(CFG$TOP_N_PURITY))
   group_purity(PCs, vars, CFG$GROUPS, CFG$TOP_N_PURITY) else NULL
if (!is.null(purity)){
   cat(sprintf("\n=== (6) share of the top %d pairs within one official FRED-QD group ===\n",
               CFG$TOP_N_PURITY))
   print(purity, row.names = FALSE) }

rb <- if (!is.null(CFG$RATE_SERIES)) rate_block(ET$draws, vars, CFG$RATE_SERIES) else NULL
if (!is.null(rb)){
   cat("\n=== (7) pooled partial correlations, rate block ===\n"); print(rb, quote = FALSE) }

conv <- convergence_table(CH$vapor, ci_edges, vars, p)
cat("\n=== (8) convergence: log-likelihood, spectral radius, every reported edge ===\n")
print(conv[order(-conv$Rhat), ], row.names = FALSE)

ppc <- ppc_table(CH$vapor, Y_train, p, CFG$PPC_N_REP, CFG$PPC_SEED)
cat(sprintf("\n=== (9) posterior predictive check, %d replicates ===\n", CFG$PPC_N_REP))
print(ppc, row.names = FALSE)

vt <- do.call(rbind, lapply(CH$vapor, function(f)
   data.frame(seed = attr(f,"seed"), minutes = round(attr(f,"minutes"), 2),
              retained_draws = dim(f$Lambda)[3])))
cat("\n=== (10) VAPOR time per chain ===\n"); print(vt, row.names = FALSE)
cat("\n=== (11) residual-precision floor hits per VAPOR chain ===\n")
print(floor_tab, row.names = FALSE)

# ============================================================================
#  FIGURES
# ============================================================================
f1 <- fig_trace(CH$vapor, ci_edges, vars, p, file.path(FIG_DIR, CFG$FIG_TRACE),
                n_edge = CFG$TRACE_EDGES)
f2 <- fig_ppc(ppc, CH$vapor, Y_train, p, file.path(FIG_DIR, CFG$FIG_PPC),
              CFG$PPC_N_REP, CFG$PPC_SEED)
cat(sprintf("\nwrote %s and %s\n", basename(f1), basename(f2)))

# ============================================================================
#  SAVE.  Raw draws per chain, every table, and the fields the existing figure
#  scripts read (vars, p, Y_train, Y_test, long$pc_vapor, long$ci_edges).
# ============================================================================
TABLES <- list(mspe_per_chain = md$per_chain, dm = md$dm,
               mspe_check = mspe_check,
               edges = ET$table, thresholds = thr, edge_ranks = rk,
               purity = purity, rate_block = rb, convergence = conv,
               ppc = ppc, vapor_time = vt, chain_log = chain_log,
               floor_hits = floor_tab, floor_iters = floor_iter_tab)
out <- list(
   config = CFG, order = ORDER,
   provenance = list(R = R.version.string,
      packages = vapply(c("VAPOR","PrecFact","bayesianVARs","glasso","posterior"),
                        function(pk) as.character(packageVersion(pk)), character(1)),
      run_at = format(Sys.time()),
      minutes_total = as.numeric(difftime(Sys.time(), t_start, units="mins"))),
   vars = vars, p = p, T_train = T_train, T_test = T_test,
   Y_train = Y_train, Y_test = Y_test, hid_sets = hid_sets,
   first_stage = list(B_hat = vr$B_hat, B_rank = rt$B, sv = rt$sv,
                      frob_kept = rt$frob_kept, resid_rr = resid_rr),
   glasso = gl, omegas = OM, losses = LOSS, tables = TABLES,
   #  the heatmap and network scripts read methods[[name]]$pcor[[1]], the shape
   #  the superseded p*_final_COMPLETE.R saved.  Keep that shape so those
   #  figures need no change; a method with no completed chain gets an empty
   #  list, which the figure scripts test for.
   #  item 4: $pcor is now ONE matrix per method under the single convention,
   #  kept in the list-of-one shape the figure scripts index with [[1]], so no
   #  panel is a single chain any more.  $pcor_per_chain keeps the old
   #  per-chain point estimates for reference; nothing reads it.
   methods = setNames(lapply(ORDER, function(nm) list(
      label = nm, omega = OM[[nm]],
      pcor = { q <- PC_FOR[[nm]]; if (is.null(q)) list() else list(q) },
      pcor_per_chain = if (length(OM[[nm]])) lapply(OM[[nm]], pcor_from_omega) else list(),
      loss = LOSS[[nm]],
      mspe = if (length(LOSS[[nm]])) vapply(LOSS[[nm]], mean, numeric(1)) else numeric(0),
      n_fits = length(LOSS[[nm]]))), ORDER),
   #  raw draws, per chain, for every stochastic method
   draws = list(vapor = CH$vapor, precfact = CH$pf_raw,
                var_precfact = CH$pf_rr, kastner = CH$kastner),
   long = list(source = "this run", rerun = TRUE,
               pc_vapor = PCs$VAPOR, pc_vapor_convention_b = ET$pc_mean,
               pc_all = PCs, pc_provenance = pc_prov,
               ci_edges = ci_edges, edge_ci = ET$E,
               diag_tab = conv, ppc = ppc, edge_ranks = rk,
               chains_kept = length(CH$vapor)))
saveRDS(out, file.path(CFG$OUT, CFG$OUT_RDS))
#  the heatmap / network / term-structure scripts read FIG_RDS by name
saveRDS(out, file.path(CFG$OUT, CFG$FIG_RDS))
saveRDS(TABLES, file.path(CFG$OUT, sub("_final", "_tables", CFG$OUT_RDS)))
cat(sprintf("\nsaved %s, %s and %s in %.1f min\n", CFG$OUT_RDS, CFG$FIG_RDS,
            sub("_final","_tables",CFG$OUT_RDS),
            as.numeric(difftime(Sys.time(), t_start, units="mins"))))
