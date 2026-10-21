# ============================================================================
#  results_summary.md, rebuilt from scratch from THIS rerun's outputs only.
#     Rscript applications/final_summary.R
#  Reads results/fits only; fits nothing.  Every section names the object it
#  came from.  Sections that cannot be rebuilt from the new fits are listed
#  at the end as dropped.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 210)
MD <- file.path(TAB_DIR, "results_summary.md")
rd <- function(f){ p <- file.path(FITS_DIR, f); if (file.exists(p)) readRDS(p) else NULL }
R29 <- rd("p29_results.rds");  R100 <- rd("p100_results.rds")
T29 <- R29$tables;             T100 <- R100$tables
SW  <- list(`29` = rd("sweeps_p29_final.rds"), `100` = rd("sweeps_p100_final.rds"))
CL  <- list(`29` = rd("chain_length_p29.rds"), `100` = rd("chain_length_p100.rds"))
PR  <- list(`29` = rd("pairs_p29.rds"),        `100` = rd("pairs_p100.rds"))
EX  <- rd("extra_tables.rds"); ES <- rd("edge_stability.rds")
PX  <- rd("p100_tables_extra.rds"); NI <- rd("nipa_p100.rds")
RK  <- rd("ranks_p29_final.rds"); RC <- rd("recheck_tables.rds"); SR <- rd("stuck_runs.rds")
DROPPED <- character(0)
sec <- function(title, x, note = NULL, src = NULL){
   if (is.null(x)){ DROPPED <<- c(DROPPED, title); return(invisible(NULL)) }
   cat(sprintf("\n## %s\n\n", title), file = MD, append = TRUE)
   if (!is.null(note)) cat(sprintf("%s\n\n", note), file = MD, append = TRUE)
   if (!is.null(src)) cat(sprintf("*Source: %s*\n\n", src), file = MD, append = TRUE)
   cat("```\n", file = MD, append = TRUE)
   cat(paste(utils::capture.output(
      if (is.data.frame(x)) print(x, row.names = FALSE) else print(x)), collapse = "\n"),
      file = MD, append = TRUE)
   cat("\n```\n", file = MD, append = TRUE) }

gh <- tryCatch(system2("git", c("-C","/Users/peterdunson/Desktop/vapor","rev-parse","--short","HEAD"),
                       stdout = TRUE)[1], error = function(e) NA)
cat(sprintf("# VAPOR application results\n\nRebuilt from scratch %s from this rerun's outputs only.\n\nVAPOR %s, git %s. DP base measure Gamma(as = 1, bs = 1); a_dir = b_dir = 0.1; residual-precision floor 1e-8.\nPanel: NA filter on the selected series only; both dimensions train 1996Q1-2020Q4 (100 rows), test 2021Q1-2025Q4 (20 rows).\nPartial-correlation summaries: mean of per-draw partial correlations pooled over four chains for every Bayesian method; point estimates for VAR + Glasso and the sample covariance. MSPE uses each method's precision point estimate, unchanged.\nMain VAR + Glasso and VAR + PrecFact use the full-rank bayesianVARs posterior-mean first stage; the \"(rank 3)\" rows are the rank-3-truncated versions.\nKastner, VAR + Glasso (both first stages), the sample covariance and the random walk are reused, not refitted.\n\nSimulation sections are absent: the simulations have not been rerun under the new base measure.\n",
            format(Sys.time(), "%Y-%m-%d %H:%M"),
            as.character(utils::packageVersion("VAPOR")), gh), file = MD)

for (d in c("29","100")){
   R <- if (d == "29") R29 else R100; TT <- if (d == "29") T29 else T100
   src <- sprintf("p%s_results.rds", d)
   cat(sprintf("\n\n# p = %s\n", d), file = MD, append = TRUE)
   sec(sprintf("p = %s (2a) MSPE per chain and mean (SD)", d), TT$mspe_per_chain, NULL, src)
   sec(sprintf("p = %s (2a) Diebold-Mariano vs VAPOR", d), TT$dm,
       "HLN corrected, h = 1, per-quarter losses averaged over chains; positive favours VAPOR.", src)
   sec(sprintf("p = %s (2a) held-out sets in which VAPOR is lower", d),
       EX[[sprintf("sets_p%s", d)]], "Of the 100 held-out sets.", "extra_tables.rds")
   sec(sprintf("p = %s (2b) VAPOR edges, pooled 95%% CI rule", d), TT$edges,
       "pcor is the pooled mean of the per-draw partial correlations over 4000 draws.", src)
   sec(sprintf("p = %s (2c) all pairs, ranks 1-10 by |pooled mean|", d),
       PR[[d]]$top10, NULL, sprintf("pairs_p%s.rds", d))
   sec(sprintf("p = %s (2d) edge counts at thresholds", d), TT$thresholds, NULL, src)
   sec(sprintf("p = %s (2e) edge recovery, per-edge ranks", d), TT$edge_ranks, NULL, src)
   sec(sprintf("p = %s (2e) edge recovery summary", d),
       RC[[sprintf("recovery_summary_p%s", d)]], NULL, "recheck_tables.rds")
   sec(sprintf("p = %s (2f) convergence", d), TT$convergence,
       "Rank-normalized split R-hat and bulk/tail ESS over the four chains.", src)
   sec(sprintf("p = %s (2g) spectral radius of B", d),
       EX[[sprintf("rho_p%s", d)]], NULL, "extra_tables.rds")
   sec(sprintf("p = %s (2h) one-step-ahead posterior predictive check", d), TT$ppc, NULL, src)
   sec(sprintf("p = %s (2i) residual-precision floor, per series", d), TT$floor_hits, NULL, src)
   sec(sprintf("p = %s (2i) residual-precision floor, iterations with any floor", d),
       TT$floor_iters, "floor_iters counts iterations, burn-in included.", src)
   sec(sprintf("p = %s (2j) sensitivity sweep", d), SW[[d]]$table,
       "One 100,000-iteration chain per setting; B stored for every fit.",
       sprintf("sweeps_p%s_final.rds", d))
   sec(sprintf("p = %s (2j) edge stability over distinct fits", d),
       ES[[sprintf("stability_p%s", d)]],
       "The reported configuration is counted once, not three times.", "edge_stability.rds")
   sec(sprintf("p = %s (2j) edges selected by each of seeds 101-106", d),
       ES[[sprintf("seeds_p%s", d)]], NULL, "edge_stability.rds")
   sec(sprintf("p = %s (2k) chain length", d), CL[[d]]$table,
       "Configs B, A and C refitted here; 100000/20000/80 reuses the main fits.",
       sprintf("chain_length_p%s.rds", d))
   sec(sprintf("p = %s (2l) minutes per chain", d), TT$chain_log, NULL, src)
   sec(sprintf("p = %s (2m) condition number of the training-sample covariance", d),
       EX[[sprintf("cond_p%s", d)]], NULL, "extra_tables.rds")
   sec(sprintf("p = %s runs of >= 5 identical consecutive draws", d), SR[[d]], NULL, "stuck_runs.rds")
}

cat("\n\n# p = 29 only\n", file = MD, append = TRUE)
sec("p = 29 (3) the six pairs among FEDFUNDS, TB3MS, GS1, GS10", PR$`29`$rates, NULL, "pairs_p29.rds")
sec("p = 29 (3) rank and k table (Table S5)", if (!is.null(RK)) RK$agg else NULL, NULL, "ranks_p29_final.rds")

cat("\n\n# p = 100 only\n", file = MD, append = TRUE)
sec("p = 100 (4) rate block", PR$`100`$rates, NULL, "pairs_p100.rds")
sec("p = 100 (4) selected negative edges", PX$negative_edges,
    "Rank of each within every method's own ordering.", "p100_tables_extra.rds")
sec("p = 100 (4) group purity of the top 100 pairs", PX$purity, NULL, "p100_tables_extra.rds")
sec("p = 100 (4) VAR + Glasso: nonzero off-diagonals and group membership",
    EX$glasso_counts, NULL, "extra_tables.rds")
sec("p = 100 (4) VAR + Glasso: 25 largest pairs VAPOR does not select",
    EX$glasso_top25, NULL, "extra_tables.rds")
sec("p = 100 (4) ranks of the selected rate edges", PX$rate_edges, NULL, "p100_tables_extra.rds")
sec("p = 100 (4) GS5--MORTGAGE30US, every method", PX$gs5_mortgage, NULL, "p100_tables_extra.rds")
sec("p = 100 (4) the p = 29 edges in the pooled p = 100 estimate", EX$p29_in_p100,
    "Cross-dimension table.", "extra_tables.rds")
sec("p = 100 (4) NIPA block, per chain", NI$nipa_per_chain, NULL, "nipa_p100.rds")
sec("p = 100 (4) NIPA block, pooled R-hat and ESS", NI$nipa_pooled, NULL, "nipa_p100.rds")
sec("p = 100 (4) NIPA accounting identity, coefficients", NI$coefficients, NULL, "nipa_p100.rds")
sec("p = 100 (4) NIPA accounting identity, fit", NI$fit, NULL, "nipa_p100.rds")
sec("p = 100 (4) NIPA block, smallest eigenvalues of the 6 x 6 correlation matrix",
    NI$smallest_three, sprintf("Condition number %.2f.", NI$condition), "nipa_p100.rds")
sec("p = 100 (4) sample correlations among the five rate series",
    if (!is.null(NI$rate_summary)) data.frame(min = round(NI$rate_summary[1],3),
       median = round(NI$rate_summary[2],3), max = round(NI$rate_summary[3],3)) else NULL,
    NULL, "nipa_p100.rds")
sec("p = 100 (4) partial-correlation convention, per method", PX$provenance, NULL,
    "p100_tables_extra.rds")

cat("\n\n# Figures\n\n```\n", file = MD, append = TRUE)
figs <- list.files(FIG_DIR, pattern = "[.]png$")
cat(paste(sprintf("%-32s %s", figs,
    format(file.mtime(file.path(FIG_DIR, figs)), "%Y-%m-%d %H:%M")), collapse = "\n"),
    file = MD, append = TRUE)
cat("\n```\n", file = MD, append = TRUE)

if (length(DROPPED)){
   cat("\n\n# Sections dropped\n\nNot rebuildable from this rerun's outputs:\n\n", file = MD, append = TRUE)
   cat(paste(sprintf("- %s", DROPPED), collapse = "\n"), file = MD, append = TRUE)
   cat("\n", file = MD, append = TRUE) }
n <- length(readLines(MD, warn = FALSE))
cat(sprintf("wrote %s: %d lines, %d sections, %d dropped\n", MD, n,
            sum(grepl("^## ", readLines(MD, warn = FALSE))), length(DROPPED)))
if (length(DROPPED)) for (x in DROPPED) cat(sprintf("   dropped: %s\n", x))
