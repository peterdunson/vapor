#!/bin/bash
# ============================================================================
#  THE entry point.  Regenerates every application and simulation table,
#  figure and results_summary.md in the paper, in dependency order.
#
#      bash applications/run_all.sh
#
#  Application code lives in applications/ and writes to applications/results/.
#  Simulation code lives in simulation/ and writes to simulation/results/.
#  results_summary.md is the one shared output and stays in
#  applications/results/tables/.
#
#  Inputs that must already exist and are NOT produced here:
#    applications/data/2026-05-qd.csv      the FRED-QD vintage (in repo)
#    applications/fredqd_groups.csv        series -> FRED-QD group
#    applications/varlist_descriptions.csv series -> description, T, G, read
#                                          by final_varlist_p100.R
#  The four results/fits/*.rds below can be refit with REBUILD_INPUTS=1; see
#  stage 0.
#    results/fits/reused_p29.rds           Kastner, VAR+Glasso (both first
#    results/fits/reused_p100.rds          stages), sample covariance, random
#                                          walk -- reused, never refitted
#    results/fits/kastner_p100_20k.rds     the p = 100 Kastner panel
#    results/fits/kastner_p100_20k_draws.rds
#
#  Each stage is nohup'd and waited on, so the log survives a closed terminal
#  and a failure stops the stages that depend on it.
#
#  NOTE: do not edit this script while it is running.  bash reads a script
#  incrementally by byte offset, so an edit shifts the remaining bytes under
#  the running shell and it resumes mid-token.
# ============================================================================
set -u
REPO=/Users/peterdunson/Desktop/vapor
RES=$REPO/applications/results
SIMRES=$REPO/simulation/results
LOGS=$RES/logs
SIMLOGS=$SIMRES/logs
mkdir -p "$LOGS" "$RES/figures" "$RES/tables" "$RES/fits" \
         "$SIMLOGS" "$SIMRES/figures" "$SIMRES/tables" "$SIMRES/fits" "$SIMRES/sim_f_parts"

_run () {          # _run <logdir> <name> <command...>
   local dir=$1 name=$2; shift 2
   echo "=== $name  $(date '+%Y-%m-%d %H:%M') ==="
   nohup "$@" > "$dir/$name.log" 2>&1 &
   wait $!
   local rc=$?
   if [ $rc -ne 0 ]; then echo "    FAILED (exit $rc), see $dir/$name.log"; fi
   return $rc
}
run    () { _run "$LOGS"    "$@"; }     # application stages
simrun () { _run "$SIMLOGS" "$@"; }     # simulation stages
cd "$REPO" || exit 1

#  0. OPTIONAL: rebuild the four input objects under results/fits/ that nothing
#     else here produces.  Skipped unless REBUILD_INPUTS=1, because the shipped
#     copies ARE the ones the paper reports and a rebuild replaces them:
#
#        REBUILD_INPUTS=1 bash applications/run_all.sh
#
#     The scripts are reconstructions, not the originals (never committed, not
#     recoverable); their headers say what they can and cannot reproduce.
if [ "${REBUILD_INPUTS:-0}" = "1" ]; then
   echo "=== REBUILD_INPUTS=1: rebuilding the four input objects ==="
   for f in reused_p29 reused_p100 kastner_p100_20k kastner_p100_20k_draws; do
      [ -f "$RES/fits/$f.rds" ] && cp "$RES/fits/$f.rds" "$RES/fits/$f.rds.prev"
   done
   run rebuild_kastner_p100 Rscript applications/rebuild_kastner_p100_20k.R || exit 1
   run rebuild_reused_p29   Rscript applications/rebuild_reused.R 29         || exit 1
   run rebuild_reused_p100  Rscript applications/rebuild_reused.R 100        || exit 1
else
   echo "=== inputs: using the shipped copies (set REBUILD_INPUTS=1 to refit) ==="
fi

#  1. the two application fits; everything downstream reads their .rds
run p29  Rscript applications/final_p29.R  || exit 1
run p100 Rscript applications/final_p100.R || exit 1

#  2. analyses that read those fits
run sweeps_p29       Rscript applications/final_sweeps.R 29        || exit 1
run sweeps_p100      Rscript applications/final_sweeps.R 100       || exit 1
run chain_length_p29 Rscript applications/final_chain_length.R 29
run chain_length_p100 Rscript applications/final_chain_length.R 100
run ranks_p29        Rscript applications/final_ranks_p29.R
run pairs_p29        Rscript applications/final_pairs.R 29
run pairs_p100       Rscript applications/final_pairs.R 100
run edge_stability   Rscript applications/final_edge_stability.R
run extra_tables     Rscript applications/final_extra_tables.R
run p100_tables      Rscript applications/final_p100_tables.R
run nipa_p100        Rscript applications/final_nipa_p100.R
run stuck_runs       Rscript applications/stuck_runs.R
run recheck          Rscript applications/recheck_ppc_recovery.R

#  3. tables and figures
run varlist_p100     Rscript applications/final_varlist_p100.R
run fig_heatmap_p29  Rscript applications/heatmap_p29.R
run fig_heatmap_p100 Rscript applications/heatmap_p100.R
run fig_network_p100 Rscript applications/network_p100_grouped.R
run fig_term_p100    Rscript applications/term_structure_p100.R
run fig_trace_p100   Rscript applications/final_trace_p100.R

#  3b. the posterior predictive check for the competing models.  Reads the
#      p<d>_results.rds / p<d>_tables.rds that stage 1 wrote and the reused
#      competitor fits.  Kastner needs B draws and a time-varying Sigma_t,
#      which no stored fit holds, so the script loads results/fits/
#      kastner_ppc_p<d>.rds when it is there and otherwise refits four chains
#      with sv_keep = "all" and caches them.  A cold p = 100 refit is ~80 min;
#      with the cache both dimensions take about two minutes.
run ppc_p29          Rscript applications/ppc_comparison.R 29
run ppc_p100         Rscript applications/ppc_comparison.R 100

#  4. the simulations.  sim_a--sim_e run in the driver's own child
#     processes; sim_f's three (p, T) settings are separate processes, combined
#     afterwards, because each is long enough to be worth running apart.
simrun simulations      Rscript simulation/final_simulations.R
simrun sim_f_p50        Rscript simulation/run_sim_f_setting.R 50 35
simrun sim_f_p100       Rscript simulation/run_sim_f_setting.R 100 70
simrun sim_f_p150       Rscript simulation/run_sim_f_setting.R 150 105
simrun sim_f_combine    Rscript simulation/combine_sim_f.R
simrun sim_tables       Rscript simulation/final_sim_tables.R
simrun fig_sim_e_heatmap Rscript simulation/final_sim_e_heatmap.R

#  5. the summaries: they read every .rds above
run summary          Rscript applications/final_summary.R
simrun sim_sections     Rscript simulation/final_sim_sections.R

#  6. last: the checks for numbers stated in the text of the paper.  This one
#     refits four p = 29 chains at the OLD Gamma(0.1, 0.1) base measure, which
#     is the only way to reproduce the claim about the idiosyncratic precisions
#     staying near zero; the shipped default never reaches the bound.  ~5 min.
run paper_checks     Rscript applications/paper_checks.R

echo "=== done  $(date '+%Y-%m-%d %H:%M') ==="
echo "application tables:  $RES/tables  (incl. paper_checks.txt)"
echo "application figures: $RES/figures"
echo "simulation tables:   $SIMRES/tables"
echo "simulation figures:  $SIMRES/figures"
echo "application logs:    $LOGS"
echo "simulation logs:     $SIMLOGS"
