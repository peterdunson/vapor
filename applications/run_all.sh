#!/bin/bash
# ============================================================================
#  THE entry point.  Regenerates every application table, figure and
#  results_summary.md in the paper, in dependency order.
#
#      bash applications/run_all.sh
#
#  Inputs that must already exist and are NOT produced here:
#    ~/Desktop/2026-05-qd.csv              the FRED-QD vintage
#    applications/fredqd_groups.csv        series -> FRED-QD group
#    applications/main-77.tex              Table S8 descriptions, read by
#                                          final_varlist_p100.R
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
LOGS=$RES/logs
mkdir -p "$LOGS" "$RES/figures" "$RES/tables" "$RES/fits"

run () {           # run <name> <command...>
   local name=$1; shift
   echo "=== $name  $(date '+%Y-%m-%d %H:%M') ==="
   nohup "$@" > "$LOGS/$name.log" 2>&1 &
   wait $!
   local rc=$?
   if [ $rc -ne 0 ]; then echo "    FAILED (exit $rc), see $LOGS/$name.log"; fi
   return $rc
}
cd "$REPO" || exit 1

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

#  4. the simulations, and the summary last: it reads every .rds above
run simulations      Rscript applications/final_simulations.R
run summary          Rscript applications/final_summary.R

echo "=== done  $(date '+%Y-%m-%d %H:%M') ==="
echo "tables:  $RES/tables"
echo "figures: $RES/figures"
echo "logs:    $LOGS"
