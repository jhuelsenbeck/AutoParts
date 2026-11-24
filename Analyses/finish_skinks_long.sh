#!/bin/bash
#
# Finish the long skink analyses (Analyses/skinks_long_runs.txt), then make the
# convergence table (long runs against the earlier 2M-cycle runs).
#
#   Analyses/finish_skinks_long.sh              4 runs at a time (the default for 16-19 cores)
#   Analyses/finish_skinks_long.sh --jobs 3     any run_all.sh option can be given
#
# Runs that finished (a .done file) are skipped. A run that was interrupted, e.g. by a
# restart of the computer, starts again from the beginning: the program does not save
# its state, and its output files are overwritten. The script can be started again as
# often as needed. The runs go from the end of the list, as before (E(k) = 6.0 and 4.0
# first, then 2.0 and 1.2).

ANALYSES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ANALYSES_DIR" || exit 1

"$ANALYSES_DIR/run_all.sh" --reverse --runs "$ANALYSES_DIR/skinks_long_runs.txt" "$@"
status=$?

# the table covers every pair finished so far, even if some runs failed
python3 "$ANALYSES_DIR/convergence_table.py" Output/empirical/skinks Output/empirical/skinks_long \
    "$ANALYSES_DIR/skinks_convergence.tsv" --subset 7 > "$ANALYSES_DIR/skinks_convergence.md"
echo "Convergence table: $ANALYSES_DIR/skinks_convergence.md (and skinks_convergence.tsv)"
exit $status
