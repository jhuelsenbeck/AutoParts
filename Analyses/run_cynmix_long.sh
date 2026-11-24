#!/bin/bash
#
# Run the long gall wasp (cynmix) analyses (Analyses/cynmix_long_runs.txt: 28 runs, each
# a 1M-cycle burn-in then 10M cycles sampled every 500), then make the convergence table
# (long runs against the earlier 2M-cycle runs in Output/empirical/cynmix).
#
#   Analyses/run_cynmix_long.sh              4 runs at a time (the default for 16-19 cores)
#   Analyses/run_cynmix_long.sh --jobs 3     any run_all.sh option can be given
#   Analyses/run_cynmix_long.sh --dry-run    show what would be run, and what is done
#
# Runs that finished (a .done file) are skipped. A run that was interrupted, e.g. by a
# restart of the computer, starts again from the beginning: the program does not save
# its state, and its output files are overwritten. The script can be started again as
# often as needed. The output goes to Analyses/Output/empirical/cynmix_long/; the table
# to Analyses/cynmix_convergence.md and .tsv.

ANALYSES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ANALYSES_DIR" || exit 1

"$ANALYSES_DIR/run_all.sh" --reverse --runs "$ANALYSES_DIR/cynmix_long_runs.txt" "$@"
status=$?
case " $* " in *" --dry-run "*) exit $status ;; esac

# the table covers every pair finished so far, even if some runs failed
python3 "$ANALYSES_DIR/convergence_table.py" Output/empirical/cynmix Output/empirical/cynmix_long \
    "$ANALYSES_DIR/cynmix_convergence.tsv" > "$ANALYSES_DIR/cynmix_convergence.md"
echo "Convergence table: $ANALYSES_DIR/cynmix_convergence.md (and cynmix_convergence.tsv)"
exit $status
