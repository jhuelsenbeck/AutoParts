#!/bin/bash
#
# Rerun all of the hyperprior analyses with the new hyperprior on the concentration parameter
# (Analyses/hyper_shape1_runs.txt: 14 gall wasp, 14 skink and 6 simulation runs), then make
# the convergence tables (new runs against the earlier hyperprior runs).
#
#   Analyses/run_hyper_shape1.sh              4 runs at a time (the default for 16-19 cores)
#   Analyses/run_hyper_shape1.sh --dry-run    show what would be run, and what is done
#   Analyses/run_hyper_shape1.sh cynmix       only the gall wasp runs (any run_all.sh option works)
#
# Runs that finished (a .done file) are skipped, so the script can be started again after an
# interruption; an interrupted run starts again from the beginning.

ANALYSES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ANALYSES_DIR" || exit 1

"$ANALYSES_DIR/run_all.sh" --runs "$ANALYSES_DIR/hyper_shape1_runs.txt" "$@"
status=$?
case " $* " in *" --dry-run "*) exit $status ;; esac

# convergence tables: earlier hyperprior runs (labelled "short") against the new ones ("long")
python3 "$ANALYSES_DIR/convergence_table.py" Output/empirical/cynmix_long Output/empirical/cynmix_hyper_shape1 \
    "$ANALYSES_DIR/cynmix_hyper_shape1_convergence.tsv" > "$ANALYSES_DIR/cynmix_hyper_shape1_convergence.md"
python3 "$ANALYSES_DIR/convergence_table.py" Output/empirical/skinks_long Output/empirical/skinks_hyper_shape1 \
    "$ANALYSES_DIR/skinks_hyper_shape1_convergence.tsv" --subset 7 > "$ANALYSES_DIR/skinks_hyper_shape1_convergence.md"
echo "Convergence tables: $ANALYSES_DIR/{cynmix,skinks}_hyper_shape1_convergence.md"
exit $status
