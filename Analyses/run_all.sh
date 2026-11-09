#!/bin/bash
#
# Run the analyses of the paper one after another (listed in Analyses/paper_runs.txt).
#
#   Analyses/run_all.sh                 build the program, then run every analysis
#   Analyses/run_all.sh empirical       only the runs whose name contains "empirical"
#   Analyses/run_all.sh cynmix_hyper    (any text in the run name selects runs)
#   Analyses/run_all.sh --dry-run       show what would be run, and what is done
#   Analyses/run_all.sh --threads 4     use 4 threads (default: the number of physical cores)
#   Analyses/run_all.sh --runs FILE     use another list of runs
#
# The script can be stopped (Ctrl-C) and started again: a run that finished has a
# <name>.done file in Analyses/Output and is skipped. A run that was interrupted, or
# that failed, is started again from the beginning. Each run writes
# Analyses/Output/<name>.tsv, .tre and .log; progress.log there lists each run's
# start, end and time. The Mac is kept awake (caffeinate) while the script runs.

set -u

ANALYSES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$ANALYSES_DIR")"
OUTPUT_DIR="$ANALYSES_DIR/Output"
BUILD_DIR="$REPO_DIR/build-runs"      # a build of its own, so it cannot clash with another build
PROGRAM="$BUILD_DIR/autoparts"
RUN_LIST="$ANALYSES_DIR/paper_runs.txt"

dryRun=0
threads=""
filter=""
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) dryRun=1 ;;
        --threads) threads="$2"; shift ;;
        --runs)    RUN_LIST="$2"; shift ;;
        -h|--help) sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)        echo "Unknown option: $1" >&2; exit 1 ;;
        *)         filter="$1" ;;
    esac
    shift
done

if [ -z "$threads" ]; then
    threads=$(sysctl -n hw.physicalcpu 2>/dev/null || echo 1)
fi

cd "$REPO_DIR" || exit 1

# build the program (a no-op if it is up to date)
mkdir -p "$OUTPUT_DIR"
if [ $dryRun -eq 0 ]; then
    echo "Building autoparts in $BUILD_DIR ..."
    if ! { cmake -S "$REPO_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release && cmake --build "$BUILD_DIR" -j "$threads"; } > "$OUTPUT_DIR/build.log" 2>&1; then
        echo "The build failed; see $OUTPUT_DIR/build.log" >&2
        exit 1
    fi
    # keep the Mac awake until this script ends
    caffeinate -i -w $$ &
fi

progress="$OUTPUT_DIR/progress.log"

formatTime() {
    local s=$1
    printf "%dh%02dm%02ds" $((s / 3600)) $(((s % 3600) / 60)) $((s % 60))
}

numSelected=0; numDone=0; numRun=0; numFailed=0
# read the list on file descriptor 3, so the program cannot consume it
while IFS= read -r line <&3; do
    case "$line" in ''|'#'*) continue ;; esac
    read -r name seed args <<< "$line"
    case "$name" in *"$filter"*) ;; *) continue ;; esac
    numSelected=$((numSelected + 1))
    out="$OUTPUT_DIR/$name"
    if [ -f "$out.done" ]; then
        numDone=$((numDone + 1))
        [ $dryRun -eq 1 ] && echo "done     $name  ($(cat "$out.done"))"
        continue
    fi
    if [ $dryRun -eq 1 ]; then
        echo "to run   $name"
        echo "         autoparts $args -o Analyses/Output/$name -seed $seed -threads $threads"
        continue
    fi

    mkdir -p "$(dirname "$out")"
    echo "$(date '+%Y-%m-%d %H:%M:%S')  start  $name" | tee -a "$progress"
    start=$(date +%s)
    # shellcheck disable=SC2086
    "$PROGRAM" $args -o "$out" -seed "$seed" -threads "$threads" > "$out.log" 2>&1 < /dev/null
    status=$?
    elapsed=$(( $(date +%s) - start ))
    if [ $status -eq 0 ]; then
        echo "finished $(date '+%Y-%m-%d %H:%M:%S') in $(formatTime $elapsed), seed $seed, $threads threads" > "$out.done"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  done   $name  ($(formatTime $elapsed))" | tee -a "$progress"
        numRun=$((numRun + 1))
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S')  FAILED $name  (exit status $status; see $out.log)" | tee -a "$progress"
        numFailed=$((numFailed + 1))
    fi
done 3< "$RUN_LIST"

if [ $dryRun -eq 1 ]; then
    echo "$numSelected runs selected: $numDone done, $((numSelected - numDone)) to run ($threads threads)"
else
    echo "$numSelected runs selected: $numDone were already done, $numRun run now, $numFailed failed"
fi
[ $numFailed -eq 0 ]
