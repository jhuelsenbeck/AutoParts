#!/bin/bash
#
# Run the analyses of the paper (listed in Analyses/paper_runs.txt).
#
#   Analyses/run_all.sh                 build the program, then run every analysis, several at
#                                       a time (default: one per 7 physical cores, e.g. 4 on 28)
#   Analyses/run_all.sh --jobs 2        run 2 analyses at a time (the cores are shared among them)
#   Analyses/run_all.sh empirical       only the runs whose name contains "empirical"
#   Analyses/run_all.sh cynmix_hyper    (any text in the run name selects runs)
#   Analyses/run_all.sh --dry-run       show what would be run, and what is done
#   Analyses/run_all.sh --threads 7     use 7 threads per analysis (default: the number of
#                                       physical cores divided by the number of jobs)
#   Analyses/run_all.sh --runs FILE     use another list of runs
#
# One analysis uses at most about a dozen threads well, so on a machine with many cores
# it is faster overall to run several analyses at once. The runs start in the order of
# the list, which puts the longest (hummingbird) runs first.
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
PROGRESS="$OUTPUT_DIR/progress.log"

formatTime() {
    local s=$1
    printf "%dh%02dm%02ds" $((s / 3600)) $(((s % 3600) / 60)) $((s % 60))
}

# Run one analysis (called by xargs below, so that several can run at once):
#   run_all.sh --run-one <threads> <name> <seed> <arguments...>
if [ "${1:-}" = "--run-one" ]; then
    threads=$2; name=$3; seed=$4
    shift 4
    out="$OUTPUT_DIR/$name"
    mkdir -p "$(dirname "$out")"
    echo "$(date '+%Y-%m-%d %H:%M:%S')  start  $name" | tee -a "$PROGRESS"
    start=$(date +%s)
    "$PROGRAM" "$@" -o "$out" -seed "$seed" -threads "$threads" > "$out.log" 2>&1 < /dev/null
    status=$?
    elapsed=$(( $(date +%s) - start ))
    if [ $status -eq 0 ]; then
        echo "finished $(date '+%Y-%m-%d %H:%M:%S') in $(formatTime $elapsed), seed $seed, $threads threads" > "$out.done"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  done   $name  ($(formatTime $elapsed))" | tee -a "$PROGRESS"
        exit 0
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S')  FAILED $name  (exit status $status; see $out.log)" | tee -a "$PROGRESS"
    exit 1
fi

dryRun=0
jobs=""
threads=""
filter=""
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) dryRun=1 ;;
        --jobs)    jobs="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --runs)    RUN_LIST="$2"; shift ;;
        -h|--help) sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)        echo "Unknown option: $1" >&2; exit 1 ;;
        *)         filter="$1" ;;
    esac
    shift
done
numCores=$(sysctl -n hw.physicalcpu 2>/dev/null || echo 1)
if [ -z "$jobs" ]; then
    jobs=$((numCores / 7))
    [ "$jobs" -lt 1 ] && jobs=1
fi
case "$jobs" in *[!0-9]*|0) echo "--jobs needs a positive whole number" >&2; exit 1 ;; esac
if [ -z "$threads" ]; then
    threads=$((numCores / jobs))
    [ "$threads" -lt 1 ] && threads=1
fi

cd "$REPO_DIR" || exit 1
mkdir -p "$OUTPUT_DIR"

# the runs selected that are not done yet, one per line: <name> <seed> <arguments...>
toRun="$OUTPUT_DIR/.to_run.$$"
trap 'rm -f "$toRun"' EXIT
: > "$toRun"
numSelected=0; numDone=0
while IFS= read -r line; do
    case "$line" in ''|'#'*) continue ;; esac
    read -r name seed args <<< "$line"
    case "$name" in *"$filter"*) ;; *) continue ;; esac
    numSelected=$((numSelected + 1))
    if [ -f "$OUTPUT_DIR/$name.done" ]; then
        numDone=$((numDone + 1))
        [ $dryRun -eq 1 ] && echo "done     $name  ($(cat "$OUTPUT_DIR/$name.done"))"
        continue
    fi
    if [ $dryRun -eq 1 ]; then
        echo "to run   $name"
        echo "         autoparts $args -o Analyses/Output/$name -seed $seed -threads $threads"
    fi
    echo "$line" >> "$toRun"
done < "$RUN_LIST"
numToRun=$((numSelected - numDone))

if [ $dryRun -eq 1 ]; then
    echo "$numSelected runs selected: $numDone done, $numToRun to run ($jobs at a time, $threads threads each)"
    exit 0
fi
if [ $numToRun -eq 0 ]; then
    echo "$numSelected runs selected: all done"
    exit 0
fi

# build the program (a no-op if it is up to date)
echo "Building autoparts in $BUILD_DIR ..."
if ! { cmake -S "$REPO_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release && cmake --build "$BUILD_DIR" -j "$numCores"; } > "$OUTPUT_DIR/build.log" 2>&1; then
    echo "The build failed; see $OUTPUT_DIR/build.log" >&2
    exit 1
fi

# keep the Mac awake until this script ends
caffeinate -i -w $$ &

echo "$numToRun runs to do ($numDone of the $numSelected selected are done), $jobs at a time, $threads threads each"
xargs -P "$jobs" -L 1 "$ANALYSES_DIR/$(basename "$0")" --run-one "$threads" < "$toRun"

numFailed=0
while IFS= read -r line; do
    read -r name rest <<< "$line"
    [ -f "$OUTPUT_DIR/$name.done" ] || numFailed=$((numFailed + 1))
done < "$toRun"
echo "$numSelected runs selected: $numDone were already done, $((numToRun - numFailed)) run now, $numFailed failed or stopped"
[ $numFailed -eq 0 ]
