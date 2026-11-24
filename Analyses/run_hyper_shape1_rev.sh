#!/bin/bash
#
# The same as run_hyper_shape1.sh, but working through Analyses/hyper_shape1_runs.txt from
# the end (hummingbirds first), for the second of two computers. All of run_all.sh's options
# work here too (e.g. --dry-run, --jobs 3).

HYPER_REVERSE=1 exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_hyper_shape1.sh" "$@"
