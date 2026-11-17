#!/bin/bash
#
# The same as run_all.sh, but working through the list of runs from the end, so that
# two computers (one running run_all.sh, the other run_all_rev.sh) can meet in the
# middle. All of run_all.sh's options work here too; see run_all.sh --help.
#
# The two computers do not see each other's output, so neither knows when they have
# met: compare their Analyses/Output/progress.log files now and again, and stop one of
# them when the runs left are the ones the other is running. Then copy one computer's
# Analyses/Output into the other's (the .done files mark which runs are finished).

exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_all.sh" --reverse "$@"
