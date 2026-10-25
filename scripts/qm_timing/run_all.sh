#!/bin/bash
# Runs the case lists one after another (Claude Generated, Oct 2026). Usage: run_all.sh case_file [case_file ...]
# Results: ${QMTIME_WORK:-~/.cache/qmtime}/results.jsonl; table: qm_timing.py table
HERE=$(cd "$(dirname "$0")" && pwd)
for q in "$@"; do
  echo "=== $q $(date) load $(cut -d' ' -f1-3 /proc/loadavg)"
  python3 "$HERE/qm_timing.py" queue "$q"
done
echo "=== done $(date)"
