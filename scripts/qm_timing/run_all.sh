#!/bin/bash
# Runs the case lists one after another (Claude Generated, Oct 2026). Usage: run_all.sh case_file [case_file ...]
# Results: ${QMTIME_WORK:-~/.cache/qmtime}/results.jsonl; table: qm_timing.py table
HERE=$(cd "$(dirname "$0")" && pwd)
for q in "$@"; do
  echo "=== $q $(date) load $(cut -d' ' -f1-3 /proc/loadavg)"
  # machine state before each batch: GPU memory held by other processes (e.g. ollama) distorts the GPU runs
  nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu,temperature.gpu --format=csv,noheader 2>/dev/null
  nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader 2>/dev/null
  pgrep -a ollama || true
  python3 "$HERE/qm_timing.py" queue "$q"
done
echo "=== done $(date)"
