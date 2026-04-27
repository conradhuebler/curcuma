#!/bin/bash
# record_baseline.sh — Measure and store performance baseline for polymer MD
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (April 2026)
#
# Usage:
#   bash record_baseline.sh [threads] [method] [note]
#
# Examples:
#   bash record_baseline.sh 4 gfnff        "pre_optimization_baseline"
#   bash record_baseline.sh 4 xtb-gfnff   "pre_optimization_baseline"
#   bash record_baseline.sh 4 gfnff-gpu   "pre_optimization_baseline"
#   bash record_baseline.sh 1 gfnff        "single_thread"
#
# Output: baseline.dat (appended, one line per run)
# baseline.dat is tracked in git — git hash ties measurement to code state.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- Parameters ---
MD_THREADS="${1:-4}"
MD_METHOD="${2:-gfnff}"
NOTE="${3:-baseline}"
MD_MAXTIME=1000
MD_SEED=42
MD_PRINT=100
TRJ_FILE="input.trj.xyz"
BASELINE_FILE="baseline.dat"

# GPU flag: gfnff-gpu maps to method=gfnff + -gpu flag
USE_GPU=false
CURCUMA_METHOD="$MD_METHOD"
if [ "$MD_METHOD" = "gfnff-gpu" ]; then
    USE_GPU=true
    CURCUMA_METHOD="gfnff"
fi

# Locate curcuma binary
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../../.." && pwd)}"
CURCUMA="$PROJECT_ROOT/release/curcuma"
if [ ! -x "$CURCUMA" ]; then
    CURCUMA="$PROJECT_ROOT/build/curcuma"
fi
if [ ! -x "$CURCUMA" ]; then
    echo "ERROR: curcuma binary not found. Set PROJECT_ROOT or build first." >&2
    exit 1
fi

# --- System info ---
HOSTNAME=$(hostname)
CPU_MODEL=$(grep -m1 "model name" /proc/cpuinfo 2>/dev/null | sed 's/.*: //' || echo "unknown")
GPU_MODEL=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo "none")
GIT_HASH=$(git -C "$PROJECT_ROOT" rev-parse --short HEAD 2>/dev/null || echo "unknown")
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

echo "============================================"
echo " GFN-FF Polymer MD Performance Baseline"
echo "============================================"
echo " System  : $HOSTNAME"
echo " CPU     : $CPU_MODEL"
echo " GPU     : $GPU_MODEL"
echo " Method  : $MD_METHOD"
echo " Threads : $MD_THREADS"
echo " Git     : $GIT_HASH"
echo " Note    : $NOTE"
echo "--------------------------------------------"

# --- Cleanup ---
rm -f "$TRJ_FILE" baseline_run.log baseline_run.err input.opt.xyz input.restart

# --- Build command ---
CMD=( "$CURCUMA"
    -md input.xyz
    -method "$CURCUMA_METHOD"
    -maxtime $MD_MAXTIME
    -threads $MD_THREADS
    -md.seed $MD_SEED
    -md.no_restart
    -md.rattle_12 false
    -md.print_frequency $MD_PRINT
)
if $USE_GPU; then
    CMD+=( -gpu cuda )
fi

echo "Command: ${CMD[*]}"
echo "Starting 1000-step MD on polymer (1410 atoms)..."

# --- Timed run ---
T_START=$(date +%s%3N)
timeout 3600 "${CMD[@]}" > baseline_run.log 2> baseline_run.err
EXIT_CODE=$?
T_END=$(date +%s%3N)

T_WALL_MS=$(( T_END - T_START ))
T_WALL_S=$(awk "BEGIN{printf \"%.2f\", $T_WALL_MS/1000}")
T_PER_STEP_MS=$(awk "BEGIN{printf \"%.1f\", $T_WALL_MS/1000}")

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: MD run failed (exit $EXIT_CODE). Check baseline_run.err." >&2
    cat baseline_run.err >&2
    exit 1
fi

echo " Done."
echo " Wall time  : ${T_WALL_S} s"
echo " ms / step  : ${T_PER_STEP_MS} ms"

# --- Sanity check ---
QUALITY="ok"
if grep -qiE "nan|inf" "$TRJ_FILE" 2>/dev/null; then
    echo "WARNING: NaN/Inf detected in trajectory." >&2
    QUALITY="NaN_detected"
fi
if grep -qi "Simulation got unstable" baseline_run.log baseline_run.err 2>/dev/null; then
    echo "WARNING: MD reported instability." >&2
    QUALITY="unstable"
fi

# --- Append to baseline.dat ---
if [ ! -f "$BASELINE_FILE" ]; then
    echo "# GFN-FF Polymer MD Baseline — 1000 steps, 1410 atoms" > "$BASELINE_FILE"
    echo "# Fields: timestamp | git_hash | method | threads | wall_s | ms_per_step | quality | note" >> "$BASELINE_FILE"
    echo "# ---" >> "$BASELINE_FILE"
fi

echo "$TIMESTAMP | $GIT_HASH | $MD_METHOD | threads=$MD_THREADS | ${T_WALL_S}s | ${T_PER_STEP_MS}ms/step | $QUALITY | $NOTE" >> "$BASELINE_FILE"

echo "--------------------------------------------"
echo " Result appended to $BASELINE_FILE"
echo "============================================"
