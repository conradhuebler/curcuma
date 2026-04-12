#!/bin/bash
# compare_baseline.sh — Compare current run against stored baseline per method
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (April 2026)
#
# Usage:
#   bash compare_baseline.sh [threads] [method] [note]
#
# Compares against the first baseline.dat entry for the same method.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MD_THREADS="${1:-4}"
MD_METHOD="${2:-gfnff}"
NOTE="${3:-comparison}"
BASELINE_FILE="baseline.dat"

USE_GPU=false
CURCUMA_METHOD="$MD_METHOD"
if [ "$MD_METHOD" = "gfnff-gpu" ]; then
    USE_GPU=true
    CURCUMA_METHOD="gfnff"
fi

if [ ! -f "$BASELINE_FILE" ]; then
    echo "ERROR: baseline.dat not found. Run record_baseline.sh first." >&2
    exit 1
fi

# Find first baseline entry for this method
BASELINE_LINE=$(grep -v "^#" "$BASELINE_FILE" | grep "| $MD_METHOD |" | head -1)
if [ -z "$BASELINE_LINE" ]; then
    echo "ERROR: No baseline entry for method '$MD_METHOD' in baseline.dat." >&2
    echo "Available methods in baseline.dat:"
    grep -v "^#" "$BASELINE_FILE" | awk -F'|' '{gsub(/ /,"",$3); print $3}' | sort -u
    exit 1
fi

BASELINE_S=$(echo "$BASELINE_LINE"   | awk -F'|' '{gsub(/[^0-9.]/, "", $5); print $5}')
BASELINE_MS=$(echo "$BASELINE_LINE"  | awk -F'|' '{gsub(/[^0-9.]/, "", $6); print $6}')
BASELINE_GIT=$(echo "$BASELINE_LINE" | awk -F'|' '{gsub(/^ +| +$/, "", $2); print $2}')
BASELINE_NOTE=$(echo "$BASELINE_LINE"| awk -F'|' '{gsub(/^ +| +$/, "", $8); print $8}')

# Locate curcuma
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../../.." && pwd)}"
CURCUMA="$PROJECT_ROOT/release/curcuma"
if [ ! -x "$CURCUMA" ]; then CURCUMA="$PROJECT_ROOT/build/curcuma"; fi

GIT_HASH=$(git -C "$PROJECT_ROOT" rev-parse --short HEAD 2>/dev/null || echo "unknown")
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

echo "============================================"
echo " GFN-FF Polymer MD Performance Comparison"
echo "============================================"
echo " Method   : $MD_METHOD"
echo " Baseline : ${BASELINE_S}s  (${BASELINE_MS} ms/step)  git=$BASELINE_GIT  \"$BASELINE_NOTE\""
echo " Running new measurement with $MD_THREADS threads..."
echo "--------------------------------------------"

rm -f input.trj.xyz compare_run.log compare_run.err input.opt.xyz input.restart

CMD=( "$CURCUMA"
    -md input.xyz
    -method "$CURCUMA_METHOD"
    -maxtime 1000
    -threads $MD_THREADS
    -md.seed 42
    -md.no_restart
    -md.rattle_12 false
    -md.print_frequency 100
)
if $USE_GPU; then CMD+=( -gpu cuda ); fi

T_START=$(date +%s%3N)
timeout 3600 "${CMD[@]}" > compare_run.log 2> compare_run.err
EXIT_CODE=$?
T_END=$(date +%s%3N)

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: MD run failed. Check compare_run.err." >&2
    exit 1
fi

T_WALL_MS=$(( T_END - T_START ))
T_WALL_S=$(awk "BEGIN{printf \"%.2f\", $T_WALL_MS/1000}")
T_PER_STEP_MS=$(awk "BEGIN{printf \"%.1f\", $T_WALL_MS/1000}")

SPEEDUP=$(awk -v b="$BASELINE_S" -v n="$T_WALL_S" \
    'BEGIN{if(n>0 && b>0) printf "%.2f", b/n; else print "N/A"}')
DELTA_MS=$(awk -v b="$BASELINE_MS" -v n="$T_PER_STEP_MS" \
    'BEGIN{printf "%.1f", b-n}')

echo " New result : ${T_WALL_S}s  (${T_PER_STEP_MS} ms/step)  git=$GIT_HASH"
echo " Speedup    : ${SPEEDUP}×"
echo " Δ ms/step  : ${DELTA_MS} ms  (positive = faster)"
echo "--------------------------------------------"

echo "$TIMESTAMP | $GIT_HASH | $MD_METHOD | threads=$MD_THREADS | ${T_WALL_S}s | ${T_PER_STEP_MS}ms/step | ok | $NOTE  [vs $BASELINE_GIT speedup=${SPEEDUP}x]" >> "$BASELINE_FILE"
echo " Appended to $BASELINE_FILE"
echo "============================================"
