#!/bin/bash
# run_all_baselines.sh — Record baseline for gfnff, xtb-gfnff, and gfnff-gpu
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (April 2026)
#
# Usage:
#   bash run_all_baselines.sh [threads] [note]
#
# Example:
#   bash run_all_baselines.sh 4 "pre_optimization_baseline"
#
# Skips gfnff-gpu if no CUDA binary / nvidia-smi found.
# Skips xtb-gfnff if xtb is not in PATH.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

THREADS="${1:-4}"
NOTE="${2:-pre_optimization_baseline}"

FAILED=0

run_method() {
    local method="$1"
    echo ""
    echo "##############################################"
    echo "  Method: $method"
    echo "##############################################"
    bash "$SCRIPT_DIR/record_baseline.sh" "$THREADS" "$method" "$NOTE"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "  → FAILED (exit $rc)" >&2
        FAILED=$(( FAILED + 1 ))
    fi
    return $rc
}

# 1. Native GFN-FF (always available)
run_method "gfnff"

# 2. XTB GFN-FF (requires xtb in PATH or USE_XTB build)
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../../.." && pwd)}"
CURCUMA="${PROJECT_ROOT}/release/curcuma"
if [ ! -x "$CURCUMA" ]; then CURCUMA="${PROJECT_ROOT}/build/curcuma"; fi

if command -v xtb &>/dev/null || "$CURCUMA" -help 2>&1 | grep -qi "xtb"; then
    run_method "xtb-gfnff" || true
else
    echo ""
    echo "  Skipping xtb-gfnff — xtb not found in PATH"
fi

# 3. GFN-FF GPU (requires CUDA build)
if nvidia-smi &>/dev/null; then
    run_method "gfnff-gpu" || true
else
    echo ""
    echo "  Skipping gfnff-gpu — no CUDA device found (nvidia-smi failed)"
fi

echo ""
echo "##############################################"
echo "  Done. Results in baseline.dat"
if [ $FAILED -gt 0 ]; then
    echo "  WARNING: $FAILED method(s) failed"
fi
echo "##############################################"
cat baseline.dat
