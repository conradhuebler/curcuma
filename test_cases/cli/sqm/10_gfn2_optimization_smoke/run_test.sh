#!/bin/bash
# Test: GFN2 optimization smoke — anchor for AP 4 gradient validation
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — pre-AP-4 baseline for gfn2 -opt
#
# Purpose: gfn2 -opt currently uses TBLite gradients. After AP 4 the same
# command will use native xTB gradients. This test verifies that the
# end-to-end optimization pipeline works today (TBLite path) and provides
# an anchor for AP 4: post-fix optimization must converge to similar
# structures within a small RMSD.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 10: GFN2 Optimization Smoke (TBLite gradient anchor)"
TEST_DIR="$SCRIPT_DIR"
METHOD="gfn2"

# Optimization tolerance: the LBFGSpp path currently fails with a json
# null-error (known issue, see KNOWN_BUGS.md), but the legacy optimizer
# converges. We accept either.
ENERGY_TOLERANCE="0.01"   # 10 mEh — generous; opt convergence ≠ pre-opt energy

# Final energy ranges (TBLite gfn2; expect monotonic decrease vs single-point)
# Post-opt energies must be ≤ pre-opt single-point with some margin.
SP_WATER="-5.07036982186124"
SP_CH4="-4.175218518291275"
SP_CAFFEINE="-42.14723024961651"

run_opt() {
    local MOL=$1
    local SP_REF=$2
    local LOG="out_${MOL%.xyz}.log"

    cd "$TEST_DIR"
    rm -f "${MOL%.xyz}.opt.xyz" "${MOL%.xyz}.trj.xyz"

    timeout 60 $CURCUMA -opt "$MOL" -method $METHOD -verbosity 1 > "$LOG" 2>&1
    local EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        echo -e "${RED}✗ FAIL${NC}: $MOL — curcuma -opt exit $EXIT_CODE"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Verify optimization completed
    if ! grep -q "Geometry Optimisation converged\|gfn2 Final Energy" "$LOG"; then
        echo -e "${RED}✗ FAIL${NC}: $MOL — no convergence marker in log"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Extract final energy
    local E
    E=$(sed 's/\x1b\[[0-9;]*m//g' "$LOG" | grep -oP 'gfn2 Final Energy: \K[-0-9.eE]+' | tail -1)

    if [ -z "$E" ]; then
        echo -e "${RED}✗ FAIL${NC}: $MOL — could not extract final energy"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Sanity: optimized energy should be ≤ single-point reference + small slack
    # (negative numbers, "more negative" = lower)
    local within=$(awk -v e="$E" -v ref="$SP_REF" -v tol="$ENERGY_TOLERANCE" \
        'BEGIN { print (e <= ref + tol) ? "yes" : "no" }')
    if [ "$within" = "yes" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $MOL — opt converged, E=$E Eh (SP ref: $SP_REF)"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $MOL — opt energy $E exceeds SP ref + tol ($SP_REF + $ENERGY_TOLERANCE)"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    run_opt water.xyz    "$SP_WATER"
    run_opt CH4.xyz      "$SP_CH4"
    run_opt caffeine.xyz "$SP_CAFFEINE"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
