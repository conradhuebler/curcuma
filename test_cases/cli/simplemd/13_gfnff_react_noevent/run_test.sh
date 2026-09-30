#!/bin/bash

# Test: gfnff react topology mode — no-event regression
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Aug 2026) - Validates:
#   1. topology_mode=react MD on NH3 at 300 K completes without crash
#   2. NO bond events fire at equilibrium (hysteresis stability)
#   3. The trajectory is BIT-IDENTICAL to a plain (auto-mode) run with the
#      same seed — react mode must not perturb non-reactive dynamics.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 13: gfnff react topology no-event regression"
TEST_DIR="$SCRIPT_DIR"

MD_ARGS="-method gfnff -maxtime 500 -md.time_step 0.5 -temperature 300 -threads 1 -md.seed 42 -md.no_restart -md.rattle_12 false"

run_test() {
    cd "$TEST_DIR"
    rm -rf runA runR
    mkdir runA runR
    cp input.xyz runA/ && cp input.xyz runR/

    (cd runA && timeout 200 $CURCUMA -md input.xyz $MD_ARGS -no_bmt > stdout.log 2> stderr.log)
    local exit_a=$?
    (cd runR && timeout 200 $CURCUMA -md input.xyz $MD_ARGS -gfnff.topology_mode react -no_bmt > stdout.log 2> stderr.log)
    local exit_r=$?
    [ $exit_a -eq 0 ] && [ $exit_r -eq 0 ]
    return $?
}

validate_results() {
    local failed=0

    # 1. No bond events at 300 K equilibrium
    TESTS_RUN=$((TESTS_RUN + 1))
    local events
    events=$(grep -ac "REACT bond" runR/stdout.log 2>/dev/null) || events=0
    if [ "${events:-0}" -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: no bond events at 300 K"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $events unexpected bond events"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. Trajectories bit-identical
    TESTS_RUN=$((TESTS_RUN + 1))
    local trjA trjR
    trjA=$(find runA -name "input.trj.xyz" | head -1)
    trjR=$(find runR -name "input.trj.xyz" | head -1)
    if [ -n "$trjA" ] && [ -n "$trjR" ] && diff -q "$trjA" "$trjR" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASS${NC}: react trajectory bit-identical to auto"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: trajectories differ (or missing)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    assert_exit_code $? 0 "both MD runs should complete"
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
