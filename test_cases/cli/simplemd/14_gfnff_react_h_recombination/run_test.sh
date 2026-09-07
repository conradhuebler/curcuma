#!/bin/bash

# Test: gfnff react topology mode — H atom recombination
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Aug 2026) - Validates:
#   1. 4 free H atoms in a tight spherical wall at 3000 K recombine:
#      at least one "REACT bond formed" event fires
#   2. every rebuild logs its dE_jump (the accepted energy discontinuity
#      is measured, not hidden)
#   3. the run stays numerically stable (no NaN/Inf, no crash)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 14: gfnff react H recombination"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    rm -f stdout.log stderr.log input.trj.xyz
    cleanup_bmt_dirs
    timeout 280 $CURCUMA -md input.xyz -method gfnff -gfnff.topology_mode react \
        -temperature 3000 -maxtime 5000 -md.time_step 0.5 -threads 1 \
        -md.seed 42 -md.no_restart -md.rattle_12 false \
        -md.wall_type spheric -md.wall_radius 2.5 \
        -no_bmt > stdout.log 2> stderr.log
    return $?
}

validate_results() {
    local failed=0

    # 1. At least one bond formation event
    TESTS_RUN=$((TESTS_RUN + 1))
    local formed
    formed=$(grep -ac "REACT bond formed" stdout.log 2>/dev/null) || formed=0
    if [ "${formed:-0}" -ge 1 ]; then
        echo -e "${GREEN}✓ PASS${NC}: $formed bond formation event(s)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: no bond formation events"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. Every rebuild logged a dE_jump
    TESTS_RUN=$((TESTS_RUN + 1))
    local rebuilds jumps
    rebuilds=$(grep -ac "REACT rebuild" stdout.log 2>/dev/null) || rebuilds=0
    jumps=$(grep -ac "dE_jump" stdout.log 2>/dev/null) || jumps=0
    if [ "${rebuilds:-0}" -ge 1 ] && [ "$rebuilds" -eq "$jumps" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $rebuilds rebuild(s), each with logged dE_jump"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: rebuilds=$rebuilds, dE_jump lines=$jumps"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 3. Numerical stability
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "Simulation got unstable|NaN/Inf velocity" stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}✗ FAIL${NC}: MD reported instability"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: no instability reported"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    assert_exit_code $? 0 "react MD should complete without crash"
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
