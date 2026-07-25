#!/bin/bash

# Test: ConfSearch with the torsion-recombination phase (Phase 3c)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026)
#
# Smoke test for the ConfGen integration: the phase must run inside the cycle, must not break the
# pipeline, and anything it adds has to survive the normal Phase 4 filters. It does NOT assert a
# conformer yield -- butane has a single rotatable torsion, so there is almost nothing to recombine;
# the yield question is answered on real molecules, not in a 60-second test.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confsearch - 02: torsion-recombination phase (3c) integration"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log
    timeout 280 "$CURCUMA" -confsearch input.xyz -method gfnff \
        -startT 500 -endT 400 -deltaT 100 -time 1200 -threads 1 \
        -confgen_phase true -confgen_max_proposals 6 \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "ConfSearch with Phase 3c completes"

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "signal (6|11)|SIGSEGV|SIGABRT|terminate called|Segmentation fault" stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: crash signature in output"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: no crash signature"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # the phase must actually run (or explicitly report why it was skipped)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "Phase 3c" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: Phase 3c ran"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        grep -oE "Phase 3c.*" stdout.log | head -2 | sed 's/^/  /'
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: Phase 3c never appeared although -confgen_phase true was given"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # the search still produces its result
    local result
    result=$(find_output_file "input.cumulative.opt.accepted.xyz")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$result" ] && [ -s "$result" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: cumulative result written ($(count_xyz_structures "$result") conformers)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: no cumulative result -- the phase broke the pipeline"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
