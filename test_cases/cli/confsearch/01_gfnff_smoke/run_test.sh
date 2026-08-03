#!/bin/bash
# Test: ConfSearch end-to-end smoke test (gfnff)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026) - regression guard for the use-after-free crash in
# ConfSearch::PerformOptimisation (commit 76e9981): the Phase 2 optimisation used to segfault in the
# CxxThreadPool destructor, which blocked ConfSearch on small inputs. This runs a short multi-cycle
# ConfSearch (strided RMSD-MTD default, ConfSearch safeguards off by default) and requires it to
# complete cleanly and produce the cumulative result.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confsearch - 01: gfnff end-to-end smoke (crash-fix regression)"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log
    timeout 280 "$CURCUMA" -confsearch input.xyz -method gfnff \
        -startT 500 -endT 400 -deltaT 100 -time 1200 -threads 1 \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    # 1. completes without a crash exit code (segfault=139, abort=134)
    assert_exit_code $RUN_EXIT 0 "ConfSearch completes (no crash)"

    # 2. no crash signature in the output (the old bug aborted with signal 6/11)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "signal (6|11)|SIGSEGV|SIGABRT|terminate called|Segmentation fault" \
        stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: crash signature in output"
        grep -iE "signal (6|11)|SIGSEGV|SIGABRT|terminate called|Segmentation fault" stdout.log stderr.log | head -3
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: no crash signature in output"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 3. cumulative result written
    local result
    result=$(find_output_file "input.cumulative.opt.accepted.xyz")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$result" ] && [ -s "$result" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: cumulative result written ($result)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: cumulative result missing (${result:-not found})"
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

main "$@"
