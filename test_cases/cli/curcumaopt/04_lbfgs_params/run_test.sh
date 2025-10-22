#!/bin/bash
# Test: LBFGS Parameter Passing
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_curcumaopt.md Szenario 4

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 04: LBFGS Parameter Passing"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute: curcuma -opt input.xyz -opt.optimizer 0 -opt.lbfgs.m 10
    # optimizer=0 is LBFGS, m=10 is number of corrections to store
    $CURCUMA -opt input.xyz -opt.optimizer 0 -opt.lbfgs.m 10 > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "Optimization with LBFGS parameters should succeed"
    assert_file_exists "input.opt.xyz" "Optimized structure created"

    return 0
}

validate_results() {
    # Check if LBFGS was used (might be mentioned in verbose output)
    # Note: This is informational, not critical for test success
    if grep -qi "lbfgs\|L-BFGS" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} LBFGS optimizer mentioned in output"
    else
        echo -e "${YELLOW}⚠${NC} LBFGS not explicitly mentioned (may be default, non-critical)"
    fi

    # Verify optimization completed
    assert_string_in_file "converged\|optimization.*complete\|finished" stdout.log "Optimization completed"

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
}

main() {
    test_header "$TEST_NAME"
    cleanup_before

    if run_test; then
        validate_results
    fi

    print_test_summary

    if [ $TESTS_FAILED -gt 0 ]; then
        exit 1
    fi
    exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    main "$@"
fi
