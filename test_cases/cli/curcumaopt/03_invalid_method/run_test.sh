#!/bin/bash
# Test: Invalid Method - Error Handling
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_curcumaopt.md Szenario 3

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 03: Invalid Method Error Handling"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute with invalid method: curcuma -opt input.xyz -opt.method non_existent_method
    # This SHOULD fail
    $CURCUMA -opt input.xyz -opt.method non_existent_method > stdout.log 2> stderr.log
    local exit_code=$?

    # We expect a non-zero exit code
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ $exit_code -ne 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: Program correctly failed with invalid method (exit code: $exit_code)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Program should have failed but succeeded (exit code: 0)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # Check for error message in stderr or stdout
    if grep -qi "error\|invalid\|unknown\|not.*recognized\|not.*found" stderr.log stdout.log; then
        echo -e "${GREEN}✓ PASS${NC}: Error message present in output"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}⚠ WARNING${NC}: No clear error message found (stderr and stdout shown below)"
        echo "=== STDERR ==="
        cat stderr.log
        echo "=== STDOUT ==="
        cat stdout.log
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    return 0
}

validate_results() {
    # No output files should be created
    assert_file_not_exists "input.opt.xyz" "No output file created on error"

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
