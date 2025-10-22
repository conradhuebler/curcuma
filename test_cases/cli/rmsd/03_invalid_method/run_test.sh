#!/bin/bash
# Test: Invalid RMSD Method - Error Handling
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_rmsd.md Szenario 3

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="rmsd - 03: Invalid Method Error Handling"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -rmsd ref.xyz target.xyz -rmsd.method non_existent > stdout.log 2> stderr.log
    local exit_code=$?

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ $exit_code -ne 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: Correctly failed with invalid method (exit: $exit_code)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Should have failed but succeeded"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    assert_string_in_file "error\|invalid\|unknown" stderr.log stdout.log "Error message present"
    return 0
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
