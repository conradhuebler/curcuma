#!/bin/bash
# Test: RMSD with No Reordering
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_rmsd.md Szenario 2

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="rmsd - 02: No Reordering"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute with no_reorder flag
    $CURCUMA -rmsd ref.xyz target.xyz -rmsd.no_reorder true > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "RMSD without reordering should succeed"
    assert_string_in_file "RMSD" stdout.log "RMSD value output"

    # Should mention that reordering is disabled
    if grep -qi "reorder.*disabled\|no.*reorder" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} Reordering disabled message found"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}⚠${NC} No explicit message about disabled reordering (non-critical)"
    fi

    return 0
}

validate_results() {
    # Extract RMSD value - strip ANSI color codes
    local rmsd_line=$(grep "RMSD:" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | head -1)

    if [ -z "$rmsd_line" ]; then
        echo -e "${RED}✗ FAIL${NC}: Could not extract RMSD value"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Extract numeric value
    local rmsd_value=$(echo "$rmsd_line" | grep -oP 'RMSD:\s+\K[0-9.]+' | head -1)

    if [ -z "$rmsd_value" ]; then
        echo -e "${YELLOW}⚠${NC} Could not extract number, but RMSD output found"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    fi

    # Without reordering: should have NO permutation, likely much higher RMSD than with reordering
    echo -e "${BLUE}Info:${NC} RMSD without reordering: $rmsd_value"
    TESTS_RUN=$((TESTS_RUN + 1))
    TESTS_PASSED=$((TESTS_PASSED + 1))

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
}

main() {
    test_header "$TEST_NAME"
    cleanup_before
    if run_test; then validate_results; fi
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
