#!/bin/bash
# Test: Standard RMSD Calculation
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_rmsd.md Szenario 1

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="rmsd - 01: Standard RMSD Calculation"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute: curcuma -rmsd ref.xyz target.xyz
    $CURCUMA -rmsd ref.xyz target.xyz > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "RMSD calculation should succeed"
    assert_string_in_file "RMSD" stdout.log "RMSD value output"

    return 0
}

validate_results() {
    # Extract RMSD value from output - strip ANSI color codes first
    # Output format: "[colors]RMSD: 2.872143"
    local rmsd_line=$(grep "RMSD:" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | head -1)

    if [ -z "$rmsd_line" ]; then
        echo -e "${RED}✗ FAIL${NC}: Could not find RMSD line"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Extract numeric value after "RMSD:"
    local rmsd_value=$(echo "$rmsd_line" | grep -oP 'RMSD:\s+\K[0-9.]+' | head -1)

    if [ -z "$rmsd_value" ]; then
        echo -e "${RED}✗ FAIL${NC}: Could not extract RMSD number from: $rmsd_line"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Validate against known reference for AAA-bGal (90 atoms)
    # Method: free (with reordering) - known golden value
    local expected_rmsd="2.872143"
    local tolerance="0.0001"

    assert_scientific_value "$expected_rmsd" "$rmsd_value" "$tolerance" "RMSD with reordering (AAA-bGal)"

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
