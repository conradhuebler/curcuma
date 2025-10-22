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
    # Extract RMSD value from output
    # Expected format: "RMSD: X.XXXX" or similar
    local rmsd_line=$(grep -i "rmsd" stdout.log | head -1)

    if [ -n "$rmsd_line" ]; then
        echo -e "${GREEN}✓${NC} RMSD output found: $rmsd_line"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))

        # Extract numeric value (simplified regex)
        local rmsd_value=$(echo "$rmsd_line" | grep -oP '\d+\.\d+' | head -1 || echo "")
        if [ -n "$rmsd_value" ]; then
            echo -e "${BLUE}Info:${NC} RMSD value extracted: $rmsd_value"
        fi
    else
        echo -e "${YELLOW}⚠${NC} RMSD value not found in expected format"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

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
