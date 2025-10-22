#!/bin/bash
# Template for Curcuma CLI tests
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated
#
# USAGE: Copy this template to your test scenario directory and modify:
# 1. Set TEST_NAME
# 2. Implement run_test() function
# 3. Implement validate_results() function

set -e  # Exit on error (but catch errors in tests)

# Import common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

# Test configuration
TEST_NAME="Template Test"
TEST_DIR="$SCRIPT_DIR"

# Main test execution
run_test() {
    # TODO: Implement test execution
    # Example:
    # $CURCUMA -i input.xyz -opt > stdout.log 2> stderr.log
    # local exit_code=$?
    #
    # assert_exit_code $exit_code 0 "Curcuma execution should succeed"
    # assert_file_exists "input.opt.xyz" "Optimized structure created"
    # assert_string_in_file "converged" stdout.log "Optimization converged"

    echo "TODO: Implement run_test() function"
    return 1
}

# Result validation
validate_results() {
    # TODO: Implement result validation
    # Example:
    # local energy=$(extract_energy_from_xyz "input.opt.xyz")
    # assert_numeric_match -123.456 $energy 0.001 "Final energy matches reference"

    echo "TODO: Implement validate_results() function"
    return 1
}

# Cleanup before test
cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
}

# Cleanup after test
cleanup_after() {
    cleanup_test_artifacts
}

# Main execution
main() {
    test_header "$TEST_NAME"

    cleanup_before

    # Run the test
    if run_test; then
        validate_results
    else
        echo -e "${RED}Test execution failed${NC}"
    fi

    # Optional: Keep artifacts for inspection
    # cleanup_after

    print_test_summary

    # Return failure if any test failed
    if [ $TESTS_FAILED -gt 0 ]; then
        exit 1
    fi
    exit 0
}

# Execute main if script is run directly (not sourced)
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    main "$@"
fi
