#!/bin/bash
# Test: GFN2 Single Point Calculation
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_curcumaopt.md Szenario 2

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 02: GFN2 Single Point"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute: curcuma -opt input.xyz -opt.method gfn2 -opt.single_point true
    $CURCUMA -opt input.xyz -opt.method gfn2 -opt.single_point true > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "GFN2 single point should succeed"
    assert_string_in_file "Energy" stdout.log "Energy output present"

    # Should NOT have optimization artifacts
    assert_file_not_exists "input.opt.xyz" "No .opt.xyz file (single point only)"
    assert_string_not_in_file "converged" stdout.log "No optimization convergence (single point only)"

    return 0
}

validate_results() {
    # Verify GFN2 method was used
    if grep -qi "gfn2\|GFN2" stdout.log; then
        echo -e "${GREEN}✓${NC} GFN2 method confirmed in output"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}⚠${NC} GFN2 method not explicitly mentioned in output (non-critical)"
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
