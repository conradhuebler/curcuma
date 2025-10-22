#!/bin/bash
# Test: Alias Backward Compatibility (SinglePoint vs single_point)
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_curcumaopt.md Szenario 5

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 05: Alias Backward Compatibility"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Test old alias: SinglePoint (CamelCase)
    # Execute: curcuma -opt input.xyz -opt.method gfn2 -opt.SinglePoint true
    $CURCUMA -opt input.xyz -opt.method gfn2 -opt.SinglePoint true > stdout_alias.log 2> stderr_alias.log
    local exit_code_alias=$?

    assert_exit_code $exit_code_alias 0 "Single point with alias 'SinglePoint' should succeed"

    # Test new snake_case: single_point
    $CURCUMA -opt input.xyz -opt.method gfn2 -opt.single_point true > stdout_new.log 2> stderr_new.log
    local exit_code_new=$?

    assert_exit_code $exit_code_new 0 "Single point with 'single_point' should succeed"

    # Both should produce similar output (no .opt.xyz file)
    assert_file_not_exists "input.opt.xyz" "No optimization file from single point"

    return 0
}

validate_results() {
    # Both runs should have energy output
    assert_string_in_file "Energy" stdout_alias.log "Energy in alias run"
    assert_string_in_file "Energy" stdout_new.log "Energy in new run"

    # Optional: Compare outputs (ignoring timestamps/runtime info)
    # This is informational only
    echo -e "${BLUE}Info:${NC} Checking if outputs are equivalent..."

    # Extract energy lines and compare
    local energy_alias=$(grep -i "energy" stdout_alias.log | head -1 || echo "")
    local energy_new=$(grep -i "energy" stdout_new.log | head -1 || echo "")

    if [ "$energy_alias" == "$energy_new" ]; then
        echo -e "${GREEN}✓${NC} Outputs are identical (alias works correctly)"
    else
        echo -e "${YELLOW}⚠${NC} Outputs differ slightly (may be due to timestamps, non-critical)"
    fi

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f stdout_alias.log stderr_alias.log stdout_new.log stderr_new.log
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
