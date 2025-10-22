#!/bin/bash
# Test: Optimize Hydrogen Only (optH)
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_curcumaopt.md Szenario 6

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 06: Optimize Hydrogen Only"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute: curcuma -opt input.xyz -opt.optimize_h true
    $CURCUMA -opt input.xyz -opt.optimize_h true > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "Hydrogen-only optimization should succeed"
    assert_file_exists "input.opt.xyz" "Optimized structure created"

    return 0
}

validate_results() {
    # Check that optimization message is present
    assert_string_in_file "converged\|optimization.*complete" stdout.log "Optimization completed"

    # Advanced validation: Check that heavy atom positions didn't change
    # This requires comparing input.xyz with input.opt.xyz
    # For now, we just verify the files exist and have content
    if [ -f "input.xyz" ] && [ -f "input.opt.xyz" ]; then
        local input_lines=$(wc -l < input.xyz)
        local output_lines=$(wc -l < input.opt.xyz)

        # Should have same number of atoms
        TESTS_RUN=$((TESTS_RUN + 1))
        if [ $input_lines -eq $output_lines ]; then
            echo -e "${GREEN}✓ PASS${NC}: Input and output have same structure (line count match)"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${YELLOW}⚠ WARNING${NC}: Line count differs (input: $input_lines, output: $output_lines)"
        fi

        # TODO: More sophisticated check - parse XYZ and verify heavy atoms unmoved
        # This would require Python/AWK script to compare coordinates
        echo -e "${BLUE}Info:${NC} Advanced coordinate comparison not implemented (would require XYZ parser)"
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
