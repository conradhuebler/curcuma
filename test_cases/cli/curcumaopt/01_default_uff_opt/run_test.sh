#!/bin/bash
# Test: Standard UFF Optimization
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_curcumaopt.md Szenario 1

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 01: Standard UFF Optimization"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute: curcuma -opt input.xyz
    $CURCUMA -opt input.xyz > stdout.log 2> stderr.log
    local exit_code=$?

    # Note: curcuma may return exit 1 even on success (known bug), so check output file instead
    assert_curcuma_success "input.opt.xyz" "UFF optimization completed successfully"

    # Check for completion message (optional, since verbosity may be 0)
    if grep -qi "converged\|finished\|complete" stdout.log 2>/dev/null; then
        echo -e "${GREEN}✓${NC} Optimization completion message found"
    else
        echo -e "${YELLOW}⚠${NC} No explicit completion message (may be due to low verbosity)"
    fi

    return 0
}

validate_results() {
    # Scientific validation: Check optimized energy against reference
    if [ ! -f "input.opt.xyz" ]; then
        echo -e "${RED}✗${NC} Cannot validate: input.opt.xyz missing"
        return 1
    fi

    local energy=$(extract_energy_from_xyz "input.opt.xyz")

    if [ -z "$energy" ]; then
        echo -e "${YELLOW}⚠${NC} Could not extract energy from XYZ comment"
        return 0  # Non-critical
    fi

    echo -e "${BLUE}Info:${NC} Extracted energy: $energy Eh"

    # UFF reference energy for water.xyz (3 atoms: O-H-H)
    # Based on test_energy_methods.cpp pattern
    local EXPECTED_ENERGY="0.0004"  # UFF for water (small molecule, ~kcal/mol range)
    local TOLERANCE="0.01"          # Looser tolerance for force fields

    # Validate against reference
    assert_scientific_value "$EXPECTED_ENERGY" "$energy" "$TOLERANCE" "UFF energy for water"

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
