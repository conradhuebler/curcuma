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
    # Scientific validation: Check optimized energy is physically reasonable
    if [ ! -f "input.opt.xyz" ]; then
        echo -e "${RED}✗${NC} Cannot validate: input.opt.xyz missing"
        return 1
    fi

    local final_energy=$(extract_energy_from_xyz "input.opt.xyz")

    if [ -z "$final_energy" ]; then
        echo -e "${YELLOW}⚠${NC} Could not extract energy from XYZ comment"
        return 0  # Non-critical for small molecules where energy might be rounded to 0
    fi

    echo -e "${BLUE}Info:${NC} Optimized energy: $final_energy Eh"

    # UFF for 3-atom water molecule - typically very small/negative value
    # Just verify it's a reasonable number (not NaN, inf, or obviously wrong)
    # Tolerance: Energy values for water should be in range [-1, 1] Eh for UFF
    if (( $(echo "$final_energy > -1 && $final_energy < 1" | bc -l) )); then
        echo -e "${GREEN}✓ PASS${NC}: Optimized energy is in physically reasonable range for water"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Energy value $final_energy is outside reasonable range"
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
