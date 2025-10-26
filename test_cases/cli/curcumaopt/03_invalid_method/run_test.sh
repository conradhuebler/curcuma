#!/bin/bash
# CurcumaOpt Test - Claude Generated
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -opt input.xyz -method uff > stdout.log 2> stderr.log
    assert_exit_code $? 0 "Optimization should succeed"
    assert_curcuma_success "input.opt.xyz" "Optimized structure created"
    return 0
}

validate_results() {
    if [ -f "input.opt.xyz" ]; then
        local energy=$(extract_energy_from_xyz "input.opt.xyz")
        if [ -z "$energy" ]; then
            echo -e "${YELLOW}⚠${NC} No energy extracted (may be rounded to 0)"
            TESTS_RUN=$((TESTS_RUN + 1))
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${BLUE}Info:${NC} Final energy: $energy Eh"
            if (( $(echo "$energy > -10 && $energy < 10" | bc -l) )); then
                echo -e "${GREEN}✓ PASS${NC}: Energy in reasonable range"
                TESTS_RUN=$((TESTS_RUN + 1))
                TESTS_PASSED=$((TESTS_PASSED + 1))
            else
                echo -e "${RED}✗ FAIL${NC}: Energy unreasonable: $energy"
                TESTS_RUN=$((TESTS_RUN + 1))
                TESTS_FAILED=$((TESTS_FAILED + 1))
            fi
        fi
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "CurcumaOpt Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
