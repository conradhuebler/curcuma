#!/bin/bash
# CurcumaOpt Test - Claude Generated (Nov 2026)
# A non-converged -opt must say so in its final convergence line and report the real gradient
# norm of the last accepted structure. Before the fix, a run that stopped before its first
# accepted step (line search below min_step at iteration 1) printed
# "gradient_norm=0.000000e+00 Eh/Bohr ... [OK]".
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

# $1 log file, $2 label: the line carries NOT CONVERGED and a gradient norm > 0
check_final_line() {
    local log=$1
    local label=$2
    assert_string_in_file "status=NOT CONVERGED" "$log" "$label: final line says NOT CONVERGED"
    local gnorm=$(grep -a "status=NOT CONVERGED" "$log" | tail -1 | sed -n 's/.*gradient_norm=\([^ ]*\) Eh\/Bohr.*/\1/p')
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$gnorm" ] && awk -v g="$gnorm" 'BEGIN { exit !(g > 0) }'; then
        echo -e "${GREEN}PASS${NC}: $label: gradient_norm=$gnorm > 0"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $label: gradient_norm='$gnorm' (expected > 0)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

run_test() {
    cd "$TEST_DIR"
    # (a) step cap reached
    $CURCUMA -opt input.xyz -method gfnff -max_iterations 2 -verbosity 2 > maxiter.log 2>&1
    assert_exit_code $? 1 "max_iterations 2: non-converged run exits non-zero"
    check_final_line maxiter.log "max_iterations 2"

    # (b) zero step at iteration 1 (no accepted step): line search min_step far too large
    $CURCUMA -opt input.xyz -method gfnff -max_iterations 3 -lbfgs_min_step 100 -verbosity 2 > zerostep.log 2>&1
    assert_exit_code $? 1 "zero step at iteration 1: non-converged run exits non-zero"
    check_final_line zerostep.log "zero step at iteration 1"
    return 0
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f maxiter.log zerostep.log; }

main() {
    test_header "CurcumaOpt non-converged final line"
    cleanup_before
    run_test
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
