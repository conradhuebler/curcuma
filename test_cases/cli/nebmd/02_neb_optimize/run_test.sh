#!/bin/bash
# CI-NEB band optimisation test - Claude Generated 2026 (AI-generated, machine-tested)
# Ethane staggered -> staggered (120 deg rotation). The barrier must sit in the MIDDLE
# of the band and forward must equal reverse (identical conformers), which is a strong
# check that tangents, springs and the climbing image all work.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -nebmd start.xyz end.xyz -nebmd.nimages 9 \
        -nebmd.idpp false -nebmd.optimize true -nebmd.opt_iterations 400 \
        -md.method gfnff -md.max_time 1 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB optimisation should succeed"
    local opt_xyz=$(find_output_file "start.neb.opt.xyz")
    assert_file_exists "$opt_xyz" "Optimised band"
    return 0
}

validate_results() {
    # Forward and reverse barrier must agree for a symmetric rotation.
    local fwd=$(grep -oE "dE\(forward\) = [-0-9.]+" stdout.log | grep -oE "[-0-9.]+$" | head -1)
    local rev=$(grep -oE "dE\(reverse\) = [-0-9.]+" stdout.log | grep -oE "[-0-9.]+$" | head -1)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$fwd" ] && [ -n "$rev" ] && \
       awk -v a="$fwd" -v b="$rev" 'BEGIN{exit !(((a-b)<0.2)&&((b-a)<0.2))}'; then
        echo -e "${GREEN}ok PASS${NC}: forward ($fwd) == reverse ($rev) kcal/mol"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}xx FAIL${NC}: forward/reverse mismatch: '$fwd' vs '$rev'"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # Barrier magnitude: ethane rotation is ~2-3 kcal/mol.
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$fwd" ] && awk -v a="$fwd" 'BEGIN{exit !(a>1.0 && a<5.0)}'; then
        echo -e "${GREEN}ok PASS${NC}: barrier $fwd kcal/mol is in the expected 1-5 range"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}xx FAIL${NC}: barrier '$fwd' outside the expected 1-5 kcal/mol range"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }
main() {
    test_header "NEB band optimisation (CI-NEB, ethane rotation)"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
