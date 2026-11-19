#!/bin/bash
# Reaction-MD example inputs (docs/REACTION_MD_EXAMPLES.md) - Claude Generated (Nov 2026)
# The JSON inputs of test_cases/reaction_md are copied here by CMake; scripts/reaction_md_check.py is run from the repository.
#  A: haber-bosch-gas (forced N-N pull, seed 1): 2 free NH3 at the end
#  B: haber-bosch-gas-control (no external potential): N-N stays below 1.5 A, no NH3
#  C: h2-pt4 (300 steps): H-H > 2.0 A, both H bound to Pt (< 1.8 A)
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

check() {  # check <description> <1 = pass>
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$2" = "1" ]; then
        echo -e "${GREEN}PASS${NC}: $1"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $1"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

run_test() {
    cd "$TEST_DIR"
    rm -rf run_a run_b run_c
    mkdir run_a run_b run_c
    cd run_a
    $CURCUMA -md ../haber.xyz -import_config ../reactions-haber-bosch-gas.json -threads 1 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "A: forced MD should succeed"
    cd ../run_b
    $CURCUMA -md ../haber.xyz -import_config ../reactions-haber-bosch-gas-control.json -threads 1 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "B: control MD should succeed"
    cd ../run_c
    $CURCUMA -md ../pt4.xyz -import_config ../reactions-h2-pt4.json -max_time 300 -threads 1 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "C: H2 on Pt4 MD should succeed"
    return 0
}

validate_results() {
    cd "$TEST_DIR"
    local a b c
    a=$(python3 -I "$PROJECT_ROOT/scripts/reaction_md_check.py" haber run_a --dt 1 --dump 10 2>&1) || a="ERR $a"
    b=$(python3 -I "$PROJECT_ROOT/scripts/reaction_md_check.py" haber run_b --dt 1 --dump 10 2>&1) || b="ERR $b"
    c=$(python3 -I "$PROJECT_ROOT/scripts/reaction_md_check.py" h2-metal run_c --dt 1 --dump 5 2>&1) || c="ERR $c"
    check "A: 2 free NH3 in the last frame" "$(echo "$a" | grep -q 'free NH3 in last frame  : 2' && echo 1 || echo 0)"
    check "A: N-N stretched above 2.5 A by the pull" "$(echo "$a" | grep -q 'N-N max > 2.5 A: True' && echo 1 || echo 0)"
    check "B: control, N2 intact (N-N max <= 2.5 A) and no NH3" "$(echo "$b" | grep -q 'NH3 (free) formed : False | N-N max > 2.5 A: False' && echo 1 || echo 0)"
    check "C: H2 dissociated on Pt4" "$(echo "$c" | grep -q 'EVENT dissociated       : True' && echo 1 || echo 0)"
}

cleanup_before() { cd "$TEST_DIR"; rm -rf run_a run_b run_c; }

main() {
    test_header "SimpleMD Reaction-MD Example Inputs Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
