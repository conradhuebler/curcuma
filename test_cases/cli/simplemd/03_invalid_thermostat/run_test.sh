#!/bin/bash
# SimpleMD Test - Claude Generated (FIXED)
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -md input.xyz -md.max_time 10 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "MD should succeed"
    assert_file_exists "input.trj.xyz" "Trajectory file"
    return 0
}

validate_results() {
    [ ! -f "input.trj.xyz" ] && return 1
    
    frames=$(count_xyz_structures "input.trj.xyz")
    if [ $frames -ge 17 ] && [ $frames -le 23 ]; then
        echo -e "${GREEN}✓ PASS${NC}: Trajectory has $frames frames (expected ~20)"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Expected ~20 frames, got $frames"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "SimpleMD Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
