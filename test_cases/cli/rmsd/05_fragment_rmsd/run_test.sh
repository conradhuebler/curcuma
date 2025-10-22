#!/bin/bash
# Test: Fragment-based RMSD
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_rmsd.md Szenario 5

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="rmsd - 05: Fragment-based RMSD"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Calculate RMSD only for specific fragments
    $CURCUMA -rmsd ref.xyz target.xyz -rmsd.fragment_reference 0 -rmsd.fragment_target 0 > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "Fragment RMSD should succeed"
    assert_string_in_file "RMSD" stdout.log "RMSD output present"

    if grep -qi "fragment" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} Fragment calculation mentioned"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return 0
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
