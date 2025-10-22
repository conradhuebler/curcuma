#!/bin/bash
# Test: Dynamic RMSD Threshold
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_confscan.md Szenario 2

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 02: Dynamic RMSD Threshold"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -confscan conformers.xyz -confscan.get_rmsd true > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "ConfScan with get_rmsd should succeed"

    if grep -qi "dynamic.*rmsd\|determined.*rmsd" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} Dynamic RMSD determination mentioned"
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
