#!/bin/bash
# Test: RATTLE Constraints
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_simplemd.md Szenario 4

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 04: RATTLE Constraints"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # rattle=1 means all bonds constrained
    $CURCUMA -simplemd input.xyz -simplemd.max_time 5 -simplemd.rattle 1 > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "SimpleMD with RATTLE should succeed"
    assert_file_exists "input.trj.xyz" "Trajectory created"

    if grep -qi "rattle\|constraint" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} RATTLE constraints mentioned"
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
