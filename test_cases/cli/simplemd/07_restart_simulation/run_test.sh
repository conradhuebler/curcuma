#!/bin/bash
# Test: Restart Simulation
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_simplemd.md Szenario 7

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 07: Restart Simulation"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Run 1: Create restart file
    $CURCUMA -simplemd input.xyz -simplemd.max_time 5 -simplemd.write_restart_frequency 5 > stdout_run1.log 2> stderr_run1.log
    local exit1=$?

    assert_exit_code $exit1 0 "Initial MD run should succeed"
    assert_file_exists "input.restart.json" "Restart file created"

    # Run 2: Restart from file
    $CURCUMA -simplemd input.xyz -simplemd.max_time 10 -simplemd.restart_file input.restart.json > stdout_run2.log 2> stderr_run2.log
    local exit2=$?

    assert_exit_code $exit2 0 "Restart MD run should succeed"

    if grep -qi "restart\|resuming\|continuing" stdout_run2.log stderr_run2.log; then
        echo -e "${GREEN}✓${NC} Restart mentioned in output"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f stdout_run1.log stderr_run1.log stdout_run2.log stderr_run2.log
}

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
