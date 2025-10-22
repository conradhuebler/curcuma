#!/bin/bash
# Test: Short NVE MD Simulation
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_simplemd.md Szenario 1

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 01: Short NVE Run"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Very short run: 10 fs
    $CURCUMA -simplemd input.xyz -simplemd.max_time 10 -simplemd.thermostat none > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "SimpleMD NVE should succeed"
    assert_file_exists "input.trj.xyz" "Trajectory file created"

    if grep -qi "STARTING MD\|MD.*simulation\|molecular dynamics" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} MD simulation started"
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
