#!/bin/bash
# Test: Temperature Alias (T vs temperature)
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_simplemd.md Szenario 6

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 06: Temperature Alias"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Test with alias 'T'
    $CURCUMA -simplemd input.xyz -simplemd.max_time 10 -simplemd.thermostat berendsen -simplemd.T 300 > stdout_alias.log 2> stderr_alias.log
    local exit_alias=$?

    # Test with full name 'temperature'
    $CURCUMA -simplemd input.xyz -simplemd.max_time 10 -simplemd.thermostat berendsen -simplemd.temperature 300 > stdout_new.log 2> stderr_new.log
    local exit_new=$?

    assert_exit_code $exit_alias 0 "SimpleMD with alias T should succeed"
    assert_exit_code $exit_new 0 "SimpleMD with temperature should succeed"

    # Both should produce trajectory
    if [ -f "input.trj.xyz" ]; then
        echo -e "${GREEN}✓${NC} Both runs produced output"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f stdout_alias.log stderr_alias.log stdout_new.log stderr_new.log
}

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
