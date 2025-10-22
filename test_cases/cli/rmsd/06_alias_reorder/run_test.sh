#!/bin/bash
# Test: Alias Compatibility (reorder vs force_reorder)
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_rmsd.md Szenario 6

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="rmsd - 06: Alias Compatibility"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Test with old alias
    $CURCUMA -rmsd ref.xyz target.xyz -rmsd.reorder true > stdout_alias.log 2> stderr_alias.log
    local exit_alias=$?

    # Test with new name
    $CURCUMA -rmsd ref.xyz target.xyz -rmsd.force_reorder true > stdout_new.log 2> stderr_new.log
    local exit_new=$?

    assert_exit_code $exit_alias 0 "RMSD with alias should succeed"
    assert_exit_code $exit_new 0 "RMSD with new name should succeed"

    # Both should produce output
    assert_string_in_file "RMSD" stdout_alias.log "RMSD in alias output"
    assert_string_in_file "RMSD" stdout_new.log "RMSD in new output"

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
