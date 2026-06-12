#!/bin/bash
# Test: confscan - 03: Unknown RMSD method falls back gracefully
# Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated
#
# An unrecognized -rmsd.method must not crash or hang: it warns and falls back
# to the recommended default (subspace), producing the same result as the
# default scan (01). Regression guard for the silent-no-output bug too.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 03: Unknown RMSD method fallback"
TEST_DIR="$SCRIPT_DIR"

EXP_ACCEPTED=14
EXP_REJECTED=30
EXP_SUMMARY="14 5 1 237"   # identical to the subspace default it falls back to

run_test() {
    cd "$TEST_DIR"
    set +e
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method non_existent_method \
        -confscan.threads 1 \
        -confscan.restart false \
        > stdout.log 2> stderr.log
    local exit_code=$?
    set -e

    assert_exit_code $exit_code 0 "Unknown method should still exit 0 (graceful fallback)"
    # Must not hang or no-op: the ensemble is read and a valid scan runs.
    assert_string_in_file "of 44 total" stdout.log "All 44 input structures were read"
    assert_file_not_empty "$(find_output_file 'conformers.accepted.xyz')" "Accepted conformers file non-empty"
    return 0
}

validate_results() {
    local accepted=$(count_xyz_structures "$(find_output_file 'conformers.accepted.xyz')")
    local rejected=$(count_xyz_structures "$(find_output_file 'conformers.rejected.xyz')")
    assert_numeric_match $EXP_ACCEPTED "$accepted" 0 "Accepted conformer count (fallback)"
    assert_numeric_match $EXP_REJECTED "$rejected" 0 "Rejected conformer count (fallback)"
    local summary=$(extract_confscan_summary stdout.log)
    assert_equals "$EXP_SUMMARY" "$summary" "Fallback result matches subspace default"
    return 0
}

cleanup_before() { cd "$TEST_DIR"; clean_dir_keep conformers.xyz run_test.sh; }

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
