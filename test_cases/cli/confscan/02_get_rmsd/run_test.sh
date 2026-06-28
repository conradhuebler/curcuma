#!/bin/bash
# Test: confscan - 02: Dynamic RMSD threshold (get_rmsd)
# Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated
#
# Exercises -confscan.get_rmsd true: the RMSD cutoff is derived from the
# ensemble instead of the fixed default, which changes the filtering outcome
# relative to the default scan (01).

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 02: Dynamic RMSD threshold (get_rmsd)"
TEST_DIR="$SCRIPT_DIR"

# Golden values measured against conformers.xyz. See GOLDEN_REFERENCES.md.
EXP_ACCEPTED=1
EXP_REJECTED=43
EXP_SUMMARY="1 0 0 112"   # accepted reorder reuse skipped

run_test() {
    cd "$TEST_DIR"
    set +e
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.get_rmsd true \
        -confscan.threads 1 \
        -confscan.restart false \
        > stdout.log 2> stderr.log
    local exit_code=$?
    set -e

    assert_exit_code $exit_code 0 "ConfScan should exit 0"
    assert_string_in_file "of 44 total" stdout.log "All 44 input structures were read"
    assert_file_not_empty "$(find_output_file 'conformers.accepted.xyz')" "Accepted conformers file non-empty"
    return 0
}

validate_results() {
    local accepted=$(count_xyz_structures "$(find_output_file 'conformers.accepted.xyz')")
    local rejected=$(count_xyz_structures "$(find_output_file 'conformers.rejected.xyz')")
    assert_numeric_match $EXP_ACCEPTED "$accepted" 0 "Accepted conformer count (get_rmsd)"
    assert_numeric_match $EXP_REJECTED "$rejected" 0 "Rejected conformer count (get_rmsd)"
    local summary=$(extract_confscan_summary stdout.log)
    assert_equals "$EXP_SUMMARY" "$summary" "Result fingerprint (accepted reorder reuse skipped)"
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
