#!/bin/bash
# Test: confscan - 07: Restart / resume
# Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated
#
# Runs ConfScan twice in the same directory (-no_bmt so the restart file persists
# in the CWD). Run 1 writes curcuma_restart.json; run 2 reads it and resumes.
# The resumed run reaches the same accepted set, but with the reorder rules reused
# (reorder count 0 instead of 5) - a distinct fingerprint proving the restart path.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 07: Restart / resume"
TEST_DIR="$SCRIPT_DIR"

EXP_ACCEPTED=14
EXP_REJECTED=30
EXP_SUMMARY_RUN1="14 5 1 237"   # fresh scan
EXP_SUMMARY_RUN2="14 0 1 237"   # resumed: reorder rules reused -> 0 new reorders

run_test() {
    cd "$TEST_DIR"
    # Fresh start: drop any restart state so run 1 is deterministic.
    rm -f curcuma_restart.json conformers.confscan_progress.json

    # NOTE: single-threaded. The threaded reorder path can intermittently stall (a
    # pre-existing concurrency issue, see AIChangelog / docs); the restart double-run
    # exposes it. threads=1 is deterministic and yields the same goldens.
    set +e
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.threads 1 \
        -no_bmt \
        > stdout.run1.log 2> stderr.run1.log
    local rc1=$?
    set -e
    assert_exit_code $rc1 0 "ConfScan run 1 should exit 0"
    assert_string_in_file "of 44 total" stdout.run1.log "Run 1 read all 44 structures"
    assert_file_exists "curcuma_restart.json" "Run 1 wrote restart file"
    local s1=$(extract_confscan_summary stdout.run1.log)
    assert_equals "$EXP_SUMMARY_RUN1" "$s1" "Run 1 fingerprint (fresh scan)"

    # Resume from the restart file.
    set +e
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.threads 1 \
        -no_bmt \
        > stdout.run2.log 2> stderr.run2.log
    local rc2=$?
    set -e
    assert_exit_code $rc2 0 "ConfScan run 2 (resume) should exit 0"
    assert_string_in_file "of 44 total" stdout.run2.log "Run 2 read all 44 structures"
    return 0
}

validate_results() {
    local accepted=$(count_xyz_structures "conformers.accepted.xyz")
    local rejected=$(count_xyz_structures "conformers.rejected.xyz")
    assert_numeric_match $EXP_ACCEPTED "$accepted" 0 "Accepted conformer count after resume"
    assert_numeric_match $EXP_REJECTED "$rejected" 0 "Rejected conformer count after resume"
    local s2=$(extract_confscan_summary stdout.run2.log)
    assert_equals "$EXP_SUMMARY_RUN2" "$s2" "Run 2 fingerprint (reorder rules reused)"
    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    clean_dir_keep conformers.xyz run_test.sh
}

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
