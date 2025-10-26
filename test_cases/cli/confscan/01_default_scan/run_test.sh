#!/bin/bash
# Test: Default ConfScan
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_confscan.md Szenario 1

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 01: Default Scan"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Claude Generated: Optimized parameters for fast execution (~10s vs 30s+ timeout)
    # Match C++ Unit-Test configuration: subspace method, 8 threads, no restart
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.threads 8 \
        -confscan.restart false \
        > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "ConfScan should succeed"
    assert_file_exists "conformers.accepted.xyz" "Accepted conformers file created"
    assert_file_exists "conformers.rejected.xyz" "Rejected conformers file created"

    return 0
}

validate_results() {
    # Scientific Validation: Count and validate conformer populations
    local accepted=$(count_xyz_structures "conformers.accepted.xyz")
    local rejected=$(count_xyz_structures "conformers.rejected.xyz")

    # Golden references from unit test runs (44 total input structures)
    assert_numeric_match 14 "$accepted" "Expected 14 accepted conformers"
    assert_numeric_match 30 "$rejected" "Expected 30 rejected conformers"

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
