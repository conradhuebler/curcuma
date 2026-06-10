#!/bin/bash
export PROJECT_ROOT='/home/conrad/src/claude_curcuma/curcuma/test_cases/cli/../..'

# Test: ConfScan Variations
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.threads 8 \
        -confscan.restart false \
        > stdout.log 2> stderr.log

    assert_exit_code $? 0 "ConfScan should succeed"

    # BMT-aware: resolve output files from CWD or BMT directory
    local accepted_file=$(find_output_file "conformers.accepted.xyz")
    local rejected_file=$(find_output_file "conformers.rejected.xyz")
    assert_file_exists "$accepted_file" "Accepted conformers"
    assert_file_exists "$rejected_file" "Rejected conformers"
    return 0
}

validate_results() {
    local accepted_file=$(find_output_file "conformers.accepted.xyz")
    local rejected_file=$(find_output_file "conformers.rejected.xyz")
    local accepted=$(count_xyz_structures "$accepted_file")
    local rejected=$(count_xyz_structures "$rejected_file")
    assert_numeric_match 14 "$accepted" "Accepted conformers"
    assert_numeric_match 30 "$rejected" "Rejected conformers"
    return 0
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "ConfScan Test"
    cleanup_before
    if run_test; then validate_results; fi
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi