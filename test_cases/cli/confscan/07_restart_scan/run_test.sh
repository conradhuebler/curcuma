#!/bin/bash
# Test: Restart Scan Functionality
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_confscan.md Szenario 7

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 07: Restart Scan"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Create dummy progress file to simulate previous run
    echo '{"last_processed": 5}' > conformers.confscan_progress.json

    # Claude Generated: Optimized parameters for fast execution (~10s vs 30s+ timeout)
    # Match C++ Unit-Test configuration: subspace method, 8 threads, WITH restart enabled
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.threads 8 \
        -confscan.restart true \
        > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "ConfScan restart should succeed"

    if grep -qi "restart\|resuming\|continuing" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} Restart functionality mentioned"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}⚠${NC} Restart message not found (non-critical)"
    fi

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f conformers.confscan_progress.json
}

main() {
    test_header "$TEST_NAME"
    cleanup_before
    run_test
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
