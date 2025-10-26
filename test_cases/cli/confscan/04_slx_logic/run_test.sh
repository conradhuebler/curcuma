#!/bin/bash
# Test: sLX Logic (Multiple Thresholds)
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_confscan.md Szenario 4

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 04: sLX Logic"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Claude Generated: Optimized parameters for fast execution (~10s vs 30s+ timeout)
    # Match C++ Unit-Test configuration: subspace method, 8 threads, no restart
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method subspace \
        -confscan.threads 8 \
        -confscan.restart false \
        -confscan.sLX "1.5,2.5" \
        > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "ConfScan with sLX should succeed"

    if grep -qi "sLX\|threshold.*1.5.*2.5\|loose.*threshold" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} sLX thresholds mentioned"
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
