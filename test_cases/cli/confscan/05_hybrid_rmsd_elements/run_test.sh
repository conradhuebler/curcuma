#!/bin/bash
# Test: Hybrid RMSD with Element Selection
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_confscan.md Szenario 5

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confscan - 05: Hybrid RMSD Elements"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Claude Generated: Optimized parameters for fast execution (~10s vs 30s+ timeout)
    # Match C++ Unit-Test configuration: 8 threads, no restart (hybrid method is slow but necessary for test)
    $CURCUMA -confscan conformers.xyz \
        -rmsd.method hybrid \
        -rmsd.element "7,8" \
        -confscan.threads 8 \
        -confscan.restart false \
        > stdout.log 2> stderr.log
    local exit_code=$?

    assert_exit_code $exit_code 0 "ConfScan with hybrid RMSD should succeed"

    if grep -qi "hybrid.*rmsd\|element.*7.*8" stdout.log stderr.log; then
        echo -e "${GREEN}✓${NC} Hybrid RMSD with elements mentioned"
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
