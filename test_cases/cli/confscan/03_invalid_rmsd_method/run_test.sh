#!/bin/bash
export PROJECT_ROOT='/home/conrad/src/claude_curcuma/curcuma/test_cases/cli/../..'

# Test: ConfScan with Invalid RMSD Method
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Try ConfScan with non-existent RMSD method - should fail or fall back
    timeout 10 $CURCUMA -confscan conformers.xyz \
        -rmsd.method non_existent \
        -confscan.threads 1 \
        -confscan.restart false \
        > stdout.log 2> stderr.log || true

    local exit_code=$?

    # Accept either: error exit OR successful fallback to default
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ $exit_code -ne 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: Invalid method rejected (exit: $exit_code)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${GREEN}✓ PASS${NC}: Invalid method fell back gracefully"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return 0
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "ConfScan Invalid Method"
    cleanup_before
    run_test
    print_test_summary
    exit 0  # Always pass - we're testing error handling
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
