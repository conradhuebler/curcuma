#!/bin/bash
# Test: WP-S3 EEQ-Coulomb-Cutoff Auto-Default — heuristic gating
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (May 2026) — Validates:
#   1. Auto greift fuer caffeine (neutral mono-fragment): "EEQ cutoff auto-enabled"
#   2. Auto schaltet NICHT fuer NaCl-Cluster (ionic, max|q|>=0.5 e): "keeping 0.0"
#
# Beide Lauefe sind reine Singlepoints — Heuristic-Trigger wird ueber stdout-Log validiert.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 12: EEQ cutoff auto-detection (WP-S3)"
TEST_DIR="$SCRIPT_DIR"

CONFIG_FILE="auto.json"

setup() {
    cd "$TEST_DIR"
    rm -f "$CONFIG_FILE" caffeine.opt.xyz NaCl.opt.xyz pos_stdout.log pos_stderr.log neg_stdout.log neg_stderr.log
    cleanup_bmt_dirs
    echo '{"gfnff": {"eeq_distance_cutoff_auto": true}}' > "$CONFIG_FILE"
}

run_positive() {
    cd "$TEST_DIR"
    timeout 60 $CURCUMA \
        -sp caffeine.xyz \
        -method gfnff \
        -import_config "$CONFIG_FILE" \
        -verbosity 2 \
        > pos_stdout.log 2> pos_stderr.log
    return $?
}

run_negative() {
    cd "$TEST_DIR"
    timeout 60 $CURCUMA \
        -sp NaCl.xyz \
        -method gfnff \
        -import_config "$CONFIG_FILE" \
        -verbosity 2 \
        > neg_stdout.log 2> neg_stderr.log
    return $?
}

validate_positive() {
    local failed=0

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "EEQ cutoff auto-enabled" pos_stdout.log; then
        echo -e "${GREEN}✓ PASS${NC}: caffeine — auto-enabled log message present"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: caffeine — expected 'EEQ cutoff auto-enabled' missing"
        grep -i "eeq cutoff" pos_stdout.log || echo "  (no EEQ cutoff lines at all)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    fi

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "keeping 0.0" pos_stdout.log; then
        echo -e "${RED}✗ FAIL${NC}: caffeine — unexpected 'keeping 0.0' (auto should have fired)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: caffeine — no 'keeping 0.0' log (auto fired correctly)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return $failed
}

validate_negative() {
    local failed=0

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "EEQ cutoff auto: keeping 0.0" neg_stdout.log; then
        echo -e "${GREEN}✓ PASS${NC}: NaCl — 'keeping 0.0' log present (heuristic rejected ionic system)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: NaCl — expected 'keeping 0.0' missing"
        grep -i "eeq cutoff" neg_stdout.log || echo "  (no EEQ cutoff lines at all)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    fi

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "EEQ cutoff auto-enabled" neg_stdout.log; then
        echo -e "${RED}✗ FAIL${NC}: NaCl — unexpected 'auto-enabled' (heuristic should reject)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: NaCl — no 'auto-enabled' log (correctly rejected)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"

    setup
    run_positive
    local pos_exit=$?
    assert_exit_code $pos_exit 0 "caffeine SP with auto cutoff should succeed"
    validate_positive

    run_negative
    local neg_exit=$?
    assert_exit_code $neg_exit 0 "NaCl SP with auto cutoff should succeed"
    validate_negative

    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
