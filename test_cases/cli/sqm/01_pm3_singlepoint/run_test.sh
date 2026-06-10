#!/bin/bash
# Test: PM3 Single Point Energies — native vs Ulysses reference
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — SQM test infrastructure for native method validation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 01: PM3 Single Point Energies"
TEST_DIR="$SCRIPT_DIR"
METHOD="pm3"

# Initial tolerance: 1 mEh (loose — documents current state vs Ulysses reference)
# TARGET: 1e-5 Eh once native implementation matches Ulysses
ENERGY_TOLERANCE="0.001"

# Reference values from ~/src/curcuma/release (Ulysses pm3)
REF_WATER="-11.930556635970035"
REF_CH4="-6.634522417383802"
REF_HCN="-10.844860409404877"
REF_CH3OH="-17.42607799450115"
REF_CH3OCH3="-22.9156162384024"
REF_C6H6="-29.50181896982359"
REF_CAFFEINE="-85.95363002544002"

run_single() {
    local MOL=$1
    local REF=$2
    local LOG="out_${MOL%.xyz}.log"

    cd "$TEST_DIR"
    $CURCUMA -sp "$MOL" -method $METHOD -verbosity 1 > "$LOG" 2>&1
    local EXIT_CODE=$?

    # Strip ANSI codes and extract energy
    local E
    E=$(sed 's/\x1b\[[0-9;]*m//g' "$LOG" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1)

    if [ -z "$E" ]; then
        echo -e "${RED}✗ FAIL${NC}: $MOL — no energy output (exit: $EXIT_CODE)"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    assert_numeric_match "$REF" "$E" "$ENERGY_TOLERANCE" "$MOL total energy"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    run_single water.xyz    "$REF_WATER"
    run_single CH4.xyz      "$REF_CH4"
    run_single HCN.xyz      "$REF_HCN"
    run_single CH3OH.xyz    "$REF_CH3OH"
    run_single CH3OCH3.xyz  "$REF_CH3OCH3"
    run_single C6H6.xyz     "$REF_C6H6"
    run_single caffeine.xyz "$REF_CAFFEINE"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
