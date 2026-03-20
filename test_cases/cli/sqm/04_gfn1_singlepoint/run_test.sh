#!/bin/bash
# Test: GFN1-xTB Single Point Energies — native vs TBLite reference
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — SQM test infrastructure for native method validation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 04: GFN1-xTB Single Point Energies"
TEST_DIR="$SCRIPT_DIR"
METHOD="gfn1"

# Initial tolerance: 1 mEh (loose — documents current state vs TBLite reference)
# TARGET: 1e-5 Eh once native implementation matches TBLite
ENERGY_TOLERANCE="0.001"

# Reference values from ~/src/curcuma/release (TBLite gfn1)
REF_WATER="-5.768641096632452"
REF_CH4="-4.274232453209084"
REF_HCN="-5.782628661866971"
REF_CH3OH="-8.961076723325403"
REF_CH3OCH3="-12.157138513370391"
REF_C6H6="-15.894349466234921"
REF_CAFFEINE="-44.50985543300134"

run_single() {
    local MOL=$1
    local REF=$2
    local LOG="out_${MOL%.xyz}.log"

    cd "$TEST_DIR"
    $CURCUMA -sp "$MOL" -method $METHOD -verbosity 1 > "$LOG" 2>&1
    local EXIT_CODE=$?

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
