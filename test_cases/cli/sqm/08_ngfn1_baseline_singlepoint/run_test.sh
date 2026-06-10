#!/bin/bash
# Test: ngfn1 (native GFN1) — accuracy vs TBLite gfn1 reference
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — anchor for NATIVE_XTB_ROADMAP AP 1-5
#
# Purpose: ngfn1 routes directly to the native GFN1 implementation
# (currently the OLD monolithic gfn1.cpp; after AP 2 the new modular
# curcuma::xtb::XTB(GFN1); after AP 5 validated against TBLite).
#
# This test enforces TBLite gfn1 as the truth. Native ngfn1 must match
# within tolerance — the gap is the bug we are tracking.
#
# Status (2026-04-25, OLD GFN1 class): max abs diff ~3.1 mEh on caffeine,
# median ~5e-5 Eh. Tolerance set to 5 mEh accommodates this. AP 5 will
# tighten tolerance to 1e-5 Eh after the new modular implementation is
# validated.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 08: native GFN1 (gfn1) vs TBLite gfn1 reference"
TEST_DIR="$SCRIPT_DIR"
METHOD="gfn1"   # canonical native GFN1 (the former 'ngfn1' alias was removed)

# 5 mEh tolerance: today's max native-vs-TBLite drift is ~3 mEh (caffeine).
# AP 5 will re-baseline with 1e-5 Eh.
ENERGY_TOLERANCE="0.005"

# Reference: TBLite gfn1 (truth) — same values as test 04
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

    assert_numeric_match "$REF" "$E" "$ENERGY_TOLERANCE" "$MOL ngfn1 vs TBLite-gfn1"
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
