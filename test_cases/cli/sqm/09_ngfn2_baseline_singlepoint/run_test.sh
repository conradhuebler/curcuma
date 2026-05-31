#!/bin/bash
# Test: ngfn2 (native GFN2) — accuracy vs TBLite gfn2 reference
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — anchor for NATIVE_XTB_ROADMAP AP 1-5
#
# Purpose: ngfn2 routes directly to the native GFN2 implementation
# (currently the OLD monolithic gfn2.cpp, documented "0/7 vs TBLite";
# after AP 2 the new modular curcuma::xtb::XTB(GFN2)).
#
# This test enforces TBLite gfn2 as the truth. The native code must reach
# this accuracy after AP 4 (gradients) and AP 5 (validation).
#
# Status (2026-04-25, OLD GFN2 class): MOST MOLECULES FAIL — caffeine and
# C6H6 are wrong-sign, CH3OCH3 returns 0. This is the documented broken
# state. The test deliberately fails today to anchor the AP roadmap goal.
# When AP 4 is complete, native ngfn2 must match TBLite within tolerance.
#
# Expected current outcome: water/CH4/HCN may still fail at 1 mEh tolerance
# (off by 100-700 mEh today); CH3OCH3/C6H6/caffeine fail catastrophically.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 09: native GFN2 (gfn2) vs TBLite gfn2 reference"
TEST_DIR="$SCRIPT_DIR"
METHOD="gfn2"   # canonical native GFN2 (the former 'ngfn2' alias was removed)

# 1 mEh tolerance — same as test 05 (gfn2 via TBLite). The native must
# eventually reach this accuracy. Currently fails until AP 4 complete.
ENERGY_TOLERANCE="0.001"

# Reference: TBLite gfn2 (truth) — same values as test 05
REF_WATER="-5.07036982186124"
REF_CH4="-4.175218518291275"
REF_HCN="-5.504066165581644"
REF_CH3OH="-8.226117414318244"
REF_CH3OCH3="-11.387712388182692"
REF_C6H6="-15.879607334047481"
REF_CAFFEINE="-42.14723024961651"

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

    assert_numeric_match "$REF" "$E" "$ENERGY_TOLERANCE" "$MOL ngfn2 baseline"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    echo -e "${YELLOW}NOTE${NC}: References are TBLite gfn2 (truth). Native ngfn2 is currently"
    echo -e "${YELLOW}      broken (0/7 vs TBLite, NATIVE_XTB_STATUS) — this test is EXPECTED"
    echo -e "${YELLOW}      to fail today and anchors the goal of AP 4 (analytical gradients)."
    echo

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
