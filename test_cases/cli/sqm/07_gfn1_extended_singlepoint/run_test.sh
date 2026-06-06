#!/bin/bash
# Test: GFN1-xTB Single Point — extended molecule set (gfnff-shared molecules)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — anchor points for NATIVE_XTB_ROADMAP AP 1-5
#
# Purpose: Reference values for molecules also used by GFN-FF tests.
# Currently routed through TBLite; after AP 3 hits native xTB.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 07: GFN1-xTB Extended Single Point (gfnff-shared molecules)"
TEST_DIR="$SCRIPT_DIR"
METHOD="gfn1"

ENERGY_TOLERANCE="0.001"

# Reference: TBLite gfn1 (generated 2026-04-25 from release/curcuma)
REF_HCl="-4.683871705040322"
REF_HH="-0.8747123347646926"
REF_OH="-5.068657264131002"
REF_ACETIC_DIMER="-31.57843891611312"
REF_C6H5COOH="-27.416918767094856"

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

    run_single HCl.xyz               "$REF_HCl"
    run_single HH.xyz                "$REF_HH"
    run_single OH.xyz                "$REF_OH"
    run_single acetic_acid_dimer.xyz "$REF_ACETIC_DIMER"
    run_single C6H5COOH.xyz          "$REF_C6H5COOH"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
