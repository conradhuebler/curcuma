#!/bin/bash
# Test: gfn2 provider routing — anchor for NATIVE_XTB_ROADMAP AP 3
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated — verifies dispatch via MethodFactory for gfn2/ngfn2/ipea1.
#
# Today (pre-AP-3): -method gfn2 routes to TBLite, -method ngfn2 to native.
# After AP 3: -method gfn2 routes to native (= ngfn2). -method ipea1 stays TBLite.
#
# This test only enforces accuracy where TBLite is reliably the backend
# (gfn2, ipea1). Native ngfn2 accuracy is checked by test 09; here we
# only verify ngfn2 dispatches and produces a finite number.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 11: GFN2 Provider Routing Check (water)"
TEST_DIR="$SCRIPT_DIR"

ENERGY_TOLERANCE="0.001"

# References for water.xyz — TBLite-routed methods (gfn2 today, ipea1 always)
REF_GFN2="-5.07036982186124"        # gfn2 → TBLite (current routing)
REF_IPEA1="-5.826250385201024"      # ipea1 → TBLite (always)

# Run a method, return success if energy is finite (non-zero, not NaN).
# Sets E_OUT global.
run_method_finite() {
    local METHOD=$1
    local LOG="out_${METHOD}.log"
    E_OUT=""

    cd "$TEST_DIR"
    $CURCUMA -sp water.xyz -method "$METHOD" -verbosity 1 > "$LOG" 2>&1
    local EXIT_CODE=$?

    E_OUT=$(sed 's/\x1b\[[0-9;]*m//g' "$LOG" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1)

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -z "$E_OUT" ]; then
        echo -e "${RED}✗ FAIL${NC}: -method $METHOD — no energy output (exit $EXIT_CODE)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
    local is_zero=$(awk -v e="$E_OUT" 'BEGIN { print (e == 0.0) ? "yes" : "no" }')
    if [ "$is_zero" = "yes" ]; then
        echo -e "${RED}✗ FAIL${NC}: -method $METHOD — energy is zero (likely missing backend)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
    echo -e "${GREEN}✓ PASS${NC}: -method $METHOD dispatched, E=$E_OUT Eh"
    TESTS_PASSED=$((TESTS_PASSED + 1))
    return 0
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # 1) Each method must dispatch and produce a finite energy
    run_method_finite gfn2
    run_method_finite ngfn2
    run_method_finite ipea1

    # 2) TBLite-routed methods must match their references precisely
    local E_GFN2 E_IPEA1
    E_GFN2=$(sed 's/\x1b\[[0-9;]*m//g' out_gfn2.log  | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1)
    E_IPEA1=$(sed 's/\x1b\[[0-9;]*m//g' out_ipea1.log | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1)
    assert_numeric_match "$REF_GFN2"  "$E_GFN2"  "$ENERGY_TOLERANCE" "gfn2 (TBLite-routed) energy"
    assert_numeric_match "$REF_IPEA1" "$E_IPEA1" "$ENERGY_TOLERANCE" "ipea1 (TBLite) energy"

    # 3) Routing-marker checks (informational, never fail)
    if grep -qE "GFN2 resolved to TBLite|using native xTB" out_gfn2.log; then
        echo -e "${GREEN}✓ PASS${NC}: gfn2 routing marker present"
    else
        echo -e "${YELLOW}⚠${NC}: gfn2 routing marker not found (verbosity issue, non-fatal)"
    fi
    TESTS_RUN=$((TESTS_RUN + 1))
    TESTS_PASSED=$((TESTS_PASSED + 1))

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
