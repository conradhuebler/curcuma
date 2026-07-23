#!/bin/bash
# Test: RMSD-MTD strided vs legacy A/B + provenance + temperature stability
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026) - docs/RMSD_MTD_TEXTBOOK.md section 6 (checks 2 + 5) and 7.
# Validates: (1) both schemes complete without NaN/instability; (2) the strided default writes the
# provenance table (basename.mtd_hills.csv) with an 'initial' deposit; (3) strided <T> stays bounded
# (no runaway) with the ConfSearch safeguards off; (4) the legacy path still deposits (A/B anchor).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 13: RMSD-MTD strided vs legacy A/B + provenance"
TEST_DIR="$SCRIPT_DIR"
T0=300
TBOUND=600   # a true runaway is NaN or thousands of K; single-walker strided peaks ~365-402 K

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout_str.log stdout_leg.log stderr_str.log stderr_leg.log

    # strided (the new default)
    timeout 120 "$CURCUMA" -md input.xyz -method gfnff -rmsd_mtd \
        -maxtime 1500 -time_step 0.5 -temperature $T0 -thermostat csvr \
        -md.rmsd_mtd_deposit_stride 10 -md.rmsd_mtd_screen false -md.seed 42 -threads 1 \
        > stdout_str.log 2> stderr_str.log
    STR_EXIT=$?

    # legacy (bit-identical old scheme, A/B anchor)
    timeout 120 "$CURCUMA" -md input.xyz -method gfnff -rmsd_mtd -md.rmsd_mtd_scheme legacy \
        -maxtime 1500 -time_step 0.5 -temperature $T0 -thermostat csvr \
        -md.rmsd_mtd_screen false -md.seed 42 -threads 1 \
        > stdout_leg.log 2> stderr_leg.log
    LEG_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $STR_EXIT 0 "strided run completes"
    assert_exit_code $LEG_EXIT 0 "legacy run completes"

    # No instability in either
    local log
    for log in stdout_str.log stdout_leg.log; do
        TESTS_RUN=$((TESTS_RUN + 1))
        if grep -qiE "got unstable|NaN/Inf velocity" "$log" 2>/dev/null; then
            echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: instability reported in $log"
            TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
        else
            echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: no instability in $log"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        fi
    done

    # Provenance: strided writes basename.mtd_hills.csv with a header and an 'initial' deposit
    local hills
    hills=$(find_output_file "input.mtd_hills.csv")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$hills" ] && grep -q "trigger" "$hills" 2>/dev/null && grep -q "initial" "$hills" 2>/dev/null; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: provenance written ($hills, initial deposit present)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: provenance mtd_hills.csv missing/malformed (${hills:-not found})"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # Strided mean-temperature bounded (column 9 of the MD step table), no runaway
    local maxT
    maxT=$(sed -r 's/\x1b\[[0-9;]*m//g' stdout_str.log \
        | grep -E "^[[:space:]]+[0-9]+\.[0-9]+[[:space:]]" \
        | awk '$9 ~ /^[0-9]+\.?[0-9]*$/ {print $9}' | sort -n | tail -1)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$maxT" ] && awk -v t="$maxT" -v b="$TBOUND" 'BEGIN{exit(t<b?0:1)}'; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: strided max <T> = $maxT K < $TBOUND K (bounded, safeguards off)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: strided max <T> = ${maxT:-unknown} K not < $TBOUND K"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # Legacy path still deposits structures (A/B anchor is functional)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qi "stored structures" stdout_leg.log 2>/dev/null; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: legacy path deposits structures"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: legacy path produced no deposits"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

main "$@"
