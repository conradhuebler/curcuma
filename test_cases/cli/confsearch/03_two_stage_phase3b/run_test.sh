#!/bin/bash

# Test: two-stage high-level re-optimisation (Phase 3b) + per-cycle output
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026)
#
# Covers three things that broke or were missing before:
#  1. -phase3b_two_stage runs crude opt -> dedup filter -> accurate opt and completes. The second
#     ConfScan instance inside one cycle used to segfault in ripser (function-local static simplex
#     enumerators holding references to a destroyed ripser object; patched at configure time).
#  2. The per-cycle ensemble files and the best-per-cycle trajectory are written for BOTH levels.
#  3. The reference energies are reported per level of theory and no cross-method delta is printed
#     (the two energies live on different potential-energy surfaces).
# uff is the accurate method here purely so the test stays fast; the code path is method-agnostic.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confsearch - 03: two-stage Phase 3b + per-cycle output"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log
    timeout 280 "$CURCUMA" -confsearch input.xyz -md_method gfnff -opt_method uff \
        -startT 500 -endT 400 -deltaT 100 -time 600 -threads 1 \
        -phase3b_two_stage true \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "ConfSearch with two-stage Phase 3b completes"

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "signal (6|11)|SIGSEGV|SIGABRT|terminate called|Segmentation fault" stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: crash signature in output"
        grep -iE "signal (6|11)|SIGSEGV|SIGABRT" stdout.log stderr.log | head -3
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: no crash signature (ripser re-entrancy)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # both stages must appear
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "Phase 3b stage 1/2" stdout.log && grep -q "Phase 3b stage 2/2" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: both optimisation stages ran"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: two-stage Phase 3b did not run"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # per-cycle ensembles on BOTH levels of theory
    local md_ens opt_ens md_best opt_best
    md_ens=$(find_output_file "input.cycle01_T500K.ensemble.gfnff.xyz")
    opt_ens=$(find_output_file "input.cycle01_T500K.ensemble.uff.xyz")
    md_best=$(find_output_file "input.best_per_cycle.gfnff.xyz")
    opt_best=$(find_output_file "input.best_per_cycle.uff.xyz")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s "$md_ens" ] && [ -s "$opt_ens" ] && [ -s "$md_best" ] && [ -s "$opt_best" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: per-cycle ensembles + best-per-cycle written on both levels"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: per-cycle output missing (md=${md_ens:-none} opt=${opt_ens:-none} best=${md_best:-none}/${opt_best:-none})"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # the per-cycle ensemble must be energy-sorted
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s "$md_ens" ] && awk '/Energy =/ {
            for (i = 1; i <= NF; i++) if ($i == "=") { e = $(i+1); break }
            if (n > 0 && e < prev - 1e-9) { bad = 1 }
            prev = e; n++
        } END { exit (bad ? 1 : 0) }' "$md_ens"; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: per-cycle ensemble is energy-sorted"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: per-cycle ensemble is not energy-sorted"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # reference energies: one line per level, and NO cross-PES delta
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "Reference energies (optimised input structure" stdout.log \
        && ! grep -qE "Initial energies:.*delta=" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: reference energies reported per level, no cross-PES delta"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: reference-energy reporting wrong"
        grep -E "Reference energies|Initial energies" stdout.log | head -3
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

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
