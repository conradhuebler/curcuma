#!/bin/bash

# Test: the per-repetition selection funnel runs on the RANKING surface (-relax_pes opt, default)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Aug 2026)
#
# What this pins, and why it matters:
#
# Until Aug 2026 a dual-method run optimised every snapshot, applied the energy window,
# deduplicated and re-scored the recombination proposals on the EXPLORATION surface. Only the
# once-per-stage REFINE saw opt_method. That is the same defect -phase3b_two_stage was introduced
# against, one level higher: structures were discarded before the ranking method ever computed
# them. Measured on a 107-atom peptide, three gfnff/gfn2 runs stayed at +28.4/+28.6/+38.4 kJ/mol
# against a GOAT reference while an all-gfn2 run reached -19.8, and the budget shows why -- the
# gfnff/gfn2 run spent 2345 QM optimisations, the all-gfn2 run 7570.
#
# Checks:
#  1. the run states which surface the funnel uses,
#  2. the funnel's files carry the RANKING method, not the exploration method,
#  3. REFINE is skipped with the reason given (its input already IS the accurate ensemble),
#  4. no file of the funnel is written on the exploration surface,
#  5. the ensemble is on one energy scale (all energies within a plausible uff range, no gfnff
#     values mixed in -- the two surfaces differ by hundreds of Hartree).
# uff is the accurate method purely so the test stays fast; the path is method-agnostic.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confsearch - 04: selection funnel on the ranking surface (-relax_pes opt)"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log
    # No -relax_pes here on purpose: this test pins the DEFAULT.
    timeout 280 "$CURCUMA" -confsearch input.xyz -md_method gfnff -opt_method uff \
        -startT 500 -endT 400 -deltaT 100 -time 600 -threads 1 \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "ConfSearch with the default funnel surface completes"

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "signal (6|11)|SIGSEGV|SIGABRT|terminate called|Segmentation fault" stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: crash signature in output"
        grep -iE "signal (6|11)|SIGSEGV|SIGABRT" stdout.log stderr.log | head -3
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: no crash signature"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 1. the choice is stated, not implied
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "selection funnel (RELAX + energy window + dedup + re-scoring) runs at uff" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: the funnel surface is reported"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: the funnel surface is not reported"
        grep -iE "selection funnel|relax_pes" stdout.log | head -3
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. + 4. the funnel writes on the ranking surface and nowhere else
    local relax_opt reduce_opt relax_md reduce_md
    relax_opt=$(find_output_file "input.cycle01_T500K_r1.s2_relax.uff.xyz")
    reduce_opt=$(find_output_file "input.cycle01_T500K_r1.s3_reduce.uff.xyz")
    [ -z "$relax_opt" ] && relax_opt=$(find_output_file "input.cycle01_T500K.s2_relax.uff.xyz")
    [ -z "$reduce_opt" ] && reduce_opt=$(find_output_file "input.cycle01_T500K.s3_reduce.uff.xyz")
    relax_md=$(find_output_file "input.cycle01_T500K_r1.s2_relax.gfnff.xyz")
    [ -z "$relax_md" ] && relax_md=$(find_output_file "input.cycle01_T500K.s2_relax.gfnff.xyz")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s "$relax_opt" ] && [ -s "$reduce_opt" ] && [ -z "$relax_md" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: RELAX and REDUCE wrote on the ranking surface only"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: funnel files wrong (relax_opt=${relax_opt:-none} reduce_opt=${reduce_opt:-none} relax_md=${relax_md:-should not exist})"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 3. REFINE is skipped, and the reason is given
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "REFINE skipped -- the funnel already ran at uff" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: REFINE skipped with its reason stated"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: REFINE skip not reported"
        grep -i "REFINE" stdout.log | head -5
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 5. the MD still runs at the exploration method -- the point of a dual-method run
    TESTS_RUN=$((TESTS_RUN + 1))
    local explore
    explore=$(find_output_file "input.cycle01_T500K_r1.s1_explore.gfnff.xyz")
    [ -z "$explore" ] && explore=$(find_output_file "input.cycle01_T500K.s1_explore.gfnff.xyz")
    if [ -s "$explore" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: the dynamics still runs at the exploration method"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: no gfnff exploration file -- the MD surface changed too"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 6. one energy scale: the funnel ensemble must not mix uff and gfnff energies. A gfnff energy
    #    for this molecule is tens of Hartree, a uff energy is a small positive number -- so a
    #    single sign/magnitude test separates them without hard-coding a golden value.
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s "$reduce_opt" ] && awk '/Energy =/ {
            for (i = 1; i <= NF; i++) if ($i == "=") { e = $(i+1) + 0; break }
            if (e < -1.0) { bad = 1 }
            n++
        } END { exit ((bad || n == 0) ? 1 : 0) }' "$reduce_opt"; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: the funnel ensemble is on one energy scale"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: the funnel ensemble mixes energy scales"
        grep -m3 "Energy =" "$reduce_opt" 2>/dev/null
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
