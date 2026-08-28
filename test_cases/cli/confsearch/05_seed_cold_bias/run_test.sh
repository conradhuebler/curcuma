#!/bin/bash

# Test: temperature-dependent seed count and bias-density-aware seed ranking
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Aug 2026)
#
# Both mechanisms answer the same measured problem: the seed ranking is purely energetic, and a
# structure that once falls below rank seed_rank can never rise again, because the pool only
# deepens. Measured on a 107-atom peptide: of the ten deepest structures per stage, the number NOT
# produced by the lowest-energy seed rises from 1/10 at 600 K to 10/10 at 300 K -- in the coldest
# stage the best seed contributed nothing, because its surroundings are saturated with hills.
#
#  -seed_rank_cold_factor  raises the seed count linearly towards endT (cheap: the cold stages
#                          produce far fewer snapshots anyway),
#  -seed_bias_penalty      demotes candidates that already carry accumulated bias, ranking only.
#
# Checks: both fire and say so, the raised count is the linear interpolation, the penalty actually
# reorders (a higher-energy candidate outranks a lower-energy one), and -- the important one --
# both are OFF by default, so a run without them is unchanged.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confsearch - 05: cold-stage seed count + bias-density ranking"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f on.log off.log stderr.log
    # startT 500 -> endT 300, deltaT 100: three stages, so the interpolation is checkable by hand
    # (factor 3 => 3 seeds at 500 K, 6 at 400 K, 9 at 300 K).
    timeout 280 "$CURCUMA" -confsearch input.xyz -method gfnff \
        -startT 500 -endT 300 -deltaT 100 -time 400 -repeat 1 -seed_rank 3 -threads 1 \
        -seed_rank_cold_factor 3.0 -seed_bias_penalty 20.0 -confgen_phase false \
        > on.log 2> stderr.log
    RUN_EXIT=$?
    cleanup_bmt_dirs
    timeout 280 "$CURCUMA" -confsearch input.xyz -method gfnff \
        -startT 500 -endT 300 -deltaT 100 -time 400 -repeat 1 -seed_rank 3 -threads 1 \
        -confgen_phase false \
        > off.log 2>> stderr.log
    RUN_EXIT_OFF=$?
}

validate_results() {
    local failed=0
    assert_exit_code $RUN_EXIT 0 "run with both mechanisms completes"
    assert_exit_code $RUN_EXIT_OFF 0 "run with the defaults completes"

    # 1. the raised seed count, and that it is the linear interpolation
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "T = 400 K -- seed_rank raised 3 -> 6" on.log \
       && grep -q "T = 300 K -- seed_rank raised 3 -> 9" on.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: seed count follows the temperature (3 -> 6 -> 9)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: cold-stage seed count wrong"
        grep -i "seed_rank raised" on.log | head -3
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. the bias density is computed and reported
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "bias density around" on.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: bias density around the candidates is measured"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: no bias-density report"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 3. the penalty actually reorders: somewhere a seed with a LARGER gap to the cycle best is
    #    ranked above one with a smaller gap. Under a purely energetic ranking that cannot happen.
    TESTS_RUN=$((TESTS_RUN + 1))
    if awk '/seed +[0-9]+:/ {
              for (i = 1; i <= NF; i++) if ($i == "kJ/mol") { d = $(i-1) + 0; break }
              if (seen && d < prev - 1e-9) { found = 1 }
              prev = d; seen = 1
            } END { exit (found ? 0 : 1) }' on.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: the bias penalty reorders the ranking"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: ranking is still purely energetic"
        grep -E "seed +[0-9]+:" on.log | head -6
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 4. and the report still measures against the LOWEST energy, not the ranking winner:
    #    the first seed of every stage must read +0.00 vs this cycle's best.
    TESTS_RUN=$((TESTS_RUN + 1))
    if ! grep -E "seed +[0-9]+:" on.log | grep -q -- "-[0-9]*\.[0-9]* kJ/mol vs this cycle's best"; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: the report still references the lowest energy"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: negative gap to 'this cycle's best' -- wrong reference"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 5. THE important one: both are off by default
    TESTS_RUN=$((TESTS_RUN + 1))
    if ! grep -q "seed_rank raised" off.log && ! grep -q "bias density around" off.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: both mechanisms are off by default"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: a mechanism fired without being asked for"
        grep -E "seed_rank raised|bias density around" off.log | head -3
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
