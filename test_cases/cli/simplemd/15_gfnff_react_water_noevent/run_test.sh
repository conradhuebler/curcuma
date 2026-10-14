#!/bin/bash

# Test: gfnff react topology mode — hydrogen-bonded water cluster stays inert
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026) - Validates:
#   1. A hydrogen-bonded water cluster at 300 K runs in react mode without a
#      single bond event: the O-H formation radius (1.60 A) lies inside the
#      hydrogen-bond contact range, so this pins that the slack-radius rule
#      keeps hydrogen bonds from being promoted to covalent bonds.
#   2. The trajectory is BIT-IDENTICAL to a plain (auto-mode) run with the same
#      seed — react mode must not perturb non-reactive dynamics.
#   3. RATTLE is refused together with react mode (constraints are frozen at
#      initialisation and cannot follow a changing topology): the guard reports
#      the refusal and not a single frame is integrated. Checked on the output
#      rather than the exit code because curcuma's -md path returns 0 even when
#      Initialise() fails (pre-existing, same reason test_utils.sh offers
#      assert_curcuma_success as a file-existence check). SimpleMD creates the
#      trajectory file during setup, so the assertion is that it stays EMPTY.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 15: gfnff react water cluster no-event"
TEST_DIR="$SCRIPT_DIR"

MD_ARGS="-method gfnff -maxtime 1000 -md.time_step 0.5 -temperature 300 -threads 1 -md.seed 42 -md.no_restart -md.rattle 0"

run_test() {
    cd "$TEST_DIR"
    rm -rf runA runR runG
    mkdir runA runR runG
    cp input.xyz runA/ && cp input.xyz runR/ && cp input.xyz runG/

    (cd runA && timeout 250 $CURCUMA -md input.xyz $MD_ARGS -no_bmt > stdout.log 2> stderr.log)
    local exit_a=$?
    (cd runR && timeout 250 $CURCUMA -md input.xyz $MD_ARGS -gfnff.topology_mode react -no_bmt > stdout.log 2> stderr.log)
    local exit_r=$?
    # Guard run: must fail fast (exit code != 0), never simulate.
    (cd runG && timeout 60 $CURCUMA -md input.xyz -method gfnff -maxtime 10 -md.rattle 1 -gfnff.topology_mode react -no_bmt > stdout.log 2> stderr.log) || true
    [ $exit_a -eq 0 ] && [ $exit_r -eq 0 ]
    return $?
}

validate_results() {
    local failed=0

    # 1. No bond events at 300 K in the hydrogen-bonded cluster
    TESTS_RUN=$((TESTS_RUN + 1))
    local events
    events=$(grep -ac "REACT bond" runR/stdout.log 2>/dev/null) || events=0
    if [ "${events:-0}" -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: no bond events in the water cluster at 300 K"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: $events unexpected bond events"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. Trajectories bit-identical
    TESTS_RUN=$((TESTS_RUN + 1))
    local trjA trjR
    trjA=$(find runA -name "input.trj.xyz" | head -1)
    trjR=$(find runR -name "input.trj.xyz" | head -1)
    if [ -n "$trjA" ] && [ -n "$trjR" ] && diff -q "$trjA" "$trjR" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASS${NC}: react trajectory bit-identical to auto"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: trajectories differ (or missing)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 3. RATTLE + react refused at initialisation: message present AND nothing simulated
    TESTS_RUN=$((TESTS_RUN + 1))
    local guard_trj guard_bytes
    guard_trj=$(find runG -name "input.trj.xyz" | head -1)
    guard_bytes=0
    [ -n "$guard_trj" ] && guard_bytes=$(wc -c < "$guard_trj")
    if grep -aq "rattle is not available with gfnff topology_mode=react" runG/stdout.log runG/stderr.log \
       && [ "$guard_bytes" -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: rattle + react refused before any step (empty trajectory)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: rattle + react was not refused ($guard_bytes bytes of trajectory)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    assert_exit_code $? 0 "both MD runs should complete"
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
