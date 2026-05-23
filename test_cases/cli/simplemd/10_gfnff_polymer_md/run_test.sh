#!/bin/bash
# Test: gfnff MD stability for polymer (1410 atoms, 1000 steps)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (April 2026) - Validates:
#   1. MD completes 1000 steps without NaN/Inf crash
#   2. Trajectory is written with expected frames
#   3. Total energy drift is bounded over 1 ps
#   4. No instability reported by MD engine
#
# Motivation: Performance and stability regression test for large systems.
#             Polymer (1410 atoms) exercises N²-scaling pair loops, thread
#             balancing, and EEQ solver robustness at scale.
# Note: Timeout is generous (3600 s) — 1000-step GFN-FF MD on 1410 atoms
#       may take 5-30 min depending on hardware.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 10: gfnff MD polymer 1000 steps (1410 atoms)"
TEST_DIR="$SCRIPT_DIR"

# MD parameters
MD_MAXTIME=1000           # 1000 fs = 1 ps (1 step = 1 fs default)
MD_SEED=42
MD_THREADS=4
MD_PRINT=100              # print every 100 steps → 10 output lines
TRJ_FILE="input.trj.xyz"

# Expected: at least 8 frames in trajectory
MIN_FRAMES=8

# Total energy drift tolerance over 1 ps with CSVR thermostat
# Larger molecule → tighter relative drift expected
ENERGY_DRIFT_TOL=0.50

run_test() {
    cd "$TEST_DIR"
    rm -f "$TRJ_FILE" input.opt.xyz input.restart stdout.log stderr.log

    echo "Running: $CURCUMA -md input.xyz -method gfnff -maxtime $MD_MAXTIME"
    echo "         -threads $MD_THREADS -md.seed $MD_SEED -md.no_restart"
    echo "         -md.rattle_12 false -md.print_frequency $MD_PRINT"
    timeout 3600 $CURCUMA \
        -md input.xyz \
        -method gfnff \
        -maxtime $MD_MAXTIME \
        -threads $MD_THREADS \
        -md.seed $MD_SEED \
        -md.no_restart \
        -md.rattle_12 false \
        -md.print_frequency $MD_PRINT \
        > stdout.log 2> stderr.log
    local exit_code=$?

    echo "Exit code: $exit_code"
    return $exit_code
}

validate_results() {
    local failed=0

    # 1. No NaN/Inf in trajectory
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "nan|inf" "$TRJ_FILE" 2>/dev/null; then
        echo -e "${RED}✗ FAIL${NC}: NaN or Inf detected in trajectory file"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: No NaN/Inf in trajectory"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 2. No NaN/Inf in MD step output after step 0
    TESTS_RUN=$((TESTS_RUN + 1))
    local nan_in_steps
    nan_in_steps=$(grep -E "^\s+[1-9][0-9]*\." stdout.log | grep -cE "[[:space:]]-?nan[[:space:]]|[[:space:]]inf[[:space:]]" 2>/dev/null) || nan_in_steps=0
    if [ "${nan_in_steps:-0}" -gt 0 ]; then
        echo -e "${RED}✗ FAIL${NC}: NaN or Inf in MD step output (t > 0)"
        grep -E "^\s+[1-9][0-9]*\." stdout.log | grep -E "nan|inf"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: No NaN/Inf in MD step output"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 3. No instability crash
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qi "Simulation got unstable\|NaN/Inf velocity" stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}✗ FAIL${NC}: MD reported instability"
        grep -i "got unstable\|NaN/Inf" stdout.log stderr.log
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: No instability reported"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 4. Trajectory frame count
    local frames
    frames=$(count_xyz_structures "$TRJ_FILE")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "${frames:-0}" -ge "$MIN_FRAMES" ]; then
        echo -e "${GREEN}✓ PASS${NC}: Trajectory has $frames frames (≥$MIN_FRAMES required)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Trajectory has only ${frames:-0} frames (need $MIN_FRAMES)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        failed=1
    fi

    # 5. Energy drift: first vs last total energy in MD step output
    # MD output columns: time  E_pot  E_pot_avg  E_kin  E_kin_avg  E_tot  ...
    local first_etot last_etot
    first_etot=$(grep -E "^\s+[1-9][0-9]*\." stdout.log | head -1 | awk '{print $6}')
    last_etot=$(grep -E "^\s+[1-9][0-9]*\." stdout.log | tail -1 | awk '{print $6}')

    if [ -n "$first_etot" ] && [ -n "$last_etot" ]; then
        drift=$(awk -v a="$first_etot" -v b="$last_etot" 'BEGIN{d=a-b; if(d<0)d=-d; print d}')
        TESTS_RUN=$((TESTS_RUN + 1))
        if awk -v d="$drift" -v t="$ENERGY_DRIFT_TOL" 'BEGIN{exit(d<t?0:1)}'; then
            echo -e "${GREEN}✓ PASS${NC}: Total energy drift $drift Eh < $ENERGY_DRIFT_TOL Eh"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}✗ FAIL${NC}: Total energy drift $drift Eh ≥ $ENERGY_DRIFT_TOL Eh"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            failed=1
        fi
    else
        echo -e "${YELLOW}⚠ WARN${NC}: Could not extract total energy values for drift check"
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"

    run_test
    local run_exit=$?

    assert_exit_code $run_exit 0 "gfnff MD polymer should complete without crash"

    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
