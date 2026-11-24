#!/bin/bash
# Test: -scf_precond frag on the device-resident ROCm GFN2/GFN1 SCF (docs/SCF_PRECOND.md, section GPU; ROCm counterpart of cli_sqm_13)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only.
#
# (a) -gpu rocm without the switch: the resident ROCm loop runs, energy equals the CPU reference
# (b) -gpu rocm -scf_precond frag: M is applied on the device, same energy (1e-7 Eh), fewer SCF iterations
# (c) the GPU frag energy equals the CPU frag energy (1e-7 Eh)
# (d) unrelaxed 197-atom cluster: the guard (default g2:4:best) fires on the device, energy unchanged
# (e) GFN1 with -gpu_strict true: M applied on the device, exit code 0, same energy, fewer SCF iterations
# (f) a second frag run gives the same gradient file (md5; -xtb.rocm_fixed_point default true)
# Skipped (exit 77) when the run in (a) reports no usable ROCm device or plugin.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 23: GFN2 SCF preconditioner on the ROCm resident loop"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.0000001"

# CPU reference (cli_sqm_12, release build of commit 1632241c, -threads 1): urea3-water8.gfn2-relaxed.
REF_UREAWATER="-83.04383863"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
iters_of()  { strip "$1" | grep -oP 'SCF converged in \K[0-9]+' | tail -1; }

run() {   # run <log> <xyz> [args...]: GPU single point, verbosity 2 (prints the SCF path and M lines)
    local LOG=$1; local XYZ=$2; shift 2
    $CURCUMA -sp "$XYZ" -method gfn2 -gpu rocm -threads 4 -verbosity 2 -no_bmt "$@" > "$LOG" 2>&1
    echo "EXIT $?" >> "$LOG"
}

check_fewer() {   # check_fewer <n_frag> <n_none> <label>
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$1" ] && [ -n "$2" ] && [ "$1" -lt "$2" ]; then
        echo -e "${GREEN}PASS${NC}: $3 ($1 < $2)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (frag '$1', none '$2')"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # (a) default GPU path; decides about the skip
    run none.log ureawater.xyz
    if ! strip none.log | grep -q "fully device-resident loop"; then
        if strip none.log | grep -qiE "no usable|GPU fallbacks during this run|requires a ROCm build|plugin"; then
            echo -e "${YELLOW}SKIP${NC}: no usable ROCm device or plugin (resident GPU loop not active)"
            exit 77
        fi
        echo -e "${RED}FAIL${NC}: -gpu rocm did not run the device-resident loop and reported no missing device"
        exit 1
    fi
    assert_string_in_file "EXIT 0" none.log "default GPU run: exit code 0"
    local E0=$(energy_of none.log); local N0=$(iters_of none.log)
    assert_numeric_match "$REF_UREAWATER" "$E0" "$ENERGY_TOLERANCE" "default GPU run: energy vs CPU reference"
    assert_string_not_in_file "SCF precond" none.log "default GPU run: no preconditioner set up"

    # (b) frag on the device
    run frag.log ureawater.xyz -scf_precond frag -dump_gradient frag1.grad
    assert_string_in_file "EXIT 0" frag.log "frag GPU run: exit code 0"
    assert_string_in_file "SCF precond: device M uploaded to the resident GPU loop" frag.log "frag GPU run: M applied on the device"
    local E1=$(energy_of frag.log); local N1=$(iters_of frag.log)
    assert_numeric_match "$E0" "$E1" "$ENERGY_TOLERANCE" "frag GPU run: energy equals the run without the switch"
    check_fewer "$N1" "$N0" "frag GPU run: fewer SCF iterations"

    # (c) CPU frag gives the same energy
    $CURCUMA -sp ureawater.xyz -method gfn2 -threads 1 -verbosity 1 -no_bmt -scf_precond frag > cpu_frag.log 2>&1
    assert_exit_code $? 0 "frag CPU run"
    assert_numeric_match "$(energy_of cpu_frag.log)" "$E1" "$ENERGY_TOLERANCE" "frag: GPU energy equals the CPU energy"

    # (d) guard on the device
    run u197_none.log ureawater197.xyz
    run u197_frag.log ureawater197.xyz -scf_precond frag
    assert_string_in_file "EXIT 0" u197_frag.log "unrelaxed 197, frag GPU run: exit code 0"
    TESTS_RUN=$((TESTS_RUN + 1))
    if strip u197_frag.log | grep -qE "guard 'g2:4:best' fired at iteration 4 .*device"; then
        echo -e "${GREEN}PASS${NC}: unrelaxed 197: device guard fired at iteration 4"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: unrelaxed 197: device guard did not fire at iteration 4"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    assert_numeric_match "$(energy_of u197_none.log)" "$(energy_of u197_frag.log)" "$ENERGY_TOLERANCE" \
        "unrelaxed 197: energy with the guard equals the run without the switch"

    # (e) GFN1: M applied on the device, no GPU fallback under -gpu_strict
    $CURCUMA -sp ureawater.xyz -method gfn1 -gpu rocm -gpu_strict true -threads 4 -verbosity 2 -no_bmt \
        -scf_precond frag > gfn1.log 2>&1
    assert_exit_code $? 0 "GFN1 GPU run with -gpu_strict true"
    assert_string_in_file "SCF precond: device M uploaded to the resident GPU loop" gfn1.log "GFN1: M applied on the device"
    assert_string_not_in_file "scf_precond 'frag' ignored" gfn1.log "GFN1: switch not ignored"
    $CURCUMA -sp ureawater.xyz -method gfn1 -gpu rocm -threads 4 -verbosity 2 -no_bmt \
        -scf_precond none > gfn1_none.log 2>&1
    assert_numeric_match "$(energy_of gfn1_none.log)" "$(energy_of gfn1.log)" "$ENERGY_TOLERANCE" \
        "GFN1: energy with -scf_precond frag equals the run without the switch"
    check_fewer "$(iters_of gfn1.log)" "$(iters_of gfn1_none.log)" "GFN1: fewer SCF iterations"

    # (f) run-to-run reproducibility of the frag gradient (fixed-point accumulation)
    run frag2.log ureawater.xyz -scf_precond frag -dump_gradient frag2.grad
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s frag1.grad ] && [ "$(md5sum < frag1.grad)" = "$(md5sum < frag2.grad)" ]; then
        echo -e "${GREEN}PASS${NC}: frag: second run gives the identical gradient file"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: frag: gradient files of two runs differ (or are missing)"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
