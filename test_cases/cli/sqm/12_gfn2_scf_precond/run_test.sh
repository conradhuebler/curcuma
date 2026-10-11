#!/bin/bash
# Test: native GFN2 SCF preconditioner (-scf_precond frag), docs/SCF_PRECOND.md
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026)
#
# (a) default (no switch): energy equals the reference of the binary before the preconditioner existed
# (b) -scf_precond frag: same energy as (a) within 1e-7 Eh and fewer SCF iterations
# (c) error paths: a charged system and a cluster of isolated ions ignore the switch with a warning
#     and give the same energy and iteration count as without it.
# (d) -scf_precond frag, 5 fs NVE MD on the unrelaxed 197-atom cluster (the guard fires there): runs without
#     error, and the first MD energy (same geometry, same fixed point) equals the run without the switch within
#     1e-7 Eh. No iteration condition: frag does not improve the cold start on this structure (measured n=1,
#     -sp -threads 8: none 24, frag 33 iterations; docs/SCF_PRECOND.md).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 12: GFN2 SCF preconditioner (scf_precond frag)"
TEST_DIR="$SCRIPT_DIR"

# Tolerance 1e-7 Eh: the preconditioner changes only the path to the SCF fixed point; the energy may
# move by the SCF threshold (1e-5 on max|dq|), which is far below this. Output resolution is 1e-8 Eh.
ENERGY_TOLERANCE="0.0000001"

# Reference values: curcuma release build of commit 1632241c (before -scf_precond existed),
# -method gfn2 -threads 1, measured 2026-10-10 (docs/SCF_PRECOND.md).
REF_UREAWATER="-83.04383863"       # urea3-water8.gfn2-relaxed, 11 SCF iterations there
REF_HYDRONIUM="-15.31425657"       # hydronium-water2, -charge 1
REF_NACL="-9.49239348"             # na2cl2-isolated-ions

energy_of() { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
iters_of()  { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'SCF converged in \K[0-9]+' | tail -1; }

run() {   # run <log> <xyz> [args...]
    local LOG=$1; local XYZ=$2; shift 2
    $CURCUMA -sp "$XYZ" -method gfn2 -threads 1 -verbosity 1 -no_bmt "$@" > "$LOG" 2>&1
    assert_exit_code $? 0 "curcuma -sp $XYZ $*"
}

first_energy() { sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'Final Energy: \K[-0-9.]+' | head -1; }

run_md() {   # run_md <name> [args...]: 5 fs NVE GFN2 MD, stdout and stderr to <name>.log
    local NAME=$1; shift
    $CURCUMA -md ureawater197.xyz -method gfn2 -max_time 5 -time_step 1 -temperature 300 \
        -thermostat none -seed 42 -threads 1 -verbosity 2 -no_bmt "$@" > "$NAME.log" 2>&1
    assert_exit_code $? 0 "curcuma -md ureawater197.xyz $*"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # (a) default path
    run none.log ureawater.xyz
    local E0=$(energy_of none.log); local N0=$(iters_of none.log)
    assert_numeric_match "$REF_UREAWATER" "$E0" "$ENERGY_TOLERANCE" "default: energy vs pre-change reference"
    assert_string_not_in_file "SCF precond" none.log "default: no preconditioner set up"

    # (b) fragment preconditioner
    run frag.log ureawater.xyz -scf_precond frag
    local E1=$(energy_of frag.log); local N1=$(iters_of frag.log)
    assert_numeric_match "$E0" "$E1" "$ENERGY_TOLERANCE" "frag: energy equals the default run"
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$N0" ] && [ -n "$N1" ] && [ "$N1" -lt "$N0" ]; then
        echo -e "${GREEN}PASS${NC}: frag: fewer SCF iterations ($N1 < $N0)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: frag: SCF iterations not reduced (frag '$N1', default '$N0')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # (c1) charged system: switch ignored with a warning
    run ion_none.log hydronium.xyz -charge 1
    run ion_frag.log hydronium.xyz -charge 1 -scf_precond frag
    assert_string_in_file "scf_precond 'frag' ignored: total charge" ion_frag.log "charged: warning printed"
    assert_numeric_match "$REF_HYDRONIUM" "$(energy_of ion_frag.log)" "$ENERGY_TOLERANCE" "charged: energy unchanged"
    assert_equals "$(iters_of ion_none.log)" "$(iters_of ion_frag.log)" "charged: same SCF iterations as without the switch"

    # (c2) isolated ions (single atoms): no block qualifies, switch ignored with a warning
    run nacl_frag.log nacl.xyz -scf_precond frag
    assert_string_in_file "no fragment qualifies" nacl_frag.log "isolated ions: warning printed"
    assert_numeric_match "$REF_NACL" "$(energy_of nacl_frag.log)" "$ENERGY_TOLERANCE" "isolated ions: energy unchanged"

    # (d) MD on a structure where the guard fires: no error, same first energy
    run_md md_none
    run_md md_frag -scf_precond frag
    assert_string_in_file "SCF precond: frag, 39 fragments" md_frag.log "MD frag: preconditioner built (39 fragments)"
    local EF0=$(first_energy md_none.log); local EF1=$(first_energy md_frag.log)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$EF0" ] && [ -n "$EF1" ]; then
        echo -e "${GREEN}PASS${NC}: MD frag: first energies found ($EF0, $EF1)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: MD frag: first energy missing (none '$EF0', frag '$EF1')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    assert_numeric_match "$EF0" "$EF1" "$ENERGY_TOLERANCE" "MD frag: first MD energy equals the run without the switch"

    # (e) CLI guard and spin/parity (Oct 2026): "-spin 1" as ONE argument is an error, not a silently
    #     ignored key; a radical without -spin gets a parity warning; with -spin 1 the energy equals
    #     xtb 6.7.1 (Cl2-, doublet: -9.321669164 Eh, identical with and without --uhf 1 in xtb).
    local REF_CL2="-9.32166916"
    $CURCUMA -sp cl2.xyz -method gfn2 -threads 1 -verbosity 1 -no_bmt -charge -1 "-spin 1" > spin_token.log 2>&1
    assert_exit_code $? 1 "flag token with whitespace (\"-spin 1\" as one argument) is refused"
    assert_string_in_file "contains whitespace" spin_token.log "flag token: error message printed"
    run cl2_nospin.log cl2.xyz -charge -1
    assert_string_in_file "inconsistent parity" cl2_nospin.log "radical without -spin: parity warning printed"
    run cl2_spin.log cl2.xyz -charge -1 -spin 1
    assert_string_not_in_file "inconsistent parity" cl2_spin.log "radical with -spin 1: no parity warning"
    assert_numeric_match "$REF_CL2" "$(energy_of cl2_spin.log)" "$ENERGY_TOLERANCE" "Cl2- with -spin 1 equals xtb 6.7.1"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
