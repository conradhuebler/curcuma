#!/bin/bash
# Test: solvent names (case, xtb aliases), errors for unusable requests, reference states of the native solvation
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only.
#
#   (a) any case and the xtb aliases give the energy of the canonical name, and that energy is not the gas phase
#   (b) an unknown solvent, a GBSA solvent without parameters, cpcm with a state other than gsolv (tblite rejects it too)
#       and an unknown state are errors (exit code != 0)
#   (c) -solvent_state (one flat argument for gfn1, gfn2 and gfnff): bar1M and reference shift the energy as xtb 6.7.1 does
#       (reference values measured with xtb 6.7.1 --alpb <solvent> <state> on methanol, identical for gfn1, gfn2, gfnff):
#       bar1M - gsolv = +3.0312463360e-3 Eh; reference - gsolv: water +1.8137894620e-3, acetone -6.2090370e-4, octanol -2.0990886180e-3.
#       One shift for every method: ideal-gas volume 24.79 L/mol as xtb, kB from CODATA 2018 (xtb: 3e-9 Eh smaller in bar1M);
#       tolerance 1e-8 Eh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 18: solvent names and reference states"
TEST_DIR="$SCRIPT_DIR"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
run() { local LOG=$1; shift; $CURCUMA -sp methanol.xyz -threads 1 -verbosity 1 -no_bmt -dump_gradient "${LOG%.log}.grad" "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }
# full-precision energy (15 decimals) from the gradient dump header; the printed energy has 8 decimals only
eprec() { grep -oP '# energy \K[-0-9.e]+' "${1%.log}.grad" | head -1; }
exit_of() { grep -oP 'EXIT \K[0-9]+' "$1" | tail -1; }   # no line anchor: curcuma may end without a newline
differs_from_gas() {  # energy $1 lower than the gas phase energy $2 by more than 1 mEh
    python3 -c "import sys; sys.exit(0 if float('$2') - float('$1') > 1e-3 else 1)"
}
expect_error() {      # log, text that must appear, description
    TESTS_RUN=$((TESTS_RUN + 1))
    local rc=$(exit_of "$1")
    if [ -n "$rc" ] && [ "$rc" != "0" ] && strip "$1" | grep -q -i "$2"; then
        echo -e "${GREEN}PASS${NC}: $3 (exit code $rc)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (exit code $rc, expected non-zero and '$2')"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}
shift_check() {       # description, energy_state, energy_gsolv, expected shift, tolerance
    local d=$(python3 -c "print(float('$2') - float('$3'))")
    TESTS_RUN=$((TESTS_RUN + 1))
    if python3 -c "import sys; sys.exit(0 if abs(($d) - ($4)) < $5 else 1)"; then
        echo -e "${GREEN}PASS${NC}: $1 shift $d Eh (xtb $4)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $1 shift $d Eh, xtb $4, tolerance $5"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    run gas.log -method gfn2
    local EGAS=$(energy_of gas.log)
    run gasff.log -method gfnff
    local EGASFF=$(energy_of gasff.log)

    # (a) case and aliases, native gfn2 (alpb) and gfnff
    run_pair() {   # canonical name, other spelling, description
        run "can_${1// /_}.log" -method gfn2 -xtb.solvent "$1"
        run "alt_${2// /_}.log" -method gfn2 -xtb.solvent "$2"
        local EC=$(energy_of "can_${1// /_}.log") EA=$(energy_of "alt_${2// /_}.log")
        assert_numeric_match "$EC" "$EA" "0.000000000001" "gfn2: '$2' = '$1'"
        TESTS_RUN=$((TESTS_RUN + 1))
        if differs_from_gas "$EA" "$EGAS"; then
            echo -e "${GREEN}PASS${NC}: gfn2 '$2' is a solvent run, not the gas phase"; TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC}: gfn2 '$2' gave the gas phase energy"; TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    }
    run_pair water Water
    run_pair water H2O
    run_pair chcl3 CHCl3
    run_pair chcl3 chloroform
    run_pair ch2cl2 dichloromethane
    run_pair ch2cl2 DCM
    run_pair ether "diethyl ether"
    run_pair thf Tetrahydrofuran
    run_pair ethylacetate "ethyl acetate"
    run_pair woctanol "octanol wet"
    run_pair cs2 "carbon disulfide"
    run_pair acetonitrile methylcyanide
    run_pair furane furan
    run_pair hexane n-hexane

    run ff_water.log -method gfnff -gfnff.solvent water
    run ff_Water.log -method gfnff -gfnff.solvent Water
    run ff_dcm.log -method gfnff -gfnff.solvent DCM
    run ff_ch2cl2.log -method gfnff -gfnff.solvent ch2cl2
    assert_numeric_match "$(energy_of ff_water.log)" "$(energy_of ff_Water.log)" "0.000000000001" "gfnff: 'Water' = 'water'"
    assert_numeric_match "$(energy_of ff_ch2cl2.log)" "$(energy_of ff_dcm.log)" "0.000000000001" "gfnff: 'DCM' = 'ch2cl2'"
    TESTS_RUN=$((TESTS_RUN + 1))
    if differs_from_gas "$(energy_of ff_Water.log)" "$EGASFF"; then
        echo -e "${GREEN}PASS${NC}: gfnff 'Water' is a solvent run"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: gfnff 'Water' gave the gas phase energy"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # 'vacuum' and 'gas' mean no solvent (docs/SOLVATION.md)
    run vac.log -method gfn2 -xtb.solvent vacuum
    assert_numeric_match "$(energy_of vac.log)" "$EGAS" "0.000000000001" "gfn2: 'vacuum' = gas phase"
    assert_string_in_file "EXIT 0" vac.log "gfn2: 'vacuum' exits cleanly"
    run vacff.log -method gfnff -gfnff.solvent Vacuum
    assert_numeric_match "$(energy_of vacff.log)" "$EGASFF" "0.000000000001" "gfnff: 'Vacuum' = gas phase"

    # (b) errors
    run e1.log -method gfn2 -xtb.solvent watr
    expect_error e1.log "unknown solvent" "gfn2: unknown solvent is an error"
    run e2.log -method gfn1 -xtb.solvent dmf -xtb.solvent_model gbsa
    expect_error e2.log "unknown solvent" "gfn1 GBSA without a dmf parameter set is an error"
    run e3.log -method gfnff -gfnff.solvent watr
    expect_error e3.log "unknown solvent" "gfnff: unknown solvent is an error"
    run e4.log -method gfn2 -xtb.solvent water -xtb.solvent_model cpcm -solvent_state bar1M
    expect_error e4.log "solvent_state" "gfn2: cpcm with solvent_state bar1M is an error (as tblite --cpcm)"
    run e5.log -method gfn2 -xtb.solvent water -solvent_state foo
    expect_error e5.log "solvent_state" "gfn2: unknown solvent_state is an error"
    run e6.log -method gfnff -gfnff.solvent water -solvent_state foo
    expect_error e6.log "solvent_state" "gfnff: unknown solvent_state is an error"

    # (c) reference states
    for METHOD in gfn2 gfn1; do
        TOL="0.00000001"
        for SOL in water acetone octanol; do
            case $SOL in water) REF="0.0018137894620";; acetone) REF="-0.00062090370";; octanol) REF="-0.002099088618";; esac
            run "st_${METHOD}_${SOL}_g.log" -method $METHOD -xtb.solvent $SOL
            run "st_${METHOD}_${SOL}_b.log" -method $METHOD -xtb.solvent $SOL -solvent_state bar1M
            run "st_${METHOD}_${SOL}_r.log" -method $METHOD -xtb.solvent $SOL -solvent_state Reference
            local G=$(eprec "st_${METHOD}_${SOL}_g.log")
            shift_check "$METHOD $SOL bar1M" "$(eprec st_${METHOD}_${SOL}_b.log)" "$G" "0.0030312463360" "$TOL"
            shift_check "$METHOD $SOL reference" "$(eprec st_${METHOD}_${SOL}_r.log)" "$G" "$REF" "$TOL"
        done
    done
    for SOL in water acetone octanol; do
        case $SOL in water) REF="0.0018137894620";; acetone) REF="-0.00062090370";; octanol) REF="-0.002099088618";; esac
        run "st_ff_${SOL}_g.log" -method gfnff -gfnff.solvent $SOL
        run "st_ff_${SOL}_b.log" -method gfnff -gfnff.solvent $SOL -solvent_state BAR1M
        run "st_ff_${SOL}_r.log" -method gfnff -gfnff.solvent $SOL -solvent_state reference
        local G=$(eprec "st_ff_${SOL}_g.log")
        shift_check "gfnff $SOL bar1M" "$(eprec st_ff_${SOL}_b.log)" "$G" "0.0030312463360" "0.00000001"
        shift_check "gfnff $SOL reference" "$(eprec st_ff_${SOL}_r.log)" "$G" "$REF" "0.00000001"
    done

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
