#!/bin/bash
# Test: spin polarisation (-xtb.spin_polarization) of native GFN1/GFN2 against tblite 0.8.0
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only (docs/SPIN_POLARISATION.md).
#
# Triplet O2 (tilted bond, so the gradient has three components) with the tblite spin-polarization container
# (Python wheel tblite 0.8.0, uhf = 2, 300 K, accuracy 0.001, spin constants of tblite v0.8.0). Checked:
#   - energy and the full gradient of GFN2 and GFN1 against the tblite numbers (energy 1e-7 Eh, gradient 1e-6 Eh/A;
#     measured at implementation: 3e-9 Eh, 1e-8 Eh/Bohr),
#   - the option off (default): -spin 2 gives the energy of the model without spin polarisation (tblite value),
#   - the option on with -spin 0 changes nothing (closed shell),
#   - the excluded SCF modes (-scf_mode diis, -scf_precond frag, -large_system_mode fragments) end with an error, not a silent fallback.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 16: spin polarisation of native GFN1/GFN2 (triplet O2) against tblite 0.8.0"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.0000001"
GRADIENT_TOLERANCE="0.000001"

# tblite 0.8.0, o2.xyz (library structure dioxygen.bent-triplet), uhf 2, 300 K, spin-polarization 1.0
REF_GFN2_SPIN="-7.929775795419304"
REF_GFN1_SPIN="-9.142703563503403"
REF_GFN2_NOSPIN="-7.90192423867321"
# gradient of atom 1 in Eh/Bohr (atom 2 is the negative)
REF_GFN2_G="-0.013069260544005215 -0.00435642018133503 -0.052712684194153786"
REF_GFN1_G="-0.01600743325862025 -0.005335811086206692 -0.06456331414310232"
BOHR_IN_ANGSTROM="0.52917721092"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }

run() {   # run <log> <grad> <method> <spin> [args...]
    local LOG=$1; local GRAD=$2; local METHOD=$3; local SPIN=$4; shift 4
    $CURCUMA -sp o2.xyz -method "$METHOD" -spin "$SPIN" -gradient -dump_gradient "$GRAD" -threads 1 \
        -scf_threshold 1e-8 -scf_max_iterations 500 -verbosity 2 -no_bmt "$@" > "$LOG" 2>&1
    echo "EXIT $?" >> "$LOG"
}

check_gradient() {   # check_gradient <curcuma grad file (Eh/A)> <ref atom 1 Eh/Bohr> <label>
    TESTS_RUN=$((TESTS_RUN + 1))
    local MAXDIFF
    MAXDIFF=$(python3 - "$1" "$2" "$BOHR_IN_ANGSTROM" <<'PY'
import sys
g = [float(x) for l in open(sys.argv[1]) if not l.startswith('#') for x in l.split()]
r = [float(x) / float(sys.argv[3]) for x in sys.argv[2].split()]
ref = r + [-x for x in r]
print(max(abs(x - y) for x, y in zip(g, ref)) if len(g) == 6 else 'nan')
PY
)
    if python3 -c "import sys; d=float('$MAXDIFF'); sys.exit(0 if d == d and d < $GRADIENT_TOLERANCE else 1)"; then
        echo -e "${GREEN}PASS${NC}: $3 (max |dgrad| $MAXDIFF Eh/A < $GRADIENT_TOLERANCE)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (max |dgrad| $MAXDIFF Eh/A)"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    run gfn2_spin.log gfn2_spin.grad gfn2 2 -xtb.spin_polarization true
    assert_string_in_file "EXIT 0" gfn2_spin.log "GFN2 triplet with spin polarisation: exit code 0"
    assert_numeric_match "$REF_GFN2_SPIN" "$(energy_of gfn2_spin.log)" "$ENERGY_TOLERANCE" "GFN2 triplet O2 energy against tblite"
    check_gradient gfn2_spin.grad "$REF_GFN2_G" "GFN2 triplet O2 gradient against tblite"

    run gfn1_spin.log gfn1_spin.grad gfn1 2 -xtb.spin_polarization true
    assert_string_in_file "EXIT 0" gfn1_spin.log "GFN1 triplet with spin polarisation: exit code 0"
    assert_numeric_match "$REF_GFN1_SPIN" "$(energy_of gfn1_spin.log)" "$ENERGY_TOLERANCE" "GFN1 triplet O2 energy against tblite"
    check_gradient gfn1_spin.grad "$REF_GFN1_G" "GFN1 triplet O2 gradient against tblite"

    # Option off: the model without spin polarisation, unchanged.
    run gfn2_off.log gfn2_off.grad gfn2 2
    assert_numeric_match "$REF_GFN2_NOSPIN" "$(energy_of gfn2_off.log)" "$ENERGY_TOLERANCE" "GFN2 triplet O2 without the option: energy of the unpolarised model"

    # Closed shell: the option has no effect.
    run gfn2_s0_off.log gfn2_s0_off.grad gfn2 0
    run gfn2_s0_on.log gfn2_s0_on.grad gfn2 0 -xtb.spin_polarization true
    assert_numeric_match "$(energy_of gfn2_s0_off.log)" "$(energy_of gfn2_s0_on.log)" "0.000000001" "GFN2 -spin 0: option on equals option off"

    # Excluded modes: a clear error, no energy.
    local NEG
    for NEG in "scf_mode diis" "scf_precond frag" "large_system_mode fragments"; do
        local TAG="neg_${NEG// /_}"
        run "${TAG}.log" "${TAG}.grad" gfn2 2 -xtb.spin_polarization true -${NEG}
        TESTS_RUN=$((TESTS_RUN + 1))
        if strip "${TAG}.log" | grep -q "spin_polarization" && ! strip "${TAG}.log" | grep -q "Single Point Energy" \
           && ! grep -q "EXIT 0" "${TAG}.log"; then
            echo -e "${GREEN}PASS${NC}: -${NEG} with spin polarisation ends with an error and no energy"; TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC}: -${NEG} with spin polarisation did not end with the expected error"; TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    done

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
