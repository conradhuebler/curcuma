#!/bin/bash
# Test: native CPCM-X (-xtb.solvent_model cpcmx) for gfn1/gfn2 against xtb 6.7.1 --cpcmx
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only.
#
#   (a) total free energies equal xtb 6.7.1 (edcfbbe, /opt/bin/xtb) `xtb mol.xyz --gfn N --cpcmx <solvent> --sp --acc 0.01
#       [--chrg q]` (TOTAL ENERGY, measured Nov 22, 2026 at the library geometries) within 5e-8 Eh; measured differences
#       <= 3.0e-8 Eh (caffeine), which are the gas-phase differences of the native GFN methods to xtb
#   (b) -xtb.cpcmx_smd_oxygen_radius oxygen equals CPCM-X 1.1.0 (built from source, same xtb.cosmo input) within 5e-8 Eh
#   (c) errors: -opt, -dump_gradient (no gradient), unknown solvent, missing database, solvent_state other than gsolv
# The solvent files in cpcmx_db/ are DB/xtb/{water,hexane,acetonitrile}.cosmo of CPCM-X 1.1.0 (LGPL-3.0, see README.md).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 21: native CPCM-X solvation (xtb --cpcmx)"
TEST_DIR="$SCRIPT_DIR"
DB="$SCRIPT_DIR/cpcmx_db"

strip() { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
exit_of() { grep -oP 'EXIT \K[0-9]+' "$1" | tail -1; }
run() { local LOG=$1 XYZ=$2; shift 2; $CURCUMA -sp "$XYZ" -threads 1 -verbosity 1 -no_bmt "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }
total() { strip "$1" | grep -oP 'total free energy = \K[-0-9.]+' | tail -1; }

expect_error() {      # log, text that must appear, description
    TESTS_RUN=$((TESTS_RUN + 1))
    local rc=$(exit_of "$1")
    if [ -n "$rc" ] && [ "$rc" != "0" ] && strip "$1" | grep -q -i "$2"; then
        echo -e "${GREEN}PASS${NC}: $3 (exit code $rc)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (exit code $rc, expected non-zero and '$2')"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    check() {   # tag, xyz, reference energy, description, args...
        local TAG=$1 XYZ=$2 REF=$3 WHAT=$4; shift 4
        run "$TAG.log" "$XYZ" -scf_threshold 1e-10 -xtb.solvent_model cpcmx -xtb.cpcmx_db "$DB" "$@"
        assert_numeric_match "$REF" "$(total $TAG.log)" "0.00000005" "$TAG = $WHAT"
    }
    # (a) xtb 6.7.1 --cpcmx
    check w_g2_water   water.xyz      -5.078368201954  "xtb 6.7.1 --cpcmx" -method gfn2 -xtb.solvent water
    check m_g2_hexane  methanol.xyz   -8.229311027620  "xtb 6.7.1 --cpcmx" -method gfn2 -xtb.solvent hexane
    check c_g2_hexane  caffeine.xyz  -42.162128307375  "xtb 6.7.1 --cpcmx" -method gfn2 -xtb.solvent hexane
    check c_g1_water   caffeine.xyz  -44.534520176408  "xtb 6.7.1 --cpcmx" -method gfn1 -xtb.solvent water
    check oh_g1_water  hydroxyl.xyz   -5.529367739674  "xtb 6.7.1 --cpcmx" -method gfn1 -xtb.solvent water -charge -1
    check n_g2_acn     ammonium.xyz   -4.610227644859  "xtb 6.7.1 --cpcmx" -method gfn2 -xtb.solvent acetonitrile -charge 1
    check g_g2_acn     glucoside.xyz -46.543788947297  "xtb 6.7.1 --cpcmx" -method gfn2 -xtb.solvent acetonitrile

    # (b) CPCM-X 1.1.0 oxygen radius rule
    check c_g2_hex_o   caffeine.xyz  -42.1610416329638 "CPCM-X 1.1.0" -method gfn2 -xtb.solvent hexane -xtb.cpcmx_smd_oxygen_radius oxygen
    check g_g2_acn_o   glucoside.xyz -46.5398944056923 "CPCM-X 1.1.0" -method gfn2 -xtb.solvent acetonitrile -xtb.cpcmx_smd_oxygen_radius oxygen

    # (c) errors
    $CURCUMA -opt water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcmx -xtb.cpcmx_db "$DB" -no_bmt > e1.log 2>&1; echo "EXIT $?" >> e1.log
    expect_error e1.log "only -sp is supported" "cpcmx: -opt is rejected"
    run e2.log water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcmx -xtb.cpcmx_db "$DB" -dump_gradient e2.grad
    expect_error e2.log "no analytic gradient" "cpcmx: -sp -dump_gradient is rejected"
    run e3.log water.xyz -method gfn2 -xtb.solvent watr -xtb.solvent_model cpcmx -xtb.cpcmx_db "$DB"
    expect_error e3.log "Unknown CPCM-X solvent" "cpcmx: unknown solvent is an error"
    run e4.log water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcmx -xtb.cpcmx_db "$SCRIPT_DIR/no_such_dir"
    expect_error e4.log "database not found" "cpcmx: missing database directory is an error"
    run e5.log water.xyz -method gfn2 -xtb.solvent dmso -xtb.solvent_model cpcmx -xtb.cpcmx_db "$DB"
    expect_error e5.log "not in the database" "cpcmx: solvent file missing from the database is an error"
    run e6.log water.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcmx -xtb.cpcmx_db "$DB" -solvent_state bar1M
    expect_error e6.log "solvent_state" "cpcmx: solvent_state bar1M is an error"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
