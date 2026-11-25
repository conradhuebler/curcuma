#!/bin/bash
# Test: -scf_criterion xtb switches the GFN2 D4 cutoffs to xtb 6.7.1 (ATM 40, pair 72 Bohr); explicit values win
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only (Known Issue #47).
#
# Structure: water8 cube (triples inside 25-40 Bohr) plus one water 40 A away (pairs at 60-79 Bohr), 27 atoms.
#   (a) no flag = explicit 25/50 (tblite default unchanged), reported as 25/50
#   (b) -scf_criterion xtb reports 40/72 and equals -scf_criterion xtb with explicit 40/72
#   (c) -scf_criterion xtb with explicit 25/50 reports 25/50 and differs from (b) (the cutoffs act on this structure)
#   (d) an explicit pair cutoff alone keeps the automatic ATM 40
#   (e) without -scf_criterion xtb, explicit 40/72 differs from the default (each cutoff reaches pairs/triples here)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 23: D4 cutoffs follow -scf_criterion xtb unless given"
TEST_DIR="$SCRIPT_DIR"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { grep -oP '^# energy \K[-0-9.e]+' "$1" | head -1; }
cut_of()    { strip "$1" | grep -oP 'D4 cutoffs: pair \K[0-9.]+ Bohr, ATM [0-9.]+' | head -1; }
run() { local TAG=$1; shift; $CURCUMA -sp water9.xyz -method gfn2 -threads 1 -verbosity 3 -no_bmt -scf_threshold 1e-9 \
            -dump_gradient "$TAG.grad" "$@" > "$TAG.log" 2>&1; echo "EXIT $?" >> "$TAG.log"; }
check() {   # check <description> <condition>
    TESTS_RUN=$((TESTS_RUN + 1))
    if eval "$2"; then echo -e "${GREEN}PASS${NC}: $1"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else echo -e "${RED}FAIL${NC}: $1"; TESTS_FAILED=$((TESTS_FAILED + 1)); fi
}
absdiff() { python3 -c "import sys; print(abs(float(sys.argv[1]) - float(sys.argv[2])))" "$1" "$2"; }
gt() { python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) > float(sys.argv[2]) else 1)" "$1" "$2"; }

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f ./*.grad ./*.log

    run base
    run base_expl -xtb.d4_atm_cutoff 25 -xtb.d4_pair_cutoff 50
    run xtb       -scf_criterion xtb
    run xtb_expl  -scf_criterion xtb -xtb.d4_atm_cutoff 40 -xtb.d4_pair_cutoff 72
    run xtb_tbl   -scf_criterion xtb -xtb.d4_atm_cutoff 25 -xtb.d4_pair_cutoff 50
    run xtb_pair  -scf_criterion xtb -xtb.d4_pair_cutoff 50
    run def_wide  -xtb.d4_atm_cutoff 40 -xtb.d4_pair_cutoff 72
    for t in base base_expl xtb xtb_expl xtb_tbl xtb_pair def_wide; do
        assert_string_in_file "EXIT 0" $t.log "$t exits cleanly"
    done
    local EB=$(energy_of base.grad) EBE=$(energy_of base_expl.grad) EX=$(energy_of xtb.grad) EXE=$(energy_of xtb_expl.grad)
    local ET=$(energy_of xtb_tbl.grad) EW=$(energy_of def_wide.grad)

    check "(a) default reports pair 50 / ATM 25" '[ "$(cut_of base.log)" = "50.0 Bohr, ATM 25.0" ]'
    check "(a) default = explicit 25/50 (bitwise energy)" '[ "$EB" = "$EBE" ]'
    check "(b) scf_criterion xtb reports pair 72 / ATM 40" '[ "$(cut_of xtb.log)" = "72.0 Bohr, ATM 40.0" ]'
    check "(b) scf_criterion xtb = explicit 40/72 (bitwise energy)" '[ "$EX" = "$EXE" ]'
    check "(c) explicit 25/50 wins over scf_criterion xtb" '[ "$(cut_of xtb_tbl.log)" = "50.0 Bohr, ATM 25.0" ]'
    check "(c) the cutoffs change the energy here (|dE| > 1e-9 Eh)" 'gt "$(absdiff "$EX" "$ET")" 1e-9'
    check "(d) explicit pair cutoff alone keeps automatic ATM 40" '[ "$(cut_of xtb_pair.log)" = "50.0 Bohr, ATM 40.0" ]'
    check "(e) explicit 40/72 without the xtb rule differs from the default" 'gt "$(absdiff "$EW" "$EB")" 1e-9'

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
