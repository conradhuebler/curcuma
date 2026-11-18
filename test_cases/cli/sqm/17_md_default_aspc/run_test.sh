#!/bin/bash
# Test: -md with native gfn1/gfn2 uses ASPC charge extrapolation by default; an explicit flag wins
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only.
#
#   (a) -md -method gfn2 / gfn1: the SCF log reports "aspc extrapolation" (default since Nov 18, 2026)
#   (b) -scf_extrapolation none: no extrapolation line (explicit choice wins, also via -xtb.scf_extrapolation)
#   (c) -scf_extrapolation gauss: the log names gauss (explicit choice wins)
#   (d) -sp has no history: no extrapolation line
#   (e) the final energy of the MD with the default equals the one with 'none' within 1e-6 Eh (same fixpoint)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 17: MD default SCC extrapolation (aspc)"
TEST_DIR="$SCRIPT_DIR"

strip() { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
md()  { local LOG=$1; shift; mkdir -p "d_$LOG"; (cd "d_$LOG" && cp ../benzene.xyz . && $CURCUMA -md benzene.xyz -threads 1 -max_time 6 -seed 3 -thermostat none -verbosity 3 -no_bmt "$@" > "../$LOG" 2>&1; echo "EXIT $?" >> "../$LOG"); }
epot_last() { strip "$1" | grep -E '^ +0\.006000 ' | head -1 | awk '{print $2}'; }

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    for METHOD in gfn2 gfn1; do
        md "${METHOD}_def.log" -method "$METHOD"
        md "${METHOD}_none.log" -method "$METHOD" -scf_extrapolation none
        md "${METHOD}_dot.log" -method "$METHOD" -xtb.scf_extrapolation none
        md "${METHOD}_gauss.log" -method "$METHOD" -scf_extrapolation gauss
        assert_string_in_file "aspc extrapolation over" ${METHOD}_def.log "$METHOD: MD uses aspc by default"
        assert_string_not_in_file "extrapolation over" ${METHOD}_none.log "$METHOD: -scf_extrapolation none switches it off"
        assert_string_not_in_file "extrapolation over" ${METHOD}_dot.log "$METHOD: -xtb.scf_extrapolation none switches it off"
        assert_string_in_file "gauss extrapolation over" ${METHOD}_gauss.log "$METHOD: -scf_extrapolation gauss is respected"
        local ED=$(epot_last ${METHOD}_def.log); local EN=$(epot_last ${METHOD}_none.log)
        assert_numeric_match "$ED" "$EN" "0.000001" "$METHOD: MD energy with the default equals 'none'"
    done

    mkdir -p d_sp; (cd d_sp && cp ../benzene.xyz . && $CURCUMA -sp benzene.xyz -method gfn2 -threads 1 -verbosity 3 -no_bmt > ../sp.log 2>&1)
    assert_string_not_in_file "extrapolation over" sp.log "-sp: no extrapolation (no history)"

    rm -rf d_*
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
