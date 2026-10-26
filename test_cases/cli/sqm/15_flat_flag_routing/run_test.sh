#!/bin/bash
# Test: a flat flag several modules register is routed by -method (-electronic_temperature)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026). AI-generated, machine-tested only.
#
# -electronic_temperature is registered by tblite, ulysses and xtb. Until Oct 16, 2026 the flat flag stayed in the
# command module with a warning, so no method read it: -method gfn2 -electronic_temperature 40000 gave the 300 K
# energy. CLI2Json now routes it to the module of the selected method (MethodFactory::parameterScopeForMethod).
#   (a) gfn2 / gfn1: the flat flag gives the same energy as -xtb.electronic_temperature
#   (b) that energy differs from the 300 K energy by more than 0.1 Eh (the smearing took effect)
#   (c) without -method (gfnff) nothing can be routed: the ambiguity warning is still printed

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 15: flat -electronic_temperature routed by -method"
TEST_DIR="$SCRIPT_DIR"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
run() { local LOG=$1; shift; $CURCUMA -sp benzene.xyz -threads 1 -verbosity 1 -no_bmt "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    for METHOD in gfn2 gfn1; do
        run "${METHOD}_base.log" -method "$METHOD"
        run "${METHOD}_flat.log" -method "$METHOD" -electronic_temperature 40000
        run "${METHOD}_dotted.log" -method "$METHOD" -xtb.electronic_temperature 40000
        local EB=$(energy_of ${METHOD}_base.log); local EF=$(energy_of ${METHOD}_flat.log); local ED=$(energy_of ${METHOD}_dotted.log)
        assert_numeric_match "$ED" "$EF" "0.000000001" "$METHOD: flat flag equals -xtb.electronic_temperature"
        TESTS_RUN=$((TESTS_RUN + 1))
        if python3 -c "import sys; sys.exit(0 if abs(float('$EF') - float('$EB')) > 0.1 else 1)"; then
            echo -e "${GREEN}PASS${NC}: $METHOD: 40000 K changes the energy ($EB -> $EF)"; TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC}: $METHOD: 40000 K did not change the energy ($EB -> $EF)"; TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    done

    run gfnff_flat.log -electronic_temperature 40000
    assert_string_in_file "registered in multiple modules" gfnff_flat.log "no -method: the ambiguity warning is still printed"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
