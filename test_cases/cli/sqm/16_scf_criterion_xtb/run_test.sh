#!/bin/bash
# Test: -scf_criterion xtb (xtb 6.7.1 convergence rule) and -scf_accuracy on native gfn2 / gfn1
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only; the rule itself is checked against the xtb
# source (scc_core.f90, scf_module.F90), not against an xtb run, by this test.
#
#   (a) -scf_criterion default gives the same energy as no flag (nothing changes by default)
#   (b) -scf_criterion xtb converges and its energy is within 1e-6 Eh of the default rule
#   (c) -scf_accuracy 0.01 (tight) needs more iterations than -scf_accuracy 1 (the loose side is bounded below by the
#       FP32 phase of mixed precision, which never accepts convergence)
#   (d) the flat flag and -xtb.scf_criterion give the same energy and iteration count
#   (e) an unknown name is an error (non-zero exit)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 16: -scf_criterion xtb and -scf_accuracy"
TEST_DIR="$SCRIPT_DIR"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }
iters_of()  { strip "$1" | grep -oP 'SCF converged in \K[0-9]+' | head -1; }
run() { local LOG=$1; shift; $CURCUMA -sp benzene.xyz -threads 1 -verbosity 3 -no_bmt "$@" > "$LOG" 2>&1; echo "EXIT $?" >> "$LOG"; }

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    for METHOD in gfn2 gfn1; do
        run "${METHOD}_base.log" -method "$METHOD"
        run "${METHOD}_def.log" -method "$METHOD" -scf_criterion default
        run "${METHOD}_xtb.log" -method "$METHOD" -scf_criterion xtb
        run "${METHOD}_tight.log" -method "$METHOD" -scf_criterion xtb -scf_accuracy 0.01
        run "${METHOD}_dotted.log" -method "$METHOD" -xtb.scf_criterion xtb
        local EB=$(energy_of ${METHOD}_base.log); local ED=$(energy_of ${METHOD}_def.log); local EX=$(energy_of ${METHOD}_xtb.log)
        local EDOT=$(energy_of ${METHOD}_dotted.log)
        local IX=$(iters_of ${METHOD}_xtb.log); local IL=$(iters_of ${METHOD}_tight.log); local IDOT=$(iters_of ${METHOD}_dotted.log)

        assert_numeric_match "$ED" "$EB" "0.000000000001" "$METHOD: -scf_criterion default = no flag"
        assert_string_in_file "EXIT 0" ${METHOD}_xtb.log "$METHOD: -scf_criterion xtb exits cleanly"
        assert_numeric_match "$EX" "$EB" "0.000001" "$METHOD: xtb rule energy within 1e-6 Eh of the default rule"
        assert_numeric_match "$EDOT" "$EX" "0.000000000001" "$METHOD: -xtb.scf_criterion = flat flag (energy)"
        TESTS_RUN=$((TESTS_RUN + 1))
        if [ -n "$IX" ] && [ "$IX" = "$IDOT" ]; then
            echo -e "${GREEN}PASS${NC}: $METHOD: -xtb.scf_criterion = flat flag (iterations $IX)"; TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC}: $METHOD: iterations flat $IX, dotted $IDOT"; TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
        TESTS_RUN=$((TESTS_RUN + 1))
        if [ -n "$IL" ] && [ -n "$IX" ] && [ "$IL" -gt "$IX" ]; then
            echo -e "${GREEN}PASS${NC}: $METHOD: -scf_accuracy 0.01 needs more iterations ($IL > $IX)"; TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC}: $METHOD: tight accuracy iterations $IL, normal $IX"; TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    done

    run bad.log -method gfn2 -scf_criterion nonsense
    TESTS_RUN=$((TESTS_RUN + 1))
    if ! grep -q "EXIT 0" bad.log; then
        echo -e "${GREEN}PASS${NC}: unknown -scf_criterion is an error"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: unknown -scf_criterion was accepted"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
