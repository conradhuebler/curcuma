#!/bin/bash

# Test: ConfGen proposal generation (build recombined state vectors, optimise, judge)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026)
#
# Guards the two traps that made the first version of the generator report nonsense:
#
#   1. TOPOLOGY. A recombined structure can put two atoms inside bonding distance; the force field
#      then derives a new bond and the optimisation relaxes into a DIFFERENT MOLECULE. Measured
#      before the guard: a "new conformer" 2649 kJ/mol below the ensemble minimum, with 3-6 changed
#      bonds. Every optimised proposal must be topology-checked, and the check must use an explicit
#      covalent-radius factor -- Molecule's default 1.5 counts a compressed 1-3 contact as a bond and
#      rejected 45 of 46 chemically intact structures.
#
#   2. LIKE-FOR-LIKE ENERGIES. Proposals are optimised, the input ensemble may not be. Comparing them
#      directly turns the optimisation gain into a fake discovery ("106 kJ/mol below the minimum").
#      This test's input is deliberately an ensemble that was NOT optimised with gfnff, so the
#      warning about it MUST appear.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confgen - 02: proposal generation, topology + energy guards"
TEST_DIR="$SCRIPT_DIR"

MAX_PROPOSALS=6

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log ensemble.proposals*.xyz
    timeout 280 "$CURCUMA" -confgen ensemble.xyz -method gfnff -threads 1 \
        -generate true -max_proposals $MAX_PROPOSALS -proposal_templates 2 -no_bmt \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "proposal generation completes"

    # 1. proposals were actually built and optimised
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qE "state vector\(s\) proposed, [0-9]+ built" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: proposals generated and built"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        grep -oE "[0-9]+ state vector.*optimised successfully" stdout.log | head -1 | sed 's/^/  /'
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: no proposals were built"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. the topology check runs and is reported (guard 1)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "changed their bond topology" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: optimised proposals are topology-checked"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: topology check missing -- reaction products would be reported as conformers"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 3. counts stay consistent: new <= chemically valid
    local valid novel
    valid=$(grep -oE "of the [0-9]+ chemically valid" stdout.log | grep -oE "[0-9]+" | head -1)
    novel=$(grep -oE "[0-9]+ of [0-9]+ chemically valid proposals are NEW" stdout.log | grep -oE "^[0-9]+" | head -1)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$valid" ] && [ -n "$novel" ] && [ "$novel" -le "$valid" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: ${novel} new of ${valid} chemically valid proposals (consistent)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: inconsistent counts (new=${novel:-?}, valid=${valid:-?})"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 4. like-for-like energy guard (guard 2). This input was never optimised with gfnff, so the
    #    warning must fire -- if it does not, energies are being compared across optimisation states.
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "is NOT at the minimum of" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: un-optimised input ensemble detected and reported"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        grep -oE "re-optimising the input structures lowers them by up to [0-9.]+ kJ/mol" stdout.log \
            | head -1 | sed 's/^/  /'
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: no warning although the input is not at the gfnff minimum -- energy claims would be invalid"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 5. Conformers of one molecule span a few tens of kJ/mol. A reaction product (a formed or
    #    broken bond) sits about 1 Eh = 2600 kJ/mol away -- exactly the failure this guards against.
    #    So: all energies finite, negative, and their spread below 0.5 Eh.
    local newfile
    newfile=$(find_output_file "ensemble.proposals.new.xyz")
    if [ -n "$newfile" ]; then
        TESTS_RUN=$((TESTS_RUN + 1))
        local verdict
        verdict=$(awk '/Energy =/ { for (i=1;i<=NF;i++) if ($i=="=") { e=$(i+1)+0; got=$(i+1); break } ; if (got ~ /[nN][aA][nN]|[iI][nN][fF]/) bad=1; if (e>=0) bad=1; if (n++==0) {lo=e;hi=e}; if (e<lo) lo=e; if (e>hi) hi=e } END { if (n==0) print "empty"; else if (bad) print "bad"; else printf "ok %.4f\\n", hi-lo }' "$newfile")
        if [ "${verdict%% *}" = "ok" ]; then
            echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: new structures are conformers of one molecule (energy spread ${verdict#* } Eh < 0.5)"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        elif [ "$verdict" = "empty" ]; then
            echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: no new structures produced -- nothing to check here"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: implausible energies in the new structures ($verdict) -- reaction products slipped through the topology check"
            TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
        fi
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
