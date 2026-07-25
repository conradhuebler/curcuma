#!/bin/bash

# Test: ConfGen ensemble analysis (torsion states + matched-pair energy decomposition)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026)
#
# Validates the analysis stage of -confgen on a real 44-structure ensemble (114 atoms):
#   1. the run completes
#   2. rotatable torsions are detected from the topology
#   3. the GFN-FF energy decomposition SUMS TO THE TOTAL ENERGY -- checked numerically from the
#      written CSV, not from a log line. This is the regression guard for the missing "Repulsion"
#      term (Jul 2026): it was absent from GFNFFComputationalMethod::getEnergyDecomposition(), which
#      left the components off by ~2781 kJ/mol on a 90-atom molecule and made every per-term
#      attribution meaningless.
#   4. matched pairs are found and the three CSV tables are written

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="confgen - 01: matched-pair analysis + decomposition completeness"
TEST_DIR="$SCRIPT_DIR"

MIN_TORSIONS=10
SUM_TOL=1e-6          # Hartree, term sum vs total energy

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log ensemble.torsions.csv ensemble.torsion_states.csv ensemble.matched_pairs.csv
    timeout 200 "$CURCUMA" -confgen ensemble.xyz -method gfnff -threads 1 -no_bmt \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "ConfGen analysis completes"

    # 1. torsions detected
    local n_torsions
    n_torsions=$(grep -oE "[0-9]+ rotatable torsion" stdout.log | grep -oE "^[0-9]+" | head -1)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "${n_torsions:-0}" -ge "$MIN_TORSIONS" ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: ${n_torsions} rotatable torsions detected (>= $MIN_TORSIONS)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: only ${n_torsions:-0} rotatable torsions detected"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 2. the three tables exist and carry data
    for f in ensemble.torsions.csv ensemble.torsion_states.csv ensemble.matched_pairs.csv; do
        local path
        path=$(find_output_file "$f")
        TESTS_RUN=$((TESTS_RUN + 1))
        if [ -n "$path" ] && [ "$(grep -cv '^#' "$path" 2>/dev/null)" -ge 1 ]; then
            echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: $f written with data rows"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: $f missing or empty"
            TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
        fi
    done

    # 3. decomposition completeness, computed from the CSV: sum of all *_Eh term columns must equal
    #    the energy_Eh column for every structure.
    local table max_dev
    table=$(find_output_file "ensemble.torsions.csv")
    max_dev=$(awk -F, '
        NR==1 {
            for (i = 1; i <= NF; i++) {
                gsub(/^[ \t#]+|[ \t]+$/, "", $i)
                if ($i == "energy_Eh") { etot = i }
                else if ($i ~ /_Eh$/)  { term[++nterm] = i }
            }
            next
        }
        {
            s = 0
            for (t = 1; t <= nterm; t++) s += $(term[t])
            d = s - $(etot); if (d < 0) d = -d
            if (d > worst) worst = d
        }
        END { printf "%.10f", worst }' "$table")
    TESTS_RUN=$((TESTS_RUN + 1))
    if awk -v d="$max_dev" -v t="$SUM_TOL" 'BEGIN{exit(d<t?0:1)}'; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: energy decomposition sums to the total energy for all structures (worst deviation $max_dev Eh)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: decomposition incomplete -- worst deviation $max_dev Eh > $SUM_TOL (a term is missing from getEnergyDecomposition)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 4. matched pairs were actually measured
    local n_pairs_rows
    n_pairs_rows=$(grep -cv '^#' "$(find_output_file ensemble.matched_pairs.csv)" 2>/dev/null)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "${n_pairs_rows:-0}" -ge 3 ]; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: ${n_pairs_rows} state transitions measured from matched pairs"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: only ${n_pairs_rows:-0} transitions found (expected >= 3)"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # 5. the additivity check is reported (it is the decision criterion for the recombination stage)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "additivity check" stdout.log; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: additivity check reported"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        grep -oE "additivity check.*" stdout.log | head -1 | sed 's/^/  /'
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: additivity check missing from the report"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
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
