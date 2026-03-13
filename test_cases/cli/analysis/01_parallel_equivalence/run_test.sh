#!/bin/bash
# Test: Analysis parallel vs sequential equivalence
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Frame-level parallelization correctness

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="analysis - 01: Parallel vs Sequential equivalence (100 frames, scattering)"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Sequential run (threads=1)
    $CURCUMA -analysis input.vtf -threads 1 -scattering_enable > seq_stdout.log 2> seq_stderr.log
    assert_exit_code $? 0 "Sequential analysis should succeed"

    # Rename sequential outputs
    mv input.general.csv seq.general.csv
    mv input.scattering_statistics.csv seq.scattering_statistics.csv

    # Parallel run (threads=4)
    $CURCUMA -analysis input.vtf -threads 4 -scattering_enable > par_stdout.log 2> par_stderr.log
    assert_exit_code $? 0 "Parallel analysis should succeed"

    # Rename parallel outputs
    mv input.general.csv par.general.csv
    mv input.scattering_statistics.csv par.scattering_statistics.csv

    return 0
}

validate_results() {
    cd "$TEST_DIR"

    # ---------------------------------------------------------------
    # Golden references (scnp_100frames.vtf, 200 CG beads, 100 frames)
    # Extracted from validated sequential run — January 2026
    # ---------------------------------------------------------------
    # general.csv: Frame 1   Gyr_u=20.0401  Rout=42.771
    # general.csv: Frame 100 Gyr_u=21.2418  Rout=38.0658
    # scattering_statistics.csv (N=100):
    #   first q-point (0.01 Å⁻¹): P_avg=0.984916  S_avg=196.973
    # ---------------------------------------------------------------

    # --- Check: output files exist ---
    assert_file_exists seq.general.csv "Sequential general.csv"
    assert_file_exists par.general.csv "Parallel general.csv"
    assert_file_exists seq.scattering_statistics.csv "Sequential scattering_statistics.csv"
    assert_file_exists par.scattering_statistics.csv "Parallel scattering_statistics.csv"

    # --- Check: exactly 101 lines (header + 100 frames) ---
    TESTS_RUN=$((TESTS_RUN + 1))
    local seq_lines=$(wc -l < seq.general.csv)
    if [ "$seq_lines" -eq 101 ]; then
        echo -e "${GREEN}✓ PASS${NC}: general.csv has 101 lines (header + 100 frames)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: general.csv has $seq_lines lines (expected 101)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # --- Check: N=100 in scattering header ---
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "N=100" seq.scattering_statistics.csv; then
        echo -e "${GREEN}✓ PASS${NC}: scattering_statistics.csv reports N=100 frames"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: scattering_statistics.csv does not report N=100"
        head -2 seq.scattering_statistics.csv
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # --- Golden: Frame 1 Gyration radius (unweighted) ---
    local gyr1=$(awk -F',' 'NR==2 {print $2}' seq.general.csv)
    assert_scientific_value "20.0401" "$gyr1" "0.001" "Frame 1 Gyr_u (golden ref)"

    # --- Golden: Frame 1 Rout ---
    local rout1=$(awk -F',' 'NR==2 {print $4}' seq.general.csv)
    assert_scientific_value "42.771" "$rout1" "0.01" "Frame 1 Rout (golden ref)"

    # --- Golden: Frame 100 Gyration radius ---
    local gyr100=$(awk -F',' 'NR==101 {print $2}' seq.general.csv)
    assert_scientific_value "21.2418" "$gyr100" "0.001" "Frame 100 Gyr_u (golden ref)"

    # --- Golden: Frame 100 Rout ---
    local rout100=$(awk -F',' 'NR==101 {print $4}' seq.general.csv)
    assert_scientific_value "38.0658" "$rout100" "0.01" "Frame 100 Rout (golden ref)"

    # --- Golden: first q-point P(q) average (q=0.01 Å⁻¹) ---
    local pavg=$(awk -F',' '/^[0-9]/{print $2; exit}' seq.scattering_statistics.csv)
    assert_scientific_value "0.984916" "$pavg" "0.001" "P(q) avg at q=0.01 (golden ref)"

    # --- Golden: first q-point S(q) average ---
    local savg=$(awk -F',' '/^[0-9]/{print $5; exit}' seq.scattering_statistics.csv)
    assert_scientific_value "196.973" "$savg" "0.1" "S(q) avg at q=0.01 (golden ref)"

    # --- Equivalence: general.csv identical seq vs par ---
    TESTS_RUN=$((TESTS_RUN + 1))
    if diff -q seq.general.csv par.general.csv > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASS${NC}: general.csv identical (seq vs par)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: general.csv differs between sequential and parallel"
        diff seq.general.csv par.general.csv | head -10
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # --- Equivalence: scattering_statistics.csv identical seq vs par ---
    TESTS_RUN=$((TESTS_RUN + 1))
    if diff -q seq.scattering_statistics.csv par.scattering_statistics.csv > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASS${NC}: scattering_statistics.csv identical (seq vs par)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: scattering_statistics.csv differs between sequential and parallel"
        diff seq.scattering_statistics.csv par.scattering_statistics.csv | head -10
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # --- Equivalence: parallel stdout confirms threaded path taken ---
    assert_string_in_file "threads" par_stdout.log "Parallel path reported in stdout"

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    rm -f seq.general.csv par.general.csv
    rm -f seq.scattering_statistics.csv par.scattering_statistics.csv
    rm -f seq_stdout.log seq_stderr.log par_stdout.log par_stderr.log
    rm -f input.general.csv input.scattering_statistics.csv
    rm -f input.scattering.gnu
}

main() {
    test_header "$TEST_NAME"
    cleanup_before

    if run_test; then
        validate_results
    fi

    print_test_summary

    if [ $TESTS_FAILED -gt 0 ]; then
        exit 1
    fi
    exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    main "$@"
fi
