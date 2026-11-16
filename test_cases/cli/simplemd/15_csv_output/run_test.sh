#!/bin/bash
# SimpleMD CSV status output test - Claude Generated (September 2026)
# Verifies that per-step energies/temperatures are written to a UTF-8 CSV
# file (<basename>.md.csv) in the BMT output directory, semicolon-delimited.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # print=100 gives several CSV rows for a 1000 fs run (dt=1 fs default)
    $CURCUMA -md input.xyz -method uff -md.max_time 1000 -md.print 100 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "MD should succeed"
    local csv_file=$(find_output_file "input.md.csv")
    assert_file_exists "$csv_file" "CSV status file in BMT dir"
    return 0
}

validate_results() {
    local csv_file=$(find_output_file "input.md.csv")
    [ -z "$csv_file" ] || [ ! -f "$csv_file" ] && return 1

    # UTF-8 BOM: first 3 bytes must be EF BB BF
    local bom_hex=$(head -c 3 "$csv_file" | od -An -tx1 | tr -d ' \n')
    if [ "$bom_hex" = "efbbbf" ]; then
        echo -e "${GREEN}✓ PASS${NC}: UTF-8 BOM present"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: UTF-8 BOM missing (got: $bom_hex)"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # Header row with semicolon delimiter (BOM precedes "step", so not anchored)
    assert_string_in_file "step;time_ps;epot;epot_avg;ekin;ekin_avg;etot;etot_avg;temperature;temperature_avg" \
        "$csv_file" "CSV header with semicolon delimiter"

    # At least one data row (line starting with an integer step + delimiter)
    if grep -qE '^[0-9]+;' "$csv_file"; then
        echo -e "${GREEN}✓ PASS${NC}: CSV has data rows"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: CSV has no data rows"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # No comma-delimited header (delimiter must be configurable, default ';')
    assert_string_not_in_file "step,time_ps,epot" "$csv_file" "CSV not comma-delimited by default"
    return 0
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "SimpleMD CSV Output Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
