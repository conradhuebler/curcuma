#!/bin/bash
# Test: Scattering analysis q-spacing, gnuplot script and statistics CSV
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026) - replaces the former loose script test_scattering_features.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="analysis - 02: Scattering q-spacing, gnuplot script, statistics CSV"
TEST_DIR="$SCRIPT_DIR"

# One frame of a 90-atom host; 20 q points between 0.01 and 2.0 A^-1.
SCATTER_ARGS=(-analysis input_aa.xyz -scattering_enable -scattering_q_min 0.01 -scattering_q_max 2.0 -scattering_q_steps 20 -no_bmt)

run_test() {
    cd "$TEST_DIR"

    # Default spacing is logarithmic
    $CURCUMA "${SCATTER_ARGS[@]}" > log_stdout.log 2> log_stderr.log
    assert_exit_code $? 0 "Scattering with default (logarithmic) spacing should succeed"
    mv input_aa.scattering_statistics.csv log.scattering_statistics.csv
    mv input_aa.scattering.gnu log.scattering.gnu

    # Linear spacing, separate argument tokens
    $CURCUMA "${SCATTER_ARGS[@]}" -scattering_q_spacing linear > lin_stdout.log 2> lin_stderr.log
    assert_exit_code $? 0 "Scattering with linear spacing should succeed"
    mv input_aa.scattering_statistics.csv lin.scattering_statistics.csv

    return 0
}

# q value of data row N (1-based, after the two header lines) of a statistics CSV
q_of_row() {
    sed -n "$((2 + $2))p" "$1" | cut -d, -f1
}

validate_results() {
    cd "$TEST_DIR"

    assert_string_in_file "logarithmic q-spacing" log_stdout.log "Default spacing reported as logarithmic"
    assert_string_in_file "linear q-spacing" lin_stdout.log "Linear spacing reported when requested"
    assert_file_exists log.scattering.gnu "Gnuplot script written"
    assert_file_exists log.scattering_statistics.csv "Statistics CSV (logarithmic)"
    assert_file_exists lin.scattering_statistics.csv "Statistics CSV (linear)"

    # 2 header lines + 20 q points
    TESTS_RUN=$((TESTS_RUN + 1))
    local nlog nlin
    nlog=$(wc -l < log.scattering_statistics.csv)
    nlin=$(wc -l < lin.scattering_statistics.csv)
    if [ "$nlog" -eq 22 ] && [ "$nlin" -eq 22 ]; then
        echo -e "${GREEN}✓ PASS${NC}: both statistics files have 20 q points"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: expected 22 lines, got $nlog (log) and $nlin (linear)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # q grid: first 0.01, last 2.0 for both spacings; second point 0.01*(200)^(1/19) = 0.013216 (log)
    # and 0.01 + 1.99/19 = 0.114737 (linear)
    assert_equals "0.010000" "$(q_of_row log.scattering_statistics.csv 1)" "First q (log)"
    assert_equals "2.000000" "$(q_of_row log.scattering_statistics.csv 20)" "Last q (log)"
    assert_equals "0.013216" "$(q_of_row log.scattering_statistics.csv 2)" "Second q (log)"
    assert_equals "0.010000" "$(q_of_row lin.scattering_statistics.csv 1)" "First q (linear)"
    assert_equals "2.000000" "$(q_of_row lin.scattering_statistics.csv 20)" "Last q (linear)"
    assert_equals "0.114737" "$(q_of_row lin.scattering_statistics.csv 2)" "Second q (linear)"

    # The script must plot the statistics file it belongs to; plotting itself needs gnuplot
    assert_string_in_file "scattering" log.scattering.gnu "Gnuplot script refers to scattering data"
    if command -v gnuplot > /dev/null 2>&1; then
        rm -f input_aa.scattering_plot.png
        cp log.scattering.gnu input_aa.scattering.gnu
        cp log.scattering_statistics.csv input_aa.scattering_statistics.csv
        TESTS_RUN=$((TESTS_RUN + 1))
        if gnuplot input_aa.scattering.gnu > gnuplot.log 2>&1 && [ -s input_aa.scattering_plot.png ]; then
            echo -e "${GREEN}✓ PASS${NC}: gnuplot renders the script"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}✗ FAIL${NC}: gnuplot did not produce input_aa.scattering_plot.png"
            TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    else
        echo "gnuplot not installed: rendering of the plot script not checked"
    fi

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    rm -f log.* lin.* log_*.log lin_*.log gnuplot.log
    rm -f input_aa.general.csv input_aa.scattering_statistics.csv input_aa.scattering.gnu input_aa.scattering_plot.png
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
