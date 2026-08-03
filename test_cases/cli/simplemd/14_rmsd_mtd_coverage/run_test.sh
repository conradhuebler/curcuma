#!/bin/bash
# Test: RMSD-MTD strided coverage / hill spacing diagnostics
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Jul 2026) - docs/RMSD_MTD_TEXTBOOK.md section 6 (check 4) and 7.3.
# Validates that the strided scheme writes basename.mtd_coverage_statistics.csv and that the
# nearest-neighbour hill spacing is consistent with r_dep: r_dep resolves to the auto FWHM(alpha=10)
# ~0.527 A, and the largest nearest-neighbour spacing stays within a few r_dep (no large gaps).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 14: RMSD-MTD coverage / hill spacing"
TEST_DIR="$SCRIPT_DIR"
RDEP_EXPECT=0.5266   # 2.35482 / sqrt(2*10)
RDEP_TOL=0.02

run_test() {
    cd "$TEST_DIR"
    cleanup_bmt_dirs
    rm -f stdout.log stderr.log
    timeout 120 "$CURCUMA" -md input.xyz -method gfnff -rmsd_mtd \
        -maxtime 2000 -time_step 0.5 -temperature 300 -thermostat csvr \
        -md.rmsd_mtd_deposit_stride 10 -md.seed 7 -threads 1 \
        > stdout.log 2> stderr.log
    RUN_EXIT=$?
}

validate_results() {
    local failed=0

    assert_exit_code $RUN_EXIT 0 "strided MD completes"

    local stats
    stats=$(find_output_file "input.mtd_coverage_statistics.csv")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -z "$stats" ] || [ ! -s "$stats" ]; then
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: coverage statistics file missing (${stats:-not found})"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
    echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: coverage statistics written ($stats)"
    TESTS_PASSED=$((TESTS_PASSED + 1))

    # Data line: n,min,mean,max,r_dep,count_above_r_dep
    local line n mn mean mx rdep
    line=$(grep -vE "^#" "$stats" | tail -1)
    n=$(echo "$line" | cut -d, -f1)
    mn=$(echo "$line" | cut -d, -f2)
    mx=$(echo "$line" | cut -d, -f4)
    rdep=$(echo "$line" | cut -d, -f5)
    echo "  coverage: n=$n min=$mn max=$mx r_dep=$rdep"

    # r_dep resolves to the auto FWHM(alpha=10)
    assert_scientific_value "$RDEP_EXPECT" "$rdep" "$RDEP_TOL" "r_dep = auto FWHM(alpha=10)"

    # No large coverage gaps: max nearest-neighbour spacing within 3*r_dep
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$mx" ] && [ -n "$rdep" ] && awk -v m="$mx" -v r="$rdep" 'BEGIN{exit(m < 3.0*r ? 0 : 1)}'; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: max nearest-neighbour spacing $mx A < 3*r_dep (no large gaps)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: max nearest-neighbour spacing ${mx:-unknown} A not < 3*r_dep=${rdep}"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    # Pool is non-trivial (more than the single initial hill)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$n" ] && [ "$n" -ge 3 ] 2>/dev/null; then
        echo -e "${GREEN}\xe2\x9c\x93 PASS${NC}: $n hills deposited"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}\xe2\x9c\x97 FAIL${NC}: only ${n:-0} hills deposited"
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

main "$@"
