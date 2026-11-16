#!/bin/bash
# NEB-MD CLI Test - Claude Generated 2026 (AI-generated, machine-tested)
# Drives a short NEB-like MD on ethane (staggered -> eclipsed) with GFN-FF and
# checks that the band runs, produces a path snapshot + energy table, and stays
# finite (no NaN/Inf). Uses BMT (default) so each run is isolated from stale
# GFN-FF parameter caches.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -nebmd start.xyz end.xyz \
        -nebmd.nimages 6 \
        -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.max_time 20 \
        -simplemd.temperature 300 -nebmd.k_spring 0.005 -nebmd.dump_frequency 10 \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD should succeed"
    local path_xyz=$(find_output_file "start.neb.path.xyz")
    assert_file_exists "$path_xyz" "NEB path snapshot"
    local energy_csv=$(find_output_file "start.neb.energy.csv")
    assert_file_exists "$energy_csv" "NEB energy table"
    return 0
}

validate_results() {
    local path_xyz=$(find_output_file "start.neb.path.xyz")
    local energy_csv=$(find_output_file "start.neb.energy.csv")
    [ -z "$path_xyz" ] || [ ! -f "$path_xyz" ] && return 1
    [ -z "$energy_csv" ] || [ ! -f "$energy_csv" ] && return 1

    # The path snapshot must contain at least one band block (6 images).
    local nblocks=$(grep -c "image=" "$path_xyz" 2>/dev/null || true)
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "${nblocks:-0}" -ge 6 ]; then
        echo -e "${GREEN}ok PASS${NC}: path.xyz has $nblocks image blocks (>= 6)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}xx FAIL${NC}: Expected >= 6 image blocks, got ${nblocks:-0}"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # Energy table must have a header + at least the step-0 rows (6 images).
    local nrows=$(wc -l < "$energy_csv")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$nrows" -ge 7 ]; then
        echo -e "${GREEN}ok PASS${NC}: energy.csv has $nrows lines (>= 7)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}xx FAIL${NC}: Expected >= 7 lines, got $nrows"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # No NaN/Inf anywhere in the outputs.
    local nan_count=$(grep -ciE "nan|inf" "$energy_csv" "$path_xyz" 2>/dev/null | awk -F: '{s+=$NF} END{print s+0}')
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "${nan_count:-0}" -eq 0 ]; then
        echo -e "${GREEN}ok PASS${NC}: no NaN/Inf in outputs"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}xx FAIL${NC}: Found $nan_count NaN/Inf occurrences"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "NEB-MD Test (gfnff, ethane)"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi