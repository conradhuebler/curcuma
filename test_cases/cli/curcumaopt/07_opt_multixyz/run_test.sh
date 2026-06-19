#!/bin/bash
# Test: Multi-XYZ optimisation with gfnff optimises all frames
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="curcumaopt - 07: Multi-XYZ gfnff optimisation"
TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"

    # Execute: optimise all 17 frames of the helicen multi-XYZ file.
    # Default max_iterations is used so the test validates converged minima.
    $CURCUMA -opt helicen.xyz -method gfnff -threads 1 -no_bmt > stdout.log 2> stderr.log
    local exit_code=$?

    # The command may return 0 even when individual frames do not converge;
    # the important result is that the output file contains all frames.
    local output_file=$(find_output_file "helicen.opt.xyz")
    assert_file_exists "$output_file" "gfnff multi-XYZ optimisation produced output file"

    return 0
}

validate_results() {
    local output_file=$(find_output_file "helicen.opt.xyz")

    if [ -z "$output_file" ] || [ ! -f "$output_file" ]; then
        echo -e "${RED}✗ FAIL${NC}: helicen.opt.xyz missing"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    local frame_count
    frame_count=$(count_xyz_structures "$output_file")

    # The input file contains 17 frames; all of them must be optimised and written out.
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$frame_count" -eq 17 ]; then
        echo -e "${GREEN}✓ PASS${NC}: Output contains all 17 optimised structures"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Expected 17 structures in helicen.opt.xyz, found $frame_count"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Verify that at least one energy value was written and looks physically reasonable.
    local first_energy
    first_energy=$(extract_energy_from_xyz "$output_file")
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -n "$first_energy" ] && (( $(echo "$first_energy < 0" | bc -l) )); then
        echo -e "${GREEN}✓ PASS${NC}: First optimised energy $first_energy Eh is negative (physically reasonable)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Could not extract a reasonable energy from $output_file"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Compare every optimised frame against the golden reference obtained by
    # optimising each frame individually. This ensures the multi-XYZ path reaches
    # the same minima as the single-structure path.
    local ref_file="$TEST_DIR/golden_energies.txt"
    if [ -f "$ref_file" ]; then
        local idx=0
        local all_match=true
        local tolerance=0.00001
        while IFS= read -r line || [ -n "$line" ]; do
            local ref_idx ref_energy
            read -r ref_idx ref_energy <<< "$line"
            local actual_energy
            actual_energy=$(extract_nth_energy_from_xyz "$output_file" "$((idx + 1))")
            TESTS_RUN=$((TESTS_RUN + 1))
            if [ -n "$actual_energy" ] && [ -n "$ref_energy" ]; then
                local diff
                diff=$(awk -v a="$actual_energy" -v b="$ref_energy" 'BEGIN { d=a-b; if (d<0) d=-d; print d }')
                local passes
                passes=$(awk -v d="$diff" -v t="$tolerance" 'BEGIN { if (d<=t) print 1; else print 0 }')
                if [ "$passes" -eq 1 ]; then
                    echo -e "${GREEN}✓ PASS${NC}: Frame $ref_idx energy matches reference ($actual_energy vs $ref_energy Eh, diff ${diff})"
                    TESTS_PASSED=$((TESTS_PASSED + 1))
                else
                    echo -e "${RED}✗ FAIL${NC}: Frame $ref_idx energy mismatch ($actual_energy vs $ref_energy Eh, diff ${diff} > $tolerance)"
                    TESTS_FAILED=$((TESTS_FAILED + 1))
                    all_match=false
                fi
            else
                echo -e "${RED}✗ FAIL${NC}: Could not extract energy for frame $ref_idx"
                TESTS_FAILED=$((TESTS_FAILED + 1))
                all_match=false
            fi
            idx=$((idx + 1))
        done < "$ref_file"

        if [ "$all_match" = true ]; then
            echo -e "${GREEN}✓ PASS${NC}: All optimised energies match individual-structure references within $tolerance Eh"
        fi
    else
        echo -e "${YELLOW}⚠ WARNING${NC}: Golden reference file $ref_file not found; skipping per-frame energy check"
    fi

    return 0
}

cleanup_before() {
    cd "$TEST_DIR"
    cleanup_test_artifacts
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
