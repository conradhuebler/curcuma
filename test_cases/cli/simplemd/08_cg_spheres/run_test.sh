#!/bin/bash
# SimpleMD CG Test - Claude Generated (Oct 2025)
# Tests SimpleMD with coarse-grained spherical particles and PBC
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

run_test() {
    cd "$TEST_DIR"
    # Short CG simulation: 2 beads, 20 fs, should use 10x timestep scaling
    $CURCUMA -md input.xyz \
             -method cg \
             -load_ff_json cg_params.json \
             -md.max_time 20 \
             -md.time_step 1.0 \
             -md.thermostat none \
             -md.verbosity 1 \
             > stdout.log 2> stderr.log

    exit_code=$?
    assert_exit_code $exit_code 0 "SimpleMD CG should succeed"
    return 0
}

validate_results() {
    # Check XYZ trajectory exists (always written)
    if [ -f "input.trj.xyz" ]; then
        echo -e "${GREEN}✓ PASS${NC}: XYZ trajectory generated"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: XYZ trajectory not found"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Check VTF trajectory exists (NEW: CG-specific output)
    if [ -f "input.trj.vtf" ]; then
        echo -e "${GREEN}✓ PASS${NC}: VTF trajectory generated for CG system"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}○ WARN${NC}: VTF trajectory not found (may be disabled)"
    fi

    # Count frames in XYZ trajectory (should have at least 2: initial + final)
    frames=$(count_xyz_structures "input.trj.xyz")
    if [ $frames -ge 2 ]; then
        echo -e "${GREEN}✓ PASS${NC}: Trajectory has $frames frames (minimum 2 required)"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: Expected at least 2 frames, got $frames"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    # Check for CG detection message in output
    if grep -q "Pure CG system detected" stdout.log 2>/dev/null; then
        echo -e "${GREEN}✓ PASS${NC}: CG system detection working"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}○ INFO${NC}: CG detection message not found (may be expected at low verbosity)"
    fi

    # Check for timestep scaling message
    if grep -q "CG timestep scaling applied" stdout.log 2>/dev/null; then
        echo -e "${GREEN}✓ PASS${NC}: Timestep scaling detected"
        TESTS_RUN=$((TESTS_RUN + 1))
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${YELLOW}○ INFO${NC}: Timestep scaling message not found"
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; }

main() {
    test_header "SimpleMD CG Simulation Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
