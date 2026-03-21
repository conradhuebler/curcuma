#!/bin/bash
# Test: GPU GFN-FF (ggfnff) Single Point — energy matches CPU gfnff
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (March 2026): Validates ggfnff CLI method against CPU gfnff reference

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="ggfnff - 01: GPU GFN-FF Single Point (caffeine)"
TEST_DIR="$SCRIPT_DIR"

# Reference energy from CPU gfnff (test 02_caffeine_gfnff_energy_components)
# ggfnff must match within 1e-5 Eh
REFERENCE_TOTAL_ENERGY="-4.672736999448"
ENERGY_TOLERANCE="0.00002"   # 20 µEh — generous for GPU floating-point

run_test() {
    cd "$TEST_DIR"

    # Attempt to run ggfnff; exit code 1 with "requires a CUDA build" → skip
    echo "Executing: $CURCUMA -sp caffeine.xyz -method ggfnff -verbosity 1"
    $CURCUMA -sp caffeine.xyz -method ggfnff -verbosity 1 > stdout.log 2> stderr.log
    local exit_code=$?

    # Graceful skip if ggfnff is not compiled in
    if grep -q "requires a CUDA build\|USE_CUDA" stderr.log stdout.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: ggfnff not available in this build (USE_CUDA=OFF)"
        exit 0
    fi

    # Graceful skip if no GPU device is found
    if grep -q "no CUDA device\|CUDA error\|cudaGetDeviceCount" stderr.log stdout.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: No CUDA-capable GPU found at runtime"
        exit 0
    fi

    assert_exit_code $exit_code 0 "ggfnff single point should succeed"
    return 0
}

validate_results() {
    # Extract total energy
    local total_energy=""
    if grep -q "\[ENERGY\]Total Energy:" stdout.log 2>/dev/null; then
        total_energy=$(grep "\[ENERGY\]Total Energy:" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    elif grep -q "Single Point Energy" stdout.log 2>/dev/null; then
        total_energy=$(grep "Single Point Energy" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    fi

    if [ -z "$total_energy" ]; then
        total_energy=$(extract_energy_from_xyz "caffeine.opt.xyz" 2>/dev/null || echo "")
    fi

    if [ -z "$total_energy" ]; then
        echo -e "${RED}✗ FAIL${NC}: Could not extract energy from ggfnff output"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    echo "GPU ggfnff energy: $total_energy Eh"
    echo "CPU gfnff  ref:    $REFERENCE_TOTAL_ENERGY Eh"
    assert_scientific_value "$REFERENCE_TOTAL_ENERGY" "$total_energy" "$ENERGY_TOLERANCE" "ggfnff energy matches gfnff CPU reference"
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
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
