#!/bin/bash
# Test: ROCm GFN-FF Single Point — energy matches CPU gfnff
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (June 2026): Validates gfnff -gpu rocm CLI against CPU gfnff reference

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="gfnff_gpu - 02: ROCm GFN-FF Single Point (caffeine)"
TEST_DIR="$SCRIPT_DIR"

# CPU gfnff reference (same value as the CUDA test 01). ROCm must match within FP.
REFERENCE_TOTAL_ENERGY="-4.672736999448"
ENERGY_TOLERANCE="0.00002"   # 20 µEh — generous for GPU floating-point reduction order

run_test() {
    cd "$TEST_DIR"

    echo "Executing: $CURCUMA -sp caffeine.xyz -method gfnff -gpu rocm -verbosity 1"
    $CURCUMA -sp caffeine.xyz -method gfnff -gpu rocm -verbosity 1 > stdout.log 2> stderr.log
    local exit_code=$?

    # Graceful skip if ROCm GFN-FF is not compiled in
    if grep -q "without ROCm GFN-FF\|USE_ROCM\|using CPU implementation" stderr.log stdout.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: gfnff ROCm not available in this build (USE_ROCM=OFF)"
        exit 0
    fi

    # Graceful skip if no HIP device is found at runtime
    if grep -q "no HIP device\|hipErrorNoDevice\|HIP error\|no usable HIP" stderr.log stdout.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: No ROCm-capable GPU found at runtime"
        exit 0
    fi

    assert_exit_code $exit_code 0 "gfnff ROCm single point should succeed"
    return 0
}

validate_results() {
    local total_energy=""
    if grep -q "\[ENERGY\]Total Energy:" stdout.log 2>/dev/null; then
        total_energy=$(grep "\[ENERGY\]Total Energy:" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    elif grep -q "Single Point Energy" stdout.log 2>/dev/null; then
        total_energy=$(grep "Single Point Energy" stdout.log | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    fi

    if [ -z "$total_energy" ]; then
        echo -e "${RED}✗ FAIL${NC}: Could not extract energy from gfnff ROCm output"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    echo "ROCm gfnff energy: $total_energy Eh"
    echo "CPU gfnff ref:     $REFERENCE_TOTAL_ENERGY Eh"
    assert_scientific_value "$REFERENCE_TOTAL_ENERGY" "$total_energy" "$ENERGY_TOLERANCE" "gfnff ROCm energy matches CPU reference"
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
