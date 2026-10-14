#!/bin/bash
# Test: ROCm GFN-FF single point of a Mn complex equals the CPU single point
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026): the D4 Gaussian weights of the device (k_gaussian_weights) used a
# `norm > 1e-10` switch that collapsed atoms with far-away reference CNs (Mn, Ru in carbonyl
# complexes) onto the first reference, so the dispersion term differed from the CPU by 2.4e-3 Eh
# for this structure (MOR41 PR05). The CPU energy is computed in the same test, no stored value.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="gfnff_gpu - 04: ROCm GFN-FF single point of a Mn complex (MOR41 PR05) against the CPU"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.00001"   # 10 uEh; the unfixed kernel differs by 2.4e-3 Eh

energy_of() {
    sed 's/\x1b\[[0-9;]*m//g' "$1" | grep "Single Point Energy" | grep -oE '[-+]?[0-9]+\.[0-9]+' | head -1
}

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -sp mol.xyz -method gfnff -gpu none -no_bmt -verbosity 1 > cpu.log 2> cpu.err
    assert_exit_code $? 0 "CPU single point should succeed"
    rm -f mol.topo.json
    $CURCUMA -sp mol.xyz -method gfnff -gpu rocm -no_bmt -verbosity 1 > rocm.log 2> rocm.err || true

    if grep -q "without ROCm GFN-FF\|USE_ROCM\|using CPU implementation\|no HIP device\|hipErrorNoDevice\|HIP error\|no usable HIP\|found no GPU plugin" rocm.log rocm.err 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: ROCm GFN-FF not available (build or device)"
        exit 0
    fi
    return 0
}

validate_results() {
    local e_cpu e_gpu
    e_cpu=$(energy_of cpu.log)
    e_gpu=$(energy_of rocm.log)
    if [ -z "$e_cpu" ] || [ -z "$e_gpu" ]; then
        echo -e "${RED}✗ FAIL${NC}: could not extract the energies (cpu='$e_cpu' rocm='$e_gpu')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
    echo "CPU:  $e_cpu Eh"
    echo "ROCm: $e_gpu Eh"
    assert_scientific_value "$e_cpu" "$e_gpu" "$ENERGY_TOLERANCE" "ROCm gfnff energy of the Mn complex matches the CPU"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    rm -f cpu.log cpu.err rocm.log rocm.err mol.topo.json
    if run_test; then
        validate_results
    fi
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
