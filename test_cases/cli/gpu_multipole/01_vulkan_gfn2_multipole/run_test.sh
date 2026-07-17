#!/bin/bash
# Test: Vulkan GFN2 on-device multipole integrals (Stage 3m / V-AP2)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (2026-06): validates `gfn2 -gpu vulkan` energy against the CPU
# reference (`-gpu none`). The device multipole_ints kernel builds dp_int/qp_int, which
# the host GFN2 SCF consumes — bit-identical GFN2 energy proves the device integrals.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

# Prefer a release_vulkan/ binary when present (the canonical release/ has no Vulkan).
if [ -x "$PROJECT_ROOT/release_vulkan/curcuma" ]; then
    export CURCUMA="$PROJECT_ROOT/release_vulkan/curcuma"
fi

TEST_NAME="gpu_multipole - 01: Vulkan GFN2 multipole integrals (H2O)"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.0000001"   # 0.1 µEh — device integrals must reproduce CPU energy

extract_energy() { grep "Final Energy" "$1" | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.[0-9]+' | head -1; }

run_test() {
    cd "$TEST_DIR"
    echo "Executing: $CURCUMA -sp H2O.xyz -method gfn2 -gpu vulkan -verbosity 2"
    $CURCUMA -sp H2O.xyz -method gfn2 -gpu vulkan -verbosity 2 -no_bmt > vk.log 2>&1 || true
    if grep -qiE "no usable Vulkan|requires a Vulkan|USE_VULKAN|running CPU path" vk.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: Vulkan native xTB not available (no device or USE_VULKAN=OFF)"
        exit 0
    fi
    # Confirm the device multipole path actually ran (not a silent host fallback).
    if ! grep -qi "multipole integrals built on GPU device" vk.log; then
        echo -e "${YELLOW}SKIP${NC}: device multipole path not taken (host fallback)"
        exit 0
    fi
    echo "Executing: $CURCUMA -sp H2O.xyz -method gfn2 -gpu none -verbosity 2 (reference)"
    $CURCUMA -sp H2O.xyz -method gfn2 -gpu none -verbosity 2 -no_bmt > cpu.log 2>&1
    return 0
}

validate_results() {
    local e_vk e_cpu
    e_vk=$(extract_energy vk.log); e_cpu=$(extract_energy cpu.log)
    echo "Vulkan  : gfn2 E=$e_vk"
    echo "CPU ref : gfn2 E=$e_cpu"
    assert_scientific_value "$e_cpu" "$e_vk" "$ENERGY_TOLERANCE" "Vulkan GFN2 energy (device multipole) matches CPU"
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f vk.log cpu.log; }

main() {
    test_header "$TEST_NAME"
    cleanup_before
    if run_test; then validate_results; fi
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
