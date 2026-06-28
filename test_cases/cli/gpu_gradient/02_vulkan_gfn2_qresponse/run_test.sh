#!/bin/bash
# Test: Vulkan GFN2 D4 EEQ charge-response gradient (Stage 5 Part A / X-AP4)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (2026-06): validates the device-resident GFN2 D4 q-response on Vulkan —
# the EEQ charges + ∂q/∂x adjoint solve now run on the GPU (d4eeq_cn/build/solve/resp shaders,
# single-workgroup no-pivot FP64 Gaussian elimination, derf erfcc). The full GFN2 D4 gradient is
# device-resident; energy + gradient norm must match the CPU reference (`-gpu none`). HCN is polar
# so the q-response is a non-trivial gradient contribution. Mirrors 01 (GFN1) + CUDA gpu_gradient.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

# Prefer a release_vulkan/ binary so this test actually exercises the device path.
if [ -x "$PROJECT_ROOT/release_vulkan/curcuma" ]; then
    export CURCUMA="$PROJECT_ROOT/release_vulkan/curcuma"
fi

TEST_NAME="gpu_gradient - 02: Vulkan GFN2 D4 q-response gradient (HCN)"
TEST_DIR="$SCRIPT_DIR"

ENERGY_TOLERANCE="0.0000001"     # 0.1 µEh — q-response does not affect the energy
GRADNORM_TOLERANCE="0.0000005"   # device q-response gradient must match CPU host model

extract_energy()   { grep "Final Energy" "$1" | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.[0-9]+' | head -1; }
extract_gradnorm() { grep "Gradient norm" "$1" | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[0-9]+\.[0-9]+[eE][-+][0-9]+' | head -1; }

run_test() {
    cd "$TEST_DIR"
    echo "Executing: $CURCUMA -sp HCN.xyz -method gfn2 -gpu vulkan -verbosity 2"
    $CURCUMA -sp HCN.xyz -method gfn2 -gpu vulkan -verbosity 2 -no_bmt > vk.log 2>&1 || true

    # Graceful skip if Vulkan is not compiled in or no FP64 device is present.
    if grep -qiE "no usable Vulkan|requires a Vulkan|USE_VULKAN|running CPU path" vk.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: Vulkan native xTB not available (no device or USE_VULKAN_XTB=OFF)"
        exit 0
    fi

    echo "Executing: $CURCUMA -sp HCN.xyz -method gfn2 -gpu none -verbosity 2 (reference)"
    $CURCUMA -sp HCN.xyz -method gfn2 -gpu none -verbosity 2 -no_bmt > cpu.log 2>&1
    return 0
}

validate_results() {
    local e_vk e_cpu g_vk g_cpu
    e_vk=$(extract_energy vk.log);    e_cpu=$(extract_energy cpu.log)
    g_vk=$(extract_gradnorm vk.log);  g_cpu=$(extract_gradnorm cpu.log)
    echo "Vulkan  : E=$e_vk  gradnorm=$g_vk"
    echo "CPU ref : E=$e_cpu  gradnorm=$g_cpu"
    assert_scientific_value "$e_cpu" "$e_vk" "$ENERGY_TOLERANCE" "Vulkan GFN2 energy matches CPU"
    assert_scientific_value "$g_cpu" "$g_vk" "$GRADNORM_TOLERANCE" "Vulkan GFN2 q-response gradient norm matches CPU"
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f vk.log cpu.log; }

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
