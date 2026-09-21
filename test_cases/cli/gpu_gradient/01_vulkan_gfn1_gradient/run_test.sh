#!/bin/bash
# Test: Vulkan GFN1 on-device nuclear gradient (Stage 4 / V-AP1)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (2026-06): validates `gfn1 -gpu vulkan` energy + gradient norm against
# the CPU reference (`-gpu none`). The device-resident GFN1 gradient (repulsion / Coulomb /
# H0-Pulay+CN gather kernels) must match the CPU path bit-for-bit. Mirrors CUDA gpu_gradient.

set -e
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

# The canonical release/ build usually has no Vulkan; prefer a release_vulkan/ binary
# when present so this test actually exercises the device path (otherwise it self-skips).
if [ -x "$PROJECT_ROOT/release_vulkan/curcuma" ]; then
    export CURCUMA="$PROJECT_ROOT/release_vulkan/curcuma"
fi

TEST_NAME="gpu_gradient - 01: Vulkan GFN1 nuclear gradient (NH3)"
TEST_DIR="$SCRIPT_DIR"

ENERGY_TOLERANCE="0.0000001"   # 0.1 µEh — device gradient must be bit-identical to CPU
GRADNORM_TOLERANCE="0.0000005"

extract_energy()  { grep "Final Energy" "$1" | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[-+]?[0-9]+\.[0-9]+' | head -1; }
extract_gradnorm() { grep "Gradient norm" "$1" | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[0-9]+\.[0-9]+[eE][-+][0-9]+' | head -1; }

run_test() {
    cd "$TEST_DIR"
    echo "Executing: $CURCUMA -sp NH3.xyz -method gfn1 -gpu vulkan -verbosity 2"
    $CURCUMA -sp NH3.xyz -method gfn1 -gpu vulkan -verbosity 2 -no_bmt > vk.log 2>&1 || true

    # Graceful skip if Vulkan is not compiled in or no FP64 device is present.
    if grep -qiE "no usable Vulkan|requires a Vulkan|USE_VULKAN|running CPU path" vk.log 2>/dev/null; then
        echo -e "${YELLOW}SKIP${NC}: Vulkan native xTB not available (no device or USE_VULKAN=OFF)"
        exit 0
    fi

    echo "Executing: $CURCUMA -sp NH3.xyz -method gfn1 -gpu none -verbosity 2 (reference)"
    $CURCUMA -sp NH3.xyz -method gfn1 -gpu none -verbosity 2 -no_bmt > cpu.log 2>&1
    return 0
}

validate_results() {
    local e_vk e_cpu g_vk g_cpu
    e_vk=$(extract_energy vk.log);    e_cpu=$(extract_energy cpu.log)
    g_vk=$(extract_gradnorm vk.log);  g_cpu=$(extract_gradnorm cpu.log)
    echo "Vulkan  : E=$e_vk  gradnorm=$g_vk"
    echo "CPU ref : E=$e_cpu  gradnorm=$g_cpu"
    assert_scientific_value "$e_cpu" "$e_vk" "$ENERGY_TOLERANCE" "Vulkan GFN1 energy matches CPU"
    assert_scientific_value "$g_cpu" "$g_vk" "$GRADNORM_TOLERANCE" "Vulkan GFN1 gradient norm matches CPU"
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
