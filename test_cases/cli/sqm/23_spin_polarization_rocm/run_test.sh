#!/bin/bash
# Test: spin polarisation (-xtb.spin_polarization) of native GFN1/GFN2 on the ROCm path against the CPU
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026). AI-generated, machine-tested only (docs/SPIN_POLARISATION.md, section GPU).
# ROCm twin of cli_sqm_22_spin_polarization_cuda (same cases, same limits).
#
# The CPU spin-polarised SCF is the reference (it matches tblite 0.8.0, cli_sqm_16). Checked on -gpu rocm:
#   - triplet O2 (tilted bond, -spin 2) GFN2 and GFN1, and Pt7 (-spin 8, 2000 K) GFN2 and GFN1, FP64 on the GPU
#     (-scf_mixed_precision false): energy within 1e-7 Eh, every gradient component within 1e-6 Eh/A of the CPU,
#   - mixed precision (FP32 steps far from convergence) on O2 and Pt7, GFN2: energy within 1e-6 Eh,
#   - the run used the spin branch of the device loop (log line); -gpu_strict true makes a GPU fallback an error,
#   - -scf_gpu_partial_diag, -scf_precond frag and -large_system_mode fragments with spin polarisation end with an
#     error, not a silent run.
# Skipped (exit 77) when the GPU run reports no usable HIP device or no ROCm plugin.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 23: spin polarisation, GFN1/GFN2 on the ROCm path against the CPU (O2, Pt7)"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.0000001"
GRADIENT_TOLERANCE="0.000001"
MIXED_TOLERANCE="0.000001"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }

run() {   # run <log> <grad> <xyz> <method> <spin> [args...]
    local LOG=$1; local GRAD=$2; local XYZ=$3; local METHOD=$4; local SPIN=$5; shift 5
    $CURCUMA -sp "$XYZ" -method "$METHOD" -spin "$SPIN" -xtb.spin_polarization true -gradient -dump_gradient "$GRAD" \
        -threads 1 -scf_threshold 1e-8 -scf_max_iterations 500 -verbosity 2 -no_bmt "$@" > "$LOG" 2>&1
    echo "EXIT $?" >> "$LOG"
}

check_gradient() {   # check_gradient <cpu grad> <gpu grad> <label>
    TESTS_RUN=$((TESTS_RUN + 1))
    local MAXDIFF
    MAXDIFF=$(python3 - "$1" "$2" <<'PY'
import sys
a = [float(x) for l in open(sys.argv[1]) if not l.startswith('#') for x in l.split()]
b = [float(x) for l in open(sys.argv[2]) if not l.startswith('#') for x in l.split()]
print(max(abs(x - y) for x, y in zip(a, b)) if a and len(a) == len(b) else 'nan')
PY
)
    if python3 -c "import sys; d=float('$MAXDIFF'); sys.exit(0 if d == d and d < $GRADIENT_TOLERANCE else 1)"; then
        echo -e "${GREEN}PASS${NC}: $3 (max |dgrad| $MAXDIFF < $GRADIENT_TOLERANCE)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (max |dgrad| $MAXDIFF)"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # <file> <spin> <extra arguments>
    local CASES=("o2.xyz|2|" "pt7.xyz|8|-xtb.electronic_temperature 2000")
    for CASE in "${CASES[@]}"; do
        IFS='|' read -r XYZ SPIN EXTRA <<< "$CASE"
        for METHOD in gfn2 gfn1; do
            local TAG="${XYZ%.xyz}_${METHOD}"
            run "${TAG}_cpu.log" "${TAG}_cpu.grad" "$XYZ" "$METHOD" "$SPIN" $EXTRA
            run "${TAG}_gpu.log" "${TAG}_gpu.grad" "$XYZ" "$METHOD" "$SPIN" -gpu rocm -gpu_strict true \
                -scf_mixed_precision false $EXTRA
            if strip "${TAG}_gpu.log" | grep -qiE "no usable HIP|no HIP device|hipErrorNoDevice|requires a ROCm build|plugin .*not (found|loaded)|cannot load"; then
                echo -e "${YELLOW}SKIP${NC}: no usable HIP device or ROCm plugin"
                exit 77
            fi
            assert_string_in_file "EXIT 0" "${TAG}_cpu.log" "$TAG CPU run: exit code 0"
            assert_string_in_file "EXIT 0" "${TAG}_gpu.log" "$TAG GPU run (FP64): exit code 0"
            assert_string_in_file "spin polarisation on the host-driven device loop" "${TAG}_gpu.log" \
                "$TAG GPU run: spin branch of the device loop"
            assert_numeric_match "$(energy_of ${TAG}_cpu.log)" "$(energy_of ${TAG}_gpu.log)" "$ENERGY_TOLERANCE" \
                "$TAG: GPU energy (FP64) equals the CPU energy"
            check_gradient "${TAG}_cpu.grad" "${TAG}_gpu.grad" "$TAG: GPU gradient (FP64) equals the CPU gradient"
            if [ "$METHOD" = gfn2 ]; then
                run "${TAG}_mixed.log" "${TAG}_mixed.grad" "$XYZ" "$METHOD" "$SPIN" -gpu rocm -gpu_strict true \
                    -scf_mixed_precision true $EXTRA
                assert_string_in_file "EXIT 0" "${TAG}_mixed.log" "$TAG GPU run (mixed precision): exit code 0"
                assert_numeric_match "$(energy_of ${TAG}_cpu.log)" "$(energy_of ${TAG}_mixed.log)" "$MIXED_TOLERANCE" \
                    "$TAG: GPU energy (mixed precision) equals the CPU energy"
            fi
        done
    done

    # Unsupported combinations on the GPU: a clear error that names the flag, and no energy. The grep requires the
    # [ERROR] line, because the echoed command line also contains both words.
    local NEG
    for NEG in "scf_gpu_partial_diag true" "scf_precond frag" "large_system_mode fragments"; do
        local TAG="neg_${NEG// /_}"
        local FLAG="${NEG%% *}"
        run "${TAG}.log" "${TAG}.grad" o2.xyz gfn2 2 -gpu rocm -${NEG}
        TESTS_RUN=$((TESTS_RUN + 1))
        if strip "${TAG}.log" | grep -q "ERROR.*spin_polarization.*${FLAG}" && ! strip "${TAG}.log" | grep -q "Single Point Energy" \
           && ! grep -q "EXIT 0" "${TAG}.log"; then
            echo -e "${GREEN}PASS${NC}: -gpu rocm -${NEG} with spin polarisation ends with an error and no energy"; TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC}: -gpu rocm -${NEG} with spin polarisation did not end with the expected error"; TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    done

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
