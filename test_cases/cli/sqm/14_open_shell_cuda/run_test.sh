#!/bin/bash
# Test: open-shell GFN1/GFN2 on the CUDA path against the CPU (Known Issue #39)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026). AI-generated, machine-tested only.
#
# Until Oct 16, 2026 the fully device-resident GFN2 loop filled the orbitals as a closed shell (energy up to
# 49 kcal/mol too low for triplets) and the device gradient of GFN1 and GFN2 built W from integer closed-shell
# occupations (gradient error up to 0.3 Eh/A). Open shells now use the host-driven loop and the host gradient.
# Triplet N2 (-spin 2) and the OH radical (-spin 1), GFN1 and GFN2, GPU against CPU:
#   energy within 1e-7 Eh, every gradient component within 1e-6 Eh/A.
# Skipped (exit 77) when the GPU run reports no usable CUDA device or plugin.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 14: open-shell GFN1/GFN2 on the CUDA path against the CPU"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.0000001"
GRADIENT_TOLERANCE="0.000001"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }

run() {   # run <log> <grad> <xyz> <method> <spin> [args...]
    local LOG=$1; local GRAD=$2; local XYZ=$3; local METHOD=$4; local SPIN=$5; shift 5
    $CURCUMA -sp "$XYZ" -method "$METHOD" -spin "$SPIN" -gradient -dump_gradient "$GRAD" -threads 1 \
        -verbosity 2 -no_bmt "$@" > "$LOG" 2>&1
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

    for CASE in "n2.xyz 2" "oh.xyz 1"; do
        set -- $CASE
        local XYZ=$1; local SPIN=$2
        for METHOD in gfn2 gfn1; do
            local TAG="${XYZ%.xyz}_${METHOD}"
            run "${TAG}_cpu.log" "${TAG}_cpu.grad" "$XYZ" "$METHOD" "$SPIN"
            run "${TAG}_gpu.log" "${TAG}_gpu.grad" "$XYZ" "$METHOD" "$SPIN" -gpu cuda
            if strip "${TAG}_gpu.log" | grep -qiE "no usable CUDA device|requires a CUDA build|plugin"; then
                echo -e "${YELLOW}SKIP${NC}: no usable CUDA device or plugin"
                exit 77
            fi
            assert_string_in_file "EXIT 0" "${TAG}_cpu.log" "$TAG CPU run: exit code 0"
            assert_string_in_file "EXIT 0" "${TAG}_gpu.log" "$TAG GPU run: exit code 0"
            assert_numeric_match "$(energy_of ${TAG}_cpu.log)" "$(energy_of ${TAG}_gpu.log)" "$ENERGY_TOLERANCE" \
                "$TAG: GPU energy equals the CPU energy"
            check_gradient "${TAG}_cpu.grad" "${TAG}_gpu.grad" "$TAG: GPU gradient equals the CPU gradient"
        done
    done

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
