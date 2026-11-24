#!/bin/bash
# Test: screened (sparse) integral storage on ROCm, -gpu_sparse_integrals on (docs/ROCM_PARITY_A7.md)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 24, 2026). AI-generated, machine-tested only.
#
# pep20 peptide (309 atoms): the distance screen keeps 73 % (GFN2) / 85 % (GFN1) of the AO pairs, so pairs are dropped.
# (a) -gpu rocm without the switch: dense storage, no storage line
# (b) GFN2 -gpu_sparse_integrals on -gpu_strict true: screened storage on the device, energy and gradient equal the
#     dense run (1e-8 Eh, 1e-8 Eh/A; measured 6e-14 Eh, 6e-15 Eh/A)
# (c) the same for GFN1
# (d) a second screened run gives the identical gradient file (md5)
# Skipped (exit 77) when the run in (a) reports no usable ROCm device or plugin.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="sqm - 24: screened integral storage on ROCm"
TEST_DIR="$SCRIPT_DIR"
ENERGY_TOLERANCE="0.00000001"

strip()     { sed 's/\x1b\[[0-9;]*m//g' "$1"; }
energy_of() { strip "$1" | grep -oP 'Single Point Energy = \K[-0-9.e]+' | head -1; }

run() {   # run <log> <method> <grad> [args...]
    local LOG=$1; local M=$2; local G=$3; shift 3
    $CURCUMA -sp pep20.xyz -method "$M" -gpu rocm -threads 4 -verbosity 2 -no_bmt -scf_threshold 1e-9 \
        -dump_gradient "$G" "$@" > "$LOG" 2>&1
    echo "EXIT $?" >> "$LOG"
}

check_grad() {   # check_grad <a> <b> <label>: max |a - b| over the gradient rows <= 1e-8 Eh/A
    TESTS_RUN=$((TESTS_RUN + 1))
    local d
    d=$(paste <(grep -v '^#' "$1") <(grep -v '^#' "$2") | awk '{for(i=1;i<=3;i++){x=$i-$(i+3); if(x<0)x=-x; if(x>m)m=x}} END{printf "%.3e", m+0}')
    if [ -s "$1" ] && [ -s "$2" ] && awk -v d="$d" 'BEGIN{exit !(d <= 1e-8)}'; then
        echo -e "${GREEN}PASS${NC}: $3 (max |dG| $d Eh/A)"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $3 (max |dG| '$d' Eh/A)"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts

    # (a) dense default; decides about the skip
    run dense2.log gfn2 dense2.grad
    if ! strip dense2.log | grep -q "fully on the device"; then
        if strip dense2.log | grep -qiE "no usable|GPU fallbacks during this run|requires a ROCm build|plugin"; then
            echo -e "${YELLOW}SKIP${NC}: no usable ROCm device or plugin"
            exit 77
        fi
    fi
    assert_string_in_file "EXIT 0" dense2.log "GFN2 dense GPU run: exit code 0"
    assert_string_not_in_file "GPU integrals:" dense2.log "GFN2 dense GPU run: no storage line (ROCm default dense)"

    # (b) GFN2 screened
    run sp2.log gfn2 sp2.grad -gpu_sparse_integrals on -gpu_strict true
    assert_string_in_file "EXIT 0" sp2.log "GFN2 screened run with -gpu_strict true: exit code 0"
    assert_string_in_file "GPU integrals: screened pair storage" sp2.log "GFN2: screened storage in use"
    assert_numeric_match "$(energy_of dense2.log)" "$(energy_of sp2.log)" "$ENERGY_TOLERANCE" "GFN2: screened energy equals dense"
    check_grad sp2.grad dense2.grad "GFN2: screened gradient equals dense"

    # (c) GFN1
    run dense1.log gfn1 dense1.grad
    run sp1.log gfn1 sp1.grad -gpu_sparse_integrals on -gpu_strict true
    assert_string_in_file "EXIT 0" sp1.log "GFN1 screened run with -gpu_strict true: exit code 0"
    assert_string_in_file "GPU integrals: screened pair storage" sp1.log "GFN1: screened storage in use"
    assert_numeric_match "$(energy_of dense1.log)" "$(energy_of sp1.log)" "$ENERGY_TOLERANCE" "GFN1: screened energy equals dense"
    check_grad sp1.grad dense1.grad "GFN1: screened gradient equals dense"

    # (d) run-to-run reproducibility
    run sp2b.log gfn2 sp2b.grad -gpu_sparse_integrals on
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s sp2.grad ] && [ "$(md5sum < sp2.grad)" = "$(md5sum < sp2b.grad)" ]; then
        echo -e "${GREEN}PASS${NC}: GFN2 screened: second run gives the identical gradient file"; TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: GFN2 screened: gradient files of two runs differ (or are missing)"; TESTS_FAILED=$((TESTS_FAILED + 1))
    fi

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
