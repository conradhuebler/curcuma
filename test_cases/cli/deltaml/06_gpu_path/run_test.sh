#!/bin/bash
# Test: Delta-ML on the GPU (descriptor, kernel, forces and fit on the device)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (October 2026) - Delta-ML GPU path, docs/DELTA_ML_GPU_2026-10.md
#
# Skips (exit 77) when the build has no CUDA plugin or the machine no device: the model then
# reports "unavailable" and stays on the CPU, which is itself asserted to be announced.
# Asserted with a device, without golden numbers:
#   - prediction: -gpu cuda equals the CPU to 1e-9 Eh and 1e-9 Eh/Angstrom on three geometries.
#     The base is overridden to UFF (no GPU path), so the difference is the correction alone
#   - FD gate of the fully assembled surface with -gpu cuda (GFN-FF base on the GPU as well):
#     -dump_gradient vs central differences of -sp energies, h 1e-3 Angstrom, tol 1e-5 Eh/A
#   - fit: a model fitted with -gpu cuda predicts what the CPU-fitted model predicts
#     (1e-7 Eh, 5e-6 Eh/A; the summation order of the normal equations, see below)
#
# Training data: the 9-atom fixture of test 01 (copied by add_cli_test into the sibling
# build directory), so no structure file is duplicated here.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="deltaml - 06: GPU path"
TEST_DIR="$SCRIPT_DIR"
INPUT="gp_input.xyz"
MODEL="gp_input.deltaml.json"

# Energy (15 decimals) and gradient rows of a -dump_gradient file.
dump_energy() { grep -oP '^# energy \K-?[0-9.]+' "$1"; }

# max |a - b| over all gradient components of two dump files
max_grad_diff() {
    [ -s "$1" ] && [ -s "$2" ] || { echo nan; return; }
    paste <(grep -v '^#' "$1") <(grep -v '^#' "$2") | awk '{for (k = 1; k <= 3; ++k) {d = $k - $(k + 3);
        if (d < 0) d = -d; if (d > m) m = d}} END {printf "%.3e", m + 0}'
}

# A missing number (failed run, absent file) yields "nan", which every "below" check rejects.
abs_diff() {
    [ -n "$1" ] && [ -n "$2" ] || { echo nan; return; }
    awk -v a="$1" -v b="$2" 'BEGIN{d = a - b; if (d < 0) d = -d; printf "%.3e", d}'
}
below() { [ "$1" != nan ] && awk -v d="$1" -v t="$2" 'BEGIN{print (d < t) ? 1 : 0}' || echo 0; }

# Write $1 with coordinate (atom $3, component $4 in 0..2) shifted by $5 Angstrom to $2.
displace() {
    awk -v a="$3" -v k="$4" -v h="$5" 'NR <= 2 {print; next}
        { if (NR - 2 == a + 1) $(k + 2) = sprintf("%.10f", $(k + 2) + h);
          printf "%-3s %18.10f %18.10f %18.10f\n", $1, $2, $3, $4 }' "$1" > "$2"
}

# frame $2 (0-based) of the multi-structure file $1 -> $3
frame() {
    awk -v f="$2" 'NR == 1 {n = $1} {k = int((NR - 1) / (n + 2)); if (k == f) print}' "$1" > "$3"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -rf *.deltaml.* probe*.xyz *.grad gpufit *.topo.json out_*.log

    cp "$SCRIPT_DIR/../01_train_and_predict/deltaml_input.xyz" "$INPUT"
    assert_file_exists "$INPUT" "training fixture copied from test 01"

    $CURCUMA -activelearn "$INPUT" -method gfnff -ref_method gfn2 -n_centers 120 -test_fraction 0.25 \
        -activelearn.energy_model false -threads 1 -no_bmt > out_train.log 2>&1
    assert_exit_code $? 0 "CPU fit exits cleanly"
    assert_file_exists "$MODEL" "CPU model written"

    frame "$INPUT" 0 probe0.xyz
    frame "$INPUT" 3 probe1.xyz
    frame "$INPUT" 7 probe2.xyz

    # --- skip without a usable device -----------------------------------------------
    $CURCUMA -sp probe0.xyz -method "deltaml:$MODEL" -deltaml.base uff -gpu cuda -no_bmt > out_probe.log 2>&1
    if grep -q "unavailable for prediction" out_probe.log; then
        echo "SKIP: no CUDA plugin or device -- the CPU fallback was announced:"
        grep "unavailable for prediction" out_probe.log | sed 's/\x1b\[[0-9;]*m//g'
        exit 77
    fi

    # --- prediction parity, base UFF on both sides ------------------------------------
    local k WORST_E=0 WORST_G=0
    for k in 0 1 2; do
        $CURCUMA -sp probe$k.xyz -method "deltaml:$MODEL" -deltaml.base uff -dump_gradient c$k.grad \
            -no_bmt > out_c$k.log 2>&1
        $CURCUMA -sp probe$k.xyz -method "deltaml:$MODEL" -deltaml.base uff -dump_gradient g$k.grad \
            -gpu cuda -no_bmt > out_g$k.log 2>&1
        local DE DG
        DE=$(abs_diff "$(dump_energy c$k.grad)" "$(dump_energy g$k.grad)")
        DG=$(max_grad_diff c$k.grad g$k.grad)
        echo "  geometry $k: |dE| $DE Eh, max |dg| $DG Eh/A"
        [ "$DE" = nan ] && WORST_E=nan
        [ "$DG" = nan ] && WORST_G=nan
        [ "$WORST_E" != nan ] && WORST_E=$(awk -v a="$WORST_E" -v b="$DE" 'BEGIN{print (b > a) ? b : a}')
        [ "$WORST_G" != nan ] && WORST_G=$(awk -v a="$WORST_G" -v b="$DG" 'BEGIN{print (b > a) ? b : a}')
    done
    assert_equals "1" "$(below "$WORST_E" 1e-9)" "GPU energy equals CPU (max |dE| $WORST_E < 1e-9 Eh)"
    assert_equals "1" "$(below "$WORST_G" 1e-9)" "GPU gradient equals CPU (max |dg| $WORST_G < 1e-9 Eh/A)"

    # --- FD gate of the assembled surface on the GPU -----------------------------------
    $CURCUMA -sp probe0.xyz -method "deltaml:$MODEL" -gpu cuda -dump_gradient fd0.grad -no_bmt \
        -gfnff.cache_topology false > out_fd0.log 2>&1
    assert_exit_code $? 0 "GPU gradient single point exits cleanly"
    local H=0.001 TOL=1e-5 MAXDEV=0 FAILED=0 spec
    for spec in "0 0" "1 1" "3 2" "5 0" "8 1"; do
        set -- $spec
        local atom=$1 comp=$2 EP EM GA FD DEV
        displace probe0.xyz probe_p.xyz "$atom" "$comp" "$H"
        displace probe0.xyz probe_m.xyz "$atom" "$comp" "-$H"
        $CURCUMA -sp probe_p.xyz -method "deltaml:$MODEL" -gpu cuda -dump_gradient p.grad -no_bmt \
            -gfnff.cache_topology false > out_fd_p.log 2>&1 || FAILED=1
        $CURCUMA -sp probe_m.xyz -method "deltaml:$MODEL" -gpu cuda -dump_gradient m.grad -no_bmt \
            -gfnff.cache_topology false > out_fd_m.log 2>&1 || FAILED=1
        EP=$(dump_energy p.grad 2>/dev/null)
        EM=$(dump_energy m.grad 2>/dev/null)
        rm -f p.grad m.grad
        GA=$(grep -v '^#' fd0.grad | awk -v r="$((atom + 1))" -v c="$((comp + 1))" 'NR == r {print $c}')
        FD=$(awk -v p="$EP" -v m="$EM" -v h="$H" 'BEGIN{printf "%.10e", (p - m) / (2 * h)}')
        [ -n "$EP" ] && [ -n "$EM" ] || FD=""
        DEV=$(abs_diff "$GA" "$FD")
        echo "  atom $atom comp $comp: analytic $GA  FD $FD  |diff| $DEV"
        [ "$DEV" = nan ] && MAXDEV=nan
        [ "$MAXDEV" != nan ] && MAXDEV=$(awk -v a="$MAXDEV" -v b="$DEV" 'BEGIN{print (b > a) ? b : a}')
    done
    assert_equals "0" "$FAILED" "all displaced GPU single points produced an energy"
    assert_equals "1" "$(below "$MAXDEV" "$TOL")" "GPU gradient matches central differences (max |diff| $MAXDEV < $TOL Eh/A)"

    # --- fit on the GPU -------------------------------------------------------------
    mkdir -p gpufit
    cp "$INPUT" gpufit/
    (cd gpufit && $CURCUMA -activelearn "$INPUT" -method gfnff -ref_method gfn2 -n_centers 120 \
        -test_fraction 0.25 -activelearn.energy_model false -threads 1 -gpu cuda -no_bmt > out_train.log 2>&1)
    assert_exit_code $? 0 "GPU fit exits cleanly"
    assert_string_not_in_file "unavailable for the fit" gpufit/out_train.log "fit ran on the device"
    assert_string_not_in_file "refitting on the CPU" gpufit/out_train.log "no device failure during the fit"
    $CURCUMA -sp probe1.xyz -method "deltaml:gpufit/$MODEL" -deltaml.base uff -dump_gradient f_gpu.grad \
        -no_bmt > out_fg.log 2>&1
    local DE DG
    DE=$(abs_diff "$(dump_energy c1.grad)" "$(dump_energy f_gpu.grad)")
    DG=$(max_grad_diff c1.grad f_gpu.grad)
    echo "  GPU-fitted vs CPU-fitted model (CPU evaluation): |dE| $DE Eh, max |dg| $DG Eh/A"
    assert_equals "1" "$(below "$DE" 1e-7)" "GPU fit predicts the CPU-fit energy (|dE| $DE < 1e-7 Eh)"
    # 5e-6, not rounding level: cuBLAS sums the normal equations in another order and the
    # system (M = 120, lambda 1e-6) amplifies that. Measured Oct 2026: the July GPU fit, which
    # built the design block on the host, differs from the CPU fit by the same 2.5e-8 Eh /
    # 7.7e-7 Eh/A as the device-built one.
    assert_equals "1" "$(below "$DG" 5e-6)" "GPU fit predicts the CPU-fit gradient (max |dg| $DG < 5e-6 Eh/A)"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
