#!/bin/bash
# Test: Delta-ML energy-only second model
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (October 2026) - Delta-ML CLI coverage, energy model
#
# '-activelearn' in train mode writes a second, energy-only model next to the main one
# (<basename>.deltaml.energy.json, "role": "energy"). Asserted here, without golden numbers:
#   - both files are written, the energy file carries role "energy", the main file "joint"
#   - the log reports both held-out blocks
#   - '-activelearn.energy_model false' writes no energy file and leaves the main model
#     byte-identical (the second fit must not touch the first)
#   - '-sp -method deltaml:<energy model>' runs; the forces warning appears only when a
#     gradient is requested (-dump_gradient, -opt), not for a plain energy
#   - the analytic gradient of the assembled energy-model surface matches central
#     differences of '-sp' energies (finite-difference gate, root CLAUDE.md): displacement
#     1e-3 Angstrom on a handful of coordinates, tolerance 1e-5 Eh/Angstrom
#
# Training data: the 9-atom fixture of test 01 (copied by add_cli_test into the sibling
# build directory), so no structure file is duplicated here.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="deltaml - 03: energy-only second model"
TEST_DIR="$SCRIPT_DIR"
INPUT="em_input.xyz"
MODEL="em_input.deltaml.json"
EMODEL="em_input.deltaml.energy.json"
PROBE="probe.xyz"
WARN_PATTERN="energy-only model"

strip_ansi() { sed 's/\x1b\[[0-9;]*m//g' "$1"; }

# -sp energy at verbosity 3 as base + correction (10 decimals each); the plain
# "Single Point Energy" line has only 8, too coarse for a 1e-3 Angstrom difference.
sp_energy() {
    local xyz=$1 log=$2
    $CURCUMA -sp "$xyz" -method "deltaml:$EMODEL" -verbosity 3 -no_bmt > "$log" 2>&1 || return 1
    local b c
    b=$(strip_ansi "$log" | grep -oP 'deltaml_base_energy: \K[-0-9.eE+]+' | head -1)
    c=$(strip_ansi "$log" | grep -oP 'deltaml_correction: \K[-0-9.eE+]+' | head -1)
    [ -n "$b" ] && [ -n "$c" ] || return 1
    awk -v b="$b" -v c="$c" 'BEGIN{printf "%.12f", b + c}'
}

# Write $1 with coordinate (atom $3, component $4 in 0..2) shifted by $5 Angstrom to $2.
displace() {
    awk -v a="$3" -v k="$4" -v h="$5" 'NR <= 2 {print; next}
        { if (NR - 2 == a + 1) $(k + 2) = sprintf("%.10f", $(k + 2) + h);
          printf "%-3s %18.10f %18.10f %18.10f\n", $1, $2, $3, $4 }' "$1" > "$2"
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -rf *.deltaml.* "$PROBE" probe_*.xyz grad.txt noem *.topo.json

    cp "$SCRIPT_DIR/../01_train_and_predict/deltaml_input.xyz" "$INPUT"
    assert_file_exists "$INPUT" "training fixture copied from test 01"

    # --- train: both models -------------------------------------------------------
    $CURCUMA -activelearn "$INPUT" -method gfnff -ref_method gfn2 \
        -n_centers 120 -test_fraction 0.25 -threads 1 -no_bmt > out_train.log 2>&1
    assert_exit_code $? 0 "training run exits cleanly"
    assert_file_exists "$MODEL" "main model written"
    assert_file_exists "$EMODEL" "energy model written next to it"
    assert_string_in_file '"role": "energy"' "$EMODEL" "energy model carries role energy"
    assert_string_in_file '"role": "joint"' "$MODEL" "main model carries role joint"
    assert_string_in_file "held-out set (energy model" out_train.log "energy model held-out block"
    assert_string_in_file "held-out set (joint model" out_train.log "joint model held-out block"

    # --- energy_model false: no second file, main model unchanged ------------------
    mkdir -p noem
    cp "$INPUT" noem/
    (cd noem && $CURCUMA -activelearn "$INPUT" -method gfnff -ref_method gfn2 \
        -n_centers 120 -test_fraction 0.25 -threads 1 -activelearn.energy_model false \
        -no_bmt > out_train.log 2>&1)
    assert_exit_code $? 0 "training without energy model exits cleanly"
    assert_file_not_exists "noem/$EMODEL" "-activelearn.energy_model false writes no energy model"
    if cmp -s "$MODEL" "noem/$MODEL"; then
        assert_equals "identical" "identical" "main model byte-identical with and without energy model"
    else
        assert_equals "identical" "different" "main model byte-identical with and without energy model"
    fi

    # --- single point and the warning ---------------------------------------------
    head -11 "$INPUT" > "$PROBE"
    $CURCUMA -sp "$PROBE" -method "deltaml:$EMODEL" -no_bmt > out_sp.log 2>&1
    assert_exit_code $? 0 "energy-model single point exits cleanly"
    assert_string_in_file "Single Point Energy" out_sp.log "energy reported"
    assert_string_not_in_file "$WARN_PATTERN" out_sp.log "no forces warning without gradient"

    $CURCUMA -sp "$PROBE" -method "deltaml:$EMODEL" -dump_gradient grad.txt -no_bmt > out_grad.log 2>&1
    assert_exit_code $? 0 "energy-model gradient single point exits cleanly"
    assert_string_in_file "$WARN_PATTERN" out_grad.log "forces warning with -dump_gradient"
    local NWARN
    NWARN=$(strip_ansi out_grad.log | grep -c "$WARN_PATTERN")
    assert_equals "1" "$NWARN" "warning printed once"

    $CURCUMA -sp "$PROBE" -method "deltaml:$MODEL" -dump_gradient grad_main.txt -no_bmt > out_grad_main.log 2>&1
    assert_string_not_in_file "$WARN_PATTERN" out_grad_main.log "no warning for the joint model"

    $CURCUMA -opt "$PROBE" -method "deltaml:$EMODEL" -opt.optimizer lbfgs -threads 1 -no_bmt \
        > out_opt.log 2>&1
    assert_string_in_file "$WARN_PATTERN" out_opt.log "forces warning with -opt"

    # --- finite-difference gate on the assembled surface ----------------------------
    # grad.txt: two comment lines, then one "gx gy gz" row per atom [Eh/Angstrom].
    local H=0.001 TOL=1e-5 MAXDEV=0 FAILED=0
    for spec in "0 0" "1 1" "3 2" "5 0" "8 1"; do
        set -- $spec
        local atom=$1 comp=$2
        displace "$PROBE" probe_p.xyz "$atom" "$comp" "$H"
        displace "$PROBE" probe_m.xyz "$atom" "$comp" "-$H"
        local EP EM GA FD DEV
        EP=$(sp_energy probe_p.xyz out_fd_p.log) || FAILED=1
        EM=$(sp_energy probe_m.xyz out_fd_m.log) || FAILED=1
        GA=$(grep -v '^#' grad.txt | awk -v r="$((atom + 1))" -v c="$((comp + 1))" 'NR == r {print $c}')
        FD=$(awk -v p="$EP" -v m="$EM" -v h="$H" 'BEGIN{printf "%.10e", (p - m) / (2 * h)}')
        DEV=$(awk -v a="$GA" -v f="$FD" 'BEGIN{d = a - f; if (d < 0) d = -d; printf "%.3e", d}')
        echo "  atom $atom comp $comp: analytic $GA  FD $FD  |diff| $DEV"
        MAXDEV=$(awk -v a="$MAXDEV" -v b="$DEV" 'BEGIN{print (b > a) ? b : a}')
    done
    assert_equals "0" "$FAILED" "all displaced single points produced an energy"
    local OK
    OK=$(awk -v d="$MAXDEV" -v t="$TOL" 'BEGIN{print (d < t) ? 1 : 0}')
    assert_equals "1" "$OK" "energy-model gradient matches central differences (max |diff| $MAXDEV < $TOL Eh/A)"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
