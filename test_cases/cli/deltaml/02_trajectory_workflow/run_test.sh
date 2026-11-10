#!/bin/bash
# Test: Delta-ML trajectory workflow -- sample, fit, run dynamics on the corrected surface
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Delta-ML CLI coverage, trajectory use case
#
# The conformer test (01) covers a model fitted to a given ensemble. This one covers the
# other use case, where the structures come from dynamics and the model is judged on forces
# rather than on energies: generate frames with a short run, fit to them, then integrate on
# the corrected surface.
#
# Deliberately no golden numbers. A 9-atom fixture says nothing about accuracy, so the test
# asserts the properties that must hold for any fit: the gradient error improves on the
# held-out split, and dynamics on the resulting surface neither blows up nor heats away.
# Judging a trajectory model by how far it moves an optimised minimum would be wrong -- a
# lower RMS force error does not constrain where the gradient vanishes.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="deltaml - 02: trajectory workflow"
TEST_DIR="$SCRIPT_DIR"
SEED_XYZ="seed.xyz"
FRAMES="frames.xyz"
MODEL="frames.deltaml.json"

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -rf "$MODEL" "$FRAMES" *.snapshots *.topo.json *.deltaml.*

    # --- generate the trajectory -------------------------------------------------
    # 300 K rather than a hot run: this should look like the dynamics the model will later
    # be used for, and a 9-atom molecule at 500 K samples a far rougher difference than any
    # realistic trajectory would.
    $CURCUMA -md "$SEED_XYZ" -method gfnff -max_time 1500 -time_step 1 -temperature 300 \
        -dump_frequency 5 -seed 4242 -print_frequency 100000 -no_bmt > out_md.log 2>&1
    assert_exit_code 0 $? "sampling run exits cleanly"

    # Legacy output mode puts the trajectory in <basename>.snapshots/.
    local TRJ
    TRJ=$(find . -name "*.trj.xyz" | head -1)
    if [ -z "$TRJ" ]; then
        assert_equals "a trajectory" "none" "sampling produced a trajectory"
        print_test_summary
        exit 1
    fi
    cp "$TRJ" "$FRAMES"
    local NFRAMES
    NFRAMES=$(count_xyz_structures "$FRAMES")
    assert_equals "1" "$(awk -v n="$NFRAMES" 'BEGIN{print (n>=10)?1:0}')" \
        "trajectory has enough frames to fit ($NFRAMES)"

    # --- fit to the frames --------------------------------------------------------
    $CURCUMA -activelearn "$FRAMES" -method gfnff -ref_method gfn2 \
        -n_centers 200 -test_fraction 0.25 -threads 1 -no_bmt > out_fit.log 2>&1
    assert_exit_code 0 $? "fit to trajectory frames exits cleanly"
    assert_file_exists "$MODEL" "model written"

    # Forces are the acceptance criterion for this use case.
    local G_MODEL G_BASE BETTER
    G_MODEL=$(sed 's/\x1b\[[0-9;]*m//g' out_fit.log | grep "gradient RMSE" | tail -1 | awk '{print $6}')
    G_BASE=$(sed 's/\x1b\[[0-9;]*m//g' out_fit.log | grep "gradient RMSE" | tail -1 | awk '{print $9}')
    BETTER=$(awk -v a="$G_MODEL" -v b="$G_BASE" 'BEGIN{print (a != "" && b != "" && a < b) ? 1 : 0}')
    assert_equals "1" "$BETTER" "held-out gradient RMSE improved (model=$G_MODEL base=$G_BASE)"

    # --- integrate on the corrected surface ---------------------------------------
    # Constant energy, so nothing external keeps the run in line: a surface whose gradient
    # does not match its own energy pumps the temperature away within a few hundred steps.
    $CURCUMA -md "$SEED_XYZ" -method "deltaml:$MODEL" -max_time 400 -time_step 0.5 \
        -temperature 300 -thermostat none -seed 99 -print_frequency 40 -write_xyz false \
        -no_bmt > out_nve.log 2>&1
    assert_exit_code 0 $? "dynamics on the corrected surface exits cleanly"

    local TMEAN
    TMEAN=$(sed 's/\x1b\[[0-9;]*m//g' out_nve.log \
        | awk 'NF>=9 && $1 ~ /^[0-9]+\.[0-9]+$/ {sum+=$8; n++} END{if(n>0) printf "%.1f", sum/n}')
    if [ -z "$TMEAN" ]; then
        assert_equals "a temperature series" "none" "dynamics produced a temperature series"
    else
        # Equipartition takes the 300 K initialisation to roughly half that; a runaway
        # surface leaves this band immediately.
        assert_equals "1" "$(awk -v t="$TMEAN" 'BEGIN{print (t>20 && t<900)?1:0}')" \
            "mean temperature stayed in a physical band (${TMEAN} K)"
    fi
    # Only the numeric columns matter: the status line carries a cosmetic "-nan" for the
    # remaining-time estimate on the very first step, before any rate can be known.
    local BADROWS
    BADROWS=$(sed 's/\x1b\[[0-9;]*m//g' out_nve.log \
        | awk 'NF>=9 && $1 ~ /^[0-9]+\.[0-9]+$/ && ($2 ~ /nan/ || $6 ~ /nan/ || $8 ~ /nan/)' | wc -l)
    assert_equals "0" "$BADROWS" "energy and temperature stayed finite"

    rm -rf *.snapshots
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
