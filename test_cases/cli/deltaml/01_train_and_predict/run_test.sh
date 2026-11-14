#!/bin/bash
# Test: Delta-ML train, predict, reload and optimise
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Delta-ML CLI coverage
#
# Asserts invariants, not accuracy. The fit quality on a 9-atom fixture with 60 randomly
# displaced structures carries no information; what must hold regardless is that the model
# trains and reloads, that the corrected energy really is base plus correction, that an
# optimisation on the corrected surface converges (which an energy/gradient pair that
# disagree would prevent), and that a missing model fails loudly instead of silently
# falling back to the bare base method.
#
# The sharp numerical check of the analytic derivatives lives in the C++ unit test
# test_deltaml, where finite differences can be taken properly.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="deltaml - 01: train, predict, reload, optimise"
TEST_DIR="$SCRIPT_DIR"
INPUT="deltaml_input.xyz"
MODEL="deltaml_input.deltaml.json"
PROBE="probe.xyz"

energy_from_log() {
    sed 's/\x1b\[[0-9;]*m//g' "$1" | grep -oP 'Single Point Energy = \K[-0-9.eE+]+' | head -1
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -f "$MODEL" "$PROBE" *.deltaml.dataset.*

    # --- train --------------------------------------------------------------------
    $CURCUMA -activelearn "$INPUT" -method gfnff -ref_method gfn2 \
        -n_centers 120 -test_fraction 0.25 -threads 1 -no_bmt > out_train.log 2>&1
    assert_exit_code 0 $? "training run exits cleanly"
    assert_file_exists "$MODEL" "model file written"
    assert_string_in_file "curcuma-deltaml-model" "$MODEL" "model carries its format tag"
    assert_string_in_file "held-out set" out_train.log "held-out evaluation reported"

    # The gradient term is what makes a delta model usable for optimisation, so the fit
    # has to improve it over the bare base method even on this crude fixture.
    local G_MODEL G_BASE BETTER
    G_MODEL=$(sed 's/\x1b\[[0-9;]*m//g' out_train.log | grep "gradient RMSE" | tail -1 | awk '{print $6}')
    G_BASE=$(sed 's/\x1b\[[0-9;]*m//g' out_train.log | grep "gradient RMSE" | tail -1 | awk '{print $9}')
    BETTER=$(awk -v a="$G_MODEL" -v b="$G_BASE" 'BEGIN{print (a != "" && b != "" && a < b) ? 1 : 0}')
    assert_equals "1" "$BETTER" "held-out gradient RMSE improved (model=$G_MODEL base=$G_BASE)"

    # --- single point: total must equal base plus correction ------------------------
    head -11 "$INPUT" > "$PROBE"

    $CURCUMA -sp "$PROBE" -method gfnff -no_bmt > out_base.log 2>&1
    assert_exit_code 0 $? "base single point exits cleanly"
    local E_BASE
    E_BASE=$(energy_from_log out_base.log)

    $CURCUMA -sp "$PROBE" -method "deltaml:$MODEL" -verbosity 3 -no_bmt > out_ml.log 2>&1
    assert_exit_code 0 $? "deltaml single point exits cleanly"
    local E_ML CORR
    E_ML=$(energy_from_log out_ml.log)
    CORR=$(sed 's/\x1b\[[0-9;]*m//g' out_ml.log \
        | grep -oP 'deltaml_correction: \K[-0-9.eE+]+' | head -1)

    if [ -z "$E_ML" ] || [ -z "$CORR" ] || [ -z "$E_BASE" ]; then
        assert_equals "energies" "missing (base=$E_BASE ml=$E_ML corr=$CORR)" \
            "could extract base energy, total energy and correction"
    else
        assert_numeric_match "$(awk -v b="$E_BASE" -v c="$CORR" 'BEGIN{printf "%.10f", b+c}')" \
            "$E_ML" "1e-7" "total energy equals base plus ML correction"
    fi

    # --- both spellings of the model reference must agree ---------------------------
    $CURCUMA -sp "$PROBE" -method deltaml -deltaml.model "$MODEL" -no_bmt > out_reload.log 2>&1
    assert_exit_code 0 $? "dotted -deltaml.model form works"
    assert_numeric_match "$E_ML" "$(energy_from_log out_reload.log)" "1e-10" \
        "shorthand and dotted form give the same energy"

    # --- optimisation on the corrected surface --------------------------------------
    # An energy and a gradient that disagree show up here: the line search stalls and the
    # run does not reach convergence.
    $CURCUMA -opt "$PROBE" -method "deltaml:$MODEL" -threads 1 -no_bmt > out_opt.log 2>&1
    assert_exit_code 0 $? "optimisation on the corrected surface exits cleanly"
    assert_file_exists "probe.opt.xyz" "optimised structure written"

    # --- a missing model must fail loudly -------------------------------------------
    $CURCUMA -sp "$PROBE" -method deltaml -no_bmt > out_nomodel.log 2>&1
    assert_string_in_file "requires a model file" out_nomodel.log \
        "missing model reports a clear error"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
