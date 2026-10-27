#!/bin/bash
# Test: active learning is reproducible, validates prospectively and stops by the chosen criterion
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (October 2026)
#
# '-activelearn -al_mode learn' samples with molecular dynamics on the current surface, labels
# a selection and refits, cycle by cycle. Asserted here:
#   (a) two runs with the same -activelearn.seed and -threads 4 write byte-identical learning
#       curves and datasets. Before October 2026 the CSVR thermostat of SimpleMD drew from one
#       process-wide engine, which the four concurrent sampling runs consumed in scheduling
#       order, so the runs diverged already in cycle 0; this assertion fails on a36d9b12.
#   (b) the learning curve carries the prospective-validation columns (val_*) and the pooled
#       ones (pool_*), finite from the second cycle on, where a model exists
#   (c) al_criterion: 'validation' with an unreachable al_target_energy runs all cycles and
#       reports "not converged"; with generous targets it stops before al_cycles and names the
#       criterion; 'both' and 'uncertainty' (indicator threshold 1, always met, patience 1)
#       stop in the first cycle with a model (val_loose/both_loose with al_window 1).
#   (d) pooled validation (al_window 3, Oct 2026): the criterion is tested only on a full pool
#       of 3 cycles x 10 frames, so generous targets with patience 2 stop in cycle 5, not 3.
#
# Input: dimethyl ether (9 atoms, structure library), GFN-FF -> GFN2; about 0.3 s per run.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="deltaml - 05: active learning reproducible, validated, criterion"
TEST_DIR="$SCRIPT_DIR"
INPUT="mol.xyz"
CURVE="mol.deltaml.learning_curve.dat"
DATASET="mol.deltaml.dataset.json"

# learn <dir> [extra args]: one active-learning run in its own directory, log without colours.
learn() {
    local dir=$1
    shift
    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$INPUT" "$dir/"
    (cd "$dir" && $CURCUMA -activelearn "$INPUT" -method gfnff -activelearn.al_mode learn \
        -activelearn.ref_method gfn2 -activelearn.al_cycles 4 -activelearn.al_md_runs 4 \
        -activelearn.al_md_time 500 -activelearn.al_dump_frequency 10 \
        -activelearn.al_candidates 6 -activelearn.al_initial 8 -activelearn.n_centers 100 \
        -activelearn.al_uncertainty 0.05 -activelearn.seed 1 -threads 4 -no_bmt "$@" \
        > raw.log 2>&1)
    local rc=$?
    sed 's/\x1b\[[0-9;]*m//g' "$dir/raw.log" > "$dir/run.log"
    return $rc
}

cycles_run() { grep -c -- "--- cycle" "$1/run.log"; }

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -rf rep1 rep2 val_strict val_loose both_loose unc_loose pool_loose

    # --- (a) reproducibility -------------------------------------------------------
    learn rep1
    assert_exit_code $? 0 "first learn run exits cleanly"
    learn rep2
    assert_exit_code $? 0 "second learn run exits cleanly"
    assert_file_exists "rep1/$CURVE" "learning curve written"
    local SAME_CURVE SAME_DATA
    SAME_CURVE=$(cmp -s "rep1/$CURVE" "rep2/$CURVE" && echo 1 || echo 0)
    SAME_DATA=$(cmp -s "rep1/$DATASET" "rep2/$DATASET" && echo 1 || echo 0)
    assert_equals "1" "$SAME_CURVE" "same seed, -threads 4: learning curves byte-identical"
    assert_equals "1" "$SAME_DATA" "same seed, -threads 4: datasets byte-identical"

    # --- (b) validation columns ------------------------------------------------------
    assert_string_in_file "val_hole_kcal" "rep1/$CURVE" "curve header names the validation columns"
    assert_string_in_file "pool_hole_q_kcal" "rep1/$CURVE" "curve header names the pool columns"
    # Rows of cycle >= 1: 15 columns, val_n > 0, pool_n > 0, columns 8-11 and 13-15 finite numbers.
    local GOOD BAD
    GOOD=$(awk '!/^#/ && $1 >= 1 && NF == 15 && $7 > 0 && $12 > 0 {ok=1; for (i=8;i<=15;i++) if (i != 12 && $i !~ /^-?[0-9.]+(e[-+]?[0-9]+)?$/) ok=0; if (ok) n++} END{print n+0}' "rep1/$CURVE")
    BAD=$(awk '!/^#/ && $1 >= 1 && !(NF == 15 && $7 > 0 && $12 > 0)' "rep1/$CURVE" | wc -l)
    assert_equals "1" "$(awk -v g="$GOOD" -v b="$BAD" 'BEGIN{print (g >= 2 && b == 0) ? 1 : 0}')" \
        "validation columns finite from cycle 2 on ($GOOD rows, $BAD bad)"
    assert_string_in_file "validation on 10 new frames" "rep1/run.log" "validation line in the log"

    # --- (c) criteria -----------------------------------------------------------------
    learn val_strict -activelearn.al_criterion validation -activelearn.al_target_energy 1e-9
    assert_exit_code $? 0 "validation criterion, unreachable target: exits cleanly"
    assert_equals "4" "$(cycles_run val_strict)" "unreachable target runs all 4 cycles"
    assert_string_in_file "not converged -- criterion 'validation'" "val_strict/run.log" \
        "unreachable target reports not converged"

    learn val_loose -activelearn.al_criterion validation -activelearn.al_target_energy 1000 \
        -activelearn.al_target_hole 1000 -activelearn.al_window 1
    assert_exit_code $? 0 "validation criterion, generous targets: exits cleanly"
    assert_equals "3" "$(cycles_run val_loose)" "generous targets stop after cycle 3 (patience 2)"
    assert_string_in_file "converged in cycle 3 by criterion validation" "val_loose/run.log" \
        "convergence names the validation criterion"
    assert_string_not_in_file "not converged" "val_loose/run.log" "no not-converged line"

    # Indicator threshold 1 is always met but also leaves nothing to select (picks need an
    # indicator above it), so these two runs use patience 1: they stop in the first cycle
    # with a model, before the empty selection would end the run.
    learn both_loose -activelearn.al_criterion both -activelearn.al_uncertainty 1.0 \
        -activelearn.al_target_energy 1000 -activelearn.al_target_hole 1000 -activelearn.al_patience 1 \
        -activelearn.al_window 1
    assert_exit_code $? 0 "both criteria, generous targets: exits cleanly"
    assert_string_in_file "converged in cycle 2 by criterion both" "both_loose/run.log" \
        "convergence names both criteria"

    learn unc_loose -activelearn.al_criterion uncertainty -activelearn.al_uncertainty 1.0 \
        -activelearn.al_patience 1 -activelearn.al_validation_frames 0
    assert_exit_code $? 0 "uncertainty criterion, threshold 1: exits cleanly"
    assert_string_in_file "converged in cycle 2 by criterion uncertainty" "unc_loose/run.log" \
        "uncertainty criterion stops at the first cycle with a model (patience 1)"

    # --- (d) pooled validation ----------------------------------------------------------
    learn pool_loose -activelearn.al_criterion validation -activelearn.al_target_energy 1000 \
        -activelearn.al_target_hole 1000 -activelearn.al_cycles 6
    assert_exit_code $? 0 "pooled validation, generous targets: exits cleanly"
    assert_string_in_file "window of 3 not yet full" "pool_loose/run.log" "pool not full in cycle 2"
    assert_string_in_file "validation pool 30 frames from 3 cycle(s) (model" "pool_loose/run.log" \
        "full pool of 3 cycles x 10 frames"
    assert_string_in_file "converged in cycle 5 by criterion validation (pool of 30 frames" \
        "pool_loose/run.log" "full pool needed: patience 2 stops in cycle 5"

    rm -rf *.snapshots
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
