#!/bin/bash
# Test: Delta-ML labels equal fresh single points, in input order
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (October 2026) - regression test for the label/deployment base-energy mismatch
#
# '-activelearn' stores, per structure, energy_base, gradient_base, energy_ref and gradient_ref.
# The deployed model adds its correction to a freshly set-up base method, so the stored base
# labels must be what '-sp' returns for the same geometry. Before Oct 2026 every thread reused
# one calculator for its whole share of the structures, and GFN-FF re-initialised on a new
# geometry kept distances of the previous one in its EEQ solver: the stored base energy was
# exact only for the first structure of each thread (docs/DELTA_ML.md, "Open: the anchor
# freezes more than it should"). Asserted here:
#   - labelling with -threads 4 writes the dataset in input order
#   - for every structure (including all that are not the first of a thread chunk):
#     |fresh - stored| <= 1e-7 Eh for energy_base and energy_ref and <= 1e-7 Eh/A for
#     gradient_base and gradient_ref, fresh = '-sp -dump_gradient -gfnff.cache_topology false'
#   - -threads 1 and -threads 4 write byte-identical datasets
#
# Input: helicene-frames17 (17 frames, 42 atoms, structure library); the frames differ in
# geometry, which is what made the old reuse visible. On bc7742e6 this test fails: dataset not
# in input order with -threads 4; with -threads 1 (order intact) energy_base is off by up to
# 1.4e-4 Eh and gradient_base by 3.9e-4 Eh/A, the reference labels are exact.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="deltaml - 04: labels equal fresh single points, input order"
TEST_DIR="$SCRIPT_DIR"
INPUT="label_set.xyz"
BASE_METHOD="gfnff"
REF_METHOD="gfn2"
TOL_E="1e-7"
TOL_G="1e-7"

label() {
    local threads=$1 dir=$2
    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$INPUT" "$dir/"
    (cd "$dir" && $CURCUMA -activelearn "$INPUT" -method "$BASE_METHOD" \
        -activelearn.ref_method "$REF_METHOD" -activelearn.n_centers 20 \
        -activelearn.energy_model false -threads "$threads" -no_bmt > label.log 2>&1)
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    cleanup_test_artifacts
    rm -rf run_t4 run_t1 fresh

    label 4 run_t4
    assert_exit_code $? 0 "labelling with -threads 4 exits cleanly"
    local DS4="run_t4/label_set.deltaml.dataset.json"
    assert_file_exists "$DS4" "dataset written (-threads 4)"

    # Fresh single points, one directory per frame, so no cache file is shared.
    local NFRAMES
    NFRAMES=$(count_xyz_structures "$INPUT")
    local NAT
    NAT=$(head -1 "$INPUT" | tr -d '[:space:]')
    local SP_FAILED=0
    for ((k = 0; k < NFRAMES; k++)); do
        local d
        d=$(printf "fresh/f%03d" "$k")
        mkdir -p "$d"
        sed -n "$((k * (NAT + 2) + 1)),$(((k + 1) * (NAT + 2)))p" "$INPUT" > "$d/x.xyz"
        (cd "$d" && $CURCUMA -sp x.xyz -method "$BASE_METHOD" -gfnff.cache_topology false \
            -dump_gradient base.grad -threads 1 -no_bmt > base.log 2>&1) || SP_FAILED=1
        (cd "$d" && $CURCUMA -sp x.xyz -method "$REF_METHOD" -gfnff.cache_topology false \
            -dump_gradient ref.grad -threads 1 -no_bmt > ref.log 2>&1) || SP_FAILED=1
    done
    assert_equals "0" "$SP_FAILED" "all $NFRAMES fresh single points ran"

    python3 -I "$SCRIPT_DIR/compare_labels.py" "$DS4" "$INPUT" fresh "$TOL_E" "$TOL_G" > compare.log 2>&1
    local CMP=$?
    cat compare.log
    assert_equals "0" "$CMP" "stored labels equal fresh single points and follow the input order"

    label 1 run_t1
    assert_exit_code $? 0 "labelling with -threads 1 exits cleanly"
    cmp -s "$DS4" "run_t1/label_set.deltaml.dataset.json"
    assert_equals "0" "$?" "-threads 1 and -threads 4 write identical datasets"

    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
