#!/bin/bash
# Test: -md.seed seeds the initial velocities
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (October 2026) - regression test for the ineffective -md.seed
#
# InitVelocities() used a static default-constructed engine, so six MD runs with the seeds 21 to
# 26 from one start structure were bit-identical. Now a positive -md.seed also seeds the initial
# velocities, while without it the velocities keep their fixed legacy sequence (golden values of
# other tests). Asserted on a short NVE run (GFN-FF, dimethyl ether), compared by the md5 of the
# MD diagnostics file:
#   - the same seed twice gives identical diagnostics
#   - two different seeds give different diagnostics
#   - two runs without a seed give identical diagnostics (legacy behaviour)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 16: -md.seed seeds the initial velocities"
TEST_DIR="$SCRIPT_DIR"

md_digest() {
    local dir=$1; shift
    rm -rf "$dir"; mkdir -p "$dir"; cp input.xyz "$dir/"
    (cd "$dir" && $CURCUMA -md input.xyz -method gfnff -md.temperature 300 -md.max_time 100 \
        -md.time_step 1 -md.thermostat none -md_diagnostics true -threads 1 -no_bmt "$@" \
        > md.log 2>&1; echo $? > rc.txt)
    cat "$dir/rc.txt"
    md5sum "$dir/input.diag.jsonl" 2>/dev/null | cut -c1-32
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    local out_a out_b out_c out_d out_e
    out_a=$(md_digest run_s21a -md.seed 21)
    out_b=$(md_digest run_s21b -md.seed 21)
    out_c=$(md_digest run_s22 -md.seed 22)
    out_d=$(md_digest run_n1)
    out_e=$(md_digest run_n2)

    assert_equals "0" "$(echo "$out_a" | head -1)" "MD with -md.seed 21 exits cleanly"
    local da db dc dd de
    da=$(echo "$out_a" | tail -1); db=$(echo "$out_b" | tail -1); dc=$(echo "$out_c" | tail -1)
    dd=$(echo "$out_d" | tail -1); de=$(echo "$out_e" | tail -1)

    assert_equals "yes" "$([ -n "$da" ] && echo yes || echo no)" "diagnostics file written"
    assert_equals "$da" "$db" "same seed twice gives identical diagnostics"
    assert_equals "different" "$([ "$da" != "$dc" ] && echo different || echo equal)" \
        "seeds 21 and 22 give different diagnostics"
    assert_equals "$dd" "$de" "runs without a seed are identical (legacy fixed velocities)"

    print_test_summary
}

main "$@"
