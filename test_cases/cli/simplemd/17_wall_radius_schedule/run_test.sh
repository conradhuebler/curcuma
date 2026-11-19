#!/bin/bash
# SimpleMD wall radius schedule and wall pressure CSV columns - Claude Generated (Nov 2026)
# He2 at 10 Angstrom (no interaction), centred on the origin by SimpleMD: both atoms sit at +-5 A.
#  A: harmonic wall R = 3 A, k = wall_temp * kB; wall force per atom k (5 - R), pressure 2 F / (4 pi R^2)
#  B: wall_radius_schedule 4 A -> 3 A over 10 of 20 steps: radius falls monotonically, is then held at 3.0, pressure rises
#  C: no wall: no wall columns in the CSV (existing files unchanged)
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

check() {  # check <description> <1 = pass>
    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$2" = "1" ]; then
        echo -e "${GREEN}PASS${NC}: $1"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}: $1"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

make_input() { printf '2\nHe2\nHe 0.0 0.0 0.0\nHe 10.0 0.0 0.0\n' > input.xyz; }

# prints the BMT-or-CWD CSV path of the current directory
csv() { find . -maxdepth 2 -name "input.md.csv" -type f | head -1; }

run_test() {
    cd "$TEST_DIR"
    rm -rf run_a run_b run_c
    mkdir run_a run_b run_c

    cd run_a; make_input
    $CURCUMA -md input.xyz -method gfn2 -threads 1 -wall_type spheric -wall_potential harmonic -wall_radius 3.0 \
        -wall_temp 1000 -temperature 1 -max_time 0.2 -dt 0.1 -print 1 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "A: MD should succeed"
    cd ..

    cd run_b; make_input
    $CURCUMA -md input.xyz -method gfn2 -threads 1 -wall_type spheric -wall_potential harmonic -wall_radius 4.0 \
        -wall_temp 1000 -temperature 1 -max_time 4 -dt 0.2 -print 1 -wall_radius_schedule "3.0:steps:10" \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "B: MD should succeed"
    cd ..

    cd run_c; make_input
    $CURCUMA -md input.xyz -method gfn2 -threads 1 -temperature 1 -max_time 0.2 -dt 0.1 -print 1 \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "C: MD should succeed"
    cd ..
    return 0
}

validate_results() {
    # A: last column triple = radius, P, P_avg. Expected P = 2 k (5 - R) / (4 pi R^2) * 4359.7447222 GPa
    cd "$TEST_DIR/run_a"
    local f=$(csv)
    assert_file_exists "$f" "A: CSV file"
    local h=0
    head -1 "$f" | grep -q "wall_radius_A;wall_pressure_GPa;wall_pressure_avg100_GPa" && h=1
    check "A: header carries wall_radius_A;wall_pressure_GPa;wall_pressure_avg100_GPa" $h
    # row of step 1: atoms moved < 1e-6 A, so within 0.2 % of the analytic value
    local ok=$(awk -F';' 'NR>1 && $1==1 { k=1000*3.1668115634556e-6; F=2*k*(5.0-3.0);
        pe=F/(4*3.14159265358979*9.0)*4359.7447222071; d=($(NF-1)-pe)/pe; if (d<0) d=-d;
        printf "%d", (d<2e-3 && $(NF-2)==3.0) }' "$f")
    check "A: wall_pressure_GPa(step 1) equals 2 k (r-R)/(4 pi R^2) (0.2 %), radius 3.0" "${ok:-0}"
    cd "$TEST_DIR/run_b"
    f=$(csv)
    assert_file_exists "$f" "B: CSV file"
    ok=$(awk -F';' 'NR>1 { r=$(NF-2); p=$(NF-1); if (NR==2) {p0=p} else if (r>rprev+1e-9) bad=1; rprev=r; last_r=r; last_p=p }
        END { d=last_r-3.0; if (d<0) d=-d; printf "%d", (!bad && d<1e-6 && last_p>p0) }' "$f")
    check "B: radius falls monotonically to 3.0, pressure rises" "${ok:-0}"
    cd "$TEST_DIR/run_c"
    f=$(csv)
    assert_file_exists "$f" "C: CSV file"
    assert_string_not_in_file "wall_radius_A" "$f" "C: no wall columns without a wall"
}

cleanup_before() { cd "$TEST_DIR"; rm -rf run_a run_b run_c; }

main() {
    test_header "SimpleMD Wall Radius Schedule Test"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
