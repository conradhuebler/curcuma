#!/bin/bash
# NEB-MD CLI Test 04 - two-stage pipeline (string -> path-CV conditional mean force)
# Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
# Claude Generated 2026 (AI-generated, machine-tested)
#
# Ethane, GFN-FF, 6 images; stage 1 = scheme 'string' for 2 ps, stage 2 = scheme
# 'independent' + umbrella + constrain on the stage-1 path for 2 ps. Checks: exit 0; both
# stages write into ONE BMT directory; start.stage2.neb.fe.csv exists; stage 2 reports that
# it read the stage-1 path with 6 images.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

pass() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_PASSED=$((TESTS_PASSED + 1)); echo -e "${GREEN}ok PASS${NC}: $1"; }
fail() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); echo -e "${RED}xx FAIL${NC}: $1"; }

run_test() {
    cd "$TEST_DIR"
    # BMT directories of earlier runs are not removed by cleanup_test_artifacts (no nebmd
    # pattern there), so the directory of THIS run is identified as the new one.
    find . -maxdepth 1 -type d -name "start.nebmd.*" | sort > .bmt_before.txt
    $CURCUMA -nebmd start.xyz end.xyz \
        -nebmd.pipeline true -nebmd.nimages 6 -nebmd.string_update 400 \
        -nebmd.string_max_time 2000 -simplemd.max_time 2000 \
        -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.temperature 300 \
        -simplemd.thermostat csvr \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD pipeline should succeed"
    return 0
}

validate_results() {
    local clean=$(sed 's/\x1b\[[0-9;]*m//g' stdout.log)
    local newdirs=$(find . -maxdepth 1 -type d -name "start.nebmd.*" | sort | comm -13 .bmt_before.txt -)
    local ndirs=$(echo "$newdirs" | grep -c "start.nebmd" || true)
    [ "$ndirs" -eq 1 ] && pass "exactly one new BMT directory for both stages" || fail "$ndirs new BMT directories, expected 1"
    local bmt=$(echo "$newdirs" | head -1)
    [ -n "$bmt" ] && [ -f "$bmt/start.neb.string.xyz" ] && pass "stage-1 string path in the BMT directory" \
        || fail "stage-1 start.neb.string.xyz missing in BMT directory"
    [ -n "$bmt" ] && [ -f "$bmt/start.stage2.neb.fe.csv" ] && pass "stage-2 fe.csv in the same BMT directory" \
        || fail "start.stage2.neb.fe.csv missing in BMT directory"
    if echo "$clean" | grep -q "reference path read from .*start.neb.string.xyz (6 images)"; then
        pass "stage 2 read the stage-1 path (6 images)"
    else
        fail "stage 2 did not report reading the stage-1 path with 6 images"
    fi
    # AP4 (Claude Generated 2026): per-point rejection fraction and the interior s range
    if grep -q "^estimator,s,dGds_Eh,G_kcal,samples,rejected_frac,dGds_gradE_Eh,dGds_metric_Eh,dGds_wall_Eh$" "$bmt/start.stage2.neb.fe.csv" 2>/dev/null; then
        pass "fe.csv has the rejected_frac column"
    else
        fail "fe.csv header without rejected_frac"
    fi
    if grep -q "^# interior s range \[" "$bmt/start.stage2.neb.fe.csv" 2>/dev/null \
        && echo "$clean" | grep -q "interior range     : s = \["; then
        pass "interior s range reported (console + fe.csv header)"
    else
        fail "interior s range missing in console or fe.csv"
    fi
    # comment lines are skipped (free-text header)
    local nonfinite=$(grep -v "^#" "$bmt/start.stage2.neb.fe.csv" 2>/dev/null | grep -ciE "nan|inf" || true)
    [ "${nonfinite:-0}" -eq 0 ] && pass "no NaN/Inf in stage-2 fe.csv" || fail "$nonfinite NaN/Inf in stage-2 fe.csv"
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f start.neb.* start.stage2.neb.* .bmt_before.txt; }

main() {
    test_header "NEB-MD pipeline string -> path CV (gfnff, ethane)"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
