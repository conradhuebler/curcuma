#!/bin/bash
# NEB-MD CLI Test 07 - string stage aborts on a dead image (string_abort_on_dead)
# Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
# Claude Generated 2026 (AI-generated, machine-tested)
#
# Ethane (staggered -> eclipsed), GFN-FF, 6 images. A dead image is forced with
# -nebmd.max_temperature 350 at a 300 K setpoint: the instantaneous temperature of an 8-atom
# image fluctuates by ~100 K, so a window crosses 350 K within the first ~20 steps.
#  (a) scheme string: exit code != 0, error mentions 'dead', no start.neb.string.xyz;
#  (b) pipeline: exit code != 0, 'dead', no start.stage2.* output and no string path;
#  (c) string_abort_on_dead false: exit 0, the run reports the frames that are NOT relaxed.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

pass() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_PASSED=$((TESTS_PASSED + 1)); echo -e "${GREEN}ok PASS${NC}: $1"; }
fail() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); echo -e "${RED}xx FAIL${NC}: $1"; }

COMMON="-nebmd.nimages 6 -nebmd.string_update 400 -nebmd.string_max_time 2000 \
 -nebmd.max_temperature 350 -simplemd.method gfnff -simplemd.time_step 0.5 \
 -simplemd.temperature 300 -simplemd.thermostat csvr"

# $1 = label, rest = extra arguments; leaves the cleaned stdout in clean_<label>.txt, returns rc
run_case() {
    local label="$1"; shift
    cd "$TEST_DIR"
    cleanup_test_artifacts; rm -f start.neb.* start.stage2.* clean_$label.txt; rm -rf start.nebmd.*
    local rc=0
    $CURCUMA -nebmd start.xyz end.xyz "$@" $COMMON > stdout_$label.txt 2>&1 || rc=$?
    sed 's/\x1b\[[0-9;]*m//g' stdout_$label.txt > clean_$label.txt
    return $rc
}

count_files() { find . -maxdepth 2 -name "$1" 2>/dev/null | wc -l; }

case_string() {
    local rc=0
    run_case string -nebmd.scheme string || rc=$?
    [ "$rc" -ne 0 ] && pass "(a) string: exit code $rc != 0" || fail "(a) string: exit code 0 despite dead image"
    grep -q "ABORTED - dead image" clean_string.txt && pass "(a) string: error message mentions 'dead'" \
        || fail "(a) string: no 'ABORTED - dead image' message"
    local n=$(count_files "start.neb.string.xyz")
    [ "$n" -eq 0 ] && pass "(a) string: no start.neb.string.xyz written" || fail "(a) string: $n string path file(s) written"
}

case_pipeline() {
    local rc=0
    run_case pipeline -nebmd.pipeline true -simplemd.max_time 400 || rc=$?
    [ "$rc" -ne 0 ] && pass "(b) pipeline: exit code $rc != 0" || fail "(b) pipeline: exit code 0 despite dead image"
    grep -q "dead" clean_pipeline.txt && pass "(b) pipeline: output mentions 'dead'" || fail "(b) pipeline: no 'dead' in output"
    grep -q "stage 2 skipped" clean_pipeline.txt && pass "(b) pipeline: stage 2 skipped" || fail "(b) pipeline: no 'stage 2 skipped'"
    local n2=$(count_files "*.stage2.*")
    [ "$n2" -eq 0 ] && pass "(b) pipeline: no *.stage2.* output" || fail "(b) pipeline: $n2 *.stage2.* file(s)"
    local n=$(count_files "start.neb.string.xyz")
    [ "$n" -eq 0 ] && pass "(b) pipeline: no start.neb.string.xyz written" || fail "(b) pipeline: $n string path file(s)"
}

case_no_abort() {
    local rc=0
    run_case noabort -nebmd.scheme string -nebmd.string_abort_on_dead false || rc=$?
    [ "$rc" -eq 0 ] && pass "(c) abort_on_dead false: exit code 0" || fail "(c) abort_on_dead false: exit code $rc"
    grep -q "dead images: path frames .* are NOT relaxed" clean_noabort.txt \
        && pass "(c) abort_on_dead false: end of stage reports frames NOT relaxed" \
        || fail "(c) abort_on_dead false: no 'NOT relaxed' report"
}

main() {
    test_header "NEB-MD string: abort on dead image (gfnff, ethane)"
    case_string
    case_pipeline
    case_no_abort
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
