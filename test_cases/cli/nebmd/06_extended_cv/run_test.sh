#!/bin/bash
# NEB-MD CLI Test 06 - extended / re-sampled path CV (cv_frames, cv_extend, cv_extend_length) + basin barrier
# Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
# Claude Generated 2026 (AI-generated, machine-tested)
#
# Ethane, GFN-FF, scheme 'independent' + umbrella + constrain, 6 windows, path CV re-sampled
# onto 12 frames and extended by 2 frames beyond each end (16 CV frames), 2 ps per window.
# Checks: exit 0; fe.csv exists; the support points reach beyond s(second path frame) and
# s(second-to-last path frame) - the outer windows sit in the former endpoint basins, which
# the default setup (fixed endpoints) never samples; 'barrier (basin' lines in console and
# fe.csv; every number in fe.csv finite; no dead window (6 estimator-I support points).
# Second run (subdir len_run/): same setup with cv_frames 24 and the extension given as a LENGTH,
# cv_extend_length 0.06 -> m = ceil(0.06 * 23) + 1 = 3 frames per end (30 CV frames). Checks:
# the 'cv_extend: 0.06 of the path length = 3 frames' log line, the CV frame line, support
# points outside [s(second), s(second-to-last)], no dead window, 6 support points.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

pass() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_PASSED=$((TESTS_PASSED + 1)); echo -e "${GREEN}ok PASS${NC}: $1"; }
fail() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); echo -e "${RED}xx FAIL${NC}: $1"; }

run_test() {
    cd "$TEST_DIR"
    $CURCUMA -nebmd start.xyz end.xyz \
        -nebmd.scheme independent -nebmd.umbrella true -nebmd.constrain true \
        -nebmd.nimages 6 -nebmd.cv_frames 12 -nebmd.cv_extend 2 \
        -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.temperature 300 \
        -simplemd.thermostat csvr -simplemd.max_time 2000 -simplemd.seed 1 \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD with extended CV should succeed"
    return 0
}

run_test_length() {
    mkdir -p "$TEST_DIR/len_run"
    cd "$TEST_DIR/len_run"
    cp ../start.xyz ../end.xyz .
    $CURCUMA -nebmd start.xyz end.xyz \
        -nebmd.scheme independent -nebmd.umbrella true -nebmd.constrain true \
        -nebmd.nimages 6 -nebmd.cv_frames 24 -nebmd.cv_extend_length 0.06 \
        -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.temperature 300 \
        -simplemd.thermostat csvr -simplemd.max_time 2000 -simplemd.seed 1 \
        > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD with cv_extend_length should succeed"
    return 0
}

# Claude Generated 2026: $1 = run directory, $2 = expected 'path CV on ...' fragment
validate_results() {
    cd "$1"
    local cvline="$2"
    local clean=$(sed 's/\x1b\[[0-9;]*m//g' stdout.log)
    local fe=$(find_output_file "start.neb.fe.csv")
    if [ -z "$fe" ] || [ ! -f "$fe" ]; then
        fail "start.neb.fe.csv not written"
        return 1
    fi
    pass "fe.csv written ($fe)"
    echo "$clean" | grep -qF "$cvline" \
        && pass "CV line: $cvline" || fail "CV frame count line missing / wrong (expected '$cvline')"
    # support points vs the original path frames on the same CV (from the fe.csv header)
    local res=$(python3 - "$fe" <<'PYEOF'
import math, re, sys
lines = open(sys.argv[1]).read().splitlines()
cv = [l for l in lines if l.startswith("# cv path:")]
m = re.search(r"s\(first\) ([-0-9.e+]+), s\(second\) ([-0-9.e+]+), s\(second-to-last\) ([-0-9.e+]+), s\(last\) ([-0-9.e+]+)", cv[0]) if cv else None
rows = [l.split(",") for l in lines if l.startswith("I,")]
s = [float(r[1]) for r in rows]
finite = all(math.isfinite(float(x)) for l in lines if not l.startswith("#") and not l.startswith("estimator")
             for x in l.split(",")[1:])
if not m or not s:
    print("nocv", len(s), int(finite)); sys.exit(0)
s1, s2 = float(m.group(2)), float(m.group(3))
print("ok", len(s), int(finite), int(min(s) < s1), int(max(s) > s2), "%.4f %.4f %.4f %.4f" % (min(s), s1, max(s), s2))
PYEOF
)
    local tag=$(echo "$res" | awk '{print $1}')
    local nI=$(echo "$res" | awk '{print $2}')
    local fin=$(echo "$res" | awk '{print $3}')
    [ "$tag" = "ok" ] && pass "fe.csv carries the '# cv path:' header" || fail "fe.csv without '# cv path:' header ($res)"
    [ "$(echo "$res" | awk '{print $4}')" = "1" ] && pass "lowest support point below s(second path frame) ($(echo "$res" | awk '{print $6, "<", $7}'))" \
        || fail "no support point below s(second path frame): $res"
    [ "$(echo "$res" | awk '{print $5}')" = "1" ] && pass "highest support point above s(second-to-last path frame) ($(echo "$res" | awk '{print $8, ">", $9}'))" \
        || fail "no support point above s(second-to-last path frame): $res"
    [ "$fin" = "1" ] && pass "every number in fe.csv is finite" || fail "non-finite value in fe.csv"
    if grep -q "^# estimator I barrier (basin, fwd/rev)" "$fe" && echo "$clean" | grep -q "barrier (basin, fwd/rev)"; then
        pass "'barrier (basin' lines in fe.csv and console"
    else
        fail "'barrier (basin' line missing in fe.csv or console"
    fi
    # basin numbers, when defined, must be finite
    local bad=$(grep "^# estimator I \(barrier (basin\|reaction dG (basin\)" "$fe" | grep -ciE "nan|inf" || true)
    [ "${bad:-0}" -eq 0 ] && pass "basin barrier values finite" || fail "non-finite basin barrier value"
    if echo "$clean" | grep -q "declared dead"; then
        fail "a window was declared dead"
    else
        pass "no dead window"
    fi
    [ "$nI" = "6" ] && pass "all 6 windows contributed a support point" || fail "$nI estimator-I support points, expected 6"
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f start.neb.*; rm -rf len_run start.nebmd.*; }

main() {
    test_header "NEB-MD extended path CV (cv_frames 12, cv_extend 2 / cv_frames 24, cv_extend_length 0.06; gfnff, ethane)"
    cleanup_before
    run_test && validate_results "$TEST_DIR" "path CV on 16 frames (path 6 frames, cv_frames 12, cv_extend 2)"
    echo "--- second run: cv_frames 24, cv_extend_length 0.06 ---"
    if run_test_length; then
        cd "$TEST_DIR/len_run"
        sed 's/\x1b\[[0-9;]*m//g' stdout.log | grep -q "cv_extend: 0.06 of the path length = 3 frames" \
            && pass "log: 'cv_extend: 0.06 of the path length = 3 frames'" \
            || fail "log line 'cv_extend: 0.06 of the path length = 3 frames' missing"
        validate_results "$TEST_DIR/len_run" "path CV on 30 frames (path 6 frames, cv_frames 24, cv_extend 3)"
    fi
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
