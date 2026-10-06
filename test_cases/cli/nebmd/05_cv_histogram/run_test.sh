#!/bin/bash
# NEB-MD CLI Test 05 - cv_histogram (reference G_ref(s) from a trajectory, no MD)
# Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
# Claude Generated 2026 (AI-generated, machine-tested)
#
# Synthetic input (pure python, fixed seed): a 6-frame reference path (linear between
# start.xyz and end.xyz) and a trajectory with every path frame 100 times plus Gaussian
# noise of 0.02 A per coordinate. Expected: the s histogram has 6 separated peaks, their
# populations give G_ref differences < 0.05 kcal/mol (equal populations), and the automatic
# z filter (on z - z_path(s)) finds no far population (unimodal) and discards <= 2 % of the
# frames (only outliers > 10 x median excess).
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_DIR="$SCRIPT_DIR"

pass() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_PASSED=$((TESTS_PASSED + 1)); echo -e "${GREEN}ok PASS${NC}: $1"; }
fail() { TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); echo -e "${RED}xx FAIL${NC}: $1"; }

make_inputs() {
    python3 - <<'PYEOF'
import random
def read(fn):
    L = open(fn).read().split("\n"); n = int(L[0])
    return [(l.split()[0], [float(x) for x in l.split()[1:4]]) for l in L[2:2 + n]]
a, b = read("start.xyz"), read("end.xyz")
N = 6
path = []
for i in range(N):
    t = i / (N - 1)
    path.append([(ea, [(1 - t) * xa[c] + t * xb[c] for c in range(3)]) for (ea, xa), (eb, xb) in zip(a, b)])
def write(f, frame, comment):
    f.write("%d\n%s\n" % (len(frame), comment))
    for e, x in frame:
        f.write("%s %.8f %.8f %.8f\n" % (e, x[0], x[1], x[2]))
with open("path.xyz", "w") as f:
    for i, fr in enumerate(path):
        write(f, fr, "path frame %d" % i)
rng = random.Random(20261001)
with open("traj.xyz", "w") as f:
    for i, fr in enumerate(path):
        for k in range(100):
            write(f, [(e, [x[c] + rng.gauss(0.0, 0.02) for c in range(3)]) for e, x in fr], "frame %d copy %d" % (i, k))
PYEOF
}

run_test() {
    cd "$TEST_DIR"
    make_inputs
    $CURCUMA -nebmd start.xyz end.xyz -nebmd.path_file path.xyz -nebmd.cv_histogram traj.xyz \
        -simplemd.method gfnff -simplemd.temperature 300 > stdout.log 2> stderr.log
    assert_exit_code $? 0 "NEB-MD cv_histogram should succeed"
    local ref=$(find_output_file "start.neb.ref.csv")
    assert_file_exists "$ref" "reference histogram"
    return 0
}

validate_results() {
    local ref=$(find_output_file "start.neb.ref.csv")
    [ -z "$ref" ] || [ ! -f "$ref" ] && return 1
    local res=$(python3 - "$ref" <<'PYEOF'
import math, sys
rows = [l.strip().split(",") for l in open(sys.argv[1]) if l[0].isdigit()]
counts = [int(r[1]) for r in rows]
# a peak = populated bins separated from the next populated bin by fewer than 3 empty bins
# (the 6 frames lie ~16-20 bins apart; single stray samples must not split a peak)
clusters, cur, gap = [], 0, 0
for c in counts + [0, 0, 0]:
    if c > 0:
        cur += c; gap = 0
    else:
        gap += 1
        if gap == 3 and cur > 0:
            clusters.append(cur); cur = 0
kT = 0.0019872043 * 300.0  # kcal/mol
G = [-kT * math.log(c / sum(clusters)) for c in clusters]
print(len(clusters), "%.6f" % (max(G) - min(G)), sum(clusters))
PYEOF
)
    local npk=$(echo "$res" | awk '{print $1}')
    local dg=$(echo "$res" | awk '{print $2}')
    local tot=$(echo "$res" | awk '{print $3}')
    [ "$npk" = "6" ] && pass "histogram has 6 separated peaks" || fail "histogram has $npk peaks, expected 6"
    if awk -v d="$dg" 'BEGIN {exit !(d < 0.05)}'; then
        pass "G_ref of the peaks agree within $dg kcal/mol (< 0.05)"
    else
        fail "G_ref peak spread $dg kcal/mol >= 0.05"
    fi
    # Pure noise around the path is unimodal: the automatic filter must not report a far
    # population and may only drop outliers > 10 x median excess. Accept up to 2 % and
    # require that the rejected count is reported and consistent.
    rej=$(grep '^# z filter:' "$ref" | sed 's/.*rejected \([0-9]*\) = .*/\1/')
    frac=$(grep '^# z filter:' "$ref" | sed 's/.*= *//')
    if [ -n "$rej" ] && [ $((tot + rej)) -eq 600 ] && awk -v f="$frac" 'BEGIN {exit !(f <= 0.02)}'; then
        pass "z filter rejected $rej of 600 frames (fraction $frac <= 0.02); $tot histogrammed"
    else
        fail "z filter: histogrammed $tot, rejected '$rej' (fraction '$frac'): $(grep '# z filter' "$ref")"
    fi
    if grep -q '^# z filter: .*unimodal' "$ref" && grep -q '^# z filter interior: ' "$ref"; then
        pass "z filter reports the unimodal rule and the interior rejection"
    else
        fail "z filter rule/interior line missing: $(grep '# z filter' "$ref")"
    fi
}

cleanup_before() { cd "$TEST_DIR"; cleanup_test_artifacts; rm -f start.neb.* path.xyz traj.xyz; }

main() {
    test_header "NEB-MD cv_histogram (synthetic ethane trajectory)"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
