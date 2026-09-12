#!/bin/bash

# Test: revgfnff reactive-mode NVE energy conservation, RELATIVE to plain gfnff
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026)
#
# Re-specified from an earlier dt^2-scaling design (renamed from
# 18_gfnff_rev_nve_dt2): plain `gfnff` itself does NOT drift as dt^2 on this
# 2 H2 / 3000 K / NVE system (measured on the final binary: 2.06e-4 / 8.3e-5
# Eh at dt 0.25/0.125 fs using the printed "Etot" table column, i.e. close to
# linear in dt, not quadratic -- the wall energy is ~0 throughout this run
# and does not change the picture). A dt^2 criterion on revgfnff alone would
# therefore be testing gfnff's own (pre-existing, unrelated) integrator
# behaviour, not the reactive stage-1b blend this test suite is about.
#
# Instead this test is RELATIVE: run the identical NVE setup with plain
# `gfnff -gfnff.topology_mode static` and with `revgfnff -gfnff.topology_mode
# react`, at two time steps, and require the two methods' energy drift to
# agree to within a factor of ~1.5 either way. If the reactive blend leaked
# energy that gfnff's own integration does not, drift(revgfnff)/drift(gfnff)
# would move far from 1; measured on the final binary it is 1.005 (dt=0.25)
# and 1.012 (dt=0.125) -- i.e. revgfnff's NVE behaviour is indistinguishable
# from plain gfnff's own (imperfect but unrelated) integrator drift here.
#
# Both dt <= the default reactive time-step cap (-md.rev_dt_cap, 0.25 fs --
# "the reactive blend is not integrable at 0.5 fs for hot X-H bonds",
# src/capabilities/simplemd.h), so revgfnff runs at the requested step
# unclamped in this test. A LARGER requested step (e.g. -md.time_step 0.5)
# would be silently clamped to 0.25 fs with a printed
# "WARNING revgfnff: requested time step ... exceeds the reactive-blend
# limit" line (src/capabilities/simplemd.cpp, docs/REV_GFNFF_STAGE1.md);
# -md.rev_dt_cap 0 disables the cap. Not exercised here -- both requested
# steps are already <= the cap, so nothing is clamped and no warning fires.
#
# The drift comparison uses the printed "Etot" table column directly, NOT a
# freshly recomputed Epot+Ekin(+Wall) sum: the two are NOT the same reading
# in the FINAL (finalizeRun) print of a run -- Etot in that one row can be
# one step stale relative to the Epot/Ekin/Wall values printed alongside it
# (a pre-existing SimpleMD print-ordering quirk, unrelated to revgfnff and
# out of scope here -- no src/ changes). Recomputing Epot+Ekin+Wall from that
# stale-adjacent row gives a spurious ~150x ratio; the plain Etot column
# (what both this test and the coordinator's own measurement use) does not
# have that problem because it is one single self-consistent reading.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 18: revgfnff NVE energy conservation vs plain gfnff"
TEST_DIR="$SCRIPT_DIR"

MIN_RATIO=0.5
MAX_RATIO=1.5
MAX_T_K=6000.0

run_one() {
    local sub="$1" method="$2" topo="$3" dt="$4"
    rm -rf "$sub"
    mkdir -p "$sub"
    cp input.xyz "$sub/input.xyz"
    ( cd "$sub" && timeout 260 $CURCUMA -md input.xyz -method "$method" -gfnff.topology_mode "$topo" \
        -temperature 3000 -maxtime 1000 -md.time_step "$dt" \
        -md.wall_radius 3.0 -md.wall_type spheric \
        -md.thermostat none -md.rattle_12 false -md.seed 42 \
        -md.no_restart -threads 1 -verbosity 1 -no_bmt \
        > stdout.log 2> stderr.log )
    return $?
}

run_test() {
    cd "$TEST_DIR"
    rm -rf gfnff_025 gfnff_0125 rev_025 rev_0125
    cleanup_bmt_dirs
    local rc=0
    run_one gfnff_025  gfnff    static 0.25  || rc=$?
    run_one gfnff_0125 gfnff    static 0.125 || rc=$?
    run_one rev_025    revgfnff react  0.25  || rc=$?
    run_one rev_0125   revgfnff react  0.125 || rc=$?
    return $rc
}

validate_results() {
    local failed=0

    # 1. Numerical stability across all four sub-runs
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "Simulation got unstable|NaN/Inf velocity" \
        gfnff_025/stdout.log gfnff_025/stderr.log \
        gfnff_0125/stdout.log gfnff_0125/stderr.log \
        rev_025/stdout.log rev_025/stderr.log \
        rev_0125/stdout.log rev_0125/stderr.log 2>/dev/null; then
        echo -e "${RED}✗ FAIL${NC}: MD reported instability in at least one run"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: no instability reported in any run"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 2. drift(revgfnff)/drift(gfnff) in [0.5, 1.5] at both dt, T_max < 6000 K
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$MIN_RATIO" "$MAX_RATIO" "$MAX_T_K" <<'PYEOF'
import re, sys
min_ratio, max_ratio, max_t = (float(x) for x in sys.argv[1:4])
ROW_RE = re.compile(r"\s+\d+\.\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+")

def load(path):
    txt = open(path).read()
    return [l.split() for l in txt.splitlines() if ROW_RE.match(l)]

def finite_row(r):
    try:
        vals = [float(r[i]) for i in (1, 3, 5, 7, 9)]
    except ValueError:
        return False
    return all(v == v for v in vals)  # NaN != NaN

stats = {}
any_nonfinite = False
for name in ["gfnff_025", "gfnff_0125", "rev_025", "rev_0125"]:
    rows = load(f"{name}/stdout.log")
    if not rows:
        print(f"FAIL: no MD table rows parsed for {name}")
        sys.exit(1)
    good = [r for r in rows if finite_row(r)]
    if len(good) < 2:
        print(f"FAIL: fewer than 2 finite MD table rows for {name}")
        sys.exit(1)
    if len(good) != len(rows):
        any_nonfinite = True
    # column 5 (0-indexed) is the printed "Etot" table column -- used as-is,
    # see header comment for why a recomputed Epot+Ekin+Wall sum is NOT used.
    drift = abs(float(good[-1][5]) - float(good[0][5]))
    tmax = max(float(r[7]) for r in good)
    stats[name] = (drift, tmax)

ratio_025 = stats["rev_025"][0] / stats["gfnff_025"][0] if stats["gfnff_025"][0] > 0 else float("inf")
ratio_0125 = stats["rev_0125"][0] / stats["gfnff_0125"][0] if stats["gfnff_0125"][0] > 0 else float("inf")
t_max_all = max(v[1] for v in stats.values())

print(f"drift gfnff(0.25)={stats['gfnff_025'][0]:.3e} rev(0.25)={stats['rev_025'][0]:.3e} ratio={ratio_025:.4f} | "
      f"gfnff(0.125)={stats['gfnff_0125'][0]:.3e} rev(0.125)={stats['rev_0125'][0]:.3e} ratio={ratio_0125:.4f} | "
      f"T_max={t_max_all:.1f} any_nonfinite_row={any_nonfinite}")

ok = True
reasons = []
if not (min_ratio <= ratio_025 <= max_ratio):
    ok = False; reasons.append(f"ratio(dt=0.25) {ratio_025:.4f} outside [{min_ratio},{max_ratio}]")
if not (min_ratio <= ratio_0125 <= max_ratio):
    ok = False; reasons.append(f"ratio(dt=0.125) {ratio_0125:.4f} outside [{min_ratio},{max_ratio}]")
if t_max_all > max_t:
    ok = False; reasons.append(f"T_max {t_max_all:.1f} > {max_t}")
if any_nonfinite:
    ok = False; reasons.append("non-finite MD table row(s) beyond the known t=0 cosmetic one")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: revgfnff/gfnff drift ratio within [${MIN_RATIO}, ${MAX_RATIO}] at both dt"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: revgfnff/gfnff drift ratio or T_max outside thresholds"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    assert_exit_code $? 0 "all four sub-runs should complete without crash"
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
