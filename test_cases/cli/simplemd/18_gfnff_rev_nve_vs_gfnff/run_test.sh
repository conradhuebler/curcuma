#!/bin/bash

# Test: revgfnff reactive-mode NVE energy conservation, RELATIVE to plain gfnff
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026)
#
# Re-specified a second time (Sep 12, 2026). The previous version of this
# test read the printed "Etot" table column at the FIRST and LAST rows and
# required the two methods' endpoint drift to agree within a factor of 1.5,
# with an explicit note that a recomputed Epot+Ekin(+Wall) sum gave "a
# spurious ~150x ratio" and was therefore avoided. That note described a
# real bug, not a property of the physics: commit 4ef8d30c fixed
# finalizeRun() printing a STALE Etot in the final table row (the value from
# the previous periodic print, up to print_frequency old) instead of the
# current Epot+Ekin. The test was calibrated against that bug -- with it
# fixed, the endpoint-drift criterion itself is unreliable for a different,
# pre-existing reason: NVE drift on a 4-atom system over ~1 ps is dominated
# by vibrational-phase noise, not by the integrator's true energy leak (see
# CLAUDE.md Known Issue #32). An endpoint difference can be accidentally
# tiny (denominator near a phase-cancellation point) or accidentally large
# for either method, unrelated to any real drift rate. Measured on the
# then-current binary: drift(revgfnff)/drift(gfnff) = 101.3 at dt=0.25 fs,
# because gfnff's OWN endpoint drift over that 1 ps window happened to be
# 3.0e-6 Eh (an accident of oscillation phase), not because revgfnff leaked
# 100x more energy (its own endpoint drift there was an unremarkable 3.0e-4
# Eh). At dt=0.125 fs the same run gave a normal-looking ratio of 1.06 --
# i.e. the metric is not just biased, it is unstable in both directions.
#
# This version measures the SLOPE of Etot(t) by ordinary least squares over
# the whole trajectory, not the two endpoints, and grounds the pass/fail
# threshold in a real calibration (below) instead of a single ratio window.
#
# CALIBRATION (2 H2 / 3000 K / NVE / seed 42, this binary, Sep 12, 2026):
# -maxtime 10000 (10 ps) at each dt, with -md.print_frequency 100 (0.1 ps)
# for enough points to fit; the fit itself de-duplicates to one row per
# 0.1 ps bucket first (SimpleMD's print condition is step-count based and
# can print a short burst of near-identical rows around each interval --
# fitting the raw stream would understate the noise). Slopes in Eh/ps:
#
#   run         slope        note
#   gfnff  dt=0.25   -1.747e-04   baseline (no reactions)
#   gfnff  dt=0.125  -8.825e-05   baseline (no reactions)
#   rev    dt=0.25   -1.292e-04   2 bonds formed+broken, sum -0.7 kJ/mol
#   rev    dt=0.125  +1.214e-04   1 bond formed+broken,  sum -0.2 kJ/mol
#
# (revgfnff's react mode legitimately forms/breaks a transient H-H bond
# once or twice in this window at 3000 K -- that is the mechanism under
# test behaving as designed, not noise; REACT summary lines in the log
# report each jump's energy and are small, 0.1-1 kJ/mol here.) Fit standard
# errors are 5.4e-05 to 8.0e-05 Eh/ps -- i.e. the same order as the slopes
# themselves on this tiny, single-seed, partly-reactive system; a run long
# enough to shrink that further (30-250 ps was tried) starts crossing into
# genuine H-exchange events that push T_max toward/over 6000 K, which is a
# real chemistry artifact of this cramped 4-atom box, not an integrator
# problem, and would make the test fail for the wrong reason. 10 ps with
# denser printing is the shortest window where the slope sign and order of
# magnitude are stable and the run stays deep inside the safe T_max regime
# (largest T_max observed here: 3452 K).
#
# THRESHOLD: |slope_rev| must not exceed 1.5x |slope_gfnff| OR an absolute
# floor of 5.0e-4 Eh/ps, whichever is larger. The floor is set to ~3x the
# largest slope magnitude (1.75e-4) and ~7x the largest fit standard error
# (8.0e-5) measured above, so ordinary seed/phase noise -- including the
# occasional small react bond-formation/breaking event -- cannot trip it,
# while a leak an order of magnitude worse than anything measured here
# still fails. This is a one-sided check: revgfnff drifting LESS than
# gfnff is not a problem (rev/dt=0.125 above has |slope_rev| < |slope_gfnff|
# only by a factor of ~1.4, comfortably inside the 1.5x margin anyway).
#
# Both dt <= the default reactive time-step cap (-md.rev_dt_cap, 0.25 fs --
# "the reactive blend is not integrable at 0.5 fs for hot X-H bonds",
# src/capabilities/simplemd.h), so revgfnff runs at the requested step
# unclamped in this test.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 18: revgfnff NVE energy conservation vs plain gfnff"
TEST_DIR="$SCRIPT_DIR"

MAXTIME_FS=10000        # 10 ps, see calibration above
PRINT_FREQUENCY_FS=100  # 0.1 ps -- dense enough for a determined slope fit
FIT_BUCKET_PS=0.1        # de-duplication bucket for the fit, matches print frequency
RATIO_FACTOR=1.5
SLOPE_FLOOR=5.0e-4       # Eh/ps, see calibration above
MAX_T_K=6000.0

run_one() {
    local sub="$1" method="$2" topo="$3" dt="$4"
    rm -rf "$sub"
    mkdir -p "$sub"
    cp input.xyz "$sub/input.xyz"
    ( cd "$sub" && timeout 260 $CURCUMA -md input.xyz -method "$method" -gfnff.topology_mode "$topo" \
        -temperature 3000 -maxtime "$MAXTIME_FS" -md.time_step "$dt" \
        -md.print_frequency "$PRINT_FREQUENCY_FS" \
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

    # 2. slope(revgfnff) <= max(RATIO_FACTOR * slope(gfnff), SLOPE_FLOOR) at both dt,
    #    no non-finite table row beyond the known cosmetic t=0 one, T_max < 6000 K.
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$RATIO_FACTOR" "$SLOPE_FLOOR" "$MAX_T_K" "$FIT_BUCKET_PS" <<'PYEOF'
import re, sys, math

ratio_factor, slope_floor, max_t, bucket_ps = (float(x) for x in sys.argv[1:5])
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

def linfit_slope(t, y):
    # Ordinary least squares slope of y vs t. Used only after de-duplicating
    # to one row per time bucket (see header) so consecutive near-identical
    # printed rows do not pose as independent evidence.
    n = len(t)
    tm = sum(t) / n
    ym = sum(y) / n
    sxx = sum((ti - tm) ** 2 for ti in t)
    if sxx == 0:
        return 0.0
    sxy = sum((ti - tm) * (yi - ym) for ti, yi in zip(t, y))
    return sxy / sxx

def thin(rows, bucket):
    seen = set()
    out = []
    for r in rows:
        b = int(float(r[0]) / bucket)
        if b not in seen:
            seen.add(b)
            out.append(r)
    return out

stats = {}
any_nonfinite = False
for name in ["gfnff_025", "gfnff_0125", "rev_025", "rev_0125"]:
    rows = load(f"{name}/stdout.log")
    if not rows:
        print(f"FAIL: no MD table rows parsed for {name}")
        sys.exit(1)
    good = [r for r in rows if finite_row(r)]
    if len(good) < 3:
        print(f"FAIL: fewer than 3 finite MD table rows for {name}")
        sys.exit(1)
    if len(good) != len(rows):
        any_nonfinite = True
    thinned = thin(good, bucket_ps)
    t = [float(r[0]) for r in thinned]
    Etot = [float(r[5]) for r in thinned]
    slope = linfit_slope(t, Etot)
    tmax = max(float(r[7]) for r in good)
    stats[name] = (slope, tmax, len(thinned))

t_max_all = max(v[1] for v in stats.values())

def threshold_check(dt_label, gfnff_key, rev_key):
    slope_g, _, n_g = stats[gfnff_key]
    slope_r, _, n_r = stats[rev_key]
    threshold = max(ratio_factor * abs(slope_g), slope_floor)
    ok = abs(slope_r) <= threshold
    print(f"dt={dt_label}: slope gfnff={slope_g:.4e} (n={n_g}) rev={slope_r:.4e} (n={n_r}) "
          f"threshold={threshold:.4e} (1.5x={ratio_factor*abs(slope_g):.4e}, floor={slope_floor:.4e}) "
          f"{'OK' if ok else 'EXCEEDED'}")
    return ok

ok_025 = threshold_check("0.25", "gfnff_025", "rev_025")
ok_0125 = threshold_check("0.125", "gfnff_0125", "rev_0125")
print(f"T_max={t_max_all:.1f} any_nonfinite_row={any_nonfinite}")

ok = True
reasons = []
if not ok_025:
    ok = False; reasons.append("slope(rev, dt=0.25) exceeds threshold")
if not ok_0125:
    ok = False; reasons.append("slope(rev, dt=0.125) exceeds threshold")
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
        echo -e "${GREEN}✓ PASS${NC}: revgfnff slope within threshold of plain gfnff's slope at both dt, T_max ok"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: revgfnff slope, T_max, or table finiteness outside thresholds"
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
