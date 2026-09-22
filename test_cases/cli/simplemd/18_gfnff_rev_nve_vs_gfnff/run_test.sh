#!/bin/bash

# Test: revgfnff reactive-mode NVE energy conservation, RELATIVE to plain gfnff
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026)
#
# Third calibration (Sep 13, 2026). History, so the numbers below are not read
# as arbitrary:
#
# 1. The test originally compared the PRINTED Etot of the first and last table
#    row. Commit 4ef8d30c fixed finalizeRun() printing a STALE Etot in the final
#    row (the value of the previous periodic print), which the criterion had
#    been calibrated against; with it fixed the endpoint ratio was shown to be
#    unstable in both directions (gfnff's denominator was accidentally ~100x
#    below its typical drift on one run). Replaced by the slope of Etot(t).
# 2. The slope version used 2 H2 (4 atoms) at 3000 K. Commit 4601be27 tightened
#    the react formation criterion (new PARAM -gfnff.rev_form_switch, default
#    "order" = join on the NARROW bond order b2_ij > rev_bo2_form = 0.1, which
#    crosses at 1.611x the covalent radius sum; "weight" = the old WIDE
#    term-weight switch at rev_bo_form = 0.05, which still reads 0.59 at a
#    hydrogen bond). Under the new criterion that system produces ZERO react
#    events -- the test then compared a revgfnff run in which the reactive
#    blend never fired against plain gfnff and measured nothing about the
#    blend. This version fixes that.
#
# WHY THE SYSTEM AND THE SETTING CHANGED (measured, not assumed):
#
# * 2 H2 in a 3.0 A wall at any temperature up to 6000 K: 0 events in 10 ps
#   NVE. The narrow switch needs the two H atoms of DIFFERENT molecules within
#   1.611 x (0.32+0.32) x 1.02^2 = 1.073 A; the measured closest approach over
#   those runs is 1.30-1.66 A. Going hotter does not help for a different
#   reason: an NVE H2 has a fixed vibrational energy, and reaching the break
#   radius (~1.07 A) needs ~45 kcal/mol in the stretch, i.e. >10000 K. Between
#   10000 K and 20000 K the events that do fire are all break/re-form cycles of
#   an EXISTING H2 bond (the fresh-pair formation path stays out of reach).
# * The spherical wall is a NON-CONSERVATIVE element of this test: it is
#   applied as a velocity impulse and its potential energy is not part of the
#   printed Epot (see the "Etot invariant" note in
#   test_cases/revgfnff/_log/NVE_TEST_STATUS.md). Measured on the 2 H2 system at
#   12000 K, 10 ps: slope(gfnff, wall) = -1.4e-3 Eh/ps but slope(gfnff, no
#   wall) = +1.4e-5 Eh/ps -- the wall contributes ~100x more drift than the
#   integrator and the reactive blend together. The test therefore runs
#   WITHOUT a wall (-md.wall_type none). The monomers then separate (the
#   reactivity under test is intramolecular H2 break/re-form, which needs no
#   confinement); nothing about the energy test depends on them staying close.
# * A single H2 pair is at a threshold and its event count is erratic (12000 K /
#   dt 0.25 gives 3416 rebuilds, dt 0.125 gives 0). A bath of 12 H2, sparsely
#   packed so that the molecules do not collide on the way out, samples the
#   threshold across 12 independent molecules and is stable instead: 16-52
#   rebuilds at BOTH time steps for every temperature from 10500 K to 12500 K.
#
# CALIBRATION (12 H2 / 11500 K / NVE / no wall / seed 42, this binary,
# Sep 13, 2026): -maxtime 10000 (10 ps) at each dt, -md.print_frequency 100
# (0.1 ps). The fit de-duplicates to one row per 0.1 ps first (SimpleMD's print
# condition is int(step*dt) % print_frequency and emits a short burst of 2-4
# near-identical rows around each interval -- fitting the raw stream would
# understate the noise). Slopes in Eh/ps:
#
#   run         slope        SE          events (formed/broken/rebuilds)
#   gfnff  dt=0.25   -9.478e-06   8.33e-05   0 / 0 / 0
#   gfnff  dt=0.125  +2.704e-05   1.78e-05   0 / 0 / 0
#   rev    dt=0.25   +3.851e-04   2.31e-04   14 / 15 / 46
#   rev    dt=0.125  +3.932e-04   2.31e-04   20 / 21 / 52
#
# T_max = 11845 K (1.03x the target). dE_jump over the 46/52 rebuilds: median
# 0.0, |max| 1.1 kJ/mol, sum -1.9 kJ/mol. Wall time 17.5 s for all four
# sub-runs on one thread (TIMEOUT 600 s is ample, unchanged).
#
# The event count was checked to be robust rather than a knife edge: the same
# four runs at 10500 / 11000 / 12000 / 12500 K give 24-50 / 20-50 / 16-44 /
# 16-38 rebuilds at the two dt values, with slopes 3.5e-4 to 5.5e-4 and T_max
# always ~1.03x the target -- no instability message anywhere.
#
# THRESHOLD: |slope_rev| must not exceed 1.5x |slope_gfnff| OR an absolute
# floor of 1.6e-3 Eh/ps, whichever is larger. The floor is set to ~4x the
# largest slope magnitude (3.93e-4) and ~7x the largest fit standard error
# (2.31e-4) measured above, so ordinary trajectory/phase noise and the
# event-count spread across the 10500-12500 K window cannot trip it, while a
# blend leak four times larger than anything measured here still fails. The
# relative term is what carries the physics: gfnff's own slope is ~1e-5, so a
# revgfnff slope an order of magnitude above its baseline fails on that term
# alone even before the floor is reached. This is a one-sided check: revgfnff
# drifting LESS than gfnff is not a problem.
#
# REBUILDS: the test additionally requires at least 20 REACT rebuilds summed
# over the two revgfnff runs (measured: 98 at the commit setting, never below
# 54 in the 10500-12500 K window). Without it the test could silently go back
# to measuring nothing -- a slope of ~0 passes trivially -- which is exactly the
# failure mode this re-specification repairs.
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

# ---------------------------------------------------------------------------
# Third re-specification (Sep 18, 2026): stage 3a made the same bath far more
# reactive, so the OPERATING POINT moved, not the criterion
# ---------------------------------------------------------------------------
# Measured with the stage-3a default (rev_budget_fix_h true, share delivered, well gauss),
# same 12-H2 bath, same 10 ps / 0.1 ps protocol, seed 42, over a temperature window
# (test_cases/revgfnff/_log/WORK_STATUS.md package 5):
#
#     T / K   dt      slope_rev Eh/ps   SE        rebuilds   slope_gfnff
#      8000   0.25    6.372e-04         3.6e-04         38   5.5e-05
#      8000   0.125   6.230e-04         3.6e-04         52   4.2e-06
#      9500   0.25    1.766e-03         7.6e-04       3274  -1.1e-05
#      9500   0.125   2.710e-03         9.0e-04        260   1.5e-06
#     10500   0.25    2.521e-03         1.2e-03        188  -7.6e-05
#     11500   0.25    4.311e-03         1.7e-03        568  -9.5e-06
#     11500   0.125   6.183e-03         2.1e-03        592   2.7e-05
#     12500   0.125   4.375e-03         1.7e-03        422   8.8e-06
#
# At the committed 11500 K the bath now produces 568-592 rebuilds where the Sep-13 calibration
# measured 16-52, and the drift scales with the event count: those events are real topology
# changes with a real dE_jump each, not integrator error. Raising the floor to cover 11500 K
# would mean 2.5e-2 Eh/ps, i.e. a test that gates almost nothing.
#
# So the TEMPERATURE moves to 8000 K, where the event count (38-52) is back in the band the
# Sep-13 floor was calibrated for, and the floor is re-derived by the same rule: ~4x the largest
# measured |slope_rev| (6.37e-04 -> 2.5e-03) and ~7x the largest fit SE (3.6e-04 -> 2.5e-03).
# The relative 1.5x term is unchanged and is still what carries the physics (gfnff's own slope
# is ~5e-05 here). MAX_T_K keeps its 1.3x-of-target meaning: 8373/8410 K measured -> 10400.
# ---------------------------------------------------------------------------
# FLAGGED, NOT RECALIBRATED (Sep 22, 2026, WORK_STATUS package 12): this test FAILS at the
# mg3 default on its dt = 0.125 arm, and the decision on what to do about it is the operator's.
# Nothing below was changed. What was measured, all on the same bath/protocol/binary:
#
#   arm                          dt 0.25              dt 0.125
#   committed calibration        6.37e-4 /  38 reb    6.23e-4 /  52 reb   (OLD MD clock)
#   HEAD pre-flip  (mg)          1.00e-3 /  70 reb    1.85e-3 / 190 reb
#   HEAD post-flip (mg3)         1.17e-3 /  74 reb    2.95e-3 / 244 reb   <- floor 2.5e-3
#
# Two separate things are visible here and they should not be conflated:
#
# 1. The committed calibration is STALE FOR A REASON THAT PREDATES THE WELL-FORM FLIP. Package
#    10's MD time-step unit fix changed what "-maxtime 10000" and "-md.time_step 0.125" mean, so
#    the same mg arm now produces 70/190 rebuilds where the Sep-18 table recorded 38/52, and its
#    dt = 0.125 slope sits at 0.74x of the floor where this file's own rule ("~4x the largest
#    measured slope") asks for 0.25x. The test only still passed pre-flip because of that 26 %
#    of headroom.
# 2. On top of that, mg3 IS reproducibly more dissipative on this particular bath. Measured with
#    8 paired replicates (a 1e-5 A displacement; -md.seed does NOT perturb this run - 7 seeds give
#    bit-identical trajectories): |slope| mg 1.854-1.986e-3 vs mg3 2.946-3.235e-3, paired
#    difference +1.04e-3 to +1.25e-3 with 8/8 positive and no overlap; mg3 is above the floor in
#    8/8 replicates, mg in 0/8. At dt = 0.25 (5 replicates) both arms stay below the floor.
#    This is a DIFFERENT statistic from the one package 11 settled (that one is the per-step
#    |dEpot| spike tail over 130 cells, where mg/mg2/mg3 are indistinguishable at n = 780/arm);
#    it neither contradicts nor is contradicted by it.
#
# Why no recalibration was applied here: neither of the two routes this file's own history
# offers is available without a judgement call. Re-deriving the floor by the documented rule
# from the mg3 numbers gives ~1.2e-2 Eh/ps, which is the "gates almost nothing" outcome the
# Sep-18 note above warns against. Moving the operating point does not work either - a scan of
# 5000-16000 K (WORK_STATUS package 12) finds no temperature where the event count returns to
# the calibrated band: below ~5200 K the bath produces zero events (MIN_REBUILDS would fail),
# 5300-5500 K is a knife edge (2-2758 rebuilds over 100 K), 6000-12000 K falls smoothly from
# 108/376 to 46/138 rebuilds, and at 13000 K the intramolecular break/re-form channel opens and
# it jumps back to 248/696. The widest margin inside the smooth stretch is at 12000 K
# (5.70e-4 / 1.28e-3, i.e. 4.4x / 2.0x below the floor), but it sits directly under that cliff.
# ---------------------------------------------------------------------------
TEMPERATURE=8000
MAXTIME_FS=10000        # 10 ps, see calibration above
PRINT_FREQUENCY_FS=100  # 0.1 ps -- dense enough for a determined slope fit
FIT_BUCKET_PS=0.1        # de-duplication bucket for the fit, matches print frequency
RATIO_FACTOR=1.5
SLOPE_FLOOR=2.5e-3       # Eh/ps, see the Sep 18, 2026 calibration above
MAX_T_K=10400.0
MIN_REBUILDS=20

run_one() {
    local sub="$1" method="$2" topo="$3" dt="$4"
    rm -rf "$sub"
    mkdir -p "$sub"
    cp input.xyz "$sub/input.xyz"
    # a stale topology cache silently changes the perceived topology (Known
    # Issue #11 in CLAUDE.md) -- always start from a clean one
    ( cd "$sub" && timeout 260 $CURCUMA -md input.xyz -method "$method" -gfnff.topology_mode "$topo" \
        -temperature "$TEMPERATURE" -maxtime "$MAXTIME_FS" -md.time_step "$dt" \
        -md.print_frequency "$PRINT_FREQUENCY_FS" \
        -md.wall_type none \
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
    #    at least MIN_REBUILDS reactive rebuilds in total, no non-finite table row,
    #    T_max < MAX_T_K.
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$RATIO_FACTOR" "$SLOPE_FLOOR" "$MAX_T_K" "$FIT_BUCKET_PS" "$MIN_REBUILDS" <<'PYEOF'
import re, sys, math

ratio_factor, slope_floor, max_t, bucket_ps, min_rebuilds = (float(x) for x in sys.argv[1:6])
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
rebuilds = 0
for name in ["gfnff_025", "gfnff_0125", "rev_025", "rev_0125"]:
    txt = open(f"{name}/stdout.log").read()
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
    if name.startswith("rev"):
        rebuilds += len(re.findall(r"REACT rebuild", txt))
        print(f"{name}: {len(re.findall(r'REACT bond formed', txt))} formed, "
              f"{len(re.findall(r'REACT bond broken', txt))} broken, "
              f"{len(re.findall(r'REACT rebuild', txt))} rebuilds")

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
print(f"T_max={t_max_all:.1f} any_nonfinite_row={any_nonfinite} rev_rebuilds={rebuilds}")

ok = True
reasons = []
if not ok_025:
    ok = False; reasons.append("slope(rev, dt=0.25) exceeds threshold")
if not ok_0125:
    ok = False; reasons.append("slope(rev, dt=0.125) exceeds threshold")
if rebuilds < min_rebuilds:
    ok = False
    reasons.append(f"only {rebuilds} reactive rebuilds in the two revgfnff runs "
                   f"(>= {min_rebuilds:.0f} required): the blend is not being exercised")
if t_max_all > max_t:
    ok = False; reasons.append(f"T_max {t_max_all:.1f} > {max_t}")
if any_nonfinite:
    ok = False; reasons.append("non-finite MD table row(s)")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: revgfnff slope within threshold of plain gfnff's slope at both dt, blend exercised, T_max ok"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: revgfnff slope, rebuild count, T_max, or table finiteness outside thresholds"
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
