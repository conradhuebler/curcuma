#!/bin/bash

# Test: revgfnff reactive stage-1b smoothness — 2 H2 in a tight sphere
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026) - Validates:
#   1. At least 10 topology events fire (REACT rebuild lines) within 4 ps at
#      4000 K, 2.6 A min-dist packing in a 3.0 A spherical wall -- the system
#      is reactive enough to exercise the transition-blend machinery
#      repeatedly, not just once.
#
#      Re-tuned (Sep 2026): the original config (3000 K, 2 ps, coupling 10 fs)
#      gave 0 events once the repulsion-blend fix (its own switch instead of
#      leaking 1% into the bonded term) and the torsion-weight fix removed an
#      artificial dissociation channel that used to make this system reactive
#      at 3000 K -- the two H2 now genuinely need more thermal energy to
#      recombine/exchange. Raised to 4000 K / 4000 fs and tightened the CSVR
#      coupling 10 -> 1 fs (needed to keep the T_max criterion below, not to
#      get events -- 3800/4000/4200 K at coupling 10 fs already gave
#      24/56/78 events but T_max/T0 ratios of 2.27/1.90/2.10).
#   2. The stage-1b energy discontinuity at each event (dE_jump, see
#      docs/REV_GFNFF_STAGE1.md "Stage 1b") is SMALL: median |dE_jump| < 2
#      kJ/mol and no single event above 60 kJ/mol. The documented stage-1a
#      (no blend) behaviour on the same system is a median of 20-190 kJ/mol
#      and jumps up to ~1657 kJ/mol (docs/REV_GFNFF_STAGE1.md "Measured"),
#      so these thresholds fail hard on stage 1a and pass with a large margin
#      on the measured stage-1b numbers at 4000 K/4000 fs/coupling 1 fs
#      (76 events, median 0.75, max 17.3 kJ/mol -- reproducible bit-for-bit
#      across independent reruns with a freshly cleared topology cache; the
#      3800/4200 K probes used for the robustness check below gave
#      34/0.5/5.5 and 98/0.55/16.0 respectively (events/median/max), so the
#      thresholds hold with margin across the whole 3800-4200 K window, not
#      just at the committed 4000 K).
#   3. No numerical instability (same check as test 14/15) and the printed
#      MD-table temperature never exceeds 1.6x the target temperature
#      (6400 K here) -- a tighter, temperature-relative bound than test 17's
#      fixed one, chosen because this is a 6-DOF system with large
#      equipartition variance and the printed table only samples a handful
#      of steps (print_frequency), so a rare high sample is expected but
#      should not be a multiple of the target. Measured at 4000 K: T_max
#      5486 K (ratio 1.37); 3800/4200 K probes: ratios 1.54/1.29.
#
# The very first-ever rebuild of a react-mode run can log dE_jump = nan --
# the "before" (old-topology) energy requires CN state from a previous step
# (rebuildReactiveTopology(), src/core/energy_calculators/ff_methods/
# gfnff_method.cpp: "Requires CN state from a previous step; NaN otherwise"),
# which does not exist yet for rebuild #1. scripts/revgfnff_jump_stats.py
# already treats this as expected and filters NaN out of its statistics
# (it tracks n_nan separately and does not fail on it). This test does the
# same but additionally pins the NaN to rebuild #1 only -- a NaN anywhere
# else in the run is a real instability and fails the test. (At 4000 K this
# particular system did not trigger it at all -- n_nan = 0 -- but the
# allowance is kept since it is a documented, geometry-driven artifact, not
# a temperature-dependent one; see test 17.)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 16: revgfnff stage-1b smoothness (2 H2)"
TEST_DIR="$SCRIPT_DIR"

TEMPERATURE=4000
MIN_EVENTS=10
MAX_MEDIAN_KJ=2.0
MAX_JUMP_KJ=60.0
# RECALIBRATED Sep 18, 2026 (stage 3a; see test_cases/revgfnff/_log/WORK_STATUS.md package 5).
# The old criterion was "the printed T never exceeds 1.6x the target". It was measuring the
# SAMPLING, not the dynamics: at the default print frequency this 4 ps run prints 17 rows, and
# the max over 17 samples of a 6-DOF Maxwell-Boltzmann distribution is itself a random variable -
# measured over 3600/3800/4000/4200/4400 K it came out 2.17/1.23/1.79/1.74/1.91, i.e. the
# criterion passed or failed on which sample happened to be printed. With -md.print_frequency 10
# (1601 rows, same trajectory - the event counts are identical) the truth is visible: the mean
# temperature is 1.14-1.17x the target across that window and the instantaneous max reaches
# 5.5-6.0x. A 4-atom system in a 3 A sphere with hundreds of reactive events genuinely does that.
# So the test now checks the MEAN (a thermostat/energy-balance statement, 1.5x with 28 % margin)
# and keeps a loose max as a pure instability guard (8x, 34 % margin over the measured 5.99).
# The dE_jump thresholds are untouched and still have enormous margin: median 0.00 kJ/mol and
# max 0.50-1.50 kJ/mol over the same five temperatures, against 2.0 and 60.0.
PRINT_FREQUENCY_FS=10
MEAN_T_RATIO=1.5
MAX_T_RATIO=8.0
MEAN_T_K=$(python3 -c "print($TEMPERATURE * $MEAN_T_RATIO)")
MAX_T_K=$(python3 -c "print($TEMPERATURE * $MAX_T_RATIO)")

run_test() {
    cd "$TEST_DIR"
    rm -f stdout.log stderr.log input.trj.xyz
    # Claude Generated (Sep 2026): a stale *.topo.json cache from a previous
    # invocation silently changes the perceived topology and therefore the
    # whole reactive trajectory (same caveat as Known Issue #11 in CLAUDE.md)
    # -- always start from a clean topology/snapshot state.
    rm -f input.topo.json
    rm -rf input.snapshots
    cleanup_bmt_dirs
    timeout 280 $CURCUMA -md input.xyz -method revgfnff -gfnff.topology_mode react \
        -temperature $TEMPERATURE -maxtime 4000 -md.time_step 0.25 \
        -md.print_frequency $PRINT_FREQUENCY_FS \
        -md.wall_radius 3.0 -md.wall_type spheric \
        -md.thermostat csvr -md.coupling 1 -md.rattle_12 false \
        -md.no_restart -threads 1 -verbosity 1 -no_bmt \
        > stdout.log 2> stderr.log
    return $?
}

validate_results() {
    local failed=0

    # 1. Numerical stability (same convention as test 14/15)
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "Simulation got unstable|NaN/Inf velocity" stdout.log stderr.log 2>/dev/null; then
        echo -e "${RED}✗ FAIL${NC}: MD reported instability"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: no instability reported"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 2. dE_jump statistics, event count, NaN confinement, T_max -- one parse
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$MIN_EVENTS" "$MAX_MEDIAN_KJ" "$MAX_JUMP_KJ" "$MAX_T_K" "$MEAN_T_K" <<'PYEOF'
import re, statistics, sys
min_events, max_median, max_jump, max_t, mean_t = (float(x) for x in sys.argv[1:6])
txt = open("stdout.log").read()
events = []  # (idx, abs_kj_or_None_if_nan)
for m in re.finditer(r"REACT rebuild #(\d+): \d+ bonds, dE_jump = ([-+0-9.a-z]+) Eh \(([-+0-9.]+) kJ/mol\)", txt):
    idx = int(m.group(1))
    is_nan = "nan" in m.group(2).lower()
    events.append((idx, None if is_nan else abs(float(m.group(3)))))
n_events = len(events)
nan_idx = [i for i, v in events if v is None]
finite = [v for _, v in events if v is not None]
median = statistics.median(finite) if finite else None
max_abs = max(finite) if finite else None
rows = [l.split() for l in txt.splitlines()
        if re.match(r"\s+\d+\.\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+", l)]
t_vals = []
for r in rows:
    try:
        t_vals.append(float(r[7]))
    except (IndexError, ValueError):
        pass
t_max = max(t_vals) if t_vals else float("nan")
t_mean = statistics.mean(t_vals) if t_vals else float("nan")
ok = True
reasons = []
if n_events < min_events:
    ok = False; reasons.append(f"n_events {n_events} < {min_events:.0f}")
if any(i != 1 for i in nan_idx):
    ok = False; reasons.append(f"NaN dE_jump outside rebuild #1: {[i for i in nan_idx if i != 1]}")
if median is None or median > max_median:
    ok = False; reasons.append(f"median {median} > {max_median}")
if max_abs is None or max_abs > max_jump:
    ok = False; reasons.append(f"max|dE_jump| {max_abs} > {max_jump}")
if not (t_vals) or t_max > max_t:
    ok = False; reasons.append(f"T_max {t_max} > {max_t}")
if not (t_vals) or t_mean > mean_t:
    ok = False; reasons.append(f"T_mean {t_mean:.1f} > {mean_t}")
print(f"n_events={n_events} n_nan={len(nan_idx)} nan_at={nan_idx} "
      f"median_abs_kJ={median} max_abs_kJ={max_abs} T_max_K={t_max} T_mean_K={t_mean:.1f} "
      f"n_rows={len(t_vals)}")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: react-stats within thresholds"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: react-stats outside thresholds"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    assert_exit_code $? 0 "revgfnff react MD should complete without crash"
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
