#!/bin/bash

# Test: revgfnff reactive stage-1b smoothness — N2 + 3 H2, no refractory period
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026) - Validates the stage-1b transition blend on a
# larger, chemically richer reactive system than test 16 (N-H and H-H bond
# formation/breaking in the same run). The revgfnff preset already ships with
# refractory = 0 and no valence cap (docs/REV_GFNFF_STAGE1.md), so this test
# exercises that default directly rather than overriding it.
#
# Pass criteria (docs/REV_GFNFF_STAGE1.md "Measured", N2 + 3 H2 row):
#   1. No numerical instability (same convention as test 14/15).
#   2. T_max (MD table, column 8) stays below 8000 K -- measured 3.97 kK,
#      well inside the documented 3.4-4.8 kK range for this system/temperature.
#   3. Median |dE_jump| < 2 kJ/mol -- stage 1a (no blend) on this system is
#      documented at median 192 kJ/mol, max 1020 kJ/mol, so this threshold
#      fails hard on stage 1a and passes with a large margin on the measured
#      stage-1b value (median 0.0 kJ/mol, reproduced across repeated runs
#      with a freshly cleared topology cache).
#   4. At least 90% of events have |dE_jump| < 5 kJ/mol -- documented stage-1b
#      share is 0.99; this run measured 1.0 (212 events).
#
# Same NaN-at-rebuild-#1 caveat as test 16: the very first ever topology
# rebuild of a react-mode run has no prior CN state to compute an old-topology
# "before" energy from (rebuildReactiveTopology(), src/core/energy_calculators/
# ff_methods/gfnff_method.cpp), so its dE_jump is nan by construction -- this
# is confirmed reproducible on this system's fixed t=0 geometry (an N-H pair
# sits right at the formation threshold from the initial packing, independent
# of the thermostat seed) and is exactly what scripts/revgfnff_jump_stats.py
# already filters out of its own statistics (n_nan, tracked separately, not a
# failure). This test allows a NaN only at rebuild #1; any later NaN fails.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 17: revgfnff stage-1b smoothness (N2 + 3 H2, no refractory)"
TEST_DIR="$SCRIPT_DIR"

MAX_T_K=8000.0
MAX_MEDIAN_KJ=2.0
MIN_FRAC_LT5=0.9

run_test() {
    cd "$TEST_DIR"
    rm -f stdout.log stderr.log input.trj.xyz
    # Claude Generated (Sep 2026): same cache caveat as test 16 -- a stale
    # *.topo.json changes the perceived topology and the whole trajectory.
    rm -f input.topo.json
    rm -rf input.snapshots
    cleanup_bmt_dirs
    timeout 280 $CURCUMA -md input.xyz -method revgfnff -gfnff.topology_mode react \
        -temperature 3500 -maxtime 3000 -md.time_step 0.25 \
        -md.wall_radius 3.5 -md.wall_type spheric \
        -md.thermostat csvr -md.coupling 10 -md.rattle_12 false \
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

    # 2. dE_jump statistics, NaN confinement to rebuild #1, T_max -- one parse
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$MAX_T_K" "$MAX_MEDIAN_KJ" "$MIN_FRAC_LT5" <<'PYEOF'
import re, statistics, sys
max_t, max_median, min_frac_lt5 = (float(x) for x in sys.argv[1:4])
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
frac_lt5 = (sum(v < 5.0 for v in finite) / len(finite)) if finite else None
rows = [l.split() for l in txt.splitlines()
        if re.match(r"\s+\d+\.\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+", l)]
t_vals = []
for r in rows:
    try:
        t_vals.append(float(r[7]))
    except (IndexError, ValueError):
        pass
t_max = max(t_vals) if t_vals else float("nan")
ok = True
reasons = []
if n_events == 0:
    ok = False; reasons.append("no REACT rebuild events at all")
if any(i != 1 for i in nan_idx):
    ok = False; reasons.append(f"NaN dE_jump outside rebuild #1: {[i for i in nan_idx if i != 1]}")
if median is None or median > max_median:
    ok = False; reasons.append(f"median {median} > {max_median}")
if frac_lt5 is None or frac_lt5 < min_frac_lt5:
    ok = False; reasons.append(f"frac_lt5 {frac_lt5} < {min_frac_lt5}")
if not t_vals or t_max > max_t:
    ok = False; reasons.append(f"T_max {t_max} > {max_t}")
print(f"n_events={n_events} n_nan={len(nan_idx)} nan_at={nan_idx} "
      f"median_abs_kJ={median} frac_lt5={frac_lt5} T_max_K={t_max}")
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
