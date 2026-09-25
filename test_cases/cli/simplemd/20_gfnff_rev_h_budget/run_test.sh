#!/bin/bash

# Test: the hydrogen valence budget of the rev-gfnff share (stage 3a(ii))
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 18, 2026)
#
# Regression test for -gfnff.rev_budget_fix_h, the default since Sep 18, 2026
# (test_cases/revgfnff/_log/{RUNAWAY,HBUDGET}_STATUS.md, FABLE_REVIEW_2 A.1/A.2).
#
# The share's effective valence is Val_i = Val_Z + softplus(settled partners - Val_Z). For a
# HYDROGEN that rule is wrong: a bridging H - a just-formed H2 still bonded to its carbon - is a
# 3c-2e bond with ONE valence, but the softplus lets its budget reach 2 as soon as the second
# partner's tight bond order crosses the settled window, and BOTH of its partial wells jump from
# half share to full share inside one 0.25 fs step with no topology event. The system then runs
# into that artificial minimum, is thrown apart, and the scan finds the pair outside its window.
#
# This is a two-armed test and the CONTROL is the point: the same trajectory with
# -gfnff.rev_budget_fix_h false must VIOLATE both bounds. Without it the test would pass on a
# build where the share is switched off entirely.
#
# System: frame 16 of test_cases/revgfnff/fit_work/c2h6_1000K.xyz - the cell in which the runaway
# was first isolated - 5 ps at 2000 K, CSVR, dt 0.25 fs, seed 42, one thread.
#
# UPDATED Sep 19, 2026 (WORK_STATUS package 6). The operator flipped two defaults in that
# session - `-gfnff.rev_share_form conserving` and `-gfnff.rev_well_form mg` - and EITHER of them
# removes this runaway on its own, so with the new defaults the negative control no longer
# violates anything and the test stopped exercising what it exists for:
#
#   * under `conserving` the hydrogen budget cannot grow at all - the excess-budget cap is
#     X = 0 for H by element (not by the fix_h flag), so -gfnff.rev_budget_fix_h is a NO-OP
#     there. Both arms measured 70.77 kJ/mol and 1.536 a0.
#   * under `mg` + `delivered` the control is 53.38 kJ/mol / 1.839 a0, i.e. also inside the
#     bounds: the deeper, wider MG well removes this particular runaway too.
#
# So the two-armed hydrogen-budget statement is now made in the mode where the budget EXISTS -
# both arms pin `-gfnff.rev_well_form gauss -gfnff.rev_share_form delivered` - and a third arm
# asserts that the SHIPPED DEFAULT is inside the same bounds. No threshold was changed.
#
# Measured (Sep 19, 2026, binary 916847ff; the first two reproduce the Sep 18 numbers exactly):
#
#   arm                                        max per-step dEpot     min r(H-H)
#   gauss+delivered, fix_h on                        59.26 kJ/mol       1.770 a0
#   gauss+delivered, fix_h false (control)         2593.60 kJ/mol       0.559 a0
#   shipped default (mg + conserving)                70.77 kJ/mol       1.536 a0
#
# Thresholds: 150 kJ/mol (2.1x margin on the shipped default, 17x below the control) and 1.0 a0
# (1.54x margin, 1.8x above the control). Two hydrogens at 0.559 a0 = 0.30 Angstrom is not
# chemistry.
#
# The per-step metric EXCLUDES intervals containing a "REACT rebuild" line: a topology change has
# its own dE_jump statistic (tests 16/17) and is not what this test is about. That exclusion is
# also why the metric is meaningful at all - the runaway carries NO topology event, which is
# exactly why the dE_jump statistics of the 22-cell grid never saw it.
#
# ============================================================================================
# KNOWN FAILING since the MD time-step unit fix (Sep 2026) - OPERATOR DECISION PENDING.
# Nothing was changed here; the failure is left visible on purpose. Claude Generated (Sep 2026).
#
# SimpleMD's step is now converted from real femtoseconds into the integrator's own time unit
# sqrt(amu*A^2/Eh) = 1.9516144 fs. Before that, "-md.time_step 0.25" above really integrated
# 0.4879 fs and "-maxtime 5000" really ran 9.76 ps; the CSVR coupling and the COM-removal
# interval are ratios/cadences in the same nominal fs and were stretched with it. The three arms
# therefore now run at a genuinely 0.25 fs step, and the negative control no longer violates:
#
#   arm                                   pre-fix (real 0.4879 fs)     now (real 0.25 fs)
#   gauss+delivered, fix_h on                59.26 kJ / 1.770 a0        25.40 kJ / 2.080 a0
#   gauss+delivered, fix_h false (control) 2593.60 kJ / 0.559 a0        32.68 kJ / 2.120 a0
#   shipped default (mg + conserving)        70.77 kJ / 1.536 a0        31.48 kJ / 1.994 a0
#
# The time-step fix is NOT implicated, proven rather than argued: restoring the OLD physical
# regime on the FIXED binary - every dt-derived setting multiplied by 1.9516144204, i.e.
#   -md.rev_dt_cap 0 -md.time_step 0.4879036051 -maxtime 9758.072102
#   -md.coupling 19.516144204 -md.remove_com_motion 195.16144204
# - reproduces all three arms EXACTLY: 59.26 / 2593.60 / 70.77 kJ and 1.770 / 0.559 / 1.536 a0.
#
# Why this is not a mechanical recalibration. Holding the discrete dynamics fixed (CSVR per-step
# ratio 0.025, COM removal every 400 steps, 20000 steps) and varying only the true step length
# gives, for the control arm (max per-step dEpot / min r(H-H)):
#
#   0.125  17.19/1.757   0.20  19.81/2.082   0.25  32.68/2.120   0.30  34.36/2.144
#   0.35   45.21/1.986   0.40  52.39/2.098   0.45  49.31/1.969   0.4879 2593.60/0.559 VIOLATES
#   0.50   65.99/1.930   0.55  91.20/1.941   0.60  641.62/0.458 VIOLATES
#
# The violation is not a threshold in dt - 0.45 and 0.50 sit on either side of 0.4879 and both
# stay inside - it is a rare event of ONE chaotic trajectory that the step length reshuffles.
# (-md.seed does not help: it does not change the initial velocities in this path, verified over
# 8 seeds, all bit-identical.) Any "recalibrated" step would therefore be a lucky draw, so the
# thresholds and flags were deliberately left untouched.
#
# The substantive question for the operator: the falsifier for -gfnff.rev_budget_fix_h now only
# exists ABOVE the shipped 0.25 fs cap, so the default's justification needs either a new cell /
# temperature that exposes the budget at 0.25 fs, or a re-scoping of what this test asserts.
# Measurement and options: test_cases/revgfnff/_log/WORK_STATUS.md package 10.
# ============================================================================================
#
# RE-SCOPED Sep 25, 2026 (feature/multi-gpu merge, operator decision): the non-exploding control
# is accepted as correct/improved behaviour. The control arm is still run and printed, but no
# longer gates the test; the fixed and shipped arms keep their 150 kJ/mol / 1.0 a0 bounds
# unchanged. Consequence, stated plainly: this test no longer falsifies -gfnff.rev_budget_fix_h.
# Measured on the merged binary (max per-step dEpot / min r(H-H)): fixed 38.61 / 2.222,
# control 25.03 / 2.064, shipped 224.44 / 1.403 -> the SHIPPED arm violates the dEpot bound.
# That violation is an energy-non-conserving event at t = 1.376-1.380 ps right after an H-H
# topology event. It was checked on 12 replicate runs (T = 1975..2030 K in 5 K steps; -md.seed
# does not change the initial velocities), and it is the known rare event of one chaotic
# trajectory, not a merge regression:
# the pre-merge binary violates on 1/12, the merged binary on 2/12, the merged binary with the
# old d4_cn_cache_threshold 0.01 on 3/12. T = 2000 K (this test) happens to be a violating draw
# now. Left FAILING on purpose for the operator: see MULTIGPU_MERGE_STATUS.md, reconciliation.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 20: rev-gfnff hydrogen valence budget"
TEST_DIR="$SCRIPT_DIR"

MAX_STEP_KJ=150.0
MIN_HH_BOHR=1.0

run_one() {
    local sub="$1"; shift
    rm -rf "$sub"
    mkdir -p "$sub"
    cp input.xyz "$sub/input.xyz"
    # a stale topology cache silently changes the perceived topology (Known Issue #11)
    ( cd "$sub" && timeout 280 $CURCUMA -md input.xyz -method revgfnff \
        -gfnff.topology_mode react -temperature 2000 -maxtime 5000 -md.time_step 0.25 \
        -md.thermostat csvr -md.coupling 10 -md.rattle_12 false -md.no_restart -md.seed 42 \
        -threads 1 -verbosity 1 -no_bmt -md.print_frequency 1 -md.dump_frequency 1 "$@" \
        > stdout.log 2> stderr.log )
    return $?
}

# the hydrogen budget only exists in the delivered share, and only the Gaussian well exposes
# this runaway - so the two-armed statement pins both (see the header). The third arm is the
# shipped default, with no flag at all.
GAUSS_DELIVERED=(-gfnff.rev_well_form gauss -gfnff.rev_share_form delivered)

run_test() {
    cd "$TEST_DIR"
    rm -rf fixed free shipped
    cleanup_bmt_dirs
    local rc=0
    run_one fixed "${GAUSS_DELIVERED[@]}" || rc=$?
    run_one free "${GAUSS_DELIVERED[@]}" -gfnff.rev_budget_fix_h false || rc=$?
    run_one shipped || rc=$?
    return $rc
}

validate_results() {
    local failed=0
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$MAX_STEP_KJ" "$MIN_HH_BOHR" <<'PYEOF'
import math, re, sys

max_step_kj, min_hh = (float(x) for x in sys.argv[1:3])
BOHR = 0.529177210903
H2KJ = 2625.4996394799
ANSI = re.compile(r"\x1b\[[0-9;]*m")
ROW = re.compile(r"\s+\d+\.\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+")


def step_max(path):
    """max |Epot(t+dt) - Epot(t)| in kJ/mol, intervals with a REACT rebuild excluded."""
    txt = ANSI.sub("", open(path, errors="replace").read())
    rows, pending = [], False
    for line in txt.splitlines():
        if "REACT rebuild" in line:
            pending = True
            continue
        if not ROW.match(line):
            continue
        f = line.split()
        try:
            ep = float(f[1])
        except (IndexError, ValueError):
            continue
        rows.append((ep, pending))
        pending = False
    out = 0.0
    for k in range(1, len(rows)):
        if rows[k][1]:
            continue
        a, b = rows[k - 1][0], rows[k][0]
        if a != a or b != b:          # NaN
            continue
        out = max(out, abs(b - a) * H2KJ)
    return out, len(rows)


def min_hh_bohr(path):
    lines = open(path).read().splitlines()
    i, best = 0, float("inf")
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        atoms = [l.split() for l in lines[i + 2:i + 2 + n]]
        i += n + 2
        h = [[float(x) for x in a[1:4]] for a in atoms if a[0] == "H"]
        for a in range(len(h)):
            for b in range(a + 1, len(h)):
                d = math.dist(h[a], h[b])
                if d < best:
                    best = d
    return best / BOHR


res = {}
for arm in ("fixed", "free", "shipped"):
    s, nrows = step_max(f"{arm}/stdout.log")
    hh = min_hh_bohr(f"{arm}/input.snapshots/input.trj.xyz")
    res[arm] = (s, hh, nrows)
    print(f"{arm:6s}: max per-step dEpot {s:9.2f} kJ/mol  min r(H-H) {hh:7.3f} a0  rows {nrows}")

ok, reasons = True, []
for arm, label in (("fixed", "gauss+delivered, fix_h on"), ("shipped", "shipped default")):
    s, hh, nrows = res[arm]
    if nrows < 1000:
        ok = False
        reasons.append(f"only {nrows} status rows parsed for the {label} arm")
    if s > max_step_kj:
        ok = False
        reasons.append(f"{label}: max per-step dEpot {s:.2f} > {max_step_kj} kJ/mol")
    if hh < min_hh:
        ok = False
        reasons.append(f"{label}: min r(H-H) {hh:.3f} < {min_hh} a0")
sf, hf, _ = res["free"]
# Sep 25, 2026 (operator decision, multi-gpu merge): the control no longer explodes and that is
# accepted as the correct/improved behaviour, so it is reported but no longer gates the test.
if sf <= max_step_kj and hf >= min_hh:
    print("INFO: negative control (gauss + delivered + -gfnff.rev_budget_fix_h false) stays inside "
          "both bounds (accepted Sep 25, 2026; no longer gating)")
else:
    print("INFO: negative control violates a bound (the pre-Sep-2026 behaviour)")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: hydrogen budget bounded in both the pinned and the shipped arm (control reported, not gating)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: hydrogen budget check"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi
    return $failed
}

main() {
    test_header "$TEST_NAME"
    cd "$TEST_DIR"
    run_test || true
    validate_results || true
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
