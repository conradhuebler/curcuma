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
# Measured (Sep 18, 2026, binary 7fef81bd):
#
#   arm                          max per-step dEpot     min r(H-H)
#   default (fix_h on)                 59.26 kJ/mol       1.770 a0
#   control (-...fix_h false)        2593.60 kJ/mol       0.559 a0
#
# Thresholds: 150 kJ/mol (2.5x margin on the default, 17x below the control) and 1.0 a0 (1.77x
# margin, 1.8x above the control). Two hydrogens at 0.559 a0 = 0.30 Angstrom is not chemistry.
#
# The per-step metric EXCLUDES intervals containing a "REACT rebuild" line: a topology change has
# its own dE_jump statistic (tests 16/17) and is not what this test is about. That exclusion is
# also why the metric is meaningful at all - the runaway carries NO topology event, which is
# exactly why the dE_jump statistics of the 22-cell grid never saw it.

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

run_test() {
    cd "$TEST_DIR"
    rm -rf fixed free
    cleanup_bmt_dirs
    local rc=0
    run_one fixed || rc=$?
    run_one free -gfnff.rev_budget_fix_h false || rc=$?
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
for arm in ("fixed", "free"):
    s, nrows = step_max(f"{arm}/stdout.log")
    hh = min_hh_bohr(f"{arm}/input.snapshots/input.trj.xyz")
    res[arm] = (s, hh, nrows)
    print(f"{arm:6s}: max per-step dEpot {s:9.2f} kJ/mol  min r(H-H) {hh:7.3f} a0  rows {nrows}")

ok, reasons = True, []
s, hh, nrows = res["fixed"]
if nrows < 1000:
    ok = False
    reasons.append(f"only {nrows} status rows parsed for the default arm")
if s > max_step_kj:
    ok = False
    reasons.append(f"default arm: max per-step dEpot {s:.2f} > {max_step_kj} kJ/mol")
if hh < min_hh:
    ok = False
    reasons.append(f"default arm: min r(H-H) {hh:.3f} < {min_hh} a0")
sf, hf, _ = res["free"]
if sf <= max_step_kj and hf >= min_hh:
    ok = False
    reasons.append("negative control (-gfnff.rev_budget_fix_h false) stayed inside BOTH bounds "
                   "- the test is no longer exercising the hydrogen budget")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: hydrogen budget bounded, control violates both bounds"
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
