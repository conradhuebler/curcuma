#!/bin/bash

# Test: revgfnff react mode must REFUSE the water-dimer hydrogen bond
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 2026)
#
# Regression test for the react formation criterion introduced by commit
# 4601be27 ("Tighten the react formation criterion, fix the dE_jump NaN at
# rebuild #1"). Before that commit a non-bonded pair joined the topology as
# soon as its WIDE term-weight switch b_ij (rev_bo_*: f_b 2.0, k -7.5) passed
# rev_bo_form = 0.05. That switch is still 0.59 at the water dimer's H...O
# distance (1.952 A) and crosses 0.05 only at 2.31x the covalent sum, i.e. it
# joins ordinary hydrogen bonds and van-der-Waals contacts. Measured on the
# Cs water dimer that meant 18 formations / 15 breaks / 33 rebuilds per ps and
# a mean Epot 2.6-4.0 kcal/mol BELOW the plain-gfnff static run -- an
# unphysical extra attraction created purely by the topology scan.
#
# The default criterion is now rev_form_switch = order, which joins on the
# NARROW bond-order switch b2_ij (rev_bo2_*: f_b 1.4, k -6.0) at
# rev_bo2_form = 0.1. That switch reads 3.9e-4 at the water dimer's H...O and
# 9.7e-7 at its O...O contact (both far below 0.1) and crosses 0.1 only at
# 1.611x the covalent sum -- the same radius as the non-rev react formation
# factor (1.6) and the bond-break radius, so formation and break are
# symmetric. The H bond must therefore be refused.
#
# What this test asserts (standard Cs water dimer, 300 K, 1 ps, CSVR,
# dt 0.25 fs, no wall, one thread):
#   1. no instability in any of the three runs;
#   2. default (order): ZERO "REACT bond formed" and ZERO "REACT rebuild" --
#      the blend must not fire on a hydrogen-bonded pair at all;
#   3. default (order): mean Epot within 0.01 kcal/mol of the plain
#      -method gfnff -gfnff.topology_mode static run of the same trajectory
#      length -- i.e. the react-mode run is energetically identical to the
#      static one, which is the whole point of refusing the join. The
#      trajectory is deterministic (test 16/18 note: -md.seed does not change
#      the initial velocities on this branch), so this is a reproducible
#      comparison, not a sampled one.
#   4. NEGATIVE CONTROL (-gfnff.rev_form_switch weight, the pre-fix
#      behaviour): at least one formation. Without this the test could pass
#      because the system went inert for an unrelated reason (e.g. the scan
#      silently stopped running) rather than because the criterion is right.
#
# Measured on the current binary (Sep 13, 2026):
#   run      formed broken rebuilds  mean Epot (Eh)    vs static
#   static        0      0        0   -0.66123056      --
#   order         0      0        0   -0.66123166      -0.0007 kcal/mol
#   weight        6      3        9   -0.67062858      -5.897  kcal/mol
#
# The weight control run also carries a mean Epot 5.9 kcal/mol below the
# static run (the STATUS log of the fix measured 2.6-4.0 kcal/mol on the same
# molecule with NVE instead of CSVR; the size of the offset depends on how
# far the softened pair drifts, so only its sign is asserted here and the
# number is quoted, not gated).
#
# 1 ps at 300 K is short on purpose: the dimer sits at a van-der-Waals
# minimum and CSVR at 300 K (coupling 10 fs, the default) is what keeps it
# together -- in NVE the two monomers drift apart within the window and the
# H...O distance stops being the thing under test. The run is deterministic,
# so 1 ps is enough to reproduce the numbers above bit-for-bit.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../test_utils.sh"

TEST_NAME="simplemd - 19: revgfnff react join refuses the water-dimer H bond"
TEST_DIR="$SCRIPT_DIR"

TEMPERATURE=300
MAXTIME_FS=1000         # 1 ps
DT=0.25
PRINT_FREQUENCY_FS=100  # 0.1 ps, enough rows for a stable block mean
EPOT_TOL_KCAL=0.01
# RECALIBRATED Sep 18, 2026 (stage 3a; test_cases/revgfnff/_log/WORK_STATUS.md package 5).
# Assertion 3 used to compare the react run against a plain `-method gfnff` static run with a
# 0.01 kcal/mol tolerance. Stage 3a (i)'s r0 fix gives revgfnff a small CONSTANT equilibrium
# offset against gfnff - measured on this very dimer, -0.5555 kcal/mol over 1 ps - so that
# comparison now fails for a reason that has nothing to do with the join. Measured with the
# stage-3a default: revgfnff REACT and revgfnff STATIC have the SAME mean Epot to 0.0000
# kcal/mol (-0.66200626 Eh both), while gfnff static is -0.66112094. The join is therefore still
# exactly free, which is what this test exists to assert - so the reference arm becomes
# `-method revgfnff -gfnff.topology_mode static`, keeping the 0.01 kcal/mol tolerance, and the
# gfnff arm stays as a loose sanity bound on the rev-vs-gfnff equilibrium offset.
#
# UPDATED Sep 19, 2026 (WORK_STATUS package 6): the operator made `-gfnff.rev_well_form mg` the
# default, and an MG well has its own fitted depth D = s |k_b|, so revgfnff's ABSOLUTE energy
# moves away from plain gfnff by design - measured on this dimer -133.5616 kcal/mol, against
# -0.5555 with the Gaussian. That offset is not what this bound is for (relative energies are
# unaffected: the 167-reaction conformer/S66 guard moves 1.0341 -> 1.0439 kcal/mol MAD), so the
# bound now compares the GAUSSIAN rev arm against gfnff, where it keeps its original meaning and
# its original 2.0 kcal/mol threshold. The mg-vs-gauss offset is printed, not gated - it is an
# AI-fitted table value and will move when the table is refitted.
GFNFF_OFFSET_TOL_KCAL=2.0
HARTREE_KCAL=627.5094740631

COMMON_ARGS="-temperature $TEMPERATURE -maxtime $MAXTIME_FS -md.time_step $DT \
    -md.print_frequency $PRINT_FREQUENCY_FS -md.thermostat csvr \
    -md.rattle_12 false -md.seed 42 -md.no_restart -threads 1 -verbosity 1 -no_bmt"

run_one() {
    local sub="$1"; shift
    rm -rf "$sub"
    mkdir -p "$sub"
    cp input.xyz "$sub/input.xyz"
    # a stale topology cache silently changes the perceived topology (Known
    # Issue #11 in CLAUDE.md) -- always start from a clean one
    ( cd "$sub" && timeout 200 $CURCUMA -md input.xyz "$@" $COMMON_ARGS \
        > stdout.log 2> stderr.log )
    return $?
}

run_test() {
    cd "$TEST_DIR"
    rm -rf gstatic static gaussstatic order weight
    cleanup_bmt_dirs
    local rc=0
    run_one gstatic -method gfnff   -gfnff.topology_mode static          || rc=$?
    run_one static -method revgfnff -gfnff.topology_mode static          || rc=$?
    run_one gaussstatic -method revgfnff -gfnff.topology_mode static \
                   -gfnff.rev_well_form gauss                            || rc=$?
    run_one order  -method revgfnff -gfnff.topology_mode react           || rc=$?
    run_one weight -method revgfnff -gfnff.topology_mode react \
                   -gfnff.rev_form_switch weight                         || rc=$?
    return $rc
}

validate_results() {
    local failed=0

    # 1. Numerical stability across all three sub-runs
    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -qiE "Simulation got unstable|NaN/Inf velocity" \
        static/stdout.log static/stderr.log \
        order/stdout.log order/stderr.log \
        weight/stdout.log weight/stderr.log 2>/dev/null; then
        echo -e "${RED}✗ FAIL${NC}: MD reported instability in at least one run"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    else
        echo -e "${GREEN}✓ PASS${NC}: no instability reported in any run"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    fi

    # 2-4. Event counts, the mean-Epot comparison, and the negative control
    TESTS_RUN=$((TESTS_RUN + 1))
    local py_out py_rc
    set +e
    py_out=$(python3 - "$EPOT_TOL_KCAL" "$HARTREE_KCAL" "$GFNFF_OFFSET_TOL_KCAL" <<'PYEOF'
import re, sys

tol_kcal, hartree_kcal, gfnff_tol = (float(x) for x in sys.argv[1:4])
ROW_RE = re.compile(r"\s+\d+\.\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+")

def counts(path):
    txt = open(path, errors="replace").read()
    return (len(re.findall(r"REACT bond formed", txt)),
            len(re.findall(r"REACT bond broken", txt)),
            len(re.findall(r"REACT rebuild", txt)))

def mean_epot(path):
    rows = [l.split() for l in open(path, errors="replace").read().splitlines()
            if ROW_RE.match(l)]
    vals = []
    for r in rows:
        try:
            vals.append(float(r[1]))
        except (IndexError, ValueError):
            pass
    return sum(vals) / len(vals) if vals else float("nan")

ok = True
reasons = []
f_o, b_o, r_o = counts("order/stdout.log")
f_w, b_w, r_w = counts("weight/stdout.log")
f_s, _, _ = counts("static/stdout.log")
f_gs, _, _ = counts("gstatic/stdout.log")

if f_s != 0:
    ok = False; reasons.append(f"revgfnff static formed {f_s} bonds (must be 0)")
if f_gs:
    ok = False; reasons.append(f"plain gfnff static formed {f_gs} bonds (must be 0)")

if f_o != 0 or r_o != 0:
    ok = False
    reasons.append(f"default criterion formed/rebuilds = {f_o}/{r_o} (must be 0/0): "
                   f"the H bond was joined")

e_static = mean_epot("static/stdout.log")
e_gstatic = mean_epot("gstatic/stdout.log")
e_gaussstatic = mean_epot("gaussstatic/stdout.log")
e_order = mean_epot("order/stdout.log")
d_kcal = (e_order - e_static) * hartree_kcal
if not (abs(d_kcal) <= tol_kcal):
    ok = False
    reasons.append(f"mean Epot {d_kcal:+.4f} kcal/mol from the revgfnff static run "
                   f"(tolerance {tol_kcal})")
# loose sanity bound on the rev-vs-gfnff equilibrium offset (stage 3a (i)'s r0 fix): it is
# -0.5555 kcal/mol here and is NOT what this test gates, but a tenfold change would mean the
# equilibrium form moved and should be looked at. Measured on the GAUSSIAN rev arm, so that a
# fitted well depth (the mg / erfmorse tables) cannot trip it - see the header.
d_gfnff = (e_gaussstatic - e_gstatic) * hartree_kcal
if not (abs(d_gfnff) <= gfnff_tol):
    ok = False
    reasons.append(f"revgfnff static is {d_gfnff:+.4f} kcal/mol from plain gfnff static "
                   f"(tolerance {gfnff_tol}) - the equilibrium offset moved")

if f_w < 1:
    ok = False
    reasons.append(f"negative control (rev_form_switch weight) formed {f_w} bonds "
                   f"-- the old criterion is no longer being exercised, so this "
                   f"test no longer proves anything")

print(f"static : formed={f_s} (revgfnff)  gfnff static formed={f_gs}")
print(f"offsets: react-revstatic={d_kcal:+.4f}  revstatic(gauss)-gfnff={d_gfnff:+.4f}  "
      f"revstatic(default well)-gfnff={(e_static - e_gstatic) * hartree_kcal:+.4f} kcal/mol")
print(f"order  : formed={f_o} broken={b_o} rebuilds={r_o}")
print(f"weight : formed={f_w} broken={b_w} rebuilds={r_w} (negative control)")
print(f"mean Epot: static={e_static:.8f} Eh order={e_order:.8f} Eh "
      f"delta={d_kcal:+.4f} kcal/mol (tol {tol_kcal})")
print(f"mean Epot weight vs static = "
      f"{(mean_epot('weight/stdout.log') - e_static) * hartree_kcal:+.3f} kcal/mol")
if not ok:
    print("FAIL reasons: " + "; ".join(reasons))
sys.exit(0 if ok else 1)
PYEOF
)
    py_rc=$?
    set -e
    echo "$py_out"
    if [ $py_rc -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: default criterion refuses the H bond, energy matches static, control fires"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: the H bond was joined, or the energy/control check failed"
        TESTS_FAILED=$((TESTS_FAILED + 1)); failed=1
    fi

    return $failed
}

main() {
    test_header "$TEST_NAME"
    run_test
    assert_exit_code $? 0 "all three sub-runs should complete without crash"
    validate_results
    print_test_summary
    [ $TESTS_FAILED -gt 0 ] && exit 1 || exit 0
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then main "$@"; fi
