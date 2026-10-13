#!/bin/bash

# Test: -gfnff.frag_charge_model ensemble with many fragments - bounded memory, index-free carrier,
#       equality with the former full enumeration, variant memory budget (docs/FRAG_CHARGE_MODEL.md)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026) - AI-generated, machine-tested; human production testing pending
#
# Until Oct 2026 the ensemble enumerated all C(|q|+G-1, |q|) charge distributions over the G
# fragment groups before filtering: 500 waters at charge +3 = 2.1e7 vectors of length 500
# (~40 GB). The selection is now a bounded dynamic programme (GFNFF::fragSelectPlacements).
# Systems are built at run time from library structures: W500 = the 500 waters of
# peo201-water500 (atoms 1411-2910, 1500 atoms, 500 fragments).
#
#   1. MEMORY   - W500 at charge +3 (default model): exit 0, peak RSS < 4 GB (os.wait4 ru_maxrss).
#   2. ORDER    - W500 + H3O+(H2O)2 (shifted 30 A away), charge +1: the carrier is unique (only the
#                 hydronium keeps every fragment closed-shell). Default energy with the cluster
#                 first equals the default energy with the cluster last (|dE| < 1e-6 kcal/mol) and
#                 the reference-rule energy with the hydronium as fragment 0; the reference rule
#                 with a water first must differ by > 50 kcal/mol (liveness).
#   3. BRUTE    - default vs -gfnff.frag_charge_enumerate true (the former full enumeration):
#                 byte-identical energy + gradient dumps for water8 at charge +3 and -2 (8 groups,
#                 cap of 6 placements active) and W500 + H3O+(H2O)2 at charge +2 (126756 distributions).
#   4. BUDGET   - W500 at charge +3 with -gfnff.frag_charge_variant_memory_mb 1: one variant only,
#                 exit 0, the warning printed exactly once, result byte-identical to
#                 -gfnff.frag_charge_max_placements 1, and a lower peak RSS than check 1.
#
# Each run gets RLIMIT_AS = 8 GB (the runs here need < 1.5 GB RSS) and a 300 s timeout, so that a
# regression to the full enumeration fails the test with an allocation error instead of exhausting
# the machine (no systemd dependence). Checked: the binary before the fix exits 1 after 3.6 s there.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 08: frag_charge_model ensemble with many fragments"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import os, resource, subprocess, sys, tempfile, shutil, threading, time
B = os.environ["CURCUMA"]
K = 627.509474
REF = ["-gfnff.frag_charge_model", "reference"]
fails = 0

def limit():
    resource.setrlimit(resource.RLIMIT_AS, (8 << 30, 8 << 30))

def sp(frame, charge, extra, verbosity=0):
    """One single point in a fresh directory; returns (exit code, energy, gradient text, stdout, peak RSS MB, wall s)."""
    d = tempfile.mkdtemp()
    try:
        with open(os.path.join(d, "s.xyz"), "w") as f:
            f.write(f"{len(frame)}\n\n" + "".join(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in frame))
        cmd = [B, "-sp", "s.xyz", "-method", "gfnff", "-gfnff.cache_topology", "false", "-charge", str(charge),
               "-threads", "4", "-verbosity", str(verbosity), "-no_bmt", "-dump_gradient", "g.txt"] + extra
        t0 = time.time()
        with open(os.path.join(d, "out.log"), "w") as o:
            p = subprocess.Popen(cmd, cwd=d, stdout=o, stderr=subprocess.STDOUT, preexec_fn=limit)
            timer = threading.Timer(300.0, p.kill)
            timer.start()
            _, status, ru = os.wait4(p.pid, 0)
            timer.cancel()
        dt = time.time() - t0
        rc = os.waitstatus_to_exitcode(status)
        out = open(os.path.join(d, "out.log"), errors="replace").read()
        g, e = None, None
        if os.path.exists(os.path.join(d, "g.txt")):
            g = open(os.path.join(d, "g.txt")).read()
            for l in g.splitlines():
                if l.startswith("# energy"):
                    e = float(l.split()[2])
        return rc, e, g, out, ru.ru_maxrss / 1024.0, dt
    finally:
        shutil.rmtree(d, ignore_errors=True)

def check(name, ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + f"{name}: {msg}", flush=True)
    if not ok: fails += 1

def rd(p):
    L = open(p).read().split("\n"); n = int(L[0])
    return [(l.split()[0], *map(float, l.split()[1:4])) for l in L[2:2 + n]]

w500 = rd("peo_water.xyz")[1410:]
assert len(w500) == 1500 and [a[0] for a in w500[:3]] == ["O", "H", "H"]
h3o = rd("h3op_h2o2.xyz")
zmax = max(a[3] for a in w500)
h3o_far = [(a[0], a[1], a[2], a[3] + zmax + 30.0) for a in h3o]
water8 = rd("water8.xyz")

# 1. memory
rc, e1, g1, out1, mb1, t1 = sp(w500, 3, [])
check("memory: W500 charge +3 finishes", rc == 0 and e1 is not None, f"exit {rc}, E = {e1}, {t1:.1f} s")
check("memory: peak RSS < 4 GB", mb1 < 4096.0, f"{mb1:.0f} MB")

# 2. order independence of the (unique) carrier
A = h3o_far + w500                       # cluster first (its first atom is a WATER oxygen)
Bo = w500 + h3o_far[6:10] + h3o_far[0:6]  # cluster last, hydronium last
Ah = h3o_far[6:10] + h3o_far[0:6] + w500  # hydronium first: the reference rule picks it
eA = sp(A, 1, [])[1]; eB = sp(Bo, 1, [])[1]; eAh_ref = sp(Ah, 1, REF)[1]; eB_ref = sp(Bo, 1, REF)[1]
ok = None not in (eA, eB, eAh_ref, eB_ref)
check("order: default, cluster first == cluster last", ok and abs(eA - eB) * K < 1e-6, f"{(eA - eB) * K:+.2e} kcal/mol" if ok else "missing energy")
check("order: default == reference with hydronium first", ok and abs(eA - eAh_ref) * K < 1e-6, f"{(eA - eAh_ref) * K:+.2e} kcal/mol" if ok else "missing energy")
check("order: reference with a water first differs (liveness)", ok and (eB_ref - eB) * K > 50.0, f"{(eB_ref - eB) * K:+.1f} kcal/mol" if ok else "missing energy")

# 3. equality with the former full enumeration
for name, frame, q in (("water8 +3", water8, 3), ("water8 -2", water8, -2), ("W500+H3O+(H2O)2 +2", A, 2)):
    rs, es, gs = sp(frame, q, [])[:3]
    re_, ee, ge = sp(frame, q, ["-gfnff.frag_charge_enumerate", "true"])[:3]
    check(f"brute force: {name}", rs == 0 and re_ == 0 and gs is not None and gs == ge,
          f"E = {es} / {ee}, dumps {'identical' if gs == ge else 'DIFFER'}")

# 4. variant memory budget
rc4, e4, g4, out4, mb4, t4 = sp(w500, 3, ["-gfnff.frag_charge_variant_memory_mb", "1"], verbosity=1)
nwarn = out4.count("exceed frag_charge_variant_memory_mb")
g1p = sp(w500, 3, ["-gfnff.frag_charge_max_placements", "1"])[2]
check("budget: runs and warns once", rc4 == 0 and nwarn == 1, f"exit {rc4}, {nwarn} warning(s)")
check("budget: equals max_placements 1", g4 is not None and g4 == g1p, "dumps identical" if g4 == g1p else "dumps DIFFER")
check("budget: lower peak RSS than the default", mb4 < mb1, f"{mb4:.0f} MB vs {mb1:.0f} MB")
sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
