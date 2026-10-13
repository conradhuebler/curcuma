#!/usr/bin/env python3
"""Count the structures of a reference set whose SCF does not converge.

    python scripts/scan_convergence.py release/curcuma gfn2 test_cases/GMTKN55-testset

Used to check whether an outlier in refset_regression.py is a real change or just
a non-converged SCF wandering (GMTKN55/gfn2 has exactly three such structures:
G21IP/b+, be+ and c+, all single-atom cations). Claude Generated (Sep 2026).
"""
import concurrent.futures as cf, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path
binary, method = sys.argv[1], sys.argv[2]
root = Path(sys.argv[3])
xyzs = sorted(root.rglob("struc.xyz")) or sorted(root.rglob("*.xyz"))
def cs(x):
    c = s = 0
    for n in (".CHRG", ".UHF"):
        f = x.parent / n
        if f.exists():
            try: v = int(f.read_text().split()[0])
            except Exception: v = 0
            if n == ".CHRG": c = v
            else: s = v
    return c, s
def run(x):
    c, s = cs(x)
    with tempfile.TemporaryDirectory() as td:
        loc = Path(td) / "struc.xyz"; shutil.copy(x, loc)
        cmd = [binary, "-sp", str(loc), "-method", method, "-threads", "1", "-no_bmt", "-verbosity", "2"]
        if c: cmd += ["-charge", str(c)]
        if s: cmd += ["-spin", str(s)]
        try: out = subprocess.run(cmd, capture_output=True, text=True, timeout=900).stdout
        except subprocess.TimeoutExpired: return (x, "timeout")
    return (x, "notconv" if "NOT converged" in out else "ok")
bad = []
with cf.ThreadPoolExecutor(max_workers=32) as ex:
    for x, st in ex.map(run, xyzs):
        if st != "ok": bad.append((str(x.relative_to(root)), st))
print(f"{method}: {len(xyzs)} structures, {len(bad)} not converged/timeout")
for n, st in bad[:20]: print("   ", st, n)
