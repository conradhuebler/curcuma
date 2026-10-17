#!/usr/bin/env python3
"""-opt -method mp2 must reach the PySCF MP2/def2-SVP minimum.

Claude Generated (Oct 2026). Runs `curcuma -opt <start> -method mp2` (tight convergence) in a temporary
directory, re-evaluates the optimised structure with `-sp`, and compares the energy and all interatomic
distances with reference_grad.json["opt"] (scripts/mp2_opt_reference.py: PySCF RHF + mp.MP2 analytic gradient,
scipy BFGS gtol 1e-7 Eh/Bohr). Plain python (no numpy).
Usage: check_mp2_opt.py --curcuma <bin> --ref reference_grad.json --mol H2O|formaldehyde
"""
import argparse, json, math, os, re, shutil, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--mol", required=True)
ap.add_argument("--tol-e", type=float, default=1e-8)
ap.add_argument("--tol-r", type=float, default=1e-4)
a = ap.parse_args()
R = json.load(open(a.ref))["opt"]["cases"][a.mol]
start = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", R["start"])
cur = os.path.abspath(a.curcuma)


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [list(map(float, l.split()[1:4])) for l in L[2:2 + n]]


with tempfile.TemporaryDirectory() as tmp:
    shutil.copy(start, os.path.join(tmp, "start.xyz"))
    p = subprocess.run([cur, "-opt", "start.xyz", "-method", "mp2", "-no_bmt", "-verbosity", "1", "-qm.threads", "4",
                        "-qm.scf_threshold", "1e-9", "-opt.convergence_preset", "tight"],
                       cwd=tmp, capture_output=True, text=True, timeout=2500)
    if p.returncode != 0 or not os.path.exists(os.path.join(tmp, "start.opt.xyz")):
        print(p.stdout[-3000:], p.stderr[-2000:])
        sys.exit("-opt -method mp2 failed")
    out = subprocess.run([cur, "-sp", "start.opt.xyz", "-method", "mp2", "-no_bmt", "-verbosity", "1",
                          "-qm.scf_threshold", "1e-10", "-qm.threads", "4"], cwd=tmp, capture_output=True, text=True,
                         timeout=600, check=True).stdout
    e = float(re.search(r"E\(mp2\)\s*=\s*(-?\d+\.\d+)", re.sub(r"\x1b\[[0-9;]*m", "", out)).group(1))
    x = read_xyz(os.path.join(tmp, "start.opt.xyz"))
d = [math.dist(x[i], x[j]) for i in range(len(x)) for j in range(i + 1, len(x))]
dr = max(abs(u - w) for u, w in zip(d, R["distances"]))
de = e - R["e_min"]
ok = abs(de) <= a.tol_e and dr <= a.tol_r
print(f"{a.mol}: E(opt) {e:.10f}  E(PySCF min) {R['e_min']:.10f}  diff {de:+.2e}  max|d r_ij| {dr:.2e} A  "
      f"{'ok' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
