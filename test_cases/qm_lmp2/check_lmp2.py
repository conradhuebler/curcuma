#!/usr/bin/env python3
"""Native local PNO-RI-MP2 (WP15, method lmp2) checks for ctest -L qm_lmp2.

Claude Generated (Oct 2026). Two modes:

  --exact      complete pair spaces (-qm.lmp2_domain_radius 1e6, -qm.lmp2_pair_cut 1e6, -qm.lmp2_pno_cut 0,
               -qm.lmp2_fock_cut 0, tight convergence) with -qm.lmp2_check true: the iterated LMP2 correlation
               energy must equal curcuma's canonical RI-MP2 with the same auxiliary basis within --tol (1e-9 Eh).
               This verifies the localisation, the PAO spaces, the RI factors and the residual equations.
  --case NAME  default thresholds; E_corr(LMP2) is compared with reference.json (produced by
               scripts/lmp2_reference.py with curcuma itself, i.e. a regression reference) within the tolerance
               stored there; the canonical RI-MP2 value of the json is used to print the local error per atom.

  --extra=...  further flags (one string);  --threads N  -qm.threads (default 4)
"""
import argparse
import json
import os
import re
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz")
ap.add_argument("--ref")
ap.add_argument("--case")
ap.add_argument("--exact", action="store_true")
ap.add_argument("--fc", action="store_true")
ap.add_argument("--extra", default="")
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--tol", type=float, default=1e-9)
a = ap.parse_args()

EXACT = ["-qm.lmp2_domain_radius", "1e6", "-qm.lmp2_pair_cut", "1e6", "-qm.lmp2_pno_cut", "0",
         "-qm.lmp2_fock_cut", "0", "-qm.lmp2_threshold", "1e-12", "-qm.lmp2_residual_threshold", "1e-10",
         "-qm.lmp2_check", "true"]


def run(xyz, flags):
    cmd = [a.curcuma, "-sp", xyz, "-method", "lmp2", "-qm.basis", "def2-SVP", "-qm.scf_threshold", "1e-10",
           "-qm.threads", str(a.threads), "-verbosity", "2", "-no_bmt", "-gpu", "none"] + flags + a.extra.split()
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    v = {}
    for k, p in {"corr": r"E\(LMP2 corr\)\s*=\s*(-?\d+\.\d+) Eh", "check": r"E\(RI-MP2 corr\)\s*=\s*(-?\d+\.\d+) Eh",
                 "hf": r"E\(HF\)\s*=\s*(-?\d+\.\d+) Eh", "total": r"E\(lmp2\)\s*=\s*(-?\d+\.\d+) Eh"}.items():
        m = re.search(p, txt)
        if m:
            v[k] = float(m.group(1))
    for line in txt.splitlines():
        if any(s in line for s in ("pairs kept", "PNOs per kept", "LMP2 iterations", "LMP2 timing", "amplitude fraction")):
            print(line.strip())
    if "corr" not in v:
        print(" ".join(cmd))
        print(txt[-3000:])
        sys.exit("no E(LMP2 corr) in the output")
    return v


if a.exact:
    flags = EXACT + (["-qm.frozen_core", "true"] if a.fc else [])
    v = run(a.xyz, flags)
    if "check" not in v:
        sys.exit("no E(RI-MP2 corr) in the output")
    d = abs(v["corr"] - v["check"])
    print(f"E_corr LMP2 (complete spaces) = {v['corr']:.12f}  RI-MP2 = {v['check']:.12f}  diff = {d:.2e}  tol = {a.tol:.0e}")
    sys.exit(0 if d <= a.tol else 1)

R = json.load(open(a.ref))
c = R["cases"][a.case]
xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", "..", c["xyz"])
v = run(xyz, c["flags"])
tol = c.get("tol", R["tolerance"])
d = abs(v["corr"] - c["e_corr_lmp2"])
print(f"E_corr LMP2 = {v['corr']:.12f}  reference = {c['e_corr_lmp2']:.12f}  diff = {d:.2e}  tol = {tol:.0e}")
err = v["corr"] - c["e_corr_rimp2"]
print(f"(info) LMP2 - canonical RI-MP2 = {err:.3e} Eh = {err / c['natoms']:.3e} Eh/atom")
sys.exit(0 if d <= tol else 1)
