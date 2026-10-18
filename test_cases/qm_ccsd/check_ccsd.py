#!/usr/bin/env python3
"""Native canonical RCCSD (WP13) against the PySCF reference in reference.json.

Claude Generated (Oct 2026). Runs `curcuma -sp <xyz> -method ccsd -verbosity 2` and compares E(HF),
the iteration-0 (MP2) correlation energy, the CCSD correlation energy, the total energy and the T1
diagnostic with the PySCF values (scripts/ccsd_reference.py: PySCF RHF conv_tol 1e-12 + cc.CCSD
conv_tol 1e-12, def2-SVP spherical, frozen=None or the chemical frozen core).

  --fc          adds -qm.frozen_core true and compares with the frozen-core reference
  --extra=...   further flags (one string)
  --threads N   -qm.threads (default 4)
Passes when every energy agrees within --tol (default 1e-8 Eh; E(HF) within --tol-hf) and the T1
diagnostic within 1e-6. The ORCA values stored in reference.json are printed for information.
"""
import argparse
import json
import os
import re
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--mol", required=True)
ap.add_argument("--fc", action="store_true")
ap.add_argument("--extra", default="")
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--tol", type=float, default=1e-8)
ap.add_argument("--tol-hf", type=float, default=1e-8)
ap.add_argument("--tol-t1", type=float, default=1e-6)
a = ap.parse_args()

R = json.load(open(a.ref))["cases"][a.mol]
ref = R["fc" if a.fc else "full"]
xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", R["xyz"])
PAT = {
    "hf": r"E\(HF\)\s*=\s*(-?\d+\.\d+) Eh",
    "mp2": r"E\(MP2 corr\)\s*=\s*(-?\d+\.\d+) Eh",
    "ccsd": r"E\(CCSD corr\)\s*=\s*(-?\d+\.\d+) Eh",
    "t1": r"T1 diagnostic\s*=\s*(\d+\.\d+)",
    "total": r"E\(ccsd\)\s*=\s*(-?\d+\.\d+) Eh",
    "iterations": r"CCSD iterations = (\d+)",
}
cmd = [a.curcuma, "-sp", xyz, "-method", "ccsd", "-qm.basis", "def2-SVP", "-qm.scf_threshold", "1e-10",
       "-qm.threads", str(a.threads), "-verbosity", "2", "-no_bmt", "-gpu", "none"]
if a.fc:
    cmd += ["-qm.frozen_core", "true"]
cmd += a.extra.split()
out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000)
txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
v = {}
for k, p in PAT.items():
    m = re.search(p, txt)
    if not m:
        print(" ".join(cmd))
        print(txt[-3000:])
        sys.exit(f"no '{k}' value in the output")
    v[k] = float(m.group(1))
for line in txt.splitlines():
    if "CCSD timing" in line or "canonical RCCSD" in line:
        print(line.strip())

checks = [
    ("E(HF)", v["hf"], R["hf"], a.tol_hf),
    ("E_corr(MP2, it 0)", v["mp2"], ref["e_mp2"], a.tol),
    ("E_corr(CCSD)", v["ccsd"], ref["e_ccsd"], a.tol),
    ("E_total(ccsd)", v["total"], R["hf"] + ref["e_ccsd"], a.tol + a.tol_hf),
    ("T1 diagnostic", v["t1"], ref["t1_diagnostic"], a.tol_t1),
]
ok = True
for name, got, want, tol in checks:
    d = abs(got - want)
    good = d <= tol
    ok &= good
    print(f"{name:18s} curcuma={got:.12f} pyscf={want:.12f} diff={d:.2e} tol={tol:.0e} {'ok' if good else 'FAIL'}")
print(f"iterations: curcuma {int(v['iterations'])}, PySCF {ref.get('pyscf_iterations')}")
orca = R.get("orca", {}).get("fc" if a.fc else "full")
if orca:
    print(f"(info) ORCA 6.1: E_corr {orca['e_ccsd']:.12f}, curcuma - ORCA = {v['ccsd'] - orca['e_ccsd']:.2e}, "
          f"E(HF) curcuma - ORCA = {v['hf'] - orca['hf']:.2e}")
sys.exit(0 if ok else 1)
