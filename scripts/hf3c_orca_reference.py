#!/usr/bin/env python3
"""ORCA 6.1 references for the native HF-3c beyond Ne (test_cases/qm_hf3c/reference_orca_z.json, ctest -L qm_hf3c_z).

Every case is one closed-shell structure of test_cases/structures, run with ORCA `! HF-3c VeryTightSCF` (gradient
cases add `EnGrad`). ORCA's `! HF-3c` is HF/MINIX + D3(BJ, HF/MINIX parameters, two-body) + gCP(HF/MINIX) + the
short-range basis term; it prints each part, and this script stores them:

  hf       'Total Energy :' of the SCF block (E(HF/MINIX), 14 decimals)
  d3       'Edisp/kcal,au:' second number (12 decimals)
  gcp_srb  'gCP+bas correction' (gCP + SRB, 9 decimals; otool_gcp computes it)
  total    'FINAL SINGLE POINT ENERGY' (12 decimals)
  gradient job.engrad (Eh/Bohr), gradient cases only

Geometry: the library xyz in Angstrom as written (ORCA converts with its own constant; curcuma uses CODATA 2018).

Optimisation references (OPT, --only opt_<name>): `! HF-3c Opt TightOpt VeryTightSCF` from the library start.

Usage: python3 scripts/hf3c_orca_reference.py --orca /path/to/orca [--only CASE|opt_CASE ...] [--work DIR]
Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import json
import os
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.normpath(os.path.join(HERE, "..", "test_cases"))
OUT = os.path.join(TC, "qm_hf3c", "reference_orca_z.json")

# case -> (library path below test_cases/structures, gradient?)
CASES = {
    "NaF": ("small/sodium-fluoride.xyz", True),
    "NaCl": ("small/sodium-chloride.xyz", False),
    "H2S": ("small/hydrogen-sulfide.d-shell-validation.xyz", True),
    "PH3": ("small/phosphine.d-shell-validation.xyz", False),
    "SiH4": ("small/silane.d-shell-validation.xyz", False),
    "HCl": ("small/hydrogen-chloride.xyz", False),
    "SO2": ("small/sulfur-dioxide.xyz", True),
    "AlCl3": ("small/aluminium-trichloride.xyz", False),
    "KCl": ("small/potassium-chloride.xyz", True),
    "ZnCl2": ("metals/zinc-dichloride.xyz", True),
    "TiCl4": ("metals/titanium-tetrachloride.xyz", False),
    "GeH4": ("small/germane.xyz", False),
    "HBr": ("small/hydrogen-bromide.xyz", True),
    "Br2": ("small/dibromine.xyz", False),
    "HI": ("small/hydrogen-iodide.xyz", True),
    "CH3I": ("small/iodomethane.xyz", True),
    "XeF2": ("small/xenon-difluoride.xyz", False),
    "CsI": ("small/caesium-iodide.xyz", False),
    "HgCl2": ("metals/mercury-dichloride.xyz", False),
    "PbH4": ("small/plumbane.xyz", True),
    "CeF4": ("metals/cerium-tetrafluoride.xyz", False),
}

# opt case -> library start structure: ORCA `! HF-3c Opt TightOpt VeryTightSCF` from the same start, final energy
OPT = {"H2S": "small/hydrogen-sulfide.d-shell-validation.xyz", "HI": "small/hydrogen-iodide.xyz"}


def run_orca_opt(orca, work, name, path):
    with open(os.path.join(TC, "structures", path)) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    kw = "HF-3c Opt TightOpt VeryTightSCF"
    d = os.path.join(work, "opt_" + name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "job.inp"), "w") as f:
        f.write("! %s\n* xyz 0 1\n%s\n*\n" % (kw, "\n".join(lines[2:2 + n])))
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(orca))
    out = subprocess.run([orca, "job.inp"], cwd=d, env=env, capture_output=True, text=True).stdout
    if "ORCA TERMINATED NORMALLY" not in out or "OPTIMIZATION RUN DONE" not in out:
        raise RuntimeError("ORCA opt failed in " + d)
    final = open(os.path.join(d, "job.xyz")).read().splitlines()[2:2 + n]
    return {"xyz": "structures/" + path, "orca_input": kw,
            "total": float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)[-1]),
            "orca_final_geometry": final}


def run_orca(orca, work, name, path, gradient):
    with open(os.path.join(TC, "structures", path)) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    kw = "HF-3c VeryTightSCF" + (" EnGrad" if gradient else "")
    inp = "! %s\n* xyz 0 1\n%s\n*\n" % (kw, "\n".join(lines[2:2 + n]))
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "job.inp"), "w") as f:
        f.write(inp)
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(orca))
    out = subprocess.run([orca, "job.inp"], cwd=d, env=env, capture_output=True, text=True).stdout
    if "ORCA TERMINATED NORMALLY" not in out:
        raise RuntimeError("ORCA failed in " + d)
    r = {"xyz": "structures/" + path, "orca_input": kw}
    r["hf"] = float(re.findall(r"Total Energy\s+:\s+(-?\d+\.\d+) Eh", out)[-1])
    r["d3"] = float(re.findall(r"Edisp/kcal,au:\s+\S+\s+(-?\d+\.\d+)", out)[-1])
    r["gcp_srb"] = float(re.findall(r"gCP\+bas correction\s+(-?\d+\.\d+)", out)[-1])
    r["total"] = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)[-1])
    if gradient:
        g = open(os.path.join(d, "job.engrad")).read().split("# The current gradient in Eh/bohr")[1]
        vals = [float(v) for v in g.split("#")[1].split()]
        r["gradient"] = [vals[i:i + 3] for i in range(0, len(vals), 3)]
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", required=True)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--work")
    a = ap.parse_args()
    work = a.work or tempfile.mkdtemp(prefix="hf3c_orca_")
    old = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res = {"_about": "ORCA 6.1.0 `! HF-3c VeryTightSCF` references, scripts/hf3c_orca_reference.py",
           "_units": "Eh, gradient Eh/Bohr", "cases": old.get("cases", {})}
    for name, (path, grad) in CASES.items():
        if a.only and name not in a.only:
            continue
        res["cases"][name] = run_orca(a.orca, work, name, path, grad)
        c = res["cases"][name]
        print("%-6s HF %.12f  D3 %.12f  gCP+SRB %.9f  total %.12f%s" % (
            name, c["hf"], c["d3"], c["gcp_srb"], c["total"], "  +gradient" if grad else ""))
    res["opt"] = old.get("opt", {})
    for name, path in OPT.items():
        if a.only and ("opt_" + name) not in a.only:
            continue
        res["opt"][name] = run_orca_opt(a.orca, work, name, path)
        print("opt %-6s total %.12f" % (name, res["opt"][name]["total"]))
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
