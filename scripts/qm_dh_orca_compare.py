#!/usr/bin/env python3
"""Native double hybrids (WP12 stage 1) against ORCA 6.1 on matched settings; writes test_cases/qm_dh/reference_orca.json.

ORCA input per case: '! <keyword> def2-SVP NoRI VeryTightSCF DefGrid3', AuxC def2/JK (= curcuma's
def2-universal-JKFIT), '%mp2 RI true' (RI-MP2 for the PT2 part, exact SCF), ORCA's default frozen core, '%method IntAcc 8
AngularGrid 7 GridPruning 0' (the "fine" grid of WP5-WP11) and, where ORCA's own definition differs from curcuma's,
explicit overrides:
  - XMuePBE 0.2195149727645171 for every PBE-exchange functional (ORCA's default mu is 0.21952, libxc's/curcuma's
    0.2195149727645171, WP6);
  - P86: ORCA has no PZ81 local correlation and uses VWN5 (printed LDAOpt VWN-5); curcuma runs with -qm.p86_lda vwn5;
  - the -d3bj rows and dsd-pbeb95-d3bj: ScalHFX/ScalDFX/ScalDFC/ScalLDAC, %mp2 PS (opposite spin) / PT (same spin) and
    D3S6/D3A1/D3S8/D3A2 set to curcuma's table values (ORCA's DSD keywords carry the no-dispersion parameter sets).
curcuma: -qm.grid 250x1454, -qm.frozen_core true, -qm.scf_threshold 1e-10, RI-PT2 with def2-universal-JKFIT.
Compared: E_KS (ORCA "Total Energy" of the SCF), E_PT2 (ORCA "RI-MP2 CORRELATION ENERGY", scaled), E_disp (ORCA
"Dispersion correction", two-body D3BJ) and the total.

Run:  python3 scripts/qm_dh_orca_compare.py [--orca DIR] [--curcuma release/curcuma] [--only H2O_b2plyp,...]
Claude Generated (Oct 2026, WP12 stage 1).
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "test_cases", "qm_dh", "reference_orca.json")
XYZ = {"H2O": "qm_1e/H2O.xyz", "NH3": "qm_1e/NH3.xyz", "CH4": "qm_1e/CH4.xyz", "HF": "qm_1e/HF.xyz",
       "Ne": "qm_1e/Ne.xyz", "benzene": "qm_hf3c/benzene.xyz"}
MU = "XMuePBE 0.2195149727645171"
# curcuma method -> (ORCA keyword, %method extra, %mp2 extra, curcuma extra flags)
DEF = {
    "b2plyp": ("B2PLYP", "", "", []),
    "b2gp-plyp": ("B2GP-PLYP", "", "", []),
    "dsd-blyp": ("DSD-BLYP", "", "", []),
    "dsd-pbep86": ("DSD-PBEP86", MU, "", ["-qm.p86_lda", "vwn5"]),
    "revdsd-pbep86": ("revDSD-PBEP86/2021", MU, "", ["-qm.p86_lda", "vwn5"]),
    "b2plyp-d3bj": ("B2PLYP D3BJ", "D3S6 0.64 D3A1 0.3065 D3S8 0.9147 D3A2 5.0570", "", []),
    "b2gp-plyp-d3bj": ("B2GP-PLYP D3BJ", "D3S6 0.56 D3A1 0.0 D3S8 0.2597 D3A2 6.3332", "", []),
    "dsd-blyp-d3bj": ("DSD-BLYP D3BJ", "ScalHFX 0.71 ScalDFX 0.29 ScalDFC 0.54 ScalLDAC 0.54 D3S6 0.57 D3A1 0.0 D3S8 0.0 "
                      "D3A2 5.4", "PS 0.47 PT 0.40", []),
    "dsd-pbep86-d3bj": ("DSD-PBEP86 D3BJ", MU + " ScalHFX 0.69 ScalDFX 0.31 ScalDFC 0.44 ScalLDAC 0.44 D3S6 0.48 D3A1 0.0 "
                        "D3S8 0.0 D3A2 5.6", "PS 0.52 PT 0.22", ["-qm.p86_lda", "vwn5"]),
    "dsd-pbeb95-d3bj": ("DSD-PBEB95 D3BJ", MU + " ScalHFX 0.66 ScalDFX 0.34 ScalDFC 0.55 ScalLDAC 0.55 D3S6 0.61 D3A1 0.0 "
                        "D3S8 0.0 D3A2 6.2", "PS 0.46 PT 0.09", []),
}
CASES = [(m, f) for f in ("b2plyp", "b2gp-plyp", "dsd-blyp", "dsd-pbep86", "revdsd-pbep86", "dsd-pbeb95-d3bj")
         for m in XYZ] + [(m, f) for f in ("b2plyp-d3bj", "b2gp-plyp-d3bj", "dsd-blyp-d3bj", "dsd-pbep86-d3bj")
                          for m in ("H2O", "NH3", "benzene")]


def orca_run(orca_dir, mol, fn):
    kw, meth, mp2, _ = DEF[fn]
    lines = open(os.path.join(ROOT, "test_cases", XYZ[mol])).read().splitlines()
    n = int(lines[0].split()[0])
    geom = "\n".join(lines[2:2 + n])
    inp = (f"! {kw} def2-SVP NoRI VeryTightSCF DefGrid3\n%basis AuxC \"def2/JK\" end\n"
           f"%method IntAcc 8 AngularGrid 7 GridPruning 0 {meth} end\n%mp2 RI true {mp2} end\n"
           f"* xyz 0 1\n{geom}\n*\n")
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "dh.inp")
        open(p, "w").write(inp)
        env = dict(os.environ, LD_LIBRARY_PATH=orca_dir, OMP_NUM_THREADS="1")
        out = subprocess.run([os.path.join(orca_dir, "orca"), p], capture_output=True, text=True, cwd=tmp, env=env).stdout
    if "ORCA TERMINATED NORMALLY" not in out:
        raise RuntimeError(f"ORCA failed for {mol} {fn}")
    get = lambda rx: float(re.findall(rx, out)[-1])
    r = {"e_ks": get(r"Total Energy\s+:\s+(-?\d+\.\d+) Eh"),
         "e_pt2": get(r"RI-MP2 CORRELATION ENERGY:\s+(-?\d+\.\d+)"),
         "total": get(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")}
    r["e_disp"] = get(r"Dispersion correction\s+(-?\d+\.\d+)") if "D3BJ" in kw else 0.0
    printed = {k: re.findall(rx, out)[-1] for k, rx in (("ScalHFX", r"ScalHFX\s+\.+\s+(\S+)"),
                                                          ("ScalDFC", r"ScalDFC\s+\.+\s+(\S+)"),
                                                          ("aabb", r"Scaling for aa/bb pairs\s+\.+\s+(\S+)"),
                                                          ("ab", r"Scaling for ab pairs\s+\.+\s+(\S+)"))
               if re.findall(rx, out)}
    r["printed"] = printed
    return r


def curcuma_run(exe, mol, fn):
    xyz = os.path.join(ROOT, "test_cases", XYZ[mol])
    args = [exe, "-sp", xyz, "-method", fn, "-no_bmt", "-qm.grid", "250x1454", "-qm.frozen_core", "true",
            "-qm.scf_threshold", "1e-10", "-qm.threads", "4"] + DEF[fn][3]
    with tempfile.TemporaryDirectory() as tmp:
        out = subprocess.run(args, capture_output=True, text=True, cwd=tmp).stdout
    get = lambda rx: float(re.findall(rx, out)[-1])
    r = {"e_ks": get(r"E\(KS hybrid\)\s+=\s+(-?\d+\.\d+)"), "e_pt2": get(r"E\(PT2\)\s+=\s+(-?\d+\.\d+)"),
         "total": get(r"E\(" + re.escape(fn) + r"\)\s+=\s+(-?\d+\.\d+)")}
    m = re.findall(r"E\(disp D3\(BJ\)\)\s+=\s+(-?\d+\.\d+)", out)
    r["e_disp"] = float(m[-1]) if m else 0.0
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default=os.environ.get("ORCA_DIR",
                    "/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2"))
    ap.add_argument("--curcuma", default=os.path.join(ROOT, "release", "curcuma"))
    ap.add_argument("--only", default="")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    cases = CASES
    if a.only:
        sel = set(a.only.split(","))
        cases = [c for c in CASES if f"{c[0]}_{c[1]}" in sel]
    res = {}
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        fo = {c: ex.submit(orca_run, a.orca, *c) for c in cases}
        fc = {c: ex.submit(curcuma_run, a.curcuma, *c) for c in cases}
        for c in cases:
            res[c] = (fo[c].result(), fc[c].result())
    print(f"{'case':28s} {'E_KS ORCA':>17s} {'dE_KS':>9s} {'dE_PT2':>9s} {'dE_disp':>9s} {'dE_total':>9s}  ORCA printed")
    out = {"_about": "ORCA 6.1.0 vs curcuma, scripts/qm_dh_orca_compare.py (settings in its docstring)", "d3": {},
           "cases": {}}
    for (mol, fn), (o, c) in res.items():
        key = f"{mol}_{fn}"
        print(f"{key:28s} {o['e_ks']:17.10f} {c['e_ks'] - o['e_ks']:+9.1e} {c['e_pt2'] - o['e_pt2']:+9.1e} "
              f"{c['e_disp'] - o['e_disp']:+9.1e} {c['total'] - o['total']:+9.1e}  {o['printed']}")
        out["cases"][key] = {"orca": o, "curcuma": c}
        if o["e_disp"] != 0.0:
            out["d3"][key] = o["e_disp"]
    if not a.only:
        json.dump(out, open(OUT, "w"), indent=1)
        print("wrote", OUT)


if __name__ == "__main__":
    main()
