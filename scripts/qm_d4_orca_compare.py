#!/usr/bin/env python3
"""Native DFT-D4 (WP12 stage 4) against ORCA 6.1; writes test_cases/qm_dh/reference_orca_d4.json.

Part 1, D4 term alone ("disp"): ORCA prints the D4 energy ("Dispersion correction", 9 decimals) and, for a KS
functional with EnGrad, the "DISPERSION GRADIENT" block (9 decimals, Eh/Bohr). The D4 term does not depend on the
orbitals, so every case runs '! PBE D4 STO-3G NoRI EnGrad' with the parameter set given explicitly
(%method D4S6 D4S8 D4A1 D4A2 D4S9 end; ORCA's dftd4 is v3.4.0 per its output) -- plus the ORCA keyword itself
('! B2PLYP D4' etc., "keyword" entries) to show that ORCA's built-in set is the same. Exception: ORCA's built-in
revDSD-PBEP86-D4/2021 set has a2 = 4.2014 (fitted to its printed energies: 5 molecules within 5e-10 Eh), not
4.201 of Mehta et al. 2023 Table 1, so its keyword entry differs from curcuma by up to 2.4e-6 Eh (benzene).

Part 2, whole double hybrid ("dh"): revDSD-PBEP86-D4 total energy and gradient on the settings of
scripts/qm_dh_orca_compare.py ('! revDSD-PBEP86-D4/2021 def2-SVP NoRI VeryTightSCF DefGrid3 EnGrad', AuxC def2/JK,
%mp2 RI true, ORCA's frozen core, IntAcc 8 AngularGrid 7 GridPruning 0, XMuePBE = libxc mu, D4A2 4.201 so that the
D4 parameters equal curcuma's; P86 on VWN5 = ORCA and curcuma default). curcuma: -sp -qm.grid 250x1454
-qm.frozen_core true -qm.scf_threshold 1e-10 -dump_gradient.

Run:  python3 scripts/qm_d4_orca_compare.py [--orca DIR] [--curcuma release/curcuma] [--part disp|dh|all]
Claude Generated (Oct 2026, WP12 stage 4).
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
OUT = os.path.join(ROOT, "test_cases", "qm_dh", "reference_orca_d4.json")
BOHR = 0.529177210903
XYZ = {"H2O": "qm_1e/H2O.xyz", "NH3": "qm_1e/NH3.xyz", "CH4": "qm_1e/CH4.xyz",
       "H2CO": "qm_grad/formaldehyde_hf3c_min_reference.xyz", "benzene": "qm_hf3c/benzene.xyz"}
# curcuma functional key -> (ORCA keyword, s6, s8, a1, a2, s9) = curcuma's d4DftParameters()
PARAMS = {
    "pbe": ("PBE D4", 1.0, 0.95948085, 0.38574991, 4.80688534, 1.0),
    "b3lyp": ("B3LYP D4", 1.0, 2.02929367, 0.40868035, 4.53807137, 1.0),
    "tpss": ("TPSS D4", 1.0, 1.76596355, 0.42822303, 4.54257102, 1.0),
    "scan": ("SCANfunc D4", 1.0, 1.46126056, 0.62930855, 6.31284039, 1.0),
    "b2plyp": ("B2PLYP D4", 0.64, 1.16888646, 0.44154604, 4.73114642, 1.0),
    "b2gp-plyp": ("B2GP-PLYP D4", 0.56, 0.94633372, 0.42907301, 5.18802602, 1.0),
    "revdsd-pbep86": ("revDSD-PBEP86-D4/2021", 0.5917, 0.0, 0.371, 4.201, 1.0),
}
MU = "XMuePBE 0.2195149727645171"


def geom(mol):
    lines = open(os.path.join(ROOT, "test_cases", XYZ[mol])).read().splitlines()
    n = int(lines[0].split()[0])
    return "\n".join(lines[2:2 + n]), n


def run_orca(orca_dir, inp):
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "job.inp")
        open(p, "w").write(inp)
        env = dict(os.environ, LD_LIBRARY_PATH=orca_dir, OMP_NUM_THREADS="1")
        out = subprocess.run([os.path.join(orca_dir, "orca"), p], capture_output=True, text=True, cwd=tmp,
                             env=env).stdout
        if "ORCA TERMINATED NORMALLY" not in out:
            raise RuntimeError(out[-2000:])
        eg = os.path.join(tmp, "job.engrad")
        grad = None
        if os.path.exists(eg):
            L = [l.strip() for l in open(eg) if l.strip() and not l.startswith("#")]
            nat = int(L[0])
            grad = [float(x) for x in L[2:2 + 3 * nat]]
    return out, grad


def disp_case(orca_dir, mol, fn, keyword):
    kw, s6, s8, a1, a2, s9 = PARAMS[fn]
    g, n = geom(mol)
    if keyword:
        if fn in ("b2plyp", "b2gp-plyp", "revdsd-pbep86"):
            inp = f"! {kw} STO-3G NoRI\n%mp2 RI false end\n* xyz 0 1\n{g}\n*\n"
        else:
            inp = f"! {kw} STO-3G NoRI\n* xyz 0 1\n{g}\n*\n"
    else:
        inp = (f"! PBE D4 STO-3G NoRI EnGrad\n%method D4S6 {s6} D4S8 {s8} D4A1 {a1} D4A2 {a2} D4S9 {s9} end\n"
               f"* xyz 0 1\n{g}\n*\n")
    out, _ = run_orca(orca_dir, inp)
    r = {"energy": float(re.findall(r"Dispersion correction\s+(-?\d+\.\d+)", out)[-1]), "input": inp}
    if not keyword:
        blk = out.split("DISPERSION GRADIENT")[-1]
        rows = re.findall(r"^\s*\d+\s+\w+\s+:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", blk, re.M)[:n]
        r["gradient"] = [float(v) for row in rows for v in row]
    return r


def dh_case(orca_dir, exe, mol):
    g, n = geom(mol)
    inp = (f"! revDSD-PBEP86-D4/2021 def2-SVP NoRI VeryTightSCF DefGrid3 EnGrad\n%basis AuxC \"def2/JK\" end\n"
           f"%method IntAcc 8 AngularGrid 7 GridPruning 0 {MU} D4A2 4.201 end\n%mp2 RI true end\n"
           f"* xyz 0 1\n{g}\n*\n")
    out, grad = run_orca(orca_dir, inp)
    o = {"total": float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)[-1]),
         "e_disp": float(re.findall(r"Dispersion correction\s+(-?\d+\.\d+)", out)[-1]), "gradient": grad,
         "input": inp}
    xyz = os.path.join(ROOT, "test_cases", XYZ[mol])
    with tempfile.TemporaryDirectory() as tmp:
        gp = os.path.join(tmp, "g.txt")
        args = [exe, "-sp", xyz, "-method", "revdsd-pbep86-d4", "-no_bmt", "-qm.grid", "250x1454",
                "-qm.frozen_core", "true", "-qm.scf_threshold", "1e-10", "-qm.threads", "4", "-dump_gradient", gp]
        txt = subprocess.run(args, capture_output=True, text=True, cwd=tmp).stdout
        rows = [l.split() for l in open(gp) if l.strip() and not l.startswith("#")]
    c = {"total": float(re.findall(r"E\(revdsd-pbep86-d4\)\s+=\s+(-?\d+\.\d+)", txt)[-1]),
         "e_disp": float(re.findall(r"E\(disp D4\)\s+=\s+(-?\d+\.\d+)", txt)[-1]),
         "gradient": [float(v) * BOHR for r in rows for v in r[:3]], "args": args[3:]}
    c["args"][0] = XYZ[mol]
    return o, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default="/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2")
    ap.add_argument("--curcuma", default=os.path.join(ROOT, "release", "curcuma"))
    ap.add_argument("--part", default="all")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    out["_about"] = "ORCA 6.1.0 (dftd4 3.4.0 inside), scripts/qm_d4_orca_compare.py (settings in its docstring)"
    if a.part in ("disp", "all"):
        jobs = [(m, f, k) for m in XYZ for f in PARAMS for k in (False, True)]
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            fut = {j: ex.submit(disp_case, a.orca, *j) for j in jobs}
            out["disp"], out["keyword"] = {}, {}
            for (m, f, k), fu in fut.items():
                (out["keyword"] if k else out["disp"])[f"{m}_{f}"] = fu.result()
        for key, r in out["disp"].items():
            print(f"{key:24s} E_D4 explicit {r['energy']:.9f}  keyword {out['keyword'][key]['energy']:.9f}")
    if a.part in ("dh", "all"):
        out["dh"] = {}
        for m in ("H2O", "NH3", "H2CO"):
            o, c = dh_case(a.orca, os.path.abspath(a.curcuma), m)
            d = max(abs(x - y) for x, y in zip(o["gradient"], c["gradient"]))
            out["dh"][f"{m}_revdsd-pbep86-d4"] = {"orca": o, "curcuma": c, "max_dg": d}
            print(f"{m:8s} revdsd-pbep86-d4 dE {c['total'] - o['total']:+.2e}  dE_D4 {c['e_disp'] - o['e_disp']:+.2e}  "
                  f"max|dg| {d:.2e} Eh/Bohr", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
