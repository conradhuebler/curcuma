#!/usr/bin/env python3
"""Double-hybrid analytic gradient (WP12 stage 2) against ORCA 6.1 on the settings of WP12 stage 1.

ORCA: '! <keyword> def2-SVP NoRI VeryTightSCF DefGrid3 EnGrad', AuxC def2/JK, '%mp2 RI true', ORCA's default frozen
core, '%method IntAcc 8 AngularGrid 7 GridPruning 0' plus the overrides of scripts/qm_dh_orca_compare.py (XMuePBE,
scaling and D3 parameters for the -d3bj rows); gradient from <base>.engrad (Eh/Bohr).
curcuma: test_qm_dh_grad grad --grid 250x1454 --fc (P86 with --p86 vwn5), grid response on.
Prints max |g_curcuma - g_ORCA| per case and writes test_cases/qm_dh_grad/reference_orca.json (energies and gradients).
Run:  python3 scripts/qm_dh_grad_orca_compare.py --exe release/test_cases/qm_dh_grad/test_qm_dh_grad [--only H2O_b2plyp]
Claude Generated (Oct 2026, WP12 stage 2).
"""
import argparse
import concurrent.futures as cf
import importlib.util
import json
import os
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "test_cases", "qm_dh_grad", "reference_orca.json")
spec = importlib.util.spec_from_file_location("qm_dh_orca_compare", os.path.join(HERE, "qm_dh_orca_compare.py"))
OC = importlib.util.module_from_spec(spec)
spec.loader.exec_module(OC)
CASES = [(m, f) for f in ("b2plyp", "dsd-blyp", "dsd-pbep86") for m in ("H2O", "NH3")] + \
        [("H2O", "b2plyp-d3bj"), ("NH3", "revdsd-pbep86")]
# Not compared: dsd-pbeb95-d3bj -- ORCA 6.1 stops with "Gradient is not yet implemented for the PWPB95, PW6B95 and
# DSD-PBEB95 funtionals" (NumGrad or its LibXC version would be needed; not run).


def orca_grad(orca_dir, mol, fn):
    kw, meth, mp2, _ = OC.DEF[fn]
    lines = open(os.path.join(ROOT, "test_cases", OC.XYZ[mol])).read().splitlines()
    n = int(lines[0].split()[0])
    geom = "\n".join(lines[2:2 + n])
    inp = (f"! {kw} def2-SVP NoRI VeryTightSCF DefGrid3 EnGrad\n%basis AuxC \"def2/JK\" end\n"
           f"%method IntAcc 8 AngularGrid 7 GridPruning 0 {meth} end\n%mp2 RI true {mp2} end\n"
           f"* xyz 0 1\n{geom}\n*\n")
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "dh.inp")
        open(p, "w").write(inp)
        env = dict(os.environ, LD_LIBRARY_PATH=orca_dir, OMP_NUM_THREADS="1")
        out = subprocess.run([os.path.join(orca_dir, "orca"), p], capture_output=True, text=True, cwd=tmp, env=env).stdout
        if "ORCA TERMINATED NORMALLY" not in out:
            raise RuntimeError(f"ORCA failed for {mol} {fn}: {out[-1500:]}")
        L = [l.strip() for l in open(os.path.join(tmp, "dh.engrad")) if l.strip() and not l.startswith("#")]
    nat = int(L[0])
    energy = float(L[1])
    grad = [float(x) for x in L[2:2 + 3 * nat]]
    return {"energy": energy, "gradient": grad, "input": inp}


def curcuma_grad(exe, mol, fn):
    extra = OC.DEF[fn][3]
    args = [exe, "grad", os.path.join(ROOT, "test_cases", OC.XYZ[mol]), "--method", fn, "--grid", "250x1454", "--fc",
            "--threads", "4"]
    if "-qm.p86_lda" in extra:
        args += ["--p86", "vwn5"]
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    out = json.loads(p.stdout[p.stdout.find("{\n"):])
    out["args"] = args[2:]
    out["args"][0] = OC.XYZ[mol]  # path relative to test_cases
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default="/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    cases = [c for c in CASES if not a.only or f"{c[0]}_{c[1]}" in a.only.split(",")]
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(orca_grad, a.orca, m, f): (m, f) for m, f in cases}
        res = {}
        for fu in cf.as_completed(futs):
            res[futs[fu]] = fu.result()
    out = {"generator": "scripts/qm_dh_grad_orca_compare.py", "orca": "6.1.0", "cases": {}}
    for m, f in cases:
        c = curcuma_grad(os.path.abspath(a.exe), m, f)
        o = res[(m, f)]
        d = max(abs(x - y) for x, y in zip(c["gradient"], o["gradient"]))
        out["cases"][f"{m}_{f}"] = {"orca_energy": o["energy"], "orca_gradient": o["gradient"],
                                    "curcuma_energy": c["energy"], "curcuma_gradient": c["gradient"], "max_dg": d,
                                    "curcuma_args": c["args"]}
        print(f"{m:4s} {f:16s} dE {c['energy'] - o['energy']:+.2e} Eh  max|g - g_ORCA| {d:.2e} Eh/Bohr  "
              f"max|g| {max(map(abs, o['gradient'])):.3e}", flush=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))


if __name__ == "__main__":
    main()
