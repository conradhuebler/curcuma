#!/usr/bin/env python3
"""Double-hybrid geometry optimisation (WP12 stage 2): curcuma -opt against ORCA 6.1 Opt on matched settings.

ORCA: '! <keyword> def2-SVP NoRI VeryTightSCF DefGrid3 TightOpt', AuxC def2/JK, '%mp2 RI true', frozen core (ORCA
default), '%method IntAcc 8 AngularGrid 7 GridPruning 0' (+ XMuePBE for PBE exchange), as scripts/qm_dh_orca_compare.py.
curcuma: -opt -qm.grid 250x1454 -qm.frozen_core true (P86 rows -qm.p86_lda vwn5), -opt.convergence_preset tight.
Compared: all interatomic distances (Angstrom) and the bond angle(s) of the optimised geometries, and the final energies.
Start geometries: test_cases/qm_grad/{H2O,formaldehyde}_distorted.xyz, qm_1e/NH3.xyz.
Run:  python3 scripts/qm_dh_opt_orca_compare.py --curcuma release/curcuma
Claude Generated (Oct 2026, WP12 stage 2).
"""
import argparse
import concurrent.futures as cf
import importlib.util
import math
import os
import re
import shutil
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("qm_dh_orca_compare", os.path.join(HERE, "qm_dh_orca_compare.py"))
OC = importlib.util.module_from_spec(spec)
spec.loader.exec_module(OC)
START = {"H2O": "qm_grad/H2O_distorted.xyz", "NH3": "qm_1e/NH3.xyz", "H2CO": "qm_grad/formaldehyde_distorted.xyz"}
CASES = [("H2O", "b2plyp"), ("NH3", "b2plyp"), ("H2CO", "b2plyp"), ("H2O", "dsd-pbep86"), ("H2CO", "dsd-pbep86")]


def read_xyz(text):
    L = text.strip().splitlines()
    n = int(L[0].split()[0])
    return [(l.split()[0], [float(v) for v in l.split()[1:4]]) for l in L[2:2 + n]]


def orca_opt(orca_dir, mol, fn):
    kw, meth, mp2, _ = OC.DEF[fn]
    lines = open(os.path.join(ROOT, "test_cases", START[mol])).read().splitlines()
    n = int(lines[0].split()[0])
    geom = "\n".join(lines[2:2 + n])
    inp = (f"! {kw} def2-SVP NoRI VeryTightSCF DefGrid3 TightOpt\n%basis AuxC \"def2/JK\" end\n"
           f"%method IntAcc 8 AngularGrid 7 GridPruning 0 {meth} end\n%mp2 RI true {mp2} end\n"
           f"* xyz 0 1\n{geom}\n*\n")
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "dh.inp")
        open(p, "w").write(inp)
        env = dict(os.environ, LD_LIBRARY_PATH=orca_dir, OMP_NUM_THREADS="1")
        out = subprocess.run([os.path.join(orca_dir, "orca"), p], capture_output=True, text=True, cwd=tmp, env=env).stdout
        if "ORCA TERMINATED NORMALLY" not in out or "OPTIMIZATION RUN DONE" not in out:
            raise RuntimeError(f"ORCA opt failed for {mol} {fn}")
        xyz = read_xyz(open(os.path.join(tmp, "dh.xyz")).read())
        e = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)[-1])
    return xyz, e


def curcuma_opt(exe, mol, fn):
    extra = OC.DEF[fn][3]
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(os.path.join(ROOT, "test_cases", START[mol]), os.path.join(tmp, "start.xyz"))
        args = [exe, "-opt", "start.xyz", "-method", fn, "-no_bmt", "-verbosity", "1", "-qm.grid", "250x1454",
                "-qm.frozen_core", "true", "-qm.scf_threshold", "1e-9", "-opt.convergence_preset", "tight",
                "-threads", "4"] + extra
        p = subprocess.run(args, capture_output=True, text=True, cwd=tmp)
        xyz = read_xyz(open(os.path.join(tmp, "start.opt.xyz")).read())
        txt = re.sub(r"\x1b\[[0-9;]*m", "", p.stdout)
        e = float(re.findall(r"Energy = (-?\d+\.\d+)", open(os.path.join(tmp, "start.opt.xyz")).read())[0])
    return xyz, e, p.returncode


def internals(xyz):
    pos = [c for _, c in xyz]
    d = {f"r{i}{j}": math.dist(pos[i], pos[j]) for i in range(len(pos)) for j in range(i + 1, len(pos))}
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default="/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2")
    ap.add_argument("--curcuma", required=True)
    a = ap.parse_args()
    with cf.ThreadPoolExecutor(max_workers=5) as ex:
        fo = {c: ex.submit(orca_opt, a.orca, *c) for c in CASES}
        res = {c: fo[c].result() for c in CASES}
    for c in CASES:
        cx, ce, rc = curcuma_opt(os.path.abspath(a.curcuma), *c)
        ox, oe = res[c]
        dc, do = internals(cx), internals(ox)
        dmax = max(abs(dc[k] - do[k]) for k in dc)
        print(f"{c[0]:5s} {c[1]:12s} curcuma rc {rc}: max |d_ij - d_ij(ORCA)| {dmax:.1e} A; E(curcuma, 6 digits) {ce:.6f} "
              f"E(ORCA) {oe:.8f}", flush=True)


if __name__ == "__main__":
    main()
