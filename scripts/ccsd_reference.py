#!/usr/bin/env python3
"""Reference CCSD energies for the native canonical RCCSD (test_cases/qm_ccsd/reference.json).

PySCF RHF (spherical def2-SVP from PySCF's own library, = curcuma's def2-SVP.dat; conv_tol 1e-12,
conv_tol_grad 1e-8) followed by cc.CCSD(mf) (canonical, exact integrals; conv_tol 1e-12,
conv_tol_normt 1e-10) with and without the chemical frozen core (frozen = number of 1s orbitals of
Li-Ne atoms, as qmint::frozenCoreOrbitals). Stored per case: E_HF, E_corr(CCSD), E_corr(MP2), the
T1 diagnostic and PySCF's iteration count, for both settings (Hartree).

With --orca the same cases are run with ORCA 6.1 (`! RHF CCSD def2-SVP ExtremeSCF [No]FrozenCore`,
%mdci STol 1e-10) and E_HF / E_corr (from job.property.txt, full precision) / T1 diagnostic are stored
under "orca" (second reference; the ctest compares with PySCF). ORCA scratch lives in
~/.cache/curcuma_ccsd_orca (never /tmp: RAM tmpfs); a watchdog kills ORCA when the disk there has less
than 20 GB free.

Run with the PySCF venv:
    ~/.venvs/qmref/bin/python scripts/ccsd_reference.py [--orca] [names...]
    env ORCA=<orca binary> (default the local 6.1 install)
Claude Generated (Oct 2026, WP13).
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time

import pyscf
from pyscf import cc, gto, mp, scf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "test_cases/qm_ccsd/reference.json")
ORCA = os.environ.get("ORCA", "/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca")
ORCA_WORK = os.path.expanduser("~/.cache/curcuma_ccsd_orca")

CASES = [  # (name, xyz relative to test_cases/)
    ("H2O", "qm_1e/H2O.xyz"),
    ("NH3", "qm_1e/NH3.xyz"),
    ("CH4", "qm_1e/CH4.xyz"),
    ("formaldehyde", "qm_hf3c/formaldehyde.xyz"),
    ("ethane", "structures/small/ethane.xyz"),
    ("water_dimer", "qm_hf3c/water_dimer.xyz"),
    ("benzene", "qm_hf3c/benzene.xyz"),
]


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]


def frozen_core(atoms):
    """Chemical frozen core as in qmint::frozenCoreOrbitals (H-He 0, Li-Ne 1)."""
    return sum(0 if gto.charge(s) <= 2 else 1 for s, _ in atoms)


def pyscf_case(atoms):
    mol = gto.M(atom=[[s, x] for s, x in atoms], basis="def2-svp", unit="Angstrom", cart=False, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-12, 1e-8, 200
    e_hf = mf.kernel()
    if not mf.converged:
        sys.exit("PySCF RHF did not converge")
    entry = {"nbf": mol.nao_nr(), "hf": e_hf}
    nfc = frozen_core(atoms)
    for tag, frz in (("full", None), ("fc", nfc if nfc > 0 else None)):
        e_mp2 = mp.MP2(mf, frozen=frz).kernel()[0]
        c = cc.CCSD(mf, frozen=frz)
        c.conv_tol, c.conv_tol_normt, c.max_cycle = 1e-12, 1e-10, 200
        niter = [0]
        c.callback = lambda env: niter.__setitem__(0, niter[0] + 1)
        t0 = time.time()
        c.kernel()
        if not c.converged:
            sys.exit("PySCF CCSD did not converge")
        entry[tag] = {"frozen": frz or 0, "e_ccsd": c.e_corr, "e_mp2": e_mp2,
                      "t1_diagnostic": float(c.get_t1_diagnostic()), "pyscf_iterations": niter[0],
                      "pyscf_seconds": round(time.time() - t0, 2)}
    return entry


def run_orca(name, atoms, fc):
    d = os.path.join(ORCA_WORK, f"{name}_{'fc' if fc else 'full'}")
    os.makedirs(d, exist_ok=True)
    outp = os.path.join(d, "job.out")
    if not (os.path.exists(outp) and "ORCA TERMINATED NORMALLY" in open(outp).read()):
        with open(os.path.join(d, "job.inp"), "w") as f:
            f.write(f"! RHF CCSD def2-SVP ExtremeSCF {'FrozenCore' if fc else 'NoFrozenCore'}\n"
                    "%mdci\n  STol 1e-10\n  MaxIter 200\nend\n%maxcore 8000\n* xyz 0 1\n")
            for s, x in atoms:
                f.write(f"{s} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}\n")
            f.write("*\n")
        env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(ORCA))
        with open(outp, "w") as fo:
            p = subprocess.Popen([ORCA, "job.inp"], cwd=d, stdout=fo, stderr=subprocess.STDOUT, env=env)
            while p.poll() is None:  # disk watchdog
                if shutil.disk_usage(d).free < 20e9:
                    p.kill()
                    sys.exit(f"ORCA killed: less than 20 GB free in {d}")
                time.sleep(2)
    txt = open(outp).read()
    if "ORCA TERMINATED NORMALLY" not in txt:
        return None
    prop = open(os.path.join(d, "job.property.txt")).read()
    val = lambda key: float(re.findall(key + r"[^\n]*\n[^\n]*\n\s*\n\s*0\s+(-?\d\.\d+e[-+]\d+)", prop)[-1])
    t1 = re.findall(r"T1 diagnostic\s+\.\.\.\s+(-?\d+\.\d+)", txt)
    for f in os.listdir(d):  # keep only input, output and property file
        if f not in ("job.inp", "job.out", "job.property.txt"):
            try:
                os.remove(os.path.join(d, f))
            except OSError:
                pass
    return {"hf": val("&refEnergy"), "e_ccsd": val("&corrEnergy"), "t1_diagnostic": float(t1[-1]) if t1 else None}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    orca = "--orca" in sys.argv
    old = json.load(open(OUT)) if os.path.exists(OUT) else {"cases": {}}
    out = {"_about": "PySCF canonical RCCSD/def2-SVP references for ctest -L qm_ccsd (scripts/ccsd_reference.py)",
           "_sources": {"pyscf": pyscf.__version__, "basis": "def2-svp (PySCF library), spherical",
                        "scf": "RHF conv_tol 1e-12, conv_tol_grad 1e-8",
                        "ccsd": "cc.CCSD(mf), conv_tol 1e-12, conv_tol_normt 1e-10, frozen=None and frozen=n_core",
                        "orca": "ORCA 6.1.0, ! RHF CCSD def2-SVP ExtremeSCF [No]FrozenCore, %mdci STol 1e-10 "
                                "(second reference, not used by the ctest)"},
           "_units": "Hartree", "cases": dict(old.get("cases", {}))}
    for name, rel in CASES:
        if args and name not in args:
            continue
        atoms = read_xyz(os.path.join(ROOT, "test_cases", rel))
        entry = pyscf_case(atoms)
        entry = {"xyz": rel, **entry}
        if orca:
            entry["orca"] = {"full": run_orca(name, atoms, False), "fc": run_orca(name, atoms, True)}
        elif name in out["cases"] and "orca" in out["cases"][name]:
            entry["orca"] = out["cases"][name]["orca"]
        out["cases"][name] = entry
        msg = (f"{name:14s} nbf={entry['nbf']:4d} E_HF={entry['hf']:.10f} CCSD={entry['full']['e_ccsd']:.12f} "
               f"CCSD(fc)={entry['fc']['e_ccsd']:.12f}")
        if "orca" in entry and entry["orca"]["full"]:
            o = entry["orca"]
            msg += (f"  ORCA-PySCF: HF {o['full']['hf'] - entry['hf']:.1e} corr {o['full']['e_ccsd'] - entry['full']['e_ccsd']:.1e}"
                    f" corr(fc) {o['fc']['e_ccsd'] - entry['fc']['e_ccsd']:.1e}")
        print(msg, flush=True)
        json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
