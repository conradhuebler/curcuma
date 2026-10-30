#!/usr/bin/env python3
"""ctest checker for the native LDA (WP5). Standard library only.

Checks (--check):
  kernel      qmxc point kernels vs libxc (reference.json "kernel", PySCF 2.14 / libxc 7.0.0) on rho = 1e-12..1e2:
              eps and v of Slater, VWN5, VWN-RPA to KERNEL_TOL relative; v vs finite difference of rho eps to 1e-8
  energy      LDA/def2-SVP SCF energy and E_xc vs PySCF RKS on the identical grid (--case), ENERGY_TOL Eh;
              electrons on the grid vs PySCF to 1e-8; RI cases against PySCF density_fit with the same aux file
  gradient    a gradient request for lda succeeds since WP8 DFT (Oct 2026): gradient_failed false, hasGradient() true;
              the values are checked in test_cases/qm_dft_grad
  invariance  atom order reversed and translation: |dE| <= 1e-9 Eh; rotation (the Lebedev grid stays fixed in space):
              |dE| reported, limit ROT_TOL (sg1)
  update      UpdateMolecule A -> B equals a fresh engine on B (1e-9 Eh): grid and V_xc cache follow the geometry
  threads     1 and 4 threads give the same energy (1e-10 Eh)
  direct      integral-direct J (scf_direct on) equals the stored-tensor J path (1e-8 Eh)
  rij         -qm.ri j runs and converges; |E(RI-J) - E(exact)| < 1e-3 Eh (a method variant, no reference value)

Claude Generated (Oct 2026, WP5).
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

KERNEL_TOL = 1e-12
ENERGY_TOL = 1e-8
ROT_TOL = 1e-4


def run(exe, args):
    p = subprocess.run([exe] + args, capture_output=True, text=True)
    if p.returncode not in (0, 2):
        print(p.stdout, p.stderr)
        raise SystemExit(f"test_qm_lda failed ({p.returncode}): {' '.join(args)}")
    return json.loads(p.stdout)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, float(x), float(y), float(z)])
    return atoms


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for s, x, y, z in atoms:
            f.write(f"{s} {x:.12f} {y:.12f} {z:.12f}\n")


def rotate(atoms, axis, angle):
    ax = [a / math.sqrt(sum(c * c for c in axis)) for a in axis]
    c, s = math.cos(angle), math.sin(angle)
    x, y, z = ax
    R = [[c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
         [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
         [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]]
    return [[a[0]] + [sum(R[i][k] * a[1 + k] for k in range(3)) for i in range(3)] for a in atoms]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--case", default="H2O_sg1")
    a = ap.parse_args()
    ref = json.load(open(a.ref))
    tc = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..")
    ok = True

    if a.check == "kernel":
        out = run(a.exe, ["kernel", a.ref])
        rho = ref["kernel"]["rho"]
        for key in ("slater", "vwn5", "vwn_rpa"):
            r = ref["kernel"][key]
            de = max(rel(e, er) for e, er in zip(out[key]["eps"], r["eps"]))
            dv = max(rel(v, vr) for v, vr in zip(out[key]["v"], r["v"]))
            dfd = max(rel(v, vf) for v, vf in zip(out[key]["v"], out[key]["v_fd"]))
            good = de <= KERNEL_TOL and dv <= KERNEL_TOL and dfd <= 1e-8
            ok &= good
            print(f"{key:8s} vs {r['libxc']:14s} rho {rho[0]:.0e}..{rho[-1]:.0e}: max rel eps {de:.2e}, v {dv:.2e} "
                  f"(tol {KERNEL_TOL:.0e}); v vs FD of rho eps {dfd:.2e} (tol 1e-8) {'OK' if good else 'FAIL'}")
    elif a.check == "energy":
        r = ref["energies"][a.case]
        out = run(a.exe, ["energy", os.path.join(tc, r["xyz"]), "--grid", r["grid"], "--vwn", r["vwn_variant"],
                          "--ri", r.get("ri", "off")])
        dE, dX, dN = out["energy"] - r["energy"], out["exc"] - r["exc"], out["nelec_grid"] - r["nelec_grid"]
        good = out["converged"] and abs(dE) <= ENERGY_TOL and abs(dX) <= ENERGY_TOL and abs(dN) <= 1e-8
        ok &= good
        print(f"{a.case}: E {out['energy']:.10f} (PySCF {r['energy']:.10f}, dE {dE:.2e}); E_xc dE {dX:.2e}; "
              f"N_grid {out['nelec_grid']:.8f} (dN {dN:.1e}); {out['iterations']} it; {out['label']} "
              f"{'OK' if good else 'FAIL'} (tol {ENERGY_TOL:.0e})")
    elif a.check == "gradient":
        out = run(a.exe, ["energy", os.path.join(tc, "qm_1e/H2O.xyz"), "--gradient"])
        good = out.get("gradient_failed") is False and out["has_gradient"] is True  # WP8 DFT (Oct 2026)
        ok &= good
        print(f"gradient request: failed={out.get('gradient_failed')} has_gradient={out['has_gradient']} "
              f"'{out.get('gradient_error')}' {'OK' if good else 'FAIL'}")
    elif a.check in ("invariance", "update", "threads", "direct", "rij"):
        mol = a.case.split("_")[0]
        xyz = os.path.join(tc, {"H2O": "qm_1e/H2O.xyz", "NH3": "qm_1e/NH3.xyz", "CH4": "qm_1e/CH4.xyz",
                                "benzene": "qm_hf3c/benzene.xyz"}[mol])
        base = run(a.exe, ["energy", xyz])["energy"]
        with tempfile.TemporaryDirectory() as td:
            atoms = read_xyz(xyz)
            if a.check == "invariance":
                tests = {"reversed": (list(reversed(atoms)), 1e-9),
                         "translated": ([[s, x + 1.3, y - 0.7, z + 2.1] for s, x, y, z in atoms], 1e-9),
                         "rotated": (rotate(atoms, (1.0, 2.0, 3.0), 0.7), ROT_TOL)}
                for name, (at, tol) in tests.items():
                    p = os.path.join(td, name + ".xyz")
                    write_xyz(p, at)
                    e = run(a.exe, ["energy", p])["energy"]
                    good = abs(e - base) <= tol
                    ok &= good
                    print(f"{mol} {name}: E {e:.10f}, dE {e - base:.2e} (tol {tol:.0e}) {'OK' if good else 'FAIL'}")
            elif a.check == "update":
                p = os.path.join(td, "moved.xyz")
                moved = [[s, x * 1.03, y * 0.98, z + 0.02 * i] for i, (s, x, y, z) in enumerate(atoms)]
                write_xyz(p, moved)
                out = run(a.exe, ["update", xyz, p])
                d = out["energy_B_updated"] - out["energy_B_fresh"]
                good = out["converged"] and abs(d) <= 1e-9 and abs(out["energy_A"] - out["energy_B_fresh"]) > 1e-6
                ok &= good
                print(f"{mol} update: E_A {out['energy_A']:.10f}, E_B updated {out['energy_B_updated']:.10f}, fresh "
                      f"{out['energy_B_fresh']:.10f}, diff {d:.2e} (warm start {out['warm_start']}) "
                      f"{'OK' if good else 'FAIL'}")
            elif a.check == "threads":
                e1 = run(a.exe, ["energy", xyz, "--threads", "1"])["energy"]
                good = abs(e1 - base) <= 1e-10
                ok &= good
                print(f"{mol} threads 1 vs 4: {e1:.12f} vs {base:.12f}, diff {e1 - base:.2e} {'OK' if good else 'FAIL'}")
            elif a.check == "direct":
                e = run(a.exe, ["energy", xyz, "--direct", "on"])["energy"]
                good = abs(e - base) <= 1e-8
                ok &= good
                print(f"{mol} integral-direct vs stored: {e:.10f} vs {base:.10f}, diff {e - base:.2e} "
                      f"{'OK' if good else 'FAIL'}")
            else:
                out = run(a.exe, ["energy", xyz, "--ri", "j"])
                d = out["energy"] - base
                good = out["converged"] and abs(d) < 1e-3
                ok &= good
                print(f"{mol} RI-J: {out['energy']:.10f}, exact {base:.10f}, fit error {d:.2e} {'OK' if good else 'FAIL'}")
    else:
        raise SystemExit("unknown check " + a.check)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
