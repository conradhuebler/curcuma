#!/usr/bin/env python3
"""Gates for the native Kohn-Sham gradient (WP8 DFT); reference.json from scripts/qm_dft_grad_pyscf_reference.py.

Checks (gradients in Eh/Bohr; standard library only, no numpy):
  pyscf --case C     curcuma gradient without grid response vs PySCF grid_response=False and with grid response vs
                     PySCF grid_response=True, identical grid (WP4); tol 1e-9 (observed <= 3.6e-11)
  fd --case M_F --fd moving|frozen
                     analytic gradient vs the central 4-point finite difference (h = 2e-3 Bohr) of curcuma's own energy:
                     moving grid (rebuilt per geometry) vs the gradient WITH grid response, frozen grid (fixed in space,
                     QMEngine::freezeXCGrid) vs the gradient WITHOUT grid response; tol 2e-9 (observed <= 3.4e-10)
  translation --case M_F
                     sum of the atomic gradients: with grid response < 1e-11 (observed 1e-14), without it the grid
                     error is reported (H2O/NH3 sg1 ~1e-5); invariance of the gradient under translation < 1e-11
  rotation --case M_F
                     gradient at a rotated geometry, rotated back, vs the original: the grid does not rotate, so this
                     is the grid error (sg1: <= 1e-4 measured); gate 3e-4 at sg1
  threads --case M_F 1 vs 4 threads bitwise (exact path)
  update             UpdateMolecule A -> B vs a fresh engine on B, max |dg| < 1e-10
  opt --case M_F     curcuma -opt -method F from qm_grad/<M>.xyz; the optimised geometry (all interatomic distances,
                     Angstrom) and the energy at it against the PySCF/scipy-BFGS minimum (same grid, grid response on);
                     the energy gate is 1e-6 because the final orientation may differ (grid not rotationally
                     invariant); in addition curcuma's energy and gradient AT the PySCF minimum (1e-9 Eh, 1e-6 Eh/Bohr)
Claude Generated (Oct 2026, WP8; --h/--tol and the meta-GGA cases of test_cases/qm_mgga_grad WP11 stage 2).
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

TC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XYZ = {"H2O": "qm_1e/H2O.xyz", "NH3": "qm_1e/NH3.xyz", "CH4": "qm_1e/CH4.xyz", "HF": "qm_1e/HF.xyz",
       "benzene": "qm_hf3c/benzene.xyz", "H2O_distorted": "qm_grad/H2O_distorted.xyz",
       "formaldehyde_distorted": "qm_grad/formaldehyde_distorted.xyz"}


def run(exe, args):
    p = subprocess.run([exe] + args, capture_output=True, text=True, timeout=3000)
    if p.returncode != 0:
        print(p.stderr)
        sys.exit(f"{exe} {' '.join(args)} failed ({p.returncode})")
    return json.loads(p.stdout)


def maxdiff(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def parse_case(case):
    """'H2O_b3lyp_vwn_rpa_rij' -> molecule, options"""
    for mol in sorted(XYZ, key=len, reverse=True):
        if case.startswith(mol + "_"):
            rest = case[len(mol) + 1:].split("_")
            break
    else:
        sys.exit(f"unknown case {case}")
    opts = ["--functional", rest[0]]
    rest = rest[1:]
    if rest[:2] == ["vwn", "rpa"]:
        opts += ["--vwn", "vwn_rpa"]
        rest = rest[2:]
    for r in rest:
        if r in ("rij", "rijk"):
            opts += ["--ri", "j" if r == "rij" else "jk"]
        elif r.startswith("sg") or "x" in r:
            opts += ["--grid", r]
    return mol, opts


def distances(xyz):
    n = len(xyz)
    return [math.dist(xyz[i], xyz[j]) for i in range(n) for j in range(i + 1, n)]


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--case", default="")
    ap.add_argument("--fd", default="moving")
    ap.add_argument("--curcuma", default="")
    ap.add_argument("--h", default="")      # FD step (Bohr), default of the driver: 2e-3 (WP11 stage 2)
    ap.add_argument("--tol", type=float, default=2e-9)  # FD gate (WP11 stage 2: meta-GGA cases)
    a = ap.parse_args()
    ref = json.load(open(a.ref))
    ok = True
    if a.check == "pyscf":
        r = ref["gradients"][a.case]
        base = ["grad", os.path.join(TC, r["xyz"]), "--functional", r["functional"], "--vwn", r["vwn"], "--grid",
                r["grid"], "--ri", r["ri"]]
        g0 = run(a.exe, base + ["--response", "0"])
        g1 = run(a.exe, base + ["--response", "1"])
        d0, d1 = maxdiff(g0["gradient"], r["gradient"]), maxdiff(g1["gradient"], r["gradient_response"])
        ok = d0 <= 1e-9 and d1 <= 1e-9
        print(f"{a.case}: E {g1['energy']:.10f} (PySCF {r['energy']:.10f}, dE {g1['energy'] - r['energy']:+.1e}); "
              f"max|g - g_pyscf| no response {d0:.1e}, with response {d1:.1e}; max|g| "
              f"{max(map(abs, r['gradient'])):.2e}; response size {maxdiff(r['gradient'], r['gradient_response']):.1e} "
              f"{'OK' if ok else 'FAIL'} (tol 1e-9)")
    elif a.check == "fd":
        mol, opts = parse_case(a.case)
        resp = "1" if a.fd == "moving" else "0"
        out = run(a.exe, ["grad", os.path.join(TC, XYZ[mol])] + opts + ["--fd", a.fd, "--response", resp]
                  + (["--h", a.h] if a.h else []))
        d = maxdiff(out["gradient"], out["fd_gradient"])
        ok = d <= a.tol
        print(f"{a.case}: analytic (grid response {'on' if resp == '1' else 'off'}) vs FD ({a.fd} grid, h "
              f"{out['fd_step_bohr']} Bohr): max diff {d:.2e}, max|g| {max(map(abs, out['gradient'])):.2e} "
              f"{'OK' if ok else 'FAIL'} (tol {a.tol:.0e})")
    elif a.check == "translation":
        mol, opts = parse_case(a.case)
        x = os.path.join(TC, XYZ[mol])
        g1 = run(a.exe, ["grad", x] + opts + ["--response", "1"])
        g0 = run(a.exe, ["grad", x] + opts + ["--response", "0"])
        rt = run(a.exe, ["rotate", x] + opts + ["--response", "1"])
        s1, s0 = max(map(abs, g1["force_sum"])), max(map(abs, g0["force_sum"]))
        ok = s1 < 1e-11 and rt["translation_max_dg"] < 1e-11
        print(f"{a.case}: |sum_A g_A| with grid response {s1:.1e} (gate 1e-11), without {s0:.1e} (grid error); "
              f"translated geometry max|dg| {rt['translation_max_dg']:.1e} {'OK' if ok else 'FAIL'}")
    elif a.check == "rotation":
        mol, opts = parse_case(a.case)
        rt = run(a.exe, ["rotate", os.path.join(TC, XYZ[mol])] + opts + ["--response", "1"])
        ok = rt["rotation_max_dg"] < 3e-4
        print(f"{a.case}: rotated (grid fixed in space): dE {rt['energy_rotated'] - rt['energy']:+.1e} Eh, max|dg| "
              f"{rt['rotation_max_dg']:.1e} Eh/Bohr (max|g| {rt['max_g']:.2e}); torque original "
              f"{max(map(abs, rt['original']['torque'])):.1e}, rotated {max(map(abs, rt['rotated']['torque'])):.1e} "
              f"{'OK' if ok else 'FAIL'} (gate 3e-4)")
    elif a.check == "threads":
        mol, opts = parse_case(a.case)
        out = run(a.exe, ["threads", os.path.join(TC, XYZ[mol])] + opts)
        ok = out["bitwise"] is True
        print(f"{a.case}: 1 vs 4 threads max|dg| {out['max_dg']:.1e} bitwise {out['bitwise']} {'OK' if ok else 'FAIL'}")
    elif a.check == "update":
        mol, opts = parse_case(a.case)
        out = run(a.exe, ["update", os.path.join(TC, XYZ[mol]), os.path.join(TC, XYZ["H2O_distorted"])] + opts)
        ok = out["max_dg"] < 1e-10 and abs(out["energy_B_updated"] - out["energy_B_fresh"]) < 1e-10
        print(f"{a.case}: UpdateMolecule vs fresh engine: dE {out['energy_B_updated'] - out['energy_B_fresh']:+.1e}, "
              f"max|dg| {out['max_dg']:.1e} {'OK' if ok else 'FAIL'}")
    elif a.check == "opt":
        r = ref["optimisations"][a.case]
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(os.path.join(TC, r["start"]), os.path.join(tmp, "start.xyz"))
            p = subprocess.run([os.path.abspath(a.curcuma), "-opt", "start.xyz", "-method", r["functional"], "-no_bmt",
                                "-verbosity", "0", "-qm.scf_threshold", "1e-9", "-opt.convergence_preset", "tight"],
                               cwd=tmp, capture_output=True, text=True, timeout=3000)
            if p.returncode != 0:
                print(p.stdout[-2000:], p.stderr[-2000:])
                sys.exit("curcuma -opt failed")
            xopt = read_xyz(os.path.join(tmp, "start.opt.xyz"))
            out = run(a.exe, ["grad", os.path.join(tmp, "start.opt.xyz"), "--functional", r["functional"], "--vwn",
                              r["vwn"], "--grid", r["grid"]])
            # the PySCF minimum must also be a stationary point of curcuma's energy (same grid, same orientation)
            with open(os.path.join(tmp, "ref.xyz"), "w") as f:
                f.write(f"{len(r['symbols'])}\nPySCF minimum\n")
                for s_, x in zip(r["symbols"], r["geometry_angstrom"]):
                    f.write(f"{s_} {x[0]:.14f} {x[1]:.14f} {x[2]:.14f}\n")
            at_ref = run(a.exe, ["grad", os.path.join(tmp, "ref.xyz"), "--functional", r["functional"], "--vwn",
                                 r["vwn"], "--grid", r["grid"]])
        dd = maxdiff(distances(xopt), distances(r["geometry_angstrom"]))
        de = out["energy"] - r["energy"]
        gmax = max(map(abs, out["gradient"]))
        de_ref, g_ref = at_ref["energy"] - r["energy"], max(map(abs, at_ref["gradient"]))
        # The optimiser may end in another orientation: the grid is not rotationally invariant (sg1, H2O: energy
        # varies by ~3e-6 Eh with the orientation), so the energy at the curcuma minimum can differ by ~1e-7 from the
        # PySCF minimum although both are minima of the same orientation-dependent surface (observed: H2O PBE 1.1e-7).
        ok = abs(de) <= 1e-6 and dd <= 2e-4 and gmax <= 1e-4 and abs(de_ref) <= 1e-9 and g_ref <= 1e-6
        print(f"{a.case}: E(opt) {out['energy']:.10f} E(PySCF min) {r['energy']:.10f} dE {de:+.1e} (gate 1e-6); max "
              f"|d_ij - d_ij(ref)| {dd:.1e} A (gate 2e-4); max|g| at the curcuma minimum {gmax:.1e} Eh/Bohr; curcuma "
              f"at the PySCF minimum: dE {de_ref:+.1e}, max|g| {g_ref:.1e} (gates 1e-9, 1e-6) {'OK' if ok else 'FAIL'}")
    else:
        sys.exit("unknown check " + a.check)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
