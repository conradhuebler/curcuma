#!/usr/bin/env python3
"""ctest checker for the native B3LYP (WP7), derived from check_qm_pbe.py. Standard library only.

Checks (--check):
  kernel      qmxc B88 and LYP point kernels vs libxc (reference.json "kernel", PySCF 2.14 / libxc 7.0.0) on 472
              (rho, sigma) pairs: eps, v_rho, v_sigma to KERNEL_TOL in a SCALED error |a - b| / scale:
                eps, v_rho: scale max(|b|, |eps_LDA-part|) with the LDA part of the kernel at that rho (Slater for B88,
                the sigma-free LYP terms for LYP); the gradient term cancels against it at large sigma
                v_sigma: max(|b|, rho |eps_ref| / sigma), eps_ref = eps of the kernel at sigma = 0
              Points libxc drops (rho below its dens_threshold, eps = 0 there) are skipped and counted.
              v_rho and v_sigma vs central finite differences of rho eps: FD_TOL in the same scaled error.
              No NaN/inf on a dense (rho 1e-20..1e2, sigma 0..1e20) scan.
  vxcfd       V_xc vs a finite difference of E_xc in the density matrix (VXCFD_TOL relative); 1 vs 4 threads equal
  energy      B3LYP/def2-SVP SCF energy and E_xc vs PySCF RKS on the identical grid (--case: VWN5 = libxc B3LYP5,
              VWN-RPA = libxc B3LYP), ENERGY_TOL Eh; electrons on the grid to 1e-8; RI cases against PySCF
              density_fit with the same aux file (rij: only_dfj = J fitted, K exact; rijk: J and K fitted);
              the exact-exchange part -(a_x/4) Tr(PK) vs PySCF's K on its density (ENERGY_TOL), a_x = 0.2
  gradient    a gradient request for b3lyp succeeds since WP8 DFT (Oct 2026): gradient_failed false, hasGradient() true;
              the values are checked in test_cases/qm_dft_grad
  invariance  atom order reversed and translation: |dE| <= 1e-9 Eh; rotation (grid fixed in space) <= ROT_TOL (sg1)
  update      UpdateMolecule A -> B equals a fresh engine on B (1e-9 Eh)
  threads     1 and 4 threads give the same energy (1e-10 Eh)
  direct      integral-direct J/K equals the stored-tensor J/K path (1e-8 Eh)

Claude Generated (Oct 2026, WP7).
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

KERNEL_TOL = 1e-12
FD_TOL = 1e-6
VXCFD_TOL = 1e-7
ENERGY_TOL = 1e-8
ROT_TOL = 1e-4


def run(exe, args):
    p = subprocess.run([exe] + args, capture_output=True, text=True)
    if p.returncode not in (0, 2):
        print(p.stdout, p.stderr)
        raise SystemExit(f"test_qm_b3lyp failed ({p.returncode}): {' '.join(args)}")
    return json.loads(p.stdout)


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


def lda_parts(rho):
    """|eps| of the sigma-free part of B88 (Slater per spin) and LYP (the F + C_F G terms) at rho."""
    cx = 1.5 * (3.0 / (4.0 * math.pi)) ** (1.0 / 3.0)
    b88 = 2.0 ** (-1.0 / 3.0) * rho ** (1.0 / 3.0) * cx
    a, bb, c, d = 0.04918, 0.132, 0.2533, 0.349
    cf = 0.3 * (3.0 * math.pi ** 2) ** (2.0 / 3.0)
    r = rho ** (-1.0 / 3.0)
    lyp = abs(-a / (1 + d * r) - a * bb * math.exp(-c * r) / (1 + d * r) * cf)
    return {"b88": b88, "lyp": lyp}


def kernel_scale(key, q, b, i, k):
    rho, sigma = k["rho"][i], k["sigma"][i]
    e0 = lda_parts(rho)[key]
    if q in ("eps", "vrho"):
        return max(abs(b), e0)
    return max(abs(b), rho * e0 / sigma if sigma > 0 else 0.0)


def check_kernel(a):
    ref = json.load(open(a.ref))
    k = ref["kernel"]
    out = run(a.exe, ["kernel", a.ref])
    n = len(k["rho"])
    ok = True
    print(f"{n} (rho, sigma) pairs, rho {min(k['rho']):.0e}..{max(k['rho']):.0e}, sigma 0..{max(k['sigma']):.1e}; "
          f"libxc sigma floors {k['sigma_floor']}, libxc density thresholds {k['dens_threshold']}")
    for key in ("b88", "lyp"):
        skipped = sum(1 for i in range(n) if k[key]["eps"][i] == 0.0)
        for q in ("eps", "vrho", "vsigma"):
            worst, worst_i, plain = 0.0, 0, 0.0
            for i in range(n):
                if k[key]["eps"][i] == 0.0:
                    continue  # dropped by libxc
                b, c = k[key][q][i], out[key][q][i]
                err = abs(c - b) / max(kernel_scale(key, q, b, i, k), 1e-300)
                if err > worst:
                    worst, worst_i = err, i
                if b != 0.0:
                    plain = max(plain, abs(c - b) / abs(b))
            good = worst <= KERNEL_TOL
            ok &= good
            print(f"{key:4s} {q:6s} vs {k[key]['libxc']:10s}: max scaled err {worst:.2e} (rho {k['rho'][worst_i]:.1e}, "
                  f"sigma {k['sigma'][worst_i]:.1e}); max plain rel {plain:.1e}; {skipped} pairs dropped by libxc "
                  f"{'OK' if good else 'FAIL'}")
        for q in ("vrho", "vsigma"):
            worst = 0.0
            for i in range(n):
                f = out[key][q + "_fd"][i]
                if f is None:
                    continue
                c = out[key][q][i]
                worst = max(worst, abs(c - f) / max(kernel_scale(key, q, c, i, k), 1e-300))
            good = worst <= FD_TOL
            ok &= good
            print(f"{key:4s} {q:6s} vs FD of rho eps: max scaled err {worst:.2e} (tol {FD_TOL:.0e}) "
                  f"{'OK' if good else 'FAIL'}")
        d = max((abs(out[key]["eps_raw"][i] - k[key]["eps"][i]) / max(abs(k[key]["eps"][i]), 1e-300)
                 for i in range(n) if k["sigma"][i] == 0.0 and k[key]["eps"][i] != 0.0), default=0.0)
        print(f"{key}: sigma = 0 without libxc's floor vs libxc: max rel eps {d:.1e} (diagnostic)")
    good = out["nonfinite"] == 0
    ok &= good
    print(f"non-finite eps/v_rho/v_sigma on {out['nonfinite_scanned']} scanned (rho, sigma, kernel) points: "
          f"{out['nonfinite']} {'OK' if good else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--case", default="H2O_sg1")
    a = ap.parse_args()
    ref = json.load(open(a.ref))
    tc = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..")
    xyzs = {"H2O": "qm_1e/H2O.xyz", "NH3": "qm_1e/NH3.xyz", "CH4": "qm_1e/CH4.xyz", "benzene": "qm_hf3c/benzene.xyz"}
    ok = True

    if a.check == "kernel":
        ok = check_kernel(a)
    elif a.check == "vxcfd":
        out = run(a.exe, ["vxcfd", os.path.join(tc, xyzs[a.case])])
        good = out["rel"] <= VXCFD_TOL and out["threads_maxdiff"] == 0.0 and out["exc_threads_diff"] == 0.0 \
            and out["asym"] == 0.0
        ok &= good
        print(f"{a.case}: Tr(V_xc D) {out['trVD']:.12f}, FD of E_xc {out['fd']:.12f}, rel {out['rel']:.2e} "
              f"(tol {VXCFD_TOL:.0e}); V_xc 1 vs 4 threads max diff {out['threads_maxdiff']:.1e}, E_xc diff "
              f"{out['exc_threads_diff']:.1e}; asym {out['asym']:.1e}; {out['npts']} points {'OK' if good else 'FAIL'}")
    elif a.check == "energy":
        r = ref["energies"][a.case]
        out = run(a.exe, ["energy", os.path.join(tc, r["xyz"]), "--grid", r["grid"], "--ri", r.get("ri", "off"),
                          "--vwn", r["vwn"]])
        dE, dX, dN = out["energy"] - r["energy"], out["exc"] - r["exc"], out["nelec_grid"] - r["nelec_grid"]
        dK = out["exact_exchange"] - r["exx"]
        good = out["converged"] and abs(dE) <= ENERGY_TOL and abs(dX) <= ENERGY_TOL and abs(dN) <= 1e-8 \
            and abs(dK) <= ENERGY_TOL and out["a_x"] == r["hybrid_coeff"]
        ok &= good
        print(f"{a.case}: E {out['energy']:.10f} (PySCF {r['energy']:.10f}, dE {dE:.2e}); E_xc dE {dX:.2e}; "
              f"E_HFX {out['exact_exchange']:.10f} (dE {dK:.2e}, a_x {out['a_x']}); "
              f"N_grid {out['nelec_grid']:.8f} (dN {dN:.1e}); {out['iterations']} it; PySCF {r['xc']}; {out['label']} "
              f"{'OK' if good else 'FAIL'} (tol {ENERGY_TOL:.0e})")
    elif a.check == "gradient":
        out = run(a.exe, ["energy", os.path.join(tc, "qm_1e/H2O.xyz"), "--gradient"])
        good = out.get("gradient_failed") is False and out["has_gradient"] is True  # WP8 DFT (Oct 2026)
        ok &= good
        print(f"gradient request: failed={out.get('gradient_failed')} has_gradient={out['has_gradient']} "
              f"'{out.get('gradient_error')}' {'OK' if good else 'FAIL'}")
    elif a.check in ("invariance", "update", "threads", "direct", "rij"):
        mol = a.case.split("_")[0]
        xyz = os.path.join(tc, xyzs[mol])
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
