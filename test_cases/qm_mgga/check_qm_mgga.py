#!/usr/bin/env python3
"""ctest checker for the native meta-GGAs TPSS and SCAN (WP11 stage 1). Standard library only.

Checks (--check):
  kernel      qmxc TPSS-X, TPSS-C, SCAN-X, SCAN-C vs libxc (reference.json "kernel", PySCF 2.14 / libxc 7.0.0) on the
              (rho, sigma, tau) table: eps, v_rho, v_sigma, v_tau to KERNEL_TOL in a SCALED error |a - b| / scale:
                eps:     exchange max(|b|, |eps_x^unif|), correlation max(|b|, |eps_c^PW92|)
                v_rho:   exchange max(|b|, |eps_x^unif|), correlation max(|b|, |eps_c^PW92|, |v_c^PW92|)
                v_sigma: max(|b|, rho |eps_ref| / sigma),  v_tau: max(|b|, rho |eps_ref| / tau)
              (eps_ref = eps_x^unif or eps_c^PW92). The scale is the size of the terms that cancel: TPSS correlation
              vanishes for z = tau_W / tau = 1 by cancellation, v_sigma sigma and v_tau tau are energy scales.
              Plain relative errors are printed too. TPSS correlation is not differentiable at z = 1 (min(z, 1),
              TPSS 2003 eq. 11): at alpha = 0 and at the "clamp" points (tau < tau_W, sigma lowered to 8 rho tau as
              libxc does, so z = 1 again) the side of the kink is decided by rounding, so its v_rho, v_sigma, v_tau
              are left out there (eps is compared). All other kernels are checked on the clamp points too.
              v_rho, v_sigma, v_tau vs central finite differences of rho eps (alpha >= 0.01, tau >= 1e-18, rho >= 1e-10,
              admissible points; points
              within 1e-4 of a SCAN cut-off are skipped, the steps keep alpha within 1e-5 max(alpha, 0.01)):
              FD_TOL in the same scaled error.
  vxcfd       V_xc vs a finite difference of E_xc in the density matrix (VXCFD_TOL relative); V_xc with 1 and 4
              threads identical (0), V_xc symmetric; --case <mol>_<functional>[_screen]
  energy      SCF energy and E_xc vs PySCF RKS on the identical grid (--case), ENERGY_TOL Eh; electrons on the grid
              vs PySCF to 1e-8; RI cases against PySCF density_fit with the same aux file
  gradient    a gradient request fails loudly (no meta-GGA gradient in stage 1): gradient_failed true,
              hasGradient() false, error text names the meta-GGA
  invariance  atom order reversed and translation: |dE| <= 1e-9 Eh; rotation (the Lebedev grid stays fixed in space):
              |dE| reported, limit ROT_TOL at the case's grid
  update      UpdateMolecule A -> B equals a fresh engine on B (1e-9 Eh)
  threads     1 and 4 threads give the same energy (1e-10 Eh)
  screen      AO screening off (-qm.xc_ao_threshold 0) vs default 1e-12: |dE| <= 1e-9 Eh
  gpu         curcuma -sp ... -gpu cuda: the meta-GGA XC part falls back to the CPU as a counted fallback (summary
              line "qm_xc"), the energy equals the CPU run to 1e-8 Eh; with -gpu_strict true the run stops (exit 3)

Claude Generated (Oct 2026, WP11 stage 1).
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile

KERNEL_TOL = 1e-10
FD_TOL = 1e-5     # secondary check (libxc is the gate); worst measured 1.8e-6: SCAN-C v_rho at rho 1e-10, alpha 1 (eps_c cancels at large t, the FD is rounding-limited), all others <= 3e-7
VXCFD_TOL = 1e-7
ENERGY_TOL = 1e-8
ROT_TOL = {"sg1": 1e-4, "sg3": 1e-5}
EXCHANGE = ("tpss_x", "scan_x")


def run(exe, args):
    p = subprocess.run([exe] + args, capture_output=True, text=True)
    if p.returncode not in (0, 2):
        print(p.stdout, p.stderr)
        raise SystemExit(f"test_qm_mgga failed ({p.returncode}): {' '.join(args)}")
    out = p.stdout
    start = out.find("{\n")  # error lines of CurcumaLogger (with colour codes) may precede the JSON
    return json.loads(out[max(start, 0):])


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


def eps_ref(key, k, i):
    if key in EXCHANGE:
        return 0.75 * (3.0 * k["rho"][i] / math.pi) ** (1.0 / 3.0)
    return abs(k["pw92"]["eps"][i])


def kernel_scale(key, q, b, i, k):
    rho, sigma, tau = k["rho"][i], max(k["sigma"][i], k["sigma_floor"]), max(k["tau"][i], 1e-20)
    e = eps_ref(key, k, i)
    if q == "eps":
        return max(abs(b), e)
    if q == "vrho":
        return max(abs(b), e) if key in EXCHANGE else max(abs(b), e, abs(k["pw92"]["vrho"][i]))
    if q == "vsigma":
        return max(abs(b), rho * e / sigma)
    return max(abs(b), rho * e / tau)


def near_cutoff(k, i):
    a = k["alpha"][i]
    return any(abs(a - c) <= 1e-4 * c for c in k["scan_cutoffs"])


def check_kernel(a):
    ref = json.load(open(a.ref))
    k = ref["kernel"]
    out = run(a.exe, ["kernel", a.ref])
    n = len(k["rho"])
    ok = True
    nclamp = sum(1 for x in k["kind"] if x == "clamp")
    print(f"{n} (rho, sigma, tau) points ({nclamp} clamp points tau < tau_W), rho {min(k['rho']):.0e}.."
          f"{max(k['rho']):.0e}, alpha 0..{max(k['alpha']):.0e}; sigma floor {k['sigma_floor']:.0e}")
    for key in ("tpss_x", "tpss_c", "scan_x", "scan_c"):
        for q in ("eps", "vrho", "vsigma", "vtau"):
            worst, wi, plain, skipped = 0.0, 0, 0.0, 0
            for i in range(n):
                if key == "tpss_c" and q != "eps" and k["alpha"][i] <= 0.0:  # z = 1: alpha = 0 and clamp points
                    skipped += 1
                    continue
                b, c = k[key][q][i], out[key][q][i]
                err = abs(c - b) / max(kernel_scale(key, q, b, i, k), 1e-300)
                if not err <= worst and not err < worst:  # NaN or larger
                    worst, wi = (err if err == err else float("inf")), i
                if b != 0.0:
                    plain = max(plain, abs(c - b) / abs(b))
            good = worst <= KERNEL_TOL
            ok &= good
            print(f"{key:6s} {q:6s} vs {k[key]['libxc']:11s}: max scaled err {worst:.2e} (rho {k['rho'][wi]:.1e}, "
                  f"sigma {k['sigma'][wi]:.1e}, alpha {k['alpha'][wi]:.3g}); max plain rel {plain:.1e}"
                  f"{f'; {skipped} kink points left out' if skipped else ''} {'OK' if good else 'FAIL'}")
        for q in ("vrho", "vsigma", "vtau"):
            worst, wi, cnt, cut = 0.0, 0, 0, 0
            for i in range(n):
                f = out[key][q + "_fd"][i]
                if f is None:
                    continue
                if near_cutoff(k, i):
                    cut += 1
                    continue
                cnt += 1
                c = out[key][q][i]
                err = abs(c - f) / max(kernel_scale(key, q, c, i, k), 1e-300)
                if not err <= worst:
                    worst, wi = (err if err == err else float("inf")), i
            good = worst <= FD_TOL
            ok &= good
            print(f"{key:6s} {q:6s} vs FD of rho eps: max scaled err {worst:.2e} over {cnt} points (rho "
                  f"{k['rho'][wi]:.1e}, alpha {k['alpha'][wi]:.3g}; {cut} next to a SCAN cut-off not used) "
                  f"(tol {FD_TOL:.0e}) {'OK' if good else 'FAIL'}")
    return ok


def cli_energy(curcuma, xyz, functional, extra):
    p = subprocess.run([curcuma, "-sp", xyz, "-method", functional, "-no_bmt", "-verbosity", "1",
                        "-qm.scf_threshold", "1e-10"] + extra, capture_output=True, text=True)
    txt = p.stdout + p.stderr
    m = re.findall(r"total [A-Z]+ energy[^0-9-]*(-?\d+\.\d+)", txt)
    return p.returncode, (float(m[-1]) if m else None), txt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--case", default="H2O_tpss_sg1")
    ap.add_argument("--curcuma", default="")
    a = ap.parse_args()
    ref = json.load(open(a.ref))
    tc = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..")
    xyzs = {"H2O": "qm_1e/H2O.xyz", "NH3": "qm_1e/NH3.xyz", "CH4": "qm_1e/CH4.xyz", "benzene": "qm_hf3c/benzene.xyz"}
    ok = True

    if a.check == "kernel":
        ok = check_kernel(a)
    elif a.check == "vxcfd":
        parts = a.case.split("_")
        mol, fn = parts[0], parts[1]
        thr = "1e-12" if len(parts) > 2 and parts[2] == "screen" else "0"
        out = run(a.exe, ["vxcfd", os.path.join(tc, xyzs[mol]), "--functional", fn, "--ao-threshold", thr])
        good = out["rel"] <= VXCFD_TOL and out["threads_maxdiff"] == 0.0 and out["exc_threads_diff"] == 0.0 \
            and out["asym"] == 0.0
        ok &= good
        print(f"{a.case}: Tr(V_xc D) {out['trVD']:.12f}, FD of E_xc {out['fd']:.12f}, rel {out['rel']:.2e} "
              f"(tol {VXCFD_TOL:.0e}); V_xc 1 vs 4 threads max diff {out['threads_maxdiff']:.1e}, E_xc diff "
              f"{out['exc_threads_diff']:.1e}; asym {out['asym']:.1e}; {out['npts']} points, screened "
              f"{out['screened']} {'OK' if good else 'FAIL'}")
    elif a.check == "energy":
        r = ref["energies"][a.case]
        out = run(a.exe, ["energy", os.path.join(tc, r["xyz"]), "--functional", r["functional"], "--grid", r["grid"],
                          "--ri", r.get("ri", "off")])
        dE, dX, dN = out["energy"] - r["energy"], out["exc"] - r["exc"], out["nelec_grid"] - r["nelec_grid"]
        good = out["converged"] and abs(dE) <= ENERGY_TOL and abs(dX) <= ENERGY_TOL and abs(dN) <= 1e-8
        ok &= good
        print(f"{a.case}: E {out['energy']:.10f} (PySCF {r['energy']:.10f}, dE {dE:.2e}); E_xc dE {dX:.2e}; "
              f"N_grid {out['nelec_grid']:.8f} (dN {dN:.1e}); {out['iterations']} it; {out['label']} "
              f"{'OK' if good else 'FAIL'} (tol {ENERGY_TOL:.0e})")
    elif a.check == "gradient":
        fn = a.case.split("_")[1]
        out = run(a.exe, ["energy", os.path.join(tc, "qm_1e/H2O.xyz"), "--functional", fn, "--gradient"])
        good = out.get("gradient_failed") is True and out["has_gradient"] is False \
            and "meta-GGA" in out.get("gradient_error", "") and out["converged"]
        ok &= good
        print(f"{fn} gradient request: failed={out.get('gradient_failed')} has_gradient={out['has_gradient']} "
              f"'{out.get('gradient_error')}' {'OK' if good else 'FAIL'}")
    elif a.check == "gpu":
        mol, fn = a.case.split("_")[:2]
        xyz = os.path.join(tc, xyzs[mol])
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)
            rc0, e_cpu, _ = cli_energy(a.curcuma, xyz, fn, ["-gpu", "none"])
            rc1, e_gpu, txt = cli_energy(a.curcuma, xyz, fn, ["-gpu", "cuda"])
            rc2, e_strict, _ = cli_energy(a.curcuma, xyz, fn, ["-gpu", "cuda", "-gpu_strict", "true"])
        counted = re.search(r"\d+ x qm_xc - first: .*meta-GGA", txt) is not None
        good = rc0 == 0 and rc1 == 0 and e_cpu is not None and e_gpu is not None and abs(e_gpu - e_cpu) <= 1e-8 \
            and counted and rc2 == 3
        ok &= good
        print(f"{mol} {fn}: CPU {e_cpu} (exit {rc0}), -gpu cuda {e_gpu} (exit {rc1}, counted qm_xc fallback "
              f"{counted}), diff {None if e_cpu is None or e_gpu is None else f'{e_gpu - e_cpu:.2e}'}; "
              f"-gpu_strict true exit {rc2} (expected 3) {'OK' if good else 'FAIL'}")
    elif a.check in ("invariance", "update", "threads", "screen"):
        parts = a.case.split("_")
        mol, fn = parts[0], parts[1]
        grid = parts[2] if len(parts) > 2 else "sg1"
        xyz = os.path.join(tc, xyzs[mol])
        opts = ["--functional", fn, "--grid", grid]
        base = run(a.exe, ["energy", xyz] + opts)["energy"]
        with tempfile.TemporaryDirectory() as td:
            atoms = read_xyz(xyz)
            if a.check == "invariance":
                tests = {"reversed": (list(reversed(atoms)), 1e-9),
                         "translated": ([[s, x + 1.3, y - 0.7, z + 2.1] for s, x, y, z in atoms], 1e-9),
                         "rotated": (rotate(atoms, (1.0, 2.0, 3.0), 0.7), ROT_TOL[grid])}
                for name, (at, tol) in tests.items():
                    p = os.path.join(td, name + ".xyz")
                    write_xyz(p, at)
                    e = run(a.exe, ["energy", p] + opts)["energy"]
                    good = abs(e - base) <= tol
                    ok &= good
                    print(f"{mol} {fn} {grid} {name}: E {e:.10f}, dE {e - base:.2e} (tol {tol:.0e}) "
                          f"{'OK' if good else 'FAIL'}")
            elif a.check == "update":
                p = os.path.join(td, "moved.xyz")
                moved = [[s, x * 1.03, y * 0.98, z + 0.02 * i] for i, (s, x, y, z) in enumerate(atoms)]
                write_xyz(p, moved)
                out = run(a.exe, ["update", xyz, p] + opts)
                d = out["energy_B_updated"] - out["energy_B_fresh"]
                good = out["converged"] and abs(d) <= 1e-9 and abs(out["energy_A"] - out["energy_B_fresh"]) > 1e-6
                ok &= good
                print(f"{mol} {fn} update: E_A {out['energy_A']:.10f}, E_B updated {out['energy_B_updated']:.10f}, "
                      f"fresh {out['energy_B_fresh']:.10f}, diff {d:.2e} (warm start {out['warm_start']}) "
                      f"{'OK' if good else 'FAIL'}")
            elif a.check == "threads":
                e1 = run(a.exe, ["energy", xyz, "--threads", "1"] + opts)["energy"]
                good = abs(e1 - base) <= 1e-10
                ok &= good
                print(f"{mol} {fn} threads 1 vs 4: {e1:.12f} vs {base:.12f}, diff {e1 - base:.2e} "
                      f"{'OK' if good else 'FAIL'}")
            else:
                e0 = run(a.exe, ["energy", xyz, "--ao-threshold", "0"] + opts)["energy"]
                good = abs(e0 - base) <= 1e-9
                ok &= good
                print(f"{mol} {fn} AO screening off vs 1e-12: {e0:.12f} vs {base:.12f}, diff {e0 - base:.2e} "
                      f"{'OK' if good else 'FAIL'}")
    else:
        raise SystemExit("unknown check " + a.check)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
