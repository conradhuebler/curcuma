#!/usr/bin/env python3
"""ctest gates for the double-hybrid analytic gradient (WP12 stage 2); reference.json from
scripts/qm_dh_grad_pyscf_reference.py. Gradients in Eh/Bohr; standard library only (no numpy).

  kernel        jet first derivatives (qm_xc_fxc_kernels.h) = production potentials (qm_xc_kernels.h) to 1e-11 relative;
                jet second derivatives vs libxc 7.0.0 (eval_xc deriv 2) <= 1e-8 relative and vs a 4-point FD of the
                production potentials <= 1e-5 relative after removing 10 x the stencil rounding noise (scale: max(|ref|,
                1e-3 max|row|), FD only where defined)
  fxc --case M_F   F_xc[D] vs FD of V_xc[P + t D] and the directional XC gradient vs FD of the XC gradient
                (basis, points, weight energy density): max diff <= 1e-7 max|value| (FD = test tool)
  fd --case ... --fd moving|frozen   analytic gradient vs 4-point FD of curcuma's own total energy (moving grid =
                with grid response, frozen grid = without), tol --tol
  pyscf --case C   analytic gradient with grid response vs FD of the PySCF double-hybrid energy (reference.json)
  orca --case C    250x1454 frozen-core gradient vs ORCA 6.1 EnGrad (reference_orca.json, --tol) and vs the stored
                   curcuma gradient (1e-9)
  translation / rotation / update   as qm_dft_grad; threads: 1 vs 4 vs 8 threads bitwise (RI and canonical)
  opt --case M_F   curcuma -opt -method F: final max|g| and the PySCF-minimum check (if stored)
Case names: <molecule>_<method>[_canonical][_fc][_vwn5][_default][_sg3]. Claude Generated (Oct 2026, WP12 stage 2).
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


def run(exe, args, timeout=3000):
    p = subprocess.run([exe] + args, capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0:
        print(p.stdout[-3000:], p.stderr[-3000:])
        sys.exit(f"{exe} {' '.join(args)} failed ({p.returncode})")
    start = p.stdout.find("{\n")
    return json.loads(p.stdout[start:])


def maxdiff(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def parse_case(case):
    # The driver's own defaults are the PySCF/libxc convention (PZ81-P86, no frozen core), NOT curcuma's double-hybrid
    # defaults (since Oct 31 2026: -qm.p86_lda vwn5, frozen core). Token "default" = curcuma's defaults (vwn5 + fc).
    for mol in sorted(XYZ, key=len, reverse=True):
        if case.startswith(mol + "_"):
            rest = case[len(mol) + 1:].split("_")
            break
    else:
        sys.exit(f"unknown case {case}")
    opts = ["--method", rest[0]]
    for r in rest[1:]:
        if r == "canonical":
            opts += ["--pt2", "canonical"]
        elif r == "fc":
            opts += ["--fc"]
        elif r == "vwn5":
            opts += ["--p86", "vwn5"]
        elif r == "default":
            opts += ["--p86", "vwn5", "--fc", "--cli-defaults"]
        elif r.startswith("sg") or "x" in r:
            opts += ["--grid", r]
    return mol, opts


def check_kernel(a, ref):
    K = ref["kernel"]
    out = run(a.exe, ["kernel", a.ref])
    names = {"rr": 0, "rs": 1, "rt": 2, "ss": 3, "st": 4, "tt": 5}
    ok = True
    for key in ("b88", "lyp", "pbex", "p86", "p86vwn", "b95"):
        o, r = out[key], K[key]
        nv = 3 if key == "b95" else 2
        wv = wl = wf = 0.0
        for i in range(len(K["rho"])):
            jv, pv = o["jet_v"][i], o["prod_v"][i]
            sc = max(abs(x) for x in pv[:nv])
            wv = max(wv, max(abs(jv[c] - pv[c]) / max(abs(pv[c]), 1e-3 * sc, 1e-300) for c in range(nv)))
            jh = o["jet_h"][i]
            hsc = max(abs(x) for x in jh)
            for comp, idx in names.items():
                if comp not in r["h"]:
                    continue
                ref_v = r["h"][comp][i]
                wl = max(wl, abs(jh[idx] - ref_v) / max(abs(ref_v), 1e-3 * hsc, 1e-300))
                if o["fd_h"][i] is not None:
                    f, nz = o["fd_h"][i][idx], o["fd_noise"][i][idx]
                    # minus 10 x the rounding noise of the stencil (tiny sigma: v changes in the last bits only)
                    wf = max(wf, max(0.0, abs(jh[idx] - f) - 10.0 * nz) / max(abs(jh[idx]), 1e-3 * hsc, 1e-300))
        good = wv <= 1e-11 and wl <= 1e-8 and wf <= 1e-5
        ok &= good
        print(f"{key:7s} jet v vs production {wv:.1e} (1e-11); jet H vs libxc {wl:.1e} (1e-8); jet H vs FD of the "
              f"production v {wf:.1e} (1e-5) {'OK' if good else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--case", default="")
    ap.add_argument("--fd", default="moving")
    ap.add_argument("--h", default="")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--curcuma", default="")
    ap.add_argument("--extra", default="")
    ap.add_argument("--orca", default="")
    a = ap.parse_args()
    ref = json.load(open(a.ref))
    ok = True
    if a.check == "kernel":
        ok = check_kernel(a, ref)
    elif a.check == "fxc":
        mol, opts = parse_case(a.case)
        o = run(a.exe, ["fxc", os.path.join(TC, XYZ[mol])] + opts)
        rel = {k: o[k + "_max_diff"] / max(o[k + "_max"], 1e-300) for k in ("fxc", "basis", "points", "edens")}
        ok = all(v <= 1e-7 for v in rel.values())
        print(f"{a.case}: F_xc vs FD {rel['fxc']:.1e}, directional XC gradient basis {rel['basis']:.1e}, points "
              f"{rel['points']:.1e}, weight energy density {rel['edens']:.1e} (relative, gate 1e-7) {'OK' if ok else 'FAIL'}")
    elif a.check == "fd":
        mol, opts = parse_case(a.case)
        resp = "1" if a.fd == "moving" else "0"
        out = run(a.exe, ["grad", os.path.join(TC, XYZ[mol])] + opts + ["--fd", a.fd, "--response", resp]
                  + (["--h", a.h] if a.h else []))
        d = maxdiff(out["gradient"], out["fd_gradient"])
        ok = d <= a.tol
        print(f"{a.case}: analytic (grid response {'on' if resp == '1' else 'off'}) vs FD ({a.fd} grid, h "
              f"{out['fd_step_bohr']} Bohr): max diff {d:.2e}, max|g| {max(map(abs, out['gradient'])):.2e}, Z-vector "
              f"{out['info']['zvector_iterations']} it. {'OK' if ok else 'FAIL'} (tol {a.tol:.0e})")
    elif a.check == "pyscf":
        r = ref["gradients"][a.case]
        fn, _, p86 = r["functional"].partition("@")
        opts = ["--method", fn, "--pt2", r["pt2"], "--grid", r["grid"], "--response", "1"] + (["--p86", p86] if p86 else [])
        if r["frozen_core"]:
            opts.append("--fc")
        out = run(a.exe, ["grad", os.path.join(TC, r["xyz"])] + opts)
        d = maxdiff(out["gradient"], r["gradient"])
        de = out["energy"] - r["energy"]
        # dispersion is not in the PySCF energy: compare without it
        if "dispersion" in out["parts"]:
            g = [x - y for x, y in zip(out["gradient"], out["parts"]["dispersion"])]
            d = maxdiff(g, r["gradient"])
            de -= out["info"]["decomposition"]["disp"]
        ok = d <= a.tol and abs(de) <= 1e-8
        print(f"{a.case}: E - E(PySCF) {de:+.1e} (1e-8); max|g - FD(PySCF energy)| {d:.2e} (tol {a.tol:.0e}), max|g| "
              f"{max(map(abs, r['gradient'])):.2e} {'OK' if ok else 'FAIL'}")
    elif a.check == "orca":
        r = json.load(open(a.orca))["cases"][a.case]
        args = list(r["curcuma_args"])
        args[0] = os.path.join(TC, args[0])
        out = run(a.exe, ["grad"] + args)
        d = maxdiff(out["gradient"], r["orca_gradient"])
        dc = maxdiff(out["gradient"], r["curcuma_gradient"])
        ok = d <= a.tol and dc <= 1e-9
        print(f"{a.case}: dE(ORCA 6.1) {out['energy'] - r['orca_energy']:+.1e} Eh; max|g - g_ORCA| {d:.2e} (tol "
              f"{a.tol:.0e}, grid-limited, see WP12 doc); vs the stored curcuma gradient {dc:.1e} (1e-9) "
              f"{'OK' if ok else 'FAIL'}")
    elif a.check == "translation":
        mol, opts = parse_case(a.case)
        x = os.path.join(TC, XYZ[mol])
        g1 = run(a.exe, ["grad", x] + opts + ["--response", "1"])
        g0 = run(a.exe, ["grad", x] + opts + ["--response", "0"])
        rt = run(a.exe, ["rotate", x] + opts + ["--response", "1"])
        s1, s0 = max(map(abs, g1["force_sum"])), max(map(abs, g0["force_sum"]))
        ok = s1 < 1e-10 and rt["translation_max_dg"] < 1e-10
        print(f"{a.case}: |sum_A g_A| with grid response {s1:.1e} (gate 1e-10), without {s0:.1e} (grid error); "
              f"translated geometry max|dg| {rt['translation_max_dg']:.1e} (1e-10) {'OK' if ok else 'FAIL'}")
    elif a.check == "rotation":
        mol, opts = parse_case(a.case)
        rt = run(a.exe, ["rotate", os.path.join(TC, XYZ[mol])] + opts + ["--response", "1"])
        ok = rt["rotation_max_dg"] < 3e-4
        print(f"{a.case}: rotated (grid fixed in space): dE {rt['energy_rotated'] - rt['energy']:+.1e} Eh, max|dg| "
              f"{rt['rotation_max_dg']:.1e} (max|g| {rt['max_g']:.2e}); torque original "
              f"{max(map(abs, rt['original']['torque'])):.1e}, rotated {max(map(abs, rt['rotated']['torque'])):.1e} "
              f"{'OK' if ok else 'FAIL'} (gate 3e-4)")
    elif a.check == "threads":
        mol, opts = parse_case(a.case)
        out = run(a.exe, ["threads", os.path.join(TC, XYZ[mol])] + opts)
        # since Oct 31 2026 mp2GradientHost uses fixed chunks, so RI and canonical PT2 are both gated bitwise
        ok = out["bitwise"] is True
        print(f"{a.case}: 1 vs 4 vs 8 threads max|dg| {out['max_dg']:.1e}, bitwise {out['bitwise']} (gate: bitwise) "
              f"{'OK' if ok else 'FAIL'}")
    elif a.check == "update":
        mol, opts = parse_case(a.case)
        out = run(a.exe, ["update", os.path.join(TC, XYZ[mol]), os.path.join(TC, XYZ["H2O_distorted"])] + opts)
        ok = out["max_dg"] < 1e-10 and abs(out["energy_B_updated"] - out["energy_B_fresh"]) < 1e-10
        print(f"{a.case}: updateGeometry vs fresh: dE {out['energy_B_updated'] - out['energy_B_fresh']:+.1e}, max|dg| "
              f"{out['max_dg']:.1e} {'OK' if ok else 'FAIL'}")
    elif a.check == "opt":
        mol, opts = parse_case(a.case)
        method = opts[1]
        # driver options -> curcuma flags, always explicit (the curcuma defaults differ from the driver's since Oct 31
        # 2026); "_default" cases pass no p86/frozen-core flag and check that the curcuma defaults are vwn5 + fc
        cli_defaults = "--cli-defaults" in opts
        opts = [x for x in opts if x != "--cli-defaults"]
        cli = []
        if "--pt2" in opts:
            cli += ["-qm.dh_pt2", opts[opts.index("--pt2") + 1]]
        if not cli_defaults:
            cli += ["-qm.p86_lda", opts[opts.index("--p86") + 1] if "--p86" in opts else "pz81"]
            cli += ["-qm.frozen_core", "true" if "--fc" in opts else "false"]
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(os.path.join(TC, XYZ[mol]), os.path.join(tmp, "start.xyz"))
            p = subprocess.run([os.path.abspath(a.curcuma), "-opt", "start.xyz", "-method", method, "-no_bmt",
                                "-verbosity", "0", "-qm.scf_threshold", "1e-9", "-opt.convergence_preset", "tight"]
                               + cli + a.extra.split(), cwd=tmp, capture_output=True, text=True, timeout=3000)
            if p.returncode != 0:
                print(p.stdout[-2000:], p.stderr[-2000:])
                sys.exit("curcuma -opt failed")
            out = run(a.exe, ["grad", os.path.join(tmp, "start.opt.xyz")] + opts)
            e0 = run(a.exe, ["grad", os.path.join(tmp, "start.xyz")] + opts)
        gmax = max(map(abs, out["gradient"]))
        ok = gmax <= 1e-4 and out["energy"] < e0["energy"]
        print(f"{a.case}: E(start) {e0['energy']:.10f} E(opt) {out['energy']:.10f} dE {out['energy'] - e0['energy']:+.2e}; "
              f"max|g| at the curcuma minimum {gmax:.1e} Eh/Bohr (gate 1e-4) {'OK' if ok else 'FAIL'}")
    else:
        sys.exit("unknown check " + a.check)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
