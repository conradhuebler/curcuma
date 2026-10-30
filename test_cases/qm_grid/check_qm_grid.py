#!/usr/bin/env python3
"""ctest checker for the native DFT grid (WP4); standard library only.

  --check lebedev   Y_lm integrals of every Lebedev grid (test_qm_grid lebedev)
  --check points    grid of --case vs PySCF (reference.json grid_ref)
  --check density   HF density integrated on the grid of --case vs PySCF (reference.json density_ref)
  --check aofd      AO gradients vs finite differences (molecule and grid of --case)

Tolerances (why): points/weights 1e-12 relative (same formulas, only the floating-point order differs);
electron count 1e-8 (both SCFs converged to ~1e-9 in the density); rho 1e-7 relative to max(rho, 1e-3), grad rho
1e-7 relative to max(rho, |grad rho|, 1e-3); weights: partition factor w/w_quad 1e-12 absolute; Lebedev 1e-12 (16-digit tables, sums of up to 5810 terms); AO FD 1e-6 relative to max(1, |grad|).
Claude Generated (Oct 2026).
"""
import argparse
import json
import os
import subprocess
import sys


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED:", " ".join(cmd), "\n", r.stderr)
        sys.exit(1)
    return json.loads(r.stdout)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--check", required=True, choices=["lebedev", "points", "density", "aofd"])
    ap.add_argument("--case")
    ap.add_argument("--cartesian", action="store_true")
    a = ap.parse_args()
    ok = True

    if a.check == "lebedev":
        res = run([a.exe, "lebedev"])
        print("%6s %6s %12s %12s %12s %12s" % ("npts", "degree", "|w|-1", "err l<=deg", "err l=deg+1", "ortho"))
        for r in res:
            o = r.get("max_err_orthonormality", float("nan"))
            print("%6d %6d %12.2e %12.2e %12.2e %12.2e" % (r["npts"], r["degree"], r["weight_sum_minus_1"],
                  r["max_err_l_le_degree"], r["max_err_l_degree_plus_1"], o))
            if abs(r["weight_sum_minus_1"]) > 1e-13 or r["max_err_l_le_degree"] > 1e-12:
                ok = False
            if r["max_err_l_degree_plus_1"] < 1e-6:  # exact one degree higher would mean a wrong degree label
                ok = False
            if "max_err_orthonormality" in r and o > 1e-12:
                ok = False
        print("PASS" if ok else "FAIL")
        return 0 if ok else 1

    ref = json.load(open(a.ref))
    case = ref["cases"][a.case]
    xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", case["xyz"])
    opts = ["--grid", case["grid"], "--radial", case["radial"], "--partition", case["partition"],
            "--adjust", case["radii_adjust"]]

    if a.check == "points":
        g = case["grid_ref"]
        out = run([a.exe, "grid", xyz] + opts + ["--stride", str(g["stride"])])
        if out["npts"] != g["npts"]:
            print("npts", out["npts"], "!=", g["npts"])
            return 1
        for k in ("weight_sum", "gauss_integral", "pz2_integral"):
            d = rel(out[k], g[k])
            print("%-15s curcuma %.15e pyscf %.15e rel %.1e" % (k, out[k], g[k], d))
            ok &= d < 1e-12
        print("gauss vs exact  rel %.1e   pz2 vs exact rel %.1e" % (rel(out["gauss_integral"], out["gauss_exact"]),
              rel(out["pz2_integral"], out["pz2_exact"])))
        for A, (x, y) in enumerate(zip(out["atom_weight_sums"], g["atom_weight_sums"])):
            ok &= rel(x, y) < 1e-12
        dx = dw = dwrel = 0.0
        for pc, pr in zip(out["points"], g["points"]):
            if pc[0] != pr[0]:
                print("index mismatch", pc[0], pr[0])
                return 1
            for k in range(3):
                dx = max(dx, abs(pc[4 + k] - pr[1 + k]) / max(1.0, abs(pr[1 + k])))
            # Gate: error of the partition factor w / w_quad (in [0, 1]), i.e. |dw| / w_quad. The plain relative
            # error is printed too: small factors (outer shells, 1 - g with g -> 1) lose relative digits to cancellation.
            dw = max(dw, abs(pc[7] - pr[4]) / pc[8])
            dwrel = max(dwrel, abs(pc[7] - pr[4]) / max(abs(pr[4]), 1e-300))
        print("sampled points %d: max coord dev %.1e, max |dw|/w_quad %.1e (max |dw|/w %.1e, not gated)" % (
            len(g["points"]), dx, dw, dwrel))
        ok &= dx < 1e-12 and dw < 1e-12
    elif a.check == "density":
        d = case["density_ref"]
        out = run([a.exe, "density", xyz] + opts + ["--stride", str(case["grid_ref"]["stride"])])
        print("E_HF curcuma %.10f pyscf %.10f" % (out["energy"], d["energy"]))
        print("N_el %d: curcuma int rho = %.12f (dev %.3e), pyscf %.12f (dev %.3e), diff %.1e" % (
            d["num_electrons"], out["nelec_grid"], out["nelec_grid"] - d["num_electrons"], d["nelec_grid"],
            d["deviation"], abs(out["nelec_grid"] - d["nelec_grid"])))
        ok &= abs(out["nelec_grid"] - d["nelec_grid"]) < 1e-8
        dr = dg = 0.0
        for pc, pr in zip(out["points"], d["points"]):
            if pc[0] != pr[0]:
                print("index mismatch")
                return 1
            # scale of the point: max(rho, |grad rho|, 1e-3); near a nucleus grad rho is a cancellation of terms of
            # the size of rho, so its error is measured against rho there
            scale = max(abs(pr[1]), max(abs(v) for v in pr[2:5]), 1e-3)
            dr = max(dr, abs(pc[1] - pr[1]) / max(abs(pr[1]), 1e-3))
            for k in range(3):
                dg = max(dg, abs(pc[2 + k] - pr[2 + k]) / scale)
        print("sampled points %d: rho max rel dev %.1e, grad rho max rel dev %.1e" % (len(d["points"]), dr, dg))
        ok &= dr < 1e-7 and dg < 1e-7
    else:
        cmd = [a.exe, "aofd", xyz] + (["--cartesian"] if a.cartesian else [])
        out = run(cmd)
        print("AO gradient vs FD: %d points, nbf %d, max rel err %.2e (max |grad| %.2e)" % (
            out["points_checked"], out["nbf"], out["max_rel_err"], out["max_abs_grad"]))
        ok &= out["max_rel_err"] < 1e-6
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
