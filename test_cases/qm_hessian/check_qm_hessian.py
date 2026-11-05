#!/usr/bin/env python3
"""Native analytic RHF Hessian gate (ctest -L qm_hess). All Eh/Bohr^2, index 3*atom + xyz.

(a) analytic Hessian vs the central finite difference of curcuma's OWN analytic gradient (the fully assembled
    QMMethod::getGradient surface): proves the Hessian is the derivative of the gradient the code computes;
(b) symmetry before symmetrisation, translational invariance (sum over atoms of every row block) and rotational
    invariance H (L r) = L g, which holds at every geometry;
(c) analytic Hessian vs PySCF (reference.json, scripts/hf_hessian_reference.py), an independent implementation:
    vs the FD of PySCF's analytic gradient (tight) and vs PySCF's analytic Hessian (looser: PySCF's own analytic
    Hessian is up to 9.0e-8 off its FD on these cases);
(d) --pieces: derivative integrals vs finite differences at fixed density (see dump_qm_hessian.cpp).
Usage: check_qm_hessian.py --dump <dump_qm_hessian> --ref reference.json --case NAME [--pieces] [--no-fd]
Claude Generated (Nov 2026).
"""
import argparse, json, os, subprocess, sys

TOL_REF = 2e-7     # vs PySCF analytic Hessian (observed <= 9.0e-8, PySCF vs its own FD alike)
TOL_REF_FD = 2e-8  # vs FD of PySCF's analytic gradient (4-point, h = 5e-4 Bohr; observed <= 9.1e-9)
TOL_FD = 1e-8      # vs FD of the own analytic gradient (h = 5e-4 Bohr, 4-point; observed <= 2.0e-9)
TOL_SYM = 1e-9     # max |H_xy - H_yx| before symmetrisation (observed <= 2.8e-11)
TOL_INV = 1e-9     # translational / rotational residual (observed <= 4.4e-11)
TOL_PIECE = 1e-7   # relative to the largest element, FD step 1e-4 Bohr (truncation ~h^2; observed <= 7.4e-9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--pieces", action="store_true")
    ap.add_argument("--no-fd", action="store_true")
    ap.add_argument("--qm", nargs=2, action="append", default=[], metavar=("KEY", "VALUE"),
                    help="extra qm-module setting for dump_qm_hessian, e.g. --qm scf_direct on")
    a = ap.parse_args()
    ref = json.load(open(a.ref))["cases"][a.case]
    xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", ref["xyz"])
    basis = {"def2-svp": "def2-SVP", "def2-tzvp": "def2-TZVP"}[ref["basis"]]
    cmd = [os.path.abspath(a.dump), xyz, "--basis", basis, "--threads", "4"]
    for key, val in a.qm:
        cmd += ["--qm", key, val]
    if not a.no_fd:
        cmd += ["--fd", "--stencil", "4", "--h", "5e-4"]
    if a.pieces:
        cmd.append("--pieces")
    got = json.loads(subprocess.run(cmd, capture_output=True, text=True, timeout=1800, check=True).stdout)

    H, Hr = got["hessian"], ref["hessian"]
    d_ref = max(abs(x - y) for x, y in zip(H, Hr))
    d_e = abs(got["energy"] - ref["energy"])
    info = got["info"]
    d_ref_fd = max(abs(x - y) for x, y in zip(H, ref["hessian_fd"]))
    ok = d_ref <= TOL_REF and d_ref_fd <= TOL_REF_FD and d_e <= 1e-8
    print(f"{a.case:28s} energy diff {d_e:.1e}  max|H - PySCF analytic| = {d_ref:.2e}  "
          f"max|H - FD(PySCF gradient)| = {d_ref_fd:.2e}  (max|H| {max(map(abs, Hr)):.3e})")
    print(f"{a.case:28s} asymmetry {info['max_asymmetry']:.2e}  translation {got['translation_residual']:.2e}  "
          f"rotation {got['rotation_residual']:.2e}  |g| {got['gradient_norm']:.2e}  "
          f"CPHF {info['cphf_iterations']} it, {info['jk_builds']} J/K, {got['hessian_ms']:.0f} ms")
    ok &= info["max_asymmetry"] <= TOL_SYM
    ok &= got["translation_residual"] <= TOL_INV and got["rotation_residual"] <= TOL_INV
    if "fd_hessian" in got:
        print(f"{a.case:28s} max|H - FD(own gradient)| = {got['max_abs_diff_fd']:.2e}  "
              f"(symmetrised FD {got['max_abs_diff_fd_sym']:.2e}, FD {got['fd_ms']:.0f} ms)")
        ok &= got["max_abs_diff_fd"] <= TOL_FD
    if "pieces" in got:
        p = got["pieces"]
        print(f"{a.case:28s} digest J-K/2 vs tensor {p['digest_vs_tensor_JK']:.2e}")
        ok &= p["digest_vs_tensor_JK"] <= 1e-11
        for key in ("S1", "H1", "G1", "hess1e", "hess2e"):
            rel = p[key + "_vs_fd"] / max(p[key + "_max"], 1e-30)
            print(f"{a.case:28s} {key:7s} vs FD: max diff {p[key + '_vs_fd']:.2e} (max {p[key + '_max']:.2e}, rel {rel:.1e})")
            ok &= rel <= TOL_PIECE
    print("ok" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
