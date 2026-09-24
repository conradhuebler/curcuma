#!/usr/bin/env python3
"""Native HF / HF-3c analytic-gradient gate (ctest -L qm_grad).

(a) analytic gradient vs the central finite difference of curcuma's OWN energy
    (proves the analytic expression is the derivative of what the code computes);
(b) analytic gradient vs PySCF's analytic RHF gradient (+ simple-dftd3 D3/gCP
    gradients for hf-3c) from reference.json -- an independent implementation.
All Eh/Bohr. Usage: check_qm_gradient.py --dump <dump_qm_gradient> --ref reference.json --case NAME
Claude Generated (Sep 2026).
"""
import argparse, json, os, subprocess, sys

TOL_REF = 1e-8   # vs PySCF / simple-dftd3 (observed <= 1.5e-10)
TOL_FD = 1e-7    # vs FD (h = 1e-4 Bohr: truncation ~h^2 g''' ~ 3e-9, plus SCF noise)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--no-fd", action="store_true")
    a = ap.parse_args()
    ref = json.load(open(a.ref))["cases"][a.case]
    xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", ref["xyz"])
    basis = "MINIX" if ref["method"] == "hf-3c" else "def2-SVP"
    cmd = [os.path.abspath(a.dump), xyz, "--method", ref["method"], "--basis", basis]
    if not a.no_fd:
        cmd.append("--fd")
    got = json.loads(subprocess.run(cmd, capture_output=True, text=True, timeout=3000, check=True).stdout)

    g, gr = got["gradient"], ref["gradient"]
    d_ref = max(abs(x - y) for x, y in zip(g, gr))
    d_e = abs(got["energy"] - ref["energy"])
    ok = d_ref <= TOL_REF and d_e <= 1e-8
    print(f"{a.case:28s} energy diff {d_e:.1e}  max|g - ref| = {d_ref:.2e}  (max|g| {max(map(abs, gr)):.3e})")
    if "fd_gradient" in got:
        d_fd = max(abs(x - y) for x, y in zip(g, got["fd_gradient"]))
        print(f"{a.case:28s} max|g - FD(own energy)| = {d_fd:.2e}")
        ok &= d_fd <= TOL_FD
    print("ok" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
