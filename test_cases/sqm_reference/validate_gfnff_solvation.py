#!/usr/bin/env python3
"""
GFN-FF self-consistent ALPB solvation validation (Claude Generated, WP5, June 2026).

Runs the real `curcuma` binary with native GFN-FF + ALPB implicit solvation and
compares the total energy against committed xtb 6.7.1 references
(test_cases/sqm_reference/reference_data/gfnff_alpb_xtb.ref.json, produced by
`xtb mol.xyz --gfnff --alpb SOLVENT`). Needs NO xtb at test time.

The self-consistent EEQ<->ALPB coupling (A_eeq += Born matrix, matching the gfnff
reference gfnff_engrad.F90:1346-1350) makes the native total match xtb to <=1e-8 Eh.

Usage:
  validate_gfnff_solvation.py --curcuma <path> --ref <ref.json> --xyz-dir <dir>
                              --mol <name> --solvent <name> [--tol <Eh>] [--quiet]
Exit 0 on pass, 1 on fail.
"""
import argparse
import json
import re
import subprocess
import sys

ANSI = re.compile(r"\x1b\[[0-9;]*m")


def run_curcuma(curcuma, xyz, solvent, quiet):
    cmd = [curcuma, "-sp", xyz, "-method", "gfnff",
           "-gfnff.solvent", solvent, "-gfnff.solvent_model", "alpb",
           "-verbosity", "1"]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        if not quiet:
            print(f"  TIMEOUT: {' '.join(cmd)}")
        return None
    text = ANSI.sub("", out.stdout + out.stderr)
    m = re.findall(r"Single Point Energy = (-?[0-9.]+)", text)
    if not m:
        if not quiet:
            print(f"  no energy parsed (exit {out.returncode})")
        return None
    return float(m[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--xyz-dir", required=True)
    ap.add_argument("--mol", required=True)
    ap.add_argument("--solvent", required=True)
    ap.add_argument("--tol", type=float, default=1e-7)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    ref = json.load(open(args.ref))
    key = f"{args.mol}_{args.solvent}"
    if key not in ref["data"]:
        print(f"FAIL {key}: not in reference {args.ref}")
        return 1
    e_ref = ref["data"][key]

    e = run_curcuma(args.curcuma, f"{args.xyz_dir}/{args.mol}.xyz", args.solvent, args.quiet)
    if e is None:
        print(f"FAIL {key}: curcuma produced no energy")
        return 1

    dE = e - e_ref
    ok = abs(dE) <= args.tol
    if not args.quiet or not ok:
        print(f"{'PASS' if ok else 'FAIL'} {key:24s} gfnff+alpb: "
              f"E_curcuma={e:.8f}  E_xtb={e_ref:.8f}  dE={dE*1000:+.5f} mEh  "
              f"tol={args.tol*1000:.4f} mEh")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
