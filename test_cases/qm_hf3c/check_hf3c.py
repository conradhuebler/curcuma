#!/usr/bin/env python3
"""Native HF-3c gate: curcuma `-method hf-3c` vs independent references.

Runs `curcuma -sp <xyz> -method hf-3c` and compares each term of
E = HF/MINIX + D3(BJ) + gCP + SRB with reference.json (PySCF + simple-dftd3,
see scripts/hf3c_reference.py). For H2O the total is also checked against the
ORCA 6.1 `! HF-3c` value stored in the same file.

Usage: check_hf3c.py --curcuma <bin> --ref reference.json --mol <name>
Claude Generated (Sep 2026).
"""
import argparse, json, os, re, subprocess, sys, tempfile

TOL = {"hf": 1e-8, "d3": 1e-9, "gcp": 1e-9, "srb": 1e-9, "total": 1e-8}
PATTERN = {"hf": r"E\(HF/MINIX\)", "d3": r"E\(D3BJ\)", "gcp": r"E\(gCP\)",
           "srb": r"E\(SRB\)", "total": r"E\(HF-3c\)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--mol", required=True)
    a = ap.parse_args()

    a.curcuma = os.path.abspath(a.curcuma)
    ref = json.load(open(a.ref))
    r = ref["molecules"][a.mol]
    xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", r["xyz"])
    with tempfile.TemporaryDirectory() as tmp:
        out = subprocess.run([a.curcuma, "-sp", xyz, "-method", "hf-3c", "-verbosity", "1",
                              "-no_bmt", "-qm.scf_threshold", "1e-8"],
                             cwd=tmp, capture_output=True, text=True, timeout=600).stdout
    got = {}
    for key, pat in PATTERN.items():
        m = re.search(pat + r"\s*=\s*(-?\d+\.\d+)", out)
        if not m:
            print(out); sys.exit(f"{a.mol}: no '{key}' line in the curcuma output")
        got[key] = float(m.group(1))

    fail = False
    for key in TOL:
        d = got[key] - r[key]
        ok = abs(d) <= TOL[key]
        fail |= not ok
        print(f"{a.mol:13s} {key:5s} curcuma {got[key]:18.12f}  ref {r[key]:18.12f}  "
              f"diff {d:+.2e}  {'ok' if ok else 'FAIL (tol %.0e)' % TOL[key]}")
    if a.mol == "H2O":
        d = got["total"] - ref["orca_h2o"]["total"]
        ok = abs(d) <= 1e-8
        fail |= not ok
        print(f"{a.mol:13s} total vs ORCA 6.1 {ref['orca_h2o']['total']:.12f}  diff {d:+.2e}  "
              f"{'ok' if ok else 'FAIL'}")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
