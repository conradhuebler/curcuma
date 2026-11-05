#!/usr/bin/env python3
"""ctest checker for the native QM SCF level shift and stability analysis (ctest -L qm_stability). Standard library only.

Checks (--check):
  ref     -qm.scf_stability check on one case of reference.json: energy |dE| <= 1e-8 Eh and every lowest eigenvalue of
          (A + B) against PySCF's dense matrix (scripts/qm_stability_reference.py), |d| <= --tol (default 1e-6 Eh)
  follow  -qm.scf_stability follow: energy and <S^2> of the solution reached by following against PySCF's followed
          solution (|dE| <= 1e-7 Eh, |dS^2| <= 1e-4), and the final analysis reports stable
  shift   -qm.scf_level_shift --shift: the converged energy (and with --gradient the analytic gradient) equals the run
          without shift (|dE| <= 1e-9 Eh, max |dg| <= 1e-8 Eh/Bohr); --extra adds -qm settings to both runs
          (e.g. scf_mode=plain); --expect-unconverged-without: the run without shift must fail (a case that needs it)
  --threads N (ref, follow): threads of dump_qm_energy (the KS cases use 4: the 1-thread V_xc is slow on this grid)

Claude Generated (Nov 2026).
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.normpath(os.path.join(HERE, ".."))
GRID = {"grid": "sg1", "vwn_variant": "vwn5"}


def run(dump, case, qm, gradient=False, expect_fail=False, threads=1):
    cmd = [dump, os.path.join(TC, case["xyz"]), case["method"], "--charge", str(case["charge"]), "--spin",
           str(case["spin"]), "--threads", str(threads)]
    q = {"basis": case["basis"], "scf_guess": "atom"}
    if case["reference"] == "uhf":
        q["reference"] = "uhf"
    if case["method"] != "hf":
        q.update(GRID)
    q.update(qm)
    for k, v in q.items():
        cmd += ["--qm", k, v if isinstance(v, str) else json.dumps(v)]
    if gradient:
        cmd.append("--gradient")
    env = dict(os.environ)
    env["CURCUMA_QM_BASIS_EMBEDDED_ONLY"] = "1"
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    out = p.stdout + p.stderr
    if expect_fail:
        return p.returncode, out
    if p.returncode != 0:
        print(out[-3000:])
        print("FAIL: dump_qm_energy exit", p.returncode)
        sys.exit(1)
    r = {"grad": {}, "stability": {}}
    for ln in out.splitlines():
        t = ln.split()
        if len(t) == 2 and t[0] in ("energy", "s2"):
            r[t[0]] = float(t[1])
        elif len(t) == 5 and t[0] == "grad":
            r["grad"][int(t[1])] = [float(v) for v in t[2:]]
        elif len(t) == 3 and t[0] == "stability":
            r["stability"][t[1]] = float(t[2])
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--check", required=True, choices=["ref", "follow", "shift"])
    ap.add_argument("--case", required=True)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--shift", type=float, default=0.5)
    ap.add_argument("--gradient", action="store_true")
    ap.add_argument("--extra", nargs="*", default=[])
    ap.add_argument("--expect-unconverged-without", action="store_true")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    with open(os.path.join(HERE, "reference.json")) as f:
        ref = json.load(f)["cases"][a.case]
    extra = {}
    for kv in a.extra:
        k, v = kv.split("=", 1)
        try:
            extra[k] = json.loads(v)
        except ValueError:
            extra[k] = v
    ok = True
    if a.check == "ref":
        r = run(a.dump, ref, dict(extra, scf_stability="check"), threads=a.threads)
        dE = r["energy"] - ref["energy"]
        print("energy %.12f  pyscf %.12f  dE %.2e" % (r["energy"], ref["energy"], dE))
        ok &= abs(dE) <= 1e-8
        for key in ("rhf_rhf", "rhf_uhf", "uhf_uhf"):
            if key not in ref:
                continue
            if key not in r["stability"]:
                print("FAIL: no stability value", key)
                ok = False
                continue
            d = r["stability"][key] - ref[key]
            print("%s %.10f  pyscf dense %.10f  d %.2e (tol %.0e)" % (key, r["stability"][key], ref[key], d, a.tol))
            ok &= abs(d) <= a.tol
    elif a.check == "follow":
        r = run(a.dump, ref, dict(extra, scf_stability="follow"), threads=a.threads)
        dE = r["energy"] - ref["followed_energy"]
        print("followed energy %.12f  pyscf %.12f  dE %.2e" % (r["energy"], ref["followed_energy"], dE))
        ok &= abs(dE) <= 1e-7
        s2 = r.get("s2", 0.0)
        print("<S^2> %.8f  pyscf %.8f" % (s2, ref["followed_s2"]))
        ok &= abs(s2 - ref["followed_s2"]) <= 1e-4
        last = [v for k, v in r["stability"].items() if k in ("rhf_rhf", "rhf_uhf", "uhf_uhf")]
        print("stability values after following:", r["stability"])
        ok &= r["stability"].get("stable", 0.0) == 1.0 and r["stability"].get("followed", 0.0) >= 1.0
    else:
        q0 = dict(extra)
        q1 = dict(extra, scf_level_shift=a.shift)
        if a.expect_unconverged_without:
            rc, out = run(a.dump, ref, q0, gradient=a.gradient, expect_fail=True)
            print("without shift: exit %d (expected != 0)" % rc)
            ok &= rc != 0
            r1 = run(a.dump, ref, q1, gradient=a.gradient)
            dE = r1["energy"] - ref["energy"]
            print("with shift %.2f: energy %.12f  pyscf %.12f  dE %.2e" % (a.shift, r1["energy"], ref["energy"], dE))
            ok &= abs(dE) <= 1e-8
        else:
            r0 = run(a.dump, ref, q0, gradient=a.gradient)
            r1 = run(a.dump, ref, q1, gradient=a.gradient)
            dE = r1["energy"] - r0["energy"]
            print("energy without %.12f  with shift %.2f %.12f  dE %.2e" % (r0["energy"], a.shift, r1["energy"], dE))
            ok &= abs(dE) <= 1e-9
            if a.gradient:
                g = max(abs(x - y) for i in r0["grad"] for x, y in zip(r0["grad"][i], r1["grad"][i]))
                print("max |dg| %.2e Eh/Bohr" % g)
                ok &= g <= 1e-8
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
