#!/usr/bin/env python3
"""Seminumerical exchange (-qm.exchange_method cosx, COSX / RIJCOSX) checks for the native QM engine.
Claude Generated (Nov 2026).

Modes (--check):
  reference   curcuma COSX total energy vs reference.json "pyscf_cosx": PySCF RHF/RKS whose K is the same COSX
              expression on curcuma's own grid (dumped with CURCUMA_QM_COSX_DUMP; AO values and the potential
              integrals int1e_grids from PySCF/libcint; scripts in docs/QM_COSX.md), --tol (default 1e-8 Eh). Also
              |E_cosx - E_exact| (curcuma, same run settings) <= --bound: the grid error must stay in its measured band.
  off         -qm.exchange_method exact gives the same 17-digit energy and gradient as no option (dump_qm_energy).
  threads     COSX energy with 1 and 4 threads bitwise equal (dump_qm_energy).
  fine        no overlap fit on a fine grid (75x302) approaches the exact exchange: |dE| <= --bound.
  rij         RIJCOSX consistency: E(-qm.ri j + cosx) - E(-qm.ri j) and E(cosx) - E(exact) (exchange error with fitted
              and with exact J) agree to --tol (the RI-J and COSX errors add up to second order).
  refuse      a request COSX does not support must fail loud (non-zero exit, message with --expect), never an energy:
              --extra holds the request (e.g. -gradient, -spin 1, -qm.ri jk).
Plain python3 (no numpy).
"""
import argparse, json, os, re, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--dump", default="", help="dump_qm_energy binary (off, threads)")
ap.add_argument("--ref", required=True)
ap.add_argument("--case", required=True)
ap.add_argument("--check", default="reference")
ap.add_argument("--tol", type=float, default=1e-8)
ap.add_argument("--bound", type=float, default=1e-3)
ap.add_argument("--extra", default="")
ap.add_argument("--expect", default="")
ap.add_argument("--xyz", default="", help="other structure (relative to the reference file) for --check refuse")
a = ap.parse_args()

here = os.path.dirname(os.path.abspath(a.ref))
ref = json.load(open(a.ref))["cases"][a.case]
xyz = os.path.normpath(os.path.join(here, a.xyz or ref["xyz"]))
base = ["-method", ref["method"], "-qm.basis", ref["basis"], "-qm.scf_threshold", "1e-10", "-qm.threads", "4"]
cosx = ["-qm.exchange_method", "cosx", "-qm.cosx_grid", ref["grid"]] + ref.get("flags", "").split()


def run(args):
    with tempfile.TemporaryDirectory() as tmp:
        out = subprocess.run([a.curcuma, "-sp", xyz, "-no_bmt"] + args + a.extra.split(), capture_output=True,
                             text=True, timeout=3000, cwd=tmp)
    return out.returncode, re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)


def energy(txt):
    ee = re.search(r"E_elec: (-?\d+\.\d+) Eh", txt)
    en = re.search(r"E_nn\s*: (-?\d+\.\d+) Eh", txt)
    if not ee or not en:
        print(txt[-3000:])
        sys.exit("no E_elec / E_nn in the output")
    return float(ee.group(1)) + float(en.group(1))


def dump(args, threads=4):
    cmd = [a.dump, xyz, ref["method"], "--qm", "basis", ref["basis"], "--qm", "scf_threshold", "1e-10",
           "--threads", str(threads), "--verbosity", "0"] + args
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000)
    lines = [l for l in out.stdout.splitlines() if l.startswith(("energy", "grad"))]
    if out.returncode != 0 or not lines:
        print(out.stdout[-2000:] + out.stderr[-2000:])
        sys.exit(f"dump_qm_energy exited with {out.returncode}")
    return lines


if a.check == "reference":
    rc, txt = run(base + cosx + ["-verbosity", "2"])
    if rc != 0:
        print(txt[-3000:])
        sys.exit(f"curcuma exited with {rc}")
    if "COSX(" not in txt:
        print(txt[-3000:])
        sys.exit("the COSX path did not run")
    e = energy(txt)
    rc, txt0 = run(base + ["-verbosity", "2"])
    e0 = energy(txt0)
    d = e - ref["pyscf_cosx"]
    g = e - e0
    ok = abs(d) <= a.tol and abs(g) <= a.bound
    print(f"{a.case}: E_cosx={e:.10f} PySCF-COSX={ref['pyscf_cosx']:.10f} diff={d:+.2e} (tol {a.tol:.0e}); "
          f"E_cosx-E_exact={g:+.3e} (bound {a.bound:.0e})" + ("  ok" if ok else "  FAIL"))
    sys.exit(0 if ok else 1)

if a.check == "off":
    plain = dump([])
    exact = dump(["--qm", "exchange_method", "exact"])
    ok = plain == exact
    print(f"{a.case}: {len(plain)} energy/gradient lines, exact == absent: {ok}")
    sys.exit(0 if ok else 1)

if a.check == "threads":
    cx = ["--qm", "exchange_method", "cosx", "--qm", "cosx_grid", ref["grid"]]
    e1 = dump(cx, 1)
    e4 = dump(cx, 4)
    ok = e1 == e4
    print(f"{a.case}: 1 thread {e1[0]}, 4 threads {e4[0]}, bitwise equal: {ok}")
    sys.exit(0 if ok else 1)

if a.check == "fine":
    rc, txt = run(base + ["-qm.exchange_method", "cosx", "-qm.cosx_grid", "75x302", "-qm.cosx_fit_overlap", "false",
                          "-verbosity", "2"])
    if rc != 0:
        print(txt[-3000:])
        sys.exit(f"curcuma exited with {rc}")
    e = energy(txt)
    e0 = energy(run(base + ["-verbosity", "2"])[1])
    ok = abs(e - e0) <= a.bound
    print(f"{a.case}: no fit, 75x302: E_cosx-E_exact={e - e0:+.3e} (bound {a.bound:.0e})" + ("  ok" if ok else "  FAIL"))
    sys.exit(0 if ok else 1)

if a.check == "rij":
    e = {}
    for key, args in (("exact", []), ("cosx", cosx), ("rij", ["-qm.ri", "j"]), ("rijcosx", ["-qm.ri", "j"] + cosx)):
        rc, txt = run(base + args + ["-verbosity", "2"])
        if rc != 0:
            print(txt[-3000:])
            sys.exit(f"curcuma exited with {rc} ({key})")
        e[key] = energy(txt)
    d1, d2 = e["cosx"] - e["exact"], e["rijcosx"] - e["rij"]
    ok = abs(d1 - d2) <= a.tol
    print(f"{a.case}: COSX error with exact J {d1:+.3e}, with RI-J {d2:+.3e}, difference {d1 - d2:+.2e} (tol {a.tol:.0e}); "
          f"RIJCOSX - exact {e['rijcosx'] - e['exact']:+.3e}" + ("  ok" if ok else "  FAIL"))
    sys.exit(0 if ok else 1)

if a.check == "refuse":
    rc, txt = run(base + cosx + ["-verbosity", "1"])
    ok = rc != 0 and all(x in txt for x in a.expect.split("|") if x)
    print(f"{a.case}: '{a.extra}' -> exit {rc}, message found: {all(x in txt for x in a.expect.split('|') if x)}"
          + ("  ok" if ok else "  FAIL"))
    if not ok:
        print(txt[-3000:])
    sys.exit(0 if ok else 1)

sys.exit("unknown --check " + a.check)
