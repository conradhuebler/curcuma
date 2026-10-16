#!/usr/bin/env python3
"""WP10 density fitting (RI-J / RI-JK) checks for the native HF engine. Claude Generated (Oct 2026).

Modes (--check):
  energy          curcuma -sp -method hf -qm.ri <ri> vs the PySCF total energy in reference.json
                  (--tol, default 1e-8 Eh; same method and auxiliary basis, so near rounding) and vs ORCA 6.1
                  (--tol-orca, default 1e-7 Eh; ORCA uses its own thresholds).
  gradient_error  -sp -gradient with -qm.ri must fail (non-zero exit, clear message), never return a number.
  hf3c            -method hf-3c ignores -qm.ri (warning) and returns the exact-integral energy.
  gpu             -gpu cuda with CURCUMA_QM_GPU_RI_CHECK=1: device J/K of the first build vs the CPU contraction
                  (<= 1e-10), energy vs PySCF, no fallback. Exit 77 (ctest SKIP) without plugin/device.
Plain python3 (no numpy): parses E_elec and E_nn (verbosity 2, 10 decimals).
"""
import argparse, json, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--mol", required=True)
ap.add_argument("--ri", default="jk")
ap.add_argument("--check", default="energy")
ap.add_argument("--set", default="pyscf", help="reference set: pyscf | pyscf_cartesian")
ap.add_argument("--tol", type=float, default=1e-8)
ap.add_argument("--tol-orca", type=float, default=1e-7)
ap.add_argument("--extra", default="", help="extra curcuma flags (one string, use --extra=...)")
a = ap.parse_args()

ref = json.load(open(a.ref))[a.mol]
xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), ref["xyz"])


def run(args, env=None):
    cmd = [a.curcuma, "-sp", xyz, "-no_bmt"] + args + a.extra.split()
    e = dict(os.environ)
    if env:
        e.update(env)
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000, env=e)
    return out.returncode, re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)


def energy(txt):
    ee = re.search(r"E_elec: (-?\d+\.\d+) Eh", txt)
    en = re.search(r"E_nn\s*: (-?\d+\.\d+) Eh", txt)
    if not ee or not en:
        print(txt[-3000:])
        sys.exit("no E_elec / E_nn in the output")
    return float(ee.group(1)) + float(en.group(1))


hf = ["-method", "hf", "-qm.scf_threshold", "1e-10", "-verbosity", "2"]

if a.check == "energy":
    rc, txt = run(hf + ["-qm.ri", a.ri])
    if rc != 0:
        print(txt[-3000:])
        sys.exit(f"curcuma exited with {rc}")
    if a.ri != "off" and f"RI-{a.ri.upper()} with" not in txt:
        print(txt[-3000:])
        sys.exit("the RI path did not run")
    e = energy(txt)
    d = e - ref[a.set][a.ri]
    ok = abs(d) <= a.tol
    line = f"{a.mol} ri={a.ri} set={a.set} E={e:.10f} ref={ref[a.set][a.ri]:.10f} diff={d:+.2e} (tol {a.tol:.0e})"
    if a.set == "pyscf":
        do = e - ref["orca"][a.ri]
        ok = ok and abs(do) <= a.tol_orca
        line += f"  ORCA diff={do:+.2e} (tol {a.tol_orca:.0e})"
    print(line + ("  ok" if ok else "  FAIL"))
    sys.exit(0 if ok else 1)

if a.check == "gradient_error":
    rc, txt = run(["-method", "hf", "-qm.ri", a.ri, "-gradient", "-verbosity", "1"])
    msg = "no analytic gradient with density fitting" in txt
    print(f"exit code {rc}, message {'found' if msg else 'missing'}")
    sys.exit(0 if (rc != 0 and msg) else 1)

if a.check == "hf3c":
    rc1, t1 = run(["-method", "hf-3c", "-qm.ri", a.ri, "-verbosity", "1"])
    rc0, t0 = run(["-method", "hf-3c", "-verbosity", "1"])
    m1 = re.search(r"Single Point Energy = (-?\d+\.\d+)", t1)
    m0 = re.search(r"Single Point Energy = (-?\d+\.\d+)", t0)
    warned = "hf-3c: -qm.ri" in t1 and "ignored" in t1
    ok = rc0 == 0 and rc1 == 0 and m0 and m1 and m0.group(1) == m1.group(1) and warned
    print(f"hf-3c with -qm.ri {a.ri}: {m1.group(1) if m1 else 'NA'}, without: {m0.group(1) if m0 else 'NA'}, "
          f"warning {'found' if warned else 'missing'}")
    sys.exit(0 if ok else 1)

if a.check == "gpu":
    rc, txt = run(hf + ["-qm.ri", a.ri, "-gpu", "cuda", "-qm.threads", "4"], env={"CURCUMA_QM_GPU_RI_CHECK": "1"})
    if "no GPU plugin" in txt or "not loadable" in txt or "cudaSetDevice failed" in txt:
        print("no usable CUDA device/plugin: skipped")
        sys.exit(77)
    m = re.search(r"QM GPU RI check: n=(\d+) naux=(\d+) max\|dJ\|=(\S+) max\|dK\|=(\S+)", txt)
    if rc != 0 or not m:
        print(txt[-3000:])
        sys.exit("device RI path did not run (no check line)")
    if "GPU fallbacks during this run" in txt:
        print(txt[-3000:])
        sys.exit("the run reported a GPU fallback")
    dJ, dK = float(m.group(3)), float(m.group(4))
    e = energy(txt)
    d = e - ref["pyscf"][a.ri]
    ok = dJ <= 1e-10 and dK <= 1e-10 and abs(d) <= a.tol
    print(f"{a.mol} ri={a.ri} device vs CPU: max|dJ|={dJ:.2e} max|dK|={dK:.2e}; E={e:.10f} diff PySCF={d:+.2e}  "
          f"{'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

sys.exit("unknown --check " + a.check)
