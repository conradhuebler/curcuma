#!/usr/bin/env python3
"""ctest checker of the ECP matrix and ECP gradient on the device (ctest -L qm_gpu_ecp). Standard library only.

One case of qm_ecp/reference.json (optionally with --qm KEY=VALUE overrides, --charge/--spin), run by dump_qm_energy:
  1. CPU path and device path (--gpu cuda, CURCUMA_QM_GPU_ECP_CHECK=1) of the same binary: energy |dE| <= --etol
     (default 1e-10 Eh), gradient max |dg| <= --gtol (default 1e-10 Eh/Bohr);
  2. the device run's own comparison lines (same density, same geometry): ECP matrix max |V_dev - V_cpu| <= --vtol
     (default floor 1e-13 Eh, else 5e-15 of max|V|), radial orders identical
     (orders_equal=1), ECP gradient max |g_dev - g_cpu| <= --egtol (default 1e-11 Eh/Bohr);
  3. no counted qm_ecp fallback (the device ECP must run) and no other fallback except those in --allow.
Exit 77 (skip) without libcurcuma_cuda.so (CURCUMA_PLUGIN_DIR or the build tree) or without a CUDA device.

Claude Generated (Nov 2026).
"""
import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.normpath(os.path.join(HERE, ".."))


def gpu_available(dump):
    cands = [os.environ.get("CURCUMA_PLUGIN_DIR", "")]
    d = os.path.dirname(os.path.abspath(dump))
    for _ in range(4):
        cands.append(d)
        d = os.path.dirname(d)
    pdir = next((c for c in cands if c and os.path.exists(os.path.join(c, "libcurcuma_cuda.so"))), "")
    ok = bool(pdir)
    if ok:
        try:
            ok = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True).returncode == 0
        except OSError:
            ok = False
    if not ok:
        print("skip: no CUDA plugin found (%s) or no CUDA device" % cands)
        sys.exit(77)
    os.environ["CURCUMA_PLUGIN_DIR"] = pdir


def run(a, xyz, method, qm, gpu):
    cmd = [a.dump, xyz, method, "--threads", str(a.threads), "--gradient"]
    for k, v in qm.items():
        cmd += ["--qm", k, json.dumps(v) if not isinstance(v, str) else v]
    if a.charge or a.spin:
        cmd += ["--charge", str(a.charge), "--spin", str(a.spin)]
    env = dict(os.environ)
    env["CURCUMA_QM_BASIS_EMBEDDED_ONLY"] = "1"
    if gpu:
        cmd += ["--gpu", "cuda"]
        env["CURCUMA_QM_GPU_ECP_CHECK"] = "1"
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    out = p.stdout + p.stderr
    m = re.search(r"^energy (\S+)$", out, re.M)
    if p.returncode != 0 or not m:
        print(out[-4000:])
        raise SystemExit("dump_qm_energy failed: " + " ".join(cmd))
    g = [[float(x) for x in l.split()[2:5]] for l in out.splitlines() if l.startswith("grad ")]
    fb = {f.group(1): int(f.group(2)) for f in re.finditer(r"^fallback (\S+) (\d+)$", out, re.M)}
    return float(m.group(1)), g, fb, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--dump", required=True)
    ap.add_argument("--reference", default=os.path.join(HERE, "reference.json"))
    ap.add_argument("--qm", nargs="*", default=[], help="engine options KEY=VALUE on top of the case")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--etol", type=float, default=1e-10)
    ap.add_argument("--gtol", type=float, default=1e-10)
    ap.add_argument("--vtol", type=float, default=1e-13)
    ap.add_argument("--egtol", type=float, default=1e-11)
    ap.add_argument("--allow", nargs="*", default=[], help="counted fallbacks that may occur (never qm_ecp)")
    a = ap.parse_args()
    gpu_available(a.dump)
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = dict(ref["qm"])
    qm.setdefault("scf_threshold", 1e-10)
    for kv in a.qm:
        k, v = kv.split("=", 1)
        try:
            qm[k] = json.loads(v)
        except ValueError:
            qm[k] = v
    e0, g0, _, _ = run(a, xyz, ref["method"], qm, False)
    e1, g1, fb, out = run(a, xyz, ref["method"], qm, True)
    ok = True
    dg = max((abs(x - y) for r1, r2 in zip(g0, g1) for x, y in zip(r1, r2)), default=0.0)
    print("%-28s %s CPU %.12f GPU %.12f dE %+.2e  max|dg| %.2e Eh/Bohr" % (a.case, " ".join(a.qm), e0, e1, e1 - e0, dg))
    if abs(e1 - e0) > a.etol or dg > a.gtol or len(g0) != len(g1) or not g0:
        ok = False
    mats = re.findall(r"QM GPU ECP check: centres=(\d+) max\|dev-cpu\|=(\S+) max\|cpu\|=(\S+) orders_equal=(\d) "
                      r"cpu_ordered=(\d+) dev_ms=(\S+) cpu_ms=(\S+)", out)
    grads = re.findall(r"QM GPU ECP gradient check: max\|dev-cpu\|=(\S+) max\|cpu\|=(\S+) dev_ms=(\S+) cpu_ms=(\S+)", out)
    if not mats or not grads:
        print(out[-3000:])
        print("FAILED: no device ECP check line (device ECP did not run)")
        ok = False
    for c, dv, vmax, eq, cpuord, dms, cms in mats:
        # absolute floor --vtol, else 5e-15 of max|V| (a few ulp of a sum of several hundred terms; measured XeF2 1.7e-15,
        # CeF4 2.2e-15 relative with the device `pow`, Nov 4 2026); energies and the ECP gradient have their own limits
        lim = max(a.vtol, 5e-15 * float(vmax))
        print("  ECP matrix: centres %s max|V_dev - V_cpu| %s (limit %.1e, max|V| %s) orders equal %s, CPU-ordered "
              "centres %s, device %s ms, CPU %s ms" % (c, dv, lim, vmax, eq, cpuord, dms, cms))
        if float(dv) > lim or eq != "1":
            ok = False
    for dv, gmax, dms, cms in grads:
        print("  ECP gradient: max|g_dev - g_cpu| %s Eh/Bohr (limit %.1e, max|g| %s) device %s ms, CPU %s ms" % (
            dv, a.egtol, gmax, dms, cms))
        if float(dv) > a.egtol:
            ok = False
    bad = {k: v for k, v in fb.items() if k == "qm_ecp" or k not in a.allow}
    print("  counted GPU fallbacks: %s (allowed: %s)" % (fb or "none", a.allow or "none"))
    if bad:
        print("FAILED: device path fell back: %s" % bad)
        ok = False
    if not ok:
        raise SystemExit("FAILED: %s" % a.case)


if __name__ == "__main__":
    main()
