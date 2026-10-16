#!/usr/bin/env python3
"""Native canonical MP2 (WP11) against the PySCF reference in reference.json.

Claude Generated (Oct 2026). Runs `curcuma -sp <xyz> -method <mp2|scs-mp2|sos-mp2> -verbosity 2`
and compares E(HF), E_OS, E_SS, E_corr (and the printed SCS/SOS values and the total of the chosen
variant) with the PySCF values (scripts/mp2_reference.py: PySCF RHF conv_tol 1e-12 + mp.MP2,
def2-SVP spherical, frozen=None or the chemical frozen core).

  --fc          adds -qm.frozen_core true and compares with the frozen-core reference
  --gpu         adds -gpu cuda; the run must report the device path and no fallback; exit 77
                (ctest SKIP) when no CUDA plugin / device is usable
  --cpu-check   also runs the CPU path (same flags without -gpu) and requires |dE_corr| <= --tol-cpu
  --extra=...   further flags (one string)
  --ri-ref F    RI-MP2 mode (methods ri-mp2, scs-ri-mp2, sos-ri-mp2): compare with F (reference_ri.json,
                scripts/rimp2_reference.py: PySCF with the same auxiliary basis files); --aux selects the
                auxiliary set (adds -qm.ri_mp2_aux), --scf rijk adds -qm.ri jk (RI-JK SCF reference)
Passes when every compared quantity agrees within --tol (default 1e-8 Eh; E(HF) within --tol-hf).
"""
import argparse, json, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--mol", required=True)
ap.add_argument("--method", default="mp2")
ap.add_argument("--fc", action="store_true")
ap.add_argument("--gpu", action="store_true")
ap.add_argument("--cpu-check", action="store_true")
ap.add_argument("--extra", default="")
ap.add_argument("--tol", type=float, default=1e-8)
ap.add_argument("--tol-hf", type=float, default=1e-8)
ap.add_argument("--tol-cpu", type=float, default=1e-10)
ap.add_argument("--ri-ref", default="")
ap.add_argument("--aux", default="def2-universal-JKFIT")
ap.add_argument("--scf", default="exact", choices=["exact", "rijk"])
a = ap.parse_args()

if a.ri_ref:
    # Claude Generated (Oct 2026, RI-MP2)
    RC = json.load(open(a.ri_ref))["cases"][a.mol]
    R = {"xyz": RC["xyz"], "hf": RC[a.scf]["hf"]}
    ref = RC[a.scf][a.aux]["fc" if a.fc else "full"]
    a.extra += " -qm.ri_mp2_aux " + a.aux + (" -qm.ri jk" if a.scf == "rijk" else "")
else:
    R = json.load(open(a.ref))["cases"][a.mol]
    ref = R["fc" if a.fc else "full"]
xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", R["xyz"])
PAT = {
    "hf": r"E\(HF\)\s*=\s*(-?\d+\.\d+) Eh",
    "os": r"E\(MP2 OS\)\s*=\s*(-?\d+\.\d+) Eh",
    "ss": r"E\(MP2 SS\)\s*=\s*(-?\d+\.\d+) Eh",
    "corr": r"E\(MP2 corr\)\s*=\s*(-?\d+\.\d+) Eh",
    "scs": r"E\(SCS-MP2 corr\)\s*=\s*(-?\d+\.\d+) Eh",
    "sos": r"E\(SOS-MP2 corr\)\s*=\s*(-?\d+\.\d+) Eh",
    "total": r"E\(" + re.escape(a.method) + r"\)\s*=\s*(-?\d+\.\d+) Eh",
}


def run(gpu):
    cmd = [a.curcuma, "-sp", xyz, "-method", a.method, "-qm.basis", "def2-SVP", "-qm.scf_threshold", "1e-10",
           "-qm.threads", "4", "-verbosity", "2", "-no_bmt"]
    if a.fc:
        cmd += ["-qm.frozen_core", "true"]
    cmd += ["-gpu", "cuda" if gpu else "none"] + a.extra.split()
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    if gpu and ("no GPU plugin" in txt or "not loadable" in txt):
        print("no usable CUDA device/plugin: skipped")
        sys.exit(77)
    val = {}
    for k, p in PAT.items():
        m = re.search(p, txt)
        if not m:
            print(" ".join(cmd))
            print(txt[-3000:])
            sys.exit(f"no '{k}' value in the output")
        val[k] = float(m.group(1))
    return val, txt


v, txt = run(a.gpu)
for line in txt.splitlines():
    if "MP2 timing" in line or "canonical MP2" in line or "RI-MP2 (aux" in line:
        print(line.strip())
if a.ri_ref and ("RI-MP2 (aux " + a.aux) not in txt:
    sys.exit("the RI-MP2 path with aux " + a.aux + " did not run")
if a.gpu:
    dev = r"RI-MP2 timing:.*energy on cuda" if a.ri_ref else r"MP2 timing:.*\(cuda"
    if "MP2 timing" not in txt or re.search(dev, txt) is None:
        sys.exit("the device MP2 path did not run")
    if "GPU fallbacks during this run" in txt:
        print(txt[-3000:])
        sys.exit("the -gpu run reported a fallback")

expect_total = R["hf"] + {"mp2": ref["e_corr"], "scs-mp2": 1.2 * ref["e_os"] + ref["e_ss"] / 3.0,
                          "sos-mp2": 1.3 * ref["e_os"]}[a.method.replace("ri-", "")]
checks = [
    ("E(HF)", v["hf"], R["hf"], a.tol_hf),
    ("E_OS", v["os"], ref["e_os"], a.tol),
    ("E_SS", v["ss"], ref["e_ss"], a.tol),
    ("E_corr", v["corr"], ref["e_corr"], a.tol),
    ("E_SCS", v["scs"], 1.2 * ref["e_os"] + ref["e_ss"] / 3.0, a.tol),
    ("E_SOS", v["sos"], 1.3 * ref["e_os"], a.tol),
    ("E_total(" + a.method + ")", v["total"], expect_total, a.tol + a.tol_hf),
]
ok = True
for name, got, want, tol in checks:
    d = abs(got - want)
    good = d <= tol
    ok &= good
    print(f"{name:18s} curcuma={got:.12f} pyscf={want:.12f} diff={d:.2e} tol={tol:.0e} {'ok' if good else 'FAIL'}")
if a.cpu_check:
    vc, _ = run(False)
    d = abs(vc["corr"] - v["corr"])
    good = d <= a.tol_cpu
    ok &= good
    print(f"{'E_corr vs CPU':18s} gpu={v['corr']:.12f} cpu={vc['corr']:.12f} diff={d:.2e} tol={a.tol_cpu:.0e} "
          f"{'ok' if good else 'FAIL'}")
sys.exit(0 if ok else 1)
