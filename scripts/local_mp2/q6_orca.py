#!/usr/bin/env python3
"""Q6 of the local-MP2 feasibility study: ORCA 6.1 DLPNO-MP2 as a qualitative external comparison.

For each system: ORCA RHF/def2-SVP, then (a) canonical RI-MP2 with def2-SVP/C, (b) DLPNO-MP2 NormalPNO,
(c) DLPNO-MP2 TightPNO, all frozen core (ORCA default) and the same RHF settings. Recovery = E_corr(DLPNO) /
E_corr(RI-MP2) within ORCA, so ORCA's aux basis (def2-SVP/C, not JKFIT) cancels. Serial runs.
Geometries are the ones of lmp2_common (written into the work directory, not into the repository).

Usage: ORCA=<path to orca> LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q6_orca.py [systems...]
Claude Generated (Oct 2026).
"""
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads

ORCA = os.environ.get("ORCA", "/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca")
JOBS = {"rimp2": "! RHF RI-MP2 def2-SVP def2-SVP/C TightSCF",
        "dlpno_normal": "! RHF DLPNO-MP2 def2-SVP def2-SVP/C NormalPNO TightSCF",
        "dlpno_tight": "! RHF DLPNO-MP2 def2-SVP def2-SVP/C TightPNO TightSCF"}


def run(name, job):
    d = os.path.join(lc.work_dir(), "orca", f"{name}_{job}")
    os.makedirs(d, exist_ok=True)
    inp = os.path.join(d, "job.inp")
    with open(inp, "w") as f:
        f.write(JOBS[job] + "\n%maxcore 8000\n* xyz 0 1\n")
        for s, x in lc.SYSTEMS[name][1]():
            f.write(f"{s} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}\n")
        f.write("*\n")
    outp = os.path.join(d, "job.out")
    if not (os.path.exists(outp) and "ORCA TERMINATED NORMALLY" in open(outp).read()):
        env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(ORCA))
        t0 = time.time()
        with open(outp, "w") as fo:
            subprocess.run([ORCA, "job.inp"], cwd=d, stdout=fo, stderr=subprocess.STDOUT, env=env, check=False)
        print(f"#  {name} {job}: {time.time() - t0:.0f}s", flush=True)
    txt = open(outp).read()
    if "ORCA TERMINATED NORMALLY" not in txt:
        sys.exit(f"ORCA failed: {outp}")
    ehf = float(re.findall(r"Total Energy\s+:\s+(-?\d+\.\d+) Eh", txt)[-1])
    m = re.findall(r"MP2 CORRELATION ENERGY\s*:?\s+(-?\d+\.\d+)", txt)
    res = {"e_hf": ehf, "e_corr": float(m[-1]), "out": os.path.relpath(outp, lc.work_dir())}
    for key, pat in (("tcutpno", r"TCutPNO\s+=\s+(\S+)"), ("tcutdo", r"TCutDO\s+\.\.\.\s+(\S+)"),
                     ("pno_mean", r"Average number of PNOs per pair\s+(\S+)"),
                     ("dE_pno", r"Semicanonical PNO correction to energy\s+\.\.\.\s+(\S+)")):
        f = re.findall(pat, txt)
        if f:
            res[key] = float(f[-1])
    return res


def main():
    names = sys.argv[1:] or ["alkane12", "water16"]
    out = []
    for name in names:
        r = {job: run(name, job) for job in JOBS}
        rec = {"system": name, **r,
               "recovery_normal_pct": 100 * r["dlpno_normal"]["e_corr"] / r["rimp2"]["e_corr"],
               "recovery_tight_pct": 100 * r["dlpno_tight"]["e_corr"] / r["rimp2"]["e_corr"],
               "err_normal": r["dlpno_normal"]["e_corr"] - r["rimp2"]["e_corr"],
               "err_tight": r["dlpno_tight"]["e_corr"] - r["rimp2"]["e_corr"]}
        out.append(rec)
        print(f"{name}: E_corr RI-MP2(def2-SVP/C) {r['rimp2']['e_corr']:.8f}  NormalPNO err {rec['err_normal']:.2e}"
              f" ({rec['recovery_normal_pct']:.3f}%)  TightPNO err {rec['err_tight']:.2e}"
              f" ({rec['recovery_tight_pct']:.3f}%)  PNOs/pair {r['dlpno_normal'].get('pno_mean')} / {r['dlpno_tight'].get('pno_mean')}",
              flush=True)
    lc.save_json("q6_orca.json", out)


if __name__ == "__main__":
    main()
