#!/usr/bin/env python3
"""ORCA 6.1 DLPNO-CCSD against ORCA canonical CCSD (qualitative comparison for the local-CCSD prototype).

For each system (geometries of ccsd_canonical.geometry, written into the work directory only): ORCA RHF/def2-SVP,
frozen core (ORCA default), (a) canonical CCSD with exact integrals, (b) DLPNO-CCSD NormalPNO, (c) DLPNO-CCSD
TightPNO, both with def2-SVP/C, (e) DLPNO-CCSD with TCutPNO = TCutPairs = 1e-12, TCutDO = TCutMKN = 0
("complete PNO", same aux basis and code: the reference that isolates the PNO/pair truncation), (d) canonical RI-MP2 with def2-SVP/C and the MP2 energy of (a), whose difference
estimates the part of the DLPNO error that comes from the auxiliary basis (DLPNO uses RI, (a) does not).
Error = E_corr(DLPNO) - E_corr(canonical CCSD), positive = correlation missing. Serial runs, %maxcore 8000.

Usage: LCCSD_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_ccsd/orca_dlpno.py [--complete water6] [systems...]
       (--complete: systems that also get the complete-PNO job; it took 19 min for (H2O)6 but ~65 min per CC iteration
       and > 180 GB scratch for naphthalene, so it is opt-in)
       env ORCA=<orca binary>, LCCSD_ORCA_DIR=<dir> (ORCA run directories; default the work directory. The
       complete-PNO job writes > 100 GB of scratch for C8H18, so put it on a disk, not on a tmpfs)
Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccsd_canonical as cc0  # first: sets the thread environment before numpy loads
import lmp2_common as lc

import json
import re
import subprocess
import time

JOBS = {"ccsd": "! RHF CCSD def2-SVP VeryTightSCF\n%mdci\n  STol 1e-8\nend\n",
        "dlpno_normal": "! RHF DLPNO-CCSD def2-SVP def2-SVP/C NormalPNO VeryTightSCF\n",
        "dlpno_tight": "! RHF DLPNO-CCSD def2-SVP def2-SVP/C TightPNO VeryTightSCF\n",
        "rimp2": "! RHF RI-MP2 def2-SVP def2-SVP/C VeryTightSCF\n",
        "dlpno_complete": "! RHF DLPNO-CCSD def2-SVP def2-SVP/C VeryTightSCF\n%mdci\n  TCutPNO 1e-12\n  TCutPairs 1e-12\n"
                          "  TCutDO 0\n  TCutMKN 0\nend\n"}


def run(name, job):
    d = os.path.join(os.environ.get("LCCSD_ORCA_DIR", lc.work_dir()), "orca_dlpno", f"{name}_{job}")
    os.makedirs(d, exist_ok=True)
    outp = os.path.join(d, "job.out")
    if not (os.path.exists(outp) and "ORCA TERMINATED NORMALLY" in open(outp).read()):
        with open(os.path.join(d, "job.inp"), "w") as f:
            f.write(JOBS[job] + "%maxcore 8000\n* xyz 0 1\n")
            for s, x in cc0.geometry(name):
                f.write(f"{s} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}\n")
            f.write("*\n")
        env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(cc0.ORCA))
        t0 = time.time()
        with open(outp, "w") as fo:
            subprocess.run([cc0.ORCA, "job.inp"], cwd=d, stdout=fo, stderr=subprocess.STDOUT, env=env, check=False)
        print(f"#  {name} {job}: {time.time() - t0:.0f}s", flush=True)
    txt = open(outp).read()
    if "ORCA TERMINATED NORMALLY" not in txt:
        sys.exit(f"ORCA failed: {outp}")
    res = {"e_hf": float(re.findall(r"Total Energy\s+:\s+(-?\d+\.\d+) Eh", txt)[-1]), "out": outp}
    if job == "rimp2":
        res["e_corr"] = float(re.findall(r"MP2 CORRELATION ENERGY\s*:?\s+(-?\d+\.\d+)", txt)[-1])
    else:
        res["e_hf"], res["e_corr"] = cc0.orca_property_energies(os.path.join(d, "job.property.txt"))
        m = re.findall(r"^\s*0\s+-\d+\.\d+\s+(-\d+\.\d+)", txt.split("COUPLED CLUSTER ITERATIONS")[-1], re.M)
        if job == "ccsd" and m:
            res["e_mp2"] = float(m[0])  # iteration 0 of the canonical CCSD = MP2 (9 dp)
    for key, pat in (("pno_mean", r"Average number of PNOs per pair\s+\.*\s*(\S+)"),
                     ("dE_pno", r"PNO truncation correction\s+\.\.\.\s+(\S+)")):
        f = re.findall(pat, txt)
        if f:
            try:
                res[key] = float(f[-1])
            except ValueError:
                pass
    return res


def main():
    args = sys.argv[1:]
    complete = ["water6"]
    if "--complete" in args:
        k = args.index("--complete")
        complete = [x for x in args[k + 1].split(",") if x]
        args = args[:k] + args[k + 2:]
    names = args or ["alkane8", "water6", "acene2"]
    out = []
    for name in names:
        r = {job: run(name, job) for job in JOBS if job != "dlpno_complete" or name in complete}
        natm = len(cc0.geometry(name))
        ref = r["ccsd"]["e_corr"]
        rec = {"system": name, "natm": natm, **r,
               "err_normal": r["dlpno_normal"]["e_corr"] - ref, "err_tight": r["dlpno_tight"]["e_corr"] - ref,
               "aux_err_mp2": r["rimp2"]["e_corr"] - r["ccsd"].get("e_mp2", float("nan")),
               "complete_minus_ccsd": None, "err_normal_vs_complete": None, "err_tight_vs_complete": None}
        if "dlpno_complete" in r:
            c = r["dlpno_complete"]["e_corr"]
            rec.update({"complete_minus_ccsd": c - ref, "err_normal_vs_complete": r["dlpno_normal"]["e_corr"] - c,
                        "err_tight_vs_complete": r["dlpno_tight"]["e_corr"] - c})
        out.append(rec)
        print(f"{name}: E_corr CCSD {ref:.8f}  NormalPNO err {rec['err_normal']:.2e} ({rec['err_normal'] / natm:.1e}/atom)"
              f"  TightPNO err {rec['err_tight']:.2e} ({rec['err_tight'] / natm:.1e}/atom)"
              f"  RI-MP2(SVP/C) - MP2 {rec['aux_err_mp2']:.2e}  complete-PNO DLPNO - CCSD {rec['complete_minus_ccsd']}"
              f"  vs complete: Normal {rec['err_normal_vs_complete']} Tight {rec['err_tight_vs_complete']}", flush=True)
        lc.save_json("orca_dlpno.json", out)


if __name__ == "__main__":
    main()
