#!/usr/bin/env python3
"""AO screening of the Kohn-Sham XC part (-qm.xc_ao_threshold) on the CPU.

Claude Generated (Oct 2026). Runs `curcuma -sp -dump_gradient` (CPU) with -qm.xc_ao_threshold 0 (no screening) and with
the default threshold, 4 threads, and once more with the default threshold on 1 thread. Passes when
  |E(default) - E(0)| <= --tol-e (1e-9 Eh) and max |g(default) - g(0)| <= --tol-g (1e-8 Eh/Bohr),
  the screened run on 1 and 4 threads gives the bitwise same gradient dump (not with --no-thread-check: the RI
  paths are not bitwise thread-independent, also without screening),
  the screened run reports "AO screening" with less than 100 % of the functions kept.
No numpy needed.
"""
import argparse, os, re, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="pbe")
ap.add_argument("--extra", default="")
ap.add_argument("--tol-e", type=float, default=1e-9)
ap.add_argument("--tol-g", type=float, default=1e-8)
ap.add_argument("--no-thread-check", action="store_true",
                help="skip 1 vs 4 threads (RI-J/RI-JK are not bitwise thread-independent, with or without screening)")
a = ap.parse_args()
work = tempfile.mkdtemp(prefix="qm_xc_screen_")


def run(tag, thr, threads):
    out = os.path.join(work, f"g_{tag}.txt")
    cmd = [os.path.abspath(a.curcuma), "-sp", os.path.abspath(a.xyz), "-method", a.method, "-gpu", "none",
           "-qm.threads", str(threads), "-qm.scf_threshold", "1e-10", "-no_bmt", "-verbosity", "2",
           "-dump_gradient", out] + (["-qm.xc_ao_threshold", thr] if thr is not None else []) + a.extra.split()
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=work, timeout=3000)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    if not os.path.exists(out):
        print(txt[-3000:])
        sys.exit(f"run {tag}: no gradient dump")
    raw = open(out).read()
    lines = raw.splitlines()
    e = float(re.search(r"energy (\S+) Eh", lines[1]).group(1))
    g = [float(x) * 0.529177210903 for l in lines[2:] for x in l.split()]  # Eh/Angstrom -> Eh/Bohr
    return e, g, raw, txt


e0, g0, _, _ = run("off", "0", 4)
e1, g1, raw4, txt = run("on4", None, 4)
raw1 = raw4 if a.no_thread_check else run("on1", None, 1)[2]
m = re.search(r"AO screening \(threshold (\S+)\): (\d+) blocks, ([\d.]+) % of the functions kept", txt)
de = abs(e1 - e0)
dg = max(abs(x - y) for x, y in zip(g0, g1))
same = raw1 == raw4
ok = m is not None and float(m.group(3)) < 100.0 and de <= a.tol_e and dg <= a.tol_g and same
print(f"{os.path.basename(a.xyz)} {a.method}: threshold {m.group(1) if m else '?'}, kept {m.group(3) if m else '?'} %; "
      f"|dE| {de:.2e} Eh (tol {a.tol_e:.0e}), max|dg| {dg:.2e} Eh/Bohr (tol {a.tol_g:.0e}); 1 vs 4 threads bitwise: "
      f"{'skipped' if a.no_thread_check else same}  {'ok' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
