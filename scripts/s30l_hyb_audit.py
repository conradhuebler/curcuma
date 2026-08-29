#!/usr/bin/env python3
"""Audit curcuma vs xtb per-atom hybridization across all S30L systems.

For each (system, part), runs curcuma -verbosity 3 and xtb --gfnff --sp, parses the
per-atom hybridization (curcuma "Pi-fragment assignments: Atom i (Z=z, hyb=h)" and
xtb "atom neighbors ... sp-hybrid ... pi" table), and reports all mismatches
(Z, hyb_cur, hyb_xtb, pi_cur, pi_xtb) aggregated across the whole set.

Usage:
    python scripts/s30l_hyb_audit.py            # all 30, all parts
    python scripts/s30l_hyb_audit.py 11 15       # subset of systems

Claude Generated (Jul 2026, F3 element audit).
"""
import os
import re
import subprocess
import sys
import tempfile
import shutil
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN = os.path.join(ROOT, "test_cases", "s30l_test_set", "_run")
CURCUMA = os.path.join(ROOT, "release", "curcuma")
XTB = os.path.expanduser("~/Downloads/xtb-6.6.1/bin/xtb")

ZSYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 16: "S", 17: "Cl", 35: "Br", 53: "I"}


def curcuma_hyb(xyz, charge):
    cmd = [CURCUMA, "-sp", xyz, "-method", "gfnff", "-charge", str(charge),
           "-verbosity", "3", "-no_bmt"]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT).stdout
    out = re.sub(r"\x1b\[[0-9;]*m", "", out)
    hyb = {}  # atom(0-based) -> (Z, hyb, pi_fragment)
    for line in out.splitlines():
        m = re.match(r"\s*Atom\s+(\d+)\s+\(Z=(\d+),\s*hyb=(\d+)\):\s*fragment\s+(\d+)", line)
        if m:
            a, z, h, frag = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
            hyb[a] = (z, h, 1 if frag > 0 else 0)
    return hyb


def xtb_hyb(xyz):
    workdir = tempfile.mkdtemp(prefix="s30l_xtb_")
    xyz_dst = os.path.join(workdir, os.path.basename(xyz))
    shutil.copy(xyz, xyz_dst)
    out = subprocess.run([XTB, os.path.basename(xyz), "--gfnff", "--sp"],
                         capture_output=True, text=True, cwd=workdir).stdout
    shutil.rmtree(workdir, ignore_errors=True)
    hyb = {}  # atom(1-based) -> (sym, hyb, pi)
    f = False
    for line in out.splitlines():
        if line.startswith("  atom   neighbors"):
            f = True
            continue
        if f:
            if line.strip() == "" or line.startswith(" #") or line.startswith(" angle"):
                f = False
                continue
            toks = line.split()
            if len(toks) >= 8 and toks[0].isdigit():
                a = int(toks[0])
                sym = toks[1]
                h = int(toks[5])  # sp-hybrid column
                pi = int(toks[7])  # pi column
                hyb[a] = (sym, h, pi)
    return hyb


def main():
    args = [int(x) for x in sys.argv[1:]] or list(range(1, 31))
    mismatches = []  # (sys, part, atom0, Z, hyb_cur, hyb_xtb, pi_cur, pi_xtb)
    for s in args:
        for part in ("A", "B", "AB"):
            xyz = os.path.join(RUN, str(s), part + ".xyz")
            if not os.path.exists(xyz):
                continue
            pchrg = os.path.join(ROOT, "test_cases", "s30l_test_set", str(s), part, ".CHRG")
            charge = 0
            if os.path.exists(pchrg):
                with open(pchrg) as f:
                    try:
                        charge = int(f.read().strip())
                    except ValueError:
                        pass
            ch = curcuma_hyb(xyz, charge)
            xh = xtb_hyb(xyz)
            for a1, (sym, hx, px) in xh.items():
                a0 = a1 - 1
                if a0 in ch:
                    z, hc, pc = ch[a0]
                    if hc != hx or pc != px:
                        mismatches.append((s, part, a0, z, hc, hx, pc, px))
    print(f"Total hybridization mismatches: {len(mismatches)}")
    # aggregate by (Z, hyb_cur, hyb_xtb)
    agg = Counter()
    pi_agg = Counter()
    for s, part, a, z, hc, hx, pc, px in mismatches:
        agg[(z, hc, hx)] += 1
        if pc != px:
            pi_agg[(z, pc, px)] += 1
    print(f"\nBy element (Z, hyb_cur, hyb_xtb) -> count:")
    for (z, hc, hx), c in sorted(agg.items()):
        print(f"  Z={z:<3}({ZSYM.get(z,'?'):<2}) cur_hyb={hc} xtb_hyb={hx}  : {c}")
    print(f"\nPi-system membership mismatches (Z, pi_cur, pi_xtb) -> count:")
    for (z, pc, px), c in sorted(pi_agg.items()):
        print(f"  Z={z:<3}({ZSYM.get(z,'?'):<2}) cur_pi={pc} xtb_pi={px}  : {c}")
    print(f"\nFirst 20 mismatches (sys/part/atom Z: cur->xtb hyb, cur->xtb pi):")
    for s, part, a, z, hc, hx, pc, px in mismatches[:20]:
        print(f"  {s}/{part} atom{a} Z={z}({ZSYM.get(z,'?')}): hyb {hc}->{hx}, pi {pc}->{px}")


if __name__ == "__main__":
    main()