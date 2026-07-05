#!/usr/bin/env python3
"""Per-bond GFN-FF parameter comparison: curcuma vs xtb 6.6.1.

Dumps curcuma per-bond (pibo, fc, r0, alpha, energy) at -verbosity 3 and xtb
per-bond (piBO, kbond, R0, alp, R) from --verbose, matches by atom pair, prints
deltas sorted by |energy delta|.

pibo recovered from fpi: fpi = 1 - 1.0*(0.315 - pibo) = 0.685 + pibo  ->  pibo = fpi - 0.685

Usage:
    python scripts/s30l_bond_compare.py 11 A
    python scripts/s30l_bond_compare.py 15 A

Claude Generated (Jul 2026, F3 diagnosis).
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN = os.path.join(ROOT, "test_cases", "s30l_test_set", "_run")
CURCUMA = os.path.join(ROOT, "release", "curcuma")
XTB = os.path.expanduser("~/Downloads/xtb-6.6.1/bin/xtb")


def curcuma_bonds(xyz, charge):
    cmd = [CURCUMA, "-sp", xyz, "-method", "gfnff", "-charge", str(charge),
           "-verbosity", "3", "-no_bmt"]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT).stdout
    out = re.sub(r"\x1b\[[0-9;]*m", "", out)
    bonds = {}  # (i,j) -> dict
    # BOND_FACTORS: i, j, bond_i=..., bond_j=..., bstr=..., fqq=..., ringf=..., fheavy=..., fpi=..., fxh=..., fcn=..., fc=..., qa1=..., qa2=...
    for line in out.splitlines():
        m = re.match(r"\s*BOND_FACTORS:\s*(\d+),\s*(\d+),\s*bond_i=([\d.]+),\s*bond_j=([\d.]+),\s*bstr=([\d.]+),\s*fqq=([\d.]+),\s*ringf=([\d.]+),\s*fheavy=([\d.]+),\s*fpi=([\d.]+),\s*fxh=([\d.]+),\s*fcn=([\d.]+),\s*fc=([-\d.]+)", line)
        if m:
            i, j = int(m.group(1)), int(m.group(2))
            key = (min(i, j), max(i, j))
            fpi = float(m.group(9))
            bonds.setdefault(key, {})
            bonds[key]["fpi"] = fpi
            bonds[key]["pibo"] = fpi - 0.685
            bonds[key]["fc"] = float(m.group(12))
            bonds[key]["fqq"] = float(m.group(6))
            bonds[key]["ringf"] = float(m.group(7))
            bonds[key]["bstr"] = float(m.group(5))
    # BOND_CSV: idx, i, j, Zi, Zj, rij, r0, fc, alpha, fqq, energy
    for line in out.splitlines():
        m = re.match(r"\s*BOND_CSV:\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*([\d.]+),\s*([\d.]+),\s*([-\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([-\d.]+)", line)
        if m:
            i, j = int(m.group(2)), int(m.group(3))
            key = (min(i, j), max(i, j))
            bonds.setdefault(key, {})
            bonds[key]["rij"] = float(m.group(6))
            bonds[key]["r0"] = float(m.group(7))
            bonds[key]["alpha"] = float(m.group(9))
            bonds[key]["energy"] = float(m.group(11))
    return bonds


def xtb_bonds(xyz):
    import tempfile, shutil
    workdir = tempfile.mkdtemp(prefix="s30l_xtb_")
    xyz_dst = os.path.join(workdir, os.path.basename(xyz))
    shutil.copy(xyz, xyz_dst)
    cmd = [XTB, os.path.basename(xyz), "--gfnff", "--sp", "--verbose"]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir).stdout
    shutil.rmtree(workdir, ignore_errors=True)
    bonds = {}
    in_table = False
    for line in out.splitlines():
        if line.startswith(" bond atoms"):
            in_table = True
            continue
        if in_table:
            if line.strip() == "" or line.startswith(" angle") or line.startswith(" #"):
                in_table = False
                continue
            # "C  C      3    1      2    5     1.419   1.286   0.500   1.023  -0.228   0.491"
            toks = line.split()
            if len(toks) < 11:
                continue
            try:
                # sym_i sym_j atom_i atom_j type ring R R0 piBO fqq kbond alp
                ai = int(toks[2])
                aj = int(toks[3])
                R = float(toks[6])
                R0 = float(toks[7])
                pibo = float(toks[8])
                fqq = float(toks[9])
                kbond = float(toks[10])
                alp = float(toks[11]) if len(toks) > 11 else float("nan")
            except (ValueError, IndexError):
                continue
            # xtb 1-based -> 0-based
            key = (min(ai - 1, aj - 1), max(ai - 1, aj - 1))
            bonds[key] = {"R": R, "R0": R0, "pibo": pibo, "fqq": fqq,
                          "kbond": kbond, "alp": alp}
    return bonds


def compare(sysid, part):
    import math
    BOHR2ANG = 1.0 / 1.8897259886  # Bohr -> Angstrom
    base = os.path.join(RUN, str(sysid), part)
    xyz = base + ".xyz"
    pchrg = os.path.join(ROOT, "test_cases", "s30l_test_set", str(sysid), part, ".CHRG")
    charge = 0
    if os.path.exists(pchrg):
        with open(pchrg) as f:
            try:
                charge = int(f.read().strip())
            except ValueError:
                charge = 0
    cb = curcuma_bonds(xyz, charge)
    xb = xtb_bonds(xyz)
    print(f"\n=== system {sysid}/{part} (charge={charge}) ===")
    print(f"  curcuma bonds: {len(cb)}   xtb bonds: {len(xb)}")
    rows = []
    tot_e_cur = 0.0
    tot_e_xtb = 0.0
    for key in sorted(set(cb) | set(xb)):
        c = cb.get(key, {})
        x = xb.get(key, {})
        # curcuma energy is printed directly (Eh). Recompute in angstrom to cross-check.
        e_cur = c.get("energy", 0.0)
        # xtb bond energy = kbond * exp(-alp*(R-R0)^2); R, R0 in Angstrom, alp in 1/angstrom^2
        R = x.get("R", 0.0)
        R0x = x.get("R0", 0.0)
        alpx = x.get("alp", 0.0)
        kbx = x.get("kbond", 0.0)
        e_xtb = kbx * math.exp(-alpx * (R - R0x) ** 2) if xb.get(key) else 0.0
        tot_e_cur += e_cur
        tot_e_xtb += e_xtb
        de = e_cur - e_xtb
        dpibo = c.get("pibo", 0.0) - x.get("pibo", 0.0)
        dfc = c.get("fc", 0.0) - x.get("kbond", 0.0)
        # convert curcuma r0 (Bohr) and alpha (1/Bohr^2) to Angstrom for comparison
        cur_r0_ang = c.get("r0", 0.0) * BOHR2ANG if c.get("r0") is not None else float("nan")
        cur_alpha_ang = c.get("alpha", 0.0) * (BOHR2ANG ** 2) if c.get("alpha") is not None else float("nan")
        dr0 = cur_r0_ang - x.get("R0", 0.0)
        dalpha = cur_alpha_ang - x.get("alp", 0.0)
        rows.append((key, de, dpibo, dfc, dr0, dalpha, c, x))
    print(f"  total bond E: cur={tot_e_cur:.6f}  xtb={tot_e_xtb:.6f}  delta={tot_e_cur-tot_e_xtb:.6f} Eh")
    only_cur = sorted(set(cb) - set(xb))
    only_xtb = sorted(set(xb) - set(cb))
    if only_cur:
        print(f"  bonds only in curcuma: {only_cur[:10]}{' ...' if len(only_cur)>10 else ''} ({len(only_cur)})")
    if only_xtb:
        print(f"  bonds only in xtb: {only_xtb[:10]}{' ...' if len(only_xtb)>10 else ''} ({len(only_xtb)})")
    rows.sort(key=lambda r: abs(r[1]), reverse=True)
    print(f"  {'atoms':<10}{'dE(Eh)':>11}{'pibo_cur':>9}{'pibo_xtb':>9}{'dPibo':>8}{'dFc':>11}{'dR0(A)':>9}{'dAlp':>9}")
    for key, de, dpibo, dfc, dr0, dalpha, c, x in rows[:25]:
        pc = c.get("pibo", float("nan"))
        px = x.get("pibo", float("nan"))
        print(f"  {str(key):<10}{de:>+11.6f}{pc:>9.4f}{px:>9.4f}{dpibo:>+8.4f}{dfc:>+11.6f}{dr0:>+9.4f}{dalpha:>+9.4f}")


def main():
    args = sys.argv[1:]
    if not args or len(args) % 2 != 0:
        print("usage: python scripts/s30l_bond_compare.py SYS PART [SYS PART ...]")
        return
    for i in range(0, len(args), 2):
        compare(int(args[i]), args[i + 1])


if __name__ == "__main__":
    main()