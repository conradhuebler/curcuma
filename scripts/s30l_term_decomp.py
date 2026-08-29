#!/usr/bin/env python3
"""Term-by-term GFN-FF energy decomposition: curcuma vs xtb 6.6.1.

For each requested (system, part) pair, runs curcuma at verbosity 2 and xtb --gfnff --sp,
parses the energy-decomposition terms, prints side-by-side with deltas (cur - xtb) in Eh.

Usage:
    python scripts/s30l_term_decomp.py 12 A
    python scripts/s30l_term_decomp.py 11 A 12 A 15 A 16 A 23 AB
    python scripts/s30l_term_decomp.py all        # all outliers: 11/12/15/16 A+AB, 23 AB, 7/8 A+AB

Output is plain ASCII (no UTF), one block per (system, part).

Claude Generated (Jul 2026).
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN = os.path.join(ROOT, "test_cases", "s30l_test_set", "_run")
CURCUMA = os.path.join(ROOT, "release", "curcuma")
XTB = os.path.expanduser("~/Downloads/xtb-6.6.1/bin/xtb")
EH2KCAL = 627.509474


def curcuma_terms(xyz, charge):
    """Return dict term->Eh from curcuma verbosity-2 output."""
    cmd = [CURCUMA, "-sp", xyz, "-method", "gfnff", "-charge", str(charge),
           "-verbosity", "2", "-no_bmt"]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT).stdout
    # strip ANSI
    out = re.sub(r"\x1b\[[0-9;]*m", "", out)
    terms = {}
    total = None
    for line in out.splitlines():
        # [RESULT]  Bond                     -19.3445269927      0.00        --
        # [RESULT]  Total                    -18.0402939797   (no trailing CPU col)
        m = re.match(r"\[RESULT\]\s*(\w[\w ()]*?)\s+([+-]?\d+\.\d+)\s*(?:\d|$)", line)
        if m:
            name = m.group(1).strip()
            val = float(m.group(2))
            # skip timing rows ("Total energy call", "Total parameter generation", etc.)
            if "ms" in line or "wall=" in line or "%" in line:
                continue
            if name.lower() == "total":
                total = val
            else:
                terms[name] = val
    # curcuma splits repulsion into bonded+nonbond; combine for comparison
    rb = terms.pop("Repulsion (bonded)", None)
    rn = terms.pop("Repulsion (nonbond)", None)
    if rb is not None or rn is not None:
        terms["Repulsion"] = (rb or 0.0) + (rn or 0.0)
    # map curcuma names to canonical (match xtb)
    canon = {
        "Bond": "Bond", "Angle": "Angle", "Dihedral": "Torsion",
        "Inversion": "Inversion", "sTors": "sTors", "Dispersion": "Dispersion",
        "Repulsion": "Repulsion", "Coulomb": "Coulomb",
        "H-bonds": "HB", "X-bonds": "XB",
        "ATM (3-body)": "ATM", "BATM": "BATM",
    }
    terms = {canon.get(k, k): v for k, v in terms.items()}
    # xtb has no separate Inversion/sTors lines -> fold into Torsion for apples-to-apples
    torsion_combined = terms.get("Torsion", 0.0) + terms.get("Inversion", 0.0) + terms.get("sTors", 0.0)
    terms.pop("Inversion", None)
    terms.pop("sTors", None)
    terms["Torsion"] = torsion_combined
    return terms, total


def xtb_terms(xyz, charge):
    """Return dict term->Eh from xtb --gfnff --sp output (clean workdir with .CHRG)."""
    import tempfile, shutil
    workdir = tempfile.mkdtemp(prefix="s30l_xtb_")
    xyz_dst = os.path.join(workdir, os.path.basename(xyz))
    shutil.copy(xyz, xyz_dst)
    with open(os.path.join(workdir, ".CHRG"), "w") as f:
        f.write(f"{charge:+d}\n")
    cmd = [XTB, os.path.basename(xyz), "--gfnff", "--sp"]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir).stdout
    shutil.rmtree(workdir, ignore_errors=True)
    terms = {}
    total = None
    for line in out.splitlines():
        # :: bond energy              -18.617443809007 Eh    ::
        m = re.match(r"\s*::\s*(.+?)\s+(-?\d+\.\d+)\s+Eh\s*::", line)
        if m:
            name = m.group(1).strip().lower()
            val = float(m.group(2))
            if "total energy" in name:
                total = val
            elif "gradient" in name or "charge" in name or "summary" in name:
                continue
            else:
                terms[name] = val
    # map xtb names to canonical
    mapping = {
        "bond energy": "Bond",
        "angle energy": "Angle",
        "torsion energy": "Torsion",
        "repulsion energy": "Repulsion",
        "electrostat energy": "Coulomb",
        "dispersion energy": "Dispersion",
        "hb energy": "HB",
        "xb energy": "XB",
        "bonded atm energy": "BATM",
    }
    mapped = {}
    for k, v in terms.items():
        if k in mapping:
            mapped[mapping[k]] = v
    return mapped, total


# curcuma canonical term order (Inversion+sTors folded into Torsion to match xtb)
CUR_ORDER = ["Bond", "Angle", "Torsion", "Dispersion",
             "Repulsion", "Coulomb", "HB", "XB", "ATM", "BATM"]
# terms only xtb has / only curcuma has -> still show


def compare(sysid, part):
    base = os.path.join(RUN, str(sysid), part)
    xyz = base + ".xyz"
    if not os.path.exists(xyz):
        print(f"  (missing {xyz})")
        return
    # charge: read from original S30L set {i}/{part}/.CHRG (same source as harness)
    charge = 0
    pchrg = os.path.join(ROOT, "test_cases", "s30l_test_set", str(sysid), part, ".CHRG")
    if os.path.exists(pchrg):
        with open(pchrg) as f:
            try:
                charge = int(f.read().strip())
            except ValueError:
                charge = 0
    cur, cur_total = curcuma_terms(xyz, charge)
    xt, xtb_total = xtb_terms(xyz, charge)
    print(f"\n=== system {sysid} / {part}  (charge={charge}) ===")
    print(f"  curcuma total = {cur_total:.10f}   xtb total = {xtb_total:.10f}   "
          f"delta = {cur_total-xtb_total:+.10f} Eh ({(cur_total-xtb_total)*EH2KCAL:+.3f} kcal/mol)")
    print(f"  {'term':<14}{'curcuma':>18}{'xtb':>18}{'delta(Eh)':>16}{'delta(kcal)':>12}")
    allterms = list(dict.fromkeys(CUR_ORDER + list(xt.keys())))
    for t in allterms:
        c = cur.get(t)
        x = xt.get(t)
        if c is None and x is None:
            continue
        cs = f"{c:.10f}" if c is not None else "   --"
        xs = f"{x:.10f}" if x is not None else "   --"
        if c is not None and x is not None:
            d = c - x
            ds = f"{d:+.10f}"
            dk = f"{d*EH2KCAL:+.4f}"
        else:
            ds = "   --"
            dk = "   --"
        print(f"  {t:<14}{cs:>18}{xs:>18}{ds:>16}{dk:>12}")
    # unaccounted (cur total - sum cur listed) vs (xtb total - sum xtb listed)
    cur_sum = sum(v for k, v in cur.items() if k != "Total")
    xt_sum = sum(v for k, v in xt.items() if k != "Total")
    print(f"  {'(listed sum)':<14}{cur_sum:>18.10f}{xt_sum:>18.10f}")
    print(f"  {'(unlisted)':<14}{(cur_total-cur_sum):>18.10f}{(xtb_total-xt_sum):>18.10f}")


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return
    if args[0] == "all":
        pairs = []
        for s in [11, 12, 15, 16]:
            pairs += [(s, "A"), (s, "AB")]
        pairs += [(23, "AB"), (7, "A"), (8, "A"), (7, "AB"), (8, "AB")]
        pairs += [(30, "AB"), (25, "AB"), (26, "AB")]
    else:
        pairs = []
        i = 0
        while i + 1 < len(args):
            pairs.append((int(args[i]), args[i + 1]))
            i += 2
    for s, p in pairs:
        compare(s, p)


if __name__ == "__main__":
    main()