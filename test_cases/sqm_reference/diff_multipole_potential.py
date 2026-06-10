#!/usr/bin/env python3
# AP6b Option(b) multipole-potential audit (Claude Generated).
#
# The integrals (dp_int/qp_int) are already proven bit-identical to tblite. This
# script compares the DOWNSTREAM multipole quantities, computed from tblite's
# converged density P:
#   - atomic moments  gf2.dpat / gf2.qpat   vs  dpat_tblite / qpat_tblite
#   - potentials      gf2.vdp  / gf2.vqp    vs  vdp_tblite  / vqp_tblite
#   - atom potential  gf2.vat_extra         vs  vat_tblite  (tblite vat is TOTAL;
#                                                 only clean for H2 where q_sh=0)
#
# For each quantity it reports max|ours-tbl| and, to expose a global sign flip,
# max|ours+tbl|. The first quantity (in the order dpat/qpat -> vdp/vqp -> vat)
# that diverges localizes the bug.
#
# Requires dumps from the patched tblite (USE_TBLITE build, e.g. release_tblite/).
# Usage:  diff_multipole_potential.py <dump.json> [...]

import json
import sys

import numpy as np


def col(d, key):
    return np.array(d[key]) if key in d else None


def cmp(name, O, T, note=""):
    if O is None or T is None:
        print(f"    {name:10s} : MISSING")
        return
    dminus = np.abs(O - T).max()
    dplus = np.abs(O + T).max()
    flip = "  <-- sign flip" if (dplus < dminus * 0.5 and dminus > 1e-9) else ""
    print(f"    {name:10s} : max|o-t|={dminus:.3e}  max|o+t|={dplus:.3e}{flip}  {note}")


def audit(path):
    d = json.load(open(path))
    mol = d.get("input_xyz", path)
    if "dpat_tblite" not in d or "gf2" not in d:
        print(f"{mol}: missing tblite potential keys (nat>6 or unpatched) -- skipped")
        return
    g = d["gf2"]
    print(f"{mol}  (nat={d['natoms']})")
    cmp("dpat", col(g, "dpat"), col(d, "dpat_tblite"))
    cmp("qpat", col(g, "qpat"), col(d, "qpat_tblite"))
    cmp("vdp", col(g, "vdp"), col(d, "vdp_tblite"))
    cmp("vqp", col(g, "vqp"), col(d, "vqp_tblite"))
    # vat: our gf2.vat_extra is the MULTIPOLE atom potential only; tblite's pot%vat
    # is multipole + the SELF-CONSISTENT D4 dispersion potential (dEdisp/dq, added
    # via dispersion%get_potential inside the SCF). The diff is therefore the D4 SCF
    # potential that native GFN2 omits (root cause of the AP6b ~1.5e-4 Eh F-shift).
    note = "(diff = self-consistent D4 dispersion potential; native GFN2 omits it)"
    cmp("vat", col(g, "vat_extra"), col(d, "vat_tblite"), note)


def main(argv):
    if len(argv) < 2:
        print("usage: diff_multipole_potential.py <dump.json> [...]")
        return 2
    for p in argv[1:]:
        audit(p)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
