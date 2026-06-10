#!/usr/bin/env python3
# AP6b multipole-integral audit (Claude Generated).
#
# Compares curcuma's native GFN2 multipole AO integrals (dp_int/qp_int, our
# production cgto_multipole path) against the GENUINE tblite integrals exposed
# by the diagnostic tblite patch (keys dipole_integral_tblite /
# quadrupole_integral_tblite).
#
# The dump_tblite_multipole tool (USE_TBLITE build, e.g. release_tblite/) must
# have been run with the patched tblite to populate the *_tblite keys. The
# integral arrays are only written for nat<=6 systems.
#
# Usage:  diff_multipole_ints.py <dump1.json> [dump2.json ...]
# Exit 0 if every molecule matches to <1e-12, else 1.

import json
import sys

import numpy as np

TOL = 1e-12


def audit(path):
    d = json.load(open(path))
    mol = d.get("input_xyz", path)
    if "dipole_integral_tblite" not in d:
        print(f"{mol}: no tblite integral keys (nat>6 or unpatched tblite) -- skipped")
        return None
    odp = np.array(d["gf2"]["dp_int"])
    oqp = np.array(d["gf2"]["qp_int"])
    tdp = np.array(d["dipole_integral_tblite"])
    tqp = np.array(d["quadrupole_integral_tblite"])
    nao = d["norbitals"]
    ao2at = [d["shell_map"][d["orbital_map"][mu]] for mu in range(nao)]
    same = np.array([[ao2at[i] == ao2at[j] for j in range(nao)] for i in range(nao)])

    def blkmax(O, T, mask):
        if not mask.any():
            return 0.0
        return max(np.abs(O[k] - T[k])[mask].max() for k in range(O.shape[0]))

    res = {
        "dp_same": blkmax(odp, tdp, same), "dp_off": blkmax(odp, tdp, ~same),
        "qp_same": blkmax(oqp, tqp, same), "qp_off": blkmax(oqp, tqp, ~same),
    }
    worst = max(res.values())
    status = "PASS" if worst < TOL else "FAIL"
    print(f"{mol:24s} nao={nao:2d}  "
          f"dp[same={res['dp_same']:.2e} off={res['dp_off']:.2e}]  "
          f"qp[same={res['qp_same']:.2e} off={res['qp_off']:.2e}]  {status}")
    return worst < TOL


def main(argv):
    if len(argv) < 2:
        print("usage: diff_multipole_ints.py <dump.json> [...]")
        return 2
    results = [audit(p) for p in argv[1:]]
    checked = [r for r in results if r is not None]
    return 0 if checked and all(checked) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
