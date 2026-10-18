#!/usr/bin/env python3
"""Regression reference for ctest -L qm_lmp2 (test_cases/qm_lmp2/reference.json).

Claude Generated (Oct 2026, WP15). The local PNO-RI-MP2 (method lmp2) has no external program that reproduces it
to 1e-8 (other localisation, domains and PNO details in ORCA DLPNO-MP2 or the Python prototype), so the reference
is curcuma itself at default thresholds (since Oct 18 2026 the local fit; the global fit and small fit radii as
extra cases): a regression reference, not a validation. Each case runs
`curcuma -sp <xyz> -method lmp2 -qm.lmp2_check true` and stores E_corr(LMP2) and the canonical RI-MP2 correlation
energy of the same run (RI-MP2 itself is validated against PySCF in ctest -L qm_mp2). Tolerance 1e-7 Eh (the
iteration stops at |dE| < 1e-9, max|R| < 1e-7; the value is independent of the thread count by construction,
measured differences between 1 and 4 threads are reported by the ctest).

Usage: python3 scripts/lmp2_reference.py --curcuma release/curcuma [--out test_cases/qm_lmp2/reference.json]
"""
import argparse
import json
import os
import re
import subprocess
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
CASES = {
    # name: (library structure relative to the repository root, extra flags)
    "butane_fc": ("test_cases/structures/medium/butane.anti.xyz", ["-qm.frozen_core", "true"]),
    "butane_boys_fc": ("test_cases/structures/medium/butane.anti.xyz", ["-qm.frozen_core", "true",
                                                                       "-qm.lmp2_localizer", "boys"]),
    "benzene_fc": ("test_cases/structures/small/benzene.ideal-d6h.xyz", ["-qm.frozen_core", "true"]),
    "acetic_acid_dimer_fc": ("test_cases/structures/clusters/acetic-acid-dimer.cyclic-hb.xyz",
                             ["-qm.frozen_core", "true", "-qm.lmp2_pair_cut", "8"]),
    # global fit from the stored RI tensor (the default before Oct 18 2026)
    "global_butane_fc": ("test_cases/structures/medium/butane.anti.xyz", ["-qm.frozen_core", "true", "-qm.lmp2_fit", "global"]),
    "global_acetic_acid_dimer_fc": ("test_cases/structures/clusters/acetic-acid-dimer.cyclic-hb.xyz",
                                    ["-qm.frozen_core", "true", "-qm.lmp2_pair_cut", "8", "-qm.lmp2_fit", "global"]),
    # local fit with integral-direct 3-index integrals at radius 4 Bohr (truncation active), Oct 2026
    "local_pair_r4_acetic_acid_dimer_fc": ("test_cases/structures/clusters/acetic-acid-dimer.cyclic-hb.xyz",
                                           ["-qm.frozen_core", "true", "-qm.lmp2_pair_cut", "8", "-qm.lmp2_fit", "local",
                                            "-qm.lmp2_fit_domain", "pair", "-qm.lmp2_fit_radius", "4", "-qm.lmp2_ao_cut", "1e-6"]),
    "local_orbital_r4_acetic_acid_dimer_fc": ("test_cases/structures/clusters/acetic-acid-dimer.cyclic-hb.xyz",
                                              ["-qm.frozen_core", "true", "-qm.lmp2_pair_cut", "8", "-qm.lmp2_fit", "local",
                                               "-qm.lmp2_fit_domain", "orbital", "-qm.lmp2_fit_radius", "4", "-qm.lmp2_ao_cut", "1e-6"]),
}
TOL = 1e-7


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "test_cases/qm_lmp2/reference.json"))
    a = ap.parse_args()
    out = {"_about": "Regression reference of the native local PNO-RI-MP2 (method lmp2), produced by "
                     "scripts/lmp2_reference.py with curcuma itself; e_corr_rimp2 = canonical RI-MP2 of the same run "
                     "(same aux basis). Not an external validation.",
           "_settings": "def2-SVP spherical, -qm.scf_threshold 1e-10, aux def2-universal-JKFIT, default -qm.lmp2_* "
                        "unless listed in flags",
           "_units": "Hartree", "tolerance": TOL, "cases": {}}
    for name, (xyz, flags) in CASES.items():
        cmd = [a.curcuma, "-sp", os.path.join(ROOT, xyz), "-method", "lmp2", "-qm.basis", "def2-SVP",
               "-qm.scf_threshold", "1e-10", "-qm.threads", "4", "-verbosity", "2", "-no_bmt", "-gpu", "none",
               "-qm.lmp2_check", "true"] + flags
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3000)
        txt = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
        m1 = re.search(r"E\(LMP2 corr\)\s*=\s*(-?\d+\.\d+)", txt)
        m2 = re.search(r"E\(RI-MP2 corr\)\s*=\s*(-?\d+\.\d+)", txt)
        m3 = re.search(r"native lmp2 .*\((\d+) atoms", txt)
        if not (m1 and m2 and m3):
            print(txt[-3000:])
            sys.exit(f"{name}: values not found")
        e, ref, nat = float(m1.group(1)), float(m2.group(1)), int(m3.group(1))
        out["cases"][name] = {"xyz": xyz, "flags": flags, "natoms": nat, "e_corr_lmp2": e, "e_corr_rimp2": ref}
        print(f"{name:22s} E_corr LMP2 {e:.12f} RI-MP2 {ref:.12f} diff {e - ref:.3e} ({(e - ref) / nat:.3e} /atom)")
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
