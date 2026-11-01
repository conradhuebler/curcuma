#!/usr/bin/env python3
"""gCP/SRB reference data of HF-3c for every element H-Rn (test_cases/qm_hf3c/reference_gcp.json).

The oracle is otool_gcp, the gCP program (H. Kruse, J. G. Brandenburg, S. Grimme, version 1.06, Sep 2014) that
ORCA 6.1 calls for `! HF-3c`; it shares no code with curcuma. For every structure it is run twice:
  -level hf3c -grad     E_gCP, E_bas (the HF-3c short-range basis term, "SRB") and the gradient of their sum
  -level hf/minix -grad E_gCP alone and its gradient (so the two terms are checked separately)
otool_gcp prints the energies with 10 decimals (Eh) and the gradient with 19 (Eh/Bohr). Its analytic gradient is
compared with a central difference (1e-3 Angstrom) of its own printed energy: "otool_grad_consistent" is false where
they differ by more than 2e-7 Eh/Bohr (FD noise of the printed energy ~5e-8). Measured Nov 2026: for pairs with an atom without virtual MINIX functions
(He, Ne: the code skips the pair before adding the other atom's term) and for pairs whose Slater exponents differ by
less than 0.1 without being the same element (Cl-O, Co/Ni/Cu/As with O or H ..., the analytic overlap derivative takes
the closed form there while the energy takes the series), otool's gradient is not the derivative of its energy.

Structures (Angstrom): for every element X = H..Rn one 4-atom cluster X, H, O, X (pairs X-H, X-O, X-X, H-O, ...;
like elements exercise the series branch of the overlap) and a few mixed heavy-element clusters (pairs between the
1s/2s/3s shell classes and the effective-number groups Rb-Xe, Cs-La, Ce-Lu, Hf-Rn). Geometries are generated, not
molecules; gCP depends only on the pair distances and element numbers.

Usage: python3 scripts/hf3c_gcp_reference.py --otool /path/to/orca/otool_gcp [--work DIR]
Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import json
import os
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "test_cases", "qm_hf3c", "reference_gcp.json")
SYM = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
       "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
       "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()


def cluster(x):
    return [(x, 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.55), ("O", 1.85, 0.35, -0.45), (x, -1.30, 2.05, 0.75)]


MIXED = {
    "mix_light_heavy": [("C", 0, 0, 0), ("Na", 0, 0, 2.3), ("Kr", 2.4, 0.2, 0.1), ("He", -1.9, 1.1, 0.4),
                        ("Li", 0.6, -2.2, 0.3)],
    "mix_rb_xe": [("Rb", 0, 0, 0), ("I", 0, 0, 3.0), ("Ag", 2.6, 0.4, 0.2), ("Sn", -2.1, 1.7, -0.6), ("H", 0.4, -1.6, 1.1)],
    "mix_cs_rn": [("Cs", 0, 0, 0), ("La", 0, 0, 3.4), ("Ce", 3.0, 0.3, 0.5), ("Pt", -2.4, 1.9, 0.2), ("Pb", 0.7, -2.6, -0.8),
                  ("Rn", -0.9, -0.4, 3.1), ("F", 1.2, 1.4, 1.6)],
    "mix_tm": [("Fe", 0, 0, 0), ("Zn", 0, 0, 2.5), ("Cu", 2.2, 0.6, 0.3), ("Br", -2.0, 1.4, -0.5), ("Ga", 0.5, -2.3, 0.9),
               ("Mo", -1.0, -1.2, 2.6), ("S", 1.3, 1.1, 1.9)],
}


FD_STEP = 1.0e-3                   # Angstrom
BOHR_IN_ANGSTROM = 0.529177210903  # CODATA 2018


def write_xyz(path, name, atoms):
    with open(path, "w") as f:
        f.write("%d\n%s\n" % (len(atoms), name))
        for s, x, y, z in atoms:
            f.write("%-2s %14.8f %14.8f %14.8f\n" % (s, x, y, z))


def run(otool, xyz, level, work):
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(otool))
    out = subprocess.run([otool, xyz, "-level", level, "-grad"], cwd=work, env=env, capture_output=True,
                         text=True, check=True).stdout
    egcp = float(re.search(r"Egcp:\s+(\S+)", out).group(1))
    m = re.search(r"Ebas:\s+(\S+)", out)
    ebas = float(m.group(1)) if m else None
    g = out.split("Ggcp", 1)[1].split("cpu-time", 1)[0].split()
    return egcp, ebas, [[float(v) for v in g[i:i + 3]] for i in range(0, len(g), 3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--otool", required=True)
    ap.add_argument("--work")
    a = ap.parse_args()
    work = a.work or tempfile.mkdtemp(prefix="gcpref_")
    structs = {"Z%02d_%s" % (i + 1, s): cluster(s) for i, s in enumerate(SYM)}
    structs.update(MIXED)
    ver = subprocess.run([a.otool, "-h"], capture_output=True, text=True,
                         env=dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(a.otool))).stdout
    version = re.search(r"Version\s+(\S+\s+\S+\s+\S+)", ver).group(1)
    res = {"_about": "gCP/SRB reference of HF-3c (otool_gcp of ORCA 6.1), scripts/hf3c_gcp_reference.py",
           "_otool_gcp": version, "_units": "energies Eh, gradient Eh/Bohr, coordinates Angstrom", "structures": {}}
    for name, atoms in structs.items():
        xyz = os.path.join(work, name + ".xyz")
        write_xyz(xyz, name, atoms)
        e1, eb, g1 = run(a.otool, xyz, "hf3c", work)
        e2, _, g2 = run(a.otool, xyz, "hf/minix", work)
        assert abs(e1 - e2) < 1e-12, (name, e1, e2)
        # otool's analytic gradient against a central difference of its own printed energy (gCP + SRB)
        fd_dev = 0.0
        for i in range(len(atoms)):
            for k in range(3):
                ep = []
                for sgn in (1, -1):
                    moved = [list(t) for t in atoms]
                    moved[i][k + 1] += sgn * FD_STEP
                    write_xyz(xyz, name, moved)
                    ea, ebb, _ = run(a.otool, xyz, "hf3c", work)
                    ep.append(ea + ebb)
                fd = (ep[0] - ep[1]) / (2.0 * FD_STEP) * BOHR_IN_ANGSTROM
                fd_dev = max(fd_dev, abs(fd - g1[i][k]))
        write_xyz(xyz, name, atoms)
        res["structures"][name] = {"atoms": [list(t) for t in atoms], "gcp": e2, "srb": eb, "grad_gcp": g2,
                                   "grad_total": g1, "otool_grad_vs_own_fd": fd_dev,
                                   "otool_grad_consistent": fd_dev < 2e-7}
        print("%-16s gCP %.10f  SRB %.10f  |G - FD| %.1e" % (name, e2, eb, fd_dev))
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
