#!/usr/bin/env python3
"""
ORCA HF/def2-SVP reference generator for the WP1 1e-integral validation.

Runs ORCA 6.1 (`! HF def2-SVP`) on an XYZ and writes <stem>.orca_ref.json next
to the xyz. The WP1-valid cross-check stored here is the **smallest overlap
eigenvalue** (parsed from the ORCA output "Diagonalization of the overlap
matrix"): it is a 1e quantity, AO-order-invariant, and extremely basis-
sensitive, so it confirms that curcuma's def2-SVP.dat == ORCA's internal
def2-SVP. diff_qm_1e.py auto-discovers this file and compares it to curcuma's
smallest S eigenvalue.

The ORCA MO (orbital) energies are ALSO stored, labelled as Fock eigenvalues
(they are the SCF-converged HF Fock spectrum Hc + 2J - K, NOT the 1e Hcore
spectrum). Curcuma at WP1 has no 2e integrals / SCF, so the full MO-spectrum
comparison vs ORCA is a WP3+ deliverable and is NOT asserted by diff_qm_1e.py.

  qm_1e_reference.py <input.xyz> [--orca /opt/orca_6_1/orca]
                 [--basis def2-SVP] [--charge 0] [--spin 1] [--out DIR]

The ORCA input is built inline (`* xyz charge mult` with the xyz atoms in
Angstrom, ORCA's native xyz unit). mult = |2S+1| derived from --spin (number of
unpaired electrons): singlet=1, doublet=2, ... Default charge 0, spin 0 (singlet).

Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse, json, os, shutil, subprocess, sys, tempfile

SYM_TO_Z = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10}


def read_xyz(path):
    with open(path) as f:
        nat = int(f.readline())
        name = f.readline().strip()
        atoms = []
        for _ in range(nat):
            p = f.readline().split()
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return name, atoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--orca", default=os.environ.get("ORCA_BIN", "/opt/orca_6_1/orca"))
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0, help="number of unpaired electrons; mult=spin+1")
    ap.add_argument("--out", default=None, help="output dir (default: xyz dir)")
    ap.add_argument("--functional", default="HF")
    args = ap.parse_args()

    if not os.path.exists(args.orca):
        sys.stderr.write("ORCA not found at %s\n" % args.orca)
        return 2
    orca_2json = os.path.join(os.path.dirname(args.orca), "orca_2json")

    name, atoms = read_xyz(args.xyz)
    stem = os.path.splitext(os.path.basename(args.xyz))[0]
    out_dir = args.out or os.path.dirname(os.path.abspath(args.xyz))
    mult = args.spin + 1  # |2S+1| with S = spin/2

    tmp = tempfile.mkdtemp(prefix="qm1e_orca_")
    try:
        inp = os.path.join(tmp, stem + ".inp")
        with open(inp, "w") as f:
            f.write("! %s %s\n" % (args.functional, args.basis))
            f.write("* xyz %d %d\n" % (args.charge, mult))
            for sym, x, y, z in atoms:
                f.write("  %s %f %f %f\n" % (sym, x, y, z))
            f.write("*\n")

        # Run ORCA (cwd must be tmp so outputs land next to the inp).
        out_path = os.path.join(tmp, stem + ".out")
        with open(out_path, "w") as out:
            r = subprocess.run([args.orca, inp], cwd=tmp, stdout=out, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            sys.stderr.write("ORCA failed (%d); see %s.out\n" % (r.returncode, stem))
            return 1

        gbw = os.path.join(tmp, stem + ".gbw")
        if not os.path.exists(gbw):
            sys.stderr.write("ORCA produced no GBW (%s)\n" % gbw)
            return 1

        # Parse the smallest overlap eigenvalue from the ORCA output (the one
        # 1e quantity ORCA reliably prints; AO-invariant, basis-sensitive).
        s_min = None
        with open(out_path) as f:
            in_diag = False
            for line in f:
                if "Diagonalization of the overlap matrix" in line:
                    in_diag = True
                elif in_diag and "Smallest eigenvalue" in line:
                    try:
                        s_min = float(line.split()[-1])
                    except (ValueError, IndexError):
                        pass
                    break
        if s_min is None:
            sys.stderr.write("could not parse ORCA smallest overlap eigenvalue\n")
            return 1

        # orca_2json writes <stem>.json next to the gbw; progress goes to stdout.
        subprocess.run([orca_2json, gbw], cwd=tmp,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        js = os.path.join(tmp, stem + ".json")
        if not os.path.exists(js):
            sys.stderr.write("orca_2json produced no JSON (%s)\n" % js)
            return 1

        data = json.load(open(js))
        mo = data["Molecule"]["MolecularOrbitals"]["MOs"]
        energies = [float(m["OrbitalEnergy"]) for m in mo]

        ref = {
            "molecule": {"name": name, "natoms": len(atoms),
                          "atoms": [{"z": SYM_TO_Z[a[0]], "symbol": a[0]} for a in atoms]},
            "method": args.functional,
            "basis": args.basis,
            "charge": args.charge,
            "multiplicity": mult,
            "nmo": len(energies),
            "s_min_eigenvalue": s_min,          # WP1-valid 1e cross-check
            "mo_energies": energies,             # Fock (HF) spectrum, WP3+ use only
            "mo_energies_note": "SCF-converged Fock eigenvalues (Hc+2J-K), NOT 1e Hcore; WP3+ comparison",
        }
        out_path = os.path.join(out_dir, stem + ".orca_ref.json")
        with open(out_path, "w") as f:
            json.dump(ref, f, indent=2)
        print("wrote %s (%d MOs)" % (out_path, len(energies)))
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())