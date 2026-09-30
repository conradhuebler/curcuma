#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Interface-level gradient and frequency contract test.

Two things every consumer of EnergyCalculator relies on, neither of which was
covered before Sep 2026:

  1. `EnergyCalculator::Gradient()` returns **Eh/Angstrom**, for every method.
     SimpleMD feeds it straight into a Verlet integrator whose coordinates are in
     Angstrom, and the finite-difference Hessian differences it against Angstrom
     displacements. GFN-FF used to return Eh/Bohr here (a factor 1/au = 1.8897),
     which made MD forces too small; the existing `gfnff_numgrad_builtin` test
     could not see it because it compares GFN-FF's INTERNAL gradient against an
     internal finite difference, both in Bohr.

  2. Vibrational frequencies come out in cm^-1 and match the reference
     implementation. The Hessian was assembled in Eh/Ang^2 but consumed as
     Eh/Bohr^2, which made every frequency too high by 1/au; an empirical
     "+47.349 cm^-1" offset in the frequency formula hid part of it.

Test 1 needs nothing external. Test 2 compares against committed xtb 6.7.1
reference values (no xtb needed at run time).

Usage: check_gradient_units.py <curcuma-binary> [--verbose]
"""
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

AU = 0.52917721092          # Angstrom per Bohr

# A geometry that is not a stationary point for any of the three methods, so every
# component is non-trivial. Angstrom.
H2O = [("O", 0.000000, 0.000000, 0.000000),
       ("H", 0.000000, 0.000000, 0.980000),
       ("H", 0.950000, 0.000000, -0.240000)]

# The geometry the committed reference frequencies were produced at (an xtb GFN-FF
# optimum; not a stationary point for gfn1/gfn2, which is fine — both codes are
# evaluated at the SAME geometry).
H2O_OPT = [("O", 0.00000000705328, 0.00000000358300, 0.00202061038004),
           ("H", 0.00000610586735, 0.77614253823073, 0.58498834016127),
           ("H", -0.00000202185815, -0.77614049630927, 0.58499104944976)]

METHODS = ("gfn1", "gfn2", "gfnff")


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for el, x, y, z in atoms:
            f.write(f"{el} {x:.12f} {y:.12f} {z:.12f}\n")


def energy(binary, xyz, method, workdir):
    out = subprocess.run([binary, "-sp", str(xyz), "-method", method,
                          "-no_bmt", "-verbosity", "0"],
                         capture_output=True, text=True, cwd=workdir).stdout
    m = re.search(r"Single Point Energy = (-?\d+\.\d+)", out)
    return float(m.group(1)) if m else None


def analytic(binary, xyz, method, workdir):
    dump = Path(workdir) / "g.dump"
    if dump.exists():
        dump.unlink()
    subprocess.run([binary, "-sp", str(xyz), "-method", method, "-gradient",
                    "-dump_gradient", str(dump), "-no_bmt", "-verbosity", "0"],
                   capture_output=True, text=True, cwd=workdir)
    if not dump.exists():
        return None
    g = []
    for line in dump.read_text().splitlines():
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) == 3:
            g.append([float(v) for v in p])
    return g or None


def test_gradient_units(binary, workdir, verbose):
    """Analytic gradient vs central FD of the energy, displacements in Angstrom."""
    ok = True
    xyz = Path(workdir) / "m.xyz"
    write_xyz(xyz, H2O)
    h = 0.002
    for method in METHODS:
        g = analytic(binary, xyz, method, workdir)
        if g is None:
            print(f"  FAIL {method}: no analytic gradient")
            ok = False
            continue
        worst = 0.0
        worst_at = None
        for i in range(len(H2O)):
            for k in range(3):
                plus = [list(a) for a in H2O]
                plus[i][k + 1] += h
                write_xyz(Path(workdir) / "p.xyz", plus)
                ep = energy(binary, Path(workdir) / "p.xyz", method, workdir)
                minus = [list(a) for a in H2O]
                minus[i][k + 1] -= h
                write_xyz(Path(workdir) / "n.xyz", minus)
                em = energy(binary, Path(workdir) / "n.xyz", method, workdir)
                if ep is None or em is None:
                    print(f"  FAIL {method}: energy evaluation failed")
                    ok = False
                    continue
                fd = (ep - em) / (2 * h)            # Eh/Angstrom
                if abs(fd) < 1e-4:                  # skip near-zero components
                    continue
                rel = abs(g[i][k] - fd) / abs(fd)
                if rel > worst:
                    worst, worst_at = rel, (i, "xyz"[k], g[i][k], fd)
        # 1e-2 is loose against FD truncation but 100x tighter than the 0.89
        # relative error a Bohr/Angstrom mix-up produces.
        good = worst < 1e-2
        ok &= good
        tag = "ok  " if good else "FAIL"
        print(f"  [{tag}] {method:6s} max relative |analytic - FD| = {worst:.2e}")
        if (verbose or not good) and worst_at:
            i, c, ga, fd = worst_at
            print(f"         worst: atom {i} {c}: analytic {ga:+.8f}  FD {fd:+.8f} Eh/Ang"
                  f"   ratio {ga / fd:.6f} (1/au = {1 / AU:.6f}, au = {AU:.6f})")
    return ok


def test_frequencies(binary, workdir, verbose):
    """Vibrational frequencies vs committed xtb 6.7.1 values."""
    ok = True
    xyz = Path(workdir) / "opt.xyz"
    write_xyz(xyz, H2O_OPT)
    for method, ref in REFERENCE.items():
        out = subprocess.run([binary, "-hessian", str(xyz), "-method", method,
                              "-no_bmt", "-verbosity", "1"],
                             capture_output=True, text=True, cwd=workdir).stdout
        nums = None
        lines = out.splitlines()
        for i, line in enumerate(lines):
            if "Vibrational Frequencies" in line:
                for cand in lines[i + 1:i + 4]:
                    vals = re.findall(r"-?\d+\.\d+", re.sub(r"\x1b\[[0-9;]*m", "", cand))
                    vals = [float(v) for v in vals if abs(float(v)) > 1.0]
                    if len(vals) >= 3:
                        nums = vals[-3:]
                        break
                break
        if nums is None:
            print(f"  FAIL {method}: no frequencies parsed")
            ok = False
            continue
        dev = max(abs(a - b) / b for a, b in zip(nums, ref)) * 100.0
        good = dev < 2.0                    # percent
        ok &= good
        tag = "ok  " if good else "FAIL"
        print(f"  [{tag}] {method:6s} max deviation from xtb 6.7.1 = {dev:.2f} %")
        if verbose or not good:
            print(f"         curcuma {['%.1f' % v for v in nums]}  xtb {ref}")
    return ok


# xtb 6.7.1 reference frequencies (cm^-1) at H2O_OPT, `xtb x.xyz --<m> --hess`.
REFERENCE = {
    "gfn1": [1490.57, 3515.28, 3626.07],
    "gfn2": [1574.99, 3486.05, 3491.64],
    "gfnff": [1632.01, 3634.77, 3637.31],
}


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: check_gradient_units.py <curcuma-binary> [--verbose]")
    binary = sys.argv[1]
    verbose = "--verbose" in sys.argv
    workdir = tempfile.mkdtemp(prefix="gradunits_")
    print("Gradient unit contract (EnergyCalculator::Gradient() must be Eh/Angstrom):")
    ok1 = test_gradient_units(binary, workdir, verbose)
    print("Vibrational frequencies vs xtb 6.7.1:")
    ok2 = test_frequencies(binary, workdir, verbose)
    if ok1 and ok2:
        print("PASS")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
