#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Contract test for the MD time axis: one reported femtosecond must be one femtosecond.

WHY THIS TEST EXISTS
--------------------
SimpleMD integrates in Angstrom / amu / Hartree. Velocities are drawn as
sqrt(kb_Eh * T / m), i.e. in sqrt(Eh/amu), and the Verlet position update is

    x[A] += dT * v  -  0.5 * g[Eh/A] * (1/m)[1/amu] * dT^2

so for both terms to produce an Angstrom the step must carry the unit

    A * sqrt(amu/Eh) = sqrt(amu * A^2 / Eh) = 1.9516144 fs,

NOT one femtosecond (CurcumaUnit::Constants::MD_TIME_UNIT_FS). Until Sep 2026 the
user's `-md.time_step X`, meant as X femtoseconds, was handed to the integrator
unconverted, so every curcuma MD really ran at X * 1.9516144 fs and `-MaxTime` was
stretched by the same factor.

NOTHING ELSE CAN CATCH IT. The integrator was internally consistent - both terms
demanded the same unit - so energies, forces, temperature and NVE energy conservation
were all correct; a consistent integrator conserves energy in ANY time unit. Only an
observable that ties the MD clock to a frequency measured OUTSIDE the MD can see it.

THE MEASUREMENT
---------------
  1. `-opt` relaxes water with the method under test, so the Hessian is taken at a
     stationary point (otherwise the "frequencies" are contaminated by the residual
     gradient and this test measures the geometry, not the clock).
  2. `-hessian` gives the two O-H stretch frequencies. That path is validated against
     xtb 6.7.1 to 0.13 % (CLAUDE.md Known Issue #28) and knows nothing about MD.
  3. For a symmetric AB2 molecule the LOCAL O-H stretch frequency is
     sqrt((w_sym^2 + w_anti^2) / 2). Displacing one O-H by 0.006 A and running NVE at
     ~0 K excites exactly that local mode.
  4. The period read off the trajectory, in the time the MD itself reports, must equal
     1 / (c * w_local).

On the pre-fix binary this fails by the factor 1.9516 for every method (measured:
gfnff 1.9536, gfn2 1.9520, gfn1 1.9511 - method-independent, because the integrator
only ever sees a gradient in Eh/Angstrom).

Usage: check_md_time_axis.py <curcuma-binary> [--verbose] [--methods gfnff,gfn2]
"""
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

VERBOSE = "--verbose" in sys.argv

# Water, near the GFN-FF minimum; every method re-optimises it before measuring.
H2O = [("O", 0.000000, 0.000000, 0.115197),
       ("H", 0.000000, 0.775927, -0.468149),
       ("H", 0.000000, -0.775927, -0.468149)]

C_CM_PER_S = 2.99792458e10
TOLERANCE_PERCENT = 2.0  # generous: the error this test hunts is a factor of two
MD_STEP_FS = 0.05
MD_TIME_FS = 60


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write("%d\n\n" % len(atoms))
        for el, x, y, z in atoms:
            f.write("%s %.10f %.10f %.10f\n" % (el, x, y, z))


def read_last_frame(path):
    lines = open(path).read().split("\n")
    n = int(lines[0].split()[0])
    atoms, block = [], None
    for start in range(0, len(lines) - n - 1, n + 2):
        rows = lines[start + 2:start + 2 + n]
        if all(len(r.split()) >= 4 for r in rows):
            block = rows
    if block is None:
        return None
    for row in block:
        p = row.split()
        atoms.append((p[0],) + tuple(float(v) for v in p[1:4]))
    return atoms


def run(binary, args, workdir):
    proc = subprocess.run([binary] + args, capture_output=True, text=True, cwd=workdir)
    if VERBOSE:
        print("  $ " + " ".join([Path(binary).name] + args))
    return re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout)


def optimise(binary, method, workdir, xyz):
    # lbfgs, not the default 'auto': on a perfect-symmetry system whose only active
    # mode is totally symmetric, 'auto' can abort with "Energy rise exceeded maximum
    # allowed" (CLAUDE.md Known Issue #28, pre-existing and unrelated to this test).
    run(binary, ["-opt", xyz, "-method", method, "-threads", "1", "-no_bmt",
                 "-opt.optimizer", "lbfgs", "-verbosity", "1"], workdir)
    out = Path(workdir) / (Path(xyz).stem + ".opt.xyz")
    return read_last_frame(out) if out.exists() else None


def stretch_frequencies(binary, method, workdir, xyz):
    """The two O-H stretch frequencies from the Hessian, in cm^-1."""
    out = run(binary, ["-hessian", xyz, "-method", method, "-threads", "1",
                       "-no_bmt", "-verbosity", "1"], workdir)
    lines = out.split("\n")
    for i, line in enumerate(lines):
        if "Vibrational Frequencies" in line:
            values = []
            for follow in lines[i + 1:i + 6]:
                values += [float(v) for v in re.findall(r"(?<![\d.])\d+\.\d+(?![\d.]*\()", follow)]
            values = sorted(v for v in values if v > 1.0)
            if len(values) >= 2:
                return values[-2], values[-1]
    return None


def md_period(binary, method, workdir, xyz):
    """Period of the O-H oscillation, in the time the MD itself reports, in fs."""
    label = Path(xyz).stem
    run(binary, ["-md", xyz, "-method", method, "-threads", "1", "-T", "0.00001",
                 "-dt", str(MD_STEP_FS), "-MaxTime", str(MD_TIME_FS),
                 "-thermostat", "none", "-seed", "1", "-no_bmt", "-verbosity", "1",
                 "-dump", "1", "-print", "1000000"], workdir)
    trj = Path(workdir) / (label + ".snapshots") / (label + ".trj.xyz")
    if not trj.exists():
        return None
    lines = open(trj).read().split("\n")
    n = int(lines[0].split()[0])
    r, t = [], []
    for start in range(0, len(lines) - n - 1, n + 2):
        coords = []
        for line in lines[start + 2:start + 2 + n]:
            p = line.split()
            if len(p) < 4:
                break
            coords.append(tuple(float(v) for v in p[1:4]))
        if len(coords) < n:
            break
        r.append(math.dist(coords[0], coords[1]))
        t.append(float(lines[start + 1].split()[0]))
    if len(r) < 20:
        return None
    mean = sum(r) / len(r)
    cross = [i for i in range(1, len(r)) if (r[i - 1] - mean) * (r[i] - mean) < 0]
    if len(cross) < 4:
        return None

    def crossing_time(i):
        f = (mean - r[i - 1]) / (r[i] - r[i - 1])
        return t[i - 1] + f * (t[i] - t[i - 1])

    ts = [crossing_time(i) for i in cross]
    # two zero crossings per period
    return sum(ts[k + 2] - ts[k] for k in range(len(ts) - 2)) / (len(ts) - 2)


def check(binary, method, workdir):
    start = Path(workdir) / ("start_%s.xyz" % method)
    write_xyz(start, H2O)
    geom = optimise(binary, method, workdir, str(start))
    if geom is None:
        return None, "could not optimise water with -method %s" % method

    ref = Path(workdir) / ("ref_%s.xyz" % method)
    write_xyz(ref, geom)
    freqs = stretch_frequencies(binary, method, workdir, str(ref))
    if freqs is None:
        return None, "could not read the vibrational frequencies from -hessian"
    w_sym, w_anti = freqs
    w_local = math.sqrt((w_sym ** 2 + w_anti ** 2) / 2.0)
    expected = 1.0e15 / (w_local * C_CM_PER_S)

    # stretch ONE O-H by 0.006 A to excite its local mode
    disp = list(geom)
    o, h = geom[0][1:], geom[1][1:]
    d = [h[k] - o[k] for k in range(3)]
    norm = math.sqrt(sum(c * c for c in d))
    disp[1] = ("H",) + tuple(h[k] + 0.006 * d[k] / norm for k in range(3))
    osc = Path(workdir) / ("osc_%s.xyz" % method)
    write_xyz(osc, disp)

    measured = md_period(binary, method, workdir, str(osc))
    if measured is None:
        return None, "could not measure an oscillation period from the trajectory"

    deviation = 100.0 * (measured / expected - 1.0)
    print("  %-7s local mode %7.1f cm-1 (%.1f / %.1f) | period expected %.4f fs, "
          "MD reports %.4f fs (%+.2f %%)"
          % (method, w_local, w_sym, w_anti, expected, measured, deviation))
    return (deviation, expected / measured), None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    binary = sys.argv[1]
    methods = ["gfnff", "gfn2"]
    for i, a in enumerate(sys.argv):
        if a == "--methods" and i + 1 < len(sys.argv):
            methods = sys.argv[i + 1].split(",")

    failures = []
    with tempfile.TemporaryDirectory() as workdir:
        for method in methods:
            result, err = check(binary, method, workdir)
            if result is None:
                print("  FAIL %-7s %s" % (method, err))
                failures.append((method, err))
                continue
            deviation, ratio = result
            if abs(deviation) > TOLERANCE_PERCENT:
                failures.append((method, "clock off by a factor %.4f" % ratio))

    if failures:
        print("\nFAILED:")
        for method, why in failures:
            print("  - %s: %s" % (method, why))
        print("\n  The integrator's own time unit is sqrt(amu*A^2/Eh) = 1.9516144 fs")
        print("  (CurcumaUnit::Constants::MD_TIME_UNIT_FS). A user-supplied step in")
        print("  femtoseconds must be multiplied by FS_TO_MD_TIME before it multiplies a")
        print("  velocity. Check SimpleMD::Verlet(), Rattle() and NoseHover().")
        return 1
    print("  OK: the MD clock agrees with the Hessian frequency for %s" % ", ".join(methods))
    return 0


if __name__ == "__main__":
    sys.exit(main())
