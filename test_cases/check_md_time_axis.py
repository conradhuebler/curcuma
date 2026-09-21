#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Contract test for the MD time axis: one reported femtosecond must be one femtosecond.

SimpleMD integrates in Angstrom / amu / Hartree. The time unit implied by that
combination is sqrt(amu * A^2 / Eh) = 1.9516 fs, NOT one femtosecond, so the step
the user asks for has to be converted before it multiplies a velocity. Until
Sep 2026 it was not, and every curcuma MD ran with a step 1.9516x larger than
requested - a default `-dt 1.0` was really 1.95 fs. Nothing caught it, because the
integrator was internally consistent: energies, forces, temperature and energy
conservation were all correct. Only the clock was wrong.

The observable that does catch it is a vibrational period, because it ties the time
axis to a frequency that is measured independently of MD:

  1. `-hessian` gives the harmonic frequencies of water. That path is validated
     against xtb 6.7.1 to 0.13 % (Known Issue #28) and knows nothing about MD.
  2. Displacing one O-H and running NVE at ~0 K makes the molecule oscillate in
     that bond's local mode, whose frequency for a symmetric AB2 molecule is
     sqrt((w_sym^2 + w_anti^2) / 2).
  3. The period measured from the trajectory, in the time the MD reports, must
     equal 1 / (c * w_local).

A mismatch here is a unit error in the integrator and nothing else. On the
pre-fix binary this test fails by a factor 1.95.

Usage: check_md_time_axis.py <curcuma-binary> [--verbose]
"""
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

VERBOSE = "--verbose" in sys.argv

# Water, near the GFN-FF minimum. Angstrom.
H2O = [("O", 0.000000, 0.000000, 0.115197),
       ("H", 0.000000, 0.775927, -0.468149),
       ("H", 0.000000, -0.775927, -0.468149)]

C_CM_PER_S = 2.99792458e10
TOLERANCE_PERCENT = 3.0   # generous: this test is looking for a factor of two


def write_xyz(path, atoms, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for el, x, y, z in atoms:
            f.write(f"{el} {x:.10f} {y:.10f} {z:.10f}\n")


def run(binary, args, workdir):
    proc = subprocess.run([binary] + args, capture_output=True, text=True, cwd=workdir)
    out = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout)
    if VERBOSE:
        print(" ".join([binary] + args))
    return out


def frequencies(binary, workdir, xyz):
    """The two O-H stretch frequencies from the Hessian, in cm^-1."""
    out = run(binary, ["-hessian", xyz, "-method", "gfnff", "-threads", "1",
                       "-no_bmt", "-verbosity", "1"], workdir)
    lines = out.split("\n")
    for i, line in enumerate(lines):
        if "Vibrational Frequencies" in line:
            values = []
            for follow in lines[i + 1:i + 5]:
                values += [float(v) for v in re.findall(r"(?<![\d.])\d+\.\d+(?![\d.]*\()", follow)]
            values = sorted(v for v in values if v > 1.0)
            if len(values) >= 2:
                return values[-2], values[-1]
    return None


def md_period(binary, workdir, xyz):
    """Period of the O-H oscillation, in the time the MD itself reports, in fs."""
    label = Path(xyz).stem
    run(binary, ["-md", xyz, "-method", "gfnff", "-threads", "1", "-T", "0.00001",
                 "-dt", "0.05", "-MaxTime", "60", "-thermostat", "none", "-seed", "1",
                 "-no_bmt", "-verbosity", "1", "-dump", "1", "-print", "100000"], workdir)
    trj = Path(workdir) / f"{label}.snapshots" / f"{label}.trj.xyz"
    if not trj.exists():
        return None
    lines = open(trj).read().split("\n")
    n = int(lines[0].split()[0])
    r, t = [], []
    for start in range(0, len(lines) - n - 1, n + 2):
        coords = []
        for line in lines[start + 2:start + 2 + n]:
            parts = line.split()
            if len(parts) < 4:
                break
            coords.append(tuple(float(v) for v in parts[1:4]))
        if len(coords) < n:
            break
        r.append(math.dist(coords[0], coords[1]))
        t.append(float(lines[start + 1].split()[0]))
    if len(r) < 20:
        return None
    mean = sum(r) / len(r)
    crossings = [i for i in range(1, len(r)) if (r[i - 1] - mean) * (r[i] - mean) < 0]
    if len(crossings) < 4:
        return None
    def crossing_time(i):
        f = (mean - r[i - 1]) / (r[i] - r[i - 1])
        return t[i - 1] + f * (t[i] - t[i - 1])
    times = [crossing_time(i) for i in crossings]
    # two crossings per period
    return sum(times[k + 2] - times[k] for k in range(len(times) - 2)) / (len(times) - 2)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    binary = sys.argv[1]

    with tempfile.TemporaryDirectory() as workdir:
        ref = Path(workdir) / "ref.xyz"
        write_xyz(ref, H2O)
        freqs = frequencies(binary, workdir, str(ref))
        if freqs is None:
            print("FAIL: could not read the vibrational frequencies from -hessian")
            return 1
        w_sym, w_anti = freqs
        w_local = math.sqrt((w_sym ** 2 + w_anti ** 2) / 2.0)
        expected = 1.0e15 / (w_local * C_CM_PER_S)

        # stretch ONE O-H by 0.006 A to excite its local mode
        stretched = list(H2O)
        ox = H2O[0][1:]
        hy = H2O[1][1:]
        d = [hy[k] - ox[k] for k in range(3)]
        norm = math.sqrt(sum(c * c for c in d))
        stretched[1] = ("H",) + tuple(hy[k] + 0.006 * d[k] / norm for k in range(3))
        md = Path(workdir) / "osc.xyz"
        write_xyz(md, stretched)

        measured = md_period(binary, workdir, str(md))
        if measured is None:
            print("FAIL: could not measure an oscillation period from the trajectory")
            return 1

        deviation = 100.0 * (measured / expected - 1.0)
        print(f"  Hessian: {w_sym:.1f} and {w_anti:.1f} cm-1, local mode {w_local:.1f} cm-1")
        print(f"  period expected {expected:.4f} fs, MD reports {measured:.4f} fs "
              f"({deviation:+.2f} %)")

        if abs(deviation) > TOLERANCE_PERCENT:
            print("\nFAILED:")
            print(f"  - the MD clock is off by a factor {expected / measured:.4f}. The integrator's "
                  "own time unit is")
            print("    sqrt(amu*A^2/Eh) = 1.9516 fs (CurcumaUnit::Constants::MD_TIME_UNIT_FS); a "
                  "user-supplied")
            print("    step in femtoseconds must be multiplied by FS_TO_MD_TIME before it "
                  "multiplies a")
            print("    velocity. Check SimpleMD::Verlet(), Rattle() and NoseHover().")
            return 1
        print("  OK")
        return 0


if __name__ == "__main__":
    sys.exit(main())
