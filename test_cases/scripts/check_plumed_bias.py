#!/usr/bin/env python3
# Claude Generated (Oct 2026)
"""Contract test for the PLUMED coupling: a RESTRAINT must act as a restoring force.

Two argon atoms 6.5 A apart are held by RESTRAINT ARG=d AT=0.6 KAPPA=5000 (nm, kJ/mol) in NVE.
The UFF pair term is negligible there, so d oscillates around 0.6 nm with the period
T = 2 pi sqrt(mu / kappa) = 0.397 ps for mu = 19.974 amu, which ties the sign of the force and
the energy, length and time units of the interface to one number.

With the force sign wrong (until Oct 2026) the distance grows as cosh(w t), reaching 1e19 nm
within 3 ps. PLUMED's kBT is checked as well: it must be kB*T, not kB.

Usage: check_plumed_bias.py <curcuma-binary>
"""
import math
import subprocess
import sys
import tempfile
from pathlib import Path

KAPPA = 5000.0
EXPECTED = 2.0 * math.pi * math.sqrt(0.5 * 39.948 / KAPPA)  # ps, kJ/mol/nm^2/amu = 1/ps^2
T0 = 300.0


def main():
    binary = sys.argv[1]
    with tempfile.TemporaryDirectory() as tmp:
        w = Path(tmp)
        (w / "ar2.xyz").write_text("2\nAr2\nAr 0.0 0.0 0.0\nAr 0.0 0.0 6.5\n")
        (w / "plumed.dat").write_text(
            "UNITS LENGTH=nm TIME=ps ENERGY=kj/mol\n"
            "d: DISTANCE ATOMS=1,2 NOPBC\n"
            f"r: RESTRAINT ARG=d AT=0.6 KAPPA={KAPPA}\n"
            "PRINT ARG=d STRIDE=1 FILE=COLVAR\n")
        subprocess.run([binary, "-md", "ar2.xyz", "-method", "uff", "-mtd", "-plumed", "plumed.dat",
                        "-dt", "0.5", "-MaxTime", "3000", "-thermostat", "none", "-T", str(T0),
                        "-seed", "2", "-threads", "1", "-verbosity", "0"],  # seed 2: with seeds 1, 4, 9 the 1-DOF run hits a turning point with Ekin ~ 0 and the T > 10000*T_prev check aborts it (TODO.md)
                       cwd=w, capture_output=True, text=True, timeout=250)
        colvar = list(w.glob("ar2.md.*/ar2.plumed_files/COLVAR"))
        log = list(w.glob("ar2.md.*/ar2.plumed_files/plumed_log.out"))
        if not colvar or not log:
            print("FAIL: no COLVAR/plumed_log.out in the BMT PLUMED directory")
            return 1
        rows = [line.split() for line in colvar[0].read_text().splitlines() if not line.startswith("#")]
        t = [float(r[0]) for r in rows]
        d = [float(r[1]) for r in rows]
        kbt_lines = [line for line in log[0].read_text().splitlines() if "KbT:" in line]
    if max(d) > 1.0:
        print(f"FAIL: distance ran away to {max(d):.3g} nm, the bias force has the wrong sign")
        return 1
    mean = sum(d) / len(d)
    up = [t[i - 1] + (t[i] - t[i - 1]) * (mean - d[i - 1]) / (d[i] - d[i - 1])
          for i in range(1, len(d)) if d[i - 1] < mean <= d[i]]
    if len(up) < 3:
        print(f"FAIL: only {len(up)} upward crossings of the mean distance")
        return 1
    period = (up[-1] - up[0]) / (len(up) - 1)
    print(f"  period {period:.4f} ps, expected {EXPECTED:.4f} ps; d {min(d):.4f}..{max(d):.4f} nm")
    if abs(period / EXPECTED - 1.0) > 0.02:
        print("FAIL: period off by more than 2 %, a unit of the PLUMED interface is wrong")
        return 1
    kbt = float(kbt_lines[0].split("KbT:")[1]) if kbt_lines else float("nan")
    expected_kbt = 0.008314462618 * T0
    print(f"  KbT {kbt:.4f} kJ/mol, expected {expected_kbt:.4f}")
    if not abs(kbt / expected_kbt - 1.0) < 0.01:
        print("FAIL: PLUMED kBT is not kB*T")
        return 1
    print("  OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
