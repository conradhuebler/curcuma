#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Contract test for the adaptive (step-rejecting) integrator.

Three things, on a SINGLE water molecule so the test costs a second:

  1. `adaptive_step` is off by default and changes nothing when it is off. The step
     wrapper (`SimpleMD::IntegratorStep`) must be a plain call to `Integrator()` in
     that case - the run has to be bit-identical to one that never knew about the
     feature.

  2. Without it, one water at 1750 K with dt = 1 fs in NVE DESTROYS itself: the
     GFN-FF O-H bond stiffens under compression (3817 cm^-1 at the equilibrium
     length, 10758 cm^-1 at 0.70 A), so the Verlet stability limit drops below the
     time step and the energy runs away. Measured: +8.78 Eh over 1 ps, with an
     average "temperature" of 1.8 million K. This is the failure mode the feature
     exists for, reproduced on three atoms.

  3. With it, the same run conserves energy (measured -0.003 Eh) without a single
     constraint and without touching the masses.

The initial velocities of a 3-atom system are deterministic here (verified over
three different seeds), so the numbers are reproducible.

See docs/MD_LARGE_SYSTEMS.md for the mechanism and the calibration measurements.

Usage: check_adaptive_step.py <curcuma-binary> [--verbose]
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# GFN-FF-optimised water, Angstrom.
H2O = [("O", 0.000000, -0.000000, 0.115197),
       ("H", 0.000000, 0.775927, -0.468149),
       ("H", 0.000000, -0.775927, -0.468149)]

VERBOSE = "--verbose" in sys.argv


def run_md(binary, workdir, label, extra):
    """One NVE run at 1750 K / dt = 1 fs. Returns (dE in Eh, average T in K)."""
    xyz = Path(workdir) / f"{label}.xyz"
    with open(xyz, "w") as f:
        f.write(f"{len(H2O)}\n\n")
        for el, x, y, z in H2O:
            f.write(f"{el} {x:.12f} {y:.12f} {z:.12f}\n")

    cmd = [binary, "-md", str(xyz), "-method", "gfnff", "-threads", "1",
           "-T", "1750", "-dt", "1.0", "-MaxTime", "1000", "-thermostat", "none",
           "-no_bmt", "-verbosity", "1", "-seed", "1"] + extra
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    if VERBOSE:
        print(" ".join(cmd))

    snaps = Path(workdir) / f"{label}.snapshots"
    first = snaps / f"{label}_step_0.json"
    last = snaps / f"{label}.final.json"
    if not last.exists():
        last = snaps / f"{label}.unstable.json"
    if not (first.exists() and last.exists()):
        print(f"FAIL: {label}: no snapshots written")
        if VERBOSE:
            print(proc.stdout[-3000:])
        return None, None

    a = json.load(open(first))["MD"]
    b = json.load(open(last))["MD"]
    dE = (b["average_Epot"] + b["average_Ekin"]) - (a["average_Epot"] + a["average_Ekin"])
    return dE, b["average_T"]


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    binary = sys.argv[1]
    failures = []

    with tempfile.TemporaryDirectory() as workdir:
        # 1 + 2: the default. Must be off, and must therefore blow up.
        dE_off, T_off = run_md(binary, workdir, "off", [])
        if dE_off is None:
            return 1
        print(f"  default (adaptive_step off): dE = {dE_off:+.4f} Eh, <T> = {T_off:.0f} K")
        if dE_off < 1.0:
            failures.append(
                f"the reference failure mode did not reproduce: dE = {dE_off:+.4f} Eh, "
                "expected > +1 Eh. Either GFN-FF changed, or adaptive_step is no longer "
                "off by default - check which, the second would be a silent default change.")

        # An explicit false must behave exactly like the default.
        dE_false, T_false = run_md(binary, workdir, "false", ["-adaptive_step", "false"])
        if dE_false is None:
            return 1
        if abs(dE_false - dE_off) > 1e-10 or abs(T_false - T_off) > 1e-8:
            failures.append(
                f"-adaptive_step false is not identical to the default: "
                f"dE {dE_false:+.6f} vs {dE_off:+.6f} Eh, <T> {T_false:.4f} vs {T_off:.4f} K")

        # 3: with the feature, the same run has to conserve energy.
        dE_on, T_on = run_md(binary, workdir, "on", ["-adaptive_step", "true"])
        if dE_on is None:
            return 1
        print(f"  adaptive_step on           : dE = {dE_on:+.4f} Eh, <T> = {T_on:.0f} K")
        if abs(dE_on) > 0.05:
            failures.append(
                f"step rejection did not conserve the energy: dE = {dE_on:+.4f} Eh, "
                "expected |dE| < 0.05 Eh")
        if T_on > 5000.0:
            failures.append(
                f"step rejection did not keep the trajectory bounded: <T> = {T_on:.0f} K")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("  OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
