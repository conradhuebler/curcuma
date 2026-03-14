#!/usr/bin/env python3
"""
Test Bond Energy Calculation

Simple test: Calculate bond stretching energy using E = k * (r - r0)^2
and compare with XTB reference values from .out files.

Formula from gfnff_engrad.F90 lines 373-399:
  E_bond = vbond(3) * (r_actual - vbond(1))^2
"""

import json
import math
from pathlib import Path

# Unit conversion
ANGSTROM_TO_BOHR = 1.889726
HARTREE_TO_EV = 27.211386

# Test data
test_cases = {
    "HH": {
        "xyz_file": "molecules/dimers/HH.xyz",
        "json_file": "molecules/dimers/HH.json",
        "bond_energy_reference": -0.164952024621,  # From HH.out
    },
    "OH": {
        "xyz_file": "molecules/dimers/OH.xyz",
        "json_file": "molecules/dimers/OH.json",
        "bond_energy_reference": -0.170910705515,  # From OH.out
    },
    "HCl": {
        "xyz_file": "molecules/dimers/HCl.xyz",
        "json_file": "molecules/dimers/HCl.json",
        "bond_energy_reference": -0.084310498873,  # From HCl.out
    },
}


def read_xyz(filename):
    """Read XYZ file and return coordinates in Bohr"""
    with open(filename) as f:
        lines = f.readlines()

    natoms = int(lines[0].strip())
    coords = []

    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        x = float(parts[1]) * ANGSTROM_TO_BOHR
        y = float(parts[2]) * ANGSTROM_TO_BOHR
        z = float(parts[3]) * ANGSTROM_TO_BOHR
        coords.append((x, y, z))

    return coords


def read_vbond_json(filename):
    """Read vbond parameters from JSON"""
    with open(filename) as f:
        data = json.load(f)
    return data["vbond"], data["blist"]


def distance(c1, c2):
    """Calculate distance between two points in Bohr"""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def calculate_bond_energy(coords, vbond_params, bond_list):
    """
    Calculate bond stretching energy - EXACT FORTRAN FORMULA

    Formula from gfnff_engrad.F90 line 699:
      E_bond = vbond(3) * exp(-vbond(2) * dr²)

    where:
      r = actual distance [Bohr]
      r0 = vbond(1) equilibrium offset [Bohr]
      dr = r - r0
      vbond(2) = exponential decay coefficient
      vbond(3) = energy amplitude [Hartree]
    """
    total_energy = 0.0
    energies = []

    for bond_idx, (i, j) in enumerate(bond_list):
        # Convert from 1-based to 0-based indexing
        i -= 1
        j -= 1

        # Get bond parameters from JSON
        # vbond(1) = r0 (equilibrium distance), vbond(2) = t8 (decay), vbond(3) = amplitude
        r0, t8, amplitude = vbond_params[bond_idx]

        # Calculate actual distance in Bohr
        r_actual = distance(coords[i], coords[j])

        # EXACT formula from gfnff_engrad.F90 line 699:
        # E_bond = vbond(3) * exp(-vbond(2) * dr²)
        dr = r_actual - r0
        E_bond = amplitude * math.exp(-t8 * dr * dr)

        energies.append({
            "bond": (i, j),
            "r_actual": r_actual,
            "r0": r0,
            "t8": t8,
            "amplitude": amplitude,
            "dr": dr,
            "E_bond": E_bond,
        })

        total_energy += E_bond

    return total_energy, energies


def main():
    print("\n" + "="*80)
    print("GFN-FF Bond Energy Test")
    print("="*80)

    results = []

    for name, data in test_cases.items():
        print(f"\n{name}:")
        print("-" * 40)

        # Read geometry and parameters
        coords = read_xyz(data["xyz_file"])
        vbond_params, bond_list = read_vbond_json(data["json_file"])

        # Calculate energy
        total_energy, energies = calculate_bond_energy(coords, vbond_params, bond_list)

        reference = data["bond_energy_reference"]
        error = abs(total_energy - reference)
        error_pct = 100.0 * error / abs(reference)

        print(f"  Calculated:  {total_energy:20.12f} Hartree")
        print(f"  Reference:   {reference:20.12f} Hartree")
        print(f"  Error:       {error:20.12f} Hartree ({error_pct:.4f}%)")
        print(f"  Match:       {'✓' if error < 0.001 else '✗'}")

        # Bond details
        for e in energies:
            print(f"\n  Bond ({e['bond'][0]},{e['bond'][1]}):")
            print(f"    r_actual  = {e['r_actual']:.8f} Bohr")
            print(f"    r0        = {e['r0']:.8f} Bohr")
            print(f"    dr        = {e['dr']:.8f} Bohr")
            print(f"    t8        = {e['t8']:.12f} (exponential decay)")
            print(f"    amplitude = {e['amplitude']:.12f} Hartree")
            print(f"    E_bond    = {e['E_bond']:.12f} Hartree")

        results.append({
            "system": name,
            "calculated": total_energy,
            "reference": reference,
            "error": error,
            "error_pct": error_pct,
        })

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\n{'System':<10} {'Calculated':>20} {'Reference':>20} {'Error %':>12}")
    print("-" * 62)

    total_error = 0
    for r in results:
        print(f"{r['system']:<10} {r['calculated']:>20.12f} {r['reference']:>20.12f} {r['error_pct']:>11.4f}%")
        total_error += r['error_pct']

    avg_error = total_error / len(results)
    print("-" * 62)
    print(f"{'Average':<10} {'':<20} {'':<20} {avg_error:>11.4f}%")

    print("\n" + "="*80 + "\n")


if __name__ == "__main__":
    main()
