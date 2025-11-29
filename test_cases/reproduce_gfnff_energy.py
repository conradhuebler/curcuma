#!/usr/bin/env python3
"""
GFN-FF Bond Energy Reproduction Script

Reproduces GFN-FF bond stretching energy from:
- XYZ geometry file (Ångström)
- JSON parameters (vbond from XTB output)

Formula: E_bond = k * (r - r0)^2
  r   = actual distance [Bohr]
  r0  = vbond(1) equilibrium offset [Bohr]
  k   = vbond(3) force constant prefactor [Hartree/Bohr²]

Source: Fortran gfnff_engrad.F90 lines 373-399

Claude Generated - Parameter extraction and energy reproduction
"""

import json
import math
import sys
from pathlib import Path


class GFNFFEnergyCalculator:
    """Calculate GFN-FF bond stretching energy from XYZ geometry and JSON parameters"""

    # Unit conversion constants
    ANGSTROM_TO_BOHR = 1.889726  # 1 Å = 1.889726 Bohr
    HARTREE_TO_EV = 27.211386   # 1 Hartree = 27.211386 eV
    HARTREE_TO_KCAL_MOL = 627.509469  # 1 Hartree = 627.509469 kcal/mol

    def __init__(self, xyz_file, json_file):
        """
        Initialize calculator from geometry and parameter files

        Args:
            xyz_file: Path to XYZ coordinate file (Ångström)
            json_file: Path to JSON parameter file (from XTB)
        """
        self.atoms = []  # List of (element, x, y, z) in Ångström
        self.coordinates_bohr = []  # Coordinates converted to Bohr
        self.bonds = []  # Bond list [(atom1, atom2), ...]
        self.vbond_params = []  # vbond parameters from JSON

        self.read_xyz(xyz_file)
        self.read_vbond_json(json_file)

        # Convert coordinates to Bohr for calculation
        self._to_bohr()

    def read_xyz(self, filename):
        """Parse XYZ coordinate file (Ångström)"""
        with open(filename, 'r') as f:
            lines = f.readlines()

        natoms = int(lines[0].strip())
        # lines[1] is comment

        for i in range(2, 2 + natoms):
            parts = lines[i].split()
            element = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            self.atoms.append((element, x, y, z))

    def read_vbond_json(self, filename):
        """Extract vbond parameters from JSON file"""
        with open(filename, 'r') as f:
            data = json.load(f)

        # Extract bond list and parameters
        if "blist" in data:
            # blist format: [[atom1, atom2], ...] with 1-based indexing
            for bond_pair in data["blist"]:
                # Convert from 1-based to 0-based indexing
                self.bonds.append((bond_pair[0] - 1, bond_pair[1] - 1))

        if "vbond" in data:
            # vbond format: [[r0_shift, k_steepness, k_strength], ...]
            self.vbond_params = data["vbond"]

    def _to_bohr(self):
        """Convert atomic coordinates from Ångström to Bohr"""
        self.coordinates_bohr = []
        for element, x, y, z in self.atoms:
            x_bohr = x * self.ANGSTROM_TO_BOHR
            y_bohr = y * self.ANGSTROM_TO_BOHR
            z_bohr = z * self.ANGSTROM_TO_BOHR
            self.coordinates_bohr.append((x_bohr, y_bohr, z_bohr))

    def _distance(self, i, j):
        """Calculate distance between atoms i and j in Bohr"""
        x1, y1, z1 = self.coordinates_bohr[i]
        x2, y2, z2 = self.coordinates_bohr[j]
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def calculate_bond_energy(self):
        """
        Calculate total bond stretching energy

        E_bond(i) = k * (r - r0)^2

        Returns:
            tuple: (total_energy_hartree, list of per-bond energies)
        """
        total_energy = 0.0
        bond_energies = []

        for bond_idx, (i, j) in enumerate(self.bonds):
            # Get vbond parameters for this bond
            if bond_idx >= len(self.vbond_params):
                raise ValueError(f"Bond {bond_idx} exceeds vbond parameter count")

            r0, k_steepness, k_strength = self.vbond_params[bond_idx]

            # Calculate actual distance
            r_actual = self._distance(i, j)

            # Calculate energy: E = k * (r - r0)^2
            dr = r_actual - r0
            E_bond = k_strength * dr * dr

            bond_energies.append({
                'bond': (i, j),
                'atoms': (self.atoms[i][0], self.atoms[j][0]),
                'r_actual': r_actual,
                'r0': r0,
                'k_strength': k_strength,
                'dr': dr,
                'energy_hartree': E_bond,
                'energy_ev': E_bond * self.HARTREE_TO_EV,
                'energy_kcal_mol': E_bond * self.HARTREE_TO_KCAL_MOL
            })

            total_energy += E_bond

        return total_energy, bond_energies

    def print_results(self, total_energy, bond_energies):
        """Print energy calculation results"""
        print(f"\n{'='*70}")
        print(f"GFN-FF Bond Energy Reproduction")
        print(f"{'='*70}")

        # Per-bond breakdown
        print(f"\nPer-Bond Energy Breakdown:")
        print(f"{'Bond':<10} {'Atoms':<15} {'r_actual':>12} {'r0':>12} {'dr':>12} {'E [Eh]':>15}")
        print(f"{'-'*86}")

        for e in bond_energies:
            atoms_str = f"{e['atoms'][0]}-{e['atoms'][1]}"
            print(f"({e['bond'][0]},{e['bond'][1]:<7})"
                  f" {atoms_str:<15}"
                  f" {e['r_actual']:>12.8f}"
                  f" {e['r0']:>12.8f}"
                  f" {e['dr']:>12.8f}"
                  f" {e['energy_hartree']:>15.10f}")

        print(f"\n{'='*70}")
        print(f"Total Bond Energy:")
        print(f"  {total_energy:20.10f} Hartree")
        print(f"  {total_energy * self.HARTREE_TO_EV:20.10f} eV")
        print(f"  {total_energy * self.HARTREE_TO_KCAL_MOL:20.10f} kcal/mol")
        print(f"{'='*70}\n")


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 reproduce_gfnff_energy.py <xyz_file> <json_file>")
        print("\nExample:")
        print("  python3 reproduce_gfnff_energy.py molecules/dimers/HH.xyz molecules/dimers/HH.json")
        sys.exit(1)

    xyz_file = Path(sys.argv[1])
    json_file = Path(sys.argv[2])

    if not xyz_file.exists():
        print(f"Error: XYZ file not found: {xyz_file}")
        sys.exit(1)

    if not json_file.exists():
        print(f"Error: JSON file not found: {json_file}")
        sys.exit(1)

    try:
        calc = GFNFFEnergyCalculator(str(xyz_file), str(json_file))
        total_energy, bond_energies = calc.calculate_bond_energy()
        calc.print_results(total_energy, bond_energies)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
