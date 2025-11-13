#!/usr/bin/env python3
"""
Example 1: Basic Molecule Operations
=====================================

Demonstrates:
- Loading molecules from files
- Accessing molecular properties
- Basic geometry manipulation
- File I/O operations

Claude Generated Example
Copyright (C) 2019 - 2025 Conrad Hübler
"""

import curcuma

def main():
    print("=" * 60)
    print("Curcuma Python Example 1: Basic Molecule Operations")
    print("=" * 60)
    print()

    # Create an empty molecule
    print("1. Creating an empty molecule:")
    mol = curcuma.Molecule()
    print(f"   {mol}")
    print()

    # Load a molecule from file (assuming test_cases directory exists)
    print("2. Loading molecule from file:")
    try:
        # You can replace this with your own XYZ file
        mol = curcuma.Molecule("../../test_cases/ethane.xyz")
        print(f"   Loaded: {mol}")
        print(f"   Number of atoms: {len(mol)}")
        print(f"   Molecular charge: {mol.charge()}")
        print(f"   Molecular mass: {mol.mass():.2f} amu")
        print()
    except Exception as e:
        print(f"   Note: Could not load test file ({e})")
        print("   Creating water molecule manually instead...")

        # Create a simple water molecule manually
        mol = curcuma.Molecule(3, 0)  # 3 atoms, charge 0
        print(f"   Created: {mol}")
        print()

    # Access geometry
    print("3. Accessing geometry:")
    geometry = mol.geometry()
    print(f"   Geometry shape: {geometry.shape}")
    print(f"   Geometry (first 3 atoms):")
    for i in range(min(3, len(mol))):
        atom = mol.get_atom(i)
        print(f"      Atom {i}: ({atom[0]:.3f}, {atom[1]:.3f}, {atom[2]:.3f}) Å")
    print()

    # Calculate distances
    if len(mol) >= 2:
        print("4. Distance calculations:")
        dist = mol.calculate_distance(0, 1)
        print(f"   Distance between atoms 0 and 1: {dist:.3f} Å")
        print()

    # Center of mass
    print("5. Center of mass:")
    com = mol.center_of_mass()
    print(f"   Center of mass: ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}) Å")
    print()

    # Save molecule
    print("6. Saving molecule:")
    try:
        mol.write_xyz("example_output.xyz")
        print("   Saved to 'example_output.xyz'")
    except Exception as e:
        print(f"   Note: Could not save file ({e})")
    print()

    # JSON export
    print("7. JSON export:")
    json_data = mol.to_json()
    print(f"   Exported to JSON (keys: {list(json_data.keys())})")
    print()

    print("=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
