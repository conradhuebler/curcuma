#!/usr/bin/env python3
"""
Example 2: Energy Calculations
===============================

Demonstrates:
- Creating energy calculators with different methods
- Single-point energy calculations
- Gradient calculations
- Accessing atomic charges and bond orders

Claude Generated Example
Copyright (C) 2019 - 2025 Conrad Hübler
"""

import curcuma

def main():
    print("=" * 60)
    print("Curcuma Python Example 2: Energy Calculations")
    print("=" * 60)
    print()

    # Load a molecule
    print("1. Loading molecule:")
    try:
        mol = curcuma.Molecule("../../test_cases/ethane.xyz")
        print(f"   Loaded: {mol}")
    except:
        # Create a simple molecule if file not found
        mol = curcuma.Molecule(3, 0)
        print(f"   Created test molecule: {mol}")
    print()

    # Universal Force Field (UFF) calculation
    print("2. UFF Energy Calculation:")
    print("   Creating UFF calculator...")
    uff_calc = curcuma.EnergyCalculator("uff", {"verbosity": 1})
    uff_calc.set_molecule(mol)

    print("   Calculating energy...")
    uff_energy = uff_calc.calculate_energy()
    print(f"   UFF Energy: {uff_energy:.6f} Hartree")
    print(f"              {uff_energy * 2625.5:.2f} kJ/mol")
    print()

    # Extended Hückel Theory (EHT) calculation
    print("3. Extended Hückel Theory (EHT) Calculation:")
    try:
        print("   Creating EHT calculator...")
        eht_calc = curcuma.EnergyCalculator("eht", {"verbosity": 1})
        eht_calc.set_molecule(mol)

        print("   Calculating energy and gradient...")
        eht_energy = eht_calc.calculate_energy(gradient=True)
        print(f"   EHT Energy: {eht_energy:.6f} Hartree")

        # Get gradient
        gradient = eht_calc.gradient()
        print(f"   Gradient shape: {gradient.shape}")
        print(f"   Max gradient component: {abs(gradient).max():.6f}")
        print()
    except Exception as e:
        print(f"   Note: EHT not available ({e})")
        print()

    # GFN2-xTB calculation (if available)
    print("4. GFN2-xTB Calculation (if available):")
    try:
        print("   Creating GFN2 calculator...")
        gfn2_calc = curcuma.EnergyCalculator("gfn2", {"verbosity": 1})
        gfn2_calc.set_molecule(mol)

        print("   Calculating energy...")
        gfn2_energy = gfn2_calc.calculate_energy()
        print(f"   GFN2 Energy: {gfn2_energy:.6f} Hartree")

        # Get atomic charges
        charges = gfn2_calc.charges()
        print(f"   Atomic charges (first 3 atoms): {charges[:3]}")
        print()
    except Exception as e:
        print(f"   Note: GFN2 not available ({e})")
        print()

    # Convenience function for quick calculations
    print("5. Using convenience function:")
    try:
        energy = curcuma.calculate_energy(mol, method="uff", verbosity=0)
        print(f"   Quick UFF energy: {energy:.6f} Hartree")
    except Exception as e:
        print(f"   Note: Could not calculate energy ({e})")
    print()

    print("=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
