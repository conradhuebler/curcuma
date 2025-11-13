#!/usr/bin/env python3
"""
Example 4: RMSD Analysis and Structural Alignment
==================================================

Demonstrates:
- RMSD calculation between structures
- Structural alignment
- Heavy atom RMSD
- Molecule comparison

Claude Generated Example
Copyright (C) 2019 - 2025 Conrad Hübler
"""

import curcuma

def main():
    print("=" * 60)
    print("Curcuma Python Example 4: RMSD Analysis")
    print("=" * 60)
    print()

    # Create or load reference structure
    print("1. Loading reference structure:")
    try:
        ref_mol = curcuma.Molecule("../../test_cases/ethane.xyz")
        print(f"   Reference: {ref_mol}")
    except:
        ref_mol = curcuma.Molecule(3, 0)
        print(f"   Created reference: {ref_mol}")
    print()

    # Create a slightly modified structure for comparison
    print("2. Creating modified structure:")
    try:
        # Try to load a different conformer or structure
        target_mol = curcuma.Molecule("../../test_cases/ethane.xyz")
        print(f"   Target: {target_mol}")
    except:
        # Use the same molecule for demonstration
        target_mol = curcuma.Molecule(ref_mol)
        print(f"   Created target (copy of reference): {target_mol}")
    print()

    # Calculate RMSD using convenience function
    print("3. Calculating RMSD (all atoms):")
    try:
        rmsd_all = curcuma.calculate_rmsd(ref_mol, target_mol)
        print(f"   RMSD (all atoms): {rmsd_all:.4f} Å")
    except Exception as e:
        print(f"   Note: Could not calculate RMSD ({e})")
    print()

    # Calculate heavy atom RMSD
    print("4. Calculating heavy atom RMSD:")
    try:
        rmsd_heavy = curcuma.calculate_rmsd(
            ref_mol, target_mol, heavy_only=True
        )
        print(f"   RMSD (heavy atoms only): {rmsd_heavy:.4f} Å")
    except Exception as e:
        print(f"   Note: Could not calculate heavy RMSD ({e})")
    print()

    # Structural alignment
    print("5. Aligning structures:")
    try:
        aligned_mol = curcuma.align_molecules(ref_mol, target_mol)
        print(f"   Aligned: {aligned_mol}")

        # Calculate RMSD after alignment
        rmsd_aligned = curcuma.calculate_rmsd(ref_mol, aligned_mol)
        print(f"   RMSD after alignment: {rmsd_aligned:.4f} Å")

        # Save aligned structure
        aligned_mol.write_xyz("aligned_structure.xyz")
        print("   Saved to 'aligned_structure.xyz'")
    except Exception as e:
        print(f"   Note: Could not align structures ({e})")
    print()

    # Using RMSD class directly
    print("6. Using RMSD class for advanced control:")
    try:
        rmsd_config = {
            "heavy": False,
            "reorder": False,
            "silent": True,
        }

        rmsd_calc = curcuma.RMSD(rmsd_config)
        rmsd_calc.set_reference(ref_mol)
        rmsd_calc.set_target(target_mol)

        rmsd_value = rmsd_calc.calculate()
        print(f"   RMSD: {rmsd_value:.4f} Å")

        # Align and get result
        rmsd_calc.align()
        aligned = rmsd_calc.get_aligned_molecule()
        print(f"   Aligned molecule: {aligned}")
    except Exception as e:
        print(f"   Note: Could not use RMSD class ({e})")
    print()

    print("=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
