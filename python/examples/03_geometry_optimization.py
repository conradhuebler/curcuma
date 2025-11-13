#!/usr/bin/env python3
"""
Example 3: Geometry Optimization
=================================

Demonstrates:
- Geometry optimization with different methods
- Convergence checking
- Trajectory output
- Energy comparison before/after optimization

Claude Generated Example
Copyright (C) 2019 - 2025 Conrad Hübler
"""

import curcuma

def main():
    print("=" * 60)
    print("Curcuma Python Example 3: Geometry Optimization")
    print("=" * 60)
    print()

    # Load a molecule
    print("1. Loading initial structure:")
    try:
        mol = curcuma.Molecule("../../test_cases/ethane.xyz")
        print(f"   Loaded: {mol}")
    except:
        mol = curcuma.Molecule(3, 0)
        print(f"   Created test molecule: {mol}")
    print()

    # Calculate initial energy
    print("2. Calculating initial energy:")
    initial_energy = curcuma.calculate_energy(mol, method="uff", verbosity=0)
    print(f"   Initial UFF energy: {initial_energy:.6f} Hartree")
    print()

    # Optimize with UFF
    print("3. Running UFF geometry optimization:")
    print("   Setting up optimizer...")

    opt_config = {
        "method": "uff",
        "max_iter": 100,
        "verbosity": 1,
        "write_traj": False,
    }

    optimizer = curcuma.Optimizer(opt_config)
    optimizer.set_molecule(mol)

    print("   Optimizing geometry...")
    optimizer.optimize()

    # Check convergence
    if optimizer.converged():
        print("   ✓ Optimization converged!")
    else:
        print("   ✗ Optimization did not converge")

    print(f"   Number of iterations: {optimizer.iterations()}")

    # Get optimized structure
    optimized_mol = optimizer.get_molecule()
    final_energy = optimizer.final_energy()
    print(f"   Final energy: {final_energy:.6f} Hartree")
    print(f"   Energy change: {(final_energy - initial_energy) * 2625.5:.2f} kJ/mol")
    print()

    # Save optimized structure
    print("4. Saving optimized structure:")
    try:
        optimized_mol.write_xyz("optimized_structure.xyz")
        print("   Saved to 'optimized_structure.xyz'")
    except Exception as e:
        print(f"   Note: Could not save file ({e})")
    print()

    # Use convenience function
    print("5. Using convenience function for optimization:")
    try:
        opt_mol_quick = curcuma.optimize_geometry(
            mol, method="uff", max_iter=50, verbosity=0
        )
        print(f"   Optimized: {opt_mol_quick}")
        quick_energy = curcuma.calculate_energy(opt_mol_quick, method="uff", verbosity=0)
        print(f"   Final energy: {quick_energy:.6f} Hartree")
    except Exception as e:
        print(f"   Note: Could not optimize ({e})")
    print()

    print("=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
