# Curcuma Python Interface - User Guide

**Claude Generated: January 2025**

## Overview

The Curcuma Python interface provides a complete, Pythonic API for computational chemistry calculations, making all capabilities of Curcuma accessible from Python with a simple and intuitive interface.

## Installation

### Requirements

- Python >= 3.7
- CMake >= 3.18
- C++14 or C++17 compiler
- NumPy >= 1.19

### Installation Steps

```bash
# Navigate to curcuma repository
cd curcuma

# Build and install Python bindings
cd python
pip install .

# Or for development (editable install)
pip install -e .
```

### Verify Installation

```python
import curcuma
curcuma.info()
```

## Quick Start Guide

### 1. Working with Molecules

```python
import curcuma

# Load molecule from file
mol = curcuma.Molecule("structure.xyz")

# Access properties
print(f"Number of atoms: {len(mol)}")
print(f"Molecular charge: {mol.charge()}")
print(f"Molecular mass: {mol.mass():.2f} amu")

# Get geometry as NumPy array
geometry = mol.geometry()  # Shape: (3, N)
print(f"Geometry shape: {geometry.shape}")

# Calculate distances
dist = mol.calculate_distance(0, 1)
print(f"Distance between atoms 0 and 1: {dist:.3f} Å")

# Save to file
mol.write_xyz("output.xyz")
```

### 2. Energy Calculations

```python
import curcuma

# Load molecule
mol = curcuma.Molecule("molecule.xyz")

# Quick energy calculation (convenience function)
energy_uff = curcuma.calculate_energy(mol, method="uff")
print(f"UFF Energy: {energy_uff:.6f} Hartree")

# Using EnergyCalculator for more control
calc = curcuma.EnergyCalculator("gfn2", {"verbosity": 1})
calc.set_molecule(mol)

# Calculate energy and gradient
energy = calc.calculate_energy(gradient=True)
gradient = calc.gradient()  # NumPy array (3, N)

# Get molecular properties
charges = calc.charges()
dipole = calc.dipole_moment()

print(f"Energy: {energy:.6f} Hartree")
print(f"Max gradient: {abs(gradient).max():.6f}")
print(f"Atomic charges: {charges}")
```

### 3. Geometry Optimization

```python
import curcuma

# Load initial structure
mol = curcuma.Molecule("initial.xyz")

# Quick optimization (convenience function)
optimized = curcuma.optimize_geometry(
    mol,
    method="gfn2",
    max_iter=100,
    verbosity=1
)
optimized.write_xyz("optimized.xyz")

# Using Optimizer for more control
opt_config = {
    "method": "uff",
    "max_iter": 100,
    "grad_conv": 1e-4,
    "energy_conv": 1e-6,
    "verbosity": 2,
    "write_traj": True,
}

optimizer = curcuma.Optimizer(opt_config)
optimizer.set_molecule(mol)
optimizer.optimize()

if optimizer.converged():
    print("✓ Optimization converged!")
    print(f"Iterations: {optimizer.iterations()}")
    print(f"Final energy: {optimizer.final_energy():.6f} Hartree")

    result = optimizer.get_molecule()
    result.write_xyz("optimized_structure.xyz")
```

### 4. RMSD Analysis

```python
import curcuma

# Load structures
reference = curcuma.Molecule("reference.xyz")
target = curcuma.Molecule("structure.xyz")

# Quick RMSD calculation
rmsd = curcuma.calculate_rmsd(reference, target)
print(f"RMSD: {rmsd:.3f} Å")

# Heavy atom RMSD only
rmsd_heavy = curcuma.calculate_rmsd(
    reference, target, heavy_only=True
)
print(f"Heavy atom RMSD: {rmsd_heavy:.3f} Å")

# Align structures
aligned = curcuma.align_molecules(reference, target)
aligned.write_xyz("aligned.xyz")

# Using RMSD class for advanced control
rmsd_calc = curcuma.RMSD({
    "heavy": True,
    "reorder": False,
    "silent": True
})

rmsd_calc.set_reference(reference)
rmsd_calc.set_target(target)
rmsd_value = rmsd_calc.calculate()
rmsd_calc.align()
aligned_mol = rmsd_calc.get_aligned_molecule()
```

### 5. Molecular Dynamics

```python
import curcuma

# Load molecule
mol = curcuma.Molecule("molecule.xyz")

# Quick MD simulation
final_structure = curcuma.run_molecular_dynamics(
    mol,
    method="uff",
    temperature=300,  # Kelvin
    steps=1000,
    timestep=0.5,  # femtoseconds
    verbosity=1
)
final_structure.write_xyz("md_final.xyz")

# Using MolecularDynamics for more control
md_config = {
    "method": "gfn2",
    "timestep": 0.5,
    "temperature": 350,
    "steps": 5000,
    "ensemble": "nvt",
    "dump": 10,  # Write every 10 steps
    "thermostat": "berendsen",
    "verbosity": 2,
}

md = curcuma.MolecularDynamics(md_config)
md.set_molecule(mol)
md.run()

final = md.get_molecule()
final.write_xyz("md_result.xyz")

# Write trajectory
md.write_trajectory("trajectory.xyz")
```

### 6. Conformational Analysis

```python
import curcuma

# Load molecule
mol = curcuma.Molecule("flexible_molecule.xyz")

# Quick dihedral scan
conformers = curcuma.scan_dihedral(
    mol,
    dihedral_atoms=[0, 1, 2, 3],  # Atom indices (0-based)
    method="gfn2",
    steps=36,  # 10-degree increments (360°/36 = 10°)
    verbosity=1
)

print(f"Generated {len(conformers)} conformers")

# Save conformer ensemble
for i, conformer in enumerate(conformers):
    conformer.write_xyz(f"conformer_{i:03d}.xyz")

# Using ConformationalScan for more control
scan_config = {
    "method": "uff",
    "dihedral": [1, 5, 8, 12],  # Atom indices
    "steps": 72,  # 5-degree increments
    "optimize": True,
    "rmsd_threshold": 0.5,  # Å
    "energy_window": 50,  # kJ/mol
    "verbosity": 2,
}

scanner = curcuma.ConformationalScan(scan_config)
scanner.set_molecule(mol)
scanner.run()

conformers = scanner.get_conformers()
scanner.write_ensemble("conformers.xyz")
```

## Available Methods

### Computational Methods

```python
# Get available methods
methods = curcuma.get_available_methods()

print("Force Fields:", methods["force_fields"])
# ['uff', 'uff-d3', 'qmdff', 'cgfnff']

print("Quantum Methods:", methods["quantum"])
# ['eht', 'gfn1', 'gfn2', 'ipea1']

print("Dispersion Corrections:", methods["dispersion"])
# ['d3', 'd4']
```

## Configuration Options

All classes accept configuration dictionaries:

```python
config = {
    "method": "gfn2",          # Calculation method
    "verbosity": 1,            # Output level (0-3)
    "threads": 4,              # Number of threads
    "max_iter": 100,           # Maximum iterations
    "grad_conv": 1e-4,         # Gradient convergence
    "energy_conv": 1e-6,       # Energy convergence
    "write_traj": True,        # Write trajectory
}
```

### Verbosity Levels

- **0**: Silent (no output)
- **1**: Minimal (final results only)
- **2**: Normal (scientific analysis)
- **3**: Debug (complete details)

## Integration with Other Tools

### NumPy Integration

All geometry and gradient arrays are NumPy arrays:

```python
import numpy as np
import curcuma

mol = curcuma.Molecule("molecule.xyz")

# Get geometry as NumPy array
geom = mol.geometry()  # Shape: (3, N)
print(f"Type: {type(geom)}")  # <class 'numpy.ndarray'>

# Manipulate with NumPy
geom_copy = geom.copy()
geom_copy += np.random.randn(*geom.shape) * 0.1  # Add noise

# Set modified geometry
mol.set_geometry(geom_copy)
```

### ASE (Atomic Simulation Environment)

```python
import curcuma
import ase
import ase.io

# Convert Curcuma Molecule to ASE Atoms
mol = curcuma.Molecule("molecule.xyz")
geom = mol.geometry()  # (3, N)
elements = [mol.get_element(i) for i in range(len(mol))]

# Create ASE Atoms object
atoms = ase.Atoms(
    numbers=elements,
    positions=geom.T  # ASE uses (N, 3) shape
)

# Use ASE visualization
ase.io.write("structure.png", atoms)
```

## Error Handling

```python
import curcuma

try:
    mol = curcuma.Molecule("nonexistent.xyz")
except Exception as e:
    print(f"Error loading molecule: {e}")

try:
    calc = curcuma.EnergyCalculator("invalid_method")
except Exception as e:
    print(f"Invalid method: {e}")
```

## Performance Tips

1. **Use silent mode for batch calculations**:
   ```python
   config = {"verbosity": 0}
   ```

2. **Reuse calculators**:
   ```python
   calc = curcuma.EnergyCalculator("gfn2")
   for mol in molecules:
       calc.set_molecule(mol)
       energy = calc.calculate_energy()
   ```

3. **Enable multi-threading**:
   ```python
   config = {"threads": 4}
   ```

## Examples

Complete working examples are available in `python/examples/`:

- `01_basic_molecule.py` - Molecule operations
- `02_energy_calculation.py` - Energy calculations
- `03_geometry_optimization.py` - Optimization workflows
- `04_rmsd_analysis.py` - RMSD and alignment

## Getting Help

### In Python

```python
import curcuma

# Module info
curcuma.info()

# Help on specific classes
help(curcuma.Molecule)
help(curcuma.EnergyCalculator)
help(curcuma.optimize_geometry)
```

### Documentation

- Main documentation: [CLAUDE.md](../CLAUDE.md)
- Implementation details: [python/PYTHON_INTERFACE.md](../python/PYTHON_INTERFACE.md)
- Python package README: [python/README.md](../python/README.md)

## Citation

If you use Curcuma in your research, please cite:

```bibtex
@software{curcuma,
  author = {Hübler, Conrad},
  title = {Curcuma: Molecular Modelling and Simulation Toolkit},
  year = {2025},
  url = {https://github.com/conradhuebler/curcuma},
  doi = {10.5281/zenodo.4302722}
}
```

---

**Created**: January 2025
**Author**: Conrad Hübler
**Implementation**: Claude AI Assistant
**License**: GPL-3.0
