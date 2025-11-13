# Curcuma Python Interface

Python bindings for the Curcuma computational chemistry toolkit, providing a simple and intuitive interface for molecular modeling, quantum chemistry calculations, and conformational analysis.

## Features

- **Molecular Structures**: Load, manipulate, and analyze molecular geometries
- **Energy Calculations**: Support for multiple QM and MM methods (UFF, GFN1/GFN2, EHT, etc.)
- **Geometry Optimization**: LBFGS and other optimization algorithms
- **Molecular Dynamics**: Simple MD simulations with various ensembles
- **Conformational Analysis**: Systematic conformational scanning and searching
- **RMSD Analysis**: Structure comparison and alignment
- **Dispersion Corrections**: D3 and D4 corrections

## Installation

### From Source (Development)

```bash
# Clone the repository
git clone https://github.com/conradhuebler/curcuma.git
cd curcuma

# Build and install Python bindings
cd python
pip install -v .

# Or for development (editable install)
pip install -e .
```

### Requirements

- Python >= 3.7
- CMake >= 3.18
- C++14 or C++17 compiler
- pybind11 >= 2.11 (automatically downloaded)
- NumPy >= 1.19

## Quick Start

### Basic Molecule Operations

```python
import curcuma

# Load a molecule from file
mol = curcuma.Molecule("structure.xyz")

# Access properties
print(f"Atoms: {len(mol)}")
print(f"Charge: {mol.charge()}")
print(f"Mass: {mol.mass():.2f} amu")

# Calculate distances
dist = mol.calculate_distance(0, 1)
print(f"Distance: {dist:.3f} Å")

# Save to file
mol.write_xyz("output.xyz")
```

### Energy Calculations

```python
import curcuma

# Load molecule
mol = curcuma.Molecule("molecule.xyz")

# UFF energy calculation
energy = curcuma.calculate_energy(mol, method="uff")
print(f"Energy: {energy:.6f} Hartree")

# GFN2-xTB calculation (if available)
calc = curcuma.EnergyCalculator("gfn2", {"verbosity": 1})
calc.set_molecule(mol)
energy = calc.calculate_energy(gradient=True)
gradient = calc.gradient()
charges = calc.charges()
```

### Geometry Optimization

```python
import curcuma

# Load initial structure
mol = curcuma.Molecule("initial.xyz")

# Optimize geometry
optimized = curcuma.optimize_geometry(
    mol,
    method="gfn2",
    max_iter=100,
    verbosity=1
)

# Save optimized structure
optimized.write_xyz("optimized.xyz")
```

### RMSD Analysis

```python
import curcuma

# Load structures
ref = curcuma.Molecule("reference.xyz")
target = curcuma.Molecule("structure.xyz")

# Calculate RMSD
rmsd = curcuma.calculate_rmsd(ref, target, heavy_only=True)
print(f"RMSD: {rmsd:.3f} Å")

# Align structures
aligned = curcuma.align_molecules(ref, target)
aligned.write_xyz("aligned.xyz")
```

### Molecular Dynamics

```python
import curcuma

# Load molecule
mol = curcuma.Molecule("molecule.xyz")

# Run MD simulation
final_structure = curcuma.run_molecular_dynamics(
    mol,
    method="uff",
    temperature=300,  # Kelvin
    steps=1000,
    timestep=0.5  # femtoseconds
)

final_structure.write_xyz("md_final.xyz")
```

### Conformational Analysis

```python
import curcuma

# Load molecule
mol = curcuma.Molecule("flexible.xyz")

# Scan dihedral angle
conformers = curcuma.scan_dihedral(
    mol,
    dihedral_atoms=[0, 1, 2, 3],  # atom indices
    method="gfn2",
    steps=36  # 10-degree increments
)

print(f"Generated {len(conformers)} conformers")
```

## Available Methods

### Force Fields
- `uff`: Universal Force Field
- `uff-d3`: UFF with D3 dispersion correction
- `qmdff`: Quantum Mechanically Derived Force Field
- `cgfnff`: Native GFN-FF implementation (work in progress)

### Quantum Methods
- `eht`: Extended Hückel Theory
- `gfn1`: GFN1-xTB (requires TBLite, XTB, or Ulysses)
- `gfn2`: GFN2-xTB (requires TBLite, Ulysses, or XTB)
- `ipea1`: iPEA1-xTB (requires TBLite)

### Dispersion Corrections
- `d3`: DFT-D3 correction
- `d4`: DFT-D4 correction

## Examples

Comprehensive examples are provided in the `examples/` directory:

- `01_basic_molecule.py` - Molecule creation and manipulation
- `02_energy_calculation.py` - Energy calculations with different methods
- `03_geometry_optimization.py` - Geometry optimization workflows
- `04_rmsd_analysis.py` - RMSD and structural alignment

Run examples:
```bash
cd python/examples
python 01_basic_molecule.py
```

## Documentation

For detailed documentation of the C++ backend, see the main [Curcuma documentation](https://github.com/conradhuebler/curcuma).

### API Reference

All Python classes and functions include detailed docstrings:

```python
import curcuma
help(curcuma.Molecule)
help(curcuma.EnergyCalculator)
help(curcuma.optimize_geometry)
```

## Development

### Building from Source

```bash
# Clone repository
git clone https://github.com/conradhuebler/curcuma.git
cd curcuma

# Build C++ library and Python bindings
mkdir build && cd build
cmake .. -DBUILD_PYTHON_BINDINGS=ON
make -j4

# Python module will be built in the build directory
```

### Running Tests

```bash
cd python
pytest tests/
```

## License

Curcuma is licensed under the GNU General Public License v3.0.

Copyright (C) 2019 - 2025 Conrad Hübler

## Citation

If you use Curcuma in your research, please cite:

```
@software{curcuma,
  author = {Hübler, Conrad},
  title = {Curcuma: Molecular Modelling and Simulation Toolkit},
  year = {2025},
  url = {https://github.com/conradhuebler/curcuma},
  doi = {10.5281/zenodo.4302722}
}
```

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues on GitHub.

## Contact

Conrad Hübler - Conrad.Huebler@gmx.net

Project Link: https://github.com/conradhuebler/curcuma

---

*Claude Generated: Python interface for Curcuma computational chemistry toolkit*
