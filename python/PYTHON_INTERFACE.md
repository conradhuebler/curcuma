# Curcuma Python Interface - Implementation Details

This document describes the Python interface implementation for Curcuma.

**Claude Generated: January 2025**

## Architecture

The Python interface uses **pybind11** to create C++ bindings for the Curcuma computational chemistry toolkit. The implementation follows an educational-first design philosophy, prioritizing clarity and ease of use.

### Directory Structure

```
curcuma/
├── python/
│   ├── pybind/                    # C++ binding source files
│   │   ├── curcuma_bindings.cpp   # Main module definition
│   │   ├── molecule_binding.cpp   # Molecule class bindings
│   │   ├── energy_calculator_binding.cpp
│   │   ├── rmsd_binding.cpp
│   │   ├── optimizer_binding.cpp
│   │   ├── simplemd_binding.cpp
│   │   └── confscan_binding.cpp
│   ├── curcuma/                   # Python package
│   │   └── __init__.py            # Package initialization
│   ├── examples/                  # Usage examples
│   │   ├── 01_basic_molecule.py
│   │   ├── 02_energy_calculation.py
│   │   ├── 03_geometry_optimization.py
│   │   └── 04_rmsd_analysis.py
│   ├── setup.py                   # Installation script
│   ├── pyproject.toml             # Modern Python packaging config
│   └── README.md                  # Python-specific documentation
```

## Components

### 1. Core Bindings (`curcuma_bindings.cpp`)

Main pybind11 module definition that ties together all submodules:
- Module documentation
- Version information
- Forward declarations for all binding functions
- Submodule registration

### 2. Molecule Bindings (`molecule_binding.cpp`)

Exposes the `Molecule` class with:
- Multiple constructors (empty, pre-allocated, from file, copy)
- File I/O (XYZ, MOL2, SDF, JSON)
- Geometry access and manipulation
- Property access (charge, mass, energy, fragments)
- Distance calculations and connectivity
- Python special methods (`__len__`, `__repr__`, `__str__`)

### 3. Energy Calculator Bindings (`energy_calculator_binding.cpp`)

Exposes the unified `EnergyCalculator` interface:
- Support for all QM/MM methods (UFF, EHT, GFN1/GFN2, etc.)
- Configuration via Python dictionaries (auto-converted to JSON)
- Energy and gradient calculations
- Property access (charges, bond orders, dipole)
- Convenience function: `calculate_energy()`

### 4. RMSD Bindings (`rmsd_binding.cpp`)

Structural analysis capabilities:
- RMSD calculation between structures
- Structural alignment
- Heavy atom filtering
- Atom reordering options
- Convenience functions: `calculate_rmsd()`, `align_molecules()`

### 5. Optimizer Bindings (`optimizer_binding.cpp`)

Geometry optimization interface:
- LBFGS and other optimization algorithms
- Convergence checking
- Trajectory output
- Convenience function: `optimize_geometry()`

### 6. SimpleMD Bindings (`simplemd_binding.cpp`)

Molecular dynamics simulations:
- NVE and NVT ensembles
- Velocity Verlet integration
- Temperature control
- Trajectory output
- Convenience function: `run_molecular_dynamics()`

### 7. ConfScan Bindings (`confscan_binding.cpp`)

Conformational analysis:
- Systematic dihedral scanning
- Energy-based filtering
- RMSD-based duplicate detection
- Ensemble output
- Convenience function: `scan_dihedral()`

## CMake Integration

The `CMakeLists.txt` includes a new section for Python bindings:

```cmake
option(BUILD_PYTHON_BINDINGS "Build Python bindings using pybind11" ON)

if(BUILD_PYTHON_BINDINGS)
    # Find Python and pybind11
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

    # Download pybind11 if not present
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v2.11.1
    )
    FetchContent_MakeAvailable(pybind11)

    # Create Python module
    pybind11_add_module(curcuma_py ${PYBIND_SRC})
    set_target_properties(curcuma_py PROPERTIES OUTPUT_NAME curcuma)
    target_link_libraries(curcuma_py PRIVATE curcuma_core)

    # Install to Python site-packages
    install(TARGETS curcuma_py LIBRARY DESTINATION ${Python3_SITEARCH})
endif()
```

## Design Principles

### 1. Educational-First API

All bindings prioritize clarity:
- Clear, descriptive names
- Comprehensive docstrings
- Intuitive parameter ordering
- Sensible default values

### 2. Pythonic Interface

Python special methods are implemented where appropriate:
- `__len__()` for atom count
- `__repr__()` for object representation
- `__str__()` for user-friendly output
- `__call__()` for callable objects

### 3. Convenience Functions

Top-level convenience functions provide simple interfaces for common tasks:
```python
# Instead of:
calc = EnergyCalculator("gfn2", config)
calc.set_molecule(mol)
energy = calc.calculate_energy()

# Users can write:
energy = calculate_energy(mol, "gfn2")
```

### 4. Configuration Handling

Python dictionaries are automatically converted to JSON for C++ consumption:
```python
config = {
    "method": "gfn2",
    "verbosity": 1,
    "max_iter": 100,
}
optimizer = Optimizer(config)
```

### 5. NumPy Integration

Geometry matrices are exposed as NumPy arrays via pybind11's Eigen integration:
```python
geometry = mol.geometry()  # Returns numpy array (3×N)
gradient = calc.gradient()  # Returns numpy array (3×N)
```

## Build and Installation

### Standard Installation

```bash
cd curcuma/python
pip install .
```

### Development Installation

```bash
cd curcuma/python
pip install -e .
```

### Build Options

Disable Python bindings in CMake:
```bash
cmake .. -DBUILD_PYTHON_BINDINGS=OFF
```

## Testing

Python tests should be added to `python/tests/`:
```python
import curcuma
import pytest

def test_molecule_creation():
    mol = curcuma.Molecule()
    assert len(mol) == 0

def test_energy_calculation():
    mol = curcuma.Molecule("test.xyz")
    energy = curcuma.calculate_energy(mol, "uff")
    assert isinstance(energy, float)
```

Run tests:
```bash
cd python
pytest tests/
```

## Future Enhancements

Potential additions to the Python interface:

1. **Trajectory Analysis**: Expose trajectory reading and analysis tools
2. **Hessian Analysis**: Expose frequency calculations and normal mode analysis
3. **QMDFF Fitting**: Expose force field parameter fitting
4. **Visualization**: Integration with ASE, RDKit, or PyMOL
5. **Type Hints**: Add complete type annotations for all functions
6. **Async Support**: Async versions of computationally expensive operations
7. **Context Managers**: Support `with` statements for resource management
8. **Serialization**: Pickle support for Molecule and other objects

## Copyright and Attribution

All Python bindings follow Curcuma's copyright policy:

```cpp
/*
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Python bindings for [component]
 */
```

Copyright remains with Conrad Hübler as project owner. Claude contributions are acknowledged in code comments.

## References

- [pybind11 Documentation](https://pybind11.readthedocs.io/)
- [Curcuma Main Repository](https://github.com/conradhuebler/curcuma)
- [Python Packaging Guide](https://packaging.python.org/)

---

**Implementation Date**: January 2025
**Implementation**: Claude AI Assistant
**Copyright**: Conrad Hübler (2019-2025)
