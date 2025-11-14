"""
Curcuma - Molecular Modelling and Simulation Toolkit
====================================================

A Python interface to the Curcuma computational chemistry toolkit,
providing access to:

- Quantum mechanical methods (Extended Hückel, GFN1/GFN2, etc.)
- Force field methods (UFF, GFN-FF, QMDFF)
- Geometry optimization (LBFGS and other algorithms)
- Molecular dynamics simulations
- Conformational analysis and searching
- RMSD calculations and structural alignment
- Dispersion corrections (D3, D4)

Example Usage
-------------
>>> import curcuma
>>> mol = curcuma.Molecule("structure.xyz")
>>> calc = curcuma.EnergyCalculator("gfn2")
>>> calc.set_molecule(mol)
>>> energy = calc.calculate_energy()
>>> print(f"Energy: {energy} Hartree")

For more information, see:
https://github.com/conradhuebler/curcuma

Claude Generated: Python package initialization
Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
"""

__version__ = "0.1.0"
__author__ = "Conrad Hübler"
__email__ = "Conrad.Huebler@gmx.net"
__license__ = "GPL-3.0"

# Import the compiled C++ extension module
try:
    from . import curcuma as _curcuma_core

    # Re-export fully implemented classes from the C++ module
    Molecule = _curcuma_core.Molecule
    EnergyCalculator = _curcuma_core.EnergyCalculator
    RMSD = _curcuma_core.RMSD

    # Fully implemented convenience functions
    calculate_energy = _curcuma_core.calculate_energy
    calculate_rmsd = _curcuma_core.calculate_rmsd
    align_molecules = _curcuma_core.align_molecules

    # Placeholder functions (not yet fully implemented - CLI-only)
    optimize_geometry = _curcuma_core.optimize_geometry
    run_molecular_dynamics = _curcuma_core.run_molecular_dynamics
    scan_dihedral = _curcuma_core.scan_dihedral

    # Module metadata
    if hasattr(_curcuma_core, "__version__"):
        __version__ = _curcuma_core.__version__

except ImportError as e:
    import warnings
    warnings.warn(
        f"Failed to import Curcuma C++ extension module: {e}\n"
        "Please ensure the module was compiled correctly.\n"
        "Try rebuilding with: pip install -v .",
        ImportWarning
    )
    raise

# Define what gets exported with "from curcuma import *"
__all__ = [
    # Fully implemented core classes
    "Molecule",
    "EnergyCalculator",
    "RMSD",

    # Fully implemented convenience functions
    "calculate_energy",
    "calculate_rmsd",
    "align_molecules",

    # Placeholder functions (CLI-only for now)
    "optimize_geometry",
    "run_molecular_dynamics",
    "scan_dihedral",

    # Metadata
    "__version__",
    "__author__",
    "__email__",
    "__license__",
]


def get_available_methods():
    """
    Get list of available computational methods

    Returns:
        dict: Dictionary of available methods by category

    Example:
        >>> methods = curcuma.get_available_methods()
        >>> print(methods["quantum"])
        ['eht', 'gfn1', 'gfn2', 'ipea1']
    """
    return {
        "force_fields": ["uff", "uff-d3", "qmdff", "cgfnff"],
        "quantum": ["eht", "gfn1", "gfn2", "ipea1"],
        "dispersion": ["d3", "d4"],
    }


def info():
    """
    Print Curcuma package information
    """
    print(f"Curcuma Python Interface v{__version__}")
    print(f"Author: {__author__}")
    print(f"License: {__license__}")
    print()
    print("Available computational methods:")
    methods = get_available_methods()
    for category, method_list in methods.items():
        print(f"  {category.replace('_', ' ').title()}:")
        for method in method_list:
            print(f"    - {method}")
    print()
    print("For documentation, see: https://github.com/conradhuebler/curcuma")
