#!/usr/bin/env python3
"""
XTB GFN-FF Energy Component Extraction Script

Extracts all 7 energy components from XTB GFN-FF output and generates
C++ test reference data for test_gfnff_unified.cpp.

Usage:
    python extract_xtb_energy_components.py <molecule>.xyz [--xtb-path /path/to/xtb]

Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
Claude Generated - December 25, 2025
"""

import os
import sys
import subprocess
import argparse
import re
from pathlib import Path
from typing import Dict, Optional

class Color:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    END = '\033[0m'

class EnergyComponents:
    """Container for GFN-FF energy components"""
    def __init__(self):
        self.total_energy: Optional[float] = None
        self.bond_energy: Optional[float] = None
        self.angle_energy: Optional[float] = None
        self.torsion_energy: Optional[float] = None
        self.repulsion_energy: Optional[float] = None
        self.electrostat_energy: Optional[float] = None
        self.dispersion_energy: Optional[float] = None
        self.batm_energy: Optional[float] = None  # H-bond correction

    def is_complete(self) -> bool:
        """Check if all energy components have been extracted"""
        return all([
            self.total_energy is not None,
            self.bond_energy is not None,
            self.angle_energy is not None,
            self.torsion_energy is not None,
            self.repulsion_energy is not None,
            self.electrostat_energy is not None,
            self.dispersion_energy is not None,
            self.batm_energy is not None,
        ])

    def to_cpp_string(self, molecule_path: str) -> str:
        """Generate C++ test data string"""
        return (
            f'{{"molecules/larger/{Path(molecule_path).name}", '
            f'{self.total_energy:.12f}, '
            f'{self.bond_energy:.12f}, '
            f'{self.angle_energy:.12f}, '
            f'{self.torsion_energy:.12f}, '
            f'{self.repulsion_energy:.12f}, '
            f'{self.electrostat_energy:.12f}, '
            f'{self.dispersion_energy:.12f}, '
            f'{self.batm_energy:.12f}}}'
        )

def run_xtb(molecule_path: str, xtb_binary: str) -> str:
    """
    Run XTB with GFN-FF on the given molecule

    Args:
        molecule_path: Path to XYZ molecule file
        xtb_binary: Path to XTB executable

    Returns:
        XTB output as string

    Raises:
        RuntimeError: If XTB execution fails
    """
    cmd = [xtb_binary, molecule_path, "--gfnff"]

    print(f"{Color.CYAN}Running: {' '.join(cmd)}{Color.END}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"XTB execution failed:\n{e.stderr}")
    except subprocess.TimeoutExpired:
        raise RuntimeError("XTB execution timed out after 60 seconds")

def parse_xtb_output(output: str) -> EnergyComponents:
    """
    Parse XTB GFN-FF output to extract all energy components

    XTB GFN-FF output format (example):

    bond             :        12.345
    angle            :         1.234
    torsion          :         0.123
    repulsion        :         2.345
    electrostat      :        -1.234
    dispersion       :        -0.567
    batm             :        -0.012
    total energy     :        14.234 Eh

    Args:
        output: XTB stdout output

    Returns:
        EnergyComponents object with extracted values
    """
    components = EnergyComponents()

    # Parse energy breakdown (from GFN-FF energy decomposition)
    # XTB format: ":: bond energy               -2.508166283882 Eh    ::"
    patterns = {
        'bond': r'::\s*bond energy\s+([-+]?\d+\.\d+)\s+Eh',
        'angle': r'::\s*angle energy\s+([-+]?\d+\.\d+)\s+Eh',
        'torsion': r'::\s*torsion energy\s+([-+]?\d+\.\d+)\s+Eh',
        'repulsion': r'::\s*repulsion energy\s+([-+]?\d+\.\d+)\s+Eh',
        'electrostat': r'::\s*electrostat energy\s+([-+]?\d+\.\d+)\s+Eh',
        'dispersion': r'::\s*dispersion energy\s+([-+]?\d+\.\d+)\s+Eh',
        'batm': r'::\s*bonded atm energy\s+([-+]?\d+\.\d+)\s+Eh',
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, output, re.IGNORECASE)
        if match:
            value = float(match.group(1))
            setattr(components, f'{key}_energy', value)
            print(f"  {Color.GREEN}✓{Color.END} {key:12s}: {value:14.8f} Eh")
        else:
            print(f"  {Color.YELLOW}⚠{Color.END} {key:12s}: NOT FOUND")

    # Parse total energy (multiple possible formats)
    total_patterns = [
        r'::\s*total energy\s+([-+]?\d+\.\d+)\s+Eh',
        r'TOTAL ENERGY\s+([-+]?\d+\.\d+)\s+Eh',
        r'\|\s*TOTAL ENERGY\s+([-+]?\d+\.\d+)\s+Eh',
    ]

    for pattern in total_patterns:
        match = re.search(pattern, output, re.IGNORECASE)
        if match:
            components.total_energy = float(match.group(1))
            print(f"  {Color.BOLD}{Color.GREEN}✓{Color.END} {'total':12s}: {components.total_energy:14.8f} Eh{Color.END}")
            break
    else:
        print(f"  {Color.RED}✗{Color.END} {'total':12s}: NOT FOUND")

    return components

def validate_components(components: EnergyComponents) -> bool:
    """
    Validate that extracted energy components sum to total energy

    Args:
        components: Extracted energy components

    Returns:
        True if validation passes
    """
    if not components.is_complete():
        print(f"{Color.RED}ERROR: Not all energy components were extracted!{Color.END}")
        return False

    # Calculate sum of components
    component_sum = (
        components.bond_energy +
        components.angle_energy +
        components.torsion_energy +
        components.repulsion_energy +
        components.electrostat_energy +
        components.dispersion_energy +
        components.batm_energy
    )

    # Check if sum matches total (within tolerance)
    diff = abs(component_sum - components.total_energy)
    tolerance = 1e-6  # Hartree

    print(f"\n{Color.BOLD}Validation:{Color.END}")
    print(f"  Component sum: {component_sum:14.8f} Eh")
    print(f"  Total energy:  {components.total_energy:14.8f} Eh")
    print(f"  Difference:    {diff:14.8f} Eh")

    if diff < tolerance:
        print(f"  {Color.GREEN}✓ Components sum correctly (diff < {tolerance}){Color.END}")
        return True
    else:
        print(f"  {Color.RED}✗ Components do NOT sum correctly (diff = {diff}){Color.END}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Extract GFN-FF energy components from XTB output'
    )
    parser.add_argument(
        'molecule',
        help='Path to XYZ molecule file'
    )
    parser.add_argument(
        '--xtb-path',
        default='/home/conrad/Downloads/xtb-6.6.1/bin/xtb',
        help='Path to XTB executable (default: /home/conrad/Downloads/xtb-6.6.1/bin/xtb)'
    )
    parser.add_argument(
        '--output-cpp',
        action='store_true',
        help='Output C++ test data format'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.molecule):
        print(f"{Color.RED}ERROR: Molecule file not found: {args.molecule}{Color.END}")
        sys.exit(1)

    if not os.path.exists(args.xtb_path):
        print(f"{Color.RED}ERROR: XTB binary not found: {args.xtb_path}{Color.END}")
        sys.exit(1)

    print(f"\n{Color.BOLD}=== XTB Energy Component Extraction ==={Color.END}")
    print(f"Molecule: {args.molecule}")
    print(f"XTB:      {args.xtb_path}\n")

    try:
        # Run XTB
        output = run_xtb(args.molecule, args.xtb_path)

        # Parse energy components
        print(f"\n{Color.BOLD}Extracted Energy Components:{Color.END}")
        components = parse_xtb_output(output)

        # Validate
        if not validate_components(components):
            print(f"\n{Color.YELLOW}WARNING: Validation failed, but continuing...{Color.END}")

        # Output C++ test data
        if args.output_cpp or components.is_complete():
            print(f"\n{Color.BOLD}C++ Test Data:{Color.END}")
            print(f"{Color.CYAN}{components.to_cpp_string(args.molecule)},{Color.END}")
            print(f"\n{Color.GREEN}Copy this line to test_gfnff_unified.cpp{Color.END}")

        print(f"\n{Color.GREEN}✓ Extraction completed successfully{Color.END}")

    except Exception as e:
        print(f"\n{Color.RED}ERROR: {e}{Color.END}")
        sys.exit(1)

if __name__ == '__main__':
    main()
