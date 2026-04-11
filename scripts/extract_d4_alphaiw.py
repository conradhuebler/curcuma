#!/usr/bin/env python3
"""
Extract alphaiw (frequency-dependent polarizabilities) from dftd4param.f90

Extracts the complete 23×7×118 alphaiw array for all elements and reference states.
This is the missing piece for accurate D4 C6 coefficient calculation.

Usage:
    python extract_d4_alphaiw.py

Output:
    test_cases/reference_data/d4_alphaiw_data.cpp

Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
Claude Generated - December 25, 2025
"""

import re
import sys
from pathlib import Path

def parse_fortran_data_line(line):
    """Extract numerical values from Fortran data statement"""
    # Remove comments
    line = re.sub(r'!.*$', '', line)

    # Find all floating point numbers
    numbers = re.findall(r'[-+]?\d+\.\d+(?:_wp)?', line)

    # Convert to Python floats
    return [float(n.replace('_wp', '')) for n in numbers]

def extract_alphaiw(fortran_file):
    """
    Extract alphaiw(23, 7, 118) array from Fortran dftd4param.f90

    Format in Fortran:
        data alphaiw(:,ref,elem) / value1, value2, ... /

    Returns:
        dict: {(elem, ref): [23 frequency values]}
    """
    alphaiw_data = {}

    with open(fortran_file, 'r') as f:
        content = f.read()

    # Pattern to match alphaiw data statements
    # Example: data alphaiw(:,1,  6) /
    pattern = r'data\s+alphaiw\s*\(\s*:\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*/([^/]+)/'

    matches = re.finditer(pattern, content, re.MULTILINE | re.DOTALL)

    for match in matches:
        ref = int(match.group(1))
        elem = int(match.group(2))
        data_block = match.group(3)

        # Parse all lines in the data block
        values = []
        for line in data_block.split('\n'):
            line = line.strip()
            if line.startswith('&'):
                line = line[1:].strip()  # Remove continuation character
            values.extend(parse_fortran_data_line(line))

        # Should have exactly 23 values (frequency points)
        if len(values) == 23:
            alphaiw_data[(elem, ref)] = values
            print(f"  ✓ Extracted alphaiw for element {elem:3d} ref {ref}: {len(values)} frequencies")
        elif len(values) > 0:
            print(f"  ⚠ Warning: element {elem} ref {ref} has {len(values)} values (expected 23)")
            alphaiw_data[(elem, ref)] = values

    return alphaiw_data

def generate_cpp_file(alphaiw_data, output_file):
    """Generate C++ source file with alphaiw data"""

    # Determine dimensions
    max_elem = max(elem for elem, ref in alphaiw_data.keys())
    max_refs_per_elem = {}
    for (elem, ref) in alphaiw_data.keys():
        max_refs_per_elem[elem] = max(max_refs_per_elem.get(elem, 0), ref)

    print(f"\nData dimensions:")
    print(f"  Elements: 1 to {max_elem}")
    print(f"  Max references per element: {max(max_refs_per_elem.values())}")
    print(f"  Frequency points: 23")
    print(f"  Total data points: {len(alphaiw_data)} reference states")

    with open(output_file, 'w') as f:
        f.write("""/*
 * D4 Frequency-Dependent Polarizabilities (alphaiw)
 * Extracted from external/gfnff/src/dftd4param.f90
 *
 * Array dimensions: alphaiw[element][reference][frequency]
 *   - 118 elements (H through Og)
 *   - Up to 7 reference states per element
 *   - 23 frequency points (imaginary frequency grid)
 *
 * Usage: C6_ij = (3/π) ∫ αᵢ(iω) × αⱼ(iω) dω
 *
 * Claude Generated - December 25, 2025
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include <vector>

// Initialize 3D array: d4_alphaiw[element-1][reference-1][frequency]
// Note: Use 0-based indexing (subtract 1 from Fortran element/reference numbers)

std::vector<std::vector<std::vector<double>>> d4_alphaiw_data;

void initialize_d4_alphaiw() {
    // Resize to [118 elements][7 references][23 frequencies]
    d4_alphaiw_data.resize(118, std::vector<std::vector<double>>(7, std::vector<double>(23, 0.0)));

""")

        # Write data grouped by element
        for elem in range(1, max_elem + 1):
            elem_refs = [(e, r) for (e, r) in alphaiw_data.keys() if e == elem]
            if not elem_refs:
                continue

            f.write(f"    // Element {elem} ({len(elem_refs)} reference states)\n")

            for (e, ref) in sorted(elem_refs):
                values = alphaiw_data[(e, ref)]
                elem_idx = e - 1  # Convert to 0-based
                ref_idx = ref - 1  # Convert to 0-based

                f.write(f"    d4_alphaiw_data[{elem_idx}][{ref_idx}] = {{\n")

                # Write in rows of 6 values for readability
                for i in range(0, len(values), 6):
                    chunk = values[i:i+6]
                    f.write("        " + ", ".join(f"{v:.7f}" for v in chunk))
                    if i + 6 < len(values):
                        f.write(",\n")
                    else:
                        f.write("\n")

                f.write("    };\n")

            f.write("\n")

        f.write("}\n")

    print(f"\n✓ Generated C++ file: {output_file}")
    print(f"  Total alphaiw reference states: {len(alphaiw_data)}")

def main():
    # Paths
    fortran_file = Path("external/gfnff/src/dftd4param.f90")
    output_file = Path("test_cases/reference_data/d4_alphaiw_data.cpp")

    if not fortran_file.exists():
        print(f"ERROR: Fortran file not found: {fortran_file}")
        print(f"  Current directory: {Path.cwd()}")
        print(f"  Please run from curcuma root directory")
        sys.exit(1)

    print("=" * 70)
    print("D4 alphaiw Data Extraction")
    print("=" * 70)
    print(f"Source: {fortran_file}")
    print(f"Output: {output_file}\n")

    # Extract data
    print("Extracting alphaiw data from Fortran...")
    alphaiw_data = extract_alphaiw(fortran_file)

    if not alphaiw_data:
        print("\nERROR: No alphaiw data found!")
        sys.exit(1)

    # Generate C++ file
    print("\nGenerating C++ source file...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    generate_cpp_file(alphaiw_data, output_file)

    print("\n" + "=" * 70)
    print("✓ Extraction completed successfully")
    print("=" * 70)
    print(f"\nNext steps:")
    print(f"  1. Include in d4param_generator.cpp: #include \"../../../{output_file}\"")
    print(f"  2. Call initialize_d4_alphaiw() in constructor")
    print(f"  3. Replace m_alpha_iw with d4_alphaiw_data")
    print(f"  4. Rebuild and test")

if __name__ == '__main__':
    main()
