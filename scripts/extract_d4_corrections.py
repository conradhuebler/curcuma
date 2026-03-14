#!/usr/bin/env python3
"""
Extract D4 polarizability correction factors from dftd4param.f90

Extracts:
- ascale(ref, elem): Atomic scaling factors (7×118)
- sscale(refsys): Reference system scaling (17 values)
- secaiw(freq, refsys): Reference system polarizabilities (23×17)
- refsys(ref, elem): Reference system mapping (7×118)

These are needed for the correction:
    alpha_corrected = ascale * (alphaiw - hcount * sscale * secaiw)

Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
Claude Generated - December 25, 2025
"""

import re
import sys
from pathlib import Path

def parse_fortran_array_init(content, array_name, indices_pattern):
    """Parse Fortran data statements for an array"""
    data = {}

    # Pattern: data array_name (indices) / values /
    pattern = rf'data\s+{array_name}\s*\({indices_pattern}\)\s*/([^/]+)/'

    matches = re.finditer(pattern, content, re.MULTILINE | re.DOTALL)

    for match in matches:
        # Extract indices
        indices_str = match.group(1)
        # Extract values
        values_str = match.group(2)

        # Parse indices
        indices = [int(x.strip()) for x in indices_str.split(',')]

        # Parse values
        values = []
        for line in values_str.split('\n'):
            line = re.sub(r'!.*$', '', line).strip()  # Remove comments
            if line.startswith('&'):
                line = line[1:].strip()
            numbers = re.findall(r'[-+]?\d+\.\d+(?:_wp)?', line)
            values.extend([float(n.replace('_wp', '')) for n in numbers])

        # Store with indices as key
        data[tuple(indices)] = values

    return data

def extract_ascale(fortran_file):
    """Extract ascale(ref, elem) - atomic scaling factors"""
    with open(fortran_file, 'r') as f:
        content = f.read()

    # Pattern: data ascale  (ref,  elem) /  value_wp /
    pattern = r'data\s+ascale\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*([-+]?\d+\.\d+_wp)\s*/'

    ascale = {}
    for match in re.finditer(pattern, content):
        ref = int(match.group(1))
        elem = int(match.group(2))
        value = float(match.group(3).replace('_wp', ''))
        ascale[(elem, ref)] = value

    print(f"  Extracted {len(ascale)} ascale values")
    return ascale

def extract_sscale_and_secaiw(fortran_file):
    """Extract sscale and secaiw for reference systems"""
    with open(fortran_file, 'r') as f:
        content = f.read()

    sscale = {}
    secaiw = {}

    # Find all SEC system definitions
    # Pattern: SEC name  - - - - - - id
    #          data secq   (id) / value / ; data sscale(id)/ value /
    #          data secaiw (:,id) / ... /

    sec_pattern = r'! SEC\s+(\w+).*?(\d{4})\s*\n.*?data\s+sscale\s*\(\s*(\d+)\s*\)\s*/\s*([-+]?\d+\.\d+_wp)'

    for match in re.finditer(sec_pattern, content, re.MULTILINE | re.DOTALL):
        name = match.group(1)
        ref_id = int(match.group(3))
        sscale_val = float(match.group(4).replace('_wp', ''))
        sscale[ref_id] = sscale_val
        print(f"  sscale[{ref_id:2d}] = {sscale_val:.6f} (SEC {name})")

    # Extract secaiw(:,refsys)
    secaiw_pattern = r'data\s+secaiw\s*\(\s*:\s*,\s*(\d+)\s*\)\s*/([^/]+)/'

    for match in re.finditer(secaiw_pattern, content, re.MULTILINE | re.DOTALL):
        ref_id = int(match.group(1))
        values_str = match.group(2)

        values = []
        for line in values_str.split('\n'):
            line = re.sub(r'!.*$', '', line).strip()
            if line.startswith('&'):
                line = line[1:].strip()
            numbers = re.findall(r'[-+]?\d+\.\d+(?:_wp)?', line)
            values.extend([float(n.replace('_wp', '')) for n in numbers])

        if len(values) == 23:
            secaiw[ref_id] = values
            print(f"  secaiw[:, {ref_id:2d}] extracted (23 frequencies)")

    return sscale, secaiw

def extract_refsys(fortran_file):
    """Extract refsys(ref, elem) - reference system mapping"""
    with open(fortran_file, 'r') as f:
        content = f.read()

    # Pattern: data refsys  (ref,  elem) / id /
    pattern = r'data\s+refsys\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*(\d+)\s*/'

    refsys = {}
    for match in re.finditer(pattern, content):
        ref = int(match.group(1))
        elem = int(match.group(2))
        sys_id = int(match.group(3))
        refsys[(elem, ref)] = sys_id

    print(f"  Extracted {len(refsys)} refsys mappings")
    return refsys

def generate_cpp_file(ascale, sscale, secaiw, refsys, output_file):
    """Generate C++ source file with correction factors"""

    # Determine dimensions
    max_elem = max(elem for elem, ref in ascale.keys())
    max_refsys = max(sscale.keys())

    with open(output_file, 'w') as f:
        f.write("""/*
 * D4 Polarizability Correction Factors
 * Extracted from external/gfnff/src/dftd4param.f90
 *
 * Correction formula: α_corrected = ascale * (αᵢⱼw - hcount * sscale * secaiw)
 *
 * Arrays:
 *   - ascale[elem][ref]: Atomic scaling factors (118×7)
 *   - sscale[refsys]: Reference system scaling (17 values)
 *   - secaiw[refsys][freq]: Reference polarizabilities (17×23)
 *   - refsys[elem][ref]: Reference system mapping (118×7)
 *
 * Claude Generated - December 25, 2025
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include <vector>
#include <map>

// Atomic scaling factors: ascale[element-1][reference-1]
std::vector<std::vector<double>> d4_ascale_data;

// Reference system scaling: sscale[refsys_id]
std::map<int, double> d4_sscale_data;

// Reference system polarizabilities: secaiw[refsys_id][frequency]
std::map<int, std::vector<double>> d4_secaiw_data;

// Reference system mapping: refsys[element-1][reference-1] -> refsys_id
std::vector<std::vector<int>> d4_refsys_data;

void initialize_d4_corrections() {
    // Resize ascale and refsys to [118 elements][7 references]
    d4_ascale_data.resize(118, std::vector<double>(7, 1.0));
    d4_refsys_data.resize(118, std::vector<int>(7, 0));

""")

        # Write ascale data
        f.write("    // Atomic scaling factors (ascale)\n")
        for (elem, ref), value in sorted(ascale.items()):
            elem_idx = elem - 1
            ref_idx = ref - 1
            f.write(f"    d4_ascale_data[{elem_idx}][{ref_idx}] = {value:.14f};\n")

        f.write("\n    // Reference system mapping (refsys)\n")
        for (elem, ref), sys_id in sorted(refsys.items()):
            elem_idx = elem - 1
            ref_idx = ref - 1
            f.write(f"    d4_refsys_data[{elem_idx}][{ref_idx}] = {sys_id};\n")

        f.write("\n    // Reference system scaling (sscale)\n")
        for sys_id, value in sorted(sscale.items()):
            f.write(f"    d4_sscale_data[{sys_id}] = {value:.14f};\n")

        f.write("\n    // Reference system polarizabilities (secaiw)\n")
        for sys_id, values in sorted(secaiw.items()):
            f.write(f"    d4_secaiw_data[{sys_id}] = {{\n")
            for i in range(0, len(values), 6):
                chunk = values[i:i+6]
                f.write("        " + ", ".join(f"{v:.7f}" for v in chunk))
                if i + 6 < len(values):
                    f.write(",\n")
                else:
                    f.write("\n")
            f.write("    };\n")

        f.write("}\n")

    print(f"\n✓ Generated C++ file: {output_file}")
    print(f"  ascale entries: {len(ascale)}")
    print(f"  sscale entries: {len(sscale)}")
    print(f"  secaiw systems: {len(secaiw)}")
    print(f"  refsys mappings: {len(refsys)}")

def main():
    fortran_file = Path("external/gfnff/src/dftd4param.f90")
    output_file = Path("test_cases/reference_data/d4_corrections_data.cpp")

    if not fortran_file.exists():
        print(f"ERROR: Fortran file not found: {fortran_file}")
        sys.exit(1)

    print("=" * 70)
    print("D4 Correction Factors Extraction")
    print("=" * 70)
    print(f"Source: {fortran_file}")
    print(f"Output: {output_file}\n")

    print("Extracting ascale...")
    ascale = extract_ascale(fortran_file)

    print("\nExtracting sscale and secaiw...")
    sscale, secaiw = extract_sscale_and_secaiw(fortran_file)

    print("\nExtracting refsys...")
    refsys = extract_refsys(fortran_file)

    print("\nGenerating C++ source file...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    generate_cpp_file(ascale, sscale, secaiw, refsys, output_file)

    print("\n" + "=" * 70)
    print("✓ Extraction completed successfully")
    print("=" * 70)

if __name__ == '__main__':
    main()
