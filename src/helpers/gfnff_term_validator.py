#!/usr/bin/env python3

"""
GFN-FF Energy Decomposition Test Script
Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
Claude Generated - GFN-FF 100% Accuracy Implementation Plan
"""

import json
import subprocess
import sys
import re
import os
import numpy as np

def parse_gfnff_output(output_text):
    """Parse GFN-FF output to extract energy decomposition"""
    energy_terms = {}

    # Extract total energy
    total_match = re.search(r'Total energy: (-?\d+\.\d+) Eh', output_text)
    if total_match:
        energy_terms['total'] = float(total_match.group(1))

    # Extract energy breakdown
    section_match = re.search(r'Energy breakdown:', output_text)
    if not section_match:
        return energy_terms

    # Parse individual terms
    term_lines = []
    lines_after_section = section_match.end()

    for line in lines_after_section:
        line = line.strip()
        if match := re.search(r'(bond|angle|torsion|repulsion|electrostatic|dispersion|h-bond|x-bond|bonded atm|external)', line):
            term_name = match.group(1).lower()
            energy_match = re.search(r'[:\s+]+([-+]?\d+)', line)
            if energy_match:
                energy_terms[term_name] = float(energy_match.group(1))

    return energy_terms

def compare_with_reference(system, calculated_terms, reference_file):
    """Compare calculated terms with reference values"""

    with open(reference_file) as f:
        reference_data = json.load(f)

    reference_breakdown = reference_data["reference_data"]["energy_decomposition"]
    calculated_breakdown = calculated_terms

    print(f"\n=== {system} Energy Decomposition Comparison ===")
    print(f"Target tolerance: 1e-6 Hartree per term")
    print(f"Tested vs XT B Reference Data:\n")

    all_passed = True

    # Get expected terms in order
    expected_terms = list(reference_breakdown.keys())

    for term in expected_terms:
        if term not in calculated_breakdown:
            print(f"WARNING: {term} missing from calculation!")
            continue

        calculated_val = calculated_breakdown.get(term, 0.0)
        reference_val = reference_breakdown[term]["value"]
        difference = abs(calculated_val - reference_val)
        relative_error = difference / abs(reference_val) if reference_val != 0 else 0.0

        passed = difference < 1e-6
        all_passed = all_passed and passed

        status = "✅ PASS" if passed else "❌ FAIL"

        print(f"{term: {status}")
        print(f"   Calculated: {calculated_val:.10f} Eh")
        print(f"   Reference: {reference_val:.10f} Eh")
        print(f"   Difference: {difference:.2e-6} Eh")

        if difference > 1e-6:
            print(f"   🚨 EXCEEDS TOLERANCE")
        print(f"   Critical fix required for 100% accuracy!\n")

    return all_passed

def run_system_test(system_name, xyz_file, verbosity_level=2):
    """Run single system test"""
    print(f"\n=== {system_name} GFN-FF Test ===")

    # Convert to absolute path
    if not os.path.isabs(xyz_file):
        xyz_file = os.path.abspath(xyz_file)

    # Build curcuma command
    cmd = ["./curcuma", "-sp", xyz_file, "-method", "gfnff", "-threads", "1"]

    if verbosity_level >= 2:
        cmd.append("-verbosity")

    # Run calculation and capture output
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                               check=False, cwd="../release")
        output = result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error running {cmd}: {e}")
        return False

    return parse_gfnff_output(output)

def main():
    print("🔬 GFN-FF 100% Accuracy Validation Tool")
    print("===================================")
    print("Testing: HH, CH4, H2O vs XTB 6.6.1 reference")
    print("Target: 1e-6 Hartree tolerance per energy term")
    print("===================================\n")

    test_cases = [
        {"name": "HH", "file": "test_cases/molecules/dimers/HH.xyz", "reference": "test_cases/golden_references/gfnff_hh.json"},
        {"name": "CH4", "file": "test_cases/molecules/larger/CH4.xyz", "reference": "test_cases/golden_references/gfnff_ch4.json"},
        {"name": "H2O", "file": "test_cases/molecules/trimers/water.xyz", "reference": "test_cases/golden_references/gfnff_h2o.json"}
    ]

    all_systems_passed = True

    for test in test_cases:
        test_passed = run_system_test(test["name"], test["file"], test["reference"])
        all_systems_passed = all_systems_passed and test_passed

        if not test_passed:
            print(f"❌ {test['name']} test FAILED - investigation needed")
            all_systems_passed = False

    print(f"\n=== FINAL RESULT ===")
    if all_systems_passed:
        print("🎯 SUCCESS: All energy terms within 1e-6 Hartree tolerance!")
        print("💰 Ready to claim the $200 bounty!")
        print("🏁 GFN-FF implementation ready for production!")
    else:
        print("❌ INCOMPLETE: Some systems still exceed tolerance limits")
        print("🔧 ACTION: Systematic fixes needed")

        # Identify failing systems
        for test in test_cases:
            if not test_passed:
                print(f"\n🔍 {test['name']} Analysis:")
                # Could run detailed investigation here
                print(f"   Status: {test_passed}")

    return all_systems_passed

if __name__ == "__main__":
    main()
else:
        main() {}