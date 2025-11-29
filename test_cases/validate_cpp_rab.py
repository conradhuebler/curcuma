#!/usr/bin/env python3
"""
Validate C++ GFN-FF Parameter Generation Against Python Reference

This script runs curcuma on test molecules and extracts the RAB_TRANSFORM debug output,
comparing it with the Python reference implementation.
"""

import subprocess
import re
import math
from pathlib import Path

# Expected values from Python test_bond_energy_exact.py
EXPECTED_R0 = {
    "HH": 1.31280068,  # Bohr
    "OH": 1.52784633,  # Bohr
    "HCl": 2.30286684, # Bohr
}

TOLERANCE = 0.0001  # Bohr

def extract_rab_transform(output):
    """Extract RAB_TRANSFORM debug lines from curcuma output"""
    lines = []
    for line in output.split('\n'):
        if 'RAB_TRANSFORM' in line:
            lines.append(line)
    return lines

def parse_r_eq(output):
    """Extract final r_eq value from RAB_TRANSFORM output"""
    for line in output.split('\n'):
        if 'r_eq =' in line:
            # Extract the last number before "Bohr"
            match = re.search(r'r_eq = .* = ([0-9.]+) Bohr', line)
            if match:
                return float(match.group(1))
    return None

def run_test(mol_name, mol_path):
    """Run curcuma on molecule and extract r0"""
    print(f"\n{'='*60}")
    print(f"Testing: {mol_name} ({mol_path})")
    print('='*60)

    try:
        # Run curcuma - note: cgfnff may not be recognized, try alternative names
        cmd = ["./release/curcuma", mol_path, "-method", "uff"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        output = result.stdout + result.stderr

        # Extract RAB_TRANSFORM lines
        rab_lines = extract_rab_transform(output)

        if not rab_lines:
            print(f"  ✗ No RAB_TRANSFORM output found")
            print(f"  (This may mean the method was not called)")
            return False

        # Print extracted debug output
        print("\n  Debug Output:")
        for line in rab_lines:
            print(f"  {line}")

        # Extract final r_eq
        r_eq = parse_r_eq('\n'.join(rab_lines))

        if r_eq is None:
            print(f"  ✗ Could not parse r_eq from output")
            return False

        # Compare with expected
        expected = EXPECTED_R0[mol_name]
        error = abs(r_eq - expected)
        error_pct = 100.0 * error / expected

        print(f"\n  Calculated r0: {r_eq:.8f} Bohr")
        print(f"  Expected r0:   {expected:.8f} Bohr")
        print(f"  Error:         {error:.8f} Bohr ({error_pct:.4f}%)")

        if error < TOLERANCE:
            print(f"  ✓ PASS")
            return True
        else:
            print(f"  ✗ FAIL (error > {TOLERANCE})")
            return False

    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout")
        return False
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False

def main():
    test_molecules = [
        ("HH", "test_cases/molecules/dimers/HH.xyz"),
        ("OH", "test_cases/molecules/dimers/OH.xyz"),
        ("HCl", "test_cases/molecules/dimers/HCl.xyz"),
    ]

    print("\n" + "="*60)
    print("C++ GFN-FF Parameter Generation Validation")
    print("="*60)

    passed = 0
    failed = 0

    for mol_name, mol_path in test_molecules:
        if run_test(mol_name, mol_path):
            passed += 1
        else:
            failed += 1

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Passed: {passed}/3")
    print(f"Failed: {failed}/3")

    if failed == 0:
        print("\n✓ All tests PASSED")
        return 0
    else:
        print("\n✗ Some tests FAILED")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
