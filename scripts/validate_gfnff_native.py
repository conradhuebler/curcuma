#!/usr/bin/env python3
"""
GFN-FF Native vs External Validation Script

Compares native C++ implementation (cgfnff) against external Fortran
library (gfnff) for comprehensive validation.

Usage:
    python validate_gfnff_native.py [--tolerance 0.5] [--test-dir path]
"""

import os
import sys
import subprocess
import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple
import re

class Color:
    """ANSI color codes"""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

def parse_energy(output: str) -> float:
    """Extract energy from curcuma output"""
    # Look for patterns like "Energy: -123.456 Hartree"
    patterns = [
        r'Energy:\s*([-+]?\d+\.\d+)',
        r'Total Energy:\s*([-+]?\d+\.\d+)',
        r'E\s*=\s*([-+]?\d+\.\d+)',
    ]

    for pattern in patterns:
        match = re.search(pattern, output)
        if match:
            return float(match.group(1))

    raise ValueError(f"Could not parse energy from output:\n{output}")

def parse_gradient(output: str) -> List[List[float]]:
    """Extract gradient from curcuma output"""
    # This is a simplified parser - adjust based on actual output format
    gradients = []
    in_gradient_section = False

    for line in output.split('\n'):
        if 'Gradient' in line or 'gradient' in line:
            in_gradient_section = True
            continue

        if in_gradient_section:
            # Look for lines with 3 floats
            match = re.findall(r'[-+]?\d+\.\d+', line)
            if len(match) == 3:
                gradients.append([float(x) for x in match])
            elif gradients and len(match) == 0:
                # End of gradient section
                break

    return gradients

def run_curcuma(xyz_file: Path, method: str, curcuma_bin: str = "./curcuma") -> Tuple[float, List[List[float]], str]:
    """Run curcuma and extract results"""
    cmd = [curcuma_bin, "-sp", str(xyz_file), "-method", method, "-grad", "true"]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
            check=False
        )

        if result.returncode != 0:
            raise RuntimeError(f"Curcuma failed: {result.stderr}")

        energy = parse_energy(result.stdout)
        gradient = parse_gradient(result.stdout)

        return energy, gradient, result.stdout

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Curcuma timeout for {xyz_file}")
    except Exception as e:
        raise RuntimeError(f"Error running curcuma: {e}")

def compare_energies(e_native: float, e_external: float, tolerance_kcal: float = 0.5) -> Tuple[bool, float]:
    """Compare energies with tolerance in kcal/mol"""
    HARTREE_TO_KCAL = 627.5094740631
    diff_kcal = abs(e_native - e_external) * HARTREE_TO_KCAL
    passed = diff_kcal <= tolerance_kcal
    return passed, diff_kcal

def compare_gradients(g_native: List[List[float]], g_external: List[List[float]],
                      tolerance_percent: float = 5.0) -> Tuple[bool, float]:
    """Compare gradients with relative tolerance"""
    if len(g_native) != len(g_external):
        return False, 999.9

    max_diff_percent = 0.0
    for gn, ge in zip(g_native, g_external):
        for gn_i, ge_i in zip(gn, ge):
            if abs(ge_i) > 1e-10:  # Avoid division by zero
                diff_percent = abs(gn_i - ge_i) / abs(ge_i) * 100
                max_diff_percent = max(max_diff_percent, diff_percent)

    passed = max_diff_percent <= tolerance_percent
    return passed, max_diff_percent

def validate_molecule(xyz_file: Path, tolerance_kcal: float, curcuma_bin: str, verbose: bool = False) -> Dict:
    """Validate single molecule"""
    mol_name = xyz_file.stem
    print(f"\n{Color.BOLD}Testing: {mol_name}{Color.END}")

    result = {
        "molecule": mol_name,
        "energy_passed": False,
        "gradient_passed": False,
        "energy_diff_kcal": 0.0,
        "gradient_diff_percent": 0.0,
        "error": None
    }

    try:
        # Run external (reference)
        print(f"  Running external GFN-FF... ", end='', flush=True)
        e_ext, g_ext, out_ext = run_curcuma(xyz_file, "gfnff", curcuma_bin)
        print(f"{Color.GREEN}✓{Color.END} E = {e_ext:.6f} Hartree")

        # Run native
        print(f"  Running native cgfnff... ", end='', flush=True)
        e_nat, g_nat, out_nat = run_curcuma(xyz_file, "cgfnff", curcuma_bin)
        print(f"{Color.GREEN}✓{Color.END} E = {e_nat:.6f} Hartree")

        # Compare energies
        energy_passed, energy_diff = compare_energies(e_nat, e_ext, tolerance_kcal)
        result["energy_passed"] = energy_passed
        result["energy_diff_kcal"] = energy_diff

        status = f"{Color.GREEN}PASS{Color.END}" if energy_passed else f"{Color.RED}FAIL{Color.END}"
        print(f"  Energy difference: {energy_diff:.3f} kcal/mol [{status}]")

        # Compare gradients (if available)
        if g_nat and g_ext:
            gradient_passed, gradient_diff = compare_gradients(g_nat, g_ext)
            result["gradient_passed"] = gradient_passed
            result["gradient_diff_percent"] = gradient_diff

            status = f"{Color.GREEN}PASS{Color.END}" if gradient_passed else f"{Color.RED}FAIL{Color.END}"
            print(f"  Gradient difference: {gradient_diff:.2f}% [{status}]")

        if verbose:
            print(f"\n{Color.BLUE}External Output:{Color.END}")
            print(out_ext[:500])
            print(f"\n{Color.BLUE}Native Output:{Color.END}")
            print(out_nat[:500])

    except Exception as e:
        result["error"] = str(e)
        print(f"  {Color.RED}ERROR: {e}{Color.END}")

    return result

def generate_report(results: List[Dict], tolerance_kcal: float) -> str:
    """Generate validation report"""
    total = len(results)
    energy_passed = sum(1 for r in results if r["energy_passed"])
    gradient_passed = sum(1 for r in results if r.get("gradient_passed", False))
    errors = sum(1 for r in results if r["error"] is not None)

    report = f"""
{'='*80}
GFN-FF NATIVE VALIDATION REPORT
{'='*80}

Configuration:
  Tolerance: {tolerance_kcal} kcal/mol (energy)
  Tolerance: 5.0% (gradient)

Summary:
  Total molecules: {total}
  Energy tests passed: {energy_passed}/{total} ({energy_passed/total*100:.1f}%)
  Gradient tests passed: {gradient_passed}/{total} ({gradient_passed/total*100:.1f}%)
  Errors: {errors}

Detailed Results:
"""

    # Sort by energy difference
    results_sorted = sorted(results, key=lambda r: r.get("energy_diff_kcal", 999.9))

    report += f"\n{'Molecule':<30} {'ΔE (kcal/mol)':<15} {'ΔGrad (%)':<12} {'Status':<10}\n"
    report += f"{'-'*80}\n"

    for r in results_sorted:
        mol = r["molecule"]
        e_diff = r.get("energy_diff_kcal", 0.0)
        g_diff = r.get("gradient_diff_percent", 0.0)

        if r["error"]:
            status = "ERROR"
        elif r["energy_passed"] and r.get("gradient_passed", True):
            status = "✅ PASS"
        elif r["energy_passed"]:
            status = "⚠️  PARTIAL"
        else:
            status = "❌ FAIL"

        report += f"{mol:<30} {e_diff:<15.3f} {g_diff:<12.2f} {status:<10}\n"

        if r["error"]:
            report += f"  Error: {r['error']}\n"

    report += f"\n{'='*80}\n"

    # Overall assessment
    pass_rate = energy_passed / total * 100
    if pass_rate >= 95:
        assessment = f"{Color.GREEN}EXCELLENT{Color.END} - Production ready"
    elif pass_rate >= 80:
        assessment = f"{Color.YELLOW}GOOD{Color.END} - Minor issues to address"
    elif pass_rate >= 50:
        assessment = f"{Color.YELLOW}FAIR{Color.END} - Significant work needed"
    else:
        assessment = f"{Color.RED}POOR{Color.END} - Major implementation issues"

    report += f"\nOverall Assessment: {assessment}\n"
    report += f"{'='*80}\n"

    return report

def main():
    parser = argparse.ArgumentParser(description="Validate native GFN-FF implementation")
    parser.add_argument("--tolerance", type=float, default=0.5,
                        help="Energy tolerance in kcal/mol (default: 0.5)")
    parser.add_argument("--test-dir", type=str, default="test_cases/gfnff_validation",
                        help="Directory containing test XYZ files")
    parser.add_argument("--curcuma", type=str, default="./build/curcuma",
                        help="Path to curcuma binary")
    parser.add_argument("--verbose", action="store_true",
                        help="Print detailed output")
    parser.add_argument("--output", type=str, default="gfnff_validation_report.txt",
                        help="Output report file")

    args = parser.parse_args()

    # Find test molecules
    test_dir = Path(args.test_dir)
    if not test_dir.exists():
        print(f"{Color.RED}Error: Test directory not found: {test_dir}{Color.END}")
        sys.exit(1)

    xyz_files = list(test_dir.glob("**/*.xyz"))
    if not xyz_files:
        print(f"{Color.RED}Error: No XYZ files found in {test_dir}{Color.END}")
        sys.exit(1)

    print(f"{Color.BOLD}GFN-FF Native Validation{Color.END}")
    print(f"Found {len(xyz_files)} test molecules")
    print(f"Tolerance: {args.tolerance} kcal/mol")

    # Validate each molecule
    results = []
    for xyz_file in xyz_files:
        result = validate_molecule(xyz_file, args.tolerance, args.curcuma, args.verbose)
        results.append(result)

    # Generate report
    report = generate_report(results, args.tolerance)
    print(report)

    # Save to file
    with open(args.output, 'w') as f:
        # Remove ANSI colors for file output
        clean_report = re.sub(r'\x1b\[[0-9;]+m', '', report)
        f.write(clean_report)

    print(f"\nReport saved to: {args.output}")

    # Exit with appropriate code
    total_passed = sum(1 for r in results if r["energy_passed"])
    if total_passed == len(results):
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
