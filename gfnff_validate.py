#!/usr/bin/env python3
"""
GFN-FF Validation Framework
Systematically compares Curcuma native GFN-FF against Fortran reference
"""

import subprocess
import re
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple

@dataclass
class EnergyBreakdown:
    """Energy component breakdown"""
    total: float = 0.0
    bond: float = 0.0
    angle: float = 0.0
    torsion: float = 0.0
    repulsion: float = 0.0
    electrostatic: float = 0.0
    dispersion: float = 0.0
    coulomb: float = 0.0  # Curcuma's name

    def dict(self):
        return {
            "total": self.total,
            "bond": self.bond,
            "angle": self.angle,
            "torsion": self.torsion,
            "repulsion": self.repulsion,
            "electrostatic": self.electrostatic,
            "dispersion": self.dispersion,
        }

    def __repr__(self):
        return json.dumps(self.dict(), indent=2)

class GFNFFValidator:
    """Validates GFN-FF implementation"""

    def __init__(self, curcuma_exe: str, reference_exe: str, test_dir: str):
        self.curcuma_exe = Path(curcuma_exe)
        self.reference_exe = Path(reference_exe)
        self.test_dir = Path(test_dir)

    def run_reference(self, xyz_file: str) -> EnergyBreakdown:
        """Run Fortran reference implementation"""
        try:
            result = subprocess.run(
                [str(self.reference_exe), str(self.test_dir / xyz_file)],
                cwd=str(self.test_dir),
                capture_output=True,
                text=True,
                timeout=30
            )

            breakdown = EnergyBreakdown()
            output = result.stdout + result.stderr

            # Extract energy values from output
            total_match = re.search(r'Total energy:\s+([-\d.]+)', output)
            if total_match:
                breakdown.total = float(total_match.group(1))

            bond_match = re.search(r'Bond energy:\s+([-\d.]+)\s+Eh', output)
            if bond_match:
                breakdown.bond = float(bond_match.group(1))

            angle_match = re.search(r'Angle energy:\s+([-\d.]+)\s+Eh', output)
            if angle_match:
                breakdown.angle = float(angle_match.group(1))

            torsion_match = re.search(r'Torsion energy:\s+([-\d.]+)\s+Eh', output)
            if torsion_match:
                breakdown.torsion = float(torsion_match.group(1))

            repulsion_match = re.search(r'Repulsion energy:\s+([-\d.]+)\s+Eh', output)
            if repulsion_match:
                breakdown.repulsion = float(repulsion_match.group(1))

            electrostatic_match = re.search(r'Electrostatic[:\s]+([-\d.]+)\s+Eh', output)
            if electrostatic_match:
                breakdown.electrostatic = float(electrostatic_match.group(1))

            dispersion_match = re.search(r'Dispersion[:\s]+([-\d.]+)\s+Eh', output)
            if dispersion_match:
                breakdown.dispersion = float(dispersion_match.group(1))

            return breakdown

        except Exception as e:
            print(f"Error running reference: {e}")
            return None

    def run_curcuma(self, xyz_file: str) -> EnergyBreakdown:
        """Run Curcuma native GFN-FF"""
        try:
            result = subprocess.run(
                [str(self.curcuma_exe), "-sp", str(self.test_dir / xyz_file),
                 "-method", "cgfnff", "-verbosity", "2"],
                capture_output=True,
                text=True,
                timeout=30
            )

            # Remove ANSI color codes
            output = result.stdout + result.stderr
            output = re.sub(r'\x1b\[[0-9;]*m', '', output)

            breakdown = EnergyBreakdown()

            # Extract total energy
            total_match = re.search(r'Single Point Energy = ([-\d.]+)\s+Eh', output)
            if total_match:
                breakdown.total = float(total_match.group(1))

            # Extract component energies
            # Pattern: [PARAM] bond_energy: -0.165008 Eh
            bond_match = re.search(r'bond_energy:\s+([-\d.]+)\s+Eh', output)
            if bond_match:
                breakdown.bond = float(bond_match.group(1))

            angle_match = re.search(r'angle_energy:\s+([-\d.]+)\s+Eh', output)
            if angle_match:
                breakdown.angle = float(angle_match.group(1))

            # Look for repulsion - try multiple patterns
            repulsion_match = re.search(r'thread_repulsion_energy:\s+([-\d.]+)\s+Eh', output)
            if repulsion_match:
                breakdown.repulsion = float(repulsion_match.group(1))
            else:
                repulsion_match = re.search(r'\w\w_repulsion:\s+([-\d.]+)\s+Eh', output)
                if repulsion_match:
                    breakdown.repulsion = float(repulsion_match.group(1))

            # Look for dispersion
            dispersion_match = re.search(r'thread_dispersion_energy:\s+([-\d.]+)\s+Eh', output)
            if dispersion_match:
                breakdown.dispersion = float(dispersion_match.group(1))
            else:
                dispersion_match = re.search(r'GFNFF_dispersion:\s+([-\d.]+)\s+Eh', output)
                if dispersion_match:
                    breakdown.dispersion = float(dispersion_match.group(1))

            # Look for Coulomb - try multiple patterns
            coulomb_match = re.search(r'thread_coulomb_energy:\s+([-\d.]+)\s+Eh', output)
            if coulomb_match:
                breakdown.coulomb = float(coulomb_match.group(1))
            else:
                # Try alternative pattern for Coulomb output
                coulomb_match = re.search(r'GFNFF_coulomb:\s+([-\d.]+)\s+Eh', output)
                if coulomb_match:
                    breakdown.coulomb = float(coulomb_match.group(1))

            return breakdown

        except Exception as e:
            print(f"Error running Curcuma: {e}")
            return None

    def compare_energies(self, ref: EnergyBreakdown, cur: EnergyBreakdown) -> Dict:
        """Compare energy components"""
        def relative_error(ref, cur):
            if abs(ref) < 1e-10:
                return 0.0 if abs(cur) < 1e-10 else float('inf')
            return abs((cur - ref) / ref) * 100.0

        comparison = {
            "total": {
                "reference": ref.total,
                "curcuma": cur.total,
                "error_abs": cur.total - ref.total,
                "error_pct": relative_error(ref.total, cur.total)
            },
            "bond": {
                "reference": ref.bond,
                "curcuma": cur.bond,
                "error_abs": cur.bond - ref.bond,
                "error_pct": relative_error(ref.bond, cur.bond)
            },
            "angle": {
                "reference": ref.angle,
                "curcuma": cur.angle,
                "error_abs": cur.angle - ref.angle,
                "error_pct": relative_error(ref.angle, cur.angle)
            },
            "repulsion": {
                "reference": ref.repulsion,
                "curcuma": cur.repulsion,
                "error_abs": cur.repulsion - ref.repulsion,
                "error_pct": relative_error(ref.repulsion, cur.repulsion)
            },
            "dispersion": {
                "reference": ref.dispersion,
                "curcuma": cur.dispersion,
                "error_abs": cur.dispersion - ref.dispersion,
                "error_pct": relative_error(ref.dispersion, cur.dispersion)
            },
            "electrostatic": {
                "reference": ref.electrostatic,
                "curcuma": cur.coulomb,
                "error_abs": cur.coulomb - ref.electrostatic,
                "error_pct": relative_error(ref.electrostatic, cur.coulomb)
            },
        }

        return comparison

    def validate_molecule(self, xyz_file: str) -> Dict:
        """Validate a single molecule"""
        print(f"\n{'='*80}")
        print(f"Validating: {xyz_file}")
        print(f"{'='*80}")

        # Run both implementations
        ref_energy = self.run_reference(xyz_file)
        cur_energy = self.run_curcuma(xyz_file)

        if not ref_energy or not cur_energy:
            print("ERROR: Failed to run implementations")
            return None

        # Compare
        comparison = self.compare_energies(ref_energy, cur_energy)

        # Print results
        self._print_results(xyz_file, comparison)

        return comparison

    def _print_results(self, molecule: str, comparison: Dict):
        """Pretty print comparison results"""
        print(f"\nMolecule: {molecule}")
        print(f"\n{'Term':<15} {'Reference':>15} {'Curcuma':>15} {'Abs Error':>15} {'%  Error':>10}")
        print("-" * 75)

        for term in ["bond", "angle", "repulsion", "dispersion", "electrostatic", "total"]:
            if term in comparison:
                c = comparison[term]
                pct = c["error_pct"]
                if pct == float('inf'):
                    pct_str = "∞"
                else:
                    pct_str = f"{pct:.2f}%"

                print(f"{term:<15} {c['reference']:>15.10f} {c['curcuma']:>15.10f} "
                      f"{c['error_abs']:>15.10f} {pct_str:>10}")

                # Color code: Green if <1%, Yellow if <5%, Red if >5%
                if pct < 1.0:
                    status = "✓ EXCELLENT"
                elif pct < 5.0:
                    status = "⚠ GOOD"
                elif pct < 10.0:
                    status = "! MODERATE"
                else:
                    status = "✗ POOR"
                print(f"         Status: {status}")

    def validate_all(self, molecules: List[str]) -> Dict:
        """Validate all molecules"""
        results = {}
        for mol in molecules:
            results[mol] = self.validate_molecule(mol)
        return results

if __name__ == "__main__":
    # Paths
    curcuma_exe = "/home/conrad/src/claude_curcuma/curcuma/release/curcuma"
    reference_exe = "/home/conrad/src/claude_curcuma/curcuma/external/gfnff/build/test/gfnff-gfnff_analyze-test"
    test_dir = "/home/conrad/src/claude_curcuma/curcuma/test_cases/molecules/dimers"

    # Test molecules
    molecules = ["HH.xyz", "HCl.xyz", "OH.xyz"]

    validator = GFNFFValidator(curcuma_exe, reference_exe, test_dir)
    results = validator.validate_all(molecules)

    # Summary
    print(f"\n{'='*80}")
    print("VALIDATION SUMMARY")
    print(f"{'='*80}")

    # Write results to file
    with open("gfnff_validation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\nResults saved to: gfnff_validation_results.json")
