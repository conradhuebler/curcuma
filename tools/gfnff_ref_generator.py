#!/usr/bin/env python3
import subprocess
import json
import re
import sys
import os
from datetime import datetime

def parse_analyzer_output(output, xyz_file):
    ref_data = {
        "molecule": {
            "name": os.path.basename(xyz_file).replace(".xyz", ""),
            "geometry_file": xyz_file,
            "atoms": []
        },
        "metadata": {
            "generator": "gfnff-analyzer-parser",
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        },
        "topology": {
            "coordination_numbers": [],
            "hybridizations": [],
            "charges": []
        },
        "bonds": [],
        "angles": [],
        "torsions": [],
        "energy_components": {},
        "gradient_decomposition": {}
    }

    lines = output.splitlines()

    # Float regex pattern: optionally signed, at least one digit before or after decimal point, optional exponent
    FLOAT_PATTERN = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"

    # 1. Parse Atoms and Topology
    atom_section = False
    for line in lines:
        if "atom   neighbors  erfCN metchar sp-hybrid imet pi  qest" in line:
            atom_section = True
            continue
        if atom_section:
            if not line.strip():
                if len(ref_data["molecule"]["atoms"]) > 0:
                    atom_section = False
                continue
            # Handle the coordinates table
            match = re.match(r"^\s*(\d+)\s+([A-Z][a-z]?)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-]?[\d.]+)\s+([-]?[\d.]+)\s+([-]?[\d.]+)\s+([-]?[\d.]+)", line)
            if match:
                idx, element, neighbors, erfCN, metchar, hybrid, imet, pi, qest, x, y, z = match.groups()
                ref_data["molecule"]["atoms"].append({
                    "element": element,
                    "x": float(x),
                    "y": float(y),
                    "z": float(z)
                })
                ref_data["topology"]["coordination_numbers"].append(float(erfCN))
                ref_data["topology"]["hybridizations"].append(int(hybrid))
                ref_data["topology"]["charges"].append(float(qest))
            else:
                if len(ref_data["molecule"]["atoms"]) > 0:
                    atom_section = False

    # 1b. Parse Phase 2 energy charges from goed_gfnff output (first block only)
    # These are the charges used for Coulomb energy evaluation (nlist%q)
    energy_charges = []
    in_goed_charges = False
    natoms = len(ref_data["molecule"]["atoms"]) if ref_data["molecule"]["atoms"] else 0
    for line in lines:
        if "Charges computed in goed_gfnff:" in line and not energy_charges:
            in_goed_charges = True
            continue
        if in_goed_charges:
            if line.strip().startswith("---") or line.strip().startswith("Atom"):
                continue
            match = re.match(r"\s*\d+\s*\|\s*([-]?[\d.]+)", line)
            if match:
                energy_charges.append(float(match.group(1)))
                if natoms > 0 and len(energy_charges) >= natoms:
                    in_goed_charges = False
            else:
                if energy_charges:
                    in_goed_charges = False
    if energy_charges:
        ref_data["energy_charges"] = energy_charges

    # 2. Parse Energy Components
    for line in lines:
        if "Total energy:" in line and "Eh" in line:
            match = re.search(r"Total energy:\s+([-]?[\d.]+)", line)
            if match:
                ref_data["energy_components"]["total"] = float(match.group(1))

        if "Bond energy:" in line:
             match = re.search(r"Bond energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["bond"] = float(match.group(1))
        if "Angle energy:" in line:
             match = re.search(r"Angle energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["angle"] = float(match.group(1))
        if "Torsion energy:" in line:
             match = re.search(r"Torsion energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["torsion"] = float(match.group(1))
        if "Repulsion energy:" in line:
             match = re.search(r"Repulsion energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["repulsion"] = float(match.group(1))
        if "Electrostatic energy:" in line:
             match = re.search(r"Electrostatic energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["electrostatic"] = float(match.group(1))
        if "Dispersion energy:" in line:
             match = re.search(r"Dispersion energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["dispersion"] = float(match.group(1))
        if "Hydrogen bond energy:" in line:
             match = re.search(r"Hydrogen bond energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["hb"] = float(match.group(1))
        if "Halogen bond energy:" in line:
             match = re.search(r"Halogen bond energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["xb"] = float(match.group(1))
        if "Bonded atm energy:" in line:
             match = re.search(r"Bonded atm energy:\s+([-]?[\d.]+)", line)
             if match: ref_data["energy_components"]["batm"] = float(match.group(1))

        if "Gradient norm:" in line:
            match = re.search(r"Gradient norm:\s+([\d.]+)", line)
            if match:
                if "gradients" not in ref_data: ref_data["gradients"] = {}
                ref_data["gradients"]["norm"] = float(match.group(1))

    # 3. Parse Detailed Parameters (Bonds, Angles, Torsions)
    for i, line in enumerate(lines):
        if "Detailed Parameter Calculation for Bond" in line:
            match = re.search(r"Bond\s+(\d+)\s+-\s+(\d+)", line)
            if match:
                current_atoms = [int(match.group(1))-1, int(match.group(2))-1]
                r0, kbond, prefactor = None, None, None
                for j in range(i+1, min(i+100, len(lines))):
                    if "Detailed Parameter Calculation" in lines[j] and j > i+1: break
                    if "r0 =" in lines[j] and "Angstrom" in lines[j]:
                        r0_match = re.search(r"r0\s*=\s*.*?" + FLOAT_PATTERN + r"\s+Angstrom", lines[j])
                        if r0_match: r0 = float(r0_match.group(1))
                    if "topo%vbond(2,i) =" in lines[j]:
                        kbond_match = re.search(r"topo%vbond\(2,i\)\s*=\s*" + FLOAT_PATTERN + r"\s*$", lines[j])
                        if kbond_match: kbond = float(kbond_match.group(1))
                    if "topo%vbond(3,i) =" in lines[j]:
                        pre_match = re.search(r"topo%vbond\(3,i\)\s*=\s*" + FLOAT_PATTERN + r"\s*$", lines[j])
                        if pre_match: prefactor = float(pre_match.group(1))

                if r0 is not None:
                    ref_data["bonds"].append({
                        "atoms": current_atoms,
                        "r0": r0,
                        "kbond": kbond,
                        "prefactor": prefactor,
                        "energy": 0.0
                    })

        elif "Detailed Parameter Calculation for Angle" in line:
            match = re.search(r"Angle\s+(\d+)\s+-\s+(\d+)\s+-\s+(\d+)", line)
            if match:
                current_atoms = [int(match.group(1))-1, int(match.group(2))-1, int(match.group(3))-1]
                phi0, fc = None, None
                for j in range(i+1, min(i+100, len(lines))):
                    if "Detailed Parameter Calculation" in lines[j] and j > i+1: break
                    if "phi0 =" in lines[j] and "radians" in lines[j]:
                        phi0_match = re.search(r"phi0\s*=\s*.*?" + FLOAT_PATTERN + r"\s+radians", lines[j])
                        if phi0_match: phi0 = float(phi0_match.group(1)) * 180.0 / 3.14159265358979323846
                    if "topo%vangl(2,i) =" in lines[j]:
                        fc_match = re.search(r"topo%vangl\(2,i\)\s*=\s*" + FLOAT_PATTERN + r"\s*$", lines[j])
                        if fc_match: fc = float(fc_match.group(1))

                if phi0 is not None:
                    ref_data["angles"].append({
                        "atoms": current_atoms,
                        "phi0": phi0,
                        "fc": fc,
                        "energy": 0.0
                    })

        elif "Detailed Parameter Calculation for Torsion" in line:
            match = re.search(r"Torsion\s+(\d+)\s+-\s+(\d+)\s+-\s+(\d+)\s+-\s+(\d+)", line)
            if match:
                current_atoms = [int(match.group(1))-1, int(match.group(2))-1, int(match.group(3))-1, int(match.group(4))-1]
                phi0, fctot, n = None, None, None
                for j in range(i+1, min(i+100, len(lines))):
                    if "Detailed Parameter Calculation" in lines[j] and j > i+1: break
                    if "phi0 =" in lines[j] and "radians" in lines[j]:
                        phi0_match = re.search(r"phi0\s*=\s*.*?" + FLOAT_PATTERN + r"\s+radians", lines[j])
                        if phi0_match: phi0 = float(phi0_match.group(1))
                    if "fctot =" in lines[j]:
                        fc_match = re.search(r"=\s*" + FLOAT_PATTERN + r"\s*$", lines[j])
                        if fc_match: fctot = float(fc_match.group(1))
                    if "where n =" in lines[j]:
                        n_match = re.search(r"where n =\s+(\d+)", lines[j])
                        if n_match: n = int(n_match.group(1))

                if phi0 is not None:
                    ref_data["torsions"].append({
                        "atoms": current_atoms,
                        "n": n,
                        "phi0": phi0,
                        "V": fctot,
                        "energy": 0.0
                    })

    # 4. Parse Gradient Decomposition
    # Fortran format uses pipe-delimited table:
    #   "    1 | C  | Bond      |      -0.00000000       0.00000000      -0.00000000  0.00000"
    #   "      |      | Angle     |      -0.00000000       0.00000001      -0.00000000  0.00000"
    # Claude Generated (February 2026): Fixed parser to match pipe-delimited Fortran output
    grad_section = False
    current_atom_idx = -1
    for line in lines:
        if "=== Force Field Component Gradient Analysis ===" in line:
            grad_section = True
            continue
        if grad_section:
            if line.strip().startswith("---") or line.strip().startswith("Atom"):
                continue
            if not line.strip() or ("======" in line) or ("GFN-FF setup done" in line):
                if current_atom_idx >= 0:
                    grad_section = False
                continue
            if "|--" in line:
                continue

            # Split by pipe delimiter
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue

            # First field: atom index (e.g. "1") or empty (continuation)
            atom_field = parts[0].strip() if len(parts) > 0 else ""
            # Second field: element type (e.g. "C") or empty
            # Third field: component name (e.g. "Bond", "Angle")
            comp_field = parts[2].strip() if len(parts) > 2 else ""
            # Fourth field: gradient values (e.g. "  -0.00000000   0.00000000  -0.00000000  0.00000")
            grad_field = parts[3].strip() if len(parts) > 3 else ""

            if not comp_field or not grad_field:
                continue

            # Skip separator and header lines
            if comp_field in ("Component", "--------"):
                continue

            # Parse atom index from first field (if present, this is a new atom)
            if atom_field and atom_field.isdigit():
                current_atom_idx = int(atom_field) - 1

            if current_atom_idx < 0:
                continue

            # Parse gradient values from grad_field
            floats = re.findall(FLOAT_PATTERN, grad_field)
            if len(floats) >= 3:
                gx, gy, gz = float(floats[0]), float(floats[1]), float(floats[2])
                # Use string keys for JSON compatibility
                atom_key = str(current_atom_idx)
                if atom_key not in ref_data["gradient_decomposition"]:
                    ref_data["gradient_decomposition"][atom_key] = {}
                ref_data["gradient_decomposition"][atom_key][comp_field] = {"x": gx, "y": gy, "z": gz}

    return ref_data

def main():
    if len(sys.argv) < 3:
        print("Usage: ./gfnff_ref_generator.py <analyzer_bin> <molecule.xyz> [output.json]")
        sys.exit(1)

    analyzer_bin = sys.argv[1]
    xyz_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(xyz_file).replace(".xyz", ".ref.json")

    # Run analyzer
    cmd = [analyzer_bin, xyz_file, "-", "2", "2"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running analyzer: {e}")
        print(e.stdout)
        print(e.stderr)
        sys.exit(1)

    ref_data = parse_analyzer_output(output, xyz_file)

    with open(output_file, 'w') as f:
        json.dump(ref_data, f, indent=2)

    print(f"Generated reference data for {ref_data['molecule']['name']} -> {output_file}")

if __name__ == "__main__":
    main()
