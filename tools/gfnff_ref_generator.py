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

    # 4. Parse Gradient Decomposition (NEW)
    grad_section = False
    current_atom_idx = -1
    for line in lines:
        if "=== Force Field Component Gradient Analysis ===" in line:
            grad_section = True
            continue
        if grad_section:
            if "======" in line or "GFN-FF setup done" in line:
                grad_section = False
                continue

            # Match atom index line
            atom_match = re.match(r"^\s*(\d+)\s+([A-Z][a-z]?)\s+([A-Za-z]+)\s+" + FLOAT_PATTERN + r"\s+" + FLOAT_PATTERN + r"\s+" + FLOAT_PATTERN, line)
            if atom_match:
                current_atom_idx = int(atom_match.group(1)) - 1
                comp = atom_match.group(3)
                gx, gy, gz = float(atom_match.group(4)), float(atom_match.group(5)), float(atom_match.group(6))

                if current_atom_idx not in ref_data["gradient_decomposition"]:
                    ref_data["gradient_decomposition"][current_atom_idx] = {}
                ref_data["gradient_decomposition"][current_atom_idx][comp] = {"x": gx, "y": gy, "z": gz}
                continue

            # Match component line for the same atom
            comp_match = re.match(r"^\s+([A-Za-z]+)\s+" + FLOAT_PATTERN + r"\s+" + FLOAT_PATTERN + r"\s+" + FLOAT_PATTERN, line)
            if comp_match and current_atom_idx != -1:
                comp = comp_match.group(1)
                gx, gy, gz = float(comp_match.group(2)), float(comp_match.group(3)), float(comp_match.group(4))
                ref_data["gradient_decomposition"][current_atom_idx][comp] = {"x": gx, "y": gy, "z": gz}

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
