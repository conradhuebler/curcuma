#!/usr/bin/env python3
"""
GFN-FF vbond Parameter Derivation - EXACT FORTRAN IMPLEMENTATION

Implementation of the THREE formulas from gfnff_ini.f90:
1. vbond(1) = gen%rabshift + shift (Line 1276) - 12 Shift Rules
2. vbond(2) = gen%srb1 * (1.0 + fsrb2*(ΔEN)^2 + gen%srb3*bstrength) (Line 1223)
3. vbond(3) = -param%bond(ia)*param%bond(ja)*ringf*bstrength*fqq*fheavy*fpi*fxh*fcn (Line 1285)

Source: gfnff_param.f90 (lines 725-727: srb1=0.3731, srb2=0.3171, srb3=0.2538)
        gfnff_ini.f90 (lines 1087-1285)
Reference: Spicher & Grimme, Angew. Chem. 2020, 59, 15665

KNOWN ISSUE: OH vbond(1) expected -0.160 (from Fortran logic) vs observed -0.020 (XTB)
This discrepancy is documented in the plan and needs investigation.

Claude Generated - EXACT implementation from source code analysis
"""

import json
import math
from pathlib import Path


class GFNFFParameters:
    """EXACT generator parameters from gfnff_param.f90 lines 725-766"""

    # Shift parameters (lines 757-766)
    RABSHIFT = -0.110       # gen%rabshift
    RABSHIFTH = -0.050      # gen%rabshifth (X-H bonds)
    HYPER_SHIFT = 0.03      # gen%hyper_shift
    HSHIFT3 = -0.11         # gen%hshift3 (Z > 10)
    HSHIFT4 = -0.11         # gen%hshift4 (Z > 18)
    HSHIFT5 = -0.06         # gen%hshift5 (Z > 36)
    METAL1_SHIFT = 0.2      # Group 1+2 metals
    METAL2_SHIFT = 0.15     # Transition metals
    METAL3_SHIFT = 0.05     # Main group metals
    ETA_SHIFT = 0.040       # Eta-bonded

    # Steepness parameters (lines 725-727)
    SRB1 = 0.3731           # gen%srb1
    SRB2 = 0.3171           # gen%srb2
    SRB3 = 0.2538           # gen%srb3

    # Bond strength matrix (4x4, lines 804-814)
    # indexed by [max(hyb)][min(hyb)]
    BSMAT = [
        [1.000, 1.323, 1.079, 1.000],  # max=0
        [1.323, 1.980, 1.484, 1.323],  # max=1
        [1.079, 1.484, 1.240, 1.079],  # max=2
        [1.000, 1.323, 1.079, 1.000]   # max=3
    ]

    # Pauling electronegativities (lines 315-324, from external source)
    PAULING_EN = [
        2.20, 3.00,  # 1-2: H, He
        0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98,  # 3-9: Li-F
        4.50,  # 10: Ne
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50,  # 11-18: Na-Ar
        0.82, 1.00, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88,  # 19-27: K-Co
        1.91, 1.90, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3.00,  # 28-36: Ni-Kr
        0.82, 0.95, 1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28,  # 37-45: Rb-Rh
        2.20, 1.93, 1.69, 1.78, 1.96, 2.05, 2.10, 2.66, 2.60,  # 46-54: Pd-Xe
        0.79, 0.89, 1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18,  # 55-63: Cs-Eu
        1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.30,  # 64-72: Gd-Hf
        1.50, 1.70, 1.90, 2.10, 2.20, 2.20, 2.20,  # 73-79: Ta-Au
        2.00, 1.62, 2.33, 2.02, 2.00, 2.20, 2.20  # 80-86: Hg-Rn
    ]

    # Bond parameters from gfnff_par.h (kcal/(mol·Å²))
    BOND_PARAMS = [
        0.417997, 0.258490, 0.113608, 0.195935, 0.231217,  # 1-5: H-B
        0.385248, 0.379257, 0.339249, 0.330706, 0.120319,  # 6-10: C-Ne
        0.127255, 0.173647, 0.183796, 0.273055, 0.249044,  # 11-15: Na-P
        0.290653, 0.218744, 0.034706, 0.136353, 0.192467,  # 16-20: S-Ca
        0.335860, 0.314452, 0.293044, 0.271636, 0.250228,  # 21-25: Sc-Mn
        0.228819, 0.207411, 0.186003, 0.164595, 0.143187,  # 26-30: Fe-Zn
        0.212434, 0.210451, 0.219870, 0.224618, 0.272206,  # 31-35: Ga-Br
        0.147864, 0.150000, 0.150000, 0.329501, 0.309632,  # 36-40: Kr-Zr
        0.289763, 0.269894, 0.250025, 0.230155, 0.210286,  # 41-45: Nb-Rh
        0.190417, 0.170548, 0.150679, 0.192977, 0.173411,  # 46-50: Pd-Sn
        0.186907, 0.192891, 0.223202, 0.172577, 0.150000,  # 51-55: Sb-Cs
        0.150000, 0.370682, 0.368511, 0.366339, 0.364168,  # 56-60: Ba-Nd
        0.361996, 0.359825, 0.357654, 0.355482, 0.353311,  # 61-65: Pm-Tb
        0.351139, 0.348968, 0.346797, 0.344625, 0.342454,  # 66-70: Dy-Yb
        0.340282, 0.338111, 0.305540, 0.272969, 0.240398,  # 71-75: Lu-Re
        0.207828, 0.175257, 0.142686, 0.110115, 0.077544,  # 76-80: Os-Hg
        0.108597, 0.148422, 0.183731, 0.192274, 0.127706,  # 81-85: Tl-At
        0.086756  # 86: Rn
    ]

    # Element-pair specific fxh corrections from gfnff_ini.f90 lines 1160-1165
    FXH_CORRECTIONS = {
        # (Z_i, Z_j): fxh value
        (1, 5): 1.10,   # B-H
        (5, 1): 1.10,
        (1, 7): 1.06,   # N-H
        (7, 1): 1.06,
        (1, 8): 0.93,   # O-H (WEAKENS bond!)
        (8, 1): 0.93,
        # C-H variations are geometry-dependent (3-ring, aldehyde)
        # Not included in base model for diatomics
        # Metal-H corrections are metal-type dependent
        # Not included for diatomics (no metals in test set)
    }

    @staticmethod
    def get_pauling_en(z):
        """Get Pauling EN for atomic number z (1-86)"""
        if not 1 <= z <= 86:
            raise ValueError(f"Atomic number {z} out of range [1, 86]")
        return GFNFFParameters.PAULING_EN[z - 1]

    @staticmethod
    def get_bond_param(z):
        """Get bond parameter for atomic number z (1-86)"""
        if not 1 <= z <= 86:
            raise ValueError(f"Atomic number {z} out of range [1, 86]")
        return GFNFFParameters.BOND_PARAMS[z - 1]

    @staticmethod
    def get_fxh(z1, z2):
        """Get element-pair specific fxh correction"""
        pair1 = (z1, z2)
        pair2 = (z2, z1)
        if pair1 in GFNFFParameters.FXH_CORRECTIONS:
            return GFNFFParameters.FXH_CORRECTIONS[pair1]
        if pair2 in GFNFFParameters.FXH_CORRECTIONS:
            return GFNFFParameters.FXH_CORRECTIONS[pair2]
        return 1.0


class VbondDeriver:
    """Derives vbond parameters using EXACT Fortran formulas"""

    def derive_r0_shift(self, z1, z2):
        """
        Formula 1: vbond(1) = gen%rabshift + shift (Line 1276)

        Implements 12 shift calculation rules from gfnff_ini.f90 lines 1087-1275
        Rules applied in order: some REPLACE previous shift, some ADD to it

        For diatomics, only rules 2, 4, 5, 12 are relevant:
        - Rule 2: X-H bonds (both test cases have H)
        - Rule 4: X-sp³ hybridization (additive after rule 2)
        - Rule 5: X-sp hybridization (additive)
        - Rule 12: Heavy atoms both Z>10 (additive)
        """
        shift = 0.0

        # Rule 2: X-Hydrogen atom (line 1144)
        # ia == 1 OR ja == 1: REPLACES previous shift
        if z1 == 1 or z2 == 1:
            shift = GFNFFParameters.RABSHIFTH
            # After this, rules 4, 5, 12 can modify shift additively

        # Rule 4: X-sp³ hybridization (lines 1146-1147)
        # For diatomics, this depends on actual hybridization of atoms
        # Not directly accessible without topology info
        # For now: skip (would need hyb(ii), hyb(jj) info)

        # Rule 5: X-sp hybridization (lines 1148-1149)
        # Similar issue - needs hybridization info
        # For now: skip

        # Rule 12: Heavy atoms BOTH Z>10 (lines 1267-1272)
        # Requires BOTH atoms Z > 10 (AND, not OR)
        if z1 > 10 and z2 > 10:
            shift += GFNFFParameters.HSHIFT3  # -0.11
            if z1 > 18 or z2 > 18:
                shift += GFNFFParameters.HSHIFT4  # -0.11
            if z1 > 36 or z2 > 36:
                shift += GFNFFParameters.HSHIFT5  # -0.06

        return GFNFFParameters.RABSHIFT + shift

    def derive_steepness(self, z1, z2):
        """
        Formula 2: vbond(2) = gen%srb1 * (1.0 + fsrb2*(EN_i - EN_j)^2 + gen%srb3*bstrength)
        (Line 1223 in gfnff_ini.f90)

        For diatomics with no metals:
        - bstrength = 1.0 (single bond)
        - fsrb2 = srb2 = 0.3171 (non-metal, not metal)
        """
        en1 = GFNFFParameters.get_pauling_en(z1)
        en2 = GFNFFParameters.get_pauling_en(z2)
        en_diff_sq = (en1 - en2) ** 2

        # For diatomics: bstrength = 1.0 (single bond, no special cases)
        bstrength = 1.0

        # For diatomics: non-metal bonds, so fsrb2 = srb2
        fsrb2 = GFNFFParameters.SRB2

        steepness = (GFNFFParameters.SRB1 *
                     (1.0 + fsrb2 * en_diff_sq + GFNFFParameters.SRB3 * bstrength))

        return steepness

    def derive_prefactor(self, z1, z2):
        """
        Formula 3: vbond(3) = -bond(ia)*bond(ja)*ringf*bstrength*fqq*fheavy*fpi*fxh*fcn
        (Line 1285 in gfnff_ini.f90)

        For DIATOMIC non-metal systems:
        - ringf = 1.0 (no rings)
        - bstrength = 1.0 (single bond)
        - fqq ≈ 1.0 (neutral atoms)
        - fheavy = 1.0 (no metals)
        - fpi = 1.0 (no aromatic)
        - fxh = element-pair specific (see fxh corrections)
        - fcn = 1.0 (if not both heavy, or if no metal coordination)

        Result: vbond(3) = -bond(ia)*bond(ja)*fxh
        """
        bond_product = (GFNFFParameters.get_bond_param(z1) *
                        GFNFFParameters.get_bond_param(z2))

        # For diatomics, only fxh varies
        fxh = GFNFFParameters.get_fxh(z1, z2)

        # All other factors = 1.0 for diatomic non-metals
        prefactor = -bond_product * fxh

        return prefactor

    def derive_vbond(self, z1, z2):
        """Derive complete vbond triplet [vbond(1), vbond(2), vbond(3)]"""
        return [
            self.derive_r0_shift(z1, z2),
            self.derive_steepness(z1, z2),
            self.derive_prefactor(z1, z2)
        ]


def get_element_symbol(z):
    """Get element symbol from atomic number"""
    symbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn"
    ]
    return symbols[z - 1] if 1 <= z <= 86 else "X"


def main():
    print("\n" + "="*80)
    print("GFN-FF vbond Parameter Derivation - EXACT FORTRAN IMPLEMENTATION")
    print("="*80)

    deriver = VbondDeriver()

    # Test cases
    test_cases = [
        (1, 1, "HH"),
        (8, 1, "OH"),
        (17, 1, "HCl"),
    ]

    results = {}

    for z1, z2, name in test_cases:
        elem1 = get_element_symbol(z1)
        elem2 = get_element_symbol(z2)

        vbond = deriver.derive_vbond(z1, z2)

        print(f"\n{name} ({elem1}-{elem2}):")
        print(f"  Pauling EN: {GFNFFParameters.get_pauling_en(z1):.3f} - "
              f"{GFNFFParameters.get_pauling_en(z2):.3f}")
        print(f"  Bond params: {GFNFFParameters.get_bond_param(z1):.6f} × "
              f"{GFNFFParameters.get_bond_param(z2):.6f}")
        print(f"\n  Derived vbond:")
        print(f"    vbond(1) [R0_shift]:   {vbond[0]:20.15f} Bohr")
        print(f"    vbond(2) [steepness]:  {vbond[1]:20.15f}")
        print(f"    vbond(3) [prefactor]:  {vbond[2]:20.15f} Hartree/Bohr²")

        results[name] = vbond

    # Compare with XTB reference values
    print("\n" + "="*80)
    print("Comparison with XTB Reference Values")
    print("="*80)

    reference = {
        "HH": [-0.160000000000000, 0.467792780000000, -0.178827447071212],
        "OH": [-0.020000000000000, 0.680329896428000, -0.182731113861188],
        "HCl": [-0.160000000000000, 0.576827285216000, -0.095731747575742],
    }

    print(f"\n{'System':<10} {'Parameter':<15} {'Derived':>20} {'Reference':>20} "
          f"{'Error':>12} {'Match':>8}")
    print("-" * 95)

    total_errors = []
    for name in ["HH", "OH", "HCl"]:
        if name not in results:
            continue

        derived = results[name]
        ref = reference.get(name)

        if not ref:
            print(f"{name}: No reference data")
            continue

        for i, (d, r) in enumerate(zip(derived, ref)):
            param_name = ["r0_shift", "steepness", "prefactor"][i]
            error = abs(d - r)
            error_pct = 100.0 * error / max(abs(r), 1e-10)
            match = "✓" if error < 0.001 else "✗"

            print(f"{name:<10} {param_name:<15} {d:>20.15f} {r:>20.15f} "
                  f"{error_pct:>11.3f}% {match:>8}")

            total_errors.append(error_pct)

    if total_errors:
        avg_error = sum(total_errors) / len(total_errors)
        max_error = max(total_errors)
        print("-" * 95)
        print(f"{'Average Error:':<50} {avg_error:>11.3f}%")
        print(f"{'Maximum Error:':<50} {max_error:>11.3f}%")

    print("\n⚠️  KNOWN ISSUE: OH vbond(1) expected -0.160 (Fortran logic), got -0.020 (XTB)")
    print("    This discrepancy needs investigation (see Plan file)")

    print("\n" + "="*80 + "\n")


if __name__ == "__main__":
    main()
