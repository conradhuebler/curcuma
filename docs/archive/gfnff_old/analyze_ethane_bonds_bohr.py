import numpy as np

# Ethane coordinates from ethane.xyz (in Angstrom)
coords_angstrom = np.array([
    [0.000000, 0.000000, 0.000000],  # C1 (atom 0)
    [0.000000, 0.000000, 1.540000],  # C2 (atom 1)
    [0.946000, 0.000000, -0.370000], # H1 (atom 2)
    [-0.473000, 0.819000, -0.370000], # H2 (atom 3)
    [-0.473000, -0.819000, -0.370000], # H3 (atom 4)
    [0.946000, 0.000000, 1.910000],  # H4 (atom 5)
    [-0.473000, 0.819000, 1.910000], # H5 (atom 6)
    [-0.473000, -0.819000, 1.910000]  # H6 (atom 7)
])

# Convert to Bohr (1 Angstrom = 1.889726 Bohr)
ANGSTROM_TO_BOHR = 1.889726
coords_bohr = coords_angstrom * ANGSTROM_TO_BOHR

# Covalent radii from gfnff_par.h (in Angstrom)
# H: 0.32, C: 0.75
covalent_radii_angstrom = {
    1: 0.32,  # H
    6: 0.75   # C
}

# Convert to Bohr
covalent_radii_bohr = {
    1: 0.32 * ANGSTROM_TO_BOHR,  # H
    6: 0.75 * ANGSTROM_TO_BOHR   # C
}

def distance(a, b):
    return np.sqrt(np.sum((a - b)**2))

print("Ethane Bond Analysis (BOHR UNITS)")
print("=================================")

# Calculate all distances in Bohr
atom_types = [6, 6, 1, 1, 1, 1, 1, 1]  # C, C, H, H, H, H, H, H

print("\nAll interatomic distances (Bohr):")
for i in range(8):
    for j in range(i+1, 8):
        dist = distance(coords_bohr[i], coords_bohr[j])
        r_cov_sum = covalent_radii_bohr[atom_types[i]] + covalent_radii_bohr[atom_types[j]]
        bond_threshold = 1.3 * r_cov_sum
        
        is_bond = dist < bond_threshold
        
        atom_names = ['C1', 'C2', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6']
        print(f"{atom_names[i]}-{atom_names[j]}: {dist:.3f} Bohr (cov sum: {r_cov_sum:.3f}, threshold: {bond_threshold:.3f}) {'[BOND]' if is_bond else ''}")

print("\nExpected bonds:")
print("- C1-C2 (bond distance ~2.907 Bohr)")
print("- C1-H1, C1-H2, C1-H3 (bond distances ~1.918 Bohr)")
print("- C2-H4, C2-H5, C2-H6 (bond distances ~1.918 Bohr)")