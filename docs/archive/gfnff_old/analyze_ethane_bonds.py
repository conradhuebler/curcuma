import numpy as np

# Ethane coordinates from ethane.xyz
coords = np.array([
    [0.000000, 0.000000, 0.000000],  # C1 (atom 0)
    [0.000000, 0.000000, 1.540000],  # C2 (atom 1)
    [0.946000, 0.000000, -0.370000], # H1 (atom 2)
    [-0.473000, 0.819000, -0.370000], # H2 (atom 3)
    [-0.473000, -0.819000, -0.370000], # H3 (atom 4)
    [0.946000, 0.000000, 1.910000],  # H4 (atom 5)
    [-0.473000, 0.819000, 1.910000], # H5 (atom 6)
    [-0.473000, -0.819000, 1.910000]  # H6 (atom 7)
])

# Covalent radii from gfnff_par.h (in Angstrom)
# H: 0.32, C: 0.75
covalent_radii = {
    1: 0.32,  # H
    6: 0.75   # C
}

def distance(a, b):
    return np.sqrt(np.sum((a - b)**2))

print("Ethane Bond Analysis")
print("====================")

# Calculate all distances
atom_types = [6, 6, 1, 1, 1, 1, 1, 1]  # C, C, H, H, H, H, H, H

print("\nAll interatomic distances (Angstrom):")
for i in range(8):
    for j in range(i+1, 8):
        dist = distance(coords[i], coords[j])
        r_cov_sum = covalent_radii[atom_types[i]] + covalent_radii[atom_types[j]]
        bond_threshold = 1.3 * r_cov_sum
        
        is_bond = dist < bond_threshold
        
        atom_names = ['C1', 'C2', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6']
        print(f"{atom_names[i]}-{atom_names[j]}: {dist:.3f} Å (cov sum: {r_cov_sum:.3f}, threshold: {bond_threshold:.3f}) {'[BOND]' if is_bond else ''}")

print("\nExpected bonds:")
print("- C1-C2 (bond distance ~1.54 Å)")
print("- C1-H1, C1-H2, C1-H3 (bond distances ~1.09 Å)")
print("- C2-H4, C2-H5, C2-H6 (bond distances ~1.09 Å)")