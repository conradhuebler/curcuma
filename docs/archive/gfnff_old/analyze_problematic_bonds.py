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

print("Problematic bond distances that should NOT be detected as bonds:")
print("================================================================")

# These are the incorrect bonds being detected
problematic_bonds = [
    (0, 5),  # C1-H4
    (0, 6),  # C1-H5
    (0, 7),  # C1-H6
    (1, 2),  # C2-H1
    (1, 3),  # C2-H2
    (1, 4),  # C2-H3
]

atom_names = ['C1', 'C2', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6']

for i, j in problematic_bonds:
    dist = distance(coords[i], coords[j])
    r_cov_sum = covalent_radii[6] + covalent_radii[1]  # C-H
    bond_threshold = 1.3 * r_cov_sum
    
    print(f"{atom_names[i]}-{atom_names[j]}: {dist:.3f} Å (cov sum: {r_cov_sum:.3f}, threshold: {bond_threshold:.3f}) {'[INCORRECTLY DETECTED AS BOND!]' if dist < bond_threshold else ''}")

print("\nActual bond distances that SHOULD be detected as bonds:")
print("======================================================")

# These are the correct bonds
correct_bonds = [
    (0, 1),  # C1-C2
    (0, 2),  # C1-H1
    (0, 3),  # C1-H2
    (0, 4),  # C1-H3
    (1, 5),  # C2-H4
    (1, 6),  # C2-H5
    (1, 7),  # C2-H6
]

for i, j in correct_bonds:
    dist = distance(coords[i], coords[j])
    if i == 0 and j == 1:  # C-C bond
        r_cov_sum = covalent_radii[6] + covalent_radii[6]
    else:  # C-H bonds
        r_cov_sum = covalent_radii[6] + covalent_radii[1]
    bond_threshold = 1.3 * r_cov_sum
    
    print(f"{atom_names[i]}-{atom_names[j]}: {dist:.3f} Å (cov sum: {r_cov_sum:.3f}, threshold: {bond_threshold:.3f}) {'[CORRECTLY DETECTED AS BOND]' if dist < bond_threshold else ''}")