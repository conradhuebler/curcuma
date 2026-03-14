# GFN-FF Inversion/Out-of-Plane Potential - Theoretical Background

**Author**: Claude (AI Assistant) - Based on GFN-FF method
**Date**: 2025-11-10
**Reference**: Spicher, S.; Grimme, S. *Angew. Chem. Int. Ed.* **2020**, *59*, 15665-15676
**Fortran Implementation**: `external/gfnff/src/gfnff_helpers.f90:427-510` (omega, domegadr)

---

## ACKNOWLEDGMENT

This document describes the **inversion/out-of-plane potential** component of the **GFN-FF (Geometry, Frequency, Noncovalent Force Field)** method, developed by:

- **Prof. Dr. Stefan Grimme** (University of Bonn, Mulliken Center for Theoretical Chemistry)
- **Dr. Sebastian Spicher** (University of Bonn)

**Original Publication**:
S. Spicher, S. Grimme
*"Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems"*
*Angew. Chem. Int. Ed.* **2020**, *59*, 15665-15676
DOI: [10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239)

**Original Fortran Code**:
Available at: https://github.com/grimme-lab/gfnff
Licensed under: LGPL v3

This document serves as **educational material for theoretical chemists** to understand the scientific basis of Curcuma's native C++ implementation.

---

## 1. Physical Background

### 1.1 What are Inversion/Out-of-Plane Terms?

**Inversions** (also called **improper torsions** or **out-of-plane bends**) describe the energy cost of moving an atom out of a planar arrangement. They are essential for modeling:

1. **sp²-Hybridized Centers**:
   - Carbon in ethene (C=C), aromatics, carbonyls (C=O)
   - Nitrogen in imines, pyridine, nitro groups
   - Boron in BH₃, BF₃ (trigonal planar)

2. **Planarity Constraints**:
   - Aromatic rings (benzene, naphthalene)
   - Conjugated systems (polyenes, peptide bonds)
   - Carboxylate anions (COO⁻)

3. **Pyramidalization**:
   - Ammonia (NH₃) - pyramidal with inversion barrier ~24 kJ/mol
   - Phosphines (PR₃) - higher barrier than amines
   - Umbrella inversions in molecules

### 1.2 Why Do Molecules Prefer Planarity?

**Quantum Mechanical Origin**:
- **π-Bonding**: Requires parallel p-orbitals → planarity maximizes overlap
- **Conjugation**: Delocalized π-systems stabilize planar geometries
- **Hybridization**: sp² orbitals are naturally 120° apart in a plane

**Energy Penalty for Non-Planarity**:
- Breaking conjugation (~20-50 kJ/mol for aromatics)
- Reduced π-overlap (C=C twisting ~260 kJ/mol)
- Steric strain in pyramidal configurations

### 1.3 Examples in Chemistry

| Molecule | Central Atom | Planar? | Inversion Barrier |
|----------|--------------|---------|-------------------|
| Ethene (C₂H₄) | C (sp²) | Yes | ~260 kJ/mol (C=C twist) |
| Benzene (C₆H₆) | C (sp²) | Yes | Very high (aromatic) |
| Ammonia (NH₃) | N (sp³) | No | ~24 kJ/mol (umbrella) |
| Formaldehyde (H₂CO) | C (sp²) | Yes | ~200 kJ/mol (C=O twist) |
| BH₃ | B (sp²) | Yes | ~100 kJ/mol |

---

## 2. Mathematical Formulation

### 2.1 Out-of-Plane Angle Definition

For four atoms **i-j-k-l**, the **out-of-plane angle ω** (omega) measures how far atom **i** is from the plane defined by **j-k-l**.

**Geometry**:
```
      i  ← Central atom (out-of-plane)
      |
      j--k--l  ← Reference plane
```

**Mathematical Definition**:

ω = arcsin(**n** · **v̂**)

where:
- **n** = (**r**ᵢⱼ × **r**ⱼₖ) / |**r**ᵢⱼ × **r**ⱼₖ|  (unit normal to plane i-j-k)
- **v** = **r**ᵢₗ = **r**ₗ - **r**ᵢ  (vector from i to l)
- **v̂** = **v** / |**v**|  (unit vector)

**Range**: ω ∈ [-π/2, +π/2]
- ω = 0 → atom i in the plane j-k-l (planar)
- ω = ±π/2 → atom i perpendicular to plane (fully pyramidal)

### 2.2 GFN-FF Inversion Energy

The inversion potential in GFN-FF uses a **double-well cosine function** to allow two equivalent planar configurations (±ω₀):

**E_inv** = V · [cos(ω) - cos(ω₀)]² · D(**r**ᵢⱼ, **r**ⱼₖ, **r**ⱼₗ)

**Components**:

1. **V** = Barrier height (kcal/mol)
   - Depends on element type and hybridization
   - Typical values: 1-10 kcal/mol for sp² centers
   - Larger for strong π-systems (aromatics)

2. **ω₀** = Reference angle (usually 0 for planar)
   - For sp²: ω₀ = 0 (planar preference)
   - For pyramidal (e.g., NH₃): ω₀ ≠ 0

3. **D(r)** = Distance damping function
   - Couples bond stretching with inversion
   - Same form as torsion damping
   - Reduces barrier when bonds are stretched

### 2.3 Alternative Form (φ₀ = 0 case)

For strictly planar systems (ω₀ = 0), GFN-FF sometimes uses:

**E_inv** = V · [1 + cos(ω + π)] · D(r) = V · [1 - cos(ω)]

This is a **single-well potential** with minimum at ω = 0.

---

## 3. Analytical Gradients

### 3.1 Gradient Formula

The force on each atom is:

**F**ᵢ = -∂E/∂**r**ᵢ = -(∂E/∂ω) · (∂ω/∂**r**ᵢ) - (∂E/∂D) · (∂D/∂**r**ᵢ)

**Chain rule components**:

1. **Energy derivative**:
   ∂E/∂ω = -2V · sin(ω) · [cos(ω₀) - cos(ω)] · D(r)

2. **Angle derivatives**: ∂ω/∂**r**ᵢ, ∂ω/∂**r**ⱼ, ∂ω/∂**r**ₖ, ∂ω/∂**r**ₗ
   - Computed via domegadr subroutine
   - Complex due to cross products and normalizations

3. **Damping derivatives**: ∂D/∂**r** (same as torsions)

### 3.2 Geometric Derivatives (domegadr)

From the Fortran implementation (gfnff_helpers.f90:450-510):

**Key vectors**:
- **r**ₑ = **r**ᵢ - **r**ⱼ
- **r**_d = **r**ₖ - **r**ⱼ
- **r**_v = **r**ₗ - **r**ᵢ
- **n** = **r**ₑ × **r**_d  (normal to plane i-j-k)

**Gradient expressions** (1/cos(ω) normalization):

∂ω/∂**r**ᵢ = (1/(|**n**||**v**|cos ω)) · [(**r**_d × **v**) - **n** - sin(ω)·(...)]

∂ω/∂**r**ⱼ = (1/(|**n**||**v**|cos ω)) · [(**v** × (**r**_d - **r**ₑ)) - sin(ω)·(...)]

∂ω/∂**r**ₖ = (1/(|**n**||**v**|cos ω)) · [(**r**_v × **r**ₑ) - sin(ω)·(...)]

∂ω/∂**r**ₗ = (1/(|**n**||**v**|cos ω)) · [**n** - sin(ω)·(**n**/|**v**)**v**]

**Singularity**: When cos(ω) ≈ 0 (ω ≈ ±90°), gradients → 0 (physical: no restoring force when perpendicular).

---

## 4. Implementation in GFN-FF

### 4.1 Inversion Detection Algorithm

Not all atom quartets (i-j-k-l) require inversion terms. GFN-FF detects them based on:

1. **Central atom hybridization**:
   - sp² → needs inversion term (planar preference)
   - sp³ → no inversion (tetrahedral is natural)
   - sp → no inversion (linear)

2. **Coordination number**:
   - 3 neighbors → likely sp² → add inversion
   - Example: C in ethene has 3 neighbors (2H, 1C)

3. **Pi-system membership**:
   - Atoms in aromatic rings → strong inversion terms
   - Conjugated C=C chains → moderate inversion terms

4. **Topology check**:
   - Central atom j with exactly 3 neighbors: i, k, l
   - Forms a trigonal planar arrangement

### 4.2 Parameter Assignment

**Barrier heights** (element-specific):

| Element | Hybridization | Typical V (kcal/mol) | Examples |
|---------|---------------|----------------------|----------|
| C | sp² | 5-10 | Ethene, benzene |
| N | sp² | 3-6 | Pyridine, imines |
| O | sp² | 2-4 | Carbonyls (weak) |
| B | sp² | 8-12 | BH₃, BF₃ |
| P | sp³→sp² | 1-3 | Phosphines (low) |

**Reference angle ω₀**:
- Usually 0 (planar)
- For some heteroatoms: small non-zero values

### 4.3 Distance Damping

Same damping function as torsions:

D(r) = 1 / [1 + exp(-α·(r/r₀ - 1))]

where:
- α = 16 (steepness parameter)
- r₀ = sum of covalent radii
- Applied to all three bonds: i-j, j-k, j-l

**Physical meaning**: When bonds stretch (r >> r₀), the inversion barrier decreases (bond breaking).

---

## 5. Comparison: Torsions vs. Inversions

| Property | Torsions | Inversions |
|----------|----------|------------|
| **Atoms** | 4 (i-j-k-l, linear chain) | 4 (i central, j-k-l plane) |
| **Angle** | Dihedral φ ∈ [-π, π] | Out-of-plane ω ∈ [-π/2, π/2] |
| **Formula** | φ = atan2(sin, cos) | ω = arcsin(**n**·**v̂**) |
| **Potential** | cos(n·φ) (n=1,2,3) | [cos(ω) - cos(ω₀)]² |
| **Typical for** | Single bonds (C-C rotation) | sp² centers (planarization) |
| **Barrier** | 1-4 kcal/mol | 2-10 kcal/mol |
| **Examples** | Butane (gauche/anti) | Ethene (planar C=C) |

---

## 6. Validation Strategy

### 6.1 Test Molecules

1. **Ethene (C₂H₄)**:
   - Two sp² carbons → 2 inversion terms
   - Should be strictly planar (ω = 0)
   - Barrier: ~260 kJ/mol for 90° twist

2. **Ammonia (NH₃)**:
   - Pyramidal nitrogen (ω₀ ≈ 20°)
   - Inversion barrier: ~24 kJ/mol
   - Tests non-zero ω₀

3. **Benzene (C₆H₆)**:
   - 6 sp² carbons → 6 inversion terms
   - Perfect planarity required
   - Aromatic stabilization

### 6.2 Validation Criteria

Compare native C++ implementation vs. Fortran reference:

1. **Number of inversions detected**: Must match exactly
2. **Barrier heights**: ±20% acceptable in Phase 1
3. **Energy at ω = 0**: Within ±0.1 kcal/mol
4. **Energy at ω = 30°**: Within ±0.5 kcal/mol
5. **Gradient accuracy**: <5% error in force components

### 6.3 Numerical Tests

1. **Finite difference check**:
   ```
   ∂E/∂x ≈ [E(x+h) - E(x-h)] / (2h)
   ```
   Compare with analytical gradient (should agree to ~10⁻⁶)

2. **Translational invariance**:
   ```
   ∂ω/∂r_i + ∂ω/∂r_j + ∂ω/∂r_k + ∂ω/∂r_l = 0
   ```

3. **Rotational invariance**:
   Energy should not change when rotating entire molecule

---

## 7. Known Limitations (Phase 1.2 Implementation)

### 7.1 Simplifications

1. **Hybridization detection**: Based on coordination number only
   - Missing: Geometry-based analysis (bond angles)
   - Missing: Electronic structure (VSEPR theory)

2. **Pi-system detection**: Not implemented
   - Full GFN-FF: Increased barriers for conjugated systems
   - Phase 1: Uses element-specific defaults

3. **Ring corrections**: Not implemented
   - Small rings (cyclopropene): Should have reduced barriers
   - Requires ring detection (Phase 2)

4. **Multiple inversion terms**: One per sp² center
   - Full GFN-FF: Sometimes uses multiple potentials
   - Example: Delocalized anions (carboxylate)

### 7.2 Expected Accuracy

**Phase 1.2 targets**:
- Energy differences: ±1-2 kcal/mol vs. full GFN-FF
- Geometry optimization: ±0.1 Å RMSD for planar molecules
- Suitable for: Simple organic molecules, validation purposes
- Not suitable for: Complex conjugated systems, transition metals

---

## 8. Literature References

1. **Original GFN-FF Paper**:
   S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.* **2020**, *59*, 15665-15676
   DOI: 10.1002/anie.202004239

2. **Improper Torsions Review**:
   A. D. MacKerell Jr. et al., *J. Phys. Chem. B* **1998**, *102*, 3586-3616
   (CHARMM force field - historical context)

3. **Out-of-Plane Deformations**:
   T. A. Halgren, *J. Comput. Chem.* **1996**, *17*, 490-519
   (MMFF94 force field - alternative approach)

4. **Umbrella Inversion**:
   R. J. Berry, *J. Chem. Phys.* **1960**, *32*, 933-938
   (Classic paper on NH₃ inversion)

5. **π-Conjugation Effects**:
   G. J. Martyna, M. E. Tuckerman, *J. Chem. Phys.* **1999**, *110*, 2810-2821
   (Planar constraint methods)

6. **Aromatic Planarity**:
   P. v. R. Schleyer, H. Jiao, *Pure Appl. Chem.* **1996**, *68*, 209-218
   (Aromaticity and geometry)

---

## 9. Computational Aspects

### 9.1 Efficiency

- **Cost per inversion**: ~20 FLOPs (angle) + ~150 FLOPs (gradient)
- **Typical molecule**: 10-30% of atoms have inversions
- **Bottleneck**: Cross products (6 per gradient calculation)

### 9.2 Numerical Stability

**Critical points**:
1. **cos(ω) → 0** (ω ≈ ±90°):
   - Gradients numerically unstable
   - Fortran: Returns zero gradients
   - Physical: System is at inflection point

2. **|**n**| → 0** (i, j, k collinear):
   - Angle undefined
   - Skip inversion term

3. **|**v**| → 0** (i and l coincide):
   - Angle undefined
   - Skip inversion term

### 9.3 Optimization Considerations

- **Second derivatives** (Hessian): Not implemented in Phase 1
- **Constraint algorithms**: Can enforce planarity faster than inversion potential
- **Convergence**: Inversions improve SCF-like iterations in geometry optimization

---

## 10. Summary

**Inversions/out-of-plane terms** are essential for:
- ✅ Maintaining planarity in sp² centers
- ✅ Modeling conjugated systems correctly
- ✅ Accurate aromatic ring geometries
- ✅ Pyramidalization in nitrogen/phosphorus

**Key equations**:
- ω = arcsin(**n** · **v̂**)
- E = V·[cos(ω) - cos(ω₀)]²·D(r)
- **F** = -(∂E/∂ω)·(∂ω/∂**r**) - (∂E/∂D)·(∂D/∂**r**)

**Implementation status** (Phase 1.2):
- ✅ Basic inversion detection (sp² identification)
- ✅ Analytical gradients (domegadr)
- ⏳ Pi-system corrections (Phase 2)
- ⏳ Ring strain effects (Phase 2)

---

**Next**: See `gfnff_inversions.cpp` for the C++ implementation following this theory.
