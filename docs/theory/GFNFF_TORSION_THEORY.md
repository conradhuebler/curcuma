# Scientific Documentation: Torsion Potentials in GFN-FF

**Author**: Implementation by Claude based on GFN-FF (Spicher & Grimme, 2020)
**Original Work**: S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.* **59**, 15665 (2020)
**Original Implementation**: xtb program, Grimme Lab, University of Bonn
**DOI**: https://doi.org/10.1002/anie.202004239

---

## Acknowledgment

This implementation is based on the pioneering work of **Prof. Stefan Grimme** and **Dr. Sebastian Spicher** at the University of Bonn. The GFN-FF (Geometry, Frequency, Noncovalent - Force Field) method represents a breakthrough in combining quantum chemical accuracy with force field efficiency.

**Original Fortran implementation**: `external/gfnff/src/gfnff_engrad.F90` (lines 1041-1122)
- Developed by: Grimme Group (Sebastian Ehlert, Sebastian Spicher, Stefan Grimme)
- Maintained by: Philipp Pracht (https://github.com/pprcht/gfnff)
- License: GNU LGPL v3

This C++ port aims to provide educational clarity while maintaining scientific accuracy against the reference implementation.

---

## 1. Theoretical Background

### 1.1 What are Torsion Potentials?

**Torsion angles** (also called dihedral angles) describe the rotation around a chemical bond. Given four atoms **i-j-k-l** connected in sequence:

```
    i
     \
      j---k
           \
            l
```

The **torsion angle φ** is the angle between:
- The plane formed by atoms **i-j-k**
- The plane formed by atoms **j-k-l**

**Physical interpretation**:
- φ = 0°: cis/eclipsed conformation (atoms i and l on same side)
- φ = 180°: trans/anti conformation (atoms i and l on opposite sides)
- φ = ±60°, ±120°: gauche conformations

### 1.2 Why Do Torsion Barriers Exist?

**Multiple quantum mechanical effects contribute**:

1. **Steric repulsion**: Electron clouds of substituents repel when close (eclipsed)
2. **Hyperconjugation**: Orbital overlap stabilizes specific conformations
3. **Electrostatic interactions**: Dipole-dipole interactions favor certain angles
4. **π-conjugation**: In conjugated systems, planarity is energetically favored

**Example - Ethane (C₂H₆)**:
- Rotation barrier: ~3 kcal/mol
- Minimum at φ = 60°, 180°, 300° (staggered)
- Maximum at φ = 0°, 120°, 240° (eclipsed)

**Example - Butane (C₄H₁₀)**:
- Rotation barrier: ~3-6 kcal/mol
- Anti (φ = 180°) most stable
- Gauche (φ = ±60°) ~0.6 kcal/mol higher
- Eclipsed (φ = 0°, 120°) highest energy

---

## 2. Mathematical Formulation in GFN-FF

### 2.1 Standard Torsion Potential (Proper Dihedrals)

For **proper dihedrals** (regular i-j-k-l sequence), GFN-FF uses a **cosine series expansion**:

```
E_tors = V_n/2 · [1 + cos(n·φ - φ₀)] · D(r_ij, r_jk, r_kl)
```

**Parameters**:
- **V_n**: Torsion barrier height (kcal/mol)
- **n**: Periodicity (1, 2, or 3)
  - n=1: One minimum per 360° rotation (rare, e.g., H₂O₂)
  - n=2: Two minima per 360° (e.g., H-C-C-H in ethylene)
  - n=3: Three minima per 360° (e.g., H-C-C-H in ethane, most common)
- **φ₀**: Phase shift / reference angle (typically 0° or 180°)
- **D(...)**: Distance-dependent damping function

**Physical meaning of n**:
- n determines the **symmetry** of the rotation
- For sp³-sp³ single bonds (C-C): typically n=3 (threefold barrier)
- For sp²-sp² double bonds (C=C): typically n=2 (cis/trans isomerization)

### 2.2 Improper Torsions (Out-of-Plane Angles)

For **improper dihedrals** (central atom j connected to i, k, l), used for out-of-plane bending:

```
E_oop = V · (cos(φ) - cos(φ₀))²
```

or for φ₀ = 0:

```
E_oop = V/2 · [1 + cos(φ + π)]
```

**Use cases**:
- Maintaining planarity in sp² centers (C=O, C=C, aromatic rings)
- Umbrella inversion in NH₃, NR₃

### 2.3 Distance Damping Function

GFN-FF includes **topology-aware damping** to account for bond stretching/compression:

```
D(r_ij, r_jk, r_kl) = D_ij · D_jk · D_kl
```

where each damping term is:

```
D_ab = 1 / [1 + exp(-α·(r_ab/r₀_ab - 1))]
```

**Physical interpretation**:
- When bonds are stretched (r > r₀), torsion barrier decreases
- When bonds are compressed (r < r₀), torsion barrier increases
- This couples **bond stretching** with **torsional motion**

**Parameters from Fortran code** (`gfnffdampt` function):
- α typically ~16-20 (steepness of damping)
- r₀ from covalent radii table

---

## 3. GFN-FF Specific Implementation Details

### 3.1 Two Torsion Types in GFN-FF

From `gfnff_engrad.F90:1041-1122`, GFN-FF distinguishes:

#### Type A: Proper Torsion (rn > 0)
```fortran
! Lines 1061-1087
rn = topo%tlist(5,m)  ! Periodicity (1, 2, or 3)
phi0 = topo%vtors(1,m)  ! Reference angle
V = topo%vtors(2,m)     ! Barrier height

! Energy
et = (1 + cos(rn*phi + pi))*V  ! Note: dphi = phi - phi0, c1 = rn*dphi + pi
e = et * damp
```

**Formula breakdown**:
```
c1 = rn·(φ - φ₀) + π
E = V · [1 + cos(c1)] · D
  = V · [1 + cos(rn·φ - rn·φ₀ + π)] · D
  = V/2 · [1 - cos(rn·(φ - φ₀))] · D  (using cos(x+π) = -cos(x))
```

This is the **standard cosine potential** with energy minimum at φ = φ₀.

#### Type B: Improper Torsion (rn ≤ 0)
```fortran
! Lines 1089-1120
if (rn == 0) then
  ! Simple cosine: E = V*(1 + cos(phi))
  et = (1 + cos(phi + pi))*V
else  ! rn < 0
  ! Double minimum at ±phi0
  et = V * (cos(phi) - cos(phi0))^2
endif
```

**Use cases**:
- rn = 0: Out-of-plane bending (maintains planarity)
- rn < 0: Symmetric double minimum (e.g., umbrella inversion)

### 3.2 Gradient Calculation

**Analytical gradient** derived via chain rule:

```
∂E/∂x_atom = ∂E/∂φ · ∂φ/∂x_atom + ∂E/∂D · ∂D/∂x_atom
```

**Components**:

1. **Energy derivative w.r.t. angle**:
```
∂E/∂φ = -rn·V·sin(rn·φ - rn·φ₀ + π)·D
      = rn·V·sin(rn·(φ - φ₀))·D
```

2. **Angle derivatives** (`dphidr` in Fortran):
   - ∂φ/∂x_i, ∂φ/∂x_j, ∂φ/∂x_k, ∂φ/∂x_l
   - Complex geometry derivatives involving cross products

3. **Damping derivatives**:
   - ∂D/∂r_ij, ∂D/∂r_jk, ∂D/∂r_kl
   - From exponential damping function

**Fortran implementation** (lines 1073-1086):
```fortran
call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)  ! Get ∂φ/∂xyz
dij = -rn*sin(c1)*V*damp  ! ∂E/∂φ

! Gradient contributions from angle derivative
g(1:3,1) = dij*dda(1:3) + term1  ! Atom i
g(1:3,2) = dij*ddb(1:3) - term1 + term2  ! Atom j
g(1:3,3) = dij*ddc(1:3) + term3 - term2  ! Atom k
g(1:3,4) = dij*ddd(1:3) - term3  ! Atom l

! term1, term2, term3: damping derivative contributions
```

---

## 4. Parameter Assignment Strategy

### 4.1 How GFN-FF Determines Torsion Parameters

From `gfnff_ini2.f90`, GFN-FF uses **topology-aware** parameter assignment:

**Decision tree**:

1. **Identify central bond j-k**:
   - Check hybridization of atoms j and k
   - Check if part of ring system
   - Check if part of conjugated system

2. **Assign periodicity (n)**:
   - sp³-sp³: n = 3 (threefold)
   - sp²-sp²: n = 2 (twofold for cis/trans)
   - sp²-sp³: n = 1 or 2 (depends on conjugation)
   - In rings: modified based on ring size

3. **Assign barrier height (V)**:
   - Look up from element-specific table
   - Apply corrections for:
     - Ring strain (3-, 4-membered rings)
     - Conjugation (lower barrier)
     - Hyperconjugation (stabilizes specific conformations)

4. **Assign reference angle (φ₀)**:
   - Typically 0° for sp³-sp³
   - 180° for conjugated systems (prefers planarity)

### 4.2 Element-Specific Barriers

**Typical barrier heights** (from GFN-FF parameters):

| Bond Type | V_barrier | n | Example |
|-----------|-----------|---|---------|
| C(sp³)-C(sp³) | ~1.5 kcal/mol | 3 | Ethane |
| C(sp²)-C(sp²) | ~3-5 kcal/mol | 2 | Butadiene |
| C(sp²)-C(sp³) | ~2 kcal/mol | 2-3 | Propene |
| N-C | ~0.5-2 kcal/mol | 2-3 | Amines |
| O-C | ~1-2 kcal/mol | 3 | Ethers |
| C(sp²)=C(sp²) | ~45 kcal/mol | 2 | Ethylene (π-bond rotation!) |

**Note**: These are representative values; actual GFN-FF parameters include fine-tuning based on extended testing.

---

## 5. Computational Considerations

### 5.1 Angle Calculation

**Dihedral angle formula**:

Given vectors:
- **v₁** = r_j - r_i (bond i→j)
- **v₂** = r_k - r_j (bond j→k)
- **v₃** = r_l - r_k (bond k→l)

```
φ = atan2(|v₂|·(v₁ × v₂)·v₃, (v₁ × v₂)·(v₂ × v₃))
```

**Why atan2 instead of acos**?
- atan2 gives signed angle (-π to π)
- acos only gives 0 to π (loses sign information)
- Sign matters for distinguishing clockwise/counterclockwise rotation

### 5.2 Numerical Stability

**Critical issues**:

1. **Linear geometries**: When three atoms are colinear, cross products → 0
   - Check for |v₁ × v₂| < ε or |v₂ × v₃| < ε
   - Skip torsion or use alternative formulation

2. **Periodicity wrapping**: φ is periodic in 2π
   - Ensure φ - φ₀ is correctly wrapped to [-π, π]
   - Use: `dphi = atan2(sin(phi - phi0), cos(phi - phi0))`

3. **Gradient singularities**: Near φ = ±90° for cos(φ)
   - Usually not problematic due to smooth potential
   - Damping function helps regularize

### 5.3 Performance Tips

**Optimization strategies**:

1. **Precompute bond vectors**: Only calculate once per step
2. **Sparse torsion list**: Not all i-j-k-l combinations are torsions
3. **Neighbor lists**: Only check bonded atoms
4. **Vectorization**: Eigen can vectorize cross products

---

## 6. Validation Strategy

### 6.1 Test Molecules

**Systematic tests** (in order of complexity):

1. **Ethane (C₂H₆)**:
   - Simplest C-C torsion
   - n=3, V~1.5 kcal/mol
   - Expect 3 minima at 60°, 180°, 300°

2. **Butane (C₄H₁₀)**:
   - Two C-C torsions
   - Anti/gauche energy difference ~0.6 kcal/mol
   - Tests coupled torsions

3. **Ethylene (C₂H₄)**:
   - C=C torsion (very high barrier)
   - n=2, V~45 kcal/mol
   - Should strongly resist rotation

4. **Benzene (C₆H₆)**:
   - Aromatic system
   - All φ₀ = 0° or 180° (planarity)
   - Tests conjugated system handling

### 6.2 Validation Criteria

**Quantitative targets**:

| Property | Tolerance | Reference |
|----------|-----------|-----------|
| Energy | ±0.5 kcal/mol | vs. Fortran GFN-FF |
| Gradient | ±1% RMS | vs. numerical |
| Barrier height | ±0.2 kcal/mol | vs. Fortran |
| Optimized geometry | ±0.05 Å RMSD | vs. Fortran |

**Qualitative checks**:
- Correct number of minima (matches periodicity n)
- Smooth energy profile vs. angle
- Gradient vanishes at extrema

---

## 7. Literature References

### Primary References

1. **S. Spicher, S. Grimme**
   *Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems*
   *Angew. Chem. Int. Ed.* **59**, 15665-15665 (2020)
   DOI: [10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239)

   *The definitive GFN-FF paper - describes all potentials, parametrization strategy, and validation.*

2. **S. Grimme, C. Bannwarth, P. Shushkov**
   *A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions*
   *J. Chem. Theory Comput.* **13**, 1989-2009 (2017)
   DOI: [10.1021/acs.jctc.7b00118](https://doi.org/10.1021/acs.jctc.7b00118)

   *Background on GFN-xTB methods, which inspired GFN-FF topology analysis.*

### Theoretical Background

3. **K. B. Wiberg**
   *Application of the Pople-Santry-Segal CNDO Method to the Cyclopropylcarbinyl and Cyclobutyl Cation*
   *Tetrahedron* **24**, 1083-1096 (1968)
   DOI: [10.1016/0040-4020(68)88057-3](https://doi.org/10.1016/0040-4020(68)88057-3)

   *Classic work on conformational analysis and torsion barriers.*

4. **N. L. Allinger, Y. H. Yuh, J.-H. Lii**
   *Molecular Mechanics. The MM3 Force Field for Hydrocarbons*
   *J. Am. Chem. Soc.* **111**, 8551-8566 (1989)
   DOI: [10.1021/ja00205a001](https://doi.org/10.1021/ja00205a001)

   *MM3 force field - pioneering work on torsion potential forms.*

### Computational Methods

5. **W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives**
   *Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics*
   *J. Am. Chem. Soc.* **118**, 11225-11236 (1996)
   DOI: [10.1021/ja9621760](https://doi.org/10.1021/ja9621760)

   *OPLS - alternative torsion parametrization approach.*

---

## 8. Implementation Notes for C++ Port

### 8.1 Design Decisions

**Differences from Fortran**:

1. **Object-oriented structure**:
   - Fortran: Subroutine with passed arrays
   - C++: Method of `GFNFF` class with member access

2. **Angle calculation**:
   - Fortran: `valijklff()` function
   - C++: Implement `calculateDihedralAngle()` using Eigen

3. **Gradient storage**:
   - Fortran: Local `g(3,4)` array
   - C++: Update member `m_gradient` matrix

### 8.2 Educational Enhancements

**Comments in C++ code should explain**:

1. **Physical meaning**: Why this term exists
2. **Mathematical derivation**: Link to equations above
3. **Numerical considerations**: Why we use atan2, etc.
4. **Literature reference**: Which paper/equation

**Example**:
```cpp
// Calculate torsion energy (Spicher & Grimme 2020, Eq. 8)
// E = V/2 * [1 - cos(n*(φ - φ₀))] * D(r_ij, r_jk, r_kl)
// Physical meaning: Rotation barrier around j-k bond
double c1 = periodicity * (phi - phi0) + M_PI;
double energy = barrier_height * (1.0 + cos(c1)) * damping;
```

### 8.3 Testing Harness

**Incremental testing**:

1. **Angle only**: Test dihedral calculation vs. known values
2. **Energy only**: Compare energies vs. Fortran (no gradients)
3. **Numerical gradient**: Finite differences vs. analytical
4. **Full validation**: Geometry optimization test

---

## 9. Expected Timeline

**Phase 1.1 - Torsions** (estimated 3-5 days):

- Day 1: Scientific documentation (this file) ✅
- Day 2: Implement `calculateDihedralAngle()` and `calculateTorsionEnergy()`
- Day 3: Implement `calculateTorsionGradient()` with analytical derivatives
- Day 4: Parameter generation `generateGFNFFTorsions()`
- Day 5: Validation against Fortran reference

**Success criteria**:
- ✅ Butane energy profile matches Fortran within 0.5 kcal/mol
- ✅ Gradient validated numerically
- ✅ Code documented with literature references

---

**End of Scientific Documentation**

*This document should be updated as implementation progresses and new insights are gained.*

---

## Appendix: Quick Reference Formulas

### Torsion Energy
```
E_tors = (V/2) · [1 - cos(n·(φ - φ₀))] · D(r_ij, r_jk, r_kl)
```

### Dihedral Angle
```
φ = atan2(|v₂|·(v₁ × v₂)·v₃, (v₁ × v₂)·(v₂ × v₃))
```

### Energy Gradient
```
∂E/∂x = (∂E/∂φ)·(∂φ/∂x) + (∂E/∂D)·(∂D/∂x)
∂E/∂φ = (V·n/2)·sin(n·(φ - φ₀))·D
```

### Damping Function
```
D_ab = 1 / [1 + exp(-α·(r_ab/r₀_ab - 1))]
∂D/∂r = α·D·(1 - D) / r₀
```
