# Collective Variables (CVs) in Curcuma

**Author:** Claude (Anthropic) & Conrad Hübler
**Date:** November 2025
**Status:** Implementation in Progress

---

## Table of Contents

1. [Introduction](#introduction)
2. [Theoretical Background](#theoretical-background)
3. [Architecture Overview](#architecture-overview)
4. [Implemented CV Types](#implemented-cv-types)
5. [Mathematical Formulations](#mathematical-formulations)
6. [Usage Examples](#usage-examples)
7. [References](#references)

---

## 1. Introduction

### What are Collective Variables?

**Collective Variables (CVs)** are low-dimensional descriptors of complex molecular configurations used in enhanced sampling methods like metadynamics. Instead of tracking all 3N atomic coordinates, CVs reduce the system to a few meaningful degrees of freedom (e.g., distance between two atoms, dihedral angle, radius of gyration).

### Why Multiple CV Types?

Different chemical processes require different CVs:
- **Protein folding:** Radius of gyration, backbone dihedrals
- **Ligand binding:** Distances, coordination numbers
- **Conformational changes:** RMSD, angles
- **Chemical reactions:** Bond distances, coordination

**Key Principle:** The choice of CV determines which regions of phase space are explored during metadynamics [1].

---

## 2. Theoretical Background

### 2.1 Enhanced Sampling with CVs

In metadynamics, a bias potential `V(s, t)` is constructed as a function of collective variables `s(r)`, where `r` are atomic coordinates:

```
V(s, t) = Σ_i w_i * exp(-(s - s_i)² / (2σ²))
```

**Key Papers:**
- **Laio & Parrinello (2002):** Original metadynamics paper [1]
- **Barducci et al. (2008):** Well-tempered metadynamics [2]
- **Branduardi et al. (2007):** Path collective variables [3]

### 2.2 Requirements for Good CVs

A good collective variable must satisfy:

1. **Differentiability:** `∇s(r)` must be well-defined (for force calculation)
2. **Low-dimensionality:** Typically 1-3 CVs (curse of dimensionality)
3. **Physical relevance:** Captures the slow degrees of freedom
4. **Computationally efficient:** Fast evaluation during MD

**Mathematical Gradient:**
```
F_bias = -∇_r V(s, t) = -(∂V/∂s) * (∂s/∂r)
```

Where:
- `∂V/∂s`: Derivative of bias potential w.r.t. CV
- `∂s/∂r`: Gradient of CV w.r.t. atomic positions (Jacobian)

---

## 3. Architecture Overview

### 3.1 Design Philosophy

**Inspired by PLUMED's architecture** [4], we implement a polymorphic CV system:

```
CollectiveVariable (abstract base class)
    ├─ CV_RMSD
    ├─ CV_Distance
    ├─ CV_Angle
    ├─ CV_Dihedral
    ├─ CV_RadiusOfGyration
    ├─ CV_Coordination
    └─ CV_PathCV (future)
```

**Key Design Decisions:**

1. **Polymorphism:** All CVs inherit from `CollectiveVariable` base class
2. **Eigen Integration:** All gradients returned as `Eigen::MatrixXd` (N_atoms × 3)
3. **Lazy Evaluation:** CVs only computed when needed
4. **SIMD-friendly:** Vectorized calculations where possible
5. **Unit-aware:** All CVs return values in standard units (Å, degrees, etc.)

### 3.2 Class Hierarchy

```cpp
/**
 * @brief Abstract base class for all collective variables
 *
 * All derived CV classes must implement:
 * - calculate(): Compute CV value from current geometry
 * - gradient(): Compute ∂s/∂r (Jacobian matrix)
 * - type(): Return CV type identifier
 *
 * References:
 * [1] Bonomi et al., Comp. Phys. Comm. 180, 1961 (2009) - PLUMED design
 */
class CollectiveVariable {
public:
    virtual double calculate(const Molecule& mol) = 0;
    virtual Geometry gradient(const Molecule& mol) = 0;
    virtual std::string type() const = 0;
    virtual std::string description() const = 0;

    // Configuration
    void setAtoms(const std::vector<int>& atoms);
    const std::vector<int>& getAtoms() const { return m_atoms; }

protected:
    std::vector<int> m_atoms;  // Atom indices involved in this CV
    std::string m_name;         // User-defined label
};
```

---

## 4. Implemented CV Types

### 4.1 RMSD (Root-Mean-Square Deviation)

**Physical Meaning:** Measures structural similarity to a reference structure.

**Formula:**
```
RMSD(r, r_ref) = sqrt(1/N * Σ_i |r_i - r_ref,i|²)
```

**Gradient (Kabsch-aligned):**
```
∂RMSD/∂r_i = (r_i - r_ref,i) / (N * RMSD)
```

**Use Cases:**
- Protein conformational changes
- Ligand docking
- Molecular recognition

**References:**
- Kabsch (1976): Optimal rotation algorithm [5]
- Grimme (2019): RMSD-MTD for conformer sampling [6]

**Implementation Notes:**
- Uses SVD-based Kabsch algorithm for optimal alignment
- Mass-weighted option available
- Fragment-aware: Can compute inter-fragment distances

---

### 4.2 Distance

**Physical Meaning:** Distance between two atoms or centers of mass.

**Formula:**
```
d(r) = |r_j - r_i|
```

**Gradient:**
```
∂d/∂r_i = -(r_j - r_i) / d
∂d/∂r_j = +(r_j - r_i) / d
```

**Use Cases:**
- Bond breaking/formation
- Ligand-protein distance
- Reaction coordinates

**Implementation Notes:**
- Periodic boundary conditions (PBC) aware
- Can use center-of-mass for groups of atoms
- Numerically stable for d → 0

---

### 4.3 Angle

**Physical Meaning:** Angle formed by three atoms (i-j-k), vertex at j.

**Formula:**
```
θ(r) = arccos((r_ij · r_kj) / (|r_ij| * |r_kj|))
```
where `r_ij = r_i - r_j`, `r_kj = r_k - r_j`

**Gradient (chain rule):**
```
∂θ/∂r_i = (1 / sin(θ)) * [(cos(θ) * r_ij / |r_ij|² ) - (r_kj / (|r_ij| * |r_kj|))]
∂θ/∂r_k = (1 / sin(θ)) * [(cos(θ) * r_kj / |r_kj|² ) - (r_ij / (|r_ij| * |r_kj|))]
∂θ/∂r_j = -(∂θ/∂r_i + ∂θ/∂r_k)
```

**Use Cases:**
- Valence angle deformations
- Hydrogen bonding geometry
- Ring puckering

**References:**
- Tribello et al. (2014): PLUMED 2 documentation [7]

**Implementation Notes:**
- Handles θ → 0° and θ → 180° singularities
- Returns angle in degrees (user-friendly) or radians (internal)

---

### 4.4 Dihedral (Torsion Angle)

**Physical Meaning:** Dihedral angle formed by four atoms (i-j-k-l).

**Formula:**
```
φ(r) = atan2((r_ij × r_jk) · (r_jk × r_kl) / |r_jk|,
              (r_ij × r_jk) · (r_kl))
```

**Gradient:**
Uses the cross-product chain rule. See implementation for full derivation.

**Use Cases:**
- Protein backbone angles (φ, ψ)
- Rotation around single bonds
- Conformer interconversion

**References:**
- Blondel & Karplus (1996): Derivation of analytical gradients [8]

**Implementation Notes:**
- Periodicity handling: φ ∈ [-180°, +180°]
- Efficient cross-product calculation
- PBC-aware

---

### 4.5 Radius of Gyration

**Physical Meaning:** Root-mean-square distance of atoms from their center of mass.

**Formula:**
```
R_g = sqrt(Σ_i m_i * |r_i - r_COM|² / Σ_i m_i)
```

where `r_COM = Σ_i m_i * r_i / Σ_i m_i`

**Gradient:**
```
∂R_g/∂r_i = (m_i / (M * R_g)) * (r_i - r_COM) * (1 - m_i/M)
```

**Use Cases:**
- Protein folding/unfolding
- Polymer collapse
- Cluster formation

**References:**
- Dill & MacCallum (2012): Protein folding review [9]

**Implementation Notes:**
- Mass-weighted by default
- Can use geometric center (unweighted)

---

### 4.6 Coordination Number

**Physical Meaning:** Number of atoms of type B within cutoff distance of atoms of type A.

**Formula (smooth switching function):**
```
CN = Σ_{i∈A} Σ_{j∈B} s(r_ij)

s(r) = (1 - (r/r_0)^n) / (1 - (r/r_0)^m)
```

where `n=6`, `m=12` (typical), `r_0` is the cutoff.

**Gradient:**
```
∂CN/∂r_i = Σ_j (∂s/∂r_ij) * (∂r_ij/∂r_i)
```

**Use Cases:**
- Solvation shells
- Ligand binding sites
- Metal coordination

**References:**
- Iannuzzi et al. (2003): Smooth coordination function [10]

**Implementation Notes:**
- Differentiable (smooth) switching function
- Rational switching function prevents discontinuities
- PBC-aware

---

## 5. Mathematical Formulations

### 5.1 Chain Rule for CV Gradients

The bias force on atom `i` due to CV `s` is:

```
F_i^bias = -(∂V/∂s) * (∂s/∂r_i)
```

**Key Observation:** We need two components:
1. `∂V/∂s`: From metadynamics algorithm (sum of Gaussians)
2. `∂s/∂r_i`: From CV implementation (Jacobian)

### 5.2 Numerical Stability

**Problem:** Some CVs have singularities (e.g., angle gradient at θ=0°).

**Solutions implemented:**
- **Distance:** Add `ε = 10^-8 Å` to denominator when d < 10^-6
- **Angle:** Switch to alternative formula near θ=0° and θ=180°
- **Dihedral:** Use `atan2()` for proper quadrant handling

### 5.3 Performance Optimization

**SIMD Vectorization:** Where possible, use Eigen's vectorized operations.

Example (Distance CV):
```cpp
// Vectorized (fast)
Eigen::Vector3d r_ij = mol.getPosition(j) - mol.getPosition(i);
double d = r_ij.norm();

// vs. component-wise (slow)
double dx = x[j] - x[i];
double dy = y[j] - y[i];
double dz = z[j] - z[i];
double d = sqrt(dx*dx + dy*dy + dz*dz);
```

---

## 6. Usage Examples

### 6.1 Single CV: Distance

```bash
./curcuma -simplemd input.xyz \
  -rmsd_mtd true \
  -cv_type distance \
  -cv_atoms "0,10" \
  -rmsd_mtd_k 0.1 \
  -rmsd_mtd_alpha 10.0
```

### 6.2 Multiple CVs: Distance + Angle

```bash
./curcuma -simplemd input.xyz \
  -rmsd_mtd true \
  -cv1_type distance -cv1_atoms "0,10" \
  -cv2_type angle -cv2_atoms "5,0,10" \
  -rmsd_mtd_k 0.1 \
  -rmsd_mtd_alpha 10.0
```

### 6.3 Protein Folding: Radius of Gyration + RMSD

```bash
./curcuma -simplemd protein.xyz \
  -rmsd_mtd true \
  -cv1_type gyration -cv1_atoms "backbone" \
  -cv2_type rmsd -cv2_atoms "all" \
  -rmsd_mtd_k 0.05 \
  -rmsd_mtd_alpha 5.0
```

---

## 7. References

[1] **Laio, A. & Parrinello, M.** (2002). *Escaping free-energy minima.* Proc. Natl. Acad. Sci. USA **99**, 12562-12566. DOI: [10.1073/pnas.202427399](https://doi.org/10.1073/pnas.202427399)

[2] **Barducci, A., Bussi, G. & Parrinello, M.** (2008). *Well-tempered metadynamics: A smoothly converging and tunable free-energy method.* Phys. Rev. Lett. **100**, 020603. DOI: [10.1103/PhysRevLett.100.020603](https://doi.org/10.1103/PhysRevLett.100.020603)

[3] **Branduardi, D., Gervasio, F. L. & Parrinello, M.** (2007). *From A to B in free energy space.* J. Chem. Phys. **126**, 054103. DOI: [10.1063/1.2432340](https://doi.org/10.1063/1.2432340)

[4] **Bonomi, M. et al.** (2009). *PLUMED: A portable plugin for free-energy calculations with molecular dynamics.* Comp. Phys. Comm. **180**, 1961-1972. DOI: [10.1016/j.cpc.2009.05.011](https://doi.org/10.1016/j.cpc.2009.05.011)

[5] **Kabsch, W.** (1976). *A solution for the best rotation to relate two sets of vectors.* Acta Cryst. A **32**, 922-923. DOI: [10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873)

[6] **Grimme, S.** (2019). *Exploration of Chemical Compound, Conformer, and Reaction Space with Meta-Dynamics Simulations.* J. Chem. Theory Comput. **15**, 2847-2862. DOI: [10.1021/acs.jctc.9b00143](https://doi.org/10.1021/acs.jctc.9b00143)

[7] **Tribello, G. A. et al.** (2014). *PLUMED 2: New feathers for an old bird.* Comp. Phys. Comm. **185**, 604-613. DOI: [10.1016/j.cpc.2013.09.018](https://doi.org/10.1016/j.cpc.2013.09.018)

[8] **Blondel, A. & Karplus, M.** (1996). *New formulation for derivatives of torsion angles and improper torsion angles in molecular mechanics.* J. Comp. Chem. **17**, 1132-1141. DOI: [10.1002/(SICI)1096-987X(19960715)17:9<1132::AID-JCC5>3.0.CO;2-T](https://doi.org/10.1002/(SICI)1096-987X(19960715)17:9<1132::AID-JCC5>3.0.CO;2-T)

[9] **Dill, K. A. & MacCallum, J. L.** (2012). *The protein-folding problem, 50 years on.* Science **338**, 1042-1046. DOI: [10.1126/science.1219021](https://doi.org/10.1126/science.1219021)

[10] **Iannuzzi, M., Laio, A. & Parrinello, M.** (2003). *Efficient exploration of reactive potential energy surfaces using Car-Parrinello molecular dynamics.* Phys. Rev. Lett. **90**, 238302. DOI: [10.1103/PhysRevLett.90.238302](https://doi.org/10.1103/PhysRevLett.90.238302)

---

## Appendix A: Comparison with PLUMED

| Feature | PLUMED | Curcuma (Planned) |
|---------|--------|-------------------|
| Number of CV types | ~200 | ~10 (core set) |
| Distance | ✅ | ✅ |
| Angle | ✅ | ✅ |
| Dihedral | ✅ | ✅ |
| RMSD | ✅ | ✅ |
| Radius of Gyration | ✅ | ✅ |
| Coordination | ✅ | ✅ |
| Path CVs | ✅ | ❌ (future) |
| Custom functions | ✅ | ❌ |
| PBC support | ✅ | ✅ |
| Grid acceleration | ✅ | 🔄 (implementing) |
| FES analysis | ✅ | 🔄 (implementing) |

**Philosophy:** Curcuma focuses on a **curated set of essential CVs** with excellent documentation, rather than exhaustive coverage.

---

## Appendix B: Testing Strategy

Each CV implementation includes:

1. **Unit tests:** Verify gradient correctness via finite differences
2. **Numerical stability tests:** Boundary cases (d→0, θ→0°, etc.)
3. **Performance benchmarks:** Compare with analytical solutions
4. **Integration tests:** Full MD runs with known free energy surfaces

Example test:
```cpp
TEST(CV_Distance, GradientCorrectness) {
    // Compare analytical gradient with numerical finite difference
    double epsilon = 1e-6;
    Geometry grad_analytical = cv.gradient(mol);
    Geometry grad_numerical = numericalGradient(cv, mol, epsilon);
    EXPECT_NEAR(grad_analytical, grad_numerical, 1e-4);
}
```

---

**End of Document**
