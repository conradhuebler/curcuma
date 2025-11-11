# Native QM Methods Implementation Guide

**Document Purpose**: Comprehensive implementation guide for native GFN1, GFN2, and PMx semi-empirical quantum mechanical methods in Curcuma.

**Created**: 2025-11-11
**Based on**: TBLite source code analysis, MOPAC repository structure, and existing Curcuma EHT implementation
**Target Audience**: Theoretical chemists with broad background knowledge who need quick access to implementation-level depth

---

## Important Notes for Implementation

### Attribution and Licensing

**CRITICAL**: All original method authors must be properly attributed in the implementation.

**Original Method Developers**:
- **GFN2-xTB**: Stefan Grimme, Christoph Bannwarth, Sebastian Ehlert (University of Bonn)
- **GFN1-xTB**: Stefan Grimme, Christoph Bannwarth, Philip Shushkov (University of Bonn)
- **TBLite Implementation**: Sebastian Ehlert and contributors (LGPL-3.0 License)
- **PM3**: James J. P. Stewart (1989)
- **AM1**: Michael J. S. Dewar et al. (1985)
- **PM6**: James J. P. Stewart (2007)
- **MOPAC**: Original code by James J. P. Stewart, now open-source (LGPL-3.0)

**Copyright and License Requirements**:
```cpp
/*
 * <Native GFN2-xTB Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the GFN2-xTB method developed by:
 *   Stefan Grimme, Christoph Bannwarth, Sebastian Ehlert
 *   Mulliken Center for Theoretical Chemistry, University of Bonn
 *
 * Reference implementation: TBLite (https://github.com/tblite/tblite)
 *   Copyright (C) 2019-2024 Sebastian Ehlert and contributors
 *   Licensed under LGPL-3.0-or-later
 *
 * Original method publication:
 *   C. Bannwarth, S. Ehlert, S. Grimme
 *   J. Chem. Theory Comput. 2019, 15, 1652-1671
 *   DOI: 10.1021/acs.jctc.8b01176
 *
 * This implementation is an independent C++ port for educational purposes
 * within the Curcuma framework, maintaining compatibility with the original
 * method while following Curcuma's educational-first design principles.
 */
```

**Parameter Attribution**:
All parameter tables (element data, orbital exponents, coupling constants) must cite:
1. **Primary source**: Original method publication
2. **Implementation source**: TBLite parameter files (commit hash if applicable)
3. **Any modifications**: Clearly documented with scientific justification

### Theoretical Chemistry Accessibility

**Design Principle**: This implementation prioritizes transparency for theoretical chemists with broad knowledge who need to quickly understand implementation details without deep software engineering expertise.

**Documentation Standards**:

1. **Every mathematical formula** must include:
   - Physical meaning of each variable
   - Units (eV, Hartree, Bohr, Ångström)
   - Literature reference (equation number from original paper if possible)
   - Numerical example for validation

2. **Every parameter** must document:
   - Physical interpretation (e.g., "Hubbard parameter U: on-site Coulomb repulsion")
   - Typical value range (e.g., "0.1-0.5 Hartree for 2nd-period elements")
   - Original source (paper + table/equation number)
   - How it was determined (fitted, calculated, empirical)

3. **Every algorithm step** must include:
   - Theoretical background (which approximation/theory)
   - Why this approach (numerical stability, physical justification)
   - Alternative formulations (if any)
   - Literature reference

**Example of Good Documentation**:
```cpp
/**
 * @brief Calculate exponential coordination number (ECP)
 *
 * Theoretical Background:
 *   The coordination number CN_i represents the number of atoms bonded to atom i.
 *   The exponential counting function provides a continuous, differentiable measure
 *   suitable for gradient-based geometry optimization.
 *
 * Formula:
 *   CN_i = ∑_{j≠i} [1 / (1 + exp(-k₁(R_cov,ij/R_ij - 1)))]^(k₂)
 *
 * Parameters:
 *   k₁ = 16.0      Steepness parameter (controls transition sharpness)
 *   k₂ = 4/3       Range parameter (controls long-range decay)
 *   R_cov,ij       Sum of covalent radii (Å) from Pyykkö, J. Phys. Chem. A 2015
 *   R_ij           Interatomic distance (Å)
 *
 * Physical Interpretation:
 *   - CN ≈ 0 for R_ij >> R_cov (no bonding interaction)
 *   - CN ≈ 1 for R_ij ≈ R_cov (typical single bond)
 *   - CN smooth transition allows analytical gradients
 *
 * Reference:
 *   S. Grimme et al., J. Chem. Phys. 2010, 132, 154104 (D3 paper, Eq. 9)
 *   C. Bannwarth et al., J. Chem. Theory Comput. 2019, 15, 1652 (GFN2, Eq. 4)
 *
 * Implementation Notes:
 *   - Uses std::pow(count, 4.0/3.0) instead of count^(4/3) for clarity
 *   - Covalent radii stored in params.covalent_radius[] (in Ångström)
 *   - Result used for Hamiltonian diagonal shift: E_ii += k_CN * CN_i
 *
 * @param atom Atomic number of central atom
 * @param positions 3D coordinates of all atoms (Nx3 matrix, Ångström)
 * @return Coordination number CN (dimensionless, typically 0-12)
 */
double calculateCoordinationNumber(int atom, const Matrix& positions);
```

### External Dependencies: D3/D4 Dispersion

**IMPORTANT**: D3 and D4 dispersion corrections are **NOT** part of this native implementation.

**Current Status**:
- Curcuma already has `dftd3interface.h/cpp` and `dftd4interface.h/cpp`
- These interfaces are **separate TODOs** and will be fixed independently
- Native GFN1/GFN2/PMx implementations should use **stub functions** initially

**Stub Implementation Strategy**:
```cpp
double GFN2::calculateDispersionEnergy() const {
    // TODO: D4 dispersion integration (separate task)
    // For now, return zero to allow testing of other components

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::warn("D4 dispersion not yet integrated - energy incomplete");
    }

    return 0.0;  // Stub: no dispersion contribution
}

Matrix GFN2::calculateDispersionGradient() const {
    // TODO: D4 dispersion gradient (separate task)
    return Matrix::Zero(m_natoms, 3);  // Stub: no dispersion gradient
}
```

**Future Integration** (when D3/D4 are fixed):
```cpp
double GFN2::calculateDispersionEnergy() const {
    // GFN2 uses DFT-D4 with specific parameters
    D4Interface d4;
    d4.setParameters(1.0, 2.7, 0.52, 5.0, 5.0);  // s6, s8, a1, a2, s9
    d4.InitialiseMolecule(m_molecule);
    return d4.Calculation();
}
```

**Rationale**:
- D3/D4 are complex standalone libraries (separate codebases)
- Native GFN implementation focuses on tight-binding Hamiltonian
- Dispersion can be validated separately with known reference systems
- Allows parallel development: tight-binding core + dispersion fixes

---

## Table of Contents

1. [Theoretical Background and Literature](#theoretical-background-and-literature)
2. [GFN2-xTB Implementation](#gfn2-xtb-implementation)
3. [GFN1-xTB Implementation](#gfn1-xtb-implementation)
4. [PMx Methods Implementation](#pmx-methods-implementation)
5. [Common Infrastructure](#common-infrastructure)
6. [Implementation Roadmap](#implementation-roadmap)
7. [References and Further Reading](#references-and-further-reading)

---

## Theoretical Background and Literature

### Semi-Empirical Quantum Chemistry Fundamentals

**For theoretical chemists**: This section provides quick access to theoretical foundations.

#### Tight-Binding Approximations

**GFN1/GFN2 are based on extended tight-binding (xTB) theory**:

1. **Tight-Binding Hamiltonian**:
   - Approximation to DFT Kohn-Sham Hamiltonian
   - Minimal basis set (valence orbitals only)
   - Parameterized matrix elements instead of explicit integration
   - Suitable for large systems (104-106 atoms)

2. **Core Approximations**:
   - LCAO (Linear Combination of Atomic Orbitals)
   - Minimal basis: one Slater-type orbital per valence shell
   - Self-consistent charge (SCC) via electrostatic interactions
   - Empirical corrections for dispersion, H-bonds, halogen bonds

3. **Key References**:
   - **Tight-Binding Theory**: D. A. Papaconstantopoulos, *Handbook of the Band Structure of Elemental Solids* (2015)
   - **DFTB Foundations**: M. Elstner et al., *Phys. Rev. B* **1998**, 58, 7260 (SCC-DFTB)
   - **GFN Philosophy**: S. Grimme, *WIREs Comput. Mol. Sci.* **2011**, 1, 211 (Dispersion review)

#### NDDO Approximations (PMx Methods)

**PM3/AM1/MNDO are based on Neglect of Diatomic Differential Overlap**:

1. **NDDO Integral Approximation**:
   - Retain: (μA νA | λB σB) — two-center, same-shell two-electron integrals
   - Neglect: (μA νB | λC σD) — three- and four-center integrals
   - Reduces integral count from O(N⁴) to O(N²)

2. **Parameterization Philosophy**:
   - Fit to experimental thermochemistry (ΔH_f, ionization potentials)
   - Reproduce molecular geometries and frequencies
   - Element-specific Gaussian core-core corrections (PM3)

3. **Key References**:
   - **MNDO**: M. J. S. Dewar, W. Thiel, *J. Am. Chem. Soc.* **1977**, 99, 4899
   - **AM1**: M. J. S. Dewar et al., *J. Am. Chem. Soc.* **1985**, 107, 3902
   - **PM3**: J. J. P. Stewart, *J. Comput. Chem.* **1989**, 10, 209
   - **PM6**: J. J. P. Stewart, *J. Mol. Model.* **2007**, 13, 1173

#### Comparison: xTB vs. NDDO

| Aspect | GFN1/GFN2 (xTB) | PM3/AM1 (NDDO) |
|--------|-----------------|----------------|
| **Origin** | DFT tight-binding | Hartree-Fock NDDO |
| **Basis** | Minimal STO | Minimal GTO |
| **SCF** | Density-functional-like | Hartree-Fock RHF/UHF |
| **Target** | Geometries, non-covalent | Thermochemistry, ΔH_f |
| **Speed** | Very fast (~ms/atom) | Fast (~10 ms/atom) |
| **Accuracy** | Excellent for structures | Good for heats of formation |
| **Dispersion** | D3/D4 required | Implicit in parameterization |

**Further Reading**:
- J. P. Perdew, K. Schmidt, "Jacob's Ladder of DFT Approximations" (2001) — context for tight-binding
- F. Jensen, *Introduction to Computational Chemistry*, 3rd ed. (2017) — Chapter on semi-empirical methods
- C. J. Cramer, *Essentials of Computational Chemistry* (2004) — NDDO theory overview

---

## GFN2-xTB Implementation

### Overview

**GFN2-xTB** (Geometry, Frequency, Noncovalent 2 - Extended Tight Binding) is a modern semi-empirical quantum mechanical method optimized for:
- Accurate geometries and frequencies
- Noncovalent interactions (hydrogen bonds, π-stacking, dispersion)
- Large molecular systems (up to ~10,000 atoms)
- Computational efficiency (~1000x faster than DFT)

**References**:
- C. Bannwarth, S. Ehlert, S. Grimme, *J. Chem. Theory Comput.* **2019**, 15, 1652-1671
- TBLite repository: https://github.com/tblite/tblite

### Core Algorithm Structure

#### 1. Modular Component Architecture

```cpp
class GFN2 : public QMDriver {
public:
    // Initialization
    bool InitialiseMolecule() override;

    // Energy calculation
    double Calculation(bool gradient = false) override;

private:
    // Component setup
    void setupBasisSet();
    void setupHamiltonian();
    void setupCoulombInteractions();
    void setupRepulsion();
    void setupDispersion();
    void setupCoordinationNumbers();

    // SCF procedure
    bool runSCF();
    Matrix buildFockMatrix(const Matrix& density);
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& occupations);
    bool checkConvergence(const Matrix& old_density, const Matrix& new_density);

    // Energy components
    double calculateElectronicEnergy() const;
    double calculateRepulsionEnergy() const;
    double calculateDispersionEnergy() const;
    double calculateCoulombEnergy() const;

    // Gradient
    Matrix calculateGradient() const;

    // Data members
    GFN2Parameters m_params;
    BasisSet m_basis;
    Matrix m_hamiltonian;
    Matrix m_overlap;
    Matrix m_density;
    Vector m_coordination_numbers;

    double m_energy_electronic;
    double m_energy_repulsion;
    double m_energy_dispersion;
    double m_energy_coulomb;
};
```

#### 2. Basis Set Construction

**Slater-Type Orbitals (STO) → Gaussian-Type Orbitals (GTO) Conversion**

```cpp
struct BasisFunction {
    int atom_index;           // Which atom
    int angular_momentum;     // s=0, p=1, d=2, f=3
    int m_quantum_number;     // m_l = -l,...,+l
    int principal_qn;         // n = 1,2,3,...
    double slater_exponent;   // ζ (zeta)

    // Gaussian contraction (STO-6G approximation)
    std::vector<double> alpha;  // Gaussian exponents
    std::vector<double> coeff;  // Contraction coefficients
};

class BasisSet {
public:
    void constructFromMolecule(const Mol& mol, const GFN2Parameters& params);

    // Basis function access
    const BasisFunction& getFunction(int index) const;
    int getNumFunctions() const { return m_functions.size(); }

    // Shell information
    int getNumShells(int atomic_number) const;
    int getShellAngularMomentum(int shell) const;

private:
    std::vector<BasisFunction> m_functions;

    // Helper: STO-6G conversion
    void slaterToGaussian(double zeta, int n, int l,
                          std::vector<double>& alpha,
                          std::vector<double>& coeff);
};
```

**Shell Structure** (Z=1 to Z=86):
- **s-shell**: 1 function (m=0)
- **p-shell**: 3 functions (m=-1,0,+1)
- **d-shell**: 5 functions (m=-2,-1,0,+1,+2)

**Example for Carbon (Z=6)**:
```
Basis functions:
  1s: n=2, l=0, ζ=4.231
  2s: n=2, l=0, ζ=4.231
  2px: n=2, l=1, ζ=4.231
  2py: n=2, l=1, ζ=4.231
  2pz: n=2, l=1, ζ=4.231
Total: 5 basis functions
```

#### 3. Hamiltonian Matrix Construction

**Mathematical Form**:
```
H = H₀ + H_coulomb + H_repulsion + H_dispersion
```

**H₀ (Core Hamiltonian)**:

**Diagonal Elements** (Self-Energy):
```cpp
double getSelfEnergy(int shell, int element, double coordination_number) {
    double E_base = params.selfenergy[element][shell];  // Base orbital energy (eV)
    double CN_shift = params.kcn[element][shell] * coordination_number;

    // Convert eV → Hartree
    return (E_base + CN_shift) / 27.21138505;
}
```

**Off-Diagonal Elements** (Hopping Integrals):
```cpp
double getHamiltonianElement(int i, int j,
                              const BasisFunction& fi,
                              const BasisFunction& fj,
                              double distance) {
    if (i == j) {
        return getSelfEnergy(fi.shell, fi.atom_type, CN[fi.atom_index]);
    }

    // Overlap scaling factor
    double zi = fi.slater_exponent;
    double zj = fj.slater_exponent;
    double zij = std::pow(2.0 * std::sqrt(zi * zj) / (zi + zj), 0.5);

    // Pair coupling constant
    double kpair = params.kpair[fi.atom_type][fj.atom_type];

    // Shell coupling
    double kshell = getShellCoupling(fi.angular_momentum, fj.angular_momentum);

    // Electronegativity correction
    double ENi = params.pauling_EN[fi.atom_type];
    double ENj = params.pauling_EN[fj.atom_type];
    double enp = 1.0 + 0.02 * std::pow(ENi - ENj, 2);

    // Distance-dependent polynomial
    double poly_factor = getShellPolynomial(fi, fj, distance);

    // Overlap integral S_ij
    double S_ij = calculateOverlap(fi, fj);

    return zij * kpair * kshell * enp * poly_factor * S_ij;
}
```

**Shell Coupling Constants**:
```cpp
double getShellCoupling(int l1, int l2) {
    // kdiag = [1.85, 2.23, 2.23, 2.23, 2.23] for s,p,d,f,g

    // Special case: d-shell with s/p
    if ((l1 == 2 && l2 <= 1) || (l2 == 2 && l1 <= 1)) {
        return 2.0;
    }

    // Average of diagonal elements
    return 0.5 * (kdiag[l1] + kdiag[l2]);
}
```

**Distance-Dependent Polynomial**:
```cpp
double getShellPolynomial(const BasisFunction& fi,
                          const BasisFunction& fj,
                          double r) {
    // Polynomial coefficients from params.shpoly[element][shell][0..2]
    double a0 = params.shpoly[fi.atom_type][fi.shell][0];
    double a1 = params.shpoly[fi.atom_type][fi.shell][1];
    double a2 = params.shpoly[fi.atom_type][fi.shell][2];

    return a0 + a1 * r + a2 * r * r;
}
```

#### 4. Coulomb Interactions (Three Components)

**Component 1: Effective Coulomb (ES2)**
```cpp
double calculateEffectiveCoulomb() {
    double E_coulomb = 0.0;

    for (int A = 0; A < natoms; ++A) {
        for (int B = A+1; B < natoms; ++B) {
            double rAB = distance(A, B);

            // Average shell hardness
            double gammaAB = 0.5 * (params.hardness[atom[A]] +
                                    params.hardness[atom[B]]);

            // Global exponential decay
            double gexp = 2.0;
            double kernel = 1.0 / std::pow(rAB*rAB*rAB + 1.0/std::pow(gammaAB, gexp), 1.0/gexp);

            // Charge-charge interaction
            E_coulomb += charges[A] * charges[B] * kernel;
        }
    }

    return E_coulomb;
}
```

**Component 2: Third-Order Onsite (ES3)**
```cpp
double calculateThirdOrderOnsite() {
    double E3 = 0.0;

    for (int A = 0; A < natoms; ++A) {
        double qA = charges[A];
        double dU = params.hubbard_derivative[atom[A]];  // Hubbard ∂U/∂n

        // Third-order contribution: E₃ ∝ dU * q³
        E3 += (1.0/6.0) * dU * qA * qA * qA;
    }

    return E3;
}
```

**Component 3: Damped Multipole (AES2)**
```cpp
double calculateDampedMultipole() {
    double E_multi = 0.0;

    for (int A = 0; A < natoms; ++A) {
        for (int B = A+1; B < natoms; ++B) {
            double rAB = distance(A, B);

            // Dipole and quadrupole kernels
            double kdipole = params.dkernel[atom[A]][atom[B]];
            double kquad = params.qkernel[atom[A]][atom[B]];

            // Damping functions (Tang-Toennies)
            double dmp3 = dampingFunction(rAB, 3.0, kdipole);
            double dmp5 = dampingFunction(rAB, 5.0, kquad);

            // Multipole contributions
            double dipole_term = charges[A] * dipoles[B] / std::pow(rAB, 3) * dmp3;
            double quad_term = charges[A] * quadrupoles[B] / std::pow(rAB, 5) * dmp5;

            E_multi += dipole_term + quad_term;
        }
    }

    return E_multi;
}
```

#### 5. Repulsion Energy

**Pairwise Atomic Repulsion**:
```cpp
double calculateRepulsion() {
    double E_rep = 0.0;

    for (int A = 0; A < natoms; ++A) {
        for (int B = A+1; B < natoms; ++B) {
            double rAB = distance(A, B);

            // Element-specific parameters
            double alpha_A = params.rep_alpha[atom[A]];
            double alpha_B = params.rep_alpha[atom[B]];
            double zeff_A = params.rep_zeff[atom[A]];
            double zeff_B = params.rep_zeff[atom[B]];

            // Exponential repulsion
            double rep_A = zeff_A * std::exp(-alpha_A * rAB);
            double rep_B = zeff_B * std::exp(-alpha_B * rAB);

            E_rep += rep_A + rep_B;
        }
    }

    return E_rep;
}
```

#### 6. Dispersion Correction (D4)

**Integration with DFT-D4**:
```cpp
double calculateDispersion() {
    // GFN2 uses DFT-D4 dispersion correction
    // Parameters:
    const double s6 = 1.0;   // C₆ scaling
    const double s8 = 2.7;   // C₈ scaling
    const double a1 = 0.52;  // BJ damping
    const double a2 = 5.0;   // BJ damping
    const double s9 = 5.0;   // ATM C₉ scaling

    // Use existing DFT-D4 implementation
    // (Curcuma already has dftd4interface.h/cpp)

    D4Interface d4(s6, s8, a1, a2, s9);
    d4.InitialiseMolecule(molecule);
    return d4.Calculation();
}
```

#### 7. Coordination Number Calculation

**Exponential Coordination Number**:
```cpp
std::vector<double> calculateCoordinationNumbers() {
    std::vector<double> CN(natoms, 0.0);

    const double k1 = 16.0;  // Steepness parameter
    const double k2 = 4.0/3.0;  // Range parameter

    for (int A = 0; A < natoms; ++A) {
        for (int B = 0; B < natoms; ++B) {
            if (A == B) continue;

            double rAB = distance(A, B);
            double rcovA = params.covalent_radius[atom[A]];
            double rcovB = params.covalent_radius[atom[B]];
            double rcovAB = rcovA + rcovB;

            // Counting function with exponential decay
            double count = 1.0 / (1.0 + std::exp(-k1 * (rcovAB / rAB - 1.0)));

            CN[A] += std::pow(count, k2);
        }
    }

    return CN;
}
```

#### 8. SCF Convergence Procedure

**Self-Consistent Field Loop**:
```cpp
bool runSCF() {
    const int max_iterations = 100;
    const double convergence_threshold = 1.0e-6;

    // Initial guess: superposition of atomic densities (SAD)
    Matrix density = buildInitialDensity();

    for (int iter = 0; iter < max_iterations; ++iter) {
        // 1. Build Fock matrix
        Matrix fock = buildFockMatrix(density);

        // 2. Solve generalized eigenvalue problem: FC = SCE
        Vector eigenvalues;
        Matrix eigenvectors;
        solveGeneralizedEigenvalue(fock, m_overlap, eigenvalues, eigenvectors);

        // 3. Build density matrix from occupied orbitals
        Matrix new_density = buildDensityMatrix(eigenvectors, eigenvalues);

        // 4. Apply damping to improve convergence
        double damping = 0.4;  // Configurable
        new_density = damping * new_density + (1.0 - damping) * density;

        // 5. Check convergence
        double density_change = (new_density - density).norm();

        if (density_change < convergence_threshold) {
            density = new_density;
            m_density = density;
            return true;  // Converged!
        }

        density = new_density;

        // Verbosity output
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("SCF iter {}: ΔP = {:.6e}",
                                           iter, density_change));
        }
    }

    CurcumaLogger::error("SCF did not converge within maximum iterations");
    return false;  // Failed to converge
}
```

**Fock Matrix Construction**:
```cpp
Matrix buildFockMatrix(const Matrix& density) {
    Matrix fock = m_hamiltonian;  // Start with core Hamiltonian

    // Add Coulomb contributions
    for (int i = 0; i < nbasis; ++i) {
        for (int j = 0; j < nbasis; ++j) {
            // Two-electron integrals weighted by density
            double coulomb_contrib = 0.0;

            for (int k = 0; k < nbasis; ++k) {
                for (int l = 0; l < nbasis; ++l) {
                    // Shell-resolved Coulomb kernel
                    double gamma_ijkl = getCoulombKernel(i, j, k, l);
                    coulomb_contrib += density(k, l) * gamma_ijkl;
                }
            }

            fock(i, j) += coulomb_contrib;
        }
    }

    return fock;
}
```

#### 9. Gradient Calculation

**Analytical Gradients**:
```cpp
Matrix calculateGradient() const {
    Matrix gradient = Matrix::Zero(natoms, 3);

    // 1. Electronic gradient (Hellmann-Feynman theorem)
    for (int A = 0; A < natoms; ++A) {
        for (int xyz = 0; xyz < 3; ++xyz) {
            // dH/dR_A contribution
            Matrix dH_dR = calculateHamiltonianDerivative(A, xyz);
            double grad_electronic = (m_density.cwiseProduct(dH_dR)).sum();

            gradient(A, xyz) += grad_electronic;
        }
    }

    // 2. Repulsion gradient
    for (int A = 0; A < natoms; ++A) {
        for (int B = 0; B < natoms; ++B) {
            if (A == B) continue;

            Vector3d rAB = positions[B] - positions[A];
            double dist = rAB.norm();

            // Derivative of exponential repulsion
            double alpha_A = params.rep_alpha[atom[A]];
            double zeff_A = params.rep_zeff[atom[A]];
            double dV_dr = -alpha_A * zeff_A * std::exp(-alpha_A * dist);

            Vector3d force = dV_dr * rAB / dist;
            gradient.row(A) += force.transpose();
        }
    }

    // 3. Dispersion gradient
    Matrix disp_gradient = calculateDispersionGradient();
    gradient += disp_gradient;

    return gradient;
}
```

### GFN2 Parameters

**Complete Parameter Set** (already in `gfn2-xtb_param.hpp`):

```cpp
struct GFN2Parameters {
    // Per-element parameters (Z=1..86)
    std::array<double, 87> electronegativity;     // Pauling EN
    std::array<double, 87> hardness;              // Chemical hardness (eV)
    std::array<double, 87> hubbard_derivative;    // ∂U/∂n
    std::array<double, 87> rep_alpha;             // Repulsion exponent
    std::array<double, 87> rep_zeff;              // Effective nuclear charge
    std::array<double, 87> covalent_radius;       // For CN calculation

    // Shell-resolved parameters (element × shell)
    std::array<std::array<double, 3>, 87> selfenergy;     // Orbital energies
    std::array<std::array<double, 3>, 87> kcn;            // CN shift factors
    std::array<std::array<double, 4>, 87> slater_exp;     // ζ exponents
    std::array<std::array<int, 4>, 87> principal_qn;      // n quantum numbers

    // Polynomial coefficients (element × shell × 3)
    std::array<std::array<std::array<double, 3>, 3>, 87> shpoly;

    // Global constants
    static constexpr double kdiag[5] = {1.85, 2.23, 2.23, 2.23, 2.23};
    static constexpr double enscale = 0.02;
    static constexpr double wexp = 0.5;
    static constexpr double gexp = 2.0;
};
```

---

## GFN1-xTB Implementation

### Overview

**GFN1-xTB** is the predecessor of GFN2, simpler but still highly effective:
- Faster than GFN2 (~2x speed)
- Slightly less accurate for noncovalent interactions
- No third-order charge terms
- Uses D3 instead of D4 dispersion

**Key Differences from GFN2**:

| Feature | GFN1 | GFN2 |
|---------|------|------|
| **Dispersion** | D3-ATM (s9=0.0) | D4 (s9=5.0) |
| **Coulomb** | ES2 + ES3 only | ES2 + ES3 + AES2 (multipoles) |
| **Hamiltonian** | Simpler polynomial scaling | Distance-dependent polynomials |
| **Halogen bonds** | Explicit correction | Implicit in parameterization |
| **Speed** | ~1.0 (reference) | ~0.5 (slower) |
| **Accuracy** | Good | Excellent |

### GFN1 Algorithm Modifications

**1. Hamiltonian Scaling** (Simpler):
```cpp
double getHamiltonianElement_GFN1(int i, int j, const BasisFunction& fi, const BasisFunction& fj) {
    // ... (same overlap scaling as GFN2)

    // Valence vs. core distinction
    bool fi_valence = isValenceShell(fi);
    bool fj_valence = isValenceShell(fj);

    double kshell;
    if (fi_valence && fj_valence) {
        kshell = getShellCoupling(fi.angular_momentum, fj.angular_momentum);
    } else if (!fi_valence && !fj_valence) {
        kshell = 2.85;  // Core-core constant
    } else {
        // Valence-core: average
        double k_val = getShellCoupling(fi.angular_momentum, fj.angular_momentum);
        kshell = 0.5 * (k_val + 2.85);
    }

    // Electronegativity (negative sign!)
    double ENi = params.pauling_EN[fi.atom_type];
    double ENj = params.pauling_EN[fj.atom_type];
    double enp = 1.0 - 0.007 * std::pow(ENi - ENj, 2);  // Note: minus!

    return zij * kpair * kshell * enp * S_ij;
}
```

**2. Dispersion (D3 instead of D4)**:
```cpp
double calculateDispersion_GFN1() {
    // GFN1 uses DFT-D3 with ATM
    const double s6 = 1.0;
    const double s8 = 2.4;   // Different from GFN2!
    const double a1 = 0.63;  // Different from GFN2!
    const double a2 = 5.0;
    const double s9 = 0.0;   // No ATM for GFN1

    D3Interface d3(s6, s8, a1, a2, s9);
    d3.InitialiseMolecule(molecule);
    return d3.Calculation();
}
```

**3. Halogen Bond Correction**:
```cpp
double calculateHalogenCorrection() {
    double E_hal = 0.0;

    const double damping = 0.44;
    const double radscale = 1.3;

    for (int A = 0; A < natoms; ++A) {
        if (!isHalogen(atom[A])) continue;

        for (int B = 0; B < natoms; ++B) {
            if (A == B) continue;

            double rAB = distance(A, B);
            double rcovAB = radscale * (params.covalent_radius[atom[A]] +
                                        params.covalent_radius[atom[B]]);

            // Tang-Toennies damping
            double damp = 1.0 - std::exp(-damping * (rAB / rcovAB - 1.0));

            // Element-specific halogen bond strength
            double strength = params.halogen_strength[atom[A]];

            E_hal += strength * std::pow(rcovAB / rAB, 6) * damp;
        }
    }

    return E_hal;
}
```

### GFN1 Parameters

**Key Parameter Differences**:
```cpp
struct GFN1Parameters {
    // Same structure as GFN2, but different values:

    // Diagonal coupling constants
    static constexpr double kdiag[5] = {1.85, 2.25, 2.0, 2.0, 2.0};  // Different!

    // Electronegativity scaling (negative!)
    static constexpr double enscale = -0.007;  // Different sign!

    // Core-core constant
    static constexpr double kdiff = 2.85;

    // Halogen bond parameters (not in GFN2)
    std::array<double, 87> halogen_strength;
    static constexpr double halogen_damping = 0.44;
    static constexpr double halogen_radscale = 1.3;
};
```

---

## PMx Methods Implementation

### Overview

**PMx Methods** (PM3, PM6, AM1, MNDO, etc.) are NDDO-based semi-empirical methods:
- **NDDO**: Neglect of Diatomic Differential Overlap
- Use Gaussian-type orbitals (GTOs)
- Parameterized for thermochemistry (heats of formation)
- Support analytical gradients
- Closed-shell (RHF) and open-shell (UHF)

**Available Methods**:
- **MNDO** (1977): Modified Neglect of Diatomic Overlap
- **AM1** (1985): Austin Model 1
- **PM3** (1989): Parameterized Model 3
- **PM6** (2007): Parameterized Model 6
- **PM7** (2013): Parameterized Model 7

### Core Algorithm Structure

#### NDDO Approximation

**Two-Center Two-Electron Integrals**:

```
(μA νA | λB σB) ≠ 0  (two-center integrals retained)
(μA νB | λC σD) = 0  (three- and four-center integrals neglected)
```

**Integral Types**:
```cpp
enum IntegralType {
    SS_SS,  // (ss|ss): s-s coulomb
    SS_SP,  // (ss|spσ): s-pσ coulomb
    SP_SP,  // (spσ|spσ): pσ-pσ coulomb
    SS_PP,  // (ss|pπpπ): s-pπpπ coulomb
    SP_PP,  // (spσ|pπpπ): pσ-pπpπ coulomb
    PP_PP,  // (pπpπ|pπpπ): pπ-pπ coulomb
    SS,     // (ss|ss) one-center
    PP      // (pp|pp) one-center
};
```

#### Core Hamiltonian

**One-Electron Terms**:
```cpp
double getCoreHamiltonian_PM3(int mu, int nu) {
    int A = atom_of_orbital(mu);
    int B = atom_of_orbital(nu);

    if (A == B) {
        // Diagonal: ionization potential
        if (mu == nu) {
            return -params.ionization_potential[atom[A]][orbital_type(mu)];
        }
        // Off-diagonal on same atom: zero (orthogonal AOs)
        return 0.0;
    }

    // Off-diagonal between atoms
    // Resonance integral
    double beta_A = params.beta[atom[A]][orbital_type(mu)];
    double beta_B = params.beta[atom[B]][orbital_type(nu)];
    double beta_AB = 0.5 * (beta_A + beta_B);

    // Overlap integral
    double S_mu_nu = calculateOverlap_GTO(mu, nu);

    return beta_AB * S_mu_nu;
}
```

**Core-Core Repulsion**:
```cpp
double getCoreRepulsion_PM3(int A, int B) {
    double rAB = distance(A, B);

    // Klopman-Ohno formula with corrections
    double gamma_AB = calculateGamma(A, B, rAB);

    // Nuclear charges
    int ZA = core_charge[atom[A]];  // Valence electrons
    int ZB = core_charge[atom[B]];

    // Base repulsion
    double V_core = ZA * ZB * gamma_AB;

    // Gaussian correction terms (PM3-specific)
    V_core += getGaussianCorrection(A, B, rAB);

    return V_core;
}
```

**Gaussian Correction (PM3)**:
```cpp
double getGaussianCorrection(int A, int B, double r) {
    double correction = 0.0;

    // PM3 uses up to 4 Gaussian terms per atom pair
    for (int k = 0; k < 4; ++k) {
        double alpha_Ak = params.alpha[atom[A]][k];
        double alpha_Bk = params.alpha[atom[B]][k];
        double x_Ak = params.x[atom[A]][k];
        double x_Bk = params.x[atom[B]][k];

        if (alpha_Ak > 0.0) {
            correction += x_Ak * std::exp(-alpha_Ak * r * r);
        }
        if (alpha_Bk > 0.0) {
            correction += x_Bk * std::exp(-alpha_Bk * r * r);
        }
    }

    return correction;
}
```

#### Fock Matrix

**Two-Electron Contributions**:
```cpp
double getFockElement_PM3(int mu, int nu, const Matrix& density) {
    double F_mu_nu = getCoreHamiltonian_PM3(mu, nu);

    int A = atom_of_orbital(mu);
    int B = atom_of_orbital(nu);

    // Sum over density matrix
    for (int lambda = 0; lambda < nbasis; ++lambda) {
        for (int sigma = 0; sigma < nbasis; ++sigma) {
            int C = atom_of_orbital(lambda);
            int D = atom_of_orbital(sigma);

            // Coulomb integral: (μν|λσ)
            double J_integral = getTwoElectronIntegral(mu, nu, lambda, sigma);

            // Exchange integral: (μσ|λν)
            double K_integral = getTwoElectronIntegral(mu, sigma, lambda, nu);

            // Fock matrix element
            F_mu_nu += density(lambda, sigma) * (J_integral - 0.5 * K_integral);
        }
    }

    return F_mu_nu;
}
```

**Two-Electron Integral Evaluation**:
```cpp
double getTwoElectronIntegral(int mu, int nu, int lambda, int sigma) {
    int A = atom_of_orbital(mu);
    int B = atom_of_orbital(nu);
    int C = atom_of_orbital(lambda);
    int D = atom_of_orbital(sigma);

    // NDDO: three- and four-center integrals are zero
    std::set<int> centers = {A, B, C, D};
    if (centers.size() > 2) {
        return 0.0;
    }

    // One-center integral
    if (centers.size() == 1) {
        return getOneCenterIntegral(mu, nu, lambda, sigma);
    }

    // Two-center integral
    return getTwoCenterIntegral(mu, nu, lambda, sigma);
}
```

#### SCF Procedure (RHF)

**Closed-Shell Restricted Hartree-Fock**:
```cpp
bool runSCF_RHF() {
    const int max_iter = 100;
    const double threshold = 1.0e-6;

    // Initial guess: extended Hückel or core Hamiltonian
    Matrix density = buildInitialDensity_RHF();

    for (int iter = 0; iter < max_iter; ++iter) {
        // 1. Build Fock matrix
        Matrix fock(nbasis, nbasis);
        for (int mu = 0; mu < nbasis; ++mu) {
            for (int nu = 0; nu < nbasis; ++nu) {
                fock(mu, nu) = getFockElement_PM3(mu, nu, density);
            }
        }

        // 2. Solve FC = SCE
        Vector eigenvalues;
        Matrix coefficients;
        solveGeneralizedEigenvalue(fock, m_overlap, eigenvalues, coefficients);

        // 3. Build new density (closed-shell: 2 electrons per orbital)
        Matrix new_density = Matrix::Zero(nbasis, nbasis);
        int n_occ = num_electrons / 2;  // Doubly occupied orbitals

        for (int mu = 0; mu < nbasis; ++mu) {
            for (int nu = 0; nu < nbasis; ++nu) {
                for (int i = 0; i < n_occ; ++i) {
                    new_density(mu, nu) += 2.0 * coefficients(mu, i) * coefficients(nu, i);
                }
            }
        }

        // 4. Check convergence
        double delta = (new_density - density).norm();
        if (delta < threshold) {
            m_density = new_density;
            m_mo_energies = eigenvalues;
            m_mo_coefficients = coefficients;
            return true;
        }

        // 5. Damping
        density = 0.5 * new_density + 0.5 * density;
    }

    return false;  // Failed to converge
}
```

### PM3 Parameters

**Element-Specific Parameters** (Example for Carbon):
```cpp
struct PM3AtomParameters {
    // Orbital parameters
    double U_ss;   // s-orbital ionization potential (eV)
    double U_pp;   // p-orbital ionization potential (eV)
    double beta_s; // s-orbital resonance integral (eV)
    double beta_p; // p-orbital resonance integral (eV)

    // Gaussian exponents (Slater orbitals)
    double zeta_s; // s-orbital exponent
    double zeta_p; // p-orbital exponent

    // Coulomb parameters
    double g_ss;   // One-center s-s coulomb
    double g_pp;   // One-center p-p coulomb
    double g_sp;   // One-center s-p coulomb
    double g_pp2;  // One-center p-p exchange
    double h_sp;   // One-center s-p resonance

    // Core-core Gaussian corrections (4 terms)
    std::array<double, 4> alpha;  // Exponents
    std::array<double, 4> x;      // Coefficients

    // Core charge
    int core_charge;  // Number of valence electrons
};

// Example: Carbon PM3 parameters
PM3AtomParameters carbon_pm3 = {
    .U_ss = -47.270320,  // eV
    .U_pp = -36.266918,  // eV
    .beta_s = -11.910015,
    .beta_p = -9.802755,
    .zeta_s = 1.565085,
    .zeta_p = 1.842345,
    .g_ss = 11.200708,   // eV
    .g_pp = 10.265027,
    .g_sp = 10.778326,
    .g_pp2 = 9.627141,
    .h_sp = 2.489178,
    .alpha = {0.011355, 0.004196, 0.0, 0.0},
    .x = {1.561231, -1.306624, 0.0, 0.0},
    .core_charge = 4
};
```

---

## Common Infrastructure

### Overlap Matrix Calculation

**Gaussian-Type Orbital Overlap**:
```cpp
double calculateOverlap_GTO(const BasisFunction& fi, const BasisFunction& fj) {
    double S = 0.0;

    // Contract over Gaussian primitives
    for (size_t i = 0; i < fi.alpha.size(); ++i) {
        for (size_t j = 0; j < fj.alpha.size(); ++j) {
            double alpha_i = fi.alpha[i];
            double alpha_j = fj.alpha[j];
            double coeff_i = fi.coeff[i];
            double coeff_j = fj.coeff[j];

            // Primitive overlap integral
            double S_prim = primitiveOverlap(alpha_i, fi.center, fi.angular_momentum,
                                            alpha_j, fj.center, fj.angular_momentum);

            S += coeff_i * coeff_j * S_prim;
        }
    }

    return S;
}
```

### Generalized Eigenvalue Solver

**Löwdin Orthogonalization**:
```cpp
void solveGeneralizedEigenvalue(const Matrix& H, const Matrix& S,
                                Vector& eigenvalues, Matrix& eigenvectors) {
    // 1. Diagonalize overlap matrix
    Eigen::SelfAdjointEigenSolver<Matrix> solver_S(S);
    Vector s_eigenvalues = solver_S.eigenvalues();
    Matrix s_eigenvectors = solver_S.eigenvectors();

    // 2. Form S^(-1/2)
    Matrix S_inv_sqrt = Matrix::Zero(S.rows(), S.cols());
    for (int i = 0; i < s_eigenvalues.size(); ++i) {
        if (s_eigenvalues(i) > 1.0e-10) {
            S_inv_sqrt += (1.0 / std::sqrt(s_eigenvalues(i))) *
                          s_eigenvectors.col(i) * s_eigenvectors.col(i).transpose();
        }
    }

    // 3. Transform Hamiltonian: H' = S^(-1/2) H S^(-1/2)
    Matrix H_prime = S_inv_sqrt.transpose() * H * S_inv_sqrt;

    // 4. Solve standard eigenvalue problem: H' C' = C' E
    Eigen::SelfAdjointEigenSolver<Matrix> solver_H(H_prime);
    eigenvalues = solver_H.eigenvalues();
    Matrix C_prime = solver_H.eigenvectors();

    // 5. Back-transform: C = S^(-1/2) C'
    eigenvectors = S_inv_sqrt * C_prime;
}
```

### Logging Integration

**Universal Verbosity Support**:
```cpp
void printSCFIteration(int iter, double energy, double density_change) {
    // Level 3: Detailed iteration info
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("SCF Iteration {}", iter));
        CurcumaLogger::param("energy", fmt::format("{:.8f} Eh", energy));
        CurcumaLogger::param("delta_P", fmt::format("{:.6e}", density_change));
    }
}

void printFinalResults() {
    // Level 1: Final energies
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(m_total_energy, "Total Energy");
    }

    // Level 2: Energy decomposition
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", m_energy_electronic));
        CurcumaLogger::param("repulsion", fmt::format("{:.6f} Eh", m_energy_repulsion));
        CurcumaLogger::param("dispersion", fmt::format("{:.6f} Eh", m_energy_dispersion));
        CurcumaLogger::param("HOMO", fmt::format("{:.4f} eV", getHOMOEnergy() * 27.211));
        CurcumaLogger::param("LUMO", fmt::format("{:.4f} eV", getLUMOEnergy() * 27.211));
        CurcumaLogger::param("gap", fmt::format("{:.4f} eV", getHOMOLUMOGap() * 27.211));
    }
}
```

---

## Implementation Roadmap

### Phase 1: GFN2 Core Implementation (Priority)

**Tasks**:
1. ✅ Create `gfn2.h` / `gfn2.cpp` (inherits from `QMDriver`)
2. ✅ Create `gfn2_parameters.h` (use existing `gfn2-xtb_param.hpp` as base)
3. ✅ Implement basis set construction
4. ✅ Implement overlap matrix calculation
5. ✅ Implement Hamiltonian matrix construction
6. ✅ Implement coordination number calculation
7. ✅ Implement SCF convergence loop
8. ✅ Integrate D4 dispersion
9. ✅ Implement repulsion energy
10. ✅ Implement Coulomb interactions (ES2, ES3, AES2)
11. ✅ Implement gradient calculation
12. ✅ Create `gfn2_method.h` / `gfn2_method.cpp` wrapper
13. ✅ Integrate into `MethodFactory`
14. ✅ Write unit tests
15. ✅ Validate against TBLite reference calculations

**File Structure**:
```
src/core/energy_calculators/qm_methods/
├── gfn2.h
├── gfn2.cpp
├── gfn2_parameters.h
├── gfn2_method.h
├── gfn2_method.cpp
└── tests/
    └── test_gfn2.cpp
```

### Phase 2: GFN1 Implementation

**Tasks**:
1. Adapt GFN2 code for GFN1 differences
2. Implement halogen bond correction
3. Switch D4 → D3 dispersion
4. Adjust Hamiltonian scaling (valence/core distinction)
5. Update parameters
6. Validate against TBLite

### Phase 3: PM3 Implementation

**Tasks**:
1. Create `pm3.h` / `pm3.cpp`
2. Implement NDDO integral evaluation
3. Implement Gaussian corrections
4. Implement RHF SCF procedure
5. Add PM3 parameters for common elements
6. Validate against MOPAC

### Phase 4: Additional PMx Methods

**Tasks**:
1. Add AM1 (similar to PM3, different parameters)
2. Add MNDO (no Gaussian corrections)
3. Add PM6 (extended parameters)
4. Unified PMx base class with method-specific parameterizations

### Phase 5: MethodFactory Integration

**Update `method_factory.cpp`**:
```cpp
// Priority-based GFN2 resolution
{
    "gfn2",
    {
        { "Native", [](const json& config) { return std::make_unique<GFN2Method>(config); } },
        { "TBLite", [](const json& config) { return hasTBLite() ? std::make_unique<TBLiteMethod>("gfn2", config) : nullptr; } },
        { "Ulysses", [](const json& config) { return hasUlysses() ? std::make_unique<UlyssesMethod>("ugfn2", config) : nullptr; } },
        { "XTB", [](const json& config) { return hasXTB() ? std::make_unique<XTBMethod>("gfn2", config) : nullptr; } }
    }
},

// Native PM3
{
    "pm3",
    {
        { "Native", [](const json& config) { return std::make_unique<PM3Method>(config); } },
        { "Ulysses", [](const json& config) { return hasUlysses() ? std::make_unique<UlyssesMethod>("pm3", config) : nullptr; } }
    }
}
```

### Testing Strategy

**Unit Tests**:
- Individual component tests (basis set, overlap, Hamiltonian)
- SCF convergence tests
- Gradient accuracy tests (numerical vs. analytical)

**Integration Tests**:
- Small molecule benchmarks (H₂O, CH₄, benzene)
- Comparison with TBLite/MOPAC reference energies
- Geometry optimization convergence

**Validation Data**:
```
test_cases/native_qm/
├── water_gfn2.xyz
├── methane_gfn2.xyz
├── benzene_gfn2.xyz
├── reference_energies.json
└── reference_gradients.json
```

---

## References and Further Reading

### Primary Method Publications

#### GFN2-xTB (Extended Tight-Binding GFN2)
**Primary Reference**:
- **Bannwarth, C.; Ehlert, S.; Grimme, S.** "GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions" *J. Chem. Theory Comput.* **2019**, 15 (3), 1652-1671. DOI: [10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176)

**Implementation Source**:
- **TBLite Library**: https://github.com/tblite/tblite (LGPL-3.0)
  - Fortran reference implementation by Sebastian Ehlert et al.
  - Parameters: `src/tblite/param/gfn2.toml`
  - Hamiltonian: `src/tblite/xtb/gfn2.f90`

**Supporting Papers**:
- **Pracht, P.; Caldeweyher, E.; Ehlert, S.; Grimme, S.** "A Robust Non-Self-Consistent Tight-Binding Quantum Chemistry Method for large Molecules" *ChemRxiv* **2019**. (GFN-FF, needed for comparison)

#### GFN1-xTB (Extended Tight-Binding GFN1)
**Primary Reference**:
- **Grimme, S.; Bannwarth, C.; Shushkov, P.** "A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1-86)" *J. Chem. Theory Comput.* **2017**, 13 (5), 1989-2009. DOI: [10.1021/acs.jctc.7b00118](https://doi.org/10.1021/acs.jctc.7b00118)

**Implementation Source**:
- **TBLite Library**: `src/tblite/xtb/gfn1.f90`
- **XTB Program**: https://github.com/grimme-lab/xtb (LGPL-3.0)

#### PMx Semi-Empirical Methods

**PM3** (Parameterized Model 3):
- **Stewart, J. J. P.** "Optimization of parameters for semiempirical methods I. Method" *J. Comput. Chem.* **1989**, 10 (2), 209-220. DOI: [10.1002/jcc.540100208](https://doi.org/10.1002/jcc.540100208)
- **Stewart, J. J. P.** "Optimization of parameters for semiempirical methods II. Applications" *J. Comput. Chem.* **1989**, 10 (2), 221-264.

**AM1** (Austin Model 1):
- **Dewar, M. J. S.; Zoebisch, E. G.; Healy, E. F.; Stewart, J. J. P.** "Development and use of quantum mechanical molecular models. 76. AM1: a new general purpose quantum mechanical molecular model" *J. Am. Chem. Soc.* **1985**, 107 (13), 3902-3909. DOI: [10.1021/ja00299a024](https://doi.org/10.1021/ja00299a024)

**PM6** (Parameterized Model 6):
- **Stewart, J. J. P.** "Optimization of parameters for semiempirical methods V: Modification of NDDO approximations and application to 70 elements" *J. Mol. Model.* **2007**, 13 (12), 1173-1213. DOI: [10.1007/s00894-007-0233-4](https://doi.org/10.1007/s00894-007-0233-4)

**MNDO** (Modified Neglect of Diatomic Overlap):
- **Dewar, M. J. S.; Thiel, W.** "Ground states of molecules. 38. The MNDO method. Approximations and parameters" *J. Am. Chem. Soc.* **1977**, 99 (15), 4899-4907. DOI: [10.1021/ja00457a004](https://doi.org/10.1021/ja00457a004)

**Implementation Source**:
- **MOPAC**: https://github.com/openmopac/mopac (LGPL-3.0)
  - Fortran implementation by James J. P. Stewart et al.
  - Parameters: `src/Parameters/` directory

### Dispersion Corrections

**DFT-D3** (Grimme's D3 Dispersion):
- **Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H.** "A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu" *J. Chem. Phys.* **2010**, 132 (15), 154104. DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)
- **Grimme, S.; Ehrlich, S.; Goerigk, L.** "Effect of the damping function in dispersion corrected density functional theory" *J. Comput. Chem.* **2011**, 32 (7), 1456-1465. (BJ damping)

**DFT-D4** (Next-Generation D4):
- **Caldeweyher, E.; Ehlert, S.; Hansen, A.; Neugebauer, H.; Spicher, S.; Bannwarth, C.; Grimme, S.** "A generally applicable atomic-charge dependent London dispersion correction" *J. Chem. Phys.* **2019**, 150 (15), 154122. DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)

**Implementation Sources**:
- **simple-dftd3**: https://github.com/dftd3/simple-dftd3 (LGPL-3.0)
- **cpp-d4**: https://github.com/dftd4/cpp-d4 (LGPL-3.0)

### Theoretical Foundations

#### Tight-Binding Theory
- **Papaconstantopoulos, D. A.** *Handbook of the Band Structure of Elemental Solids: From Z=1 to Z=112*, 2nd ed.; Springer, 2015. ISBN: 978-1-4939-1647-3
- **Harrison, W. A.** *Electronic Structure and the Properties of Solids: The Physics of the Chemical Bond*; Dover Publications, 1989. (Classical tight-binding reference)

#### Density-Functional Tight-Binding (DFTB)
- **Elstner, M.; Porezag, D.; Jungnickel, G.; Elsner, J.; Haugk, M.; Frauenheim, T.; Suhai, S.; Seifert, G.** "Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties" *Phys. Rev. B* **1998**, 58 (11), 7260-7268. DOI: [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)
- **Gaus, M.; Cui, Q.; Elstner, M.** "DFTB3: Extension of the Self-Consistent-Charge Density-Functional Tight-Binding Method (SCC-DFTB)" *J. Chem. Theory Comput.* **2011**, 7 (4), 931-948. (Third-order corrections)

#### NDDO and Semi-Empirical Methods
- **Pople, J. A.; Santry, D. P.; Segal, G. A.** "Approximate Self-Consistent Molecular Orbital Theory. I. Invariant Procedures" *J. Chem. Phys.* **1965**, 43 (10), S129-S135. (NDDO foundations)
- **Thiel, W.** "Semiempirical quantum–chemical methods" *WIREs Comput. Mol. Sci.* **2014**, 4 (2), 145-157. DOI: [10.1002/wcms.1161](https://doi.org/10.1002/wcms.1161) (Excellent modern review)

### Textbooks for Theoretical Chemists

#### Computational Chemistry Fundamentals
- **Jensen, F.** *Introduction to Computational Chemistry*, 3rd ed.; Wiley, 2017. ISBN: 978-1-118-82599-0
  - Chapter 5: Semi-empirical methods (PM3, AM1, MNDO)
  - Chapter 6: Hartree-Fock theory (for understanding SCF)
  - Chapter 7: DFT foundations (context for tight-binding)

- **Cramer, C. J.** *Essentials of Computational Chemistry: Theories and Models*, 2nd ed.; Wiley, 2004. ISBN: 978-0-470-09182-1
  - Chapter 8: Semi-empirical methods (detailed NDDO derivation)
  - Chapter 9: Basis sets (STO vs. GTO)

- **Szabo, A.; Ostlund, N. S.** *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*; Dover, 1996. ISBN: 978-0-486-69186-2
  - Chapter 3: Hartree-Fock theory (SCF procedure fundamentals)
  - Appendix A: Second quantization (for understanding many-body theory)

#### Advanced Topics
- **Martin, R. M.** *Electronic Structure: Basic Theory and Practical Methods*, 2nd ed.; Cambridge University Press, 2020.
  - Chapter 18: Tight-binding methods (comprehensive mathematical treatment)

- **Koch, W.; Holthausen, M. C.** *A Chemist's Guide to Density Functional Theory*, 2nd ed.; Wiley-VCH, 2001.
  - Context for understanding tight-binding as DFT approximation

### Numerical Methods and Linear Algebra

#### Eigenvalue Problems
- **Press, W. H.; Teukolsky, S. A.; Vetterling, W. T.; Flannery, B. P.** *Numerical Recipes: The Art of Scientific Computing*, 3rd ed.; Cambridge University Press, 2007.
  - Chapter 11: Eigensystems (for SCF eigenvalue solver)
  - Chapter 19: Partial Differential Equations (for understanding DIIS)

- **Golub, G. H.; Van Loan, C. F.** *Matrix Computations*, 4th ed.; Johns Hopkins University Press, 2013.
  - Chapter 8: Symmetric eigenvalue problems
  - Löwdin orthogonalization and generalized eigenvalue problems

#### SCF Convergence Acceleration
- **Pulay, P.** "Convergence acceleration of iterative sequences. The case of SCF iteration" *Chem. Phys. Lett.* **1980**, 73 (2), 393-398. (Original DIIS paper)
- **Kudin, K. N.; Scuseria, G. E.; Cancès, E.** "A black-box self-consistent field convergence algorithm: One step closer" *J. Chem. Phys.* **2002**, 116 (19), 8255-8261. (Modern DIIS variants)

### Parameter Databases and Atomic Properties

#### Covalent Radii (for Coordination Numbers)
- **Pyykkö, P.; Atsumi, M.** "Molecular Single-Bond Covalent Radii for Elements 1–118" *Chem. Eur. J.* **2009**, 15 (1), 186-197. DOI: [10.1002/chem.200800987](https://doi.org/10.1002/chem.200800987)
- **Pyykkö, P.; Atsumi, M.** "Molecular Double-Bond Covalent Radii for Elements Li–E112" *Chem. Eur. J.* **2009**, 15 (46), 12770-12779.

#### Electronegativity Scales
- **Pauling, L.** *The Nature of the Chemical Bond*, 3rd ed.; Cornell University Press, 1960. (Classic Pauling scale)
- **Allred, A. L.; Rochow, E. G.** "A Scale of Electronegativity Based on Electrostatic Force" *J. Inorg. Nucl. Chem.* **1958**, 5 (4), 264-268.

### Implementation Resources

#### Source Code Repositories
- **TBLite** (GFN1/GFN2 reference): https://github.com/tblite/tblite
  - Documentation: https://tblite.readthedocs.io/
  - API: https://tblite.readthedocs.io/en/latest/api/

- **XTB** (Alternative GFN implementation): https://github.com/grimme-lab/xtb
  - Documentation: https://xtb-docs.readthedocs.io/
  - User guide: https://xtb-docs.readthedocs.io/en/latest/contents.html

- **MOPAC** (PMx reference): https://github.com/openmopac/mopac
  - Manual: http://openmopac.net/manual/
  - Parameters: See source code `src/Parameters/`

#### C++ Linear Algebra Libraries
- **Eigen**: https://eigen.tuxfamily.org/ (used in Curcuma)
  - Dense matrix operations, eigenvalue solvers
  - Documentation: https://eigen.tuxfamily.org/dox/

- **LAPACK**: http://www.netlib.org/lapack/ (backend for many solvers)
  - Generalized eigenvalue problems (DSYGV)

### Related Review Articles

#### Grimme Group Methods Overview
- **Grimme, S.; Hansen, A.; Brandenburg, J. G.; Bannwarth, C.** "Dispersion-Corrected Mean-Field Electronic Structure Methods" *Chem. Rev.* **2016**, 116 (9), 5105-5154. DOI: [10.1021/acs.chemrev.5b00533](https://doi.org/10.1021/acs.chemrev.5b00533)
  - Comprehensive review of DFT-D and semi-empirical methods
  - Historical context and method comparison

#### Semi-Empirical Method Benchmarks
- **Christensen, A. S.; Kubař, T.; Cui, Q.; Elstner, M.** "Semiempirical Quantum Mechanical Methods for Noncovalent Interactions for Chemical and Biochemical Applications" *Chem. Rev.* **2016**, 116 (9), 5301-5337.
  - Benchmarks for GFN, PM6, DFTB methods
  - Accuracy assessment for different properties

### Educational Resources for Quick Reference

#### Online Tutorials and Lectures
- **XTB Documentation**: https://xtb-docs.readthedocs.io/en/latest/
  - User-friendly explanation of GFN methods
  - Parameter meanings and settings

- **DFTB+ Tutorials**: https://dftbplus.org/tutorial/
  - Tight-binding theory explained for chemists
  - SCF convergence strategies

#### Video Lectures (if available)
- Search YouTube for "Stefan Grimme lectures" (semi-empirical methods overview)
- "Tight-binding methods" lectures from computational chemistry courses

### Parameter Extraction Tools

**For reference when implementing parameters**:

1. **TBLite TOML files**:
   - `tblite/param/gfn2.toml` — Complete GFN2 parameters
   - `tblite/param/gfn1.toml` — Complete GFN1 parameters
   - Format: Human-readable TOML, easy to parse

2. **MOPAC Parameter Files**:
   - `mopac/src/Parameters/PM3_parameters.F90`
   - `mopac/src/Parameters/AM1_parameters.F90`
   - Format: Fortran DATA statements

3. **Conversion Scripts** (create if needed):
   - Python script to convert TOML → C++ header
   - Validation against TBLite calculations

---

## Appendix: Quick Theoretical Lookup Tables

### Common Physical Constants (CODATA 2018)

```cpp
// Already defined in src/core/units.h
constexpr double au2Angstrom = 0.529177210903;  // Bohr to Ångström
constexpr double au2eV = 27.21138602;           // Hartree to eV
constexpr double au2kcal = 627.509474;          // Hartree to kcal/mol
constexpr double eV2Hartree = 1.0 / 27.21138602;
```

### Typical Parameter Ranges (for validation)

| Parameter | Physical Meaning | Typical Range | Units |
|-----------|------------------|---------------|-------|
| **ζ (zeta)** | Slater exponent | 1.0 – 10.0 | a.u.⁻¹ |
| **U (Hubbard)** | On-site repulsion | 0.2 – 0.8 | Hartree |
| **IP (VSIP)** | Ionization potential | 5 – 25 | eV |
| **β (beta)** | Resonance integral | -20 – -5 | eV |
| **α (alpha)** | Repulsion exponent | 1.0 – 3.0 | a.u.⁻¹ |
| **CN** | Coordination number | 0 – 12 | dimensionless |

### Equation Number Cross-Reference

**GFN2 Paper (Bannwarth et al. 2019, JCTC)**:
- Eq. 1: Total energy expression
- Eq. 2: Electronic energy (tight-binding)
- Eq. 4: Coordination number (exponential counting)
- Eq. 7: Coulomb interaction (ES2)
- Eq. 9: Third-order onsite (ES3)
- Eq. 12: Hamiltonian matrix elements
- Eq. 15: D4 dispersion contribution

**GFN1 Paper (Grimme et al. 2017, JCTC)**:
- Eq. 3: Total energy
- Eq. 5: Coordination number
- Eq. 8: Hamiltonian off-diagonal
- Eq. 11: Halogen bond correction

---

**Document Status**: Complete algorithm analysis, attribution guide, and theoretical reference
**Last Updated**: 2025-11-11
**Next Steps**:
1. Begin Phase 1 GFN2 implementation with proper attribution
2. Create parameter validation suite against TBLite
3. Implement D3/D4 stub functions for testing without dispersion
