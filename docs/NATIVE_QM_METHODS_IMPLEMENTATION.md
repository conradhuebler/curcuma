# Native QM Methods Implementation Guide

**Document Purpose**: Comprehensive implementation guide for native GFN1, GFN2, and PMx semi-empirical quantum mechanical methods in Curcuma.

**Created**: 2025-11-11
**Based on**: TBLite source code analysis, MOPAC repository structure, and existing Curcuma EHT implementation

---

## Table of Contents

1. [GFN2-xTB Implementation](#gfn2-xtb-implementation)
2. [GFN1-xTB Implementation](#gfn1-xtb-implementation)
3. [PMx Methods Implementation](#pmx-methods-implementation)
4. [Common Infrastructure](#common-infrastructure)
5. [Implementation Roadmap](#implementation-roadmap)

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

## References

### GFN2-xTB
- Bannwarth, C.; Ehlert, S.; Grimme, S. *J. Chem. Theory Comput.* **2019**, 15, 1652-1671. DOI: 10.1021/acs.jctc.8b01176
- TBLite repository: https://github.com/tblite/tblite

### GFN1-xTB
- Grimme, S.; Bannwarth, C.; Shushkov, P. *J. Chem. Theory Comput.* **2017**, 13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118

### PMx Methods
- Stewart, J. J. P. *J. Mol. Model.* **2007**, 13, 1173-1213. (PM6)
- Stewart, J. J. P. *J. Comput. Chem.* **1989**, 10, 209-220. (PM3)
- Dewar, M. J. S.; et al. *J. Am. Chem. Soc.* **1985**, 107, 3902-3909. (AM1)
- MOPAC repository: https://github.com/openmopac/mopac

### DFT-D3/D4
- Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. *J. Chem. Phys.* **2010**, 132, 154104. (D3)
- Caldeweyher, E.; et al. *J. Chem. Phys.* **2019**, 150, 154122. (D4)

---

**Document Status**: Complete algorithm analysis and implementation guide
**Next Steps**: Begin Phase 1 GFN2 implementation or request clarifications
