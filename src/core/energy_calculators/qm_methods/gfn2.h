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
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "src/core/global.h"
#include "gfn2-xtb_param.hpp"
#include "gfn2_params_loader.h"
#include "qm_driver.h"
#include "STOIntegrals.hpp"
#include "src/core/solvation/gbsa.h"

#include "src/core/energy_calculators/ff_methods/d4param_generator.h"

#ifdef USE_D4
#include "dftd4interface.h"
#endif

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <memory>
#include <vector>

/**
 * @brief Native implementation of GFN2-xTB (Extended Tight-Binding GFN2)
 *
 * Theoretical Background:
 *   GFN2-xTB is a semi-empirical quantum mechanical method based on tight-binding
 *   density functional theory. It provides excellent accuracy for molecular geometries,
 *   vibrational frequencies, and noncovalent interactions at computational cost
 *   ~1000x lower than DFT.
 *
 * Core Approximations:
 *   - Minimal valence basis (Slater-type orbitals)
 *   - Parameterized Hamiltonian matrix elements
 *   - Self-consistent charge (SCC) treatment via electrostatics
 *   - D4 dispersion correction (external, see dispersion stub)
 *   - Third-order charge corrections for improved accuracy
 *
 * Energy Expression:
 *   E_total = E_electronic + E_repulsion + E_coulomb + E_dispersion
 *
 * Reference:
 *   C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 2019, 15, 1652
 *   Equations referenced: Eq. 1 (total energy), Eq. 2 (electronic), Eq. 4 (CN),
 *                         Eq. 7 (Coulomb), Eq. 9 (3rd-order), Eq. 12 (Hamiltonian)
 *
 * Implementation Philosophy:
 *   - Educational transparency: algorithms clearly visible
 *   - Theoretical chemistry accessibility: formulas with physical meaning
 *   - Modular design: each energy component separable for testing
 *   - Performance: vectorized operations where possible
 *
 * Claude Generated: Full native GFN2 implementation based on TBLite source analysis
 */
class GFN2 : public QMDriver {
public:
    /**
     * @brief Default constructor
     * Initializes GFN2 with default parameters from ArrayParameters database
     */
    GFN2();

    /**
     * @brief Constructor with custom parameters (for testing/research)
     * @param params Custom GFN2 parameter set
     */
    explicit GFN2(const ArrayParameters& params);

    /**
     * @brief Destructor
     */
    virtual ~GFN2() = default;

    // =================================================================================
    // QMDriver Interface Implementation
    // =================================================================================

    /**
     * @brief Initialize molecule structure for GFN2 calculation
     *
     * Sets up:
     *   - Atomic positions and types
     *   - Charge and spin state
     *   - Basis set construction
     *   - Coordination numbers
     *
     * @return true if initialization successful
     */
    virtual bool InitialiseMolecule() override;

    /**
     * @brief Perform GFN2 energy and gradient calculation
     *
     * Algorithm:
     *   1. Calculate coordination numbers (Eq. 4 in ref.)
     *   2. Build overlap matrix S
     *   3. Build Hamiltonian matrix H₀ (Eq. 12)
     *   4. SCF convergence loop (self-consistent charges)
     *   5. Calculate repulsion energy (pairwise exponential)
     *   6. Calculate Coulomb energy (ES2 + ES3 + AES2, Eq. 7-11)
     *   7. Calculate dispersion energy (D4 stub, Eq. 15)
     *   8. Optionally: calculate analytical gradients
     *
     * @param gradient If true, calculate analytical gradients
     * @return Total energy in Hartree
     */
    virtual double Calculation(bool gradient = false) override;

    // =================================================================================
    // GFN2-Specific Methods
    // =================================================================================

    /**
     * @brief Get current coordination numbers for all atoms
     *
     * Coordination numbers are calculated using exponential counting function:
     *   CN_i = ∑_{j≠i} [1 / (1 + exp(-k₁(R_cov,ij/R_ij - 1)))]^(k₂)
     *
     * Parameters (from GFN2 paper):
     *   k₁ = 16.0    (steepness)
     *   k₂ = 4/3     (range decay)
     *
     * @return Vector of coordination numbers (size = natoms)
     */
    Vector getCoordinationNumbers() const { return m_coordination_numbers; }

    /**
     * @brief Get energy decomposition for analysis
     *
     * Useful for:
     *   - Debugging energy components
     *   - Understanding contributions
     *   - Comparing with TBLite reference
     *
     * @return JSON object with energy breakdown
     */
    json getEnergyDecomposition() const;

    /**
     * @brief Get HOMO-LUMO gap
     * @return Gap in eV
     */
    double getHOMOLUMOGap() const;

    /**
     * @brief Get HOMO energy
     * @return HOMO energy in eV
     */
    double getHOMOEnergy() const;

    /**
     * @brief Get LUMO energy
     * @return LUMO energy in eV
     */
    double getLUMOEnergy() const;

    /**
     * @brief Get partial atomic charges (Mulliken)
     * @return Vector of atomic charges (size = natoms)
     */
    Vector getPartialCharges() const { return m_charges; }

    /**
     * @brief Enable implicit solvation with GBSA model
     *
     * Activates native GBSA (Generalized Born + Surface Area) solvation.
     * Must be called before Calculation() to include solvation effects.
     *
     * @param solvent Solvent name (e.g., "water", "dmso", "acetone", "none")
     *
     * Example:
     *   GFN2 gfn2;
     *   gfn2.InitialiseMolecule(molecule);
     *   gfn2.setSolvent("water");  // Enable water solvation
     *   double energy = gfn2.Calculation();
     *
     * @note Claude Generated: Native GBSA integration (November 2025)
     */
    void setSolvent(const std::string& solvent) {
        m_solvent = solvent;
        if (solvent != "none" && solvent != "" && solvent != "vacuum") {
            m_solvation = std::make_unique<Curcuma::Solvation::GBSA>(solvent, "GFN2");
        } else {
            m_solvation.reset();
        }
    }

private:
    // =================================================================================
    // Basis Set Construction
    // =================================================================================

    /**
     * @brief Construct minimal valence basis set
     *
     * GFN2 uses Slater-type orbitals (STO) for each valence shell:
     *   - s-shell: 1 function
     *   - p-shell: 3 functions (px, py, pz)
     *   - d-shell: 5 functions (d-orbitals)
     *
     * STOs are contracted to Gaussians (STO-6G) for integral evaluation.
     *
     * Reference: TBLite `add_basis` function in gfn2.f90
     *
     * @return Number of basis functions created
     */
    int buildBasisSet();

    /**
     * @brief Build overlap matrix S
     *
     * Overlap integral between STOs:
     *   S_ij = <φ_i | φ_j>
     *
     * Uses analytical STO overlap formulas from STOIntegrals.hpp
     *
     * @param basisset Basis function descriptors
     * @return Overlap matrix (nbasis × nbasis, symmetric, positive definite)
     */
    Matrix MakeOverlap(std::vector<STO::Orbital>& basisset) override;

    // =================================================================================
    // Hamiltonian Construction
    // =================================================================================

    /**
     * @brief Build core Hamiltonian matrix H₀
     *
     * Diagonal elements (self-energy, GFN2 Eq. 12a):
     *   H_ii = E_shell(element, shell) + k_CN * CN_i
     *
     * Off-diagonal elements (hopping integral, GFN2 Eq. 12b):
     *   H_ij = z_ij * k_pair * k_shell * (1 + 0.02 * ΔEN²) * poly(r_ij) * S_ij
     *
     * where:
     *   z_ij = (2√(ζ_i·ζ_j)/(ζ_i+ζ_j))^0.5   orbital overlap scaling
     *   k_pair = pairwise coupling constant
     *   k_shell = shell-dependent coupling
     *   ΔEN = electronegativity difference (Pauling scale)
     *   poly(r) = distance-dependent polynomial correction
     *   S_ij = overlap integral
     *
     * Reference: TBLite `gfn2_h0spec` derived type, functions:
     *   - get_selfenergy
     *   - get_hscale
     *   - get_cnshift
     *   - get_shpoly
     *
     * @param S Overlap matrix
     * @param basisset Basis function descriptors
     * @return Hamiltonian matrix (nbasis × nbasis, symmetric)
     */
    Matrix MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset) override;

    /**
     * @brief Calculate atomic dipole moments from density matrix
     * @return Vector of atomic dipoles (natoms x 3)
     *
     * Claude Generated: Extracted from calculateGradient for modularity
     */
    std::vector<Eigen::Vector3d> calculateAtomicDipoles() const;

    /**
     * @brief Calculate atomic quadrupole moments from density matrix
     * @return Vector of atomic quadrupoles (natoms x 6 in xx, xy, xz, yy, yz, zz order)
     *
     * Claude Generated: Based on TBLite multipole.f90 implementation
     */
    std::vector<Eigen::VectorXd> calculateAtomicQuadrupoles() const;

    /**
     * @brief Calculate AES2 anisotropic electrostatic energy
     * @param dipoles Atomic dipole moments
     * @param quadrupoles Atomic quadrupole moments
     * @return AES2 energy contribution
     *
     * Reference: TBLite multipole.f90 get_energy (lines 188-215)
     */
    double calculateAES2Energy(const std::vector<Eigen::Vector3d>& dipoles,
                               const std::vector<Eigen::VectorXd>& quadrupoles) const;

    /**
     * @brief Calculate AES2 charge-dipole interaction matrix
     * @param rad Multipole damping radii
     * @param kdmp3 Damping exponent for inverse quadratic
     * @return Matrix (3 x natoms x natoms) for charge-dipole interactions
     *
     * Reference: TBLite multipole.f90 get_multipole_matrix_0d
     */
    Matrix calculateAES2ChargeDipoleMatrix(const Vector& rad, double kdmp3) const;

    /**
     * @brief Calculate AES2 dipole-dipole interaction matrix
     * @param rad Multipole damping radii
     * @param kdmp5 Damping exponent for inverse cubic
     * @return Vector of matrices for dipole-dipole interactions (9 x natoms x natoms)
     *
     * Note: Implementation uses vector of matrices for 4D dipole-dipole interactions
     * Reference: TBLite multipole.f90 get_multipole_matrix_0d
     */
    std::vector<Matrix> calculateAES2DipoleDipoleMatrix(const Vector& rad, double kdmp5) const;

    /**
     * @brief Calculate AES2 charge-quadrupole interaction matrix
     * @param rad Multipole damping radii
     * @param kdmp5 Damping exponent for inverse cubic
     * @return Matrix (6 x natoms x natoms) for charge-quadrupole interactions
     *
     * Reference: TBLite multipole.f90 get_multipole_matrix_0d
     */
    Matrix calculateAES2ChargeQuadrupoleMatrix(const Vector& rad, double kdmp5) const;

    /**
     * @brief Calculate self-energy for given shell
     *
     * Theoretical Background:
     *   Self-energy represents the energy of an electron in an isolated atomic orbital.
     *   In tight-binding theory, this is parameterized from valence state ionization
     *   potentials (VSIP).
     *
     * GFN2 Enhancement:
     *   Self-energy is shifted by coordination number to account for chemical environment:
     *   E_ii = E_base + k_CN * CN_i
     *
     * Parameters (from ArrayParameters):
     *   - selfenergy[element][shell]: base orbital energy (eV)
     *   - kcn[element][shell]: CN shift parameter
     *
     * Units: Returns energy in Hartree (converted from eV)
     *
     * @param element Atomic number
     * @param shell Shell index (0=s, 1=p, 2=d)
     * @param CN Coordination number of atom
     * @return Self-energy in Hartree
     */
    double getSelfEnergy(int element, int shell, double CN) const;

    /**
     * @brief Get shell-resolved Hubbard parameter for Coulomb calculations
     *
     * The Hubbard parameter γ (gamma) represents the chemical hardness or
     * self-Coulomb repulsion of an electron in a specific shell. GFN2 uses
     * shell-resolved Hubbard parameters:
     *
     * Formula: γ_shell = HUBBARD_PARAMETER[Z] * SHELL_HUBBARD_CORR[Z][shell]
     *
     * Parameters from TBLite gfn2.f90:
     *   - HUBBARD_PARAMETER[Z]: Base chemical hardness (86 elements)
     *   - SHELL_HUBBARD_CORR[Z][shell]: Shell-specific correction (s, p, d)
     *
     * Reference: TBLite gfn2.f90, effective.f90 (Coulomb kernel)
     *
     * Claude Generated (January 2025): Shell-resolved Hubbard for GFN2 SCC
     *
     * @param Z Atomic number (1-86)
     * @param shell Shell index (0=s, 1=p, 2=d)
     * @return Shell-resolved Hubbard parameter in Hartree
     */
    double getShellHubbard(int Z, int shell) const;

    /**
     * @brief Calculate Mulliken atomic charges from density matrix
     *
     * Mulliken population analysis: q_A = Z_A - Σ_{μ∈A} (P·S)_μμ
     *
     * Claude Generated (January 2025)
     *
     * @param density Density matrix P
     * @return Vector of atomic charges (size = natoms)
     */
    Vector calculateMullikenCharges(const Matrix& density) const;

    /**
     * @brief Get valence electron count for GFN2 minimal basis
     *
     * Returns the number of valence electrons for an element based on its
     * position in the periodic table. This is used for Mulliken charge analysis
     * in GFN2's minimal valence basis where core electrons are not included.
     *
     * Claude Generated (January 2025)
     *
     * @param Z Atomic number (1-86)
     * @return Number of valence electrons
     */
    int getValenceElectronCount(int Z) const;

    /**
     * @brief Calculate Hamiltonian scaling factor for off-diagonal element
     *
     * Implements GFN2 Eq. 12b hopping integral formula.
     *
     * @param fi Basis function i
     * @param fj Basis function j
     * @param distance Interatomic distance in Bohr
     * @return Scaling factor (dimensionless)
     */
    double getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const;

    // =================================================================================
    // Coordination Numbers
    // =================================================================================

    /**
     * @brief Calculate exponential coordination number (ECP)
     *
     * Theoretical Background:
     *   The coordination number CN_i represents the number of atoms bonded to atom i.
     *   The exponential counting function provides a continuous, differentiable measure
     *   suitable for gradient-based geometry optimization.
     *
     * Formula (GFN2 Eq. 4):
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
     * @return Vector of coordination numbers (dimensionless, typically 0-12)
     */
    Vector calculateCoordinationNumbers();

    // =================================================================================
    // SCF Procedure
    // =================================================================================

    /**
     * @brief Run self-consistent field (SCF) convergence loop
     *
     * Algorithm:
     *   1. Initial guess: superposition of atomic densities (SAD)
     *   2. Build Fock matrix F = H₀ + G(P)
     *      where G(P) contains Coulomb contributions
     *   3. Solve generalized eigenvalue problem: FC = SCE
     *   4. Build density matrix from occupied orbitals
     *   5. Check convergence: ||P_new - P_old|| < threshold
     *   6. Apply damping: P = α·P_new + (1-α)·P_old
     *   7. Repeat until converged
     *
     * Convergence criteria:
     *   - Density matrix change: ||ΔP|| < 1.0e-6
     *   - Energy change: |ΔE| < 1.0e-7 Hartree
     *   - Maximum iterations: 100
     *
     * Reference: Standard SCF procedure (e.g., Szabo/Ostlund, Modern Quantum Chemistry)
     *
     * @return true if SCF converged
     */
    bool runSCF();

    /**
     * @brief Build Fock matrix from current density
     *
     * F = H₀ + G(P)
     *
     * where G(P) contains:
     *   - Coulomb repulsion (shell-resolved)
     *   - Exchange correlation (implicit in tight-binding parameterization)
     *
     * @param density Current density matrix
     * @return Fock matrix (nbasis × nbasis)
     */
    Matrix buildFockMatrix(const Matrix& density);

    /**
     * @brief Calculate Coulomb kernel γ_AB using Klopman-Ohno-Mataga-Nishimoto
     *
     * Formula (TBLite effective.f90:271-277):
     *   γ_AB = 1 / [R_AB^gexp + (1/γ_avg)^gexp]^(1/gexp)
     *   where γ_avg = (γ_A + γ_B) / 2, gexp = 2.0 for GFN2
     *
     * Claude Generated (January 2025)
     *
     * @param gamma_A Hubbard parameter for atom A (shell-resolved)
     * @param gamma_B Hubbard parameter for atom B (shell-resolved)
     * @param R_AB Distance between atoms in Bohr
     * @return Coulomb kernel γ_AB
     */
    double calculateCoulombKernel(double gamma_A, double gamma_B, double R_AB) const;

    /**
     * @brief Build density matrix from molecular orbital coefficients
     *
     * Closed-shell formula:
     *   P_μν = 2 * ∑_{i=1}^{N_occ} C_μi * C_νi
     *
     * where N_occ = N_electrons / 2 (number of doubly occupied orbitals)
     *
     * @param mo_coefficients MO coefficient matrix (nbasis × nbasis)
     * @param mo_energies Orbital energies (nbasis)
     * @return Density matrix (nbasis × nbasis, symmetric)
     */
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies);

    // =================================================================================
    // Energy Components
    // =================================================================================

    /**
     * @brief Calculate electronic energy from density and Fock matrices
     *
     * E_elec = ∑_μν P_μν * (H_μν + F_μν) / 2
     *
     * @return Electronic energy in Hartree
     */
    double calculateElectronicEnergy() const;

    /**
     * @brief Calculate pairwise repulsion energy
     *
     * Repulsion between atomic cores (GFN2 uses exponential form):
     *   E_rep = ∑_{A<B} Z_eff,A * exp(-α_A * R_AB) + Z_eff,B * exp(-α_B * R_AB)
     *
     * Parameters (element-specific, from ArrayParameters):
     *   - rep_alpha[element]: exponential decay parameter
     *   - rep_zeff[element]: effective nuclear charge
     *
     * Units: Returns energy in Hartree
     *
     * @return Repulsion energy in Hartree
     */
    double calculateRepulsionEnergy() const;

    /**
     * @brief Calculate Coulomb interaction energy (ES2 + ES3 + AES2)
     *
     * Three components (GFN2 Eq. 7-11):
     *
     * 1. ES2 (effective Coulomb, Eq. 7):
     *    E_ES2 = ∑_{A<B} q_A * q_B * γ_AB(R_AB)
     *    where γ_AB is shell-resolved Coulomb kernel
     *
     * 2. ES3 (third-order onsite, Eq. 9):
     *    E_ES3 = (1/6) * ∑_A dU_A * q_A³
     *    where dU is Hubbard derivative
     *
     * 3. AES2 (damped multipole, Eq. 10-11):
     *    E_AES2 includes dipole and quadrupole contributions
     *    with Tang-Toennies damping
     *
     * Reference: GFN2 paper section 2.2 "Electrostatics"
     *
     * @return Coulomb energy in Hartree
     */
    double calculateCoulombEnergy() const;

    /**
     * @brief Calculate dispersion energy (D4)
     *
     * GFN2 uses DFT-D4 dispersion correction with parameters:
     *   s6 = 1.0, s8 = 2.7, a1 = 0.52, a2 = 5.0
     *
     * Claude Generated (March 2026): Integrated D4 dispersion
     * - When USE_D4 is defined: Full D4 calculation via DFTD4Interface
     * - When USE_D4 is not defined: Returns 0.0 with warning
     * - Expected contribution: ~0.01-0.1 Hartree for organic molecules
     *
     * Reference: E. Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122
     *
     * @return Dispersion energy in Hartree
     */
    double calculateDispersionEnergy() const;

    // =================================================================================
    // Gradient Calculation
    // =================================================================================

    /**
     * @brief Calculate analytical gradients (force = -gradient)
     *
     * Energy gradient with respect to nuclear coordinates:
     *   dE/dR_A = dE_elec/dR_A + dE_rep/dR_A + dE_disp/dR_A
     *
     * Uses Hellmann-Feynman theorem for electronic part.
     *
     * @return Gradient matrix (natoms × 3) in Hartree/Bohr
     */
    Matrix calculateGradient() const;

    /**
     * @brief Calculate numerical gradient using finite differences (for testing)
     *
     * Uses central finite differences with step size h:
     *   dE/dx ≈ (E(x+h) - E(x-h)) / (2h)
     *
     * @param h Step size in Angstrom (default: 0.001 Å)
     * @return Numerical gradient matrix (natoms × 3) in Hartree/Bohr
     *
     * Claude Generated (March 2026): For validating analytical gradients
     */
    Matrix calculateNumericalGradient(double h = 0.001) const;

    // =================================================================================
    // Data Members
    // =================================================================================

    ArrayParameters m_params;                  ///< Legacy GFN2 parameters (for compatibility)
    GFN2Params::ParameterDatabase m_param_db;  ///< Shell-resolved parameter database (November 2025)
    std::vector<STO::Orbital> m_basis;        ///< Basis set (STOs)
    int m_nbasis;                              ///< Number of basis functions

    Vector m_coordination_numbers;             ///< Coordination numbers (natoms)
    Vector m_charges;                          ///< Mulliken partial charges (natoms)

    Matrix m_overlap;                          ///< Overlap matrix S (nbasis × nbasis)
    Matrix m_hamiltonian;                      ///< Core Hamiltonian H₀ (nbasis × nbasis)
    Matrix m_fock;                             ///< Fock matrix F (nbasis × nbasis)
    Matrix m_density;                          ///< Density matrix P (nbasis × nbasis)
    Matrix m_mo;                               ///< Molecular orbital coefficients (nbasis × nbasis)
    Vector m_energies;                         ///< Orbital energies (nbasis)
    Matrix m_S_inv_sqrt;                       ///< Cached S^(-1/2) for SCF (computed once)

    // Energy components (for decomposition analysis)
    double m_energy_electronic;                ///< Electronic energy
    double m_energy_repulsion;                 ///< Core-core repulsion
    double m_energy_coulomb;                   ///< Coulomb interaction
    double m_energy_aes2;                      ///< Anisotropic electrostatic energy (AES2)
    double m_energy_dispersion;                ///< D4 dispersion (stub)
    double m_energy_solvation;                 ///< GBSA solvation energy (native)

    // Solvation model (optional, native implementation)
    std::unique_ptr<Curcuma::Solvation::GBSA> m_solvation;  ///< GBSA solvation model
    std::string m_solvent;                     ///< Solvent name (e.g., "water", "none")

#ifdef USE_D4
    std::unique_ptr<DFTD4Interface> m_d4;      ///< D4 dispersion calculator (external cpp-d4)
#endif
    mutable std::unique_ptr<D4ParameterGenerator> m_d4_native;  ///< Native D4 fallback (CN-only weighting)

    // SCF convergence parameters
    int m_scf_max_iterations;                  ///< Maximum SCF iterations (default: 100)
    double m_scf_threshold;                    ///< Convergence threshold (default: 1.0e-6)
    double m_scf_damping;                      ///< Density damping factor (default: 0.4)

    bool m_scf_converged;                      ///< SCF convergence status
};
