/*
 * <Native PM6 Implementation>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the PM6 method developed by:
 *   James J. P. Stewart
 *   Stewart Computational Chemistry, Colorado Springs, CO, USA, 2007
 *
 * Reference implementation: MOPAC (http://openmopac.net/)
 *   Copyright (C) James J. P. Stewart
 *   Licensed under LGPL-3.0
 *
 * Original method publication:
 *   J. J. P. Stewart
 *   J. Mol. Model. 2007, 13, 1173-1213
 *   DOI: 10.1007/s00894-007-0233-4
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
#include "qm_driver.h"
#include "STOIntegrals.hpp"

#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <map>

/**
 * @brief Native implementation of PM6 (Parametric Method 6)
 *
 * Theoretical Background:
 *   PM6 is the most modern and comprehensive NDDO-based semi-empirical method.
 *   It extends the PM3/AM1 framework with improved parameterization covering
 *   70+ elements (H through Pu) with high accuracy for organic and inorganic
 *   chemistry, biochemistry, and materials science.
 *
 * NDDO Approximation:
 *   - Same two-center integral framework as PM3/AM1/MNDO
 *   - Enhanced Gaussian core-core repulsion corrections
 *   - Extensively re-parameterized against large molecular database
 *   - Improved treatment of metals and heavy elements
 *
 * Energy Expression:
 *   E_total = E_electronic + E_core-core
 *   E_elec = ∑_μν P_μν H_μν + (1/2) ∑_μνλσ P_μν P_λσ (μν|λσ)
 *   E_core = ∑_A<B [Z_A Z_B (ss|ss)_AB + V_Gauss_AB]
 *
 * Gaussian Core Correction:
 *   V_Gauss = Z_A Z_B (ss|ss) * ∑_k [a_k exp(-b_k (R - c_k)²)]
 *   PM6 uses 2-4 Gaussian functions per element pair
 *
 * PM6 vs PM3 vs AM1 vs MNDO:
 *   - MNDO (1977): 4 elements, simple core repulsion
 *   - AM1 (1985): 12 elements, Gaussian corrections
 *   - PM3 (1989): 18 elements, thermochemistry focus
 *   - PM6 (2007): 70+ elements, broad applicability, best accuracy
 *
 * PM6 Improvements:
 *   - Extended element coverage (H-Pu)
 *   - Improved accuracy for heats of formation
 *   - Better geometries and barrier heights
 *   - Parameterized for solid-state and biological systems
 *   - Enhanced treatment of transition metals
 *
 * Reference:
 *   J. J. P. Stewart, J. Mol. Model. 2007, 13, 1173-1213
 *   MOPAC manual: http://openmopac.net/manual/
 *
 * Implementation Philosophy:
 *   - Educational transparency: NDDO approximations clearly visible
 *   - Most modern native semi-empirical method in Curcuma
 *   - Direct PM6 implementation with authentic parameters
 *
 * Claude Generated: Native PM6 implementation for Curcuma
 */
class PM6 : public QMDriver {
public:
    PM6();
    virtual ~PM6() = default;

    // QMDriver Interface
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false) override;

    // PM6-Specific Methods
    json getEnergyDecomposition() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getPartialCharges() const { return m_charges; }
    double getHeatOfFormation() const;

private:
    // Basis set (minimal valence only)
    int buildBasisSet();
    Matrix MakeOverlap(std::vector<STO::Orbital>& basisset) override;
    Matrix MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset) override;

    // PM6 Hamiltonian construction
    double getCoreIntegral(int element, int orbital_type) const;
    double getResonanceIntegral(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const;

    // Two-electron integrals (NDDO approximation with full MNDO multipole expansion)
    double calculateTwoElectronIntegral(int mu, int nu, int lambda, int sigma) const;
    double getGammaAB(int Z_A, int Z_B, double R_AB) const;  // (ss|ss) integral

    // SCF procedure
    bool runSCF();
    Matrix buildFockMatrix(const Matrix& density);
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies);

    // Energy components
    double calculateElectronicEnergy() const;
    double calculateCoreRepulsionEnergy() const;  // PM6: Enhanced Gaussian-corrected core repulsion

    // Gradients
    Matrix calculateGradient() const;

    // PM6 Parameters
    struct PM6Params {
        double U_ss;  // One-center, one-electron integral (s)
        double U_pp;  // One-center, one-electron integral (p)
        double beta_s;  // Resonance integral parameter (s)
        double beta_p;  // Resonance integral parameter (p)
        double zeta_s;  // Slater exponent (s)
        double zeta_p;  // Slater exponent (p)
        double alpha;  // Core repulsion parameter
        double D1;    // Dipole expansion parameter (for MNDO multipole integrals)
        double D2;    // Quadrupole expansion parameter (for d-orbitals, =0 for sp)
        double rho_s; // Orbital exponent for s-type ERIs (ρ_s)
        double rho_p; // Orbital exponent for p-type ERIs (ρ_p)
        std::vector<double> gauss_a;  // Gaussian expansion coefficients
        std::vector<double> gauss_b;  // Gaussian expansion exponents
        std::vector<double> gauss_c;  // Gaussian expansion centers
    };

    std::map<int, PM6Params> m_pm6_params;  // Element-specific parameters
    void initializePM6Parameters();

    // Data members
    std::vector<STO::Orbital> m_basis;
    int m_nbasis;

    Vector m_charges;

    Matrix m_overlap;
    Matrix m_hamiltonian;
    Matrix m_fock;
    Matrix m_density;

    double m_energy_electronic;
    double m_energy_core_repulsion;

    int m_scf_max_iterations;
    double m_scf_threshold;
    double m_scf_damping;
    bool m_scf_converged;
};
