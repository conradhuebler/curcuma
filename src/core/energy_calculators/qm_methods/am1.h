/*
 * <Native AM1 Implementation>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the AM1 (Austin Model 1) method developed by:
 *   Michael J. S. Dewar and coworkers
 *   University of Texas at Austin, 1985
 *
 * Reference implementation: MOPAC (http://openmopac.net/)
 *   Copyright (C) James J. P. Stewart
 *   Licensed under LGPL-3.0
 *
 * Original method publication:
 *   M. J. S. Dewar, E. G. Zoebisch, E. F. Healy, J. J. P. Stewart
 *   J. Am. Chem. Soc. 1985, 107, 3902-3909
 *   DOI: 10.1021/ja00299a024
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
 * @brief Native implementation of AM1 (Austin Model 1)
 *
 * Theoretical Background:
 *   AM1 is an improved semi-empirical method based on MNDO with Gaussian
 *   core-core repulsion corrections. It addresses MNDO's over-repulsion of
 *   proximate atoms, improving hydrogen bonding and conformational energies.
 *
 * NDDO Approximation:
 *   - Same two-center integral framework as MNDO
 *   - Enhanced core-core repulsion with Gaussian expansions
 *   - Re-parameterized for improved accuracy over MNDO
 *   - Better treatment of hydrogen bonds and barriers
 *
 * Energy Expression:
 *   E_total = E_electronic + E_core-core
 *   E_elec = ∑_μν P_μν H_μν + (1/2) ∑_μνλσ P_μν P_λσ (μν|λσ)
 *   E_core = ∑_A<B [Z_A Z_B (ss|ss)_AB + V_Gauss_AB]
 *
 * Gaussian Core Correction:
 *   V_Gauss = Z_A Z_B (ss|ss) * ∑_k [a_k exp(-b_k (R - c_k)²)]
 *   Typically 2-4 Gaussian functions per element
 *
 * AM1 vs MNDO vs PM3:
 *   - MNDO (1977): Original, simple core repulsion
 *   - AM1 (1985): Gaussian core corrections, better H-bonds
 *   - PM3 (1989): Re-parameterized, thermochemistry focus
 *
 * Reference:
 *   M. J. S. Dewar et al., J. Am. Chem. Soc. 1985, 107, 3902
 *   MOPAC manual: http://openmopac.net/manual/
 *
 * Implementation Philosophy:
 *   - Educational transparency: NDDO approximations clearly visible
 *   - Foundation for understanding semi-empirical methods
 *   - Direct AM1 implementation with Gaussian corrections
 *
 * Claude Generated: Native AM1 implementation for Curcuma
 */
class AM1 : public QMDriver {
public:
    AM1();
    virtual ~AM1() = default;

    // QMDriver Interface
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false) override;

    // AM1-Specific Methods
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

    // AM1 Hamiltonian construction
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
    double calculateCoreRepulsionEnergy() const;  // AM1: Gaussian-corrected core repulsion

    // Gradients
    Matrix calculateGradient() const;

    // AM1 Parameters
    struct AM1Params {
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

    std::map<int, AM1Params> m_am1_params;  // Element-specific parameters
    void initializeAM1Parameters();

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
