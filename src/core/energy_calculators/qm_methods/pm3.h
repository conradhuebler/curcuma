/*
 * <Native PM3 Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the PM3 method developed by:
 *   James J. P. Stewart
 *   Stewart Computational Chemistry
 *
 * Reference implementation: MOPAC (http://openmopac.net/)
 *   Copyright (C) James J. P. Stewart
 *   Licensed under LGPL-3.0
 *
 * Original method publication:
 *   J. J. P. Stewart
 *   J. Comput. Chem. 1989, 10, 209-220
 *   DOI: 10.1002/jcc.540100208
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
 * @brief Native implementation of PM3 (Parametric Method 3)
 *
 * Theoretical Background:
 *   PM3 is a semi-empirical quantum mechanical method based on the NDDO
 *   (Neglect of Diatomic Differential Overlap) approximation. It uses
 *   minimal valence basis sets and neglects three- and four-center integrals.
 *
 * NDDO Approximation:
 *   - Only two-center integrals retained
 *   - Core-core repulsion with Gaussian expansions
 *   - Resonance integrals from Slater-Condon parameters
 *   - Element-specific parametrization for thermochemistry
 *
 * Energy Expression:
 *   E_total = E_electronic + E_core-core
 *   E_elec = ∑_μν P_μν H_μν + (1/2) ∑_μνλσ P_μν P_λσ (μν|λσ)
 *
 * PM3 vs Other Methods:
 *   - PM3: Optimized for organic molecules, heats of formation
 *   - AM1: Alternative parameter set, similar accuracy
 *   - MNDO: Predecessor, less accurate for H-bonds
 *   - PM6: Modern successor, broader element coverage
 *
 * Reference:
 *   J. J. P. Stewart, J. Comput. Chem. 1989, 10, 209
 *   MOPAC manual: http://openmopac.net/manual/
 *
 * Implementation Philosophy:
 *   - Educational transparency: NDDO approximations visible
 *   - Simplified parameter access (full PM3 has 100+ params per element)
 *   - Modular design for testing AM1, MNDO variants
 *
 * Claude Generated: Native PM3 implementation for Curcuma
 */
class PM3 : public QMDriver {
public:
    PM3();
    virtual ~PM3() = default;

    // QMDriver Interface
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false) override;

    // PM3-Specific Methods
    json getEnergyDecomposition() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getPartialCharges() const { return m_charges; }
    double getHeatOfFormation() const;  // PM3-specific: ΔH_f

private:
    // Basis set (minimal valence only)
    int buildBasisSet();
    Matrix MakeOverlap(std::vector<STO::Orbital>& basisset) override;
    Matrix MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset) override;

    // PM3 Hamiltonian construction
    double getCoreIntegral(int element, int orbital_type) const;
    double getResonanceIntegral(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const;

    // Two-electron integrals (NDDO approximation)
    double calculateTwoElectronIntegral(int mu, int nu, int lambda, int sigma) const;
    double getGammaAB(int Z_A, int Z_B, double R_AB) const;  // (ss|ss) integral

    // SCF procedure (similar to GFN methods)
    bool runSCF();
    Matrix buildFockMatrix(const Matrix& density);
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies);

    // Energy components
    double calculateElectronicEnergy() const;
    double calculateCoreRepulsionEnergy() const;  // PM3 Gaussian expansions

    // Gradients
    Matrix calculateGradient() const;

    // PM3 Parameters (simplified - real PM3 has many more)
    struct PM3Params {
        double U_ss;  // One-center, one-electron integral (s)
        double U_pp;  // One-center, one-electron integral (p)
        double beta_s;  // Resonance integral parameter (s)
        double beta_p;  // Resonance integral parameter (p)
        double zeta_s;  // Slater exponent (s)
        double zeta_p;  // Slater exponent (p)
        double alpha;  // Core repulsion parameter
        std::vector<double> gauss_a;  // Gaussian expansion coefficients
        std::vector<double> gauss_b;  // Gaussian expansion exponents
        std::vector<double> gauss_c;  // Gaussian expansion widths
    };

    std::map<int, PM3Params> m_pm3_params;  // Element-specific parameters
    void initializePM3Parameters();

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
