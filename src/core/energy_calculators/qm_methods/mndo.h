/*
 * <Native MNDO Implementation>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the MNDO method developed by:
 *   Michael J. S. Dewar and Walter Thiel
 *   University of Texas at Austin
 *
 * Reference implementation: MOPAC (http://openmopac.net/)
 *   Copyright (C) James J. P. Stewart
 *   Licensed under LGPL-3.0
 *
 * Original method publications:
 *   M. J. S. Dewar, W. Thiel
 *   J. Am. Chem. Soc. 1977, 99, 4899-4907
 *   DOI: 10.1021/ja00457a004
 *
 *   M. J. S. Dewar, W. Thiel
 *   Theor. Chim. Acta 1977, 46, 89-104
 *   DOI: 10.1007/BF00548085
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
 * @brief Native implementation of MNDO (Modified Neglect of Diatomic Overlap)
 *
 * Theoretical Background:
 *   MNDO is the original semi-empirical quantum mechanical method based on the
 *   NDDO (Neglect of Diatomic Differential Overlap) approximation. It was the
 *   predecessor to AM1 and PM3.
 *
 * NDDO Approximation:
 *   - Only two-center integrals retained
 *   - Core-core repulsion without Gaussian expansions (unlike PM3/AM1)
 *   - Resonance integrals from Slater-Condon parameters
 *   - Element-specific parametrization for ground state properties
 *
 * Energy Expression:
 *   E_total = E_electronic + E_core-core
 *   E_elec = ∑_μν P_μν H_μν + (1/2) ∑_μνλσ P_μν P_λσ (μν|λσ)
 *   E_core = ∑_A<B (Z_A Z_B / R_AB + core-core repulsion terms)
 *
 * MNDO vs Later Methods:
 *   - MNDO (1977): Original parametrization, simple core repulsion
 *   - AM1 (1985): Improved with Gaussian core corrections
 *   - PM3 (1989): Re-parameterized for thermochemistry
 *   - PM6 (2007): Extended to more elements, better accuracy
 *
 * Reference:
 *   M. J. S. Dewar, W. Thiel, J. Am. Chem. Soc. 1977, 99, 4899
 *   MOPAC manual: http://openmopac.net/manual/
 *
 * Implementation Philosophy:
 *   - Educational transparency: NDDO approximations clearly visible
 *   - Direct MNDO implementation without Gaussian corrections
 *   - Foundation for AM1/PM3 implementations
 *
 * Claude Generated: Native MNDO implementation for Curcuma
 */
class MNDO : public QMDriver {
public:
    MNDO();
    virtual ~MNDO() = default;

    // QMDriver Interface
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false) override;

    // MNDO-Specific Methods
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

    // MNDO Hamiltonian construction
    double getCoreIntegral(int element, int orbital_type) const;
    double getResonanceIntegral(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const;

    // Two-electron integrals (NDDO approximation with full MNDO multipole expansion)
    double calculateTwoElectronIntegral(int mu, int nu, int lambda, int sigma) const;

    // SCF procedure
    bool runSCF();
    Matrix buildFockMatrix(const Matrix& density);
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies);

    // Energy components
    double calculateElectronicEnergy() const;
    double calculateCoreRepulsionEnergy() const;  // MNDO: Simple formula without Gaussians

    // Gradients
    Matrix calculateGradient() const;

    // MNDO Parameters
    struct MNDOParams {
        double U_ss;  // One-center, one-electron integral (s)
        double U_pp;  // One-center, one-electron integral (p)
        double beta_s;  // Resonance integral parameter (s)
        double beta_p;  // Resonance integral parameter (p)
        double zeta_s;  // Slater exponent (s)
        double zeta_p;  // Slater exponent (p)
        double alpha;  // Core repulsion parameter (MNDO-specific)
        double D1;    // Dipole expansion parameter (for MNDO multipole integrals)
        double D2;    // Quadrupole expansion parameter (for d-orbitals, =0 for sp)
        double rho_s; // Orbital exponent for s-type ERIs (ρ_s)
        double rho_p; // Orbital exponent for p-type ERIs (ρ_p)
        // MNDO does NOT use Gaussian expansions (unlike PM3/AM1)
    };

    std::map<int, MNDOParams> m_mndo_params;  // Element-specific parameters
    void initializeMNDOParameters();

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
