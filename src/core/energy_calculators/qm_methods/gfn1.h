/*
 * <Native GFN1-xTB Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the GFN1-xTB method developed by:
 *   Stefan Grimme, Christoph Bannwarth, Philip Shushkov
 *   Mulliken Center for Theoretical Chemistry, University of Bonn
 *
 * Reference implementation: TBLite (https://github.com/tblite/tblite)
 *   Copyright (C) 2019-2024 Sebastian Ehlert and contributors
 *   Licensed under LGPL-3.0-or-later
 *
 * Original method publication:
 *   S. Grimme, C. Bannwarth, P. Shushkov
 *   J. Chem. Theory Comput. 2017, 13, 1989-2009
 *   DOI: 10.1021/acs.jctc.7b00118
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
#include "gfn2-xtb_param.hpp"  // Reuse parameter structure
#include "qm_driver.h"
#include "STOIntegrals.hpp"

#include <Eigen/Dense>
#include <memory>
#include <vector>

/**
 * @brief Native implementation of GFN1-xTB (Extended Tight-Binding GFN1)
 *
 * Theoretical Background:
 *   GFN1-xTB is the predecessor of GFN2-xTB, offering a good balance between
 *   accuracy and computational efficiency. It's particularly well-suited for
 *   organometallic systems and includes explicit halogen bond corrections.
 *
 * Key Differences from GFN2:
 *   - D3 dispersion correction (instead of D4)
 *   - Halogen bond (XB) correction term
 *   - Simpler valence shell structure
 *   - Different parametrization philosophy
 *
 * Energy Expression:
 *   E_total = E_electronic + E_repulsion + E_coulomb + E_disp(D3) + E_XB
 *
 * Reference:
 *   S. Grimme et al., J. Chem. Theory Comput. 2017, 13, 1989
 *   TBLite implementation: src/tblite/xtb/gfn1.f90
 *
 * Implementation Philosophy:
 *   - Educational transparency: simplified compared to full TBLite
 *   - Chemically reasonable approximations for missing parameters
 *   - Modular design for easy testing and extension
 *
 * Claude Generated: Native GFN1 implementation for Curcuma
 */
class GFN1 : public QMDriver {
public:
    GFN1();
    explicit GFN1(const ArrayParameters& params);
    virtual ~GFN1() = default;

    // QMDriver Interface
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false) override;

    // GFN1-Specific Methods
    Vector getCoordinationNumbers() const { return m_coordination_numbers; }
    json getEnergyDecomposition() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getPartialCharges() const { return m_charges; }

private:
    // Basis set
    int buildBasisSet();
    Matrix MakeOverlap(const std::vector<STO::Orbital>& basisset) override;
    Matrix MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset) override;

    // Hamiltonian construction (similar to GFN2 but simpler)
    double getSelfEnergy(int element, int shell, double CN) const;
    double getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const;

    // Coordination numbers (same as GFN2)
    Vector calculateCoordinationNumbers();

    // SCF procedure
    bool runSCF();
    Matrix buildFockMatrix(const Matrix& density);
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies);

    // Energy components
    double calculateElectronicEnergy() const;
    double calculateRepulsionEnergy() const;
    double calculateCoulombEnergy() const;
    double calculateDispersionEnergy() const;  // D3 stub
    double calculateHalogenBondCorrection() const;  // GFN1-specific XB correction

    // Gradients
    Matrix calculateGradient() const;

    // Data members
    ArrayParameters m_params;
    std::vector<STO::Orbital> m_basis;
    int m_nbasis;

    Vector m_coordination_numbers;
    Vector m_charges;

    Matrix m_overlap;
    Matrix m_hamiltonian;
    Matrix m_fock;
    Matrix m_density;

    // Energy components
    double m_energy_electronic;
    double m_energy_repulsion;
    double m_energy_coulomb;
    double m_energy_dispersion;   // D3 stub
    double m_energy_halogen_bond;  // XB correction

    // SCF parameters
    int m_scf_max_iterations;
    double m_scf_threshold;
    double m_scf_damping;
    bool m_scf_converged;
};
