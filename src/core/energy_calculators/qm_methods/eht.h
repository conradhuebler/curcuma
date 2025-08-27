/*
 * <Extended Hückel Theory Implementation in Curcuma>
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Extended Hückel Theory (EHT) is a semi-empirical quantum chemical method
 * developed by Roald Hoffmann. This implementation provides:
 * - Parameterized atomic orbital energies (VSIP values)
 * - Slater-type orbital overlap and Hamiltonian matrix construction
 * - Molecular orbital analysis and electronic structure calculations
 * - Support for common organic elements (H, C, N, O, F)
 */

#pragma once
#include "src/core/global.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <memory>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"
#include "eht_parameters.h"

#include "interface/abstract_interface.h"
#include "json.hpp"
#include "qm_driver.h"

typedef std::vector<STO::Orbital> Basisset;

/**
 * @brief Extended Hückel Theory (EHT) Implementation
 *
 * This class implements the Extended Hückel Theory method for semi-empirical
 * quantum chemical calculations. EHT is particularly useful for studying
 * molecular orbital interactions and electronic structure of organic molecules.
 *
 * Key features:
 * - Parameterized atomic orbital energies based on VSIP values
 * - Slater-type orbital basis functions
 * - Wolfsberg-Helmholz approximation for off-diagonal Hamiltonian elements
 * - Support for H, C, N, O, F atoms
 * - Molecular orbital analysis and orbital energy calculation
 *
 * References:
 * - Hoffmann, R. "An Extended Hückel Theory" J. Chem. Phys. 39, 1397 (1963)
 * - Hoffmann, R. "Solids and Surfaces: A Chemist's View of Bonding" (1988)
 */
class EHT : public QMDriver {
public:
    /**
     * @brief Default constructor
     * Initializes EHT with default Wolfsberg-Helmholz constant K=1.75
     */
    EHT();

    /**
     * @brief Constructor with custom Wolfsberg-Helmholz constant
     * @param K Wolfsberg-Helmholz constant (typically 1.75)
     */
    explicit EHT(double K);

    /**
     * @brief Initialize molecule from Mol object
     * @param molecule Molecule object containing geometry and atom types
     * @return true if initialization successful
     */
    // virtual bool InitialiseMolecule(const Mol& molecule) override;

    /**
     * @brief Initialize molecule with current stored molecule data
     * @return true if initialization successful
     */
    virtual bool InitialiseMolecule() override;

    /**
     * @brief Perform EHT calculation
     * @param gradient Currently not supported for EHT (will print warning)
     * @return Total electronic energy (currently returns 0, as EHT focuses on orbital energies)
     */
    virtual double Calculation(bool gradient = false) override;

    /**
     * @brief Get molecular orbital coefficients matrix
     * @return Matrix where columns represent molecular orbitals
     */
    Matrix MolecularOrbitals() const { return m_mo; }

    /**
     * @brief Get orbital energies in ascending order
     * @return Vector of orbital energies in eV
     */
    Vector Energies() const { return m_energies; }

    /**
     * @brief Get total number of valence electrons
     * @return Number of electrons in the system
     */
    int NumElectrons() const { return m_num_electrons; }

    /**
     * @brief Get HOMO-LUMO gap
     * @return Energy gap between HOMO and LUMO in eV
     */
    double getHOMOLUMOGap() const;

    /**
     * @brief Get HOMO energy
     * @return Highest Occupied Molecular Orbital energy in eV
     */
    double getHOMOEnergy() const;

    /**
     * @brief Get LUMO energy
     * @return Lowest Unoccupied Molecular Orbital energy in eV
     */
    double getLUMOEnergy() const;

    /**
     * @brief Set Wolfsberg-Helmholz constant
     * @param K New K value (typically between 1.5 and 2.0)
     */
    void setWolfsbergHelmholzConstant(double K) { m_K = K; }

    /**
     * @brief Get current Wolfsberg-Helmholz constant
     * @return Current K value
     */
    double getWolfsbergHelmholzConstant() const { return m_K; }

    /**
     * @brief Print orbital analysis summary
     * @param num_orbitals_around_gap Number of orbitals to show around HOMO-LUMO gap
     */
    void printOrbitalAnalysis(int num_orbitals_around_gap = 5) const;

    /**
     * @brief Print orbital analysis using CurcumaLogger for Level 2 verbosity
     * Structured output with HOMO/LUMO analysis and molecular properties
     * Claude Generated
     */
    void printOrbitalAnalysisVerbose() const;

private:
    /**
     * @brief Construct basis set from molecular geometry and atom types
     * @return Vector of STO::Orbital objects representing the basis set
     */
    Basisset MakeBasis();

    /**
     * @brief Construct overlap matrix S
     * @param basisset Basis set of atomic orbitals
     * @return Overlap matrix S_ij = <φ_i|φ_j>
     */
    Matrix MakeOverlap(Basisset& basisset) override;

    /**
     * @brief Construct Hamiltonian matrix using Wolfsberg-Helmholz approximation
     * @param S Overlap matrix
     * @param basisset Basis set of atomic orbitals
     * @return Hamiltonian matrix H_ij
     *
     * Uses the approximation:
     * H_ii = VSIP_i (Valence State Ionization Potential)
     * H_ij = K * S_ij * (VSIP_i + VSIP_j) / 2  for i ≠ j
     */
    Matrix MakeH(const Matrix& S, const Basisset& basisset) override;

    /**
     * @brief Validate that all atoms in molecule are supported
     * @return true if all atoms have EHT parameters available
     */
    bool validateMolecule() const;

    /**
     * @brief Calculate total electronic energy (sum of occupied orbital energies)
     * @return Total electronic energy in eV
     */
    double calculateElectronicEnergy() const;

private:
    double m_K; ///< Wolfsberg-Helmholz constant (typically 1.75)
};
