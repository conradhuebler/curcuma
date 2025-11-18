/*
 * <GFN-FF GBSA Solvation Model Header. >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

/**
 * @file gbsa.h
 * @brief Generalized Born Surface Area (GBSA) implicit solvation model for GFN-FF
 *
 * Implementation of GBSA solvation model based on:
 * - Generalized Born approximation for electrostatic solvation
 * - Solvent Accessible Surface Area (SASA) for non-polar contribution
 * - Reference: gfnff Fortran implementation (src/gbsa/)
 *
 * Energy formula:
 *   E_solv = E_Born + E_SASA + E_shift
 *   E_Born = 0.5 * q^T * B * q (+ E_HB if hydrogen bonding enabled)
 *   E_SASA = gamma * SASA
 *
 * Where:
 *   - q = atomic charges (from EEQ)
 *   - B = Born matrix (from Born radii and dielectric kernel)
 *   - gamma = surface tension parameter (element-specific)
 *   - SASA = solvent accessible surface area
 *
 * Claude Generated (2025): C++ implementation following Fortran reference
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "src/tools/general.h"

namespace GBSA {

/**
 * @brief GBSA parameters for a specific solvent
 *
 * Reference: gfnff/include/param_gbsa_*.fh
 */
struct GBSAParameters {
    std::string solvent_name;

    // Solvent properties
    double dielectric_constant;  ///< Dielectric constant (epsilon)
    double molar_mass;            ///< Molar mass in g/mol
    double density;               ///< Density in g/cm³

    // Born radii parameters
    double born_scale;           ///< c1 - scaling factor for Born radii
    double probe_radius;         ///< rprobe - probe radius for SASA integration (Å)
    double born_offset;          ///< soset - offset parameter for Born radii

    // Energy shifts
    double free_energy_shift;    ///< gshift - free energy shift (kcal/mol)

    // ALPB parameter (0.0 for GBSA, non-zero for ALPB)
    double alpha;

    // Element-specific parameters (94 elements, indexed by atomic number - 1)
    std::vector<double> surface_tension;  ///< gamscale[94] - surface tension scaling
    std::vector<double> descreening;      ///< sx[94] - descreening parameters
    std::vector<double> hbond_strength;   ///< tmp[94] - hydrogen bond strength

    /**
     * @brief Default constructor with zeros
     */
    GBSAParameters()
        : dielectric_constant(0.0)
        , molar_mass(0.0)
        , density(0.0)
        , born_scale(0.0)
        , probe_radius(0.0)
        , born_offset(0.0)
        , free_energy_shift(0.0)
        , alpha(0.0)
        , surface_tension(94, 0.0)
        , descreening(94, 0.0)
        , hbond_strength(94, 0.0)
    {
    }
};

/**
 * @brief GBSA Solvation Model Class
 *
 * Implements generalized Born solvation with surface area correction.
 * Based on gfnff/src/gbsa/gbsa.f90 implementation.
 *
 * Claude Generated (2025): Educational implementation with clear scientific documentation
 */
class GBSASolvation {
public:
    /**
     * @brief Constructor
     * @param atomcount Number of atoms
     * @param atoms Atomic numbers (element types)
     * @param geometry Atomic coordinates (atomcount x 3)
     * @param solvent Solvent name ("water", "acetone", etc.)
     */
    GBSASolvation(int atomcount, const std::vector<int>& atoms,
        const Matrix& geometry, const std::string& solvent = "water");

    /**
     * @brief Calculate GBSA solvation energy
     * @param charges Atomic charges from EEQ (Vector of size atomcount)
     * @param gradient If true, also calculate gradients
     * @return Total solvation energy in Hartree
     *
     * Formula: E_solv = 0.5 * q^T * B * q + E_SASA + gshift
     */
    double calculateEnergy(const Vector& charges, bool gradient = false);

    /**
     * @brief Get gradient of solvation energy
     * @return Gradient matrix (atomcount x 3)
     */
    Matrix getGradient() const { return m_gradient; }

    /**
     * @brief Get Born radii
     * @return Vector of Born radii (Å)
     */
    Vector getBornRadii() const { return m_born_radii; }

    /**
     * @brief Get SASA values
     * @return Vector of atomic SASA values (Ų)
     */
    Vector getSASA() const { return m_sasa; }

    /**
     * @brief Get energy components
     * @param e_born Born electrostatic energy (output)
     * @param e_sasa Surface area energy (output)
     * @param e_hb Hydrogen bonding energy (output)
     * @param e_shift Free energy shift (output)
     */
    void getEnergyComponents(double& e_born, double& e_sasa,
        double& e_hb, double& e_shift) const;

    /**
     * @brief Update geometry for new calculation
     * @param geometry New atomic coordinates
     */
    void updateGeometry(const Matrix& geometry);

    /**
     * @brief Enable/disable hydrogen bonding correction
     * @param enable True to enable H-bond correction
     */
    void setHydrogenBondingCorrection(bool enable) { m_use_hbond = enable; }

private:
    // Molecular data
    int m_atomcount;
    std::vector<int> m_atoms;
    Matrix m_geometry;

    // Solvation parameters
    GBSAParameters m_params;
    std::string m_solvent;

    // Computed quantities
    Vector m_born_radii;          ///< Born radii for each atom
    Matrix m_born_matrix;         ///< Born interaction matrix
    Vector m_sasa;                ///< Solvent accessible surface area per atom
    Matrix m_gradient;            ///< Energy gradient
    double m_gsasa;               ///< Total SASA energy
    double m_gborn;               ///< Total Born energy
    double m_ghb;                 ///< Total H-bond energy

    // Van der Waals radii
    Vector m_vdw_radii;           ///< vdW radii for each atom
    Vector m_vdw_radii_offset;    ///< vdW radii with offset (svdw)
    Vector m_descreening_radii;   ///< Descreening radii (rho)

    // Flags
    bool m_use_hbond;             ///< Use hydrogen bonding correction
    bool m_initialized;           ///< Initialization flag

    /**
     * @brief Load GBSA parameters for specified solvent
     * @param solvent Solvent name
     * @return true if parameters loaded successfully
     */
    bool loadSolventParameters(const std::string& solvent);

    /**
     * @brief Calculate Born radii using GBOBC integrator
     *
     * Reference: gfnff/src/gbsa/born.f90 - compute_bornr()
     * Formula: Born radius calculation with tanh correction
     *   R_Born,i = c1 / (1/svdw_i - 1/vdw_i * tanh(arg))
     *   arg = psi * (alpha + psi * (gamma * psi - beta))
     *
     * GBOBCII parameters: alpha=1.0, beta=0.8, gamma=4.85
     */
    void calculateBornRadii();

    /**
     * @brief Calculate psi integral for Born radii
     *
     * Reference: gfnff/src/gbsa/born.f90 - compute_psi()
     * Psi integral accounts for neighbor overlap contributions
     */
    void calculatePsiIntegral();

    /**
     * @brief Construct Born interaction matrix
     *
     * Reference: gfnff/src/gbsa/kernel.f90 - addBornMatStill()
     * Formula: B_ij = f(r_ij, R_i, R_j, epsilon)
     *
     * Still kernel:
     *   f_ij = 1/sqrt(r_ij² + R_i*R_j*exp(-r_ij²/(4*R_i*R_j)))
     * with dielectric screening: (1 - 1/epsilon)
     */
    void constructBornMatrix();

    /**
     * @brief Calculate solvent accessible surface area
     *
     * Reference: gfnff/src/gbsa/sasa.f90 - compute_numsa()
     * Numerical integration using Lebedev angular grids
     */
    void calculateSASA();

    /**
     * @brief Get van der Waals radius for element
     * @param atomic_number Atomic number (1-94)
     * @return vdW radius in Å
     *
     * Reference: gfnff/src/gbsa/vdwrad.f90
     */
    double getVdWRadius(int atomic_number) const;

    /**
     * @brief Initialize van der Waals radii for all atoms
     */
    void initializeVdWRadii();
};

/**
 * @brief Get GBSA parameters for water (GFN2 parametrization)
 *
 * Reference: gfnff/include/param_gbsa_h2o.fh - gfn2_h2o
 */
GBSAParameters getWaterParameters();

/**
 * @brief Get GBSA parameters for acetone
 *
 * Reference: gfnff/include/param_gbsa_acetone.fh
 */
GBSAParameters getAcetoneParameters();

/**
 * @brief Get GBSA parameters for DMSO
 *
 * Reference: gfnff/include/param_gbsa_dmso.fh
 */
GBSAParameters getDMSOParameters();

/**
 * @brief Get GBSA parameters by solvent name
 * @param solvent Solvent name (case-insensitive)
 * @return GBSA parameters, or water parameters if not found
 */
GBSAParameters getSolventParameters(const std::string& solvent);

} // namespace GBSA
