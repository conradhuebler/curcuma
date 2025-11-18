/*
 * <GFN-FF GBSA Solvation Model Implementation. >
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
 * @file gbsa.cpp
 * @brief Implementation of GBSA solvation model for GFN-FF
 *
 * Claude Generated (2025): C++ implementation following Fortran reference
 * Reference: gfnff/src/gbsa/ (Fortran implementation)
 */

#include "gbsa.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "src/core/curcuma_logger.h"

namespace GBSA {

// GBOBCII constants (from born.f90)
constexpr double ALPHA_GBOBC = 1.0;  ///< alp parameter
constexpr double BETA_GBOBC = 0.8;   ///< bet parameter
constexpr double GAMMA_GBOBC = 4.85; ///< gam parameter

// Van der Waals radii (Bondi radii in Å)
// Reference: gfnff/src/gbsa/vdwrad.f90 - vanDerWaalsRadD3
static const double VDW_RADII_D3[94] = {
    1.09155, 0.86735, 1.74780, 1.54910,  // H-Be
    1.60800, 1.45515, 1.31125, 1.24085,  // B-O
    1.14980, 1.06870, 1.85410, 1.74195,  // F-Mg
    2.00530, 1.89585, 1.75085, 1.65535,  // Al-S
    1.55230, 1.45740, 2.12055, 2.05175,  // Cl-Ca
    1.94515, 1.88210, 1.86055, 1.72070,  // Sc-Cr
    1.77310, 1.72105, 1.71635, 1.67310,  // Mn-Ni
    1.65040, 1.61545, 1.97895, 1.93095,  // Cu-Ge
    1.83125, 1.76340, 1.68310, 1.60480,  // As-Kr
    2.30880, 2.23820, 2.10980, 2.02985,  // Rb-Zr
    1.92980, 1.87715, 1.78450, 1.73115,  // Nb-Ru
    1.69875, 1.67625, 1.66540, 1.73100,  // Rh-Cd
    2.13115, 2.09370, 2.00750, 1.94505,  // In-Te
    1.86900, 1.79445, 2.52835, 2.59070,  // I-Ba
    2.31305, 2.31005, 2.28510, 2.26355,  // La-Nd
    2.24515, 2.22575, 2.21170, 2.06215,  // Pm-Gd
    2.12135, 2.07705, 2.13970, 2.12250,  // Tb-Er
    2.11040, 2.09930, 2.00650, 2.12250,  // Tm-Hf
    2.04900, 1.99275, 1.94775, 1.87450,  // Ta-Ir
    1.72280, 1.67625, 1.62820, 1.67995,  // Pt-Hg
    2.15635, 2.13820, 2.05875, 2.00270,  // Tl-Po
    1.93220, 1.86080, 2.53980, 2.46470,  // At-Ra
    2.35215, 2.21260, 2.22970, 2.19785,  // Ac-U
    2.17695, 2.21705                      // Np-Pu
};

GBSASolvation::GBSASolvation(int atomcount, const std::vector<int>& atoms,
    const Matrix& geometry, const std::string& solvent)
    : m_atomcount(atomcount)
    , m_atoms(atoms)
    , m_geometry(geometry)
    , m_solvent(solvent)
    , m_use_hbond(false)
    , m_initialized(false)
{
    // Allocate arrays
    m_born_radii = Vector::Zero(m_atomcount);
    m_born_matrix = Matrix::Zero(m_atomcount, m_atomcount);
    m_sasa = Vector::Zero(m_atomcount);
    m_gradient = Matrix::Zero(m_atomcount, 3);
    m_vdw_radii = Vector::Zero(m_atomcount);
    m_vdw_radii_offset = Vector::Zero(m_atomcount);
    m_descreening_radii = Vector::Zero(m_atomcount);

    // Load solvent parameters
    if (!loadSolventParameters(solvent)) {
        CurcumaLogger::warn("Failed to load solvent parameters for: " + solvent + ", using water");
        loadSolventParameters("water");
    }

    // Initialize van der Waals radii
    initializeVdWRadii();

    m_initialized = true;
}

bool GBSASolvation::loadSolventParameters(const std::string& solvent)
{
    m_params = getSolventParameters(solvent);
    return m_params.dielectric_constant > 0.0;
}

void GBSASolvation::initializeVdWRadii()
{
    // Initialize van der Waals radii for all atoms
    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];
        m_vdw_radii[i] = getVdWRadius(z);

        // svdw = vdw + born_offset
        m_vdw_radii_offset[i] = m_vdw_radii[i] + m_params.born_offset;

        // rho = descreening parameter * vdw
        double descreening = (z >= 1 && z <= 94) ? m_params.descreening[z - 1] : 1.0;
        m_descreening_radii[i] = descreening * m_vdw_radii[i];
    }
}

double GBSASolvation::getVdWRadius(int atomic_number) const
{
    if (atomic_number < 1 || atomic_number > 94) {
        return 2.0; // Default radius for unknown elements
    }
    return VDW_RADII_D3[atomic_number - 1];
}

void GBSASolvation::updateGeometry(const Matrix& geometry)
{
    m_geometry = geometry;
}

double GBSASolvation::calculateEnergy(const Vector& charges, bool gradient)
{
    if (!m_initialized) {
        CurcumaLogger::error("GBSA not initialized");
        return 0.0;
    }

    if (charges.size() != m_atomcount) {
        CurcumaLogger::error("Charge vector size mismatch in GBSA");
        return 0.0;
    }

    // Reset energy components
    m_gborn = 0.0;
    m_gsasa = 0.0;
    m_ghb = 0.0;

    // Step 1: Calculate Born radii
    calculateBornRadii();

    // Step 2: Construct Born interaction matrix
    constructBornMatrix();

    // Step 3: Calculate SASA
    calculateSASA();

    // Step 4: Calculate Born energy: E_Born = 0.5 * q^T * B * q
    Vector tmp = 0.5 * (m_born_matrix * charges);
    m_gborn = charges.dot(tmp);

    // Step 5: Calculate SASA energy: E_SASA = sum_i gamma_i * SASA_i
    m_gsasa = 0.0;
    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];
        double gamma = (z >= 1 && z <= 94) ? m_params.surface_tension[z - 1] : 0.0;
        m_gsasa += gamma * m_sasa[i];
    }

    // Step 6: Hydrogen bonding correction (optional)
    if (m_use_hbond) {
        for (int i = 0; i < m_atomcount; ++i) {
            int z = m_atoms[i];
            double hb_strength = (z >= 1 && z <= 94) ? m_params.hbond_strength[z - 1] : 0.0;
            m_ghb += hb_strength * charges[i] * charges[i];
        }
        m_gborn -= m_ghb; // Subtract H-bond from Born energy
    }

    // Step 7: Total solvation energy
    double e_total = m_gborn + m_gsasa + m_params.free_energy_shift;

    // Step 8: Calculate gradients if requested
    if (gradient) {
        // TODO: Implement gradient calculation
        // Reference: gfnff/src/gbsa/gbsa.f90 - addGradient()
        m_gradient.setZero();
    }

    return e_total;
}

void GBSASolvation::calculateBornRadii()
{
    /**
     * Born radii calculation using GBOBC integrator
     * Reference: gfnff/src/gbsa/born.f90 - compute_bornr()
     *
     * Formula:
     *   1. Calculate psi integral (neighbor overlap)
     *   2. Apply tanh correction:
     *      R_Born,i = c1 / (1/svdw_i - 1/vdw_i * tanh(arg))
     *      arg = psi * (alpha + psi * (gamma * psi - beta))
     */

    // Step 1: Calculate psi integral for each atom
    Vector psi = Vector::Zero(m_atomcount);

    for (int i = 0; i < m_atomcount; ++i) {
        double psi_i = 0.0;
        double rhoi = m_descreening_radii[i];
        double vdwi = m_vdw_radii[i];

        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j)
                continue;

            Vector ri = m_geometry.row(i);
            Vector rj = m_geometry.row(j);
            double r_ij = (ri - rj).norm();

            double rhoj = m_descreening_radii[j];
            double vdwj = m_vdw_radii[j];

            // GBOBC pairwise descreening integral
            // Reference: gfnff/src/gbsa/born.f90 - compute_psi()

            // Check overlap conditions
            double rh = rhoi + rhoj;
            double r1 = 1.0 / r_ij;

            // Case 1: Complete overlap (j inside i's sphere)
            if (r_ij < rhoi - rhoj) {
                // gi = 1/(rhoi - rhoj)
                double gi = 1.0 / (rhoi - rhoj);
                psi_i += gi;
                continue;
            }

            // Case 2: Complete overlap (i inside j's sphere)
            if (r_ij < rhoj - rhoi) {
                // gj = 1/(rhoj - rhoi)
                double gj = 1.0 / (rhoj - rhoi);
                psi_i += gj;
                continue;
            }

            // Case 3: Partial overlap or separated
            if (r_ij > rh) {
                // No overlap - negligible contribution
                continue;
            }

            // Calculate overlap integral
            double ap = r_ij + rhoj;
            double am = r_ij - rhoj;
            double lnab = std::log(ap / am);

            double gi = (1.0 / (r_ij - rhoj) - 1.0 / (r_ij + rhoj) + 0.25 * r1 * lnab) / rhoi;
            psi_i += gi;
        }

        psi[i] = psi_i;
    }

    // Step 2: Convert psi to Born radii with tanh correction
    for (int i = 0; i < m_atomcount; ++i) {
        double svdwi = m_vdw_radii_offset[i];
        double vdwi = m_vdw_radii[i];
        double s1 = 1.0 / svdwi;
        double v1 = 1.0 / vdwi;
        double s2 = 0.5 * svdwi;

        // Scale psi by svdw/2
        double br = psi[i] * s2;

        // GBOBCII tanh correction
        // arg = br * (alpha + br * (gamma * br - beta))
        double arg2 = br * (GAMMA_GBOBC * br - BETA_GBOBC);
        double arg = br * (ALPHA_GBOBC + arg2);

        double th = std::tanh(arg);

        // Born radius: R = c1 / (1/svdw - 1/vdw * tanh(arg))
        double br_final = m_params.born_scale / (s1 - v1 * th);

        m_born_radii[i] = br_final;
    }
}

void GBSASolvation::constructBornMatrix()
{
    /**
     * Construct Born interaction matrix
     * Reference: gfnff/src/gbsa/kernel.f90 - addBornMatStill()
     *
     * Still kernel:
     *   f_ij = 1/sqrt(r_ij² + R_i*R_j*exp(-r_ij²/(4*R_i*R_j)))
     * with dielectric screening: (1 - 1/epsilon) * f_ij
     */

    double keps = 1.0 - 1.0 / m_params.dielectric_constant;

    m_born_matrix.setZero();

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = 0; j <= i; ++j) {
            double bij = 0.0;

            if (i == j) {
                // Self-interaction: B_ii = 1/R_i
                bij = 1.0 / m_born_radii[i];
            } else {
                // Pairwise interaction using Still kernel
                Vector ri = m_geometry.row(i);
                Vector rj = m_geometry.row(j);
                double r_ij = (ri - rj).norm();

                double Ri = m_born_radii[i];
                double Rj = m_born_radii[j];

                // Still kernel: f_ij = 1/sqrt(r_ij² + Ri*Rj*exp(-r_ij²/(4*Ri*Rj)))
                double aiaj = Ri * Rj;
                double arg = -r_ij * r_ij / (4.0 * aiaj);
                double expterm = std::exp(arg);
                double fij = 1.0 / std::sqrt(r_ij * r_ij + aiaj * expterm);

                bij = fij;
            }

            // Apply dielectric screening
            bij *= keps;

            // Fill symmetric matrix
            m_born_matrix(i, j) = bij;
            m_born_matrix(j, i) = bij;
        }
    }
}

void GBSASolvation::calculateSASA()
{
    /**
     * Calculate Solvent Accessible Surface Area
     * Reference: gfnff/src/gbsa/sasa.f90 - compute_numsa()
     *
     * Simplified implementation using probe sphere rolling
     * Full implementation would use Lebedev angular grids for numerical integration
     */

    m_sasa.setZero();

    double probe_radius = m_params.probe_radius;

    for (int i = 0; i < m_atomcount; ++i) {
        double ri_sasa = m_vdw_radii[i] + probe_radius;

        // Simple approximation: SASA = 4*pi*r_sasa²
        // TODO: Implement proper numerical integration with neighbor exclusion
        // Reference: gfnff/src/gbsa/lebedev.f90 for angular grids
        double sasa_i = 4.0 * M_PI * ri_sasa * ri_sasa;

        // Neighbor correction (approximate)
        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j)
                continue;

            Vector rvec_i = m_geometry.row(i);
            Vector rvec_j = m_geometry.row(j);
            double r_ij = (rvec_i - rvec_j).norm();

            double rj_sasa = m_vdw_radii[j] + probe_radius;
            double cutoff = ri_sasa + rj_sasa;

            if (r_ij < cutoff) {
                // Simple overlap reduction
                double overlap_factor = 1.0 - (cutoff - r_ij) / cutoff;
                sasa_i *= overlap_factor;
            }
        }

        m_sasa[i] = sasa_i;
    }
}

void GBSASolvation::getEnergyComponents(double& e_born, double& e_sasa,
    double& e_hb, double& e_shift) const
{
    e_born = m_gborn;
    e_sasa = m_gsasa;
    e_hb = m_ghb;
    e_shift = m_params.free_energy_shift;
}

// ============================================================================
// Solvent Parameter Database
// ============================================================================

GBSAParameters getWaterParameters()
{
    /**
     * GFN2-xTB GBSA parameters for water
     * Reference: gfnff/include/param_gbsa_h2o.fh - gfn2_h2o
     */
    GBSAParameters params;
    params.solvent_name = "water";
    params.dielectric_constant = 80.2;
    params.molar_mass = 18.0;
    params.density = 1.0;
    params.born_scale = 1.55243817;
    params.probe_radius = 0.97531379;
    params.free_energy_shift = 1.16556316;
    params.born_offset = 0.13032636;
    params.alpha = 0.0; // GBSA (not ALPB)

    // Surface tension scaling (gamscale) for elements 1-94
    params.surface_tension = {
        -0.33498306, 1.26871861, -4.02682203, -6.61323872, -0.68114203,
        -0.74769065, -2.31291292, 0.91797911, 1.25543348, 0.85760021,
        -6.74101323, -1.19096758, -0.09427821, 1.53981552, -0.28566074,
        -0.48555420, -0.10208045, 0.27015708, -1.34530799, -1.82638776,
        -2.33181969, 1.94895791, -4.03488678, -3.69916894, -1.76072219,
        6.56392322, 2.00051154, -0.51490455, 0.01915301, -0.47581759,
        -0.68235927, 0.47915022, -0.53992351, 0.68840219, -0.28044414,
        0.17421825, -1.20561618, -0.83130281, -4.57650868, 0.56423064,
        -4.64819627, -4.63737691, -0.84355651, -2.92994729, -1.08086259,
        -0.58039850, -0.33243759, -0.03930933, -0.40836368, 0.19337851,
        0.09829368, -0.32995482, -0.83107400, -0.29130139, -1.25860685,
        -0.45587905, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.58524764, -0.54932016, 4.54316358, -2.41569825,
        0.37534414, -0.20293106, 0.79627337, 0.49884735, 0.27154747,
        0.59002423, 0.60254474, 1.17063879, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    // Descreening parameters (sx) for elements 1-94
    params.descreening = {
        0.71893869, 2.06028647, 0.40608033, 0.92184652, 0.77279023,
        0.74298311, 0.90261230, 0.75369019, 0.61625017, 1.84783975,
        1.20080317, 0.65802384, 0.72103171, 0.83603641, 0.66440982,
        0.89558269, 0.62985587, 1.27858461, 0.65299061, 0.62472181,
        0.98563134, 0.98178333, 1.09523331, 0.91537500, 0.84181077,
        0.78712610, 0.82166926, 0.77862239, 0.75925377, 1.06720895,
        0.78280891, 0.76954265, 0.64491049, 0.79964594, 0.63761709,
        0.54781361, 0.3, 0.53156870, 1.26880121, 1.12551818,
        1.25661787, 0.58174619, 1.21090637, 1.13642998, 1.17199113,
        1.09945430, 0.76967747, 0.48847954, 0.73965346, 0.71276027,
        0.67600887, 0.71055175, 0.41233944, 0.92879335, 0.62896377,
        0.68883366, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 1.13198750, 0.80446560, 0.59316114, 0.88315904,
        0.98359339, 0.88309289, 1.08473432, 1.09367781, 0.64535600,
        0.76250916, 0.41983434, 1.02977992, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8
    };

    // Hydrogen bond strength (tmp) for elements 1-94
    params.hbond_strength = {
        6.70894947, 0.0, 0.0, 0.0, 0.0,
        1.26459036, 3.52206160, 2.30440543, 1.98829409, 0.0,
        0.0, 0.0, 0.0, 0.0, 2.68116653,
        0.38262428, 1.02948365, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 2.37162766, 2.02275702,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 2.61716906, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    return params;
}

GBSAParameters getAcetoneParameters()
{
    // TODO: Implement acetone parameters from param_gbsa_acetone.fh
    GBSAParameters params;
    params.solvent_name = "acetone";
    params.dielectric_constant = 20.7;
    params.molar_mass = 58.08;
    params.density = 0.7845;
    // ... other parameters
    return params;
}

GBSAParameters getDMSOParameters()
{
    // TODO: Implement DMSO parameters from param_gbsa_dmso.fh
    GBSAParameters params;
    params.solvent_name = "dmso";
    params.dielectric_constant = 46.7;
    params.molar_mass = 78.13;
    params.density = 1.1;
    // ... other parameters
    return params;
}

GBSAParameters getSolventParameters(const std::string& solvent)
{
    std::string solv_lower = solvent;
    std::transform(solv_lower.begin(), solv_lower.end(), solv_lower.begin(), ::tolower);

    if (solv_lower == "water" || solv_lower == "h2o") {
        return getWaterParameters();
    } else if (solv_lower == "acetone") {
        return getAcetoneParameters();
    } else if (solv_lower == "dmso") {
        return getDMSOParameters();
    }

    // Default to water if unknown solvent
    CurcumaLogger::warn("Unknown solvent: " + solvent + ", defaulting to water");
    return getWaterParameters();
}

} // namespace GBSA
