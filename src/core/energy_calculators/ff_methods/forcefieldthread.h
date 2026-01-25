/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#pragma once

#include "src/core/global.h"

#include "src/core/hbonds.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#ifdef USE_D3
#include "src/core/energy_calculators/qm_methods/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/energy_calculators/qm_methods/dftd4interface.h"
#endif

#include "qmdff_par.h"
#include "uff_par.h"
#include "gfnff_geometry.h"  // Claude Generated (2025): GFN-FF geometry functions

#include <functional>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

struct Bond {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, distance = 0;
    double fc = 0, exponent = 0, r0_ij = 0, r0_ik = 0;
    double rabshift = 0.0;  // Claude Generated (Dec 2025): GFN-FF rabshift (vbond(1)) for validation
    double fqq = 1.0;       // Claude Generated (Jan 7, 2026): GFN-FF charge-dependent force constant factor

    // Claude Generated (Jan 18, 2026): Dynamic r0 calculation parameters
    // Reference: Fortran gfnff_rab.f90:147-153 - r0 recalculated at each Calculate()
    // Formula: r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff
    int z_i = 0, z_j = 0;           // Atomic numbers for parameter lookup
    double r0_base_i = 0.0;          // r0_gfnff[z_i-1] (Bohr)
    double r0_base_j = 0.0;          // r0_gfnff[z_j-1] (Bohr)
    double cnfak_i = 0.0;            // cnfak_gfnff[z_i-1]
    double cnfak_j = 0.0;            // cnfak_gfnff[z_j-1]
    double ff = 1.0;                 // EN-correction: 1 - k1*|ΔEN| - k2*ΔEN²

    // Claude Generated (Jan 24, 2026): Hydrogen bridge bond modulation (egbond_hb)
    // Reference: Fortran gfnff_engrad.F90:449-453, 919-994
    // For X-H bonds participating in HB: alpha_modified = (1 - 0.1*hb_cn_H) * alpha
    int nr_hb = 0;           // Number of HB interactions this bond participates in
    double hb_cn_H = 0.0;    // HB coordination number for hydrogen atom (used if nr_hb >= 1)
};

struct Angle {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0;
    double fc = 0, r0_ij = 0, r0_ik = 0, theta0_ijk = 0;
    double C0 = 0, C1 = 0, C2 = 0;
};

struct Dihedral {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, l = 0;
    double V = 0, n = 0, phi0 = 0;
    bool is_extra = false;  // Claude Generated (Jan 1, 2026): GFN-FF extra sp3-sp3 torsion flag
    bool is_nci = false;   // Claude Generated (Jan 13, 2026): NCI (Non-Covalent Interaction) torsion flag
};

struct Inversion {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, l = 0;
    double fc = 0, C0 = 0, C1 = 0, C2 = 0;
};

struct vdW {
    // === Existing members ===
    int type = 1; // 1 = UFF, 2 = QMDFF, 3 = CG
    int i = 0, j = 0;
    double C_ij = 0, r0_ij = 0;

    // === NEW: CG-specific parameters (only used when type=3) ===
    // Complete ellipsoid support structure (for future extensibility)
    Eigen::Vector3d shape_i = Eigen::Vector3d(2.0, 2.0, 2.0); // (x,y,z)-radii for atom i
    Eigen::Vector3d shape_j = Eigen::Vector3d(2.0, 2.0, 2.0); // (x,y,z)-radii for atom j
    Eigen::Vector3d orient_i = Eigen::Vector3d(0.0, 0.0, 0.0); // Euler angles for atom i
    Eigen::Vector3d orient_j = Eigen::Vector3d(0.0, 0.0, 0.0); // Euler angles for atom j

    // CG potential parameters
    double sigma = 4.0, epsilon = 0.0; // LJ parameters
    int cg_potential_type = 1; // 1=LJ_1612, 2=LJ_612, 3=tabulated
};

struct EQ {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0;
    double q_i = 0, q_j = 0, epsilon = 1;
};

/**
 * @brief GFN-FF specific non-bonded interaction structures
 *
 * Claude Generated (2025): Phase 4 pairwise parallelizable non-bonded terms
 * These are computed for all atom pairs and parallelized like UFF vdW terms.
 *
 * Design Philosophy:
 * - NOT add-on corrections (integrated into main energy loop)
 * - Pairwise parallelizable (like UFF vdW approach)
 * - Each pair (i,j) stored once with full parameters
 * - ForceFieldThread handles parallelization automatically
 */

/**
 * @brief D3/D4 Dispersion pairwise term
 *
 * Reference: gfnff_engrad.F90 (D4 integration)
 * Formula: E_disp = -Σ_ij f_damp(r_ij) * (C6_ij/r_ij^6 + C8_ij/r_ij^8 + ...)
 */
struct GFNFFDispersion {
    int i = 0, j = 0;           ///< Atom pair indices
    double C6 = 0.0;            ///< C6 dispersion coefficient
    double C8 = 0.0;            ///< C8 dispersion coefficient (optional)
    double r_cut = 50.0;        ///< Cutoff radius (Bohr)
    double s6 = 1.0;            ///< Scaling factor for C6
    double s8 = 1.0;            ///< Scaling factor for C8
    double a1 = 0.0;            ///< BJ damping parameter 1
    double a2 = 0.0;            ///< BJ damping parameter 2
};

/**
 * @brief GFN-FF Repulsion pairwise term
 *
 * Reference: gfnff_engrad.F90:407-439 (bonded repulsion)
 * Formula: E_rep = repab * exp(-α*r^1.5) / r
 */
struct GFNFFRepulsion {
    int i = 0, j = 0;           ///< Atom pair indices
    double alpha = 0.0;         ///< Repulsion exponent (sqrt(repa_i * repa_j))
    double repab = 0.0;         ///< Repulsion prefactor (repz_i * repz_j * scale)
    double r_cut = 50.0;        ///< Cutoff radius (Bohr)
};

/**
 * @brief EEQ-based Coulomb electrostatics pairwise term
 *
 * Reference: Phase 3 EEQ charges
 * Formula: E_coul = q_i * q_j * erf(γ_ij * r_ij) / r_ij
 */
struct GFNFFCoulomb {
    int i = 0, j = 0;           ///< Atom pair indices
    double q_i = 0.0;           ///< EEQ charge on atom i
    double q_j = 0.0;           ///< EEQ charge on atom j
    double gamma_ij = 0.0;      ///< Damping parameter (1/sqrt(α_i + α_j))
    double chi_i = 0.0;         ///< Electronegativity of atom i (for self-energy term)
    double chi_j = 0.0;         ///< Electronegativity of atom j (for self-energy term)
    double gam_i = 0.0;         ///< Chemical hardness γ_i (for self-interaction term)
    double gam_j = 0.0;         ///< Chemical hardness γ_j (for self-interaction term)
    double alp_i = 0.0;         ///< Chemical softness α_i (for self-interaction term)
    double alp_j = 0.0;         ///< Chemical softness α_j (for self-interaction term)
    double r_cut = 50.0;        ///< Cutoff radius (Bohr)
};

/**
 * @brief GFN-FF Hydrogen Bond (HB) three-body term
 *
 * Reference: gfnff_engrad.F90 - abhgfnff_eg1() and abhgfnff_eg2new()
 * Structure: A-H...B (donor A, hydrogen H, acceptor B)
 *
 * Case 1: Simple A...H...B geometry
 * Case 2: A-H...B with neighbor orientation damping
 *
 * Formula: E_HB = B_AH × C_AH × C_B × (-C_acidity × R_damp × Q_H^outl)
 * Damping: Combined short-range, long-range, and out-of-line damping
 *
 * Claude Generated (2025): Phase 1 - HB/XB Implementation
 */
struct GFNFFHydrogenBond {
    int i = 0;              ///< Donor atom A index (bonded to H)
    int j = 0;              ///< Hydrogen atom H index
    int k = 0;              ///< Acceptor atom B index

    // Parameter storage
    double basicity_A = 0.0;    ///< HB basicity of atom A (xhbas[A])
    double basicity_B = 0.0;    ///< HB basicity of atom B (xhbas[B])
    double acidity_A = 0.0;     ///< HB acidity of atom A (xhaci[A])
    double acidity_B = 0.0;     ///< HB acidity of atom B (xhaci[B])

    // Pre-computed charge factors (calculated once from EEQ charges)
    double q_H = 0.0;           ///< EEQ charge on H
    double q_A = 0.0;           ///< EEQ charge on A
    double q_B = 0.0;           ///< EEQ charge on B

    // Geometry thresholds
    double r_cut = 50.0;        ///< Distance cutoff (Bohr)

    // Case 2/3 specific (neighbor orientation)
    int case_type = 1;          ///< 1 = simple, 2 = case 2 (oriented), 3 = case 3 (carbonyl/nitro)
    std::vector<int> neighbors_A;  ///< Neighbor indices of donor A (for Case 2/3)
    std::vector<int> neighbors_B;  ///< Neighbor indices of acceptor B (for Case 2/3)
    int acceptor_parent_index = -1; ///< Parent of acceptor (e.g., C in C=O) for Case 3
};

/**
 * @brief GFN-FF Halogen Bond (XB) three-body term
 *
 * Reference: gfnff_engrad.F90 - rbxgfnff_eg()
 * Structure: A-X...B (donor A, halogen X, acceptor B)
 *
 * Formula: E_XB = -R_damp × Q_outl × C_B × Q_B × C_X × Q_X
 * Damping: Uses XB-specific parameters (XB_BACUT, XB_SCUT, XB_LONGCUT)
 *
 * Claude Generated (2025): Phase 1 - HB/XB Implementation
 */
struct GFNFFHalogenBond {
    int i = 0;              ///< Donor atom A index (bonded to X)
    int j = 0;              ///< Halogen atom X index
    int k = 0;              ///< Acceptor atom B index

    // Parameter storage
    double basicity_B = 0.0;    ///< XB basicity of B (xhbas[B])
    double acidity_X = 0.0;     ///< XB acidity of X (xbaci[X])

    // Pre-computed charge factors
    double q_X = 0.0;           ///< EEQ charge on X
    double q_B = 0.0;           ///< EEQ charge on B

    // Geometry thresholds
    double r_cut = 50.0;        ///< Distance cutoff (Bohr)
};

/**
 * @brief ATM (Axilrod-Teller-Muto) three-body dispersion term
 *
 * Reference: external/cpp-d4/src/damping/atm.cpp
 * Formula: E_ATM = sum_{i<j<k} ang * fdmp * C9 / 3 * triple_scale
 *
 * Used by both D3 and D4 (same formula, different C6 sources)
 * - D3: C6 from CN-weighted interpolation (D3ParameterGenerator)
 * - D4: C6 from CN+charge-weighted Casimir-Polder (D4ParameterGenerator)
 *
 * Claude Generated (2025): ATM three-body dispersion
 */
struct ATMTriple {
    int i = 0;              ///< First atom index
    int j = 0;              ///< Second atom index (i < j)
    int k = 0;              ///< Third atom index (j < k)

    // Pre-computed C6 coefficients (from D3/D4 generator)
    double C6_ij = 0.0;     ///< Pairwise C6 coefficient (i,j)
    double C6_ik = 0.0;     ///< Pairwise C6 coefficient (i,k)
    double C6_jk = 0.0;     ///< Pairwise C6 coefficient (j,k)

    // Damping parameters
    double s9 = 1.0;        ///< Global scaling factor for C9
    double a1 = 0.0;        ///< BJ damping parameter a1
    double a2 = 0.0;        ///< BJ damping parameter a2 (Bohr)
    double alp = 14.0;      ///< Damping exponent

    // Method identifier
    std::string atm_method = "d3";  ///< "d3" or "d4"

    // Symmetry factor (pre-computed during generation)
    double triple_scale = 1.0;  ///< 1.0 (ijk), 0.5 (iij/ijj), 1/6 (iii)
};

/**
 * @brief GFN-FF Bonded ATM (batm) three-body term
 *
 * Reference: external/gfnff/src/gfnff_ini.f90:745-779, gfnff_engrad.F90:562-603, batmgfnff_eg:3267-3334
 * Formula: E_batm = sum_{1,4-pairs} c9 * (ang + 1.0) / rav3
 * where:
 *   c9 = ff * zb3atm_i * zb3atm_j * zb3atm_k
 *   ff = (1 - fqq*q_i) * (1 - fqq*q_j) * (1 - fqq*q_k) with fqq=3.0
 *   ang = 0.375 * (r_ik*r_jk - r_ij^2) * (r_ij*r_ik - r_jk^2) * (r_ij*r_jk - r_ik^2) / r_ijk^3
 *   rav3 = (r_ij*r_jk*r_ik)^1.5
 *
 * CLAUDE GENERATED (January 17, 2026): GFN-FF bonded ATM implementation
 * Key difference from general ATM: Only for 1,4-pairs (bpair==3)
 * - batm: O(N_bonds) terms - restricted topology
 * - D3/D4 ATM: O(N³) terms - general all-atom combinations
 */
struct GFNFFBatmTriple {
    int i = 0;              ///< 1,4-pair first atom index
    int j = 0;              ///< 1,4-pair second atom index (bpair[i][j] == 3)
    int k = 0;              ///< Center atom neighbor (either neighbor of i or j)

    // Pre-computed element parameters
    double zb3atm_i = 0.0;  ///< zb3atm parameter for atom i
    double zb3atm_j = 0.0;  ///< zb3atm parameter for atom j
    double zb3atm_k = 0.0;  ///< zb3atm parameter for atom k
};

class ForceFieldThread : public CxxThread {

public:
    ForceFieldThread(int thread, int threads);
    virtual int execute() override;
    virtual int Type() const { return m_method; }  // Claude Generated (2025): Return actual method type, not hardcoded 1!
    void addBond(const Bond& bonds);
    void addAngle(const Angle& angles);
    void addDihedral(const Dihedral& dihedrals);
    void addInversion(const Inversion& inversions);
    void addvdW(const vdW& vdWs);
    void addEQ(const EQ& EQs);

    void addGFNFFBond(const Bond& bonds);
    void addGFNFFAngle(const Angle& angles);
    void addGFNFFDihedral(const Dihedral& dihedrals);
    void addGFNFFExtraTorsion(const Dihedral& extra_torsion);  // Claude Generated (Jan 2, 2026): Extra sp3-sp3 gauche torsions
    void addGFNFFInversion(const Inversion& inversions);
    void addGFNFFvdW(const vdW& vdWs);

    // Phase 4: GFN-FF pairwise non-bonded addition methods (Claude Generated 2025)
    void addGFNFFDispersion(const GFNFFDispersion& dispersion);
    void addGFNFFBondedRepulsion(const GFNFFRepulsion& repulsion);
    void addGFNFFNonbondedRepulsion(const GFNFFRepulsion& repulsion);
    void addGFNFFCoulomb(const GFNFFCoulomb& coulomb);

    // D3/D4 parameter integration methods
    void addD3Dispersion(const GFNFFDispersion& d3_dispersion);
    void addD4Dispersion(const GFNFFDispersion& d4_dispersion);
    void CalculateD3DispersionContribution();
    void CalculateD4DispersionContribution();  // Claude Generated - Dec 25, 2025

    // Phase 1.2: GFN-FF hydrogen bond and halogen bond addition methods (Claude Generated 2025)
    void addGFNFFHydrogenBond(const GFNFFHydrogenBond& hbond);
    void addGFNFFHalogenBond(const GFNFFHalogenBond& xbond);

    // ATM three-body dispersion (D3/D4)
    void addATMTriple(const ATMTriple& triple);

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    // GFN-FF bonded ATM (batm) three-body terms for 1,4-pairs
    void addGFNFFBatmTriple(const GFNFFBatmTriple& batm_triple);
    void CalculateGFNFFBatmContribution();

    // Phase 6: Assign atoms for self-energy calculation (Claude Generated Dec 2025)
    // Thread-safe: Each thread calculates self-energy only for its assigned atoms
    void assignAtomsForSelfEnergy(const std::vector<int>& atom_indices);

    // Phase 3: Initialize atom types for covalent radius calculations in GFN-FF
    void Initialise(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
    }

    // Phase 5A: Set EEQ charges for fqq calculation (Claude Generated Nov 2025)
    void setEEQCharges(const Vector& charges)
    {
        m_eeq_charges = charges;
    }

    // Claude Generated (Jan 18, 2026): Set D3 coordination numbers for dynamic r0 calculation
    // Reference: Fortran gfnff_engrad.F90:432 - CN recalculated at each energy evaluation
    void setD3CN(const Vector& d3_cn)
    {
        m_d3_cn = d3_cn;
    }

    inline void UpdateGeometry(const Matrix& geometry, bool gradient)
    {
        m_geometry = geometry;
        m_calculate_gradient = gradient;
        m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
    }

    inline void setGeometry(const Matrix& geometry, bool gradient)
    {
        m_geometry = geometry;
        m_calculate_gradient = gradient;
        m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
    }

    inline void setMethod(int method)
    {
        m_method = method;
    }
    double BondEnergy() { return m_bond_energy; }
    double AngleEnergy() { return m_angle_energy; }
    double DihedralEnergy() { return m_dihedral_energy; }
    double InversionEnergy() { return m_inversion_energy; }
    double VdWEnergy() { return m_vdw_energy; }
    double RepEnergy() { return m_rep_energy; }
    double EQEnergy() { return m_eq_energy; }

    // Phase 4: GFN-FF pairwise non-bonded energy components (Claude Generated 2025)
    double DispersionEnergy() { return m_dispersion_energy; }
    double CoulombEnergy() { return m_coulomb_energy; }

    // Claude Generated 2025: Native D3/D4 energy components
    double D3Energy() { return m_d3_energy; }
    double D4Energy() { return m_d4_energy; }

    // Phase 1.2: HB/XB energy components (Claude Generated 2025)
    double HydrogenBondEnergy() { return m_energy_hbond; }
    double HalogenBondEnergy() { return m_energy_xbond; }

    // Claude Generated (December 2025): ATM three-body dispersion energy
    double ATMEnergy() { return m_atm_energy; }

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    double BatmEnergy() { return m_batm_energy; }

    // Phase 2: GFN-FF parameter flag setters (Claude Generated Dec 2025)
    void setDispersionEnabled(bool enabled) { m_dispersion_enabled = enabled; }
    void setHydrogenBondEnabled(bool enabled) { m_hbond_enabled = enabled; }
    void setRepulsionEnabled(bool enabled) { m_repulsion_enabled = enabled; }
    void setCoulombEnabled(bool enabled) { m_coulomb_enabled = enabled; }

    Matrix Gradient() const { return m_gradient; }

private:
    void CalculateUFFBondContribution();
    void CalculateUFFAngleContribution();
    void CalculateUFFDihedralContribution();
    void CalculateUFFInversionContribution();
    void CalculateUFFvdWContribution();

    void CalculateQMDFFBondContribution();
    void CalculateQMDFFAngleContribution();
    void CalculateQMDFFDihedralContribution();
    void CalculateQMDFFEspContribution();
    void CalculateESPContribution();

    void CalculateGFNFFBondContribution();
    void CalculateGFNFFAngleContribution();
    void CalculateGFNFFDihedralContribution();
    void CalculateGFNFFExtraTorsionContribution();  // Claude Generated (Jan 2, 2026): Extra sp3-sp3 gauche torsions
    void CalculateGFNFFInversionContribution();
    void CalculateGFNFFvdWContribution();

    // Phase 4: GFN-FF pairwise non-bonded calculation functions (Claude Generated 2025)
    void CalculateGFNFFDispersionContribution();
    void CalculateGFNFFBondedRepulsionContribution();
    void CalculateGFNFFNonbondedRepulsionContribution();
    void CalculateGFNFFCoulombContribution();

    // Phase 5: GFN-FF hydrogen bond and halogen bond calculation functions (Claude Generated 2025)
    void CalculateGFNFFHydrogenBondContribution();
    void CalculateGFNFFHalogenBondContribution();

    // ATM three-body dispersion calculation functions (Claude Generated 2025)
    void CalculateATMContribution();
    void CalculateATMGradient();  // Claude Generated (2025): Analytical ATM gradients

    // double HarmonicBondStretching();

    // double LJBondStretching();

    std::function<double()> CalculateBondContribution;
    std::function<double()> CalculateAngleContribution;
    std::function<double()> CalculateTorsionContribution;
    std::function<double()> CalculateInversionContribution;
    std::function<double()> CalculateVdWContribution;
    // std::function<double()> CalculateESPContribution;
    std::function<double()> CalculateHBondContribution;

    std::vector<Bond> m_uff_bonds;
    std::vector<Angle> m_uff_angles;
    std::vector<Dihedral> m_uff_dihedrals, m_qmdff_dihedrals;
    std::vector<Inversion> m_uff_inversions, m_qmdff_inversions;
    std::vector<vdW> m_uff_vdWs;
    std::vector<EQ> m_EQs;

    std::vector<Bond> m_gfnff_bonds;
    std::vector<Angle> m_gfnff_angles;
    std::vector<Dihedral> m_gfnff_dihedrals;        // Primary torsions (n=3, n=2, etc.)
    std::vector<Dihedral> m_gfnff_extra_torsions;   // Extra sp3-sp3 gauche torsions (n=1) - Claude Generated (Jan 2, 2026)
    std::vector<Inversion> m_gfnff_inversions;
    std::vector<vdW> m_gfnff_vdWs;  // Legacy (will be replaced by pairwise terms below)

    // Phase 4: GFN-FF pairwise parallelizable non-bonded terms
    std::vector<GFNFFDispersion> m_gfnff_dispersions;     // D3/D4 dispersion
    std::vector<GFNFFRepulsion> m_gfnff_bonded_repulsions;    // GFN-FF bonded repulsion
    std::vector<GFNFFRepulsion> m_gfnff_nonbonded_repulsions; // GFN-FF non-bonded repulsion
    std::vector<GFNFFCoulomb> m_gfnff_coulombs;         // EEQ Coulomb electrostatics

    // Phase 6: Atom assignment for self-energy calculation (Claude Generated Dec 2025)
    // Each thread gets a subset of atoms to calculate self-energy terms
    // This prevents duplicate self-energy when same atom appears in multiple threads
    std::vector<int> m_assigned_atoms_for_self_energy;

    // D3/D4 native dispersion pairs
    std::vector<GFNFFDispersion> m_d3_dispersions;  // Native D3 parameters
    std::vector<GFNFFDispersion> m_d4_dispersions;  // Native D4 parameters

    // Phase 1.2: GFN-FF hydrogen bond and halogen bond terms (Claude Generated 2025)
    std::vector<GFNFFHydrogenBond> m_gfnff_hbonds;     // Hydrogen bonds (HB)
    std::vector<GFNFFHalogenBond> m_gfnff_xbonds;      // Halogen bonds (XB)

    // ATM three-body dispersion (D3/D4)
    std::vector<ATMTriple> m_atm_triples;  // ATM three-body terms

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    // GFN-FF bonded ATM (batm) three-body terms for 1,4-pairs
    std::vector<GFNFFBatmTriple> m_gfnff_batms;  // Batm triples

protected:
    Matrix m_geometry, m_gradient;
    double m_energy = 0, m_bond_energy = 0.0, m_angle_energy = 0.0, m_dihedral_energy = 0.0, m_inversion_energy = 0.0, m_vdw_energy = 0.0, m_rep_energy = 0.0, m_eq_energy = 0.0;

    // Phase 4: Separate energy components for GFN-FF non-bonded terms
    double m_dispersion_energy = 0.0;  // D3/D4 dispersion
    double m_coulomb_energy = 0.0;     // EEQ Coulomb electrostatics
    double m_d3_energy = 0.0;          // Native D3 dispersion
    double m_d4_energy = 0.0;          // Native D4 dispersion

    // Phase 1.2: HB/XB energy components (Claude Generated 2025)
    double m_energy_hbond = 0.0;       // Hydrogen bond energy
    double m_energy_xbond = 0.0;       // Halogen bond energy
    double m_atm_energy = 0.0;         // ATM three-body dispersion energy (Claude Generated December 2025)

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    double m_batm_energy = 0.0;        // Batm three-body dispersion energy

    double m_final_factor = 1;
    double m_bond_scaling = 1, m_angle_scaling = 1, m_dihedral_scaling = 1, m_inversion_scaling = 1, m_vdw_scaling = 1, m_rep_scaling = 1;
    double m_au = 1;
    double m_d = 1e-3;
    int m_calc_gradient = 1;
    int m_thread = 0, m_threads = 0, m_method = 1;
    bool m_calculate_gradient = true;

    // Phase 3: Atom types for covalent radius calculations in GFN-FF
    std::vector<int> m_atom_types;

    // Phase 5A: EEQ charges for fqq angle correction (Claude Generated Nov 2025)
    Vector m_eeq_charges;

    // Claude Generated (Jan 18, 2026): D3 coordination numbers for dynamic r0 calculation
    // These are recalculated from current geometry at each Calculate() call
    Vector m_d3_cn;

    // Phase 1.2: Cached bonded pairs for fast lookup in repulsion calculation (Claude Generated - Dec 2025)
    // Built once in execute() to avoid O(N_bonds × log(N_bonds)) overhead per energy call
    std::set<std::pair<int, int>> m_bonded_pairs;
    bool m_bonded_pairs_cached = false;

    // Phase 2: Parameter flags for GFN-FF term control (Claude Generated Dec 2025)
    // These control which energy terms are calculated - saves CPU time for disabled terms
    bool m_dispersion_enabled = true;
    bool m_hbond_enabled = true;
    bool m_repulsion_enabled = true;
    bool m_coulomb_enabled = true;
};

class H4Thread : public ForceFieldThread {

public:
    H4Thread(int thread, int threads);
    ~H4Thread();
    virtual int Type() const { return 3; }

    void setParamater(const json& parameter)
    {
        m_h4correction.set_OH_O(parameter["h4_oh_o"].get<double>());
        m_h4correction.set_OH_N(parameter["h4_oh_n"].get<double>());
        m_h4correction.set_NH_O(parameter["h4_nh_o"].get<double>());
        m_h4correction.set_NH_N(parameter["h4_nh_n"].get<double>());

        m_h4correction.set_WH_O(parameter["h4_wh_o"].get<double>());
        m_h4correction.set_NH4(parameter["h4_nh4"].get<double>());
        m_h4correction.set_COO(parameter["h4_coo"].get<double>());
        m_h4correction.set_HH_Rep_K(parameter["hh_rep_k"].get<double>());
        m_h4correction.set_HH_Rep_E(parameter["hh_rep_e"].get<double>());
        m_h4correction.set_HH_Rep_R0(parameter["hh_rep_r0"].get<double>());
    }

    void Initialise(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
        m_h4correction.allocate(m_atom_types.size());
    }
    virtual int execute() override;

private:
    hbonds4::H4Correction m_h4correction;
    std::vector<int> m_atom_types;
};
