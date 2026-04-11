/*
 * <GFN-FF Parameter Structures for Curcuma>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (March 2026): Extracted GFN-FF specific parameter structs
 * from forcefieldthread.h for clean separation of concerns.
 *
 * This header defines all GFN-FF specific interaction parameter structures
 * and the aggregate GFNFFParameterSet used for native in-memory parameter
 * transfer between GFNFF (generator) and ForceField (calculator).
 *
 * Generic structs (Bond, Angle, Dihedral, Inversion, vdW, EQ) remain in
 * forcefieldthread.h as they are shared by UFF/QMDFF/GFN-FF.
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "json.hpp"
using json = nlohmann::json;

// Forward-declare generic structs from forcefieldthread.h
struct Bond;
struct Angle;
struct Dihedral;
struct Inversion;
struct vdW;

/**
 * @brief Force field method type for unified parameter set
 *
 * Claude Generated (March 2026): Distinguishes UFF, QMDFF, and GFN-FF within
 * the shared ForceFieldParameterSet / FFWorkspace architecture.
 */
enum class FFMethodType { UFF = 1, QMDFF = 2, GFN_FF = 3 };

/**
 * @brief Bond-HB mapping entry for dncoord_erf calculation
 *
 * Claude Generated (Feb 21, 2026): Stores A-H...B mapping for HB coordination number
 * Reference: Fortran gfnff_data_types.f90:88,118-120 (bond_hb_AH, bond_hb_B, bond_hb_Bn)
 */
struct BondHBEntry {
    int A = 0;                     ///< Donor atom index (bonded to H, must be N or O)
    int H = 0;                     ///< Hydrogen atom index
    std::vector<int> B_atoms;      ///< Acceptor atom indices (N or O atoms)
};

/**
 * @brief HB coordination number gradient entry for egbond_hb chain-rule gradient
 *
 * Claude Generated (Feb 22, 2026): Stores d(hb_cn_H)/dr for each (H, B) pair.
 */
struct HBGradEntry {
    int H_atom = 0;                ///< Hydrogen atom index
    int B_atom = 0;                ///< Acceptor atom index
    Eigen::Vector3d dCN_dH;        ///< d(hb_cn_H)/dr_H for this (H, B) pair
    Eigen::Vector3d dCN_dB;        ///< d(hb_cn_H)/dr_B for this (H, B) pair (= -dCN_dH)
};

/**
 * @brief D3/D4 Dispersion pairwise term (GFN-FF modified formula)
 *
 * Reference: gfnff_gdisp0.f90:365-377
 * GFN-FF uses MODIFIED BJ damping: E = -0.5 * C6 * (t6 + 2*r4r2ij*t8)
 *
 * P1c (Apr 2026): Removed legacy D3 fields (C8, s6, s8, a1, a2).
 * 88 bytes → 48 bytes per pair. D3 uses separate D3DispersionPair struct.
 */
struct GFNFFDispersion {
    int i = 0, j = 0;           ///< Atom pair indices
    double C6 = 0.0;            ///< C6 dispersion coefficient (CN-weighted)
    double r4r2ij = 0.0;        ///< Product: 3 * sqrtZr4r2_i * sqrtZr4r2_j
    double r0_squared = 0.0;    ///< Pre-computed R0^2 = (a1*sqrt(r4r2ij) + a2)^2
    double r_cut = 50.0;        ///< Cutoff radius (Bohr)
    double zetac6 = 1.0;        ///< Zeta charge scaling
};

/**
 * @brief D3 Dispersion pairwise term (standard BJ damping)
 *
 * D3 uses standard BJ damping: E = -s6*C6/(r^6+R0^6) - s8*C8/(r^8+R0^8)
 * with R0 = a1*sqrt(C8/C6) + a2
 *
 * Claude Generated (Apr 2026): P1c — Separated from GFNFFDispersion to enable
 * struct size optimization (88 → 48 bytes for GFN-FF/D4 hot path).
 */
struct D3DispersionPair {
    int i = 0, j = 0;
    double C6 = 0.0;
    double C8 = 0.0;
    double s6 = 1.0;
    double s8 = 1.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double r_cut = 50.0;
};

/**
 * @brief GFN-FF Repulsion pairwise term
 *
 * Reference: gfnff_engrad.F90:407-439
 * Formula: E_rep = repab * exp(-alpha*r^1.5) / r
 */
struct GFNFFRepulsion {
    int i = 0, j = 0;
    double alpha = 0.0;
    double repab = 0.0;
    double r_cut = 50.0;
};

/**
 * @brief EEQ-based Coulomb electrostatics pairwise term
 *
 * Reference: Phase 3 EEQ charges
 * Formula: E_coul = q_i * q_j * erf(gamma_ij * r_ij) / r_ij
 *
 * P3a (Apr 2026): Shrunk from 120→56 bytes.
 * Removed redundant per-atom fields now stored in per-atom vectors:
 *   q_i/q_j → m_eeq_charges (dynamic EEQ charges)
 *   chi_i/chi_j → computed as chi_base + cnf*sqrt(cn) at runtime
 *   gam_i/gam_j → m_coulomb_gam per-atom vector
 *   alp_i/alp_j → m_coulomb_alp per-atom vector
 */
struct GFNFFCoulomb {
    int i = 0, j = 0;
    double gamma_ij = 0.0;
    double chi_base_i = 0.0, chi_base_j = 0.0;
    double cnf_i = 0.0, cnf_j = 0.0;
    double r_cut = 50.0;
};

/**
 * @brief GFN-FF Hydrogen Bond (HB) three-body term
 *
 * Reference: gfnff_engrad.F90 - abhgfnff_eg1() and abhgfnff_eg2new()
 * Structure: A-H...B (donor A, hydrogen H, acceptor B)
 */
struct GFNFFHydrogenBond {
    int i = 0;              ///< Donor atom A
    int j = 0;              ///< Hydrogen atom H
    int k = 0;              ///< Acceptor atom B

    double basicity_A = 0.0, basicity_B = 0.0;
    double acidity_A = 0.0, acidity_B = 0.0;
    double q_H = 0.0, q_A = 0.0, q_B = 0.0;
    double r_cut = 50.0;

    int case_type = 1;
    std::vector<int> neighbors_A;
    std::vector<int> neighbors_B;
    int acceptor_parent_index = -1;
    std::vector<int> neighbors_C;

    mutable double strength = 1.0;
    double sigmoid_width = 0.5;
};

/**
 * @brief GFN-FF Halogen Bond (XB) three-body term
 *
 * Reference: gfnff_engrad.F90 - rbxgfnff_eg()
 * Structure: A-X...B (donor A, halogen X, acceptor B)
 */
struct GFNFFHalogenBond {
    int i = 0;              ///< Donor atom A
    int j = 0;              ///< Halogen atom X
    int k = 0;              ///< Acceptor atom B

    double basicity_B = 0.0;
    double acidity_X = 0.0;
    double q_X = 0.0, q_B = 0.0;
    double r_cut = 50.0;

    mutable double strength = 1.0;
    double sigmoid_width = 0.5;
};

/**
 * @brief GFN-FF Triple Bond Torsion (sTors_eg)
 *
 * Reference: gfnff_engrad.F90:3454
 * Potential: E = erefhalf * (1 - cos(2*phi))
 */
struct GFNFFSTorsion {
    int i = 0, j = 0, k = 0, l = 0;
    double erefhalf = 3.75e-4;
};

/**
 * @brief ATM (Axilrod-Teller-Muto) three-body dispersion term
 *
 * Reference: external/cpp-d4/src/damping/atm.cpp
 * Used by both D3 and D4 (same formula, different C6 sources)
 */
struct ATMTriple {
    int i = 0, j = 0, k = 0;
    double C6_ij = 0.0, C6_ik = 0.0, C6_jk = 0.0;
    double s9 = 1.0, a1 = 0.0, a2 = 0.0, alp = 14.0;
    std::string atm_method = "d3";
    double triple_scale = 1.0;
};

/**
 * @brief GFN-FF Bonded ATM (batm) three-body term
 *
 * Reference: gfnff_ini.f90:745-779, gfnff_engrad.F90:562-603
 * Only for 1,4-pairs (bpair==3)
 */
struct GFNFFBatmTriple {
    int i = 0, j = 0, k = 0;
    double zb3atm_i = 0.0, zb3atm_j = 0.0, zb3atm_k = 0.0;
};

/**
 * @brief Complete GFN-FF parameter set for native in-memory transfer
 *
 * Claude Generated (March 2026): Aggregate of all GFN-FF interaction parameters.
 * Replaces JSON as internal data transfer between GFNFF (generator) and ForceField (calculator).
 * JSON serialization preserved only for file-based parameter caching.
 *
 * Usage:
 *   // In GFNFF: generate native parameter set
 *   GFNFFParameterSet params = generateGFNFFParameterSet();
 *   // In ForceField: consume directly (no JSON round-trip)
 *   m_forcefield->setGFNFFParameters(params);
 *   // For disk cache: serialize when needed
 *   json j = params.toJSON();
 */
struct GFNFFParameterSet {
    // Bonded terms (use generic structs from forcefieldthread.h)
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::vector<Dihedral> extra_dihedrals;
    std::vector<Inversion> inversions;
    std::vector<GFNFFSTorsion> storsions;

    // Non-bonded pairwise terms
    std::vector<GFNFFDispersion> dispersions;
    std::vector<GFNFFRepulsion> bonded_repulsions;
    std::vector<GFNFFRepulsion> nonbonded_repulsions;
    std::vector<GFNFFCoulomb> coulombs;

    // Three-body terms
    std::vector<GFNFFHydrogenBond> hbonds;
    std::vector<GFNFFHalogenBond> xbonds;
    std::vector<ATMTriple> atm_triples;
    std::vector<GFNFFBatmTriple> batm_triples;

    // Bond-HB mapping data
    std::vector<BondHBEntry> bond_hb_data;

    // UFF/QMDFF specific non-bonded pairs (empty for GFN-FF)
    std::vector<vdW> vdws;

    // Charges
    Eigen::VectorXd eeq_charges;
    Eigen::VectorXd topology_charges;

    // P3a (Apr 2026): Per-atom Coulomb self-energy parameters (extracted from pairs)
    Eigen::VectorXd eeq_gam;    ///< Per-atom chemical hardness (charge-corrected)
    Eigen::VectorXd eeq_alp;    ///< Per-atom alpha (charge-corrected)

    // Metadata
    std::string dispersion_method = "d4";
    double e0 = 0.0;

    // Method type (distinguishes UFF, QMDFF, GFN-FF in FFWorkspace)
    FFMethodType method_type = FFMethodType::GFN_FF;

    // Term enable flags
    bool dispersion_enabled = true;
    bool hbond_enabled = true;
    bool repulsion_enabled = true;
    bool coulomb_enabled = true;

    /**
     * @brief Serialize to JSON for file-based parameter caching
     * @return JSON object with all parameters
     */
    json toJSON() const;

    /**
     * @brief Deserialize from JSON (for loading cached parameters from disk)
     * @param j JSON object with parameters
     * @return GFNFFParameterSet
     */
    static GFNFFParameterSet fromJSON(const json& j);
};

/// Alias for backward compatibility — GFNFFParameterSet IS ForceFieldParameterSet
using ForceFieldParameterSet = GFNFFParameterSet;
