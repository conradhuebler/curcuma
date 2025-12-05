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

    // Pre-computed charge factors (calculated once from EEQ charges)
    double q_H = 0.0;           ///< EEQ charge on H
    double q_A = 0.0;           ///< EEQ charge on A
    double q_B = 0.0;           ///< EEQ charge on B

    // Geometry thresholds
    double r_cut = 50.0;        ///< Distance cutoff (Bohr)

    // Case 2 specific (neighbor orientation)
    int case_type = 1;          ///< 1 = simple, 2 = with orientation
    std::vector<int> neighbors_B;  ///< Neighbor indices of B (for Case 2)
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
    void addGFNFFInversion(const Inversion& inversions);
    void addGFNFFvdW(const vdW& vdWs);

    // Phase 4: GFN-FF pairwise non-bonded addition methods (Claude Generated 2025)
    void addGFNFFDispersion(const GFNFFDispersion& dispersion);
    void addGFNFFRepulsion(const GFNFFRepulsion& repulsion);
    void addGFNFFCoulomb(const GFNFFCoulomb& coulomb);

    // Phase 1.2: GFN-FF hydrogen bond and halogen bond addition methods (Claude Generated 2025)
    void addGFNFFHydrogenBond(const GFNFFHydrogenBond& hbond);
    void addGFNFFHalogenBond(const GFNFFHalogenBond& xbond);

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

    // Phase 1.2: HB/XB energy components (Claude Generated 2025)
    double HydrogenBondEnergy() { return m_energy_hbond; }
    double HalogenBondEnergy() { return m_energy_xbond; }

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
    void CalculateGFNFFInversionContribution();
    void CalculateGFNFFvdWContribution();

    // Phase 4: GFN-FF pairwise non-bonded calculation functions (Claude Generated 2025)
    void CalculateGFNFFDispersionContribution();
    void CalculateGFNFFRepulsionContribution();
    void CalculateGFNFFCoulombContribution();

    // Phase 5: GFN-FF hydrogen bond and halogen bond calculation functions (Claude Generated 2025)
    void CalculateGFNFFHydrogenBondContribution();
    void CalculateGFNFFHalogenBondContribution();

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
    std::vector<Dihedral> m_gfnff_dihedrals;
    std::vector<Inversion> m_gfnff_inversions;
    std::vector<vdW> m_gfnff_vdWs;  // Legacy (will be replaced by pairwise terms below)

    // Phase 4: GFN-FF pairwise parallelizable non-bonded terms
    std::vector<GFNFFDispersion> m_gfnff_dispersions;  // D3/D4 dispersion
    std::vector<GFNFFRepulsion> m_gfnff_repulsions;    // GFN-FF repulsion
    std::vector<GFNFFCoulomb> m_gfnff_coulombs;        // EEQ Coulomb electrostatics

    // Phase 1.2: GFN-FF hydrogen bond and halogen bond terms (Claude Generated 2025)
    std::vector<GFNFFHydrogenBond> m_gfnff_hbonds;     // Hydrogen bonds (HB)
    std::vector<GFNFFHalogenBond> m_gfnff_xbonds;      // Halogen bonds (XB)

protected:
    Matrix m_geometry, m_gradient;
    double m_energy = 0, m_bond_energy = 0.0, m_angle_energy = 0.0, m_dihedral_energy = 0.0, m_inversion_energy = 0.0, m_vdw_energy = 0.0, m_rep_energy = 0.0, m_eq_energy = 0.0;

    // Phase 4: Separate energy components for GFN-FF non-bonded terms
    double m_dispersion_energy = 0.0;  // D3/D4 dispersion
    double m_coulomb_energy = 0.0;     // EEQ Coulomb electrostatics

    // Phase 1.2: HB/XB energy components (Claude Generated 2025)
    double m_energy_hbond = 0.0;       // Hydrogen bond energy
    double m_energy_xbond = 0.0;       // Halogen bond energy

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
};

#ifdef USE_D3
class D3Thread : public ForceFieldThread {

public:
    D3Thread(int thread, int threads);
    ~D3Thread();
    virtual int Type() const { return 2; }

    void setParamater(const json& parameter)
    {
        m_d3->UpdateParametersD3(parameter);
    }

    void Initialise(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
        m_d3->InitialiseMolecule(m_atom_types);
    }
    virtual int execute() override;

private:
    DFTD3Interface* m_d3;
    std::vector<int> m_atom_types;
};
#endif

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
