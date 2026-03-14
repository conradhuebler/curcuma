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

#include "gfnff_parameters.h"  // Claude Generated (March 2026): GFN-FF specific parameter structs
#include "qmdff_par.h"
#include "uff_par.h"
#include "gfnff_geometry.h"  // Claude Generated (2025): GFN-FF geometry functions

#include <functional>
#include <set>
#include <vector>
#include <unordered_map>  // Claude Generated (February 2026): For term timing storage
#include <chrono>         // Claude Generated (February 2026): For timing measurements

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

// BondHBEntry, HBGradEntry, GFNFFDispersion, GFNFFRepulsion, GFNFFCoulomb,
// GFNFFHydrogenBond, GFNFFHalogenBond, GFNFFSTorsion, ATMTriple, GFNFFBatmTriple
// are defined in gfnff_parameters.h (included above)

struct Bond {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0;
    double distance = 0.0;
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
    int type = 1; // 1 = UFF, 2 = QMDFF, 3 = GFN-FF
    int i = 0, j = 0, k = 0, l = 0;
    double fc = 0, C0 = 0, C1 = 0, C2 = 0;
    // GFN-FF specific: potential_type 0 = V*(1-cos(omega))*damp, -1 = V*(cos(omega)-cos(omega0))^2*damp
    int potential_type = 0;
    double omega0 = 0.0; // Equilibrium out-of-plane angle (radians), used for potential_type=-1
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

// GFN-FF specific structs (GFNFFDispersion, GFNFFRepulsion, GFNFFCoulomb,
// GFNFFHydrogenBond, GFNFFHalogenBond, GFNFFSTorsion, ATMTriple, GFNFFBatmTriple,
// BondHBEntry, HBGradEntry) are defined in gfnff_parameters.h

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
    void addGFNFFSTorsion(const GFNFFSTorsion& storsion);  // Claude Generated (March 2026): Triple bond torsions
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

    // Claude Generated (Feb 15, 2026): HB/XB bulk setters for MD dynamic updates
    // Reference: Fortran gfnff_engrad.F90:246-260 - rebuild HB/XB lists at each energy call
    void setGFNFFHBonds(const std::vector<GFNFFHydrogenBond>& hbonds);
    void setGFNFFHalogenBonds(const std::vector<GFNFFHalogenBond>& xbonds);

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

    // Claude Generated (Feb 21, 2026): Set Phase-1 topology charges for BATM
    // Reference: Fortran gfnff_engrad.F90:620 uses topo%qa (Phase-1, fixed) for BATM
    // CRITICAL: BATM must use Phase-1 charges, not Phase-2 EEQ charges
    void setTopologyCharges(const Vector& charges)
    {
        m_topology_charges = charges;
    }

    // Claude Generated (Jan 18, 2026): Set D3 coordination numbers for dynamic r0 calculation
    // Reference: Fortran gfnff_engrad.F90:432 - CN recalculated at each energy evaluation
    void setD3CN(const Vector& d3_cn)
    {
        m_d3_cn = d3_cn;
    }

    // Claude Generated (Feb 1, 2026): Set CN and CNF for Coulomb charge derivative gradients
    // Reference: Fortran gfnff_engrad.F90:418-422 - qtmp(i) = q(i)*cnf(i)/(2*sqrt(cn(i)))
    void setCN(const Vector& cn)
    {
        m_cn = cn;
    }

    void setCNF(const Vector& cnf)
    {
        m_cnf = cnf;
    }

    // Claude Generated (Feb 1, 2026): Set CN derivatives for Coulomb charge derivative gradients
    // dcn[dim](i,j) = dCN(j)/dr(i,dim) - how atom i's position affects atom j's CN
    // Reference: Fortran gfnff_engrad.F90:422 - call gemv(dcn,qtmp,g,alpha=-1.0,beta=1.0)
    void setCNDerivatives(const std::vector<Matrix>& dcn)
    {
        m_dcn = dcn;
    }

    // Claude Generated (Feb 21, 2026): Bond-HB mapping for dncoord_erf calculation
    // Reference: Fortran gfnff_engrad.F90:1069-1120
    // These entries define which H atoms participate in HBs and which B atoms to count
    void setBondHBData(const std::vector<BondHBEntry>& bond_hb_data) {
        m_bond_hb_data = bond_hb_data;
    }

    // Claude Generated (Feb 15, 2026): dEdcn accumulator for CN chain-rule gradient terms
    // Accumulates dE/dCN contributions from bond (dr0/dCN) and dispersion (dC6/dCN) terms
    // Applied after thread completion via dcn chain rule: gradient += dcn * dEdcn
    const Vector& getDEdcn() const { return m_dEdcn; }

    // Claude Generated (Mar 2026): Bond-specific dEdcn for per-component gradient attribution
    // Separates bond dr0/dCN from dispersion dC6/dCN so CN corrections go to the right component
    const Vector& getDEdcnBond() const { return m_dEdcn_bond; }

    // Claude Generated (Feb 15, 2026): Set dc6dcn matrix for dispersion CN gradient
    // Reference: Fortran gfnff_gdisp0.f90:262-305 - dc6dcn(i,j) = dC6(i,j)/dCN(i)
    void setDispersionDC6DCN(const Matrix& dc6dcn) { m_dc6dcn = dc6dcn; }

    inline void UpdateGeometry(const Matrix& geometry, bool gradient)
    {
        m_geometry = geometry;
        m_calculate_gradient = gradient;
        m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
        m_dEdcn = Vector::Zero(m_geometry.rows());
        m_dEdcn_bond = Vector::Zero(m_geometry.rows());
        if (m_store_gradient_components) {
            initGradientComponents(m_geometry.rows());
        }
    }

    inline void setGeometry(const Matrix& geometry, bool gradient)
    {
        m_geometry = geometry;
        m_calculate_gradient = gradient;
        m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
        m_dEdcn = Vector::Zero(m_geometry.rows());
        m_dEdcn_bond = Vector::Zero(m_geometry.rows());
        if (m_store_gradient_components) {
            initGradientComponents(m_geometry.rows());
        }
    }

    inline void setMethod(int method)
    {
        m_method = method;
        // Claude Generated (Feb 1, 2026): Update distance multiplier when method changes
        // Method 3 (GFN-FF) and 5 (D3-only) use Bohr coordinates internally
        if (m_method == 3 || m_method == 5) {
            m_au = 1.0;
        } else {
            m_au = 1.889726125; // Angstrom to Bohr for UFF-based non-bonded terms
        }
    }
    double BondEnergy() { return m_bond_energy; }
    double AngleEnergy() { return m_angle_energy; }
    double DihedralEnergy() { return m_dihedral_energy; }
    double InversionEnergy() { return m_inversion_energy; }
    double VdWEnergy() { return m_vdw_energy; }
    double RepEnergy() { return m_rep_energy; }
    double BondedRepEnergy() { return m_bonded_rep_energy; }
    double NonbondedRepEnergy() { return m_nonbonded_rep_energy; }
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

    // Claude Generated (February 2026): Accessor for individual term timings
    inline const std::unordered_map<std::string, long long>& getTermTimings() const {
        return m_term_timings;
    }

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    double BatmEnergy() { return m_batm_energy; }

    // Claude Generated (March 2026): Triple bond torsions
    double STorsEnergy() { return m_stors_energy; }

    // Phase 2: GFN-FF parameter flag setters (Claude Generated Dec 2025)
    void setDispersionEnabled(bool enabled) { m_dispersion_enabled = enabled; }
    void setHydrogenBondEnabled(bool enabled) { m_hbond_enabled = enabled; }
    void setRepulsionEnabled(bool enabled) { m_repulsion_enabled = enabled; }
    void setCoulombEnabled(bool enabled) { m_coulomb_enabled = enabled; }

    Matrix Gradient() const { return m_gradient; }

    // Claude Generated (February 2026): Per-component gradient storage for validation
    // Activated by setStoreGradientComponents(true) before calculation
    void setStoreGradientComponents(bool store) { m_store_gradient_components = store; }
    bool storeGradientComponents() const { return m_store_gradient_components; }

    // Per-component gradient getters (only valid if m_store_gradient_components == true)
    const Matrix& GradientBond() const { return m_gradient_bond; }
    const Matrix& GradientAngle() const { return m_gradient_angle; }
    const Matrix& GradientTorsion() const { return m_gradient_torsion; }
    const Matrix& GradientRepulsion() const { return m_gradient_repulsion; }
    const Matrix& GradientCoulomb() const { return m_gradient_coulomb; }
    const Matrix& GradientDispersion() const { return m_gradient_dispersion; }
    const Matrix& GradientHB() const { return m_gradient_hb; }
    const Matrix& GradientXB() const { return m_gradient_xb; }
    const Matrix& GradientBATM() const { return m_gradient_batm; }
    const Matrix& GradientATM() const { return m_gradient_atm; }   ///< ATM three-body dispersion gradient (Claude Generated Mar 2026)

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
    void computeHBCoordinationNumbers();   // Claude Generated (Feb 21, 2026): dncoord_erf for bond-HB CN
    void CalculateGFNFFAngleContribution();
    void CalculateGFNFFDihedralContribution();
    void CalculateGFNFFExtraTorsionContribution();  // Claude Generated (Jan 2, 2026): Extra sp3-sp3 gauche torsions
    void CalculateGFNFFInversionContribution();
    void CalculateGFNFFSTorsionContribution();  // Claude Generated (March 2026): Triple bond torsions
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
    std::vector<GFNFFSTorsion> m_gfnff_storsions;  // Triple bond torsions
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

    // Claude Generated (Feb 21, 2026): Bond-HB mapping for dncoord_erf calculation
    std::vector<BondHBEntry> m_bond_hb_data;

    // Claude Generated (Feb 22, 2026): HB gradient entries for chain-rule gradient
    // Computed in computeHBCoordinationNumbers(), applied in CalculateGFNFFBondContribution()
    // Reference: Fortran gfnff_engrad.F90:1054-1063 (hb_dcn)
    std::vector<HBGradEntry> m_hb_grad_entries;

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
    double m_stors_energy = 0.0;       // Triple bond torsion energy (sTors_eg)
    double m_atm_energy = 0.0;         // ATM three-body dispersion energy (Claude Generated December 2025)

    // Claude Generated (March 2026): Separate bonded/non-bonded repulsion for diagnostics
    double m_bonded_rep_energy = 0.0;    // GFN-FF bonded repulsion (REPSCALB=1.7583)
    double m_nonbonded_rep_energy = 0.0; // GFN-FF non-bonded repulsion (REPSCALN=0.4270)

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

    // Claude Generated (Feb 21, 2026): Phase-1 topology charges for BATM
    // Reference: Fortran gfnff_engrad.F90:620 passes topo%qa (Phase-1) to batmgfnff_eg
    // These are FIXED charges computed once at initialization, unlike m_eeq_charges which
    // are geometry-dependent. Using Phase-2 charges in BATM causes energy/gradient inconsistency
    // because the BATM gradient has no dq/dx term → energy drift in MD.
    Vector m_topology_charges;

    // Claude Generated (Jan 18, 2026): D3 coordination numbers for dynamic r0 calculation
    // These are recalculated from current geometry at each Calculate() call
    Vector m_d3_cn;

    // Claude Generated (Feb 1, 2026): CN, CNF, and CN derivatives for Coulomb charge derivative gradients
    // Reference: Fortran gfnff_engrad.F90:418-422 - charge derivative via CN
    Vector m_cn;                              // Coordination numbers per atom
    Vector m_cnf;                             // CNF parameters per atom (cnf_eeq from gfnff_par.h)
    std::vector<Matrix> m_dcn;                // CN derivatives: dcn[dim](i,j) = dCN(j)/dr(i,dim)

    // Claude Generated (Feb 15, 2026): dE/dCN accumulator for bond dr0/dCN and dispersion dC6/dCN
    // After thread completion, this is summed across threads and dcn chain rule applied
    Vector m_dEdcn;

    // Claude Generated (Mar 2026): Bond-only dE/dCN for per-component gradient attribution
    Vector m_dEdcn_bond;

    // Claude Generated (Feb 15, 2026): dc6dcn matrix for dispersion CN gradient
    // dc6dcn(i,j) = dC6(i,j)/dCN(i) - set from D4ParameterGenerator
    Matrix m_dc6dcn;

    // Claude Generated (February 2026): Individual energy term timing
    std::unordered_map<std::string, long long> m_term_timings;  // milliseconds

    // Claude Generated (February 2026): Helper to time energy term calculation
    template<typename Func>
    void timeEnergyTerm(const std::string& term_name, Func calculation) {
        auto start = std::chrono::high_resolution_clock::now();
        calculation();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        m_term_timings[term_name] = duration.count();
    }

    // Claude Generated (February 2026): Per-component gradient decomposition for validation
    bool m_store_gradient_components = false;
    Matrix m_gradient_bond, m_gradient_angle, m_gradient_torsion;
    Matrix m_gradient_repulsion, m_gradient_coulomb, m_gradient_dispersion;
    Matrix m_gradient_hb, m_gradient_xb;
    Matrix m_gradient_batm;   ///< BATM three-body gradient component (Claude Generated Mar 2026)
    Matrix m_gradient_atm;    ///< ATM three-body dispersion gradient component (Claude Generated Mar 2026)

    void initGradientComponents(int natoms) {
        m_gradient_bond = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_angle = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_torsion = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_repulsion = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_coulomb = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_dispersion = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_hb = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_xb = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_batm = Eigen::MatrixXd::Zero(natoms, 3);
        m_gradient_atm = Eigen::MatrixXd::Zero(natoms, 3);
    }

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
    virtual int Type() const override { return 3; }

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
