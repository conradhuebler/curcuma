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

#include "forcefieldthread.h"

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

#include <functional>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

static const json FFJson = {
    { "threads", 1 },
    { "gradient", 1 }
};

class ForceField {

public:
    ForceField(const json& controller);
    ~ForceField();

    // inline void setAtomCount(int atom) { m_natoms = atom; }
    inline void setAtomTypes(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
        m_natoms = atom_types.size();
    }

    // Claude Generated: Temporary method for EnergyCalculator compatibility
    // TODO: Eventually merge QMInterface and ForceField into unified interface
    void setMolecule(const Mol& mol);
    void UpdateGeometry(const Matrix& geometry);
    inline void UpdateGeometry(const double* coord);
    inline void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);

    double Calculate(bool gradient = true);

    Matrix Gradient() const { return m_gradient; }

    // Claude Generated: Energy component getters for regression testing (Nov 2025)
    inline double BondEnergy() const { return m_bond_energy; }
    inline double AngleEnergy() const { return m_angle_energy; }
    inline double DihedralEnergy() const { return m_dihedral_energy; }
    inline double InversionEnergy() const { return m_inversion_energy; }
    inline double VdWEnergy() const { return m_vdw_energy; }
    inline double RepulsionEnergy() const { return m_rep_energy; }
    inline double HHEnergy() const { return m_gfnff_repulsion; }  // Claude Generated (Dec 2025): GFN-FF repulsion energy
    inline double DispersionEnergy() const { return m_dispersion_energy; }
    inline double D3Energy() const { return m_d3_energy; }  // Claude Generated (Jan 2, 2026): D3 dispersion energy
    inline double D4Energy() const { return m_d4_energy; }  // Claude Generated (Jan 2, 2026): D4 dispersion energy
    inline double CoulombEnergy() const { return m_coulomb_energy; }
    inline double ElectrostatEnergy() const { return m_eq_energy; }
    inline double HydrogenBondEnergy() const { return m_energy_hbond; }   // Claude Generated (2025): Phase 5
    inline double HalogenBondEnergy() const { return m_energy_xbond; }    // Claude Generated (2025): Phase 5
    inline double ATMEnergy() const { return m_atm_energy; }        // Claude Generated (December 2025): ATM energy
    inline double BatmEnergy() const { return m_batm_energy; }       // Claude Generated (Jan 17, 2026): Batm energy

    void setParameter(const json& parameter);
    void setParameterFile(const std::string& file);

    // Parameter caching functions for all FF methods (UFF, GFN-FF, QMDFF, etc.)
    bool saveParametersToFile(const std::string& filename) const;
    bool loadParametersFromFile(const std::string& filename);
    json exportCurrentParameters() const;
    bool hasParameters() const { return !m_parameters.empty(); }

    // Auto-parameter file management: input.xyz -> input.param.json
    bool tryLoadAutoParameters(const std::string& method);
    bool autoSaveParameters() const;
    static std::string generateParameterFileName(const std::string& geometry_file);
    void setParameterCaching(bool enable) { m_enable_caching = enable; }

    // Phase 5A: Distribute EEQ charges to all threads for fqq calculation (Claude Generated Nov 2025)
    void distributeEEQCharges(const Vector& charges);

    // Get cached EEQ charges for GFN-FF cache restoration (Claude Generated Dec 2025)
    const Vector& getCachedEEQCharges() const { return m_eeq_charges; }

    // Claude Generated (Jan 18, 2026): Distribute D3 CN to all threads for dynamic r0 calculation
    // Reference: Fortran gfnff_engrad.F90:432 - CN recalculated at each energy evaluation
    void distributeD3CN(const Vector& d3_cn);

    // Get cached D3 CN for validation (Claude Generated Jan 18, 2026)
    const Vector& getCachedD3CN() const { return m_d3_cn; }

    // Claude Generated (Feb 1, 2026): Distribute CN, CNF, and CN derivatives for Coulomb gradients
    // Reference: Fortran gfnff_engrad.F90:418-422 - charge derivative via CN
    void distributeCNandDerivatives(const Vector& cn, const Vector& cnf,
                                     const std::vector<Matrix>& dcn);

    Eigen::MatrixXd NumGrad();

    // Claude Generated: Parameter analysis functionality
    void printParameterSummary() const;

private:
    void AutoRanges();
    void setBonds(const json& bonds);
    void setAngles(const json& angles);
    void setDihedrals(const json& dihedrals);
    void setInversions(const json& inversions);
    void setESPs(const json& esps);

    // Claude Generated: CG parameter generation methods
    void generateCGParameters(const json& cg_config);
    Eigen::Vector3d getCGShapeForAtom(int atom_index, const json& config);
    Eigen::Vector3d getCGOrientationForAtom(int atom_index, const json& config);

    std::vector<ForceFieldThread*> m_stored_threads;
    CxxThreadPool* m_threadpool;
    void setvdWs(const json& vdws);

    // Phase 4.2: GFN-FF pairwise non-bonded parameter setters (Claude Generated 2025)
    void setGFNFFDispersions(const json& dispersions);
    void setD4Dispersions(const json& dispersions);  // Claude Generated - Dec 25, 2025: Native D4 dispersion
    void setGFNFFBondedRepulsions(const json& repulsions);
    void setGFNFFNonbondedRepulsions(const json& repulsions);
    void setGFNFFCoulombs(const json& coulombs);

    // Phase 3: GFN-FF hydrogen bond and halogen bond parameter setters (Claude Generated 2025)
    void setGFNFFHydrogenBonds(const json& hbonds);
    void setGFNFFHalogenBonds(const json& xbonds);

    // ATM three-body dispersion parameter setter (Claude Generated 2025)
    void setATMTriples(const json& triples);

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    // GFN-FF bonded ATM (batm) parameter setter for 1,4-pairs
    void setGFNFFBatms(const json& batms);

    // Claude Generated: Energy component storage for regression testing (Nov 2025)
    double m_bond_energy = 0.0;
    double m_angle_energy = 0.0;
    double m_dihedral_energy = 0.0;
    double m_inversion_energy = 0.0;
    double m_vdw_energy = 0.0;
    double m_rep_energy = 0.0;
    double m_gfnff_repulsion = 0.0;  // Claude Generated (Dec 2025): GFN-FF repulsion energy (standard exponential repulsion)
    double m_eq_energy = 0.0;
    double m_dispersion_energy = 0.0;
    double m_coulomb_energy = 0.0;
    double m_energy_hbond = 0.0;    // Claude Generated (2025): Phase 5 - Hydrogen bond energy
    double m_energy_xbond = 0.0;    // Claude Generated (2025): Phase 5 - Halogen bond energy
    double m_atm_energy = 0.0;      // Claude Generated (December 2025): ATM three-body dispersion energy
    double m_batm_energy = 0.0;     // Claude Generated (January 17, 2026): Batm three-body dispersion energy
    double m_d3_energy = 0.0;       // Claude Generated (Jan 2, 2026): D3 dispersion energy (for UFF-D3 and GFN-FF)
    double m_d4_energy = 0.0;       // Claude Generated (Jan 2, 2026): D4 dispersion energy (for GFN-FF)

    Matrix m_geometry, m_gradient;
    std::vector<int> m_atom_types;
    std::string m_method = "uff";
    StringList m_uff_methods = { "uff", "uff-d3" };  // Claude Generated (December 19, 2025): uff-d3 already included
    StringList m_qmdff_methods = { "qmdff", "quff" };
    double m_e0 = 0;
    int m_natoms = 0;
    int m_threads = 1;
    int m_gradient_type = 1;
    std::vector<Bond> m_bonds;
    std::vector<Angle> m_angles;
    std::vector<Dihedral> m_dihedrals;        // Primary torsions (n=3, n=2, etc.)
    std::vector<Dihedral> m_extra_dihedrals;   // Extra sp3-sp3 gauche torsions (n=1) - Claude Generated (Jan 2, 2026)
    std::vector<Inversion> m_inversions;
    std::vector<vdW> m_vdWs;
    std::vector<EQ> m_EQs;

    // Phase 4.2: GFN-FF pairwise non-bonded storage (Claude Generated 2025)
    std::vector<GFNFFDispersion> m_gfnff_dispersions;
    std::vector<GFNFFDispersion> m_d4_dispersions;  // Claude Generated - Dec 25, 2025: Native D4 dispersion pairs
    std::vector<GFNFFRepulsion> m_gfnff_bonded_repulsions;
    std::vector<GFNFFRepulsion> m_gfnff_nonbonded_repulsions;
    std::vector<GFNFFCoulomb> m_gfnff_coulombs;

    // Phase 3: GFN-FF hydrogen bond and halogen bond storage (Claude Generated 2025)
    std::vector<GFNFFHydrogenBond> m_gfnff_hbonds;
    std::vector<GFNFFHalogenBond> m_gfnff_xbonds;

    // ATM three-body dispersion storage (D3/D4)
    std::vector<ATMTriple> m_atm_triples;

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    // GFN-FF bonded ATM (batm) parameters for 1,4-pairs
    std::vector<GFNFFBatmTriple> m_gfnff_batms;

    // EEQ charges for GFN-FF (cached with parameters - Claude Generated Dec 2025)
    Vector m_eeq_charges;

    // Claude Generated (Jan 18, 2026): D3 coordination numbers for dynamic r0 calculation
    // Recalculated from current geometry at each Calculate() call for cgfnff
    Vector m_d3_cn;

    // Claude Generated (Feb 1, 2026): CN, CNF, and CN derivatives for Coulomb charge derivative gradients
    // Reference: Fortran gfnff_engrad.F90:418-422 - qtmp(i) = q(i)*cnf(i)/(2*sqrt(cn(i)))
    Vector m_cn;                    // Coordination numbers per atom
    Vector m_cnf;                   // CNF parameters per atom (for qtmp calculation)
    std::vector<Matrix> m_dcn;      // CN derivatives: dcn[dim](i,j) = dCN(j)/dr(i,dim)

    json m_parameters;
    std::string m_auto_param_file; // Auto-detected parameter file path
    bool m_enable_caching = true; // Can be disabled for multi-threading
    bool m_in_setParameter = false; // Claude Generated: Recursive guard for setParameter()
};
