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
    inline double HHEnergy() const { return m_hh_energy; }  // Claude Generated (Dec 2025, Session 9): GFN-FF HH repulsion
    inline double DispersionEnergy() const { return m_dispersion_energy; }
    inline double CoulombEnergy() const { return m_coulomb_energy; }
    inline double ElectrostatEnergy() const { return m_eq_energy; }
    inline double HydrogenBondEnergy() const { return m_energy_hbond; }   // Claude Generated (2025): Phase 5
    inline double HalogenBondEnergy() const { return m_energy_xbond; }    // Claude Generated (2025): Phase 5

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
    void setGFNFFRepulsions(const json& repulsions);
    void setGFNFFCoulombs(const json& coulombs);

    // Phase 3: GFN-FF hydrogen bond and halogen bond parameter setters (Claude Generated 2025)
    void setGFNFFHydrogenBonds(const json& hbonds);
    void setGFNFFHalogenBonds(const json& xbonds);

    // Claude Generated: Energy component storage for regression testing (Nov 2025)
    double m_bond_energy = 0.0;
    double m_angle_energy = 0.0;
    double m_dihedral_energy = 0.0;
    double m_inversion_energy = 0.0;
    double m_vdw_energy = 0.0;
    double m_rep_energy = 0.0;
    double m_hh_energy = 0.0;      // Claude Generated (Dec 2025, Session 9): GFN-FF HH repulsion energy
    double m_eq_energy = 0.0;
    double m_dispersion_energy = 0.0;
    double m_coulomb_energy = 0.0;
    double m_energy_hbond = 0.0;    // Claude Generated (2025): Phase 5 - Hydrogen bond energy
    double m_energy_xbond = 0.0;    // Claude Generated (2025): Phase 5 - Halogen bond energy

    Matrix m_geometry, m_gradient;
    std::vector<int> m_atom_types;
    std::string m_method = "uff";
    StringList m_uff_methods = { "uff", "uff-d3" };
    StringList m_qmdff_methods = { "qmdff", "quff" };
    double m_e0 = 0;
    int m_natoms = 0;
    int m_threads = 1;
    int m_gradient_type = 1;
    std::vector<Bond> m_bonds;
    std::vector<Angle> m_angles;
    std::vector<Dihedral> m_dihedrals;
    std::vector<Inversion> m_inversions;
    std::vector<vdW> m_vdWs;
    std::vector<EQ> m_EQs;

    // Phase 4.2: GFN-FF pairwise non-bonded storage (Claude Generated 2025)
    std::vector<GFNFFDispersion> m_gfnff_dispersions;
    std::vector<GFNFFRepulsion> m_gfnff_repulsions;
    std::vector<GFNFFCoulomb> m_gfnff_coulombs;

    // Phase 3: GFN-FF hydrogen bond and halogen bond storage (Claude Generated 2025)
    std::vector<GFNFFHydrogenBond> m_gfnff_hbonds;
    std::vector<GFNFFHalogenBond> m_gfnff_xbonds;

    json m_parameters;
    std::string m_auto_param_file; // Auto-detected parameter file path
    bool m_enable_caching = true; // Can be disabled for multi-threading
    bool m_in_setParameter = false; // Claude Generated: Recursive guard for setParameter()
};
