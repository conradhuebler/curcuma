/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "ff_terms.h"
#include "ff_workspace.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

/**
 * @brief UFF / UFF-D3 / QMDFF force field front-end.
 *
 * ForceField parses the interaction lists produced by ForceFieldGenerator
 * (bonds, angles, dihedrals, inversions, vdW pairs and, for uff-d3, the D3
 * dispersion pairs), owns the parameter cache (input.xyz -> input.param.json)
 * and evaluates energy and gradient through FFWorkspace. GFN-FF has its own
 * engine (class GFNFF) and does not pass through this class.
 *
 * Claude Generated (Sep 2026): the legacy ForceFieldThread engine, the CG pair
 * loop and the GFN-FF parameter intake were removed. All of it had been dead
 * since GFNFF moved to FFWorkspace and the MethodFactory stopped routing the
 * "cg" / "d3" method names here.
 */
class ForceField {

public:
    ForceField(const json& controller);
    ~ForceField();

    inline void setAtomTypes(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
        m_natoms = atom_types.size();
    }

    // Claude Generated: Temporary method for EnergyCalculator compatibility
    // TODO: Eventually merge QMInterface and ForceField into unified interface
    void setMolecule(const Mol& mol);
    void UpdateGeometry(const Matrix& geometry);
    void UpdateGeometry(const double* coord);
    void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);

    /// Energy in Hartree (and gradient in Hartree/Bohr if requested) from FFWorkspace.
    double Calculate(bool gradient = true);

    Matrix Gradient() const { return m_gradient; }

    // Claude Generated: Energy component getters for regression testing (Nov 2025)
    inline double BondEnergy() const { return m_bond_energy; }
    inline double AngleEnergy() const { return m_angle_energy; }
    inline double DihedralEnergy() const { return m_dihedral_energy; }
    inline double InversionEnergy() const { return m_inversion_energy; }
    inline double VdWEnergy() const { return m_vdw_energy; }
    inline double RepulsionEnergy() const { return m_rep_energy; }
    inline double DispersionEnergy() const { return m_dispersion_energy; }
    inline double D3Energy() const { return m_d3_energy; }  // Claude Generated (Jan 2, 2026): D3 dispersion energy (uff-d3)

    void setParameter(const json& parameter);
    void setParameterFile(const std::string& file);

    // Parameter caching functions (UFF, UFF-D3, QMDFF)
    bool saveParametersToFile(const std::string& filename) const;
    bool loadParametersFromFile(const std::string& filename);
    json exportCurrentParameters() const;
    bool hasParameters() const { return !m_parameters.empty(); }

    // Auto-parameter file management: input.xyz -> input.param.json
    bool tryLoadAutoParameters(const std::string& method);
    bool autoSaveParameters() const;
    static std::string generateParameterFileName(const std::string& geometry_file);
    void setParameterCaching(bool enable) { m_enable_caching = enable; }

    Eigen::MatrixXd NumGrad();

    // Claude Generated: Parameter analysis functionality
    void printParameterSummary() const;

private:
    // Coarse-grained method (Sep 2026, restored from the removed thread engine): one vdW
    // entry of type 3 per CG-CG pair, parameters from cg_default / cg_per_atom /
    // pair_interactions of the controller (e.g. -load_ff_json FILE).
    void generateCGParameters(const json& cg_config);
    Eigen::Vector3d getCGShapeForAtom(int atom_index, const json& config) const;
    Eigen::Vector3d getCGOrientationForAtom(int atom_index, const json& config) const;
    void setBonds(const json& bonds);
    void setAngles(const json& angles);
    void setDihedrals(const json& dihedrals);
    void setInversions(const json& inversions);
    void setvdWs(const json& vdws);
    void setESPs(const json& esps);

    /// Claude Generated (Sep 2026): parse "d3_dispersion_pairs" (D3ParameterGenerator output,
    /// forwarded by ForceFieldGenerator for uff-d3) into the GFNFFDispersion pair layout that
    /// FFWorkspace evaluates with dispersion_method = "d3".
    void setD3DispersionPairs(const json& pairs);

    /// Claude Generated (Sep 2026): build the FFWorkspace from the parsed interaction lists.
    void buildWorkspace();

    CxxThreadPool* m_threadpool;

    // Claude Generated: Energy component storage for regression testing (Nov 2025)
    double m_bond_energy = 0.0;
    double m_angle_energy = 0.0;
    double m_dihedral_energy = 0.0;
    double m_inversion_energy = 0.0;
    double m_vdw_energy = 0.0;
    double m_rep_energy = 0.0;
    double m_dispersion_energy = 0.0;
    double m_d3_energy = 0.0;       // Claude Generated (Jan 2, 2026): D3 dispersion energy (uff-d3)

    GeoGradMatrix m_geometry, m_gradient;  // WP-G: RowMajor N×3 hot data
    std::vector<int> m_atom_types;
    std::string m_method = "uff";
    double m_e0 = 0;
    int m_natoms = 0;
    int m_threads = 1;
    std::vector<Bond> m_bonds;
    std::vector<Angle> m_angles;
    std::vector<Dihedral> m_dihedrals;
    std::vector<Inversion> m_inversions;
    std::vector<vdW> m_vdWs;
    std::vector<EQ> m_EQs;

    /// uff-d3 pairwise D3 dispersion terms (GFNFFDispersion layout consumed by FFWorkspace)
    std::vector<GFNFFDispersion> m_d3_dispersions;

    json m_parameters;
    std::string m_auto_param_file; // Auto-detected parameter file path
    bool m_enable_caching = true; // Can be disabled for multi-threading
    bool m_in_setParameter = false; // Claude Generated: Recursive guard for setParameter()

    // Claude Generated (March 2026): FFWorkspace evaluates all UFF/QMDFF terms
    std::unique_ptr<FFWorkspace> m_workspace;
};
