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

#include <set>

#include "src/core/molecule.h"
#include "src/core/qmdff_par.h"
#include "src/core/uff_par.h"

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

static json BondJson{
    { "type", 1 },
    { "i", 0 },
    { "j", 0 },
    { "k", 0 },
    { "distance", 1 },
    { "fc", 0 },
    { "exponent", 0 },
    { "r0_ij", 0 },
    { "r0_ik", 0 }
};

static json AngleJson{
    { "type", 1 },
    { "i", 0 },
    { "j", 0 },
    { "k", 0 },
    { "fc", 0 },
    { "r0_ij", 0 },
    { "r0_ik", 0 },
    { "theta0_ijk", 0 },
    { "C0", 0 },
    { "C1", 0 },
    { "C2", 0 }
};

static json DihedralJson{
    { "type", 1 },
    { "i", 0 },
    { "j", 0 },
    { "k", 0 },
    { "l", 0 },
    { "fc", 0 },
    { "r0_ij", 0 },
    { "r0_ik", 0 },
    { "theta0_ijk", 0 },
    { "V", 0 },
    { "n", 0 },
    { "phi0", 0 }
};

static json InversionJson{
    { "type", 1 },
    { "i", 0 },
    { "j", 0 },
    { "k", 0 },
    { "l", 0 },
    { "fc", 0 },
    { "C0", 0 },
    { "C1", 0 },
    { "C2", 0 }
};

static json vdWJson{
    { "type", 1 },
    { "i", 0 },
    { "j", 0 },
    { "C_ij", 1 },
    { "r0_ij", 0 },
};

static json EQJson{
    { "type", 1 },
    { "i", 0 },
    { "j", 0 },
    { "C_ij", 1 },
    { "r0_ij", 0 },
};

const json FFGenerator{
    { "method", "uff" },
    { "d3", 0 },
    { "d4", 0 },
    { "d3_s6", 1.0 },
    { "d3_s8", 2.7 },
    { "d3_s9", 1.0 },
    { "d3_a1", 0.45 },
    { "d3_a2", 4.0 },
    { "d3_alp", 1 },
    { "bond_scaling", 1 },
    { "angle_scaling", 1 },
    { "inversion_scaling", 1 },
    { "vdw_scaling", 1 },
    { "rep_scaling", 1 },
    { "dihedral_scaling", 1 },
    { "coulomb_scaling", 1 },
    { "h4", 0 },
    { "hh", 0 },
    { "h4_oh_o", 2.32 },
    { "h4_oh_n", 3.10 },
    { "h4_nh_o", 1.07 },
    { "h4_nh_n", 2.01 },
    { "h4_wh_o", 0.42 },
    { "h4_nh4", 3.61 },
    { "h4_coo", 1.41 },
    { "hh_rep_k", 0.42 },
    { "hh_rep_e", 12.7 },
    { "hh_rep_r0", 2.3 },
    { "bond_force", 1.0584 / 7.25 * 1.8 },
    { "angle_force", 1.0584 / 7.25 },
    { "torsion_force", 1 / 627.503 },
    { "inversion_force", 1 / 627.503 },
    { "vdw_force", 1 / 627.503 },
    { "h4_scaling", 0 },
    { "hh_scaling", 0 },
    { "e0", 0 }
};

class ForceFieldGenerator {
public:
    ForceFieldGenerator(const json& controller);

    void setMolecule(const Molecule& molecule);
    void Generate(const std::vector<std::pair<int, int>>& formed_bonds = std::vector<std::pair<int, int>>());
    json getParameter();

private:
    double UFFBondRestLength(int i, int j, double order);
    void AssignUffAtomTypes();

    void setBonds(const TContainer& bonds);
    void setAngles();
    void setDihedrals();
    void setInversions();
    void setvdWs();

    json Bonds() const;
    json Angles() const;
    json Dihedrals() const;
    json Inversions() const;
    json vdWs() const;
    json writeUFF();

    Molecule m_molecule;
    Matrix m_topo, m_geometry, m_distance;

    std::vector<std::vector<int>> m_stored_bonds;
    std::vector<std::vector<int>> m_identified_rings;
    std::vector<int> m_atom_types, m_coordination;
    std::vector<std::set<int>> m_ignored_vdw;
    std::vector<json> m_bonds, m_angles, m_dihedrals, m_inversions, m_vdws, m_eqs;
    double m_uff_bond_force = 1.0584 /* in Eh kcal/mol = 664.12 */, m_uff_angle_force = 1.0584 /* in Eh kcal/mol = 664.12 */, m_uff_dihedral_force = 1, m_uff_inversion_force = 1, m_vdw_force = 1, m_scaling = 1.4;

    double m_au = 1;

    int m_ff_type = 1;
    std::string m_method = "uff";
    json m_parameter;
};
