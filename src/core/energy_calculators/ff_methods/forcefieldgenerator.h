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
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include <set>

#include "d3param_generator.h"
#include "d4param_generator.h"
#include "qmdff_par.h"
#include "src/core/molecule.h"
#include "uff_par.h"

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
    { "q_i", 0 },
    { "q_j", 0 },
    { "epsilon", 1 }
};

// Claude Generated 2025: ForceField Generator Parameter Registry - replaces const FFGenerator JSON
BEGIN_PARAMETER_DEFINITION(forcefield)
    // Method Selection
    PARAM(method, String, "uff", "Force field method (uff, uff-d3, d3, qmdff).", "Method", {})
    PARAM(d3_preset, String, "pbe0", "D3 functional preset for d3-only method (pbe0, blyp, b3lyp, tpss, pbe, bp86, gfnff).", "D3-Only", {})

    // Dispersion Corrections (D3/D4)
    // NOTE (Claude Generated December 21, 2025): D3 parameters removed - UFF-D3 uses D3ParameterGenerator directly
    // D3ParameterGenerator has factory methods with correct parameter sets for all functionals
    PARAM(d3, Int, 0, "Enable DFT-D3 dispersion correction (0=off, 1=on).", "Dispersion", {})
    PARAM(d4, Int, 0, "Enable DFT-D4 dispersion correction (0=off, 1=on).", "Dispersion", {})

    // Force Field Term Scaling
    PARAM(bond_scaling, Double, 1.0, "Scaling factor for bond stretch terms.", "Scaling", {})
    PARAM(angle_scaling, Double, 1.0, "Scaling factor for angle bend terms.", "Scaling", {})
    PARAM(dihedral_scaling, Double, 1.0, "Scaling factor for torsion/dihedral terms.", "Scaling", {})
    PARAM(inversion_scaling, Double, 1.0, "Scaling factor for inversion/improper terms.", "Scaling", {})
    PARAM(vdw_scaling, Double, 1.0, "Scaling factor for van der Waals terms.", "Scaling", {})
    PARAM(rep_scaling, Double, 1.0, "Scaling factor for repulsion terms.", "Scaling", {})
    PARAM(coulomb_scaling, Double, 1.0, "Scaling factor for electrostatic/Coulomb terms.", "Scaling", {})

    // H4 Hydrogen Bonding Corrections
    PARAM(h4, Int, 1, "Enable H4 hydrogen bonding correction (0=off, 1=on).", "HBond", {})
    PARAM(hh, Int, 1, "Enable H-H repulsion correction (0=off, 1=on).", "HBond", {})
    PARAM(h4_oh_o, Double, 2.32, "H4 correction for O-H...O hydrogen bonds (kJ/mol).", "HBond", {})
    PARAM(h4_oh_n, Double, 3.10, "H4 correction for O-H...N hydrogen bonds (kJ/mol).", "HBond", {})
    PARAM(h4_nh_o, Double, 1.07, "H4 correction for N-H...O hydrogen bonds (kJ/mol).", "HBond", {})
    PARAM(h4_nh_n, Double, 2.01, "H4 correction for N-H...N hydrogen bonds (kJ/mol).", "HBond", {})
    PARAM(h4_wh_o, Double, 0.42, "H4 correction for water H...O interactions (kJ/mol).", "HBond", {})
    PARAM(h4_nh4, Double, 3.61, "H4 correction for NH4+ interactions (kJ/mol).", "HBond", {})
    PARAM(h4_coo, Double, 1.41, "H4 correction for carboxylate COO- interactions (kJ/mol).", "HBond", {})
    PARAM(hh_rep_k, Double, 0.42, "H-H repulsion force constant k.", "HBond", {})
    PARAM(hh_rep_e, Double, 12.7, "H-H repulsion well depth epsilon (kJ/mol).", "HBond", {})
    PARAM(hh_rep_r0, Double, 2.3, "H-H repulsion equilibrium distance r0 (Angstrom).", "HBond", {})
    PARAM(h4_scaling, Double, 0.0, "Global scaling factor for all H4 corrections.", "HBond", {})
    PARAM(hh_scaling, Double, 0.0, "Global scaling factor for H-H repulsion.", "HBond", {})

    // Force Constants (unit conversions for UFF)
    PARAM(bond_force, Double, 0.26306, "Bond stretch force constant (Hartree/Angstrom^2).", "ForceConstants", {})
    PARAM(angle_force, Double, 0.14595, "Angle bend force constant (Hartree/radian^2).", "ForceConstants", {})
    PARAM(torsion_force, Double, 0.001593, "Torsion/dihedral force constant (Hartree).", "ForceConstants", {})
    PARAM(inversion_force, Double, 0.001593, "Inversion/improper force constant (Hartree).", "ForceConstants", {})
    PARAM(vdw_force, Double, 0.001593, "Van der Waals force constant (Hartree).", "ForceConstants", {})

    // Other Parameters
    PARAM(energy_offset, Double, 0.0, "Constant energy offset to add to calculated energy (Hartree).", "General", {"e0"})
END_PARAMETER_DEFINITION

class ForceFieldGenerator {
public:
    ForceFieldGenerator(const ConfigManager& config);

    void setMolecule(const Mol& mol);
    void Generate(const std::vector<std::pair<int, int>>& formed_bonds = std::vector<std::pair<int, int>>());
    json getParameter();

    // D3/D4 parameter generation
    void GenerateD3Parameters();
    void GenerateD4Parameters();

    // Claude Generated (December 19, 2025): UFF-D3 hybrid method
    json GenerateUFFD3Parameters();

    // Claude Generated (December 21, 2025): D3-only dispersion method
    json GenerateD3OnlyParameters(const std::string& preset = "pbe0");

private:
    double UFFBondRestLength(int i, int j, double order);
    void AssignUffAtomTypes();

    // Phase 3 LITE: GFN-FF parameter generation (Claude Generated Nov 2025)
    void GenerateGFNFF();

    void setBonds(const TContainer& bonds);
    void setAngles();
    void setDihedrals();
    void setInversions();
    void setNCI();
    // void setESP();

    json Bonds() const;
    json Angles() const;
    json Dihedrals() const;
    json Inversions() const;
    json vdWs() const;
    json ESPs() const;

    json writeUFF();

    Mol m_mol;
    int m_atoms;

    Matrix m_topo, m_geometry, m_distance;

    StringList m_uff_methods = { "uff", "uff-d3" };
    StringList m_qmdff_methods = { "qmdff", "quff" };
    Vector m_partial_charges;

    std::vector<std::vector<int>> m_stored_bonds;
    std::vector<std::vector<int>> m_identified_rings;
    std::vector<int> m_atom_types, m_ff_atom_types, m_coordination;
    std::vector<std::set<int>> m_ignored_vdw, m_1_4_charges;
    std::vector<json> m_bonds, m_angles, m_dihedrals, m_inversions, m_vdws, m_esps;
    double m_uff_bond_force = 1.0584 /* in Eh kcal/mol = 664.12 */, m_uff_angle_force = 1.0584 /* in Eh kcal/mol = 664.12 */, m_uff_dihedral_force = 1, m_uff_inversion_force = 1, m_vdw_force = 1, m_scaling = 1.4;

    double m_au = 1;

    int m_ff_type = 1;
    std::string m_method = "uff";
    json m_parameter;
};
