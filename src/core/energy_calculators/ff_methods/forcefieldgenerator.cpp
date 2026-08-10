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

#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "src/core/topology.h"
#include "src/tools/general.h"

#include "forcefieldfunctions.h"

// Claude Generated (Aug 2026): QMDFF term generation (port of xtb src/qmdff.f90)
#include "d3param_generator.h"
#include "qmdff_data.h"
#include "qmdff_terms.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"

#include <chrono>
#include <map>
#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>

#include "json.hpp"
using json = nlohmann::json;

#include "forcefieldgenerator.h"

ForceFieldGenerator::ForceFieldGenerator(const ConfigManager& config)
{
    // Claude Generated: Type-safe parameter access via ConfigManager (Phase 3B)
    m_method = config.get<std::string>("method", "uff");

    // Build parameter json from ConfigManager for backward compatibility
    m_parameter = json();
    m_parameter["method"] = config.get<std::string>("method", "uff");
    m_parameter["d3"] = config.get<int>("d3", 0);
    m_parameter["d4"] = config.get<int>("d4", 0);
    m_parameter["d3_s6"] = config.get<double>("d3_s6", 1.0);
    m_parameter["d3_s8"] = config.get<double>("d3_s8", 2.7);
    m_parameter["d3_s9"] = config.get<double>("d3_s9", 1.0);
    m_parameter["d3_a1"] = config.get<double>("d3_a1", 0.45);
    m_parameter["d3_a2"] = config.get<double>("d3_a2", 4.0);
    m_parameter["d3_alp"] = config.get<double>("d3_alp", 1.0);
    m_parameter["d3method"] = config.get<std::string>("d3method", "none"); // Claude Generated: Added d3method parameter support
    m_parameter["bond_scaling"] = config.get<double>("bond_scaling", 1.0);
    m_parameter["angle_scaling"] = config.get<double>("angle_scaling", 1.0);
    m_parameter["dihedral_scaling"] = config.get<double>("dihedral_scaling", 1.0);
    m_parameter["inversion_scaling"] = config.get<double>("inversion_scaling", 1.0);
    m_parameter["vdw_scaling"] = config.get<double>("vdw_scaling", 1.0);
    m_parameter["rep_scaling"] = config.get<double>("rep_scaling", 1.0);
    m_parameter["coulomb_scaling"] = config.get<double>("coulomb_scaling", 1.0);
    m_parameter["h4"] = config.get<int>("h4", 1);
    m_parameter["hh"] = config.get<int>("hh", 1);
    m_parameter["h4_oh_o"] = config.get<double>("h4_oh_o", 2.32);
    m_parameter["h4_oh_n"] = config.get<double>("h4_oh_n", 3.10);
    m_parameter["h4_nh_o"] = config.get<double>("h4_nh_o", 1.07);
    m_parameter["h4_nh_n"] = config.get<double>("h4_nh_n", 2.01);
    m_parameter["h4_wh_o"] = config.get<double>("h4_wh_o", 0.42);
    m_parameter["h4_nh4"] = config.get<double>("h4_nh4", 3.61);
    m_parameter["h4_coo"] = config.get<double>("h4_coo", 1.41);
    m_parameter["hh_rep_k"] = config.get<double>("hh_rep_k", 0.42);
    m_parameter["hh_rep_e"] = config.get<double>("hh_rep_e", 12.7);
    m_parameter["hh_rep_r0"] = config.get<double>("hh_rep_r0", 2.3);
    m_parameter["bond_force"] = config.get<double>("bond_force", 1.0584 / 7.25 * 1.8);
    m_parameter["angle_force"] = config.get<double>("angle_force", 1.0584 / 7.25);
    m_parameter["torsion_force"] = config.get<double>("torsion_force", 1.0 / 627.503);
    m_parameter["inversion_force"] = config.get<double>("inversion_force", 1.0 / 627.503);
    m_parameter["vdw_force"] = config.get<double>("vdw_force", 1.0 / 627.503);
    m_parameter["h4_scaling"] = config.get<double>("h4_scaling", 0.0);
    m_parameter["hh_scaling"] = config.get<double>("hh_scaling", 0.0);
    m_parameter["e0"] = config.get<double>("energy_offset", 0.0);
    m_parameter["verbosity"] = config.get<int>("verbosity", 0);  // Claude Generated (2025): Store verbosity for ForceFieldGenerator logging
}

void ForceFieldGenerator::setMolecule(const Mol& mol)
{
    m_mol = mol;
    m_geometry = mol.m_geometry;
    m_atoms = mol.m_number_atoms;
    m_atom_types = mol.m_atoms;
    m_ff_atom_types = mol.m_atoms;
    m_partial_charges = mol.m_partial_charges;
}

// Claude Generated: Enhanced Generate method with detailed timing and status output
void ForceFieldGenerator::Generate(const std::vector<std::pair<int, int>>& formed_bonds)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Level 1+: Parameter generation start
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Starting force field parameter generation");
        CurcumaLogger::param("atoms", m_atoms);
        CurcumaLogger::param("method", m_method);
    }

    m_coordination = std::vector<int>(m_atoms, 0);
    m_topo = Eigen::MatrixXd::Zero(m_atoms, m_atoms);
    m_distance = Eigen::MatrixXd::Zero(m_atoms, m_atoms);

    auto topology_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Topology initialization
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing topology and distance matrices");
    }

    TContainer bonds;
    m_scaling = 1.4;
    if (m_mol.m_has_topology && formed_bonds.empty()) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Using persistent topology matrix for force field generation");
        }
        for (int i = 0; i < m_atoms; ++i) {
            m_stored_bonds.push_back(std::vector<int>());
            m_ignored_vdw.push_back(std::set<int>({ i }));
            m_1_4_charges.push_back(std::set<int>({ i }));
        }
        for (int i = 0; i < m_atoms; ++i) {
            for (int j = i + 1; j < m_atoms; ++j) {
                double distance = (m_geometry.row(i) - m_geometry.row(j)).norm() * m_au;
                m_distance(i, j) = distance;
                m_distance(j, i) = distance;
                if (m_mol.m_topology(i, j) == 1) {
                    if (bonds.insert({ i, j })) {
                        m_coordination[i]++;
                        m_coordination[j]++;
                        m_stored_bonds[i].push_back(j);
                        m_stored_bonds[j].push_back(i);
                        m_ignored_vdw[i].insert(j);
                        m_ignored_vdw[j].insert(i);
                    }
                    m_topo(i, j) = 1;
                    m_topo(j, i) = 1;
                }
            }
        }
    } else if (formed_bonds.size() == 0) {
        for (int i = 0; i < m_atoms; ++i) {
            m_stored_bonds.push_back(std::vector<int>());
            m_ignored_vdw.push_back(std::set<int>({ i }));
            m_1_4_charges.push_back(std::set<int>({ i }));
            for (int j = 0; j < m_atoms && m_stored_bonds[i].size() < CoordinationNumber[m_mol.m_atoms[i]]; ++j) {
                if (i == j)
                    continue;
                double x_i = m_geometry(i, 0) * m_au;
                double x_j = m_geometry(j, 0) * m_au;

                double y_i = m_geometry(i, 1) * m_au;
                double y_j = m_geometry(j, 1) * m_au;

                double z_i = m_geometry(i, 2) * m_au;
                double z_j = m_geometry(j, 2) * m_au;

                double r_ij = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));
                m_distance(i, j) = r_ij;
                m_distance(j, i) = r_ij;
                if (r_ij <= (Elements::CovalentRadius[m_mol.m_atoms[i]] + Elements::CovalentRadius[m_mol.m_atoms[j]]) * m_scaling * m_au) {
                    if (bonds.insert({ std::min(i, j), std::max(i, j) })) {
                        m_coordination[i]++;
                        m_stored_bonds[i].push_back(j);
                        m_ignored_vdw[i].insert(j);
                    }
                    m_topo(i, j) = 1;
                    m_topo(j, i) = 1;
                }
            }
        }
    } else {
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_stored_bonds.push_back(std::vector<int>());
            m_ignored_vdw.push_back(std::set<int>({ i }));
        }
        for (const std::pair<int, int>& bond : formed_bonds) {

            int i = bond.first - 1;
            int j = bond.second - 1;
            if (bonds.insert({ std::min(i, j), std::max(i, j) })) {
                m_coordination[i]++;
                m_coordination[j]++;

                m_stored_bonds[i].push_back(j);
                m_stored_bonds[j].push_back(i);

                m_ignored_vdw[i].insert(j);
                m_ignored_vdw[j].insert(i);
            }
            m_topo(i, j) = 1;
            m_topo(j, i) = 1;
        }
    }

    auto topology_end = std::chrono::high_resolution_clock::now();
    auto topology_duration = std::chrono::duration_cast<std::chrono::milliseconds>(topology_end - topology_start);
    // Level 2+: Bond detection results
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("Bond detection completed");
        CurcumaLogger::param("bonds_found", static_cast<int>(bonds.Storage().size()));
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("topology_time", fmt::format("{} ms", topology_duration.count()));
        }
    }

    auto atomtype_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Atom type assignment
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Assigning UFF atom types");
    }
    AssignUffAtomTypes();
    // if (m_rings)
    // std::cout << "Crude ring finding method ... " << std::endl;
    int maxsize = m_atom_types.size();
    int maxcache = m_atom_types.size();
    if (m_atom_types.size() > 1000)
        maxcache = 5;
    m_identified_rings = Topology::FindRings(m_stored_bonds, m_atom_types.size(), maxsize, maxcache);

    auto atomtype_end = std::chrono::high_resolution_clock::now();
    auto atomtype_duration = std::chrono::duration_cast<std::chrono::milliseconds>(atomtype_end - atomtype_start);
    // Level 2+: Atom type assignment results
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("Atom type assignment and ring detection completed");
        CurcumaLogger::param("rings_found", static_cast<int>(m_identified_rings.size()));
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("atomtype_time", fmt::format("{} ms", atomtype_duration.count()));
        }
    }

    // Level 2+: Force field method setup
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Setting up force field method: {}", m_method));
    }
    if (m_method.compare("uff") == 0) {
        m_ff_type = 1;
        // Claude Generated: Check for d3method parameter to enable D3/D4 dispersion
        std::string d3method = m_parameter["d3method"].get<std::string>();
        if (d3method != "none" && d3method != "free") {
            m_parameter["d3"] = 1;
            m_parameter["vdw_scaling"] = 0;
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info("UFF D3/D4 dispersion enabled via d3method parameter: " + d3method);
            }
        }
    } else if (m_method.compare("uff-d3") == 0) {
        m_ff_type = 1;
        m_parameter["d3"] = 1;
        m_parameter["vdw_scaling"] = 0;
    } else if (m_method.compare("quff") == 0) {
        m_ff_type = 3;

    } else if (m_method.compare("qmdff") == 0) {
        m_ff_type = 2;
        // Claude Generated (Aug 2026): NO separate D3/D4 pass here. QMDFF carries its own
        // dispersion inside the non-covalent block (setQMDFFTerms), and the FFWorkspace
        // runs QMDFF with dispersion_enabled=false — the former `m_parameter["d3"] = 1`
        // therefore generated a full D3 parameter set that was then discarded.
        m_parameter["vdw_scaling"] = 0;
    } else if (m_method.compare("gfnff") == 0) {
        m_ff_type = 4;  // GFN-FF type (Phase 3 integration)
        // GFN-FF uses built-in dispersion/repulsion instead of D3/D4
        m_parameter["vdw_scaling"] = 0;
    } else if (m_method.compare("d3") == 0) {
        m_ff_type = 5;  // D3-only type (Claude Generated December 21, 2025)
        m_parameter["d3"] = 1;
        m_parameter["vdw_scaling"] = 0;

        // D3-only: Skip ALL bonded terms
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::info("D3-only method: skipping bonded term generation");
        }

        // Get preset from m_parameter (default: pbe0)
        std::string preset = m_parameter.value("d3_preset", "pbe0");

        // Generate D3 parameters only
        json d3_only_params = GenerateD3OnlyParameters(preset);
        m_parameter = MergeJson(m_parameter, d3_only_params);

        // Early return - no bonded terms needed
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success("D3-only parameter generation complete");
        }
        return;  // Skip bonds, angles, dihedrals, inversions
    }

    m_uff_bond_force = m_parameter["bond_force"];
    m_uff_angle_force = m_parameter["angle_force"];
    m_uff_dihedral_force = m_parameter["torsion_force"];
    m_uff_inversion_force = m_parameter["inversion_force"];
    m_vdw_force = m_parameter["vdw_force"];

    // Generate force field parameters with timing
    auto bonds_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Bond parameter generation
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating bond parameters");
    }
    setBonds(bonds);
    auto bonds_end = std::chrono::high_resolution_clock::now();
    auto bonds_duration = std::chrono::duration_cast<std::chrono::milliseconds>(bonds_end - bonds_start);

    auto angles_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Angle parameter generation
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating angle parameters");
    }
    setAngles();
    auto angles_end = std::chrono::high_resolution_clock::now();
    auto angles_duration = std::chrono::duration_cast<std::chrono::milliseconds>(angles_end - angles_start);

    auto dihedrals_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Dihedral parameter generation
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating dihedral parameters");
    }
    setDihedrals();
    auto dihedrals_end = std::chrono::high_resolution_clock::now();
    auto dihedrals_duration = std::chrono::duration_cast<std::chrono::milliseconds>(dihedrals_end - dihedrals_start);

    auto inversions_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Inversion parameter generation
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating inversion parameters");
    }
    setInversions();
    auto inversions_end = std::chrono::high_resolution_clock::now();
    auto inversions_duration = std::chrono::duration_cast<std::chrono::milliseconds>(inversions_end - inversions_start);

    auto nci_start = std::chrono::high_resolution_clock::now();
    // Level 2+: Non-covalent parameter generation
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating non-covalent interaction parameters");
    }
    setNCI();
    auto nci_end = std::chrono::high_resolution_clock::now();
    auto nci_duration = std::chrono::duration_cast<std::chrono::milliseconds>(nci_end - nci_start);

    // D3/D4 parameter generation if enabled
    if (m_parameter["d3"].get<int>() == 1) {
        GenerateD3Parameters();
    }
    if (m_parameter["d4"].get<int>() == 1) {
        GenerateD4Parameters();
    }

    // Claude Generated (Aug 2026): QMDFF-specific terms (torsions in QMDFF's Fourier
    // form, non-covalent pair list, hydrogen/halogen bonds). Only meaningful for the
    // QMDFF methods; see setQMDFFTerms() for the mapping and its limitations.
    if (m_ff_type == 2) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Generating QMDFF torsion / non-covalent / HB parameters");
        }
        setQMDFFTerms();
    }

    // Summary output with timing
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - start_time);

    // Level 1+: Completion summary
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::success("Force field parameter generation completed");
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("total_time", fmt::format("{} ms", total_duration.count()));
            CurcumaLogger::param("bonds_generated", static_cast<int>(m_bonds.size()));
            CurcumaLogger::param("angles_generated", static_cast<int>(m_angles.size()));
            CurcumaLogger::param("dihedrals_generated", static_cast<int>(m_dihedrals.size()));
            CurcumaLogger::param("inversions_generated", static_cast<int>(m_inversions.size()));
            CurcumaLogger::param("nci_pairs_generated", static_cast<int>(m_vdws.size()));
        }
    }

    // Level 3+: Detailed timing breakdown
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Detailed timing breakdown:");
        CurcumaLogger::param("topology_time", fmt::format("{} ms", topology_duration.count()));
        CurcumaLogger::param("atomtype_time", fmt::format("{} ms", atomtype_duration.count()));
        CurcumaLogger::param("bonds_time", fmt::format("{} ms", bonds_duration.count()));
        CurcumaLogger::param("angles_time", fmt::format("{} ms", angles_duration.count()));
        CurcumaLogger::param("dihedrals_time", fmt::format("{} ms", dihedrals_duration.count()));
        CurcumaLogger::param("inversions_time", fmt::format("{} ms", inversions_duration.count()));
        CurcumaLogger::param("nci_time", fmt::format("{} ms", nci_duration.count()));
    }
}

double ForceFieldGenerator::UFFBondRestLength(int i, int j, double n)
{
    double cRi = UFFParameters[m_atom_types[i]][cR];
    double cRj = UFFParameters[m_atom_types[j]][cR];
    double cXii = UFFParameters[m_atom_types[i]][cXi];
    double cXij = UFFParameters[m_atom_types[j]][cXi];

    double lambda = 0.13332;
    double r_BO = -lambda * (cRi + cRj) * log(n);
    double r_EN = cRi * cRj * (sqrt(cXii) - sqrt(cXij)) * (sqrt(cXii) - sqrt(cXij)) / (cRi * cXii + cRj * cXij);
    double r_0 = cRi + cRj;
    return (r_0 + r_BO - r_EN) * m_au;
}

void ForceFieldGenerator::setBonds(const TContainer& bonds)
{
    std::vector<std::pair<int, int>> bonded;
    std::vector<int> atom_types = m_mol.m_atoms;
    for (const auto& bond : bonds.Storage()) {
        if (std::find(bonded.begin(), bonded.end(), std::pair<int, int>(std::min(bond[0], bond[1]), std::max(bond[0], bond[1]))) != bonded.end())
            continue;
        bonded.push_back(std::pair<int, int>(std::min(bond[0], bond[1]), std::max(bond[0], bond[1])));
        json uffbond = BondJson;

        uffbond["i"] = bond[0];
        uffbond["j"] = bond[1];
        uffbond["type"] = m_ff_type;
        int bond_order = 1;

        if (std::find(Conjugated.cbegin(), Conjugated.cend(), m_atom_types[bond[0]]) != Conjugated.cend() && std::find(Conjugated.cbegin(), Conjugated.cend(), m_atom_types[bond[1]]) != Conjugated.cend())
            bond_order = 2;
        else if (std::find(Triples.cbegin(), Triples.cend(), m_atom_types[bond[0]]) != Triples.cend() || std::find(Triples.cbegin(), Triples.cend(), m_atom_types[bond[1]]) != Triples.cend())
            bond_order = 3;
        else
            bond_order = 1;
        double r0_ij = 0;
        if (m_ff_type == 1)
            r0_ij = UFFBondRestLength(bond[0], bond[1], bond_order);
        else
            r0_ij = m_distance(bond[0], bond[1]);
        uffbond["r0_ij"] = r0_ij;
        double cZi = UFFParameters[m_atom_types[bond[0]]][cZ];
        double cZj = UFFParameters[m_atom_types[bond[1]]][cZ];
        if (m_ff_type == 2)
            uffbond["fc"] = 0.01;
        else
            uffbond["fc"] = m_uff_bond_force * cZi * cZj / (r0_ij * r0_ij * r0_ij);

        double exp = (ka(atom_types[bond[0]]) * ka(atom_types[bond[1]]) + kEN * std::abs((AllenEN(atom_types[bond[0]]) - AllenEN(atom_types[bond[1]])) * (AllenEN(atom_types[bond[0]]) - AllenEN(atom_types[bond[1]]))));
        exp /= 4.75;
        uffbond["exponent"] = exp;
        m_bonds.push_back(uffbond);

        int i = bond[0];
        int j = bond[1];

        std::vector<int> k_bodies;
        for (auto t : m_stored_bonds[i]) {
            k_bodies.push_back(t);

            if (t == j)
                continue;
            json uffangle = AngleJson;
            uffangle["j"] = i;
            uffangle["i"] = std::min(t, j);
            uffangle["k"] = std::max(j, t);

            if (std::find(m_angles.begin(), m_angles.end(), uffangle) == m_angles.end())
                m_angles.push_back(uffangle);
            m_ignored_vdw[j].insert(i);
            m_ignored_vdw[j].insert(t);
        }

        std::vector<int> l_bodies;
        for (auto t : m_stored_bonds[j]) {
            l_bodies.push_back(t);

            if (t == i)
                continue;
            json uffangle = AngleJson;
            uffangle["j"] = j;
            uffangle["i"] = std::min(t, i);
            uffangle["k"] = std::max(i, t);
            if (std::find(m_angles.begin(), m_angles.end(), uffangle) == m_angles.end())
                m_angles.push_back(uffangle);

            m_ignored_vdw[j].insert(i);
            m_ignored_vdw[j].insert(t);
        }

        for (int k : k_bodies) {
            for (int l : l_bodies) {
                if (k == i || k == j || k == l || i == j || i == l || j == l)
                    continue;
                json uffdihedral = DihedralJson;
                uffdihedral["i"] = k;
                uffdihedral["j"] = i;
                uffdihedral["k"] = j;
                uffdihedral["l"] = l;
                if (std::find(m_dihedrals.begin(), m_dihedrals.end(), uffdihedral) == m_dihedrals.end()) {
                    m_dihedrals.push_back(uffdihedral);
                    m_ignored_vdw[i].insert(k);
                    m_ignored_vdw[i].insert(l);
                    m_ignored_vdw[j].insert(k);
                    m_ignored_vdw[j].insert(l);
                    m_ignored_vdw[k].insert(l);
                    m_ignored_vdw[l].insert(k);
                    m_1_4_charges[i].insert(k);
                }
            }
        }
        if (m_stored_bonds[i].size() == 3) {
            json inversion = InversionJson;
            inversion["i"] = i;
            inversion["j"] = m_stored_bonds[i][0];
            inversion["k"] = m_stored_bonds[i][1];
            inversion["l"] = m_stored_bonds[i][2];
            if (std::find(m_inversions.begin(), m_inversions.end(), inversion) == m_inversions.end()) {
                m_inversions.push_back(inversion);
            }
        }
        if (m_stored_bonds[j].size() == 3) {
            json inversion = InversionJson;
            inversion["i"] = j;
            inversion["j"] = m_stored_bonds[j][0];
            inversion["k"] = m_stored_bonds[j][1];
            inversion["l"] = m_stored_bonds[j][2];
            if (std::find(m_inversions.begin(), m_inversions.end(), inversion) == m_inversions.end()) {
                m_inversions.push_back(inversion);
            }
        }
    }
}

void ForceFieldGenerator::setAngles()
{
    std::vector<int> atom_types = m_atom_types;

    for (int index = 0; index < m_angles.size(); ++index) {
        int i = m_angles[index]["i"];
        int j = m_angles[index]["j"];
        int k = m_angles[index]["k"];
        if (i == j || i == k || j == k) {
            m_angles[index]["type"] = 0; // this will be set to be removed
            continue;
        }
        double f = pi / 180.0;
        double r0_ij = UFFBondRestLength(i, j, 1);
        double r0_jk = UFFBondRestLength(j, k, 1);
        double Theta0 = UFFParameters[m_atom_types[j]][cTheta0];
        double cosTheta0 = cos(Theta0 * f);
        double r0_ik = sqrt(r0_ij * r0_ij + r0_jk * r0_jk - 2. * r0_ij * r0_jk * cosTheta0);
        double param = m_uff_angle_force;
        double beta = 2.0 * param / (r0_ij * r0_jk);
        double preFactor = beta * UFFParameters[m_atom_types[j]][cZ] * UFFParameters[m_atom_types[k]][cZ] / (r0_ik * r0_ik * r0_ik * r0_ik * r0_ik);
        double rTerm = r0_ij * r0_jk;
        double inner = 3.0 * rTerm * (1.0 - cosTheta0 * cosTheta0) - r0_ik * r0_ik * cosTheta0;
        m_angles[index]["fc"] = preFactor * rTerm * inner;
        double C2 = 1 / (4 * std::max(sin(Theta0 * f) * sin(Theta0 * f), 1e-4));
        double C1 = -4 * C2 * cosTheta0;
        double C0 = C2 * (2 * cosTheta0 * cosTheta0 + 1);
        m_angles[index]["C0"] = C0;
        m_angles[index]["C1"] = C1;
        m_angles[index]["C2"] = C2;
#pragma message("TODO: Check UFF angle bending")
        Eigen::Vector3d rij = m_geometry.row(i) - m_geometry.row(j);
        Eigen::Vector3d nij = rij.normalized();
        Eigen::Vector3d rkj = m_geometry.row(k) - m_geometry.row(j);
        Eigen::Vector3d nkj = rkj.normalized();
        double costheta = nij.dot(nkj);
        m_angles[index]["theta0_ijk"] = acos(costheta); // m_molecule.CalculateAngle(i, j, k) / f; // UFF::AngleBending(m_geometry.row(j), m_geometry.row(i), m_geometry.row(k), derivate, false);
        /*
                    if(std::find(m_qmdff_methods.begin(), m_qmdff_methods.end(), m_method) != m_qmdff_methods.end())
                    {
                        json uffbond = BondJson;

                        uffbond["i"] = i;
                        uffbond["j"] = k;
                        uffbond["k"] = k;

                        uffbond["type"] = m_ff_type;

                        r0_ij = m_distance(i, j);
                        r0_ik = m_distance(i, k);

                        uffbond["r0_ij"] = r0_ik;
                        uffbond["r0_ik"] = r0_ik;


                        uffbond["fc"] = 0.01;
                        double exp =ka13 + kb13* ka(atom_types[i]) * ka(atom_types[k]);
                        bool found = false;
                        for(const auto &ring : m_identified_rings)
                        {
                            for(const auto &atom : ring)
                            {
                                if((atom == i || atom == j || atom == k) && !found)
                                {
                                    exp += k13r;
                                    found = true;
                                }
                            }
                        }
                       // exp /= 3.0;
                        uffbond["exponent"] = exp;
                        uffbond["distance"] = 2;
                        m_bonds.push_back(uffbond);
                    }
    */
    }
}

void ForceFieldGenerator::setDihedrals()
{
    for (int index = 0; index < m_dihedrals.size(); ++index) {
        int i = m_dihedrals[index]["i"];
        int j = m_dihedrals[index]["j"];
        int k = m_dihedrals[index]["k"];
        int l = m_dihedrals[index]["l"];

        m_dihedrals[index]["n"] = 2;
        double f = pi / 180.0;
        double bond_order = 1;
        m_dihedrals[index]["V"] = 2;
        m_dihedrals[index]["n"] = 3;
        m_dihedrals[index]["phi0"] = 180 * f;

        if (std::find(Conjugated.cbegin(), Conjugated.cend(), m_atom_types[k]) != Conjugated.cend() && std::find(Conjugated.cbegin(), Conjugated.cend(), m_atom_types[j]) != Conjugated.cend())
            bond_order = 2;
        else if (std::find(Triples.cbegin(), Triples.cend(), m_atom_types[k]) != Triples.cend() || std::find(Triples.cbegin(), Triples.cend(), m_atom_types[j]) != Triples.cend())
            bond_order = 3;
        else
            bond_order = 1;

        if (m_coordination[j] == 4 && m_coordination[k] == 4) // 2*sp3
        {
            m_dihedrals[index]["V"] = sqrt(UFFParameters[m_atom_types[j]][cV] * UFFParameters[m_atom_types[k]][cV]) * m_uff_dihedral_force;
            m_dihedrals[index]["phi0"] = 180 * f;
            m_dihedrals[index]["n"] = 3;
        }
        if (m_coordination[j] == 3 && m_coordination[k] == 3) // 2*sp2
        {
            m_dihedrals[index]["V"] = 5 * sqrt(UFFParameters[m_atom_types[j]][cU] * UFFParameters[m_atom_types[k]][cU]) * (1 + 4.18 * log(bond_order)) * m_uff_dihedral_force;
            m_dihedrals[index]["phi0"] = 180 * f;
            m_dihedrals[index]["n"] = 2;
        } else if ((m_coordination[j] == 4 && m_coordination[k] == 3) || (m_coordination[j] == 3 && m_coordination[k] == 4)) {
            m_dihedrals[index]["V"] = sqrt(UFFParameters[m_atom_types[j]][cV] * UFFParameters[m_atom_types[k]][cV]) * m_uff_dihedral_force;
            m_dihedrals[index]["phi0"] = 0 * f;
            m_dihedrals[index]["n"] = 6;

        } else {
            m_dihedrals[index]["V"] = 5 * sqrt(UFFParameters[m_atom_types[j]][cU] * UFFParameters[m_atom_types[k]][cU]) * (1 + 4.18 * log(bond_order)) * m_uff_dihedral_force;
            m_dihedrals[index]["phi0"] = 90 * f;
        }
    }
}

void ForceFieldGenerator::setInversions()
{
    for (int index = 0; index < m_inversions.size(); ++index) {
        const int i = m_inversions[index]["i"];
        if (m_coordination[i] != 3)
            continue;

        int j = m_inversions[index]["j"];
        int k = m_inversions[index]["k"];
        int l = m_inversions[index]["l"];

        double C0 = 0.0;
        double C1 = 0.0;
        double C2 = 0.0;
        double f = pi / 180.0;
        double kijkl = 0;

#pragma message("TODO: Check UFF inversion")
        if (6 <= m_mol.m_atoms[i] && m_mol.m_atoms[i] <= 8) {
            C0 = 1.0;
            C1 = -1.0;
            C2 = 0.0;
            kijkl = 6 * m_uff_inversion_force;
            if (m_mol.m_atoms[j] == 8 || m_mol.m_atoms[k] == 8 || m_mol.m_atoms[l] == 8)
                kijkl = 50 * m_uff_inversion_force;
        } else {
            double w0 = pi / 180.0;
            switch (m_mol.m_atoms[i]) {
            // if the central atom is phosphorous
            case 15:
                w0 *= 84.4339;
                break;

                // if the central atom is arsenic
            case 33:
                w0 *= 86.9735;
                break;

                // if the central atom is antimonium
            case 51:
                w0 *= 87.7047;
                break;

                // if the central atom is bismuth
            case 83:
                w0 *= 90.0;
                break;
            }
            C2 = 1.0;
            C1 = -4.0 * cos(w0 * f);
            C0 = -(C1 * cos(w0 * f) + C2 * cos(2.0 * w0 * f));
            kijkl = 22.0 / (C0 + C1 + C2) * m_uff_inversion_force;
        }
        m_inversions[index]["C0"] = C0;
        m_inversions[index]["C1"] = C1;
        m_inversions[index]["C2"] = C2;
        m_inversions[index]["fc"] = kijkl;
    }
}

void ForceFieldGenerator::setNCI()
{
    // std::vector<double> charges = m_mol.m_partial_charges // m_molecule.getPartialCharges();
    double q_thresh = 1e-5;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        for (int j = i + 1; j < m_atom_types.size(); ++j) {
            if (std::find(m_ignored_vdw[i].begin(), m_ignored_vdw[i].end(), j) != m_ignored_vdw[i].end() || std::find(m_ignored_vdw[j].begin(), m_ignored_vdw[j].end(), i) != m_ignored_vdw[j].end())
                continue;
            json vdW = vdWJson;
            double cDi = UFFParameters[m_atom_types[i]][cD];
            double cDj = UFFParameters[m_atom_types[j]][cD];
            double cxi = UFFParameters[m_atom_types[i]][cx];
            double cxj = UFFParameters[m_atom_types[j]][cx];
            vdW["C_ij"] = sqrt(cDi * cDj) * m_vdw_force;
            vdW["i"] = i;
            vdW["j"] = j;
            vdW["r0_ij"] = sqrt(cxi * cxj);
            m_vdws.push_back(vdW);

            if (m_mol.m_partial_charges.size() <= i || m_mol.m_partial_charges.size() <= j)
                continue;
            if (std::abs(m_mol.m_partial_charges[i] < q_thresh) || std::abs(m_mol.m_partial_charges[j] < q_thresh))
                continue;
            json esp = EQJson;
            esp["i"] = i;
            esp["j"] = j;
            esp["q_i"] = m_mol.m_partial_charges[i];
            esp["q_j"] = m_mol.m_partial_charges[j];
            m_esps.push_back(esp);
        }
    }
    /*
    for(int i = 0; i < m_1_4_charges.size(); ++i)
    {
        for(const auto j : m_1_4_charges[i])
        {
            if(i == j)
                continue;
            json esp = EQJson;
            esp["i"] = i;
            esp["j"] = j;
            esp["q_i"] = charges[i];
            esp["q_j"] = charges[j];
            esp["epsilon"] = 0.85;
            m_esps.push_back(esp);
        }
    }
*/
}

void ForceFieldGenerator::AssignUffAtomTypes()
{
    std::vector<int> atom_types = m_atom_types;
    for (int i = 0; i < atom_types.size(); ++i) {
        switch (atom_types[i]) {
        case 1: // Hydrogen
            if (m_stored_bonds[i].size() == 2)
                m_atom_types[i] = 3; // Bridging Hydrogen
            else
                m_atom_types[i] = 1;
            break;
        case 2: // Helium
            m_atom_types[i] = 4;
            break;
        case 3: // Li
            m_atom_types[i] = 5;
            break;
        case 4: // Be
            m_atom_types[i] = 6;
            break;
        case 5: // B
            m_atom_types[i] = 7;
            break;
        case 6: // C
            if (m_coordination[i] == 4)
                m_atom_types[i] = 9;
            else if (m_coordination[i] == 3)
                m_atom_types[i] = 10;
            else // if (coordination == 2)
                m_atom_types[i] = 12;
            break;
        case 7: // N
            if (m_coordination[i] == 3)
                m_atom_types[i] = 13;
            else if (m_coordination[i] == 2)
                m_atom_types[i] = 14;
            else // if (coordination == 2)
                m_atom_types[i] = 15;
            break;
        case 8: // O
            if (m_coordination[i] == 3)
                m_atom_types[i] = 17;
            else if (m_coordination[i] == 2)
                m_atom_types[i] = 19;
            else // if (coordination == 2)
                m_atom_types[i] = 21;
            break;
        case 9: // F
            m_atom_types[i] = 22;
            break;
        case 10: // Ne
            m_atom_types[i] = 23;
            break;
        case 11: // Na
            m_atom_types[i] = 24;
            break;
        case 12: // Mg
            m_atom_types[i] = 25;
            break;
        case 13: // Al
            m_atom_types[i] = 26;
            break;
        case 14: // Si
            m_atom_types[i] = 27;
            break;
        case 15: // P
#pragma message("maybe add organometallic phosphorous (28)")
            m_atom_types[i] = 29;
            break;
        case 16: // S
            if (m_coordination[i] == 2)
                m_atom_types[i] = 31;
            else // ok, currently we do not discriminate between SO2 and SO3, just because there is H2SO3 and H2SO4
                m_atom_types[i] = 32;
#pragma message("we have to add organic S")
            break;
        case 17: // Cl
            m_atom_types[i] = 36;
            break;
        case 18: // Ar
            m_atom_types[i] = 37;
            break;
        case 19: // K
            m_atom_types[i] = 38;
            break;
        case 20: // Ca
            m_atom_types[i] = 39;
            break;
        case 21: // Sc
            m_atom_types[i] = 40;
            break;
        case 22: // Ti
            if (m_coordination[i] == 6)
                m_atom_types[i] = 41;
            else
                m_atom_types[i] = 42;
            break;
        case 23: // Va
            m_atom_types[i] = 43;
            break;
        case 24: // Cr
            m_atom_types[i] = 44;
            break;
        case 25: // Mn
            m_atom_types[i] = 45;
            break;
        case 26: // Fe
            if (m_coordination[i] == 6)
                m_atom_types[i] = 46;
            else
                m_atom_types[i] = 47;
            break;
        case 27: // Co
            m_atom_types[i] = 48;
            break;
        case 28: // Ni
            m_atom_types[i] = 49;
            break;
        case 29: // Cu
            m_atom_types[i] = 50;
            break;
        case 30: // Zn
            m_atom_types[i] = 51;
            break;
        case 31: // Ga
            m_atom_types[i] = 52;
            break;
        case 32: // Ge
            m_atom_types[i] = 53;
            break;
        case 33: // As
            m_atom_types[i] = 54;
            break;
        case 34: // Se
            m_atom_types[i] = 55;
            break;
        case 35: // Br
            m_atom_types[i] = 56;
            break;
        case 36: // Kr
            m_atom_types[i] = 57;
            break;
        case 37: // Rb
            m_atom_types[i] = 58;
            break;
        case 38: // Sr
            m_atom_types[i] = 59;
            break;
        case 39: // Y
            m_atom_types[i] = 60;
            break;
        case 40: // Zr
            m_atom_types[i] = 61;
            break;
        case 41: // Nb
            m_atom_types[i] = 62;
            break;
        case 42: // Mo
            if (m_coordination[i] == 6)
                m_atom_types[i] = 63;
            else
                m_atom_types[i] = 64;
            break;
        case 43: // Tc
            m_atom_types[i] = 65;
            break;
        case 44: // Ru
            m_atom_types[i] = 66;
            break;
        case 45: // Rh
            m_atom_types[i] = 67;
            break;
        case 46: // Pd
            m_atom_types[i] = 68;
            break;
        case 47: // Ag
            m_atom_types[i] = 69;
            break;
        case 48: // Cd
            m_atom_types[i] = 70;
            break;
        case 49: // In
            m_atom_types[i] = 71;
            break;
        case 50: // Sn
            m_atom_types[i] = 72;
            break;
        case 51: // Sb
            m_atom_types[i] = 73;
            break;
        case 52: // Te
            m_atom_types[i] = 74;
            break;
        case 53: // I
            m_atom_types[i] = 75;
            break;
        case 54: // Xe
            m_atom_types[i] = 76;
            break;
        default:
            m_atom_types[i] = 0;
        };
    }
}

// =============================================================================
// QMDFF term generation — Claude Generated (August 2026)
//
// Reference: xtb b7dbd36^:src/qmdff.f90 (module xtb_qmdff), deleted from xtb in
// commit b7dbd36 "Refactoring of external drivers (#568)".
//
// IMPORTANT SCOPE NOTE
// --------------------
// The genuine QMDFF parametrisation is a fit to a QM Hessian performed by
// Grimme's standalone QMDFF program, which writes a `solvent` parameter file;
// xtb's qmdff.f90 only ever *evaluated* such a file. Curcuma derives its own
// parameters from the UFF-style topology (and optionally refines the bond/angle
// force constants with `-qmdfffit`). The three lists built here therefore differ
// in origin from a real QMDFF parameter set:
//
//   torsions : an EXACT re-expression of curcuma's UFF torsion terms in the
//              QMDFF Fourier form, so that the QMDFF dissociation damping
//              applies. The barriers themselves are still UFF-derived.
//   NCI      : essentially parameter-free — C6 from the validated native D3
//              tables, R0 from r2r4, Z_A from valelff, alpha from the D3 R0ab
//              table. Only the atomic charges are external input.
//   HB/XB    : topological detection plus the reference's charge-dependent
//              scaling (hbpara). The QMDFF program derives the same list from
//              the QM topology.
//
// In other words: the FUNCTIONAL FORM is now the reference one, the PARAMETERS
// are curcuma's. That distinction matters when comparing against published
// QMDFF numbers — see docs/QMDFF.md.
// =============================================================================

void ForceFieldGenerator::setQMDFFTerms()
{
    m_qmdff_torsions.clear();
    m_qmdff_ncis.clear();
    m_qmdff_hbonds.clear();

    const int natoms = static_cast<int>(m_mol.m_atoms.size());
    if (natoms == 0)
        return;

    // -------------------------------------------------------------------------
    // 1) Torsions: rewrite E_UFF = V/2 (1 - cos(n phi0) cos(n phi)) as two QMDFF
    //    Fourier terms. With the QMDFF reference angle set to pi both erf
    //    branches coincide (dphi1 == dphi2 == phi - pi), so
    //      v1 (1 + cos(n phi)) + v2 (1 - cos(n phi))
    //        = (v1 + v2) + (v1 - v2) cos(n phi)
    //    reproduces the UFF term exactly for
    //      v1 = V/4 (1 - cos(n phi0)),  v2 = V/4 (1 + cos(n phi0)).
    //    The phases n*pi and n*pi + pi undo the phi -> phi - pi shift.
    // -------------------------------------------------------------------------
    for (const auto& d : m_dihedrals) {
        if (d.value("type", 1) == 0)
            continue;
        const double V = d.value("V", 0.0);
        const double n = d.value("n", 0.0);
        const double phi0 = d.value("phi0", 0.0);
        if (std::abs(V) < 1e-14 || n <= 0.0)
            continue;

        const double cn0 = std::cos(n * phi0);
        const double v1 = 0.25 * V * (1.0 - cn0);
        const double v2 = 0.25 * V * (1.0 + cn0);

        json terms = json::array();
        terms.push_back({ { "n", n }, { "phase", n * pi }, { "v", v1 } });
        terms.push_back({ { "n", n }, { "phase", n * pi + pi }, { "v", v2 } });

        json t;
        t["i"] = d["i"];
        t["j"] = d["j"];
        t["k"] = d["k"];
        t["l"] = d["l"];
        t["out_of_plane"] = false;
        t["phi0"] = pi;
        t["scale"] = 1.0;
        t["terms"] = terms;
        m_qmdff_torsions.push_back(t);
    }

    // -------------------------------------------------------------------------
    // 1b) Out-of-plane terms. QMDFF keeps them in the torsion list with
    //     tors(6,·) == 2 (qmdff.f90:466-504) instead of a separate improper term,
    //     so curcuma's UFF impropers are re-expressed in that form. Index order
    //     differs: curcuma's Inversion has the CENTRAL atom first, QMDFF has it
    //     second (its damping runs over the bonds j-i, j-k, j-l), hence the swap.
    //
    //     The single-minimum branch (rn > 0) is used:
    //       E = v (1 + cos(omega - omega0 + pi)) = v (1 - cos(omega - omega0))
    //     with omega0 taken from the reference structure, consistent with how the
    //     bond and angle references are taken. The barrier v is curcuma's UFF
    //     improper force constant — the reference derives it from a QM Hessian,
    //     so this magnitude is an approximation, not a ported value.
    // -------------------------------------------------------------------------
    for (const auto& inv : m_inversions) {
        if (inv.value("type", 1) == 0)
            continue;
        const double fc = inv.value("fc", 0.0);
        if (std::abs(fc) < 1e-14)
            continue;
        const int centre = inv["i"];
        const int s1 = inv["j"];
        const int s2 = inv["k"];
        const int s3 = inv["l"];
        if (centre >= natoms || s1 >= natoms || s2 >= natoms || s3 >= natoms)
            continue;

        Eigen::Matrix<double, 4, 3> dummy;
        const double omega0 = QMDFFTerms::outOfPlaneAngle(
            m_geometry.row(s1).transpose(), m_geometry.row(centre).transpose(),
            m_geometry.row(s2).transpose(), m_geometry.row(s3).transpose(),
            dummy, false);

        json terms = json::array();
        terms.push_back({ { "n", 1.0 }, { "phase", 0.0 }, { "v", 0.0 } });

        json t;
        t["i"] = s1;
        t["j"] = centre;
        t["k"] = s2;
        t["l"] = s3;
        t["out_of_plane"] = true;
        t["phi0"] = omega0;
        t["scale"] = fc;
        t["terms"] = terms;
        m_qmdff_torsions.push_back(t);
    }

    // -------------------------------------------------------------------------
    // 2) Topological distance classes for the non-covalent list.
    //    BFS over the bond graph; the QMDFF class nk is the bond count:
    //      1,2 -> 1   1,3 -> 2   1,4 -> 3   1,5 -> 4   further/disconnected -> 5
    //    Classes 1 and 2 have eps1 = eps2 = 0 in the reference and are omitted.
    // -------------------------------------------------------------------------
    constexpr int kUnreachable = 99;
    std::vector<std::vector<int>> bond_dist(natoms, std::vector<int>(natoms, kUnreachable));
    for (int start = 0; start < natoms; ++start) {
        bond_dist[start][start] = 0;
        std::vector<int> frontier{ start };
        int depth = 0;
        while (!frontier.empty() && depth < 5) {
            std::vector<int> next;
            ++depth;
            for (int a : frontier) {
                if (a < 0 || a >= static_cast<int>(m_stored_bonds.size()))
                    continue;
                for (int b : m_stored_bonds[a]) {
                    if (b < 0 || b >= natoms)
                        continue;
                    if (bond_dist[start][b] > depth) {
                        bond_dist[start][b] = depth;
                        next.push_back(b);
                    }
                }
            }
            frontier.swap(next);
        }
    }

    // -------------------------------------------------------------------------
    // 3) Non-covalent pairs. The functional form is ff_ini's (qmdff.f90:132-155):
    //      sr42 = 3 s8 r2r4_i r2r4_j = s8 * (C8/C6)
    //      R0   = a1 sqrt(3 r2r4_i r2r4_j) + a2 = a1 sqrt(C8/C6) + a2
    //
    //    Claude Generated (Aug 2026): the C6 source is **D4**, not the D3 of the 2014
    //    reference. D4 gives a charge- AND coordination-number-dependent C6, and we feed
    //    it the same QM charges the electrostatic term uses, so dispersion and
    //    electrostatics are consistent. Consequently the damping constants are D4's own
    //    (s6=1.0, s8=2.0, a1=0.58, a2=4.80 from the `d4param` registry) rather than
    //    QMDFF's D3-era 2.70/0.45/4.00 — mixing D4 C6 with D3-fitted damping would be
    //    worse than either. This changes every `-method qmdff` energy; see docs/QMDFF.md.
    // -------------------------------------------------------------------------
    const bool have_charges = (m_mol.m_partial_charges.size() == natoms);

    ConfigManager d4_config("d4param", m_parameter);
    const double d4_s6 = d4_config.get<double>("d4_s6", 1.0);
    const double d4_s8 = d4_config.get<double>("d4_s8", 2.0);
    const double d4_a1 = d4_config.get<double>("d4_a1", 0.58);
    const double d4_a2 = d4_config.get<double>("d4_a2", 4.80);

    D4ParameterGenerator d4(d4_config);
    d4.setBuildPairLists(false); // we own the pair loop; the JSON list would be discarded
    if (have_charges) {
        // Must precede GenerateParameters: reuse the QM charges instead of an internal
        // EEQ solve, so zeta(q) and the point-charge term see the same charges.
        d4.setTopologyCharges(m_partial_charges);
    }
    {
        const Matrix geometry_bohr = m_geometry * 1.889726125; // D4 expects Bohr
        d4.GenerateParameters(m_mol.m_atoms, geometry_bohr);
    }

    for (int i = 0; i < natoms; ++i) {
        for (int j = i + 1; j < natoms; ++j) {
            const int d = bond_dist[i][j];
            int nk = 5;
            if (d == 1)
                nk = 1;
            else if (d == 2)
                nk = 2;
            else if (d == 3)
                nk = 3;
            else if (d == 4)
                nk = 4;
            if (nk <= 2)
                continue; // eps1 = eps2 = 0 — the reference contributes nothing here

            const int zi = m_mol.m_atoms[i];
            const int zj = m_mol.m_atoms[j];

            const double c6 = d4_s6 * d4.getChargeWeightedC6(zi, zj, i, j);
            // C8/C6 = 3 * sqrtZr4r2_i * sqrtZr4r2_j (the D3/D4 convention)
            const double c8_over_c6 = 3.0 * d4.getSqrtZr4r2(zi) * d4.getSqrtZr4r2(zj);

            json n;
            n["i"] = i;
            n["j"] = j;
            n["nk"] = nk;
            n["c6"] = c6;
            n["r0_bj"] = d4_a1 * std::sqrt(c8_over_c6) + d4_a2;
            n["sr42"] = d4_s8 * c8_over_c6;
            n["zab"] = QMDFFData::valenceElectrons(zi) * QMDFFData::valenceElectrons(zj);
            n["alpha"] = QMDFFData::repulsionAlpha(zi, zj);
            n["qq"] = have_charges ? (m_mol.m_partial_charges[i] * m_mol.m_partial_charges[j]) : 0.0;
            m_qmdff_ncis.push_back(n);
        }
    }

    // -------------------------------------------------------------------------
    // 4) Hydrogen and halogen bonds. QMDFF stores triples (A, B, H) where A and B
    //    are the heavy partners; the strength is charge dependent (qmdff.f90:244-254):
    //      HB: c_A = hbpara(10, 5, q_A) * scalehb(Z_A), same for B
    //      XB: c   = hbpara(-6.5, 1, q_X) * scalexb(Z_X)
    //    Without partial charges the scaling degenerates, so the lists stay empty.
    // -------------------------------------------------------------------------
    if (have_charges) {
        auto isBonded = [&](int a, int b) {
            if (a < 0 || a >= static_cast<int>(m_stored_bonds.size()))
                return false;
            return std::find(m_stored_bonds[a].begin(), m_stored_bonds[a].end(), b) != m_stored_bonds[a].end();
        };

        // --- hydrogen bonds: H covalently bound to a donor A, acceptor B not bonded to H
        for (int h = 0; h < natoms; ++h) {
            if (m_mol.m_atoms[h] != 1)
                continue;
            for (int a = 0; a < natoms; ++a) {
                if (a == h || !isBonded(h, a))
                    continue;
                if (QMDFFData::scaleHB(m_mol.m_atoms[a]) <= 0.0)
                    continue;
                for (int b = 0; b < natoms; ++b) {
                    if (b == a || b == h)
                        continue;
                    if (QMDFFData::scaleHB(m_mol.m_atoms[b]) <= 0.0)
                        continue;
                    if (isBonded(h, b))
                        continue;
                    if (bond_dist[a][b] <= 2)
                        continue; // same functional group, not a hydrogen bridge

                    json e;
                    e["a"] = a;
                    e["b"] = b;
                    e["h"] = h;
                    e["c1"] = QMDFFTerms::hbpara(10.0, 5.0, m_mol.m_partial_charges[a])
                        * QMDFFData::scaleHB(m_mol.m_atoms[a]);
                    e["c2"] = QMDFFTerms::hbpara(10.0, 5.0, m_mol.m_partial_charges[b])
                        * QMDFFData::scaleHB(m_mol.m_atoms[b]);
                    e["halogen"] = false;
                    m_qmdff_hbonds.push_back(e);
                }
            }
        }

        // --- halogen bonds: X covalently bound to A, acceptor B elsewhere
        for (int x = 0; x < natoms; ++x) {
            const double sxb = QMDFFData::scaleXB(m_mol.m_atoms[x]);
            if (sxb <= 0.0)
                continue;
            for (int a = 0; a < natoms; ++a) {
                if (a == x || !isBonded(x, a))
                    continue;
                for (int b = 0; b < natoms; ++b) {
                    if (b == a || b == x)
                        continue;
                    if (QMDFFData::scaleHB(m_mol.m_atoms[b]) <= 0.0)
                        continue;
                    if (bond_dist[x][b] <= 2)
                        continue;

                    json e;
                    e["a"] = a;
                    e["b"] = b;
                    e["h"] = x;
                    e["c1"] = sxb * QMDFFTerms::hbpara(-6.5, 1.0, m_mol.m_partial_charges[x]);
                    e["c2"] = 0.0;
                    e["halogen"] = true;
                    m_qmdff_hbonds.push_back(e);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // 5) LJ -> Morse promotion (qmdff.f90:994-1019, subroutine setmorse):
    //    bonds whose well depth exceeds morsethr dissociate properly with a Morse
    //    potential; weaker ones stay Lennard-Jones.
    // -------------------------------------------------------------------------
    const double morsethr = m_parameter.value("qmdff_morsethr", 99.0);
    int morse_count = 0;
    for (auto& b : m_bonds) {
        const bool morse = (b.value("fc", 0.0) > morsethr);
        b["qmdff_potential"] = morse ? 1 : 0;
        if (morse)
            ++morse_count;
    }

    // Visible at verbosity >= 1: silently dropping two of the five non-bonded
    // contributions would otherwise look like a converged, complete calculation.
    if (!have_charges && CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::warn("QMDFF: no partial charges available — electrostatics and HB/XB "
                            "terms are switched off. Use -qmdfffit to derive charges from a "
                            "QM calculation.");
    }
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("qmdff_torsions", static_cast<int>(m_qmdff_torsions.size()));
        CurcumaLogger::param("qmdff_nci_pairs", static_cast<int>(m_qmdff_ncis.size()));
        CurcumaLogger::param("qmdff_hb_xb_terms", static_cast<int>(m_qmdff_hbonds.size()));
        CurcumaLogger::param("qmdff_morse_bonds", fmt::format("{} (threshold {:.3f})", morse_count, morsethr));
    }
}

json ForceFieldGenerator::getParameter()
{
    json parameters = m_parameter;
    parameters["bonds"] = Bonds();
    parameters["angles"] = Angles();
    parameters["dihedrals"] = Dihedrals();
    parameters["inversions"] = Inversions();
    parameters["vdws"] = vdWs();
    parameters["esps"] = ESPs();

    // Claude Generated (Aug 2026): QMDFF-specific lists (empty for UFF/D3)
    if (!m_qmdff_torsions.empty())
        parameters["qmdff_torsions"] = m_qmdff_torsions;
    if (!m_qmdff_ncis.empty())
        parameters["qmdff_ncis"] = m_qmdff_ncis;
    if (!m_qmdff_hbonds.empty())
        parameters["qmdff_hbonds"] = m_qmdff_hbonds;

    return parameters;
}

json ForceFieldGenerator::Bonds() const
{
    json bonds;
    int index = 0;

    for (int i = 0; i < m_bonds.size(); ++i) {
        if (m_bonds[i]["type"] != 0) {
            bonds[index] = m_bonds[i];
            index++;
        }
    }
    return bonds;
}

json ForceFieldGenerator::Angles() const
{
    json angles;
    int index = 0;

    for (int i = 0; i < m_angles.size(); ++i) {
        if (m_angles[i]["type"] != 0) {
            angles[index] = m_angles[i];
            index++;
        }
    }
    return angles;
}

json ForceFieldGenerator::Dihedrals() const
{
    json dihedrals;
    int index = 0;

    for (int i = 0; i < m_dihedrals.size(); ++i) {
        if (m_dihedrals[i]["type"] != 0) {
            dihedrals[index] = m_dihedrals[i];
            index++;
        }
    }
    return dihedrals;
}

json ForceFieldGenerator::Inversions() const
{
    json inversions;
    int index = 0;
    for (int i = 0; i < m_inversions.size(); ++i) {
        if (m_inversions[i]["type"] != 0 && m_inversions[i]["fc"] != 0) {
            inversions[index] = m_inversions[i];
            index++;
        }
    }
    return inversions;
}

json ForceFieldGenerator::vdWs() const
{
    json vdws;
    int index = 0;

    for (int i = 0; i < m_vdws.size(); ++i) {
        if (m_vdws[i]["type"] != 0) {
            vdws[index] = m_vdws[i];
            index++;
        }
    }
    return vdws;
}

json ForceFieldGenerator::ESPs() const
{
    json esps;
    int index = 0;

    for (int i = 0; i < m_esps.size(); ++i) {
        if (m_esps[i]["type"] != 0) {
            esps[index] = m_esps[i];
            index++;
        }
    }
    return esps;
}

// Claude Generated 2025: D3 parameter generation
void ForceFieldGenerator::GenerateD3Parameters()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating D3 dispersion parameters");
    }

    // Claude Generated (2025): Ensure d3_s9=1.0 for ATM generation
    double d3_s9 = m_parameter.value("d3_s9", 1.0);
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("d3_s9_in_config", d3_s9);
    }

    // Create D3 parameter generator with force field configuration + explicit s9=1.0
    json d3_config_json = m_parameter;
    d3_config_json["d3_s9"] = d3_s9;  // Claude Generated (2025): Force ATM generation

    ConfigManager d3_config("d3param", d3_config_json);
    D3ParameterGenerator d3_gen(d3_config);

    // Extract atom types from molecule
    std::vector<int> atoms = m_mol.m_atoms;

    // Generate D3 parameters
    d3_gen.GenerateParameters(atoms, m_geometry);
    json d3_params = d3_gen.getParameters();

    // Merge D3 parameters into main parameter set
    if (d3_params.contains("d3_dispersion_pairs")) {
        m_parameter["d3_dispersion_pairs"] = d3_params["d3_dispersion_pairs"];
        m_parameter["d3_damping"] = d3_params["d3_damping"];
        m_parameter["d3_enabled"] = true;

        // Claude Generated (2025): Extract ATM triples for three-body dispersion
        if (d3_params.contains("atm_triples") && !d3_params["atm_triples"].is_null()) {
            m_parameter["atm_triples"] = d3_params["atm_triples"];
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::param("atm_triples", static_cast<int>(d3_params["atm_triples"].size()));
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D3 parameter generation completed");
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("d3_time", fmt::format("{} ms", duration.count()));
            if (d3_params.contains("d3_dispersion_pairs")) {
                CurcumaLogger::param("d3_pairs", static_cast<int>(d3_params["d3_dispersion_pairs"].size()));
            }
        }
    }
}

// Claude Generated 2025: D4 parameter generation
void ForceFieldGenerator::GenerateD4Parameters()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Generating D4 dispersion parameters");
    }

    // Create D4 parameter generator with force field configuration
    ConfigManager d4_config("d4param", m_parameter);
    D4ParameterGenerator d4_gen(d4_config);

    // Extract atom types from molecule
    std::vector<int> atoms = m_mol.m_atoms;

    // Generate D4 parameters with EEQ charges from geometry (Dec 2025 - Phase 2)
    d4_gen.GenerateParameters(atoms, m_geometry);
    json d4_params = d4_gen.getParameters();

    // Merge D4 parameters into main parameter set
    if (d4_params.contains("d4_dispersion_pairs")) {
        m_parameter["d4_dispersion_pairs"] = d4_params["d4_dispersion_pairs"];
        m_parameter["d4_damping"] = d4_params["d4_damping"];
        m_parameter["d4_enabled"] = true;

        // Claude Generated (2025): Extract ATM triples for three-body dispersion
        if (d4_params.contains("atm_triples") && !d4_params["atm_triples"].is_null()) {
            m_parameter["atm_triples"] = d4_params["atm_triples"];
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::param("atm_triples", static_cast<int>(d4_params["atm_triples"].size()));
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D4 parameter generation completed");
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("d4_time", fmt::format("{} ms", duration.count()));
            if (d4_params.contains("d4_dispersion_pairs")) {
                CurcumaLogger::param("d4_pairs", static_cast<int>(d4_params["d4_dispersion_pairs"].size()));
            }
        }
    }
}

// ============================================================================
// Claude Generated (December 19, 2025): UFF-D3 Hybrid Method
// ============================================================================

json ForceFieldGenerator::GenerateUFFD3Parameters()
{
    /**
     * @brief Generate UFF bonded parameters + native D3 dispersion
     *
     * This method combines:
     * 1. UFF bonded terms (bonds, angles, dihedrals, inversions, vdW)
     * 2. Native D3 dispersion (validated 10/11 molecules <1% error)
     *
     * UFF-D3 provides a fast, accurate hybrid method for molecular mechanics
     * with geometry-dependent dispersion correction.
     *
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010) [D3-BJ]
     *
     * Claude Generated: December 19, 2025
     */

    auto start_time = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Generating UFF-D3 parameters");
        CurcumaLogger::param("atoms", m_atoms);
    }

    // Step 1: Generate UFF bonded parameters (bonds, angles, dihedrals, inversions, vdW)
    Generate();  // Fills m_bonds, m_angles, m_dihedrals, m_inversions, m_vdws

    // Step 2: Get UFF parameters as JSON
    json uff_params = getParameter();

    // Step 3: Add D3 dispersion correction
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Adding native D3 dispersion correction");
    }

    try {
        // D3 configuration with PBE0/BJ parameters (recommended for UFF-D3)
        json d3_config_json = {
            {"d3_a1", 0.4145},   // PBE0/BJ damping parameter 1
            {"d3_a2", 4.8593},   // PBE0/BJ damping parameter 2 (Bohr)
            {"d3_alp", 14.0},    // Alpha parameter
            {"d3_s6", 1.0},      // C6 scaling factor
            {"d3_s8", 1.2177},    // C8 scaling factor (PBE0)
            {"d3_s9", 1.0}       // ATM three-body scaling factor
        };

        ConfigManager d3_config("d3param", d3_config_json);
        D3ParameterGenerator d3_gen(d3_config);

        // Generate D3 parameters with geometry-dependent CN calculation
        d3_gen.GenerateParameters(m_mol.m_atoms, m_geometry);
        json d3_params = d3_gen.getParameters();

        // Merge D3 parameters into UFF parameter set
        if (d3_params.contains("d3_dispersion_pairs")) {
            uff_params["d3_dispersion_pairs"] = d3_params["d3_dispersion_pairs"];
            uff_params["d3_damping"] = d3_params["d3_damping"];
            uff_params["d3_enabled"] = true;

            // Claude Generated (2025): Add ATM triples for three-body dispersion
            if (d3_params.contains("atm_triples") && !d3_params["atm_triples"].is_null()) {
                uff_params["atm_triples"] = d3_params["atm_triples"];
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::param("atm_triples", static_cast<int>(d3_params["atm_triples"].size()));
                }
            }

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::param("d3_pairs", static_cast<int>(d3_params["d3_dispersion_pairs"].size()));
            }
        } else {
            if (CurcumaLogger::get_verbosity() >= 1) {
                CurcumaLogger::warn("D3 generated no dispersion pairs");
            }
            uff_params["d3_enabled"] = false;
        }

    } catch (const std::exception& e) {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::error(fmt::format("D3 parameter generation failed: {}", e.what()));
        }
        uff_params["d3_enabled"] = false;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::success("UFF-D3 parameter generation completed");
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("total_time", fmt::format("{} ms", duration.count()));
        }
    }

    return uff_params;
}

json ForceFieldGenerator::GenerateD3OnlyParameters(const std::string& preset)
{
    /**
     * @brief Generate pure D3 dispersion parameters without bonded terms
     *
     * Supports 7 functional presets:
     * - pbe0 (default), blyp, b3lyp, tpss, pbe, bp86, gfnff
     *
     * Returns JSON with ONLY:
     * - d3_dispersion_pairs
     * - d3_damping
     * - d3_enabled = true
     *
     * NO bonded terms (bonds, angles, dihedrals, inversions, vdw)
     *
     * Claude Generated: December 21, 2025
     */

    auto start_time = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Generating D3-only dispersion parameters");
        CurcumaLogger::param("preset", preset);
        CurcumaLogger::param("atoms", m_atoms);
    }

    json d3_params;

    try {
        // Create D3ParameterGenerator with preset
        D3ParameterGenerator d3_gen = D3ParameterGenerator::createForMethod(preset);

        // Generate D3 parameters
        d3_gen.GenerateParameters(m_mol.m_atoms, m_geometry);
        json generated = d3_gen.getParameters();

        // Extract ONLY D3 data
        if (generated.contains("d3_dispersion_pairs")) {
            d3_params["d3_dispersion_pairs"] = generated["d3_dispersion_pairs"];
            d3_params["d3_damping"] = generated["d3_damping"];
            d3_params["d3_enabled"] = true;

            // Claude Generated (2025): Extract ATM triples for three-body dispersion
            if (generated.contains("atm_triples") && !generated["atm_triples"].is_null()) {
                d3_params["atm_triples"] = generated["atm_triples"];
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::param("atm_triples", static_cast<int>(generated["atm_triples"].size()));
                }
            }

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::param("d3_pairs", static_cast<int>(generated["d3_dispersion_pairs"].size()));
            }
        } else {
            CurcumaLogger::warn("D3 generated no dispersion pairs");
            d3_params["d3_enabled"] = false;
        }

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("D3 generation failed: {}", e.what()));
        d3_params["d3_enabled"] = false;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::success("D3-only parameter generation completed");
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("total_time", fmt::format("{} ms", duration.count()));
        }
    }

    return d3_params;
}
