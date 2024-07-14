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

#include "src/core/global.h"
#include "src/core/topology.h"
#include "src/tools/general.h"

#include "json.hpp"
using json = nlohmann::json;

#include "forcefieldgenerator.h"

ForceFieldGenerator::ForceFieldGenerator(const json& controller)
{
    m_parameter = MergeJson(FFGenerator, controller);
    m_method = m_parameter["method"];
}

void ForceFieldGenerator::setMolecule(const Molecule& molecule)
{
    m_geometry = molecule.getGeometry();
    m_molecule = molecule;
}

void ForceFieldGenerator::Generate(const std::vector<std::pair<int, int>>& formed_bonds)
{
    m_atom_types = std::vector<int>(m_molecule.Atoms().size(), 0);
    m_coordination = std::vector<int>(m_molecule.Atoms().size(), 0);
    m_topo = Eigen::MatrixXd::Zero(m_molecule.Atoms().size(), m_molecule.Atoms().size());
    TContainer bonds;
    m_scaling = 1.4;
    if (formed_bonds.size() == 0) {
        for (int i = 0; i < m_molecule.Atoms().size(); ++i) {
            m_stored_bonds.push_back(std::vector<int>());
            m_ignored_vdw.push_back(std::set<int>({ i }));
            for (int j = 0; j < m_molecule.Atoms().size() && m_stored_bonds[i].size() < CoordinationNumber[m_molecule.Atoms()[i]]; ++j) {
                if (i == j)
                    continue;
                double x_i = m_geometry(i, 0) * m_au;
                double x_j = m_geometry(j, 0) * m_au;

                double y_i = m_geometry(i, 1) * m_au;
                double y_j = m_geometry(j, 1) * m_au;

                double z_i = m_geometry(i, 2) * m_au;
                double z_j = m_geometry(j, 2) * m_au;

                double r_ij = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));

                if (r_ij <= (Elements::CovalentRadius[m_molecule.Atoms()[i]] + Elements::CovalentRadius[m_molecule.Atoms()[j]]) * m_scaling * m_au) {
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
    AssignUffAtomTypes();
    m_glob_scale = 7.25;
    // if (m_rings)
    // std::cout << "Crude ring finding method ... " << std::endl;
    int maxsize = m_atom_types.size();
    int maxcache = m_atom_types.size();
    if (m_atom_types.size() > 1000)
        maxcache = 5;
    m_identified_rings = Topology::FindRings(m_stored_bonds, m_atom_types.size(), maxsize, maxcache);
    // std::cout << "... done!" << std::endl;

    if (m_method.compare("uff-d3") == 0) {
        m_parameter["d3"] = 1;
        m_parameter["vdw_scaling"] = 0;
    } else if (m_method.compare("qmdff") == 0) {
        m_ff_type = 2;
        m_parameter["d3"] = 1;
        m_parameter["vdw_scaling"] = 0;
    }

    m_uff_bond_force = m_parameter["bond_force"];
    m_uff_angle_force = m_parameter["angle_force"];
    m_uff_dihedral_force = m_parameter["torsion_force"];
    m_uff_inversion_force = m_parameter["inversion_force"];
    m_vdw_force = m_parameter["vdw_force"];

    // std::cout << m_parameter << std::endl;

    setBonds(bonds);

    setAngles();

    setDihedrals();

    setInversions();

    setvdWs();
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
    for (const auto& bond : bonds.Storage()) {
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
        double r0_ij = UFFBondRestLength(bond[0], bond[1], bond_order);
        uffbond["r0_ij"] = r0_ij;
        double cZi = UFFParameters[m_atom_types[bond[0]]][cZ];
        double cZj = UFFParameters[m_atom_types[bond[1]]][cZ];
        uffbond["fc"] = m_uff_bond_force * cZi * cZj / (r0_ij * r0_ij * r0_ij);
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
            inversion["i"] = i;
            inversion["j"] = m_stored_bonds[i][0];
            inversion["k"] = m_stored_bonds[i][1];
            inversion["l"] = m_stored_bonds[i][2];
            if (std::find(m_inversions.begin(), m_inversions.end(), inversion) == m_inversions.end()) {
                m_inversions.push_back(inversion);
            }
        }
    }
}

void ForceFieldGenerator::setAngles()
{
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
        if (6 <= m_molecule.Atoms()[i] && m_molecule.Atoms()[i] <= 8) {
            C0 = 1.0;
            C1 = -1.0;
            C2 = 0.0;
            kijkl = 6 * m_uff_inversion_force;
            if (m_molecule.Atoms()[j] == 8 || m_molecule.Atoms()[k] == 8 || m_molecule.Atoms()[l] == 8)
                kijkl = 50 * m_uff_inversion_force;
        } else {
            double w0 = pi / 180.0;
            switch (m_molecule.Atoms()[i]) {
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

void ForceFieldGenerator::setvdWs()
{
    double vdw_scale = 7.5;
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
        }
    }
}

void ForceFieldGenerator::AssignUffAtomTypes()
{
    std::vector<int> atom_types = m_molecule.Atoms();
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
        /* if (m_verbose) {
             std::cout << i << " " << m_atom_types[i] << " " << m_stored_bonds[i].size() << " " << m_atom_types[i] << std::endl;
         }*/
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
