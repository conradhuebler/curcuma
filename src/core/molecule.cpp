/*
 * <Some globale definition for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "elements.h"

#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include <Eigen/Dense>

#include <fmt/core.h>
#include <fmt/format.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <map>
#include <sstream>

#include "molecule.h"

Molecule::Molecule(int n, int q)
{
    m_charge = q;
    InitialiseEmptyGeometry(n);
}

Molecule::Molecule() = default;

Molecule::Molecule(const Molecule& other)
{
    m_geometry = other.m_geometry;
    m_charge = other.m_charge;
    m_fragments = other.m_fragments;
    m_atoms = other.m_atoms;
    m_name = other.m_name;
    m_energy = other.m_energy;
    m_spin = other.m_spin;
    m_bonds = other.m_bonds;
}
/*
Molecule& Molecule::operator=(const Molecule& other)
{
    m_geometry = other.m_geometry;
    m_charge = other.m_charge;
    m_fragments = other.m_fragments;
    m_atoms = other.m_atoms;
    m_name = other.m_name;
    m_energy = other.m_energy;
    m_spin = other.m_spin;
    return *this;
}
*/
Molecule::Molecule(const Molecule* other)
{
    m_geometry = other->m_geometry;
    m_charge = other->m_charge;
    m_fragments = other->m_fragments;
    m_atoms = other->m_atoms;
    m_name = other->m_name;
    m_energy = other->m_energy;
    m_spin = other->m_spin;
    m_bonds = other->m_bonds;
}
/*
Molecule& Molecule::operator=(const Molecule* other)
{
    m_geometry = other->m_geometry;
    m_charge = other->m_charge;
    m_fragments = other->m_fragments;
    m_atoms = other->m_atoms;
    m_name = other->m_name;
    m_energy = other->m_energy;
    m_spin = other->m_spin;
    return *this;
}
*/
Molecule::Molecule(const Mol& other)
{
    setXYZComment(other.m_commentline);
    m_geometry = other.m_geometry;
    m_charge = other.m_charge;
    m_atoms = other.m_atoms;
    m_energy = other.m_energy;
    m_spin = other.m_spin;
    m_bonds = other.m_bonds;
}

Molecule::Molecule(const Mol* other)
{
    setXYZComment(other->m_commentline);
    m_geometry = other->m_geometry;
    m_charge = other->m_charge;
    m_atoms = other->m_atoms;
    m_energy = other->m_energy;
    m_spin = other->m_spin;
    m_bonds = other->m_bonds;
}

Molecule::Molecule(const std::string& file)
{
    /*
    std::vector<std::string> lines;
    int atoms = 0;
    int index = 0;
    int i = 0;
    bool xyzfile = std::string(file).find(".xyz") != std::string::npos || std::string(file).find(".trj") != std::string::npos;
    if (xyzfile) {
        auto m_file = new std::ifstream(file);
        for (std::string line; getline(*m_file, line);) {
            if (index == 0 && xyzfile) {
                try {
                    atoms = stoi(line);
                } catch (const std::invalid_argument& arg) {
                    atoms = 0;
                }
                InitialiseEmptyGeometry(atoms);
            }
            if (xyzfile) {
                if (i == 1)
                    setXYZComment(line);
                if (i > 1) {
                    setXYZ(line, i - 2);
                }
                if (i - 1 == atoms) {
                    break;
                }
                ++i;
            } else {
                setAtom(line, i);
            }
            index++;
        }
    }*/
    LoadMolecule(Files::LoadFile(file));
}

Molecule::~Molecule()
{
}

json Molecule::ExportJson() const
{
    json structure;
    structure["atoms"] = m_atoms.size();
    structure["elements"] = Tools::Vector2String(m_atoms);
    structure["name"] = m_name;
    for (int i = 0; i < m_atoms.size(); ++i) {
        structure["atom" + std::to_string(i)] = Tools::DoubleVector2String({ m_geometry[i][0], m_geometry[i][1], m_geometry[i][2] });
    }
    structure["charge"] = m_charge;
    return structure;
}

void Molecule::WriteJsonFile(const std::string& filename)
{
    std::ofstream input;
    input.open(filename, std::ios::out);
    input << ExportJson();
    input.close();
}

void Molecule::ImportJson(const std::string& jsonfile)
{
    json molecule;
    std::ifstream file(jsonfile);
    file >> molecule;
    ImportJson(molecule);
}

void Molecule::ImportJson(const json& molecule)
{
    int atoms = molecule["atoms"];
    m_name = molecule["name"];
    m_atoms = Tools::String2Vector(molecule["elements"]);
    m_charge = molecule["charge"];
    for (int i = 0; i < atoms; ++i) {
        auto position = Tools::String2DoubleVec(molecule["atom" + std::to_string(i)], "|");
        m_geometry.push_back(std::array<double, 3>({ position[0], position[1], position[2] }));
    }
}
void Molecule::Initialise(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    m_charge = charge;
    m_spin = spin;
    Molecule mol;
    for (int i = 0; i < natoms; ++i) {
        m_geometry.push_back({ coord[3 * i], coord[3 * i + 1], coord[3 * i + 2] });
        m_atoms.push_back(attyp[i]);
        m_mass += Elements::AtomicMass[attyp[i]];
    }
}

void Molecule::ApplyReorderRule(const std::vector<int>& rule)
{
    Molecule mol;
    for (auto& i : rule)
        mol.addPair(Atom(i));
    m_geometry = mol.m_geometry;
    m_atoms = mol.m_atoms;
}

void Molecule::print_geom(bool moreinfo) const
{
    std::cout << AtomCount() << std::endl;
    std::cout << Name() << " " << std::setprecision(12) << Energy() << std::endl;
    for (int i = 0; i < AtomCount(); i++) {
        printf("%s %8.5f %8.5f %8.5f\n", Elements::ElementAbbr[m_atoms[i]].c_str(), m_geometry[i][0], m_geometry[i][1], m_geometry[i][2]);
    }
    if (moreinfo) {
        std::cout << std::endl
                  << std::endl
                  << "***********************************************************" << std::endl;
        std::cout << "**         Center = " << Centroid().transpose() << std::endl;
        std::cout << "**         Number of Fragments = " << GetFragments().size() << std::endl;
        std::cout << "**         Ia = " << Ia() << std::endl;
        std::cout << "**         Ib = " << Ib() << std::endl;
        std::cout << "**         Ic = " << Ic() << std::endl;
        std::cout << "***********************************************************" << std::endl;
    }
}

int Molecule::Check() const
{
    return 0;

    for (int i = 0; i < AtomCount(); ++i) {
        if (std::isnan(m_geometry[i][0]) || std::isnan(m_geometry[i][1]) || std::isnan(m_geometry[i][2]))
            return 2;
        for (int j = 0; j < i; ++j) {
            if (CalculateDistance(i, j) < 1e-1)
                return 1;
        }
    }
    return 0;
}

void Molecule::printFragmente()
{
    if (m_fragments.size() == 0)
        GetFragments();

    std::cout << std::endl
              << std::endl
              << "***********************************************************" << std::endl;
    std::cout << "**         Center = " << Centroid().transpose() << std::endl;
    std::cout << "**         Number of Fragments = " << GetFragments().size() << std::endl;
    std::cout << "**         Ia = " << Ia() << std::endl;
    std::cout << "**         Ib = " << Ib() << std::endl;
    std::cout << "**         Ic = " << Ic() << std::endl;
    std::cout << "***********************************************************" << std::endl;

    for (std::size_t i = 0; i < m_fragments.size(); ++i) {
        for (const auto& atom : m_fragments[i]) {
            printf("%s(%i) %8.5f %8.5f %8.5f\n", Elements::ElementAbbr[m_atoms[atom]].c_str(), i + 1, m_geometry[atom][0], m_geometry[atom][1], m_geometry[atom][2]);
        }
    }
}

void Molecule::printAtom(int i) const
{
    if (i < AtomCount())
        printf("%s %8.5f %8.5f %8.5f", Elements::ElementAbbr[m_atoms[i]].c_str(), m_geometry[i][0], m_geometry[i][1], m_geometry[i][2]);
}


void Molecule::InitialiseEmptyGeometry(int atoms)
{
    for (int i = 0; i < atoms; i++) {
        std::array<double, 3> atom = { 0, 0, 0 };
        m_geometry.push_back(atom);
    }
    m_dirty = true;
}

bool Molecule::addPair(const std::pair<int, Position>& atom)
{
    bool exist = true;
    // const std::array<double, 3> at = { atom.second(0),  atom.second(1),  atom.second(2)};
    m_geometry.push_back({ atom.second(0), atom.second(1), atom.second(2) });
    m_atoms.push_back(atom.first);
    m_mass += Elements::AtomicMass[atom.first];

    for (std::size_t i = 0; i < AtomCount(); ++i)
        for (std::size_t j = i + 1; j < AtomCount(); ++j)
            if (CalculateDistance(i, j) < 1e-6)
                exist = false;

    m_dirty = true;

    return exist;
}

double Molecule::CalculateMass()
{
    double mass = 0;
    for (int atom : m_atoms)
        mass += Elements::AtomicMass[atom];
    m_mass = mass;
    return mass;
}

void Molecule::setXYZComment(const std::string& comment)
{
    StringList list = Tools::SplitString(comment);
    if (comment.find("Curcuma") != std::string::npos && list.size() >= 8) {
        try {
            setEnergy(std::stod(list[4]));
            setCharge(std::stod(list[9]));
        } catch (const std::invalid_argument& what) {
            try {
                setEnergy(std::stod(list[3]));
                setCharge(std::stod(list[8]));
            } catch (const std::invalid_argument& what) {
            }
        }
    } else {
        if (list.size() == 7) {
            setXYZComment_7(list);
        } else if (list.size() == 6) {
            setXYZComment_6(list);
        } else if (list.size() == 4) {
            setXYZComment_4(list);
        } else if (list.size() == 1) {
            if (list[0].compare("") == 0) {
                // Ignore empty comment line
            } else {
                try {
                    setEnergy(std::stod(list[0]));
                } catch (const std::string& what_arg) {
                } catch (const std::invalid_argument& arg) {
                }
            }
        } else {
            for (const std::string& s : list) {
                double energy = 0;
                if (Tools::isDouble(s)) {
                    energy = std::stod(s);
                    setEnergy(energy);
                    break;
                }
            }
        }
    }
}

bool Molecule::setXYZComment_0(const StringList& list)
{
    for (const std::string& s : list) {
        double energy = 0;
        if (Tools::isDouble(s)) {
            energy = std::stod(s);
            setEnergy(energy);
            break;
        }
    }

    return true;
}

bool Molecule::setXYZComment_1(const StringList& list)
{
    for (const std::string& s : list) {
        double energy = 0;
        if (Tools::isDouble(s)) {
            energy = std::stod(s);
            setEnergy(energy);
            break;
        }
    }
    return true;
}

bool Molecule::setXYZComment_2(const StringList& list)
{
    for (const std::string& s : list) {
        double energy = 0;
        if (Tools::isDouble(s)) {
            energy = std::stod(s);
            setEnergy(energy);
            break;
        }
    }
    return true;
}

bool Molecule::setXYZComment_3(const StringList& list)
{
    for (const std::string& s : list) {
        double energy = 0;
        if (Tools::isDouble(s)) {
            energy = std::stod(s);
            setEnergy(energy);
            break;
        }
    }
    return true;
}

bool Molecule::setXYZComment_4(const StringList& list)
{
    if (list[0].compare("SCF") == 0 && list[1].compare("done") == 0) {
        try {
            setEnergy(std::stod((list[2])));
        } catch (const std::string& what_arg) {
            setEnergy(0);
        } catch (const std::invalid_argument& argument) {
            setEnergy(0);
        }
    } else {
        setName(list[0]);
        if (list[3] == "")
            setEnergy(0);
        else {
            try {
                setEnergy(std::stod((list[3])));
            } catch (const std::string& what_arg) {
                setEnergy(0);
            } catch (const std::invalid_argument& what) {
                setEnergy(0);
            }
        }
    }
    return true;
}

bool Molecule::setXYZComment_5(const StringList& list)
{
    return true;
}

bool Molecule::setXYZComment_6(const StringList& list)
{
    try {
        setEnergy(std::stod((list[3])));
    } catch (const std::string& what_arg) {
        setEnergy(0);
    } catch (const std::invalid_argument& ia) {
        setEnergy(0);
    }

    return true;
}

bool Molecule::setXYZComment_7(const StringList& list)
{
    if (list[2].compare("gnorm:") == 0) {
        try {
            setEnergy(std::stod((list[1])));
        } catch (const std::string& what_arg) {
            setEnergy(0);
        } catch (const std::invalid_argument& argument) {
            setEnergy(0);
        }
    } else {
        try {
            setEnergy(std::stod((list[4])));
        } catch (const std::string& what_arg) {
            setEnergy(0);
        } catch (const std::invalid_argument& argument) {
            setEnergy(0);
        }
    }
    return true;
}

bool Molecule::setXYZComment_8(const StringList& list)
{
    return true;
}

bool Molecule::Contains(const std::pair<int, Position>& atom)
{
    for (std::size_t i = 0; i < AtomCount(); ++i) {
        if (GeometryTools::Distance(Atom(i).second, atom.second) < 1e-6)
            return true;
    }

    return false;
}

double Molecule::CalculateDistance(int i, int j) const
{
    if (i >= AtomCount() || j >= AtomCount())
        return 0;

    double x_i = m_geometry[i][0];
    double x_j = m_geometry[j][0];

    double y_i = m_geometry[i][1];
    double y_j = m_geometry[j][1];

    double z_i = m_geometry[i][2];
    double z_j = m_geometry[j][2];

    return sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));
}

double Molecule::DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const
{
    return pos1[0]*pos2[0]+pos1[1]*pos2[1]+pos1[2]*pos2[2];
}

double Molecule::CalculateAngle(int atom1, int atom2, int atom3) const
{
    std::array<double, 3> atom_0 = { m_geometry[atom2] }; // Proton
    std::array<double, 3> atom_1 = { m_geometry[atom3] }; // Acceptor
    std::array<double, 3> atom_2 = { m_geometry[atom1] }; // Donor

    std::array<double, 3> vec_1 = { atom_0[0]-atom_1[0], atom_0[1]-atom_1[1], atom_0[2]-atom_1[2] };
    std::array<double, 3> vec_2 = { atom_0[0]-atom_2[0], atom_0[1]-atom_2[1], atom_0[2]-atom_2[2] };
    
    return acos(DotProduct(vec_1,vec_2)/(sqrt(DotProduct(vec_1,vec_1)*DotProduct(vec_2,vec_2))))*360/2.0/pi;
}

void Molecule::ParseString(const std::string& internal, std::vector<std::string>& elements)
{
    std::string element;
    const char *delim = " ";
    for (const char & c : internal)
    {
        if(*delim != c)
            element.push_back(c);
        else
        {
            if (element.size()) {
                elements.push_back(element);
            }
            element.clear();
        }
    }

    elements.push_back(element);
}

void Molecule::setAtom(const std::string& internal, int i)
{
    std::vector<std::string> elements;
    ParseString(internal, elements);

    int atom = 0;
    if (Tools::isInt(elements[0]))
        atom = std::stoi(elements[0]);
    else
        atom = Elements::String2Element(elements[0]);
    m_atoms.push_back(atom);

    if(elements.size() == 7)
    {
        int atom_1 = stoi(elements[1]);
        int atom_2 = stoi(elements[2]);
        int atom_3 = stoi(elements[3]);

        double r_ij = stod(elements[4]);
        double omega = stod(elements[5]);
        double theta = stod(elements[6]);

        double r_ik = CalculateDistance(atom_2, atom_1);

        double x = r_ik + r_ij*cos(omega);
        double y = r_ij*sin(omega);
        double z = 0;
        m_geometry[i][0] = x;
        m_geometry[i][1] = y;
        m_geometry[i][2] = z;
    }
    m_dirty = true;
}

void Molecule::setXYZ(const std::string& internal, int i)
{
    std::vector<std::string > elements;
    ParseString(internal, elements);
    int atom = 0;
    if (Tools::isInt(elements[0]))
        atom = std::stoi(elements[0]);
    else
        atom = Elements::String2Element(elements[0]);
    m_atoms.push_back(atom);

    if (elements.size() >= 4) {
        double x = stod(elements[1]);
        double y = stod(elements[2]);
        double z = stod(elements[3]);

        m_geometry[i][0] = x;
        m_geometry[i][1] = y;
        m_geometry[i][2] = z;
    }

    m_dirty = true;
}

void Molecule::clear()
{
    m_atoms.clear();
    m_geometry.clear();
    m_dirty = true;
}

void Molecule::LoadMolecule(const Molecule& molecule)
{
    clear();
    m_charge = molecule.Charge();
    m_atoms = molecule.Atoms();
    m_energy = molecule.Energy();
    m_bonds = molecule.m_bonds;
    InitialiseEmptyGeometry(molecule.AtomCount());
    setGeometry(molecule.getGeometry());
}

void Molecule::LoadMolecule(const Molecule* molecule)
{
    clear();
    m_charge = molecule->Charge();
    m_atoms = molecule->Atoms();
    m_bonds = molecule->m_bonds;

    InitialiseEmptyGeometry(molecule->AtomCount());
    setGeometry(molecule->getGeometry());
}

void Molecule::LoadMolecule(const Mol& molecule)
{
    clear();

    setXYZComment(molecule.m_commentline);
    m_charge = molecule.m_charge;
    m_atoms = molecule.m_atoms;
    m_energy = molecule.m_energy;
    InitialiseEmptyGeometry(AtomCount());
    m_geometry = (molecule.m_geometry);
    m_bonds = molecule.m_bonds;
}

void Molecule::LoadMolecule(const Mol* molecule)
{
    clear();

    setXYZComment(molecule->m_commentline);
    m_charge = molecule->m_charge;
    m_atoms = molecule->m_atoms;
    m_energy = molecule->m_energy;
    InitialiseEmptyGeometry(AtomCount());
    m_geometry = (molecule->m_geometry);
    m_bonds = molecule->m_bonds;
}

Molecule Molecule::getFragmentMolecule(int fragment) const
{
    // Lets make that one day faster, but not today ...
    Molecule result;
    auto atoms = GetFragments()[fragment];
    for (auto atom : atoms) {
        result.addPair(Atom(atom));
    }
    return result;
}
Geometry Molecule::getGeometry(bool protons) const
{
    if (protons) {
        Geometry geometry(m_geometry.size(), 3);
        for (int i = 0; i < m_geometry.size(); ++i) {
            geometry(i, 0) = m_geometry[i][0];
            geometry(i, 1) = m_geometry[i][1];
            geometry(i, 2) = m_geometry[i][2];
        }
        return geometry;
    } else {
        std::vector<int> indicies;
        for (int i = 0; i < m_geometry.size(); ++i) {
            if (m_atoms[i] != 1) {
                indicies.push_back(i);
            }
        }
        Geometry geometry(indicies.size(), 3);
        int index = 0;
        for (int i : indicies) {
            geometry(index, 0) = m_geometry[i][0];
            geometry(index, 1) = m_geometry[i][1];
            geometry(index, 2) = m_geometry[i][2];
            index++;
        }
        return geometry;
    }
}

Geometry Molecule::getGeometry(const IntPair& pair, bool protons) const
{
    int start = pair.first;
    int end = pair.second;

    if (start < 0 || start >= m_geometry.size())
        start = 0;

    if (end < 0 || end >= m_geometry.size())
        end = m_geometry.size();

    Geometry geometry(m_geometry.size(), 3);
    int index = 0;

    if (protons) {
        for (int i = start; i < end; ++i) {
            geometry(index, 0) = m_geometry[i][0];
            geometry(index, 1) = m_geometry[i][1];
            geometry(index, 2) = m_geometry[i][2];
            index++;
        }
    } else {
        for (int i = start; i < end; ++i) {
            if (m_atoms[i] != 1) {
                geometry(index, 0) = m_geometry[i][0];
                geometry(index, 1) = m_geometry[i][1];
                geometry(index, 2) = m_geometry[i][2];
                index++;
            }
        }
    }
    return geometry.block(0, 0, index, 3);
}

Geometry Molecule::getGeometry(std::vector<int> atoms, bool protons) const
{
    Geometry geometry(m_geometry.size(), 3);
    int index = 0;
    if (protons) {
        for (int i : atoms) {
            geometry(index, 0) = m_geometry[i][0];
            geometry(index, 1) = m_geometry[i][1];
            geometry(index, 2) = m_geometry[i][2];
            index++;
        }
    } else {
        for (int i : atoms) {
            if (m_atoms[i] != 1) {
                geometry(index, 0) = m_geometry[i][0];
                geometry(index, 1) = m_geometry[i][1];
                geometry(index, 2) = m_geometry[i][2];
                index++;
            }
        }
    }
    return geometry.block(0, 0, index, 3);
}

Geometry Molecule::getGeometryByFragment(int fragment, bool protons) const
{
    if (fragment == -1)
        return getGeometry(protons);
    else
        return getGeometry(m_fragments[fragment], protons);
}

bool Molecule::setGeometry(const Geometry &geometry)
{
    if (geometry.rows() != m_geometry.size())
        return false;

    m_dirty = true;

    for (int i = 0; i < m_geometry.size(); ++i) {
        m_geometry[i][0] = geometry(i, 0);
        m_geometry[i][1] = geometry(i, 1);
        m_geometry[i][2] = geometry(i, 2);
    }

    return true;
}

bool Molecule::setGeometryByFragment(const Geometry& geometry, int fragment, bool protons)
{
    if (fragment >= GetFragments().size())
        return false;

    std::vector<int> frag = m_fragments[fragment];
    int index = 0;
    if (protons) {
        for (int i : frag) {
            m_geometry[i][0] = geometry(index, 0);
            m_geometry[i][1] = geometry(index, 1);
            m_geometry[i][2] = geometry(index, 2);
            index++;
        }
    } else {
        for (int i : frag) {
            if (Atom(i).first == 1)
                continue;
            m_geometry[i][0] = geometry(index, 0);
            m_geometry[i][1] = geometry(index, 1);
            m_geometry[i][2] = geometry(index, 2);
            index++;
        }
    }
    return true;
}

std::string Molecule::LowerDistanceMatrix() const
{
    std::ostringstream stream;
    for (int i = 0; i < AtomCount(); ++i) {
        for (int j = 0; j <= i; ++j) {
            stream << std::to_string(CalculateDistance(i, j)) + ",";
        }
        stream << std::endl;
    }
    std::string matrix;
    matrix = stream.str();
    return matrix;
}

std::vector<float> Molecule::LowerDistanceVector() const
{
    std::vector<float> vector;
    for (int i = 0; i < AtomCount(); ++i) {
        for (int j = 0; j < i; ++j) {
            vector.push_back((CalculateDistance(i, j)));
        }
    }
    return vector;
}

std::vector<double> Molecule::DeltaEN() const
{
    std::vector<double> vector;
    for (int i = 0; i < AtomCount(); ++i) {
        for (int j = 0; j < i; ++j) {
            vector.push_back(std::abs(Elements::PaulingEN[m_atoms[i]] - Elements::PaulingEN[m_atoms[j]]));
        }
    }
    return vector;
}

Position Molecule::Centroid(bool protons, int fragment) const
{
    return GeometryTools::Centroid(getGeometryByFragment(fragment, protons));
}

Position Molecule::MassCentroid(bool protons, int fragment) const
{
    Position pos = { 0, 0, 0 };
    double mass = 0.0;
    for (int i = 0; i < m_atoms.size(); ++i) {
        pos += Atom(i).second;
        mass += Elements::AtomicMass[Atom(i).first];
    }
    pos /= mass;
    return pos;
}

Eigen::Vector3d Molecule::COM(bool protons, int fragment)
{
    // todo implement heavy and fragements ...
    if (m_mass < 1)
        CalculateMass();
    Eigen::Vector3d com = { 0, 0, 0 };
    for (int i = 0; i < m_geometry.size(); ++i) {
        double mass = Elements::AtomicMass[m_atoms[i]];
        com(0) += mass * m_geometry[i][0];
        com(1) += mass * m_geometry[i][1];
        com(2) += mass * m_geometry[i][2];
    }
    com(0) /= m_mass;
    com(1) /= m_mass;
    com(2) /= m_mass;
    return com;
}

std::vector<Position> Molecule::CalculateDipoleMoments(const std::vector<double>& scaling) const
{
    std::vector<Position> dipole_moments;

    if (m_charges.size() != m_geometry.size()) {
        std::cout << "No partial charges available" << std::endl;
        return dipole_moments;
    }

    for (int f = 0; f < GetFragments().size(); ++f) {
        Position pos = { 0, 0, 0 }, dipole = { 0, 0, 0 };
        double mass = 0;
        for (int i : m_fragments[f]) {
            double m = Elements::AtomicMass[m_atoms[i]];
            mass += m;
            pos(0) += m * m_geometry[i][0];
            pos(1) += m * m_geometry[i][1];
            pos(2) += m * m_geometry[i][2];
        }
        pos(0) /= mass;
        pos(1) /= mass;
        pos(2) /= mass;
        for (int i : m_fragments[f]) {
            double scale = 3;
            if (scaling.size() > i)
                scale = scaling[i];
            dipole(0) += m_charges[i] * (m_geometry[i][0] - pos(0)) * scale;
            dipole(1) += m_charges[i] * (m_geometry[i][1] - pos(1)) * scale;
            dipole(2) += m_charges[i] * (m_geometry[i][2] - pos(2)) * scale;
            //  std::cout << scale << " ";
        }
        // std::cout << std::endl;
        dipole_moments.push_back(dipole);
    }
    return dipole_moments;
}

Position Molecule::CalculateDipoleMoment(const std::vector<double>& scaling) const
{
    double mass = 0;

    Position pos = { 0, 0, 0 }, dipole = { 0, 0, 0 };
    if (m_charges.size() != m_geometry.size()) {
        std::cout << "No partial charges available" << std::endl;
        return dipole;
    }
    for (int i = 0; i < m_geometry.size(); ++i) {
        double m = Elements::AtomicMass[m_atoms[i]];
        mass += m;
        pos(0) += m * m_geometry[i][0];
        pos(1) += m * m_geometry[i][1];
        pos(2) += m * m_geometry[i][2];
    }
    pos(0) /= mass;
    pos(1) /= mass;
    pos(2) /= mass;
    for (int i = 0; i < m_geometry.size(); ++i) {
        double scale = 3;
        if (scaling.size() > i)
            scale = scaling[i];
        dipole(0) += m_charges[i] * (m_geometry[i][0] - pos(0)) * scale;
        dipole(1) += m_charges[i] * (m_geometry[i][1] - pos(1)) * scale;
        dipole(2) += m_charges[i] * (m_geometry[i][2] - pos(2)) * scale;
        // std::cout << scale << " ";
    }
    // std::cout << std::endl;
    return dipole;
}

std::pair<double, double> Molecule::GyrationRadius(bool protons, int fragment)
{
    Eigen::Vector3d com = COM(protons, fragment);
    double gyr = 0, gyr_mass = 0;
    for (int i = 0; i < m_geometry.size(); ++i) {
        gyr += ((com(0) - m_geometry[i][0]) * (com(0) - m_geometry[i][0]) + (com(1) - m_geometry[i][1]) * (com(1) - m_geometry[i][1]) + (com(2) - m_geometry[i][2]) * (com(2) - m_geometry[i][2]));
        gyr_mass += Elements::AtomicMass[m_atoms[i]] * ((com(0) - m_geometry[i][0]) * (com(0) - m_geometry[i][0]) + (com(1) - m_geometry[i][1]) * (com(1) - m_geometry[i][1]) + (com(2) - m_geometry[i][2]) * (com(2) - m_geometry[i][2]));
    }
    gyr /= double(m_geometry.size());
    gyr_mass /= double(m_mass);
    return std::pair<double, double>(gyr, gyr_mass);
}

std::pair<int, Position> Molecule::Atom(int i) const
{
    return std::pair<int, Position>(m_atoms[i], { m_geometry[i][0], m_geometry[i][1], m_geometry[i][2] });
}

void Molecule::writeXYZFile(const std::string& filename) const
{
    std::ofstream input;
    input.open(filename, std::ios::out);
    input << XYZString();
    input.close();
}

void Molecule::writeXYZFile(const std::string& filename, const  std::vector<int> &order) const
{
    std::ofstream input;
    input.open(filename, std::ios::out);
    input << XYZString(order);
    input.close();
}

std::string Molecule::Header() const
{
#ifdef GCC
    return fmt::format("{} ** Energy = {:10f} Eh ** Charge = {} ** Spin = {} ** Curcuma {} ({})\n", m_name, Energy(), Charge(), Spin(), qint_version, git_tag);
#else
    return fmt::format("{} ** Energy = {:} Eh ** Charge = {} ** Spin = {} ** Curcuma {} ({})\n", m_name, Energy(), Charge(), Spin(), qint_version, git_tag);
#endif
}

std::string Molecule::Atom2String(int i) const
{
#ifdef GCC
    return fmt::format("{}  {:f}    {:f}    {:f}\n", Elements::ElementAbbr[m_atoms[i]].c_str(), m_geometry[i][0], m_geometry[i][1], m_geometry[i][2]);
#else
    return fmt::format("{}  {:}    {:}    {:}\n", Elements::ElementAbbr[m_atoms[i]].c_str(), m_geometry[i][0], m_geometry[i][1], m_geometry[i][2]);
#endif
}

void Molecule::writeXYZFragments(const std::string& filename) const
{
    /* For now we write only single fragments and pairs, more complex structures (quarternary etc) will follow some time later */
    auto fragments = GetFragments();
    for (int frag = 0; frag < fragments.size(); ++frag) {
        auto fragment = fragments[frag];
        std::ofstream input;
        input.open(filename + "_F" + std::to_string(frag + 1) + ".xyz", std::ios::out);

        std::string output;
        output += Header();
        int atoms = 0;
        for (int j = 0; j < fragment.size(); ++j) {
            int i = fragment[j];
            atoms++;
            output += Atom2String(i);
        }
        input << atoms << std::endl;
        input << output;

        input.close();
        for (int frag2 = frag + 1; frag2 < fragments.size(); ++frag2) {
            auto fragment2 = fragments[frag2];
            for (auto a : fragment)
                fragment2.push_back(a);

            std::ofstream input;
            input.open(filename + "_F" + std::to_string(frag + 1) + "_F" + std::to_string(frag2 + 1) + ".xyz", std::ios::out);

            std::string output;
            output += Header();
            int atoms = 0;
            for (int j = 0; j < fragment2.size(); ++j) {
                int i = fragment2[j];
                atoms++;
                output += Atom2String(i);
            }
            input << atoms << std::endl;
            input << output;

            input.close();
        }
    }
}

void Molecule::appendXYZFile(const std::string& filename) const
{
    std::string output;
    output += fmt::format("{}\n", AtomCount());
    output += Header();
    for (int i = 0; i < AtomCount(); ++i) {
        output += Atom2String(i);
    }
    std::ofstream input;
    input.open(filename, std::ios_base::app);
    //std::cout << output << std::endl;
    input << output;
    input.close();
}

std::string Molecule::XYZString() const
{
    std::string output;
    output += fmt::format("{}\n", AtomCount());
    output += Header();
    for (int i = 0; i < AtomCount(); ++i) {
        output += Atom2String(i);
    }
    return output;
}

std::string Molecule::XYZString(const std::vector<int> &order) const
{
    std::string output;
    output += fmt::format("{}\n", AtomCount());
    output += Header();
    for (int i : order) {
        std::cout << i << " ";
        output += Atom2String(i);
    }
    return output;
}

std::vector<int> Molecule::BoundHydrogens(int atom, double scaling) const
{
    std::vector<int> result;
    if (atom >= AtomCount() || Atom(atom).first == 1)
        return result;

    for (int i = 0; i < AtomCount(); ++i) {
        // std::cout << Atom(i).first << std::endl;
        if (atom == i || Atom(i).first != 1)
            continue;
        double distance = CalculateDistance(i, atom);
        // std::cout << i << " " << atom << " " << distance << std::endl;
        if (distance < (Elements::CovalentRadius[Atom(atom).first] + Elements::CovalentRadius[1]) * scaling)
            result.push_back(i);
    }
    return result;
}

void Molecule::Center(bool mass)
{
    if (!mass)
        setGeometry(GeometryTools::TranslateGeometry(getGeometry(), GeometryTools::Centroid(getGeometry()), Position{ 0, 0, 0 }));
    else
        setGeometry(GeometryTools::TranslateGeometry(getGeometry(), MassCentroid(), Position{ 0, 0, 0 }));
}

void Molecule::CalculateRotationalConstants()
{
    m_Ia = 0, m_Ib = 0, m_Ic = 0;
    double mass = 0;
    Position pos = { 0, 0, 0 };
    for (int i = 0; i < AtomCount(); ++i) {
        double m = Elements::AtomicMass[m_atoms[i]];
        mass += m;
        pos(0) += m * m_geometry[i][0];
        pos(1) += m * m_geometry[i][1];
        pos(2) += m * m_geometry[i][2];
    }
    pos(0) /= mass;
    pos(1) /= mass;
    pos(2) /= mass;

    Geometry geom = GeometryTools::TranslateGeometry(getGeometry(), pos, { 0, 0, 0 });
    Geometry matrix = Geometry::Zero(3, 3);
    for (int i = 0; i < AtomCount(); ++i) {
        double m = Elements::AtomicMass[m_atoms[i]];
        double x = geom(i, 0);
        double y = geom(i, 1);
        double z = geom(i, 2);
        double x2 = x * x;
        double y2 = y * y;
        double z2 = z * z;
        matrix(0, 0) += m * (y2 + z2);
        matrix(1, 1) += m * (x2 + z2);
        matrix(2, 2) += m * (x2 + y2);
        matrix(0, 1) -= m * x * y;
        matrix(0, 2) -= m * x * z;
        matrix(1, 2) -= m * y * z;
    }
    matrix(1, 0) = matrix(0, 1);
    matrix(2, 0) = matrix(0, 2);
    matrix(2, 1) = matrix(1, 2);

    Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
    diag_I.compute(matrix);
    double conv = 1.6605402E-24 * 10E-10 * 10E-10 * 10;
    double conv2 = 6.6260755E-34 / pi / pi / 8;

    m_alignmentAxes = diag_I.eigenvectors();

    m_Ia = conv2 / (diag_I.eigenvalues()(0) * conv);
    m_Ib = conv2 / (diag_I.eigenvalues()(1) * conv);
    m_Ic = conv2 / (diag_I.eigenvalues()(2) * conv);
}

void Molecule::AlignAxis()
{
    CalculateRotationalConstants();
    Geometry coordinates = getGeometry(true);
    {
        Eigen::Vector3d alignmentAxis = m_alignmentAxes.col(0);
        Eigen::Vector3d referenceAxis = Eigen::Vector3d::UnitZ(); // Specify the desired reference axis (e.g., z-axis)
        Eigen::Quaterniond alignmentQuaternion;
        alignmentQuaternion.setFromTwoVectors(alignmentAxis, referenceAxis);
        coordinates = (alignmentQuaternion.toRotationMatrix() * getGeometry().transpose()).transpose();
    }
    {
        Eigen::Vector3d alignmentAxis = m_alignmentAxes.col(1);
        Eigen::Vector3d referenceAxis = Eigen::Vector3d::UnitY(); // Specify the desired reference axis (e.g., z-axis)
        Eigen::Quaterniond alignmentQuaternion;
        alignmentQuaternion.setFromTwoVectors(alignmentAxis, referenceAxis);
        coordinates = (alignmentQuaternion.toRotationMatrix() * getGeometry().transpose()).transpose();
    }
    {
        Eigen::Vector3d alignmentAxis = m_alignmentAxes.col(2);
        Eigen::Vector3d referenceAxis = Eigen::Vector3d::UnitX(); // Specify the desired reference axis (e.g., z-axis)
        Eigen::Quaterniond alignmentQuaternion;
        alignmentQuaternion.setFromTwoVectors(alignmentAxis, referenceAxis);
        coordinates = (alignmentQuaternion.toRotationMatrix() * getGeometry().transpose()).transpose();
    }
    setGeometry(coordinates);
}

std::map<int, std::vector<int>> Molecule::getConnectivtiy(double scaling, int latest) const
{
    if (latest == -1 || latest > AtomCount())
        latest = AtomCount();
    std::map<int, std::vector<int>> connections;
    for (int i = 0; i < AtomCount(); ++i) {
        std::vector<int> connect = BoundHydrogens(i, scaling);
        if (connect.size() == 0)
            continue;

        connections.insert(std::pair<int, std::vector<int>>(i, connect));
    }
    return connections;
}

void Molecule::PrintConnectivitiy(double scaling) const
{
    print_geom();

    for (std::size_t i = 0; i < AtomCount(); ++i) {
        std::cout << Atom(i).first << " ... " << Atom(i).second.transpose() << std::endl;
        if (Atom(i).first != 1) {
            for (std::size_t j = i + 1; j < AtomCount(); ++j) {
                if (Atom(j).first == 1) {
                    double distance = CalculateDistance(i, j);
                    if (distance < (Elements::CovalentRadius[Atom(i).first] + Elements::CovalentRadius[Atom(j).first]) * scaling)
                        std::cout << "Atom " << i << " and Atom " << j << ": Distance = " << distance << " - Cov Rad: " << (Elements::CovalentRadius[Atom(i).first] + Elements::CovalentRadius[Atom(j).first]) << std::endl;
                }
            }
        }
    }
}

void Molecule::AnalyseIntermoleculeDistance() const
{
    double cutoff = 2.5;
    for (std::size_t i = 0; i < m_fragments.size(); ++i) {
        for (std::size_t j = i + 1; j < m_fragments.size(); ++j) {
            for (int a : m_fragments[i]) {
                for (int b : m_fragments[j]) {
                    double distance = CalculateDistance(a, b);
                    if (distance < cutoff) {
                        if (Atom(a).first == 1 && Atom(b).first != 1 || Atom(a).first != 1 && Atom(b).first == 1)
                            std::cout << std::setprecision(6) << distance << " " << a << "(" << Atom(a).first << ") - " << b << "(" << Atom(b).first << ")" << std::endl;
                    }
                }
            }
        }
    }
}

std::vector<std::vector<int>> Molecule::GetFragments(double scaling) const
{
    if (scaling != m_scaling)
        m_dirty = true;
    if (m_fragments.size() > 0 && !m_dirty)
        return m_fragments;
    m_mass_fragments.clear();
    m_scaling = scaling;
    std::multimap<double, std::vector<int>> ordered_list;

    std::vector<int> fragment, atoms, cached;
    double mass = 0;
    for (std::size_t i = 0; i < m_atoms.size(); ++i)
        atoms.push_back(i);

    m_fragments.clear();

    while (atoms.size()) {
        fragment.push_back(atoms.at(0));
        atoms.erase(atoms.begin());
        for (std::size_t i = 0; i < fragment.size(); ++i) {
            for (std::size_t j = 0; j < atoms.size(); ++j) {
                double distance = CalculateDistance(fragment[i], atoms[j]);
                //std::cout << "Fragment No: " << m_fragments.size() + 1 << " ("<< fragment.size() <<")  ##  Atom " << fragment[i] << " (Index " << i << ") and Atom " << atoms[j] << " - Distance: " << distance << " Thresh " << (Elements::CovalentRadius[Atom(fragment[i]).first] + Elements::CovalentRadius[Atom(atoms[j]).first]);
                if (distance < (Elements::CovalentRadius[Atom(fragment[i]).first] + Elements::CovalentRadius[Atom(atoms[j]).first]) * m_scaling) {
                    //std::cout << " ... taken " << std::endl;
                    fragment.push_back(atoms.at(j));
                    mass += Elements::AtomicMass[Atom(atoms.at(j)).first];
                    atoms.erase(atoms.begin() + j);
                    --j; // We have to decrease the atom vector counter since it got one element removed and we dont want to skip single elements
                } //else
                //   std::cout  << std::endl;
            }
        }
        std::sort(fragment.begin(), fragment.end());
        //m_fragments.push_back(fragment);
        mass *= -1; // I think, that is an easy way to ***
        ordered_list.insert(std::pair<double, std::vector<int>>(mass, fragment));
        fragment.clear();
        mass = 0;
    }
    for (const auto& entry : ordered_list) {
        //m_mass_fragments.push_back(-1 * entry.first); // *** make the std::map container sort in reverse order :-)
        m_fragments.push_back(entry.second);
    }
    for (int i = 0; i < m_fragments.size(); ++i) {
        double mass = 0;
        for (auto atom : m_fragments[i]) {
            mass += Elements::AtomicMass[Atom(atom).first];
            m_fragment_assignment.insert(std::pair<int, int>(atom, i));
        }
        m_mass_fragments.push_back(mass);
    }

    m_dirty = false;
    return m_fragments;
}

void Molecule::InitialiseConnectedMass(double scaling, bool protons)
{
    for (int i = 0; i < AtomCount(); ++i) {
        int mass = 0;
        auto atom_i = Atom(i);
        if (!protons && atom_i.first == 1)
            continue;
        for (int j = i + 1; j < AtomCount(); ++j) {
            auto atom_j = Atom(j);
            double distance = CalculateDistance(i, j);
            if (distance < (Elements::CovalentRadius[atom_i.first] + Elements::CovalentRadius[atom_j.first]) * scaling) {
                mass += atom_j.first; //Elements::AtomicMass[atom_j.first - 1];
            }
        }
        m_connect_mass.push_back(mass);
    }
}

std::vector<int> Molecule::WhiteListProtons() const
{
    std::vector<int> whitelist_proton;
    for (std::size_t i = 0; i < AtomCount(); ++i) {
        if (Atom(i).first != 1)
            continue;
        double distance = 1e8;
        int element = 0;
        for (std::size_t j = 0; j < AtomCount(); ++j) {
            if (i == j)
                continue;

            double d = CalculateDistance(i, j);
            if (d < distance)
                element = Atom(j).first;
            distance = std::min(d, distance);
        }
        if (element == 7 || element == 8) {
            whitelist_proton.push_back(i);
        }
    }
    return whitelist_proton;
}

void Molecule::MapHydrogenBonds()
{
    std::vector<int> whitelist_proton = WhiteListProtons();
    m_HydrogenBondMap = Matrix::Zero(AtomCount(), AtomCount());
    double h_radius = Elements::CovalentRadius[1];
    for (int hydrogen : whitelist_proton) {
        int accepted_donor = -1;
        double accepted_distance = m_hbond_cutoff;
        for (int donor = 0; donor < AtomCount(); ++donor) {
            int element = Atom(donor).first;
            if (!(element == 8 || element == 7))
                continue;

            double distance = CalculateDistance(donor, hydrogen);
            if (distance > (Elements::CovalentRadius[element] + h_radius) * m_scaling && distance < accepted_distance)
                accepted_donor = donor;
        }
        if (accepted_donor > -1) {
            m_HydrogenBondMap(accepted_donor, hydrogen) = 1;
            m_HydrogenBondMap(hydrogen, accepted_donor) = 1;
        }
    }
}

Matrix Molecule::HydrogenBondMatrix(int f1, int f2)
{
    if (m_HydrogenBondMap.rows() != AtomCount())
        MapHydrogenBonds();

    if (f1 == f2 && f1 == -1)
        return m_HydrogenBondMap;

    Matrix HydrogenBondMap = Matrix::Zero(AtomCount(), AtomCount());

    for (int i = 0; i < AtomCount(); ++i) {
        for (int j = i + 1; j < AtomCount(); ++j) {
            HydrogenBondMap(i, j) = m_HydrogenBondMap(i, j) * // has to be identified as hydrogen bond before
                int((m_fragment_assignment[i] == f1) || (m_fragment_assignment[j] == f1) || (f1 == -1)) * // f1 is either i, j or -1
                int((m_fragment_assignment[i] == f2) || (m_fragment_assignment[j] == f2) || (f2 == -1)); // f2 is either i, j or -1
            HydrogenBondMap(j, i) = HydrogenBondMap(i, j);
        }
    }
    return HydrogenBondMap;
}

Molecule Molecule::ElementsRemoved(const std::vector<int>& elements)
{
    Molecule mol;
    for (int i = 0; i < AtomCount(); ++i) {
        auto pair = Atom(i);
        if (std::find(elements.begin(), elements.end(), pair.first) == elements.end())
            mol.addPair(pair);
    }
    return mol;
}

Molecule Molecule::AtomsRemoved(const std::vector<int>& atoms)
{
    Molecule mol;
    for (int i = 0; i < AtomCount(); ++i) {
        auto pair = Atom(i);
        if (std::find(atoms.begin(), atoms.end(), i) == atoms.end())
            mol.addPair(pair);
    }
    return mol;
}

std::pair<Matrix, Matrix> Molecule::DistanceMatrix() const
{
    Matrix distance = Eigen::MatrixXd::Zero(AtomCount(), AtomCount());
    Matrix topo = distance;
    for (int i = 0; i < AtomCount(); ++i) {
        for (int j = 0; j < i; ++j) {
            distance(i, j) = CalculateDistance(i, j);
            distance(j, i) = distance(i, j);
            topo(i, j) = distance(i, j) <= (Elements::CovalentRadius[Atom(i).first] + Elements::CovalentRadius[Atom(j).first]) * m_scaling;
            topo(j, i) = topo(i, j);
        }
    }
    return std::pair<Matrix, Matrix>(distance, topo);
}

std::vector<double> Molecule::GetBox() const
{
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;

    for (int i = 0; i < AtomCount(); ++i) {
        auto pair = Atom(i);
        x_min = std::min(pair.second[0], x_min);
        x_min = std::max(pair.second[0], x_max);

        y_min = std::min(pair.second[1], y_min);
        y_min = std::max(pair.second[1], y_max);

        z_min = std::min(pair.second[2], z_min);
        z_min = std::max(pair.second[2], z_max);
    }
    std::vector<double> box;
    box.push_back(x_max - x_min);
    box.push_back(y_max - y_min);
    box.push_back(z_max - z_min);
    return box;
}
