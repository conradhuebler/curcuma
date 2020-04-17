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
}

Molecule::Molecule(const Molecule* other)
{
    m_geometry = other->m_geometry;
    m_charge = other->m_charge;
    m_fragments = other->m_fragments;
    m_atoms = other->m_atoms;
    m_name = other->m_name;
    m_energy = other->m_energy;
}

Molecule::~Molecule()
{
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
        printf("%s %8.5f %8.5f %8.5f", m_atoms[i], m_geometry[i][0], m_geometry[i][1], m_geometry[i][2]);
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

    for (std::size_t i = 0; i < AtomCount(); ++i)
        for (std::size_t j = i + 1; j < AtomCount(); ++j)
            if (Distance(i, j) < 1e-6)
                exist = false;

    m_dirty = true;

    return exist;
}

void Molecule::setXYZComment(const std::string& comment)
{
    StringList list = Tools::SplitString(comment);
    if (list.size() == 7) {
        if (list[2].compare("gnorm:") == 0) {
            //setEnergy(std::stod((list[1])));
        } else
            try {
                setEnergy(std::stod((list[4])));
            } catch (const std::string& what_arg) {
                setEnergy(0);
            }
    } else if (list.size() == 4) {
        if (list[0].compare("SCF") == 0 && list[1].compare("done") == 0) {
            try {
                setEnergy(std::stod((list[2])));
            } catch (const std::string& what_arg) {
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
                }
            }
        }
    } else if (list.size() == 1) {
        try {
            setEnergy(std::stod((list[0])));
        } catch (const std::string& what_arg) {
        }
    } else {
        for (const string& s : list) {
            double energy = 0;
            if (Tools::isDouble(s)) {
                energy = std::stod(s);
                setEnergy(energy);
                break;
            }
        }
    }
}

bool Molecule::Contains(const std::pair<int, Position>& atom)
{

    for (std::size_t i = 0; i < AtomCount(); ++i) {
        if (GeometryTools::Distance(Atom(i).second, atom.second) < 1e-6)
            return true;
    }

    return false;
}

double Molecule::Distance(int i, int j) const
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


double Molecule::angle(int atom1, int atom2, int atom3) const
{

    std::array<double, 3> atom_0 = { m_geometry[atom2 - 1] }; // Proton
    std::array<double, 3> atom_1 = { m_geometry[atom3 - 1] }; // Acceptor
    std::array<double, 3> atom_2 = { m_geometry[atom1 - 1] }; // Donor

    std::array<double, 3> vec_1 = { atom_0[0]-atom_1[0], atom_0[1]-atom_1[1], atom_0[2]-atom_1[2] };
    std::array<double, 3> vec_2 = { atom_0[0]-atom_2[0], atom_0[1]-atom_2[1], atom_0[2]-atom_2[2] };
    
    return acos(DotProduct(vec_1,vec_2)/(sqrt(DotProduct(vec_1,vec_1)*DotProduct(vec_2,vec_2))))*360/2.0/pi;
}



void Molecule::setAtom(const std::string& internal, int i)
{
    std::vector<std::string > elements;
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
    m_atoms.push_back(Elements::String2Element(elements[0]));

    if(elements.size() == 7)
    {
        
        int atom_1 = stoi(elements[1]);
        int atom_2 = stoi(elements[2]);
        int atom_3 = stoi(elements[3]);
        
        double r_ij = stod(elements[4]);
        double omega = stod(elements[5]);
        double theta = stod(elements[6]);
        
        double r_ik = Distance(atom_2, atom_1);
        
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
    m_atoms.push_back(Elements::String2Element(elements[0]));

    if(elements.size() >= 4)
    {
        
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
    InitialiseEmptyGeometry(molecule.AtomCount());
    setGeometry(molecule.getGeometry());
}

void Molecule::LoadMolecule(const Molecule* molecule)
{
    clear();
    m_charge = molecule->Charge();
    m_atoms = molecule->Atoms();
    InitialiseEmptyGeometry(molecule->AtomCount());
    setGeometry(molecule->getGeometry());
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
    ;
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
    ;
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

Position Molecule::Centroid(bool protons, int fragment) const
{
    return GeometryTools::Centroid(getGeometryByFragment(fragment, protons));
}

std::pair<int, Position> Molecule::Atom(int i) const
{
    return std::pair<int, Position>(m_atoms[i], { m_geometry[i][0], m_geometry[i][1], m_geometry[i][2] });
}

void Molecule::writeXYZFile(const std::string& filename) const
{
    std::ofstream input;
    input.open(filename, ios::out);
    input << AtomCount() << std::endl
          << Name() << " ** Energy = " << std::setprecision(12) << Energy() << " Eh **"
          << std::endl;
    for (int i = 0; i < AtomCount(); ++i) {
        input << Elements::ElementAbbr[m_atoms[i]].c_str() << "      " << m_geometry[i][0] << "      " << m_geometry[i][1] << "      " << m_geometry[i][2] << std::endl;
    }
    input.close();
}

void Molecule::appendXYZFile(const std::string& filename) const
{
    std::ofstream input;
    input.open(filename, std::ios_base::app);
    input << AtomCount() << std::endl
          << Name() << " ** Energy = " << std::setprecision(12) << Energy() << " Eh **"
          << std::endl;
    for (int i = 0; i < AtomCount(); ++i) {
        input << Elements::ElementAbbr[m_atoms[i]].c_str() << "      " << m_geometry[i][0] << "      " << m_geometry[i][1] << "      " << m_geometry[i][2] << std::endl;
    }
    input.close();
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
        double distance = Distance(i, atom);
        // std::cout << i << " " << atom << " " << distance << std::endl;
        if (distance < (Elements::CovalentRadius[Atom(atom).first] + Elements::CovalentRadius[1]) * scaling)
            result.push_back(i);
    }
    return result;
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
    // std::cout << Centroid().transpose() << std::endl << pos.transpose() << std::endl;
    Geometry geom = GeometryTools::TranslateGeometry(getGeometry(), Centroid(), pos);
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
        matrix(0, 1) += m * x * y;
        matrix(0, 2) += m * x * z;
        matrix(1, 2) += m * y * z;
    }
    matrix(1, 0) = matrix(0, 1);
    matrix(2, 0) = matrix(0, 2);
    matrix(2, 1) = matrix(1, 2);

    Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
    diag_I.compute(matrix);

    //std::cout << diag_I.eigenvalues().transpose() << std::endl;
    m_Ia = diag_I.eigenvalues()(0);
    m_Ib = diag_I.eigenvalues()(1);
    m_Ic = diag_I.eigenvalues()(2);
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
                    double distance = Distance(i, j);
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
                    double distance = Distance(a, b);
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
                double distance = Distance(fragment[i], atoms[j]);
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
        m_mass_fragments.push_back(-1 * entry.first); // *** make the std::map container sort in reverse order :-)
        m_fragments.push_back(entry.second);
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
            double distance = Distance(i, j);
            if (distance < (Elements::CovalentRadius[atom_i.first - 1] + Elements::CovalentRadius[atom_j.first - 1]) * scaling) {
                //      std::cout << atom_i.first - 1<< " "  << atom_j.first - 1<< " " << Elements::AtomicMass[atom_j.first - 1] << std::endl;
                mass += atom_j.first; //Elements::AtomicMass[atom_j.first - 1];
            }
        }
        m_connect_mass.push_back(mass);
        //std::cout << i << " " << mass << std::endl;
    }
}
