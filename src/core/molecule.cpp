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

#include <cstdio>
#include <fstream>
#include <sstream>
#include <istream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <array>

#include "molecule.h"

void Molecule::print_geom() const
{
    std::cout << AtomCount() << endl
              << endl;
    for (int i = 0; i < AtomCount(); i++) {
        printf("%s %8.5f %8.5f %8.5f\n", Elements::ElementAbbr[m_atoms[i]].c_str(), geom[i][0], geom[i][1], geom[i][2]);
    }
}

void Molecule::printAtom(int i) const
{
    if (i < AtomCount())
        printf("%s %8.5f %8.5f %8.5f", m_atoms[i], geom[i][0], geom[i][1], geom[i][2]);
}


void Molecule::translate(double x, double y, double z)
{
    for (int i = 0; i < AtomCount(); i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

Molecule::Molecule(int n, int q)
{
    m_charge = q;
    InitialiseEmptyGeometry(n);
}

Molecule::Molecule()
{
}

Molecule::~Molecule()
{ 

}

void Molecule::InitialiseEmptyGeometry(int atoms)
{
    for (int i = 0; i < atoms; i++) {
        std::array<double, 3> atom = { 0, 0, 0 };
        geom.push_back(atom);
    }
}

bool Molecule::addPair(const std::pair<int, Position>& atom)
{
    bool exist = true;
    // const std::array<double, 3> at = { atom.second(0),  atom.second(1),  atom.second(2)};
    geom.push_back({ atom.second(0), atom.second(1), atom.second(2) });
    m_atoms.push_back(atom.first);

    for (int i = 0; i < AtomCount(); ++i)
        for (int j = i + 1; j < AtomCount(); ++j)
            if (Distance(i, j) < 1e-6)
                exist = false;
    return exist;
}

double Molecule::Distance(int i, int j) const
{
    if (i >= AtomCount() || j >= AtomCount())
        return 0;
    
    double x_i = geom[i][0];
    double x_j = geom[j][0];
    
    double y_i = geom[i][1];
    double y_j = geom[j][1];
    
    double z_i = geom[i][2];
    double z_j = geom[j][2];
    
    return sqrt( ( ((x_i-x_j)*(x_i-x_j))+((y_i-y_j)*(y_i-y_j))+((z_i-z_j)*(z_i-z_j)) ));
    
}

double Molecule::DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const
{
    return pos1[0]*pos2[0]+pos1[1]*pos2[1]+pos1[2]*pos2[2];
}


double Molecule::angle(int atom1, int atom2, int atom3) const
{

    std::array<double, 3> atom_0 = { geom[atom2 - 1] }; // Proton
    std::array<double, 3> atom_1 = { geom[atom3 - 1] }; // Acceptor
    std::array<double, 3> atom_2 = { geom[atom1 - 1] }; // Donor
    
    
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
        geom[i][0] = x;
        geom[i][1] = y;
        geom[i][2] = z;
    }
        
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
        
        geom[i][0] = x;
        geom[i][1] = y;
        geom[i][2] = z;
    }
        
}

void Molecule::clear()
{
    m_atoms.clear();
    geom.clear();
}

void Molecule::LoadMolecule(const Molecule& molecule)
{
    clear();
    m_charge = molecule.Charge();
    m_atoms = molecule.Atoms();
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

Geometry Molecule::getGeometry() const
{
    Geometry geometry(geom.size(), 3);

    for(int i = 0; i < geom.size(); ++i)
    {
        geometry(i, 0) = geom[i][0];
        geometry(i, 1) = geom[i][1];
        geometry(i, 2) = geom[i][2];
    }
    return geometry;
}

bool Molecule::setGeometry(const Geometry &geometry)
{
    if(geometry.rows() != geom.size())
        return false;

    for(int i = 0; i < geom.size(); ++i)
    {
        geom[i][0] = geometry(i, 0);
        geom[i][1] = geometry(i, 1);
        geom[i][2] = geometry(i, 2);
    }

    return true;
}

Position  Molecule::Centroid(bool hydrogen) const
{
    Position position = {0 , 0 , 0};
    Geometry geom = getGeometry();
    if(hydrogen)
    {
        for(int i = 0; i < geom.rows(); ++i)
        {
            position += geom.row(i);
        }

        position /= double(geom.rows());
    }else
    {
        throw -1;
    }
    return position;
}

std::pair<int, Position> Molecule::Atom(int i) const
{
    return std::pair<int, Position>(m_atoms[i], { geom[i][0], geom[i][1], geom[i][2] });
}

void Molecule::writeXYZFile(const std::string& filename)
{
    std::ofstream input;
    input.open(filename, ios::out);
    input << AtomCount() << std::endl
          << std::endl;
    for (int i = 0; i < AtomCount(); ++i) {
        input << Elements::ElementAbbr[m_atoms[i]].c_str() << "\t" << geom[i][0] << "\t" << geom[i][1] << "\t" << geom[i][2] << std::endl;
    }
    input.close();
}
