/*
 * <Internal Coordinate Handler for chemical structures.>
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

#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "src/core/global.h"

using namespace std;
 
class Molecule
{
  public:
    Molecule(int n, int q = 0);
    Molecule();
    ~Molecule();

    void print_geom() const;
    void printAtom(int i) const;
    inline int Charge() const { return m_charge; }
    void setCharge(int charge) { m_charge = charge; }

    //     void rotate(double phi);
    void translate(double x, double y, double z);
//     double bond(int atom1, int atom2) const;
    double angle(int atom1, int atom2, int atom3) const;
    double DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const;
//     double torsion(int atom1, int atom2, int atom3, int atom4);

    void clear();

    void LoadMolecule(const Molecule& molecule);
    void LoadMolecule(const Molecule* molecule);

    void setAtom(const std::string &internal, int i);
    void setXYZ(const std::string &coord, int i);
    bool addPair(const std::pair<int, Position>& atom);
    double Distance(int i, int j) const;

    Geometry getGeometry() const;
    bool setGeometry(const Geometry &geometry);
    Position Centroid(bool hydrogen = true) const;
    inline int AtomCount() const { return m_atoms.size(); }
    std::vector<int> Atoms() const { return m_atoms; }
    std::pair<int, Position> Atom(int i) const;

    void writeXYZFile(const std::string& filename);

private:
    void InitialiseEmptyGeometry(int atoms);

    int m_charge = 0;
    int *zvals;
    std::vector< std::array<double, 3> > geom;
    std::vector<int> m_atoms;

    string point_group;
};
