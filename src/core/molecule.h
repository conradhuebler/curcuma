/*
 * <Internal Coordinate Handler for chemical structures.>
 * Copyright (C) 2019  Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <string>
#include <vector> 
#include <array>

#include <Eigen/Dense>

#include "src/core/global.h"

using namespace std;
 
class Molecule
{
  public:
    Molecule(int n, int q = 0);
    ~Molecule();

    void print_geom() const;
    void printAtom(int i) const;
//     void rotate(double phi);
    void translate(double x, double y, double z);
//     double bond(int atom1, int atom2) const;
    double angle(int atom1, int atom2, int atom3) const;
    double DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const;
//     double torsion(int atom1, int atom2, int atom3, int atom4);
 


    void setAtom(const std::string &internal, int i);
    void setXYZ(const std::string &coord, int i);
    double Distance(int i, int j) const;

    Geometry getGeometry() const;
    bool setGeometry(const Geometry &geometry);
    Position Centroid(bool hydrogen = true) const;

private:
    int natom;
    int charge;
    int *zvals;
    std::vector<std::string > atoms;
    std::vector< std::array<double, 3> > geom;
    std::vector< int > m_elements;

//     double **geom;
    string point_group;
};
