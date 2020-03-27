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
#include <map>
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
    Molecule(const Molecule* other);
    Molecule(const Molecule& other);
    Molecule();
    ~Molecule();

    void print_geom(bool moreinfo = true) const;
    void printFragmente();
    void printAtom(int i) const;
    inline int Charge() const { return m_charge; }
    void setCharge(int charge) { m_charge = charge; }

    double angle(int atom1, int atom2, int atom3) const;
    double DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const;
//     double torsion(int atom1, int atom2, int atom3, int atom4);

    void clear();

    void LoadMolecule(const Molecule& molecule);
    void LoadMolecule(const Molecule* molecule);

    void setAtom(const std::string &internal, int i);
    void setXYZ(const std::string &coord, int i);
    bool addPair(const std::pair<int, Position>& atom);
    bool Contains(const std::pair<int, Position>& atom);

    double Distance(int i, int j) const;
    Geometry getGeometry(const IntPair& pair, bool protons = true) const;
    Geometry getGeometry(std::vector<int> atoms, bool protons = true) const;
    Geometry getGeometryByFragment(int fragment, bool protons = true) const;
    inline Geometry getGeometry(bool protons = true) const { return getGeometry(IntPair(0, -1), protons); }

    bool setGeometry(const Geometry &geometry);
    bool setGeometryByFragment(const Geometry& geometry, int fragment, bool protons = true);

    Position Centroid(bool hydrogen = true, int fragment = -1) const;
    inline std::size_t AtomCount() const { return m_atoms.size(); }
    std::vector<int> Atoms() const { return m_atoms; }
    std::pair<int, Position> Atom(int i) const;

    void writeXYZFile(const std::string& filename) const;
    inline void writeXYZFile() const { writeXYZFile(Name() + ".xyz"); }

    void appendXYZFile(const std::string& filename) const;
    inline void appendXYZFile() const { appendXYZFile(Name() + ".xyz"); }

    std::vector<int> BoundHydrogens(int atom, double scaling = 1.5) const;
    std::map<int, std::vector<int>> getConnectivtiy(double scaling = 1.5, int latest = -1) const;

    void PrintConnectivitiy(double scaling = 1.5) const;
    inline void PrintConnectivitiy() const { PrintConnectivitiy(m_scaling); }

    /*! \brief Return Fragments, will be determined on first call, then stored */
    std::vector<std::vector<int>> GetFragments(double scaling) const;
    inline std::vector<std::vector<int>> GetFragments() const { return GetFragments(m_scaling); }

    inline void setName(const std::string& name) { m_name = name; }

    inline void setEnergy(double energy) { m_energy = energy; }

    inline double Energy() const { return m_energy; }

    inline std::string Name() const { return m_name; }
    void CalculateRotationalConstants();

    inline double Ia() const { return m_Ia; }
    inline double Ib() const { return m_Ib; }
    inline double Ic() const { return m_Ic; }

    void AnalyseIntermoleculeDistance() const;

    inline void setScaling(double scaling) { m_scaling = scaling; }

private:
    void InitialiseEmptyGeometry(int atoms);

    int m_charge = 0;
    std::vector<std::array<double, 3>> m_geometry;
    std::vector<int> m_atoms;
    mutable std::vector<std::vector<int>> m_fragments;
    mutable std::vector<double> m_mass_fragments;
    string point_group;
    mutable bool m_dirty = true;
    std::string m_name;
    double m_energy = 0, m_Ia = 0, m_Ib = 0, m_Ic = 0;
    mutable double m_scaling = 1.5;
};
