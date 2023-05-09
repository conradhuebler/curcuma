/*
 * <Internal Coordinate Handler for chemical structures.>
 * Copyright (C) 2019 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

typedef std::pair<int, Position> AtomDef;

struct Mol {
    double m_energy;
    double m_spin;

    int m_number_atoms;
    int m_charge;

    std::string m_commentline;

    std::vector<std::array<double, 3>> m_geometry;
    std::vector<int> m_atoms;
};

class Molecule
{
  public:
    Molecule(int n, int q = 0);
    Molecule(const Molecule* other);
    Molecule(const Molecule& other);
    Molecule(const std::string& file);
    Molecule(const Mol& mol);
    Molecule(const Mol* mol);

    Molecule();
    ~Molecule();

    /* Molecule& operator=(const Molecule& molecule);
     Molecule& operator=(const Molecule* molecule);
 */
    void print_geom(bool moreinfo = true) const;
    void printFragmente();
    void printAtom(int i) const;
    inline int Charge() const { return m_charge; }
    void setCharge(int charge) { m_charge = charge; }

    inline double Mass() const { return m_mass; }
    double CalculateMass();
    std::vector<double> FragmentMass() const { return m_mass_fragments; }

    void InitialiseConnectedMass(double scaling = 1.3, bool protons = true);
    inline double ConnectedMass(int atom) const { return m_connect_mass[atom]; }
    double CalculateAngle(int atom1, int atom2, int atom3) const;
    double DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const;

    void clear();

    void LoadMolecule(const Molecule& molecule);
    void LoadMolecule(const Molecule* molecule);

    void LoadMolecule(const Mol& molecule);
    void LoadMolecule(const Mol* molecule);

    void setAtom(const std::string &internal, int i);
    void setXYZ(const std::string &coord, int i);

    void setXYZComment(const std::string& comment);

    bool addPair(const std::pair<int, Position>& atom);
    bool Contains(const std::pair<int, Position>& atom);

    Molecule getFragmentMolecule(int fragment) const;

    double CalculateDistance(int i, int j) const;
    Geometry getGeometry(const IntPair& pair, bool protons = true) const;
    Geometry getGeometry(std::vector<int> atoms, bool protons = true) const;
    Geometry getGeometryByFragment(int fragment, bool protons = true) const;
    Geometry getGeometry(bool protons = true) const;

    std::string LowerDistanceMatrix() const;
    std::vector<float> LowerDistanceVector() const;
    std::vector<double> DeltaEN() const;

    bool setGeometry(const Geometry &geometry);
    bool setGeometryByFragment(const Geometry& geometry, int fragment, bool protons = true);

    Position Centroid(bool hydrogen = true, int fragment = -1) const;
    inline std::size_t AtomCount() const { return m_atoms.size(); }
    std::vector<int> Atoms() const { return m_atoms; }
    std::pair<int, Position> Atom(int i) const;

    void writeXYZFile(const std::string& filename) const;
    void writeXYZFile(const std::string& filename, const  std::vector<int> &order) const;

    inline void writeXYZFile() const { writeXYZFile(Name() + ".xyz"); }

    void appendXYZFile(const std::string& filename) const;
    inline void appendXYZFile() const { appendXYZFile(Name() + ".xyz"); }

    std::string XYZString() const;
    std::string XYZString(const std::vector<int> &order) const;


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

    std::string Atom2String(int i) const;
    std::string Header() const;

    void CalculateRotationalConstants();

    inline double Ia() const { return m_Ia; }
    inline double Ib() const { return m_Ib; }
    inline double Ic() const { return m_Ic; }

    void AnalyseIntermoleculeDistance() const;

    inline void setScaling(double scaling) { m_scaling = scaling; }

    void MapHydrogenBonds();
    Matrix HydrogenBondMatrix(int f1, int f2);
    void writeXYZFragments(const std::string& basename) const;

    int Check() const;

    inline void setSpin(int spin) { m_spin = spin; }
    inline int Spin() const { return m_spin; }

    Molecule ElementsRemoved(const std::vector<int>& elements);
    Molecule AtomsRemoved(const std::vector<int>& atoms);

    int CountElement(int element)
    {
        return std::count(m_atoms.cbegin(), m_atoms.cend(), element);
    }

    void setPersisentImage(const Eigen::MatrixXd image)
    {
        m_persistentImage = image;
    }

    Eigen::MatrixXd getPersisentImage() const
    {
        return m_persistentImage;
    }

    void Center();

    std::pair<Matrix, Matrix> DistanceMatrix() const;

    std::vector<std::array<double, 3>> Coords() const { return m_geometry; }

    Matrix RotationMatrix() const { return m_rotation_matrix; }

private:
    void ParseString(const std::string& internal, std::vector<std::string>& elements);

    bool setXYZComment_0(const StringList& list);
    bool setXYZComment_1(const StringList& list);
    bool setXYZComment_2(const StringList& list);
    bool setXYZComment_3(const StringList& list);
    bool setXYZComment_4(const StringList& list);
    bool setXYZComment_5(const StringList& list);
    bool setXYZComment_6(const StringList& list);
    bool setXYZComment_7(const StringList& list);
    bool setXYZComment_8(const StringList& list);
    bool setXYZComment_10(const StringList& list);

    std::vector<int> WhiteListProtons() const;

    void InitialiseEmptyGeometry(int atoms);

    int m_charge = 0, m_spin = 0;
    std::vector<std::array<double, 3>> m_geometry;
    std::vector<int> m_atoms;

    std::vector<int> m_connect_mass;
    Matrix m_HydrogenBondMap;
    Eigen::MatrixXd m_persistentImage, m_rotation_matrix;
    mutable std::vector<std::vector<int>> m_fragments;
    mutable std::map<int, int> m_fragment_assignment;

    mutable std::vector<double> m_mass_fragments;
    mutable bool m_dirty = true;
    std::string m_name;
    double m_energy = 0, m_Ia = 0, m_Ib = 0, m_Ic = 0, m_mass = 0, m_hbond_cutoff = 3;
    mutable double m_scaling = 1.5;
};
