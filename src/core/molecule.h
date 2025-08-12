/*
 * <Internal Coordinate Handler for chemical structures.>
 * Copyright (C) 2019 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
 *               2024 Gerd Gehrisch
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

#include "json.hpp"
using json = nlohmann::json;

typedef std::pair<int, Position> AtomDef;
namespace curcuma {
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

    Molecule& operator=(const Molecule& molecule) = default;
    Mol getMolInfo() const;

    json ExportJson() const;
    void WriteJsonFile(const std::string& filename);

    void ImportJson(const std::string& jsonfile);
    void ImportJson(const json& molecule);

    void Initialise(const int* attyp, const double* coord, const int natoms, const double charge, const int spin);
    void ApplyReorderRule(const std::vector<int>& rule);

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
    // double DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const;

    void clear();

    void LoadMolecule(const Molecule& molecule);
    void LoadMolecule(const Molecule* molecule);

    void LoadMolecule(const Mol& molecule);
    void LoadMolecule(const Mol* molecule);

    void setAtom(const std::string &internal, int i);
    void setXYZ(const std::string &coord, int i);

    void setXYZComment(const std::string& comment);

    bool addPair(const std::pair<int, Position>& atom);
    inline bool addAtom(const std::pair<int, Position>& atom){ return addPair(atom); }

    bool Contains(const std::pair<int, Position>& atom);

    std::vector<double> GetBox() const;

    Molecule getFragmentMolecule(int fragment) const;
    Molecule getFragmentMolecule(const std::vector<int>& atoms) const;

    double CalculateDistance(int i, int j) const;
    std::pair<double, double> GyrationRadius(double hmass = 1, bool hydrogen = true, int fragment = -1);

    /*! \brief Methode to get the geometry of the molecule as array of vectors
     *
     */
    Geometry getGeometry(const IntPair& pair, bool protons = true) const;
    Geometry getGeometry(std::vector<int> atoms, bool protons = true) const;
    Geometry getGeometryByFragment(int fragment, bool protons = true) const;
    Geometry getGeometry(bool protons = true) const;

    std::string LowerDistanceMatrix(bool exclude_bonds = false, bool print_elements = false) const;
    std::string DistanceMatrixString(bool exclude_bonds = false, bool print_elements = false) const;

    std::vector<float> LowerDistanceVector(bool exclude_hydrogen = false) const;
    std::vector<double> DeltaEN() const;

    bool setGeometry(const Geometry &geometry);
    bool setGeometryByFragment(const Geometry& geometry, int fragment, bool protons = true);

    Position Centroid(bool hydrogen = true, int fragment = -1) const;

    /*! \brief Method to calc the center of Mass
     *
     */
    Position MassCentroid(bool hydrogen = true, int fragment = -1) const;
    Eigen::Vector3d COM(bool hydrogen = true, int fragment = -1);

    /*! \brief Methode to get number of atoms
     *
     * @return size of the atoms
     */
    inline std::size_t AtomCount() const { return m_atoms.size(); }

    /*! \brief Methode to get array of all the atoms
     *
     * @return array of the atomnumber
     */
    std::vector<int> Atoms() const { return m_atoms; }
    std::vector<int> FragString2Indicies(const std::string& string) const;

    /*! \brief Methode to get atom number and XYZ from index
     *
     * @param i: index of the atom
     * @return pair of atom number and the xyz position
     */

    std::pair<int, Position> Atom(int i) const;

    void writeXYZFile(const std::string& filename) const;
    void writeXYZFile(const std::string& filename, const  std::vector<int> &order) const;

    inline void writeXYZFile() const { writeXYZFile(Name() + ".xyz"); }

    void appendXYZFile(const std::string& filename) const;
    inline void appendXYZFile() const { appendXYZFile(Name() + ".xyz"); }

    void appendDipoleFile(const std::string& filename) const;
    inline void appendDipoleFile() const { appendDipoleFile(Name() + ".dip"); }

    std::string XYZString() const;
    std::string XYZString(const std::vector<int> &order) const;

    /*! \brief Methode to calculate the dipole moments of single molecules
     *
     */
    std::vector<Position> CalculateDipoleMoments(const std::vector<double>& scaling = std::vector<double>(), const std::vector<std::vector<int>>& fragments = std::vector<std::vector<int>>()) const;

    /*! \brief Methode to calculate the dipole moments of whole structure
     * unit of dipol is electron times angstron
     */
    Position CalculateDipoleMoment(const Vector& scaling = Vector(), bool bond = false) const;
    Position CalculateDipoleMoment(const std::vector<double>& scaling = std::vector<double>(), bool bond = false) const;

    Geometry ChargeDistribution() const;

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
    // Claude Generated: Get molecular formula
    std::string Formula() const;

    /*! \brief Methode to get the atom name with position
     *
     * Input i is the number of the atom,
     * Output is a string as in .xyz file
     * */
    std::string Atom2String(int i) const;

    /*! \brief Methode to get the header of the xyz file
     *
     * */
    std::string Header() const;

    void CalculateRotationalConstants();
    void AlignAxis(const std::vector<int>& axisPermutation, const std::vector<int>& axisOrientation);

    inline double Ia() const { return m_Ia; }
    inline double Ib() const { return m_Ib; }
    inline double Ic() const { return m_Ic; }

    void AnalyseIntermoleculeDistance() const;

    inline void setScaling(double scaling) { m_scaling = scaling; }

    void MapHydrogenBonds();
    Matrix HydrogenBondMatrix(int f1, int f2);
    void writeXYZFragments(const std::string& basename) const;

    /*! no use at the moment
     *
     */
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

    /*! Methode to translate the coord to the centroid
     *
     */
    void Center(bool mass = false);

    std::pair<Matrix, Matrix> DistanceMatrix() const;

    Geometry Coords() const { return m_geometry; }

    Matrix AlignmentAxes() const { return m_alignmentAxes; }

    void setPartialCharges(const Vector& charges) { m_charges = charges; }

    void setPartialCharges(const std::vector<double>& charges)
    {
        m_charges = Vector::Zero(charges.size());
        for (int i = 0; i < charges.size(); ++i)
            m_charges(i) = charges[i];
    }

    inline Vector getPartialCharges() const { return m_charges; }

    void setDipole(const Position& dipole) { m_dipole = dipole; }

    inline Position getDipole() const { return m_dipole; }

    inline std::vector<std::pair<int, int>> Bonds() const { return m_bonds; }

    inline void addBorderPoint(const Position& point) { m_borders.push_back(point); }

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
    Position m_dipole;
    Geometry m_geometry;
    std::vector<int> m_atoms;
    Vector m_charges;
    std::vector<Position> m_borders;
    std::vector<int> m_connect_mass;
    Matrix m_HydrogenBondMap;
    Eigen::MatrixXd m_persistentImage, m_alignmentAxes;
    mutable std::vector<std::vector<int>> m_fragments;
    mutable std::map<int, int> m_fragment_assignment;

    mutable std::vector<double> m_mass_fragments;
    std::vector<std::pair<int, int>> m_bonds;

    mutable bool m_dirty = true;
    std::string m_name;
    double m_energy = 0, m_Ia = 0, m_Ib = 0, m_Ic = 0, m_mass = 0, m_hbond_cutoff = 3;
    mutable double m_scaling = 1.5;
};
} // namespace curcuma