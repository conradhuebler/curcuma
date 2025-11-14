/*
 * <Simple UFF implementation for Cucuma. >
 * Copyright (C) 2022 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

/* H4 Correction taken from
 * https://www.rezacovi.cz/science/sw/h_bonds4.c
 * Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)
 *            http://dx.doi.org/10.1021/ct200751e
 *
 */

#include "src/core/global.h"

#include "src/core/hbonds.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#ifdef USE_D3
#include "src/core/energy_calculators/qm_methods/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/energy_calculators/qm_methods/dftd4interface.h"
#endif

#include "uff_par.h"
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

class UFFThread : public CxxThread {
public:
    UFFThread(int thread, int threads)
        : m_thread(thread)
        , m_threads(threads)
    {
        setAutoDelete(false);
        m_final_factor = 1 / 2625.15 * 4.19;
        // m_d = parameters["differential"].get<double>();
        m_d = 1e-7;
    }
    //~UFFThread();

    virtual int execute() override;

    void UpdateGeometry(const double* coord);
    void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);
    void UpdateGeometry(Matrix* geometry)
    {
        m_geometry = geometry;
        m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    }

    Matrix Gradient() const { return m_gradient; }

    void setMolecule(const std::vector<int>& atom_types, Matrix* geometry)
    {
        m_atom_types = atom_types;
        m_geometry = geometry;
        m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    }
    void readUFF(const json& parameters);

    inline double Energy() const { return m_energy; }
    inline double BondEnergy() const { return m_bond_energy; }
    inline double AngleEnergy() const { return m_angle_energy; }
    inline double DihedralEnergy() const { return m_dihedral_energy; }
    inline double InversionEnergy() const { return m_inversion_energy; }
    inline double VdWEnergy() const { return m_vdw_energy; }

    void AddBond(const UFFBond& bond)
    {
        m_uffbonds.push_back(bond);
    }
    void AddAngle(const UFFAngle& angle)
    {
        m_uffangle.push_back(angle);
    }
    void AddDihedral(const UFFDihedral& dihedral)
    {
        m_uffdihedral.push_back(dihedral);
    }
    void AddInversion(const UFFInversion& inversion)
    {
        m_uffinversion.push_back(inversion);
    }
    void AddvdW(const UFFvdW& vdw)
    {
        m_uffvdwaals.push_back(vdw);
    }

private:
    double CalculateAngle(int atom1, int atom2, int atom3) const
    {
        auto atom_0 = m_geometry->row(atom1);
        auto atom_1 = m_geometry->row(atom2);
        auto atom_2 = m_geometry->row(atom3);

        auto vec_1 = atom_0 - atom_1;
        auto vec_2 = atom_0 - atom_2;

        return acos(vec_1.dot(vec_2) / sqrt(vec_1.dot(vec_1)) * sqrt(vec_2.dot(vec_2))) * 360 / 2.0 / pi;
    }

    Eigen::Vector3d NormalVector(const Eigen::Vector3d& aa, const Eigen::Vector3d& ab, const Eigen::Vector3d& ac)
    {
        Eigen::Vector3d aba = ab - aa;
        Eigen::Vector3d abc = ab - ac;
        return aba.cross(abc);
    }

    Eigen::Vector3d NormalVector(int i, int j, int k)
    {
        auto aa = m_geometry->row(i);
        auto ab = m_geometry->row(j);
        auto ac = m_geometry->row(k);

        Eigen::Vector3d aba = ab - aa;
        Eigen::Vector3d abc = ab - ac;
        return aba.cross(abc);
    }

    double BondEnergy(double distance, double r, double k_ij, double D_ij = 70);
    double CalculateBondStretching();

    double AngleBend(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, double kijk, double C0, double C1, double C2);
    double CalculateAngleBending();

    double Dihedral(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double V, double n, double phi0);
    double CalculateDihedral();

    double FullInversion(const int& i, const int& j, const int& k, const int& l, double d_forceConstant, double C0, double C1, double C2);
    double Inversion(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double k_ijkl, double C0, double C1, double C2);
    double CalculateInversion();

    double NonBonds(const Eigen::Vector3d& i, const Eigen::Vector3d& j, double Dij, double xij);
    double CalculateNonBonds();
    double CalculateElectrostatic();

    double Distance(double x1, double x2, double y1, double y2, double z1, double z2) const;
    double DotProduct(double x1, double x2, double y1, double y2, double z1, double z2) const;

    inline Eigen::Vector3d Position(int pos) const { return m_geometry->row(pos); }
    std::vector<int> m_atom_types, m_uff_atom_types, m_coordination;
    std::vector<std::vector<int>> m_stored_bonds;
    std::vector<std::vector<int>> m_identified_rings;

    Matrix *m_geometry, m_gradient;

    std::vector<UFFBond> m_uffbonds;
    int m_uff_bond_start = 0, m_uff_bond_end = 0;

    std::vector<UFFAngle> m_uffangle;
    int m_uff_angle_start = 0, m_uff_angle_end = 0;

    std::vector<UFFDihedral> m_uffdihedral;
    int m_uff_dihedral_start = 0, m_uff_dihedral_end = 0;

    std::vector<UFFInversion> m_uffinversion;
    int m_uff_inv_start = 0, m_uff_inv_end = 0;

    std::vector<UFFvdW> m_uffvdwaals;
    int m_uff_vdw_start = 0, m_uff_vdw_end = 0;

    double m_scaling = 1.15;
    Matrix m_topo;
    bool m_CalculateGradient = true, m_initialised = false;
    std::string m_writeparam = "none", m_writeuff = "none";
    double m_d = 1e-6;
    double m_au = 1;
    double m_h4_scaling = 1, m_hh_scaling = 1;
    double m_final_factor = 1;
    double m_bond_force = 664.12;
    double m_angle_force = 664.12;
    double m_energy = 0.0;
    double m_bond_energy = 0.0, m_angle_energy = 0.0, m_dihedral_energy = 0.0, m_inversion_energy = 0.0, m_vdw_energy = 0.0, m_d4_energy = 0.0, m_d3_energy = 0.0, m_energy_h4 = 0.0, m_energy_hh = 0.0;

    bool m_use_d3 = false;
    bool m_use_d4 = false;
    bool m_rings = false;
    int m_calc_gradient = 0;

    double m_bond_scaling = 1, m_angle_scaling = 1, m_dihedral_scaling = 1, m_inversion_scaling = 1, m_vdw_scaling = 1, m_rep_scaling = 1, m_coulmob_scaling = 1;
    int m_thread = 0, m_threads = 0;
};

class eigenUFF {
public:
    eigenUFF(const json& controller);
    ~eigenUFF();

    void UpdateGeometry(const double* coord);
    void UpdateGeometry(const Matrix& geometry);

    void setMolecule(const std::vector<int>& atom_types, const Matrix& geometry)
    {
        m_atom_types = atom_types;
        m_geometry = geometry;
    }

    void Initialise(const std::vector<std::pair<int, int>>& formed_bonds = std::vector<std::pair<int, int>>());

    double Calculate(bool gradient = true, bool verbose = false);

    Matrix Gradient() const { return m_gradient; }
    void Gradient(double* gradient) const;

    void NumGrad(double* gradient);
    Eigen::MatrixXd NumGrad();

    void writeParameterFile(const std::string& file) const;
    void writeUFFFile(const std::string& file) const;

    // json writeParameter() const;
    json writeUFF() const;

    json Bonds() const;
    json Angles() const;
    json Dihedrals() const;
    json Inversions() const;
    json vdWs() const;

    void setBonds(const json& bonds);
    void setBonds(const TContainer& bonds, std::vector<std::set<int>>& ignored_vdw, TContainer& angels, TContainer& dihedrals, TContainer& inversions);
    void setBonds(const std::vector<std::pair<int, int>>& bonds);

    void setAngles(const json& angles);
    void setAngles(const TContainer& angles, const std::vector<std::set<int>>& ignored_vdw);

    void setDihedrals(const json& dihedrals);
    void setDihedrals(const TContainer& dihedrals);

    void setInversions(const json& inversions);
    void setInversions(const TContainer& inversions);

    void setvdWs(const json& vdws);
    void setvdWs(const std::vector<std::set<int>>& ignored_vdw);

    void readParameterFile(const std::string& file);
    void readUFFFile(const std::string& file);

    void readParameter(const json& parameters);
    void readUFF(const json& parameters);

    void setInitialisation(bool initialised) { m_initialised = initialised; }

    void AutoRanges();
    void setBondRanges(int start, int end)
    {
        m_uff_bond_start = start;
        m_uff_bond_end = end;
    }

    void seAngleRanges(int start, int end)
    {
        m_uff_angle_start = start;
        m_uff_angle_end = end;
    }

    void setDihedralRanges(int start, int end)
    {
        m_uff_dihedral_start = start;
        m_uff_dihedral_end = end;
    }

    void setInversionRanges(int start, int end)
    {
        m_uff_inv_start = start;
        m_uff_inv_end = end;
    }

    void setVDWRanges(int start, int end)
    {
        m_uff_vdw_start = start;
        m_uff_vdw_end = end;
    }

private:
    void AssignUffAtomTypes();

    double BondRestLength(int i, int j, double order);

    std::vector<int> m_atom_types, m_uff_atom_types, m_coordination;
    std::vector<std::vector<int>> m_stored_bonds;
    std::vector<std::vector<int>> m_identified_rings;

    Matrix m_geometry, m_gradient;

    std::vector<UFFBond> m_uffbonds;
    int m_uff_bond_start = 0, m_uff_bond_end = 0;

    std::vector<UFFAngle> m_uffangle;
    int m_uff_angle_start = 0, m_uff_angle_end = 0;

    std::vector<UFFDihedral> m_uffdihedral;
    int m_uff_dihedral_start = 0, m_uff_dihedral_end = 0;

    std::vector<UFFInversion> m_uffinversion;
    int m_uff_inv_start = 0, m_uff_inv_end = 0;

    std::vector<UFFvdW> m_uffvdwaals;
    int m_uff_vdw_start = 0, m_uff_vdw_end = 0;

    double m_scaling = 1.15;
    Matrix m_topo;
    bool m_CalculateGradient = true, m_initialised = false;
    std::string m_writeparam = "none", m_writeuff = "none";
    double m_d = 1e-6;
    double m_au = 1;
    double m_h4_scaling = 1, m_hh_scaling = 1;
    double m_final_factor = 1;
    double m_bond_force = 664.12;
    double m_angle_force = 664.12;

    bool m_use_d3 = false;
    bool m_use_d4 = false;
    bool m_rings = false;
    double m_bond_scaling = 1, m_angle_scaling = 1, m_dihedral_scaling = 1, m_inversion_scaling = 1, m_vdw_scaling = 1, m_rep_scaling = 1, m_coulmob_scaling = 1;
    bool m_numtorsion = false;
    int m_calc_gradient = 0;

    int m_threads = 1;
    hbonds4::H4Correction m_h4correction;
    std::vector<UFFThread*> m_stored_threads;
    CxxThreadPool* m_threadpool;

#ifdef USE_D3
    DFTD3Interface* m_d3;
#endif
#ifdef USE_D4
    DFTD4Interface* m_d4;
#endif
};
