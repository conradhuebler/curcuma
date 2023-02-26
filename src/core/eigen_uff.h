/*
 * <Simple UFF implementation for Cucuma. >
 * Copyright (C) 2022 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "hbonds.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#ifdef USE_D3
#include "src/core/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/dftd4interface.h"
#endif

#include "src/core/uff_par.h"
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

/*
class UFFThread : public CxxThread
{
public:
    UFFThread();
    ~UFFThread();

    virtual int execute() override
    {
        m_mddriver = new SimpleMD(m_controller, false);
        m_mddriver->setMolecule(m_molecule);
        m_mddriver->setBaseName("thread" + std::to_string(m_thread));
        m_mddriver->Initialise();
        m_mddriver->start();
        return 0;
    }
private:
};
*/

class eigenUFF {
public:
    eigenUFF(const json& controller);

    void UpdateGeometry(const double* coord);
    void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);

    void setMolecule(const std::vector<int>& atom_types, const std::vector<std::array<double, 3>>& geometry)
    {
        m_atom_types = atom_types;
        m_geometry = Eigen::MatrixXd::Zero(m_atom_types.size(), 3); // geometry;
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_geometry(i, 0) = geometry[i][0];
            m_geometry(i, 1) = geometry[i][1];
            m_geometry(i, 2) = geometry[i][2];
        }
    }

    void Initialise();

    double Calculate(bool gradient = true, bool verbose = false);

    //   std::vector<std::array<double, 3>> Gradient() const { return m_gradient; }
    void Gradient(double* gradient) const;

    void NumGrad(double* gradient);
    std::vector<std::array<double, 3>> NumGrad();

    void writeParameterFile(const std::string& file) const;
    void writeUFFFile(const std::string& file) const;

    json writeParameter() const;
    json writeUFF() const;

    void readParameterFile(const std::string& file);
    void readUFFFile(const std::string& file);

    void readParameter(const json& parameters);
    void readUFF(const json& parameters);

    void setInitialisation(bool initialised) { m_initialised = initialised; }

private:
    void AssignUffAtomTypes();
    void FindRings();

    double BondRestLength(int i, int j, double order);

    double DotProduct(const v& pos1, const v& pos2) const
    {
        return pos1[0] * pos2[0] + pos1[1] * pos2[1] + pos1[2] * pos2[2];
    }
    double Norm(const v& pos1) const
    {
        return sqrt(DotProduct(pos1, pos1));
    }
    double CalculateAngle(int atom1, int atom2, int atom3) const
    {
        v atom_0 = { m_geometry[atom1] };
        v atom_1 = { m_geometry[atom2] };
        v atom_2 = { m_geometry[atom3] };

        v vec_1 = { atom_0[0] - atom_1[0], atom_0[1] - atom_1[1], atom_0[2] - atom_1[2] };
        v vec_2 = { atom_0[0] - atom_2[0], atom_0[1] - atom_2[1], atom_0[2] - atom_2[2] };

        return acos(DotProduct(vec_1, vec_2) / (sqrt(DotProduct(vec_1, vec_1) * DotProduct(vec_2, vec_2)))) * 360 / 2.0 / pi;
    }

    v NormalVector(const v& i, const v& j, const v& k)
    {
        Eigen::Vector3d aa = Eigen::Vector3d{ i[0], i[1], i[2] };
        Eigen::Vector3d ab = Eigen::Vector3d{ j[0], j[1], j[2] };
        Eigen::Vector3d ac = Eigen::Vector3d{ k[0], k[1], k[2] };

        Eigen::Vector3d aba = ab - aa;
        Eigen::Vector3d abc = ab - ac;
        auto p = aba.cross(abc);
        return v{ p[0], p[1], p[2] };
    }

    Eigen::Vector3d NormalVector(int i, int j, int k)
    {
        Eigen::Vector3d aa = Eigen::Vector3d{ m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2) };
        Eigen::Vector3d ab = Eigen::Vector3d{ m_geometry(j, 0), m_geometry(j, 1), m_geometry(j, 2) };
        Eigen::Vector3d ac = Eigen::Vector3d{ m_geometry(k, 0), m_geometry(k, 1), m_geometry(k, 2) };

        Eigen::Vector3d aba = ab - aa;
        Eigen::Vector3d abc = ab - ac;
        return aba.cross(abc);
    }
    double CalculateBondStretching();

    double AngleBend(const std::array<double, 3>& i, const std::array<double, 3>& j, const std::array<double, 3>& k, double kijk, double C0, double C1, double C2);
    double CalculateAngleBending();

    double Dihedral(const v& i, const v& j, const v& k, const v& l, double V, double n, double phi0);
    double CalculateDihedral();

    double FullInversion(const int& i, const int& j, const int& k, const int& l, double d_forceConstant, double C0, double C1, double C2);
    double Inversion(const v& i, const v& j, const v& k, const v& l, double k_ijkl, double C0, double C1, double C2);
    double CalculateInversion();

    double NonBonds(const v& i, const v& j, double Dij, double xij);
    double CalculateNonBonds();
    double CalculateElectrostatic();

    double Distance(double x1, double x2, double y1, double y2, double z1, double z2) const;
    double DotProduct(double x1, double x2, double y1, double y2, double z1, double z2) const;

    double BondEnergy(double distance, double r, double k_ij, double D_ij = 70);

    std::vector<int> m_atom_types, m_uff_atom_types, m_coordination;
    std::vector<std::vector<int>> m_stored_bonds;
    std::vector<std::vector<int>> m_identified_rings;

    Matrix m_geometry, m_gradient;

    std::vector<UFFBond> m_uffbonds;
    std::vector<UFFAngle> m_uffangle;
    std::vector<UFFDihedral> m_uffdihedral;
    std::vector<UFFInversion> m_uffinversion;
    std::vector<UFFvdW> m_uffvdwaals;

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
    bool m_verbose = false;
    bool m_rings = false;
    double m_bond_scaling = 1, m_angle_scaling = 1, m_dihedral_scaling = 1, m_inversion_scaling = 1, m_vdw_scaling = 1, m_rep_scaling = 1, m_coulmob_scaling = 1;

    hbonds4::H4Correction m_h4correction;
#ifdef USE_D3
    DFTD3Interface* m_d3;
#endif
#ifdef USE_D4
    DFTD4Interface* m_d4;
#endif
};
