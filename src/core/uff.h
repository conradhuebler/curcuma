/*
 * <Simple UFF implementation for Cucuma. >
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "hbonds.h"

#include <vector>

#include <Eigen/Dense>

#include "src/core/global.h"

/* this will be the basic organics stuff the uff */
/* directly taken from
 * https://github.com/openbabel/openbabel/blob/master/data/UFF.prm
 * and
 * https://github.com/openbabel/openbabel/blob/master/src/forcefields/forcefielduff.cpp
 * and
 * https://github.com/openbabel/openbabel/blob/master/src/forcefields/forcefielduff.h
 *
 * originally published at: J. Am. Chem. Soc. (1992) 114(25) p. 10024-10035.
 */

const std::vector<double> Dummy = { 0.01, 180, 0.4, 5000, 12, 10.0, 0, 0, 9.66, 14.92, 0.7 };

const std::vector<double> H = { 0.354, 180, 2.886, 0.044, 12, 0.712, 0, 0, 4.528, 6.9452, 0.371 };

const std::vector<double> C3 = { 0.757, 109.47, 3.851, 0.105, 12.73, 1.912, 2.119, 2, 5.343, 5.063, 0.759 };
const std::vector<double> CR = { 0.729, 120, 3.851, 0.105, 12.73, 1.912, 0, 2, 5.343, 5.063, 0.759 };
const std::vector<double> C2 = { 0.732, 120, 3.851, 0.105, 12.73, 1.912, 0, 2, 5.343, 5.063, 0.759 };
const std::vector<double> C1 = { 0.706, 180, 3.851, 0.105, 12.73, 1.912, 0, 2, 5.343, 5.063, 0.759 };

const std::vector<double> N3 = { 0.7, 106.7, 3.66, 0.069, 13.407, 2.544, 0.45, 2, 6.899, 5.88, 0.715 };
const std::vector<double> NR = { 0.699, 120, 3.66, 0.069, 13.407, 2.544, 0, 2, 6.899, 5.88, 0.715 };
const std::vector<double> N2 = { 0.685, 111.2, 3.66, 0.069, 13.407, 2.544, 0, 2, 6.899, 5.88, 0.715 };
const std::vector<double> N1 = { 0.656, 180, 3.66, 0.069, 13.407, 2.544, 0, 2, 6.899, 5.88, 0.715 };

const std::vector<double> O3 = { 0.658, 104.51, 3.5, 0.06, 14.085, 2.3, 0.018, 2, 8.741, 6.682, 0.669 };
const std::vector<double> O3z = { 0.528, 146, 3.5, 0.06, 14.085, 2.3, 0.018, 2, 8.741, 6.682, 0.669 };
const std::vector<double> OR = { 0.68, 110, 3.5, 0.06, 14.085, 2.3, 0, 2, 8.741, 6.682, 0.669 };

const std::vector<double> O2 = { 0.634, 120, 3.5, 0.06, 14.085, 2.3, 0, 2, 8.741, 6.682, 0.669 };
const std::vector<double> O1 = { 0.639, 180, 3.5, 0.06, 14.085, 2.3, 0, 2, 8.741, 6.682, 0.669 };

const int cR = 0;
const int cTheta0 = 1;
const int cx = 2;
const int cD = 3;
const int cZeta = 4;
const int cZ = 5;
const int cV = 6;
const int cU = 7;
const int cXi = 8;
const int cHard = 9;
const int cRadius = 10;

typedef std::array<double, 3> v;

inline std::array<double, 3> AddVector(const v& x, const v& y)
{
    return std::array<double, 3>{ x[0] + y[0], x[1] + y[1], x[2] + y[2] };
}
inline std::array<double, 3> SubVector(const v& x, const v& y)
{
    return std::array<double, 3>{ x[0] - y[0], x[1] - y[1], x[2] - y[2] };
}
class TContainer {
public:
    TContainer() = default;
    bool insert(std::vector<int> vector)
    {
        std::vector<int> cache = vector;
        std::sort(cache.begin(), cache.end());
        if (std::find(m_sorted.begin(), m_sorted.end(), cache) != m_sorted.end())
            return false;
        m_storage.push_back(vector);
        m_sorted.push_back(cache);
        return true;
    }
    const std::vector<std::vector<int>>& Storage() const { return m_storage; }

private:
    std::vector<std::vector<int>> m_storage, m_sorted;
};

class UFF {
public:
    UFF();

    void UpdateGeometry(const double* coord);
    void setMolecule(const std::vector<int>& atom_types, const std::vector<std::array<double, 3>>& geometry)
    {
        m_atom_types = atom_types;
        m_geometry = geometry;
    }

    void Initialise();
    double Calculate(bool gradient = true);

    std::vector<std::array<double, 3>> Gradient() const { return m_gradient; }
    void Gradient(double* gradient) const;

    void NumGrad(double* gradient);

private:
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
        Eigen::Vector3d aa = Eigen::Vector3d{ m_geometry[i][0], m_geometry[i][1], m_geometry[i][2] };
        Eigen::Vector3d ab = Eigen::Vector3d{ m_geometry[j][0], m_geometry[j][1], m_geometry[j][2] };
        Eigen::Vector3d ac = Eigen::Vector3d{ m_geometry[k][0], m_geometry[k][1], m_geometry[k][2] };

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

    double BondEnergy(double distance, double r, double k_ij, double D_ij);

    std::vector<int> m_atom_types, m_uff_atom_types, m_coordination;
    std::vector<std::array<double, 3>> m_geometry, m_gradient;
    std::vector<std::vector<double>> m_parameter;
    std::vector<std::pair<int, int>> m_bonds, m_non_bonds;
    std::vector<std::array<int, 3>> m_bond_angle;
    std::vector<std::array<int, 4>> m_dihedrals, m_inversions;
    std::vector<std::vector<int>> m_bonds_2, m_non_bonds_2, m_angles_2, m_dihedrals_2, m_inversions_2;

    double m_scaling = 1.15;
    Matrix m_topo;
    bool m_CalculateGradient = true, m_initialised = false;
    double m_d = 1e-6;
    double m_au = 1;
    double h_e1 = 1, h_e2 = 1;
    double m_final_factor = 1;
    ;
};
