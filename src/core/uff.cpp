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

#include "src/tools/general.h"

#include <Eigen/Dense>

#include "uff.h"

UFF::UFF()
{
    m_parameter = {
        Dummy,
        H,
        Dummy, // D
        Dummy, // H_b
        Dummy, // He
        Dummy, // Li
        Dummy, // Be3
        Dummy, // B3
        Dummy, // B2
        C3,
        CR,
        C2,
        C1,
        N3,
        NR,
        N2,
        N1,
        O3,
        O3z,
        OR,
        O2,
        O1,
    };
}

void UFF::Initialise()
{
    m_uff_atom_types = std::vector<int>(m_atom_types.size(), 0);
    m_bonds.clear();
    m_topo = Eigen::MatrixXd::Zero(m_atom_types.size(), m_atom_types.size());
    for (int i = 0; i < m_atom_types.size(); ++i) {
        int coordination = 0;
        for (int j = 0; j < i; ++j) {
            double x_i = m_geometry[i][0];
            double x_j = m_geometry[j][0];

            double y_i = m_geometry[i][1];
            double y_j = m_geometry[j][1];

            double z_i = m_geometry[i][2];
            double z_j = m_geometry[j][2];

            double distance = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));

            coordination += distance <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling;
            if (coordination) {
                m_bonds.push_back(std::pair<int, int>(i, j));
                m_topo(i, j) = 1;
                m_topo(j, i) = 1;
            }
        }

        if (m_atom_types[i] == 1) {
            m_atom_types[i] = 1;
        } else if (m_atom_types[i] == 6) {
            if (coordination == 3)
                m_atom_types[i] == 9;
            else if (coordination == 2)
                m_atom_types[i] == 10;
            else if (coordination == 1)
                m_atom_types[i] == 12;
        } else if (m_atom_types[i] == 7) {
            if (coordination == 3)
                m_atom_types[i] == 13;
            else if (coordination == 2)
                m_atom_types[i] == 14;
            else if (coordination == 1)
                m_atom_types[i] == 16;
        } else if (m_atom_types[i] == 8) {
            if (coordination == 3)
                m_atom_types[i] == 17;
            else if (coordination == 2)
                m_atom_types[i] == 19;
            else if (coordination == 1)
                m_atom_types[i] == 21;
        }
    }
}

double UFF::Calculate()
{
    double energy = 0.0;
    energy += CalculateBondStretching() + CalculateAngleBending() + CalculateDihedral() + CalculateInversion() + CalculateNonBonds() + CalculateElectrostatic();

    return energy;
}

double UFF::CalculateBondStretching()
{
    double energy = 0.0;
    for (const std::pair<int, int>& bond : m_bonds) {
        double x_i = m_geometry[bond.first][0];
        double x_j = m_geometry[bond.second][0];

        double y_i = m_geometry[bond.first][1];
        double y_j = m_geometry[bond.second][1];

        double z_i = m_geometry[bond.first][2];
        double z_j = m_geometry[bond.second][2];

        double distance = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));
        double r = m_parameter[bond.first][0] + m_parameter[bond.second][0];
        double k_ij = 664.12 * m_parameter[bond.first][5] * m_parameter[bond.second][5] / (distance * distance * distance);
        double D_ij = 70;
        double alpha = sqrt(k_ij / (2 * D_ij));
        double exp_ij = exp(-1 * alpha * (r - distance) - 1);
        energy += D_ij * (exp_ij * exp_ij);
    }
    return energy;
}

double UFF::CalculateAngleBending()
{
    double energy = 0.0;

    return energy;
}

double UFF::CalculateDihedral()
{
    double energy = 0.0;

    return energy;
}

double UFF::CalculateInversion()
{
    double energy = 0.0;

    return energy;
}
double UFF::CalculateNonBonds()
{
    double energy = 0.0;

    return energy;
}
double UFF::CalculateElectrostatic()
{
    double energy = 0.0;

    return energy;
}
