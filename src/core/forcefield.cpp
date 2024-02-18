/*
 * < Generic force field class for curcuma . >
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

#include "src/core/forcefieldderivaties.h"
#include "src/core/qmdff_par.h"
#include "src/core/uff_par.h"

#include "forcefield.h"

ForceFieldThread::ForceFieldThread()
{
}

double ForceFieldThread::LJBondStretching()
{
    double factor = m_final_factor * m_bond_scaling;
    double energy = 0.0;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };

    for (int index = 0; index < m_bonds.size(); ++index) {
        const auto& bond = m_bonds[index];

        Vector i = Position(bond.i);
        Vector j = Position(bond.j);

        // Matrix derivate;
        energy += StretchEnergy(i - j, bond.r0_ij, bond.fc, bond.exponent) * factor;
        if (m_calculate_gradient) {
            /*if (m_calc_gradient == 0) {
                double diff = 0;
                m_gradient.row(a) += diff * derivate.row(0);
                m_gradient.row(b) += diff * derivate.row(1);

            } else if (m_calc_gradient == 1) {*/
            double d_x = (StretchEnergy(i + dx - j, bond.r0_ij, bond.fc, bond.exponent) - StretchEnergy(i - dx - j, bond.r0_ij, bond.fc, bond.exponent)) / (2 * m_d) * factor;
            double d_y = (StretchEnergy(i + dy - j, bond.r0_ij, bond.fc, bond.exponent) - StretchEnergy(i - dy - j, bond.r0_ij, bond.fc, bond.exponent)) / (2 * m_d) * factor;
            double d_z = (StretchEnergy(i + dz - j, bond.r0_ij, bond.fc, bond.exponent) - StretchEnergy(i - dz - j, bond.r0_ij, bond.fc, bond.exponent)) / (2 * m_d) * factor;
            m_gradient(bond.a, 0) += d_x;
            m_gradient(bond.a, 1) += d_y;
            m_gradient(bond.a, 2) += d_z;

            m_gradient(bond.b, 0) -= d_x;
            m_gradient(bond.b, 1) -= d_y;
            m_gradient(bond.b, 2) -= d_x;
            //}
        }
    }

    return energy;
}

double ForceFieldThread::HarmonicBondStretching()
{
    double factor = m_final_factor * m_bond_scaling;
    double energy = 0.0;

    for (int index = 0; index < m_bonds.size(); ++index) {
        const auto& bond = m_bonds[index];
        const int i = bond.i;
        const int j = bond.j;

        Vector x = Position(i);
        Vector y = Position(j);
        Matrix derivate;
        double r0 = BondStretching(x, y, derivate, m_calculate_gradient);

        energy += (0.5 * bond.fc * (r0 - bond.r0_ij) * (r0 - bond.r0_ij)) * factor;
        if (m_calculate_gradient) {
            double diff = (bond.fc) * (r0 - bond.r0_ij) * factor;
            m_gradient.row(i) += diff * derivate.row(0);
            m_gradient.row(j) += diff * derivate.row(1);
        }
    }

    return energy;
}

double QMDFFThread::CalculateAngleBending()
{
    double threshold = 1e-2;
    double energy = 0.0;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_qmdffangle.size(); ++index) {
        const auto& angle = m_qmdffangle[index];
        const int a = angle.a;
        const int b = angle.b;
        const int c = angle.c;

        auto atom_a = Position(a);
        auto atom_b = Position(b);
        auto atom_c = Position(c);
        Matrix derivate;
        double costheta = AngleBending(atom_a, atom_b, atom_c, derivate, m_CalculateGradient);
        std::function<double(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double, double, double, double)> angle_function;
        if (std::abs(costheta - pi) < threshold)

            angle_function = [this](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, double thetae, double kabc, double reAB, double reAC) -> double {
                double val = this->LinearAngleBend(a, b, c, thetae, kabc, reAB, reAC);
                if (std::isnan(val))
                    return 0;
                else
                    return val;
            };
        else
            angle_function = [this](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, double thetae, double kabc, double reAB, double reAC) -> double {
                double val = this->AngleBend(a, b, c, thetae, kabc, reAB, reAC);
                if (std::isnan(val))
                    return 0;
                else
                    return val;
            };

        double e = angle_function(atom_a, atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC); //(angle.kijk * (angle.C0 + angle.C1 * costheta + angle.C2 * (2 * costheta * costheta - 1))) * m_final_factor * m_angle_scaling;

        if (m_CalculateGradient) {
            if (m_calc_gradient == 0) {
                /*
                double sintheta = sin(acos(costheta));
                double dEdtheta = -angle.kijk * sintheta * (angle.C1 + 4 * angle.C2 * costheta) * m_final_factor * m_angle_scaling;
                m_gradient.row(i) += dEdtheta * derivate.row(0);
                m_gradient.row(j) += dEdtheta * derivate.row(1);
                m_gradient.row(k) += dEdtheta * derivate.row(2);
                */
            } else if (m_calc_gradient == 1) {
                m_gradient(a, 0) += (angle_function(AddVector(atom_a, dx), atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(SubVector(atom_a, dx), atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
                m_gradient(a, 1) += (angle_function(AddVector(atom_a, dy), atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(SubVector(atom_a, dy), atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
                m_gradient(a, 2) += (angle_function(AddVector(atom_a, dz), atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(SubVector(atom_a, dz), atom_b, atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);

                m_gradient(b, 0) += (angle_function(atom_a, AddVector(atom_b, dx), atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(atom_a, SubVector(atom_b, dx), atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
                m_gradient(b, 1) += (angle_function(atom_a, AddVector(atom_b, dy), atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(atom_a, SubVector(atom_b, dy), atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
                m_gradient(b, 2) += (angle_function(atom_a, AddVector(atom_b, dz), atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(atom_a, SubVector(atom_b, dz), atom_c, angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);

                m_gradient(c, 0) += (angle_function(atom_a, atom_b, AddVector(atom_c, dx), angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(atom_a, atom_b, SubVector(atom_c, dx), angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
                m_gradient(c, 1) += (angle_function(atom_a, atom_b, AddVector(atom_c, dy), angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(atom_a, atom_b, SubVector(atom_c, dy), angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
                m_gradient(c, 2) += (angle_function(atom_a, atom_b, AddVector(atom_c, dz), angle.thetae, angle.kabc, angle.reAB, angle.reAC) - angle_function(atom_a, atom_b, SubVector(atom_c, dz), angle.thetae, angle.kabc, angle.reAB, angle.reAC)) / (2 * m_d);
            }
        }
    }
    return energy;
}

ForceField::ForceField(const json& controller)
{
}

void ForceField::UpdateGeometry(const Matrix& geometry)
{
    m_geometry = geometry;
}

void ForceField::UpdateGeometry(const double* coord)
{
#pragma message("replace with raw data")
    for (int i = 0; i < m_natoms; ++i) {
        m_geometry(i, 0) = coord[3 * i + 0];
        m_geometry(i, 1) = coord[3 * i + 1];
        m_geometry(i, 2) = coord[3 * i + 2];
    }
}

void ForceField::UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);
{
#pragma message("replace with raw data")
    for (int i = 0; i < m_natoms; ++i) {
        m_geometry(i, 0) = geometry[i][0];
        m_geometry(i, 1) = geometry[i][1];
        m_geometry(i, 2) = geometry[i][2];
    }
}

void ForceField::setParameter(const json& parameter)
{
}
