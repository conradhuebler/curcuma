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

#include "src/core/global.h"
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
    // m_au = au;
}

void UFF::Initialise()
{
    m_uff_atom_types = std::vector<int>(m_atom_types.size(), 0);
    m_coordination = std::vector<int>(m_atom_types.size(), 0);

    m_bonds.clear();
    m_topo = Eigen::MatrixXd::Zero(m_atom_types.size(), m_atom_types.size());
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_gradient.push_back({ 0, 0, 0 });
        for (int j = i + 1; j < m_atom_types.size(); ++j) {
            double x_i = m_geometry[i][0] * m_au;
            double x_j = m_geometry[j][0] * m_au;

            double y_i = m_geometry[i][1] * m_au;
            double y_j = m_geometry[j][1] * m_au;

            double z_i = m_geometry[i][2] * m_au;
            double z_j = m_geometry[j][2] * m_au;

            double r_ij = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));

            if (r_ij <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling * m_au) {
                m_coordination[i]++;
                m_coordination[j]++;
                m_bonds.push_back(std::pair<int, int>(i, j));
                m_topo(i, j) = 1;
                m_topo(j, i) = 1;
                for (int k = j; k < m_atom_types.size(); ++k) {
                    if (i == k || j == k)
                        continue;

                    double x_k = m_geometry[k][0] * m_au;
                    double y_k = m_geometry[k][1] * m_au;
                    double z_k = m_geometry[k][2] * m_au;
                    double r_ik = sqrt((((x_i - x_k) * (x_i - x_k)) + ((y_i - y_k) * (y_i - y_k)) + ((z_i - z_k) * (z_i - z_k))));
                    if (r_ik <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[k]]) * m_scaling * m_au) {
                        m_bond_angle.push_back(std::array<int, 3>{ i, j, k });

                        for (int l = k; l < m_atom_types.size(); ++l) {
                            if (i == l || j == l || k == l)
                                continue;

                            double x_l = m_geometry[l][0] * m_au;
                            double y_l = m_geometry[l][1] * m_au;
                            double z_l = m_geometry[l][2] * m_au;
                            double r_kl = sqrt((((x_l - x_k) * (x_l - x_k)) + ((y_l - y_k) * (y_l - y_k)) + ((z_l - z_k) * (z_l - z_k))));
                            double r_jl = sqrt((((x_l - x_j) * (x_l - x_j)) + ((y_l - y_j) * (y_l - y_j)) + ((z_l - z_j) * (z_l - z_j))));
                            if (r_kl <= (Elements::CovalentRadius[m_atom_types[l]] + Elements::CovalentRadius[m_atom_types[k]]) * m_scaling * m_au) {
                                m_dihedrals.push_back(std::array<int, 4>{ j, i, k, l });
                            }
                            if (r_jl <= (Elements::CovalentRadius[m_atom_types[l]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling) {
                                m_dihedrals.push_back(std::array<int, 4>{ l, j, i, k });
                            }
                        }
                    }
                }
            }
        }
        if (m_atom_types[i] == 1) {
            m_uff_atom_types[i] = 1;
        } else if (m_atom_types[i] == 6) {
            if (m_coordination[i] == 4)
                m_uff_atom_types[i] = 9;
            else if (m_coordination[i] == 3)
                m_uff_atom_types[i] = 10;
            else // if (coordination == 2)
                m_uff_atom_types[i] = 12;
        } else if (m_atom_types[i] == 7) {
            if (m_coordination[i] == 3)
                m_uff_atom_types[i] = 13;
            else if (m_coordination[i] == 2)
                m_uff_atom_types[i] = 14;
            else // (coordination == 1)
                m_uff_atom_types[i] = 15;
        } else if (m_atom_types[i] == 8) {
            if (m_coordination[i] == 3)
                m_uff_atom_types[i] = 17;
            else if (m_coordination[i] == 2)
                m_uff_atom_types[i] = 19;
            else // (coordination == 1)
                m_uff_atom_types[i] = 21;
        }
    }
    m_initialised = true;
}

void UFF::UpdateGeometry(const double* coord)
{
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry[i][0] = coord[3 * i] * au;
        m_geometry[i][1] = coord[3 * i + 1] * au;
        m_geometry[i][2] = coord[3 * i + 2] * au;

        m_gradient[i] = { 0, 0, 0 };
    }
}

void UFF::Gradient(double* grad) const
{
    // std::cout << std::endl;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        grad[3 * i] = m_gradient[i][0];
        grad[3 * i + 1] = m_gradient[i][1];
        grad[3 * i + 2] = m_gradient[i][2];
        //  std::cout << m_gradient[i][0] << " " << m_gradient[i][1] << " " <<m_gradient[i][0] << std::endl;
    }
    // std::cout << std::endl;
}

void UFF::NumGrad(double* grad)
{
    double dx = 1e-6;
    bool g = m_CalculateGradient;
    m_CalculateGradient = false;
    double E1, E2;
    //  std::cout << std::endl;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry[i][j] += dx;
            E1 = Calculate();
            m_geometry[i][j] -= 2 * dx;
            E2 = Calculate();
            grad[3 * i + j] = (E1 - E2) / (2 * dx);
            m_geometry[i][j] += dx;
            // std::cout << E1 << " " << E2 << " " << (E1-E2)/(2*dx) << std::endl;
        }
        // std::cout << grad[3*i] << " " << grad[3*i + 1] << " " <<grad[3*i + 2] << std::endl;
        /*
        m_geometry[i][1] +=dx;
        E1 = Calculate();
        m_geometry[i][1] -= 2*dx;
        E2 = Calculate();
        grad[3*i + 1 ] = (E1-E2)/(2*dx);
        m_geometry[i][1] +=dx;

        m_geometry[i][2] +=dx;
        E1 = Calculate();
        m_geometry[i][2] -= 2*dx;
        E2 = Calculate();
        grad[3*i + 2 ] = (E1-E2)/(2*dx);
        m_geometry[i][2] +=dx;
        */
        //  std::cout << m_gradient[i][0] << " " << m_gradient[i][1] << " " <<m_gradient[i][0] << std::endl;
    }
    //  std::cout << std::endl;
    m_CalculateGradient = g;
}

double UFF::BondRestLength(int i, int j, double n)
{
    double lambda = 0.13332;
    double r_BO = -lambda * (m_parameter[m_uff_atom_types[i]][0] + m_parameter[m_uff_atom_types[j]][0]) * log(n);
    double r_EN; // = m_parameter[ m_uff_atom_types[i] ][0]*m_parameter[ m_uff_atom_types[j] ][0]*(sqrt(m_parameter[ m_uff_atom_types[i] ][7])-sqrt(m_parameter[ m_uff_atom_types[j] ][7]))*(sqrt(m_parameter[ m_uff_atom_types[i] ][7])-sqrt(m_parameter[ m_uff_atom_types[j] ][7]))/(m_parameter[ m_uff_atom_types[i] ][0]*m_parameter[ m_uff_atom_types[i] ][7]+m_parameter[ m_uff_atom_types[j] ][0]*m_parameter[ m_uff_atom_types[j] ][7]);
    return (m_parameter[m_uff_atom_types[i]][0] + m_parameter[m_uff_atom_types[j]][0] + r_BO + r_EN) * m_au;
}

double UFF::Calculate()
{
    double energy = 0.0;
    energy = (CalculateBondStretching() + CalculateAngleBending() + CalculateDihedral() + CalculateInversion() + CalculateNonBonds() + CalculateElectrostatic()) / 2625.15 * 4.19;

    return energy;
}

double UFF::Distance(double x1, double x2, double y1, double y2, double z1, double z2) const
{
    return sqrt((((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2))));
}

double UFF::DotProduct(double x1, double x2, double y1, double y2, double z1, double z2) const
{
    return x1 * x2 + y1 * y2 + z1 * z2;
}

double UFF::BondEnergy(double distance, double r, double kij, double D_ij)
{
    // double k_ij = k/(r*r*r);
    // std::cout << k_ij << std::endl;
    return 0.5 * kij * (distance - r) * (distance - r);
    double alpha = sqrt(kij / (2 * D_ij));
    double exp_ij = exp(-1 * alpha * (r - distance) - 1);
    return D_ij * (exp_ij * exp_ij);
}

double UFF::CalculateBondStretching()
{
    m_d = 1e-5;
    double factor = 1;
    double energy = 0.0;
    for (const std::pair<int, int>& bond : m_bonds) {
        const int i = bond.first;
        const int j = bond.second;
        double xi = m_geometry[i][0] * m_au;
        double xj = m_geometry[j][0] * m_au;

        double yi = m_geometry[i][1] * m_au;
        double yj = m_geometry[j][1] * m_au;

        double zi = m_geometry[i][2] * m_au;
        double zj = m_geometry[j][2] * m_au;
        double D_ij = 70;
        double rij = BondRestLength(i, j, 1);
        // if(m_atom_types[j] == 7 && m_atom_types[i] == 1 || m_atom_types[i] == 7 && m_atom_types[j] == 1 )
        // {
        //  std::cout << rij << " " << Distance(xi, xj, yi, yj, zi, zj) << " " << m_atom_types[i] << " (" << i << ")" << m_atom_types[j] << "(" <<j <<")"<< std::endl;
        //     BondRestLength(i,j,1);
        // }
        double kij = 664.12 * m_parameter[m_uff_atom_types[i]][5] * m_parameter[m_uff_atom_types[j]][5] / (rij * rij * rij);
        double benergy = BondEnergy(Distance(xi, xj, yi, yj, zi, zj), rij, kij, D_ij);
        energy += benergy;
        // std::cout << Distance(x_i, x_j, y_i, y_j, z_i, z_j) << " (" << r - Distance(x_i, x_j, y_i, y_j, z_i, z_j) << ") ";
        if (m_CalculateGradient) {
            double dE_p = BondEnergy(Distance(xi + m_d, xj, yi, yj, zi, zj), rij, kij, D_ij);
            double dE_m = BondEnergy(Distance(xi - m_d, xj, yi, yj, zi, zj), rij, kij, D_ij);
            //   std::cout << dE_p << " " << dE_m << " " << factor*(dE_p - dE_m)/(2.0*m_d) <<  " " <<  factor*(dE_p - dE_m) << " "<< (2.0*m_d)<< std::endl;
            /*
                        m_gradient[i][0] += factor*(dE_p - dE_m)/(2.0*m_d);
                        m_gradient[j][0] -= factor*(dE_p - dE_m)/(2.0*m_d);
                        m_gradient[i][1] += factor*(dE_p - dE_m)/(2.0*m_d);
                        m_gradient[j][1] -= factor*(dE_p - dE_m)/(2.0*m_d);
                        m_gradient[i][2] += factor*(dE_p - dE_m)/(2.0*m_d);
                        m_gradient[j][2] -= factor*(dE_p - dE_m)/(2.0*m_d);
            */

            m_gradient[i][0] += (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[j][0] -= (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[i][1] += (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[j][1] -= (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[i][2] += (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[j][2] -= (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));

            /*
                        dE_p = BondEnergy(Distance(x_i, x_j + m_d, y_i, y_j, z_i, z_j), r, k, D_ij);
                        dE_m = BondEnergy(Distance(x_i, x_j - m_d, y_i, y_j, z_i, z_j), r, k, D_ij);
                        m_gradient[j][0] += factor*(dE_p - dE_m)/(2.0*m_d);

                        dE_p = BondEnergy(Distance(x_i, x_j, y_i + m_d, y_j, z_i, z_j), r, k, D_ij);
                        dE_m = BondEnergy(Distance(x_i, x_j, y_i - m_d, y_j, z_i, z_j), r, k, D_ij);
                        m_gradient[i][1] += factor*(dE_p - dE_m)/(2.0*m_d);

                        dE_p = BondEnergy(Distance(x_i, x_j, y_i, y_j + m_d, z_i, z_j), r, k, D_ij);
                        dE_m = BondEnergy(Distance(x_i, x_j, y_i, y_j - m_d, z_i, z_j), r, k, D_ij);
                        m_gradient[j][1] += factor*(dE_p - dE_m)/(2.0*m_d);

                        dE_p = BondEnergy(Distance(x_i, x_j, y_i, y_j, z_i + m_d, z_j), r, k, D_ij);
                        dE_m = BondEnergy(Distance(x_i, x_j, y_i, y_j, z_i - m_d, z_j), r, k, D_ij);
                        m_gradient[i][2] += factor*(dE_p - dE_m)/(2.0*m_d);

                        dE_p = BondEnergy(Distance(x_i, x_j, y_i, y_j, z_i, z_j + m_d), r, k, D_ij);
                        dE_m = BondEnergy(Distance(x_i, x_j, y_i, y_j, z_i, z_j - m_d), r, k, D_ij);
                        m_gradient[j][2] += factor*(dE_p - dE_m)/(2.0*m_d);
            */
        }
        //    std::cout << energy << " ";
    }
    // std::cout << std::endl;
    return energy;
}

double UFF::CalculateAngleBending()
{
    double energy = 0.0;
    for (const auto& bond : m_bond_angle) {
        const int i = bond[0];
        const int j = bond[1];
        const int k = bond[2];
        /*
                double xi = m_geometry[i][0] * au;
                double xj = m_geometry[j][0]* au;
                double xk = m_geometry[k][0]* au;

                double yi = m_geometry[i][1]* au;
                double yj = m_geometry[j][1]* au;
                double yk = m_geometry[k][1]* au;

                double zi = m_geometry[i][2]* au;
                double zj = m_geometry[j][2]* au;
                double zk = m_geometry[k][2]* au;
        */

        double rij = BondRestLength(i, j, 1);
        double rjk = BondRestLength(i, k, 1);
        double cosTheta0 = cos(m_parameter[m_uff_atom_types[i]][1]);
        double rik = sqrt(rij * rij + rjk * rjk - 2. * rij * rjk * cosTheta0);
        double param = 664.12;
        double beta = 2.0 * param / (rij * rjk);
        double preFactor = beta * m_parameter[m_uff_atom_types[i]][5] * m_parameter[m_uff_atom_types[k]][5] / (rik * rik * rik * rik * rik);
        double rTerm = rij * rjk;
        double inner = 3.0 * rTerm * (1.0 - cosTheta0 * cosTheta0) - rik * rik * cosTheta0;
        double k_ijk = preFactor * rTerm * inner;
        double bondtype = 1;
        double theta = CalculateAngle(i, j, k);
        /*std::cout << m_coordination[i]  << " " << theta
                  << " " << m_parameter[ m_uff_atom_types[i] ][1]
                  << " " << m_parameter[ m_uff_atom_types[j] ][1]
                  << " " << m_parameter[ m_uff_atom_types[k] ][1]
                    << std::endl;
                    */
        energy += k_ijk / bondtype * bondtype * (1 - cos(bondtype * theta));
    }
    return energy;
}

double UFF::CalculateDihedral()
{
    double energy = 0.0;
    /*for(const auto &bond : m_dihedrals)
    {
        std::cout << bond[0] << " " << bond[1] << " " << bond[2] << " " << bond[3] << std::endl;
    }*/
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
