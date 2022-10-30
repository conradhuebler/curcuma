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

/* H4 Correction taken from
 * https://www.rezacovi.cz/science/sw/h_bonds4.c
 * Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)
 *            http://dx.doi.org/10.1021/ct200751e
 *
 */

#include "hbonds.h"

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
    TContainer bonds, nonbonds, angels, dihedrals, inversions;
    for (int i = 0; i < m_atom_types.size(); ++i) {

        m_gradient.push_back({ 0, 0, 0 });
        for (int j = 0; j < m_atom_types.size(); ++j) {
            if (i == j)
                continue;
            double x_i = m_geometry[i][0] * m_au;
            double x_j = m_geometry[j][0] * m_au;

            double y_i = m_geometry[i][1] * m_au;
            double y_j = m_geometry[j][1] * m_au;

            double z_i = m_geometry[i][2] * m_au;
            double z_j = m_geometry[j][2] * m_au;

            double r_ij = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));

            if (r_ij <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling * m_au) {
                if (bonds.insert({ i, j })) {
                    m_coordination[i]++;
                    m_coordination[j]++;
                }
                m_bonds.push_back(std::pair<int, int>(i, j));

                m_topo(i, j) = 1;
                m_topo(j, i) = 1;
                for (int k = 0; k < m_atom_types.size(); ++k) {
                    if (i == k || j == k)
                        continue;

                    double x_k = m_geometry[k][0] * m_au;
                    double y_k = m_geometry[k][1] * m_au;
                    double z_k = m_geometry[k][2] * m_au;
                    double r_ik = sqrt((((x_i - x_k) * (x_i - x_k)) + ((y_i - y_k) * (y_i - y_k)) + ((z_i - z_k) * (z_i - z_k))));
                    if (r_ik <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[k]]) * m_scaling * m_au) {
                        m_bond_angle.push_back(std::array<int, 3>{ i, j, k });
                        angels.insert({ i, j, k });

                        for (int l = 0; l < m_atom_types.size(); ++l) {
                            if (i == l || j == l || k == l)
                                continue;

                            double x_l = m_geometry[l][0] * m_au;
                            double y_l = m_geometry[l][1] * m_au;
                            double z_l = m_geometry[l][2] * m_au;
                            double r_kl = sqrt((((x_l - x_k) * (x_l - x_k)) + ((y_l - y_k) * (y_l - y_k)) + ((z_l - z_k) * (z_l - z_k))));
                            double r_jl = sqrt((((x_l - x_j) * (x_l - x_j)) + ((y_l - y_j) * (y_l - y_j)) + ((z_l - z_j) * (z_l - z_j))));
                            double r_il = sqrt((((x_l - x_i) * (x_l - x_i)) + ((y_l - y_i) * (y_l - y_i)) + ((z_l - z_i) * (z_l - z_i))));
                            if (r_kl <= (Elements::CovalentRadius[m_atom_types[l]] + Elements::CovalentRadius[m_atom_types[k]]) * m_scaling * m_au) {
                                m_dihedrals.push_back(std::array<int, 4>{ j, i, k, l });
                                dihedrals.insert({ j, i, k, l });
                            }
                            if (r_jl <= (Elements::CovalentRadius[m_atom_types[l]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling) {
                                m_dihedrals.push_back(std::array<int, 4>{ l, j, i, k });
                                dihedrals.insert({ l, j, i, k });
                            }
                            if (r_il <= (Elements::CovalentRadius[m_atom_types[l]] + Elements::CovalentRadius[m_atom_types[i]]) * m_scaling) {
                                m_inversions.push_back(std::array<int, 4>{ i, j, k, l });
                                inversions.insert({ i, j, k, l });
                            }
                        }
                    }
                }
            } else {
                m_non_bonds.push_back(std::pair<int, int>(i, j));
                nonbonds.insert({ i, j });
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
    m_bonds_2 = bonds.Storage();
    m_non_bonds_2 = nonbonds.Storage();
    m_angles_2 = angels.Storage();
    m_dihedrals_2 = dihedrals.Storage();
    m_inversions_2 = inversions.Storage();

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
    for (int i = 0; i < m_atom_types.size(); ++i) {
        grad[3 * i] = m_gradient[i][0] + h_e1 * h_gradient[i].x + h_e2 * h_gradient2[i].x;
        grad[3 * i + 1] = m_gradient[i][1] + h_e1 * h_gradient[i].y + h_e2 * h_gradient2[i].y;
        grad[3 * i + 2] = m_gradient[i][2] + h_e1 * h_gradient[i].z + h_e2 * h_gradient2[i].z;
    }
}

void UFF::NumGrad(double* grad)
{
    double dx = 1e-5;
    bool g = m_CalculateGradient;
    m_CalculateGradient = false;
    double E1, E2;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry[i][j] += dx;
            E1 = Calculate();
            m_geometry[i][j] -= 2 * dx;
            E2 = Calculate();
            grad[3 * i + j] = (E1 - E2) / (2 * dx);
            m_geometry[i][j] += dx;
            if (j == 0)
                m_geometry[i][j] += h_e1 * h_gradient[i].x + h_e2 * h_gradient2[i].x;
            else if (j == 1)
                m_geometry[i][j] += h_e1 * h_gradient[i].y + h_e2 * h_gradient2[i].y;
            else if (j == 2)
                m_geometry[i][j] += h_e1 * h_gradient[i].z + h_e2 * h_gradient2[i].z;
        }
    }
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

    hbonds4::atom_t geometry[m_atom_types.size()];
    hbonds4::coord_t* gradient;
    hbonds4::coord_t* gradient2;

    for (int i = 0; i < m_atom_types.size(); ++i) {
        geometry[i].x = m_geometry[i][0] * m_au;
        geometry[i].y = m_geometry[i][1] * m_au;
        geometry[i].z = m_geometry[i][2] * m_au;
        geometry[i].e = m_atom_types[i];
    }
    hbonds4::gradient_allocate(m_atom_types.size(), &gradient); // Allocate memory for H4 gradient
    hbonds4::gradient_allocate(m_atom_types.size(), &gradient2); // Allocate memory for HH repulsion gradient

    energy = (CalculateBondStretching() + CalculateAngleBending() + CalculateDihedral() + CalculateInversion() + CalculateNonBonds() + CalculateElectrostatic()) / 2625.15 * 4.19;

    double energy_h4;
    // energy_h4 = hbonds4::energy_corr_h4(m_atom_types.size(), geometry, gradient);
    double energy_hh;
    // energy_hh = hbonds4::energy_corr_hh_rep(m_atom_types.size(), geometry, gradient2);
    energy += h_e1 * energy_h4 + h_e2 * energy_hh;
    delete h_gradient;
    delete h_gradient2;
    h_gradient = gradient;
    h_gradient2 = gradient2;
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
    for (const auto& bond : m_bonds_2) {
        const int i = bond[0];
        const int j = bond[1];
        double xi = m_geometry[i][0] * m_au;
        double xj = m_geometry[j][0] * m_au;

        double yi = m_geometry[i][1] * m_au;
        double yj = m_geometry[j][1] * m_au;

        double zi = m_geometry[i][2] * m_au;
        double zj = m_geometry[j][2] * m_au;
        double D_ij = 70;
        double rij = BondRestLength(i, j, 1);
        double kij = 664.12 * m_parameter[m_uff_atom_types[i]][5] * m_parameter[m_uff_atom_types[j]][5] / (rij * rij * rij);
        double benergy = BondEnergy(Distance(xi, xj, yi, yj, zi, zj), rij, kij, D_ij);
        energy += benergy;
        if (m_CalculateGradient) {
            double dE_p = BondEnergy(Distance(xi + m_d, xj, yi, yj, zi, zj), rij, kij, D_ij);
            double dE_m = BondEnergy(Distance(xi - m_d, xj, yi, yj, zi, zj), rij, kij, D_ij);
            m_gradient[i][0] += (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[j][0] -= (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[i][1] += (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[j][1] -= (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[i][2] += (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[j][2] -= (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
        }
    }
    return energy;
}

double UFF::CalculateAngleBending()
{
    double energy = 0.0;
    for (const auto& bond : m_angles_2) {
        const int i = bond[0];
        const int j = bond[1];
        const int k = bond[2];

        double f = pi / 180.0;
        double rij = BondRestLength(i, j, 1);
        double rjk = BondRestLength(i, k, 1);
        double Theta0 = m_parameter[m_uff_atom_types[i]][1];
        double cosTheta0 = cos(Theta0 * f);
        double rik = sqrt(rij * rij + rjk * rjk - 2. * rij * rjk * cosTheta0);
        double param = 664.12;
        double beta = 2.0 * param / (rij * rjk);
        double preFactor = beta * m_parameter[m_uff_atom_types[i]][5] * m_parameter[m_uff_atom_types[k]][5] / (rik * rik * rik * rik * rik);
        double rTerm = rij * rjk;
        double inner = 3.0 * rTerm * (1.0 - cosTheta0 * cosTheta0) - rik * rik * cosTheta0;
        double k_ijk = preFactor * rTerm * inner;
        double bondtype = 1;
        double theta = (CalculateAngle(i, j, k));
        double costheta = cos(theta * f);

        double C2 = 1 / (4 * sin(Theta0 * f) * sin(Theta0 * f));
        double C1 = -4 * C2 * cosTheta0;
        double C0 = C2 * (2 * cosTheta0 * cosTheta0 + 1);
        //   if(m_atom_types[i] == 7)
        //       energy += k_ijk*(C0+C1*cosTheta0+C2*cos(2*m_parameter[m_uff_atom_types[i]][1]));
        //   else
        energy += k_ijk * (C0 + C1 * costheta + C2 * cos(2 * theta * f));
        // energy += k_ijk / (bondtype * bondtype) * (1 - cos(bondtype * theta));
    }
    return energy;
}

double UFF::CalculateDihedral()
{
    double energy = 0.0;
    for (const auto& bond : m_dihedrals_2) {
        const int i = bond[0];
        const int j = bond[1];
        const int k = bond[2];
        const int l = bond[3];
        /*
            auto normal = [this](int a, int b, int c)->Eigen::Vector3d
            {
                Eigen::Vector3d aa = Eigen::Vector3d{m_geometry[a][0] , m_geometry[a][1] , m_geometry[a][2] };
                Eigen::Vector3d ab = Eigen::Vector3d{m_geometry[b][0] , m_geometry[b][1] , m_geometry[b][2] };
                Eigen::Vector3d ac = Eigen::Vector3d{m_geometry[c][0] , m_geometry[c][1] , m_geometry[c][2] };

                Eigen::Vector3d aba = ab-aa;
                Eigen::Vector3d abc = ab-ac;
                return aba.cross(abc);
            };*/
        /*
                Eigen::Vector3d ai = Eigen::Vector3d{m_geometry[i][0] * au, m_geometry[i][1] * au, m_geometry[i][2] * au};
                Eigen::Vector3d aj = Eigen::Vector3d{m_geometry[j][0] * au, m_geometry[j][1] * au, m_geometry[j][2] * au};
                Eigen::Vector3d ak = Eigen::Vector3d{m_geometry[k][0] * au, m_geometry[k][1] * au, m_geometry[k][2] * au};
                Eigen::Vector3d al = Eigen::Vector3d{m_geometry[l][0] * au, m_geometry[l][1] * au, m_geometry[l][2] * au};
        */
        Eigen::Vector3d nabc = NormalVector(i, j, k);
        Eigen::Vector3d nbcd = NormalVector(j, k, l);
        double central_bond_order = 1;

        double phi = acos(nabc.dot(nbcd) / (nabc.norm() * nbcd.norm())) * 360 / 2.0 / pi;
        int n = 3;
        double V = 5 * sqrt(m_parameter[m_uff_atom_types[j]][7] * m_parameter[m_uff_atom_types[k]][7]) * (1 + 4.18 * log(central_bond_order));
        double phi0 = 180;
        double f = pi / 180.0;

        energy += 1 / 2 * V * (1 - cos(n * phi0 * f) * cos(n * phi * f));
        // std::cout << i << " " << j << " " << k << " " << l  << " " << " " <<" "<< theta << " " << std::endl;
    }
    return energy;
}

double UFF::CalculateInversion()
{
    double energy = 0.0;
    // return energy;
    for (const auto& inversion : m_inversions_2) {
        const int i = inversion[0];
        if (m_coordination[i] != 3)
            continue;
        double f = pi / 180.0;

        const int j = inversion[1];
        const int k = inversion[2];
        const int l = inversion[3];
        Eigen::Vector3d ai = Eigen::Vector3d{ m_geometry[i][0], m_geometry[i][1], m_geometry[i][2] };
        Eigen::Vector3d al = Eigen::Vector3d{ m_geometry[l][0], m_geometry[l][1], m_geometry[l][2] };
        Eigen::Vector3d ail = ai - al;

        Eigen::Vector3d nbcd = NormalVector(i, j, k);

        double C0 = 0.0;
        double C1 = 0.0;
        double C2 = 0.0;
        double d_forceConstant = 0;
        if (6 <= m_atom_types[i] && m_atom_types[i] <= 8) {
            C0 = 1.0;
            C1 = -1.0;
            C2 = 0.0;
            d_forceConstant = 6;
            if (m_atom_types[j] == 8 || m_atom_types[k] == 8 || m_atom_types[l] == 8)
                d_forceConstant = 50;
        } else {
            double w0 = pi / 180.0;
            switch (m_atom_types[i]) {
            // if the central atom is phosphorous
            case 15:
                w0 *= 84.4339;
                break;

            // if the central atom is arsenic
            case 33:
                w0 *= 86.9735;
                break;

            // if the central atom is antimonium
            case 51:
                w0 *= 87.7047;
                break;

            // if the central atom is bismuth
            case 83:
                w0 *= 90.0;
                break;
            }
            C2 = 1.0;
            C1 = -4.0 * cos(w0 * f);
            C0 = -(C1 * cos(w0 * f) + C2 * cos(2.0 * w0 * f));
            d_forceConstant = 22.0 / (C0 + C1 + C2);
        }
        // res /= 3.0;
        double cosY = (nbcd.dot(ail) / (nbcd.norm() * ail.norm())); //* 360 / 2.0 / pi;
        double sinYSq = 1.0 - cosY * cosY;
        double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
        // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
        double cos2W = 2.0 * sinY * sinY - 1.0;
        double res = d_forceConstant * (C0 + C1 * sinY + C2 * cos2W);
        energy += res;
        // std::cout << i << " (" << m_coordination[i]<< ") " << j << " " << k << " " << l  << " " << " " <<phi << std::endl;
    }
    return energy;
}
double UFF::CalculateNonBonds()
{
    double energy = 0.0;

    for (const auto& bond : m_non_bonds_2) {
        const int i = bond[0];
        const int j = bond[1];
        /*
    for(int i = 0; i < m_atom_types.size(); ++i)
    {
        for(int j = i + 1; j < m_atom_types.size(); ++j)
        {*/
        double xi = m_geometry[i][0];
        double xj = m_geometry[j][0];

        double yi = m_geometry[i][1];
        double yj = m_geometry[j][1];

        double zi = m_geometry[i][2];
        double zj = m_geometry[j][2];

        double x = Distance(xi, xj, yi, yj, zi, zj) * m_au;
        double D_IJ = sqrt(m_parameter[m_uff_atom_types[i]][3] * m_parameter[m_uff_atom_types[j]][3]);
        double x_IJ = sqrt(m_parameter[m_uff_atom_types[i]][2] * m_parameter[m_uff_atom_types[j]][2]);
        // std::cout << x << " (" << D_IJ*(-2*pow((x_IJ/x),6) + pow((x_IJ/x),12)) << ") ";
        energy += D_IJ * (-2 * pow((x_IJ / x), 6) + pow((x_IJ / x), 12));
        // }
    }
    // std::cout << std::endl;
    return energy;
}
double UFF::CalculateElectrostatic()
{
    double energy = 0.0;

    return energy;
}
