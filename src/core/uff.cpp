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
    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = 1e-6;
    // m_au = au;
}

void UFF::Initialise()
{
    m_uff_atom_types = std::vector<int>(m_atom_types.size(), 0);
    m_coordination = std::vector<int>(m_atom_types.size(), 0);

    m_bonds.clear();
    /*
        hbonds4::atom_t geometry[m_atom_types.size()];
        hbonds4::coord_t gradient[m_atom_types.size()];
        hbonds4::coord_t gradient2[m_atom_types.size()];
    */
    m_topo = Eigen::MatrixXd::Zero(m_atom_types.size(), m_atom_types.size());
    TContainer bonds, nonbonds, angels, dihedrals, inversions;
    for (int i = 0; i < m_atom_types.size(); ++i) {

        m_gradient.push_back({ 0, 0, 0 });
        /*
        gradient[i].x = 0;
        gradient[i].y = 0;
        gradient[i].z = 0;

        gradient2[i].x = 0;
        gradient2[i].y = 0;
        gradient2[i].z = 0;
        */
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
    }
    for (int i = 0; i < m_atom_types.size(); ++i) {
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
    /*
    h_geometry = geometry;
    h_gradient = gradient;
    h_gradient2 = gradient2;
    */
    m_initialised = true;
}

void UFF::UpdateGeometry(const double* coord)
{
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry[i][0] = coord[3 * i] * au;
        m_geometry[i][1] = coord[3 * i + 1] * au;
        m_geometry[i][2] = coord[3 * i + 2] * au;
        /*
                h_geometry[i].x = coord[3 * i] * au;
                h_geometry[i].y = coord[3 * i + 1] * au;
                h_geometry[i].z = coord[3 * i + 2] * au;
                h_geometry[i].e = m_atom_types[i];
        */
        m_gradient[i] = { 0, 0, 0 };
        /*
                h_gradient[i].x = 0;
                h_gradient[i].y = 0;
                h_gradient[i].z = 0;

                h_gradient2[i].x = 0;
                h_gradient2[i].y = 0;
                h_gradient2[i].z = 0;
                */
    }
}

void UFF::Gradient(double* grad) const
{
    double factor = 1;
    //  std::cout << std::endl;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        grad[3 * i] = m_gradient[i][0] * factor /* + m_final_factor * h_gradient[i].x + m_final_factor * h_gradient2[i].x */;
        grad[3 * i + 1] = m_gradient[i][1] * factor /*+  m_final_factor * h_gradient[i].y + m_final_factor * h_gradient2[i].y*/;
        grad[3 * i + 2] = m_gradient[i][2] * factor /*+  m_final_factor * h_gradient[i].z + m_final_factor * h_gradient2[i].z*/;
        //      std::cout << grad[3 * i] << " " << grad[3 * i + 1] << " " << grad[3 * i + 2] << std::endl;
    }
    //  std::cout << std::endl;
}

void UFF::NumGrad(double* grad)
{
    // std::cout << std::endl;
    double dx = m_d;
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
            /*
            if (j == 0)
                m_geometry[i][j] += m_final_factor * h_gradient[i].x + m_final_factor * h_gradient2[i].x;
            else if (j == 1)
                m_geometry[i][j] += m_final_factor * h_gradient[i].y + m_final_factor * h_gradient2[i].y;
            else if (j == 2)
                m_geometry[i][j] += m_final_factor * h_gradient[i].z + m_final_factor * h_gradient2[i].z;
                */
        }
        //    std::cout << grad[3 * i] << " " << grad[3 * i + 1] << " " << grad[3 * i + 2] << std::endl;
    }
    //  std::cout << std::endl;

    m_CalculateGradient = g;
}

double UFF::BondRestLength(int i, int j, double n)
{
    const double cRi = m_parameter[m_uff_atom_types[i]][cR];
    const double cRj = m_parameter[m_uff_atom_types[j]][cR];
    const double cXii = m_parameter[m_uff_atom_types[i]][cXi];
    const double cXij = m_parameter[m_uff_atom_types[j]][cXi];

    double lambda = 0.13332;
    double r_BO = -lambda * (cRi + cRj) * log(n);
    double r_EN = cRi * cRj * (sqrt(cXii) - sqrt(cXij)) * (sqrt(cXii) - sqrt(cXij)) / (cRi * cXii + cRj * cXij);
    double r_0 = cRi + cRj;
    return (r_0 + r_BO + r_EN) * m_au;
}

double UFF::Calculate(bool grd)
{
    m_CalculateGradient = grd;
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

    energy = (CalculateBondStretching() + CalculateAngleBending() + CalculateDihedral() + CalculateInversion() + CalculateNonBonds() + CalculateElectrostatic()); // ;

    double energy_h4 = hbonds4::energy_corr_h4(m_atom_types.size(), geometry, gradient);
    double energy_hh = hbonds4::energy_corr_hh_rep(m_atom_types.size(), geometry, gradient2);
    energy += m_final_factor * h_e1 * energy_h4 + m_final_factor * h_e2 * energy_hh;

    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_gradient[i][0] += m_final_factor * h_e1 * gradient[i].x + m_final_factor * h_e2 * gradient2[i].x;
        m_gradient[i][1] += m_final_factor * h_e1 * gradient[i].y + m_final_factor * h_e2 * gradient2[i].y;
        m_gradient[i][2] += m_final_factor * h_e1 * gradient[i].z + m_final_factor * h_e2 * gradient2[i].z;
    }

    delete gradient;
    delete gradient2;

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
    double energy = (0.5 * kij * (distance - r) * (distance - r)) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;

    double alpha = sqrt(kij / (2 * D_ij));
    double exp_ij = exp(-1 * alpha * (r - distance) - 1);
    return D_ij * (exp_ij * exp_ij);
}

double UFF::CalculateBondStretching()
{
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
        const double cZi = m_parameter[m_uff_atom_types[i]][cZ];
        const double cZj = m_parameter[m_uff_atom_types[j]][cZ];
        double kij = 664.12 * cZi * cZj / (rij * rij * rij);
        double benergy = BondEnergy(Distance(xi, xj, yi, yj, zi, zj), rij, kij, D_ij);
        // std::cout << Distance(xi, xj, yi, yj, zi, zj) << " " << rij << std::endl;
        energy += benergy;
        if (m_CalculateGradient) {

            m_gradient[i][0] += (BondEnergy(Distance(xi + m_d, xj, yi, yj, zi, zj), rij, kij, D_ij) - BondEnergy(Distance(xi - m_d, xj, yi, yj, zi, zj), rij, kij, D_ij)) / (2 * m_d);
            m_gradient[i][1] += (BondEnergy(Distance(xi, xj, yi + m_d, yj, zi, zj), rij, kij, D_ij) - BondEnergy(Distance(xi, xj, yi - m_d, yj, zi, zj), rij, kij, D_ij)) / (2 * m_d);
            m_gradient[i][2] += (BondEnergy(Distance(xi, xj, yi, yj, zi + m_d, zj), rij, kij, D_ij) - BondEnergy(Distance(xi, xj, yi, yj, zi - m_d, zj), rij, kij, D_ij)) / (2 * m_d);

            m_gradient[j][0] += (BondEnergy(Distance(xi, xj + m_d, yi, yj, zi, zj), rij, kij, D_ij) - BondEnergy(Distance(xi, xj - m_d, yi, yj, zi, zj), rij, kij, D_ij)) / (2 * m_d);
            m_gradient[j][1] += (BondEnergy(Distance(xi, xj, yi, yj + m_d, zi, zj), rij, kij, D_ij) - BondEnergy(Distance(xi, xj, yi, yj - m_d, zi, zj), rij, kij, D_ij)) / (2 * m_d);
            m_gradient[j][2] += (BondEnergy(Distance(xi, xj, yi, yj, zi, zj + m_d), rij, kij, D_ij) - BondEnergy(Distance(xi, xj, yi, yj, zi, zj - m_d), rij, kij, D_ij)) / (2 * m_d);

            /*
            double diff = (1.0 * kij * (xi - xj) * (sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj)) - rij)) / sqrt((zi - zj) * (zi - zj) + (yi - yj) * (yi - yj) + (xi - xj) * (xi - xj));
            m_gradient[i][0] += diff;
            m_gradient[j][0] -= diff;

            m_gradient[i][1] += diff;
            m_gradient[j][1] -= diff;

            m_gradient[i][2] += diff;
            m_gradient[j][2] -= diff;
            */
        }
    }
    return energy;
}
double UFF::AngleBend(const std::array<double, 3>& i, const std::array<double, 3>& j, const std::array<double, 3>& k, double kijk, double C0, double C1, double C2)
{
    std::array<double, 3> vec_1 = { i[0] - j[0], i[1] - j[1], i[2] - j[2] };
    std::array<double, 3> vec_2 = { i[0] - k[0], i[1] - k[1], i[2] - k[2] };

    double costheta = (DotProduct(vec_1, vec_2) / (sqrt(DotProduct(vec_1, vec_1) * DotProduct(vec_2, vec_2))));
    double energy = (kijk * (C0 + C1 * costheta + C2 * (2 * costheta * costheta - 1))) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFF::CalculateAngleBending()
{
    double energy = 0.0;
    std::array<double, 3> dx = { m_d, 0, 0 };
    std::array<double, 3> dy = { 0, m_d, 0 };
    std::array<double, 3> dz = { 0, 0, m_d };
    for (const auto& bond : m_angles_2) {
        const int i = bond[0];
        const int j = bond[1];
        const int k = bond[2];

        std::array<double, 3> atom_i = { m_geometry[i] };
        std::array<double, 3> atom_j = { m_geometry[j] };
        std::array<double, 3> atom_k = { m_geometry[k] };

        double f = pi / 180.0;
        double rij = BondRestLength(i, j, 1);
        double rjk = BondRestLength(j, k, 1);
        double Theta0 = m_parameter[m_uff_atom_types[i]][cTheta0];
        double cosTheta0 = cos(Theta0 * f);
        double rik = sqrt(rij * rij + rjk * rjk - 2. * rij * rjk * cosTheta0);
        double param = 664.12;
        double beta = 2.0 * param / (rij * rjk);
        double preFactor = beta * m_parameter[m_uff_atom_types[i]][cZ] * m_parameter[m_uff_atom_types[k]][cZ] / (rik * rik * rik * rik * rik);
        double rTerm = rij * rjk;
        double inner = 3.0 * rTerm * (1.0 - cosTheta0 * cosTheta0) - rik * rik * cosTheta0;
        double k_ijk = preFactor * rTerm * inner;
        double bondtype = 1;
        double theta = CalculateAngle(i, j, k);
        double costheta = cos(theta * f);
        double C2 = 1 / (4 * sin(Theta0 * f) * sin(Theta0 * f));
        double C1 = -4 * C2 * cosTheta0;
        double C0 = C2 * (2 * cosTheta0 * cosTheta0 + 1);
        // energy += k_ijk * (C0 + C1 * costheta + C2 * (2*costheta*costheta -1));
        energy += AngleBend(atom_i, atom_j, atom_k, k_ijk, C0, C1, C2);
        // std::cout << k_ijk * (C0 + C1 * costheta + C2 * (2*costheta*costheta -1)) - AngleBend(atom_i,atom_j,atom_k, k_ijk, C0,C1,C2) << std::endl;
        if (m_CalculateGradient) {
            // double dE = -k_ijk * (C1*sin(theta*f) + 2.0 * C2*sin(2 * theta*f));

            m_gradient[i][0] += (AngleBend(AddVector(atom_i, dx), atom_j, atom_k, k_ijk, C0, C1, C2) - AngleBend(SubVector(atom_i, dx), atom_j, atom_k, k_ijk, C0, C1, C2)) / (2 * m_d);
            m_gradient[i][1] += (AngleBend(AddVector(atom_i, dy), atom_j, atom_k, k_ijk, C0, C1, C2) - AngleBend(SubVector(atom_i, dy), atom_j, atom_k, k_ijk, C0, C1, C2)) / (2 * m_d);
            m_gradient[i][2] += (AngleBend(AddVector(atom_i, dz), atom_j, atom_k, k_ijk, C0, C1, C2) - AngleBend(SubVector(atom_i, dz), atom_j, atom_k, k_ijk, C0, C1, C2)) / (2 * m_d);

            m_gradient[j][0] += (AngleBend(atom_i, AddVector(atom_j, dx), atom_k, k_ijk, C0, C1, C2) - AngleBend(atom_i, SubVector(atom_j, dx), atom_k, k_ijk, C0, C1, C2)) / (2 * m_d);
            m_gradient[j][1] += (AngleBend(atom_i, AddVector(atom_j, dy), atom_k, k_ijk, C0, C1, C2) - AngleBend(atom_i, SubVector(atom_j, dy), atom_k, k_ijk, C0, C1, C2)) / (2 * m_d);
            m_gradient[j][2] += (AngleBend(atom_i, AddVector(atom_j, dz), atom_k, k_ijk, C0, C1, C2) - AngleBend(atom_i, SubVector(atom_j, dz), atom_k, k_ijk, C0, C1, C2)) / (2 * m_d);

            m_gradient[k][0] += (AngleBend(atom_i, atom_j, AddVector(atom_k, dx), k_ijk, C0, C1, C2) - AngleBend(atom_i, atom_j, SubVector(atom_k, dx), k_ijk, C0, C1, C2)) / (2 * m_d);
            m_gradient[k][1] += (AngleBend(atom_i, atom_j, AddVector(atom_k, dy), k_ijk, C0, C1, C2) - AngleBend(atom_i, atom_j, SubVector(atom_k, dy), k_ijk, C0, C1, C2)) / (2 * m_d);
            m_gradient[k][2] += (AngleBend(atom_i, atom_j, AddVector(atom_k, dz), k_ijk, C0, C1, C2) - AngleBend(atom_i, atom_j, SubVector(atom_k, dz), k_ijk, C0, C1, C2)) / (2 * m_d);
        }
        // energy += k_ijk / (bondtype * bondtype) * (1 - cos(bondtype * theta * f));
    }
    return energy;
}

double UFF::Dihedral(const v& i, const v& j, const v& k, const v& l, double V, double n, double phi0)
{
    v nabc = NormalVector(i, j, k);
    v nbcd = NormalVector(j, k, l);
    double n_abc = Norm(nabc);
    double n_bcd = Norm(nbcd);
    double dotpr = DotProduct(nabc, nbcd);
    double phi = acos(dotpr / (n_abc * n_bcd)) * 360 / 2.0 / pi;
    double f = pi / 180.0;

    double energy = (1 / 2.0 * V * (1 - cos(n * phi0) * cos(n * phi * f))) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFF::CalculateDihedral()
{
    double energy = 0.0;
    std::array<double, 3> dx = { m_d, 0, 0 };
    std::array<double, 3> dy = { 0, m_d, 0 };
    std::array<double, 3> dz = { 0, 0, m_d };

    for (const auto& bond : m_dihedrals_2) {
        const int i = bond[0];
        const int j = bond[1];
        const int k = bond[2];
        const int l = bond[3];
        v atom_i = { m_geometry[i] };
        v atom_j = { m_geometry[j] };
        v atom_k = { m_geometry[k] };
        v atom_l = { m_geometry[l] };

        Eigen::Vector3d nabc = NormalVector(i, j, k);
        Eigen::Vector3d nbcd = NormalVector(j, k, l);
        double central_bond_order = 1;

        double phi = acos(nabc.dot(nbcd) / (nabc.norm() * nbcd.norm())) * 360 / 2.0 / pi;
        int n = 3;
        double V;
        double phi0 = 0;
        double f = pi / 180.0;

        if (m_coordination[j] == 4 && m_coordination[k] == 4) {
            V = sqrt(m_parameter[m_uff_atom_types[j]][cV] * m_parameter[m_uff_atom_types[k]][cV]);
            phi0 = 180 * f;
        } else {
            V = 5 * sqrt(m_parameter[m_uff_atom_types[j]][cU] * m_parameter[m_uff_atom_types[k]][cU]) * (1 + 4.18 * log(central_bond_order));
            phi0 = 90 * f;
        }

        // std::cout << 1 / 2 * V * (1 - cos(n * phi0) * cos(n * phi * f)) << std::endl;
        // energy += 1 / 2 * V * (1 - cos(n * phi0) * cos(n * phi * f));
        energy += Dihedral(atom_i, atom_j, atom_k, atom_l, V, n, phi0);
        if (m_CalculateGradient) {
            // double dE = -k_ijk * (C1*sin(theta*f) + 2.0 * C2*sin(2 * theta*f));

            m_gradient[i][0] += (Dihedral(AddVector(atom_i, dx), atom_j, atom_k, atom_l, V, n, phi0) - Dihedral(SubVector(atom_i, dx), atom_j, atom_k, atom_l, V, n, phi0)) / (2 * m_d);
            m_gradient[i][1] += (Dihedral(AddVector(atom_i, dy), atom_j, atom_k, atom_l, V, n, phi0) - Dihedral(SubVector(atom_i, dy), atom_j, atom_k, atom_l, V, n, phi0)) / (2 * m_d);
            m_gradient[i][2] += (Dihedral(AddVector(atom_i, dz), atom_j, atom_k, atom_l, V, n, phi0) - Dihedral(SubVector(atom_i, dz), atom_j, atom_k, atom_l, V, n, phi0)) / (2 * m_d);

            m_gradient[j][0] += (Dihedral(atom_i, AddVector(atom_j, dx), atom_k, atom_l, V, n, phi0) - Dihedral(atom_i, SubVector(atom_j, dx), atom_k, atom_l, V, n, phi0)) / (2 * m_d);
            m_gradient[j][1] += (Dihedral(atom_i, AddVector(atom_j, dy), atom_k, atom_l, V, n, phi0) - Dihedral(atom_i, SubVector(atom_j, dy), atom_k, atom_l, V, n, phi0)) / (2 * m_d);
            m_gradient[j][2] += (Dihedral(atom_i, AddVector(atom_j, dz), atom_k, atom_l, V, n, phi0) - Dihedral(atom_i, SubVector(atom_j, dz), atom_k, atom_l, V, n, phi0)) / (2 * m_d);

            m_gradient[k][0] += (Dihedral(atom_i, atom_j, AddVector(atom_k, dx), atom_l, V, n, phi0) - Dihedral(atom_i, atom_j, SubVector(atom_k, dx), atom_l, V, n, phi0)) / (2 * m_d);
            m_gradient[k][1] += (Dihedral(atom_i, atom_j, AddVector(atom_k, dy), atom_l, V, n, phi0) - Dihedral(atom_i, atom_j, SubVector(atom_k, dy), atom_l, V, n, phi0)) / (2 * m_d);
            m_gradient[k][2] += (Dihedral(atom_i, atom_j, AddVector(atom_k, dz), atom_l, V, n, phi0) - Dihedral(atom_i, atom_j, SubVector(atom_k, dz), atom_l, V, n, phi0)) / (2 * m_d);

            m_gradient[l][0] += (Dihedral(atom_i, atom_j, atom_k, AddVector(atom_l, dx), V, n, phi0) - Dihedral(atom_i, atom_j, atom_k, SubVector(atom_l, dx), V, n, phi0)) / (2 * m_d);
            m_gradient[l][1] += (Dihedral(atom_i, atom_j, atom_k, AddVector(atom_l, dy), V, n, phi0) - Dihedral(atom_i, atom_j, atom_k, SubVector(atom_l, dy), V, n, phi0)) / (2 * m_d);
            m_gradient[l][2] += (Dihedral(atom_i, atom_j, atom_k, AddVector(atom_l, dz), V, n, phi0) - Dihedral(atom_i, atom_j, atom_k, SubVector(atom_l, dz), V, n, phi0)) / (2 * m_d);
        }

        // std::cout << 1 / 2 * V * (1 - cos(n * phi0) * cos(n * phi * f)) - Dihedral(atom_i, atom_j, atom_k, atom_l, V, n, phi0) << std::endl;
        // std::cout << i << " " << j << " " << k << " " << l  << " " << " " <<" "<< theta << " " << std::endl;
    }
    return energy;
}
double UFF::Inversion(const v& i, const v& j, const v& k, const v& l, double k_ijkl, double C0, double C1, double C2)
{
    v ail = SubVector(i, l);
    v nbcd = NormalVector(i, j, k);

    double cosY = (DotProduct(nbcd, ail) / (Norm(nbcd) * Norm(ail))); //* 360 / 2.0 / pi;
    // std::cout << cosY* 360 / 2.0 / pi << " ";

    double sinYSq = 1.0 - cosY * cosY;
    double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
    // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
    double cos2W = 2.0 * sinY * sinY - 1.0;
    double energy = (k_ijkl * (C0 + C1 * sinY + C2 * cos2W)) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFF::FullInversion(const int& i, const int& j, const int& k, const int& l, double d_forceConstant, double C0, double C1, double C2)
{
    double energy;
    std::array<double, 3> dx = { m_d, 0, 0 };
    std::array<double, 3> dy = { 0, m_d, 0 };
    std::array<double, 3> dz = { 0, 0, m_d };
    v atom_i = { m_geometry[i] };
    v atom_j = { m_geometry[j] };
    v atom_k = { m_geometry[k] };
    v atom_l = { m_geometry[l] };
    energy += Inversion(atom_i, atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2);
    if (m_CalculateGradient) {
        // double dE = -k_ijk * (C1*sin(theta*f) + 2.0 * C2*sin(2 * theta*f));

        m_gradient[i][0] += (Inversion(AddVector(atom_i, dx), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dx), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[i][1] += (Inversion(AddVector(atom_i, dy), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dy), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[i][2] += (Inversion(AddVector(atom_i, dz), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dz), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);

        m_gradient[j][0] += (Inversion(atom_i, AddVector(atom_j, dx), atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dx), atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[j][1] += (Inversion(atom_i, AddVector(atom_j, dy), atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dy), atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[j][2] += (Inversion(atom_i, AddVector(atom_j, dz), atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dz), atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);

        m_gradient[k][0] += (Inversion(atom_i, atom_j, AddVector(atom_k, dx), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dx), atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[k][1] += (Inversion(atom_i, atom_j, AddVector(atom_k, dy), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dy), atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[k][2] += (Inversion(atom_i, atom_j, AddVector(atom_k, dz), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dz), atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);

        m_gradient[l][0] += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dx), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dx), d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[l][1] += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dy), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dy), d_forceConstant, C0, C1, C2)) / (2 * m_d);
        m_gradient[l][2] += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dz), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dz), d_forceConstant, C0, C1, C2)) / (2 * m_d);
    }
    return energy;
}
double UFF::CalculateInversion()
{
    double energy = 0.0;
    // return energy;
    std::array<double, 3> dx = { m_d, 0, 0 };
    std::array<double, 3> dy = { 0, m_d, 0 };
    std::array<double, 3> dz = { 0, 0, m_d };
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
        v atom_i = { m_geometry[i] };
        v atom_j = { m_geometry[j] };
        v atom_k = { m_geometry[k] };
        v atom_l = { m_geometry[l] };
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
        // std::cout << res - Inversion(atom_i, atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) << std::endl;
        energy += FullInversion(i, j, k, l, d_forceConstant, C0, C1, C2);
        // energy += FullInversion(i,j,l,k, d_forceConstant, C0, C1, C2);
        // energy += FullInversion(i,l,j,k, d_forceConstant, C0, C1, C2);

        /*
        energy += Inversion(atom_i, atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2);
        if(m_CalculateGradient)
        {
            //double dE = -k_ijk * (C1*sin(theta*f) + 2.0 * C2*sin(2 * theta*f));

            m_gradient[i][0] += (Inversion(AddVector(atom_i, dx), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dx), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[i][1] += (Inversion(AddVector(atom_i, dy), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dy), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[i][2] += (Inversion(AddVector(atom_i, dz), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dz), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);

            m_gradient[j][0] += (Inversion(atom_i, AddVector(atom_j, dx),  atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dx),  atom_k, atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[j][1] += (Inversion(atom_i, AddVector(atom_j, dy),  atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dy),  atom_k, atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[j][2] += (Inversion(atom_i, AddVector(atom_j, dz),  atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dz),  atom_k, atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);

            m_gradient[k][0] += (Inversion(atom_i, atom_j, AddVector(atom_k, dx), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dx), atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[k][1] += (Inversion(atom_i, atom_j, AddVector(atom_k, dy), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dy), atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[k][2] += (Inversion(atom_i, atom_j, AddVector(atom_k, dz), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dz), atom_l, d_forceConstant, C0, C1, C2))/(2*m_d);

            m_gradient[l][0] += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dx), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dx), d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[l][1] += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dy), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dy), d_forceConstant, C0, C1, C2))/(2*m_d);
            m_gradient[l][2] += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dz), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dz), d_forceConstant, C0, C1, C2))/(2*m_d);
        }
        */
        // std::cout << i << " (" << m_coordination[i]<< ") " << j << " " << k << " " << l  << " " << " " <<phi << std::endl;
    }
    //      std::cout << std::endl;
    return energy;
}

double UFF::NonBonds(const v& i, const v& j, double Dij, double xij)
{
    double r = Distance(i[0], j[0], i[1], j[1], i[2], j[2]) * m_au;
    double pow6 = pow((xij / r), 6);
    double energy = Dij * (-2 * pow6 + pow6 * pow6) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFF::CalculateNonBonds()
{
    double energy = 0.0;
    std::array<double, 3> dx = { m_d, 0, 0 };
    std::array<double, 3> dy = { 0, m_d, 0 };
    std::array<double, 3> dz = { 0, 0, m_d };
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

        v atom_i = { m_geometry[i] };
        v atom_j = { m_geometry[j] };
        double cDi = m_parameter[m_uff_atom_types[i]][cD];
        double cDj = m_parameter[m_uff_atom_types[j]][cD];
        double cxi = m_parameter[m_uff_atom_types[i]][cx];
        double cxj = m_parameter[m_uff_atom_types[j]][cx];
        double D_IJ = sqrt(cDi * cDj);

        double x_IJ = sqrt(cxi * cxj);

        // std::cout << x << " (" << D_IJ*(-2*pow((x_IJ/x),6) + pow((x_IJ/x),12)) << ") ";
        energy += NonBonds(atom_i, atom_j, D_IJ, x_IJ); // D_IJ * (-2 * pow((x_IJ / x), 6) + pow((x_IJ / x), 12));
        //     std::cout << D_IJ * (-2 * pow((x_IJ / x), 6) + pow((x_IJ / x), 12)) - NonBonds(atom_i, atom_j, D_IJ, x_IJ)<< std::endl;
        if (m_CalculateGradient) {
            // double dE = -k_ijk * (C1*sin(theta*f) + 2.0 * C2*sin(2 * theta*f));

            m_gradient[i][0] += (NonBonds(AddVector(atom_i, dx), atom_j, D_IJ, x_IJ) - NonBonds(SubVector(atom_i, dx), atom_j, D_IJ, x_IJ)) / (2 * m_d);
            m_gradient[i][1] += (NonBonds(AddVector(atom_i, dy), atom_j, D_IJ, x_IJ) - NonBonds(SubVector(atom_i, dy), atom_j, D_IJ, x_IJ)) / (2 * m_d);
            m_gradient[i][2] += (NonBonds(AddVector(atom_i, dz), atom_j, D_IJ, x_IJ) - NonBonds(SubVector(atom_i, dz), atom_j, D_IJ, x_IJ)) / (2 * m_d);

            m_gradient[j][0] += (NonBonds(atom_i, AddVector(atom_j, dx), D_IJ, x_IJ) - NonBonds(atom_i, SubVector(atom_j, dx), D_IJ, x_IJ)) / (2 * m_d);
            m_gradient[j][1] += (NonBonds(atom_i, AddVector(atom_j, dy), D_IJ, x_IJ) - NonBonds(atom_i, SubVector(atom_j, dy), D_IJ, x_IJ)) / (2 * m_d);
            m_gradient[j][2] += (NonBonds(atom_i, AddVector(atom_j, dz), D_IJ, x_IJ) - NonBonds(atom_i, SubVector(atom_j, dz), D_IJ, x_IJ)) / (2 * m_d);
        }
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
