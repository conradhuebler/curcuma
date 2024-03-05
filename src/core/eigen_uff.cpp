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

#include "src/core/global.h"
#include "src/tools/general.h"

/* H4 Correction taken from
 * https://www.rezacovi.cz/science/sw/h_bonds4.c
 * Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)
 *            http://dx.doi.org/10.1021/ct200751e
 *
 */

#include "hbonds.h"
#ifdef USE_D4
#include "src/core/dftd4interface.h"
#endif

#include <Eigen/Dense>

#include "src/core/forcefieldderivaties.h"
#include "src/core/topology.h"

#include "eigen_uff.h"

#include "json.hpp"
using json = nlohmann::json;

int UFFThread::execute()
{
    //    m_CalculateGradient = grd;
    m_d4_energy = 0;
    m_d3_energy = 0;
    m_bond_energy = CalculateBondStretching();
    m_angle_energy = CalculateAngleBending();
    m_dihedral_energy = CalculateDihedral();
    m_inversion_energy = CalculateInversion();
    m_vdw_energy = CalculateNonBonds();
    /* + CalculateElectrostatic(); */
    m_energy = m_bond_energy + m_angle_energy + m_dihedral_energy + m_inversion_energy + m_vdw_energy;
    return 0;
}

void UFFThread::readUFF(const json& parameter)
{
    //  json parameter = MergeJson(UFFParameterJson, parameters);

    m_d = parameter["differential"].get<double>();

    m_bond_scaling = parameter["bond_scaling"].get<double>();
    m_angle_scaling = parameter["angle_scaling"].get<double>();
    m_dihedral_scaling = parameter["dihedral_scaling"].get<double>();
    m_inversion_scaling = parameter["inversion_scaling"].get<double>();
    m_vdw_scaling = parameter["vdw_scaling"].get<double>();
    m_rep_scaling = parameter["rep_scaling"].get<double>();

    m_coulmob_scaling = parameter["coulomb_scaling"].get<double>();

    m_bond_force = parameter["bond_force"].get<double>();
    m_angle_force = parameter["angle_force"].get<double>();
    m_calc_gradient = parameter["gradient"].get<int>();
}

double UFFThread::Distance(double x1, double x2, double y1, double y2, double z1, double z2) const
{
    return sqrt((((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2))));
}

double UFFThread::BondEnergy(double distance, double r, double kij, double D_ij)
{
    double energy = (0.5 * kij * (distance - r) * (distance - r)) * m_final_factor * m_bond_scaling;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFFThread::CalculateBondStretching()
{
    double factor = 1;
    double energy = 0.0;

    for (int index = 0; index < m_uffbonds.size(); ++index) {
        const auto& bond = m_uffbonds[index];
        const int i = bond.i;
        const int j = bond.j;

        Vector x = Position(i);
        Vector y = Position(j);
        Matrix derivate;
        double r0 = BondStretching(x, y, derivate, m_CalculateGradient);

        energy += (0.5 * bond.kij * (r0 - bond.r0) * (r0 - bond.r0)) * m_final_factor * m_bond_scaling;
        if (m_CalculateGradient) {
            if (m_calc_gradient == 0) {
                double diff = (bond.kij) * (r0 - bond.r0) * m_final_factor * m_bond_scaling;
                m_gradient.row(i) += diff * derivate.row(0);
                m_gradient.row(j) += diff * derivate.row(1);

            } else if (m_calc_gradient == 1) {
                double xi = (*m_geometry)(i, 0) * m_au;
                double xj = (*m_geometry)(j, 0) * m_au;

                double yi = (*m_geometry)(i, 1) * m_au;
                double yj = (*m_geometry)(j, 1) * m_au;

                double zi = (*m_geometry)(i, 2) * m_au;
                double zj = (*m_geometry)(j, 2) * m_au;
                m_gradient(i, 0) += (BondEnergy(Distance(xi + m_d, xj, yi, yj, zi, zj), bond.r0, bond.kij) - BondEnergy(Distance(xi - m_d, xj, yi, yj, zi, zj), bond.r0, bond.kij)) / (2 * m_d);
                m_gradient(i, 1) += (BondEnergy(Distance(xi, xj, yi + m_d, yj, zi, zj), bond.r0, bond.kij) - BondEnergy(Distance(xi, xj, yi - m_d, yj, zi, zj), bond.r0, bond.kij)) / (2 * m_d);
                m_gradient(i, 2) += (BondEnergy(Distance(xi, xj, yi, yj, zi + m_d, zj), bond.r0, bond.kij) - BondEnergy(Distance(xi, xj, yi, yj, zi - m_d, zj), bond.r0, bond.kij)) / (2 * m_d);

                m_gradient(j, 0) += (BondEnergy(Distance(xi, xj + m_d, yi, yj, zi, zj), bond.r0, bond.kij) - BondEnergy(Distance(xi, xj - m_d, yi, yj, zi, zj), bond.r0, bond.kij)) / (2 * m_d);
                m_gradient(j, 1) += (BondEnergy(Distance(xi, xj, yi, yj + m_d, zi, zj), bond.r0, bond.kij) - BondEnergy(Distance(xi, xj, yi, yj - m_d, zi, zj), bond.r0, bond.kij)) / (2 * m_d);
                m_gradient(j, 2) += (BondEnergy(Distance(xi, xj, yi, yj, zi, zj + m_d), bond.r0, bond.kij) - BondEnergy(Distance(xi, xj, yi, yj, zi, zj - m_d), bond.r0, bond.kij)) / (2 * m_d);
            }
        }
    }

    return energy;
}

double UFFThread::AngleBend(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, double kijk, double C0, double C1, double C2)
{
    Eigen::Vector3d vec_1 = j - i;
    Eigen::Vector3d vec_2 = j - k;

    double costheta = (vec_1.dot(vec_2) / (sqrt(vec_1.dot(vec_1) * vec_2.dot(vec_2))));
    double energy = (kijk * (C0 + C1 * costheta + C2 * (2 * costheta * costheta - 1))) * m_final_factor * m_angle_scaling;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFFThread::CalculateAngleBending()
{
    double energy = 0.0;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_uffangle.size(); ++index) {
        const auto& angle = m_uffangle[index];
        const int i = angle.i;
        const int j = angle.j;
        const int k = angle.k;

        auto atom_i = Position(i);
        auto atom_j = Position(j);
        auto atom_k = Position(k);
        Matrix derivate;
        double costheta = AngleBending(atom_i, atom_j, atom_k, derivate, m_CalculateGradient);
        double e = (angle.kijk * (angle.C0 + angle.C1 * costheta + angle.C2 * (2 * costheta * costheta - 1))) * m_final_factor * m_angle_scaling;
        energy += e;
        if (m_CalculateGradient) {
            if (m_calc_gradient == 0) {
                double sintheta = sin(acos(costheta));
                double dEdtheta = -angle.kijk * sintheta * (angle.C1 + 4 * angle.C2 * costheta) * m_final_factor * m_angle_scaling;
                m_gradient.row(i) += dEdtheta * derivate.row(0);
                m_gradient.row(j) += dEdtheta * derivate.row(1);
                m_gradient.row(k) += dEdtheta * derivate.row(2);

            } else if (m_calc_gradient == 1) {
                m_gradient(i, 0) += (AngleBend(AddVector(atom_i, dx), atom_j, atom_k, angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(SubVector(atom_i, dx), atom_j, atom_k, angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
                m_gradient(i, 1) += (AngleBend(AddVector(atom_i, dy), atom_j, atom_k, angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(SubVector(atom_i, dy), atom_j, atom_k, angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
                m_gradient(i, 2) += (AngleBend(AddVector(atom_i, dz), atom_j, atom_k, angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(SubVector(atom_i, dz), atom_j, atom_k, angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);

                m_gradient(j, 0) += (AngleBend(atom_i, AddVector(atom_j, dx), atom_k, angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(atom_i, SubVector(atom_j, dx), atom_k, angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
                m_gradient(j, 1) += (AngleBend(atom_i, AddVector(atom_j, dy), atom_k, angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(atom_i, SubVector(atom_j, dy), atom_k, angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
                m_gradient(j, 2) += (AngleBend(atom_i, AddVector(atom_j, dz), atom_k, angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(atom_i, SubVector(atom_j, dz), atom_k, angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);

                m_gradient(k, 0) += (AngleBend(atom_i, atom_j, AddVector(atom_k, dx), angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(atom_i, atom_j, SubVector(atom_k, dx), angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
                m_gradient(k, 1) += (AngleBend(atom_i, atom_j, AddVector(atom_k, dy), angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(atom_i, atom_j, SubVector(atom_k, dy), angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
                m_gradient(k, 2) += (AngleBend(atom_i, atom_j, AddVector(atom_k, dz), angle.kijk, angle.C0, angle.C1, angle.C2) - AngleBend(atom_i, atom_j, SubVector(atom_k, dz), angle.kijk, angle.C0, angle.C1, angle.C2)) / (2 * m_d);
            }
        }
    }
    return energy;
}

double UFFThread::Dihedral(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double V, double n, double phi0)
{
    Eigen::Vector3d nabc = NormalVector(i, j, k);
    Eigen::Vector3d nbcd = NormalVector(j, k, l);
    double n_abc = (nabc).norm();
    double n_bcd = (nbcd).norm();
    double dotpr = nabc.dot(nbcd);
    Eigen::Vector3d rji = j - i;
    double sign = (-1 * rji).dot(nbcd) < 0 ? -1 : 1;
    double phi = pi + sign * acos(dotpr / (n_abc * n_bcd));
    double energy = (1 / 2.0 * V * (1 - cos(n * phi0) * cos(n * phi))) * m_final_factor * m_dihedral_scaling;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFFThread::CalculateDihedral()
{
    double energy = 0.0;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_uffdihedral.size(); ++index) {
        const auto& dihedral = m_uffdihedral[index];
        const int i = dihedral.i;
        const int j = dihedral.j;
        const int k = dihedral.k;
        const int l = dihedral.l;
        Eigen::Vector3d atom_i = Position(i);
        Eigen::Vector3d atom_j = Position(j);
        Eigen::Vector3d atom_k = Position(k);
        Eigen::Vector3d atom_l = Position(l);
        Eigen::Vector3d eins = { 1.0, 1.0, 1.0 };
        Eigen::Vector3d nabc = NormalVector(atom_i, atom_j, atom_k);
        Eigen::Vector3d nbcd = NormalVector(atom_j, atom_k, atom_l);

        double n_abc = (nabc).norm();
        double n_bcd = (nbcd).norm();
        double dotpr = nabc.dot(nbcd);
        Eigen::Vector3d m = nabc;
        Eigen::Vector3d n = nbcd;
        Eigen::Vector3d rji = atom_j - atom_i;

        double sign = (-1 * rji).dot(n) < 0 ? -1 : 1;
        double phi = pi + sign * acos(dotpr / (n_abc * n_bcd));
        double sinphi = sin(phi);
        double e = Dihedral(atom_i, atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0);
        if (isnan(e))
            continue;
        energy += e;
        if (m_CalculateGradient) {
            if (m_calc_gradient == 0) {
                Eigen::Vector3d rkj = atom_k - atom_j;
                Eigen::Vector3d rkl = atom_k - atom_l;
                double tmp = 0;
                double dEdphi = (1 / 2.0 * dihedral.V * dihedral.n * (cos(dihedral.n * dihedral.phi0) * sin(dihedral.n * phi))) * m_final_factor * m_dihedral_scaling;
                if (isnan(dEdphi))
                    continue;
                /* some old code stored */
                //  Eigen::Vector3d dEdi = -1*dEdphi * rkj.norm()*m/(m.norm()*m.norm());
                //  Eigen::Vector3d dEdi = dEdphi * (rji/rji.norm()).cross(rkj/(rkj.norm()))/(sinphi*sinphi*rji.norm());
                //  Eigen::Vector3d dEdi = dEdphi * (atom_k - atom_j).norm()/(nabc.norm()*nabc.norm())*nabc;

                //  Eigen::Vector3d dEdl = dEdphi * rkj.norm()*n/(n.norm()*n.norm());
                //  Eigen::Vector3d dEdl = dEdphi * (-1*rkj/rkj.norm()).cross(rkl/(rkl.norm()))/(sinphi*sinphi*rkl.norm());

                //  Eigen::Vector3d dEdj = -1 * dEdphi * (-1 * (rji.dot(rkj) / (rkj.norm() * rkj.norm()) - 1) * (rkj.norm() / nabc.norm()) * nabc - rkj.dot(rkl) / (rkj.norm() * rkj.norm()) * nbcd / (nbcd.norm() * nbcd.norm()));
                //  Eigen::Vector3d dEdj = dEdi*-1*rji.norm()/(rkj.norm() * -1*cos(phi)) + dEdl*rkl.norm()/(rkj.norm() * -1*cos(phi));

                //  Eigen::Vector3d dEdk = -1 * dEdphi * ((rkj.dot(rkl) / (rkj.norm() * rkj.norm()) - 1) * (rkj.norm() / nbcd.norm()) * nbcd - rji.dot(rkj) / (rkj.norm() * rkj.norm()) * nabc / (nabc.norm() * nabc.norm()));
                //  Eigen::Vector3d dEdk = -1*dEdl - ((-1*rji).dot(rkj)/(rkj.norm()*rkj.norm())*dEdi) + (rkl.dot(rkj)/(rkj.norm()*rkj.norm())*dEdl);

                Eigen::Vector3d dEdi = dEdphi * rkj.norm() / (nabc.norm() * nabc.norm()) * nabc;
                Eigen::Vector3d dEdl = -1 * dEdphi * rkj.norm() / (nbcd.norm() * nbcd.norm()) * nbcd;
                Eigen::Vector3d dEdj = -1 * dEdi + ((-1 * rji).dot(rkj) / (rkj.norm() * rkj.norm()) * dEdi) - (rkl.dot(rkj) / (rkj.norm() * rkj.norm()) * dEdl);
                Eigen::Vector3d dEdk = -1 * (dEdi + dEdj + dEdl);

                if (isnan(dEdi.sum()) || isnan(dEdj.sum()) || isnan(dEdk.sum()) || isnan(dEdl.sum()))
                    continue;
                m_gradient(i, 0) += dEdi(0);
                m_gradient(i, 1) += dEdi(1);
                m_gradient(i, 2) += dEdi(2);

                m_gradient(l, 0) += dEdl(0);
                m_gradient(l, 1) += dEdl(1);
                m_gradient(l, 2) += dEdl(2);

                m_gradient(j, 0) += dEdj(0);
                m_gradient(j, 1) += dEdj(1);
                m_gradient(j, 2) += dEdj(2);

                m_gradient(k, 0) += dEdk(0);
                m_gradient(k, 1) += dEdk(1);
                m_gradient(k, 2) += dEdk(2);
            } else if (m_calc_gradient == 1) {
                e *= 1000;
                double tmp = 0.0;

                tmp = (Dihedral(AddVector(atom_i, dx), atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(SubVector(atom_i, dx), atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(i, 0) += 0;
                else
                    m_gradient(i, 0) += tmp;

                tmp = (Dihedral(AddVector(atom_i, dy), atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(SubVector(atom_i, dy), atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(i, 1) += 0;
                else
                    m_gradient(i, 1) += tmp;

                tmp = (Dihedral(AddVector(atom_i, dz), atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(SubVector(atom_i, dz), atom_j, atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(i, 2) += 0;
                else
                    m_gradient(i, 2) += tmp;

                tmp = (Dihedral(atom_i, AddVector(atom_j, dx), atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, SubVector(atom_j, dx), atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(j, 0) += 0;
                else
                    m_gradient(j, 0) += tmp;

                tmp = (Dihedral(atom_i, AddVector(atom_j, dy), atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, SubVector(atom_j, dy), atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(j, 1) += 0;
                else
                    m_gradient(j, 1) += tmp;

                tmp = (Dihedral(atom_i, AddVector(atom_j, dz), atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, SubVector(atom_j, dz), atom_k, atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(j, 2) += 0;
                else
                    m_gradient(j, 2) += tmp;

                tmp = (Dihedral(atom_i, atom_j, AddVector(atom_k, dx), atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, atom_j, SubVector(atom_k, dx), atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(k, 0) += 0;
                else
                    m_gradient(k, 0) += tmp;

                tmp = (Dihedral(atom_i, atom_j, AddVector(atom_k, dy), atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, atom_j, SubVector(atom_k, dy), atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(k, 1) += 0;
                else
                    m_gradient(k, 1) += tmp;

                tmp = (Dihedral(atom_i, atom_j, AddVector(atom_k, dz), atom_l, dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, atom_j, SubVector(atom_k, dz), atom_l, dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(k, 2) += 0;
                else
                    m_gradient(k, 2) += tmp;

                tmp = (Dihedral(atom_i, atom_j, atom_k, AddVector(atom_l, dx), dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, atom_j, atom_k, SubVector(atom_l, dx), dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(l, 0) += 0;
                else
                    m_gradient(l, 0) += tmp;

                tmp = (Dihedral(atom_i, atom_j, atom_k, AddVector(atom_l, dy), dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, atom_j, atom_k, SubVector(atom_l, dy), dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(l, 1) += 0;
                else
                    m_gradient(l, 1) += tmp;

                tmp = (Dihedral(atom_i, atom_j, atom_k, AddVector(atom_l, dz), dihedral.V, dihedral.n, dihedral.phi0) - Dihedral(atom_i, atom_j, atom_k, SubVector(atom_l, dz), dihedral.V, dihedral.n, dihedral.phi0)) / (2 * m_d);
                if (std::abs(tmp) > std::abs(e))
                    m_gradient(l, 2) += 0;
                else
                    m_gradient(l, 2) += tmp;
            }
        }
    }
    return energy;
}
double UFFThread::Inversion(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double k_ijkl, double C0, double C1, double C2)
{
    Eigen::Vector3d ail = SubVector(i, l);
    Eigen::Vector3d nijk = NormalVector(i, j, k);

    double cosY = (nijk.dot(ail) / ((nijk).norm() * (ail).norm()));

    double sinYSq = 1.0 - cosY * cosY;
    double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
    double cos2Y = sinY * sinY - 1.0;
    double energy = (k_ijkl * (C0 + C1 * sinY + C2 * cos2Y)) * m_final_factor * m_inversion_scaling;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFFThread::FullInversion(const int& i, const int& j, const int& k, const int& l, double d_forceConstant, double C0, double C1, double C2)
{
    double energy;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    Eigen::Vector3d atom_i = Position(i);
    Eigen::Vector3d atom_j = Position(j);
    Eigen::Vector3d atom_k = Position(k);
    Eigen::Vector3d atom_l = Position(l);
    energy += Inversion(atom_i, atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2);
    if (m_CalculateGradient) {
        if (m_calc_gradient == 0) {
            Eigen::Vector3d rji = (atom_j - atom_i);
            Eigen::Vector3d rjk = (atom_k - atom_i);
            Eigen::Vector3d rjl = (atom_l - atom_i);

            if (rji.norm() < 1e-5 || rjk.norm() < 1e-5 || rjl.norm() < 1e-5)
                return 0;

            double dji = rji.norm();
            double djk = rjk.norm();
            double djl = rjl.norm();
            rji /= dji;
            rjk /= djk;
            rjl /= djl;

            Eigen::Vector3d nijk = rji.cross(rjk);
            nijk /= nijk.norm();

            double cosY = (nijk.dot(rjl));
            double sinYSq = 1.0 - cosY * cosY;
            double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
            double cosTheta = (rji.dot(rjk));
            double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
            double sinTheta = std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

            double dEdY = -1 * (d_forceConstant * (C1 * cosY - 4 * C2 * cosY * sinY)) * m_final_factor * m_inversion_scaling;

            Eigen::Vector3d p1 = rji.cross(rjk);
            Eigen::Vector3d p2 = rjk.cross(rjl);
            Eigen::Vector3d p3 = rjl.cross(rji);

            const double sin_dl = p1.dot(rjl) / sinTheta;

            // the wilson angle:
            const double dll = asin(sin_dl);
            const double cos_dl = cos(dll);
            /*
                    double term1 = cosY * sinTheta;
                    double term2 = cosY / (sinY * sinThetaSq);

                    Eigen::Vector3d dYdi = (p1/term1 - (rji-rjk*cosTheta)*term2)/dji;
                    Eigen::Vector3d dYdk = (p2/term1 - (rjk-rji*cosTheta)*term2)/djk;
                    Eigen::Vector3d dYdl = (p3/term1 - (rjl*cosY/sinY))/djl;

             */
            Eigen::Vector3d dYdl = (p1 / sinTheta - (rjl * sin_dl)) / djl;
            Eigen::Vector3d dYdi = ((p2 + (((-rji + rjk * cosTheta) * sin_dl) / sinTheta)) / dji) / sinTheta;
            Eigen::Vector3d dYdk = ((p3 + (((-rjk + rji * cosTheta) * sin_dl) / sinTheta)) / djk) / sinTheta;
            Eigen::Vector3d dYdj = -1 * (dYdi + dYdk + dYdl);

            m_gradient(i, 0) += dEdY * dYdj(0);
            m_gradient(i, 1) += dEdY * dYdj(1);
            m_gradient(i, 2) += dEdY * dYdj(2);

            m_gradient(j, 0) += dEdY * dYdi(0);
            m_gradient(j, 1) += dEdY * dYdi(1);
            m_gradient(j, 2) += dEdY * dYdi(2);

            m_gradient(k, 0) += dEdY * dYdk(0);
            m_gradient(k, 1) += dEdY * dYdk(1);
            m_gradient(k, 2) += dEdY * dYdk(2);

            m_gradient(l, 0) += dEdY * dYdl(0);
            m_gradient(l, 1) += dEdY * dYdl(1);
            m_gradient(l, 2) += dEdY * dYdl(2);
        } else if (m_calc_gradient == 1) {
            m_gradient(i, 0) += (Inversion(AddVector(atom_i, dx), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dx), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(i, 1) += (Inversion(AddVector(atom_i, dy), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dy), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(i, 2) += (Inversion(AddVector(atom_i, dz), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(SubVector(atom_i, dz), atom_j, atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);

            m_gradient(j, 0) += (Inversion(atom_i, AddVector(atom_j, dx), atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dx), atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(j, 1) += (Inversion(atom_i, AddVector(atom_j, dy), atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dy), atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(j, 2) += (Inversion(atom_i, AddVector(atom_j, dz), atom_k, atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, SubVector(atom_j, dz), atom_k, atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);

            m_gradient(k, 0) += (Inversion(atom_i, atom_j, AddVector(atom_k, dx), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dx), atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(k, 1) += (Inversion(atom_i, atom_j, AddVector(atom_k, dy), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dy), atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(k, 2) += (Inversion(atom_i, atom_j, AddVector(atom_k, dz), atom_l, d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, SubVector(atom_k, dz), atom_l, d_forceConstant, C0, C1, C2)) / (2 * m_d);

            m_gradient(l, 0) += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dx), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dx), d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(l, 1) += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dy), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dy), d_forceConstant, C0, C1, C2)) / (2 * m_d);
            m_gradient(l, 2) += (Inversion(atom_i, atom_j, atom_k, AddVector(atom_l, dz), d_forceConstant, C0, C1, C2) - Inversion(atom_i, atom_j, atom_k, SubVector(atom_l, dz), d_forceConstant, C0, C1, C2)) / (2 * m_d);
        }
    }
    return energy;
}

double UFFThread::CalculateInversion()
{
    double energy = 0.0;
    for (int index = 0; index < m_uffinversion.size(); ++index) {
        const auto& inversion = m_uffinversion[index];
        const int i = inversion.i;
        const int j = inversion.j;
        const int k = inversion.k;
        const int l = inversion.l;
        energy += FullInversion(i, j, k, l, inversion.kijkl, inversion.C0, inversion.C1, inversion.C2);
    }
    return energy;
}

double UFFThread::NonBonds(const Eigen::Vector3d& i, const Eigen::Vector3d& j, double Dij, double xij)
{
    double r = (i - j).norm() * m_au;
    double pow6 = pow((xij / r), 6);
    double energy = Dij * (-2 * pow6 * m_vdw_scaling + pow6 * pow6 * m_rep_scaling) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double UFFThread::CalculateNonBonds()
{
    double energy = 0.0;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_uffvdwaals.size(); ++index) {
        const auto& vdw = m_uffvdwaals[index];
        const int i = vdw.i;
        const int j = vdw.j;
        Eigen::Vector3d atom_i = Position(i);
        Eigen::Vector3d atom_j = Position(j);
        double r = (atom_i - atom_j).norm() * m_au;
        double pow6 = pow((vdw.xij / r), 6);

        energy += vdw.Dij * (-2 * pow6 * m_vdw_scaling + pow6 * pow6 * m_rep_scaling) * m_final_factor;
        if (m_CalculateGradient) {
            if (m_calc_gradient == 0) {
                double diff = 12 * vdw.Dij * (pow6 * m_vdw_scaling - pow6 * pow6 * m_rep_scaling) / (r * r) * m_final_factor;
                m_gradient(i, 0) += diff * (atom_i(0) - atom_j(0));
                m_gradient(i, 1) += diff * (atom_i(1) - atom_j(1));
                m_gradient(i, 2) += diff * (atom_i(2) - atom_j(2));

                m_gradient(j, 0) -= diff * (atom_i(0) - atom_j(0));
                m_gradient(j, 1) -= diff * (atom_i(1) - atom_j(1));
                m_gradient(j, 2) -= diff * (atom_i(2) - atom_j(2));
            } else if (m_calc_gradient == 1) {
                m_gradient(i, 0) += (NonBonds(AddVector(atom_i, dx), atom_j, vdw.Dij, vdw.xij) - NonBonds(SubVector(atom_i, dx), atom_j, vdw.Dij, vdw.xij)) / (2 * m_d);
                m_gradient(i, 1) += (NonBonds(AddVector(atom_i, dy), atom_j, vdw.Dij, vdw.xij) - NonBonds(SubVector(atom_i, dy), atom_j, vdw.Dij, vdw.xij)) / (2 * m_d);
                m_gradient(i, 2) += (NonBonds(AddVector(atom_i, dz), atom_j, vdw.Dij, vdw.xij) - NonBonds(SubVector(atom_i, dz), atom_j, vdw.Dij, vdw.xij)) / (2 * m_d);

                m_gradient(j, 0) += (NonBonds(atom_i, AddVector(atom_j, dx), vdw.Dij, vdw.xij) - NonBonds(atom_i, SubVector(atom_j, dx), vdw.Dij, vdw.xij)) / (2 * m_d);
                m_gradient(j, 1) += (NonBonds(atom_i, AddVector(atom_j, dy), vdw.Dij, vdw.xij) - NonBonds(atom_i, SubVector(atom_j, dy), vdw.Dij, vdw.xij)) / (2 * m_d);
                m_gradient(j, 2) += (NonBonds(atom_i, AddVector(atom_j, dz), vdw.Dij, vdw.xij) - NonBonds(atom_i, SubVector(atom_j, dz), vdw.Dij, vdw.xij)) / (2 * m_d);
            }
        }
    }
    return energy;
}

double UFFThread::CalculateElectrostatic()
{
    double energy = 0.0;

    return energy;
}

eigenUFF::eigenUFF(const json& controller)
{
    json parameter = MergeJson(UFFParameterJson, controller);
    m_threadpool = new CxxThreadPool();
    m_threadpool->setProgressBar(CxxThreadPool::ProgressBarType::None);
#ifdef USE_D3
    m_use_d3 = parameter["d3"].get<int>();
    if (m_use_d3)
        m_d3 = new DFTD3Interface(controller);
#endif

#ifdef USE_D4
    m_use_d4 = parameter["d4"].get<int>();
    if (m_use_d4)
        m_d4 = new DFTD4Interface(controller);
#endif

    std::string param_file = parameter["param_file"];
    std::string uff_file = parameter["uff_file"];
    if (param_file.compare("none") != 0) {
        readParameterFile(param_file);
    }

    if (uff_file.compare("none") != 0) {
        readUFFFile(uff_file);
    }

    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = parameter["differential"].get<double>();

    readUFF(parameter);

    m_writeparam = parameter["writeparam"];
    m_writeuff = parameter["writeuff"];
    m_verbose = parameter["verbose"];
    m_rings = parameter["rings"];
    m_threads = parameter["threads"];
    m_scaling = 1.4;
    // m_au = au;
}

eigenUFF::~eigenUFF()
{
    delete m_threadpool;
    for (int i = 0; i < m_stored_threads.size(); ++i)
        delete m_stored_threads[i];
}

void eigenUFF::Initialise(const std::vector<std::pair<int, int>>& formed_bonds)
{
    if (m_initialised)
        return;

    m_uff_atom_types = std::vector<int>(m_atom_types.size(), 0);
    m_coordination = std::vector<int>(m_atom_types.size(), 0);
    std::vector<std::set<int>> ignored_vdw;
    m_topo = Eigen::MatrixXd::Zero(m_atom_types.size(), m_atom_types.size());
    TContainer bonds, nonbonds, angles, dihedrals, inversions;
    m_scaling = 1.4;
    m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    if (formed_bonds.size() == 0) {
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_stored_bonds.push_back(std::vector<int>());
            ignored_vdw.push_back(std::set<int>({ i }));
            for (int j = 0; j < m_atom_types.size() && m_stored_bonds[i].size() < CoordinationNumber[m_atom_types[i]]; ++j) {
                if (i == j)
                    continue;
                double x_i = m_geometry(i, 0) * m_au;
                double x_j = m_geometry(j, 0) * m_au;

                double y_i = m_geometry(i, 1) * m_au;
                double y_j = m_geometry(j, 1) * m_au;

                double z_i = m_geometry(i, 2) * m_au;
                double z_j = m_geometry(j, 2) * m_au;

                double r_ij = sqrt((((x_i - x_j) * (x_i - x_j)) + ((y_i - y_j) * (y_i - y_j)) + ((z_i - z_j) * (z_i - z_j))));

                if (r_ij <= (Elements::CovalentRadius[m_atom_types[i]] + Elements::CovalentRadius[m_atom_types[j]]) * m_scaling * m_au) {
                    if (bonds.insert({ std::min(i, j), std::max(i, j) })) {
                        m_coordination[i]++;
                        m_stored_bonds[i].push_back(j);
                        ignored_vdw[i].insert(j);
                    }
                    m_topo(i, j) = 1;
                    m_topo(j, i) = 1;
                }
            }
        }
    } else {
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_stored_bonds.push_back(std::vector<int>());
            ignored_vdw.push_back(std::set<int>({ i }));
        }
        for (const std::pair<int, int>& bond : formed_bonds) {

            int i = bond.first - 1;
            int j = bond.second - 1;

            if (bonds.insert({ std::min(i, j), std::max(i, j) })) {
                m_coordination[i]++;
                m_coordination[j]++;

                m_stored_bonds[i].push_back(j);
                m_stored_bonds[j].push_back(i);

                ignored_vdw[i].insert(j);
                ignored_vdw[j].insert(i);
            }
            m_topo(i, j) = 1;
            m_topo(j, i) = 1;
        }
    }
    AssignUffAtomTypes();
    if (m_rings)
        m_identified_rings = Topology::FindRings(m_stored_bonds, m_atom_types.size());

    bonds.clean();
    setBonds(bonds, ignored_vdw, angles, dihedrals, inversions);

    angles.clean();
    setAngles(angles, ignored_vdw);

    dihedrals.clean();
    setDihedrals(dihedrals);

    inversions.clean();
    setInversions(inversions);

    nonbonds.clean();
    setvdWs(ignored_vdw);

    m_h4correction.allocate(m_atom_types.size());

#ifdef USE_D3
    if (m_use_d3)
        m_d3->InitialiseMolecule(m_atom_types);
#endif

#ifdef USE_D4
    if (m_use_d4)
        m_d4->InitialiseMolecule(m_atom_types);
#endif

    if (m_writeparam.compare("none") != 0)
        writeParameterFile(m_writeparam + ".json");

    if (m_writeuff.compare("none") != 0)
        writeUFFFile(m_writeuff + ".json");

    AutoRanges();
    m_initialised = true;
}

void eigenUFF::setBonds(const TContainer& bonds, std::vector<std::set<int>>& ignored_vdw, TContainer& angels, TContainer& dihedrals, TContainer& inversions)
{
    for (const auto& bond : bonds.Storage()) {
        UFFBond b;

        b.i = bond[0];
        b.j = bond[1];
        int bond_order = 1;

        if (std::find(Conjugated.cbegin(), Conjugated.cend(), m_uff_atom_types[b.i]) != Conjugated.cend() && std::find(Conjugated.cbegin(), Conjugated.cend(), m_uff_atom_types[b.j]) != Conjugated.cend())
            bond_order = 2;
        else if (std::find(Triples.cbegin(), Triples.cend(), m_uff_atom_types[b.i]) != Triples.cend() || std::find(Triples.cbegin(), Triples.cend(), m_uff_atom_types[b.j]) != Triples.cend())
            bond_order = 3;
        else
            bond_order = 1;

        b.r0 = BondRestLength(b.i, b.j, bond_order);
        double cZi = UFFParameters[m_uff_atom_types[b.i]][cZ];
        double cZj = UFFParameters[m_uff_atom_types[b.j]][cZ];
        b.kij = 0.5 * m_bond_force * cZi * cZj / (b.r0 * b.r0 * b.r0);

        m_uffbonds.push_back(b);

        int i = bond[0];
        int j = bond[1];

        std::vector<int> k_bodies;
        for (auto t : m_stored_bonds[i]) {
            k_bodies.push_back(t);

            if (t == j)
                continue;
            angels.insert({ std::min(t, j), i, std::max(j, t) });
            ignored_vdw[i].insert(t);
        }

        std::vector<int> l_bodies;
        for (auto t : m_stored_bonds[j]) {
            l_bodies.push_back(t);

            if (t == i)
                continue;
            angels.insert({ std::min(i, t), j, std::max(t, i) });
            ignored_vdw[j].insert(t);
        }

        for (int k : k_bodies) {
            for (int l : l_bodies) {
                if (k == i || k == j || k == l || i == j || i == l || j == l)
                    continue;
                dihedrals.insert({ k, i, j, l });
                ignored_vdw[i].insert(k);
                ignored_vdw[i].insert(l);
                ignored_vdw[j].insert(k);
                ignored_vdw[j].insert(l);
                ignored_vdw[k].insert(l);
                ignored_vdw[l].insert(k);
            }
        }
        if (m_stored_bonds[i].size() == 3) {
            inversions.insert({ i, m_stored_bonds[i][0], m_stored_bonds[i][1], m_stored_bonds[i][2] });
        }
        if (m_stored_bonds[j].size() == 3) {
            inversions.insert({ j, m_stored_bonds[j][0], m_stored_bonds[j][1], m_stored_bonds[j][2] });
        }
    }
}

void eigenUFF::setAngles(const TContainer& angles, const std::vector<std::set<int>>& ignored_vdw)
{
    for (const auto& angle : angles.Storage()) {
        UFFAngle a;

        a.i = angle[0];
        a.j = angle[1];
        a.k = angle[2];
        if (a.i == a.j || a.i == a.k || a.j == a.k)
            continue;
        double f = pi / 180.0;
        double rij = BondRestLength(a.i, a.j, 1);
        double rjk = BondRestLength(a.j, a.k, 1);
        double Theta0 = UFFParameters[m_uff_atom_types[a.j]][cTheta0];
        double cosTheta0 = cos(Theta0 * f);
        double rik = sqrt(rij * rij + rjk * rjk - 2. * rij * rjk * cosTheta0);
        double param = m_angle_force;
        double beta = 2.0 * param / (rij * rjk);
        double preFactor = beta * UFFParameters[m_uff_atom_types[a.j]][cZ] * UFFParameters[m_uff_atom_types[a.k]][cZ] / (rik * rik * rik * rik * rik);
        double rTerm = rij * rjk;
        double inner = 3.0 * rTerm * (1.0 - cosTheta0 * cosTheta0) - rik * rik * cosTheta0;
        a.kijk = preFactor * rTerm * inner;
        a.C2 = 1 / (4 * std::max(sin(Theta0 * f) * sin(Theta0 * f), 1e-4));
        a.C1 = -4 * a.C2 * cosTheta0;
        a.C0 = a.C2 * (2 * cosTheta0 * cosTheta0 + 1);
        m_uffangle.push_back(a);
    }
}

void eigenUFF::setDihedrals(const TContainer& dihedrals)
{
    for (const auto& dihedral : dihedrals.Storage()) {
        UFFDihedral d;
        d.i = dihedral[0];
        d.j = dihedral[1];
        d.k = dihedral[2];
        d.l = dihedral[3];

        d.n = 2;
        double f = pi / 180.0;
        double bond_order = 1;
        d.V = 2;
        d.n = 3;
        d.phi0 = 180 * f;

        if (std::find(Conjugated.cbegin(), Conjugated.cend(), m_uff_atom_types[d.k]) != Conjugated.cend() && std::find(Conjugated.cbegin(), Conjugated.cend(), m_uff_atom_types[d.j]) != Conjugated.cend())
            bond_order = 2;
        else if (std::find(Triples.cbegin(), Triples.cend(), m_uff_atom_types[d.k]) != Triples.cend() || std::find(Triples.cbegin(), Triples.cend(), m_uff_atom_types[d.j]) != Triples.cend())
            bond_order = 3;
        else
            bond_order = 1;

        if (m_coordination[d.j] == 4 && m_coordination[d.k] == 4) // 2*sp3
        {
            d.V = sqrt(UFFParameters[m_uff_atom_types[d.j]][cV] * UFFParameters[m_uff_atom_types[d.k]][cV]);
            d.phi0 = 180 * f;
            d.n = 3;
        }
        if (m_coordination[d.j] == 3 && m_coordination[d.k] == 3) // 2*sp2
        {
            d.V = 5 * sqrt(UFFParameters[m_uff_atom_types[d.j]][cU] * UFFParameters[m_uff_atom_types[d.k]][cU]) * (1 + 4.18 * log(bond_order));
            d.phi0 = 180 * f;
            d.n = 2;
        } else if ((m_coordination[d.j] == 4 && m_coordination[d.k] == 3) || (m_coordination[d.j] == 3 && m_coordination[d.k] == 4)) {
            d.V = sqrt(UFFParameters[m_uff_atom_types[d.j]][cV] * UFFParameters[m_uff_atom_types[d.k]][cV]);
            d.phi0 = 0 * f;
            d.n = 6;

        } else {
            d.V = 5 * sqrt(UFFParameters[m_uff_atom_types[d.j]][cU] * UFFParameters[m_uff_atom_types[d.k]][cU]) * (1 + 4.18 * log(bond_order));
            d.phi0 = 90 * f;
        }

        m_uffdihedral.push_back(d);
    }
}

void eigenUFF::setInversions(const TContainer& inversions)
{
    for (const auto& inversion : inversions.Storage()) {
        const int i = inversion[0];
        if (m_coordination[i] != 3)
            continue;

        UFFInversion inv;
        inv.i = i;
        inv.j = inversion[1];
        inv.k = inversion[2];
        inv.l = inversion[3];

        double C0 = 0.0;
        double C1 = 0.0;
        double C2 = 0.0;
        double f = pi / 180.0;
        double kijkl = 0;
        if (6 <= m_atom_types[i] && m_atom_types[i] <= 8) {
            C0 = 1.0;
            C1 = -1.0;
            C2 = 0.0;
            kijkl = 6;
            if (m_atom_types[inv.j] == 8 || m_atom_types[inv.k] == 8 || m_atom_types[inv.l] == 8)
                kijkl = 50;
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
            kijkl = 22.0 / (C0 + C1 + C2);
        }
        inv.C0 = C0;
        inv.C1 = C1;
        inv.C2 = C2;
        inv.kijkl = kijkl;
        m_uffinversion.push_back(inv);
    }
}

void eigenUFF::setvdWs(const std::vector<std::set<int>>& ignored_vdw)
{
    for (int i = 0; i < m_atom_types.size(); ++i) {
        for (int j = i + 1; j < m_atom_types.size(); ++j) {
            if (std::find(ignored_vdw[i].begin(), ignored_vdw[i].end(), j) != ignored_vdw[i].end() || std::find(ignored_vdw[j].begin(), ignored_vdw[j].end(), i) != ignored_vdw[j].end())
                continue;
            UFFvdW v;
            v.i = i;
            v.j = j;

            double cDi = UFFParameters[m_uff_atom_types[v.i]][cD];
            double cDj = UFFParameters[m_uff_atom_types[v.j]][cD];
            double cxi = UFFParameters[m_uff_atom_types[v.i]][cx];
            double cxj = UFFParameters[m_uff_atom_types[v.j]][cx];
            v.Dij = sqrt(cDi * cDj) * 2;

            v.xij = sqrt(cxi * cxj);

            m_uffvdwaals.push_back(v);
        }
    }
}

void eigenUFF::AssignUffAtomTypes()
{
    for (int i = 0; i < m_atom_types.size(); ++i) {
        switch (m_atom_types[i]) {
        case 1: // Hydrogen
            if (m_stored_bonds[i].size() == 2)
                m_uff_atom_types[i] = 3; // Bridging Hydrogen
            else
                m_uff_atom_types[i] = 1;
            break;
        case 2: // Helium
            m_uff_atom_types[i] = 4;
            break;
        case 3: // Li
            m_uff_atom_types[i] = 5;
            break;
        case 4: // Be
            m_uff_atom_types[i] = 6;
            break;
        case 5: // B
            m_uff_atom_types[i] = 7;
            break;
        case 6: // C
            if (m_coordination[i] == 4)
                m_uff_atom_types[i] = 9;
            else if (m_coordination[i] == 3)
                m_uff_atom_types[i] = 10;
            else // if (coordination == 2)
                m_uff_atom_types[i] = 12;
            break;
        case 7: // N
            if (m_coordination[i] == 3)
                m_uff_atom_types[i] = 13;
            else if (m_coordination[i] == 2)
                m_uff_atom_types[i] = 14;
            else // if (coordination == 2)
                m_uff_atom_types[i] = 15;
            break;
        case 8: // O
            if (m_coordination[i] == 3)
                m_uff_atom_types[i] = 17;
            else if (m_coordination[i] == 2)
                m_uff_atom_types[i] = 19;
            else // if (coordination == 2)
                m_uff_atom_types[i] = 21;
            break;
        case 9: // F
            m_uff_atom_types[i] = 22;
            break;
        case 10: // Ne
            m_uff_atom_types[i] = 23;
            break;
        case 11: // Na
            m_uff_atom_types[i] = 24;
            break;
        case 12: // Mg
            m_uff_atom_types[i] = 25;
            break;
        case 13: // Al
            m_uff_atom_types[i] = 26;
            break;
        case 14: // Si
            m_uff_atom_types[i] = 27;
            break;
        case 15: // P
#pragma message("maybe add organometallic phosphorous (28)")
            m_uff_atom_types[i] = 29;
            break;
        case 16: // S
            if (m_coordination[i] == 2)
                m_uff_atom_types[i] = 31;
            else // ok, currently we do not discriminate between SO2 and SO3, just because there is H2SO3 and H2SO4
                m_uff_atom_types[i] = 32;
#pragma message("we have to add organic S")
            break;
        case 17: // Cl
            m_uff_atom_types[i] = 36;
            break;
        case 18: // Ar
            m_uff_atom_types[i] = 37;
            break;
        case 19: // K
            m_uff_atom_types[i] = 38;
            break;
        case 20: // Ca
            m_uff_atom_types[i] = 39;
            break;
        case 21: // Sc
            m_uff_atom_types[i] = 40;
            break;
        case 22: // Ti
            if (m_coordination[i] == 6)
                m_uff_atom_types[i] = 41;
            else
                m_uff_atom_types[i] = 42;
            break;
        case 23: // Va
            m_uff_atom_types[i] = 43;
            break;
        case 24: // Cr
            m_uff_atom_types[i] = 44;
            break;
        case 25: // Mn
            m_uff_atom_types[i] = 45;
            break;
        case 26: // Fe
            if (m_coordination[i] == 6)
                m_uff_atom_types[i] = 46;
            else
                m_uff_atom_types[i] = 47;
            break;
        case 27: // Co
            m_uff_atom_types[i] = 48;
            break;
        case 28: // Ni
            m_uff_atom_types[i] = 49;
            break;
        case 29: // Cu
            m_uff_atom_types[i] = 50;
            break;
        case 30: // Zn
            m_uff_atom_types[i] = 51;
            break;
        case 31: // Ga
            m_uff_atom_types[i] = 52;
            break;
        case 32: // Ge
            m_uff_atom_types[i] = 53;
            break;
        case 33: // As
            m_uff_atom_types[i] = 54;
            break;
        case 34: // Se
            m_uff_atom_types[i] = 55;
            break;
        case 35: // Br
            m_uff_atom_types[i] = 56;
            break;
        case 36: // Kr
            m_uff_atom_types[i] = 57;
            break;
        case 37: // Rb
            m_uff_atom_types[i] = 58;
            break;
        case 38: // Sr
            m_uff_atom_types[i] = 59;
            break;
        case 39: // Y
            m_uff_atom_types[i] = 60;
            break;
        case 40: // Zr
            m_uff_atom_types[i] = 61;
            break;
        case 41: // Nb
            m_uff_atom_types[i] = 62;
            break;
        case 42: // Mo
            if (m_coordination[i] == 6)
                m_uff_atom_types[i] = 63;
            else
                m_uff_atom_types[i] = 64;
            break;
        case 43: // Tc
            m_uff_atom_types[i] = 65;
            break;
        case 44: // Ru
            m_uff_atom_types[i] = 66;
            break;
        case 45: // Rh
            m_uff_atom_types[i] = 67;
            break;
        case 46: // Pd
            m_uff_atom_types[i] = 68;
            break;
        case 47: // Ag
            m_uff_atom_types[i] = 69;
            break;
        case 48: // Cd
            m_uff_atom_types[i] = 70;
            break;
        case 49: // In
            m_uff_atom_types[i] = 71;
            break;
        case 50: // Sn
            m_uff_atom_types[i] = 72;
            break;
        case 51: // Sb
            m_uff_atom_types[i] = 73;
            break;
        case 52: // Te
            m_uff_atom_types[i] = 74;
            break;
        case 53: // I
            m_uff_atom_types[i] = 75;
            break;
        case 54: // Xe
            m_uff_atom_types[i] = 76;
            break;
        default:
            m_uff_atom_types[i] = 0;
        };
        if (m_verbose) {
            std::cout << i << " " << m_atom_types[i] << " " << m_stored_bonds[i].size() << " " << m_uff_atom_types[i] << std::endl;
        }
    }
}

void eigenUFF::writeParameterFile(const std::string& file) const
{
    json parameters = writeUFF();
    parameters["bonds"] = Bonds();
    parameters["angles"] = Angles();
    parameters["dihedrals"] = Dihedrals();
    parameters["inversions"] = Inversions();
    parameters["vdws"] = vdWs();
    std::ofstream parameterfile(file);
    parameterfile << parameters;
}

void eigenUFF::writeUFFFile(const std::string& file) const
{
    json parameters = writeUFF();
    std::ofstream parameterfile(file);
    parameterfile << writeUFF();
}

json eigenUFF::Bonds() const
{
    json bonds;
    for (int i = 0; i < m_uffbonds.size(); ++i) {
        json bond;
        bond["i"] = m_uffbonds[i].i;
        bond["j"] = m_uffbonds[i].j;
        bond["r0"] = m_uffbonds[i].r0;
        bond["kij"] = m_uffbonds[i].kij;
        bonds[i] = bond;
    }
    return bonds;
}

json eigenUFF::Angles() const
{
    json angles;

    for (int i = 0; i < m_uffangle.size(); ++i) {
        json angle;
        angle["i"] = m_uffangle[i].i;
        angle["j"] = m_uffangle[i].j;
        angle["k"] = m_uffangle[i].k;

        angle["kijk"] = m_uffangle[i].kijk;
        angle["C0"] = m_uffangle[i].C0;
        angle["C1"] = m_uffangle[i].C1;
        angle["C2"] = m_uffangle[i].C2;

        angles[i] = angle;
    }
    return angles;
}

json eigenUFF::Dihedrals() const
{
    json dihedrals;
    for (int i = 0; i < m_uffdihedral.size(); ++i) {
        json dihedral;
        dihedral["i"] = m_uffdihedral[i].i;
        dihedral["j"] = m_uffdihedral[i].j;
        dihedral["k"] = m_uffdihedral[i].k;
        dihedral["l"] = m_uffdihedral[i].l;
        dihedral["V"] = m_uffdihedral[i].V;
        dihedral["n"] = m_uffdihedral[i].n;
        dihedral["phi0"] = m_uffdihedral[i].phi0;

        dihedrals[i] = dihedral;
    }
    return dihedrals;
}

json eigenUFF::Inversions() const
{
    json inversions;
    for (int i = 0; i < m_uffinversion.size(); ++i) {
        json inversion;
        inversion["i"] = m_uffinversion[i].i;
        inversion["j"] = m_uffinversion[i].j;
        inversion["k"] = m_uffinversion[i].k;
        inversion["l"] = m_uffinversion[i].l;
        inversion["kijkl"] = m_uffinversion[i].kijkl;
        inversion["C0"] = m_uffinversion[i].C0;
        inversion["C1"] = m_uffinversion[i].C1;
        inversion["C2"] = m_uffinversion[i].C2;

        inversions[i] = inversion;
    }
    return inversions;
}

json eigenUFF::vdWs() const
{
    json vdws;
    for (int i = 0; i < m_uffvdwaals.size(); ++i) {
        json vdw;
        vdw["i"] = m_uffvdwaals[i].i;
        vdw["j"] = m_uffvdwaals[i].j;
        vdw["Dij"] = m_uffvdwaals[i].Dij;
        vdw["xij"] = m_uffvdwaals[i].xij;
        vdws[i] = vdw;
    }
    return vdws;
}
/*
json eigenUFF::writeParameter() const
{
    json parameters;

    parameters["bonds"] = Bonds();
    parameters["angles"] = Angles();
    parameters["dihedrals"] = Dihedrals();
    parameters["inversions"] = Inversions();
    parameters["vdws"] = vdWs();

    parameters["bond_scaling"] = m_bond_scaling;
    parameters["angle_scaling"] = m_angle_scaling;
    parameters["inversion_scaling"] = m_inversion_scaling;
    parameters["vdw_scaling"] = m_vdw_scaling;
    parameters["rep_scaling"] = m_rep_scaling;
    parameters["dihedral_scaling"] = m_dihedral_scaling;
    parameters["gradient"] = m_calc_gradient;

    parameters["coulomb_scaling"] = m_coulmob_scaling;

    parameters["bond_force"] = m_bond_force;
    parameters["angle_force"] = m_angle_force;

    parameters["h4_scaling"] = m_h4_scaling;
    parameters["hh_scaling"] = m_hh_scaling;

    parameters["h4_oh_o"] = m_h4correction.get_OH_O();
    parameters["h4_oh_n"] = m_h4correction.get_OH_N();
    parameters["h4_nh_o"] = m_h4correction.get_NH_O();
    parameters["h4_nh_n"] = m_h4correction.get_NH_N();

    parameters["h4_wh_o"] = m_h4correction.get_WH_O();
    parameters["h4_nh4"] = m_h4correction.get_NH4();
    parameters["h4_coo"] = m_h4correction.get_COO();
    parameters["hh_rep_k"] = m_h4correction.get_HH_Rep_K();
    parameters["hh_rep_e"] = m_h4correction.get_HH_Rep_E();
    parameters["hh_rep_r0"] = m_h4correction.get_HH_Rep_R0();

#ifdef USE_D3
    if (m_use_d3) {
        parameters["d_s6"] = m_d3->ParameterS6();
        parameters["d_s8"] = m_d3->ParameterS8();
        parameters["d_s9"] = m_d3->ParameterS9();
        parameters["d_a1"] = m_d3->ParameterA1();
        parameters["d_a2"] = m_d3->ParameterA2();
    }
#endif

#ifdef USE_D4
    if (m_use_d4) {
        parameters["d_s6"] = m_d4->Parameter().s6;
        parameters["d_s8"] = m_d4->Parameter().s8;
        parameters["d_s10"] = m_d4->Parameter().s10;
        parameters["d_s9"] = m_d4->Parameter().s9;
        parameters["d_a1"] = m_d4->Parameter().a1;
        parameters["d_a2"] = m_d4->Parameter().a2;
    }
#endif
    return parameters;
}
*/
json eigenUFF::writeUFF() const
{
    json parameters;
    parameters["differential"] = m_d;
    parameters["bond_scaling"] = m_bond_scaling;
    parameters["angle_scaling"] = m_angle_scaling;
    parameters["inversion_scaling"] = m_inversion_scaling;
    parameters["vdw_scaling"] = m_vdw_scaling;
    parameters["rep_scaling"] = m_rep_scaling;
    parameters["dihedral_scaling"] = m_dihedral_scaling;
    parameters["gradient"] = m_calc_gradient;

    parameters["coulomb_scaling"] = m_coulmob_scaling;

    parameters["bond_force"] = m_bond_force;
    parameters["angle_force"] = m_angle_force;

    parameters["h4_scaling"] = m_h4_scaling;
    parameters["hh_scaling"] = m_hh_scaling;

    parameters["h4_oh_o"] = m_h4correction.get_OH_O();
    parameters["h4_oh_n"] = m_h4correction.get_OH_N();
    parameters["h4_nh_o"] = m_h4correction.get_NH_O();
    parameters["h4_nh_n"] = m_h4correction.get_NH_N();

    parameters["h4_wh_o"] = m_h4correction.get_WH_O();
    parameters["h4_nh4"] = m_h4correction.get_NH4();
    parameters["h4_coo"] = m_h4correction.get_COO();
    parameters["hh_rep_k"] = m_h4correction.get_HH_Rep_K();
    parameters["hh_rep_e"] = m_h4correction.get_HH_Rep_E();
    parameters["hh_rep_r0"] = m_h4correction.get_HH_Rep_R0();

#ifdef USE_D3
    if (m_use_d3) {
        parameters["d_s6"] = m_d3->ParameterS6();
        parameters["d_s8"] = m_d3->ParameterS8();
        parameters["d_s9"] = m_d3->ParameterS9();
        parameters["d_a1"] = m_d3->ParameterA1();
        parameters["d_a2"] = m_d3->ParameterA2();
    }
#endif

#ifdef USE_D4
    if (m_use_d4) {
        parameters["d4_s6"] = m_d4->Parameter().s6;
        parameters["d4_s8"] = m_d4->Parameter().s8;
        parameters["d4_s10"] = m_d4->Parameter().s10;
        parameters["d4_s9"] = m_d4->Parameter().s9;
        parameters["d4_a1"] = m_d4->Parameter().a1;
        parameters["d4_a2"] = m_d4->Parameter().a2;
    }
#endif
    return parameters;
}
void eigenUFF::readUFF(const json& parameters)
{
    json parameter = MergeJson(UFFParameterJson, parameters);

#ifdef USE_D3
    if (m_use_d3)
        m_d3->UpdateParameters(parameter);
#endif

#ifdef USE_D4
    if (m_use_d4)
        m_d4->UpdateParameters(parameter);
#endif

    m_d = parameter["differential"].get<double>();

    m_bond_scaling = parameter["bond_scaling"].get<double>();
    m_angle_scaling = parameter["angle_scaling"].get<double>();
    m_dihedral_scaling = parameter["dihedral_scaling"].get<double>();
    m_inversion_scaling = parameter["inversion_scaling"].get<double>();
    m_vdw_scaling = parameter["vdw_scaling"].get<double>();
    m_rep_scaling = parameter["rep_scaling"].get<double>();

    m_coulmob_scaling = parameter["coulomb_scaling"].get<double>();

    m_bond_force = parameter["bond_force"].get<double>();
    m_angle_force = parameter["angle_force"].get<double>();
    m_calc_gradient = parameter["gradient"].get<int>();
    m_h4_scaling = parameter["h4_scaling"].get<double>();
    m_hh_scaling = parameter["hh_scaling"].get<double>();

    m_h4correction.set_OH_O(parameter["h4_oh_o"].get<double>());
    m_h4correction.set_OH_N(parameter["h4_oh_n"].get<double>());
    m_h4correction.set_NH_O(parameter["h4_nh_o"].get<double>());
    m_h4correction.set_NH_N(parameter["h4_nh_n"].get<double>());

    m_h4correction.set_WH_O(parameter["h4_wh_o"].get<double>());
    m_h4correction.set_NH4(parameter["h4_nh4"].get<double>());
    m_h4correction.set_COO(parameter["h4_coo"].get<double>());
    m_h4correction.set_HH_Rep_K(parameter["hh_rep_k"].get<double>());
    m_h4correction.set_HH_Rep_E(parameter["hh_rep_e"].get<double>());
    m_h4correction.set_HH_Rep_R0(parameter["hh_rep_r0"].get<double>());
}

void eigenUFF::readParameter(const json& parameters)
{
    m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    readUFF(parameters);
    // while (m_gradient.size() < m_atom_types.size())
    //     m_gradient.push_back({ 0, 0, 0 });

    //  m_d = parameters["differential"].get<double>();
    /*
    #ifdef USE_D3
        if (m_use_d3)
            m_d3->UpdateParameters(parameters);
    #endif

    #ifdef USE_D4
        if (m_use_d4)
            m_d4->UpdateParameters(parameters);
    #endif

        m_bond_scaling = parameters["bond_scaling"].get<double>();
        m_angle_scaling = parameters["angle_scaling"].get<double>();
        m_dihedral_scaling = parameters["dihedral_scaling"].get<double>();
        m_inversion_scaling = parameters["inversion_scaling"].get<double>();
        m_vdw_scaling = parameters["vdw_scaling"].get<double>();
        m_rep_scaling = parameters["rep_scaling"].get<double>();

        m_coulmob_scaling = parameters["coulomb_scaling"].get<double>();

        m_bond_force = parameters["bond_force"].get<double>();
        m_angle_force = parameters["angle_force"].get<double>();
        m_calc_gradient = parameters["gradient"].get<int>();

        m_h4_scaling = parameters["h4_scaling"].get<double>();
        m_hh_scaling = parameters["hh_scaling"].get<double>();

        m_h4correction.set_OH_O(parameters["h4_oh_o"].get<double>());
        m_h4correction.set_OH_N(parameters["h4_oh_n"].get<double>());
        m_h4correction.set_NH_O(parameters["h4_nh_o"].get<double>());
        m_h4correction.set_NH_N(parameters["h4_nh_n"].get<double>());

        m_h4correction.set_WH_O(parameters["h4_wh_o"].get<double>());
        m_h4correction.set_NH4(parameters["h4_nh4"].get<double>());
        m_h4correction.set_COO(parameters["h4_coo"].get<double>());
        m_h4correction.set_HH_Rep_K(parameters["hh_rep_k"].get<double>());
        m_h4correction.set_HH_Rep_E(parameters["hh_rep_e"].get<double>());
        m_h4correction.set_HH_Rep_R0(parameters["hh_rep_r0"].get<double>());
        */
    setBonds(parameters["bonds"]);
    setAngles(parameters["angles"]);
    setDihedrals(parameters["dihedrals"]);
    setInversions(parameters["inversions"]);
    setvdWs(parameters["vdws"]);

    AutoRanges();
    m_initialised = true;
}

void eigenUFF::setBonds(const json& bonds)
{
    m_uffbonds.clear();
    for (int i = 0; i < bonds.size(); ++i) {
        json bond = bonds[i].get<json>();
        UFFBond b;

        b.i = bond["i"].get<int>();
        b.j = bond["j"].get<int>();
        b.r0 = bond["r0"].get<double>();
        b.kij = bond["kij"].get<double>();
        m_uffbonds.push_back(b);
    }
}

void eigenUFF::setAngles(const json& angles)
{
    m_uffangle.clear();
    for (int i = 0; i < angles.size(); ++i) {
        json angle = angles[i].get<json>();
        UFFAngle a;

        a.i = angle["i"].get<int>();
        a.j = angle["j"].get<int>();
        a.k = angle["k"].get<int>();
        a.C0 = angle["C0"].get<double>();
        a.C1 = angle["C1"].get<double>();
        a.C2 = angle["C2"].get<double>();
        a.kijk = angle["kijk"].get<double>();
        m_uffangle.push_back(a);
    }
}

void eigenUFF::setDihedrals(const json& dihedrals)
{
    m_uffdihedral.clear();
    for (int i = 0; i < dihedrals.size(); ++i) {
        json dihedral = dihedrals[i].get<json>();
        UFFDihedral d;

        d.i = dihedral["i"].get<int>();
        d.j = dihedral["j"].get<int>();
        d.k = dihedral["k"].get<int>();
        d.l = dihedral["l"].get<int>();
        d.V = dihedral["V"].get<double>();
        d.n = dihedral["n"].get<double>();
        d.phi0 = dihedral["phi0"].get<double>();
        m_uffdihedral.push_back(d);
    }
}

void eigenUFF::setInversions(const json& inversions)
{
    m_uffinversion.clear();
    for (int i = 0; i < inversions.size(); ++i) {
        json inversion = inversions[i].get<json>();
        UFFInversion inv;

        inv.i = inversion["i"].get<int>();
        inv.j = inversion["j"].get<int>();
        inv.k = inversion["k"].get<int>();
        inv.l = inversion["l"].get<int>();
        inv.kijkl = inversion["kijkl"].get<double>();
        inv.C0 = inversion["C0"].get<double>();
        inv.C1 = inversion["C1"].get<double>();
        inv.C2 = inversion["C2"].get<double>();

        m_uffinversion.push_back(inv);
    }
}

void eigenUFF::setvdWs(const json& vdws)
{
    m_uffvdwaals.clear();
    for (int i = 0; i < vdws.size(); ++i) {
        json vdw = vdws[i].get<json>();
        UFFvdW v;

        v.i = vdw["i"].get<int>();
        v.j = vdw["j"].get<int>();
        v.Dij = vdw["Dij"].get<double>();
        v.xij = vdw["xij"].get<double>();

        m_uffvdwaals.push_back(v);
    }
}

void eigenUFF::readUFFFile(const std::string& file)
{
    nlohmann::json parameters;
    std::ifstream parameterfile(file);
    try {
        parameterfile >> parameters;
    } catch (nlohmann::json::type_error& e) {
    } catch (nlohmann::json::parse_error& e) {
    }
    readUFF(parameters);
}

void eigenUFF::readParameterFile(const std::string& file)
{
    nlohmann::json parameters;
    std::ifstream parameterfile(file);
    try {
        parameterfile >> parameters;
    } catch (nlohmann::json::type_error& e) {
    } catch (nlohmann::json::parse_error& e) {
    }
    readParameter(parameters);
}

void eigenUFF::AutoRanges()
{
    for (int i = 0; i < m_threads; ++i) {
        UFFThread* thread = new UFFThread(i, m_threads);
        thread->readUFF(writeUFF());
        thread->setMolecule(m_atom_types, &m_geometry);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
        for (int j = int(i * m_uffbonds.size() / double(m_threads)); j < int((i + 1) * m_uffbonds.size() / double(m_threads)); ++j)
            thread->AddBond(m_uffbonds[j]);

        for (int j = int(i * m_uffangle.size() / double(m_threads)); j < int((i + 1) * m_uffangle.size() / double(m_threads)); ++j)
            thread->AddAngle(m_uffangle[j]);

        for (int j = int(i * m_uffdihedral.size() / double(m_threads)); j < int((i + 1) * m_uffdihedral.size() / double(m_threads)); ++j)
            thread->AddDihedral(m_uffdihedral[j]);

        for (int j = int(i * m_uffinversion.size() / double(m_threads)); j < int((i + 1) * m_uffinversion.size() / double(m_threads)); ++j)
            thread->AddInversion(m_uffinversion[j]);

        for (int j = int(i * m_uffvdwaals.size() / double(m_threads)); j < int((i + 1) * m_uffvdwaals.size() / double(m_threads)); ++j)
            thread->AddvdW(m_uffvdwaals[j]);
    }

    m_uff_bond_end = m_uffbonds.size();
    m_uff_angle_end = m_uffangle.size();
    m_uff_dihedral_end = m_uffdihedral.size();
    m_uff_inv_end = m_uffinversion.size();
    m_uff_vdw_end = m_uffvdwaals.size();
}

void eigenUFF::UpdateGeometry(const double* coord)
{
    if (m_gradient.rows() != m_atom_types.size())
        m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);

    if (m_gradient.size() != m_atom_types.size()) {
        m_h4correction.allocate(m_atom_types.size());
    }

    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry(i, 0) = coord[3 * i + 0] * au;
        m_geometry(i, 1) = coord[3 * i + 1] * au;
        m_geometry(i, 2) = coord[3 * i + 2] * au;

        m_gradient(i, 0) = 0;
        m_gradient(i, 1) = 0;
        m_gradient(i, 2) = 0;
    }
}

void eigenUFF::UpdateGeometry(const Matrix& geometry)
{
    if (m_gradient.rows() != m_atom_types.size())
        m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);

    if (m_gradient.size() != m_atom_types.size()) {
        m_h4correction.allocate(m_atom_types.size());
    }
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry(i, 0) = geometry(i, 0);
        m_geometry(i, 1) = geometry(i, 1);
        m_geometry(i, 2) = geometry(i, 2);

        m_gradient(i, 0) = 0;
        m_gradient(i, 1) = 0;
        m_gradient(i, 2) = 0;
    }
}

void eigenUFF::Gradient(double* grad) const
{
    double factor = 1;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        grad[3 * i + 0] = m_gradient(i, 0) * factor;
        grad[3 * i + 1] = m_gradient(i, 1) * factor;
        grad[3 * i + 2] = m_gradient(i, 2) * factor;
    }
}

void eigenUFF::NumGrad(double* grad)
{
    double dx = m_d;
    bool g = m_CalculateGradient;
    m_CalculateGradient = false;
    double E1, E2;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry(i, j) += dx;
            E1 = Calculate();
            m_geometry(i, j) -= 2 * dx;
            E2 = Calculate();
            grad[3 * i + j] = (E1 - E2) / (2 * dx);
            m_geometry(i, j) += dx;
        }
    }
    m_CalculateGradient = g;
}

Eigen::MatrixXd eigenUFF::NumGrad()
{
    Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);

    double dx = m_d;
    bool g = m_CalculateGradient;
    m_CalculateGradient = false;
    double E1, E2;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry(i, j) += dx;
            E1 = Calculate();
            m_geometry(i, j) -= 2 * dx;
            E2 = Calculate();
            gradient(i, j) = (E1 - E2) / (2 * dx);
            m_geometry(i, j) += dx;
        }
    }
    m_CalculateGradient = g;
    return gradient;
}

double eigenUFF::BondRestLength(int i, int j, double n)
{
    double cRi = UFFParameters[m_uff_atom_types[i]][cR];
    double cRj = UFFParameters[m_uff_atom_types[j]][cR];
    double cXii = UFFParameters[m_uff_atom_types[i]][cXi];
    double cXij = UFFParameters[m_uff_atom_types[j]][cXi];

    double lambda = 0.13332;
    double r_BO = -lambda * (cRi + cRj) * log(n);
    double r_EN = cRi * cRj * (sqrt(cXii) - sqrt(cXij)) * (sqrt(cXii) - sqrt(cXij)) / (cRi * cXii + cRj * cXij);
    double r_0 = cRi + cRj;
    return (r_0 + r_BO - r_EN) * m_au;
}

double eigenUFF::Calculate(bool grd, bool verbose)
{
    m_CalculateGradient = grd;
    hbonds4::atom_t geometry[m_atom_types.size()];
    for (int i = 0; i < m_atom_types.size(); ++i) {
        geometry[i].x = m_geometry(i,0) * m_au;
        geometry[i].y = m_geometry(i,1) * m_au;
        geometry[i].z = m_geometry(i,2) * m_au;
        geometry[i].e = m_atom_types[i];
        m_h4correction.GradientH4()[i].x = 0;
        m_h4correction.GradientH4()[i].y = 0;
        m_h4correction.GradientH4()[i].z = 0;

        m_h4correction.GradientHH()[i].x = 0;
        m_h4correction.GradientHH()[i].y = 0;
        m_h4correction.GradientHH()[i].z = 0;

#ifdef USE_D4
        if (m_use_d4)
            m_d4->UpdateAtom(i, m_geometry(i, 0) / au, m_geometry(i, 1) / au, m_geometry(i, 2) / au);
#endif

#ifdef USE_D3
        if (m_use_d3)
            m_d3->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
#endif
    }
    double energy = 0.0;
    double d4_energy = 0;
    double d3_energy = 0;
    double bond_energy = 0.0;
    double angle_energy = 0.0;
    double dihedral_energy = 0.0;
    double inversion_energy = 0.0;
    double vdw_energy = 0.0;

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_stored_threads[i]->UpdateGeometry(&m_geometry);
    }

    m_threadpool->Reset();
    m_threadpool->setActiveThreadCount(m_threads);

    m_threadpool->StartAndWait();
    m_threadpool->setWakeUp(m_threadpool->WakeUp() / 2.0);
    for (int i = 0; i < m_stored_threads.size(); ++i) {
        bond_energy += m_stored_threads[i]->BondEnergy();
        angle_energy += m_stored_threads[i]->AngleEnergy();
        dihedral_energy += m_stored_threads[i]->DihedralEnergy();
        inversion_energy += m_stored_threads[i]->InversionEnergy();
        vdw_energy += m_stored_threads[i]->VdWEnergy();
        m_gradient += m_stored_threads[i]->Gradient();
    }
    /* + CalculateElectrostatic(); */
    energy = bond_energy + angle_energy + dihedral_energy + inversion_energy + vdw_energy;
    /*
    if(m_calc_gradient == 2 && m_CalculateGradient)
    {
        m_gradient = NumGrad();
    }
    */
#ifdef USE_D3
    if (m_use_d3) {
        if (grd) {
            double grad[3 * m_atom_types.size()];
            d3_energy = m_d3->DFTD3Calculation(grad);
            for (int i = 0; i < m_atom_types.size(); ++i) {
                m_gradient(i, 0) += grad[3 * i + 0] * au;
                m_gradient(i, 1) += grad[3 * i + 1] * au;
                m_gradient(i, 2) += grad[3 * i + 2] * au;
            }
        } else
            d3_energy = m_d3->DFTD3Calculation(0);
    }
#endif

#ifdef USE_D4
    if (m_use_d4) {
        if (grd) {
            double grad[3 * m_atom_types.size()];
            d4_energy = m_d4->DFTD4Calculation(grad);
            for (int i = 0; i < m_atom_types.size(); ++i) {
                m_gradient(i, 0) += grad[3 * i + 0] * au;
                m_gradient(i, 1) += grad[3 * i + 1] * au;
                m_gradient(i, 2) += grad[3 * i + 2] * au;
            }
        } else
            d4_energy = m_d4->DFTD4Calculation(0);
    }
#endif

    double energy_h4 = 0;
    if (m_h4_scaling > 1e-8)
        energy_h4 = m_h4correction.energy_corr_h4(m_atom_types.size(), geometry);
    double energy_hh = 0;
    if (m_hh_scaling > 1e-8)
        m_h4correction.energy_corr_hh_rep(m_atom_types.size(), geometry);
    energy += m_final_factor * m_h4_scaling * energy_h4 + m_final_factor * m_hh_scaling * energy_hh + d3_energy + d4_energy;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_gradient(i, 0) += m_final_factor * m_h4_scaling * m_h4correction.GradientH4()[i].x + m_final_factor * m_hh_scaling * m_h4correction.GradientHH()[i].x;
        m_gradient(i, 1) += m_final_factor * m_h4_scaling * m_h4correction.GradientH4()[i].y + m_final_factor * m_hh_scaling * m_h4correction.GradientHH()[i].y;
        m_gradient(i, 2) += m_final_factor * m_h4_scaling * m_h4correction.GradientH4()[i].z + m_final_factor * m_hh_scaling * m_h4correction.GradientHH()[i].z;
    }

    if (verbose) {
        std::cout << "Total energy " << energy << " Eh. Sum of " << std::endl
                  << "Bond Energy " << bond_energy << " Eh" << std::endl
                  << "Angle Energy " << angle_energy << " Eh" << std::endl
                  << "Dihedral Energy " << dihedral_energy << " Eh" << std::endl
                  << "Inversion Energy " << inversion_energy << " Eh" << std::endl
                  << "Nonbonded Energy " << vdw_energy << " Eh" << std::endl
                  << "D3 Energy " << d3_energy << " Eh" << std::endl
                  << "D4 Energy " << d4_energy << " Eh" << std::endl
                  << "HBondCorrection " << m_final_factor * m_h4_scaling * energy_h4 << " Eh" << std::endl
                  << "HHRepCorrection " << m_final_factor * m_hh_scaling * energy_hh << " Eh" << std::endl
                  << std::endl;

        for (int i = 0; i < m_atom_types.size(); ++i) {
            std::cout << m_gradient(i, 0) << " " << m_gradient(i, 1) << " " << m_gradient(i, 2) << std::endl;
        }
    }
    return energy;
}
