/*
 * <Simple QMDFF implementation for Cucuma. >
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
#include "src/core/dftd3interface.h"
#include "src/core/dftd4interface.h"

#include "src/core/qmdff_par.h"
#include "src/core/topology.h"

#include <Eigen/Dense>

#include "src/core/forcefieldderivaties.h"

#include "qmdff.h"

#include "json.hpp"
using json = nlohmann::json;

int QMDFFThread::execute()
{
    m_calc_gradient = 1;
    //    m_CalculateGradient = grd;
    m_d4_energy = 0;
    m_d3_energy = 0;
    m_bond_energy = CalculateStretchEnergy();
    m_angle_energy = CalculateAngleBending();
    m_dihedral_energy = CalculateDihedral();
    m_inversion_energy = CalculateInversion();
    m_vdw_energy = CalculateNonBonds();
    /* + CalculateElectrostatic(); */
    m_energy = m_bond_energy + m_angle_energy + m_dihedral_energy + m_inversion_energy + m_vdw_energy;
    return 0;
}

void QMDFFThread::readFF(const json& parameter)
{
    //  json parameter = MergeJson(UFFParameterJson, parameters);
    /*
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
        */
}

double QMDFFThread::Distance(double x1, double x2, double y1, double y2, double z1, double z2) const
{
    return sqrt((((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2))));
}

double QMDFFThread::StretchEnergy(double rAB, double reAB, double kAB, double a)
{
    double ratio = reAB / rAB;
    double energy = kAB * (1 + pow(ratio, a) - 2 * pow(ratio, a * 0.5));
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double QMDFFThread::CalculateStretchEnergy()
{
    double factor = 1;
    double energy = 0.0;

    for (int index = 0; index < m_qmdffbonds.size(); ++index) {
        const auto& bond = m_qmdffbonds[index];
        const int a = bond.a;
        const int b = bond.b;
        double xa = (*m_geometry)(a, 0) * m_au;
        double xb = (*m_geometry)(b, 0) * m_au;

        double ya = (*m_geometry)(a, 1) * m_au;
        double yb = (*m_geometry)(b, 1) * m_au;

        double za = (*m_geometry)(a, 2) * m_au;
        double zb = (*m_geometry)(b, 2) * m_au;

        Vector x = Position(a);
        Vector y = Position(b);
        Matrix derivate;
        energy += StretchEnergy(Distance(xa, xb, ya, yb, za, zb), bond.reAB, bond.kAB, bond.exponA);
        if (m_CalculateGradient) {
            if (m_calc_gradient == 0) {
                double diff = 0;
                m_gradient.row(a) += diff * derivate.row(0);
                m_gradient.row(b) += diff * derivate.row(1);

            } else if (m_calc_gradient == 1) {
                double xa = (*m_geometry)(a, 0) * m_au;
                double xb = (*m_geometry)(b, 0) * m_au;

                double ya = (*m_geometry)(a, 1) * m_au;
                double yb = (*m_geometry)(b, 1) * m_au;

                double za = (*m_geometry)(a, 2) * m_au;
                double zb = (*m_geometry)(b, 2) * m_au;

                double d_x = (StretchEnergy(Distance(xa + m_d, xb, ya, yb, za, zb), bond.reAB, bond.kAB, bond.exponA) - StretchEnergy(Distance(xa - m_d, xb, ya, yb, za, zb), bond.reAB, bond.kAB, bond.exponA)) / (2 * m_d);
                double d_y = (StretchEnergy(Distance(xa, xb, ya + m_d, yb, za, zb), bond.reAB, bond.kAB, bond.exponA) - StretchEnergy(Distance(xa, xb, ya - m_d, yb, za, zb), bond.reAB, bond.kAB, bond.exponA)) / (2 * m_d);
                double d_z = (StretchEnergy(Distance(xa, xb, ya, yb, za + m_d, zb), bond.reAB, bond.kAB, bond.exponA) - StretchEnergy(Distance(xa, xb, ya, yb, za - m_d, zb), bond.reAB, bond.kAB, bond.exponA)) / (2 * m_d);
                m_gradient(a, 0) += d_x;
                m_gradient(a, 1) += d_y;
                m_gradient(a, 2) += d_z;

                m_gradient(b, 0) -= d_x; //(StretchEnergy(Distance(xa, xb + m_d, ya, yb, za, zb), bond.reAB, bond.kAB, bond.exponA) - StretchEnergy(Distance(xa, xb - m_d, ya, yb, za, zb), bond.reAB, bond.kAB, bond.exponA)) / (2 * m_d);
                m_gradient(b, 1) -= d_y; //(StretchEnergy(Distance(xa, xb, ya, yb + m_d, za, zb), bond.reAB, bond.kAB, bond.exponA) - StretchEnergy(Distance(xa, xb, ya, yb - m_d, za, zb), bond.reAB, bond.kAB, bond.exponA)) / (2 * m_d);
                m_gradient(b, 2) -= d_z; //(StretchEnergy(Distance(xa, xb, ya, yb, za, zb + m_d), bond.reAB, bond.kAB, bond.exponA) - StretchEnergy(Distance(xa, xb, ya, yb, za, zb - m_d), bond.reAB, bond.kAB, bond.exponA)) / (2 * m_d);
            }
        }
    }

    return energy;
}

double QMDFFThread::AngleDamping(double rAB, double rAC, double reAB, double reAC)
{
    double kdamp = 1;
    double ratioAB = rAB / reAB;
    double ratioAC = rAC / reAC;

    double fABinv = 1 + kdamp * pow(ratioAB, 4);
    double fACinv = 1 + kdamp * pow(ratioAC, 4);
    double finv = fABinv * fACinv;
    return 1 / finv;
}

double QMDFFThread::AngleBend(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, double thetae, double kabc, double reAB, double reAC)
{
    Eigen::Vector3d vec_1 = a - b;
    Eigen::Vector3d vec_2 = a - c;
    double damp = AngleDamping(Distance(a(0), b(0), a(1), b(1), a(2), b(2)), Distance(a(0), c(0), a(1), c(1), a(2), c(2)), reAB, reAC);
    double costheta = (vec_1.dot(vec_2) / (sqrt(vec_1.dot(vec_1) * vec_2.dot(vec_2))));
    double costhetae = cos(thetae);
    double energy = (kabc * damp * (costhetae - costheta) * (costhetae - costheta)) * m_angle_scaling;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double QMDFFThread::LinearAngleBend(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, double thetae, double kabc, double reAB, double reAC)
{
    if (kabc < 0)
        return 0;
    Eigen::Vector3d vec_1 = a - b;
    Eigen::Vector3d vec_2 = a - c;
    double damp = AngleDamping(Distance(a(0), b(0), a(1), b(1), a(2), b(2)), Distance(a(0), c(0), a(1), c(1), a(2), c(2)), reAB, reAC);
    double theta = acos(vec_1.dot(vec_2) / (sqrt(vec_1.dot(vec_1) * vec_2.dot(vec_2))));
    double energy = (kabc * damp * (thetae - theta) * (thetae - theta)) * m_angle_scaling;
    if (isnan(energy))
        return 0;
    else
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

double QMDFFThread::TorsionDamping(double rCA, double rAB, double rBD, double reCA, double reAB, double reBD)
{
    double kdamp = 1;
    double ratioCA = rCA / reCA;
    double ratioAB = rAB / reAB;
    double ratioBD = rAB / reBD;

    double fCAinv = 1 + kdamp * pow(ratioCA, 4);
    double fABinv = 1 + kdamp * pow(ratioAB, 4);
    double fBDinv = 1 + kdamp * pow(ratioBD, 4);

    double finv = fCAinv * fABinv * fBDinv;
    return 1 / finv;
}

double QMDFFThread::Dihedral(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double V, double n, double phi0)
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

double QMDFFThread::CalculateDihedral()
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
double QMDFFThread::Inversion(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double k_ijkl, double C0, double C1, double C2)
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

double QMDFFThread::FullInversion(const int& i, const int& j, const int& k, const int& l, double d_forceConstant, double C0, double C1, double C2)
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

double QMDFFThread::CalculateInversion()
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

double QMDFFThread::NonBonds(const Eigen::Vector3d& i, const Eigen::Vector3d& j, double Dij, double xij)
{
    double r = (i - j).norm() * m_au;
    double pow6 = pow((xij / r), 6);
    double energy = Dij * (-2 * pow6 * m_vdw_scaling + pow6 * pow6 * m_rep_scaling) * m_final_factor;
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double QMDFFThread::CalculateNonBonds()
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

double QMDFFThread::CalculateElectrostatic()
{
    double energy = 0.0;

    return energy;
}

QMDFF::QMDFF(const json& controller)
{
    json parameter = MergeJson(QMDFFParameterJson, controller);
    // std::cout << parameter << std::endl;

    m_threadpool = new CxxThreadPool();
    m_threadpool->setProgressBar(CxxThreadPool::ProgressBarType::None);

    std::string param_file = parameter["param_file"];
    std::string uff_file = parameter["ff_file"];
    if (param_file.compare("none") != 0) {
        readParameterFile(param_file);
    }

    if (uff_file.compare("none") != 0) {
        readFFFile(uff_file);
    }

    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = parameter["differential"].get<double>();
    json d3settings = DFTD3Settings;
    // the following parameter are simply taken from the HF calculation with def2-SVP basis set (later should be not important)
    d3settings["d_s6"] = 1.0;
    d3settings["d_s8"] = s8;
    d3settings["d_s9"] = 1.0;

    d3settings["d_a1"] = a1;
    d3settings["d_a2"] = a2;
    d3settings["d_alp"] = 1;

    m_d3 = new DFTD3Interface(d3settings);
    m_h4correction.allocate(m_atom_types.size());

    readFF(parameter);

    m_writeparam = parameter["writeparam"];
    m_writeuff = parameter["writeff"];
    m_verbose = parameter["verbose"];
    m_rings = parameter["rings"];
    m_threads = parameter["threads"];
    m_const = parameter["const"];
    m_variable = parameter["variable"];
    // std::cout << "Calculating variable " << int(m_variable) << " and constants terms " << int(m_const) << std::endl;
    m_scaling = 1.4;
    // m_au = au;
}

QMDFF::~QMDFF()
{
    delete m_d3;
    delete m_threadpool;
    for (int i = 0; i < m_stored_threads.size(); ++i)
        delete m_stored_threads[i];
}

void QMDFF::setParameter(const json& parameter)
{
    m_threadpool->clear();
    m_stored_threads.clear();

    if (parameter.contains("bonds"))
        setBonds(parameter["bonds"]);
    if (parameter.contains("angles"))
        setAngles(parameter["angles"]);

    AutoRanges();
}

json QMDFF::writeParameter() const
{
    json parameters = writeFF();
    parameters["bonds"] = Bonds();
    // parameters["angles"] = Angles();
    // parameters["dihedrals"] = Dihedrals();
    // parameters["inversions"] = Inversions();
    // parameters["vdws"] = vdWs();
    return parameters;
}
void QMDFF::Initialise()
{
    if (m_initialised)
        return;
    // std::cout << "Initialising QMDFF (see S. Grimmme, J. Chem. Theory Comput. 2014, 10, 10, 4497–4514 [10.1021/ct500573f]) for the original publication!" << std::endl;

    m_uff_atom_types = std::vector<int>(m_atom_types.size(), 0);
    m_coordination = std::vector<int>(m_atom_types.size(), 0);
    std::vector<std::set<int>> ignored_vdw;
    m_topo = Eigen::MatrixXd::Zero(m_atom_types.size(), m_atom_types.size());
    TContainer bonds, nonbonds, angles, dihedrals, inversions;
    m_scaling = 1.4;
    m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
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

    if (m_rings)
        m_identified_rings = Topology::FindRings(m_stored_bonds, m_atom_types.size());
    /*
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
        */
    if (m_writeparam.compare("none") != 0)
        writeParameterFile(m_writeparam + ".json");

    if (m_writeuff.compare("none") != 0)
        writeFFFile(m_writeuff + ".json");

    AutoRanges();
    m_initialised = true;
}

void QMDFF::setMolecule(const std::vector<int>& atom_types, const std::vector<std::array<double, 3>>& geometry)
{
    m_atom_types = atom_types;
    m_geometry = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry(i, 0) = geometry[i][0];
        m_geometry(i, 1) = geometry[i][1];
        m_geometry(i, 2) = geometry[i][2];
    }
    m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    m_d3->InitialiseMolecule(m_atom_types);
    m_h4correction.allocate(m_atom_types.size());
}

void QMDFF::setMolecule(const std::vector<int>& atom_types, const Matrix& geometry)
{
    m_atom_types = atom_types;
    m_geometry = geometry;
    m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    m_d3->InitialiseMolecule(m_atom_types);
}

void QMDFF::writeParameterFile(const std::string& file) const
{
    json parameters = writeFF();
    parameters["bonds"] = Bonds();
    parameters["angles"] = Angles();
    parameters["dihedrals"] = Dihedrals();
    parameters["inversions"] = Inversions();
    parameters["vdws"] = vdWs();
    std::ofstream parameterfile(file);
    parameterfile << parameters;
}

void QMDFF::writeFFFile(const std::string& file) const
{
    json parameters = writeFF();
    std::ofstream parameterfile(file);
    parameterfile << writeFF();
}

json QMDFF::Bonds() const
{
    json bonds;
    for (int i = 0; i < m_qmdffbonds.size(); ++i) {
        json bond;
        bond["i"] = m_qmdffbonds[i].a;
        bond["j"] = m_qmdffbonds[i].b;
        bond["reAB"] = m_qmdffbonds[i].reAB;
        bond["kAB"] = m_qmdffbonds[i].kAB;
        bond["exponA"] = m_qmdffbonds[i].exponA;
        bond["distance"] = m_qmdffbonds[i].distance;
        bonds[i] = bond;
    }
    return bonds;
}

json QMDFF::Angles() const
{
    json angles;

    for (int i = 0; i < m_qmdffangle.size(); ++i) {
        json angle;
        angle["a"] = m_qmdffangle[i].a;
        angle["b"] = m_qmdffangle[i].b;
        angle["c"] = m_qmdffangle[i].c;

        angle["kabc"] = m_qmdffangle[i].kabc;
        angle["thetae"] = m_qmdffangle[i].thetae;
        angle["reAB"] = m_qmdffangle[i].reAB;
        angle["reAC"] = m_qmdffangle[i].reAC;

        angles[i] = angle;
    }
    return angles;
}

json QMDFF::Dihedrals() const
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

json QMDFF::Inversions() const
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

json QMDFF::vdWs() const
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

json QMDFF::writeFF() const
{
    json parameters;

    return parameters;
}
void QMDFF::readFF(const json& parameters)
{
    json parameter = MergeJson(UFFParameterJson, parameters);
    setParameter(parameter);
}

void QMDFF::readParameter(const json& parameters)
{
    m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);
    readFF(parameters);

    AutoRanges();
    m_initialised = true;
}

void QMDFF::setBonds(const json& bonds)
{
    m_qmdffbonds.clear();
    for (int i = 0; i < bonds.size(); ++i) {
        json bond = bonds[i].get<json>();
        QMDFFBond b;

        b.a = bond["a"].get<int>();
        b.b = bond["b"].get<int>();
        b.reAB = bond["reAB"].get<double>();
        b.kAB = bond["kAB"].get<double>();
        b.exponA = bond["exponA"].get<double>();
        b.distance = bond["distance"].get<int>();
        m_qmdffbonds.push_back(b);
    }
}

void QMDFF::setAngles(const json& angles)
{
    m_qmdffangle.clear();
    for (int i = 0; i < angles.size(); ++i) {
        json angle = angles[i].get<json>();
        QMDFFAngle a;

        a.a = angle["a"].get<int>();
        a.b = angle["b"].get<int>();
        a.c = angle["c"].get<int>();
        a.thetae = angle["thetae"].get<double>();
        a.reAB = angle["reAB"].get<double>();
        a.reAC = angle["reAC"].get<double>();
        a.kabc = angle["kabc"].get<double>();
        m_qmdffangle.push_back(a);
    }
}

void QMDFF::setDihedrals(const json& dihedrals)
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

void QMDFF::setInversions(const json& inversions)
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
/*
void QMDFF::setvdWs(const json& vdws)
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
*/
void QMDFF::readFFFile(const std::string& file)
{
    nlohmann::json parameters;
    std::ifstream parameterfile(file);
    try {
        parameterfile >> parameters;
    } catch (nlohmann::json::type_error& e) {
    } catch (nlohmann::json::parse_error& e) {
    }
    readFF(parameters);
}

void QMDFF::readParameterFile(const std::string& file)
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

void QMDFF::AutoRanges()
{
    for (int i = 0; i < m_threads; ++i) {
        QMDFFThread* thread = new QMDFFThread(i, m_threads);
        thread->readFF(writeFF());
        thread->setMolecule(m_atom_types, &m_geometry);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
        if (m_variable)
            for (int j = int(i * m_qmdffbonds.size() / double(m_threads)); j < int((i + 1) * m_qmdffbonds.size() / double(m_threads)); ++j)
                thread->AddBond(m_qmdffbonds[j]);
        if (m_variable)
            for (int j = int(i * m_qmdffangle.size() / double(m_threads)); j < int((i + 1) * m_qmdffangle.size() / double(m_threads)); ++j)
                thread->AddAngle(m_qmdffangle[j]);

        if (m_const)
            for (int j = int(i * m_uffdihedral.size() / double(m_threads)); j < int((i + 1) * m_uffdihedral.size() / double(m_threads)); ++j)
                thread->AddDihedral(m_uffdihedral[j]);

        if (m_const)
            for (int j = int(i * m_uffinversion.size() / double(m_threads)); j < int((i + 1) * m_uffinversion.size() / double(m_threads)); ++j)
                thread->AddInversion(m_uffinversion[j]);
        /*
                for (int j = int(i * m_uffvdwaals.size() / double(m_threads)); j < int((i + 1) * m_uffvdwaals.size() / double(m_threads)); ++j)
                    thread->AddvdW(m_uffvdwaals[j]);
                */
    }

    m_uff_bond_end = m_qmdffbonds.size();
    m_uff_angle_end = m_qmdffangle.size();
    m_uff_dihedral_end = m_uffdihedral.size();
    m_uff_inv_end = m_uffinversion.size();
    // m_uff_vdw_end = m_uffvdwaals.size();
}

void QMDFF::UpdateGeometry(const double* coord)
{
    if (m_gradient.rows() != m_atom_types.size())
        m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);

    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry(i, 0) = coord[3 * i + 0] * au;
        m_geometry(i, 1) = coord[3 * i + 1] * au;
        m_geometry(i, 2) = coord[3 * i + 2] * au;

        m_gradient(i, 0) = 0;
        m_gradient(i, 1) = 0;
        m_gradient(i, 2) = 0;
    }
}

void QMDFF::UpdateGeometry(const std::vector<std::array<double, 3>>& geometry)
{
    if (m_gradient.rows() != m_atom_types.size())
        m_gradient = Eigen::MatrixXd::Zero(m_atom_types.size(), 3);

    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_geometry(i, 0) = geometry[i][0];
        m_geometry(i, 1) = geometry[i][1];
        m_geometry(i, 2) = geometry[i][2];

        m_gradient(i, 0) = 0;
        m_gradient(i, 1) = 0;
        m_gradient(i, 2) = 0;
    }
}

void QMDFF::Gradient(double* grad) const
{
    double factor = 1;
    for (int i = 0; i < m_atom_types.size(); ++i) {
        grad[3 * i + 0] = m_gradient(i, 0) * factor;
        grad[3 * i + 1] = m_gradient(i, 1) * factor;
        grad[3 * i + 2] = m_gradient(i, 2) * factor;
    }
}

void QMDFF::NumGrad(double* grad)
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

Eigen::MatrixXd QMDFF::NumGrad()
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

double QMDFF::Calculate(bool grd, bool verbose)
{
    m_CalculateGradient = grd;
    hbonds4::atom_t geometry[m_atom_types.size()];
    for (int i = 0; i < m_atom_types.size(); ++i) {
        geometry[i].x = m_geometry(i, 0) * m_au;
        geometry[i].y = m_geometry(i, 1) * m_au;
        geometry[i].z = m_geometry(i, 2) * m_au;
        geometry[i].e = m_atom_types[i];
        m_h4correction.GradientH4()[i].x = 0;
        m_h4correction.GradientH4()[i].y = 0;
        m_h4correction.GradientH4()[i].z = 0;

        m_h4correction.GradientHH()[i].x = 0;
        m_h4correction.GradientHH()[i].y = 0;
        m_h4correction.GradientHH()[i].z = 0;
        m_d3->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
    }

    double energy = 0.0;

    double bond_energy = 0.0;
    double angle_energy = 0.0;
    double dihedral_energy = 0.0;
    double inversion_energy = 0.0;
    // double vdw_energy = 0.0;
    double d3_energy = 0;

    m_threadpool->setActiveThreadCount(m_threads);

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_stored_threads[i]->UpdateGeometry(&m_geometry);
    }

    m_threadpool->Reset();
    m_threadpool->StartAndWait();
    m_threadpool->setWakeUp(m_threadpool->WakeUp() / 2.0);
    for (int i = 0; i < m_stored_threads.size(); ++i) {
        bond_energy += m_stored_threads[i]->BondEnergy();
        angle_energy += m_stored_threads[i]->AngleEnergy();
        dihedral_energy += m_stored_threads[i]->DihedralEnergy();
        inversion_energy += m_stored_threads[i]->InversionEnergy();
        // vdw_energy += m_stored_threads[i]->VdWEnergy();
        m_gradient += m_stored_threads[i]->Gradient();
    }
    /* + CalculateElectrostatic(); */

    double energy_h4 = 0;
    double energy_hh = 0;

    if (m_const) {

        if (m_h4_scaling > 1e-8)
            energy_h4 = m_h4correction.energy_corr_h4(m_atom_types.size(), geometry);
        if (m_hh_scaling > 1e-8)
            m_h4correction.energy_corr_hh_rep(m_atom_types.size(), geometry);
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_gradient(i, 0) += m_final_factor * m_h4_scaling * m_h4correction.GradientH4()[i].x + m_final_factor * m_hh_scaling * m_h4correction.GradientHH()[i].x;
            m_gradient(i, 1) += m_final_factor * m_h4_scaling * m_h4correction.GradientH4()[i].y + m_final_factor * m_hh_scaling * m_h4correction.GradientHH()[i].y;
            m_gradient(i, 2) += m_final_factor * m_h4_scaling * m_h4correction.GradientH4()[i].z + m_final_factor * m_hh_scaling * m_h4correction.GradientHH()[i].z;
        }

        if (grd) {
            double grad[3 * m_atom_types.size()];
            d3_energy = m_d3->DFTD3Calculation(grad);
            for (int i = 0; i < m_atom_types.size(); ++i) {
                double val = grad[3 * i + 0] * au;
                if (!std::isnan(val) && std::abs(val) < 1e10)
                    m_gradient(i, 0) += val;
                val = grad[3 * i + 1] * au;
                if (!std::isnan(val) && std::abs(val) < 1e10)
                    m_gradient(i, 1) += val;
                val = grad[3 * i + 2] * au;
                if (!std::isnan(val) && std::abs(val) < 1e10)
                    m_gradient(i, 2) += val;
            }
        } else
            d3_energy = m_d3->DFTD3Calculation(0);
        energy += m_final_factor * m_h4_scaling * energy_h4 + m_final_factor * m_hh_scaling * energy_hh + d3_energy;
    }
    energy = bond_energy + angle_energy + dihedral_energy + inversion_energy; // + vdw_energy;
    if (verbose) {
        std::cout << "Total energy " << energy << " Eh. Sum of " << std::endl
                  << "Bond Energy " << bond_energy << " Eh" << std::endl
                  << "Angle Energy " << angle_energy << " Eh" << std::endl
                  << "Dihedral Energy " << dihedral_energy << " Eh" << std::endl
                  << "Inversion Energy " << inversion_energy << " Eh" << std::endl
                  << "Nonbonded Energy " << d3_energy + m_final_factor * m_h4_scaling * energy_h4 + m_final_factor * m_hh_scaling * energy_hh << " Eh" << std::endl
                  << std::endl;

        for (int i = 0; i < m_atom_types.size(); ++i) {
            std::cout << m_gradient(i, 0) << " " << m_gradient(i, 1) << " " << m_gradient(i, 2) << std::endl;
        }
    }
    return energy;
}
