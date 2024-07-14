/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "forcefieldfunctions.h"

#include "forcefield.h"

ForceFieldThread::ForceFieldThread(int thread, int threads)
    : m_thread(thread)
    , m_threads(threads)
{
    setAutoDelete(false);
    m_final_factor = 1; // / 2625.15 * 4.19;
    // m_d = parameters["differential"].get<double>();
    m_d = 1e-7;
}

int ForceFieldThread::execute()
{
    m_angle_energy = 0.0;
    m_bond_energy = 0.0;
    m_vdw_energy = 0;
    m_rep_energy = 0;
    m_inversion_energy = 0;
    m_dihedral_energy = 0;
    m_angle_energy = 0;
    m_bond_energy = 0.0;

    CalculateUFFBondContribution();
    CalculateUFFAngleContribution();
    CalculateUFFDihedralContribution();
    CalculateUFFInversionContribution();
    CalculateUFFvdWContribution();

    /*
    CalculateQMDFFBondContribution();
    CalculateQMDFFAngleContribution();
    CalculateQMDFFDihedralContribution();
    */
    return 0;
}

void ForceFieldThread::addBond(const Bond& bonds)
{
    if (bonds.type == 1)
        m_uff_bonds.push_back(bonds);
    else if (bonds.type == 2)
        m_qmdff_bonds.push_back(bonds);
}

void ForceFieldThread::addAngle(const Angle& angles)
{
    if (angles.type == 1)
        m_uff_angles.push_back(angles);
    else if (angles.type == 2)
        m_qmdff_angles.push_back(angles);
}

void ForceFieldThread::addDihedral(const Dihedral& dihedrals)
{
    if (dihedrals.type == 1)
        m_uff_dihedrals.push_back(dihedrals);
    else if (dihedrals.type == 2)
        m_qmdff_dihedrals.push_back(dihedrals);
}

void ForceFieldThread::addInversion(const Inversion& inversions)
{
    if (inversions.type == 1)
        m_uff_inversions.push_back(inversions);
    else if (inversions.type == 2)
        m_qmdff_inversions.push_back(inversions);
}

void ForceFieldThread::addvdW(const vdW& vdWs)
{
    if (vdWs.type == 1)
        m_uff_vdWs.push_back(vdWs);
}

void ForceFieldThread::addEQ(const EQ& EQs)
{
    if (EQs.type == 2)
        m_qmdff_EQs.push_back(EQs);
}

void ForceFieldThread::CalculateUFFBondContribution()
{
    double factor = m_final_factor * m_bond_scaling;

    for (int index = 0; index < m_uff_bonds.size(); ++index) {
        const auto& bond = m_uff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(i, j, derivate, m_calculate_gradient);

        m_bond_energy += (0.5 * bond.fc * (rij - bond.r0_ij) * (rij - bond.r0_ij)) * factor;
        if (m_calculate_gradient) {
            double diff = (bond.fc) * (rij - bond.r0_ij) * factor;
            m_gradient.row(bond.i) += diff * derivate.row(0);
            m_gradient.row(bond.j) += diff * derivate.row(1);
        }
    }
}

void ForceFieldThread::CalculateUFFAngleContribution()
{
    for (int index = 0; index < m_uff_angles.size(); ++index) {
        const auto& angle = m_uff_angles[index];
        auto i = m_geometry.row(angle.i);
        auto j = m_geometry.row(angle.j);
        auto k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_calculate_gradient);
        m_angle_energy += (angle.fc * (angle.C0 + angle.C1 * costheta + angle.C2 * (2 * costheta * costheta - 1))) * m_final_factor * m_angle_scaling;
        if (m_calculate_gradient) {
            double sintheta = sin(acos(costheta));
            double dEdtheta = -angle.fc * sintheta * (angle.C1 + 4 * angle.C2 * costheta) * m_final_factor * m_angle_scaling;
            m_gradient.row(angle.i) += dEdtheta * derivate.row(0);
            m_gradient.row(angle.j) += dEdtheta * derivate.row(1);
            m_gradient.row(angle.k) += dEdtheta * derivate.row(2);
        }
    }
}

void ForceFieldThread::CalculateUFFDihedralContribution()
{
    for (int index = 0; index < m_uff_dihedrals.size(); ++index) {
        const auto& dihedral = m_uff_dihedrals[index];
        Eigen::Vector3d i = m_geometry.row(dihedral.i);
        Eigen::Vector3d j = m_geometry.row(dihedral.j);
        Eigen::Vector3d k = m_geometry.row(dihedral.k);
        Eigen::Vector3d l = m_geometry.row(dihedral.l);
        Eigen::Vector3d nijk = UFF::NormalVector(i, j, k);
        Eigen::Vector3d njkl = UFF::NormalVector(j, k, l);
        double n_ijk = (nijk).norm();
        double n_jkl = (njkl).norm();
        double dotpr = nijk.dot(njkl);
        Eigen::Vector3d ji = j - i;

        double sign = (-1 * ji).dot(njkl) < 0 ? -1 : 1;
        double phi = pi + sign * acos(dotpr / (n_ijk * n_jkl));
        double sinphi = sin(phi);
        double tmp_energy = (1 / 2.0 * dihedral.V * (1 - cos(dihedral.n * dihedral.phi0) * cos(dihedral.n * phi))) * m_final_factor * m_dihedral_scaling;
        if (isnan(tmp_energy))
            continue;
        m_dihedral_energy += tmp_energy;
        if (m_calculate_gradient) {

            Eigen::Vector3d kj = k - j;
            Eigen::Vector3d kl = k - l;
            double tmp = 0;
            double dEdphi = (1 / 2.0 * dihedral.V * dihedral.n * (cos(dihedral.n * dihedral.phi0) * sin(dihedral.n * phi))) * m_final_factor * m_dihedral_scaling;
            if (isnan(dEdphi))
                continue;

            Eigen::Vector3d dEdi = dEdphi * kj.norm() / (nijk.norm() * nijk.norm()) * nijk;
            Eigen::Vector3d dEdl = -1 * dEdphi * kj.norm() / (njkl.norm() * njkl.norm()) * njkl;
            Eigen::Vector3d dEdj = -1 * dEdi + ((-1 * ji).dot(kj) / (kj.norm() * kj.norm()) * dEdi) - (kl.dot(kj) / (kj.norm() * kj.norm()) * dEdl);
            Eigen::Vector3d dEdk = -1 * (dEdi + dEdj + dEdl);

            if (isnan(dEdi.sum()) || isnan(dEdj.sum()) || isnan(dEdk.sum()) || isnan(dEdl.sum()))
                continue;
            m_gradient(dihedral.i, 0) += dEdi(0);
            m_gradient(dihedral.i, 1) += dEdi(1);
            m_gradient(dihedral.i, 2) += dEdi(2);

            m_gradient(dihedral.l, 0) += dEdl(0);
            m_gradient(dihedral.l, 1) += dEdl(1);
            m_gradient(dihedral.l, 2) += dEdl(2);

            m_gradient(dihedral.j, 0) += dEdj(0);
            m_gradient(dihedral.j, 1) += dEdj(1);
            m_gradient(dihedral.j, 2) += dEdj(2);

            m_gradient(dihedral.k, 0) += dEdk(0);
            m_gradient(dihedral.k, 1) += dEdk(1);
            m_gradient(dihedral.k, 2) += dEdk(2);
        }
    }
}

void ForceFieldThread::CalculateUFFInversionContribution()
{

    for (int index = 0; index < m_uff_inversions.size(); ++index) {
        const auto& inversion = m_uff_inversions[index];

        Eigen::Vector3d i = m_geometry.row(inversion.i);
        Eigen::Vector3d j = m_geometry.row(inversion.j);
        Eigen::Vector3d k = m_geometry.row(inversion.k);
        Eigen::Vector3d l = m_geometry.row(inversion.l);

        Eigen::Vector3d ail = SubVector(i, l);
        Eigen::Vector3d nijk = UFF::NormalVector(i, j, k);

        double cosY = (nijk.dot(ail) / ((nijk).norm() * (ail).norm()));

        double sinYSq = 1.0 - cosY * cosY;
        double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
        double cos2Y = sinY * sinY - 1.0;

        double tmp_energy = (inversion.fc * (inversion.C0 + inversion.C1 * sinY + inversion.C2 * cos2Y)) * m_final_factor * m_inversion_scaling;
        if (std::isnan(tmp_energy))
            continue;
        m_inversion_energy += tmp_energy;

        // energy += Inversion(i, j, k, l, d_forceConstant, C0, C1, C2);
        if (m_calculate_gradient) {
            Eigen::Vector3d ji = (j - i);
            Eigen::Vector3d jk = (k - i);
            Eigen::Vector3d jl = (l - i);

            if (ji.norm() < 1e-5 || jk.norm() < 1e-5 || jl.norm() < 1e-5)
                continue;

            double dji = ji.norm();
            double djk = jk.norm();
            double djl = jl.norm();
            ji /= dji;
            jk /= djk;
            jl /= djl;

            Eigen::Vector3d nijk = ji.cross(jk);
            nijk /= nijk.norm();

            double cosY = (nijk.dot(jl));
            double sinYSq = 1.0 - cosY * cosY;
            double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
            double cosTheta = (ji.dot(jk));
            double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
            double sinTheta = std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

            double dEdY = -1 * (inversion.fc * (inversion.C1 * cosY - 4 * inversion.C2 * cosY * sinY)) * m_final_factor * m_inversion_scaling;

            Eigen::Vector3d p1 = ji.cross(jk);
            Eigen::Vector3d p2 = jk.cross(jl);
            Eigen::Vector3d p3 = jl.cross(ji);

            const double sin_dl = p1.dot(jl) / sinTheta;

            // the wilson angle:
            const double dll = asin(sin_dl);
            const double cos_dl = cos(dll);

            Eigen::Vector3d dYdl = (p1 / sinTheta - (jl * sin_dl)) / djl;
            Eigen::Vector3d dYdi = ((p2 + (((-ji + jk * cosTheta) * sin_dl) / sinTheta)) / dji) / sinTheta;
            Eigen::Vector3d dYdk = ((p3 + (((-jk + ji * cosTheta) * sin_dl) / sinTheta)) / djk) / sinTheta;
            Eigen::Vector3d dYdj = -1 * (dYdi + dYdk + dYdl);

            m_gradient(inversion.i, 0) += dEdY * dYdj(0);
            m_gradient(inversion.i, 1) += dEdY * dYdj(1);
            m_gradient(inversion.i, 2) += dEdY * dYdj(2);

            m_gradient(inversion.j, 0) += dEdY * dYdi(0);
            m_gradient(inversion.j, 1) += dEdY * dYdi(1);
            m_gradient(inversion.j, 2) += dEdY * dYdi(2);

            m_gradient(inversion.k, 0) += dEdY * dYdk(0);
            m_gradient(inversion.k, 1) += dEdY * dYdk(1);
            m_gradient(inversion.k, 2) += dEdY * dYdk(2);

            m_gradient(inversion.l, 0) += dEdY * dYdl(0);
            m_gradient(inversion.l, 1) += dEdY * dYdl(1);
            m_gradient(inversion.l, 2) += dEdY * dYdl(2);
        }
    }
}

void ForceFieldThread::CalculateUFFvdWContribution()
{

    for (int index = 0; index < m_uff_vdWs.size(); ++index) {
        const auto& vdw = m_uff_vdWs[index];
        Eigen::Vector3d i = m_geometry.row(vdw.i);
        Eigen::Vector3d j = m_geometry.row(vdw.j);
        double ij = (i - j).norm() * m_au;
        double pow6 = pow((vdw.r0_ij / ij), 6);

        m_vdw_energy += vdw.C_ij * (-2 * pow6 * m_vdw_scaling) * m_final_factor;
        m_rep_energy += vdw.C_ij * (pow6 * pow6 * m_rep_scaling) * m_final_factor;
        if (m_calculate_gradient) {
            double diff = 12 * vdw.C_ij * (pow6 * m_vdw_scaling - pow6 * pow6 * m_rep_scaling) / (ij * ij) * m_final_factor;
            m_gradient(vdw.i, 0) += diff * (i(0) - j(0));
            m_gradient(vdw.i, 1) += diff * (i(1) - j(1));
            m_gradient(vdw.i, 2) += diff * (i(2) - j(2));

            m_gradient(vdw.j, 0) -= diff * (i(0) - j(0));
            m_gradient(vdw.j, 1) -= diff * (i(1) - j(1));
            m_gradient(vdw.j, 2) -= diff * (i(2) - j(2));
        }
    }
}

void ForceFieldThread::CalculateQMDFFBondContribution()
{
    double factor = m_final_factor * m_bond_scaling;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };

    for (int index = 0; index < m_qmdff_bonds.size(); ++index) {
        const auto& bond = m_qmdff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);

        // Matrix derivate;
        m_bond_energy += QMDFF::LJStretchEnergy((i - j).norm(), bond.r0_ij, bond.fc, bond.exponent) * factor;
        if (m_calculate_gradient) {
            /*if (m_calc_gradient == 0) {
                double diff = 0;
                m_gradient.row(a) += diff * derivate.row(0);
                m_gradient.row(b) += diff * derivate.row(1);

            } else if (m_calc_gradient == 1) {*/
            double d_x = (QMDFF::LJStretchEnergy((i + dx - j).norm(), bond.r0_ij, bond.fc, bond.exponent) - QMDFF::LJStretchEnergy((i - dx - j).norm(), bond.r0_ij, bond.fc, bond.exponent)) / (2 * m_d) * factor;
            double d_y = (QMDFF::LJStretchEnergy((i + dy - j).norm(), bond.r0_ij, bond.fc, bond.exponent) - QMDFF::LJStretchEnergy((i - dy - j).norm(), bond.r0_ij, bond.fc, bond.exponent)) / (2 * m_d) * factor;
            double d_z = (QMDFF::LJStretchEnergy((i + dz - j).norm(), bond.r0_ij, bond.fc, bond.exponent) - QMDFF::LJStretchEnergy((i - dz - j).norm(), bond.r0_ij, bond.fc, bond.exponent)) / (2 * m_d) * factor;
            m_gradient(bond.i, 0) += d_x;
            m_gradient(bond.i, 1) += d_y;
            m_gradient(bond.i, 2) += d_z;

            m_gradient(bond.j, 0) -= d_x;
            m_gradient(bond.j, 1) -= d_y;
            m_gradient(bond.j, 2) -= d_x;
            //}
        }
    }
}
/*
double ForceFieldThread::HarmonicBondStretching()
{
    double energy = 0;

    return energy;
}
*/
void ForceFieldThread::CalculateQMDFFAngleContribution()
{
    double threshold = 1e-2;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_qmdff_angles.size(); ++index) {
        const auto& angle = m_qmdff_angles[index];
        /*
        const int a = angle.a;
        const int b = angle.b;
        const int c = angle.c;
*/
        Vector i = m_geometry.row(angle.i);
        Vector j = m_geometry.row(angle.j);
        Vector k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = AngleBending(i, j, k, derivate, m_calculate_gradient);
        std::function<double(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double, double, double, double)> angle_function;
        if (std::abs(costheta - pi) < threshold)

            angle_function = [this](const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, double thetae, double fc, double r0_ij, double r0_ik) -> double {
                double val = QMDFF::LinearAngleBend(i, j, k, thetae, fc, r0_ij, r0_ik);
                if (std::isnan(val))
                    return 0;
                else
                    return val;
            };
        else
            angle_function = [this](const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, double thetae, double fc, double r0_ij, double r0_ik) -> double {
                double val = QMDFF::AngleBend(i, j, k, thetae, fc, r0_ij, r0_ik);
                if (std::isnan(val))
                    return 0;
                else
                    return val;
            };

        m_angle_energy += angle_function(i, j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik); //(angle.kijk * (angle.C0 + angle.C1 * costheta + angle.C2 * (2 * costheta * costheta - 1))) * m_final_factor * m_angle_scaling;

        if (m_calculate_gradient) {
            if (m_calc_gradient == 0) {

            } else if (m_calc_gradient == 1) {
                m_gradient(angle.i, 0) += (angle_function((i + dx), j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function((i - dx), j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
                m_gradient(angle.i, 1) += (angle_function((i + dy), j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function((i - dy), j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
                m_gradient(angle.i, 2) += (angle_function((i + dz), j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function((i - dz), j, k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);

                m_gradient(angle.j, 0) += (angle_function(i, (j + dx), k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function(i, (j - dx), k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
                m_gradient(angle.j, 1) += (angle_function(i, (j + dy), k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function(i, (j - dy), k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
                m_gradient(angle.j, 2) += (angle_function(i, (j + dz), k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function(i, (j - dz), k, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);

                m_gradient(angle.k, 0) += (angle_function(i, j, k + dx, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function(i, j, k - dx, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
                m_gradient(angle.k, 1) += (angle_function(i, j, k + dy, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function(i, j, k - dy, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
                m_gradient(angle.k, 2) += (angle_function(i, j, k + dz, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik) - angle_function(i, j, k - dz, angle.theta0_ijk, angle.fc, angle.r0_ij, angle.r0_ik)) / (2 * m_d);
            }
        }
    }
}

D3Thread::D3Thread(int thread, int threads)
    : ForceFieldThread(thread, threads)
{
    setAutoDelete(false);
    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = 1e-7;
#ifdef USE_D3
    m_d3 = new DFTD3Interface();
#endif
}

D3Thread::~D3Thread()
{
#ifdef USE_D3
    delete m_d3;
#endif
}

int D3Thread::execute()
{
#ifdef USE_D3
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_d3->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
    }

    if (m_calculate_gradient) {
        double grad[3 * m_atom_types.size()];
        m_vdw_energy = m_d3->DFTD3Calculation(grad);
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_gradient(i, 0) += grad[3 * i + 0] * au;
            m_gradient(i, 1) += grad[3 * i + 1] * au;
            m_gradient(i, 2) += grad[3 * i + 2] * au;
        }
    } else
        m_vdw_energy = m_d3->DFTD3Calculation(0);
#else
    std::cerr << "D3 is not included, sorry for that" << std::endl;
    exit(1);
#endif
    return 0;
}

H4Thread::H4Thread(int thread, int threads)
    : ForceFieldThread(thread, threads)
{
    setAutoDelete(false);
    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = 1e-7;
}

H4Thread::~H4Thread()
{
}

int H4Thread::execute()
{
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
    }

    m_vdw_energy = m_h4correction.energy_corr_h4(m_atom_types.size(), geometry) * m_vdw_scaling * m_final_factor;
    m_rep_energy = m_h4correction.energy_corr_hh_rep(m_atom_types.size(), geometry) * m_rep_scaling * m_final_factor;

    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_gradient(i, 0) += m_final_factor * m_vdw_scaling * m_h4correction.GradientH4()[i].x + m_final_factor * m_rep_scaling * m_h4correction.GradientHH()[i].x;
        m_gradient(i, 1) += m_final_factor * m_vdw_scaling * m_h4correction.GradientH4()[i].y + m_final_factor * m_rep_scaling * m_h4correction.GradientHH()[i].y;
        m_gradient(i, 2) += m_final_factor * m_vdw_scaling * m_h4correction.GradientH4()[i].z + m_final_factor * m_rep_scaling * m_h4correction.GradientHH()[i].z;
    }
    return 0;
}
