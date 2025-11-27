/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "forcefieldderivaties.h"
#include "qmdff_par.h"
#include "uff_par.h"

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
    std::cout << "Forcefield threads" << std::endl;
    if (m_method == 1) {
        CalculateUFFBondContribution();
        CalculateUFFAngleContribution();

    } else if (m_method == 2) {
        CalculateQMDFFBondContribution();
        // CalculateUFFBondContribution();
        CalculateQMDFFAngleContribution();
    } else if (m_method == 3) {
        std::cout << "GFN-FF Energy Calculation Started in Thread " << m_thread << std::endl;
        // GFN-FF bonded terms
        CalculateGFNFFBondContribution();
        CalculateGFNFFAngleContribution();
        CalculateGFNFFDihedralContribution();
        CalculateGFNFFInversionContribution();

        // GFN-FF non-bonded pairwise parallelizable terms (Phase 4)
        CalculateGFNFFDispersionContribution();  // D3/D4 dispersion
        CalculateGFNFFRepulsionContribution();   // GFN-FF repulsion
        CalculateGFNFFCoulombContribution();     // EEQ Coulomb electrostatics

        // Legacy vdW (will be replaced by pairwise terms above)
        CalculateGFNFFvdWContribution();
    }

    if (m_method != 3) {
        CalculateUFFDihedralContribution();
        CalculateUFFInversionContribution();
        CalculateUFFvdWContribution();
    }
    CalculateESPContribution();
    /*
    CalculateQMDFFDihedralContribution();
    */
    return 0;
}

void ForceFieldThread::addBond(const Bond& bonds)
{
    // if (bonds.type == 1)
    m_uff_bonds.push_back(bonds);
    // else if (bonds.type == 2)
    //     m_qmdff_bonds.push_back(bonds);
}

void ForceFieldThread::addAngle(const Angle& angles)
{
    // if (angles.type == 1)
    m_uff_angles.push_back(angles);
    // else if (angles.type == 2)
    //     m_qmdff_angles.push_back(angles);
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
    // if (EQs.type == 2)
    m_EQs.push_back(EQs);
}

void ForceFieldThread::addGFNFFBond(const Bond& bonds)
{
    m_gfnff_bonds.push_back(bonds);
}

void ForceFieldThread::addGFNFFAngle(const Angle& angles)
{
    m_gfnff_angles.push_back(angles);
}

void ForceFieldThread::addGFNFFDihedral(const Dihedral& dihedrals)
{
    m_gfnff_dihedrals.push_back(dihedrals);
}

void ForceFieldThread::addGFNFFInversion(const Inversion& inversions)
{
    m_gfnff_inversions.push_back(inversions);
}

void ForceFieldThread::addGFNFFvdW(const vdW& vdWs)
{
    m_gfnff_vdWs.push_back(vdWs);
}

// Phase 4: GFN-FF pairwise non-bonded addition methods (Claude Generated 2025)

void ForceFieldThread::addGFNFFDispersion(const GFNFFDispersion& dispersion)
{
    m_gfnff_dispersions.push_back(dispersion);
}

void ForceFieldThread::addGFNFFRepulsion(const GFNFFRepulsion& repulsion)
{
    m_gfnff_repulsions.push_back(repulsion);
}

void ForceFieldThread::addGFNFFCoulomb(const GFNFFCoulomb& coulomb)
{
    m_gfnff_coulombs.push_back(coulomb);
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

        m_vdw_energy += vdw.C_ij * (-2 * pow6 * m_vdw_scaling) * m_final_factor / 100;
        m_rep_energy += vdw.C_ij * (pow6 * pow6 * m_rep_scaling) * m_final_factor / 100;
        if (m_calculate_gradient) {
            double diff = 12 * vdw.C_ij * (pow6 * m_vdw_scaling - pow6 * pow6 * m_rep_scaling) / (ij * ij) * m_final_factor / 100;
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
    m_d = 1e-5;
    double factor = m_final_factor * m_bond_scaling;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_uff_bonds.size(); ++index) {
        const auto& bond = m_uff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);
        Vector ij = i - j;
        double distance = (ij).norm();

        // Matrix derivate;
        double fc = bond.r0_ij;
        const double ratio = bond.r0_ij / distance;
        m_bond_energy += fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
        if (m_calculate_gradient) {

            double diff = 1 * fc * (-1 * bond.exponent * pow(ratio, bond.exponent - 1) + 2 * bond.exponent * 0.5 * pow(ratio, bond.exponent * 0.5 - 1));
            /*
                        Vector ijx = i+ dx - j;
                        double distancex1 = (ijx).norm();
                        double ratio = bond.r0_ij / distancex1;

                        double dx_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
                        ijx = i - dx - j;
                        distancex1 = (ijx).norm();
                        ratio = bond.r0_ij / distancex1;
                        double dx_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));


                        Vector ijy = i+ dy - j;
                        double distancey1 = (ijy).norm();
                        ratio = bond.r0_ij / distancey1;

                        double dy_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
                        ijy = i - dy - j;
                        distancey1 = (ijy).norm();
                        ratio = bond.r0_ij / distancey1;
                        double dy_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));


                        Vector ijz = i+ dz - j;
                        double distancez1 = (ijz).norm();
                        ratio = bond.r0_ij / distancez1;

                        double dz_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
                        ijz = i - dz - j;
                        distancez1 = (ijz).norm();
                        ratio = bond.r0_ij / distancez1;
                        double dz_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
            */
            // std::cout << (dx_p - dx_m)/(2.0*m_d) << " " << (dy_p - dy_m)/(2.0*m_d) <<" "<< (dz_p - dz_m)/(2.0*m_d) << " :: " << (diff * ij / (distance)).transpose() << std::endl;
            /*
            m_gradient(bond.i, 0) += (dx_p - dx_m)/(2*m_d);
            m_gradient(bond.i, 1) += (dy_p - dy_m)/(2*m_d);
            m_gradient(bond.i, 2) += (dz_p - dz_m)/(2*m_d);

            m_gradient(bond.j, 0) -= (dx_p - dx_m)/(2*m_d);
            m_gradient(bond.j, 1) -= (dy_p - dy_m)/(2*m_d);
            m_gradient(bond.j, 2) -= (dz_p - dz_m)/(2*m_d);
            */
            m_gradient.row(bond.i) += diff * ij / (distance);
            m_gradient.row(bond.j) -= diff * ij / (distance);
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
    for (int index = 0; index < m_uff_angles.size(); ++index) {
        const auto& angle = m_uff_angles[index];

        Vector i = m_geometry.row(angle.i);
        Vector j = m_geometry.row(angle.j);
        Vector k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta0_ijk = cos(angle.theta0_ijk * pi / 180.0);
        double costheta = 0;
        double energy = 0;
        double dEdTheta = 0;

        if (std::abs(costheta0_ijk + 1) < threshold) {
            costheta = UFF::AngleBending(i, j, k, derivate, true);
            energy = angle.fc * (costheta - costheta0_ijk) * (costheta - costheta0_ijk);
            dEdTheta = 2 * angle.fc * (costheta - costheta0_ijk);
        } else {
            costheta = UFF::AngleBending(i, j, k, derivate, true);
            energy = angle.fc * (costheta - costheta0_ijk) * (costheta - costheta0_ijk);
            dEdTheta = 2 * angle.fc * (costheta - costheta0_ijk);
        }

        m_angle_energy += energy;
        derivate *= dEdTheta;
        if (m_calculate_gradient) {
            m_gradient.row(angle.i) -= derivate.row(0);
            m_gradient.row(angle.j) -= derivate.row(1);
            m_gradient.row(angle.k) -= derivate.row(2);
        }
    }
}

void ForceFieldThread::CalculateESPContribution()
{

    for (int index = 0; index < m_EQs.size(); ++index) {
        const auto& eq = m_EQs[index];
        Vector i = m_geometry.row(eq.i);
        Vector j = m_geometry.row(eq.j);
        Vector ij = i - j;
        double distance = (ij).norm();
        m_eq_energy += eq.epsilon * eq.q_i * eq.q_j / (distance);
        if (m_calculate_gradient) {
            double diff = -1 * eq.epsilon * eq.q_i * eq.q_j / (distance * distance);
            m_gradient.row(eq.i) += diff * ij / distance;
            m_gradient.row(eq.j) -= diff * ij / distance;
        }
    }
}

#ifdef USE_D3
D3Thread::D3Thread(int thread, int threads)
    : ForceFieldThread(thread, threads)
{
    setAutoDelete(false);
    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = 1e-7;
    m_d3 = new DFTD3Interface();
}

D3Thread::~D3Thread()
{
    delete m_d3;
}

int D3Thread::execute()
{
    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_d3->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
    }

    if (m_calculate_gradient) {
        double grad[3 * m_atom_types.size()];
        m_vdw_energy = m_d3->Calculation(grad);
        for (int i = 0; i < m_atom_types.size(); ++i) {
            m_gradient(i, 0) += grad[3 * i + 0] * au;
            m_gradient(i, 1) += grad[3 * i + 1] * au;
            m_gradient(i, 2) += grad[3 * i + 2] * au;
        }
    } else
        m_vdw_energy = m_d3->Calculation(0);
    return 0;
}
#endif

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

void ForceFieldThread::CalculateGFNFFBondContribution()
{
    // Phase 1.3: Correct GFN-FF exponential bond potential
    // Formula from Fortran gfnff_engrad.F90:675-721
    // E_bond = k_b * exp(-α * (r - r₀)²)
    
    double factor = m_final_factor * m_bond_scaling;

    for (int index = 0; index < m_gfnff_bonds.size(); ++index) {
        const auto& bond = m_gfnff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(i, j, derivate, m_calculate_gradient);

        // GFN-FF exponential bond stretching: E = k_b * exp(-α * (r-r₀)²)
        // Note: k_b is stored as NEGATIVE in gfnff.cpp:1358 for attractive bonds
        double dr = rij - bond.r0_ij;           // r - r₀
        double alpha = bond.exponent;            // α stored in exponent field
        double k_b = bond.fc;                    // Force constant (already negative!)
        double exp_term = std::exp(-alpha * dr * dr);
        double energy = k_b * exp_term;          // k_b already contains sign
        std::cout << "Bond " << bond.i << "-" << bond.j << ": "
                  << "k_b=" << k_b << ", alpha=" << alpha << ", "
                  << "r=" << rij << ", r0=" << bond.r0_ij << ", "
                  << "dr=" << dr << ", exp=" << exp_term << ", "
                  << "E=" << energy << ", factor=" << factor << std::endl;
        m_bond_energy += energy * factor;

        if (m_calculate_gradient) {
            // dE/dr = -2*α*dr*E (chain rule)
            double dEdr = -2.0 * alpha * dr * energy;
            m_gradient.row(bond.i) += dEdr * factor * derivate.row(0);
            m_gradient.row(bond.j) += dEdr * factor * derivate.row(1);
        }
    }
}

void ForceFieldThread::CalculateGFNFFAngleContribution()
{
    // Phase 3: GFN-FF angle bending potential WITH distance-damping
    // Formula from Fortran gfnff_engrad.F90:857-916 (subroutine egbend)
    //
    // CRITICAL: GFN-FF applies distance damping to angle energies!
    // Damping function: damp = 1.0 / (1.0 + ((r² / rcut)²))
    // where rcut = atcuta * (rcov_i + rcov_j)²
    // Reference: gfnff_engrad.F90:1223-1232 (subroutine gfnffdampa)
    //
    // Full formula:
    //   E = k * (cosθ - cosθ₀)² * damp_ij * damp_jk
    //   where damp_ij and damp_jk are calculated from bond distances
    //
    // For linear angles (θ₀ ≈ π):
    //   E = k * (θ - θ₀)² * damp_ij * damp_jk (harmonic in θ)

    const double factor = m_final_factor * m_angle_scaling;
    const double pi = 3.14159265358979323846;
    const double linear_threshold = 1.0e-6;  // Fortran: pi-c0 .lt. 1.d-6
    const double atcuta = 0.595;  // GFN-FF angle damping parameter (gfnff_param.f90:463)

    for (int index = 0; index < m_gfnff_angles.size(); ++index) {
        const auto& angle = m_gfnff_angles[index];
        auto i = m_geometry.row(angle.i);
        auto j = m_geometry.row(angle.j);
        auto k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_calculate_gradient);

        // Clamp cosine to valid range [-1, 1] for numerical stability
        costheta = std::max(-1.0, std::min(1.0, costheta));

        double theta = std::acos(costheta);    // Current angle [0, π]
        double theta0 = angle.theta0_ijk;      // Equilibrium angle [0, π]
        double k_ijk = angle.fc;               // Force constant

        // ===== PHASE 5A: fqq charge-dependent angle correction =====
        // Reference: Fortran gfnff_ini.f90:1426-1430
        // Formula: fqq = 1.0 - (qa_j*qa_i + qa_j*qa_k) * qfacBEN
        // Parameter: qfacBEN = -0.54 (gfnff_param.f90:741)
        // Claude Generated (Nov 2025)
        const double qfacBEN = -0.54;
        double fqq = 1.0;  // Default (no correction)

        if (angle.i < m_eeq_charges.size() &&
            angle.j < m_eeq_charges.size() &&
            angle.k < m_eeq_charges.size()) {

            double qa_j = m_eeq_charges[angle.j];  // Center atom
            double qa_i = m_eeq_charges[angle.i];
            double qa_k = m_eeq_charges[angle.k];

            // Safety check: charges should be reasonable ([-1, 1] range)
            if (std::abs(qa_j) < 1.0 && std::abs(qa_i) < 1.0 && std::abs(qa_k) < 1.0) {
                double charge_product = qa_j * qa_i + qa_j * qa_k;
                fqq = 1.0 - charge_product * qfacBEN;

                // TODO Phase 5B (LOW PRIORITY): Metal case uses 2.5x stronger correction
                // Requires is_metal[] to be passed to threads
            }
        }

        // Apply fqq correction to force constant
        k_ijk *= fqq;

        double energy, dedtheta;

        // Check if equilibrium angle is linear (θ₀ ≈ π)
        if (std::abs(pi - theta0) < linear_threshold) {
            // Linear case: E = k*(θ - θ₀)² (harmonic in θ)
            // Fortran gfnff_engrad.F90:895-897
            double dtheta = theta - theta0;
            energy = k_ijk * dtheta * dtheta;
            dedtheta = 2.0 * k_ijk * dtheta;
        } else {
            // Normal case: E = k*(cosθ - cosθ₀)² (cosine-based)
            // Fortran gfnff_engrad.F90:899-900
            double costheta0 = std::cos(theta0);
            double dcostheta = costheta - costheta0;
            energy = k_ijk * dcostheta * dcostheta;

            // dE/dθ = dE/d(cosθ) * d(cosθ)/dθ
            //       = 2*k*(cosθ - cosθ₀) * (-sinθ)
            //       = 2*k*sinθ*(cosθ₀ - cosθ)
            double sintheta = std::sin(theta);
            dedtheta = 2.0 * k_ijk * sintheta * (costheta0 - costheta);
        }

        // ===== PHASE 3 ADDITION: Distance-dependent damping =====
        // Calculate damping factors for bonds i-j and j-k
        // Formula: damp = 1.0 / (1.0 + ((r² / rcut)²))
        // where rcut = atcuta * (rcov_i + rcov_j)²

        double r_ij_sq = (i - j).squaredNorm();  // r_ij²
        double r_jk_sq = (k - j).squaredNorm();  // r_jk²

        // Get covalent radii from atom types (using GFN-FF D3-style covalent radii in Bohr)
        // CRITICAL: GFN-FF uses D3-style covalent radii, NOT the angle radii!
        // Reference: gfnff_param.f90:381-404 (covalentRadD3 array)
        // Format: Angström values pre-converted to Bohr with factor aatoau * 4.0 / 3.0
        // where aatoau = 1.8897261246257702 (Angström to Bohr)
        static const std::vector<double> gfnff_d3_cov_radii_bohr = {
            0.32*1.88972612462*4.0/3.0, 0.46*1.88972612462*4.0/3.0,  // H,He
            1.20*1.88972612462*4.0/3.0, 0.94*1.88972612462*4.0/3.0, 0.77*1.88972612462*4.0/3.0, 0.75*1.88972612462*4.0/3.0, 0.71*1.88972612462*4.0/3.0, 0.63*1.88972612462*4.0/3.0, 0.64*1.88972612462*4.0/3.0, 0.67*1.88972612462*4.0/3.0,  // Li-Ne
            1.40*1.88972612462*4.0/3.0, 1.25*1.88972612462*4.0/3.0, 1.13*1.88972612462*4.0/3.0, 1.04*1.88972612462*4.0/3.0, 1.10*1.88972612462*4.0/3.0, 1.02*1.88972612462*4.0/3.0, 0.99*1.88972612462*4.0/3.0, 0.96*1.88972612462*4.0/3.0,  // Na-Ar
            1.76*1.88972612462*4.0/3.0, 1.54*1.88972612462*4.0/3.0,  // K,Ca
            1.33*1.88972612462*4.0/3.0, 1.22*1.88972612462*4.0/3.0, 1.21*1.88972612462*4.0/3.0, 1.10*1.88972612462*4.0/3.0, 1.07*1.88972612462*4.0/3.0,  // Sc-
            1.04*1.88972612462*4.0/3.0, 1.00*1.88972612462*4.0/3.0, 0.99*1.88972612462*4.0/3.0, 1.01*1.88972612462*4.0/3.0, 1.09*1.88972612462*4.0/3.0,  // -Zn
            1.12*1.88972612462*4.0/3.0, 1.09*1.88972612462*4.0/3.0, 1.15*1.88972612462*4.0/3.0, 1.10*1.88972612462*4.0/3.0, 1.14*1.88972612462*4.0/3.0, 1.17*1.88972612462*4.0/3.0,  // Ga-Kr
            1.89*1.88972612462*4.0/3.0, 1.67*1.88972612462*4.0/3.0,  // Rb,Sr
            1.47*1.88972612462*4.0/3.0, 1.39*1.88972612462*4.0/3.0, 1.32*1.88972612462*4.0/3.0, 1.24*1.88972612462*4.0/3.0, 1.15*1.88972612462*4.0/3.0,  // Y-
            1.13*1.88972612462*4.0/3.0, 1.13*1.88972612462*4.0/3.0, 1.08*1.88972612462*4.0/3.0, 1.15*1.88972612462*4.0/3.0, 1.23*1.88972612462*4.0/3.0,  // -Cd
            1.28*1.88972612462*4.0/3.0, 1.26*1.88972612462*4.0/3.0, 1.26*1.88972612462*4.0/3.0, 1.23*1.88972612462*4.0/3.0, 1.32*1.88972612462*4.0/3.0, 1.31*1.88972612462*4.0/3.0  // In-Xe
        };

        auto get_rcov_bohr = [&](int atomic_number) -> double {
            if (atomic_number >= 1 && atomic_number <= static_cast<int>(gfnff_d3_cov_radii_bohr.size())) {
                return gfnff_d3_cov_radii_bohr[atomic_number - 1];
            }
            return 1.0 * 1.88972612462 * 4.0 / 3.0;  // Fallback for unknown elements (1.0 Å converted to Bohr)
        };

        double rcov_i = get_rcov_bohr(m_atom_types[angle.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[angle.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[angle.k]);

        // Calculate rcut values
        double rcut_ij_sq = atcuta * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk_sq = atcuta * (rcov_j + rcov_k) * (rcov_j + rcov_k);

        // Calculate damping factors: damp = 1.0 / (1.0 + (r²/rcut²)²)
        double rr_ij = (r_ij_sq / rcut_ij_sq);
        double rr_jk = (r_jk_sq / rcut_jk_sq);
        rr_ij = rr_ij * rr_ij;  // (r²/rcut²)²
        rr_jk = rr_jk * rr_jk;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);
        double damp = damp_ij * damp_jk;

        // ===== PHASE 6: Damping derivatives for gradient =====
        // Claude Generated (Nov 2025): Complete angle gradients with damping derivatives
        // Reference: Fortran gfnff_engrad.F90:890-892, 911-915, 1223-1232
        //
        // Damping derivative: ∂damp/∂(r²) = -4*rr / (r² * (1 + rr)²)
        // where rr = (r²/rcut²)²
        //
        // Formula (Fortran gfnff_engrad.F90:1231):
        //   ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)
        double damp2ij = -2.0 * 2.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij));
        double damp2jk = -2.0 * 2.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk));

        // Phase 3: Apply distance-dependent damping to energy
        // Fortran formula: e = ea*damp where damp = damp_ij*damp_jk
        m_angle_energy += energy * damp * factor;

        if (m_calculate_gradient) {
            // ===== COMPLETE GRADIENT WITH DAMPING DERIVATIVES =====
            // Claude Generated (Nov 2025): Full GFN-FF angle gradients
            // Reference: Fortran gfnff_engrad.F90:904-915
            //
            // Complete gradient formula:
            //   ∂E/∂x = (∂E/∂θ * damp) * (∂θ/∂x) + (∂E/∂damp) * (∂damp/∂x)
            //         = dedθ * damp * derivate + ea * (∂damp/∂x)
            //
            // where ∂damp/∂x has contributions from both damp_ij and damp_jk

            // Distance vectors (NOT normalized - Fortran lines 880-881)
            // vab = xyz(:,i) - xyz(:,j)  →  vector from j to i
            // vcb = xyz(:,k) - xyz(:,j)  →  vector from j to k
            Vector vab = i - j;
            Vector vcb = k - j;

            // Damping gradient contributions (Fortran lines 911-912)
            // term1 = ea * damp2ij * dampjk * vab
            // term2 = ea * damp2jk * dampij * vcb
            Vector term1 = energy * damp2ij * damp_jk * vab;
            Vector term2 = energy * damp2jk * damp_ij * vcb;

            // Angle gradient contributions (already computed by UFF::AngleBending)
            // Note: derivate matrix contains ∂θ/∂x for all 3 atoms
            Vector grad_angle_i = dedtheta * damp * factor * derivate.row(0);
            Vector grad_angle_j = dedtheta * damp * factor * derivate.row(1);
            Vector grad_angle_k = dedtheta * damp * factor * derivate.row(2);

            // Complete gradients with damping terms (Fortran lines 913-915)
            // g(:,1) = -dedb*damp - term1 - term2  (Atom i: negative signs!)
            // g(:,2) =  deda*damp + term1          (Atom j: center)
            // g(:,3) =  dedc*damp + term2          (Atom k)
            m_gradient.row(angle.i) += grad_angle_i - term1 - term2;
            m_gradient.row(angle.j) += grad_angle_j + term1;
            m_gradient.row(angle.k) += grad_angle_k + term2;
        }
    }
}

void ForceFieldThread::CalculateGFNFFDihedralContribution()
{
    // GFN-FF torsion: E = V*(1 + cos(n*φ - φ0))

    for (int index = 0; index < m_gfnff_dihedrals.size(); ++index) {
        const auto& dihedral = m_gfnff_dihedrals[index];

        auto i = m_geometry.row(dihedral.i);
        auto j = m_geometry.row(dihedral.j);
        auto k = m_geometry.row(dihedral.k);
        auto l = m_geometry.row(dihedral.l);

        Matrix derivate;
        double phi = UFF::Torsion(i, j, k, l, derivate, m_calculate_gradient);

        double V = dihedral.V;
        double n = dihedral.n;
        double phi0 = dihedral.phi0;

        double energy = V * (1 + cos(n * phi - phi0));
        m_dihedral_energy += energy * m_final_factor * m_dihedral_scaling;

        if (m_calculate_gradient) {
            double dEdphi = -V * n * sin(n * phi - phi0) * m_final_factor * m_dihedral_scaling;
            m_gradient.row(dihedral.i) += dEdphi * derivate.row(0);
            m_gradient.row(dihedral.j) += dEdphi * derivate.row(1);
            m_gradient.row(dihedral.k) += dEdphi * derivate.row(2);
            m_gradient.row(dihedral.l) += dEdphi * derivate.row(3);
        }
    }
}

void ForceFieldThread::CalculateGFNFFInversionContribution()
{
    // GFN-FF inversion: E = k*(C0 + C1*cos(θ) + C2*cos(2θ))

    for (int index = 0; index < m_gfnff_inversions.size(); ++index) {
        const auto& inversion = m_gfnff_inversions[index];

        auto i = m_geometry.row(inversion.i); // out-of-plane atom
        auto j = m_geometry.row(inversion.j); // plane atom 1
        auto k = m_geometry.row(inversion.k); // plane atom 2
        auto l = m_geometry.row(inversion.l); // central atom

        Matrix derivate;
        double theta = UFF::OutOfPlane(i, j, k, l, derivate, m_calculate_gradient);

        double energy = inversion.fc * (inversion.C0 + inversion.C1 * cos(theta) + inversion.C2 * cos(2 * theta));
        m_inversion_energy += energy * m_final_factor * m_inversion_scaling;

        if (m_calculate_gradient) {
            double dEdtheta = inversion.fc * (-inversion.C1 * sin(theta) - 2 * inversion.C2 * sin(2 * theta)) * m_final_factor * m_inversion_scaling;
            m_gradient.row(inversion.i) += dEdtheta * derivate.row(0);
            m_gradient.row(inversion.j) += dEdtheta * derivate.row(1);
            m_gradient.row(inversion.k) += dEdtheta * derivate.row(2);
            m_gradient.row(inversion.l) += dEdtheta * derivate.row(3);
        }
    }
}

// DEPRECATED: Legacy vdW calculation replaced by Phase 4 pairwise terms
// (CalculateGFNFFDispersionContribution + CalculateGFNFFRepulsionContribution)
// TODO: Remove after verifying all tests pass without this
void ForceFieldThread::CalculateGFNFFvdWContribution()
{
#pragma message("TODO (DEPRECATED - use Phase 4 pairwise): Implement proper GFN-FF dispersion with D4-like correction")
    // TODO: GFN-FF uses sophisticated dispersion correction (similar to D4)
    // This is a simplified implementation for testing

    for (int index = 0; index < m_gfnff_vdWs.size(); ++index) {
        const auto& vdw = m_gfnff_vdWs[index];
        auto i = m_geometry.row(vdw.i);
        auto j = m_geometry.row(vdw.j);

        Vector rij = i - j;
        double distance = rij.norm();

        if (distance > 0) {
            // Simplified van der Waals interaction
            // TODO (DEPRECATED): Replace with proper GFN-FF dispersion calculation (use Phase 4 pairwise instead)
            double C6 = vdw.C_ij;
            double r0 = vdw.r0_ij;
            double ratio = r0 / distance;

            // Simple Lennard-Jones-like term (placeholder)
            double energy = -C6 * pow(ratio, 6);

            m_vdw_energy += energy * m_final_factor * m_vdw_scaling;

            if (m_calculate_gradient) {
                double dEdr = 6.0 * C6 * pow(ratio, 6) / distance * m_final_factor * m_vdw_scaling;
                Vector grad = dEdr * rij / distance;
                m_gradient.row(vdw.i) += grad.transpose();
                m_gradient.row(vdw.j) -= grad.transpose();
            }
        }
    }
}

// ============================================================================
// Phase 4: GFN-FF Pairwise Non-Bonded Terms (Claude Generated 2025)
// ============================================================================
// These functions implement GFN-FF non-bonded interactions as pairwise
// parallelizable terms following the UFF vdW pattern (NOT as add-on corrections).
// Each pair (i,j) is pre-computed with full parameters and stored in vectors.

void ForceFieldThread::CalculateGFNFFDispersionContribution()
{
    /**
     * @brief D3/D4 Dispersion with Becke-Johnson damping (pairwise parallelizable)
     *
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010) [D3-BJ damping]
     * Formula: E_disp = -Σ_ij f_damp(r_ij) * (s6*C6/r^6 + s8*C8/r^8)
     * BJ damping: f_damp = r^n / (r^n + (a1*sqrt(C8/C6) + a2)^n)
     *
     * Claude Generated (2025): Phase 4 pairwise non-bonded implementation
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_dispersions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread calculating {} dispersion pairs", m_gfnff_dispersions.size()));
    }

    for (int index = 0; index < m_gfnff_dispersions.size(); ++index) {
        const auto& disp = m_gfnff_dispersions[index];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;  // Convert to atomic units if needed

        if (rij > disp.r_cut || rij < 1e-10) continue;  // Skip if beyond cutoff or too close

        // Becke-Johnson damping function (order n=6 for C6, n=8 for C8)
        double r_crit = disp.a1 * std::sqrt(disp.C8 / (disp.C6 + 1e-14)) + disp.a2;

        // C6 term: -s6*C6/r^6 with BJ damping
        double r6 = std::pow(rij, 6);
        double damp6 = std::pow(r_crit, 6);
        double f_damp6 = r6 / (r6 + damp6);
        double E_C6 = -disp.s6 * disp.C6 * f_damp6 / r6;

        // C8 term: -s8*C8/r^8 with BJ damping
        double r8 = r6 * rij * rij;
        double damp8 = std::pow(r_crit, 8);
        double f_damp8 = r8 / (r8 + damp8);
        double E_C8 = -disp.s8 * disp.C8 * f_damp8 / r8;

        double energy = (E_C6 + E_C8) * m_final_factor;
        m_dispersion_energy += energy;

        if (m_calculate_gradient) {
            // Analytical gradient: dE/dr = dE_C6/dr + dE_C8/dr
            // d/dr[f_damp * C_n / r^n] = -n*f_damp*C_n/r^(n+1) + C_n/r^n * df_damp/dr

            // C6 gradient
            double df_damp6_dr = 6.0 * r6 * damp6 / (std::pow(r6 + damp6, 2) * rij);
            double dE_C6_dr = -6.0 * E_C6 / rij + (-disp.s6 * disp.C6 / r6) * df_damp6_dr;

            // C8 gradient
            double df_damp8_dr = 8.0 * r8 * damp8 / (std::pow(r8 + damp8, 2) * rij);
            double dE_C8_dr = -8.0 * E_C8 / rij + (-disp.s8 * disp.C8 / r8) * df_damp8_dr;

            double dEdr = (dE_C6_dr + dE_C8_dr) * m_final_factor;
            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(disp.i) += grad.transpose();
            m_gradient.row(disp.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_dispersions.size() > 0) {
        CurcumaLogger::param("thread_dispersion_energy", fmt::format("{:.6f} Eh", m_dispersion_energy));
    }
}

void ForceFieldThread::CalculateGFNFFRepulsionContribution()
{
    /**
     * @brief GFN-FF Repulsion term (pairwise parallelizable)
     *
     * Reference: Fortran gfnff_engrad.F90:407-439 (bonded repulsion with repscalb=1.7583)
     *                      gfnff_engrad.F90:228 (non-bonded repulsion with repscaln=0.4270)
     * Formula: E_rep = repab * exp(-α*r^β) / r
     * GFN-FF uses β=1.5 (r^1.5 exponent)
     *
     * Claude Generated (2025): Phase 4 pairwise non-bonded implementation
     * Feb 2025: Added bonded/non-bonded scaling factors from Fortran gfnff_param.f90:373-374
     */

    // Scaling factors from Fortran (gfnff_param.f90:373-374)
    static const double REPSCALB = 1.7583;  // Bonded repulsion scaling
    static const double REPSCALN = 0.4270;  // Non-bonded repulsion scaling

    // Build set of bonded pairs for fast lookup
    std::set<std::pair<int, int>> bonded_pairs;
    for (const auto& bond : m_gfnff_bonds) {
        bonded_pairs.insert({bond.i, bond.j});
        bonded_pairs.insert({bond.j, bond.i});  // symmetric
    }

    for (int index = 0; index < m_gfnff_repulsions.size(); ++index) {
        const auto& rep = m_gfnff_repulsions[index];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        if (rij > rep.r_cut || rij < 1e-10) continue;

        // Determine bonded status (used only for reference; scaling already embedded in repab)
        (void)bonded_pairs; // avoid unused warning if not needed elsewhere

        // Base GFN‑FF repulsion energy (without additional scaling factors)
        // E = repab * exp(-α * r^1.5) / r
        double r_1_5 = std::pow(rij, 1.5);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;

        // Apply global energy scaling and final factor
        m_rep_energy += base_energy * m_final_factor * m_rep_scaling;

        if (m_calculate_gradient) {
            // Derivative of base_energy with respect to r
            // dE/dr = -base_energy / r - 1.5 * α * sqrt(r) * base_energy
            double dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
            // Include the same scaling factors as used for the energy
            dEdr *= m_final_factor * m_rep_scaling;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            m_gradient.row(rep.i) += grad.transpose();
            m_gradient.row(rep.j) -= grad.transpose();
        }
    }
}

void ForceFieldThread::CalculateGFNFFCoulombContribution()
{
    /**
     * @brief EEQ-based Coulomb electrostatics with erf damping (pairwise parallelizable)
     *
     * Reference: Phase 3 EEQ charges + Fortran gfnff_engrad.F90:1378-1389
     * Formula: E_coul = q_i * q_j * erf(γ_ij * r_ij) / r_ij
     * Error function damping avoids 1/r singularity at short range
     *
     * Claude Generated (2025): Phase 4 pairwise non-bonded implementation
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_coulombs.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread calculating {} Coulomb pairs", m_gfnff_coulombs.size()));
    }

    for (int index = 0; index < m_gfnff_coulombs.size(); ++index) {
        const auto& coul = m_gfnff_coulombs[index];

        Eigen::Vector3d ri = m_geometry.row(coul.i);
        Eigen::Vector3d rj = m_geometry.row(coul.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        if (rij > coul.r_cut || rij < 1e-10) continue;

        // EEQ Coulomb with erf damping: E = q_i * q_j * erf(γ*r) / r
        double gamma_r = coul.gamma_ij * rij;
        double erf_term = std::erf(gamma_r);
        double energy = coul.q_i * coul.q_j * erf_term / rij;

        m_coulomb_energy += energy * m_final_factor;

        if (m_calculate_gradient) {
            // dE/dr = q_i*q_j * [d/dr(erf(γ*r)/r)]
            //       = q_i*q_j * [γ*exp(-(γ*r)^2)*(2/sqrt(π))/r - erf(γ*r)/r^2]
            const double sqrt_pi = 1.772453850905516;  // sqrt(π)
            double exp_term = std::exp(-gamma_r * gamma_r);
            double derf_dr = coul.gamma_ij * exp_term * (2.0 / sqrt_pi);
            double dEdr = coul.q_i * coul.q_j * (derf_dr / rij - erf_term / (rij * rij)) * m_final_factor;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(coul.i) += grad.transpose();
            m_gradient.row(coul.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_coulombs.size() > 0) {
        CurcumaLogger::param("thread_coulomb_energy", fmt::format("{:.6f} Eh", m_coulomb_energy));
    }
}
