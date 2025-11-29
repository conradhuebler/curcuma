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
#include "gfnff_par.h"

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

        // GFN-FF hydrogen bond and halogen bond terms (Phase 5)
        CalculateGFNFFHydrogenBondContribution();  // HB three-body terms
        CalculateGFNFFHalogenBondContribution();   // XB three-body terms

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

// Phase 4: GFN-FF hydrogen bond and halogen bond addition methods (Claude Generated 2025)

void ForceFieldThread::addGFNFFHydrogenBond(const GFNFFHydrogenBond& hbond)
{
    m_gfnff_hbonds.push_back(hbond);
}

void ForceFieldThread::addGFNFFHalogenBond(const GFNFFHalogenBond& xbond)
{
    m_gfnff_xbonds.push_back(xbond);
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
    // Claude Generated (2025): Now using correct GFNFF_Geometry::calculateDihedralAngle

    for (int index = 0; index < m_gfnff_dihedrals.size(); ++index) {
        const auto& dihedral = m_gfnff_dihedrals[index];

        // Extract atom positions as Eigen::Vector3d
        Eigen::Vector3d r_i = m_geometry.row(dihedral.i).head<3>();
        Eigen::Vector3d r_j = m_geometry.row(dihedral.j).head<3>();
        Eigen::Vector3d r_k = m_geometry.row(dihedral.k).head<3>();
        Eigen::Vector3d r_l = m_geometry.row(dihedral.l).head<3>();

        // Use GFN-FF dihedral angle calculation (not UFF!)
        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);

        double V = dihedral.V;
        double n = dihedral.n;
        double phi0 = dihedral.phi0;

        // GFN-FF energy formula
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
    // Claude Generated (2025): Now using correct GFNFF_Geometry::calculateOutOfPlaneAngle

    for (int index = 0; index < m_gfnff_inversions.size(); ++index) {
        const auto& inversion = m_gfnff_inversions[index];

        // Extract atom positions as Eigen::Vector3d
        Eigen::Vector3d r_i = m_geometry.row(inversion.i).head<3>();  // out-of-plane atom
        Eigen::Vector3d r_j = m_geometry.row(inversion.j).head<3>();  // plane atom 1
        Eigen::Vector3d r_k = m_geometry.row(inversion.k).head<3>();  // plane atom 2
        Eigen::Vector3d r_l = m_geometry.row(inversion.l).head<3>();  // central atom

        // Use GFN-FF out-of-plane angle calculation (not UFF!)
        Matrix derivate;
        double theta = GFNFF_Geometry::calculateOutOfPlaneAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);

        // GFN-FF inversion energy formula
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

// ============================================================================
// Phase 5: GFN-FF HB/XB Damping Helper Functions
// ============================================================================
// Claude Generated (2025): Direct translation from gfnff_engrad.F90
// Reference: abhgfnff_eg1(), rbxgfnff_eg()

/**
 * @brief Out-of-line damping for HB/XB geometry
 *
 * Penalizes non-linear A-H...B or A-X...B geometries
 * Formula: outl = 2 / (1 + exp(bacut/r_AB × (r_AH + r_HB)/r_AB - 1))
 *
 * @param r_AH Distance A-H (or A-X) in Bohr
 * @param r_HB Distance H-B (or X-B) in Bohr
 * @param r_AB Distance A-B in Bohr
 * @param bacut Angle cut-off parameter (HB_BACUT or XB_BACUT)
 */
inline double damping_out_of_line(double r_AH, double r_HB, double r_AB, double bacut)
{
    double ratio = (r_AH + r_HB) / r_AB;
    double exponent = bacut / r_AB * (ratio - 1.0);
    return 2.0 / (1.0 + std::exp(exponent));
}

/**
 * @brief Short-range damping for HB/XB
 *
 * Prevents unphysical close-contact overbinding
 * Formula: damp_s = 1 / (1 + (scut × r_vdw / r²)^alp)
 *
 * @param r Distance between atoms in Bohr
 * @param r_vdw Combined van der Waals radius in Bohr
 * @param scut Short-range cut-off parameter (HB_SCUT or XB_SCUT)
 * @param alp Damping exponent (HB_ALP)
 */
inline double damping_short_range(double r, double r_vdw, double scut, double alp)
{
    double ratio = scut * r_vdw / (r * r);
    return 1.0 / (1.0 + std::pow(ratio, alp));
}

/**
 * @brief Long-range damping for HB/XB
 *
 * Smooth cut-off at large distances
 * Formula: damp_l = 1 / (1 + (r² / longcut)^alp)
 *
 * @param r Distance between atoms in Bohr
 * @param longcut Long-range cut-off (HB_LONGCUT or HB_LONGCUT_XB)
 * @param alp Damping exponent (HB_ALP)
 */
inline double damping_long_range(double r, double longcut, double alp)
{
    double r_sq = r * r;
    return 1.0 / (1.0 + std::pow(r_sq / longcut, alp));
}

/**
 * @brief Charge-dependent scaling factor
 *
 * EEQ charge modulation for HB/XB strength
 * Formula: Q = exp(st × q) / (exp(st × q) + sf)
 *
 * @param q EEQ charge on atom
 * @param st Scaling strength (HB_ST or XB_ST)
 * @param sf Scaling offset (HB_SF or XB_SF)
 */
inline double charge_scaling(double q, double st, double sf)
{
    double exp_term = std::exp(st * q);
    return exp_term / (exp_term + sf);
}

// ============================================================================
// Phase 5: GFN-FF HB/XB Energy Calculation Methods
// ============================================================================

void ForceFieldThread::CalculateGFNFFHydrogenBondContribution()
{
    /**
     * @brief Calculate GFN-FF hydrogen bond energy contribution
     *
     * Reference: gfnff_engrad.F90 - abhgfnff_eg1() (Case 1) and abhgfnff_eg2new() (Case 2)
     * Formula: E_HB = B_AH × C_AH × C_B × (-C_acidity × R_damp × Q_H^outl)
     *
     * Claude Generated (2025): Direct translation from Fortran implementation
     */

    using namespace GFNFFParameters;

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_hbonds.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} hydrogen bonds", m_thread, m_gfnff_hbonds.size()));
    }

    for (const auto& hb : m_gfnff_hbonds) {
        // Get atom positions
        Eigen::Vector3d pos_A = m_geometry.row(hb.i).transpose();
        Eigen::Vector3d pos_H = m_geometry.row(hb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(hb.k).transpose();

        // Calculate distance vectors
        Eigen::Vector3d r_AH_vec = pos_H - pos_A;
        Eigen::Vector3d r_HB_vec = pos_B - pos_H;
        Eigen::Vector3d r_AB_vec = pos_B - pos_A;

        // Calculate distances (convert to Bohr via m_au)
        double r_AH = r_AH_vec.norm() * m_au;
        double r_HB = r_HB_vec.norm() * m_au;
        double r_AB = r_AB_vec.norm() * m_au;

        // Distance cutoff check
        if (r_AB > hb.r_cut) continue;

        // --- Donor-Acceptor Strength Term B_AH ---
        // Reference: gfnff_engrad.F90:1445-1448
        // B = (C_A × r_AH^4 + C_B × r_HB^4) / (r_AH^4 + r_HB^4)
        double r_AH_4 = r_AH * r_AH * r_AH * r_AH;
        double r_HB_4 = r_HB * r_HB * r_HB * r_HB;
        double B_AH = (hb.basicity_A * r_AH_4 + hb.basicity_B * r_HB_4) / (r_AH_4 + r_HB_4);

        // --- Charge Scaling Factors ---
        double Q_H = charge_scaling(hb.q_H, HB_ST, HB_SF);
        double Q_A = charge_scaling(hb.q_A, HB_ST, HB_SF);
        double Q_B = charge_scaling(hb.q_B, HB_ST, HB_SF);

        // --- Combined Damping R_damp ---
        // Use covalent radii as approximation for vdW radius
        double r_vdw_AB = covalent_radii[hb.i] + covalent_radii[hb.k];  // Angström -> need conversion
        r_vdw_AB *= 1.88973;  // Convert Å to Bohr

        double damp_short = damping_short_range(r_AB, r_vdw_AB, HB_SCUT, HB_ALP);
        double damp_long = damping_long_range(r_AB, HB_LONGCUT, HB_ALP);
        double damp_outl = damping_out_of_line(r_AH, r_HB, r_AB, HB_BACUT);

        double R_damp = damp_short * damp_long * damp_outl / (r_AB * r_AB * r_AB);

        // --- Acidity Term C_acidity ---
        // Reference: hbonds() subroutine - complex mixing of basicity/acidity
        double C_acidity = hb.acidity_A;  // Simplified for Case 1

        // --- Total HB Energy ---
        // Note: Negative sign convention from Fortran
        double E_HB = B_AH * Q_A * Q_B * (-C_acidity) * R_damp * Q_H;

        m_energy_hbond += E_HB * m_final_factor;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  HB({}-{}-{}): r_AB={:.3f} Bohr, E={:.6f} Eh",
                hb.i, hb.j, hb.k, r_AB, E_HB * m_final_factor));
        }

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Reference: gfnff_engrad.F90 - abhgfnff_eg1() lines 1501-1557
        // Claude Generated (2025): Phase 6 - HB analytical gradients
        if (m_calculate_gradient) {
            // --- 1. Compute Energy Derivative Components ---

            // ∂B_AH/∂r_AH and ∂B_AH/∂r_HB (donor-acceptor mixing term derivatives)
            double r_AH_3 = r_AH * r_AH * r_AH;
            double r_HB_3 = r_HB * r_HB * r_HB;
            double denominator = r_AH_4 + r_HB_4;
            double denominator_sq = denominator * denominator;

            double dB_drAH = 4.0 * r_AH_3 * (hb.basicity_A * r_HB_4 - hb.basicity_B * r_AH_4) / denominator_sq;
            double dB_drHB = 4.0 * r_HB_3 * (hb.basicity_B * r_AH_4 - hb.basicity_A * r_HB_4) / denominator_sq;

            // ∂(damping)/∂r derivatives
            // Short-range damping: damp_s = 1 / (1 + (scut × r_vdw / r²)^alp)
            double ratio_short = HB_SCUT * r_vdw_AB / (r_AB * r_AB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term / (r_AB * (1.0 + damp_short_term));

            // Long-range damping: damp_l = 1 / (1 + (r² / longcut)^alp)
            double ratio_long = (r_AB * r_AB) / HB_LONGCUT;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_AB * damp_long * damp_long_term / (HB_LONGCUT * (1.0 + damp_long_term));

            // Out-of-line damping: outl = 2 / (1 + exp(bacut/r_AB × (r_AH + r_HB)/r_AB - 1))
            double ratio = (r_AH + r_HB) / r_AB;
            double exponent = HB_BACUT / r_AB * (ratio - 1.0);
            double exp_term = std::exp(exponent);
            double denom_outl = 1.0 + exp_term;

            // Derivatives of out-of-line damping w.r.t. each distance
            double ddamp_outl_drAH = -2.0 * exp_term * HB_BACUT / (r_AB * r_AB * denom_outl * denom_outl);
            double ddamp_outl_drHB = ddamp_outl_drAH;  // Same formula
            double ddamp_outl_drAB = 2.0 * exp_term * HB_BACUT * (r_AH + r_HB) / (r_AB * r_AB * r_AB * denom_outl * denom_outl);

            // Combined damping derivative: R_damp = damp_short * damp_long * damp_outl / r_AB³
            double R_damp_no_r3 = damp_short * damp_long * damp_outl;

            double dRdamp_drAB = (ddamp_short_dr * damp_long * damp_outl +
                                  damp_short * ddamp_long_dr * damp_outl +
                                  damp_short * damp_long * ddamp_outl_drAB) / (r_AB * r_AB * r_AB)
                                - 3.0 * R_damp / r_AB;

            double dRdamp_drAH = damp_short * damp_long * ddamp_outl_drAH / (r_AB * r_AB * r_AB);
            double dRdamp_drHB = damp_short * damp_long * ddamp_outl_drHB / (r_AB * r_AB * r_AB);

            // --- 2. Energy Derivatives w.r.t. Distances ---
            // E_HB = B_AH × Q_A × Q_B × (-C_acidity) × R_damp × Q_H
            double E_prefactor = Q_A * Q_B * (-C_acidity) * Q_H;

            double dE_drAH = E_prefactor * (dB_drAH * R_damp + B_AH * dRdamp_drAH);
            double dE_drHB = E_prefactor * (dB_drHB * R_damp + B_AH * dRdamp_drHB);
            double dE_drAB = E_prefactor * B_AH * dRdamp_drAB;

            // --- 3. Chain Rule: Position Vector Derivatives ---
            // ∇_atom r_ij = (r_j - r_i) / |r_j - r_i|
            Eigen::Vector3d grad_rAH_unit = r_AH_vec / r_AH;  // Direction A → H
            Eigen::Vector3d grad_rHB_unit = r_HB_vec / r_HB;  // Direction H → B
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / r_AB;  // Direction A → B

            // --- 4. Accumulate Gradients on Each Atom ---
            // Atom A: influenced by r_AH (negative direction) and r_AB (negative direction)
            Eigen::Vector3d grad_A = (-dE_drAH * grad_rAH_unit - dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(hb.i) += grad_A.transpose();

            // Atom H: influenced by r_AH (positive direction) and r_HB (negative direction)
            Eigen::Vector3d grad_H = (dE_drAH * grad_rAH_unit - dE_drHB * grad_rHB_unit) * m_final_factor;
            m_gradient.row(hb.j) += grad_H.transpose();

            // Atom B: influenced by r_HB (positive direction) and r_AB (positive direction)
            Eigen::Vector3d grad_B = (dE_drHB * grad_rHB_unit + dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(hb.k) += grad_B.transpose();

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "    HB Gradient: |∇A|={:.4f} |∇H|={:.4f} |∇B|={:.4f} Eh/Bohr",
                    grad_A.norm(), grad_H.norm(), grad_B.norm()));
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_hbonds.size() > 0) {
        CurcumaLogger::param("thread_hbond_energy", fmt::format("{:.6f} Eh", m_energy_hbond));
    }
}

void ForceFieldThread::CalculateGFNFFHalogenBondContribution()
{
    /**
     * @brief Calculate GFN-FF halogen bond energy contribution
     *
     * Reference: gfnff_engrad.F90 - rbxgfnff_eg() (lines 2928-3043)
     * Formula: E_XB = -R_damp × Q_outl × C_B × Q_B × C_X × Q_X
     *
     * Claude Generated (2025): XB uses different damping parameters than HB
     */

    using namespace GFNFFParameters;

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_xbonds.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} halogen bonds", m_thread, m_gfnff_xbonds.size()));
    }

    for (const auto& xb : m_gfnff_xbonds) {
        // Get atom positions
        Eigen::Vector3d pos_A = m_geometry.row(xb.i).transpose();
        Eigen::Vector3d pos_X = m_geometry.row(xb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(xb.k).transpose();

        // Calculate distance vectors
        Eigen::Vector3d r_AX_vec = pos_X - pos_A;
        Eigen::Vector3d r_XB_vec = pos_B - pos_X;
        Eigen::Vector3d r_AB_vec = pos_B - pos_A;

        // Calculate distances (convert to Bohr via m_au)
        double r_AX = r_AX_vec.norm() * m_au;
        double r_XB = r_XB_vec.norm() * m_au;
        double r_AB = r_AB_vec.norm() * m_au;

        // Distance cutoff check
        if (r_XB > xb.r_cut) continue;

        // --- Charge Scaling Factors ---
        // XB uses different parameters: XB_ST=15, XB_SF=0.03 (much weaker offset)
        double Q_X = charge_scaling(xb.q_X, XB_ST, XB_SF);
        double Q_B = charge_scaling(xb.q_B, XB_ST, XB_SF);

        // --- Combined Damping R_damp ---
        double r_vdw_XB = covalent_radii[xb.j] + covalent_radii[xb.k];  // Angström
        r_vdw_XB *= 1.88973;  // Convert Å to Bohr

        // XB uses different cutoff parameters
        double damp_short = damping_short_range(r_XB, r_vdw_XB, XB_SCUT, HB_ALP);
        double damp_long = damping_long_range(r_XB, HB_LONGCUT_XB, HB_ALP);
        double damp_outl = damping_out_of_line(r_AX, r_XB, r_AB, XB_BACUT);

        double R_damp = damp_short * damp_long * damp_outl / (r_XB * r_XB * r_XB);

        // --- Total XB Energy ---
        // C_B = 1.0 (fixed), C_X = xbaci parameter
        double C_B = 1.0;
        double C_X = xb.acidity_X;

        double E_XB = -R_damp * C_B * Q_B * C_X * Q_X;

        m_energy_xbond += E_XB * m_final_factor;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  XB({}-{}-{}): r_XB={:.3f} Bohr, E={:.6f} Eh",
                xb.i, xb.j, xb.k, r_XB, E_XB * m_final_factor));
        }

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Reference: gfnff_engrad.F90 - rbxgfnff_eg() gradient section
        // Claude Generated (2025): Phase 6 - XB analytical gradients
        if (m_calculate_gradient) {
            // XB gradients are simpler than HB (no donor-acceptor mixing term)

            // --- 1. Compute Damping Derivatives ---
            // Short-range damping: damp_s = 1 / (1 + (scut × r_vdw / r²)^alp)
            double ratio_short = XB_SCUT * r_vdw_XB / (r_XB * r_XB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term / (r_XB * (1.0 + damp_short_term));

            // Long-range damping: damp_l = 1 / (1 + (r² / longcut)^alp)
            double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_XB * damp_long * damp_long_term / (HB_LONGCUT_XB * (1.0 + damp_long_term));

            // Out-of-line damping derivatives
            double ratio = (r_AX + r_XB) / r_AB;
            double exponent = XB_BACUT / r_AB * (ratio - 1.0);
            double exp_term = std::exp(exponent);
            double denom_outl = 1.0 + exp_term;

            // Derivatives of out-of-line damping w.r.t. each distance
            double ddamp_outl_drAX = -2.0 * exp_term * XB_BACUT / (r_AB * r_AB * denom_outl * denom_outl);
            double ddamp_outl_drXB = ddamp_outl_drAX;  // Same formula
            double ddamp_outl_drAB = 2.0 * exp_term * XB_BACUT * (r_AX + r_XB) / (r_AB * r_AB * r_AB * denom_outl * denom_outl);

            // Combined damping derivative: R_damp = damp_short * damp_long * damp_outl / r_XB³
            // Note: For XB, primary distance is r_XB, not r_AB as in HB

            double dRdamp_drXB = (ddamp_short_dr * damp_long * damp_outl +
                                  damp_short * ddamp_long_dr * damp_outl +
                                  damp_short * damp_long * ddamp_outl_drXB) / (r_XB * r_XB * r_XB)
                                - 3.0 * R_damp / r_XB;

            double dRdamp_drAX = damp_short * damp_long * ddamp_outl_drAX / (r_XB * r_XB * r_XB);
            double dRdamp_drAB = damp_short * damp_long * ddamp_outl_drAB / (r_XB * r_XB * r_XB);

            // --- 2. Energy Derivatives w.r.t. Distances ---
            // E_XB = -R_damp × C_B × Q_B × C_X × Q_X
            double E_prefactor = -C_B * Q_B * C_X * Q_X;

            double dE_drXB = E_prefactor * dRdamp_drXB;
            double dE_drAX = E_prefactor * dRdamp_drAX;
            double dE_drAB = E_prefactor * dRdamp_drAB;

            // --- 3. Chain Rule: Position Vector Derivatives ---
            Eigen::Vector3d grad_rAX_unit = r_AX_vec / r_AX;  // Direction A → X
            Eigen::Vector3d grad_rXB_unit = r_XB_vec / r_XB;  // Direction X → B
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / r_AB;  // Direction A → B

            // --- 4. Accumulate Gradients on Each Atom ---
            // Atom A: influenced by r_AX (negative direction) and r_AB (negative direction)
            Eigen::Vector3d grad_A = (-dE_drAX * grad_rAX_unit - dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(xb.i) += grad_A.transpose();

            // Atom X: influenced by r_AX (positive direction) and r_XB (negative direction)
            Eigen::Vector3d grad_X = (dE_drAX * grad_rAX_unit - dE_drXB * grad_rXB_unit) * m_final_factor;
            m_gradient.row(xb.j) += grad_X.transpose();

            // Atom B: influenced by r_XB (positive direction) and r_AB (positive direction)
            Eigen::Vector3d grad_B = (dE_drXB * grad_rXB_unit + dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(xb.k) += grad_B.transpose();

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "    XB Gradient: |∇A|={:.4f} |∇X|={:.4f} |∇B|={:.4f} Eh/Bohr",
                    grad_A.norm(), grad_X.norm(), grad_B.norm()));
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_xbonds.size() > 0) {
        CurcumaLogger::param("thread_xbond_energy", fmt::format("{:.6f} Eh", m_energy_xbond));
    }
}
