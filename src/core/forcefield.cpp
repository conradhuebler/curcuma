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
#include "forcefieldthread.h"

ForceField::ForceField(const json& controller)
{
    json parameter = MergeJson(UFFParameterJson, controller);

    m_threadpool = new CxxThreadPool();
    m_threadpool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    m_threads = parameter["threads"];
    m_gradient_type = parameter["gradient"];
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

void ForceField::UpdateGeometry(const std::vector<std::array<double, 3>>& geometry)
{
#pragma message("replace with raw data")
    for (int i = 0; i < m_natoms; ++i) {
        m_geometry(i, 0) = geometry[i][0];
        m_geometry(i, 1) = geometry[i][1];
        m_geometry(i, 2) = geometry[i][2];
    }
}

void ForceField::setParameter(const json& parameters)
{
    if (parameters.contains("bonds"))
        setBonds(parameters["bonds"]);
    if (parameters.contains("angles"))
        setAngles(parameters["angles"]);
    if (parameters.contains("dihedrals"))
        setDihedrals(parameters["dihedrals"]);
    if (parameters.contains("inversions"))
        setInversions(parameters["inversions"]);
    if (parameters.contains("vdws"))
        setvdWs(parameters["vdws"]);
    m_parameters = parameters;
    m_method = m_parameters["method"];
    if (m_parameters.contains("e0"))
        m_e0 = m_parameters["e0"];
    AutoRanges();
}

void ForceField::setBonds(const json& bonds)
{
    m_bonds.clear();
    for (int i = 0; i < bonds.size(); ++i) {
        json bond = bonds[i].get<json>();
        Bond b;
        b.type = bond["type"];
        b.i = bond["i"];
        b.j = bond["j"];
        b.k = bond["k"];
        b.distance = bond["distance"];
        b.exponent = bond["exponent"];

        b.r0_ij = bond["r0_ij"];
        b.r0_ik = bond["r0_ik"];

        b.fc = bond["fc"];
        m_bonds.push_back(b);
    }
}

void ForceField::setAngles(const json& angles)
{
    m_angles.clear();
    for (int i = 0; i < angles.size(); ++i) {
        json angle = angles[i].get<json>();
        Angle a;

        a.type = angle["type"];

        a.i = angle["i"];
        a.j = angle["j"];
        a.k = angle["k"];
        a.C0 = angle["C0"];
        a.C1 = angle["C1"];
        a.C2 = angle["C2"];
        a.fc = angle["fc"];
        a.r0_ij = angle["r0_ij"];
        a.r0_ik = angle["r0_ik"];
        a.theta0_ijk = angle["theta0_ijk"];
        m_angles.push_back(a);
    }
}

void ForceField::setDihedrals(const json& dihedrals)
{
    m_dihedrals.clear();
    for (int i = 0; i < dihedrals.size(); ++i) {
        json dihedral = dihedrals[i].get<json>();
        Dihedral d;
        d.type = dihedral["type"];

        d.i = dihedral["i"];
        d.j = dihedral["j"];
        d.k = dihedral["k"];
        d.l = dihedral["l"];
        d.V = dihedral["V"];
        d.n = dihedral["n"];
        d.phi0 = dihedral["phi0"];
        m_dihedrals.push_back(d);
    }
}

void ForceField::setInversions(const json& inversions)
{
    m_inversions.clear();
    for (int i = 0; i < inversions.size(); ++i) {
        json inversion = inversions[i].get<json>();
        Inversion inv;
        inv.type = inversion["type"];

        inv.i = inversion["i"];
        inv.j = inversion["j"];
        inv.k = inversion["k"];
        inv.l = inversion["l"];
        inv.fc = inversion["fc"];
        inv.C0 = inversion["C0"];
        inv.C1 = inversion["C1"];
        inv.C2 = inversion["C2"];

        m_inversions.push_back(inv);
    }
}

void ForceField::setvdWs(const json& vdws)
{
    m_vdWs.clear();
    for (int i = 0; i < vdws.size(); ++i) {
        json vdw = vdws[i].get<json>();
        vdW v;
        v.type = vdw["type"];

        v.i = vdw["i"];
        v.j = vdw["j"];
        v.C_ij = vdw["C_ij"];
        v.r0_ij = vdw["r0_ij"];

        m_vdWs.push_back(v);
    }
}

void ForceField::AutoRanges()
{
    int free_threads = m_threads;
    int d3 = m_parameters["d3"];
    if (d3) {
        if (free_threads > 1)
            free_threads--;
        D3Thread* thread = new D3Thread(m_threads - 1, free_threads);
        thread->setParamater(m_parameters);
        thread->Initialise(m_atom_types);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
    }
    int h4 = m_parameters["h4"];
    if (h4) {
        if (free_threads > 1)
            free_threads--;
        H4Thread* thread = new H4Thread(m_threads - 1, free_threads);
        thread->setParamater(m_parameters);
        thread->Initialise(m_atom_types);
        if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) {
            thread->setMethod(1);
        } else if (std::find(m_qmdff_methods.begin(), m_qmdff_methods.end(), m_method) != m_qmdff_methods.end()) {
            thread->setMethod(2);
        }
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
    }
    for (int i = 0; i < free_threads; ++i) {
        ForceFieldThread* thread = new ForceFieldThread(i, free_threads);
        thread->setGeometry(m_geometry, false);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
        thread->setMethod(2);
        for (int j = int(i * m_bonds.size() / double(free_threads)); j < int((i + 1) * m_bonds.size() / double(free_threads)); ++j)
            thread->addBond(m_bonds[j]);

        for (int j = int(i * m_angles.size() / double(free_threads)); j < int((i + 1) * m_angles.size() / double(free_threads)); ++j)
            thread->addAngle(m_angles[j]);

        for (int j = int(i * m_dihedrals.size() / double(free_threads)); j < int((i + 1) * m_dihedrals.size() / double(free_threads)); ++j)
            thread->addDihedral(m_dihedrals[j]);

        for (int j = int(i * m_inversions.size() / double(free_threads)); j < int((i + 1) * m_inversions.size() / double(free_threads)); ++j)
            thread->addInversion(m_inversions[j]);

        for (int j = int(i * m_vdWs.size() / double(free_threads)); j < int((i + 1) * m_vdWs.size() / double(free_threads)); ++j)
            thread->addvdW(m_vdWs[j]);

        for (int j = int(i * m_EQs.size() / double(free_threads)); j < int((i + 1) * m_EQs.size() / double(free_threads)); ++j)
            thread->addEQ(m_EQs[j]);
    }
}

Eigen::MatrixXd ForceField::NumGrad()
{
    Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(m_natoms, 3);

    double dx = 1e-6; // m_d;
    // bool g = m_CalculateGradient;
    // m_CalculateGradient = false;
    double E1, E2;
    for (int i = 0; i < m_natoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            m_geometry(i, j) += dx;
            E1 = Calculate(false, false);
            m_geometry(i, j) -= 2 * dx;
            E2 = Calculate(false, false);
            gradient(i, j) = (E1 - E2) / (2 * dx);
            m_geometry(i, j) += dx;
        }
    }
    // m_CalculateGradient = g;
    return gradient;
}

double ForceField::Calculate(bool gradient, bool verbose)
{
    m_gradient = Eigen::MatrixXd::Zero(m_geometry.rows(), 3);
    double energy = 0.0;
    double d4_energy = 0;
    double d3_energy = 0;
    double bond_energy = 0.0;
    double angle_energy = 0.0;
    double dihedral_energy = 0.0;
    double inversion_energy = 0.0;
    double vdw_energy = 0.0;
    double rep_energy = 0.0;
    double eq_energy = 0.0;
    double h4_energy = 0.0;
    double hh_energy = 0.0;

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        m_stored_threads[i]->UpdateGeometry(m_geometry, gradient);
    }

    m_threadpool->Reset();
    m_threadpool->setActiveThreadCount(m_threads);

    m_threadpool->StartAndWait();
    m_threadpool->setWakeUp(m_threadpool->WakeUp() / 2);

    for (int i = 0; i < m_stored_threads.size(); ++i) {
        bond_energy += m_stored_threads[i]->BondEnergy();
        angle_energy += m_stored_threads[i]->AngleEnergy();
        dihedral_energy += m_stored_threads[i]->DihedralEnergy();
        inversion_energy += m_stored_threads[i]->InversionEnergy();
        if (m_stored_threads[i]->Type() != 3)
            vdw_energy += m_stored_threads[i]->VdWEnergy();
        else
            h4_energy += m_stored_threads[i]->VdWEnergy();
        if (m_stored_threads[i]->Type() != 3)
            rep_energy += m_stored_threads[i]->RepEnergy();
        else
            hh_energy += m_stored_threads[i]->RepEnergy();
        // eq_energy += m_stored_threads[i]->RepEnergy();

        m_gradient += m_stored_threads[i]->Gradient();
    }

    energy = m_e0 + bond_energy + angle_energy + dihedral_energy + inversion_energy + vdw_energy + rep_energy + eq_energy + h4_energy + hh_energy;
    if (verbose) {
        std::cout << "Total energy " << energy << " Eh. Sum of " << std::endl
                  << "E0 (from QMDFF) " << m_e0 << " Eh" << std::endl
                  << "Bond Energy " << bond_energy << " Eh" << std::endl
                  << "Angle Energy " << angle_energy << " Eh" << std::endl
                  << "Dihedral Energy " << dihedral_energy << " Eh" << std::endl
                  << "Inversion Energy " << inversion_energy << " Eh" << std::endl
                  << "Nonbonded Energy " << vdw_energy + rep_energy << " Eh" << std::endl
                  << "D3 Energy " << d3_energy << " Eh" << std::endl
                  << "D4 Energy " << d4_energy << " Eh" << std::endl
                  << "HBondCorrection " << h4_energy << " Eh" << std::endl
                  << "HHRepCorrection " << hh_energy << " Eh" << std::endl
                  << std::endl;

        // for (int i = 0; i < m_natoms; ++i) {
        //     std::cout << m_gradient(i, 0) << " " << m_gradient(i, 1) << " " << m_gradient(i, 2) << std::endl;
        // }
    }
    return energy;
}
