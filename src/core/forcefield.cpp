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

void ForceField::setParameter(const json& parameter)
{
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
    for (int i = 0; i < m_threads; ++i) {
        ForceFieldThread* thread = new ForceFieldThread(i, m_threads);
        thread->setGeometry(m_geometry);
        m_threadpool->addThread(thread);
        m_stored_threads.push_back(thread);
        for (int j = int(i * m_bonds.size() / double(m_threads)); j < int((i + 1) * m_bonds.size() / double(m_threads)); ++j)
            thread->addBond(m_bonds[j]);

        for (int j = int(i * m_angles.size() / double(m_threads)); j < int((i + 1) * m_angles.size() / double(m_threads)); ++j)
            thread->addAngle(m_angles[j]);

        for (int j = int(i * m_dihedrals.size() / double(m_threads)); j < int((i + 1) * m_dihedrals.size() / double(m_threads)); ++j)
            thread->addDihedral(m_dihedrals[j]);

        for (int j = int(i * m_inversions.size() / double(m_threads)); j < int((i + 1) * m_inversions.size() / double(m_threads)); ++j)
            thread->addInversion(m_inversions[j]);

        for (int j = int(i * m_vdWs.size() / double(m_threads)); j < int((i + 1) * m_vdWs.size() / double(m_threads)); ++j)
            thread->addvdW(m_vdWs[j]);

        for (int j = int(i * m_EQs.size() / double(m_threads)); j < int((i + 1) * m_EQs.size() / double(m_threads)); ++j)
            thread->addEQ(m_EQs[j]);
    }
}
