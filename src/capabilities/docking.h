/*
 * <Docking tool for structures. >
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#pragma once

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "src/tools/geometry.h"

#include <map>

class Docking {
public:
    Docking();

    /*! \brief Set the structure of the host molecule */
    void setHostStructure(const Molecule& molecule)
    {
        m_host_structure = molecule;
        m_initial_anchor = m_host_structure.Centroid();
    }

    /*! \brief Set the structure of the guest molecule */
    void setGuestStructure(const Molecule& molecule) { m_guest_structure = molecule; }

    /*! \brief Initial Position, where the guest structure has to be placed */
    void setAnchorPosition(const Position& position) { m_initial_anchor = position; }

    /*! \brief Set the rotation on the x axis */
    void setXRotation(int xxx) { m_xxx_rotation = xxx; }

    /*! \brief Set the rotation on the y axis */
    void setYRotation(int yyy) { m_yyy_rotation = yyy; }

    /*! \brief Set the rotation on the z axis */
    void setZRotation(int zzz) { m_zzz_rotation = zzz; }

    void PerformDocking();

    Molecule getMolecule() const { return m_supramol; }

    void setCheck(bool check) { m_check = check; }

private:
    Molecule m_host_structure, m_guest_structure, m_supramol;
    Position m_initial_anchor = Position{ 0, 0, 0 };
    int m_xxx_rotation = 1, m_yyy_rotation = 1, m_zzz_rotation = 1;
    std::map<double, Vector> m_docking_list;
    std::map<double, Molecule*> m_result_list;
    std::vector<Position> m_anchor_accepted, m_rotation_accepted;
    bool m_check = false;
};
