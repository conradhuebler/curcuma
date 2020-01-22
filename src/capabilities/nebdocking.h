/*
 * <Prepare and Optimise Structures for Nudge-Elastic-Band calculation. >
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

#include "src/core/global.h"
#include "src/core/molecule.h"

class NEBDocking {
public:
    NEBDocking();

    inline void setStructures(const Molecule& start, const Molecule& end)
    {
        m_start = start;
        m_end = end;
    }
    void Prepare();

    /*! \brief Number of Proton changes allowed */
    inline int ProtonTransfer() const { return m_pt; }

    /*! \brief Set number of allowed proton transfer */
    inline void setProtonTransfer(int pt) { m_pt = pt; }

private:
    std::pair<Molecule, Molecule> DockForNEB(const Molecule& first, const Molecule& second);
    Molecule m_start, m_end;
    int m_pt = 0;
};
