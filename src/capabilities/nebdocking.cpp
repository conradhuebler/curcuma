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

#include "src/capabilities/LevMarNEBDocking.h"

#include "src/capabilities/rmsd.h"

#include "nebdocking.h"

NEBDocking::NEBDocking()
{
}

void NEBDocking::Prepare()
{
    /* First reorder structures and contrain proton connectivity */
    RMSDDriver* driver = new RMSDDriver(m_start, m_end);
    driver->setForceReorder(true);
    driver->setCheckConnections(true);
    driver->setProtonTransfer(ProtonTransfer());
    driver->AutoPilot();

    Molecule end = driver->TargetAligned();
    Molecule first, second;
    if (m_start.GetFragments().size() >= m_end.GetFragments().size()) {
        second = DockForNEB(m_start, end);
    } else {
        second = DockForNEB(end, m_start);
    }
}

Molecule NEBDocking::DockForNEB(const Molecule& first, const Molecule& second)
{
    Molecule result;
    result = OptimiseAtoms(&first, second);
    result.writeXYZFile("neb_ende.xyz");
    return result;
}
