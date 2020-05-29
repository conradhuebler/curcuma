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

#include "src/capabilities/optimiser/LevMarNEBDocking.h"

#include "src/capabilities/rmsd.h"

#include "nebdocking.h"

NEBDocking::NEBDocking()
{
}

void NEBDocking::Prepare()
{
    /* First reorder structures and contrain proton connectivity */
    RMSDDriver* driver = new RMSDDriver();
    driver->setReference(m_start);
    driver->setTarget(m_end);
    driver->setForceReorder(true);
    driver->setCheckConnections(true);
    driver->setProtonTransfer(ProtonTransfer());
    driver->AutoPilot();

    Molecule end = driver->TargetAligned();
    end.writeXYZFile("neb_pre_ende.xyz");

    Molecule first, second;
    std::pair<Molecule, Molecule> result;

    result = DockForNEB(m_start, end);

    driver->setReference(result.first);
    driver->setTarget(result.second);
    driver->setForceReorder(false);

    driver->AutoPilot();
    end = driver->TargetAligned();
    end.writeXYZFile("neb_ende.xyz");
}

std::pair<Molecule, Molecule> NEBDocking::DockForNEB(const Molecule& first, const Molecule& second)
{
    std::pair<Molecule, Molecule> result;
    result = OptimiseAtoms(first, second);

    result.first.writeXYZFile("neb_first.xyz");
    result.second.writeXYZFile("neb_second.xyz");

    return result;
}
