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

#include "src/core/pseudoff.h"

#include "src/capabilities/LevMarDocking.h"

#include <iostream>

#include "docking.h"

Docking::Docking()
{
}

void Docking::PerformDocking()
{
    Molecule result = m_host_structure;

    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;
    std::cout << PseudoFF::LennardJones(m_host_structure, m_guest_structure) << std::endl;

    Geometry destination = GeometryTools::TranslateMolecule(m_guest_structure, m_guest_structure.Centroid(), m_initial_anchor);
    m_guest_structure.setGeometry(destination);

    std::cout << PseudoFF::LennardJones(m_host_structure, m_guest_structure) << std::endl;
    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;

    Position anchor = OptimiseAnchor(&m_host_structure, m_guest_structure, m_initial_anchor);

    destination = GeometryTools::TranslateMolecule(m_guest_structure, m_guest_structure.Centroid(), anchor);
    m_guest_structure.setGeometry(destination);

    std::cout << PseudoFF::LennardJones(m_host_structure, m_guest_structure) << std::endl;
    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;

    for (int i = 0; i < m_guest_structure.AtomCount(); ++i) {
        result.addPair(m_guest_structure.Atom(i));
    }
    m_supramol = result;
}
