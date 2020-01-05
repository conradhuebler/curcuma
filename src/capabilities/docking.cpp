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

/* Junk code


    for(int x = 0; x < 1; ++x)
    {
        Geometry rotation = GeometryTools::RotationX(x*10);
        for(int y = 0; y < 36; ++y)
        {
            rotation *= GeometryTools::RotationY(y*10);

            for(int z = 0; z < 1; ++z)
            {
                rotation *= GeometryTools::RotationZ(z*10).transpose();
                m_guest_structure.setGeometry(destination*rotation);
                m_guest_structure.writeXYZFile(std::to_string(x*10) + "_" + std::to_string(y*10)+ "_" + std::to_string(z*10)+ "_rot.xyz");
            }
        }
    }


 */
void Docking::PerformDocking()
{
    Molecule guest = m_guest_structure;
    Molecule result = m_host_structure;

    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;
    std::cout << PseudoFF::LennardJones(m_host_structure, m_guest_structure) << std::endl;
    Geometry stored_guest = m_guest_structure.getGeometry();
    Position initial_centroid = m_guest_structure.Centroid();

    /*
    Geometry destination = GeometryTools::TranslateMolecule(guest,initial_centroid, m_initial_anchor);

    std::cout << PseudoFF::LennardJones(m_host_structure, m_guest_structure) << std::endl;
    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;
    */

    for (int x = 0; x < 36; ++x) {
        for (int y = 0; y < 36; ++y) {
            for (int z = 0; z < 36; ++z) {
                std::pair<Position, Position> pair = OptimiseAnchor(&m_host_structure, guest, m_initial_anchor, Position{ x * 10, y * 10, z * 10 });
                Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, pair.first, pair.second);
                guest.setGeometry(destination);
                double energy = PseudoFF::LennardJones(m_host_structure, guest);
                m_docking_list.insert(std::pair<double, Vector>(energy, PositionPair2Vector(pair)));
                guest = m_guest_structure;
            }
        }
    }
    guest = m_guest_structure;

    int index = 0;
    for (const auto& element : m_docking_list) {
        if (element.first > 0)
            break;

        std::pair<Position, Position> pair = Vector2PositionPair(element.second);
        std::cout << pair.first.transpose() << " " << pair.second.transpose() << std::endl;
        Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, pair.first, pair.second);
        guest.setGeometry(destination);
        std::cout << guest.Centroid().transpose() << " = Centroid of Guest ### Psudeo FF Energy = " << element.first << std::endl;
        if (index == 0) {
            for (int i = 0; i < guest.AtomCount(); ++i) {
                result.addPair(guest.Atom(i));
            }
            m_supramol = result;
        }
        index++;
    }
}
