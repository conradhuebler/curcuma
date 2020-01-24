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

#include "src/capabilities/optimiser/LevMarDocking.h"

#include <iostream>

#include "docking.h"

Docking::Docking()
{
}


void Docking::PerformDocking()
{
    Molecule guest = m_guest_structure;
    Molecule result = m_host_structure;

    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;

    Geometry stored_guest = m_guest_structure.getGeometry();
    Position initial_centroid = m_guest_structure.Centroid();

    /*
    Geometry destination = GeometryTools::TranslateMolecule(guest,initial_centroid, m_initial_anchor);

    std::cout << PseudoFF::LennardJones(m_host_structure, m_guest_structure) << std::endl;
    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;
    */
    std::cout << std::endl
              << "** Docking Phase 0 - Starting **" << std::endl
              << std::endl;
    int X = 36;
    int Y = 36;
    int Z = 36;

    for (int x = 0; x < X; ++x) {
        for (int y = 0; y < Y; ++y) {
            for (int z = 0; z < Z; ++z) {
                Molecule* molecule = new Molecule(m_host_structure);
                guest = m_guest_structure;

                std::pair<Position, Position> pair = OptimiseAnchor(&m_host_structure, guest, m_initial_anchor, Position{ x * 10, y * 10, z * 10 });

                bool accept = true;

                if (GeometryTools::Distance(m_initial_anchor, pair.first) > 1e5) {
                    accept = false;
                    continue;
                }

                for (std::size_t i = 0; i < m_anchor_accepted.size(); ++i) {
                    Position anchor = m_anchor_accepted[i];
                    Position rotation = m_rotation_accepted[i];
                    if (GeometryTools::Distance(anchor, pair.first) < 1e-1 || GeometryTools::Distance(rotation, pair.second) < 1e-1)
                        accept = false;
                    //std::cout << GeometryTools::Distance(anchor, pair.first) << " --- "  << GeometryTools::Distance(rotation, pair.second) << std::endl;
                }
                if (accept == false)
                    continue;
                m_anchor_accepted.push_back(pair.first);
                m_rotation_accepted.push_back(pair.second);

                Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, pair.first, pair.second);

                // std::cout << pair.first.transpose() << " " << pair.second.transpose() << std::endl;

                guest.setGeometry(destination);
                double energy = PseudoFF::LennardJones(m_host_structure, guest);
                m_docking_list.insert(std::pair<double, Vector>(energy, PositionPair2Vector(pair)));
                result = m_host_structure;
                for (std::size_t i = 0; i < guest.AtomCount(); ++i) {
                    molecule->addPair(guest.Atom(i));
                }
                m_result_list.insert(std::pair<double, Molecule*>(energy, molecule));
            }
        }
        std::cout << (x / double(X)) * 100 << "% - " << m_anchor_accepted.size() << " stored structures." << std::endl;
    }
    guest = m_guest_structure;
    std::cout << std::endl
              << "** Docking Phase 0 - Finished **" << std::endl;
    int index = 0;

    for (const auto& pair : m_result_list) {
        ++index;

        std::cout << pair.first << std::endl;
        pair.second->appendXYZFile("docked_structures_" + std::to_string(pair.second->GetFragments(1.3).size()) + ".xyz");
        pair.second->print_geom();
        delete pair.second;
    }
}
