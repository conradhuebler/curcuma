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

#include "src/tools/general.h"

#include "src/capabilities/optimiser/LBFGSInterface.h"

#include "src/capabilities/optimiser/LevMarDocking.h"
#include "src/capabilities/optimiser/XTBDocking.h"

#include <algorithm>
#include <future>
#include <iostream>
#include <thread>

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

    Geometry geometry = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, m_initial_anchor, Position{ 0, 0, 0 });
    /*
    Molecule temp_guest = guest;
    temp_guest.setGeometry(geometry);
    Geometry geom2 = PrepareHost(&m_host_structure, &temp_guest);
    Molecule host2 = m_host_structure;
    host2.setGeometry(geom2);


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

    int max = 10;
    if (m_check) {
        X = 10;
        Y = 10;
        Z = 10;
        max = 36;
        std::cout << "Fewer sampling " << std::endl;
    }
    int excluded = 0, all = 0;
    for (int x = 0; x < X; ++x) {
        for (int y = 0; y < Y; ++y) {
            for (int z = 0; z < Z; ++z) {
                if (m_result_list.size() >= 5)
                    continue;
                ++all;
                Molecule* molecule = new Molecule(m_host_structure);
                guest = m_guest_structure;

                std::pair<Position, Position> pair = OptimiseAnchor(&m_host_structure, guest, m_initial_anchor, Position{ x * max, y * max, z * max });

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
                }
                if (accept == false) {
                    excluded++;
                    continue;
                }
                m_anchor_accepted.push_back(pair.first);
                m_rotation_accepted.push_back(pair.second);

                Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, pair.first, pair.second);

                guest.setGeometry(destination);
                double distance = GeometryTools::Distance(pair.first, m_host_structure.Centroid());
                m_sum_distance += distance;
                m_docking_list.insert(std::pair<double, Vector>(distance, PositionPair2Vector(pair)));
                result = m_host_structure;
                for (std::size_t i = 0; i < guest.AtomCount(); ++i) {
                    molecule->addPair(guest.Atom(i));
                }
                molecule->setEnergy(distance);
                m_result_list.insert(std::pair<double, Molecule*>(distance, molecule));
            }
        }
        std::cout << (x / double(X)) * 100 << "% - " << m_anchor_accepted.size() << " stored structures. " << excluded << " structures were skipped, due to being duplicate!" << all << " checked!" << std::endl;
    }

    std::cout << m_anchor_accepted.size() << " stored structures. " << excluded << " structures were skipped, due to being duplicate!" << all << " checked!" << std::endl;

    guest = m_guest_structure;
    std::cout << std::endl
              << "** Docking Phase 0 - Finished **" << std::endl;
    int index = 0;
    std::vector<int> frags = { 0, 0 };
    for (const auto& pair : m_result_list) {
        ++index;
        frags[pair.second->GetFragments(1.3).size()]++;
        //std::cout << pair.first << std::endl;
        const std::string name = "Docking_B" + std::to_string(frags[pair.second->GetFragments(1.3).size()] / 100 + 1) + "_F" + std::to_string(pair.second->GetFragments(1.3).size()) + ".xyz";
        pair.second->appendXYZFile(name);
        if (!std::binary_search(m_files.begin(), m_files.end(), name))
            m_files.push_back(name);
    }
    PostOptimise();
}

void Docking::PostOptimise()
{
    int threads = 1; // xtb is not thread-safe yet
    std::map<double, Molecule*> result_list;
    auto iter = m_result_list.begin();
    while (iter != m_result_list.end()) {
        std::vector<Thread*> thread_block;
        std::cout << "Batch calculation started! " << std::endl;
        for (int i = 0; i < threads; ++i) {
            if (iter == m_result_list.end())
                continue;

            auto pair = *iter;

            Thread* th = new Thread;
            th->setMolecule(pair.second);
            th->start();
            thread_block.push_back(th);

            ++iter;
        }
        std::cout << "Batch evaluation ... " << std::endl;
        for (auto thread : thread_block) {
            thread->wait();

            Molecule* mol2 = new Molecule(thread->getMolecule());
            result_list.insert(std::pair<double, Molecule*>(mol2->Energy(), mol2));
            delete thread;
        }
        std::cout << "Done!" << std::endl;
    }
    std::cout << "** Docking Phase 2 - Finished **" << std::endl;
    /*
    for (const auto& pair : m_result_list) {
        Molecule *mol2 = new Molecule(OptimiseGeometry(pair.second, false, true, 1e4, 0.5));
        result_list.insert(std::pair<double, Molecule *>(mol2->Energy(), mol2));
        delete pair.second;
    }
    */
    std::vector<int> frags = { 0, 0 };
    int index = 0;
    for (const auto& pair : result_list) {
        ++index;
        frags[pair.second->GetFragments(1.3).size()]++;
        //std::cout << pair.first << std::endl;
        const std::string name = "Optimise_B" + std::to_string(frags[pair.second->GetFragments(1.3).size()] / 100 + 1) + "_F" + std::to_string(pair.second->GetFragments(1.3).size()) + ".xyz";
        pair.second->appendXYZFile(name);
        if (!std::binary_search(m_files.begin(), m_files.end(), name))
            m_files.push_back(name);

        delete pair.second;
    }
}
