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

#include "src/capabilities/confscan.h"

#include <algorithm>
#include <future>
#include <iostream>
#include <thread>

#include "json.hpp"
using json = nlohmann::json;

#include "docking.h"

Docking::Docking(const json& controller)
    : CurcumaMethod(DockingJson)
{
    UpdateController(controller);
    LoadControlJson();
}

void Docking::LoadControlJson()
{
    std::cout << m_controller << std::endl;
    double Pos_X = Json2KeyWord<double>(m_controller, "Pos_X");
    double Pos_Y = Json2KeyWord<double>(m_controller, "Pos_Y");
    double Pos_Z = Json2KeyWord<double>(m_controller, "Pos_Z");
    m_initial_anchor = Position{ Pos_X, Pos_Y, Pos_Z };
    m_step_X = Json2KeyWord<int>(m_controller, "Step_X");
    m_step_Y = Json2KeyWord<int>(m_controller, "Step_Y");
    m_step_Z = Json2KeyWord<int>(m_controller, "Step_Z");
    m_AutoPos = Json2KeyWord<bool>(m_controller, "AutoPos");
    m_PostFilter = Json2KeyWord<bool>(m_controller, "Filter");
    m_PostOptimise = Json2KeyWord<bool>(m_controller, "PostOpt");
    m_host = Json2KeyWord<std::string>(m_controller, "host");
    m_guest = Json2KeyWord<std::string>(m_controller, "guest");
    m_complex = Json2KeyWord<std::string>(m_controller, "complex");
}

bool Docking::Initialise()
{
    if ((m_host.compare("none") == 0 && m_guest.compare("none") == 0) && (m_complex.compare("none") == 0)) {
        AppendError("Neither host nor guest or a complex structure given. Sorry");
        return false;
    }
    Molecule host, guest, complex;
    bool host_loaded = true, guest_loaded = true, complex_loaded = true;
    try {
        host = Tools::LoadFile(m_host);
    } catch (int error) {
        host_loaded = false;
    }
    try {
        guest = Tools::LoadFile(m_guest);
    } catch (int error) {
        guest_loaded = false;
    }
    try {
        complex = Tools::LoadFile(m_complex);
    } catch (int error) {
        complex_loaded = false;
    }
    if (!host_loaded && !guest_loaded && !complex_loaded) {
        AppendError("Neither host nor guest or a complex structure where loaded. Sorry");
        return false;
    }

    if ((host.AtomCount() == 0 && guest.AtomCount() == 0) && (complex.AtomCount() == 0)) {
        AppendError("Host, guest or complex structure are empty. Nothing I can do now.");
        return false;
    }
    bool useComplex = false;
    if (complex.AtomCount()) {
        auto fragments = complex.GetFragments(m_scaling);
        if (fragments.size() == 2) {
            std::cout << "Complex structure is used." << std::endl;
            useComplex = true;

            m_host_structure.LoadMolecule(complex.getFragmentMolecule(0));
            m_guest_structure.LoadMolecule(complex.getFragmentMolecule(1));
            m_initial_anchor = m_guest_structure.Centroid(true);
            return true;
        }
    }
    if (!useComplex && host.AtomCount() && guest.AtomCount()) {
        m_host_structure.LoadMolecule(host);
        m_guest_structure.LoadMolecule(guest);
        return true;
    }
    AppendError("Hmm, something was missing. Please check your input structures.");
    return false;
}

void Docking::PerformDocking()
{
    double frag_scaling = 1.2;

    Molecule guest = m_guest_structure;
    Molecule result = m_host_structure;
    guest.CalculateMass();
    result.CalculateMass();
    m_fragments_mass.push_back(result.Mass());
    m_fragments_mass.push_back(guest.Mass());

    /* for( auto a : m_fragments_mass)
        std::cout << a << " " ;
    */
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

    int max_X = 360 / double(m_step_X);
    int max_Y = 360 / double(m_step_Y);
    int max_Z = 360 / double(m_step_Z);

    int excluded = 0, all = 0;
    for (int x = 0; x < m_step_X; ++x) {
        for (int y = 0; y < m_step_Y; ++y) {
            for (int z = 0; z < m_step_Z; ++z) {
                ++all;
                Molecule* molecule = new Molecule(m_host_structure);
                guest = m_guest_structure;

                std::pair<Position, Position> pair = OptimiseAnchor(&m_host_structure, guest, m_initial_anchor, Position{ x * max_X, y * max_Y, z * max_Z });

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
        std::cout << (x / double(m_step_X)) * 100 << "% - " << m_anchor_accepted.size() << " stored structures. " << excluded << " structures were skipped, due to being duplicate! " << all << " checked all together!" << std::endl;
    }

    std::cout << m_anchor_accepted.size() << " stored structures. " << excluded << " structures were skipped, due to being duplicate!" << all << " checked!" << std::endl;

    guest = m_guest_structure;
    std::cout << std::endl
              << "** Docking Phase 0 - Finished **" << std::endl;
    int index = 0;
    // Will be removed some day
    std::vector<int> frags = { 0, 0, 0, 0, 0 };
    for (const auto& pair : m_result_list) {
        ++index;
        frags[pair.second->GetFragments(frag_scaling).size()]++;
        //std::cout << pair.first << std::endl;
        const std::string name = "Docking_B" + std::to_string(frags[pair.second->GetFragments(frag_scaling).size()] / 100 + 1) + "_F" + std::to_string(pair.second->GetFragments(frag_scaling).size()) + ".xyz";
        pair.second->appendXYZFile(name);
        if (!std::binary_search(m_files.begin(), m_files.end(), name))
            m_files.push_back(name);
    }
    PostOptimise();
}

void Docking::PostOptimise()
{
    double frag_scaling = 1.5;
    int threads = 1; // xtb is not thread-safe yet
    std::map<double, Molecule*> result_list, final_results;
    auto iter = m_result_list.begin();
    int index = 0;
    json opt = OptJson;
    opt["dE"] = 50;
    opt["dRMSD"] = 0.1;
    while (iter != m_result_list.end()) {
        std::vector<Thread*> thread_block;
        std::cout << "Batch calculation started! " << std::endl;
        for (int i = 0; i < threads; ++i) {
            if (iter == m_result_list.end())
                continue;

            auto pair = *iter;

            Thread* th = new Thread;
            th->setMolecule(pair.second);
            th->setController(opt);
            th->start();
            thread_block.push_back(th);

            ++iter;
            index++;
        }
        std::cout << "Batch evaluation ... " << std::endl;
        for (auto thread : thread_block) {
            thread->wait();

            Molecule* mol2 = new Molecule(thread->getMolecule());
            result_list.insert(std::pair<double, Molecule*>(mol2->Energy(), mol2));
            delete thread;
        }
        std::cout << "Done! " << index << " of " << m_result_list.size() << " (" << index / double(m_result_list.size()) * 100 << " %)" << std::endl;
    }
    std::cout << "** Docking Phase 2 - Finished **" << std::endl;

    {
        // Will be removed some day
        std::vector<int> frags = { 0, 0, 0, 0, 0 };
        int index = 0;
        for (const auto& pair : result_list) {
            ++index;
            frags[pair.second->GetFragments(frag_scaling).size()]++;
            const std::string name = "Optimise_B" + std::to_string(frags[pair.second->GetFragments(frag_scaling).size()] / 100 + 1) + "_F" + std::to_string(pair.second->GetFragments(frag_scaling).size()) + ".xyz";
            pair.second->appendXYZFile(name);
            if (!std::binary_search(m_files.begin(), m_files.end(), name))
                m_files.push_back(name);

            if (pair.second->GetFragments(frag_scaling).size() == 2) {
                double sum = 0;
                std::vector<double> fragments = pair.second->FragmentMass();

                for (int i = 0; i < fragments.size(); ++i)
                    sum += abs(m_fragments_mass[i] - fragments[i]);
                for (auto a : fragments)
                    std::cout << a << " ";
                std::cout << " Difference " << sum;
                if (sum < 1e-3) {
                    final_results.insert(pair);
                    std::cout << " taking! E = " << pair.first << " Eh" << std::endl;
                } else {
                    std::cout << " skipping! E = " << pair.first << " Eh" << std::endl;
                }
            }
        }
    }

    std::cout << "** Docking Phase 3 - Fast Filtering structures with correct fragments **" << std::endl;

    json confscan = ConfScanJson;
    confscan["heavy"] = true;
    confscan["maxrank"] = -1;
    confscan["writeXYZ"] = false;
    confscan["ForceReorder"] = false;
    confscan["check"] = false;
    confscan["energy"] = 1.0;
    confscan["noname"] = true;
    confscan["preventReorder"] = true;
    confscan["maxenergy"] = -1;

    ConfScan* scan = new ConfScan(confscan);
    scan->setMolecules(final_results);
    scan->scan();

    std::vector<Molecule*> result = scan->Result();

    std::cout << "** Docking Phase 4 - Writing structures to Final_Result.xyz **" << std::endl;

    {
        int index = 0;
        for (const auto mol : result) {
            ++index;
            const std::string name = "Final_Result.xyz";
            mol->appendXYZFile(name);
            delete mol;
        }
    }
}
