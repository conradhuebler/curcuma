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

#include "src/tools/general.h"

#include "src/capabilities/confscan.h"
#include "src/capabilities/curcumaopt.h"
#include "src/core/parameter_registry.h"  // Claude Generated 2025
#include "src/capabilities/optimiser/LevMarDocking.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <algorithm>
#include <iostream>
#include <map>

#include "json.hpp"
using json = nlohmann::json;

#include "docking.h"

Docking::Docking(const json& controller, bool silent)
    : Docking(ConfigManager("docking", controller), silent)
{
}

Docking::Docking(const ConfigManager& config, bool silent)
    : CurcumaMethod(json{}, config.exportConfig(), silent)
{
    UpdateController(config.exportConfig());
}

void Docking::LoadControlJson()
{
    double Pos_X = Json2KeyWord<double>(m_defaults, "Pos_X");
    double Pos_Y = Json2KeyWord<double>(m_defaults, "Pos_Y");
    double Pos_Z = Json2KeyWord<double>(m_defaults, "Pos_Z");
    m_initial_anchor = { Position{ Pos_X, Pos_Y, Pos_Z } };
    m_step_X = Json2KeyWord<int>(m_defaults, "Step_X");
    m_step_Y = Json2KeyWord<int>(m_defaults, "Step_Y");
    m_step_Z = Json2KeyWord<int>(m_defaults, "Step_Z");
    m_AutoPos = Json2KeyWord<bool>(m_defaults, "AutoPos");
    m_PostFilter = Json2KeyWord<bool>(m_defaults, "Filter");
    m_PostOptimise = Json2KeyWord<bool>(m_defaults, "PostOpt");
    m_NoOpt = Json2KeyWord<bool>(m_defaults, "NoOpt");
    m_host = Json2KeyWord<std::string>(m_defaults, "host");
    m_guest = Json2KeyWord<std::string>(m_defaults, "guest");
    m_complex = Json2KeyWord<std::string>(m_defaults, "complex");
    m_centroid_max_distance = Json2KeyWord<double>(m_defaults, "CentroidMaxDistance");
    m_centroid_tol_distance = Json2KeyWord<double>(m_defaults, "CentroidTolDis");
    m_centroid_rot_distance = Json2KeyWord<double>(m_defaults, "RotationTolDis");
    m_energy_threshold = Json2KeyWord<double>(m_defaults, "EnergyThreshold");

    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_docking_threads = Json2KeyWord<int>(m_defaults, "DockingThreads");
    m_charge = Json2KeyWord<int>(m_defaults, "Charge");
    m_cycles = Json2KeyWord<int>(m_defaults, "Cycles");

    m_RMSDthreads = Json2KeyWord<int>(m_defaults, "RMSDThreads");
    m_RMSDElement = Json2KeyWord<int>(m_defaults, "RMSDElement");
    m_RMSDmethod = Json2KeyWord<std::string>(m_defaults, "RMSDMethod");
}

bool Docking::Initialise()
{
    if ((m_host.compare("none") == 0 && m_guest.compare("none") == 0) && (m_complex.compare("none") == 0)) {
        AppendError("Neither host nor guest or a complex structure given. Sorry");
        return false;
    }
    Molecule host, guest, complex;
    bool host_loaded = true, guest_loaded = true, complex_loaded = true;
    if (!(m_host.compare("none") == 0)) {
        try {
            host = Files::LoadFile(m_host);
        } catch (int error) {
            host_loaded = false;
        }
    } else
        host_loaded = false;
    if (!(m_guest.compare("none") == 0)) {
        try {
            guest = Files::LoadFile(m_guest);
        } catch (int error) {
            guest_loaded = false;
        }
    } else
        guest_loaded = false;

    if (!(m_complex.compare("none") == 0)) {
        try {
            complex = Files::LoadFile(m_complex);
        } catch (int error) {
            complex_loaded = false;
        }
    } else
        complex_loaded = false;

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
            useComplex = true;

            m_host_structure.LoadMolecule(complex.getFragmentMolecule(0));
            m_guest_structure.LoadMolecule(complex.getFragmentMolecule(1));
            m_initial_anchor = { m_guest_structure.Centroid(true) };
        }
    }
    if (!useComplex && host.AtomCount() && guest.AtomCount()) {
        m_host_structure.LoadMolecule(host);
        m_guest_structure.LoadMolecule(guest);
        return true;
    } else if (useComplex && !host_loaded && guest_loaded && guest.AtomCount()) {
        std::cout << "Complex structure is used, with substrat replace by -guest argument!" << std::endl;
        m_guest_structure.LoadMolecule(guest);
        return true;
    } else if (useComplex && !host_loaded && !guest_loaded) {
        std::cout << "Complex structure is used." << std::endl;
        return true;
    }
    AppendError("Hmm, something was missing. Please check your input structures.");
    return false;
}

void Docking::PerformDocking()
{
#ifndef USE_XTB
    m_PostOptimise = false;
#endif

    double frag_scaling = 1.2;

    Molecule guest = m_guest_structure;
    Molecule result = m_host_structure;
    guest.CalculateMass();
    result.CalculateMass();
    m_fragments_mass.push_back(result.Mass());
    m_fragments_mass.push_back(guest.Mass());

    std::cout << m_guest_structure.Centroid().transpose() << " = Centroid of Guest" << std::endl;

    if (m_AutoPos)
        m_initial_anchor = { m_host_structure.Centroid() };

    Geometry stored_guest = m_guest_structure.getGeometry();
    Position initial_centroid = m_guest_structure.Centroid();

    // Geometry geometry = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, m_initial_anchor, Position{ 0, 0, 0 });

    std::cout << std::endl
              << "** Docking Phase 0 - Starting **" << std::endl
              << std::endl;

    int max_X = 360 / double(m_step_X);
    int max_Y = 360 / double(m_step_Y);
    int max_Z = 360 / double(m_step_Z);

    int excluded = 0, all = 0, distance = 0;
    while (m_current_cycle < m_cycles) {
        std::cout << m_initial_anchor.size() << " initial anchors" << std::endl;
        if (m_NoOpt) // Loop unrolled
        {
            for (int x = 0; x < m_step_X; ++x) {
                for (int y = 0; y < m_step_Y; ++y) {
                    for (int z = 0; z < m_step_Z; ++z) {
                        ++all;
                        Molecule* molecule = new Molecule(m_host_structure);
                        guest = m_guest_structure;
                        for (const Position& anchor : m_initial_anchor) {
                            Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, anchor, Position{ x * max_X, y * max_Y, z * max_Z });
                            guest.setGeometry(destination);
                            for (std::size_t i = 0; i < guest.AtomCount(); ++i) {
                                molecule->addPair(guest.Atom(i));
                            }
                            m_docking_result.insert(std::pair<double, Molecule*>(all, molecule));
                        }
                    }
                }
                std::cout << (x / double(m_step_X)) * 100 << "% done." << std::endl;
            }
        } else {
            std::vector<DockThread*> threads;
            CxxThreadPool* pool = new CxxThreadPool;
            pool->setActiveThreadCount(m_threads);
            for (int x = 0; x < m_step_X; ++x) {
                for (int y = 0; y < m_step_Y; ++y) {
                    for (int z = 0; z < m_step_Z; ++z) {
                        for (const Position& anchor : m_initial_anchor) {
                            DockThread* thread = new DockThread(m_host_structure, guest);
                            thread->setPosition(anchor);
                            thread->setRotation(Position{ x * max_X, y * max_Y, z * max_Z });
                            pool->addThread(thread);
                            threads.push_back(thread);
                        }
                    }
                }
            }
            pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
            pool->DynamicPool();
            pool->StartAndWait();
            std::cout //<< std::endl
                << "** Docking Phase 0 - Finished - Now collection results **" << std::endl
                << std::endl;
            for (const auto* t : pool->getFinishedThreads()) {
                const DockThread* thread = static_cast<const DockThread*>(t);
                ++all;
                Molecule* molecule = new Molecule(m_host_structure);
                guest = m_guest_structure;

                bool accept = true;

                if (GeometryTools::Distance(thread->InitialPosition(), thread->LastPosition()) > m_centroid_max_distance) {
                    accept = false;
                    distance++;
                    m_initial_list.push_back(std::pair<Position, Position>(thread->InitialPosition(), thread->InitialRotation()));
                    continue;
                }

                for (std::size_t i = 0; i < m_anchor_accepted.size(); ++i) {
                    Position anchor = m_anchor_accepted[i];
                    Position rotation = m_rotation_accepted[i];
                    if (GeometryTools::Distance(anchor, thread->LastPosition()) < m_centroid_tol_distance || GeometryTools::Distance(rotation, thread->LastRotation()) < m_centroid_rot_distance)
                        accept = false;
                }
                if (accept == false) {
                    excluded++;
                    continue;
                }
                m_anchor_accepted.push_back(thread->LastPosition());
                m_rotation_accepted.push_back(thread->LastRotation());

                Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, thread->LastPosition(), thread->LastRotation());

                guest.setGeometry(destination);
                double distance = GeometryTools::Distance(thread->LastPosition(), m_host_structure.Centroid());
                m_sum_distance += distance;
                m_docking_list.insert(std::pair<double, Vector>(distance, PositionPair2Vector(std::pair<Position, Position>(thread->LastPosition(), thread->LastRotation()))));
                result = m_host_structure;
                for (std::size_t i = 0; i < guest.AtomCount(); ++i) {
                    molecule->addPair(guest.Atom(i));
                }
                molecule->setEnergy(distance);
                molecule->setCharge(m_charge);
                m_docking_result.insert(std::pair<double, Molecule*>(all, molecule));
                if (molecule->GetFragments(frag_scaling).size() == 2)
                    m_reuse_anchor.push_back(thread->LastPosition());
            }
            delete pool;
        }
        for (auto& init : m_initial_list) {
            Geometry destination = GeometryTools::TranslateAndRotate(stored_guest, initial_centroid, init.first, init.second);
            Molecule molecule = Molecule(m_host_structure);
            guest.setGeometry(destination);
            for (std::size_t i = 0; i < guest.AtomCount(); ++i) {
                molecule.addPair(guest.Atom(i));
            }
            molecule.appendXYZFile("Docking_Failed.xyz");
        }

        std::cout << m_anchor_accepted.size() << " stored structures. " << std::endl
                  << excluded << " structures were skipped, due to being duplicate!" << std::endl
                  << all << " checked!" << std::endl;
        guest = m_guest_structure;

        m_current_cycle++;
        m_initial_anchor = m_reuse_anchor;
    }
    std::cout << std::endl
              << "** Docking Phase 0 - Finished **" << std::endl;
    int index = 0;
    // Will be removed some day
    for (const auto& pair : m_docking_result) {
        ++index;
        std::string name;
        if (!m_NoOpt)
            name = "Docking_F" + std::to_string(pair.second->GetFragments(frag_scaling).size()) + ".xyz";
        else
            name = "Docking.xyz";
        pair.second->appendXYZFile(name);
        if (!std::binary_search(m_files.begin(), m_files.end(), name))
            m_files.push_back(name);
    }

    if (!m_NoOpt && m_PostOptimise) {
#ifdef USE_XTB
        PostOptimise();
#else
        std::cerr << "xtb support was not included into the binary. Sorry for that!" << std::endl;
#endif
    }
}

void Docking::PostOptimise()
{
    OptimiseBatch();

    CollectStructures();

    if (!m_PostFilter)
        return;

    std::cout << "** Docking Phase 3 - Filtering structures with hybrid reordering**" << std::endl;
    FilterStructures();
}

void Docking::OptimiseBatch()
{
    double frag_scaling = 1.5;
    /*
        json GFNFF = CurcumaOptJson;
        GFNFF["printOutput"] = false;
        GFNFF["threads"] = 1; // m_threads;
        GFNFF["gfn"] = 66;
        GFNFF["dE"] = 1;
        GFNFF["dRMSD"] = 0.1;
        GFNFF["GradNorm"] = 0.01;
        GFNFF["ConvCount"] = 3;
    */
    json GFN2_crude = CurcumaOptJson;
    GFN2_crude["SinglePoint"] = false;
    GFN2_crude["gfn"] = 2;
    GFN2_crude["threads"] = m_threads;
    GFN2_crude["dE"] = 10;
    GFN2_crude["dRMSD"] = 0.1;
    GFN2_crude["GradNorm"] = 0.01;
    GFN2_crude["ConvCount"] = 3;

    m_optimise = new CurcumaOpt(GFN2_crude, false);
    m_singlepoint = new CurcumaOpt(GFN2_crude, false);

    /* Divide list of structures in those which already have 2 fragments (supramolecular complex) and those which have atoms to close to each other*/
    auto iter = m_docking_result.begin();
    while (iter != m_docking_result.end()) {
        auto pair = *iter;
        if (pair.second->GetFragments(frag_scaling).size() == 2)
            m_optimise->addMolecule(pair.second);
        else if (pair.second->GetFragments(frag_scaling).size() == 1)
            m_singlepoint->addMolecule(pair.second);
        ++iter;
    }
    /* Optimisation */
    std::cout << "** Crude preoptimisation of " << m_optimise->Molecules()->size() << " complexes **" << std::endl;

    m_optimise->overrideBasename("Optimise_Crude_F2");
    m_optimise->start();

    /* Optimisation of the structures with incorrect fragments */
    std::cout << "** Crude preoptimisation of " << m_singlepoint->Molecules()->size() << " structures with wrong connectivity **" << std::endl;

    m_singlepoint->overrideBasename("Optimise_Crude_FX");
    m_singlepoint->start();

    /* Update Optimiser to perform single point calculation (GFN2) of the at GFN-FF level optimised structures */

    /*
        m_optimise->UpdateController(GFN2_crude);
        m_optimise->start();
    */
    /* Estimate the approximate energy of the correct complexes */
    double e0 = -1e7;
    for (const auto& m : *(m_optimise->Molecules())) {
        e0 = std::max(m.Energy(), e0);
    }
    std::cout << " *** energy threshold " << e0 << " Eh ***" << std::endl;

    // m_optimise->clear();

    const std::string exclude_energy = "AboveThreshold.xyz";
    int added = 0, excluded = 0;
    for (const auto& m : *(m_singlepoint->Molecules())) {
        if (abs(m.Energy() - e0) * 2625.5 < m_energy_threshold) {
            m_optimise->addMolecule(m);
            added++;
            /*
            Molecule* mol1 = new Molecule(m);
            m_temp_results.insert(std::pair<double, Molecule*>(m.Energy(), mol1));
            */
        } else {
            excluded++;
            m.appendXYZFile(exclude_energy);
        }
    }
    std::cout << " ***" << added << " complexes to optimisation batch ***" << std::endl;
    std::cout << " ***" << excluded << " structures dropped ***" << std::endl;

    json GFN2 = CurcumaOptJson;
    GFN2["gfn"] = 2;
    GFN2["threads"] = m_threads;
    GFN2["printOutput"] = false;

    std::cout << "** Standard optimsation of " << m_optimise->Molecules()->size() << " complexes **" << std::endl;

    m_optimise->UpdateController(GFN2);
    m_optimise->start();
}

void Docking::CollectStructures()
{
    double frag_scaling = 1.5;
    const std::string name = "Final_Result.xyz";
    const std::string excluded = "Excluded_Result.xyz";
    int added = 0, dropped = 0;

    for (const auto& m : *(m_optimise->Molecules())) {
        m.GetFragments(frag_scaling).size();
        if (m.GetFragments(frag_scaling).size() == 2) {
            double sum = 0;
            std::vector<double> fragments = m.FragmentMass();
            for (unsigned int i = 0; i < fragments.size(); ++i)
                sum += abs(m_fragments_mass[i] - fragments[i]);
            if (sum < 1e-3) {
                added++;
                m_optimisation_result.insert(std::pair<double, Molecule*>(m.Energy(), new Molecule(m)));
                // std::cout << " adding new Molecule! E = " << m.Energy() << " Eh\n\n"
                //           << std::endl;
                m.appendXYZFile(name);
            } else
                m.appendXYZFile(excluded);
            dropped++;
        }
    }
    std::cout << " ***" << added << " complexes to optimisation batch and written to" << name << " ***" << std::endl;
    std::cout << " ***" << dropped << " structures dropped  - written to " << excluded << " ***" << std::endl;
}

void Docking::FilterStructures()
{
    // Claude Generated 2025: Use ParameterRegistry instead of static JSON
    json confscan = ParameterRegistry::getInstance().getDefaultJson("confscan");
    confscan["forceReorder"] = true;
    confscan["silent"] = true;
    confscan["RMSDMethod"] = m_RMSDmethod;
    confscan["RMSDThreads"] = m_threads;
    confscan["RMSDElement"] = m_RMSDElement;
    // json controller;
    // controller["ConfScan"] = confscan;
    ConfScan* scan = new ConfScan(confscan, true);
    // scan->setMolecules(m_optimisation_result);
    scan->setFileName("Final_Result.xyz");
    scan->start();

    std::vector<Molecule*> result = scan->Result();
    std::cout << "** Docking Phase 4 - Structures were written to Final_Result.accepted.xyz **" << std::endl;
}
