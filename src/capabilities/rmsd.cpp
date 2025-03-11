/*
 * <RMSD calculator for chemical structures.>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "rmsd_functions.h"

extern "C" {
#include "src/capabilities/c_code/interface.h"
}

#include "munkres.h"

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include <cstdlib>
#include <stdio.h>

#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include <LBFGS.h>

#include "json.hpp"
using json = nlohmann::json;

#include "rmsd.h"
RMSDThread::RMSDThread(const Molecule& reference_molecule, const Molecule& target, const Geometry& reference, const Matrix& reference_topology, const std::vector<int> intermediate, double connected_mass, int element, int topo)
    : m_reference_molecule(reference_molecule)
    , m_target(target)
    , m_reference(reference)
    , m_reference_topology(reference_topology)
    , m_intermediate(intermediate)
    , m_connected_mass(connected_mass)
    , m_element(element)
    , m_topo(topo)
{
    if (m_topo == 0) {
        m_evaluator = [this](const Molecule& target_local) -> double {
            auto t = RMSDFunctions::getAligned(m_reference, target_local.getGeometry(), 1);
            return RMSDFunctions::getRMSD(m_reference, t);
        };
    } else if (m_topo == 1) {
        m_evaluator = [this](const Molecule& target_local) -> double {
            auto target_topo = target_local.DistanceMatrix().second;
            return (target_topo - m_reference_topology).cwiseAbs().sum();
        };
    } else {
        m_evaluator = [this](const Molecule& target_local) -> double {
            Molecule tar;
            Matrix reference_matrix = m_reference_molecule.DistanceMatrix().second;
            Geometry reference_geometry = m_reference_molecule.getGeometry();
            Geometry target_geometry = target_local.getGeometry();
            Geometry step = (reference_geometry - target_geometry) / m_topo;
            int topo = 0;
            for (int j = 1; j <= m_topo; ++j) {
                target_geometry += step;
                topo += CompareTopoMatrix(reference_matrix, tar.DistanceMatrix().second);
            }
            return topo;
        };
    }
}

int RMSDThread::execute()
{
    Molecule reference;
    Molecule target;

    for (int i = 0; i < m_intermediate.size(); i++) {
        target.addPair(m_target.Atom(m_intermediate[i]));
    }

    const int i = reference.AtomCount();
    Molecule reference_local(reference);
    std::map<double, int> match;
    bool found_none = true;

    for (int j = 0; j < m_target.AtomCount(); ++j) {
        if (m_target.Atoms()[j] == m_element) {
            found_none = false;
            Molecule target_local(target);
            if (target_local.addPair(m_target.Atom(j))) {
                /*
                //const auto t = RMSDFunctions::getAligned(m_reference, GeometryTools::TranslateGeometry(target_local.getGeometry(), target_local.Centroid(true), Position{ 0, 0, 0 }), 1);
                const auto t = RMSDFunctions::getAligned(m_reference, target_local.getGeometry(), 1);
                double value = 0;
                if(!m_topo)
                    value = RMSDFunctions::getRMSD(m_reference, t);
                else{
                    Matrix target_topo = target_local.DistanceMatrix().second;
                    value = (target_topo - m_reference_topology).cwiseAbs().sum();
                }
                */
                double value = m_evaluator(target_local);
                m_calculations++;
                if (target_local.AtomCount() <= m_target.AtomCount()) {
                    {
                        match.insert(std::pair<double, int>(value, j));
                    }
                }
            }
        }
    }
    m_match = match.size();
    for (const auto& element : match) {
        std::vector<int> temp = m_intermediate;
        temp.push_back(element.second);
        m_shelf.insert(std::pair<double, std::vector<int>>(element.first, temp));
    }

    return 0;
}

RMSDDriver::RMSDDriver(const json& controller, bool silent)
    : CurcumaMethod(RMSDJson, controller, silent)
{
    UpdateController(controller);
}

RMSDDriver::~RMSDDriver()
{
}

void RMSDDriver::LoadControlJson()
{
    m_fragment_reference = Json2KeyWord<int>(m_defaults, "fragment_reference");
    m_fragment_target = Json2KeyWord<int>(m_defaults, "fragment_target");
    m_fragment = Json2KeyWord<int>(m_defaults, "fragment");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_initial_fragment = Json2KeyWord<int>(m_defaults, "init");
    m_pt = Json2KeyWord<int>(m_defaults, "pt");
    m_molaligntol = Json2KeyWord<int>(m_defaults, "molaligntol");

    m_force_reorder = Json2KeyWord<bool>(m_defaults, "reorder");
    m_protons = !Json2KeyWord<bool>(m_defaults, "heavy");
    if (!m_verbose)
        m_silent = Json2KeyWord<bool>(m_defaults, "silent");
    m_intermedia_storage = Json2KeyWord<double>(m_defaults, "storage");
    m_dynamic_center = Json2KeyWord<bool>(m_defaults, "DynamicCenter");
    m_topo = Json2KeyWord<int>(m_defaults, "topo");
    m_write = Json2KeyWord<int>(m_defaults, "write");
    m_noreorder = Json2KeyWord<bool>(m_defaults, "noreorder");
    m_update_rotation = Json2KeyWord<bool>(m_defaults, "update-rotation");
    m_split = Json2KeyWord<bool>(m_defaults, "split");
    m_nofree = Json2KeyWord<bool>(m_defaults, "nofree");
    m_maxtrial = Json2KeyWord<int>(m_defaults, "maxtrial");
    m_kmstat = Json2KeyWord<bool>(m_defaults, "kmstat");
    m_km_convergence = Json2KeyWord<double>(m_defaults, "km_conv");

#pragma message("these hacks to overcome the json stuff are not nice, TODO!")
    try {
        std::string element = m_defaults["Element"].get<std::string>();
        StringList elements = Tools::SplitString(element, ",");
        for (const std::string& str : elements) {
            try {
                m_element_templates.push_back(std::stod(str));
            } catch (const std::invalid_argument& arg) {
                continue;
            }
        }
        if (m_element_templates.size())
            m_element = m_element_templates[0];
    } catch (const nlohmann::detail::type_error& error) {
        m_element = Json2KeyWord<int>(m_defaults, "element");
        m_element_templates.push_back(m_element);
    }

    m_check = Json2KeyWord<bool>(m_defaults, "check");
    m_damping = Json2KeyWord<double>(m_defaults, "damping");
    m_molalign = Json2KeyWord<std::string>(m_defaults, "molalignbin");

    std::string method = Json2KeyWord<std::string>(m_defaults, "method");
    m_munkress_cycle = 1;
    m_limit = 0;
    if (method.compare("template") == 0)
        m_method = 2;
    else if (method.compare("incr") == 0)
        m_method = 1;
    else if (method.compare("hybrid0") == 0)
        m_method = 3;
    else if (method.compare("hybrid") == 0 || method.compare("subspace") == 0) {
        m_method = 4;
        m_limit = Json2KeyWord<int>(m_defaults, "limit");
    } else if (method.compare("free") == 0) {
        m_method = 5;
    } else if (method.compare("molalign") == 0)
        m_method = 6;
    else if (method.compare("dtemplate") == 0) {
        m_method = 7;
        m_limit = Json2KeyWord<int>(m_defaults, "limit");
    } else
        m_method = 1;
    if (!m_silent) {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nPermutation of atomic indices performed according to {0} \n\n", m_method);
    }
    m_costmatrix = Json2KeyWord<int>(m_defaults, "costmatrix");
    std::string order = Json2KeyWord<std::string>(m_defaults, "order");
    int cycles = Json2KeyWord<int>(m_defaults, "cycles");
    if (cycles != -1)
        m_munkress_cycle = cycles;

    if (!order.empty() && std::filesystem::exists(order)) {
        std::ifstream file(order);
        if (file.is_open()) {
            std::string line;
            std::string fileContent;
            while (std::getline(file, line)) {
                fileContent += line;
            }
            file.close();

            // Parse den Inhalt der Datei
            // m_reorder_rules = Tools::ParseStringToVector(fileContent);
            std::vector<int> vector = Tools::String2Vector(fileContent);
            std::cout << vector.size() << std::endl;
            if (vector.size() != 0) {
                m_reorder_rules = vector;
                m_method = 10;
            }
        }
    }

    std::cout << m_defaults["reference_atoms"] << std::endl;
    std::string reference_atoms = Json2KeyWord<std::string>(m_defaults, "reference_atoms");
    std::string target_atoms = Json2KeyWord<std::string>(m_defaults, "target_atoms");
    if (reference_atoms.size() > 0)
        m_reference_atoms = Tools::ParseStringToVector(reference_atoms);

    if (target_atoms.size() > 0) {
        m_target_atoms = Tools::ParseStringToVector(target_atoms);
    }
    if (m_reference_atoms.size() != m_target_atoms.size()) {
        std::cout << "Number of reference and target atoms must be the same, exiting ...";
        exit(1);
    } else if (m_reference_atoms.size() > 0) {
        std::cout << "Aligning molecules using atoms: ";
        for (int i = 0; i < m_reference_atoms.size(); ++i)
            std::cout << m_reference_atoms[i] << " " << m_target_atoms[i] << "\n";
    }
    if (!m_silent) {
        fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold, "\nCurrent Configuration:\n");
        fmt::print("Fragment Reference: {}\n", m_fragment_reference);
        fmt::print("Fragment Target: {}\n", m_fragment_target);
        fmt::print("Fragment: {}\n", m_fragment);
        fmt::print("Threads: {}\n", m_threads);
        fmt::print("Initial Fragment: {}\n", m_initial_fragment);
        fmt::print("PT: {}\n", m_pt);
        fmt::print("Molalign Tolerance: {}\n", m_molaligntol);
        fmt::print("Force Reorder: {}\n", m_force_reorder);
        fmt::print("Protons: {}\n", m_protons);
        fmt::print("Silent: {}\n", m_silent);
        fmt::print("Intermedia Storage: {}\n", m_intermedia_storage);
        fmt::print("Dynamic Center: {}\n", m_dynamic_center);
        fmt::print("Topo: {}\n", m_topo);
        fmt::print("Write: {}\n", m_write);
        fmt::print("No Reorder: {}\n", m_noreorder);
        fmt::print("Update Rotation: {}\n", m_update_rotation);
        fmt::print("Split: {}\n", m_split);
        fmt::print("No Free: {}\n", m_nofree);
        fmt::print("Max Trial: {}\n", m_maxtrial);
        fmt::print("KM Stat: {}\n", m_kmstat);
        fmt::print("KM Convergence: {}\n", m_km_convergence);
        fmt::print("Element: {}\n", m_element);
        fmt::print("Check: {}\n", m_check);
        fmt::print("Damping: {}\n", m_damping);
        fmt::print("Molalign Bin: {}\n", m_molalign);
        fmt::print("Method: {}\n", m_method);
        fmt::print("Cost Matrix: {}\n", m_costmatrix);
        fmt::print("Order: {}\n", order);
        fmt::print("Cycles: {}\n", m_munkress_cycle);
        fmt::print("Reference Atoms: {}\n", fmt::join(m_reference_atoms, ", "));
        fmt::print("Target Atoms: {}\n", fmt::join(m_target_atoms, ", "));
    }
}

void RMSDDriver::setMatchingAtoms(const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms)
{
    m_reference_atoms = reference_atoms;
    m_target_atoms = target_atoms;
#pragma message("TODO: remove notreorder in the future")
    m_noreorder = true;
}

void RMSDDriver::start()
{
    if (m_reference.AtomCount() == 0 || m_target.AtomCount() == 0) {
        std::cout << "At least one structure is empty, exiting ...";
        return;
    }
    RunTimer timer(false);
    clear();
    bool rmsd_calculated = false;

    if (m_reference_atoms.size() == m_target_atoms.size() && m_reference_atoms.size() > 0) {
        m_target_original = m_target;
        m_reference_original = m_reference;
        m_target.clear();
        m_reference.clear();
        for (int i = 0; i < m_reference_atoms.size(); ++i) {
            m_reference.addPair(m_reference_original.Atom(m_reference_atoms[i]));
            m_target.addPair(m_target_original.Atom(m_target_atoms[i]));
        }
    }

    if (m_reference.AtomCount() < m_target.AtomCount()) {
        m_swap = true;
        Molecule tmp = m_reference;
        m_reference = m_target;
        m_target = tmp;
    }

    if (m_initial_fragment != -1 && m_initial.size() == 0)
        m_initial = m_reference.GetFragments()[m_initial_fragment];

    if (m_initial.size())
        InitialiseOrder();

    if (m_fragment != -1) {
        m_fragment_reference = m_fragment;
        m_fragment_target = m_fragment;
    }

    if (m_protons == false)
        ProtonDepleted();

    int reference_fragments = m_reference.GetFragments(m_scaling).size();
    int target_fragments = m_target.GetFragments(m_scaling).size();

    m_reference.InitialiseConnectedMass(1.5, m_protons);
    m_target.InitialiseConnectedMass(1.5, m_protons);

    m_reference_centered = m_reference;
    m_target_centered = m_target;

    m_reference.Center();
    m_target.Center();

    m_target_aligned = m_target;
    m_reference_aligned.LoadMolecule(m_reference);
    if (m_reference.Atoms() != m_target.Atoms() || m_force_reorder) {
        if (!m_noreorder || m_method == 10) {
            // std::cout << "Reordering Atoms" << std::endl;
            ReorderMolecule();
            rmsd_calculated = true;
        }
    }
    Molecule temp_ref, temp_tar;
    if (m_target_atoms.size() == 0) {
        int consent = true;
        for (int i = 0; i < m_reference.AtomCount() && i < m_target.AtomCount(); ++i) {
            if (m_reference.Atom(i).first == m_target.Atom(i).first) {
                temp_ref.addPair(m_reference.Atom(i));
                temp_tar.addPair(m_target.Atom(i));
            } else
                consent = false;
        }
        if (consent) {
            if (m_fragment_reference != -1 && m_fragment_target != -1) {
                m_rmsd = CustomRotation();
            } else if (!rmsd_calculated)
                m_rmsd = BestFitRMSD();
        } else {
            if (!m_silent)
                fmt::print("Partial RMSD is calculated, only from those atoms, that match each other.\n\n\n");
            m_rmsd = PartialRMSD(temp_ref, temp_tar);
        }
    } else {
        m_rmsd = BestFitRMSD();
        Eigen::Matrix3d R = m_rotation;
        auto reference = CenterMolecule(m_reference_original.getGeometry());
        auto target = CenterMolecule(m_target_original.getGeometry());
        m_target.clear();
        m_reference.clear();
        m_reference_aligned = m_reference_original;
        m_target_aligned = m_target_original;
        const auto t = RMSDFunctions::applyRotation(target, R);
        m_reference_aligned.setGeometry(reference);
        m_target_aligned.setGeometry(t);
        m_reference.setGeometry(reference);
        m_target.setGeometry(t);
        //        rmsd = RMSDFunctions::getRMSD(reference, t);
    }

    m_htopo_diff = CompareTopoMatrix(m_reference_aligned.HydrogenBondMatrix(-1, -1), m_target_aligned.HydrogenBondMatrix(-1, -1));
    if (!m_silent) {
        std::cout << std::endl
                  << "RMSD calculation took " << timer.Elapsed() << " msecs." << std::endl;
        std::cout << "Difference in Topological Hydrogen Bond Matrix is " << m_htopo_diff << std::endl;
        std::cout << std::endl
                  << std::endl;
        CheckTopology();
    }
    if (m_swap) {
        Molecule reference = m_reference;
        Molecule reference_aligned = m_reference_aligned;
        m_reference = m_target;
        m_reference_aligned = m_target_aligned;
        m_target = reference;
        m_target_aligned = reference_aligned;
    }
}

double RMSDDriver::SimpleRMSD()
{
    double rmsd = 0;
    rmsd = RMSDFunctions::getRMSD(m_reference.getGeometry(), m_target.getGeometry());
    return rmsd;
}

double RMSDDriver::BestFitRMSD()
{
    double rmsd = 0;
    auto reference = CenterMolecule(m_reference.getGeometry());
    auto target = CenterMolecule(m_target.getGeometry());
    // const auto t = RMSDFunctions::getAligned(reference, target, 1);
    Eigen::Matrix3d R = RMSDFunctions::BestFitRotation(reference, target, 1);
    const auto t = RMSDFunctions::applyRotation(target, R);
    m_reference_aligned.setGeometry(reference);
    m_target_aligned.setGeometry(t);
    m_reference.setGeometry(reference);
    m_target.setGeometry(t);
    rmsd = RMSDFunctions::getRMSD(reference, t);
    m_rmsd = rmsd;
    m_rotation = R;
    return rmsd;
}

double RMSDDriver::PartialRMSD(const Molecule& ref, const Molecule& tar)
{
    double rmsd = 0;
    {
        auto operators = GetOperateVectors(ref, tar);
        Eigen::Matrix3d R = operators.first;
        auto reference = GeometryTools::TranslateGeometry(m_reference.getGeometry(), ref.Centroid(), Position{ 0, 0, 0 });
        auto target = GeometryTools::TranslateGeometry(m_target.getGeometry(), tar.Centroid(), Position{ 0, 0, 0 });
        m_reference_aligned.setGeometry(reference);
        m_target_aligned.setGeometry(RMSDFunctions::applyRotation(target, R));
    }
    {
        auto reference = CenterMolecule(ref.getGeometry());
        auto target = CenterMolecule(tar.getGeometry());
        const auto t = RMSDFunctions::getAligned(reference, target, 1);
        rmsd = RMSDFunctions::getRMSD(reference, t);
    }
    return rmsd;
}

double RMSDDriver::CustomRotation()
{
    double rmsd = 0;
    int fragment_reference = m_fragment_reference;
    int fragment_target = m_fragment_target;
    auto reference_frag = CenterMolecule(m_reference.getGeometryByFragment(fragment_reference));
    auto target_frag = CenterMolecule(m_target.getGeometryByFragment(fragment_target));

    Eigen::Matrix3d rotation = RMSDFunctions::BestFitRotation(reference_frag, target_frag, 1);

    auto reference = GeometryTools::TranslateGeometry(m_reference.getGeometry(), GeometryTools::Centroid(m_reference.getGeometryByFragment(fragment_reference)), Position{ 0, 0, 0 }); // CenterMolecule(reference_mol);
    auto target = GeometryTools::TranslateGeometry(m_target.getGeometry(), GeometryTools::Centroid(m_target.getGeometryByFragment(fragment_target)), Position{ 0, 0, 0 }); // CenterMolecule(target_mol);
    const auto t = RMSDFunctions::applyRotation(target, rotation);
    m_reference_aligned.setGeometry(reference);
    m_target_aligned.setGeometry(t);
    rmsd = RMSDFunctions::getRMSD(reference, t);
    m_rotation = rotation;
    return rmsd;
}

void RMSDDriver::ReorderIncremental()
{
    int inter_size = m_reference.AtomCount() * (m_reference.AtomCount() - 1) * m_intermedia_storage;

    m_reorder_reference = m_reference;
    m_reorder_target = m_target;

    m_reorder_reference.setGeometry(CenterMolecule(m_reference.getGeometry()));
    m_reorder_target.setGeometry(CenterMolecule(m_target.getGeometry()));

    std::pair<Molecule, LimitedStorage> result = InitialisePair();
    Molecule ref = result.first;
    LimitedStorage storage_shelf = result.second;

    int reference_reordered = 0;
    int reference_not_reorordered = 0;
    int max = std::min(m_reference.AtomCount(), m_target.AtomCount());
    int combinations = 0;
    int wake_up = 100;
    CxxThreadPool* pool = new CxxThreadPool;
    if (m_silent)
        pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    else
        pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
    pool->setActiveThreadCount(m_threads);
    std::vector<AtomDef> atoms;
    while (
        m_reorder_reference_geometry.rows() < m_reorder_reference.AtomCount() && m_reorder_reference_geometry.rows() < m_reorder_target.AtomCount() && ((reference_reordered + reference_not_reorordered) <= m_reference.AtomCount())) {
        int thread_count = 0;

        Molecule reference = ref;
        int i = reference.AtomCount();
        double mass = m_reference.ConnectedMass(i);
        auto atom = m_reorder_reference.Atom(i);
        if (!m_silent) {
            std::cout << int((reference_reordered + reference_not_reorordered) / double(max) * 100) << " % " << std::endl;
        }
        int element = atom.first;
        reference.addPair(atom);
        if (m_dynamic_center)
            m_reorder_reference_geometry = GeometryTools::TranslateGeometry(reference.getGeometry(), reference.Centroid(true), Position{ 0, 0, 0 });
        else
            m_reorder_reference_geometry = reference.getGeometry();
        std::vector<RMSDThread*> threads;
        for (const auto& e : *storage_shelf.data()) {
            RMSDThread* thread = new RMSDThread(reference, m_reorder_target, m_reorder_reference_geometry, reference.DistanceMatrix().second, e.second, mass, element, m_topo);
            pool->addThread(thread);
            threads.push_back(thread);
            thread_count++;
        }
        pool->StaticPool();
        pool->setWakeUp(wake_up);
        int match = 0;
        /* For now, lets just dont start the threads if the current element can not be found in target */
        if (std::find(m_target.Atoms().begin(), m_target.Atoms().end(), element) != m_target.Atoms().end()) {
            pool->StartAndWait();
        }
        LimitedStorage storage_shelf_next(inter_size);
        for (const auto t : threads) {
            RMSDThread* thread = static_cast<RMSDThread*>(t);
            for (const auto& item : (*thread->data())) {
                if (thread->Match()) {
                    storage_shelf_next.addItem(item.first, item.second);
                }
                match += thread->Match();
            }
            combinations += thread->Calculations();
        }
        if (match == 0) {
            Molecule ref_0;
            auto atom = m_reorder_reference.Atom(i);
            for (std::size_t index = 0; index < m_reference.AtomCount(); index++) {
                if (i != index)
                    ref_0.addPair(m_reference.Atom(index));
            }
            if (std::find(atoms.begin(), atoms.end(), m_reorder_reference.Atom(i)) == atoms.end()) {
                ref_0.addPair(m_reorder_reference.Atom(i));
                atoms.push_back(m_reorder_reference.Atom(i));
                m_reorder_reference = ref_0;
                fmt::print("Atom order of reference molecule was altered!\nThe current atom (formerly {}) was pushed to the end of the list!\n", i + reference_reordered);
                fmt::print("Element {} : ({:f} {:f} {:f})\n", atom.first, atom.second[0], atom.second[1], atom.second[2]);
                reference_reordered++;
            } else
                reference_reordered = m_reference.AtomCount();
        } else {
            ref = reference;
            storage_shelf = storage_shelf_next;
            int index = 0;
            for (const auto& i : *storage_shelf.data()) {
                m_intermedia_rules.push_back(i.second);
                if (index > m_limit)
                    break;
                index++;
            }
            reference_not_reorordered++;
        }
        wake_up = 2 * pool->WakeUp();
        pool->clear();
    }
    delete pool;

    int count = 0;
    for (const auto& e : *storage_shelf.data()) {
        if (count > max * (max - 1) * m_intermedia_storage)
            continue;
        std::vector<int> rule = FillMissing(m_reference, e.second);
        if (std::find(m_stored_rules.begin(), m_stored_rules.end(), rule) == m_stored_rules.end()) {
            m_stored_rules.push_back(rule);
            m_intermedia_rules.push_back(rule);
            m_intermediate_cost_matrices.insert(std::pair<double, std::vector<int>>(e.first, rule));
            count++;
        }
    }
    // m_reorder_rules = m_results.begin()->second;
    if (m_stored_rules.size() == 0) {
        if (!m_silent)
            std::cout << "No new solution found, sorry" << std::endl;
        for (int i = 0; i < m_reference.AtomCount(); ++i)
            m_reorder_rules.push_back(i);
    } else {
        m_reorder_rules = m_stored_rules[0];
    }
    m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
    m_target = m_target_reordered;
    m_target_aligned = m_target;
    if (!m_silent)
        std::cout << "Overall " << combinations << " where evaluated!" << std::endl;
}

std::vector<int> RMSDDriver::FillMissing(const Molecule& molecule, const std::vector<int>& order)
{
    std::vector<int> result = order;
    for (int i = 0; i < molecule.AtomCount(); ++i) {
        if (std::find(result.begin(), result.end(), i) == result.end())
            result.push_back(i);
    }
    return result;
}

double RMSDDriver::Rules2RMSD(const std::vector<int> rules, int fragment)
{
    Molecule target;
    for (auto i : rules)
        target.addPair(m_target.Atom(i));
    int tmp_ref = m_fragment_reference;
    int tmp_tar = m_fragment_target;

    m_fragment_target = fragment;
    m_fragment_reference = fragment;
    Molecule ref;
    Molecule tar;
    double rmsd = CalculateRMSD(m_reference, target, &ref, &tar);
    m_htopo_diff = CompareTopoMatrix(ref.HydrogenBondMatrix(-1, -1), tar.HydrogenBondMatrix(-1, -1));

    m_fragment_reference = tmp_ref;
    m_fragment_target = tmp_tar;
    return rmsd;
}

StructComp RMSDDriver::Rule2RMSD(const std::vector<int> rules, int fragment)
{
    StructComp result;
    Molecule target;
    for (auto i : rules)
        target.addPair(m_target.Atom(i));
    int tmp_ref = m_fragment_reference;
    int tmp_tar = m_fragment_target;

    m_fragment_target = fragment;
    m_fragment_reference = fragment;
    Molecule ref;
    Molecule tar;
    Matrix topo_reference = m_reference.DistanceMatrix().second;
    Matrix topo_target = target.DistanceMatrix().second;
    result.rmsd = CalculateRMSD(m_reference, target, &ref, &tar);
    result.diff_hydrogen_bonds = CompareTopoMatrix(ref.HydrogenBondMatrix(-1, -1), tar.HydrogenBondMatrix(-1, -1));
    result.diff_topology = (topo_target - topo_reference).cwiseAbs().sum();
    m_fragment_reference = tmp_ref;
    m_fragment_target = tmp_tar;
    return result;
}

void RMSDDriver::clear()
{
    m_results.clear();
    m_stored_rules.clear();
    m_prepared_cost_matrices.clear();
    m_intermediate_cost_matrices.clear();
    m_intermedia_rules.clear();

    m_rmsd = 0.0;
}

void RMSDDriver::CheckTopology()
{
    Matrix reference_matrix = m_reference.DistanceMatrix().second;
    int index = 0;
    int best_topo = reference_matrix.size() * 2;
    int best_index = 0;
    std::vector<int> best_rule;
    for (const auto& rule : m_stored_rules) {
        if (index >= m_write)
            break;

        index++;
        if (index == 1)
            continue;

        Molecule target;

        for (auto i : rule)
            target.addPair(m_target_original.Atom(i));
        Molecule ref;
        Molecule tar;
        double rmsd = CalculateRMSD(m_reference, target, &ref, &tar);
        tar.writeXYZFile("target." + std::to_string(index) + ".reordered.xyz");
        int topo0 = CompareTopoMatrix(reference_matrix, tar.DistanceMatrix().second);

        Geometry reference_geometry = ref.getGeometry();
        Geometry target_geometry = tar.getGeometry();
        Geometry step = (reference_geometry - target_geometry) * 0.1;
        int topo = 0;
        for (int j = 1; j <= 10; ++j) {
            target_geometry += step;
            tar.setGeometry(target_geometry);
            // tar.appendXYZFile(std::to_string(index) + ".test.xyz");
            topo += CompareTopoMatrix(reference_matrix, tar.DistanceMatrix().second);
        }
        // ref.appendXYZFile(std::to_string(index) + ".test.xyz");

        std::cout << index << " " << rmsd << " " << topo0 << " " << topo << std::endl;
        /*  if (topo < best_topo) {
              best_topo = topo;
              best_rule = rule;
              best_index = index;
          }*/
    }
    /*
    if(index > 1)
        std::cout << best_index << " " << best_topo << std::endl;
        */
}

void RMSDDriver::ProtonDepleted()
{
    if (!m_silent)
        std::cerr << "Will perform calculation on proton depleted structure." << std::endl;

    Molecule reference;
    for (std::size_t i = 0; i < m_reference.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_reference.Atom(i);
        if (atom.first != 1)
            reference.addPair(atom);
    }

    Molecule target;
    for (std::size_t i = 0; i < m_target.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_target.Atom(i);
        if (atom.first != 1)
            target.addPair(atom);
    }
    m_reference = reference;
    m_target = target;
    m_init_count = m_heavy_init;
}

double RMSDDriver::CalculateRMSD(const Molecule& reference_mol, const Molecule& target_mol, Molecule* ret_ref, Molecule* ret_tar, int factor) const
{
    double rmsd = 0;
    auto reference = CenterMolecule(reference_mol.getGeometry());
    auto target = CenterMolecule(target_mol.getGeometry());
    const auto t = RMSDFunctions::getAligned(reference, target, 1);
    if (ret_ref != NULL) {
        ret_ref->LoadMolecule(reference_mol);
        ret_ref->setGeometry(reference);
    }
    if (ret_tar != NULL) {
        ret_tar->LoadMolecule(target_mol);
        ret_tar->setGeometry(t);
    }
    rmsd = RMSDFunctions::getRMSD(reference, t);
    return rmsd;
}

std::vector<double> RMSDDriver::IndivRMSD(const Molecule& reference_mol, const Molecule& target_mol, int factor) const
{
    Eigen::Matrix3d R = RMSDFunctions::BestFitRotation(reference_mol, target_mol, factor);

    std::vector<double> terms;

    Geometry reference = CenterMolecule(reference_mol, m_fragment_reference);
    Geometry target = CenterMolecule(target_mol, m_fragment_target);
    Geometry rotated = target * R;

    for (int i = 0; i < rotated.rows(); ++i) {
        terms.push_back((rotated(i, 0) - reference(i, 0)) * (rotated(i, 0) - reference(i, 0)) + (rotated(i, 1) - reference(i, 1)) * (rotated(i, 1) - reference(i, 1)) + (rotated(i, 2) - reference(i, 2)) * (rotated(i, 2) - reference(i, 2)));
    }
    return terms;
}

double RMSDDriver::CalculateRMSD()
{
    Molecule *reference = new Molecule, *target = new Molecule;
    double rmsd = CalculateRMSD(m_reference, m_target, reference, target);

    m_reference_aligned.LoadMolecule(reference);
    m_target_aligned.LoadMolecule(target);

    return rmsd;
}

Geometry RMSDDriver::CenterMolecule(const Molecule& mol, int fragment) const
{
    const Geometry cached = mol.getGeometryByFragment(fragment, m_protons);
    return GeometryTools::TranslateGeometry(cached, GeometryTools::Centroid(cached), Position{ 0, 0, 0 });
}

Geometry RMSDDriver::CenterMolecule(const Geometry& mol) const
{
    return GeometryTools::TranslateGeometry(mol, GeometryTools::Centroid(mol), Position{ 0, 0, 0 });
}

void RMSDDriver::InitialiseOrder()
{
    Molecule target, reference;

    for (int a : m_initial) {
        m_heavy_init += m_reference.Atom(a).first != 1;
        reference.addPair(m_reference.Atom(a));
        target.addPair(m_target.Atom(a));
    }

    /* As they may be different, we have two loops working here */
    for (int i = 0; i < m_reference.AtomCount(); ++i)
        if (!reference.Contains(m_reference.Atom(i)))
            reference.addPair(m_reference.Atom(i));

    for (int i = 0; i < m_target.AtomCount(); ++i)
        if (!target.Contains(m_target.Atom(i)))
            target.addPair(m_target.Atom(i));

    m_reference = reference;
    m_target = target;
    m_init_count = m_initial.size();
}

std::pair<Molecule, LimitedStorage> RMSDDriver::InitialisePair()
{
    Molecule reference;
    LimitedStorage storage(m_reference.AtomCount() * (m_reference.AtomCount() - 1) * m_intermedia_storage);
    int index = 0;
    std::vector<int> elements = m_reorder_reference.Atoms();
    if (m_initial.size() == 0) {
        for (int i = 0; i < m_reorder_reference.AtomCount() && index < 2; i++) {
            reference.addPair(m_reorder_reference.Atom(i));
            index++;
        }
        m_reorder_reference_geometry = reference.getGeometry();

        std::vector<int> elements_target = m_reorder_target.Atoms();
        std::vector<int> tmp_reference = reference.Atoms();

        for (int i = 0; i < m_reorder_target.AtomCount(); ++i) {
            for (int j = i + 1; j < m_reorder_target.AtomCount(); ++j) {
                if (tmp_reference[0] == elements_target[i] && tmp_reference[1] == elements_target[j])
                    storage.addItem(storage.size(), { i, j });
                if (tmp_reference[0] == elements_target[j] && tmp_reference[1] == elements_target[i])
                    storage.addItem(storage.size(), { j, i });
            }
        }
    } else {
        std::vector<int> start;
        for (int i = 0; i < m_init_count; ++i)
            start.push_back(i);
        m_intermediate_results.push_back(start);
    }
    return std::pair<Molecule, LimitedStorage>(reference, storage);
}

void RMSDDriver::ReorderMolecule()
{
    auto R = GetOperateVectors(m_reference, m_target);
    Eigen::Matrix3d rotation = R.first;

    auto blob = MakeCostMatrix(OptimiseRotation(rotation));

    // if (m_method == 4)
    //     blob.first = 1;
    if (m_method != 6)
        InsertRotation(blob);
    if (m_method == 1)
        ReorderIncremental();
    else if (m_method == 2)
        TemplateReorder();
    else if (m_method == 3)
        HeavyTemplate();
    else if (m_method == 4)
        AtomTemplate();
    else if (m_method == 5)
        TemplateFree();
    else if (m_method == 6) {
        if (!MolAlignLib()) {
            InsertRotation(blob);
            FinaliseTemplate();

            m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
            m_target_aligned = m_target;
            m_reorder_rules = m_results.begin()->second;
            m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
            m_target = m_target_reordered;
            m_rmsd = m_results.begin()->first;
            return;
        }
    } else if (m_method == 7)
        DistanceTemplate();
    if (m_method != 6 && m_method != 10) {
        FinaliseTemplate();

        m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
        m_target_aligned = m_target;
        m_reorder_rules = m_results.begin()->second;
        m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
        m_target = m_target_reordered;
        m_rmsd = m_results.begin()->first;
    } else if (m_method == 10) {
        m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
        m_rmsd = Rules2RMSD(m_reorder_rules);
        std::cout << m_rmsd << std::endl;
    }
}

void RMSDDriver::FinaliseTemplate()
{
    Molecule target = m_target;
    std::vector<std::vector<int>> rules = m_stored_rules;
    double rmsd_prev = 10;
    int eq_counter = 0, incr_counter = 0;
    int iter = 0;
    RunTimer time;
    double prev_cost = -1;
    std::ofstream result_file;
    if (m_kmstat)
        result_file.open("kmstat.dat", std::ios_base::app);
    for (auto permutation : m_prepared_cost_matrices) {
        iter++;

        if (eq_counter > m_maxtrial || incr_counter > m_maxtrial || iter > m_maxtrial)
            break;
        if (prev_cost != -1 && prev_cost < permutation.first / 10) {
            break;
        }

        prev_cost = permutation.first;
        auto result = SolveCostMatrix(permutation.second);

        m_target_reordered = ApplyOrder(result, target);
        double rmsd = Rules2RMSD(result);
        if (!m_silent)
            std::cout << permutation.first << " " << rmsd << " " << eq_counter << " " << time.Elapsed() << std::endl;
        if (m_kmstat)
            result_file << permutation.first << " " << rmsd << " " << std::endl;

        m_results.insert(std::pair<double, std::vector<int>>(rmsd, result));
        time.Reset();

        eq_counter += std::abs(rmsd - rmsd_prev) < m_km_convergence;
        incr_counter += rmsd > rmsd_prev;
        rmsd_prev = rmsd;
    }
    if (m_kmstat)
        result_file.close();

    m_reorder_rules = m_results.begin()->second;
}

void RMSDDriver::DistanceTemplate()
{
    if (!m_silent)
        std::cout << "Prepare template structure on atom distances:" << std::endl;

    auto pairs = PrepareDistanceTemplate();
}

void RMSDDriver::AtomTemplate()
{
    auto pairs = PrepareAtomTemplate(m_element_templates);
    if (pairs.first.size() == 0 || pairs.second.size() == 0) {
        std::cout << "Templates are empty, maybe try different elements" << std::endl;
        return;
    }
}

void RMSDDriver::HeavyTemplate()
{
    if (!m_silent)
        std::cout << "Prepare heavy atom template structure:" << std::endl;

    PrepareHeavyTemplate();
}
Matrix RMSDDriver::OptimiseRotation(const Eigen::Matrix3d& rotation)
{
    Eigen::Vector3d euler = rotation.eulerAngles(0, 1, 2);

    LBFGSRotation opt(3);
    opt.m_reference = m_reference.getGeometry();
    opt.m_target = m_target.getGeometry();
    opt.m_reference_atoms = m_reference.Atoms();
    opt.m_target_atoms = m_target.Atoms();
    opt.m_costmatrix = m_costmatrix;

    int iteration;
    int maxiter = 1000;
    LBFGSParam<double> param;
    LBFGSSolver<double> solver(param);
    Vector vector(3);
    for (int i = 0; i < 3; ++i)
        vector(i) = euler(i);

    double fx;
    int converged = solver.InitializeSingleSteps(opt, vector, fx);
    Vector v_old = vector;
    for (iteration = 1; iteration <= maxiter; ++iteration) {
        try {
            solver.SingleStep(opt, vector, fx);
        } catch (const std::logic_error& error_result) {
            break;
        } catch (const std::runtime_error& error_result) {
            break;
        }

        if ((vector - v_old).norm() < 1e-5)
            break;
        v_old = vector;
    }
    Eigen::Matrix3d n;
    n = Eigen::AngleAxisd(vector[0], Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(vector[1], Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(vector[2], Eigen::Vector3d::UnitZ());

    return n;
}

void RMSDDriver::TemplateFree()
{
    auto R = GetOperateVectors(m_reference, m_target);
    Eigen::Matrix3d rotation = R.first;

    // Eigen::Matrix3d n = OptimiseRotation(R.first);
    for (int i = 0; i < 360; i += 60) {
        Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
        m(0, 0) = 1;
        m(1, 1) = 1;
        m(2, 2) = 1;
        m = m * Eigen::AngleAxisd(i, Eigen::Vector3d::UnitX());
        // std::cout << m << std::endl;

        auto blob = MakeCostMatrix(OptimiseRotation(m));
        //      std::cout << blob.first << " " << std::endl;
        m_cost_limit = blob.first;

        //    if (m_method == 4)
        //        blob.first = 1;
        InsertRotation(blob);
    }
    for (int j = 0; j < 360; j += 60) {
        Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
        m(0, 0) = 1;
        m(1, 1) = 1;
        m(2, 2) = 1;
        m = m * Eigen::AngleAxisd(j, Eigen::Vector3d::UnitY());
        // std::cout << m << std::endl;

        auto blob = MakeCostMatrix(OptimiseRotation(m));
        //     std::cout << blob.first << " " << std::endl ;

        m_cost_limit = blob.first;

        //   if (m_method == 4)
        //       blob.first = 1;
        InsertRotation(blob);
    }
    for (int k = 0; k < 360; k += 60) {
        Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
        m(0, 0) = 1;
        m(1, 1) = 1;
        m(2, 2) = 1;
        m = m * Eigen::AngleAxisd(k, Eigen::Vector3d::UnitZ());
        // std::cout << m << std::endl;

        auto blob = MakeCostMatrix(OptimiseRotation(m));
        //    std::cout << blob.first << " " << std::endl ;

        m_cost_limit = blob.first;

        //   if (m_method == 4)
        //       blob.first = 1;
        InsertRotation(blob);
    }
}

void RMSDDriver::InsertRotation(std::pair<double, Matrix>& rotation)
{
    for (auto i : m_prepared_cost_matrices) {
        if ((rotation.second - i.second).cwiseAbs().sum() / double(i.second.rows() * i.second.rows()) < 10)
            return;
    }
    m_prepared_cost_matrices.insert(rotation);
}

bool RMSDDriver::TemplateReorder()
{
    auto fragments = CheckFragments();

    if (fragments.first == -1 || fragments.second == -1)
        return false;

    m_fragment = -1;
    m_fragment_target = -1;
    m_fragment_reference = -1;

    Molecule reference_mol = m_reference.getFragmentMolecule(fragments.first);
    Molecule target_mol = m_target.getFragmentMolecule(fragments.second);

    auto operators = GetOperateVectors(fragments.first, fragments.second);
    Eigen::Matrix3d R = operators.first;

    Geometry cached_reference = m_reference.getGeometryByFragment(fragments.first, m_protons);
    Geometry cached_target = m_target.getGeometryByFragment(fragments.second, m_protons);

    Geometry ref = GeometryTools::TranslateMolecule(m_reference, GeometryTools::Centroid(cached_reference), Position{ 0, 0, 0 });
    Geometry tget = GeometryTools::TranslateMolecule(m_target, GeometryTools::Centroid(cached_target), Position{ 0, 0, 0 });

    // Eigen::MatrixXd tar = tget.transpose();

    Geometry rotated = tget * R;

    Molecule ref_mol = m_reference;
    ref_mol.setGeometry(ref);

    Molecule tar_mol = m_target;
    tar_mol.setGeometry(rotated);

    auto blob = MakeCostMatrix(cached_reference, rotated);
    m_cost_limit = blob.first;
    // if (!m_silent)
    // if (m_method == 4)
    //    blob.first = 1;

    InsertRotation(blob);
    return true;
}

std::pair<std::vector<int>, std::vector<int>> RMSDDriver::PrepareHeavyTemplate()
{
    Molecule reference;
    std::vector<int> reference_indicies, target_indicies;
    for (std::size_t i = 0; i < m_reference.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_reference.Atom(i);
        if (atom.first != 1) {
            reference.addPair(atom);
            reference_indicies.push_back(i);
        }
    }

    Molecule target;
    for (std::size_t i = 0; i < m_target.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_target.Atom(i);
        if (atom.first != 1) {
            target.addPair(atom);
            target_indicies.push_back(i);
        }
    }

    Molecule cached_reference_mol = m_reference;
    Molecule cached_target_mol = m_target;

    m_reference = reference;
    m_target = target;

    m_init_count = m_heavy_init;

    ReorderIncremental();
    std::vector<std::vector<int>> transformed_rules;
    for (int i = 0; i < m_stored_rules.size(); ++i) {
        std::vector<int> tmp;
        for (auto index : m_stored_rules[i])
            tmp.push_back(target_indicies[index]);
        transformed_rules.push_back(tmp);
    }
    m_stored_rules = transformed_rules;
    std::vector<int> target_indices = m_reorder_rules;

    m_reference = cached_reference_mol;
    m_target = cached_target_mol;

    return std::pair<std::vector<int>, std::vector<int>>(reference_indicies, target_indices);
}

std::pair<std::vector<int>, std::vector<int>> RMSDDriver::PrepareAtomTemplate(int templateatom)
{
    Molecule reference;
    std::vector<int> reference_indicies, target_indicies;
    for (std::size_t i = 0; i < m_reference.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_reference.Atom(i);
        if (atom.first == templateatom) {
            reference.addPair(atom);
            reference_indicies.push_back(i);
        }
    }

    Molecule target;
    for (std::size_t i = 0; i < m_target.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_target.Atom(i);
        if (atom.first == templateatom) {
            target.addPair(atom);
            target_indicies.push_back(i);
        }
    }
    if (target.AtomCount() == 0 || reference.AtomCount() == 0) {
        std::cout << " Template list is empty, try different elements please" << std::endl;
        exit(0);
    }
    Molecule cached_reference_mol = m_reference;
    Molecule cached_target_mol = m_target;

    m_reference = reference;
    m_target = target;

    m_init_count = m_heavy_init;

    ReorderIncremental();
    m_reference = cached_reference_mol;
    m_target = cached_target_mol;
    for (int i = 0; i < m_stored_rules.size(); ++i) {
        std::vector<int> tmp;
        for (auto index : m_stored_rules[i])
            tmp.push_back(target_indicies[index]);

        auto OperateVectors = GetOperateVectors(reference_indicies, tmp);
        auto result = MakeCostMatrix(OperateVectors.first);
        // result.first = m_prepared_cost_matrices.size();

        InsertRotation(result);
    }

    for (const auto& indices : m_intermedia_rules) {
        auto result = MakeCostMatrix(reference_indicies, indices);
        // result.first = m_prepared_cost_matrices.size();
        InsertRotation(result);
    }

    return std::pair<std::vector<int>, std::vector<int>>(reference_indicies, reference_indicies);
}

std::pair<std::vector<int>, std::vector<int>> RMSDDriver::PrepareAtomTemplate(const std::vector<int>& templateatom)
{
    Molecule reference;
    std::vector<int> reference_indicies, target_indicies;
    for (std::size_t i = 0; i < m_reference.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_reference.Atom(i);
        if (std::find(templateatom.begin(), templateatom.end(), atom.first) != templateatom.end()) {
            reference.addPair(atom);
            reference_indicies.push_back(i);
        }
    }

    Molecule target;
    for (std::size_t i = 0; i < m_target.AtomCount(); ++i) {
        std::pair<int, Position> atom = m_target.Atom(i);
        if (std::find(templateatom.begin(), templateatom.end(), atom.first) != templateatom.end()) {
            target.addPair(atom);
            target_indicies.push_back(i);
        }
    }

    Molecule cached_reference_mol = m_reference;
    Molecule cached_target_mol = m_target;

    m_reference = reference;
    m_target = target;

    m_init_count = m_heavy_init;

    ReorderIncremental();

    m_reference = cached_reference_mol;
    m_target = cached_target_mol;

    std::vector<std::vector<int>> transformed_rules;
    for (int i = 0; i < m_stored_rules.size(); ++i) {
        std::vector<int> tmp;
        for (auto index : m_stored_rules[i])
            tmp.push_back(target_indicies[index]);
        transformed_rules.push_back(tmp);
        m_intermedia_rules.push_back(tmp);
        auto OperateVectors = GetOperateVectors(reference_indicies, tmp);
        auto result = MakeCostMatrix(OperateVectors.first);
        // result.first = m_prepared_cost_matrices.size() + 1;

        InsertRotation(result);
    }

    for (const auto& indices : m_intermedia_rules) {
        auto result = MakeCostMatrix(reference_indicies, indices);
        // result.first = m_prepared_cost_matrices.size();
        InsertRotation(result);
    }
    m_stored_rules = transformed_rules;
    std::vector<int> target_indices = m_reorder_rules;

    return std::pair<std::vector<int>, std::vector<int>>(reference_indicies, target_indices);
}

std::pair<std::vector<int>, std::vector<int>> RMSDDriver::PrepareDistanceTemplate() // const std::vector<int>& templateatom)
{
    if (!m_silent)
        std::cout << "Start Prepare Template" << std::endl;
    RunTimer time;
    std::map<double, std::pair<int, int>> m_distance_reference, m_distance_target;

    Position ref_centroid = m_reference.Centroid(true);
    Position tar_centroid = m_target.Centroid(true);

    Geometry reference_geometry = GeometryTools::TranslateMolecule(m_reference, m_reference.Centroid(true), Position{ 0, 0, 0 });
    Geometry target_geometry = GeometryTools::TranslateMolecule(m_target, m_target.Centroid(true), Position{ 0, 0, 0 });

    std::map<double, int> distance_reference, distance_target;
    Molecule reference;
    std::vector<int> reference_indicies, target_indicies;
    for (std::size_t i = 0; i < m_reference.AtomCount() && i < m_target.AtomCount(); ++i) {

        double distance = (m_reference.Atom(i).second - ref_centroid).norm();
        distance_reference[distance] = i;

        distance = (m_target.Atom(i).second - tar_centroid).norm();
        distance_target[distance] = i;

        for (int j = i + 1; j < m_reference.AtomCount() && j < m_target.AtomCount(); ++j) {
            double distance1 = (m_reference.Atom(i).second - m_reference.Atom(j).second).norm();
            m_distance_reference[distance1] = std::pair<int, int>(i, j);

            double distance2 = (m_target.Atom(i).second - m_target.Atom(j).second).norm();
            m_distance_target[distance2] = std::pair<int, int>(i, j);
        }
    }

    Molecule target;
    int i = 0;

    auto ref_end = m_distance_reference.cend();
    auto tar_end = m_distance_target.cend();

    while (reference_indicies.size() < m_limit) {
        ref_end--;
        tar_end--;
        std::pair<int, Position> atom_r1 = m_reference.Atom(ref_end->second.first);
        std::pair<int, Position> atom_r2 = m_reference.Atom(ref_end->second.second);

        std::pair<int, Position> atom_t1 = m_target.Atom(tar_end->second.first);
        std::pair<int, Position> atom_t2 = m_target.Atom(tar_end->second.second);

        if (atom_r1.first == atom_t1.first && std::find(reference_indicies.begin(), reference_indicies.end(), ref_end->second.first) == reference_indicies.end() && std::find(target_indicies.begin(), target_indicies.end(), tar_end->second.first) == target_indicies.end()) {
            {
                reference.addPair(atom_r1);
                reference_indicies.push_back(ref_end->second.first);

                target.addPair(atom_t1);
                target_indicies.push_back(tar_end->second.first);
            }
        }
        if (atom_r2.first == atom_t2.first && std::find(reference_indicies.begin(), reference_indicies.end(), ref_end->second.second) == reference_indicies.end() && std::find(target_indicies.begin(), target_indicies.end(), tar_end->second.second) == target_indicies.end()) {
            {
                reference.addPair(atom_r2);
                reference_indicies.push_back(ref_end->second.second);

                target.addPair(atom_t2);
                target_indicies.push_back(tar_end->second.second);
            }
        }
        ++i;
    }

    Molecule cached_reference_mol = m_reference;
    Molecule cached_target_mol = m_target;

    m_reference = reference;
    m_target = target;

    m_init_count = m_heavy_init;

    ReorderIncremental();

    m_reference = cached_reference_mol;
    m_target = cached_target_mol;
    for (const auto& indices : m_intermedia_rules) {
        auto result = MakeCostMatrix(reference_indicies, indices);
        InsertRotation(result);
    }
    m_intermedia_rules.clear();

    std::map<double, std::vector<int>> order_rules;
    std::vector<std::vector<int>> transformed_rules;
    for (int i = 0; i < m_stored_rules.size(); ++i) {
        std::vector<int> tmp;
        for (auto index : m_stored_rules[i])
            tmp.push_back(target_indicies[index]);

        auto OperateVectors = GetOperateVectors(reference_indicies, tmp);
        auto result = MakeCostMatrix(OperateVectors.first);
        InsertRotation(result);
    }

    m_stored_rules.clear();
    std::vector<int> target_indices = m_reorder_rules;
    return std::pair<std::vector<int>, std::vector<int>>(reference_indicies, target_indices);
}

std::pair<double, Matrix> RMSDDriver::MakeCostMatrix(const std::vector<int>& permuation)
{
    std::vector<int> first(permuation.size(), 0);
    for (int i = 0; i < permuation.size(); ++i)
        first[i] = i;
    return MakeCostMatrix(std::pair<std::vector<int>, std::vector<int>>(first, permuation));
}

std::pair<double, Matrix> RMSDDriver::MakeCostMatrix(const std::vector<int>& reference, const std::vector<int>& target)
{
    auto operators = GetOperateVectors(reference, target);
    Eigen::Matrix3d R = operators.first;

    Geometry cached_reference = m_reference_centered.getGeometry();
    Geometry cached_target = m_target_centered.getGeometry();

    //  Eigen::MatrixXd tar = cached_target.transpose();

    Geometry rotated = cached_target * R;

    Molecule ref_mol = m_reference;
    ref_mol.setGeometry(cached_reference);

    Molecule tar_mol = m_target;
    tar_mol.setGeometry(rotated);
    return MakeCostMatrix(cached_reference, rotated);
}

std::pair<double, Matrix> RMSDDriver::MakeCostMatrix(const std::pair<std::vector<int>, std::vector<int>>& pair)
{
    return MakeCostMatrix(pair.first, pair.second);
}

std::pair<double, Matrix> RMSDDriver::MakeCostMatrix(const Matrix& rotation)
{
    // Geometry target = m_target.getGeometry().transpose();
    Geometry rotated = m_target.getGeometry() * rotation;
    Geometry reference = m_reference.getGeometry();
    return MakeCostMatrix(reference, rotated);
}

std::pair<double, Matrix> RMSDDriver::MakeCostMatrix(const Geometry& reference, const Geometry& target /*, const std::vector<int> reference_atoms, const std::vector<int> target_atoms*/)
{
    double penalty = 1e23;

    Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(m_reference.AtomCount(), m_reference.AtomCount());
    double sum = 0;
    for (int i = 0; i < m_reference.AtomCount(); ++i) {
        double min = penalty;
        int row_count = 0;
        for (int j = 0; j < m_reference.AtomCount(); ++j) {
            double d = (target.row(j) - reference.row(i)).norm();
            double norm = (target.row(j).norm() - reference.row(i).norm());
            distance(i, j) = Cost(d, norm, m_costmatrix);
            distance(i, j) += penalty * (m_reference.Atom(i).first != m_target.Atom(j).first);
            min = std::min(min, distance(i, j));
        }
        sum += min;
    }
    return std::pair<double, Matrix>(sum, distance);
}

std::pair<double, Matrix> RMSDDriver::MakeCostMatrix(const Geometry& reference, const Geometry& target, const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms, int costmatrix)
{
    double penalty = 1e23;
    Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(reference_atoms.size(), reference_atoms.size());
    double sum = 0;
    for (int i = 0; i < reference_atoms.size(); ++i) {
        double min = penalty;
        int row_count = 0;
        for (int j = 0; j < reference_atoms.size(); ++j) {
            double d = (target.row(j) - reference.row(i)).norm();
            double norm = (target.row(j).norm() - reference.row(i).norm());
            distance(i, j) = Cost(d, norm, costmatrix);

            distance(i, j) += penalty * (reference_atoms[i] != target_atoms[j]);
            min = std::min(min, distance(i, j));
        }
        sum += min;
    }
    return std::pair<double, Matrix>(sum, distance);
}

std::vector<int> RMSDDriver::SolveCostMatrix(Matrix& distance)
{
    std::vector<int> new_order;
    new_order.resize(distance.rows());
    double difference = distance.cwiseAbs().sum();
    int iter = 0;
    for (iter = 0; iter < 10 && difference != 0; ++iter) {
        bool use_c = false;
        int dim = distance.rows();
        if (use_c) {
            double* table = new double[dim * dim];
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    table[i * dim + j] = distance(i, j);
                }
            }
            int* order = new int[dim];
            assign(dim, table, order);
            for (int i = 0; i < dim; ++i)
                new_order[i] = order[i];
            delete (table);
            delete (order);
        } else {
            auto result = MunkressAssign(distance);
            for (int i = 0; i < result.cols(); ++i) {
                for (int j = 0; j < result.rows(); ++j) {
                    if (result(i, j) == 1) {
                        new_order[i] = j;
                        break;
                    }
                }
            }
        }
        auto pair = MakeCostMatrix(new_order);
        difference = (distance - pair.second).cwiseAbs().sum();
        distance = pair.second;
    }
    if (!m_silent)
        std::cout << iter << std::endl;
    return new_order;
}

Molecule RMSDDriver::ApplyOrder(const std::vector<int>& order, const Molecule& mol)
{
    Molecule result;
    for (auto i : order) {
        if (result.AtomCount() >= mol.AtomCount())
            break;
        result.addPair(mol.Atom(i));
    }
    return result;
}

std::pair<int, int> RMSDDriver::CheckFragments()
{
    if (m_fragment != -1)
        return std::pair<int, int>(m_fragment, m_fragment);

    auto ref_mass = m_reference.FragmentMass();
    auto tar_mass = m_target.FragmentMass();

    for (int i = ref_mass.size() - 1; i >= 0; --i) {
        for (int j = tar_mass.size() - 1; j >= 0; --j) {
            if (std::abs(ref_mass[i] - tar_mass[j]) < 1e-5)
                return std::pair<int, int>(i, j);
        }
    }
    return std::pair<int, int>(-1, -1);
}

std::pair<Matrix, Position> RMSDDriver::GetOperateVectors(int fragment_reference, int fragment_target)
{
    Molecule reference_mol = m_reference.getFragmentMolecule(fragment_reference);
    Molecule target_mol = m_target.getFragmentMolecule(fragment_target);
    return GetOperateVectors(reference_mol, target_mol);
}

std::pair<Matrix, Position> RMSDDriver::GetOperateVectors(const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms)
{
    Molecule reference_mol;
    Molecule target_mol;
    if (reference_atoms.size() == target_atoms.size()) {
        for (int i = 0; i < reference_atoms.size(); ++i) {
            reference_mol.addPair(m_reference.Atom(reference_atoms[i]));
            target_mol.addPair(m_target.Atom(target_atoms[i]));
        }
    } else {
        for (int i = 0; i < target_atoms.size(); ++i) {
            reference_mol.addPair(m_reference.Atom(i));
            target_mol.addPair(m_target.Atom(target_atoms[i]));
        }
    }
    return GetOperateVectors(reference_mol, target_mol);
}

std::pair<Matrix, Position> RMSDDriver::GetOperateVectors(const Molecule& reference, const Molecule& target)
{
    Eigen::Matrix3d R = RMSDFunctions::BestFitRotation(reference, target, 1);

    Geometry cached_reference = reference.getGeometry();
    Geometry cached_target = target.getGeometry();

    Position translate = GeometryTools::Centroid(cached_reference) - GeometryTools::Centroid(cached_target);

    return std::pair<Matrix, Position>(R, translate);
}

bool RMSDDriver::MolAlignLib()
{
    m_reference.writeXYZFile("molaign_ref.xyz");
    m_target.writeXYZFile("molalign_tar.xyz");

    FILE* FileOpen;
    std::string command = m_molalign + m_molalignarg + " molaign_ref.xyz   molalign_tar.xyz " + " 2>&1";
    if (!m_silent)
        std::cout << command << std::endl;
    FileOpen = popen(command.c_str(), "r");
    bool ok = true;
    bool rndm = false;
    bool error = false;
    char line[130];

    while (fgets(line, sizeof line, FileOpen)) {
        //  ok = std::string(line).find("RMSD") != std::string::npos;
        rndm = std::string(line).find("random") != std::string::npos;
        error = std::string(line).find("Error") != std::string::npos;
        printf("%s", line);
    }
    pclose(FileOpen);

    if (std::filesystem::exists("aligned.xyz") and !rndm) {
        if (!m_silent) {
            fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nPlease cite the follow research report!\nJ. Chem. Inf. Model. 2023, 63, 4, 1157–1165 - DOI: 10.1021/acs.jcim.2c01187\n\n");
        }
        FileIterator file("aligned.xyz", true);
        m_reference_centered = file.Next();
        m_target_reordered = file.Next();
        m_target_aligned = m_target_reordered;
        m_target = m_target_reordered;
        m_rmsd = CalculateRMSD();
        std::filesystem::remove("aligned.xyz");
    } else {
        if (!rndm && !error) {
            fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "Molalign was not found. Consider getting it from\nhttps://github.com/qcuaeh/molalignlib\nEither adding the location of the binary to the path for executables or append\n-molalignbin /yourpath/molalign\nto your argument list!\n");
            return false;
        }
    }
    if (rndm) {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "molalign has trouble with random numbers, falling back to plain Kuhn-Munkres ...\n");
        return false;
    }
    if (error) {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "molalign has trouble with finding the solution, try to increase tolerance via -molaligntol XX \n");
        return false;
    }
    return true;
}

Geometry RMSDDriver::Gradient() const
{
    return (m_reference.getGeometry() - m_target.getGeometry()) / (RMSD() * m_target.getGeometry().rows());
}
