/*
 * <RMSD calculator for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/molecule.h"
#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "json.hpp"
using json = nlohmann::json;

#include "rmsd.h"

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
                //const auto t = RMSDFunctions::getAligned(m_reference, GeometryTools::TranslateGeometry(target_local.getGeometry(), target_local.Centroid(true), Position{ 0, 0, 0 }), 1);
                const auto t = RMSDFunctions::getAligned(m_reference, target_local.getGeometry(), 1);

                double rmsd_local = RMSDFunctions::getRMSD(m_reference, t);

                if (target_local.AtomCount() <= m_target.AtomCount()) {
                    {
                        match.insert(std::pair<double, int>(rmsd_local, j));
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
    m_force_reorder = Json2KeyWord<bool>(m_defaults, "reorder");
    m_protons = !Json2KeyWord<bool>(m_defaults, "heavy");
    m_silent = Json2KeyWord<bool>(m_defaults, "silent");
    m_intermedia_storage = Json2KeyWord<double>(m_defaults, "storage");
    m_dynamic_center = Json2KeyWord<bool>(m_defaults, "DynamicCenter");

    m_noreorder = Json2KeyWord<bool>(m_defaults, "noreorder");
    m_element = Json2KeyWord<int>(m_defaults, "element");

    std::string method = Json2KeyWord<std::string>(m_defaults, "method");

    if (method.compare("template") == 0)
        m_method = 2;
    else if (method.compare("incr") == 0)
        m_method = 1;
    else if (method.compare("hybrid0") == 0)
        m_method = 3;
    else if (method.compare("hybrid") == 0)
        m_method = 4;
    else if (method.compare("free") == 0)
        m_method = 5;
    else
        m_method = 1;

    std::string order = Json2KeyWord<std::string>(m_defaults, "order");

    std::vector<int> vector = Tools::String2Vector(order);
    if (vector.size() != 0)
        m_reorder_rules = vector;
}

void RMSDDriver::start()
{
    RunTimer timer(false);

    clear();

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


    /*
    if (std::abs(m_reference.Mass() - m_target.Mass()) > 1e-4) {
        bool stop = false;
        if (!m_silent)
            fmt::print("Atom count diveres for both molecules, check for supramolecular fragments, that may match!\nFor now, will take only one pair.");

        auto ref_fragments = m_reference.GetFragments();
        auto tar_fragments = m_target.GetFragments();
        for (int i = 0; i < ref_fragments.size() && !stop; ++i) {
            for (int j = 0; j < tar_fragments.size() && !stop; ++j) {
                if (abs(m_reference.getFragmentMolecule(i).Mass() - m_target.getFragmentMolecule(j).Mass()) < 1e-4) {
                    if (!m_silent) {
                        fmt::print("\n\nOverwriting reference molecule (mass = {0:f}) with fragment {1} (mass = {2:f}).\n", m_reference.Mass(), i, m_reference.getFragmentMolecule(i).Mass());
                        fmt::print("Overwriting target molecule (mass = {0:f}) with fragment {1} (mass = {2:f}).\n\n\n", m_target.Mass(), j, m_target.getFragmentMolecule(j).Mass());
                    }
                    m_reference = m_reference.getFragmentMolecule(i);
                    m_target = m_target.getFragmentMolecule(j);
                    stop = true;*/
    /*
                     * https://stackoverflow.com/questions/9695902/how-to-break-out-of-nested-loops - NoGo2
                     */
    /*}
            }
        }
    }*/

    m_target_aligned = m_target;
    m_reference_aligned.LoadMolecule(m_reference);

    if (m_reference.Atoms() != m_target.Atoms() || m_force_reorder) {
        if (!m_noreorder)
            ReorderMolecule();
    }
    Molecule temp_ref, temp_tar;
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
        } else
            m_rmsd = BestFitRMSD();
    } else {
        if (!m_silent)
            fmt::print("Partial RMSD is calculated, only from those atoms, that match each other.\n\n\n");
        m_rmsd = PartialRMSD(temp_ref, temp_tar);
    }
    /*
    auto terms = IndivRMSD(m_reference, m_target);
    for(const auto  &i:terms)
        std::cout << i << std::endl;
    */
    m_htopo_diff = CompareTopoMatrix(m_reference_aligned.HydrogenBondMatrix(-1, -1), m_target_aligned.HydrogenBondMatrix(-1, -1));
    if (!m_silent) {
        std::cout << "RMSD calculation took " << timer.Elapsed() << " msecs." << std::endl;
        std::cout << "Difference in Topological Hydrogen Bond Matrix is " << m_htopo_diff << std::endl;
    }
    if (m_swap) {
        Molecule reference = m_reference;
        //Molecule reference_reorder = m_reference_reordered;
        Molecule reference_aligned = m_reference_aligned;
        m_reference = m_target;
        m_reference_aligned = m_target_aligned;
        // m_reference_reordered = m_target_reordered;
        m_target = reference;
        //m_target_reordered = reference_reorder;
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
    const auto t = RMSDFunctions::getAligned(reference, target, 1);
    m_reference_aligned.setGeometry(reference);
    m_target_aligned.setGeometry(t);
    rmsd = RMSDFunctions::getRMSD(reference, t);
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
    auto target = GeometryTools::TranslateGeometry(m_target.getGeometry(), GeometryTools::Centroid(m_target.getGeometryByFragment(fragment_target)), Position{ 0, 0, 0 }); //CenterMolecule(target_mol);
    const auto t = RMSDFunctions::applyRotation(target, rotation);
    m_reference_aligned.setGeometry(reference);
    //m_reference.writeXYZFile("ref.xyz");
    m_target_aligned.setGeometry(t);
    //m_target_aligned.writeXYZFile("tar.xyz");
    rmsd = RMSDFunctions::getRMSD(reference, t);

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
    int wake_up = 100;
    std::vector<AtomDef> atoms;
    while (
        m_reorder_reference_geometry.rows() < m_reorder_reference.AtomCount() && m_reorder_reference_geometry.rows() < m_reorder_target.AtomCount() && ((reference_reordered + reference_not_reorordered) <= m_reference.AtomCount())) {
        int thread_count = 0;

        CxxThreadPool* pool = new CxxThreadPool;
        pool->setActiveThreadCount(m_threads);
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

        for (const auto& e : *storage_shelf.data()) {
            RMSDThread* thread = new RMSDThread(m_reorder_target, m_reorder_reference_geometry, e.second, mass, element);
            pool->addThread(thread);
            thread_count++;
        }
        pool->StaticPool();
        pool->setWakeUp(wake_up);
        pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
        int match = 0;
        /* For now, lets just dont start the threads if the current element can not be found in target */
        if (std::find(m_target.Atoms().begin(), m_target.Atoms().end(), element) != m_target.Atoms().end()) {
            pool->StartAndWait();
        }
        LimitedStorage storage_shelf_next(inter_size);
        for (const auto t : pool->Finished()) {
            RMSDThread* thread = static_cast<RMSDThread*>(t);
            for (const auto& item : (*thread->data())) {
                if (thread->Match()) {
                    storage_shelf_next.addItem(item.second, item.first);
                }
                match += thread->Match();
            }
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
            reference_not_reorordered++;
        }
        wake_up = 2 * pool->WakeUp();
        delete pool;
    }
    int count = 0;
    for (const auto& e : *storage_shelf.data()) {
        if (count > 50)
            continue;
        m_stored_rules.push_back(FillMissing(m_reference, e.second));
        count++;
    }
    m_reorder_rules = m_stored_rules[0];
    m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
    m_target = m_target_reordered;
    m_target_aligned = m_target;
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

    double rmsd = CalculateRMSD(m_reference, target);
    m_fragment_reference = tmp_ref;
    m_fragment_target = tmp_tar;
    return rmsd;
}

void RMSDDriver::clear()
{
    m_results.clear();
    m_connectivity.clear();
    m_stored_rules.clear();
}

void RMSDDriver::ProtonDepleted()
{
    if (!m_silent)
        std::cerr << "Will perform calculation on proton depleted structure." << std::endl;

    Molecule reference;
    for(std::size_t i = 0; i < m_reference.AtomCount(); ++i)
    {
        std::pair<int, Position> atom = m_reference.Atom(i);
        if(atom.first != 1)
            reference.addPair(atom);
    }

    Molecule target;
    for(std::size_t i = 0; i < m_target.AtomCount(); ++i)
    {
        std::pair<int, Position> atom = m_target.Atom(i);
        if(atom.first != 1)
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
    Eigen::MatrixXd tar = target.transpose();

    Geometry rotated = tar.transpose() * R;

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
/*
int RMSDDriver::CheckConnectivitiy(const Molecule& mol1, const Molecule& mol2) const
{
    auto connect_1 = mol1.getConnectivtiy(m_scaling);
    auto connect_2 = mol2.getConnectivtiy(m_scaling);

    if (connect_1.size() != connect_2.size())
        return -1;
    int difference = 0;
    int match = 0;
    if (!m_silent)
        if (m_print_intermediate)
            std::cout << std::endl
                      << std::endl
                      << " *** Connectivitiy check will be performed ***" << std::endl
                      << std::endl;
    for (std::size_t i = 0; i < connect_1.size(); ++i) {
        auto target = connect_1[i];
        auto reference = connect_2[i];

        difference += Tools::VectorDifference(reference, target);
        if (!m_silent)
            if (m_print_intermediate)
                std::cout << "vector difference " << Tools::VectorDifference(reference, target) << std::endl;
        if (reference == target) {
            if (!m_silent)
                if (m_print_intermediate) {
                    std::cout << i << " matches. Fine!" << std::endl;
                    for (const auto& i : reference)
                        std::cout << " " << i;
                    std::cout << std::endl;

                    for (const auto& i : target)
                        std::cout << " " << i;
                    std::cout << std::endl;
                }
            match++;
        } else {
            if (!m_silent)
                if (m_print_intermediate) {
                    std::cout << "No match for " << i << std::endl;
                    for (const auto& i : reference)
                        std::cout << " " << i;
                    std::cout << std::endl;

                    for (const auto& i : target)
                        std::cout << " " << i;
                    std::cout << std::endl;
                }
        }
    }
    if (!m_silent)
        if (m_print_intermediate)
            std::cout << std::endl
                      << std::endl
                      << " *** Connectivitiy check done! ***" << std::endl
                      << std::endl;

    return difference;
}

int RMSDDriver::CheckConnectivitiy(const Molecule& mol1) const
{
    auto connect = mol1.getConnectivtiy(m_scaling);

    if (m_connectivity.size() != connect.size())
        return -1;
    if (!m_silent)
        if (m_print_intermediate)
            std::cout << std::endl
                      << std::endl
                      << " *** Connectivitiy check will be performed ***" << std::endl
                      << std::endl;

    int match = 0;
    int difference = 0;

    for (std::size_t i = 0; i < connect.size(); ++i) {
        auto target = connect[i];

        auto reference = m_connectivity.at(i);
        difference += Tools::VectorDifference(reference, target);

        if (!m_silent)
            if (m_print_intermediate)
                std::cout << "vector difference " << Tools::VectorDifference(reference, target) << std::endl;

        if (reference == target) {
            if (m_print_intermediate) {
                std::cout << i << " matches. Fine!" << std::endl;
                for (const auto& i : reference)
                    std::cout << " " << i;
                std::cout << std::endl;

                for (const auto& i : target)
                    std::cout << " " << i;
                std::cout << std::endl;
            }
            match++;
        } else {
            if (m_print_intermediate) {
                std::cout << "No match for " << i << std::endl;
                for (const auto& i : reference)
                    std::cout << " " << i;
                std::cout << std::endl;

                for (const auto& i : target)
                    std::cout << " " << i;
                std::cout << std::endl;
            }
        }
    }
    if (m_print_intermediate)
        std::cout << std::endl
                  << std::endl
                  << " *** Connectivitiy check done! ***" << std::endl
                  << std::endl;

    return difference;
}
*/
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
                    storage.addItem({ i, j }, storage.size());
                if (tmp_reference[0] == elements_target[j] && tmp_reference[1] == elements_target[i])
                    storage.addItem({ j, i }, storage.size());
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
    if (m_reorder_rules.size() != 0) {
        m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
        m_target = m_target_reordered;
        m_target_aligned = m_target;
        return;
    }

    double scaling = 1.5;
    m_connectivity = m_reference.getConnectivtiy(scaling);

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
}

void RMSDDriver::AtomTemplate()
{
    std::cout << "Prepare atom template structure:" << std::endl;

    auto pairs = PrepareAtomTemplate(m_element);
    FinaliseTemplate(pairs);

    m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
    m_target = m_target_reordered;
    m_target_aligned = m_target;
}

void RMSDDriver::HeavyTemplate()
{
    std::cout << "Prepare heavy atom template structure:" << std::endl;

    auto pairs = PrepareHeavyTemplate();
    FinaliseTemplate(pairs);

    m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
    m_target = m_target_reordered;
    m_target_aligned = m_target;
}

void RMSDDriver::TemplateFree()
{
    auto operators = GetOperateVectors(m_reference, m_target);
    Eigen::Matrix3d R = operators.first;

    Geometry cached_reference = m_reference.getGeometry();
    Geometry cached_target = m_target.getGeometry();

    Geometry ref = GeometryTools::TranslateMolecule(m_reference, m_reference.Centroid(true), Position{ 0, 0, 0 });
    Geometry tget = GeometryTools::TranslateMolecule(m_target, m_target.Centroid(true), Position{ 0, 0, 0 });

    Eigen::MatrixXd tar = tget.transpose();

    Geometry rotated = tar.transpose() * R;

    Molecule ref_mol = m_reference;
    ref_mol.setGeometry(ref);
    Molecule tar_mol = m_target;
    tar_mol.setGeometry(rotated);

    std::vector<int> new_order = DistanceReorder(ref_mol, tar_mol);
    m_target_reordered = ApplyOrder(new_order, m_target);
    m_target = m_target_reordered;
    m_target_aligned = m_target;
}

void RMSDDriver::FinaliseTemplate(std::pair<std::vector<int>, std::vector<int>> pairs)
{
    std::vector<int> tmp;
    for (int j = 0; j < m_reorder_rules.size(); ++j)
        tmp.push_back(j);

    Molecule target = m_target;
    std::map<double, std::vector<int>> local_results;
    for (int outer = 0; outer < m_stored_rules.size() && outer < 10; ++outer) {
        pairs.second = m_stored_rules[outer];

        for (int i = 0; i < 5; ++i) {
            auto result = AlignByVectorPair(pairs);
            m_reorder_rules = result;
            m_target_reordered = ApplyOrder(m_reorder_rules, target);
            local_results.insert(std::pair<double, std::vector<int>>(Rules2RMSD(m_reorder_rules), m_reorder_rules));

            result = AlignByVectorPair(tmp, m_reorder_rules);
            if (m_reorder_rules == result)
                break;
            m_reorder_rules = result;
            m_target_reordered = ApplyOrder(m_reorder_rules, target);
            local_results.insert(std::pair<double, std::vector<int>>(Rules2RMSD(m_reorder_rules), m_reorder_rules));
        }
    }
    m_stored_rules.clear();
    for (const auto& i : local_results) {
        m_stored_rules.push_back(i.second);
    }
    m_reorder_rules = local_results.begin()->second;
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
    m_reference = cached_reference_mol;
    m_target = cached_target_mol;

    return std::pair<std::vector<int>, std::vector<int>>(reference_indicies, target_indices);
}

std::vector<int> RMSDDriver::AlignByVectorPair(std::vector<int> first, std::vector<int> second)
{
    auto operators = GetOperateVectors(first, second);
    Eigen::Matrix3d R = operators.first;

    Geometry cached_reference = m_reference.getGeometry(first, m_protons);
    Geometry cached_target = m_target.getGeometry(second, m_protons);

    Geometry ref = GeometryTools::TranslateMolecule(m_reference, m_reference.Centroid(true), Position{ 0, 0, 0 });
    Geometry tget = GeometryTools::TranslateMolecule(m_target, m_target.Centroid(true), Position{ 0, 0, 0 });

    Eigen::MatrixXd tar = tget.transpose();

    Geometry rotated = tar.transpose() * R;

    Molecule ref_mol = m_reference;
    ref_mol.setGeometry(ref);

    Molecule tar_mol = m_target;
    tar_mol.setGeometry(rotated);

    return DistanceReorder(ref_mol, tar_mol);
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
    for (int i = 0; i < reference_atoms.size(); ++i) {
        reference_mol.addPair(m_reference.Atom(reference_atoms[i]));
        target_mol.addPair(m_target.Atom(target_atoms[i]));
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

    Eigen::MatrixXd tar = tget.transpose();

    Geometry rotated = tar.transpose() * R;

    Molecule ref_mol = m_reference;
    ref_mol.setGeometry(ref);

    Molecule tar_mol = m_target;
    tar_mol.setGeometry(rotated);

    m_reorder_rules = DistanceReorder(ref_mol, tar_mol);
    m_target_reordered = ApplyOrder(m_reorder_rules, m_target);
    m_target = m_target_reordered;
    return true;
}

std::vector<int> RMSDDriver::DistanceReorder(const Molecule& reference, const Molecule& target)
{
    std::vector<int> new_order(m_reference.AtomCount(), -1), done_ref, done_tar;

    while (done_ref.size() < m_reference.AtomCount()) {
        double distance = 1e10;
        int match_reference = 0;
        int match_target = 0;
        for (int i = 0; i < target.AtomCount(); ++i) {
            if (std::find(done_tar.begin(), done_tar.end(), i) != done_tar.end())
                continue;
            for (int j = 0; j < reference.AtomCount(); ++j) {
                if (std::find(done_ref.begin(), done_ref.end(), j) != done_ref.end())
                    continue;

                if (target.Atom(i).first != reference.Atom(j).first)
                    continue;

                const double local_distance = GeometryTools::Distance(target.Atom(i).second, reference.Atom(j).second);
                if (local_distance <= distance) {
                    distance = local_distance;
                    match_target = i;
                    match_reference = j;
                }
            }
        }
        new_order[match_target] = match_reference;
        done_tar.push_back(match_target);
        done_ref.push_back(match_reference);
    }
    return new_order;
}
