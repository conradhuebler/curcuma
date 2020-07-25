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

#include "src/core/molecule.h"
#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "json.hpp"
using json = nlohmann::json;

#include "rmsd.h"

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
    m_fragment_reference = Json2KeyWord<int>(m_defaults, "fragment");
    m_fragment_target = Json2KeyWord<int>(m_defaults, "fragment");
    m_fragment = Json2KeyWord<int>(m_defaults, "fragment");
    m_initial_fragment = Json2KeyWord<int>(m_defaults, "init");
    m_pt = Json2KeyWord<int>(m_defaults, "pt");
    m_force_reorder = Json2KeyWord<bool>(m_defaults, "reorder");
    m_protons = !Json2KeyWord<bool>(m_defaults, "heavy");
    m_silent = Json2KeyWord<bool>(m_defaults, "silent");
    m_intermedia_storage = Json2KeyWord<double>(m_defaults, "storage");
    m_noreorder = Json2KeyWord<bool>(m_defaults, "noreorder");
    std::string method = Json2KeyWord<std::string>(m_defaults, "method");

    if (method.compare("template") == 0)
        m_method = 2;
    else
        m_method = 1;
}

void RMSDDriver::start()
{
    RunTimer timer(false);

    if (m_target.AtomCount() > m_reference.AtomCount()) {
        Molecule reference = m_reference;
        m_reference = m_target;
        m_target = reference;
    }
    m_intermedia_storage = 1;
    clear();

    if (m_initial_fragment != -1 && m_initial.size() == 0)
        m_initial = m_reference.GetFragments()[m_initial_fragment];

    if (m_initial.size())
        InitialiseOrder();

    m_reference.InitialiseConnectedMass(1.5, m_protons);
    m_target.InitialiseConnectedMass(1.5, m_protons);

    if(m_protons == false)
        ProtonDepleted();

    if (m_fragment_reference < -1 || m_fragment_reference > m_reference.GetFragments(m_scaling).size()) {
        m_fragment_reference = -1;
    }
    if (m_fragment_target < -1 || m_fragment_target > m_target.GetFragments(m_scaling).size()) {
        m_fragment_target = -1;
    }

    if ((m_reference.AtomCount() != m_target.AtomCount()) && (m_fragment_target != -1) && (m_fragment_reference != -1)) {
        if (m_reference.getFragmentMolecule(m_fragment_target).AtomCount() == m_target.getFragmentMolecule(m_fragment_target).AtomCount())
            m_partial_rmsd = true;
        else {
            std::cout << "Nothing fitting at all." << std::endl;
            m_noreorder = false;
        }
    }

    if (m_reference.Atoms() == m_target.Atoms())
        m_rmsd_raw = CalculateRMSD(m_reference, m_target);

    if ((m_reference.Atoms() != m_target.Atoms() || ForceReorder()) && m_noreorder == false) {
        if (!m_silent)
            std::cerr << "Molecules are different. What now?\n";

        if (m_reference.AtomCount() == m_target.AtomCount()) {
            if (!m_silent) {
                std::cerr << "Try to reorder the structures?\n";
                std::cerr << "Initial RMSD is " << m_rmsd_raw << std::endl;
            }

            ReorderMolecule();
        } else {
            ReorderMolecule();
            if (!m_silent)
                std::cout << "RMSD calculation took " << timer.Elapsed() << " msecs." << std::endl;
            return;
        }
    } else
        m_target_reordered = m_target;

    Molecule *reference = new Molecule, *target = new Molecule;
    //m_target_reordered.writeXYZFile("tar.xyz");
    //m_reference.writeXYZFile("ref.xyz");
    m_rmsd = CalculateRMSD(m_reference, m_target_reordered, reference, target);

    m_reference_aligned.LoadMolecule(reference);
    m_target_aligned.LoadMolecule(target);
    //m_target_aligned.writeXYZFile("last_align.xyz");
    m_htopo_diff = CompareTopoMatrix(m_reference_aligned.HydrogenBondMatrix(-1, -1), m_target_aligned.HydrogenBondMatrix(-1, -1));
    if (!m_silent) {
        std::cout << "RMSD calculation took " << timer.Elapsed() << " msecs." << std::endl;
        std::cout << "Difference in Topological Hydrogen Bond Matrix is " << m_htopo_diff << std::endl;
    }
    delete reference;
    delete target;
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
    m_storage.clear();
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

Eigen::Matrix3d RMSDDriver::BestFitRotationShort(const Geometry& reference_mol, const Geometry& target_mol) const
{
    /* The rmsd kabsch algorithmn was adopted from here:
 * https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
 * The specific git commit was
 * https://github.com/oleg-alexandrov/projects/blob/e7b1eb7a4d83d41af563c24859072e4ddd9b730b/eigen/Kabsch.cpp
 */

    Geometry target = CenterMolecule(target_mol);

    Eigen::MatrixXd ref = reference_mol.transpose();
    Eigen::MatrixXd tar = target.transpose();

    Eigen::MatrixXd Cov = ref * tar.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
        d = 1 * 1.0;
    else
        d = 1 * -1.0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;

    return svd.matrixV() * I * svd.matrixU().transpose();
}

Eigen::Matrix3d RMSDDriver::BestFitRotation(const Geometry& reference_mol, const Geometry& target_mol, int factor) const
{
    /* The rmsd kabsch algorithmn was adopted from here:
 * https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
 * The specific git commit was
 * https://github.com/oleg-alexandrov/projects/blob/e7b1eb7a4d83d41af563c24859072e4ddd9b730b/eigen/Kabsch.cpp
 */

    Geometry reference = CenterMolecule(reference_mol);
    Geometry target = CenterMolecule(target_mol);

    Eigen::MatrixXd ref = reference.transpose();
    Eigen::MatrixXd tar = target.transpose();

    Eigen::MatrixXd Cov = ref * tar.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
        d = factor * 1.0;
    else
        d = factor * -1.0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;

    return svd.matrixV() * I * svd.matrixU().transpose();
}

Eigen::Matrix3d RMSDDriver::BestFitRotation(const Molecule& reference_mol, const Molecule& target_mol, int factor) const
{
    const Geometry ref = reference_mol.getGeometryByFragment(m_fragment_reference, m_protons);
    const Geometry tar = target_mol.getGeometryByFragment(m_fragment_target, m_protons);
    return BestFitRotation(ref, tar, factor);
}

double RMSDDriver::CalculateShortRMSD(const Geometry& reference_mol, const Molecule& target_mol) const
{
    const Geometry tar = CenterMolecule(target_mol.getGeometry());

    Eigen::Matrix3d R = BestFitRotationShort(reference_mol, tar);

    double rmsd = 0;

    Eigen::MatrixXd target_transposed = tar.transpose();

    Geometry rotated = target_transposed.transpose() * R;
    for (int i = 0; i < rotated.rows(); ++i) {
        rmsd += (rotated(i, 0) - reference_mol(i, 0)) * (rotated(i, 0) - reference_mol(i, 0)) + (rotated(i, 1) - reference_mol(i, 1)) * (rotated(i, 1) - reference_mol(i, 1)) + (rotated(i, 2) - reference_mol(i, 2)) * (rotated(i, 2) - reference_mol(i, 2));
    }
    rmsd /= double(rotated.rows());

    return sqrt(rmsd);
}

double RMSDDriver::CalculateRMSD(const Molecule& reference_mol, const Molecule& target_mol, Molecule* ret_ref, Molecule* ret_tar, int factor) const
{
    if (reference_mol.AtomCount() != target_mol.AtomCount() && m_partial_rmsd == false)
        return 10;

    Eigen::Matrix3d R = BestFitRotation(reference_mol, target_mol, factor);

    double rmsd = 0;
    // int fragment_reference = m_fragment_reference;
    //  int fragment_target = m_fragment_target;

    Geometry reference;
    Geometry target;

    if (!m_partial_rmsd) {
        reference = GeometryTools::TranslateGeometry(reference_mol.getGeometry(), GeometryTools::Centroid(m_reference.getGeometryByFragment(m_fragment_reference)), Position{ 0, 0, 0 }); // CenterMolecule(reference_mol);
        target = GeometryTools::TranslateGeometry(target_mol.getGeometry(), GeometryTools::Centroid(target_mol.getGeometryByFragment(m_fragment_target)), Position{ 0, 0, 0 }); //CenterMolecule(target_mol);
    } else {
        reference = GeometryTools::TranslateGeometry(reference_mol.getGeometryByFragment(m_fragment_reference), GeometryTools::Centroid(m_reference.getGeometryByFragment(m_fragment_reference)), Position{ 0, 0, 0 }); // CenterMolecule(reference_mol);
        target = GeometryTools::TranslateGeometry(target_mol.getGeometryByFragment(m_fragment_target), GeometryTools::Centroid(target_mol.getGeometryByFragment(m_fragment_target)), Position{ 0, 0, 0 }); //CenterMolecule(target_mol);
    }
    Eigen::MatrixXd tar = target.transpose();


    Geometry rotated = tar.transpose()*R;
    for(int i = 0; i < rotated.rows(); ++i)
    {
        rmsd += (rotated(i, 0) - reference(i, 0))*(rotated(i, 0) - reference(i, 0)) +
                      (rotated(i, 1) - reference(i, 1))*(rotated(i, 1) - reference(i, 1)) +
                      (rotated(i, 2) - reference(i, 2))*(rotated(i, 2) - reference(i, 2));
    }
    rmsd /= double(rotated.rows());

    if (ret_tar != nullptr) {
        ret_tar->LoadMolecule(target_mol);
        ret_tar->setGeometry(rotated);
        // ret_tar->writeXYZFile("tar2.xyz");
    }

    if (ret_ref != nullptr) {
        ret_ref->LoadMolecule(reference_mol);
        ret_ref->setGeometry(reference);
        // ret_ref->writeXYZFile("ref2.xyz");
    }
    //m_fragment_reference = fragment_reference;
    return sqrt(rmsd);
}

std::vector<double> RMSDDriver::IndivRMSD(const Molecule& reference_mol, const Molecule& target_mol, int factor) const
{
    Eigen::Matrix3d R = BestFitRotation(reference_mol, target_mol, factor);

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

void RMSDDriver::InitialisePair()
{
    Molecule reference;
    int index = 0;
    std::vector<int> elements = m_reference.Atoms();
    /*
    std::vector<std::vector<int>> fragments = m_reference.GetFragments();
    for(int i = 0; i < fragments.size(); ++i)
    {
        int index = 0;
        for(int j = 0; j < fragments[i].size() && index < 1; ++j)
        {
            if (elements[fragments[i][j]] != 1) {
                reference.addPair(m_reference.Atom(fragments[i][j]));
                index++;
            }
        }
    }
    */

    if (m_initial.size() == 0) {
        for (int i = 0; i < m_reference.AtomCount() && index < 2; i++) {
            if (elements[i] != 1) {
                reference.addPair(m_reference.Atom(i));
                index++;
            }
        }

        std::vector<int> elements_target = m_target.Atoms();
        std::vector<int> tmp_reference = reference.Atoms();

        for (int i = 0; i < m_target.AtomCount(); ++i) {
            if (m_target.Atom(i).first == 1) // Skip first atom if Proton
                continue;
            for (int j = i + 1; j < m_target.AtomCount(); ++j) {
                if (m_target.Atom(j).first == 1) // Skip second atom if Proton
                    continue;
                if (tmp_reference[0] == elements_target[i] && tmp_reference[1] == elements_target[j])
                    m_intermediate_results.push({ i, j });
                if (tmp_reference[0] == elements_target[j] && tmp_reference[1] == elements_target[i])
                    m_intermediate_results.push({ j, i });
            }
        }
    } else {
        std::vector<int> start;
        for (int i = 0; i < m_init_count; ++i)
            start.push_back(i);
        m_intermediate_results.push(start);
    }
}

void RMSDDriver::ReorderMolecule()
{
    double scaling = 1.5;
    m_connectivity = m_reference.getConnectivtiy(scaling);

    if (m_method == 1)
        ReorderStraight();
    else if (m_method == 2)
        TemplateReorder();

    FinaliseReorder();
}

void RMSDDriver::FinaliseReorder()
{
    if (m_results.size()) {
        m_reorder_rules = m_results.begin()->second;
    } else
        m_stored_rules.push_back(m_reorder_rules);

    Molecule mol2 = ApplyOrder(m_reorder_rules, m_target);
    m_rmsd = CalculateRMSD(m_reference, mol2);
    m_target_reordered = mol2;
    // m_target_reordered.writeXYZFile("lalala3.xyz");
    for (const auto& element : m_results) {
        m_stored_rules.push_back(element.second);
    }
    /*
        if (CheckConnectivitiy(m_reference, mol2) == m_pt) {
            if (!m_silent)
                std::cout << "Found fitting solution, taking ... " << std::endl;
            m_target_reordered = mol2;

        }*/

    /*
    for (const auto& element : m_results) {
        Molecule mol2;
        for (int i = 0; i < element.second.size(); i++) {
            mol2.addPair(m_target.Atom(element.second[i]));
        }
        m_reorder_rules = element.second;
        m_rmsd = CalculateRMSD(m_reference, mol2);
        if (count == 0) // store the best result in any way
            m_target_reordered = mol2;

        if (CheckConnectivitiy(m_reference, mol2) == m_pt) {
            if (!m_silent)
                std::cout << "Found fitting solution, taking ... " << std::endl;
            m_target_reordered = mol2;
            break;
        }
        count++;
        break;
    }
    */
}

void RMSDDriver::ReorderStraight()
{
    int inter_size = m_reference.AtomCount() * (m_reference.AtomCount() - 1) * m_intermedia_storage;
    m_storage = std::vector<IntermediateStorage>(m_reference.AtomCount() - 1, IntermediateStorage(inter_size));

    double scaling = 1.5;

    /* Lets initialise a molecule with two atom from the reference */
    InitialisePair();

    m_rmsd = CalculateRMSD(m_reference, m_target) / double(m_reference.AtomCount());

    while (m_intermediate_results.size()) {
        std::vector<int> inter = m_intermediate_results.front();
        //std::cout << inter.size() << std::endl;
        //if(inter.size() <  m_target.AtomCount())
        SolveIntermediate(inter);
        m_intermediate_results.pop();
    }
    int i = 0;
    int next = 0;
    for (; i < m_storage.size(); ++i) {
        if (!m_silent) {
            /*int ct = 0;
             for (const auto& element : (*m_storage[i].data())) {
                {
                    ct++;
                    std::cout << " " << element.first << " ";
                    if(ct == 10)
                        break;
                }
            }*/
            std::cout << double(i) / double(m_reference.AtomCount()) * 100 << " % done " << std::endl;
        }
        // std::cout << (*m_storage[i].data()).size() << " size of block ";
        int counter = 0;
        // std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        for (const auto& element : (*m_storage[i].data())) {
            //   std::cout << element.first << "(" << m_threshold << ")" << std::endl;
            if (element.first < m_threshold) {
                /*if(i > 34)
                {
                    if(counter < 20)
                        std::cout << element.first << " " << Tools::Vector2String(element.second) << std::endl;
                }*/
                counter++;
                if (!SolveIntermediate(element.second) && !next)
                    next = i;
            }
        }
        //std::cout << " " << counter << " results used " << (*m_storage[i].data()).size() - counter << " ignored (" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() << ")" << std::endl;

        //std::sort(m_last_rmsd.begin(), m_last_rmsd.end());
        //m_threshold = (2*Tools::median(m_last_rmsd) +  m_last_rmsd[0])/3.0;
        m_threshold = Tools::mean(m_last_rmsd);
        // std::cout << m_threshold << std::endl;
        if ((*m_storage[i].data()).size()) {
            m_last_rmsd.clear();
        }
    }
    if (next)
        i = next;
    for (; i < m_storage.size(); ++i) {
        //std::cout << "next " << i << std::endl;
        if (!m_silent) {
            /*int ct = 0;
             for (const auto& element : (*m_storage[i].data())) {
                {
                    ct++;
                    std::cout << " " << element.first << " ";
                    if(ct == 10)
                        break;
                }
            }*/
            //std::cout << double(i) / double(m_reference.AtomCount()) * 100 << " % done " << std::endl;
        }
        for (const auto& element : (*m_storage[i].data())) {
            //if(element.first < m_threshold)
            SolveIntermediate(element.second);
        }

        m_threshold = Tools::median(m_last_rmsd);
        //std::cout << m_threshold << std::endl;
        if ((*m_storage[i].data()).size()) {
            m_last_rmsd.clear();
        }
    }

    std::vector<int> intermediate;
    m_print_intermediate = true;

    if (m_results.size() == 0) {
        for (std::size_t index = m_storage.size() - 1; index >= 0; --index) {

            if (m_storage[index].data()->size()) {
                int count = 0;
                for (const auto& element : (*m_storage[index].data())) {
                    std::cout << element.first;
                    m_rmsd = element.first;
                    for (int l : element.second)
                        std::cout << " " << l;
                    std::cout << std::endl;
                    ReconstructTarget(element.second);
                    return;
                    count++;
                    //if(count)
                    //    break;
                }
                break;
            }
        }
    }
    m_print_intermediate = false;
    return;
}

bool RMSDDriver::SolveIntermediate(std::vector<int> intermediate, bool fast)
{
    if (m_reference_reordered + intermediate.size() >= m_reference.AtomCount())
        return false;
    //std::cout << "Reference reordered " << m_reference_reordered << std::endl;
    Molecule reference;
    Molecule target;
    //std::vector<int> elements_target = m_target.Atoms();

    for (int i = 0; i < intermediate.size(); i++) {
        reference.addPair(m_reference.Atom(i));
        target.addPair(m_target.Atom(intermediate[i]));
    }

    const int i = reference.AtomCount();
    //for (int i = reference.AtomCount(); i < m_reference.AtomCount(); ++i) {
    Molecule reference_local(reference);
    const int element_local = m_reference.Atom(i).first;
    // std::cout << "Current element to check is " << element_local << " at reference position " << i << std::endl;
    const Position blob = GeometryTools::Centroid(CenterMolecule(reference_local.getGeometry()));
    const Position blob1 = GeometryTools::Centroid(CenterMolecule(target.getGeometry()));

    //std::cout << GeometryTools::Distance(blob, m_reference.Atom(i).second) << std::endl;
    reference_local.addPair(m_reference.Atom(i));
    //double rmsd = 1e22;
    std::map<double, int> match;
    const Geometry ref = CenterMolecule(reference_local.getGeometry());

    bool found_none = true;

    for (int j = 0; j < m_target.AtomCount(); ++j) {
        if (m_target.Atoms()[j] == element_local) {
            if ((m_reference.ConnectedMass(i) != m_target.ConnectedMass(j) || GeometryTools::Distance(blob1, m_target.Atom(j).second) > 1.5 * GeometryTools::Distance(blob, m_reference.Atom(i).second)) && fast) {
                continue;
            }

            found_none = false;
            Molecule target_local(target);
            if (target_local.addPair(m_target.Atom(j))) {
                double rmsd_local = CalculateShortRMSD(ref, target_local);

                if (target_local.AtomCount() < m_target.AtomCount()) {
                    //std::cout << "count smaller " << m_target.AtomCount() << " ";
                    //rmsd = rmsd_local;
                    if (CheckConnections()) {
                        int difference = CheckConnectivitiy(reference_local, target_local);
                        if (difference <= m_pt)
                            match.insert(std::pair<double, int>(rmsd_local, j));
                    } else {
                        match.insert(std::pair<double, int>(rmsd_local, j));
                        m_last_rmsd.push_back(rmsd_local);
                        // std::cout << target_local.AtomCount() <<  " added  (" << match.size() << ")" << std::endl;
                    }
                } else {
                    std::vector<int> inter = intermediate;
                    inter.push_back(j);
                    m_results.insert(std::pair<double, std::vector<int>>(rmsd_local, inter));
                }
            }
        }
    }

    if (match.size() == 0) {
        Molecule ref;
        for (std::size_t index = 0; index < m_reference.AtomCount(); index++) {
            if (i != index)
                ref.addPair(m_reference.Atom(index));
        }
        ref.addPair(m_reference.Atom(i));
        m_reference = ref;
        m_reference_reordered++;
        SolveIntermediate(intermediate);
        return false;
    }

    for (const auto& element : match) {
        std::vector<int> temp = intermediate;
        temp.push_back(element.second);
        m_storage[temp.size() - 1].addItem(temp, element.first);
    }
    return true;
}

void RMSDDriver::ReconstructTarget(const std::vector<int>& atoms)
{
    std::vector<int> reorder_rules;
    Molecule target_cached, target, reference;
    int index = 0;
    for (int atom : atoms) {
        target_cached.addPair(m_target.Atom(atom));
        reorder_rules.push_back(atom);
        reference.addPair(m_reference.Atom(index));
        index++;
    }

    target.LoadMolecule(target_cached);
    for (std::size_t i = 0; i < m_target.AtomCount(); ++i) {
        if (target.Contains(m_target.Atom(i)))
            continue;
        target.addPair(m_target.Atom(i));
        reorder_rules.push_back(i);
    }
    m_reorder_rules = reorder_rules;
    m_target_reordered.LoadMolecule(target);

    if (m_postprocess) {
        double mean_0 = 0;
        bool allow_loop = true;
        //std::cout << CalculateRMSD(reference, target_cached) << std::endl;

        while (allow_loop) {
            Molecule target2, reference2;
            auto terms = IndivRMSD(reference, target_cached);
            double mean = Tools::median(terms);
            if (abs(mean_0 - mean) < 0.1) {
                allow_loop = false;
            }
            for (std::size_t index = 0; index < terms.size(); ++index) {
                if (terms[index] < 1.5 * mean) {
                    target2.addPair(target_cached.Atom(index));
                    reference2.addPair(reference.Atom(index));
                } // else
                //std::cout << index << " " << terms[index] << " " << mean << std::endl;
            }

            target_cached = target2;
            reference = reference2;
            //std::cout << CalculateRMSD(reference, target_cached) << " " << reference.AtomCount() << std::endl;
            mean_0 = mean;
        }
    }

    //std::cout << CalculateRMSD(reference, target_cached) << std::endl;
    Eigen::Matrix3d R = BestFitRotation(reference, target_cached);

    Geometry reference_geom, target_geom;

    reference_geom = GeometryTools::TranslateGeometry(m_reference.getGeometry(), GeometryTools::Centroid(reference.getGeometry()), Position{ 0, 0, 0 });
    target_geom = GeometryTools::TranslateGeometry(target.getGeometry(), GeometryTools::Centroid(target_cached.getGeometry()), Position{ 0, 0, 0 });

    Eigen::MatrixXd tar = target_geom.transpose();

    Geometry rotated = tar.transpose() * R;
    target.setGeometry(rotated);

    m_reference_aligned.LoadMolecule(m_reference);
    m_reference_aligned.setGeometry(reference_geom);
    m_target_aligned.LoadMolecule(target);
}

Molecule RMSDDriver::ApplyOrder(const std::vector<int>& order, const Molecule& mol)
{
    Molecule result;
    for (auto i : order)
        result.addPair(mol.Atom(i));
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
            if (ref_mass[i] == tar_mass[j])
                return std::pair<int, int>(i, j);
        }
    }
    return std::pair<int, int>(-1, -1);
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

    Eigen::Matrix3d R = BestFitRotation(reference_mol, target_mol, 1);

    std::vector<double> terms;

    Geometry reference = CenterMolecule(reference_mol, m_fragment_reference);
    Geometry target = CenterMolecule(target_mol, m_fragment_target);

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

    std::vector<int> new_order(m_reference.AtomCount(), -1), done_ref, done_tar;

    while (done_ref.size() < m_reference.AtomCount()) {
        double distance = 1e10;
        int match_reference = 0;
        int match_target = 0;
        for (int i = 0; i < tar_mol.AtomCount(); ++i) {
            if (std::find(done_tar.begin(), done_tar.end(), i) != done_tar.end())
                continue;
            for (int j = 0; j < ref_mol.AtomCount(); ++j) {
                if (std::find(done_ref.begin(), done_ref.end(), j) != done_ref.end())
                    continue;

                if (tar_mol.Atom(i).first != ref_mol.Atom(j).first)
                    continue;

                const double local_distance = GeometryTools::Distance(tar_mol.Atom(i).second, ref_mol.Atom(j).second);
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
    m_reorder_rules = new_order;
    return true;
}
