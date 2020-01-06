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
#include "src/tools/geometry.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "rmsd.h"

RMSDDriver::RMSDDriver(const Molecule &reference, const Molecule &target)
    : m_reference(reference)
    , m_target(target)
{

}

void RMSDDriver::AutoPilot()
{
    m_rmsd_raw = CalculateRMSD(m_reference, m_target);

    if (m_reference.Atoms() != m_target.Atoms() || ForceReorder()) {
        std::cout << "Molecules are different. What now?\n";

        if (m_reference.AtomCount() == m_target.AtomCount()) {
            std::cout << "Try to reorder the structures?\n";
            std::cout << "Initial RMSD is " << m_rmsd_raw << std::endl;
            ReorderMolecule();
        }
    }

    Molecule *reference = new Molecule, *target = new Molecule;
    m_rmsd = CalculateRMSD(m_reference, m_target_reordered, reference, target);

    m_reference_aligned.LoadMolecule(reference);
    m_target_aligned.LoadMolecule(target);
}

Eigen::Matrix3d RMSDDriver::BestFitRotation(const Molecule& reference_mol, const Molecule& target_mol, int factor) const
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

double RMSDDriver::CalculateRMSD(const Molecule& reference_mol, const Molecule& target_mol, Molecule* ret_ref, Molecule* ret_tar, int factor) const
{
    Eigen::Matrix3d R = BestFitRotation(reference_mol, target_mol, factor);

    double rmsd = 0;

    Geometry reference = CenterMolecule(reference_mol);
    Geometry target = CenterMolecule(target_mol);
    Eigen::MatrixXd tar = target.transpose();

    // std::cout << std::endl << R << std::endl << std::endl;

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
        rotated.transposeInPlace();
        ret_tar->setGeometry(rotated);
    }

    if (ret_ref != nullptr) {
        ret_ref->LoadMolecule(reference_mol);
        reference.transposeInPlace();
        ret_ref->setGeometry(reference);
    }
    return sqrt(rmsd);
}
std::vector<double> RMSDDriver::IndivRMSD(const Molecule& reference_mol, const Molecule& target_mol, int factor) const
{
    Eigen::Matrix3d R = BestFitRotation(reference_mol, target_mol, factor);

    std::vector<double> terms;

    Geometry reference = CenterMolecule(reference_mol);
    Geometry target = CenterMolecule(target_mol);
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

Geometry RMSDDriver::CenterMolecule(const Molecule &mol) const
{
    return GeometryTools::TranslateMolecule(mol, mol.Centroid(), Position{0, 0, 0});
}

void RMSDDriver::ReorderMolecule()
{
    int inter_size = m_reference.AtomCount() * (m_reference.AtomCount() - 1);
    m_storage = std::vector<IntermediateStorage>(m_reference.AtomCount() - 1, IntermediateStorage(inter_size));

    /* Lets initialise a molecule with two atom from the reference */
    Molecule reference;
    int index = 0;
    std::vector<int> elements = m_reference.Atoms();
    for (int i = 0; i < m_reference.AtomCount() && index < 2; i++) {
        if (elements[i] != 1) {
            reference.addPair(m_reference.Atom(i));
            index++;
        }
    }

    std::vector<int> elements_target = m_target.Atoms();
    std::vector<int> tmp_reference = reference.Atoms();
    for (int i = 0; i < m_target.AtomCount(); ++i)
        for (int j = i + 1; j < m_target.AtomCount(); ++j) {
            if (tmp_reference[0] == elements_target[i] && tmp_reference[1] == elements_target[j])
                m_intermediate_results.push({ i, j });
            if (tmp_reference[0] == elements_target[j] && tmp_reference[1] == elements_target[i])
                m_intermediate_results.push({ j, i });
        }

    m_rmsd = CalculateRMSD(m_reference, m_target) / double(m_reference.AtomCount());

    while (m_intermediate_results.size()) {
        std::vector<int> inter = m_intermediate_results.front();
        //std::cout << inter.size() << std::endl;
        //if(inter.size() <  m_target.AtomCount())
        SolveIntermediate(inter);
        m_intermediate_results.pop();
    }

    for (int i = 0; i < m_storage.size(); ++i) {
        std::cout << double(i) / double(m_reference.AtomCount()) * 100 << " % done " << std::endl;

        for (const auto& element : (*m_storage[i].data())) {
            SolveIntermediate(element.second);
        }
    }

    int count = 0;
    std::vector<int> intermediate;
    for (const auto& element : m_results) {
        for (int i = 0; i < element.second.size(); i++) {
            m_target_reordered.addPair(m_target.Atom(element.second[i]));
        }
        m_rmsd = CalculateRMSD(m_reference, m_target_reordered);

        count++;
        if (count == 1)
            break;
    }

    return;
}

void RMSDDriver::SolveIntermediate(std::vector<int> intermediate)
{
    Molecule reference;
    Molecule target;
    std::vector<int> elements_target = m_target.Atoms();

    for (int i = 0; i < intermediate.size(); i++) {
        reference.addPair(m_reference.Atom(i));
        target.addPair(m_target.Atom(intermediate[i]));
    }

    int i = reference.AtomCount();
    //for (int i = reference.AtomCount(); i < m_reference.AtomCount(); ++i) {
    Molecule reference_local(reference);
    int element_local = m_reference.Atom(i).first;
    reference_local.addPair(m_reference.Atom(i));
    double rmsd = 1e22;
    std::map<double, int> match;

    for (int j = 0; j < m_target.AtomCount(); ++j) {
        if (elements_target[j] == element_local) {
            Molecule target_local(target);
            if (target_local.addPair(m_target.Atom(j))) {
                double rmsd_local = CalculateRMSD(reference_local, target_local);
                if (target_local.AtomCount() < m_target.AtomCount()) {
                    rmsd = rmsd_local;
                    match.insert(std::pair<double, int>(rmsd_local, j));
                } else {
                    std::vector<int> inter = intermediate;
                    inter.push_back(j);
                    m_results.insert(std::pair<double, std::vector<int>>(rmsd_local, inter));
                }
            }
        }
    }
    for (const auto& element : match) {
        std::vector<int> temp = intermediate;
        temp.push_back(element.second);
        m_storage[temp.size() - 1].addItem(temp, element.first);
    }
}
