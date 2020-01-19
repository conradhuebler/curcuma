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

#pragma once

#include "src/core/molecule.h"
#include "src/core/global.h"

#include <map>
#include <queue>

class IntermediateStorage {
public:
    inline IntermediateStorage(unsigned int size)
        : m_size(size)
    {
    }

    inline void addItem(const std::vector<int>& vector, double rmsd)
    {

        m_shelf.insert(std::pair<double, std::vector<int>>(rmsd, vector));
        if (m_shelf.size() >= m_size)
            m_shelf.erase(--m_shelf.end());
    }

    const std::map<double, std::vector<int>>* data() const { return &m_shelf; }

private:
    unsigned int m_size;
    std::map<double, std::vector<int>> m_shelf;
};

class RMSDDriver{

public:
    RMSDDriver(const Molecule& reference, const Molecule& target);
    RMSDDriver(const Molecule* reference, const Molecule* target);

    /*! \brief Use the AutoPilot to automatically perform everything, results are stored as long the object exsist */
    void AutoPilot();

    double CalculateRMSD();
    double CalculateRMSD(const Molecule& reference, const Molecule& target, Molecule* ret_ref = nullptr, Molecule* ret_tar = nullptr, int factor = 1) const;

    void ProtonDepleted();

    std::vector<double> IndivRMSD(const Molecule& reference, const Molecule& target, int factor = 1) const;

    void ReorderMolecule();

    /*! \brief Return the reference molecule centered */
    inline Molecule ReferenceAligned() const { return m_reference_aligned; }

    /*! \brief Return the target molecule centered and aligned to the reference molecule */
    inline Molecule TargetAligned() const { return m_target_aligned; }

    /*! \brief Return the target molecule reorderd but remaining at the original position */
    inline Molecule TargetReorderd() const { return m_target_reordered; }

    /*! \brief Return best-fit reordered RMSD */
    inline double RMSD() const { return m_rmsd; }

    /*! \brief Return best-fit RMSD with reordering */
    inline double RMSDRaw() const { return m_rmsd_raw; }

    /*! \brief Force Reordering, even the sequence of elements are equal */
    inline void setForceReorder(bool reorder) { m_force_reorder = reorder; }

    /*! \brief Check, if Reordering is forced */
    inline bool ForceReorder() const { return m_force_reorder; }

    /*! \brief Get n'th/rd best fit result */
    Molecule getFitIndex(int index);

    /*! \brief Set the index of the fragment that is used for rmsd calculation/atom reordering */
    inline void setFragment(int fragment) { m_fragment = fragment; }

    /*! \brief Set to use protons (true = default or false) */
    inline void setProtons(bool protons) { m_protons = protons; }

    /*! \brief Set Connectivitiy Check forced (true or false = default) */
    inline void setCheckConnections(bool check) { m_check_connections = check; }

    /*! \brief Force Connectivitiy Check */
    inline bool CheckConnections() const { return m_check_connections; }

    /*! \brief Number of Proton changes allowed */
    inline int ProtonTransfer() const { return m_pt; }

    /*! \brief Set number of allowed proton transfer */
    inline void setProtonTransfer(int pt) { m_pt = pt; }

    /*! \brief Set silent */
    inline void setSilent(bool silent) { m_silent = silent; }

private:
    void ReorderStraight();

    void InitialisePair();

    void SolveIntermediate(std::vector<int> intermediate);

    bool CheckConnectivitiy(const Molecule& mol1, const Molecule& mol2) const;
    bool CheckConnectivitiy(const Molecule& mol1) const;

    Eigen::Matrix3d BestFitRotation(const Molecule& reference, const Molecule& target, int factor = 1) const;
    Eigen::Matrix3d BestFitRotation(const Geometry& reference, const Geometry& target, int factor = 1) const;

    Geometry CenterMolecule(const Molecule &mol) const;
    Geometry CenterMolecule(const Geometry& mol) const;

    Molecule m_reference, m_target, m_reference_aligned, m_target_aligned, m_target_reordered;
    bool m_force_reorder = false, m_protons = true, m_print_intermediate = false, m_silent = false;
    std::queue<std::vector<int>> m_intermediate_results;
    std::map<double, std::vector<int>> m_results;
    std::map<int, std::vector<int>> m_connectivity;
    std::vector<IntermediateStorage> m_storage;
    double m_rmsd = 0, m_rmsd_raw = 0, m_scaling = 1.5;
    bool m_check_connections = false;
    int m_hit = 1, m_pt = 0;
    mutable int m_fragment = -1;
};
