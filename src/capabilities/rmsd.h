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

#include "rmsd_functions.h"

#include "src/core/molecule.h"
#include "src/core/global.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <chrono>
#include <functional>
#include <map>
#include <queue>

#include "json.hpp"
using json = nlohmann::json;

#include "curcumamethod.h"

struct StructComp {
    double rmsd = 0;
    double diff_hydrogen_bonds = 0;
    double diff_topology = 0;
};

class RMSDThread : public CxxThread {
public:
    RMSDThread(const Molecule& reference_molecule, const Molecule& target, const Geometry& reference, const Matrix& reference_topology, const std::vector<int> intermediate, double connected_mass, int element, int topo);
    inline virtual ~RMSDThread() = default;

    int execute() override;

    const std::map<double, std::vector<int>>* data() const { return &m_shelf; }
    inline int Match() const { return m_match; }

private:
    Molecule m_target;
    Molecule m_reference_molecule;
    Geometry m_reference;
    Matrix m_reference_topology;
    std::map<double, std::vector<int>> m_shelf;
    std::vector<int> m_intermediate;
    double m_connected_mass = 0;
    int m_element = -1;
    int m_match;
    int m_topo = 0;
    std::function<double(const Molecule&)> m_evaluator;
};

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

static const json RMSDJson = {
    { "reorder", false },
    { "check", false },
    { "heavy", false },
    { "fragment", -1 },
    { "fragment_reference", -1 },
    { "fragment_target", -1 },
    { "init", -1 },
    { "pt", 0 },
    { "silent", false },
    { "storage", 1.0 },
    { "method", "incr" },
    { "noreorder", false },
    { "threads", 1 },
    { "Element", 7 },
    { "DynamicCenter", false },
    { "order", "" },
    { "check", false },
    { "topo", 0 },
    { "write", 0 },
    { "moi", false },
    { "update-rotation", false },
    { "damping", 0.8 },
    { "split", false }
};

class RMSDDriver : public CurcumaMethod {
public:
    RMSDDriver(const json& controller = RMSDJson, bool silent = true);

    virtual ~RMSDDriver();

    inline void setReference(const Molecule& reference) { m_reference = reference; }
    inline void setTarget(const Molecule& target)
    {
        m_target = target;
        m_target_original = target;
    }

    double Rules2RMSD(const std::vector<int> rules, int fragment = -1);
    StructComp Rule2RMSD(const std::vector<int> rules, int fragment = 1);

    double CalculateRMSD();
    double CalculateRMSD(const Molecule& reference, const Molecule& target, Molecule* ret_ref = nullptr, Molecule* ret_tar = nullptr, int factor = 1) const;

    void ProtonDepleted();

    std::vector<double> IndivRMSD(const Molecule& reference, const Molecule& target, int factor = 1) const;

    void ReorderMolecule();

    /*! \brief Return the reference molecule centered */
    inline Molecule ReferenceAligned() const { return m_reference_aligned; }

    /*! \brief Return the reference molecule centered */
    inline const Molecule* ReferenceAlignedReference() const { return &m_reference_aligned; }

    /*! \brief Return the target molecule centered and aligned to the reference molecule */
    inline Molecule TargetAligned() const { return m_target_aligned; }

    /*! \brief Return the target molecule centered and aligned to the reference molecule */
    inline const Molecule* TargetAlignedReference() const { return &m_target_aligned; }

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
    inline void setFragment(int fragment)
    {
        m_fragment = fragment;
        m_fragment_reference = fragment;
        m_fragment_target = fragment;
    }

    /*! \brief Set the index of the fragment that is used for rmsd calculation/atom reordering */
    inline void setFragmentTarget(int fragment) { m_fragment_target = fragment; }

    /*! \brief Set the index of the fragment that is used for rmsd calculation/atom reordering */
    inline void setFragmentReference(int fragment) { m_fragment_reference = fragment; }

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

    /*! \brief Set silent */
    inline void setPartialRMSD(bool partial_rmsd) { m_partial_rmsd = partial_rmsd; }

    void setScaling(double scaling) { m_scaling = scaling; }

    inline void setIntermediateStorage(double storage) { m_intermedia_storage = storage; }

    inline std::vector<int> ReorderRules() const { return m_reorder_rules; }

    inline void setInitial(std::vector<int> initial) { m_initial = initial; }
    inline void setInitialFragment(int fragment) { m_initial_fragment = fragment; }

    void start() override;

    Molecule ApplyOrder(const std::vector<int>& order, const Molecule& mol);

    std::vector<std::vector<int>> StoredRules() const { return m_stored_rules; }

    inline int HBondTopoDifference() const { return m_htopo_diff; }

    double SimpleRMSD();

    double BestFitRMSD();

    double CustomRotation();

    double PartialRMSD(const Molecule& ref, const Molecule& tar);

    void clear();

    void setThreads(int threads) { m_threads = threads; }

private:
    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("RMSD") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    void ReorderIncremental();

    void HeavyTemplate();

    void AtomTemplate();

    void TemplateFree();

    void CheckTopology();

    std::pair<std::vector<int>, std::vector<int>> PrepareHeavyTemplate();

    std::pair<std::vector<int>, std::vector<int>> PrepareAtomTemplate(int templateatom);
    std::pair<std::vector<int>, std::vector<int>> PrepareAtomTemplate(const std::vector<int>& templateatom);

    void FinaliseTemplate(std::pair<std::vector<int>, std::vector<int>> initial);

    std::vector<int> DistanceReorder(const Molecule& reference, const Molecule& target);

    std::vector<int> DistanceReorderV1(const Molecule& reference, const Molecule& target);
    std::vector<int> DistanceReorderV2(const Molecule& reference, const Molecule& target);
    std::pair<std::vector<int>, std::vector<int>> DistanceReorderV3(const Molecule& reference, const Molecule& target);
    std::vector<int> FillOrder(const Molecule& reference, const Molecule& target, const std::vector<int>& order);

    std::vector<int> AlignByVectorPair(std::vector<int> first, std::vector<int> second);
    inline std::vector<int> AlignByVectorPair(std::pair<std::vector<int>, std::vector<int>> pair)
    {
        return AlignByVectorPair(pair.first, pair.second);
    }

    std::vector<int> FillMissing(const Molecule& molecule, const std::vector<int>& order);

    void InitialiseOrder();
    std::pair<Molecule, LimitedStorage> InitialisePair();
    /*
    int CheckConnectivitiy(const Molecule& mol1, const Molecule& mol2) const;
    int CheckConnectivitiy(const Molecule& mol1) const;
*/
    bool TemplateReorder();
    std::pair<int, int> CheckFragments();

    Geometry CenterMolecule(const Molecule& mol, int fragment) const;
    Geometry CenterMolecule(const Geometry& molt) const;

    std::pair<Matrix, Position> GetOperateVectors(int fragment_reference, int fragment_target);
    std::pair<Matrix, Position> GetOperateVectors(const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms);
    std::pair<Matrix, Position> GetOperateVectors(const Molecule& reference, const Molecule& target);

    Molecule m_reference, m_target, m_target_original, m_reference_aligned, m_target_aligned, m_target_reordered, m_reorder_reference, m_reorder_target;
    Geometry m_reorder_reference_geometry;
    bool m_force_reorder = false, m_protons = true, m_print_intermediate = false, m_silent = false, m_moi = false;
    std::vector<std::vector<int>> m_intermediate_results;
    std::map<double, std::vector<int>> m_results;
    std::vector<double> m_last_rmsd;
    std::vector<int> m_reorder_rules;
    std::vector<std::vector<int>> m_stored_rules;
    std::map<int, std::vector<int>> m_connectivity;
    double m_rmsd = 0, m_rmsd_raw = 0, m_scaling = 1.5, m_intermedia_storage = 1, m_threshold = 99, m_damping = 0.8;
    bool m_check = false;
    bool m_check_connections = false, m_postprocess = true, m_noreorder = false, m_swap = false, m_dynamic_center = false;
    bool m_update_rotation = false, m_split = false;
    int m_hit = 1, m_pt = 0, m_reference_reordered = 0, m_heavy_init = 0, m_init_count = 0, m_initial_fragment = -1, m_method = 1, m_htopo_diff = -1, m_partial_rmsd = -1, m_threads = 1, m_element = 7, m_write = 0, m_topo = 0;
    mutable int m_fragment = -1, m_fragment_reference = -1, m_fragment_target = -1;
    std::vector<int> m_initial, m_element_templates;
};
