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

#pragma once

#include "rmsd/rmsd_functions.h"
#include "rmsd/rmsd_strategies.h" // For AlignmentMethod enum

#include "src/core/molecule.h"
#include "src/core/global.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <queue>

#include <LBFGS.h>

#include "json.hpp"
using json = nlohmann::json;

#include "curcumamethod.h"
#include "src/core/parameter_macros.h"  // Claude Generated - For PARAM macro definitions
#include "src/core/config_manager.h"    // Claude Generated - Modern parameter access layer

// Claude Generated - Forward declarations for strategy pattern
class AlignmentStrategy;
struct AlignmentResult;

// Claude Generated - AlignmentConfig needs full definition for member variable
struct AlignmentConfig {
    AlignmentMethod method = AlignmentMethod::INCREMENTAL; // Method identifier
    int limit = 10; // Limit for template methods
    int element = 7; // Element for template methods
    std::vector<int> element_templates; // Multiple elements for templates
    bool force_reorder = false; // Force reordering
    bool update_rotation = false; // Update rotation optimization
    int threads = 1; // Number of threads
    std::string molalign_bin = "molalign"; // MolAlign binary path
    std::string molalign_args = " -remap -fast -tol 10"; // MolAlign arguments
    int molalign_tolerance = 10; // MolAlign tolerance

    AlignmentConfig() = default;
};

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
    inline int Calculations() const { return m_calculations; }

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
    int m_calculations = 0;
    std::function<double(const Molecule&)> m_evaluator;
};

// Claude Generated 2025: RMSDJson removed - replaced by ParameterRegistry + ConfigManager
// Legacy static JSON object removed - all defaults now managed through Parameter Registry System

class RMSDDriver : public CurcumaMethod {
    // Claude Generated - Friend classes for Strategy Pattern access to private members
    friend class IncrementalAlignmentStrategy;
    friend class TemplateAlignmentStrategy;
    friend class HeavyTemplateStrategy;
    friend class AtomTemplateStrategy;
    friend class InertiaAlignmentStrategy;
    friend class MolAlignStrategy;
    friend class DistanceTemplateStrategy;
    friend class PredefinedOrderStrategy;

public:
    RMSDDriver(const json& controller = json(), bool silent = true);  // Claude Generated 2025: Default to empty JSON, ParameterRegistry provides defaults

    virtual ~RMSDDriver();

    // Claude Generated - Handle move semantics due to std::unique_ptr member
    RMSDDriver(const RMSDDriver&) = delete;
    RMSDDriver& operator=(const RMSDDriver&) = delete;
    RMSDDriver(RMSDDriver&&) = default;
    RMSDDriver& operator=(RMSDDriver&&) = default;

    inline void setReference(const Molecule& reference) { m_reference = reference; }
    inline void setTarget(const Molecule& target)
    {
        m_target = target;
        m_target_original = target;
    }

    void setMatchingAtoms(const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms);

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

    int getKuhnMunkresIterations() const { return m_kuhn_munkres_iterations; }

    void clear();
    void reset();

    void setThreads(int threads) { m_threads = threads; }

    bool MolAlignLib();
    static std::pair<double, Matrix> MakeCostMatrix(const Geometry& reference, const Geometry& target, const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms, int costmatrix);

    Geometry Gradient() const;

    // Claude Generated - Public helper methods for Strategy Pattern access
    Geometry CenterMolecule(const Molecule& mol, int fragment) const;
    Geometry CenterMolecule(const Geometry& geom) const;
    std::pair<Molecule, LimitedStorage> InitialisePair();
    std::vector<int> FillMissing(const Molecule& molecule, const std::vector<int>& order);
    void InsertRotation(std::pair<double, Matrix>& rotation);
    std::pair<int, int> CheckFragments();
    std::pair<Matrix, Position> GetOperateVectors(int fragment_reference, int fragment_target);
    std::pair<Matrix, Position> GetOperateVectors(const std::vector<int>& reference_atoms, const std::vector<int>& target_atoms);
    std::pair<double, Matrix> MakeCostMatrix(const Geometry& reference, const Geometry& target);
    std::pair<double, Matrix> MakeCostMatrix(const std::vector<int>& reference, const std::vector<int>& target);
    std::pair<double, Matrix> MakeCostMatrix(const Matrix& rotation);
    std::pair<std::vector<int>, std::vector<int>> PrepareHeavyTemplate();
    std::pair<std::vector<int>, std::vector<int>> PrepareAtomTemplate(const std::vector<int>& templateatom);

private:
    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    // Claude Generated - Thematic parameter loading methods for better organization
    void LoadFragmentAndThreadingParameters();
    void LoadAlignmentMethodParameters();
    void LoadElementTemplateParameters();
    void LoadCostMatrixParameters();
    void LoadAtomSelectionParameters();
    void LoadFileOrderParameters();
    void DisplayConfigurationSummary();

    // Claude Generated - Strategy pattern methods
    void InitializeAlignmentStrategy();
    AlignmentConfig CreateAlignmentConfig() const;

    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("RMSD") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    void CheckTopology();

    Matrix OptimiseRotation(const Eigen::Matrix3d& rotation);

    std::pair<std::vector<int>, std::vector<int>> PrepareDistanceTemplate();

    std::pair<std::vector<int>, std::vector<int>> PrepareAtomTemplate(int templateatom);

    void FinaliseTemplate();

    std::vector<int> DistanceReorder(const Molecule& reference, const Molecule& target, int max = 2);

    // std::vector<int> FillOrder(const Molecule& reference, const Molecule& target, const std::vector<int>& order);
    std::vector<int> Munkress(const Molecule& reference, const Molecule& target);

    std::vector<int> AlignByVectorPair(std::vector<int> first, std::vector<int> second);
    inline std::vector<int> AlignByVectorPair(std::pair<std::vector<int>, std::vector<int>> pair)
    {
        return AlignByVectorPair(pair.first, pair.second);
    }

    static inline double Cost(double distance, double norm, int costmatrix)
    {
        if (costmatrix == 1)
            return distance * distance;
        else if (costmatrix == 2)
            return distance;
        else if (costmatrix == 3)
            return distance + norm;
        else if (costmatrix == 4)
            return distance * distance + norm * norm;
        else if (costmatrix == 5)
            return distance * norm;
        else if (costmatrix == 6)
            return distance * distance * norm * norm;
        else
            return distance * distance;
    }

    void InitialiseOrder();
    /*
    int CheckConnectivitiy(const Molecule& mol1, const Molecule& mol2) const;
    int CheckConnectivitiy(const Molecule& mol1) const;
    */

    std::pair<double, Matrix> MakeCostMatrix(const std::vector<int>& permutation);
    std::pair<double, Matrix> MakeCostMatrix(const std::pair<std::vector<int>, std::vector<int>>& pair);

    std::vector<int> SolveCostMatrix(Matrix& distance);

    std::pair<Matrix, Position> GetOperateVectors(const Molecule& reference, const Molecule& target);

    Molecule m_reference, m_target, m_target_original, m_reference_aligned, m_reference_original, m_target_aligned, m_target_reordered, m_reorder_reference, m_reorder_target, m_reference_centered, m_target_centered;
    Geometry m_reorder_reference_geometry;
    bool m_force_reorder = false, m_protons = true, m_print_intermediate = false, m_silent = false;
    ConfigManager m_config;  // Claude Generated - Modern type-safe parameter access
    std::vector<std::vector<int>> m_intermediate_results;
    std::map<double, std::vector<int>> m_results, m_intermediate_cost_matrices;
    std::vector<double> m_last_rmsd;
    std::vector<int> m_reorder_rules;
    std::vector<std::vector<int>> m_stored_rules, m_intermedia_rules;
    std::vector<double> m_tmp_rmsd;
    double m_rmsd = 0, m_rmsd_raw = 0, m_scaling = 1.5, m_intermedia_storage = 1, m_threshold = 99, m_damping = 0.8, m_km_convergence = 1e-3;
    bool m_check = false;
    bool m_check_connections = false, m_postprocess = true, m_noreorder = false, m_swap = false, m_dynamic_center = false;
    bool m_update_rotation = false, m_split = false, m_nofree = false;
    bool m_kmstat = false;

    int m_hit = 1, m_pt = 0, m_reference_reordered = 0, m_heavy_init = 0, m_init_count = 0, m_initial_fragment = -1, m_method = 1, m_htopo_diff = -1, m_partial_rmsd = -1, m_threads = 1, m_element = 7, m_write = 0, m_topo = 0;
    int m_munkress_cycle = 1;
    int m_molaligntol = 10;
    int m_limit = 10;
    int m_costmatrix = 1;
    int m_kuhn_munkres_max_iterations = 2;
    int m_kuhn_munkres_iterations = 0;
    double m_cost_limit = 0, m_target_rmsd = 0.0;
    mutable int m_fragment = -1, m_fragment_reference = -1, m_fragment_target = -1;
    std::vector<int> m_initial, m_element_templates;
    std::vector<int> m_reference_atoms, m_target_atoms;
    Eigen::Matrix3d m_rotation;
    std::string m_molalign = "molalign", m_molalignarg = " -remap -fast -tol 10";
    std::map<double, Matrix> m_prepared_cost_matrices;

    // Claude Generated - Strategy pattern members
    std::unique_ptr<AlignmentStrategy> m_alignment_strategy;
    AlignmentConfig m_alignment_config;

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    // Claude Generated - Parameter Registry Integration (October 2025)
    BEGIN_PARAMETER_DEFINITION(rmsd)

    // --- General & Threading ---
    PARAM(threads, Int, 1, "Number of threads for parallel execution.", "Performance", {})
    PARAM(protons, Bool, true, "Include protons in the calculation (opposite of 'heavy').", "General", {"heavy"})
    PARAM(force_reorder, Bool, false, "Force reordering even if atom counts match.", "General", {"reorder"})
    PARAM(no_reorder, Bool, false, "Disable all reordering logic.", "General", {"noreorder"})

    // --- Alignment Method ---
    PARAM(method, String, "incr", "Alignment method: hungarian|incr|template|hybrid|subspace|inertia|molalign|dtemplate|predefined.", "Method", {"RMSDmethod", "rmsdmethod"})
    PARAM(limit, Int, 0, "Limit for subspace and dtemplate methods.", "Method", {})
    PARAM(element, String, "7", "Element(s) for template methods (e.g., \"7,8\").", "Method", {"Element"})
    PARAM(order_file, String, "", "Path to a file with a predefined atom order.", "Method", {"order"})

    // --- Cost Matrix & Kuhn-Munkres ---
    PARAM(cost_matrix_type, Int, 0, "Type of cost matrix to use.", "Cost Matrix", {"costmatrix"})
    PARAM(km_cycles, Int, -1, "Number of Kuhn-Munkres cycles (-1 for auto).", "Cost Matrix", {"cycles"})
    PARAM(km_max_iterations, Int, 1000, "Max iterations for the KM algorithm.", "Cost Matrix", {"km_maxiterations"})
    PARAM(km_convergence, Double, 1e-6, "Convergence threshold for the KM algorithm.", "Cost Matrix", {"km_conv"})

    // --- Atom Selection ---
    PARAM(fragment, Int, -1, "Use only a specific fragment index for both molecules.", "Selection", {})
    PARAM(fragment_reference, Int, -1, "Fragment index for the reference molecule.", "Selection", {})
    PARAM(fragment_target, Int, -1, "Fragment index for the target molecule.", "Selection", {})
    PARAM(reference_atoms, String, "", "Specific atom indices for the reference molecule.", "Selection", {})
    PARAM(target_atoms, String, "", "Specific atom indices for the target molecule.", "Selection", {})

    // --- MolAlign Integration ---
    PARAM(molalign_bin, String, "molalign", "Path to the molalign binary.", "MolAlign", {"molalignbin"})
    PARAM(molalign_args, String, "", "Additional arguments for the molalign binary.", "MolAlign", {})
    PARAM(molalign_tolerance, Int, 100, "Tolerance for molalign.", "MolAlign", {"molaligntol"})

    // --- Advanced Options ---
    PARAM(init, Int, -1, "Initial fragment index.", "Advanced", {})
    PARAM(pt, Int, 0, "Number of allowed proton transfers.", "Advanced", {})
    PARAM(storage, Double, 1.0, "Intermediate storage factor.", "Advanced", {})
    PARAM(dynamic_center, Bool, false, "Use dynamic centering.", "Advanced", {"DynamicCenter"})
    PARAM(topo, Int, 0, "Topology comparison mode.", "Advanced", {})
    PARAM(write, Int, 0, "Write intermediate results.", "Advanced", {})
    PARAM(update_rotation, Bool, false, "Update rotation optimization.", "Advanced", {"update-rotation"})
    PARAM(damping, Double, 0.8, "Damping factor for optimization.", "Advanced", {})
    PARAM(split, Bool, false, "Split calculation into fragments.", "Advanced", {})
    PARAM(kmstat, Bool, false, "Output Kuhn-Munkres statistics.", "Advanced", {})
    PARAM(target_rmsd, Double, 0.0, "Target RMSD value for early stopping.", "Advanced", {})
    PARAM(check, Bool, false, "Enable connectivity checking.", "Advanced", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};

using namespace LBFGSpp;

class LBFGSRotation {
public:
    LBFGSRotation(int n_)
    {
    }
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
        Eigen::Matrix3d n;
        n = Eigen::AngleAxisd(x[0], Eigen::Vector3d::UnitX())
            * Eigen::AngleAxisd(x[1], Eigen::Vector3d::UnitY())
            * Eigen::AngleAxisd(x[2], Eigen::Vector3d::UnitZ());

        // Eigen::MatrixXd tar = m_target.transpose();

        // Geometry rotated = tar.transpose() * n;
        double fx = 0.0;
        double dx = 1e-5;

        auto result = RMSDDriver::MakeCostMatrix(m_reference, m_target * n, m_reference_atoms, m_target_atoms, m_costmatrix);
        // std::cout << result.first << std::endl;
        Eigen::VectorXd tmp = x;
        for (int i = 0; i < 3; ++i) {
            tmp[i] += dx;
            // std::cout << tmp.transpose() << std::endl;
            Eigen::Matrix3d n;
            n = Eigen::AngleAxisd(tmp[0], Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(tmp[1], Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(tmp[2], Eigen::Vector3d::UnitZ());

            // Geometry rotated = tar.transpose() * n;

            auto p = RMSDDriver::MakeCostMatrix(m_reference, m_target * n, m_reference_atoms, m_target_atoms, m_costmatrix).first;

            tmp[i] -= 2 * dx;
            // std::cout << tmp.transpose() << std::endl;

            n = Eigen::AngleAxisd(tmp[0], Eigen::Vector3d::UnitX())
                * Eigen::AngleAxisd(tmp[1], Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(tmp[2], Eigen::Vector3d::UnitZ());

            // rotated = tar.transpose() * n;

            auto m = RMSDDriver::MakeCostMatrix(m_reference, m_target * n, m_reference_atoms, m_target_atoms, m_costmatrix).first;
            // std::cout << p << " " << m << " " << (p-m)/(2*dx) << std::endl << std::endl;
            grad[i] = (p - m) / (2 * dx);
        }

        return result.first;
    }

    Vector Parameter() const { return m_parameter; }

    Geometry m_reference, m_target;
    std::vector<int> m_reference_atoms, m_target_atoms;
    int m_costmatrix;

private:
    Vector m_parameter;
};
