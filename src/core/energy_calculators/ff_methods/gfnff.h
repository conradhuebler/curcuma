/*
 * <GFN-FF Implementation for Curcuma>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * GFN-FF (Geometry, Frequency, Noncovalent - Force Field) is a fully
 * automated, quantum chemistry-based force field for the accurate
 * description of structures and dynamics of large molecular systems.
 */

#pragma once

#include "src/core/units.h"
#include "json.hpp"
#include "src/core/config_manager.h"
#include "src/core/parameter_macros.h"
#include "src/core/energy_calculators/ff_methods/forcefield.h"

#include "src/core/energy_calculators/ff_methods/ff_workspace.h"  // Claude Generated (Mar 2026): Unified workspace
#include "src/core/energy_calculators/ff_methods/eeq_solver.h"  // EEQ charge calculation (Dec 2025 - Phase 3)
#include "src/core/energy_calculators/ff_methods/huckel_solver.h"  // Full Hückel calculation (Jan 2026 - Phase 1)
#include "src/core/energy_calculators/dispersion/d4param_generator.h"  // Claude Generated (Feb 15, 2026): D4 for dc6dcn gradient
#include "src/core/energy_calculators/ff_methods/alpb_solvation.h"  // Claude Generated (Mar 2026): ALPB solvation
#include "src/core/global.h"
#include "src/core/functional_groups.h"
#include "src/core/periodic_table.h"
#include <functional>
#include <algorithm>
#include <tuple>
#include <utility>
#include <limits>
#include <optional>
#include <vector>
#include <memory>

using json = nlohmann::json;

/**
 * @brief GFN-FF Implementation as Standalone Force Field
 *
 * GFN-FF combines quantum chemical accuracy with force field efficiency.
 * It provides:
 * - Automatic parametrization based on extended tight-binding (xTB)
 * - Accurate treatment of non-covalent interactions
 * - Gradients for geometry optimization and dynamics
 * - Coverage of the periodic table up to Z=86
 *
 * References:
 * - Spicher, S.; Grimme, S. "Robust Atomistic Modeling of Materials,
 *   Organometallic, and Biochemical Systems" Angew. Chem. Int. Ed. 59, 15665 (2020)
 *
 * Note: GFNFF is semantically a force field, not a quantum method.
 * It was previously inheriting from QMInterface but has been refactored
 * to standalone class for semantic correctness (Phase 2, November 2025).
 */

/**
 * @brief Triangular indexing function for symmetric matrices
 *
 * Claude Generated (January 10, 2026) - Port from Fortran gfnff_helpers.f90:416-423
 * Converts atom pair (i,j) to linear index for upper-triangular storage
 *
 * @param i First atom index (0-based)
 * @param j Second atom index (0-based)
 * @return Linear index in triangular array
 *
 * Formula (C++ 0-based): lin = min(i,j) + max(i,j)*(max(i,j)+1)/2
 */
inline int lin(int i, int j) {
    int imax = std::max(i, j);
    int imin = std::min(i, j);
    return imin + imax * (imax + 1) / 2;
}

/**
 * @brief Sparse N x N integer table for short-range topological pair data
 *
 * Claude Generated (Sep 2026). Holds the topological pair tables of GFN-FF (the
 * reference's topo%bpair and curcuma's BFS bond-count distance) which are only
 * informative for pairs a few bonds apart: every other off-diagonal entry has one
 * and the same "far" value. Storing only the near pairs makes both tables
 * O(N * k) instead of O(N^2) (two dense tables were ~1.7 GB at N = 14640).
 *
 * Row i holds (j, value) sorted by ascending j, so iterating a row visits the
 * partners in the same order a dense `for (j ...)` loop would.
 */
struct SparseTopoTable {
    int n = 0;
    int diag_value = 0;       ///< value for i == j
    int default_value = 0;    ///< value for every pair not stored in a row
    std::vector<std::vector<std::pair<int, int>>> rows;

    void reset(int size, int diag, int far)
    {
        n = size;
        diag_value = diag;
        default_value = far;
        rows.assign(size, {});
    }

    bool empty() const { return n == 0; }

    /// Value of pair (i, j); binary search in row i (rows hold a few dozen entries).
    int get(int i, int j) const
    {
        if (i == j) return diag_value;
        const auto& r = rows[i];
        auto it = std::lower_bound(r.begin(), r.end(), j,
            [](const std::pair<int, int>& e, int key) { return e.first < key; });
        return (it != r.end() && it->first == j) ? it->second : default_value;
    }

    const std::vector<std::pair<int, int>>& row(int i) const { return rows[i]; }

    /// Number of stored (off-diagonal, non-default) entries, both directions counted.
    size_t storedEntries() const
    {
        size_t s = 0;
        for (const auto& r : rows) s += r.size();
        return s;
    }
};

/**
 * @brief Check if an atom is classified as a metal
 *
 * Claude Generated (January 2026) - Phase 4: Metal scaling
 * Reference: gfnff_method.cpp:4382-4386
 *
 * @param atomic_number Atomic number (Z)
 * @return true if atom is a transition metal or lanthanide/actinide
 *
 * Metal classification: Transition metals (Sc-Zn, Y-Cd, La-Hg) + Lanthanides/Actinides
 */
inline bool isMetalAtom(int atomic_number) {
    return (atomic_number >= 21 && atomic_number <= 30) ||   // Sc-Zn
           (atomic_number >= 39 && atomic_number <= 48) ||   // Y-Cd
           (atomic_number >= 57 && atomic_number <= 80) ||   // La-Hg
           (atomic_number >= 89 && atomic_number <= 103);    // Ac-Lr
}

/**
 * @brief Unified GFN-FF energy + timing report (CPU and GPU paths share this format)
 *
 * Claude Generated (May 2026): Single source of truth for the verbosity-2 output that
 * was previously fragmented between forcefield.cpp (CPU) and gfnff_gpu_method.cpp (GPU).
 * Populated by the wrapper after the energy calculation, then printed by
 * printGFNFFEnergyReport() for a uniform user-facing display.
 *
 * Conventions:
 *   - All times in milliseconds. -1.0 means "not applicable for this path".
 *   - cpu_sum: sum of per-thread CPU time (NOT wall-clock). Wall-clock is in t_pool_wall.
 *   - gpu: GPU phase wall-clock (proxy for kernel time; kernels are async, no per-kernel CUDA events).
 */
struct GFNFFEnergyReport {
    // Energy components (Hartree)
    double bond = 0, angle = 0, dihedral = 0, inversion = 0, stors = 0;
    double dispersion = 0, bonded_rep = 0, nonbonded_rep = 0;
    double coulomb = 0, hbond = 0, xbond = 0, atm = 0, batm = 0;
    // Claude Generated (May 2026, HB-investigation): per-case HB diagnostic split.
    // hbond = case1 + case2 + case3 + case4. Counts compare against Fortran nhb1/nhb2.
    double hbond_case1 = 0, hbond_case2 = 0, hbond_case3 = 0, hbond_case4 = 0;
    int hbond_case1_count = 0, hbond_case2_count = 0, hbond_case3_count = 0, hbond_case4_count = 0;
    double total = 0;

    // Per-term timing: cpu_sum (across all threads), gpu (phase wall-clock); -1 = N/A
    struct TermTiming { double cpu_sum = -1.0, gpu = -1.0; };
    TermTiming t_bond, t_angle, t_dihedral, t_inversion, t_stors;
    TermTiming t_dispersion, t_bonded_rep, t_nonbonded_rep;
    TermTiming t_coulomb, t_hbond, t_xbond, t_atm, t_batm;

    // Gradient
    double gradient_norm = -1.0;
    TermTiming t_gradient;

    // Phase summary — wall-clock per phase
    double t_wall = 0.0;          // total wall-clock of Calculation()
    int    n_cpu_threads = 1;
    double t_pool_wall = -1.0;    // wall-clock of thread pool (ForceField)
    double t_pool_cpu_sum = -1.0; // sum of all per-term CPU timings
    double t_cn_eeq_cpu = -1.0;   // CN + EEQ wall-clock (CPU serial)
    double t_hbxb = -1.0;         // HB/XB re-detection wall-clock
    double t_gradient_cpu = -1.0; // chain-rule gradient (CPU serial after pool)

    // GPU-only phase timings (-1 on CPU path)
    double t_gpu_cn = -1.0;       // GPU CN kernel phase wall-clock
    double t_cpu_eeq_gpu_path = -1.0; // CPU EEQ overlapping with GPU charge-indep kernels
    double t_gpu_phase2 = -1.0;   // Phase 2: Coulomb + DMA

    // One-time init costs (carried forward from param-gen phase so every
    // energy report shows the full story: setup + calculation)
    double t_topology = -1.0;  // ms, -1 if not measured / cached
    double t_param_gen = -1.0; // ms, -1 if not measured / cached
    double t_gpu_upload = -1.0; // ms, -1 on CPU path or if not measured

    bool is_gpu = false;
};

/**
 * @brief Print GFN-FF energy decomposition + timing breakdown at verbosity >= 2.
 *
 * Claude Generated (May 2026): Same format for CPU and GPU paths. Uses CurcumaLogger::result().
 */
void printGFNFFEnergyReport(const GFNFFEnergyReport& r);

/**
 * @brief One-time parameter generation profiling report (printed at verbosity >= 2)
 *
 * Claude Generated (May 2026): Captures timings for both topology construction
 * (`calculateTopologyInfo()`) and parameter generation (`generateGFNFFParameterSet()`).
 * On topology-cache hit, topology phase fields are 0 (skipped) but the table still
 * prints so the layout is stable.
 *
 * Future Phase B will fill `t_parallel_block_*` for OpenMP / CxxThreadPool comparison.
 */
struct GFNFFParamGenReport {
    enum Backend { Sequential, CxxThreadPool, OpenMPSections };
    Backend backend = Sequential;
    int n_threads = 1;
    int n_atoms = 0;
    bool topology_cached = false;     // true if .topo.json hit (most topo phases skipped)

    // Topology sub-phases (-1 = not measured this run)
    double t_distance_matrix    = -1;
    double t_cn_hyb_pi_rings    = -1;
    double t_eeq_phase1         = -1;
    double t_eeq_phase1_corr    = -1;
    double t_eeq_phase2         = -1;
    // A0 (Jul 2026): the old single `t_pi_bond_orders` bucket bracketed three
    // unrelated pieces of work, which made it useless for guiding optimisation.
    // Split into the ipis per-pi-system EEQ re-solves, the FT-HMO solve itself,
    // and the bond-type classification loop.
    double t_pi_charges_eeq     = -1;
    double t_huckel             = -1;
    double t_bond_types         = -1;
    double t_topo_distances     = -1;
    double t_topology_total     = -1;
    int    n_pi_systems         = -1;  // number of pi-systems (drives the two above)

    // Parameter generation phases — all measured in generateGFNFFParameterSet
    double t_bonds       = -1;
    double t_angles      = -1;
    double t_torsions    = -1;
    double t_inversions  = -1;
    double t_storsions   = -1;
    double t_coulomb     = -1;
    double t_repulsion   = -1;
    double t_dispersion  = -1;
    double t_batm        = -1;
    double t_hbxb        = -1;
    double t_crossref    = -1;

    // Parallel-block summary (cxxthreadpool path; -1 if sequential or backend=Sequential)
    double t_parallel_block_wall    = -1;
    double t_parallel_block_cpu_sum = -1;

    double t_param_gen_total = -1;
};

void printGFNFFParamGenReport(const GFNFFParamGenReport& r);

// P2b (Apr 2026): CN cutoff parameters — configurable via CLI
// Three modes:
//   cn_cutoff_bohr > 0: Neighbor-list mode (default 10.0 Bohr, fast O(N*k))
//   cn_cutoff_bohr = 0, cn_accuracy > 0: Fortran accuracy-based threshold (cnthr = 100 - log10(acc)*50)
//   cn_cutoff_bohr = 0, cn_accuracy = 0: Full O(N²) reference mode (no cutoff)
// FIX (Jul 23, 2026): default raised 6.0 -> 10.0 Bohr. The reference cnthr
// (gfnff_param.f90:551, accuracy=1) is 100 Bohr^2 = 10 Bohr; the old 6.0 Bohr was
// TIGHTER than the reference and truncated erf-CN contributions for heavy/metal atoms
// (large covalent radii push the erf transition past 6 Bohr). That gave a wrong
// dynamic-bond-r0 CN in the FFWorkspace energy path — e.g. PR23 (Ir complex) bond
// energy +1.63 kcal. At 10 Bohr the CN is converged and matches the reference; the SP
// bond energy is bit-identical to the legacy per-bond path (which used the full CN).
BEGIN_PARAMETER_DEFINITION(gfnff)
// ---- EEQ charges ----
// The PARAMs of the former registry module eeq_solver live here since Oct 2026 (EEQSolver is
// used by GFN-FF only). EEQSolver still reads the JSON scope "eeq_solver": -eeq_solver.<name>
// lands there and wins over -gfnff.<name>; the flat -<name> and -gfnff.<name> reach it through
// the gfnff scope (GFNFF::forwardEEQSolverParams, ConfigManager flat keys, fixEEQSolverChoice in
// method_factory.cpp, eeq_value() in gfnff_gpu_method_impl.h). ParameterRegistry::canonicalModule
// maps lookups by "eeq_solver" to this module.
PARAM(accuracy, String, "normal", "Accepted for compatibility, but no code reads it, so it has no effect. The tolerances are set by eeq_ppcg_tol, pcg_tolerance, cn_accuracy and hb_accuracy.", "EEQ", {},
        "enum=loose|normal|medium|high")
PARAM(allow_unconverged_charges, Bool, false, "A failed EEQ solve always falls back to Phase-1 or uniform charges; true reports this as a warning instead of an error. The GPU path ignores it.", "EEQ", {})
PARAM(skip_phase2, Bool, false, "Use the Phase-1 topology charges as final charges and skip the geometry-dependent Phase-2 EEQ solve. Faster, less accurate.", "EEQ", {})
PARAM(solve_method, String, "auto", "EEQ linear solver: auto, cholesky, ppcg, pcg, lu, ldlt or batched (inexact). auto takes ppcg from eeq_ppcg_min_atoms and eeq_ppcg_min_nfrag on, else cholesky. See EEQ_SOLVER_CONTROL.md.", "EEQ", {})
PARAM(solve, String, "auto", "Older name of solve_method; a non-auto solve_method wins when both are given.", "EEQ", {})
PARAM(eeq_ppcg_min_atoms, Int, 500, "auto takes ppcg only from this many atoms on; smaller systems use the exact cholesky solve.", "EEQ", {})
PARAM(eeq_ppcg_min_nfrag, Int, 1, "auto takes ppcg only with at least this many fragments. 0: auto never takes it; solve_method ppcg still runs it.", "EEQ", {})
PARAM(eeq_ppcg_tol, Double, 1e-12, "Relative residual tolerance of the CPU ppcg solve, |P r| <= tol*(|b|+1). The GPU ppcg uses pcg_tolerance as an absolute tolerance instead.", "EEQ", {})
PARAM(eeq_ppcg_max_iter, Int, 500, "Iteration cap of ppcg on CPU and GPU; a CPU solve that does not converge uses the exact cholesky solve for that step.", "EEQ", {})
PARAM(eeq_ppcg_precond, String, "localinv", "ppcg preconditioner: 'localinv' (sparse local inverse, fewer iterations, costs a build) or 'diag' (Jacobi). If localinv is not positive definite the solve uses diag with a warning.", "EEQ", {"ppcg_precond"})
PARAM(eeq_localinv_radius_bohr, Double, 8.0, "localinv: keep the inverse entries of atom pairs within this distance (Bohr).", "EEQ", {})
PARAM(eeq_localinv_ext_bohr, Double, 10.0, "localinv: radius (Bohr) of the block inverted per atom, at least eeq_localinv_radius_bohr. The build cost grows with the cube of the atoms inside it.", "EEQ", {})
PARAM(eeq_localinv_rebuild_bohr, Double, 1.0, "localinv: rebuild once an atom moved more than this (Bohr) since the last build; 0 rebuilds at every solve. A stale preconditioner changes the iteration count, not the charges.", "EEQ", {})
PARAM(max_pcg_iterations, Int, 100, "Iteration cap of the pcg solver on the CPU; the GPU pcg uses 200 unless this is given.", "EEQ", {})
PARAM(pcg_tolerance, Double, 1e-10, "Residual tolerance of the pcg solver; on the GPU also the absolute tolerance of ppcg.", "EEQ", {})
PARAM(pcg_large_threshold, Int, 500, "Atom count above which the pcg solver uses the large-system cap and tolerance (pcg_large_system_*).", "EEQ", {})
PARAM(pcg_large_system_iterations, Int, 100, "pcg above pcg_large_threshold atoms: the cap becomes max(max_pcg_iterations, min(pcg_large_system_scaling*N, this)).", "EEQ", {})
PARAM(pcg_large_system_scaling, Int, 10, "pcg above pcg_large_threshold atoms: iterations per atom in the cap, see pcg_large_system_iterations.", "EEQ", {})
PARAM(pcg_large_system_tol_factor, Double, 1e-6, "pcg above pcg_large_threshold atoms: the tolerance becomes max(pcg_tolerance, factor*(|rhs|+1)).", "EEQ", {})
PARAM(eeq_max_iterations, Int, -1, "Older name: a value > 0 sets max_pcg_iterations; -1 leaves it alone.", "EEQ", {})
PARAM(eeq_tolerance, Double, -1.0, "Older name: a value > 0 sets pcg_tolerance; -1 leaves it alone.", "EEQ", {})
PARAM(eeq_accuracy, Double, 1e-6, "Older name of pcg_large_system_tol_factor, which it sets; larger is faster and less accurate.", "EEQ", {})
PARAM(eeq_batched_max_error_eh, Double, 1e-6, "batched solver: warn when the estimated energy error of the dropped inter-fragment coupling exceeds this (Eh). 0 skips the estimate.", "EEQ", {})
PARAM(eeq_distance_cutoff, Double, 0.0, "Drop EEQ matrix elements beyond this distance (Bohr); 0 keeps all, as the reference. A cutoff makes the charges inconsistent with the full Coulomb energy and spoils MD energy conservation.", "EEQ", {})
PARAM(eeq_distance_cutoff_auto, Bool, false, "After Phase 1, set the EEQ cutoff to 30 Bohr when the system is one fragment with all |q| < 0.5 e. An explicit eeq_distance_cutoff wins.", "EEQ", {})
PARAM(eeq_refactor_eps_bohr, Double, 0.0, "Reuse the Cholesky factor of the EEQ matrix until an atom moved more than this (Bohr), corrected by eeq_refine_iters steps. 0 factorises at every call (exact).", "EEQ", {})
PARAM(eeq_refactor_force_every, Int, 0, "With eeq_refactor_eps_bohr > 0: refactorise at least every N calls; 0 = on displacement only.", "EEQ", {})
PARAM(eeq_refine_iters, Int, 1, "With a reused Cholesky factor: iterative-refinement steps per solve, O(N^2) each; more steps leave a smaller charge error. 0 = none.", "EEQ", {})
PARAM(eeq_matrix_rebuild_eps_bohr, Double, 0.0, "Reuse the off-diagonal EEQ Coulomb matrix while atoms moved less than this (Bohr) and the CN changed less than 0.05. 0 rebuilds it at every call.", "EEQ", {})
PARAM(eeq_extrapolation, String, "none", "Start guess of the pcg solver from earlier steps: 'none' (last solution), 'aspc' (fixed-step MD) or 'gauss' (polynomial fit). The converged charges are unchanged.", "EEQ", {})
PARAM(eeq_extrapolation_order, Int, 3, "Order of eeq_extrapolation: for 'aspc' the predictor order k (k+2 steps), for 'gauss' the polynomial degree.", "EEQ", {})
PARAM(use_iterative_refinement, Bool, false, "Repeat the Phase-2 EEQ solve until the charges change less than convergence_threshold or max_iterations is reached. false: one solve.", "EEQ", {})
PARAM(max_iterations, Int, 50, "Most Phase-2 repetitions with use_iterative_refinement; no effect otherwise.", "EEQ", {})
PARAM(convergence_threshold, Double, 1e-6, "use_iterative_refinement: stop once the largest charge change between two Phase-2 solves is below this (e).", "EEQ", {})
PARAM(calculate_cn, Bool, true, "Compute the coordination numbers inside the EEQ solver when the caller passes none.", "EEQ", {})
PARAM(static_charges, Bool, false, "Keep the EEQ charges of the first call for all later steps. Only for equilibrium MD; wrong for charge transfer or ionic dynamics.", "EEQ", {})
PARAM(verbosity, Int, 0, "Print level of the EEQ solver alone (0-3). GFN-FF sets it from -verbosity unless -eeq_solver.verbosity is given.", "Output", {})
PARAM(dump_charges, Bool, false, "Write the Phase-1 and Phase-2 charges of each final-charge solve to charges_dump_N<atoms>.json in the working directory, for diagnosis.", "Output", {})
// ---- Net charge of multi-fragment systems ----
PARAM(frag_charge_model, String, "ensemble", "Net charge of a charged multi-fragment system: 'reference' puts it on the fragment of atom 1 (pprcht/xtb), 'ensemble' weights the placements continuously; CPU only. See FRAG_CHARGE_MODEL.md.", "Fragment charge", {})
PARAM(frag_charge_s_max, Double, 1.0, "ensemble: end of the contact window in units of the bond threshold r_thr; between 1 and s_max close fragments are blended. 1.0 = no window, placement rule only.", "Fragment charge", {})
PARAM(frag_charge_tau, Double, 1.0, "ensemble: temperature (kcal/mol) of the Boltzmann average over chemically identical charge placements; keeps energy and gradient smooth at degeneracies.", "Fragment charge", {})
PARAM(frag_charge_sigma, Double, 0.05, "ensemble: width (e) of the weighting by free Phase-1 EEQ charge between different charge carriers that the parity rule leaves tied.", "Fragment charge", {})
PARAM(frag_charge_max_edges, Int, 4, "ensemble: most fragment contacts blended at once (2^n corners, at most 10); further contacts count as fully separated, with a warning.", "Fragment charge", {})
PARAM(frag_charge_max_placements, Int, 6, "ensemble: most integer charge placements evaluated per corner, preselected by free Phase-1 EEQ charge.", "Fragment charge", {})
PARAM(frag_charge_variant_memory_mb, Int, 8192, "ensemble: memory budget (MB) for the charge variants held at once; above it fewer placements and corners are kept, with a warning. See FRAG_CHARGE_MODEL.md.", "Fragment charge", {})
PARAM(frag_charge_enumerate, Bool, false, "ensemble, validation only: enumerate all charge distributions instead of the bounded selection, which gives the same placements. Refused above 2 million distributions.", "Fragment charge", {})
PARAM(frag_charge_autodetect, Bool, false, "Charged system of exactly two fragments: try both charge placements and keep the lower EEQ energy. Off because the reference never does this and it can pick the wrong fragment.", "Fragment charge", {})
// ---- Topology and coordination number ----
// ORIGIN NOTE (do not duplicate on merge): nh_linear_fix and its guard in
// determineHybridizationFortran() were first added on the `confsearch` branch (51830efa,
// Aug 29 2026) and ported here in Sep 2026. The guard has since been REFINED here (it also
// requires the nitrogen's heavy partner to be branched) and the help text was shortened in
// Oct 2026, so the two branches differ: THIS version is the newer one, take it on merge.
PARAM(topology_mode, String, "constant", "'constant' freezes the topology after setup, as the reference; 'auto' recomputes it once an atom moved > 0.5 Bohr; 'react' re-detects bonds during MD. 'default' = constant.", "Topology", {})
PARAM(cn_cutoff_bohr, Double, 10.0, "Neighbour-list cutoff (Bohr) of the coordination number; 10 matches the reference. 0 uses the cn_accuracy threshold instead.", "Topology", {})
PARAM(cn_accuracy, Double, 1.0, "With cn_cutoff_bohr 0: CN threshold as in the reference, 100 - 50*log10(acc) Bohr^2. 0 takes all pairs.", "Topology", {})
PARAM(static_cn, Bool, false, "Keep CN, its derivatives and the D4 weights and C6 of the first call for all later steps. Only for equilibrium MD; wrong when hybridisation changes.", "Topology", {})
PARAM(static_all, Bool, false, "Shorthand for static_charges and static_cn together.", "Topology", {})
// ---- Reactive topology (topology_mode react, docs/GFNFF_REACT_TOPOLOGY.md) ----
PARAM(react_bond_form_factor, Double, 1.6, "react: a pair becomes bonded when r < factor * covalent-radius sum (element-scaled). Must stay below react_bond_break_factor.", "Reactive", {})
PARAM(react_bond_break_factor, Double, 2.6, "react: a bond breaks when r > factor * covalent-radius sum (element-scaled); the gap to react_bond_form_factor is the hysteresis against flicker.", "Reactive", {})
PARAM(react_check_every, Int, 5, "react: run the O(N^2) bond scan every N energy calls; 0 = on displacement only.", "Reactive", {})
PARAM(react_check_disp_bohr, Double, 0.25, "react: also scan when an atom moved more than this (Bohr) since the last scan; 0 disables this trigger.", "Reactive", {})
PARAM(react_refractory_scans, Int, 10, "react: a pair whose bond just broke may not re-form for this many scans, which stops form/break cycles. 0 disables.", "Reactive", {})
PARAM(react_valence_cap, Bool, true, "react: refuse a new bond while an atom uses its valence plus one exchange slot, counting bond orders. false allows unconstrained formation.", "Reactive", {})
PARAM(react_exchange_scans, Int, 20, "react: an atom may exceed its nominal valence for at most this many scans, then its weakest bond breaks. 0 disables.", "Reactive", {})
PARAM(react_slack_form_factor, Double, 1.2, "react: tighter formation factor for bonds that push an atom into the exchange slot above its sigma valence.", "Reactive", {})
// ---- Hydrogen and halogen bonds ----
PARAM(hb_accuracy, Double, 0.1, "Accuracy of the HB search thresholds, 200 and 400 - 50*log10(acc) Bohr^2; smaller finds more pairs. 0.1 is the reference value.", "HB/XB", {})
PARAM(hb_thr1_bohr2, Double, 0.0, "Override hbthr1, the A-B distance-squared cutoff of the HB search (Bohr^2); 0 derives it from hb_accuracy.", "HB/XB", {})
PARAM(hb_thr2_bohr2, Double, 0.0, "Override hbthr2, the A-H-B sum-of-squares cutoff of the HB search (Bohr^2); 0 derives it from hb_accuracy.", "HB/XB", {})
PARAM(hb_update_rmsd_bohr, Double, 0.3, "Rebuild the HB/XB lists once sqrt(sum d^2)/N since the last build exceeds this (Bohr), as the reference; 0.3 is the reference value.", "HB/XB", {})
PARAM(hb_update_force_every, Int, 10, "Also rebuild the HB/XB lists every N energy calls, since the displacement trigger rarely fires in large systems. 0 = trigger only, 1 = every call.", "HB/XB", {})
PARAM(hb_min_pair_energy_eh, Double, 1e-9, "Drop unbound (case 1) H-bond candidates with exact |E_HB| below this (Eh); donor-bonded cases and the serial small-system search keep all. 0 keeps all.", "HB/XB", {})
PARAM(xb_min_pair_energy_eh, Double, 1e-9, "Drop halogen-bond candidates whose exact |E_XB| is below this (Eh). 0 keeps all.", "HB/XB", {})
PARAM(hb_parallel_min_pairs, Int, 500, "Run the H-bond detection on the thread pool above this many A-B pairs; 0 = always, -1 = never.", "HB/XB", {})
// ---- Dispersion ----
PARAM(dispersion_cutoff_bohr, Double, 0.0, "Cutoff (Bohr) of the D4 dispersion pair list; 0 keeps all pairs, as the reference. A cutoff also widens the CN-derivative range to match.", "Dispersion", {})
PARAM(dispersion_c6_update, Bool, true, "Recompute the D4 pair C6 from the current CN at each new geometry, so MD and optimisation agree with single points. false keeps C6 of the setup geometry.", "Dispersion", {})
PARAM(disp_half_contraction, Bool, true, "Faster per-atom contraction for the D4 C6 and dC6/dCN; it changes the summation order (energies agree to about 1e-10 Eh). false for the plain contraction.", "Dispersion", {})
// ---- Electrostatics ----
PARAM(coulomb_r_cut, Double, 0.0, "Hard per-pair cutoff (Bohr) of the GFN-FF electrostatics; 0 = none, as the reference. A cutoff makes the energy discontinuous and truncates the 1/r sum.", "Coulomb", {})
PARAM(coulomb_implicit, Bool, true, "CPU: evaluate the Coulomb pairs on the fly from per-atom data instead of storing an N^2/2 pair list; same energies. Not used with eeq_distance_cutoff > 0.", "Coulomb", {})
PARAM(farfield_shortcut, Bool, true, "CPU Coulomb energy and EEQ matrix: use 1/r for erf(g r)/r where g r >= 6.5, where erf is exactly 1.0 in FP64. Bit-identical in the cases checked; false evaluates erf everywhere.", "Coulomb", {})
// ---- Pair lists ----
PARAM(nb_cell_list_min_atoms, Int, 800, "From this atom count on, non-bonded neighbours (repulsion, HB/XB, cutoff Coulomb) come from a cell list; below, all pairs are scanned. 0 = always the cell list.", "Pair lists", {"hb_cell_list_min_atoms"})
PARAM(nonbonded_skin_bohr, Double, 2.0, "Verlet skin (Bohr) of the repulsion and cutoff Coulomb lists, rebuilt once an atom moved more than skin/2 (exact). Active from nb_cell_list_min_atoms on; 0 uses nonbonded_rebuild_every.", "Pair lists", {})
PARAM(nonbonded_rebuild_every, Int, 1, "With nonbonded_skin_bohr 0: rebuild the repulsion and cutoff Coulomb lists every N energy calls; 1 is always safe. See GFNFF_PAIR_LIST_REFRESH.md.", "Pair lists", {})
// ---- Implicit solvation (use the dotted -gfnff.solvent; the flat -solvent is shared) ----
PARAM(solvent, String, "none", "Implicit solvent for native GFN-FF ALPB, such as water or dmso; 'none' = gas phase, and the EEQ charges feel the solvent. Use -gfnff.solvent, since the flat -solvent is shared.", "Solvation", {})
PARAM(solvent_model, String, "alpb", "Solvation model. GFN-FF has ALPB only, as the reference: 'gbsa' (or code 2) runs ALPB with a warning.", "Solvation", {})
// ---- Deliberate deviations from the reference (docs/KNOWN_ISSUES_ARCHIVE.md switch table) ----
PARAM(nh_linear_fix, Bool, true, "Keep a 2-coordinate N-H nitrogen from becoming sp by the angle-only rule (> 160 deg), which lets stretched N-H geometries optimise into a too deep artefact. false = reference rule.", "Reference", {})
PARAM(storsion_reference_loop_bug, Bool, false, "Reproduce the reference loop bug that evaluates only the last C#C triple-bond torsion; curcuma sums all of them. true only to match reference totals.", "Reference", {})
PARAM(hh_repulsion_bpair, Bool, true, "Classify H...H 1,3/1,4 repulsion from the bond-pair table, as the reference. false uses a BFS bond count, which can depend on the atom order.", "Reference", {})
PARAM(amideh_acidity_order_bug, Bool, false, "Reproduce the reference atom-order bug in the amide-H donor acidity, which makes the H-bond depend on whether N precedes its H. true only to match reference totals.", "Reference", {})
PARAM(dispersion_atm, Bool, false, "Add a D4 three-body (ATM) term over bonded triples. The reference dispersion is pairwise, so it is off.", "Reference", {})
// ---- GPU ----
PARAM(gpu_block_size, Int, 0, "GPU kernel block size: 0 adaptive, else 32, 64, 128, 256 or 512.", "GPU", {})
PARAM(gpu_cn_pair_regen, Bool, true, "GPU: rebuild the CN-derivative pair list when the topology displacement check fires. false builds it once.", "GPU", {})
PARAM(gpu_cn_pair_regen_every, Int, 0, "GPU: also rebuild the CN-derivative pair list every N gradient steps; 0 = on the displacement check only.", "GPU", {})
PARAM(gpu_cn_pair_cutoff_factor, Double, 2.5, "GPU CN-derivative pair list cutoff, factor*(rcov_i+rcov_j). Larger is slower and more complete, smaller is faster and truncates more.", "GPU", {})
PARAM(gpu_cn_pair_skin_bohr, Double, 1.0, "GPU: Verlet skin (Bohr) of the CN pair list, rebuilt only once atoms moved far enough to lose a term above 1e-29. Negative rebuilds at every topology check.", "GPU", {})
PARAM(gpu_eeq_overlap, Bool, true, "CUDA: the device EEQ solve waits only for its own inputs and overlaps with the dispersion and charge-independent kernels. false waits for all earlier work.", "GPU", {})
PARAM(gpu_dispersion_overlap, Bool, true, "GPU: run the D4 dispersion gradient kernel parallel to the EEQ solve, which it does not need; changes only the summation order. false runs it after.", "GPU", {})
PARAM(gpu_hb_detect_on_device, Bool, true, "GPU: build the H-bond triples on the device at each HB/XB refresh (same set as the host, other order; halogen bonds stay on the host). false = host detection.", "GPU", {})
PARAM(gpu_hb_detect_verify, Bool, false, "With gpu_hb_detect_on_device: also run the host detection at each refresh and report differences (HBVERIFY lines on stderr). Diagnostic only.", "GPU", {})
PARAM(gpu_staged_list_upload, Bool, true, "GPU MD: convert rebuilt pair lists on host threads into pinned memory and copy them asynchronously (CUDA: D4, repulsion, H-bond; ROCm: repulsion). false = serial copies.", "GPU", {})
PARAM(gpu_coulomb_implicit, Bool, true, "GPU: enumerate the Coulomb pairs on the device instead of uploading an N^2/2 host list. Not used with eeq_distance_cutoff > 0.", "GPU", {})
PARAM(gpu_disp_pairs_on_device, Bool, true, "GPU: build the D4 dispersion pair list on the device, in the host order. false builds it on the host.", "GPU", {})
PARAM(gpu_rep_pairs_on_device, Bool, true, "GPU: rebuild the repulsion pair list on the device when the skin rule fires; only for distance-filtered lists (nb_cell_list_min_atoms). false builds it on the host.", "GPU", {})
PARAM(eeq_mixed_precision, Bool, false, "GPU EEQ: factor the matrix in FP32 and refine to FP64 accuracy, for FP64-weak cards. Opt-in; few-fragment paths only.", "GPU", {})
PARAM(eeq_mixed_precision_iters, Int, 2, "Refinement steps of eeq_mixed_precision, at least 1.", "GPU", {})
PARAM(gpu_block_jacobi_max_frag_atoms, Int, 300, "GPU pcg: fragments larger than this get the diagonal instead of the dense block-Jacobi preconditioner, whose build grows steeply with size. See EEQ_SOLVER_CONTROL.md.", "GPU", {})
PARAM(gpu_block_jacobi_max_nfrag, Int, 400, "GPU pcg: above this many fragments no block-Jacobi preconditioner is built (one host round trip per fragment); the diagonal Jacobi is used.", "GPU", {})
PARAM(gpu_split_devices, String, "auto", "CUDA: GPUs to split one molecule over (Coulomb tiles, ppcg EEQ). 'auto' = the -gpu_devices pool or all, unless pinned or a batch worker; a list like '0,1'; 'none'.", "GPU", {})
PARAM(gpu_split_min_atoms, Int, 1000, "CUDA: use the multi-GPU split from this many atoms on.", "GPU", {})
PARAM(gpu_split_eeq_min_atoms, Int, 4000, "CUDA: split the EEQ over several GPUs from this many atoms on; below, the peer latency outweighs the gain. See MULTI_GPU.md.", "GPU", {})
PARAM(gpu_eeq_mixed_precision, Bool, true, "CUDA multi-GPU EEQ: PCG matvec on an FP32 copy, restarted from the FP64 residual until the FP64 tolerance is met. false (or CURCUMA_EEQ_MIXED=0) = FP64.", "GPU", {})
PARAM(gpu_eeq_mixed_rho, Double, 1e-5, "gpu_eeq_mixed_precision: each FP32 round stops at this fraction of its starting residual. Larger values need more FP64 restarts.", "GPU", {})
PARAM(rocm_dc6dcn_tables, Bool, false, "ROCm: per-pair C6 and dC6/dCN from per-atom partner tables (the CUDA scheme) instead of the fused gather, which is faster. CURCUMA_DC6DCN_TABLES=0|1 overrides it.", "GPU", {})
PARAM(rocm_device_schur, Bool, true, "ROCm: do the Schur step of the dense EEQ solve on the device, so only nfrag + nfrag^2 numbers reach the host. false does it on the host.", "GPU", {})
PARAM(rocm_fixed_point, Bool, true, "ROCm: sum forces, dE/dCN, energies and CN in 64-bit fixed point, so repeated runs are bitwise identical. false uses FP64 atomics.", "GPU", {})
PARAM(eeq_rocm_cpu_fragment_threshold, Int, 16, "ROCm: from this fragment count on, a stalled device EEQ falls back to the CPU; giving the key explicitly forces the exact CPU EEQ there. 0 = always the device.", "GPU", {})
END_PARAMETER_DEFINITION

class GFNFF {
public:
    /// Test hook: #EEQ PCG solves that used the multi-step warm-start extrapolation
    /// (eeq_extrapolation). 0 with the default 'none'. Claude Generated.
    long eeqPcgExtrapolationCount() const {
        return m_eeq_solver ? m_eeq_solver->pcgExtrapolationCount() : -1;
    }

    /// Test hook: drop the EEQ Cholesky/matrix caches, exactly as getCachedTopology()
    /// does before a full topology rebuild. Lets the re-entrancy regression test
    /// reproduce the MD/opt rebuild path in-process. Claude Generated (Jul 2026).
    void invalidateEEQCachesForTest() {
        if (m_eeq_solver) {
            m_eeq_solver->invalidateCholeskyCache();
            m_eeq_solver->invalidateMatrixCache();
        }
    }

    /**
     * @brief Static topology data — computed once at initialization, never changes
     *
     * Claude Generated (March 2026): Separated from TopologyInfo for clean
     * distinction between one-time topology and per-step dynamic data.
     *
     * Contains: connectivity, hybridization, ring membership, Phase-1 EEQ charges,
     * correction parameters, BATM topology, cached EEQ element parameters.
     */
    struct GFNFFTopology {
        // Atom classification
        Vector neighbor_counts;                                  // Simple neighbor counts (integer CN)
        std::vector<int> hybridization;                          // 0=none/octahedral, 1=sp, 2=sp2, 3=sp3, 5=hypervalent (determineHybridizationFortran)
        std::vector<int> pi_fragments;                           // Pi fragment assignment per atom
        std::vector<int> pi_atoms_final;                         // Fortran post-Hueckel piadr (gfnff_ini.f90:1016, "piadr = itmp"): 1 iff the atom ends a bond inside a SOLVED pi-system. Stricter than pi_fragments and the array every consumer after the Hueckel section tests - Claude Generated Sep 2026
        std::vector<int> itag;                                   // -1 iff atom is eta-coordinated to a metal (Fortran itag; gfnff_ini2.f90:170-198) - Claude Generated Jul 2026
        std::vector<int> pi_system_charge;                       // ipis: charge per pi-system (subtract from nelpi) - Claude Generated Jul 2026
        std::vector<int> ring_sizes;                             // Smallest ring containing each atom
        std::vector<bool> is_metal;                              // Metal atom flags
        std::vector<bool> is_aromatic;                           // Aromatic atom flags
        std::vector<bool> is_amide_h;                            // Amide hydrogen flags (Coulomb chi correction)

        // Connectivity
        std::vector<std::vector<int>> neighbor_lists;            // Full neighbor connectivity
        std::vector<std::vector<int>> adjacency_list;            // Per-atom bonded neighbor list
        SparseTopoTable topo_distances;                          // BFS shortest-path bond counts, stored up to 5 bonds (0 = self, 999 = further/unconnected) - sparse since Sep 2026

        // Fortran multi-list neighbour construction (gfnff_ini2.f90:128-130, 197-202).
        // Fortran keeps four lists and assigns hybridization from the metal-reduced,
        // eta-aware mixture rather than from the full connectivity. Claude Generated (Jul 2026).
        std::vector<std::vector<int>> nb_full;                   // nbf: getnb(icase=1), no filtering
        std::vector<std::vector<int>> nb_hc;                     // topo%nb DURING the hyb loop: getnb(icase=2), drops all bonds of highly-coordinated atoms
        std::vector<std::vector<int>> nb_nometal;                // nbm: getnb(icase=3), metals + unusually coordinated heavy atoms removed
        std::vector<double> metallic_character;                  // mchar (gfnff_ini.f90:249), gates the nbm metal filter

        // Functional groups
        std::vector<FunctionalGroupType> functional_groups;      // Per-atom classification

        // Phase-1 EEQ charges and corrections (computed once, fixed)
        Vector topology_charges;                                 // Phase 1 EEQ charges (qa)
        Vector eeq_charges;                                      // Phase 1 EEQ charges (alias)
        Vector dxi;                                              // Electronegativity corrections
        Vector dgam;                                             // Hardness corrections
        Vector dalpha;                                           // Polarizability corrections
        Vector alpeeq;                                           // Charge-corrected alpha² values

        // Bond classification
        std::vector<int> bond_types;                             // Per-bond type (1-7)
        std::vector<double> pi_bond_orders;                      // Triangular lin(i,j) format

        // Cached EEQ element parameters (performance optimization)
        std::vector<double> eeq_chi;                             // Electronegativity per atom
        std::vector<double> eeq_gam;                             // Chemical hardness per atom
        std::vector<double> eeq_alp;                             // Damping parameter (squared) per atom
        std::vector<double> eeq_cnf;                             // CN correction factor per atom

        // BATM topology
        SparseTopoTable bpair;                                   // Reference topo%bpair (nbondmat): 1/2/3 stored, 5 = further, 0 = self - sparse since Sep 2026
        std::vector<std::tuple<int,int,int>> b3list;             // Batm triples (i,j,k)
        int nbatm = 0;

        // Ring enumeration
        std::vector<std::vector<int>> rings;                     // rings[ring_id] = {atom0, ...}
        std::vector<std::vector<int>> atom_to_rings;             // atom_to_rings[atom] = {ring_id0, ...}

        // Molecular fragments
        int nfrag = 1;
        std::vector<int> fraglist;                               // Fragment ID per atom (1-indexed)
        std::vector<double> qfrag;                               // Target charge per fragment
    };

    /**
     * @brief Per-geometry-step dynamic data — recalculated when geometry changes
     *
     * Claude Generated (March 2026): Separated from TopologyInfo for clean
     * distinction between one-time topology and per-step dynamic data.
     *
     * Contains: coordination numbers, distance matrices.
     * Phase-2 EEQ charges are NOT stored here — they are managed by GFNFF::Calculation().
     */
    struct GFNFFDynamicState {
        Vector coordination_numbers;                             // D3 CN (geometry-dependent)
        Eigen::MatrixXd distance_matrix;                         // N×N distances in Bohr (only for initial topology, not per-step)
    };

    /**
     * @brief Combined topology + dynamic state (backward-compatible wrapper)
     *
     * Inherits from both GFNFFTopology and GFNFFDynamicState so that all existing
     * code accessing `topo_info.hybridization`, `topo_info.coordination_numbers`, etc.
     * continues to work unchanged. New code can use the base types directly when
     * only static or dynamic data is needed.
     *
     * Claude Generated (March 2026): Architecture cleanup — inheritance-based split
     */
    struct TopologyInfo : public GFNFFTopology, public GFNFFDynamicState {
    };

    /**
     * @brief GFN-FF Results structure (matches Fortran gfnff_results)
     *
     * Claude Generated (Mar 2026): Unified energy and gradient decomposition
     * Reference: external/gfnff/src/gfnff_engrad.F90:35-63 (gfnff_results type)
     *
     * Provides complete decomposition of energy and gradient components for:
     * - Debugging and validation
     * - Energy component export for external analysis
     * - Restart file compatibility
     * - Solvation energy decomposition (when ALPB active)
     */
    struct GFNFFResults {
        // Total energy and gradient norm
        double e_total = 0.0;      ///< Total energy in Hartree
        double gnorm = 0.0;         ///< Gradient norm in Hartree/Bohr

        // Bonded energy components (Hartree)
        double e_bond = 0.0;        ///< Bond stretching
        double e_angle = 0.0;      ///< Angle bending
        double e_torsion = 0.0;    ///< Dihedral torsion (primary + extra)
        double e_inversion = 0.0; ///< Out-of-plane bending
        double e_storsion = 0.0;  ///< Triple bond torsions (sTors_eg)

        // Non-bonded energy components (Hartree)
        double e_repulsion = 0.0;     ///< Non-bonded repulsion (exponential)
        double e_bonded_repulsion = 0.0; ///< Bonded repulsion (1,2/1,3 scaling)
        double e_coulomb = 0.0;       ///< Electrostatic (EEQ charges)
        double e_dispersion = 0.0;    ///< D3/D4 dispersion
        double e_hb = 0.0;            ///< Hydrogen bonds
        double e_xb = 0.0;            ///< Halogen bonds

        // Three-body dispersion (Hartree)
        double e_atm = 0.0;      ///< D3/D4 ATM three-body dispersion
        double e_batm = 0.0;     ///< Bonded ATM (1,4-pairs)

        // External contributions (future)
        double e_ext = 0.0;      ///< External field (reserved)

        // Solvation energy components (when ALPB active, in Hartree)
        double g_born = 0.0;     ///< Born electrostatic solvation
        double g_sasa = 0.0;     ///< Non-polar surface area
        double g_hb_solv = 0.0;  ///< HB solvation correction
        double g_shift = 0.0;    ///< Free energy shift
        double g_solv = 0.0;     ///< Total solvation energy

        // Dipole moment (Debye)
        Eigen::Vector3d dipole = Eigen::Vector3d::Zero();

        // Gradient components (Hartree/Bohr, 3×N matrices)
        // Only populated when gradient calculation was requested
        Eigen::MatrixXd g_bond;       ///< Bond gradient (3×N)
        Eigen::MatrixXd g_angle;      ///< Angle gradient (3×N)
        Eigen::MatrixXd g_torsion;    ///< Torsion gradient (3×N)
        Eigen::MatrixXd g_repulsion;  ///< Repulsion gradient (3×N)
        Eigen::MatrixXd g_coulomb;    ///< Coulomb gradient (3×N)
        Eigen::MatrixXd g_dispersion; ///< Dispersion gradient (3×N)
        Eigen::MatrixXd g_hb;         ///< Hydrogen bond gradient (3×N)
        Eigen::MatrixXd g_xb;         ///< Halogen bond gradient (3×N)
        Eigen::MatrixXd g_atm;        ///< ATM gradient (3×N)
        Eigen::MatrixXd g_batm;       ///< BATM gradient (3×N)
        Eigen::MatrixXd g_total;      ///< Total gradient (3×N)

        // Atomic charges (EEQ)
        Eigen::VectorXd charges;  ///< Phase-2 EEQ charges

        /**
         * @brief Export to JSON for serialization
         * @return JSON object with all energy components
         */
        json toJSON() const;

        /**
         * @brief Import from JSON for restart
         * @param j JSON object with results
         */
        void fromJSON(const json& j);
    };

    /**
     * @brief Default constructor
     */
    GFNFF();

    /**
     * @brief Constructor with custom parameters
     * @param parameters JSON configuration for GFN-FF
     */
    explicit GFNFF(const json& parameters);

    /**
     * @brief Destructor
     */
    virtual ~GFNFF();

    /**
     * @brief Initialize molecule for GFN-FF calculation from Mol object
     * @param molecule Molecule to initialize
     * @return true if initialization successful
     */
    bool InitialiseMolecule(const Mol& molecule);

    /**
     * @brief Initialize molecule for GFN-FF calculation (parameterless)
     * @return true if initialization successful
     */
    bool InitialiseMolecule();

    /**
     * @brief Update molecular geometry
     * @param geometry New geometry matrix
     * @return true if update successful
     */
    bool UpdateMolecule(const Matrix& geometry);

    /**
     * @brief Update molecular geometry (parameterless)
     * @return true if update successful
     */
    bool UpdateMolecule();

    /**
     * @brief Perform GFN-FF calculation
     * @param gradient Calculate gradients if true
     * @return Total energy in Hartree
     */
    double Calculation(bool gradient = false);

    /**
     * @brief Get analytical gradients
     * @return Gradient matrix (N_atoms x 3) in Hartree/Bohr
     */
    Geometry Gradient() const { return m_gradient; }

    /**
     * @brief Check if gradients are available
     * @return true (GFN-FF always provides gradients)
     */
    bool hasGradient() const { return true; }

    /**
     * @brief Calculate numerical gradient via finite differences
     * @param dx Displacement step size (default 1e-5 Bohr)
     * @return Gradient matrix (N_atoms x 3) in Hartree/Bohr
     *
     * Claude Generated (Feb 21, 2026): Diagnosis infrastructure for gradient validation.
     * Perturbs each atom coordinate and recalculates energy including EEQ charge update.
     * This captures the full geometric dependence (including dq/dx terms) for comparison
     * with analytical gradient.
     *
     * Reference: Plan unified-baking-gizmo.md Step 1a
     */
    Matrix NumGrad(double dx = 1e-5);

    /**
     * @brief Numerical gradient with FIXED charges and CN (no EEQ/CN recalculation)
     * @param dx Displacement step size (default 1e-5 Bohr)
     * @return Gradient matrix (N_atoms x 3) in Hartree/Bohr
     *
     * Claude Generated (Feb 23, 2026): Isolates gradient formula bugs from missing dq/dx.
     * Only updates geometry in ForceField, does NOT recalculate CN or EEQ charges.
     * Comparison with analytical gradient tests ONLY the direct gradient formulas
     * and CN chain-rule terms. Any deviation indicates a real gradient bug.
     */
    Matrix NumGradFixedCharges(double dx = 1e-5);

    /**
     * @brief Get atomic partial charges (Phase 2 energy charges - nlist%q)
     * @return Vector of atomic charges (final energy charges)
     *
     * Returns Phase 2 EEQ charges used for energy calculation.
     * These are the charges used in gradient calculations and for Coulomb energy.
     */
    Vector Charges() const;

    /**
     * @brief Get topology charges (Phase 1 - topo%qa)
     * @return Vector of topology charges
     *
     * Returns Phase 1 EEQ topology charges used for parameter generation.
     * These charges use integer neighbor count and are used to compute
     * corrections for bonds, angles, and other topological terms.
     *
     * Claude Generated (January 4, 2026)
     */
    Vector getTopologyCharges() const;

    /**
     * @brief Get energy charges (Phase 2 - nlist%q)
     * @return Vector of energy charges
     *
     * Alias for Charges() - returns Phase 2 EEQ charges.
     * Provided for clarity when comparing both charge types.
     *
     * Claude Generated (January 4, 2026)
     */
    Vector getEnergyCharges() const { return Charges(); }

    /**
     * @brief Get bond orders (Wiberg bond orders)
     * @return Vector of bond orders
     */
    Vector BondOrders() const;

    /**
     * @brief Get complete GFN-FF results structure
     * @return GFNFFResults with all energy and gradient components
     *
     * Claude Generated (Mar 2026): Unified energy decomposition
     * Reference: Fortran gfnff_engrad.F90:35-63 (gfnff_results type)
     *
     * Provides all energy components matching Fortran output:
     * - Bonded: bond, angle, torsion, inversion, storsion
     * - Non-bonded: repulsion, coulomb, dispersion, hb, xb
     * - Three-body: atm, batm
     * - Solvation: g_born, g_sasa, g_hb, g_shift (if ALPB active)
     *
     * Gradient components are populated if gradient=true in Calculation().
     */
    GFNFFResults getResults() const;

    /**
     * @brief Set calculation parameters
     * @param parameters JSON configuration
     */
    void setParameters(const json& parameters);

    /**
     * @brief Set thread count used for CN/EEQ/D4 phases and parameter generation
     *
     * Claude Generated (WP1, May 2026): replaces ad-hoc `m_parameters.value("threads", 1)`
     * reads scattered through Calculation()/prepareCNAndEEQ()/parameter generation. The
     * QM-wrapper (GFNFFComputationalMethod::setThreadCount) forwards into here so that
     * post-construction thread count changes from EnergyCalculator are honored.
     */
    void setThreadCount(int threads) { m_threads = (threads > 0 ? threads : 1); m_parameters["threads"] = m_threads; }

    /**
     * @brief Get currently configured thread count
     */
    int threadCount() const { return m_threads; }

    /**
     * @brief Fix the EEQ fragment grouping and the integer group charges of this instance
     *
     * Used by the frag_charge_model ensemble: each charge variant is a GFNFF instance whose
     * Phase-1/Phase-2 EEQ constraints are given here instead of being perceived + placed by the
     * reference rule. fraglist is 1-based per atom (group id), qfrag[g-1] the group charge.
     * Must be called before InitialiseMolecule(). Claude Generated (Sep 2026).
     */
    void setFragmentOverride(const std::vector<int>& fraglist, int nfrag, const std::vector<double>& qfrag);

    /// frag_charge_model ensemble: number of charge variants evaluated in the last call (0 = inactive)
    int fragEnsembleVariantCount() const { return m_frag_last_nvariants; }

    /// frag_charge_model resolved to ensemble (false: reference rule, or react topology mode). Claude Generated (Oct 2026).
    bool fragEnsembleRequested() const { return m_frag_ensemble; }
    /// True when an implicit solvent (ALPB) is active; the GPU wrappers do not apply it. Claude Generated (Oct 2026).
    bool solvationActive() const { return m_solvation != nullptr; }

    /**
     * @brief Charge placements of the frag_charge_model ensemble for one corner (bounded selection)
     *
     * A placement distributes nq = |charge| units (sign `unit`, +1 or -1) over G groups. Rank per
     * placement = sum over groups of [odd electron count] + (G+1)*[no electron left], with
     * nelec[g] the electrons of the NEUTRAL group; score S = sum_g (unit*c_g) * qfree_g[g] / nq
     * (summed in group order, as the full enumeration did). Returns the minimum-rank placements:
     * all of them in the old enumeration order (lexicographically descending) if there are at
     * most max_placements, else the max_placements best by S (ties: enumeration order).
     * Dynamic programming over (group, units used) with pruning by safe dominance:
     * O(G * nq^2 * max_placements) time, O(G * nq) memory; bit-identical to
     * fragEnumeratePlacements() (equivalence test: test_frag_charge_selection).
     * exact_tiebreak is set false only if a pathological cluster of near-equal scores (within
     * rounding) had to be truncated, where the tie order may then differ from the enumeration.
     * Claude Generated (Oct 2026).
     */
    static std::vector<std::vector<int>> fragSelectPlacements(int nq, double unit, const std::vector<long>& nelec,
        const std::vector<double>& qfree_g, int max_placements, bool* exact_tiebreak = nullptr);

    /**
     * @brief Reference for fragSelectPlacements(): the pre-Oct-2026 full enumeration
     *
     * Enumerates all C(nq+G-1, nq) compositions (memory grows like G^nq), keeps the minimum
     * rank, stable-sorts by S if more than max_placements remain and truncates. Validation and
     * tests only (-gfnff.frag_charge_enumerate). Claude Generated (Oct 2026), code moved from
     * fragEnsembleBlend().
     */
    static std::vector<std::vector<int>> fragEnumeratePlacements(int nq, double unit, const std::vector<long>& nelec,
        const std::vector<double>& qfree_g, int max_placements);

    /**
     * @brief Retrieve cached topology information, computing it once if needed
     * @return Reference to topology information
     *
     * Provided for validation and testing purposes.
     */
    const TopologyInfo& getTopologyInfo() const { return getCachedTopology(); }

    /// WP-A (Jun 2026): set by the GPU method when gpu_disp_pairs_on_device is on, so
    /// the host D4 generator computes only CN + Gaussian weights and skips its O(N^2)
    /// pair loop (the GPU builds the pair list on device). Must be set before
    /// InitialiseMolecule. No effect on the CPU path.
    void setSkipHostDispPairs(bool v) { m_skip_host_disp_pairs = v; }

    /// Claude Generated (Sep 2026): set by the GPU method so the host does not build the
    /// N^2/2 Coulomb pair list (the device enumerates the pairs). Only honoured without an
    /// EEQ distance cutoff. Must be set before InitialiseMolecule. No effect on the CPU path.
    void setImplicitCoulombPairs(bool v) { m_implicit_coulomb_pairs = v; }

    /// Claude Generated (Oct 2026): set by a GPU method that never runs this object's CPU
    /// workspace (ROCm wrapper): updateHBXBIfNeeded() then re-detects the HB/XB lists for
    /// getLastHBonds()/getLastXBonds() but does not copy them into m_workspace nor rebuild its
    /// bond-HB data (12 + 16 ms per re-detection at 582k H bonds). No effect on the CPU path.
    void setSkipHostHBWorkspaceUpdate(bool v) { m_skip_host_hb_workspace = v; }
    /// Claude Generated (Oct 2026): GPU wrappers read the repulsion lists from getLast*Repulsions() and never run the CPU workspace; skips the 1.5 M-entry copy into it per rebuild.
    void setSkipHostRepulsionWorkspaceUpdate(bool v) { m_skip_host_rep_workspace = v; }

    /**
     * @brief Export topology information for restart/topology I/O
     * @return JSON object with topology data (fragments, charges, hybridization, CN)
     *
     * Claude Generated (Mar 2026): Phase 3 - Restart/Topology I/O
     * Reference: Fortran gfnff_restart.f90 (write_restart_gff)
     *
     * Exports topology information needed for MD restart:
     * - Fragment information (nfrag, fraglist, qfrag)
     * - Topology charges (Phase-1 EEQ)
     * - Hybridization states
     * - Coordination numbers
     * - Ring membership
     */
    json exportTopology() const;

    /**
     * @brief Import topology information from restart file
     * @param topo_json JSON object with topology data
     * @return true if topology was successfully imported
     *
     * Claude Generated (Mar 2026): Phase 3 - Restart/Topology I/O
     */
    bool importTopology(const json& topo_json);

    /**
     * @brief Compute fingerprint for topology cache validation
     * @return Hash string based on atom count, types, and bond list
     *
     * Claude Generated (March 2026): Fingerprint ensures topology cache is
     * invalidated when molecular connectivity changes.
     */
    std::string computeTopologyFingerprint() const;

    /**
     * @brief Calculate full topology information for advanced parametrization
     * @return Complete topology information
     */
    /// Full topology build. Runs calculateTopologyInfoOnce() twice, mirroring Fortran's
    /// q-loop (gfnff_ini.f90:258-263): pass 1 with qa=0, pass 2 with the pass-1 charges
    /// shrinking the bond radii. Claude Generated (Jul 2026).
    TopologyInfo calculateTopologyInfo() const;

    /// One pass of the topology build (bond list -> four neighbour lists -> hybridization
    /// -> rings/pi/EEQ). Uses m_bond_qa for the getnb radius shrink.
    TopologyInfo calculateTopologyInfoOnce() const;

    /**
     * @brief Generate GFN-FF parameters as native C++ structs (no JSON)
     * @return GFNFFParameterSet with all interaction parameters
     *
     * Claude Generated (March 2026): Primary parameter generation path.
     * Called after InitialiseMolecule() when topology is available.
     */
    GFNFFParameterSet generateGFNFFParameterSet();


    /**
     * @brief Consume cached parameter set for external use.
     *
     * Claude Generated (March 2026): initializeForceField() stores a heap copy of the
     * parameter set. External consumers (GPU wrapper, etc.) call consumeCachedParameterSet()
     * once to take ownership; subsequent calls return nullptr.
     */
    /**
     * @brief Consume cached parameter set (avoids extra generateGFNFFParameterSet call).
     *
     * Claude Generated (March 2026): initializeForceField() stores a heap copy of the
     * parameter set. External callers (e.g. GPU wrapper) call this once to take ownership;
     * subsequent calls return nullptr.
     */
    std::unique_ptr<GFNFFParameterSet> consumeCachedParameterSet() { return std::move(m_cached_parameter_set); }

    // === GPU orchestration helpers (Claude Generated March 2026) ===
    // These expose internal CN/EEQ computation so that GGFNFFComputationalMethod
    // can orchestrate GPU + CPU-residual without duplicating logic.

    /**
     * @brief Timing breakdown for prepareCNAndEEQ sub-steps.
     * Claude Generated (April 2026): Collected for consolidated summary in Calculation().
     */
    struct PrepTiming {
        double total = 0.0;
        double cn = 0.0;
        double eeq_topo = 0.0;
        double cnf = 0.0;
        double dcn = 0.0;
        double d4_gw = 0.0;
        double eeq_solve = 0.0;
        double charge_dist = 0.0;
    };

    /**
     * @brief Compute CN, EEQ charges, and (if gradient) CN derivatives for current geometry.
     * Results are stored internally and distributed to m_workspace.
     * Call getters below to retrieve results for external workspaces.
     * @param gradient  If true, also compute gradient-related data (cnf, dc6dcn)
     * @param gpu_only  If true, skip sparse dcn matrix build and CPU forcefield/workspace
     *                  distribution (GPU has its own k_cn_chainrule kernel)
     * @param external_cn  If non-null, use these CN values instead of computing on CPU.
     *                     Claude Generated (March 2026): Enables GPU CN bypass.
     * @param out_timing  If non-null, filled with per-sub-step durations (ms).
     */
    void prepareCNAndEEQ(bool gradient, bool gpu_only = false, const Vector* external_cn = nullptr, bool skip_eeq = false, PrepTiming* out_timing = nullptr);

    /**
     * @brief Extract EEQ parameters for GPU solver (O(N) CPU work only)
     *
     * Claude Generated (March 2026): GPU EEQ Phase 7 — parameter extraction.
     * Computes dxi, dgam, alpha_corrected, gam_corrected, rhs_atoms from
     * cached topology and current CN. No matrix build or solve — that's done
     * on GPU by EEQSolverGPU.
     *
     * Must be called AFTER prepareCNAndEEQ(gradient, gpu_only=true, cn, skip_eeq=true).
     *
     * @param cn Current coordination numbers (from GPU or CPU)
     * @return EEQGPUParams struct with all arrays for GPU solver
     */
    struct EEQGPUParams {
        std::vector<double> alpha_corrected;  ///< [N] charge-corrected alpha² values
        std::vector<double> gam_corrected;    ///< [N] corrected hardness
        std::vector<double> rhs_atoms;        ///< [N] electronegativity RHS: -chi + dxi + cnf*sqrt(cn)
        std::vector<double> rhs_constraints;  ///< [nfrag] target charges per fragment
        std::vector<int>    fraglist;         ///< [N] fragment ID per atom (1-indexed)
        int nfrag = 1;                        ///< Number of molecular fragments
        // WP2: topology-constant RHS components for GPU kernel k_build_eeq_rhs
        std::vector<double> chi_corrected_static; ///< [N] -chi + dxi + amide_corr (no CN term)
        std::vector<double> cnf;                  ///< [N] cnf_eeq per atom
    };
    EEQGPUParams prepareEEQParametersForGPU(const Vector& cn) const;

    /**
     * @brief Get EEQ distance cutoff from parameters (Bohr).
     * Used by GPU path to match CPU matrix conditioning behavior.
     */
    double getEEQDistanceCutoff() const;

    /**
     * @brief Forward gfnff-level solver parameters to eeq_solver sub-config.
     * Handles: solve, eeq_max_iterations, eeq_tolerance, eeq_accuracy, eeq_distance_cutoff.
     */
    void forwardEEQSolverParams(json& eeq_params);

    /**
     * @brief WP-S3 (May 2026): apply eeq_distance_cutoff_auto heuristic after Phase-1.
     * Pre-condition: m_cached_topology contains valid topology_charges and nfrag.
     * Sets EEQSolver cutoff to 30 Bohr if nfrag==1 and max|q|<0.5 e, else clears the
     * override. Honours an explicit user-set eeq_distance_cutoff>0 (manual wins).
     */
    void applyEEQCutoffAutoIfRequested();

    /**
     * @brief Re-detect HB/XB pairs if geometry has changed enough (RMSD > 0.3 Bohr).
     * Updates the given workspace with new HB/XB interaction lists.
     * @param ws External workspace to update (e.g. CPU residual workspace for GPU path)
     */
    void updateHBXBIfNeeded(FFWorkspace* ws);

    /// True if updateHBXBIfNeeded() ran and changed lists since last call
    bool consumeHBXBUpdate() { bool r = m_hbxb_updated; m_hbxb_updated = false; return r; }

    // ---- Device-side H-bond detection (Claude Generated, Oct 2026; CUDA and ROCm, -gfnff.gpu_hb_detect_on_device) ----

    /// Per-atom inputs of the H-bond detection in flat arrays, the part of detectHydrogenBondsNative() that
    /// is O(N): basicity/acidity with the functional-group overrides, the HB hydrogens, the donor/acceptor
    /// atom filter, the thresholds, and the neighbour / bond / hydrogen-on-atom lists as CSR. The device
    /// enumerates the (A, H, B) triples from them (FFWorkspaceGPU::detectHBondsOnDevice()).
    using HBDetectInputs = GFNFFHBDetectInputs;
    HBDetectInputs prepareHBDetectInputs(const Vector& charges) const;

    /// The trigger part of updateHBXBIfNeeded() (periodic forced rebuild or RMSD test), with its bookkeeping.
    /// True when the lists have to be rebuilt now.
    bool hbxbRefreshDue();

    /// Halogen bonds from the host detection (they stay on the host: a few entries).
    std::vector<GFNFFHalogenBond> detectXBForDevice() const;

    /// Host detection of the H-bond list as the reference of the device list (verification only).
    std::vector<GFNFFHydrogenBond> detectHBHostReference() const;

    /// Bookkeeping after a device-side re-detection: stores the XB list, the counts and the reference geometry,
    /// invalidates the host H-bond list and flags the update (consumeHBXBUpdate()).
    void commitHBXBDeviceUpdate(int n_hb, std::vector<GFNFFHalogenBond>&& xbonds);

    /**
     * @brief Periodically rebuild the bonded/non-bonded repulsion pair lists from the
     * current geometry.
     *
     * Claude Generated (Sep 2026): generateRepulsionPairsNative() builds the non-bonded
     * repulsion list ONCE, at InitialiseMolecule() time, using a hard 20 Bohr distance
     * cutoff (NB_REP_RCUT in gfnff_method.cpp) — a pair further apart than that at t=0
     * contributes exactly zero repulsion for the rest of the run, even after the two
     * atoms diffuse within bonding distance of each other. Unlike the HB/XB list
     * (updateHBXBIfNeeded()) there was no periodic refresh at all, so an atom pair that
     * starts outside the cutoff (common in a large, loosely packed system — see
     * docs/MD_LARGE_SYSTEMS.md, the isolated-water-in-a-cavity case) can pass straight
     * through the repulsive wall in MD, producing an unphysical near-zero contact
     * distance and, eventually, an EEQ/energy blow-up.
     *
     * Fix (design 1 of the two considered — see the task write-up): re-derive the pair
     * list from the CURRENT geometry every `nonbonded_rebuild_every` energy evaluations
     * (default 1 = every step). generateRepulsionPairsNative() is side-effect-free (reads
     * the already-cached bond list + topology, returns fresh vectors) and does not touch
     * generateGFNFFParameterSet() (the "third call causes heap corruption" function noted
     * at its call site) — it regenerates only the two repulsion vectors, nothing else.
     * Rebuilding unconditionally at the default interval is simplest and provably correct
     * (no pair can be missed for longer than one rebuild interval); a Verlet-style skin
     * list (rebuild at a wider cutoff, trigger on max displacement) was the documented
     * fallback if the per-call cost turned out to be prohibitive — measured cheap enough
     * on the reference case that the fallback was not needed.
     *
     * @param ws External workspace to update as well (e.g. CPU residual workspace for the
     *           GPU path, mirroring updateHBXBIfNeeded()'s `ws` parameter).
     */
    void updateNonbondedRepulsionIfNeeded(FFWorkspace* ws);

    /// True if updateNonbondedRepulsionIfNeeded() rebuilt the lists since last check
    bool consumeNonbondedRepulsionUpdate() { bool r = m_nb_rep_updated; m_nb_rep_updated = false; return r; }

    /**
     * @brief Per-atom and per-pair tables from which a device rebuilds the non-bonded repulsion list.
     *
     * Claude Generated (Oct 2026): the pair term of generateRepulsionPairsNative() is
     * alpha = sqrt(dum_i * dum_j) * ff, repab = repz_i * repz_j * REPSCALN, r_cut = 20 Bohr, with
     * dum_i = repan[Z_i] * (1 + q_i * QREPSCAL) * (1 + NREPSCAL / (1 + cn_i^2)) per atom and ff a
     * function of the element pair (H-H: HHFAC, times HH13REP / HH14REP for topological distance
     * 2 / 3). Everything except the distance test depends on the topology only, so these tables
     * stay valid until the next topology rebuild.
     */
    struct RepulsionDeviceInputs {
        std::vector<int> valid;             ///< [N] 1 if the element has repulsion parameters
        std::vector<int> metal;             ///< [N] 1 if PeriodicTable::getMetalType(Z) > 0
        std::vector<double> dum;            ///< [N] per-atom alpha factor
        std::vector<double> repz;           ///< [N] effective nuclear charge factor
        std::vector<long long> bond_keys;   ///< sorted i * N + j of bond-list entries with i < j (excluded pairs)
        std::vector<long long> hh_keys;     ///< sorted i * N + j (i < j) of H-H pairs at topological distance 2 or 3
        std::vector<int> hh_dist;           ///< distance (2 or 3) belonging to hh_keys
        double repscaln = 0.0, hhfac = 0.0, hh13 = 0.0, hh14 = 0.0;
        double build_radius_bohr = 0.0;     ///< 20 Bohr kernel cutoff + Verlet skin
        double r_cut_bohr = 20.0;           ///< value stored as GFNFFRepulsion::r_cut
    };

    /// Fills the tables above; false if the topology cannot be mapped (the caller keeps the host build).
    bool repulsionDeviceInputs(RepulsionDeviceInputs& out) const;

    /// Claude Generated (Oct 2026): the GPU wrapper rebuilds the non-bonded repulsion list on the device;
    /// updateNonbondedRepulsionIfNeeded() then only decides WHEN (skin rule, counters) and flags the update.
    void setSkipHostRepulsionPairs(bool v) { m_skip_host_rep_pairs = v; }

    /// Last rebuilt bonded/non-bonded repulsion pair lists (for GPU SoA re-upload)
    const std::vector<GFNFFRepulsion>& getLastBondedRepulsions() const { return m_last_bonded_reps; }
    const std::vector<GFNFFRepulsion>& getLastNonbondedRepulsions() const { return m_last_nonbonded_reps; }

    /**
     * @brief Rebuild the D4 dispersion pair list when an atom may have crossed into its
     * evaluation cutoff.
     *
     * Claude Generated (Sep 2026): the D4 pair list (D4ParameterGenerator::
     * GenerateDispersionPairsNative()) is built ONCE, at InitialiseMolecule() time, from every
     * pair closer than R_build = 60 Bohr, and the kernel evaluates each stored pair up to
     * r_cut = R_eval = 50 Bohr. A pair further apart than 60 Bohr at t=0 was absent for the
     * rest of the run, even after its atoms diffused together — the same one-shot architecture
     * as the repulsion list (updateNonbondedRepulsionIfNeeded()); the HB/XB update log even
     * labelled the dispersion count "(static)".
     *
     * Unlike repulsion, this list already carries a skin of R_build - R_eval = 10 Bohr. A pair
     * that was outside R_build at the last build can only come inside R_eval after the two
     * atoms approached each other by the full skin, which needs at least one of them to move
     * half of it (|dr_ij| <= |d_i| + |d_j| <= 2 max_k |d_k|). Rebuilding whenever the largest
     * single-atom displacement since the last build exceeds skin/2 = 5 Bohr is therefore
     * exact — the classic Verlet neighbour-list argument (L. Verlet, Phys. Rev. 159, 98
     * (1967)) — and costs one O(N) displacement scan per step. A step-count trigger like
     * the repulsion one was rejected on measurement: the build takes ~670 ms at 7320 atoms
     * (9.5 M pairs), about half of an MD step, while the displacement trigger fires a few
     * times per picosecond at most.
     *
     * The size-dependent RMSD formula of shouldUpdateHBXB() is deliberately NOT used (it is a
     * port-fidelity question of its own, see TODO.md); the trigger here is a true per-atom
     * maximum. When `dispersion_cutoff_bohr` shrinks the skin to zero (cutoff <= 50 Bohr)
     * the list falls back to the `nonbonded_rebuild_every` step count.
     *
     * The generator is reused, NOT recreated: the workspace holds a pointer into its dC6/dCN
     * matrix (setDC6DCNPtr()). Must run before prepareCNAndEEQ() so that the per-step Gaussian
     * weights, dC6/dCN and C6 refresh of that step already see the new pair list.
     *
     * @param ws External workspace to update as well (mirrors updateHBXBIfNeeded()).
     */
    void updateDispersionPairsIfNeeded(FFWorkspace* ws);

    /// True if updateDispersionPairsIfNeeded() rebuilt the list since the last check
    bool consumeDispersionPairsUpdate() { bool r = m_disp_pairs_updated; m_disp_pairs_updated = false; return r; }

    /**
     * @brief Recompute every stored D4 pair's C6 from the current-step coordination numbers.
     *
     * Claude Generated (Sep 2026): C6_ij = sum_ab W_i^a(CN_i) W_j^b(CN_j) C6ref_ab depends on
     * the geometry through the CN-Gaussian weights W (D4: Caldeweyher et al., J. Chem. Phys.
     * 150, 154122 (2019); GFN-FF CN-only weighting, gfnff_gdisp0.f90:405). The pair list
     * stored C6 at the SETUP geometry and nothing ever updated it, while the per-step
     * dispersion gradient already used dC6/dCN at the CURRENT CN (updateCNValuesForGradient()
     * -> dc6dcn). After any geometry change the energy and its own gradient therefore
     * described two different functions. Measured (Sep 2026, triose, 66 atoms): after 200 fs of
     * 800 K MD the MD dispersion energy was 2.83e-4 Eh (0.18 kcal/mol) away from a fresh
     * single point at the same geometry, and a geometry optimisation converged 3.0e-4 Eh away
     * from the single point at its own final geometry — in both cases the entire difference
     * sat in the dispersion term. Single points are unaffected: there the setup geometry is
     * the evaluation geometry.
     *
     * Reads the Gaussian weights / half-contraction the generator holds for this step, so it
     * must run right after updateCNValuesForGradient(). It inherits that function's
     * `d4_cn_cache_threshold` skip: when no CN changed by more than the threshold, the weights
     * (and therefore C6 and dC6/dCN) are left at the previous step, consistently for both.
     * Controlled by the `dispersion_c6_update` PARAM (default on).
     * @return true if the C6 of the installed D4 pair list were rewritten from the generator's
     *         current weights (stale geometry, non-empty list); false if nothing was done.
     */
    bool refreshDispersionC6();

    /// Claude Generated (Sep 2026): true if the stored D4 C6 were computed at a different
    /// geometry than the current one (and the C6 refresh is enabled). See refreshDispersionC6().
    bool dispersionC6Stale() const;

    /**
     * @brief Periodically rebuild the explicit Coulomb pair list (distance-truncated EEQ only).
     *
     * Claude Generated (Sep 2026): with `eeq_distance_cutoff > 0` the Coulomb term is evaluated
     * from an explicit pair list built ONCE (generateCoulombPairsNative(), cell list at the
     * cutoff) whose kernel cutoff equals its build radius — zero skin, so any pair that moves
     * inside the cutoff after setup is simply never evaluated. The default implicit path
     * (`coulomb_implicit`, and the GPU `gpu_coulomb_implicit`) enumerates every pair every step
     * and has no list to go stale, and the explicit list with no distance cutoff contains all
     * N(N-1)/2 pairs, so neither needs this. Same schedule as the repulsion list
     * (`nonbonded_rebuild_every`, default every step).
     *
     * @param ws External workspace to update as well (mirrors updateHBXBIfNeeded()).
     */
    void updateCoulombPairsIfNeeded(FFWorkspace* ws);

    /// True if updateCoulombPairsIfNeeded() rebuilt the list since the last check
    bool consumeCoulombPairsUpdate() { bool r = m_coul_pairs_updated; m_coul_pairs_updated = false; return r; }

    /**
     * @brief Largest single-atom displacement (3D norm, same unit as the inputs) between two
     * geometries of equal shape; +infinity if the shapes differ.
     *
     * Claude Generated (Sep 2026): the size-independent quantity a neighbour-list skin is
     * compared against (see updateDispersionPairsIfNeeded()).
     */
    static double maxAtomDisplacement(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b);

    // Claude Generated (Apr 2026): Timing accessors for GPU orchestrator
    double getParamGenTimeMs() const { return m_param_gen_time_ms; }
    double getTopologyTimeMs() const { return m_topology_time_ms; }

    /// Set external topology decision from GPU displacement check (Claude Generated March 2026).
    /// If set, needsFullTopologyUpdate() uses this value instead of CPU computation.
    void setExternalTopologyDecision(bool needs_full_update) const {
        m_external_topology_decision = needs_full_update;
    }

    /// Returns true if the last getCachedTopology() triggered a full topology recalculation.
    /// One-shot: resets to false after read. Used by GPU path to know when to update ref geometry.
    bool consumeFullTopologyUpdate() const {
        bool r = m_full_topology_recalculated;
        m_full_topology_recalculated = false;
        return r;
    }

    // === React topology mode (Claude Generated Aug 2026) ===
    // Event-driven reactive bond topology: the bond list is re-detected with a distance
    // hysteresis during MD and all bonded terms (bonds, angles, torsions, inversions)
    // plus the bonded/non-bonded repulsion partition are rebuilt when it changes.
    // Active only for topology_mode == "react". See docs/GFNFF_REACT_TOPOLOGY.md.

    /**
     * @brief Run the react-mode bond scan and rebuild all bonded terms if the bond set changed.
     *
     * No-op unless topology_mode == "react". Called at the start of Calculation() so
     * EEQ constraints, HB/XB detection and all force-field terms see the new topology
     * within the same step.
     */
    void updateReactiveTopologyIfNeeded();

    /// One-shot: true if updateReactiveTopologyIfNeeded() rebuilt the topology since the
    /// last call. Consumed by the GPU/HIP wrappers to trigger a device workspace rebuild.
    bool consumeReactRebuild() { bool r = m_react_rebuilt; m_react_rebuilt = false; return r; }

    /// Current authoritative react-mode bond set (canonical i<j pairs). Empty unless react mode.
    const std::vector<std::pair<int,int>>& reactiveBonds() const { return m_react_bonds; }

    /// Number of bonded-term rebuilds since initialisation (react mode).
    int reactiveRebuildCount() const { return m_react_rebuild_count; }

    /// Bond orders parallel to reactiveBonds(): 1 + round(Hueckel pi order), clamped to
    /// 1..3, refreshed at every rebuild (react mode). Cheap to read; never recomputes.
    const std::vector<int>& reactiveBondOrders() const { return m_react_bond_orders; }

    /// Topology mode in effect after parameter parsing: "auto", "constant" or "react".
    const std::string& topologyMode() const { return m_topology_mode; }

    /**
     * @brief One react-mode topology change event (Claude Generated Sep 2026).
     *
     * Recorded by detectReactiveBondChanges(); de_jump_eh is filled by the rebuild:
     * E(new topology with its EEQ charges) - E(old topology with its charges) at the
     * same geometry, i.e. the potential-energy discontinuity the dynamics experiences.
     * NaN when it could not be measured (no CN state yet, first energy call).
     */
    struct ReactEvent {
        long call = 0;                                    ///< energy-call counter at the scan
        std::vector<std::pair<int, int>> formed;          ///< canonical i<j pairs
        std::vector<std::pair<int, int>> broken;          ///< canonical i<j pairs (incl. exchange resolutions)
        double de_jump_eh = std::numeric_limits<double>::quiet_NaN();
    };

    /// Move out all events recorded since the last call (react mode; empty otherwise).
    /// A run that never consumes them (the CLI) keeps only the most recent
    /// kReactEventLimit; the log line of every event is written regardless.
    std::vector<ReactEvent> consumeReactEvents() { return std::exchange(m_react_events, {}); }
    static constexpr size_t kReactEventLimit = 4096;

    const std::vector<GFNFFHydrogenBond>& getLastHBonds() const { return m_last_hbonds; }
    const std::vector<GFNFFHalogenBond>& getLastXBonds() const { return m_last_xbonds; }

    /**
     * @brief Result of rebuildBondHBData() — per-bond HB metadata for GPU upload.
     *
     * Claude Generated (Apr 2026): GPU path needs bond nr_hb / hb_H_atom arrays and
     * the flat BondHBEntry list to be rebuilt after each dynamic HB re-detection.
     */
    struct BondHBRebuildResult {
        std::vector<BondHBEntry> bond_hb_data;   ///< Flat (A,H,B) entries for GPU HB-alpha pairs
        std::vector<int> bond_nr_hb;             ///< Per-bond HB count (size = nb)
        std::vector<int> bond_hb_H_atom;         ///< Per-bond H atom index, -1 if not an HB bond (size = nb)
    };

    /**
     * @brief Rebuild bond-HB cross-reference data from freshly re-detected HB list.
     *
     * Claude Generated (Apr 2026): Called by GFNFFGPUMethod after consumeHBXBUpdate()
     * to propagate new HB topology into the GPU bond SoA (nr_hb, hb_H_atom) and the
     * HB-alpha pair list.  Mirrors the cross-referencing logic in generateGFNFFParameterSet().
     *
     * @param hbonds  New HB list from getLastHBonds()
     * @param bonds   Static bond list from GFNFFParameterSet (m_gpu_params->bonds)
     * @return BondHBRebuildResult with updated bond_hb_data, per-bond nr_hb, per-bond hb_H_atom
     */
    /// Same as above from plain (A, H, B) arrays (device detection, which keeps no host struct list)
    /// Shared body of the two rebuilds: n entries, get(m, A, H, B) returns entry m
    BondHBRebuildResult rebuildBondHBDataCore(size_t n, const std::function<void(size_t, int&, int&, int&)>& get,
                                              const std::vector<Bond>& bonds) const;
    BondHBRebuildResult rebuildBondHBDataFromTriples(const std::vector<int>& hb_i, const std::vector<int>& hb_j,
                                                     const std::vector<int>& hb_k,
                                                     const std::vector<Bond>& bonds) const;
    BondHBRebuildResult rebuildBondHBData(const std::vector<GFNFFHydrogenBond>& hbonds,
                                           const std::vector<Bond>& bonds) const;

    // Getters for CN/EEQ results (valid after prepareCNAndEEQ)
    const Vector& getLastCN() const { return m_last_cn; }
    const Vector& getLastCharges() const { return m_charges; }

    // F-Q4 (Claude Generated): true iff the most recent InitialiseMolecule/Calculation
    // had the EEQ solver fall back to uniform/placeholder charges. The wrapper refuses
    // the result instead of returning a Coulomb energy built on wrong charges.
    // m_eeq_solve_failed is the sticky init-time flag (cleared per molecule); the live
    // solver flag catches per-step (opt/MD) fallbacks (cleared at each Calculation).
    bool eeqSolveFailed() const {
        return m_eeq_solve_failed || (m_eeq_solver && m_eeq_solver->lastSolveFailed());
    }
    // WP-G fix (May 2026): m_geometry_bohr is GeoGradMatrix (RowMajor). Returning
    // `const Matrix&` was creating a dangling reference to a temporary produced by
    // Eigen's implicit storage-order conversion — segfaulted in the GPU path
    // (GFNFFGPUComputationalMethod::calculateEnergy).
    const GeoGradMatrix& getGeometryBohr() const { return m_geometry_bohr; }

    /**
     * @brief True if m_charges are valid for the supplied geometry.
     * Mirrors the eeq_charges_current check in prepareCNAndEEQ
     * (gfnff_method.cpp:761-765). Used by the GPU path to short-circuit
     * a redundant Phase-2 EEQ when the geometry has not moved since the
     * last solve (typical for SP at the topology-time geometry).
     * Claude Generated (May 2026).
     */
    bool areEEQChargesCurrent(const Matrix& geom_bohr) const
    {
        return (m_charges.size() == m_atomcount)
            && (m_last_eeq_geometry.rows() == geom_bohr.rows())
            && (m_last_eeq_geometry.cols() == geom_bohr.cols())
            && (m_last_eeq_geometry == geom_bohr);
    }

    /**
     * @brief Store GPU-computed charges via memcpy (no Eigen heap alloc, no ForceField distribution).
     * Claude Generated (March 2026): Safe for GPU path where CUDA corrupts heap metadata.
     * Requires preAllocateForGPUPath() called first (m_charges pre-sized).
     * @param data  Pointer to N doubles
     * @param n     Number of atoms (must match m_charges.size())
     */
    void storeChargesFromGPU(const double* data, int n)
    {
        if (n == m_charges.size())
            std::memcpy(m_charges.data(), data, n * sizeof(double));
    }
    // Claude Generated (WP4, May 2026): CNDerivStore replaces std::vector<SpMatrix>
    const CNDerivStore& getLastCNDerivatives() const { return m_last_dcn; }

    /// WP-D Stage D (May 2026): fused CN + DCN single-pass result.
    struct CNAndDerivResult {
        Vector cn_values;                         ///< post-log CN (size N)
        Vector cn_raw;                            ///< pre-log erf-sum (size N)
        CNDerivStore dcn_store;                   ///< gradient pair-list + diagonal (dlogdcn applied)
        std::vector<std::vector<int>> neighbors;  ///< symmetric neighbor list (reusable by D4/EEQ)
    };
    const Vector& getLastCNF() const { return m_last_cnf; }
    FFWorkspace* getWorkspace() const { return m_workspace.get(); }
    /// GPU wrappers: keep the FULL parameter set for consumeCachedParameterSet() (default: bonded terms only).
    void setKeepFullParameterSet(bool keep) { m_keep_full_parameter_set = keep; }
    /// Shared CxxThreadPool (null before initializeForceField()).
    CxxThreadPool* threadPool() const;

    // Static-Mode (WP-S1, May 2026): expose frozen-state flags so the GPU method can
    // propagate them to FFWorkspaceGPU before launching kernels.
    bool staticCNFrozen() const { return m_static_cn && m_static_state_captured; }
    bool staticChargesFrozen() const { return m_static_charges && m_static_state_captured; }

    // WP-P1 (May 2026): last per-phase timings (ms) for MD diagnostics JSONL dump.
    const PrepTiming& getLastPrepTiming() const { return m_last_prep_timing; }

    // WP-D (May 2026): raw CN (pre log-squash) — populated by prepareCNAndEEQ.
    const Vector& getLastCNRaw() const { return m_last_cn_raw; }

    /// WP-P1 (May 2026): force per-phase chrono collection regardless of verbosity level.
    /// Set by SimpleMD when md_diagnostics_timing=true so the JSONL gets non-zero values.
    void setForcePhaseTiming(bool on) { m_force_phase_timing = on; }
    bool forcePhaseTiming() const { return m_force_phase_timing; }

    // Claude Generated (March 2026): Phase 2 GPU dc6dcn — expose D4 internals
    D4ParameterGenerator* getD4Generator() { return m_d4_generator.get(); }

    // Claude Generated (Apr 2026): P1a — Delegate CN-change threshold check to D4ParameterGenerator
    bool canSkipD4GaussianWeightsUpdate(const std::vector<double>& cn) const;
    void recordD4CNValues(const std::vector<double>& cn);

    /**
     * @brief Set external CN values (from GPU computation).
     * Claude Generated (March 2026): Phase 1 GPU CN migration.
     * Allows GPU-computed CN to replace CPU CN before EEQ calculation.
     * @param cn External CN values (size N)
     */
    void setLastCN(const Vector& cn) { m_last_cn = cn; }

    /**
     * @brief Pre-allocate per-step Eigen Vectors to correct size.
     * Claude Generated (March 2026): Must be called BEFORE CUDA init to avoid
     * heap corruption — CUDA allocations corrupt adjacent heap metadata,
     * making subsequent Eigen Vector resizes crash.
     * After this call, prepareCNAndEEQ() uses memcpy instead of Eigen assignment.
     * @param natoms Number of atoms
     */
    void preAllocateForGPUPath(int natoms)
    {
        m_last_cn  = Vector::Zero(natoms);
        m_last_cnf = Vector::Zero(natoms);
        if (m_charges.size() != natoms)
            m_charges = Vector::Zero(natoms);
        m_gpu_path_preallocated = true;
    }

    /// True if preAllocateForGPUPath() was called (enables memcpy path in prepareCNAndEEQ)
    bool isGPUPathPreallocated() const { return m_gpu_path_preallocated; }

    // ===== WP-FF-DistMatrix-Sharing (May 2026) =====
    // XTB gfnff_engrad.F90:175-195 pattern: compute packed-triangular sqrab/srab
    // ONCE per energy call, then every term/sub-routine indexes into them.
    // Curcuma had 95 .norm() calls in forcefieldthread.cpp + 53 in eeq_solver.cpp,
    // each recomputing the same distances. This API centralizes the work.

    /// Compute the shared packed-triangular distance array m_shared_srab (size N(N+1)/2) from
    /// m_geometry_bohr and the atoms-too-close flag (any r_ij < 1e-10 Bohr) for the EEQ solver.
    /// Only the off-diagonal slots (j < i) are written; nothing reads the diagonal slots.
    /// Thread-pool parallelized via threadPool().
    /// Called from GFNFF::Calculation() before prepareCNAndEEQ.
    void computeSharedDistances() const;

    /// Packed lower-triangular index for atom pair (i, j) with i != j.
    /// 0-based: idx = i*(i+1)/2 + j  for i > j.
    /// Caller responsible for i != j; diagonal not stored in packed form.
    static inline int triIdx(int i, int j) noexcept {
        if (j > i) { int t = i; i = j; j = t; }
        return i * (i + 1) / 2 + j;
    }

private:
    /**
     * @brief Initialize GFN-FF force field and generate parameters
     * @return true if successful
     */
    bool initializeForceField();

    /**
     * @brief Calculate topology and connectivity for GFN-FF
     * @return true if successful
     */
    bool calculateTopology();

    /**
     * @brief Retrieve cached topology information, computing it once if needed
     */
    const TopologyInfo& getCachedTopology() const;

    /**
     * @brief Retrieve cached bond list, computing it once if needed
     */
    const std::vector<std::pair<int,int>>& getCachedBondList() const;

    /**
     * @brief Calculate topological distances (bond counts) between all atom pairs using BFS
     * @param adjacency_list Per-atom neighbor connectivity
     * @return Sparse table of shortest path lengths up to 5 bonds (0=same, 1=bonded,
     *         2=1,3-pair, 3=1,4-pair, ...); every pair further apart or unconnected reads 999
     *
     * Claude Generated (Dec 24, 2025): Breadth-First Search for 1,3/1,4 topology factors.
     * Sparse storage since Sep 2026 (was a dense N x N matrix).
     */
    SparseTopoTable calculateTopologyDistances(const std::vector<std::vector<int>>& adjacency_list) const;

    /**
     * @brief Verbatim port of the reference's nbondmat (gfnff_ini2.f90:1280-1357).
     *
     * Produces topo%bpair: 1 for a direct bond as recorded in EITHER direction, 2 and 3
     * for pairs that reach each other SYMMETRICALLY within that many bonds, 5 for
     * everything else. The symmetry requirement (pairsbond's `dai .and. daj`,
     * gfnff_ini2.f90:1380) is what stops an eta bond — stored only on the metal's side —
     * from bridging a longer path, while the level-1 pass still records it as a bond.
     * Curcuma previously approximated this with a plain BFS plus an "eta-free" variant,
     * which got the two halves right separately but never together.
     *
     * @param nb Per-atom neighbour list; the reference passes topo%nb, i.e. the nbdum
     *           mixture that curcuma keeps in TopologyInfo::adjacency_list.
     * @return Sparse table holding the tags 1/2/3; every other pair reads 5, i == j reads 0.
     *         Sparse since Sep 2026 (was a dense N x N matrix plus a dense N x N
     *         membership matrix during construction).
     */
    SparseTopoTable computeBpairNbondmat(const std::vector<std::vector<int>>& nb) const;

    /**
     * @brief Detect molecular fragments (connected components)
     * @param adjacency_list Per-atom neighbor connectivity
     * @return Pair of (nfrag, fraglist)
     *
     * Claude Generated (Jan 31, 2026) - Ported from Fortran gfnff_helpers.f90:49-78 (mrecgff)
     */
    std::pair<int, std::vector<int>> detectMolecularFragments(const std::vector<std::vector<int>>& adjacency_list) const;

    /**
     * @brief Classify bond type according to GFN-FF topology rules
     * @param atom_i First atom index
     * @param atom_j Second atom index
     * @param hyb_i Hybridization of atom i (0=sp3, 1=sp, 2=sp2, 3=terminal, 5=hypervalent)
     * @param hyb_j Hybridization of atom j
     * @param is_metal_i True if atom i is a metal
     * @param is_metal_j True if atom j is a metal
     * @return Bond type (btyp): 1=single, 2=pi, 3=sp/linear, 4=hypervalent, 5=metal, 6=eta, 7=TM-TM
     *
     * Claude Generated (Jan 2, 2026): Ported from Fortran gfnff_ini.f90:1131-1148
     * Used for extra torsion filtering (btyp < 5 excludes metal bonds)
     */
    int classifyBondType(int atom_i, int atom_j, int hyb_i, int hyb_j,
                         bool is_metal_i, bool is_metal_j) const;

    /**
     * @brief Validate molecular structure for GFN-FF
     * @return true if molecule is valid
     */
    bool validateMolecule() const;

    /**
     * @brief Generate GFN-FF torsion parameters from topology
     *
     * Claude Generated (2025): Ported from Grimme Lab GFN-FF (Spicher & Grimme 2020)
     * Reference: external/gfnff/src/gfnff_engrad.F90:1041-1122 (egtors subroutine)
     *
     * Implements proper and improper torsion potentials for molecular rotations.
     * See docs/theory/GFNFF_TORSION_THEORY.md for scientific background.
     *
     * @return JSON array of torsion parameters
     */
    json generateGFNFFTorsions() const;

    /**
     * @brief Generate GFN-FF triple bond torsion (sTors_eg) parameters
     *
     * Claude Generated (March 2026): specialized torsion for sp-sp systems
     * Reference: external/gfnff/src/gfnff_ini.f90:2308 (specialTorsList)
     *
     * @return JSON array of triple bond torsion parameters
     */
    json generateGFNFFSTorsions() const;

    /**
     * @brief Generate GFN-FF inversion/out-of-plane parameters
     *
     * Claude Generated (2025): Inversion term implementation
     * Reference: external/gfnff/src/gfnff_helpers.f90:427-510 (omega, domegadr)
     *
     * Implements out-of-plane bending potentials for sp² centers (planarity constraints).
     * See docs/theory/GFNFF_INVERSION_THEORY.md for scientific background.
     *
     * @return JSON array of inversion parameters
     */
    json generateGFNFFInversions() const;

    // =================================================================================
    // NATIVE STRUCT GENERATORS (March 2026 — Phase 2 architecture cleanup)
    // Return std::vector<StructType> directly, bypassing JSON serialization.
    // The JSON generators above are kept as wrappers for file cache compatibility.
    // =================================================================================

    /// Generate bond parameters as native Bond structs
    std::vector<Bond> generateBondsNative(const TopologyInfo& topo_info) const;

    /// Generate angle parameters as native Angle structs
    std::vector<Angle> generateAnglesNative(const TopologyInfo& topo_info) const;

    /// Generate torsion parameters as native Dihedral structs (dihedrals + extra_dihedrals)
    std::pair<std::vector<Dihedral>, std::vector<Dihedral>> generateTorsionsNative() const;

    /// Generate inversion parameters as native Inversion structs
    std::vector<Inversion> generateInversionsNative() const;

    /// Generate triple bond torsion parameters as native GFNFFSTorsion structs
    std::vector<GFNFFSTorsion> generateSTorsionsNative() const;

    /// Generate Coulomb pair parameters as native GFNFFCoulomb structs
    std::vector<GFNFFCoulomb> generateCoulombPairsNative() const;

    /// Per-atom EEQ Coulomb self-energy inputs (chi_base/gam/alp/cnf/chi_static),
    /// independent of the pair list above — see generateCoulombSelfEnergyNative().
    struct CoulombSelfEnergy {
        Eigen::VectorXd chi_base, gam, alp, cnf, chi_static;
    };

    /// Generate per-atom Coulomb self-energy parameters (Claude Generated Sep 2026).
    /// Mirrors the per-atom half of generateCoulombPairsNative()'s fillPair(), but
    /// runs unconditionally for every atom (no pairing), matching the Fortran
    /// reference (gfnff_engrad.F90:1378-1389: the self-energy statement executes
    /// for every atom i regardless of whether the inner j<i pairwise loop has any
    /// iterations). Needed so a single isolated atom — where the pair list is
    /// structurally empty — still gets a nonzero EEQ self-energy.
    CoulombSelfEnergy generateCoulombSelfEnergyNative() const;

    /// Generate repulsion pair parameters as native GFNFFRepulsion structs (bonded + nonbonded)
    /// @param reuse_nonbonded optional storage of the previous non-bonded list (Oct 2026): its memory is taken over
    ///        and overwritten, so a rebuild during MD does not fault in a fresh 1.5 M-entry vector. Same list as without.
    std::pair<std::vector<GFNFFRepulsion>, std::vector<GFNFFRepulsion>> generateRepulsionPairsNative(
        std::vector<GFNFFRepulsion>* reuse_nonbonded = nullptr) const;

    /// Generate dispersion pair parameters as native GFNFFDispersion structs + ATM triples + method name
    std::tuple<std::vector<GFNFFDispersion>, std::vector<ATMTriple>, std::string> generateDispersionPairsNative() const;

    /// Claude Generated (Sep 2026): user `dispersion_cutoff_bohr` (0 = none), top-level or gfnff scope
    double dispersionCutoffBohr() const;
    /// Claude Generated (Sep 2026): the WP-Disp distance filter on a D4 pair list (no-op without a cutoff)
    void applyDispersionCutoff(std::vector<GFNFFDispersion>& dispersions, bool report) const;
    /// Claude Generated (Sep 2026): build radius minus evaluation radius of the D4 pair list (Bohr)
    double dispersionSkinBohr() const;
    /// Claude Generated (Sep 2026): Verlet skin of the repulsion / explicit-Coulomb lists (nonbonded_skin_bohr)
    double nonbondedSkinBohr() const;
    /// Claude Generated (Sep 2026): electrostatics pair cutoff in Bohr; coulomb_r_cut <= 0 (default) = none, as the reference
    double coulombRCutBohr() const;
    /// Claude Generated (Sep 2026): true if the repulsion list is cell-list built (distance-filtered),
    /// false if the O(N^2) build stores every pair (N < nb_cell_list_min_atoms) and so cannot go stale
    bool repulsionListIsDistanceFiltered() const;

    /// Detect hydrogen bonds as native GFNFFHydrogenBond structs
    /// @param storage  Optional previous list (Claude Generated, Oct 2026): taken over, its entries are reused as
    ///                 the slots of the result in the parallel cell-list path (same result as without it)
    std::vector<GFNFFHydrogenBond> detectHydrogenBondsNative(const Vector& charges,
                                                             std::vector<GFNFFHydrogenBond>* storage = nullptr) const;

    /// Detect halogen bonds as native GFNFFHalogenBond structs
    std::vector<GFNFFHalogenBond> detectHalogenBondsNative(const Vector& charges) const;

    // Sep 2026 (docs/GFNFF_PERFORMANCE_LEVERS.md lever #1): cheap, EXACT pre-struct HB/XB
    // strength estimates, used by detectHydrogenBondsNative/detectHalogenBondsNative to skip
    // GFNFFHydrogenBond/GFNFFHalogenBond allocation for candidates whose |E| the energy kernel
    // (ff_workspace_gfnff.cpp calcHydrogenBonds/calcHalogenBonds) would compute as negligible.
    // Each reuses the exact same GFNFFParameters damping primitives as that kernel — this is a
    // reordering of the real formula, not an approximation.
    // ONLY case 1 HB (below) and XB are pruned this way. Case 2/3/4 HB are intentionally
    // NEVER pruned by |E_HB|: their acceptor B, when N/O, feeds bond_hb_data / hb_cn_H
    // (ff_workspace_gfnff.cpp computeHBCoordinationNumbers), a purely GEOMETRIC erf-based
    // count that rescales the donor-H BOND term (egbond_hb) and does not correlate with
    // |E_HB| — an earlier attempt at pruning case 2/4 by |E_HB| shifted the bond term by
    // ~0.18 kcal/mol on a 66-atom test system despite the HB term itself moving by <1e-8 Eh.
    // Case 1 has no such coupling: a case-1 H is by construction not bonded to either
    // flanking atom, so it can never match the bond_hb_data lookup key.
    double estimateHBStrengthCase1(int A, int H, int B,
                                    double basicity_A, double basicity_B,
                                    double acidity_A, double acidity_B,
                                    double q_H, double q_A, double q_B) const;
    double estimateXBStrength(int A, int X, int B,
                               double acidity_X, double q_X, double q_B) const;

    /// Generate BATM triple parameters as native GFNFFBatmTriple structs
    std::vector<GFNFFBatmTriple> generateBatmTriplesNative(const TopologyInfo& topo_info) const;

    // Phase 4.2: GFN-FF pairwise non-bonded parameter generation (Claude Generated 2025)

    /**
     * @brief Extract D3/D4 configuration from main GFN-FF config
     * @param method Dispersion method: "d3" or "d4"
     * @return ConfigManager for D3/D4 parameter generator
     *
     * Claude Generated (December 2025): Configuration helper
     * Extracts relevant parameters (s6, s8, a1, a2) from main config
     * and creates ConfigManager for D3ParameterGenerator or D4ParameterGenerator.
     */
    ConfigManager extractDispersionConfig(const std::string& method) const;

    /**
     * @brief Get covalent radius for element
     * @param atomic_number Element atomic number
     * @return Covalent radius in Angstrom
     */
    double getCovalentRadius(int atomic_number) const;

    // GFN-FF parameter structures
    struct GFNFFBondParams {
        double force_constant;        // k_b in Fortran (energy scale)
        double equilibrium_distance;  // r₀ reference bond length
        double alpha;                 // α exponential decay parameter (was: anharmonic_factor)
        double rabshift;              // Claude Generated (Dec 2025): vbond(1) = gen%rabshift + shift
        double fqq;                   // Claude Generated (Jan 7, 2026): charge-dependent force constant factor

        // Claude Generated (Jan 18, 2026): Dynamic r0 calculation parameters
        // Reference: Fortran gfnff_rab.f90:147-153 - r0 recalculated at each Calculate()
        // Formula: r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff
        int z_i = 0, z_j = 0;           // Atomic numbers for parameter lookup
        double r0_base_i = 0.0;          // r0_gfnff[z_i-1] (Bohr)
        double r0_base_j = 0.0;          // r0_gfnff[z_j-1] (Bohr)
        double cnfak_i = 0.0;            // cnfak_gfnff[z_i-1]
        double cnfak_j = 0.0;            // cnfak_gfnff[z_j-1]
        double ff = 1.0;                 // EN-correction: 1 - k1*|ΔEN| - k2*ΔEN²
    };

    struct GFNFFAngleParams {
        double force_constant;     // k_ijk in Fortran
        double equilibrium_angle;  // θ₀ reference angle
        // Phase 1.3: Removed c0,c1,c2 Fourier coefficients (were dummy values)
        // GFN-FF uses simple angle bending, not Fourier expansion
    };

    /**
     * @brief GFN-FF torsion parameters
     *
     * Claude Generated (2025): Based on GFN-FF method (Spicher & Grimme 2020)
     *
     * Torsion potential: E = V/2 * [1 - cos(n*(φ - φ₀))] * D(r_ij, r_jk, r_kl)
     *
     * Scientific background:
     * - Describes rotation around central bond j-k in i-j-k-l sequence
     * - Periodicity n determines symmetry (1, 2, or 3)
     * - Barrier height V controls rotation difficulty
     * - Distance damping D couples stretching with rotation
     *
     * Reference: docs/theory/GFNFF_TORSION_THEORY.md
     */
    struct GFNFFTorsionParams {
        double barrier_height;     ///< V_n: Energy barrier in Hartree (Spicher & Grimme Eq. 8)
        int periodicity;            ///< n: Rotational symmetry (1, 2, or 3)
        double phase_shift;         ///< φ₀: Reference angle in radians
        bool is_improper;          ///< True for out-of-plane/improper torsions
        double fij_corrected;      ///< Central bond factor after all corrections (H-count, amide, alphaCO, CN, N-sat, hypervalent)
        double fkl_corrected;      ///< Outer atom factor after all corrections (N-reduction, hypervalent, CN)
        double fqq;                ///< Charge correction factor: 1 + |qa_j*qa_k| * qfacTOR
    };

    /**
     * @brief GFN-FF inversion/out-of-plane parameters
     *
     * Claude Generated (2025): Based on GFN-FF method (Spicher & Grimme 2020)
     *
     * Inversion potential: E = V * [cos(ω) - cos(ω₀)]² * D(r_ij, r_jk, r_jl)
     *
     * Scientific background:
     * - Describes out-of-plane bending for atom i relative to plane j-k-l
     * - Enforces planarity at sp² centers (aromatics, C=C, C=O)
     * - ω = out-of-plane angle ∈ [-π/2, +π/2]
     * - Double-well potential allows ±ω₀ equivalence
     *
     * Reference: docs/theory/GFNFF_INVERSION_THEORY.md
     */
    struct GFNFFInversionParams {
        double barrier_height;     ///< V: Energy barrier in kcal/mol
        double reference_angle;     ///< ω₀: Reference angle in radians (usually 0 for planar)
        int potential_type;        ///< 0: double-well [cos(ω)-cos(ω₀)]², 1: single-well [1-cos(ω)]
    };

    /**
     * @brief EEQ (Electronegativity Equalization) parameters
     *
     * Claude Generated (2025): Phase 3 EEQ charge calculation
     *
     * Parameters from gfnff_param.f90 (angewChem2020 parameter set)
     * Reference: S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673
     *
     * Scientific background:
     * - EEQ method solves linear system A·q = b for atomic charges
     * - chi: atomic electronegativity (controls charge distribution)
     * - gam: chemical hardness (resistance to charge transfer)
     * - alp: damping parameter for Coulomb interaction
     * - cnf: coordination number correction factor
     */
    struct EEQParameters {
        double chi;  ///< Electronegativity (angewChem2020)
        double gam;  ///< Chemical hardness (angewChem2020)
        double alp;  ///< Damping parameter (angewChem2020)
        double cnf;  ///< CN correction factor (angewChem2020)
        double xi_corr;  ///< Environment correction (topology-dependent, optional)
    };

    /**
     * @brief Get GFN-FF bond parameters for element pair
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param distance Current bond distance
     * @return GFN-FF bond parameters
     */
    /**
     * @brief Get GFN-FF bond parameters with full topology corrections (Phase 9)
     * @param atom1 First atom index
     * @param atom2 Second atom index
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param distance Current bond distance
     * @param topo Topology information (CN, hyb, charges, rings)
     * @return GFN-FF bond parameters with all corrections
     */
    GFNFFBondParams getGFNFFBondParameters(int atom1, int atom2, int z1, int z2,
                                            double distance, const TopologyInfo& topo) const;

    /**
     * @brief Get EEQ parameters for an element
     *
     * Claude Generated (2025): Phase 3 EEQ implementation
     * Reference: gfnff_param.f90 (chi/gam/alp/cnf_angewChem2020)
     *
     * Returns EEQ parameters from angewChem2020 parameter set:
     * - chi: atomic electronegativity
     * - gam: chemical hardness
     * - alp: damping parameter for erf(γ*r)/r Coulomb interaction
     * - cnf: coordination number correction factor
     *
     * @param atomic_number Element atomic number (1-86)
     * @return EEQ parameters structure
     */
    EEQParameters getEEQParameters(int atomic_number) const;

    /**
     * @brief Get GFN-FF angle parameters for angle triplet
     * @param atom_i First atom index (needed for hybridization lookup)
     * @param atom_j Center atom index (determines equilibrium angle via hybridization)
     * @param atom_k Third atom index
     * @param current_angle Current angle in radians (used for geometry-dependent overrides)
     * @param coord_numbers Pre-computed coordination numbers for all atoms (Claude Generated February 2026)
     * @return GFN-FF angle parameters with topology-based equilibrium angles
     *
     * Claude Generated (Nov 2025): Phase 2 implementation uses topology-aware parameters
     * including charge-dependent corrections (fqq), coordination number scaling (fn),
     * element-specific corrections (f2), and small-angle corrections (fbsmall).
     *
     * Claude Generated (February 2026): Added coord_numbers parameter to eliminate
     * 2,614 redundant CN calculations (26 seconds → 0.01 seconds, 2600× speedup!)
     */
    GFNFFAngleParams getGFNFFAngleParameters(int atom_i, int atom_j, int atom_k,
                                              double current_angle, const TopologyInfo& topo_info,
                                              const Vector& coord_numbers) const;

    /**
     * @brief Get GFN-FF torsion parameters for atom quartet
     *
     * Claude Generated (2025): Topology-aware parameter assignment
     * Reference: external/gfnff/src/gfnff_ini2.f90 (torsion setup)
     *
     * Assigns torsion parameters based on:
     * - Hybridization of central atoms j and k
     * - Ring membership (strain corrections)
     * - Conjugation status (planarity preferences)
     *
     * @param z_i Atomic number of first atom
     * @param z_j Atomic number of second atom (central bond)
     * @param z_k Atomic number of third atom (central bond)
     * @param z_l Atomic number of fourth atom
     * @param hyb_j Hybridization of atom j (1=sp, 2=sp2, 3=sp3)
     * @param hyb_k Hybridization of atom k
     * @return GFN-FF torsion parameters
     */
    GFNFFTorsionParams getGFNFFTorsionParameters(int z_i, int z_j, int z_k, int z_l,
                                                  int hyb_j, int hyb_k,
                                                  double qa_j = 0.0, double qa_k = 0.0,
                                                  double cn_i = 2.0, double cn_l = 2.0,
                                                  bool in_ring = false, int ring_size = 0,
                                                  int i_atom_idx = -1, int j_atom_idx = -1,
                                                  int k_atom_idx = -1, int l_atom_idx = -1,
                                                  int bond_type = 1) const;

    /**
     * @brief Calculate dihedral angle for four atoms
     *
     * Claude Generated (2025): Standard computational chemistry formula
     * Reference: Allen & Tildesley "Computer Simulation of Liquids" (1987)
     *
     * Computes signed dihedral angle φ ∈ [-π, π] between planes i-j-k and j-k-l.
     * Uses atan2 for proper sign handling (crucial for gradient calculation).
     *
     * Physical interpretation:
     * - φ = 0°: cis/eclipsed (atoms i and l on same side)
     * - φ = 180°: trans/anti (atoms i and l on opposite sides)
     *
     * @param i Index of first atom
     * @param j Index of second atom
     * @param k Index of third atom
     * @param l Index of fourth atom
     * @return Dihedral angle in radians [-π, π]
     */
    double calculateDihedralAngle(int i, int j, int k, int l) const;

    /**
     * @brief Calculate derivatives of dihedral angle w.r.t. atomic positions
     *
     * Claude Generated (2025): Analytical gradient for torsions
     * Reference: external/gfnff/src/gfnff_engrad.F90 (dphidr subroutine)
     *
     * Computes ∂φ/∂x_i, ∂φ/∂x_j, ∂φ/∂x_k, ∂φ/∂x_l for chain rule in gradient.
     * Critical for analytical force calculation in molecular dynamics.
     *
     * @param i Index of first atom
     * @param j Index of second atom
     * @param k Index of third atom
     * @param l Index of fourth atom
     * @param phi Current dihedral angle (from calculateDihedralAngle)
     * @param dda Output: ∂φ/∂x_i (3D vector)
     * @param ddb Output: ∂φ/∂x_j (3D vector)
     * @param ddc Output: ∂φ/∂x_k (3D vector)
     * @param ddd Output: ∂φ/∂x_l (3D vector)
     */
    void calculateDihedralGradient(int i, int j, int k, int l, double phi,
                                     Vector& dda, Vector& ddb,
                                     Vector& ddc, Vector& ddd) const;

    /**
     * @brief Calculate damping function for torsion potential
     *
     * Claude Generated (2025): Distance-dependent damping
     * Reference: external/gfnff/src/gfnff_engrad.F90 (gfnffdampt function)
     *
     * Damping function couples bond stretching with torsional motion:
     * D(r) = 1 / [1 + exp(-α*(r/r₀ - 1))]
     *
     * Physical meaning:
     * - Stretched bonds → weaker torsion barrier (easier rotation)
     * - Compressed bonds → stronger torsion barrier
     *
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param r_squared Squared distance r² in Bohr²
     * @param damp Output: Damping value D(r)
     * @param damp_deriv Output: Derivative ∂D/∂r
     */
    void calculateTorsionDamping(int z1, int z2, double r_squared,
                                  double& damp, double& damp_deriv) const;

    // =================================================================================
    // INVERSION/OUT-OF-PLANE HELPER FUNCTIONS (Phase 1.2)
    // =================================================================================

    /**
     * @brief Calculate out-of-plane angle (omega) for atom i relative to plane j-k-l
     *
     * Claude Generated (2025): Inversion angle calculation
     * Reference: external/gfnff/src/gfnff_helpers.f90:427-448 (omega function)
     *
     * Computes ω ∈ [-π/2, +π/2] measuring deviation from planarity:
     *   ω = arcsin(n · v̂)
     * where:
     *   n = (r_ij × r_jk) / |r_ij × r_jk|  (normal to plane i-j-k)
     *   v = r_il  (vector from i to l)
     *
     * Physical interpretation:
     * - ω = 0: atom i in plane j-k-l (planar, typical for sp²)
     * - ω = ±π/2: atom i perpendicular to plane (pyramidal)
     *
     * @param i Index of central atom (out-of-plane)
     * @param j Index of first plane atom
     * @param k Index of second plane atom
     * @param l Index of third plane atom
     * @return Out-of-plane angle in radians [-π/2, π/2]
     */
    double calculateOutOfPlaneAngle(int i, int j, int k, int l) const;

    /**
     * @brief Calculate derivatives of out-of-plane angle w.r.t. atomic positions
     *
     * Claude Generated (2025): Analytical gradient for inversions
     * Reference: external/gfnff/src/gfnff_helpers.f90:450-510 (domegadr subroutine)
     *
     * Computes ∂ω/∂x_i, ∂ω/∂x_j, ∂ω/∂x_k, ∂ω/∂x_l for chain rule in gradient.
     * Critical for analytical force calculation in geometry optimization.
     *
     * @param i Index of central atom
     * @param j Index of first plane atom
     * @param k Index of second plane atom
     * @param l Index of third plane atom
     * @param omega Current out-of-plane angle (from calculateOutOfPlaneAngle)
     * @param grad_i Output: ∂ω/∂x_i (3D vector)
     * @param grad_j Output: ∂ω/∂x_j (3D vector)
     * @param grad_k Output: ∂ω/∂x_k (3D vector)
     * @param grad_l Output: ∂ω/∂x_l (3D vector)
     */
    void calculateInversionGradient(int i, int j, int k, int l, double omega,
                                     Vector& grad_i, Vector& grad_j,
                                     Vector& grad_k, Vector& grad_l) const;

    /**
     * @brief Get GFN-FF inversion parameters for atom quartet
     *
     * Claude Generated (2025): Topology-aware inversion parameter assignment
     *
     * Assigns inversion parameters based on:
     * - Hybridization of central atom i (sp² → needs inversion)
     * - Element type (C, N, O, B different barriers)
     * - Pi-system membership (aromatics → higher barriers)
     *
     * @param z_i Atomic number of central atom (out-of-plane)
     * @param z_j Atomic number of plane atom j
     * @param z_k Atomic number of plane atom k
     * @param z_l Atomic number of plane atom l
     * @param hyb_i Hybridization of central atom i (1=sp, 2=sp², 3=sp³)
     * @return GFN-FF inversion parameters
     */
    GFNFFInversionParams getGFNFFInversionParameters(int z_i, int z_j, int z_k, int z_l,
                                                      int hyb_i) const;

    // =================================================================================
    // Advanced GFN-FF Parameter Generation (for future implementation)
    // =================================================================================

    /**
     * @brief Calculate coordination number derivatives for gradients
     * @param cn Coordination numbers
     * @param threshold Coordination number threshold (squared distance in Bohr²)
     * @return 3D tensor of CN derivatives (3 x natoms x natoms)
     */
    // Claude Generated (WP4, May 2026): returns CNDerivStore (pair-list + diag) instead of std::vector<SpMatrix>
    // Eliminates ~1000 ms triplet+setFromTriplets cost on mixture.xyz N=6200 (74 % of CN+EEQ phase per WP1).
    /// WP-D (May 2026): main implementation. When `cn_raw_in.size() == m_atomcount`
    /// the function skips the redundant N²-erf loop in step 1 and uses `cn_raw_in`
    /// directly for dlogdcn. Otherwise (empty vec) falls back to internal recompute.
    /// WP-D Stage C (May 2026): when `neighbors != nullptr`, the dcn step-3 loop iterates
    /// the pre-filtered neighbor list instead of the full N² triangle, eliminating O(N²)
    /// threshold checks. Pass nullptr to use the original N² fallback.
    CNDerivStore calculateCoordinationNumberDerivatives(const Vector& cn, const Vector& cn_raw_in,
                                                       double threshold = 1600.0,
                                                       CxxThreadPool* pool = nullptr,
                                                       int num_threads = 1,
                                                       const std::vector<std::vector<int>>* neighbors = nullptr) const;

    /// Backward-compat wrapper — legacy callers without cn_raw access.
    inline CNDerivStore calculateCoordinationNumberDerivatives(
        const Vector& cn, double threshold = 1600.0,
        CxxThreadPool* pool = nullptr, int num_threads = 1) const
    {
        return calculateCoordinationNumberDerivatives(cn, Vector{}, threshold, pool, num_threads);
    }

    /// WP-D Stage D (May 2026): fused CN + DCN in a single O(N²) pair pass.
    /// Eliminates the redundant geometry-read pass that calculateCoordinationNumberDerivatives
    /// performs after calculateGFNFFCNWithNeighbors. Only called when GFNFF_CN_DCN_FUSION is
    /// defined and gradient=true and cn_cutoff_bohr > 0 and !gpu_only and !reuse_cn.
    CNAndDerivResult computeCNAndDerivativesFused(
        double cn_cutoff_bohr,
        CxxThreadPool* pool,
        int num_threads) const;

    /**
     * @brief Determine hybridization states for all atoms
     * @return Vector of hybridization states (1=sp, 2=sp2, 3=sp3, 4=sp3d, 5=sp3d2)
     */
    /**
     * @brief Determine hybridization states (PHASE 2 OPTIMIZED)
     * @param adjacency_list Pre-computed bond connectivity (eliminates O(N²) loop)
     * @return Hybridization states (1=sp, 2=sp2, 3=sp3, etc.)
     */
    std::vector<int> determineHybridization(const std::vector<std::vector<int>>& adjacency_list) const;

    /**
     * @brief Detect pi-systems and conjugated fragments (PHASE 2 OPTIMIZED)
     * @param hyb Hybridization states
     * @param adjacency_list Pre-computed bond connectivity (eliminates O(N²) loop)
     * @param nb_full Full (unfiltered) neighbour list, for the N/S pi-veto below
     * @return Vector mapping atoms to pi-fragment IDs (0 = no pi-system)
     */
    std::vector<int> detectPiSystems(const std::vector<int>& hyb,
                                     const std::vector<std::vector<int>>& adjacency_list,
                                     const std::vector<std::vector<int>>& nb_full) const;

    /**
     * @brief Find smallest ring size for each atom and enumerate all rings
     * @param adjacency_list Pre-computed bond connectivity (eliminates O(N²) loop)
     * @param topo_info TopologyInfo to populate with rings and atom_to_rings
     * @return Vector of smallest ring sizes (0 = not in ring)
     *
     * Claude Generated (Feb 9, 2026): Also populates topo_info.rings and topo_info.atom_to_rings
     */
    std::vector<int> findSmallestRings(const std::vector<std::vector<int>>& adjacency_list,
                                       TopologyInfo& topo_info) const;

    /**
     * @brief Check if two atoms are in the same ring
     * @param i First atom index
     * @param j Second atom index
     * @param ring_size Output: size of the smallest ring they share (0 if not in same ring)
     * @return true if atoms are in the same ring
     *
     * Claude Generated (Feb 9, 2026): Rewritten to use pre-computed ring membership
     * instead of broken ad-hoc path-finding. O(1) lookup via atom_to_rings.
     */
    bool areAtomsInSameRing(int i, int j, int& ring_size) const;

    /**
     * @brief Find smallest ring containing all four atoms of a torsion quartet
     * @param i First atom (terminal)
     * @param j Second atom (central bond)
     * @param k Third atom (central bond)
     * @param l Fourth atom (terminal)
     * @return Size of smallest ring containing all four atoms, or 0 if none
     *
     * Claude Generated (Feb 9, 2026): Equivalent to Fortran ringstors/rings4
     * Reference: gfnff_ini.f90:1846 (rings4 = ringstors(ii,jj,kk,ll,...))
     */
    int smallestRingContainingAll(int i, int j, int k, int l) const;

    /**
     * @brief Find smallest ring containing all three atoms of an angle bend
     * @param i Center atom
     * @param j First neighbor
     * @param k Second neighbor
     * @return Size of smallest ring containing all three atoms, or 0 if none
     *
     * Claude Generated (Feb 11, 2026): Equivalent to Fortran ringsbend(i,j,k,...)
     * Reference: gfnff_ini2.f90:503-543
     */
    int smallestRingContainingBend(int i, int j, int k) const;

    /**
     * @brief Find largest ring containing all four atoms of a torsion quartet
     * @param i First atom (terminal)
     * @param j Second atom (central bond)
     * @param k Third atom (central bond)
     * @param l Fourth atom (terminal)
     * @return Size of largest ring containing all four atoms, or 0 if none
     *
     * Claude Generated (Feb 9, 2026): For Fortran ringl == rings4 check
     */
    int largestRingContainingAll(int i, int j, int k, int l) const;

    /**
     * @brief Calculate topology-dependent electronegativity corrections (dxi)
     *
     * Claude Generated (2025): Element-specific topology corrections for EEQ
     * Reference: gfnff_ini.f90:247-297 (dxicalc subroutine)
     *
     * Implements the missing dxi topology corrections that significantly
     * improve EEQ charge accuracy, especially for heteroatoms.
     *
     * dxi corrections include:
     * - Boron hydrogen-dependent corrections
     * - Carbene carbon corrections
     * - Oxygen nitro group and water corrections
     * - Group 6 (S, Se, etc.) overcoordination corrections
     * - Group 7 (Cl, Br, I) metal-dependent corrections
     *
     * @param atoms Atomic numbers
     * @param coordination_numbers Coordination numbers
     * @param neighbor_list Bond connectivity information
     * @param topo_info Topology information (pi-systems, etc.)
     * @return Vector of dxi corrections for each atom
     */
    Vector calculateDXI(const std::vector<int>& atoms,
                        const Vector& coordination_numbers,
                        const std::vector<std::vector<int>>& neighbor_list,
                        const TopologyInfo& topo_info) const;

    /**
     * @brief Calculate EEQ charges using extended electronegativity equalization
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param rings Ring information
     * @return Vector of EEQ charges
     */
    Vector calculateEEQCharges(const Vector& cn, const std::vector<int>& hyb, const std::vector<int>& rings) const;
    /// EEQ electrostatic energy for a charge set (fragment charge-assignment trial).
    double calculateEEQEnergy(const Vector& charges, const Vector& cn) const;

    /**
     * @brief Calculate dgam (charge-dependent hardness) corrections
     *
     * Claude Generated (December 2025, Session 6): Extracted from calculateEEQCharges()
     * Reference: external/gfnff/src/gfnff_ini.f90:677-688
     *
     * Calculates charge-dependent gamma corrections that refine the EEQ hardness matrix
     * based on computed atomic charges and element type.
     *
     * @param qa_charges Base EEQ charges (from Phase 3.3 of calculateEEQCharges)
     * @param hybridization Hybridization state per atom (1=sp, 2=sp2, 3=sp3)
     * @param ring_sizes Smallest ring size per atom (0 if not in ring)
     * @return dgam corrections: Delta-gamma values to add to hardness matrix diagonal
     */
    Vector calculateDgam(const Vector& qa_charges,
                        const std::vector<int>& hybridization,
                        const std::vector<int>& ring_sizes) const;

    /**
     * @brief Detect eta(η)-coordinated atoms (metal-alkene/alkyne/Cp side-on bonding)
     *
     * Claude Generated (July 2026). Faithful port of the Fortran etacoord logic
     * (external/gfnff/src/gfnff_ini2.f90:170-198). Sets itag[i] = -1 for genuine
     * η-coordinated carbons (only ati<=10, in practice C), -1 otherwise 0. Used by
     * the angle-bending feta metal correction so it fires ONLY for real η ligands
     * (metal-alkene/alkyne/cyclopentadienyl), NOT σ-bonded π ligands like CO.
     *
     * @param neighbor_lists Full bonded adjacency (== Fortran nbf, includes metals)
     * @return itag vector (size m_atomcount): -1 if η-coordinated, else 0
     */
    /// @param neighbor_lists nbf, the full list (getnb icase=1)
    /// @param nbm            the metal-filtered list (getnb icase=3) — the reference's
    ///                       nbm(20,i) is the SIZE OF THAT LIST, not "nbf minus metals"
    std::vector<int> computeEtaCoordination(const std::vector<std::vector<int>>& neighbor_lists,
                                            const std::vector<std::vector<int>>& nbm) const;

    /**
     * @brief Estimate per-atom "metallic character" mchar
     *
     * Claude Generated (July 2026). Port of external/gfnff/src/gfnff_ini.f90:243-250:
     *   mchar(i) = exp(-0.005 * en(Z_i)^8) * dum2 / (cn(i) + 1)
     * where dum2 = sum_j || d logCN_i / d R_j || over the GFN-FF logistic CN
     * (gfnff_cn.f90:gfnff_dlogcoord). Used only to gate the metal filter of the
     * metal-reduced neighbour list nbm (getnb icase=3: mchar > 0.25 => drop).
     *
     * The exp(-0.005*en^8) factor is extremely stiff (en=2.0 -> 0.28, en=3.0 -> ~0),
     * so it acts as a near-hard electronegativity switch around en ~ 2.4.
     *
     * NOTE: uses GFNFFParameters::gfnff_en (Fortran param%en), NOT the rab_en table.
     *
     * @return mchar per atom (size m_atomcount)
     */
    std::vector<double> computeMetallicCharacter() const;

    /**
     * @brief Build the Fortran four-list neighbour set (nbf / topo%nb / nbm / nbdum)
     *
     * Claude Generated (July 2026). Port of external/gfnff/src/gfnff_ini2.f90:128-130
     * and 197-202. Fortran derives three lists from the same distance criterion and
     * then assigns hybridization from a per-atom mixture:
     *   nbf   (icase=1) full connectivity, no filtering
     *   nb_hc (icase=2) drops ALL bonds of highly-coordinated atoms
     *                   (hc_crit = 4 if group <= 2 - which includes every transition
     *                    metal, since periodic_group is negative for the d-block - else 6)
     *   nbm   (icase=3) drops metals (mchar > 0.25 or metal_type > 0) and heavy atoms
     *                   above their normal CN (nbf > normcn and Z > 10)
     *   nbdum          per-atom mixture: nbm for eta-coordinated atoms, else nbf
     *
     * Fills topo.nb_full / nb_hc / nb_nometal / metallic_character / itag. The mixture
     * is what Fortran finally stores as topo%nb (gfnff_ini2.f90:335).
     *
     * @param topo Topology to populate (nb_full etc. are overwritten)
     * @param[out] nbdum The per-atom eta-aware mixture (Fortran's final topo%nb)
     */
    void buildNeighborListSet(GFNFFTopology& topo, std::vector<std::vector<int>>& nbdum) const;

    /**
     * @brief CN of the nearest non-metal neighbour of an atom
     *
     * Claude Generated (July 2026). Port of gfnff_ini2.f90:431-450 (nn_nearest_noM).
     * Fortran calls this with topo%nb, i.e. the HC-filtered list nb_hc.
     *
     * @param ii Atom index
     * @param nb Neighbour list to search (pass nb_hc for Fortran parity)
     * @param distance_matrix Interatomic distances
     * @return CN of the closest non-metal neighbour, or 0 if there is none
     */
    int nnNearestNoM(int ii, const std::vector<std::vector<int>>& nb,
                     const Eigen::MatrixXd& distance_matrix) const;

    /**
     * @brief Fortran-faithful hybridization assignment
     *
     * Claude Generated (July 2026). Transcription of gfnff_ini2.f90:211-333, replacing
     * the geometry-first heuristic in determineHybridization(). Reads the four-list
     * neighbour set: nb20i from the nbdum mixture, nbdiff/nbmdiff from nbf vs nb_hc/nbm,
     * and indexes the NO2/B-N/N-SO2 and CO->sp rules into nb_hc (which is what Fortran's
     * topo%nb still holds at this point in the initialization).
     *
     * Also WRITES itag: sets +1 for carbenes and NO2 nitrogen (consumed by Hueckel/HB),
     * on top of the -1 eta tags already set by computeEtaCoordination().
     *
     * @param topo Topology holding nb_full / nb_hc / nb_nometal / distance_matrix
     * @param nbdum The eta-aware per-atom mixture
     * @param itag In/out tag vector (-1 eta in, +1 carbene/NO2 out)
     * @return Hybridization per atom (0 none/octahedral, 1 sp, 2 sp2, 3 sp3, 5 hypervalent)
     */
    std::vector<int> determineHybridizationFortran(const GFNFFTopology& topo,
                                                   const std::vector<std::vector<int>>& nbdum,
                                                   const Eigen::MatrixXd& distance_matrix,
                                                   std::vector<int>& itag) const;

    /**
     * @brief Calculate simplified π-bond orders for all atom pairs
     *
     * Claude Generated (January 10, 2026) - Phase 2C: π-bond order approximation
     *
     * Claude Generated (January 14, 2026) - Updated for Phase 1: Full Hückel implementation
     *
     * Two modes available (controlled by m_use_full_huckel):
     *
     * **Full Hückel mode (default, m_use_full_huckel=true)**:
     * Uses iterative self-consistent Hückel method from gfnff_ini.f90:928-1062.
     * - P-dependent off-diagonal coupling (prevents over-delocalization)
     * - Charge-dependent diagonal elements
     * - Fermi smearing at 4000K for biradical handling
     * - Exact π-bond orders from density matrix elements
     *
     * **Simplified mode (m_use_full_huckel=false)**:
     * Approximation based on hybridization and bond types:
     * - Single bonds (sp3-sp3): pbo = 0.0
     * - π bonds (sp2-sp2, sp2-sp): pbo = 0.5-1.0
     * - sp bonds (sp-sp): pbo = 1.5
     *
     * @param bond_list Vector of bonded atom pairs
     * @param hybridization Hybridization state per atom
     * @param pi_fragments Pi-system fragment IDs
     * @param charges EEQ atomic charges (needed for full Hückel)
     * @param geometry_bohr N×3 geometry matrix in Bohr (P2a: replaces distance matrix)
     * @return Vector of π-bond orders in triangular format (access via lin(i,j))
     */
    std::vector<double> calculatePiBondOrders(
        const std::vector<std::pair<int,int>>& bond_list,
        const std::vector<int>& hybridization,
        const std::vector<int>& pi_fragments,
        const std::vector<double>& charges = {},
        const Eigen::MatrixXd& geometry_bohr = Eigen::MatrixXd(),
        const std::vector<int>& pi_system_charge = {},
        const std::vector<int>& itag = {},
        std::vector<int>* pi_atoms_final = nullptr) const;

    // Advanced parameter structures (EEQParameters already defined above at line 298)
    // TopologyInfo now defined at line 51 (public section) for use in function signatures

    /**
     * @brief Get EEQ parameters for specific atom with environment corrections
     * @param atom_idx Atom index
     * @param topo_info Topology information
     * @return EEQ parameters for this atom
     */
    EEQParameters getEEQParameters(int atom_idx, const TopologyInfo& topo_info) const;

    // =================================================================================
    // Two-Phase EEQ System Methods (Session 5, December 2025)
    // =================================================================================

    /**
     * @brief Phase 1: Calculate topology charges using base EEQ parameters
     *
     * Solves EEQ using ONLY base parameters without any corrections:
     * - chi = -chi_base (NO dxi corrections)
     * - gamma = gam_base (NO dgam corrections)
     * - alpha = alp_base^2 (NO dalpha corrections)
     *
     * These topology charges (qa) are then used in Phase 2 to calculate
     * the correction terms (dxi, dgam, dalpha).
     *
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param rings Ring information
     * @return Topology charges (qa) from base EEQ
     *
     * Reference: Fortran gfnff_ini.f90:405-421
     */
    Vector calculateTopologyCharges(const Vector& cn, const std::vector<int>& hyb,
                                     const std::vector<int>& rings) const;

    /**
     * @brief Calculate dxi (electronegativity) corrections
     *
     * Applies the full cascade logic from Fortran gfnff_ini.f90:361-403
     * with 30+ lines of element-specific, group-specific, and neighbor-dependent corrections.
     *
     * Corrections include:
     * - Boron: +nh*0.015 (hydrogen neighbors)
     * - Carbon: carbene (-0.15), free CO (+0.15)
     * - Oxygen: nitro (+0.05), water (-0.02), overcoordination (+nn*0.005)
     * - Group 6: overcoordination correction
     * - Group 7: polyvalent halogen corrections
     *
     * @param topo Topology information with neighbor lists and functional groups
     * @param qa_charges Topology charges from Phase 1
     * @return Electronegativity corrections per atom
     *
     * Reference: Fortran gfnff_ini.f90:361-403
     */
    Vector calculateDxi(const TopologyInfo& topo, const Vector& qa_charges) const;

    /**
     * @brief Calculate dalpha (polarizability) corrections
     *
     * Applies charge-dependent polarizability corrections:
     * alpeeq = (alp_base + ff * qa)^2
     *
     * where ff depends on element and group:
     * - C: +0.09, N: -0.21
     * - Group 6: -0.03, Group 7: +0.50
     * - Main-group metals: +0.3, Transition metals: -0.1
     *
     * @param qa_charges Topology charges from Phase 1
     * @return Polarizability corrections per atom (NOT yet squared)
     *
     * Reference: Fortran gfnff_ini.f90:694-707
     */
    Vector calculateDalpha(const Vector& qa_charges) const;

    /**
     * @brief Phase 2: Calculate final charges with all corrections applied
     *
     * Solves second EEQ with corrected parameters:
     * - chi = -chi_base + dxi [+ amide correction]
     * - gamma = gam_base + dgam
     * - alpha = (alp_base + dalpha)^2
     *
     * Also applies special amide hydrogen correction: chi -= 0.02
     *
     * @param topo Topology information with all corrections calculated
     * @return Final EEQ charges (q) to be used for force field calculations
     *
     * Reference: Fortran gfnff_ini.f90:694-707
     */
    Vector calculateFinalCharges(const TopologyInfo& topo) const;

public:
    // =================================================================================
    // Energy Component Access (Claude Generated November 2025)
    // =================================================================================

    /**
     * @brief Get bond energy component
     * @return Bond stretching energy or 0 if not calculated
     */
    double BondEnergy() const;

    /**
     * @brief Get angle energy component
     * @return Angle bending energy or 0 if not calculated
     */
    double AngleEnergy() const;

    /**
     * @brief Get dihedral energy component
     * @return Dihedral torsion energy or 0 if not calculated
     */
    double DihedralEnergy() const;

    /**
     * @brief Get inversion energy component
     * @return Inversion/out-of-plane energy or 0 if not calculated
     */
    double InversionEnergy() const;
    double STorsEnergy() const;

    /**
     * @brief Get repulsion energy component
     * @return Core-core repulsion energy or 0 if not calculated
     */
    double RepulsionEnergy() const;
    double BondedRepulsionEnergy() const;
    double NonbondedRepulsionEnergy() const;

    /**
     * @brief Get dispersion energy component
     * @return Dispersion correction energy or 0 if not calculated
     */
    double DispersionEnergy() const;

    /**
     * @brief Get Coulomb electrostatic energy component
     * @return Electrostatic energy or 0 if not calculated
     */
    double CoulombEnergy() const;

    /**
     * @brief Get D4 dispersion energy component
     * @return D4 dispersion energy or 0 if not calculated
     * Claude Generated (Jan 2, 2026): D4 dispersion energy accessor
     */
    double D4Energy() const;

    /**
     * @brief Get batm (bonded ATM) energy component
     * @return Batm energy or 0 if not calculated
     * Claude Generated (Jan 17, 2026): Batm energy accessor for 1,4-pairs
     */
    double BatmEnergy() const;

    // Claude Generated (April 2026): PBC accessors for GPU path
    bool hasPBC() const { return m_has_pbc; }
    Eigen::Matrix3d getUnitCellBohr() const {
        constexpr double ANG2BOHR = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
        return m_unit_cell * ANG2BOHR;
    }

    /**
     * @brief Get hydrogen bond energy component
     * @return Hydrogen bond energy or 0 if not calculated
     * Claude Generated (Dec 2025): H-bond energy accessor via ForceField
     */
    double HydrogenBondEnergy() const;

    /**
     * @brief Get halogen bond energy component
     * @return Halogen bond energy or 0 if not calculated
     * Claude Generated (Dec 2025): X-bond energy accessor via ForceField
     */
    double HalogenBondEnergy() const;

    /**
     * @brief Get atm (three-body dispersion) energy component
     * @return ATM energy or 0 if not calculated
     * Claude Generated (Dec 2025): ATM energy accessor via ForceField
     */
    double ATMEnergy() const;

    // =================================================================================
    // Per-Component Gradient Decomposition (Claude Generated February 2026)
    // =================================================================================

    /**
     * @brief Enable per-component gradient storage for validation
     * @param store true to activate component gradient accumulation
     *
     * When enabled, each energy term's gradient contribution is stored separately
     * in addition to the total gradient. Adds ~15% overhead due to matrix snapshots.
     * Only use for validation, not for production MD/optimization.
     */
    void setStoreGradientComponents(bool store);

    /// Per-component gradient getters (Bohr units, same as ForceField internal)
    Matrix GradientBond() const;
    Matrix GradientAngle() const;
    Matrix GradientTorsion() const;
    Matrix GradientRepulsion() const;
    Matrix GradientCoulomb() const;
    Matrix GradientDispersion() const;
    Matrix GradientHB() const;
    Matrix GradientXB() const;
    Matrix GradientBATM() const;
    Matrix GradientATM() const;   ///< ATM three-body dispersion gradient (Claude Generated Mar 2026)

    /// Claude Generated (Mar 2026): Expose dispersion CN chain-rule correction for diagnostics
    Matrix getDispCNCorrection() const;

    // =================================================================================
    // vbond Parameter Access for Verification (Claude Generated November 2025)
    // =================================================================================

    // =================================================================================
    // TWO-PHASE EEQ SYSTEM (Claude Generated November 2025, Session 5)
    // =================================================================================

    /**
     * @brief Phase 1: Calculate topology-aware base charges (qa) via EEQ
     *
     * The two-phase EEQ system separates charge calculation into:
     * 1. TOPOLOGY PHASE: Base charges from atomic properties (electronegativity, hardness)
     *    using coordination-dependent parameters
     * 2. CORRECTION PHASE: Apply dxi, dgam, dalpha corrections for refined accuracy
     *
     * Reference: angewChem 2020, GFN-FF parameter set
     *   - chi[z]: Electronegativity for element z
     *   - gam[z]: Chemical hardness for element z
     *   - alp[z]: Damping parameter for erf(γ*r)/r Coulomb
     *   - cnf[z]: Coordination number correction factor
     *
     * Extended Hückel Theory (EHT) approximation:
     *   qa_i = -χ_i - J_ii + Σ_j(1/(2*J_ij) - 1/(2*r_ij))
     *
     * @param topo_info TopologyInfo structure with coordination numbers, hybridization
     * @return true if Phase 1 charges calculated successfully
     *
     * Output written to: topo_info.topology_charges (qa in Hartree)
     */
    bool calculateTopologyCharges(TopologyInfo& topo_info) const;

    /**
     * @brief Calculate dxi (electronegativity) corrections for Phase 2
     *
     * dxi corrects electronegativity based on:
     * - Local environment (neighbor count, hybridization)
     * - Bonding context (pi-systems, heteroatom effects)
     * - Functional group classification
     *
     * Physical meaning: Electronegativity is NOT constant - it depends on
     * chemical context. Atoms in electron-withdrawing groups become more
     * electronegative.
     *
     * @param topo_info TopologyInfo with topology_charges and hybrid classifications
     * @return true if dxi corrections calculated
     *
     * Output written to: topo_info.dxi (corrections to chi in Hartree)
     */
    bool calculateDxi(TopologyInfo& topo_info) const;

    /**
     * @brief Calculate dalpha (polarizability) corrections for Phase 2
     *
     * dalpha corrects the damping parameter (alpha) based on:
     * - Atomic size changes (coordination-dependent)
     * - Electronic environment (hybridization, charge state)
     * - Pi-system participation
     *
     * Physical meaning: Polarizability (and hence Coulomb damping) adapts to
     * local electronic density. More polarizable atoms in electron-rich
     * environments use different damping.
     *
     * @param topo_info TopologyInfo with coordination numbers and charges
     * @return true if dalpha corrections calculated
     *
     * Output written to: topo_info.dalpha (corrections to alpha)
     */
    bool calculateDalpha(TopologyInfo& topo_info) const;

    /**
     * @brief Calculate charge-dependent alpha (alpeeq) for EEQ
     *
     * Claude Generated (January 2026)
     *
     * Computes charge-dependent alpha values used in EEQ matrix construction.
     * This implements the formula from Fortran gfnff_ini.f90:718-725:
     *
     *   alpeeq(i) = (alpha_base + ff*qa(i))²
     *
     * where ff is element-specific:
     *   - Carbon (Z=6): ff = 0.09
     *   - Nitrogen (Z=7): ff = -0.21
     *   - Group 6 (O,S,Se): ff = -0.03
     *   - Group 7 (Halogens): ff = 0.50
     *   - Main group metals: ff = 0.3
     *   - Transition metals: ff = -0.1
     *
     * Physical meaning: Charge state modifies atomic polarizability (Gaussian width).
     * Positive charges increase alpha (softer, more diffuse) for C/main-group metals,
     * negative charges increase alpha for N/halogens/transition metals.
     *
     * CRITICAL: This must be called AFTER topology_charges are computed,
     * and the resulting alpeeq values are used UNCHANGED in all subsequent
     * EEQ calculations (no iteration).
     *
     * Reference: Fortran gfnff_data_types.f90:128 - "atomic alpha for EEQ, squared"
     *
     * @param topo_info TopologyInfo with topology_charges (Phase 1 charges qa)
     * @return true if alpeeq calculated successfully
     *
     * Output written to: topo_info.alpeeq (squared alpha values)
     */
    bool calculateAlpeeq(TopologyInfo& topo_info) const;

    /**
     * @brief Get generated ForceField parameters
     *
     * Returns the full parameter set generated by the ForceField engine,
     * including bonds, angles, dihedrals, electrostatics, dispersion, etc.
     *
     * This is different from getParameters() which only returns the input JSON.
     *
     * @return JSON object containing all force field parameters
     * @return Empty JSON if ForceField is not initialized
     *
     * Claude Generated - December 27, 2025
     */
    /// JSON view of the native parameter set (bonds/angles/dihedrals/extra_dihedrals/inversions).
    json getForceFieldParameters() const {
        json out;
        out["bonds"] = getBondParameters();
        out["angles"] = getAngleParameters();
        json tors = getTorsionParameters();
        out["dihedrals"] = tors.value("primary", json::array());
        out["extra_dihedrals"] = tors.value("extra", json::array());
        out["inversions"] = getInversionParameters();
        return out;
    }

    /**
     * @brief Phase 2: Calculate final refined charges by solving corrected EEQ
     *
     * Iteratively solves EEQ with corrections:
     *   qa_i_final = qa_i + dxi_i + (dgam correction)
     *   Then re-solve EEQ with modified parameters
     *
     * The correction application is:
     * 1. Modify electronegativity: χ'_i = χ_i + dxi_i
     * 2. Recalculate Coulomb matrix with corrected polarizabilities: α'_i = α_i + dalpha_i
     * 3. Re-solve the linear EEQ system to consistency
     *
     * Convergence: Typically 1-2 iterations for <0.01 e change per atom
     *
     * @param topo_info TopologyInfo with topology_charges, dxi, dalpha corrections
     * @param max_iterations Maximum iterations for EEQ convergence (default: 10)
     * @param convergence_threshold Threshold for charge change (default: 1e-5 Hartree)
     * @return true if Phase 2 refinement successful
     *
     * Output written to: topo_info.eeq_charges (final qa in Hartree)
     */
    bool calculateFinalCharges(TopologyInfo& topo_info, int max_iterations = 10,
                               double convergence_threshold = 1e-5) const;

    /**
     * @brief Set external charges (for testing/validation)
     * @param charges Atomic partial charges to use instead of EEQ-calculated charges
     *
     * Claude Generated (December 2025): Testing utility for charge-dependent validation
     * Bypasses EEQ calculation and uses provided reference charges directly.
     * Call AFTER InitialiseMolecule() but BEFORE Calculation().
     * This method is intended for validation purposes to isolate energy calculation
     * errors from EEQ charge calculation errors.
     */
    void setCharges(const Vector& charges);

    /**
     * @brief Skip Phase-2 EEQ charge recalculation in Calculation()
     * Claude Generated (March 2026): Diagnostic for charge injection tests
     * When true, Calculation() uses whatever charges are currently set (via setCharges())
     * instead of recalculating from EEQ solver. Restore to false after diagnostic.
     */
    void setSkipEEQRecalc(bool skip) { m_skip_eeq_recalc = skip; }
    void setRepDiag(bool diag) { m_rep_diag = diag; }

    /**
     * @brief Get current bond parameters (for validation)
     *
     * Claude Generated (December 2025): Phase 3 - Parameter validation infrastructure
     * Returns the current bond parameters as JSON for systematic validation.
     * Used by Test 5 to compare generated parameters against XTB reference.
     *
     * @return JSON array of bond parameters with fields: i, j, distance, k, fc, etc.
     */
    json getBondParameters() const;

    /**
     * @brief Get current angle parameters (for validation)
     *
     * Claude Generated (December 2025): Phase 3 - Parameter validation infrastructure
     *
     * @return JSON array of angle parameters
     */
    json getAngleParameters() const;

    /**
     * @brief Get current torsion parameters (primary + extra) for validation
     *
     * Claude Generated (March 2026): Per-torsion diagnostic infrastructure
     * Returns JSON with "primary" (from m_dihedrals) and "extra" (from m_extra_dihedrals) arrays.
     *
     * @return JSON object with "primary" and "extra" arrays of torsion parameters
     */
    json getTorsionParameters() const;

    /**
     * @brief Get current inversion parameters for validation
     *
     * Claude Generated (March 2026): Per-torsion diagnostic infrastructure
     *
     * @return JSON array of inversion parameters including potential_type and omega0
     */
    json getInversionParameters() const;

private:
    // Molecular structure (formerly from QMInterface base class)
    int m_atomcount = 0; ///< Number of atoms
    GeoGradMatrix m_geometry; ///< Molecular geometry in Angström — WP-G: RowMajor
    GeoGradMatrix m_gradient; ///< Gradient in Hartree/Bohr — WP-G: RowMajor
    std::vector<int> m_atoms; ///< Atomic numbers (Z values)
    int m_charge = 0; ///< Total molecular charge
    int m_spin = 0; ///< Spin multiplicity

    // GFN-FF specific
    json m_parameters; ///< GFN-FF parameters
    int m_threads = 1; ///< Claude Generated (WP1, May 2026): cached thread count, kept in sync with m_parameters["threads"]
    std::unique_ptr<CxxThreadPool> m_pool; ///< Shared worker pool (topology setup, EEQ, workspace kernels)
    std::unique_ptr<FFWorkspace> m_workspace; ///< Claude Generated (Mar 2026): Unified workspace (replaces ForceField path)

    GeoGradMatrix m_geometry_bohr; ///< Geometry in Bohr (GFN-FF parameters are in Bohr) — WP-G: RowMajor

    // EEQ charge calculation (Dec 2025 - Phase 3: Extraction and delegation)
    std::unique_ptr<EEQSolver> m_eeq_solver; ///< Standalone EEQ solver (replaces embedded EEQ code)
    mutable bool m_eeq_solve_failed = false; ///< F-Q4: EEQ fell back to placeholder charges (fail-loud)

    // Hückel solver for π-bond orders (Jan 2026 - Phase 1: Full Hückel implementation)
    std::unique_ptr<HuckelSolver> m_huckel_solver; ///< Full iterative Hückel solver
    bool m_use_full_huckel = true; ///< Use full Hückel calculation (default: true, set to false for simplified approximation)

    // ATM three-body dispersion terms (extracted from D3/D4 - Claude Generated Jan 2025)

    // Claude Generated (Feb 15, 2026): D4ParameterGenerator kept alive for runtime dc6dcn computation
    // Reference: Fortran gfnff_gdisp0.f90:382-395 - dc6dcn needed for dispersion CN gradient
    mutable std::unique_ptr<D4ParameterGenerator> m_d4_generator;

    // Claude Generated (Mar 2026): ALPB solvation model
    // Reference: Fortran external/gfnff/src/gbsa/gbsa.f90
    // Initialized when solvent != "none", called in Calculation() after EEQ charges
    std::unique_ptr<ALPBSolvation> m_solvation;
    std::string m_solvent = "none";  ///< Solvent name ("none" = gas phase)
    std::string m_solvent_model_label = "ALPB";  ///< "ALPB" | "GBSA" for logging (WP5)

    /**
     * @brief HB/XB dynamic update support for MD simulations
     *
     * Claude Generated (February 2026): Reference geometry tracking for HB/XB list updates
     * Reference: Fortran gfnff_engrad.F90:246-260, gfnff_ini2.f90:715-717
     *
     * During MD simulations, HB/XB pairs can form or break as geometry changes.
     * This struct tracks the reference geometry and triggers list rebuilds when
     * per-atom RMSD exceeds 0.3 Bohr threshold.
     */
    struct HBReferenceGeometry {
        Eigen::MatrixXd reference_positions;  ///< Positions when lists were built (Bohr)
        int nhb_count = 0;    ///< Number of hydrogen bonds
        int nxb_count = 0;    ///< Number of halogen bonds
        bool needs_update = true;  ///< Force update on first call
    };

    mutable std::optional<HBReferenceGeometry> m_hb_reference;  ///< HB/XB reference geometry tracker

    /**
     * @brief Check if HB/XB lists need updating based on geometry change
     *
     * Claude Generated (February 2026)
     * Reference: Fortran gfnff_ini2.f90:715-717
     *
     * Uses per-atom RMSD threshold of 0.3 Bohr to trigger list rebuild.
     * This matches Fortran behavior where lists are rebuilt when geometry
     * changes significantly during MD simulations.
     *
     * @param current_geometry Current geometry in Bohr
     * @return true if lists should be rebuilt
     */
    bool shouldUpdateHBXB(const Eigen::MatrixXd& current_geometry) const;

    // Geometry change detection for intelligent caching
    // WP-G (May 2026): aligned with m_geometry_bohr RowMajor type
    class GeometryChangeDetector {
    private:
        GeoGradMatrix m_last_geometry;
        double m_change_threshold = 1e-6;

    public:
        bool geometryChanged(const GeoGradMatrix& new_geometry) const {
            if (m_last_geometry.rows() != new_geometry.rows() ||
                m_last_geometry.cols() != new_geometry.cols()) {
                return true;
            }
            // Only invalidate cache if change exceeds threshold
            return (m_last_geometry - new_geometry).array().abs().maxCoeff() > m_change_threshold;
        }

        void updateGeometry(const GeoGradMatrix& new_geometry) {
            m_last_geometry = new_geometry;
        }

        void reset() {
            m_last_geometry = GeoGradMatrix();
        }
    };

    bool m_initialized; ///< Initialization status
    bool m_skip_eeq_recalc = false; ///< Skip Phase-2 EEQ recalculation (for charge injection diagnostic)
    bool m_rep_diag = false; ///< Dump repulsion alphanb diagnostic

    // Static-Mode (WP-S1): freeze CN/charges across MD steps to avoid recompute
    bool m_static_charges = false;          ///< If true, skip Phase-2 EEQ after first successful call
    bool m_static_cn = false;               ///< If true, skip CN/dcn/D4-weight recompute after first call
    bool m_static_state_captured = false;   ///< Becomes true once initial CN/charges have been captured

    // WP-S3 (May 2026): runtime state of the EEQ cutoff auto-detection
    bool m_eeq_cutoff_auto_active = false;  ///< true if applyEEQCutoffAutoIfRequested set 30 Bohr

    // WP-P1 (May 2026): cached after each prepareCNAndEEQ() call; consumed by MD diagnostics dump.
    mutable PrepTiming m_last_prep_timing{};
    bool m_force_phase_timing = false;  ///< If true, prepareCNAndEEQ collects per-phase timings even at verbosity < 2

    // WP-D (May 2026): cached raw CN from CNCalculator — reused by calculateCoordinationNumberDerivatives
    // to avoid recomputing the N²-erf loop in dcn step 1.
    mutable Vector m_last_cn_raw{};
    // WP-D Stage C (May 2026): symmetric CN neighbor list — reused by dcn to avoid N² pair scan.
    // Populated by prepareCNAndEEQ when cn_cutoff_bohr > 0; empty otherwise (→ N² fallback).
    mutable std::vector<std::vector<int>> m_last_cn_neighbors{};

    // Claude Generated (April 2026): Periodic Boundary Conditions
    bool m_has_pbc = false;                                              ///< PBC active flag
    Eigen::Matrix3d m_unit_cell = Eigen::Matrix3d::Identity();          ///< Unit cell (Angstrom, from Mol)

    double m_energy_total; ///< Total energy in Hartree
    Vector m_charges; ///< Atomic partial charges
    Vector m_bond_orders; ///< Wiberg bond orders

    // Geometry change detector for intelligent caching
    mutable GeometryChangeDetector m_geometry_tracker;

    // Two-tier topology caching (March 2026)
    // Tier 1: Static topology (bonds, rings, hybridization) - only invalidates on large geometry change
    // Tier 2: Dynamic state (CN, distances) - invalidates on small geometry change
    mutable std::optional<TopologyInfo> m_cached_topology;
    mutable Eigen::MatrixXd m_last_topology_geometry;  // Geometry when topology was last calculated
    mutable Eigen::MatrixXd m_last_eeq_geometry;       // Geometry for which m_charges are currently valid
    mutable bool m_static_topology_valid = false;       // True if static topology is current
    mutable bool m_full_topology_recalculated = false;  // Set by getCachedTopology() on full update
    mutable std::optional<bool> m_external_topology_decision; ///< GPU displacement check result
    mutable std::optional<std::vector<std::pair<int,int>>> m_cached_bond_list;

    std::vector<std::pair<int,int>> m_forced_bonds; ///< External bonds merged with geometric detection

    // React topology mode state (Claude Generated Aug 2026). See docs/GFNFF_REACT_TOPOLOGY.md.
    std::vector<std::pair<int,int>> m_react_bonds; ///< Authoritative bond set (canonical i<j), owns m_forced_bonds in react mode
    Eigen::MatrixXd m_react_ref_geometry; ///< Geometry (Bohr) at the last hysteresis scan
    long m_react_calls = 0; ///< Energy calls since init (drives the scan cadence)
    int m_react_rebuild_count = 0; ///< Bonded-term rebuilds so far
    bool m_react_rebuilt = false; ///< One-shot flag consumed by the GPU/HIP wrappers
    double m_react_form_factor = 1.6; ///< Bond-formation threshold factor (optimistic)
    double m_react_break_factor = 2.6; ///< Bond-keeping threshold factor (conservative)
    int m_react_check_every = 5; ///< Scan every N energy calls (0 = displacement only)
    double m_react_check_disp = 0.25; ///< Scan when any atom moved more than this (Bohr)
    int m_react_refractory_scans = 10; ///< Scans a broken pair must wait before re-forming
    bool m_react_valence_cap = true; ///< Enforce the bond-order-aware valence cap on formation
    int m_react_exchange_scans = 20; ///< Max scans an atom may stay over nominal valence
    double m_react_slack_form_factor = 1.2; ///< Tighter formation radius for slack-consuming bonds
    std::map<int, int> m_react_overvalence_streak; ///< Atom -> consecutive over-valent scans
    std::map<std::pair<int, int>, int> m_react_refractory; ///< Pair -> remaining blocked scans
    bool m_react_owns_bonds = false; ///< react mode: m_forced_bonds is authoritative even when empty
    std::vector<int> m_react_bond_orders; ///< Parallel to m_react_bonds, see reactiveBondOrders()
    std::vector<ReactEvent> m_react_events; ///< Events since the last consumeReactEvents()

    /// Refresh m_react_bond_orders from the cached Hueckel pi orders (after a rebuild / init).
    void refreshReactBondOrders();

    /**
     * @brief React mode: O(N^2) hysteresis scan over all atom pairs; updates m_react_bonds.
     * An existing bond survives while r < break_factor*thr; a new pair becomes a bond
     * at r < form_factor*thr, with thr = (rcov_i+rcov_j)*fat_i*fat_j.
     * @return true if the bond set changed
     */
    bool detectReactiveBondChanges();

    /**
     * @brief React mode: regenerate ALL bonded terms, the repulsion partition, Coulomb
     * parameters and HB/XB lists from the current m_react_bonds; push into both engines.
     * @return true on success
     */
    bool rebuildReactiveTopology();

    // Topology caching mode: "auto" (two-tier caching) or "constant" (never recalculate)
    std::string m_topology_mode = "constant";

    // Claude Generated (March 2026): Topology persistence in param.json
    bool m_cache_topology = true;   ///< Cache Phase-1 EEQ topology in param.json (opt-out)
    /// Claude Generated (Oct 2026): per-chunk buffers of generateRepulsionPairsNative() kept between MD rebuilds (only
    /// when the caller passes reuse storage): no 70 MB free + page-faulting re-allocation per refresh at 6200 atoms.
    mutable std::vector<std::vector<GFNFFRepulsion>> m_rep_part_cache;
    bool m_skip_host_rep_workspace = false;  ///< Oct 2026: GPU wrapper owns the repulsion lists (setSkipHostRepulsionWorkspaceUpdate)
    bool m_skip_host_hb_workspace = false;  ///< Oct 2026: GPU wrapper owns the HB/XB lists (setSkipHostHBWorkspaceUpdate)
    bool m_skip_host_rep_pairs = false;   ///< Oct 2026: the GPU wrapper builds the non-bonded repulsion list
    bool m_skip_host_disp_pairs = false;  ///< WP-A: GPU builds D4 pairs; skip host O(N^2) loop
    bool m_implicit_coulomb_pairs = false; ///< GPU enumerates Coulomb pairs; skip host list
    bool m_print_timing = true;     ///< Print init timing summary at verbosity >= 1

    // Claude Generated (April 2026): Timing for consolidated summary
    double m_param_gen_time_ms = 0.0;       ///< Parameter generation time (generateGFNFFParameterSet) for summary
    mutable double m_topology_time_ms = 0.0; ///< Topology generation time (calculateTopologyInfo) for summary (mutable: set in const method)

    // Claude Generated (May 2026): Profiling report for verbosity-2 param-gen breakdown.
    // mutable so calculateTopologyInfo() (const) can populate sub-phase timings.
    mutable GFNFFParamGenReport m_param_gen_report;

    // Check if geometry change warrants full topology recalculation (vs just dynamic state)
    bool needsFullTopologyUpdate(const Eigen::MatrixXd& geometry_bohr) const;

    // Dynamic state update for Tier 2 caching
    void updateDynamicState(TopologyInfo& topo) const;

    // Claude Generated (March 2026): Heap-stored parameter copy for external consumers.
    // Set in initializeForceField(), consumed once via consumeCachedParameterSet().
    std::unique_ptr<GFNFFParameterSet> m_cached_parameter_set;
    bool m_keep_full_parameter_set = false; ///< GPU wrappers need the full pair lists; CPU keeps bonded terms only
    std::unique_ptr<GFNFFParameterSet> makeParameterSetCache(const GFNFFParameterSet& p) const;
    // EEQ Phase-2 topology input, rebuilt only when the topology changes (B3, Sep 2026)
    std::optional<EEQSolver::TopologyInput> m_eeq_topo_cache;
    unsigned m_eeq_topo_cache_version = 0;
    mutable unsigned m_topology_version = 0; ///< bumped whenever m_cached_topology is (re)assigned (also from const getCachedTopology)

    // Claude Generated (March 2026): Last re-detected HB/XB lists from updateHBXBIfNeeded()
    std::vector<GFNFFHydrogenBond> m_last_hbonds;
    std::vector<GFNFFHalogenBond> m_last_xbonds;
    bool m_hbxb_updated = false;  ///< True if updateHBXBIfNeeded() ran since last check
    bool m_hbxb_fresh = false;    ///< True if HB/XB lists were freshly built during init and geometry is unchanged
    long m_hbxb_update_calls = 0; ///< Task #11: call counter for hb_update_force_every periodic rebuild

    // Claude Generated (Sep 2026): Last rebuilt repulsion pair lists from
    // updateNonbondedRepulsionIfNeeded() — see gfnff_method.cpp and the declaration above.
    std::vector<GFNFFRepulsion> m_last_bonded_reps;
    std::vector<GFNFFRepulsion> m_last_nonbonded_reps;
    bool m_nb_rep_updated = false;    ///< True if updateNonbondedRepulsionIfNeeded() rebuilt the lists since last check
    long m_nb_rep_update_calls = 0;   ///< Call counter for nonbonded_rebuild_every periodic rebuild

    // Claude Generated (Sep 2026): D4 pair-list skin tracking (updateDispersionPairsIfNeeded()).
    Eigen::MatrixXd m_disp_list_ref_geometry; ///< Geometry (Bohr) the current D4 pair list was built at
    Eigen::MatrixXd m_disp_c6_geometry;        ///< Geometry (Bohr) the stored D4 C6 belong to (refreshDispersionC6())
    bool m_disp_pairs_updated = false;         ///< True if the D4 list was rebuilt since last check
    long m_disp_update_calls = 0;              ///< Call counter for the zero-skin step-count fallback
    long m_disp_rebuild_count = 0;             ///< Number of displacement-triggered rebuilds (diagnostic)

    // Claude Generated (Sep 2026): explicit Coulomb list refresh (updateCoulombPairsIfNeeded()).
    bool m_coul_pairs_updated = false;         ///< True if the Coulomb list was rebuilt since last check
    Eigen::MatrixXd m_coul_list_ref_geometry;  ///< Geometry (Bohr) the explicit Coulomb list was built at
    Eigen::MatrixXd m_rep_list_ref_geometry;   ///< Geometry (Bohr) the repulsion list was built at (skin mode)
    long m_rep_rebuild_count = 0;              ///< Number of repulsion-list rebuilds (diagnostic)
    long m_coul_update_calls = 0;              ///< Call counter for nonbonded_rebuild_every

    // Claude Generated (March 2026): State from last prepareCNAndEEQ() call
    Vector m_last_cn;    ///< Coordination numbers
    Vector m_last_cnf;   ///< CN-dependent EEQ factors per atom
    bool m_gpu_path_preallocated = false; ///< True after preAllocateForGPUPath()

    /// Use the Fortran-faithful hybridization (determineHybridizationFortran) instead of
    /// the legacy geometry-first heuristic. Set false to fall back to the pre-Jul-2026
    /// behaviour for bisection. Claude Generated (Jul 2026).
    bool m_use_fortran_hyb = true;

    /// Topology charges (Fortran topo%qa) fed back into the getnb bond-radius shrink
    /// (`rtmp -= qa*fq`, gfnff_ini2.f90:122). Empty on the first pass, which is exactly
    /// Fortran's pass 1 (gfnff_ini.f90:258 sets qa=0 before the q-loop). Filled by the
    /// second q-loop pass. Claude Generated (Jul 2026).
    mutable std::vector<double> m_bond_qa;
    // q-loop pass-2 carry-over of the fragmentation. The reference gates its whole fragment
    // block on `if (topo%nfrag <= 1)` (gfnff_ini.f90:467), so the second pass KEEPS the
    // fragmentation and qfrag found in pass 1 even when the charge-shrunk radii have since
    // merged two fragments into one. Empty nfrag (0) means "detect normally".
    mutable int m_frag_carry_nfrag = 0;
    mutable std::vector<int> m_frag_carry_list;
    mutable std::vector<double> m_frag_carry_qfrag;

    // ===== frag_charge_model ensemble (Claude Generated, Sep 24, 2026) =====
    // See docs/FRAG_CHARGE_MODEL.md. A charge VARIANT is a complete GFNFF instance whose
    // EEQ fragment groups and integer group charges are fixed by m_frag_override; the master
    // blends the variants' energies/gradients: E = sum_c W_c sum_p omega_cp E_cp.
    struct FragOverride {
        bool active = false;
        int nfrag = 0;
        std::vector<int> fraglist;   // 1-based group id per atom
        std::vector<double> qfrag;   // integer charge per group
    };
    FragOverride m_frag_override;
    bool m_frag_is_variant = false;      ///< this instance is a variant (never nests)
    mutable int m_frag_split_pass = 0;   ///< q-loop pass that split the fragments (1 or 2), set by calculateTopologyInfo
    bool m_frag_ensemble = false;        ///< frag_charge_model == ensemble
    double m_frag_s_max = 1.1;
    double m_frag_tau_eh = 1.0 / 627.5094740631;
    double m_frag_sigma = 0.05;          ///< e, width of the free-charge softmax over carrier classes
    int m_frag_max_edges = 4;
    int m_frag_max_placements = 6;
    int m_frag_variant_memory_mb = 8192; ///< budget for the held variants (Claude Generated, Oct 2026)
    bool m_frag_enumerate = false;       ///< validation: full enumeration of the placements (pre-Oct-2026 path)
    bool m_frag_warned_memory = false;
    bool m_frag_warned_select = false;
    struct FragContact { int i, j; double thr, L, dL; };   // dL = dL/ds
    struct FragEdge { int f, g; double lambda = 1.0; std::vector<FragContact> contacts; };
    struct FragVariant {
        std::string key;
        std::unique_ptr<GFNFF> ff;
        double energy = 0.0;
        Matrix gradient;
        Vector charges;
        FFEnergyComponents comp;
        bool ok = false;
    };
    std::vector<FragVariant> m_frag_variants;          // cache, keyed by FragVariant::key
    unsigned m_frag_variants_topo_version = ~0u;
    Vector m_frag_free_q;                               // free single-constraint Phase-1 charges (placement pre-selection)
    bool m_frag_blend_valid = false;
    FFEnergyComponents m_frag_blend_comp;
    int m_frag_last_nvariants = 0;
    bool m_frag_warned_edges = false;
    double fragPass1Threshold(int i, int j) const;     ///< getnb bond threshold of the pass that split the fragments, Bohr
    std::vector<FragEdge> fragWindowEdges(const std::vector<int>& fraglist, int nfrag) const;
    GFNFF* fragVariant(const std::string& key, const std::vector<int>& group_of_atom, int ngroups,
                       const std::vector<double>& qgroup);
    double fragEnsembleBlend(bool gradient, double e_master);
    double calculationSingle(bool gradient);            ///< the reference (single-variant) Calculation()
    CNDerivStore m_last_dcn; ///< CN derivatives (gradient only). Claude Generated (WP4, May 2026): pair-list replaces std::vector<SpMatrix>

    // WP-FF-DistMatrix-Sharing (May 2026): shared packed-triangular distance array.
    // Filled by computeSharedDistances() once per energy call; its only reader is the EEQSolver
    // phase 2 (setExternalDistances). Layout: lower-triangular, index via
    // GFNFF::triIdx(i, j) = i*(i+1)/2 + j (i > j); the diagonal slots are allocated, never written
    // or read. Oct 2026: the squared-distance twin m_shared_sqrab had no reader and is gone.
    mutable Eigen::VectorXd m_shared_srab;
    mutable int             m_shared_dist_N = 0;
    mutable bool            m_shared_dist_too_close = false; ///< any r_ij < 1e-10 Bohr in m_shared_srab

    // Conversion factors
    static constexpr double HARTREE_TO_KCAL = 627.5094740631;
    // Sep 2026: were 0.5291772105638411 (~CODATA 2014) while m_geometry_bohr is built with
    // CurcumaUnit (CODATA 2018), so the numerical gradient's Bohr->Angstrom round trip was
    // off by 6.4e-10 relative. One constant for geometry now - Claude Generated.
    static constexpr double BOHR_TO_ANGSTROM = CurcumaUnit::Length::bohr_radius_or_legacy(0.5291772105638411);
    static constexpr double KCAL_TO_HARTREE = 1.0 / 627.5094740631;
    static constexpr double ANGSTROM_TO_BOHR = CurcumaUnit::Length::angstrom_to_bohr_or_legacy(1.0 / 0.5291772105638411);

    /**
     * @brief Element-specific radius scaling factors (fat array from gfnff_ini2.f90:76-97)
     *
     * Claude Generated (January 2026) - Phase 3: Element-specific neighbor detection
     * Applied in bond detection: threshold = 1.3 * (rcov_i + rcov_j) * fat[Z_i] * fat[Z_j]
     *
     * Default: 1.0 for all elements
     * Special adjustments for specific elements to improve bond detection accuracy
     */
    static constexpr double fat[87] = {
        0.0,   // 0: placeholder (atom numbers start at 1)
        1.02,  // 1: H
        1.00,  // 2: He
        1.00,  // 3: Li
        1.03,  // 4: Be
        1.02,  // 5: B
        1.00,  // 6: C
        1.00,  // 7: N
        1.02,  // 8: O
        1.05,  // 9: F
        1.10,  // 10: Ne
        1.01,  // 11: Na
        1.02,  // 12: Mg
        1.00,  // 13: Al
        1.00,  // 14: Si
        0.97,  // 15: P
        1.00,  // 16: S
        1.00,  // 17: Cl
        1.10,  // 18: Ar
        1.02,  // 19: K
        1.02,  // 20: Ca
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  // 21-30: Sc-Zn
        1.00, 1.00, 1.00, 0.99,  // 31-34: Ga-Se (34: Se = 0.99)
        1.00, 1.00, 1.00,        // 35-37: Br-Rb
        1.02,  // 38: Sr
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  // 39-49: Y-In
        1.01,  // 50: Sn
        0.99,  // 51: Sb
        0.95,  // 52: Te
        0.98,  // 53: I
        1.00, 1.00,  // 54-55: Xe-Cs
        1.02,  // 56: Ba
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  // 57-66: La-Dy
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,        // 67-75: Ho-Re
        1.02,  // 76: Os
        1.00, 1.00, 1.00, 1.00, 1.00,  // 77-81: Ir-Tl
        1.06,  // 82: Pb
        0.95,  // 83: Bi
        1.00, 1.00, 1.00   // 84-86: Po-Rn
    };
};

/**
 * @brief GFN-FF Architecture: Two-Phase Implementation
 *
 * PARAMETER GENERATION (this class):
 * - generateTopologyAwareBonds()     → JSON["bonds"]
 * - generateTopologyAwareAngles()    → JSON["angles"]
 * - generateGFNFFTorsions()          → JSON["dihedrals"]
 * - generateGFNFFInversions()        → JSON["inversions"]
 * - generateGFNFFDispersionPairs()   → JSON["gfnff_dispersions"]
 * - generateGFNFFRepulsionPairs()    → JSON["gfnff_repulsions"]
 * - generateGFNFFCoulombPairs()      → JSON["gfnff_coulombs"]
 *
 * TERM CALCULATION (ForceFieldThread):
 * - CalculateGFNFFBondContribution()
 * - CalculateGFNFFAngleContribution()
 * - CalculateGFNFFDihedralContribution()
 * - CalculateGFNFFInversionContribution()
 * - CalculateGFNFFDispersionContribution()
 * - CalculateGFNFFRepulsionContribution()
 * - CalculateGFNFFCoulombContribution()
 *
 * To add new terms: Modify BOTH this class (generation) AND ForceFieldThread (calculation)
 * See: src/core/energy_calculators/ff_methods/CLAUDE.md for detailed checklist
 */
