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

#include "json.hpp"
#include <map>  // rev-gfnff per-element overrides (Sep 2026)
#include "src/core/config_manager.h"
#include "src/core/parameter_macros.h"
#include "src/core/energy_calculators/ff_methods/forcefield.h"

#include "src/core/energy_calculators/ff_methods/ff_workspace.h"  // Claude Generated (Mar 2026): Unified workspace
#include "src/core/energy_calculators/ff_methods/eeq_solver.h"  // EEQ charge calculation (Dec 2025 - Phase 3)
#include "src/core/energy_calculators/ff_methods/huckel_solver.h"  // Full Hückel calculation (Jan 2026 - Phase 1)
#include "src/core/energy_calculators/ff_methods/gfnff_param_tables.h"  // Claude Generated (Sep 2026): runtime parameter tables (rev-gfnff WP1a)
#include "src/core/energy_calculators/dispersion/d4param_generator.h"  // Claude Generated (Feb 15, 2026): D4 for dc6dcn gradient
#include "src/core/energy_calculators/ff_methods/alpb_solvation.h"  // Claude Generated (Mar 2026): ALPB solvation
#include "src/core/global.h"
#include "src/core/functional_groups.h"
#include "src/core/periodic_table.h"
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
    double over_coord = 0;   // rev-gfnff stage 1 (Sep 2026)
    double sqe_hardness = 0; // rev-gfnff stage 2 (Sep 2026)
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
PARAM(accuracy, String, "normal", "Accuracy profile: loose|normal|medium|high. Maps to EEQ and CN parameters.", "Basic", {},
        "enum=loose|normal|medium|high")
PARAM(allow_unconverged_charges, Bool, false, "Allow calculation to continue with unconverged EEQ charges (warn instead of abort).", "Advanced", {})
PARAM(skip_phase2, Bool, false, "Skip Phase 2 EEQ refinement and use Phase 1 topology charges directly. Faster but less accurate.", "Advanced", {})
// ORIGIN NOTE (do not duplicate on merge): nh_linear_fix and its guard in
// determineHybridizationFortran() were first added on the `confsearch` branch (51830efa,
// Aug 29 2026), where the artefact was found, and ported here in Sep 2026. The guard has
// since been REFINED here (it now also requires the nitrogen's heavy partner to be
// branched), so the two branches are no longer identical: THIS version is the newer one,
// take it on merge. The PARAM text below is still byte-identical to confsearch's.
PARAM(nh_linear_fix, Bool, true, "Do not let the angle-only GEODEP rule (input angle > 160 deg -> sp) promote a 2-coordinate nitrogen that carries a hydrogen to sp hybridisation. Why: the reference rule (gfnff_ini.f90, gen%linthr) declares ANY near-linear input angle linear-by-design, so a thermally stretched =N-H (measured: a guanidine imine N-H at 179 deg in a hot MD snapshot) is re-perceived as sp, gets theta0=180, and the distortion becomes its own equilibrium -- the structure optimises INTO the artefact and appears ~160 kJ/mol too deep (xtb 6.7.1 reproduces this with -276 kJ/mol, so it is an inherited method defect, not a port bug). In a conformer search, where every snapshot optimisation derives its own topology, one such event founds a self-reinforcing family (measured: 75 percent of a WEKLQ pool within three temperature stages). The genuine sp cases of an N-H nitrogen (H-N=C isocyanide-like, R-N=N terminal, metal nitriles, azides) are all caught by the STRUCTURAL rules that run before the angle fallback and are unaffected by this guard. Set false for bit-faithful reference (xtb/pprcht) behaviour, e.g. for validation against the Fortran implementations.", "Advanced", {})
PARAM(frag_charge_autodetect, Bool, false, "For a CHARGED molecule that falls into exactly TWO fragments, try both placements of the net charge and keep the one with the lower EEQ electrostatic energy. Off by default because that is NOT what the reference does: its auto-detection block (gfnff_ini.f90, nfrag==2 branch) is gated on sum(qfrag(2:nfrag)) > 999 while qfrag is pre-initialised to [charge, 0, ...], so the block is dead code in both pprcht and xtb and the effective rule is 'whole charge on fragment 0'. Enabling the trial changes which fragment carries the charge and can be very wrong: GMTKN55 AHB21/21 (formate ... HF) then puts the -1 on the two-atom HF fragment, giving its hydrogen a charge of -0.52 and shifting the Coulomb term by 237 kcal/mol. Enable only to reproduce curcuma's pre-Sep-2026 behaviour or to experiment with the placement rule.", "Advanced", {})
PARAM(frag_charge_model, String, "ensemble", "How the net charge of a CHARGED molecule that GFN-FF perceives as several fragments is placed. reference = the pprcht/xtb rule: the whole charge on fragment 0, i.e. on the fragment that contains atom 1, one integer placement, identical to every release so far. ensemble = chemistry-aware and continuous: every integer placement of the charge over the fragments is a complete GFN-FF evaluation with parameters consistent with that placement, placements are chosen by electron-count parity - fewest odd-electron fragments - then weighted by their free Phase-1 EEQ charge, and chemically identical placements by their energy with temperature frag_charge_tau, and fragments that are close to their pass-1 bond threshold are blended continuously with a merged corner that treats them as one EEQ group - weight 1-lambda, lambda = smootherstep of r/r_thr between 1 and frag_charge_s_max - so the energy is continuous where the fragment count changes. Engaged only for charged systems with two or more fragments; CPU; not in react topology mode. See docs/FRAG_CHARGE_MODEL.md, test_cases/revgfnff/_log/FRAG_CHARGE_STATUS.md.", "Advanced", {})
PARAM(frag_charge_s_max, Double, 1.0, "frag_charge_model ensemble: end of the contact window in units of the pass-1 bond threshold r_thr of each atom pair. Between s = r/r_thr = 1, where pass 1 splits the pair, and s_max two fragments are partly merged; beyond s_max they are fully separated and only integer placements remain. 1.0 switches the window off: only the placement rule acts, and the energy keeps the reference rule's step at the split.", "Advanced", {})
PARAM(frag_charge_tau, Double, 1.0, "frag_charge_model ensemble: temperature in kcal/mol of the Boltzmann average over integer charge placements. Small = the lowest-energy placement wins; placements closer than a few tau are averaged, which keeps the energy and gradient smooth at degeneracies such as a symmetric X...X- pair.", "Advanced", {})
PARAM(frag_charge_sigma, Double, 0.05, "frag_charge_model ensemble: width in elementary charges of the weighting between chemically DIFFERENT charge carriers that the electron-count parity rule leaves tied, by the free single-constraint Phase-1 EEQ charge each carrier holds. Chemically identical carriers, e.g. the two ends of a symmetric X...X- pair, are weighted by their energies with frag_charge_tau instead.", "Advanced", {})
PARAM(frag_charge_atomic_ea, Bool, false, "frag_charge_model ensemble, opt-in (Sep 2026, test_cases/revgfnff/_log/I2_CLF_STATUS.md): for a net charge of -1 whose parity-tied carrier candidates are all SINGLE atoms, weight the chemically different carriers by the experimental atomic electron affinity instead of the free Phase-1 EEQ charge. The asymptote gap E(A-) + E(B) - E(A) - E(B-) is exactly EA(B) - EA(A), while the EEQ free charge measures electronegativity, and the two order F and Cl oppositely: ClF- then dissociates to Cl- + F (the lower state) instead of Cl + F-. Falls back to the free-charge rule when any candidate is polyatomic or its element has no tabulated value.", "Advanced", {})
PARAM(frag_charge_ea_sigma, Double, 0.02, "frag_charge_atomic_ea: width in eV of the softmax over the atomic electron affinities of the candidate carriers (0.02 eV: a hard choice except at near-ties such as Br / F, 0.04 eV apart).", "Advanced", {})
PARAM(frag_charge_max_edges, Int, 4, "frag_charge_model ensemble: maximum number of fragment-fragment contacts inside the window that are blended - 2^n corners are evaluated. Contacts beyond this count, the most separated first, are treated as fully separated and a warning is printed.", "Advanced", {})
PARAM(frag_charge_max_placements, Int, 6, "frag_charge_model ensemble: maximum number of integer charge placements evaluated per corner. With more groups the candidates are pre-selected by the free single-constraint Phase-1 EEQ charge of each group.", "Advanced", {})
PARAM(cn_cutoff_bohr, Double, 10.0, "CN neighbor list cutoff radius in Bohr (reference cnthr=100 Bohr^2=10 Bohr). 0 = use accuracy-based threshold instead.", "Advanced", {})
PARAM(cn_accuracy, Double, 1.0, "CN accuracy for threshold calculation (cnthr = 100 - log10(acc)*50). Only used when cn_cutoff_bohr = 0. Set to 0 for full O(N^2) reference mode.", "Advanced", {})
PARAM(solve, String, "auto",
      "EEQ solver method: lu, schur_cholesky, pcg, auto. Passed to eeq_solver.", "Algorithm", {})
PARAM(eeq_max_iterations, Int, -1,
      "Max EEQ PCG iterations. -1 = use adaptive default (200 small, 5000 large). Passed to eeq_solver.", "Algorithm", {})
PARAM(eeq_tolerance, Double, -1.0,
      "EEQ PCG tolerance. -1 = use adaptive default (1e-10 small, 1e-6*||rhs|| large). Passed to eeq_solver.", "Algorithm", {})
PARAM(eeq_accuracy, Double, 1e-6,
      "Target EEQ charge accuracy (e). Sets pcg_large_system_tol_factor. Higher = faster but less accurate.", "Algorithm", {})
// WP-C (May 2026): canonical PARAM definition lives in eeq_solver.h:965 with default 0.0
// (matches Fortran goed_gfnff). Forwarded to eeq_solver sub-config via
// gfnff_method.cpp::forwardEEQSolverParams. Removing the duplicate here eliminates the
// 30.0-vs-0.0 default discrepancy that overrode the Fortran-matching behaviour.
// PARAM(eeq_distance_cutoff, Double, 30.0, ...) — REMOVED, see eeq_solver.h:965
PARAM(gpu_block_size, Int, 0,
      "GPU kernel block size (0 = adaptive, 32/64/128/256/512). 512 = max occupancy. Passed to ff_workspace_gpu.", "Advanced", {})
PARAM(nb_cell_list_min_atoms, Int, 800,
      "Min atom count to use SpatialCellList for non-bonded neighbor detection (HB/XB and Coulomb when eeq_distance_cutoff>0). Below this falls back to O(N²) loop. 0 = always use cell list.", "Advanced", {"hb_cell_list_min_atoms"})
PARAM(hb_parallel_min_pairs, Int, 500,
      "Min AB-pair count to parallelise HB detection via CxxThreadPool. 0 = always parallel (if pool available). -1 = never parallel.", "Advanced", {})
PARAM(static_charges, Bool, false,
      "Skip Phase-2 EEQ refinement after initialisation; reuse initial charges. Saves ~30 ms/step. Invalid for charge-transfer or ionic dynamics.", "Performance", {})
PARAM(static_cn, Bool, false,
      "Cache CN, dcn, CNF and D4 Gaussian-weights/dc6dcn after first call. Saves ~25 ms/step CPU, ~5-10 ms/step GPU. Invalid for sp2/sp3 changes.", "Performance", {})
PARAM(static_all, Bool, false,
      "Shorthand: enables static_charges=true AND static_cn=true. Use only for stable production NVT/NPT in equilibrium regime.", "Performance", {})
PARAM(eeq_distance_cutoff_auto, Bool, false,
      "Auto-enable eeq_distance_cutoff=30 Bohr after Phase-1 when nfrag==1 and max|q|<0.5 e. Saves ~12 ms/step polymer. Falls back to 0.0 for ionic/multi-fragment systems.", "Performance", {})
PARAM(dispersion_cutoff_bohr, Double, 0.0, "Cutoff (Bohr) for D4 dispersion pair-list. 0 = full O(N^2) (Fortran-parity). Recommended for large systems: 15.0. Energy drift < 1 muEh at 15 Bohr. When active, CN-derivative stencil is extended to cover the cutoff range.", "Performance", {})
PARAM(dispersion_c6_update, Bool, true, "Recompute every stored D4 pair C6 from the current coordination numbers on every evaluation (stale-CN fix B, Known Issue #32). Before this fix C6 stayed at the setup geometry while the dispersion gradient already used dC6/dCN at the current CN, so a reused calculator (optimisation, MD, finite differences, batch reuse) evaluated a frozen-C6 energy that did not match a single point at the same geometry. Single points are unaffected. false restores the old frozen-C6 behaviour for comparisons.", "Algorithm", {})
PARAM(disp_half_contraction, Bool, true,
      "Lever 3 Opt B: per-atom half-contraction fast path for the D4 dispersion C6 and dc6dcn build. About 7x faster inner contraction on large systems; reassociates the FP sum at ~1e-16 so energy matches to ~1e-10 Eh and gradient to ~1e-7. Set false for strictly bit-identical reproductions.", "Performance", {})
PARAM(eeq_refactor_eps_bohr, Double, 0.05,
      "WP-EEQ-Cache: EEQ Cholesky refactorization threshold (max atom displacement, Bohr). "
      "Skips O(N^3) factorization when geometry change below this. "
      "Set to 0.0 (or negative) to disable the cache entirely — bit-identical to pre-WP. "
      "Forwarded to eeq_solver.eeq_refactor_eps_bohr.", "Performance", {})
PARAM(eeq_refactor_force_every, Int, 0,
      "WP-EEQ-Cache: Force EEQ Cholesky refactorization every N steps. "
      "0 = geometry-triggered only. Recommended: 100 for long MD. "
      "Forwarded to eeq_solver.eeq_refactor_force_every.", "Performance", {})
PARAM(eeq_refine_iters, Int, 1, "A4: iterative-refinement steps when the EEQ solve reuses a cached Cholesky factor. Keeps charges exact for the current geometry at O(N^2) cost, so a loose refactor threshold does not corrupt the gradient. 0 disables. Forwarded to eeq_solver.eeq_refine_iters.", "Performance", {})
// NOTE: each PARAM is kept on a SINGLE line on purpose. The param_parser clears its
// buffer on the first ')' it sees, so a multi-line PARAM whose help text contains '(...)'
// is silently dropped from the registry (see eeq_refactor_* above). Single-line is safe.
PARAM(gpu_cn_pair_regen, Bool, true, "Task 10: regenerate the GPU CN-derivative pair list when the topology-displacement check fires (atoms moved past 0.5 Bohr). true keeps the list fresh during opt/MD; false reverts to the legacy build-once-and-latch behaviour.", "Performance", {})
PARAM(gpu_cn_pair_regen_every, Int, 0, "Task 10: force GPU CN-derivative pair-list regeneration every N gradient steps (0 = topology-triggered only). MD safety net for slow drift below the displacement threshold.", "Performance", {})
PARAM(gpu_cn_pair_cutoff_factor, Double, 2.5, "Task 10: GPU CN-derivative pair-list cutoff = factor*(rcov_i+rcov_j). Default 2.5 (erf-CN derivative ~exp(-126) beyond this). Larger = more pairs toward the CPU 40 Bohr reach and slower; smaller = faster with more truncation error.", "Performance", {})
PARAM(hb_accuracy, Double, 0.1, "Task 11: HB-list accuracy driving the Fortran thresholds hbthr1 = 200 - log10(acc)*50, hbthr2 = 400 - log10(acc)*50 (Bohr^2). Smaller = larger cutoffs = more HB pairs = more accurate but slower. 0.1 reproduces the gfnff reference.", "Performance", {})
PARAM(hb_thr1_bohr2, Double, 0.0, "Task 11: direct override of hbthr1, the A-B distance-squared cutoff for nhb2 detection (Bohr^2). 0 = derive from hb_accuracy.", "Performance", {})
PARAM(hb_thr2_bohr2, Double, 0.0, "Task 11: direct override of hbthr2, the A-H-B sum-of-squares cutoff for nhb1 detection (Bohr^2). 0 = derive from hb_accuracy.", "Performance", {})
PARAM(hb_update_rmsd_bohr, Double, 0.3, "Task 11: per-atom RMSD (Bohr) that triggers an HB/XB list rebuild. 0.3 reproduces the gfnff reference (gfnff_ini2.f90:717). Smaller = rebuild more often = less near-threshold staleness in MD but slower.", "Performance", {})
PARAM(hb_update_force_every, Int, 10, "Force an HB/XB list rebuild every N energy evaluations, in addition to the RMSD trigger (0 = RMSD-triggered only). Default 10 since Sep 2026: the RMSD trigger computes sqrt(sum d^2)/N (faithful to gfnff_ini2.f90:717, kept unchanged pending a reference check - see TODO.md), which is sqrt(N) smaller than a per-atom RMSD, so it practically never fires beyond a few dozen atoms (triose, 66 atoms, 200 fs at 800 K: per-atom RMSD 2.2 Bohr, trigger value 0.27, no rebuild; the stale list held 1602 H-bond triples against 1488 in a fresh build, 6e-6 Eh apart). One rebuild costs ~100 ms on polymer_2x (7320 atoms, 150k triples, ~20-50 triples change per 0.5 fs step), so every 10 steps is ~0.7 percent of an MD step. 1 = every step (exact classification, ~7 percent there).", "Performance", {})
PARAM(hb_min_pair_energy_eh, Double, 1e-9, "Sep 2026 (docs/GFNFF_PERFORMANCE_LEVERS.md lever #1): skip allocating a GFNFFHydrogenBond for a case-1 (unbound A...H...B) candidate when the EXACT |E_HB| the energy kernel would compute is below this (Eh), evaluated at detection time from the same formula calcHydrogenBonds uses - not an approximation, and not the rejected raw-distance cut from the same doc (which moved energy 0.79 Eh). Case 2/3/4 (donor-bonded H) are NEVER pruned by this, regardless of value: their acceptor feeds bond_hb_data/hb_cn_H, a geometric quantity uncorrelated with |E_HB| - pruning them shifted the bond term by 0.18 kcal/mol on a test system even though the HB term itself barely moved. 0 disables. On a 3000-water/14640-atom system the case-1 candidate list alone was over 1M entries dominated by long-range, damping-suppressed near-zero contributors.", "Performance", {})
PARAM(xb_min_pair_energy_eh, Double, 1e-9, "Sep 2026 (docs/GFNFF_PERFORMANCE_LEVERS.md lever #1, XB counterpart): skip allocating a GFNFFHalogenBond when the EXACT |E_XB| the energy kernel would compute is below this (Eh) - the XB formula has no case-dependent branch, so this is exact, not an approximation. 0 disables.", "Performance", {})
PARAM(nonbonded_rebuild_every, Int, 1, "The non-bonded repulsion pair list is built from a hard 20 Bohr distance cutoff; a pair that starts beyond it and diffuses closer during MD was never re-evaluated (fixed Sep 2026 - see GFNFF::updateNonbondedRepulsionIfNeeded). This rebuilds it from the current geometry every N energy evaluations. 1 (default) = every step, unconditionally correct - no pair can cross the cutoff undetected. The same schedule drives the explicit Coulomb list of eeq_distance_cutoff > 0 and the D4 list when dispersion_cutoff_bohr <= 50 leaves it without a skin. Raise only after confirming the rebuild cost matters for your system size; a stale list can let two atoms pass through the repulsive wall with zero force, which is far more expensive to debug than the rebuild.", "Performance", {})
PARAM(nonbonded_skin_bohr, Double, 0.0, "Verlet skin (Bohr) for the non-bonded repulsion list and the explicit Coulomb list of eeq_distance_cutoff > 0. 0 (default) = rebuild on the nonbonded_rebuild_every step count. > 0 = build the lists that much wider than their kernel cutoff (20 Bohr repulsion, eeq_distance_cutoff Coulomb) and rebuild only once some atom moved more than skin/2 since the last build - exact, since no pair can then cross the kernel cutoff unseen (Verlet 1967). Energies are unchanged up to floating-point summation order (the longer list shifts the thread partition); only the rebuild schedule and the list length change. See docs/GFNFF_PAIR_LIST_REFRESH.md for measured costs.", "Performance", {})
PARAM(dispersion_c6_update, Bool, true, "Recompute every stored D4 pair C6 from the current coordination numbers on every evaluation (CPU: the 'Stale-CN fix B' block at the end of GFNFF prepare, which also covers stored rev-gfnff corner lists; GPU: the device-side refresh). Before Sep 2026 C6 stayed at the setup geometry while the dispersion gradient already used dC6/dCN at the current CN, so MD and optimisation energies drifted away from a single point at the same geometry (triose: 0.18 kcal/mol after 200 fs at 800 K). Single points are unaffected. false restores the old frozen-C6 behaviour for comparisons.", "Algorithm", {})
PARAM(eeq_mixed_precision, Bool, false, "WP-B GPU only: factor the EEQ Coulomb matrix in FP32 then refine the solution with the FP64 residual, dsposv-style, for full FP64 accuracy at a fraction of the FP64-factor cost on FP64-weak GPUs. Opt-in on CUDA and ROCm (default OFF; enable per card after measuring). Applies to the factor-dominated few-fragment solve paths; the many-fragment general path stays FP64.", "Performance", {})
PARAM(eeq_mixed_precision_iters, Int, 2, "WP-B GPU only: number of FP64-residual / FP32-correction refinement steps for eeq_mixed_precision. Minimum 1. Two steps reach FP64 accuracy on the validation set.", "Performance", {})
PARAM(coulomb_implicit, Bool, true, "CPU: evaluate the N^2/2 Coulomb pairs on the fly from the per-atom EEQ charges and alpeeq instead of building and storing a pair list. The stored list costs 128 bytes per pair - 3.4 GB and ~0.5 s of pure write bandwidth at 7320 atoms, which threading does not remove (measured). Energies agree with the stored path to rounding; set false for the stored list (e.g. to compare). Not used with eeq_distance_cutoff > 0, where the list is already short. DEFAULT TRUE since Sep 18, 2026. Two consequences, both measured: -gfnff.dump_params no longer contains a Coulomb list, so its md5 changes by construction (the energies do not), and the partition is by ATOM instead of by pair, so the reduction order can change at -threads > 1. Set false for the stored list.", "Performance", {})
PARAM(coulomb_r_cut, Double, 100.0, "GFN-FF electrostatics: per-pair distance cutoff in Bohr. The reference (Fortran goed_gfnff) has NO cutoff; 100 Bohr was chosen as an 'effective no-cutoff' because no pair of the validation sets reaches it. That assumption breaks for a system wider than ~53 Angstrom: the cutoff is HARD (no switching), so a pair crossing it changes the energy discontinuously - measured on a 500-water cluster, moving one oxygen by 0.0005 Angstrom jumped the energy by 5 kJ/mol and made the analytic gradient wrong by 0.96 Eh/Angstrom at that atom. In MD such crossings inject energy. Raise it (or set a very large value) for systems above ~50 Angstrom; the price on polymer_2x (7320 atoms) is a single point 3.6 -> 5.4 s. Below ~53 Angstrom nothing changes, so every reference set is unaffected. See docs/MD_LARGE_SYSTEMS.md.", "Performance", {})
PARAM(gpu_coulomb_implicit, Bool, true, "GPU only: the device enumerates all Coulomb atom pairs itself (per-atom gather, gamma_ij from per-atom alpeeq) instead of reading an N^2/2 pair list built on the host. Saves the host list (24.5 M pairs / 2.7 GB and ~1.5 s at 7320 atoms). Not used with eeq_distance_cutoff > 0. Set false for the stored pair list.", "Performance", {})
PARAM(gpu_disp_pairs_on_device, Bool, false, "WP-A GPU only: build the D4 dispersion pair list on the device via a two-pass enumeration plus per-pair C6 contraction, replacing the host O(N^2) GenerateDispersionPairsNative loop and the per-build H2D upload. Default OFF keeps the proven host build. Bit-identical to the host list up to the FP order of the device Gaussian weights. Measured Sep 2026 on polymer_2x (7320 atoms): pair build 630 -> 32 ms, SP wall 6.8 -> 5.6 s, energy identical; see docs/GPU_TUNING.md.", "Performance", {})
PARAM(eeq_rocm_cpu_fragment_threshold, Int, 16, "ROCm GFN-FF only: fragment count at or above which the device EEQ solve is replaced by the exact CPU PCG block-Jacobi warm-start solver, whose O(N^2 k) cost beats the device dense N x N Cholesky O(N^3) for solvent boxes and keeps ROCm charges identical to the CPU path. Set 0 to always use the device solve.", "Performance", {})
// Implicit solvation (WP5, Claude Generated June 2026). Registering these here is
// what makes -gfnff.solvent reach GFNFF::InitialiseMolecule (the value was silently
// ignored before, since the gfnff module declared no solvent PARAM). Use the dotted
// -gfnff.solvent form; the flat -solvent is ambiguous across providers.
PARAM(solvent, String, "none",
      "Implicit solvent for the native GFN-FF ALPB/GBSA model (e.g. 'water', 'dmso', "
      "'acetone', 'chloroform'). 'none' (default) runs gas phase. Born electrostatics "
      "at the EEQ charges + CDS surface term + state shift; the reaction field is a "
      "post-hoc add-on (the EEQ charges do not yet feel the solvent). Use the dotted "
      "-gfnff.solvent (the flat -solvent is ambiguous across providers).", "Solvation", {})
PARAM(solvent_model, String, "alpb",
      "GFN-FF implicit solvation model: 'alpb' (default, P16 Born kernel) or 'gbsa' "
      "(Still kernel, no shape term). CPCM is not implemented natively. Legacy numeric "
      "codes (2=gbsa, 3=alpb) are also accepted.", "Solvation", {})
// React topology mode (Claude Generated Aug 2026): event-driven reactive bond topology.
// Bonds may form and break during MD; all bonded terms are rebuilt at change events.
// See docs/GFNFF_REACT_TOPOLOGY.md. PARAMs stay single-line, see note above.
PARAM(topology_mode, String, "auto", "Topology mode: auto = adaptive two-tier caching, constant = frozen after init, react = bond topology is re-detected with hysteresis during MD and all bonded terms are rebuilt at change events. default is accepted as an alias for auto.", "Basic", {})
PARAM(reuse_topology_check, Bool, false, "Re-validate a carried-over force-field topology against the geometry it is used on: the bond graph the interaction lists (bonds, angles, torsions, repulsion partition, EEQ fragments) were built from is compared with the one the current geometry yields, and they are rebuilt for this frame - with a warning - when the two differ. Enabled automatically by -batch_reuse_topology true, which is what makes calculator reuse safe for bond-stretch scans, dissociation curves, conformer series and multi-molecule batches; a homogeneous series (frames of one MD trajectory of one molecule) never trips it, so its numbers and cost are unchanged. false (the default, and the pre-Sep-2026 behaviour) trusts the first frame's topology unconditionally, which is silently wrong for structurally different frames by 18-117 kcal/mol at one geometry (test_cases/revgfnff/_log/OUTLIER_STATUS.md section F); pass -gfnff.reuse_topology_check false to opt back into that. MD and geometry optimisation freeze the topology on purpose (one molecule, one trajectory) and never enable it.", "Basic", {})
PARAM(react_bond_form_factor, Double, 1.6, "React mode: a non-bonded pair becomes a bond when r < factor * covalent-radius sum * element fat scaling. Optimistic on purpose: the Gaussian bond well is weak at this distance and formation is expected mid-collision. Must stay below react_bond_break_factor and below typical hydrogen-bond contact distances.", "Reactive", {})
PARAM(react_bond_break_factor, Double, 2.6, "React mode: an existing bond is removed when r > factor * covalent-radius sum * element fat scaling. Conservative on purpose: the bond is kept until its Gaussian well has largely decayed, so removal causes only a small energy jump. The wide gap to react_bond_form_factor is the hysteresis that prevents flicker.", "Reactive", {})
PARAM(react_check_every, Int, 5, "React mode: run the O N^2 hysteresis bond scan every N energy calls. 0 = displacement-triggered only.", "Reactive", {})
PARAM(react_check_disp_bohr, Double, 0.25, "React mode: also run the bond scan when any atom moved more than this distance in Bohr since the last scan. 0 disables the displacement trigger.", "Reactive", {})
PARAM(react_refractory_scans, Int, 10, "React mode: a pair whose bond just broke may not re-form for this many scans. Interrupts the form/break cycle that otherwise pumps the recombination energy through the thermostat over and over. 0 disables.", "Reactive", {})
PARAM(react_valence_cap, Bool, true, "React mode: refuse a new bond while an atom already uses its element valence plus one exchange slack, counting bond orders so multiple bonds consume valence. Prevents unphysical agglomerates; disable to sample unconstrained formation. Refused formations are logged at verbosity 2.", "Reactive", {})
PARAM(react_exchange_scans, Int, 20, "React mode: an atom may stay above its nominal valence for at most this many scans, then its weakest bond is broken. Forces exchange intermediates like a hydrogen bridging two heavy atoms to resolve instead of staying geometrically locked. 0 disables.", "Reactive", {})
PARAM(react_slack_form_factor, Double, 1.2, "React mode: tighter formation radius factor for bonds that push an atom above its nominal sigma valence into the exchange slack. A genuine exchange intermediate has the extra partner near bond distance; the ordinary optimistic factor would re-create bridges endlessly.", "Reactive", {})
PARAM(rev_enabled, Bool, false, "rev-gfnff stage 1: continuous bond order on every bonded term, blended bonded/non-bonded repulsion and an over-coordination energy (set by -method revgfnff).", "Reactive", {})
PARAM(rev_bond_weight, Bool, true, "rev-gfnff: multiply the bond well by the continuous bond order b_ij(r).", "Reactive", {})
PARAM(rev_term_weights, Bool, true, "rev-gfnff: multiply angle, torsion and inversion terms by the product of their bond weights.", "Reactive", {})
PARAM(rev_blend_repulsion, Bool, true, "rev-gfnff: evaluate the repulsion of every pair as b E_bonded + (1-b) E_nonbonded.", "Reactive", {})
PARAM(rev_over_coord, Bool, true, "rev-gfnff: add the over-coordination energy p_Z softplus(sum_j b_ij BO_ij - Val_Z)^2 (replaces the react valence cap).", "Reactive", {})
PARAM(rev_h_not_sp, Bool, true, "rev-gfnff stage 3a(ii): an sp hydrogen is not a BRIDGING atom (the reference scales a bridging H bond to 0.30 of its strength, which leaves the H-H wells of an exchange transition state 3.3x too shallow). rev-only.", "Reactive", {})
PARAM(rev_valence_share, Bool, true, "rev-gfnff stage 3a(ii): share the bond well over the valence the two ends still have free, c_ij = 1/2 (f_i + f_j) with f_i = clip((Val_i - sum_{k != j} w_ik)/w_ij, 0, 1), w the pair's term weight and Val_i = Val_Z(i) + softplus(sum_k shareClip(b_ik) - Val_Z(i)) the hypervalent-correct effective valence (a nominal valence plus a saturating settled-partner count). Two wells that share one valence then add to one well's worth (an exchange transition state) instead of two; a lone bond, a saturated equilibrium and a genuine hypervalent atom (ammonium, hydronium, perchlorate) all read c = 1 exactly.", "Reactive", {})
PARAM(rev_share_onethree, Bool, false, "rev-gfnff stage 3a(ii): the SMOOTH 1,3 proxy of the valence share - a pair claims no valence for the bond order that leaks onto it from a SETTLED shared partner, g_p = clip(1 - sum_k sigma_ik sigma_jk) with sigma = the settled weight of the tight bond order, and its own well is multiplied by c = 1 - g (1 - (f_i + f_j)/2) instead of (f_i + f_j)/2. Continuous by construction (a function of the existing bond orders, no threshold and no count), so it does not introduce the topology switch that a literal 1,3 test is. Its purpose: a compact polyhedron's perception carries 1,3 contacts as bonds whose WIDE term weight reads ~1, and without the proxy they claim the full valence of both ends (BF4- at B-F = 1.143 A, +569.7 kcal/mol), while the migrating pair of an exchange transition state, geometrically indistinguishable, must claim it. rev-only. DEFAULT OFF (measured, test_cases/revgfnff/_log/PROXY_STATUS.md): it fixes the BF4- probe exactly and is FD-exact and bit-identical at every equilibrium, but it gives a 1,3 contact pair the full well of its pair where the plain share suppresses it, and hot react MD then shows max |dE_jump| 3331 vs 471 kJ/mol and 24 vs 3 events >= 50 kJ/mol on the 22-cell grid. Opt-in experiment for the c_ij design decision.", "Reactive", {})
PARAM(rev_bo_center, Double, 2.0, "rev-gfnff: switching radius factor f_b of the bond order, R = f_b (rcov_i + rcov_j) fat_i fat_j.", "Reactive", {})
PARAM(rev_bo_width, Double, -7.5, "rev-gfnff: steepness k of the erf bond-order switch b = 0.5 (1 + erf(k (r - R)/R)); negative so b -> 1 inside R.", "Reactive", {})
PARAM(rev_bo2_center, Double, 1.4, "rev-gfnff: switching radius factor of the BOND ORDER used by the over-coordination sum and the repulsion blend (1 at the bond, 0 at 1,3 and hydrogen-bond distances).", "Reactive", {})
PARAM(rev_bo2_width, Double, -6.0, "rev-gfnff: steepness of the bond-order switch.", "Reactive", {})
PARAM(rev_over_shift, Double, 0.5, "rev-gfnff: the penalty argument is bo_sum - valence - shift, so a saturated atom pays nothing.", "Reactive", {})
PARAM(rev_over_preset, String, "stage1a", "rev-gfnff: built-in set of the over-coordination parameters (per-element p_over, over_shift, sigma valences). stage1a (default) = the placeholder set: p_over 0.3 for every element, over_shift 0.5, nominal sigma valences. fit2026-09-12 = the Levenberg-Marquardt fit against reference class C plus the GMTKN55 barriers, mirroring test_cases/revgfnff/params/rev_over_fit_2026-09-12.json; it stays opt-in until stage 3 re-measures the barriers. An explicit rev_over_p / rev_over_shift and the rev section of a -gfnff.param_file / -gfnff.param_json override document win over the preset.", "Reactive", {})
PARAM(rev_blend, Bool, true, "rev-gfnff stage 1b: blend the bonded terms of the old and the new topology while the transition pair crosses its weight window instead of swapping them at the event.", "Reactive", {})
PARAM(rev_bo3_center, Double, 1.6, "rev-gfnff stage 1b: centre of the TRANSITION coordinate switch in units of the covalent sum (x fat_i fat_j). The re-parametrisation of the neighbours blends in over rev_tr_begin..rev_tr_end of this switch, i.e. between 1.63x (1,3 pairs read 0.02) and 1.31x (equilibrium bonds read above 0.95); the pair's own well keeps the wide weight.", "Reactive", {})
PARAM(rev_bo3_width, Double, -8.0, "rev-gfnff stage 1b: steepness k of the transition coordinate switch (negative: 1 inside the centre).", "Reactive", {})
PARAM(rev_bo4_center, Double, 1.7, "rev-gfnff: centre of the REPULSION BLEND switch (b E_bonded + (1-b) E_nonbonded) in units of the covalent sum; 1 at a bond (0.99992 at 1.1x), ~0 at 1,3 and H-bond distances. Separate from the E_over switch, which a fit may soften.", "Reactive", {})
PARAM(rev_bo4_width, Double, -12.0, "rev-gfnff: steepness k of the repulsion blend switch.", "Reactive", {})
PARAM(rev_bo5_center, Double, 1.3, "rev-gfnff: centre of the repulsion blend switch of NON-BONDED-list pairs (units of the covalent sum): ~0 at 1,4 and H-bond distances, so nothing leaks at equilibrium; bonded-list pairs use rev_bo4_*.", "Reactive", {})
PARAM(rev_bo5_width, Double, -12.0, "rev-gfnff: steepness of the non-bonded repulsion blend switch.", "Reactive", {})
PARAM(rev_tr_begin, Double, 0.02, "rev-gfnff stage 1b: transition coordinate at which a forming pair starts its re-parametrisation (s = 0) and at which a breaking pair has finished it (s = 1); also the formation threshold of a fading well or of a 1,3 pair.", "Reactive", {})
PARAM(rev_tr_end, Double, 0.8, "rev-gfnff stage 1b: transition coordinate at which a formation is complete (s = 1) and below which a bond starts breaking (s = 0).", "Reactive", {})
PARAM(rev_tr_revert, Double, 0.75, "rev-gfnff stage 1b: a breaking transition still at s = 0 is undone once the coordinate climbs back above this value (hysteresis against vibrational chatter).", "Reactive", {})
PARAM(rev_tr_prebreak, Double, 0.5, "rev-gfnff stage 1b: a topology bond starts its breaking transition (s = 0) once its transition coordinate falls below this value; lower than rev_tr_end so that a hot bond vibrating around the formation end does not chatter.", "Reactive", {})
PARAM(rev_demote_cooldown, Int, 0, "rev-gfnff stage 1b: a transition interrupted below s = 0.5 is demoted (old topology restored, jump s dE) and its pair may not start again for this many scans; above 0.5 it is promoted (jump (1 - s) dE).", "Reactive", {})
PARAM(rev_bo13_form, Double, 0.1, "rev-gfnff stage 1b: a topological 1,3 pair becomes a bond (ring closure) once its bond order on the E_over switch exceeds this value; its transition then runs on that switch from here to 0.9.", "Reactive", {})
PARAM(rev_bo13_ordinary_join, Bool, false, "rev-gfnff stage 3a(ii) diagnostic: let a topological 1,3 pair join the bond list on the ORDINARY formation criterion rev_bo2_form with the ordinary transition window, instead of the special 1,3 window rev_bo13_form .. 0.9. The special window exists to stop a geminal pair from stealing valence early; it is keyed on the discrete graph bit 'shares a SETTLED neighbour', which can flip while the pair is already inside its window and then applies the formation as a hard swap at s = 1. Set this true to measure how much of the hot-MD jump tail that classification carries. DEFAULT OFF: the delivered scan behaviour is unchanged.", "Reactive", {})
PARAM(rev_well_form, String, "mg3", "rev-gfnff stage 3a(iii): the shape of the bond well. 'gauss' is the delivered GFN-FF form k_b exp(-alpha (r - r0)^2) multiplied by the reactive term weight - bit-identical to the state before this option. 'mg' and 'erfmorse' are the two curvature-pinned two-parameter forms of FABLE_REVIEW_2 B: both are E = -D (2y - y^2) with y = exp(-(a x + beta x^2)) resp. y = erfc((x - u)/sigma)/erfc(-u/sigma), x = r - r0, and both reproduce the delivered r_min AND the delivered force constant 2 alpha |k_b| BY CONSTRUCTION (a is a closed form, u a bisection), so only the depth scale s = D/|k_b| and the tail are fitted - per element pair, on the class-A r2SCAN-3c bond scans, against a charge-frozen rest (src/core/energy_calculators/ff_methods/rev_well_table.h, generated by scripts/revgfnff_wellfit.py). The delivered Gaussian is 19.2 kcal/mol RMS off the reference on the break side (median over 32 bond types); both new forms are 2.1. The term weight is NOT applied to these wells: they decay by themselves, so multiplying by it would truncate the tail that was just fitted. An element pair with no class-A data keeps the Gaussian and says so at verbosity 2. AI-FITTED, machine-tested only. DEFAULT mg from Sep 19 to Sep 22, 2026, when mg3 took over (see below); 'mg' remains available and reproduces that behaviour bit-for-bit. That first flip was made because mg and erfmorse are indistinguishable on the data and mg is the cheaper one (its a is a closed form where erfmorse needs a per-bond bisection, 1.42x the react-MD wall time before that bisection is cached). The flip buys class-A median rms 24.49 -> 19.50, median dev D_e -25.53 -> -12.74 kcal/mol and median dev r90 -0.330 -> -0.058 A, at an equilibrium bond-length shift of up to 0.0063 A and a pooled conformer/S66 guard MAD of 1.0341 -> 1.0439 kcal/mol. It does NOT touch plain -method gfnff (the whole rev path is gated on rev_enabled). The full benefit needs the stage-3b bond-order-resolved table - the delivered one is keyed on the ELEMENT PAIR, which is why the harness median lands at 19.50 and not at the 2.07 the per-system fit reaches. 'gauss' reproduces the pre-Sep-19 behaviour. OPT-IN, Sep 20, 2026: 'mg2' and 'mg3' are the same MG well with the curvature and the minimum FREED - four parameters per key (depth scale s, curvature scale ca, tail beta, and dr0, an offset added to the model's own dynamic r0), fitted the same way on the same class-A curves (rev_well_table_v2.h). ca = 1, dr0 = 0 is exactly 'mg'. 'mg2' is keyed on the element pair like 'mg'; 'mg3' is additionally keyed on the BOND ORDER and interpolates linearly between the single/double/triple fits of that pair at the pair's continuous order 1 + pibo * (both ends sp ? 2 : 1) - a topology constant inside one energy call, so it adds no geometry derivative and no switch. Neither changes what gauss|mg|erfmorse compute. DEFAULT mg3 since Sep 22, 2026 (operator decision, WORK_STATUS package 12): against mg it buys class-A harness median rms 22.15 -> 13.22, median dev D_e -14.32 -> -6.86 kcal/mol, median |b_model - b_r2SCAN-3c| over the 32 class-A bonds 0.0246 -> 0.0044 A, class D dE_MAD 5.152 -> 4.555 and grad_RMS 14.428 -> 11.313, rkt06 path rms 2.3788 -> 2.2665, at a pooled conformer/S66 guard MAD of 1.0439 -> 1.0547 kcal/mol (+0.011) and an equilibrium bond-length shift of 0.0212 A vs gauss (which moves TOWARDS the reference). The react-MD smoothness concern raised against mg3 in Sep 2026 was measured away in package 11: over 11700 paired-replicate trajectories (780 per arm and time step) neither mg2 nor mg3 is distinguishable from mg on the per-step tail at any time step. 'mg2' is the same accuracy family without the bond-order dimension and is equal-or-worse than mg3 on every measured row at the same guard cost. AI-FITTED, machine-tested only.", "Reactive", {})
PARAM(rev_well_order_override, Double, -1.0, "rev-gfnff stage 3b DIAGNOSTIC: force every bond's continuous bond order to this value instead of the one the FT-HMO pi order gives, so the bond-order-resolved well table (rev_well_form mg3) can be scanned in its own interpolation variable. Negative (the default) leaves the computed order alone and the run is bit-identical to not passing the flag. It exists because the order is a TOPOLOGY quantity and cannot otherwise be varied continuously from outside: with it, dE/d(order) is measurable by finite differences and the claim that the well responds to a drifting pi order proportionally (never in a step) is a measurement rather than an argument. Only mg3 reads it. Not for production use.", "Reactive", {})
PARAM(rev_share_form, String, "conserving", "rev-gfnff stage 3a(ii): which valence-share formula multiplies the bond well. 'delivered' is the left-over rule c_ij = 1/2 (f_i + f_j) with f_i = clip((Val_i - sum_{k != j} w_ik)/w_ij): every pair asks how much of atom i's valence the OTHER partners left over, so with a term weight of ~1 on every partner ONE partner too many makes every pair of that atom see nothing left and the atom hands out 0 of its valences - the measured cause of the artificial radical adducts (a free H that hits a saturated C, N or O lands 54-107 kcal/mol below the r2SCAN-3c reference). 'conserving' is FABLE_REVIEW_2 A.5: f_i = min(1, Val_i / S_i) per ATOM with S_i the atom's whole claim sum, c_ij = f_i f_j (the product, not the mean - with the mean the rkt06 exchange transition state gets c = 0.75 and the path breaks), so sum_j f_i w_ij = min(Val_i, S_i) exactly and an over-coordinated atom SPREADS its valence instead of forfeiting it. In that mode the excess budget is also granted by CHARGE instead of by element: Val_i = Val_Z + min(G(S_i - Val_Z), X_i) with X_i = 0 for H and F, 1 for group 13, 6 - Val_Z for period >= 3 groups 15-17, and the clipped positive topological charge of the atom plus its H partners otherwise - what separates NH4+ from NH3 + H is the electron count, and the charge is the only electron-count information a force field has, plus the donor rule of -gfnff.rev_share_donor_rule. DEFAULT conserving since Sep 19, 2026 (operator decision): with the donor rule the dative/ylide regression that held the flip back is 0.00 kcal/mol, and the mode removes the artificial radical adducts (class-C dev min -87..-107 -> -1.5..+0.0 kcal/mol) and takes the 20-cell react-MD grid from 487 to 72 per-step events >= 50 kJ/mol. 'delivered' reproduces the pre-Sep-19 behaviour.", "Reactive", {})
PARAM(rev_share_min_width, Double, 0.1, "rev-gfnff stage 3a(ii), rev_share_form conserving only: half-width (in units of Val/S) of the C1 smooth min that replaces min(1, x). The min is EXACTLY 1 for x >= 1 - so a saturated equilibrium atom keeps a share factor of literally 1.0 and its term stays bit-identical - and EXACTLY x for x <= 1 - width, so valence conservation is exact where the share bites; only the corner between the two is rounded by a cubic. Smaller values follow min(1, x) more closely at the price of a stiffer force through the join.", "Reactive", {})
PARAM(rev_share_donor_rule, Bool, true, "rev-gfnff stage 3a(ii), rev_share_form conserving only: grant the excess budget X_i >= 1 to an atom that DONATES into a dative bond. The charge rule cannot see one - a dative bond puts a whole valence into the acceptor's empty orbital while the donor's EEQ charge stays near +0.2 - so without this rule an amine borane's nitrogen gets Val ~ 3.2 against S ~ 4 and all four of its wells are scaled by ~0.8, measured as +73 to +110 kcal/mol on H3N-BH3, H3N-O and H3N-CH2 where the delivered share is inert. An atom qualifies if, in the corner being evaluated, it has a partner that is either a group-13 element (B, Al: an empty p orbital, which no bond count can reveal) or an atom carrying FEWER partners than its own nominal sigma valence, i.e. a free coordination site (the amine oxide's one-coordinate O, the ylide's three-coordinate C). The grant is a maximum against the charge rule, never a replacement, and it applies only where the charge rule decides - group 13 and the period >= 3 octet expansion already carry a larger cap. Both tests read the corner's own bond list, so the cap stays a per-corner constant with no chain rule. Costs nothing at equilibrium (the cap is only reachable once an atom claims MORE than its nominal valence) and nothing on a radical adduct (a one-coordinate H is not deficient). DEFAULT ON; false is the ablation arm.", "Reactive", {})
PARAM(rev_budget_fix_h, Bool, true, "rev-gfnff stage 3a(ii): hydrogen keeps its nominal valence 1 in the share budget - never hypervalent. The alternative rule (false) gives every atom Val_i = Val_Z + softplus(settled count - Val_Z), so an H that sits between two partners (a just-formed H2 still bonded to its carbon) gets Val_H -> 2 as soon as the H-H tight bond order crosses the settled window, and BOTH its wells jump from half share to full share in one step - the measured origin of the hot react-MD blow-ups. A hydrogen has ONE valence; a bridging H is a 3c-2e bond whose two partial wells must SHARE it. With the flag on, Val_H = 1 exactly and its derivative channel is zero; every other element keeps the softplus budget. DEFAULT ON since Sep 18, 2026: on the 20-cell react-MD grid it removes both runaways (max per-step dEpot 2593.6 -> 59.3 kJ/mol, T_max 1.6e8 -> 8306 K, hard swaps 5 -> 0) and no falsifier moves (equilibria, hypervalent ions, BF4-, rkt06 all bit-identical). false reproduces the pre-Sep-18 behaviour.", "Reactive", {})
PARAM(rev_pair_validity, Bool, false, "rev-gfnff pair-validity gate (FABLE_BOND_STATE_2.md sec 2.1-rev): a per-corner, geometry-free veto on a listed bonded pair (i,j) that is a closed-shell repulsion rather than a bond - true unless one end has a free valence slot (n_other < Val_Z + cap, the SAME per-atom budget cap the conserving valence share computes, FFWorkspace::shareCapForAtom), or is a metal, or is a proton bridging two lone-pair atoms (an X-H-Y 3c-4e bridge), or the pair is a doubly-bridged M...M diagonal alongside >= 2 pure H bridges, or the pair shares a metal/deficient neighbour (a sigma complex, e.g. Kubas eta2-H2), or the pair's local shell carries >= 0.5 e of Phase-1 topological charge (a cationic 3c-2e bond, H3+/CH5+). An invalid pair is removed from that corner's own topology and the corner is REGENERATED from the reduced bond list (hybridisation, angles, torsions, pi-systems, Phase-1 EEQ all freshly derived, exactly as if the pair had never been perceived) - not merely zeroed in the bond energy, since the reference test (a compressed BF4- probe) requires the reduced corner's hybridisation and EEQ to match a naturally 4-bonded topology bit-for-bit. Fixes the compressed-BF4- artefact (a fresh 10-bond perception, six spurious F...F contacts alongside the four genuine B-F bonds, scores hundreds of kcal/mol above the same force field's own 4-bond evaluation) and the geminal H...H well of a hot react-MD trajectory (an artefact bond between two already-saturated hydrogens of the same carbon, which re-parametrises the sibling C-H bonds). DEFAULT OFF: a rev-gfnff research mechanism, not a port-fidelity change, bit-identical when off (verified over GMTKN55+MOR41+S30L-CI, 2647 structures, 0 moved) and on every already-valid reference structure with it on. See docs/REV_GFNFF_STAGE3A.md and test_cases/revgfnff/_log/PAIR_VALIDITY_IMPL_STATUS.md.", "Reactive", {})
PARAM(rev_h_scope, Bool, false, "rev-gfnff Q5 (FABLE_BOND_STATE_2.md sec 7.6, 'H-scope'): master switch for the hydrogen-perception rule set. A bridging hydrogen (an X-H-Y 3c-4e/3c-2e bond, a migrating H) is given hyb=1 ('sp') by the topology perception today, which lets it trigger four rules meant for a genuine sp centre: the bsmat[1][*] bond-strength column instead of the terminal-H bsmat[hyb_X][0] (H1), 3-ring/fxh membership through the bridge triangle (R1), an angle centred on the bridging H with theta0=180 (H2), and counting as an sp/sp2 picon neighbour for pi-conjugation (P1, a structural consequence of H1 with no separate switch). All four are per-corner constants (element + partner count + the corner's own ring/pi membership), so a change is carried by the existing s-blend with no new derivative. Subsumes rev_h_not_sp for real hydrogen (once hyb(H)=0 always, rev_h_not_sp's is_bridge/0.30-scaling test becomes structurally unreachable for Z==1); rev_h_not_sp itself is kept, unchanged, as the narrower mechanism reachable when rev_h_scope_h1 is off. Sub-switches rev_h_scope_h1/h2/r1 are ablation arms, read only when this is true. DEFAULT OFF: bit-identical over GMTKN55+MOR41+S30L-CI (2647 structures) when off; touches exactly the 80 structures with a genuinely 2+-coordinate hydrogen when on, none other. See docs/REV_GFNFF_STAGE3A.md and test_cases/revgfnff/_log/H_SCOPE_IMPL_STATUS.md.", "Reactive", {})
PARAM(rev_h_scope_h1, Bool, true, "rev-gfnff Q5 'H1' (only read when rev_h_scope is true): every atom with Z==1 gets hyb=0, whatever its partner count (tests the ELEMENT, not the periodic group - the grp==1 branch that also covers Li/Na/K stays reachable for those). X-H bond strength = bsmat[hyb_X][0], the same value a terminal H on that X gets; H-H is unchanged (already a pure Z==1&&Z==1 special case, bstren[1]=1.00).", "Reactive", {})
PARAM(rev_h_scope_h2, Bool, true, "rev-gfnff Q5 'H2' (only read when rev_h_scope is true): no GFN-FF angle term is ever centred on a Z==1 atom. Needed together with H1 - H1 alone gives a spurious tetrahedral theta0=109.5 deg at a bridging H (tried and rejected once before this Q5 work, see the note in determineHybridizationFortran).", "Reactive", {})
PARAM(rev_h_scope_r1, Bool, true, "rev-gfnff Q5 'R1' (only read when rev_h_scope is true): ring enumeration excludes every Z==1 atom entirely - no ringf and no 3-ring fxh correction reaches any bond of a hydrogen-bridged ring, including the bridging bonds themselves (fxh is keyed on the heavy atom's own ring membership, so excluding H from ring perception clears it for every bond of that heavy atom at once, bridging and sibling alike).", "Reactive", {})
PARAM(rev_max_transitions, Int, 4, "rev-gfnff stage 1b: transitions blended at the same time (2^k topology corners are evaluated per step); a further event snaps the transition closest to either end of its window.", "Reactive", {})
PARAM(rev_bo_form, Double, 0.05, "rev-gfnff react scan: a non-bonded pair (never a 1,3 pair) joins the bond list once its term WEIGHT exceeds this value, i.e. where its terms are still ~0.", "Reactive", {})
PARAM(rev_bo_break, Double, 0.02, "rev-gfnff react scan: a bond leaves the list once its term weight falls below this value.", "Reactive", {})
PARAM(rev_form_switch, String, "order", "rev-gfnff react scan: which switch decides that a non-bonded pair BECOMES a bond. order (default, Sep 12, 2026) = the narrow bond-order switch rev_bo2_* at the threshold rev_bo2_form, i.e. the switch that already reads ~0 at 1,3 and hydrogen-bond distances; the new pair's own well then enters through the transition blend (s) instead of at once, so the join stays energy-neutral. weight = the previous behaviour: the WIDE term-weight switch rev_bo_* at rev_bo_form, which still reads above 0.05 out to 2.31x the covalent sum and therefore joins hydrogen bonds and van-der-Waals contacts (water dimer: 18 formations/ps, Epot 4 kcal/mol below the static run).", "Reactive", {})
PARAM(rev_bo2_form, Double, 0.1, "rev-gfnff react scan: formation threshold on the NARROW bond order (rev_bo2_*) for rev_form_switch = order. 0.1 crosses at 1.611x the covalent sum, which is the non-rev react formation factor (1.6) and the radius at which a bond starts breaking (transition coordinate 0.5 at 1.600x), so formation and break are symmetric; the water dimer reads 3.9e-4 at its H...O and 9.7e-7 at its O...O contact.", "Reactive", {})
PARAM(rev_over_p, Double, 0.3, "rev-gfnff: default over-coordination prefactor p_Z in Eh (per-element values via the rev.p_over override).", "Reactive", {})
PARAM(rev_over_k, Double, 10.0, "rev-gfnff: softplus steepness of the over-coordination penalty.", "Reactive", {})
// ---- rev-gfnff stage 2 (Claude Generated, Sep 2026): split-charge (SQE) model ---------------
// docs/REV_GFNFF_STAGE2.md. q_i = q0_i + sum_j p_ij with a per-bond hardness
// kappa_ij = kappa0_ij / b_ij; kappa0_ij = 1/2 (kappa_Z(i) + kappa_Z(j)). kappa -> 0 on a
// connected bond graph reproduces the constrained EEQ minimum exactly, b -> 0 pins the charge
// on the separating fragment, so the model interpolates continuously between the two limits
// that today's discrete fragment count brackets.
PARAM(rev_charge_model, String, "eeq", "rev-gfnff stage 2: charge model. eeq = today's per-fragment constrained EEQ; sqe = split charges on the bond graph with a per-bond hardness kappa_ij = kappa0_ij / b_ij (no fragment constraint). Default eeq until the kappa_Z fit is done.", "Reactive", {})
PARAM(rev_sqe_kappa_H, Double, 0.0, "rev-gfnff stage 2: bond-hardness parameter kappa_Z of hydrogen (Eh). 0 = the pair costs nothing to polarise, i.e. the connected-graph EEQ limit.", "Reactive", {})
PARAM(rev_sqe_kappa_C, Double, 0.0, "rev-gfnff stage 2: bond-hardness parameter kappa_Z of carbon (Eh).", "Reactive", {})
PARAM(rev_sqe_kappa_N, Double, 0.0, "rev-gfnff stage 2: bond-hardness parameter kappa_Z of nitrogen (Eh).", "Reactive", {})
PARAM(rev_sqe_kappa_O, Double, 0.0, "rev-gfnff stage 2: bond-hardness parameter kappa_Z of oxygen (Eh).", "Reactive", {})
PARAM(rev_sqe_kappa_F, Double, 0.0, "rev-gfnff stage 2: bond-hardness parameter kappa_Z of fluorine (Eh).", "Reactive", {})
PARAM(rev_sqe_kappa_Cl, Double, 0.0, "rev-gfnff stage 2: bond-hardness parameter kappa_Z of chlorine (Eh).", "Reactive", {})
PARAM(rev_sqe_bmin, Double, 1e-3, "rev-gfnff stage 2: a pair whose bond order falls below this is dropped from the split-charge system (its hardness kappa0/b would exceed 1000 kappa0, i.e. it is rigid).", "Reactive", {})
PARAM(rev_sqe_kappa_form, String, "inverse", "rev-gfnff stage 2 B2: functional form of the bond hardness kappa(b). inverse = kappa0/b (the original, factor 1.6 only between b=0.97 and b=0.6); power = kappa0/b^n with n = rev_sqe_kappa_exponent (steeper, still finite at b=1); vanishing = kappa0 (1-b)/b (exactly 0 at b=1, so only a weakened bond is penalised). Default inverse keeps every existing result bit-identical.", "Reactive", {})
PARAM(rev_sqe_kappa_exponent, Double, 3.0, "rev-gfnff stage 2 B2: exponent n of the power form of kappa(b) = kappa0/b^n. Ignored by the other forms.", "Reactive", {})
PARAM(rev_sqe_q0_rule, String, "mu", "rev-gfnff stage 2 B2: where a charged fragment's integer charge sits at p = 0 (the static/cold-start reference charges q0). uniform = spread flat over the fragment (the original rule; for a symmetric anion that IS already the EEQ minimum, so kappa has no lever); mu = localised on the fragment atoms with the lowest EEQ chemical potential for an electron, highest for a hole. At kappa = 0 the two are provably identical (the increments p reach the same EEQ minimum from either start).", "Reactive", {})
PARAM(rev_sqe_q0_mu_tau, Double, 1.0, "rev-gfnff stage 2, rev_sqe_q0_rule mu: temperature in kcal/mol (per elementary charge) of the smooth placement. Instead of the hard argmin over the EEQ chemical potential, the energy is the Boltzmann average sum_p w_p E_p over the whole-unit placements p, w_p ~ exp(mu.q0_p/tau), each E_p a full split-charge solve; the analytic gradient includes dw_p/dx. Placements more than 34 tau above the best carry no weight, so away from a mu crossing the result is the hard rule to the last bit, and at an exact symmetric tie (Cl2-, a C2v carboxylate) the equal branches give the hard energy too; only the force cusp at a crossing is removed. Costs one extra split-charge solve per near-tied placement. 0 = the original hard rule (kept only to reproduce old numbers). See test_cases/revgfnff/_log/MU_CUSP_STATUS.md.", "Reactive", {})
PARAM(rev_sqe_base_q0_keep, Bool, true, "rev-gfnff stage 2 react mode (Sep 29, 2026, needs rev_charge_model sqe): when the first transition of a set starts, the old-topology base corner - which carries the full blend weight at s = 0 - keeps the q0 the slot corner used a moment before (the frozen corner q0, the P2 Phase-1 placement, or the q0 rule's placement), instead of re-deriving it by rounding the converged charges (and, in P3 flat mode, re-localising it). Re-deriving took a discrete decision at full weight and made the energy jump at every transition start when any pair hardness is > 0 (formate C-H at kappa 0.5: 39.9 kcal/mol; Cl2- at kappa_Cl 0.85: 9.25). Inert when every pair hardness is 0 (kappa 0, harris). false = the old capture, kept only to reproduce old numbers. See test_cases/revgfnff/_log/MU_CUSP_STATUS.md section 10.", "Reactive", {})
PARAM(rev_sqe_phase1, Bool, false, "rev-gfnff P2 (Sep 2026, opt-in, needs rev_charge_model sqe): solve the Phase-1 topology charges qa with the same split-charge model as the Phase-2 charges (pairs = the topology bonds at their TOPOLOGICAL bond order 1, same kappa_Z and q0 rule), so the charge-dependent hardness dgam(qa)/alpeeq(qa) of the Coulomb self-energy localises together with the charge instead of staying at the delocalised constrained-EEQ qa. qa stays a function of the topology alone. At kappa = 0 on a connected graph identical to the constrained Phase 1. See test_cases/revgfnff/_log/P2P3_STATUS.md.", "Reactive", {})
PARAM(rev_sqe_group_pairs_only, Bool, false, "rev-gfnff stage 2, Sep 25, 2026, opt-in, needs rev_charge_model sqe, use together with rev_sqe_virtual_pairs: the Phase-2 split-charge solve drops every pair whose two atoms lie in DIFFERENT EEQ constraint groups - the Phase-2 analogue of the P2 Phase-1 restriction. Such a pair is a pass-2 bond between two pass-1 fragments, Known Issue 17, e.g. the C-X bond of an [X-CH3-X]- SN2 transition state or an X2- between the pass-1 split and the static bond cutoff; without the flag it moves charge across the group border, so SQE at kappa 0 differs from the constrained EEQ by up to 200 kcal/mol. With both flags SQE at kappa 0 reproduces the constrained EEQ exactly. CHANGES the recommended X2- setting in that distance band, see test_cases/revgfnff/_log/SQE_INVARIANT_STATUS.md.", "Reactive", {})
PARAM(rev_sqe_virtual_pairs, Bool, false, "rev-gfnff stage 2 (Sep 24, 2026, opt-in, needs rev_charge_model sqe): the Phase-2 split-charge solve chains every bond-graph component that shares one EEQ constraint group with zero-hardness VIRTUAL pairs, the Phase-2 analogue of the P2 Phase-1 virtual pairs. Without it charge cannot move between two unbonded atoms of the same constraint group and stays at the integer q0 placement - which the frag_charge_model ensemble merged corner creates just past the bond cutoff (X2- label dependence, SQE(kappa 0) != EEQ by up to 108 kcal/mol). At kappa = 0 the reachable charge space is then the constrained EEQ one. See test_cases/revgfnff/_log/X2_SCOPE_STATUS.md.", "Reactive", {})
PARAM(rev_excess_electron, Bool, false, "rev-gfnff P3 (Sep 2026, opt-in, needs rev_charge_model sqe): perceive excess electrons with no bonding slot left (an anionic fragment whose atoms are valence-saturated, e.g. Cl2-, F2-: x = max(0, -Q_f - free slots)), spread x over the fragment's bonds that have a calibrated half-order well row, lower their bond order by x/2 (the mg3 well table then reads the half-order row) and add x*rev_excess_kappa to their split-charge hardness (the resonance of the 2c-3e bond lives in the well, not in the Coulomb term). A per-topology constant; inert for every neutral fragment. See test_cases/revgfnff/_log/P2P3_STATUS.md.", "Reactive", {})
PARAM(rev_pi_excess_electron, Bool, false, "rev-gfnff P3 pi* extension (Sep 25, 2026, opt-in PROTOTYPE, needs rev_excess_electron): perceive the pi* excess electron of a DIATOMIC radical anion (O2-, S2-), which the sigma-slot budget of rev_excess_electron cannot see. Fires only for an isolated two-atom bond-graph component with a pi component (continuous order > 1) and fragment charge exactly -1, on a pair with a calibrated pi-excess well row (rev_well_table_v2.h kPiExcessEntries); the bond's mg3 well is then replaced by that row and the pair's split-charge excess (flat/frac/harris) counts y = 1. No polyatomic pi system, no dianion, no aromatic. See test_cases/revgfnff/_log/PI_STAR_STATUS.md.", "Reactive", {})
PARAM(rev_excess_bond_extend, Double, 1.0, "rev-gfnff P3 (Sep 27, 2026, opt-in, needs rev_excess_electron; the pi* pairs also rev_pi_excess_electron): factor on the ordinary getnb bond threshold for a 2c-3e candidate pair, so the X-X bond - and with it the excess count x, the half-order / pi-excess well and the harris x*g term - survives past the ordinary cutoff (~1.05 r_min) out to where the reference binding has gone. A candidate is a pair of atoms that the ordinary perception leaves WITHOUT ANY bond, whose element pair has a half-order (or, with rev_pi_excess_electron, a pi-excess) well row, in a system of negative net charge, and that are each other's only such candidate. Nothing else is ever bonded by it. The frag_charge_model ensemble window moves with it ([f, f*frag_charge_s_max] times the ordinary threshold). 1.0 = off (bit-identical). See test_cases/revgfnff/_log/X2_COMPRESSED_SURVEY_STATUS.md section 8.", "Reactive", {})
PARAM(rev_excess_kappa, Double, 100.0, "rev-gfnff P3: flat split-charge hardness (Eh) per excess electron on a perceived 2c-3e pair (x*rev_excess_kappa, no bond-order dependence). Large = the charge stays on its q0 atom, i.e. the Coulomb term carries no delocalisation energy for that pair.", "Reactive", {})
PARAM(rev_excess_mode, String, "flat", "rev-gfnff P3 (Sep 2026, opt-in alternative, see test_cases/revgfnff/_log/P2P3_ALTERNATIVES_STATUS.md): how the Coulomb term is kept from double-counting a perceived 2c-3e pair. flat = the shipped x*rev_excess_kappa split-charge hardness (localises the charge on its q0 atom); frac = fractional-charge correction E_x = 1/2 x c K_ij(r) q_i q_j with K_ij the pair EEQ curvature and c = rev_excess_frac_c: removes the fraction c of the pair delocalisation energy while the charges stay symmetric, polarisable and independent of q0; harris = no extra hardness (charges free, as flat at rev_excess_kappa 0) plus a non-self-consistent energy x*g(r) per perceived pair, g fitted to DLPNO-CCSD(T) (rev_harris_table.h, test_cases/revgfnff/_log/P2P3_HARRIS_STATUS.md).", "Reactive", {})
PARAM(rev_excess_frac_c, Double, 0.9, "rev-gfnff P3, rev_excess_mode frac: fraction c in [0, 0.99] of the pair EEQ curvature removed (1 - c of the EEQ delocalisation energy is kept; the axial polarisability scales as 1/(1 - c)).", "Reactive", {})
PARAM(rev_excess_react_consistent, Bool, true, "rev-gfnff P3 flat mode (Sep 2026, repair, default on; only acts when rev_excess_electron is on AND in react-mode topology corners, see test_cases/revgfnff/_log/P2P3_ALTERNATIVES_STATUS.md): keep the flat excess-electron hardness consistent across react-mode topology corners. (a) a corner that perceives an excess-electron pair gets its q0 re-localised on that pair (on the atom its Phase-1 qa holds the charge on, then lowest EEQ mu, then index), because the flat hardness freezes the charge at q0 and a delocalised corner q0 would carry the full EEQ delocalisation energy; (b) a pair in flight that is not a bond in one corner inherits the largest x*kappa_x any corner assigns it (the corner that has the bond decides).", "Reactive", {})
PARAM(storsion_reference_loop_bug, Bool, false, "Reproduce the reference implementation's triple-bond-torsion (sTors) loop bug bit-for-bit. Both pprcht/gfnff and xtb 6.7.1 call sTors_eg(m,...) with the array SIZE m instead of the loop index, so they evaluate only the LAST detected C-triplebond-C torsion, m times, and drop all others (and give exactly zero whenever the last slot was never filled). Curcuma sums every detected torsion, which is what the term is meant to do - its erefhalf is a DLPNO-CCSD(T) diphenylacetylene reference value, not a fitted parameter. Enable only to reproduce reference totals exactly.", "Advanced", {})
PARAM(param_file, String, "", "rev-gfnff: JSON file with sparse parameter overrides ({gen:{...}, tables:{name:{Z:value}}, rev:{...}}) deep-merged over the built-in GFN-FF tables; unknown keys abort.", "Advanced", {})
PARAM(param_json, String, "", "rev-gfnff: the same override document given inline as a JSON string; merged after param_file.", "Advanced", {})
PARAM(dump_params, String, "", "rev-gfnff: write the generated per-term force-field parameters (bonds fc/r0/alpha, angles, torsions, ...) of the current molecule to this JSON file after initialisation.", "Advanced", {})
PARAM(hh_repulsion_bpair, Bool, true, "Take the 1,3/1,4 classification of the non-bonded H...H repulsion factors hh13rep/hh14rep from the reference pair table topo%bpair, i.e. nbondmat, as in gfnff_ini.f90:755-756. false restores the former plain BFS bond count, which differs across neighbour entries stored on one side only - eta bonds, main-group metals - and there depends on the atom order of the input. Default since Sep 2026.", "Advanced", {})
// OPERATOR NOTE (Sep 25, 2026): default false adopted from feature/multi-gpu; RE-EVALUATE when rev-gfnff stage 3 (the H/C/N/O/F/Cl element-table refit, docs/REV_GFNFF_ROADMAP.md WP5) happens - not settled permanently.
PARAM(dispersion_atm, Bool, false, "Add a D4 Axilrod-Teller-Muto three-body dispersion term over BONDED i-j-k triples only, s9=1, zero damping. The GFN-FF reference has no such term - its dispersion is pairwise, gfnff_gdisp0.f90 d3_gradient - and bonded triples are the ones the damping suppresses, so the term is at most 0.025 kcal/mol on MOR41+GMTKN55. Default off since Sep 2026; true reproduces the earlier curcuma totals.", "Advanced", {})
PARAM(amideh_acidity_order_bug, Bool, false, "Reproduce the reference atom-order bug in the amide-H donor acidity: gfnff_ini.f90:793-798 resets hbaci of atom i and scales hbaci of its first neighbour in the same loop, so the 0.8 amide factor on the nitrogen survives only when the nitrogen precedes its H in the input. The H-bond then depends on the atom numbering and is 1.25x stronger otherwise. curcuma applies the factor regardless of order; enable only to reproduce reference totals for arbitrary atom orders.", "Advanced", {})
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
        /// rev-gfnff P3 (Claude Generated, Sep 23, 2026): perceived excess electrons per bond,
        /// keyed on (min, max) atom index; only pairs with a nonzero x are stored. Empty unless
        /// rev_excess_electron is on. See GFNFF::revExcessElectrons().
        std::map<std::pair<int, int>, double> rev_excess;
        /// rev-gfnff P3 pi* prototype (Claude Generated, Sep 25, 2026): perceived pi* excess
        /// electrons y of a diatomic radical anion, same key; empty unless rev_pi_excess_electron.
        /// See GFNFF::revPiExcessElectrons() and _log/PI_STAR_STATUS.md.
        std::map<std::pair<int, int>, double> rev_pi_excess;
        /// rev-gfnff P2: the q0 the Phase-1 split-charge solve used (empty unless
        /// rev_sqe_phase1). The static Phase-2 q0 reuses it, so both phases localise a charge
        /// on the SAME atom (their mu probes use different matrices and could disagree).
        Vector rev_sqe_q0;
        /// rev-gfnff P2: the hybridisation Phase 1C computed dgam with (the GEODEP sp2 -> sp3
        /// promotion after the Hueckel section changes topo.hybridization later, so re-deriving
        /// dgam from the final array would change dgam even where qa did not change).
        std::vector<int> rev_hyb_eeq;
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
     *
     * Claude Generated (Sep 2026): now a thin gate wrapper (gfnff_pair_validity.cpp) around
     * generateGFNFFParameterSetImpl() - see that function's declaration below. Off by default
     * (-gfnff.rev_pair_validity), in which case it is exactly the one call it always was.
     */
    GFNFFParameterSet generateGFNFFParameterSet();

    /// Claude Generated (Sep 2026): the actual generator, renamed out of the way of the gate
    /// wrapper above. Every existing call site keeps calling generateGFNFFParameterSet(); this
    /// is what that now calls (once, or twice under the pair-validity gate).
    GFNFFParameterSet generateGFNFFParameterSetImpl();

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

    // === Reused-topology invalidation (Claude Generated Sep 2026) ===
    // A force-field topology is perceived once (initializeForceField ->
    // generateGFNFFParameterSet) and then reused. needsFullTopologyUpdate() + the geometry
    // tracker already refresh the *TopologyInfo* on a >0.5 Bohr displacement, but that never
    // touched the interaction lists the energy is actually built from, so a calculator reused
    // across structurally different frames (-batch_reuse_topology true) silently ran every
    // frame on the first frame's bond graph (18-117 kcal/mol at the same geometry;
    // OUTLIER_STATUS.md section F). Calculation() now compares the perceived bond graph with
    // the one the lists were built for and rebuilds when they differ; the PARAM
    // reuse_topology_check turns that off.

    /// Rebuild the whole force-field interaction list from the current geometry/topology.
    /// Returns false (and leaves the previous lists in place) if parameter generation failed.
    bool rebuildForceFieldForCurrentGeometry();

    /// Geometric bond perception on the CURRENT geometry, uncached. The geometric half of
    /// getCachedBondList(), shared so the reuse check needs no second implementation.
    std::vector<std::pair<int,int>> perceiveGeometricBonds() const;
    /// rev_excess_bond_extend (Claude Generated, Sep 27, 2026): true if (i, j) is an eligible
    /// 2c-3e element pair in a negatively charged system (flags on, rows present). Distance-free.
    bool revX2PairElements(int i, int j) const;
    /// rev_excess_bond_extend: the extra bonds of the 2c-3e candidate rule on top of \p bonds
    /// (the ordinary list of the same pass); \p qshift = the pass's charge shrink per atom.
    std::vector<std::pair<int,int>> revX2ExtendBonds(const std::vector<std::pair<int,int>>& bonds,
        const std::vector<double>& qshift, const std::vector<double>& fm_atom) const;
    /// ordinary getnb threshold (Bohr) of pass 1 (qa = 0) - the isolation test of the rule above
    double getnbThresholdPass1(int i, int j) const;
    /// rev_excess_bond_extend: the full, distance-free eligibility of (i, j) - element pair, net
    /// charge -1, both atoms without any ordinary pass-1 bond, no third atom inside the pair's
    /// bond ellipsoid. One function for the perception AND the ensemble window.
    bool revX2PairExtendable(int i, int j) const;

    /// Canonical i<j bond graph of a neighbour list, for the reuse comparison.
    static std::vector<std::pair<int, int>>
    canonicalBondGraph(const std::vector<std::vector<int>>& neighbor_lists);

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
    /// Number of events not yet consumed (diagnostics). Claude Generated (Sep 2026).
    size_t pendingReactEvents() const { return m_react_events.size(); }
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
     * @param bonds   Static bond list from GFNFFParameterSet (m_gpu_params_leaked->bonds)
     * @return BondHBRebuildResult with updated bond_hb_data, per-bond nr_hb, per-bond hb_H_atom
     */
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
    const Matrix* getDC6DCNPtr() const { return m_d4_generator ? &m_d4_generator->getDC6DCN() : nullptr; }
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

    /// Compute shared packed-triangular distance arrays from m_geometry_bohr.
    /// Allocates/refreshes m_shared_sqrab and m_shared_srab (size N(N+1)/2).
    /// Thread-pool parallelized via threadPool().
    /// Called from GFNFF::Calculation() after prepareCNAndEEQ, before m_workspace->calculate().
    void computeSharedDistances() const;

    const Eigen::VectorXd& sharedSqrab() const { return m_shared_sqrab; }
    const Eigen::VectorXd& sharedSrab()  const { return m_shared_srab; }

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
     * @brief rev-gfnff pair-validity gate (Claude Generated, Sep 2026; gfnff_pair_validity.cpp;
     * FABLE_BOND_STATE_2.md sec 2.1-rev): every listed bonded pair of `topo` that VALID() rejects
     * - a closed-shell repulsion between two saturated, lone-pair-free, non-metal centres, with
     * no shared metal/deficient neighbour and no nearby antibonding excess charge. Reads only
     * `topo` (neighbour lists, is_metal, Phase-1 topology_charges) plus GFNFF::revValence() and
     * FFWorkspace::shareCapForAtom() - no geometry, no FFWorkspace instance. Canonical (i < j)
     * pairs, in no particular order. Called by generateGFNFFParameterSet() only when
     * -gfnff.rev_pair_validity is on.
     */
    std::vector<std::pair<int, int>> findInvalidPairValidityPairs(const TopologyInfo& topo) const;

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
    std::pair<std::vector<GFNFFRepulsion>, std::vector<GFNFFRepulsion>> generateRepulsionPairsNative() const;

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
    /// Claude Generated (Sep 2026): true if the repulsion list is cell-list built (distance-filtered),
    /// false if the O(N^2) build stores every pair (N < nb_cell_list_min_atoms) and so cannot go stale
    bool repulsionListIsDistanceFiltered() const;

    /// Detect hydrogen bonds as native GFNFFHydrogenBond structs
    std::vector<GFNFFHydrogenBond> detectHydrogenBondsNative(const Vector& charges) const;

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
    /// rev-gfnff stage 1: over-coordination energy of the last calculation (0 unless enabled)
    double OverCoordEnergy() const { return m_frag_blend_valid ? m_frag_blend_comp.over_coord : (m_workspace ? m_workspace->energyComponents().over_coord : 0.0); }
    /// rev-gfnff stage 2 (Sep 2026): bond-hardness energy of the split-charge model
    double SqeHardnessEnergy() const { return m_frag_blend_valid ? m_frag_blend_comp.sqe_hardness : (m_workspace ? m_workspace->energyComponents().sqe_hardness : 0.0); }

    // Claude Generated (April 2026): PBC accessors for GPU path
    bool hasPBC() const { return m_has_pbc; }
    Eigen::Matrix3d getUnitCellBohr() const {
        constexpr double ANG2BOHR = 1.0 / 0.529177210903;
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
    std::shared_ptr<const GFNFFTables> m_tables = GFNFFTables::defaults(); ///< Claude Generated (Sep 2026): runtime parameter tables (defaults or -gfnff.param_file overrides)
    const GFNFFTables& T() const { return *m_tables; } ///< the tables this instance evaluates with
    RevSettings m_rev_settings;                 ///< rev-gfnff stage 1 (Sep 2026): kernel switches and per-atom data
    double m_rev_bo_form = 0.05, m_rev_bo_break = 0.02; ///< rev-gfnff react scan thresholds on the term weight
    bool m_rev_form_order = true;   ///< rev_form_switch == "order": the narrow bond order decides a formation (Sep 12, 2026)
    double m_rev_bo2_form = 0.1;    ///< formation threshold on the narrow bond order (rev_form_switch == "order")
    double m_rev_tr_begin = 0.02, m_rev_tr_end = 0.8, m_rev_tr_revert = 0.75;
    double m_rev_tr_prebreak = 0.5;    ///< stage 1b: a bond starts breaking below this coordinate (below the formation end: a formed bond sticks)
    double m_rev_bo13_form = 0.1;      ///< stage 1b: a 1,3 pair closes a ring once its bond order (E_over switch) exceeds this
    bool m_rev_bo13_ordinary_join = false; ///< stage 3a(ii) diagnostic: drop the special 1,3 window, join on the ordinary criterion
    std::string m_rev_well_form = "mg3";    ///< stage 3a(iii): gauss | mg | erfmorse | mg2 | mg3 bond well (Sep 18, 2026; DEFAULT mg3 since Sep 22, 2026, was mg)
    std::string m_rev_share_form = "conserving"; ///< stage 3a(ii): delivered | conserving share formula (Sep 18, 2026; DEFAULT conserving since Sep 19, 2026)
    bool m_rev_share_donor_rule = true; ///< stage 3a(ii) conserving: a dative donor gets X_i >= 1 (Claude Generated, Sep 19, 2026)
    bool m_rev_budget_fix_h = true;    ///< stage 3a(ii): hydrogen keeps Val = 1 in the share budget (never hypervalent); Claude Generated Sep 15, 2026, DEFAULT ON since Sep 18, 2026
    bool m_rev_pair_validity = false;  ///< pair-validity gate (FABLE_BOND_STATE_2.md sec 2.1-rev); Claude Generated Sep 2026, DEFAULT OFF (gfnff_pair_validity.cpp)
    std::vector<Bond> m_rev_fading; ///< stage 1b: wells of broken bonds, kept in the bond list (weight w) until w < rev_bo_break
    std::map<std::pair<int, int>, long> m_rev_cooldown; ///< stage 1b: demoted pair -> first scan call at which it may start a transition again
    int m_rev_demote_cooldown = 0;                      ///< stage 1b: scans a demoted pair has to wait
    // stage 1b multi-transition state (Sep 12, 2026): corner bond lists / EEQ inputs by mask
    /// soft mu q0 rule (Claude Generated, Sep 24, 2026; _log/MU_CUSP_STATUS.md): the candidate
    /// integer placements of the mu rule and their weights w_p = exp(mu.q0_p / tau) / Z. The
    /// energy is the weighted average sum_p w_p E_p over the placements' full SQE energies.
    struct RevQ0Blend {
        std::vector<Vector> q0;              ///< q0 of each placement (whole units)
        std::vector<double> w;               ///< weights, sum 1, descending
        Vector q0u;                          ///< the uniform charges the probe mu was evaluated at
        EEQSolver::ProbeGeometryTerms terms; ///< alpha / cnf / cutoff of the probe matrix
        double tau = 0.0;                    ///< Eh
    };
    struct CornerEEQ {
        Vector topology_charges;
        std::vector<int> hybridization;
        std::optional<EEQSolver::TopologyInput> topo;
        std::optional<Vector> alpeeq;
        // rev-gfnff stage 2 (Sep 2026): the corner's split-charge model. `q0` is the ONE
        // discrete decision left (docs/REV_GFNFF_STAGE2.md) and is frozen when the corner is
        // created, at s = 0, so no energy jump can come from it; `sqe_pairs` is the corner's
        // bond list (fading wells included) plus the pairs of the transitions in flight, the
        // bond order b of each pair being taken at the current geometry every step.
        Vector q0;
        std::vector<std::pair<int, int>> sqe_pairs;
        /// rev-gfnff P3 (Sep 23, 2026): the flat excess-electron hardness per pair of sqe_pairs
        /// (same order), a corner constant like q0. Empty = all zero.
        std::vector<double> sqe_kappa_x;
        /// P3 alternative "frac" (Sep 23, 2026): fractional-charge correction factor per pair
        /// of sqe_pairs (same order), a corner constant. Empty = all zero.
        std::vector<double> sqe_frac_c;
        /// P3 alternative "harris" (Sep 24, 2026): topological excess-electron count x per pair
        /// of sqe_pairs (same order), a corner constant. Empty = all zero.
        std::vector<double> sqe_harris_x;
        /// soft mu q0 rule (Sep 24, 2026, _log/MU_CUSP_STATUS.md): set only when q0 came from the
        /// geometry-dependent mu rule at THIS step (slot corner, no transition in flight) and more
        /// than one placement carries weight. Frozen corners leave it empty.
        std::shared_ptr<RevQ0Blend> q0_blend;
    };
    bool m_rev_pending = false;                    ///< the next rebuild starts m_rev_pending_tr
    RevTransition m_rev_pending_tr;
    int m_rev_max_transitions = 4;                 ///< 2^k corners at most
    std::vector<std::vector<std::pair<int, int>>> m_rev_corner_bonds; ///< bond list of every corner (index = mask)
    std::vector<CornerEEQ> m_rev_corner_eeq;       ///< per-corner EEQ inputs for the per-step solve
    std::vector<std::pair<int, int>> m_rev_base_bonds; ///< topology before the first transition of a set
    CornerEEQ m_rev_base_eeq;
    void loadParameterOverrides();                 ///< param_file / param_json tables + rev settings (ctor and setParameters)
    /// EEQ inputs of the currently cached topology. `corner_generation` picks the q0 rule of
    /// the design: true = round the current fragment charge sums (a corner created by a
    /// topology event), false = the initialisation rule (integer fragment charges of qfrag
    /// spread uniformly). Falls back to the initialisation rule when no charges exist yet.
    CornerEEQ captureCornerEEQ(bool corner_generation = true);
    bool prepareTransitionCorners(std::vector<GFNFFParameterSet>& out); ///< generate the new corners except the all-ones one
    void appendFadingWells(GFNFFParameterSet& params) const;
    void installCornerPrepare();                   ///< per-corner EEQ solve callback into the workspace
    void dropTransitionCorners(int t, bool keep_new);
    void finishTransition(int t, bool keep_new, const char* why); ///< complete / revert / snap one transition (measured jump)
    void snapTransition(int t, const char* why, std::vector<std::pair<int, int>>* next);
    std::map<int, double> m_rev_p_over, m_rev_valence;  ///< per-element overrides from the rev section of the tables
    // ---- rev-gfnff stage 2 (Claude Generated, Sep 2026): split-charge model ------------------
    bool m_rev_sqe = false;                     ///< rev_charge_model == "sqe"
    double m_rev_sqe_bmin = 1e-3;               ///< pairs below this bond order are rigid
    std::map<int, double> m_rev_sqe_kappa;      ///< Z -> kappa_Z (Eh); absent = 0
    int m_rev_sqe_kappa_form = 0;               ///< B2: EEQSolver::SqeKappaForm as an int
    double m_rev_sqe_kappa_exponent = 3.0;      ///< B2: n of the power form
    bool m_rev_sqe_q0_mu = true;                ///< B2: q0 localised by chemical potential (else uniform)
    double m_rev_sqe_q0_mu_tau = 1.0 / 627.5094740631; ///< soft mu rule temperature (Eh); 0 = hard rule
    bool m_rev_sqe_base_q0_keep = true; ///< react: base corner inherits the slot's q0 at a transition start (Sep 29, 2026)
    /// soft mu q0 rule, this step's slot-corner blend (null = single placement, nothing to add)
    std::shared_ptr<RevQ0Blend> m_rev_q0_blend;
    std::vector<Vector> m_rev_q0_blend_q;         ///< SQE charges of every placement
    std::vector<Vector> m_rev_q0_blend_p;         ///< split charges of every placement (pair order)
    std::vector<double> m_rev_q0_blend_e;         ///< model energy of every placement (Eh)
    std::vector<EEQSolver::SqePair> m_rev_q0_blend_pairs; ///< the slot corner's pair list
    int m_rev_q0_blend_ref = 0;                   ///< placement the workspace evaluated
    double m_rev_q0_blend_de = 0.0;               ///< sum_p w_p E_p - E_ref, added to the energy
    /// dE/dx of the blend on top of the reference placement's workspace gradient (Eh/Bohr)
    Matrix revSqeQ0BlendGradient() const;
    /// explicit geometry gradient of the SQE model energy at fixed (q, p): Coulomb pairs,
    /// the CN part of chi and the pair hardness (Eh/Bohr)
    Matrix revSqeModelGradient(const Vector& q, const Vector& p, const EEQSolver::ProbeGeometryTerms& terms) const;
    /// sum_k v_k dmu_k/dx for the probe mu_k = x_k(CN) - sum_j A_kj(r) q0u_j (Eh/Bohr)
    Matrix revSqeMuDerivative(const Vector& v, const Vector& q0u, const EEQSolver::ProbeGeometryTerms& terms) const;
    // ---- rev-gfnff P2/P3 (Claude Generated, Sep 23, 2026; _log/P2P3_STATUS.md) ---------------
    bool m_rev_sqe_phase1 = false;              ///< P2: Phase-1 qa solved with the SQE model
    bool m_rev_sqe_virtual = false;             ///< Phase-2 virtual pairs across bond-graph components of one constraint group (X2_SCOPE_STATUS.md)
    bool m_rev_sqe_group_pairs = false;         ///< Phase-2 split-charge pairs only inside one constraint group (SQE_INVARIANT_STATUS.md)
    bool m_rev_excess = false;                  ///< P3: excess-electron perception on
    double m_rev_excess_kappa = 100.0;          ///< P3: flat hardness per excess electron (Eh)
    bool m_rev_excess_frac = false;             ///< P3 alternative: rev_excess_mode == "frac"
    double m_rev_excess_frac_c = 0.9;           ///< P3 alternative: c of the fractional-charge correction
    bool m_rev_excess_harris = false;           ///< P3 alternative: rev_excess_mode == "harris" (_log/P2P3_HARRIS_STATUS.md)
    bool m_rev_excess_react_consistent = true;  ///< P3 repair: q0 + kappa_x consistent over react corners
    /// P3 repair (a): re-localise q0 on every perceived excess-electron pair of a corner
    void revLocaliseExcessQ0(CornerEEQ& ce, const TopologyInfo& topo) const;
    /// P3: per-bond excess electrons x_ij of a topology (see the PARAM rev_excess_electron).
    std::map<std::pair<int, int>, double> revExcessElectrons(const TopologyInfo& topo) const;
    bool m_rev_pi_excess = false;               ///< P3 pi* prototype (_log/PI_STAR_STATUS.md)
    double m_rev_bond_extend = 1.0;             ///< rev_excess_bond_extend: 2c-3e candidate bond threshold factor (1 = off)
    /// P3 pi* prototype: per-bond pi* excess electrons y_ij of a diatomic radical anion
    std::map<std::pair<int, int>, double> revPiExcessElectrons(const TopologyInfo& topo) const;
    /// sigma x_ij + pi* y_ij of one pair (0 if not perceived)
    double revExcessTotal(const TopologyInfo& topo, int i, int j) const;
    /// P2: re-solve topo.topology_charges with the split-charge model on the Phase-1 matrix and
    /// refresh alpeeq/dgam from them; also fills topo.rev_excess (P3). Runs at the end of every
    /// calculateTopologyInfoOnce() pass. No-op unless m_rev_sqe && (m_rev_sqe_phase1 || m_rev_excess).
    void revApplyPhase1Sqe(TopologyInfo& topo) const;
    /// P3: x_ij kappa_x for one pair of a topology (0 if not perceived)
    double revExcessKappa(const TopologyInfo& topo, int i, int j) const;
    /// P3 alternative "frac": min(x_ij, 1) c for one pair of a topology (0 if not perceived / flat mode)
    double revExcessFracC(const TopologyInfo& topo, int i, int j) const;
    /// P3 alternative "harris": x_ij for one pair of a topology (0 if not perceived / not harris mode)
    double revExcessHarrisX(const TopologyInfo& topo, int i, int j) const;
    double revSqeKappa(int Z) const { auto it = m_rev_sqe_kappa.find(Z); return it == m_rev_sqe_kappa.end() ? 0.0 : it->second; }
    /// the pair set of a corner: its bond graph + the fading wells + the pairs in transition
    std::vector<std::pair<int, int>> revSqePairs(const std::vector<std::vector<int>>& neighbor_lists) const;
    /// q0 of the initialisation (static / cold-start) rule for the integer fragment charges of
    /// qfrag. `uniform` spreads them flat; `mu` (B2, default) localises them on the extreme-
    /// chemical-potential atoms of the fragment, which is what gives kappa a lever at all.
    /// The mu variant needs the EEQ inputs of the corner the q0 belongs to.
    Vector revSqeQ0Fragments(const EEQSolver::TopologyInput& ti,
                             const Vector& topology_charges,
                             const std::vector<int>& hybridization,
                             const std::optional<Vector>& alpeeq,
                             std::shared_ptr<RevQ0Blend>* blend_out = nullptr) const;
    /// soft mu rule (Sep 24, 2026): the integer placements of the mu rule within 34 tau of the
    /// best one and their Boltzmann weights over sum_i e_i; returns the best placement's q0
    Vector revSqeQ0MuBlend(const Vector& q0u, const Vector& mu,
                           const EEQSolver::ProbeGeometryTerms& terms, const std::vector<double>& qf,
                           const std::vector<int>& count, const std::vector<int>& frag,
                           std::shared_ptr<RevQ0Blend>* blend_out) const;
    /// q0 of the corner-generation rule: rounded fragment sums of `q_now`, residual to the
    /// fragment with the lowest (highest) EEQ chemical potential when an electron (a hole) is left over
    Vector revSqeQ0Rounded(const CornerEEQ& ce, const Vector& q_now) const;
    /// the split-charge data of the SLOT corner (the expected topology): its pair set comes
    /// from the current cached topology, its q0 from the stored all-ones corner if a transition
    /// set exists (frozen at creation) and otherwise from the initialisation rule
    CornerEEQ revSlotCorner(const TopologyInfo& topo) const;
    /// solve the split-charge system of one corner and push its p values into the workspace slot
    Vector revSolveSplitCharges(const CornerEEQ& ce, const Vector& topology_charges,
                                const std::vector<int>& hybridization,
                                const std::optional<EEQSolver::TopologyInput>& topo,
                                const std::optional<Vector>& alpeeq,
                                CxxThreadPool* pool, int threads);
    void setupRevSettings();                    ///< fill m_rev_settings from PARAMs + tables (ctor) 
    void fillRevPerAtom();                      ///< per-atom rcov/fat/p/valence (after the atoms are known)
    double revValence(int Z) const;             ///< nominal sigma valence of an element
    double revOverP(int Z) const;               ///< over-coordination prefactor p_Z of an element (preset / rev section / rev_over_p)
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

    // Reused-topology invalidation (Claude Generated Sep 2026). m_ff_bond_graph is the
    // canonical bond graph of the topology the CURRENT interaction lists were built from
    // (recorded in generateGFNFFParameterSet); m_reuse_seen_bonds caches the perception of
    // the current geometry that is compared against it each energy call.
    std::vector<std::pair<int,int>> m_ff_bond_graph;
    std::vector<std::pair<int,int>> m_reuse_seen_bonds;
    Eigen::MatrixXd m_reuse_seen_geometry;
    bool m_reuse_topology_check = false; ///< PARAM reuse_topology_check (parsed in the ctor; false = trust the first frame)

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
    /// rev-gfnff stage 3b (Claude Generated, Sep 2026): the CONTINUOUS bond order of a pair, the
    /// key of the bond-order-resolved well table (Bond::rev_order). See the definition for why
    /// this is the same expression as refreshReactBondOrders() without its rounding/threshold.
    double continuousBondOrder(int i, int j, const TopologyInfo& topo) const;

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
    std::string m_topology_mode = "auto";

    // Claude Generated (March 2026): Topology persistence in param.json
    bool m_cache_topology = true;   ///< Cache Phase-1 EEQ topology in param.json (opt-out)
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
    // See docs/FRAG_CHARGE_MODEL.md, FRAG_CHARGE_STATUS.md section 1. A charge VARIANT is a complete GFNFF instance whose
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
    bool m_frag_atomic_ea = false;       ///< frag_charge_atomic_ea: single-atom anion carriers ranked by atomic EA (I2_CLF_STATUS.md)
    double m_frag_ea_sigma = 0.02;       ///< eV, width of the EA softmax
    int m_frag_max_edges = 4;
    int m_frag_max_placements = 6;
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

    // WP-FF-DistMatrix-Sharing (May 2026): shared packed-triangular distance arrays.
    // Filled by computeSharedDistances() once per energy call; consumed by
    // ForceFieldThread term loops and EEQSolver phase-2.  Layout: lower-triangular,
    // index via GFNFF::triIdx(i, j) = i*(i+1)/2 + j (i > j).
    mutable Eigen::VectorXd m_shared_sqrab;
    mutable Eigen::VectorXd m_shared_srab;
    mutable int             m_shared_dist_N = 0;

    // Conversion factors
    static constexpr double HARTREE_TO_KCAL = 627.5094740631;
    static constexpr double BOHR_TO_ANGSTROM = 0.5291772105638411;
    static constexpr double KCAL_TO_HARTREE = 1.0 / 627.5094740631;
    static constexpr double ANGSTROM_TO_BOHR = 1.0 / 0.5291772105638411;

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
