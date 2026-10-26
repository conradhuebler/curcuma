/*
 * < EEQ (Electronegativity Equalization) Charge Solver >
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
 * This file contains a standalone EEQ (Electronegativity Equalization) solver
 * extracted from GFN-FF implementation for reusability across force field methods.
 *
 * References:
 * - S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673
 *   (GFN-FF with EEQ charge model)
 * - R. T. Sanderson, Chemical Bonds and Bond Energy, Academic Press, 1976
 *   (Original electronegativity equalization concept)
 */

#pragma once

#include "src/core/global.h"
#include "src/core/config_manager.h"
#include "src/core/parameter_macros.h"
#include "src/core/charge_extrapolation.h"

#include <Eigen/Dense>
#include <deque>
#include <functional>
#include <string>
#include <vector>
#include <memory>
#include <optional>

class CxxThreadPool;  // Forward declaration for pool-based parallelisation

/**
 * @brief Distance mode for EEQ matrix construction
 *
 * Controls which type of interatomic distances are used in the EEQ Coulomb matrix:
 * - Topological: Floyd-Warshall bond path distances (for Phase 1 topology charges)
 * - Geometric: Euclidean xyz distances (for Phase 2 energy charges)
 *
 * CRITICAL FIX (January 17, 2026): The Fortran reference uses DIFFERENT distance modes
 * for Phase 1 vs Phase 2:
 * - Phase 1 (goedeckera): Topological distances (`pair(ij)` from Floyd-Warshall)
 * - Phase 2 (goed_gfnff): Geometric distances (`r(ij)` from xyz coordinates)
 *
 * Claude Generated - January 17, 2026
 */
enum class EEQDistanceMode {
    Topological,  ///< Phase 1: Floyd-Warshall bond path distances
    Geometric     ///< Phase 2: Euclidean xyz distances
};

/**
 * @brief Linear solve method for EEQ system
 *
 * Controls the algorithm used to solve the augmented EEQ linear system:
 * - LU: PartialPivLU on full augmented matrix (baseline, always works)
 * - SchurCholesky: Exploit SPD structure of NxN Coulomb sub-matrix via Cholesky,
 *   then Schur complement for constraint rows (~2x faster than LU)
 * - PCG: Preconditioned Conjugate Gradient with Schur complement constraint
 *   handling. O(N²·k) where k << N for warm-started MD/optimization (~10-30x for N>500)
 * - Auto: First-call benchmark chooses fastest method (SchurCholesky or PCG),
 *   then caches decision for subsequent calls. Default for optimal performance.
 *
 * Claude Generated - March 2026 (Performance optimization)
 */
enum class EEQSolveMethod {
    LU,             ///< PartialPivLU on full augmented system (baseline)
    SchurCholesky,  ///< Cholesky on NxN SPD block + Schur complement for constraints
    LDLT,           ///< Bunch-Kaufman LDL^T on NxN block (indefinite-capable) + Schur; LLT-fail fallback
    Batched,        ///< Per-fragment Cholesky (GPU only); CPU falls back to SchurCholesky
    PCG,            ///< Preconditioned Conjugate Gradient with warm start
    ProjectedPCG,   ///< Constraint-projected PCG on A_nn: ONE iterative solve for any number of fragments
    Auto            ///< Threshold rule: ProjectedPCG for N >= eeq_ppcg_min_atoms and nfrag >= eeq_ppcg_min_nfrag, else SchurCholesky. Resolved by EEQSolver::resolveMethod (CPU and GPU share it).
};

/**
 * @brief Local-inverse preconditioner for the projected PCG of the EEQ system (opt-in,
 *        `-eeq_solver.eeq_ppcg_precond localinv`). Claude Generated (Oct 2026).
 *
 * **Problem.** The fragment-constrained EEQ system is the KKT system
 *     A q + C^T lambda = b,   C q = d,
 * A = A_nn (N x N, SPD: hardness diagonal + erf(r_ij/sqrt(a_i+a_j))/r_ij), C the 0/1 fragment
 * membership (nfrag x N). The projected CG (EEQSolver::solveWithProjectedPCG) iterates on the
 * tangent space null(C): residuals are projected with P = I - C^T (C C^T)^-1 C.
 *
 * **Constraint preconditioner** (Gould, Hribar, Nocedal, SIAM J. Sci. Comput. 23 (2001) 1376;
 * Keller, Gould, Wathen, SIAM J. Matrix Anal. Appl. 21 (2000) 1300). For any SPD G ~ A^-1 the
 * preconditioned residual is the solution z of
 *     [G^-1 C^T; C 0] [z; mu] = [r; 0]   =>   z = P_G r,
 *     P_G = G - G C^T S^-1 C G,   S = C G C^T  (nfrag x nfrag, SPD).
 * Properties used by the CG: (i) C P_G = 0, so z and every search direction stay in null(C) and
 * q keeps satisfying C q = d; (ii) P_G is symmetric, P_G C^T = 0, and for r in null(C), r != 0,
 *     r^T P_G r = |G^1/2 r|^2 - |G^1/2 C^T S^-1 C G r|^2 > 0
 * (Cauchy-Schwarz in the G inner product; equality only for r in range(C^T), which is orthogonal
 * to null(C)) - P_G is positive definite on the tangent space exactly when G is SPD. (iii) For
 * G = A^-1, P_G = (P A P)^+ on null(C): the ideal preconditioner. The Jacobi path is the special
 * case G = diag(A)^-1, where S is diagonal (fragments are disjoint).
 *
 * **G = local inverse** (a sparse approximate inverse built from local dense inverses; cf.
 * Benzi, J. Comput. Phys. 182 (2002) 418, section 5). A is long-ranged, A^-1 is short-ranged
 * (screening by the hardness diagonal). For each atom i: E_i = atoms within R_ext of i, invert
 * the dense block A[E_i,E_i] (here: Cholesky + one solve for the column of i), keep the entries
 * of row i with r_ij <= R. G_ij = (Gl_ij + Gl_ji)/2. Symmetrising does not guarantee
 * definiteness, so G is checked with a sparse Cholesky; S is factorised the same way. Any
 * failure (or rz <= 0 / pAp <= 0 during the CG) means Jacobi is used instead - never a wrong M.
 *
 * Cost: build O(N n_ext^3 / 3) with n_ext ~ atoms within R_ext (independent per atom, parallel,
 * the result does not depend on the thread count); per CG iteration two sparse products and
 * one nfrag-size sparse triangular solve pair, against one dense N^2 matvec.
 */
struct LocalInversePreconditioner {
    bool ok = false;           ///< false: build failed (G or S not SPD) - use Jacobi
    std::string why;           ///< reason of a failed build
    int natoms = 0;
    int nfrag = 0;
    std::vector<int> frag;     ///< fragment index of every atom (from C)
    Matrix geometry;           ///< geometry (Bohr) of the build, for the rebuild criterion
    Eigen::SparseMatrix<double, Eigen::RowMajor> G;  ///< symmetric local inverse, N x N
    Eigen::SparseMatrix<double, Eigen::RowMajor> H;  ///< G C^T, N x nfrag
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> S_llt;  ///< factor of S = C G C^T
    Matrix S_inv;              ///< dense S^-1 (nfrag x nfrag) when the GPU applies it (else empty)
    double nnz_per_row = 0.0;  ///< mean stored entries per row of G
    double mean_block = 0.0;   ///< mean dense block size n_ext
    double build_ms = 0.0;
    double blocks_ms = 0.0;    ///< part of build_ms: neighbour lists + dense block inverses
    double check_ms = 0.0;     ///< part of build_ms: sparse Cholesky of G (definiteness check)
    double sinv_ms = 0.0;      ///< part of build_ms: dense S^-1 (GPU path only)

    /// z = P_G r = G r - H S^-1 C G r. Deterministic (fixed summation order).
    void apply(const Vector& r, Vector& z, Vector& t) const;
};

/**
 * @brief EEQ (Electronegativity Equalization) charge solver
 *
 * Standalone utility for calculating atomic partial charges using the EEQ method.
 * Implements the two-phase algorithm from GFN-FF:
 * - Phase 1: Topology charges (quick initial charges from augmented linear system)
 * - Phase 2: Final charges (iterative refinement with environmental corrections)
 *
 * Usage:
 *   ConfigManager config("eeq_solver", json{});
 *   EEQSolver solver(config);
 *   Vector charges = solver.calculateCharges(atoms, geometry_bohr, total_charge);
 *
 * Claude Generated - December 2025
 */
class EEQSolver {
public:
    /// Read-only view of a (row-major) Matrix block, e.g. A.topLeftCorner(n, n), without copying.
    using MatrixCRef = Eigen::Ref<const Matrix, 0, Eigen::OuterStride<>>;

    /**
     * @brief Topology information for EEQ Phase 1
     *
     * Used to compute topological distances via Floyd-Warshall algorithm
     * instead of geometric (line-of-sight) distances.
     *
     * Topological distances represent shortest paths through the bond graph,
     * where each bond has length = sum of covalent radii. This matches the
     * XTB implementation in gfnff_ini.f90:431-461.
     *
     * Claude Generated December 2025
     */
    struct TopologyInput {
        std::vector<std::vector<int>> neighbor_lists;  // neighbor_lists[i] = atoms bonded to i
        std::vector<double> covalent_radii;            // Covalent radii in ANGSTROM (NOT Bohr!)
        int nfrag = 1;                                 // Number of molecular fragments
        std::vector<int> fraglist;                     // fraglist[i] = fragment ID for atom i (1-indexed)
        std::vector<double> qfrag;                     // qfrag[f] = target charge for fragment f
        // Fortran itag (gfnff_ini2.f90): +1 carbene C / NO2 N, -1 eta-coordinated, 0 otherwise.
        // Supplied so the dxi corrections can test the REAL tag instead of re-deriving it from
        // geometry: the stored tag has already been through the aryne rule and the qa < -0.4
        // override, and re-deriving it silently skips both. Empty = fall back to the heuristic.
        std::vector<int> itag;
        // Fortran piadr as the EEQ section sees it: the PRE-Hückel pi-candidate list
        // (gfnff_ini.f90:312-336), i.e. GFNFF::detectPiSystems()'s output, nonzero for a
        // candidate. Supplied so the dgam / amide rules use the SAME array as the rest of
        // the force field instead of EEQSolver::detectPiSystem()'s own older inference,
        // which has neither the NR3-X / SO3 vetoes nor B/Cl/S in its element lists.
        // Empty = fall back to that inference.
        std::vector<char> is_pi;
    };

    /**
     * @brief Construct EEQ solver from configuration
     * @param config ConfigManager with eeq_solver parameters
     */
    explicit EEQSolver(const ConfigManager& config);

    /**
     * @brief Default destructor
     */
    ~EEQSolver() = default;

    // ===== Convergence Statistics =====

    /**
     * @brief Print summary of PCG convergence statistics since last reset.
     *
     * Prints one line: "EEQ PCG: X/Y converged (worst |r|=Z.ZeZZ)" if any
     * calls did not converge, or nothing if all converged.
     * Resets statistics after printing.
     */
    void printConvergenceSummary();

    /**
     * @brief Reset convergence statistics counters.
     */
    void resetConvergenceStats();

    // ===== Main API =====

    /**
     * @brief Calculate atomic partial charges using two-phase EEQ algorithm
     *
     * @param atoms Atomic numbers (Z) for each atom
     * @param geometry_bohr Cartesian coordinates in Bohr (natoms x 3)
     * @param total_charge Total molecular charge (default 0)
     * @param cn_hint Optional pre-calculated coordination numbers (nullptr to calculate internally)
     * @param hyb_hint Optional pre-calculated hybridization states (nullptr to use defaults)
     * @param topology Optional topology information for Floyd-Warshall topological distances
     *                 If provided, uses topological distances (more accurate, matches XTB)
     *                 If omitted, falls back to geometric distances (deprecated)
     * @return Vector of atomic partial charges in elementary charge units (e)
     *
     * Automatically uses two-phase algorithm for better accuracy:
     * - Phase 1: Solve augmented EEQ system for topology charges
     * - Phase 2: Iterative refinement with dxi, dgam, dalpha corrections
     *
     * @note CRITICAL FIX (Dec 28, 2025): Added topology parameter to enable Floyd-Warshall
     *       topological distances in Phase 1, matching XTB reference implementation.
     *       Without topology, geometric distances produce charges 4-5× too large!
     */
    Vector calculateCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge = 0,
        const Vector* cn_hint = nullptr,
        const std::vector<int>* hyb_hint = nullptr,
        const std::optional<TopologyInput>& topology = std::nullopt
    );

    // ===== Advanced API (two-phase workflow) =====

    /**
     * @brief Phase 1: Calculate initial topology charges
     *
     * Solves augmented EEQ linear system:
     * - A_{ii} = gam_i + sqrt(2/π)/sqrt(alpha_i)
     * - A_{ij} = erf(γ_ij * r_ij) / r_ij
     * - Constraint row/column for charge conservation
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param total_charge Total molecular charge
     * @param cn Coordination numbers (for chi corrections)
     * @param topology Optional topology information (neighbor lists, covalent radii)
     *                 If provided, uses Floyd-Warshall topological distances
     *                 If omitted, falls back to geometric distances (deprecated)
     * @return Vector of topology charges (qa)
     *
     * @note Topological distances computed via Floyd-Warshall algorithm match XTB reference
     *       and produce charges ~4-5× smaller than geometric distances.
     *       Reference: XTB gfnff_ini.f90:431-461
     */
    Vector calculateTopologyCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt,
        bool use_corrections = false,  // CRITICAL FIX (Jan 4, 2026): default false to match gfnff_final.cpp
        CxxThreadPool* pool = nullptr,  // Claude Generated (WP2, May 2026): pool forwarded to dispatchSolve for Stage-4 batched parallelisation
        int num_threads = 1
    );

    /**
     * @brief Phase-1 topology charges for several fragment-charge assignments (A1, Jul 2026)
     *
     * The augmented Phase-1 matrix A depends on dxi, chi/gam/alpha, the topological
     * distances and `fraglist` — but NOT on `qfrag`, which enters only the constraint
     * RHS `x(natoms+f)`. The GFN-FF ipis block (gfnff_method.cpp) re-solves Phase 1
     * once per pi-system with a different qfrag each time, so it was rebuilding a
     * bit-identical matrix (Dijkstra + N x N erf fill) for every solve.
     *
     * This builds the system once and then solves it for each supplied qfrag,
     * patching only the constraint rows of the RHS. Each solve calls the same
     * dispatchSolve() with the same A and the same x as the per-call path, so the
     * returned charges are bit-identical to calling calculateTopologyCharges() in a
     * loop.
     *
     * @param qfrag_variants One qfrag vector per desired solve (size nfrag each).
     * @return One charge vector per variant, in order.
     */
    std::vector<Vector> calculateTopologyChargesMultiRHS(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& cn,
        const std::optional<TopologyInput>& topology,
        const std::vector<std::vector<double>>& qfrag_variants,
        bool use_corrections = false,
        CxxThreadPool* pool = nullptr,
        int num_threads = 1
    );

    /**
     * @brief Phase 2: Calculate final EEQ charges with environmental corrections
     *
     * Single linear solve (matches XTB gfnff_ini.f90:699-706) with corrected parameters:
     * - chi' = chi + dxi(environment)
     * - gam' = gam + dgam(charge-dependent, uses topology_charges)
     * - alpha' = alpha_base² (BASE VALUE, NOT charge-corrected!)
     *
     * CRITICAL (Jan 7, 2026): Alpha uses BASE values in Phase 2, unlike Phase 1c which
     * applies alpha' = (alpha_base + ff*qa)². This is the key difference between topology
     * and energy charge calculations. The dxi/dgam corrections are applied, but alpha
     * remains constant to maintain numerical stability.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param total_charge Total molecular charge
     * @param topology_charges Phase-1 charges (qa) used for dgam calculation (NOT for alpha!)
     * @param cn Coordination numbers
     * @param hybridization Hybridization states (1=sp, 2=sp2, 3=sp3)
     * @param topology Optional topology information for integer neighbor counts in CNF calculation
     * @param alpeeq Optional charge-dependent alpha values (squared) from Phase 1
     *               If provided, these pre-computed alpeeq values are used instead of base alpha
     *               Formula: alpeeq(i) = (alpha_base + ff*qa(i))²
     *               Reference: Fortran gfnff_ini.f90:718-725
     * @return Vector of final EEQ charges
     */
    Vector calculateFinalCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& topology_charges,
        const Vector& cn,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology = std::nullopt,
        bool use_corrections = false,  // CRITICAL FIX (Jan 4, 2026): default false to match gfnff_final.cpp
        const std::optional<Vector>& alpeeq = std::nullopt,  // Claude Generated (January 2026): Charge-dependent alpha
        CxxThreadPool* pool = nullptr,  // Claude Generated (Mar 2026): Pool-based parallelisation
        int num_threads = 1
    );

    /**
     * @brief Calculate dgam corrections with full pi-system and amide detection
     *
     * Claude Generated (March 2026): Public interface for computing dgam with
     * complete N-specific corrections (ff=-0.14 for pi-N, ff=-0.16 for amide N).
     * Ensures Coulomb energy parameters are consistent with EEQ solver parameters.
     * Reference: Fortran gfnff_ini.f90:697-724 computes dgam ONCE for both.
     *
     * @param atoms Atomic numbers
     * @param topology_charges Phase-1 charges (qa)
     * @param hybridization Hybridization states
     * @param cn Coordination numbers
     * @param topology Topology information (neighbor lists, etc.)
     * @return Vector of dgam corrections
     */
    Vector calculateDgamFull(
        const std::vector<int>& atoms,
        const Vector& topology_charges,
        const std::vector<int>& hybridization,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    );

    /**
     * @brief Calculate dxi corrections with full environment detection
     *
     * Claude Generated (March 2026): Public interface for computing dxi with
     * complete Fortran-matching corrections (carbene C, free CO, nitro O,
     * polyvalent halogens, pi-system neighbor EN averaging, etc.).
     * Ensures Coulomb energy dxi is consistent with EEQ solver dxi.
     * Reference: Fortran gfnff_ini.f90:358-403
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @param topology Topology information (neighbor lists, etc.)
     * @return Vector of dxi corrections
     */
    Vector calculateDxiFull(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    );

    /**
     * @brief Detect amide hydrogens (public wrapper)
     *
     * Claude Generated (March 2026): Public access to amide hydrogen detection
     * for use in Coulomb chi_base correction.
     * Reference: Fortran gfnff_ini.f90:717 - amideH correction to chieeq
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param cn Coordination numbers
     * @param topology Topology information
     * @return Vector of amide hydrogen flags
     */
    std::vector<bool> detectAmideHydrogensFull(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    ) const;

    /**
     * @brief Parse solve method string to enum (public so GPU dispatch can reuse).
     * Claude Generated - March 2026; moved to public for WP7-B (May 2026).
     *
     * Canonical names: auto, cholesky, ldlt, lu, pcg, ppcg, batched. Aliases: schur_cholesky
     * (= cholesky), projected_pcg (= ppcg). Oct 2026: an unknown name throws
     * std::invalid_argument listing the valid names (it used to fall back to auto silently).
     */
    static EEQSolveMethod parseSolveMethod(const std::string& method_str);

    /// Canonical CLI name of a method ("cholesky", "ppcg", ...).
    static const char* methodName(EEQSolveMethod m);

    /// Print an "EEQ solver: ..." line at verbosity 1, once per distinct line process-wide
    /// (several calculator instances, e.g. batch workers, must not repeat it). Thread-safe.
    static void reportSolverLine(const std::string& line);

    /// Comma separated list of the canonical names (for messages).
    static std::string validMethodNames();

    /// True when the GPU back-end (`cuda`, `rocm`, `vulkan`) has a device implementation of the method
    /// (auto resolves to a supported one). Oct 2026: lu (dense augmented LU, cuSOLVER getrf) and ldlt
    /// (cusolverDnDsytrf + Xsytrs on A_nn, Schur) run on CUDA only - requesting an unsupported one with
    /// -gpu is an error, not a silent substitute. ROCm: only auto and cholesky (no device ppcg/pcg/batched).
    /// Vulkan: GFN-FF runs on the CPU, so every name is accepted.
    static bool gpuSupports(EEQSolveMethod m, const std::string& backend = "cuda");

    /// Comma separated names a back-end can run (for the error message). Claude Generated (Oct 2026).
    static std::string gpuMethodNames(const std::string& backend);

    /**
     * @brief Energy error of the batched solve, estimated to second order. Claude Generated (Oct 2026).
     *
     * The batched solver minimises the EEQ functional E(q) = 1/2 q^T A q - x^T q with the
     * inter-fragment blocks of A set to zero. Its charges q_b leave the residual r = A q_b - x of
     * the full system; with the per-fragment charge constraints the exact energy excess is
     * E(q_b) - E(q*) = 1/2 r^T (P A P)^+ r. Replacing A by its fragment-diagonal part gives
     *   dE ~ 1/2 sum_f r_f^T d_f,   [A_ff 1; 1^T 0] [d_f; mu] = [r_f; 0],
     * i.e. one block-Jacobi step with the same blocks the batched solve factorised. Exact when
     * the coupling is weak; it only uses the residual (one O(N^2) matrix-vector product) and
     * the fragment blocks.
     *
     * @param residual        r = A q_b - x over all atoms (A_nn and the atom right-hand side)
     * @param fragment_atoms  atom indices of each fragment
     * @param A_elem          element (i, j) of A_nn for atoms of the same fragment
     * @return estimated E(q_b) - E(q*) in Eh (>= 0)
     */
    static double estimateBatchedError(const Vector& residual,
                                       const std::vector<std::vector<int>>& fragment_atoms,
                                       const std::function<double(int, int)>& A_elem);

    /// Last batched-error estimate of this solver (Eh), -1 when none was made. Claude Generated (Oct 2026).
    double lastBatchedErrorEstimate() const { return m_last_batched_error; }

    /**
     * @brief Resolve the requested method to the one that runs. Claude Generated (Oct 2026).
     *
     * ONE rule for CPU and GPU. Every explicit method is returned unchanged (the user's choice
     * is never overridden by a size threshold). Auto becomes ProjectedPCG when
     * natoms >= min_atoms and min_nfrag > 0 and nfrag >= min_nfrag (eeq_ppcg_min_atoms /
     * eeq_ppcg_min_nfrag, default 500 / 1), otherwise SchurCholesky. Auto never picks an
     * inexact method (batched).
     */
    static EEQSolveMethod resolveMethod(EEQSolveMethod requested, int natoms, int nfrag,
                                        int min_atoms, int min_nfrag);

    /// WP-S3 (May 2026): post-init override for `eeq_distance_cutoff` (Bohr).
    /// Set by GFNFF after Phase-1 when the auto-detection heuristic triggers
    /// (nfrag==1 and max|topology_charges| < 0.5 e). Pass a negative value
    /// to clear the override and fall back to the ConfigManager default.
    void setEEQDistanceCutoff(double cutoff) { m_eeq_distance_cutoff_override = cutoff; }

    /// WP-EEQ-Cache (May 2026): invalidate cached Cholesky factor.
    /// Call after setMolecule() or topology rebuild to force refactorization on next step.
    ///
    /// Also clears the Phase-2 pending buffers (Jul 2026). solveWithSchurCholesky
    /// treats "pending buffers non-empty" as "I am Phase 2" (:2113-2119, :2166-2169);
    /// resetting the factor while leaving them populated let the next Phase 1 persist
    /// its topological factor into the freshly-emptied cache, which Phase 2 then
    /// consumed. Dropping the factor and its key together keeps the two consistent.
    void invalidateCholeskyCache()
    {
        m_chol_cache.reset();
        m_pending_geometry.resize(0, 0);
        m_pending_cn.resize(0);
    }

    /// WP-EEQ-Matrix-Cache (May 2026): invalidate cached A_nn off-diagonal.
    /// Call after topology rebuild to force matrix re-build on next Phase 2.
    void invalidateMatrixCache() { m_A_cache_valid = false; }

    /// WP-FF-DistMatrix-Sharing (May 2026): consume externally-computed packed
    /// lower-triangular distance array (srab[i*(i+1)/2 + j] for i > j, Bohr).
    /// When set, calculateFinalCharges skips its own O(N^2) distance loop and
    /// fills m_phase2_distances from this view. Lifetime managed by caller (GFNFF).
    /// too_close (optional, Oct 2026): flag the caller computed while writing srab (any off-diagonal
    /// r_ij < 1e-10); when given, calculateFinalCharges skips its own atoms-too-close scan over srab.
    void setExternalDistances(const Eigen::VectorXd* srab, const bool* too_close = nullptr)
    {
        m_external_srab = srab;
        m_external_srab_too_close = too_close;
    }

    /// WP5 (June 2026): self-consistent implicit-solvation reaction field. When set,
    /// calculateFinalCharges adds this symmetric nat×nat Born matrix B into the
    /// top-left block of the augmented EEQ matrix before the linear solve, so the
    /// charges feel the solvent (matches the gfnff reference A_eeq += gbsa%bornMat,
    /// gfnff_engrad.F90:1346-1350). nullptr = gas phase. Lifetime managed by caller
    /// (GFNFF holds the ALPBSolvation that owns the matrix). Claude Generated.
    void setReactionField(const Eigen::MatrixXd* born_matrix) { m_reaction_field = born_matrix; }

    /// WP-S3 (May 2026): effective Coulomb-matrix cutoff used by the
    /// sparsification (override wins, otherwise the ConfigManager value).
    double getEEQDistanceCutoffEffective() const {
        return (m_eeq_distance_cutoff_override >= 0.0)
                   ? m_eeq_distance_cutoff_override
                   : m_config.get<double>("eeq_distance_cutoff", 0.0);
    }

private:
    /**
     * @brief Element-specific EEQ parameters
     *
     * From Spicher & Grimme Angew. Chem. Int. Ed. 2020, 59, 15665 (angewChem2020 set)
     */
    struct EEQParameters {
        double chi;    ///< Electronegativity (Hartree)
        double gam;    ///< Chemical hardness (Hartree)
        double alp;    ///< Damping parameter (Bohr^-1, SQUARED!)
        double cnf;    ///< CN correction factor (dimensionless)
    };

    /**
     * @brief Get EEQ parameters for element Z
     *
     * @param Z Atomic number (1-86)
     * @param cn Coordination number (for chi correction)
     * @return EEQParameters structure
     *
     * CRITICAL: alpha value is SQUARED (alpha_eeq[Z-1]²) before storage!
     */
    EEQParameters getParameters(int Z, double cn = 0.0) const;

    /**
     * @brief Generate uniform fallback charges when EEQ solver fails
     * Claude Generated (March 2026): Graceful degradation instead of hard crash
     *
     * @param natoms Number of atoms
     * @param total_charge Total molecular charge
     * @param context Description of failure context (for warning message)
     * @return Vector with q_i = total_charge / natoms
     */
    Vector generateFallbackCharges(int natoms, int total_charge, const std::string& context) const;

    // ===== Correction Terms =====

    /**
     * @brief Calculate electronegativity corrections (dxi)
     *
     * Environment-dependent chi corrections based on coordination number
     * and neighbor electronegativity.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @return Vector of dxi corrections (Hartree)
     */
    Vector calculateDxi(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    );

    /**
     * @brief Calculate hardness corrections (dgam)
     *
     * Charge-dependent gamma corrections with element-specific factors.
     * Formula: dgam(i) = qa(i) * ff
     * where ff depends on element and hybridization.
     *
     * @param atoms Atomic numbers
     * @param charges Current charges (from Phase 1)
     * @param hybridization Hybridization states
     * @return Vector of dgam corrections (Hartree)
     */
    Vector calculateDgam(
        const std::vector<int>& atoms,
        const Vector& charges,
        const std::vector<int>& hybridization,
        const std::vector<bool>& is_pi_atom = {},
        const std::vector<bool>& is_amide = {}
    );

    // NOTE: calculateDalpha() removed - alpha now calculated inline with charge-dependent formula

    // ===== NEW Single-Solve Helper Functions (Jan 4, 2026) =====

    /**
     * @brief Build EEQ Coulomb matrix with ALL charge-dependent corrections
     *
     * Constructs the augmented EEQ matrix with:
     * - Diagonal: gam_corrected + sqrt(2/π)/sqrt(alpha_corrected)
     * - Off-diagonal: Coulomb terms with corrected alpha
     * - Constraint row/column for charge conservation
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @param current_charges Current charge estimate (used for charge-dependent dgam/alpha)
     * @param dxi Electronegativity corrections
     * @param dgam Hardness corrections (gam - qa*ff)
     * @param hybridization Hybridization states
     * @param topology Optional topology for topological distances
     * @param distance_mode Distance calculation mode (Topological for Phase 1, Geometric for Phase 2)
     *                      CRITICAL (Jan 17, 2026): Fortran uses topological for Phase 1 (goedeckera),
     *                      geometric for Phase 2 (goed_gfnff). Default=Topological for backward compat.
     * @return Augmented EEQ matrix (natoms+1)×(natoms+1)
     */
    Matrix buildCorrectedEEQMatrix(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const Vector& current_charges,
        const Vector& dxi,
        const Vector& dgam,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology,
        EEQDistanceMode distance_mode = EEQDistanceMode::Topological
    );

    /**
     * @brief Solve augmented EEQ linear system with corrected parameters
     *
     * Sets up RHS: x(i) = -chi + dxi [+ CNF*sqrt(nb) if use_cnf_term=true]
     * Solves the augmented system: A*[q; lambda] = [x; total_charge]
     *
     * CRITICAL (Jan 4, 2026): CNF term ONLY in Phase 1 (topology charges)!
     * Phase 2 (final charges with dgam) does NOT use CNF term!
     * Reference: XTB gfnff_ini.f90 lines 563-570 vs 696-707
     *
     * @param A Augmented EEQ matrix from buildCorrectedEEQMatrix()
     * @param atoms Atomic numbers
     * @param cn Coordination numbers
     * @param dxi Electronegativity corrections
     * @param total_charge Total molecular charge
     * @param topology Optional topology for integer neighbor counts in CNF
     * @param use_cnf_term If true, adds CNF*sqrt(nb) to RHS (Phase 1 only!)
     * @return Vector of atomic charges (natoms elements)
     */
    Vector solveEEQ(
        const Matrix& A,
        const std::vector<int>& atoms,
        const Vector& cn,
        const Vector& dxi,
        int total_charge,
        const std::optional<TopologyInput>& topology,
        bool use_cnf_term = true,  // Default true for backward compatibility
        bool use_integer_nb = true // Phase 1 uses integer nb, Phase 2 uses fractional cn
    );

    // ===== Helper Functions =====

    /**
     * @brief Calculate coordination numbers from geometry
     *
     * Uses covalent radii with exponential distance weighting.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @return Vector of coordination numbers
     */
    Vector calculateCoordinationNumbers(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr
    ) const;

    /**
     * @brief Detect hybridization states from geometry
     *
     * Simple geometric analysis: 1=sp, 2=sp2, 3=sp3
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @return Vector of hybridization states
     */
    std::vector<int> detectHybridization(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    ) const;

    /**
     * @brief Detect atoms in pi-systems
     *
     * Simplified detection for EEQ corrections: (sp or sp2) AND (C,N,O,F,S)
     * plus picon (sp3 N,O,F bonded to sp/sp2 atom).
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param topology Topology information
     * @return Vector of pi-system membership flags
     */
    std::vector<bool> detectPiSystem(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology
    ) const;

    /**
     * @brief Detect amide nitrogens
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param is_pi_atom Pi-system membership flags
     * @param topology Topology information
     * @param cn Coordination numbers
     * @param exact_pi_membership Use true is_pi_atom[] membership for the
     *        neighbour pi-checks instead of the Fortran piadr index-cutoff
     *        replication (see .cpp for why both variants must exist)
     * @return Vector of amide nitrogen flags
     */
    std::vector<bool> detectAmideNitrogens(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::vector<bool>& is_pi_atom,
        const std::optional<TopologyInput>& topology,
        const Vector& cn,
        bool exact_pi_membership = false
    ) const;

    /**
     * @brief Detect amide hydrogens
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param is_amide Vector of amide nitrogen flags
     * @param topology Topology information
     * @return Vector of amide hydrogen flags
     *
     * Claude Generated January 2026 - Ported from Fortran gfnff_ini2.f90:1566
     */
    std::vector<bool> detectAmideHydrogens(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::vector<bool>& is_amide,
        const std::optional<TopologyInput>& topology
    ) const;

    /**
     * @brief Compute topological distances via multi-source Dijkstra (O(N·E·logN))
     *
     * Performance replacement for Floyd-Warshall (O(N³)). Uses per-atom Dijkstra
     * on the sparse bond graph with early termination at TDIST_THR cutoff.
     * OpenMP-parallelized over source atoms. Uses float32 arithmetic to match
     * Fortran real(sp) rounding behavior.
     *
     * Claude Generated - March 2026 (Performance optimization)
     *
     * @param atoms Atomic numbers
     * @param topology Topology information (neighbor lists, covalent radii)
     * @return Matrix of topological distances in Bohr (identical to Floyd-Warshall output)
     */
    Matrix computeTopologicalDistancesSparse(
        const std::vector<int>& atoms,
        const TopologyInput& topology
    ) const;

    // ===== Solve Methods =====

    /**
     * @brief Dispatch solve of augmented EEQ system using configured solver method
     * Claude Generated (March 2026): Unified solver dispatch for Phase 1 and Phase 2
     *
     * Routes to SchurCholesky, PCG, LU, or Auto based on m_solve_method.
     * Returns atomic charges (natoms elements), or fallback charges on failure.
     *
     * @param A Full augmented matrix (natoms+nfrag) × (natoms+nfrag)
     * @param x Full RHS vector (natoms+nfrag)
     * @param natoms Number of atoms
     * @param nfrag Number of fragments
     * @param total_charge Total molecular charge (for fallback)
     * @return Vector of atomic charges (natoms elements)
     */
    Vector dispatchSolve(
        const Matrix& A,
        const Vector& x,
        int natoms,
        int nfrag,
        int total_charge,
        CxxThreadPool* pool = nullptr,   // Claude Generated (WP2, May 2026): pool for Stage-4 batched parallelisation
        int num_threads = 1,              // Default 1 keeps existing callers (Phase 1, solveEEQ) serial
        const char* phase = ""            // label for the verbosity-1 "EEQ solver:" line
    );

    /**
     * @brief Solve augmented EEQ system via Schur complement + Cholesky
     *
     * Exploits the SPD structure of the NxN Coulomb sub-matrix A:
     * 1. Cholesky factorize A (O(N³/6) vs O(N³/3) for LU)
     * 2. Solve A·z₁ = b and A·z₂ = Cᵀ via back-substitution
     * 3. Form Schur complement S = C·z₂ (nfrag × nfrag)
     * 4. Solve constraint: λ = S⁻¹·(C·z₁ - d)
     * 5. Final charges: q = z₁ - z₂·λ
     *
     * Falls back to LU if Cholesky fails (matrix not SPD).
     *
     * Claude Generated - March 2026 (Performance optimization)
     *
     * @param A_nn NxN Coulomb+hardness sub-matrix (must be SPD)
     * @param rhs_atoms RHS vector for atoms (N elements)
     * @param C Constraint matrix (nfrag × N)
     * @param rhs_constraints RHS for constraints (nfrag elements)
     * @param natoms Number of atoms
     * @param nfrag Number of fragments
     * @return Charge vector (N elements), or empty vector on failure
     */
    /**
     * @brief Constraint-projected preconditioned CG for the fragment-constrained EEQ system.
     *
     * Solves  min 1/2 q^T A_nn q - b^T q  s.t.  C q = qfrag  (KKT: A_nn q + C^T lambda = b) with
     * a single CG run on the constraint tangent space: the residual and the Jacobi-preconditioned
     * direction are projected with P = I - C^T (C C^T)^{-1} C, which for the 0/1 fragment
     * membership matrix C is just "subtract the per-fragment mean" (O(N)). Cost per iteration is
     * one dense matvec, independent of nfrag — the Schur route needs nfrag+1 solves (O(N^2 nfrag)),
     * which is what made many-fragment boxes (solvent) expensive.
     *
     * Not bit-identical to the direct solve: charges are converged to |P r| <= tol (|b|+1),
     * i.e. an approximation controlled by eeq_ppcg_tol (default 1e-10, energy error ~1e-10 Eh).
     * Warm-started from the previous step's charges (re-projected onto the current qfrag).
     * Returns an empty vector if not converged within eeq_ppcg_max_iter (caller falls back).
     * Claude Generated (Sep 2026), machine-tested on a 3000-atom water box.
     */
    Vector solveWithProjectedPCG(
        const MatrixCRef& A_nn,
        const Vector& rhs_atoms,
        const MatrixCRef& C,
        const Vector& rhs_constraints,
        int natoms,
        int nfrag);

public:
    /**
     * @brief Build the local-inverse preconditioner (see LocalInversePreconditioner).
     * @param geometry_bohr N x 3 coordinates in Bohr (neighbourhoods)
     * @param frag          fragment index per atom, 0..nfrag-1
     * @param fill_block    fills the dense block A[idx, idx] (lower triangle suffices) into
     *                      the column-major n x n buffer `out`; called concurrently from
     *                      several threads, must be read-only on shared data
     * @param radius_bohr   R: entries kept per row; ext_bohr: R_ext, block that is inverted
     * @return always non-null; ->ok tells whether it may be used
     * Claude Generated (Oct 2026). Result is bit-identical for any thread count.
     */
    static std::shared_ptr<LocalInversePreconditioner> buildLocalInverse(
        const Matrix& geometry_bohr, const std::vector<int>& frag, int nfrag,
        const std::function<void(const std::vector<int>& idx, double* out)>& fill_block,
        double radius_bohr, double ext_bohr, CxxThreadPool* pool, int num_threads,
        bool want_dense_S_inverse = false);

    /// Statistics of the local-inverse preconditioner (test and benchmark hooks).
    long localInverseBuilds() const { return m_localinv_builds; }
    long localInverseBuildFailures() const { return m_localinv_build_failures; }
    long localInverseRuntimeFallbacks() const { return m_localinv_runtime_fallbacks; }

private:

    Vector solveWithSchurCholesky(
        const MatrixCRef& A_nn,
        const Vector& rhs_atoms,
        const MatrixCRef& C,
        const Vector& rhs_constraints,
        int natoms,
        int nfrag
    );

    /**
     * @brief Solve A·x = b via Preconditioned Conjugate Gradient
     *
     * Iterative O(N²·k) solver where k is iteration count (typically 10-50 with warm start).
     * Uses Jacobi (diagonal) preconditioner.
     *
     * Claude Generated - March 2026 (Performance optimization)
     *
     * @param A NxN SPD matrix
     * @param b RHS vector
     * @param x0 Initial guess (warm start from previous step)
     * @param max_iter Maximum CG iterations
     * @param tol Convergence tolerance on residual norm
     * @param block_pc Optional block-Jacobi preconditioner (Stage 2, May 2026).
     *                 When non-null, replaces the diagonal Jacobi with per-fragment
     *                 Cholesky inversion — typically reduces k from 30-100 to 2-5.
     * @return Solution vector
     */
    struct BlockJacobiPC;  // Forward decl — defined below
    Vector solveWithPCG(
        const MatrixCRef& A,
        const Vector& b,
        const Vector& x0,
        int max_iter,
        double tol,
        const BlockJacobiPC* block_pc = nullptr
    );

    /**
     * @brief Block-Jacobi preconditioner for multi-fragment EEQ PCG.
     *
     * Stage 2 (May 2026): Replaces the diagonal Jacobi `M_inv = diag(A)^-1` with
     * a block-diagonal `M^-1 = diag(A_ff^-1)`, where A_ff is the sub-block of A_nn
     * for fragment f's atoms. Each block is factorized via dense Cholesky once.
     *
     * Cost:
     *   - Build:  O(Σ atoms_f³)  ≤ O(N · max_frag²)  — negligible for small fragments.
     *   - Apply:  O(Σ atoms_f²)  ≤ O(N · max_frag)   — also negligible.
     *
     * Convergence: PCG sees only the inter-fragment Coulomb coupling as residual,
     * since intra-fragment is exactly inverted. For weakly coupled fragments
     * (e.g. solvated clusters) k typically drops to 2-5 iterations.
     *
     * Falls back to identity (no-op) if any per-fragment block is not SPD.
     */
    struct BlockJacobiPC {
        std::vector<Eigen::LLT<Matrix>> blocks;
        std::vector<std::vector<int>> atom_ids;  // atom_ids[f] = global indices of fragment f's atoms
        bool valid = false;

        /// Apply z = M^-1 r blockwise. Atoms not in any fragment fall back to identity (z = r).
        Vector apply(const Vector& r) const;

        /// Apply Z = M^-1 R to all columns of R at once (one multi-RHS triangular solve per fragment).
        /// Not bit-identical to apply() per column (blocked kernel). Claude Generated (Oct 2026).
        Matrix applyColumns(const Matrix& R) const;
    };

    /// Build BlockJacobiPC from A_nn and the fragment constraint matrix C (nfrag × natoms).
    /// Returns a default-constructed (invalid) PC if any fragment block fails Cholesky.
    static BlockJacobiPC buildBlockJacobi(const MatrixCRef& A_nn, const MatrixCRef& C);

    /**
     * @brief Multi-RHS Block-PCG solver — solves A·X = B with B ∈ R^{N × m} simultaneously.
     *
     * Stage 3 (May 2026): Replaces the per-fragment PCG loop in dispatchSolve. The single
     * matvec `A · P` becomes a GEMM (BLAS-3) over an N × m matrix instead of m independent
     * BLAS-2 matvecs, sharing memory bandwidth and benefiting from SIMD/cache reuse.
     *
     * Each column iterates with its own α, β, residual; convergence is checked per column.
     * No deflation — once any column has converged below tolerance, it continues to be
     * updated harmlessly (its residual stays small) until all columns finish or max_iter
     * is reached. Eigen handles the column-wise scalar broadcasting via .array() ops.
     *
     * @param A        N×N SPD matrix (A_nn for EEQ)
     * @param B        N×m matrix of right-hand-sides (column 0 = atomic RHS, cols 1..m-1 = C^T)
     * @param X0       N×m initial guess (zeros or warm-start matrix)
     * @param max_iter Maximum CG iterations
     * @param tol      Convergence tolerance on max-column residual norm
     * @param block_pc Optional block-Jacobi preconditioner (applied column-wise)
     * @return         N×m solution matrix X with A·X ≈ B
     */
    Matrix solveWithPCG_multiRHS(
        const MatrixCRef& A,
        const Matrix& B,
        const Matrix& X0,
        int max_iter,
        double tol,
        const BlockJacobiPC* block_pc = nullptr
    );

    // ===== Configuration =====

    ConfigManager m_config;           ///< Configuration manager
    int m_max_iterations;             ///< Maximum iterations for Phase 2 refinement
    double m_convergence_threshold;   ///< Convergence threshold for charge changes (e)
    int m_verbosity;                  ///< Verbosity level (0-3)
    bool m_calculate_cn;              ///< Auto-calculate CN if not provided
    bool m_allow_unconverged;         ///< Allow continuing with unconverged charges (Claude Generated Apr 2026)
    bool m_skip_phase2;               ///< Skip Phase 2 and use Phase 1 topology charges directly (Claude Generated Apr 2026)
    EEQSolveMethod m_solve_method;    ///< Linear solve algorithm selection
    double m_eeq_distance_cutoff_override = -1.0;  ///< WP-S3: post-init cutoff override (-1 = use config)
    double m_refactor_eps = 0.0;        ///< WP-EEQ-Cache: max displacement (Bohr) before re-factorizing; 0 = cache off (default since Sep 27, 2026)
    int    m_refactor_force_every = 0;  ///< WP-EEQ-Cache: force refactorization every N steps (0 = disabled)
    int    m_refine_iters = 1;          ///< A4: iterative-refinement steps on a cached-factor solve (0 = off)
    // Projected PCG (many-fragment path, Sep 2026): see solveWithProjectedPCG()
    int    m_ppcg_min_nfrag = 1;        ///< auto-select PPCG when nfrag >= this (0 = never auto)
    int    m_ppcg_min_atoms = 500;      ///< ... and natoms >= this
    double m_ppcg_tol = 1e-12;          ///< relative residual tolerance (|P r| <= tol * (|b|+1))
    int    m_ppcg_max_iter = 500;       ///< iteration cap; non-convergence falls back to Schur-Cholesky
    Vector m_ppcg_last_q;               ///< warm start (previous step's charges)
    // Local-inverse preconditioner of the projected PCG (Claude Generated, Oct 2026), see
    // LocalInversePreconditioner. On by default (eeq_ppcg_precond = localinv), diag restores Jacobi.
    bool   m_farfield_shortcut = true;  ///< from farfield_shortcut (default true), see eeq_solver.cpp
    bool   m_ppcg_localinv = false;     ///< set from eeq_ppcg_precond in the constructor (default localinv)
    double m_localinv_radius = 8.0;     ///< eeq_localinv_radius_bohr (kept entries per row)
    double m_localinv_ext = 10.0;       ///< eeq_localinv_ext_bohr (block that is inverted)
    double m_localinv_rebuild = 1.0;    ///< eeq_localinv_rebuild_bohr (0 = rebuild every solve)
    std::shared_ptr<LocalInversePreconditioner> m_localinv;  ///< last build (also a failed one)
    const Matrix* m_ppcg_geometry = nullptr;  ///< Phase-2 geometry (Bohr) during dispatchSolve; null in Phase 1
    const std::vector<int>* m_ppcg_atom_frag = nullptr;  ///< Phase-2 atom -> fragment (0-based, -1 none) during dispatchSolve
    CxxThreadPool* m_ppcg_pool = nullptr;     ///< pool for the parallel build (dispatchSolve)
    int    m_ppcg_threads = 1;
    long   m_localinv_builds = 0;            ///< successful builds
    long   m_localinv_build_failures = 0;    ///< builds rejected (G or S not positive definite)
    long   m_localinv_runtime_fallbacks = 0; ///< solves restarted with Jacobi (rz <= 0 or pAp <= 0)
    double m_localinv_build_ms_total = 0.0;  ///< wall time of all build attempts
    double m_matrix_rebuild_eps = 0.0;  ///< WP-EEQ-Matrix-Cache: max displacement before A_nn off-diag rebuild (0 = disabled)

    // ===== Cached Data for Energy Calculation =====

    // EEQ Solver intelligent caching for performance optimization

    mutable Vector m_dxi_stored;      ///< Stored dxi corrections from last calculateCharges() call
    // PCG warm-start cache for iterative EEQ solve
    // Claude Generated - March 2026 (Performance optimization)
    // May 2026: m_pcg_last_Z2 promoted from Vector to Matrix so all nfrag constraint
    // solves warm-start from their previous result, not just fragment 0.
    mutable Vector m_pcg_last_z1;  ///< Previous A⁻¹·b solution (warm start for PCG)
    mutable Matrix m_pcg_last_Z2;  ///< Previous A⁻¹·C^T columns (N × nfrag), one warm start per fragment
    mutable bool m_pcg_cache_valid = false;  ///< Whether PCG warm-start cache is usable

    // Multi-step PCG warm-start extrapolation (opt-in; generalises the 1-step cache above).
    // History of the converged solution matrix X = [z1 | Z2] (N × (nfrag+1)), most-recent at
    // front(); the PCG initial guess X0 is extrapolated from it (charge_extrapolation.h). The
    // converged charges are unchanged (PCG converges regardless of X0) — only iterations drop.
    // Claude Generated.
    std::string m_eeq_extrapolation = "none";   ///< none | aspc | gauss
    int         m_eeq_extrap_order  = 3;         ///< ASPC order / Gauss degree
    mutable std::deque<Matrix> m_pcg_x_history;  ///< past [z1|Z2] solutions, front=newest
    mutable bool m_eeq_extrap_cited = false;     ///< one-shot citation guard
    mutable long m_pcg_extrap_count = 0;         ///< #PCG solves that used the multi-step warm-start

public:
    /// Number of PCG solves whose initial guess was built by multi-step extrapolation
    /// (0 when eeq_extrapolation=none or before the history is long enough). Test hook.
    long pcgExtrapolationCount() const { return m_pcg_extrap_count; }

    /// Fail-loud status (F-Q4, Claude Generated): true iff the most recent solve fell
    /// back to uniform/placeholder charges (generateFallbackCharges). The caller must
    /// clearSolveStatus() before a solve and check lastSolveFailed() after, so a wrong
    /// charge set is never silently propagated into the Coulomb energy/gradient.
    bool lastSolveFailed() const { return m_last_solve_failed; }
    void clearSolveStatus() const { m_last_solve_failed = false; }

private:


    // WP7-B (May 2026): warn-once flag for CPU "batched" → cholesky fallback
    mutable bool m_batched_cpu_warned = false;

    // F-Q4 (Claude Generated): set by generateFallbackCharges when a solve degrades to
    // uniform/placeholder charges; read via lastSolveFailed(), reset via clearSolveStatus().
    mutable bool m_last_solve_failed = false;

    // Oct 2026: energy criterion for the explicit batched solver (replaces the old minimum
    // inter-fragment distance, eeq_batched_min_distance, and its per-step O(N^2) contact scan).
    double m_batched_max_error = 1e-6;                 ///< eeq_batched_max_error_eh (0 = no check)
    mutable std::vector<std::string> m_batched_checked; ///< "<phase>:<natoms>:<nfrag>" already estimated
    mutable double m_last_batched_error = -1.0;         ///< last estimate (Eh), -1 = none

    // Convergence statistics — accumulate across calls, print summary instead of per-call spam
    // Claude Generated (April 2026)
    mutable int m_pcg_total_calls = 0;       ///< Total PCG solve calls
    mutable int m_pcg_nonconv_calls = 0;    ///< PCG calls that did not converge
    mutable int m_pcg_total_iters = 0;      ///< Total PCG iterations across all calls
    mutable double m_pcg_worst_residual = 0.0; ///< Worst |r| among non-converged calls

    // Claude Generated (Apr 2026, corrected Jul 2026): hint that Phase 2 has previously
    // produced implausible charges for this system size.
    //
    // This used to be a permanent freeze to the Phase-1 charges ("skip Phase 2 forever").
    // That was removed in Jun 2026 (see docs/GFNFF_POLARIZATION_AUDIT.md) — the flag is
    // now TRANSIENT: it is cleared on the next plausible Phase-2 solve (:3959), and its
    // only remaining effect is to steer dispatchSolve to the Batched solver, and only for
    // nfrag >= 8 (:1468-1478). It no longer suppresses Phase 2.
    mutable bool m_phase2_historically_implausible = false;
    mutable int m_phase2_implausible_natoms = 0;

    // ===== Pre-allocated Buffers for calculateFinalCharges (Claude Generated Mar 2026) =====
    // Avoid ~26 MB alloc+free per gradient step for large molecules (N=1280).
    // Buffers are allocated once and reused when atom count stays the same.
    //
    // WP-FF-Packed-Triangular (May 2026): m_phase2_distances is packed lower-triangular
    //   (N(N+1)/2 doubles, half of full N×N). Index: data[i*(i+1)/2 + j] for i > j.
    //   Diagonal slot exists but is never written (atom self-distance is meaningless;
    //   Coulomb diagonal is built from gam_corrected + TSQRT2PI/sqrt(alpha) separately).
    //   Matches the convention used by GFNFF::triIdx / m_shared_srab / m_external_srab.
    mutable Eigen::VectorXd m_phase2_distances;  ///< packed lower-triangular distance buffer
    mutable Matrix m_phase2_A;           ///< (N+nfrag)×(N+nfrag) augmented matrix buffer
    mutable Vector m_phase2_rhs;         ///< (N+nfrag) RHS vector buffer
    mutable int m_phase2_buf_natoms = 0; ///< Atom count for current buffer size
    mutable int m_phase2_buf_nfrag = 0;  ///< Fragment count for current buffer size

    // WP-EEQ-Cache (May 2026): Cholesky factor cache across MD/opt steps.
    // Avoids O(N^3/6) refactorization when geometry displacement < m_refactor_eps.
    struct EEQCholeskyCache {
        Eigen::MatrixXd chol_factor;  ///< Cached Cholesky factor of A_nn (column-major, threaded LAPACK dpotrf; see eeq_solver.cpp)
        Matrix Z2;                 ///< A_nn^{-1}·C^T  (N x nfrag) — reused across steps
        Matrix S;                  ///< Schur complement C·Z2 (nfrag x nfrag)
        Matrix last_geometry;      ///< Geometry (N x 3) at factorization time
        Vector last_cn;            ///< CN vector at factorization time
        int steps_since_refactor = 0;
        int cached_natoms = 0;
        int cached_nfrag  = 0;
        bool valid = false;
        void reset() { valid = false; steps_since_refactor = 0; }
    };
    mutable EEQCholeskyCache m_chol_cache;
    mutable Matrix m_pending_geometry;  ///< Set by calculateFinalCharges() before dispatchSolve
    mutable Vector m_pending_cn;        ///< Set by calculateFinalCharges() before dispatchSolve

    // WP-EEQ-Matrix-Cache (May 2026): A_nn off-diagonal cache across MD steps.
    // The Coulomb off-diagonal A_ij = erf(γ_ij·r_ij)/r_ij depends on r_ij and on
    // alpha_corrected (a function of topology_charges + CN). When geometry change
    // < m_matrix_rebuild_eps AND CN change < 0.05, the off-diagonal is reusable;
    // only the diagonal needs rebuilding each step (cheap O(N)). Default disabled
    // (eps=0) for safety; user opt-in for performance.
    mutable Matrix m_cached_A_offdiag;     ///< N×N off-diagonal of A_nn at cache time
    mutable Matrix m_cached_A_geometry;    ///< Geometry (N×3) at cache time
    mutable Vector m_cached_A_cn;          ///< CN vector at cache time
    mutable int    m_cached_A_natoms = 0;  ///< Atom count for cache validity check
    mutable bool   m_A_cache_valid = false;

    /// WP-FF-DistMatrix-Sharing (May 2026): external packed-triangular srab from GFNFF.
    /// When non-null, calculateFinalCharges fills m_phase2_distances from this instead of
    /// recomputing the O(N^2) sqrt loop.
    const Eigen::VectorXd* m_external_srab = nullptr;
    const bool* m_external_srab_too_close = nullptr; ///< atoms-too-close flag of m_external_srab (owner GFNFF)

    /// WP5 (June 2026): implicit-solvation reaction-field Born matrix (nat×nat),
    /// added to the EEQ matrix in calculateFinalCharges. nullptr = gas phase.
    /// Column-major to match ALPBSolvation::bornMatrix(); symmetric, so the storage
    /// order vs the row-major augmented EEQ matrix is irrelevant.
    const Eigen::MatrixXd* m_reaction_field = nullptr;

    /// Last truly successful charges (Phase 2 if available, otherwise Phase 1).
    /// Initialized from topology_charges on first call, then overwritten with the Phase 2
    /// result after a successful solve (line 2919). NOT overwritten on plausibility-check
    /// failure — the cache preserves the last good result. Used as fallback when the solver
    /// fails on subsequent steps. Validated on read: allFinite(), |q|_max < 50, charge sum
    /// within tolerance, to avoid using a corrupted cache.
    mutable Vector m_last_successful_charges;

    /// Validate cached charges are plausible for use as fallback
    bool isValidChargeCache(const Vector& charges, int natoms, int total_charge) const {
        return charges.size() == natoms
            && charges.allFinite()
            && std::abs(charges.sum() - total_charge) < 1.0
            && charges.cwiseAbs().maxCoeff() < 50.0;
    }


    /// Ensure buffers are large enough. Only reallocates if size changed.
    void ensurePhase2Buffers(int natoms, int nfrag) const {
        if (natoms != m_phase2_buf_natoms || nfrag != m_phase2_buf_nfrag) {
            int m = natoms + nfrag;
            // WP-FF-Packed-Triangular: packed lower-triangular = N(N+1)/2 doubles.
            m_phase2_distances.resize(static_cast<Eigen::Index>(natoms) * (natoms + 1) / 2);
            m_phase2_A.resize(m, m);
            m_phase2_rhs.resize(m);
            m_phase2_buf_natoms = natoms;
            m_phase2_buf_nfrag = nfrag;
        }
    }
};

// ===== Parameter Definitions =====
// Claude Generated (Oct 2026): the PARAMs of this solver live in the gfnff block of gfnff.h
// (category EEQ); the former registry module eeq_solver was merged into gfnff. EEQSolver still
// reads the JSON scope "eeq_solver" (ConfigManager("eeq_solver", ...)), and
// ParameterRegistry::canonicalModule() answers registry lookups by that name from gfnff, so
// defaults, aliases and the dotted -eeq_solver.<name> spelling are unchanged.
