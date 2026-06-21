/*
 * <EEQSolverHip - CUDA-Accelerated EEQ Charge Solver>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU EEQ solver using cuSOLVER Cholesky.
 * Builds the N×N Coulomb matrix on GPU via k_eeq_build_matrix kernel,
 * then solves via cusolverDnDpotrf (Cholesky) + cusolverDnDpotrs (triangular solve).
 *
 * Interface is CUDA-agnostic (Pimpl pattern) — can be forward-declared
 * from non-CUDA translation units.
 *
 * Minimum compute capability: 6.0 (Pascal) — native double atomicAdd.
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF EEQ)
 */

#pragma once

#ifdef USE_ROCM_GFNFF

#include <memory>
#include <vector>

struct EEQSolverHipImpl;

/**
 * @brief GPU-accelerated EEQ Coulomb matrix build + Cholesky solve.
 *
 * Usage:
 *   EEQSolverHip eeq_gpu(max_natoms);
 *   bool ok = eeq_gpu.solve(N, nfrag, cx, cy, cz, alpha, gam, fraglist,
 *                            rhs_atoms, rhs_constr, z1, Z2);  // cx/cy/cz: SoA device ptrs
 *   if (!ok) { ... CPU fallback ... }
 *   // CPU Schur complement: q = z1 - Z2 * S^{-1} * (C*z1 - d)
 */
class EEQSolverHip {
public:
    explicit EEQSolverHip(int max_natoms);
    ~EEQSolverHip();

    EEQSolverHip(const EEQSolverHip&) = delete;
    EEQSolverHip& operator=(const EEQSolverHip&) = delete;

    /**
     * @brief Build N×N Coulomb matrix on GPU + Cholesky solve via cuSOLVER.
     *
     * @param natoms           Number of atoms N
     * @param nfrag            Number of fragments (typically 1)
     * @param cx               GPU pointer to [N] x-coordinates (Bohr, SoA)
     * @param cy               GPU pointer to [N] y-coordinates (Bohr, SoA)
     * @param cz               GPU pointer to [N] z-coordinates (Bohr, SoA)
     * @param alpha_corrected  Host [N] charge-corrected alpha² values (already squared!)
     * @param gam_corrected    Host [N] corrected hardness values
     * @param fraglist         Host [N] fragment ID per atom (1-indexed), empty for single-fragment
     * @param rhs_atoms        Host [N] RHS vector for atoms
     * @param rhs_constraints  Host [nfrag] target charges per fragment
     * @param out_z1           Host [N] output: A⁻¹ · b_atoms
     * @param out_Z2           Host [N*nfrag] output: A⁻¹ · C^T columns (column-major)
     * @return true on success, false if Cholesky fails (not SPD → caller falls back to CPU)
     */
    /**
     * @param force_refactor   If true: always rebuild matrix + dpotrf (full solve).
     *                         If false: reuse cached Cholesky factor L (lazy solve).
     *                         Caller decides based on geometry RMSD from last refactorization.
     *                         Default: true (always refactorize = matches reference).
     */
    bool solve(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* alpha_corrected,
        const double* gam_corrected,
        const std::vector<int>& fraglist,
        const double* rhs_atoms,
        const double* rhs_constraints,
        double* out_z1,
        double* out_Z2,
        double cutoff_sq = 0.0,
        bool force_refactor = true
    );

    /**
     * @brief WP2 variant: alpha, gam, rhs_atoms already on GPU.
     *
     * Eliminates per-step H2D upload of alpha_corrected, gam_corrected, and rhs_atoms.
     * alpha and gam come from uploadEEQTopologyParams (topology-constant);
     * rhs_atoms comes from k_build_eeq_rhs (written on main stream each step).
     * Caller must synchronize the main workspace stream before calling
     * (finalizeCNForCPU() satisfies this requirement).
     *
     * @param d_alpha_corrected  Device [N] alpha² (from FFWorkspaceGPU::getDeviceAlphaPtr())
     * @param d_gam_corrected    Device [N] gam+dgam (from FFWorkspaceGPU::getDeviceGamPtr())
     * @param d_rhs_atoms        Device [N] RHS from k_build_eeq_rhs
     * @param d_rhs_constraints  Device [nfrag] target charges (unused in GPU solve, for WP4)
     */
    bool solveWithDeviceRHS(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* d_alpha_corrected,
        const double* d_gam_corrected,
        const double* d_rhs_atoms,
        const double* d_rhs_constraints,
        const std::vector<int>& fraglist,
        double* out_z1,
        double* out_Z2,
        double cutoff_sq = 0.0,
        bool force_refactor = true
    );

    /// ⚠️ NOT USED in production (kept for reference). ~4 s slower than solve() + CPU Schur.
    bool solveAndComputeCharges(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* alpha_corrected,
        const double* gam_corrected,
        const std::vector<int>& fraglist,
        const double* rhs_atoms,
        const double* rhs_constraints,
        double* out_charges,
        double cutoff_sq = 0.0
    );

    /**
     * @brief WP5-A: WP2 solve + GPU-side Schur complement (no D2H for z1/Z2).
     *
     * Replaces the 22 KB D2H + CPU O(N) Schur loop + 11 KB H2D with:
     *   - k_eeq_reduce_sums (sum z1, Z2 on GPU)
     *   - 2-scalar D2H (16 bytes) for lambda
     *   - k_eeq_schur_nfrag1 (charges stored in d_rhs[0..N-1])
     *
     * On success: charges available via getDeviceChargesPtr().
     * Caller must NOT call setEEQCharges() afterwards — use setEEQDeviceCharges() instead.
     *
     * Only supported for nfrag == 1 (fast path). Returns false for nfrag > 1
     * or if Cholesky fails; caller must fall back to CPU solve + setEEQCharges().
     *
     * @param rhs_c0         m_eeq_rhs_constraints[0] — target fragment charge (scalar)
     * @param force_refactor See solveWithDeviceRHS()
     */
    bool solveWithDeviceRHSAndGPUSchur(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* d_alpha_corrected,
        const double* d_gam_corrected,
        const double* d_rhs_atoms,
        const std::vector<int>& fraglist,
        double rhs_c0,
        double cutoff_sq = 0.0,
        bool force_refactor = true
    );

    /// WP5-A: Device pointer to EEQ charges after successful solveWithDeviceRHSAndGPUSchur().
    /// Points to d_rhs[0..N-1] (valid until next solve call).
    /// EEQ stream is fully synced before return — safe to use immediately on any stream.
    double* getDeviceChargesPtr();

    /// Step B (May 2026): true if the last successful refactorization used the
    /// LU fallback (cusolverDnDgetrf) because Cholesky failed on an indefinite
    /// matrix. Persists across lazy solves until the next refactorization.
    /// Used by gfnff_gpu_method.cpp to log a one-line warning matching the
    /// CPU dispatcher's "ill-conditioned matrix → using LU" message.
    bool isUsingLUFallback() const;

    /// Step B (May 2026): leading minor index returned by cusolverDnDpotrf when
    /// the most recent Cholesky attempt failed. 0 if Cholesky succeeded.
    int getLastCholInfo() const;

    /**
     * @brief WP7-A: General Schur for nfrag > 1 (May 2026), full GPU-resident path.
     *
     * Generalizes solveWithDeviceRHSAndGPUSchur() to multi-fragment systems:
     *   - Single N×N Cholesky (cross-fragment Coulomb correctly retained).
     *   - dpotrs solves [b_atoms | C^T] in one call (nfrag+1 RHS columns).
     *   - GPU reduction yields Cz1[nfrag] + S[nfrag×nfrag] (8·nfrag·(nfrag+1) bytes D2H).
     *   - Tiny CPU Gauss-elim for λ; H2D upload (8·nfrag bytes).
     *   - One apply kernel produces charges in d_rhs[0..N-1].
     *
     * uploadFragmentTopology() must have been called first (provides d_atom_frag and
     * d_rhs_constraint_cols). Returns false on Cholesky failure or invalid topo —
     * caller must fall back to solveWithDeviceRHS() + CPU Schur.
     *
     * On success: charges available via getDeviceChargesPtr() (same contract as WP5-A).
     */
    bool solveWithDeviceRHSAndGPUSchurGeneral(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* d_alpha_corrected,
        const double* d_gam_corrected,
        const double* d_rhs_atoms,
        const std::vector<int>& fraglist,
        const std::vector<double>& rhs_constraints,  ///< [nfrag] target charge per fragment
        double cutoff_sq = 0.0,
        bool force_refactor = true
    );

    // ── WP6: Batched per-fragment Cholesky (nfrag > 1) ──────────────────────

    /**
     * @brief Upload fragment topology for the batched Cholesky path (nfrag > 1).
     *
     * Claude Generated (May 2026): Topology-constant; call once per topology build
     * before any solveWithDeviceRHSAndGPUSchurBatched() call.
     * Computes fragment sizes, A-block offsets, sorted-atom permutation on CPU,
     * uploads all index arrays to GPU, pre-fills constraint cols (all-ones) in
     * d_rhs_blocks, and queries the cuSOLVER workspace for the largest fragment.
     *
     * @param nfrag    Number of fragments
     * @param fraglist [N] fragment ID per atom (1-indexed)
     * @param natoms   Total atom count N
     */
    void uploadFragmentTopology(int nfrag,
                                const std::vector<int>& fraglist,
                                int natoms);

    /// Returns true after a successful uploadFragmentTopology() call.
    bool isFragmentTopoValid() const;

    /**
     * @brief WP7-B: cache the minimum atom-atom distance between different fragments.
     *
     * Claude Generated (May 2026): O(N²) scan over the fragment-sorted atom map
     * (uploadFragmentTopology() must have been called first). Computes the minimum
     * |r_i - r_j|² over all pairs (i,j) with different fragment IDs and caches the
     * squared result. Used by the dispatch in gfnff_gpu_method.cpp to decide whether
     * to warn before running the batched (cross-fragment-Coulomb-dropping) solver.
     *
     * Coordinates are CPU-side Bohr; this scan happens once per topology build, not
     * per MD step. If never called, getMinFragmentDistanceSq() returns -1.
     */
    void updateMinFragmentDistance(const double* host_x,
                                   const double* host_y,
                                   const double* host_z,
                                   int natoms);

    /// WP7-B: cached squared min inter-fragment distance (Bohr²). -1 if not yet computed.
    double getMinFragmentDistanceSq() const;

    /**
     * @brief WP7-C: GPU PCG solver for nfrag>=1 (May 2026).
     *
     * Iterative O(k·N²) replacement for the WP7-A Cholesky path. Reuses WP7-A's
     * Schur-complement infrastructure: the (1+nfrag) PCG runs solve A·z1 = b and
     * A·Z2[:,f] = e_f independently, then k_eeq_reduce_fragment_sums + CPU Gauss-elim
     * + k_eeq_schur_general apply the constraint exactly. Warm-start uses
     * d_z1_persistent / d_Z2_persistent (kept across calls until topology changes).
     *
     * Returns false on PCG stall (max_iter reached without |r|<tol) — caller must
     * fall back to solveWithDeviceRHSAndGPUSchurGeneral (WP7-A).
     *
     * @param max_iter        Per-PCG-call iteration cap (e.g. 200).
     * @param tol             Convergence tolerance on |r| (e.g. 1e-10).
     * @param force_refactor  If true: rebuild A and re-extract Jacobi precond.
     */
    bool solveWithDeviceRHSAndGPUPCG(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* d_alpha_corrected,
        const double* d_gam_corrected,
        const double* d_rhs_atoms,
        const std::vector<int>& fraglist,
        const std::vector<double>& rhs_constraints,
        int    max_iter,
        double tol,
        double cutoff_sq    = 0.0,
        bool   force_refactor = true
    );

    /**
     * @brief Batched per-fragment EEQ: independent N_f×N_f Cholesky + host Schur per fragment.
     *
     * Claude Generated (May 2026): Replaces single N×N Cholesky for nfrag > 1 systems.
     * Each fragment is solved independently (cross-fragment Coulomb = 0).
     * On success: charges for all N atoms available via getDeviceChargesPtr()
     *   (in global atom order, in d_rhs[0..N-1]).
     * Returns false if any fragment's Cholesky fails → caller falls back to CPU.
     *
     * uploadFragmentTopology() must have been called first.
     *
     * @param natoms            Total atom count N
     * @param nfrag             Number of fragments
     * @param cx/cy/cz          Device SoA coordinate pointers (global atom order)
     * @param d_alpha           Device [N] alpha² (from uploadEEQTopologyParams)
     * @param d_gam             Device [N] gam_corrected (from uploadEEQTopologyParams)
     * @param d_rhs_atoms       Device [N] RHS from k_build_eeq_rhs
     * @param rhs_constraints   Host [nfrag] target charge per fragment
     * @param cutoff_sq         Distance cutoff² (0 = no cutoff)
     * @param force_refactor    If false: reuse cached Cholesky factors (lazy solve)
     */
    bool solveWithDeviceRHSAndGPUSchurBatched(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* d_alpha,
        const double* d_gam,
        const double* d_rhs_atoms,
        const std::vector<double>& rhs_constraints,
        double cutoff_sq = 0.0,
        bool force_refactor = true);

private:
    std::unique_ptr<EEQSolverHipImpl> m_impl;
    int m_max_natoms;
    int m_last_N = 0;  ///< cached N for workspace reuse
};

#endif // USE_ROCM_GFNFF
