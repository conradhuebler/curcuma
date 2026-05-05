/*
 * <EEQSolverGPU - CUDA-Accelerated EEQ Charge Solver>
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

#ifdef USE_CUDA

#include <memory>
#include <vector>

struct EEQSolverGPUImpl;

/**
 * @brief GPU-accelerated EEQ Coulomb matrix build + Cholesky solve.
 *
 * Usage:
 *   EEQSolverGPU eeq_gpu(max_natoms);
 *   bool ok = eeq_gpu.solve(N, nfrag, cx, cy, cz, alpha, gam, fraglist,
 *                            rhs_atoms, rhs_constr, z1, Z2);  // cx/cy/cz: SoA device ptrs
 *   if (!ok) { ... CPU fallback ... }
 *   // CPU Schur complement: q = z1 - Z2 * S^{-1} * (C*z1 - d)
 */
class EEQSolverGPU {
public:
    explicit EEQSolverGPU(int max_natoms);
    ~EEQSolverGPU();

    EEQSolverGPU(const EEQSolverGPU&) = delete;
    EEQSolverGPU& operator=(const EEQSolverGPU&) = delete;

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
        bool force_refactor = true
    );

    /// WP5-A: Device pointer to EEQ charges after successful solveWithDeviceRHSAndGPUSchur().
    /// Points to d_rhs[0..N-1] (valid until next solve call).
    /// EEQ stream is fully synced before return — safe to use immediately on any stream.
    double* getDeviceChargesPtr();

private:
    std::unique_ptr<EEQSolverGPUImpl> m_impl;
    int m_max_natoms;
    int m_last_N = 0;  ///< cached N for workspace reuse
};

#endif // USE_CUDA
