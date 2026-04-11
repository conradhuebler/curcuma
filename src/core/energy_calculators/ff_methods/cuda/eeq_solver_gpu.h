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
    bool solve(
        int natoms, int nfrag,
        const double* cx, const double* cy, const double* cz,
        const double* alpha_corrected,
        const double* gam_corrected,
        const std::vector<int>& fraglist,
        const double* rhs_atoms,
        const double* rhs_constraints,
        double* out_z1,
        double* out_Z2
    );

private:
    std::unique_ptr<EEQSolverGPUImpl> m_impl;
    int m_max_natoms;
    int m_last_N = 0;  ///< cached N for workspace reuse
};

#endif // USE_CUDA
