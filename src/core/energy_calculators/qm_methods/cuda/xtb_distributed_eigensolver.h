/*
 * < Multi-GPU dense symmetric eigensolver for the native xTB SCF >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

// Claude Generated (Sep 2026, multi-GPU step 3).
//
// The per-iteration eigensolve of the (Cholesky-reduced) Fock matrix is the O(n^3) wall of a
// large GFN1/GFN2 SCF. This class spreads that one step over several GPUs; everything else in the
// device-resident SCF stays on the context's own device. The caller hands over the reduced matrix
// on its device, the solver distributes the columns, solves, and writes the eigenvectors back into
// the same buffer and the eigenvalues into a device vector - a drop-in for cusolverDn?syevd.
//
// Backends:
//   mp - cuSOLVERMp (NCCL). ScaLAPACK-style: one rank per GPU; curcuma is a single process, so
//        every call runs one host thread per GPU on communicators from ncclCommInitAll.
//        Measured n = 15444 on 4x RTX A4500 (PCIe): FP64 50.2 -> 20.7 s, FP32 7.6 -> 5.2 s.
//   mg - cusolverMg (single-process multi-GPU, part of the CUDA toolkit, deprecated in CUDA 13).
//        Fallback: FP64 50.2 -> 26.4 s on 4 GPUs, but FP32 SLOWER than one GPU (12.6 vs 7.6 s).
//
// Data layout (both): a 1 x ndev process grid, column blocks of `block` columns dealt out
// cyclically (block b lives on device b % ndev), every device holds all n rows of its columns.

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace curcuma {
namespace xtb {
namespace gpu {

class DistributedEigensolver {
public:
    virtual ~DistributedEigensolver() = default;

    /**
     * @brief Create a solver over `devices`.
     * @param backend "auto" (mp if built with cuSOLVERMp, else mg), "mp" or "mg"
     * @param block   column block size (128 was fastest at n = 15444)
     * @param why     reason when nullptr is returned
     */
    static std::unique_ptr<DistributedEigensolver> create(const std::string& backend,
                                                          const std::vector<int>& devices,
                                                          int block, std::string& why);

    /// Backend names compiled into this plugin, e.g. "mp,mg".
    static std::string availableBackends();

    virtual const char* name() const = 0;
    virtual int deviceCount() const = 0;

    /**
     * @brief Solve A x = lambda x in place.
     * @param n          matrix order
     * @param A          n x n column-major matrix (double*, or float* when fp32) on src_device;
     *                   the lower triangle is read, all n^2 elements are overwritten by the
     *                   eigenvectors (column k = eigenvector k)
     * @param eig        length-n device vector on src_device for the ascending eigenvalues
     * @param fp32       single precision
     * @param src_device device holding A and eig; its stream must be synchronised
     * @return false on any failure (the caller then uses the single-GPU solver)
     */
    virtual bool solve(int n, void* A, void* eig, bool fp32, int src_device) = 0;

    /// Free all per-device matrices and workspaces (handles and communicators stay); the next
    /// solve re-creates them. Used before the gradient to lower the peak memory.
    virtual void releaseBuffers() {}

    /// True when solveGeneralized() is available (cuSOLVERMp backend).
    virtual bool supportsGeneralized() const { return false; }

    /**
     * @brief Solve the generalized problem F C = S C e with S = L L^T, entirely on the devices:
     *        reduction A = L^-1 F L^-T (sygst), eigensolve (syevd), back-transform C = L^-T Q (trsm).
     * @param A            F (lower triangle read) on src_device; overwritten by C (column k = vector k)
     * @param L            lower Cholesky factor of S on src_device, same precision as A. Only read
     *                     when `l_generation` differs from the previous call (or the precision or n
     *                     changed): the factor stays distributed across calls. May be nullptr when
     *                     the caller knows the factor is current.
     * @param l_generation changes whenever the caller's L changes (new geometry)
     * @return false on failure; inputIntact() tells whether A is still F
     */
    virtual bool solveGeneralized(int n, void* A, const void* L, long l_generation, void* eig, bool fp32,
                                  int src_device)
    {
        (void)n; (void)A; (void)L; (void)l_generation; (void)eig; (void)fp32; (void)src_device;
        return false;
    }

    /// True when the distributed factor for (n, fp32, l_generation) is already in place.
    virtual bool hasMetric(int n, bool fp32, long l_generation) const
    {
        (void)n; (void)fp32; (void)l_generation;
        return false;
    }

    /// Wall-clock split of the last solve() in ms: column scatter, collective solve, gather.
    /// For solveGeneralized() the solve part contains reduction + syevd + back-transform, and
    /// reduce_ms / back_ms give those two separately (0 for solve()).
    void lastTimings(double& scatter_ms, double& solve_ms, double& gather_ms) const
    {
        scatter_ms = m_scatter_ms; solve_ms = m_solve_ms; gather_ms = m_gather_ms;
    }
    void lastGeneralizedTimings(double& reduce_ms, double& back_ms) const
    {
        reduce_ms = m_reduce_ms; back_ms = m_back_ms;
    }

    /// After a failed solve(): true when A was not modified (the caller can still solve it).
    bool inputIntact() const { return m_input_intact; }

protected:
    bool m_input_intact = true;
    double m_scatter_ms = 0.0, m_solve_ms = 0.0, m_gather_ms = 0.0;
    double m_reduce_ms = 0.0, m_back_ms = 0.0;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma
