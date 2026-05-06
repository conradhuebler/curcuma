/*
 * <EEQSolverGPU - CUDA-Accelerated EEQ Charge Solver>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU EEQ solver using cuSOLVER Cholesky.
 * Builds the N×N Coulomb matrix on GPU via k_eeq_build_matrix kernel,
 * then solves via cusolverDnDpotrf (Cholesky) + cusolverDnDpotrs (triangular solve).
 *
 * Matrix build: O(N²) with N*(N+1)/2 threads (lower triangle, erf() per pair)
 * Cholesky solve: O(N³/6) via cuSOLVER (2× faster than dense LU)
 * CPU Schur complement: O(nfrag) — typically scalar division for nfrag=1
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF EEQ)
 */

#include "eeq_solver_gpu.h"
#include "gfnff_soa.h"
#include "gfnff_kernels.cuh"

#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <string>
#include <vector>

// ============================================================================
// Error-checking helpers
// ============================================================================

static void checkCudaEEQ(cudaError_t err, const char* msg)
{
    if (err != cudaSuccess)
        throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
}

static void checkCusolverEEQ(cusolverStatus_t status, const char* msg)
{
    if (status != CUSOLVER_STATUS_SUCCESS)
        throw std::runtime_error(std::string(msg) + ": cusolver error " + std::to_string(static_cast<int>(status)));
}

// ============================================================================
// Pimpl implementation struct
// ============================================================================

struct EEQSolverGPUImpl {
    cusolverDnHandle_t cusolver_handle = nullptr;
    cudaStream_t stream = nullptr;

    // Device buffers
    CudaBuffer<double> d_alpha;      ///< [N] alpha_corrected
    CudaBuffer<double> d_gam;        ///< [N] gam_corrected
    CudaBuffer<double> d_A;          ///< [N*N] Coulomb matrix; after dpotrf: Cholesky factor L
    CudaBuffer<double> d_rhs;        ///< [N*(nfrag+1)] RHS matrix (column-major)
    CudaBuffer<double> d_workspace;  ///< cuSOLVER workspace
    CudaBuffer<int>    d_info;       ///< [1] cuSOLVER info
    CudaBuffer<double> d_sums;       ///< [2] GPU reduction: [Cz1, S] for Schur complement

    // Pinned host buffers for async D2H transfer
    double* h_result = nullptr;      ///< [N*(nfrag+1)] downloaded solution
    int*    h_info = nullptr;        ///< [1] cusolver status

    int workspace_size = 0;

    // Lazy Cholesky state: d_A holds the cached factor L after refactorization.
    // force_refactor flag (passed per solve() call) decides whether to rebuild.
    // Caller (GFNFFGPUComputationalMethod) decides based on geometry RMSD.
    bool m_has_cached_factor = false;

    // Cached host RHS buffer to avoid per-step heap allocation
    std::vector<double> m_h_rhs_cache;
};

// ============================================================================
// EEQ Matrix Build Kernel
// ============================================================================

/// sqrt(2/pi) constant — matches EEQSolver CPU implementation
__device__ constexpr double TSQRT2PI_GPU = 0.797884560802866;

/**
 * @brief Build N×N EEQ Coulomb matrix on GPU (column-major for cuSOLVER)
 *
 * Each thread handles one (i,j) pair in the lower triangle (i >= j).
 * Thread mapping: tid → (i,j) via inverse triangular number formula.
 *
 * Diagonal:   A(i,i) = gam[i] + sqrt(2/pi) / sqrt(alpha[i])
 * Off-diag:   gamma_ij = 1/sqrt(alpha[i] + alpha[j])
 *             A(i,j) = erf(gamma_ij * r_ij) / r_ij
 *
 * Column-major layout: A[col*N + row].
 * Symmetric matrix → writes both lower and upper triangle for cuSOLVER potrs.
 *
 * Reference: Spicher/Grimme JCTC 2020 — GFN-FF EEQ Coulomb matrix
 */
__global__ __launch_bounds__(256, 4)
void k_eeq_build_matrix(
    int N,
    const double* __restrict__ cx,       // [N] x-coordinates (Bohr, SoA)
    const double* __restrict__ cy,       // [N] y-coordinates (Bohr, SoA)
    const double* __restrict__ cz,       // [N] z-coordinates (Bohr, SoA)
    const double* __restrict__ alpha,     // [N] alpha² = (alpha_base + ff*qa)², already squared
    const double* __restrict__ gam,       // [N] gam_corrected
    double* __restrict__ A,               // [N*N] output column-major
    double cutoff_sq)                     // distance cutoff squared (0.0 = no cutoff)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int n_lower = N * (N + 1) / 2;
    if (tid >= n_lower) return;

    // Map flat tid to (i,j) in lower triangle: i >= j
    // tid = j + i*(i+1)/2  →  i = floor((-1 + sqrt(1 + 8*tid)) / 2)
    int i = static_cast<int>((-1.0 + sqrt(1.0 + 8.0 * static_cast<double>(tid))) * 0.5);
    // Correct for floating-point rounding
    while (i * (i + 1) / 2 > tid) --i;
    while ((i + 1) * (i + 2) / 2 <= tid) ++i;
    int j = tid - i * (i + 1) / 2;

    if (i == j) {
        // Diagonal: gam + sqrt(2/pi) / sqrt(alpha)
        A[i * N + i] = gam[i] + TSQRT2PI_GPU / sqrt(alpha[i]);
    } else {
        // Off-diagonal: compute interatomic distance (SoA layout)
        double dx = cx[i] - cx[j];
        double dy = cy[i] - cy[j];
        double dz = cz[i] - cz[j];
        double r_sq = dx * dx + dy * dy + dz * dz;

        // Distance cutoff: zero out near-zero long-range elements
        // Matches CPU EEQ behavior for matrix conditioning (same cutoff_sq logic)
        if (cutoff_sq > 0.0 && r_sq > cutoff_sq) {
            A[j * N + i] = 0.0;
            A[i * N + j] = 0.0;
        } else {
            double r = sqrt(r_sq);
            double gamma_ij = 1.0 / sqrt(alpha[i] + alpha[j]);
            double val = erf(gamma_ij * r) / r;

            // Column-major: element (row, col) at col*N + row
            // Write both triangles (cuSOLVER potrs reads full matrix)
            A[j * N + i] = val;   // lower triangle
            A[i * N + j] = val;   // upper triangle (symmetric)
        }
    }
}

// ============================================================================
// EEQSolverGPU Implementation
// ============================================================================

EEQSolverGPU::EEQSolverGPU(int max_natoms)
    : m_max_natoms(max_natoms)
{
    m_impl = std::make_unique<EEQSolverGPUImpl>();

    // Create cuSOLVER handle + dedicated stream
    checkCusolverEEQ(cusolverDnCreate(&m_impl->cusolver_handle), "cusolverDnCreate");
    checkCudaEEQ(cudaStreamCreate(&m_impl->stream), "cudaStreamCreate EEQ");
    checkCusolverEEQ(cusolverDnSetStream(m_impl->cusolver_handle, m_impl->stream),
                     "cusolverDnSetStream");

    // Pre-allocate device buffers for max atom count
    m_impl->d_alpha.alloc(max_natoms);
    m_impl->d_gam.alloc(max_natoms);
    m_impl->d_A.alloc(max_natoms * max_natoms);
    m_impl->d_info.alloc(1);

    // Pre-allocate d_rhs for max nrhs = max_nfrag + 1 (assume max_nfrag = 8)
    constexpr int MAX_NFRAG = 8;
    m_impl->d_rhs.alloc(max_natoms * (MAX_NFRAG + 1));

    // Pinned host buffers for async download
    checkCudaEEQ(cudaMallocHost(&m_impl->h_result,
                                 max_natoms * (MAX_NFRAG + 1) * sizeof(double)),
                 "pinned h_result");
    checkCudaEEQ(cudaMallocHost(&m_impl->h_info, sizeof(int)), "pinned h_info");
}

EEQSolverGPU::~EEQSolverGPU()
{
    if (m_impl) {
        if (m_impl->h_result) cudaFreeHost(m_impl->h_result);
        if (m_impl->h_info) cudaFreeHost(m_impl->h_info);
        if (m_impl->cusolver_handle) cusolverDnDestroy(m_impl->cusolver_handle);
        if (m_impl->stream) cudaStreamDestroy(m_impl->stream);
    }
}

bool EEQSolverGPU::solve(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* alpha_corrected,
    const double* gam_corrected,
    const std::vector<int>& fraglist,
    const double* rhs_atoms,
    const double* rhs_constraints,
    double* out_z1,
    double* out_Z2,
    double cutoff_sq,
    bool force_refactor)
{
    const int N = natoms;
    const int nrhs = nfrag + 1;  // b_atoms + nfrag constraint columns

    // Refactorize if: first call, N changed, or caller says geometry moved enough.
    // Caller (GFNFFGPUComputationalMethod) tracks RMSD from last refactorization.
    const bool do_refactor = force_refactor
                          || !m_impl->m_has_cached_factor
                          || (N != m_last_N);

    if (do_refactor) {
        // --- 1. Upload alpha, gam to device ---
        m_impl->d_alpha.upload(alpha_corrected, N);
        m_impl->d_gam.upload(gam_corrected, N);

        // --- 2. Build N×N Coulomb matrix on GPU ---
        {
            int n_lower = N * (N + 1) / 2;
            int block = 256;
            int grid = (n_lower + block - 1) / block;
            k_eeq_build_matrix<<<grid, block, 0, m_impl->stream>>>(
                N, cx, cy, cz, m_impl->d_alpha.ptr, m_impl->d_gam.ptr, m_impl->d_A.ptr, cutoff_sq);
        }

        // --- 3. cuSOLVER workspace query (cached for same N) ---
        if (N != m_last_N) {
            int lwork = 0;
            checkCusolverEEQ(
                cusolverDnDpotrf_bufferSize(m_impl->cusolver_handle,
                                             CUBLAS_FILL_MODE_LOWER,
                                             N, m_impl->d_A.ptr, N, &lwork),
                "potrf_bufferSize");
            m_impl->workspace_size = lwork;
            if (m_impl->d_workspace.n < lwork)
                m_impl->d_workspace.alloc(lwork);
            m_last_N = N;
        }

        // --- 4. Cholesky factorization: A = L·L^T (overwrites d_A with L) ---
        checkCusolverEEQ(
            cusolverDnDpotrf(m_impl->cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, m_impl->d_A.ptr, N,
                              m_impl->d_workspace.ptr, m_impl->workspace_size,
                              m_impl->d_info.ptr),
            "cusolverDnDpotrf");

        // --- 5. Check factorization success ---
        checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                      cudaMemcpyDeviceToHost, m_impl->stream),
                     "download d_info");
        checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrf");

        if (*m_impl->h_info != 0) {
            // Cholesky failed (matrix not SPD) → caller falls back to CPU
            m_impl->m_has_cached_factor = false;
            return false;
        }

        m_impl->m_has_cached_factor = true;
    } else {
        // Non-refactor step: reuse cached L from d_A, skip matrix build + dpotrf
    }

    // --- 6. Build RHS matrix on CPU and upload (always fresh — captures CN-dependent χ) ---
    // Column-major layout: [b_atoms | C^T_1 | ... | C^T_nfrag], each column N doubles.
    {
        const int rhs_size = N * nrhs;
        if (m_impl->d_rhs.n < rhs_size)
            m_impl->d_rhs.alloc(rhs_size);

        // Reuse cached buffer to avoid per-step heap allocation
        m_impl->m_h_rhs_cache.assign(rhs_size, 0.0);
        double* h_rhs = m_impl->m_h_rhs_cache.data();

        // Column 0: atom electronegativity RHS
        std::memcpy(h_rhs, rhs_atoms, N * sizeof(double));

        // Columns 1..nfrag: C^T columns (C(f,j)=1 if atom j in fragment f)
        for (int f = 0; f < nfrag; ++f) {
            double* col = h_rhs + (f + 1) * N;
            if (fraglist.empty()) {
                if (f == 0) {
                    for (int j = 0; j < N; ++j)
                        col[j] = 1.0;
                }
            } else {
                for (int j = 0; j < N; ++j) {
                    if (fraglist[j] == f + 1)
                        col[j] = 1.0;
                }
            }
        }

        m_impl->d_rhs.upload(h_rhs, rhs_size);
    }

    // --- 7. Triangular solve: L·L^T · X = B (in-place, overwrites d_rhs) ---
    checkCusolverEEQ(
        cusolverDnDpotrs(m_impl->cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, nrhs,
                          m_impl->d_A.ptr, N,
                          m_impl->d_rhs.ptr, N,
                          m_impl->d_info.ptr),
        "cusolverDnDpotrs");

    // --- 8. Download solution columns (z1 and Z2) ---
    const int result_size = N * nrhs;
    checkCudaEEQ(cudaMemcpyAsync(m_impl->h_result, m_impl->d_rhs.ptr,
                                  result_size * sizeof(double),
                                  cudaMemcpyDeviceToHost, m_impl->stream),
                 "download solution");
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrs");

    // --- 9. Copy to caller buffers ---
    std::memcpy(out_z1, m_impl->h_result, N * sizeof(double));
    if (nfrag > 0) {
        std::memcpy(out_Z2, m_impl->h_result + N, N * nfrag * sizeof(double));
    }

    return true;
}


// ============================================================================
// solveWithDeviceRHS (WP2) — alpha/gam/rhs already on GPU, skip H2D uploads
// Claude Generated (May 2026): Eliminates 3×N double H2D per MD step.
// ============================================================================

bool EEQSolverGPU::solveWithDeviceRHS(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha_corrected,
    const double* d_gam_corrected,
    const double* d_rhs_atoms,
    const double* /* d_rhs_constraints — reserved for WP4 GPU Schur */,
    const std::vector<int>& fraglist,
    double* out_z1,
    double* out_Z2,
    double cutoff_sq,
    bool force_refactor)
{
    const int N = natoms;
    const int nrhs = nfrag + 1;

    const bool do_refactor = force_refactor
                          || !m_impl->m_has_cached_factor
                          || (N != m_last_N);

    if (do_refactor) {
        // alpha/gam already on GPU — pass device pointers directly to matrix builder
        {
            int n_lower = N * (N + 1) / 2;
            int block = 256;
            int grid = (n_lower + block - 1) / block;
            k_eeq_build_matrix<<<grid, block, 0, m_impl->stream>>>(
                N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, m_impl->d_A.ptr, cutoff_sq);
        }

        if (N != m_last_N) {
            int lwork = 0;
            checkCusolverEEQ(
                cusolverDnDpotrf_bufferSize(m_impl->cusolver_handle,
                                             CUBLAS_FILL_MODE_LOWER,
                                             N, m_impl->d_A.ptr, N, &lwork),
                "potrf_bufferSize");
            m_impl->workspace_size = lwork;
            if (m_impl->d_workspace.n < lwork)
                m_impl->d_workspace.alloc(lwork);
            m_last_N = N;
        }

        checkCusolverEEQ(
            cusolverDnDpotrf(m_impl->cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, m_impl->d_A.ptr, N,
                              m_impl->d_workspace.ptr, m_impl->workspace_size,
                              m_impl->d_info.ptr),
            "cusolverDnDpotrf");

        checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                      cudaMemcpyDeviceToHost, m_impl->stream),
                     "download d_info");
        checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrf");

        if (*m_impl->h_info != 0) {
            m_impl->m_has_cached_factor = false;
            return false;
        }
        m_impl->m_has_cached_factor = true;
    }

    // Build RHS matrix:
    //   Column 0: D2D copy from d_rhs_atoms (avoids H2D for CN-dependent χ term)
    //   Columns 1..nfrag: topology-constant constraint columns (build on CPU, upload H2D)
    {
        const int rhs_size = N * nrhs;
        if (m_impl->d_rhs.n < rhs_size)
            m_impl->d_rhs.alloc(rhs_size);

        // Column 0: device-to-device (d_rhs_atoms is ready after main stream sync)
        checkCudaEEQ(cudaMemcpy(m_impl->d_rhs.ptr, d_rhs_atoms, N * sizeof(double),
                                 cudaMemcpyDeviceToDevice),
                     "D2D rhs_atoms col0");

        // Columns 1..nfrag: binary fragment indicator (topology-constant)
        if (nfrag > 0) {
            m_impl->m_h_rhs_cache.assign(N * nfrag, 0.0);
            for (int f = 0; f < nfrag; ++f) {
                double* col = m_impl->m_h_rhs_cache.data() + f * N;
                if (fraglist.empty()) {
                    if (f == 0)
                        for (int j = 0; j < N; ++j) col[j] = 1.0;
                } else {
                    for (int j = 0; j < N; ++j)
                        if (fraglist[j] == f + 1) col[j] = 1.0;
                }
            }
            checkCudaEEQ(cudaMemcpy(m_impl->d_rhs.ptr + N, m_impl->m_h_rhs_cache.data(),
                                     N * nfrag * sizeof(double), cudaMemcpyHostToDevice),
                          "H2D constraint cols");
        }
    }

    checkCusolverEEQ(
        cusolverDnDpotrs(m_impl->cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, nrhs,
                          m_impl->d_A.ptr, N,
                          m_impl->d_rhs.ptr, N,
                          m_impl->d_info.ptr),
        "cusolverDnDpotrs");

    const int result_size = N * nrhs;
    checkCudaEEQ(cudaMemcpyAsync(m_impl->h_result, m_impl->d_rhs.ptr,
                                  result_size * sizeof(double),
                                  cudaMemcpyDeviceToHost, m_impl->stream),
                 "download solution");
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrs");

    std::memcpy(out_z1, m_impl->h_result, N * sizeof(double));
    if (nfrag > 0)
        std::memcpy(out_Z2, m_impl->h_result + N, N * nfrag * sizeof(double));

    return true;
}

// ============================================================================
// solveAndComputeCharges — GPU Schur complement for nfrag == 1
// ============================================================================

bool EEQSolverGPU::solveAndComputeCharges(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* alpha_corrected,
    const double* gam_corrected,
    const std::vector<int>& fraglist,
    const double* rhs_atoms,
    const double* rhs_constraints,
    double* out_charges,
    double cutoff_sq)
{
    const int N = natoms;
    const int nrhs = nfrag + 1;

    // --- 1. Upload alpha, gam to device ---
    m_impl->d_alpha.upload(alpha_corrected, N);
    m_impl->d_gam.upload(gam_corrected, N);

    // --- 2. Build N×N Coulomb matrix on GPU ---
    {
        int n_lower = N * (N + 1) / 2;
        int block = 256;
        int grid = (n_lower + block - 1) / block;

        k_eeq_build_matrix<<<grid, block, 0, m_impl->stream>>>(
            N, cx, cy, cz, m_impl->d_alpha.ptr, m_impl->d_gam.ptr, m_impl->d_A.ptr, cutoff_sq);
    }

    // --- 3. Build RHS matrix on CPU and upload ---
    {
        const int rhs_size = N * nrhs;
        if (m_impl->d_rhs.n < rhs_size)
            m_impl->d_rhs.alloc(rhs_size);

        std::vector<double> h_rhs(rhs_size, 0.0);
        std::memcpy(h_rhs.data(), rhs_atoms, N * sizeof(double));
        for (int f = 0; f < nfrag; ++f) {
            double* col = h_rhs.data() + (f + 1) * N;
            if (fraglist.empty()) {
                if (f == 0) {
                    for (int j = 0; j < N; ++j)
                        col[j] = 1.0;
                }
            } else {
                for (int j = 0; j < N; ++j) {
                    if (fraglist[j] == f + 1)
                        col[j] = 1.0;
                }
            }
        }
        m_impl->d_rhs.upload(h_rhs.data(), rhs_size);
    }

    // --- 4. cuSOLVER workspace query ---
    if (N != m_last_N) {
        int lwork = 0;
        checkCusolverEEQ(
            cusolverDnDpotrf_bufferSize(m_impl->cusolver_handle,
                                         CUBLAS_FILL_MODE_LOWER,
                                         N, m_impl->d_A.ptr, N, &lwork),
            "potrf_bufferSize");
        m_impl->workspace_size = lwork;
        if (m_impl->d_workspace.n < lwork)
            m_impl->d_workspace.alloc(lwork);
        m_last_N = N;
    }

    // --- 5. Cholesky factorization ---
    checkCusolverEEQ(
        cusolverDnDpotrf(m_impl->cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, m_impl->d_A.ptr, N,
                          m_impl->d_workspace.ptr, m_impl->workspace_size,
                          m_impl->d_info.ptr),
        "cusolverDnDpotrf");

    checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                  cudaMemcpyDeviceToHost, m_impl->stream),
                 "download d_info");
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrf");

    if (*m_impl->h_info != 0) {
        return false;  // not SPD → caller falls back to CPU
    }

    // --- 6. Triangular solve (in-place in d_rhs) ---
    checkCusolverEEQ(
        cusolverDnDpotrs(m_impl->cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, nrhs,
                          m_impl->d_A.ptr, N,
                          m_impl->d_rhs.ptr, N,
                          m_impl->d_info.ptr),
        "cusolverDnDpotrs");

    // --- 7. Schur complement on GPU (nfrag == 1 fast path) ---
    if (nfrag == 1) {
        // Allocate/zero d_sums [2]
        if (m_impl->d_sums.n < 2)
            m_impl->d_sums.alloc(2);
        cudaMemsetAsync(m_impl->d_sums.ptr, 0, 2 * sizeof(double), m_impl->stream);

        // Reduce Cz1 = sum(z1) and S = sum(Z2_col0)
        int reduce_block = 256;
        int reduce_grid = (N + reduce_block - 1) / reduce_block;
        k_eeq_reduce_sums<<<reduce_grid, reduce_block,
                             2 * reduce_block * sizeof(double), m_impl->stream>>>(
            N, m_impl->d_rhs.ptr, m_impl->d_sums.ptr);

        // Download 2 scalars
        double h_sums[2] = {0.0, 0.0};
        checkCudaEEQ(cudaMemcpyAsync(h_sums, m_impl->d_sums.ptr, 2 * sizeof(double),
                                      cudaMemcpyDeviceToHost, m_impl->stream),
                     "download sums");
        checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after reduce");

        double Cz1 = h_sums[0];
        double S   = h_sums[1];
        double lambda = (Cz1 - rhs_constraints[0]) / S;

        // Compute charges on GPU
        int block = 256;
        int grid = (N + block - 1) / block;
        k_eeq_schur_nfrag1<<<grid, block, 0, m_impl->stream>>>(
            N, m_impl->d_rhs.ptr, lambda, m_impl->d_rhs.ptr);  // reuse d_rhs as output

        // Download charges only (N doubles)
        checkCudaEEQ(cudaMemcpyAsync(out_charges, m_impl->d_rhs.ptr,
                                      N * sizeof(double),
                                      cudaMemcpyDeviceToHost, m_impl->stream),
                     "download charges");
        checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after schur");
    } else {
        // nfrag > 1: not handled by this fast path — caller should use solve() + CPU Schur
        return false;
    }

    return true;
}

// ============================================================================
// solveWithDeviceRHSAndGPUSchur (WP5-A) — WP2 solve + GPU Schur complement
// Claude Generated (May 2026): Eliminates 22 KB D2H (z1+Z2) + O(N) CPU loop
// + 11 KB H2D (charges). Only the Schur scalar lambda crosses CPU boundary.
// ============================================================================

bool EEQSolverGPU::solveWithDeviceRHSAndGPUSchur(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha_corrected,
    const double* d_gam_corrected,
    const double* d_rhs_atoms,
    const std::vector<int>& fraglist,
    double rhs_c0,
    double cutoff_sq,
    bool force_refactor)
{
    // Only nfrag == 1 supported; multi-fragment needs batched Cholesky (future work)
    if (nfrag != 1) return false;

    const int N = natoms;
    const int nrhs = nfrag + 1;  // col0: b_atoms, col1: fragment indicator

    const bool do_refactor = force_refactor
                          || !m_impl->m_has_cached_factor
                          || (N != m_last_N);

    if (do_refactor) {
        // Build N×N Coulomb matrix on GPU (alpha/gam already on device via WP2)
        {
            int n_lower = N * (N + 1) / 2;
            int block = 256;
            int grid = (n_lower + block - 1) / block;
            k_eeq_build_matrix<<<grid, block, 0, m_impl->stream>>>(
                N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, m_impl->d_A.ptr, cutoff_sq);
        }

        if (N != m_last_N) {
            int lwork = 0;
            checkCusolverEEQ(
                cusolverDnDpotrf_bufferSize(m_impl->cusolver_handle,
                                             CUBLAS_FILL_MODE_LOWER,
                                             N, m_impl->d_A.ptr, N, &lwork),
                "potrf_bufferSize");
            m_impl->workspace_size = lwork;
            if (m_impl->d_workspace.n < lwork)
                m_impl->d_workspace.alloc(lwork);
            m_last_N = N;
        }

        checkCusolverEEQ(
            cusolverDnDpotrf(m_impl->cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, m_impl->d_A.ptr, N,
                              m_impl->d_workspace.ptr, m_impl->workspace_size,
                              m_impl->d_info.ptr),
            "cusolverDnDpotrf");

        checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                      cudaMemcpyDeviceToHost, m_impl->stream),
                     "download d_info");
        checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrf");

        if (*m_impl->h_info != 0) {
            m_impl->m_has_cached_factor = false;
            return false;
        }
        m_impl->m_has_cached_factor = true;
    }

    // Build RHS: col0 = d_rhs_atoms (D2D), col1 = fragment indicator (H2D, topology-constant)
    {
        const int rhs_size = N * nrhs;
        if (m_impl->d_rhs.n < rhs_size)
            m_impl->d_rhs.alloc(rhs_size);

        checkCudaEEQ(cudaMemcpy(m_impl->d_rhs.ptr, d_rhs_atoms, N * sizeof(double),
                                 cudaMemcpyDeviceToDevice),
                     "D2D rhs_atoms col0");

        // Column 1: binary fragment indicator (topology-constant for nfrag=1)
        if (m_impl->m_h_rhs_cache.size() < (size_t)N)
            m_impl->m_h_rhs_cache.assign(N, 0.0);
        else
            std::fill(m_impl->m_h_rhs_cache.begin(), m_impl->m_h_rhs_cache.begin() + N, 0.0);

        if (fraglist.empty()) {
            std::fill(m_impl->m_h_rhs_cache.begin(), m_impl->m_h_rhs_cache.begin() + N, 1.0);
        } else {
            for (int j = 0; j < N; ++j)
                if (fraglist[j] == 1) m_impl->m_h_rhs_cache[j] = 1.0;
        }
        checkCudaEEQ(cudaMemcpy(m_impl->d_rhs.ptr + N, m_impl->m_h_rhs_cache.data(),
                                 N * sizeof(double), cudaMemcpyHostToDevice),
                     "H2D constraint col1");
    }

    // Triangular solve: L·L^T · X = B (in-place, d_rhs[0..N-1]=z1, d_rhs[N..2N-1]=Z2)
    checkCusolverEEQ(
        cusolverDnDpotrs(m_impl->cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, nrhs,
                          m_impl->d_A.ptr, N,
                          m_impl->d_rhs.ptr, N,
                          m_impl->d_info.ptr),
        "cusolverDnDpotrs");

    // GPU Schur: q[i] = z1[i] - Z2[i] * lambda
    //   lambda = (sum(z1) - rhs_c0) / sum(Z2)
    if (m_impl->d_sums.n < 2)
        m_impl->d_sums.alloc(2);
    cudaMemsetAsync(m_impl->d_sums.ptr, 0, 2 * sizeof(double), m_impl->stream);

    int rb = 256, rg = (N + rb - 1) / rb;
    k_eeq_reduce_sums<<<rg, rb, 2 * rb * sizeof(double), m_impl->stream>>>(
        N, m_impl->d_rhs.ptr, m_impl->d_sums.ptr);

    // Download only 2 scalars (16 bytes) for lambda — not the full N-vector
    double h_sums[2] = {0.0, 0.0};
    checkCudaEEQ(cudaMemcpyAsync(h_sums, m_impl->d_sums.ptr, 2 * sizeof(double),
                                  cudaMemcpyDeviceToHost, m_impl->stream),
                 "download Schur sums");
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after reduce");

    double Cz1    = h_sums[0];
    double S      = h_sums[1];
    double lambda = (Cz1 - rhs_c0) / S;

    // Write charges into d_rhs[0..N-1] (overwrites z1 in-place)
    int b = 256, g = (N + b - 1) / b;
    k_eeq_schur_nfrag1<<<g, b, 0, m_impl->stream>>>(
        N, m_impl->d_rhs.ptr, lambda, m_impl->d_rhs.ptr);

    // Sync so caller can safely D2D-copy from getDeviceChargesPtr()
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after schur_nfrag1");

    return true;
}

double* EEQSolverGPU::getDeviceChargesPtr()
{
    return m_impl->d_rhs.ptr;
}
