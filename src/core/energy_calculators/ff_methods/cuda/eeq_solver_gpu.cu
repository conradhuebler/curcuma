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
    CudaBuffer<double> d_A;          ///< [N*N] column-major Coulomb matrix
    CudaBuffer<double> d_rhs;        ///< [N*(nfrag+1)] RHS matrix (column-major)
    CudaBuffer<double> d_workspace;  ///< cuSOLVER workspace
    CudaBuffer<int>    d_info;       ///< [1] cuSOLVER info

    // Pinned host buffers for async D2H transfer
    double* h_result = nullptr;      ///< [N*(nfrag+1)] downloaded solution
    int*    h_info = nullptr;        ///< [1] cusolver status

    int workspace_size = 0;
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
    const double* __restrict__ coords,   // [N*3] row-major AoS: x₀y₀z₀x₁y₁z₁... (Bohr)
    const double* __restrict__ alpha,     // [N] alpha² = (alpha_base + ff*qa)², already squared
    const double* __restrict__ gam,       // [N] gam_corrected
    double* __restrict__ A)               // [N*N] output column-major
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
        // Off-diagonal: compute interatomic distance
        double dx = coords[i * 3 + 0] - coords[j * 3 + 0];
        double dy = coords[i * 3 + 1] - coords[j * 3 + 1];
        double dz = coords[i * 3 + 2] - coords[j * 3 + 2];
        double r = sqrt(dx * dx + dy * dy + dz * dz);

        double gamma_ij = 1.0 / sqrt(alpha[i] + alpha[j]);
        double val = erf(gamma_ij * r) / r;

        // Column-major: element (row, col) at col*N + row
        // Write both triangles (cuSOLVER potrs reads full matrix)
        A[j * N + i] = val;   // lower triangle
        A[i * N + j] = val;   // upper triangle (symmetric)
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
    const double* d_coords,
    const double* alpha_corrected,
    const double* gam_corrected,
    const std::vector<int>& fraglist,
    const double* rhs_atoms,
    const double* rhs_constraints,
    double* out_z1,
    double* out_Z2)
{
    const int N = natoms;
    const int nrhs = nfrag + 1;  // b_atoms + nfrag constraint columns

    // --- 1. Upload alpha, gam to device ---
    m_impl->d_alpha.upload(alpha_corrected, N);
    m_impl->d_gam.upload(gam_corrected, N);

    // --- 2. Build N×N Coulomb matrix on GPU ---
    {
        int n_lower = N * (N + 1) / 2;
        int block = 256;
        int grid = (n_lower + block - 1) / block;

        k_eeq_build_matrix<<<grid, block, 0, m_impl->stream>>>(
            N, d_coords, m_impl->d_alpha.ptr, m_impl->d_gam.ptr, m_impl->d_A.ptr);
    }

    // --- 3. Build RHS matrix on CPU and upload ---
    // Column-major layout: [b_atoms | C^T_1 | ... | C^T_nfrag]
    // Each column is N doubles.
    {
        const int rhs_size = N * nrhs;
        if (m_impl->d_rhs.n < rhs_size)
            m_impl->d_rhs.alloc(rhs_size);

        std::vector<double> h_rhs(rhs_size, 0.0);

        // Column 0: atom electronegativity RHS
        std::memcpy(h_rhs.data(), rhs_atoms, N * sizeof(double));

        // Columns 1..nfrag: C^T columns (constraint matrix transpose)
        // C(f, j) = 1 if atom j belongs to fragment f (1-indexed fraglist)
        for (int f = 0; f < nfrag; ++f) {
            double* col = h_rhs.data() + (f + 1) * N;
            if (fraglist.empty()) {
                // Single fragment: all atoms in fragment 0
                if (f == 0) {
                    for (int j = 0; j < N; ++j)
                        col[j] = 1.0;
                }
            } else {
                for (int j = 0; j < N; ++j) {
                    if (fraglist[j] == f + 1)  // fraglist is 1-indexed
                        col[j] = 1.0;
                }
            }
        }

        m_impl->d_rhs.upload(h_rhs.data(), rhs_size);
    }

    // --- 4. cuSOLVER workspace query (cache for same N) ---
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

    // --- 5. Cholesky factorization: A = L·L^T ---
    checkCusolverEEQ(
        cusolverDnDpotrf(m_impl->cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, m_impl->d_A.ptr, N,
                          m_impl->d_workspace.ptr, m_impl->workspace_size,
                          m_impl->d_info.ptr),
        "cusolverDnDpotrf");

    // --- 6. Check factorization success ---
    checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                  cudaMemcpyDeviceToHost, m_impl->stream),
                 "download d_info");
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after potrf");

    if (*m_impl->h_info != 0) {
        // Cholesky failed (matrix not SPD) → caller falls back to CPU
        return false;
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
    // Column 0 → z1 (A⁻¹ · b_atoms)
    std::memcpy(out_z1, m_impl->h_result, N * sizeof(double));

    // Columns 1..nfrag → Z2 (A⁻¹ · C^T, column-major)
    if (nfrag > 0) {
        std::memcpy(out_Z2, m_impl->h_result + N, N * nfrag * sizeof(double));
    }

    return true;
}
