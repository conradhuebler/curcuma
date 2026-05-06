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
    double* h_result      = nullptr; ///< [N*(nfrag+1)] downloaded solution — resizable
    int     m_h_result_size = 0;     ///< current allocation size in doubles
    int*    h_info        = nullptr; ///< [1] cusolver status

    int workspace_size = 0;

    // Lazy Cholesky state: d_A holds the cached factor L after refactorization.
    // force_refactor flag (passed per solve() call) decides whether to rebuild.
    // Caller (GFNFFGPUComputationalMethod) decides based on geometry RMSD.
    bool m_has_cached_factor = false;

    // Cached host RHS buffer to avoid per-step heap allocation
    std::vector<double> m_h_rhs_cache;

    // ── WP6: Batched per-fragment Cholesky (nfrag > 1) ──────────────────────
    // Fragment layout (topology-constant, uploaded once per topology build):
    //   d_frag_sizes[f]        = N_f
    //   d_frag_offsets_A[f]    = prefix sum of N_f² (start in d_A_blocks)
    //   d_frag_offsets_rhs[f]  = prefix sum of N_f*2 (start in d_rhs_blocks)
    //   d_frag_offsets_pair[f] = prefix sum of N_f*(N_f+1)/2 (thread dispatch)
    //   d_frag_atom_offsets[f] = prefix sum of N_f (atom start per fragment)
    //   d_frag_atom_map[k]     = global atom index for sorted position k
    CudaBuffer<double> d_A_blocks;           ///< [sum N_f²] per-fragment Coulomb blocks
    CudaBuffer<double> d_rhs_blocks;         ///< [sum N_f*2] per-fragment RHS
    CudaBuffer<int>    d_frag_sizes;         ///< [nfrag]
    CudaBuffer<int>    d_frag_offsets_A;     ///< [nfrag]
    CudaBuffer<int>    d_frag_offsets_rhs;   ///< [nfrag]
    CudaBuffer<int>    d_frag_offsets_pair;  ///< [nfrag+1]
    CudaBuffer<int>    d_frag_atom_offsets;  ///< [nfrag+1]
    CudaBuffer<int>    d_frag_atom_map;      ///< [N]
    CudaBuffer<double> d_frag_workspace;     ///< [max_frag_ws] reused per-fragment dpotrf
    CudaBuffer<int>    d_frag_info;          ///< [1] dpotrf info per fragment

    // CPU mirrors (computed at uploadFragmentTopology, read during batched solve)
    std::vector<int>    h_frag_sizes;
    std::vector<int>    h_frag_offsets_A;
    std::vector<int>    h_frag_offsets_rhs;
    std::vector<int>    h_frag_offsets_pair;
    std::vector<int>    h_frag_atom_offsets;
    std::vector<int>    h_frag_atom_map;
    std::vector<double> h_charges_global;    ///< [max_natoms] scatter buffer for Schur output

    int  m_nfrag_batched   = 0;
    int  m_max_frag_N      = 0;
    int  m_frag_ws_size    = 0;
    bool m_frag_topo_valid = false;
    bool m_frag_refactored = false;   ///< cached Cholesky factors valid (lazy solve)

    // ── WP7-A: General Schur for nfrag > 1 (May 2026) ────────────────────────
    // Topology-constant device buffers populated by uploadFragmentTopology():
    //   d_atom_frag[i]            = fraglist[i] - 1 (0-indexed fragment id)
    //   d_rhs_constraint_cols[g·N+i] = (fraglist[i] == g+1) ? 1.0 : 0.0
    // Per-step buffers (sized to nfrag at upload):
    //   d_Cz1_general[f]          = Σ_{i∈frag_f} z1[i]   — reduction output
    //   d_S_general[f*nfrag+g]    = Σ_{i∈frag_f} Z2[i,g] — reduction output
    //   d_lambda_general[f]       = Lagrange multipliers (H2D after CPU Schur)
    CudaBuffer<int>     d_atom_frag;
    CudaBuffer<double>  d_rhs_constraint_cols;
    CudaBuffer<double>  d_Cz1_general;
    CudaBuffer<double>  d_S_general;
    CudaBuffer<double>  d_lambda_general;
    std::vector<double> h_Cz1_general;
    std::vector<double> h_S_general;
    std::vector<double> h_lambda_general;
};

// ============================================================================
// h_result resize helper
// ============================================================================

/// Ensure h_result pinned buffer holds at least `needed` doubles.
/// Reallocates (cudaFreeHost + cudaMallocHost) if current size is insufficient.
static void ensureHResult(EEQSolverGPUImpl& impl, int needed)
{
    if (needed <= impl.m_h_result_size) return;
    if (impl.h_result) {
        if (cudaFreeHost(impl.h_result) != cudaSuccess)
            throw std::runtime_error("ensureHResult: cudaFreeHost failed");
        impl.h_result = nullptr;
    }
    if (cudaMallocHost(&impl.h_result, needed * sizeof(double)) != cudaSuccess)
        throw std::runtime_error("ensureHResult: cudaMallocHost failed for "
                                 + std::to_string(needed) + " doubles");
    impl.m_h_result_size = needed;
}

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

    // Pre-allocate d_rhs and h_result for up to MAX_NFRAG fragments.
    // ensureHResult() grows h_result when nfrag > MAX_NFRAG.
    constexpr int MAX_NFRAG = 8;
    m_impl->d_rhs.alloc(max_natoms * (MAX_NFRAG + 1));

    const int initial_h_result = max_natoms * (MAX_NFRAG + 1);
    checkCudaEEQ(cudaMallocHost(&m_impl->h_result, initial_h_result * sizeof(double)),
                 "pinned h_result");
    m_impl->m_h_result_size = initial_h_result;
    checkCudaEEQ(cudaMallocHost(&m_impl->h_info, sizeof(int)), "pinned h_info");

    // WP6: pre-size CPU scatter buffer for batched Schur (host vector, not pinned)
    m_impl->h_charges_global.resize(max_natoms, 0.0);
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
    ensureHResult(*m_impl, result_size);
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
    ensureHResult(*m_impl, result_size);
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

// ============================================================================
// WP7-A: solveWithDeviceRHSAndGPUSchurGeneral — full N×N Cholesky + GPU Schur
//        for nfrag > 1. Replaces D2H of z1+Z2 (~N·(1+nfrag) doubles) with
//        D2H of Cz1+S (nfrag + nfrag² doubles).
// Claude Generated (May 2026)
// Reference: docs/GPU_WP7_EEQ_LARGE_SYSTEMS.md (Strategie A)
// ============================================================================

/// Solve a small dense linear system S·λ = rhs (in-place; nfrag×nfrag) on CPU
/// via Gauss-elim with partial pivoting. Identical to the routine inlined in
/// gfnff_gpu_method.cpp:682-711 (kept duplicate to avoid cross-TU helper for now).
static void solveSchurSystemInPlace(double* S, double* rhs, double* lambda, int nfrag)
{
    for (int col = 0; col < nfrag; ++col) {
        int    pivot = col;
        double maxv  = std::abs(S[col * nfrag + col]);
        for (int row = col + 1; row < nfrag; ++row) {
            double v = std::abs(S[row * nfrag + col]);
            if (v > maxv) { maxv = v; pivot = row; }
        }
        if (pivot != col) {
            for (int k = 0; k < nfrag; ++k)
                std::swap(S[col * nfrag + k], S[pivot * nfrag + k]);
            std::swap(rhs[col], rhs[pivot]);
        }
        double inv_d = (std::abs(S[col * nfrag + col]) > 1e-15)
                           ? 1.0 / S[col * nfrag + col] : 0.0;
        for (int row = col + 1; row < nfrag; ++row) {
            double fac = S[row * nfrag + col] * inv_d;
            for (int k = col + 1; k < nfrag; ++k)
                S[row * nfrag + k] -= fac * S[col * nfrag + k];
            rhs[row] -= fac * rhs[col];
            S[row * nfrag + col] = 0.0;
        }
    }
    for (int i = nfrag - 1; i >= 0; --i) {
        lambda[i] = rhs[i];
        for (int j = i + 1; j < nfrag; ++j)
            lambda[i] -= S[i * nfrag + j] * lambda[j];
        double diag = S[i * nfrag + i];
        lambda[i] = (std::abs(diag) > 1e-15) ? lambda[i] / diag : 0.0;
    }
}

bool EEQSolverGPU::solveWithDeviceRHSAndGPUSchurGeneral(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha_corrected,
    const double* d_gam_corrected,
    const double* d_rhs_atoms,
    const std::vector<int>& /*fraglist*/,
    const std::vector<double>& rhs_constraints,
    double cutoff_sq,
    bool force_refactor)
{
    auto& impl = *m_impl;
    if (nfrag < 1) return false;
    if (!impl.m_frag_topo_valid) return false;
    if (impl.m_nfrag_batched != nfrag) return false;
    if (impl.d_atom_frag.n < natoms) return false;
    if (impl.d_rhs_constraint_cols.n < natoms * nfrag) return false;

    const int N    = natoms;
    const int nrhs = nfrag + 1;

    const bool do_refactor = force_refactor
                          || !impl.m_has_cached_factor
                          || (N != m_last_N);

    if (do_refactor) {
        // Build full N×N Coulomb matrix (lower-triangle threads) — same as WP5-A.
        {
            int n_lower = N * (N + 1) / 2;
            int block   = 256;
            int grid    = (n_lower + block - 1) / block;
            k_eeq_build_matrix<<<grid, block, 0, impl.stream>>>(
                N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, impl.d_A.ptr, cutoff_sq);
        }
        if (N != m_last_N) {
            int lwork = 0;
            checkCusolverEEQ(
                cusolverDnDpotrf_bufferSize(impl.cusolver_handle,
                                             CUBLAS_FILL_MODE_LOWER,
                                             N, impl.d_A.ptr, N, &lwork),
                "potrf_bufferSize (general)");
            impl.workspace_size = lwork;
            if (impl.d_workspace.n < lwork)
                impl.d_workspace.alloc(lwork);
            m_last_N = N;
        }
        checkCusolverEEQ(
            cusolverDnDpotrf(impl.cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, impl.d_A.ptr, N,
                              impl.d_workspace.ptr, impl.workspace_size,
                              impl.d_info.ptr),
            "cusolverDnDpotrf (general)");
        checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_info.ptr, sizeof(int),
                                      cudaMemcpyDeviceToHost, impl.stream),
                     "download d_info (general)");
        checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after potrf (general)");

        if (*impl.h_info != 0) {
            impl.m_has_cached_factor = false;
            return false;
        }
        impl.m_has_cached_factor = true;
    }

    // Assemble (nfrag+1)-column RHS on device:
    //   col 0       = d_rhs_atoms (D2D)
    //   cols 1..nfrag = topology-constant indicator columns (D2D from cache)
    {
        const int rhs_size = N * nrhs;
        if (impl.d_rhs.n < rhs_size)
            impl.d_rhs.alloc(rhs_size);
        checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, d_rhs_atoms,
                                      N * sizeof(double),
                                      cudaMemcpyDeviceToDevice, impl.stream),
                     "D2D rhs_atoms col0 (general)");
        checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr + N,
                                      impl.d_rhs_constraint_cols.ptr,
                                      (size_t)N * nfrag * sizeof(double),
                                      cudaMemcpyDeviceToDevice, impl.stream),
                     "D2D rhs_constraint_cols (general)");
    }

    // dpotrs: in-place solve A·X = [b_atoms | C^T]. After this:
    //   d_rhs[0..N-1]                   = z1
    //   d_rhs[(g+1)·N .. (g+2)·N-1]    = Z2[:,g]   for g = 0..nfrag-1
    checkCusolverEEQ(
        cusolverDnDpotrs(impl.cusolver_handle,
                          CUBLAS_FILL_MODE_LOWER,
                          N, nrhs,
                          impl.d_A.ptr, N,
                          impl.d_rhs.ptr, N,
                          impl.d_info.ptr),
        "cusolverDnDpotrs (general)");

    // GPU reduction: Cz1[f] = Σ_{i∈frag_f} z1[i], S[f,g] = Σ_{i∈frag_f} Z2[i,g].
    cudaMemsetAsync(impl.d_Cz1_general.ptr, 0, nfrag * sizeof(double), impl.stream);
    cudaMemsetAsync(impl.d_S_general.ptr,   0, (size_t)nfrag * nfrag * sizeof(double), impl.stream);
    {
        int rb = 256;
        int rg = (N + rb - 1) / rb;
        k_eeq_reduce_fragment_sums<<<rg, rb, 0, impl.stream>>>(
            N, nfrag,
            impl.d_rhs.ptr,
            impl.d_atom_frag.ptr,
            impl.d_Cz1_general.ptr,
            impl.d_S_general.ptr);
    }

    // D2H the small reduction outputs.
    checkCudaEEQ(cudaMemcpyAsync(impl.h_Cz1_general.data(), impl.d_Cz1_general.ptr,
                                  nfrag * sizeof(double),
                                  cudaMemcpyDeviceToHost, impl.stream),
                 "D2H Cz1 (general)");
    checkCudaEEQ(cudaMemcpyAsync(impl.h_S_general.data(), impl.d_S_general.ptr,
                                  (size_t)nfrag * nfrag * sizeof(double),
                                  cudaMemcpyDeviceToHost, impl.stream),
                 "D2H S (general)");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after reduce (general)");

    // CPU Schur solve: S·λ = Cz1 - Q_frag (Gauss-elim with partial pivoting).
    std::vector<double> schur_rhs(nfrag, 0.0);
    for (int f = 0; f < nfrag; ++f) {
        double q_target = (f < (int)rhs_constraints.size()) ? rhs_constraints[f] : 0.0;
        schur_rhs[f] = impl.h_Cz1_general[f] - q_target;
    }
    solveSchurSystemInPlace(impl.h_S_general.data(),
                             schur_rhs.data(),
                             impl.h_lambda_general.data(),
                             nfrag);

    // H2D the tiny λ vector.
    checkCudaEEQ(cudaMemcpyAsync(impl.d_lambda_general.ptr,
                                  impl.h_lambda_general.data(),
                                  nfrag * sizeof(double),
                                  cudaMemcpyHostToDevice, impl.stream),
                 "H2D lambda (general)");

    // Apply: q[i] = z1[i] − Σ_g Z2[i,g]·λ[g]; write in-place into d_rhs[0..N-1].
    {
        int b = 256;
        int g = (N + b - 1) / b;
        k_eeq_schur_general<<<g, b, 0, impl.stream>>>(
            N, nfrag,
            impl.d_rhs.ptr,
            impl.d_lambda_general.ptr,
            impl.d_rhs.ptr);
    }

    // Sync so caller can D2D-copy from getDeviceChargesPtr() immediately.
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after schur_general");
    return true;
}

// ============================================================================
// WP6: uploadFragmentTopology — one-time setup for batched Cholesky (nfrag > 1)
// Claude Generated (May 2026)
// ============================================================================

void EEQSolverGPU::uploadFragmentTopology(int nfrag,
                                           const std::vector<int>& fraglist,
                                           int natoms)
{
    const int N = natoms;
    auto& impl = *m_impl;

    impl.m_frag_topo_valid  = false;
    impl.m_frag_refactored  = false;
    impl.m_nfrag_batched    = nfrag;

    // --- 1. Count atoms per fragment ---
    impl.h_frag_sizes.assign(nfrag, 0);
    for (int i = 0; i < N; ++i) {
        int f = fraglist[i] - 1;  // 1-indexed → 0-indexed
        if (f >= 0 && f < nfrag)
            impl.h_frag_sizes[f]++;
    }

    // --- 2. Build sorted atom map: atoms sorted by fragment ID (stable) ---
    impl.h_frag_atom_map.resize(N);
    {
        std::vector<int> pos(nfrag, 0);
        // Compute starting positions for each fragment
        impl.h_frag_atom_offsets.resize(nfrag + 1, 0);
        for (int f = 0; f < nfrag; ++f)
            impl.h_frag_atom_offsets[f + 1] = impl.h_frag_atom_offsets[f] + impl.h_frag_sizes[f];
        // Scatter atoms into sorted positions
        for (int i = 0; i < N; ++i) {
            int f = fraglist[i] - 1;
            int slot = impl.h_frag_atom_offsets[f] + pos[f];
            impl.h_frag_atom_map[slot] = i;
            pos[f]++;
        }
    }

    // --- 3. Compute prefix sums ---
    int max_frag_N = 0;
    int total_A    = 0;
    int total_rhs  = 0;
    int total_pair = 0;
    impl.h_frag_offsets_A.resize(nfrag);
    impl.h_frag_offsets_rhs.resize(nfrag);
    impl.h_frag_offsets_pair.resize(nfrag + 1);
    impl.h_frag_offsets_pair[0] = 0;

    for (int f = 0; f < nfrag; ++f) {
        int Nf = impl.h_frag_sizes[f];
        impl.h_frag_offsets_A[f]       = total_A;
        impl.h_frag_offsets_rhs[f]     = total_rhs;
        impl.h_frag_offsets_pair[f + 1] = total_pair + Nf * (Nf + 1) / 2;
        total_A    += Nf * Nf;
        total_rhs  += Nf * 2;
        total_pair += Nf * (Nf + 1) / 2;
        if (Nf > max_frag_N) max_frag_N = Nf;
    }
    impl.m_max_frag_N = max_frag_N;

    // --- 4. Allocate / upload device buffers ---
    if (total_A > 0)   impl.d_A_blocks.alloc(total_A);
    if (total_rhs > 0) impl.d_rhs_blocks.alloc(total_rhs);
    impl.d_frag_sizes.alloc(nfrag);
    impl.d_frag_offsets_A.alloc(nfrag);
    impl.d_frag_offsets_rhs.alloc(nfrag);
    impl.d_frag_offsets_pair.alloc(nfrag + 1);
    impl.d_frag_atom_offsets.alloc(nfrag + 1);
    impl.d_frag_atom_map.alloc(N);
    impl.d_frag_info.alloc(1);

    impl.d_frag_sizes.upload(impl.h_frag_sizes.data(), nfrag);
    impl.d_frag_offsets_A.upload(impl.h_frag_offsets_A.data(), nfrag);
    impl.d_frag_offsets_rhs.upload(impl.h_frag_offsets_rhs.data(), nfrag);
    impl.d_frag_offsets_pair.upload(impl.h_frag_offsets_pair.data(), nfrag + 1);
    impl.d_frag_atom_offsets.upload(impl.h_frag_atom_offsets.data(), nfrag + 1);
    impl.d_frag_atom_map.upload(impl.h_frag_atom_map.data(), N);

    // --- 5. Pre-fill col1 (constraint = all-ones per fragment) in d_rhs_blocks ---
    // Layout: each fragment block is [col0: N_f doubles | col1: N_f doubles]
    // col1 is the binary indicator (1.0 for every atom in the fragment).
    {
        std::vector<double> h_rhs_init(total_rhs, 0.0);
        for (int f = 0; f < nfrag; ++f) {
            int Nf     = impl.h_frag_sizes[f];
            int offset = impl.h_frag_offsets_rhs[f];
            // col1 starts at offset + Nf
            for (int k = 0; k < Nf; ++k)
                h_rhs_init[offset + Nf + k] = 1.0;
        }
        checkCudaEEQ(cudaMemcpy(impl.d_rhs_blocks.ptr, h_rhs_init.data(),
                                 total_rhs * sizeof(double), cudaMemcpyHostToDevice),
                     "H2D rhs_blocks init");
    }

    // --- 6. Query cuSOLVER workspace size for the largest fragment ---
    if (max_frag_N > 0) {
        int lwork = 0;
        checkCusolverEEQ(
            cusolverDnDpotrf_bufferSize(impl.cusolver_handle,
                                        CUBLAS_FILL_MODE_LOWER,
                                        max_frag_N,
                                        impl.d_A_blocks.ptr,  // just a valid device ptr
                                        max_frag_N, &lwork),
            "frag potrf_bufferSize");
        impl.m_frag_ws_size = lwork;
        if (lwork > 0)
            impl.d_frag_workspace.alloc(lwork);
    }

    // --- 7. Resize host scatter buffer if needed ---
    if ((int)impl.h_charges_global.size() < N)
        impl.h_charges_global.resize(N, 0.0);

    // --- 8. WP7-A: per-atom fragment id + topology-constant indicator columns ---
    // d_atom_frag[i]               = fraglist[i] - 1 (0-indexed)
    // d_rhs_constraint_cols[g·N+i] = 1.0 if fraglist[i]==g+1, else 0.0
    const int rhs_cols_size = N * nfrag;
    {
        std::vector<int>    h_atom_frag(N, 0);
        std::vector<double> h_constraint_cols(rhs_cols_size, 0.0);
        for (int i = 0; i < N; ++i) {
            int f = fraglist[i] - 1;
            h_atom_frag[i] = f;
            if (f >= 0 && f < nfrag)
                h_constraint_cols[f * N + i] = 1.0;
        }
        impl.d_atom_frag.alloc(N);
        impl.d_rhs_constraint_cols.alloc(rhs_cols_size);
        impl.d_atom_frag.upload(h_atom_frag.data(), N);
        checkCudaEEQ(cudaMemcpy(impl.d_rhs_constraint_cols.ptr,
                                 h_constraint_cols.data(),
                                 (size_t)rhs_cols_size * sizeof(double),
                                 cudaMemcpyHostToDevice),
                     "H2D d_rhs_constraint_cols");
    }
    impl.d_Cz1_general.alloc(nfrag);
    impl.d_S_general.alloc(nfrag * nfrag);
    impl.d_lambda_general.alloc(nfrag);
    impl.h_Cz1_general.assign(nfrag, 0.0);
    impl.h_S_general.assign((size_t)nfrag * (size_t)nfrag, 0.0);
    impl.h_lambda_general.assign(nfrag, 0.0);

    impl.m_frag_topo_valid = true;
}

bool EEQSolverGPU::isFragmentTopoValid() const
{
    return m_impl->m_frag_topo_valid;
}

// ============================================================================
// WP6: solveWithDeviceRHSAndGPUSchurBatched — batched Cholesky for nfrag > 1
// Claude Generated (May 2026)
// ============================================================================

bool EEQSolverGPU::solveWithDeviceRHSAndGPUSchurBatched(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha,
    const double* d_gam,
    const double* d_rhs_atoms,
    const std::vector<double>& rhs_constraints,
    double cutoff_sq,
    bool force_refactor)
{
    auto& impl = *m_impl;
    if (!impl.m_frag_topo_valid) return false;
    if (nfrag != impl.m_nfrag_batched) return false;

    const int N           = natoms;
    const int total_pairs = impl.h_frag_offsets_pair[nfrag];
    const int total_rhs   = impl.h_frag_offsets_rhs[nfrag - 1]
                          + impl.h_frag_sizes[nfrag - 1] * 2;

    const bool do_refactor = force_refactor || !impl.m_frag_refactored;

    // --- 1. Build per-fragment Coulomb blocks (if refactorizing) ---
    if (do_refactor && total_pairs > 0) {
        int block = 256;
        int grid  = (total_pairs + block - 1) / block;
        k_eeq_build_fragment_matrices<<<grid, block, 0, impl.stream>>>(
            total_pairs, cx, cy, cz, d_alpha, d_gam,
            impl.d_frag_sizes.ptr, impl.d_frag_offsets_A.ptr,
            impl.d_frag_offsets_pair.ptr, impl.d_frag_atom_offsets.ptr,
            impl.d_frag_atom_map.ptr, impl.d_A_blocks.ptr, nfrag, cutoff_sq);
    }

    // --- 2. Gather per-step RHS (col0) from global d_rhs_atoms ---
    {
        int block = 256;
        int grid  = (N + block - 1) / block;
        k_eeq_gather_rhs_fragments<<<grid, block, 0, impl.stream>>>(
            N, d_rhs_atoms,
            impl.d_frag_atom_map.ptr, impl.d_frag_atom_offsets.ptr,
            impl.d_frag_offsets_rhs.ptr, impl.d_frag_sizes.ptr,
            impl.d_rhs_blocks.ptr, nfrag);
    }

    // --- 3. Per-fragment dpotrf + dpotrs (sequential on EEQ stream) ---
    for (int f = 0; f < nfrag; ++f) {
        int    Nf      = impl.h_frag_sizes[f];
        int    A_off   = impl.h_frag_offsets_A[f];
        int    rhs_off = impl.h_frag_offsets_rhs[f];
        double* d_Af   = impl.d_A_blocks.ptr + A_off;
        double* d_Bf   = impl.d_rhs_blocks.ptr + rhs_off;

        if (do_refactor) {
            checkCusolverEEQ(
                cusolverDnDpotrf(impl.cusolver_handle,
                                  CUBLAS_FILL_MODE_LOWER,
                                  Nf, d_Af, Nf,
                                  impl.d_frag_workspace.ptr, impl.m_frag_ws_size,
                                  impl.d_frag_info.ptr),
                "frag dpotrf");

            checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_frag_info.ptr, sizeof(int),
                                          cudaMemcpyDeviceToHost, impl.stream),
                         "frag download d_info");
            checkCudaEEQ(cudaStreamSynchronize(impl.stream), "frag sync after potrf");

            if (*impl.h_info != 0) {
                impl.m_frag_refactored = false;
                return false;  // fragment not SPD → CPU fallback
            }
        }

        // Triangular solve: L·L^T · [z1_f | Z2_f] = [b_f | ones_f]  (2 RHS columns)
        checkCusolverEEQ(
            cusolverDnDpotrs(impl.cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              Nf, 2, d_Af, Nf,
                              d_Bf, Nf,
                              impl.d_frag_info.ptr),
            "frag dpotrs");
    }

    if (do_refactor)
        impl.m_frag_refactored = true;

    // --- 4. Download all per-fragment solutions (2*N doubles total) ---
    ensureHResult(impl, total_rhs);
    checkCudaEEQ(cudaMemcpyAsync(impl.h_result, impl.d_rhs_blocks.ptr,
                                  total_rhs * sizeof(double),
                                  cudaMemcpyDeviceToHost, impl.stream),
                 "frag D2H rhs_blocks");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "frag sync after D2H");

    // --- 5. Host Schur complement per fragment + scatter to global order ---
    // Each fragment: lambda_f = (sum(z1_f) - rhs_constraints[f]) / sum(Z2_f)
    //                q_f[k]  = z1_f[k] - Z2_f[k] * lambda_f
    for (int f = 0; f < nfrag; ++f) {
        int    Nf      = impl.h_frag_sizes[f];
        int    rhs_off = impl.h_frag_offsets_rhs[f];
        double* z1_f   = impl.h_result + rhs_off;
        double* Z2_f   = impl.h_result + rhs_off + Nf;

        double Cz1 = 0.0, S = 0.0;
        for (int k = 0; k < Nf; ++k) {
            Cz1 += z1_f[k];
            S   += Z2_f[k];
        }

        double qfrag_f  = (f < (int)rhs_constraints.size()) ? rhs_constraints[f] : 0.0;
        double lambda_f = (S != 0.0) ? (Cz1 - qfrag_f) / S : 0.0;

        int atom_off = impl.h_frag_atom_offsets[f];
        for (int k = 0; k < Nf; ++k) {
            int gi = impl.h_frag_atom_map[atom_off + k];
            impl.h_charges_global[gi] = z1_f[k] - Z2_f[k] * lambda_f;
        }
    }

    // --- 6. Upload global charges to d_rhs[0..N-1] (getDeviceChargesPtr() target) ---
    if (impl.d_rhs.n < N)
        impl.d_rhs.alloc(N);

    checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, impl.h_charges_global.data(),
                                  N * sizeof(double),
                                  cudaMemcpyHostToDevice, impl.stream),
                 "frag H2D charges");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "frag sync after H2D charges");

    return true;
}
