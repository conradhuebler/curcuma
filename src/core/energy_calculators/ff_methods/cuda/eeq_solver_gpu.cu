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
#include "src/core/curcuma_logger.h"

#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <limits>
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

static void checkCublasEEQ(cublasStatus_t status, const char* msg)
{
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(msg) + ": cublas error " + std::to_string(static_cast<int>(status)));
}

// ============================================================================
// Pimpl implementation struct
// ============================================================================

struct EEQSolverGPUImpl {
    cusolverDnHandle_t cusolver_handle = nullptr;
    cublasHandle_t     cublas_handle   = nullptr;  ///< WP7-C: PCG matvec/dot/axpy
    cudaStream_t stream = nullptr;
    cudaEvent_t  ev_legacy_barrier = nullptr;  ///< F-12: orders impl.stream after the workspace stream
    cudaEvent_t  ev_input_ready    = nullptr;  ///< K3 (Oct 2026): workspace event, not owned; replaces the barrier
    bool         stream_high_priority = false; ///< K2/K3 (Oct 2026): stream created with the greatest priority

    /// Claude Generated (Oct 2026, CUDA/ROCm parity item 10): orders `stream` after the work the FFWorkspaceGPU
    /// streams did for this step (coordinates, corrected alpha/gam, RHS). All curcuma streams are blocking
    /// streams: they are ordered against the legacy default stream but not against each other. A device-RHS
    /// solve that launches its matrix build as the first operation on `stream` therefore needs this wait; the
    /// paths that start with a blocking cudaMemcpy get it implicitly. Measured before this call was added to
    /// the lu and batched paths: water dimer MD with a second GPU job on the same card gave 29 and 22 different
    /// trajectories in 40 runs for `lu` and `batched`, 1 for the others.
    /// `allow_input_event`: use the workspace's input-ready event (K3) instead of the legacy barrier, only valid
    /// for the paths that need nothing the event does not cover (the projected PCG).
    void orderAfterWorkspace(bool allow_input_event)
    {
        if (allow_input_event && ev_input_ready) {
            cudaStreamWaitEvent(stream, ev_input_ready, 0);
            return;
        }
        if (!ev_legacy_barrier)
            cudaEventCreateWithFlags(&ev_legacy_barrier, cudaEventDisableTiming);
        cudaEventRecord(ev_legacy_barrier, cudaStreamLegacy);
        cudaStreamWaitEvent(stream, ev_legacy_barrier, 0);
    }

    // Claude Generated (Sep 2026, F-3): multi-GPU split of the WP7-E projected PCG. Each helper
    // device owns a column block [c0, c1) of the symmetric EEQ matrix; this device owns [0, c1_own).
    struct SplitPart {
        int dev = -1, c0 = 0, c1 = 0;
        cudaStream_t   stream = nullptr;
        cublasHandle_t blas   = nullptr;
        cudaEvent_t    done   = nullptr;
        CudaBuffer<double> A, x, y, z, alpha, gam, p, yblk;
        CudaBuffer<float> Af, xf, yf;   ///< Oct 2026: FP32 block, its input slice and result (mixed-precision matvec)
    };
    int own_device = -1;
    std::vector<int> split_devs;                 ///< helper devices (configured)
    int split_min_atoms = 0;
    std::vector<std::unique_ptr<SplitPart>> split_parts;
    int split_N = 0;
    int split_c1_own = 0;
    CudaBuffer<double> d_A_own;                  ///< [N * c1_own] own column block
    CudaBuffer<float> d_A_own_f, d_xf, d_yf;     ///< Oct 2026: FP32 block of this device, input and result vectors
    double ppcg_true_sq = -1.0;                  ///< Oct 2026: |P(b-Aq)|^2 of the last mixed solve (measured with the FP64 matvec), <0 if unknown
    bool split_mixed_req = true;                 ///< Oct 2026: mixed precision requested (gpu_eeq_mixed_precision)
    double mixed_rho = 1e-5;                     ///< round reduction of the mixed solve (gpu_eeq_mixed_rho)
    bool split_mixed = false;                    ///< Oct 2026: matvec of the split runs in FP32 (CURCUMA_EEQ_MIXED)
    CudaBuffer<double> d_ppcg_scal;              ///< PCG device scalars (Sep 2026, one sync per iteration)
    cudaEvent_t ev_split = nullptr;              ///< recorded on stream, waited on by helpers
    bool split_failed = false;
    std::string split_status;
    long split_matvecs = 0, split_builds = 0;
    void releaseSplit() {
        for (auto& pp : split_parts) {
            if (!pp) continue;
            cudaSetDevice(pp->dev);
            pp->Af.free(); pp->xf.free(); pp->yf.free();
            pp->A.free(); pp->x.free(); pp->y.free(); pp->z.free(); pp->alpha.free();
            pp->gam.free(); pp->p.free(); pp->yblk.free();
            if (pp->blas)   cublasDestroy(pp->blas);
            if (pp->done)   cudaEventDestroy(pp->done);
            if (pp->stream) cudaStreamDestroy(pp->stream);
        }
        split_parts.clear();
        if (own_device >= 0) cudaSetDevice(own_device);
        d_A_own.free(); d_A_own_f.free(); d_xf.free(); d_yf.free();
        split_N = 0;
    }

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

    // ── Step B (May 2026): LU fallback for indefinite EEQ matrices ──────────
    // cusolverDnDpotrf requires SPD; the EEQ Coulomb matrix can be indefinite
    // (the polymer test case has min_gersh ≈ -60.7, κ_est = inf). When dpotrf
    // returns info != 0, rebuild d_A and refactor with cusolverDnDgetrf (LU
    // with partial pivoting), then use dgetrs in place of dpotrs. m_using_lu
    // persists with m_has_cached_factor: lazy solves on subsequent calls reuse
    // the cached LU factor and pick dgetrs accordingly. Mirrors the CPU
    // dispatcher rule (eeq_solver.cpp:1391-1430) that switches PCG → LU when
    // the Gershgorin bound is non-positive.
    CudaBuffer<int>    d_lu_pivots;       ///< [N] partial-pivot indices from dgetrf
    CudaBuffer<double> d_lu_workspace;    ///< cuSOLVER dgetrf workspace
    int  lu_workspace_size = 0;
    bool m_using_lu        = false;       ///< true: cached factor is LU, use dgetrs
    int  m_last_chol_info  = 0;           ///< first failed leading minor (for diagnostics)

    // ── WP-B (Jun 2026): FP32-factor + FP64-refine (mixed-precision dsposv) ──────
    // When m_mixed_precision is on, d_A keeps the FP64 *matrix* (not factored in
    // place); d_A_f32 holds the FP32 Cholesky factor. The per-step solve does an
    // FP32 substitution + FP64 residual refinement → full FP64 accuracy at a
    // fraction of the FP64-factor cost on FP64-weak cards. m_cached_factor_mixed
    // records the cached-factor kind so lazy (non-refactor) solves dispatch right.
    bool   m_mixed_precision      = false;
    int    m_refine_iters         = 2;     ///< max FP64-residual / FP32-correction steps (>=1)
    bool   m_cached_factor_mixed  = false; ///< cached factor lives in d_A_f32 (FP32)
    CudaBuffer<float>  d_A_f32;            ///< [N*N] FP32 Cholesky factor
    CudaBuffer<float>  d_ws_f32;           ///< cuSOLVER spotrf workspace (FP32)
    int    ws_f32_size            = 0;
    CudaBuffer<float>  d_rhs_f32;          ///< [N*nrhs] FP32 solve scratch
    CudaBuffer<double> d_x_ref;            ///< [N*nrhs] FP64 running solution
    CudaBuffer<double> d_b_ref;            ///< [N*nrhs] preserved FP64 RHS
    CudaBuffer<double> d_r_ref;            ///< [N*nrhs] FP64 residual

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
    // Claude Generated (Oct 2026): pristine copy of d_rhs_blocks as uploaded (col0 = 0,
    // col1 = 1). dpotrs overwrites BOTH columns of d_rhs_blocks with the solution, so the
    // all-ones constraint column must be restored before every solve; it used to be set only
    // at topology upload, and from the second call on the solver used Z2 of the previous
    // step as right-hand side (wrong lambda, MD blow-up on polymer_2x at step 2).
    CudaBuffer<double> d_rhs_blocks_init;    ///< [sum N_f*2] col0 = 0, col1 = 1 per fragment
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

    int    m_nfrag_batched      = 0;
    int    m_max_frag_N         = 0;
    int    m_frag_ws_size       = 0;
    bool   m_frag_topo_valid    = false;
    bool   m_frag_refactored    = false;  ///< cached Cholesky factors valid (lazy solve)

    // ── WP7-C: GPU PCG (May 2026) ────────────────────────────────────────────
    // Persistent warm-start buffers (survive across solve calls; invalidated on
    // topology change or large geometry jump):
    CudaBuffer<double>  d_z1_persistent;       ///< [N] last z1 solution
    CudaBuffer<double>  d_Z2_persistent;       ///< [N·nfrag] last Z2 columns
    bool                m_pcg_warm_valid = false;
    bool                m_pcg_M_inv_valid = false;  ///< Jacobi precond cache (re-extract on refactor)
    // WP7-D (Jun 2026): block-Jacobi preconditioner. When valid, the PCG applies
    // z = blockdiag(A_ff^-1)·r (exact per-fragment inverse, stored in d_A_blocks) instead
    // of the diagonal Jacobi — far fewer iterations for many-fragment systems. Falls back
    // to the diagonal (always extracted) when a fragment block is not SPD or too large.
    bool                m_pcg_block_jacobi_valid = false;
    // Per-call scratch (size N — re-allocated when N grows):
    CudaBuffer<double>  d_pcg_M_inv;           ///< [N] 1/A[i,i]
    CudaBuffer<double>  d_pcg_r;               ///< [N] residual
    CudaBuffer<double>  d_pcg_z;               ///< [N] M_inv·r
    CudaBuffer<double>  d_pcg_p;               ///< [N] search direction
    CudaBuffer<double>  d_pcg_Ap;              ///< [N] A·p
    CudaBuffer<double>  d_pcg_dot_scratch;     ///< [2] device-pointer scalars (rz, pAp)
    CudaBuffer<double>  d_pcg_rnorm_scratch;   ///< [1] |r|² for convergence check
    int                 m_pcg_total_iters     = 0;
    int                 m_pcg_total_calls     = 0;
    int                 m_pcg_nonconv_calls   = 0;

    // ── WP7-A: General Schur for nfrag > 1 (May 2026) ────────────────────────
    // Topology-constant device buffers populated by uploadFragmentTopology():
    //   d_atom_frag[i]            = fraglist[i] - 1 (0-indexed fragment id)
    //   d_rhs_constraint_cols[g·N+i] = (fraglist[i] == g+1) ? 1.0 : 0.0
    // Per-step buffers (sized to nfrag at upload):
    //   d_Cz1_general[f]          = Σ_{i∈frag_f} z1[i]   — reduction output
    //   d_S_general[f*nfrag+g]    = Σ_{i∈frag_f} Z2[i,g] — reduction output
    //   d_lambda_general[f]       = Lagrange multipliers (H2D after CPU Schur)
    CudaBuffer<int>     d_atom_frag;
    CudaBuffer<int>     d_frag_off;   // CSR of the atoms of each fragment (ascending), for deterministic fragment sums
    CudaBuffer<int>     d_frag_idx;
    CudaBuffer<double>  d_rhs_constraint_cols;
    CudaBuffer<double>  d_Cz1_general;
    CudaBuffer<double>  d_S_general;
    CudaBuffer<double>  d_lambda_general;
    std::vector<double> h_Cz1_general;
    std::vector<double> h_S_general;
    std::vector<double> h_lambda_general;

    // ── WP7-E: GPU projected PCG (Sep 2026) ──────────────────────────────────
    // Single-solve port of the CPU EEQSolver::solveWithProjectedPCG: instead of
    // (nfrag+1) independent PCG solves + a CPU Schur reduction (WP7-C above), one
    // PCG runs directly in the constraint tangent space (Σ_{i∈frag} v_i = 0),
    // returning the final charges without any Schur step. Existing d_pcg_r/z/p/Ap
    // and d_pcg_M_inv (WP7-C, above) are reused as per-call N-sized scratch; the
    // buffers here are the small (nfrag-sized) additions the projection needs.
    //   d_frag_inv_atoms[f]  = 1 / N_f                      (topology-constant)
    //   d_frag_sM[f]         = Σ_{i∈frag_f} M_inv[i]         (valid with m_pcg_M_inv_valid)
    //   d_frag_sum_scratch   = per-call Σ_{i∈frag_f} v[i] reduction target
    //   d_frag_delta_scratch = per-call per-fragment correction to scatter back
    //   d_rhs_constraints_dev = per-call upload of the (small) host rhs_constraints
    //   d_ppcg_q_persistent  = warm-start charges, independent of WP7-C's z1/Z2
    CudaBuffer<double>  d_frag_inv_atoms;
    CudaBuffer<double>  d_frag_sM;
    CudaBuffer<double>  d_frag_sum_scratch;
    CudaBuffer<double>  d_frag_delta_scratch;
    CudaBuffer<double>  d_rhs_constraints_dev;
    CudaBuffer<double>  d_ppcg_q_persistent;
    bool                m_ppcg_warm_valid = false;
    // Claude Generated (Oct 2026): local-inverse preconditioner of WP7-E (eeq_ppcg_precond
    // localinv), uploaded from the host: CSR of G and H = G C^T, dense S^-1 (nfrag x nfrag).
    CudaBuffer<int>     d_li_G_rowptr, d_li_G_col, d_li_H_rowptr, d_li_H_col;
    CudaBuffer<double>  d_li_G_val, d_li_H_val, d_li_S_inv, d_li_t, d_li_mu, d_li_q0;
    bool                m_li_valid = false;
    int                 m_li_N = 0, m_li_nfrag = 0;
    long                m_li_runtime_fallbacks = 0;
    // Claude Generated (Oct 2026): device build of the local inverse (buildLocalInverseOnDevice).
    // m_li_device = true: G was built here and is applied as z = P G r (no H, no S^-1, see there).
    bool                m_li_device = false;
    bool                m_li_build_pending = false;   ///< requestLocalInverseDeviceBuild() was called
    double              m_li_radius = 8.0, m_li_ext = 10.0;
    CudaBuffer<int>     d_li_nblk, d_li_nkeep, d_li_fail, d_li_atoms;
    CudaBuffer<float4>  d_li_xf;                      ///< float copy of the coordinates (prefilter)
    CudaBuffer<unsigned int> d_li_maxabs;
    CudaBuffer<double>  d_li_Gl_val, d_li_scratch;
    int                 m_li_max_smem = -1;           ///< opt-in shared memory per block (bytes), queried once
    EEQSolverGPU::LocalInverseDeviceStats m_li_stats;
    // Claude Generated (Oct 2026): dense LU of the augmented (N+nfrag) system
    // (solve_method lu, solveWithDeviceRHSAndGPUAugmentedLU).
    CudaBuffer<double>  d_aug_A;       ///< [M*M] augmented matrix, then its LU factor (column-major)
    CudaBuffer<double>  d_aug_rhs;     ///< [M] [b_atoms; Q_frag], then the solution [q; lambda]
    CudaBuffer<int>     d_aug_piv;     ///< [M] getrf pivots
    CudaBuffer<double>  d_aug_ws;      ///< getrf workspace
    int                 m_aug_ws_size = 0;
    // Claude Generated (Oct 2026): Bunch-Kaufman LDL^T of A_nn (solve_method ldlt): cusolverDnDsytrf
    // takes 32-bit pivots, cusolverDnXsytrs 64-bit ones.
    CudaBuffer<int>           d_ldlt_piv;     ///< [N] sytrf pivots (LAPACK convention)
    CudaBuffer<long long>     d_ldlt_piv64;   ///< [N] the same as int64 for Xsytrs
    CudaBuffer<double>        d_ldlt_ws;      ///< sytrf workspace (doubles)
    CudaBuffer<unsigned char> d_ldlt_xs;      ///< Xsytrs device workspace (bytes)
    std::vector<unsigned char> h_ldlt_xs;     ///< Xsytrs host workspace (bytes)
    int                 m_ppcg_total_iters   = 0;
    int                 m_ppcg_total_calls   = 0;
    int                 m_ppcg_nonconv_calls = 0;
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
 *             (= 1 / r_ij exactly when r_ij^2 >= 36 (alpha_i + alpha_j), no erf call there)
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
            const double asum = alpha[i] + alpha[j];
            double val;
            if (r_sq >= GFNFF_COUL_FARFIELD_X2 * asum) {
                // Far field (Claude Generated, Oct 2026): erf(gamma_ij r) == 1.0 exactly for
                // gamma_ij r >= 6 (GFNFF_COUL_FARFIELD_X2), so 1/r is bit-identical to erf()/r.
                val = 1.0 / r;
            } else {
                double gamma_ij = 1.0 / sqrt(asum);
                val = erf(gamma_ij * r) / r;
            }

            // Column-major: element (row, col) at col*N + row
            // Write both triangles (cuSOLVER potrs reads full matrix)
            A[j * N + i] = val;   // lower triangle
            A[i * N + j] = val;   // upper triangle (symmetric)
        }
    }
}

/**
 * @brief Column block [c0, c0+ncols) of the EEQ matrix, column-major (N rows), for the
 * multi-GPU split. Claude Generated (Sep 2026). Same element formula and operand order as
 * k_eeq_build_matrix (i = max(row,col), j = min(row,col)), so the values are bit-identical.
 * Far field uses 1/r like k_eeq_build_matrix (erf saturates to 1.0, GFNFF_COUL_FARFIELD_X2).
 */
__global__ __launch_bounds__(256, 4)
void k_eeq_build_cols(int N, int c0, int ncols,
                      const double* __restrict__ cx, const double* __restrict__ cy,
                      const double* __restrict__ cz, const double* __restrict__ alpha,
                      const double* __restrict__ gam, double* __restrict__ A_blk, double cutoff_sq,
                      float* __restrict__ A_blk32 = nullptr)
{
    // Oct 2026: 2-D grid (x = rows, y = columns, grid-stride in y) instead of one linear index: the 64-bit
    // tid % N and tid / N cost more than the element itself in the far field. Same element, same value.
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= N) return;
    for (int lc = blockIdx.y; lc < ncols; lc += gridDim.y) {
    const long long tid = static_cast<long long>(lc) * N + row;
    const int col = c0 + lc;
    double v;
    if (row == col) {
        v = gam[row] + TSQRT2PI_GPU / sqrt(alpha[row]);
    } else {
        const int i = row > col ? row : col;
        const int j = row > col ? col : row;
        const double dx = cx[i] - cx[j];
        const double dy = cy[i] - cy[j];
        const double dz = cz[i] - cz[j];
        const double r_sq = dx * dx + dy * dy + dz * dz;
        if (cutoff_sq > 0.0 && r_sq > cutoff_sq) {
            v = 0.0;
        } else {
            const double asum = alpha[i] + alpha[j];
            if (r_sq >= GFNFF_COUL_FARFIELD_X2 * asum) {
                // Far field (Claude Generated, Oct 2026): erf == 1.0 as in k_eeq_build_matrix; rsqrt (1 ulp)
                // replaces the FP64 sqrt and division (r_sq >= 36 asum > 0, so no zero distance here)
                v = rsqrt(r_sq);
            } else {
                const double r = sqrt(r_sq);
                const double gamma_ij = 1.0 / sqrt(asum);
                v = erf(gamma_ij * r) / r;
            }
        }
    }
    A_blk[tid] = v;
    if (A_blk32) A_blk32[tid] = static_cast<float>(v);   // Oct 2026: FP32 copy for the mixed-precision matvec
    }
}

/// 1 / A_ii straight from alpha/gam (= k_pcg_extract_diag_inv on the full matrix). Claude Generated.
__global__ void k_eeq_diag_inv(int N, const double* __restrict__ alpha,
                               const double* __restrict__ gam, double* __restrict__ M_inv)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    const double diag = gam[i] + TSQRT2PI_GPU / sqrt(alpha[i]);
    M_inv[i] = (diag != 0.0) ? 1.0 / diag : 0.0;
}

// ============================================================================
// WP-B (Jun 2026): mixed-precision cast/axpy helpers (FP32 factor + FP64 refine)
// ============================================================================

__global__ void k_cast_d2f(const double* __restrict__ in, float* __restrict__ out, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = static_cast<float>(in[i]);
}

__global__ void k_cast_f2d(const float* __restrict__ in, double* __restrict__ out, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = static_cast<double>(in[i]);
}

// x += (double)dx — accumulate the FP32 correction into the FP64 running solution.
__global__ void k_axpy_f2d(double* __restrict__ x, const float* __restrict__ dx, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) x[i] += static_cast<double>(dx[i]);
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
    cudaGetDevice(&m_impl->own_device);
    checkCusolverEEQ(cusolverDnSetStream(m_impl->cusolver_handle, m_impl->stream),
                     "cusolverDnSetStream");

    // WP7-C: cuBLAS handle (PCG matvec/dot/axpy). Same stream as cusolver — no sync needed.
    checkCublasEEQ(cublasCreate(&m_impl->cublas_handle), "cublasCreate");
    checkCublasEEQ(cublasSetStream(m_impl->cublas_handle, m_impl->stream), "cublasSetStream");

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
        if (m_impl->cublas_handle) cublasDestroy(m_impl->cublas_handle);
        if (m_impl->cusolver_handle) cusolverDnDestroy(m_impl->cusolver_handle);
        m_impl->releaseSplit();
        if (m_impl->ev_split) cudaEventDestroy(m_impl->ev_split);
        if (m_impl->ev_legacy_barrier) cudaEventDestroy(m_impl->ev_legacy_barrier);
        if (m_impl->stream) cudaStreamDestroy(m_impl->stream);
    }
}

// ============================================================================
// WP-B (Jun 2026): mixed-precision FP32-factor + FP64-refine solve (dsposv-style)
// ============================================================================

// Claude Generated (Sep 2026, F-3): see eeq_solver_gpu.h.
void EEQSolverGPU::setSplitDevices(const std::vector<int>& helper_devices, int min_atoms)
{
    m_impl->releaseSplit();
    m_impl->split_devs.clear();
    for (int d : helper_devices)
        if (d != m_impl->own_device) m_impl->split_devs.push_back(d);
    m_impl->split_min_atoms = min_atoms;
    m_impl->split_failed = false;
    m_impl->split_status.clear();
}

// Claude Generated (Oct 2026): see eeq_solver_gpu.h.
void EEQSolverGPU::setSplitMixedPrecision(bool on, double rho)
{
    m_impl->releaseSplit();   // the FP32 blocks are allocated with the split
    m_impl->split_mixed_req = on;
    if (rho > 0.0) m_impl->mixed_rho = rho;
}

std::string EEQSolverGPU::splitStatus() const
{
    const auto& I = *m_impl;
    if (I.split_devs.empty()) return {};
    if (!I.split_status.empty()) return I.split_status;
    if (I.split_N == 0)
        return "configured, not used yet (projected-PCG path needs nfrag > 1 and N >= "
               + std::to_string(I.split_min_atoms) + ")";
    std::string list = std::to_string(I.own_device);
    for (int d : I.split_devs) list += "," + std::to_string(d);
    return "EEQ matrix + PCG matvec on devices [" + list + "], " + std::to_string(I.split_builds)
        + " builds, " + std::to_string(I.split_matvecs) + " matvecs";
}

void EEQSolverGPU::setMixedPrecision(bool enabled, int refine_iters)
{
    m_impl->m_mixed_precision = enabled;
    m_impl->m_refine_iters    = (refine_iters < 1) ? 1 : refine_iters;
}

// FP32 Cholesky of the FP64 matrix currently in d_A. d_A is left intact (the
// residual refinement needs the unfactored matrix); the FP32 factor lands in
// d_A_f32. Returns the cuSOLVER info (0 = SPD success). On failure the caller
// falls back to the FP64 dpotrf/LU path, which sees the still-intact d_A.
int EEQSolverGPU::mixedFactor(int N)
{
    const int n2 = N * N;
    if (m_impl->d_A_f32.n < n2) m_impl->d_A_f32.alloc(n2);

    int block = 256;
    int grid  = (n2 + block - 1) / block;
    k_cast_d2f<<<grid, block, 0, m_impl->stream>>>(m_impl->d_A.ptr, m_impl->d_A_f32.ptr, n2);

    int lwork = 0;
    checkCusolverEEQ(
        cusolverDnSpotrf_bufferSize(m_impl->cusolver_handle, CUBLAS_FILL_MODE_LOWER,
                                    N, m_impl->d_A_f32.ptr, N, &lwork),
        "Spotrf_bufferSize");
    if (m_impl->d_ws_f32.n < lwork) m_impl->d_ws_f32.alloc(lwork);
    m_impl->ws_f32_size = lwork;

    checkCusolverEEQ(
        cusolverDnSpotrf(m_impl->cusolver_handle, CUBLAS_FILL_MODE_LOWER,
                         N, m_impl->d_A_f32.ptr, N,
                         m_impl->d_ws_f32.ptr, m_impl->ws_f32_size, m_impl->d_info.ptr),
        "cusolverDnSpotrf");
    checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                 cudaMemcpyDeviceToHost, m_impl->stream),
                 "download Spotrf info");
    checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after Spotrf");
    return *m_impl->h_info;
}

// Solve A·X = B (B in d_rhs, nrhs columns) using the FP32 factor in d_A_f32 plus
// FP64 iterative refinement against the FP64 matrix in d_A. The FP64 solution is
// written back into d_rhs, matching the dpotrs contract the callers expect.
void EEQSolverGPU::mixedSolveRefine(int N, int nrhs)
{
    const int sz = N * nrhs;
    if (m_impl->d_rhs_f32.n < sz) m_impl->d_rhs_f32.alloc(sz);
    if (m_impl->d_x_ref.n   < sz) m_impl->d_x_ref.alloc(sz);
    if (m_impl->d_b_ref.n   < sz) m_impl->d_b_ref.alloc(sz);
    if (m_impl->d_r_ref.n   < sz) m_impl->d_r_ref.alloc(sz);

    const int block = 256;
    const int grid  = (sz + block - 1) / block;

    // Preserve the FP64 RHS B (d_rhs is overwritten with x at the end).
    checkCudaEEQ(cudaMemcpyAsync(m_impl->d_b_ref.ptr, m_impl->d_rhs.ptr, sz * sizeof(double),
                                 cudaMemcpyDeviceToDevice, m_impl->stream), "copy B");

    // x0 = spotrs(L_f32, (float)B)
    k_cast_d2f<<<grid, block, 0, m_impl->stream>>>(m_impl->d_b_ref.ptr, m_impl->d_rhs_f32.ptr, sz);
    checkCusolverEEQ(
        cusolverDnSpotrs(m_impl->cusolver_handle, CUBLAS_FILL_MODE_LOWER, N, nrhs,
                         m_impl->d_A_f32.ptr, N, m_impl->d_rhs_f32.ptr, N, m_impl->d_info.ptr),
        "cusolverDnSpotrs x0");
    k_cast_f2d<<<grid, block, 0, m_impl->stream>>>(m_impl->d_rhs_f32.ptr, m_impl->d_x_ref.ptr, sz);

    // Iterative refinement: r = b − A·x (FP64), dx = A⁻¹·r (FP32), x += dx.
    const double done = 1.0, dnegone = -1.0;
    const int iters = (m_impl->m_refine_iters < 1) ? 1 : m_impl->m_refine_iters;
    for (int it = 0; it < iters; ++it) {
        checkCudaEEQ(cudaMemcpyAsync(m_impl->d_r_ref.ptr, m_impl->d_b_ref.ptr, sz * sizeof(double),
                                     cudaMemcpyDeviceToDevice, m_impl->stream), "copy B->R");
        // R = -A·x + R  (A symmetric, lower triangle; column-major)
        checkCublasEEQ(
            cublasDsymm(m_impl->cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                        N, nrhs, &dnegone, m_impl->d_A.ptr, N, m_impl->d_x_ref.ptr, N,
                        &done, m_impl->d_r_ref.ptr, N),
            "cublasDsymm residual");
        k_cast_d2f<<<grid, block, 0, m_impl->stream>>>(m_impl->d_r_ref.ptr, m_impl->d_rhs_f32.ptr, sz);
        checkCusolverEEQ(
            cusolverDnSpotrs(m_impl->cusolver_handle, CUBLAS_FILL_MODE_LOWER, N, nrhs,
                             m_impl->d_A_f32.ptr, N, m_impl->d_rhs_f32.ptr, N, m_impl->d_info.ptr),
            "cusolverDnSpotrs refine");
        k_axpy_f2d<<<grid, block, 0, m_impl->stream>>>(m_impl->d_x_ref.ptr, m_impl->d_rhs_f32.ptr, sz);
    }

    // Solution back into d_rhs (callers read z1/Z2 from d_rhs).
    checkCudaEEQ(cudaMemcpyAsync(m_impl->d_rhs.ptr, m_impl->d_x_ref.ptr, sz * sizeof(double),
                                 cudaMemcpyDeviceToDevice, m_impl->stream), "copy x->d_rhs");
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

        // --- 3+4. Factor: FP32 mixed-precision (WP-B) first, else FP64 Cholesky ---
        bool mixed_done = false;
        if (m_impl->m_mixed_precision) {
            if (mixedFactor(N) == 0) {  // FP32 Cholesky SPD; d_A kept as the FP64 matrix
                mixed_done = true;
                m_impl->m_cached_factor_mixed = true;
                m_impl->m_has_cached_factor   = true;
                m_last_N = N;
            }
            // else: FP32 factor not SPD → FP64 path below (mixedFactor left d_A intact)
        }

        if (!mixed_done) {
            m_impl->m_cached_factor_mixed = false;
            // cuSOLVER FP64 workspace query (cached for same N)
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

            // Cholesky factorization: A = L·L^T (overwrites d_A with L)
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
                // Cholesky failed (matrix not SPD) → caller falls back to CPU
                m_impl->m_has_cached_factor = false;
                return false;
            }

            m_impl->m_has_cached_factor = true;
        }
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

    // --- 7. Triangular solve: mixed (FP32 factor + FP64 refine) or FP64 dpotrs ---
    if (m_impl->m_cached_factor_mixed) {
        mixedSolveRefine(N, nrhs);
    } else {
        checkCusolverEEQ(
            cusolverDnDpotrs(m_impl->cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, nrhs,
                              m_impl->d_A.ptr, N,
                              m_impl->d_rhs.ptr, N,
                              m_impl->d_info.ptr),
            "cusolverDnDpotrs");
    }

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

    m_impl->orderAfterWorkspace(false);   // parity item 10: cx/cy/cz, alpha, gam, rhs come from the workspace streams

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

        // Factor: FP32 mixed-precision (WP-B) first, else FP64 Cholesky.
        bool mixed_done = false;
        if (m_impl->m_mixed_precision) {
            if (mixedFactor(N) == 0) {  // FP32 Cholesky SPD; d_A kept as the FP64 matrix
                mixed_done = true;
                m_impl->m_cached_factor_mixed = true;
                m_impl->m_has_cached_factor   = true;
                m_last_N = N;
            }
        }

        if (!mixed_done) {
            m_impl->m_cached_factor_mixed = false;
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

    if (m_impl->m_cached_factor_mixed) {
        mixedSolveRefine(N, nrhs);
    } else {
        checkCusolverEEQ(
            cusolverDnDpotrs(m_impl->cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, nrhs,
                              m_impl->d_A.ptr, N,
                              m_impl->d_rhs.ptr, N,
                              m_impl->d_info.ptr),
            "cusolverDnDpotrs");
    }

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
        k_eeq_reduce_sums_det<<<1, reduce_block,
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

    m_impl->orderAfterWorkspace(false);   // parity item 10

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

        // Factor: FP32 mixed-precision (WP-B) first, else FP64 Cholesky (+ LU fallback).
        bool mixed_done = false;
        if (m_impl->m_mixed_precision) {
            if (mixedFactor(N) == 0) {  // FP32 Cholesky SPD; d_A kept as the FP64 matrix
                mixed_done = true;
                m_impl->m_cached_factor_mixed = true;
                m_impl->m_has_cached_factor   = true;
                m_impl->m_using_lu            = false;
                m_last_N = N;
            }
            // else: FP32 not SPD → FP64 dpotrf/LU below (mixedFactor left d_A intact)
        }

        if (!mixed_done) {
        m_impl->m_cached_factor_mixed = false;
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
            // Step B (Claude Generated, May 2026): Cholesky failed, matrix is
            // not SPD. Mirror the CPU dispatcher (eeq_solver.cpp:1391-1430)
            // and fall back to LU with partial pivoting (cusolverDnDgetrf),
            // which handles indefinite matrices. Without this fallback, an
            // indefinite EEQ matrix (e.g. the polymer test case with
            // min_gersh ≈ -60.7) would silently produce wrong charges and
            // a ~0.98 Eh Coulomb energy error.
            m_impl->m_last_chol_info = *m_impl->h_info;

            // dpotrf clobbered d_A during the partial factorization — rebuild.
            {
                int n_lower = N * (N + 1) / 2;
                int block = 256;
                int grid = (n_lower + block - 1) / block;
                k_eeq_build_matrix<<<grid, block, 0, m_impl->stream>>>(
                    N, cx, cy, cz, d_alpha_corrected, d_gam_corrected,
                    m_impl->d_A.ptr, cutoff_sq);
            }

            if (m_impl->d_lu_pivots.n < N)
                m_impl->d_lu_pivots.alloc(N);
            int lwork_lu = 0;
            checkCusolverEEQ(
                cusolverDnDgetrf_bufferSize(m_impl->cusolver_handle, N, N,
                                             m_impl->d_A.ptr, N, &lwork_lu),
                "getrf_bufferSize");
            if (m_impl->d_lu_workspace.n < lwork_lu) {
                m_impl->d_lu_workspace.alloc(lwork_lu);
                m_impl->lu_workspace_size = lwork_lu;
            }

            checkCusolverEEQ(
                cusolverDnDgetrf(m_impl->cusolver_handle, N, N,
                                  m_impl->d_A.ptr, N,
                                  m_impl->d_lu_workspace.ptr,
                                  m_impl->d_lu_pivots.ptr,
                                  m_impl->d_info.ptr),
                "cusolverDnDgetrf");
            checkCudaEEQ(cudaMemcpyAsync(m_impl->h_info, m_impl->d_info.ptr, sizeof(int),
                                          cudaMemcpyDeviceToHost, m_impl->stream),
                         "download d_info (LU)");
            checkCudaEEQ(cudaStreamSynchronize(m_impl->stream), "sync after getrf");

            if (*m_impl->h_info != 0) {
                m_impl->m_has_cached_factor = false;
                m_impl->m_using_lu = false;
                return false;
            }
            m_impl->m_has_cached_factor = true;
            m_impl->m_using_lu = true;
        } else {
            m_impl->m_has_cached_factor = true;
            m_impl->m_using_lu = false;
        }
        }  // end if (!mixed_done)
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

    // Triangular solve: mixed (FP32 factor + FP64 refine), or FP64 Cholesky / LU.
    // Cached factor type drives dispatch; m_using_lu / m_cached_factor_mixed persist
    // with m_has_cached_factor so lazy-solve calls pick the right routine.
    if (m_impl->m_cached_factor_mixed) {
        mixedSolveRefine(N, nrhs);
    } else if (m_impl->m_using_lu) {
        checkCusolverEEQ(
            cusolverDnDgetrs(m_impl->cusolver_handle,
                              CUBLAS_OP_N, N, nrhs,
                              m_impl->d_A.ptr, N,
                              m_impl->d_lu_pivots.ptr,
                              m_impl->d_rhs.ptr, N,
                              m_impl->d_info.ptr),
            "cusolverDnDgetrs");
    } else {
        checkCusolverEEQ(
            cusolverDnDpotrs(m_impl->cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, nrhs,
                              m_impl->d_A.ptr, N,
                              m_impl->d_rhs.ptr, N,
                              m_impl->d_info.ptr),
            "cusolverDnDpotrs");
    }

    // GPU Schur: q[i] = z1[i] - Z2[i] * lambda
    //   lambda = (sum(z1) - rhs_c0) / sum(Z2)
    if (m_impl->d_sums.n < 2)
        m_impl->d_sums.alloc(2);
    cudaMemsetAsync(m_impl->d_sums.ptr, 0, 2 * sizeof(double), m_impl->stream);

    int rb = 256;
    k_eeq_reduce_sums_det<<<1, rb, 2 * rb * sizeof(double), m_impl->stream>>>(
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

bool EEQSolverGPU::isUsingLUFallback() const
{
    return m_impl->m_using_lu;
}

int EEQSolverGPU::getLastCholInfo() const
{
    return m_impl->m_last_chol_info;
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

// Claude Generated (Oct 2026): the Schur tail shared by WP7-A (Cholesky/LU) and the device ldlt.
// d_rhs holds [z1 | Z2] (N x (nfrag+1), column-major) on entry and the charges in d_rhs[0..N-1]
// on return; synchronizes the EEQ stream. Moved verbatim out of solveWithDeviceRHSAndGPUSchurGeneral.
static void finishGeneralSchur(EEQSolverGPUImpl& impl, int N, int nfrag,
                               const std::vector<double>& rhs_constraints)
{
    // GPU reduction: Cz1[f] = Σ_{i∈frag_f} z1[i], S[f,g] = Σ_{i∈frag_f} Z2[i,g].
    cudaMemsetAsync(impl.d_Cz1_general.ptr, 0, nfrag * sizeof(double), impl.stream);
    cudaMemsetAsync(impl.d_S_general.ptr,   0, (size_t)nfrag * nfrag * sizeof(double), impl.stream);
    {
        int rb = 256;
        int rg = (int)(((long long)nfrag * (nfrag + 1) + rb - 1) / rb);
        k_eeq_reduce_fragment_sums_csr<<<rg, rb, 0, impl.stream>>>(
            N, nfrag,
            impl.d_rhs.ptr,
            impl.d_frag_off.ptr, impl.d_frag_idx.ptr,
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
    if (nfrag < 1) {
        CurcumaLogger::warn("EEQ WP7-A: nfrag < 1, aborting");
        return false;
    }
    if (!impl.m_frag_topo_valid) {
        CurcumaLogger::warn("EEQ WP7-A: fragment topo invalid, aborting");
        return false;
    }
    if (impl.m_nfrag_batched != nfrag) {
        CurcumaLogger::warn(fmt::format(
            "EEQ WP7-A: nfrag mismatch (batched={}, requested={}), aborting",
            impl.m_nfrag_batched, nfrag));
        return false;
    }
    if (impl.d_atom_frag.n < natoms) {
        CurcumaLogger::warn(fmt::format(
            "EEQ WP7-A: atom_frag buffer too small ({} < {}), aborting",
            impl.d_atom_frag.n, natoms));
        return false;
    }
    if (impl.d_rhs_constraint_cols.n < natoms * nfrag) {
        CurcumaLogger::warn(fmt::format(
            "EEQ WP7-A: rhs_constraint_cols buffer too small ({} < {}), aborting",
            impl.d_rhs_constraint_cols.n, natoms * nfrag));
        return false;
    }

    const int N    = natoms;
    const int nrhs = nfrag + 1;
    impl.orderAfterWorkspace(false);   // parity item 10

    const bool do_refactor = force_refactor
                          || !impl.m_has_cached_factor
                          || (N != m_last_N);

    if (do_refactor) {
        // WP-B (Jun 2026): the nfrag>1 general path keeps the FP64 Cholesky/LU. With
        // nrhs = nfrag+1 ≈ N for many fragments, the FP64 residual GEMM of the mixed-
        // precision refinement costs as much as the solve itself (the FP32-factor win
        // only pays off when the O(N³) factor dominates, i.e. few RHS), and it doubles
        // the N×nrhs scratch. Mixed precision is applied only to the factor-dominated
        // few-fragment paths (solve / solveWithDeviceRHS / GPU-Schur nfrag=1).
        impl.m_cached_factor_mixed = false;
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
            CurcumaLogger::warn(fmt::format(
                "EEQ WP7-A: Cholesky factorization failed (info={}), attempting LU fallback",
                *impl.h_info));
            // Step B (Claude Generated, May 2026): Cholesky failed → LU fallback.
            // Same rationale as in solveWithDeviceRHSAndGPUSchur (WP5-A).
            impl.m_last_chol_info = *impl.h_info;
            {
                int n_lower = N * (N + 1) / 2;
                int block   = 256;
                int grid    = (n_lower + block - 1) / block;
                k_eeq_build_matrix<<<grid, block, 0, impl.stream>>>(
                    N, cx, cy, cz, d_alpha_corrected, d_gam_corrected,
                    impl.d_A.ptr, cutoff_sq);
            }

            if (impl.d_lu_pivots.n < N)
                impl.d_lu_pivots.alloc(N);
            int lwork_lu = 0;
            checkCusolverEEQ(
                cusolverDnDgetrf_bufferSize(impl.cusolver_handle, N, N,
                                             impl.d_A.ptr, N, &lwork_lu),
                "getrf_bufferSize (general)");
            if (impl.d_lu_workspace.n < lwork_lu) {
                impl.d_lu_workspace.alloc(lwork_lu);
                impl.lu_workspace_size = lwork_lu;
            }

            checkCusolverEEQ(
                cusolverDnDgetrf(impl.cusolver_handle, N, N,
                                  impl.d_A.ptr, N,
                                  impl.d_lu_workspace.ptr,
                                  impl.d_lu_pivots.ptr,
                                  impl.d_info.ptr),
                "cusolverDnDgetrf (general)");
            checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_info.ptr, sizeof(int),
                                          cudaMemcpyDeviceToHost, impl.stream),
                         "download d_info (general LU)");
            checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after getrf (general)");

            if (*impl.h_info != 0) {
                CurcumaLogger::warn(fmt::format(
                    "EEQ WP7-A: LU factorization failed (info={}), aborting",
                    *impl.h_info));
                impl.m_has_cached_factor = false;
                impl.m_using_lu = false;
                return false;
            }
            impl.m_has_cached_factor = true;
            impl.m_using_lu = true;
        } else {
            impl.m_has_cached_factor = true;
            impl.m_using_lu = false;
        }
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

    // In-place solve A·X = [b_atoms | C^T]. After this:
    //   d_rhs[0..N-1]                = z1
    //   d_rhs[(g+1)·N .. (g+2)·N-1]  = Z2[:,g]   for g = 0..nfrag-1
    // Cholesky (dpotrs) or LU (dgetrs) depending on cached factor type.
    if (impl.m_using_lu) {
        checkCusolverEEQ(
            cusolverDnDgetrs(impl.cusolver_handle,
                              CUBLAS_OP_N, N, nrhs,
                              impl.d_A.ptr, N,
                              impl.d_lu_pivots.ptr,
                              impl.d_rhs.ptr, N,
                              impl.d_info.ptr),
            "cusolverDnDgetrs (general)");
    } else {
        checkCusolverEEQ(
            cusolverDnDpotrs(impl.cusolver_handle,
                              CUBLAS_FILL_MODE_LOWER,
                              N, nrhs,
                              impl.d_A.ptr, N,
                              impl.d_rhs.ptr, N,
                              impl.d_info.ptr),
            "cusolverDnDpotrs (general)");
    }

    finishGeneralSchur(impl, N, nfrag, rhs_constraints);
    return true;
}

// ============================================================================
// WP7-C: GPU PCG (May 2026) — iterative O(k·N²) replacement for Cholesky.
// Reuses WP7-A reduce + Schur kernels for constraint handling.
// ============================================================================

// Run a single preconditioned CG solve A·x = b with warm-start in d_x.
// d_b and d_x must be N device doubles. Returns true on convergence (|r|<tol),
// false on stall. A and M_inv are read from impl.d_A / impl.d_pcg_M_inv.
//
// Uses cuBLAS with POINTER_MODE_DEVICE for inner products; one D2H per iter for
// |r|² convergence check. All work is on impl.stream — no extra sync until exit.
// WP7-D (Jun 2026, Claude Generated): largest fragment a GPU block-Jacobi block may have.
// Bounds the dynamic shared memory (max_frag_N doubles staged per block) and the O(N_f²)
// per-thread inner GEMV. Larger fragments fall back to the diagonal Jacobi.
// Sep 2026: 2048 was found (curcuma EEQ-Löser-Benchmark) to be far too permissive — a single
// 1410-atom outlier fragment (e.g. one polymer chain among hundreds of water molecules) still
// passes it, and its dense potrf+potri costs as much as the exact solver PCG exists to avoid.
// Kept only as the hard ceiling; the tuned default lives in the PARAM
// gpu_block_jacobi_max_frag_atoms (eeq_solver.h) and is threaded in per call — see
// buildBlockJacobiFactors's max_frag_atoms/max_nfrag parameters below.
static constexpr int GPU_BLOCK_JACOBI_MAX_NF = 2048;

// Build the per-fragment block-Jacobi inverse blocks into impl.d_A_blocks (port of the CPU
// EEQSolver::buildBlockJacobi). Builds each fragment's Coulomb block, Cholesky-factorizes it
// (potrf), inverts it (potri), and symmetrizes. Sets impl.m_pcg_block_jacobi_valid. On any
// non-SPD block or a too-large fragment, marks the PC invalid and the caller keeps the
// (always-extracted) diagonal Jacobi. Runs only at PCG refactor; amortized across MD steps.
//
// Sep 2026 (curcuma EEQ-Löser-Benchmark): this function's cost is dominated by TWO independent
// factors, benchmarked separately — max_frag_atoms guards the first (a single large fragment's
// dense potrf+potri), max_nfrag guards the second (the serial host loop's two
// cudaStreamSynchronize calls per fragment, catastrophic by itself at nfrag~1500+ even with
// uniformly tiny fragments). Both must be checked; neither alone covers a polymer_2x-scale
// solvated system. See docs/GPU_TUNING.md and eeq_solver.h's PARAM help text for measurements.
static void buildBlockJacobiFactors(EEQSolverGPUImpl& impl, int nfrag,
                                     const double* cx, const double* cy, const double* cz,
                                     const double* d_alpha, const double* d_gam,
                                     double cutoff_sq,
                                     int max_frag_atoms, int max_nfrag)
{
    impl.m_pcg_block_jacobi_valid = false;
    if (nfrag < 2 || !impl.m_frag_topo_valid) return;
    if (nfrag > max_nfrag) return;
    const int frag_atom_cap = (max_frag_atoms < GPU_BLOCK_JACOBI_MAX_NF) ? max_frag_atoms : GPU_BLOCK_JACOBI_MAX_NF;
    if (impl.m_max_frag_N <= 0 || impl.m_max_frag_N > frag_atom_cap) return;

    const int total_pairs = impl.h_frag_offsets_pair[nfrag];
    if (total_pairs <= 0) return;

    // 1. Build per-fragment Coulomb blocks (column-major) into d_A_blocks.
    //    Overwrites the WP6 batched factor cache → invalidate it.
    impl.m_frag_refactored = false;
    {
        int block = 256;
        int grid  = (total_pairs + block - 1) / block;
        k_eeq_build_fragment_matrices<<<grid, block, 0, impl.stream>>>(
            total_pairs, cx, cy, cz, d_alpha, d_gam,
            impl.d_frag_sizes.ptr, impl.d_frag_offsets_A.ptr,
            impl.d_frag_offsets_pair.ptr, impl.d_frag_atom_offsets.ptr,
            impl.d_frag_atom_map.ptr, impl.d_A_blocks.ptr, nfrag, cutoff_sq);
    }

    // 2. Ensure workspace covers both potrf and potri for the largest fragment.
    {
        int lw_tri = 0, lw_inv = 0;
        checkCusolverEEQ(
            cusolverDnDpotrf_bufferSize(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER,
                                        impl.m_max_frag_N, impl.d_A_blocks.ptr,
                                        impl.m_max_frag_N, &lw_tri),
            "bj potrf_bufferSize");
        checkCusolverEEQ(
            cusolverDnDpotri_bufferSize(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER,
                                        impl.m_max_frag_N, impl.d_A_blocks.ptr,
                                        impl.m_max_frag_N, &lw_inv),
            "bj potri_bufferSize");
        int lw = (lw_tri > lw_inv) ? lw_tri : lw_inv;
        if (impl.d_frag_workspace.n < lw) { impl.d_frag_workspace.alloc(lw); impl.m_frag_ws_size = lw; }
        else if (impl.m_frag_ws_size < lw) impl.m_frag_ws_size = lw;
    }

    // 3. Per-fragment potrf + potri → explicit inverse in the lower triangle.
    for (int f = 0; f < nfrag; ++f) {
        int     Nf    = impl.h_frag_sizes[f];
        if (Nf <= 0) continue;
        double* d_Af  = impl.d_A_blocks.ptr + impl.h_frag_offsets_A[f];

        checkCusolverEEQ(
            cusolverDnDpotrf(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER, Nf,
                             d_Af, Nf, impl.d_frag_workspace.ptr, impl.m_frag_ws_size,
                             impl.d_frag_info.ptr),
            "bj potrf");
        checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_frag_info.ptr, sizeof(int),
                                      cudaMemcpyDeviceToHost, impl.stream), "bj D2H info potrf");
        checkCudaEEQ(cudaStreamSynchronize(impl.stream), "bj sync potrf");
        if (*impl.h_info != 0) return;   // not SPD → keep diagonal Jacobi

        checkCusolverEEQ(
            cusolverDnDpotri(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER, Nf,
                             d_Af, Nf, impl.d_frag_workspace.ptr, impl.m_frag_ws_size,
                             impl.d_frag_info.ptr),
            "bj potri");
        checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_frag_info.ptr, sizeof(int),
                                      cudaMemcpyDeviceToHost, impl.stream), "bj D2H info potri");
        checkCudaEEQ(cudaStreamSynchronize(impl.stream), "bj sync potri");
        if (*impl.h_info != 0) return;   // singular → keep diagonal Jacobi
    }

    // 4. Mirror lower → upper so the apply can do a full symmetric GEMV.
    k_eeq_symmetrize_blocks<<<nfrag, 256, 0, impl.stream>>>(
        nfrag, impl.d_A_blocks.ptr, impl.d_frag_sizes.ptr, impl.d_frag_offsets_A.ptr);

    impl.m_pcg_block_jacobi_valid = true;
}

// WP7-D: apply the PCG preconditioner — block-Jacobi when valid, else diagonal Jacobi.
static inline void applyPrecondPCG(EEQSolverGPUImpl& impl, int N,
                                    const double* d_r, double* d_z)
{
    if (impl.m_pcg_block_jacobi_valid) {
        size_t shmem = (size_t)impl.m_max_frag_N * sizeof(double);
        k_eeq_block_jacobi_apply<<<impl.m_nfrag_batched, 256, shmem, impl.stream>>>(
            impl.m_nfrag_batched, impl.d_A_blocks.ptr,
            impl.d_frag_sizes.ptr, impl.d_frag_offsets_A.ptr,
            impl.d_frag_atom_offsets.ptr, impl.d_frag_atom_map.ptr, d_r, d_z);
    } else {
        int blk = 256, grd = (N + blk - 1) / blk;
        k_pcg_apply_precond<<<grd, blk, 0, impl.stream>>>(N, impl.d_pcg_M_inv.ptr, d_r, d_z);
    }
}

static bool runSinglePCG(EEQSolverGPUImpl& impl, int N,
                          const double* d_b, double* d_x,
                          int max_iter, double tol)
{
    // Resize per-call scratch if N grew.
    if (impl.d_pcg_r.n  < N) impl.d_pcg_r.alloc(N);
    if (impl.d_pcg_z.n  < N) impl.d_pcg_z.alloc(N);
    if (impl.d_pcg_p.n  < N) impl.d_pcg_p.alloc(N);
    if (impl.d_pcg_Ap.n < N) impl.d_pcg_Ap.alloc(N);
    if (impl.d_pcg_dot_scratch.n   < 2) impl.d_pcg_dot_scratch.alloc(2);
    if (impl.d_pcg_rnorm_scratch.n < 1) impl.d_pcg_rnorm_scratch.alloc(1);

    cublasHandle_t blas = impl.cublas_handle;
    cudaStream_t   s    = impl.stream;
    const double  one   = 1.0;
    const double  zero  = 0.0;

    // --- Init: r = b − A·x, z = M_inv·r, p = z, rz = r·z ---
    // Compute A·x → d_pcg_Ap (reuse Ap scratch as Ax holder for init).
    checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "set ptr mode host");
    checkCublasEEQ(
        cublasDsymv(blas, CUBLAS_FILL_MODE_LOWER, N,
                    &one, impl.d_A.ptr, N, d_x, 1, &zero, impl.d_pcg_Ap.ptr, 1),
        "cublasDsymv init Ax");

    {
        int blk = 256, grd = (N + blk - 1) / blk;
        k_pcg_init_residual<<<grd, blk, 0, s>>>(N, d_b, impl.d_pcg_Ap.ptr, impl.d_pcg_r.ptr);
        applyPrecondPCG(impl, N, impl.d_pcg_r.ptr, impl.d_pcg_z.ptr);  // WP7-D: block-Jacobi or diagonal
        // Initial p ← z (use cublasDcopy for simplicity).
        checkCublasEEQ(cublasDcopy(blas, N, impl.d_pcg_z.ptr, 1, impl.d_pcg_p.ptr, 1),
                       "cublasDcopy p=z");
    }

    // rz = r·z (device pointer) ; |r|² (device pointer)
    checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "set ptr mode device");
    checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_z.ptr, 1,
                              impl.d_pcg_dot_scratch.ptr /*[0]=rz*/),
                   "cublasDdot rz init");

    // Pull rz to host once for the loop (cheap; 8 bytes).
    double h_rz = 0.0;
    checkCudaEEQ(cudaMemcpyAsync(&h_rz, impl.d_pcg_dot_scratch.ptr, sizeof(double),
                                  cudaMemcpyDeviceToHost, s),
                 "D2H rz init");
    checkCudaEEQ(cudaStreamSynchronize(s), "sync rz init");

    bool converged = false;
    int  iters     = 0;
    for (int k = 0; k < max_iter; ++k) {
        iters = k + 1;

        // Ap = A·p
        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "ptr host");
        checkCublasEEQ(
            cublasDsymv(blas, CUBLAS_FILL_MODE_LOWER, N,
                        &one, impl.d_A.ptr, N, impl.d_pcg_p.ptr, 1,
                        &zero, impl.d_pcg_Ap.ptr, 1),
            "cublasDsymv Ap");

        // pAp = p·Ap (device pointer)
        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev");
        checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_p.ptr, 1, impl.d_pcg_Ap.ptr, 1,
                                  impl.d_pcg_dot_scratch.ptr + 1 /*[1]=pAp*/),
                       "cublasDdot pAp");

        // D2H pAp (8 bytes).
        double h_pAp = 0.0;
        checkCudaEEQ(cudaMemcpyAsync(&h_pAp, impl.d_pcg_dot_scratch.ptr + 1, sizeof(double),
                                      cudaMemcpyDeviceToHost, s),
                     "D2H pAp");
        checkCudaEEQ(cudaStreamSynchronize(s), "sync pAp");

        if (std::abs(h_pAp) < 1e-30) break;  // degenerate direction
        double alpha = h_rz / h_pAp;
        double neg_alpha = -alpha;

        // x += α·p ; r -= α·Ap
        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "ptr host");
        checkCublasEEQ(cublasDaxpy(blas, N, &alpha, impl.d_pcg_p.ptr,  1, d_x, 1),
                       "Daxpy x+=ap");
        checkCublasEEQ(cublasDaxpy(blas, N, &neg_alpha, impl.d_pcg_Ap.ptr, 1,
                                    impl.d_pcg_r.ptr, 1),
                       "Daxpy r-=aAp");

        // |r|² = r·r (device pointer)
        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev");
        checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_r.ptr, 1,
                                  impl.d_pcg_rnorm_scratch.ptr),
                       "Ddot rnorm");
        double h_rnorm_sq = 0.0;
        checkCudaEEQ(cudaMemcpyAsync(&h_rnorm_sq, impl.d_pcg_rnorm_scratch.ptr,
                                      sizeof(double), cudaMemcpyDeviceToHost, s),
                     "D2H rnorm");
        checkCudaEEQ(cudaStreamSynchronize(s), "sync rnorm");

        if (h_rnorm_sq < tol * tol) { converged = true; break; }

        // z_new = M_inv · r ; rz_new = r · z_new   (WP7-D: block-Jacobi or diagonal)
        applyPrecondPCG(impl, N, impl.d_pcg_r.ptr, impl.d_pcg_z.ptr);
        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev");
        checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_z.ptr, 1,
                                  impl.d_pcg_dot_scratch.ptr /*[0]=rz_new*/),
                       "Ddot rz_new");
        double h_rz_new = 0.0;
        checkCudaEEQ(cudaMemcpyAsync(&h_rz_new, impl.d_pcg_dot_scratch.ptr, sizeof(double),
                                      cudaMemcpyDeviceToHost, s),
                     "D2H rz_new");
        checkCudaEEQ(cudaStreamSynchronize(s), "sync rz_new");

        if (std::abs(h_rz) < 1e-30) break;  // degenerate
        double beta = h_rz_new / h_rz;
        h_rz = h_rz_new;

        // p = z + β·p (in-place via dedicated kernel — avoids cublasDcopy + axpy).
        {
            int blk = 256, grd = (N + blk - 1) / blk;
            k_pcg_dir_update<<<grd, blk, 0, s>>>(N, impl.d_pcg_z.ptr, beta,
                                                  impl.d_pcg_p.ptr, impl.d_pcg_p.ptr);
        }
    }

    impl.m_pcg_total_calls += 1;
    impl.m_pcg_total_iters += iters;
    if (!converged) impl.m_pcg_nonconv_calls += 1;
    return converged;
}

bool EEQSolverGPU::solveWithDeviceRHSAndGPUPCG(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha_corrected,
    const double* d_gam_corrected,
    const double* d_rhs_atoms,
    const std::vector<int>& /*fraglist*/,
    const std::vector<double>& rhs_constraints,
    int    max_iter,
    double tol,
    double cutoff_sq,
    bool   force_refactor,
    int    block_jacobi_max_frag_atoms,
    int    block_jacobi_max_nfrag)
{
    auto& impl = *m_impl;
    if (nfrag < 1) return false;
    if (!impl.m_frag_topo_valid) return false;
    if (impl.m_nfrag_batched != nfrag) return false;
    if (impl.d_atom_frag.n < natoms) return false;
    if (impl.d_rhs_constraint_cols.n < natoms * nfrag) return false;

    const int N    = natoms;
    const int nrhs = nfrag + 1;
    impl.orderAfterWorkspace(false);   // parity item 10

    const bool do_refactor = force_refactor
                          || !impl.m_pcg_M_inv_valid
                          || (N != m_last_N);

    // Claude Generated (Oct 2026): per-call breakdown for the verbosity-2 diagnostic below
    // (setup = matrix + preconditioner, solves = the nfrag+1 PCG runs). Host wall clock; the
    // setup is synchronous through buildBlockJacobiFactors or the first PCG's D2H.
    using pcg_clock = std::chrono::steady_clock;
    const auto pcg_t0 = pcg_clock::now();
    const int pcg_iters0 = impl.m_pcg_total_iters;
    double pcg_setup_ms = 0.0;
    int pcg_iter_z1 = 0, pcg_iter_max = 0;

    if (do_refactor) {
        // (Re)build A on GPU — same kernel WP7-A uses for the Cholesky path.
        {
            int n_lower = N * (N + 1) / 2;
            int block   = 256;
            int grid    = (n_lower + block - 1) / block;
            k_eeq_build_matrix<<<grid, block, 0, impl.stream>>>(
                N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, impl.d_A.ptr, cutoff_sq);
        }
        if (impl.d_pcg_M_inv.n < N) impl.d_pcg_M_inv.alloc(N);
        {
            int blk = 256, grd = (N + blk - 1) / blk;
            k_pcg_extract_diag_inv<<<grd, blk, 0, impl.stream>>>(
                N, impl.d_A.ptr, impl.d_pcg_M_inv.ptr);
        }
        impl.m_pcg_M_inv_valid = true;
        // WP7-D (Jun 2026): build the per-fragment block-Jacobi inverse (nfrag>=2). On
        // success the PCG applies the exact per-fragment inverse instead of the diagonal
        // Jacobi (sets m_pcg_block_jacobi_valid); on failure the diagonal above is kept.
        // Sep 2026: size-gated (see buildBlockJacobiFactors) — skipped entirely above
        // block_jacobi_max_nfrag or when the largest fragment exceeds
        // block_jacobi_max_frag_atoms, falling back to the diagonal Jacobi (always correct,
        // just more PCG iterations). See eeq_solver.h PARAM docs for the measurements.
        buildBlockJacobiFactors(impl, nfrag, cx, cy, cz,
                                d_alpha_corrected, d_gam_corrected, cutoff_sq,
                                block_jacobi_max_frag_atoms, block_jacobi_max_nfrag);
        if (impl.m_pcg_block_jacobi_valid) {
            CurcumaLogger::info("EEQ GPU PCG: block-Jacobi preconditioner active");
        } else if (nfrag >= 2) {
            CurcumaLogger::info(fmt::format(
                "EEQ GPU PCG: block-Jacobi skipped (nfrag={}, max_frag_N={}) — "
                "using diagonal Jacobi (gpu_block_jacobi_max_nfrag={}, "
                "gpu_block_jacobi_max_frag_atoms={})",
                nfrag, impl.m_max_frag_N, block_jacobi_max_nfrag, block_jacobi_max_frag_atoms));
        }
        checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after pcg setup");
        pcg_setup_ms = std::chrono::duration<double, std::milli>(pcg_clock::now() - pcg_t0).count();
        // Note: do NOT touch m_last_N here — it is owned by the cuSOLVER paths
        // (WP5-A / WP7-A) and tracks their workspace allocation. PCG bypasses
        // cuSOLVER entirely, so writing m_last_N here would falsely tell the
        // Cholesky paths their workspace is still valid for the current N.
        // Geometry change ⇒ warm-start may be far from new optimum, but PCG still
        // converges from any x0. Keep the cached z1/Z2 anyway — Z2 in particular
        // is very stable across MD steps (matches CPU m_pcg_last_z2 strategy).
    }

    // Allocate d_rhs (column-major [N × (nfrag+1)]) — first column receives z1
    // after PCG, columns 1..nfrag receive Z2. Same layout as WP7-A.
    {
        const int rhs_size = N * nrhs;
        if (impl.d_rhs.n < rhs_size) impl.d_rhs.alloc(rhs_size);
    }

    // --- PCG #1: A · z1 = b_atoms.   x0 ← d_z1_persistent on warm; else zero. ---
    if (impl.m_pcg_warm_valid) {
        // d_rhs[0..N-1] ← d_z1_persistent  (D2D copy, on stream)
        checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, impl.d_z1_persistent.ptr,
                                      N * sizeof(double),
                                      cudaMemcpyDeviceToDevice, impl.stream),
                     "D2D z1 warm-start");
    } else {
        cudaMemsetAsync(impl.d_rhs.ptr, 0, N * sizeof(double), impl.stream);
    }
    bool ok = runSinglePCG(impl, N, d_rhs_atoms, impl.d_rhs.ptr, max_iter, tol);
    if (!ok) return false;
    pcg_iter_z1 = impl.m_pcg_total_iters - pcg_iters0;

    // --- PCG #2..nfrag+1: A · Z2[:,f] = e_f for each fragment. ---
    for (int f = 0; f < nfrag; ++f) {
        double* d_x_col = impl.d_rhs.ptr + (f + 1) * N;
        const double* d_b_col = impl.d_rhs_constraint_cols.ptr + f * N;
        if (impl.m_pcg_warm_valid && impl.d_Z2_persistent.n >= (f + 1) * N) {
            checkCudaEEQ(cudaMemcpyAsync(d_x_col, impl.d_Z2_persistent.ptr + f * N,
                                          N * sizeof(double),
                                          cudaMemcpyDeviceToDevice, impl.stream),
                         "D2D Z2 warm-start");
        } else {
            cudaMemsetAsync(d_x_col, 0, N * sizeof(double), impl.stream);
        }
        const int it_before = impl.m_pcg_total_iters;
        ok = runSinglePCG(impl, N, d_b_col, d_x_col, max_iter, tol);
        if (!ok) return false;
        pcg_iter_max = std::max(pcg_iter_max, impl.m_pcg_total_iters - it_before);
    }
    if (CurcumaLogger::get_verbosity() >= 2) {
        const double total_ms = std::chrono::duration<double, std::milli>(pcg_clock::now() - pcg_t0).count();
        const int iters = impl.m_pcg_total_iters - pcg_iters0;
        CurcumaLogger::info(fmt::format(
            "EEQ GPU PCG (WP7-C): {} solves, {} iterations (z1 {}, max per Z2 column {}, mean {:.1f}), "
            "setup {:.1f} ms, solves {:.1f} ms, preconditioner {}",
            nrhs, iters, pcg_iter_z1, pcg_iter_max, double(iters) / nrhs, pcg_setup_ms,
            total_ms - pcg_setup_ms, impl.m_pcg_block_jacobi_valid ? "block-Jacobi" : "diagonal"));
    }

    // --- Persist warm-start BEFORE the Schur apply overwrites d_rhs[0..N-1]. ---
    checkCudaEEQ(cudaMemcpyAsync(impl.d_z1_persistent.ptr, impl.d_rhs.ptr,
                                  N * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream),
                 "D2D persist z1");
    if (impl.d_Z2_persistent.n < N * nfrag) impl.d_Z2_persistent.alloc(N * nfrag);
    checkCudaEEQ(cudaMemcpyAsync(impl.d_Z2_persistent.ptr, impl.d_rhs.ptr + N,
                                  N * nfrag * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream),
                 "D2D persist Z2");
    impl.m_pcg_warm_valid = true;

    // --- Schur reduction (reuse WP7-A kernels): Cz1[nfrag] + S[nfrag×nfrag] ---
    cudaMemsetAsync(impl.d_Cz1_general.ptr, 0, nfrag * sizeof(double), impl.stream);
    cudaMemsetAsync(impl.d_S_general.ptr,   0, (size_t)nfrag * nfrag * sizeof(double), impl.stream);
    {
        int rb = 256, rg = (int)(((long long)nfrag * (nfrag + 1) + rb - 1) / rb);
        k_eeq_reduce_fragment_sums_csr<<<rg, rb, 0, impl.stream>>>(
            N, nfrag, impl.d_rhs.ptr, impl.d_frag_off.ptr, impl.d_frag_idx.ptr,
            impl.d_Cz1_general.ptr, impl.d_S_general.ptr);
    }

    checkCudaEEQ(cudaMemcpyAsync(impl.h_Cz1_general.data(), impl.d_Cz1_general.ptr,
                                  nfrag * sizeof(double),
                                  cudaMemcpyDeviceToHost, impl.stream),
                 "D2H Cz1 (pcg)");
    checkCudaEEQ(cudaMemcpyAsync(impl.h_S_general.data(), impl.d_S_general.ptr,
                                  (size_t)nfrag * nfrag * sizeof(double),
                                  cudaMemcpyDeviceToHost, impl.stream),
                 "D2H S (pcg)");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after reduce (pcg)");

    // CPU Schur: S·λ = Cz1 − Q_frag.
    std::vector<double> schur_rhs(nfrag, 0.0);
    for (int f = 0; f < nfrag; ++f) {
        double q_target = (f < (int)rhs_constraints.size()) ? rhs_constraints[f] : 0.0;
        schur_rhs[f] = impl.h_Cz1_general[f] - q_target;
    }
    solveSchurSystemInPlace(impl.h_S_general.data(), schur_rhs.data(),
                             impl.h_lambda_general.data(), nfrag);

    checkCudaEEQ(cudaMemcpyAsync(impl.d_lambda_general.ptr,
                                  impl.h_lambda_general.data(),
                                  nfrag * sizeof(double),
                                  cudaMemcpyHostToDevice, impl.stream),
                 "H2D lambda (pcg)");

    // Apply: q[i] = z1[i] − Σ_g Z2[i,g]·λ[g]; in-place into d_rhs[0..N-1].
    {
        int b = 256, g = (N + b - 1) / b;
        k_eeq_schur_general<<<g, b, 0, impl.stream>>>(
            N, nfrag, impl.d_rhs.ptr, impl.d_lambda_general.ptr, impl.d_rhs.ptr);
    }

    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after schur_general (pcg)");
    return true;
}

// ============================================================================
// WP7-E: GPU projected PCG (Sep 2026)
// Single-solve port of the CPU EEQSolver::solveWithProjectedPCG: instead of
// (nfrag+1) independent PCG solves + a CPU Schur reduction (WP7-C above), one
// PCG runs directly in the constraint tangent space (Σ_{i∈frag} v_i = 0) and
// returns the final charges without any Schur step. Motivated by the fact
// that WP7-C's per-fragment loop (runSinglePCG called nfrag+1 times) and its
// block-Jacobi setup (buildBlockJacobiFactors) both scale with fragment count
// / size in ways that make WP7-C intractable for many-fragment, large systems
// (e.g. a solvated polymer with nfrag~1500) — see the "curcuma EEQ-Löser-
// Benchmark" project (docs / vault). This path uses only the diagonal Jacobi,
// exactly mirroring the CPU reference; there is no block-Jacobi cost to skip.
// Claude Generated.
// ============================================================================

// Deterministic fragment sum out[f] = sum_{i in f} v[i] (fixed order; replaces the atomicAdd k_frag_sum).
static inline void fragSumDet(EEQSolverGPUImpl& impl, int nfrag, const double* d_v, double* d_out, cudaStream_t s)
{
    // many small fragments: one warp each (8 fragments per block); few large ones: one block each
    if (static_cast<long long>(impl.d_frag_idx.n) <= 64LL * nfrag)
        k_frag_sum_warp<<<(nfrag + 7) / 8, 256, 0, s>>>(nfrag, d_v, impl.d_frag_off.ptr, impl.d_frag_idx.ptr, d_out);
    else
        k_frag_sum_csr<<<nfrag, 256, 256 * sizeof(double), s>>>(nfrag, d_v, impl.d_frag_off.ptr, impl.d_frag_idx.ptr, d_out);
}

// r ← P·r, the projection onto the constraint tangent space: subtract the
// per-fragment mean (mirrors the CPU `project` lambda).
static inline void ppcgProjectResidual(EEQSolverGPUImpl& impl, int N, int nfrag, double* d_v)
{
    cudaMemsetAsync(impl.d_frag_sum_scratch.ptr, 0, (size_t)nfrag * sizeof(double), impl.stream);
    {
        fragSumDet(impl, nfrag, d_v, impl.d_frag_sum_scratch.ptr, impl.stream);
    }
    {
        int blk = 256, grd = (nfrag + blk - 1) / blk;
        k_frag_project_delta<<<grd, blk, 0, impl.stream>>>(
            nfrag, impl.d_frag_sum_scratch.ptr, impl.d_frag_inv_atoms.ptr,
            impl.d_frag_delta_scratch.ptr);
    }
    {
        int blk = 256, grd = (N + blk - 1) / blk;
        k_frag_scatter_add<<<grd, blk, 0, impl.stream>>>(
            N, d_v, impl.d_atom_frag.ptr, impl.d_frag_delta_scratch.ptr);
    }
}

// z ← M_inv·r, then an M-weighted per-fragment correction so that C·z = 0
// (mirrors the CPU `precondition` lambda).
// Claude Generated (Oct 2026): kernels of the local-inverse preconditioner (z = G r - H S^-1 C G r).
// One thread per row, fixed summation order inside a row.
__global__ void k_li_csr_spmv(int nrows, const int* __restrict__ rp, const int* __restrict__ col,
                              const double* __restrict__ val, const double* __restrict__ x,
                              double* __restrict__ y, int subtract)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nrows) return;
    double s = 0.0;
    for (int e = rp[i]; e < rp[i + 1]; ++e) s += val[e] * x[col[e]];
    if (subtract) y[i] -= s; else y[i] = s;
}
// mu = S^-1 t, S^-1 dense column-major (symmetric), one thread per row.
__global__ void k_li_dense_gemv(int n, const double* __restrict__ Sinv, const double* __restrict__ t,
                                double* __restrict__ mu)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double s = 0.0;
    for (int j = 0; j < n; ++j) s += Sinv[static_cast<size_t>(j) * n + i] * t[j];
    mu[i] = s;
}
// Breakdown flag for the local inverse: rz <= 0 means the preconditioner is not positive definite.
__global__ void k_li_check_rz(const double* __restrict__ rz, double* flag)
{
    if (blockIdx.x == 0 && threadIdx.x == 0 && !(*rz > 0.0)) *flag = 1.0;
}

static inline void ppcgPrecondition(EEQSolverGPUImpl& impl, int N, int nfrag,
                                     const double* d_r, double* d_z)
{
    int blk = 256, grd = (N + blk - 1) / blk;
    if (impl.m_li_valid && impl.m_li_N == N && impl.m_li_nfrag == nfrag && impl.m_li_device) {
        // Device build (Oct 2026): z = P G r, P = I - C^T (C C^T)^-1 C (subtract the fragment mean).
        // Symmetric and positive definite on null(C) whenever G is; C z = 0 afterwards.
        k_li_csr_spmv<<<grd, blk, 0, impl.stream>>>(N, impl.d_li_G_rowptr.ptr, impl.d_li_G_col.ptr,
                                                    impl.d_li_G_val.ptr, d_r, d_z, 0);
        ppcgProjectResidual(impl, N, nfrag, d_z);
        return;
    }
    if (impl.m_li_valid && impl.m_li_N == N && impl.m_li_nfrag == nfrag) {
        // z = G r ; t = C z ; mu = S^-1 t ; z -= H mu  (C z = 0 afterwards)
        k_li_csr_spmv<<<grd, blk, 0, impl.stream>>>(N, impl.d_li_G_rowptr.ptr, impl.d_li_G_col.ptr,
                                                    impl.d_li_G_val.ptr, d_r, d_z, 0);
        cudaMemsetAsync(impl.d_li_t.ptr, 0, (size_t)nfrag * sizeof(double), impl.stream);
        fragSumDet(impl, nfrag, d_z, impl.d_li_t.ptr, impl.stream);
        const int grd_f = (nfrag + blk - 1) / blk;
        k_li_dense_gemv<<<grd_f, blk, 0, impl.stream>>>(nfrag, impl.d_li_S_inv.ptr, impl.d_li_t.ptr, impl.d_li_mu.ptr);
        k_li_csr_spmv<<<grd, blk, 0, impl.stream>>>(N, impl.d_li_H_rowptr.ptr, impl.d_li_H_col.ptr,
                                                    impl.d_li_H_val.ptr, impl.d_li_mu.ptr, d_z, 1);
        return;
    }
    k_pcg_apply_precond<<<grd, blk, 0, impl.stream>>>(N, impl.d_pcg_M_inv.ptr, d_r, d_z);
    cudaMemsetAsync(impl.d_frag_sum_scratch.ptr, 0, (size_t)nfrag * sizeof(double), impl.stream);
    fragSumDet(impl, nfrag, d_z, impl.d_frag_sum_scratch.ptr, impl.stream);
    k_frag_precond_correct<<<grd, blk, 0, impl.stream>>>(
        N, impl.d_pcg_M_inv.ptr, impl.d_atom_frag.ptr,
        impl.d_frag_sum_scratch.ptr, impl.d_frag_sM.ptr, d_z);
}

// ============================================================================
// Multi-GPU split of the projected-PCG EEQ (Claude Generated, Sep 2026, F-3)
// ============================================================================

static bool splitWanted(const EEQSolverGPUImpl& impl, int N)
{
    return !impl.split_devs.empty() && !impl.split_failed && N >= impl.split_min_atoms;
}

/// Allocate the column blocks for N atoms (once per N). false = split unusable (status set).
static bool splitPrepare(EEQSolverGPUImpl& impl, int N)
{
    if (impl.split_N == N && !impl.split_parts.empty()) return true;
    impl.releaseSplit();
    { const char* e = std::getenv("CURCUMA_EEQ_MIXED"); impl.split_mixed = impl.split_mixed_req && !(e && *e == '0'); }
    const int nparts = static_cast<int>(impl.split_devs.size()) + 1;
    // Equal column shares; this device takes the first block.
    std::vector<int> bounds(nparts + 1, 0);
    for (int k = 0; k <= nparts; ++k)
        bounds[k] = static_cast<int>((static_cast<long long>(N) * k) / nparts);
    try {
        for (int h = 0; h < static_cast<int>(impl.split_devs.size()); ++h) {
            auto part = std::make_unique<EEQSolverGPUImpl::SplitPart>();
            part->dev = impl.split_devs[h];
            part->c0 = bounds[h + 1];
            part->c1 = bounds[h + 2];
            const int nc = part->c1 - part->c0;
            if (cudaSetDevice(part->dev) != cudaSuccess)
                throw std::runtime_error("cudaSetDevice(" + std::to_string(part->dev) + ") failed");
            int can = 0;
            cudaDeviceCanAccessPeer(&can, part->dev, impl.own_device);
            if (can) cudaDeviceEnablePeerAccess(impl.own_device, 0);
            cudaGetLastError();   // "already enabled" is fine
            // Oct 2026: highest priority: the matvec chain of the PCG is the critical path and must overtake
            // the dispersion helper kernels queued on the same device.
            int prio_lo = 0, prio_hi = 0;
            cudaDeviceGetStreamPriorityRange(&prio_lo, &prio_hi);
            if (cudaStreamCreateWithPriority(&part->stream, cudaStreamDefault, prio_hi) != cudaSuccess
                || cudaEventCreateWithFlags(&part->done, cudaEventDisableTiming) != cudaSuccess
                || cublasCreate(&part->blas) != CUBLAS_STATUS_SUCCESS
                || cublasSetStream(part->blas, part->stream) != CUBLAS_STATUS_SUCCESS)
                throw std::runtime_error("stream/event/cuBLAS setup failed on device " + std::to_string(part->dev));
            part->A.alloc(N * nc);
            part->x.alloc(N); part->y.alloc(N); part->z.alloc(N);
            part->alpha.alloc(N); part->gam.alloc(N);
            part->p.alloc(N); part->yblk.alloc(nc);
            if (impl.split_mixed) { part->Af.alloc(N * nc); part->xf.alloc(N); part->yf.alloc(nc); }
            impl.split_parts.push_back(std::move(part));
        }
        cudaSetDevice(impl.own_device);
        for (int d : impl.split_devs) {
            int can = 0;
            cudaDeviceCanAccessPeer(&can, impl.own_device, d);
            if (can) cudaDeviceEnablePeerAccess(d, 0);
        }
        cudaGetLastError();
        impl.split_c1_own = bounds[1];
        impl.d_A_own.alloc(N * impl.split_c1_own);
        if (impl.split_mixed) { impl.d_A_own_f.alloc(N * impl.split_c1_own); impl.d_xf.alloc(N); impl.d_yf.alloc(N); }
        if (!impl.ev_split) cudaEventCreateWithFlags(&impl.ev_split, cudaEventDisableTiming);
        impl.split_N = N;
    } catch (const std::exception& e) {
        impl.releaseSplit();
        cudaSetDevice(impl.own_device);
        impl.split_failed = true;
        impl.split_status = std::string("unavailable: ") + e.what();
        return false;
    }
    return true;
}

/// Build every column block (this device's on impl.stream, the helpers' on their streams).
static void splitBuild(EEQSolverGPUImpl& impl, int N, const double* cx, const double* cy,
                       const double* cz, const double* alpha, const double* gam, double cutoff_sq)
{
    const size_t bytes = sizeof(double) * static_cast<size_t>(N);
    cudaEventRecord(impl.ev_split, impl.stream);
    for (auto& pp : impl.split_parts) {
        auto& P = *pp;
        cudaSetDevice(P.dev);
        cudaStreamWaitEvent(P.stream, impl.ev_split, 0);
        cudaMemcpyPeerAsync(P.x.ptr, P.dev, cx, impl.own_device, bytes, P.stream);
        cudaMemcpyPeerAsync(P.y.ptr, P.dev, cy, impl.own_device, bytes, P.stream);
        cudaMemcpyPeerAsync(P.z.ptr, P.dev, cz, impl.own_device, bytes, P.stream);
        cudaMemcpyPeerAsync(P.alpha.ptr, P.dev, alpha, impl.own_device, bytes, P.stream);
        cudaMemcpyPeerAsync(P.gam.ptr, P.dev, gam, impl.own_device, bytes, P.stream);
        const int nc = P.c1 - P.c0;
        const int blk = 256;
        const dim3 grd((N + blk - 1) / blk, std::min(nc, 65535));
        k_eeq_build_cols<<<grd, blk, 0, P.stream>>>(N, P.c0, nc, P.x.ptr, P.y.ptr, P.z.ptr,
                                                     P.alpha.ptr, P.gam.ptr, P.A.ptr, cutoff_sq,
                                                     impl.split_mixed ? P.Af.ptr : nullptr);
    }
    cudaSetDevice(impl.own_device);
    const int blk = 256;
    const dim3 grd((N + blk - 1) / blk, std::min(impl.split_c1_own, 65535));
    k_eeq_build_cols<<<grd, blk, 0, impl.stream>>>(N, 0, impl.split_c1_own, cx, cy, cz,
                                                   alpha, gam, impl.d_A_own.ptr, cutoff_sq,
                                                   impl.split_mixed ? impl.d_A_own_f.ptr : nullptr);
    ++impl.split_builds;
}

/// Column dot products in a fixed order (Claude Generated, Oct 2026): y[c] = sum_i A[i + c N] x[i], one block per
/// column, thread t sums i = t, t+256, ... and a fixed tree adds the 256 partials. cuBLAS gemv picks its algorithm
/// from the block shape (nc), so the result of a column could depend on the number of GPUs; this kernel does not.
__global__ void k_split_coldot_f64(int N, const double* __restrict__ A, const double* __restrict__ x, double* __restrict__ y)
{
    __shared__ double sh[256];
    const double* col = A + static_cast<size_t>(blockIdx.x) * N;
    double acc = 0.0;
    for (int i = threadIdx.x; i < N; i += 256) acc = fma(col[i], x[i], acc);
    sh[threadIdx.x] = acc;
    __syncthreads();
    for (int s = 128; s > 0; s >>= 1) {
        if (threadIdx.x < s) sh[threadIdx.x] += sh[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) y[blockIdx.x] = sh[0];
}

/// FP32 matrix and vector, float partials per thread, the tree in double (the sum of 256 values).
__global__ void k_split_coldot_f32(int N, const float* __restrict__ A, const float* __restrict__ x, float* __restrict__ y)
{
    __shared__ double sh[256];
    const float* col = A + static_cast<size_t>(blockIdx.x) * N;
    float acc = 0.0f;
    for (int i = threadIdx.x; i < N; i += 256) acc = fmaf(col[i], x[i], acc);
    sh[threadIdx.x] = static_cast<double>(acc);
    __syncthreads();
    for (int s = 128; s > 0; s >>= 1) {
        if (threadIdx.x < s) sh[threadIdx.x] += sh[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) y[blockIdx.x] = static_cast<float>(sh[0]);
}

/// y = A x with the column blocks spread over the devices. x and y live on this device.
static void splitMatvecF32(EEQSolverGPUImpl& impl, int N, const double* x, double* y)
{
    // Oct 2026: FP32 matrix, FP32 vectors, FP64 result vector. Halves the matrix traffic and the peer copies.
    const int blk = 256, grd = (N + blk - 1) / blk;
    k_cast_d2f<<<grd, blk, 0, impl.stream>>>(x, impl.d_xf.ptr, N);
    cudaEventRecord(impl.ev_split, impl.stream);
    for (auto& pp : impl.split_parts) {
        auto& P = *pp;
        const int nc = P.c1 - P.c0;
        cudaSetDevice(P.dev);
        cudaStreamWaitEvent(P.stream, impl.ev_split, 0);
        cudaMemcpyPeerAsync(P.xf.ptr, P.dev, impl.d_xf.ptr, impl.own_device, sizeof(float) * N, P.stream);
        k_split_coldot_f32<<<nc, 256, 0, P.stream>>>(N, P.Af.ptr, P.xf.ptr, P.yf.ptr);
        cudaMemcpyPeerAsync(impl.d_yf.ptr + P.c0, impl.own_device, P.yf.ptr, P.dev, sizeof(float) * nc, P.stream);
        cudaEventRecord(P.done, P.stream);
    }
    cudaSetDevice(impl.own_device);
    k_split_coldot_f32<<<impl.split_c1_own, 256, 0, impl.stream>>>(N, impl.d_A_own_f.ptr, impl.d_xf.ptr, impl.d_yf.ptr);
    for (auto& pp : impl.split_parts)
        cudaStreamWaitEvent(impl.stream, pp->done, 0);
    k_cast_f2d<<<grd, blk, 0, impl.stream>>>(impl.d_yf.ptr, y, N);
    ++impl.split_matvecs;
}

static void splitMatvec(EEQSolverGPUImpl& impl, int N, const double* x, double* y)
{
    cudaEventRecord(impl.ev_split, impl.stream);
    for (auto& pp : impl.split_parts) {
        auto& P = *pp;
        const int nc = P.c1 - P.c0;
        cudaSetDevice(P.dev);
        cudaStreamWaitEvent(P.stream, impl.ev_split, 0);
        cudaMemcpyPeerAsync(P.p.ptr, P.dev, x, impl.own_device, sizeof(double) * N, P.stream);
        k_split_coldot_f64<<<nc, 256, 0, P.stream>>>(N, P.A.ptr, P.p.ptr, P.yblk.ptr);
        cudaMemcpyPeerAsync(y + P.c0, impl.own_device, P.yblk.ptr, P.dev, sizeof(double) * nc, P.stream);
        cudaEventRecord(P.done, P.stream);
    }
    cudaSetDevice(impl.own_device);
    k_split_coldot_f64<<<impl.split_c1_own, 256, 0, impl.stream>>>(N, impl.d_A_own.ptr, x, y);
    for (auto& pp : impl.split_parts)
        cudaStreamWaitEvent(impl.stream, pp->done, 0);
    ++impl.split_matvecs;
}

/// A x through the split when it is set up for N, else the single-device dense symv.
static void eeqMatvec(EEQSolverGPUImpl& impl, int N, const double* x, double* y, bool exact = false)
{
    if (impl.split_N == N && !impl.split_parts.empty()) {
        if (impl.split_mixed && !exact) { splitMatvecF32(impl, N, x, y); return; }
        splitMatvec(impl, N, x, y);
        return;
    }
    const double one = 1.0, zero = 0.0;
    cublasSetPointerMode(impl.cublas_handle, CUBLAS_POINTER_MODE_HOST);
    checkCublasEEQ(cublasDsymv(impl.cublas_handle, CUBLAS_FILL_MODE_LOWER, N,
                               &one, impl.d_A.ptr, N, x, 1, &zero, y, 1),
                   "cublasDsymv (ppcg)");
}

/// q += alpha p, r -= alpha Ap with alpha = rz / pAp read on the device (Claude Generated, Sep 2026).
/// pAp <= 0 (non-descent) sets *flag and leaves q, r unchanged, as the host loop used to break.
__global__ void k_ppcg_update_qr(int N, const double* __restrict__ rz, const double* __restrict__ pAp,
                                 const double* __restrict__ p, const double* __restrict__ Ap,
                                 double* __restrict__ q, double* __restrict__ r, double* flag)
{
    const double den = *pAp;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!(den > 0.0)) {
        if (i == 0) *flag = 1.0;
        return;
    }
    if (i >= N) return;
    const double alpha = *rz / den;
    q[i] += alpha * p[i];
    r[i] -= alpha * Ap[i];
}

/// p = z + beta p with beta = rz_new / rz_old on the device; |rz_old| < 1e-30 sets *flag.
__global__ void k_ppcg_update_p(int N, const double* __restrict__ rz_new, const double* __restrict__ rz_old,
                                const double* __restrict__ z, double* __restrict__ p, double* flag)
{
    const double old = *rz_old;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (fabs(old) < 1e-30) {
        if (i == 0) *flag = 1.0;
        return;
    }
    if (i >= N) return;
    const double beta = *rz_new / old;
    p[i] = z[i] + beta * p[i];
}

// One projected PCG solve. d_q is both the (feasibility-shifted, warm-start)
// input and the final-charges output, in place. Port of
// EEQSolver::solveWithProjectedPCG (eeq_solver.cpp).
static bool runProjectedPCG(EEQSolverGPUImpl& impl, int N, int nfrag,
                             const double* d_rhs_atoms, double* d_q,
                             int max_iter, double tol)
{
    if (impl.d_pcg_r.n  < N) impl.d_pcg_r.alloc(N);
    if (impl.d_pcg_z.n  < N) impl.d_pcg_z.alloc(N);
    if (impl.d_pcg_p.n  < N) impl.d_pcg_p.alloc(N);
    if (impl.d_pcg_Ap.n < N) impl.d_pcg_Ap.alloc(N);
    if (impl.d_pcg_dot_scratch.n   < 2) impl.d_pcg_dot_scratch.alloc(2);
    impl.ppcg_true_sq = -1.0;
    if (impl.d_pcg_rnorm_scratch.n < 1) impl.d_pcg_rnorm_scratch.alloc(1);

    cublasHandle_t blas = impl.cublas_handle;
    cudaStream_t   s    = impl.stream;
    const double  one   = 1.0;
    const double  zero  = 0.0;

    // --- Feasible start: q[i] += (rhs_constraints[frag[i]] − Σ_{j∈frag} q[j]) / N_frag. ---
    cudaMemsetAsync(impl.d_frag_sum_scratch.ptr, 0, (size_t)nfrag * sizeof(double), s);
    {
        fragSumDet(impl, nfrag, d_q, impl.d_frag_sum_scratch.ptr, s);
    }
    {
        int blk = 256, grd = (nfrag + blk - 1) / blk;
        k_frag_feasibility_delta<<<grd, blk, 0, s>>>(
            nfrag, impl.d_frag_sum_scratch.ptr, impl.d_rhs_constraints_dev.ptr,
            impl.d_frag_inv_atoms.ptr, impl.d_frag_delta_scratch.ptr);
    }
    {
        int blk = 256, grd = (N + blk - 1) / blk;
        k_frag_scatter_add<<<grd, blk, 0, s>>>(N, d_q, impl.d_atom_frag.ptr, impl.d_frag_delta_scratch.ptr);
    }

    // --- r = b − A·q ; project(r) ; z = precondition(r) ; p = z ; rz = r·z ---
    checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "ptr host (ppcg init)");
    eeqMatvec(impl, N, d_q, impl.d_pcg_Ap.ptr, true);   // dense symv or the multi-GPU split
    {
        int blk = 256, grd = (N + blk - 1) / blk;
        k_pcg_init_residual<<<grd, blk, 0, s>>>(N, d_rhs_atoms, impl.d_pcg_Ap.ptr, impl.d_pcg_r.ptr);
    }
    ppcgProjectResidual(impl, N, nfrag, impl.d_pcg_r.ptr);
    ppcgPrecondition(impl, N, nfrag, impl.d_pcg_r.ptr, impl.d_pcg_z.ptr);
    checkCublasEEQ(cublasDcopy(blas, N, impl.d_pcg_z.ptr, 1, impl.d_pcg_p.ptr, 1), "Dcopy p=z (ppcg)");

    checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev (ppcg)");
    checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_z.ptr, 1,
                              impl.d_pcg_dot_scratch.ptr), "Ddot rz init (ppcg)");
    double h_rz = 0.0;
    checkCudaEEQ(cudaMemcpyAsync(&h_rz, impl.d_pcg_dot_scratch.ptr, sizeof(double),
                                  cudaMemcpyDeviceToHost, s), "D2H rz init (ppcg)");

    checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_r.ptr, 1,
                              impl.d_pcg_rnorm_scratch.ptr), "Ddot rnorm init (ppcg)");
    double h_rnorm_sq = 0.0;
    checkCudaEEQ(cudaMemcpyAsync(&h_rnorm_sq, impl.d_pcg_rnorm_scratch.ptr, sizeof(double),
                                  cudaMemcpyDeviceToHost, s), "D2H rnorm init (ppcg)");
    checkCudaEEQ(cudaStreamSynchronize(s), "sync init (ppcg)");

    bool converged = (h_rnorm_sq <= tol * tol);
    int  iters     = 0;
    // Oct 2026: mixed precision (CURCUMA_EEQ_MIXED): the iteration matvec runs on the FP32 copy of the matrix,
    // everything else in FP64. The FP32 matvec perturbs A by ~1e-7 relative, so the recursion residual drifts
    // from b - A q; a round therefore stops at max(tol, rho * |r at its start|), the true residual is recomputed
    // with the FP64 matvec and the Krylov space is restarted from it (defect correction).
    const bool mixed = impl.split_mixed && impl.split_N == N && !impl.split_parts.empty();
    const double kMixedRho = impl.mixed_rho;
    double round_tol = tol;
    if (mixed) round_tol = std::max(tol, kMixedRho * std::sqrt(h_rnorm_sq));
    int rounds = 1;
    // Oct 2026: with the local inverse, rz <= 0 means it is not positive definite: stop (the
    // caller redoes the solve with Jacobi).
    const bool li_on = impl.m_li_valid && impl.m_li_N == N && impl.m_li_nfrag == nfrag;
    if (li_on && !converged && !(h_rz > 0.0)) return false;

    // Claude Generated (Sep 2026): one host sync per iteration instead of three. alpha and beta
    // are formed on the device from device scalars (the divisions are the same IEEE operations the
    // host did); only |r|^2 and a breakdown flag come back for the convergence test. rz lives in
    // two slots used alternately, so no kernel overwrites a scalar other threads still read.
    // Breakdown rules as before: pAp <= 0 or |rz| < 1e-30 stop the iteration (flag -> alpha 0).
    if (impl.d_ppcg_scal.n < 6) impl.d_ppcg_scal.alloc(6);   // [rz0, rz1, pAp, rnorm, flag, -]
    double* d_scal = impl.d_ppcg_scal.ptr;
    cudaMemsetAsync(d_scal, 0, 6 * sizeof(double), s);
    checkCudaEEQ(cudaMemcpyAsync(d_scal, &h_rz, sizeof(double), cudaMemcpyHostToDevice, s),
                 "H2D rz (ppcg)");
    int cur = 0;   // slot holding the current rz
    const int blk_n = 256, grd_n = (N + blk_n - 1) / blk_n;
    double h_back[2] = { 0.0, 0.0 };

    for (int k = 0; k < max_iter && !converged; ++k) {
        iters = k + 1;

        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "ptr host (ppcg)");
        eeqMatvec(impl, N, impl.d_pcg_p.ptr, impl.d_pcg_Ap.ptr);   // dense symv or the multi-GPU split

        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev (ppcg)");
        checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_p.ptr, 1, impl.d_pcg_Ap.ptr, 1, d_scal + 2),
                       "Ddot pAp (ppcg)");
        k_ppcg_update_qr<<<grd_n, blk_n, 0, s>>>(N, d_scal + cur, d_scal + 2, impl.d_pcg_p.ptr,
                                                 impl.d_pcg_Ap.ptr, d_q, impl.d_pcg_r.ptr, d_scal + 4);

        ppcgProjectResidual(impl, N, nfrag, impl.d_pcg_r.ptr);

        checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_r.ptr, 1, d_scal + 3),
                       "Ddot rnorm (ppcg)");
        checkCudaEEQ(cudaMemcpyAsync(h_back, d_scal + 3, 2 * sizeof(double),
                                      cudaMemcpyDeviceToHost, s), "D2H rnorm+flag (ppcg)");
        checkCudaEEQ(cudaStreamSynchronize(s), "sync rnorm (ppcg)");
        h_rnorm_sq = h_back[0];
        if (h_back[1] != 0.0) break;   // pAp <= 0 or degenerate rz (see k_ppcg_update_*)

        if (h_rnorm_sq <= round_tol * round_tol) {
            if (!mixed) { converged = true; break; }
            // end of a round: true residual with the FP64 matvec, restart from it
            if (++rounds > 12) break;
            checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_HOST), "ptr host (ppcg restart)");
            eeqMatvec(impl, N, d_q, impl.d_pcg_Ap.ptr, true);
            k_pcg_init_residual<<<grd_n, blk_n, 0, s>>>(N, d_rhs_atoms, impl.d_pcg_Ap.ptr, impl.d_pcg_r.ptr);
            ppcgProjectResidual(impl, N, nfrag, impl.d_pcg_r.ptr);
            ppcgPrecondition(impl, N, nfrag, impl.d_pcg_r.ptr, impl.d_pcg_z.ptr);
            checkCublasEEQ(cublasDcopy(blas, N, impl.d_pcg_z.ptr, 1, impl.d_pcg_p.ptr, 1), "Dcopy p=z (ppcg restart)");
            checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev (ppcg restart)");
            checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_z.ptr, 1, d_scal + 0), "Ddot rz (ppcg restart)");
            checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_r.ptr, 1, d_scal + 3), "Ddot rnorm (ppcg restart)");
            cudaMemsetAsync(d_scal + 4, 0, sizeof(double), s);
            double h4[4] = { 0.0, 0.0, 0.0, 0.0 };
            checkCudaEEQ(cudaMemcpyAsync(h4, d_scal, 4 * sizeof(double), cudaMemcpyDeviceToHost, s), "D2H restart (ppcg)");
            checkCudaEEQ(cudaStreamSynchronize(s), "sync restart (ppcg)");
            cur = 0;
            h_rnorm_sq = h4[3];
            if (h_rnorm_sq <= tol * tol) { converged = true; impl.ppcg_true_sq = h_rnorm_sq; break; }
            if (li_on && !(h4[0] > 0.0)) break;
            round_tol = std::max(tol, kMixedRho * std::sqrt(h_rnorm_sq));
            continue;
        }

        ppcgPrecondition(impl, N, nfrag, impl.d_pcg_r.ptr, impl.d_pcg_z.ptr);
        checkCublasEEQ(cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE), "ptr dev (ppcg)");
        checkCublasEEQ(cublasDdot(blas, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_z.ptr, 1, d_scal + (1 - cur)),
                       "Ddot rz_new (ppcg)");
        if (li_on) k_li_check_rz<<<1, 1, 0, s>>>(d_scal + (1 - cur), d_scal + 4);
        k_ppcg_update_p<<<grd_n, blk_n, 0, s>>>(N, d_scal + (1 - cur), d_scal + cur,
                                                impl.d_pcg_z.ptr, impl.d_pcg_p.ptr, d_scal + 4);
        cur = 1 - cur;
    }

    impl.m_ppcg_total_calls += 1;
    impl.m_ppcg_total_iters += iters;
    // Claude Generated (Sep 2026): CURCUMA_PPCG_ITERS=1 prints the iteration count of every solve.
    static const bool print_iters = [] { const char* e = std::getenv("CURCUMA_PPCG_ITERS"); return e && *e && *e != '0'; }();
    if (print_iters)
        fprintf(stderr, "PPCG N=%d nfrag=%d iters=%d converged=%d |r|^2=%.3e rounds=%d%s\n",
                N, nfrag, iters, converged ? 1 : 0, h_rnorm_sq, mixed ? rounds : 0,
                li_on ? (impl.m_li_device ? " precond=localinv gpu build=device" : " precond=localinv gpu build=host") : "");
    if (!converged) impl.m_ppcg_nonconv_calls += 1;
    return converged;
}

// ============================================================================
// Device build of the local-inverse preconditioner (Claude Generated, Oct 2026)
// ----------------------------------------------------------------------------
// Same construction as the host EEQSolver::buildLocalInverse (maths: eeq_solver.h,
// LocalInversePreconditioner): for every atom i the dense block B = A[E_i, E_i] over the atoms
// within R_ext (ascending atom index, i included) is Cholesky-factorised, the column of i of B^-1
// is solved, and its entries for atoms within R are row i of Gl; G = (Gl + Gl^T) / 2.
// Differences to the host build:
//   * B is gathered from the device matrix d_A (k_eeq_build_matrix), so G is the local inverse of
//     exactly the matrix the CG multiplies with;
//   * neighbours by a brute-force distance scan (one thread block per atom, N^2 distances, ~0.1 ms
//     at N = 7320) instead of a cell list;
//   * no global definiteness check (the host's sparse Cholesky of G, 17 of 47 ms): every local
//     pivot must be > 0 (a cheap necessary condition), the rest is checked by the CG itself
//     (rz > 0, pAp > 0) and by a true residual after convergence (solveWithDeviceRHSAndGPUProjectedPCG);
//   * applied as z = P G r (orthogonal fragment projection) instead of G - H S^-1 H^T: no
//     S = C G C^T, whose dense factorisation costs 4.2 ms (potrf) + 9.4 ms (potri) for nfrag = 1502
//     on an RTX 5080. On the CPU this variant needs 0-2 more iterations (SP: polymer_2x 17 vs 16,
//     peo201-water500 15 vs 14, urea400-water1000 11 vs 9, peo201-chain 14 vs 14).
// Reproducible: every atom row is computed by one thread block in a fixed order, no atomics on
// floating-point values.
// ============================================================================

/// Squared distance with explicit round-to-nearest operations (no FMA contraction), so that
/// li_r2(i, j) == li_r2(j, i) bitwise and the kept pattern of G is exactly symmetric.
__device__ __forceinline__ double li_r2(const double* __restrict__ cx, const double* __restrict__ cy,
                                        const double* __restrict__ cz, int i, int j)
{
    const double dx = cx[i] - cx[j], dy = cy[i] - cy[j], dz = cz[i] - cz[j];
    return __dadd_rn(__dadd_rn(__dmul_rn(dx, dx), __dmul_rn(dy, dy)), __dmul_rn(dz, dz));
}

/// Start of column c in a packed column-major lower triangle of order n.
__device__ __forceinline__ int li_col(int c, int n) { return c * n - (c * (c - 1)) / 2; }

/// Single-precision copy of the coordinates, shifted by atom 0 (keeps the magnitudes at the size of
/// the system), for a cheap distance PRE-filter; w = 0. maxabs[0] receives max |shifted coordinate|
/// as float bits (non-negative floats order like ints).
__global__ void k_li_pack_f(int N, const double* __restrict__ cx, const double* __restrict__ cy,
                            const double* __restrict__ cz, float4* __restrict__ xf, unsigned int* maxabs)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= N) return;
    const float x = static_cast<float>(cx[j] - cx[0]);
    const float y = static_cast<float>(cy[j] - cy[0]);
    const float z = static_cast<float>(cz[j] - cz[0]);
    xf[j] = make_float4(x, y, z, 0.0f);
    const float m = fmaxf(fabsf(x), fmaxf(fabsf(y), fabsf(z)));
    atomicMax(maxabs, __float_as_uint(m));
}

/// Squared-distance bound for the single-precision prefilter: every pair with exact r^2 <= E2 has
/// float r^2 <= this. Error per coordinate difference <= 2 eps_f (|x_i| + |x_j|) + rounding of the
/// shift, taken generously as 8 eps_f maxabs; the decision itself is always made with li_r2.
__device__ __forceinline__ float li_prefilter_bound(double E2, const unsigned int* maxabs)
{
    const double delta = 8.0 * 1.2e-7 * static_cast<double>(__uint_as_float(*maxabs)) + 1e-6;
    const double r = sqrt(E2) + 2.0 * sqrt(3.0) * delta;
    return static_cast<float>(r * r * (1.0 + 1e-5));
}

/// Per atom: block size n_i (atoms within R_ext, i included) and kept entries k_i (within R).
__global__ void k_li_count(int N, const double* __restrict__ cx, const double* __restrict__ cy,
                           const double* __restrict__ cz, const float4* __restrict__ xf,
                           const unsigned int* __restrict__ maxabs, double E2, double R2,
                           int* __restrict__ nblk, int* __restrict__ nkeep)
{
    __shared__ int s_n, s_k;
    const float E2f = li_prefilter_bound(E2, maxabs);
    for (int i = blockIdx.x; i < N; i += gridDim.x) {
        if (threadIdx.x == 0) { s_n = 0; s_k = 0; }
        __syncthreads();
        const float4 pi = xf[i];
        int n = 0, k = 0;
        for (int j = threadIdx.x; j < N; j += blockDim.x) {
            const float4 pj = xf[j];
            const float fx = pi.x - pj.x, fy = pi.y - pj.y, fz = pi.z - pj.z;
            if (fx * fx + fy * fy + fz * fz > E2f) continue;   // certainly outside R_ext
            const double r2 = li_r2(cx, cy, cz, i, j);
            if (r2 <= E2) { ++n; if (r2 <= R2) ++k; }
        }
        atomicAdd(&s_n, n);   // integer sums: order-independent
        atomicAdd(&s_k, k);
        __syncthreads();
        if (threadIdx.x == 0) { nblk[i] = s_n; nkeep[i] = s_k; }
        __syncthreads();
    }
}

/// Row i of Gl for the atoms atoms[0..count) (one thread block per atom, grid-stride). The atoms
/// are launched in buckets of similar block size n_i, so that the dynamic shared memory (packed
/// lower triangle + solved column + atom list) is sized for the bucket and not for the largest
/// block of the system; a bucket whose largest block does not fit uses a per-block slice of global
/// scratch instead (use_global). fail: 1 = a pivot <= 0 (block not positive definite), 2 =
/// neighbour or kept count differs from k_li_count.
__global__ void k_li_blocks(int N, const int* __restrict__ atoms, int count,
                            const double* __restrict__ cx, const double* __restrict__ cy,
                            const double* __restrict__ cz, const float4* __restrict__ xf,
                            const unsigned int* __restrict__ maxabs, double E2, double R2,
                            const double* __restrict__ A, int lda,
                            const int* __restrict__ nblk, const int* __restrict__ rowptr, int n_max,
                            int* __restrict__ Gcol, double* __restrict__ Glval, int* fail,
                            int use_global, double* gscratch, long long gstride)
{
    extern __shared__ double li_smem[];
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    const int bd = blockDim.x;
    const long long npk = static_cast<long long>(n_max) * (n_max + 1) / 2;
    double* L;
    int* idx;
    if (use_global) {
        L = gscratch + static_cast<long long>(blockIdx.x) * gstride;
        idx = reinterpret_cast<int*>(li_smem);
    } else {
        L = li_smem;
        idx = reinterpret_cast<int*>(li_smem + npk + n_max);
    }
    double* x = L + npk;
    unsigned char* keep = reinterpret_cast<unsigned char*>(idx + n_max);
    int* wlist = reinterpret_cast<int*>(keep + ((n_max + 15) / 16) * 16);   // nwarps * n_max: (j << 1) | kept

    __shared__ int s_k, s_fail;
    __shared__ int s_wcnt[32];
    const float E2f = li_prefilter_bound(E2, maxabs);

    for (int b = blockIdx.x; b < count; b += gridDim.x) {
        const int i = atoms[b];
        if (tid == 0) { s_k = -1; s_fail = 0; }
        const float4 pi = xf[i];
        // (a) neighbours within R_ext in ascending atom order: warp w scans the w-th contiguous
        //     segment of atom indices (ballots keep the order inside a warp), then the per-warp
        //     lists are concatenated in warp order.
        const int seg0 = static_cast<int>(static_cast<long long>(N) * warp / nwarps);
        const int seg1 = static_cast<int>(static_cast<long long>(N) * (warp + 1) / nwarps);
        int* wl = wlist + warp * n_max;
        int cnt = 0;
        for (int j0 = seg0; j0 < seg1; j0 += 32) {
            const int j = j0 + lane;
            bool in = false, kept = false;
            if (j < seg1) {
                const float4 pj = xf[j];
                const float fx = pi.x - pj.x, fy = pi.y - pj.y, fz = pi.z - pj.z;
                if (fx * fx + fy * fy + fz * fz <= E2f) {
                    const double r2 = li_r2(cx, cy, cz, i, j);
                    in = (r2 <= E2);
                    kept = (r2 <= R2);
                }
            }
            const unsigned m = __ballot_sync(0xffffffffu, in);
            if (in) {
                const int p = cnt + __popc(m & ((1u << lane) - 1u));
                if (p < n_max) wl[p] = (j << 1) | (kept ? 1 : 0);
            }
            cnt += __popc(m);
        }
        if (lane == 0) s_wcnt[warp] = cnt;
        __syncthreads();
        int off = 0, n = 0;
        for (int w = 0; w < nwarps; ++w) { if (w < warp) off += s_wcnt[w]; n += s_wcnt[w]; }
        if (n <= n_max) {
            for (int p = lane; p < cnt; p += 32) {
                const int v = wl[p], j = v >> 1;
                idx[off + p] = j;
                keep[off + p] = static_cast<unsigned char>(v & 1);
                if (j == i) s_k = off + p;
            }
        }
        __syncthreads();
        const int k = s_k;
        if (n > n_max || n != nblk[i] || k < 0) {
            if (tid == 0) atomicExch(fail, 2);
            return;
        }
        // (b) packed lower triangle of B = A[idx, idx] from the device matrix (column-major, lda)
        for (int c = warp; c < n; c += nwarps) {
            const double* Acol = A + static_cast<long long>(idx[c]) * lda;
            double* Lc = L + li_col(c, n);
            for (int r = c + lane; r < n; r += 32) Lc[r - c] = Acol[idx[r]];
        }
        __syncthreads();
        // (c) right-looking Cholesky B = L L^T (same operation order per element as the host
        //     localCholeskyLower). Pivot <= 0: the block is not positive definite.
        for (int kk = 0; kk < n; ++kk) {
            double* Lk = L + li_col(kk, n);
            if (tid == 0) {
                const double p = Lk[0];
                if (!(p > 0.0)) s_fail = 1; else Lk[0] = sqrt(p);
            }
            __syncthreads();
            if (s_fail) break;
            const double inv = 1.0 / Lk[0];
            for (int r = kk + 1 + tid; r < n; r += bd) Lk[r - kk] *= inv;
            __syncthreads();
            for (int c = kk + 1 + warp; c < n; c += nwarps) {
                const double l = Lk[c - kk];
                double* Lc = L + li_col(c, n);
                for (int r = c + lane; r < n; r += 32) Lc[r - c] -= l * Lk[r - kk];
            }
            __syncthreads();
        }
        if (s_fail) {
            if (tid == 0) atomicExch(fail, 1);
            return;
        }
        // (d) L y = e_k (y_m = 0 for m < k)
        for (int m = tid; m < n; m += bd) x[m] = (m == k) ? 1.0 : 0.0;
        __syncthreads();
        for (int c = k; c < n; ++c) {
            const double* Lc = L + li_col(c, n);
            if (tid == 0) x[c] /= Lc[0];
            __syncthreads();
            const double xc = x[c];
            for (int r = c + 1 + tid; r < n; r += bd) x[r] -= Lc[r - c] * xc;
            __syncthreads();
        }
        // (e) L^T x = y, column-oriented: x_c final, then x_r -= L(c, r) x_c for r < c
        for (int c = n - 1; c >= 0; --c) {
            if (tid == 0) x[c] /= L[li_col(c, n)];
            __syncthreads();
            const double xc = x[c];
            for (int r = tid; r < c; r += bd) x[r] -= L[li_col(r, n) + (c - r)] * xc;
            __syncthreads();
        }
        // (f) kept entries of row i, ascending columns
        if (tid == 0) {
            int p = rowptr[i];
            const int pend = rowptr[i + 1];
            for (int m = 0; m < n; ++m) {
                if (!keep[m]) continue;
                if (p < pend) { Gcol[p] = idx[m]; Glval[p] = x[m]; }
                ++p;
            }
            if (p != pend) atomicExch(fail, 2);
        }
        __syncthreads();
    }
}

/// G_ij = (Gl_ij + Gl_ji) / 2 on the (symmetric) pattern; Gl_ji by binary search in row j.
__global__ void k_li_symmetrize(int N, const int* __restrict__ rowptr, const int* __restrict__ col,
                                const double* __restrict__ Glval, double* __restrict__ Gval, int* fail)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    for (int e = rowptr[i]; e < rowptr[i + 1]; ++e) {
        const int j = col[e];
        int lo = rowptr[j], hi = rowptr[j + 1] - 1, pos = -1;
        while (lo <= hi) {
            const int mid = (lo + hi) >> 1;
            const int c = col[mid];
            if (c == i) { pos = mid; break; }
            if (c < i) lo = mid + 1; else hi = mid - 1;
        }
        if (pos < 0) { atomicExch(fail, 3); return; }
        Gval[e] = 0.5 * Glval[e] + 0.5 * Glval[pos];
    }
}

/// Test hook CURCUMA_LOCALINV_TEST_INDEFINITE: 1 = G -> -G (negative definite), 2 = G_ii -> -G_ii
/// for every 10th atom (indefinite), 3 = G unchanged but the true-residual check fails every time.
/// Only for exercising the fallbacks; 0 / unset = off.
static int liTestMode()
{
    static const int mode = [] { const char* e = std::getenv("CURCUMA_LOCALINV_TEST_INDEFINITE"); return e ? std::atoi(e) : 0; }();
    return mode;
}

/// Kernel of the test hook (modes 1 and 2).
__global__ void k_li_test_perturb(int N, const int* __restrict__ rowptr, const int* __restrict__ col,
                                  double* __restrict__ Gval, int mode)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    for (int e = rowptr[i]; e < rowptr[i + 1]; ++e)
        if (mode == 1 || (mode == 2 && col[e] == i && i % 10 == 0)) Gval[e] = -Gval[e];
}

/// Builds G on the device from impl.d_A (must hold the current EEQ matrix, full, column-major).
/// On success the preconditioner is active in device mode (m_li_device). Returns false with a
/// reason; the caller then runs Jacobi. Claude Generated (Oct 2026).
static bool buildLocalInverseOnDevice(EEQSolverGPUImpl& impl, int N, int nfrag, const double* cx,
                                      const double* cy, const double* cz, std::string& why)
{
    cudaStream_t s = impl.stream;
    // The matrix assembly queued before is not part of the build time: wait for it first (the
    // count download below synchronises anyway, so this costs nothing extra).
    checkCudaEEQ(cudaStreamSynchronize(s), "sync before li build");
    const auto t0 = std::chrono::steady_clock::now();
    impl.m_li_valid = false;
    const double R2 = impl.m_li_radius * impl.m_li_radius;
    const double E2 = impl.m_li_ext * impl.m_li_ext;
    if (impl.d_li_nblk.n < N) impl.d_li_nblk.alloc(N);
    if (impl.d_li_nkeep.n < N) impl.d_li_nkeep.alloc(N);
    if (impl.d_li_fail.n < 1) impl.d_li_fail.alloc(1);
    if (impl.d_li_G_rowptr.n < N + 1) impl.d_li_G_rowptr.alloc(N + 1);
    if (impl.d_li_xf.n < N) impl.d_li_xf.alloc(N);
    if (impl.d_li_maxabs.n < 1) impl.d_li_maxabs.alloc(1);
    cudaMemsetAsync(impl.d_li_maxabs.ptr, 0, sizeof(unsigned int), s);
    k_li_pack_f<<<(N + 255) / 256, 256, 0, s>>>(N, cx, cy, cz, impl.d_li_xf.ptr, impl.d_li_maxabs.ptr);
    constexpr int kThreads = 128;
    k_li_count<<<N, kThreads, 0, s>>>(N, cx, cy, cz, impl.d_li_xf.ptr, impl.d_li_maxabs.ptr, E2, R2,
                                      impl.d_li_nblk.ptr, impl.d_li_nkeep.ptr);
    std::vector<int> nblk(N), nkeep(N);
    checkCudaEEQ(cudaMemcpyAsync(nblk.data(), impl.d_li_nblk.ptr, N * sizeof(int), cudaMemcpyDeviceToHost, s), "D2H li nblk");
    checkCudaEEQ(cudaMemcpyAsync(nkeep.data(), impl.d_li_nkeep.ptr, N * sizeof(int), cudaMemcpyDeviceToHost, s), "D2H li nkeep");
    checkCudaEEQ(cudaStreamSynchronize(s), "sync li count");
    const double count_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

    std::vector<int> rowptr(N + 1, 0);
    int n_max = 0;
    double sum_block = 0.0;
    for (int i = 0; i < N; ++i) {
        rowptr[i + 1] = rowptr[i] + nkeep[i];
        n_max = std::max(n_max, nblk[i]);
        sum_block += nblk[i];
    }
    const int nnz = rowptr[N];
    if (impl.d_li_G_col.n < nnz) impl.d_li_G_col.alloc(nnz);
    if (impl.d_li_G_val.n < nnz) impl.d_li_G_val.alloc(nnz);
    if (impl.d_li_Gl_val.n < nnz) impl.d_li_Gl_val.alloc(nnz);
    checkCudaEEQ(cudaMemcpyAsync(impl.d_li_G_rowptr.ptr, rowptr.data(), (N + 1) * sizeof(int),
                                 cudaMemcpyHostToDevice, s), "H2D li rowptr");
    cudaMemsetAsync(impl.d_li_fail.ptr, 0, sizeof(int), s);

    // Shared memory: packed lower triangle + solved column (doubles), atom list (int), kept flags.
    if (impl.m_li_max_smem < 0) {
        int dev = 0, v = 0;
        cudaGetDevice(&dev);
        cudaDeviceGetAttribute(&v, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev);
        impl.m_li_max_smem = v;
    }
    // Buckets of block size (<= 32, <= 64, <= 128, larger), each launch with the shared memory of its
    // largest block: most blocks are small (mean 31 atoms on polymer_2x, largest 86), and memory for
    // the largest block everywhere would leave few blocks resident per multiprocessor.
    const int bucket_limit[4] = { 32, 64, 128, std::max(n_max, 129) };
    std::vector<int> order;
    order.reserve(N);
    int bucket_begin[5] = { 0, 0, 0, 0, 0 };
    int bucket_nmax[4] = { 0, 0, 0, 0 };
    for (int bkt = 0; bkt < 4; ++bkt) {
        bucket_begin[bkt] = static_cast<int>(order.size());
        const int lo = bkt == 0 ? 0 : bucket_limit[bkt - 1];
        for (int i = 0; i < N; ++i)
            if (nblk[i] > lo && nblk[i] <= bucket_limit[bkt]) { order.push_back(i); bucket_nmax[bkt] = std::max(bucket_nmax[bkt], nblk[i]); }
    }
    bucket_begin[4] = static_cast<int>(order.size());
    if (impl.d_li_atoms.n < N) impl.d_li_atoms.alloc(N);
    checkCudaEEQ(cudaMemcpyAsync(impl.d_li_atoms.ptr, order.data(), order.size() * sizeof(int),
                                 cudaMemcpyHostToDevice, s), "H2D li atom order");
    const int nwarps = kThreads / 32;
    bool any_global = false;
    for (int bkt = 0; bkt < 4; ++bkt) {
        const int count = bucket_begin[bkt + 1] - bucket_begin[bkt];
        if (count == 0) continue;
        const int nb = bucket_nmax[bkt];
        const long long npk = static_cast<long long>(nb) * (nb + 1) / 2;
        const size_t small_bytes = static_cast<size_t>(nb) * (sizeof(int) + 1) + 16
                                 + static_cast<size_t>(nwarps) * nb * sizeof(int) + 16;
        const size_t full_bytes = static_cast<size_t>(npk + nb) * sizeof(double) + small_bytes;
        const size_t static_reserve = 1024;   // s_k, s_fail, s_wcnt and alignment
        const bool use_global = full_bytes + static_reserve > static_cast<size_t>(impl.m_li_max_smem);
        const size_t smem = use_global ? small_bytes : full_bytes;
        int grid = count;
        long long gstride = 0;
        if (use_global) {
            any_global = true;
            gstride = ((npk + nb + 31) / 32) * 32;
            const long long budget = 256LL * 1024 * 1024 / sizeof(double);   // 256 MB of scratch at most
            grid = static_cast<int>(std::max(1LL, std::min<long long>(count, budget / gstride)));
            if (impl.d_li_scratch.n < grid * gstride) impl.d_li_scratch.alloc(static_cast<int>(grid * gstride));
        }
        if (smem > 48 * 1024)
            checkCudaEEQ(cudaFuncSetAttribute(k_li_blocks, cudaFuncAttributeMaxDynamicSharedMemorySize,
                                              static_cast<int>(smem)), "li smem attribute");
        k_li_blocks<<<grid, kThreads, smem, s>>>(N, impl.d_li_atoms.ptr + bucket_begin[bkt], count, cx, cy, cz,
                                                 impl.d_li_xf.ptr, impl.d_li_maxabs.ptr, E2, R2,
                                                 impl.d_A.ptr, N, impl.d_li_nblk.ptr,
                                                 impl.d_li_G_rowptr.ptr, nb, impl.d_li_G_col.ptr,
                                                 impl.d_li_Gl_val.ptr, impl.d_li_fail.ptr, use_global ? 1 : 0,
                                                 use_global ? impl.d_li_scratch.ptr : nullptr, gstride);
    }
    checkCudaEEQ(cudaGetLastError(), "launch k_li_blocks");
    {
        const int blk = 256, grd = (N + blk - 1) / blk;
        k_li_symmetrize<<<grd, blk, 0, s>>>(N, impl.d_li_G_rowptr.ptr, impl.d_li_G_col.ptr,
                                            impl.d_li_Gl_val.ptr, impl.d_li_G_val.ptr, impl.d_li_fail.ptr);
        const int test_mode = liTestMode();
        if (test_mode == 1 || test_mode == 2)
            k_li_test_perturb<<<grd, blk, 0, s>>>(N, impl.d_li_G_rowptr.ptr, impl.d_li_G_col.ptr,
                                                  impl.d_li_G_val.ptr, test_mode);
    }
    int h_fail = 0;
    checkCudaEEQ(cudaMemcpyAsync(&h_fail, impl.d_li_fail.ptr, sizeof(int), cudaMemcpyDeviceToHost, s), "D2H li fail");
    checkCudaEEQ(cudaStreamSynchronize(s), "sync li build");
    const double build_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    impl.m_li_stats.last_build_ms = build_ms;
    impl.m_li_stats.total_build_ms += build_ms;

    // Test hook shared with the host build: CURCUMA_LOCALINV_TEST_REJECT=1 rejects every build.
    static const bool test_reject = [] { const char* e = std::getenv("CURCUMA_LOCALINV_TEST_REJECT"); return e && *e && *e != '0'; }();
    if (h_fail == 1) why = "a dense block A[E,E] is not positive definite (local Cholesky pivot <= 0)";
    else if (h_fail == 2) why = "neighbour list of the device build inconsistent";
    else if (h_fail == 3) why = "pattern of the device-built G is not symmetric";
    else if (test_reject) why = "rejected by CURCUMA_LOCALINV_TEST_REJECT";
    const bool ok = (h_fail == 0) && !test_reject;

    static const bool print_iters = [] { const char* e = std::getenv("CURCUMA_PPCG_ITERS"); return e && *e && *e != '0'; }();
    if (print_iters)
        fprintf(stderr, "LOCALINV GPU device build ok=%d build_ms=%.2f count_ms=%.2f nnz_row=%.1f block_mean=%.1f block_max=%d mem=%s\n",
                ok ? 1 : 0, build_ms, count_ms, static_cast<double>(nnz) / N, sum_block / N, n_max,
                any_global ? "global" : "shared");
    if (!ok) return false;
    impl.m_li_N = N;
    impl.m_li_nfrag = nfrag;
    impl.m_li_device = true;
    impl.m_li_valid = true;
    return true;
}

bool EEQSolverGPU::solveWithDeviceRHSAndGPUProjectedPCG(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha_corrected,
    const double* d_gam_corrected,
    const double* d_rhs_atoms,
    const std::vector<int>& /*fraglist*/,
    const std::vector<double>& rhs_constraints,
    int    max_iter,
    double tol,
    double cutoff_sq,
    bool   force_refactor)
{
    auto& impl = *m_impl;
    if (nfrag < 1) return false;
    if (!impl.m_frag_topo_valid) return false;
    if (impl.m_nfrag_batched != nfrag) return false;
    if (impl.d_atom_frag.n < natoms) return false;
    if (impl.d_frag_inv_atoms.n < nfrag) return false;

    const int N = natoms;

    // Claude Generated (Sep 2026, MULTI_GPU_GAPS F-12): d_rhs_atoms and the coordinates are
    // written on the FFWorkspaceGPU main stream, which this solver does not own. The other
    // device-RHS paths start with a blocking cudaMemcpy, which orders them implicitly through
    // the legacy default stream; this path had no such operation since WP5-D removed the sync.
    // An event recorded on the legacy stream waits for all prior work in every blocking stream
    // (all curcuma streams are created with cudaStreamCreate), so impl.stream waits for it.
    // K3 (Oct 2026, Claude Generated): with an input-ready event from the workspace, wait only for
    // the main-stream point after k_build_eeq_rhs; dc6dcn and the Phase-1 worker streams run in
    // parallel to the matrix build. Without it, the legacy barrier as before.
    impl.orderAfterWorkspace(true);

    const bool do_refactor = force_refactor
                          || !impl.m_pcg_M_inv_valid
                          || (N != m_last_N);

    // Claude Generated (Oct 2026): CURCUMA_PPCG_GPU_TIMING=1 synchronises after each phase and prints
    // matrix / local-inverse build / CG / true-residual wall times (measurement hook, adds syncs).
    static const bool phase_timing = [] { const char* e = std::getenv("CURCUMA_PPCG_GPU_TIMING"); return e && *e && *e != '0'; }();
    auto phase_clock = [&]() {
        if (phase_timing) cudaStreamSynchronize(impl.stream);
        return std::chrono::steady_clock::now();
    };
    const auto tp0 = phase_clock();

    // Claude Generated (Sep 2026, F-3): multi-GPU split when configured and N is large enough.
    // A split that fails to set up (status set) falls back to the single-device path for good.
    const bool use_split = splitWanted(impl, N) && splitPrepare(impl, N);
    if (!use_split && impl.split_N != 0) impl.releaseSplit();

    if (do_refactor || (use_split && impl.split_builds == 0)) {
        if (impl.d_pcg_M_inv.n < N) impl.d_pcg_M_inv.alloc(N);
        if (use_split) {
            splitBuild(impl, N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, cutoff_sq);
            int blk = 256, grd = (N + blk - 1) / blk;
            k_eeq_diag_inv<<<grd, blk, 0, impl.stream>>>(N, d_alpha_corrected, d_gam_corrected,
                                                          impl.d_pcg_M_inv.ptr);
        } else {
            {
                int n_lower = N * (N + 1) / 2;
                int block   = 256;
                int grid    = (n_lower + block - 1) / block;
                k_eeq_build_matrix<<<grid, block, 0, impl.stream>>>(
                    N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, impl.d_A.ptr, cutoff_sq);
            }
            int blk = 256, grd = (N + blk - 1) / blk;
            k_pcg_extract_diag_inv<<<grd, blk, 0, impl.stream>>>(
                N, impl.d_A.ptr, impl.d_pcg_M_inv.ptr);
        }
        impl.m_pcg_M_inv_valid = true;

        // Fragment-weighted diagonal sum sM[f] = Σ_{i∈f} M_inv[i] — needed by the
        // projected preconditioner; rebuilt whenever M_inv itself is rebuilt.
        cudaMemsetAsync(impl.d_frag_sM.ptr, 0, (size_t)nfrag * sizeof(double), impl.stream);
        {
            fragSumDet(impl, nfrag, impl.d_pcg_M_inv.ptr, impl.d_frag_sM.ptr, impl.stream);
        }
        // No block-Jacobi here (unlike WP7-C): the CPU algorithm this ports uses only the
        // diagonal Jacobi. Measured Sep 29, 2026: an exact block inverse for every small fragment
        // (1500 waters of polymer_2x) left the iteration count unchanged (73/66/68 vs 75/68/68) -
        // the inter-fragment Coulomb coupling limits convergence, not the intra-fragment one.
    }

    const auto tp1 = phase_clock();
    // Oct 2026: local inverse built on the device from the matrix just assembled (requested by
    // the wrapper via requestLocalInverseDeviceBuild when its rebuild criterion fired).
    if (impl.m_li_build_pending) {
        impl.m_li_build_pending = false;
        std::string why;
        bool built = false;
        if (use_split) {
            // Oct 2026: the blocks need the full dense matrix. The split keeps only column blocks, so
            // assemble the dense matrix once on the own device for this rebuild (about every 8 MD steps);
            // a host build of the local inverse cost ~70 ms per rebuild on polymer_2x.
            impl.m_li_valid = false;
            try {
                impl.d_A.ensure(N * N);
                const int n_lower = N * (N + 1) / 2;
                k_eeq_build_matrix<<<(n_lower + 255) / 256, 256, 0, impl.stream>>>(
                    N, cx, cy, cz, d_alpha_corrected, d_gam_corrected, impl.d_A.ptr, cutoff_sq);
                built = buildLocalInverseOnDevice(impl, N, nfrag, cx, cy, cz, why);
            } catch (const std::exception& e) {
                why = std::string("dense matrix for the split rebuild: ") + e.what();
            }
        } else {
            built = buildLocalInverseOnDevice(impl, N, nfrag, cx, cy, cz, why);
        }
        if (built) {
            impl.m_li_stats.builds += 1;
        } else {
            impl.m_li_stats.build_failures += 1;
            impl.m_li_stats.last_why = why;
        }
    }

    if (impl.d_rhs.n < N) impl.d_rhs.alloc(N);

    // Warm-start q ← persistent cache (D2D), else zero. (A linear extrapolation from the last two
    // solutions was tried Sep 29, 2026: polymer_2x MD, 20 steps, 1-2 of ~68 iterations saved -
    // the iteration count is set by the conditioning, not by the start vector. Not kept.)
    if (impl.m_ppcg_warm_valid && impl.d_ppcg_q_persistent.n >= N) {
        checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, impl.d_ppcg_q_persistent.ptr,
                                      N * sizeof(double), cudaMemcpyDeviceToDevice, impl.stream),
                     "D2D ppcg warm-start");
    } else {
        cudaMemsetAsync(impl.d_rhs.ptr, 0, N * sizeof(double), impl.stream);
    }

    // Upload the (small, nfrag-sized) per-call constraint RHS.
    {
        std::vector<double> h_rc(nfrag, 0.0);
        for (int f = 0; f < nfrag && f < (int)rhs_constraints.size(); ++f)
            h_rc[f] = rhs_constraints[f];
        checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs_constraints_dev.ptr, h_rc.data(),
                                      nfrag * sizeof(double), cudaMemcpyHostToDevice, impl.stream),
                     "H2D rhs_constraints (ppcg)");
    }

    // Oct 2026: with the local-inverse preconditioner keep the start vector, so a breakdown
    // (rz <= 0, pAp <= 0) or a stall can be redone with the diagonal Jacobi from the same start.
    const bool li_on = impl.m_li_valid && impl.m_li_N == N && impl.m_li_nfrag == nfrag;
    if (li_on) {
        if (impl.d_li_q0.n < N) impl.d_li_q0.alloc(N);
        checkCudaEEQ(cudaMemcpyAsync(impl.d_li_q0.ptr, impl.d_rhs.ptr, N * sizeof(double),
                                      cudaMemcpyDeviceToDevice, impl.stream), "D2D localinv q0");
    }
    const auto tp2 = phase_clock();
    bool ok = runProjectedPCG(impl, N, nfrag, d_rhs_atoms, impl.d_rhs.ptr, max_iter, tol);
    const auto tp3 = phase_clock();
    bool residual_failed = false;
    if (ok && li_on && impl.m_li_device) {
        // Acceptance of a device-built G (Claude Generated, Oct 2026): the device build has no
        // global definiteness check, and an indefinite M can let the RECURSIVE residual of the CG
        // drift away from the true one. One extra matvec: |P(b - A q)| must be <= 2 tol (the
        // recursion gap from rounding is far below tol for a positive definite M).
        constexpr double kTrueResidualFactor = 2.0;
        double h_true_sq = impl.ppcg_true_sq;   // mixed solve: the last restart measured exactly this quantity
        if (h_true_sq < 0.0) {
            eeqMatvec(impl, N, impl.d_rhs.ptr, impl.d_pcg_Ap.ptr, true);
            {
                const int blk = 256, grd = (N + blk - 1) / blk;
                k_pcg_init_residual<<<grd, blk, 0, impl.stream>>>(N, d_rhs_atoms, impl.d_pcg_Ap.ptr, impl.d_pcg_r.ptr);
            }
            ppcgProjectResidual(impl, N, nfrag, impl.d_pcg_r.ptr);
            checkCublasEEQ(cublasSetPointerMode(impl.cublas_handle, CUBLAS_POINTER_MODE_HOST), "ptr host (li true res)");
            checkCublasEEQ(cublasDdot(impl.cublas_handle, N, impl.d_pcg_r.ptr, 1, impl.d_pcg_r.ptr, 1, &h_true_sq),
                           "Ddot true residual (li)");
        }
        const double ratio = (tol > 0.0) ? std::sqrt(h_true_sq) / tol : 0.0;
        static const bool print_iters = [] { const char* e = std::getenv("CURCUMA_PPCG_ITERS"); return e && *e && *e != '0'; }();
        if (print_iters) fprintf(stderr, "LOCALINV GPU true residual |P(b-Aq)|/tol=%.3f\n", ratio);
        if (liTestMode() == 3) h_true_sq = std::numeric_limits<double>::infinity();   // test hook
        if (!(h_true_sq <= kTrueResidualFactor * kTrueResidualFactor * tol * tol)) {
            ok = false;
            residual_failed = true;
            impl.m_li_stats.residual_fallbacks += 1;
        } else {
            impl.m_li_stats.last_true_residual_ratio = ratio;
        }
    }
    if (phase_timing) {
        const auto tp4 = phase_clock();
        auto ms = [](std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a).count(); };
        fprintf(stderr, "PPCG GPU timing matrix_ms=%.3f li_build_ms=%.3f cg_ms=%.3f true_res_ms=%.3f\n",
                ms(tp0, tp1), ms(tp1, tp2), ms(tp2, tp3), ms(tp3, tp4));
    }
    if (!ok && li_on) {
        impl.m_li_valid = false;   // never continue with this M; the host rebuilds or stays on Jacobi
        impl.m_li_runtime_fallbacks += 1;
        impl.m_li_stats.runtime_fallbacks += 1;
        impl.m_li_stats.last_why = residual_failed ? "true residual above 2 tol" : "breakdown (rz <= 0 / pAp <= 0) or stall";
        CurcumaLogger::warn(fmt::format("EEQ GPU ppcg localinv: {}, solve repeated with diagonal "
                                        "Jacobi (fallbacks so far: {})", impl.m_li_stats.last_why,
                                        impl.m_li_runtime_fallbacks));
        checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, impl.d_li_q0.ptr, N * sizeof(double),
                                      cudaMemcpyDeviceToDevice, impl.stream), "D2D localinv q0 restore");
        ok = runProjectedPCG(impl, N, nfrag, d_rhs_atoms, impl.d_rhs.ptr, max_iter, tol);
    }
    if (!ok) return false;

    // Persist warm-start (charges are already at d_rhs[0..N-1] — the shared
    // getDeviceChargesPtr() contract — so just cache a copy for next call).
    if (impl.d_ppcg_q_persistent.n < N) impl.d_ppcg_q_persistent.alloc(N);
    checkCudaEEQ(cudaMemcpyAsync(impl.d_ppcg_q_persistent.ptr, impl.d_rhs.ptr,
                                  N * sizeof(double), cudaMemcpyDeviceToDevice, impl.stream),
                 "D2D persist ppcg q");
    impl.m_ppcg_warm_valid = true;

    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after ppcg");
    return true;
}

// ============================================================================
// Local-inverse preconditioner upload (Claude Generated, Oct 2026)
// ============================================================================

bool EEQSolverGPU::setLocalInversePreconditioner(int N, int nfrag,
    const std::vector<int>& G_rowptr, const std::vector<int>& G_col, const std::vector<double>& G_val,
    const std::vector<int>& H_rowptr, const std::vector<int>& H_col, const std::vector<double>& H_val,
    const std::vector<double>& S_inv)
{
    auto& impl = *m_impl;
    impl.m_li_valid = false;
    if (static_cast<int>(G_rowptr.size()) != N + 1 || static_cast<int>(H_rowptr.size()) != N + 1
        || static_cast<long long>(S_inv.size()) != static_cast<long long>(nfrag) * nfrag)
        return false;
    auto up_i = [&](CudaBuffer<int>& b, const std::vector<int>& v) {
        if (b.n < static_cast<int>(v.size())) b.alloc(static_cast<int>(v.size()));
        return cudaMemcpyAsync(b.ptr, v.data(), v.size() * sizeof(int), cudaMemcpyHostToDevice, impl.stream) == cudaSuccess;
    };
    auto up_d = [&](CudaBuffer<double>& b, const std::vector<double>& v) {
        if (b.n < static_cast<int>(v.size())) b.alloc(static_cast<int>(v.size()));
        return cudaMemcpyAsync(b.ptr, v.data(), v.size() * sizeof(double), cudaMemcpyHostToDevice, impl.stream) == cudaSuccess;
    };
    bool ok = up_i(impl.d_li_G_rowptr, G_rowptr) && up_i(impl.d_li_G_col, G_col) && up_d(impl.d_li_G_val, G_val)
           && up_i(impl.d_li_H_rowptr, H_rowptr) && up_i(impl.d_li_H_col, H_col) && up_d(impl.d_li_H_val, H_val)
           && up_d(impl.d_li_S_inv, S_inv);
    if (impl.d_li_t.n < nfrag) impl.d_li_t.alloc(nfrag);
    if (impl.d_li_mu.n < nfrag) impl.d_li_mu.alloc(nfrag);
    ok = ok && cudaStreamSynchronize(impl.stream) == cudaSuccess;
    if (!ok) return false;
    impl.m_li_N = N;
    impl.m_li_nfrag = nfrag;
    impl.m_li_device = false;
    impl.m_li_valid = true;
    return true;
}

void EEQSolverGPU::clearLocalInversePreconditioner() { m_impl->m_li_valid = false; }
bool EEQSolverGPU::localInversePreconditionerActive() const { return m_impl->m_li_valid; }
long EEQSolverGPU::localInverseRuntimeFallbacks() const { return m_impl->m_li_runtime_fallbacks; }

void EEQSolverGPU::requestLocalInverseDeviceBuild(double radius_bohr, double ext_bohr)
{
    auto& impl = *m_impl;
    impl.m_li_radius = radius_bohr;
    impl.m_li_ext = std::max(ext_bohr, radius_bohr);
    impl.m_li_build_pending = true;
    impl.m_li_valid = false;   // the previous G must not be used with the new geometry's request pending
}

// K3 (Oct 2026, Claude Generated): see header.
void EEQSolverGPU::setInputReadyEvent(void* cuda_event)
{
    m_impl->ev_input_ready = static_cast<cudaEvent_t>(cuda_event);
}

// K2/K3 (Oct 2026, Claude Generated): see header. Recreates the solver stream with the requested
// priority and rebinds the cuSOLVER/cuBLAS handles; numerics are unchanged (same kernels, order).
bool EEQSolverGPU::setStreamPriorityHigh(bool high)
{
    auto& impl = *m_impl;
    if (impl.stream_high_priority == high) return true;
    int least = 0, greatest = 0;
    if (cudaDeviceGetStreamPriorityRange(&least, &greatest) != cudaSuccess) return false;
    cudaStream_t s = nullptr;
    if (cudaStreamCreateWithPriority(&s, cudaStreamDefault, high ? greatest : least) != cudaSuccess)
        return false;
    cudaStreamSynchronize(impl.stream);
    if (cusolverDnSetStream(impl.cusolver_handle, s) != CUSOLVER_STATUS_SUCCESS
        || cublasSetStream(impl.cublas_handle, s) != CUBLAS_STATUS_SUCCESS) {
        cusolverDnSetStream(impl.cusolver_handle, impl.stream);
        cublasSetStream(impl.cublas_handle, impl.stream);
        cudaStreamDestroy(s);
        return false;
    }
    cudaStreamDestroy(impl.stream);
    impl.stream = s;
    impl.stream_high_priority = high;
    return true;
}

bool EEQSolverGPU::localInverseDeviceBuildPossible(int N) const
{
    // Split: the build needs the dense N x N matrix on the own device (8 N^2 bytes); keep it to sizes
    // that fit comfortably, larger systems keep the host build.
    return !splitWanted(*m_impl, N) || N <= 20000;
}

EEQSolverGPU::LocalInverseDeviceStats EEQSolverGPU::localInverseDeviceStats() const { return m_impl->m_li_stats; }

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

    impl.m_frag_topo_valid       = false;
    impl.m_frag_refactored       = false;
    impl.m_pcg_warm_valid        = false; // WP7-C: warm-start invalid after topology change
    impl.m_pcg_M_inv_valid       = false; // WP7-C: re-extract M_inv on next solve
    impl.m_pcg_block_jacobi_valid = false; // WP7-D: rebuild block-Jacobi after topology change
    impl.m_nfrag_batched         = nfrag;

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
    // Claude Generated (Sep 2026, F-4): one fragment is uploaded only for the projected PCG
    // (WP7-E); its per-fragment dense block would be a second N x N matrix (429 MB at 7320 atoms)
    // that no path reads for nfrag == 1 (the batched solvers run only for nfrag > 1).
    if (total_A > 0 && nfrag > 1) impl.d_A_blocks.alloc(total_A);
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
        if (total_rhs > 0 && nfrag > 1) {
            impl.d_rhs_blocks_init.alloc(total_rhs);
            checkCudaEEQ(cudaMemcpy(impl.d_rhs_blocks_init.ptr, h_rhs_init.data(),
                                     total_rhs * sizeof(double), cudaMemcpyHostToDevice),
                         "H2D rhs_blocks_init");
        }
    }

    // --- 6. Query cuSOLVER workspace size for the largest fragment ---
    if (max_frag_N > 0 && nfrag > 1) {
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
        {   // CSR (Oct 2026): fixed-order fragment sums, see k_frag_sum_csr
            std::vector<int> off(nfrag + 1, 0), idx;
            for (int i = 0; i < N; ++i)
                if (h_atom_frag[i] >= 0 && h_atom_frag[i] < nfrag) ++off[h_atom_frag[i] + 1];
            for (int f = 0; f < nfrag; ++f) off[f + 1] += off[f];
            idx.assign(off[nfrag], 0);
            std::vector<int> fill(off.begin(), off.end() - 1);
            for (int i = 0; i < N; ++i)
                if (h_atom_frag[i] >= 0 && h_atom_frag[i] < nfrag) idx[fill[h_atom_frag[i]]++] = i;
            impl.d_frag_off.alloc(nfrag + 1);
            impl.d_frag_off.upload(off.data(), nfrag + 1);
            impl.d_frag_idx.alloc(std::max<size_t>(idx.size(), 1));
            if (!idx.empty()) impl.d_frag_idx.upload(idx.data(), idx.size());
        }
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

    // WP7-C: persistent warm-start buffers (sized to current topology).
    // Per-call PCG scratch is grown lazily inside solveSinglePCG.
    impl.d_z1_persistent.alloc(N);
    impl.d_Z2_persistent.alloc(N * nfrag);

    // WP7-E: projected-PCG topology-constant buffers (see EEQSolverGPUImpl for layout).
    {
        std::vector<double> h_inv_atoms(nfrag, 0.0);
        for (int f = 0; f < nfrag; ++f)
            h_inv_atoms[f] = (impl.h_frag_sizes[f] > 0) ? (1.0 / (double)impl.h_frag_sizes[f]) : 0.0;
        impl.d_frag_inv_atoms.alloc(nfrag);
        impl.d_frag_inv_atoms.upload(h_inv_atoms.data(), nfrag);
    }
    impl.d_frag_sM.alloc(nfrag);
    impl.d_frag_sum_scratch.alloc(nfrag);

    impl.d_frag_delta_scratch.alloc(nfrag);
    impl.d_rhs_constraints_dev.alloc(nfrag);
    impl.d_ppcg_q_persistent.alloc(N);
    impl.m_ppcg_warm_valid = false;  // WP7-E: warm-start invalid after topology change

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

    impl.orderAfterWorkspace(false);   // parity item 10: no blocking copy before the first kernel on this stream

    // --- 1. Build per-fragment Coulomb blocks (if refactorizing) ---
    if (do_refactor && total_pairs > 0) {
        // WP7-D: the batched path overwrites d_A_blocks with Cholesky factors, so any
        // cached PCG block-Jacobi inverse (also in d_A_blocks) is now stale.
        impl.m_pcg_block_jacobi_valid = false;
        int block = 256;
        int grid  = (total_pairs + block - 1) / block;
        k_eeq_build_fragment_matrices<<<grid, block, 0, impl.stream>>>(
            total_pairs, cx, cy, cz, d_alpha, d_gam,
            impl.d_frag_sizes.ptr, impl.d_frag_offsets_A.ptr,
            impl.d_frag_offsets_pair.ptr, impl.d_frag_atom_offsets.ptr,
            impl.d_frag_atom_map.ptr, impl.d_A_blocks.ptr, nfrag, cutoff_sq);
    }

    // --- 2. Restore col1 = 1 (dpotrs of the previous call overwrote it with Z2), then
    //        gather the per-step RHS (col0) from global d_rhs_atoms ---
    if (impl.d_rhs_blocks_init.n < total_rhs) return false;
    checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs_blocks.ptr, impl.d_rhs_blocks_init.ptr,
                                  (size_t)total_rhs * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream),
                 "frag D2D rhs_blocks reset");
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

// ============================================================================
// solve_method lu on the device: dense LU of the augmented EEQ system
// Claude Generated (Oct 2026)
// ============================================================================

/**
 * @brief Augmented EEQ matrix [A C^T; C 0] (M = N + nfrag rows), column-major, one thread per element.
 *
 * Atom block: same element formula and operand order as k_eeq_build_cols / k_eeq_build_matrix
 * (i = max(row,col), j = min(row,col)), so A_nn is bit-identical to the Cholesky paths' matrix.
 * Constraint block: C(f, a) = 1 when atom a belongs to fragment f (atom_frag, 0-indexed), as in
 * the CPU EEQSolver augmented matrix; the nfrag x nfrag corner is zero.
 */
__global__ __launch_bounds__(256, 4)
void k_eeq_build_augmented(int N, int M,
                           const double* __restrict__ cx, const double* __restrict__ cy,
                           const double* __restrict__ cz, const double* __restrict__ alpha,
                           const double* __restrict__ gam, const int* __restrict__ atom_frag,
                           double* __restrict__ Aug, double cutoff_sq)
{
    const long long tid = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long long tot = static_cast<long long>(M) * M;
    if (tid >= tot) return;
    const int row = static_cast<int>(tid % M);
    const int col = static_cast<int>(tid / M);
    double v;
    if (row < N && col < N) {
        if (row == col) {
            v = gam[row] + TSQRT2PI_GPU / sqrt(alpha[row]);
        } else {
            const int i = row > col ? row : col;
            const int j = row > col ? col : row;
            const double dx = cx[i] - cx[j];
            const double dy = cy[i] - cy[j];
            const double dz = cz[i] - cz[j];
            const double r_sq = dx * dx + dy * dy + dz * dz;
            if (cutoff_sq > 0.0 && r_sq > cutoff_sq) {
                v = 0.0;
            } else {
                const double r = sqrt(r_sq);
                const double asum = alpha[i] + alpha[j];
                // Far field (Claude Generated, Oct 2026): bit-identical, see k_eeq_build_matrix
                if (r_sq >= GFNFF_COUL_FARFIELD_X2 * asum) {
                    v = 1.0 / r;
                } else {
                    const double gamma_ij = 1.0 / sqrt(asum);
                    v = erf(gamma_ij * r) / r;
                }
            }
        }
    } else if (row >= N && col >= N) {
        v = 0.0;
    } else {
        const int atom = row < N ? row : col;
        const int frag = row < N ? col - N : row - N;
        v = (atom_frag[atom] == frag) ? 1.0 : 0.0;
    }
    Aug[tid] = v;
}

bool EEQSolverGPU::solveWithDeviceRHSAndGPUAugmentedLU(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha,
    const double* d_gam,
    const double* d_rhs_atoms,
    const std::vector<double>& rhs_constraints,
    double cutoff_sq)
{
    auto& impl = *m_impl;
    if (nfrag < 1 || !impl.m_frag_topo_valid || impl.m_nfrag_batched != nfrag
        || impl.d_atom_frag.n < natoms) {
        CurcumaLogger::warn("EEQ GPU lu: fragment topology not on the device");
        return false;
    }
    const int N = natoms;
    const long long M_ll = static_cast<long long>(N) + nfrag;
    if (M_ll * M_ll > static_cast<long long>(std::numeric_limits<int>::max())) {
        CurcumaLogger::warn(fmt::format("EEQ GPU lu: augmented system too large ({} rows)", M_ll));
        return false;
    }
    const int M = static_cast<int>(M_ll);

    impl.orderAfterWorkspace(false);   // parity item 10

    impl.d_aug_A.ensure(M * M);
    impl.d_aug_rhs.ensure(M);
    impl.d_aug_piv.ensure(M);

    // 1. Augmented matrix on the device.
    {
        const long long tot = static_cast<long long>(M) * M;
        const int block = 256;
        const int grid = static_cast<int>((tot + block - 1) / block);
        k_eeq_build_augmented<<<grid, block, 0, impl.stream>>>(
            N, M, cx, cy, cz, d_alpha, d_gam, impl.d_atom_frag.ptr, impl.d_aug_A.ptr, cutoff_sq);
        checkCudaEEQ(cudaGetLastError(), "k_eeq_build_augmented launch");
    }

    // 2. RHS [b_atoms; Q_frag].
    std::vector<double> h_q(nfrag, 0.0);
    for (int f = 0; f < nfrag && f < static_cast<int>(rhs_constraints.size()); ++f)
        h_q[f] = rhs_constraints[f];
    checkCudaEEQ(cudaMemcpyAsync(impl.d_aug_rhs.ptr, d_rhs_atoms, N * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream), "D2D lu rhs atoms");
    checkCudaEEQ(cudaMemcpyAsync(impl.d_aug_rhs.ptr + N, h_q.data(), nfrag * sizeof(double),
                                  cudaMemcpyHostToDevice, impl.stream), "H2D lu rhs constraints");

    // 3. getrf (partial pivoting) + getrs, one RHS column.
    int lwork = 0;
    checkCusolverEEQ(cusolverDnDgetrf_bufferSize(impl.cusolver_handle, M, M, impl.d_aug_A.ptr, M, &lwork),
                     "getrf_bufferSize (augmented lu)");
    if (impl.m_aug_ws_size < lwork || impl.d_aug_ws.n < lwork) {
        impl.d_aug_ws.alloc(lwork);
        impl.m_aug_ws_size = lwork;
    }
    checkCusolverEEQ(cusolverDnDgetrf(impl.cusolver_handle, M, M, impl.d_aug_A.ptr, M,
                                      impl.d_aug_ws.ptr, impl.d_aug_piv.ptr, impl.d_info.ptr),
                     "cusolverDnDgetrf (augmented lu)");
    checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_info.ptr, sizeof(int),
                                  cudaMemcpyDeviceToHost, impl.stream), "D2H info (augmented lu)");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after getrf (augmented lu)");
    if (*impl.h_info != 0) {
        CurcumaLogger::warn(fmt::format("EEQ GPU lu: getrf failed (info={})", *impl.h_info));
        return false;
    }
    checkCusolverEEQ(cusolverDnDgetrs(impl.cusolver_handle, CUBLAS_OP_N, M, 1,
                                      impl.d_aug_A.ptr, M, impl.d_aug_piv.ptr,
                                      impl.d_aug_rhs.ptr, M, impl.d_info.ptr),
                     "cusolverDnDgetrs (augmented lu)");

    // 4. Charges = first N entries of the solution -> d_rhs[0..N-1] (getDeviceChargesPtr()).
    if (impl.d_rhs.n < N) impl.d_rhs.alloc(N);
    checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, impl.d_aug_rhs.ptr, N * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream), "D2D lu charges");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after augmented lu");
    return true;
}

// ============================================================================
// solve_method ldlt on the device: Bunch-Kaufman LDL^T of A_nn + Schur complement
// Claude Generated (Oct 2026)
// ============================================================================

__global__ void k_int_to_int64(int n, const int* __restrict__ in, long long* __restrict__ out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = static_cast<long long>(in[i]);
}

bool EEQSolverGPU::solveWithDeviceRHSAndGPUSchurLDLT(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha,
    const double* d_gam,
    const double* d_rhs_atoms,
    const std::vector<double>& rhs_constraints,
    double cutoff_sq)
{
    auto& impl = *m_impl;
    if (nfrag < 1 || !impl.m_frag_topo_valid || impl.m_nfrag_batched != nfrag
        || impl.d_atom_frag.n < natoms || impl.d_rhs_constraint_cols.n < natoms * nfrag) {
        CurcumaLogger::warn("EEQ GPU ldlt: fragment topology not on the device");
        return false;
    }
    const int N = natoms;
    const int nrhs = nfrag + 1;
    if (static_cast<long long>(N) * N > static_cast<long long>(std::numeric_limits<int>::max())) return false;

    impl.orderAfterWorkspace(false);   // parity item 10

    // d_A is shared with the Cholesky/LU paths: their cached factor is gone after this call.
    impl.m_has_cached_factor  = false;
    impl.m_using_lu           = false;
    impl.m_cached_factor_mixed = false;
    impl.d_A.ensure(N * N);
    impl.d_ldlt_piv.ensure(N);
    impl.d_ldlt_piv64.ensure(N);
    impl.d_rhs.ensure(N * nrhs);

    // 1. A_nn (same kernel and values as WP5-A/WP7-A; both triangles written, sytrf reads the lower).
    {
        const int n_lower = N * (N + 1) / 2;
        const int block = 256;
        k_eeq_build_matrix<<<(n_lower + block - 1) / block, block, 0, impl.stream>>>(
            N, cx, cy, cz, d_alpha, d_gam, impl.d_A.ptr, cutoff_sq);
    }

    // 2. LDL^T (Bunch-Kaufman), lower triangle.
    int lwork = 0;
    checkCusolverEEQ(cusolverDnDsytrf_bufferSize(impl.cusolver_handle, N, impl.d_A.ptr, N, &lwork),
                     "sytrf_bufferSize");
    impl.d_ldlt_ws.ensure(std::max(1, lwork));
    checkCusolverEEQ(cusolverDnDsytrf(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER, N, impl.d_A.ptr, N,
                                      impl.d_ldlt_piv.ptr, impl.d_ldlt_ws.ptr, std::max(1, lwork),
                                      impl.d_info.ptr),
                     "cusolverDnDsytrf");
    checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_info.ptr, sizeof(int),
                                  cudaMemcpyDeviceToHost, impl.stream), "D2H info (sytrf)");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after sytrf");
    if (*impl.h_info != 0) {
        CurcumaLogger::warn(fmt::format("EEQ GPU ldlt: sytrf failed (info={})", *impl.h_info));
        return false;
    }
    k_int_to_int64<<<(N + 255) / 256, 256, 0, impl.stream>>>(N, impl.d_ldlt_piv.ptr, impl.d_ldlt_piv64.ptr);

    // 3. RHS [b_atoms | C^T], as in WP7-A.
    checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr, d_rhs_atoms, N * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream), "D2D rhs_atoms col0 (ldlt)");
    checkCudaEEQ(cudaMemcpyAsync(impl.d_rhs.ptr + N, impl.d_rhs_constraint_cols.ptr,
                                  (size_t)N * nfrag * sizeof(double),
                                  cudaMemcpyDeviceToDevice, impl.stream), "D2D constraint cols (ldlt)");

    // 4. Solve A_nn [z1 | Z2] = [b | C^T], in blocks of at most 1024 right-hand sides:
    //    cusolverDnXsytrs of cuSOLVER 12.2 (CUDA 13) returns wrong columns, with info = 0, for
    //    nrhs > 1024 (standalone check, SPD and indefinite, n = 500..7320: residual 2.5e2 at
    //    nrhs = 1025, 1e-15 at 1024 and for any chunked call).
    constexpr int kSytrsMaxRhs = 1024;
    size_t ws_dev = 0, ws_host = 0;
    const int64_t n64 = N;
    checkCusolverEEQ(cusolverDnXsytrs_bufferSize(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER, n64,
                                                 static_cast<int64_t>(std::min(nrhs, kSytrsMaxRhs)),
                                                 CUDA_R_64F, impl.d_A.ptr, n64,
                                                 reinterpret_cast<const int64_t*>(impl.d_ldlt_piv64.ptr),
                                                 CUDA_R_64F, impl.d_rhs.ptr, n64, &ws_dev, &ws_host),
                     "Xsytrs_bufferSize");
    impl.d_ldlt_xs.ensure(static_cast<int>(std::max<size_t>(1, ws_dev)));
    if (impl.h_ldlt_xs.size() < std::max<size_t>(1, ws_host)) impl.h_ldlt_xs.resize(std::max<size_t>(1, ws_host));
    for (int c0 = 0; c0 < nrhs; c0 += kSytrsMaxRhs) {
        const int64_t nc = std::min(kSytrsMaxRhs, nrhs - c0);
        checkCusolverEEQ(cusolverDnXsytrs(impl.cusolver_handle, CUBLAS_FILL_MODE_LOWER, n64, nc,
                                          CUDA_R_64F, impl.d_A.ptr, n64,
                                          reinterpret_cast<const int64_t*>(impl.d_ldlt_piv64.ptr),
                                          CUDA_R_64F, impl.d_rhs.ptr + static_cast<size_t>(c0) * N, n64,
                                          impl.d_ldlt_xs.ptr, ws_dev, impl.h_ldlt_xs.data(), ws_host,
                                          impl.d_info.ptr),
                         "cusolverDnXsytrs");
    }
    checkCudaEEQ(cudaMemcpyAsync(impl.h_info, impl.d_info.ptr, sizeof(int),
                                  cudaMemcpyDeviceToHost, impl.stream), "D2H info (sytrs)");
    checkCudaEEQ(cudaStreamSynchronize(impl.stream), "sync after sytrs");
    if (*impl.h_info != 0) {
        CurcumaLogger::warn(fmt::format("EEQ GPU ldlt: sytrs failed (info={})", *impl.h_info));
        return false;
    }

    // 5. Schur complement over the fragment constraints -> charges in d_rhs[0..N-1].
    finishGeneralSchur(impl, N, nfrag, rhs_constraints);
    return true;
}
