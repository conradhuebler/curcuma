/*
 * <Native xTB GPU Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): cuSOLVER/cuBLAS handle + stream lifecycle. Kept
 * free of CurcumaLogger so nvcc never sees host-only logging headers; the host
 * wrapper does all user-facing logging based on ok()/deviceName().
 */

#ifdef USE_CUDA

#include "xtb_gpu_context.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <set>
#include <thread>
#include <vector>

// CudaBuffer<T> (RAII cudaMalloc/cudaFree) — the project's established device
// allocation helper; routing every allocation through it is the mitigation for
// the GFN-FF GPU heap-corruption class of bug.
#include "../../ff_methods/cuda/gfnff_soa.h"

// Stage 3 device integral functions (CN counting; later overlap/multipole) and
// the host element-parameter tables that seed __constant__ memory.
#include "xtb_gpu_integrals_device.cuh"
#include "../parameters/xtb_params_extra.hpp"
#include "xtb_distributed_eigensolver.h"

namespace curcuma {
namespace xtb {
namespace gpu {

// Claude Generated (Sep 2026): used by Impl::distSolve below, defined with the other kernels.
__global__ void k_probe_fill(float* v, int n, unsigned seed);
__global__ void k_scale_by(float* x, const float* s, int n);
__global__ void k_probe_filld(double* v, int n, unsigned seed);
__global__ void k_scale_byd(double* x, const double* s, int n);

struct XtbGpuContext::Impl {
    cudaStream_t       stream   = nullptr;
    cublasHandle_t     cublas   = nullptr;
    cusolverDnHandle_t cusolver = nullptr;
    int                device   = -1;
    std::string        name;
    bool               ok       = false;

    // Device-resident GFN1 SCF state (Stage 2). Allocated once per geometry in
    // residentBegin and reused across SCF iterations: H0/S/L are uploaded once,
    // C/P are produced on the device and never leave until residentFinalize.
    int                resident_n = 0;
    CudaBuffer<double> dH0, dS, dL;      // geometry-constant (uploaded once)
    CudaBuffer<double> dC, dP, dCw;      // eigenvectors / density / scaled cols
    CudaBuffer<double> dEps, dVao, dOcc, dPop;  // length-n working vectors
    CudaBuffer<double> dWork;            // cuSOLVER dsyevd workspace
    CudaBuffer<int>    dInfo;            // cuSOLVER devInfo
    int                lwork = 0;

    // Mixed-precision (FP32) eigensolve buffers — far-from-convergence SCF
    // iterations solve in FP32 (≫ faster on consumer FP64-limited GPUs), the
    // FP64 path takes over near convergence so the converged energy is FP64.
    CudaBuffer<float>  dCf, dLf, dWorkf; // F/eigenvectors, lower L, ssyevd workspace
    CudaBuffer<float>  dEpsf;            // FP32 eigenvalues
    int                lwork_f32 = 0;

    // GFN2 multipole (Stage 2b). dDpInt holds the 3 dipole AO-integral matrices
    // contiguously (3·n·n), dQpInt the 6 quadrupole matrices (6·n·n); both
    // geometry-constant. dVdp/dVqp are the per-iteration multipole potentials
    // (3·nat / 6·nat), dDpAt/dQpAt the per-iteration atomic moments.
    int                resident_nat = 0;
    CudaBuffer<double> dDpInt, dQpInt;   // 3·nn / 6·nn (uploaded once)
    CudaBuffer<int>    dAo2at;           // AO→atom map, length n (uploaded once)
    CudaBuffer<double> dVdp, dVqp;       // 3·nat / 6·nat (per iteration)
    CudaBuffer<double> dDpAt, dQpAt;     // 3·nat / 6·nat (per iteration)

    // Stage 3: device-side integral build. The flattened basis is molecule-
    // constant (uploaded once in beginBasis); only dXyz changes per geometry.
    int                basis_nat = 0;    // atoms in the uploaded basis
    int                basis_nsh = 0;    // shells in the uploaded basis
    int                basis_nao = 0;    // AOs in the uploaded basis
    int                basis_is_gfn2 = 0;
    CudaBuffer<int>    dZ;               // atomic numbers, length nat
    CudaBuffer<int>    dSh2at;           // shell → atom map, length nsh
    CudaBuffer<int>    dAng;             // angular momentum per shell, length nsh
    CudaBuffer<int>    dIaoSh;           // first-AO offset per shell, length nsh
    CudaBuffer<int>    dNaoSh;           // AOs per shell, length nsh
    CudaBuffer<int>    dShNprim;         // primitives per shell, length nsh
    CudaBuffer<int>    dShPrimOff;       // prim offset per shell, length nsh
    CudaBuffer<int>    dAo2sh;           // AO → shell map, length nao (GFN2 multipole)
    CudaBuffer<int>    dValence;         // GFN1 valence flags per shell, length nsh
    CudaBuffer<double> dPrimAlpha;       // flattened primitive exponents
    CudaBuffer<double> dPrimCoeff;       // flattened primitive coefficients
    CudaBuffer<double> dShZeta;          // slater_exp per shell, length nsh
    CudaBuffer<double> dSelfE0;          // raw shell self-energies, length nsh
    CudaBuffer<double> dKcn;             // CN coefficients, length nsh
    CudaBuffer<double> dShpoly;          // distance-polynomial coeff, length nsh
    CudaBuffer<double> dHardness;        // Coulomb per-shell hardness, length nsh
    CudaBuffer<double> dXyz;             // geometry in bohr, length 3·nat (per geometry)
    CudaBuffer<double> dCN;              // coordination numbers, length nat (resident)
    CudaBuffer<double> dSE;              // CN-shifted self-energies, length nsh (resident)
    CudaBuffer<double> dGamma;           // Coulomb γ matrix, nsh² col-major (resident)
    CudaBuffer<double> dSdR;             // overlap derivative dS/dR_A, 3·nn (Stage 4)
    CudaBuffer<double> dPotrfWork;       // cuSOLVER potrf workspace (device chol)
    int                potrf_lwork = 0;

    // Stage 4 gradient: per-atom repulsion params (molecule-constant) + the
    // gradient work buffers. dGrad layout [3*i+k] Eh/Bohr.
    CudaBuffer<double> dRepAlpha, dRepZeff;  // per-atom, length nat
    CudaBuffer<double> dW;                    // energy-weighted density, nao²
    CudaBuffer<double> dQsh;                  // shell charges, length nsh
    CudaBuffer<double> dGrad;                 // gradient, 3·nat
    CudaBuffer<double> dEdcn;                 // CN coupling, length nat

    // Stage 5 (Part A): single-shot D4 EEQ charge model. Self-contained — does
    // NOT depend on beginBasis; eeqCharges allocates/uploads its own length-N
    // params + (N+1)×(N+1) augmented matrix and keeps the LU factor + CN resident
    // so eeqChargeResponseGradient reuses them for the adjoint solve. Claude Generated.
    int                eeq_n = 0;             // N for which the LU factor is valid
    CudaBuffer<double> dEeqXyz;               // 3·N geometry (Bohr)
    CudaBuffer<double> dEeqChi, dEeqGam;      // per-atom χ, γ (length N)
    CudaBuffer<double> dEeqAlp, dEeqCnf;      // per-atom α² , κ (length N)
    CudaBuffer<double> dEeqRcov;              // per-atom 4/3·rcov·Å→Bohr (length N)
    CudaBuffer<double> dEeqCn, dEeqCnRaw;     // log-compressed CN + raw CN (length N)
    CudaBuffer<double> dEeqM;                 // augmented matrix / LU factor, (N+1)²
    CudaBuffer<double> dEeqRhs;               // RHS [b;Q] → [q;λ], length N+1
    CudaBuffer<double> dEeqQ;                 // atomic charges q (length N, resident)
    CudaBuffer<double> dEeqAdjRhs;            // adjoint RHS [dEdq;0] → [z;…], length N+1
    CudaBuffer<double> dEeqDedq;              // uploaded dE_D4/dq (length N)
    CudaBuffer<double> dEeqU;                 // per-atom CN-response weight (length N)
    CudaBuffer<double> dEeqGrad;              // response gradient [3·N], [3a+k]
    CudaBuffer<double> dEeqWork;              // cuSOLVER getrf workspace
    CudaBuffer<int>    dEeqIpiv;              // LU pivots, length N+1
    int                eeq_lwork = 0;

    // Stage 5 (Part B1): device atomic Mulliken charges from the resident density.
    // dQat(A) = n0_at(A) − Σ_{μ∈A} pop_ao(μ); reduced from the resident dPop via
    // the resident dAo2at map, kept resident for the in-SCF D4 potential (B2).
    CudaBuffer<double> dN0at;                 // reference atom occupations (nat)
    CudaBuffer<double> dQat;                  // atomic charges q_at (nat, resident)
    // Stage 6 (S6.2): reference shell occupations for the resident q_sh reduction
    // (q_sh = n0_sh − Σ_{μ∈s} pop_ao(μ), scattered via dAo2sh into the resident
    // dQsh). Mirrors dN0at/dQat for shells. Claude Generated.
    CudaBuffer<double> dN0sh;                 // reference shell occupations (nsh)

    // Stage 5 (Part B2): in-SCF GFN2 D4 atom-potential dE_D4/dq on the device.
    // The CN-Gaussian + zeta weights are built on the host (buildRefWFlat) and
    // uploaded per iteration as W/dWq (nat·MAX_REF); the device runs the O(N²)
    // 7×7 reference contraction × BJ disp_sum → dEdq(A). Reference data
    // (c6_flat, sqrtZr4r2, refn, Z) + geometry + BJ params upload once/geometry.
    static constexpr int D4_MAX_REF = 7;      // mirrors D4ParameterGenerator::MAX_REF
    int                d4_nat = 0;
    double             d4_s6 = 0, d4_s8 = 0, d4_a1 = 0, d4_a2 = 0, d4_cut = 0;
    bool               d4_c6_uploaded = false;  // c6_flat is element data → upload once/process
    CudaBuffer<int>    dD4Z;                   // atomic numbers (nat)
    CudaBuffer<int>    dD4Nref;                // reference count per atom (nat)
    CudaBuffer<double> dD4Sqrt;                // sqrtZr4r2 per atom (nat)
    CudaBuffer<double> dD4Xyz;                 // geometry, Bohr (3·nat)
    CudaBuffer<double> dD4C6Flat;              // reference C6 block (MAX_ELEM²·MAX_REF²)
    CudaBuffer<double> dD4W, dD4dWq;           // per-iter weights (nat·MAX_REF)
    CudaBuffer<double> dD4Dedq;                // output dE_D4/dq (nat)
    // Claude Generated (Sep 2026): post-SCF 2-body gradient + ATM on the device.
    CudaBuffer<double> dD4dWc, dD4Eat, dD4Grad, dD4Dcn;   // nat*7 / nat / 3 nat / nat
    CudaBuffer<double> dD4AtmC6, dD4AtmDc6;               // nat^2 each (q=0 reference)
    CudaBuffer<int>    dD4NbPtr, dD4Nb;                   // ATM neighbour list
    std::vector<double> h_d4_xyz;                         // host copy of the D4 geometry

    // Stage 6 (S6.2b): q-independent per-atom reference data so the device rebuilds
    // W/dWq from the resident SCF charges (k_d4_build_refw), removing the host
    // buildRefWFlat + W/dWq upload from the loop. Uploaded once per geometry by
    // beginDispersionWeights. Claude Generated.
    CudaBuffer<double> dD4Cn;                  // geometry-fixed CN per atom (nat)
    CudaBuffer<double> dD4Gi;                  // eta·gc per atom (nat)
    CudaBuffer<double> dD4Zeff;                // effective nuclear charge per atom (nat)
    CudaBuffer<double> dD4Refcn;              // ngw-bucketing reference CN (nat·MAX_REF)
    CudaBuffer<double> dD4Refcovcn;          // CN-Gaussian reference covCN (nat·MAX_REF)
    CudaBuffer<double> dD4Refq;              // reference charges (nat·MAX_REF)

    // Stage 5 (Part B3/B4): full device GFN2 potential build. The geometry-fixed
    // multipole interaction matrices (amat_*, nat²) + the per-shell third-order
    // hardness gamma3 upload once per geometry (beginPotential); residentSolvePotential
    // then builds v_sh (γ·q_sh + third-order) + the multipole v_dp/v_qp/v_at scalar
    // shift + the resident D4 dEdq, expands v_ao, and folds into the Fock+eigensolve —
    // so the host uploads only q_sh/dp_at/qp_at (+ the host-built D4 W/dWq) per iter.
    int                pot_nsh = 0;
    CudaBuffer<double> dMpAmatSD;   // charge-dipole, 3·nat² (col-major blocks [k])
    CudaBuffer<double> dMpAmatDD;   // dipole-dipole, 9·nat² (block [a*3+b])
    CudaBuffer<double> dMpAmatSQ;   // charge-quadrupole, 6·nat²
    // Claude Generated (Sep 2026): on-the-fly interaction (no stored nat² matrices).
    bool               mp_otf = false;
    double             mp_dmp3 = 3.0, mp_dmp5 = 4.0;
    CudaBuffer<double> dMpXyz;      // geometry (Bohr), 3·nat
    CudaBuffer<double> dMpRad;      // CN-dependent damping radii, nat
    CudaBuffer<double> dMpDkernel;  // on-site dipole XC kernel, nat
    CudaBuffer<double> dMpQkernel;  // on-site quadrupole XC kernel, nat
    CudaBuffer<double> dGamma3;     // per-shell third-order hardness Γ_s, nsh
    CudaBuffer<double> dPotQsh;     // uploaded shell charges q_sh, nsh
    CudaBuffer<double> dInDpAt;     // uploaded atomic dipoles dp_at, 3·nat
    CudaBuffer<double> dInQpAt;     // uploaded atomic quadrupoles qp_at, 6·nat
    CudaBuffer<double> dVsh;        // shell potential v_sh, nsh
    CudaBuffer<double> dVat;        // atom potential v_at, nat

    // WP4b (Claude Generated June 2026): in-SCF implicit-solvation reaction field.
    // dSolvB is the nat×nat Born interaction matrix B (keps-scaled, incl. self-energy
    // + HB diagonal; symmetric, column-major). When solv_active, the device potential
    // build adds v_at += B·q_at (GFN2 Mulliken charges) so the SCF feels the solvent.
    CudaBuffer<double> dSolvB;      // Born matrix B, nat·nat
    bool               solv_active = false;

    // Stage 6 (S6.1): device occupation. k_occupations fills the resident dOcc
    // from the resident dEps (no host round-trip in the loop); dOccMu/dOccNcol
    // carry the converged chemical potential + the occupied-column count for the
    // component test. Claude Generated.
    CudaBuffer<double> dOccMu;      // chemical potential µ (1 double)
    CudaBuffer<int>    dOccNcol;    // last column with occ>1e-12 (1 int)

    // Stage 6 (S6.3): device SCC energy. dEScratch holds γ·q_sh (nsh) for the
    // Coulomb dot; dESca accumulates the third-order + multipole scalars (2 doubles,
    // atomicAdd targets). Claude Generated.
    CudaBuffer<double> dEScratch;   // γ·q_sh, nsh
    CudaBuffer<double> dESca;       // [E_third, E_multipole]

    // Stage 6 (S6.4): device Broyden mixer state. The packed SCC vector has length
    // broyden_N (= nsh + 9·nat for GFN2). The history dF/u live as the first M
    // columns of N×max_hist column-major matrices (a ring buffer over the slots —
    // vnext is invariant to the column order, so dropping the oldest slot matches
    // the host FIFO); vin_last/F_last + the M×M Gram solve scratch stay resident
    // across SCF iterations. Claude Generated.
    int    broyden_N = 0, broyden_iter = 0, broyden_push = 0, broyden_maxhist = 20;
    double broyden_alpha = 0.25, broyden_w0 = 0.01;
    CudaBuffer<double> dBroyVin, dBroyVout, dBroyVnext;          // scratch (test upload)
    CudaBuffer<double> dBroyF, dBroyFLast, dBroyVinLast, dBroyDFtmp;  // length N
    CudaBuffer<double> dBroyDFmat, dBroyUmat;                    // N × max_hist (col-major)
    CudaBuffer<double> dBroyGram, dBroyC, dBroyGamma;            // M×M / M / M

    // Stage 6 (S6.5): fused device-driven loop. The geometry-fixed scalars + the
    // convergence dq scratch; everything else reuses the resident Stage 2-5/6
    // buffers. The host polls only dq + the 4 energy scalars (O(1)) per step.
    int    loop_nsh = 0, loop_nat = 0, loop_nao = 0, loop_nocc_pairs = 0;
    double loop_Tele = 0.0, loop_nelec = 0.0;
    CudaBuffer<double> dDq;         // max|q_sh_out − q_sh_in| (1 double)

    // Claude Generated (Sep 2026): pre-allocation memory check (beginBasis).
    bool        memory_check = true;
    std::string last_error;

    // Claude Generated (Sep 2026): phase profiler for the device path, enabled by the
    // environment variable CURCUMA_GPU_PROFILE. Each mark synchronises the stream (kernels are
    // asynchronous, so host clocks are meaningless without it) and accumulates the wall time
    // since the previous mark under a phase name. Disabled = no synchronisation, no cost.
    bool prof = std::getenv("CURCUMA_GPU_PROFILE") != nullptr;
    // Claude Generated (Sep 29, 2026): the thread-per-atom originals of the warp-per-atom pair kernels.
    bool thread_per_atom = std::getenv("CURCUMA_GPU_THREAD_PER_ATOM") != nullptr;
    std::vector<std::string> prof_names;
    std::vector<double>      prof_ms;
    std::vector<int>         prof_calls;
    std::vector<double>      prof_mem_mib;   // largest device memory in use at a mark of that phase
    std::chrono::steady_clock::time_point prof_t;
    void profStart()
    {
        if (!prof) return;
        cudaStreamSynchronize(stream);
        prof_t = std::chrono::steady_clock::now();
    }
    void profMark(const char* name)
    {
        if (!prof) return;
        cudaStreamSynchronize(stream);
        const auto t = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t - prof_t).count();
        prof_t = t;
        profAdd(name, ms);
    }
    /// Add a timing without moving the mark (sub-phases measured elsewhere, e.g. inside the
    /// multi-GPU eigensolver). Also records the device memory in use on this context's device.
    void profAdd(const char* name, double ms)
    {
        if (!prof) return;
        size_t free_b = 0, total_b = 0;
        double used = 0.0;
        if (cudaMemGetInfo(&free_b, &total_b) == cudaSuccess)
            used = static_cast<double>(total_b - free_b) / (1024.0 * 1024.0);
        for (size_t i = 0; i < prof_names.size(); ++i) {
            if (prof_names[i] == name) {
                prof_ms[i] += ms; ++prof_calls[i];
                prof_mem_mib[i] = std::max(prof_mem_mib[i], used);
                return;
            }
        }
        prof_names.emplace_back(name);
        prof_ms.push_back(ms);
        prof_calls.push_back(1);
        prof_mem_mib.push_back(used);
    }

    // Claude Generated (Sep 2026): screened (sparse) AO-pair storage for S, H0, dp_int, qp_int.
    // See the "Screened (sparse) AO-pair storage" kernel block for the layout.
    // Claude Generated (Sep 2026, multi-GPU): the pattern density P(r,c) = sum_k Cw(r,k) C(c,k)
    // is a sum over the occupied columns, so it splits over GPUs without any approximation: every
    // device evaluates the SAME pattern over its own slice of k and the partial results are added.
    // Only a column slice of C (n x kn) travels per step, not the whole matrix.
    struct DensityHelper {
        int          device = -1;
        cudaStream_t stream = nullptr;
        cudaEvent_t  done   = nullptr;
        double*      dC     = nullptr;   // n x kn slice of the eigenvectors
        double*      dCw    = nullptr;   // the same slice scaled by the occupations
        double*      dOcc   = nullptr;   // kn occupations
        double*      dPsp   = nullptr;   // partial density on the pattern (nnz)
        int*         dRow   = nullptr;   // pattern (geometry-constant)
        int*         dCol   = nullptr;
        int          cap_cols = 0;
        int          cap_nnz  = 0;
        long         pattern_gen = -1;
        double*      stage  = nullptr;   // nnz staging buffer on the PRIMARY device
    };
    std::vector<int>           dens_devices;    // helper devices (never the context's own)
    std::vector<DensityHelper> dens_helpers;
    bool                       dens_failed = false;
    int                        dens_min_nao = 4000;
    std::string                dens_error;
    std::string                dens_skip_reason;   ///< G2-6: why the split did not run (last cause)
    int                        dens_steps = 0;
    long                       pattern_generation = 0;   // bumped when the screened pattern changes

    void releaseDensityHelpers()
    {
        for (DensityHelper& h : dens_helpers) {
            cudaSetDevice(h.device);
            for (void** p : { reinterpret_cast<void**>(&h.dC), reinterpret_cast<void**>(&h.dCw),
                              reinterpret_cast<void**>(&h.dOcc), reinterpret_cast<void**>(&h.dPsp),
                              reinterpret_cast<void**>(&h.dRow), reinterpret_cast<void**>(&h.dCol) })
                if (*p) { cudaFree(*p); *p = nullptr; }
            if (h.done) { cudaEventDestroy(h.done); h.done = nullptr; }
            if (h.stream) { cudaStreamDestroy(h.stream); h.stream = nullptr; }
            cudaSetDevice(device);
            if (h.stage) { cudaFree(h.stage); h.stage = nullptr; }
            h.cap_cols = h.cap_nnz = 0;
            h.pattern_gen = -1;
        }
        dens_helpers.clear();
        cudaSetDevice(device);
    }

    // Claude Generated (Sep 2026, multi-GPU step 3): optional multi-GPU eigensolve, created on
    // first use. dist_failed latches after a failure so the SCF does not retry every iteration.
    std::vector<int> dist_devices;
    std::string      dist_backend = "auto";
    int              dist_block = 128;
    int              dist_min_nao = 4000;
    int              dist_potrf = 0;          // Cholesky factorisations done on the devices
    void*            dist_backup = nullptr;   // FP64 verification copy of F on a helper device
    int              dist_backup_dev = -1;
    size_t           dist_backup_bytes = 0;
    std::string      dist_potrf_status;       // why one was not (empty: never refused)
    bool             dist_fp32 = true;
    bool             dist_failed = false;
    // Claude Generated (Sep 2026): -1 = this backend's FP32 solve was verified WRONG on this
    // machine, 1 = verified correct, 0 = not tested yet. cusolverMg FP32 failed here (4x A4500,
    // CUDA 13.3, n = 15444: the SCF diverged from the first iteration while the eigenvalue sum
    // still matched, so a trace check does not catch it) - but that is one library version on one
    // machine, so it is measured per run instead of hard-coded.
    int              dist_fp32_state = 0;
    int              dist_fp64_state = 0;   // same, for the FP64 solves (see verifyDistributed)
    bool             dist_verify = true;    // -gpu_eigensolver_verify
    CudaBuffer<float>  dProbeV, dProbeY, dProbeT, dProbeA;   // FP32 probe vectors + copy of A
    CudaBuffer<float>  dProbeU;                               // scratch for the back-transformed check
    CudaBuffer<double> dProbeVd, dProbeYd, dProbeTd, dProbeAd, dProbeSd;
    int              dist_solves = 0;
    long             l_generation = 0;   // bumped whenever dL (Cholesky factor of S) is rewritten
    std::string      dist_status;
    std::string      dist_degrade;   ///< G2-7: part of the solve stays on this device (see accessor)
    std::unique_ptr<DistributedEigensolver> dist;

    /// 1 = solved on several GPUs, 0 = not used (input intact, run the single-GPU solver),
    /// -1 = failed after the input was overwritten.
    /// With back_L != nullptr (FP32 only) the backend also back-transforms the eigenvectors,
    /// C = L^-T Q, on the devices when it can; *back_done then tells whether it did.
    int distSolve(int n, void* A, void* eig, bool fp32, const void* back_L = nullptr,
                  bool* back_done = nullptr)
    {
        if (back_done) *back_done = false;
        if (!distReady(n, fp32)) return 0;
        if (cudaStreamSynchronize(stream) != cudaSuccess) return 0;
        // First FP32 solve on this machine: keep a copy of the input so the returned eigenpairs
        // can be checked against it, and so this call can be redone on this device if they are
        // wrong. One n x n FP32 copy, once per run; if it does not fit, FP32 is not distributed
        // (an unverified FP32 backend could return silently wrong eigenvectors).
        // CURCUMA_GPU_EIG_VERIFY_ALWAYS=1 verifies EVERY distributed FP32 solve, not just the
        // first. Diagnostic: cusolverMg passed the first check here and still made the SCF
        // diverge, i.e. a later call went wrong. Claude Generated (Sep 2026).
        // EVERY distributed solve is verified, not just the first: cusolverMg passed the first
        // check here and degraded later (solve 1 residual 1.4e-6, solve 2 8.7e-2), which a one-shot
        // gate cannot catch - the run then converged to an energy 1.65 kcal/mol off. The check is
        // three matrix-vector products plus one n x n copy, i.e. well under a percent of a solve.
        bool verify = fp32 && dist_verify && dist_fp32_state >= 0;
        bool verify64 = !fp32 && dist_verify && dist_fp64_state >= 0;
        if (verify64) {
            try {
                dProbeAd.ensure(static_cast<int>(static_cast<size_t>(n) * n));
                dProbeVd.ensure(n); dProbeYd.ensure(n); dProbeTd.ensure(n);
            } catch (...) {
                dProbeAd.free(); dProbeVd.free(); dProbeYd.free(); dProbeTd.free();
                dist_fp64_state = -1;
                dist_status = std::string(dist->name()) + ": FP64 verification needs one more n x n "
                    "buffer than fits on device " + std::to_string(device)
                    + "; the eigensolve stays on that device";
                return 0;
            }
            const int b = 256;
            k_probe_filld<<<(n + b - 1) / b, b, 0, stream>>>(dProbeVd.ptr, n, 0x85EBCA6Bu);
            if (cudaGetLastError() != cudaSuccess
                || cudaMemcpyAsync(dProbeAd.ptr, A, sizeof(double) * static_cast<size_t>(n) * n,
                                   cudaMemcpyDeviceToDevice, stream) != cudaSuccess
                || cudaStreamSynchronize(stream) != cudaSuccess)
                verify64 = false;
        }
        if (verify) {
            try {
                dProbeA.ensure(static_cast<int>(static_cast<size_t>(n) * n));
                dProbeV.ensure(n); dProbeY.ensure(n); dProbeT.ensure(n);
            } catch (...) {
                dProbeA.free(); dProbeV.free(); dProbeY.free(); dProbeT.free();
                dist_fp32_state = -1;
                dist_status = std::string(dist->name()) + ": FP32 verification needs one more "
                    "n x n buffer than fits on device " + std::to_string(device)
                    + "; only FP64 is distributed";
                return 0;
            }
            const int b = 256;
            k_probe_fill<<<(n + b - 1) / b, b, 0, stream>>>(dProbeV.ptr, n, 0x9E3779B9u);
            if (cudaGetLastError() != cudaSuccess
                || cudaMemcpyAsync(dProbeA.ptr, A, sizeof(float) * static_cast<size_t>(n) * n,
                                   cudaMemcpyDeviceToDevice, stream) != cudaSuccess
                || cudaStreamSynchronize(stream) != cudaSuccess)
                verify = false;
        }
        // Claude Generated (Sep 29, 2026): FP32 back-transform on the devices (column split)
        // instead of one n x n strsm on this device. CURCUMA_GPU_EIG_LOCAL_BACKXF=1 keeps it here.
        const bool back_on_devices = fp32 && back_L && dist->supportsBackTransform()
                                     && !std::getenv("CURCUMA_GPU_EIG_LOCAL_BACKXF");
        bool solved = back_on_devices
            ? dist->solveBackTransformed(n, A, back_L, l_generation, eig, fp32, device)
            : dist->solve(n, A, eig, fp32, device);
        // Claude Generated (Sep 2026): CURCUMA_GPU_EIG_CORRUPT=1 deliberately damages the returned
        // eigenvectors, to check that verifyDistributedFp32() actually notices. A detector that is
        // never tested against a known-bad input is not a detector.
        if (solved && verify && std::getenv("CURCUMA_GPU_EIG_CORRUPT")) {
            const size_t ncorrupt = static_cast<size_t>(n) * std::max(1, n / 20);
            cudaMemsetAsync(A, 0, sizeof(float) * ncorrupt, stream);
            cudaStreamSynchronize(stream);
        }
        cudaSetDevice(device);
        if (prof) {
            double sc = 0.0, so = 0.0, ga = 0.0;
            dist->lastTimings(sc, so, ga);
            profAdd(fp32 ? "  multi-GPU FP32: scatter" : "  multi-GPU FP64: scatter", sc);
            profAdd(fp32 ? "  multi-GPU FP32: solve" : "  multi-GPU FP64: solve", so);
            profAdd(fp32 ? "  multi-GPU FP32: gather" : "  multi-GPU FP64: gather", ga);
            if (back_on_devices) {
                double re = 0.0, ba = 0.0;
                dist->lastGeneralizedTimings(re, ba);
                profAdd("  multi-GPU FP32: back-transform (in solve)", ba);
            }
        }
        if (solved) {
            ++dist_solves;
            bool rejected = false;
            if (verify && verifyDistributedFp32(n, static_cast<const float*>(A),
                                                static_cast<const float*>(eig),
                                                back_on_devices ? static_cast<const float*>(back_L) : nullptr)
                && dist_fp32_state < 0) {
                // Wrong eigenpairs: put the input back and let the caller solve on this device.
                rejected = (cudaMemcpyAsync(A, dProbeA.ptr, sizeof(float) * static_cast<size_t>(n) * n,
                                            cudaMemcpyDeviceToDevice, stream) == cudaSuccess
                            && cudaStreamSynchronize(stream) == cudaSuccess);
                --dist_solves;
            }
            if (verify64 && verifyDistributedFp64(n, static_cast<const double*>(A),
                                                  static_cast<const double*>(eig), false)
                && dist_fp64_state < 0) {
                rejected = (cudaMemcpyAsync(A, dProbeAd.ptr, sizeof(double) * static_cast<size_t>(n) * n,
                                            cudaMemcpyDeviceToDevice, stream) == cudaSuccess
                            && cudaStreamSynchronize(stream) == cudaSuccess);
                --dist_solves;
            }
            // The probe buffers stay allocated: re-allocating ~1 GB per solve cost 38 s of the
            // 222 s polymer_2x run (measured), while keeping them costs one n x n buffer.
            // releaseEigenWorkspaces() frees them when the SCF is over.
            if (rejected) return 0;            // input restored -> single-GPU path
            if (back_done) *back_done = back_on_devices;
            // Only THIS call's precision decides: an earlier FP32 rejection must not fail a
            // perfectly good FP64 solve (it did, and aborted the SCF at the first FP64 iteration).
            if ((fp32 && dist_fp32_state < 0) || (!fp32 && dist_fp64_state < 0)) return -1;
            return 1;
        }
        dist_failed = true;
        dist_status = std::string(dist->name()) + " failed (n = " + std::to_string(n)
            + (fp32 ? ", FP32" : ", FP64") + "); single-GPU solver from here on";
        return dist->inputIntact() ? 0 : -1;
    }

    /// Claude Generated (Sep 2026): drop every eigensolver workspace and FP32 copy on this device
    /// and the multi-GPU solver's per-device buffers. Called once the SCF has converged: the
    /// post-SCF phase (dense P rebuild, D4, gradient W) otherwise stacks on top of them (polymer_2x:
    /// the device peak sat between SCF end and gradient start while ~3.5 GB per GPU were idle).
    // Stage 2 pseudo-diagonalisation state (Claude Generated, Sep 29, 2026).
    bool   pseudo_enabled = false;
    double pseudo_max_dq = 0.05;
    int    pseudo_max_steps = 8;
    int    pseudo_consecutive = 0;   // pseudo steps since the last full solve
    int    pseudo_taken = 0, pseudo_rejected = 0;
    double last_dq = 1e300;          // max|dq| of the previous resident step
    bool   in_resident_step = false; // pseudo steps only inside residentScfStep
    bool   cprev_valid = false;      // dCprev32 holds the previous step's eigenvectors (AO basis)
    CudaBuffer<float>  dCprev32, dPsT, dPsGov, dPsWork;
    CudaBuffer<double> dPsG, dPsStats;
    CudaBuffer<unsigned int> dPsMax;
    // FP64 pseudo steps (-scf_pseudo_diag_fp64): the same step in FP64 (stage 3 of the cost plan,
    // and the whole SCF on GPUs whose FP64 is not throttled).
    bool   pseudo_fp64 = false;
    bool   cprev64_valid = false;
    int    pseudo_consecutive64 = 0;
    float  pseudo_last_max_theta = 0.0f;
    bool   c_exact_eigenvectors = true;   // dC/dEps come from a full eigensolve (false after a pseudo step)
    CudaBuffer<double> dCprev64, dPsT64, dPsGov64, dPsWork64;
    void releasePseudoDiag()
    {
        dCprev32.free(); dPsT.free(); dPsGov.free(); dPsWork.free(); dPsG.free(); dPsStats.free(); dPsMax.free();
        dCprev64.free(); dPsT64.free(); dPsGov64.free(); dPsWork64.free();
        cprev_valid = cprev64_valid = false;
        pseudo_consecutive = pseudo_consecutive64 = 0;
    }
    bool pseudoDiagStep(int n);
    int  pseudoDiagStep64(int n);
    template <typename R>
    int pseudoRotateOrthonormalise(int n, int no, const R* F, const R* L, R* C, R* Tb, R* Gov, R* M,
                                   CudaBuffer<R>& work, double* g, const char* tag);

    void releaseEigenWorkspaces()
    {
        releasePseudoDiag();
        dWork.free(); lwork = 0;
        dCf.free(); dLf.free(); dWorkf.free(); lwork_f32 = 0;
        dProbeA.free(); dProbeV.free(); dProbeY.free(); dProbeT.free(); dProbeU.free();
        dProbeAd.free(); dProbeVd.free(); dProbeYd.free(); dProbeTd.free(); dProbeSd.free();
        releaseRemoteBackup();
        if (dist) { dist->releaseBuffers(); cudaSetDevice(device); }
    }

    /**
     * @brief Check that a distributed FP32 solve returned eigenpairs of the matrix it was given.
     *
     * A random vector v is pushed through the matrix twice: y = A v from a copy of the input
     * (symv, lower triangle) and y2 = Q (eps .* (Q^T v)) from what the solver returned. They must
     * agree. Three matrix-vector products plus one n x n copy, done ONCE per run: cusolverMg's
     * FP32 eigenvectors were wrong here (4x A4500, CUDA 13.3, n = 15444 - the SCF diverged from
     * the first iteration while the eigenvalue SUM still matched, so a trace check does not catch
     * it), but that is one library version on one machine, so it is measured rather than assumed.
     * @return true when the verdict was reached (dist_fp32_state set), false when it could not run
     */
    bool verifyDistributedFp32(int n, const float* Q, const float* eps_dev, const float* L = nullptr)
    {
        if (dProbeA.n < n || dProbeV.n < n) return false;
        const float one = 1.0f, zero = 0.0f, minus = -1.0f;
        const int b = 256;
        if (cublasSsymv(cublas, CUBLAS_FILL_MODE_LOWER, n, &one, dProbeA.ptr, n, dProbeV.ptr, 1,
                        &zero, dProbeY.ptr, 1) != CUBLAS_STATUS_SUCCESS)                 // y = A v
            return false;
        if (L) {
            // Claude Generated (Sep 29, 2026): the solver returned C = L^-T Q, not Q. With
            // Q = L^T C the same identity reads  A v = L^T C (eps .* (C^T (L v))), i.e. two extra
            // triangular matrix-vector products; the residual and its reference are unchanged.
            try { dProbeU.ensure(n); } catch (...) { return false; }
            if (cudaMemcpyAsync(dProbeU.ptr, dProbeV.ptr, sizeof(float) * n, cudaMemcpyDeviceToDevice,
                                stream) != cudaSuccess
                || cublasStrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, L, n,
                               dProbeU.ptr, 1) != CUBLAS_STATUS_SUCCESS                  // u = L v
                || cublasSgemv(cublas, CUBLAS_OP_T, n, n, &one, Q, n, dProbeU.ptr, 1, &zero,
                               dProbeT.ptr, 1) != CUBLAS_STATUS_SUCCESS)                 // t = C^T u
                return false;
            k_scale_by<<<(n + b - 1) / b, b, 0, stream>>>(dProbeT.ptr, eps_dev, n);      // t *= eps
            if (cudaGetLastError() != cudaSuccess
                || cublasSgemv(cublas, CUBLAS_OP_N, n, n, &one, Q, n, dProbeT.ptr, 1, &zero,
                               dProbeU.ptr, 1) != CUBLAS_STATUS_SUCCESS                  // u = C t
                || cublasStrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, L, n,
                               dProbeU.ptr, 1) != CUBLAS_STATUS_SUCCESS                  // u = L^T u
                || cublasSscal(cublas, n, &minus, dProbeY.ptr, 1) != CUBLAS_STATUS_SUCCESS
                || cublasSaxpy(cublas, n, &one, dProbeU.ptr, 1, dProbeY.ptr, 1)
                       != CUBLAS_STATUS_SUCCESS)                                         // y := u - y
                return false;
        } else {
            if (cublasSgemv(cublas, CUBLAS_OP_T, n, n, &one, Q, n, dProbeV.ptr, 1, &zero,
                            dProbeT.ptr, 1) != CUBLAS_STATUS_SUCCESS)                    // t = Q^T v
                return false;
            k_scale_by<<<(n + b - 1) / b, b, 0, stream>>>(dProbeT.ptr, eps_dev, n);      // t *= eps
            if (cudaGetLastError() != cudaSuccess) return false;
            if (cublasSgemv(cublas, CUBLAS_OP_N, n, n, &one, Q, n, dProbeT.ptr, 1, &minus,
                            dProbeY.ptr, 1) != CUBLAS_STATUS_SUCCESS)                    // y := Q t - y
                return false;
        }
        float res = 0.0f, ref = 0.0f;
        if (cublasSnrm2(cublas, n, dProbeY.ptr, 1, &res) != CUBLAS_STATUS_SUCCESS
            || cublasSnrm2(cublas, n, dProbeT.ptr, 1, &ref) != CUBLAS_STATUS_SUCCESS
            || cudaStreamSynchronize(stream) != cudaSuccess)
            return false;
        const double rel = (ref > 0.0f) ? static_cast<double>(res) / static_cast<double>(ref) : 1.0;
        // A = Q L Q^T alone is NOT enough: it is invariant under a consistent permutation of the
        // eigenpairs, while the SCF fills the LEADING columns and therefore needs the eigenvalues
        // ascending. Measured here (polymer, nao 3222, cusolverMg FP32): the residual check passed
        // at 1e-6 while the SCF diverged, which is exactly what a permuted spectrum looks like.
        std::vector<float> eps_host(n);
        if (cudaMemcpy(eps_host.data(), eps_dev, sizeof(float) * n, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        // Only count inversions that MATTER: in FP32 two nearly degenerate orbitals can come back
        // swapped by ~1e-7 Eh, which changes nothing for the occupation, while a genuinely permuted
        // spectrum shows up as large steps in the wrong direction. Measured: cuSOLVERMp at
        // nao = 15444 returns exactly one such harmless swap, and rejecting it cost 36 s.
        const double span = std::fabs(static_cast<double>(eps_host[n - 1]) - static_cast<double>(eps_host[0]));
        const double ord_tol = std::max(1.0e-6, 1.0e-6 * span);
        int inversions = 0;
        for (int i = 1; i < n; ++i)
            if (static_cast<double>(eps_host[i - 1]) - static_cast<double>(eps_host[i]) > ord_tol)
                ++inversions;
        // FP32 over n ~ 1e4 accumulations lands well below 1e-3 when the solve is correct; the
        // observed failure was O(1), so the threshold does not need to be tight.
        if (std::getenv("CURCUMA_GPU_EIG_VERIFY_ALWAYS"))
            std::fprintf(stderr, "[eig verify] solve %d: relative residual %.3e, %d inversions\n",
                         dist_solves, rel, inversions);
        dist_fp32_state = (rel < 1.0e-3 && inversions == 0) ? 1 : -1;
        if (dist_fp32_state < 0)
            dist_status = std::string(dist->name()) + ": its FP32 eigenpairs failed verification "
                "here (relative residual " + std::to_string(rel) + ", " + std::to_string(inversions)
                + " eigenvalues out of ascending order); FP32 iterations stay on device "
                + std::to_string(device) + " and only FP64 is distributed";
        else
            dist_status = std::string(dist->name()) + ": FP32 eigenpairs verified (relative "
                "residual " + std::to_string(rel) + ", spectrum ascending)";
        return true;
    }

    /**
     * @brief FP64 twin of verifyDistributedFp32, for both the plain and the generalized solve.
     *
     * plain:       A v  ==  C (eps .* (C^T v))
     * generalized: F v  ==  S C (eps .* (C^T (S v)))  with S = L L^T applied by two trmv calls,
     *              which is the identity F = S C diag(eps) C^T S for  F C = S C diag(eps),
     *              C^T S C = I.
     * dProbeAd holds the copy of A (plain) or F (generalized) taken before the solve.
     */
    bool verifyDistributedFp64(int n, const double* C, const double* eps_dev, bool generalized,
                               bool y_ready = false)
    {
        if ((!y_ready && dProbeAd.n < n) || dProbeVd.n < n || dProbeYd.n < n) return false;
        const double one = 1.0, zero = 0.0, minus = -1.0;
        const int b = 256;
        // y = M v   (M = A or F, symmetric, lower triangle); y_ready: computed before the solve
        if (!y_ready
            && cublasDsymv(cublas, CUBLAS_FILL_MODE_LOWER, n, &one, dProbeAd.ptr, n, dProbeVd.ptr, 1,
                           &zero, dProbeYd.ptr, 1) != CUBLAS_STATUS_SUCCESS)
            return false;
        // u = v (plain) or u = S v = L (L^T v) (generalized), kept in dProbeTd
        if (cudaMemcpyAsync(dProbeTd.ptr, dProbeVd.ptr, sizeof(double) * n,
                            cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
            return false;
        if (generalized) {
            if (cublasDtrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n,
                            dL.ptr, n, dProbeTd.ptr, 1) != CUBLAS_STATUS_SUCCESS
                || cublasDtrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n,
                               dL.ptr, n, dProbeTd.ptr, 1) != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        // t = C^T u ; t *= eps ; w = C t   (w reuses dProbeTd)
        // Scratch for C^T u; kept across calls with the other probe buffers.
        try { dProbeSd.ensure(n); } catch (...) { return false; }
        CudaBuffer<double>& tmp = dProbeSd;
        if (cublasDgemv(cublas, CUBLAS_OP_T, n, n, &one, C, n, dProbeTd.ptr, 1, &zero,
                        tmp.ptr, 1) != CUBLAS_STATUS_SUCCESS)
            return false;
        k_scale_byd<<<(n + b - 1) / b, b, 0, stream>>>(tmp.ptr, eps_dev, n);
        if (cudaGetLastError() != cudaSuccess) return false;
        if (cublasDgemv(cublas, CUBLAS_OP_N, n, n, &one, C, n, tmp.ptr, 1, &zero,
                        dProbeTd.ptr, 1) != CUBLAS_STATUS_SUCCESS)
            return false;
        if (generalized) {   // y2 = S w
            if (cublasDtrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n,
                            dL.ptr, n, dProbeTd.ptr, 1) != CUBLAS_STATUS_SUCCESS
                || cublasDtrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n,
                               dL.ptr, n, dProbeTd.ptr, 1) != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        double res = 0.0, ref = 0.0;
        if (cublasDaxpy(cublas, n, &minus, dProbeYd.ptr, 1, dProbeTd.ptr, 1) != CUBLAS_STATUS_SUCCESS
            || cublasDnrm2(cublas, n, dProbeTd.ptr, 1, &res) != CUBLAS_STATUS_SUCCESS
            || cublasDnrm2(cublas, n, dProbeYd.ptr, 1, &ref) != CUBLAS_STATUS_SUCCESS
            || cudaStreamSynchronize(stream) != cudaSuccess)
            return false;
        std::vector<double> eps_host(n);
        if (cudaMemcpy(eps_host.data(), eps_dev, sizeof(double) * n, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        // FP64 rounding is ~1e-16, so anything above a tiny tolerance is a real permutation.
        const double span = std::fabs(eps_host[n - 1] - eps_host[0]);
        const double ord_tol = std::max(1.0e-12, 1.0e-12 * span);
        int inversions = 0;
        for (int i = 1; i < n; ++i)
            if (eps_host[i - 1] - eps_host[i] > ord_tol) ++inversions;
        const double rel = (ref > 0.0) ? res / ref : 1.0;
        if (std::getenv("CURCUMA_GPU_EIG_VERIFY_ALWAYS"))
            std::fprintf(stderr, "[eig verify] FP64 solve %d: relative residual %.3e, %d inversions\n",
                         dist_solves, rel, inversions);
        // FP64 has ~1e-16 per operation, so a correct solve lands far below 1e-8 even at n = 15444.
        dist_fp64_state = (rel < 1.0e-8 && inversions == 0) ? 1 : -1;
        if (dist_fp64_state < 0)
            dist_status = std::string(dist->name()) + ": its FP64 eigenpairs failed verification "
                "here (relative residual " + std::to_string(rel) + ", " + std::to_string(inversions)
                + " eigenvalues out of ascending order); the eigensolve stays on device "
                + std::to_string(device);
        return true;
    }

    /// Buffer of `bytes` on the eigensolver device (other than this one) with the most free
    /// memory, kept across calls; false when none has room. Claude Generated (Sep 29, 2026).
    bool ensureRemoteBackup(size_t bytes)
    {
        if (std::getenv("CURCUMA_GPU_EIG_LOCAL_BACKUP")) return false;
        if (dist_backup && dist_backup_bytes >= bytes) return true;
        releaseRemoteBackup();
        int best = -1;
        size_t best_free = 0;
        for (int d : dist_devices) {
            if (d == device) continue;
            size_t fr = 0, tot = 0;
            cudaSetDevice(d);
            if (cudaMemGetInfo(&fr, &tot) == cudaSuccess && fr > best_free) { best_free = fr; best = d; }
        }
        bool ok = false;
        // keep 1 GB head room on the helper for its own solver buffers
        if (best >= 0 && best_free > bytes + (size_t(1) << 30)) {
            cudaSetDevice(best);
            if (cudaMalloc(&dist_backup, bytes) == cudaSuccess) {
                dist_backup_dev = best;
                dist_backup_bytes = bytes;
                ok = true;
            } else {
                cudaGetLastError();
                dist_backup = nullptr;
            }
        }
        cudaSetDevice(device);
        return ok;
    }

    void releaseRemoteBackup()
    {
        if (dist_backup) {
            cudaSetDevice(dist_backup_dev);
            cudaFree(dist_backup);
            cudaSetDevice(device);
        }
        dist_backup = nullptr;
        dist_backup_dev = -1;
        dist_backup_bytes = 0;
    }

    /// Create the multi-GPU solver on first use; false when it is not to be used for this call.
    bool distReady(int n, bool fp32)
    {
        if (dist_devices.size() < 2 || dist_failed || n < dist_min_nao || (fp32 && !dist_fp32))
            return false;
        if (!dist) {
            std::string why;
            dist = DistributedEigensolver::create(dist_backend, dist_devices, dist_block, why);
            cudaSetDevice(device);
            if (!dist) {
                dist_failed = true;
                dist_status = "unavailable: " + why;
                return false;
            }
        }
        // Claude Generated (Sep 2026, docs/MULTI_GPU_GAPS.md G2-8): cusolverMg is FP64-only by
        // design (its FP32 eigenpairs degraded after the first call, see the verification note
        // above). It used to be excluded only by the per-solve check, i.e. not at all with
        // -gpu_eigensolver_verify false. FP32 iterations stay on this device.
        if (fp32 && std::string(dist->name()) == "cusolverMg")
            return false;
        // FP32 correctness is a per-machine question (see dist_fp32_state): the first FP32 solve
        // is verified, and only a failed verification turns FP32 off for the rest of the run.
        // the wrapper reports dist_status; a backend that failed verification is not used again
        if (fp32 && dist_fp32_state < 0) return false;
        if (!fp32 && dist_fp64_state < 0) return false;
        return true;
    }

    /**
     * @brief Cholesky of the overlap in dL (S in the lower triangle) on the eigensolver's devices
     *        (Claude Generated, Sep 29, 2026).
     *
     * Same gates as the FP64 eigensolve (distReady), plus a check of the result: a copy of S is
     * kept, and S v == L (L^T v) is tested with a random v (relative residual < 1e-10; a correct
     * FP64 factorisation lands near 1e-15). On a failed check S is put back and the caller runs
     * the single-GPU potrf. The copy is one n x n buffer during setup, far below the SCF peak.
     * CURCUMA_GPU_DIST_POTRF=0 keeps the factorisation on this device.
     * @return 1 = dL holds L, 0 = not done (dL still holds S), -1 = dL destroyed
     */
    int distCholesky(int n)
    {
        if (!distReady(n, false) || !dist->supportsCholesky()) return 0;
        if (const char* e = std::getenv("CURCUMA_GPU_DIST_POTRF"); e && std::string(e) == "0") return 0;
        const size_t nn = static_cast<size_t>(n) * n;
        CudaBuffer<double> copy, v, y, t;
        try { copy.alloc(static_cast<int>(nn)); v.alloc(n); y.alloc(n); t.alloc(n); }
        catch (...) { return 0; }
        const int b = 256;
        k_probe_filld<<<(n + b - 1) / b, b, 0, stream>>>(v.ptr, n, 0x27D4EB2Fu);
        if (cudaGetLastError() != cudaSuccess
            || cudaMemcpyAsync(copy.ptr, dL.ptr, sizeof(double) * nn, cudaMemcpyDeviceToDevice, stream) != cudaSuccess
            || cudaStreamSynchronize(stream) != cudaSuccess)
            return 0;
        const bool done = dist->cholesky(n, dL.ptr, device);
        cudaSetDevice(device);
        if (prof) {
            double sc = 0.0, so = 0.0, ga = 0.0;
            dist->lastTimings(sc, so, ga);
            profAdd("  multi-GPU potrf: scatter", sc);
            profAdd("  multi-GPU potrf: factorise", so);
            profAdd("  multi-GPU potrf: gather", ga);
        }
        bool good = false;
        if (done) {
            const double one = 1.0, zero = 0.0, minus = -1.0;
            double res = 0.0, ref = 0.0;
            good = cublasDsymv(cublas, CUBLAS_FILL_MODE_LOWER, n, &one, copy.ptr, n, v.ptr, 1, &zero,
                               y.ptr, 1) == CUBLAS_STATUS_SUCCESS                       // y = S v
                && cudaMemcpyAsync(t.ptr, v.ptr, sizeof(double) * n, cudaMemcpyDeviceToDevice, stream) == cudaSuccess
                && cublasDtrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, dL.ptr, n,
                               t.ptr, 1) == CUBLAS_STATUS_SUCCESS                         // t = L^T v
                && cublasDtrmv(cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, dL.ptr, n,
                               t.ptr, 1) == CUBLAS_STATUS_SUCCESS                         // t = L t
                && cublasDnrm2(cublas, n, y.ptr, 1, &ref) == CUBLAS_STATUS_SUCCESS
                && cublasDaxpy(cublas, n, &minus, y.ptr, 1, t.ptr, 1) == CUBLAS_STATUS_SUCCESS
                && cublasDnrm2(cublas, n, t.ptr, 1, &res) == CUBLAS_STATUS_SUCCESS
                && cudaStreamSynchronize(stream) == cudaSuccess
                && ref > 0.0 && res / ref < 1.0e-10;
            if (std::getenv("CURCUMA_GPU_EIG_VERIFY_ALWAYS"))
                std::fprintf(stderr, "[chol verify] relative residual %.3e\n", ref > 0.0 ? res / ref : -1.0);
        }
        if (good) {
            ++dist_potrf;
            return 1;
        }
        // Not factorised or failed the check: put S back for the single-GPU potrf.
        if (cudaMemcpyAsync(dL.ptr, copy.ptr, sizeof(double) * nn, cudaMemcpyDeviceToDevice, stream) != cudaSuccess
            || cudaStreamSynchronize(stream) != cudaSuccess)
            return -1;
        if (dist_potrf_status.empty())
            dist_potrf_status = done ? "multi-GPU Cholesky of S failed its check; factorised on device "
                                           + std::to_string(device)
                                     : "multi-GPU Cholesky of S failed; factorised on device " + std::to_string(device);
        return 0;
    }

    /// Generalized variant (reduction + solve + back-transform on the devices). Same return codes
    /// as distSolve(); 0 also when the backend has no generalized path.
    int distSolveGeneralized(int n, void* A, const void* L, void* eig, bool fp32)
    {
        if (!distReady(n, fp32) || !dist->supportsGeneralized()) return 0;
        if (cudaStreamSynchronize(stream) != cudaSuccess) return 0;
        // Same verification as the plain path, on the generalized residual F c = eps S c with
        // S = L L^T. A backend that passes once is trusted for the rest of the run; one that fails
        // is dropped and the solve returns to this device. CURCUMA_GPU_EIG_VERIFY_ALWAYS=1 checks
        // every call - which is how the cusolverMg reuse defect was found (first call exact, second
        // call wrong).
        bool verify64 = !fp32 && dist_verify && dist_fp64_state >= 0;
        // Claude Generated (Sep 29, 2026): the verification needs y = F v and, only to redo a
        // rejected solve, a copy of F. y is now taken before the solve (F is still intact), and
        // the copy lives on one of the helper devices when it fits there: it was the largest
        // single buffer on this device during the FP64 solve (1.9 GB at n = 15444; with 2 GPUs the
        // device peaked 2.3 GB above the single-GPU run). On a helper-device failure the copy
        // stays here as before.
        bool backup_remote = false;
        if (verify64) {
            try {
                dProbeVd.ensure(n); dProbeYd.ensure(n); dProbeTd.ensure(n);
            } catch (...) {
                dProbeVd.free(); dProbeYd.free(); dProbeTd.free();
                verify64 = false;
            }
        }
        if (verify64) {
            const size_t bytes = sizeof(double) * static_cast<size_t>(n) * n;
            backup_remote = ensureRemoteBackup(bytes)
                && cudaMemcpyPeer(dist_backup, dist_backup_dev, A, device, bytes) == cudaSuccess;
            cudaSetDevice(device);
            if (!backup_remote) {
                try {
                    dProbeAd.ensure(static_cast<int>(static_cast<size_t>(n) * n));
                } catch (...) {
                    dProbeAd.free(); dProbeVd.free(); dProbeYd.free(); dProbeTd.free();
                    dist_fp64_state = -1;
                    dist_status = std::string(dist->name()) + ": FP64 verification needs one more n x n "
                        "buffer than fits on device " + std::to_string(device)
                        + "; the eigensolve stays on that device";
                    return 0;
                }
            }
            const double one = 1.0, zero = 0.0;
            const int b = 256;
            k_probe_filld<<<(n + b - 1) / b, b, 0, stream>>>(dProbeVd.ptr, n, 0x85EBCA6Bu);
            if (cudaGetLastError() != cudaSuccess
                || cublasDsymv(cublas, CUBLAS_FILL_MODE_LOWER, n, &one, static_cast<const double*>(A), n,
                               dProbeVd.ptr, 1, &zero, dProbeYd.ptr, 1) != CUBLAS_STATUS_SUCCESS   // y = F v
                || (!backup_remote
                    && cudaMemcpyAsync(dProbeAd.ptr, A, bytes, cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
                || cudaStreamSynchronize(stream) != cudaSuccess)
                verify64 = false;
        }
        const bool solved = dist->solveGeneralized(n, A, L, l_generation, eig, fp32, device);
        cudaSetDevice(device);
        if (prof) {
            double sc = 0.0, so = 0.0, ga = 0.0, re = 0.0, ba = 0.0;
            dist->lastTimings(sc, so, ga);
            dist->lastGeneralizedTimings(re, ba);
            const std::string p = fp32 ? "  multi-GPU FP32: " : "  multi-GPU FP64: ";
            profAdd((p + "scatter F (+L)").c_str(), sc);
            profAdd((p + "sygst").c_str(), re);
            profAdd((p + "syevd").c_str(), so - re - ba);
            profAdd((p + "trsm back-transform").c_str(), ba);
            profAdd((p + "gather C").c_str(), ga);
        }
        if (solved) {
            ++dist_solves;
            bool rejected = false;
            if (verify64 && verifyDistributedFp64(n, static_cast<const double*>(A),
                                                  static_cast<const double*>(eig), true, true)
                && dist_fp64_state < 0) {
                const size_t bytes = sizeof(double) * static_cast<size_t>(n) * n;
                rejected = backup_remote
                    ? cudaMemcpyPeer(A, device, dist_backup, dist_backup_dev, bytes) == cudaSuccess
                    : (cudaMemcpyAsync(A, dProbeAd.ptr, bytes, cudaMemcpyDeviceToDevice, stream) == cudaSuccess
                       && cudaStreamSynchronize(stream) == cudaSuccess);
                cudaSetDevice(device);
                --dist_solves;
            }
            if (rejected) return 0;
            if (dist_fp64_state < 0) return -1;
            return 1;
        }
        dist_failed = true;
        dist_status = std::string(dist->name()) + " generalized solve failed (n = " + std::to_string(n)
            + (fp32 ? ", FP32" : ", FP64") + "); single-GPU solver from here on";
        return dist->inputIntact() ? 0 : -1;
    }

    int    sparse_mode = 1;          // 0 = always dense, 1 = auto, 2 = always sparse
    bool   sparse = false;           // storage used for the current geometry
    double sparse_eps = 1.0e-20;     // integrals below this are dropped
    int    sp_nnz = 0;
    double sp_fraction = 1.0;        // nnz / nao^2
    double sp_rmax = 0.0;            // largest element-pair cutoff (Bohr)
    CudaBuffer<int>    dSpRow, dSpCol, dSpColPtr, dSpPerm;
    // Claude Generated (Sep 29, 2026): atom-pair screen for the dense S/H0 build in screened
    // mode. Element kind per atom and the squared element-pair cutoff (the same table the pair
    // pattern is built from, slightly widened), so k_overlap_h0 skips the atom pairs whose
    // integrals the SCF never stores. screen_nk == 0: no screen (dense mode, or disabled).
    CudaBuffer<int>    dAtKind;
    CudaBuffer<double> dScreenCut2;
    int                screen_nk = 0;
    CudaBuffer<double> dSpS, dSpH0, dSpDp, dSpQp, dSpTmp;
    // Density: the resident screened loop keeps P only on the pattern (dSpP) and rebuilds the
    // dense dP on demand (gradient, host download). p_dense_valid says which one is current.
    CudaBuffer<double> dSpP;
    bool               p_dense_valid = true;
    int                last_ncol = 0;
    // Point 2 (Claude Generated, Sep 2026): the device gradient's H0-Pulay and CN-onsite
    // kernels read P and W ONLY at stored screened pairs (verified exhaustively - see
    // computeGradient) and at the AO diagonal (always in-pattern by construction, see
    // buildScreenedPairs). dSpW is the pattern-only twin of dSpP for the energy-weighted
    // density W (built by the *_spP gradient kernels below via k_density_sp, weight 2*eps
    // instead of occ); dSpDiagIdx[mu] is the pair-storage index e of the (mu,mu) entry,
    // computed once in buildScreenedPairs (same loop that builds dSpRow/dSpCol, no extra
    // scan). Neither is needed when the geometry is small enough to stay dense.
    CudaBuffer<double> dSpW;
    CudaBuffer<int>    dSpDiagIdx;
    std::vector<int>   h_sp_row, h_sp_col;       // host copy (dense downloads)
    // Per-atom screening data captured in beginBasis (geometry independent).
    std::vector<int>    h_z, h_at_ao0, h_at_nao;
    std::vector<double> h_at_amin, h_at_cmax;
    bool                h_ao_contiguous = false;

    /// Free every nao²- and nat²-sized buffer. Called when a basis cannot be set up, so a
    /// refused or half-allocated basis does not keep gigabytes of device memory pinned for the
    /// rest of the process (the CPU fallback then runs with a clean device).
    void releaseLarge()
    {
        for (CudaBuffer<double>* b : { &dH0, &dS, &dL, &dC, &dP, &dCw, &dWork, &dDpInt, &dQpInt,
                                       &dSdR, &dPotrfWork, &dW, &dGamma, &dMpAmatSD, &dMpAmatDD,
                                       &dMpAmatSQ, &dSolvB, &dEeqM, &dEeqWork })
            b->free();
        for (CudaBuffer<float>* b : { &dCf, &dLf, &dWorkf })
            b->free();
        for (CudaBuffer<double>* b : { &dSpS, &dSpH0, &dSpDp, &dSpQp, &dSpTmp, &dSpP, &dSpW })
            b->free();
        releaseDensityHelpers();
        mp_otf = false;
        p_dense_valid = true;
        for (CudaBuffer<int>* b : { &dSpRow, &dSpCol, &dSpColPtr, &dSpPerm, &dSpDiagIdx })
            b->free();
        h_sp_row.clear(); h_sp_row.shrink_to_fit();
        h_sp_col.clear(); h_sp_col.shrink_to_fit();
        sparse = false;
        sp_nnz = 0;
        lwork = 0;
        lwork_f32 = 0;
        potrf_lwork = 0;
        resident_n = 0;
        basis_nao = 0;
    }
};

// ---- Stage 3 element parameter tables in __constant__ memory --------------
// Only the D3 covalent radii (au) are needed for the CN kernel; more tables
// (hubbard, shell_hubbard, pauling_en, …) join here in 3b/3c/3d. The 86-element
// arrays are far under the 64 KB constant limit. Seeded once via ensureStage3Constants.
namespace {
__constant__ double c_covrad_d3[86];   // covalent_rad_d3_au(z), index z-1
__constant__ double c_pauling[86];     // pauling_en[z-1]
__constant__ double c_atomic_rad[86];  // atomic_rad_au(z), index z-1

// Claude Generated (Sep 2026, multi-GPU): __constant__ memory belongs to ONE device, so the
// tables must be uploaded once PER DEVICE, not once per process. The old `static bool`
// left every device but the first with uninitialised tables (wrong CN/EEQ, no error)
// and was a data race when several workers built contexts at the same time.
void ensureStage3Constants()
{
    static std::mutex mtx;
    static std::set<int> loaded_devices;
    int dev = 0;
    if (cudaGetDevice(&dev) != cudaSuccess) return;
    std::lock_guard<std::mutex> lock(mtx);
    if (loaded_devices.count(dev)) return;
    double covrad[86], pauling[86], arad[86];
    for (int z = 1; z <= 86; ++z) {
        covrad[z - 1]  = curcuma::xtb::covalent_rad_d3_au(z);
        pauling[z - 1] = curcuma::xtb::pauling_en[z - 1];
        arad[z - 1]    = curcuma::xtb::atomic_rad_au(z);
    }
    cudaMemcpyToSymbol(c_covrad_d3, covrad, sizeof(covrad));
    cudaMemcpyToSymbol(c_pauling, pauling, sizeof(pauling));
    cudaMemcpyToSymbol(c_atomic_rad, arad, sizeof(arad));
    loaded_devices.insert(dev);
}
} // namespace

// ---- Stage 3 device kernels ----------------------------------------------

// Coordination numbers: one thread per atom i sums c_ij over all j≠i. Mirrors
// cn_exp (GFN1) / cn_gfn (GFN2) (xtb_params_extra.hpp:153/183). The CPU loops
// unique pairs j<i and adds to both; summing all j per atom is identical.
// cutoff = 25 au, r²<1e-12 guards coincident atoms.
__global__ void k_cn(int nat, const double* __restrict__ xyz, const int* __restrict__ z,
                     int is_gfn2, double* __restrict__ cn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const double xi = xyz[3 * i + 0], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
    const double rci = c_covrad_d3[z[i] - 1];
    const double cutoff2 = 25.0 * 25.0;
    double sum = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const double dx = xi - xyz[3 * j + 0];
        const double dy = yi - xyz[3 * j + 1];
        const double dz = zi - xyz[3 * j + 2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cutoff2 || r2 < 1.0e-12) continue;
        const double r  = sqrt(r2);
        const double rc = rci + c_covrad_d3[z[j] - 1];
        sum += d_cn_pair(r, rc, is_gfn2 != 0);
    }
    cn[i] = sum;
}

// CN-shifted shell self-energies: se[s] = selfenergy[s] − kcn[s]·CN[sh2at[s]]
// (mirrors XTB::getSelfEnergies, xtb_h0.cpp:79).
__global__ void k_self_energy(int nsh, const double* __restrict__ selfenergy,
                              const double* __restrict__ kcn, const int* __restrict__ sh2at,
                              const double* __restrict__ cn, double* __restrict__ se)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    se[s] = selfenergy[s] - kcn[s] * cn[sh2at[s]];
}

// Overlap S and bare Hamiltonian H0: one thread per shell-pair (ish_a, ish_b)
// computes the shell-pair h_factor once, then loops AO pairs writing S and H0.
// Mirrors XTB::getHamiltonianH0 (xtb_h0.cpp:140-231) index-for-index, incl. the
// on-atom same-orbital S=1 special case and the {py,pz,px} AO ordering. Outputs
// column-major nao×nao (S, H0 symmetric → layout-safe). se is the CN-shifted
// self-energy from k_self_energy.
__global__ void k_overlap_h0(
    int nsh, int nao, int is_gfn2,
    const int* __restrict__ sh2at, const int* __restrict__ ang_sh,
    const int* __restrict__ iao_sh, const int* __restrict__ nao_sh,
    const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly,
    const double* __restrict__ se, const int* __restrict__ z,
    const int* __restrict__ valence, const double* __restrict__ xyz,
    double* __restrict__ S, double* __restrict__ H0,
    const int* __restrict__ at_kind, const double* __restrict__ cut2, int nk)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    const int b = blockIdx.y * blockDim.y + threadIdx.y;
    if (a >= nsh || b >= nsh) return;

    const int iat = sh2at[a], jat = sh2at[b];
    const int la = ang_sh[a], lb = ang_sh[b];
    const double xa = xyz[3 * iat + 0], ya = xyz[3 * iat + 1], za = xyz[3 * iat + 2];
    const double xb = xyz[3 * jat + 0], yb = xyz[3 * jat + 1], zb = xyz[3 * jat + 2];
    // Claude Generated (Sep 29, 2026): screened storage keeps only the atom pairs inside the
    // element-pair cutoff (every integral beyond it is below sparse_eps). The caller zeroed S and
    // H0, so a pair outside the (slightly widened) cutoff is skipped instead of computed.
    if (nk > 0 && iat != jat) {
        const double sx = xa - xb, sy = ya - yb, sz = za - zb;
        if (sx * sx + sy * sy + sz * sz > cut2[at_kind[iat] * nk + at_kind[jat]]) return;
    }
    const double avg_eps = 0.5 * (se[a] + se[b]);

    double h_factor;
    if (iat == jat) {
        h_factor = 1.0;
    } else {
        const int zi = z[iat], zj = z[jat];
        const double dx = xa - xb, dy = ya - yb, dz = za - zb;
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double rr = sqrt(sqrt(r2) / (c_atomic_rad[zi - 1] + c_atomic_rad[zj - 1]));
        const double pi_ij = (1.0 + shpoly[a] * rr) * (1.0 + shpoly[b] * rr);
        double hs;
        if (is_gfn2 == 0) {
            const bool vi = valence[a] != 0, vj = valence[b] != 0;
            if (vi && vj) {
                double den = c_pauling[zi - 1] - c_pauling[zj - 1]; den *= den;
                hs = d_kpair_gfn1(zi, zj) * d_kshell_gfn1(la, lb) * (1.0 + (-7.0e-3) * den);
            } else if (vi && !vj) {
                hs = 0.5 * (d_kshell_gfn1(la, la) + 2.85);
            } else if (!vi && vj) {
                hs = 0.5 * (d_kshell_gfn1(lb, lb) + 2.85);
            } else {
                hs = 2.85;
            }
        } else {
            double den = c_pauling[zi - 1] - c_pauling[zj - 1]; den *= den;
            const double enp = 1.0 + 2.0e-2 * den;
            const double km  = d_kshell_gfn2(la, lb) * enp;  // kpair=1 for GFN2
            const double za_ = sh_zeta[a], zb_ = sh_zeta[b];
            const double zij = pow(2.0 * sqrt(za_ * zb_) / (za_ + zb_), 0.5);  // wexp=0.5
            hs = zij * km;
        }
        h_factor = hs * pi_ij;
    }

    const int ia_start = iao_sh[a], ia_nao = nao_sh[a];
    const int jb_start = iao_sh[b], jb_nao = nao_sh[b];
    const double* pa_alpha = prim_alpha + sh_prim_off[a];
    const double* pa_coeff = prim_coeff + sh_prim_off[a];
    const int npa = sh_nprim[a];
    const double* pb_alpha = prim_alpha + sh_prim_off[b];
    const double* pb_coeff = prim_coeff + sh_prim_off[b];
    const int npb = sh_nprim[b];

    // X-I1: d-touching shell pairs use the cartesian->spherical element function;
    // pure s/p pairs keep the scalar d_cgto_overlap path (byte-identical).
    const bool dpair = (la >= 2 || lb >= 2);
    for (int ia = 0; ia < ia_nao; ++ia) {
        const int mu = ia_start + ia;
        for (int jb = 0; jb < jb_nao; ++jb) {
            const int nu = jb_start + jb;
            double s_ab;
            if (!dpair) {
                const int ta = d_ao_to_type(la, ia);
                const int tb = d_ao_to_type(lb, jb);
                if (ta < 0 || tb < 0) continue;
                s_ab = (iat == jat && a == b && ta == tb)
                    ? 1.0
                    : d_cgto_overlap(pa_alpha, pa_coeff, npa, pb_alpha, pb_coeff, npb,
                                     xa, ya, za, xb, yb, zb, ta, tb);
            } else {
                s_ab = d_overlap_elem(la, ia, lb, jb, pa_alpha, pa_coeff, npa,
                                      pb_alpha, pb_coeff, npb, xa, ya, za, xb, yb, zb);
                if (iat == jat && a == b && ia == jb) s_ab = 1.0;
            }
            const size_t idx = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
            S[idx]  = s_ab;
            H0[idx] = avg_eps * h_factor * s_ab;
        }
    }
}

// Coulomb γ matrix: one thread per shell-pair (is, js). Mirrors
// build_gamma_matrix (xtb_coulomb.hpp:98), gexp=2 for both methods, per-shell
// hardness g precomputed on the host. Column-major nsh×nsh (symmetric).
__global__ void k_gamma(int nsh, int is_gfn2, const int* __restrict__ sh2at,
                        const double* __restrict__ g, const double* __restrict__ xyz,
                        double* __restrict__ gamma)
{
    const int is = blockIdx.x * blockDim.x + threadIdx.x;
    const int js = blockIdx.y * blockDim.y + threadIdx.y;
    if (is >= nsh || js >= nsh) return;
    const int iat = sh2at[is], jat = sh2at[js];
    double val;
    if (is == js) {
        val = g[is];                                  // diagonal: raw hardness
    } else if (iat == jat) {
        val = d_coulomb_average(g[is], g[js], is_gfn2 != 0);  // on-atom cross-shell
    } else {
        const double dx = xyz[3 * iat + 0] - xyz[3 * jat + 0];
        const double dy = xyz[3 * iat + 1] - xyz[3 * jat + 1];
        const double dz = xyz[3 * iat + 2] - xyz[3 * jat + 2];
        const double r1 = sqrt(dx * dx + dy * dy + dz * dz);
        const double r1g = pow(r1, 2.0);              // gexp=2
        const double gam = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
        val = pow(r1g + pow(gam, -2.0), -0.5);        // (R² + γ̄⁻²)^(−1/2)
    }
    gamma[static_cast<size_t>(is) + static_cast<size_t>(js) * nsh] = val;
}

// GFN2 multipole integrals: one thread per AO pair (mu, nu). Computes the
// global-origin dipole/quadrupole (d_cgto_multipole), then the per-column origin
// shift (origin = atom of the column AO nu) + traceless transform using the
// resident overlap S. Mirrors setupMultipole (xtb_multipole.cpp:84-159).
// dp_int contiguous 3·nn, qp_int 6·nn, column-major (mu,nu) at mu+nu*nao —
// matching k_add_fock_multipole / k_multipole_moments. dp/qp NOT symmetric.
__global__ void k_multipole_ints(
    int nao, const int* __restrict__ ao2sh, const int* __restrict__ ao2at,
    const int* __restrict__ iao_sh, const int* __restrict__ ang_sh,
    const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ xyz, const double* __restrict__ S,
    double* __restrict__ dp_int, double* __restrict__ qp_int)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu >= nao || nu >= nao) return;
    const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
    const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;

    const int isha = ao2sh[mu], iat = ao2at[mu];
    const int ishb = ao2sh[nu], jat = ao2at[nu];
    const int la = ang_sh[isha], lb = ang_sh[ishb];
    const int sa = mu - iao_sh[isha], sb = nu - iao_sh[ishb];
    const bool dpair = (la >= 2 || lb >= 2);   // X-I1

    for (int k = 0; k < 3; ++k) dp_int[static_cast<size_t>(k) * nn + mn] = 0.0;
    for (int k = 0; k < 6; ++k) qp_int[static_cast<size_t>(k) * nn + mn] = 0.0;

    const double* aA = prim_alpha + sh_prim_off[isha];
    const double* cA = prim_coeff + sh_prim_off[isha];
    const int npa = sh_nprim[isha];
    const double* aB = prim_alpha + sh_prim_off[ishb];
    const double* cB = prim_coeff + sh_prim_off[ishb];
    const int npb = sh_nprim[ishb];

    double Sx, D[3], Q[6];
    if (!dpair) {
        const int ta = d_ao_to_type(la, sa);
        const int tb = d_ao_to_type(lb, sb);
        if (ta < 0 || tb < 0) return;
        d_cgto_multipole(aA, cA, npa, aB, cB, npb,
                         xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                         xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], ta, tb, Sx, D, Q);
    } else {
        d_multipole_elem(la, sa, lb, sb, aA, cA, npa, aB, cB, npb,
                         xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                         xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], Sx, D, Q);
    }

    // Origin shift to the column atom (nu) + traceless transform with overlap S.
    const double Rx = xyz[3*jat+0], Ry = xyz[3*jat+1], Rz = xyz[3*jat+2];
    const double Smn = S[mn];
    const double dx = D[0], dy = D[1], dz = D[2];
    dp_int[0*nn + mn] = dx - Rx * Smn;
    dp_int[1*nn + mn] = dy - Ry * Smn;
    dp_int[2*nn + mn] = dz - Rz * Smn;

    const double qxx = Q[0] - 2*Rx*dx + Rx*Rx*Smn;
    const double qxy = Q[1] - Rx*dy - Ry*dx + Rx*Ry*Smn;
    const double qyy = Q[2] - 2*Ry*dy + Ry*Ry*Smn;
    const double qxz = Q[3] - Rx*dz - Rz*dx + Rx*Rz*Smn;
    const double qyz = Q[4] - Ry*dz - Rz*dy + Ry*Rz*Smn;
    const double qzz = Q[5] - 2*Rz*dz + Rz*Rz*Smn;
    const double tr = 0.5 * (qxx + qyy + qzz);
    qp_int[0*nn + mn] = 1.5 * qxx - tr;
    qp_int[1*nn + mn] = 1.5 * qxy;
    qp_int[2*nn + mn] = 1.5 * qyy - tr;
    qp_int[3*nn + mn] = 1.5 * qxz;
    qp_int[4*nn + mn] = 1.5 * qyz;
    qp_int[5*nn + mn] = 1.5 * qzz - tr;
}

// Overlap derivative dS_μν/dR_{atom(μ)} for every AO pair: one thread per (μ,ν),
// 3 components written to dSdR (contiguous 3·nn, column-major (μ,ν) at mu+nu*nao).
// Crux primitive of the Stage-4 H0/Pulay gradient; validated standalone first.
__global__ void k_overlap_grad(
    int nao, const int* __restrict__ ao2sh, const int* __restrict__ ao2at,
    const int* __restrict__ iao_sh, const int* __restrict__ ang_sh,
    const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ xyz, double* __restrict__ dSdR)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu >= nao || nu >= nao) return;
    const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
    const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;

    const int isha = ao2sh[mu], iat = ao2at[mu];
    const int ishb = ao2sh[nu], jat = ao2at[nu];
    const int la = ang_sh[isha], lb = ang_sh[ishb];      // X-I1
    const int sa = mu - iao_sh[isha], sb = nu - iao_sh[ishb];

    double g[3] = {0.0, 0.0, 0.0};
    if (la >= 2 || lb >= 2) {
        d_overlap_grad_elem(la, sa, lb, sb,
            prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
            prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
            xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
            xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], g);
    } else {
        const int ta = d_ao_to_type(la, sa);
        const int tb = d_ao_to_type(lb, sb);
        if (ta >= 0 && tb >= 0) {
            d_cgto_overlap_grad(
                prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
                prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
                xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], ta, tb, g);
        }
    }
    dSdR[0*nn + mn] = g[0];
    dSdR[1*nn + mn] = g[1];
    dSdR[2*nn + mn] = g[2];
}

// ---- Stage 4 gradient kernels (grad layout [3*i+k], Eh/Bohr) ---------------

// Section 1: repulsion gradient. One thread per atom i, inner loop j<i.
__global__ void k_grad_repulsion(int nat, int is_gfn2, const int* __restrict__ z,
                                 const double* __restrict__ xyz, const double* __restrict__ rep_alpha,
                                 const double* __restrict__ rep_zeff, double kexp, double rexp,
                                 double kexp_light, double* __restrict__ grad)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int zi = z[i];
    const double alfi = rep_alpha[i], zeffi = rep_zeff[i];
    const double xi = xyz[3*i+0], yi = xyz[3*i+1], zit = xyz[3*i+2];
    for (int j = 0; j < i; ++j) {
        const double dx = xi - xyz[3*j+0], dy = yi - xyz[3*j+1], dz = zit - xyz[3*j+2];
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1.0e-12) continue;
        const double r = sqrt(r2);
        const double alpha_pair = sqrt(alfi * rep_alpha[j]);
        double kexp_pair = kexp;
        if (is_gfn2 && zi <= 2 && z[j] <= 2) kexp_pair = kexp_light;
        const double r_kexp = pow(r, kexp_pair);
        const double E_pair = zeffi * rep_zeff[j] / pow(r, rexp) * exp(-alpha_pair * r_kexp);
        const double dEdr = -(rexp / r + alpha_pair * kexp_pair * pow(r, kexp_pair - 1.0)) * E_pair;
        const double fx = dEdr*dx/r, fy = dEdr*dy/r, fz = dEdr*dz/r;
        atomicAdd(&grad[3*i+0], fx); atomicAdd(&grad[3*j+0], -fx);
        atomicAdd(&grad[3*i+1], fy); atomicAdd(&grad[3*j+1], -fy);
        atomicAdd(&grad[3*i+2], fz); atomicAdd(&grad[3*j+2], -fz);
    }
}

// Section 2a: on-site CN coupling dEdcn[iat] += −kcn[ish]·P(μ,μ). One thread per AO μ.
__global__ void k_grad_cn_onsite(int nao, const double* __restrict__ P,
                                 const double* __restrict__ kcn, const int* __restrict__ ao2sh,
                                 const int* __restrict__ ao2at, double* __restrict__ dEdcn)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    const double Pmm = P[static_cast<size_t>(mu) + static_cast<size_t>(mu) * nao];
    atomicAdd(&dEdcn[ao2at[mu]], (-kcn[ao2sh[mu]]) * Pmm);
}

// Point 2 (Claude Generated, Sep 2026): pattern-only twin - P(mu,mu) read straight out of the
// compact dSpP array via the diagonal index dSpDiagIdx[mu] (built once in buildScreenedPairs;
// the AO diagonal is always in the screened pattern - an atom is always its own neighbour).
__global__ void k_grad_cn_onsite_spP(int nao, const double* __restrict__ Psp,
                                     const int* __restrict__ diag_idx,
                                     const double* __restrict__ kcn, const int* __restrict__ ao2sh,
                                     const int* __restrict__ ao2at, double* __restrict__ dEdcn)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    const double Pmm = Psp[diag_idx[mu]];
    atomicAdd(&dEdcn[ao2at[mu]], (-kcn[ao2sh[mu]]) * Pmm);
}

// Section 2b: H0/Pulay off-site gradient. One thread per AO pair (μ,ν) with iat<jat.
// Mirrors xtb_gradient.cpp:228-438 (GFN1/GFN2 isotropic part; the GFN2 multipole
// integral Pulay block is Stage 4b and handled separately).
// Pair body of the H0/Pulay gradient, shared by the dense kernel (every (mu,nu)), the
// screened-pair kernel (stored pairs only; S and H0 come from the sparse arrays, P/W from
// the dense buffers at the pattern position) and the pattern-only screened-pair kernel
// (S/H0/P/W all come from the sparse arrays - Point 2, Sep 2026).
// Claude Generated (Sep 2026): factored out of k_grad_h0_pulay unchanged. Pmn/Wmn are now
// passed in by the caller (was a dense P[mn]/W[mn] lookup here) so the same body serves a
// dense buffer, a dense buffer read at a pattern index, or a pattern-only array - the caller
// decides where Pmn/Wmn come from, this function no longer knows or cares.
__device__ __forceinline__ void d_grad_h0_pulay_pair(
    int mu, int nu, double Smn, double H0mn, double Pmn, double Wmn,
    int is_gfn2,
    const int* __restrict__ ao2sh, const int* __restrict__ ao2at, const int* __restrict__ ang,
    const int* __restrict__ iao_sh, const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly, const double* __restrict__ kcn,
    const int* __restrict__ valence, const int* __restrict__ z, const double* __restrict__ se,
    const double* __restrict__ xyz, const double* __restrict__ v_ao,
    const double* __restrict__ v_dp, const double* __restrict__ v_qp,
    double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int iat = ao2at[mu], jat = ao2at[nu];
    if (iat >= jat) return;  // unique atom pairs, off-site only

    const int isha = ao2sh[mu], ishb = ao2sh[nu];
    const int la = ang[isha], lb = ang[ishb];
    const int sa = mu - iao_sh[isha], sb = nu - iao_sh[ishb];
    const bool dpair = (la >= 2 || lb >= 2);   // X-I1
    int ta = -1, tb = -1;
    if (!dpair) {
        ta = d_ao_to_type(la, sa);
        tb = d_ao_to_type(lb, sb);
        if (ta < 0 || tb < 0) return;
    }

    const double xa = xyz[3*iat+0], ya = xyz[3*iat+1], za = xyz[3*iat+2];
    const double xb = xyz[3*jat+0], yb = xyz[3*jat+1], zb = xyz[3*jat+2];
    const double dxij = xa - xb, dyij = ya - yb, dzij = za - zb;
    const double r2 = dxij*dxij + dyij*dyij + dzij*dzij;
    if (r2 < 1.0e-12) return;
    const double r = sqrt(r2);
    const int zi = z[iat], zj = z[jat];
    const double rad_sum = c_atomic_rad[zi-1] + c_atomic_rad[zj-1];
    const double rr = sqrt(r / rad_sum);
    const double pi_a = 1.0 + shpoly[isha] * rr, pi_b = 1.0 + shpoly[ishb] * rr;

    double hs;
    if (is_gfn2 == 0) {
        const bool vi = valence[isha] != 0, vj = valence[ishb] != 0;
        if (vi && vj) {
            double den = c_pauling[zi-1] - c_pauling[zj-1]; den *= den;
            hs = d_kpair_gfn1(zi, zj) * d_kshell_gfn1(la, lb) * (1.0 + (-7.0e-3) * den);
        } else if (vi && !vj) hs = 0.5 * (d_kshell_gfn1(la, la) + 2.85);
        else if (!vi && vj)   hs = 0.5 * (d_kshell_gfn1(lb, lb) + 2.85);
        else                  hs = 2.85;
    } else {
        double den = c_pauling[zi-1] - c_pauling[zj-1]; den *= den;
        const double enp = 1.0 + 2.0e-2 * den;
        const double km = d_kshell_gfn2(la, lb) * enp;
        const double za_ = sh_zeta[isha], zb_ = sh_zeta[ishb];
        hs = pow(2.0 * sqrt(za_ * zb_) / (za_ + zb_), 0.5) * km;
    }
    const double h_factor = hs * pi_a * pi_b;
    const double h_av = 0.5 * (se[isha] + se[ishb]) * h_factor;
    const double dlog_pi_dr_r = (shpoly[isha] / pi_a + shpoly[ishb] / pi_b) * rr / (2.0 * r2);

    double dS[3];
    if (dpair) {
        d_overlap_grad_elem(la, sa, lb, sb,
            prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
            prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
            xa, ya, za, xb, yb, zb, dS);
    } else {
        d_cgto_overlap_grad(
            prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
            prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
            xa, ya, za, xb, yb, zb, ta, tb, dS);
    }

    const double sval = 2.0*Pmn*h_av - 2.0*Wmn - Pmn*(v_ao[mu] + v_ao[nu]);
    const double shp  = 2.0*Pmn*H0mn*dlog_pi_dr_r;
    double Gx = sval*dS[0] + shp*dxij;
    double Gy = sval*dS[1] + shp*dyij;
    double Gz = sval*dS[2] + shp*dzij;

    // GFN2 multipole-integral Pulay term (xtb_gradient.cpp:350-414). The
    // transformed dp_int/qp_int derivatives are contracted with the converged
    // multipole potential v_dp/v_qp; G_sval[l] -= Pmn·term[l]. dp_int origin is
    // the column atom jat, so dR = R_jat − R_iat.
    if (is_gfn2 && v_dp) {
        double D_mp[3], dD_dA[3][3], dQ_dA[3][6];
        if (dpair) {
            d_multipole_grad_elem(la, sa, lb, sb,
                prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
                prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
                xa, ya, za, xb, yb, zb, D_mp, dD_dA, dQ_dA);
        } else {
            d_cgto_multipole_grad_transformed(
                prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
                prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
                xa, ya, za, xb, yb, zb, ta, tb, D_mp, dD_dA, dQ_dA);
        }
        const double dR[3] = { xb - xa, yb - ya, zb - za };
        const int qa6[6] = {0,0,1,0,1,2}, qb6[6] = {0,1,1,2,2,2};
        double term[3];
        for (int l = 0; l < 3; ++l) {
            double t = 0.0;
            for (int k = 0; k < 3; ++k) {
                t += dD_dA[l][k] * v_dp[k + jat*3];
                const double dDiat = dD_dA[l][k] + dR[k]*dS[l] - (k==l ? Smn : 0.0);
                t += dDiat * v_dp[k + iat*3];
            }
            double dqr[6];
            for (int q = 0; q < 6; ++q) {
                const int a = qa6[q], b = qb6[q];
                dqr[q] = -(b==l ? D_mp[a] : 0.0) + dR[b]*dD_dA[l][a]
                       -  (a==l ? D_mp[b] : 0.0) + dR[a]*dD_dA[l][b]
                       + (-(a==l ? dR[b] : 0.0) - (b==l ? dR[a] : 0.0))*Smn
                       + dR[a]*dR[b]*dS[l];
            }
            const double dtr_c = 0.5*(dqr[0] + dqr[2] + dqr[5]);
            for (int q = 0; q < 6; ++q) {
                const bool is_diag = (qa6[q] == qb6[q]);
                t += dQ_dA[l][q] * v_qp[q + jat*6];
                const double dQiat = dQ_dA[l][q] + 1.5*dqr[q] - (is_diag ? dtr_c : 0.0);
                t += dQiat * v_qp[q + iat*6];
            }
            term[l] = t;
        }
        Gx -= Pmn * term[0];
        Gy -= Pmn * term[1];
        Gz -= Pmn * term[2];
    }

    atomicAdd(&grad[3*iat+0], Gx); atomicAdd(&grad[3*jat+0], -Gx);
    atomicAdd(&grad[3*iat+1], Gy); atomicAdd(&grad[3*jat+1], -Gy);
    atomicAdd(&grad[3*iat+2], Gz); atomicAdd(&grad[3*jat+2], -Gz);

    const double cn_c = h_factor * Pmn * Smn;
    atomicAdd(&dEdcn[iat], (-kcn[isha]) * cn_c);
    atomicAdd(&dEdcn[jat], (-kcn[ishb]) * cn_c);
}

__global__ void k_grad_h0_pulay(
    int nao, int is_gfn2,
    const int* __restrict__ ao2sh, const int* __restrict__ ao2at, const int* __restrict__ ang,
    const int* __restrict__ iao_sh, const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly, const double* __restrict__ kcn,
    const int* __restrict__ valence, const int* __restrict__ z, const double* __restrict__ se,
    const double* __restrict__ xyz, const double* __restrict__ P, const double* __restrict__ S,
    const double* __restrict__ H0, const double* __restrict__ W, const double* __restrict__ v_ao,
    const double* __restrict__ v_dp, const double* __restrict__ v_qp,
    double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu >= nao || nu >= nao) return;
    if (ao2at[mu] >= ao2at[nu]) return;  // unique atom pairs, off-site only
    const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
    d_grad_h0_pulay_pair(mu, nu, S[mn], H0[mn], P[mn], W[mn], is_gfn2, ao2sh, ao2at, ang, iao_sh,
                         sh_nprim, sh_prim_off, prim_alpha, prim_coeff, sh_zeta, shpoly, kcn,
                         valence, z, se, xyz, v_ao, v_dp, v_qp, grad, dEdcn);
}

// Screened-pair twin: one thread per stored pair e = (row, col). P and W are dense (a full
// nao x nao buffer) but only ever read at the pattern position (mu,nu) - kept as the
// reference/fallback for callers that upload P/C explicitly (pc_resident=false in
// computeGradient) rather than reuse the resident screened-pattern density.
__global__ void k_grad_h0_pulay_sp(
    int nnz, const int* __restrict__ row, const int* __restrict__ col,
    const double* __restrict__ Ssp, const double* __restrict__ H0sp,
    int nao, int is_gfn2,
    const int* __restrict__ ao2sh, const int* __restrict__ ao2at, const int* __restrict__ ang,
    const int* __restrict__ iao_sh, const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly, const double* __restrict__ kcn,
    const int* __restrict__ valence, const int* __restrict__ z, const double* __restrict__ se,
    const double* __restrict__ xyz, const double* __restrict__ P, const double* __restrict__ W,
    const double* __restrict__ v_ao, const double* __restrict__ v_dp, const double* __restrict__ v_qp,
    double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    const int mu = row[e], nu = col[e];
    if (ao2at[mu] >= ao2at[nu]) return;
    const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
    d_grad_h0_pulay_pair(mu, nu, Ssp[e], H0sp[e], P[mn], W[mn], is_gfn2, ao2sh, ao2at, ang, iao_sh,
                         sh_nprim, sh_prim_off, prim_alpha, prim_coeff, sh_zeta, shpoly, kcn,
                         valence, z, se, xyz, v_ao, v_dp, v_qp, grad, dEdcn);
}

// Point 2 (Claude Generated, Sep 2026): pattern-only twin of k_grad_h0_pulay_sp - P and W
// come from the compact nnz-sized dSpP/dSpW arrays (same pair index e that already indexes
// Ssp/H0sp), no dense nao x nao buffer at all. Used on the device-resident path
// (pc_resident=true) where dSpP is already the converged density from the SCF loop's last
// iteration and dSpW is built once per gradient call (see computeGradient).
__global__ void k_grad_h0_pulay_spP(
    int nnz, const int* __restrict__ row, const int* __restrict__ col,
    const double* __restrict__ Ssp, const double* __restrict__ H0sp,
    const double* __restrict__ Psp, const double* __restrict__ Wsp,
    int is_gfn2,
    const int* __restrict__ ao2sh, const int* __restrict__ ao2at, const int* __restrict__ ang,
    const int* __restrict__ iao_sh, const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly, const double* __restrict__ kcn,
    const int* __restrict__ valence, const int* __restrict__ z, const double* __restrict__ se,
    const double* __restrict__ xyz,
    const double* __restrict__ v_ao, const double* __restrict__ v_dp, const double* __restrict__ v_qp,
    double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    const int mu = row[e], nu = col[e];
    if (ao2at[mu] >= ao2at[nu]) return;
    d_grad_h0_pulay_pair(mu, nu, Ssp[e], H0sp[e], Psp[e], Wsp[e], is_gfn2, ao2sh, ao2at, ang, iao_sh,
                         sh_nprim, sh_prim_off, prim_alpha, prim_coeff, sh_zeta, shpoly, kcn,
                         valence, z, se, xyz, v_ao, v_dp, v_qp, grad, dEdcn);
}

// Section 3: isotropic Coulomb gradient. One thread per shell is, inner js<is.
//
// Point 5a (Claude Generated, Sep 2026): the only call site (computeGradient, below) always
// passes gexp=2.0 (xtb_coulomb.hpp:98, Klopman-Ohno kernel, both GFN1/GFN2 - the energy-side
// gamma-matrix kernel a few hundred lines up already hardcodes the same 2.0 for the same
// reason). The general path costs 4 FP64 pow() per shell pair; on this card (Ampere-class)
// FP64 transcendentals run at a fraction of the FP64 add/mul rate, so this loop was
// essentially the whole measured 234 ms Coulomb-gradient phase on polymer_2x (7320 atoms).
//
// For gexp==2 exactly, let X = r^2 + gam_bar^-2, so gamma = X^(-1/2). The code's "dgamma_dr"
// is actually d(gamma)/dr * (1/r) (that 1/r is what lets grad += force*dx work directly on
// the raw coordinate difference, no separate normalisation below) - algebraically
// dgamma/dr * (1/r) = -r^(gexp-2)*gamma^(gexp+1), and r^(gexp-2) = r^0 = 1 identically for
// gexp=2, leaving -gamma^3 with no r-dependence left to evaluate at all. Verified against the
// general formula symbolically (both reduce to the same closed form) and numerically at
// several (r, gam_bar) points before wiring in. NOT bit-identical to the general path -
// gamma*gamma*gamma rounds one ULP differently than pow(gamma,3.0) at some inputs - see the
// noise-floor/validation numbers in the report; treat this the same way the shell-pair-blocked
// integral kernels' FMA-contraction note is treated (docs/SQM_PERFORMANCE.md).
__global__ void k_grad_coulomb(int nsh, int is_gfn2, const int* __restrict__ sh2at,
                               const double* __restrict__ g, const double* __restrict__ q_sh,
                               const double* __restrict__ xyz, double gexp, double* __restrict__ grad)
{
    const int is = blockIdx.x * blockDim.x + threadIdx.x;
    if (is >= nsh) return;
    const int iat = sh2at[is];
    const bool fast_gexp2 = (gexp == 2.0);
    for (int js = 0; js < is; ++js) {
        const int jat = sh2at[js];
        if (iat == jat) continue;
        const double dx = xyz[3*iat+0] - xyz[3*jat+0];
        const double dy = xyz[3*iat+1] - xyz[3*jat+1];
        const double dz = xyz[3*iat+2] - xyz[3*jat+2];
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1.0e-12) continue;
        const double gam_bar = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
        double dgamma_dr;
        if (fast_gexp2) {
            const double gamma = 1.0 / sqrt(r2 + 1.0 / (gam_bar * gam_bar));
            dgamma_dr = -(gamma * gamma * gamma);
        } else {
            const double r1 = sqrt(r2);
            const double gamma = pow(pow(r1, gexp) + pow(gam_bar, -gexp), -1.0 / gexp);
            dgamma_dr = -pow(r1, gexp - 2.0) * pow(gamma, gexp + 1.0);
        }
        const double force = q_sh[is] * q_sh[js] * dgamma_dr;
        atomicAdd(&grad[3*iat+0], force*dx); atomicAdd(&grad[3*jat+0], -force*dx);
        atomicAdd(&grad[3*iat+1], force*dy); atomicAdd(&grad[3*jat+1], -force*dy);
        atomicAdd(&grad[3*iat+2], force*dz); atomicAdd(&grad[3*jat+2], -force*dz);
    }
}

// ---- Device kernels (Stage 2 resident SCF) -------------------------------
// All matrices column-major n×n: element (i,j) at index i + j*n.

// F = H0 − ½·S·(v_ao(i) + v_ao(j))   (GFN1 isotropic Fock; H0, S symmetric).
__global__ void k_build_fock_iso(double* F, const double* H0, const double* S,
                                 const double* vao, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < n && j < n) {
        const size_t idx = static_cast<size_t>(i) + static_cast<size_t>(j) * n;
        F[idx] = H0[idx] - 0.5 * S[idx] * (vao[i] + vao[j]);
    }
}

// FP64↔FP32 element-wise conversion (mixed-precision eigensolve).
__global__ void k_d2f(const double* __restrict__ in, float* __restrict__ out, size_t n)
{
    const size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = static_cast<float>(in[i]);
}
// Claude Generated (Sep 2026): deterministic pseudo-random probe vector and elementwise scale,
// used by verifyDistributedFp32() to check that a distributed eigensolve really returned
// eigenpairs of the matrix it was given.
__global__ void k_probe_fill(float* __restrict__ v, int n, unsigned seed)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    unsigned h = seed ^ (static_cast<unsigned>(i) * 2654435761u);
    h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
    v[i] = (static_cast<float>(h & 0xFFFFFu) / 524288.0f) - 1.0f;   // in [-1, 1)
}

__global__ void k_scale_by(float* __restrict__ x, const float* __restrict__ s, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) x[i] *= s[i];
}

__global__ void k_probe_filld(double* __restrict__ v, int n, unsigned seed)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    unsigned h = seed ^ (static_cast<unsigned>(i) * 2654435761u);
    h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
    v[i] = (static_cast<double>(h & 0xFFFFFu) / 524288.0) - 1.0;
}

__global__ void k_scale_byd(double* __restrict__ x, const double* __restrict__ s, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) x[i] *= s[i];
}

__global__ void k_f2d(const float* __restrict__ in, double* __restrict__ out, size_t n)
{
    const size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = static_cast<double>(in[i]);
}

// Scale the leading ncol columns of C by the occupation: Cw(:,k) = occ[k]·C(:,k).
__global__ void k_scale_cols(double* Cw, const double* C, const double* occ,
                             int n, int ncol)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int k = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < n && k < ncol) {
        const size_t idx = static_cast<size_t>(i) + static_cast<size_t>(k) * n;
        Cw[idx] = C[idx] * occ[k];
    }
}

// Mulliken AO populations pop(μ) = Σ_ν P(μ,ν)·S(μ,ν)  (S symmetric → = Σ P_μν S_νμ).
__global__ void k_pop_ao(double* pop, const double* P, const double* S, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu < n) {
        double s = 0.0;
        for (int nu = 0; nu < n; ++nu) {
            const size_t idx = static_cast<size_t>(mu) + static_cast<size_t>(nu) * n;
            s += P[idx] * S[idx];
        }
        pop[mu] = s;
    }
}

// GFN2 anisotropic Fock contribution (mirrors the multipole branch of XTB::buildFock):
//   F(μ,ν) −= ½·[ Σ_k dp_int[k](μ,ν)·v_dp(k,jat) + dp_int[k](ν,μ)·v_dp(k,iat)
//              + Σ_k qp_int[k](μ,ν)·v_qp(k,jat) + qp_int[k](ν,μ)·v_qp(k,iat) ]
// with iat=ao2at[μ], jat=ao2at[ν]. dp_int contiguous 3·nn, qp_int 6·nn (col-major).
__global__ void k_add_fock_multipole(double* F, const double* dp_int, const double* qp_int,
                                     const double* v_dp, const double* v_qp,
                                     const int* ao2at, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu < n && nu < n) {
        const size_t nn  = static_cast<size_t>(n) * static_cast<size_t>(n);
        const size_t mn  = static_cast<size_t>(mu) + static_cast<size_t>(nu) * n; // (μ,ν)
        const size_t nm  = static_cast<size_t>(nu) + static_cast<size_t>(mu) * n; // (ν,μ)
        const int iat = ao2at[mu];
        const int jat = ao2at[nu];
        double dd = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double* dk = dp_int + static_cast<size_t>(k) * nn;
            dd += dk[mn] * v_dp[k + jat * 3] + dk[nm] * v_dp[k + iat * 3];
        }
        double qq = 0.0;
        for (int k = 0; k < 6; ++k) {
            const double* qk = qp_int + static_cast<size_t>(k) * nn;
            qq += qk[mn] * v_qp[k + jat * 6] + qk[nm] * v_qp[k + iat * 6];
        }
        F[mn] -= 0.5 * (dd + qq);
    }
}

// GFN2 atomic multipole moments (mirrors the multipole block of updatePopulations):
//   dp_at(k,iat) −= Σ_ν P(ν,μ)·dp_int[k](ν,μ)   summed over μ∈iat   (3×nat)
//   qp_at(k,iat) −= Σ_ν P(ν,μ)·qp_int[k](ν,μ)   summed over μ∈iat   (6×nat)
// One thread per μ; column-μ dot products, atomic-scattered into the owning atom.
// dp_at/qp_at must be zeroed before launch.
__global__ void k_multipole_moments(double* dp_at, double* qp_at,
                                    const double* P, const double* dp_int, const double* qp_int,
                                    const int* ao2at, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu < n) {
        const size_t nn   = static_cast<size_t>(n) * static_cast<size_t>(n);
        const size_t col  = static_cast<size_t>(mu) * n;   // start of column μ
        const int    iat  = ao2at[mu];
        for (int k = 0; k < 3; ++k) {
            const double* dk = dp_int + static_cast<size_t>(k) * nn + col;
            double acc = 0.0;
            for (int nu = 0; nu < n; ++nu) acc += P[col + nu] * dk[nu];
            atomicAdd(&dp_at[k + iat * 3], -acc);
        }
        for (int k = 0; k < 6; ++k) {
            const double* qk = qp_int + static_cast<size_t>(k) * nn + col;
            double acc = 0.0;
            for (int nu = 0; nu < n; ++nu) acc += P[col + nu] * qk[nu];
            atomicAdd(&qp_at[k + iat * 6], -acc);
        }
    }
}

// ====================================================================== *
//  Screened (sparse) AO-pair storage for S, H0 and the GFN2 multipole integrals.
//  Claude Generated (Sep 2026, large systems).
//
//  The Gaussian integrals decay as exp(-a_i a_j/(a_i+a_j) R^2) with the distance R of
//  the two atoms, so beyond a basis-derived cutoff every element of S, H0, dp_int and
//  qp_int is below 1e-20 and contributes nothing representable to the SCF. Storing only
//  the AO pairs of atom pairs within that cutoff turns eleven nao^2 matrices into eleven
//  nnz-length arrays (polymer_2x: 16 % of the pairs within 40 Bohr).
//
//  Layout: entry e is the AO pair (row[e], col[e]), ordered column-major (col ascending,
//  then row ascending) - the same order as the dense buffers - and colptr[mu] is the first
//  entry of column mu. perm[e] is the entry of the transposed pair (col[e], row[e]); the
//  atom-distance screen is symmetric, so it always exists. Every kernel below performs the
//  same arithmetic, in the same order, as its dense counterpart, so with nothing screened
//  away (forced mode on a small molecule) the results are identical to the dense path.
// ====================================================================== *

// values[e] = dense[row[e] + col[e]*n]
__global__ void k_sp_gather(const double* __restrict__ dense, const int* __restrict__ row,
                            const int* __restrict__ col, int n, int nnz, double* __restrict__ values)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    values[e] = dense[static_cast<size_t>(row[e]) + static_cast<size_t>(col[e]) * n];
}

// dense[row[e] + col[e]*n] = values[e]   (dense must be zeroed before)
__global__ void k_sp_scatter(const double* __restrict__ values, const int* __restrict__ row,
                             const int* __restrict__ col, int n, int nnz, double* __restrict__ dense)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    dense[static_cast<size_t>(row[e]) + static_cast<size_t>(col[e]) * n] = values[e];
}

// Sparse twin of k_multipole_ints: one thread per stored pair, S taken from the screened
// overlap values of the same entry. dp/qp are 3*nnz / 6*nnz, component k at k*nnz + e.
__global__ void k_multipole_ints_sp(
    int nnz, const int* __restrict__ row, const int* __restrict__ col,
    const int* __restrict__ ao2sh, const int* __restrict__ ao2at,
    const int* __restrict__ iao_sh, const int* __restrict__ ang_sh,
    const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ xyz, const double* __restrict__ Ssp,
    double* __restrict__ dp_sp, double* __restrict__ qp_sp)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    const size_t ne = static_cast<size_t>(nnz);
    const size_t ee = static_cast<size_t>(e);
    const int mu = row[e], nu = col[e];

    const int isha = ao2sh[mu], iat = ao2at[mu];
    const int ishb = ao2sh[nu], jat = ao2at[nu];
    const int la = ang_sh[isha], lb = ang_sh[ishb];
    const int sa = mu - iao_sh[isha], sb = nu - iao_sh[ishb];
    const bool dpair = (la >= 2 || lb >= 2);

    for (int k = 0; k < 3; ++k) dp_sp[static_cast<size_t>(k) * ne + ee] = 0.0;
    for (int k = 0; k < 6; ++k) qp_sp[static_cast<size_t>(k) * ne + ee] = 0.0;

    const double* aA = prim_alpha + sh_prim_off[isha];
    const double* cA = prim_coeff + sh_prim_off[isha];
    const int npa = sh_nprim[isha];
    const double* aB = prim_alpha + sh_prim_off[ishb];
    const double* cB = prim_coeff + sh_prim_off[ishb];
    const int npb = sh_nprim[ishb];

    double Sx, D[3], Q[6];
    if (!dpair) {
        const int ta = d_ao_to_type(la, sa);
        const int tb = d_ao_to_type(lb, sb);
        if (ta < 0 || tb < 0) return;
        d_cgto_multipole(aA, cA, npa, aB, cB, npb,
                         xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                         xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], ta, tb, Sx, D, Q);
    } else {
        d_multipole_elem(la, sa, lb, sb, aA, cA, npa, aB, cB, npb,
                         xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                         xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], Sx, D, Q);
    }

    const double Rx = xyz[3*jat+0], Ry = xyz[3*jat+1], Rz = xyz[3*jat+2];
    const double Smn = Ssp[e];
    const double dx = D[0], dy = D[1], dz = D[2];
    dp_sp[0*ne + ee] = dx - Rx * Smn;
    dp_sp[1*ne + ee] = dy - Ry * Smn;
    dp_sp[2*ne + ee] = dz - Rz * Smn;

    const double qxx = Q[0] - 2*Rx*dx + Rx*Rx*Smn;
    const double qxy = Q[1] - Rx*dy - Ry*dx + Rx*Ry*Smn;
    const double qyy = Q[2] - 2*Ry*dy + Ry*Ry*Smn;
    const double qxz = Q[3] - Rx*dz - Rz*dx + Rx*Rz*Smn;
    const double qyz = Q[4] - Ry*dz - Rz*dy + Ry*Rz*Smn;
    const double qzz = Q[5] - 2*Rz*dz + Rz*Rz*Smn;
    const double tr = 0.5 * (qxx + qyy + qzz);
    qp_sp[0*ne + ee] = 1.5 * qxx - tr;
    qp_sp[1*ne + ee] = 1.5 * qxy;
    qp_sp[2*ne + ee] = 1.5 * qyy - tr;
    qp_sp[3*ne + ee] = 1.5 * qxz;
    qp_sp[4*ne + ee] = 1.5 * qyz;
    qp_sp[5*ne + ee] = 1.5 * qzz - tr;
}

// Sparse twin of k_build_fock_iso (+ optional k_add_fock_multipole): F must be zeroed first.
// F(mu,nu) = H0 - 1/2 S (v_mu + v_nu) - 1/2 [dp/qp(mu,nu).v(jat) + dp/qp(nu,mu).v(iat)].
__global__ void k_build_fock_sp(double* __restrict__ F, int n, int nnz,
                                const int* __restrict__ row, const int* __restrict__ col,
                                const int* __restrict__ perm,
                                const double* __restrict__ H0sp, const double* __restrict__ Ssp,
                                const double* __restrict__ vao,
                                const double* __restrict__ dp_sp, const double* __restrict__ qp_sp,
                                const double* __restrict__ v_dp, const double* __restrict__ v_qp,
                                const int* __restrict__ ao2at)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    const int mu = row[e], nu = col[e];
    const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * n;
    double f = H0sp[e] - 0.5 * Ssp[e] * (vao[mu] + vao[nu]);
    if (dp_sp) {
        const size_t ne = static_cast<size_t>(nnz);
        const size_t em = static_cast<size_t>(e), et = static_cast<size_t>(perm[e]);
        const int iat = ao2at[mu], jat = ao2at[nu];
        double dd = 0.0;
        for (int k = 0; k < 3; ++k) {
            const size_t off = static_cast<size_t>(k) * ne;
            dd += dp_sp[off + em] * v_dp[k + jat * 3] + dp_sp[off + et] * v_dp[k + iat * 3];
        }
        double qq = 0.0;
        for (int k = 0; k < 6; ++k) {
            const size_t off = static_cast<size_t>(k) * ne;
            qq += qp_sp[off + em] * v_qp[k + jat * 6] + qp_sp[off + et] * v_qp[k + iat * 6];
        }
        f -= 0.5 * (dd + qq);
    }
    F[mn] = f;
}

// Sparse twin of k_pop_ao: pop(mu) = sum_nu P(mu,nu) S(mu,nu), nu ascending. Column mu of
// the pattern lists exactly the nu with S(nu,mu) stored; S(mu,nu) is at perm[e].
__global__ void k_pop_ao_sp(double* __restrict__ pop, const double* __restrict__ P,
                            const double* __restrict__ Ssp, const int* __restrict__ row,
                            const int* __restrict__ colptr, const int* __restrict__ perm, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= n) return;
    double s = 0.0;
    for (int e = colptr[mu]; e < colptr[mu + 1]; ++e) {
        const size_t idx = static_cast<size_t>(mu) + static_cast<size_t>(row[e]) * n;
        s += P[idx] * Ssp[perm[e]];
    }
    pop[mu] = s;
}

// Sparse twin of k_multipole_moments (column mu, rows ascending).
__global__ void k_multipole_moments_sp(double* dp_at, double* qp_at, const double* __restrict__ P,
                                       const double* __restrict__ dp_sp, const double* __restrict__ qp_sp,
                                       const int* __restrict__ row, const int* __restrict__ colptr,
                                       const int* __restrict__ ao2at, int n, int nnz)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= n) return;
    const size_t ne  = static_cast<size_t>(nnz);
    const size_t col = static_cast<size_t>(mu) * n;
    const int iat = ao2at[mu];
    for (int k = 0; k < 3; ++k) {
        const double* dk = dp_sp + static_cast<size_t>(k) * ne;
        double acc = 0.0;
        for (int e = colptr[mu]; e < colptr[mu + 1]; ++e) acc += P[col + row[e]] * dk[e];
        atomicAdd(&dp_at[k + iat * 3], -acc);
    }
    for (int k = 0; k < 6; ++k) {
        const double* qk = qp_sp + static_cast<size_t>(k) * ne;
        double acc = 0.0;
        for (int e = colptr[mu]; e < colptr[mu + 1]; ++e) acc += P[col + row[e]] * qk[e];
        atomicAdd(&qp_at[k + iat * 6], -acc);
    }
}

// Density on the screened pattern only (Claude Generated, Sep 2026). Populations, multipole
// moments and the band energy read P exclusively at stored pairs, so the SCF loop never needs
// the dense nao x nao product: P(row,col) = sum_k Cw(row,k) C(col,k) with Cw = C diag(occ),
// the same sum the dense GEMM P = Cw C^T performs for that element. polymer_2x: 6 % of the
// elements (the dense GEMM was 26 % of the SCF time) and no dense P buffer during the loop.
__global__ void k_density_sp(int nnz, const int* __restrict__ row, const int* __restrict__ col,
                             const double* __restrict__ Cw, const double* __restrict__ C,
                             int n, int ncol, double* __restrict__ Psp)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nnz) return;
    const size_t r = static_cast<size_t>(row[e]), c = static_cast<size_t>(col[e]);
    double acc = 0.0;
    for (int k = 0; k < ncol; ++k) {
        const size_t off = static_cast<size_t>(k) * n;
        acc += Cw[r + off] * C[c + off];
    }
    Psp[e] = acc;
}

// pop(mu) = sum_nu P(mu,nu) S(mu,nu) from pattern-only P (both at the transposed entry).
__global__ void k_pop_ao_spP(double* __restrict__ pop, const double* __restrict__ Psp,
                             const double* __restrict__ Ssp, const int* __restrict__ colptr,
                             const int* __restrict__ perm, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= n) return;
    double s = 0.0;
    for (int e = colptr[mu]; e < colptr[mu + 1]; ++e) {
        const int t = perm[e];
        s += Psp[t] * Ssp[t];
    }
    pop[mu] = s;
}

// Multipole moments from pattern-only P: column mu, P(row,mu) = Psp[e].
__global__ void k_multipole_moments_spP(double* dp_at, double* qp_at, const double* __restrict__ Psp,
                                        const double* __restrict__ dp_sp, const double* __restrict__ qp_sp,
                                        const int* __restrict__ colptr, const int* __restrict__ ao2at,
                                        int n, int nnz)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= n) return;
    const size_t ne = static_cast<size_t>(nnz);
    const int iat = ao2at[mu];
    for (int k = 0; k < 3; ++k) {
        const double* dk = dp_sp + static_cast<size_t>(k) * ne;
        double acc = 0.0;
        for (int e = colptr[mu]; e < colptr[mu + 1]; ++e) acc += Psp[e] * dk[e];
        atomicAdd(&dp_at[k + iat * 3], -acc);
    }
    for (int k = 0; k < 6; ++k) {
        const double* qk = qp_sp + static_cast<size_t>(k) * ne;
        double acc = 0.0;
        for (int e = colptr[mu]; e < colptr[mu + 1]; ++e) acc += Psp[e] * qk[e];
        atomicAdd(&qp_at[k + iat * 6], -acc);
    }
}

// Post-SCF: the whole 2-body D4 (energy + nuclear gradient + dE/dCN + dE/dq) in one per-atom
// GATHER — a superset of k_d4_dedq. No atomics: the pair force is antisymmetric (atom j's own
// thread accumulates the mirror term) and the energy is per-atom (host halves Σ e_atom). Port of
// D4Evaluator::computeEnergyAndGradient's per-reference 2-body path: the direct radial gradient is
// dE_dr = -C6·ddisp_dr2·(R_i-R_j) with ddisp_dr2 = s6·(-6r⁴t6²) + s8·r4r2·(-8r⁶t8²) (the C6·ζc6
// rescale cancels in the per-reference path), dEdcn_i = -dc6dcn_i·disp_sum, dEdq_i = -dc6dq_i·disp_sum.
// Claude Generated. CUDA port of the ROCm kernel (Sep 2026), unchanged.
__global__ void k_d4_grad(int nat, int max_elem, int max_ref,
                          const int* __restrict__ Z, const int* __restrict__ nref,
                          const double* __restrict__ sqrtZr4r2, const double* __restrict__ xyz,
                          const double* __restrict__ c6_flat,
                          const double* __restrict__ W, const double* __restrict__ dWq,
                          const double* __restrict__ dWc,
                          double s6, double s8, double a1, double a2, double cut2,
                          double* __restrict__ e_atom, double* __restrict__ grad,
                          double* __restrict__ dEdcn, double* __restrict__ dEdq)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int ei = Z[i] - 1;
    const int nri = nref[i];
    if (ei < 0 || ei >= max_elem || nri <= 0) {
        e_atom[i] = 0.0; grad[3*i+0] = grad[3*i+1] = grad[3*i+2] = 0.0; dEdcn[i] = 0.0; dEdq[i] = 0.0;
        return;
    }
    const double xi = xyz[3 * i + 0], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
    const double sq_i = sqrtZr4r2[i];
    const double* Wi   = W   + static_cast<size_t>(i) * max_ref;
    const double* dWqi = dWq + static_cast<size_t>(i) * max_ref;
    const double* dWci = dWc + static_cast<size_t>(i) * max_ref;

    double eacc = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, cnacc = 0.0, qacc = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const double dx = xi - xyz[3 * j + 0];
        const double dy = yi - xyz[3 * j + 1];
        const double dz = zi - xyz[3 * j + 2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cut2 || r2 < 1.0e-20) continue;
        const int ej = Z[j] - 1;
        const int nrj = nref[j];
        if (ej < 0 || ej >= max_elem || nrj <= 0) continue;

        const double r4r2ij = 3.0 * sq_i * sqrtZr4r2[j];
        const double r0 = a1 * sqrt(r4r2ij) + a2;
        const double r0_2 = r0 * r0;
        const double r0_6 = r0_2 * r0_2 * r0_2;
        const double r0_8 = r0_6 * r0_2;
        const double r6 = r2 * r2 * r2;
        const double r8 = r6 * r2;
        const double t6 = 1.0 / (r6 + r0_6);
        const double t8 = 1.0 / (r8 + r0_8);
        const double disp_sum = s6 * t6 + s8 * r4r2ij * t8;
        const double d6 = -6.0 * r2 * r2 * t6 * t6;            // -6·r⁴·t6²
        const double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;       // -8·r⁶·t8²
        const double ddisp_dr2 = s6 * d6 + s8 * r4r2ij * d8;

        const double* Wj = W + static_cast<size_t>(j) * max_ref;
        const size_t base = (static_cast<size_t>(ei) * max_elem + ej)
                          * static_cast<size_t>(max_ref) * max_ref;
        double C6 = 0.0, dc6dcni = 0.0, dc6dqi = 0.0;
        for (int a = 0; a < nri; ++a) {
            const size_t basea = base + static_cast<size_t>(a) * max_ref;
            double sW = 0.0;
            for (int b = 0; b < nrj; ++b) sW += Wj[b] * c6_flat[basea + b];
            C6      += Wi[a]   * sW;
            dc6dcni += dWci[a] * sW;
            dc6dqi  += dWqi[a] * sW;
        }
        eacc  += -C6 * disp_sum;                 // per-atom energy (host: E = ½·Σ e_atom)
        const double f = -C6 * ddisp_dr2;        // dE_dr = f·(R_i-R_j)
        gx += f * dx; gy += f * dy; gz += f * dz;
        cnacc += -dc6dcni * disp_sum;            // CN chain (host distributes via ∂CN/∂R)
        qacc  += -dc6dqi  * disp_sum;            // q-response first half (host folds ∂q/∂R)
    }
    e_atom[i] = eacc;
    grad[3 * i + 0] = gx; grad[3 * i + 1] = gy; grad[3 * i + 2] = gz;
    dEdcn[i] = cnacc;
    dEdq[i]  = qacc;
}

// D4 ATM 3-body (energy + nuclear gradient + dE/dCN), CUDA port of the ROCm k_d4_atm (Claude
// Generated, Sep 2026) with one change: the pair loops run over a host-built neighbour list of
// atoms within the ATM cutoff (nb[nbptr[a]..nbptr[a+1]), ascending) instead of all atoms. Every
// skipped pair failed the r2 > cut2 test in the ROCm kernel and contributed nothing, and the
// remaining pairs are visited in the same (ascending) order, so the sums are the same. At 7320
// atoms this turns ~2e11 distance tests into ~5e8.
__global__ void k_d4_atm_nl(int nat, const double* __restrict__ xyz, const double* __restrict__ r4r2,
                            const double* __restrict__ c6, const double* __restrict__ dc6dcn,
                            const int* __restrict__ nbptr, const int* __restrict__ nb,
                            double s9, double a1, double a2, double alp, double cut2,
                            double* __restrict__ e_atom, double* __restrict__ grad,
                            double* __restrict__ dEdcn)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    const double eps = 2.220446049250313e-16;
    const double xa = xyz[3*a+0], ya = xyz[3*a+1], za = xyz[3*a+2];
    const double r4r2a = r4r2[a];
    const size_t nat_s = static_cast<size_t>(nat);
    const bool alp16 = (alp == 16.0);
    double eacc = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, cnacc = 0.0;
    for (int ix = nbptr[a]; ix < nbptr[a + 1]; ++ix) {
        const int x = nb[ix];
        const double vaxx = xyz[3*x+0]-xa, vaxy = xyz[3*x+1]-ya, vaxz = xyz[3*x+2]-za;
        const double r2ax = vaxx*vaxx + vaxy*vaxy + vaxz*vaxz;
        if (r2ax > cut2 || r2ax < eps) continue;
        const double c6ax = c6[a*nat_s + x];
        const double r0ax = a1*sqrt(3.0*r4r2a*r4r2[x]) + a2;
        for (int iy = nbptr[a]; iy < ix; ++iy) {        // y < x (list ascending)
            const int y = nb[iy];
            const double vayx = xyz[3*y+0]-xa, vayy = xyz[3*y+1]-ya, vayz = xyz[3*y+2]-za;
            const double r2ay = vayx*vayx + vayy*vayy + vayz*vayz;
            if (r2ay > cut2 || r2ay < eps) continue;
            const double dxyx = xyz[3*y+0]-xyz[3*x+0], dxyy = xyz[3*y+1]-xyz[3*x+1], dxyz_ = xyz[3*y+2]-xyz[3*x+2];
            const double r2xy = dxyx*dxyx + dxyy*dxyy + dxyz_*dxyz_;
            if (r2xy > cut2 || r2xy < eps) continue;

            const double c6ay = c6[a*nat_s + y];
            const double c6xy = c6[x*nat_s + y];
            const double c9 = -s9 * sqrt(fabs(c6ax*c6ay*c6xy));
            const double r0ay = a1*sqrt(3.0*r4r2a*r4r2[y]) + a2;
            const double r0xy = a1*sqrt(3.0*r4r2[x]*r4r2[y]) + a2;
            const double r0 = r0ax*r0ay*r0xy;
            const double r2 = r2ax*r2ay*r2xy;
            const double r1 = sqrt(r2);
            const double r3 = r2*r1;
            const double r5 = r3*r2;
            // alp = 16 (GFN2): (r0/r1)^(16/3) = q^5 * cbrt(q), far cheaper than pow on the device.
            const double q = r0/r1;
            const double pw = alp16 ? (q*q*q*q*q) * cbrt(q) : pow(q, alp/3.0);
            const double fdmp = 1.0/(1.0 + 6.0*pw);
            const double ang = 0.375*(r2ax + r2xy - r2ay)*(r2ax - r2xy + r2ay)*(-r2ax + r2xy + r2ay)/r5 + 1.0/r3;
            const double dE = ang*fdmp*c9;
            eacc += -dE;

            const double dfdmp = -2.0*alp*pw*fdmp*fdmp;
            const double dang_ax = -0.375*(r2ax*r2ax*r2ax + r2ax*r2ax*(r2xy+r2ay)
                          + r2ax*(3.0*r2xy*r2xy + 2.0*r2xy*r2ay + 3.0*r2ay*r2ay)
                          - 5.0*(r2xy-r2ay)*(r2xy-r2ay)*(r2xy+r2ay))/r5;
            const double gcax = c9*(-dang_ax*fdmp + ang*dfdmp)/r2ax;
            gx += -gcax*vaxx; gy += -gcax*vaxy; gz += -gcax*vaxz;
            const double dang_ay = -0.375*(r2ay*r2ay*r2ay + r2ay*r2ay*(r2xy+r2ax)
                          + r2ay*(3.0*r2xy*r2xy + 2.0*r2xy*r2ax + 3.0*r2ax*r2ax)
                          - 5.0*(r2xy-r2ax)*(r2xy-r2ax)*(r2xy+r2ax))/r5;
            const double gcay = c9*(-dang_ay*fdmp + ang*dfdmp)/r2ay;
            gx += -gcay*vayx; gy += -gcay*vayy; gz += -gcay*vayz;
            cnacc += -dE*0.5*(dc6dcn[a*nat_s+x]/c6ax + dc6dcn[a*nat_s+y]/c6ay);
        }
    }
    e_atom[a] = eacc;
    grad[3*a+0] = gx; grad[3*a+1] = gy; grad[3*a+2] = gz;
    dEdcn[a] = cnacc;
}

// ====================================================================== *
//  Stage 5 (Part A): single-shot D4 EEQ charge model device kernels.
//  Verbatim port of curcuma::dispersion::D4ChargeModel (d4_charge_model.cpp):
//  one smooth augmented linear system [[A,1],[1,0]]·[q;λ]=[b;Q], solved via LU,
//  plus the analytic ∂q/∂x charge-response (adjoint + closed-form pair loop).
//  Constants mirror the CPU file exactly. Claude Generated (2026-06).
// ====================================================================== *
namespace {
constexpr double D4EEQ_TSQRT2PI        = 0.797884560802866;   // sqrt(2/π)
constexpr double D4EEQ_TWO_OVER_SQRTPI = 1.1283791670955126;  // 2/sqrt(π)
constexpr double D4EEQ_KN     = -7.5;       // GFN-FF erf-CN steepness
constexpr double D4EEQ_CNMAX  = 4.4;        // CN log-compression cap
constexpr double D4EEQ_CN_EPS = 1.0e-10;    // guard for 1/sqrt(CN)
} // namespace

// CN (GFN-FF log-compressed erf form) + raw CN. One thread per atom i. rcov_bohr
// is the pre-scaled (4/3·rcov·Å→Bohr) covalent radius; rcov==0 → CN 0 (skipped).
// Mirrors D4ChargeModel::computeCharges CN loop (d4_charge_model.cpp:58-73).
__global__ void k_d4eeq_cn(int N, const double* __restrict__ xyz,
                           const double* __restrict__ rcov_bohr,
                           double* __restrict__ cn, double* __restrict__ cn_raw)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    const double rci = rcov_bohr[i];
    if (rci == 0.0) { cn[i] = 0.0; cn_raw[i] = 0.0; return; }
    const double xi = xyz[3 * i + 0], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
    double raw = 0.0;
    for (int j = 0; j < N; ++j) {
        if (j == i) continue;
        const double rcj = rcov_bohr[j];
        if (rcj == 0.0) continue;
        const double dx = xi - xyz[3 * j + 0];
        const double dy = yi - xyz[3 * j + 1];
        const double dz = zi - xyz[3 * j + 2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);
        const double rcij = rci + rcj;
        const double dr = (r - rcij) / rcij;
        raw += 0.5 * (1.0 + erf(D4EEQ_KN * dr));
    }
    cn_raw[i] = raw;
    const double log1p_ecnmax = log(1.0 + exp(D4EEQ_CNMAX));
    cn[i] = log1p_ecnmax - log(1.0 + exp(D4EEQ_CNMAX - raw));
}

// Augmented EEQ matrix M = [[A,1],[1ᵀ,0]], column-major (N+1)×(N+1). One thread
// per (row r, col c). A_ii = γ_i + √(2/π)/√α_i² ; A_ij = erf(γ_ij·r)/r with
// γ_ij = 1/√(α_i²+α_j²); border = 1; corner = 0. (d4_charge_model.cpp:79-92).
__global__ void k_d4eeq_build(int N, const double* __restrict__ xyz,
                              const double* __restrict__ alpha_sq,
                              const double* __restrict__ gam,
                              double* __restrict__ M)
{
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    const int c = blockIdx.y * blockDim.y + threadIdx.y;
    const int m = N + 1;
    if (r >= m || c >= m) return;
    double val;
    if (r == N && c == N) {
        val = 0.0;
    } else if (r == N || c == N) {
        val = 1.0;
    } else if (r == c) {
        val = gam[r] + D4EEQ_TSQRT2PI / sqrt(alpha_sq[r]);
    } else {
        const double dx = xyz[3 * r + 0] - xyz[3 * c + 0];
        const double dy = xyz[3 * r + 1] - xyz[3 * c + 1];
        const double dz = xyz[3 * r + 2] - xyz[3 * c + 2];
        const double rr = sqrt(dx * dx + dy * dy + dz * dz);
        const double gammij = 1.0 / sqrt(alpha_sq[r] + alpha_sq[c]);
        val = erf(gammij * rr) / rr;
    }
    M[static_cast<size_t>(r) + static_cast<size_t>(c) * m] = val;
}

// RHS c = [b; Q], length N+1. b_i = -χ_i + κ_i·√max(CN_i,0); c(N)=total_charge.
__global__ void k_d4eeq_rhs(int N, const double* __restrict__ chi,
                            const double* __restrict__ cnf, const double* __restrict__ cn,
                            double total_charge, double* __restrict__ c_out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > N) return;
    if (i == N) { c_out[N] = total_charge; return; }
    c_out[i] = -chi[i] + cnf[i] * sqrt(fmax(cn[i], 0.0));
}

// Adjoint RHS [dEdq; 0], length N+1 (for M·z = [dEdq;0]).
__global__ void k_d4eeq_adjoint_rhs(int N, const double* __restrict__ dEdq,
                                    double* __restrict__ rhs)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > N) return;
    rhs[i] = (i < N) ? dEdq[i] : 0.0;
}

// Per-atom b-term response weight u_i = z_q(i)·κ_i/(2√CN_i)·g_i, g_i the
// log-compression factor 1/(1+e^(cn_raw_i−cnmax)). (d4_charge_model.cpp:116-122).
__global__ void k_d4eeq_u(int N, const double* __restrict__ zq, const double* __restrict__ cnf,
                          const double* __restrict__ cn, const double* __restrict__ cn_raw,
                          double* __restrict__ u)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    const double cni = cn[i];
    if (cni <= D4EEQ_CN_EPS || cnf[i] == 0.0) { u[i] = 0.0; return; }
    const double g = 1.0 / (1.0 + exp(cn_raw[i] - D4EEQ_CNMAX));
    u[i] = zq[i] * cnf[i] / (2.0 * sqrt(cni)) * g;
}

// Charge-response gradient, one thread per atom a (no atomicAdd): the CPU pair
// loop (b<a) gives a += g_pair and b −= g_pair; summing the full b≠a loop per
// atom reproduces both halves since coeff is pair-symmetric and û flips sign.
// grad layout [3a+k], Eh/Bohr. Mirrors addChargeResponseGradient (lines 124-152).
__global__ void k_d4eeq_response(int N, const double* __restrict__ xyz,
                                 const double* __restrict__ alpha_sq,
                                 const double* __restrict__ rcov_bohr,
                                 const double* __restrict__ q, const double* __restrict__ zq,
                                 const double* __restrict__ u, double* __restrict__ grad)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= N) return;
    const double xa = xyz[3 * a + 0], ya = xyz[3 * a + 1], za = xyz[3 * a + 2];
    const double alp_a = alpha_sq[a], rc_a = rcov_bohr[a];
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int b = 0; b < N; ++b) {
        if (b == a) continue;
        const double dx = xa - xyz[3 * b + 0];
        const double dy = ya - xyz[3 * b + 1];
        const double dz = za - xyz[3 * b + 2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);
        if (r < 1.0e-10) continue;
        const double inv_r = 1.0 / r;
        // A-term: −c_ab·A'(r)·û,  c_ab = z_q(a)q_b + z_q(b)q_a.
        const double gammij = 1.0 / sqrt(alp_a + alpha_sq[b]);
        const double gr = gammij * r;
        const double Aprime = gammij * D4EEQ_TWO_OVER_SQRTPI * exp(-gr * gr) * inv_r
                              - erf(gr) * inv_r * inv_r;
        const double c_ab = zq[a] * q[b] + zq[b] * q[a];
        // b-term (CN): (u_a+u_b)·d(cn_raw)/dr·û.
        double Draw = 0.0;
        const double rc_b = rcov_bohr[b];
        if (rc_a > 0.0 && rc_b > 0.0) {
            const double rcij = rc_a + rc_b;
            const double dr = (r - rcij) / rcij;
            const double earg = D4EEQ_KN * dr;
            Draw = 0.5 * D4EEQ_TWO_OVER_SQRTPI * exp(-earg * earg) * D4EEQ_KN / rcij;
        }
        const double coeff = (u[a] + u[b]) * Draw - c_ab * Aprime;
        gx += coeff * dx * inv_r;
        gy += coeff * dy * inv_r;
        gz += coeff * dz * inv_r;
    }
    grad[3 * a + 0] = gx;
    grad[3 * a + 1] = gy;
    grad[3 * a + 2] = gz;
}

// Stage 5 (Part B1): atomic Mulliken charges from the resident AO populations.
// q_at is pre-seeded to n0_at; one thread per AO subtracts its population into
// its atom bin. Mirrors updatePopulationsFromPopAo (xtb_scf.cpp): n_at(A) =
// Σ_{μ∈A} pop_ao(μ), q_at(A) = n0_at(A) − n_at(A). Claude Generated.
__global__ void k_qat_scatter(int nao, const double* __restrict__ pop,
                              const int* __restrict__ ao2at, double* __restrict__ q_at)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    atomicAdd(&q_at[ao2at[mu]], -pop[mu]);
}

// Stage 6 (S6.2): shell Mulliken charges from the resident AO populations. q_sh is
// pre-seeded to n0_sh; one thread per AO subtracts its population into its shell
// bin. The shell half of updatePopulationsFromPopAo (xtb_scf.cpp): n_sh(s) =
// Σ_{μ∈s} pop_ao(μ), q_sh(s) = n0_sh(s) − n_sh(s). Claude Generated.
__global__ void k_qsh_scatter(int nao, const double* __restrict__ pop,
                              const int* __restrict__ ao2sh, double* __restrict__ q_sh)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    atomicAdd(&q_sh[ao2sh[mu]], -pop[mu]);
}

// Stage 6 (S6.2b): device ports of d4_zeta / d4_dzeta (d4_charge_scaling.h).
__device__ __forceinline__ double d4_zeta_dev(double a, double c, double qref, double qmod)
{
    if (qmod < 0.0) return exp(a);
    return exp(a * (1.0 - exp(c * (1.0 - qref / qmod))));
}
__device__ __forceinline__ double d4_dzeta_dev(double a, double c, double qref, double qmod)
{
    if (qmod < 0.0) return 0.0;
    const double g = exp(c * (1.0 - qref / qmod));
    const double z = exp(a * (1.0 - g));
    return -a * c * g * z * qref / (qmod * qmod);
}

// Stage 6 (S6.2b): rebuild the per-atom D4 reference weights W = gwk(CN)·ζ(q) and
// dWq = ∂W/∂q from the SCF charges on the device — the exact port of
// D4ParameterGenerator::buildAtomRefW (want_grad). One thread per atom; the
// q-independent tables (cn, gi=eta·gc, zeff, refcn, refcovcn, refq, nref) are
// uploaded once per geometry, q is read from the resident charges, and the
// outputs feed k_d4_dedq with no host round-trip. ga=3, wf=6, MAXCN=19 (dftd4
// defaults); the ngw bucketing uses refcn, the CN-Gaussian uses refcovcn.
// Claude Generated.
__global__ void k_d4_build_refw(int nat, int max_ref, const double* __restrict__ q,
                                const double* __restrict__ cn, const double* __restrict__ gi,
                                const double* __restrict__ zeff, const int* __restrict__ nref,
                                const double* __restrict__ refcn, const double* __restrict__ refcovcn,
                                const double* __restrict__ refq,
                                double* __restrict__ W, double* __restrict__ dWq)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    constexpr double ga = 3.0, wf = 6.0;
    constexpr int MAXCN = 19;
    const int MR = max_ref;                       // 7
    const size_t base = static_cast<size_t>(a) * MR;
    for (int ir = 0; ir < MR; ++ir) { W[base + ir] = 0.0; dWq[base + ir] = 0.0; }
    int nr = nref[a];
    if (nr <= 0) return;
    if (nr > MR) nr = MR;
    const double cna = cn[a];
    const double gia = gi[a];
    const double zef = zeff[a];

    // ngw[ir] from refcn (dftd4 set_refgw): count refs sharing a rounded CN bucket.
    int cnc[MAXCN + 1];
    for (int k = 0; k <= MAXCN; ++k) cnc[k] = 0;
    cnc[0] = 1;
    for (int ir = 0; ir < nr; ++ir) {
        int icn = static_cast<int>(lround(refcn[base + ir]));
        if (icn < 0) icn = 0; if (icn > MAXCN) icn = MAXCN;
        cnc[icn] += 1;
    }
    // CN-Gaussian weights gwk = expw/norm (uses refcovcn, dftd4 weight_references).
    double expw[7];
    double norm = 0.0;
    for (int ir = 0; ir < nr; ++ir) {
        int icn = static_cast<int>(lround(refcn[base + ir]));
        if (icn < 0) icn = 0; if (icn > MAXCN) icn = MAXCN;
        const int k = cnc[icn];
        const int ngw = k * (k + 1) / 2;
        const double dcov = cna - refcovcn[base + ir];
        double ew = 0.0;
        for (int igw = 1; igw <= ngw; ++igw) {
            const double wfe = igw * wf;
            ew += exp(-wfe * dcov * dcov);
        }
        expw[ir] = ew;
        norm += ew;
    }
    const double ninv = (norm > 0.0) ? 1.0 / norm : 0.0;
    const double qmod = q[a] + zef;
    for (int ir = 0; ir < nr; ++ir) {
        double gwk = expw[ir] * ninv;
        if (!isfinite(gwk)) gwk = 0.0;
        const double qref = refq[base + ir] + zef;
        W[base + ir]   = gwk * d4_zeta_dev(ga, gia, qref, qmod);
        dWq[base + ir] = gwk * d4_dzeta_dev(ga, gia, qref, qmod);
    }
}

// Stage 5 (Part B2): in-SCF GFN2 D4 atom-potential dE_D4/dq, one thread per atom
// i (no atomicAdd — like k_d4eeq_response, each atom sums its own half of the
// symmetric pair contributions). For every j≠i within the 50-Bohr cutoff:
//   r4r2ij = 3·sqrtZr4r2_i·sqrtZr4r2_j ; R0 = a1·√r4r2ij + a2
//   disp_sum = s6/(r⁶+R0⁶) + s8·r4r2ij/(r⁸+R0⁸)              (BJ, geometry-fixed)
//   dc6/dq_i = ΣΣ dWq_i[a]·W_j[b]·c6ref(ei,ej,a,b)            (7×7 contraction)
//   dEdq(i) += −dc6/dq_i · disp_sum
// Mirrors D4Evaluator::computeEnergyAndGradient's per-reference path + contractC6Gfn2
// (the c6 cache is element-symmetric, so atom j's thread yields the CPU's dc6dqj).
// MAX_ELEM=118, MAX_REF=7. Claude Generated.
__global__ void k_d4_dedq(int nat, int max_elem, int max_ref,
                          const int* __restrict__ Z, const int* __restrict__ nref,
                          const double* __restrict__ sqrtZr4r2, const double* __restrict__ xyz,
                          const double* __restrict__ c6_flat,
                          const double* __restrict__ W, const double* __restrict__ dWq,
                          double s6, double s8, double a1, double a2, double cut2,
                          double* __restrict__ dEdq)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int ei = Z[i] - 1;
    const int nri = nref[i];
    if (ei < 0 || ei >= max_elem || nri <= 0) { dEdq[i] = 0.0; return; }
    const double xi = xyz[3 * i + 0], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
    const double sq_i = sqrtZr4r2[i];
    const double* dWqi = dWq + static_cast<size_t>(i) * max_ref;

    double acc = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const double dx = xi - xyz[3 * j + 0];
        const double dy = yi - xyz[3 * j + 1];
        const double dz = zi - xyz[3 * j + 2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cut2 || r2 < 1.0e-20) continue;
        const int ej = Z[j] - 1;
        const int nrj = nref[j];
        if (ej < 0 || ej >= max_elem || nrj <= 0) continue;

        const double r4r2ij = 3.0 * sq_i * sqrtZr4r2[j];
        const double r0 = a1 * sqrt(r4r2ij) + a2;
        const double r0_2 = r0 * r0;
        const double r0_6 = r0_2 * r0_2 * r0_2;
        const double r0_8 = r0_6 * r0_2;
        const double r6 = r2 * r2 * r2;
        const double r8 = r6 * r2;
        const double t6 = 1.0 / (r6 + r0_6);
        const double t8 = 1.0 / (r8 + r0_8);
        const double disp_sum = s6 * t6 + s8 * r4r2ij * t8;

        const double* Wj = W + static_cast<size_t>(j) * max_ref;
        const size_t base = (static_cast<size_t>(ei) * max_elem + ej)
                          * static_cast<size_t>(max_ref) * max_ref;
        double dc6dqi = 0.0;
        for (int a = 0; a < nri; ++a) {
            const double dwa = dWqi[a];
            if (dwa == 0.0) continue;
            const size_t basea = base + static_cast<size_t>(a) * max_ref;
            double s = 0.0;
            for (int b = 0; b < nrj; ++b) s += Wj[b] * c6_flat[basea + b];
            dc6dqi += dwa * s;
        }
        acc += -dc6dqi * disp_sum;
    }
    dEdq[i] = acc;
}

// ---- Stage 5 (Part B3/B4): full device GFN2 potential build kernels ---------

// q_at(A) = Σ_{s∈A} q_sh(s). One thread per shell, atomicAdd into the (zeroed)
// atom bin. (q_at = n0_at − n_at = Σ_{s∈A}(n0_sh − n_sh) = Σ_{s∈A} q_sh.)
__global__ void k_qsh_to_qat(int nsh, const double* __restrict__ q_sh,
                             const int* __restrict__ sh2at, double* __restrict__ q_at)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    atomicAdd(&q_at[sh2at[s]], q_sh[s]);
}

// GFN2 shell third-order: v_sh(s) += q_sh(s)²·Γ_s (added onto the γ·q_sh gemv
// result already in v_sh). Mirrors addThirdOrderPotential (shell-resolved).
__global__ void k_vsh_third(int nsh, const double* __restrict__ q_sh,
                            const double* __restrict__ gamma3, double* __restrict__ v_sh)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    v_sh[s] += q_sh[s] * q_sh[s] * gamma3[s];
}

// GFN2 multipole potential (one thread per atom i) — port of addMultipolePotential
// (xtb_multipole.cpp). amat_* are column-major nat×nat blocks; dp_at/qp_at are the
// uploaded mixed atomic moments ([k+j*3] / [k+j*6]); q_at the device shell→atom sum.
//   vdp(k,i) = Σ_j amat_sd[k](i,j)·q_at(j) + Σ_a amat_dd[k][a](i,j)·dp_at(a,j)
//              + 2·dkernel(i)·dp_at(k,i)
//   vqp(k,i) = Σ_j amat_sq[k](i,j)·q_at(j) + 2·qkernel(i)·qp_at(k,i)·mpscale_q[k]
//   v_at(i)  = Σ_j [Σ_k amat_sd[k](j,i)·dp_at(k,j) + Σ_k amat_sq[k](j,i)·qp_at(k,j)]
__global__ void k_multipole_potential(
    int nat, const double* __restrict__ amat_sd, const double* __restrict__ amat_dd,
    const double* __restrict__ amat_sq, const double* __restrict__ dkernel,
    const double* __restrict__ qkernel, const double* __restrict__ q_at,
    const double* __restrict__ dp_at, const double* __restrict__ qp_at,
    double* __restrict__ v_dp, double* __restrict__ v_qp, double* __restrict__ v_at)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
    const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
    double vd[3] = {0.0, 0.0, 0.0};
    double vq[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double vat = 0.0;
    for (int j = 0; j < nat; ++j) {
        const double qj = q_at[j];
        const size_t ij = static_cast<size_t>(i) + static_cast<size_t>(j) * nat;  // (i,j)
        const size_t ji = static_cast<size_t>(j) + static_cast<size_t>(i) * nat;  // (j,i)
        const double dpj0 = dp_at[0 + j * 3], dpj1 = dp_at[1 + j * 3], dpj2 = dp_at[2 + j * 3];
        for (int k = 0; k < 3; ++k) {
            vd[k] += amat_sd[static_cast<size_t>(k) * nn + ij] * qj
                   + amat_dd[(static_cast<size_t>(k) * 3 + 0) * nn + ij] * dpj0
                   + amat_dd[(static_cast<size_t>(k) * 3 + 1) * nn + ij] * dpj1
                   + amat_dd[(static_cast<size_t>(k) * 3 + 2) * nn + ij] * dpj2;
            vat += amat_sd[static_cast<size_t>(k) * nn + ji] * dp_at[k + j * 3];
        }
        for (int k = 0; k < 6; ++k) {
            vq[k] += amat_sq[static_cast<size_t>(k) * nn + ij] * qj;
            vat += amat_sq[static_cast<size_t>(k) * nn + ji] * qp_at[k + j * 6];
        }
    }
    for (int k = 0; k < 3; ++k) v_dp[k + i * 3] = vd[k] + 2.0 * dkernel[i] * dp_at[k + i * 3];
    for (int k = 0; k < 6; ++k)
        v_qp[k + i * 6] = vq[k] + 2.0 * qkernel[i] * qp_at[k + i * 6] * mpscale_q[k];
    v_at[i] = vat;
}

// Claude Generated (Sep 2026, large systems): on-the-fly twins of k_multipole_potential and
// k_energy_multipole. The GFN2 multipole interaction matrices are 18 nat^2 doubles (7.7 GB at
// 7320 atoms) but each element is a closed-form function of one atom pair (xtb_multipole.cpp
// step 5), so the kernels rebuild the elements they need instead of reading stored matrices.
// amat(row t, col s) uses v = x_s - x_t, r = |v|, rr = (mrad_t + mrad_s)/2 / r:
//   sd[k]    = v_k g3 fdmp3
//   dd[a][b] = delta_ab g3 fdmp5 - v_a v_b 3 g5 fdmp5
//   sq       = (vx^2, 2 vx vy, vy^2, 2 vx vz, 2 vy vz, vz^2) g5 fdmp5
__device__ __forceinline__ void d_mp_amat_pair(double vx, double vy, double vz, double rad_sum_half,
                                         double dmp3, double dmp5,
                                         double sd[3], double dd[9], double sq[6])
{
    const double r1 = sqrt(vx*vx + vy*vy + vz*vz);
    const double g1 = 1.0 / r1;
    const double g3 = g1 * g1 * g1;
    const double g5 = g3 * g1 * g1;
    const double rr = rad_sum_half * g1;
    const double fdmp3 = 1.0 / (1.0 + 6.0 * pow(rr, dmp3));
    const double fdmp5 = 1.0 / (1.0 + 6.0 * pow(rr, dmp5));
    sd[0] = vx * g3 * fdmp3; sd[1] = vy * g3 * fdmp3; sd[2] = vz * g3 * fdmp3;
    const double dd_iso = g3 * fdmp5, dd_anis = 3.0 * g5 * fdmp5;
    const double v[3] = {vx, vy, vz};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            dd[a * 3 + b] = ((a == b) ? dd_iso : 0.0) - v[a] * v[b] * dd_anis;
    sq[0] = vx * vx * g5 * fdmp5;
    sq[1] = 2.0 * vx * vy * g5 * fdmp5;
    sq[2] = vy * vy * g5 * fdmp5;
    sq[3] = 2.0 * vx * vz * g5 * fdmp5;
    sq[4] = 2.0 * vy * vz * g5 * fdmp5;
    sq[5] = vz * vz * g5 * fdmp5;
}

__global__ void k_multipole_potential_otf(
    int nat, const double* __restrict__ xyz, const double* __restrict__ mrad, double dmp3, double dmp5,
    const double* __restrict__ dkernel, const double* __restrict__ qkernel,
    const double* __restrict__ q_at, const double* __restrict__ dp_at, const double* __restrict__ qp_at,
    double* __restrict__ v_dp, double* __restrict__ v_qp, double* __restrict__ v_at)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
    double vd[3] = {0.0, 0.0, 0.0};
    double vq[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double vat = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;   // amat diagonal is zero
        const double qj = q_at[j];
        const double dpj0 = dp_at[0 + j * 3], dpj1 = dp_at[1 + j * 3], dpj2 = dp_at[2 + j * 3];
        const double half = 0.5 * (mrad[i] + mrad[j]);
        double sd_ij[3], dd_ij[9], sq_ij[6], sd_ji[3], dd_ji[9], sq_ji[6];
        // (row i, col j): v = x_j - x_i ; (row j, col i): v = x_i - x_j
        d_mp_amat_pair(xyz[3*j] - xyz[3*i], xyz[3*j+1] - xyz[3*i+1], xyz[3*j+2] - xyz[3*i+2],
                       half, dmp3, dmp5, sd_ij, dd_ij, sq_ij);
        d_mp_amat_pair(xyz[3*i] - xyz[3*j], xyz[3*i+1] - xyz[3*j+1], xyz[3*i+2] - xyz[3*j+2],
                       half, dmp3, dmp5, sd_ji, dd_ji, sq_ji);
        for (int k = 0; k < 3; ++k) {
            vd[k] += sd_ij[k] * qj
                   + dd_ij[k * 3 + 0] * dpj0
                   + dd_ij[k * 3 + 1] * dpj1
                   + dd_ij[k * 3 + 2] * dpj2;
            vat += sd_ji[k] * dp_at[k + j * 3];
        }
        for (int k = 0; k < 6; ++k) {
            vq[k] += sq_ij[k] * qj;
            vat += sq_ji[k] * qp_at[k + j * 6];
        }
    }
    for (int k = 0; k < 3; ++k) v_dp[k + i * 3] = vd[k] + 2.0 * dkernel[i] * dp_at[k + i * 3];
    for (int k = 0; k < 6; ++k)
        v_qp[k + i * 6] = vq[k] + 2.0 * qkernel[i] * qp_at[k + i * 6] * mpscale_q[k];
    v_at[i] = vat;
}

__global__ void k_energy_multipole_otf(int nat, const double* __restrict__ xyz, const double* __restrict__ mrad,
                                       double dmp3, double dmp5,
                                       const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                       const double* __restrict__ dp_at, const double* __restrict__ qp_at,
                                       const double* __restrict__ q_at, double* e_out)
{
    extern __shared__ double sdata[];
    const int tid = threadIdx.x;
    const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
    const int i = blockIdx.x * blockDim.x + tid;
    double ei = 0.0;
    if (i < nat) {
        double dpi[3], qpi[6];
        for (int k = 0; k < 3; ++k) dpi[k] = dp_at[k + static_cast<size_t>(i) * 3];
        for (int k = 0; k < 6; ++k) qpi[k] = qp_at[k + static_cast<size_t>(i) * 6];
        for (int j = 0; j < nat; ++j) {
            if (j == i) continue;
            const double qj = q_at[j];
            double sd[3], dd[9], sq[6];
            d_mp_amat_pair(xyz[3*j] - xyz[3*i], xyz[3*j+1] - xyz[3*i+1], xyz[3*j+2] - xyz[3*i+2],
                           0.5 * (mrad[i] + mrad[j]), dmp3, dmp5, sd, dd, sq);
            for (int k = 0; k < 3; ++k)
                ei += dpi[k] * sd[k] * qj;
            for (int a = 0; a < 3; ++a) {
                const double dpia = dpi[a];
                for (int b = 0; b < 3; ++b)
                    ei += 0.5 * dpia * dd[a * 3 + b] * dp_at[b + static_cast<size_t>(j) * 3];
            }
            for (int k = 0; k < 6; ++k)
                ei += qpi[k] * sq[k] * qj;
        }
        const double dk = dkernel[i], qk = qkernel[i];
        for (int k = 0; k < 3; ++k) ei += dk * dpi[k] * dpi[k];
        for (int k = 0; k < 6; ++k) ei += qk * qpi[k] * qpi[k] * mpscale_q[k];
    }
    sdata[tid] = ei; __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) { if (tid < st) sdata[tid] += sdata[tid + st]; __syncthreads(); }
    if (tid == 0) atomicAdd(e_out, sdata[0]);
}

// ---------------------------------------------------------------------------------------------
// Claude Generated (Sep 29, 2026): warp-per-atom variants of the three O(N^2)/O(N^3) atom-pair
// kernels. The thread-per-atom originals launch only nat threads (polymer_2x: 7320 threads, i.e.
// ~5 warps per SM on an A4500), so the FP64 pair loops ran latency-bound: k_d4_atm_nl 2.97 s,
// k_multipole_potential_otf 0.29 s per SCF iteration. Here a warp owns atom i, its 32 lanes stride
// over the partner atoms and a fixed shuffle tree reduces the partial sums - deterministic, but in
// a different summation order than the originals (last-bit differences). The multipole potential
// additionally evaluates the pair interaction once: for v -> -v the damping (|v| only) and the
// dipole-dipole and quadrupole tensors are unchanged and sd flips sign, all exactly in floating
// point, so the (j,i) terms are taken from the (i,j) evaluation bit for bit.
// CURCUMA_GPU_THREAD_PER_ATOM=1 selects the originals (for comparison).
// ---------------------------------------------------------------------------------------------
__device__ __forceinline__ double d_warp_sum(double v)
{
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

__global__ void k_multipole_potential_otf_w(
    int nat, const double* __restrict__ xyz, const double* __restrict__ mrad, double dmp3, double dmp5,
    const double* __restrict__ dkernel, const double* __restrict__ qkernel,
    const double* __restrict__ q_at, const double* __restrict__ dp_at, const double* __restrict__ qp_at,
    double* __restrict__ v_dp, double* __restrict__ v_qp, double* __restrict__ v_at)
{
    const int i = static_cast<int>((static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x) >> 5);
    const int lane = threadIdx.x & 31;
    if (i >= nat) return;   // uniform per warp
    const double xi = xyz[3*i], yi = xyz[3*i+1], zi = xyz[3*i+2];
    const double mi = mrad[i];
    double vd[3] = {0.0, 0.0, 0.0};
    double vq[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double vat = 0.0;
    for (int j = lane; j < nat; j += 32) {
        if (j == i) continue;   // amat diagonal is zero
        const double qj = q_at[j];
        const double dpj0 = dp_at[0 + j * 3], dpj1 = dp_at[1 + j * 3], dpj2 = dp_at[2 + j * 3];
        double sd[3], dd[9], sq[6];
        d_mp_amat_pair(xyz[3*j] - xi, xyz[3*j+1] - yi, xyz[3*j+2] - zi,
                       0.5 * (mi + mrad[j]), dmp3, dmp5, sd, dd, sq);
        for (int k = 0; k < 3; ++k) {
            vd[k] += sd[k] * qj + dd[k * 3 + 0] * dpj0 + dd[k * 3 + 1] * dpj1 + dd[k * 3 + 2] * dpj2;
            vat += (-sd[k]) * dp_at[k + j * 3];          // sd(j,i) = -sd(i,j)
        }
        for (int k = 0; k < 6; ++k) {
            vq[k] += sq[k] * qj;
            vat += sq[k] * qp_at[k + j * 6];             // sq(j,i) = sq(i,j)
        }
    }
    for (int k = 0; k < 3; ++k) vd[k] = d_warp_sum(vd[k]);
    for (int k = 0; k < 6; ++k) vq[k] = d_warp_sum(vq[k]);
    vat = d_warp_sum(vat);
    if (lane != 0) return;
    const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
    for (int k = 0; k < 3; ++k) v_dp[k + i * 3] = vd[k] + 2.0 * dkernel[i] * dp_at[k + i * 3];
    for (int k = 0; k < 6; ++k)
        v_qp[k + i * 6] = vq[k] + 2.0 * qkernel[i] * qp_at[k + i * 6] * mpscale_q[k];
    v_at[i] = vat;
}

// Energy twin: one warp per atom i, then the block's warp sums are added and one atomicAdd per
// block goes to e_out (as in the original). blockDim must be a multiple of 32, shared = warps.
__global__ void k_energy_multipole_otf_w(int nat, const double* __restrict__ xyz, const double* __restrict__ mrad,
                                         double dmp3, double dmp5,
                                         const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                         const double* __restrict__ dp_at, const double* __restrict__ qp_at,
                                         const double* __restrict__ q_at, double* e_out)
{
    extern __shared__ double sdata[];
    const int i = static_cast<int>((static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x) >> 5);
    const int lane = threadIdx.x & 31, wib = threadIdx.x >> 5;
    double ei = 0.0;
    if (i < nat) {
        const double xi = xyz[3*i], yi = xyz[3*i+1], zi = xyz[3*i+2];
        double dpi[3], qpi[6];
        for (int k = 0; k < 3; ++k) dpi[k] = dp_at[k + static_cast<size_t>(i) * 3];
        for (int k = 0; k < 6; ++k) qpi[k] = qp_at[k + static_cast<size_t>(i) * 6];
        for (int j = lane; j < nat; j += 32) {
            if (j == i) continue;
            const double qj = q_at[j];
            double sd[3], dd[9], sq[6];
            d_mp_amat_pair(xyz[3*j] - xi, xyz[3*j+1] - yi, xyz[3*j+2] - zi,
                           0.5 * (mrad[i] + mrad[j]), dmp3, dmp5, sd, dd, sq);
            for (int k = 0; k < 3; ++k)
                ei += dpi[k] * sd[k] * qj;
            for (int a = 0; a < 3; ++a) {
                const double dpia = dpi[a];
                for (int b = 0; b < 3; ++b)
                    ei += 0.5 * dpia * dd[a * 3 + b] * dp_at[b + static_cast<size_t>(j) * 3];
            }
            for (int k = 0; k < 6; ++k)
                ei += qpi[k] * sq[k] * qj;
        }
        ei = d_warp_sum(ei);
        if (lane == 0) {
            const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
            const double dk = dkernel[i], qk = qkernel[i];
            for (int k = 0; k < 3; ++k) ei += dk * dpi[k] * dpi[k];
            for (int k = 0; k < 6; ++k) ei += qk * qpi[k] * qpi[k] * mpscale_q[k];
        }
    }
    if (lane == 0) sdata[wib] = ei;
    __syncthreads();
    if (threadIdx.x == 0) {
        double s = 0.0;
        for (int w = 0; w < static_cast<int>(blockDim.x >> 5); ++w) s += sdata[w];
        atomicAdd(e_out, s);
    }
}

// Warp-per-atom twin of k_d4_atm_nl (see the note at k_multipole_potential_otf_w): lanes stride
// over the neighbour index ix, each keeps the inner y < x loop of the original.
__global__ void k_d4_atm_nl_w(int nat, const double* __restrict__ xyz, const double* __restrict__ r4r2,
                              const double* __restrict__ c6, const double* __restrict__ dc6dcn,
                              const int* __restrict__ nbptr, const int* __restrict__ nb,
                              double s9, double a1, double a2, double alp, double cut2,
                              double* __restrict__ e_atom, double* __restrict__ grad,
                              double* __restrict__ dEdcn)
{
    const int a = static_cast<int>((static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x) >> 5);
    const int lane = threadIdx.x & 31;
    if (a >= nat) return;   // uniform per warp
    const double eps = 2.220446049250313e-16;
    const double xa = xyz[3*a+0], ya = xyz[3*a+1], za = xyz[3*a+2];
    const double r4r2a = r4r2[a];
    const size_t nat_s = static_cast<size_t>(nat);
    const bool alp16 = (alp == 16.0);
    double eacc = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, cnacc = 0.0;
    for (int ix = nbptr[a] + lane; ix < nbptr[a + 1]; ix += 32) {
        const int x = nb[ix];
        const double vaxx = xyz[3*x+0]-xa, vaxy = xyz[3*x+1]-ya, vaxz = xyz[3*x+2]-za;
        const double r2ax = vaxx*vaxx + vaxy*vaxy + vaxz*vaxz;
        if (r2ax > cut2 || r2ax < eps) continue;
        const double c6ax = c6[a*nat_s + x];
        const double r0ax = a1*sqrt(3.0*r4r2a*r4r2[x]) + a2;
        for (int iy = nbptr[a]; iy < ix; ++iy) {        // y < x (list ascending)
            const int y = nb[iy];
            const double vayx = xyz[3*y+0]-xa, vayy = xyz[3*y+1]-ya, vayz = xyz[3*y+2]-za;
            const double r2ay = vayx*vayx + vayy*vayy + vayz*vayz;
            if (r2ay > cut2 || r2ay < eps) continue;
            const double dxyx = xyz[3*y+0]-xyz[3*x+0], dxyy = xyz[3*y+1]-xyz[3*x+1], dxyz_ = xyz[3*y+2]-xyz[3*x+2];
            const double r2xy = dxyx*dxyx + dxyy*dxyy + dxyz_*dxyz_;
            if (r2xy > cut2 || r2xy < eps) continue;

            const double c6ay = c6[a*nat_s + y];
            const double c6xy = c6[x*nat_s + y];
            const double c9 = -s9 * sqrt(fabs(c6ax*c6ay*c6xy));
            const double r0ay = a1*sqrt(3.0*r4r2a*r4r2[y]) + a2;
            const double r0xy = a1*sqrt(3.0*r4r2[x]*r4r2[y]) + a2;
            const double r0 = r0ax*r0ay*r0xy;
            const double r2 = r2ax*r2ay*r2xy;
            const double r1 = sqrt(r2);
            const double r3 = r2*r1;
            const double r5 = r3*r2;
            const double q = r0/r1;
            const double pw = alp16 ? (q*q*q*q*q) * cbrt(q) : pow(q, alp/3.0);
            const double fdmp = 1.0/(1.0 + 6.0*pw);
            const double ang = 0.375*(r2ax + r2xy - r2ay)*(r2ax - r2xy + r2ay)*(-r2ax + r2xy + r2ay)/r5 + 1.0/r3;
            const double dE = ang*fdmp*c9;
            eacc += -dE;

            const double dfdmp = -2.0*alp*pw*fdmp*fdmp;
            const double dang_ax = -0.375*(r2ax*r2ax*r2ax + r2ax*r2ax*(r2xy+r2ay)
                          + r2ax*(3.0*r2xy*r2xy + 2.0*r2xy*r2ay + 3.0*r2ay*r2ay)
                          - 5.0*(r2xy-r2ay)*(r2xy-r2ay)*(r2xy+r2ay))/r5;
            const double gcax = c9*(-dang_ax*fdmp + ang*dfdmp)/r2ax;
            gx += -gcax*vaxx; gy += -gcax*vaxy; gz += -gcax*vaxz;
            const double dang_ay = -0.375*(r2ay*r2ay*r2ay + r2ay*r2ay*(r2xy+r2ax)
                          + r2ay*(3.0*r2xy*r2xy + 2.0*r2xy*r2ax + 3.0*r2ax*r2ax)
                          - 5.0*(r2xy-r2ax)*(r2xy-r2ax)*(r2xy+r2ax))/r5;
            const double gcay = c9*(-dang_ay*fdmp + ang*dfdmp)/r2ay;
            gx += -gcay*vayx; gy += -gcay*vayy; gz += -gcay*vayz;
            cnacc += -dE*0.5*(dc6dcn[a*nat_s+x]/c6ax + dc6dcn[a*nat_s+y]/c6ay);
        }
    }
    eacc = d_warp_sum(eacc);
    gx = d_warp_sum(gx); gy = d_warp_sum(gy); gz = d_warp_sum(gz);
    cnacc = d_warp_sum(cnacc);
    if (lane != 0) return;
    e_atom[a] = eacc;
    grad[3*a+0] = gx; grad[3*a+1] = gy; grad[3*a+2] = gz;
    dEdcn[a] = cnacc;
}

// v_at(A) += dE_D4/dq(A) (resident from k_d4_dedq). One thread per atom.
__global__ void k_vat_add_d4(int nat, const double* __restrict__ d4_dedq, double* __restrict__ v_at)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    v_at[a] += d4_dedq[a];
}

// Expand the shell+atom potential to AO resolution: v_ao(μ)=v_sh(ao2sh[μ])+v_at(ao2at[μ]).
__global__ void k_expand_vao(int nao, const double* __restrict__ v_sh, const double* __restrict__ v_at,
                             const int* __restrict__ ao2sh, const int* __restrict__ ao2at,
                             double* __restrict__ v_ao)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    v_ao[mu] = v_sh[ao2sh[mu]] + v_at[ao2at[mu]];
}

// ---- Stage 6 (S6.1) occupation kernel -------------------------------------
// Single-block port of XTB::occupationsFromEps (xtb_scf.cpp). One block, blockDim
// a power of two; eps/occ are length n (grid-stride within the block). For Tele=0
// the closed-shell integer fill (2.0 per occupied pair) is exact; for Tele>0 the
// Fermi level is found by bisection over [eps_min-1, eps_max+1] — 100 iterations,
// the same x<=500 exp clamp and 1e-14 width break as the host. The block-tree
// electron-count sum differs from the host's sequential sum only in rounding, so
// µ and occ agree to ~1e-13 (the bisection is self-correcting). mu_out/ncol_out
// optional (component test). Claude Generated.
__global__ void k_occupations(const double* __restrict__ eps, double* __restrict__ occ,
                              int n, double kT, double n_elec, int nocc_pairs,
                              int use_fermi, double* mu_out, int* ncol_out)
{
    extern __shared__ double sdata[];   // blockDim doubles
    const int tid = threadIdx.x;
    const int nthreads = blockDim.x;

    if (!use_fermi) {
        for (int i = tid; i < n; i += nthreads)
            occ[i] = (i < nocc_pairs) ? 2.0 : 0.0;
        if (tid == 0) { if (mu_out) *mu_out = 0.0; if (ncol_out) *ncol_out = nocc_pairs; }
        return;
    }

    // eps_min / eps_max via block reduction.
    double vmin = 1e300, vmax = -1e300;
    for (int i = tid; i < n; i += nthreads) {
        const double e = eps[i];
        vmin = fmin(vmin, e); vmax = fmax(vmax, e);
    }
    sdata[tid] = vmin; __syncthreads();
    for (int s = nthreads >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmin(sdata[tid], sdata[tid + s]); __syncthreads(); }
    const double eps_min = sdata[0]; __syncthreads();
    sdata[tid] = vmax; __syncthreads();
    for (int s = nthreads >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    const double eps_max = sdata[0]; __syncthreads();

    __shared__ double mu_lo_s, mu_hi_s;
    if (tid == 0) { mu_lo_s = eps_min - 1.0; mu_hi_s = eps_max + 1.0; }
    __syncthreads();

    for (int bisect = 0; bisect < 100; ++bisect) {
        const double mu = 0.5 * (mu_lo_s + mu_hi_s);
        double psum = 0.0;
        for (int i = tid; i < n; i += nthreads) {
            const double x = (eps[i] - mu) / kT;
            psum += 2.0 / (1.0 + exp(fmin(x, 500.0)));
        }
        sdata[tid] = psum; __syncthreads();
        for (int s = nthreads >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] += sdata[tid + s]; __syncthreads(); }
        const double n_sum = sdata[0];
        __syncthreads();
        if (tid == 0) { if (n_sum > n_elec) mu_hi_s = mu; else mu_lo_s = mu; }
        __syncthreads();
        if (mu_hi_s - mu_lo_s < 1e-14) break;
    }
    const double mu_f = 0.5 * (mu_lo_s + mu_hi_s);
    int local_ncol = 0;
    for (int i = tid; i < n; i += nthreads) {
        const double x = (eps[i] - mu_f) / kT;
        const double o = 2.0 / (1.0 + exp(fmin(x, 500.0)));
        occ[i] = o;
        if (o > 1.0e-12) local_ncol = i + 1;
    }
    __syncthreads();
    sdata[tid] = static_cast<double>(local_ncol); __syncthreads();
    for (int s = nthreads >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    if (tid == 0) { if (mu_out) *mu_out = mu_f; if (ncol_out) *ncol_out = static_cast<int>(sdata[0]); }
}

// ---- Stage 6 (S6.3) SCC energy kernels ------------------------------------
// GFN2 shell third-order energy E = Σ_s q_sh(s)³·Γ_s/3 (coulomb::energy_third_order,
// GFN2 branch). Grid-stride block reduction → atomicAdd into e_out (pre-zeroed).
// Claude Generated.
__global__ void k_energy_third_order_shell(int nsh, const double* __restrict__ q_sh,
                                           const double* __restrict__ gamma3, double* e_out)
{
    extern __shared__ double sdata[];
    const int tid = threadIdx.x;
    double e = 0.0;
    for (int s = blockIdx.x * blockDim.x + tid; s < nsh; s += gridDim.x * blockDim.x) {
        const double q = q_sh[s];
        e += q * q * q * gamma3[s] / 3.0;
    }
    sdata[tid] = e; __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) { if (tid < st) sdata[tid] += sdata[tid + st]; __syncthreads(); }
    if (tid == 0) atomicAdd(e_out, sdata[0]);
}

// GFN2 multipole energy E (XTB::energyMultipole): one thread per atom i sums its
// SD/DD/SQ row contractions over j plus the on-site dipole/quadrupole XC, then a
// block reduction → atomicAdd into e_out (pre-zeroed). amat layouts mirror
// k_multipole_potential exactly (block a*3+b for DD; col-major (i,j)=i+j·nat;
// moments dp_at[k+i·3], qp_at[k+i·6]). Uses the OUTPUT moments/charges. Claude Generated.
__global__ void k_energy_multipole(int nat, const double* __restrict__ amat_sd,
                                   const double* __restrict__ amat_dd, const double* __restrict__ amat_sq,
                                   const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                   const double* __restrict__ dp_at, const double* __restrict__ qp_at,
                                   const double* __restrict__ q_at, double* e_out)
{
    extern __shared__ double sdata[];
    const int tid = threadIdx.x;
    const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
    const int i = blockIdx.x * blockDim.x + tid;
    double ei = 0.0;
    if (i < nat) {
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        double dpi[3], qpi[6];
        for (int k = 0; k < 3; ++k) dpi[k] = dp_at[k + static_cast<size_t>(i) * 3];
        for (int k = 0; k < 6; ++k) qpi[k] = qp_at[k + static_cast<size_t>(i) * 6];
        for (int j = 0; j < nat; ++j) {
            const double qj = q_at[j];
            const size_t ij = static_cast<size_t>(i) + static_cast<size_t>(j) * nat;
            for (int k = 0; k < 3; ++k)
                ei += dpi[k] * amat_sd[static_cast<size_t>(k) * nn + ij] * qj;
            for (int a = 0; a < 3; ++a) {
                const double dpia = dpi[a];
                for (int b = 0; b < 3; ++b)
                    ei += 0.5 * dpia * amat_dd[(static_cast<size_t>(a) * 3 + b) * nn + ij]
                              * dp_at[b + static_cast<size_t>(j) * 3];
            }
            for (int k = 0; k < 6; ++k)
                ei += qpi[k] * amat_sq[static_cast<size_t>(k) * nn + ij] * qj;
        }
        const double dk = dkernel[i], qk = qkernel[i];
        for (int k = 0; k < 3; ++k) ei += dk * dpi[k] * dpi[k];
        for (int k = 0; k < 6; ++k) ei += qk * qpi[k] * qpi[k] * mpscale_q[k];
    }
    sdata[tid] = ei; __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) { if (tid < st) sdata[tid] += sdata[tid + st]; __syncthreads(); }
    if (tid == 0) atomicAdd(e_out, sdata[0]);
}

// ---- Stage 6 (S6.4) Broyden mixer kernels ---------------------------------
__global__ void k_vec_sub(double* __restrict__ out, const double* __restrict__ a,
                          const double* __restrict__ b, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = a[i] - b[i];
}

// Normalised residual change dF = dFraw/norm and the Broyden direction
// u = alpha·dF + (vin − vin_last)/norm, written into the history slot columns.
// Mirrors BroydenMixer::update lines 106-107. Claude Generated.
__global__ void k_broyden_dfu(double* __restrict__ dfcol, double* __restrict__ ucol,
                              const double* __restrict__ dfraw, const double* __restrict__ vin,
                              const double* __restrict__ vinlast, double alpha, double inv_norm, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double dfn = dfraw[i] * inv_norm;
    dfcol[i] = dfn;
    ucol[i]  = alpha * dfn + (vin[i] - vinlast[i]) * inv_norm;
}

// Single-thread solve of (w0²I + a)·gamma = c (M ≤ 20; a col-major M×M, symmetric):
// Gaussian elimination with partial pivoting → gamma = (w0²I+a)^{-1}·c, the host
// beta·c (BroydenMixer::update lines 133-135). Claude Generated.
__global__ void k_broyden_solve(int M, const double* __restrict__ a, const double* __restrict__ c,
                                double reg, double* __restrict__ gamma)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    const int MM = 20;
    double A[MM][MM];
    double b[MM];
    for (int i = 0; i < M; ++i) {
        b[i] = c[i];
        for (int j = 0; j < M; ++j) A[i][j] = a[i + j * M] + (i == j ? reg : 0.0);
    }
    for (int k = 0; k < M; ++k) {
        int piv = k; double mx = fabs(A[k][k]);
        for (int i = k + 1; i < M; ++i) { const double v = fabs(A[i][k]); if (v > mx) { mx = v; piv = i; } }
        if (piv != k) {
            for (int j = 0; j < M; ++j) { const double t = A[k][j]; A[k][j] = A[piv][j]; A[piv][j] = t; }
            const double t = b[k]; b[k] = b[piv]; b[piv] = t;
        }
        const double akk = A[k][k];
        for (int i = k + 1; i < M; ++i) {
            const double f = A[i][k] / akk;
            for (int j = k; j < M; ++j) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
    }
    for (int i = M - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < M; ++j) s -= A[i][j] * gamma[j];
        gamma[i] = s / A[i][i];
    }
}

// Stage 6 (S6.5): max|a[i] − b[i]| over n (the SCF convergence dq on q_sh).
// One block, grid-stride, tree max-reduction → out (1 double). Claude Generated.
__global__ void k_maxabsdiff(const double* __restrict__ a, const double* __restrict__ b,
                             int n, double* out)
{
    extern __shared__ double sdata[];
    const int tid = threadIdx.x;
    double m = 0.0;
    for (int i = tid; i < n; i += blockDim.x) m = fmax(m, fabs(a[i] - b[i]));
    sdata[tid] = m; __syncthreads();
    for (int s = blockDim.x >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    if (tid == 0) *out = sdata[0];
}

XtbGpuContext::XtbGpuContext(int device)
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
        return; // no device — caller falls back to CPU

    // Claude Generated (Sep 2026, multi-GPU): device < 0 keeps the historical behaviour
    // (whatever device is current on this thread, normally 0). An explicit index binds the
    // calling thread to that device BEFORE the stream and the cuBLAS/cuSOLVER handles are
    // created, because CUDA ties all of them to the device that is current at creation.
    if (device >= count)
        return; // invalid index — caller falls back to CPU (the adapter warns)
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess)
        return;
    if (cudaGetDevice(&m_impl->device) != cudaSuccess)
        return;

    cudaDeviceProp prop{};
    if (cudaGetDeviceProperties(&prop, m_impl->device) == cudaSuccess)
        m_impl->name = prop.name;

    if (cudaStreamCreate(&m_impl->stream) != cudaSuccess)
        return;

    if (cublasCreate(&m_impl->cublas) != CUBLAS_STATUS_SUCCESS)
        return;
    cublasSetStream(m_impl->cublas, m_impl->stream);

    if (cusolverDnCreate(&m_impl->cusolver) != CUSOLVER_STATUS_SUCCESS)
        return;
    cusolverDnSetStream(m_impl->cusolver, m_impl->stream);

    m_impl->ok = true;
}

XtbGpuContext::~XtbGpuContext()
{
    if (!m_impl)
        return;
    // The multi-GPU solver switches devices while it tears down; drop it first.
    m_impl->dist.reset();
    m_impl->releaseDensityHelpers();
    m_impl->releaseRemoteBackup();
    // Free the handles and every CudaBuffer (destroyed with m_impl) on OUR device.
    bindDevice();
    if (m_impl->cusolver) cusolverDnDestroy(m_impl->cusolver);
    if (m_impl->cublas)   cublasDestroy(m_impl->cublas);
    if (m_impl->stream)   cudaStreamDestroy(m_impl->stream);
}

bool XtbGpuContext::ok() const { return m_impl && m_impl->ok; }

bool XtbGpuContext::deviceHasFastFp64() const
{
    if (!m_impl || m_impl->device < 0) return false;
    cudaDeviceProp prop{};
    if (cudaGetDeviceProperties(&prop, m_impl->device) != cudaSuccess) return false;
    // FP64:FP32 is 1:2 on the datacenter parts and 1:32 / 1:64 elsewhere. By compute capability:
    // 6.0 P100, 7.0 V100, 8.0 A100, 9.x H100/H200/GH200, 10.x B200 - all 1:2. Consumer/workstation
    // (7.5, 8.6, 8.9, 12.x incl. the RTX PRO Blackwell parts) are not.
    const int cc = prop.major * 10 + prop.minor;
    return cc == 60 || cc == 70 || cc == 80 || prop.major == 9 || prop.major == 10;
}

std::string XtbGpuContext::deviceName() const
{
    return m_impl ? m_impl->name : std::string();
}

int XtbGpuContext::deviceId() const { return m_impl ? m_impl->device : -1; }

bool XtbGpuContext::bindDevice() const
{
    if (!m_impl || m_impl->device < 0) return false;
    int cur = -1;
    if (cudaGetDevice(&cur) == cudaSuccess && cur == m_impl->device) return true;
    return cudaSetDevice(m_impl->device) == cudaSuccess;
}

bool XtbGpuContext::deviceAvailable()
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

int XtbGpuContext::deviceCount()
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess ? count : 0;
}

bool XtbGpuContext::solveGeneralizedEigenF64(const double* F, const double* L, int n,
                                             double* C, double* eps)
{
    if (!ok() || n <= 0 || !F || !L || !C || !eps)
        return false;

    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);

    CudaBuffer<double> dA, dL, dW;
    try {
        dA.alloc(static_cast<int>(nn)); // Fock → Ã → C̃ → C (in place)
        dL.alloc(static_cast<int>(nn)); // lower Cholesky factor L (constant per geometry)
        dW.alloc(n);                    // eigenvalues
    } catch (...) {
        return false;
    }

    if (cudaMemcpyAsync(dA.ptr, F, sizeof(double) * nn, cudaMemcpyHostToDevice, stream) != cudaSuccess)
        return false;
    if (cudaMemcpyAsync(dL.ptr, L, sizeof(double) * nn, cudaMemcpyHostToDevice, stream) != cudaSuccess)
        return false;

    const double one = 1.0;

    // Step 1: dA <- F · L⁻ᵀ      (solve X·Lᵀ = F)
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one, dL.ptr, n, dA.ptr, n)
        != CUBLAS_STATUS_SUCCESS)
        return false;
    // Step 2: dA <- L⁻¹ · dA     (solve L·X = dA)  ⇒  Ã = L⁻¹·F·L⁻ᵀ
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &one, dL.ptr, n, dA.ptr, n)
        != CUBLAS_STATUS_SUCCESS)
        return false;

    // Standard symmetric eigensolve of Ã (lower): eigenvectors overwrite dA
    // (column-major), ascending eigenvalues into dW.
    int lwork = 0;
    if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                    CUBLAS_FILL_MODE_LOWER, n, dA.ptr, n, dW.ptr, &lwork)
        != CUSOLVER_STATUS_SUCCESS)
        return false;

    CudaBuffer<double> dWork;
    CudaBuffer<int>    dInfo;
    try {
        dWork.alloc(lwork > 0 ? lwork : 1);
        dInfo.alloc(1);
    } catch (...) {
        return false;
    }

    if (cusolverDnDsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                         n, dA.ptr, n, dW.ptr, dWork.ptr, lwork, dInfo.ptr)
        != CUSOLVER_STATUS_SUCCESS)
        return false;

    // Step 3: dA <- L⁻ᵀ · C̃      (solve Lᵀ·C = C̃)  ⇒  generalized eigenvectors C
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one, dL.ptr, n, dA.ptr, n)
        != CUBLAS_STATUS_SUCCESS)
        return false;

    int info = 1;
    if (cudaMemcpyAsync(&info, dInfo.ptr, sizeof(int), cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaMemcpyAsync(C, dA.ptr, sizeof(double) * nn, cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaMemcpyAsync(eps, dW.ptr, sizeof(double) * n, cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaStreamSynchronize(stream) != cudaSuccess)
        return false;

    return info == 0; // cusolver devInfo: 0 = success
}

/* ====================================================================== *
 *  Device-resident GFN1 SCF (Stage 2). H0/S/L stay on the GPU for the
 *  whole loop; only length-n vectors cross the bus per iteration.
 * ====================================================================== */

bool XtbGpuContext::residentBegin(const double* H0, const double* S,
                                  const double* L, int n)
{
    if (!ok() || n <= 0 || !H0 || !S || !L)
        return false;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    try {
        // Geometry-constant matrices, uploaded exactly once for this geometry.
        m_impl->sparse = false;   // host-uploaded integrals are always dense
        m_impl->dH0.upload(H0, static_cast<int>(nn), m_impl->stream);
        m_impl->dS.upload(S, static_cast<int>(nn), m_impl->stream);
        m_impl->dL.upload(L, static_cast<int>(nn), m_impl->stream);
        ++m_impl->l_generation;
        // Resident work buffers, reused across iterations.
        m_impl->dC.alloc(static_cast<int>(nn));
        m_impl->dP.alloc(static_cast<int>(nn));
        m_impl->dCw.alloc(static_cast<int>(nn));
        m_impl->dEps.alloc(n);
        m_impl->dVao.alloc(n);
        m_impl->dOcc.alloc(n);
        m_impl->dPop.alloc(n);
        // cuSOLVER dsyevd workspace (size is geometry-constant for fixed n).
        int lwork = 0;
        if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                        CUBLAS_FILL_MODE_LOWER, n, m_impl->dC.ptr, n,
                                        m_impl->dEps.ptr, &lwork)
            != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->lwork = lwork;
        m_impl->dWork.alloc(lwork > 0 ? lwork : 1);
        m_impl->dInfo.alloc(1);
    } catch (...) {
        return false;
    }
    m_impl->resident_n = n;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

// AP1 (Claude Generated): pad the host eigenvalue array beyond the partial window
// [0,neig) with an ascending sentinel safely above the highest computed eigenvalue.
// Those orbitals then carry occ≈0 in the host Fermi/integer occupation, so the
// density (lowest ncol≤neig columns) and the band/Fermi search are bit-faithful to
// a full solve whenever neig covers the occupied(+kT) window (validated by the caller).
static inline void fillEpsSentinel(double* eps_out, int neig, int n)
{
    const double sentinel = eps_out[neig - 1] + 1.0e3;  // Hartree; >> any valence eps
    for (int i = neig; i < n; ++i) eps_out[i] = sentinel;
}

// Reduce the Fock now in dC to standard form with the cached L, solve (dsyevd),
// back-transform → generalized eigenvectors in dC, eigenvalues → eps_out. The
// device analogue of the CPU dsygst+dsyevd+dtrsm path, reusing the resident L (no
// per-iteration L upload). Shared by residentSolve and residentSolveMultipole.
// ---- Stage 2: pseudo-diagonalisation (docs/GFN2_GPU_COST_PLAN.md) ----------------------------
// Claude Generated (Sep 29, 2026). Between two SCF steps the Fock matrix changes little, so the
// previous eigenvectors C (AO basis, C^T S C = I) are already close. One first-order Jacobi
// rotation between the occupied (o) and virtual (v) blocks removes the coupling
//     G = C^T F C,   theta_ai = G_ai / (G_ii - G_aa),   C_o += C_v theta,   C_v -= C_o theta^T
// (Stewart, Csaszar, Pulay, J. Comput. Chem. 3 (1982) 227), and a Cholesky re-orthonormalisation
// in the S metric (C^T S C = (L^T C)^T (L^T C) = R^T R, C <- C R^-1) restores C^T S C = I. The
// Cholesky order keeps the occupied span: column i only mixes with columns <= i. Only GEMMs, so it
// replaces an eigensolve whose tridiagonalisation is memory-bound. The orbital energies are the
// Rayleigh quotients G_pp of the rotated-from vectors (occupied block first, virtual block after -
// the occupation kernel does not need them sorted).

// g[p] = sum_mu C(mu,p) T(mu,p), accumulated in FP64; one block per column.
template <typename R>
__global__ void k_pseudo_coldot(const R* __restrict__ C, const R* __restrict__ T, int n,
                                double* __restrict__ g)
{
    extern __shared__ double sred[];
    const int p = blockIdx.x;
    const size_t off = static_cast<size_t>(p) * n;
    double acc = 0.0;
    for (int mu = threadIdx.x; mu < n; mu += blockDim.x)
        acc += static_cast<double>(C[off + mu]) * static_cast<double>(T[off + mu]);
    sred[threadIdx.x] = acc;
    __syncthreads();
    for (int s = blockDim.x >> 1; s > 0; s >>= 1) {
        if (threadIdx.x < s) sred[threadIdx.x] += sred[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) g[p] = sred[0];
}

// stats[0] = max_i g_i (occupied), stats[1] = min_a g_a (virtual). One block.
__global__ void k_pseudo_gap(const double* __restrict__ g, int no, int n, double* __restrict__ stats)
{
    extern __shared__ double sred[];
    double hi = -1e300, lo = 1e300;
    for (int p = threadIdx.x; p < n; p += blockDim.x) {
        if (p < no) hi = fmax(hi, g[p]); else lo = fmin(lo, g[p]);
    }
    double* smax = sred;
    double* smin = sred + blockDim.x;
    smax[threadIdx.x] = hi; smin[threadIdx.x] = lo;
    __syncthreads();
    for (int s = blockDim.x >> 1; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            smax[threadIdx.x] = fmax(smax[threadIdx.x], smax[threadIdx.x + s]);
            smin[threadIdx.x] = fmin(smin[threadIdx.x], smin[threadIdx.x + s]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) { stats[0] = smax[0]; stats[1] = smin[0]; }
}

// theta_ai = G_ai / (g_i - g_a) in place (G_ov is nv x no, column-major), and the largest
// |theta| as the bit pattern of a non-negative float (monotone as unsigned int).
template <typename R>
__global__ void k_pseudo_theta(R* __restrict__ gov, const double* __restrict__ g, int no, int nv,
                               unsigned int* __restrict__ maxbits)
{
    const size_t total = static_cast<size_t>(nv) * no;
    unsigned int local = 0u;
    for (size_t e = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x; e < total;
         e += static_cast<size_t>(gridDim.x) * blockDim.x) {
        const int a = static_cast<int>(e % nv), i = static_cast<int>(e / nv);
        const double t = static_cast<double>(gov[e]) / (g[i] - g[no + a]);
        gov[e] = static_cast<R>(t);
        const unsigned int bits = __float_as_uint(fabsf(static_cast<float>(t)));
        local = bits > local ? bits : local;
    }
    atomicMax(maxbits, local);
}

// Diagonal of the S-metric Gram matrix, M_pp = sum_mu B(mu,p)^2, accumulated in FP64 and written
// over the syrk diagonal. Measured on polymer: the FP32 syrk diagonal of these sums of ~3000
// squares is systematically ~1e-6 too small per column (sum over the occupied block -1.8e-3), so
// the Cholesky normalised every orbital slightly too long and the electron count drifted by
// ~4e-3 (SCC energy 1.9 mEh off). The off-diagonal elements keep their syrk values.
template <typename R>
__global__ void k_pseudo_gram_diag(const R* __restrict__ B, int n, R* __restrict__ M)
{
    extern __shared__ double sred[];
    const int p = blockIdx.x;
    const size_t off = static_cast<size_t>(p) * n;
    double acc = 0.0;
    for (int mu = threadIdx.x; mu < n; mu += blockDim.x) {
        const double x = static_cast<double>(B[off + mu]);
        acc += x * x;
    }
    sred[threadIdx.x] = acc;
    __syncthreads();
    for (int s = blockDim.x >> 1; s > 0; s >>= 1) {
        if (threadIdx.x < s) sred[threadIdx.x] += sred[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) M[static_cast<size_t>(p) + off] = static_cast<R>(sred[0]);
}

// Precision overloads of the cuBLAS/cuSOLVER calls of the pseudo step.
namespace pseudo_blas {
inline cublasStatus_t gemm(cublasHandle_t h, cublasOperation_t a, cublasOperation_t b, int m, int n, int k,
                           const float* al, const float* A, int lda, const float* B, int ldb,
                           const float* be, float* C, int ldc)
{ return cublasSgemm(h, a, b, m, n, k, al, A, lda, B, ldb, be, C, ldc); }
inline cublasStatus_t gemm(cublasHandle_t h, cublasOperation_t a, cublasOperation_t b, int m, int n, int k,
                           const double* al, const double* A, int lda, const double* B, int ldb,
                           const double* be, double* C, int ldc)
{ return cublasDgemm(h, a, b, m, n, k, al, A, lda, B, ldb, be, C, ldc); }
inline cublasStatus_t trmm(cublasHandle_t h, int n, const float* al, const float* L, const float* C, float* B)
{ return cublasStrmm(h, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, al, L, n, C, n, B, n); }
inline cublasStatus_t trmm(cublasHandle_t h, int n, const double* al, const double* L, const double* C, double* B)
{ return cublasDtrmm(h, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, al, L, n, C, n, B, n); }
inline cublasStatus_t syrk(cublasHandle_t h, int n, const float* al, const float* B, const float* be, float* M)
{ return cublasSsyrk(h, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, n, n, al, B, n, be, M, n); }
inline cublasStatus_t syrk(cublasHandle_t h, int n, const double* al, const double* B, const double* be, double* M)
{ return cublasDsyrk(h, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, n, n, al, B, n, be, M, n); }
inline cublasStatus_t trsm(cublasHandle_t h, int n, const float* al, const float* R, float* C)
{ return cublasStrsm(h, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, al, R, n, C, n); }
inline cublasStatus_t trsm(cublasHandle_t h, int n, const double* al, const double* R, double* C)
{ return cublasDtrsm(h, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, al, R, n, C, n); }
inline cusolverStatus_t potrfSize(cusolverDnHandle_t h, int n, float* M, int* lw)
{ return cusolverDnSpotrf_bufferSize(h, CUBLAS_FILL_MODE_UPPER, n, M, n, lw); }
inline cusolverStatus_t potrfSize(cusolverDnHandle_t h, int n, double* M, int* lw)
{ return cusolverDnDpotrf_bufferSize(h, CUBLAS_FILL_MODE_UPPER, n, M, n, lw); }
inline cusolverStatus_t potrf(cusolverDnHandle_t h, int n, float* M, float* w, int lw, int* info)
{ return cusolverDnSpotrf(h, CUBLAS_FILL_MODE_UPPER, n, M, n, w, lw, info); }
inline cusolverStatus_t potrf(cusolverDnHandle_t h, int n, double* M, double* w, int lw, int* info)
{ return cusolverDnDpotrf(h, CUBLAS_FILL_MODE_UPPER, n, M, n, w, lw, info); }
} // namespace pseudo_blas

/**
 * One pseudo-diagonalisation step in precision R (float or double), Claude Generated (Sep 29, 2026).
 * In:  F (n x n Fock, AO basis), L (lower Cholesky factor of S), C (previous eigenvectors, C^T S C = I,
 *      occupied block first). Scratch: Tb (n x n), Gov (nv x no), M (n x n; may alias F - F is not read
 *      after step 1), work (potrf).
 * Out: C rotated and S-orthonormalised, g = Rayleigh quotients of the rotated-from vectors.
 * Returns 1 on success, 0 when the step was rejected before C was touched (small gap or a large
 * rotation), -1 on a failure after C was modified.
 */
template <typename R>
int XtbGpuContext::Impl::pseudoRotateOrthonormalise(int n, int no, const R* F, const R* L, R* C, R* Tb,
                                                    R* Gov, R* M, CudaBuffer<R>& work, double* g, const char* tag)
{
    XtbGpuContext::Impl& I = *this;
    namespace pb = pseudo_blas;
    const int nv = n - no;
    const int b = 256;
    const R one = 1, zero = 0, mone = -1;
    cudaStream_t stream = I.stream;
    R* Cv = C + static_cast<size_t>(no) * n;
    auto mark = [&](const char* what) { I.profMark((std::string(tag) + what).c_str()); };

    // 1. T = F C, g = diag(C^T F C), G_ov = C_v^T (F C_o).
    if (pb::gemm(I.cublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, F, n, C, n, &zero, Tb, n) != CUBLAS_STATUS_SUCCESS)
        return 0;
    k_pseudo_coldot<R><<<n, b, b * sizeof(double), stream>>>(C, Tb, n, g);
    if (pb::gemm(I.cublas, CUBLAS_OP_T, CUBLAS_OP_N, nv, no, n, &one, Cv, n, Tb, n, &zero, Gov, nv) != CUBLAS_STATUS_SUCCESS)
        return 0;
    mark("F*C, diag(G), G_ov");

    // 2. Is one rotation safe? The occupied/virtual gap must be large against kT (occupations stay
    //    0/2 and the two blocks keep their order) and every rotation angle small.
    k_pseudo_gap<<<1, b, 2 * b * sizeof(double), stream>>>(g, no, n, I.dPsStats.ptr);
    I.dPsMax.zero(1, stream);
    k_pseudo_theta<R><<<1024, b, 0, stream>>>(Gov, g, no, nv, I.dPsMax.ptr);
    if (cudaGetLastError() != cudaSuccess) return 0;
    double stats[2] = {0.0, 0.0};
    unsigned int maxbits = 0u;
    I.dPsStats.download(stats, 2, stream);
    I.dPsMax.download(&maxbits, 1, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess) return 0;
    float max_theta = 0.0f;
    std::memcpy(&max_theta, &maxbits, sizeof(float));
    I.pseudo_last_max_theta = max_theta;
    const double kT = I.loop_Tele * 3.166808e-6;
    if (!(stats[1] - stats[0] > std::max(30.0 * kT, 1.0e-3)) || !(max_theta < 0.1f)) return 0;
    mark("gap/theta check");

    // 3. Rotate: keep the old C_v in Tb, then C_v -= C_o theta^T, C_o += C_v(old) theta.
    if (cudaMemcpyAsync(Tb, Cv, sizeof(R) * static_cast<size_t>(n) * nv, cudaMemcpyDeviceToDevice, stream)
            != cudaSuccess
        || pb::gemm(I.cublas, CUBLAS_OP_N, CUBLAS_OP_T, n, nv, no, &mone, C, n, Gov, nv, &one, Cv, n)
            != CUBLAS_STATUS_SUCCESS
        || pb::gemm(I.cublas, CUBLAS_OP_N, CUBLAS_OP_N, n, no, nv, &one, Tb, n, Gov, nv, &one, C, n)
            != CUBLAS_STATUS_SUCCESS)
        return -1;
    mark("rotation");

    // 4. S-orthonormalise: B = L^T C (into Tb), M = B^T B (upper, diagonal re-summed in FP64),
    //    M = R^T R, C <- C R^-1. The Cholesky order keeps the occupied span (column i mixes only
    //    with columns <= i).
    if (pb::trmm(I.cublas, n, &one, L, C, Tb) != CUBLAS_STATUS_SUCCESS) return -1;
    mark("orth trmm L^T C");
    if (pb::syrk(I.cublas, n, &one, Tb, &zero, M) != CUBLAS_STATUS_SUCCESS) return -1;
    k_pseudo_gram_diag<R><<<n, b, b * sizeof(double), stream>>>(Tb, n, M);
    if (cudaGetLastError() != cudaSuccess) return -1;
    mark("orth syrk + FP64 diagonal");
    int lw = 0;
    if (pb::potrfSize(I.cusolver, n, M, &lw) != CUSOLVER_STATUS_SUCCESS) return -1;
    try { work.ensure(std::max(1, lw)); I.dInfo.ensure(1); } catch (...) { return -1; }
    if (pb::potrf(I.cusolver, n, M, work.ptr, lw, I.dInfo.ptr) != CUSOLVER_STATUS_SUCCESS) return -1;
    int info = 0;
    I.dInfo.download(&info, 1, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess || info != 0) return -1;
    mark("orth potrf");
    if (pb::trsm(I.cublas, n, &one, M, C) != CUBLAS_STATUS_SUCCESS) return -1;
    mark("orth trsm C R^-1");
    return 1;
}

bool XtbGpuContext::Impl::pseudoDiagStep(int n)
{
    const int no = loop_nocc_pairs, nv = n - no;
    if (no <= 0 || nv <= 0) return false;
    const size_t nn = static_cast<size_t>(n) * n;
    try {
        dCf.ensure(static_cast<int>(nn));
        dLf.ensure(static_cast<int>(nn));
        dPsT.ensure(static_cast<int>(nn));
        dPsGov.ensure(static_cast<int>(static_cast<size_t>(nv) * no));
        dPsG.ensure(n);
        dPsStats.ensure(2);
        dPsMax.ensure(1);
    } catch (...) { return false; }
    const int b = 256;
    const int gnn = static_cast<int>((nn + b - 1) / b);
    k_d2f<<<gnn, b, 0, stream>>>(dC.ptr, dCf.ptr, nn);   // F -> FP32
    k_d2f<<<gnn, b, 0, stream>>>(dL.ptr, dLf.ptr, nn);   // L -> FP32
    if (cudaGetLastError() != cudaSuccess) return false;
    // M aliases the FP32 F (dCf): F is only read in step 1.
    const int rc = pseudoRotateOrthonormalise<float>(n, no, dCf.ptr, dLf.ptr, dCprev32.ptr, dPsT.ptr,
                                                     dPsGov.ptr, dCf.ptr, dPsWork, dPsG.ptr, "eig FP32 pseudo: ");
    if (rc < 0) cprev_valid = false;
    if (rc != 1) return false;
    // Hand the result to the density: C in FP64, orbital energies = the Rayleigh quotients.
    k_f2d<<<gnn, b, 0, stream>>>(dCprev32.ptr, dC.ptr, nn);
    if (cudaGetLastError() != cudaSuccess
        || cudaMemcpyAsync(dEps.ptr, dPsG.ptr, sizeof(double) * n, cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return false;
    c_exact_eigenvectors = false;
    return true;
}

int XtbGpuContext::Impl::pseudoDiagStep64(int n)
{
    const int no = loop_nocc_pairs, nv = n - no;
    if (no <= 0 || nv <= 0) return 0;
    const size_t nn = static_cast<size_t>(n) * n;
    // The FP32 scratch is idle in an FP64 step; release it before the FP64 buffers grow.
    dCf.free(); dLf.free(); dPsT.free(); dPsGov.free(); dPsWork.free();
    try {
        dCprev64.ensure(static_cast<int>(nn));
        dPsT64.ensure(static_cast<int>(nn));
        dPsGov64.ensure(static_cast<int>(static_cast<size_t>(nv) * no));
        dPsG.ensure(n);
        dPsStats.ensure(2);
        dPsMax.ensure(1);
    } catch (...) { return 0; }
    const int b = 256;
    const int gnn = static_cast<int>((nn + b - 1) / b);
    if (!cprev64_valid) {   // start from the FP32 eigenvectors of the previous step
        k_f2d<<<gnn, b, 0, stream>>>(dCprev32.ptr, dCprev64.ptr, nn);
        if (cudaGetLastError() != cudaSuccess) return 0;
        cprev64_valid = true;
    }
    // M aliases F (dC): F is only read in step 1; C is copied back into dC afterwards.
    const int rc = pseudoRotateOrthonormalise<double>(n, no, dC.ptr, dL.ptr, dCprev64.ptr, dPsT64.ptr,
                                                      dPsGov64.ptr, dC.ptr, dPsWork64, dPsG.ptr, "eig FP64 pseudo: ");
    // rc < 0: F (dC) may already hold the Gram matrix - no silent fallback to a full solve on it.
    if (rc < 0) { cprev64_valid = false; return -1; }
    if (rc != 1) return 0;
    if (cudaMemcpyAsync(dC.ptr, dCprev64.ptr, sizeof(double) * nn, cudaMemcpyDeviceToDevice, stream) != cudaSuccess
        || cudaMemcpyAsync(dEps.ptr, dPsG.ptr, sizeof(double) * n, cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return -1;
    // The FP32 copy follows, so a later FP32 step starts from these vectors.
    if (!dCprev32.empty() && dCprev32.n >= static_cast<int>(nn)) {
        k_d2f<<<gnn, b, 0, stream>>>(dCprev64.ptr, dCprev32.ptr, nn);
        cprev_valid = cudaGetLastError() == cudaSuccess;
    }
    c_exact_eigenvectors = false;
    return 1;
}

void XtbGpuContext::setPseudoDiagonalisation(bool on, double max_dq, int max_steps, bool fp64)
{
    if (!m_impl) return;
    m_impl->pseudo_enabled = on;
    m_impl->pseudo_max_dq = max_dq;
    m_impl->pseudo_max_steps = max_steps;
    m_impl->pseudo_fp64 = on && fp64;
    if (!on) m_impl->releasePseudoDiag();
}

void XtbGpuContext::pseudoDiagonalisationCounts(int& taken, int& rejected) const
{
    taken = m_impl ? m_impl->pseudo_taken : 0;
    rejected = m_impl ? m_impl->pseudo_rejected : 0;
}

bool XtbGpuContext::eigensolveResidentFock(double* eps_out, bool fp32, int n_eig, bool download_eps)
{
    Impl& I = *m_impl;
    const int n = I.resident_n;
    const bool partial = (n_eig > 0 && n_eig < n);
    const bool closed_shell = std::fabs(I.loop_nelec - 2.0 * I.loop_nocc_pairs) < 1e-9;
    const bool pseudo_ok = I.pseudo_enabled && I.in_resident_step && !partial && closed_shell
        && I.loop_nao == n && I.last_dq < I.pseudo_max_dq;
    auto finish_pseudo = [&](const char* what) {
        ++I.pseudo_taken;
        I.profMark(what);
        if (download_eps && eps_out) I.dEps.download(eps_out, n, I.stream);
        return cudaStreamSynchronize(I.stream) == cudaSuccess;
    };
    if (pseudo_ok && fp32 && I.cprev_valid && I.pseudo_consecutive < I.pseudo_max_steps) {
        if (I.pseudoDiagStep(n)) {
            ++I.pseudo_consecutive;
            I.cprev64_valid = false;   // the FP64 copy is now one step behind
            return finish_pseudo("eig FP32 pseudo: hand-over to the density");
        }
        ++I.pseudo_rejected;
        // dCprev32 may be half-rotated; the full solve below refreshes it.
    } else if (pseudo_ok && !fp32 && I.pseudo_fp64 && (I.cprev64_valid || I.cprev_valid)
               && I.pseudo_consecutive64 < I.pseudo_max_steps) {
        const int rc = I.pseudoDiagStep64(n);
        if (rc < 0) {
            I.last_error = "FP64 pseudo-diagonalisation step failed after the Fock matrix was consumed";
            return false;
        }
        if (rc == 1) {
            ++I.pseudo_consecutive64;
            return finish_pseudo("eig FP64 pseudo: hand-over to the density");
        }
        ++I.pseudo_rejected;
    }
    const bool ok = eigensolveResidentFockFull(eps_out, fp32, n_eig, download_eps);
    if (ok) I.c_exact_eigenvectors = true;
    if (ok && I.pseudo_enabled && !partial) {
        const size_t nn = static_cast<size_t>(n) * n;
        const int b = 256;
        bool stored = false;
        try {
            I.dCprev32.ensure(static_cast<int>(nn));
            k_d2f<<<static_cast<int>((nn + b - 1) / b), b, 0, I.stream>>>(I.dC.ptr, I.dCprev32.ptr, nn);
            stored = cudaGetLastError() == cudaSuccess;
        } catch (...) { stored = false; }
        I.cprev_valid = stored;
        I.pseudo_consecutive = 0;
        I.cprev64_valid = false;
        if (I.pseudo_fp64 && !fp32) {   // an FP64 solve seeds the FP64 pseudo steps directly
            try {
                I.dCprev64.ensure(static_cast<int>(nn));
                I.cprev64_valid = cudaMemcpyAsync(I.dCprev64.ptr, I.dC.ptr, sizeof(double) * nn,
                                                  cudaMemcpyDeviceToDevice, I.stream) == cudaSuccess;
            } catch (...) { I.cprev64_valid = false; }
            I.pseudo_consecutive64 = 0;
        }
    }
    return ok;
}

bool XtbGpuContext::eigensolveResidentFockFull(double* eps_out, bool fp32, int n_eig,
                                               bool download_eps)
{
    const int n = m_impl->resident_n;
    cudaStream_t stream = m_impl->stream;
    const double one = 1.0;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);

    // AP1 (Claude Generated): partial diagonalisation. The density only needs the
    // occupied(+buffer) columns, so when 0 < n_eig < n we compute the lowest n_eig
    // eigenpairs (cusolverDnDsyevdx / Ssyevdx, range il=1..neig). The reduction
    // (full n×n) and the FP32/FP64 split are unchanged; only the eigensolver and the
    // back-transform (neig columns) shrink. eps_out[neig..n) gets an ascending
    // sentinel so the host occupation/Fermi logic is bit-faithful (those orbitals
    // carry occ≈0). neig==n is the full spectrum (cusolverDnDsyevd).
    const bool partial = (n_eig > 0 && n_eig < n);
    const int neig = partial ? n_eig : n;

    // Claude Generated (Sep 2026, multi-GPU step 3): FP64 iterations run the whole generalized
    // solve on several GPUs (sygst + syevd + trsm; L is distributed once per geometry). This device
    // only sends F and receives C, and allocates no cuSOLVER workspace. FP32 iterations distribute
    // only the eigensolve (distSolve in the FP32 branch below): measured on polymer_2x (n = 15444,
    // 4x A4500) the distributed reduction wins in FP64 (sygst 3.7 s vs two single-GPU trsm 21 s)
    // but not in FP32 (0.56 vs 0.57 s), and the cuBLASMp FP32 back-transform is slower than the
    // single-GPU trsm (1.74 vs 0.38 s).
    // Claude Generated (Sep 2026, G2-7): name the two cases in which the multi-GPU solver is set up
    // but this FP64 solve still runs (partly) on this device.
    if (m_impl->dist_devices.size() > 1 && n >= m_impl->dist_min_nao && m_impl->dist_degrade.empty()) {
        if (partial)
            m_impl->dist_degrade = "-scf_gpu_partial_diag runs the eigensolve on one device (the split is off)";
        else if (!fp32 && m_impl->distReady(n, false) && !m_impl->dist->supportsGeneralized())
            m_impl->dist_degrade = "no cuBLASMp: FP64 reduction and back-transform stay on one device, only syevd is split";
    }
    if (!partial && !fp32 && m_impl->distReady(n, false) && m_impl->dist->supportsGeneralized()) {
        // No single-GPU workspaces or FP32 copies while the FP64 solve is distributed.
        m_impl->dCf.free(); m_impl->dLf.free(); m_impl->dWorkf.free(); m_impl->lwork_f32 = 0;
        if (!m_impl->dWork.empty()) { m_impl->dWork.free(); m_impl->lwork = 0; }
        const int rc = m_impl->distSolveGeneralized(n, m_impl->dC.ptr, m_impl->dL.ptr, m_impl->dEps.ptr, false);
        if (rc < 0) return false;
        if (rc == 1) {
            m_impl->profMark("eig FP64: multi-GPU sygst+syevd+trsm");
            if (download_eps && eps_out) m_impl->dEps.download(eps_out, n, stream);
            return cudaStreamSynchronize(stream) == cudaSuccess;
        }
        // rc == 0: F is intact, continue with the single-GPU path below.
    }

    if (fp32) {
        // Mixed precision: reduce + diagonalise + back-transform in FP32, then
        // convert the eigenvectors/values back to FP64 so the resident density
        // (FP64) is unchanged. ~5–10× faster than FP64 on consumer GPUs.
        try {
            // Claude Generated (Sep 2026): the FP64 syevd workspace (~3.6 nao^2 doubles) is idle
            // while the SCF runs in FP32; release it so both precisions never coexist.
            if (!m_impl->dWork.empty()) { m_impl->dWork.free(); m_impl->lwork = 0; }
            if (m_impl->dCf.n   < static_cast<int>(nn)) m_impl->dCf.alloc(static_cast<int>(nn));
            if (m_impl->dLf.n   < static_cast<int>(nn)) m_impl->dLf.alloc(static_cast<int>(nn));
            if (m_impl->dEpsf.n < n)                    m_impl->dEpsf.alloc(n);
        } catch (...) { return false; }
        const int b = 256;
        const int gnn = static_cast<int>((nn + b - 1) / b);
        k_d2f<<<gnn, b, 0, stream>>>(m_impl->dC.ptr, m_impl->dCf.ptr, nn);  // F → FP32
        k_d2f<<<gnn, b, 0, stream>>>(m_impl->dL.ptr, m_impl->dLf.ptr, nn);  // L → FP32
        if (cudaGetLastError() != cudaSuccess) return false;
        const float onef = 1.0f;
        if (cublasStrsm(m_impl->cublas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                        CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &onef,
                        m_impl->dLf.ptr, n, m_impl->dCf.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
        if (cublasStrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                        CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &onef,
                        m_impl->dLf.ptr, n, m_impl->dCf.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
        m_impl->profMark("eig FP32: copy + reduce");
        bool back_done = false;   // eigenvectors already back-transformed on the devices
        const int dist_rc = partial ? 0 : m_impl->distSolve(n, m_impl->dCf.ptr, m_impl->dEpsf.ptr, true,
                                                            m_impl->dLf.ptr, &back_done);
        if (dist_rc < 0) return false;
        int lwork = 0;
        if (dist_rc == 1) {
            m_impl->profMark("eig FP32: syevd (multi-GPU)");
        } else if (partial) {
            const float vl = 0.0f, vu = 0.0f; int h_meig = 0;
            if (cusolverDnSsyevdx_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                             CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_LOWER, n,
                                             m_impl->dCf.ptr, n, vl, vu, 1, neig, &h_meig,
                                             m_impl->dEpsf.ptr, &lwork) != CUSOLVER_STATUS_SUCCESS)
                return false;
            if (m_impl->lwork_f32 < lwork) {
                try { m_impl->dWorkf.alloc(lwork > 0 ? lwork : 1); } catch (...) { return false; }
                m_impl->lwork_f32 = lwork;
            }
            if (cusolverDnSsyevdx(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                  CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_LOWER, n,
                                  m_impl->dCf.ptr, n, vl, vu, 1, neig, &h_meig,
                                  m_impl->dEpsf.ptr, m_impl->dWorkf.ptr, m_impl->lwork_f32,
                                  m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
                return false;
            m_impl->profMark("eig FP32: syevd");
        } else {
            if (cusolverDnSsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                            CUBLAS_FILL_MODE_LOWER, n, m_impl->dCf.ptr, n,
                                            m_impl->dEpsf.ptr, &lwork) != CUSOLVER_STATUS_SUCCESS)
                return false;
            if (m_impl->lwork_f32 < lwork) {
                try { m_impl->dWorkf.alloc(lwork > 0 ? lwork : 1); } catch (...) { return false; }
                m_impl->lwork_f32 = lwork;
            }
            if (cusolverDnSsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                 n, m_impl->dCf.ptr, n, m_impl->dEpsf.ptr, m_impl->dWorkf.ptr,
                                 m_impl->lwork_f32, m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
                return false;
            m_impl->profMark("eig FP32: syevd");
        }
        // Back-transform only the neig computed eigenvectors (columns 0..neig-1).
        if (!back_done
            && cublasStrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                           CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, neig, &onef,
                           m_impl->dLf.ptr, n, m_impl->dCf.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
        // Convert the neig columns (column-major: first neig·n contiguous) + neig eps.
        const size_t conv = static_cast<size_t>(neig) * static_cast<size_t>(n);
        k_f2d<<<static_cast<int>((conv + b - 1) / b), b, 0, stream>>>(m_impl->dCf.ptr, m_impl->dC.ptr, conv);
        k_f2d<<<(neig + b - 1) / b, b, 0, stream>>>(m_impl->dEpsf.ptr, m_impl->dEps.ptr, neig);
        if (cudaGetLastError() != cudaSuccess) return false;
        m_impl->profMark("eig FP32: back-transform + to FP64");

        int info = 0;
        if (dist_rc != 1
            && cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
                               cudaMemcpyDeviceToHost, stream) != cudaSuccess)
            return false;
        if (download_eps && eps_out) m_impl->dEps.download(eps_out, neig, stream);
        if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
        if (download_eps && eps_out && partial) fillEpsSentinel(eps_out, neig, n);
        return info == 0;
    }

    // dC holds F → Ã = L⁻¹·F·L⁻ᵀ → C̃ → C.
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    m_impl->profMark("eig FP64: reduce");
    const int dist_rc = partial ? 0 : m_impl->distSolve(n, m_impl->dC.ptr, m_impl->dEps.ptr, false);
    if (dist_rc < 0) return false;
    if (dist_rc == 1) {
        m_impl->profMark("eig FP64: syevd (multi-GPU)");
    } else if (partial) {
        const double vl = 0.0, vu = 0.0; int h_meig = 0; int lwork_dx = 0;
        if (cusolverDnDsyevdx_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                         CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_LOWER, n,
                                         m_impl->dC.ptr, n, vl, vu, 1, neig, &h_meig,
                                         m_impl->dEps.ptr, &lwork_dx) != CUSOLVER_STATUS_SUCCESS)
            return false;
        if (m_impl->lwork < lwork_dx) {   // only grow; dWork stays ≥ both syevd and syevdx
            try { m_impl->dWork.alloc(lwork_dx > 0 ? lwork_dx : 1); } catch (...) { return false; }
            m_impl->lwork = lwork_dx;
        }
        if (cusolverDnDsyevdx(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                              CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_LOWER, n,
                              m_impl->dC.ptr, n, vl, vu, 1, neig, &h_meig,
                              m_impl->dEps.ptr, m_impl->dWork.ptr, m_impl->lwork,
                              m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->profMark("eig FP64: syevd");
    } else {
        // Claude Generated (Sep 2026): FP64 workspace on demand; drop the FP32 copies first.
        if (m_impl->lwork <= 0 || m_impl->dWork.empty()) {
            m_impl->dCf.free(); m_impl->dLf.free(); m_impl->dWorkf.free(); m_impl->lwork_f32 = 0;
            int lwork = 0;
            if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                            CUBLAS_FILL_MODE_LOWER, n, m_impl->dC.ptr, n,
                                            m_impl->dEps.ptr, &lwork) != CUSOLVER_STATUS_SUCCESS)
                return false;
            try { m_impl->dWork.alloc(lwork > 0 ? lwork : 1); } catch (...) { return false; }
            m_impl->lwork = lwork;
        }
        if (cusolverDnDsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                             n, m_impl->dC.ptr, n, m_impl->dEps.ptr, m_impl->dWork.ptr,
                             m_impl->lwork, m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->profMark("eig FP64: syevd");
    }
    // Back-transform only the neig computed eigenvectors (columns 0..neig-1).
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, neig, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    m_impl->profMark("eig FP64: back-transform");

    int info = 0;
    if (dist_rc != 1
        && cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
                           cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (download_eps && eps_out) m_impl->dEps.download(eps_out, neig, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess)
        return false;
    if (download_eps && eps_out && partial) fillEpsSentinel(eps_out, neig, n);
    return info == 0;
}

bool XtbGpuContext::residentSolve(const double* v_ao, int n, double* eps_out, bool fp32,
                                  int n_eig)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !v_ao || !eps_out)
        return false;
    cudaStream_t stream = m_impl->stream;

    // Upload the AO potential (the only matrix-sized input is already resident).
    m_impl->dVao.upload(v_ao, n, stream);

    // F = H0 − ½·S·(v_ao⊕v_ao), built straight into the eigenvector buffer dC.
    if (!buildFockIntoC(n, /*multipole=*/false))
        return false;

    return eigensolveResidentFock(eps_out, fp32, n_eig);
}

bool XtbGpuContext::residentDensity(const double* occ, int ncol, int n,
                                    double* pop_ao_out, double* band_out)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !pop_ao_out || !band_out)
        return false;
    cudaStream_t stream = m_impl->stream;
    const double one = 1.0, zero = 0.0;

    try { m_impl->dP.ensure(static_cast<int>(static_cast<size_t>(n) * n)); } catch (...) { return false; }
    m_impl->p_dense_valid = true;
    if (ncol > 0) {
        m_impl->dOcc.upload(occ, ncol, stream);
        // Cw(:,k) = occ[k]·C(:,k)  for k < ncol  (n x ncol; Claude Generated Sep 2026: sized to
        // the occupied columns instead of n x n).
        try { m_impl->dCw.ensure(n * ncol); } catch (...) { return false; }
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (ncol + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, stream>>>(m_impl->dCw.ptr, m_impl->dC.ptr,
                                                 m_impl->dOcc.ptr, n, ncol);
        if (cudaGetLastError() != cudaSuccess)
            return false;
        // P = Cw[:, :ncol] · C[:, :ncol]ᵀ  ⇒  P_ij = Σ_k occ_k C_ik C_jk.
        if (cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, ncol,
                        &one, m_impl->dCw.ptr, n, m_impl->dC.ptr, n,
                        &zero, m_impl->dP.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
    } else {
        // No occupied orbitals — zero density.
        if (cudaMemsetAsync(m_impl->dP.ptr, 0, sizeof(double) * static_cast<size_t>(n) * n,
                            stream) != cudaSuccess)
            return false;
    }

    // pop_ao(μ) = Σ_ν P(μ,ν)·S(μ,ν); band energy = Σ_μν P_μν·H0_μν.
    if (!populationsAndBand(n, band_out))
        return false;

    m_impl->dPop.download(pop_ao_out, n, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

// Stage 6: density + Mulliken-AO from the RESIDENT occupations (dOcc, set by
// k_occupations) — no occ upload, no pop_ao download (dPop stays resident for the
// device q_sh/q_at reductions). band_out = Σ P⊙H0 (host scalar). Claude Generated.
// Claude Generated (Sep 2026, multi-GPU): pattern density over several GPUs.
//
// P(r,c) = sum_k Cw(r,k) C(c,k) is a plain sum over the occupied columns, so each device can
// evaluate the full pattern over its own contiguous slice of k and the partials are added on the
// primary device. Per SCF step only the column slices of C (n x kn each) and the partial pattern
// arrays travel; the pattern indices are geometry-constant and uploaded once.
//
// Returns false when the split is not set up or fails, and leaves the target buffer untouched
// in that case, so the caller falls back to the single-device kernel.
//
// Point 6 (Claude Generated, Sep 2026): for_weighted_density retargets the same column-split
// SDDMM at the gradient's W = C_occ*diag(2eps)*C_occ^T instead of the SCF loop's P. Both are
// the identical sum-over-occupied-columns pattern (k_density_sp doesn't know or care what the
// weight vector means), so nothing below needs duplicating - only which CudaBuffer the result
// lands in changes. The caller (computeGradient) uploads its own weight vector (2*eps) into
// the SAME dOcc buffer the SCF loop uses before calling this, so dOcc itself needs no
// generalisation either.
bool XtbGpuContext::densityPatternDistributed(int n, int ncol, bool for_weighted_density)
{
    Impl& I = *m_impl;
    // Claude Generated (Sep 2026, G2-6): remember WHY the split did not run, so the status names
    // the real cause instead of always "nao below min".
    if (!I.dens_devices.empty() && !I.dens_failed) {
        if (!I.sparse || I.sp_nnz <= 0)
            I.dens_skip_reason = "dense integral storage; the split needs screened storage (-gpu_sparse_integrals on)";
        else if (n < I.dens_min_nao)
            I.dens_skip_reason = "nao " + std::to_string(n) + " below gpu_density_min_nao = " + std::to_string(I.dens_min_nao);
    }
    if (I.dens_failed || I.dens_devices.empty() || !I.sparse || I.sp_nnz <= 0 || ncol <= 0
        || n < I.dens_min_nao)
        return false;   // below ~3000 basis functions the transfers cost more than the split saves

    const int nhelp = static_cast<int>(I.dens_devices.size());
    const int nworker = nhelp + 1;               // the helpers plus this device
    if (ncol < 16 * nworker) return false;       // too few columns to be worth splitting
    CudaBuffer<double>& target = for_weighted_density ? I.dSpW : I.dSpP;
    // The pattern density buffer is normally allocated by the single-device branch, which this
    // function replaces.
    try { target.ensure(I.sp_nnz); I.dCw.ensure(n * ncol); } catch (...) { return false; }

    // Create the helpers (streams, events, pattern buffers) on first use.
    if (I.dens_helpers.empty()) {
        for (int d : I.dens_devices) {
            Impl::DensityHelper h;
            h.device = d;
            cudaSetDevice(d);
            int can = 0;
            if (cudaDeviceCanAccessPeer(&can, d, I.device) == cudaSuccess && can)
                cudaDeviceEnablePeerAccess(I.device, 0);
            cudaGetLastError();   // "already enabled" is expected and must not linger
            if (cudaStreamCreate(&h.stream) != cudaSuccess
                || cudaEventCreateWithFlags(&h.done, cudaEventDisableTiming) != cudaSuccess) {
                cudaSetDevice(I.device);
                I.dens_failed = true;
                I.dens_error = "stream/event creation failed";
                I.releaseDensityHelpers();
                return false;
            }
            I.dens_helpers.push_back(h);
        }
        cudaSetDevice(I.device);
        for (int d : I.dens_devices) {
            int can = 0;
            if (cudaDeviceCanAccessPeer(&can, I.device, d) == cudaSuccess && can)
                cudaDeviceEnablePeerAccess(d, 0);
            cudaGetLastError();
        }
    }

    // Column split: contiguous slices, the primary device takes the first one.
    const int per = (ncol + nworker - 1) / nworker;
    const int own_cols = std::min(per, ncol);

    auto fail = [&](const char* where, cudaError_t err = cudaGetLastError()) {
        cudaSetDevice(I.device);
        I.dens_failed = true;
        I.dens_error = std::string(where) + ": " + cudaGetErrorString(err);
        I.releaseDensityHelpers();
        return false;
    };

    // Launch the helpers first so they overlap with this device's own slice.
    for (int i = 0; i < nhelp; ++i) {
        Impl::DensityHelper& h = I.dens_helpers[i];
        const int k0 = std::min(ncol, (i + 1) * per);
        const int kn = std::min(per, ncol - k0);
        if (kn <= 0) continue;
        cudaSetDevice(h.device);
        // Buffers: column slice (grown to the largest slice seen) and the pattern.
        if (h.cap_cols < kn) {
            if (h.dC) cudaFree(h.dC);
            if (h.dCw) cudaFree(h.dCw);
            if (h.dOcc) cudaFree(h.dOcc);
            const size_t bytes = sizeof(double) * static_cast<size_t>(n) * kn;
            if (cudaMalloc(reinterpret_cast<void**>(&h.dC), bytes) != cudaSuccess
                || cudaMalloc(reinterpret_cast<void**>(&h.dCw), bytes) != cudaSuccess
                || cudaMalloc(reinterpret_cast<void**>(&h.dOcc), sizeof(double) * kn) != cudaSuccess)
                return fail("alloc column slice");
            h.cap_cols = kn;
        }
        if (h.cap_nnz < I.sp_nnz) {
            if (h.dPsp) cudaFree(h.dPsp);
            if (h.dRow) cudaFree(h.dRow);
            if (h.dCol) cudaFree(h.dCol);
            if (cudaMalloc(reinterpret_cast<void**>(&h.dPsp), sizeof(double) * I.sp_nnz) != cudaSuccess
                || cudaMalloc(reinterpret_cast<void**>(&h.dRow), sizeof(int) * I.sp_nnz) != cudaSuccess
                || cudaMalloc(reinterpret_cast<void**>(&h.dCol), sizeof(int) * I.sp_nnz) != cudaSuccess)
                return fail("alloc pattern");
            h.cap_nnz = I.sp_nnz;
            h.pattern_gen = -1;
            cudaSetDevice(I.device);
            if (h.stage) cudaFree(h.stage);
            if (cudaMalloc(reinterpret_cast<void**>(&h.stage), sizeof(double) * I.sp_nnz) != cudaSuccess)
                return fail("alloc staging");
            cudaSetDevice(h.device);
        }
        if (h.pattern_gen != I.pattern_generation) {
            if (cudaMemcpyPeerAsync(h.dRow, h.device, I.dSpRow.ptr, I.device,
                                    sizeof(int) * I.sp_nnz, h.stream) != cudaSuccess
                || cudaMemcpyPeerAsync(h.dCol, h.device, I.dSpCol.ptr, I.device,
                                       sizeof(int) * I.sp_nnz, h.stream) != cudaSuccess)
                return fail("copy pattern");
            h.pattern_gen = I.pattern_generation;
        }
        // Column slice of C and the matching occupations.
        if (cudaMemcpyPeerAsync(h.dC, h.device, I.dC.ptr + static_cast<size_t>(k0) * n, I.device,
                                sizeof(double) * static_cast<size_t>(n) * kn, h.stream) != cudaSuccess
            || cudaMemcpyPeerAsync(h.dOcc, h.device, I.dOcc.ptr + k0, I.device,
                                   sizeof(double) * kn, h.stream) != cudaSuccess)
            return fail("copy column slice");
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (kn + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, h.stream>>>(h.dCw, h.dC, h.dOcc, n, kn);
        const int bs = 256;
        k_density_sp<<<(I.sp_nnz + bs - 1) / bs, bs, 0, h.stream>>>(
            I.sp_nnz, h.dRow, h.dCol, h.dCw, h.dC, n, kn, h.dPsp);
        if (const cudaError_t e = cudaGetLastError(); e != cudaSuccess) return fail("helper kernels", e);
        if (cudaMemcpyPeerAsync(h.stage, I.device, h.dPsp, h.device,
                                sizeof(double) * I.sp_nnz, h.stream) != cudaSuccess
            || cudaEventRecord(h.done, h.stream) != cudaSuccess)
            return fail("copy partial back");
    }

    // This device's own slice, into the target pattern buffer (dSpP or dSpW).
    cudaSetDevice(I.device);
    {
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (own_cols + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, I.stream>>>(I.dCw.ptr, I.dC.ptr, I.dOcc.ptr, n, own_cols);
        const int bs = 256;
        k_density_sp<<<(I.sp_nnz + bs - 1) / bs, bs, 0, I.stream>>>(
            I.sp_nnz, I.dSpRow.ptr, I.dSpCol.ptr, I.dCw.ptr, I.dC.ptr, n, own_cols, target.ptr);
        if (const cudaError_t e = cudaGetLastError(); e != cudaSuccess) return fail("own slice", e);
    }

    // Add the partials once each helper has written its staging buffer.
    const double one = 1.0;
    for (int i = 0; i < nhelp; ++i) {
        Impl::DensityHelper& h = I.dens_helpers[i];
        const int k0 = std::min(ncol, (i + 1) * per);
        if (std::min(per, ncol - k0) <= 0) continue;
        if (cudaStreamWaitEvent(I.stream, h.done, 0) != cudaSuccess) return fail("wait event");
        if (cublasDaxpy(I.cublas, I.sp_nnz, &one, h.stage, 1, target.ptr, 1) != CUBLAS_STATUS_SUCCESS)
            return fail("reduce partials");
    }
    ++I.dens_steps;
    return true;
}

bool XtbGpuContext::residentDensityResident(int n, int ncol, double* band_out)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !band_out) return false;
    Impl& I = *m_impl;
    cudaStream_t stream = I.stream;
    const double one = 1.0, zero = 0.0;
    I.last_ncol = ncol;
    if (ncol > 0) {
        try { I.dCw.ensure(n * ncol); } catch (...) { return false; }   // n x ncol
    }
    // Claude Generated (Sep 2026): pattern density split over several GPUs (see
    // densityPatternDistributed); falls through to the single-device kernels below when it is
    // not configured or fails.
    const bool split_done = ncol > 0 && densityPatternDistributed(n, ncol);
    if (ncol > 0 && !split_done) {
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (ncol + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, stream>>>(I.dCw.ptr, I.dC.ptr, I.dOcc.ptr, n, ncol);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    if (split_done) {
        I.p_dense_valid = false;
    } else if (I.sparse) {
        // Claude Generated (Sep 2026): pattern-only density (see k_density_sp). A row-major
        // variant (contiguous rows, occ applied in the kernel) was measured 8x slower on polymer
        // (4.69 vs 0.60 s over 12 SCF steps) and dropped.
        try { I.dSpP.ensure(I.sp_nnz); } catch (...) { return false; }
        if (ncol > 0) {
            const int bs = 256;
            k_density_sp<<<(I.sp_nnz + bs - 1) / bs, bs, 0, stream>>>(
                I.sp_nnz, I.dSpRow.ptr, I.dSpCol.ptr, I.dCw.ptr, I.dC.ptr, n, ncol, I.dSpP.ptr);
            if (cudaGetLastError() != cudaSuccess) return false;
        } else {
            I.dSpP.zero(I.sp_nnz, stream);
        }
        I.p_dense_valid = false;
    } else if (ncol > 0) {
        if (cublasDgemm(I.cublas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, ncol,
                        &one, I.dCw.ptr, n, I.dC.ptr, n,
                        &zero, I.dP.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
        I.p_dense_valid = true;
    } else {
        if (cudaMemsetAsync(I.dP.ptr, 0, sizeof(double) * static_cast<size_t>(n) * n, stream) != cudaSuccess)
            return false;
        I.p_dense_valid = true;
    }
    return populationsAndBand(n, band_out);
}

// Claude Generated (Sep 2026): rebuild the dense density from the resident eigenvectors and
// occupations of the last step when the loop only kept it on the screened pattern.
bool XtbGpuContext::ensureDenseDensity(int n)
{
    Impl& I = *m_impl;
    if (!I.sparse || I.p_dense_valid) return true;
    cudaStream_t stream = I.stream;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    try { I.dP.ensure(static_cast<int>(nn)); } catch (...) { return false; }
    const int ncol = I.last_ncol;
    if (ncol > 0) {
        try { I.dCw.ensure(n * ncol); } catch (...) { return false; }
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (ncol + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, stream>>>(I.dCw.ptr, I.dC.ptr, I.dOcc.ptr, n, ncol);
        if (cudaGetLastError() != cudaSuccess) return false;
        const double one = 1.0, zero = 0.0;
        if (cublasDgemm(I.cublas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, ncol,
                        &one, I.dCw.ptr, n, I.dC.ptr, n, &zero, I.dP.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
    } else if (cudaMemsetAsync(I.dP.ptr, 0, sizeof(double) * nn, stream) != cudaSuccess) {
        return false;
    }
    I.p_dense_valid = true;
    return true;
}

bool XtbGpuContext::residentFinalize(double* P_colmajor, double* C_colmajor, int n)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !P_colmajor || !C_colmajor)
        return false;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    m_impl->releaseEigenWorkspaces();   // SCF finished; see releaseEigenWorkspaces
    if (!ensureDenseDensity(n)) return false;
    m_impl->dP.download(P_colmajor, static_cast<int>(nn), stream);
    m_impl->dC.download(C_colmajor, static_cast<int>(nn), stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Device-resident GFN2 multipole (Stage 2b). Layered on residentBegin.
 * ====================================================================== */

bool XtbGpuContext::residentBeginMultipole(const double* dp_int3, const double* qp_int6,
                                           const int* ao2at, int n, int nat)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || nat <= 0
        || !dp_int3 || !qp_int6 || !ao2at)
        return false;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    try {
        // Geometry-constant multipole integrals (3·nn / 6·nn) + AO→atom map.
        m_impl->dDpInt.upload(dp_int3, static_cast<int>(3 * nn), m_impl->stream);
        m_impl->dQpInt.upload(qp_int6, static_cast<int>(6 * nn), m_impl->stream);
        m_impl->dAo2at.upload(ao2at, n, m_impl->stream);
        // Per-iteration multipole potentials / atomic moments.
        m_impl->dVdp.ensure(3 * nat);
        m_impl->dVqp.ensure(6 * nat);
        m_impl->dDpAt.ensure(3 * nat);
        m_impl->dQpAt.ensure(6 * nat);
    } catch (...) {
        return false;
    }
    m_impl->resident_nat = nat;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::residentSolveMultipole(const double* v_ao, const double* v_dp,
                                           const double* v_qp, int n, double* eps_out,
                                           bool fp32, int n_eig)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || m_impl->resident_nat <= 0
        || !v_ao || !v_dp || !v_qp || !eps_out)
        return false;
    cudaStream_t stream = m_impl->stream;
    const int nat = m_impl->resident_nat;

    // Upload the iteration's potentials (isotropic AO + anisotropic multipole).
    m_impl->dVao.upload(v_ao, n, stream);
    m_impl->dVdp.upload(v_dp, 3 * nat, stream);
    m_impl->dVqp.upload(v_qp, 6 * nat, stream);

    // F = H0 − ½·S·(v_ao⊕v_ao) into dC, then add the GFN2 multipole contribution.
    if (!buildFockIntoC(n, /*multipole=*/true))
        return false;

    return eigensolveResidentFock(eps_out, fp32, n_eig);
}

bool XtbGpuContext::residentMultipoleMoments(double* dp_at3, double* qp_at6, int n, int nat)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || nat != m_impl->resident_nat
        || !dp_at3 || !qp_at6)
        return false;
    cudaStream_t stream = m_impl->stream;

    if (!multipoleMomentsResident(n, nat))
        return false;

    m_impl->dDpAt.download(dp_at3, 3 * nat, stream);
    m_impl->dQpAt.download(qp_at6, 6 * nat, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Device-side integral build (Stage 3).
 * ====================================================================== */

size_t XtbGpuContext::estimateResidentBytes(int nat, int nsh, int nao, bool is_gfn2) const
{
    return estimateStorageBytes(nat, nsh, nao, is_gfn2, 0.0);
}

size_t XtbGpuContext::estimateStorageBytes(int nat, int nsh, int nao, bool is_gfn2, double nnz) const
{
    if (!ok() || nao <= 0) return 0;
    const double nn = static_cast<double>(nao) * nao;
    // cuSOLVER workspace sizes for this n (queries only; the array arguments are not read).
    int lw_syevd = 0, lw_ssyevd = 0, lw_potrf = 0;
    cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                nao, nullptr, nao, nullptr, &lw_syevd);
    cusolverDnSsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                nao, nullptr, nao, nullptr, &lw_ssyevd);
    cusolverDnDpotrf_bufferSize(m_impl->cusolver, CUBLAS_FILL_MODE_LOWER, nao, nullptr, nao, &lw_potrf);

    double doubles = 0.0;
    // Claude Generated (Sep 2026): peak of the resident SCF (the gradient frees the eigensolver
    // workspaces before it allocates W, so it stays below this).
    if (nnz > 0.0) {
        doubles += nn;                                    // L (dense Cholesky factor)
        doubles += 3.0 * nnz + 1.5 * nnz;                 // S, H0, scratch + row/col/perm indices
    } else {
        doubles += 3.0 * nn;                              // S, H0, L
    }
    doubles += 3.0 * nn;                                  // C (holds F), P, Cw (<= n x n)
    // Only one precision's eigensolver workspace exists at a time.
    const double fp64_work = static_cast<double>(lw_syevd);
    const double fp32_work = nn + 0.5 * lw_ssyevd;        // F + L copies (2 x nn/2) + ssyevd work
    doubles += std::max(fp64_work, fp32_work) + lw_potrf;
    doubles += static_cast<double>(nsh) * nsh;            // Coulomb gamma
    doubles += 2.0 * (static_cast<double>(nat) + 1.0) * (nat + 1.0);  // D4-EEQ LU + getrf work
    if (is_gfn2) {
        doubles += 9.0 * (nnz > 0.0 ? nnz : nn);          // dipole (3) + quadrupole (6) integrals
        const double amat = 18.0 * static_cast<double>(nat) * nat;
        if (amat * sizeof(double) <= 1.0e9) doubles += amat;  // else rebuilt on the fly
    }
    doubles += 64.0 * (nao + nsh + nat) + 32.0 * 1024 * 1024;  // vectors + slack
    return static_cast<size_t>(doubles * sizeof(double));
}

std::string XtbGpuContext::profileReport() const
{
    if (!m_impl || !m_impl->prof || m_impl->prof_names.empty()) return {};
    std::string out = "GPU phase profile (CURCUMA_GPU_PROFILE, stream-synchronised; indented rows are\n"
                      "sub-phases already counted in the row above; mem = device memory in use):\n";
    double total = 0.0;
    for (size_t i = 0; i < m_impl->prof_ms.size(); ++i)
        if (m_impl->prof_names[i].rfind("  ", 0) != 0) total += m_impl->prof_ms[i];
    char line[200];
    for (size_t i = 0; i < m_impl->prof_names.size(); ++i) {
        std::snprintf(line, sizeof(line), "  %-44s %11.1f ms %6.1f %%  %5d calls  mem %7.0f MiB\n",
                      m_impl->prof_names[i].c_str(), m_impl->prof_ms[i],
                      total > 0 ? 100.0 * m_impl->prof_ms[i] / total : 0.0, m_impl->prof_calls[i],
                      m_impl->prof_mem_mib[i]);
        out += line;
    }
    std::snprintf(line, sizeof(line), "  %-44s %11.1f ms\n", "sum", total);
    out += line;
    return out;
}

void XtbGpuContext::setSparseIntegrals(int mode)
{
    if (m_impl) m_impl->sparse_mode = std::max(0, std::min(2, mode));
}
bool XtbGpuContext::sparseIntegrals() const { return m_impl && m_impl->sparse; }
double XtbGpuContext::sparseFraction() const { return m_impl ? m_impl->sp_fraction : 1.0; }
double XtbGpuContext::sparseCutoffBohr() const { return m_impl ? m_impl->sp_rmax : 0.0; }

void XtbGpuContext::setDistributedEigensolver(const std::vector<int>& devices,
                                              const std::string& backend, int block, int min_nao,
                                              bool fp32, bool verify)
{
    if (!m_impl) return;
    m_impl->dist.reset();
    bindDevice();
    m_impl->dist_devices.clear();
    int count = 0;
    cudaGetDeviceCount(&count);
    for (int d : devices) {
        if (d < 0) {   // "all"
            m_impl->dist_devices.clear();
            for (int i = 0; i < count; ++i) m_impl->dist_devices.push_back(i);
            break;
        }
        if (d < count && std::find(m_impl->dist_devices.begin(), m_impl->dist_devices.end(), d)
                             == m_impl->dist_devices.end())
            m_impl->dist_devices.push_back(d);
    }
    m_impl->dist_backend = backend.empty() ? std::string("auto") : backend;
    m_impl->dist_block = block > 0 ? block : 128;
    m_impl->dist_min_nao = std::max(0, min_nao);
    m_impl->dist_fp32 = fp32;
    m_impl->dist_verify = verify;
    m_impl->dist_fp32_state = 0;
    m_impl->dist_fp64_state = 0;
    m_impl->dist_failed = false;
    m_impl->dist_solves = 0;
    m_impl->dist_status.clear();
}

void XtbGpuContext::setDensityDevices(const std::vector<int>& devices, int min_nao)
{
    if (!m_impl) return;
    bindDevice();
    m_impl->releaseDensityHelpers();
    m_impl->dens_devices.clear();
    m_impl->dens_failed = false;
    m_impl->dens_steps = 0;
    m_impl->dens_min_nao = std::max(0, min_nao);
    int count = 0;
    cudaGetDeviceCount(&count);
    for (int d : devices) {
        if (d < 0) {   // "all"
            m_impl->dens_devices.clear();
            for (int i = 0; i < count; ++i)
                if (i != m_impl->device) m_impl->dens_devices.push_back(i);
            return;
        }
        if (d >= 0 && d < count && d != m_impl->device
            && std::find(m_impl->dens_devices.begin(), m_impl->dens_devices.end(), d)
                   == m_impl->dens_devices.end())
            m_impl->dens_devices.push_back(d);
    }
}

std::string XtbGpuContext::densityDevicesStatus() const
{
    if (!m_impl || m_impl->dens_devices.empty()) return {};
    if (m_impl->dens_steps == 0 && !m_impl->dens_failed)
        return "not used (" + (m_impl->dens_skip_reason.empty()
            ? std::string("never reached: only the GFN2 resident SCF loop splits the density, GFN1 has no such path")
            : m_impl->dens_skip_reason) + ")";
    if (m_impl->dens_failed)
        return "failed (" + m_impl->dens_error + "); single-device density from here on";
    std::string list = std::to_string(m_impl->device);
    for (int d : m_impl->dens_devices) list += "," + std::to_string(d);
    return "pattern density on devices [" + list + "], " + std::to_string(m_impl->dens_steps) + " steps";
}

std::string XtbGpuContext::distributedEigensolverDegradation() const
{
    return m_impl ? m_impl->dist_degrade : std::string();
}

std::string XtbGpuContext::distributedEigensolverStatus() const
{
    if (!m_impl || m_impl->dist_devices.size() < 2) return {};
    if (!m_impl->dist_status.empty()) return m_impl->dist_status;
    if (!m_impl->dist)
        return "not used (nao below gpu_eigensolver_min_nao = " + std::to_string(m_impl->dist_min_nao) + ")";
    return std::string(m_impl->dist->name()) + " on " + std::to_string(m_impl->dist->deviceCount())
        + " GPUs, " + std::to_string(m_impl->dist_solves) + " solves"
        + (m_impl->dist_verify ? " (each verified)" : " (verification off)")
        + (m_impl->dist_potrf > 0 ? ", " + std::to_string(m_impl->dist_potrf) + " Cholesky factorisation(s) of S (checked)"
                                  : std::string())
        + (m_impl->dist_potrf_status.empty() ? std::string() : "; " + m_impl->dist_potrf_status);
}

void XtbGpuContext::setMemoryCheck(bool on)
{
    if (m_impl) m_impl->memory_check = on;
}

std::string XtbGpuContext::lastError() const
{
    return m_impl ? m_impl->last_error : std::string();
}

bool XtbGpuContext::beginBasis(const XtbGpuBasisData& b)
{
    if (!ok() || b.nat <= 0 || b.nsh <= 0 || b.nao <= 0 || b.nprim_total <= 0
        || !b.z || !b.sh2at || !b.ang_sh || !b.iao_sh || !b.nao_sh
        || !b.sh_nprim || !b.sh_prim_off || !b.prim_alpha || !b.prim_coeff
        || !b.sh_zeta || !b.selfenergy || !b.kcn || !b.shpoly)
        return false;
    m_impl->last_error.clear();

    // A new basis size invalidates every large buffer (and the memory check, which now runs
    // in computeIntegrals once the geometry - hence the screened pair count - is known).
    if (b.nao != m_impl->basis_nao)
        m_impl->releaseLarge();

    ensureStage3Constants();
    try {
        // Molecule-constant uploads (synchronous memcpy inside CudaBuffer::upload).
        m_impl->dZ.upload(b.z, b.nat, m_impl->stream);
        m_impl->dSh2at.upload(b.sh2at, b.nsh, m_impl->stream);
        m_impl->dAng.upload(b.ang_sh, b.nsh, m_impl->stream);
        m_impl->dIaoSh.upload(b.iao_sh, b.nsh, m_impl->stream);
        m_impl->dNaoSh.upload(b.nao_sh, b.nsh, m_impl->stream);
        m_impl->dShNprim.upload(b.sh_nprim, b.nsh, m_impl->stream);
        m_impl->dShPrimOff.upload(b.sh_prim_off, b.nsh, m_impl->stream);
        m_impl->dPrimAlpha.upload(b.prim_alpha, b.nprim_total, m_impl->stream);
        m_impl->dPrimCoeff.upload(b.prim_coeff, b.nprim_total, m_impl->stream);
        m_impl->dShZeta.upload(b.sh_zeta, b.nsh, m_impl->stream);
        m_impl->dSelfE0.upload(b.selfenergy, b.nsh, m_impl->stream);
        m_impl->dKcn.upload(b.kcn, b.nsh, m_impl->stream);
        m_impl->dShpoly.upload(b.shpoly, b.nsh, m_impl->stream);
        if (b.shell_hardness) m_impl->dHardness.upload(b.shell_hardness, b.nsh, m_impl->stream);
        if (b.rep_alpha) m_impl->dRepAlpha.upload(b.rep_alpha, b.nat, m_impl->stream);
        if (b.rep_zeff)  m_impl->dRepZeff.upload(b.rep_zeff, b.nat, m_impl->stream);
        // Valence flags (GFN1). For GFN2 the kernel ignores them; upload zeros so
        // the device pointer is always valid.
        if (b.valence) {
            m_impl->dValence.upload(b.valence, b.nsh, m_impl->stream);
        } else {
            std::vector<int> zeros(b.nsh, 0);
            m_impl->dValence.upload(zeros.data(), b.nsh, m_impl->stream);
        }
        // Resident integral outputs + per-geometry geometry buffer.
        m_impl->dCN.ensure(b.nat);
        m_impl->dSE.ensure(b.nsh);
        m_impl->dXyz.ensure(3 * b.nat);
        m_impl->dGamma.ensure(b.nsh * b.nsh);
        // AO→atom / AO→shell maps (used by the multipole integrals and the
        // overlap-derivative kernel; cheap, uploaded for both methods).
        if (b.ao2at) m_impl->dAo2at.upload(b.ao2at, b.nao, m_impl->stream);
        if (b.ao2sh) m_impl->dAo2sh.upload(b.ao2sh, b.nao, m_impl->stream);
        if (m_impl->dInfo.n < 1) m_impl->dInfo.alloc(1);
        // The nao²-sized integral storage (dense S/H0/L/dp/qp or the screened pair arrays)
        // is allocated by computeIntegrals, which knows the geometry.
    } catch (const std::exception& e) {
        m_impl->last_error = std::string("device allocation failed: ") + e.what();
        m_impl->releaseLarge();
        return false;
    }
    // Claude Generated (Sep 2026): per-atom screening data for the sparse integral storage.
    // amin = smallest primitive exponent on the atom, cmax = largest |contraction coefficient|;
    // together they bound every integral between two atoms (see screenRadius). The AOs of an
    // atom must form one contiguous, ascending block (the pair layout relies on it).
    m_impl->h_z.assign(b.z, b.z + b.nat);
    m_impl->h_at_amin.assign(b.nat, 1.0e300);
    m_impl->h_at_cmax.assign(b.nat, 0.0);
    for (int ish = 0; ish < b.nsh; ++ish) {
        const int at = b.sh2at[ish];
        for (int p = 0; p < b.sh_nprim[ish]; ++p) {
            const int ip = b.sh_prim_off[ish] + p;
            m_impl->h_at_amin[at] = std::min(m_impl->h_at_amin[at], b.prim_alpha[ip]);
            m_impl->h_at_cmax[at] = std::max(m_impl->h_at_cmax[at], std::fabs(b.prim_coeff[ip]));
        }
    }
    m_impl->h_at_ao0.assign(b.nat, -1);
    m_impl->h_at_nao.assign(b.nat, 0);
    m_impl->h_ao_contiguous = (b.ao2at != nullptr);
    if (b.ao2at) {
        for (int mu = 0; mu < b.nao; ++mu) {
            const int at = b.ao2at[mu];
            if (mu > 0 && at < b.ao2at[mu - 1]) { m_impl->h_ao_contiguous = false; break; }
            if (m_impl->h_at_ao0[at] < 0) m_impl->h_at_ao0[at] = mu;
            ++m_impl->h_at_nao[at];
        }
    }

    m_impl->basis_nat     = b.nat;
    m_impl->basis_nsh     = b.nsh;
    m_impl->basis_nao     = b.nao;
    m_impl->basis_is_gfn2 = b.is_gfn2;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

namespace {

// Radius (Bohr) beyond which every overlap-type integral between two atoms with smallest
// primitive exponents a1/a2 and largest contraction coefficients c1/c2 is below eps.
// A primitive product decays as c1 c2 (pi/g)^1.5 exp(-a1 a2/g R^2), g = a1 + a2; the angular
// and multipole factors add at most a polynomial in R, bounded here by (1 + R)^6, and the
// origin shift of the multipole integrals multiplies by up to (1 + |r|max)^2. The bound is
// deliberately loose: a larger radius costs memory, never accuracy. Claude Generated.
double screenRadius(double a1, double c1, double a2, double c2, double rabs, double eps)
{
    const double g = a1 + a2;
    const double mu = a1 * a2 / g;
    const double logc = std::log(c1 * c2 * std::pow(M_PI / g, 1.5) + 1.0e-300)
        + 2.0 * std::log1p(rabs);
    const double target = std::log(eps);
    auto significant = [&](double R) { return logc + 6.0 * std::log1p(R) - mu * R * R > target; };
    double lo = 0.0, hi = 8.0;
    while (significant(hi) && hi < 1.0e6) hi *= 2.0;
    for (int it = 0; it < 80; ++it) {
        const double mid = 0.5 * (lo + hi);
        (significant(mid) ? lo : hi) = mid;
    }
    return hi;
}

} // namespace

bool XtbGpuContext::buildScreenedPairs(const double* xyz_bohr, bool& use_sparse)
{
    use_sparse = false;
    Impl& I = *m_impl;
    const int nat = I.basis_nat, nao = I.basis_nao;
    I.sp_fraction = 1.0;
    I.screen_nk = 0;
    if (I.sparse_mode == 0 || !I.h_ao_contiguous || nat <= 0) return true;

    // Element-pair cutoffs (atoms of one element share amin/cmax).
    double rabs = 0.0;
    for (int a = 0; a < 3 * nat; ++a) rabs = std::max(rabs, std::fabs(xyz_bohr[a]));
    std::map<int, int> kind_of_z;
    std::vector<int> kind(nat);
    std::vector<double> k_amin, k_cmax;
    for (int a = 0; a < nat; ++a) {
        auto it = kind_of_z.find(I.h_z[a]);
        if (it == kind_of_z.end()) {
            it = kind_of_z.emplace(I.h_z[a], static_cast<int>(k_amin.size())).first;
            k_amin.push_back(I.h_at_amin[a]);
            k_cmax.push_back(I.h_at_cmax[a]);
        }
        kind[a] = it->second;
    }
    const int nk = static_cast<int>(k_amin.size());
    std::vector<double> cut2(static_cast<size_t>(nk) * nk);
    I.sp_rmax = 0.0;
    for (int i = 0; i < nk; ++i)
        for (int j = 0; j < nk; ++j) {
            const double r = screenRadius(k_amin[i], k_cmax[i], k_amin[j], k_cmax[j], rabs, I.sparse_eps);
            cut2[static_cast<size_t>(i) * nk + j] = r * r;
            I.sp_rmax = std::max(I.sp_rmax, r);
        }

    // Neighbour atoms per atom (symmetric by construction), ascending.
    const unsigned hw = std::max(1u, std::min(std::thread::hardware_concurrency(), 32u));
    std::vector<std::vector<int>> nb(nat);
    {
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < hw; ++t)
            pool.emplace_back([&, t]() {
                for (int b = static_cast<int>(t); b < nat; b += static_cast<int>(hw)) {
                    const double xb = xyz_bohr[3*b], yb = xyz_bohr[3*b+1], zb = xyz_bohr[3*b+2];
                    for (int a = 0; a < nat; ++a) {
                        const double dx = xyz_bohr[3*a] - xb, dy = xyz_bohr[3*a+1] - yb, dz = xyz_bohr[3*a+2] - zb;
                        if (dx*dx + dy*dy + dz*dz < cut2[static_cast<size_t>(kind[a]) * nk + kind[b]])
                            nb[b].push_back(a);
                    }
                }
            });
        for (auto& th : pool) th.join();
    }

    // Column sizes -> nnz.
    std::vector<long long> col_atom_nnz(nat, 0);
    long long nnz = 0;
    for (int b = 0; b < nat; ++b) {
        for (int a : nb[b]) col_atom_nnz[b] += I.h_at_nao[a];
        nnz += col_atom_nnz[b] * I.h_at_nao[b];
    }
    I.sp_fraction = static_cast<double>(nnz) / (static_cast<double>(nao) * nao);
    const bool want = (I.sparse_mode == 2) || (I.sparse_mode == 1 && I.sp_fraction < 0.5);
    // The device buffers are int-indexed; the quadrupole block is 6 * nnz.
    if (!want || nnz <= 0 || nnz > INT_MAX / 6) return true;

    // colptr per AO (column-major), then row/col/perm.
    std::vector<int> colptr(nao + 1, 0);
    {
        // AOs are contiguous and ascending by atom, so the columns of atom b follow those of b-1.
        long long acc = 0;
        for (int b = 0; b < nat; ++b)
            for (int k = 0; k < I.h_at_nao[b]; ++k) {
                colptr[I.h_at_ao0[b] + k] = static_cast<int>(acc);
                acc += col_atom_nnz[b];
            }
        colptr[nao] = static_cast<int>(acc);
    }
    // AO offset of each neighbour atom block inside a column of atom b.
    std::vector<std::vector<int>> nb_off(nat);
    for (int b = 0; b < nat; ++b) {
        nb_off[b].resize(nb[b].size());
        int off = 0;
        for (size_t i = 0; i < nb[b].size(); ++i) { nb_off[b][i] = off; off += I.h_at_nao[nb[b][i]]; }
    }
    I.h_sp_row.assign(static_cast<size_t>(nnz), 0);
    I.h_sp_col.assign(static_cast<size_t>(nnz), 0);
    std::vector<int> perm(static_cast<size_t>(nnz), 0);
    // Point 2 (Claude Generated, Sep 2026): diag_idx[mu] = the pair-storage index e of the
    // (mu,mu) entry, so the gradient's CN-onsite kernel can read P(mu,mu) straight out of the
    // pattern-only dSpP instead of needing a dense nao x nao buffer. The diagonal AO block is
    // always in the pattern (an atom is its own neighbour, cut2[same-kind] > 0), and the loop
    // below already visits it - when the neighbour atom a equals the column atom b, r==k is
    // exactly mu==nu. Captured in the SAME pass, no extra scan.
    std::vector<int> diag_idx(nao, -1);
    {
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < hw; ++t)
            pool.emplace_back([&, t]() {
                for (int b = static_cast<int>(t); b < nat; b += static_cast<int>(hw)) {
                    for (int k = 0; k < I.h_at_nao[b]; ++k) {
                        const int nu = I.h_at_ao0[b] + k;      // column AO (atom b)
                        long long e = colptr[nu];
                        for (size_t ia = 0; ia < nb[b].size(); ++ia) {
                            const int a = nb[b][ia];
                            // position of atom b inside the neighbour list of atom a
                            const auto& la = nb[a];
                            const size_t pos = static_cast<size_t>(std::lower_bound(la.begin(), la.end(), b) - la.begin());
                            const int boff = nb_off[a][pos];
                            for (int r = 0; r < I.h_at_nao[a]; ++r, ++e) {
                                const int mu = I.h_at_ao0[a] + r;   // row AO (atom a)
                                I.h_sp_row[e] = mu;
                                I.h_sp_col[e] = nu;
                                // transpose (nu, mu): column mu (atom a), row nu (atom b)
                                perm[e] = colptr[mu] + boff + k;
                                if (a == b && r == k) diag_idx[nu] = static_cast<int>(e);
                            }
                        }
                    }
                }
            });
        for (auto& th : pool) th.join();
    }

    try {
        I.dSpRow.upload(I.h_sp_row.data(), static_cast<int>(nnz), I.stream);
        I.dSpCol.upload(I.h_sp_col.data(), static_cast<int>(nnz), I.stream);
        I.dSpPerm.upload(perm.data(), static_cast<int>(nnz), I.stream);
        I.dSpColPtr.upload(colptr.data(), nao + 1, I.stream);
        I.dSpDiagIdx.upload(diag_idx.data(), nao, I.stream);
    } catch (const std::exception& ex) {
        I.last_error = std::string("screened pair upload failed: ") + ex.what();
        return false;
    }
    I.sp_nnz = static_cast<int>(nnz);
    ++I.pattern_generation;
    use_sparse = true;
    // Claude Generated (Sep 29, 2026): the same cutoff table for the dense S/H0 build. Widened by
    // 1e-6 relative so a pair the host put into the pattern can never be skipped on the device
    // because of a different rounding of the distance (FMA contraction); a borderline pair
    // outside the pattern that is computed anyway only adds a value below sparse_eps to the
    // dense S that is factorised. CURCUMA_GPU_FULL_DENSE_S=1 computes every pair (old path).
    I.screen_nk = 0;
    if (!std::getenv("CURCUMA_GPU_FULL_DENSE_S")) {
        std::vector<double> cut2w(cut2.size());
        for (size_t i = 0; i < cut2.size(); ++i) cut2w[i] = cut2[i] * (1.0 + 1.0e-6);
        try {
            I.dAtKind.upload(kind.data(), nat, I.stream);
            I.dScreenCut2.upload(cut2w.data(), nk * nk, I.stream);
            // upload() is asynchronous from pageable host memory; the vectors die on return.
            if (cudaStreamSynchronize(I.stream) != cudaSuccess) return false;
            I.screen_nk = nk;
        } catch (const std::exception&) {
            I.screen_nk = 0;   // unscreened dense build, same result
        }
    }
    return true;
}

bool XtbGpuContext::computeCnSelfEnergy(const double* xyz_bohr)
{
    if (!ok() || m_impl->basis_nat <= 0 || m_impl->basis_nsh <= 0 || !xyz_bohr)
        return false;
    const int nat = m_impl->basis_nat;
    const int nsh = m_impl->basis_nsh;
    cudaStream_t stream = m_impl->stream;

    m_impl->dXyz.upload(xyz_bohr, 3 * nat, stream);

    const int b = 128;
    k_cn<<<(nat + b - 1) / b, b, 0, stream>>>(nat, m_impl->dXyz.ptr, m_impl->dZ.ptr,
                                              m_impl->basis_is_gfn2, m_impl->dCN.ptr);
    if (cudaGetLastError() != cudaSuccess)
        return false;
    k_self_energy<<<(nsh + b - 1) / b, b, 0, stream>>>(nsh, m_impl->dSelfE0.ptr,
                                                       m_impl->dKcn.ptr, m_impl->dSh2at.ptr,
                                                       m_impl->dCN.ptr, m_impl->dSE.ptr);
    if (cudaGetLastError() != cudaSuccess)
        return false;
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::computeIntegrals(const double* xyz_bohr)
{
    if (!ok() || m_impl->basis_nao <= 0 || !xyz_bohr) return false;
    const int nat = m_impl->basis_nat;
    const int nsh = m_impl->basis_nsh;
    const int nao = m_impl->basis_nao;
    cudaStream_t stream = m_impl->stream;
    Impl& I = *m_impl;
    const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);

    // Claude Generated (Sep 2026): choose the storage for this geometry (dense nao^2 or the
    // screened pair list), check that it fits, then allocate. The check runs once per basis
    // and whenever the storage kind changes; later geometries of the same molecule reuse the
    // buffers they already hold (the free-memory figure no longer contains them).
    const bool was_sparse = I.sparse;
    const bool had_storage = !I.dL.empty();
    bool use_sparse = false;
    if (!buildScreenedPairs(xyz_bohr, use_sparse)) return false;
    if (I.memory_check && (!had_storage || use_sparse != was_sparse)) {
        const double nnz = use_sparse ? static_cast<double>(I.sp_nnz) : 0.0;
        if (had_storage) { I.releaseLarge(); I.basis_nao = nao; }
        const size_t need = estimateStorageBytes(nat, nsh, nao, I.basis_is_gfn2 != 0, nnz);
        size_t free_b = 0, total_b = 0;
        if (cudaMemGetInfo(&free_b, &total_b) == cudaSuccess && need > free_b) {
            char msg[400];
            std::snprintf(msg, sizeof(msg),
                          "needs about %.1f GB of device memory for nao=%d, nat=%d (%s storage, "
                          "%.1f %% of AO pairs), but only %.1f of %.1f GB are free on device %d",
                          need / 1073741824.0, nao, nat, use_sparse ? "screened" : "dense",
                          100.0 * I.sp_fraction, free_b / 1073741824.0, total_b / 1073741824.0, I.device);
            I.last_error = msg;
            I.releaseLarge();
            I.basis_nao = nao;
            return false;
        }
        if (had_storage) {
            // releaseLarge dropped the pair arrays built above; rebuild them.
            if (!buildScreenedPairs(xyz_bohr, use_sparse)) return false;
        }
    }
    try {
        I.dL.ensure(static_cast<int>(nn));
        if (use_sparse) {
            // Dense S is built straight into dL (then factorised in place) and dense H0 into
            // dC, which the SCF reuses for the Fock matrix: no extra nao^2 buffer.
            for (CudaBuffer<double>* buf : { &I.dS, &I.dH0, &I.dDpInt, &I.dQpInt }) buf->free();
            I.dC.ensure(static_cast<int>(nn));
            I.dSpS.ensure(I.sp_nnz);
            I.dSpH0.ensure(I.sp_nnz);
            I.dSpTmp.ensure(I.sp_nnz);
            if (I.basis_is_gfn2) {
                I.dSpDp.ensure(3 * I.sp_nnz);
                I.dSpQp.ensure(6 * I.sp_nnz);
            }
        } else {
            for (CudaBuffer<double>* buf : { &I.dSpS, &I.dSpH0, &I.dSpTmp, &I.dSpDp, &I.dSpQp }) buf->free();
            for (CudaBuffer<int>* buf : { &I.dSpRow, &I.dSpCol, &I.dSpColPtr, &I.dSpPerm }) buf->free();
            I.h_sp_row.clear(); I.h_sp_col.clear();
            I.dS.ensure(static_cast<int>(nn));
            I.dH0.ensure(static_cast<int>(nn));
            if (I.basis_is_gfn2) {
                I.dDpInt.ensure(static_cast<int>(3 * nn));
                I.dQpInt.ensure(static_cast<int>(6 * nn));
            }
        }
        int lwork = 0;
        if (cusolverDnDpotrf_bufferSize(I.cusolver, CUBLAS_FILL_MODE_LOWER, nao, I.dL.ptr, nao, &lwork)
            != CUSOLVER_STATUS_SUCCESS)
            return false;
        I.potrf_lwork = lwork;
        I.dPotrfWork.ensure(lwork > 0 ? lwork : 1);
    } catch (const std::exception& ex) {
        I.last_error = std::string("device allocation failed: ") + ex.what();
        I.releaseLarge();
        I.basis_nao = nao;
        return false;
    }
    I.sparse = use_sparse;

    // CN + self-energies (uploads xyz, runs k_cn + k_self_energy, syncs).
    I.profStart();
    if (!computeCnSelfEnergy(xyz_bohr)) return false;
    I.profMark("integrals: CN + self-energies");

    // Overlap S + bare Hamiltonian H0 (one thread per shell-pair), dense.
    const double* S_dense_dst = use_sparse ? I.dL.ptr : I.dS.ptr;
    const double* H_dense_dst = use_sparse ? I.dC.ptr : I.dH0.ptr;
    const dim3 block(16, 16);
    const dim3 grid((nsh + block.x - 1) / block.x, (nsh + block.y - 1) / block.y);
    const int screen_nk = use_sparse ? I.screen_nk : 0;
    if (screen_nk > 0
        && (cudaMemsetAsync(const_cast<double*>(S_dense_dst), 0, sizeof(double) * nn, stream) != cudaSuccess
            || cudaMemsetAsync(const_cast<double*>(H_dense_dst), 0, sizeof(double) * nn, stream) != cudaSuccess))
        return false;
    k_overlap_h0<<<grid, block, 0, stream>>>(
        nsh, nao, I.basis_is_gfn2,
        I.dSh2at.ptr, I.dAng.ptr, I.dIaoSh.ptr, I.dNaoSh.ptr,
        I.dShNprim.ptr, I.dShPrimOff.ptr, I.dPrimAlpha.ptr,
        I.dPrimCoeff.ptr, I.dShZeta.ptr, I.dShpoly.ptr,
        I.dSE.ptr, I.dZ.ptr, I.dValence.ptr, I.dXyz.ptr,
        const_cast<double*>(S_dense_dst), const_cast<double*>(H_dense_dst),
        screen_nk > 0 ? I.dAtKind.ptr : nullptr, screen_nk > 0 ? I.dScreenCut2.ptr : nullptr, screen_nk);
    if (cudaGetLastError() != cudaSuccess) return false;
    if (use_sparse) {
        const int b1 = 256;
        const int g1 = (I.sp_nnz + b1 - 1) / b1;
        k_sp_gather<<<g1, b1, 0, stream>>>(I.dL.ptr, I.dSpRow.ptr, I.dSpCol.ptr, nao, I.sp_nnz, I.dSpS.ptr);
        k_sp_gather<<<g1, b1, 0, stream>>>(I.dC.ptr, I.dSpRow.ptr, I.dSpCol.ptr, nao, I.sp_nnz, I.dSpH0.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    I.profMark("integrals: overlap + H0");

    // Coulomb γ matrix (independent of S; one thread per shell-pair).
    if (!I.dHardness.empty()) {
        const dim3 gblock(16, 16);
        const dim3 ggrid((nsh + gblock.x - 1) / gblock.x, (nsh + gblock.y - 1) / gblock.y);
        k_gamma<<<ggrid, gblock, 0, stream>>>(nsh, I.basis_is_gfn2, I.dSh2at.ptr,
                                              I.dHardness.ptr, I.dXyz.ptr, I.dGamma.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
    }

    // GFN2 multipole integrals (dp_int/qp_int); need S for the origin shift.
    if (I.basis_is_gfn2) {
        if (use_sparse) {
            const int b1 = 256;
            k_multipole_ints_sp<<<(I.sp_nnz + b1 - 1) / b1, b1, 0, stream>>>(
                I.sp_nnz, I.dSpRow.ptr, I.dSpCol.ptr, I.dAo2sh.ptr, I.dAo2at.ptr, I.dIaoSh.ptr, I.dAng.ptr,
                I.dShNprim.ptr, I.dShPrimOff.ptr, I.dPrimAlpha.ptr, I.dPrimCoeff.ptr, I.dXyz.ptr,
                I.dSpS.ptr, I.dSpDp.ptr, I.dSpQp.ptr);
            if (cudaGetLastError() != cudaSuccess) return false;
        } else if (!I.dDpInt.empty()) {
            const dim3 mblock(16, 16);
            const dim3 mgrid((nao + mblock.x - 1) / mblock.x, (nao + mblock.y - 1) / mblock.y);
            k_multipole_ints<<<mgrid, mblock, 0, stream>>>(
                nao, I.dAo2sh.ptr, I.dAo2at.ptr, I.dIaoSh.ptr, I.dAng.ptr,
                I.dShNprim.ptr, I.dShPrimOff.ptr, I.dPrimAlpha.ptr,
                I.dPrimCoeff.ptr, I.dXyz.ptr, I.dS.ptr,
                I.dDpInt.ptr, I.dQpInt.ptr);
            if (cudaGetLastError() != cudaSuccess) return false;
        }
    }
    I.profMark("integrals: gamma + multipole");

    // L = chol(S): factor the lower triangle in place (matches the CPU Eigen matrixL()
    // convention; the trsm path reads fill=LOWER only). Dense storage copies S -> L first.
    if (!use_sparse) {
        if (cudaMemcpyAsync(I.dL.ptr, I.dS.ptr, sizeof(double) * nn,
                            cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
            return false;
    }
    const int dist_chol = I.distCholesky(nao);
    if (dist_chol < 0) return false;
    int info = 1;
    if (dist_chol == 1) {
        info = 0;
        I.profMark("integrals: Cholesky of S (multi-GPU)");
    } else {
        if (cusolverDnDpotrf(I.cusolver, CUBLAS_FILL_MODE_LOWER, nao,
                             I.dL.ptr, nao, I.dPotrfWork.ptr,
                             I.potrf_lwork, I.dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
            return false;
        if (cudaMemcpyAsync(&info, I.dInfo.ptr, sizeof(int),
                            cudaMemcpyDeviceToHost, stream) != cudaSuccess)
            return false;
        if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
        I.profMark("integrals: Cholesky of S");
    }
    ++I.l_generation;
    (void)nat;
    return info == 0;
}

// ---- Storage-independent SCF building blocks (Claude Generated, Sep 2026) ----------------
// Each dispatches to the dense nao^2 kernel or its screened-pair twin; the arithmetic per
// matrix element is the same in both.

bool XtbGpuContext::buildFockIntoC(int n, bool multipole)
{
    Impl& I = *m_impl;
    cudaStream_t stream = I.stream;
    if (I.sparse) {
        const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
        if (cudaMemsetAsync(I.dC.ptr, 0, sizeof(double) * nn, stream) != cudaSuccess) return false;
        const bool mp = multipole && !I.dSpDp.empty();
        const int bs = 256;
        k_build_fock_sp<<<(I.sp_nnz + bs - 1) / bs, bs, 0, stream>>>(
            I.dC.ptr, n, I.sp_nnz, I.dSpRow.ptr, I.dSpCol.ptr, I.dSpPerm.ptr,
            I.dSpH0.ptr, I.dSpS.ptr, I.dVao.ptr,
            mp ? I.dSpDp.ptr : nullptr, mp ? I.dSpQp.ptr : nullptr,
            mp ? I.dVdp.ptr : nullptr, mp ? I.dVqp.ptr : nullptr, I.dAo2at.ptr);
        return cudaGetLastError() == cudaSuccess;
    }
    const dim3 block(16, 16);
    const dim3 grid((n + block.x - 1) / block.x, (n + block.y - 1) / block.y);
    k_build_fock_iso<<<grid, block, 0, stream>>>(I.dC.ptr, I.dH0.ptr, I.dS.ptr, I.dVao.ptr, n);
    if (cudaGetLastError() != cudaSuccess) return false;
    if (multipole) {
        k_add_fock_multipole<<<grid, block, 0, stream>>>(I.dC.ptr, I.dDpInt.ptr, I.dQpInt.ptr,
                                                         I.dVdp.ptr, I.dVqp.ptr, I.dAo2at.ptr, n);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    return true;
}

bool XtbGpuContext::populationsAndBand(int n, double* band_out)
{
    Impl& I = *m_impl;
    cudaStream_t stream = I.stream;
    const int b1 = 128;
    if (I.sparse && !I.p_dense_valid) {
        k_pop_ao_spP<<<(n + b1 - 1) / b1, b1, 0, stream>>>(I.dPop.ptr, I.dSpP.ptr, I.dSpS.ptr,
                                                           I.dSpColPtr.ptr, I.dSpPerm.ptr, n);
        if (cudaGetLastError() != cudaSuccess) return false;
        return cublasDdot(I.cublas, I.sp_nnz, I.dSpP.ptr, 1, I.dSpH0.ptr, 1, band_out)
            == CUBLAS_STATUS_SUCCESS;
    }
    if (I.sparse) {
        k_pop_ao_sp<<<(n + b1 - 1) / b1, b1, 0, stream>>>(I.dPop.ptr, I.dP.ptr, I.dSpS.ptr,
                                                          I.dSpRow.ptr, I.dSpColPtr.ptr, I.dSpPerm.ptr, n);
        if (cudaGetLastError() != cudaSuccess) return false;
        // Band energy over the stored pairs: gather P, then dot with H0.
        const int bs = 256;
        k_sp_gather<<<(I.sp_nnz + bs - 1) / bs, bs, 0, stream>>>(I.dP.ptr, I.dSpRow.ptr, I.dSpCol.ptr,
                                                                 n, I.sp_nnz, I.dSpTmp.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
        return cublasDdot(I.cublas, I.sp_nnz, I.dSpTmp.ptr, 1, I.dSpH0.ptr, 1, band_out)
            == CUBLAS_STATUS_SUCCESS;
    }
    k_pop_ao<<<(n + b1 - 1) / b1, b1, 0, stream>>>(I.dPop.ptr, I.dP.ptr, I.dS.ptr, n);
    if (cudaGetLastError() != cudaSuccess) return false;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    return cublasDdot(I.cublas, static_cast<int>(nn), I.dP.ptr, 1, I.dH0.ptr, 1, band_out)
        == CUBLAS_STATUS_SUCCESS;
}

bool XtbGpuContext::multipoleMomentsResident(int n, int nat)
{
    Impl& I = *m_impl;
    cudaStream_t stream = I.stream;
    // Accumulator-scatter kernels need zeroed targets.
    if (cudaMemsetAsync(I.dDpAt.ptr, 0, sizeof(double) * 3 * nat, stream) != cudaSuccess) return false;
    if (cudaMemsetAsync(I.dQpAt.ptr, 0, sizeof(double) * 6 * nat, stream) != cudaSuccess) return false;
    const int b1 = 128;
    if (I.sparse && !I.p_dense_valid) {
        k_multipole_moments_spP<<<(n + b1 - 1) / b1, b1, 0, stream>>>(
            I.dDpAt.ptr, I.dQpAt.ptr, I.dSpP.ptr, I.dSpDp.ptr, I.dSpQp.ptr,
            I.dSpColPtr.ptr, I.dAo2at.ptr, n, I.sp_nnz);
    } else if (I.sparse) {
        k_multipole_moments_sp<<<(n + b1 - 1) / b1, b1, 0, stream>>>(
            I.dDpAt.ptr, I.dQpAt.ptr, I.dP.ptr, I.dSpDp.ptr, I.dSpQp.ptr,
            I.dSpRow.ptr, I.dSpColPtr.ptr, I.dAo2at.ptr, n, I.sp_nnz);
    } else {
        k_multipole_moments<<<(n + b1 - 1) / b1, b1, 0, stream>>>(
            I.dDpAt.ptr, I.dQpAt.ptr, I.dP.ptr, I.dDpInt.ptr, I.dQpInt.ptr, I.dAo2at.ptr, n);
    }
    return cudaGetLastError() == cudaSuccess;
}

bool XtbGpuContext::residentBeginComputed()
{
    if (!ok() || m_impl->basis_nao <= 0) return false;
    const int n = m_impl->basis_nao;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    try {
        // dH0/dS/dL already hold the device-computed integrals (computeIntegrals).
        // Allocate the resident SCF work buffers, like residentBegin but without
        // uploading any matrix. ensure() reuses the allocation across MD/opt steps.
        // Claude Generated (Sep 2026): Cw (n x ncol) and the FP64/FP32 eigensolver workspaces
        // are allocated on first use by the density build and eigensolveResidentFock, which
        // keep only one precision's workspace at a time.
        m_impl->dC.ensure(static_cast<int>(nn));
        // Screened storage keeps P on the pattern during the loop; dense P only on demand.
        if (!m_impl->sparse) m_impl->dP.ensure(static_cast<int>(nn));
        else m_impl->dP.free();
        m_impl->dEps.ensure(n);
        m_impl->dVao.ensure(n);
        m_impl->dOcc.ensure(n);
        m_impl->dPop.ensure(n);
        if (m_impl->dInfo.n < 1) m_impl->dInfo.alloc(1);
    } catch (const std::exception& e) {
        m_impl->last_error = std::string("device allocation failed: ") + e.what();
        m_impl->releaseLarge();
        return false;
    }
    m_impl->resident_n = n;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadCn(double* cn_out)
{
    if (!ok() || m_impl->basis_nat <= 0 || !cn_out) return false;
    m_impl->dCN.download(cn_out, m_impl->basis_nat, m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadSelfEnergy(double* se_out)
{
    if (!ok() || m_impl->basis_nsh <= 0 || !se_out) return false;
    m_impl->dSE.download(se_out, m_impl->basis_nsh, m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

// Screened storage: scatter the stored pairs into a zeroed host dense matrix (the host keeps
// dense S/H0 for properties and the CPU fallback). Claude Generated (Sep 2026).
static bool scatterToHost(CudaBuffer<double>& values, const std::vector<int>& row,
                          const std::vector<int>& col, int n, int ncomp, int nnz,
                          double* out, cudaStream_t stream)
{
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    std::vector<double> v(static_cast<size_t>(ncomp) * nnz);
    values.download(v.data(), ncomp * nnz, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
    // Claude Generated (Sep 29, 2026): zero fill and scatter split over host threads. Every
    // stored pair has its own target element, so the ranges never collide and the result is
    // identical to the serial loop (polymer_2x: 313 ms per matrix serial).
    const unsigned nt = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
    auto parallel = [nt](size_t count, const auto& body) {
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < nt; ++t)
            pool.emplace_back([&, t]() { body(count * t / nt, count * (t + 1) / nt); });
        for (auto& th : pool) th.join();
    };
    parallel(static_cast<size_t>(ncomp) * nn,
             [&](size_t a, size_t b) { std::fill(out + a, out + b, 0.0); });
    for (int k = 0; k < ncomp; ++k) {
        const double* src = v.data() + static_cast<size_t>(k) * nnz;
        double* dst = out + static_cast<size_t>(k) * nn;
        parallel(static_cast<size_t>(nnz), [&](size_t a, size_t b) {
            for (size_t e = a; e < b; ++e)
                dst[static_cast<size_t>(row[e]) + static_cast<size_t>(col[e]) * n] = src[e];
        });
    }
    return true;
}

bool XtbGpuContext::downloadOverlap(double* S_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !S_out) return false;
    if (m_impl->sparse)
        return scatterToHost(m_impl->dSpS, m_impl->h_sp_row, m_impl->h_sp_col, m_impl->basis_nao, 1,
                             m_impl->sp_nnz, S_out, m_impl->stream);
    const size_t nn = static_cast<size_t>(m_impl->basis_nao) * m_impl->basis_nao;
    m_impl->dS.download(S_out, static_cast<int>(nn), m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadH0(double* H0_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !H0_out) return false;
    if (m_impl->sparse)
        return scatterToHost(m_impl->dSpH0, m_impl->h_sp_row, m_impl->h_sp_col, m_impl->basis_nao, 1,
                             m_impl->sp_nnz, H0_out, m_impl->stream);
    const size_t nn = static_cast<size_t>(m_impl->basis_nao) * m_impl->basis_nao;
    m_impl->dH0.download(H0_out, static_cast<int>(nn), m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadCholesky(double* L_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !L_out) return false;
    const size_t nn = static_cast<size_t>(m_impl->basis_nao) * m_impl->basis_nao;
    m_impl->dL.download(L_out, static_cast<int>(nn), m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadGamma(double* gamma_out)
{
    if (!ok() || m_impl->basis_nsh <= 0 || m_impl->dGamma.empty() || !gamma_out) return false;
    const size_t nn = static_cast<size_t>(m_impl->basis_nsh) * m_impl->basis_nsh;
    m_impl->dGamma.download(gamma_out, static_cast<int>(nn), m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadMultipoleInts(double* dp_int3, double* qp_int6)
{
    if (!ok() || m_impl->basis_nao <= 0 || !dp_int3 || !qp_int6)
        return false;
    if (m_impl->sparse) {
        if (m_impl->dSpDp.empty()) return false;
        return scatterToHost(m_impl->dSpDp, m_impl->h_sp_row, m_impl->h_sp_col, m_impl->basis_nao, 3,
                             m_impl->sp_nnz, dp_int3, m_impl->stream)
            && scatterToHost(m_impl->dSpQp, m_impl->h_sp_row, m_impl->h_sp_col, m_impl->basis_nao, 6,
                             m_impl->sp_nnz, qp_int6, m_impl->stream);
    }
    if (m_impl->dDpInt.empty()) return false;
    const size_t nn = static_cast<size_t>(m_impl->basis_nao) * m_impl->basis_nao;
    m_impl->dDpInt.download(dp_int3, static_cast<int>(3 * nn), m_impl->stream);
    m_impl->dQpInt.download(qp_int6, static_cast<int>(6 * nn), m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Stage 5 (Part A): single-shot D4 EEQ charge model on the device.
 *  Self-contained (no beginBasis): builds CN + the (N+1) augmented matrix +
 *  RHS on the device, factors with cusolverDnDgetrf (partial-pivot LU — the
 *  augmented system is symmetric *indefinite*, so LU, not Cholesky), solves
 *  with getrs. The LU factor, CN and per-atom params stay resident so
 *  eeqChargeResponseGradient reuses them. Claude Generated.
 * ====================================================================== */
bool XtbGpuContext::eeqCharges(int N, const double* xyz_bohr,
                               const double* chi, const double* gam,
                               const double* alpha_sq, const double* cnf,
                               const double* rcov_bohr, double total_charge,
                               double* q_out)
{
    if (!ok() || N <= 0 || !xyz_bohr || !chi || !gam || !alpha_sq || !cnf
        || !rcov_bohr || !q_out)
        return false;
    const int m = N + 1;
    cudaStream_t stream = m_impl->stream;
    try {
        m_impl->dEeqXyz.ensure(3 * N);
        m_impl->dEeqChi.ensure(N);
        m_impl->dEeqGam.ensure(N);
        m_impl->dEeqAlp.ensure(N);
        m_impl->dEeqCnf.ensure(N);
        m_impl->dEeqRcov.ensure(N);
        m_impl->dEeqCn.ensure(N);
        m_impl->dEeqCnRaw.ensure(N);
        m_impl->dEeqM.ensure(m * m);
        m_impl->dEeqRhs.ensure(m);
        m_impl->dEeqQ.ensure(N);
        m_impl->dEeqIpiv.ensure(m);
        if (m_impl->dInfo.n < 1) m_impl->dInfo.alloc(1);
    } catch (...) {
        return false;
    }

    m_impl->dEeqXyz.upload(xyz_bohr, 3 * N, stream);
    m_impl->dEeqChi.upload(chi, N, stream);
    m_impl->dEeqGam.upload(gam, N, stream);
    m_impl->dEeqAlp.upload(alpha_sq, N, stream);
    m_impl->dEeqCnf.upload(cnf, N, stream);
    m_impl->dEeqRcov.upload(rcov_bohr, N, stream);

    const int b = 128;
    k_d4eeq_cn<<<(N + b - 1) / b, b, 0, stream>>>(N, m_impl->dEeqXyz.ptr,
                                                  m_impl->dEeqRcov.ptr,
                                                  m_impl->dEeqCn.ptr, m_impl->dEeqCnRaw.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    const dim3 blk(16, 16);
    const dim3 grd((m + blk.x - 1) / blk.x, (m + blk.y - 1) / blk.y);
    k_d4eeq_build<<<grd, blk, 0, stream>>>(N, m_impl->dEeqXyz.ptr, m_impl->dEeqAlp.ptr,
                                           m_impl->dEeqGam.ptr, m_impl->dEeqM.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    k_d4eeq_rhs<<<(m + b - 1) / b, b, 0, stream>>>(N, m_impl->dEeqChi.ptr, m_impl->dEeqCnf.ptr,
                                                   m_impl->dEeqCn.ptr, total_charge,
                                                   m_impl->dEeqRhs.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    // LU factor M (getrf) then solve M·[q;λ] = [b;Q] (getrs). Workspace query is
    // size-stable for fixed m, so the buffer is reused across geometry steps.
    int lwork = 0;
    if (cusolverDnDgetrf_bufferSize(m_impl->cusolver, m, m, m_impl->dEeqM.ptr, m, &lwork)
        != CUSOLVER_STATUS_SUCCESS)
        return false;
    try {
        m_impl->eeq_lwork = lwork;
        m_impl->dEeqWork.ensure(lwork > 0 ? lwork : 1);
    } catch (...) {
        return false;
    }
    if (cusolverDnDgetrf(m_impl->cusolver, m, m, m_impl->dEeqM.ptr, m,
                         m_impl->dEeqWork.ptr, m_impl->dEeqIpiv.ptr, m_impl->dInfo.ptr)
        != CUSOLVER_STATUS_SUCCESS)
        return false;
    if (cusolverDnDgetrs(m_impl->cusolver, CUBLAS_OP_N, m, 1, m_impl->dEeqM.ptr, m,
                         m_impl->dEeqIpiv.ptr, m_impl->dEeqRhs.ptr, m, m_impl->dInfo.ptr)
        != CUSOLVER_STATUS_SUCCESS)
        return false;
    // rhs head N holds the atomic charges q → keep resident for the response.
    if (cudaMemcpyAsync(m_impl->dEeqQ.ptr, m_impl->dEeqRhs.ptr, sizeof(double) * N,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return false;
    int info = 1;
    if (cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
                        cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    m_impl->dEeqQ.download(q_out, N, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
    m_impl->eeq_n = (info == 0) ? N : 0;
    return info == 0;
}

bool XtbGpuContext::eeqChargeResponseGradient(int N, const double* dEdq, double* grad_add)
{
    if (!ok() || N <= 0 || N != m_impl->eeq_n || !dEdq || !grad_add)
        return false;
    const int m = N + 1;
    cudaStream_t stream = m_impl->stream;
    if (N <= 1) {                       // no pairwise geometry dependence
        for (int i = 0; i < 3 * N; ++i) grad_add[i] = 0.0;
        return true;
    }
    try {
        m_impl->dEeqDedq.ensure(N);
        m_impl->dEeqAdjRhs.ensure(m);
        m_impl->dEeqU.ensure(N);
        m_impl->dEeqGrad.ensure(3 * N);
    } catch (...) {
        return false;
    }
    m_impl->dEeqDedq.upload(dEdq, N, stream);

    const int b = 128;
    k_d4eeq_adjoint_rhs<<<(m + b - 1) / b, b, 0, stream>>>(N, m_impl->dEeqDedq.ptr,
                                                           m_impl->dEeqAdjRhs.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    // Adjoint: M·z = [dEdq;0], reusing the LU factor + pivots from eeqCharges.
    if (cusolverDnDgetrs(m_impl->cusolver, CUBLAS_OP_N, m, 1, m_impl->dEeqM.ptr, m,
                         m_impl->dEeqIpiv.ptr, m_impl->dEeqAdjRhs.ptr, m, m_impl->dInfo.ptr)
        != CUSOLVER_STATUS_SUCCESS)
        return false;
    // dEeqAdjRhs head N = z_q. Build the per-atom CN-response weight then the pairs.
    k_d4eeq_u<<<(N + b - 1) / b, b, 0, stream>>>(N, m_impl->dEeqAdjRhs.ptr, m_impl->dEeqCnf.ptr,
                                                 m_impl->dEeqCn.ptr, m_impl->dEeqCnRaw.ptr,
                                                 m_impl->dEeqU.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    k_d4eeq_response<<<(N + b - 1) / b, b, 0, stream>>>(N, m_impl->dEeqXyz.ptr, m_impl->dEeqAlp.ptr,
                                                        m_impl->dEeqRcov.ptr, m_impl->dEeqQ.ptr,
                                                        m_impl->dEeqAdjRhs.ptr, m_impl->dEeqU.ptr,
                                                        m_impl->dEeqGrad.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    // Download the N×3 [3a+k] response contribution; the host adds it to its accumulator.
    m_impl->dEeqGrad.download(grad_add, 3 * N, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Stage 5 (Part B1): atomic Mulliken charges from the resident density.
 *  Reduces the resident dPop (= pop_ao, set by residentDensity) into the
 *  resident dQat via the resident AO→atom map; q_at_out optionally downloads
 *  the result for validation. dQat stays resident for the in-SCF D4 potential.
 * ====================================================================== */
bool XtbGpuContext::residentAtomicCharges(const double* n0_at, int nat, double* q_at_out)
{
    if (!ok() || nat <= 0 || m_impl->resident_n <= 0 || m_impl->dPop.empty()
        || m_impl->dAo2at.empty() || !n0_at)
        return false;
    const int nao = m_impl->resident_n;
    cudaStream_t stream = m_impl->stream;
    try {
        m_impl->dN0at.ensure(nat);
        m_impl->dQat.ensure(nat);
    } catch (...) {
        return false;
    }
    // q_at ← n0_at, then subtract each AO population into its atom bin.
    m_impl->dN0at.upload(n0_at, nat, stream);
    if (cudaMemcpyAsync(m_impl->dQat.ptr, m_impl->dN0at.ptr, sizeof(double) * nat,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return false;
    const int b = 128;
    k_qat_scatter<<<(nao + b - 1) / b, b, 0, stream>>>(nao, m_impl->dPop.ptr,
                                                       m_impl->dAo2at.ptr, m_impl->dQat.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    if (q_at_out) m_impl->dQat.download(q_at_out, nat, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

// Stage 6 (S6.2): shell Mulliken charges from the resident density populations.
// q_sh ← n0_sh, then subtract each AO population into its shell bin (dAo2sh). The
// shell half of updatePopulationsFromPopAo, mirroring residentAtomicCharges. The
// result stays resident in dQsh for the device SCC energy (S6.3) + Broyden (S6.4);
// q_sh_out optionally downloads it for validation. Requires the resident dAo2sh
// (Stage-3 device-integral path). Claude Generated.
bool XtbGpuContext::residentShellCharges(const double* n0_sh, int nsh, double* q_sh_out)
{
    if (!ok() || nsh <= 0 || m_impl->resident_n <= 0 || m_impl->dPop.empty()
        || m_impl->dAo2sh.empty() || !n0_sh)
        return false;
    const int nao = m_impl->resident_n;
    cudaStream_t stream = m_impl->stream;
    try {
        m_impl->dN0sh.ensure(nsh);
        m_impl->dQsh.ensure(nsh);
    } catch (...) {
        return false;
    }
    m_impl->dN0sh.upload(n0_sh, nsh, stream);
    if (cudaMemcpyAsync(m_impl->dQsh.ptr, m_impl->dN0sh.ptr, sizeof(double) * nsh,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return false;
    const int b = 128;
    k_qsh_scatter<<<(nao + b - 1) / b, b, 0, stream>>>(nao, m_impl->dPop.ptr,
                                                       m_impl->dAo2sh.ptr, m_impl->dQsh.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    if (q_sh_out) m_impl->dQsh.download(q_sh_out, nsh, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Stage 5 (Part B2): in-SCF GFN2 D4 atom-potential dE_D4/dq on the device.
 *  beginDispersion uploads the geometry-fixed reference data once per geometry
 *  (c6_flat is element data → uploaded once per process); dispersionDedq runs
 *  the per-iteration O(N²) contraction + BJ disp_sum from the host-built W/dWq.
 * ====================================================================== */
bool XtbGpuContext::beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                                    const int* nref, const double* xyz_bohr,
                                    const double* c6_flat, int c6_flat_len,
                                    double s6, double s8, double a1, double a2, double cutoff)
{
    if (!ok() || nat <= 0 || !Z || !sqrtZr4r2 || !nref || !xyz_bohr
        || !c6_flat || c6_flat_len <= 0)
        return false;
    cudaStream_t stream = m_impl->stream;
    try {
        m_impl->dD4Z.ensure(nat);
        m_impl->dD4Nref.ensure(nat);
        m_impl->dD4Sqrt.ensure(nat);
        m_impl->dD4Xyz.ensure(3 * nat);
        m_impl->dD4W.ensure(nat * Impl::D4_MAX_REF);
        m_impl->dD4dWq.ensure(nat * Impl::D4_MAX_REF);
        m_impl->dD4Dedq.ensure(nat);
        if (!m_impl->d4_c6_uploaded || m_impl->dD4C6Flat.n < c6_flat_len)
            m_impl->dD4C6Flat.ensure(c6_flat_len);
    } catch (...) {
        return false;
    }
    m_impl->dD4Z.upload(Z, nat, stream);
    m_impl->dD4Nref.upload(nref, nat, stream);
    m_impl->dD4Sqrt.upload(sqrtZr4r2, nat, stream);
    m_impl->dD4Xyz.upload(xyz_bohr, 3 * nat, stream);
    m_impl->h_d4_xyz.assign(xyz_bohr, xyz_bohr + 3 * nat);
    // The reference C6 block is element data (geometry- AND molecule-independent):
    // upload it only once per process. ensure() keeps the allocation across steps.
    if (!m_impl->d4_c6_uploaded) {
        m_impl->dD4C6Flat.upload(c6_flat, c6_flat_len, stream);
        m_impl->d4_c6_uploaded = true;
    }
    m_impl->d4_nat = nat;
    m_impl->d4_s6 = s6; m_impl->d4_s8 = s8; m_impl->d4_a1 = a1; m_impl->d4_a2 = a2;
    m_impl->d4_cut = cutoff;
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out)
{
    if (!ok() || nat <= 0 || nat != m_impl->d4_nat || !W || !dWq || !dEdq_out
        || m_impl->dD4C6Flat.empty())
        return false;
    cudaStream_t stream = m_impl->stream;
    m_impl->dD4W.upload(W, nat * Impl::D4_MAX_REF, stream);
    m_impl->dD4dWq.upload(dWq, nat * Impl::D4_MAX_REF, stream);
    const double cut2 = m_impl->d4_cut * m_impl->d4_cut;
    const int b = 128;
    k_d4_dedq<<<(nat + b - 1) / b, b, 0, stream>>>(
        nat, 118, Impl::D4_MAX_REF, m_impl->dD4Z.ptr, m_impl->dD4Nref.ptr,
        m_impl->dD4Sqrt.ptr, m_impl->dD4Xyz.ptr, m_impl->dD4C6Flat.ptr,
        m_impl->dD4W.ptr, m_impl->dD4dWq.ptr,
        m_impl->d4_s6, m_impl->d4_s8, m_impl->d4_a1, m_impl->d4_a2, cut2,
        m_impl->dD4Dedq.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    m_impl->dD4Dedq.download(dEdq_out, nat, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

// Claude Generated (Sep 2026): post-SCF D4 on the device (CUDA twins of the ROCm path). For a
// 7k-atom GFN2 run the host ATM term alone was ~90 % of the 40 s post-SCF phase.
bool XtbGpuContext::dispersionGradient(int nat, const double* W, const double* dWq, const double* dWc,
                                       double* e_atom_out, double* grad_out,
                                       double* dEdcn_out, double* dEdq_out)
{
    if (!ok() || nat <= 0 || nat != m_impl->d4_nat || !W || !dWq || !dWc || !e_atom_out
        || !grad_out || !dEdcn_out || !dEdq_out || m_impl->dD4C6Flat.empty())
        return false;
    Impl& I = *m_impl;
    cudaStream_t stream = I.stream;
    const int wlen = nat * Impl::D4_MAX_REF;
    try {
        I.dD4dWc.ensure(wlen);
        I.dD4Eat.ensure(nat);
        I.dD4Grad.ensure(3 * nat);
        I.dD4Dcn.ensure(nat);
    } catch (...) {
        return false;
    }
    I.dD4W.upload(W, wlen, stream);
    I.dD4dWq.upload(dWq, wlen, stream);
    I.dD4dWc.upload(dWc, wlen, stream);
    const double cut2 = I.d4_cut * I.d4_cut;
    const int b = 128;
    k_d4_grad<<<(nat + b - 1) / b, b, 0, stream>>>(
        nat, 118, Impl::D4_MAX_REF, I.dD4Z.ptr, I.dD4Nref.ptr, I.dD4Sqrt.ptr, I.dD4Xyz.ptr,
        I.dD4C6Flat.ptr, I.dD4W.ptr, I.dD4dWq.ptr, I.dD4dWc.ptr,
        I.d4_s6, I.d4_s8, I.d4_a1, I.d4_a2, cut2,
        I.dD4Eat.ptr, I.dD4Grad.ptr, I.dD4Dcn.ptr, I.dD4Dedq.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    I.dD4Eat.download(e_atom_out, nat, stream);
    I.dD4Grad.download(grad_out, 3 * nat, stream);
    I.dD4Dcn.download(dEdcn_out, nat, stream);
    I.dD4Dedq.download(dEdq_out, nat, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::dispersionATM(int nat, const double* c6, const double* dc6dcn,
                                  double s9, double a1, double a2, double alp, double cutoff,
                                  double* e_atom_out, double* grad_out, double* dEdcn_out)
{
    Impl& I = *m_impl;
    if (!ok() || nat <= 0 || nat != I.d4_nat || !c6 || !dc6dcn || !e_atom_out || !grad_out
        || !dEdcn_out || static_cast<int>(I.h_d4_xyz.size()) != 3 * nat)
        return false;
    cudaStream_t stream = I.stream;
    const double cut2 = cutoff * cutoff;

    // Neighbour list (host, threaded): atoms within the cutoff, ascending, self excluded.
    std::vector<std::vector<int>> nbl(nat);
    {
        const double* xyz = I.h_d4_xyz.data();
        const unsigned hw = std::max(1u, std::min(std::thread::hardware_concurrency(), 32u));
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < hw; ++t)
            pool.emplace_back([&, t]() {
                for (int a = static_cast<int>(t); a < nat; a += static_cast<int>(hw))
                    for (int x = 0; x < nat; ++x) {
                        if (x == a) continue;
                        const double dx = xyz[3*x] - xyz[3*a], dy = xyz[3*x+1] - xyz[3*a+1], dz = xyz[3*x+2] - xyz[3*a+2];
                        if (dx*dx + dy*dy + dz*dz <= cut2) nbl[a].push_back(x);
                    }
            });
        for (auto& th : pool) th.join();
    }
    std::vector<int> nbptr(nat + 1, 0), nb;
    size_t total = 0;
    for (int a = 0; a < nat; ++a) total += nbl[a].size();
    if (total > static_cast<size_t>(INT_MAX)) return false;
    nb.reserve(total);
    for (int a = 0; a < nat; ++a) {
        nbptr[a] = static_cast<int>(nb.size());
        nb.insert(nb.end(), nbl[a].begin(), nbl[a].end());
    }
    nbptr[nat] = static_cast<int>(nb.size());

    const size_t nn = static_cast<size_t>(nat) * nat;
    if (nn > static_cast<size_t>(INT_MAX)) return false;
    try {
        I.dD4AtmC6.ensure(static_cast<int>(nn));
        I.dD4AtmDc6.ensure(static_cast<int>(nn));
        I.dD4Eat.ensure(nat);
        I.dD4Grad.ensure(3 * nat);
        I.dD4Dcn.ensure(nat);
        I.dD4NbPtr.upload(nbptr.data(), nat + 1, stream);
        I.dD4Nb.upload(nb.empty() ? nbptr.data() : nb.data(), std::max<int>(1, static_cast<int>(nb.size())), stream);
    } catch (...) {
        return false;
    }
    I.dD4AtmC6.upload(c6, static_cast<int>(nn), stream);
    I.dD4AtmDc6.upload(dc6dcn, static_cast<int>(nn), stream);
    if (I.thread_per_atom) {
        const int b = 64;
        k_d4_atm_nl<<<(nat + b - 1) / b, b, 0, stream>>>(
            nat, I.dD4Xyz.ptr, I.dD4Sqrt.ptr, I.dD4AtmC6.ptr, I.dD4AtmDc6.ptr, I.dD4NbPtr.ptr, I.dD4Nb.ptr,
            s9, a1, a2, alp, cut2, I.dD4Eat.ptr, I.dD4Grad.ptr, I.dD4Dcn.ptr);
    } else {
        const int b = 128;                                   // 4 atoms (warps) per block
        const long long threads = 32LL * nat;
        k_d4_atm_nl_w<<<static_cast<unsigned>((threads + b - 1) / b), b, 0, stream>>>(
            nat, I.dD4Xyz.ptr, I.dD4Sqrt.ptr, I.dD4AtmC6.ptr, I.dD4AtmDc6.ptr, I.dD4NbPtr.ptr, I.dD4Nb.ptr,
            s9, a1, a2, alp, cut2, I.dD4Eat.ptr, I.dD4Grad.ptr, I.dD4Dcn.ptr);
    }
    if (cudaGetLastError() != cudaSuccess) return false;
    I.dD4Eat.download(e_atom_out, nat, stream);
    I.dD4Grad.download(grad_out, 3 * nat, stream);
    I.dD4Dcn.download(dEdcn_out, nat, stream);
    const bool okk = cudaStreamSynchronize(stream) == cudaSuccess;
    // The nat^2 reference matrices are only needed for this call.
    I.dD4AtmC6.free();
    I.dD4AtmDc6.free();
    return okk;
}

/* ====================================================================== *
 *  Stage 6 (S6.2b): device rebuild of the per-atom D4 reference weights W/dWq
 *  from the SCF charges (k_d4_build_refw). beginDispersionWeights uploads the
 *  q-independent reference tables once per geometry; dispersionBuildRefW (test
 *  entry) uploads a frozen q and downloads W/dWq. The device-driven loop (S6.5)
 *  launches k_d4_build_refw directly on the resident charges + dD4W/dD4dWq.
 * ====================================================================== */
bool XtbGpuContext::beginDispersionWeights(int nat, const double* cn, const double* gi,
                                           const double* zeff, const double* refcn,
                                           const double* refcovcn, const double* refq,
                                           const int* nref)
{
    if (!ok() || nat <= 0 || !cn || !gi || !zeff || !refcn || !refcovcn || !refq || !nref)
        return false;
    cudaStream_t stream = m_impl->stream;
    const int MR = Impl::D4_MAX_REF;
    try {
        m_impl->dD4Cn.ensure(nat);
        m_impl->dD4Gi.ensure(nat);
        m_impl->dD4Zeff.ensure(nat);
        m_impl->dD4Nref.ensure(nat);
        m_impl->dD4Refcn.ensure(nat * MR);
        m_impl->dD4Refcovcn.ensure(nat * MR);
        m_impl->dD4Refq.ensure(nat * MR);
        m_impl->dD4W.ensure(nat * MR);
        m_impl->dD4dWq.ensure(nat * MR);
    } catch (...) {
        return false;
    }
    m_impl->dD4Cn.upload(cn, nat, stream);
    m_impl->dD4Gi.upload(gi, nat, stream);
    m_impl->dD4Zeff.upload(zeff, nat, stream);
    m_impl->dD4Nref.upload(nref, nat, stream);
    m_impl->dD4Refcn.upload(refcn, nat * MR, stream);
    m_impl->dD4Refcovcn.upload(refcovcn, nat * MR, stream);
    m_impl->dD4Refq.upload(refq, nat * MR, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::dispersionBuildRefW(int nat, const double* q, double* W_out, double* dWq_out)
{
    if (!ok() || nat <= 0 || !q || m_impl->dD4Cn.empty()) return false;
    cudaStream_t stream = m_impl->stream;
    const int MR = Impl::D4_MAX_REF;
    m_impl->dQat.ensure(nat);
    m_impl->dQat.upload(q, nat, stream);
    const int b = 128;
    k_d4_build_refw<<<(nat + b - 1) / b, b, 0, stream>>>(
        nat, MR, m_impl->dQat.ptr, m_impl->dD4Cn.ptr, m_impl->dD4Gi.ptr,
        m_impl->dD4Zeff.ptr, m_impl->dD4Nref.ptr, m_impl->dD4Refcn.ptr,
        m_impl->dD4Refcovcn.ptr, m_impl->dD4Refq.ptr,
        m_impl->dD4W.ptr, m_impl->dD4dWq.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    if (W_out)   m_impl->dD4W.download(W_out, nat * MR, stream);
    if (dWq_out) m_impl->dD4dWq.download(dWq_out, nat * MR, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Stage 6 (S6.3): device SCC energy from the resident OUTPUT charges/moments.
 *  E_coulomb = ½ q_shᵀ γ q_sh (cuBLAS gemv+dot on the resident dGamma/dQsh);
 *  E_third = Σ q_sh³·Γ_s/3 (GFN2 shell, dGamma3) and E_multipole (SD/DD/SQ +
 *  on-site) via block reductions. The band energy stays Σ P⊙H0 (residentDensity).
 *  Requires the resident dQsh (residentShellCharges) + dQat (residentAtomicCharges)
 *  + dDpAt/dQpAt (residentMultipoleMoments) + dGamma/dGamma3/amat (Stage 3/5).
 * ====================================================================== */
bool XtbGpuContext::sccEnergy(int nat, int nsh, double* e_coulomb, double* e_third,
                              double* e_multipole)
{
    if (!ok() || nat <= 0 || nsh <= 0 || m_impl->dGamma.empty() || m_impl->dQsh.empty())
        return false;
    cudaStream_t stream = m_impl->stream;
    const double one = 1.0, zero = 0.0;

    // E_coulomb = ½ q_shᵀ (γ q_sh).
    m_impl->dEScratch.ensure(nsh);
    if (cublasDgemv(m_impl->cublas, CUBLAS_OP_N, nsh, nsh, &one, m_impl->dGamma.ptr, nsh,
                    m_impl->dQsh.ptr, 1, &zero, m_impl->dEScratch.ptr, 1) != CUBLAS_STATUS_SUCCESS)
        return false;
    double ec = 0.0;
    if (cublasDdot(m_impl->cublas, nsh, m_impl->dQsh.ptr, 1, m_impl->dEScratch.ptr, 1, &ec)
        != CUBLAS_STATUS_SUCCESS)
        return false;

    // E_third (GFN2 shell) + E_multipole via reductions into dESca[0..1].
    m_impl->dESca.ensure(2);
    if (cudaMemsetAsync(m_impl->dESca.ptr, 0, sizeof(double) * 2, stream) != cudaSuccess)
        return false;
    const int block = 256;
    if (m_impl->dGamma3.n >= nsh && e_third) {
        const int grid = (nsh + block - 1) / block;
        k_energy_third_order_shell<<<grid, block, block * sizeof(double), stream>>>(
            nsh, m_impl->dQsh.ptr, m_impl->dGamma3.ptr, m_impl->dESca.ptr + 0);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    if (m_impl->mp_otf && e_multipole) {
        if (m_impl->thread_per_atom) {
            const int grid = (nat + block - 1) / block;
            k_energy_multipole_otf<<<grid, block, block * sizeof(double), stream>>>(
                nat, m_impl->dMpXyz.ptr, m_impl->dMpRad.ptr, m_impl->mp_dmp3, m_impl->mp_dmp5,
                m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dDpAt.ptr, m_impl->dQpAt.ptr,
                m_impl->dQat.ptr, m_impl->dESca.ptr + 1);
        } else {
            const int wb = 128;                              // 4 atoms (warps) per block
            const long long threads = 32LL * nat;
            k_energy_multipole_otf_w<<<static_cast<unsigned>((threads + wb - 1) / wb), wb,
                                       (wb / 32) * sizeof(double), stream>>>(
                nat, m_impl->dMpXyz.ptr, m_impl->dMpRad.ptr, m_impl->mp_dmp3, m_impl->mp_dmp5,
                m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dDpAt.ptr, m_impl->dQpAt.ptr,
                m_impl->dQat.ptr, m_impl->dESca.ptr + 1);
        }
        if (cudaGetLastError() != cudaSuccess) return false;
    } else if (!m_impl->dMpAmatSD.empty() && e_multipole) {
        const int grid = (nat + block - 1) / block;
        k_energy_multipole<<<grid, block, block * sizeof(double), stream>>>(
            nat, m_impl->dMpAmatSD.ptr, m_impl->dMpAmatDD.ptr, m_impl->dMpAmatSQ.ptr,
            m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dDpAt.ptr, m_impl->dQpAt.ptr,
            m_impl->dQat.ptr, m_impl->dESca.ptr + 1);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    double esca[2] = {0.0, 0.0};
    m_impl->dESca.download(esca, 2, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
    if (e_coulomb)   *e_coulomb   = 0.5 * ec;
    if (e_third)     *e_third     = esca[0];
    if (e_multipole) *e_multipole = esca[1];
    return true;
}

/* ====================================================================== *
 *  Stage 6 (S6.4): device Broyden mixer (port of BroydenMixer::update). The
 *  vector ops + the M dot products (Gram via cuBLAS gemv/gemm) + the tiny M×M
 *  regularised solve all run on the device; the history + vin_last/F_last stay
 *  resident, so the mixed SCC vector never leaves the GPU in the fused loop.
 * ====================================================================== */
bool XtbGpuContext::broydenBegin(int N, double alpha, int max_hist, double w0)
{
    if (!ok() || N <= 0 || max_hist <= 0 || max_hist > 20) return false;
    try {
        m_impl->dBroyVin.ensure(N); m_impl->dBroyVout.ensure(N); m_impl->dBroyVnext.ensure(N);
        m_impl->dBroyF.ensure(N); m_impl->dBroyFLast.ensure(N);
        m_impl->dBroyVinLast.ensure(N); m_impl->dBroyDFtmp.ensure(N);
        m_impl->dBroyDFmat.ensure(N * max_hist); m_impl->dBroyUmat.ensure(N * max_hist);
        m_impl->dBroyGram.ensure(max_hist * max_hist);
        m_impl->dBroyC.ensure(max_hist); m_impl->dBroyGamma.ensure(max_hist);
    } catch (...) {
        return false;
    }
    m_impl->broyden_N = N; m_impl->broyden_iter = 0; m_impl->broyden_push = 0;
    m_impl->broyden_maxhist = max_hist; m_impl->broyden_alpha = alpha; m_impl->broyden_w0 = w0;
    return true;
}

// Device-pointer core (resident loop + the test wrapper). dvin/dvout/dvnext are
// device pointers of length broyden_N. Queues all work on the stream; the only
// host sync is the cuBLAS nrm2 (the norm<1e-14 branch decision). Claude Generated.
bool XtbGpuContext::runBroydenUpdate(const double* dvin, const double* dvout, double* dvnext)
{
    if (m_impl->broyden_N <= 0) return false;
    cudaStream_t stream = m_impl->stream;
    const int N = m_impl->broyden_N;
    const int maxh = m_impl->broyden_maxhist;
    const double alpha = m_impl->broyden_alpha;
    const size_t bytesN = sizeof(double) * static_cast<size_t>(N);
    const int b1 = 256;
    const int grid = (N + b1 - 1) / b1;

    // F = vout − vin.
    k_vec_sub<<<grid, b1, 0, stream>>>(m_impl->dBroyF.ptr, dvout, dvin, N);
    if (cudaGetLastError() != cudaSuccess) return false;
    ++m_impl->broyden_iter;

    auto linear_step = [&]() -> bool {
        // vnext = vin + alpha·F; store vin_last, F_last.
        if (cudaMemcpyAsync(dvnext, dvin, bytesN, cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
            return false;
        if (cublasDaxpy(m_impl->cublas, N, &alpha, m_impl->dBroyF.ptr, 1, dvnext, 1) != CUBLAS_STATUS_SUCCESS)
            return false;
        cudaMemcpyAsync(m_impl->dBroyVinLast.ptr, dvin, bytesN, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(m_impl->dBroyFLast.ptr, m_impl->dBroyF.ptr, bytesN, cudaMemcpyDeviceToDevice, stream);
        return true;
    };

    if (m_impl->broyden_iter == 1)
        return linear_step();

    // dFraw = F − F_last; norm = ‖dFraw‖.
    k_vec_sub<<<grid, b1, 0, stream>>>(m_impl->dBroyDFtmp.ptr, m_impl->dBroyF.ptr, m_impl->dBroyFLast.ptr, N);
    if (cudaGetLastError() != cudaSuccess) return false;
    double norm = 0.0;
    if (cublasDnrm2(m_impl->cublas, N, m_impl->dBroyDFtmp.ptr, 1, &norm) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (norm < 1.0e-14)
        return linear_step();

    const double inv_norm = 1.0 / norm;
    const int slot = m_impl->broyden_push % maxh;
    k_broyden_dfu<<<grid, b1, 0, stream>>>(
        m_impl->dBroyDFmat.ptr + static_cast<size_t>(slot) * N,
        m_impl->dBroyUmat.ptr + static_cast<size_t>(slot) * N,
        m_impl->dBroyDFtmp.ptr, dvin, m_impl->dBroyVinLast.ptr, alpha, inv_norm, N);
    if (cudaGetLastError() != cudaSuccess) return false;
    ++m_impl->broyden_push;
    const int M = (m_impl->broyden_push < maxh) ? m_impl->broyden_push : maxh;

    // c = DFmatᵀ·F  (length M);  a = DFmatᵀ·DFmat  (M×M, col-major).
    const double one = 1.0, zero = 0.0, neg = -1.0;
    if (cublasDgemv(m_impl->cublas, CUBLAS_OP_T, N, M, &one, m_impl->dBroyDFmat.ptr, N,
                    m_impl->dBroyF.ptr, 1, &zero, m_impl->dBroyC.ptr, 1) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cublasDgemm(m_impl->cublas, CUBLAS_OP_T, CUBLAS_OP_N, M, M, N, &one,
                    m_impl->dBroyDFmat.ptr, N, m_impl->dBroyDFmat.ptr, N, &zero,
                    m_impl->dBroyGram.ptr, M) != CUBLAS_STATUS_SUCCESS)
        return false;
    // (w0²I + a)·gamma = c  (one thread, M ≤ 20).
    const double reg = m_impl->broyden_w0 * m_impl->broyden_w0;
    k_broyden_solve<<<1, 1, 0, stream>>>(M, m_impl->dBroyGram.ptr, m_impl->dBroyC.ptr, reg,
                                         m_impl->dBroyGamma.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    // vnext = vin + alpha·F − Umat·gamma.
    if (cudaMemcpyAsync(dvnext, dvin, bytesN, cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return false;
    if (cublasDaxpy(m_impl->cublas, N, &alpha, m_impl->dBroyF.ptr, 1, dvnext, 1) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cublasDgemv(m_impl->cublas, CUBLAS_OP_N, N, M, &neg, m_impl->dBroyUmat.ptr, N,
                    m_impl->dBroyGamma.ptr, 1, &one, dvnext, 1) != CUBLAS_STATUS_SUCCESS)
        return false;
    // store vin_last = vin, F_last = F.
    cudaMemcpyAsync(m_impl->dBroyVinLast.ptr, dvin, bytesN, cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(m_impl->dBroyFLast.ptr, m_impl->dBroyF.ptr, bytesN, cudaMemcpyDeviceToDevice, stream);
    return true;
}

// Component-test entry: upload vin/vout, run the device update, download vnext.
bool XtbGpuContext::broydenUpdate(int N, const double* vin, const double* vout, double* vnext)
{
    if (!ok() || N != m_impl->broyden_N || !vin || !vout || !vnext) return false;
    cudaStream_t stream = m_impl->stream;
    m_impl->dBroyVin.upload(vin, N, stream);
    m_impl->dBroyVout.upload(vout, N, stream);
    if (!runBroydenUpdate(m_impl->dBroyVin.ptr, m_impl->dBroyVout.ptr, m_impl->dBroyVnext.ptr))
        return false;
    m_impl->dBroyVnext.download(vnext, N, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Stage 6 (S6.5): fully device-resident GFN2 SCF loop (host polls only O(1)).
 *  beginResidentLoop uploads the initial SCC guess + the EEQ-guess q_at once and
 *  sets up the device Broyden; residentScfStep runs ONE fused iteration —
 *    D4 refw(q_at) → potential build + Fock + eigensolve → occupation → density →
 *    q_sh/q_at/moments → SCC energy → Broyden mix → unpack (next input) —
 *  entirely on the device, returning only dq + the 4 energy scalars. The eps /
 *  occ / pop_ao / moments / q_sh never cross the bus. residentLoopCharges
 *  downloads the converged charges into the host wavefunction once at the end, so
 *  the existing post-SCF energy + gradient path runs unchanged. Claude Generated.
 * ====================================================================== */
bool XtbGpuContext::beginResidentLoop(int nsh, int nat, int nao, double Tele, double n_elec,
                                      int nocc_pairs, const double* q_sh0, const double* dp_at0,
                                      const double* qp_at0, const double* q_at0,
                                      const double* n0_sh, const double* n0_at,
                                      double alpha, int max_hist, double w0)
{
    if (!ok() || nsh <= 0 || nat <= 0 || nao <= 0 || nao != m_impl->resident_n
        || m_impl->pot_nsh != nsh || m_impl->resident_nat != nat
        || !q_sh0 || !dp_at0 || !qp_at0 || !q_at0 || !n0_sh || !n0_at)
        return false;
    cudaStream_t stream = m_impl->stream;
    try {
        m_impl->dQat.ensure(nat); m_impl->dQsh.ensure(nsh);
        m_impl->dN0sh.ensure(nsh); m_impl->dN0at.ensure(nat);
        m_impl->dDq.ensure(1);
        m_impl->dOccMu.ensure(1); m_impl->dOccNcol.ensure(1);   // occupied-column count (density)
    } catch (...) {
        return false;
    }
    // Initial mixed SCC input + the EEQ-guess q_at (drives the first D4 weights) +
    // the reference shell/atom occupations (scatter seeds for the q_sh/q_at output).
    m_impl->dPotQsh.upload(q_sh0, nsh, stream);
    m_impl->dInDpAt.upload(dp_at0, 3 * nat, stream);
    m_impl->dInQpAt.upload(qp_at0, 6 * nat, stream);
    m_impl->dQat.upload(q_at0, nat, stream);
    m_impl->dN0sh.upload(n0_sh, nsh, stream);
    m_impl->dN0at.upload(n0_at, nat, stream);
    // Device Broyden over the packed [q_sh; dp_at; qp_at] vector.
    if (!broydenBegin(nsh + 9 * nat, alpha, max_hist, w0)) return false;
    m_impl->loop_nsh = nsh; m_impl->loop_nat = nat; m_impl->loop_nao = nao;
    m_impl->loop_Tele = Tele; m_impl->loop_nelec = n_elec; m_impl->loop_nocc_pairs = nocc_pairs;
    m_impl->last_dq = 1e300;             // a new SCF starts without usable previous eigenvectors
    m_impl->cprev_valid = m_impl->cprev64_valid = false;
    m_impl->pseudo_consecutive = m_impl->pseudo_consecutive64 = 0;
    m_impl->pseudo_taken = m_impl->pseudo_rejected = 0;
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::residentScfStep(bool fp32, double* dq_out, double* e_band,
                                    double* e_coulomb, double* e_third, double* e_multipole)
{
    if (!ok() || m_impl->loop_nao <= 0 || !dq_out || !e_band) return false;
    // Claude Generated (Sep 2026): every failure below names itself - the host used to report
    // "GPU resident SCF step failed at iteration N" with no reason. Nested calls may set a more
    // specific message first; fail() only fills an empty one.
    m_impl->last_error.clear();
    auto fail = [&](const char* what) {
        if (m_impl->last_error.empty()) m_impl->last_error = what;
        return false;
    };
    cudaStream_t stream = m_impl->stream;
    const int nsh = m_impl->loop_nsh, nat = m_impl->loop_nat, nao = m_impl->loop_nao;
    const int b = 128;
    m_impl->profStart();

    // 1. D4 reference weights from the current (previous-output) resident q_at.
    k_d4_build_refw<<<(nat + b - 1) / b, b, 0, stream>>>(
        nat, Impl::D4_MAX_REF, m_impl->dQat.ptr, m_impl->dD4Cn.ptr, m_impl->dD4Gi.ptr,
        m_impl->dD4Zeff.ptr, m_impl->dD4Nref.ptr, m_impl->dD4Refcn.ptr,
        m_impl->dD4Refcovcn.ptr, m_impl->dD4Refq.ptr, m_impl->dD4W.ptr, m_impl->dD4dWq.ptr);
    if (cudaGetLastError() != cudaSuccess) return fail("D4 reference-weight kernel failed");
    m_impl->profMark("scf: D4 reference weights");

    // 2. Potential build (overwrites dQat with the input-derived q_at) + Fock +
    //    eigensolve; eps stays resident (no download).
    m_impl->in_resident_step = true;
    const bool solved = buildDevicePotentialAndSolve(nao, fp32, /*n_eig=*/0, /*eps_out=*/nullptr,
                                                     /*download_eps=*/false);
    m_impl->in_resident_step = false;
    if (!solved) return fail("potential build / Fock / eigensolve failed");

    // 3. Occupation on the device (resident eps → resident occ).
    const double kT = m_impl->loop_Tele * 3.166808e-6;
    const int use_fermi = (m_impl->loop_Tele > 0.0) ? 1 : 0;
    const int blk = 256;
    k_occupations<<<1, blk, blk * sizeof(double), stream>>>(
        m_impl->dEps.ptr, m_impl->dOcc.ptr, nao, kT, m_impl->loop_nelec,
        m_impl->loop_nocc_pairs, use_fermi, m_impl->dOccMu.ptr, m_impl->dOccNcol.ptr);
    if (cudaGetLastError() != cudaSuccess) return fail("occupation kernel failed");
    // Claude Generated (Sep 2026): build P from the occupied columns only (last column with
    // occ > 1e-12, reported by the kernel), exactly like the CPU density (xtb_scf.cpp,
    // leftCols(ncol)). The dropped columns weigh < 1e-12; the GEMM shrinks from nao^3 to
    // nao^2 * ncol (polymer: 26 % of the SCF time was this product over all nao columns).
    int ncol = nao;
    m_impl->dOccNcol.download(&ncol, 1, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess) return fail("stream synchronisation after occupations failed");
    if (ncol <= 0 || ncol > nao) ncol = nao;
    m_impl->profMark("scf: occupations");

    // 4. Density + Mulliken-AO over the occupied columns.
    if (!residentDensityResident(nao, ncol, e_band)) return fail("density / Mulliken populations failed");
    m_impl->profMark("scf: density P + populations");

    // 5. Output charges/moments from the resident density.
    // q_sh = n0_sh − Σ_{μ∈s} pop_ao; q_at = n0_at − Σ_{μ∈A} pop_ao.
    if (cudaMemcpyAsync(m_impl->dQsh.ptr, m_impl->dN0sh.ptr, sizeof(double) * nsh,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess) return fail("shell-charge copy failed");
    k_qsh_scatter<<<(nao + b - 1) / b, b, 0, stream>>>(nao, m_impl->dPop.ptr,
                                                       m_impl->dAo2sh.ptr, m_impl->dQsh.ptr);
    if (cudaMemcpyAsync(m_impl->dQat.ptr, m_impl->dN0at.ptr, sizeof(double) * nat,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess) return fail("atomic-charge copy failed");
    k_qat_scatter<<<(nao + b - 1) / b, b, 0, stream>>>(nao, m_impl->dPop.ptr,
                                                       m_impl->dAo2at.ptr, m_impl->dQat.ptr);
    // Atomic multipole moments dp_at/qp_at from the resident density (resident).
    if (!multipoleMomentsResident(nao, nat)) return fail("atomic multipole moments failed");
    m_impl->profMark("scf: charges + multipole moments");

    // 6. SCC energy components (resident charges/moments).
    if (!sccEnergy(nat, nsh, e_coulomb, e_third, e_multipole)) return fail("SCC energy failed");
    m_impl->profMark("scf: SCC energy");

    // 7. Convergence dq = max|q_sh_out − q_sh_in| (q_sh_in = current dPotQsh).
    k_maxabsdiff<<<1, blk, blk * sizeof(double), stream>>>(m_impl->dQsh.ptr, m_impl->dPotQsh.ptr,
                                                           nsh, m_impl->dDq.ptr);
    if (cudaGetLastError() != cudaSuccess) return fail("convergence (max|dq|) kernel failed");
    m_impl->dDq.download(dq_out, 1, stream);

    // 8. Broyden: pack input [dPotQsh; dInDpAt; dInQpAt] + output [dQsh; dDpAt; dQpAt],
    //    mix, and unpack the next input back into the resident dPotQsh/dInDpAt/dInQpAt.
    {
        double* vin = m_impl->dBroyVin.ptr;
        double* vout = m_impl->dBroyVout.ptr;
        const size_t off_dp = static_cast<size_t>(nsh);
        const size_t off_qp = static_cast<size_t>(nsh) + 3 * nat;
        cudaMemcpyAsync(vin, m_impl->dPotQsh.ptr, sizeof(double) * nsh, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(vin + off_dp, m_impl->dInDpAt.ptr, sizeof(double) * 3 * nat, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(vin + off_qp, m_impl->dInQpAt.ptr, sizeof(double) * 6 * nat, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(vout, m_impl->dQsh.ptr, sizeof(double) * nsh, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(vout + off_dp, m_impl->dDpAt.ptr, sizeof(double) * 3 * nat, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(vout + off_qp, m_impl->dQpAt.ptr, sizeof(double) * 6 * nat, cudaMemcpyDeviceToDevice, stream);
        if (!runBroydenUpdate(vin, vout, m_impl->dBroyVnext.ptr)) return fail("device Broyden update failed");
        cudaMemcpyAsync(m_impl->dPotQsh.ptr, m_impl->dBroyVnext.ptr, sizeof(double) * nsh, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(m_impl->dInDpAt.ptr, m_impl->dBroyVnext.ptr + off_dp, sizeof(double) * 3 * nat, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(m_impl->dInQpAt.ptr, m_impl->dBroyVnext.ptr + off_qp, sizeof(double) * 6 * nat, cudaMemcpyDeviceToDevice, stream);
    }
    m_impl->profMark("scf: dq + Broyden");
    if (cudaStreamSynchronize(stream) != cudaSuccess) return fail("stream synchronisation at end of step failed");
    m_impl->last_dq = *dq_out;   // gates the next step's pseudo-diagonalisation
    return true;
}

bool XtbGpuContext::residentLoopCharges(double* q_sh, double* q_at, double* dp_at,
                                        double* qp_at, double* eps)
{
    if (!ok() || m_impl->loop_nao <= 0) return false;
    m_impl->releaseEigenWorkspaces();   // the resident loop is over (called once after it)
    cudaStream_t stream = m_impl->stream;
    const int nsh = m_impl->loop_nsh, nat = m_impl->loop_nat, nao = m_impl->loop_nao;
    if (q_sh)  m_impl->dQsh.download(q_sh, nsh, stream);
    if (q_at)  m_impl->dQat.download(q_at, nat, stream);
    if (dp_at) m_impl->dDpAt.download(dp_at, 3 * nat, stream);
    if (qp_at) m_impl->dQpAt.download(qp_at, 6 * nat, stream);
    if (eps)   m_impl->dEps.download(eps, nao, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Stage 5 (Part B3/B4): full device GFN2 potential build + resident solve.
 *  beginPotential uploads the geometry-fixed multipole interaction matrices +
 *  the per-shell third-order hardness once per geometry. residentSolvePotential
 *  builds v_sh (γ·q_sh + third-order) + the multipole v_dp/v_qp/v_at scalar shift
 *  + the resident D4 dE/dq on the device, expands v_ao, and folds into the same
 *  Fock build + eigensolve as residentSolveMultipole — so the host SCF loop
 *  uploads only q_sh/dp_at/qp_at (+ the host-built D4 reference weights).
 * ====================================================================== */
bool XtbGpuContext::beginPotential(int nat, int nsh,
                                   const double* amat_sd, const double* amat_dd,
                                   const double* amat_sq, const double* dkernel,
                                   const double* qkernel, const double* gamma3)
{
    if (!ok() || nat <= 0 || nsh <= 0 || !amat_sd || !amat_dd || !amat_sq
        || !dkernel || !qkernel || !gamma3)
        return false;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
    // Claude Generated (Sep 2026): the STORED multipole path holds 18 nat^2 doubles on the
    // device - 7.7 GB at 7320 atoms - on top of everything the SCF already has resident.
    // Without this check the allocation (or a later one starved by it) failed deep inside the
    // SCF and the user saw only "GPU resident SCF step failed at iteration 0". Measured on
    // polymer_2x with `-gpu_multipole_otf off` on a 20 GB A4500.
    {
        const size_t need = 18 * nn * sizeof(double);
        size_t free_b = 0, total_b = 0;
        if (m_impl->memory_check && cudaMemGetInfo(&free_b, &total_b) == cudaSuccess
            && need + (64ull << 20) > free_b) {
            char msg[420];
            std::snprintf(msg, sizeof(msg),
                          "the stored GFN2 multipole interaction matrices need %.1f GB on device %d "
                          "for nat=%d, but only %.1f of %.1f GB are free - rebuild them per "
                          "iteration instead (-gpu_multipole_otf auto, the default above ~2700 atoms)",
                          need / 1073741824.0, m_impl->device, nat,
                          free_b / 1073741824.0, total_b / 1073741824.0);
            m_impl->last_error = msg;
            return false;
        }
    }
    try {
        m_impl->dMpAmatSD.ensure(static_cast<int>(3 * nn));
        m_impl->dMpAmatDD.ensure(static_cast<int>(9 * nn));
        m_impl->dMpAmatSQ.ensure(static_cast<int>(6 * nn));
        m_impl->dMpDkernel.ensure(nat);
        m_impl->dMpQkernel.ensure(nat);
        m_impl->dGamma3.ensure(nsh);
        m_impl->dPotQsh.ensure(nsh);
        m_impl->dInDpAt.ensure(3 * nat);
        m_impl->dInQpAt.ensure(6 * nat);
        m_impl->dVsh.ensure(nsh);
        m_impl->dVat.ensure(nat);
        m_impl->dQat.ensure(nat);
    } catch (...) {
        m_impl->last_error = "device allocation for the stored GFN2 multipole interaction "
                             "matrices failed (try -gpu_multipole_otf auto)";
        return false;
    }
    m_impl->mp_otf = false;
    m_impl->dMpXyz.free();
    m_impl->dMpRad.free();
    m_impl->dMpAmatSD.upload(amat_sd, static_cast<int>(3 * nn), stream);
    m_impl->dMpAmatDD.upload(amat_dd, static_cast<int>(9 * nn), stream);
    m_impl->dMpAmatSQ.upload(amat_sq, static_cast<int>(6 * nn), stream);
    m_impl->dMpDkernel.upload(dkernel, nat, stream);
    m_impl->dMpQkernel.upload(qkernel, nat, stream);
    m_impl->dGamma3.upload(gamma3, nsh, stream);
    m_impl->pot_nsh = nsh;
    // WP4b: solvation is opt-in per geometry via beginSolvation() (called after this);
    // reset here so a solvent-free geometry never re-uses a stale Born matrix.
    m_impl->solv_active = false;
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

// Claude Generated (Sep 2026): device potential without stored interaction matrices. Uploads
// the geometry and the CN-dependent damping radii (nat + 3 nat doubles instead of 18 nat^2);
// k_multipole_potential_otf / k_energy_multipole_otf rebuild the matrix elements per iteration.
bool XtbGpuContext::beginPotentialOnTheFly(int nat, int nsh, const double* xyz_bohr,
                                           const double* mrad, double dmp3, double dmp5,
                                           const double* dkernel, const double* qkernel,
                                           const double* gamma3)
{
    if (!ok() || nat <= 0 || nsh <= 0 || !xyz_bohr || !mrad || !dkernel || !qkernel || !gamma3)
        return false;
    cudaStream_t stream = m_impl->stream;
    try {
        m_impl->dMpAmatSD.free();
        m_impl->dMpAmatDD.free();
        m_impl->dMpAmatSQ.free();
        m_impl->dMpXyz.ensure(3 * nat);
        m_impl->dMpRad.ensure(nat);
        m_impl->dMpDkernel.ensure(nat);
        m_impl->dMpQkernel.ensure(nat);
        m_impl->dGamma3.ensure(nsh);
        m_impl->dPotQsh.ensure(nsh);
        m_impl->dInDpAt.ensure(3 * nat);
        m_impl->dInQpAt.ensure(6 * nat);
        m_impl->dVsh.ensure(nsh);
        m_impl->dVat.ensure(nat);
        m_impl->dQat.ensure(nat);
    } catch (...) {
        return false;
    }
    m_impl->dMpXyz.upload(xyz_bohr, 3 * nat, stream);
    m_impl->dMpRad.upload(mrad, nat, stream);
    m_impl->dMpDkernel.upload(dkernel, nat, stream);
    m_impl->dMpQkernel.upload(qkernel, nat, stream);
    m_impl->dGamma3.upload(gamma3, nsh, stream);
    m_impl->mp_dmp3 = dmp3;
    m_impl->mp_dmp5 = dmp5;
    m_impl->mp_otf = true;
    m_impl->pot_nsh = nsh;
    m_impl->solv_active = false;
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

// WP4b (Claude Generated June 2026): upload the nat×nat Born interaction matrix B so
// the device potential build adds the in-SCF reaction field v_at += B·q_at. Call once
// per geometry after beginPotential() (GFN2 Mulliken path only). born_mat is the
// symmetric, keps-scaled m_born_mat (column-major == row-major for a symmetric matrix).
bool XtbGpuContext::beginSolvation(int nat, const double* born_mat)
{
    if (!ok() || nat <= 0 || !born_mat) return false;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
    try {
        m_impl->dSolvB.ensure(static_cast<int>(nn));
    } catch (...) {
        return false;
    }
    m_impl->dSolvB.upload(born_mat, static_cast<int>(nn), stream);
    m_impl->solv_active = true;
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::residentSolvePotential(const double* q_sh, const double* dp_at,
                                           const double* qp_at, const double* W,
                                           const double* dWq, int n, double* eps_out,
                                           bool fp32, int n_eig)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || m_impl->resident_nat <= 0
        || m_impl->pot_nsh <= 0 || m_impl->dD4C6Flat.empty() || m_impl->dGamma.empty()
        || !q_sh || !dp_at || !qp_at || !W || !dWq || !eps_out)
        return false;
    cudaStream_t stream = m_impl->stream;
    const int nat = m_impl->resident_nat;
    const int nsh = m_impl->pot_nsh;

    // Upload the iteration's mixed SCC quantities + the host-built D4 weights.
    m_impl->dPotQsh.upload(q_sh, nsh, stream);
    m_impl->dInDpAt.upload(dp_at, 3 * nat, stream);
    m_impl->dInQpAt.upload(qp_at, 6 * nat, stream);
    m_impl->dD4W.upload(W, nat * Impl::D4_MAX_REF, stream);
    m_impl->dD4dWq.upload(dWq, nat * Impl::D4_MAX_REF, stream);

    return buildDevicePotentialAndSolve(n, fp32, n_eig, eps_out, /*download_eps=*/true);
}

// Stage 6 core: build the full GFN2 potential from the RESIDENT mixed SCC inputs
// (dPotQsh/dInDpAt/dInQpAt + dD4W/dD4dWq, set by an upload or by the fused step)
// → Fock → eigensolve. Identical to residentSolvePotential's body minus the host
// uploads; download_eps=false keeps the eigenvalues resident for the device loop.
bool XtbGpuContext::buildDevicePotentialAndSolve(int n, bool fp32, int n_eig,
                                                 double* eps_out, bool download_eps)
{
    cudaStream_t stream = m_impl->stream;
    const int nat = m_impl->resident_nat;
    const int nsh = m_impl->pot_nsh;
    const int b = 128;
    // D4 atom-potential dE/dq (B2 kernel; resident, no download).
    {
        const double cut2 = m_impl->d4_cut * m_impl->d4_cut;
        k_d4_dedq<<<(nat + b - 1) / b, b, 0, stream>>>(
            nat, 118, Impl::D4_MAX_REF, m_impl->dD4Z.ptr, m_impl->dD4Nref.ptr,
            m_impl->dD4Sqrt.ptr, m_impl->dD4Xyz.ptr, m_impl->dD4C6Flat.ptr,
            m_impl->dD4W.ptr, m_impl->dD4dWq.ptr,
            m_impl->d4_s6, m_impl->d4_s8, m_impl->d4_a1, m_impl->d4_a2, cut2,
            m_impl->dD4Dedq.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    // q_at = Σ_{s∈A} q_sh(s) (input-derived, for the multipole potential).
    m_impl->dQat.zero(nat, stream);
    k_qsh_to_qat<<<(nsh + b - 1) / b, b, 0, stream>>>(nsh, m_impl->dPotQsh.ptr,
                                                      m_impl->dSh2at.ptr, m_impl->dQat.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    // v_sh = γ·q_sh (cuBLAS gemv on the resident symmetric γ) + shell third-order.
    {
        const double one = 1.0, zero = 0.0;
        if (cublasDgemv(m_impl->cublas, CUBLAS_OP_N, nsh, nsh, &one, m_impl->dGamma.ptr, nsh,
                        m_impl->dPotQsh.ptr, 1, &zero, m_impl->dVsh.ptr, 1) != CUBLAS_STATUS_SUCCESS)
            return false;
        k_vsh_third<<<(nsh + b - 1) / b, b, 0, stream>>>(nsh, m_impl->dPotQsh.ptr,
                                                         m_impl->dGamma3.ptr, m_impl->dVsh.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
    }
    // Multipole potential v_dp/v_qp + v_at scalar shift, then v_at += D4.
    if (m_impl->mp_otf && m_impl->thread_per_atom) {
        k_multipole_potential_otf<<<(nat + b - 1) / b, b, 0, stream>>>(
            nat, m_impl->dMpXyz.ptr, m_impl->dMpRad.ptr, m_impl->mp_dmp3, m_impl->mp_dmp5,
            m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dQat.ptr,
            m_impl->dInDpAt.ptr, m_impl->dInQpAt.ptr,
            m_impl->dVdp.ptr, m_impl->dVqp.ptr, m_impl->dVat.ptr);
    } else if (m_impl->mp_otf) {
        const long long threads = 32LL * nat;
        k_multipole_potential_otf_w<<<static_cast<unsigned>((threads + b - 1) / b), b, 0, stream>>>(
            nat, m_impl->dMpXyz.ptr, m_impl->dMpRad.ptr, m_impl->mp_dmp3, m_impl->mp_dmp5,
            m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dQat.ptr,
            m_impl->dInDpAt.ptr, m_impl->dInQpAt.ptr,
            m_impl->dVdp.ptr, m_impl->dVqp.ptr, m_impl->dVat.ptr);
    } else {
        k_multipole_potential<<<(nat + b - 1) / b, b, 0, stream>>>(
            nat, m_impl->dMpAmatSD.ptr, m_impl->dMpAmatDD.ptr, m_impl->dMpAmatSQ.ptr,
            m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dQat.ptr,
            m_impl->dInDpAt.ptr, m_impl->dInQpAt.ptr,
            m_impl->dVdp.ptr, m_impl->dVqp.ptr, m_impl->dVat.ptr);
    }
    if (cudaGetLastError() != cudaSuccess) return false;
    k_vat_add_d4<<<(nat + b - 1) / b, b, 0, stream>>>(nat, m_impl->dD4Dedq.ptr, m_impl->dVat.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    // WP4b: in-SCF implicit solvation reaction field v_at += B·q_at (GFN2 Mulliken).
    // B is symmetric so column-major (uploaded) vs row-major is identical; beta=1
    // accumulates onto the multipole+D4 v_at. dQat is the input-derived atomic charge.
    if (m_impl->solv_active) {
        const double one = 1.0;
        if (cublasDgemv(m_impl->cublas, CUBLAS_OP_N, nat, nat, &one, m_impl->dSolvB.ptr, nat,
                        m_impl->dQat.ptr, 1, &one, m_impl->dVat.ptr, 1) != CUBLAS_STATUS_SUCCESS)
            return false;
    }
    // Expand to AO: v_ao(μ) = v_sh(ao2sh[μ]) + v_at(ao2at[μ]).
    k_expand_vao<<<(n + b - 1) / b, b, 0, stream>>>(n, m_impl->dVsh.ptr, m_impl->dVat.ptr,
                                                    m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr,
                                                    m_impl->dVao.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    m_impl->profMark("scf: potential (D4 dE/dq, gamma, multipole)");
    // F = H0 − ½·S·(v_ao⊕v_ao) + GFN2 multipole; then eigensolve (shared path).
    if (!buildFockIntoC(n, /*multipole=*/true)) return false;
    m_impl->profMark("scf: Fock build");

    return eigensolveResidentFock(eps_out, fp32, n_eig, download_eps);
}

bool XtbGpuContext::computeOverlapGrad(const double* xyz_bohr, double* dSdR_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !xyz_bohr || !dSdR_out) return false;
    const int nao = m_impl->basis_nao;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
    try {
        if (m_impl->dSdR.n < static_cast<int>(3 * nn)) m_impl->dSdR.alloc(static_cast<int>(3 * nn));
    } catch (...) { return false; }
    m_impl->dXyz.upload(xyz_bohr, 3 * m_impl->basis_nat, stream);

    const dim3 block(16, 16);
    const dim3 grid((nao + block.x - 1) / block.x, (nao + block.y - 1) / block.y);
    k_overlap_grad<<<grid, block, 0, stream>>>(
        nao, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dIaoSh.ptr, m_impl->dAng.ptr,
        m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
        m_impl->dPrimCoeff.ptr, m_impl->dXyz.ptr, m_impl->dSdR.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    m_impl->dSdR.download(dSdR_out, static_cast<int>(3 * nn), stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::computeGradient(const double* P, const double* C, const double* eps,
                                    int nocc_orbs, const double* v_ao, const double* q_sh,
                                    const double* v_dp, const double* v_qp,
                                    double* grad_out, double* dEdcn_out, bool pc_resident)
{
    // AP8 (Claude Generated): pc_resident=true reuses the resident density dP and MO
    // coefficients dC the device-resident SCF already left on the GPU (downloaded
    // once by residentFinalize but not cleared), skipping the two nao²-sized P/C
    // host→device uploads (~2 ms/step on complex). The caller (calculateGradientGpu,
    // only invoked on the device-resident path) guarantees they are the converged
    // values; the gradient only reads dP/dC, so reuse is bit-identical to uploading.
    if (!ok() || m_impl->basis_nao <= 0 || !eps || !v_ao || !q_sh
        || !grad_out || !dEdcn_out || (!pc_resident && (!P || !C)))
        return false;
    const int nat = m_impl->basis_nat;
    const int nsh = m_impl->basis_nsh;
    const int nao = m_impl->basis_nao;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
    const double one = 1.0, zero = 0.0;

    // Points 1+2 (Claude Generated, Sep 2026): on the device-resident sparse path
    // (pc_resident=true, screened pattern in use) P and W are only ever read at stored
    // pattern positions by the kernels below (k_grad_h0_pulay_spP, k_grad_cn_onsite_spP -
    // verified exhaustively, see the *_spP kernels above) - so neither needs the dense
    // nao x nao buffer at all. use_pattern_pw selects that fast path; pc_resident=false
    // (explicit P/C upload) and the small/dense-geometry case keep the original dense
    // buffers unchanged as the reference/fallback.
    const bool use_pattern_pw = pc_resident && m_impl->sparse;

    // Claude Generated (Sep 2026): phase marks for the device gradient, same
    // profMark/CURCUMA_GPU_PROFILE mechanism as the integral build and the resident
    // SCF loop (env-gated; no cost when unset — see profMark()).
    m_impl->profStart();

    try {
        // Claude Generated (Sep 2026): the eigensolver workspaces are not needed by the gradient;
        // release them BEFORE allocating W so the gradient phase does not add to the SCF peak
        // (the allocation used to come first).
        m_impl->releaseEigenWorkspaces();
        // Point 2: dW/dP (dense, nao x nao) are only needed by the dense-buffer paths.
        if (!use_pattern_pw) {
            if (m_impl->dW.n < static_cast<int>(nn)) m_impl->dW.alloc(static_cast<int>(nn));
            if (m_impl->dP.n  < static_cast<int>(nn)) m_impl->dP.alloc(static_cast<int>(nn));
        } else {
            m_impl->dSpW.ensure(m_impl->sp_nnz);
        }
        // dCw only ever holds nao x nocc_orbs (k_scale_cols/k_density_sp write no more than
        // that) - was allocated at the full nao x nao here although residentDensityResident
        // already sizes it correctly; fixed while touching this block (memory only, no
        // measured time effect - flagged by the operator, not a separate optimisation).
        const size_t cw_n = nocc_orbs > 0 ? static_cast<size_t>(nao) * nocc_orbs : 0;
        if (cw_n > 0 && m_impl->dCw.n < static_cast<int>(cw_n)) m_impl->dCw.alloc(static_cast<int>(cw_n));
        if (m_impl->dC.n  < static_cast<int>(nn)) m_impl->dC.alloc(static_cast<int>(nn));
        if (m_impl->dVao.n < nao) m_impl->dVao.alloc(nao);
        if (m_impl->dOcc.n < nao) m_impl->dOcc.alloc(nao);
        m_impl->dQsh.ensure(nsh);
        m_impl->dGrad.ensure(3 * nat);
        m_impl->dEdcn.ensure(nat);
    } catch (...) { return false; }
    m_impl->profMark("grad: workspace alloc");

    // Upload the converged SCF state (P/C symmetric-or-column-major from host).
    // AP8: skip the nao²-sized P/C uploads when they are already resident.
    if (!pc_resident) {
        m_impl->dP.upload(P, static_cast<int>(nn), stream);
        m_impl->dC.upload(C, static_cast<int>(nn), stream);
        m_impl->p_dense_valid = true;
    } else if (!use_pattern_pw) {
        if (!ensureDenseDensity(nao)) return false;   // dense geometry: rebuild if stale
    }
    // use_pattern_pw: nothing to do - dSpP already holds the converged pattern density
    // from the SCF loop's last iteration (built by the SAME k_density_sp call the loop
    // already ran for the population/band energy, from the SAME dCw/dC this call reuses
    // via pc_resident). Point 1+2 together: no dense P, no download, no DGEMM rebuild.
    m_impl->profMark("grad: dense density rebuild (ensureDenseDensity)");
    m_impl->dVao.upload(v_ao, nao, stream);
    m_impl->dQsh.upload(q_sh, nsh, stream);
    // GFN2 multipole potentials (converged) for the multipole-integral Pulay term.
    const bool with_mp = (m_impl->basis_is_gfn2 && v_dp && v_qp);
    if (with_mp) {
        try {
            if (m_impl->dVdp.n < 3 * nat) m_impl->dVdp.alloc(3 * nat);
            if (m_impl->dVqp.n < 6 * nat) m_impl->dVqp.alloc(6 * nat);
        } catch (...) { return false; }
        m_impl->dVdp.upload(v_dp, 3 * nat, stream);
        m_impl->dVqp.upload(v_qp, 6 * nat, stream);
    }
    m_impl->profMark("grad: upload v_ao/q_sh/v_dp/v_qp");

    // Energy-weighted density W = C_occ · diag(2·ε_occ) · C_occᵀ.
    // Point 2 (Claude Generated, Sep 2026): on the pattern path, build W ONLY at the stored
    // screened pairs with the same SDDMM kernel (k_density_sp) the SCF loop already uses for
    // P, weight swapped for 2*eps - no nao x nao DGEMM. polymer_2x: 12.75 s -> measured below.
    // Point 6 (Claude Generated, Sep 2026): W is this same SDDMM, so it can go through
    // densityPatternDistributed exactly like the SCF loop's P does - see its
    // for_weighted_density flag. dOcc already holds 2*eps (uploaded just above) and the
    // full-nao dCw computed here is kept as the fallback the single-device k_density_sp below
    // still needs when the split is not configured or the geometry is too small to split.
    // Claude Generated (Sep 29, 2026): after a pseudo-diagonalisation step the columns of C span the
    // right occupied space but do not diagonalise the occupied block of the Fock matrix, so
    // C diag(2 eps) C^T misses the off-diagonal G_ij (polymer: 8.6e-4 Eh/A after FP64 pseudo steps
    // started from FP32 vectors). The basis-invariant form W = C_o (2 G_oo) C_o^T with
    // G_oo = C_o^T F C_o is exact for integer occupations (the pseudo step's gap gate). F is rebuilt
    // from the potentials uploaded above, i.e. the converged ones.
    const bool w_from_goo = nocc_orbs > 0 && !m_impl->c_exact_eigenvectors;
    if (w_from_goo) {
        const int no = nocc_orbs;
        CudaBuffer<double> Cbak, Goo;
        try {
            Cbak.ensure(static_cast<int>(nn));
            Goo.ensure(static_cast<int>(static_cast<size_t>(no) * no));
            m_impl->dCw.ensure(static_cast<int>(static_cast<size_t>(nao) * no));
        } catch (...) { return false; }
        const double two = 2.0;
        if (cudaMemcpyAsync(Cbak.ptr, m_impl->dC.ptr, sizeof(double) * nn, cudaMemcpyDeviceToDevice, stream)
                != cudaSuccess
            || !buildFockIntoC(nao, with_mp)
            || cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_N, nao, no, nao, &one, m_impl->dC.ptr, nao,
                           Cbak.ptr, nao, &zero, m_impl->dCw.ptr, nao) != CUBLAS_STATUS_SUCCESS          // T = F C_o
            || cublasDgemm(m_impl->cublas, CUBLAS_OP_T, CUBLAS_OP_N, no, no, nao, &one, Cbak.ptr, nao,
                           m_impl->dCw.ptr, nao, &zero, Goo.ptr, no) != CUBLAS_STATUS_SUCCESS            // G_oo
            || cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_N, nao, no, no, &two, Cbak.ptr, nao,
                           Goo.ptr, no, &zero, m_impl->dCw.ptr, nao) != CUBLAS_STATUS_SUCCESS            // C_o 2G_oo
            || cudaMemcpyAsync(m_impl->dC.ptr, Cbak.ptr, sizeof(double) * nn, cudaMemcpyDeviceToDevice, stream)
                != cudaSuccess)
            return false;
        if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
        m_impl->profMark("grad: W from G_oo (after a pseudo-diagonalisation step)");
    }
    if (nocc_orbs > 0) {
        if (!w_from_goo) {
            std::vector<double> occ2(nocc_orbs);
            for (int k = 0; k < nocc_orbs; ++k) occ2[k] = 2.0 * eps[k];
            m_impl->dOcc.upload(occ2.data(), nocc_orbs, stream);
            const dim3 block(16, 16);
            const dim3 grid((nao + block.x - 1) / block.x, (nocc_orbs + block.y - 1) / block.y);
            k_scale_cols<<<grid, block, 0, stream>>>(m_impl->dCw.ptr, m_impl->dC.ptr,
                                                     m_impl->dOcc.ptr, nao, nocc_orbs);
            if (cudaGetLastError() != cudaSuccess) return false;
        }
        if (use_pattern_pw) {
            // The split rebuilds Cw from C and the occupations on each helper, which the G_oo form
            // does not fit; that case stays on this device.
            const bool w_split = !w_from_goo
                && densityPatternDistributed(nao, nocc_orbs, /*for_weighted_density=*/true);
            if (!w_split) {
                const int bs = 256;
                k_density_sp<<<(m_impl->sp_nnz + bs - 1) / bs, bs, 0, stream>>>(
                    m_impl->sp_nnz, m_impl->dSpRow.ptr, m_impl->dSpCol.ptr,
                    m_impl->dCw.ptr, m_impl->dC.ptr, nao, nocc_orbs, m_impl->dSpW.ptr);
                if (cudaGetLastError() != cudaSuccess) return false;
            }
        } else if (cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_T, nao, nao, nocc_orbs,
                        &one, m_impl->dCw.ptr, nao, m_impl->dC.ptr, nao,
                        &zero, m_impl->dW.ptr, nao) != CUBLAS_STATUS_SUCCESS)
            return false;
    } else if (use_pattern_pw) {
        m_impl->dSpW.zero(m_impl->sp_nnz, stream);
    } else {
        if (cudaMemsetAsync(m_impl->dW.ptr, 0, sizeof(double) * nn, stream) != cudaSuccess)
            return false;
    }
    m_impl->profMark("grad: W = C_occ*diag(2eps)*C_occ^T (DGEMM)");

    if (cudaMemsetAsync(m_impl->dGrad.ptr, 0, sizeof(double) * 3 * nat, stream) != cudaSuccess)
        return false;
    if (cudaMemsetAsync(m_impl->dEdcn.ptr, 0, sizeof(double) * nat, stream) != cudaSuccess)
        return false;

    // Repulsion scalars (rep_kexp=1.5, rep_rexp=1.0 both; rep_kexp_light: GFN1 1.5, GFN2 1.0).
    const double kexp = 1.5, rexp = 1.0;
    const double kexp_light = m_impl->basis_is_gfn2 ? 1.0 : 1.5;

    const int b1 = 128;
    k_grad_repulsion<<<(nat + b1 - 1) / b1, b1, 0, stream>>>(
        nat, m_impl->basis_is_gfn2, m_impl->dZ.ptr, m_impl->dXyz.ptr,
        m_impl->dRepAlpha.ptr, m_impl->dRepZeff.ptr, kexp, rexp, kexp_light, m_impl->dGrad.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    // Point 2 (Claude Generated, Sep 2026): CN-onsite reads only P(mu,mu); on the pattern
    // path it comes from dSpP via the diagonal index instead of the dense buffer.
    if (use_pattern_pw) {
        k_grad_cn_onsite_spP<<<(nao + b1 - 1) / b1, b1, 0, stream>>>(
            nao, m_impl->dSpP.ptr, m_impl->dSpDiagIdx.ptr, m_impl->dKcn.ptr,
            m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dEdcn.ptr);
    } else {
        k_grad_cn_onsite<<<(nao + b1 - 1) / b1, b1, 0, stream>>>(
            nao, m_impl->dP.ptr, m_impl->dKcn.ptr, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr,
            m_impl->dEdcn.ptr);
    }
    if (cudaGetLastError() != cudaSuccess) return false;
    m_impl->profMark("grad: repulsion + CN onsite");

    const dim3 block(16, 16);
    const dim3 grid((nao + block.x - 1) / block.x, (nao + block.y - 1) / block.y);
    if (use_pattern_pw) {
        const int bs = 256;
        k_grad_h0_pulay_spP<<<(m_impl->sp_nnz + bs - 1) / bs, bs, 0, stream>>>(
            m_impl->sp_nnz, m_impl->dSpRow.ptr, m_impl->dSpCol.ptr, m_impl->dSpS.ptr, m_impl->dSpH0.ptr,
            m_impl->dSpP.ptr, m_impl->dSpW.ptr,
            m_impl->basis_is_gfn2, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dAng.ptr,
            m_impl->dIaoSh.ptr, m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
            m_impl->dPrimCoeff.ptr, m_impl->dShZeta.ptr, m_impl->dShpoly.ptr, m_impl->dKcn.ptr,
            m_impl->dValence.ptr, m_impl->dZ.ptr, m_impl->dSE.ptr, m_impl->dXyz.ptr,
            m_impl->dVao.ptr,
            with_mp ? m_impl->dVdp.ptr : nullptr, with_mp ? m_impl->dVqp.ptr : nullptr,
            m_impl->dGrad.ptr, m_impl->dEdcn.ptr);
    } else if (m_impl->sparse) {
        const int bs = 256;
        k_grad_h0_pulay_sp<<<(m_impl->sp_nnz + bs - 1) / bs, bs, 0, stream>>>(
            m_impl->sp_nnz, m_impl->dSpRow.ptr, m_impl->dSpCol.ptr, m_impl->dSpS.ptr, m_impl->dSpH0.ptr,
            nao, m_impl->basis_is_gfn2, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dAng.ptr,
            m_impl->dIaoSh.ptr, m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
            m_impl->dPrimCoeff.ptr, m_impl->dShZeta.ptr, m_impl->dShpoly.ptr, m_impl->dKcn.ptr,
            m_impl->dValence.ptr, m_impl->dZ.ptr, m_impl->dSE.ptr, m_impl->dXyz.ptr, m_impl->dP.ptr,
            m_impl->dW.ptr, m_impl->dVao.ptr,
            with_mp ? m_impl->dVdp.ptr : nullptr, with_mp ? m_impl->dVqp.ptr : nullptr,
            m_impl->dGrad.ptr, m_impl->dEdcn.ptr);
    } else {
        k_grad_h0_pulay<<<grid, block, 0, stream>>>(
            nao, m_impl->basis_is_gfn2, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dAng.ptr,
            m_impl->dIaoSh.ptr, m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
            m_impl->dPrimCoeff.ptr, m_impl->dShZeta.ptr, m_impl->dShpoly.ptr, m_impl->dKcn.ptr,
            m_impl->dValence.ptr, m_impl->dZ.ptr, m_impl->dSE.ptr, m_impl->dXyz.ptr, m_impl->dP.ptr,
            m_impl->dS.ptr, m_impl->dH0.ptr, m_impl->dW.ptr, m_impl->dVao.ptr,
            with_mp ? m_impl->dVdp.ptr : nullptr, with_mp ? m_impl->dVqp.ptr : nullptr,
            m_impl->dGrad.ptr, m_impl->dEdcn.ptr);
    }
    if (cudaGetLastError() != cudaSuccess) return false;
    m_impl->profMark(use_pattern_pw ? "grad: H0 Pulay (screened pair kernel, pattern-only P/W)"
                      : m_impl->sparse ? "grad: H0 Pulay (screened pair kernel)"
                                       : "grad: H0 Pulay (dense)");

    k_grad_coulomb<<<(nsh + b1 - 1) / b1, b1, 0, stream>>>(
        nsh, m_impl->basis_is_gfn2, m_impl->dSh2at.ptr, m_impl->dHardness.ptr,
        m_impl->dQsh.ptr, m_impl->dXyz.ptr, 2.0, m_impl->dGrad.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
    m_impl->profMark("grad: Coulomb");

    m_impl->dGrad.download(grad_out, 3 * nat, stream);
    m_impl->dEdcn.download(dEdcn_out, nat, stream);
    const bool ok_sync = cudaStreamSynchronize(stream) == cudaSuccess;
    m_impl->profMark("grad: download + sync");
    return ok_sync;
}

// Stage 6 (S6.1): device occupation. Component-test entry — uploads a frozen eps
// into the resident dEps, runs the single-block k_occupations, downloads occ (+
// µ, ncol). The device-driven loop (S6.5) calls the same kernel on the resident
// dEps with no upload/download. Claude Generated.
bool XtbGpuContext::occupations(const double* eps, int n, double Tele, double n_elec,
                                double* occ_out, int* ncol_out, double* mu_out)
{
    if (!ok() || n <= 0 || !eps || !occ_out) return false;
    cudaStream_t stream = m_impl->stream;
    m_impl->dEps.ensure(n);
    m_impl->dOcc.ensure(n);
    m_impl->dOccMu.ensure(1);
    m_impl->dOccNcol.ensure(1);
    m_impl->dEps.upload(eps, n, stream);

    const double kT = Tele * 3.166808e-6;            // K → Hartree (host constant)
    const int nocc_pairs = static_cast<int>(std::floor(n_elec / 2.0));
    const int use_fermi = (Tele > 0.0) ? 1 : 0;
    const int block = 256;                            // power of two for the tree reduction
    k_occupations<<<1, block, block * sizeof(double), stream>>>(
        m_impl->dEps.ptr, m_impl->dOcc.ptr, n, kT, n_elec, nocc_pairs, use_fermi,
        m_impl->dOccMu.ptr, m_impl->dOccNcol.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    m_impl->dOcc.download(occ_out, n, stream);
    int ncol_h = 0; double mu_h = 0.0;
    m_impl->dOccNcol.download(&ncol_h, 1, stream);
    m_impl->dOccMu.download(&mu_h, 1, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
    if (ncol_out) *ncol_out = ncol_h;
    if (mu_out)   *mu_out   = mu_h;
    return true;
}

bool XtbGpuContext::residentBeginMultipoleComputed()
{
    const bool have_mp = m_impl && (m_impl->sparse ? !m_impl->dSpDp.empty() : !m_impl->dDpInt.empty());
    if (!ok() || m_impl->basis_nat <= 0 || !have_mp) return false;
    const int nat = m_impl->basis_nat;
    try {
        // dp_int/qp_int are already resident (computeIntegrals); dAo2at uploaded in
        // beginBasis. Allocate only the per-iteration potential / moment buffers.
        m_impl->dVdp.ensure(3 * nat);
        m_impl->dVqp.ensure(6 * nat);
        m_impl->dDpAt.ensure(3 * nat);
        m_impl->dQpAt.ensure(6 * nat);
    } catch (...) {
        return false;
    }
    m_impl->resident_nat = nat;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_CUDA
