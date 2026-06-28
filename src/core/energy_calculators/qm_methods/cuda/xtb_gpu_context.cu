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

#ifdef USE_CUDA_XTB

#include "xtb_gpu_context.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <cmath>

// CudaBuffer<T> (RAII cudaMalloc/cudaFree) — the project's established device
// allocation helper; routing every allocation through it is the mitigation for
// the GFN-FF GPU heap-corruption class of bug.
#include "../../ff_methods/cuda/gfnff_soa.h"

// Stage 3 device integral functions (CN counting; later overlap/multipole) and
// the host element-parameter tables that seed __constant__ memory.
#include "xtb_gpu_integrals_device.cuh"
#include "../parameters/xtb_params_extra.hpp"

namespace curcuma {
namespace xtb {
namespace gpu {

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
};

// ---- Stage 3 element parameter tables in __constant__ memory --------------
// Only the D3 covalent radii (au) are needed for the CN kernel; more tables
// (hubbard, shell_hubbard, pauling_en, …) join here in 3b/3c/3d. The 86-element
// arrays are far under the 64 KB constant limit. Seeded once via ensureStage3Constants.
namespace {
__constant__ double c_covrad_d3[86];   // covalent_rad_d3_au(z), index z-1
__constant__ double c_pauling[86];     // pauling_en[z-1]
__constant__ double c_atomic_rad[86];  // atomic_rad_au(z), index z-1

void ensureStage3Constants()
{
    static bool loaded = false;
    if (loaded) return;
    double covrad[86], pauling[86], arad[86];
    for (int z = 1; z <= 86; ++z) {
        covrad[z - 1]  = curcuma::xtb::covalent_rad_d3_au(z);
        pauling[z - 1] = curcuma::xtb::pauling_en[z - 1];
        arad[z - 1]    = curcuma::xtb::atomic_rad_au(z);
    }
    cudaMemcpyToSymbol(c_covrad_d3, covrad, sizeof(covrad));
    cudaMemcpyToSymbol(c_pauling, pauling, sizeof(pauling));
    cudaMemcpyToSymbol(c_atomic_rad, arad, sizeof(arad));
    loaded = true;
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
    double* __restrict__ S, double* __restrict__ H0)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    const int b = blockIdx.y * blockDim.y + threadIdx.y;
    if (a >= nsh || b >= nsh) return;

    const int iat = sh2at[a], jat = sh2at[b];
    const int la = ang_sh[a], lb = ang_sh[b];
    const double xa = xyz[3 * iat + 0], ya = xyz[3 * iat + 1], za = xyz[3 * iat + 2];
    const double xb = xyz[3 * jat + 0], yb = xyz[3 * jat + 1], zb = xyz[3 * jat + 2];
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

// Section 2b: H0/Pulay off-site gradient. One thread per AO pair (μ,ν) with iat<jat.
// Mirrors xtb_gradient.cpp:228-438 (GFN1/GFN2 isotropic part; the GFN2 multipole
// integral Pulay block is Stage 4b and handled separately).
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

    const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
    const double Pmn = P[mn], Smn = S[mn], H0mn = H0[mn], Wmn = W[mn];
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

// Section 3: isotropic Coulomb gradient. One thread per shell is, inner js<is.
__global__ void k_grad_coulomb(int nsh, int is_gfn2, const int* __restrict__ sh2at,
                               const double* __restrict__ g, const double* __restrict__ q_sh,
                               const double* __restrict__ xyz, double gexp, double* __restrict__ grad)
{
    const int is = blockIdx.x * blockDim.x + threadIdx.x;
    if (is >= nsh) return;
    const int iat = sh2at[is];
    for (int js = 0; js < is; ++js) {
        const int jat = sh2at[js];
        if (iat == jat) continue;
        const double dx = xyz[3*iat+0] - xyz[3*jat+0];
        const double dy = xyz[3*iat+1] - xyz[3*jat+1];
        const double dz = xyz[3*iat+2] - xyz[3*jat+2];
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1.0e-12) continue;
        const double r1 = sqrt(r2);
        const double gam_bar = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
        const double gamma = pow(pow(r1, gexp) + pow(gam_bar, -gexp), -1.0 / gexp);
        const double dgamma_dr = -pow(r1, gexp - 2.0) * pow(gamma, gexp + 1.0);
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

XtbGpuContext::XtbGpuContext()
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
        return; // no device — caller falls back to CPU

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
    if (m_impl->cusolver) cusolverDnDestroy(m_impl->cusolver);
    if (m_impl->cublas)   cublasDestroy(m_impl->cublas);
    if (m_impl->stream)   cudaStreamDestroy(m_impl->stream);
}

bool XtbGpuContext::ok() const { return m_impl && m_impl->ok; }

std::string XtbGpuContext::deviceName() const
{
    return m_impl ? m_impl->name : std::string();
}

int XtbGpuContext::deviceId() const { return m_impl ? m_impl->device : -1; }

bool XtbGpuContext::deviceAvailable()
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
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
        m_impl->dH0.upload(H0, static_cast<int>(nn), m_impl->stream);
        m_impl->dS.upload(S, static_cast<int>(nn), m_impl->stream);
        m_impl->dL.upload(L, static_cast<int>(nn), m_impl->stream);
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
bool XtbGpuContext::eigensolveResidentFock(double* eps_out, bool fp32, int n_eig,
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

    if (fp32) {
        // Mixed precision: reduce + diagonalise + back-transform in FP32, then
        // convert the eigenvectors/values back to FP64 so the resident density
        // (FP64) is unchanged. ~5–10× faster than FP64 on consumer GPUs.
        try {
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
        int lwork = 0;
        if (partial) {
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
        }
        // Back-transform only the neig computed eigenvectors (columns 0..neig-1).
        if (cublasStrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                        CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, neig, &onef,
                        m_impl->dLf.ptr, n, m_impl->dCf.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
        // Convert the neig columns (column-major: first neig·n contiguous) + neig eps.
        const size_t conv = static_cast<size_t>(neig) * static_cast<size_t>(n);
        k_f2d<<<static_cast<int>((conv + b - 1) / b), b, 0, stream>>>(m_impl->dCf.ptr, m_impl->dC.ptr, conv);
        k_f2d<<<(neig + b - 1) / b, b, 0, stream>>>(m_impl->dEpsf.ptr, m_impl->dEps.ptr, neig);
        if (cudaGetLastError() != cudaSuccess) return false;

        int info = 1;
        if (cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
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
    if (partial) {
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
    } else {
        if (cusolverDnDsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                             n, m_impl->dC.ptr, n, m_impl->dEps.ptr, m_impl->dWork.ptr,
                             m_impl->lwork, m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
            return false;
    }
    // Back-transform only the neig computed eigenvectors (columns 0..neig-1).
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, neig, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;

    int info = 1;
    if (cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
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
    const dim3 block(16, 16);
    const dim3 grid((n + block.x - 1) / block.x, (n + block.y - 1) / block.y);
    k_build_fock_iso<<<grid, block, 0, stream>>>(m_impl->dC.ptr, m_impl->dH0.ptr,
                                                 m_impl->dS.ptr, m_impl->dVao.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
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

    if (ncol > 0) {
        m_impl->dOcc.upload(occ, ncol, stream);
        // Cw(:,k) = occ[k]·C(:,k)  for k < ncol.
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

    // pop_ao(μ) = Σ_ν P(μ,ν)·S(μ,ν).
    const int b1 = 128;
    k_pop_ao<<<(n + b1 - 1) / b1, b1, 0, stream>>>(m_impl->dPop.ptr, m_impl->dP.ptr,
                                                   m_impl->dS.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
        return false;

    // Band energy = Σ_μν P_μν·H0_μν (host-pointer dot is blocking; sync to be safe).
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    if (cublasDdot(m_impl->cublas, static_cast<int>(nn), m_impl->dP.ptr, 1,
                   m_impl->dH0.ptr, 1, band_out) != CUBLAS_STATUS_SUCCESS)
        return false;

    m_impl->dPop.download(pop_ao_out, n, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

// Stage 6: density + Mulliken-AO from the RESIDENT occupations (dOcc, set by
// k_occupations) — no occ upload, no pop_ao download (dPop stays resident for the
// device q_sh/q_at reductions). band_out = Σ P⊙H0 (host scalar). Claude Generated.
bool XtbGpuContext::residentDensityResident(int n, int ncol, double* band_out)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !band_out) return false;
    cudaStream_t stream = m_impl->stream;
    const double one = 1.0, zero = 0.0;
    if (ncol > 0) {
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (ncol + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, stream>>>(m_impl->dCw.ptr, m_impl->dC.ptr,
                                                 m_impl->dOcc.ptr, n, ncol);
        if (cudaGetLastError() != cudaSuccess) return false;
        if (cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, ncol,
                        &one, m_impl->dCw.ptr, n, m_impl->dC.ptr, n,
                        &zero, m_impl->dP.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
    } else if (cudaMemsetAsync(m_impl->dP.ptr, 0, sizeof(double) * static_cast<size_t>(n) * n,
                               stream) != cudaSuccess) {
        return false;
    }
    const int b1 = 128;
    k_pop_ao<<<(n + b1 - 1) / b1, b1, 0, stream>>>(m_impl->dPop.ptr, m_impl->dP.ptr,
                                                   m_impl->dS.ptr, n);
    if (cudaGetLastError() != cudaSuccess) return false;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    return cublasDdot(m_impl->cublas, static_cast<int>(nn), m_impl->dP.ptr, 1,
                      m_impl->dH0.ptr, 1, band_out) == CUBLAS_STATUS_SUCCESS;
}

bool XtbGpuContext::residentFinalize(double* P_colmajor, double* C_colmajor, int n)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !P_colmajor || !C_colmajor)
        return false;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
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
    const dim3 block(16, 16);
    const dim3 grid((n + block.x - 1) / block.x, (n + block.y - 1) / block.y);
    k_build_fock_iso<<<grid, block, 0, stream>>>(m_impl->dC.ptr, m_impl->dH0.ptr,
                                                 m_impl->dS.ptr, m_impl->dVao.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
        return false;
    k_add_fock_multipole<<<grid, block, 0, stream>>>(m_impl->dC.ptr, m_impl->dDpInt.ptr,
                                                     m_impl->dQpInt.ptr, m_impl->dVdp.ptr,
                                                     m_impl->dVqp.ptr, m_impl->dAo2at.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
        return false;

    return eigensolveResidentFock(eps_out, fp32, n_eig);
}

bool XtbGpuContext::residentMultipoleMoments(double* dp_at3, double* qp_at6, int n, int nat)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || nat != m_impl->resident_nat
        || !dp_at3 || !qp_at6)
        return false;
    cudaStream_t stream = m_impl->stream;

    // Accumulator-scatter kernel needs zeroed targets.
    if (cudaMemsetAsync(m_impl->dDpAt.ptr, 0, sizeof(double) * 3 * nat, stream) != cudaSuccess)
        return false;
    if (cudaMemsetAsync(m_impl->dQpAt.ptr, 0, sizeof(double) * 6 * nat, stream) != cudaSuccess)
        return false;

    const int b1 = 128;
    k_multipole_moments<<<(n + b1 - 1) / b1, b1, 0, stream>>>(
        m_impl->dDpAt.ptr, m_impl->dQpAt.ptr, m_impl->dP.ptr,
        m_impl->dDpInt.ptr, m_impl->dQpInt.ptr, m_impl->dAo2at.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
        return false;

    m_impl->dDpAt.download(dp_at3, 3 * nat, stream);
    m_impl->dQpAt.download(qp_at6, 6 * nat, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

/* ====================================================================== *
 *  Device-side integral build (Stage 3).
 * ====================================================================== */

bool XtbGpuContext::beginBasis(const XtbGpuBasisData& b)
{
    if (!ok() || b.nat <= 0 || b.nsh <= 0 || b.nao <= 0 || b.nprim_total <= 0
        || !b.z || !b.sh2at || !b.ang_sh || !b.iao_sh || !b.nao_sh
        || !b.sh_nprim || !b.sh_prim_off || !b.prim_alpha || !b.prim_coeff
        || !b.sh_zeta || !b.selfenergy || !b.kcn || !b.shpoly)
        return false;
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
        const size_t nn = static_cast<size_t>(b.nao) * static_cast<size_t>(b.nao);
        m_impl->dS.ensure(static_cast<int>(nn));
        m_impl->dH0.ensure(static_cast<int>(nn));
        m_impl->dL.ensure(static_cast<int>(nn));
        // AO→atom / AO→shell maps (used by the multipole integrals and the
        // overlap-derivative kernel; cheap, uploaded for both methods).
        if (b.ao2at) m_impl->dAo2at.upload(b.ao2at, b.nao, m_impl->stream);
        if (b.ao2sh) m_impl->dAo2sh.upload(b.ao2sh, b.nao, m_impl->stream);
        // GFN2: the resident multipole integral buffers (computed by
        // computeIntegrals; no upload). dp_int 3·nn, qp_int 6·nn col-major.
        if (b.is_gfn2) {
            m_impl->dDpInt.ensure(static_cast<int>(3 * nn));
            m_impl->dQpInt.ensure(static_cast<int>(6 * nn));
        }
        // Device Cholesky workspace (size is geometry-constant for fixed nao).
        int lwork = 0;
        if (cusolverDnDpotrf_bufferSize(m_impl->cusolver, CUBLAS_FILL_MODE_LOWER,
                                        b.nao, m_impl->dL.ptr, b.nao, &lwork)
            != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->potrf_lwork = lwork;
        m_impl->dPotrfWork.ensure(lwork > 0 ? lwork : 1);
        if (m_impl->dInfo.n < 1) m_impl->dInfo.alloc(1);
    } catch (...) {
        return false;
    }
    m_impl->basis_nat     = b.nat;
    m_impl->basis_nsh     = b.nsh;
    m_impl->basis_nao     = b.nao;
    m_impl->basis_is_gfn2 = b.is_gfn2;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
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

    // CN + self-energies (uploads xyz, runs k_cn + k_self_energy, syncs).
    if (!computeCnSelfEnergy(xyz_bohr)) return false;

    // Overlap S + bare Hamiltonian H0 (one thread per shell-pair).
    const dim3 block(16, 16);
    const dim3 grid((nsh + block.x - 1) / block.x, (nsh + block.y - 1) / block.y);
    k_overlap_h0<<<grid, block, 0, stream>>>(
        nsh, nao, m_impl->basis_is_gfn2,
        m_impl->dSh2at.ptr, m_impl->dAng.ptr, m_impl->dIaoSh.ptr, m_impl->dNaoSh.ptr,
        m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
        m_impl->dPrimCoeff.ptr, m_impl->dShZeta.ptr, m_impl->dShpoly.ptr,
        m_impl->dSE.ptr, m_impl->dZ.ptr, m_impl->dValence.ptr, m_impl->dXyz.ptr,
        m_impl->dS.ptr, m_impl->dH0.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    // Coulomb γ matrix (independent of S; one thread per shell-pair).
    if (!m_impl->dHardness.empty()) {
        const dim3 gblock(16, 16);
        const dim3 ggrid((nsh + gblock.x - 1) / gblock.x, (nsh + gblock.y - 1) / gblock.y);
        k_gamma<<<ggrid, gblock, 0, stream>>>(nsh, m_impl->basis_is_gfn2, m_impl->dSh2at.ptr,
                                              m_impl->dHardness.ptr, m_impl->dXyz.ptr,
                                              m_impl->dGamma.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
    }

    // GFN2 multipole integrals (dp_int/qp_int), one thread per AO pair. Needs the
    // resident overlap S for the origin shift, so it runs after k_overlap_h0.
    if (m_impl->basis_is_gfn2 && !m_impl->dDpInt.empty()) {
        const dim3 mblock(16, 16);
        const dim3 mgrid((nao + mblock.x - 1) / mblock.x, (nao + mblock.y - 1) / mblock.y);
        k_multipole_ints<<<mgrid, mblock, 0, stream>>>(
            nao, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dIaoSh.ptr, m_impl->dAng.ptr,
            m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
            m_impl->dPrimCoeff.ptr, m_impl->dXyz.ptr, m_impl->dS.ptr,
            m_impl->dDpInt.ptr, m_impl->dQpInt.ptr);
        if (cudaGetLastError() != cudaSuccess) return false;
    }

    // L = chol(S): copy S → L, factor the lower triangle in place (matches the
    // CPU Eigen matrixL() convention; the trsm path reads fill=LOWER only).
    const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
    if (cudaMemcpyAsync(m_impl->dL.ptr, m_impl->dS.ptr, sizeof(double) * nn,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
        return false;
    if (cusolverDnDpotrf(m_impl->cusolver, CUBLAS_FILL_MODE_LOWER, nao,
                         m_impl->dL.ptr, nao, m_impl->dPotrfWork.ptr,
                         m_impl->potrf_lwork, m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
        return false;
    int info = 1;
    if (cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
                        cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaStreamSynchronize(stream) != cudaSuccess) return false;
    return info == 0;
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
        m_impl->dC.ensure(static_cast<int>(nn));
        m_impl->dP.ensure(static_cast<int>(nn));
        m_impl->dCw.ensure(static_cast<int>(nn));
        m_impl->dEps.ensure(n);
        m_impl->dVao.ensure(n);
        m_impl->dOcc.ensure(n);
        m_impl->dPop.ensure(n);
        int lwork = 0;
        if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                        CUBLAS_FILL_MODE_LOWER, n, m_impl->dC.ptr, n,
                                        m_impl->dEps.ptr, &lwork) != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->lwork = lwork;
        m_impl->dWork.ensure(lwork > 0 ? lwork : 1);
        if (m_impl->dInfo.n < 1) m_impl->dInfo.alloc(1);
    } catch (...) {
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

bool XtbGpuContext::downloadOverlap(double* S_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !S_out) return false;
    const size_t nn = static_cast<size_t>(m_impl->basis_nao) * m_impl->basis_nao;
    m_impl->dS.download(S_out, static_cast<int>(nn), m_impl->stream);
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::downloadH0(double* H0_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !H0_out) return false;
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
    if (!ok() || m_impl->basis_nao <= 0 || m_impl->dDpInt.empty() || !dp_int3 || !qp_int6)
        return false;
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
    if (!m_impl->dMpAmatSD.empty() && e_multipole) {
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
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::residentScfStep(bool fp32, double* dq_out, double* e_band,
                                    double* e_coulomb, double* e_third, double* e_multipole)
{
    if (!ok() || m_impl->loop_nao <= 0 || !dq_out || !e_band) return false;
    cudaStream_t stream = m_impl->stream;
    const int nsh = m_impl->loop_nsh, nat = m_impl->loop_nat, nao = m_impl->loop_nao;
    const int b = 128;

    // 1. D4 reference weights from the current (previous-output) resident q_at.
    k_d4_build_refw<<<(nat + b - 1) / b, b, 0, stream>>>(
        nat, Impl::D4_MAX_REF, m_impl->dQat.ptr, m_impl->dD4Cn.ptr, m_impl->dD4Gi.ptr,
        m_impl->dD4Zeff.ptr, m_impl->dD4Nref.ptr, m_impl->dD4Refcn.ptr,
        m_impl->dD4Refcovcn.ptr, m_impl->dD4Refq.ptr, m_impl->dD4W.ptr, m_impl->dD4dWq.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    // 2. Potential build (overwrites dQat with the input-derived q_at) + Fock +
    //    eigensolve; eps stays resident (no download).
    if (!buildDevicePotentialAndSolve(nao, fp32, /*n_eig=*/0, /*eps_out=*/nullptr,
                                      /*download_eps=*/false))
        return false;

    // 3. Occupation on the device (resident eps → resident occ).
    const double kT = m_impl->loop_Tele * 3.166808e-6;
    const int use_fermi = (m_impl->loop_Tele > 0.0) ? 1 : 0;
    const int blk = 256;
    k_occupations<<<1, blk, blk * sizeof(double), stream>>>(
        m_impl->dEps.ptr, m_impl->dOcc.ptr, nao, kT, m_impl->loop_nelec,
        m_impl->loop_nocc_pairs, use_fermi, nullptr, nullptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    // 4. Density + Mulliken-AO (full columns; occ is 0 past the Fermi window).
    if (!residentDensityResident(nao, nao, e_band)) return false;

    // 5. Output charges/moments from the resident density.
    // q_sh = n0_sh − Σ_{μ∈s} pop_ao; q_at = n0_at − Σ_{μ∈A} pop_ao.
    if (cudaMemcpyAsync(m_impl->dQsh.ptr, m_impl->dN0sh.ptr, sizeof(double) * nsh,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess) return false;
    k_qsh_scatter<<<(nao + b - 1) / b, b, 0, stream>>>(nao, m_impl->dPop.ptr,
                                                       m_impl->dAo2sh.ptr, m_impl->dQsh.ptr);
    if (cudaMemcpyAsync(m_impl->dQat.ptr, m_impl->dN0at.ptr, sizeof(double) * nat,
                        cudaMemcpyDeviceToDevice, stream) != cudaSuccess) return false;
    k_qat_scatter<<<(nao + b - 1) / b, b, 0, stream>>>(nao, m_impl->dPop.ptr,
                                                       m_impl->dAo2at.ptr, m_impl->dQat.ptr);
    // Atomic multipole moments dp_at/qp_at from the resident density (resident;
    // the scatter kernel needs zeroed targets).
    if (cudaMemsetAsync(m_impl->dDpAt.ptr, 0, sizeof(double) * 3 * nat, stream) != cudaSuccess) return false;
    if (cudaMemsetAsync(m_impl->dQpAt.ptr, 0, sizeof(double) * 6 * nat, stream) != cudaSuccess) return false;
    k_multipole_moments<<<(nao + b - 1) / b, b, 0, stream>>>(
        m_impl->dDpAt.ptr, m_impl->dQpAt.ptr, m_impl->dP.ptr,
        m_impl->dDpInt.ptr, m_impl->dQpInt.ptr, m_impl->dAo2at.ptr, nao);
    if (cudaGetLastError() != cudaSuccess) return false;

    // 6. SCC energy components (resident charges/moments).
    if (!sccEnergy(nat, nsh, e_coulomb, e_third, e_multipole)) return false;

    // 7. Convergence dq = max|q_sh_out − q_sh_in| (q_sh_in = current dPotQsh).
    k_maxabsdiff<<<1, blk, blk * sizeof(double), stream>>>(m_impl->dQsh.ptr, m_impl->dPotQsh.ptr,
                                                           nsh, m_impl->dDq.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;
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
        if (!runBroydenUpdate(vin, vout, m_impl->dBroyVnext.ptr)) return false;
        cudaMemcpyAsync(m_impl->dPotQsh.ptr, m_impl->dBroyVnext.ptr, sizeof(double) * nsh, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(m_impl->dInDpAt.ptr, m_impl->dBroyVnext.ptr + off_dp, sizeof(double) * 3 * nat, cudaMemcpyDeviceToDevice, stream);
        cudaMemcpyAsync(m_impl->dInQpAt.ptr, m_impl->dBroyVnext.ptr + off_qp, sizeof(double) * 6 * nat, cudaMemcpyDeviceToDevice, stream);
    }
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::residentLoopCharges(double* q_sh, double* q_at, double* dp_at,
                                        double* qp_at, double* eps)
{
    if (!ok() || m_impl->loop_nao <= 0) return false;
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
        return false;
    }
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
    k_multipole_potential<<<(nat + b - 1) / b, b, 0, stream>>>(
        nat, m_impl->dMpAmatSD.ptr, m_impl->dMpAmatDD.ptr, m_impl->dMpAmatSQ.ptr,
        m_impl->dMpDkernel.ptr, m_impl->dMpQkernel.ptr, m_impl->dQat.ptr,
        m_impl->dInDpAt.ptr, m_impl->dInQpAt.ptr,
        m_impl->dVdp.ptr, m_impl->dVqp.ptr, m_impl->dVat.ptr);
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

    // F = H0 − ½·S·(v_ao⊕v_ao) + GFN2 multipole; then eigensolve (shared path).
    const dim3 block(16, 16);
    const dim3 grid((n + block.x - 1) / block.x, (n + block.y - 1) / block.y);
    k_build_fock_iso<<<grid, block, 0, stream>>>(m_impl->dC.ptr, m_impl->dH0.ptr,
                                                 m_impl->dS.ptr, m_impl->dVao.ptr, n);
    if (cudaGetLastError() != cudaSuccess) return false;
    k_add_fock_multipole<<<grid, block, 0, stream>>>(m_impl->dC.ptr, m_impl->dDpInt.ptr,
                                                     m_impl->dQpInt.ptr, m_impl->dVdp.ptr,
                                                     m_impl->dVqp.ptr, m_impl->dAo2at.ptr, n);
    if (cudaGetLastError() != cudaSuccess) return false;

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

    try {
        if (m_impl->dW.n < static_cast<int>(nn)) m_impl->dW.alloc(static_cast<int>(nn));
        if (m_impl->dCw.n < static_cast<int>(nn)) m_impl->dCw.alloc(static_cast<int>(nn));
        if (m_impl->dP.n  < static_cast<int>(nn)) m_impl->dP.alloc(static_cast<int>(nn));
        if (m_impl->dC.n  < static_cast<int>(nn)) m_impl->dC.alloc(static_cast<int>(nn));
        if (m_impl->dVao.n < nao) m_impl->dVao.alloc(nao);
        if (m_impl->dOcc.n < nao) m_impl->dOcc.alloc(nao);
        m_impl->dQsh.ensure(nsh);
        m_impl->dGrad.ensure(3 * nat);
        m_impl->dEdcn.ensure(nat);
    } catch (...) { return false; }

    // Upload the converged SCF state (P/C symmetric-or-column-major from host).
    // AP8: skip the nao²-sized P/C uploads when they are already resident.
    if (!pc_resident) {
        m_impl->dP.upload(P, static_cast<int>(nn), stream);
        m_impl->dC.upload(C, static_cast<int>(nn), stream);
    }
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

    // Energy-weighted density W = C_occ · diag(2·ε_occ) · C_occᵀ.
    if (nocc_orbs > 0) {
        std::vector<double> occ2(nocc_orbs);
        for (int k = 0; k < nocc_orbs; ++k) occ2[k] = 2.0 * eps[k];
        m_impl->dOcc.upload(occ2.data(), nocc_orbs, stream);
        const dim3 block(16, 16);
        const dim3 grid((nao + block.x - 1) / block.x, (nocc_orbs + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, stream>>>(m_impl->dCw.ptr, m_impl->dC.ptr,
                                                 m_impl->dOcc.ptr, nao, nocc_orbs);
        if (cudaGetLastError() != cudaSuccess) return false;
        if (cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_T, nao, nao, nocc_orbs,
                        &one, m_impl->dCw.ptr, nao, m_impl->dC.ptr, nao,
                        &zero, m_impl->dW.ptr, nao) != CUBLAS_STATUS_SUCCESS)
            return false;
    } else {
        if (cudaMemsetAsync(m_impl->dW.ptr, 0, sizeof(double) * nn, stream) != cudaSuccess)
            return false;
    }

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

    k_grad_cn_onsite<<<(nao + b1 - 1) / b1, b1, 0, stream>>>(
        nao, m_impl->dP.ptr, m_impl->dKcn.ptr, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr,
        m_impl->dEdcn.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    const dim3 block(16, 16);
    const dim3 grid((nao + block.x - 1) / block.x, (nao + block.y - 1) / block.y);
    k_grad_h0_pulay<<<grid, block, 0, stream>>>(
        nao, m_impl->basis_is_gfn2, m_impl->dAo2sh.ptr, m_impl->dAo2at.ptr, m_impl->dAng.ptr,
        m_impl->dIaoSh.ptr, m_impl->dShNprim.ptr, m_impl->dShPrimOff.ptr, m_impl->dPrimAlpha.ptr,
        m_impl->dPrimCoeff.ptr, m_impl->dShZeta.ptr, m_impl->dShpoly.ptr, m_impl->dKcn.ptr,
        m_impl->dValence.ptr, m_impl->dZ.ptr, m_impl->dSE.ptr, m_impl->dXyz.ptr, m_impl->dP.ptr,
        m_impl->dS.ptr, m_impl->dH0.ptr, m_impl->dW.ptr, m_impl->dVao.ptr,
        with_mp ? m_impl->dVdp.ptr : nullptr, with_mp ? m_impl->dVqp.ptr : nullptr,
        m_impl->dGrad.ptr, m_impl->dEdcn.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    k_grad_coulomb<<<(nsh + b1 - 1) / b1, b1, 0, stream>>>(
        nsh, m_impl->basis_is_gfn2, m_impl->dSh2at.ptr, m_impl->dHardness.ptr,
        m_impl->dQsh.ptr, m_impl->dXyz.ptr, 2.0, m_impl->dGrad.ptr);
    if (cudaGetLastError() != cudaSuccess) return false;

    m_impl->dGrad.download(grad_out, 3 * nat, stream);
    m_impl->dEdcn.download(dEdcn_out, nat, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
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
    if (!ok() || m_impl->basis_nat <= 0 || m_impl->dDpInt.empty()) return false;
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

#endif // USE_CUDA_XTB
