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

    for (int ia = 0; ia < ia_nao; ++ia) {
        const int mu = ia_start + ia;
        const int ta = d_ao_to_type(la, ia);
        if (ta < 0) continue;
        for (int jb = 0; jb < jb_nao; ++jb) {
            const int nu = jb_start + jb;
            const int tb = d_ao_to_type(lb, jb);
            if (tb < 0) continue;
            const double s_ab = (iat == jat && a == b && ta == tb)
                ? 1.0
                : d_cgto_overlap(pa_alpha, pa_coeff, npa, pb_alpha, pb_coeff, npb,
                                 xa, ya, za, xb, yb, zb, ta, tb);
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
    const int ta = d_ao_to_type(ang_sh[isha], mu - iao_sh[isha]);
    const int tb = d_ao_to_type(ang_sh[ishb], nu - iao_sh[ishb]);

    for (int k = 0; k < 3; ++k) dp_int[static_cast<size_t>(k) * nn + mn] = 0.0;
    for (int k = 0; k < 6; ++k) qp_int[static_cast<size_t>(k) * nn + mn] = 0.0;
    if (ta < 0 || tb < 0) return;

    const double* aA = prim_alpha + sh_prim_off[isha];
    const double* cA = prim_coeff + sh_prim_off[isha];
    const int npa = sh_nprim[isha];
    const double* aB = prim_alpha + sh_prim_off[ishb];
    const double* cB = prim_coeff + sh_prim_off[ishb];
    const int npb = sh_nprim[ishb];

    double Sx, D[3], Q[6];
    d_cgto_multipole(aA, cA, npa, aB, cB, npb,
                     xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                     xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], ta, tb, Sx, D, Q);

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
    const int ta = d_ao_to_type(ang_sh[isha], mu - iao_sh[isha]);
    const int tb = d_ao_to_type(ang_sh[ishb], nu - iao_sh[ishb]);

    double g[3] = {0.0, 0.0, 0.0};
    if (ta >= 0 && tb >= 0) {
        d_cgto_overlap_grad(
            prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
            prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
            xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
            xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], ta, tb, g);
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
    double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu >= nao || nu >= nao) return;
    const int iat = ao2at[mu], jat = ao2at[nu];
    if (iat >= jat) return;  // unique atom pairs, off-site only

    const int isha = ao2sh[mu], ishb = ao2sh[nu];
    const int la = ang[isha], lb = ang[ishb];
    const int ta = d_ao_to_type(la, mu - iao_sh[isha]);
    const int tb = d_ao_to_type(lb, nu - iao_sh[ishb]);
    if (ta < 0 || tb < 0) return;

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
    d_cgto_overlap_grad(
        prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
        prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
        xa, ya, za, xb, yb, zb, ta, tb, dS);

    const double sval = 2.0*Pmn*h_av - 2.0*Wmn - Pmn*(v_ao[mu] + v_ao[nu]);
    const double shp  = 2.0*Pmn*H0mn*dlog_pi_dr_r;
    const double Gx = sval*dS[0] + shp*dxij;
    const double Gy = sval*dS[1] + shp*dyij;
    const double Gz = sval*dS[2] + shp*dzij;
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

// Reduce the Fock now in dC to standard form with the cached L, solve (dsyevd),
// back-transform → generalized eigenvectors in dC, eigenvalues → eps_out. The
// device analogue of the CPU dsygst+dsyevd+dtrsm path, reusing the resident L (no
// per-iteration L upload). Shared by residentSolve and residentSolveMultipole.
bool XtbGpuContext::eigensolveResidentFock(double* eps_out)
{
    const int n = m_impl->resident_n;
    cudaStream_t stream = m_impl->stream;
    const double one = 1.0;

    // dC holds F → Ã = L⁻¹·F·L⁻ᵀ → C̃ → C.
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cusolverDnDsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                         n, m_impl->dC.ptr, n, m_impl->dEps.ptr, m_impl->dWork.ptr,
                         m_impl->lwork, m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
        return false;
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;

    int info = 1;
    if (cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
                        cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    m_impl->dEps.download(eps_out, n, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess)
        return false;
    return info == 0;
}

bool XtbGpuContext::residentSolve(const double* v_ao, int n, double* eps_out)
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

    return eigensolveResidentFock(eps_out);
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
        m_impl->dVdp.alloc(3 * nat);
        m_impl->dVqp.alloc(6 * nat);
        m_impl->dDpAt.alloc(3 * nat);
        m_impl->dQpAt.alloc(6 * nat);
    } catch (...) {
        return false;
    }
    m_impl->resident_nat = nat;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::residentSolveMultipole(const double* v_ao, const double* v_dp,
                                           const double* v_qp, int n, double* eps_out)
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

    return eigensolveResidentFock(eps_out);
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
        m_impl->dCN.alloc(b.nat);
        m_impl->dSE.alloc(b.nsh);
        m_impl->dXyz.alloc(3 * b.nat);
        m_impl->dGamma.alloc(b.nsh * b.nsh);
        const size_t nn = static_cast<size_t>(b.nao) * static_cast<size_t>(b.nao);
        m_impl->dS.alloc(static_cast<int>(nn));
        m_impl->dH0.alloc(static_cast<int>(nn));
        m_impl->dL.alloc(static_cast<int>(nn));
        // AO→atom / AO→shell maps (used by the multipole integrals and the
        // overlap-derivative kernel; cheap, uploaded for both methods).
        if (b.ao2at) m_impl->dAo2at.upload(b.ao2at, b.nao, m_impl->stream);
        if (b.ao2sh) m_impl->dAo2sh.upload(b.ao2sh, b.nao, m_impl->stream);
        // GFN2: the resident multipole integral buffers (computed by
        // computeIntegrals; no upload). dp_int 3·nn, qp_int 6·nn col-major.
        if (b.is_gfn2) {
            m_impl->dDpInt.alloc(static_cast<int>(3 * nn));
            m_impl->dQpInt.alloc(static_cast<int>(6 * nn));
        }
        // Device Cholesky workspace (size is geometry-constant for fixed nao).
        int lwork = 0;
        if (cusolverDnDpotrf_bufferSize(m_impl->cusolver, CUBLAS_FILL_MODE_LOWER,
                                        b.nao, m_impl->dL.ptr, b.nao, &lwork)
            != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->potrf_lwork = lwork;
        m_impl->dPotrfWork.alloc(lwork > 0 ? lwork : 1);
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
        // uploading any matrix.
        m_impl->dC.alloc(static_cast<int>(nn));
        m_impl->dP.alloc(static_cast<int>(nn));
        m_impl->dCw.alloc(static_cast<int>(nn));
        m_impl->dEps.alloc(n);
        m_impl->dVao.alloc(n);
        m_impl->dOcc.alloc(n);
        m_impl->dPop.alloc(n);
        int lwork = 0;
        if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                        CUBLAS_FILL_MODE_LOWER, n, m_impl->dC.ptr, n,
                                        m_impl->dEps.ptr, &lwork) != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->lwork = lwork;
        m_impl->dWork.alloc(lwork > 0 ? lwork : 1);
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
                                    double* grad_out, double* dEdcn_out)
{
    if (!ok() || m_impl->basis_nao <= 0 || !P || !C || !eps || !v_ao || !q_sh
        || !grad_out || !dEdcn_out)
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
        m_impl->dQsh.alloc(nsh);
        m_impl->dGrad.alloc(3 * nat);
        m_impl->dEdcn.alloc(nat);
    } catch (...) { return false; }

    // Upload the converged SCF state (P/C symmetric-or-column-major from host).
    m_impl->dP.upload(P, static_cast<int>(nn), stream);
    m_impl->dC.upload(C, static_cast<int>(nn), stream);
    m_impl->dVao.upload(v_ao, nao, stream);
    m_impl->dQsh.upload(q_sh, nsh, stream);

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

bool XtbGpuContext::residentBeginMultipoleComputed()
{
    if (!ok() || m_impl->basis_nat <= 0 || m_impl->dDpInt.empty()) return false;
    const int nat = m_impl->basis_nat;
    try {
        // dp_int/qp_int are already resident (computeIntegrals); dAo2at uploaded in
        // beginBasis. Allocate only the per-iteration potential / moment buffers.
        m_impl->dVdp.alloc(3 * nat);
        m_impl->dVqp.alloc(6 * nat);
        m_impl->dDpAt.alloc(3 * nat);
        m_impl->dQpAt.alloc(6 * nat);
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
