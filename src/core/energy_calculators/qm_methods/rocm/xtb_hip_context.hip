/*
 * <Native xTB ROCm/HIP Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): native GFN1/GFN2 ROCm device engine.
 *   Stage 0: device handshake (query + stream).
 *   Stage 1: solveGeneralized — the per-iteration F C = S C ε via rocSOLVER dsygvd.
 *   Stage 2: device-resident GFN1 SCF (residentBegin/Solve/Density/Finalize) — Fock
 *            build + populations via the __global__ kernels below, density via rocBLAS,
 *            eigensolve via rocSOLVER; H0/S stay resident.
 * Compiled by hipcc into a plain object (CMake add_custom_command), then linked into
 * curcuma_core with g++ — so the device kernels build, but the HIP-only --offload-arch /
 * --hip-link flags never reach the GNU link (which would drop libgomp). Links libamdhip64
 * + librocsolver + librocblas.
 */

#ifdef USE_ROCM

#include "xtb_hip_context.h"

#include <hip/hip_runtime.h>
#ifdef HAVE_ROCSOLVER
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include "xtb_hip_integrals.hiph"                 // CN/overlap/H0/gamma device helpers
#include "../parameters/xtb_params_extra.hpp"     // covalent_rad_d3_au / pauling_en / atomic_rad_au
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>

namespace curcuma {
namespace xtb {
namespace gpu {

#ifdef HAVE_ROCSOLVER
namespace {
// GFN1 isotropic Fock: F[i,j] = H0[i,j] − ½·S[i,j]·(v[i] + v[j]). Column-major n×n.
__global__ void k_fock(const double* H0, const double* S, const double* v, double* F, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;  // row
    const int j = blockIdx.y * blockDim.y + threadIdx.y;  // col
    if (i >= n || j >= n) return;
    const int idx = i + j * n;
    F[idx] = H0[idx] - 0.5 * S[idx] * (v[i] + v[j]);
}
// Column scaling: B[i,j] = A[i,j]·d[j] (forms Cw = C·diag(occ)).
__global__ void k_scale_cols(const double* A, const double* d, double* B, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n || j >= n) return;
    const int idx = i + j * n;
    B[idx] = A[idx] * d[j];
}
/**
 * @brief Copy the lower triangle of a column-major n x n matrix onto its upper triangle.
 *
 * rocblas_dsyrk writes only the requested (lower) triangle; the population/band and
 * gradient kernels read the full symmetric density, so P(i,j) = P(j,i) for i < j.
 * Claude Generated (Oct 2026).
 */
__global__ void k_mirror_lower(double* P, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;   // row
    const int j = blockIdx.y * blockDim.y + threadIdx.y;   // column
    if (i >= n || j >= n || i >= j) return;
    P[i + static_cast<size_t>(j) * n] = P[j + static_cast<size_t>(i) * n];
}
/**
 * @brief Probe kernels of the factor check: t = L^T v, y = L x, y = S v (S symmetric, lower stored).
 *
 * One thread per output row/column, O(n) work each; used by XtbHipContext::factorCheck to test
 * L L^T v = S v and L (Linv v) = v with one fixed vector v, so a corrupt Cholesky factor or
 * inverse is detected at O(n^2) cost per geometry. Claude Generated (Oct 2026).
 */
__global__ void k_chk_ltv(const double* L, const double* v, double* t, int n)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    double a = 0.0;
    for (int j = k; j < n; ++j) a += L[j + static_cast<size_t>(k) * n] * v[j];
    t[k] = a;
}
__global__ void k_chk_lv(const double* L, const double* x, double* y, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double a = 0.0;
    for (int k = 0; k <= i; ++k) a += L[i + static_cast<size_t>(k) * n] * x[k];
    y[i] = a;
}
__global__ void k_chk_sv(const double* S, const double* v, double* y, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double a = 0.0;
    for (int j = 0; j < n; ++j)
        a += (i >= j ? S[i + static_cast<size_t>(j) * n] : S[j + static_cast<size_t>(i) * n]) * v[j];
    y[i] = a;
}
// Mulliken populations + band-energy partials from the resident density P (symmetric):
//   out[i]   = Σ_j P[i,j]·S[i,j],   out[n+i] = Σ_j P[i,j]·H0[i,j].
__global__ void k_popband(const double* P, const double* S, const double* H0, double* out, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double pop = 0.0, bnd = 0.0;
    for (int j = 0; j < n; ++j) {
        const double p = P[i + j * n];
        pop += p * S[i + j * n];
        bnd += p * H0[i + j * n];
    }
    out[i] = pop; out[n + i] = bnd;
}

// ---- Stage 3 integral build kernels (ported from cuda/xtb_gpu_context.cu) ----
// Element parameter tables (covrad/pauling/arad, length 86) are passed as device
// arrays rather than __constant__ symbols.

// Coordination numbers: one thread per atom (cutoff 25 bohr).
__global__ void k_cn(int nat, const double* __restrict__ xyz, const int* __restrict__ z,
                     int is_gfn2, const double* __restrict__ covrad, double* __restrict__ cn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const double xi = xyz[3 * i + 0], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
    const double rci = covrad[z[i] - 1];
    const double cutoff2 = 25.0 * 25.0;
    double sum = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const double dx = xi - xyz[3 * j + 0], dy = yi - xyz[3 * j + 1], dz = zi - xyz[3 * j + 2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cutoff2 || r2 < 1.0e-12) continue;
        const double r  = sqrt(r2);
        const double rc = rci + covrad[z[j] - 1];
        sum += d_cn_pair(r, rc, is_gfn2 != 0);
    }
    cn[i] = sum;
}

// CN-shifted shell self-energies: se[s] = selfenergy[s] − kcn[s]·CN[sh2at[s]].
__global__ void k_self_energy(int nsh, const double* __restrict__ selfenergy,
                              const double* __restrict__ kcn, const int* __restrict__ sh2at,
                              const double* __restrict__ cn, double* __restrict__ se)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    se[s] = selfenergy[s] - kcn[s] * cn[sh2at[s]];
}

// Overlap S + bare Hamiltonian H0: one thread per shell-pair (a,b).
__global__ void k_overlap_h0(
    int nsh, int nao, int is_gfn2,
    const int* __restrict__ sh2at, const int* __restrict__ ang_sh,
    const int* __restrict__ iao_sh, const int* __restrict__ nao_sh,
    const int* __restrict__ sh_nprim, const int* __restrict__ sh_prim_off,
    const double* __restrict__ prim_alpha, const double* __restrict__ prim_coeff,
    const double* __restrict__ sh_zeta, const double* __restrict__ shpoly,
    const double* __restrict__ se, const int* __restrict__ z,
    const int* __restrict__ valence, const double* __restrict__ xyz,
    const double* __restrict__ atomic_rad, const double* __restrict__ pauling,
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
        const double rr = sqrt(sqrt(r2) / (atomic_rad[zi - 1] + atomic_rad[zj - 1]));
        const double pi_ij = (1.0 + shpoly[a] * rr) * (1.0 + shpoly[b] * rr);
        double hs;
        if (is_gfn2 == 0) {
            const bool vi = valence[a] != 0, vj = valence[b] != 0;
            if (vi && vj) {
                double den = pauling[zi - 1] - pauling[zj - 1]; den *= den;
                hs = d_kpair_gfn1(zi, zj) * d_kshell_gfn1(la, lb) * (1.0 + (-7.0e-3) * den);
            } else if (vi && !vj) {
                hs = 0.5 * (d_kshell_gfn1(la, la) + 2.85);
            } else if (!vi && vj) {
                hs = 0.5 * (d_kshell_gfn1(lb, lb) + 2.85);
            } else {
                hs = 2.85;
            }
        } else {
            double den = pauling[zi - 1] - pauling[zj - 1]; den *= den;
            const double enp = 1.0 + 2.0e-2 * den;
            const double km  = d_kshell_gfn2(la, lb) * enp;
            const double za_ = sh_zeta[a], zb_ = sh_zeta[b];
            const double zij = pow(2.0 * sqrt(za_ * zb_) / (za_ + zb_), 0.5);
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

// Coulomb γ matrix: one thread per shell-pair (gexp=2; per-shell hardness g host-built).
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
        val = g[is];
    } else if (iat == jat) {
        val = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
    } else {
        const double dx = xyz[3 * iat + 0] - xyz[3 * jat + 0];
        const double dy = xyz[3 * iat + 1] - xyz[3 * jat + 1];
        const double dz = xyz[3 * iat + 2] - xyz[3 * jat + 2];
        const double r1g = dx * dx + dy * dy + dz * dz;          // gexp=2 → R²
        const double gam = d_coulomb_average(g[is], g[js], is_gfn2 != 0);
        val = pow(r1g + pow(gam, -2.0), -0.5);
    }
    gamma[static_cast<size_t>(is) + static_cast<size_t>(js) * nsh] = val;
}

// ---- Stage 4 nuclear-gradient kernels (GFN1; grad layout [3*i+k], Eh/Bohr) ---

// WP10 (Claude Generated, Oct 2026; port of the CUDA gfix_add, xtb_gpu_context.cu, docs/GPU_REPRODUCIBILITY.md):
// with fix != 0 the gradient, dE/dCN and moment accumulators take 64-bit integers (multiples of 2^-52) added
// into the same 8 bytes, so the sum does not depend on the arrival order; k_xfix_to_double converts in place
// before the download. fix == 0 keeps the FP64 atomicAdd of the original kernels. Range +-2048 per element;
// a contribution outside it (or not finite) sets xfix_bad and the host treats the gradient as failed.
constexpr double XTB_HIP_FIX_SCALE = 4503599627370496.0;   // 2^52, as CUDA GFIX
__device__ int xfix_bad;

__device__ __forceinline__ void xfix_add(double* addr, double v)
{
    const double av = fabs(v);
    if (!(av * XTB_HIP_FIX_SCALE < 9.2e18)) {   // outside the int64 range or not finite
        atomicOr(&xfix_bad, 1);
        return;
    }
    atomicAdd(reinterpret_cast<unsigned long long*>(addr),
              static_cast<unsigned long long>(llrint(v * XTB_HIP_FIX_SCALE)));
}

/// One accumulation: fixed point (fix != 0) or the FP64 atomicAdd of the original kernels.
__device__ __forceinline__ void xacc_add(double* addr, double v, int fix)
{
    if (fix) xfix_add(addr, v);
    else     atomicAdd(addr, v);
}

/// In-place conversion of a fixed-point accumulator (see xfix_add) into doubles.
__global__ void k_xfix_to_double(int n, double* __restrict__ buf)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    buf[i] = static_cast<double>(reinterpret_cast<const long long*>(buf)[i]) * (1.0 / XTB_HIP_FIX_SCALE);
}

// Section 1: repulsion gradient. One thread per atom i, inner loop j<i.
__global__ void k_grad_repulsion(int nat, int is_gfn2, const int* __restrict__ z,
                                 const double* __restrict__ xyz, const double* __restrict__ rep_alpha,
                                 const double* __restrict__ rep_zeff, double kexp, double rexp,
                                 double kexp_light, double* __restrict__ grad, int fix)
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
        xacc_add(&grad[3*i+0], fx, fix); xacc_add(&grad[3*j+0], -fx, fix);
        xacc_add(&grad[3*i+1], fy, fix); xacc_add(&grad[3*j+1], -fy, fix);
        xacc_add(&grad[3*i+2], fz, fix); xacc_add(&grad[3*j+2], -fz, fix);
    }
}

// Section 2a: on-site CN coupling dEdcn[iat] += −kcn[ish]·P(μ,μ). One thread per AO μ.
__global__ void k_grad_cn_onsite(int nao, const double* __restrict__ P,
                                 const double* __restrict__ kcn, const int* __restrict__ ao2sh,
                                 const int* __restrict__ ao2at, double* __restrict__ dEdcn, int fix)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    const double Pmm = P[static_cast<size_t>(mu) + static_cast<size_t>(mu) * nao];
    xacc_add(&dEdcn[ao2at[mu]], (-kcn[ao2sh[mu]]) * Pmm, fix);
}

// Section 2b: H0/Pulay off-site gradient (GFN1 isotropic; no multipole). One thread
// per AO pair (μ,ν) with iat<jat.
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
    const double* __restrict__ atomic_rad, const double* __restrict__ pauling,
    double* __restrict__ grad, double* __restrict__ dEdcn, int fix)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu >= nao || nu >= nao) return;
    const int iat = ao2at[mu], jat = ao2at[nu];
    if (iat >= jat) return;

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
    const double rad_sum = atomic_rad[zi-1] + atomic_rad[zj-1];
    const double rr = sqrt(r / rad_sum);
    const double pi_a = 1.0 + shpoly[isha] * rr, pi_b = 1.0 + shpoly[ishb] * rr;

    double hs;
    if (is_gfn2 == 0) {
        const bool vi = valence[isha] != 0, vj = valence[ishb] != 0;
        if (vi && vj) {
            double den = pauling[zi-1] - pauling[zj-1]; den *= den;
            hs = d_kpair_gfn1(zi, zj) * d_kshell_gfn1(la, lb) * (1.0 + (-7.0e-3) * den);
        } else if (vi && !vj) hs = 0.5 * (d_kshell_gfn1(la, la) + 2.85);
        else if (!vi && vj)   hs = 0.5 * (d_kshell_gfn1(lb, lb) + 2.85);
        else                  hs = 2.85;
    } else {
        double den = pauling[zi-1] - pauling[zj-1]; den *= den;
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

    // GFN2 multipole-integral Pulay term (R-AP3; verbatim port of cuda k_grad_h0_pulay
    // :623-663). The transformed dp_int/qp_int derivatives are contracted with the
    // converged multipole potential v_dp/v_qp; G[l] -= Pmn·term[l]. dp_int origin is the
    // column atom jat, so dR = R_jat - R_iat.
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

    xacc_add(&grad[3*iat+0], Gx, fix); xacc_add(&grad[3*jat+0], -Gx, fix);
    xacc_add(&grad[3*iat+1], Gy, fix); xacc_add(&grad[3*jat+1], -Gy, fix);
    xacc_add(&grad[3*iat+2], Gz, fix); xacc_add(&grad[3*jat+2], -Gz, fix);

    const double cn_c = h_factor * Pmn * Smn;
    xacc_add(&dEdcn[iat], (-kcn[isha]) * cn_c, fix);
    xacc_add(&dEdcn[jat], (-kcn[ishb]) * cn_c, fix);
}

// Section 3: isotropic Coulomb gradient. One thread per shell is, inner js<is.
__global__ void k_grad_coulomb(int nsh, int is_gfn2, const int* __restrict__ sh2at,
                               const double* __restrict__ g, const double* __restrict__ q_sh,
                               const double* __restrict__ xyz, double gexp, double* __restrict__ grad, int fix)
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
        xacc_add(&grad[3*iat+0], force*dx, fix); xacc_add(&grad[3*jat+0], -force*dx, fix);
        xacc_add(&grad[3*iat+1], force*dy, fix); xacc_add(&grad[3*jat+1], -force*dy, fix);
        xacc_add(&grad[3*iat+2], force*dz, fix); xacc_add(&grad[3*jat+2], -force*dz, fix);
    }
}

// Stage 3m (R-AP1): GFN2 AO multipole integrals dp_int(3)/qp_int(6). One thread per AO
// pair (mu,nu): global-origin dipole/quadrupole (d_cgto_multipole), then the per-column
// origin shift (origin = atom of column AO nu) + traceless transform with the resident S.
// Verbatim port of CUDA k_multipole_ints. dp_int 3·nn, qp_int 6·nn, column-major (mu,nu)
// at mu+nu*nao; NOT symmetric.
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

    for (int k = 0; k < 3; ++k) dp_int[static_cast<size_t>(k) * nn + mn] = 0.0;
    for (int k = 0; k < 6; ++k) qp_int[static_cast<size_t>(k) * nn + mn] = 0.0;

    const int isha = ao2sh[mu], iat = ao2at[mu];
    const int ishb = ao2sh[nu], jat = ao2at[nu];
    const int la = ang_sh[isha], lb = ang_sh[ishb];
    const int sa = mu - iao_sh[isha], sb = nu - iao_sh[ishb];
    const bool dpair = (la >= 2 || lb >= 2);   // X-I1

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

// Stage 2b (R-AP2): GFN2 anisotropic Fock contribution, added in place to the resident
// isotropic Fock. Verbatim port of CUDA k_add_fock_multipole. One thread per AO pair.
//   F(mu,nu) -= ½·[ Σ_k dp_int[k](mu,nu)·v_dp(k,jat) + dp_int[k](nu,mu)·v_dp(k,iat)
//                 + Σ_k qp_int[k](mu,nu)·v_qp(k,jat) + qp_int[k](nu,mu)·v_qp(k,iat) ]
__global__ void k_add_fock_multipole(double* F, const double* dp_int, const double* qp_int,
                                     const double* v_dp, const double* v_qp,
                                     const int* ao2at, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu < n && nu < n) {
        const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
        const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * n;
        const size_t nm = static_cast<size_t>(nu) + static_cast<size_t>(mu) * n;
        const int iat = ao2at[mu], jat = ao2at[nu];
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

// Stage 2b: GFN2 atomic dipole/quadrupole moments from the resident density. Verbatim
// port of CUDA k_multipole_moments (HIP has atomicAdd(double), so the scatter is direct):
//   dp_at(k,iat) -= Σ_ν P(ν,μ)·dp_int[k](ν,μ) over μ∈iat. One thread per μ. Targets zeroed.
__global__ void k_multipole_moments(double* dp_at, double* qp_at,
                                    const double* P, const double* dp_int, const double* qp_int,
                                    const int* ao2at, int n, int fix)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu < n) {
        const size_t nn  = static_cast<size_t>(n) * static_cast<size_t>(n);
        const size_t col = static_cast<size_t>(mu) * n;
        const int    iat = ao2at[mu];
        for (int k = 0; k < 3; ++k) {
            const double* dk = dp_int + static_cast<size_t>(k) * nn + col;
            double acc = 0.0;
            for (int nu = 0; nu < n; ++nu) acc += P[col + nu] * dk[nu];
            xacc_add(&dp_at[k + iat * 3], -acc, fix);
        }
        for (int k = 0; k < 6; ++k) {
            const double* qk = qp_int + static_cast<size_t>(k) * nn + col;
            double acc = 0.0;
            for (int nu = 0; nu < n; ++nu) acc += P[col + nu] * qk[nu];
            xacc_add(&qp_at[k + iat * 6], -acc, fix);
        }
    }
}

// Claude Generated (Oct 2026), -xtb.rocm_device_multipole: the GFN2 multipole potential and
// energy of the host SCF loop (XTB::addMultipolePotential / XTB::energyMultipole,
// xtb_multipole.cpp) on the device, from the host-built interaction matrices.
// amat holds 15 column-major nat x nat blocks, uploaded verbatim from the host:
//   [0..2]  amat_sd[k]                    (exactly antisymmetric)
//   [3..8]  amat_dd[0][0],[0][1],[0][2],[1][1],[1][2],[2][2]   (symmetric, [a][b] == [b][a])
//   [9..14] amat_sq[k]  packed (xx, xy, yy, xz, yz, zz)       (symmetric)
// Column i is contiguous in j, so every term is read from column i as on the host:
//   amat_sd[k](i,j) = -sd_k(j,i),  amat_dd[a][b](i,j) = dd_ab(j,i),  amat_sq[k](i,j) = sq_k(j,i).
// One block per target atom i; the threads stride over j and the partial sums are combined
// by a shared-memory tree reduction with __syncthreads (no assumption on the wavefront size).
// The block size must be a power of two (kMpBlock).
//   vdp(k,i) = sum_j [ -sd_k(j,i) q_j + sum_a dd_ka(j,i) dp(a,j) ] + 2 dkernel_i dp(k,i)
//   vqp(k,i) = sum_j sq_k(j,i) q_j + 2 qkernel_i qp(k,i) mpscale_q[k]
//   v_at(i)  = sum_j [ sum_k sd_k(j,i) dp(k,j) + sum_k sq_k(j,i) qp(k,j) ]
// Moments: dp[k + 3 j], qp[k + 6 j] (column-major 3 x nat / 6 x nat as Eigen on the host).
constexpr int kMpBlock = 256;
__global__ void k_mp_potential_stored(int nat, const double* __restrict__ amat,
                                      const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                      const double* __restrict__ q, const double* __restrict__ dp,
                                      const double* __restrict__ qp,
                                      double* __restrict__ v_dp, double* __restrict__ v_qp,
                                      double* __restrict__ v_at)
{
    __shared__ double red[10][kMpBlock];
    const int i = blockIdx.x;
    const int t = threadIdx.x;
    if (i >= nat) return;   // uniform per block
    const size_t nn  = static_cast<size_t>(nat) * static_cast<size_t>(nat);
    const size_t off = static_cast<size_t>(i) * static_cast<size_t>(nat);   // column i
    const double* sd0 = amat + 0 * nn + off;
    const double* sd1 = amat + 1 * nn + off;
    const double* sd2 = amat + 2 * nn + off;
    const double* d00 = amat + 3 * nn + off;
    const double* d01 = amat + 4 * nn + off;
    const double* d02 = amat + 5 * nn + off;
    const double* d11 = amat + 6 * nn + off;
    const double* d12 = amat + 7 * nn + off;
    const double* d22 = amat + 8 * nn + off;
    const double* sq  = amat + 9 * nn + off;   // block k at sq + k*nn
    double vd0 = 0.0, vd1 = 0.0, vd2 = 0.0, acc = 0.0;
    double vq[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int j = t; j < nat; j += blockDim.x) {
        const double qj = q[j];
        const double p0 = dp[3 * j + 0], p1 = dp[3 * j + 1], p2 = dp[3 * j + 2];
        const double s0 = sd0[j], s1 = sd1[j], s2 = sd2[j];
        vd0 += -s0 * qj + d00[j] * p0 + d01[j] * p1 + d02[j] * p2;
        vd1 += -s1 * qj + d01[j] * p0 + d11[j] * p1 + d12[j] * p2;
        vd2 += -s2 * qj + d02[j] * p0 + d12[j] * p1 + d22[j] * p2;
        acc += s0 * p0 + s1 * p1 + s2 * p2;
        for (int k = 0; k < 6; ++k) {
            const double sqk = sq[static_cast<size_t>(k) * nn + j];
            vq[k] += sqk * qj;
            acc   += sqk * qp[6 * j + k];
        }
    }
    red[0][t] = vd0; red[1][t] = vd1; red[2][t] = vd2;
    for (int k = 0; k < 6; ++k) red[3 + k][t] = vq[k];
    red[9][t] = acc;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) {
        if (t < st)
            for (int c = 0; c < 10; ++c) red[c][t] += red[c][t + st];
        __syncthreads();
    }
    if (t == 0) {
        const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
        for (int k = 0; k < 3; ++k)
            v_dp[3 * i + k] = red[k][0] + 2.0 * dkernel[i] * dp[3 * i + k];
        for (int k = 0; k < 6; ++k)
            v_qp[6 * i + k] = red[3 + k][0] + 2.0 * qkernel[i] * qp[6 * i + k] * mpscale_q[k];
        v_at[i] = red[9][0];
    }
}

// Claude Generated (Oct 2026): per-atom share of the GFN2 multipole energy (XTB::energyMultipole),
// same matrix layout and reduction as k_mp_potential_stored. e_at[i] receives
//   sum_j [ -sum_k dp(k,i) sd_k(j,i) q_j + 0.5 sum_ab dp(a,i) dd_ab(j,i) dp(b,j)
//           + sum_k qp(k,i) sq_k(j,i) q_j ]
//   + dkernel_i sum_k dp(k,i)^2 + qkernel_i sum_k qp(k,i)^2 mpscale_q[k];
// the host adds the nat shares in atom order (deterministic).
__global__ void k_mp_energy_stored(int nat, const double* __restrict__ amat,
                                   const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                   const double* __restrict__ q, const double* __restrict__ dp,
                                   const double* __restrict__ qp, double* __restrict__ e_at)
{
    __shared__ double red[kMpBlock];
    const int i = blockIdx.x;
    const int t = threadIdx.x;
    if (i >= nat) return;
    const size_t nn  = static_cast<size_t>(nat) * static_cast<size_t>(nat);
    const size_t off = static_cast<size_t>(i) * static_cast<size_t>(nat);
    const double* sd0 = amat + 0 * nn + off;
    const double* sd1 = amat + 1 * nn + off;
    const double* sd2 = amat + 2 * nn + off;
    const double* d00 = amat + 3 * nn + off;
    const double* d01 = amat + 4 * nn + off;
    const double* d02 = amat + 5 * nn + off;
    const double* d11 = amat + 6 * nn + off;
    const double* d12 = amat + 7 * nn + off;
    const double* d22 = amat + 8 * nn + off;
    const double* sq  = amat + 9 * nn + off;
    const double a0 = dp[3 * i + 0], a1 = dp[3 * i + 1], a2 = dp[3 * i + 2];
    double qi[6];
    for (int k = 0; k < 6; ++k) qi[k] = qp[6 * i + k];
    double e = 0.0;
    for (int j = t; j < nat; j += blockDim.x) {
        const double qj = q[j];
        const double p0 = dp[3 * j + 0], p1 = dp[3 * j + 1], p2 = dp[3 * j + 2];
        e -= (a0 * sd0[j] + a1 * sd1[j] + a2 * sd2[j]) * qj;
        e += 0.5 * (a0 * (d00[j] * p0 + d01[j] * p1 + d02[j] * p2)
                  + a1 * (d01[j] * p0 + d11[j] * p1 + d12[j] * p2)
                  + a2 * (d02[j] * p0 + d12[j] * p1 + d22[j] * p2));
        double sqq = 0.0;
        for (int k = 0; k < 6; ++k) sqq += qi[k] * sq[static_cast<size_t>(k) * nn + j];
        e += sqq * qj;
    }
    red[t] = e;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) {
        if (t < st) red[t] += red[t + st];
        __syncthreads();
    }
    if (t == 0) {
        const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
        double es = red[0];
        for (int k = 0; k < 3; ++k) es += dkernel[i] * dp[3 * i + k] * dp[3 * i + k];
        for (int k = 0; k < 6; ++k) es += qkernel[i] * qi[k] * qi[k] * mpscale_q[k];
        e_at[i] = es;
    }
}

// ---- WP1 (Claude Generated, Oct 2026): GFN2 potential of the host SCF loop built on the device ----
// Ports of the CUDA kernels of buildDevicePotentialOnly (cuda/xtb_gpu_context.cu). All element-wise,
// one thread per shell / atom / AO / matrix element; no wavefront-size assumption.

/// Stored multipole interaction matrices from the geometry (CUDA k_mp_build_amat, same element formula as
/// the host step 5 of XTB::setupMultipole). Element (row t, col c) uses v = x_c - x_t, rr = (mrad_t+mrad_c)/2/r:
///   sd[k] = v_k g3 fdmp3,  dd[a][b] = delta_ab g3 fdmp5 - 3 v_a v_b g5 fdmp5,
///   sq    = (vx^2, 2 vx vy, vy^2, 2 vx vz, 2 vy vz, vz^2) g5 fdmp5,  fdmpN = 1/(1 + 6 rr^dmpN).
/// Written into the 15-block layout of k_mp_potential_stored (dd packed 00,01,02,11,12,22), column-major
/// (index t + c nat); the diagonal is written as zero.
__global__ void k_mp_build_amat15(int nat, const double* __restrict__ xyz, const double* __restrict__ mrad,
                                  double dmp3, double dmp5, double* __restrict__ amat)
{
    const long long tid = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    const long long nn = static_cast<long long>(nat) * nat;
    if (tid >= nn) return;
    const int t = static_cast<int>(tid % nat), c = static_cast<int>(tid / nat);
    if (t == c) { for (int b = 0; b < 15; ++b) amat[b * nn + tid] = 0.0; return; }
    const double vx = xyz[3 * c] - xyz[3 * t], vy = xyz[3 * c + 1] - xyz[3 * t + 1], vz = xyz[3 * c + 2] - xyz[3 * t + 2];
    const double r1 = sqrt(vx * vx + vy * vy + vz * vz);
    const double g1 = 1.0 / r1;
    const double g3 = g1 * g1 * g1;
    const double g5 = g3 * g1 * g1;
    const double rr = 0.5 * (mrad[t] + mrad[c]) * g1;
    const double fdmp3 = 1.0 / (1.0 + 6.0 * pow(rr, dmp3));
    const double fdmp5 = 1.0 / (1.0 + 6.0 * pow(rr, dmp5));
    amat[0 * nn + tid] = vx * g3 * fdmp3;
    amat[1 * nn + tid] = vy * g3 * fdmp3;
    amat[2 * nn + tid] = vz * g3 * fdmp3;
    const double dd_iso = g3 * fdmp5, dd_anis = 3.0 * g5 * fdmp5;
    amat[3 * nn + tid] = dd_iso - vx * vx * dd_anis;
    amat[4 * nn + tid] = -vx * vy * dd_anis;
    amat[5 * nn + tid] = -vx * vz * dd_anis;
    amat[6 * nn + tid] = dd_iso - vy * vy * dd_anis;
    amat[7 * nn + tid] = -vy * vz * dd_anis;
    amat[8 * nn + tid] = dd_iso - vz * vz * dd_anis;
    amat[9 * nn + tid]  = vx * vx * g5 * fdmp5;
    amat[10 * nn + tid] = 2.0 * vx * vy * g5 * fdmp5;
    amat[11 * nn + tid] = vy * vy * g5 * fdmp5;
    amat[12 * nn + tid] = 2.0 * vx * vz * g5 * fdmp5;
    amat[13 * nn + tid] = 2.0 * vy * vz * g5 * fdmp5;
    amat[14 * nn + tid] = vz * vz * g5 * fdmp5;
}

// ---- WP11 (Claude Generated, Oct 2026): GFN2 multipole interaction on the fly ----
// Port of the CUDA k_multipole_potential_otf / k_energy_multipole_otf (cuda/xtb_gpu_context.cu): above the 4 GB
// limit of the 15 stored nat^2 blocks (nat > ~5770), or when forced (-gpu_multipole_otf on), the kernels rebuild
// each element from the geometry instead of reading dMpAmat. d_mp_elems15 is the element formula of
// k_mp_build_amat15 (same expressions in the same order, element (row t, col c) of the 15-block layout, zero on
// the diagonal), and k_mp_potential_otf / k_mp_energy_otf are k_mp_potential_stored / k_mp_energy_stored with the
// reads replaced by it: one block per atom i, threads stride over j, the same fixed shared-memory tree. So the
// result follows the stored path to the last bit as long as the compiler evaluates the element expressions alike
// in both kernels (checked by comparison, docs/ROCM_RESIDENT_LOOP_WP.md WP11). CUDA uses a thread or a warp per
// atom with a different summation order. No atomics (WP10 fixed-point switch not involved). FP64.
__device__ __forceinline__ void d_mp_elems15(const double* __restrict__ xyz, const double* __restrict__ mrad,
                                             int t, int c, double dmp3, double dmp5, double e[15])
{
    if (t == c) { for (int b = 0; b < 15; ++b) e[b] = 0.0; return; }
    const double vx = xyz[3 * c] - xyz[3 * t], vy = xyz[3 * c + 1] - xyz[3 * t + 1], vz = xyz[3 * c + 2] - xyz[3 * t + 2];
    const double r1 = sqrt(vx * vx + vy * vy + vz * vz);
    const double g1 = 1.0 / r1;
    const double g3 = g1 * g1 * g1;
    const double g5 = g3 * g1 * g1;
    const double rr = 0.5 * (mrad[t] + mrad[c]) * g1;
    const double fdmp3 = 1.0 / (1.0 + 6.0 * pow(rr, dmp3));
    const double fdmp5 = 1.0 / (1.0 + 6.0 * pow(rr, dmp5));
    e[0] = vx * g3 * fdmp3;
    e[1] = vy * g3 * fdmp3;
    e[2] = vz * g3 * fdmp3;
    const double dd_iso = g3 * fdmp5, dd_anis = 3.0 * g5 * fdmp5;
    e[3] = dd_iso - vx * vx * dd_anis;
    e[4] = -vx * vy * dd_anis;
    e[5] = -vx * vz * dd_anis;
    e[6] = dd_iso - vy * vy * dd_anis;
    e[7] = -vy * vz * dd_anis;
    e[8] = dd_iso - vz * vz * dd_anis;
    e[9]  = vx * vx * g5 * fdmp5;
    e[10] = 2.0 * vx * vy * g5 * fdmp5;
    e[11] = vy * vy * g5 * fdmp5;
    e[12] = 2.0 * vx * vz * g5 * fdmp5;
    e[13] = 2.0 * vy * vz * g5 * fdmp5;
    e[14] = vz * vz * g5 * fdmp5;
}

/// On-the-fly twin of k_mp_potential_stored (WP11); geo = xyz (3 nat, Bohr) | mrad (nat).
__global__ void k_mp_potential_otf(int nat, const double* __restrict__ geo, double dmp3, double dmp5,
                                   const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                   const double* __restrict__ q, const double* __restrict__ dp,
                                   const double* __restrict__ qp,
                                   double* __restrict__ v_dp, double* __restrict__ v_qp,
                                   double* __restrict__ v_at)
{
    __shared__ double red[10][kMpBlock];
    const int i = blockIdx.x;
    const int t = threadIdx.x;
    if (i >= nat) return;   // uniform per block
    const double* xyz = geo;
    const double* mrad = geo + 3 * static_cast<size_t>(nat);
    double vd0 = 0.0, vd1 = 0.0, vd2 = 0.0, acc = 0.0;
    double vq[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int j = t; j < nat; j += blockDim.x) {
        double e[15];
        d_mp_elems15(xyz, mrad, j, i, dmp3, dmp5, e);   // column i, row j
        const double qj = q[j];
        const double p0 = dp[3 * j + 0], p1 = dp[3 * j + 1], p2 = dp[3 * j + 2];
        const double s0 = e[0], s1 = e[1], s2 = e[2];
        vd0 += -s0 * qj + e[3] * p0 + e[4] * p1 + e[5] * p2;
        vd1 += -s1 * qj + e[4] * p0 + e[6] * p1 + e[7] * p2;
        vd2 += -s2 * qj + e[5] * p0 + e[7] * p1 + e[8] * p2;
        acc += s0 * p0 + s1 * p1 + s2 * p2;
        for (int k = 0; k < 6; ++k) {
            const double sqk = e[9 + k];
            vq[k] += sqk * qj;
            acc   += sqk * qp[6 * j + k];
        }
    }
    red[0][t] = vd0; red[1][t] = vd1; red[2][t] = vd2;
    for (int k = 0; k < 6; ++k) red[3 + k][t] = vq[k];
    red[9][t] = acc;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) {
        if (t < st)
            for (int c = 0; c < 10; ++c) red[c][t] += red[c][t + st];
        __syncthreads();
    }
    if (t == 0) {
        const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
        for (int k = 0; k < 3; ++k)
            v_dp[3 * i + k] = red[k][0] + 2.0 * dkernel[i] * dp[3 * i + k];
        for (int k = 0; k < 6; ++k)
            v_qp[6 * i + k] = red[3 + k][0] + 2.0 * qkernel[i] * qp[6 * i + k] * mpscale_q[k];
        v_at[i] = red[9][0];
    }
}

/// On-the-fly twin of k_mp_energy_stored (WP11); per-atom shares e_at[i], summed in atom order by the caller.
__global__ void k_mp_energy_otf(int nat, const double* __restrict__ geo, double dmp3, double dmp5,
                                const double* __restrict__ dkernel, const double* __restrict__ qkernel,
                                const double* __restrict__ q, const double* __restrict__ dp,
                                const double* __restrict__ qp, double* __restrict__ e_at)
{
    __shared__ double red[kMpBlock];
    const int i = blockIdx.x;
    const int t = threadIdx.x;
    if (i >= nat) return;
    const double* xyz = geo;
    const double* mrad = geo + 3 * static_cast<size_t>(nat);
    const double a0 = dp[3 * i + 0], a1 = dp[3 * i + 1], a2 = dp[3 * i + 2];
    double qi[6];
    for (int k = 0; k < 6; ++k) qi[k] = qp[6 * i + k];
    double e = 0.0;
    for (int j = t; j < nat; j += blockDim.x) {
        double m[15];
        d_mp_elems15(xyz, mrad, j, i, dmp3, dmp5, m);
        const double qj = q[j];
        const double p0 = dp[3 * j + 0], p1 = dp[3 * j + 1], p2 = dp[3 * j + 2];
        e -= (a0 * m[0] + a1 * m[1] + a2 * m[2]) * qj;
        e += 0.5 * (a0 * (m[3] * p0 + m[4] * p1 + m[5] * p2)
                  + a1 * (m[4] * p0 + m[6] * p1 + m[7] * p2)
                  + a2 * (m[5] * p0 + m[7] * p1 + m[8] * p2));
        double sqq = 0.0;
        for (int k = 0; k < 6; ++k) sqq += qi[k] * m[9 + k];
        e += sqq * qj;
    }
    red[t] = e;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1) {
        if (t < st) red[t] += red[t + st];
        __syncthreads();
    }
    if (t == 0) {
        const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
        double es = red[0];
        for (int k = 0; k < 3; ++k) es += dkernel[i] * dp[3 * i + k] * dp[3 * i + k];
        for (int k = 0; k < 6; ++k) es += qkernel[i] * qi[k] * qi[k] * mpscale_q[k];
        e_at[i] = es;
    }
}

/// GFN2 shell third order, v_sh(s) += q_sh(s)^2 Gamma_s, on top of gamma q_sh (XTB::addThirdOrderPotential).
/// GFN1 (WP3) launches it over the atoms: v_at(A) += q_A^2 Gamma_A.
__global__ void k_vsh_third(int nsh, const double* __restrict__ q_sh, const double* __restrict__ gamma3,
                            double* __restrict__ v_sh)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    v_sh[s] += q_sh[s] * q_sh[s] * gamma3[s];
}

/// v_at(A) += dE_D4/dq(A) (XTB::addDispersionPotential), on the multipole scalar shift.
__global__ void k_vat_add_d4(int nat, const double* __restrict__ d4_dedq, double* __restrict__ v_at)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    v_at[a] += d4_dedq[a];
}

/// v_ao(mu) = v_sh(ao2sh[mu]) + v_at(ao2at[mu]) (expand_potential, xtb_scf.cpp).
__global__ void k_expand_vao(int nao, const double* __restrict__ v_sh, const double* __restrict__ v_at,
                             const int* __restrict__ ao2sh, const int* __restrict__ ao2at,
                             double* __restrict__ v_ao)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= nao) return;
    v_ao[mu] = v_sh[ao2sh[mu]] + v_at[ao2at[mu]];
}

// ---- WP2 (Claude Generated, Oct 2026): fully device-resident GFN2 SCF loop, closed shell ----
// Ports of the CUDA resident-loop kernels (cuda/xtb_gpu_context.cu, Stage 6) with three deliberate
// differences that keep the ROCm loop on the path of the host-driven loop (xtb_native.cpp):
//  - the SCC vector is packed in the host order of packSCC (q_sh | dp component-major | qp
//    component-major), the D4 reference weights are built from the mixed INPUT charges (CUDA: the
//    previous output charges), so the device Broyden sees the same sequence of vectors as BroydenMixer;
//  - charges and moments are fixed-order gathers (no atomics), in the summation order of
//    updatePopulationsFromPopAo and of the host unpackSCC;
//  - the D4 weight kernel is the current host buildAtomRefW incl. the underflow rule (all Gaussians 0:
//    the whole weight on the reference with the largest covalent CN; CUDA puts gwk = 0 there).
// Reductions use one block with __syncthreads tree steps (power-of-two block size), no wavefront-size
// assumption.

/// Unpack the mixed SCC vector X (host packing) into the per-atom inputs of the potential: q_at(A) = sum of
/// q_sh over the shells of A in ascending shell order (host unpackSCC), dp/qp atom-major (q | dp | qp in mpin).
/// moments = 0 (GFN1, WP3): X holds the shell charges alone and only q_at is written.
__global__ void k_res_unpack_atoms(int nsh, int nat, const double* __restrict__ X, const int* __restrict__ at_sh_off,
                                   double* __restrict__ mpin, int moments)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    double q = 0.0;
    for (int s = at_sh_off[a]; s < at_sh_off[a + 1]; ++s) q += X[s];
    mpin[a] = q;
    if (!moments) return;
    for (int k = 0; k < 3; ++k) mpin[nat + 3 * a + k] = X[nsh + k * nat + a];
    for (int k = 0; k < 6; ++k) mpin[4 * nat + 6 * a + k] = X[nsh + 3 * nat + k * nat + a];
}

/// Pack the output charges and atom-major moments (outmp = q_at | dp | qp) with the shell charges into the
/// host packing of packSCC.
__global__ void k_res_pack_out(int nsh, int nat, const double* __restrict__ q_sh, const double* __restrict__ outmp,
                               double* __restrict__ X)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nsh) { X[i] = q_sh[i]; return; }
    const int r = i - nsh;
    if (r < 3 * nat) { const int k = r / nat, a = r % nat; X[i] = outmp[nat + 3 * a + k]; return; }
    const int r2 = r - 3 * nat;
    if (r2 < 6 * nat) { const int k = r2 / nat, a = r2 % nat; X[i] = outmp[4 * nat + 6 * a + k]; }
}

/// D4 reference weights W = gwk(CN) zeta(q) and dW/dq per atom from the device charges q: port of
/// D4ParameterGenerator::buildAtomRefW (d4param_generator.cpp) with want_grad, incl. its underflow rule.
/// ga = 3, wf = 6, MAXCN = 19 (dftd4 defaults); the ngw bucketing uses refcn, the Gaussian refcovcn.
/// WP8a (Claude Generated, Oct 2026): dWc = d(gwk)/dCN zeta when dWc != nullptr (post-SCF D4), in the host
/// order expd/norm - gwk (dnorm/norm), non-finite -> 0; W and dWq are unchanged by it.
__global__ void k_res_d4_refw(int nat, int max_ref, const double* __restrict__ q, const double* __restrict__ cn,
                              const double* __restrict__ gi, const double* __restrict__ zeff,
                              const int* __restrict__ nref, const double* __restrict__ refcn,
                              const double* __restrict__ refcovcn, const double* __restrict__ refq,
                              double* __restrict__ W, double* __restrict__ dWq, double* __restrict__ dWc)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    constexpr double ga = 3.0, wf = 6.0;
    constexpr int MAXCN = 19;
    const size_t base = static_cast<size_t>(a) * max_ref;
    for (int ir = 0; ir < max_ref; ++ir) { W[base + ir] = 0.0; dWq[base + ir] = 0.0; if (dWc) dWc[base + ir] = 0.0; }
    int nr = nref[a];
    if (nr <= 0) return;
    if (nr > max_ref) nr = max_ref;
    if (nr > 7) nr = 7;
    const double cna = cn[a], gia = gi[a], zef = zeff[a];
    int cnc[MAXCN + 1];
    for (int k = 0; k <= MAXCN; ++k) cnc[k] = 0;
    cnc[0] = 1;
    for (int ir = 0; ir < nr; ++ir) {
        int icn = static_cast<int>(lround(refcn[base + ir]));
        if (icn < 0) icn = 0;
        if (icn > MAXCN) icn = MAXCN;
        cnc[icn] += 1;
    }
    double expw[7], expd[7];
    double norm = 0.0, dnorm = 0.0, maxcn = 0.0;
    for (int ir = 0; ir < nr; ++ir) {
        int icn = static_cast<int>(lround(refcn[base + ir]));
        if (icn < 0) icn = 0;
        if (icn > MAXCN) icn = MAXCN;
        const int k = cnc[icn];
        const int ngw = k * (k + 1) / 2;
        double ew = 0.0, ed = 0.0;
        for (int igw = 1; igw <= ngw; ++igw) {
            const double wfe = igw * wf;
            const double d = cna - refcovcn[base + ir];
            const double gw = exp(-wfe * d * d);
            ew += gw;
            if (dWc) ed += 2.0 * wfe * (refcovcn[base + ir] - cna) * gw;
        }
        expw[ir] = ew; expd[ir] = ed;
        norm += ew; dnorm += ed;
        maxcn = fmax(maxcn, refcovcn[base + ir]);
    }
    const double ninv = 1.0 / norm;
    const double dnorm_n = dnorm * ninv;
    const double qmod = q[a] + zef;
    for (int ir = 0; ir < nr; ++ir) {
        double gwk = expw[ir] * ninv;
        double dgwk = expd[ir] * ninv - gwk * dnorm_n;
        if (!isfinite(gwk)) gwk = (refcovcn[base + ir] == maxcn) ? 1.0 : 0.0;
        if (!isfinite(dgwk)) dgwk = 0.0;
        const double qref = refq[base + ir] + zef;
        double z, dz;
        if (qmod < 0.0) { z = exp(ga); dz = 0.0; }
        else {
            const double g = exp(gia * (1.0 - qref / qmod));
            z = exp(ga * (1.0 - g));
            dz = -ga * gia * g * z * qref / (qmod * qmod);
        }
        W[base + ir] = gwk * z;
        dWq[base + ir] = gwk * dz;
        if (dWc) dWc[base + ir] = dgwk * z;
    }
}

/// Closed-shell occupations of the resident eigenvalues (XTB::occupationsFromEps): Fermi smearing by
/// bisection on mu (100 steps, stop below 1e-14) or the integer filling of nocc_pairs orbitals at T = 0.
/// Port of CUDA k_occupations; one block, sdata = blockDim doubles. ncol = last orbital with occ > 1e-12.
__global__ void k_res_occupations(const double* __restrict__ eps, double* __restrict__ occ, int n, double kT,
                                  double n_elec, int nocc_pairs, int use_fermi, int* ncol_out)
{
    extern __shared__ double sdata[];
    __shared__ double mu_lo_s, mu_hi_s;
    const int tid = threadIdx.x;
    const int nt = blockDim.x;
    if (!use_fermi) {
        const int m = max(0, min(nocc_pairs, n));
        for (int i = tid; i < n; i += nt) occ[i] = (i < m) ? 2.0 : 0.0;
        if (tid == 0) *ncol_out = m;
        return;
    }
    double vmin = 1e300, vmax = -1e300;
    for (int i = tid; i < n; i += nt) { vmin = fmin(vmin, eps[i]); vmax = fmax(vmax, eps[i]); }
    sdata[tid] = vmin; __syncthreads();
    for (int s = nt >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmin(sdata[tid], sdata[tid + s]); __syncthreads(); }
    const double eps_min = sdata[0]; __syncthreads();
    sdata[tid] = vmax; __syncthreads();
    for (int s = nt >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    const double eps_max = sdata[0]; __syncthreads();
    if (tid == 0) { mu_lo_s = eps_min - 1.0; mu_hi_s = eps_max + 1.0; }
    __syncthreads();
    for (int bisect = 0; bisect < 100; ++bisect) {
        const double mu = 0.5 * (mu_lo_s + mu_hi_s);
        double psum = 0.0;
        for (int i = tid; i < n; i += nt) psum += 2.0 / (1.0 + exp(fmin((eps[i] - mu) / kT, 500.0)));
        sdata[tid] = psum; __syncthreads();
        for (int s = nt >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] += sdata[tid + s]; __syncthreads(); }
        const double n_sum = sdata[0];
        __syncthreads();
        if (tid == 0) { if (n_sum > n_elec) mu_hi_s = mu; else mu_lo_s = mu; }
        __syncthreads();
        if (mu_hi_s - mu_lo_s < 1e-14) break;
    }
    const double mu_f = 0.5 * (mu_lo_s + mu_hi_s);
    int local_ncol = 0;
    for (int i = tid; i < n; i += nt) {
        const double o = 2.0 / (1.0 + exp(fmin((eps[i] - mu_f) / kT, 500.0)));
        occ[i] = o;
        if (o > 1.0e-12) local_ncol = i + 1;
    }
    __syncthreads();
    sdata[tid] = static_cast<double>(local_ncol); __syncthreads();
    for (int s = nt >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    if (tid == 0) *ncol_out = static_cast<int>(sdata[0]);
}

/**
 * @brief One spin channel of the open-shell filling (WP4): aufbau filling with nel electrons, then (kT > 0)
 *        Newton iterations on the Fermi level towards the integer target homo. Port of XTB::fermiFillingChannel
 *        (xtb_scf.cpp; tblite get_aufbau_filling + get_fermi_filling_) and of the CUDA occupation_channel: same
 *        start value, x < 50 cut, 200 cycles, stop on total_dfermi == 0 or |occt - total| <= sqrt(eps).
 *        Unlike CUDA the number and derivative sums run in ascending orbital order on thread 0, the order of the
 *        host loop. dfs = scratch (n). Claude Generated (Oct 2026).
 */
__device__ void res_fermi_channel(const double* __restrict__ eps, double* __restrict__ out, double* __restrict__ dfs,
                                  int n, double nel, double kT)
{
    __shared__ double s_ef, s_tot, s_dtot;
    __shared__ int s_stop;
    const int tid = threadIdx.x, nt = blockDim.x;
    for (int i = tid; i < n; i += nt) out[i] = 0.0;
    __syncthreads();
    if (n == 0 || nel <= 0.0) return;
    int homo = static_cast<int>(floor(nel));
    for (int i = tid; i < min(homo, n); i += nt) out[i] = 1.0;
    const double frac = nel - floor(nel);
    if (tid == 0 && homo < n) out[homo] = frac;
    if (frac > 0.5) ++homo;
    __syncthreads();
    if (homo <= 0 || kT <= 0.0) return;   // T = 0: the aufbau filling is the answer
    if (tid == 0) s_ef = 0.5 * (eps[max(homo, 1) - 1] + eps[min(homo + 1, n) - 1]);
    __syncthreads();
    const double occt = static_cast<double>(homo);
    const double thr = 1.4901161193847656e-08;   // sqrt(epsilon(double))
    for (int cycle = 0; cycle < 200; ++cycle) {
        const double ef = s_ef;
        for (int i = tid; i < n; i += nt) {
            double f = 0.0, df = 0.0;
            const double x = (eps[i] - ef) / kT;
            if (x < 50.0) {
                const double ex = exp(x);
                f = 1.0 / (ex + 1.0);
                df = ex / (kT * (ex + 1.0) * (ex + 1.0));
            }
            out[i] = f;
            dfs[i] = df;
        }
        __syncthreads();
        if (tid == 0) {
            double tn = 0.0, td = 0.0;
            for (int i = 0; i < n; ++i) { tn += out[i]; td += dfs[i]; }
            s_tot = tn; s_dtot = td;
            s_stop = (td == 0.0) ? 1 : 0;
            if (!s_stop) {
                s_ef = ef + (occt - tn) / td;
                if (fabs(occt - tn) <= thr) s_stop = 1;
            }
        }
        __syncthreads();
        if (s_stop) break;
    }
    __syncthreads();
}

/// Open-shell occupations of the resident eigenvalues (WP4): two fillings (n_alpha, n_beta) over the same
/// orbitals, occ = occ_a + occ_b (XTB::occupationsFromEps, openShell branch; CUDA k_occupations_uhf).
/// One block; ncol = last orbital with occ > 1e-12. Claude Generated (Oct 2026).
__global__ void k_res_occupations_uhf(const double* __restrict__ eps, double* __restrict__ occ,
                                      double* __restrict__ occ_a, double* __restrict__ occ_b,
                                      double* __restrict__ dfs, int n, double kT, double n_alpha, double n_beta,
                                      int* ncol_out)
{
    extern __shared__ double sdata[];
    const int tid = threadIdx.x, nt = blockDim.x;
    res_fermi_channel(eps, occ_a, dfs, n, n_alpha, kT);
    res_fermi_channel(eps, occ_b, dfs, n, n_beta, kT);
    int local_ncol = 0;
    for (int i = tid; i < n; i += nt) {
        const double o = occ_a[i] + occ_b[i];
        occ[i] = o;
        if (o > 1.0e-12) local_ncol = i + 1;
    }
    __syncthreads();
    sdata[tid] = static_cast<double>(local_ncol); __syncthreads();
    for (int s = nt >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    if (tid == 0) *ncol_out = static_cast<int>(sdata[0]);
}

/// Scale vector of the density product: occ (GEMM) or sqrt(occ) (SYRK).
__global__ void k_res_occ_scale(int n, const double* __restrict__ occ, int use_sqrt, double* __restrict__ out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = use_sqrt ? sqrt(fmax(occ[i], 0.0)) : occ[i];
}

/// Shell populations n_sh(s) = sum of pop_ao over the AOs of s (ascending), q_sh = n0_sh - n_sh
/// (updatePopulationsFromPopAo, xtb_scf.cpp).
__global__ void k_res_qsh(int nsh, const double* __restrict__ pop, const int* __restrict__ iao_sh,
                          const int* __restrict__ nao_sh, const double* __restrict__ n0_sh,
                          double* __restrict__ n_sh, double* __restrict__ q_sh)
{
    const int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nsh) return;
    double acc = 0.0;
    for (int mu = iao_sh[s]; mu < iao_sh[s] + nao_sh[s]; ++mu) acc += pop[mu];
    n_sh[s] = acc;
    q_sh[s] = n0_sh[s] - acc;
}

/// Atomic charges q_at(A) = n0_at(A) - sum of n_sh over the shells of A (ascending), as on the host.
__global__ void k_res_qat(int nat, const int* __restrict__ at_sh_off, const double* __restrict__ n_sh,
                          const double* __restrict__ n0_at, double* __restrict__ q_at)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    double acc = 0.0;
    for (int s = at_sh_off[a]; s < at_sh_off[a + 1]; ++s) acc += n_sh[s];
    q_at[a] = n0_at[a] - acc;
}

/// Per-AO partials of the atomic moments (k_multipole_moments without the scatter): part[9 mu + k] =
/// sum_nu P(nu, mu) int_k(nu, mu), k = 0..2 dipole, 3..8 quadrupole.
__global__ void k_res_mom_part(int n, const double* __restrict__ P, const double* __restrict__ dp_int,
                               const double* __restrict__ qp_int, double* __restrict__ part)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu >= n) return;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    const size_t col = static_cast<size_t>(mu) * n;
    for (int k = 0; k < 3; ++k) {
        const double* dk = dp_int + static_cast<size_t>(k) * nn + col;
        double acc = 0.0;
        for (int nu = 0; nu < n; ++nu) acc += P[col + nu] * dk[nu];
        part[9 * static_cast<size_t>(mu) + k] = acc;
    }
    for (int k = 0; k < 6; ++k) {
        const double* qk = qp_int + static_cast<size_t>(k) * nn + col;
        double acc = 0.0;
        for (int nu = 0; nu < n; ++nu) acc += P[col + nu] * qk[nu];
        part[9 * static_cast<size_t>(mu) + 3 + k] = acc;
    }
}

/// Atomic moments from the per-AO partials in ascending AO order: dp(k, A) = -sum part, atom-major.
__global__ void k_res_mom_gather(int nat, const int* __restrict__ at_ao_off, const double* __restrict__ part,
                                 double* __restrict__ dp, double* __restrict__ qp)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    double d[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int mu = at_ao_off[a]; mu < at_ao_off[a + 1]; ++mu)
        for (int k = 0; k < 9; ++k) d[k] -= part[9 * static_cast<size_t>(mu) + k];
    for (int k = 0; k < 3; ++k) dp[3 * a + k] = d[k];
    for (int k = 0; k < 6; ++k) qp[6 * a + k] = d[3 + k];
}

/// Sequential sums in index order (one thread; the vectors are length nao, nsh or nat):
///   mode 0: out = sum x;  mode 1: out = sum x*y;  mode 2: out = sum x^3 y / 3 (third order: GFN2 over the
///   shells, GFN1 over the atoms, coulomb::energy_third_order).
__global__ void k_res_seqsum(int n, const double* __restrict__ x, const double* __restrict__ y, int mode,
                             double* __restrict__ out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    double s = 0.0;
    if (mode == 0)      for (int i = 0; i < n; ++i) s += x[i];
    else if (mode == 1) for (int i = 0; i < n; ++i) s += x[i] * y[i];
    else                for (int i = 0; i < n; ++i) s += x[i] * x[i] * x[i] * y[i] / 3.0;
    *out = s;
}

/// max_i |a_i - b_i| (the SCF convergence max|dq| over the whole SCC vector); one block, tree max.
__global__ void k_res_maxabsdiff(const double* __restrict__ a, const double* __restrict__ b, int n,
                                 double* __restrict__ out)
{
    extern __shared__ double sdata[];
    const int tid = threadIdx.x;
    double m = 0.0;
    for (int i = tid; i < n; i += blockDim.x) m = fmax(m, fabs(a[i] - b[i]));
    sdata[tid] = m; __syncthreads();
    for (int s = blockDim.x >> 1; s > 0; s >>= 1) { if (tid < s) sdata[tid] = fmax(sdata[tid], sdata[tid + s]); __syncthreads(); }
    if (tid == 0) *out = sdata[0];
}

/// Broyden mixer kernels (BroydenMixer::update, broyden_mixer.h): out = a - b; the normalised residual change
/// dF = dFraw / |dFraw| and u = alpha dF + (vin - vin_last) / |dFraw| into the history columns.
__global__ void k_res_vec_sub(double* __restrict__ out, const double* __restrict__ a, const double* __restrict__ b, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = a[i] - b[i];
}
__global__ void k_res_broyden_dfu(double* __restrict__ dfcol, double* __restrict__ ucol,
                                  const double* __restrict__ dfraw, const double* __restrict__ vin,
                                  const double* __restrict__ vinlast, double alpha, double inv_norm, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double dfn = dfraw[i] * inv_norm;
    dfcol[i] = dfn;
    ucol[i] = alpha * dfn + (vin[i] - vinlast[i]) * inv_norm;
}
/// (w0^2 I + a) gamma = c by Gaussian elimination with partial pivoting, one thread (M <= 20; the host
/// forms the inverse, so the two agree to rounding).
__global__ void k_res_broyden_solve(int M, const double* __restrict__ a, const double* __restrict__ c,
                                    double reg, double* __restrict__ gamma)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    double A[20][20];
    double b[20];
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

// ---- WP5: Broyden seed preconditioner M (Claude Generated, Oct 2026; port of the CUDA k_pc_diag, k_pc_blocks,
// k_broyden_u_pc, k_broyden_u_unprecond; docs/SCF_PRECOND.md) ----
/// out = C v, C = fallback (alpha on indices outside every Schwarz block, 0 inside).
__global__ void k_res_pc_diag(double* __restrict__ out, const double* __restrict__ fb, const double* __restrict__ v,
                              int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = fb[i] * v[i];
}
/// out += R_b^T B_b R_b v, one thread block per Schwarz block b: v_b = R_b v gathered into shared memory, one
/// thread per row of B_b (row-major), the row sum in column order (no warp-size assumption). The own rows of the
/// blocks are disjoint and the rows of non-own atoms of a segment block are stored as exact zeros, so every out[i]
/// gets at most one non-zero add and the atomic add does not make the result order dependent.
__global__ void k_res_pc_blocks(double* __restrict__ out, const double* __restrict__ v, const double* __restrict__ B,
                                const int* __restrict__ boff, const int* __restrict__ ioff,
                                const int* __restrict__ idx)
{
    extern __shared__ double vl[];
    const int b = blockIdx.x;
    const int i0 = ioff[b];
    const int n = ioff[b + 1] - i0;
    const double* Bb = B + boff[b];
    for (int k = threadIdx.x; k < n; k += blockDim.x) vl[k] = v[idx[i0 + k]];
    __syncthreads();
    for (int r = threadIdx.x; r < n; r += blockDim.x) {
        const double* row = Bb + static_cast<size_t>(r) * n;
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += row[j] * vl[j];
        atomicAdd(&out[idx[i0 + r]], s);
    }
}
/// u = M dF + (vin - vin_last) / |dF| (mdf = M dF already formed): BroydenMixer::update with a step hook.
__global__ void k_res_broyden_u_pc(double* __restrict__ ucol, const double* __restrict__ mdf,
                                   const double* __restrict__ vin, const double* __restrict__ vinlast,
                                   double inv_norm, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) ucol[i] = mdf[i] + (vin[i] - vinlast[i]) * inv_norm;
}
/// Guard hit with "keep": u_k += alpha dF_k - M dF_k (BroydenMixer::dropStepHook(true)).
__global__ void k_res_broyden_u_unprecond(double* __restrict__ ucol, const double* __restrict__ dfcol,
                                          const double* __restrict__ mdf, double alpha, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) ucol[i] += alpha * dfcol[i] - mdf[i];
}

// FP32 mixed-precision (X-AP3): element-wise FP64<->FP32 conversion for the optional
// rocsolver_ssygvd eigensolve path (far-from-convergence iterations of -scf_mixed_precision).
__global__ void k_d2f(const double* __restrict__ src, float* __restrict__ dst, size_t n)
{ size_t i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) dst[i] = static_cast<float>(src[i]); }
__global__ void k_f2d(const float* __restrict__ src, double* __restrict__ dst, size_t n)
{ size_t i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) dst[i] = static_cast<double>(src[i]); }

// ---- WP7: pseudo-diagonalisation on the resident loop (Claude Generated, Oct 2026) ----------------
// Port of the CUDA kernels k_pseudo_coldot / k_pseudo_gap / k_pseudo_theta / k_pseudo_gram_diag
// (cuda/xtb_gpu_context.cu, stage 2 of docs/GFN2_GPU_COST_PLAN.md). Between two SCF steps the previous
// eigenvectors C (C^T S C = I, occupied block first) are close to the new ones; one first-order Jacobi
// rotation between the occupied (o) and virtual (v) blocks,
//     G = C^T F C,   theta_ai = G_ai / (G_ii - G_aa),   C_o += C_v theta,   C_v -= C_o theta^T
// (Stewart, Csaszar, Pulay, J. Comput. Chem. 3 (1982) 227), and a Cholesky re-orthonormalisation in the
// S metric replace the eigensolve. All block reductions are trees over blockDim (a power of two) with
// __syncthreads, no warp-size assumption (Known Issue F-17).

/// g[p] = sum_mu C(mu,p) T(mu,p) in FP64; one block per column, sred = blockDim doubles.
template <typename R>
__global__ void k_pseudo_coldot(const R* __restrict__ C, const R* __restrict__ T, int n, double* __restrict__ g)
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

/// stats[0] = max_i g_i (occupied), stats[1] = min_a g_a (virtual); one block, sred = 2 blockDim doubles.
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

/// theta_ai = G_ai / (g_i - g_a) in place (G_ov nv x no, column-major) and the largest |theta| as the
/// bit pattern of a non-negative float (monotone as unsigned int, so atomicMax orders it).
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

/// Diagonal of the S-metric Gram matrix M_pp = sum_mu B(mu,p)^2 in FP64, written over the syrk diagonal
/// (the FP32 syrk diagonal is about 1e-6 too small per column; CUDA measurement, see the CUDA kernel).
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

/// Precision overloads of the rocBLAS/rocSOLVER calls of the pseudo step (CUDA pseudo_blas).
namespace pseudo_blas {
inline rocblas_status gemm(rocblas_handle h, rocblas_operation a, rocblas_operation b, int m, int n, int k,
                           const float* al, const float* A, int lda, const float* B, int ldb,
                           const float* be, float* C, int ldc)
{ return rocblas_sgemm(h, a, b, m, n, k, al, A, lda, B, ldb, be, C, ldc); }
inline rocblas_status gemm(rocblas_handle h, rocblas_operation a, rocblas_operation b, int m, int n, int k,
                           const double* al, const double* A, int lda, const double* B, int ldb,
                           const double* be, double* C, int ldc)
{ return rocblas_dgemm(h, a, b, m, n, k, al, A, lda, B, ldb, be, C, ldc); }
// B = L^T C (out of place, L lower)
inline rocblas_status trmm(rocblas_handle h, int n, const float* al, const float* L, const float* C, float* B)
{ return rocblas_strmm(h, rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose, rocblas_diagonal_non_unit,
                       n, n, al, L, n, C, n, B, n); }
inline rocblas_status trmm(rocblas_handle h, int n, const double* al, const double* L, const double* C, double* B)
{ return rocblas_dtrmm(h, rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose, rocblas_diagonal_non_unit,
                       n, n, al, L, n, C, n, B, n); }
// M = B^T B (upper)
inline rocblas_status syrk(rocblas_handle h, int n, const float* al, const float* B, const float* be, float* M)
{ return rocblas_ssyrk(h, rocblas_fill_upper, rocblas_operation_transpose, n, n, al, B, n, be, M, n); }
inline rocblas_status syrk(rocblas_handle h, int n, const double* al, const double* B, const double* be, double* M)
{ return rocblas_dsyrk(h, rocblas_fill_upper, rocblas_operation_transpose, n, n, al, B, n, be, M, n); }
// C <- C R^-1 (R upper)
inline rocblas_status trsm(rocblas_handle h, int n, const float* al, const float* R, float* C)
{ return rocblas_strsm(h, rocblas_side_right, rocblas_fill_upper, rocblas_operation_none, rocblas_diagonal_non_unit,
                       n, n, al, R, n, C, n); }
inline rocblas_status trsm(rocblas_handle h, int n, const double* al, const double* R, double* C)
{ return rocblas_dtrsm(h, rocblas_side_right, rocblas_fill_upper, rocblas_operation_none, rocblas_diagonal_non_unit,
                       n, n, al, R, n, C, n); }
// M = R^T R (upper)
inline rocblas_status potrf(rocblas_handle h, int n, float* M, rocblas_int* info)
{ return rocsolver_spotrf(h, rocblas_fill_upper, n, M, n, info); }
inline rocblas_status potrf(rocblas_handle h, int n, double* M, rocblas_int* info)
{ return rocsolver_dpotrf(h, rocblas_fill_upper, n, M, n, info); }
} // namespace pseudo_blas

// Stage 5 (Part B2): in-SCF GFN2 D4 atom-potential dE_D4/dq, one thread per atom
// i (no atomicAdd — each atom sums its own half of the symmetric pair contributions).
// Verbatim port of the CUDA k_d4_dedq (cuda/xtb_gpu_context.cu). For every j≠i within
// the cutoff:
//   r4r2ij = 3·sqrtZr4r2_i·sqrtZr4r2_j ; R0 = a1·√r4r2ij + a2
//   disp_sum = s6/(r⁶+R0⁶) + s8·r4r2ij/(r⁸+R0⁸)              (BJ, geometry-fixed)
//   dc6/dq_i = ΣΣ dWq_i[a]·W_j[b]·c6ref(ei,ej,a,b)            (7×7 contraction)
//   dEdq(i) += −dc6/dq_i · disp_sum
// Mirrors D4Evaluator::computeEnergyAndGradient's per-reference path + contractC6Gfn2.
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

// Post-SCF: the whole 2-body D4 (energy + nuclear gradient + dE/dCN + dE/dq) in one per-atom
// GATHER — a superset of k_d4_dedq. No atomics: the pair force is antisymmetric (atom j's own
// thread accumulates the mirror term) and the energy is per-atom (host halves Σ e_atom). Port of
// D4Evaluator::computeEnergyAndGradient's per-reference 2-body path: the direct radial gradient is
// dE_dr = -C6·ddisp_dr2·(R_i-R_j) with ddisp_dr2 = s6·(-6r⁴t6²) + s8·r4r2·(-8r⁶t8²) (the C6·ζc6
// rescale cancels in the per-reference path), dEdcn_i = -dc6dcn_i·disp_sum, dEdq_i = -dc6dq_i·disp_sum.
// Claude Generated.
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

// D4 ATM 3-body (energy + nuclear gradient + dE/dCN) — per-atom GATHER over unordered pairs
// {x,y} (the triple {a,x,y}). The triple-edge force is antisymmetric (dG_ba = -dG_ab; the dftd4
// dang is symmetric in the edge endpoints), so atom `a` accumulates only ITS member's force from
// each triple → no atomics, each triple's energy counted once per member (host: E = ⅓·Σ e_atom).
// Verbatim port of D4Evaluator::computeATM (dftd4 atm.f90). q=0 reference C6/dC6dCN supplied as
// nat² matrices (c6 symmetric; dc6dcn[a·nat+b] = ∂C6(a,b)/∂CN_a). Claude Generated.
__global__ void k_d4_atm(int nat, const double* __restrict__ xyz, const double* __restrict__ r4r2,
                         const double* __restrict__ c6, const double* __restrict__ dc6dcn,
                         double s9, double a1, double a2, double alp, double cut2,
                         double* __restrict__ e_atom, double* __restrict__ grad,
                         double* __restrict__ dEdcn)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    const double eps = 2.220446049250313e-16;   // std::numeric_limits<double>::epsilon()
    const double xa = xyz[3*a+0], ya = xyz[3*a+1], za = xyz[3*a+2];
    const double r4r2a = r4r2[a];
    double eacc = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, cnacc = 0.0;
    for (int x = 0; x < nat; ++x) {
        if (x == a) continue;
        const double vaxx = xyz[3*x+0]-xa, vaxy = xyz[3*x+1]-ya, vaxz = xyz[3*x+2]-za;  // v_ax = r_x-r_a
        const double r2ax = vaxx*vaxx + vaxy*vaxy + vaxz*vaxz;
        if (r2ax > cut2 || r2ax < eps) continue;
        const double c6ax = c6[(size_t)a*nat + x];
        const double r0ax = a1*sqrt(3.0*r4r2a*r4r2[x]) + a2;
        for (int y = 0; y < x; ++y) {
            if (y == a) continue;
            const double vayx = xyz[3*y+0]-xa, vayy = xyz[3*y+1]-ya, vayz = xyz[3*y+2]-za;  // v_ay = r_y-r_a
            const double r2ay = vayx*vayx + vayy*vayy + vayz*vayz;
            if (r2ay > cut2 || r2ay < eps) continue;
            const double dxyx = xyz[3*y+0]-xyz[3*x+0], dxyy = xyz[3*y+1]-xyz[3*x+1], dxyz_ = xyz[3*y+2]-xyz[3*x+2];
            const double r2xy = dxyx*dxyx + dxyy*dxyy + dxyz_*dxyz_;
            if (r2xy > cut2 || r2xy < eps) continue;

            const double c6ay = c6[(size_t)a*nat + y];
            const double c6xy = c6[(size_t)x*nat + y];
            const double c9 = -s9 * sqrt(fabs(c6ax*c6ay*c6xy));
            const double r0ay = a1*sqrt(3.0*r4r2a*r4r2[y]) + a2;
            const double r0xy = a1*sqrt(3.0*r4r2[x]*r4r2[y]) + a2;
            const double r0 = r0ax*r0ay*r0xy;
            const double r2 = r2ax*r2ay*r2xy;
            const double r1 = sqrt(r2);
            const double r3 = r2*r1;
            const double r5 = r3*r2;
            const double pw = pow(r0/r1, alp/3.0);
            const double fdmp = 1.0/(1.0 + 6.0*pw);
            // ang symmetric; (r2ij,r2jk,r2ik) = (r2ax, r2xy, r2ay)  [i=a,j=x,k=y]
            const double ang = 0.375*(r2ax + r2xy - r2ay)*(r2ax - r2xy + r2ay)*(-r2ax + r2xy + r2ay)/r5 + 1.0/r3;
            const double dE = ang*fdmp*c9;       // triple = 1 (distinct atoms)
            eacc += -dE;

            const double dfdmp = -2.0*alp*pw*fdmp*fdmp;
            // edge a-x (third y): dang(r2ax; {r2xy,r2ay})
            const double dang_ax = -0.375*(r2ax*r2ax*r2ax + r2ax*r2ax*(r2xy+r2ay)
                          + r2ax*(3.0*r2xy*r2xy + 2.0*r2xy*r2ay + 3.0*r2ay*r2ay)
                          - 5.0*(r2xy-r2ay)*(r2xy-r2ay)*(r2xy+r2ay))/r5;
            const double gcax = c9*(-dang_ax*fdmp + ang*dfdmp)/r2ax;   // dG_ax = gcax·v_ax
            gx += -gcax*vaxx; gy += -gcax*vaxy; gz += -gcax*vaxz;        // a's force = -dG_ax
            // edge a-y (third x): dang(r2ay; {r2xy,r2ax})
            const double dang_ay = -0.375*(r2ay*r2ay*r2ay + r2ay*r2ay*(r2xy+r2ax)
                          + r2ay*(3.0*r2xy*r2xy + 2.0*r2xy*r2ax + 3.0*r2ax*r2ax)
                          - 5.0*(r2xy-r2ax)*(r2xy-r2ax)*(r2xy+r2ax))/r5;
            const double gcay = c9*(-dang_ay*fdmp + ang*dfdmp)/r2ay;
            gx += -gcay*vayx; gy += -gcay*vayy; gz += -gcay*vayz;
            // CN chain for a: -dE·0.5·(dc6dcn(a,x)/c6ax + dc6dcn(a,y)/c6ay)
            cnacc += -dE*0.5*(dc6dcn[(size_t)a*nat+x]/c6ax + dc6dcn[(size_t)a*nat+y]/c6ay);
        }
    }
    e_atom[a] = eacc;
    grad[3*a+0] = gx; grad[3*a+1] = gy; grad[3*a+2] = gz;
    dEdcn[a] = cnacc;
}

/* ====================================================================== *
 *  WP8a (Claude Generated, Oct 2026): post-SCF dispersion and CN on the device.
 *  Ports of the CUDA kernels of XtbGpuContext::d4CovalentCnDevice / d4PostScfDevice /
 *  d3DispersionDevice (cuda/xtb_gpu_context.cu). Same formulas as the host reference:
 *  curcuma::dispersion::computeD4CovalentCN / addD4CovalentCNGradient (d4_ncoord.cpp),
 *  D4ParameterGenerator::buildAtmC6Flat (d4param_generator.cpp), D4Evaluator::computeATM,
 *  D3ParameterGenerator::prepareForEnergyGradient / getEnergyAndGradient / computeDC6DCN and
 *  CNCalculator::calculateD3CN / addD3CNGradient. Every atom gathers over its partners (no
 *  atomics, fixed order), FP64 throughout. Thread per atom: no warp-size assumption.
 * ====================================================================== */

/// EN-weighted covalent CN of dftd4 (computeD4CovalentCN): CN_i = sum_j k4 exp(-(|EN_i-EN_j|+k5)^2/k6)
/// * 0.5 (1 + erf(-kn (r/(R_i+R_j) - 1))), r <= cutoff (Bohr). rad/en indexed by Z (0..118).
__global__ void k_d4_cn(int nat, const int* __restrict__ Z, const double* __restrict__ xyz,
                        const double* __restrict__ rad, const double* __restrict__ en, double kn, double k4,
                        double k5, double k6, double cutoff, double* __restrict__ cn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int zi = Z[i];
    double c = 0.0;
    if (zi >= 1 && zi <= 118) {
        for (int j = 0; j < nat; ++j) {
            if (j == i) continue;
            const int zj = Z[j];
            if (zj < 1 || zj > 118) continue;
            const double dx = xyz[3*i] - xyz[3*j], dy = xyz[3*i+1] - xyz[3*j+1], dz = xyz[3*i+2] - xyz[3*j+2];
            const double r = sqrt(dx*dx + dy*dy + dz*dz);
            if (r > cutoff) continue;
            const double rr = r / (rad[zi] + rad[zj]);
            const double t = fabs(en[zi] - en[zj]) + k5;
            const double den = k4 * exp(-(t * t) / k6);
            c += den * (0.5 * (1.0 + erf(-kn * (rr - 1.0))));
        }
    }
    cn[i] = c;
}

/// grad_i += sum_j (dEdcn_i + dEdcn_j) dcount_ij/dr (R_i - R_j)/r with the dftd4 covalent CN
/// (addD4CovalentCNGradient), gathered per atom.
__global__ void k_d4_cn_fold(int nat, const int* __restrict__ Z, const double* __restrict__ xyz,
                             const double* __restrict__ dEdcn, const double* __restrict__ rad,
                             const double* __restrict__ en, double kn, double k4, double k5, double k6,
                             double cutoff, double* __restrict__ grad)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int zi = Z[i];
    if (zi < 1 || zi > 118) return;
    const double hlfosqrtpi = 1.0 / 1.7724538509055159;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const int zj = Z[j];
        if (zj < 1 || zj > 118) continue;
        const double rx = xyz[3*i] - xyz[3*j], ry = xyz[3*i+1] - xyz[3*j+1], rz = xyz[3*i+2] - xyz[3*j+2];
        const double r = sqrt(rx*rx + ry*ry + rz*rz);
        if (r > cutoff || r < 1e-12) continue;
        const double rcovij = rad[zi] + rad[zj];
        const double rr = r / rcovij;
        const double t = fabs(en[zi] - en[zj]) + k5;
        const double den = k4 * exp(-(t * t) / k6);
        const double u = kn * (rr - 1.0);
        const double dcountf = den * (-kn * hlfosqrtpi * exp(-(u * u))) / rcovij;
        const double factor = (dEdcn[i] + dEdcn[j]) * dcountf / r;
        gx += factor * rx; gy += factor * ry; gz += factor * rz;
    }
    grad[3*i] += gx; grad[3*i+1] += gy; grad[3*i+2] += gz;
}

/// q = 0 reference C6 and dC6/dCN matrices of the ATM term (buildAtmC6Flat): c6[a nat + b] = c6[b nat + a],
/// dc6[a nat + b] = dC6(a,b)/dCN_a, dc6[b nat + a] = dC6(a,b)/dCN_b. One thread per (a, b <= a).
__global__ void k_d4_atm_c6mat(int nat, const int* __restrict__ Z, const int* __restrict__ nref,
                               const double* __restrict__ W0, const double* __restrict__ dWc0,
                               const double* __restrict__ c6flat,
                               double* __restrict__ c6m, double* __restrict__ dc6m)
{
    const long long t = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (t >= static_cast<long long>(nat) * nat) return;
    const int a = static_cast<int>(t / nat), b = static_cast<int>(t % nat);
    if (b > a) return;
    constexpr int MR = 7, ME = 118;
    const int ea = Z[a] - 1, eb = Z[b] - 1;
    double c6 = 0.0, di = 0.0, dj = 0.0;
    if (ea >= 0 && ea < ME && eb >= 0 && eb < ME) {
        const int nra = min(nref[a], MR), nrb = min(nref[b], MR);
        const size_t base = static_cast<size_t>(ea) * ME * MR * MR + static_cast<size_t>(eb) * MR * MR;
        const double* Wa = W0 + static_cast<size_t>(a) * MR;
        const double* Wb = W0 + static_cast<size_t>(b) * MR;
        const double* dWa = dWc0 + static_cast<size_t>(a) * MR;
        const double* dWb = dWc0 + static_cast<size_t>(b) * MR;
        for (int ra = 0; ra < nra; ++ra)
            for (int rb = 0; rb < nrb; ++rb) {
                const double c6ref = c6flat[base + static_cast<size_t>(ra) * MR + rb];
                c6 += Wa[ra] * Wb[rb] * c6ref;
                di += dWa[ra] * Wb[rb] * c6ref;
                dj += Wa[ra] * dWb[rb] * c6ref;
            }
    }
    const size_t nat_s = static_cast<size_t>(nat);
    c6m[a * nat_s + b] = c6;  c6m[b * nat_s + a] = c6;
    dc6m[a * nat_s + b] = di; dc6m[b * nat_s + a] = dj;
}

/// ATM neighbour list: atoms x != a with |r_x - r_a|^2 <= cut2, ascending (count, then fill at the offsets).
__global__ void k_d4_nb_count(int nat, const double* __restrict__ xyz, double cut2, int* __restrict__ cnt)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    int c = 0;
    for (int x = 0; x < nat; ++x) {
        if (x == a) continue;
        const double dx = xyz[3*x] - xyz[3*a], dy = xyz[3*x+1] - xyz[3*a+1], dz = xyz[3*x+2] - xyz[3*a+2];
        if (dx*dx + dy*dy + dz*dz <= cut2) ++c;
    }
    cnt[a] = c;
}
__global__ void k_d4_nb_fill(int nat, const double* __restrict__ xyz, double cut2,
                             const int* __restrict__ off, int* __restrict__ nb)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    int p = off[a];
    for (int x = 0; x < nat; ++x) {
        if (x == a) continue;
        const double dx = xyz[3*x] - xyz[3*a], dy = xyz[3*x+1] - xyz[3*a+1], dz = xyz[3*x+2] - xyz[3*a+2];
        if (dx*dx + dy*dy + dz*dz <= cut2) nb[p++] = x;
    }
}

/// ATM three-body over the neighbour list: k_d4_atm with the pair loops restricted to nb[nbptr[a]..nbptr[a+1])
/// (ascending, the same pairs in the same order, so the sums equal those of k_d4_atm bit for bit; pow as the
/// host, not the cbrt shortcut of the CUDA k_d4_atm_nl).
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
            const double pw = pow(r0/r1, alp/3.0);
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

/// Block-per-atom twin of k_d4_atm_nl (port of the CUDA warp-per-atom k_d4_atm_nl_w without a warp-size
/// assumption): the threads of one block stride over the neighbour index ix of atom a = blockIdx.x, each keeps
/// the inner y < x loop; the five partial sums are reduced by a fixed shared-memory tree (deterministic, but a
/// different summation order than the host loop, so the result agrees with it up to rounding).
/// blockDim.x must be a power of two <= 256; shared memory 5 * blockDim.x doubles.
__global__ void k_d4_atm_nl_b(int nat, const double* __restrict__ xyz, const double* __restrict__ r4r2,
                              const double* __restrict__ c6, const double* __restrict__ dc6dcn,
                              const int* __restrict__ nbptr, const int* __restrict__ nb,
                              double s9, double a1, double a2, double alp, double cut2,
                              double* __restrict__ e_atom, double* __restrict__ grad,
                              double* __restrict__ dEdcn)
{
    extern __shared__ double sred[];
    const int a = blockIdx.x;
    const int tid = threadIdx.x, nt = blockDim.x;
    if (a >= nat) return;   // uniform per block
    const double eps = 2.220446049250313e-16;
    const double xa = xyz[3*a+0], ya = xyz[3*a+1], za = xyz[3*a+2];
    const double r4r2a = r4r2[a];
    const size_t nat_s = static_cast<size_t>(nat);
    double eacc = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, cnacc = 0.0;
    for (int ix = nbptr[a] + tid; ix < nbptr[a + 1]; ix += nt) {
        const int x = nb[ix];
        const double vaxx = xyz[3*x+0]-xa, vaxy = xyz[3*x+1]-ya, vaxz = xyz[3*x+2]-za;
        const double r2ax = vaxx*vaxx + vaxy*vaxy + vaxz*vaxz;
        if (r2ax > cut2 || r2ax < eps) continue;
        const double c6ax = c6[a*nat_s + x];
        const double r0ax = a1*sqrt(3.0*r4r2a*r4r2[x]) + a2;
        for (int iy = nbptr[a]; iy < ix; ++iy) {
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
            const double pw = pow(r0/r1, alp/3.0);
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
    sred[tid] = eacc; sred[nt + tid] = gx; sred[2*nt + tid] = gy; sred[3*nt + tid] = gz; sred[4*nt + tid] = cnacc;
    __syncthreads();
    for (int s = nt / 2; s > 0; s >>= 1) {
        if (tid < s)
            for (int k = 0; k < 5; ++k) sred[k*nt + tid] += sred[k*nt + tid + s];
        __syncthreads();
    }
    if (tid != 0) return;
    e_atom[a] = sred[0];
    grad[3*a+0] = sred[nt]; grad[3*a+1] = sred[2*nt]; grad[3*a+2] = sred[3*nt];
    dEdcn[a] = sred[4*nt];
}

/// out[i] = a[i] + b[i]
__global__ void k_d4_add2(int n, const double* __restrict__ a, const double* __restrict__ b, double* __restrict__ out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = a[i] + b[i];
}

// ---- GFN1 D3(BJ): CN, Gaussian weights, pair energy / gradient / dE/dCN, CN chain rule (geometry in Angstrom) ----
constexpr int D3_MR = 7, D3_ME = 103;   // D3ParameterGenerator::MAX_REF / MAX_ELEM

/// D3 CN (CNCalculator::calculateD3CN): CN_i = sum_j 1/(1 + exp(-k1 (k2 (R_i+R_j)/r - 1))), r <= cut_ang.
__global__ void k_d3_cn(int nat, const int* __restrict__ Z, const double* __restrict__ xyz,
                        const double* __restrict__ rcov, int nrcov, double k1, double k2, double cut_ang,
                        double* __restrict__ cn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int ei = Z[i] - 1;
    double c = 0.0;
    if (ei >= 0 && ei < nrcov) {
        const double rci = rcov[ei];
        for (int j = 0; j < nat; ++j) {
            if (j == i) continue;
            const int ej = Z[j] - 1;
            if (ej < 0 || ej >= nrcov) continue;
            const double dx = xyz[3*i] - xyz[3*j], dy = xyz[3*i+1] - xyz[3*j+1], dz = xyz[3*i+2] - xyz[3*j+2];
            const double r = sqrt(dx*dx + dy*dy + dz*dz);
            if (cut_ang > 0.0 && r > cut_ang) continue;
            const double arg = -k1 * (k2 * (rci + rcov[ej]) / r - 1.0);
            c += 1.0 / (1.0 + exp(arg));
        }
    }
    cn[i] = c;
}

/// D3 Gaussian weights and dW/dCN (precomputeGaussianWeights / computeDC6DCN, wf = 4): normalised by the sum,
/// non-finite weight -> 1 on the reference with the largest CN (s-dftd3 is_exceptional), non-finite derivative -> 0.
__global__ void k_d3_weights(int nat, const int* __restrict__ Z, const double* __restrict__ cn,
                             const int* __restrict__ nref_tab, const double* __restrict__ refcn_tab,
                             double* __restrict__ gw, double* __restrict__ dgw)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    for (int r = 0; r < D3_MR; ++r) { gw[i * D3_MR + r] = 0.0; dgw[i * D3_MR + r] = 0.0; }
    const int e = Z[i] - 1;
    if (e < 0 || e >= D3_ME) return;
    const int nref = min(nref_tab[e], D3_MR);
    if (nref <= 0) return;
    const double wf = 4.0;
    double raw[D3_MR], rawd[D3_MR];
    double sum = 0.0, dsum = 0.0, max_cn_ref = refcn_tab[e * D3_MR];
    for (int r = 0; r < nref; ++r) {
        const double cref = refcn_tab[e * D3_MR + r];
        const double diff = cn[i] - cref;
        raw[r] = exp(-wf * diff * diff);
        rawd[r] = -2.0 * wf * diff * raw[r];
        sum += raw[r]; dsum += rawd[r];
        if (r > 0) max_cn_ref = fmax(max_cn_ref, cref);
    }
    const double inv = 1.0 / sum;
    for (int r = 0; r < nref; ++r) {
        double g = raw[r] * inv;
        if (!isfinite(g)) g = (refcn_tab[e * D3_MR + r] == max_cn_ref) ? 1.0 : 0.0;
        gw[i * D3_MR + r] = g;
        double d = (rawd[r] * sum - raw[r] * dsum) / (sum * sum);
        if (!isfinite(d)) d = 0.0;
        dgw[i * D3_MR + r] = d;
    }
}

/// D3(BJ) pair terms gathered over all partners j of atom i (getEnergyAndGradient): e_atom (pair energy, each
/// pair seen from both atoms, so E = 1/2 sum), direct gradient (Eh/Bohr from Angstrom distances), dE/dCN_i.
__global__ void k_d3_pairs(int nat, int nt, const int* __restrict__ Z, const int* __restrict__ atype,
                           const int* __restrict__ nref_tab, const double* __restrict__ xyz,
                           const double* __restrict__ gw, const double* __restrict__ dgw,
                           const double* __restrict__ c6tab, const double* __restrict__ r6tab,
                           double s6, double s8, double a1, double a2, double bohr_per_ang, int need_grad,
                           double* __restrict__ e_atom, double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    e_atom[i] = 0.0; grad[3*i] = grad[3*i+1] = grad[3*i+2] = 0.0; dEdcn[i] = 0.0;
    const int ti = atype[i];
    if (ti < 0) return;
    const int ni = min(nref_tab[Z[i] - 1], D3_MR);
    if (ni <= 0) return;
    double eacc = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, dcnacc = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const int tj = atype[j];
        if (tj < 0) continue;
        const int nj = min(nref_tab[Z[j] - 1], D3_MR);
        if (nj <= 0) continue;
        const double* T = c6tab + (static_cast<size_t>(ti) * nt + tj) * (D3_MR * D3_MR);
        double c6 = 0.0, dc6i = 0.0;
        for (int ra = 0; ra < ni; ++ra) {
            const double wi = gw[i * D3_MR + ra], dwi = dgw[i * D3_MR + ra];
            for (int rb = 0; rb < nj; ++rb) {
                const double c6ref = T[ra * D3_MR + rb];
                c6 += wi * gw[j * D3_MR + rb] * c6ref;
                if (need_grad) dc6i += dwi * gw[j * D3_MR + rb] * c6ref;
            }
        }
        if (c6 <= 0.0) continue;
        const double r6_ratio = r6tab[static_cast<size_t>(ti) * nt + tj];
        const double c8 = c6 * r6_ratio;
        const double dx = xyz[3*i] - xyz[3*j], dy = xyz[3*i+1] - xyz[3*j+1], dz = xyz[3*i+2] - xyz[3*j+2];
        const double r_ang = sqrt(dx*dx + dy*dy + dz*dz);
        if (r_ang < 1e-6) continue;
        const double r = r_ang * bohr_per_ang;
        const double r2 = r * r;
        const double r0 = a1 * sqrt(r6_ratio) + a2;
        const double r6_0 = pow(r0, 6.0), r8_0 = pow(r0, 8.0);
        const double t6 = 1.0 / (pow(r, 6.0) + r6_0);
        const double t8 = 1.0 / (pow(r, 8.0) + r8_0);
        eacc += (-s6 * c6 * t6) + (-s8 * c8 * t8);
        if (!need_grad) continue;
        const double d6 = -6.0 * r2 * r2 * t6 * t6;
        const double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
        const double gdisp = s6 * d6 + s8 * r6_ratio * d8;
        const double fac = -c6 * gdisp;
        gx += fac * (dx * bohr_per_ang); gy += fac * (dy * bohr_per_ang); gz += fac * (dz * bohr_per_ang);
        const double edisp = s6 * t6 + s8 * r6_ratio * t8;
        dcnacc -= dc6i * edisp;
    }
    e_atom[i] = eacc;
    grad[3*i] = gx; grad[3*i+1] = gy; grad[3*i+2] = gz;
    dEdcn[i] = dcnacc;
}

/// D3 CN chain rule (addD3CNGradient): grad_i += sum_j (dEdcn_i + dEdcn_j) dcount/dr (R_i - R_j)/r * unit.
__global__ void k_d3_cn_fold(int nat, const int* __restrict__ Z, const double* __restrict__ xyz,
                             const double* __restrict__ dEdcn, const double* __restrict__ rcov, int nrcov,
                             double k1, double k2, double unit, double cut_ang, double* __restrict__ grad)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int ei = Z[i] - 1;
    if (ei < 0 || ei >= nrcov) return;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int j = 0; j < nat; ++j) {
        if (j == i) continue;
        const int ej = Z[j] - 1;
        if (ej < 0 || ej >= nrcov) continue;
        const double rx = xyz[3*i] - xyz[3*j], ry = xyz[3*i+1] - xyz[3*j+1], rz = xyz[3*i+2] - xyz[3*j+2];
        const double r = sqrt(rx*rx + ry*ry + rz*rz);
        if (r < 1e-10) continue;
        if (cut_ang > 0.0 && r > cut_ang) continue;
        const double r0 = rcov[ei] + rcov[ej];
        const double count = 1.0 / (1.0 + exp(-k1 * (k2 * r0 / r - 1.0)));
        const double dcount_dr = -count * (1.0 - count) * k1 * k2 * r0 / (r * r);
        const double factor = (dEdcn[i] + dEdcn[j]) * dcount_dr / r;
        gx += factor * rx * unit; gy += factor * ry * unit; gz += factor * rz * unit;
    }
    grad[3*i] += gx; grad[3*i+1] += gy; grad[3*i+2] += gz;
}

// ===== Single-shot D4 EEQ charge model (Stage 5, Part A) — verbatim port of the CUDA
// k_d4eeq_* kernels (cuda/xtb_gpu_context.cu). The host D4ChargeModel does one (N+1) augmented
// LU solve for the charges + an adjoint solve for ∂q/∂x; here the matrix build / RHS / CN / the
// adjoint response are HIP kernels and the LU is rocSOLVER (dgetrf/dgetrs). Claude Generated.
constexpr double D4EEQ_TSQRT2PI        = 0.797884560802866;   // sqrt(2/π)
constexpr double D4EEQ_TWO_OVER_SQRTPI = 1.1283791670955126;  // 2/sqrt(π)
constexpr double D4EEQ_KN     = -7.5;       // GFN-FF erf-CN steepness
constexpr double D4EEQ_CNMAX  = 4.4;        // CN log-compression cap
constexpr double D4EEQ_CN_EPS = 1.0e-10;    // guard for 1/sqrt(CN)

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

// Augmented EEQ matrix M = [[A,1],[1ᵀ,0]], column-major (N+1)×(N+1).
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

// Per-atom b-term response weight u_i = z_q(i)·κ_i/(2√CN_i)·g_i.
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

// Charge-response gradient, one thread per atom a (no atomicAdd): the pair coeff is symmetric
// and û flips sign, so the full b≠a loop per atom reproduces both halves. grad [3a+k], Eh/Bohr.
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
        const double gammij = 1.0 / sqrt(alp_a + alpha_sq[b]);
        const double gr = gammij * r;
        const double Aprime = gammij * D4EEQ_TWO_OVER_SQRTPI * exp(-gr * gr) * inv_r
                              - erf(gr) * inv_r * inv_r;
        const double c_ab = zq[a] * q[b] + zq[b] * q[a];
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

// ---------------------------------------------------------------------------------------------
// WP8b (Claude Generated, Oct 2026): repulsion energy and the O(nat^2) tail of the GPU gradient on the device.
// Port of the CUDA k_rep_energy, d_mp_grad_pair, k_grad_tail_mp, k_grad_tail_dcn, k_grad_tail_cn
// (cuda/xtb_gpu_context.cu), same formulas as XTB::calcRepulsionEnergy (xtb_h0.cpp) and XTB::calculateGradientGpu
// sections 5 and 4 (xtb_gradient.cpp). Differences from the host: the repulsion is summed per atom i over j < i and
// the per-atom sums are added on the host in atom order (host: one running sum); the gradient kernels are per-atom
// gathers that evaluate each pair from both sides (host: each pair once, added to both atoms, threaded stripes).
// One thread per atom, no warp-size assumption.
// ---------------------------------------------------------------------------------------------

/// e_atom[i] = sum_{j<i} Z_i Z_j / R^rexp exp(-sqrt(a_i a_j) R^kexp), kexp_light for H/He pairs in GFN2.
__global__ void k_rep_energy(int nat, int is_gfn2, const int* __restrict__ z, const double* __restrict__ xyz,
                             const double* __restrict__ alpha, const double* __restrict__ zeff,
                             double kexp, double rexp, double kexp_light, double* __restrict__ e_atom)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const int zi = z[i];
    double e = 0.0;
    for (int j = 0; j < i; ++j) {
        const double dx = xyz[3*i] - xyz[3*j], dy = xyz[3*i+1] - xyz[3*j+1], dz = xyz[3*i+2] - xyz[3*j+2];
        const double r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1.0e-12) continue;
        const double alpha_pair = sqrt(alpha[i] * alpha[j]);
        double kexp_pair = kexp;
        if (is_gfn2 && zi <= 2 && z[j] <= 2) kexp_pair = kexp_light;
        e += (zeff[i] * zeff[j] / pow(r, rexp)) * exp(-alpha_pair * pow(r, kexp_pair));
    }
    e_atom[i] = e;
}

/// One pair (iat, jat) of the GFN2 direct multipole gradient as the host loop of section 5 evaluates it:
/// g (added to iat, subtracted from jat) and fddr = dE/dr_mp of the pair (added to both). dp [k + 3 atom],
/// qp [k + 6 atom].
__device__ __forceinline__ void d_mp_grad_pair(int iat, int jat, const double* __restrict__ xyz,
                                               const double* __restrict__ mrad, double dmp3, double dmp5,
                                               const double* __restrict__ q, const double* __restrict__ dp,
                                               const double* __restrict__ qp, double& gx, double& gy,
                                               double& gz, double& fddr, bool& valid)
{
    const double vx = xyz[3*jat+0] - xyz[3*iat+0];
    const double vy = xyz[3*jat+1] - xyz[3*iat+1];
    const double vz = xyz[3*jat+2] - xyz[3*iat+2];
    const double r2 = vx*vx + vy*vy + vz*vz;
    valid = !(r2 < 1.0e-12);
    gx = gy = gz = fddr = 0.0;
    if (!valid) return;
    const double r1 = sqrt(r2);
    const double g1 = 1.0 / r1;
    const double g3 = g1*g1*g1, g5 = g3*g1*g1, g7 = g5*g1*g1;
    const double rr_dmp = 0.5 * (mrad[iat] + mrad[jat]);
    const double fdmp3 = 1.0 / (1.0 + 6.0 * pow(rr_dmp * g1, dmp3));
    const double fdmp5 = 1.0 / (1.0 + 6.0 * pow(rr_dmp * g1, dmp5));
    const double ddmp3 = g5 * fdmp3 * (-3.0 + dmp3 * (1.0 - fdmp3));
    const double ddmp5 = fdmp5 * (-5.0 - dmp5 * (fdmp5 - 1.0));
    const double* dpi = dp + 3 * iat; const double* dpj = dp + 3 * jat;
    const double* qpi = qp + 6 * iat; const double* qpj = qp + 6 * jat;
    const double qi = q[iat], qj = q[jat];

    const double dpiqj = (vx*dpi[0]+vy*dpi[1]+vz*dpi[2]) * qj;
    const double qidpj = (vx*dpj[0]+vy*dpj[1]+vz*dpj[2]) * qi;
    const double diff_sd = dpiqj - qidpj;
    gx = -ddmp3*vx*diff_sd + fdmp3*g3*(qi*dpj[0] - qj*dpi[0]);
    gy = -ddmp3*vy*diff_sd + fdmp3*g3*(qi*dpj[1] - qj*dpi[1]);
    gz = -ddmp3*vz*diff_sd + fdmp3*g3*(qi*dpj[2] - qj*dpi[2]);
    const double fddr_sd = 3.0*diff_sd*dmp3*fdmp3*g3*(fdmp3/rr_dmp)*pow(rr_dmp*g1, dmp3);

    const double dpidpj = dpi[0]*dpj[0] + dpi[1]*dpj[1] + dpi[2]*dpj[2];
    const double dpiv = vx*dpi[0] + vy*dpi[1] + vz*dpi[2];
    const double dpjv = vx*dpj[0] + vy*dpj[1] + vz*dpj[2];
    const double edd = dpidpj*r2 - 3.0*dpjv*dpiv;
    gx += -2.0*fdmp5*g5*dpidpj*vx + 3.0*fdmp5*g5*(dpiv*dpj[0]+dpjv*dpi[0]) - edd*ddmp5*g7*vx;
    gy += -2.0*fdmp5*g5*dpidpj*vy + 3.0*fdmp5*g5*(dpiv*dpj[1]+dpjv*dpi[1]) - edd*ddmp5*g7*vy;
    gz += -2.0*fdmp5*g5*dpidpj*vz + 3.0*fdmp5*g5*(dpiv*dpj[2]+dpjv*dpi[2]) - edd*ddmp5*g7*vz;
    const double fddr_dd = 3.0*edd*dmp5*fdmp5*g5*(fdmp5/rr_dmp)*pow(rr_dmp*g1, dmp5);

    const double eq =
          (qj*qpi[0]+qpj[0]*qi)*vx*vx
        + 2.0*(qj*qpi[1]+qpj[1]*qi)*vx*vy
        + (qj*qpi[2]+qpj[2]*qi)*vy*vy
        + 2.0*(qj*qpi[3]+qpj[3]*qi)*vx*vz
        + 2.0*(qj*qpi[4]+qpj[4]*qi)*vy*vz
        + (qj*qpi[5]+qpj[5]*qi)*vz*vz;
    const double tivx = qpi[0]*vx+qpi[1]*vy+qpi[3]*vz;
    const double tivy = qpi[1]*vx+qpi[2]*vy+qpi[4]*vz;
    const double tivz = qpi[3]*vx+qpi[4]*vy+qpi[5]*vz;
    const double tjvx = qpj[0]*vx+qpj[1]*vy+qpj[3]*vz;
    const double tjvy = qpj[1]*vx+qpj[2]*vy+qpj[4]*vz;
    const double tjvz = qpj[3]*vx+qpj[4]*vy+qpj[5]*vz;
    gx += -eq*ddmp5*g7*vx - 2.0*fdmp5*g5*(qi*tjvx+qj*tivx);
    gy += -eq*ddmp5*g7*vy - 2.0*fdmp5*g5*(qi*tjvy+qj*tivy);
    gz += -eq*ddmp5*g7*vz - 2.0*fdmp5*g5*(qi*tjvz+qj*tivz);
    const double fddr_sq = eq*3.0*dmp5*fdmp5*g5*(fdmp5/rr_dmp)*pow(rr_dmp*g1, dmp5);
    fddr = fddr_sd + fddr_dd + fddr_sq;
}

/// Section 5 gather: grad[3a+k] = multipole gradient on atom a, dEdr[a] = dE/dr_mp(a).
__global__ void k_grad_tail_mp(int nat, const double* __restrict__ xyz, const double* __restrict__ mrad,
                               double dmp3, double dmp5, const double* __restrict__ q,
                               const double* __restrict__ dp, const double* __restrict__ qp,
                               double* __restrict__ grad, double* __restrict__ dEdr)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    double gx = 0.0, gy = 0.0, gz = 0.0, fd = 0.0;
    for (int b = 0; b < nat; ++b) {
        if (b == a) continue;
        double px, py, pz, pf; bool ok;
        if (b < a) {                    // host pair (iat = a, jat = b): +g on a
            d_mp_grad_pair(a, b, xyz, mrad, dmp3, dmp5, q, dp, qp, px, py, pz, pf, ok);
            if (!ok) continue;
            gx += px; gy += py; gz += pz;
        } else {                        // host pair (iat = b, jat = a): -g on a
            d_mp_grad_pair(b, a, xyz, mrad, dmp3, dmp5, q, dp, qp, px, py, pz, pf, ok);
            if (!ok) continue;
            gx -= px; gy -= py; gz -= pz;
        }
        fd += pf;
    }
    grad[3*a] = gx; grad[3*a+1] = gy; grad[3*a+2] = gz;
    dEdr[a] = fd;
}

/// dE/dCN += dE/dr_mp * dmrad/dCN (host: t2 = (rmax - rad) / (1 + t1), dmrad/dCN = -t2 kexp t1 / (1 + t1)).
__global__ void k_grad_tail_dcn(int nat, const double* __restrict__ dEdr, const double* __restrict__ cn,
                                const double* __restrict__ rad, const double* __restrict__ vcn,
                                double shift, double kexp, double rmax, double* __restrict__ dEdcn)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nat) return;
    const double arg = cn[i] - vcn[i] - shift;
    const double t1 = exp(-kexp * arg);
    const double t2 = (rmax - rad[i]) / (1.0 + t1);
    dEdcn[i] += dEdr[i] * (-t2 * kexp * t1 / (1.0 + t1));
}

/// Section 4 gather: CN chain rule grad[3a+k] += sum_b (dEdcn_a + dEdcn_b) dcount/dr (r_a - r_b)_k / r.
__global__ void k_grad_tail_cn(int nat, int is_gfn1, const double* __restrict__ xyz,
                               const double* __restrict__ rcov, const double* __restrict__ dEdcn,
                               double* __restrict__ grad)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nat) return;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int b = 0; b < nat; ++b) {
        if (b == a) continue;
        const double dx = xyz[3*a] - xyz[3*b], dy = xyz[3*a+1] - xyz[3*b+1], dz = xyz[3*a+2] - xyz[3*b+2];
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1.0e-12) continue;
        const double r = sqrt(r2);
        const double rc = rcov[a] + rcov[b];
        double dcndr;
        if (is_gfn1) {
            const double k = 16.0;
            const double c = 1.0 / (1.0 + exp(-k * (rc / r - 1.0)));
            dcndr = -c * (1.0 - c) * k * rc / r2;
        } else {
            const double ka = 10.0, kb = 20.0, rs = 2.0;
            const double c1 = 1.0 / (1.0 + exp(-ka * (rc / r - 1.0)));
            const double c2 = 1.0 / (1.0 + exp(-kb * ((rc + rs) / r - 1.0)));
            const double dc1dr = -c1 * (1.0 - c1) * ka * rc / r2;
            const double dc2dr = -c2 * (1.0 - c2) * kb * (rc + rs) / r2;
            dcndr = dc1dr * c2 + c1 * dc2dr;
        }
        const double factor = (dEdcn[a] + dEdcn[b]) * dcndr / r;
        gx += factor * dx; gy += factor * dy; gz += factor * dz;
    }
    grad[3*a] += gx; grad[3*a+1] += gy; grad[3*a+2] += gz;
}
} // namespace
#endif // HAVE_ROCSOLVER

struct XtbHipContext::Impl {
    bool        ok = false;
    int         device = -1;
    std::string name;
    hipStream_t stream = nullptr;

#ifdef HAVE_ROCSOLVER
    rocblas_handle handle = nullptr;   // rocSOLVER reuses the rocBLAS handle
    bool    handle_ok = false;
    int     solver_n  = 0;             // cached device-buffer size
    double* dF   = nullptr;            // F → eigenvectors (n×n)
    double* dS   = nullptr;            // S (overwritten by its Cholesky factor)
    double* dEps = nullptr;            // eigenvalues (n)
    double* dE   = nullptr;            // rocSOLVER workspace (n)
    rocblas_int* dInfo = nullptr;      // convergence flag
    // FP32 mixed-precision (X-AP3): allocated lazily on the first fp32 solve.
    float*  dF32 = nullptr; float* dS32 = nullptr; float* dEps32 = nullptr; float* dE32 = nullptr;
    int     fp32_n = 0;

    bool ensureFp32(int n) {
        if (n == fp32_n && dF32) return true;
        for (void* p : { (void*)dF32, (void*)dS32, (void*)dEps32, (void*)dE32 }) if (p) hipFree(p);
        dF32 = dS32 = dEps32 = dE32 = nullptr; fp32_n = 0;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dF32, nn * sizeof(float)) != hipSuccess) return false;
        if (hipMalloc(&dS32, nn * sizeof(float)) != hipSuccess) return false;
        if (hipMalloc(&dEps32, n * sizeof(float)) != hipSuccess) return false;
        if (hipMalloc(&dE32, n * sizeof(float)) != hipSuccess) return false;
        fp32_n = n;
        return true;
    }
    // ----- cached Cholesky factor of the resident overlap (Claude Generated, Oct 2026) -----
    // sygvd = potrf(S) + sygst + syevd + trsm. S (dS0) is fixed for one geometry, so the
    // factor L (S = L L^T, lower) is computed once per S and the per-iteration solve runs
    // sygst -> syevd -> trsm with it. The cache key is (n, s_version); s_version is bumped
    // whenever dS0 is (re)written or reallocated. computeIntegrals already factors dS0 into
    // dL; that factor is reused when it belongs to the current S (dL_ver == s_version), else
    // the private buffer dLc holds it. dLc32 is the FP32 copy of the FP64 factor.
    bool    opt_cached_chol  = false;   ///< -xtb.rocm_cached_cholesky (default false since Oct 17, 2026)
    bool    opt_density_syrk = false;   ///< -xtb.rocm_density_syrk (default false)
    unsigned long long s_version = 1;   ///< bumped on every write/realloc of dS0
    unsigned long long dL_ver = 0;      ///< s_version for which dL = chol(dS0) is valid (0 = none)
    double* dLc = nullptr;   int chol_n = 0;   unsigned long long chol_ver = 0;
    float*  dLc32 = nullptr; int chol32_n = 0; unsigned long long chol32_ver = 0;
    const double* curL = nullptr;       ///< FP64 factor of the current S (dL or dLc), set by cholFp64
    bool    chol_alloc_failed = false;  ///< allocation failed once: stay on sygvd (no retry per step)
    // ----- explicit inverse factor Linv = L^-1 (Claude Generated, Oct 2026) -----
    // -xtb.rocm_explicit_linv (default false): Linv is formed once per S (rocsolver_dtrtri on a copy of
    // the cached FP64 factor; dLinv32 is its FP32 rounding) and each step reduces with two
    // rocblas_trmm, F' = Linv F Linv^T, and back-transforms with C = Linv^T Y (trmm) instead of
    // sygst/trsm. Same cache key (n, s_version) as the factor. Any failure -> sygst/trsm path.
    bool    opt_explicit_linv = false;  ///< -xtb.rocm_explicit_linv (default false since Oct 17, 2026)
    double* dLinv = nullptr;   int linv_n = 0;   unsigned long long linv_ver = 0;
    float*  dLinv32 = nullptr; int linv32_n = 0; unsigned long long linv32_ver = 0;
    bool    linv_alloc_failed = false;  ///< allocation failed once: no retry per step
    unsigned long long linv_bad_ver = 0; ///< s_version for which trtri reported a singular L
    unsigned long long chol_bad_ver = 0; ///< s_version for which potrf failed (no retry per step)
    unsigned long long chol_chk_ver = 0; ///< s_version for which the dL from computeIntegrals passed factorCheck
    // ----- factor check (Claude Generated, Oct 2026) -----
    // rocSOLVER potrf/sygvd were observed to return a wrong result in about 0.05 percent of the
    // calls when several processes share the GPU (docs/ROCM_GFN_POTENTIAL_2026-10-14.md,
    // scripts/rocsolver_concurrency_check.cpp). With the cached factor one wrong L would persist
    // through the whole SCF, so every new L (and Linv) is probed once per geometry with a fixed
    // vector: |L L^T v - S v| and |L (Linv v) - v| must stay at rounding level, else the factor is
    // recomputed (twice at most) and the sygvd path takes over if it stays wrong.
    double* dChk = nullptr; int chk_n = 0;   ///< 4 n doubles: v, t, y, s
    int     chk_force_fail = 0;              ///< fault injection: number of checks still forced to fail
    // ----- fallback events (Claude Generated, Oct 2026) -----
    // Bits of XtbHipContext::SolverFallbackEvent, read by the method wrapper, which prints each
    // new bit once (this TU has no logger). fail_inject is the diagnostic fault injection of
    // CURCUMA_ROCM_FAIL_LINV (0 none, 1 alloc, 2 trtri, 3 chol), read once per process.
    unsigned fallback_events = 0;
    int      fail_inject = failInjectMode();
    static int failInjectMode() {
        static const int mode = [] {
            const char* e = std::getenv("CURCUMA_ROCM_FAIL_LINV");
            if (!e) return 0;
            const std::string v(e);
            return v == "alloc" ? 1 : v == "trtri" ? 2 : v == "chol" ? 3
                 : v == "factor" ? 4 : v == "factorall" ? 5 : 0;
        }();
        return mode;
    }

    void freeCholCache() {
        if (dChk) (void)hipFree(dChk);
        dChk = nullptr; chk_n = 0; chol_chk_ver = 0;
        if (dLc) (void)hipFree(dLc);
        if (dLc32) (void)hipFree(dLc32);
        dLc = nullptr; dLc32 = nullptr; chol_n = chol32_n = 0; chol_ver = chol32_ver = 0; curL = nullptr;
        if (dLinv) (void)hipFree(dLinv);
        if (dLinv32) (void)hipFree(dLinv32);
        dLinv = nullptr; dLinv32 = nullptr; linv_n = linv32_n = 0; linv_ver = linv32_ver = 0;
    }
    /**
     * @brief Probe a factor: L L^T v = S v (Linv == nullptr) or L (Linv v) = v (Linv given).
     * @param L the lower Cholesky factor of dS0
     * @param Linv its inverse, or nullptr to test the factor itself
     * @return true when the probe agrees to 1e-10 relative; Claude Generated (Oct 2026).
     * Returns true when the probe buffers cannot be allocated (the check is skipped, not failed).
     * CURCUMA_ROCM_FAIL_LINV=factor makes the first check of a process fail, =factorall all of them.
     * Claude Generated (Oct 2026).
     */
    bool factorCheck(const double* L, const double* Linv, int n) {
        if (fail_inject == 4 || fail_inject == 5) {
            if (fail_inject == 5 || chk_force_fail == 0) { if (fail_inject == 4) chk_force_fail = -1; return false; }
        }
        if (!dChk || chk_n != n) {
            if (dChk) (void)hipFree(dChk);
            dChk = nullptr; chk_n = 0;
            if (hipMalloc(&dChk, 4 * static_cast<size_t>(n) * sizeof(double)) != hipSuccess) {
                dChk = nullptr; (void)hipGetLastError(); return true;
            }
            std::vector<double> v(n);
            for (int i = 0; i < n; ++i) v[i] = 1.0 + 0.5 * std::sin(1.7 * i + 0.3);
            if (hipMemcpy(dChk, v.data(), n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
                (void)hipFree(dChk); dChk = nullptr; (void)hipGetLastError(); return true;
            }
            chk_n = n;
        }
        double* dv = dChk; double* dt = dChk + n; double* dy = dChk + 2 * static_cast<size_t>(n); double* ds = dChk + 3 * static_cast<size_t>(n);
        const dim3 b(64), g((n + 63) / 64);
        if (!Linv) {
            hipLaunchKernelGGL(k_chk_ltv, g, b, 0, stream, L, dv, dt, n);
            hipLaunchKernelGGL(k_chk_lv, g, b, 0, stream, L, dt, dy, n);
            hipLaunchKernelGGL(k_chk_sv, g, b, 0, stream, dS0, dv, ds, n);
        } else {
            hipLaunchKernelGGL(k_chk_lv, g, b, 0, stream, Linv, dv, dt, n);
            hipLaunchKernelGGL(k_chk_lv, g, b, 0, stream, L, dt, dy, n);
            hipLaunchKernelGGL(k_chk_lv, g, b, 0, stream, L, dv, ds, n);   // unused scratch (keeps one code path)
        }
        if (hipStreamSynchronize(stream) != hipSuccess) return true;
        std::vector<double> y(n), t(n), v(n);
        if (hipMemcpy(y.data(), dy, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return true;
        if (hipMemcpy(t.data(), ds, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return true;
        if (hipMemcpy(v.data(), dv, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return true;
        const std::vector<double>& ref = Linv ? v : t;   // S v for the factor, v itself for L (Linv v)
        double err = 0.0, sc = 1.0;
        for (int i = 0; i < n; ++i) { err = std::max(err, std::fabs(y[i] - ref[i])); sc = std::max(sc, std::fabs(ref[i])); }
        return err <= 1.0e-10 * sc;
    }
    /// Recompute dLc = chol(dS0) once; true when potrf succeeded with info 0.
    bool factorDLc(int n) {
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpyAsync(dLc, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (rocsolver_dpotrf(handle, rocblas_fill_lower, n, dLc, n, dInfo) != rocblas_status_success) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return info == 0;
    }
    /**
     * @brief Lower Cholesky factor L of the resident overlap dS0 (S = L L^T), FP64.
     * @param n matrix dimension (must equal res_n)
     * @return device pointer to L (column-major, lower triangle valid) or nullptr when the
     *         cache is off or allocation/factorisation failed; the caller then uses sygvd.
     * Called before dF is touched, so a nullptr return leaves the sygvd fallback intact.
     */
    const double* cholFp64(int n) {
        curL = nullptr;
        if (!opt_cached_chol || chol_alloc_failed || n <= 0 || n != res_n || !dS0) return nullptr;
        if (chol_bad_ver == s_version) return nullptr;   // potrf failed for this S already
        if (fail_inject == 3) {   // CURCUMA_ROCM_FAIL_LINV=chol: behave as if potrf failed
            chol_bad_ver = s_version; fallback_events |= XtbHipContext::kFallbackCholFailed; return nullptr;
        }
        if (dL && bnao == n && dL_ver == s_version) {   // from computeIntegrals
            if (chol_chk_ver == s_version) { curL = dL; return curL; }
            if (factorCheck(dL, nullptr, n)) { chol_chk_ver = s_version; curL = dL; return curL; }
            dL_ver = 0; fallback_events |= XtbHipContext::kFallbackFactorRetried;   // recompute into dLc below
        }
        if (dLc && chol_n == n && chol_ver == s_version) { curL = dLc; return curL; }
        const size_t nn = static_cast<size_t>(n) * n;
        if (chol_n != n || !dLc) {
            if (dLc) (void)hipFree(dLc);
            dLc = nullptr; chol_n = 0; chol_ver = 0;
            if (hipMalloc(&dLc, nn * sizeof(double)) != hipSuccess) {
                dLc = nullptr; chol_alloc_failed = true; (void)hipGetLastError();
                fallback_events |= XtbHipContext::kFallbackCholFailed; return nullptr;
            }
            chol_n = n;
        }
        rocblas_int info = 0;
        if (hipMemcpyAsync(dLc, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return nullptr;
        if (rocsolver_dpotrf(handle, rocblas_fill_lower, n, dLc, n, dInfo) != rocblas_status_success) return nullptr;
        if (hipStreamSynchronize(stream) != hipSuccess) return nullptr;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return nullptr;
        if (info != 0) {   // S not positive definite: sygvd reports it the usual way
            chol_bad_ver = s_version; fallback_events |= XtbHipContext::kFallbackCholFailed; return nullptr;
        }
        // factor check (see factorCheck): recompute a wrong factor at most twice, then give up
        for (int attempt = 0; !factorCheck(dLc, nullptr, n); ++attempt) {
            fallback_events |= XtbHipContext::kFallbackFactorRetried;
            if (attempt >= 2 || !factorDLc(n)) {
                chol_bad_ver = s_version; fallback_events |= XtbHipContext::kFallbackCholFailed; return nullptr;
            }
        }
        chol_ver = s_version;
        curL = dLc;
        return curL;
    }
    /**
     * @brief FP32 copy of the FP64 factor from cholFp64 (rounded once per S), or nullptr.
     */
    const float* cholFp32(int n) {
        const double* L = cholFp64(n);
        if (!L) return nullptr;
        if (dLc32 && chol32_n == n && chol32_ver == s_version) return dLc32;
        const size_t nn = static_cast<size_t>(n) * n;
        if (chol32_n != n || !dLc32) {
            if (dLc32) (void)hipFree(dLc32);
            dLc32 = nullptr; chol32_n = 0; chol32_ver = 0;
            if (hipMalloc(&dLc32, nn * sizeof(float)) != hipSuccess) {
                dLc32 = nullptr; (void)hipGetLastError();
                fallback_events |= XtbHipContext::kFallbackCholFailed; return nullptr;
            }
            chol32_n = n;
        }
        const int blk = 256; const int g = static_cast<int>((nn + blk - 1) / blk);
        hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, L, dLc32, nn);
        chol32_ver = s_version;
        return dLc32;
    }

    /**
     * @brief Explicit inverse Linv = L^-1 of the cached FP64 factor of dS0 (lower triangle valid).
     * @param n matrix dimension (must equal res_n)
     * @return device pointer or nullptr when the switch is off, the factor is unavailable, the
     *         buffer cannot be allocated or trtri reports a singular L; the caller then keeps the
     *         sygst/trsm path. Claude Generated (Oct 2026).
     */
    const double* linvFp64(int n) {
        if (!opt_explicit_linv || linv_alloc_failed || linv_bad_ver == s_version) return nullptr;
        const double* L = cholFp64(n);
        if (!L) return nullptr;
        if (dLinv && linv_n == n && linv_ver == s_version) return dLinv;
        const size_t nn = static_cast<size_t>(n) * n;
        if (linv_n != n || !dLinv) {
            if (dLinv) (void)hipFree(dLinv);
            dLinv = nullptr; linv_n = 0; linv_ver = 0;
            if (fail_inject == 1 || hipMalloc(&dLinv, nn * sizeof(double)) != hipSuccess) {
                dLinv = nullptr; linv_alloc_failed = true; (void)hipGetLastError();
                fallback_events |= XtbHipContext::kFallbackLinvAlloc; return nullptr;
            }
            linv_n = n;
        }
        for (int attempt = 0;; ++attempt) {   // Linv = trtri(L), probed by L (Linv v) = v; one recompute
            if (hipMemcpyAsync(dLinv, L, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return nullptr;
            if (rocsolver_dtrtri(handle, rocblas_fill_lower, rocblas_diagonal_non_unit, n, dLinv, n, dInfo)
                != rocblas_status_success) {
                linv_bad_ver = s_version; fallback_events |= XtbHipContext::kFallbackLinvTrtri; return nullptr;
            }
            if (hipStreamSynchronize(stream) != hipSuccess) return nullptr;
            rocblas_int info = 0;
            if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return nullptr;
            if (fail_inject == 2) info = 1;   // CURCUMA_ROCM_FAIL_LINV=trtri: behave as if L were singular
            if (info != 0) {
                linv_bad_ver = s_version; fallback_events |= XtbHipContext::kFallbackLinvTrtri; return nullptr;
            }
            if (factorCheck(L, dLinv, n)) break;
            fallback_events |= XtbHipContext::kFallbackFactorRetried;
            if (attempt >= 1) { linv_bad_ver = s_version; fallback_events |= XtbHipContext::kFallbackLinvTrtri; return nullptr; }
        }
        linv_ver = s_version;
        return dLinv;
    }
    /**
     * @brief FP32 rounding of linvFp64 (once per S), or nullptr. Claude Generated (Oct 2026).
     */
    const float* linvFp32(int n) {
        const double* Li = linvFp64(n);
        if (!Li) return nullptr;
        if (dLinv32 && linv32_n == n && linv32_ver == s_version) return dLinv32;
        const size_t nn = static_cast<size_t>(n) * n;
        if (linv32_n != n || !dLinv32) {
            if (dLinv32) (void)hipFree(dLinv32);
            dLinv32 = nullptr; linv32_n = 0; linv32_ver = 0;
            if (hipMalloc(&dLinv32, nn * sizeof(float)) != hipSuccess) {
                dLinv32 = nullptr; linv_alloc_failed = true; (void)hipGetLastError();
                fallback_events |= XtbHipContext::kFallbackLinvAlloc; return nullptr;
            }
            linv32_n = n;
        }
        const int blk = 256; const int g = static_cast<int>((nn + blk - 1) / blk);
        hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, Li, dLinv32, nn);
        linv32_ver = s_version;
        return dLinv32;
    }

    /**
     * @brief FP64 generalized eigensolve F C = S C eps on the resident matrices.
     *
     * dF holds F on entry and C (ascending eigenvalues) on exit, dEps the eigenvalues.
     * With the cached factor: dsygst (itype 1, lower) -> dsyevd -> dtrsm (C = L^-T Y), the
     * same sequence rocsolver_dsygvd runs after its own potrf. Otherwise dsygvd on a copy
     * of dS0. Claude Generated (Oct 2026, the sygvd branch moved here from resSolve/resSolveMP).
     */
    bool eigensolveFp64(int n) {
        c_exact = true;   // WP7: dF is overwritten by exact eigenvectors (or the step fails)
        const size_t nn = static_cast<size_t>(n) * n;
        const double* Li = linvFp64(n);
        const double* L = Li ? nullptr : cholFp64(n);
        if (Li) {
            // Claude Generated (Oct 2026), -xtb.rocm_explicit_linv: F' = Linv F Linv^T with two
            // out-of-place trmm (dS is the work buffer of the sygvd branch, free here), then
            // syevd on the lower triangle and C = Linv^T Y in place. F is mirrored lower ->
            // upper first so that, like sygst, only the lower triangle of F enters.
            const double one = 1.0;
            dim3 b2(16, 16), g2((n + 15) / 16, (n + 15) / 16);
            hipLaunchKernelGGL(k_mirror_lower, g2, b2, 0, stream, dF, n);
            if (rocblas_dtrmm(handle, rocblas_side_left, rocblas_fill_lower, rocblas_operation_none,
                              rocblas_diagonal_non_unit, n, n, &one, Li, n, dF, n, dS, n) != rocblas_status_success) return false;
            if (rocblas_dtrmm(handle, rocblas_side_right, rocblas_fill_lower, rocblas_operation_transpose,
                              rocblas_diagonal_non_unit, n, n, &one, Li, n, dS, n, dF, n) != rocblas_status_success) return false;
            if (rocsolver_dsyevd(handle, rocblas_evect_original, rocblas_fill_lower, n, dF, n,
                                 dEps, dE, dInfo) != rocblas_status_success) return false;
            if (rocblas_dtrmm(handle, rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose,
                              rocblas_diagonal_non_unit, n, n, &one, Li, n, dF, n, dF, n) != rocblas_status_success) return false;
        } else if (L) {
            double* Lw = const_cast<double*>(L);   // non-const in the API; not written (checked)
            if (rocsolver_dsygst(handle, rocblas_eform_ax, rocblas_fill_lower, n, dF, n, Lw, n)
                != rocblas_status_success) return false;
            if (rocsolver_dsyevd(handle, rocblas_evect_original, rocblas_fill_lower, n, dF, n,
                                 dEps, dE, dInfo) != rocblas_status_success) return false;
            const double one = 1.0;
            if (rocblas_dtrsm(handle, rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose,
                              rocblas_diagonal_non_unit, n, n, &one, Lw, n, dF, n) != rocblas_status_success) return false;
        } else {
            if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
            if (hipStreamSynchronize(stream) != hipSuccess) return false;
            // F C = S C eps on the device; dF <- eigenvectors (ascending), dS destroyed.
            if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                                 n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
        }
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return info == 0;
    }

    // FP32 generalized eigensolve: dF/dS0 (FP64) → ssygvd → dF (C, FP64) + dEps (FP64).
    // dS0 stays intact (resident overlap). Returns false on any error → caller can retry FP64.
    // With the cached factor (opt_cached_chol): ssygst -> ssyevd -> strsm on the FP32 copy
    // of the FP64 factor instead of ssygvd (which would refactor the FP32-rounded S).
    bool eigensolveFp32(int n) {
        c_exact = true;   // WP7: dF is overwritten by exact eigenvectors (or the step fails)
        if (!ensureFp32(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        const int blk = 256; const int g = static_cast<int>((nn + blk - 1) / blk);
        const float* Li32 = linvFp32(n);
        const float* L32 = Li32 ? nullptr : cholFp32(n);
        bool c_in_s32 = false;   // explicit-Linv path leaves C in dS32
        if (Li32) {   // F mirrored lower -> upper in FP64 (see eigensolveFp64), Claude Generated (Oct 2026)
            dim3 b2(16, 16), g2((n + 15) / 16, (n + 15) / 16);
            hipLaunchKernelGGL(k_mirror_lower, g2, b2, 0, stream, dF, n);
        }
        hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, dF, dF32, nn);    // F → FP32
        if (Li32) {
            // -xtb.rocm_explicit_linv: F' = Linv F Linv^T (two out-of-place strmm via dS32),
            // ssyevd, C = Linv^T Y out of place into dS32 (converted to FP64 from there).
            const float onef = 1.0f;
            if (rocblas_strmm(handle, rocblas_side_left, rocblas_fill_lower, rocblas_operation_none,
                              rocblas_diagonal_non_unit, n, n, &onef, Li32, n, dF32, n, dS32, n) != rocblas_status_success) return false;
            if (rocblas_strmm(handle, rocblas_side_right, rocblas_fill_lower, rocblas_operation_transpose,
                              rocblas_diagonal_non_unit, n, n, &onef, Li32, n, dS32, n, dF32, n) != rocblas_status_success) return false;
            if (rocsolver_ssyevd(handle, rocblas_evect_original, rocblas_fill_lower, n, dF32, n,
                                 dEps32, dE32, dInfo) != rocblas_status_success) return false;
            if (rocblas_strmm(handle, rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose,
                              rocblas_diagonal_non_unit, n, n, &onef, Li32, n, dF32, n, dS32, n) != rocblas_status_success) return false;
            c_in_s32 = true;
        } else if (L32) {
            float* Lw = const_cast<float*>(L32);   // rocSOLVER/rocBLAS take non-const B; not written
            if (rocsolver_ssygst(handle, rocblas_eform_ax, rocblas_fill_lower, n, dF32, n, Lw, n)
                != rocblas_status_success) return false;
            if (rocsolver_ssyevd(handle, rocblas_evect_original, rocblas_fill_lower, n, dF32, n,
                                 dEps32, dE32, dInfo) != rocblas_status_success) return false;
            const float onef = 1.0f;   // C = L^-T Y (back-transform, as in sygvd itype 1 lower)
            if (rocblas_strsm(handle, rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose,
                              rocblas_diagonal_non_unit, n, n, &onef, Lw, n, dF32, n) != rocblas_status_success) return false;
        } else {
            hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, dS0, dS32, nn);   // S → FP32 (work)
            if (rocsolver_ssygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                                 n, dF32, n, dS32, n, dEps32, dE32, dInfo) != rocblas_status_success) return false;
        }
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;
        const int gn = static_cast<int>((nn + blk - 1) / blk), gn1 = (n + blk - 1) / blk;
        hipLaunchKernelGGL(k_f2d, dim3(gn), dim3(blk), 0, stream, c_in_s32 ? dS32 : dF32, dF, nn);   // C → FP64 (dF)
        hipLaunchKernelGGL(k_f2d, dim3(gn1), dim3(blk), 0, stream, dEps32, dEps, n); // eps → FP64
        return hipStreamSynchronize(stream) == hipSuccess;
    }

    bool ensureSolver(int n) {
        if (!handle_ok) {
            if (rocblas_create_handle(&handle) != rocblas_status_success) return false;
            rocblas_set_stream(handle, stream);
            handle_ok = true;
        }
        if (n == solver_n && dF) return true;
        freeSolver();
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dF, nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dS, nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dEps, n * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dE, n * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dInfo, sizeof(rocblas_int)) != hipSuccess) return false;
        solver_n = n;
        return true;
    }
    void freeSolver() {
        for (void* p : { (void*)dF, (void*)dS, (void*)dEps, (void*)dE, (void*)dInfo })
            if (p) hipFree(p);
        dF = dS = dEps = dE = nullptr; dInfo = nullptr; solver_n = 0;
        for (void* p : { (void*)dF32, (void*)dS32, (void*)dEps32, (void*)dE32 }) if (p) hipFree(p);
        dF32 = dS32 = dEps32 = dE32 = nullptr; fp32_n = 0;
    }

    // ----- resident GFN1 SCF (Stage 2) ------------------------------------
    // Reuses dF (eigenvectors C), dS (Cholesky work copy), dEps/dE/dInfo from the
    // solver cache; adds the geometry-resident H0/S + density work buffers.
    int     res_n = 0;
    double* dH0 = nullptr;   // resident H0 (n×n)
    double* dS0 = nullptr;   // resident S  (n×n, untouched; dS is the work copy)
    double* dvao = nullptr;  // v_ao (n)
    double* dP  = nullptr;   // density (n×n)
    double* dCw = nullptr;   // C·diag(occ) (n×n)
    double* dOcc = nullptr;  // occupations (n)
    double* dPB = nullptr;   // pop + bandpart (2n)

    bool ensureResident(int n) {
        if (!ensureSolver(n)) return false;
        if (n == res_n && dH0) return true;
        freeResident();
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dH0,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dS0,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dvao, n  * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dP,   nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dCw,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dOcc, n  * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dPB,  2 * n * sizeof(double)) != hipSuccess) return false;
        res_n = n;
        return true;
    }
    void freeResident() {
        for (void* p : { (void*)dH0, (void*)dS0, (void*)dvao, (void*)dP, (void*)dCw, (void*)dOcc, (void*)dPB })
            if (p) hipFree(p);
        dH0 = dS0 = dvao = dP = dCw = dOcc = dPB = nullptr; res_n = 0;
        freeCholCache();   // the factor belongs to dS0 (Claude Generated, Oct 2026)
        ++s_version;
    }

    bool resBegin(const double* H0, const double* S, int n) {
        if (!ok || n <= 0 || !ensureResident(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dH0, H0, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        ++s_version;   // new S: invalidates the cached Cholesky factor (Claude Generated, Oct 2026)
        if (hipMemcpy(dS0, S,  nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        return true;
    }
    // fp32=true (X-AP3): solve this step in FP32 (ssygvd, or the cached-factor chain) — far-
    // from-convergence -scf_mixed_precision iterations. dF holds C (FP64) + dEps the eigenvalues
    // either way.
    bool resSolve(const double* v_ao, double* eps_out, bool fp32) {
        const int n = res_n; if (n <= 0) return false;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);   // F → dF
        if (fp32) {
            if (!eigensolveFp32(n)) return false;   // dF ← C (FP64), dEps ← eps (FP64)
        } else {
            if (!eigensolveFp64(n)) return false;   // dF <- C (ascending), dEps <- eps
        }
        if (hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;   // dF holds C resident, ascending order
    }
    bool resDensity(const double* occ, int ncol, double* pop_ao, double* band) {
        const int n = res_n; if (n <= 0 || ncol < 0 || ncol > n) return false;
        // Claude Generated (Oct 2026), opt-in -xtb.rocm_density_syrk: P = Cw Cw^T with
        // Cw = C(:, :ncol) diag(sqrt(occ)) via rocblas_dsyrk (lower triangle, n^2 ncol flops
        // instead of 2 n^2 ncol), then mirrored to the full matrix the kernels read. Needs
        // occ >= 0; any negative occupation keeps the GEMM path.
        bool syrk = opt_density_syrk && ncol > 0;
        for (int r = 0; syrk && r < ncol; ++r) if (!(occ[r] >= 0.0)) syrk = false;
        std::vector<double> occ_full(n, 0.0);          // ascending order matches C columns
        for (int r = 0; r < ncol; ++r) occ_full[r] = syrk ? std::sqrt(occ[r]) : occ[r];
        if (hipMemcpy(dOcc, occ_full.data(), n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_scale_cols, grd, blk, 0, stream, dF, dOcc, dCw, n);   // Cw = C·diag(occ or sqrt(occ))
        const double one = 1.0, zero = 0.0;
        // P = Cw · Cᵀ
        // Claude Generated (Oct 2026): only the first ncol columns carry weight, so the inner
        // dimension is ncol instead of n (same sum, the dropped terms were exact zeros).
        if (syrk) {
            if (rocblas_dsyrk(handle, rocblas_fill_lower, rocblas_operation_none, n, ncol,
                              &one, dCw, n, &zero, dP, n) != rocblas_status_success) return false;
            hipLaunchKernelGGL(k_mirror_lower, grd, blk, 0, stream, dP, n);   // upper <- lower
        } else if (ncol > 0) {
            if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                              n, n, ncol, &one, dCw, n, dF, n, &zero, dP, n) != rocblas_status_success) return false;
        } else if (hipMemsetAsync(dP, 0, static_cast<size_t>(n) * n * sizeof(double), stream) != hipSuccess) return false;
        dim3 blk1(64), grd1((n + 63) / 64);
        hipLaunchKernelGGL(k_popband, grd1, blk1, 0, stream, dP, dS0, dH0, dPB, n);
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        std::vector<double> pb(static_cast<size_t>(2) * n);
        if (hipMemcpy(pb.data(), dPB, 2 * n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        double b = 0.0;
        for (int i = 0; i < n; ++i) { pop_ao[i] = pb[i]; b += pb[n + i]; }
        *band = b;
        return true;
    }
    bool resFinalize(double* P, double* C) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(P, dP, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(C, dF, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;
    }

    // ----- Stage 3: device integral build ---------------------------------
    bool basis_ready = false;
    int  bnat = 0, bnsh = 0, bnao = 0, bnprim = 0, bis_gfn2 = 0, bgeom_nat = 0;
    int*    dz = nullptr; int* dsh2at = nullptr; int* dang_sh = nullptr; int* diao_sh = nullptr;
    int*    dnao_sh = nullptr; int* dsh_nprim = nullptr; int* dsh_prim_off = nullptr; int* dvalence = nullptr;
    double* dprim_alpha = nullptr; double* dprim_coeff = nullptr; double* dsh_zeta = nullptr;
    double* dselfenergy = nullptr; double* dkcn = nullptr; double* dshpoly = nullptr; double* dghard = nullptr;
    double* dcovrad = nullptr; double* dpauling = nullptr; double* darad = nullptr;   // element tables (86)
    double* dxyz = nullptr; double* dcn = nullptr; double* dse = nullptr; double* dgamma = nullptr; double* dL = nullptr;
    // Stage 4 (gradient) buffers:
    int*    dao2sh = nullptr; int* dao2at = nullptr;
    double* drep_alpha = nullptr; double* drep_zeff = nullptr;
    double* dW = nullptr; double* dVao = nullptr; double* dQsh = nullptr;
    double* dGrad = nullptr; double* dEdcn_buf = nullptr;
    // Stage 3m (GFN2 multipole integrals): dp_int 3·nn / qp_int 6·nn.
    double* dDpInt = nullptr; double* dQpInt = nullptr;
    bool    mp_computed = false;
    // Stage 2b (resident multipole SCF): per-iteration v_dp/v_qp (up) / dp_at/qp_at (down).
    double* dVdp = nullptr; double* dVqp = nullptr; double* dDpAt = nullptr; double* dQpAt = nullptr;

    // ----- Stage 5 (Part B2): in-SCF GFN2 D4 dE_D4/dq on device ------------
    // Reference data uploaded once per geometry (begDispersion); the element C6 block
    // (geometry- and molecule-independent) only once per process (d4_c6_uploaded).
    // dispDedq then runs the per-iteration O(N²) contraction + BJ disp_sum.
    static constexpr int D4_MAX_REF = 7;
    int     d4_nat = 0, d4_c6_len = 0;
    double  d4_s6 = 0.0, d4_s8 = 0.0, d4_a1 = 0.0, d4_a2 = 0.0, d4_cut = 0.0;
    bool    d4_c6_uploaded = false;
    int*    dD4Z = nullptr;     // nat
    int*    dD4Nref = nullptr;  // nat
    double* dD4Sqrt = nullptr;  // nat (√(½·Z·⟨r⁴⟩/⟨r²⟩))
    double* dD4Xyz = nullptr;   // 3·nat (Bohr)
    double* dD4W = nullptr;     // nat·D4_MAX_REF (per-atom reference weights)
    double* dD4dWq = nullptr;   // nat·D4_MAX_REF (∂W/∂q)
    double* dD4Dedq = nullptr;  // nat (output)
    double* dD4C6 = nullptr;    // c6_flat (118²·7² element references)
    // 2-body D4 gradient (dispDispGrad): the extra ∂W/∂CN weights + the gather outputs.
    double* dD4dWc = nullptr;   // nat·D4_MAX_REF (∂W/∂CN)
    double* dD4Eat = nullptr;   // nat (per-atom energy; host halves+sums)
    double* dD4Grad = nullptr;  // 3·nat (nuclear gradient, [3*i+k])
    double* dD4Dcn = nullptr;   // nat (dE/dCN)
    // ATM 3-body (dispATM): the q=0 reference C6 / ∂C6/∂CN matrices (reuses dD4Eat/Grad/Dcn out).
    double* dD4c6 = nullptr;    // nat·nat (C6, symmetric)
    double* dD4dc6 = nullptr;   // nat·nat (∂C6(a,b)/∂CN_a)

    // ----- Single-shot D4 EEQ (Stage 5, Part A): charges + ∂q/∂x response -----
    int     eeq_n = 0;          // cached N (LU valid for the response when ==N)
    double* dEeqXyz = nullptr;     // 3N
    double* dEeqChi = nullptr;     // N
    double* dEeqGam = nullptr;     // N
    double* dEeqAlp = nullptr;     // N
    double* dEeqCnf = nullptr;     // N
    double* dEeqRcov = nullptr;    // N
    double* dEeqCn = nullptr;      // N
    double* dEeqCnRaw = nullptr;   // N
    double* dEeqM = nullptr;       // m·m (LU factor, reused by the adjoint)
    double* dEeqRhs = nullptr;     // m
    double* dEeqQ = nullptr;       // N (EEQ charges)
    rocblas_int* dEeqIpiv = nullptr; // m (LU pivots)
    double* dEeqDedq = nullptr;    // N
    double* dEeqAdjRhs = nullptr;  // m
    double* dEeqU = nullptr;       // N
    double* dEeqGrad = nullptr;    // 3N

    template <class T> bool up(T*& dptr, const T* host, size_t count) {
        if (dptr) hipFree(dptr);
        dptr = nullptr;
        if (hipMalloc(&dptr, count * sizeof(T)) != hipSuccess) return false;
        return hipMemcpy(dptr, host, count * sizeof(T), hipMemcpyHostToDevice) == hipSuccess;
    }

    bool begBasis(const XtbHipBasisData& bd) {
        if (!ok || bd.nao <= 0 || !ensureResident(bd.nao)) return false;  // allocates dH0/dS0/...
        bnat = bd.nat; bnsh = bd.nsh; bnao = bd.nao; bnprim = bd.nprim_total; bis_gfn2 = bd.is_gfn2;
        ++s_version; dL_ver = 0;   // dS0/dL contents are stale until computeIntegrals (Oct 2026)
        bool g = true;
        g &= up(dz, bd.z, bnat);
        g &= up(dsh2at, bd.sh2at, bnsh);    g &= up(dang_sh, bd.ang_sh, bnsh);
        g &= up(diao_sh, bd.iao_sh, bnsh);  g &= up(dnao_sh, bd.nao_sh, bnsh);
        g &= up(dsh_nprim, bd.sh_nprim, bnsh); g &= up(dsh_prim_off, bd.sh_prim_off, bnsh);
        if (bd.valence) g &= up(dvalence, bd.valence, bnsh);   // GFN1 only; null for GFN2
        else { if (dvalence) hipFree(dvalence); dvalence = nullptr; }
        g &= up(dprim_alpha, bd.prim_alpha, bnprim); g &= up(dprim_coeff, bd.prim_coeff, bnprim);
        g &= up(dsh_zeta, bd.sh_zeta, bnsh); g &= up(dselfenergy, bd.selfenergy, bnsh);
        g &= up(dkcn, bd.kcn, bnsh);         g &= up(dshpoly, bd.shpoly, bnsh);
        g &= up(dghard, bd.shell_hardness, bnsh);
        loop_ready = false;                   // WP2: a new basis ends any resident loop
        g &= buildAtomCsr(bd);
        double covrad[86], pauling[86], arad[86];
        for (int z = 1; z <= 86; ++z) {
            covrad[z - 1]  = curcuma::xtb::covalent_rad_d3_au(z);
            pauling[z - 1] = curcuma::xtb::pauling_en[z - 1];
            arad[z - 1]    = curcuma::xtb::atomic_rad_au(z);
        }
        g &= up(dcovrad, covrad, 86); g &= up(dpauling, pauling, 86); g &= up(darad, arad, 86);
        // Stage 4 gradient basis arrays (skip gracefully if not provided).
        if (bd.ao2sh) g &= up(dao2sh, bd.ao2sh, bnao);
        if (bd.ao2at) g &= up(dao2at, bd.ao2at, bnao);
        if (bd.rep_alpha) g &= up(drep_alpha, bd.rep_alpha, bnat);
        if (bd.rep_zeff)  g &= up(drep_zeff, bd.rep_zeff, bnat);
        if (bnat != bgeom_nat || !dgamma) {
            for (void* p : { (void*)dxyz, (void*)dcn, (void*)dse, (void*)dgamma, (void*)dL,
                             (void*)dW, (void*)dVao, (void*)dQsh, (void*)dGrad, (void*)dEdcn_buf }) if (p) hipFree(p);
            dxyz = dcn = dse = dgamma = dL = dW = dVao = dQsh = dGrad = dEdcn_buf = nullptr;
            g &= hipMalloc(&dxyz, 3 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dcn, static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dse, static_cast<size_t>(bnsh) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dgamma, static_cast<size_t>(bnsh) * bnsh * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dL, static_cast<size_t>(bnao) * bnao * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dW, static_cast<size_t>(bnao) * bnao * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dVao, static_cast<size_t>(bnao) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dQsh, static_cast<size_t>(bnsh) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dGrad, 3 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dEdcn_buf, static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            bgeom_nat = bnat;
        }
        // Stage 3m + 2b: GFN2 multipole-integral buffers + per-iteration potentials/moments.
        for (void* p : { (void*)dDpInt, (void*)dQpInt, (void*)dVdp, (void*)dVqp, (void*)dDpAt, (void*)dQpAt })
            if (p) hipFree(p);
        dDpInt = dQpInt = dVdp = dVqp = dDpAt = dQpAt = nullptr;
        mp_computed = false;
        mp_pot_valid = false;   // dVdp/dVqp reallocated (Claude Generated, Oct 2026)
        if (bis_gfn2) {
            const size_t nn = static_cast<size_t>(bnao) * bnao;
            g &= hipMalloc(&dDpInt, 3 * nn * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dQpInt, 6 * nn * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dVdp, 3 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dVqp, 6 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dDpAt, 3 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
            g &= hipMalloc(&dQpAt, 6 * static_cast<size_t>(bnat) * sizeof(double)) == hipSuccess;
        }
        basis_ready = g;
        return g;
    }

    /// WP8b (Claude Generated, Oct 2026): resident = true takes v_ao (dVao), q_sh (dQsh) and for GFN2 v_dp/v_qp
    /// (dVdp/dVqp) from convPotentials() instead of the host arrays (which may then be nullptr).
    bool gradient(const double* eps, int nocc_orbs, const double* v_ao, const double* q_sh,
                  const double* v_dp, const double* v_qp,
                  double* grad_out, double* dEdcn_out, bool resident = false) {
        if (!basis_ready || bnao <= 0 || !dP) return false;
        if (resident && !conv_pot_valid) { grad_error = "converged potentials are not on the device"; return false; }
        if (!resident) conv_pot_valid = false;   // dVao/dQsh/dVdp/dVqp are overwritten below
        const int nao = bnao, nat = bnat, nsh = bnsh;
        if (!resident) {
        if (hipMemcpy(dVao, v_ao, nao * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dQsh, q_sh, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        }
        // GFN2 (R-AP3): upload the converged multipole potentials for the multipole-integral
        // Pulay term inside k_grad_h0_pulay (column-major 3×nat / 6×nat, matching v_dp[k+iat*3]).
        const bool mp_grad = bis_gfn2 && dVdp && dVqp && (resident || (v_dp && v_qp));
        if (mp_grad && !resident) {
            if (hipMemcpy(dVdp, v_dp, 3 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
            if (hipMemcpy(dVqp, v_qp, 6 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        }
        const double one = 1.0, zero = 0.0;
        // Claude Generated (Oct 2026): inner dimension = number of weighted columns (rest are zero).
        int kocc = std::min(nocc_orbs, nao);
        // WP7 (Claude Generated, Oct 2026; port of the CUDA w_from_goo branch): after a pseudo-diagonalisation
        // step (-scf_pseudo_diag_fp64 on the converged step) the columns of C span the occupied space but do
        // not diagonalise the occupied block of F, so C diag(2 eps) C^T misses the off-diagonal G_ij. The
        // basis-invariant W = C_o (2 G_oo) C_o^T with G_oo = C_o^T F C_o is exact for integer occupations (the
        // pseudo step's gap gate); F is rebuilt into dS from the converged potentials uploaded above.
        const bool w_from_goo = !c_exact && ps_no_last > 0 && ps_no_last <= nao;
        if (w_from_goo) {
            const int no = ps_no_last;
            if (!dS || solver_n != nao || res_n != nao || !dH0 || !dS0 || (bis_gfn2 && (!mp_grad || !dDpInt || !dQpInt)))
                return false;
            const double two = 2.0;
            dim3 b(16, 16), gr((nao + 15) / 16, (nao + 15) / 16);
            hipLaunchKernelGGL(k_fock, gr, b, 0, stream, dH0, dS0, dVao, dS, nao);
            if (bis_gfn2)
                hipLaunchKernelGGL(k_add_fock_multipole, gr, b, 0, stream, dS, dDpInt, dQpInt, dVdp, dVqp, dao2at, nao);
            if (hipGetLastError() != hipSuccess
                || rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_none, nao, no, nao, &one, dS, nao,
                                 dF, nao, &zero, dW, nao) != rocblas_status_success           // T = F C_o (dW scratch)
                || rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, no, no, nao, &one, dF, nao,
                                 dW, nao, &zero, dS, no) != rocblas_status_success            // G_oo (dS scratch)
                || rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_none, nao, no, no, &two, dF, nao,
                                 dS, no, &zero, dCw, nao) != rocblas_status_success)          // Cw = C_o 2 G_oo
                return false;
            kocc = no;
        } else {
        // Energy-weighted density W = C·diag(2·ε_occ)·Cᵀ (reuse dOcc/dCw; dF = C resident).
        std::vector<double> occ2(nao, 0.0);
        for (int k = 0; k < nocc_orbs && k < nao; ++k) occ2[k] = 2.0 * eps[k];
        if (hipMemcpy(dOcc, occ2.data(), nao * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        { dim3 b(16, 16), gr((nao + 15) / 16, (nao + 15) / 16);
          hipLaunchKernelGGL(k_scale_cols, gr, b, 0, stream, dF, dOcc, dCw, nao); }
        }
        if (kocc > 0) {
            if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                              nao, nao, kocc, &one, dCw, nao, dF, nao, &zero, dW, nao) != rocblas_status_success) return false;
        } else if (hipMemsetAsync(dW, 0, static_cast<size_t>(nao) * nao * sizeof(double), stream) != hipSuccess) return false;
        if (hipMemsetAsync(dGrad, 0, 3 * nat * sizeof(double), stream) != hipSuccess) return false;
        if (hipMemsetAsync(dEdcn_buf, 0, nat * sizeof(double), stream) != hipSuccess) return false;
        const double kexp = 1.5, rexp = 1.0, kexp_light = bis_gfn2 ? 1.0 : 1.5;
        const int fixi = fixed_point ? 1 : 0;   // WP10: -xtb.rocm_fixed_point
        { dim3 b(128), gr((nat + 127) / 128);
          hipLaunchKernelGGL(k_grad_repulsion, gr, b, 0, stream, nat, bis_gfn2, dz, dxyz,
                             drep_alpha, drep_zeff, kexp, rexp, kexp_light, dGrad, fixi); }
        { dim3 b(128), gr((nao + 127) / 128);
          hipLaunchKernelGGL(k_grad_cn_onsite, gr, b, 0, stream, nao, dP, dkcn, dao2sh, dao2at, dEdcn_buf, fixi); }
        { dim3 b(16, 16), gr((nao + 15) / 16, (nao + 15) / 16);
          hipLaunchKernelGGL(k_grad_h0_pulay, gr, b, 0, stream, nao, bis_gfn2,
                             dao2sh, dao2at, dang_sh, diao_sh, dsh_nprim, dsh_prim_off,
                             dprim_alpha, dprim_coeff, dsh_zeta, dshpoly, dkcn, dvalence, dz, dse,
                             dxyz, dP, dS0, dH0, dW, dVao,
                             mp_grad ? dVdp : nullptr, mp_grad ? dVqp : nullptr,
                             darad, dpauling, dGrad, dEdcn_buf, fixi); }
        { dim3 b(128), gr((nsh + 127) / 128);
          hipLaunchKernelGGL(k_grad_coulomb, gr, b, 0, stream, nsh, bis_gfn2, dsh2at, dghard,
                             dQsh, dxyz, 2.0, dGrad, fixi); }
        if (fixed_point) {   // WP10: integer accumulators -> doubles before the download
            hipLaunchKernelGGL(k_xfix_to_double, dim3((3 * nat + 255) / 256), dim3(256), 0, stream, 3 * nat, dGrad);
            hipLaunchKernelGGL(k_xfix_to_double, dim3((nat + 255) / 256), dim3(256), 0, stream, nat, dEdcn_buf);
        }
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (fixed_point && !fixFlagClear("gradient")) return false;
        if (hipMemcpy(grad_out, dGrad, 3 * nat * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(dEdcn_out, dEdcn_buf, nat * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;
    }

    bool computeIntegrals(const double* xyz) {
        cn_ok = false;   // WP8a: set again once this build succeeds
        conv_pot_valid = false;   // WP8b: the converged potentials belong to the previous geometry
        if (!basis_ready || bnao <= 0) return false;
        ++s_version; dL_ver = 0;   // dS0 is rebuilt below; dL matches it once potrf succeeds
        if (hipMemcpy(dxyz, xyz, 3 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        { dim3 b(64), gr((bnat + 63) / 64);
          hipLaunchKernelGGL(k_cn, gr, b, 0, stream, bnat, dxyz, dz, bis_gfn2, dcovrad, dcn); }
        { dim3 b(64), gr((bnsh + 63) / 64);
          hipLaunchKernelGGL(k_self_energy, gr, b, 0, stream, bnsh, dselfenergy, dkcn, dsh2at, dcn, dse); }
        { dim3 b(16, 16), gr((bnsh + 15) / 16, (bnsh + 15) / 16);
          hipLaunchKernelGGL(k_overlap_h0, gr, b, 0, stream, bnsh, bnao, bis_gfn2,
                             dsh2at, dang_sh, diao_sh, dnao_sh, dsh_nprim, dsh_prim_off,
                             dprim_alpha, dprim_coeff, dsh_zeta, dshpoly, dse, dz, dvalence, dxyz,
                             darad, dpauling, dS0, dH0); }   // S → dS0, H0 → dH0 (resident)
        { dim3 b(16, 16), gr((bnsh + 15) / 16, (bnsh + 15) / 16);
          hipLaunchKernelGGL(k_gamma, gr, b, 0, stream, bnsh, bis_gfn2, dsh2at, dghard, dxyz, dgamma); }
        // L = chol(S): work on a copy of S0 (resident S must survive).
        if (hipMemcpyAsync(dL, dS0, static_cast<size_t>(bnao) * bnao * sizeof(double),
                           hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (rocsolver_dpotrf(handle, rocblas_fill_lower, bnao, dL, bnao, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;   // dS0/dH0 are now the device-built integrals (resident)
        dL_ver = s_version;            // dL = chol(dS0): reused by the cached-Cholesky solve
        // GFN2: build dp_int/qp_int on the device too (resident overlap dS0 just built),
        // so they are resident for the multipole Fock/moments + downloadable for the host gradient.
        if (bis_gfn2 && !computeMultipole()) return false;
        cn_ok = true;
        return true;
    }

    bool dlMat(double* out, const double* src, size_t count) {
        return hipMemcpy(out, src, count * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
    }

    // Stage 3m: GFN2 AO multipole integrals on device (needs dS0 = device overlap).
    bool computeMultipole() {
        if (!basis_ready || !bis_gfn2 || bnao <= 0 || !dDpInt || !dQpInt) return false;
        dim3 b(16, 16), gr((bnao + 15) / 16, (bnao + 15) / 16);
        hipLaunchKernelGGL(k_multipole_ints, gr, b, 0, stream, bnao,
                           dao2sh, dao2at, diao_sh, dang_sh, dsh_nprim, dsh_prim_off,
                           dprim_alpha, dprim_coeff, dxyz, dS0, dDpInt, dQpInt);
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        mp_computed = true;
        return true;
    }
    bool downloadMP(double* dp3, double* qp6) {
        if (!mp_computed || !dDpInt || !dQpInt) return false;
        const size_t nn = static_cast<size_t>(bnao) * bnao;
        if (hipMemcpy(dp3, dDpInt, 3 * nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(qp6, dQpInt, 6 * nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;
    }

    // ----- Stage 2b: device-resident GFN2 multipole SCF --------------------
    bool beginResidentMP() { return basis_ready && bis_gfn2 && mp_computed && dF && dP; }
    // One resident SCF step with the GFN2 anisotropic Fock term: F = H0 − ½S(v⊕v) (dF),
    // then F −= ½·multipole(v_dp,v_qp); solve F C = S C ε via rocSOLVER (dF ← C resident).
    bool resSolveMP(const double* v_ao, const double* v_dp, const double* v_qp, double* eps_out, bool fp32) {
        conv_pot_valid = false;   // WP8b: dVdp/dVqp are overwritten
        const int n = res_n; if (n <= 0 || !beginResidentMP()) return false;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dVdp, v_dp, 3 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dVqp, v_qp, 6 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);              // F → dF
        hipLaunchKernelGGL(k_add_fock_multipole, grd, blk, 0, stream, dF, dDpInt, dQpInt,
                           dVdp, dVqp, dao2at, n);                                            // F −= ½ multipole
        if (fp32) {
            if (!eigensolveFp32(n)) return false;   // dF ← C (FP64), dEps ← eps (FP64)
        } else {
            if (!eigensolveFp64(n)) return false;   // dF <- C (ascending), dEps <- eps
        }
        return hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
    }
    // Atomic dp_at(3×nat)/qp_at(6×nat) from the resident density dP (k_multipole_moments).
    bool multipoleMomentsImpl(double* dp_at, double* qp_at) {
        if (!beginResidentMP() || bnat <= 0) return false;
        if (hipMemsetAsync(dDpAt, 0, 3 * static_cast<size_t>(bnat) * sizeof(double), stream) != hipSuccess) return false;
        if (hipMemsetAsync(dQpAt, 0, 6 * static_cast<size_t>(bnat) * sizeof(double), stream) != hipSuccess) return false;
        dim3 b(128), gr((bnao + 127) / 128);
        hipLaunchKernelGGL(k_multipole_moments, gr, b, 0, stream, dDpAt, dQpAt, dP, dDpInt, dQpInt, dao2at, bnao,
                           fixed_point ? 1 : 0);
        if (fixed_point) {   // WP10: integer accumulators -> doubles before the download
            hipLaunchKernelGGL(k_xfix_to_double, dim3((3 * bnat + 255) / 256), dim3(256), 0, stream, 3 * bnat, dDpAt);
            hipLaunchKernelGGL(k_xfix_to_double, dim3((6 * bnat + 255) / 256), dim3(256), 0, stream, 6 * bnat, dQpAt);
        }
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (fixed_point && !fixFlagClear("multipole moments")) return false;
        if (hipMemcpy(dp_at, dDpAt, 3 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return hipMemcpy(qp_at, dQpAt, 6 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
    }

    // ----- device multipole potential of the host SCF loop (Claude Generated, Oct 2026) -----
    // -xtb.rocm_device_multipole: the 15 nat^2 interaction matrices (k_mp_potential_stored) are
    // uploaded once per host setupMultipole (the caller passes its version counter), the
    // per-iteration potential writes dVdp/dVqp in place, so solveMultipole needs no upload.
    double* dMpAmat = nullptr;   ///< 15 nat^2 (layout: k_mp_potential_stored)
    double* dMpKern = nullptr;   ///< dkernel (nat) | qkernel (nat)
    double* dMpIn   = nullptr;   ///< q (nat) | dp (3 nat) | qp (6 nat)
    double* dMpOut  = nullptr;   ///< v_at (nat) or the per-atom energy shares (nat)
    int     mp_nat  = 0;
    unsigned long long mp_ver = 0;   ///< host version of the uploaded matrices (0 = none)
    bool    mp_pot_valid = false;    ///< dVdp/dVqp hold a device-built potential (cleared by beginBasis)
    std::string mp_error;            ///< why the last upload failed

    // WP11 (Claude Generated, Oct 2026): on-the-fly mode of the WP1 device potential (beginPotentialOnTheFly). dMpAmat
    // stays unallocated; the potential and energy kernels rebuild the elements from dMpGeo.
    bool    mp_otf  = false;         ///< the multipole interaction of the device potential is rebuilt per call
    double* dMpGeo  = nullptr;       ///< xyz (3 nat, Bohr) | mrad (nat) for the on-the-fly kernels
    double  mp_dmp3 = 3.0, mp_dmp5 = 4.0;

    void freeMpPot() {
        for (void* p : { (void*)dMpAmat, (void*)dMpKern, (void*)dMpIn, (void*)dMpOut, (void*)dMpGeo }) if (p) (void)hipFree(p);
        dMpAmat = dMpKern = dMpIn = dMpOut = dMpGeo = nullptr;
        mp_nat = 0; mp_ver = 0; mp_pot_valid = false; mp_otf = false;
    }
    /// WP11: multipole potential of the device potential (stored or on the fly) from dMpIn into dVdp/dVqp/dMpOut.
    void launchMpPotential(int nat) {
        const size_t un = static_cast<size_t>(nat);
        if (mp_otf)
            hipLaunchKernelGGL(k_mp_potential_otf, dim3(nat), dim3(kMpBlock), 0, stream, nat, dMpGeo, mp_dmp3, mp_dmp5,
                               dMpKern, dMpKern + nat, dMpIn, dMpIn + un, dMpIn + 4 * un, dVdp, dVqp, dMpOut);
        else
            hipLaunchKernelGGL(k_mp_potential_stored, dim3(nat), dim3(kMpBlock), 0, stream, nat,
                               dMpAmat, dMpKern, dMpKern + nat, dMpIn, dMpIn + un, dMpIn + 4 * un,
                               dVdp, dVqp, dMpOut);
    }
    /// WP11: per-atom multipole energy shares (stored or on the fly) for the charges q | dp | qp (10 nat, packed).
    void launchMpEnergy(int nat, const double* qin, double* e_at) {
        const size_t un = static_cast<size_t>(nat);
        if (mp_otf)
            hipLaunchKernelGGL(k_mp_energy_otf, dim3(nat), dim3(kMpBlock), 0, stream, nat, dMpGeo, mp_dmp3, mp_dmp5,
                               dMpKern, dMpKern + nat, qin, qin + un, qin + 4 * un, e_at);
        else
            hipLaunchKernelGGL(k_mp_energy_stored, dim3(nat), dim3(kMpBlock), 0, stream, nat, dMpAmat, dMpKern,
                               dMpKern + nat, qin, qin + un, qin + 4 * un, e_at);
    }
    // Diagnostic fault injection for the device multipole path, CURCUMA_ROCM_FAIL_MP=
    // limit|alloc|upload|potential|energy|solve (read once per process). "potential",
    // "energy" and "solve" fail from the third call on, i.e. in the middle of the SCF.
    static const std::string& mpFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_MP"); return std::string(e ? e : ""); }();
        return m;
    }
    bool mpInject(const char* what) {
        if (mpFailMode() != what) return false;
        return ++mp_inject_calls[what[0] == 'p' ? 0 : what[0] == 'e' ? 1 : 2] >= 3;
    }
    int mp_inject_calls[3] = {0, 0, 0};
    bool mpPotUpload(int nat, const double* const* mats, const double* dk, const double* qk,
                     unsigned long long version, double max_bytes) {
        if (!ok || nat <= 0 || version == 0) return false;
        if (nat == mp_nat && version == mp_ver && dMpAmat) return true;   // already resident
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        const double bytes = 15.0 * static_cast<double>(nn) * sizeof(double);
        if (bytes > max_bytes || mpFailMode() == "limit") {
            freeMpPot();
            mp_error = "interaction matrices need " + std::to_string(bytes / 1.0e9).substr(0, 5)
                     + " GB, above the limit of " + std::to_string(max_bytes / 1.0e9).substr(0, 5) + " GB";
            return false;
        }
        if (nat != mp_nat || !dMpAmat) {
            freeMpPot();
            if (mpFailMode() == "alloc" || hipMalloc(&dMpAmat, 15 * nn * sizeof(double)) != hipSuccess
                || hipMalloc(&dMpKern, 2 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                || hipMalloc(&dMpIn, 10 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                || hipMalloc(&dMpOut, static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) {
                (void)hipGetLastError();
                freeMpPot();
                mp_error = "device allocation of " + std::to_string(bytes / 1.0e9).substr(0, 5) + " GB failed";
                return false;
            }
            mp_nat = nat;
        }
        mp_ver = 0;   // stays 0 until the whole upload succeeded
        if (mpFailMode() == "upload") { mp_error = "matrix upload failed"; return false; }
        for (int b = 0; b < 15; ++b)
            if (hipMemcpy(dMpAmat + static_cast<size_t>(b) * nn, mats[b], nn * sizeof(double),
                          hipMemcpyHostToDevice) != hipSuccess) { mp_error = "matrix upload failed"; return false; }
        if (hipMemcpy(dMpKern, dk, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpKern + nat, qk, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            mp_error = "kernel upload failed"; return false;
        }
        mp_ver = version;
        mp_error.clear();
        return true;
    }
    bool mpUploadMoments(const double* q, const double* dp, const double* qp) {
        const size_t n = static_cast<size_t>(mp_nat);
        return hipMemcpy(dMpIn, q, n * sizeof(double), hipMemcpyHostToDevice) == hipSuccess
            && hipMemcpy(dMpIn + n, dp, 3 * n * sizeof(double), hipMemcpyHostToDevice) == hipSuccess
            && hipMemcpy(dMpIn + 4 * n, qp, 6 * n * sizeof(double), hipMemcpyHostToDevice) == hipSuccess;
    }
    // v_dp/v_qp land in dVdp/dVqp (resident, read by resSolveMPResident) and are also copied
    // to the host, so the host Potential stays complete; v_at is returned.
    bool mpPotential(const double* q, const double* dp, const double* qp,
                     double* v_at, double* v_dp, double* v_qp) {
        mp_pot_valid = false;
        conv_pot_valid = false;   // WP8b: dVdp/dVqp are overwritten
        if (mpInject("potential")) return false;
        if (!dMpAmat || mp_ver == 0 || mp_nat != bnat || !dVdp || !dVqp) return false;
        if (!mpUploadMoments(q, dp, qp)) return false;
        const int nat = mp_nat;
        const size_t n = static_cast<size_t>(nat);
        hipLaunchKernelGGL(k_mp_potential_stored, dim3(nat), dim3(kMpBlock), 0, stream, nat,
                           dMpAmat, dMpKern, dMpKern + nat, dMpIn, dMpIn + n, dMpIn + 4 * n,
                           dVdp, dVqp, dMpOut);
        if (hipGetLastError() != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (hipMemcpy(v_at, dMpOut, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(v_dp, dVdp, 3 * n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(v_qp, dVqp, 6 * n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        mp_pot_valid = true;
        return true;
    }
    bool mpEnergy(const double* q, const double* dp, const double* qp, double* e_out) {
        if (!dMpAmat || mp_ver == 0 || mp_nat <= 0 || mpInject("energy")) return false;
        if (!mpUploadMoments(q, dp, qp)) return false;
        const int nat = mp_nat;
        const size_t n = static_cast<size_t>(nat);
        hipLaunchKernelGGL(k_mp_energy_stored, dim3(nat), dim3(kMpBlock), 0, stream, nat,
                           dMpAmat, dMpKern, dMpKern + nat, dMpIn, dMpIn + n, dMpIn + 4 * n, dMpOut);
        if (hipGetLastError() != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        std::vector<double> e_at(n);
        if (hipMemcpy(e_at.data(), dMpOut, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        double e = 0.0;
        for (double v : e_at) e += v;
        *e_out = e;
        return true;
    }
    // resSolveMP without the v_dp/v_qp upload: uses the potential mpPotential left in dVdp/dVqp.
    bool resSolveMPResident(const double* v_ao, double* eps_out, bool fp32) {
        const int n = res_n; if (n <= 0 || mpInject("solve") || !beginResidentMP() || !mp_pot_valid) return false;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);
        hipLaunchKernelGGL(k_add_fock_multipole, grd, blk, 0, stream, dF, dDpInt, dQpInt,
                           dVdp, dVqp, dao2at, n);
        if (fp32) {
            if (!eigensolveFp32(n)) return false;
        } else {
            if (!eigensolveFp64(n)) return false;
        }
        return hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
    }

    // ----- WP1: GFN2 potential of the host SCF loop on the device (Claude Generated, Oct 2026) -----
    // -xtb.rocm_device_potential. Port of the CUDA beginPotential / beginPotentialStored /
    // residentSolvePotential / buildDevicePotentialOnly. The multipole interaction matrices share
    // dMpAmat/dMpKern with the -xtb.rocm_device_multipole path (15-block layout of
    // k_mp_potential_stored); both paths are never active in the same SCF. Per iteration the host
    // uploads q_sh, q_at, dp_at, qp_at and the D4 reference weights W/dWq; the device builds
    //   v_sh = gamma q_sh + q_sh^2 Gamma_s,  v_at = multipole shift + dE_D4/dq,  v_dp, v_qp,
    // expands v_ao, adds the Fock and solves. FP64 throughout (the eigensolve follows fp32).
    int     pot_nsh = 0;                 ///< > 0 once a beginPotential* succeeded for this basis
    double* dPotQsh = nullptr;           ///< q_sh (nsh)
    double* dVsh    = nullptr;           ///< v_sh (nsh)
    double* dGamma3 = nullptr;           ///< Gamma_s (nsh)
    unsigned long long pot_ver_next = (1ull << 62);   ///< private version range of the matrices (host versions are small)
    std::string pot_error;               ///< why the last begin/solve of this path failed
    int     pot_inject_calls = 0;
    int     pot_alloc_nsh = 0;           ///< size of dPotQsh/dVsh/dGamma3
    /// Diagnostic fault injection CURCUMA_ROCM_FAIL_POT=begin|solve (read once per process): "begin"
    /// fails beginPotential and beginPotentialStored, "solve" fails resSolvePotential from its third
    /// call on (in the middle of the SCF).
    static const std::string& potFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_POT"); return std::string(e ? e : ""); }();
        return m;
    }
    void freePot() {
        for (void* p : { (void*)dPotQsh, (void*)dVsh, (void*)dGamma3 }) if (p) (void)hipFree(p);
        dPotQsh = dVsh = dGamma3 = nullptr; pot_nsh = 0; pot_alloc_nsh = 0;
    }
    /// Shared part of both begin variants: per-shell buffers, Gamma_s, the kernels, the checks that the
    /// resident GFN2 state (integrals, gamma, multipole integrals, D4 reference data) is in place.
    bool potBeginCommon(int nat, int nsh, const double* dk, const double* qk, const double* g3,
                        double max_bytes, bool otf = false) {
        pot_nsh = 0;
        mp_otf = false;   // WP11: set again by potBeginOtf once its begin succeeded
        solv_active = false; solv_failed = false;   // WP6: the Born matrix is per geometry (beginSolvation after this begin)
        if (potFailMode() == "begin") { pot_error = "forced by CURCUMA_ROCM_FAIL_POT=begin"; return false; }
        if (!ok || !basis_ready || !bis_gfn2 || nat != bnat || nsh != bnsh || !dgamma || !dao2sh || !dao2at
            || !beginResidentMP()) { pot_error = "device integrals of this geometry are not resident"; return false; }
        if (nat != d4_nat || !dD4C6 || !dD4W || !dD4dWq || !dD4Dedq) {
            pot_error = "D4 reference data not on the device"; return false;
        }
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        const double bytes = 15.0 * static_cast<double>(nn) * sizeof(double);
        if (otf) {
            // WP11: no stored matrices; a resident dMpAmat of an earlier stored begin is released (the memory is the
            // point of this mode).
            if (nat != mp_nat || dMpAmat || !dMpGeo) {
                freeMpPot();
                if (hipMalloc(&dMpKern, 2 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                    || hipMalloc(&dMpIn, 10 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                    || hipMalloc(&dMpOut, static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                    || hipMalloc(&dMpGeo, 4 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) {
                    (void)hipGetLastError();
                    freeMpPot();
                    pot_error = "device allocation of the on-the-fly multipole buffers failed";
                    return false;
                }
                mp_nat = nat;
            }
        } else if (bytes > max_bytes) {
            pot_error = "multipole interaction matrices need " + std::to_string(bytes / 1.0e9).substr(0, 5)
                      + " GB, above the limit of " + std::to_string(max_bytes / 1.0e9).substr(0, 5) + " GB";
            return false;
        }
        if (!otf && (nat != mp_nat || !dMpAmat)) {
            freeMpPot();
            if (hipMalloc(&dMpAmat, 15 * nn * sizeof(double)) != hipSuccess
                || hipMalloc(&dMpKern, 2 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                || hipMalloc(&dMpIn, 10 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess
                || hipMalloc(&dMpOut, static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) {
                (void)hipGetLastError();
                freeMpPot();
                pot_error = "device allocation of " + std::to_string(bytes / 1.0e9).substr(0, 5) + " GB failed";
                return false;
            }
            mp_nat = nat;
        }
        mp_ver = 0;   // the matrices are rewritten by the caller
        if (!dPotQsh || !dVsh || !dGamma3 || pot_alloc_nsh != nsh) {
            freePot();
            if (hipMalloc(&dPotQsh, nsh * sizeof(double)) != hipSuccess
                || hipMalloc(&dVsh, nsh * sizeof(double)) != hipSuccess
                || hipMalloc(&dGamma3, nsh * sizeof(double)) != hipSuccess) {
                (void)hipGetLastError(); freePot(); pot_error = "device allocation of the shell buffers failed"; return false;
            }
            pot_alloc_nsh = nsh;
        }
        if (hipMemcpy(dMpKern, dk, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpKern + nat, qk, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dGamma3, g3, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            pot_error = "kernel upload failed"; return false;
        }
        return true;
    }
    /// beginPotential: the host interaction matrices (3 nat^2 sd, 9 nat^2 dd with block a*3+b, 6 nat^2 sq).
    bool potBegin(int nat, int nsh, const double* sd, const double* dd, const double* sq,
                  const double* dk, const double* qk, const double* g3, double max_bytes) {
        if (!potBeginCommon(nat, nsh, dk, qk, g3, max_bytes)) return false;
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        const double* src[15] = { sd, sd + nn, sd + 2 * nn,
                                  dd, dd + nn, dd + 2 * nn, dd + 4 * nn, dd + 5 * nn, dd + 8 * nn,
                                  sq, sq + nn, sq + 2 * nn, sq + 3 * nn, sq + 4 * nn, sq + 5 * nn };
        for (int b = 0; b < 15; ++b)
            if (hipMemcpy(dMpAmat + static_cast<size_t>(b) * nn, src[b], nn * sizeof(double),
                          hipMemcpyHostToDevice) != hipSuccess) { pot_error = "matrix upload failed"; return false; }
        mp_ver = pot_ver_next++;
        pot_nsh = nsh;
        return true;
    }
    /// beginPotentialStored: the matrices built on the device from the geometry and the damping radii.
    bool potBeginStored(int nat, int nsh, const double* xyz, const double* mrad, double dmp3, double dmp5,
                        const double* dk, const double* qk, const double* g3, double max_bytes) {
        if (!potBeginCommon(nat, nsh, dk, qk, g3, max_bytes)) return false;
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        // dMpIn (10 nat) is free until the first iteration: geometry (3 nat) and radii (nat) go there.
        if (hipMemcpy(dMpIn, xyz, 3 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpIn + 3 * nat, mrad, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            pot_error = "geometry upload failed"; return false;
        }
        const unsigned blocks = static_cast<unsigned>((nn + 255) / 256);
        hipLaunchKernelGGL(k_mp_build_amat15, dim3(blocks), dim3(256), 0, stream, nat, dMpIn, dMpIn + 3 * nat,
                           dmp3, dmp5, dMpAmat);
        if (hipGetLastError() != hipSuccess || hipStreamSynchronize(stream) != hipSuccess) {
            pot_error = "matrix build kernel failed"; return false;
        }
        mp_ver = pot_ver_next++;
        pot_nsh = nsh;
        return true;
    }
    /// WP11 (Claude Generated, Oct 2026): beginPotentialOnTheFly, port of the CUDA XtbGpuContext::beginPotentialOnTheFly.
    /// Geometry and damping radii stay resident in dMpGeo; the potential and energy kernels rebuild the elements.
    bool potBeginOtf(int nat, int nsh, const double* xyz, const double* mrad, double dmp3, double dmp5,
                     const double* dk, const double* qk, const double* g3) {
        if (!potBeginCommon(nat, nsh, dk, qk, g3, 0.0, true)) return false;
        if (hipMemcpy(dMpGeo, xyz, 3 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpGeo + 3 * nat, mrad, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            pot_error = "geometry upload failed"; return false;
        }
        mp_dmp3 = dmp3;
        mp_dmp5 = dmp5;
        mp_otf = true;
        mp_ver = pot_ver_next++;
        pot_nsh = nsh;
        return true;
    }
    /// One SCF step: potential from the uploaded SCC input -> Fock -> eigensolve; eps downloaded, C in dF.
    bool resSolvePotential(const double* q_sh, const double* q_at, const double* dp, const double* qp,
                           const double* W, const double* dWq, double* eps_out, bool fp32) {
        const int n = res_n, nat = bnat, nsh = pot_nsh;
        conv_pot_valid = false;   // WP8b: dVdp/dVqp are overwritten
        if (potFailMode() == "solve" && ++pot_inject_calls >= 3) {
            pot_error = "forced by CURCUMA_ROCM_FAIL_POT=solve"; return false;
        }
        if (n <= 0 || nsh <= 0 || nsh != bnsh || (!dMpAmat && !mp_otf) || mp_ver < (1ull << 62) || mp_nat != nat
            || nat != d4_nat || !beginResidentMP()) { pot_error = "device potential not set up"; return false; }
        const size_t un = static_cast<size_t>(nat);
        const size_t wlen = un * D4_MAX_REF;
        if (hipMemcpy(dPotQsh, q_sh, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpIn, q_at, un * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpIn + un, dp, 3 * un * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dMpIn + 4 * un, qp, 6 * un * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dD4W, W, wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dD4dWq, dWq, wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            pot_error = "upload of the SCC input failed"; return false;
        }
        if (const char* err = potStepCore(fp32)) { pot_error = err; return false; }
        mp_pot_valid = false;   // dVdp/dVqp now belong to this path
        if (hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) {
            pot_error = "eigenvalue download failed"; return false;
        }
        return true;
    }
    /**
     * @brief Device half of one SCF step from the resident inputs dPotQsh (q_sh), dMpIn (q_at | dp | qp,
     *        atom-major) and dD4W/dD4dWq: D4 dE/dq, v_sh = gamma q_sh + q_sh^2 Gamma_s, multipole v_dp/v_qp and
     *        v_at shift, v_at += D4, v_ao, Fock with the multipole term, eigensolve (dF <- C, dEps <- eps).
     * @return nullptr on success, else a short reason. Shared by resSolvePotential (WP1) and the resident
     *         loop (WP2). Claude Generated (Oct 2026).
     */
    const char* potStepCore(bool fp32) {
        if (const char* err = potBuildCore()) return err;
        return fockAndSolve(true, fp32);
    }
    /**
     * @brief Potential half of potStepCore (D4 dE/dq, v_sh, multipole v_dp/v_qp and shift, solvation, v_ao into
     *        dvao), without the Fock build and the eigensolve. Split off for convPotentials (WP8b); the kernel
     *        sequence of potStepCore is unchanged. Claude Generated (Oct 2026).
     * @return nullptr on success, else a short reason
     */
    const char* potBuildCore() {
        const int n = res_n, nat = bnat, nsh = pot_nsh;
        const int b = 128;
        // D4 dE/dq (resident, no download)
        hipLaunchKernelGGL(k_d4_dedq, dim3((nat + b - 1) / b), dim3(b), 0, stream, nat, 118, D4_MAX_REF,
                           dD4Z, dD4Nref, dD4Sqrt, dD4Xyz, dD4C6, dD4W, dD4dWq,
                           d4_s6, d4_s8, d4_a1, d4_a2, d4_cut * d4_cut, dD4Dedq);
        // v_sh = gamma q_sh (gamma symmetric, column-major nsh x nsh) + third order
        const double one = 1.0, zero = 0.0;
        if (rocblas_dgemv(handle, rocblas_operation_none, nsh, nsh, &one, dgamma, nsh, dPotQsh, 1,
                          &zero, dVsh, 1) != rocblas_status_success) return "gamma q_sh gemv failed";
        hipLaunchKernelGGL(k_vsh_third, dim3((nsh + b - 1) / b), dim3(b), 0, stream, nsh, dPotQsh, dGamma3, dVsh);
        // multipole v_dp/v_qp (resident for the Fock) and the v_at shift, then v_at += D4
        // (WP11: stored matrices or rebuilt on the fly)
        launchMpPotential(nat);
        hipLaunchKernelGGL(k_vat_add_d4, dim3((nat + b - 1) / b), dim3(b), 0, stream, nat, dD4Dedq, dMpOut);
        if (const char* err = solvApply(nat, dMpIn, dMpOut)) return err;   // WP6: v_at += B q_at
        hipLaunchKernelGGL(k_expand_vao, dim3((n + b - 1) / b), dim3(b), 0, stream, n, dVsh, dMpOut,
                           dao2sh, dao2at, dvao);
        return nullptr;
    }
    /**
     * @brief Fock F = H0 - S (v_ao + v_ao^T) / 2 from the resident dvao (plus the GFN2 multipole term from
     *        dVdp/dVqp when multipole is set) and the eigensolve (dF <- C, dEps <- eps).
     * @return nullptr on success, else a short reason. Shared by potStepCore (GFN2) and potStepCoreGfn1.
     *         Claude Generated (Oct 2026, WP3).
     */
    const char* fockAndSolve(bool multipole, bool fp32) {
        const int n = res_n;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);
        if (multipole)
            hipLaunchKernelGGL(k_add_fock_multipole, grd, blk, 0, stream, dF, dDpInt, dQpInt, dVdp, dVqp, dao2at, n);
        if (hipGetLastError() != hipSuccess) return "potential kernel launch failed";
        return eigensolveStep(n, fp32);   // WP7: pseudo step or full solve
    }

    // ----- WP7: pseudo-diagonalisation on the resident loop (Claude Generated, Oct 2026) -----
    // Port of the CUDA XtbGpuContext::setPseudoDiagonalisation / eigensolveResidentFock pseudo branch /
    // pseudoDiagStep / pseudoDiagStep64 / pseudoRotateOrthonormalise (cuda/xtb_gpu_context.cu). Same gates:
    // only inside a resident step (resStep), closed shell, previous max|dq| < pseudo_max_dq, at most
    // pseudo_max_steps in a row, HOMO-LUMO estimate > max(30 kT, 1e-3) and every |theta| < 0.1; FP32 steps,
    // and FP64 steps only with -scf_pseudo_diag_fp64. Differences from CUDA: the previous eigenvectors are
    // stored only by full solves inside the resident loop (the host-driven loop never takes a pseudo step);
    // the FP64 step keeps F intact (Tb = dS, M = dCw), so a failure after the rotation falls back to the full
    // solve instead of failing the SCF step; a failure (rc < 0) or an allocation failure turns the pseudo
    // path off for the rest of the run with one notice (pseudoDiagonalisationNotice). L is the lower Cholesky
    // factor of S of this geometry (the probed dL of the device integral build, else a private factor).
    bool   ps_enabled = false;
    double ps_max_dq = 0.3;
    int    ps_max_steps = 12;
    bool   ps_fp64 = false;
    bool   ps_disabled = false;          ///< a failure turned the path off until the next setPseudoDiagonalisation
    std::string ps_notice;               ///< why (read once by the method wrapper)
    int    ps_consecutive = 0, ps_consecutive64 = 0;
    int    ps_taken = 0, ps_rejected = 0;
    double ps_last_dq = 1e300;           ///< max|dq| of the previous resident step
    bool   ps_in_step = false;           ///< set by resStep around the eigensolve
    bool   cprev_valid = false, cprev64_valid = false;
    bool   c_exact = true;               ///< dF/dEps come from a full eigensolve (false after a pseudo step)
    int    ps_no_last = 0;               ///< occupied block of the last pseudo step (gradient W from G_oo)
    int    ps_inject_calls = 0;
    float*  dCprev32 = nullptr; size_t cap_cprev32 = 0;
    float*  dPsGov = nullptr;   size_t cap_psgov = 0;
    float*  dPsL32 = nullptr;   size_t cap_psl32 = 0; unsigned long long psl32_ver = 0;
    double* dCprev64 = nullptr; size_t cap_cprev64 = 0;
    double* dPsGov64 = nullptr; size_t cap_psgov64 = 0;
    double* dPsG = nullptr;     size_t cap_psg = 0;     ///< Rayleigh quotients (n) + stats (2) at the end
    unsigned int* dPsMax = nullptr;

    /// Diagnostic fault injection CURCUMA_ROCM_FAIL_PSEUDO=begin|step (read once per process): "begin" fails the
    /// first storage of the previous eigenvectors (as an allocation failure), "step" fails a pseudo step after the
    /// rotation from its third attempt on.
    static const std::string& psFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_PSEUDO"); return std::string(e ? e : ""); }();
        return m;
    }
    static bool psTrace() { static const bool t = std::getenv("CURCUMA_PSEUDO_TRACE") != nullptr; return t; }
    /// Keep a device buffer of at least n elements (reallocated only when it grows).
    template <class T> bool psEnsure(T*& p, size_t& cap, size_t n) {
        if (p && cap >= n) return true;
        if (p) (void)hipFree(p);
        p = nullptr; cap = 0;
        if (hipMalloc(&p, std::max<size_t>(n, 1) * sizeof(T)) != hipSuccess) { (void)hipGetLastError(); p = nullptr; return false; }
        cap = n;
        return true;
    }
    /// Per-SCF reset at every loop begin (CUDA beginResidentLoop*): no previous eigenvectors, counts at 0.
    void psResetLoop() {
        cprev_valid = cprev64_valid = false;
        ps_consecutive = ps_consecutive64 = 0;
        ps_taken = ps_rejected = 0;
        ps_last_dq = 1e300;
        c_exact = true;
    }
    void freePseudo() {
        for (void* p : { (void*)dCprev32, (void*)dPsGov, (void*)dPsL32, (void*)dCprev64, (void*)dPsGov64, (void*)dPsG,
                         (void*)dPsMax })
            if (p) (void)hipFree(p);
        dCprev32 = dPsGov = dPsL32 = nullptr; dCprev64 = dPsGov64 = dPsG = nullptr; dPsMax = nullptr;
        cap_cprev32 = cap_psgov = cap_psl32 = cap_cprev64 = cap_psgov64 = cap_psg = 0; psl32_ver = 0;
        cprev_valid = cprev64_valid = false;
        ps_consecutive = ps_consecutive64 = 0;
    }
    /// Turn the path off for the rest of the run (one notice through the wrapper) and release its buffers.
    void psDisable(const std::string& why) {
        if (!ps_disabled && ps_notice.empty()) ps_notice = why;
        ps_disabled = true;
        freePseudo();
    }
    /// Lower Cholesky factor of dS0 (FP64) for the S metric of the pseudo step, independent of
    /// -xtb.rocm_cached_cholesky (cholFp64 with the cache forced on; its fallback bits stay untouched).
    const double* pseudoL64(int n) {
        const bool save = opt_cached_chol;
        const unsigned ev = fallback_events;
        opt_cached_chol = true;
        const double* L = cholFp64(n);
        opt_cached_chol = save;
        if (!save) { fallback_events = ev; chol_alloc_failed = false; }
        return L;
    }
    /// FP32 rounding of pseudoL64, once per S.
    const float* pseudoL32(int n) {
        const double* L = pseudoL64(n);
        if (!L) return nullptr;
        if (dPsL32 && psl32_ver == s_version && cap_psl32 >= static_cast<size_t>(n) * n) return dPsL32;
        const size_t nn = static_cast<size_t>(n) * n;
        if (!psEnsure(dPsL32, cap_psl32, nn)) return nullptr;
        const int blk = 256; const int g = static_cast<int>((nn + blk - 1) / blk);
        hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, L, dPsL32, nn);
        psl32_ver = s_version;
        return dPsL32;
    }
    /**
     * @brief One pseudo-diagonalisation step in precision R (float or double). Port of the CUDA
     *        pseudoRotateOrthonormalise, Claude Generated (Oct 2026).
     * In: F (n x n Fock, AO basis), L (lower factor of S), C (previous eigenvectors, occupied block first).
     * Scratch: Tb (n x n), Gov (nv x no), M (n x n; may alias F, which is only read in step 1).
     * Out: C rotated and S-orthonormalised, g = Rayleigh quotients of the rotated-from vectors (dPsG).
     * @return 1 on success, 0 when rejected before C was touched (gate or BLAS error), -1 after C was modified
     */
    template <typename R>
    int pseudoRotateOrthonormalise(int n, int no, const R* F, const R* L, R* C, R* Tb, R* Gov, R* M) {
        namespace pb = pseudo_blas;
        const int nv = n - no;
        const int b = 256;
        const R one = 1, zero = 0, mone = -1;
        R* Cv = C + static_cast<size_t>(no) * n;
        double* g = dPsG;
        double* stats = dPsG + n;
        // 1. T = F C, g = diag(C^T F C), G_ov = C_v^T (F C_o)
        if (pb::gemm(handle, rocblas_operation_none, rocblas_operation_none, n, n, n, &one, F, n, C, n, &zero, Tb, n)
            != rocblas_status_success) return 0;
        hipLaunchKernelGGL(k_pseudo_coldot<R>, dim3(n), dim3(b), b * sizeof(double), stream, C, Tb, n, g);
        if (pb::gemm(handle, rocblas_operation_transpose, rocblas_operation_none, nv, no, n, &one, Cv, n, Tb, n, &zero,
                     Gov, nv) != rocblas_status_success) return 0;
        // 2. gate: occupied/virtual gap large against kT, every rotation angle small
        hipLaunchKernelGGL(k_pseudo_gap, dim3(1), dim3(b), 2 * b * sizeof(double), stream, g, no, n, stats);
        if (hipMemsetAsync(dPsMax, 0, sizeof(unsigned int), stream) != hipSuccess) return 0;
        hipLaunchKernelGGL(k_pseudo_theta<R>, dim3(1024), dim3(b), 0, stream, Gov, g, no, nv, dPsMax);
        if (hipGetLastError() != hipSuccess) return 0;
        double st[2] = {0.0, 0.0};
        unsigned int maxbits = 0u;
        if (hipMemcpyAsync(st, stats, 2 * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipMemcpyAsync(&maxbits, dPsMax, sizeof(unsigned int), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess) return 0;
        float max_theta = 0.0f;
        std::memcpy(&max_theta, &maxbits, sizeof(float));
        const double kT = loop_Tele * 3.166808e-6;
        if (psTrace())
            fprintf(stderr, "PSTRACE gate: homo %.17g lumo %.17g gap %.17g max_theta %.9g\n", st[0], st[1],
                    st[1] - st[0], static_cast<double>(max_theta));
        if (!(st[1] - st[0] > std::max(30.0 * kT, 1.0e-3)) || !(max_theta < 0.1f)) return 0;
        // 3. rotate: old C_v into Tb, C_v -= C_o theta^T, C_o += C_v(old) theta
        if (hipMemcpyAsync(Tb, Cv, sizeof(R) * static_cast<size_t>(n) * nv, hipMemcpyDeviceToDevice, stream) != hipSuccess
            || pb::gemm(handle, rocblas_operation_none, rocblas_operation_transpose, n, nv, no, &mone, C, n, Gov, nv, &one,
                        Cv, n) != rocblas_status_success
            || pb::gemm(handle, rocblas_operation_none, rocblas_operation_none, n, no, nv, &one, Tb, n, Gov, nv, &one,
                        C, n) != rocblas_status_success)
            return -1;
        if (psFailMode() == "step" && ++ps_inject_calls >= 3) return -1;   // diagnostic fault injection
        // 4. S-orthonormalise: B = L^T C (Tb), M = B^T B (upper, diagonal re-summed in FP64), M = R^T R,
        //    C <- C R^-1 (column i mixes only with columns <= i, so the occupied span is kept)
        if (pb::trmm(handle, n, &one, L, C, Tb) != rocblas_status_success) return -1;
        if (pb::syrk(handle, n, &one, Tb, &zero, M) != rocblas_status_success) return -1;
        hipLaunchKernelGGL(k_pseudo_gram_diag<R>, dim3(n), dim3(b), b * sizeof(double), stream, Tb, n, M);
        if (hipGetLastError() != hipSuccess) return -1;
        if (pb::potrf(handle, n, M, dInfo) != rocblas_status_success) return -1;
        rocblas_int info = 0;
        if (hipMemcpyAsync(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess || info != 0) return -1;
        if (pb::trsm(handle, n, &one, M, C) != rocblas_status_success) return -1;
        return 1;
    }
    /// FP32 pseudo step from dCprev32; on success C (FP64) in dF and the Rayleigh quotients in dEps.
    /// @return 1 taken, 0 rejected (dF = F intact), -1 failed after the rotation (dF = F intact, dCprev32 invalid)
    int pseudoDiagStep(int n) {
        const int no = loop_nocc_pairs, nv = n - no;
        if (no <= 0 || nv <= 0) return 0;
        const size_t nn = static_cast<size_t>(n) * n;
        if (!ensureFp32(n) || !psEnsure(dPsGov, cap_psgov, static_cast<size_t>(nv) * no)
            || !psEnsure(dPsG, cap_psg, static_cast<size_t>(n) + 2) || !dPsMax) return -2;
        const float* L32 = pseudoL32(n);
        if (!L32) return -2;
        const int b = 256, gnn = static_cast<int>((nn + b - 1) / b);
        hipLaunchKernelGGL(k_d2f, dim3(gnn), dim3(b), 0, stream, dF, dF32, nn);   // F -> FP32 (dF keeps F)
        if (hipGetLastError() != hipSuccess) return 0;
        // M aliases the FP32 F (dF32): F is only read in step 1; Tb = dS32.
        const int rc = pseudoRotateOrthonormalise<float>(n, no, dF32, L32, dCprev32, dS32, dPsGov, dF32);
        if (rc < 0) { cprev_valid = false; return -1; }
        if (rc != 1) return 0;
        hipLaunchKernelGGL(k_f2d, dim3(gnn), dim3(b), 0, stream, dCprev32, dF, nn);
        if (hipGetLastError() != hipSuccess
            || hipMemcpyAsync(dEps, dPsG, sizeof(double) * n, hipMemcpyDeviceToDevice, stream) != hipSuccess)
            return -3;   // dF may be half converted: F is lost
        c_exact = false; ps_no_last = no;
        return 1;
    }
    /// FP64 pseudo step from dCprev64 (or the FP32 copy); F stays in dF until success. Return codes as above.
    int pseudoDiagStep64(int n) {
        const int no = loop_nocc_pairs, nv = n - no;
        if (no <= 0 || nv <= 0) return 0;
        const size_t nn = static_cast<size_t>(n) * n;
        if (!dS || !dCw || !psEnsure(dCprev64, cap_cprev64, nn) || !psEnsure(dPsGov64, cap_psgov64, static_cast<size_t>(nv) * no)
            || !psEnsure(dPsG, cap_psg, static_cast<size_t>(n) + 2) || !dPsMax) return -2;
        const double* L = pseudoL64(n);
        if (!L) return -2;
        const int b = 256, gnn = static_cast<int>((nn + b - 1) / b);
        if (!cprev64_valid) {   // start from the FP32 eigenvectors of the previous step
            hipLaunchKernelGGL(k_f2d, dim3(gnn), dim3(b), 0, stream, dCprev32, dCprev64, nn);
            if (hipGetLastError() != hipSuccess) return 0;
            cprev64_valid = true;
        }
        // Tb = dS (work copy of the sygvd branch), M = dCw (density scratch, written later in the step).
        const int rc = pseudoRotateOrthonormalise<double>(n, no, dF, L, dCprev64, dS, dPsGov64, dCw);
        if (rc < 0) { cprev64_valid = false; return -1; }
        if (rc != 1) return 0;
        if (hipMemcpyAsync(dF, dCprev64, sizeof(double) * nn, hipMemcpyDeviceToDevice, stream) != hipSuccess
            || hipMemcpyAsync(dEps, dPsG, sizeof(double) * n, hipMemcpyDeviceToDevice, stream) != hipSuccess)
            return -3;
        if (dCprev32 && cap_cprev32 >= nn) {   // the FP32 copy follows, so a later FP32 step starts here
            hipLaunchKernelGGL(k_d2f, dim3(gnn), dim3(b), 0, stream, dCprev64, dCprev32, nn);
            cprev_valid = hipGetLastError() == hipSuccess;
        }
        c_exact = false; ps_no_last = no;
        return 1;
    }
    /**
     * @brief Eigensolve of one SCF step: the pseudo step when the gates allow it, else the full solve
     *        (eigensolveFp32/Fp64), after which the eigenvectors are stored for the next pseudo step.
     *        dF holds F on entry and C on exit, dEps the eigenvalues (Rayleigh quotients after a pseudo step).
     *        Port of the CUDA eigensolveResidentFock. Claude Generated (Oct 2026, WP7).
     * @return nullptr on success, else a short reason
     */
    const char* eigensolveStep(int n, bool fp32) {
        const bool closed_shell = !loop_open && std::fabs(loop_nelec - 2.0 * loop_nocc_pairs) < 1e-9;
        const bool active = ps_enabled && !ps_disabled && ps_in_step && loop_ready;
        const bool pseudo_ok = active && closed_shell && loop_nao == n && ps_last_dq < ps_max_dq;
        if (psTrace() && ps_in_step)
            fprintf(stderr, "PSTRACE call: fp32 %d pseudo_ok %d cprev %d consecutive %d last_dq %.17g\n", fp32 ? 1 : 0,
                    pseudo_ok ? 1 : 0, cprev_valid ? 1 : 0, ps_consecutive, ps_last_dq);
        // Return codes of pseudoDiagStep*: 1 taken, 0 rejected by the gate, -1 failed after the rotation with F
        // still in dF, -2 buffers not allocated, -3 failed after F was overwritten. Only -3 fails the SCF step
        // (the host then restarts the SCF on the host-driven loop, WP2 fallback); the others run the full solve.
        auto on_failure = [&](int rc, const char* prec) -> const char* {
            if (rc == -2) psDisable(std::string("device memory for the ") + prec
                                    + " pseudo-diagonalisation buffers could not be allocated");
            else psDisable(std::string("an ") + prec + " pseudo-diagonalisation step failed after the rotation"
                           + (psFailMode() == "step" ? " (forced by CURCUMA_ROCM_FAIL_PSEUDO=step)" : ""));
            return rc == -3 ? "pseudo-diagonalisation step failed after the Fock matrix was overwritten" : nullptr;
        };
        if (pseudo_ok && fp32 && cprev_valid && ps_consecutive < ps_max_steps) {
            const int rc = pseudoDiagStep(n);
            if (rc == 1) {
                ++ps_consecutive; ++ps_taken;
                cprev64_valid = false;   // the FP64 copy is now one step behind
                if (psTrace()) fprintf(stderr, "PSTRACE pseudo step taken\n");
                return hipStreamSynchronize(stream) == hipSuccess ? nullptr : "FP32 pseudo step failed";
            }
            if (rc < 0) { if (const char* e = on_failure(rc, "FP32")) return e; }
            else { ++ps_rejected; if (psTrace()) fprintf(stderr, "PSTRACE rejected (gate failed): full solve follows\n"); }
        } else if (pseudo_ok && !fp32 && ps_fp64 && (cprev64_valid || cprev_valid) && ps_consecutive64 < ps_max_steps) {
            const int rc = pseudoDiagStep64(n);
            if (rc == 1) {
                ++ps_consecutive64; ++ps_taken;
                return hipStreamSynchronize(stream) == hipSuccess ? nullptr : "FP64 pseudo step failed";
            }
            if (rc < 0) { if (const char* e = on_failure(rc, "FP64")) return e; }
            else ++ps_rejected;
        }
        if (fp32) {
            if (!eigensolveFp32(n)) return "FP32 eigensolve failed";
        } else {
            if (!eigensolveFp64(n)) return "FP64 eigensolve failed";
        }
        c_exact = true;
        if (active && !ps_disabled) {   // store the eigenvectors for the next pseudo step
            const size_t nn = static_cast<size_t>(n) * n;
            const int b = 256, gnn = static_cast<int>((nn + b - 1) / b);
            const bool inject = psFailMode() == "begin";
            bool have_max = dPsMax != nullptr;
            if (!have_max && !inject) {
                have_max = hipMalloc(&dPsMax, sizeof(unsigned int)) == hipSuccess;
                if (!have_max) { (void)hipGetLastError(); dPsMax = nullptr; }
            }
            if (inject || !have_max || !psEnsure(dCprev32, cap_cprev32, nn)) {
                psDisable(inject ? "the previous eigenvectors were not stored (forced by CURCUMA_ROCM_FAIL_PSEUDO=begin)"
                                 : "device memory for the previous eigenvectors could not be allocated");
            } else {
                hipLaunchKernelGGL(k_d2f, dim3(gnn), dim3(b), 0, stream, dF, dCprev32, nn);
                cprev_valid = hipGetLastError() == hipSuccess;
                ps_consecutive = 0;
                cprev64_valid = false;
                if (ps_fp64 && !fp32) {   // an FP64 solve seeds the FP64 pseudo steps directly
                    cprev64_valid = psEnsure(dCprev64, cap_cprev64, nn)
                        && hipMemcpyAsync(dCprev64, dF, sizeof(double) * nn, hipMemcpyDeviceToDevice, stream) == hipSuccess;
                    ps_consecutive64 = 0;
                }
            }
        }
        return nullptr;
    }
    /**
     * @brief GFN1 device half of one resident SCF step (WP3): v_sh = gamma q_sh from the mixed input dResX,
     *        v_at = Gamma_A q_A^2 (atom-resolved third order, coulomb::potential_third_order_atom) from q_at_in,
     *        v_ao, Fock without multipole term, eigensolve. Port of the CUDA buildDevicePotentialAndSolveGfn1.
     * @param q_at_in device q_at (nat) of the input charges
     * @return nullptr on success, else a short reason. Claude Generated (Oct 2026).
     */
    const char* potStepCoreGfn1(const double* q_at_in, bool fp32) {
        if (const char* err = potBuildGfn1(dResX, q_at_in)) return err;
        return fockAndSolve(false, fp32);
    }
    /**
     * @brief Potential half of potStepCoreGfn1 for the device shell charges q_sh (nsh) and atom charges q_at_in
     *        (nat): v_sh = gamma q_sh, v_at = Gamma_A q_A^2 (+ solvation), v_ao into dvao. Split off for
     *        convPotentials (WP8b); the kernel sequence of potStepCoreGfn1 is unchanged. Claude Generated (Oct 2026).
     * @return nullptr on success, else a short reason
     */
    const char* potBuildGfn1(const double* q_sh, const double* q_at_in) {
        const int n = res_n, nat = loop_nat, nsh = loop_nsh;
        const int b = 128;
        const double one = 1.0, zero = 0.0;
        if (rocblas_dgemv(handle, rocblas_operation_none, nsh, nsh, &one, dgamma, nsh, q_sh, 1,
                          &zero, dResVsh, 1) != rocblas_status_success) return "gamma q_sh gemv failed";
        if (hipMemsetAsync(dResVat, 0, static_cast<size_t>(nat) * sizeof(double), stream) != hipSuccess)
            return "v_at reset failed";
        hipLaunchKernelGGL(k_vsh_third, dim3((nat + b - 1) / b), dim3(b), 0, stream, nat, q_at_in, dResG3at, dResVat);
        if (const char* err = solvApply(nat, q_at_in, dResVat)) return err;   // WP6: v_at += B q_at + B cm5
        hipLaunchKernelGGL(k_expand_vao, dim3((n + b - 1) / b), dim3(b), 0, stream, n, dResVsh, dResVat,
                           dao2sh, dao2at, dvao);
        return nullptr;
    }

    // ----- WP6: implicit solvation (native ALPB/GBSA) on the device (Claude Generated, Oct 2026) -----
    // -xtb.rocm_device_solvation. Port of the CUDA beginSolvation and its use in buildDevicePotentialOnly: the
    // nat x nat Born matrix B of ALPBSolvation (keps-scaled, incl. Born self-energy, HB diagonal and ALPB shape
    // term; symmetric) is uploaded once per geometry, and the device potential adds the in-SCF reaction field
    //   v_at += B q_at (+ v_shift),   v_shift = B cm5 for GFN1 (q_solute = q_Mulliken + cm5), absent for GFN2,
    // which is ALPBSolvation::addPotential (alpb_solvation.cpp) with the same q_at the host loop uses. The
    // solvation energy E = 1/2 q_s^T B q_s + G_sasa + G_shift is not part of the SCC energy (neither on the CPU),
    // so the host evaluates it post-SCF from the converged charges. GFN2: potStepCore (WP1 host-driven step and
    // WP2 resident loop); GFN1: potStepCoreGfn1 (WP3 resident loop).
    double* dSolvB = nullptr;            ///< Born matrix B (nat x nat, column-major == row-major)
    double* dSolvShift = nullptr;        ///< constant v_at shift (nat), GFN1 B cm5
    int     solv_nat = 0;                ///< size dSolvB / dSolvShift were allocated for
    bool    solv_active = false;         ///< B of the current geometry is resident and applied
    bool    solv_shift = false;          ///< dSolvShift is applied
    bool    solv_failed = false;         ///< a solvation step failed; sticky until the next begin
    int     solv_inject_calls = 0;
    std::string solv_error;              ///< why the last solvation begin/step failed
    /// Diagnostic fault injection CURCUMA_ROCM_FAIL_SOLV=begin|step (read once per process): "begin" fails
    /// beginSolvation, "step" fails the solvation part of the potential build from its third call on.
    static const std::string& solvFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_SOLV"); return std::string(e ? e : ""); }();
        return m;
    }
    void freeSolv() {
        for (void* p : { (void*)dSolvB, (void*)dSolvShift }) if (p) (void)hipFree(p);
        dSolvB = dSolvShift = nullptr; solv_nat = 0; solv_active = false; solv_shift = false;
    }
    /**
     * @brief Upload B (and the optional shift) for the current geometry and switch the reaction field on.
     * @param max_bytes limit for B (nat^2 doubles); above it the begin fails and the host keeps solvation
     * @return false (solv_error says why) on a limit, allocation or copy failure; the field stays off
     * Claude Generated (Oct 2026, WP6).
     */
    bool solvBegin(int nat, const double* born, const double* shift, double max_bytes) {
        solv_active = false; solv_shift = false; solv_failed = false;
        solv_error.clear();
        if (solvFailMode() == "begin") { solv_error = "forced by CURCUMA_ROCM_FAIL_SOLV=begin"; return false; }
        if (!ok || nat <= 0 || !born || nat != bnat) { solv_error = "atom count does not match the device basis"; return false; }
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        const double bytes = static_cast<double>(nn) * sizeof(double);
        if (bytes > max_bytes) {
            solv_error = "Born matrix needs " + std::to_string(bytes / 1.0e9).substr(0, 5) + " GB, above the limit of "
                       + std::to_string(max_bytes / 1.0e9).substr(0, 5) + " GB";
            return false;
        }
        if (nat != solv_nat || !dSolvB || !dSolvShift) {
            freeSolv();
            if (hipMalloc(&dSolvB, nn * sizeof(double)) != hipSuccess
                || hipMalloc(&dSolvShift, static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) {
                (void)hipGetLastError(); freeSolv();
                solv_error = "device allocation of the Born matrix (" + std::to_string(bytes / 1.0e9).substr(0, 5)
                           + " GB) failed";
                return false;
            }
            solv_nat = nat;
        }
        if (hipMemcpy(dSolvB, born, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || (shift && hipMemcpy(dSolvShift, shift, static_cast<size_t>(nat) * sizeof(double),
                                   hipMemcpyHostToDevice) != hipSuccess)) {
            solv_error = "Born matrix upload failed"; return false;
        }
        solv_shift = (shift != nullptr);
        solv_active = true;
        return true;
    }
    /**
     * @brief v_at += B q_at (+ shift) on the stream; no-op when no solvent is active.
     * @return nullptr on success, else a short reason (also in solv_error; solv_failed is set)
     * Claude Generated (Oct 2026, WP6).
     */
    const char* solvApply(int nat, const double* q_at, double* v_at) {
        if (!solv_active) return nullptr;
        auto fail = [&](const char* why) -> const char* { solv_failed = true; solv_error = why; return why; };
        if (solv_failed) return "device solvation failed earlier for this geometry";
        if (nat != solv_nat || !dSolvB) return fail("device solvation not set up for this geometry");
        if (solvFailMode() == "step" && ++solv_inject_calls >= 3) return fail("forced by CURCUMA_ROCM_FAIL_SOLV=step");
        const double one = 1.0;
        if (rocblas_dgemv(handle, rocblas_operation_none, nat, nat, &one, dSolvB, nat, q_at, 1, &one, v_at, 1)
            != rocblas_status_success) return fail("Born matrix gemv failed");
        if (solv_shift && rocblas_daxpy(handle, nat, &one, dSolvShift, 1, v_at, 1) != rocblas_status_success)
            return fail("solvation shift axpy failed");
        return nullptr;
    }

    // ----- WP2: fully device-resident GFN2 SCF loop, closed shell (Claude Generated, Oct 2026) -----
    // -xtb.rocm_resident_loop. Port of the CUDA beginResidentLoop / residentScfStep / residentLoopCharges /
    // beginDispersionWeights / runBroydenUpdate (cuda/xtb_gpu_context.cu, Stage 6) onto the WP1 potential
    // (potStepCore). One step: unpack the mixed SCC vector -> D4 weights at the input charges -> potential,
    // Fock, eigensolve -> occupations -> density -> output charges and moments -> SCC energies -> max|dq|
    // -> Broyden. Per step the host receives ncol (one int) and five scalars. Differences from CUDA: see the
    // kernel block above (k_res_*).
    int*    dAtShOff = nullptr;          ///< atom -> first shell (nat + 1), from beginBasis
    int*    dAtAoOff = nullptr;          ///< atom -> first AO (nat + 1), from beginBasis
    bool    at_csr_ok = false;           ///< shells and AOs are atom-contiguous in ascending atom order
    // D4 reference tables (beginDispersionWeights)
    int     rw_nat = 0;
    double* dRwCn = nullptr; double* dRwGi = nullptr; double* dRwZeff = nullptr;
    double* dRwRefcn = nullptr; double* dRwRefcov = nullptr; double* dRwRefq = nullptr;
    int*    dRwNref = nullptr;
    // loop state (beginResidentLoop)
    int     loop_nsh = 0, loop_nat = 0, loop_nao = 0, loop_N = 0, loop_nocc_pairs = 0, loop_step = 0;
    double  loop_Tele = 0.0, loop_nelec = 0.0;
    bool    loop_ready = false;
    bool    loop_gfn1 = false;           ///< WP3: GFN1 loop (SCC vector = q_sh, no D4 weights, no moments)
    size_t  loop_cap[5] = {0, 0, 0, 0, 0};   ///< sizes the loop buffers were allocated for (N, nsh, nat, nao, max_hist)
    double* dResX = nullptr;             ///< mixed SCC input (N, host packing)
    double* dResXout = nullptr;          ///< output SCC vector (N)
    double* dResXnext = nullptr;         ///< Broyden result (N)
    double* dResQat0 = nullptr;          ///< q_at of the initial guess (D4 weights of step 0)
    double* dResN0sh = nullptr; double* dResN0at = nullptr;
    double* dResNsh = nullptr;           ///< shell populations n_sh
    double* dResQsh = nullptr;           ///< output q_sh
    double* dResOutMp = nullptr;         ///< output q_at | dp | qp (atom-major, 10 nat)
    double* dResPart = nullptr;          ///< per-AO moment partials (9 nao)
    double* dResEat = nullptr;           ///< per-atom multipole energy shares (nat)
    double* dResTmp = nullptr;           ///< scratch (max(nsh, nao))
    double* dResScal = nullptr;          ///< dq, band, q.gamma.q, E3, Emp
    int*    dResNcol = nullptr;
    double* dResVsh = nullptr;           ///< GFN1 v_sh (nsh)
    double* dResVat = nullptr;           ///< GFN1 v_at (nat)
    double* dResG3at = nullptr;          ///< GFN1 Gamma_A (nat)
    double* dResQin = nullptr;           ///< GFN1 input q_at (nat)
    // WP4 open shell (setResidentOpenShell): two spin fillings over the same eigenvectors
    bool    loop_open = false;
    double  loop_nalpha = 0.0, loop_nbeta = 0.0;
    double* dResOccA = nullptr; double* dResOccB = nullptr;   ///< alpha / beta occupations (nao)
    double* dResDf = nullptr;            ///< Fermi derivative scratch of the Newton step (nao)
    // device Broyden (BroydenMixer of broyden_mixer.h)
    int     b_maxh = 0, b_iter = 0, b_push = 0;
    double  b_alpha = 0.0, b_w0 = 0.0;
    double* dBF = nullptr; double* dBFlast = nullptr; double* dBVinLast = nullptr; double* dBDFtmp = nullptr;
    double* dBDF = nullptr; double* dBU = nullptr;          ///< history columns (N x maxh each)
    double* dBGram = nullptr; double* dBC = nullptr; double* dBGamma = nullptr;
    std::string res_error;               ///< why the last begin/step of the loop failed
    int     res_inject_calls = 0;
    // WP5 (Claude Generated, Oct 2026): Broyden seed preconditioner M (broydenSetPreconditioner). The blocks
    // stay resident under pc_tag across SCFs; pc_active, the guard and b_step are per SCF (reset by a begin).
    int     b_step = 0;                  ///< Broyden updates since the begin (guard index n)
    unsigned long long pc_tag = 0;
    bool    pc_active = false, pc_keep = false, pc_best = false;
    int     pc_N = 0, pc_nblocks = 0, pc_maxdim = 0, pc_fired_at = -1, pc_scratch_N = 0;
    XtbHipContext::BroydenGuard pc_guard;
    double* dPcFb = nullptr; double* dPcB = nullptr; int* dPcIoff = nullptr; int* dPcBoff = nullptr;
    int*    dPcIdx = nullptr;
    double* dPcTmp = nullptr; double* dPcXbest = nullptr; double* dPcFbest = nullptr;   ///< scratch (pc_scratch_N)
    void freePc() {
        for (void* p : { (void*)dPcFb, (void*)dPcB, (void*)dPcIoff, (void*)dPcBoff, (void*)dPcIdx, (void*)dPcTmp,
                         (void*)dPcXbest, (void*)dPcFbest })
            if (p) (void)hipFree(p);
        dPcFb = dPcB = dPcTmp = dPcXbest = dPcFbest = nullptr; dPcIoff = dPcBoff = dPcIdx = nullptr;
        pc_tag = 0; pc_active = false; pc_N = pc_nblocks = pc_maxdim = pc_scratch_N = 0; pc_fired_at = -1;
        pc_guard = nullptr;
    }
    /// Upload or reactivate M for the SCF begun last (see XtbHipContext::broydenSetPreconditioner).
    bool pcSet(const XtbHipContext::BroydenPrecond* pc, unsigned long long tag, XtbHipContext::BroydenGuard guard,
               bool keep, bool best) {
        pc_active = false;
        res_error.clear();
        if (tag == 0) { res_error = "preconditioner tag 0"; return false; }
        if (!loop_ready || loop_N <= 0) { res_error = "device Broyden not initialised"; return false; }
        const int N = loop_N;
        if (pc) {
            if (pc->N != N || pc->nblocks < 0 || !pc->fallback
                || (pc->nblocks > 0 && (!pc->idx_off || !pc->idx || !pc->b_off || !pc->B))) {
                res_error = "preconditioner size does not match the device Broyden vector"; return false;
            }
            // The gathered block vector lives in shared memory (64 KB LDS on gfx1200; 6144 doubles as on CUDA).
            if (pc->max_dim > 6144) {
                res_error = "a fragment block exceeds 6144 SCC entries (shared-memory limit of the device kernel)";
                return false;
            }
            pc_tag = 0;   // invalid until the upload below completed
            const int nidx = pc->nblocks > 0 ? pc->idx_off[pc->nblocks] : 0;
            const int nB = pc->nblocks > 0 ? pc->b_off[pc->nblocks] : 0;
            bool g = grow(dPcFb, N) && hipMemcpy(dPcFb, pc->fallback, N * sizeof(double), hipMemcpyHostToDevice) == hipSuccess;
            if (g && pc->nblocks > 0)
                g = grow(dPcIoff, pc->nblocks + 1) && grow(dPcBoff, pc->nblocks + 1) && grow(dPcIdx, nidx) && grow(dPcB, nB)
                    && hipMemcpy(dPcIoff, pc->idx_off, (pc->nblocks + 1) * sizeof(int), hipMemcpyHostToDevice) == hipSuccess
                    && hipMemcpy(dPcBoff, pc->b_off, (pc->nblocks + 1) * sizeof(int), hipMemcpyHostToDevice) == hipSuccess
                    && hipMemcpy(dPcIdx, pc->idx, static_cast<size_t>(nidx) * sizeof(int), hipMemcpyHostToDevice) == hipSuccess
                    && hipMemcpy(dPcB, pc->B, static_cast<size_t>(nB) * sizeof(double), hipMemcpyHostToDevice) == hipSuccess;
            if (!g) { res_error = "preconditioner upload failed"; return false; }
            pc_N = N; pc_nblocks = pc->nblocks; pc_maxdim = pc->max_dim;
            pc_tag = tag;
        } else if (pc_tag != tag || pc_N != N) {
            res_error = "no resident preconditioner with this tag"; return false;
        }
        if (pc_scratch_N != N) {
            pc_scratch_N = 0;
            if (!grow(dPcTmp, N) || !grow(dPcXbest, N) || !grow(dPcFbest, N)) {
                res_error = "preconditioner scratch allocation failed"; return false;
            }
            pc_scratch_N = N;
        }
        pc_guard = std::move(guard);
        pc_keep = keep; pc_best = best;
        pc_fired_at = -1;
        pc_active = true;
        return true;
    }
    /// dst = M src (device vectors of length pc_N).
    bool pcApply(const double* src, double* dst) {
        const int N = pc_N, b1 = 256;
        hipLaunchKernelGGL(k_res_pc_diag, dim3((N + b1 - 1) / b1), dim3(b1), 0, stream, dst, dPcFb, src, N);
        if (pc_nblocks > 0) {
            const size_t smem = sizeof(double) * static_cast<size_t>(std::max(1, pc_maxdim));
            hipLaunchKernelGGL(k_res_pc_blocks, dim3(pc_nblocks), dim3(128), smem, stream, dst, src, dPcB, dPcBoff,
                               dPcIoff, dPcIdx);
        }
        return hipGetLastError() == hipSuccess;
    }
    /// Guard hit: drop M as BroydenMixer::dropStepHook(keep) (keep: u_k rewritten to the scalar form, the ring
    /// order is irrelevant; otherwise reset, so the update in progress becomes the first one).
    bool pcDrop(bool keep) {
        if (!pc_active) return true;
        const int N = loop_N;
        if (keep) {
            const int M = (b_push < b_maxh) ? b_push : b_maxh;
            const int b1 = 256, grid = (N + b1 - 1) / b1;
            for (int k = 0; k < M; ++k) {
                double* df = dBDF + static_cast<size_t>(k) * N;
                double* u = dBU + static_cast<size_t>(k) * N;
                if (!pcApply(df, dPcTmp)) return false;
                hipLaunchKernelGGL(k_res_broyden_u_unprecond, dim3(grid), dim3(b1), 0, stream, u, df, dPcTmp, b_alpha, N);
                if (hipGetLastError() != hipSuccess) return false;
            }
        } else {
            b_iter = 1;
            b_push = 0;
        }
        pc_active = false;
        return true;
    }

    /// Diagnostic fault injection CURCUMA_ROCM_FAIL_RES=begin|step (read once per process): "begin" fails
    /// beginResidentLoop, "step" fails residentScfStep from its third call on.
    static const std::string& resFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_RES"); return std::string(e ? e : ""); }();
        return m;
    }
    void freeResLoop() {
        freeResLoopBuffers();
        freePc();   // WP5
        b_maxh = 0;
        for (void* p : { (void*)dRwCn, (void*)dRwGi, (void*)dRwZeff, (void*)dRwRefcn, (void*)dRwRefcov,
                         (void*)dRwRefq, (void*)dRwNref, (void*)dAtShOff, (void*)dAtAoOff })
            if (p) (void)hipFree(p);
        dRwCn = dRwGi = dRwZeff = dRwRefcn = dRwRefcov = dRwRefq = nullptr; dRwNref = nullptr; rw_nat = 0;
        dAtShOff = dAtAoOff = nullptr; at_csr_ok = false;
    }
    /// Atom offsets of the shells and AOs (once per basis, from beginBasis); false if they are not contiguous.
    bool buildAtomCsr(const XtbHipBasisData& bd) {
        at_csr_ok = false;
        if (!bd.sh2at || !bd.ao2at || bd.nat <= 0) return true;   // no resident loop for this basis, not an error
        std::vector<int> sh(bd.nat + 1, 0), ao(bd.nat + 1, 0);
        for (int s = 0; s < bd.nsh; ++s) {
            if (bd.sh2at[s] < 0 || bd.sh2at[s] >= bd.nat || (s > 0 && bd.sh2at[s] < bd.sh2at[s - 1])) return true;
            ++sh[bd.sh2at[s] + 1];
        }
        for (int m = 0; m < bd.nao; ++m) {
            if (bd.ao2at[m] < 0 || bd.ao2at[m] >= bd.nat || (m > 0 && bd.ao2at[m] < bd.ao2at[m - 1])) return true;
            ++ao[bd.ao2at[m] + 1];
        }
        for (int a = 0; a < bd.nat; ++a) { sh[a + 1] += sh[a]; ao[a + 1] += ao[a]; }
        if (!up(dAtShOff, sh.data(), sh.size()) || !up(dAtAoOff, ao.data(), ao.size())) return false;
        at_csr_ok = true;
        return true;
    }
    template <class T> bool grow(T*& p, size_t n) {
        if (p) (void)hipFree(p);
        p = nullptr;
        if (hipMalloc(&p, std::max<size_t>(n, 1) * sizeof(T)) != hipSuccess) { (void)hipGetLastError(); p = nullptr; return false; }
        return true;
    }
    bool resBeginWeights(int nat, const double* cn, const double* gi, const double* zeff, const double* refcn,
                         const double* refcov, const double* refq, const int* nref) {
        rw_nat = 0;
        rw_geom_ok = false;
        const size_t un = static_cast<size_t>(nat), wl = un * D4_MAX_REF;
        if (!ok || nat <= 0 || !up(dRwCn, cn, un) || !up(dRwGi, gi, un) || !up(dRwZeff, zeff, un)
            || !up(dRwNref, nref, un) || !up(dRwRefcn, refcn, wl) || !up(dRwRefcov, refcov, wl)
            || !up(dRwRefq, refq, wl)) { res_error = "upload of the D4 reference tables failed"; return false; }
        rw_nat = nat;
        rw_geom_ok = true;   // WP8a: tables of the geometry of the last begDispersion
        return true;
    }
    bool resBeginLoop(int nsh, int nat, int nao, double Tele, double n_elec, int nocc_pairs,
                      const double* q_sh0, const double* dp0, const double* qp0, const double* q_at0,
                      const double* n0_sh, const double* n0_at, double alpha, int max_hist, double w0) {
        loop_ready = false;
        conv_pot_valid = false;   // WP8b
        res_error.clear();
        if (resFailMode() == "begin") { res_error = "forced by CURCUMA_ROCM_FAIL_RES=begin"; return false; }
        if (nao != res_n || nsh != pot_nsh || nat != bnat || nat != d4_nat || nat != rw_nat || (!dMpAmat && !mp_otf)
            || mp_ver < (1ull << 62) || mp_nat != nat || !beginResidentMP()) {
            res_error = "device potential of this geometry is not set up"; return false;
        }
        if (!at_csr_ok) { res_error = "shells or AOs are not atom-contiguous"; return false; }
        if (max_hist <= 0 || max_hist > 20) { res_error = "Broyden history outside 1..20"; return false; }
        const int N = nsh + 9 * nat;
        if (!resAllocLoop(N, nsh, nat, nao, max_hist)) return false;
        // X = packSCC(q_sh0, dp0, qp0) in the host order (dp0/qp0 are column-major 3 x nat / 6 x nat).
        std::vector<double> X(N);
        for (int s = 0; s < nsh; ++s) X[s] = q_sh0[s];
        for (int k = 0; k < 3; ++k) for (int a = 0; a < nat; ++a) X[nsh + k * nat + a] = dp0[k + 3 * a];
        for (int k = 0; k < 6; ++k) for (int a = 0; a < nat; ++a) X[nsh + 3 * nat + k * nat + a] = qp0[k + 6 * a];
        if (hipMemcpy(dResX, X.data(), N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResQat0, q_at0, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResN0sh, n0_sh, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResN0at, n0_at, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            res_error = "upload of the initial guess failed"; return false;
        }
        loop_nsh = nsh; loop_nat = nat; loop_nao = nao; loop_N = N;
        loop_Tele = Tele; loop_nelec = n_elec; loop_nocc_pairs = nocc_pairs; loop_step = 0;
        b_maxh = max_hist; b_iter = 0; b_push = 0; b_alpha = alpha; b_w0 = w0;
        b_step = 0; pc_active = false; pc_fired_at = -1; pc_guard = nullptr;   // WP5: M is per SCF
        loop_gfn1 = false;
        loop_open = false; loop_nalpha = loop_nbeta = 0.0;
        psResetLoop();   // WP7: a new SCF starts without usable previous eigenvectors
        loop_ready = true;
        return true;
    }
    /**
     * @brief Start a resident GFN1 SCF (WP3): SCC vector = the shell charges alone (host packSCC for GFN1),
     *        no D4 weights and no multipole state. Needs the device integrals (gamma, H0, S) of this geometry.
     *        Port of the CUDA beginResidentLoopGfn1. Claude Generated (Oct 2026).
     * @param gamma3_at Gamma_A per atom (coulomb::atom_hubbard_deriv_gfn1); q_at0 the q_at of the initial
     *        guess (third order of the first step, as the host loop)
     */
    bool resBeginLoopGfn1(int nsh, int nat, int nao, double Tele, double n_elec, int nocc_pairs,
                          const double* q_sh0, const double* q_at0, const double* n0_sh, const double* n0_at,
                          const double* gamma3_at, double alpha, int max_hist, double w0) {
        loop_ready = false;
        conv_pot_valid = false;   // WP8b
        solv_active = false; solv_failed = false;   // WP6: the GFN1 Born matrix follows this begin (beginSolvation)
        res_error.clear();
        if (resFailMode() == "begin") { res_error = "forced by CURCUMA_ROCM_FAIL_RES=begin"; return false; }
        if (!basis_ready || bis_gfn2 || nao != res_n || nao != bnao || nsh != bnsh || nat != bnat || !dgamma
            || !dao2sh || !dao2at || !dF || !dP) {
            res_error = "device integrals of this geometry are not resident"; return false;
        }
        if (!at_csr_ok) { res_error = "shells or AOs are not atom-contiguous"; return false; }
        if (max_hist <= 0 || max_hist > 20) { res_error = "Broyden history outside 1..20"; return false; }
        if (!resAllocLoop(nsh, nsh, nat, nao, max_hist)) return false;
        if (hipMemcpy(dResX, q_sh0, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResQat0, q_at0, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResN0sh, n0_sh, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResN0at, n0_at, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess
            || hipMemcpy(dResG3at, gamma3_at, nat * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) {
            res_error = "upload of the initial guess failed"; return false;
        }
        loop_nsh = nsh; loop_nat = nat; loop_nao = nao; loop_N = nsh;
        loop_Tele = Tele; loop_nelec = n_elec; loop_nocc_pairs = nocc_pairs; loop_step = 0;
        b_maxh = max_hist; b_iter = 0; b_push = 0; b_alpha = alpha; b_w0 = w0;
        b_step = 0; pc_active = false; pc_fired_at = -1; pc_guard = nullptr;   // WP5: M is per SCF
        loop_gfn1 = true;
        loop_open = false; loop_nalpha = loop_nbeta = 0.0;
        psResetLoop();   // WP7: a new SCF starts without usable previous eigenvectors
        loop_ready = true;
        return true;
    }
    /// Loop buffers for an SCC vector of length N (kept while the sizes do not change). Claude Generated (Oct 2026).
    bool resAllocLoop(int N, int nsh, int nat, int nao, int max_hist) {
        const size_t want[5] = { static_cast<size_t>(N), static_cast<size_t>(nsh), static_cast<size_t>(nat),
                                 static_cast<size_t>(nao), static_cast<size_t>(max_hist) };
        if (std::equal(want, want + 5, loop_cap)) return true;
        freeResLoopBuffers();
        const size_t uN = N, un = nat, uao = nao, ush = nsh, mh = max_hist;
        bool g = grow(dResX, uN) && grow(dResXout, uN) && grow(dResXnext, uN) && grow(dResQat0, un)
              && grow(dResN0sh, ush) && grow(dResN0at, un) && grow(dResNsh, ush) && grow(dResQsh, ush)
              && grow(dResOutMp, 10 * un) && grow(dResPart, 9 * uao) && grow(dResEat, un)
              && grow(dResTmp, std::max(ush, uao)) && grow(dResScal, 8) && grow(dResNcol, 1)
              && grow(dResVsh, ush) && grow(dResVat, un) && grow(dResG3at, un) && grow(dResQin, un)
              && grow(dResOccA, uao) && grow(dResOccB, uao) && grow(dResDf, uao)
              && grow(dBF, uN) && grow(dBFlast, uN) && grow(dBVinLast, uN) && grow(dBDFtmp, uN)
              && grow(dBDF, uN * mh) && grow(dBU, uN * mh) && grow(dBGram, mh * mh) && grow(dBC, mh)
              && grow(dBGamma, mh);
        if (!g) { freeResLoopBuffers(); res_error = "device allocation for the resident loop failed"; return false; }
        std::copy(want, want + 5, loop_cap);
        return true;
    }
    void freeResLoopBuffers() {
        for (void* p : { (void*)dResX, (void*)dResXout, (void*)dResXnext, (void*)dResQat0, (void*)dResN0sh,
                         (void*)dResN0at, (void*)dResNsh, (void*)dResQsh, (void*)dResOutMp, (void*)dResPart,
                         (void*)dResEat, (void*)dResTmp, (void*)dResScal, (void*)dResNcol,
                         (void*)dResVsh, (void*)dResVat, (void*)dResG3at, (void*)dResQin,
                         (void*)dResOccA, (void*)dResOccB, (void*)dResDf,
                         (void*)dBF, (void*)dBFlast, (void*)dBVinLast, (void*)dBDFtmp, (void*)dBDF, (void*)dBU,
                         (void*)dBGram, (void*)dBC, (void*)dBGamma })
            if (p) (void)hipFree(p);
        dResX = dResXout = dResXnext = dResQat0 = dResN0sh = dResN0at = dResNsh = dResQsh = nullptr;
        dResOutMp = dResPart = dResEat = dResTmp = dResScal = nullptr; dResNcol = nullptr;
        dResVsh = dResVat = dResG3at = dResQin = nullptr;
        dResOccA = dResOccB = dResDf = nullptr;
        dBF = dBFlast = dBVinLast = dBDFtmp = dBDF = dBU = dBGram = dBC = dBGamma = nullptr;
        std::fill(loop_cap, loop_cap + 5, size_t(0)); loop_ready = false;
    }
    /**
     * @brief Device Broyden update x_next = BroydenMixer::update(vin, vout) (broyden_mixer.h): linear step
     *        vin + alpha F in the first update and when |F - F_last| < 1e-14, otherwise the modified Broyden
     *        step with the ring history (weights 1, regularisation w0^2). Port of CUDA runBroydenUpdate.
     *        WP5 (Claude Generated, Oct 2026): with M active (pcSet) alpha*F becomes M*F and alpha*dF becomes
     *        M*dF, as with BroydenMixer::setStepHook; the guard gets ||F||_2 (one scalar to the host) and on a
     *        hit drops M (pcDrop), with "best" from the stored x_best / F_best. Without M only the norm of dF
     *        crosses the bus and the kernels are those of WP2.
     */
    const char* resBroyden(const double* vin, const double* vout, double* vnext) {
        const int N = loop_N, maxh = b_maxh;
        const size_t bytes = sizeof(double) * static_cast<size_t>(N);
        const int b1 = 256, grid = (N + b1 - 1) / b1;
        const double alpha = b_alpha;
        hipLaunchKernelGGL(k_res_vec_sub, dim3(grid), dim3(b1), 0, stream, dBF, vout, vin, N);
        ++b_iter;
        const int step = b_step++;
        const double* vin_eff = vin;   // x_best after a guard hit with "best"
        if (pc_active && pc_guard) {
            double fnorm = 0.0;
            if (rocblas_dnrm2(handle, N, dBF, 1, &fnorm) != rocblas_status_success) return "Broyden nrm2 failed";
            bool store_best = false;
            const bool hit = pc_guard(step, fnorm, &store_best);
            if (store_best
                && (hipMemcpyAsync(dPcXbest, vin, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess
                    || hipMemcpyAsync(dPcFbest, dBF, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess))
                return "Broyden copy failed";
            if (hit) {
                pc_fired_at = step;
                if (!pcDrop(pc_keep)) return "preconditioner drop failed";
                if (pc_best) {
                    vin_eff = dPcXbest;
                    if (hipMemcpyAsync(dBF, dPcFbest, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess)
                        return "Broyden copy failed";
                }
            }
        }
        const bool pc = pc_active;
        const double one = 1.0, zero = 0.0, neg = -1.0;
        // vnext += alpha F (scalar seed) or M F (preconditioned seed)
        auto add_seed = [&]() -> const char* {
            if (pc) {
                if (!pcApply(dBF, dPcTmp)) return "preconditioner kernel failed";
                if (rocblas_daxpy(handle, N, &one, dPcTmp, 1, vnext, 1) != rocblas_status_success) return "Broyden axpy failed";
            } else if (rocblas_daxpy(handle, N, &alpha, dBF, 1, vnext, 1) != rocblas_status_success) return "Broyden axpy failed";
            return nullptr;
        };
        auto linear = [&]() -> const char* {
            if (hipMemcpyAsync(vnext, vin_eff, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess) return "Broyden copy failed";
            if (const char* e = add_seed()) return e;
            if (hipMemcpyAsync(dBVinLast, vin_eff, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess
                || hipMemcpyAsync(dBFlast, dBF, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess) return "Broyden copy failed";
            return nullptr;
        };
        if (b_iter == 1) return linear();
        hipLaunchKernelGGL(k_res_vec_sub, dim3(grid), dim3(b1), 0, stream, dBDFtmp, dBF, dBFlast, N);
        double norm = 0.0;
        if (rocblas_dnrm2(handle, N, dBDFtmp, 1, &norm) != rocblas_status_success) return "Broyden nrm2 failed";
        if (norm < 1.0e-14) return linear();
        const double inv_norm = 1.0 / norm;
        const int slot = b_push % maxh;
        double* dfcol = dBDF + static_cast<size_t>(slot) * N;
        double* ucol = dBU + static_cast<size_t>(slot) * N;
        hipLaunchKernelGGL(k_res_broyden_dfu, dim3(grid), dim3(b1), 0, stream,
                           dfcol, ucol, dBDFtmp, vin_eff, dBVinLast, alpha, inv_norm, N);
        if (pc) {   // u = M dF + dx/|dF| overwrites the scalar u of k_res_broyden_dfu
            if (!pcApply(dfcol, dPcTmp)) return "preconditioner kernel failed";
            hipLaunchKernelGGL(k_res_broyden_u_pc, dim3(grid), dim3(b1), 0, stream, ucol, dPcTmp, vin_eff, dBVinLast,
                               inv_norm, N);
        }
        ++b_push;
        const int M = (b_push < maxh) ? b_push : maxh;
        // c = DF^T F (M), a = DF^T DF (M x M, column-major)
        if (rocblas_dgemv(handle, rocblas_operation_transpose, N, M, &one, dBDF, N, dBF, 1, &zero, dBC, 1)
            != rocblas_status_success) return "Broyden gemv failed";
        if (rocblas_dgemm(handle, rocblas_operation_transpose, rocblas_operation_none, M, M, N, &one,
                          dBDF, N, dBDF, N, &zero, dBGram, M) != rocblas_status_success) return "Broyden gemm failed";
        hipLaunchKernelGGL(k_res_broyden_solve, dim3(1), dim3(1), 0, stream, M, dBGram, dBC, b_w0 * b_w0, dBGamma);
        // vnext = vin + alpha F - U gamma (vin + M F - U gamma with M)
        if (hipMemcpyAsync(vnext, vin_eff, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess) return "Broyden copy failed";
        if (const char* e = add_seed()) return e;
        if (rocblas_dgemv(handle, rocblas_operation_none, N, M, &neg, dBU, N, dBGamma, 1, &one, vnext, 1)
            != rocblas_status_success) return "Broyden gemv failed";
        if (hipMemcpyAsync(dBVinLast, vin_eff, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess
            || hipMemcpyAsync(dBFlast, dBF, bytes, hipMemcpyDeviceToDevice, stream) != hipSuccess) return "Broyden copy failed";
        if (hipGetLastError() != hipSuccess) return "Broyden kernel launch failed";
        return nullptr;
    }
    /// One fused SCF iteration on the device (see the block comment above). Outputs: max|dq|, band energy
    /// Tr(P H0), Coulomb, third-order and multipole energies of the output charges.
    bool resStep(bool fp32, double* dq, double* e_band, double* e_coul, double* e_third, double* e_mp) {
        res_error.clear();
        conv_pot_valid = false;   // WP8b: the potential buffers are overwritten by this step
        auto fail = [&](const char* why) { if (res_error.empty()) res_error = why; return false; };
        if (!loop_ready) return fail("resident loop not set up");
        if (resFailMode() == "step" && ++res_inject_calls >= 3) return fail("forced by CURCUMA_ROCM_FAIL_RES=step");
        const int nsh = loop_nsh, nat = loop_nat, n = loop_nao, N = loop_N;
        // GFN1 (WP3, Claude Generated Oct 2026): SCC vector = q_sh, no D4 weights, no moments, no multipole energy.
        const bool g1 = loop_gfn1;
        if (g1) {
            if (n != res_n || n != bnao || nsh != bnsh || nat != bnat || bis_gfn2 || !basis_ready || !at_csr_ok)
                return fail("device state changed since the loop began");
        } else if (n != res_n || nsh != pot_nsh || nat != bnat || nat != d4_nat || nat != rw_nat || mp_nat != nat
                   || mp_ver < (1ull << 62) || !at_csr_ok) return fail("device state changed since the loop began");
        const size_t un = static_cast<size_t>(nat);
        const int b = 128, gnat = (nat + b - 1) / b, gsh = (nsh + b - 1) / b, gao = (n + b - 1) / b;
        if (g1) {
            // 1-3 (GFN1). q_at of the mixed input (host unpackSCC; step 0: the q_at of the initial guess, which
            // the host third order reads before the first unpack), then potential, Fock, eigensolve.
            hipLaunchKernelGGL(k_res_unpack_atoms, dim3(gnat), dim3(b), 0, stream, nsh, nat, dResX, dAtShOff,
                               dResQin, 0);
            if (hipGetLastError() != hipSuccess) return fail("input kernel failed");
            ps_in_step = true;   // WP7: the eigensolve of this step may be a pseudo step
            const char* err = potStepCoreGfn1(loop_step == 0 ? dResQat0 : dResQin, fp32);
            ps_in_step = false;
            if (err) return fail(err);
        } else {
        // 1. unpack the mixed input: q_sh -> dPotQsh, q_at (sum over shells) | dp | qp -> dMpIn
        if (hipMemcpyAsync(dPotQsh, dResX, nsh * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess)
            return fail("input copy failed");
        hipLaunchKernelGGL(k_res_unpack_atoms, dim3(gnat), dim3(b), 0, stream, nsh, nat, dResX, dAtShOff, dMpIn, 1);
        // 2. D4 reference weights at the input q_at (step 0: the q_at of the initial guess, as the host loop)
        hipLaunchKernelGGL(k_res_d4_refw, dim3(gnat), dim3(b), 0, stream, nat, D4_MAX_REF,
                           loop_step == 0 ? dResQat0 : dMpIn, dRwCn, dRwGi, dRwZeff, dRwNref, dRwRefcn, dRwRefcov,
                           dRwRefq, dD4W, dD4dWq, static_cast<double*>(nullptr));
        if (hipGetLastError() != hipSuccess) return fail("input / D4 weight kernel failed");
        // 3. potential, Fock, eigensolve
        ps_in_step = true;   // WP7: the eigensolve of this step may be a pseudo step
        const char* err = potStepCore(fp32);
        ps_in_step = false;
        if (err) return fail(err);
        mp_pot_valid = false;
        }
        // 4. occupations (resident eps -> dOcc), ncol to the host for the GEMM size
        const int blk = 256;
        if (loop_open)   // WP4: n_alpha / n_beta fillings over the same eigenvectors, summed
            hipLaunchKernelGGL(k_res_occupations_uhf, dim3(1), dim3(blk), blk * sizeof(double), stream, dEps, dOcc,
                               dResOccA, dResOccB, dResDf, n, loop_Tele * 3.166808e-6, loop_nalpha, loop_nbeta,
                               dResNcol);
        else
            hipLaunchKernelGGL(k_res_occupations, dim3(1), dim3(blk), blk * sizeof(double), stream, dEps, dOcc, n,
                               loop_Tele * 3.166808e-6, loop_nelec, loop_nocc_pairs, loop_Tele > 0.0 ? 1 : 0,
                               dResNcol);
        int ncol = n;
        if (hipGetLastError() != hipSuccess
            || hipMemcpy(&ncol, dResNcol, sizeof(int), hipMemcpyDeviceToHost) != hipSuccess)
            return fail("occupation kernel failed");
        if (ncol < 0 || ncol > n) ncol = n;
        // 5. density P = C diag(occ) C^T over the occupied columns (as resDensity), populations, band
        const bool syrk = opt_density_syrk && ncol > 0;
        hipLaunchKernelGGL(k_res_occ_scale, dim3(gao), dim3(b), 0, stream, n, dOcc, syrk ? 1 : 0, dResTmp);
        dim3 blk2(16, 16), grd2((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_scale_cols, grd2, blk2, 0, stream, dF, dResTmp, dCw, n);
        const double one = 1.0, zero = 0.0;
        if (syrk) {
            if (rocblas_dsyrk(handle, rocblas_fill_lower, rocblas_operation_none, n, ncol,
                              &one, dCw, n, &zero, dP, n) != rocblas_status_success) return fail("density syrk failed");
            hipLaunchKernelGGL(k_mirror_lower, grd2, blk2, 0, stream, dP, n);
        } else if (ncol > 0) {
            if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                              n, n, ncol, &one, dCw, n, dF, n, &zero, dP, n) != rocblas_status_success)
                return fail("density gemm failed");
        } else if (hipMemsetAsync(dP, 0, static_cast<size_t>(n) * n * sizeof(double), stream) != hipSuccess)
            return fail("density reset failed");
        hipLaunchKernelGGL(k_popband, dim3((n + 63) / 64), dim3(64), 0, stream, dP, dS0, dH0, dPB, n);
        hipLaunchKernelGGL(k_res_seqsum, dim3(1), dim3(1), 0, stream, n, dPB + n, nullptr, 0, dResScal + 1);
        // 6. output charges and moments (fixed-order gathers)
        hipLaunchKernelGGL(k_res_qsh, dim3(gsh), dim3(b), 0, stream, nsh, dPB, diao_sh, dnao_sh, dResN0sh,
                           dResNsh, dResQsh);
        hipLaunchKernelGGL(k_res_qat, dim3(gnat), dim3(b), 0, stream, nat, dAtShOff, dResNsh, dResN0at, dResOutMp);
        if (!g1) {
            hipLaunchKernelGGL(k_res_mom_part, dim3(gao), dim3(b), 0, stream, n, dP, dDpInt, dQpInt, dResPart);
            hipLaunchKernelGGL(k_res_mom_gather, dim3(gnat), dim3(b), 0, stream, nat, dAtAoOff, dResPart,
                               dResOutMp + un, dResOutMp + 4 * un);
        }
        if (hipGetLastError() != hipSuccess) return fail("density / charge kernel failed");
        // 7. SCC energies of the output: q.gamma.q (gemv + ordered dot), third order, multipole
        //    (GFN1: atom-resolved third order sum q_A^3 Gamma_A / 3, no multipole)
        if (rocblas_dgemv(handle, rocblas_operation_none, nsh, nsh, &one, dgamma, nsh, dResQsh, 1,
                          &zero, dResTmp, 1) != rocblas_status_success) return fail("Coulomb energy gemv failed");
        hipLaunchKernelGGL(k_res_seqsum, dim3(1), dim3(1), 0, stream, nsh, dResQsh, dResTmp, 1, dResScal + 2);
        if (g1) {
            hipLaunchKernelGGL(k_res_seqsum, dim3(1), dim3(1), 0, stream, nat, dResOutMp, dResG3at, 2, dResScal + 3);
            if (hipMemsetAsync(dResScal + 4, 0, sizeof(double), stream) != hipSuccess) return fail("energy reset failed");
        } else {
        hipLaunchKernelGGL(k_res_seqsum, dim3(1), dim3(1), 0, stream, nsh, dResQsh, dGamma3, 2, dResScal + 3);
        launchMpEnergy(nat, dResOutMp, dResEat);   // WP11: stored matrices or on the fly
        hipLaunchKernelGGL(k_res_seqsum, dim3(1), dim3(1), 0, stream, nat, dResEat, nullptr, 0, dResScal + 4);
        }
        // 8. pack the output, max|x_out - x_in|, Broyden, next input
        hipLaunchKernelGGL(k_res_pack_out, dim3((N + 255) / 256), dim3(256), 0, stream, nsh, nat, dResQsh,
                           dResOutMp, dResXout);
        hipLaunchKernelGGL(k_res_maxabsdiff, dim3(1), dim3(blk), blk * sizeof(double), stream, dResXout, dResX, N,
                           dResScal + 0);
        if (hipGetLastError() != hipSuccess) return fail("energy / pack kernel failed");
        if (const char* err = resBroyden(dResX, dResXout, dResXnext)) return fail(err);
        double sc[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        if (hipMemcpy(sc, dResScal, 5 * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess) return fail("scalar download failed");
        std::swap(dResX, dResXnext);   // the mixed input of the next step
        ++loop_step;
        ps_last_dq = sc[0];            // WP7: gates the next step's pseudo-diagonalisation
        *dq = sc[0]; *e_band = sc[1]; *e_coul = 0.5 * sc[2]; *e_third = sc[3]; *e_mp = sc[4];
        return true;
    }
    /**
     * @brief Open shell on the resident loop (WP4): the following steps fill the eigenvectors with n_alpha and
     *        n_beta electrons per spin channel (k_res_occupations_uhf) instead of the closed-shell filling of
     *        n_elec. Call after resBeginLoop / resBeginLoopGfn1 (a new begin resets it). Claude Generated (Oct 2026).
     */
    bool resSetOpenShell(double n_alpha, double n_beta) {
        if (!loop_ready || loop_nao <= 0 || !(n_alpha >= 0.0) || !(n_beta >= 0.0)) {
            res_error = "open shell needs a resident loop that has begun"; return false;
        }
        loop_open = true; loop_nalpha = n_alpha; loop_nbeta = n_beta; loop_nelec = n_alpha + n_beta;
        return true;
    }
    /// Converged output charges/moments and eigenvalues of the last step to the host.
    bool resCharges(double* q_sh, double* q_at, double* dp, double* qp, double* eps) {
        if (!loop_ready) { res_error = "resident loop not set up"; return false; }
        const size_t un = static_cast<size_t>(loop_nat);
        bool g = true;
        if (q_sh) g = g && hipMemcpy(q_sh, dResQsh, loop_nsh * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
        if (q_at) g = g && hipMemcpy(q_at, dResOutMp, un * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
        if (dp && !loop_gfn1)   // GFN1 has no moments: the host arrays stay as passed (zero)
            g = g && hipMemcpy(dp, dResOutMp + un, 3 * un * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
        if (qp && !loop_gfn1)
            g = g && hipMemcpy(qp, dResOutMp + 4 * un, 6 * un * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
        if (eps)  g = g && hipMemcpy(eps, dEps, loop_nao * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
        if (!g) res_error = "download of the converged charges failed";
        return g;
    }

    // ----- Stage 5 (Part B2): in-SCF GFN2 D4 dE_D4/dq ---------------------
    // begDispersion uploads the geometry-fixed reference data once per geometry; the
    // element C6 block (c6_flat) is uploaded only on the first call (process lifetime).
    bool begDispersion(int nat, const int* Z, const double* sqrtZr4r2, const int* nref,
                       const double* xyz, const double* c6_flat, int c6_len,
                       double s6, double s8, double a1, double a2, double cutoff) {
        if (!ok || nat <= 0 || c6_len <= 0) return false;
        rw_geom_ok = false;    // WP8a: the D4 reference tables of the previous geometry are stale now
        if (nat != d4_nat) {   // atom count changed (new molecule) → drop the nat-sized buffers
            for (void* p : { (void*)dD4W, (void*)dD4dWq, (void*)dD4Dedq,
                             (void*)dD4dWc, (void*)dD4Eat, (void*)dD4Grad, (void*)dD4Dcn,
                             (void*)dD4c6, (void*)dD4dc6 }) if (p) hipFree(p);
            dD4W = dD4dWq = dD4Dedq = nullptr;
            dD4dWc = dD4Eat = dD4Grad = dD4Dcn = nullptr;
            dD4c6 = dD4dc6 = nullptr;
        }
        bool g = true;
        g &= up(dD4Z, Z, static_cast<size_t>(nat));
        g &= up(dD4Nref, nref, static_cast<size_t>(nat));
        g &= up(dD4Sqrt, sqrtZr4r2, static_cast<size_t>(nat));
        g &= up(dD4Xyz, xyz, 3 * static_cast<size_t>(nat));
        if (!dD4W)    g &= hipMalloc(&dD4W,    static_cast<size_t>(nat) * D4_MAX_REF * sizeof(double)) == hipSuccess;
        if (!dD4dWq)  g &= hipMalloc(&dD4dWq,  static_cast<size_t>(nat) * D4_MAX_REF * sizeof(double)) == hipSuccess;
        if (!dD4Dedq) g &= hipMalloc(&dD4Dedq, static_cast<size_t>(nat) * sizeof(double)) == hipSuccess;
        // The reference C6 block is element data (geometry- AND molecule-independent):
        // upload it only once per process.
        if (!d4_c6_uploaded || d4_c6_len < c6_len) {
            if (dD4C6) hipFree(dD4C6);
            dD4C6 = nullptr;
            g &= hipMalloc(&dD4C6, static_cast<size_t>(c6_len) * sizeof(double)) == hipSuccess;
            if (g) g &= hipMemcpy(dD4C6, c6_flat, static_cast<size_t>(c6_len) * sizeof(double),
                                  hipMemcpyHostToDevice) == hipSuccess;
            d4_c6_uploaded = g;
            d4_c6_len = c6_len;
        }
        d4_nat = nat;
        d4_s6 = s6; d4_s8 = s8; d4_a1 = a1; d4_a2 = a2; d4_cut = cutoff;
        return g && hipStreamSynchronize(stream) == hipSuccess;
    }
    // Per SCF iteration: upload W/dWq (the q-dependent reference weights), run the O(N²)
    // contraction + BJ disp_sum kernel, download dEdq. Mirrors XTB::addDispersionPotential.
    bool dispDedq(int nat, const double* W, const double* dWq, double* dEdq_out) {
        if (!ok || nat <= 0 || nat != d4_nat || !dD4C6) return false;
        const size_t wlen = static_cast<size_t>(nat) * D4_MAX_REF;
        if (hipMemcpy(dD4W,   W,   wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dD4dWq, dWq, wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        const double cut2 = d4_cut * d4_cut;
        dim3 b(128), gr((nat + 127) / 128);
        hipLaunchKernelGGL(k_d4_dedq, gr, b, 0, stream, nat, 118, D4_MAX_REF,
                           dD4Z, dD4Nref, dD4Sqrt, dD4Xyz, dD4C6, dD4W, dD4dWq,
                           d4_s6, d4_s8, d4_a1, d4_a2, cut2, dD4Dedq);
        if (hipGetLastError() != hipSuccess) return false;
        if (hipMemcpy(dEdq_out, dD4Dedq, static_cast<size_t>(nat) * sizeof(double),
                      hipMemcpyDeviceToHost) != hipSuccess) return false;
        return hipStreamSynchronize(stream) == hipSuccess;
    }
    // Post-SCF: the whole 2-body D4 (energy + gradient + dE/dCN + dE/dq) in one gather.
    bool dispDispGrad(int nat, const double* W, const double* dWq, const double* dWc,
                      double* e_atom, double* grad, double* dEdcn, double* dEdq_out) {
        if (!ok || nat <= 0 || nat != d4_nat || !dD4C6) return false;
        const size_t wlen = static_cast<size_t>(nat) * D4_MAX_REF;
        if (!dD4dWc  && hipMalloc(&dD4dWc,  wlen * sizeof(double)) != hipSuccess) return false;
        if (!dD4Eat  && hipMalloc(&dD4Eat,  static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) return false;
        if (!dD4Grad && hipMalloc(&dD4Grad, 3 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) return false;
        if (!dD4Dcn  && hipMalloc(&dD4Dcn,  static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) return false;
        if (hipMemcpy(dD4W,   W,   wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dD4dWq, dWq, wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dD4dWc, dWc, wlen * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        const double cut2 = d4_cut * d4_cut;
        dim3 b(128), gr((nat + 127) / 128);
        hipLaunchKernelGGL(k_d4_grad, gr, b, 0, stream, nat, 118, D4_MAX_REF,
                           dD4Z, dD4Nref, dD4Sqrt, dD4Xyz, dD4C6, dD4W, dD4dWq, dD4dWc,
                           d4_s6, d4_s8, d4_a1, d4_a2, cut2,
                           dD4Eat, dD4Grad, dD4Dcn, dD4Dedq);
        if (hipGetLastError() != hipSuccess) return false;
        if (hipMemcpy(e_atom, dD4Eat, static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(grad, dD4Grad, 3 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(dEdcn, dD4Dcn, static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(dEdq_out, dD4Dedq, static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return hipStreamSynchronize(stream) == hipSuccess;
    }
    // Post-SCF: the D4 ATM 3-body (energy + gradient + dE/dCN) in one gather, reusing the
    // resident geometry + √r4r2. c6/dc6dcn are nat² (q=0 reference). Outputs reuse dD4Eat/Grad/Dcn.
    bool dispATM(int nat, const double* c6, const double* dc6dcn,
                 double s9, double a1, double a2, double alp, double cutoff,
                 double* e_atom, double* grad, double* dEdcn) {
        if (!ok || nat <= 0 || nat != d4_nat || !dD4Xyz || !dD4Sqrt) return false;
        const size_t nn = static_cast<size_t>(nat) * nat;
        // nat² matrices: (re)allocate when the size changes (begDispersion resets d4_nat).
        if (!dD4c6  && hipMalloc(&dD4c6,  nn * sizeof(double)) != hipSuccess) return false;
        if (!dD4dc6 && hipMalloc(&dD4dc6, nn * sizeof(double)) != hipSuccess) return false;
        if (!dD4Eat  && hipMalloc(&dD4Eat,  static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) return false;
        if (!dD4Grad && hipMalloc(&dD4Grad, 3 * static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) return false;
        if (!dD4Dcn  && hipMalloc(&dD4Dcn,  static_cast<size_t>(nat) * sizeof(double)) != hipSuccess) return false;
        if (hipMemcpy(dD4c6,  c6,      nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dD4dc6, dc6dcn,  nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        const double cut2 = cutoff * cutoff;
        dim3 b(128), gr((nat + 127) / 128);
        hipLaunchKernelGGL(k_d4_atm, gr, b, 0, stream, nat, dD4Xyz, dD4Sqrt, dD4c6, dD4dc6,
                           s9, a1, a2, alp, cut2, dD4Eat, dD4Grad, dD4Dcn);
        if (hipGetLastError() != hipSuccess) return false;
        if (hipMemcpy(e_atom, dD4Eat, static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(grad, dD4Grad, 3 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(dEdcn, dD4Dcn, static_cast<size_t>(nat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return hipStreamSynchronize(stream) == hipSuccess;
    }
    // ----- WP8a (Claude Generated, Oct 2026): post-SCF D4, GFN1 D3(BJ), D4 covalent CN, Hamiltonian CN -----
    // Port of the CUDA d4PostScfDevice / d3DispersionDevice / d4CovalentCnDevice / downloadCn. The post-SCF D4
    // reuses the begDispersion data (Z, nref, sqrt r4r2, xyz, element C6 block) and the resBeginWeights tables
    // of the same geometry (rw_geom_ok); the weights come from k_res_d4_refw (host buildAtomRefW incl. its
    // underflow rule) and the 2-body term from k_d4_grad. Fault injection CURCUMA_ROCM_FAIL_POSTSCF=d4|d3|cn.
    template <class T> struct DBuf {
        T* p = nullptr; size_t cap = 0;
        bool ensure(size_t n) {
            if (p && cap >= n) return true;
            if (p) (void)hipFree(p);
            p = nullptr; cap = 0;
            if (hipMalloc(&p, std::max<size_t>(n, 1) * sizeof(T)) != hipSuccess) { (void)hipGetLastError(); p = nullptr; return false; }
            cap = n;
            return true;
        }
        void release() { if (p) (void)hipFree(p); p = nullptr; cap = 0; }
    };
    bool rw_geom_ok = false;          ///< resBeginWeights ran after the last begDispersion (tables of this geometry)
    bool cn_ok = false;               ///< dcn holds the Hamiltonian CN of the last successful computeIntegrals
    std::string pst_error;             ///< why the last post-SCF call failed or declined
    bool pst_unavailable = false;      ///< the last call declined (precondition), not a device failure
    double pst_atm_ms = 0.0;           ///< ATM kernel time of the last d4PostScf (hipEvent), -1 if not run
    DBuf<double> psQ, psQ0, psRad, psEn, psW0, psdWq0, psdWc0, psAtmEat, psAtmGrad, psAtmDcn, psGtot, psDtot;
    DBuf<int> psNbCnt, psNbPtr, psNb;
    DBuf<int> cnZ; DBuf<double> cnXyz, cnOut;
    DBuf<int> d3Z, d3Type, d3Nref; DBuf<double> d3Xyz, d3Cn, d3Gw, d3Dgw, d3Eat, d3Grad, d3Dcn, d3C6, d3R6, d3Refcn, d3Rcov;

    static const std::string& pstFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_POSTSCF"); return std::string(e ? e : ""); }();
        return m;
    }
    bool pstDecline(const char* why) { pst_error = why; pst_unavailable = true; return false; }
    bool pstFail(const char* why) { pst_error = why; pst_unavailable = false; (void)hipGetLastError(); return false; }

    /// Whole post-SCF D4 of GFN2: weights at q_at, 2-body (energy, gradient, dE/dCN, dE/dq), ATM at q = 0 (C6
    /// matrices, neighbour list, three-body kernel), dE/dCN fold into the gradient with the dftd4 covalent CN.
    bool d4PostScf(int nat, const double* q_at, bool with_atm, double s9, double a1, double a2, double alp,
                   double atm_cutoff, double cn_cutoff, const double* rad, const double* en, double kn, double k4,
                   double k5, double k6, double* energy_out, double* grad_out, double* dEdq_out) {
        pst_error.clear(); pst_unavailable = false; pst_atm_ms = -1.0;
        if (pstFailMode() == "d4") return pstFail("forced by CURCUMA_ROCM_FAIL_POSTSCF=d4");
        if (!ok || nat <= 0 || nat != d4_nat || !dD4C6 || !dD4Xyz || !dD4W || !dD4dWq || !dD4Dedq)
            return pstDecline("no device D4 data of this geometry");
        if (!rw_geom_ok || nat != rw_nat)
            return pstDecline("no device D4 reference tables of this geometry (host-driven SCF loop)");
        const size_t un = static_cast<size_t>(nat), wl = un * D4_MAX_REF, nn = un * un;
        bool g = psQ.ensure(un) && psRad.ensure(119) && psEn.ensure(119) && psGtot.ensure(3 * un)
              && psDtot.ensure(un);
        if (g && !dD4dWc)  g = hipMalloc(&dD4dWc,  wl * sizeof(double)) == hipSuccess;
        if (g && !dD4Eat)  g = hipMalloc(&dD4Eat,  un * sizeof(double)) == hipSuccess;
        if (g && !dD4Grad) g = hipMalloc(&dD4Grad, 3 * un * sizeof(double)) == hipSuccess;
        if (g && !dD4Dcn)  g = hipMalloc(&dD4Dcn,  un * sizeof(double)) == hipSuccess;
        if (g && with_atm) {
            g = psQ0.ensure(un) && psW0.ensure(wl) && psdWq0.ensure(wl) && psdWc0.ensure(wl) && psAtmEat.ensure(un)
                && psAtmGrad.ensure(3 * un) && psAtmDcn.ensure(un) && psNbCnt.ensure(un) && psNbPtr.ensure(un + 1);
            if (g && !dD4c6)  g = hipMalloc(&dD4c6,  nn * sizeof(double)) == hipSuccess;
            if (g && !dD4dc6) g = hipMalloc(&dD4dc6, nn * sizeof(double)) == hipSuccess;
        }
        if (!g) return pstFail("device allocation for the post-SCF D4 failed");
        if (hipMemcpyAsync(psQ.p, q_at, un * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(psRad.p, rad, 119 * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(psEn.p, en, 119 * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess)
            return pstFail("upload of the converged charges failed");
        const int b = 128, gr = (nat + b - 1) / b;
        const double cut2 = d4_cut * d4_cut;
        // 1. weights at the converged charges, 2-body energy / gradient / dE/dCN / dE/dq
        hipLaunchKernelGGL(k_res_d4_refw, dim3(gr), dim3(b), 0, stream, nat, D4_MAX_REF, psQ.p, dRwCn, dRwGi, dRwZeff,
                           dRwNref, dRwRefcn, dRwRefcov, dRwRefq, dD4W, dD4dWq, dD4dWc);
        hipLaunchKernelGGL(k_d4_grad, dim3(gr), dim3(b), 0, stream, nat, 118, D4_MAX_REF, dD4Z, dD4Nref, dD4Sqrt,
                           dD4Xyz, dD4C6, dD4W, dD4dWq, dD4dWc, d4_s6, d4_s8, d4_a1, d4_a2, cut2,
                           dD4Eat, dD4Grad, dD4Dcn, dD4Dedq);
        if (hipGetLastError() != hipSuccess) return pstFail("D4 weight / 2-body kernel failed");
        std::vector<double> e_atm;
        if (with_atm) {
            // 2. ATM: weights at q = 0, C6 matrices, neighbour list, three-body kernel
            if (hipMemsetAsync(psQ0.p, 0, un * sizeof(double), stream) != hipSuccess) return pstFail("memset failed");
            hipLaunchKernelGGL(k_res_d4_refw, dim3(gr), dim3(b), 0, stream, nat, D4_MAX_REF, psQ0.p, dRwCn, dRwGi,
                               dRwZeff, dRwNref, dRwRefcn, dRwRefcov, dRwRefq, psW0.p, psdWq0.p, psdWc0.p);
            hipLaunchKernelGGL(k_d4_atm_c6mat, dim3(static_cast<unsigned>((nn + 255) / 256)), dim3(256), 0, stream,
                               nat, dD4Z, dD4Nref, psW0.p, psdWc0.p, dD4C6, dD4c6, dD4dc6);
            const double atm_cut2 = atm_cutoff * atm_cutoff;
            hipLaunchKernelGGL(k_d4_nb_count, dim3(gr), dim3(b), 0, stream, nat, dD4Xyz, atm_cut2, psNbCnt.p);
            if (hipGetLastError() != hipSuccess) return pstFail("ATM C6 / neighbour count kernel failed");
            std::vector<int> cnt(un), ptr(un + 1, 0);
            if (hipMemcpyAsync(cnt.data(), psNbCnt.p, un * sizeof(int), hipMemcpyDeviceToHost, stream) != hipSuccess
                || hipStreamSynchronize(stream) != hipSuccess) return pstFail("download of the neighbour counts failed");
            for (size_t a = 0; a < un; ++a) ptr[a + 1] = ptr[a] + cnt[a];
            if (!psNb.ensure(static_cast<size_t>(std::max(1, ptr[un]))))
                return pstFail("device allocation for the ATM neighbour list failed");
            if (hipMemcpyAsync(psNbPtr.p, ptr.data(), (un + 1) * sizeof(int), hipMemcpyHostToDevice, stream) != hipSuccess)
                return pstFail("upload of the neighbour offsets failed");
            hipLaunchKernelGGL(k_d4_nb_fill, dim3(gr), dim3(b), 0, stream, nat, dD4Xyz, atm_cut2, psNbPtr.p, psNb.p);
            hipEvent_t e0 = nullptr, e1 = nullptr;
            const bool ev = hipEventCreate(&e0) == hipSuccess && hipEventCreate(&e1) == hipSuccess;
            if (ev) (void)hipEventRecord(e0, stream);
            // Block per atom (default) or thread per atom (CURCUMA_ROCM_ATM_THREAD=1: the summation order of the
            // host loop and of k_d4_atm, slower on gfx1200 by divergence of the per-lane neighbour lists).
            static const bool atm_thread = [] { const char* e = std::getenv("CURCUMA_ROCM_ATM_THREAD"); return e && e[0] == '1'; }();
            if (atm_thread) {
                const int bt = 64;
                hipLaunchKernelGGL(k_d4_atm_nl, dim3((nat + bt - 1) / bt), dim3(bt), 0, stream, nat, dD4Xyz, dD4Sqrt,
                                   dD4c6, dD4dc6, psNbPtr.p, psNb.p, s9, a1, a2, alp, atm_cut2,
                                   psAtmEat.p, psAtmGrad.p, psAtmDcn.p);
            } else {
                const int bb = 64;
                hipLaunchKernelGGL(k_d4_atm_nl_b, dim3(nat), dim3(bb), 5 * bb * sizeof(double), stream, nat, dD4Xyz,
                                   dD4Sqrt, dD4c6, dD4dc6, psNbPtr.p, psNb.p, s9, a1, a2, alp, atm_cut2,
                                   psAtmEat.p, psAtmGrad.p, psAtmDcn.p);
            }
            if (ev) (void)hipEventRecord(e1, stream);
            if (hipGetLastError() != hipSuccess) {
                if (e0) (void)hipEventDestroy(e0);
                if (e1) (void)hipEventDestroy(e1);
                return pstFail("ATM kernel failed");
            }
            hipLaunchKernelGGL(k_d4_add2, dim3((3 * nat + 255) / 256), dim3(256), 0, stream, 3 * nat, dD4Grad,
                               psAtmGrad.p, psGtot.p);
            hipLaunchKernelGGL(k_d4_add2, dim3((nat + 255) / 256), dim3(256), 0, stream, nat, dD4Dcn, psAtmDcn.p,
                               psDtot.p);
            e_atm.resize(un);
            if (hipMemcpyAsync(e_atm.data(), psAtmEat.p, un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess)
                return pstFail("download of the ATM energies failed");
            if (ev && hipEventSynchronize(e1) == hipSuccess) {
                float ms = 0.0f;
                if (hipEventElapsedTime(&ms, e0, e1) == hipSuccess) pst_atm_ms = ms;
            }
            if (e0) (void)hipEventDestroy(e0);
            if (e1) (void)hipEventDestroy(e1);
        } else {
            if (hipMemcpyAsync(psGtot.p, dD4Grad, 3 * un * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess
                || hipMemcpyAsync(psDtot.p, dD4Dcn, un * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess)
                return pstFail("device copy failed");
        }
        // 3. CN chain rule of dE/dCN (2-body + ATM) into the gradient
        hipLaunchKernelGGL(k_d4_cn_fold, dim3(gr), dim3(b), 0, stream, nat, dD4Z, dD4Xyz, psDtot.p, psRad.p, psEn.p,
                           kn, k4, k5, k6, cn_cutoff, psGtot.p);
        if (hipGetLastError() != hipSuccess) return pstFail("D4 CN fold kernel failed");
        std::vector<double> e2(un);
        if (hipMemcpyAsync(e2.data(), dD4Eat, un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipMemcpyAsync(grad_out, psGtot.p, 3 * un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipMemcpyAsync(dEdq_out, dD4Dedq, un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess)
            return pstFail("download of the D4 results failed");
        double E = 0.0;
        for (size_t a = 0; a < un; ++a) E += e2[a];
        E *= 0.5;                                   // the 2-body gather counts every pair twice (host: 0.5 sum)
        if (with_atm) {
            double ea = 0.0;
            for (size_t a = 0; a < un; ++a) ea += e_atm[a];
            E += ea / 3.0;                          // every triple counted by its three members
        }
        *energy_out = E;
        return true;
    }

    /// GFN1 D3(BJ) on the device: CN, weights, pair energy / gradient / dE/dCN, CN chain rule (Angstrom in,
    /// gradient in Eh/Bohr with fold_unit = Bohr per Angstrom).
    bool d3Dispersion(int nat, const int* Z, const double* xyz_ang, int nt, const int* atom_type,
                      const double* c6tab, const double* r6tab, const int* nref_tab, const double* refcn_tab,
                      const double* rcov_ang, int nrcov, double s6, double s8, double a1, double a2, double k1,
                      double k2, double cn_cut_ang, double bohr_per_ang, double fold_unit, bool need_gradient,
                      double* energy_out, double* grad_out) {
        pst_error.clear(); pst_unavailable = false;
        if (pstFailMode() == "d3") return pstFail("forced by CURCUMA_ROCM_FAIL_POSTSCF=d3");
        if (!ok || nat <= 0 || nt <= 0) return pstDecline("no atoms or no D3 types");
        const size_t un = static_cast<size_t>(nat), nt2 = static_cast<size_t>(nt) * nt;
        if (!(d3Z.ensure(un) && d3Type.ensure(un) && d3Xyz.ensure(3 * un) && d3Cn.ensure(un)
              && d3Gw.ensure(un * D3_MR) && d3Dgw.ensure(un * D3_MR) && d3Eat.ensure(un) && d3Grad.ensure(3 * un)
              && d3Dcn.ensure(un) && d3C6.ensure(nt2 * D3_MR * D3_MR) && d3R6.ensure(nt2) && d3Nref.ensure(D3_ME)
              && d3Refcn.ensure(static_cast<size_t>(D3_ME) * D3_MR) && d3Rcov.ensure(static_cast<size_t>(nrcov))))
            return pstFail("device allocation for the GFN1 D3 dispersion failed");
        auto upA = [this](void* d, const void* h, size_t bytes) {
            return hipMemcpyAsync(d, h, bytes, hipMemcpyHostToDevice, stream) == hipSuccess;
        };
        if (!(upA(d3Z.p, Z, un * sizeof(int)) && upA(d3Type.p, atom_type, un * sizeof(int))
              && upA(d3Xyz.p, xyz_ang, 3 * un * sizeof(double))
              && upA(d3C6.p, c6tab, nt2 * D3_MR * D3_MR * sizeof(double)) && upA(d3R6.p, r6tab, nt2 * sizeof(double))
              && upA(d3Nref.p, nref_tab, D3_ME * sizeof(int))
              && upA(d3Refcn.p, refcn_tab, static_cast<size_t>(D3_ME) * D3_MR * sizeof(double))
              && upA(d3Rcov.p, rcov_ang, static_cast<size_t>(nrcov) * sizeof(double))))
            return pstFail("upload of the D3 tables failed");
        const int b = 128, gr = (nat + b - 1) / b;
        hipLaunchKernelGGL(k_d3_cn, dim3(gr), dim3(b), 0, stream, nat, d3Z.p, d3Xyz.p, d3Rcov.p, nrcov, k1, k2,
                           cn_cut_ang, d3Cn.p);
        hipLaunchKernelGGL(k_d3_weights, dim3(gr), dim3(b), 0, stream, nat, d3Z.p, d3Cn.p, d3Nref.p, d3Refcn.p,
                           d3Gw.p, d3Dgw.p);
        hipLaunchKernelGGL(k_d3_pairs, dim3(gr), dim3(b), 0, stream, nat, nt, d3Z.p, d3Type.p, d3Nref.p, d3Xyz.p,
                           d3Gw.p, d3Dgw.p, d3C6.p, d3R6.p, s6, s8, a1, a2, bohr_per_ang, need_gradient ? 1 : 0,
                           d3Eat.p, d3Grad.p, d3Dcn.p);
        if (need_gradient)
            hipLaunchKernelGGL(k_d3_cn_fold, dim3(gr), dim3(b), 0, stream, nat, d3Z.p, d3Xyz.p, d3Dcn.p, d3Rcov.p,
                               nrcov, k1, k2, fold_unit, cn_cut_ang, d3Grad.p);
        if (hipGetLastError() != hipSuccess) return pstFail("D3 kernel failed");
        std::vector<double> e(un);
        if (hipMemcpyAsync(e.data(), d3Eat.p, un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || (need_gradient && hipMemcpyAsync(grad_out, d3Grad.p, 3 * un * sizeof(double), hipMemcpyDeviceToHost,
                                                stream) != hipSuccess)
            || hipStreamSynchronize(stream) != hipSuccess)
            return pstFail("download of the D3 results failed");
        double E = 0.0;
        for (size_t a = 0; a < un; ++a) E += e[a];
        *energy_out = 0.5 * E;   // every pair is gathered from both of its atoms
        return true;
    }

    /// dftd4 covalent CN of the geometry (Bohr) for the next D4 GenerateParameters.
    bool d4CovalentCn(int nat, const int* Z, const double* xyz_bohr, const double* rad, const double* en,
                      double kn, double k4, double k5, double k6, double cutoff, double* cn_out) {
        pst_error.clear(); pst_unavailable = false;
        if (pstFailMode() == "cn") return pstFail("forced by CURCUMA_ROCM_FAIL_POSTSCF=cn");
        if (!ok || nat <= 0) return pstDecline("no atoms");
        const size_t un = static_cast<size_t>(nat);
        if (!(cnZ.ensure(un) && cnXyz.ensure(3 * un) && cnOut.ensure(un) && psRad.ensure(119) && psEn.ensure(119)))
            return pstFail("device allocation for the D4 covalent CN failed");
        if (hipMemcpyAsync(cnZ.p, Z, un * sizeof(int), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(cnXyz.p, xyz_bohr, 3 * un * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(psRad.p, rad, 119 * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(psEn.p, en, 119 * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess)
            return pstFail("upload for the D4 covalent CN failed");
        const int b = 128;
        hipLaunchKernelGGL(k_d4_cn, dim3((nat + b - 1) / b), dim3(b), 0, stream, nat, cnZ.p, cnXyz.p, psRad.p, psEn.p,
                           kn, k4, k5, k6, cutoff, cnOut.p);
        if (hipGetLastError() != hipSuccess) return pstFail("D4 covalent CN kernel failed");
        if (hipMemcpyAsync(cn_out, cnOut.p, un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess)
            return pstFail("download of the D4 covalent CN failed");
        return true;
    }

    /// Hamiltonian CN of the last device integral build (k_cn in computeIntegrals).
    bool dlCn(double* cn_out) {
        pst_error.clear(); pst_unavailable = false;
        if (pstFailMode() == "cn") return pstFail("forced by CURCUMA_ROCM_FAIL_POSTSCF=cn");
        if (!ok || !cn_ok || bnat <= 0 || !dcn) return pstDecline("no device integral build of this geometry");
        if (hipMemcpy(cn_out, dcn, static_cast<size_t>(bnat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess)
            return pstFail("download of the Hamiltonian CN failed");
        return true;
    }
    // ----- WP8b (Claude Generated, Oct 2026): gradient potentials, gradient tail and repulsion on the device -----
    // Ports of the CUDA convergedPotentials (with computeGradient(..., resident_potentials)), gradientTailDevice and
    // repulsionEnergyDevice (cuda/xtb_gpu_context.cu). convPotentials builds the potentials of the gradient at the
    // converged OUTPUT charges of the resident loop (dResQsh, dResOutMp) with the potential half of the loop step
    // (potBuildCore / potBuildGfn1; GFN2: D4 reference weights at the output q_at as the host addDispersionPotential),
    // and leaves v_ao in dVao, q_sh in dQsh and v_dp/v_qp in dVdp/dVqp for gradient(..., resident = true). The host
    // rebuilt the same potentials from the downloaded charges (XTB::addCoulombShellPotential, addThirdOrderPotential,
    // addMultipolePotential, addDispersionPotential). Fault injection CURCUMA_ROCM_FAIL_GRAD=resident|tail (every call
    // of that path fails; "defer" is read by HipScfBackend for the lazy host copies).
    bool conv_pot_valid = false;      ///< dVao/dQsh/dVdp/dVqp hold the converged potentials of the current loop
    std::string grad_error;           ///< why the last WP8b call failed
    /// WP10 (Claude Generated, Oct 2026): -xtb.rocm_fixed_point, order-independent fixed-point accumulation in
    /// k_grad_repulsion, k_grad_cn_onsite, k_grad_h0_pulay, k_grad_coulomb and k_multipole_moments (xfix_add).
    bool fixed_point = false;
    /// Reads and resets xfix_bad after a synchronised fixed-point accumulation; false (with grad_error) when a
    /// contribution was outside the fixed-point range or not finite, so the caller fails over as for any device error.
    bool fixFlagClear(const char* what) {
        int bad = 0;
        if (hipMemcpyFromSymbol(&bad, HIP_SYMBOL(xfix_bad), sizeof(int)) != hipSuccess) return false;
        if (bad == 0) return true;
        const int zero = 0;
        (void)hipMemcpyToSymbol(HIP_SYMBOL(xfix_bad), &zero, sizeof(int));
        grad_error = std::string("fixed-point accumulation out of range or not finite (") + what + ")";
        return false;
    }
    DBuf<double> tlXyz, tlDcn, tlRcov, tlGrad, tlQ, tlDp, tlQp, tlMrad, tlCn, tlRad, tlVcn, tlDer, tlGmp;
    DBuf<double> rpXyz, rpAlpha, rpZeff, rpEat;
    DBuf<int> rpZ;
    static const std::string& gradFailMode() {
        static const std::string m = [] { const char* e = std::getenv("CURCUMA_ROCM_FAIL_GRAD"); return std::string(e ? e : ""); }();
        return m;
    }
    bool gradFail(const char* why) { grad_error = why; (void)hipGetLastError(); return false; }
    /**
     * @brief Potentials of the gradient at the converged output charges of the resident loop (GFN2: gfn2 = true).
     * @return false (grad_error) without a finished resident loop of the matching method, with an active device
     *         solvent (the host keeps that potential) or on a device failure
     */
    bool convPotentials(bool gfn2) {
        conv_pot_valid = false;
        grad_error.clear();
        if (gradFailMode() == "resident") return gradFail("forced by CURCUMA_ROCM_FAIL_GRAD=resident");
        if (!ok || !loop_ready || loop_step <= 0 || loop_gfn1 == gfn2 || loop_nao != res_n || loop_nao != bnao
            || loop_nsh != bnsh || loop_nat != bnat || !dVao || !dQsh || !dvao)
            return gradFail("no finished resident loop of this geometry");
        if (solv_active) return gradFail("device solvent active");
        const int nsh = loop_nsh, nat = loop_nat, n = loop_nao;
        const size_t un = static_cast<size_t>(nat);
        if (gfn2) {
            if (nsh != pot_nsh || nat != d4_nat || nat != rw_nat || mp_nat != nat || mp_ver < (1ull << 62)
                || !dVdp || !dVqp || !beginResidentMP()) return gradFail("device potential of this geometry is not set up");
            if (hipMemcpyAsync(dPotQsh, dResQsh, nsh * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess
                || hipMemcpyAsync(dMpIn, dResOutMp, 10 * un * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess)
                return gradFail("copy of the converged charges failed");
            const int b = 128;
            hipLaunchKernelGGL(k_res_d4_refw, dim3((nat + b - 1) / b), dim3(b), 0, stream, nat, D4_MAX_REF, dMpIn,
                               dRwCn, dRwGi, dRwZeff, dRwNref, dRwRefcn, dRwRefcov, dRwRefq, dD4W, dD4dWq,
                               static_cast<double*>(nullptr));
            if (hipGetLastError() != hipSuccess) return gradFail("D4 weight kernel failed");
            if (const char* err = potBuildCore()) return gradFail(err);
            mp_pot_valid = false;   // dVdp/dVqp now hold the converged potential
        } else {
            if (const char* err = potBuildGfn1(dResQsh, dResOutMp)) return gradFail(err);
        }
        if (hipGetLastError() != hipSuccess) return gradFail("potential kernel launch failed");
        if (hipMemcpyAsync(dVao, dvao, n * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess
            || hipMemcpyAsync(dQsh, dResQsh, nsh * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess) return gradFail("potential copy failed");
        conv_pot_valid = true;
        return true;
    }
    /// Sections 5 (GFN2 with_mp) and 4 of the GPU gradient (see XtbHipContext::gradientTailDevice).
    bool gradTail(int nat, bool is_gfn2, bool with_mp, const double* xyz, const double* dEdcn, const double* rcov,
                  const double* q, const double* dp, const double* qp, const double* mrad, const double* cn,
                  const double* rad, const double* vcn, double dmp3, double dmp5, double shift, double kexp,
                  double rmax, double* grad_out) {
        grad_error.clear();
        if (gradFailMode() == "tail") return gradFail("forced by CURCUMA_ROCM_FAIL_GRAD=tail");
        if (!ok || nat <= 0) return gradFail("no device");
        const size_t un = static_cast<size_t>(nat);
        bool g = tlXyz.ensure(3 * un) && tlDcn.ensure(un) && tlRcov.ensure(un) && tlGrad.ensure(3 * un);
        if (with_mp)
            g = g && tlQ.ensure(un) && tlDp.ensure(3 * un) && tlQp.ensure(6 * un) && tlMrad.ensure(un)
                  && tlCn.ensure(un) && tlRad.ensure(un) && tlVcn.ensure(un) && tlDer.ensure(un) && tlGmp.ensure(3 * un);
        if (!g) return gradFail("device allocation for the gradient tail failed");
        auto upd = [&](double* d, const double* h, size_t c) {
            return hipMemcpyAsync(d, h, c * sizeof(double), hipMemcpyHostToDevice, stream) == hipSuccess;
        };
        g = upd(tlXyz.p, xyz, 3 * un) && upd(tlDcn.p, dEdcn, un) && upd(tlRcov.p, rcov, un)
            && hipMemsetAsync(tlGrad.p, 0, 3 * un * sizeof(double), stream) == hipSuccess;
        const int b = 128, gr = (nat + b - 1) / b;
        if (g && with_mp) {
            g = upd(tlQ.p, q, un) && upd(tlDp.p, dp, 3 * un) && upd(tlQp.p, qp, 6 * un) && upd(tlMrad.p, mrad, un)
                && upd(tlCn.p, cn, un) && upd(tlRad.p, rad, un) && upd(tlVcn.p, vcn, un);
            if (g) {
                hipLaunchKernelGGL(k_grad_tail_mp, dim3(gr), dim3(b), 0, stream, nat, tlXyz.p, tlMrad.p, dmp3, dmp5,
                                   tlQ.p, tlDp.p, tlQp.p, tlGmp.p, tlDer.p);
                hipLaunchKernelGGL(k_grad_tail_dcn, dim3(gr), dim3(b), 0, stream, nat, tlDer.p, tlCn.p, tlRad.p,
                                   tlVcn.p, shift, kexp, rmax, tlDcn.p);
                g = hipMemcpyAsync(tlGrad.p, tlGmp.p, 3 * un * sizeof(double), hipMemcpyDeviceToDevice, stream)
                    == hipSuccess;
            }
        }
        if (!g) return gradFail("upload of the gradient tail inputs failed");
        hipLaunchKernelGGL(k_grad_tail_cn, dim3(gr), dim3(b), 0, stream, nat, is_gfn2 ? 0 : 1, tlXyz.p, tlRcov.p,
                           tlDcn.p, tlGrad.p);
        if (hipGetLastError() != hipSuccess) return gradFail("gradient tail kernel failed");
        if (hipMemcpyAsync(grad_out, tlGrad.p, 3 * un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess) return gradFail("download of the gradient tail failed");
        return true;
    }
    /// Repulsion energy (see XtbHipContext::repulsionEnergyDevice); per-atom sums added on the host in atom order.
    bool repEnergy(int nat, bool is_gfn2, const int* z, const double* xyz, const double* alpha, const double* zeff,
                   double kexp, double rexp, double kexp_light, double* energy_out) {
        grad_error.clear();
        if (gradFailMode() == "tail") return gradFail("forced by CURCUMA_ROCM_FAIL_GRAD=tail");
        if (!ok || nat <= 0) return gradFail("no device");
        const size_t un = static_cast<size_t>(nat);
        if (!rpZ.ensure(un) || !rpXyz.ensure(3 * un) || !rpAlpha.ensure(un) || !rpZeff.ensure(un) || !rpEat.ensure(un))
            return gradFail("device allocation for the repulsion energy failed");
        if (hipMemcpyAsync(rpZ.p, z, un * sizeof(int), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(rpXyz.p, xyz, 3 * un * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(rpAlpha.p, alpha, un * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess
            || hipMemcpyAsync(rpZeff.p, zeff, un * sizeof(double), hipMemcpyHostToDevice, stream) != hipSuccess)
            return gradFail("upload of the repulsion inputs failed");
        const int b = 128;
        hipLaunchKernelGGL(k_rep_energy, dim3((nat + b - 1) / b), dim3(b), 0, stream, nat, is_gfn2 ? 1 : 0, rpZ.p,
                           rpXyz.p, rpAlpha.p, rpZeff.p, kexp, rexp, kexp_light, rpEat.p);
        if (hipGetLastError() != hipSuccess) return gradFail("repulsion kernel failed");
        std::vector<double> e(un);
        if (hipMemcpyAsync(e.data(), rpEat.p, un * sizeof(double), hipMemcpyDeviceToHost, stream) != hipSuccess
            || hipStreamSynchronize(stream) != hipSuccess) return gradFail("download of the repulsion energy failed");
        double E = 0.0;
        for (size_t a = 0; a < un; ++a) E += e[a];
        *energy_out = E;
        return true;
    }
    void freeGradTail() {
        for (auto* b : { &tlXyz, &tlDcn, &tlRcov, &tlGrad, &tlQ, &tlDp, &tlQp, &tlMrad, &tlCn, &tlRad, &tlVcn, &tlDer,
                         &tlGmp, &rpXyz, &rpAlpha, &rpZeff, &rpEat })
            b->release();
        rpZ.release();
    }
    void freePostScf() {
        for (auto* b : { &psQ, &psQ0, &psRad, &psEn, &psW0, &psdWq0, &psdWc0, &psAtmEat, &psAtmGrad, &psAtmDcn,
                         &psGtot, &psDtot, &cnXyz, &cnOut, &d3Xyz, &d3Cn, &d3Gw, &d3Dgw, &d3Eat, &d3Grad, &d3Dcn,
                         &d3C6, &d3R6, &d3Refcn, &d3Rcov })
            b->release();
        for (auto* b : { &psNbCnt, &psNbPtr, &psNb, &cnZ, &d3Z, &d3Type, &d3Nref }) b->release();
    }
    void freeDispersion() {
        for (void* p : { (void*)dD4Z, (void*)dD4Nref, (void*)dD4Sqrt, (void*)dD4Xyz,
                         (void*)dD4W, (void*)dD4dWq, (void*)dD4Dedq, (void*)dD4C6,
                         (void*)dD4dWc, (void*)dD4Eat, (void*)dD4Grad, (void*)dD4Dcn,
                         (void*)dD4c6, (void*)dD4dc6 })
            if (p) hipFree(p);
        dD4Z = dD4Nref = nullptr; dD4Sqrt = dD4Xyz = dD4W = dD4dWq = dD4Dedq = dD4C6 = nullptr;
        dD4dWc = dD4Eat = dD4Grad = dD4Dcn = nullptr;
        dD4c6 = dD4dc6 = nullptr;
        d4_nat = d4_c6_len = 0; d4_c6_uploaded = false;
    }

    // ----- Single-shot D4 EEQ: charges + ∂q/∂x response -------------------
    int eeq_alloc_n = 0;
    void freeEeq() {
        for (void* p : { (void*)dEeqXyz,(void*)dEeqChi,(void*)dEeqGam,(void*)dEeqAlp,(void*)dEeqCnf,
                         (void*)dEeqRcov,(void*)dEeqCn,(void*)dEeqCnRaw,(void*)dEeqM,(void*)dEeqRhs,
                         (void*)dEeqQ,(void*)dEeqIpiv,(void*)dEeqDedq,(void*)dEeqAdjRhs,(void*)dEeqU,
                         (void*)dEeqGrad }) if (p) hipFree(p);
        dEeqXyz = dEeqChi = dEeqGam = dEeqAlp = dEeqCnf = dEeqRcov = dEeqCn = dEeqCnRaw = nullptr;
        dEeqM = dEeqRhs = dEeqQ = dEeqDedq = dEeqAdjRhs = dEeqU = dEeqGrad = nullptr;
        dEeqIpiv = nullptr; eeq_alloc_n = 0; eeq_n = 0;
    }
    bool ensureEeq(int N) {
        if (!handle_ok) {
            if (rocblas_create_handle(&handle) != rocblas_status_success) return false;
            rocblas_set_stream(handle, stream); handle_ok = true;
        }
        if (!dInfo && hipMalloc(&dInfo, sizeof(rocblas_int)) != hipSuccess) return false;
        if (N == eeq_alloc_n && dEeqM) return true;
        freeEeq();
        const size_t m = static_cast<size_t>(N) + 1;
        bool g = true;
        g &= hipMalloc(&dEeqXyz, 3 * static_cast<size_t>(N) * sizeof(double)) == hipSuccess;
        for (double** p : { &dEeqChi,&dEeqGam,&dEeqAlp,&dEeqCnf,&dEeqRcov,&dEeqCn,&dEeqCnRaw,
                            &dEeqQ,&dEeqDedq,&dEeqU })
            g &= hipMalloc(p, static_cast<size_t>(N) * sizeof(double)) == hipSuccess;
        g &= hipMalloc(&dEeqM, m * m * sizeof(double)) == hipSuccess;
        g &= hipMalloc(&dEeqRhs, m * sizeof(double)) == hipSuccess;
        g &= hipMalloc(&dEeqAdjRhs, m * sizeof(double)) == hipSuccess;
        g &= hipMalloc(&dEeqGrad, 3 * static_cast<size_t>(N) * sizeof(double)) == hipSuccess;
        g &= hipMalloc(&dEeqIpiv, m * sizeof(rocblas_int)) == hipSuccess;
        if (!g) { freeEeq(); return false; }
        eeq_alloc_n = N;
        return true;
    }
    // Build the augmented EEQ matrix, LU-factor it (kept for the adjoint), solve for the charges.
    bool eeqChargesImpl(int N, const double* xyz, const double* chi, const double* gam,
                        const double* alpha_sq, const double* cnf, const double* rcov,
                        double total_charge, double* q_out) {
        if (!ok || N <= 0 || !ensureEeq(N)) return false;
        const int m = N + 1;
        if (hipMemcpy(dEeqXyz, xyz, 3 * static_cast<size_t>(N) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dEeqChi, chi, N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dEeqGam, gam, N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dEeqAlp, alpha_sq, N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dEeqCnf, cnf, N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dEeqRcov, rcov, N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        hipLaunchKernelGGL(k_d4eeq_cn, dim3((N + 127) / 128), dim3(128), 0, stream, N, dEeqXyz, dEeqRcov, dEeqCn, dEeqCnRaw);
        dim3 blk(16, 16), grd((m + 15) / 16, (m + 15) / 16);
        hipLaunchKernelGGL(k_d4eeq_build, grd, blk, 0, stream, N, dEeqXyz, dEeqAlp, dEeqGam, dEeqM);
        hipLaunchKernelGGL(k_d4eeq_rhs, dim3((m + 127) / 128), dim3(128), 0, stream, N, dEeqChi, dEeqCnf, dEeqCn, total_charge, dEeqRhs);
        if (hipGetLastError() != hipSuccess) return false;
        if (rocsolver_dgetrf(handle, m, m, dEeqM, m, dEeqIpiv, dInfo) != rocblas_status_success) return false;
        if (rocsolver_dgetrs(handle, rocblas_operation_none, m, 1, dEeqM, m, dEeqIpiv, dEeqRhs, m) != rocblas_status_success) return false;
        if (hipMemcpyAsync(dEeqQ, dEeqRhs, N * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        rocblas_int info = 1;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(q_out, dEeqQ, N * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        eeq_n = (info == 0) ? N : 0;
        return info == 0;
    }
    // ∂q/∂x response: adjoint solve M·z = [dEdq;0] (reuse the LU), then the pairwise gradient.
    bool eeqRespImpl(int N, const double* dEdq, double* grad_add) {
        if (!ok || N <= 0 || N != eeq_n || !dEeqM) return false;
        const int m = N + 1;
        if (N <= 1) { for (int i = 0; i < 3 * N; ++i) grad_add[i] = 0.0; return true; }
        if (hipMemcpy(dEeqDedq, dEdq, N * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        hipLaunchKernelGGL(k_d4eeq_adjoint_rhs, dim3((m + 127) / 128), dim3(128), 0, stream, N, dEeqDedq, dEeqAdjRhs);
        if (hipGetLastError() != hipSuccess) return false;
        if (rocsolver_dgetrs(handle, rocblas_operation_none, m, 1, dEeqM, m, dEeqIpiv, dEeqAdjRhs, m) != rocblas_status_success) return false;
        hipLaunchKernelGGL(k_d4eeq_u, dim3((N + 127) / 128), dim3(128), 0, stream, N, dEeqAdjRhs, dEeqCnf, dEeqCn, dEeqCnRaw, dEeqU);
        hipLaunchKernelGGL(k_d4eeq_response, dim3((N + 127) / 128), dim3(128), 0, stream, N, dEeqXyz, dEeqAlp, dEeqRcov, dEeqQ, dEeqAdjRhs, dEeqU, dEeqGrad);
        if (hipGetLastError() != hipSuccess) return false;
        if (hipMemcpy(grad_add, dEeqGrad, 3 * static_cast<size_t>(N) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return hipStreamSynchronize(stream) == hipSuccess;
    }

    void freeBasis() {
        for (void* p : { (void*)dz,(void*)dsh2at,(void*)dang_sh,(void*)diao_sh,(void*)dnao_sh,
                         (void*)dsh_nprim,(void*)dsh_prim_off,(void*)dvalence,(void*)dprim_alpha,
                         (void*)dprim_coeff,(void*)dsh_zeta,(void*)dselfenergy,(void*)dkcn,
                         (void*)dshpoly,(void*)dghard,(void*)dcovrad,(void*)dpauling,(void*)darad,
                         (void*)dxyz,(void*)dcn,(void*)dse,(void*)dgamma,(void*)dL,
                         (void*)dao2sh,(void*)dao2at,(void*)drep_alpha,(void*)drep_zeff,
                         (void*)dW,(void*)dVao,(void*)dQsh,(void*)dGrad,(void*)dEdcn_buf,
                         (void*)dDpInt,(void*)dQpInt,(void*)dVdp,(void*)dVqp,(void*)dDpAt,(void*)dQpAt })
            if (p) hipFree(p);
        dz = dsh2at = dang_sh = diao_sh = dnao_sh = dsh_nprim = dsh_prim_off = dvalence = nullptr;
        dprim_alpha = dprim_coeff = dsh_zeta = dselfenergy = dkcn = dshpoly = dghard = nullptr;
        dcovrad = dpauling = darad = dxyz = dcn = dse = dgamma = dL = nullptr;
        dao2sh = dao2at = nullptr; drep_alpha = drep_zeff = nullptr;
        dW = dVao = dQsh = dGrad = dEdcn_buf = nullptr;
        dDpInt = dQpInt = dVdp = dVqp = dDpAt = dQpAt = nullptr; mp_computed = false;
        basis_ready = false; bgeom_nat = 0;
    }
#endif // HAVE_ROCSOLVER
};

XtbHipContext::XtbHipContext(int device)
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (hipGetDeviceCount(&count) != hipSuccess || count <= 0)
        return;
    const int dev = device < 0 ? 0 : device;   // Claude Generated (Sep 2026, multi-GPU)
    if (dev >= count || hipSetDevice(dev) != hipSuccess)
        return;

    hipDeviceProp_t prop;
    if (hipGetDeviceProperties(&prop, dev) != hipSuccess)
        return;
    if (hipStreamCreate(&m_impl->stream) != hipSuccess)
        return;

    m_impl->device = dev;
    m_impl->name   = prop.name;
    m_impl->ok     = true;
}

XtbHipContext::~XtbHipContext()
{
    if (!m_impl) return;
    bindDevice();   // free device memory on the device that owns it
#ifdef HAVE_ROCSOLVER
    m_impl->freeMpPot();   // Claude Generated (Oct 2026)
    m_impl->freePot();     // WP1 (Claude Generated, Oct 2026)
    m_impl->freeResLoop(); // WP2 (Claude Generated, Oct 2026)
    m_impl->freePseudo();  // WP7 (Claude Generated, Oct 2026)
    m_impl->freeSolv();    // WP6 (Claude Generated, Oct 2026)
    m_impl->freeEeq();
    m_impl->freePostScf(); // WP8a (Claude Generated, Oct 2026)
    m_impl->freeGradTail(); // WP8b (Claude Generated, Oct 2026)
    m_impl->freeDispersion();
    m_impl->freeBasis();
    m_impl->freeResident();
    m_impl->freeSolver();
    if (m_impl->handle_ok) rocblas_destroy_handle(m_impl->handle);
#endif
    if (m_impl->stream)
        hipStreamDestroy(m_impl->stream);
}

bool XtbHipContext::ok() const { return m_impl && m_impl->ok; }

std::string XtbHipContext::deviceName() const
{
    return m_impl ? m_impl->name : std::string();
}

int XtbHipContext::deviceId() const { return m_impl ? m_impl->device : -1; }

bool XtbHipContext::bindDevice() const
{
    if (!m_impl || m_impl->device < 0) return false;
    return hipSetDevice(m_impl->device) == hipSuccess;
}

bool XtbHipContext::deviceAvailable()
{
    int count = 0;
    return hipGetDeviceCount(&count) == hipSuccess && count > 0;
}

bool XtbHipContext::solveGeneralized(const double* F, const double* S, int n,
                                     double* eps_out, double* C_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)F; (void)S; (void)n; (void)eps_out; (void)C_colmajor;
    return false;   // built without rocSOLVER → caller uses the CPU eigensolver
#else
    if (!m_impl || !m_impl->ok || n <= 0 || !F || !S || !eps_out || !C_colmajor)
        return false;
    Impl& d = *m_impl;
    if (!d.ensureSolver(n)) return false;
    const size_t nn = static_cast<size_t>(n) * n;

    // F C = S C ε  (itype = AX), eigenvectors of the ORIGINAL problem (CᵀSC = I).
    // dF is overwritten with the eigenvectors, dS with its Cholesky factor.
    if (hipMemcpy(d.dF, F, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
    if (hipMemcpy(d.dS, S, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;

    const rocblas_status st = rocsolver_dsygvd(
        d.handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
        n, d.dF, n, d.dS, n, d.dEps, d.dE, d.dInfo);
    if (st != rocblas_status_success) return false;
    if (hipDeviceSynchronize() != hipSuccess) return false;

    rocblas_int info = 0;
    if (hipMemcpy(&info, d.dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
    if (info != 0) return false;   // S not positive-definite or no convergence

    if (hipMemcpy(C_colmajor, d.dF, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
    if (hipMemcpy(eps_out, d.dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
    return true;
#endif
}

// Claude Generated (Oct 2026): resident-solver switches, see the header.
void XtbHipContext::setSolverOptions(bool cached_cholesky, bool density_syrk, bool explicit_linv)
{
#ifndef HAVE_ROCSOLVER
    (void)cached_cholesky; (void)density_syrk; (void)explicit_linv;
#else
    if (!m_impl) return;
    m_impl->opt_explicit_linv = explicit_linv;
    m_impl->opt_cached_chol  = cached_cholesky;
    m_impl->opt_density_syrk = density_syrk;
#endif
}

// WP10 (Claude Generated, Oct 2026): fixed-point gradient/moment accumulation, see the header.
void XtbHipContext::setFixedPoint(bool on)
{
#ifndef HAVE_ROCSOLVER
    (void)on;
#else
    if (m_impl) m_impl->fixed_point = on;
#endif
}

// Claude Generated (Oct 2026): solver fallback bitmask, see the header.
unsigned XtbHipContext::solverFallbackEvents() const
{
#ifndef HAVE_ROCSOLVER
    return 0u;
#else
    if (!m_impl) return 0u;
    return m_impl->fallback_events | (m_impl->fail_inject ? kFallbackInjected : 0u);
#endif
}

bool XtbHipContext::residentBegin(const double* H0, const double* S, int n)
{
#ifndef HAVE_ROCSOLVER
    (void)H0; (void)S; (void)n; return false;
#else
    return m_impl && m_impl->ok && H0 && S && n > 0 && m_impl->resBegin(H0, S, n);
#endif
}
bool XtbHipContext::residentSolve(const double* v_ao, double* eps_out, bool fp32)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)eps_out; (void)fp32; return false;
#else
    return m_impl && m_impl->ok && v_ao && eps_out && m_impl->resSolve(v_ao, eps_out, fp32);
#endif
}
bool XtbHipContext::residentDensity(const double* occ, int ncol, double* pop_ao, double* band)
{
#ifndef HAVE_ROCSOLVER
    (void)occ; (void)ncol; (void)pop_ao; (void)band; return false;
#else
    return m_impl && m_impl->ok && occ && pop_ao && band && m_impl->resDensity(occ, ncol, pop_ao, band);
#endif
}
bool XtbHipContext::residentFinalize(double* P_colmajor, double* C_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)P_colmajor; (void)C_colmajor; return false;
#else
    return m_impl && m_impl->ok && P_colmajor && C_colmajor && m_impl->resFinalize(P_colmajor, C_colmajor);
#endif
}

bool XtbHipContext::beginBasis(const XtbHipBasisData& basis)
{
#ifndef HAVE_ROCSOLVER
    (void)basis; return false;
#else
    return m_impl && m_impl->ok && m_impl->begBasis(basis);
#endif
}
bool XtbHipContext::beginComputed(const double* xyz_bohr)
{
#ifndef HAVE_ROCSOLVER
    (void)xyz_bohr; return false;
#else
    return m_impl && m_impl->ok && xyz_bohr && m_impl->computeIntegrals(xyz_bohr);
#endif
}
bool XtbHipContext::downloadOverlap(double* S_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)S_colmajor; return false;
#else
    return m_impl && m_impl->ok && S_colmajor && m_impl->basis_ready
        && m_impl->dlMat(S_colmajor, m_impl->dS0, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
#endif
}
bool XtbHipContext::downloadH0(double* H0_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)H0_colmajor; return false;
#else
    return m_impl && m_impl->ok && H0_colmajor && m_impl->basis_ready
        && m_impl->dlMat(H0_colmajor, m_impl->dH0, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
#endif
}
bool XtbHipContext::downloadCholesky(double* L_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)L_colmajor; return false;
#else
    return m_impl && m_impl->ok && L_colmajor && m_impl->basis_ready
        && m_impl->dlMat(L_colmajor, m_impl->dL, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
#endif
}
bool XtbHipContext::downloadGamma(double* gamma_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)gamma_colmajor; return false;
#else
    return m_impl && m_impl->ok && gamma_colmajor && m_impl->basis_ready
        && m_impl->dlMat(gamma_colmajor, m_impl->dgamma, static_cast<size_t>(m_impl->bnsh) * m_impl->bnsh);
#endif
}
bool XtbHipContext::gradient(const double* eps, int nocc_orbs, const double* v_ao,
                             const double* q_sh, const double* v_dp, const double* v_qp,
                             double* grad_out, double* dEdcn_out)
{
#ifndef HAVE_ROCSOLVER
    (void)eps; (void)nocc_orbs; (void)v_ao; (void)q_sh; (void)v_dp; (void)v_qp;
    (void)grad_out; (void)dEdcn_out; return false;
#else
    return m_impl && m_impl->ok && eps && v_ao && q_sh && grad_out && dEdcn_out
        && m_impl->gradient(eps, nocc_orbs, v_ao, q_sh, v_dp, v_qp, grad_out, dEdcn_out);
#endif
}
bool XtbHipContext::beginMultipoleComputed()
{
#ifndef HAVE_ROCSOLVER
    return false;
#else
    return m_impl && m_impl->ok && m_impl->beginResidentMP();
#endif
}
bool XtbHipContext::downloadMultipoleInts(double* dp_int3, double* qp_int6)
{
#ifndef HAVE_ROCSOLVER
    (void)dp_int3; (void)qp_int6; return false;
#else
    return m_impl && m_impl->ok && dp_int3 && qp_int6 && m_impl->downloadMP(dp_int3, qp_int6);
#endif
}
bool XtbHipContext::solveMultipole(const double* v_ao, const double* v_dp,
                                   const double* v_qp, double* eps_out, bool fp32)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)v_dp; (void)v_qp; (void)eps_out; (void)fp32; return false;
#else
    return m_impl && m_impl->ok && v_ao && v_dp && v_qp && eps_out
        && m_impl->resSolveMP(v_ao, v_dp, v_qp, eps_out, fp32);
#endif
}
// Claude Generated (Oct 2026): -xtb.rocm_device_multipole entry points (see the header).
bool XtbHipContext::beginMultipolePotential(int nat, const double* const* mats, const double* dkernel,
                                            const double* qkernel, unsigned long long version,
                                            double max_bytes)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)mats; (void)dkernel; (void)qkernel; (void)version; (void)max_bytes; return false;
#else
    if (!m_impl || !m_impl->ok || !mats || !dkernel || !qkernel) return false;
    for (int b = 0; b < 15; ++b) if (!mats[b]) return false;
    return m_impl->mpPotUpload(nat, mats, dkernel, qkernel, version, max_bytes);
#endif
}
std::string XtbHipContext::multipolePotentialError() const
{
#ifndef HAVE_ROCSOLVER
    return "build without rocSOLVER";
#else
    return m_impl ? m_impl->mp_error : std::string();
#endif
}
bool XtbHipContext::multipolePotential(const double* q_at, const double* dp_at, const double* qp_at,
                                       double* v_at, double* v_dp, double* v_qp)
{
#ifndef HAVE_ROCSOLVER
    (void)q_at; (void)dp_at; (void)qp_at; (void)v_at; (void)v_dp; (void)v_qp; return false;
#else
    return m_impl && m_impl->ok && q_at && dp_at && qp_at && v_at && v_dp && v_qp
        && m_impl->mpPotential(q_at, dp_at, qp_at, v_at, v_dp, v_qp);
#endif
}
bool XtbHipContext::multipoleEnergy(const double* q_at, const double* dp_at, const double* qp_at,
                                    double* e_out)
{
#ifndef HAVE_ROCSOLVER
    (void)q_at; (void)dp_at; (void)qp_at; (void)e_out; return false;
#else
    return m_impl && m_impl->ok && q_at && dp_at && qp_at && e_out
        && m_impl->mpEnergy(q_at, dp_at, qp_at, e_out);
#endif
}
bool XtbHipContext::solveMultipoleResident(const double* v_ao, double* eps_out, bool fp32)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)eps_out; (void)fp32; return false;
#else
    return m_impl && m_impl->ok && v_ao && eps_out && m_impl->resSolveMPResident(v_ao, eps_out, fp32);
#endif
}
// WP1 entry points, -xtb.rocm_device_potential (Claude Generated, Oct 2026; see the header).
bool XtbHipContext::beginPotential(int nat, int nsh, const double* amat_sd, const double* amat_dd,
                                   const double* amat_sq, const double* dkernel, const double* qkernel,
                                   const double* gamma3, double max_bytes)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)nsh; (void)amat_sd; (void)amat_dd; (void)amat_sq; (void)dkernel; (void)qkernel;
    (void)gamma3; (void)max_bytes; return false;
#else
    if (!m_impl || !m_impl->ok || !amat_sd || !amat_dd || !amat_sq || !dkernel || !qkernel || !gamma3) return false;
    m_impl->pot_error.clear();
    return m_impl->potBegin(nat, nsh, amat_sd, amat_dd, amat_sq, dkernel, qkernel, gamma3, max_bytes);
#endif
}
bool XtbHipContext::beginPotentialStored(int nat, int nsh, const double* xyz_bohr, const double* mrad,
                                         double dmp3, double dmp5, const double* dkernel,
                                         const double* qkernel, const double* gamma3, double max_bytes)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)nsh; (void)xyz_bohr; (void)mrad; (void)dmp3; (void)dmp5; (void)dkernel; (void)qkernel;
    (void)gamma3; (void)max_bytes; return false;
#else
    if (!m_impl || !m_impl->ok || !xyz_bohr || !mrad || !dkernel || !qkernel || !gamma3) return false;
    m_impl->pot_error.clear();
    return m_impl->potBeginStored(nat, nsh, xyz_bohr, mrad, dmp3, dmp5, dkernel, qkernel, gamma3, max_bytes);
#endif
}

// WP11 (Claude Generated, Oct 2026)
bool XtbHipContext::beginPotentialOnTheFly(int nat, int nsh, const double* xyz_bohr, const double* mrad,
                                           double dmp3, double dmp5, const double* dkernel,
                                           const double* qkernel, const double* gamma3)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)nsh; (void)xyz_bohr; (void)mrad; (void)dmp3; (void)dmp5; (void)dkernel; (void)qkernel;
    (void)gamma3; return false;
#else
    if (!m_impl || !m_impl->ok || !xyz_bohr || !mrad || !dkernel || !qkernel || !gamma3) return false;
    m_impl->pot_error.clear();
    return m_impl->potBeginOtf(nat, nsh, xyz_bohr, mrad, dmp3, dmp5, dkernel, qkernel, gamma3);
#endif
}

bool XtbHipContext::onTheFlyMultipoleActive() const
{
#ifndef HAVE_ROCSOLVER
    return false;
#else
    return m_impl && m_impl->mp_otf && m_impl->pot_nsh > 0;
#endif
}
bool XtbHipContext::solvePotential(const double* q_sh, const double* q_at, const double* dp_at,
                                   const double* qp_at, const double* W, const double* dWq,
                                   double* eps_out, bool fp32)
{
#ifndef HAVE_ROCSOLVER
    (void)q_sh; (void)q_at; (void)dp_at; (void)qp_at; (void)W; (void)dWq; (void)eps_out; (void)fp32; return false;
#else
    if (!m_impl || !m_impl->ok || !q_sh || !q_at || !dp_at || !qp_at || !W || !dWq || !eps_out) return false;
    return m_impl->resSolvePotential(q_sh, q_at, dp_at, qp_at, W, dWq, eps_out, fp32);
#endif
}
std::string XtbHipContext::devicePotentialError() const
{
#ifndef HAVE_ROCSOLVER
    return "build without rocSOLVER";
#else
    return m_impl ? m_impl->pot_error : std::string();
#endif
}
// ---- WP6: implicit solvation on the device (Claude Generated, Oct 2026) ----
bool XtbHipContext::beginSolvation(int nat, const double* born_mat, const double* v_shift, double max_bytes)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)born_mat; (void)v_shift; (void)max_bytes; return false;
#else
    if (!m_impl || !m_impl->ok || !born_mat) return false;
    return m_impl->solvBegin(nat, born_mat, v_shift, max_bytes);
#endif
}
std::string XtbHipContext::solvationError() const
{
#ifndef HAVE_ROCSOLVER
    return "build without rocSOLVER";
#else
    return m_impl ? m_impl->solv_error : std::string();
#endif
}
bool XtbHipContext::solvationFailed() const
{
#ifndef HAVE_ROCSOLVER
    return false;
#else
    return m_impl && m_impl->solv_failed;
#endif
}
// ---- WP2: fully device-resident GFN2 SCF loop (Claude Generated, Oct 2026) ----
bool XtbHipContext::beginDispersionWeights(int nat, const double* cn, const double* gi, const double* zeff,
                                           const double* refcn, const double* refcovcn, const double* refq,
                                           const int* nref)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)cn; (void)gi; (void)zeff; (void)refcn; (void)refcovcn; (void)refq; (void)nref; return false;
#else
    if (!m_impl || !m_impl->ok || !cn || !gi || !zeff || !refcn || !refcovcn || !refq || !nref) return false;
    return m_impl->resBeginWeights(nat, cn, gi, zeff, refcn, refcovcn, refq, nref);
#endif
}
bool XtbHipContext::beginResidentLoop(int nsh, int nat, int nao, double Tele, double n_elec, int nocc_pairs,
                                      const double* q_sh0, const double* dp_at0, const double* qp_at0,
                                      const double* q_at0, const double* n0_sh, const double* n0_at,
                                      double alpha, int max_hist, double w0)
{
#ifndef HAVE_ROCSOLVER
    (void)nsh; (void)nat; (void)nao; (void)Tele; (void)n_elec; (void)nocc_pairs; (void)q_sh0; (void)dp_at0;
    (void)qp_at0; (void)q_at0; (void)n0_sh; (void)n0_at; (void)alpha; (void)max_hist; (void)w0; return false;
#else
    if (!m_impl || !m_impl->ok || nsh <= 0 || nat <= 0 || nao <= 0 || !q_sh0 || !dp_at0 || !qp_at0 || !q_at0
        || !n0_sh || !n0_at) return false;
    return m_impl->resBeginLoop(nsh, nat, nao, Tele, n_elec, nocc_pairs, q_sh0, dp_at0, qp_at0, q_at0, n0_sh,
                                n0_at, alpha, max_hist, w0);
#endif
}
// ---- WP3: fully device-resident GFN1 SCF loop (Claude Generated, Oct 2026) ----
bool XtbHipContext::beginResidentLoopGfn1(int nsh, int nat, int nao, double Tele, double n_elec, int nocc_pairs,
                                          const double* q_sh0, const double* q_at0, const double* n0_sh,
                                          const double* n0_at, const double* gamma3_at, double alpha, int max_hist,
                                          double w0)
{
#ifndef HAVE_ROCSOLVER
    (void)nsh; (void)nat; (void)nao; (void)Tele; (void)n_elec; (void)nocc_pairs; (void)q_sh0; (void)q_at0;
    (void)n0_sh; (void)n0_at; (void)gamma3_at; (void)alpha; (void)max_hist; (void)w0; return false;
#else
    if (!m_impl || !m_impl->ok || nsh <= 0 || nat <= 0 || nao <= 0 || !q_sh0 || !q_at0 || !n0_sh || !n0_at
        || !gamma3_at) return false;
    return m_impl->resBeginLoopGfn1(nsh, nat, nao, Tele, n_elec, nocc_pairs, q_sh0, q_at0, n0_sh, n0_at,
                                    gamma3_at, alpha, max_hist, w0);
#endif
}
bool XtbHipContext::residentScfStep(bool fp32, double* dq, double* e_band, double* e_coulomb, double* e_third,
                                    double* e_multipole)
{
#ifndef HAVE_ROCSOLVER
    (void)fp32; (void)dq; (void)e_band; (void)e_coulomb; (void)e_third; (void)e_multipole; return false;
#else
    if (!m_impl || !m_impl->ok || !dq || !e_band || !e_coulomb || !e_third || !e_multipole) return false;
    return m_impl->resStep(fp32, dq, e_band, e_coulomb, e_third, e_multipole);
#endif
}
// WP7 (Claude Generated, Oct 2026): pseudo-diagonalisation switches and counts, see the header.
void XtbHipContext::setPseudoDiagonalisation(bool on, double max_dq, int max_steps, bool fp64)
{
#ifdef HAVE_ROCSOLVER
    if (!m_impl) return;
    m_impl->ps_enabled = on;
    m_impl->ps_max_dq = max_dq;
    m_impl->ps_max_steps = max_steps;
    m_impl->ps_fp64 = on && fp64;
    m_impl->ps_disabled = false;
    m_impl->ps_notice.clear();
    if (!on) m_impl->freePseudo();
#else
    (void)on; (void)max_dq; (void)max_steps; (void)fp64;
#endif
}

void XtbHipContext::pseudoDiagonalisationCounts(int& taken, int& rejected) const
{
#ifdef HAVE_ROCSOLVER
    taken = m_impl ? m_impl->ps_taken : 0;
    rejected = m_impl ? m_impl->ps_rejected : 0;
#else
    taken = rejected = 0;
#endif
}

std::string XtbHipContext::pseudoDiagonalisationNotice() const
{
#ifdef HAVE_ROCSOLVER
    return m_impl ? m_impl->ps_notice : std::string();
#else
    return {};
#endif
}

bool XtbHipContext::residentLoopCharges(double* q_sh, double* q_at, double* dp_at, double* qp_at, double* eps)
{
#ifndef HAVE_ROCSOLVER
    (void)q_sh; (void)q_at; (void)dp_at; (void)qp_at; (void)eps; return false;
#else
    return m_impl && m_impl->ok && m_impl->resCharges(q_sh, q_at, dp_at, qp_at, eps);
#endif
}
// ---- WP4: open shell on the resident loop (Claude Generated, Oct 2026) ----
bool XtbHipContext::setResidentOpenShell(double n_alpha, double n_beta)
{
#ifndef HAVE_ROCSOLVER
    (void)n_alpha; (void)n_beta; return false;
#else
    return m_impl && m_impl->ok && m_impl->resSetOpenShell(n_alpha, n_beta);
#endif
}
std::string XtbHipContext::residentLoopError() const
{
#ifndef HAVE_ROCSOLVER
    return "build without rocSOLVER";
#else
    return m_impl ? m_impl->res_error : std::string();
#endif
}
// ---- WP5: Broyden seed preconditioner on the resident loop (Claude Generated, Oct 2026) ----
bool XtbHipContext::broydenSetPreconditioner(const BroydenPrecond* pc, unsigned long long tag, BroydenGuard guard,
                                             bool keep, bool best)
{
#ifndef HAVE_ROCSOLVER
    (void)pc; (void)tag; (void)guard; (void)keep; (void)best; return false;
#else
    if (!m_impl || !m_impl->ok) return false;
    bindDevice();
    return m_impl->pcSet(pc, tag, std::move(guard), keep, best);
#endif
}
unsigned long long XtbHipContext::broydenPreconditionerTag() const
{
#ifndef HAVE_ROCSOLVER
    return 0;
#else
    return m_impl ? m_impl->pc_tag : 0;
#endif
}
bool XtbHipContext::broydenPreconditionerActive() const
{
#ifndef HAVE_ROCSOLVER
    return false;
#else
    return m_impl && m_impl->pc_active;
#endif
}
int XtbHipContext::broydenGuardFiredAt() const
{
#ifndef HAVE_ROCSOLVER
    return -1;
#else
    return m_impl ? m_impl->pc_fired_at : -1;
#endif
}
bool XtbHipContext::multipoleMoments(double* dp_at, double* qp_at)
{
#ifndef HAVE_ROCSOLVER
    (void)dp_at; (void)qp_at; return false;
#else
    return m_impl && m_impl->ok && dp_at && qp_at && m_impl->multipoleMomentsImpl(dp_at, qp_at);
#endif
}
bool XtbHipContext::beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                                    const int* nref, const double* xyz_bohr,
                                    const double* c6_flat, int c6_flat_len,
                                    double s6, double s8, double a1, double a2, double cutoff)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)Z; (void)sqrtZr4r2; (void)nref; (void)xyz_bohr; (void)c6_flat;
    (void)c6_flat_len; (void)s6; (void)s8; (void)a1; (void)a2; (void)cutoff; return false;
#else
    return m_impl && m_impl->ok && Z && sqrtZr4r2 && nref && xyz_bohr && c6_flat
        && m_impl->begDispersion(nat, Z, sqrtZr4r2, nref, xyz_bohr, c6_flat, c6_flat_len,
                                 s6, s8, a1, a2, cutoff);
#endif
}
bool XtbHipContext::dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)W; (void)dWq; (void)dEdq_out; return false;
#else
    return m_impl && m_impl->ok && W && dWq && dEdq_out && m_impl->dispDedq(nat, W, dWq, dEdq_out);
#endif
}
bool XtbHipContext::dispersionGradient(int nat, const double* W, const double* dWq, const double* dWc,
                                       double* e_atom_out, double* grad_out,
                                       double* dEdcn_out, double* dEdq_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)W; (void)dWq; (void)dWc; (void)e_atom_out; (void)grad_out;
    (void)dEdcn_out; (void)dEdq_out; return false;
#else
    return m_impl && m_impl->ok && W && dWq && dWc && e_atom_out && grad_out && dEdcn_out && dEdq_out
        && m_impl->dispDispGrad(nat, W, dWq, dWc, e_atom_out, grad_out, dEdcn_out, dEdq_out);
#endif
}
bool XtbHipContext::dispersionATM(int nat, const double* c6, const double* dc6dcn,
                                  double s9, double a1, double a2, double alp, double cutoff,
                                  double* e_atom_out, double* grad_out, double* dEdcn_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)c6; (void)dc6dcn; (void)s9; (void)a1; (void)a2; (void)alp; (void)cutoff;
    (void)e_atom_out; (void)grad_out; (void)dEdcn_out; return false;
#else
    return m_impl && m_impl->ok && c6 && dc6dcn && e_atom_out && grad_out && dEdcn_out
        && m_impl->dispATM(nat, c6, dc6dcn, s9, a1, a2, alp, cutoff, e_atom_out, grad_out, dEdcn_out);
#endif
}
// ---- WP8a: post-SCF dispersion and CN on the device (Claude Generated, Oct 2026) ----
bool XtbHipContext::d4PostScfDevice(int nat, const double* q_at, bool with_atm, double s9, double a1, double a2,
                                    double alp, double atm_cutoff, double cn_cutoff, const double* rad,
                                    const double* en, double kn, double k4, double k5, double k6,
                                    double* energy_out, double* grad_out, double* dEdq_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)q_at; (void)with_atm; (void)s9; (void)a1; (void)a2; (void)alp; (void)atm_cutoff;
    (void)cn_cutoff; (void)rad; (void)en; (void)kn; (void)k4; (void)k5; (void)k6; (void)energy_out;
    (void)grad_out; (void)dEdq_out; return false;
#else
    if (!m_impl || !q_at || !rad || !en || !energy_out || !grad_out || !dEdq_out) return false;
    return m_impl->d4PostScf(nat, q_at, with_atm, s9, a1, a2, alp, atm_cutoff, cn_cutoff, rad, en, kn, k4, k5, k6,
                             energy_out, grad_out, dEdq_out);
#endif
}
bool XtbHipContext::d3DispersionDevice(int nat, const int* Z, const double* xyz_ang, int nt, const int* atom_type,
                                       const double* c6tab, const double* r6tab, const int* nref_tab,
                                       const double* refcn_tab, const double* rcov_ang, int nrcov, double s6,
                                       double s8, double a1, double a2, double k1, double k2, double cn_cut_ang,
                                       double bohr_per_ang, double fold_unit, bool need_gradient,
                                       double* energy_out, double* grad_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)Z; (void)xyz_ang; (void)nt; (void)atom_type; (void)c6tab; (void)r6tab; (void)nref_tab;
    (void)refcn_tab; (void)rcov_ang; (void)nrcov; (void)s6; (void)s8; (void)a1; (void)a2; (void)k1; (void)k2;
    (void)cn_cut_ang; (void)bohr_per_ang; (void)fold_unit; (void)need_gradient; (void)energy_out; (void)grad_out;
    return false;
#else
    if (!m_impl || !Z || !xyz_ang || !atom_type || !c6tab || !r6tab || !nref_tab || !refcn_tab || !rcov_ang
        || nrcov <= 0 || !energy_out || (need_gradient && !grad_out)) return false;
    return m_impl->d3Dispersion(nat, Z, xyz_ang, nt, atom_type, c6tab, r6tab, nref_tab, refcn_tab, rcov_ang, nrcov,
                                s6, s8, a1, a2, k1, k2, cn_cut_ang, bohr_per_ang, fold_unit, need_gradient,
                                energy_out, grad_out);
#endif
}
bool XtbHipContext::d4CovalentCnDevice(int nat, const int* Z, const double* xyz_bohr, const double* rad,
                                       const double* en, double kn, double k4, double k5, double k6, double cutoff,
                                       double* cn_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)Z; (void)xyz_bohr; (void)rad; (void)en; (void)kn; (void)k4; (void)k5; (void)k6; (void)cutoff;
    (void)cn_out; return false;
#else
    if (!m_impl || !Z || !xyz_bohr || !rad || !en || !cn_out) return false;
    return m_impl->d4CovalentCn(nat, Z, xyz_bohr, rad, en, kn, k4, k5, k6, cutoff, cn_out);
#endif
}
bool XtbHipContext::downloadCn(double* cn_out)
{
#ifndef HAVE_ROCSOLVER
    (void)cn_out; return false;
#else
    return m_impl && cn_out && m_impl->dlCn(cn_out);
#endif
}
std::string XtbHipContext::postScfError() const
{
#ifndef HAVE_ROCSOLVER
    return std::string();
#else
    return m_impl ? m_impl->pst_error : std::string();
#endif
}
bool XtbHipContext::postScfUnavailable() const
{
#ifndef HAVE_ROCSOLVER
    return true;
#else
    return m_impl && m_impl->pst_unavailable;
#endif
}
double XtbHipContext::postScfAtmMs() const
{
#ifndef HAVE_ROCSOLVER
    return -1.0;
#else
    return m_impl ? m_impl->pst_atm_ms : -1.0;
#endif
}
// ---- WP8b: gradient potentials, gradient tail and repulsion on the device (Claude Generated, Oct 2026) ----
bool XtbHipContext::convergedPotentials(bool gfn2)
{
#ifndef HAVE_ROCSOLVER
    (void)gfn2; return false;
#else
    return m_impl && m_impl->convPotentials(gfn2);
#endif
}
bool XtbHipContext::gradientResident(const double* eps, int nocc_orbs, double* grad_out, double* dEdcn_out)
{
#ifndef HAVE_ROCSOLVER
    (void)eps; (void)nocc_orbs; (void)grad_out; (void)dEdcn_out; return false;
#else
    if (!m_impl || !m_impl->ok || !eps || !grad_out || !dEdcn_out) return false;
    if (m_impl->gradient(eps, nocc_orbs, nullptr, nullptr, nullptr, nullptr, grad_out, dEdcn_out, true)) return true;
    if (m_impl->grad_error.empty()) m_impl->grad_error = "device gradient with the resident potentials failed";
    return false;
#endif
}
bool XtbHipContext::gradientTailDevice(int nat, bool is_gfn2, bool with_mp, const double* xyz_bohr,
                                       const double* dEdcn, const double* rcov_au, const double* q_at,
                                       const double* dp_at, const double* qp_at, const double* mrad,
                                       const double* cn, const double* rad, const double* vcn, double dmp3,
                                       double dmp5, double shift, double kexp, double rmax, double* grad_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)is_gfn2; (void)with_mp; (void)xyz_bohr; (void)dEdcn; (void)rcov_au; (void)q_at; (void)dp_at;
    (void)qp_at; (void)mrad; (void)cn; (void)rad; (void)vcn; (void)dmp3; (void)dmp5; (void)shift; (void)kexp;
    (void)rmax; (void)grad_out; return false;
#else
    if (!m_impl || !xyz_bohr || !dEdcn || !rcov_au || !grad_out
        || (with_mp && (!q_at || !dp_at || !qp_at || !mrad || !cn || !rad || !vcn))) return false;
    return m_impl->gradTail(nat, is_gfn2, with_mp, xyz_bohr, dEdcn, rcov_au, q_at, dp_at, qp_at, mrad, cn, rad, vcn,
                            dmp3, dmp5, shift, kexp, rmax, grad_out);
#endif
}
bool XtbHipContext::repulsionEnergyDevice(int nat, bool is_gfn2, const int* z, const double* xyz_bohr,
                                          const double* alpha, const double* zeff, double kexp, double rexp,
                                          double kexp_light, double* energy_out)
{
#ifndef HAVE_ROCSOLVER
    (void)nat; (void)is_gfn2; (void)z; (void)xyz_bohr; (void)alpha; (void)zeff; (void)kexp; (void)rexp;
    (void)kexp_light; (void)energy_out; return false;
#else
    if (!m_impl || !z || !xyz_bohr || !alpha || !zeff || !energy_out) return false;
    return m_impl->repEnergy(nat, is_gfn2, z, xyz_bohr, alpha, zeff, kexp, rexp, kexp_light, energy_out);
#endif
}
std::string XtbHipContext::gradientPathError() const
{
#ifndef HAVE_ROCSOLVER
    return std::string();
#else
    return m_impl ? m_impl->grad_error : std::string();
#endif
}
bool XtbHipContext::eeqCharges(int N, const double* xyz_bohr, const double* chi, const double* gam,
                              const double* alpha_sq, const double* cnf, const double* rcov_bohr,
                              double total_charge, double* q_out)
{
#ifndef HAVE_ROCSOLVER
    (void)N; (void)xyz_bohr; (void)chi; (void)gam; (void)alpha_sq; (void)cnf; (void)rcov_bohr;
    (void)total_charge; (void)q_out; return false;
#else
    return m_impl && m_impl->ok && xyz_bohr && chi && gam && alpha_sq && cnf && rcov_bohr && q_out
        && m_impl->eeqChargesImpl(N, xyz_bohr, chi, gam, alpha_sq, cnf, rcov_bohr, total_charge, q_out);
#endif
}
bool XtbHipContext::eeqChargeResponseGradient(int N, const double* dEdq, double* grad_add)
{
#ifndef HAVE_ROCSOLVER
    (void)N; (void)dEdq; (void)grad_add; return false;
#else
    return m_impl && m_impl->ok && dEdq && grad_add && m_impl->eeqRespImpl(N, dEdq, grad_add);
#endif
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM
