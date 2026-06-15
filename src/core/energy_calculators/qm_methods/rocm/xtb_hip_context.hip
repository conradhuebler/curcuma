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

#ifdef USE_ROCM_XTB

#include "xtb_hip_context.h"

#include <hip/hip_runtime.h>
#ifdef HAVE_ROCSOLVER
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include "xtb_hip_integrals.hiph"                 // CN/overlap/H0/gamma device helpers
#include "../parameters/xtb_params_extra.hpp"     // covalent_rad_d3_au / pauling_en / atomic_rad_au
#endif

#include <memory>
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
    const double* __restrict__ atomic_rad, const double* __restrict__ pauling,
    double* __restrict__ grad, double* __restrict__ dEdcn)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    const int nu = blockIdx.y * blockDim.y + threadIdx.y;
    if (mu >= nao || nu >= nao) return;
    const int iat = ao2at[mu], jat = ao2at[nu];
    if (iat >= jat) return;

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
    const int ta = d_ao_to_type(ang_sh[isha], mu - iao_sh[isha]);
    const int tb = d_ao_to_type(ang_sh[ishb], nu - iao_sh[ishb]);
    if (ta < 0 || tb < 0) return;

    double Sx, D[3], Q[6];
    d_cgto_multipole(prim_alpha + sh_prim_off[isha], prim_coeff + sh_prim_off[isha], sh_nprim[isha],
                     prim_alpha + sh_prim_off[ishb], prim_coeff + sh_prim_off[ishb], sh_nprim[ishb],
                     xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                     xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], ta, tb, Sx, D, Q);

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
                                    const int* ao2at, int n)
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
    }

    bool resBegin(const double* H0, const double* S, int n) {
        if (!ok || n <= 0 || !ensureResident(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dH0, H0, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dS0, S,  nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        return true;
    }
    bool resSolve(const double* v_ao, double* eps_out) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);   // F → dF
        if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        // F C = S C ε on the device; dF ← eigenvectors (ascending), dS destroyed.
        if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                             n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;
        if (hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;   // dF holds C resident, ascending order
    }
    bool resDensity(const double* occ, int ncol, double* pop_ao, double* band) {
        const int n = res_n; if (n <= 0 || ncol < 0 || ncol > n) return false;
        std::vector<double> occ_full(n, 0.0);          // ascending order matches C columns
        for (int r = 0; r < ncol; ++r) occ_full[r] = occ[r];
        if (hipMemcpy(dOcc, occ_full.data(), n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_scale_cols, grd, blk, 0, stream, dF, dOcc, dCw, n);   // Cw = C·diag(occ)
        const double one = 1.0, zero = 0.0;
        // P = Cw · Cᵀ
        if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                          n, n, n, &one, dCw, n, dF, n, &zero, dP, n) != rocblas_status_success) return false;
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

    template <class T> bool up(T*& dptr, const T* host, size_t count) {
        if (dptr) hipFree(dptr);
        dptr = nullptr;
        if (hipMalloc(&dptr, count * sizeof(T)) != hipSuccess) return false;
        return hipMemcpy(dptr, host, count * sizeof(T), hipMemcpyHostToDevice) == hipSuccess;
    }

    bool begBasis(const XtbHipBasisData& bd) {
        if (!ok || bd.nao <= 0 || !ensureResident(bd.nao)) return false;  // allocates dH0/dS0/...
        bnat = bd.nat; bnsh = bd.nsh; bnao = bd.nao; bnprim = bd.nprim_total; bis_gfn2 = bd.is_gfn2;
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

    bool gradient(const double* eps, int nocc_orbs, const double* v_ao, const double* q_sh,
                  double* grad_out, double* dEdcn_out) {
        if (!basis_ready || bnao <= 0 || !dP) return false;
        const int nao = bnao, nat = bnat, nsh = bnsh;
        if (hipMemcpy(dVao, v_ao, nao * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dQsh, q_sh, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        // Energy-weighted density W = C·diag(2·ε_occ)·Cᵀ (reuse dOcc/dCw; dF = C resident).
        std::vector<double> occ2(nao, 0.0);
        for (int k = 0; k < nocc_orbs && k < nao; ++k) occ2[k] = 2.0 * eps[k];
        if (hipMemcpy(dOcc, occ2.data(), nao * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        { dim3 b(16, 16), gr((nao + 15) / 16, (nao + 15) / 16);
          hipLaunchKernelGGL(k_scale_cols, gr, b, 0, stream, dF, dOcc, dCw, nao); }
        const double one = 1.0, zero = 0.0;
        if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                          nao, nao, nao, &one, dCw, nao, dF, nao, &zero, dW, nao) != rocblas_status_success) return false;
        if (hipMemsetAsync(dGrad, 0, 3 * nat * sizeof(double), stream) != hipSuccess) return false;
        if (hipMemsetAsync(dEdcn_buf, 0, nat * sizeof(double), stream) != hipSuccess) return false;
        const double kexp = 1.5, rexp = 1.0, kexp_light = bis_gfn2 ? 1.0 : 1.5;
        { dim3 b(128), gr((nat + 127) / 128);
          hipLaunchKernelGGL(k_grad_repulsion, gr, b, 0, stream, nat, bis_gfn2, dz, dxyz,
                             drep_alpha, drep_zeff, kexp, rexp, kexp_light, dGrad); }
        { dim3 b(128), gr((nao + 127) / 128);
          hipLaunchKernelGGL(k_grad_cn_onsite, gr, b, 0, stream, nao, dP, dkcn, dao2sh, dao2at, dEdcn_buf); }
        { dim3 b(16, 16), gr((nao + 15) / 16, (nao + 15) / 16);
          hipLaunchKernelGGL(k_grad_h0_pulay, gr, b, 0, stream, nao, bis_gfn2,
                             dao2sh, dao2at, dang_sh, diao_sh, dsh_nprim, dsh_prim_off,
                             dprim_alpha, dprim_coeff, dsh_zeta, dshpoly, dkcn, dvalence, dz, dse,
                             dxyz, dP, dS0, dH0, dW, dVao, darad, dpauling, dGrad, dEdcn_buf); }
        { dim3 b(128), gr((nsh + 127) / 128);
          hipLaunchKernelGGL(k_grad_coulomb, gr, b, 0, stream, nsh, bis_gfn2, dsh2at, dghard,
                             dQsh, dxyz, 2.0, dGrad); }
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (hipMemcpy(grad_out, dGrad, 3 * nat * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(dEdcn_out, dEdcn_buf, nat * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;
    }

    bool computeIntegrals(const double* xyz) {
        if (!basis_ready || bnao <= 0) return false;
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
        // GFN2: build dp_int/qp_int on the device too (resident overlap dS0 just built),
        // so they are resident for the multipole Fock/moments + downloadable for the host gradient.
        if (bis_gfn2 && !computeMultipole()) return false;
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
    bool resSolveMP(const double* v_ao, const double* v_dp, const double* v_qp, double* eps_out) {
        const int n = res_n; if (n <= 0 || !beginResidentMP()) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dVdp, v_dp, 3 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dVqp, v_qp, 6 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);              // F → dF
        hipLaunchKernelGGL(k_add_fock_multipole, grd, blk, 0, stream, dF, dDpInt, dQpInt,
                           dVdp, dVqp, dao2at, n);                                            // F −= ½ multipole
        if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                             n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;
        return hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
    }
    // Atomic dp_at(3×nat)/qp_at(6×nat) from the resident density dP (k_multipole_moments).
    bool multipoleMomentsImpl(double* dp_at, double* qp_at) {
        if (!beginResidentMP() || bnat <= 0) return false;
        if (hipMemsetAsync(dDpAt, 0, 3 * static_cast<size_t>(bnat) * sizeof(double), stream) != hipSuccess) return false;
        if (hipMemsetAsync(dQpAt, 0, 6 * static_cast<size_t>(bnat) * sizeof(double), stream) != hipSuccess) return false;
        dim3 b(128), gr((bnao + 127) / 128);
        hipLaunchKernelGGL(k_multipole_moments, gr, b, 0, stream, dDpAt, dQpAt, dP, dDpInt, dQpInt, dao2at, bnao);
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        if (hipMemcpy(dp_at, dDpAt, 3 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return hipMemcpy(qp_at, dQpAt, 6 * static_cast<size_t>(bnat) * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess;
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

XtbHipContext::XtbHipContext()
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (hipGetDeviceCount(&count) != hipSuccess || count <= 0)
        return;
    if (hipSetDevice(0) != hipSuccess)
        return;

    hipDeviceProp_t prop;
    if (hipGetDeviceProperties(&prop, 0) != hipSuccess)
        return;
    if (hipStreamCreate(&m_impl->stream) != hipSuccess)
        return;

    m_impl->device = 0;
    m_impl->name   = prop.name;
    m_impl->ok     = true;
}

XtbHipContext::~XtbHipContext()
{
    if (!m_impl) return;
#ifdef HAVE_ROCSOLVER
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

bool XtbHipContext::residentBegin(const double* H0, const double* S, int n)
{
#ifndef HAVE_ROCSOLVER
    (void)H0; (void)S; (void)n; return false;
#else
    return m_impl && m_impl->ok && H0 && S && n > 0 && m_impl->resBegin(H0, S, n);
#endif
}
bool XtbHipContext::residentSolve(const double* v_ao, double* eps_out)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)eps_out; return false;
#else
    return m_impl && m_impl->ok && v_ao && eps_out && m_impl->resSolve(v_ao, eps_out);
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
                             const double* q_sh, double* grad_out, double* dEdcn_out)
{
#ifndef HAVE_ROCSOLVER
    (void)eps; (void)nocc_orbs; (void)v_ao; (void)q_sh; (void)grad_out; (void)dEdcn_out; return false;
#else
    return m_impl && m_impl->ok && eps && v_ao && q_sh && grad_out && dEdcn_out
        && m_impl->gradient(eps, nocc_orbs, v_ao, q_sh, grad_out, dEdcn_out);
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
                                   const double* v_qp, double* eps_out)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)v_dp; (void)v_qp; (void)eps_out; return false;
#else
    return m_impl && m_impl->ok && v_ao && v_dp && v_qp && eps_out
        && m_impl->resSolveMP(v_ao, v_dp, v_qp, eps_out);
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

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM_XTB
