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
    const double* __restrict__ v_dp, const double* __restrict__ v_qp,
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

// FP32 mixed-precision (X-AP3): element-wise FP64<->FP32 conversion for the optional
// rocsolver_ssygvd eigensolve path (far-from-convergence iterations of -scf_mixed_precision).
__global__ void k_d2f(const double* __restrict__ src, float* __restrict__ dst, size_t n)
{ size_t i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) dst[i] = static_cast<float>(src[i]); }
__global__ void k_f2d(const float* __restrict__ src, double* __restrict__ dst, size_t n)
{ size_t i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) dst[i] = static_cast<double>(src[i]); }

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
    // FP32 generalized eigensolve: dF/dS0 (FP64) → ssygvd → dF (C, FP64) + dEps (FP64).
    // dS0 stays intact (resident overlap). Returns false on any error → caller can retry FP64.
    bool eigensolveFp32(int n) {
        if (!ensureFp32(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        const int blk = 256; const int g = static_cast<int>((nn + blk - 1) / blk);
        hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, dF, dF32, nn);    // F → FP32
        hipLaunchKernelGGL(k_d2f, dim3(g), dim3(blk), 0, stream, dS0, dS32, nn);   // S → FP32 (work)
        if (rocsolver_ssygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                             n, dF32, n, dS32, n, dEps32, dE32, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;
        const int gn = static_cast<int>((nn + blk - 1) / blk), gn1 = (n + blk - 1) / blk;
        hipLaunchKernelGGL(k_f2d, dim3(gn), dim3(blk), 0, stream, dF32, dF, nn);   // C → FP64 (dF)
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
    }

    bool resBegin(const double* H0, const double* S, int n) {
        if (!ok || n <= 0 || !ensureResident(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dH0, H0, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dS0, S,  nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        return true;
    }
    // fp32=true (X-AP3): solve this step with rocsolver_ssygvd (FP32) — far-from-convergence
    // -scf_mixed_precision iterations. dF holds C (FP64) + dEps the eigenvalues either way.
    bool resSolve(const double* v_ao, double* eps_out, bool fp32) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);   // F → dF
        if (fp32) {
            if (!eigensolveFp32(n)) return false;   // dF ← C (FP64), dEps ← eps (FP64)
        } else {
            if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
            if (hipStreamSynchronize(stream) != hipSuccess) return false;
            // F C = S C ε on the device; dF ← eigenvectors (ascending), dS destroyed.
            if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                                 n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
            if (hipDeviceSynchronize() != hipSuccess) return false;
            rocblas_int info = 0;
            if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
            if (info != 0) return false;
        }
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
                  const double* v_dp, const double* v_qp,
                  double* grad_out, double* dEdcn_out) {
        if (!basis_ready || bnao <= 0 || !dP) return false;
        const int nao = bnao, nat = bnat, nsh = bnsh;
        if (hipMemcpy(dVao, v_ao, nao * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dQsh, q_sh, nsh * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        // GFN2 (R-AP3): upload the converged multipole potentials for the multipole-integral
        // Pulay term inside k_grad_h0_pulay (column-major 3×nat / 6×nat, matching v_dp[k+iat*3]).
        const bool mp_grad = bis_gfn2 && v_dp && v_qp && dVdp && dVqp;
        if (mp_grad) {
            if (hipMemcpy(dVdp, v_dp, 3 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
            if (hipMemcpy(dVqp, v_qp, 6 * static_cast<size_t>(nat) * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        }
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
                             dxyz, dP, dS0, dH0, dW, dVao,
                             mp_grad ? dVdp : nullptr, mp_grad ? dVqp : nullptr,
                             darad, dpauling, dGrad, dEdcn_buf); }
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
    bool resSolveMP(const double* v_ao, const double* v_dp, const double* v_qp, double* eps_out, bool fp32) {
        const int n = res_n; if (n <= 0 || !beginResidentMP()) return false;
        const size_t nn = static_cast<size_t>(n) * n;
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
            if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
            if (hipStreamSynchronize(stream) != hipSuccess) return false;
            if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                                 n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
            if (hipDeviceSynchronize() != hipSuccess) return false;
            rocblas_int info = 0;
            if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
            if (info != 0) return false;
        }
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

    // ----- Stage 5 (Part B2): in-SCF GFN2 D4 dE_D4/dq ---------------------
    // begDispersion uploads the geometry-fixed reference data once per geometry; the
    // element C6 block (c6_flat) is uploaded only on the first call (process lifetime).
    bool begDispersion(int nat, const int* Z, const double* sqrtZr4r2, const int* nref,
                       const double* xyz, const double* c6_flat, int c6_len,
                       double s6, double s8, double a1, double a2, double cutoff) {
        if (!ok || nat <= 0 || c6_len <= 0) return false;
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
    m_impl->freeEeq();
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

#endif // USE_ROCM_XTB
