/*
 * < CUDA 4-centre ERI build of the native QM engine (Q-G2) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, Q-G2 of docs/QM_GPU_ROADMAP.md).
//
// Builds the stored ERI tensor (mu nu | lam sig) of the native QM engine directly in device
// memory and hands it to the Q-G1 J/K builder (QmGpuJkCuda), so the tensor never crosses
// the PCIe bus. Reference: qmint::buildERI in qm_integrals.cpp; the host flattens the very
// same shell-pair tables (qmint::buildEriGpuTables), so shells, primitive-pair order,
// Hermite coefficients and both Schwarz screens are identical to the CPU.
//
// McMurchie-Davidson (Helgaker, Jorgensen, Olsen, ch. 9.9; McMurchie, Davidson,
// J. Comput. Phys. 26, 218 (1978)):
//   (ab|cd) = sum_prim c_a c_b c_c c_d 2 pi^(5/2) / (p q sqrt(p+q))
//             sum_{tuv} E^ab_t E^ab_u E^ab_v sum_{tau nu phi} (-1)^(tau+nu+phi)
//             E^cd_tau E^cd_nu E^cd_phi R^0_{t+tau, u+nu, v+phi}(rho, P-Q)
//
// Work decomposition (deliberately simple, correctness first):
//   * one thread block per canonical, Schwarz-surviving shell quartet; one launch per
//     angular-momentum class (LA,LB,LC,LD), so the dynamic shared memory is sized per class;
//   * per primitive quartet: thread 0 evaluates the Boys function (same grid + Taylor table
//     as the CPU), the block builds R^0 on the simplex t+u+v <= L level by level in shared
//     memory, then every thread accumulates its own cartesian component quartets in FP64
//     (each component owned by exactly one thread: no atomics, deterministic);
//   * the cartesian -> spherical 6d->5d transform is applied per quartet in shared memory,
//     one axis at a time (as transformAxis on the CPU);
//   * every tensor element is written by exactly one component of exactly one quartet
//     (canonical representative of its 8-fold orbit), then copied to all 8 permutations.
// Limits: s, p, d shells (L <= 2), checked on the host; FP64 throughout.

#include "qm_gpu_jk_cuda.h"

#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

constexpr double kPI = 3.14159265358979323846264338327950288;
constexpr int kNTaylor = 7;  // must match BoysTable::NTAYLOR (checked at creation)

struct DevTables {
    int n_out;
    double prim_thr, pref;
    const int *sh_ncomp, *sh_nprim, *sh_lmn_off, *sh_coef_off, *sh_out_first, *sh_T_off;
    const int* lmn;
    const double *coef, *T;
    const int *pr_A, *pr_B, *pr_prim_off, *pr_nprim, *pp_ia, *pp_ib;
    const double *pp_p, *pp_P, *pp_schwarz;
    const long long* pp_E_off;
    const double* E;
    const int *q_bra, *q_ket;
    const double* boys_f;
    int boys_ntab;
    double boys_dT, boys_Tmax;
};

// Boys function F_0..F_maxN, the same arithmetic as boysArrayInto (qm_integrals.cpp).
__device__ void boysDev(int maxN, double T, double* F, const DevTables& d)
{
    if (T >= d.boys_Tmax) {
        const double eT = exp(-T);
        F[0] = 0.5 * sqrt(kPI / T);
        for (int n = 0; n < maxN; ++n) F[n + 1] = ((2.0 * n + 1.0) * F[n] - eT) / (2.0 * T);
        return;
    }
    const int i = (int)(T / d.boys_dT + 0.5);
    const double dd = i * d.boys_dT - T;
    const double* fi = d.boys_f + (size_t)i * (d.boys_ntab + 1) + maxN;
    const double inv[kNTaylor] = { 1.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4, 1.0 / 5, 1.0 / 6 };
    double acc = fi[kNTaylor - 1];
    for (int k = kNTaylor - 2; k >= 0; --k) acc = fi[k] + dd * inv[k + 1] * acc;
    F[maxN] = acc;
    if (maxN > 0) {
        const double eT = exp(-T), T2 = 2.0 * T;
        for (int n = maxN - 1; n >= 0; --n) F[n] = (T2 * F[n + 1] + eT) / (2.0 * n + 1.0);
    }
}

__device__ inline void set8(double* out, size_t n, int a, int b, int c, int e, double v)
{
    const size_t n2 = n * n;
    out[((size_t)a * n + b) * n2 + (size_t)c * n + e] = v;
    out[((size_t)b * n + a) * n2 + (size_t)c * n + e] = v;
    out[((size_t)a * n + b) * n2 + (size_t)e * n + c] = v;
    out[((size_t)b * n + a) * n2 + (size_t)e * n + c] = v;
    out[((size_t)c * n + e) * n2 + (size_t)a * n + b] = v;
    out[((size_t)e * n + c) * n2 + (size_t)a * n + b] = v;
    out[((size_t)c * n + e) * n2 + (size_t)b * n + a] = v;
    out[((size_t)e * n + c) * n2 + (size_t)b * n + a] = v;
}

// One block per shell quartet of one class. Dynamic shared memory:
//   2 * s^3 (R levels, s = L+1) + 2 * ncart (component block + transform scratch) doubles.
__global__ void eriClassKernel(DevTables d, long long q0, int LA, int LB, int LC, int LD,
                               double* __restrict__ out)
{
    extern __shared__ double smem[];
    __shared__ double sF[16], sBp[16], sW[3], sPref;

    const long long q = q0 + blockIdx.x;
    const int i = d.q_bra[q], j = d.q_ket[q];
    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const int nA = d.sh_ncomp[A], nB = d.sh_ncomp[B], nC = d.sh_ncomp[C], nD = d.sh_ncomp[D];
    const int L = LA + LB + LC + LD;
    const int s = L + 1, s2 = s * s, s3 = s2 * s;
    const int ncart = nA * nB * nC * nD;
    double* Rb[2] = { smem, smem + s3 };
    double* blk = smem + 2 * s3;
    double* tmp = blk + ncart;
    const int tid = threadIdx.x, bs = blockDim.x;

    for (int k = tid; k < ncart; k += bs) blk[k] = 0.0;

    const int sAB_t = (LA + 1) * (LB + 1), sAB_i = LB + 1, szAB = (LA + LB + 1) * sAB_t;
    const int sCD_t = (LC + 1) * (LD + 1), sCD_i = LD + 1, szCD = (LC + LD + 1) * sCD_t;
    const int pb0 = d.pr_prim_off[i], nbp = d.pr_nprim[i];
    const int pk0 = d.pr_prim_off[j], nkp = d.pr_nprim[j];
    const double ket_max = nkp > 0 ? d.pp_schwarz[pk0] : 0.0;
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B], npC = d.sh_nprim[C], npD = d.sh_nprim[D];

    for (int ib = 0; ib < nbp; ++ib) {
        const int pb = pb0 + ib;
        const double sb = d.pp_schwarz[pb];
        if (sb * ket_max < d.prim_thr) break;  // same primitive screen as shellQuartet
        for (int ik = 0; ik < nkp; ++ik) {
            const int pk = pk0 + ik;
            if (sb * d.pp_schwarz[pk] < d.prim_thr) break;
            __syncthreads();  // previous component loop is done with R
            if (tid == 0) {
                const double p = d.pp_p[pb], qq = d.pp_p[pk];
                const double rho = p * qq / (p + qq);
                const double Wx = d.pp_P[3 * pb] - d.pp_P[3 * pk];
                const double Wy = d.pp_P[3 * pb + 1] - d.pp_P[3 * pk + 1];
                const double Wz = d.pp_P[3 * pb + 2] - d.pp_P[3 * pk + 2];
                boysDev(L, rho * (Wx * Wx + Wy * Wy + Wz * Wz), sF, d);
                const double b = -2.0 * rho;
                sBp[0] = 1.0;
                for (int n = 1; n <= L; ++n) sBp[n] = sBp[n - 1] * b;
                sW[0] = Wx; sW[1] = Wy; sW[2] = Wz;
                sPref = d.pref / (p * qq * sqrt(p + qq));
                Rb[L & 1][0] = sBp[L] * sF[L];
            }
            __syncthreads();
            // R^n_tuv on the simplex t+u+v <= L-n, level n kept in Rb[n & 1]
            // (same recursion and direction order as buildR0Simplex: z, then y, then x).
            for (int n = L - 1; n >= 0; --n) {
                const double* c = Rb[(n + 1) & 1];
                double* x = Rb[n & 1];
                const int top = L - n;
                for (int idx = tid; idx < s3; idx += bs) {
                    const int t = idx / s2, u = (idx / s) % s, v = idx % s;
                    if (t + u + v > top) continue;
                    double r;
                    if (t + u + v == 0)
                        r = sBp[n] * sF[n];
                    else if (v > 0)
                        r = (v >= 2 ? (v - 1) * c[idx - 2] : 0.0) + sW[2] * c[idx - 1];
                    else if (u > 0)
                        r = (u >= 2 ? (u - 1) * c[idx - 2 * s] : 0.0) + sW[1] * c[idx - s];
                    else
                        r = (t >= 2 ? (t - 1) * c[idx - 2 * s2] : 0.0) + sW[0] * c[idx - s2];
                    x[idx] = r;
                }
                __syncthreads();
            }
            const double* R = Rb[0];
            const double pref = sPref;
            const int ia = d.pp_ia[pb], ibb = d.pp_ib[pb], ic = d.pp_ia[pk], idd = d.pp_ib[pk];
            const double* Eab = d.E + d.pp_E_off[pb];
            const double* Ecd = d.E + d.pp_E_off[pk];

            for (int k = tid; k < ncart; k += bs) {
                const int e = k % nD, c = (k / nD) % nC, b = (k / (nD * nC)) % nB, a = k / (nD * nC * nB);
                const double w = pref * d.coef[d.sh_coef_off[A] + a * npA + ia] * d.coef[d.sh_coef_off[B] + b * npB + ibb]
                               * d.coef[d.sh_coef_off[C] + c * npC + ic] * d.coef[d.sh_coef_off[D] + e * npD + idd];
                if (w == 0.0) continue;
                const int* la = d.lmn + d.sh_lmn_off[A] + 3 * a;
                const int* lb = d.lmn + d.sh_lmn_off[B] + 3 * b;
                const int* lc = d.lmn + d.sh_lmn_off[C] + 3 * c;
                const int* ld = d.lmn + d.sh_lmn_off[D] + 3 * e;
                const int tx = la[0] + lb[0], ty = la[1] + lb[1], tz = la[2] + lb[2];
                const int kx = lc[0] + ld[0], ky = lc[1] + ld[1], kz = lc[2] + ld[2];
                double ex[5], ey[5], ez[5], fx[5], fy[5], fz[5];
                for (int m = 0; m <= tx; ++m) ex[m] = Eab[m * sAB_t + la[0] * sAB_i + lb[0]];
                for (int m = 0; m <= ty; ++m) ey[m] = Eab[szAB + m * sAB_t + la[1] * sAB_i + lb[1]];
                for (int m = 0; m <= tz; ++m) ez[m] = Eab[2 * szAB + m * sAB_t + la[2] * sAB_i + lb[2]];
                for (int m = 0; m <= kx; ++m) fx[m] = ((m & 1) ? -1.0 : 1.0) * Ecd[m * sCD_t + lc[0] * sCD_i + ld[0]];
                for (int m = 0; m <= ky; ++m) fy[m] = ((m & 1) ? -1.0 : 1.0) * Ecd[szCD + m * sCD_t + lc[1] * sCD_i + ld[1]];
                for (int m = 0; m <= kz; ++m) fz[m] = ((m & 1) ? -1.0 : 1.0) * Ecd[2 * szCD + m * sCD_t + lc[2] * sCD_i + ld[2]];
                double val = 0.0;
                for (int t = 0; t <= tx; ++t) {
                    if (ex[t] == 0.0) continue;
                    for (int u = 0; u <= ty; ++u) {
                        const double etu = ex[t] * ey[u];
                        if (etu == 0.0) continue;
                        for (int v = 0; v <= tz; ++v) {
                            const double eb = etu * ez[v];
                            if (eb == 0.0) continue;
                            double kv = 0.0;
                            for (int tau = 0; tau <= kx; ++tau)
                                for (int nu = 0; nu <= ky; ++nu) {
                                    const double ftn = fx[tau] * fy[nu];
                                    const double* Rr = R + ((t + tau) * s + (u + nu)) * s + v;
                                    for (int ph = 0; ph <= kz; ++ph) kv += ftn * fz[ph] * Rr[ph];
                                }
                            val += eb * kv;
                        }
                    }
                }
                blk[k] += w * val;
            }
        }
    }
    __syncthreads();

    // cartesian -> active (spherical) basis, one axis at a time (transformAxis)
    int dim[4] = { nA, nB, nC, nD };
    const int s4[4] = { A, B, C, D };
    double* x = blk;
    double* y = tmp;
    for (int ax = 0; ax < 4; ++ax) {
        const int toff = d.sh_T_off[s4[ax]];
        if (toff < 0) continue;
        const int nin = dim[ax];
        // only d shells carry a T (6 -> 5, checked at creation); s/p are the identity
        const int nout = (nin == 6) ? 5 : nin;
        int outer = 1, inner = 1;
        for (int m = 0; m < ax; ++m) outer *= dim[m];
        for (int m = ax + 1; m < 4; ++m) inner *= dim[m];
        const int total = outer * nout * inner;
        for (int idx = tid; idx < total; idx += bs) {
            const int r = idx % inner, jj = (idx / inner) % nout, o = idx / (inner * nout);
            double acc = 0.0;
            for (int a = 0; a < nin; ++a) {
                const double tt = d.T[toff + a * nout + jj];
                if (tt == 0.0) continue;
                acc += tt * x[((size_t)o * nin + a) * inner + r];
            }
            y[idx] = acc;
        }
        __syncthreads();
        double* sw = x; x = y; y = sw;
        dim[ax] = nout;
    }

    // canonical representative of every 8-fold orbit inside this quartet -> set8
    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B], oC = d.sh_out_first[C], oD = d.sh_out_first[D];
    const int total = dim[0] * dim[1] * dim[2] * dim[3];
    const size_t n = (size_t)d.n_out;
    for (int idx = tid; idx < total; idx += bs) {
        const int e = idx % dim[3], c = (idx / dim[3]) % dim[2], b = (idx / (dim[3] * dim[2])) % dim[1],
                  a = idx / (dim[3] * dim[2] * dim[1]);
        if (A == B && a < b) continue;
        if (C == D && c < e) continue;
        if (i == j && a * dim[1] + b < c * dim[3] + e) continue;
        set8(out, n, oA + a, oB + b, oC + c, oD + e, x[idx]);
    }
}

// ---------------------------------------------------------------------------------------------
// Thread-per-quartet variant for the small classes (Claude Generated, Oct 2026).
// The block kernel above leaves most lanes idle when a quartet has only a handful of
// cartesian components (an (ss|ss) quartet: 1 of 32) and synchronises the block for every
// primitive quartet. Small classes are the bulk of all quartets, so here one THREAD owns one
// shell quartet: everything is compile-time sized (template on LA..LD), R^n_tuv is kept in a
// single simplex array updated in place from the highest t+u+v down (R^n_tuv only reads
// entries with a smaller t+u+v of level n+1), and the component sums are fully unrolled.
// Same arithmetic and the same primitive screening as the block kernel.
// ---------------------------------------------------------------------------------------------

__host__ __device__ constexpr int ncomp(int L) { return (L + 1) * (L + 2) / 2; }
// cartesian power on axis ax of component c of a shell with momentum L, in the order of the
// basis parser: p = x, y, z; d = xx, yy, zz, xy, xz, yz (checked on the host before use)
__host__ __device__ constexpr int cpow(int L, int c, int ax)
{
    return L == 0 ? 0
        : L == 1 ? (c == ax ? 1 : 0)
        : (c < 3 ? (c == ax ? 2 : 0) : c == 3 ? (ax < 2 ? 1 : 0) : c == 4 ? (ax != 1 ? 1 : 0) : (ax > 0 ? 1 : 0));
}
// index of (t,u,v) in a simplex t+u+v <= L, layout by t, then u, then v
__host__ __device__ constexpr int sidx(int L, int t, int u, int v)
{
    int k = 0;
    for (int a = 0; a < t; ++a) k += (L - a + 1) * (L - a + 2) / 2;
    for (int b = 0; b < u; ++b) k += (L - t - b + 1);
    return k + v;
}

// One axis of the 6d -> 5d transform on a register block with compile-time dims.
template <int D0, int D1, int D2, int D3, int AX>
__device__ inline void transformAxisReg(const double* in, double* out, const double* T)
{
    constexpr int dims[4] = { D0, D1, D2, D3 };
    constexpr int nin = dims[AX], nout = 5;
    constexpr int outer = (AX > 0 ? D0 : 1) * (AX > 1 ? D1 : 1) * (AX > 2 ? D2 : 1);
    constexpr int inner = (AX < 3 ? D3 : 1) * (AX < 2 ? D2 : 1) * (AX < 1 ? D1 : 1);
#pragma unroll
    for (int o = 0; o < outer; ++o)
#pragma unroll
        for (int jj = 0; jj < nout; ++jj)
#pragma unroll
            for (int r = 0; r < inner; ++r) {
                double acc = 0.0;
#pragma unroll
                for (int a = 0; a < nin; ++a) {
                    const double tt = T[a * nout + jj];
                    if (tt != 0.0) acc += tt * in[(o * nin + a) * inner + r];
                }
                out[(o * nout + jj) * inner + r] = acc;
            }
}

template <int LA, int LB, int LC, int LD, bool SPH>
__global__ void __launch_bounds__(64) eriThreadKernel(DevTables d, long long q0, long long nq,
                                                       double* __restrict__ out)
{
    const long long qi = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (qi >= nq) return;
    constexpr int nA = ncomp(LA), nB = ncomp(LB), nC = ncomp(LC), nD = ncomp(LD);
    constexpr int L = LA + LB + LC + LD;
    constexpr int NC = nA * nB * nC * nD;
    constexpr int NS = (L + 1) * (L + 2) * (L + 3) / 6;  // simplex size
    constexpr int sAB_t = (LA + 1) * (LB + 1), sAB_i = LB + 1, szAB = (LA + LB + 1) * sAB_t;
    constexpr int sCD_t = (LC + 1) * (LD + 1), sCD_i = LD + 1, szCD = (LC + LD + 1) * sCD_t;

    const int i = d.q_bra[q0 + qi], j = d.q_ket[q0 + qi];
    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const int pb0 = d.pr_prim_off[i], nbp = d.pr_nprim[i];
    const int pk0 = d.pr_prim_off[j], nkp = d.pr_nprim[j];
    const double ket_max = nkp > 0 ? d.pp_schwarz[pk0] : 0.0;
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B], npC = d.sh_nprim[C], npD = d.sh_nprim[D];
    const double* cA = d.coef + d.sh_coef_off[A];
    const double* cB = d.coef + d.sh_coef_off[B];
    const double* cC = d.coef + d.sh_coef_off[C];
    const double* cD = d.coef + d.sh_coef_off[D];

    double acc[NC];
#pragma unroll
    for (int k = 0; k < NC; ++k) acc[k] = 0.0;

    for (int ib = 0; ib < nbp; ++ib) {
        const int pb = pb0 + ib;
        const double sb = d.pp_schwarz[pb];
        if (sb * ket_max < d.prim_thr) break;
        const double p = d.pp_p[pb];
        const double Px = d.pp_P[3 * pb], Py = d.pp_P[3 * pb + 1], Pz = d.pp_P[3 * pb + 2];
        const double* Eab = d.E + d.pp_E_off[pb];
        const int ia = d.pp_ia[pb], ibb = d.pp_ib[pb];
        for (int ik = 0; ik < nkp; ++ik) {
            const int pk = pk0 + ik;
            if (sb * d.pp_schwarz[pk] < d.prim_thr) break;
            const double qq = d.pp_p[pk];
            const double rho = p * qq / (p + qq);
            const double Wx = Px - d.pp_P[3 * pk], Wy = Py - d.pp_P[3 * pk + 1], Wz = Pz - d.pp_P[3 * pk + 2];
            double F[L + 1];
            boysDev(L, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F, d);
            double bp[L + 1];
            bp[0] = 1.0;
#pragma unroll
            for (int n = 1; n <= L; ++n) bp[n] = bp[n - 1] * (-2.0 * rho);
            // R^n on the simplex, in place: level n overwrites entries with t+u+v <= L-n,
            // largest t+u+v first (reads only smaller sums, still level n+1).
            double R[NS];
            R[0] = bp[L] * F[L];
#pragma unroll
            for (int n = L - 1; n >= 0; --n) {
#pragma unroll
                for (int sum = L - n; sum >= 1; --sum)
#pragma unroll
                    for (int t = 0; t <= sum; ++t)
#pragma unroll
                        for (int u = 0; u <= sum - t; ++u) {
                            const int v = sum - t - u;
                            double r;
                            if (v > 0)
                                r = (v >= 2 ? (v - 1) * R[sidx(L, t, u, v - 2)] : 0.0) + Wz * R[sidx(L, t, u, v - 1)];
                            else if (u > 0)
                                r = (u >= 2 ? (u - 1) * R[sidx(L, t, u - 2, v)] : 0.0) + Wy * R[sidx(L, t, u - 1, v)];
                            else
                                r = (t >= 2 ? (t - 1) * R[sidx(L, t - 2, u, v)] : 0.0) + Wx * R[sidx(L, t - 1, u, v)];
                            R[sidx(L, t, u, v)] = r;
                        }
                R[0] = bp[n] * F[n];
            }
            const double pref = d.pref / (p * qq * sqrt(p + qq));
            const double* Ecd = d.E + d.pp_E_off[pk];
            const int ic = d.pp_ia[pk], idd = d.pp_ib[pk];

#pragma unroll
            for (int a = 0; a < nA; ++a)
#pragma unroll
                for (int b = 0; b < nB; ++b) {
                    const double wab = pref * cA[a * npA + ia] * cB[b * npB + ibb];
#pragma unroll
                    for (int c = 0; c < nC; ++c)
#pragma unroll
                        for (int e = 0; e < nD; ++e) {
                            const double w = wab * cC[c * npC + ic] * cD[e * npD + idd];
                            const int la0 = cpow(LA, a, 0), la1 = cpow(LA, a, 1), la2 = cpow(LA, a, 2);
                            const int lb0 = cpow(LB, b, 0), lb1 = cpow(LB, b, 1), lb2 = cpow(LB, b, 2);
                            const int lc0 = cpow(LC, c, 0), lc1 = cpow(LC, c, 1), lc2 = cpow(LC, c, 2);
                            const int ld0 = cpow(LD, e, 0), ld1 = cpow(LD, e, 1), ld2 = cpow(LD, e, 2);
                            const int tx = la0 + lb0, ty = la1 + lb1, tz = la2 + lb2;
                            const int kx = lc0 + ld0, ky = lc1 + ld1, kz = lc2 + ld2;
                            double val = 0.0;
#pragma unroll
                            for (int t = 0; t <= tx; ++t)
#pragma unroll
                                for (int u = 0; u <= ty; ++u)
#pragma unroll
                                    for (int v = 0; v <= tz; ++v) {
                                        const double eb = Eab[t * sAB_t + la0 * sAB_i + lb0]
                                                        * Eab[szAB + u * sAB_t + la1 * sAB_i + lb1]
                                                        * Eab[2 * szAB + v * sAB_t + la2 * sAB_i + lb2];
                                        double kv = 0.0;
#pragma unroll
                                        for (int tau = 0; tau <= kx; ++tau)
#pragma unroll
                                            for (int nu = 0; nu <= ky; ++nu)
#pragma unroll
                                                for (int ph = 0; ph <= kz; ++ph) {
                                                    const double f = (((tau + nu + ph) & 1) ? -1.0 : 1.0)
                                                        * Ecd[tau * sCD_t + lc0 * sCD_i + ld0]
                                                        * Ecd[szCD + nu * sCD_t + lc1 * sCD_i + ld1]
                                                        * Ecd[2 * szCD + ph * sCD_t + lc2 * sCD_i + ld2];
                                                    kv += f * R[sidx(L, t + tau, u + nu, v + ph)];
                                                }
                                        val += eb * kv;
                                    }
                            acc[((a * nB + b) * nC + c) * nD + e] += w * val;
                        }
                }
        }
    }

    // 6d -> 5d per axis (register blocks; SPH = spherical basis), then the canonical writes
    constexpr bool tA2 = SPH && LA == 2, tB2 = SPH && LB == 2, tC2 = SPH && LC == 2, tD2 = SPH && LD == 2;
    constexpr int oA_ = tA2 ? 5 : nA, oB_ = tB2 ? 5 : nB, oC_ = tC2 ? 5 : nC, oD_ = tD2 ? 5 : nD;
    double s1[oA_ * nB * nC * nD], s2[oA_ * oB_ * nC * nD], s3[oA_ * oB_ * oC_ * nD], s4[oA_ * oB_ * oC_ * oD_];
    const int tA = d.sh_T_off[A], tB = d.sh_T_off[B], tC = d.sh_T_off[C], tD = d.sh_T_off[D];
    if constexpr (tA2) { transformAxisReg<nA, nB, nC, nD, 0>(acc, s1, d.T + tA); }
    else {
#pragma unroll
        for (int k = 0; k < NC; ++k) s1[k] = acc[k];
    }
    if constexpr (tB2) { transformAxisReg<oA_, nB, nC, nD, 1>(s1, s2, d.T + tB); }
    else {
#pragma unroll
        for (int k = 0; k < oA_ * nB * nC * nD; ++k) s2[k] = s1[k];
    }
    if constexpr (tC2) { transformAxisReg<oA_, oB_, nC, nD, 2>(s2, s3, d.T + tC); }
    else {
#pragma unroll
        for (int k = 0; k < oA_ * oB_ * nC * nD; ++k) s3[k] = s2[k];
    }
    if constexpr (tD2) { transformAxisReg<oA_, oB_, oC_, nD, 3>(s3, s4, d.T + tD); }
    else {
#pragma unroll
        for (int k = 0; k < oA_ * oB_ * oC_ * oD_; ++k) s4[k] = s3[k];
    }

    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B], oC = d.sh_out_first[C], oD = d.sh_out_first[D];
    const size_t n = (size_t)d.n_out;
#pragma unroll
    for (int a = 0; a < oA_; ++a)
#pragma unroll
        for (int b = 0; b < oB_; ++b)
#pragma unroll
            for (int c = 0; c < oC_; ++c)
#pragma unroll
                for (int e = 0; e < oD_; ++e) {
                    if (A == B && a < b) continue;
                    if (C == D && c < e) continue;
                    if (i == j && a * oB_ + b < c * oD_ + e) continue;
                    set8(out, n, oA + a, oB + b, oC + c, oD + e, s4[((a * oB_ + b) * oC_ + c) * oD_ + e]);
                }
}

// Largest class (product of cartesian component counts) run thread-per-quartet; the rest
// use the block kernel. Chosen by measurement on the RTX 5090 (docs/QM_GPU_ROADMAP.md).
constexpr int kThreadKernelMaxNC = 27;

template <int LA, int LB, int LC, int LD>
bool launchThreadKernel(const DevTables& d, long long q0, long long nq, double* out, cudaStream_t st, bool sph)
{
    if constexpr (ncomp(LA) * ncomp(LB) * ncomp(LC) * ncomp(LD) <= kThreadKernelMaxNC) {
        constexpr int bs = 64;
        const long long nb = (nq + bs - 1) / bs;
        if (sph)
            eriThreadKernel<LA, LB, LC, LD, true><<<(unsigned)nb, bs, 0, st>>>(d, q0, nq, out);
        else
            eriThreadKernel<LA, LB, LC, LD, false><<<(unsigned)nb, bs, 0, st>>>(d, q0, nq, out);
        return true;
    } else {
        return false;
    }
}

// Dispatch the 81 (L <= 2) classes to their template instance; false = use the block kernel.
bool launchThreadClass(int LA, int LB, int LC, int LD, const DevTables& d, long long q0, long long nq,
                       double* out, cudaStream_t st, bool sph)
{
#define QG2_CASE(a, b, c, e) \
    case ((a * 3 + b) * 3 + c) * 3 + e: return launchThreadKernel<a, b, c, e>(d, q0, nq, out, st, sph);
#define QG2_C3(a, b, c) QG2_CASE(a, b, c, 0) QG2_CASE(a, b, c, 1) QG2_CASE(a, b, c, 2)
#define QG2_C2(a, b) QG2_C3(a, b, 0) QG2_C3(a, b, 1) QG2_C3(a, b, 2)
#define QG2_C1(a) QG2_C2(a, 0) QG2_C2(a, 1) QG2_C2(a, 2)
    switch (((LA * 3 + LB) * 3 + LC) * 3 + LD) {
        QG2_C1(0)
        QG2_C1(1)
        QG2_C1(2)
    default: return false;
    }
#undef QG2_C1
#undef QG2_C2
#undef QG2_C3
#undef QG2_CASE
}

template <typename T>
bool upload(T** dptr, const T* host, size_t count, std::vector<void*>& owned)
{
    *dptr = nullptr;
    if (count == 0) return true;
    if (cudaMalloc(dptr, count * sizeof(T)) != cudaSuccess) { *dptr = nullptr; return false; }
    owned.push_back(*dptr);
    return cudaMemcpy(*dptr, host, count * sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess;
}

void setWhy(char* why, int len, const std::string& s)
{
    if (why && len > 0) {
        std::strncpy(why, s.c_str(), (size_t)len - 1);
        why[len - 1] = '\0';
    }
}

int ncompOf(int L) { return (L + 1) * (L + 2) / 2; }

} // namespace

// Build the ERI tensor on the device and return a J/K builder that owns it (QMGpuJK*), or
// nullptr with a reason in `why` (unsupported basis, out of device memory, CUDA error).
extern "C" void* curcuma_cuda_create_qm_eri_jk(const qmint::QMGpuEriTables* t, int device, char* why, int whylen)
{
    using clk = std::chrono::steady_clock;
    auto ms = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const auto t0 = clk::now();
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return nullptr; }
    cudaFree(nullptr);  // context creation, timed separately
    {
        // lazy module loading (CUDA_MODULE_LOADING=LAZY): load both kernels now, so the
        // kernel timing below is kernel work only
        cudaFuncAttributes a{};
        cudaFuncGetAttributes(&a, eriClassKernel);
        qmgpu::warmJkKernel();
    }
    const auto t1 = clk::now();

    if (t->boys_ntaylor != kNTaylor) { setWhy(why, whylen, "Boys table layout mismatch"); return nullptr; }
    for (int s = 0; s < t->nshell; ++s)
        if (t->sh_L[s] > 2 || t->sh_ncomp[s] != ncompOf(t->sh_L[s])) {
            setWhy(why, whylen, "shell with L > 2 (device kernel handles s, p, d)");
            return nullptr;
        }
    for (int s = 0; s < t->nshell; ++s)
        if (t->sh_T_off[s] >= 0 && !(t->sh_L[s] == 2 && t->sh_out_width[s] == 5)) {
            setWhy(why, whylen, "unexpected spherical transform block");
            return nullptr;
        }

    const size_t n = (size_t)t->n_out;
    const size_t tensor_bytes = n * n * n * n * sizeof(double);
    const size_t table_bytes = (size_t)t->n_E * 8 + (size_t)t->nquartet * 8 + (size_t)t->nprimpair * 64 + (1u << 20);
    size_t free_b = 0, total_b = 0;
    cudaMemGetInfo(&free_b, &total_b);
    const size_t margin = 256ull << 20;
    if (tensor_bytes + table_bytes + 3 * n * n * 8 + margin > free_b) {
        char buf[256];
        std::snprintf(buf, sizeof buf, "ERI tensor %.1f MB does not fit into free device memory %.1f MB",
                      tensor_bytes / 1.0e6, free_b / 1.0e6);
        setWhy(why, whylen, buf);
        return nullptr;
    }

    std::vector<void*> owned;
    auto cleanup = [&owned]() { for (void* p : owned) cudaFree(p); owned.clear(); };
    DevTables d{};
    d.n_out = t->n_out;
    d.prim_thr = t->prim_thr;
    d.pref = t->pref;
    d.boys_ntab = t->boys_ntab;
    d.boys_dT = t->boys_dT;
    d.boys_Tmax = t->boys_Tmax;
    bool ok = upload(const_cast<int**>(&d.sh_ncomp), t->sh_ncomp, t->nshell, owned)
        && upload(const_cast<int**>(&d.sh_nprim), t->sh_nprim, t->nshell, owned)
        && upload(const_cast<int**>(&d.sh_lmn_off), t->sh_lmn_off, t->nshell, owned)
        && upload(const_cast<int**>(&d.sh_coef_off), t->sh_coef_off, t->nshell, owned)
        && upload(const_cast<int**>(&d.sh_out_first), t->sh_out_first, t->nshell, owned)
        && upload(const_cast<int**>(&d.sh_T_off), t->sh_T_off, t->nshell, owned)
        && upload(const_cast<int**>(&d.lmn), t->lmn, t->n_lmn, owned)
        && upload(const_cast<double**>(&d.coef), t->coef, t->n_coef, owned)
        && upload(const_cast<double**>(&d.T), t->T, t->n_T, owned)
        && upload(const_cast<int**>(&d.pr_A), t->pr_A, t->npair, owned)
        && upload(const_cast<int**>(&d.pr_B), t->pr_B, t->npair, owned)
        && upload(const_cast<int**>(&d.pr_prim_off), t->pr_prim_off, t->npair, owned)
        && upload(const_cast<int**>(&d.pr_nprim), t->pr_nprim, t->npair, owned)
        && upload(const_cast<int**>(&d.pp_ia), t->pp_ia, t->nprimpair, owned)
        && upload(const_cast<int**>(&d.pp_ib), t->pp_ib, t->nprimpair, owned)
        && upload(const_cast<double**>(&d.pp_p), t->pp_p, t->nprimpair, owned)
        && upload(const_cast<double**>(&d.pp_P), t->pp_P, 3 * (size_t)t->nprimpair, owned)
        && upload(const_cast<double**>(&d.pp_schwarz), t->pp_schwarz, t->nprimpair, owned)
        && upload(const_cast<long long**>(&d.pp_E_off), t->pp_E_off, t->nprimpair, owned)
        && upload(const_cast<double**>(&d.E), t->E, (size_t)t->n_E, owned)
        && upload(const_cast<int**>(&d.q_bra), t->q_bra, (size_t)t->nquartet, owned)
        && upload(const_cast<int**>(&d.q_ket), t->q_ket, (size_t)t->nquartet, owned)
        && upload(const_cast<double**>(&d.boys_f), t->boys_f, (size_t)t->boys_npts * (t->boys_ntab + 1), owned);
    if (!ok) { cleanup(); setWhy(why, whylen, "upload of the ERI tables failed"); return nullptr; }

    double* dE = nullptr;
    if (cudaMalloc(&dE, tensor_bytes) != cudaSuccess) {
        cleanup();
        setWhy(why, whylen, "cudaMalloc of the ERI tensor failed");
        return nullptr;
    }
    cudaMemset(dE, 0, tensor_bytes);
    cudaDeviceSynchronize();
    const auto t2 = clk::now();

    cudaEvent_t e0, e1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    // CURCUMA_QM_GPU_ERI_PROFILE=1: per-class kernel times on stderr (serialises the launches)
    const char* penv = std::getenv("CURCUMA_QM_GPU_ERI_PROFILE");
    const bool prof = penv && penv[0] == '1';
    cudaEvent_t p0, p1;
    cudaEventCreate(&p0);
    cudaEventCreate(&p1);
    // the thread kernel hard-codes the parser's component order; verify it, and whether the
    // d shells are spherical (all carry a T) or cartesian (none does)
    bool order_ok = true, spherical = false, cartesian_d = false;
    for (int sh = 0; sh < t->nshell; ++sh) {
        const int Ls = t->sh_L[sh];
        for (int c = 0; c < t->sh_ncomp[sh]; ++c)
            for (int ax = 0; ax < 3; ++ax)
                if (t->lmn[t->sh_lmn_off[sh] + 3 * c + ax] != cpow(Ls, c, ax)) order_ok = false;
        if (Ls == 2) (t->sh_T_off[sh] >= 0 ? spherical : cartesian_d) = true;
    }
    if (spherical && cartesian_d) order_ok = false;
    cudaEventRecord(e0);
    cudaEventRecord(p0);
    // Classes are independent (disjoint tensor elements), so they run on several streams:
    // for small molecules the per-class latency, not the work, sets the time.
    constexpr int kStreams = 8;
    cudaStream_t streams[kStreams];
    for (int k = 0; k < kStreams; ++k) cudaStreamCreateWithFlags(&streams[k], cudaStreamNonBlocking);
    const char* tenv = std::getenv("CURCUMA_QM_GPU_ERI_BLOCK_ONLY");  // diagnostic: block kernel everywhere
    const bool thread_ok = order_ok && !(tenv && tenv[0] == '1');
    int nthread_classes = 0;
    for (int c = 0; c < t->nclass; ++c) {
        const long long q0 = t->class_off[c], nq = t->class_off[c + 1] - q0;
        if (nq <= 0) continue;
        const int LA = t->class_L[4 * c], LB = t->class_L[4 * c + 1], LC = t->class_L[4 * c + 2], LD = t->class_L[4 * c + 3];
        const int L = LA + LB + LC + LD, s3 = (L + 1) * (L + 1) * (L + 1);
        const int ncart = ncompOf(LA) * ncompOf(LB) * ncompOf(LC) * ncompOf(LD);
        cudaStream_t st = prof ? streams[0] : streams[c % kStreams];
        int bs = 32;
        bool threaded = thread_ok && launchThreadClass(LA, LB, LC, LD, d, q0, nq, dE, st, spherical);
        if (threaded) {
            ++nthread_classes;
            bs = 1;
        } else {
            const size_t shm = (size_t)(2 * s3 + 2 * ncart) * sizeof(double);
            while (bs < ncart && bs < 128) bs *= 2;
            // grid.x limit 2^31-1: chunk the class if ever needed
            for (long long off = 0; off < nq; off += 2000000000LL) {
                const long long cnt = (nq - off < 2000000000LL) ? nq - off : 2000000000LL;
                eriClassKernel<<<(unsigned)cnt, bs, shm, st>>>(d, q0 + off, LA, LB, LC, LD, dE);
            }
        }
        if (cudaGetLastError() != cudaSuccess) break;
        if (prof) {
            cudaEventRecord(p1, st);
            cudaEventSynchronize(p1);
            float pm = 0.0f;
            cudaEventElapsedTime(&pm, p0, p1);
            std::fprintf(stderr, "qm_gpu_eri class %d%d%d%d quartets %lld ncart %d %s: %.3f ms\n", LA, LB, LC, LD,
                         nq, ncart, threaded ? "thread" : "block", (double)pm);
            cudaEventRecord(p0, st);
        }
    }
    for (int k = 0; k < kStreams; ++k) {
        cudaStreamSynchronize(streams[k]);
        cudaStreamDestroy(streams[k]);
    }
    cudaEventRecord(e1);
    const cudaError_t kerr = cudaEventSynchronize(e1) == cudaSuccess ? cudaGetLastError() : cudaErrorUnknown;
    float kms = 0.0f;
    cudaEventElapsedTime(&kms, e0, e1);
    cudaEventDestroy(e0);
    cudaEventDestroy(e1);
    cudaEventDestroy(p0);
    cudaEventDestroy(p1);
    cleanup();
    if (kerr != cudaSuccess) {
        cudaFree(dE);
        setWhy(why, whylen, std::string("ERI kernel failed: ") + cudaGetErrorString(kerr));
        return nullptr;
    }
    const auto t3 = clk::now();

    auto* jk = new qmgpu::QmGpuJkCuda(dE, (int)n, device, true);
    if (!jk->ok()) {
        delete jk;
        setWhy(why, whylen, "allocation of the J/K work buffers failed");
        return nullptr;
    }
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "context+module load %.1f ms, tables+tensor alloc %.1f ms, ERI kernels %.1f ms (%lld quartets, %d classes, %d thread-per-quartet), total %.1f ms",
                  ms(t0, t1), ms(t1, t2), (double)kms, t->nquartet, t->nclass, nthread_classes, ms(t0, t3));
    jk->setStats(buf);
    return static_cast<qmint::QMGpuJK*>(jk);
}
