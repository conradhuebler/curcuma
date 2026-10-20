/*
 * < CUDA integral-direct J/K: thread-per-quartet McMurchie-Davidson kernels (templates) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 20 2026, docs/QM_GPU_ROADMAP.md section 2 "Large systems").
//
// Kernel templates shared by qm_gpu_direct_md.cu (FP64 instances, the exact builds) and
// qm_gpu_direct_fp32.cu (FP32 instances, early SCF phase of -qm.gpu_eri_precision mixed).
//
// Why: the block kernels (eriClassKernel, eriHermiteKernel) ran one thread block per shell quartet
// and synchronised the block several times per primitive quartet (Boys by one thread, R level by
// level, ket sums, bra sums); on the RTX 5090 the classes with more than 27 cartesian components
// reached only 2-24 quartets/us and took ~80 % of the J/K time. Here, as in the thread kernel of
// qm_gpu_eri_kernels.cuh, one THREAD owns one shell quartet (or a fixed slice of its ket components,
// see MODE below), everything is sized at compile time, and no barrier is needed.
//
// Same McMurchie-Davidson integrals (Helgaker, Jorgensen, Olsen ch. 9.9), reorganised so that the
// work per primitive quartet is small:
//   1. setup (mdInit): per primitive pair k the Hermite expansion of every ACTIVE component pair,
//        F[k][a'b'][h] = sum_{a,b} T_A(a,a') T_B(b,b') c_a c_b E^x_t E^y_u E^z_v,
//      contraction coefficients and the 6d->5d transform folded in (a linear map, so it may act on the
//      Hermite coefficients instead of the finished integrals);
//   2. per bra primitive pair pb, contracted over the ket primitive pairs pk:
//        Kc[h_ab][c'e'] = sum_pk pref(pb,pk) sum_{h_cd} (-1)^{|h_cd|} F[pk][c'e'][h_cd] R_{h_ab + h_cd}(pb,pk)
//      (R^0_tuv from the Boys function as everywhere else; pref = 2 pi^(5/2) / (p q sqrt(p+q)));
//   3. (a'b'|c'e') += sum_{h_ab} F[pb][a'b'][h_ab] Kc[h_ab][c'e'];
//   4. digest into Jt, Kt with the weights of digestBlock (FP64 atomicAdd).
// Entries of F that vanish for every primitive pair (Hermite index beyond the cartesian powers of
// the components, known at compile time from the component order and the fixed pattern of the 5d
// combinations) are skipped at compile time. All loops over components and Hermite indices are
// compile-time loops (sfor), so the front end, not the optimiser, removes the zero terms (an
// earlier #pragma unroll form needed ~25 min of nvcc for 8 classes).
//
// MODE (chosen at compile time from the register state R + Kc + acc, mdMode): 0 = one thread per
// quartet; 1 = one thread per (quartet, c'); 2 = per (quartet, e'); 3 = per (quartet, c', e'). The
// slice index is blockIdx.y, so it is uniform per block and its components are template constants
// (no divergence between lanes of different slices). Partial sums of a slice go straight to the
// atomics. Screening (quartetOf, primitive break) is the one of the other kernels, so the quartet
// count is unchanged; the FP64 integrals differ from eriHermiteKernel only in the summation order.
//
// RT = float (FP32 phase): R, Boys, Kc, acc and the digest sums in FP32 from FP32 copies of F, P and
// the Boys table; the product-centre difference W = P - Q is formed in FP64 before the cast (large
// coordinates would lose the short distances), the atomics into Jt, Kt stay FP64.

#include "qm_gpu_direct_md.cuh"

#include <cuda_runtime.h>

#include <type_traits>
#include <utility>

namespace {

using qmgpu::MdTables;

__host__ __device__ constexpr int nOut(int L) { return L == 2 ? 5 : ncomp(L); }
__host__ __device__ constexpr int nSim(int L) { return (L + 1) * (L + 2) * (L + 3) / 6; }

/// (t, u, v) of simplex index h (layout of sidx: t, then u, then v).
struct Tuv { int t, u, v; };
__host__ __device__ constexpr Tuv tuvOf(int L, int h)
{
    for (int t = 0; t <= L; ++t)
        for (int u = 0; u <= L - t; ++u)
            for (int v = 0; v <= L - t - u; ++v) {
                if (h == 0) return Tuv{ t, u, v };
                --h;
            }
    return Tuv{ 0, 0, 0 };
}

/// Is cartesian component k (parser order) part of active component m of a shell L?
/// s, p: identity; d (always 5d here): z2 {xx,yy,zz}, xz, yz, x2-y2 {xx,yy}, xy -- the pattern of
/// buildSphericalTransform (qm_integrals.cpp), verified against the uploaded T in mdInit.
__host__ __device__ constexpr bool inSupp(int L, int m, int k)
{
    return L < 2 ? (m == k) : (m == 0 ? k < 3 : m == 1 ? k == 4 : m == 2 ? k == 5 : m == 3 ? k < 2 : k == 3);
}

/// Can F[k][m1 m2][h=(t,u,v)] be nonzero? E^x_t(lx1, lx2) = 0 for t > lx1 + lx2.
__host__ __device__ constexpr bool hnz(int L1, int m1, int L2, int m2, int t, int u, int v)
{
    for (int a = 0; a < ncomp(L1); ++a) {
        if (!inSupp(L1, m1, a)) continue;
        for (int b = 0; b < ncomp(L2); ++b) {
            if (!inSupp(L2, m2, b)) continue;
            if (t <= cpow(L1, a, 0) + cpow(L2, b, 0) && u <= cpow(L1, a, 1) + cpow(L2, b, 1)
                && v <= cpow(L1, a, 2) + cpow(L2, b, 2))
                return true;
        }
    }
    return false;
}

/// Register state (values) of one thread for a mode; see the file header.
__host__ __device__ constexpr int mdState(int LA, int LB, int LC, int LD, int mode)
{
    const int sc = (mode & 1) ? 1 : nOut(LC), sd = (mode & 2) ? 1 : nOut(LD);
    return nSim(LA + LB + LC + LD) + (nSim(LA + LB) + nOut(LA) * nOut(LB)) * sc * sd;
}

/// Smallest split whose state fits (fewer lanes first); thresholds by measurement (RTX 5090).
__host__ __device__ constexpr int mdMode(int LA, int LB, int LC, int LD)
{
    if (mdState(LA, LB, LC, LD, 0) <= 100) return 0;
    const int nC = nOut(LC), nD = nOut(LD);
    const bool c_ok = nC > 1 && mdState(LA, LB, LC, LD, 1) <= 160;
    const bool d_ok = nD > 1 && mdState(LA, LB, LC, LD, 2) <= 160;
    if (c_ok && d_ok) return nC <= nD ? 1 : 2;
    if (c_ok) return 1;
    if (d_ok) return 2;
    return 3;
}

__host__ __device__ constexpr int mdSplit(int LA, int LB, int LC, int LD)
{
    const int mode = mdMode(LA, LB, LC, LD);
    return ((mode & 1) ? nOut(LC) : 1) * ((mode & 2) ? nOut(LD) : 1);
}

/// Compile-time loop: f(std::integral_constant<int, I>) for I = 0..N-1.
template <typename F, int... I>
__device__ __forceinline__ void sforImpl(F&& f, std::integer_sequence<int, I...>)
{
    (f(std::integral_constant<int, I>{}), ...);
}
template <int N, typename F>
__device__ __forceinline__ void sfor(F&& f)
{
    sforImpl(f, std::make_integer_sequence<int, N>{});
}

/// FP32 Boys function F_0..F_MAXN: same grid, Taylor table and recursions as boysDev, in FP32.
/// The Taylor offset i*dT - T uses dT = dT_hi + dT_lo and one fma (0.05 is not a float: the plain
/// i * 0.05f - T is off by up to ~2e-6 at the end of the grid). Claude Generated (Oct 20 2026).
template <int MAXN>
__device__ __forceinline__ void boysF32(float T, float* F, const float* __restrict__ tab, int ntab, float dT_hi,
                                        float dT_lo, float Tmax)
{
    if (T >= Tmax) {
        const float eT = __expf(-T);
        F[0] = 0.5f * sqrtf(3.14159265358979f / T);
#pragma unroll
        for (int n = 0; n < MAXN; ++n) F[n + 1] = ((2.0f * n + 1.0f) * F[n] - eT) / (2.0f * T);
        return;
    }
    const int i = (int)(T / dT_hi + 0.5f);
    const float dd = fmaf((float)i, dT_hi, -T) + (float)i * dT_lo;
    const float* fi = tab + (size_t)i * (ntab + 1) + MAXN;
    const float inv[kNTaylor] = { 1.0f, 1.0f, 1.0f / 2, 1.0f / 3, 1.0f / 4, 1.0f / 5, 1.0f / 6 };
    float acc = fi[kNTaylor - 1];
#pragma unroll
    for (int k = kNTaylor - 2; k >= 0; --k) acc = fi[k] + dd * inv[k + 1] * acc;
    F[MAXN] = acc;
    if (MAXN > 0) {
        const float eT = __expf(-T), T2 = 2.0f * T;
#pragma unroll
        for (int n = MAXN - 1; n >= 0; --n) F[n] = (T2 * F[n + 1] + eT) / (2.0f * n + 1.0f);
    }
}

/// One thread, one quartet (i, j), ket slice PART. Claude Generated (Oct 20 2026).
template <typename RT, int LA, int LB, int LC, int LD, int PART>
__device__ __forceinline__ void mdBody(const DevTables& d, const MdTables& m, int i, int j, double pthr)
{
    constexpr bool F32 = std::is_same<RT, float>::value;
    constexpr int MODE = mdMode(LA, LB, LC, LD);
    constexpr int nA = nOut(LA), nB = nOut(LB), nC = nOut(LC), nD = nOut(LD);
    constexpr int SC = (MODE & 1) ? 1 : nC, SD = (MODE & 2) ? 1 : nD;  // ket components per thread
    constexpr int nDgrp = nD / SD;
    constexpr int c0 = (MODE & 1) ? PART / nDgrp : 0;
    constexpr int e0 = (MODE & 2) ? PART % nDgrp : 0;
    constexpr int L = LA + LB + LC + LD, Lab = LA + LB, Lcd = LC + LD;
    constexpr int NS = nSim(L), NHAB = nSim(Lab), NHCD = nSim(Lcd);
    constexpr int NAB = nA * nB, NG = SC * SD;

    const RT* __restrict__ Ftab = nullptr;
    if constexpr (F32) Ftab = m.F32;
    else Ftab = m.F;

    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const int pb0 = d.pr_prim_off[i], nbp = d.pr_nprim[i];
    const int pk0 = d.pr_prim_off[j], nkp = d.pr_nprim[j];
    const double ket_max = nkp > 0 ? d.pp_schwarz[pk0] : 0.0;

    RT acc[NAB * NG];
#pragma unroll
    for (int k = 0; k < NAB * NG; ++k) acc[k] = RT(0);

    for (int ib = 0; ib < nbp; ++ib) {
        const int pb = pb0 + ib;
        const double sb = d.pp_schwarz[pb];
        if (sb * ket_max < pthr) break;  // same primitive screen as the other kernels
        const double p = d.pp_p[pb];
        const double Px = d.pp_P[3 * pb], Py = d.pp_P[3 * pb + 1], Pz = d.pp_P[3 * pb + 2];
        RT kc[NHAB * NG];
#pragma unroll
        for (int k = 0; k < NHAB * NG; ++k) kc[k] = RT(0);
        for (int ik = 0; ik < nkp; ++ik) {
            const int pk = pk0 + ik;
            if (sb * d.pp_schwarz[pk] < pthr) break;
            const double qq = d.pp_p[pk];
            RT Wx, Wy, Wz, rho, pref;
            RT Fb[L + 1];
            if constexpr (F32) {
                Wx = (float)(Px - d.pp_P[3 * pk]);
                Wy = (float)(Py - d.pp_P[3 * pk + 1]);
                Wz = (float)(Pz - d.pp_P[3 * pk + 2]);
                const float pf = (float)p, qf = (float)qq;
                rho = pf * qf / (pf + qf);
                boysF32<L>(rho * (Wx * Wx + Wy * Wy + Wz * Wz), Fb, m.boys_f32, d.boys_ntab, (float)d.boys_dT,
                           (float)(d.boys_dT - (double)(float)d.boys_dT), (float)d.boys_Tmax);
                pref = (float)d.pref / (pf * qf * sqrtf(pf + qf));
            } else {
                rho = p * qq / (p + qq);
                Wx = Px - d.pp_P[3 * pk];
                Wy = Py - d.pp_P[3 * pk + 1];
                Wz = Pz - d.pp_P[3 * pk + 2];
                boysDev(L, rho * (Wx * Wx + Wy * Wy + Wz * Wz), Fb, d);
                pref = d.pref / (p * qq * sqrt(p + qq));
            }
            RT bp[L + 1];
            bp[0] = RT(1);
#pragma unroll
            for (int n = 1; n <= L; ++n) bp[n] = bp[n - 1] * (RT(-2) * rho);
            // R^n on the simplex, in place: level n overwrites the entries with 1 <= t+u+v <= L-n in
            // descending simplex index; every read (one coordinate lowered by 1 or 2) has a smaller
            // index and is still level n+1 (same recursion and direction order as buildR0Simplex).
            RT R[NS];
            R[0] = bp[L] * Fb[L];
            sfor<L>([&](auto NI) {
                constexpr int n = L - 1 - decltype(NI)::value;
                sfor<NS - 1>([&](auto K) {
                    constexpr int h = NS - 1 - decltype(K)::value;
                    constexpr Tuv x = tuvOf(L, h);
                    if constexpr (x.t + x.u + x.v <= L - (L - 1 - decltype(NI)::value)) {
                        if constexpr (x.v > 0) {
                            RT r = Wz * R[sidx(L, x.t, x.u, x.v - 1)];
                            if constexpr (x.v >= 2) r += RT(x.v - 1) * R[sidx(L, x.t, x.u, x.v - 2)];
                            R[h] = r;
                        } else if constexpr (x.u > 0) {
                            RT r = Wy * R[sidx(L, x.t, x.u - 1, x.v)];
                            if constexpr (x.u >= 2) r += RT(x.u - 1) * R[sidx(L, x.t, x.u - 2, x.v)];
                            R[h] = r;
                        } else {
                            RT r = Wx * R[sidx(L, x.t - 1, x.u, x.v)];
                            if constexpr (x.t >= 2) r += RT(x.t - 1) * R[sidx(L, x.t - 2, x.u, x.v)];
                            R[h] = r;
                        }
                    }
                });
                R[0] = bp[n] * Fb[n];
            });
            const RT* Fcd = Ftab + m.F_off[pk];
            // step 2: ket Hermite sums, contracted over the ket primitives
            sfor<SC>([&](auto CI) {
                constexpr int c = decltype(CI)::value, cc = c0 + c;
                sfor<SD>([&](auto EI) {
                    constexpr int e = decltype(EI)::value, ee = e0 + e;
                    sfor<NHCD>([&](auto HC) {
                        constexpr int hc = decltype(HC)::value;
                        constexpr Tuv y = tuvOf(Lcd, hc);
                        if constexpr (hnz(LC, cc, LD, ee, y.t, y.u, y.v)) {
                            const RT ev = (((y.t + y.u + y.v) & 1) ? -pref : pref) * __ldg(Fcd + (cc * nD + ee) * NHCD + hc);
                            sfor<NHAB>([&](auto HA) {
                                constexpr int ha = decltype(HA)::value;
                                constexpr Tuv x = tuvOf(Lab, ha);
                                kc[ha * NG + c * SD + e] += ev * R[sidx(L, x.t + y.t, x.u + y.u, x.v + y.v)];
                            });
                        }
                    });
                });
            });
        }
        // step 3: bra Hermite sums
        const RT* Fab = Ftab + m.F_off[pb];
        sfor<NAB>([&](auto ABI) {
            constexpr int ab = decltype(ABI)::value, a = ab / nB, b = ab % nB;
            sfor<NHAB>([&](auto HA) {
                constexpr int ha = decltype(HA)::value;
                constexpr Tuv x = tuvOf(Lab, ha);
                if constexpr (hnz(LA, a, LB, b, x.t, x.u, x.v)) {
                    const RT ev = __ldg(Fab + ab * NHAB + ha);
#pragma unroll
                    for (int g = 0; g < NG; ++g) acc[ab * NG + g] += ev * kc[ha * NG + g];
                }
            });
        });
    }

    // step 4: digest (weights of digestBlock; partial sums of a slice go straight to the FP64 atomics)
    const size_t n = (size_t)d.n_out;
    const RT* __restrict__ P = nullptr;
    if constexpr (F32) P = m.P32;
    else P = d.P;
    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B], oC = d.sh_out_first[C] + c0, oD = d.sh_out_first[D] + e0;
    const double s4 = 4.0 * digestScale(A, B, C, D, i, j), s2 = 0.5 * s4;
#define QG3MD_V(a, b, c, e) acc[((a) * nB + (b)) * NG + (c) * SD + (e)]
#pragma unroll
    for (int a = 0; a < nA; ++a)
#pragma unroll
        for (int b = 0; b < nB; ++b) {  // Jt_ab += 4 w P_cd
            RT s = RT(0);
#pragma unroll
            for (int c = 0; c < SC; ++c)
#pragma unroll
                for (int e = 0; e < SD; ++e) s += QG3MD_V(a, b, c, e) * P[(oC + c) * n + oD + e];
            atomicAdd(&d.J[(oA + a) * n + oB + b], s4 * (double)s);
        }
#pragma unroll
    for (int c = 0; c < SC; ++c)
#pragma unroll
        for (int e = 0; e < SD; ++e) {  // Jt_cd += 4 w P_ab
            RT s = RT(0);
#pragma unroll
            for (int a = 0; a < nA; ++a)
#pragma unroll
                for (int b = 0; b < nB; ++b) s += QG3MD_V(a, b, c, e) * P[(oA + a) * n + oB + b];
            atomicAdd(&d.J[(oC + c) * n + oD + e], s4 * (double)s);
        }
#pragma unroll
    for (int a = 0; a < nA; ++a)
#pragma unroll
        for (int c = 0; c < SC; ++c) {  // Kt_ac += 2 w P_bd
            RT s = RT(0);
#pragma unroll
            for (int b = 0; b < nB; ++b)
#pragma unroll
                for (int e = 0; e < SD; ++e) s += QG3MD_V(a, b, c, e) * P[(oB + b) * n + oD + e];
            atomicAdd(&d.K[(oA + a) * n + oC + c], s2 * (double)s);
        }
#pragma unroll
    for (int b = 0; b < nB; ++b)
#pragma unroll
        for (int e = 0; e < SD; ++e) {  // Kt_bd += 2 w P_ac
            RT s = RT(0);
#pragma unroll
            for (int a = 0; a < nA; ++a)
#pragma unroll
                for (int c = 0; c < SC; ++c) s += QG3MD_V(a, b, c, e) * P[(oA + a) * n + oC + c];
            atomicAdd(&d.K[(oB + b) * n + oD + e], s2 * (double)s);
        }
#pragma unroll
    for (int a = 0; a < nA; ++a)
#pragma unroll
        for (int e = 0; e < SD; ++e) {  // Kt_ad += 2 w P_bc
            RT s = RT(0);
#pragma unroll
            for (int b = 0; b < nB; ++b)
#pragma unroll
                for (int c = 0; c < SC; ++c) s += QG3MD_V(a, b, c, e) * P[(oB + b) * n + oC + c];
            atomicAdd(&d.K[(oA + a) * n + oD + e], s2 * (double)s);
        }
#pragma unroll
    for (int b = 0; b < nB; ++b)
#pragma unroll
        for (int c = 0; c < SC; ++c) {  // Kt_bc += 2 w P_ad
            RT s = RT(0);
#pragma unroll
            for (int a = 0; a < nA; ++a)
#pragma unroll
                for (int e = 0; e < SD; ++e) s += QG3MD_V(a, b, c, e) * P[(oA + a) * n + oD + e];
            atomicAdd(&d.K[(oB + b) * n + oC + c], s2 * (double)s);
        }
#undef QG3MD_V
}

template <typename RT, int LA, int LB, int LC, int LD, int PART>
__device__ __forceinline__ void mdDispatch(int part, const DevTables& d, const MdTables& m, int i, int j, double pthr)
{
    if constexpr (PART < mdSplit(LA, LB, LC, LD)) {
        if (part == PART) mdBody<RT, LA, LB, LC, LD, PART>(d, m, i, j, pthr);
        else mdDispatch<RT, LA, LB, LC, LD, PART + 1>(part, d, m, i, j, pthr);
    }
}

constexpr int kMdBlock = 32;

template <typename RT, int LA, int LB, int LC, int LD>
__global__ void __launch_bounds__(kMdBlock) eriMdKernel(DevTables d, MdTables m, QSrc src, long long nq)
{
    const long long qi = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (qi >= nq) return;
    int i, j;
    double pthr;
    if (!quartetOf<true>(d, src, qi, i, j, pthr)) return;
    if (blockIdx.y == 0) countQuartet(d.counter);
    mdDispatch<RT, LA, LB, LC, LD, 0>((int)blockIdx.y, d, m, i, j, pthr);
}

template <typename RT, int LA, int LB, int LC, int LD>
bool mdLaunchT(const DevTables& d, const MdTables& m, const QSrc& src, long long nq, cudaStream_t st)
{
    constexpr int ns = mdSplit(LA, LB, LC, LD);
    const long long nb = (nq + kMdBlock - 1) / kMdBlock;
    eriMdKernel<RT, LA, LB, LC, LD><<<dim3((unsigned)nb, ns, 1), kMdBlock, 0, st>>>(d, m, src, nq);
    return true;
}

} // namespace
