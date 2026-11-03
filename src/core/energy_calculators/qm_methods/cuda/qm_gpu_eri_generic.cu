/*
 * < CUDA ERI kernels for quartet classes with f and g shells (stage 6b) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Nov 2026, stage 6b of docs/QM_GPU_ROADMAP.md). AI-generated, machine-tested.
//
// One generic (runtime-L) McMurchie-Davidson kernel for every shell quartet class with an f or g shell,
// used by the stored device tensor (Q-G2, DIRECT = false) and the integral-direct J/K (Q-G3, DIRECT =
// true). The integrals are those of the CPU kernel mdQuartet / mdQuartetT (qm_integrals.cpp) in the
// Hermite-matrix form (Helgaker, Jorgensen, Olsen ch. 9.9.4; McMurchie, Davidson, J. Comput. Phys. 26,
// 218 (1978)):
//   (a'b'|c'e') = sum_pb sum_h F[pb][a'b'][h] Kc[h][c'e'],
//   Kc[h][c'e'] = sum_pk pref(pb,pk) sum_h' (-1)^{|h'|} F[pk][c'e'][h'] R^0_{h + h'}(pb, pk),
//   pref = 2 pi^(5/2) / (p q sqrt(p+q)), R^0 from the Boys function (grid + Taylor table of the CPU),
// F = Hermite expansion of the active component pairs (contraction and spherical transform folded in,
// genInit). Screening: quartetOf (the Q-G2 list or the Q-G3 density screen) and the primitive break of
// the other kernels, so the quartet set is that of the s/p/d kernels and of the CPU.
//
// Work decomposition (correctness first): one warp (small classes) or one block of 128 threads (large classes) per
// shell quartet, grid-stride over the quartets of a launch; the pair with the higher L is computed on the bra side
// (swap, decided per class); R^0 level by level in shared memory (degree-ordered simplex, t+u+v <= 16);
// Kc and the quartet block in a per-group slot of global memory (an (gg|gg) quartet has 81 x 81 active
// components, too large for shared memory). The stored tensor writes each element once (canonical orbit
// representative, set8); the direct digest is digestBlock of qm_gpu_eri_kernels.cuh (FP64 atomics).

#include "qm_gpu_eri_generic.cuh"

#include <cuda_runtime.h>

#include <algorithm>
#include <string>
#include <vector>

namespace {

constexpr int kGenMaxL = 4;                     // g shells
constexpr int kGenMaxLtot = 4 * kGenMaxL;       // (gg|gg)
constexpr int kGenNS = (kGenMaxLtot + 1) * (kGenMaxLtot + 2) * (kGenMaxLtot + 3) / 6;  // 969
constexpr int kGenBlock = 128;
constexpr long long kGenWarpWork = 40000;  // active component quartets x Hermite simplex, warp vs block per quartet

__host__ __device__ inline int gSim(int L) { return L < 0 ? 0 : (L + 1) * (L + 2) * (L + 3) / 6; }
/// Degree-ordered simplex index of (t,u,v): all entries of degree < t+u+v first, inside a degree by t, then u.
__host__ __device__ inline int gIdx(int t, int u, int v)
{
    const int dg = t + u + v;
    return gSim(dg - 1) + t * (dg + 1) - t * (t - 1) / 2 + u;
}

/// F tables: one thread per (primitive pair k, active component pair a'b'). Same expansion as mdBuildFKernel
/// (qm_gpu_direct_md.cu) with runtime L, a general transform block and the degree order of h.
__global__ void genBuildFKernel(DevTables d, const int* __restrict__ ppA, const int* __restrict__ ppB,
                                const long long* __restrict__ F_off, const int* __restrict__ sh_L,
                                const int* __restrict__ width, double* __restrict__ F, int npp, int maxab)
{
    const long long gidx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (gidx >= (long long)npp * maxab) return;
    const int k = (int)(gidx / maxab), ab = (int)(gidx % maxab);
    const int A = ppA[k], B = ppB[k];
    const int LA = sh_L[A], LB = sh_L[B];
    const int ncA = d.sh_ncomp[A], ncB = d.sh_ncomp[B];
    const int oa = width[A], ob = width[B];
    if (ab >= oa * ob) return;
    const int am = ab / ob, bm = ab % ob;
    const int Lab = LA + LB, nh = gSim(Lab);
    const int sAB_t = (LA + 1) * (LB + 1), sAB_i = LB + 1, szAB = (Lab + 1) * sAB_t;
    const double* E = d.E + d.pp_E_off[k];
    const int ia = d.pp_ia[k], ib = d.pp_ib[k];
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B];
    const int tA = d.sh_T_off[A], tB = d.sh_T_off[B];
    double* out = F + F_off[k] + (long long)ab * nh;
    for (int h = 0; h < nh; ++h) out[h] = 0.0;
    for (int a = 0; a < ncA; ++a) {
        const double ta = tA < 0 ? (a == am ? 1.0 : 0.0) : d.T[tA + a * oa + am];
        if (ta == 0.0) continue;
        const int* la = d.lmn + d.sh_lmn_off[A] + 3 * a;
        for (int b = 0; b < ncB; ++b) {
            const double tb = tB < 0 ? (b == bm ? 1.0 : 0.0) : d.T[tB + b * ob + bm];
            if (tb == 0.0) continue;
            const int* lb = d.lmn + d.sh_lmn_off[B] + 3 * b;
            const double w = ta * tb * d.coef[d.sh_coef_off[A] + a * npA + ia] * d.coef[d.sh_coef_off[B] + b * npB + ib];
            for (int t = 0; t <= la[0] + lb[0]; ++t) {
                const double ex = E[t * sAB_t + la[0] * sAB_i + lb[0]];
                for (int u = 0; u <= la[1] + lb[1]; ++u) {
                    const double exy = ex * E[szAB + u * sAB_t + la[1] * sAB_i + lb[1]];
                    for (int v = 0; v <= la[2] + lb[2]; ++v)
                        out[gIdx(t, u, v)] += w * exy * E[2 * szAB + v * sAB_t + la[2] * sAB_i + lb[2]];
                }
            }
        }
    }
}

template <int GS>
__device__ __forceinline__ void genSync()
{
    if constexpr (GS > 32) __syncthreads();
    else __syncwarp();
}

/// GS threads per shell quartet (GS = kGenBlock: one block, the large classes; GS = 32: one warp, kGenBlock / 32
/// quartets per block, the many small f-type classes), grid-stride, see the file header. Dynamic shared memory per
/// group: two R level buffers of gSim(L) + 2 (L + 1) + 4 doubles; the (t,u,v) table is shared by the block.
/// Claude Generated (Nov 2026).
template <bool DIRECT, int GS>
__global__ void __launch_bounds__(kGenBlock) eriGenericKernel(DevTables d, qmgpu::GenTables g, QSrc src, long long nq, int LA,
                                                              int LB, int LC, int LD, double* __restrict__ out, int swap)
{
    extern __shared__ double gsh[];
    __shared__ int sM[kGenNS];  // (t,u,v) of degree-ordered index k, packed t | u << 8 | v << 16
    constexpr int NG = kGenBlock / GS;  // quartet groups per block
    const int grp = threadIdx.x / GS, tid = threadIdx.x % GS, bs = GS;
    // swap (host decision per class): compute the quartet as (CD|AB), the pair with the higher L on the bra side, so
    // that the ket Hermite step per primitive quartet stays small; J/K digest and set8 are symmetric in the two pairs.
    if (swap) {
        int t = LA; LA = LC; LC = t;
        t = LB; LB = LD; LD = t;
    }
    const int L = LA + LB + LC + LD, Lab = LA + LB, Lcd = LC + LD;
    const int NS = gSim(L), NHAB = gSim(Lab), NHCD = gSim(Lcd);
    for (int dg = 0; dg <= L; ++dg)
        for (int x = threadIdx.x; x < (dg + 1) * (dg + 2) / 2; x += blockDim.x) {
            int t = 0, r = x;  // x-th (t, u) pair of degree dg, by t then u
            while (r > dg - t) { r -= dg - t + 1; ++t; }
            const int u = r, v = dg - t - u;
            sM[gIdx(t, u, v)] = t | (u << 8) | (v << 16);
        }
    __syncthreads();
    double* base = gsh + (size_t)grp * (2 * NS + 2 * (L + 1) + 4);
    double* sR[2] = { base, base + NS };
    double* sF = base + 2 * NS;
    double* sBp = sF + (L + 1);
    double* sW = sBp + (L + 1);
    double& sPref = sW[3];
    const long long gslot = (long long)blockIdx.x * NG + grp;
    double* Kc = g.scratch + gslot * g.slot;

    for (long long qi = gslot; qi < nq; qi += (long long)gridDim.x * NG) {
        int i, j;
        double pthr;
        if (!quartetOf<DIRECT>(d, src, qi, i, j, pthr)) continue;  // uniform per group
        if (swap) { const int t = i; i = j; j = t; }
        if constexpr (DIRECT) {
            if (tid == 0) atomicAdd(d.counter, 1ull);
        }
        const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
        const int nA = g.width[A], nB = g.width[B], nC = g.width[C], nD = g.width[D];
        const int nab = nA * nB, ncd = nC * nD;
        double* acc = Kc + (long long)NHAB * ncd;
        for (int k = tid; k < nab * ncd; k += bs) acc[k] = 0.0;
        const int pb0 = d.pr_prim_off[i], nbp = d.pr_nprim[i];
        const int pk0 = d.pr_prim_off[j], nkp = d.pr_nprim[j];
        const double ket_max = nkp > 0 ? d.pp_schwarz[pk0] : 0.0;
        for (int ib = 0; ib < nbp; ++ib) {
            const int pb = pb0 + ib;
            const double sb = d.pp_schwarz[pb];
            if (sb * ket_max < pthr) break;  // same primitive screen as the other kernels
            genSync<GS>();                   // previous bra step is done with Kc
            for (int k = tid; k < NHAB * ncd; k += bs) Kc[k] = 0.0;
            for (int ik = 0; ik < nkp; ++ik) {
                const int pk = pk0 + ik;
                if (sb * d.pp_schwarz[pk] < pthr) break;
                genSync<GS>();  // previous ket step is done with R
                if (tid == 0) {
                    const double p = d.pp_p[pb], qq = d.pp_p[pk];
                    const double rho = p * qq / (p + qq);
                    const double Wx = d.pp_P[3 * pb] - d.pp_P[3 * pk];
                    const double Wy = d.pp_P[3 * pb + 1] - d.pp_P[3 * pk + 1];
                    const double Wz = d.pp_P[3 * pb + 2] - d.pp_P[3 * pk + 2];
                    boysDev(L, rho * (Wx * Wx + Wy * Wy + Wz * Wz), sF, d);
                    const double bb = -2.0 * rho;
                    sBp[0] = 1.0;
                    for (int n = 1; n <= L; ++n) sBp[n] = sBp[n - 1] * bb;
                    sW[0] = Wx; sW[1] = Wy; sW[2] = Wz;
                    sPref = d.pref / (p * qq * sqrt(p + qq));
                    sR[L & 1][0] = sBp[L] * sF[L];
                }
                genSync<GS>();
                // R^n_tuv on the simplex t+u+v <= L-n from level n+1 (recursion of buildR0Simplex: z, y, x)
                for (int n = L - 1; n >= 0; --n) {
                    const double* c = sR[(n + 1) & 1];
                    double* x = sR[n & 1];
                    const int top = gSim(L - n);
                    for (int k = tid; k < top; k += bs) {
                        double r;
                        if (k == 0) {
                            r = sBp[n] * sF[n];
                        } else {
                            const int tv = sM[k];
                            const int t = tv & 255, u = (tv >> 8) & 255, v = tv >> 16;
                            if (v > 0)
                                r = (v >= 2 ? (v - 1) * c[gIdx(t, u, v - 2)] : 0.0) + sW[2] * c[gIdx(t, u, v - 1)];
                            else if (u > 0)
                                r = (u >= 2 ? (u - 1) * c[gIdx(t, u - 2, v)] : 0.0) + sW[1] * c[gIdx(t, u - 1, v)];
                            else
                                r = (t >= 2 ? (t - 1) * c[gIdx(t - 2, u, v)] : 0.0) + sW[0] * c[gIdx(t - 1, u, v)];
                        }
                        x[k] = r;
                    }
                    genSync<GS>();
                }
                const double* R = sR[0];
                const double pref = sPref;
                const double* Fcd = g.F + g.F_off[pk];
                // ket Hermite sums, contracted over the ket primitives (each entry owned by one thread)
                for (int e = tid; e < NHAB * ncd; e += bs) {
                    const int ha = e / ncd, cd = e % ncd;
                    const int ta = sM[ha];
                    const int t0 = ta & 255, u0 = (ta >> 8) & 255, v0 = ta >> 16;
                    const double* f = Fcd + (long long)cd * NHCD;
                    double s = 0.0;
                    for (int hc = 0; hc < NHCD; ++hc) {
                        const double fv = f[hc];
                        if (fv == 0.0) continue;
                        const int tc = sM[hc];
                        const int t1 = tc & 255, u1 = (tc >> 8) & 255, v1 = tc >> 16;
                        const double rv = R[gIdx(t0 + t1, u0 + u1, v0 + v1)];
                        s += ((t1 + u1 + v1) & 1) ? -fv * rv : fv * rv;
                    }
                    Kc[e] += pref * s;
                }
            }
            genSync<GS>();
            // bra Hermite sums
            const double* Fab = g.F + g.F_off[pb];
            for (int e = tid; e < nab * ncd; e += bs) {
                const int ab = e / ncd, cd = e % ncd;
                const double* f = Fab + (long long)ab * NHAB;
                double s = 0.0;
                for (int ha = 0; ha < NHAB; ++ha) {
                    const double fv = f[ha];
                    if (fv != 0.0) s += fv * Kc[(long long)ha * ncd + cd];
                }
                acc[e] += s;
            }
        }
        genSync<GS>();
        const int oA = d.sh_out_first[A], oB = d.sh_out_first[B], oC = d.sh_out_first[C], oD = d.sh_out_first[D];
        if constexpr (DIRECT) {
            const int dim[4] = { nA, nB, nC, nD };
            digestBlock(acc, dim, oA, oB, oC, oD, digestScale(A, B, C, D, i, j), d, tid, bs);
        } else {
            // canonical representative of every 8-fold orbit inside this quartet -> set8 (as eriClassKernel)
            const size_t n = (size_t)d.n_out;
            for (int idx = tid; idx < nab * ncd; idx += bs) {
                const int e = idx % nD, c = (idx / nD) % nC, b = (idx / (nD * nC)) % nB, a = idx / (nD * nC * nB);
                if (A == B && a < b) continue;
                if (C == D && c < e) continue;
                if (i == j && a * nB + b < c * nD + e) continue;
                set8(out, n, oA + a, oB + b, oC + c, oD + e, acc[idx]);
            }
        }
        genSync<GS>();  // the next quartet of this group reuses acc and Kc
    }
}

} // namespace

namespace qmgpu {

bool genInit(const qmint::QMGpuEriTables* t, const DevTables& d, GenTables& g, std::vector<void*>& owned,
             std::string& why)
{
    g = GenTables{};
    for (int s = 0; s < t->nshell; ++s)
        if (t->sh_L[s] > kGenMaxL) {
            why = "shell with L > 4 (device kernels handle s-g)";
            return false;
        }
    const int npp = t->nprimpair;
    std::vector<int> ppA((size_t)npp), ppB((size_t)npp);
    std::vector<long long> off((size_t)npp + 1);
    long long tot = 0;
    int maxab = 1, maxh = 1;
    for (int i = 0; i < t->npair; ++i) {
        const int A = t->pr_A[i], B = t->pr_B[i];
        const int wab = t->sh_out_width[A] * t->sh_out_width[B];
        const int nh = gSim(t->sh_L[A] + t->sh_L[B]);
        maxab = std::max(maxab, wab);
        maxh = std::max(maxh, nh);
        for (int k = t->pr_prim_off[i]; k < t->pr_prim_off[i] + t->pr_nprim[i]; ++k) {
            ppA[k] = A;
            ppB[k] = B;
            off[k] = tot;
            tot += (long long)wab * nh;
        }
    }
    off[npp] = tot;
    double* Fd = nullptr;
    if (tot > 0 && cudaMalloc(&Fd, (size_t)tot * sizeof(double)) != cudaSuccess) {
        cudaGetLastError();
        why = "device memory for the f/g Hermite tables (" + std::to_string(tot * 8 / 1000000) + " MB)";
        return false;
    }
    if (Fd) owned.push_back(Fd);
    int *dA = nullptr, *dB = nullptr, *dL = nullptr, *dW = nullptr;
    long long* doff = nullptr;
    if (!upload(&dA, ppA.data(), (size_t)npp, owned) || !upload(&dB, ppB.data(), (size_t)npp, owned)
        || !upload(&doff, off.data(), (size_t)npp, owned) || !upload(&dL, t->sh_L, (size_t)t->nshell, owned)
        || !upload(&dW, t->sh_out_width, (size_t)t->nshell, owned)) {
        why = "upload of the f/g Hermite table offsets failed";
        return false;
    }
    if (npp > 0) {
        const long long nthr = (long long)npp * maxab;
        genBuildFKernel<<<(unsigned)((nthr + 127) / 128), 128>>>(d, dA, dB, doff, dL, dW, Fd, npp, maxab);
        if (cudaDeviceSynchronize() != cudaSuccess || cudaGetLastError() != cudaSuccess) {
            why = "f/g Hermite table kernel failed";
            return false;
        }
    }
    // scratch: Kc (Hermite x ket components) + the quartet block, per resident block
    const long long slot = (long long)maxh * maxab + (long long)maxab * maxab;
    int dev = 0, sms = 16;
    cudaGetDevice(&dev);
    cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev);
    size_t free_b = 0, total_b = 0;
    cudaMemGetInfo(&free_b, &total_b);
    const double budget = std::min(256.0e6, 0.25 * (double)free_b);
    // slots = quartet groups: 4 blocks per SM of kGenBlock / 32 warp groups each
    int nslots = std::max(kGenBlock / 32, std::min(4 * sms * (kGenBlock / 32), (int)(budget / ((double)slot * sizeof(double)))));
    nslots -= nslots % (kGenBlock / 32);
    double* sc = nullptr;
    if (cudaMalloc(&sc, (size_t)nslots * slot * sizeof(double)) != cudaSuccess) {
        cudaGetLastError();
        why = "device memory for the f/g ERI scratch";
        return false;
    }
    owned.push_back(sc);
    g.F = Fd;
    g.F_off = doff;
    g.sh_L = dL;
    g.width = dW;
    g.scratch = sc;
    g.slot = slot;
    g.nslots = nslots;
    g.bytes = (size_t)tot * sizeof(double) + (size_t)nslots * slot * sizeof(double);
    return true;
}

bool genLaunch(int LA, int LB, int LC, int LD, const DevTables& d, const GenTables& g, const QSrc& src, long long nq,
               bool direct, double* out, cudaStream_t st)
{
    if (nq <= 0) return true;
    if (LA > kGenMaxL || LB > kGenMaxL || LC > kGenMaxL || LD > kGenMaxL || !g.scratch) return false;
    // one warp per quartet for the classes with few active component quartets (most f-type classes), one block of
    // kGenBlock threads per quartet for the large ones (threshold by measurement, RTX 5090)
    const int swap = (LC + LD > LA + LB) ? 1 : 0;
    const int L = LA + LB + LC + LD;
    const long long ncomp4 = (long long)(2 * LA + 1) * (2 * LB + 1) * (2 * LC + 1) * (2 * LD + 1);
    const bool warp = ncomp4 * gSim(L) <= kGenWarpWork;
    const int ngrp = warp ? kGenBlock / 32 : 1;
    const size_t shm = (size_t)ngrp * (2 * gSim(L) + 2 * (L + 1) + 4) * sizeof(double);
    const unsigned nb = (unsigned)std::max<long long>(1, std::min<long long>((nq + ngrp - 1) / ngrp, g.nslots / ngrp));
    if (direct) {
        if (warp) eriGenericKernel<true, 32><<<nb, kGenBlock, shm, st>>>(d, g, src, nq, LA, LB, LC, LD, out, swap);
        else eriGenericKernel<true, kGenBlock><<<nb, kGenBlock, shm, st>>>(d, g, src, nq, LA, LB, LC, LD, out, swap);
    } else {
        if (warp) eriGenericKernel<false, 32><<<nb, kGenBlock, shm, st>>>(d, g, src, nq, LA, LB, LC, LD, out, swap);
        else eriGenericKernel<false, kGenBlock><<<nb, kGenBlock, shm, st>>>(d, g, src, nq, LA, LB, LC, LD, out, swap);
    }
    return cudaGetLastError() == cudaSuccess;
}

} // namespace qmgpu
