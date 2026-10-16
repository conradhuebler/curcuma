/*
 * < CUDA 3-index part of the RI (density-fitting) gradient of the native QM engine >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, RI gradient). AI-generated, machine-tested.
//
// Device form of qmri::RIJK::threeIndexGradient / derivBlock (qm_ri.cpp, section 6):
//   for every canonical orbital shell pair (A,B) (one CUDA block) and every auxiliary shell C of its
//   screened task list: for every bra primitive pair, the ket Hermite sum
//     G_c(t,u,v) = sum_ket-prims pref c_C sum_tau,nu,phi (-1)^(tau+nu+phi) E^C_tau E^C_nu E^C_phi R_{t+tau,u+nu,v+phi}
//   (McMurchie-Davidson, Helgaker 9.9.33, R from the Boys function, Helgaker 9.9.13), then per bra
//   component pair H = sum_c W_abc G_c and the centre derivatives
//     d/dA_k (ab|c) = 2 alpha (a+1_k b|c) - a_k (a-1_k b|c)   (same for B),
//   accumulated with the host weights W (3-index density, cartesian). d/dC = -(d/dA + d/dB).
// One block per bra pair, one warp per task (auxiliary shell): the gradient on the atoms of A and B is
// summed in registers over all tasks, the auxiliary-centre part goes to global memory with atomicAdd
// (order not deterministic: the last bits vary between runs, as for the other device gradients). Shared
// memory per warp: R and its recursion buffer on the compact Hermite simplex t+u+v <= Ltot <= 9, G
// (<= 15 aux components x 56). The host checks the angular momentum limits.

#include "../qm_gpu_ri_grad_backend.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr int kMaxLtot = 9;                    // orbital d raised (2+2+1) + auxiliary g (4)
constexpr int kS = kMaxLtot + 1;               // R cube edge
constexpr int kNSimp = (kMaxLtot + 1) * (kMaxLtot + 2) * (kMaxLtot + 3) / 6;  // 220 Hermite indices t+u+v <= 9
constexpr int kNSimpB = 56;                    // t+u+v <= 5 (bra: orbital d raised by one)
constexpr int kMaxNC = 15;                     // cartesian components of a g shell
constexpr int kThreads = 128;

struct Dev {
    int natoms, nshell, aux_base, npairs;
    const int *sh_L, *sh_atom, *sh_ncomp, *sh_lmn_off, *lmn, *sh_nprim, *sh_coef_off, *sh_exp_off;
    const double *coef, *exps;
    const int *pair_A, *pair_B, *pair_pp_off, *pair_npp, *aux_pp_off, *aux_npp;
    const int *pp_ia, *pp_ib;
    const double* pp_geom;
    const long long* pp_E_off;
    const double* E;
    const long long* task_start;
    const int* task_aux;
    const long long* task_woff;
    const double* W;
    const double* boys;
};

// Boys function F_0..F_mmax(T): the scheme of qm_integrals.cpp boysArrayInto -- 7-term Taylor expansion
// of F_mmax around the nearest grid point of a table (dT = 0.05, T < 36, built on the host in long
// double by downward recursion), downward recursion to F_0; T >= 36: asymptotic F_0, upward recursion.
constexpr double kBoysDT = 0.05, kBoysTmax = 36.0;
constexpr int kBoysNTAB = kMaxLtot + 7;  // n = 0..kBoysNTAB tabulated
__device__ void boysRI(int mmax, double T, double* F, const double* tab)
{
    const double eT = exp(-T);
    if (T >= kBoysTmax) {
        F[0] = 0.5 * sqrt(3.14159265358979323846 / T);
        for (int m = 0; m < mmax; ++m) F[m + 1] = ((2 * m + 1) * F[m] - eT) / (2.0 * T);
        return;
    }
    const int i = (int)(T / kBoysDT + 0.5);
    const double dd = i * kBoysDT - T;
    const double* fi = tab + (size_t)i * (kBoysNTAB + 1) + mmax;
    const double inv[7] = { 1.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4, 1.0 / 5, 1.0 / 6 };
    double acc = fi[6];
    for (int k = 5; k >= 0; --k) acc = fi[k] + dd * inv[k + 1] * acc;
    F[mmax] = acc;
    for (int m = mmax - 1; m >= 0; --m) F[m] = (2.0 * T * F[m + 1] + eT) / (2 * m + 1);
}

__global__ void __launch_bounds__(kThreads) riGrad3cKernel(Dev d, double* grad)
{
    const int k = blockIdx.x;
    const long long t0 = d.task_start[k], t1 = d.task_start[k + 1];
    if (t0 == t1) return;
    constexpr int kWarps = kThreads / 32;
    // Hermite indices (t,u,v) with t+u+v <= 9 stored compactly, ordered by t+u+v (then t, u): the entries
    // with t+u+v <= m are a prefix of length sLev[m], for R (m = Ltot) and for G (m = LAB) alike.
    // Every warp works on its own task (auxiliary shell) with its own R, recursion buffer and G.
    __shared__ double sR[kWarps][kNSimp], sL[kWarps][kNSimp], sG[kWarps][kMaxNC * kNSimpB];
    __shared__ double sBF[kWarps][kMaxLtot + 2];
    __shared__ short sR3[3 * kNSimp];     // (t,u,v) of compact index x
    __shared__ short sIdx[kS * kS * kS];  // compact index of (t,u,v), cube stride kS
    __shared__ int sLev[kMaxLtot + 2];
    __shared__ double sRed[6][kWarps];
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;

    const int A = d.pair_A[k], B = d.pair_B[k];
    const int LA = d.sh_L[A], LB = d.sh_L[B], nA = d.sh_ncomp[A], nB = d.sh_ncomp[B];
    const int aA = d.sh_atom[A], aB = d.sh_atom[B];
    const double wpair = (A != B) ? 2.0 : 1.0;
    const int LAB = LA + LB + 1;
    const int sAB_t = (LA + 2) * (LB + 2), sAB_i = LB + 2;
    const int* lmnA = d.lmn + d.sh_lmn_off[A];
    const int* lmnB = d.lmn + d.sh_lmn_off[B];
    const double* coefA = d.coef + d.sh_coef_off[A];
    const double* coefB = d.coef + d.sh_coef_off[B];
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B];
    const double* expA = d.exps + d.sh_exp_off[A];
    const double* expB = d.exps + d.sh_exp_off[B];

    if (tid == 0) {
        int c = 0;
        for (int m = 0; m <= kMaxLtot; ++m) {
            for (int t = 0; t <= m; ++t)
                for (int u = 0; u <= m - t; ++u) {
                    const int v = m - t - u;
                    sR3[3 * c] = (short)t;
                    sR3[3 * c + 1] = (short)u;
                    sR3[3 * c + 2] = (short)v;
                    sIdx[(t * kS + u) * kS + v] = (short)c;
                    ++c;
                }
            sLev[m] = c;
        }
    }
    __syncthreads();
    const int nsimp = sLev[LAB];  // bra simplex size, also the stride of G
    double* wR = sR[warp];
    double* wL = sL[warp];
    double* wG = sG[warp];
    double* wBF = sBF[warp];
    // per-lane accumulators of the 6 derivative components (A_xyz, B_xyz) over the warp's tasks
    double accAll[6] = { 0, 0, 0, 0, 0, 0 };

    for (long long ti = t0 + warp; ti < t1; ti += kWarps) {
        const int C = d.aux_base + d.task_aux[ti];
        const int LC = d.sh_L[C], nC = d.sh_ncomp[C], aC = d.sh_atom[C];
        const int* lmnC = d.lmn + d.sh_lmn_off[C];
        const double* coefC = d.coef + d.sh_coef_off[C];
        const int npC = d.sh_nprim[C];
        const double* Wt = d.W + d.task_woff[ti];
        const int Ltot = LAB + LC;
        const int sCD_t = LC + 1;  // (LC+1)*(0+1); sCD_i = 1
        const int kp0 = d.aux_pp_off[C - d.aux_base], nkp = d.aux_npp[C - d.aux_base];
        double acc[6] = { 0, 0, 0, 0, 0, 0 };

        for (int pbi = 0; pbi < d.pair_npp[k]; ++pbi) {
            const int pb = d.pair_pp_off[k] + pbi;
            const int ia = d.pp_ia[pb], ib = d.pp_ib[pb];
            const double p = d.pp_geom[4 * pb], Px = d.pp_geom[4 * pb + 1], Py = d.pp_geom[4 * pb + 2], Pz = d.pp_geom[4 * pb + 3];
            const long long esz = (long long)(LA + LB + 3) * sAB_t;
            const double* Ex = d.E + d.pp_E_off[pb];
            const double* Ey = Ex + esz;
            const double* Ez = Ey + esz;
            __syncwarp();
            for (int x = lane; x < nC * nsimp; x += 32) wG[x] = 0.0;
            for (int pki = 0; pki < nkp; ++pki) {
                const int pk = kp0 + pki;
                const int ka = d.pp_ia[pk];
                const double q = d.pp_geom[4 * pk];
                const double rho = p * q / (p + q);
                const double Wx = Px - d.pp_geom[4 * pk + 1], Wy = Py - d.pp_geom[4 * pk + 2], Wz = Pz - d.pp_geom[4 * pk + 3];
                __syncwarp();
                if (lane == 0) {
                    double F[kMaxLtot + 1];
                    boysRI(Ltot, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F, d.boys);
                    double bp = 1.0;
                    for (int n = 0; n <= Ltot; ++n) {
                        wBF[n] = bp * F[n];
                        bp *= -2.0 * rho;
                    }
                    wR[0] = wBF[Ltot];
                }
                __syncwarp();
                // R^0_tuv by levels n = Ltot-1 .. 0 (Helgaker 9.9.18-9.9.20), as qmri coulombR0; level n
                // fills the simplex t+u+v <= Ltot-n from the previous level
                double* cur = wR;
                double* nxt = wL;
                for (int n = Ltot - 1; n >= 0; --n) {
                    const int cnt = sLev[Ltot - n];
                    for (int x = lane; x < cnt; x += 32) {
                        const int t = sR3[3 * x], u = sR3[3 * x + 1], v = sR3[3 * x + 2];
                        const int base = (t * kS + u) * kS + v;
                        double r;
                        if (x == 0) r = wBF[n];
                        else if (v > 0) r = (v >= 2 ? (v - 1) * cur[sIdx[base - 2]] : 0.0) + Wz * cur[sIdx[base - 1]];
                        else if (u > 0) r = (u >= 2 ? (u - 1) * cur[sIdx[base - 2 * kS]] : 0.0) + Wy * cur[sIdx[base - kS]];
                        else r = (t >= 2 ? (t - 1) * cur[sIdx[base - 2 * kS * kS]] : 0.0) + Wx * cur[sIdx[base - kS * kS]];
                        nxt[x] = r;
                    }
                    __syncwarp();
                    double* tmp = cur;
                    cur = nxt;
                    nxt = tmp;
                }
                const double pref = 2.0 * pow(3.14159265358979323846, 2.5) / (p * q * sqrt(p + q));
                const long long ksz = (long long)(LC + 1) * sCD_t;
                const double* Kx = d.E + d.pp_E_off[pk];
                const double* Ky = Kx + ksz;
                const double* Kz = Ky + ksz;
                for (int x = lane; x < nC * nsimp; x += 32) {
                    const int cc = x / nsimp, j = x % nsimp;
                    const int t = sR3[3 * j], u = sR3[3 * j + 1], v = sR3[3 * j + 2];
                    const int l0 = lmnC[3 * cc], l1 = lmnC[3 * cc + 1], l2 = lmnC[3 * cc + 2];
                    const double w = pref * coefC[cc * npC + ka];
                    double a1 = 0.0;
                    for (int k1 = 0; k1 <= l0; ++k1) {
                        const double e1 = ((k1 & 1) ? -1.0 : 1.0) * Kx[k1 * sCD_t + l0];
                        if (e1 == 0.0) continue;
                        double a2 = 0.0;
                        for (int k2 = 0; k2 <= l1; ++k2) {
                            const double e2 = ((k2 & 1) ? -1.0 : 1.0) * Ky[k2 * sCD_t + l1];
                            if (e2 == 0.0) continue;
                            const short* Ir = sIdx + ((t + k1) * kS + u + k2) * kS + v;
                            double a3 = 0.0;
                            for (int k3 = 0; k3 <= l2; ++k3) a3 += ((k3 & 1) ? -1.0 : 1.0) * Kz[k3 * sCD_t + l2] * cur[Ir[k3]];
                            a2 += e2 * a3;
                        }
                        a1 += e1 * a2;
                    }
                    wG[cc * nsimp + j] += w * a1;
                }
            }
            __syncwarp();
            // centre derivatives: one work item per (bra component pair, derivative component)
            const double alpha = expA[ia], beta = expB[ib];
            for (int x = lane; x < nA * nB * 6; x += 32) {
                const int role = x % 6, cab = x / 6, ca = cab / nB, cb = cab % nB;
                const double c_ab = coefA[ca * npA + ia] * coefB[cb * npB + ib];
                if (c_ab == 0.0) continue;
                const double* w = Wt + (size_t)cab * nC;
                int a[3] = { lmnA[3 * ca], lmnA[3 * ca + 1], lmnA[3 * ca + 2] };
                int b[3] = { lmnB[3 * cb], lmnB[3 * cb + 1], lmnB[3 * cb + 2] };
                const bool onB = role >= 3;
                const int axis = role % 3;
                int* sh = onB ? b : a;
                const int l = sh[axis];
                auto val = [&]() {
                    double sum = 0.0;
                    for (int t = 0; t <= a[0] + b[0]; ++t) {
                        const double et = Ex[t * sAB_t + a[0] * sAB_i + b[0]];
                        if (et == 0.0) continue;
                        for (int u = 0; u <= a[1] + b[1]; ++u) {
                            const double etu = et * Ey[u * sAB_t + a[1] * sAB_i + b[1]];
                            if (etu == 0.0) continue;
                            double sv = 0.0;
                            for (int v = 0; v <= a[2] + b[2]; ++v) {
                                const double ez = Ez[v * sAB_t + a[2] * sAB_i + b[2]];
                                if (ez == 0.0) continue;
                                const double* gp = wG + sIdx[(t * kS + u) * kS + v];
                                double h = 0.0;
                                for (int cc = 0; cc < nC; ++cc) h += w[cc] * gp[cc * nsimp];
                                sv += ez * h;
                            }
                            sum += etu * sv;
                        }
                    }
                    return sum;
                };
                sh[axis] = l + 1;
                double dv = 2.0 * (onB ? beta : alpha) * val();
                if (l > 0) {
                    sh[axis] = l - 1;
                    dv -= l * val();
                }
                acc[role] += c_ab * dv;
            }
        }
        // task total of the 6 components -> the auxiliary centre gets -(dA + dB) (warp reduction)
        double gc[3] = { 0, 0, 0 };
        for (int r = 0; r < 6; ++r) {
            double v = acc[r];
            for (int off = 16; off > 0; off >>= 1) v += __shfl_down_sync(0xffffffffu, v, off);
            gc[r % 3] += v;
            accAll[r] += acc[r];
        }
        if (lane == 0)
            for (int x = 0; x < 3; ++x) atomicAdd(&grad[3 * aC + x], -wpair * gc[x]);
    }
    // block totals on the atoms of A and B
    for (int r = 0; r < 6; ++r) {
        double v = accAll[r];
        for (int off = 16; off > 0; off >>= 1) v += __shfl_down_sync(0xffffffffu, v, off);
        if (lane == 0) sRed[r][warp] = v;
    }
    __syncthreads();
    if (tid < 6) {
        double tot = 0.0;
        for (int w = 0; w < kWarps; ++w) tot += sRed[tid][w];
        atomicAdd(&grad[3 * (tid >= 3 ? aB : aA) + tid % 3], wpair * tot);
    }
}

template <typename T, typename U>
bool upload(std::vector<void*>& owned, T** dst, const U* src, size_t count)
{
    static_assert(sizeof(T) == sizeof(U), "upload: element size mismatch");
    void* p = nullptr;
    if (cudaMalloc(&p, sizeof(T) * std::max<size_t>(1, count)) != cudaSuccess) return false;
    owned.push_back(p);
    if (count && cudaMemcpy(p, src, sizeof(T) * count, cudaMemcpyHostToDevice) != cudaSuccess) return false;
    *dst = static_cast<T*>(p);
    return true;
}

void setWhyRiGrad(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

} // namespace

// 3-index part of the RI gradient on the device (natoms x 3, Eh/Bohr, added into grad which the caller
// zeroes). Returns 1 on success, 0 with a reason. Claude Generated (Oct 2026, RI gradient).
extern "C" int curcuma_cuda_qm_ri_grad(const qmint::QMGpuRIGradInput* in, double* grad, qmint::QMGpuRIGradOutput* out,
                                       int device, char* why, int whylen)
{
    if (!in || !grad || in->npairs <= 0) { setWhyRiGrad(why, whylen, "empty RI gradient input"); return 0; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyRiGrad(why, whylen, "cudaSetDevice failed"); return 0; }
    for (int s = 0; s < in->nshell; ++s) {
        const int L = in->sh_L[s];
        if ((s < in->aux_base && L > 2) || L > 4 || in->sh_ncomp[s] > kMaxNC) {
            setWhyRiGrad(why, whylen, "angular momentum above the device limit (orbital d, auxiliary g)");
            return 0;
        }
    }
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<void*> owned;
    auto cleanup = [&]() {
        for (void* p : owned) cudaFree(p);
        cudaGetLastError();
    };
    Dev d{};
    d.natoms = in->natoms;
    d.nshell = in->nshell;
    d.aux_base = in->aux_base;
    d.npairs = in->npairs;
    const int nsh = in->nshell, naux = in->nshell - in->aux_base;
    double* dgrad = nullptr;
    std::vector<double> zero(3 * (size_t)in->natoms, 0.0);
    static const std::vector<double> boys = [] {
        // same construction as qm_integrals.cpp BoysTable: long-double downward recursion from far above
        const int npts = (int)std::lround(kBoysTmax / kBoysDT) + 2, M = kBoysNTAB + 200;
        std::vector<double> f((size_t)npts * (kBoysNTAB + 1));
        std::vector<long double> G(M + 2);
        for (int i = 0; i < npts; ++i) {
            const long double T = (long double)i * kBoysDT, eT = expl(-T);
            G[M + 1] = 0.0L;
            for (int n = M; n >= 0; --n) G[n] = (2.0L * T * G[n + 1] + eT) / (2.0L * n + 1.0L);
            for (int n = 0; n <= kBoysNTAB; ++n) f[(size_t)i * (kBoysNTAB + 1) + n] = (double)G[n];
        }
        return f;
    }();
    bool ok = upload(owned, &d.sh_L, in->sh_L, nsh) && upload(owned, &d.sh_atom, in->sh_atom, nsh)
        && upload(owned, &d.sh_ncomp, in->sh_ncomp, nsh) && upload(owned, &d.sh_lmn_off, in->sh_lmn_off, nsh)
        && upload(owned, &d.lmn, in->lmn, (size_t)in->nlmn) && upload(owned, &d.sh_nprim, in->sh_nprim, nsh)
        && upload(owned, &d.sh_coef_off, in->sh_coef_off, nsh) && upload(owned, &d.coef, in->coef, (size_t)in->ncoef)
        && upload(owned, &d.sh_exp_off, in->sh_exp_off, nsh) && upload(owned, &d.exps, in->exps, (size_t)in->nexp)
        && upload(owned, &d.pair_A, in->pair_A, in->npairs) && upload(owned, &d.pair_B, in->pair_B, in->npairs)
        && upload(owned, &d.pair_pp_off, in->pair_pp_off, in->npairs) && upload(owned, &d.pair_npp, in->pair_npp, in->npairs)
        && upload(owned, &d.aux_pp_off, in->aux_pp_off, naux) && upload(owned, &d.aux_npp, in->aux_npp, naux)
        && upload(owned, &d.pp_ia, in->pp_ia, in->npp) && upload(owned, &d.pp_ib, in->pp_ib, in->npp)
        && upload(owned, &d.pp_geom, in->pp_geom, 4 * (size_t)in->npp) && upload(owned, &d.pp_E_off, in->pp_E_off, in->npp)
        && upload(owned, &d.E, in->E, (size_t)in->nE) && upload(owned, &d.task_start, in->task_start, (size_t)in->npairs + 1)
        && upload(owned, &d.task_aux, in->task_aux, (size_t)in->ntask) && upload(owned, &d.task_woff, in->task_woff, (size_t)in->ntask)
        && upload(owned, &d.W, in->W, (size_t)in->nW) && upload(owned, &d.boys, boys.data(), boys.size()) && upload(owned, &dgrad, zero.data(), zero.size());
    if (!ok) {
        cleanup();
        setWhyRiGrad(why, whylen, "device allocation/upload failed (RI gradient tables and weights)");
        return 0;
    }
    // load the kernel module before the timed launch (lazy loading, and the PTX JIT on devices newer than
    // the compiled architectures, would otherwise be counted as kernel time)
    cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr, riGrad3cKernel);
    cudaDeviceSynchronize();
    const auto t1 = std::chrono::steady_clock::now();
    riGrad3cKernel<<<in->npairs, kThreads>>>(d, dgrad);
    cudaError_t e = cudaDeviceSynchronize();
    if (e == cudaSuccess) e = cudaGetLastError();
    if (e != cudaSuccess) {
        cleanup();
        setWhyRiGrad(why, whylen, std::string("RI gradient kernel failed: ") + cudaGetErrorString(e));
        return 0;
    }
    const auto t2 = std::chrono::steady_clock::now();
    if (cudaMemcpy(zero.data(), dgrad, sizeof(double) * zero.size(), cudaMemcpyDeviceToHost) != cudaSuccess) {
        cleanup();
        setWhyRiGrad(why, whylen, "RI gradient download failed");
        return 0;
    }
    for (size_t i = 0; i < zero.size(); ++i) grad[i] += zero[i];
    cleanup();
    if (out) {
        out->ms_upload = std::chrono::duration<double, std::milli>(t1 - t0).count();
        out->ms_kernel = std::chrono::duration<double, std::milli>(t2 - t1).count();
    }
    return 1;
}
