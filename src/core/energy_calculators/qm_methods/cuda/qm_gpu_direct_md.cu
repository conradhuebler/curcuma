/*
 * < CUDA integral-direct J/K: thread-per-quartet kernels for the large classes (Q-G3) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 20 2026, docs/QM_GPU_ROADMAP.md section 2 "Large systems").
// FP64 instances of the thread-per-quartet McMurchie-Davidson kernels (algorithm, modes and
// screening: qm_gpu_direct_md_impl.cuh), the setup of their Hermite tables (mdInit) and the class
// dispatch (mdLaunch).

#include "qm_gpu_direct_md_impl.cuh"

#include <cmath>
#include <string>
#include <vector>

namespace {

/// F tables: one thread per (primitive pair, a'b'). Claude Generated (Oct 20 2026).
__global__ void mdBuildFKernel(DevTables d, const int* __restrict__ ppA, const int* __restrict__ ppB,
                               const long long* __restrict__ F_off, double* __restrict__ F, int npp, int maxab)
{
    const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= (long long)npp * maxab) return;
    const int k = (int)(g / maxab), ab = (int)(g % maxab);
    const int A = ppA[k], B = ppB[k];
    const int ncA = d.sh_ncomp[A], ncB = d.sh_ncomp[B];
    if (ncA > 6 || ncB > 6) return;  // stage 6b: f/g pairs have no md table (generic kernel, qm_gpu_eri_generic.cu)
    const int LA = ncA == 1 ? 0 : ncA == 3 ? 1 : 2, LB = ncB == 1 ? 0 : ncB == 3 ? 1 : 2;
    const int oa = nOut(LA), ob = nOut(LB);
    if (ab >= oa * ob) return;
    const int am = ab / ob, bm = ab % ob;
    const int Lab = LA + LB, nh = nSim(Lab);
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
            int h = 0;
            for (int t = 0; t <= Lab; ++t)
                for (int u = 0; u <= Lab - t; ++u)
                    for (int v = 0; v <= Lab - t - u; ++v, ++h) {
                        if (t > la[0] + lb[0] || u > la[1] + lb[1] || v > la[2] + lb[2]) continue;
                        out[h] += w * E[t * sAB_t + la[0] * sAB_i + lb[0]] * E[szAB + u * sAB_t + la[1] * sAB_i + lb[1]]
                                * E[2 * szAB + v * sAB_t + la[2] * sAB_i + lb[2]];
                    }
        }
    }
}

} // namespace

namespace qmgpu {

bool mdInit(const qmint::QMGpuEriTables* t, const DevTables& d, MdTables& m, std::vector<void*>& owned,
            std::string& why, size_t& bytes, double* dropped)
{
    m = MdTables{};
    bytes = 0;
    if (dropped) *dropped = 0.0;
    // component order (compile-time pattern) and 5d d shells with the transform pattern of inSupp
    for (int sh = 0; sh < t->nshell; ++sh) {
        const int Ls = t->sh_L[sh];
        if (Ls > 2) continue;  // stage 6b: f/g shells are never in an md class
        for (int c = 0; c < t->sh_ncomp[sh]; ++c)
            for (int ax = 0; ax < 3; ++ax)
                if (t->lmn[t->sh_lmn_off[sh] + 3 * c + ax] != cpow(Ls, c, ax)) { why = "component order"; return false; }
        if (Ls == 2) {
            if (t->sh_T_off[sh] < 0 || t->sh_out_width[sh] != 5) { why = "cartesian d shell"; return false; }
            // Stage 6b (Nov 2026): the S-metric orthonormalisation of buildSphericalTransform leaves rounding-level
            // entries outside the 5d pattern for some contracted d shells (seen in def2-TZVP/QZVPP; exactly zero in
            // the s/p/d sets checked before). Up to 1e-14 they are accepted (the md kernels drop them, the generic and
            // CPU kernels keep them); anything larger keeps the old kernels.
            for (int a = 0; a < 6; ++a)
                for (int mm = 0; mm < 5; ++mm) {
                    const double v = std::fabs(t->T[t->sh_T_off[sh] + a * 5 + mm]);
                    if (v == 0.0 || inSupp(2, mm, a)) continue;
                    if (v > 1.0e-14) { why = "5d transform pattern"; return false; }
                    if (dropped && v > *dropped) *dropped = v;
                }
        } else if (t->sh_T_off[sh] >= 0) {
            why = "transform on an s/p shell";
            return false;
        }
    }
    const int npp = t->nprimpair;
    std::vector<int> ppA((size_t)npp), ppB((size_t)npp);
    std::vector<long long> off((size_t)npp + 1);
    long long tot = 0;
    for (int i = 0; i < t->npair; ++i) {
        const int A = t->pr_A[i], B = t->pr_B[i];
        const int LA = t->sh_L[A], LB = t->sh_L[B];
        const long long sz = (LA > 2 || LB > 2) ? 0 : (long long)nOut(LA) * nOut(LB) * nSim(LA + LB);
        for (int k = t->pr_prim_off[i]; k < t->pr_prim_off[i] + t->pr_nprim[i]; ++k) {
            ppA[k] = A;
            ppB[k] = B;
            off[k] = tot;
            tot += sz;
        }
    }
    off[npp] = tot;
    double* Fd = nullptr;
    if (tot > 0 && cudaMalloc(&Fd, (size_t)tot * sizeof(double)) != cudaSuccess) {
        cudaGetLastError();
        why = "device memory for the Hermite tables (" + std::to_string(tot * 8 / 1000000) + " MB)";
        return false;
    }
    if (Fd) owned.push_back(Fd);
    int *dA = nullptr, *dB = nullptr;
    long long* doff = nullptr;
    if (!upload(&dA, ppA.data(), (size_t)npp, owned) || !upload(&dB, ppB.data(), (size_t)npp, owned)
        || !upload(&doff, off.data(), (size_t)npp, owned)) {
        why = "upload of the Hermite table offsets failed";
        return false;
    }
    if (npp > 0) {
        const int maxab = 25;
        const long long nthr = (long long)npp * maxab;
        mdBuildFKernel<<<(unsigned)((nthr + 127) / 128), 128>>>(d, dA, dB, doff, Fd, npp, maxab);
        if (cudaDeviceSynchronize() != cudaSuccess || cudaGetLastError() != cudaSuccess) {
            why = "Hermite table kernel failed";
            return false;
        }
    }
    m.F = Fd;
    m.F_off = doff;
    m.nF = tot;
    m.nBoys = (long long)t->boys_npts * (t->boys_ntab + 1);
    bytes = (size_t)tot * sizeof(double) + (size_t)npp * (2 * sizeof(int) + sizeof(long long));
    return true;
}

bool mdLaunch(int LA, int LB, int LC, int LD, const DevTables& d, const MdTables& m, const QSrc& src, long long nq,
              cudaStream_t st)
{
    if (nq <= 0) return true;
#define QG3MD_CASE(a, b, c, e) \
    case ((a * 3 + b) * 3 + c) * 3 + e: return mdLaunchT<double, a, b, c, e>(d, m, src, nq, st);
    switch (((LA * 3 + LB) * 3 + LC) * 3 + LD) {
        // every class with bra pair class X <= ket pair class Y (the launches of qm_gpu_direct.cu)
        QG3MD_CASE(0, 0, 0, 0) QG3MD_CASE(0, 0, 0, 1) QG3MD_CASE(0, 0, 0, 2) QG3MD_CASE(0, 0, 1, 0) QG3MD_CASE(0, 0, 1, 1) QG3MD_CASE(0, 0, 1, 2) QG3MD_CASE(0, 0, 2, 0) QG3MD_CASE(0, 0, 2, 1) QG3MD_CASE(0, 0, 2, 2)
        QG3MD_CASE(0, 1, 0, 1) QG3MD_CASE(0, 1, 0, 2) QG3MD_CASE(0, 1, 1, 0) QG3MD_CASE(0, 1, 1, 1) QG3MD_CASE(0, 1, 1, 2) QG3MD_CASE(0, 1, 2, 0) QG3MD_CASE(0, 1, 2, 1) QG3MD_CASE(0, 1, 2, 2)
        QG3MD_CASE(0, 2, 0, 2) QG3MD_CASE(0, 2, 1, 0) QG3MD_CASE(0, 2, 1, 1) QG3MD_CASE(0, 2, 1, 2) QG3MD_CASE(0, 2, 2, 0) QG3MD_CASE(0, 2, 2, 1) QG3MD_CASE(0, 2, 2, 2)
        QG3MD_CASE(1, 0, 1, 0) QG3MD_CASE(1, 0, 1, 1) QG3MD_CASE(1, 0, 1, 2) QG3MD_CASE(1, 0, 2, 0) QG3MD_CASE(1, 0, 2, 1) QG3MD_CASE(1, 0, 2, 2)
        QG3MD_CASE(1, 1, 1, 1) QG3MD_CASE(1, 1, 1, 2) QG3MD_CASE(1, 1, 2, 0) QG3MD_CASE(1, 1, 2, 1) QG3MD_CASE(1, 1, 2, 2)
        QG3MD_CASE(1, 2, 1, 2) QG3MD_CASE(1, 2, 2, 0) QG3MD_CASE(1, 2, 2, 1) QG3MD_CASE(1, 2, 2, 2)
        QG3MD_CASE(2, 0, 2, 0) QG3MD_CASE(2, 0, 2, 1) QG3MD_CASE(2, 0, 2, 2)
        QG3MD_CASE(2, 1, 2, 1) QG3MD_CASE(2, 1, 2, 2)
        QG3MD_CASE(2, 2, 2, 2)
    default: return false;
    }
#undef QG3MD_CASE
}

} // namespace qmgpu
