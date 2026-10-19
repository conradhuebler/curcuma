/*
 * < CUDA integral-direct J/K of the native QM engine (Q-G3) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, Q-G3 of docs/QM_GPU_ROADMAP.md).
//
// Integral-direct Coulomb/exchange build on the device: no stored ERI tensor, device memory
// O(n^2) (P, Pmax, J, K) plus the O(npair) shell-pair tables. Device counterpart of
// qmint::DirectJK::build (qm_integrals.cpp; Almloef, Faegri, Korsell, J. Comput. Chem. 3, 385
// (1982); Haeser, Ahlrichs, J. Comput. Chem. 10, 104 (1989)):
//   1. Pmax(A,B) = max |P| over the shell block (A,B), Pall = max Pmax  (kernel pmaxKernel)
//   2. every canonical shell quartet with
//        Q_i Q_j max(Pmax over AB, CD, AC, AD, BC, BD) >= thr
//      is computed with the Q-G2 kernels (qm_gpu_eri_kernels.cuh, DIRECT = true) and digested
//      at once into the unsymmetrised Jt, Kt (FP64 atomicAdd, weight deg/8);
//   3. J = (Jt + Jt^T)/2, K = (Kt + Kt^T)/2.
// Quartet enumeration (Oct 19 2026, replaces the pairs(X) x pairs(Y) rectangle whose O(npair^2)
// launch cost dominated beyond ~1000 basis functions): the pairs of each pair class are sorted by
// Schwarz factor, largest first. Per build and class pair X <= Y the host walks the bra list and keeps
// for each bra i only the ket prefix with Q_i Q_j Pall >= thr (a superset of the exact screen, since
// pq <= Pall); for X == Y the kets start at the bra's own position. Every unordered pair {i, j} is
// enumerated once, so launches cover ~ the Schwarz-surviving quartets and nothing else.
// The incremental-on-dP scheme and the full-rebuild cadence are the host SCF's (qm_scf.cpp),
// unchanged: this builder is linear in P. Atomic accumulation makes the last bits of J/K depend
// on the thread schedule (not run-to-run bit-identical).

#include "qm_gpu_eri_kernels.cuh"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

/// Pmax(A,B) = Pmax(B,A) = max |P| over the shell block, for A <= B (as DirectJK).
/// Also Pall = max over all blocks (atomicMax on the bit pattern: non-negative doubles order like
/// their unsigned 64-bit images).
__global__ void pmaxKernel(const double* __restrict__ P, const int* __restrict__ first, const int* __restrict__ width,
                           int ns, int n, double* __restrict__ Pmax, unsigned long long* __restrict__ pall_bits)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (long long)ns * ns) return;
    const int A = (int)(k / ns), B = (int)(k % ns);
    if (A > B) return;
    double m = 0.0;
    for (int a = 0; a < width[A]; ++a)
        for (int b = 0; b < width[B]; ++b) m = fmax(m, fabs(P[(size_t)(first[A] + a) * n + first[B] + b]));
    Pmax[(size_t)A * ns + B] = m;
    Pmax[(size_t)B * ns + A] = m;
    atomicMax(pall_bits, (unsigned long long)__double_as_longlong(m));
}

/// Q-G4: J = (Jt + Jt^T)/2, K = (Kt + Kt^T)/2 into separate output buffers (device side of the
/// host loop in QmGpuDirectJkCuda::build). Claude Generated (Oct 2026).
__global__ void symmetriseKernel(const double* __restrict__ Jt, const double* __restrict__ Kt,
                                 double* __restrict__ J, double* __restrict__ K, int n)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (size_t)n * n) return;
    const size_t a = k / n, b = k % n, kt = b * n + a;
    J[k] = 0.5 * (Jt[k] + Jt[kt]);
    K[k] = 0.5 * (Kt[k] + Kt[kt]);
}

/// Two-step McMurchie-Davidson variant of eriClassKernel for the large classes of the
/// integral-direct build (Claude Generated, Oct 19 2026). Profile (C20H42, RTX 5090): the classes
/// with more than 27 cartesian components (d shells, (pp|pp)) were ~13 % of the quartets but ~83 %
/// of the J/K time, because eriClassKernel evaluates, per primitive quartet and per component
/// quartet, the full double Hermite sum sum_tuv E^ab_tuv sum_tau,nu,phi (+-) E^cd R_{t+tau,..}.
/// Here the ket sum is done once per ket component pair (c, e) and bra Hermite index h = (t,u,v),
/// t+u+v <= LA+LB:
///   Kh[h][ce] = sum_{tau,nu,phi} (-1)^(tau+nu+phi) E^cd_tau E^cd_nu E^cd_phi R_{t+tau,u+nu,v+phi}
/// and then (ab|ce) += w sum_{t<=tx,u<=ty,v<=tz} E^ab_t E^ab_u E^ab_v Kh[h(t,u,v)][ce]
/// (Helgaker, Jorgensen, Olsen ch. 9.9.4; the same two-step order as the CPU gradient kernel).
/// Same integrals as eriClassKernel up to the summation order (rounding level); screening,
/// Boys, R recursion, 6d->5d and the digest are unchanged. One block per quartet. (Tried and
/// dropped, C20H42 on the RTX 5090: one warp per quartet, and warps over the primitive quartets
/// with shared atomics -- both slower than this block-wide form.)
__global__ void eriHermiteKernel(DevTables d, QSrc src, int LA, int LB, int LC, int LD)
{
    extern __shared__ double smem[];
    __shared__ double sF[16], sBp[16], sW[3], sPref;
    __shared__ int hmap[125];   // (t,u,v) -> bra Hermite index, t,u,v <= 4
    __shared__ int hlist[35 * 3];

    int i, j;
    double pthr;
    if (!quartetOf<true>(d, src, (long long)blockIdx.x, i, j, pthr)) return;  // uniform per block
    if (threadIdx.x == 0) atomicAdd(d.counter, 1ull);
    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const int nA = d.sh_ncomp[A], nB = d.sh_ncomp[B], nC = d.sh_ncomp[C], nD = d.sh_ncomp[D];
    const int L = LA + LB + LC + LD, Lab = LA + LB;
    const int s = L + 1, s2 = s * s, s3 = s2 * s;
    const int nab = nA * nB, ncd = nC * nD, ncart = nab * ncd;
    const int nh = (Lab + 1) * (Lab + 2) * (Lab + 3) / 6;
    const int tid = threadIdx.x, bs = blockDim.x;
    double* Rb[2] = { smem, smem + s3 };
    double* blk = smem + 2 * s3;
    double* tmp = blk + ncart;
    double* Kh = tmp + ncart;  // nh * ncd

    if (tid == 0) {
        int h = 0;
        for (int t = 0; t <= Lab; ++t)
            for (int u = 0; u <= Lab - t; ++u)
                for (int v = 0; v <= Lab - t - u; ++v) {
                    hmap[(t * 5 + u) * 5 + v] = h;
                    hlist[3 * h] = t; hlist[3 * h + 1] = u; hlist[3 * h + 2] = v;
                    ++h;
                }
    }
    for (int k = tid; k < ncart; k += bs) blk[k] = 0.0;

    const int sAB_t = (LA + 1) * (LB + 1), sAB_i = LB + 1, szAB = (LA + LB + 1) * sAB_t;
    const int sCD_t = (LC + 1) * (LD + 1), sCD_i = LD + 1, szCD = (LC + LD + 1) * sCD_t;
    const int pb0 = d.pr_prim_off[i], nbp = d.pr_nprim[i];
    const int pk0 = d.pr_prim_off[j], nkp = d.pr_nprim[j];
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B], npC = d.sh_nprim[C], npD = d.sh_nprim[D];

    const double ket_max = nkp > 0 ? d.pp_schwarz[pk0] : 0.0;
    for (int ib = 0; ib < nbp; ++ib) {
        const int pb = pb0 + ib;
        const double sb = d.pp_schwarz[pb];
        if (sb * ket_max < pthr) break;  // same primitive screen as eriClassKernel
        for (int ik = 0; ik < nkp; ++ik) {
            const int pk = pk0 + ik;
            if (sb * d.pp_schwarz[pk] < pthr) break;
            __syncthreads();
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
                Rb[L & 1][0] = sBp[L] * sF[L];
            }
            __syncthreads();
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
            const int ia = d.pp_ia[pb], ibb = d.pp_ib[pb], ic = d.pp_ia[pk], idd = d.pp_ib[pk];
            const double* Eab = d.E + d.pp_E_off[pb];
            const double* Ecd = d.E + d.pp_E_off[pk];

            // step 1: ket Hermite sums per (h, ce)
            for (int k = tid; k < nh * ncd; k += bs) {
                const int h = k / ncd, ce = k % ncd, c = ce / nD, e = ce % nD;
                const int t = hlist[3 * h], u = hlist[3 * h + 1], v = hlist[3 * h + 2];
                const int* lc = d.lmn + d.sh_lmn_off[C] + 3 * c;
                const int* ld = d.lmn + d.sh_lmn_off[D] + 3 * e;
                const int kx = lc[0] + ld[0], ky = lc[1] + ld[1], kz = lc[2] + ld[2];
                double kv = 0.0;
                for (int tau = 0; tau <= kx; ++tau) {
                    const double fx = ((tau & 1) ? -1.0 : 1.0) * Ecd[tau * sCD_t + lc[0] * sCD_i + ld[0]];
                    if (fx == 0.0) continue;
                    for (int nu = 0; nu <= ky; ++nu) {
                        const double ftn = fx * ((nu & 1) ? -1.0 : 1.0) * Ecd[szCD + nu * sCD_t + lc[1] * sCD_i + ld[1]];
                        if (ftn == 0.0) continue;
                        const double* Rr = R + ((t + tau) * s + (u + nu)) * s + v;
                        for (int ph = 0; ph <= kz; ++ph)
                            kv += ftn * ((ph & 1) ? -1.0 : 1.0) * Ecd[2 * szCD + ph * sCD_t + lc[2] * sCD_i + ld[2]] * Rr[ph];
                    }
                }
                Kh[k] = kv;
            }
            __syncthreads();
            // step 2: bra Hermite sums per component quartet
            const double pref = sPref;
            for (int k = tid; k < ncart; k += bs) {
                const int ab = k / ncd, ce = k % ncd;
                const int a = ab / nB, b = ab % nB, c = ce / nD, e = ce % nD;
                const double w = pref * d.coef[d.sh_coef_off[A] + a * npA + ia] * d.coef[d.sh_coef_off[B] + b * npB + ibb]
                               * d.coef[d.sh_coef_off[C] + c * npC + ic] * d.coef[d.sh_coef_off[D] + e * npD + idd];
                if (w == 0.0) continue;
                const int* la = d.lmn + d.sh_lmn_off[A] + 3 * a;
                const int* lb = d.lmn + d.sh_lmn_off[B] + 3 * b;
                const int tx = la[0] + lb[0], ty = la[1] + lb[1], tz = la[2] + lb[2];
                double val = 0.0;
                for (int t = 0; t <= tx; ++t) {
                    const double ex = Eab[t * sAB_t + la[0] * sAB_i + lb[0]];
                    if (ex == 0.0) continue;
                    for (int u = 0; u <= ty; ++u) {
                        const double etu = ex * Eab[szAB + u * sAB_t + la[1] * sAB_i + lb[1]];
                        if (etu == 0.0) continue;
                        for (int v = 0; v <= tz; ++v)
                            val += etu * Eab[2 * szAB + v * sAB_t + la[2] * sAB_i + lb[2]] * Kh[hmap[(t * 5 + u) * 5 + v] * ncd + ce];
                    }
                }
                blk[k] += w * val;
            }
        }
    }
    __syncthreads();

    int dim[4] = { nA, nB, nC, nD };
    const int s4[4] = { A, B, C, D };
    double* x = blk;
    double* y = tmp;
    for (int ax = 0; ax < 4; ++ax) {
        const int toff = d.sh_T_off[s4[ax]];
        if (toff < 0) continue;
        const int nin = dim[ax];
        const int nout = (nin == 6) ? 5 : nin;
        int outer = 1, inner = 1;
        for (int m = 0; m < ax; ++m) outer *= dim[m];
        for (int m = ax + 1; m < 4; ++m) inner *= dim[m];
        const int total = outer * nout * inner;
        for (int idx = tid; idx < total; idx += bs) {
            const int r = idx % inner, jj = (idx / inner) % nout, o = idx / (inner * nout);
            double acc = 0.0;
            for (int a2 = 0; a2 < nin; ++a2) {
                const double tt = d.T[toff + a2 * nout + jj];
                if (tt == 0.0) continue;
                acc += tt * x[((size_t)o * nin + a2) * inner + r];
            }
            y[idx] = acc;
        }
        __syncthreads();
        double* sw = x; x = y; y = sw;
        dim[ax] = nout;
    }
    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B], oC = d.sh_out_first[C], oD = d.sh_out_first[D];
    digestBlock(x, dim, oA, oB, oC, oD, digestScale(A, B, C, D, i, j), d, tid, bs);
}

/// Launch eriHermiteKernel for one class (chunked at the grid limit).
void launchHermiteClass(int LA, int LB, int LC, int LD, const DevTables& d, const QSrc& src, long long nq, cudaStream_t st)
{
    const int L = LA + LB + LC + LD, s3 = (L + 1) * (L + 1) * (L + 1), Lab = LA + LB;
    const int ncd = ncomp(LC) * ncomp(LD);
    const int ncart = ncomp(LA) * ncomp(LB) * ncd;
    const int nh = (Lab + 1) * (Lab + 2) * (Lab + 3) / 6;
    const size_t shm = (size_t)(2 * s3 + 2 * ncart + nh * ncd) * sizeof(double);
    int bs = 32;
    while (bs < std::max(ncart, nh * ncd) && bs < 128) bs *= 2;
    const long long chunk = 2000000000LL;
    for (long long off = 0; off < nq; off += chunk) {
        const long long cnt = (nq - off < chunk) ? nq - off : chunk;
        QSrc s = src;
        s.q0 = src.q0 + off;
        eriHermiteKernel<<<(unsigned)cnt, bs, shm, st>>>(d, s, LA, LB, LC, LD);
    }
}

class QmGpuDirectJkCuda final : public qmint::QMGpuJK {
public:
    QmGpuDirectJkCuda() = default;
    ~QmGpuDirectJkCuda() override
    {
        for (void* p : m_owned) cudaFree(p);
        if (m_rows_dev_off) { cudaFree(m_rows_dev_off); cudaFree(m_rows_dev_bra); cudaFree(m_rows_dev_ket0); }
        for (auto& s : m_streams)
            if (s) cudaStreamDestroy(s);
        if (m_e0) cudaEventDestroy(m_e0);
        if (m_e1) cudaEventDestroy(m_e1);
    }

    bool init(const qmint::QMGpuEriTables* t, int device, std::string& why)
    {
        if (!t->pc_off || !t->pr_schwarz) { why = "tables without pair classes (built with list_quartets)"; return false; }
        if (!uploadCommonTables(t, m_d, m_owned, why)) return false;
        m_n = t->n_out;
        m_ns = t->nshell;
        m_d.thr = t->screening;
        m_spherical = false;
        const char* tenv = std::getenv("CURCUMA_QM_GPU_ERI_BLOCK_ONLY");  // diagnostic: block kernel everywhere
        m_thread_ok = threadKernelUsable(t, m_spherical) && !(tenv && tenv[0] == '1');
        m_pc_off.assign(t->pc_off, t->pc_off + 10);
        // pairs of each class sorted by Schwarz factor, largest first (stable: CPU order among equals)
        m_sorted.assign(t->pc_pairs, t->pc_pairs + t->npair);
        m_sorted_q.resize((size_t)t->npair);
        for (int X = 0; X < 9; ++X) {
            std::stable_sort(m_sorted.begin() + m_pc_off[X], m_sorted.begin() + m_pc_off[X + 1],
                             [&](int a, int b) { return t->pr_schwarz[a] > t->pr_schwarz[b]; });
            for (long long k = m_pc_off[X]; k < m_pc_off[X + 1]; ++k) m_sorted_q[k] = t->pr_schwarz[m_sorted[k]];
        }
        const char* penv = std::getenv("CURCUMA_QM_GPU_DIRECT_PROFILE");
        m_profile = penv && penv[0] == '1';
        const char* benv = std::getenv("CURCUMA_QM_GPU_DIRECT_BLOCK_KERNEL");  // diagnostic: old block kernel
        m_block_only = benv && benv[0] == '1';
        const size_t n2 = (size_t)m_n * m_n;
        bool ok = upload(const_cast<double**>(&m_d.pr_schwarz), t->pr_schwarz, (size_t)t->npair, m_owned)
            && upload(&m_sorted_dev, m_sorted.data(), (size_t)t->npair, m_owned)
            && upload(&m_width, t->sh_out_width, (size_t)t->nshell, m_owned)
            && alloc(&m_P, n2) && alloc(&m_J, n2) && alloc(&m_K, n2) && alloc(&m_Pmax, (size_t)m_ns * m_ns)
            && alloc(&m_counter, 1) && alloc(&m_pall, 1);
        if (!ok) { why = "device allocation / upload for the direct J/K failed"; return false; }
        m_d.P = m_P;
        m_d.Pmax = m_Pmax;
        m_d.J = m_J;
        m_d.K = m_K;
        m_d.counter = m_counter;
        m_first = const_cast<int*>(m_d.sh_out_first);
        for (auto& s : m_streams)
            if (cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking) != cudaSuccess) { why = "stream creation failed"; return false; }
        cudaEventCreate(&m_e0);
        cudaEventCreate(&m_e1);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        m_bytes = 0;
        return cudaGetLastError() == cudaSuccess;
    }

    bool build(const double* P, double* J, double* K) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n2 = (size_t)m_n * m_n;
        if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        float kms = 0.0f;
        unsigned long long cnt = 0;
        if (!runKernels(kms)) return false;
        if (cudaMemcpy(J, m_J, n2 * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
            || cudaMemcpy(K, m_K, n2 * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
            || cudaMemcpy(&cnt, m_counter, sizeof cnt, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        // J = (Jt + Jt^T)/2, K = (Kt + Kt^T)/2, as DirectJK::build
        const int n = m_n;
        for (int a = 0; a < n; ++a)
            for (int b = a; b < n; ++b) {
                const double js = 0.5 * (J[(size_t)a * n + b] + J[(size_t)b * n + a]);
                const double ks = 0.5 * (K[(size_t)a * n + b] + K[(size_t)b * n + a]);
                J[(size_t)a * n + b] = J[(size_t)b * n + a] = js;
                K[(size_t)a * n + b] = K[(size_t)b * n + a] = ks;
            }
        m_last_quartets = (long long)cnt;
        ++m_builds;
        m_kernel_ms += kms;
        m_total_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    // Q-G4 (Claude Generated, Oct 2026): device-pointer variant for the resident SCF loop. The
    // symmetrisation runs on the device with the same expression as the host loop above
    // (0.5 * (a + b), commutative, so bit-identical to build() for the same Jt, Kt).
    bool buildDevice(const double* dP, double* dJ, double* dK) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n2 = (size_t)m_n * m_n;
        if (cudaMemcpy(m_P, dP, n2 * sizeof(double), cudaMemcpyDeviceToDevice) != cudaSuccess) return false;
        float kms = 0.0f;
        if (!runKernels(kms)) return false;
        symmetriseKernel<<<(unsigned)((n2 + 255) / 256), 256>>>(m_J, m_K, dJ, dK, m_n);
        unsigned long long cnt = 0;
        if (cudaGetLastError() != cudaSuccess
            || cudaMemcpy(&cnt, m_counter, sizeof cnt, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        m_last_quartets = (long long)cnt;
        ++m_builds;
        m_kernel_ms += kms;
        m_total_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    /// Q-G3 screening threshold of the following builds (host SCF: looser early phase); the
    /// primitive threshold follows as thr * 1e-3, as buildEriGpuTables. Claude Generated (Oct 2026).
    bool setScreening(double thr, bool prim_density) override
    {
        m_d.thr = thr;
        m_d.prim_thr = thr * 1.0e-3;
        m_d.prim_density = prim_density ? 1 : 0;
        return true;
    }

    /// Pmax + all quartet-class launches on m_P into the unsymmetrised m_J, m_K; synchronous.
    bool runKernels(float& kms)
    {
        const size_t n2 = (size_t)m_n * m_n;
        cudaMemset(m_J, 0, n2 * sizeof(double));
        cudaMemset(m_K, 0, n2 * sizeof(double));
        cudaMemset(m_counter, 0, sizeof(unsigned long long));
        cudaMemset(m_pall, 0, sizeof(unsigned long long));
        cudaEventRecord(m_e0);
        const long long nsn = (long long)m_ns * m_ns;
        pmaxKernel<<<(unsigned)((nsn + 127) / 128), 128>>>(m_P, m_first, m_width, m_ns, m_n, m_Pmax, m_pall);
        unsigned long long pall_bits = 0;
        if (cudaMemcpy(&pall_bits, m_pall, sizeof pall_bits, cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        double Pall;
        std::memcpy(&Pall, &pall_bits, sizeof Pall);

        // Ragged rows per class pair X <= Y (see the file header). The predicate !(Qi*Qj*Pall < thr)
        // is monotone in Qj (sorted lists), and Qi*Qj*pq <= Qi*Qj*Pall, so the prefix is a superset.
        const auto th0 = std::chrono::steady_clock::now();
        const double thr = m_d.thr;
        m_row_off.clear();
        m_row_bra.clear();
        m_row_ket0.clear();
        struct Launch { int X, Y; size_t row0; int nrow; long long nq; };
        std::vector<Launch> launches;
        long long launched = 0;
        for (int X = 0; X < 9; ++X) {
            const long long x0 = m_pc_off[X], nX = m_pc_off[X + 1] - x0;
            if (nX == 0) continue;
            for (int Y = X; Y < 9; ++Y) {
                const long long y0 = m_pc_off[Y], nY = m_pc_off[Y + 1] - y0;
                if (nY == 0) continue;
                const double* qX = m_sorted_q.data() + x0;
                const double* qY = m_sorted_q.data() + y0;
                Launch L{ X, Y, m_row_off.size(), 0, 0 };
                long long k = nY;  // ket prefix length, non-increasing along the bra list
                for (long long b = 0; b < nX; ++b) {
                    if (qX[b] * qY[0] * Pall < thr) break;
                    while (k > 0 && qX[b] * qY[k - 1] * Pall < thr) --k;
                    const long long start = (X == Y) ? b : 0;
                    if (k <= start) {
                        if (X == Y) break;  // k only shrinks while b grows
                        continue;
                    }
                    m_row_off.push_back(L.nq);
                    m_row_bra.push_back(m_sorted[x0 + b]);
                    m_row_ket0.push_back((int)start);
                    L.nq += k - start;
                }
                L.nrow = (int)(m_row_off.size() - L.row0);
                if (L.nrow > 0) {
                    launches.push_back(L);
                    launched += L.nq;
                }
            }
        }
        const size_t nrows = m_row_off.size();
        if (nrows > m_rows_cap) {
            if (m_rows_dev_off) { cudaFree(m_rows_dev_off); cudaFree(m_rows_dev_bra); cudaFree(m_rows_dev_ket0); }
            m_rows_cap = nrows + nrows / 4 + 1024;
            if (cudaMalloc(&m_rows_dev_off, m_rows_cap * sizeof(long long)) != cudaSuccess
                || cudaMalloc(&m_rows_dev_bra, m_rows_cap * sizeof(int)) != cudaSuccess
                || cudaMalloc(&m_rows_dev_ket0, m_rows_cap * sizeof(int)) != cudaSuccess) {
                m_rows_dev_off = nullptr; m_rows_dev_bra = m_rows_dev_ket0 = nullptr; m_rows_cap = 0;
                return false;
            }
        }
        if (nrows > 0
            && (cudaMemcpy(m_rows_dev_off, m_row_off.data(), nrows * sizeof(long long), cudaMemcpyHostToDevice) != cudaSuccess
                || cudaMemcpy(m_rows_dev_bra, m_row_bra.data(), nrows * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess
                || cudaMemcpy(m_rows_dev_ket0, m_row_ket0.data(), nrows * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess))
            return false;
        m_rows_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - th0).count();
        m_launched += launched;

        cudaEventRecord(m_e1);  // later launches on the streams wait for Pmax, the memsets and the rows
        for (auto& s : m_streams) cudaStreamWaitEvent(s, m_e1, 0);

        int launch = 0;
        for (const Launch& L : launches) {
            const int LA = L.X / 3, LB = L.X % 3, LC = L.Y / 3, LD = L.Y % 3;
            const QSrc src{ 0, m_rows_dev_off + L.row0, m_rows_dev_bra + L.row0, m_rows_dev_ket0 + L.row0,
                            m_sorted_dev + m_pc_off[L.Y], L.nrow };
            cudaStream_t st = m_profile ? (cudaStream_t)0 : m_streams[launch++ % kStreams];
            unsigned long long c0 = 0;
            cudaEvent_t pe0 = nullptr, pe1 = nullptr;
            if (m_profile) {
                cudaDeviceSynchronize();
                cudaMemcpy(&c0, m_counter, sizeof c0, cudaMemcpyDeviceToHost);
                cudaEventCreate(&pe0);
                cudaEventCreate(&pe1);
                cudaEventRecord(pe0, st);
            }
            const bool threaded = m_thread_ok && launchThreadClass<true>(LA, LB, LC, LD, m_d, src, L.nq, nullptr, st, m_spherical);
            if (!threaded) {
                if (m_block_only)
                    launchBlockClass<true>(LA, LB, LC, LD, m_d, src, L.nq, nullptr, st);
                else
                    launchHermiteClass(LA, LB, LC, LD, m_d, src, L.nq, st);
            }
            if (cudaGetLastError() != cudaSuccess) return false;
            if (m_profile) {
                cudaEventRecord(pe1, st);
                cudaEventSynchronize(pe1);
                float ms = 0.0f;
                cudaEventElapsedTime(&ms, pe0, pe1);
                unsigned long long c1 = 0;
                cudaMemcpy(&c1, m_counter, sizeof c1, cudaMemcpyDeviceToHost);
                ProfileEntry& e = m_prof[L.X * 9 + L.Y];
                e.ms += ms;
                e.launched += L.nq;
                e.computed += (long long)(c1 - c0);
                e.threaded = threaded;
                cudaEventDestroy(pe0);
                cudaEventDestroy(pe1);
            }
        }
        for (auto& s : m_streams)
            if (cudaStreamSynchronize(s) != cudaSuccess) return false;
        cudaEventRecord(m_e1);
        if (cudaEventSynchronize(m_e1) != cudaSuccess || cudaGetLastError() != cudaSuccess) return false;
        cudaEventElapsedTime(&kms, m_e0, m_e1);
        return true;
    }

    std::string name() const override { return m_name; }
    int n() const override { return m_n; }
    long long lastQuartets() const override { return m_last_quartets; }
    std::string stats() const override
    {
        char buf[512];
        std::snprintf(buf, sizeof buf,
                      "%d direct J/K builds, kernels %.1f ms, incl. transfers %.1f ms, quartet rows (host) %.1f ms, "
                      "launched %lld quartets, device memory %.1f MB",
                      m_builds, m_kernel_ms, m_total_ms, m_rows_ms, m_launched, m_bytes / 1.0e6);
        std::string out = buf;
        if (m_profile) {  // CURCUMA_QM_GPU_DIRECT_PROFILE=1: per pair class (serialised launches)
            static const char* lab = "spd";
            for (int k = 0; k < 81; ++k) {
                const ProfileEntry& e = m_prof[k];
                if (e.launched == 0) continue;
                const int X = k / 9, Y = k % 9;
                std::snprintf(buf, sizeof buf, "\n  class (%c%c|%c%c) %s: %.1f ms, launched %lld, computed %lld, %.1f computed/us",
                              lab[X / 3], lab[X % 3], lab[Y / 3], lab[Y % 3], e.threaded ? "thread" : (m_block_only ? "block" : "hermite"), e.ms,
                              e.launched, e.computed, e.ms > 0 ? e.computed / (1000.0 * e.ms) : 0.0);
                out += buf;
            }
        }
        return out;
    }

private:
    template <typename T>
    bool alloc(T** p, size_t count)
    {
        if (cudaMalloc(p, count * sizeof(T)) != cudaSuccess) { *p = nullptr; return false; }
        m_owned.push_back(*p);
        return true;
    }

    static constexpr int kStreams = 8;
    DevTables m_d{};
    std::vector<void*> m_owned;
    std::vector<long long> m_pc_off;
    std::vector<int> m_sorted;          // pair indices, per class sorted by Schwarz factor (desc)
    std::vector<double> m_sorted_q;     // their Schwarz factors
    int* m_sorted_dev = nullptr;
    std::vector<long long> m_row_off;   // ragged quartet rows of the current build (host)
    std::vector<int> m_row_bra, m_row_ket0;
    long long* m_rows_dev_off = nullptr;
    int *m_rows_dev_bra = nullptr, *m_rows_dev_ket0 = nullptr;
    size_t m_rows_cap = 0;
    unsigned long long* m_pall = nullptr;
    double m_rows_ms = 0.0;
    long long m_launched = 0;
    bool m_profile = false;
    bool m_block_only = false;
    struct ProfileEntry { double ms = 0.0; long long launched = 0, computed = 0; bool threaded = false; };
    ProfileEntry m_prof[81];
    int* m_width = nullptr;
    int* m_first = nullptr;
    double *m_P = nullptr, *m_J = nullptr, *m_K = nullptr, *m_Pmax = nullptr;
    unsigned long long* m_counter = nullptr;
    cudaStream_t m_streams[kStreams] = {};
    cudaEvent_t m_e0 = nullptr, m_e1 = nullptr;
    int m_n = 0, m_ns = 0, m_builds = 0;
    bool m_spherical = false, m_thread_ok = false;
    long long m_last_quartets = 0;
    double m_kernel_ms = 0.0, m_total_ms = 0.0;
    size_t m_bytes = 0;
    std::string m_name = "cuda";

public:
    void setBytes(size_t b) { m_bytes = b; }
};

} // namespace

// Integral-direct J/K builder (QMGpuJK*), or nullptr with a reason in `why`. Device memory:
// the shell/pair tables and 4 n^2 + nshell^2 doubles; no n^4 tensor. Claude Generated (Q-G3).
extern "C" void* curcuma_cuda_create_qm_direct_jk(const qmint::QMGpuEriTables* t, int device, char* why, int whylen)
{
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return nullptr; }
    cudaFree(nullptr);  // context creation
    size_t free0 = 0, total = 0;
    cudaMemGetInfo(&free0, &total);
    auto b = std::make_unique<QmGpuDirectJkCuda>();
    std::string err;
    if (!b->init(t, device, err)) {
        cudaGetLastError();  // clear a sticky-free error state for the CPU fallback
        setWhy(why, whylen, err);
        return nullptr;
    }
    size_t free1 = 0;
    cudaMemGetInfo(&free1, &total);
    b->setBytes(free0 > free1 ? free0 - free1 : 0);
    return static_cast<qmint::QMGpuJK*>(b.release());
}
