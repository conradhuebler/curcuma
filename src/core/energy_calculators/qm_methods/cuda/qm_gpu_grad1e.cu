/*
 * < CUDA one-electron part of the analytic RHF gradient of the native QM engine (Q-G8) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, Q-G8 of docs/QM_GPU_ROADMAP.md).
//
// Device counterpart of qmint::gradientOneElectron (qm_integrals.cpp), Eh/Bohr:
//
//   g_X = sum_mn P_mn d(T_mn + V_mn)/dX - sum_mn W_mn dS_mn/dX,
//   V_mn = - sum_C Z_C <m| 1/r_C |n>,   W = 2 sum_i^occ eps_i C_i C_i^T,
//
// (Pople, Krishnan, Schlegel, Binkley, Int. J. Quantum Chem. S13, 225 (1979); Helgaker,
// Jorgensen, Olsen, Molecular Electronic-Structure Theory, ch. 9). The nuclear repulsion gradient
// is not part of it (QMEngine adds it on the host, as for the CPU path).
//
// Summation scheme. The CPU loops over every ORDERED basis-function pair (m, n) and differentiates
// the first function only. Here every CANONICAL shell pair (A <= B) is visited once with the weight
// f = 2 (A != B) or 1 (A == B), and both function centres are differentiated:
//   f sum_{a in A, b in B} P_ab [ d/dA <a|O|b> -> atom(A),  d/dB <a|O|b> -> atom(B) ],
// which is the same sum because P and the operators are symmetric. Translational invariance then
// gives the rest:
//   S, T (two-centre): d/dB = -d/dA, so only d/dA is evaluated (atom(B) gets the negative); a pair
//         whose two shells sit on one atom contributes nothing and is skipped;
//   V (operator part, nucleus C moves): d/dC <a|1/r_C|b> = -(d/dA + d/dB) <a|1/r_C|b>, exactly the
//         Hellmann-Feynman expression of the CPU; a pair on the atom of C itself is skipped (its
//         basis and operator parts cancel).
//
// Derivative integrals (centre derivative of a cartesian Gaussian, as everywhere in the WP8 block):
//   d/dA_x [x_A^i e^{-a x_A^2}] = 2a x_A^{i+1} e^{..} - i x_A^{i-1} e^{..}.
// The Hermite tables E^{ij}_t of the Q-G7 gradient tables (qmint::buildEriGradGpuTables, built one
// angular-momentum step higher on both shells, K_AB folded into the x table) supply everything:
//   overlap      S_ab = (pi/p)^{3/2} s_x s_y s_z,  s(i,j) = E^{ij}_0 (the 1D overlap);
//   kinetic      T_ab = 1/2 sum_k <d_k a|d_k b>, per direction
//                t(i,j) = 4ab s(i+1,j+1) - 2a j s(i+1,j-1) - 2b i s(i-1,j+1) + i j s(i-1,j-1)
//                (the identity of primitiveKinetic); its A derivative needs s(LA+2, j), one step
//                beyond the table, from the McMurchie-Davidson recursion at t = 0:
//                E^{i+1,j}_0 = X_PA E^{ij}_0 + E^{ij}_1;
//   attraction   <a|1/r_C|b> = sum_tuv E_t E_u E_v R^0_tuv(p, P - C),  R^n_000 = (2 pi/p) (-2p)^n F_n(p|PC|^2)
//                (Helgaker 9.9.10-9.9.14, the recursion of buildRblock: z before y before x); the
//                A/B derivatives replace one 1D coefficient by 2a E^{i+1,j}_t - i E^{i-1,j}_t
//                (resp. 2b E^{i,j+1}_t - j E^{i,j-1}_t), Hermite order up to LA+LB+1.
// Boys function: boysDev, the grid + Taylor table of the CPU (qm_gpu_eri_kernels.cuh, unchanged).
//
// Work decomposition: one thread per (canonical shell pair, slot), slot = nucleus C (attraction)
// or one extra slot per pair for S and T; the thread loops over the primitive pairs and the
// cartesian component pairs and contracts with P (and W) at once, so only 3 natoms numbers leave
// the device. Pairs are taken in pair-class order (shells sorted by L), so a warp mostly works on
// one class. Accumulation per block in shared memory (FP64 shared atomics) up to 400 atoms, then
// FP64 atomicAdd into the global gradient (last bits depend on the schedule). No screening (the CPU
// 1e gradient has none). Limits: s-g shells (stage 6b, Nov 2026: instance MAXL = 9 for f/g bases, MAXL = 5 for s/p/d),
// FP64 throughout.

#include "qm_gpu_eri_kernels.cuh"   // DevTables, boysDev, uploadCommonTables, upload, setWhy (unchanged)
#include "../qm_gpu_grad1e_backend.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

// Hermite order LA + LB + 1: 5 for two d shells (the s/p/d instance, code of before), 9 for two g shells (stage 6b,
// Nov 2026: the kernels are templates on MAXL; a basis with f/g shells runs the MAXL = 9 instance)
constexpr int kG1MaxL = 5;                 // LA + LB + 1 for two d shells
constexpr int kG1MaxLfg = 9;               // LA + LB + 1 for two g shells
__host__ __device__ constexpr int g1SimplexC(int L) { return L < 0 ? 0 : (L + 1) * (L + 2) * (L + 3) / 6; }
constexpr int kG1SharedMaxAtoms = 400;     // 9.6 KB of shared memory; above: global atomics directly
constexpr int kG1Block = 128;

__host__ __device__ inline int g1Simplex(int L) { return L < 0 ? 0 : (L + 1) * (L + 2) * (L + 3) / 6; }
/// Position of R_tuv in the simplex array (ordered by degree t+u+v, then t, then u).
__device__ inline int g1Idx(int t, int u, int v)
{
    const int d = t + u + v;
    return g1Simplex(d - 1) + t * (d + 1) - t * (t - 1) / 2 + u;
}

struct G1Args {
    const int* pairs;        // shell-pair indices in pair-class order (pc_pairs)
    long long npairs;
    const int* sh_atom;
    const double* sh_xyz;    // shell centres (Bohr)
    const double* atom_xyz;  // nuclei (Bohr)
    const double* atom_Z;
    const double *pp_ea, *pp_eb;
    const double *P, *W;     // cartesian density and energy-weighted density, n x n
    int natoms;
    double* grad;            // natoms x 3
    bool shared_grad;
};

__device__ inline int lOf(int ncomp) { return ncomp == 1 ? 0 : (ncomp == 3 ? 1 : 2); }
/// stage 6b: any L <= 4 (ncomp = 1, 3, 6, 10, 15)
__device__ inline int lOfFG(int ncomp) { return ncomp == 1 ? 0 : ncomp == 3 ? 1 : ncomp == 6 ? 2 : ncomp == 10 ? 3 : 4; }
template <int MAXL>
__device__ inline int lOfT(int ncomp) { if constexpr (MAXL <= kG1MaxL) return lOf(ncomp); else return lOfFG(ncomp); }

/// R^0_tuv (t+u+v <= L) of the attraction integral for one primitive pair and one nucleus, in place
/// in R (higher degrees first, as buildR<1> of the 2e gradient), base (2 pi / p) (-2p)^n F_n(T).
template <int MAXL>
__device__ void buildR1e(const DevTables& d, double p, double PCx, double PCy, double PCz, int L, double* R)
{
    double F[MAXL + 1], base[MAXL + 1];
    boysDev(L, p * (PCx * PCx + PCy * PCy + PCz * PCz), F, d);
    double fac = 2.0 * kPI / p;
    for (int n = 0; n <= L; ++n) {
        base[n] = fac * F[n];
        fac *= -2.0 * p;
    }
    R[0] = base[L];
    for (int lv = L - 1; lv >= 0; --lv) {
        for (int deg = L - lv; deg >= 1; --deg)          // reads degrees deg-1, deg-2 (still level lv+1)
            for (int t = 0; t <= deg; ++t)
                for (int u = 0; u <= deg - t; ++u) {
                    const int v = deg - t - u;
                    double r;
                    if (v > 0)
                        r = (v >= 2 ? (v - 1) * R[g1Idx(t, u, v - 2)] : 0.0) + PCz * R[g1Idx(t, u, v - 1)];
                    else if (u > 0)
                        r = (u >= 2 ? (u - 1) * R[g1Idx(t, u - 2, v)] : 0.0) + PCy * R[g1Idx(t, u - 1, v)];
                    else
                        r = (t >= 2 ? (t - 1) * R[g1Idx(t - 2, u, v)] : 0.0) + PCx * R[g1Idx(t - 1, u, v)];
                    R[g1Idx(t, u, v)] = r;
                }
        R[0] = base[lv];
    }
}

/// Attraction part of shell pair pi for nucleus C: gA += f sum P_ab d/dA <a|-Z_C/r_C|b>, gB likewise
/// for the B centre (the operator part -(gA + gB) goes to C at the caller).
template <int MAXL>
__device__ void attractionPair(const DevTables& d, const G1Args& g, int pi, int C, double gA[3], double gB[3])
{
    const int A = d.pr_A[pi], B = d.pr_B[pi];
    const int nA = d.sh_ncomp[A], nB = d.sh_ncomp[B];
    const int LA = lOfT<MAXL>(nA), LB = lOfT<MAXL>(nB), L = LA + LB + 1;
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B];
    const double *cA = d.coef + d.sh_coef_off[A], *cB = d.coef + d.sh_coef_off[B];
    const int *lA = d.lmn + d.sh_lmn_off[A], *lB = d.lmn + d.sh_lmn_off[B];
    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B];
    const size_t n = (size_t)d.n_out;
    const int sI = LB + 2, sT = (LA + 2) * (LB + 2), sz = (LA + LB + 3) * sT;
    const double f = (A != B ? 2.0 : 1.0) * (-g.atom_Z[C]);
    const double Cx = g.atom_xyz[3 * C], Cy = g.atom_xyz[3 * C + 1], Cz = g.atom_xyz[3 * C + 2];
    double R[g1SimplexC(MAXL)];
    for (int k = d.pr_prim_off[pi]; k < d.pr_prim_off[pi] + d.pr_nprim[pi]; ++k) {
        const double p = d.pp_p[k];
        const double a2 = 2.0 * g.pp_ea[k], b2 = 2.0 * g.pp_eb[k];
        const int ia = d.pp_ia[k], ib = d.pp_ib[k];
        const double* E = d.E + d.pp_E_off[k];
        buildR1e<MAXL>(d, p, d.pp_P[3 * k] - Cx, d.pp_P[3 * k + 1] - Cy, d.pp_P[3 * k + 2] - Cz, L, R);
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b) {
                const double w = f * g.P[(oA + a) * n + oB + b] * cA[a * npA + ia] * cB[b * npB + ib];
                if (w == 0.0) continue;
                const int* la = lA + 3 * a;
                const int* lb = lB + 3 * b;
                // 1D Hermite coefficients per direction, t <= i+j+1 (E^{ij}_t = 0 for t > i+j)
                double e[3][MAXL + 1], da[3][MAXL + 1], db[3][MAXL + 1];
                int nt[3], lab = 0;
                for (int x = 0; x < 3; ++x) {
                    const double* Ex = E + x * sz;
                    const int i = la[x], j = lb[x];
                    nt[x] = i + j + 1;
                    lab += i + j;
                    for (int t = 0; t <= nt[x]; ++t) {
                        const double* Et = Ex + t * sT;
                        e[x][t] = Et[i * sI + j];
                        da[x][t] = a2 * Et[(i + 1) * sI + j] - (i > 0 ? i * Et[(i - 1) * sI + j] : 0.0);
                        db[x][t] = b2 * Et[i * sI + j + 1] - (j > 0 ? j * Et[i * sI + j - 1] : 0.0);
                    }
                }
                double sA0 = 0, sA1 = 0, sA2 = 0, sB0 = 0, sB1 = 0, sB2 = 0;
                for (int t = 0; t <= nt[0]; ++t)
                    for (int u = 0; u <= nt[1] && t + u <= lab + 1; ++u)
                        for (int v = 0; v <= nt[2] && t + u + v <= lab + 1; ++v) {
                            // terms with two raised orders vanish (E^{ij}_t = 0 for t > i+j)
                            const double r = R[g1Idx(t, u, v)];
                            const double eyz = e[1][u] * e[2][v], exz = e[0][t] * e[2][v], exy = e[0][t] * e[1][u];
                            sA0 += da[0][t] * eyz * r;
                            sA1 += da[1][u] * exz * r;
                            sA2 += da[2][v] * exy * r;
                            sB0 += db[0][t] * eyz * r;
                            sB1 += db[1][u] * exz * r;
                            sB2 += db[2][v] * exy * r;
                        }
                gA[0] += w * sA0;
                gA[1] += w * sA1;
                gA[2] += w * sA2;
                gB[0] += w * sB0;
                gB[1] += w * sB1;
                gB[2] += w * sB2;
            }
    }
}

/// Overlap and kinetic part of shell pair pi: X = f sum_ab [P_ab dT_ab/dA - W_ab dS_ab/dA].
template <int MAXL>
__device__ void overlapKineticPair(const DevTables& d, const G1Args& g, int pi, double X[3])
{
    const int A = d.pr_A[pi], B = d.pr_B[pi];
    const int nA = d.sh_ncomp[A], nB = d.sh_ncomp[B];
    const int LA = lOfT<MAXL>(nA), LB = lOfT<MAXL>(nB);
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B];
    const double *cA = d.coef + d.sh_coef_off[A], *cB = d.coef + d.sh_coef_off[B];
    const int *lA = d.lmn + d.sh_lmn_off[A], *lB = d.lmn + d.sh_lmn_off[B];
    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B];
    const size_t n = (size_t)d.n_out;
    const int sI = LB + 2, sT = (LA + 2) * (LB + 2), sz = (LA + LB + 3) * sT;
    const double f = (A != B ? 2.0 : 1.0);
    for (int k = d.pr_prim_off[pi]; k < d.pr_prim_off[pi] + d.pr_nprim[pi]; ++k) {
        const double p = d.pp_p[k];
        const double al = g.pp_ea[k], be = g.pp_eb[k];
        const int ia = d.pp_ia[k], ib = d.pp_ib[k];
        const double* E = d.E + d.pp_E_off[k];
        const double pref = (kPI / p) * sqrt(kPI / p);
        double PA[3];
        for (int x = 0; x < 3; ++x) PA[x] = d.pp_P[3 * k + x] - g.sh_xyz[3 * A + x];
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b) {
                const double cc = cA[a * npA + ia] * cB[b * npB + ib];
                const double wP = f * cc * pref * g.P[(oA + a) * n + oB + b];
                const double wW = f * cc * pref * g.W[(oA + a) * n + oB + b];
                if (wP == 0.0 && wW == 0.0) continue;
                const int* la = lA + 3 * a;
                const int* lb = lB + 3 * b;
                double s[3], ds[3], tk[3], dt[3];
                for (int x = 0; x < 3; ++x) {
                    const double* Ex = E + x * sz;
                    // 1D overlap s(i,j) = E^{ij}_0; i = LA+2 from E^{i,j}_0 = X_PA E^{i-1,j}_0 + E^{i-1,j}_1
                    auto s1 = [&](int i, int j) -> double {
                        if (i < 0 || j < 0) return 0.0;
                        if (i <= LA + 1) return Ex[i * sI + j];
                        return PA[x] * Ex[(i - 1) * sI + j] + Ex[sT + (i - 1) * sI + j];
                    };
                    // 1D kinetic-type <d a|d b> (identity of primitiveKinetic)
                    auto t1 = [&](int i, int j) -> double {
                        if (i < 0) return 0.0;
                        return 4.0 * al * be * s1(i + 1, j + 1) - 2.0 * al * j * s1(i + 1, j - 1)
                             - 2.0 * be * i * s1(i - 1, j + 1) + (double)(i * j) * s1(i - 1, j - 1);
                    };
                    const int i = la[x], j = lb[x];
                    s[x] = s1(i, j);
                    ds[x] = 2.0 * al * s1(i + 1, j) - i * s1(i - 1, j);
                    tk[x] = t1(i, j);
                    dt[x] = 2.0 * al * t1(i + 1, j) - i * t1(i - 1, j);
                }
                const double dS[3] = { ds[0] * s[1] * s[2], s[0] * ds[1] * s[2], s[0] * s[1] * ds[2] };
                const double dT[3] = { 0.5 * (dt[0] * s[1] * s[2] + ds[0] * (tk[1] * s[2] + s[1] * tk[2])),
                                       0.5 * (dt[1] * s[0] * s[2] + ds[1] * (tk[0] * s[2] + s[0] * tk[2])),
                                       0.5 * (dt[2] * s[0] * s[1] + ds[2] * (tk[0] * s[1] + s[0] * tk[1])) };
                for (int x = 0; x < 3; ++x) X[x] += wP * dT[x] - wW * dS[x];
            }
    }
}

/// One thread per (shell pair, slot); slot < natoms: attraction of nucleus slot, slot == natoms: S, T.
template <int MAXL>
__global__ void __launch_bounds__(kG1Block) grad1eKernel(DevTables d, G1Args g)
{
    extern __shared__ double sgrad[];
    if (g.shared_grad) {
        for (int k = threadIdx.x; k < 3 * g.natoms; k += blockDim.x) sgrad[k] = 0.0;
        __syncthreads();
    }
    double* dst = g.shared_grad ? sgrad : g.grad;
    const long long nslot = g.natoms + 1, total = g.npairs * nslot;
    for (long long e = (long long)blockIdx.x * blockDim.x + threadIdx.x; e < total; e += (long long)gridDim.x * blockDim.x) {
        const int pi = g.pairs[e / nslot];
        const int slot = (int)(e % nslot);
        const int aA = g.sh_atom[d.pr_A[pi]], aB = g.sh_atom[d.pr_B[pi]];
        if (slot == g.natoms) {
            if (aA == aB) continue;  // two-centre terms on one atom: no contribution
            double X[3] = { 0.0, 0.0, 0.0 };
            overlapKineticPair<MAXL>(d, g, pi, X);
            for (int x = 0; x < 3; ++x) {
                atomicAdd(&dst[3 * aA + x], X[x]);
                atomicAdd(&dst[3 * aB + x], -X[x]);
            }
        } else {
            const int C = slot;
            if (aA == aB && aA == C) continue;  // basis and operator parts cancel
            double gA[3] = { 0.0, 0.0, 0.0 }, gB[3] = { 0.0, 0.0, 0.0 };
            attractionPair<MAXL>(d, g, pi, C, gA, gB);
            for (int x = 0; x < 3; ++x) {
                atomicAdd(&dst[3 * aA + x], gA[x]);
                atomicAdd(&dst[3 * aB + x], gB[x]);
                atomicAdd(&dst[3 * C + x], -(gA[x] + gB[x]));
            }
        }
    }
    if (g.shared_grad) {
        __syncthreads();
        for (int k = threadIdx.x; k < 3 * g.natoms; k += blockDim.x)
            if (sgrad[k] != 0.0) atomicAdd(&g.grad[k], sgrad[k]);
    }
}

class QmGpuGrad1eCuda final : public qmint::QMGpuGrad1e {
public:
    ~QmGpuGrad1eCuda() override
    {
        for (void* p : m_owned) cudaFree(p);
        if (m_e0) cudaEventDestroy(m_e0);
        if (m_e1) cudaEventDestroy(m_e1);
    }

    bool init(const qmint::QMGpuEriTables* t, const qmint::QMGpuGrad1eExtra* x, int device, std::string& why)
    {
        if (!t->pc_off || !t->pc_pairs) { why = "tables without pair classes"; return false; }
        if (t->n_T != 0) { why = "gradient tables must be cartesian (no 6d->5d transform)"; return false; }
        if (!x || x->base.nprimpair != t->nprimpair || x->base.natoms <= 0 || !x->sh_xyz || !x->atom_xyz || !x->atom_Z) {
            why = "1e gradient extra data inconsistent";
            return false;
        }
        if (!uploadCommonTables(t, m_d, m_owned, why, 4)) return false;  // also refuses L > 4 (stage 6b: f/g)
        int maxL = 0;
        for (int s = 0; s < t->nshell; ++s) maxL = std::max(maxL, t->sh_L[s]);
        m_fg = maxL > 2;
        m_n = t->n_out;
        m_natoms = x->base.natoms;
        m_npairs = t->npair;
        const size_t n2 = (size_t)m_n * m_n;
        const bool ok = upload(&m_pairs, t->pc_pairs, (size_t)t->npair, m_owned)
            && upload(&m_sh_atom, x->base.sh_atom, (size_t)t->nshell, m_owned)
            && upload(&m_pp_ea, x->base.pp_ea, (size_t)t->nprimpair, m_owned)
            && upload(&m_pp_eb, x->base.pp_eb, (size_t)t->nprimpair, m_owned)
            && upload(&m_sh_xyz, x->sh_xyz, 3 * (size_t)t->nshell, m_owned)
            && upload(&m_atom_xyz, x->atom_xyz, 3 * (size_t)m_natoms, m_owned)
            && upload(&m_atom_Z, x->atom_Z, (size_t)m_natoms, m_owned)
            && alloc(&m_P, n2) && alloc(&m_W, n2) && alloc(&m_grad, 3 * (size_t)m_natoms);
        if (!ok) { why = "device allocation / upload for the 1e gradient failed"; return false; }
        cudaEventCreate(&m_e0);
        cudaEventCreate(&m_e1);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        m_sms = prop.multiProcessorCount > 0 ? prop.multiProcessorCount : 16;
        return cudaGetLastError() == cudaSuccess;
    }

    bool compute(const double* P, const double* W, double* grad) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n2 = (size_t)m_n * m_n;
        if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess
            || cudaMemcpy(m_W, W, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
            return false;
        cudaMemset(m_grad, 0, 3 * (size_t)m_natoms * sizeof(double));
        G1Args g{};
        g.pairs = m_pairs;
        g.npairs = m_npairs;
        g.sh_atom = m_sh_atom;
        g.sh_xyz = m_sh_xyz;
        g.atom_xyz = m_atom_xyz;
        g.atom_Z = m_atom_Z;
        g.pp_ea = m_pp_ea;
        g.pp_eb = m_pp_eb;
        g.P = m_P;
        g.W = m_W;
        g.natoms = m_natoms;
        g.grad = m_grad;
        g.shared_grad = m_natoms <= kG1SharedMaxAtoms;
        const long long total = (long long)m_npairs * (m_natoms + 1);
        const long long nb = std::min<long long>((total + kG1Block - 1) / kG1Block, 16LL * m_sms);
        const size_t shm = g.shared_grad ? 3 * (size_t)m_natoms * sizeof(double) : 0;
        cudaEventRecord(m_e0);
        if (m_fg)  // stage 6b: Hermite order up to 9
            grad1eKernel<kG1MaxLfg><<<(unsigned)std::max<long long>(nb, 1), kG1Block, shm>>>(m_d, g);
        else
            grad1eKernel<kG1MaxL><<<(unsigned)std::max<long long>(nb, 1), kG1Block, shm>>>(m_d, g);
        if (cudaGetLastError() != cudaSuccess) return false;
        cudaEventRecord(m_e1);
        if (cudaEventSynchronize(m_e1) != cudaSuccess || cudaGetLastError() != cudaSuccess) return false;
        float kms = 0.0f;
        cudaEventElapsedTime(&kms, m_e0, m_e1);
        if (cudaMemcpy(grad, m_grad, 3 * (size_t)m_natoms * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        m_kernel_ms = kms;
        m_total_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    std::string name() const override { return m_name; }
    std::string stats() const override
    {
        char buf[200];
        std::snprintf(buf, sizeof buf, "%d shell pairs x %d slots, kernel %.2f ms, incl. transfers %.2f ms, device memory %.1f MB",
                      m_npairs, m_natoms + 1, m_kernel_ms, m_total_ms, m_bytes / 1.0e6);
        return buf;
    }
    void setBytes(size_t b) { m_bytes = b; }

private:
    template <typename T>
    bool alloc(T** p, size_t count)
    {
        if (cudaMalloc(p, count * sizeof(T)) != cudaSuccess) { *p = nullptr; return false; }
        m_owned.push_back(*p);
        return true;
    }

    DevTables m_d{};
    std::vector<void*> m_owned;
    int *m_pairs = nullptr, *m_sh_atom = nullptr;
    double *m_pp_ea = nullptr, *m_pp_eb = nullptr, *m_sh_xyz = nullptr, *m_atom_xyz = nullptr, *m_atom_Z = nullptr;
    double *m_P = nullptr, *m_W = nullptr, *m_grad = nullptr;
    cudaEvent_t m_e0 = nullptr, m_e1 = nullptr;
    int m_n = 0, m_natoms = 0, m_npairs = 0, m_sms = 16;
    bool m_fg = false;  // stage 6b: basis with f/g shells
    double m_kernel_ms = 0.0, m_total_ms = 0.0;
    size_t m_bytes = 0;
    std::string m_name = "cuda";
};

} // namespace

// 1e-gradient builder (QMGpuGrad1e*), or nullptr with a reason in `why`. `t` from
// qmint::buildEriGradGpuTables (cartesian, raised Hermite tables), the same tables the 2e gradient
// (Q-G7) reads. Claude Generated (Oct 2026, Q-G8).
extern "C" void* curcuma_cuda_create_qm_grad1e(const qmint::QMGpuEriTables* t, const qmint::QMGpuGrad1eExtra* x, int device,
                                               char* why, int whylen)
{
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return nullptr; }
    cudaFree(nullptr);  // context creation
    size_t free0 = 0, total = 0;
    cudaMemGetInfo(&free0, &total);
    auto b = std::make_unique<QmGpuGrad1eCuda>();
    std::string err;
    if (!b->init(t, x, device, err)) {
        cudaGetLastError();
        setWhy(why, whylen, err);
        return nullptr;
    }
    size_t free1 = 0;
    cudaMemGetInfo(&free1, &total);
    b->setBytes(free0 > free1 ? free0 - free1 : 0);
    return static_cast<qmint::QMGpuGrad1e*>(b.release());
}
