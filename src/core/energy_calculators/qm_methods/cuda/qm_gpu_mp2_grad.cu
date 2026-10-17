/*
 * < CUDA two-electron derivative-integral contraction of the analytic MP2 gradient >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, MP2 gradient). Device counterpart of qmint::gradientTwoElectronMP2
// (qm_integrals.cpp, declared in qm_mp2.h):
//
//   g_X = 1/2 sum_{mnls} D_mnls d(mn|ls)/dX,
//   D_mnls = P_mn P_ls - 1/4 (P_ml P_ns + P_ms P_nl)                                         (RHF)
//          + P2_mn P_ls + P_mn P2_ls - 1/4 (P2_ml P_ns + P_ml P2_ns + P2_ms P_nl + P_ms P2_nl)  (P2 x P)
//          + G_mnls                                                     (non-separable MP2 2PDM)
//
// with P, P2 and G in the CARTESIAN basis (G: ncart^4 doubles on the device, 8-fold symmetric; the host
// takes the active-basis 2PDM to the cartesian one). Every quartet is evaluated (no screening), as on the CPU.
//
// The kernel machinery below (Hermite maps, R^0 recursion, orientation derivatives, warp kernel, class
// launches) is a COPY of cuda/qm_gpu_grad.cu (Q-G7), kept in a separate translation unit so that the
// validated Q-G7 file is not touched; only step 1 (the quartet weights D) and the argument block differ.
// A change of the derivative scheme in one file must go into the other (invariant in qm_methods/CLAUDE.md).
//
// Original Q-G7 header, for the scheme:
// Device counterpart of qmint::gradientTwoElectron (qm_integrals.cpp):
//
//   g_X = 1/2 sum_{mnls} D_mnls d(mn|ls)/dX,   D_mnls = P_mn P_ls - 1/4 (P_ml P_ns + P_ms P_nl)
//
// over canonical shell quartets (pair i <= pair j) only, each weighted with
// f = 1/2 (A!=B ? 2:1)(C!=D ? 2:1)(i!=j ? 2:1); per quartet the derivatives on the centres of A, B
// and C are evaluated, the one on D follows from translational invariance d/dD = -(d/dA+d/dB+d/dC);
// quartets whose four shells sit on one atom contribute nothing and are skipped. All of this is the
// CPU scheme. Integral-direct: every derivative integral is contracted with D at once, nothing of
// size n^4 exists; device memory is P and Pmax (O(n^2)) plus the pair tables.
//
// Derivative integrals (Helgaker, Jorgensen, Olsen ch. 9; McMurchie, Davidson, J. Comput. Phys. 26,
// 218 (1978)): the centre derivative of a cartesian Gaussian is
//   d/dA_x [x_A^l e^{-a r_A^2}] = 2a x_A^{l+1} e^{..} - l x_A^{l-1} e^{..},
// so with Hermite tables built one step higher on both shells of a pair (host:
// qmint::buildEriGradGpuTables) a derivative is an ordinary McMurchie-Davidson sum with a modified
// one-dimensional bra Hermite coefficient, e.g. for d/dA_x
//   E'_t = 2 alpha E_t^{i+1,j} - i E_t^{i-1,j}      (t <= i+j+1, y and z unchanged).
//
// Contraction order (the density is folded in as early as possible, cf. the ket-first contraction of
// the CPU ERI kernel shellQuartet):
//   1. per shell quartet: D_abcd for all component quartets;
//   2. per ket primitive pair: Dket_ab(tau,nu,phi) = sum_cd D_abcd c_c c_d (-1)^(tau+nu+phi)
//      E^cd_tau E^cd_nu E^cd_phi   (ket Hermite simplex tau+nu+phi <= LC+LD);
//   3. per primitive quartet: R^0 on the simplex t+u+v <= L = LA+LB+LC+LD+1, and for every bra
//      component pair ab and bra Hermite index (t,u,v) with t+u+v <= LA+LB+1
//        K_ab(t,u,v) = sum_{tau nu phi} Dket_ab(tau,nu,phi) R_{t+tau,u+nu,v+phi},
//      contracted with c_a c_b and the modified bra coefficients into dE/dA and dE/dB;
//   4. the same routine with bra and ket swapped, (CD|AB), gives dE/dC.
// The CPU builds h(c', d) blocks for every reachable ket power instead; same sums, regrouped.
//
// Work decomposition: one launch per pair-class combination (X, Y) over the rectangle
// pairs(X) x pairs(Y), the non-canonical half returns at once (as Q-G3), grid-stride loops, one warp
// per shell quartet (D, Dket and the R^n level buffers in shared memory, warp-shuffle reduction).
// Per orientation the warp works in one of two modes:
//   entry mode      primitive quartets in sequence, lanes over the Dket and the bra entries;
//   primitive mode  when the bra has fewer than 32 Hermite entries (s and sp-type bras): Dket of
//                   every ket primitive pair at once, then lanes over the primitive quartets, each
//                   lane with its own R^0 in thread-local memory (benzene 233 -> 173 ms kernels).
// (A thread-per-quartet kernel for the small classes was slower than both and was removed.)
// Accumulation per block in shared memory (FP64 shared atomics) up to 400 atoms, then FP64 atomicAdd
// into the global gradient; the last bits depend on the schedule. Screening: density-weighted
// shell-quartet screen with the threshold passed by the caller (0 = every quartet, i.e. exactly the
// CPU quartet set); no primitive screening (the CPU gradient has none). Limits: s, p, d shells
// (L <= 2; the derivatives reach f-type orders through the raised tables). Diagnostic:
// CURCUMA_QM_GPU_GRAD_PROFILE=1 serialises the class launches and appends per-class times to stats().

#include "qm_gpu_eri_kernels.cuh"   // DevTables, boysDev, uploadCommonTables, upload, setWhy (unchanged)
#include "../qm_gpu_mp2_grad_backend.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

namespace {

__host__ __device__ inline int nsimplex(int L) { return L < 0 ? 0 : (L + 1) * (L + 2) * (L + 3) / 6; }

// primitive mode (lanes over primitive quartets): thread-local R^0 simplex, L <= 5
constexpr int kPrimNS = 56;
constexpr int kMaxWarpsPerBlock = 4;
constexpr int kSharedGradMaxAtoms = 400;          // 9.6 KB; above: global atomics directly
constexpr size_t kSharedBudget = 46 * 1024;       // stay below the 48 KB default (no opt-in)

struct GradArgs {
    const int* sh_atom;
    const double *pp_ea, *pp_eb;  // primitive exponents of the two shells of a primitive pair
    const double* P;              // cartesian RHF density, n x n
    const double* P2;             // cartesian MP2 relaxed density correction, n x n
    const double* G;              // cartesian non-separable 2PDM, n^4
    const double* Pmax;           // unused (no screening), kept for the copied quartet test
    double thr;                   // 0: every quartet
    int natoms;
    double* grad;                 // natoms x 3
    unsigned long long* counter;
    const int *rb, *rk;           // pair lists of the two classes
    long long nk, nq;             // rectangle: nq = |rb| * nk
    int LA, LB, LC, LD, L;        // L = LA+LB+LC+LD+1
    int NS, ndq, ndket, per_warp; // simplex(L), D block, Dket, doubles per warp
    bool prim1, prim2;            // orientation (AB|CD) / (CD|AB) in primitive mode
    bool shared_grad;
};

template <int GS>
__device__ inline void gsync()
{
    if constexpr (GS > 1) __syncwarp();
}

/// Hermite index tables for degree <= L, built once per block. M[k] = (t,u,v) packed t|u<<8|v<<16,
/// ordered by degree t+u+v (so the simplex of degree <= m is the prefix of length nsimplex(m));
/// map[(t*s+u)*s+v] = k (s = L+1).
__device__ void buildHermiteMaps(int L, short* map, int* M)
{
    const int s = L + 1, s3 = s * s * s;
    for (int idx = threadIdx.x; idx < s3; idx += blockDim.x) {
        const int t = idx / (s * s), u = (idx / s) % s, v = idx % s;
        const int deg = t + u + v;
        if (deg > L) { map[idx] = -1; continue; }
        const int k = nsimplex(deg - 1) + t * (deg + 1) - t * (t - 1) / 2 + u;
        map[idx] = (short)k;
        M[k] = t | (u << 8) | (v << 16);
    }
}

/// Shell quartet of rectangle element q, or false (non-canonical half, one-atom quartet, screened).
__device__ inline bool gradQuartet(const DevTables& d, const GradArgs& g, long long q, int& i, int& j)
{
    i = g.rb[q / g.nk];
    j = g.rk[q % g.nk];
    if (i > j) return false;
    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const int aA = g.sh_atom[A];
    if (aA == g.sh_atom[B] && aA == g.sh_atom[C] && aA == g.sh_atom[D]) return false;
    if (g.thr > 0.0) {
        const int ns = d.nshell;
        const double bound = g.Pmax[A * ns + B] * g.Pmax[C * ns + D]
            + 0.25 * (g.Pmax[A * ns + C] * g.Pmax[B * ns + D] + g.Pmax[A * ns + D] * g.Pmax[B * ns + C]);
        if (d.pr_schwarz[i] * d.pr_schwarz[j] * bound < g.thr) return false;
    }
    return true;
}

/// R^0_{tuv} (t+u+v <= L, stored at map positions) of one primitive quartet (bra pair pb, ket pair
/// pk), same recursion as buildR0Simplex (base (-2 rho)^n F_n(T), z before y before x). GS = 1:
/// in place in Rw (higher degrees first); GS = 32: two level buffers Rw, Rw + NS, result in Rw.
/// sF: F[0..L], (-2 rho)^n at sF + 10, W at sF + 20..22, prefactor at sF + 23.
template <int GS>
__device__ void buildR(const DevTables& d, int pb, int pk, int L, int NS, int lane, const int* M, const short* map,
                       double* Rw, double* sF)
{
    const int s = L + 1;
    gsync<GS>();  // every lane is done with the previous R and sF
    if (lane == 0) {
        const double p = d.pp_p[pb], q = d.pp_p[pk];
        const double rho = p * q / (p + q);
        const double Wx = d.pp_P[3 * pb] - d.pp_P[3 * pk];
        const double Wy = d.pp_P[3 * pb + 1] - d.pp_P[3 * pk + 1];
        const double Wz = d.pp_P[3 * pb + 2] - d.pp_P[3 * pk + 2];
        boysDev(L, rho * (Wx * Wx + Wy * Wy + Wz * Wz), sF, d);
        double* bp = sF + 10;
        bp[0] = 1.0;
        for (int m = 1; m <= L; ++m) bp[m] = bp[m - 1] * (-2.0 * rho);
        sF[20] = Wx;
        sF[21] = Wy;
        sF[22] = Wz;
        sF[23] = d.pref / (p * q * sqrt(p + q));
        if constexpr (GS > 1) Rw[(L & 1) * NS] = bp[L] * sF[L];
        else Rw[0] = bp[L] * sF[L];
    }
    gsync<GS>();
    const double Wx = sF[20], Wy = sF[21], Wz = sF[22];
    auto rec = [&](const double* c, int k) {
        const int tv = M[k];
        const int t = tv & 255, u = (tv >> 8) & 255, v = tv >> 16;
        if (v > 0) return (v >= 2 ? (v - 1) * c[map[(t * s + u) * s + v - 2]] : 0.0) + Wz * c[map[(t * s + u) * s + v - 1]];
        if (u > 0) return (u >= 2 ? (u - 1) * c[map[(t * s + u - 2) * s + v]] : 0.0) + Wy * c[map[(t * s + u - 1) * s + v]];
        return (t >= 2 ? (t - 1) * c[map[((t - 2) * s + u) * s + v]] : 0.0) + Wx * c[map[((t - 1) * s + u) * s + v]];
    };
    if constexpr (GS == 1) {
        for (int lv = L - 1; lv >= 0; --lv) {
            for (int k = nsimplex(L - lv) - 1; k >= 1; --k) Rw[k] = rec(Rw, k);  // reads lower degrees (level lv+1)
            Rw[0] = sF[10 + lv] * sF[lv];
        }
    } else {
        for (int lv = L - 1; lv >= 0; --lv) {
            const double* c = Rw + ((lv + 1) & 1) * NS;
            double* x = Rw + (lv & 1) * NS;
            const int top = nsimplex(L - lv);
            for (int k = lane; k < top; k += GS) x[k] = (k == 0) ? sF[10 + lv] * sF[lv] : rec(c, k);
            gsync<GS>();
        }
    }
}

/// Bra-centre derivatives of one orientation (bra pair pi, ket pair pj), steps 2 and 3 of the
/// header: acc[0..2] += sum D d(ab|cd)/dA_k, and with BOTH acc[3..5] += ... /dB_k (A, B = shells of
/// pi). Dq holds D for (AB|CD) of the quartet as [ab][cd]; `transposed` = this orientation is (CD|AB).
template <int GS, bool BOTH>
__device__ void orientationGrad(const DevTables& d, const GradArgs& g, int pi, int pj, const double* Dq, bool transposed,
                                bool prim_mode, const int* M, const short* map, double* Rw, double* sF, double* Dket,
                                int lane, double acc[6])
{
    const int A = d.pr_A[pi], B = d.pr_B[pi], C = d.pr_A[pj], D = d.pr_B[pj];
    const int La = transposed ? g.LC : g.LA, Lb = transposed ? g.LD : g.LB;
    const int Lc = transposed ? g.LA : g.LC, Ld = transposed ? g.LB : g.LD;
    const int nB = d.sh_ncomp[B], nD = d.sh_ncomp[D];
    const int nBra = d.sh_ncomp[A] * nB, nKet = d.sh_ncomp[C] * nD;
    const int npA = d.sh_nprim[A], npB = d.sh_nprim[B], npC = d.sh_nprim[C], npD = d.sh_nprim[D];
    const double *cA = d.coef + d.sh_coef_off[A], *cB = d.coef + d.sh_coef_off[B];
    const double *cC = d.coef + d.sh_coef_off[C], *cD = d.coef + d.sh_coef_off[D];
    const int *lA = d.lmn + d.sh_lmn_off[A], *lB = d.lmn + d.sh_lmn_off[B];
    const int *lC = d.lmn + d.sh_lmn_off[C], *lD = d.lmn + d.sh_lmn_off[D];
    // raised Hermite tables: E[t][i][j], i <= L+1, j <= L+1, t <= Li+Lj+2
    const int sAB_t = (La + 2) * (Lb + 2), sAB_i = Lb + 2, szAB = (La + Lb + 3) * sAB_t;
    const int sCD_t = (Lc + 2) * (Ld + 2), sCD_i = Ld + 2, szCD = (Lc + Ld + 3) * sCD_t;
    const int L = g.L, s = L + 1;
    const int NSK = nsimplex(Lc + Ld), NSB1 = nsimplex(La + Lb + 1);
    const int pb0 = d.pr_prim_off[pi], nbp = d.pr_nprim[pi];
    const int pk0 = d.pr_prim_off[pj], nkp = d.pr_nprim[pj];

    // step 2 for one entry: Dket_ab(h) = sum_cd D_abcd c_c c_d (-1)^|h| E^cd_tau E^cd_nu E^cd_phi
    auto dketEntry = [&](int pk, int ab, int k) {
        const int ic = d.pp_ia[pk], id = d.pp_ib[pk];
        const double* Ecd = d.E + d.pp_E_off[pk];
        const int tv = M[k];
        const int tau = tv & 255, nu = (tv >> 8) & 255, ph = tv >> 16;
        double sum = 0.0;
        for (int cd = 0; cd < nKet; ++cd) {
            const int c = cd / nD, dd = cd % nD;
            const int* lc = lC + 3 * c;
            const int* ld = lD + 3 * dd;
            if (tau > lc[0] + ld[0] || nu > lc[1] + ld[1] || ph > lc[2] + ld[2]) continue;
            const double dv = transposed ? Dq[cd * nBra + ab] : Dq[ab * nKet + cd];
            sum += dv * cC[c * npC + ic] * cD[dd * npD + id] * Ecd[tau * sCD_t + lc[0] * sCD_i + ld[0]]
                 * Ecd[szCD + nu * sCD_t + lc[1] * sCD_i + ld[1]] * Ecd[2 * szCD + ph * sCD_t + lc[2] * sCD_i + ld[2]];
        }
        return ((tau + nu + ph) & 1) ? -sum : sum;
    };
    // step 3 for one bra entry (ab, (t,u,v) = M[k1]) of one primitive quartet
    auto braEntry = [&](int pb, const double* R, double pref, const double* Dk, int ab, int k1) {
        const double alpha2 = 2.0 * g.pp_ea[pb], beta2 = 2.0 * g.pp_eb[pb];
        const int ia = d.pp_ia[pb], ibb = d.pp_ib[pb];
        const double* Eab = d.E + d.pp_E_off[pb];
        const int a = ab / nB, b = ab % nB;
        const int tv = M[k1];
        const int tuv[3] = { tv & 255, (tv >> 8) & 255, tv >> 16 };
        const int* la = lA + 3 * a;
        const int* lb = lB + 3 * b;
        double ev[3], dA[3], dB[3];
#pragma unroll
        for (int x = 0; x < 3; ++x) {
            const double* Ex = Eab + x * szAB + tuv[x] * sAB_t;
            const int i = la[x], j = lb[x];
            ev[x] = Ex[i * sAB_i + j];  // zero for t > i+j (forward recursion)
            dA[x] = alpha2 * Ex[(i + 1) * sAB_i + j] - (i > 0 ? i * Ex[(i - 1) * sAB_i + j] : 0.0);
            dB[x] = BOTH ? beta2 * Ex[i * sAB_i + j + 1] - (j > 0 ? j * Ex[i * sAB_i + j - 1] : 0.0) : 0.0;
        }
        const double pA0 = dA[0] * ev[1] * ev[2], pA1 = ev[0] * dA[1] * ev[2], pA2 = ev[0] * ev[1] * dA[2];
        const double pB0 = dB[0] * ev[1] * ev[2], pB1 = ev[0] * dB[1] * ev[2], pB2 = ev[0] * ev[1] * dB[2];
        if (pA0 == 0.0 && pA1 == 0.0 && pA2 == 0.0 && pB0 == 0.0 && pB1 == 0.0 && pB2 == 0.0) return;
        double K = 0.0;
        for (int k = 0; k < NSK; ++k) {
            const int hv = M[k];
            K += Dk[k] * R[map[((tuv[0] + (hv & 255)) * s + tuv[1] + ((hv >> 8) & 255)) * s + tuv[2] + (hv >> 16)]];
        }
        const double w = pref * cA[a * npA + ia] * cB[b * npB + ibb] * K;
        acc[0] += w * pA0;
        acc[1] += w * pA1;
        acc[2] += w * pA2;
        if (BOTH) {
            acc[3] += w * pB0;
            acc[4] += w * pB1;
            acc[5] += w * pB2;
        }
    };

    if (prim_mode) {
        // Primitive mode (few bra entries): Dket of every ket primitive pair at once, then lanes over
        // the primitive quartets, each with its own R^0 in thread-local memory.
        const int per = nBra * NSK;
        gsync<GS>();  // previous use of Dket is finished
        for (int e = lane; e < nkp * per; e += GS) Dket[e] = dketEntry(pk0 + e / per, (e % per) / NSK, e % NSK);
        gsync<GS>();
        double Rl[kPrimNS], sFl[24];
        for (int e = lane; e < nbp * nkp; e += GS) {
            const int pb = pb0 + e / nkp, ik = e % nkp;
            buildR<1>(d, pb, pk0 + ik, L, g.NS, 0, M, map, Rl, sFl);
            for (int ab = 0; ab < nBra; ++ab)
                for (int k1 = 0; k1 < NSB1; ++k1) braEntry(pb, Rl, sFl[23], Dket + ik * per + ab * NSK, ab, k1);
        }
        return;
    }
    // Entry mode: primitive quartets in sequence, lanes over the Dket and bra entries.
    for (int ik = 0; ik < nkp; ++ik) {
        const int pk = pk0 + ik;
        gsync<GS>();  // previous K loop is done with Dket
        for (int e = lane; e < nBra * NSK; e += GS) Dket[e] = dketEntry(pk, e / NSK, e % NSK);
        gsync<GS>();
        for (int ib = 0; ib < nbp; ++ib) {
            const int pb = pb0 + ib;
            buildR<GS>(d, pb, pk, L, g.NS, lane, M, map, Rw, sF);
            for (int e = lane; e < nBra * NSB1; e += GS) braEntry(pb, Rw, sF[23], Dket + (e / NSB1) * NSK, e / NSB1, e % NSB1);
        }
    }
}

/// d/dA, d/dB, d/dC (x,y,z each) of one canonical shell quartet, summed over its components with the
/// density weights (before the degeneracy factor). Lane-partial sums for GS = 32.
template <int GS>
__device__ void quartetGrad(const DevTables& d, const GradArgs& g, int i, int j, const int* M, const short* map,
                            double* Rw, double* sF, double* Dq, double* Dket, int lane, double acc[9])
{
    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const int nB = d.sh_ncomp[B], nC = d.sh_ncomp[C], nD = d.sh_ncomp[D];
    const int nCD = nC * nD, nq = d.sh_ncomp[A] * nB * nCD;
    const int oA = d.sh_out_first[A], oB = d.sh_out_first[B], oC = d.sh_out_first[C], oD = d.sh_out_first[D];
    const size_t n = (size_t)d.n_out;
    const double* P = g.P;
    gsync<GS>();  // previous quartet is done with Dq
    const double* P2 = g.P2;
    const double* G = g.G;
    // step 1: D_abcd as gradientTwoElectronMP2 (RHF + P2 x P + G)
    for (int e = lane; e < nq; e += GS) {
        const int ab = e / nCD, cd = e % nCD;
        const size_t m = oA + ab / nB, nn = oB + ab % nB, l = oC + cd / nD, sg = oD + cd % nD;
        Dq[e] = P[m * n + nn] * P[l * n + sg] - 0.25 * (P[m * n + l] * P[nn * n + sg] + P[m * n + sg] * P[nn * n + l])
              + (P2[m * n + nn] * P[l * n + sg] + P[m * n + nn] * P2[l * n + sg])
              - 0.25 * (P2[m * n + l] * P[nn * n + sg] + P[m * n + l] * P2[nn * n + sg] + P2[m * n + sg] * P[nn * n + l]
                        + P[m * n + sg] * P2[nn * n + l])
              + G[((m * n + nn) * n + l) * n + sg];
    }
    gsync<GS>();
    double o1[6] = { 0, 0, 0, 0, 0, 0 }, o2[6] = { 0, 0, 0, 0, 0, 0 };
    orientationGrad<GS, true>(d, g, i, j, Dq, false, g.prim1, M, map, Rw, sF, Dket, lane, o1);
    orientationGrad<GS, false>(d, g, j, i, Dq, true, g.prim2, M, map, Rw, sF, Dket, lane, o2);
#pragma unroll
    for (int k = 0; k < 6; ++k) acc[k] = o1[k];
#pragma unroll
    for (int k = 0; k < 3; ++k) acc[6 + k] = o2[k];
}

/// Add f * (gA, gB, gC, -(gA+gB+gC)) of one quartet to the atoms of its four shells.
__device__ inline void scatterQuartet(const DevTables& d, const GradArgs& g, double* sgrad, int i, int j,
                                      const double acc[9])
{
    const int A = d.pr_A[i], B = d.pr_B[i], C = d.pr_A[j], D = d.pr_B[j];
    const double f = 0.5 * (A != B ? 2.0 : 1.0) * (C != D ? 2.0 : 1.0) * (i != j ? 2.0 : 1.0);
    const int at[4] = { g.sh_atom[A], g.sh_atom[B], g.sh_atom[C], g.sh_atom[D] };
    double* dst = g.shared_grad ? sgrad : g.grad;
#pragma unroll
    for (int k = 0; k < 3; ++k) {
        const double a = f * acc[k], b = f * acc[3 + k], c = f * acc[6 + k];
        atomicAdd(&dst[3 * at[0] + k], a);
        atomicAdd(&dst[3 * at[1] + k], b);
        atomicAdd(&dst[3 * at[2] + k], c);
        atomicAdd(&dst[3 * at[3] + k], -(a + b + c));
    }
}

/// Dynamic shared memory: [shared gradient][per-warp scratch (warp kernel)][map][M].
__device__ inline void blockSetup(const GradArgs& g, double* smem, int nwarp, double*& sgrad, double*& wbase,
                                  short*& map, int*& M)
{
    const int s = g.L + 1, s3 = s * s * s;
    sgrad = smem;
    wbase = smem + (g.shared_grad ? 3 * g.natoms : 0);
    map = reinterpret_cast<short*>(wbase + (size_t)nwarp * g.per_warp);
    M = reinterpret_cast<int*>(map + ((s3 + 3) / 4) * 4);
    buildHermiteMaps(g.L, map, M);
    if (g.shared_grad)
        for (int k = threadIdx.x; k < 3 * g.natoms; k += blockDim.x) sgrad[k] = 0.0;
    __syncthreads();
}

__device__ inline void flushSharedGrad(const GradArgs& g, const double* sgrad)
{
    __syncthreads();
    if (!g.shared_grad) return;
    for (int k = threadIdx.x; k < 3 * g.natoms; k += blockDim.x)
        if (sgrad[k] != 0.0) atomicAdd(&g.grad[k], sgrad[k]);
}

/// One warp per shell quartet: R^n, D and Dket in shared memory, lanes over their entries.
__global__ void __launch_bounds__(128) grad2eWarpKernel(DevTables d, GradArgs g)
{
    extern __shared__ double gsm[];
    const int nwarp = blockDim.x / 32, warp = threadIdx.x / 32, lane = threadIdx.x % 32;
    double *sgrad, *wbase;
    short* map;
    int* M;
    blockSetup(g, gsm, nwarp, sgrad, wbase, map, M);
    double* Rw = wbase + (size_t)warp * g.per_warp;  // 2 NS
    double* sF = Rw + 2 * g.NS;                      // 24
    double* Dq = sF + 24;                            // ndq
    double* Dket = Dq + g.ndq;                       // ndket
    const long long nw = (long long)gridDim.x * nwarp;
    for (long long q = (long long)blockIdx.x * nwarp + warp; q < g.nq; q += nw) {
        int i, j;
        if (!gradQuartet(d, g, q, i, j)) continue;  // uniform per warp
        if (lane == 0) atomicAdd(g.counter, 1ull);
        double acc[9];
        quartetGrad<32>(d, g, i, j, M, map, Rw, sF, Dq, Dket, lane, acc);
#pragma unroll
        for (int m = 0; m < 9; ++m)
            for (int off = 16; off > 0; off >>= 1) acc[m] += __shfl_down_sync(0xffffffffu, acc[m], off);
        if (lane == 0) scatterQuartet(d, g, sgrad, i, j, acc);
    }
    flushSharedGrad(g, sgrad);
}

/// Pmax(A,B) = max |P| over the cartesian shell block (A, B).
__global__ void gradPmaxKernel(const double* __restrict__ P, const int* __restrict__ first, const int* __restrict__ width,
                               int ns, int n, double* __restrict__ Pmax)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= ns * ns) return;
    const int A = k / ns, B = k % ns;
    double m = 0.0;
    for (int a = 0; a < width[A]; ++a)
        for (int b = 0; b < width[B]; ++b) m = fmax(m, fabs(P[(size_t)(first[A] + a) * n + first[B] + b]));
    Pmax[k] = m;
}

class QmGpuMp2Grad2eCuda final : public qmint::QMGpuMP2Grad2e {
public:
    ~QmGpuMp2Grad2eCuda() override
    {
        for (void* p : m_owned) cudaFree(p);
        for (auto& st : m_streams)
            if (st) cudaStreamDestroy(st);
        if (m_e0) cudaEventDestroy(m_e0);
        if (m_e1) cudaEventDestroy(m_e1);
    }

    bool init(const qmint::QMGpuEriTables* t, const qmint::QMGpuGradExtra* x, int device, std::string& why)
    {
        if (!t->pc_off || !t->pr_schwarz) { why = "tables without pair classes"; return false; }
        if (t->n_T != 0) { why = "gradient tables must be cartesian (no 6d->5d transform)"; return false; }
        if (!x || x->nprimpair != t->nprimpair || x->natoms <= 0) { why = "gradient extra data inconsistent"; return false; }
        if (!uploadCommonTables(t, m_d, m_owned, why)) return false;
        m_n = t->n_out;
        m_ns = t->nshell;
        m_natoms = x->natoms;
        m_pc_off.assign(t->pc_off, t->pc_off + 10);
        // largest number of primitive pairs of a shell pair, per pair class (primitive-mode Dket size)
        for (int X = 0; X < 9; ++X) {
            m_maxprim[X] = 0;
            for (int k = t->pc_off[X]; k < t->pc_off[X + 1]; ++k) m_maxprim[X] = std::max(m_maxprim[X], t->pr_nprim[t->pc_pairs[k]]);
        }
        const size_t n2 = (size_t)m_n * m_n;
        bool ok = upload(const_cast<double**>(&m_d.pr_schwarz), t->pr_schwarz, (size_t)t->npair, m_owned)
            && upload(&m_pc_pairs, t->pc_pairs, (size_t)t->npair, m_owned)
            && upload(&m_width, t->sh_out_width, (size_t)t->nshell, m_owned)
            && upload(&m_sh_atom, x->sh_atom, (size_t)t->nshell, m_owned)
            && upload(&m_pp_ea, x->pp_ea, (size_t)t->nprimpair, m_owned)
            && upload(&m_pp_eb, x->pp_eb, (size_t)t->nprimpair, m_owned)
            && alloc(&m_P, n2) && alloc(&m_P2, n2) && alloc(&m_Pmax, (size_t)m_ns * m_ns) && alloc(&m_grad, 3 * (size_t)m_natoms)
            && alloc(&m_counter, 1);
        if (!ok) { why = "device allocation / upload for the 2e gradient failed"; return false; }
        for (auto& st : m_streams)
            if (cudaStreamCreateWithFlags(&st, cudaStreamNonBlocking) != cudaSuccess) { why = "stream creation failed"; return false; }
        cudaEventCreate(&m_e0);
        cudaEventCreate(&m_e1);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        m_sms = prop.multiProcessorCount > 0 ? prop.multiProcessorCount : 16;
        return cudaGetLastError() == cudaSuccess;
    }

    bool compute(const double* P, const double* P2, const double* G, double* grad) override
    {
        const double threshold = 0.0;  // every quartet, as qmint::gradientTwoElectronMP2
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n2 = (size_t)m_n * m_n;
        if (!m_G && !alloc(&m_G, n2 * n2)) return false;
        if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess
            || cudaMemcpy(m_P2, P2, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess
            || cudaMemcpy(m_G, G, n2 * n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
            return false;
        m_upload_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        cudaMemset(m_grad, 0, 3 * (size_t)m_natoms * sizeof(double));
        cudaMemset(m_counter, 0, sizeof(unsigned long long));
        cudaEventRecord(m_e0);
        const int nsn = m_ns * m_ns;
        gradPmaxKernel<<<(nsn + 127) / 128, 128>>>(m_P, m_d.sh_out_first, m_width, m_ns, m_n, m_Pmax);
        cudaEventRecord(m_e1);
        for (auto& st : m_streams) cudaStreamWaitEvent(st, m_e1, 0);

        GradArgs g{};
        g.sh_atom = m_sh_atom;
        g.pp_ea = m_pp_ea;
        g.pp_eb = m_pp_eb;
        g.P = m_P;
        g.P2 = m_P2;
        g.G = m_G;
        g.Pmax = m_Pmax;
        g.thr = threshold;
        g.natoms = m_natoms;
        g.grad = m_grad;
        g.counter = m_counter;
        const bool prof = std::getenv("CURCUMA_QM_GPU_GRAD_PROFILE") != nullptr;  // per-class timing, serialised
        m_prof.clear();
        m_tprof = std::chrono::steady_clock::now();
        int launch = 0;
        for (int X = 0; X < 9; ++X) {
            const long long nX = m_pc_off[X + 1] - m_pc_off[X];
            if (nX == 0) continue;
            for (int Y = 0; Y < 9; ++Y) {
                const long long nY = m_pc_off[Y + 1] - m_pc_off[Y];
                if (nY == 0) continue;
                g.LA = X / 3; g.LB = X % 3; g.LC = Y / 3; g.LD = Y % 3;
                g.L = g.LA + g.LB + g.LC + g.LD + 1;
                g.rb = m_pc_pairs + m_pc_off[X];
                g.rk = m_pc_pairs + m_pc_off[Y];
                g.nk = nY;
                g.nq = nX * nY;
                const int nAB = ncompOf(g.LA) * ncompOf(g.LB), nCD = ncompOf(g.LC) * ncompOf(g.LD);
                g.NS = nsimplex(g.L);
                g.ndq = nAB * nCD;
                // Primitive mode for an orientation whose bra has fewer Hermite entries than lanes (and
                // whose R^0 fits the thread-local array); its Dket then holds every ket primitive pair.
                const int e1 = nAB * nsimplex(g.LA + g.LB + 1), e2 = nCD * nsimplex(g.LC + g.LD + 1);
                g.prim1 = e1 < 32 && g.NS <= kPrimNS;
                g.prim2 = e2 < 32 && g.NS <= kPrimNS;
                const int dk1 = nAB * nsimplex(g.LC + g.LD) * (g.prim1 ? m_maxprim[Y] : 1);
                const int dk2 = nCD * nsimplex(g.LA + g.LB) * (g.prim2 ? m_maxprim[X] : 1);
                g.ndket = std::max(dk1, dk2);
                g.per_warp = 2 * g.NS + 24 + g.ndq + g.ndket;
                const int s = g.L + 1, s3 = s * s * s;
                const size_t maps = (size_t)((s3 + 3) / 4) * 4 * sizeof(short) + (size_t)g.NS * sizeof(int);
                g.shared_grad = m_natoms <= kSharedGradMaxAtoms;
                const size_t sgrad = g.shared_grad ? 3 * (size_t)m_natoms * sizeof(double) : 0;
                cudaStream_t st = m_streams[launch++ % kStreams];
                const long long maxBlocks = 8LL * m_sms;
                const size_t warp_bytes = (size_t)g.per_warp * sizeof(double);
                const int wpb = (int)std::min<size_t>(kMaxWarpsPerBlock, (kSharedBudget - sgrad - maps) / warp_bytes);
                if (wpb < 1) return false;  // does not happen for s, p, d shells and contraction depth <= 30
                const long long nb = std::min((g.nq + wpb - 1) / wpb, maxBlocks);
                grad2eWarpKernel<<<(unsigned)nb, 32 * wpb, sgrad + maps + (size_t)wpb * warp_bytes, st>>>(m_d, g);
                if (cudaGetLastError() != cudaSuccess) return false;
                if (prof) {
                    if (cudaDeviceSynchronize() != cudaSuccess) return false;
                    const auto tn = std::chrono::steady_clock::now();
                    m_prof += fmt_class(g, std::chrono::duration<double, std::milli>(tn - m_tprof).count());
                    m_tprof = tn;
                }
            }
        }
        for (auto& st : m_streams)
            if (cudaStreamSynchronize(st) != cudaSuccess) return false;
        cudaEventRecord(m_e1);
        if (cudaEventSynchronize(m_e1) != cudaSuccess || cudaGetLastError() != cudaSuccess) return false;
        float kms = 0.0f;
        cudaEventElapsedTime(&kms, m_e0, m_e1);
        unsigned long long cnt = 0;
        if (cudaMemcpy(grad, m_grad, 3 * (size_t)m_natoms * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
            || cudaMemcpy(&cnt, m_counter, sizeof cnt, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        m_last_quartets = (long long)cnt;
        m_kernel_ms = kms;
        m_total_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    std::string name() const override { return m_name; }
    long long lastQuartets() const override { return m_last_quartets; }
    std::string stats() const override
    {
        char buf[256];
        std::snprintf(buf, sizeof buf, "%lld quartets, kernels %.1f ms, 2PDM upload %.1f ms, total %.1f ms, device memory %.1f MB",
                      m_last_quartets, m_kernel_ms, m_upload_ms, m_total_ms, (m_bytes + (m_G ? (double)m_n * m_n * m_n * m_n * 8.0 : 0.0)) / 1.0e6);
        return std::string(buf) + m_prof;
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
    static std::string fmt_class(const GradArgs& g, double ms)
    {
        char buf[128];
        std::snprintf(buf, sizeof buf, "; (%d%d|%d%d) %s%s nq=%lld %.2f ms", g.LA, g.LB, g.LC, g.LD, g.prim1 ? "P" : "E",
                      g.prim2 ? "P" : "E", g.nq, ms);
        return buf;
    }

    static constexpr int kStreams = 8;
    DevTables m_d{};
    std::vector<void*> m_owned;
    std::vector<long long> m_pc_off;
    int *m_pc_pairs = nullptr, *m_width = nullptr, *m_sh_atom = nullptr;
    double *m_pp_ea = nullptr, *m_pp_eb = nullptr, *m_P = nullptr, *m_Pmax = nullptr, *m_grad = nullptr;
    double *m_P2 = nullptr, *m_G = nullptr;
    double m_upload_ms = 0.0;
    unsigned long long* m_counter = nullptr;
    cudaStream_t m_streams[kStreams] = {};
    cudaEvent_t m_e0 = nullptr, m_e1 = nullptr;
    int m_n = 0, m_ns = 0, m_natoms = 0, m_sms = 16;
    int m_maxprim[9] = {};
    long long m_last_quartets = 0;
    double m_kernel_ms = 0.0, m_total_ms = 0.0;
    size_t m_bytes = 0;
    std::string m_name = "cuda", m_prof;
    std::chrono::steady_clock::time_point m_tprof = std::chrono::steady_clock::now();
};

} // namespace

// MP2 2e-gradient builder (QMGpuMP2Grad2e*), or nullptr with a reason in `why`. `t` from
// qmint::buildEriGradGpuTables (cartesian, raised Hermite tables). Claude Generated (Oct 2026, MP2 gradient).
extern "C" void* curcuma_cuda_create_qm_mp2_grad2e(const qmint::QMGpuEriTables* t, const qmint::QMGpuGradExtra* x, int device,
                                                   char* why, int whylen)
{
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return nullptr; }
    cudaFree(nullptr);  // context creation
    size_t free0 = 0, total = 0;
    cudaMemGetInfo(&free0, &total);
    const double need = (double)t->n_out * t->n_out * t->n_out * t->n_out * 8.0;
    if (need > 0.9 * (double)free0) {
        setWhy(why, whylen, "cartesian 2PDM (" + std::to_string((long long)(need / 1.0e6)) + " MB) does not fit into free device memory");
        return nullptr;
    }
    auto b = std::make_unique<QmGpuMp2Grad2eCuda>();
    std::string err;
    if (!b->init(t, x, device, err)) {
        cudaGetLastError();
        setWhy(why, whylen, err);
        return nullptr;
    }
    size_t free1 = 0;
    cudaMemGetInfo(&free1, &total);
    b->setBytes(free0 > free1 ? free0 - free1 : 0);
    return static_cast<qmint::QMGpuMP2Grad2e*>(b.release());
}
