/*
 * <Analytic nuclear gradient of canonical closed-shell MP2 (native RHF)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Contents:
 *   - qmint::mp2GradientHost         amplitudes, unrelaxed density, Lagrangian intermediates and the
 *                                    back-transformed amplitudes (AO 2PDM) on the CPU
 *   - QMEngine::mp2G4                C^T (4 J - 2 K)(C X C^T) C with the exact J/K builder
 *   - QMEngine::computeMP2Gradient   Z-vector, relaxed and energy-weighted densities, derivative integrals
 *
 * Claude Generated (Oct 2026). AI-generated, machine-tested (see the WP11 section of
 * docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md for what was compared).
 *
 * This program is free software under GPL-3.0
 */

// =================================================================================================
// Derivation (closed shell, real canonical RHF orbitals; written out here and checked numerically
// against PySCF and against the finite difference of curcuma's own energy, not copied from a paper).
// References for the structure of the result: J. A. Pople, R. Krishnan, H. B. Schlegel, J. S. Binkley,
// Int. J. Quantum Chem. Symp. 13, 225 (1979); N. C. Handy, H. F. Schaefer, J. Chem. Phys. 81, 5031
// (1984) (Z-vector); C. M. Aikens, S. P. Webb, R. L. Bell, G. D. Fletcher, M. W. Schmidt, M. S. Gordon,
// Theor. Chem. Acc. 110, 233 (2003) (frozen core, closed and open shell).
//
// Indices: I,J frozen core, i,j,k active occupied, a,b,c virtual, p,q,r,s any MO, k,l in sums over
// "occ" run over ALL doubly occupied (frozen + active). Chemists' notation (pq|rs).
//
// Energy:  E = E_HF + E2,  E2 = sum_{ijab} (ia|jb) T~_ijab,
//   t_ijab = (ia|jb)/D_ijab,  D_ijab = e_i + e_j - e_a - e_b,
//   T~_ijab = (c_os + c_ss) t_ijab - c_ss t_ijba        (MP2: c_os = c_ss = 1 -> 2t - t^ba;
//                                                         SCS: 1.2, 1/3; SOS: 1.3, 0).
// E2 = sum K (c1 K - c_ss K^ba)/D is a quadratic form in K = (ia|jb), so dE2/dK_iajb = 2 T~_ijab and
// dE2/dD_ijab = -t_ijab T~_ijab. Written with a general Fock matrix (Hylleraas functional, stationary in
// t), the derivative at fixed t is
//   dE2 = 2 sum T~ dK + sum_ij P_ij df_ij + sum_ab P_ab df_ab,
//   P_ij = -2 sum_kab T~_ikab t_jkab,   P_ab = 2 sum_ijc T~_ijac t_ijbc      (unrelaxed density).
//
// Orbital response: C^x = C U, U_pq + U_qp = -S^(x)_pq (S^(x) = C^T S^x C, skeleton derivative).
//   K^x_iajb = (ia|jb)^(x) + sum_p [U_pi (pa|jb) + U_pa (ip|jb) + U_pj (ia|pb) + U_pb (ia|jp)],
//   f^x_pq   = f^(x)_pq + U_qp e_q + U_pq e_p + sum_{r, k occ} U_rk A_pq,rk,
//   A_pq,rs  = 4 (pq|rs) - (pr|qs) - (ps|qr),   sum_pq X_pq A_pq,rs = [C^T (4 J - 2 K)(C X C^T) C]_rs =: G4(X)_rs,
// with f^(x)_pq = h^(x)_pq + sum_k [2 (pq|kk)^(x) - (pk|qk)^(x)] (AO derivative integrals, fixed C).
// Using T~_ijab = T~_jiba, all U terms collect into sum_rs U_rs Y_rs with
//   Y_pi = L1_pi + 2 P_pi e_p [p active] + G4(P)_pi,   L1_pi = 4 sum_jab T~_ijab (pa|jb)   (i active)
//   Y_pI = G4(P)_pI                                                                            (I frozen)
//   Y_pa = L2_pa + 2 P_pa e_p [p virtual],              L2_pa = 4 sum_ijb T~_ijab (ip|jb)
// Pairs within an invariant block (active-active, virtual-virtual, frozen-frozen) have a symmetric Y and
// only need U_rs + U_sr = -S_rs. The two non-redundant pair types:
//  * frozen-active (I,j): f_Ij = 0 for canonical orbitals gives
//      U_Ij (e_I - e_j) = -[f^(x)_Ij - S_Ij e_j - 1/2 sum_kl S_kl A_Ij,kl + sum_bl U_bl A_Ij,bl],
//    so with w_Ij = (Y_Ij - Y_jI)/(e_j - e_I) = L1_Ij/(e_j - e_I) the pair contributes
//      w_Ij [f^(x)_Ij - S_Ij e_j - 1/2 sum_kl S_kl A_Ij,kl + sum_bl U_bl A_Ij,bl] - S_Ij Y_jI;
//  * virtual-occupied (a,k), k frozen or active: the CPHF equations (Delta + A) U = B,
//      (e_a - e_k) U_ak + sum_bl A_ak,bl U_bl = -[f^(x)_ak - S_ak e_k - 1/2 sum_ml S_ml A_ak,ml].
//    Z-vector (Handy-Schaefer): with the orbital gradient
//      L_ak = Y_ak - Y_ka + sum_Ij w_Ij A_Ij,ak = L1_ak [k active] - L2_ka + G4(P + W_fz)_ak
//    (W_fz: symmetric, w_Ij/2 in the frozen-active blocks) solve (e_a - e_k) z_ak + G4(Z)_ak = -L_ak
//    (Z symmetric, z_ak/2 in the virtual-occupied blocks); then sum U_ak L_ak = sum z_ak [f^(x)_ak - ...].
//
// Result. Relaxed MP2 density (MO, symmetric):  P2 = P (oo, vv) + w/2 (frozen-active) + z/2 (vir-occ);
// energy-weighted correction Ws (MO, symmetric) such that E2^x contains + sum_pq Ws_pq S^(x)_pq:
//   occ-occ (all k,l)        -1/2 G4(P2)_kl
//   active-active (i,j)      -1/4 (Y'_ij + Y'_ji),  Y'_ij = L1_ij + 2 P_ij e_i
//   virtual-virtual (a,b)    -1/4 (Y'_ab + Y'_ba),  Y'_ab = L2_ab + 2 P_ab e_a
//   virtual-occupied (a,k)   1/2 (-L2_ka - z_ak e_k)          (both (a,k) and (k,a))
//   frozen-active (I,j)      -1/2 w_Ij e_j                     (both (I,j) and (j,I))
// and the gradient
//   dE/dX = sum (P_HF + P2)_mn h^X_mn - sum (W_HF - Ws)_mn S^X_mn + 1/2 sum D_mnls (mn|ls)^X + dE_nn/dX,
//   D = P_HF P_HF - 1/4 (P_HF P_HF + P_HF P_HF)                 (RHF, as qmint::gradientTwoElectron)
//     + P2_mn P_ls + P_mn P2_ls - 1/4 (P2_ml P_ns + P_ml P2_ns + P2_ms P_nl + P_ms P2_nl)
//     + 4 Gamma^sym,  Gamma_mnls = sum_ijab T~_ijab C_mi C_na C_lj C_sb  (sym: 8-fold average),
// all densities and Ws taken to the AO basis with C (. ) C^T. The separable P2 x P_HF term is the skeleton
// part of sum P2_pq f^(x)_pq, the 4 Gamma^sym term is 2 sum T~ (ia|jb)^(x).
//
// Evaluation (o active occupied, v virtual, n basis functions; the energy code's batching over i):
//   pass 1 per i batch: (i nu|lam sig) (first quarter), M_i = (i nu|jb), K_i = (ia|jb);
//   amplitudes, energies, P_ij, P_ab, Y2_{nu a} = sum_ijb (i nu|jb) T~_ijab  (L2 = 4 C^T Y2);
//   pass 2 per i batch: R_i^{nu lam sig} = sum_jab T~_ijab C_nu a C_lam j C_sig b,
//     Y1_{kappa i} = sum (kappa nu|lam sig) R_i^{nu lam sig} (L1 = 4 C^T Y1), Gamma += C_mu i R_i.
// Memory: the n^4 AO tensor plus the n^4 Gamma plus O(o v o v + o n o v + batch n^3). The heavy GEMMs use
// qmint::gemmABt (the linked BLAS is single-threaded MKL on the development machine) threaded with OpenMP.
// =================================================================================================

#include "qm_citations.h"
#include "qm_engine.h"
#include "qm_mp2.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <vector>

namespace {
double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

/// out (cols x rows) = in (rows x cols)^T, both row-major.
void transposeInto(const double* in, size_t rows, size_t cols, double* out)
{
    for (size_t r = 0; r < rows; ++r)
        for (size_t c = 0; c < cols; ++c) out[c * rows + r] = in[r * cols + c];
}
}  // namespace

namespace qmint {

// Claude Generated (Oct 2026): see the derivation block above.
void mp2GradientHost(const double* eri, int n, const Matrix& C, const Vector& eps, int n_frozen, int n_occ,
                     double c_os, double c_ss, int threads, double batch_mb, double* Gamma, MP2GradIntermediates& out,
                     MP2Result& res)
{
    const int o = n_occ - n_frozen;
    const int v = n - n_occ;
    threads = std::max(1, threads);
    const size_t n2 = (size_t)n * n, n3 = n2 * n, n4 = n3 * n;
    const size_t ov = (size_t)o * v, vov = (size_t)v * ov, ovov = ov * ov;
    res.e_os = res.e_ss = 0.0;
    out.Poo = Matrix::Zero(o, o);
    out.Pvv = Matrix::Zero(v, v);
    out.Y1 = Matrix::Zero(n, o);
    out.Y2 = Matrix::Zero(n, v);
    {  // parallel first touch / zeroing of the n^4 output
        const size_t chunk = (size_t)1 << 20;
        const long long nch = (long long)((n4 + chunk - 1) / chunk);
#ifdef _OPENMP
#pragma omp parallel for num_threads(std::max(1, threads)) schedule(static)
#endif
        for (long long c = 0; c < nch; ++c) {
            const size_t x0 = (size_t)c * chunk;
            std::memset(Gamma + x0, 0, std::min(chunk, n4 - x0) * sizeof(double));
        }
    }
    if (o <= 0 || v <= 0) return;
    const auto t_amp0 = std::chrono::steady_clock::now();

    // Row-major MO blocks and their transposes.
    std::vector<double> Co((size_t)n * o), Cv((size_t)n * v), CoT((size_t)o * n), CvT((size_t)v * n);
    for (int m = 0; m < n; ++m) {
        for (int i = 0; i < o; ++i) Co[(size_t)m * o + i] = CoT[(size_t)i * n + m] = C(m, n_frozen + i);
        for (int a = 0; a < v; ++a) Cv[(size_t)m * v + a] = CvT[(size_t)a * n + m] = C(m, n_occ + a);
    }
    const double per_i_mb = (double)n3 * 8.0 / 1.0e6;
    const int batch = std::max(1, std::min(o, (int)std::floor(batch_mb / per_i_mb)));

    // ---- pass 1: (i nu|jb) and (ia|jb) for every active i ----
    std::vector<double> Kall(ovov);           // [i][a][j][b] = (ia|jb), later reused for T~
    std::vector<double> Mall((size_t)o * n * ov);  // [i][nu][j][b] = (i nu|jb)
    {
        std::vector<double> X1((size_t)batch * n3);
        const Eigen::Map<const Eigen::MatrixXd> Ecm(eri, (Eigen::Index)n3, n);  // column-major view, as mp2EnergyHost
        Eigen::MatrixXd CoCm(n, o);
        for (int m = 0; m < n; ++m)
            for (int i = 0; i < o; ++i) CoCm(m, i) = Co[(size_t)m * o + i];
        for (int i0 = 0; i0 < o; i0 += batch) {
            const int ib = std::min(batch, o - i0);
            // (i nu|lam sig) = sum_mu C_mu,i (mu nu|lam sig): tall-skinny E_cm (n^3 x n) * Co, OpenMP over row chunks.
            // Oct 31 2026: the chunk size is FIXED (8192 rows), not derived from the thread count: the GEMM of a row
            // block may round differently with another block height, so a thread-dependent split made the MP2 and
            // canonical double-hybrid gradients differ in the last bits between thread counts (now 1/4/8 bitwise;
            // the bits differ from a06631c0 at the 1e-15 level). mp2EnergyHost (qm_mp2.cpp) keeps its split.
            {
                const Eigen::Index rows = (Eigen::Index)n3;
                const Eigen::Index chunk = 8192;
                const Eigen::Index nchunk = (rows + chunk - 1) / chunk;
                Eigen::Map<Eigen::MatrixXd> X(X1.data(), rows, ib);
                const Eigen::MatrixXd Cb = CoCm.middleCols(i0, ib);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
                for (Eigen::Index c = 0; c < nchunk; ++c) {
                    const Eigen::Index r0 = c * chunk, len = std::min(chunk, rows - r0);
                    X.middleRows(r0, len).noalias() = Ecm.middleRows(r0, len) * Cb;
                }
            }
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
            {
                std::vector<double> T((size_t)o * n), Mt(ov * n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (int ii = 0; ii < ib; ++ii) {
                    const int i = i0 + ii;
                    const double* X1i = X1.data() + (size_t)ii * n3;
                    double* Mi = Mall.data() + (size_t)i * n * ov;
                    for (int nu = 0; nu < n; ++nu) {
                        // X1i[nu] is symmetric in (lam, sig): T[j][sig] = sum_lam C_lam,j X1i[nu][lam][sig]
                        gemmABt(CoT.data(), o, X1i + (size_t)nu * n2, n, n, T.data());
                        gemmABt(T.data(), o, CvT.data(), v, n, Mi + (size_t)nu * ov);  // [j][b]
                    }
                    transposeInto(Mi, n, ov, Mt.data());
                    gemmABt(CvT.data(), v, Mt.data(), (int)ov, n, Kall.data() + (size_t)i * vov);  // [a][jb]
                }
            }
        }
    }

    // ---- amplitudes, energies; T~ overwrites K ----
    std::vector<double> tall(ovov);
    std::vector<double> os_i(o, 0.0), ss_i(o, 0.0);
    const double c1 = c_os + c_ss;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (int i = 0; i < o; ++i) {
        const double ei = eps(n_frozen + i);
        const double* Ki = Kall.data() + (size_t)i * vov;
        double* ti = tall.data() + (size_t)i * vov;
        double os = 0.0, ss = 0.0;
        for (int a = 0; a < v; ++a) {
            const double ea = eps(n_occ + a);
            for (int j = 0; j < o; ++j) {
                const double eija = ei + eps(n_frozen + j) - ea;
                for (int b = 0; b < v; ++b) {
                    const size_t x = ((size_t)a * o + j) * v + b;
                    const double k = Ki[x];
                    const double kx = Ki[((size_t)b * o + j) * v + a];
                    const double inv_d = 1.0 / (eija - eps(n_occ + b));
                    ti[x] = k * inv_d;
                    os += k * k * inv_d;
                    ss += k * (k - kx) * inv_d;
                }
            }
        }
        os_i[i] = os;
        ss_i[i] = ss;
        double* Tti = Kall.data() + (size_t)i * vov;  // T~ in place of K (K is no longer needed for this i)
        for (int a = 0; a < v; ++a)
            for (int j = 0; j < o; ++j)
                for (int b = 0; b < v; ++b)
                    Tti[((size_t)a * o + j) * v + b] = c1 * ti[((size_t)a * o + j) * v + b] - c_ss * ti[((size_t)b * o + j) * v + a];
    }
    for (int i = 0; i < o; ++i) {
        res.e_os += os_i[i];
        res.e_ss += ss_i[i];
    }
    const double* Tt = Kall.data();

    // ---- unrelaxed density and Y2 (per-i parts stored, summed in fixed order) ----
    {
        std::vector<double> pv_i((size_t)o * v * v), y2_i((size_t)o * n * v);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
        for (int i = 0; i < o; ++i) {
            gemmABt(Tt + (size_t)i * vov, v, tall.data() + (size_t)i * vov, v, (int)ov, pv_i.data() + (size_t)i * v * v);
            gemmABt(Mall.data() + (size_t)i * n * ov, n, Tt + (size_t)i * vov, v, (int)ov, y2_i.data() + (size_t)i * n * v);
        }
        for (int i = 0; i < o; ++i)
            for (int a = 0; a < v; ++a)
                for (int b = 0; b < v; ++b) out.Pvv(a, b) += 2.0 * pv_i[((size_t)i * v + a) * v + b];
        for (int i = 0; i < o; ++i)
            for (int m = 0; m < n; ++m)
                for (int a = 0; a < v; ++a) out.Y2(m, a) += y2_i[((size_t)i * n + m) * v + a];
        // P_ij = -2 sum_{akb} T~[i][akb] t[j][akb]: one gemmABt row block per thread item
        std::vector<double> po((size_t)o * o);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
        for (int i = 0; i < o; ++i)
            gemmABt(Tt + (size_t)i * vov, 1, tall.data(), o, (int)vov, po.data() + (size_t)i * o);
        for (int i = 0; i < o; ++i)
            for (int j = 0; j < o; ++j) out.Poo(i, j) = -2.0 * po[(size_t)i * o + j];
    }
    std::vector<double>().swap(Mall);
    std::vector<double>().swap(tall);
    res.t_grad_amplitudes_ms = msSince(t_amp0);

    // ---- pass 2: R_i, Y1 and Gamma ----
    const auto t_bt0 = std::chrono::steady_clock::now();
    std::vector<double> Rb((size_t)batch * n3);
    std::vector<double> y1((size_t)n * o, 0.0);
    for (int i0 = 0; i0 < o; i0 += batch) {
        const int ib = std::min(batch, o - i0);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
        {
            std::vector<double> Ttr(ov * v), Z((size_t)n * ov), V((size_t)o * n), Vt((size_t)n * o);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int ii = 0; ii < ib; ++ii) {
                const int i = i0 + ii;
                transposeInto(Tt + (size_t)i * vov, v, ov, Ttr.data());            // [jb][a]
                gemmABt(Cv.data(), n, Ttr.data(), (int)ov, v, Z.data());             // Z[nu][j][b] = sum_a C_nu,a T~_ijab
                double* Ri = Rb.data() + (size_t)ii * n3;
                for (int nu = 0; nu < n; ++nu) {
                    gemmABt(Z.data() + (size_t)nu * ov, o, Cv.data(), n, v, V.data());  // V[j][sig] = sum_b Z[nu][j][b] C_sig,b
                    transposeInto(V.data(), o, n, Vt.data());
                    gemmABt(Co.data(), n, Vt.data(), n, o, Ri + (size_t)nu * n2);    // R[nu][lam][sig]
                }
            }
        }
        // Y1[kappa][i] = sum_{x} (kappa | x) R_i[x], x = (nu lam sig): rows of the tensor in blocks of 4
        {
            const int nblk = (n + 3) / 4;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
            {
                std::vector<double> tmp((size_t)4 * ib);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (int b = 0; b < nblk; ++b) {
                    const int k0 = 4 * b, rows = std::min(4, n - k0);
                    gemmABt(eri + (size_t)k0 * n3, rows, Rb.data(), ib, (int)n3, tmp.data());
                    for (int r = 0; r < rows; ++r)
                        for (int ii = 0; ii < ib; ++ii) y1[(size_t)(k0 + r) * o + i0 + ii] = tmp[(size_t)r * ib + ii];
                }
            }
        }
        // Gamma[mu][x] += sum_i C_mu,i R_i[x], chunks of x (disjoint writes)
        {
            const size_t chunk = 512;
            const size_t nchunk = (n3 + chunk - 1) / chunk;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
            {
                std::vector<double> acc(chunk);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
                for (size_t c = 0; c < nchunk; ++c) {
                    const size_t x0 = c * chunk, len = std::min(chunk, n3 - x0);
                    for (int mu = 0; mu < n; ++mu) {
                        std::fill(acc.begin(), acc.begin() + len, 0.0);
                        for (int ii = 0; ii < ib; ++ii) {
                            const double cmi = Co[(size_t)mu * o + i0 + ii];
                            const double* r = Rb.data() + (size_t)ii * n3 + x0;
#pragma omp simd
                            for (size_t x = 0; x < len; ++x) acc[x] += cmi * r[x];
                        }
                        double* g = Gamma + (size_t)mu * n3 + x0;
                        for (size_t x = 0; x < len; ++x) g[x] += acc[x];
                    }
                }
            }
        }
    }
    for (int m = 0; m < n; ++m)
        for (int i = 0; i < o; ++i) out.Y1(m, i) = y1[(size_t)m * o + i];
    res.t_grad_backtrans_ms = msSince(t_bt0);
}

// Claude Generated (Oct 2026): Gamma_cart = (Q x Q x Q x Q) Gamma_act in ONE pass over the output (memory traffic
// n_act^4 read + n_cart^4 written; the first version did four one-index passes over full-size intermediates and
// took 2.7 s for naphthalene). Q has one to three nonzeros per cartesian row (block-diagonal per shell):
//   out[c1 c2 c3 c4] = sum_{s4 in row c4} Q(c4,s4) tmp[s4],
//   tmp[s4]          = sum_{s1,s2,s3 in rows c1,c2,c3} Q(c1,s1) Q(c2,s2) Q(c3,s3) in[s1 s2 s3 s4]   (contiguous in s4).
void mp2ActiveToCartesian(const double* in, const Matrix& Q, std::unique_ptr<double[]>& out, int threads)
{
    const int nc = (int)Q.rows(), na = (int)Q.cols();
    std::vector<std::vector<std::pair<int, double>>> row(nc);
    for (int c = 0; c < nc; ++c)
        for (int s = 0; s < na; ++s)
            if (Q(c, s) != 0.0) row[c].push_back({ s, Q(c, s) });
    const size_t NA = (size_t)na, NC = (size_t)nc;
    out.reset(new double[NC * NC * NC * NC]);
    double* o = out.get();
#ifdef _OPENMP
#pragma omp parallel num_threads(std::max(1, threads))
#endif
    {
        std::vector<double> tmp(na);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (long long c12 = 0; c12 < (long long)nc * nc; ++c12) {
            const int c1 = (int)(c12 / nc), c2 = (int)(c12 % nc);
            for (int c3 = 0; c3 < nc; ++c3) {
                std::fill(tmp.begin(), tmp.end(), 0.0);
                for (const auto& r1 : row[c1])
                    for (const auto& r2 : row[c2])
                        for (const auto& r3 : row[c3]) {
                            const double q = r1.second * r2.second * r3.second;
                            const double* x = in + ((r1.first * NA + r2.first) * NA + r3.first) * NA;
#pragma omp simd
                            for (int s4 = 0; s4 < na; ++s4) tmp[s4] += q * x[s4];
                        }
                double* y = o + (((size_t)c1 * NC + c2) * NC + c3) * NC;
                for (int c4 = 0; c4 < nc; ++c4) {
                    double v = 0.0;
                    for (const auto& r4 : row[c4]) v += r4.second * tmp[r4.first];
                    y[c4] = v;
                }
            }
        }
    }
}

}  // namespace qmint

// =================================================================================================
// QMEngine part: Lagrangian, Z-vector, densities, derivative integrals
// =================================================================================================

// Claude Generated (Oct 2026): G4(X) = C^T (4 J(C X C^T) - 2 K(C X C^T)) C, X symmetric (MO basis).
Matrix QMEngine::mp2G4(const Matrix& X_mo, const qmint::ERITensor* eri) const
{
    const Matrix Xao = m_mo * X_mo * m_mo.transpose();
    Matrix J, K;
    if (eri) {
        J = qmint::buildCoulomb(*eri, Xao, m_threads);
        K = qmint::buildExchange(*eri, Xao, m_threads);
    } else {
        buildJKExact(Xao, J, K);
    }
    return m_mo.transpose() * (4.0 * J - 2.0 * K) * m_mo;
}

// Claude Generated (Oct 2026). See the derivation block at the top of this file.
bool QMEngine::computeMP2Gradient(const qmint::MP2Options& opt, double c_os, double c_ss, qmint::MP2Result& res)
{
    res = qmint::MP2Result();
    res.gradient = true;
    qmcite::mp2("mp2", true);   // Z-vector method (Handy & Schaefer 1984)
    m_gradient_failed = false;
    m_gradient_error.clear();
    auto fail = [&](const std::string& msg) {
        m_gradient = Matrix::Zero(m_atomcount, 3);
        m_gradient_failed = true;
        m_gradient_error = msg;
        CurcumaLogger::error(msg);
        return false;
    };
    if (m_functional != QMFunctional::HF) return fail("MP2 gradient: needs the HF reference");
    if (!m_scf_converged) return fail("MP2 gradient: the RHF SCF is not converged");
    if (riActive())
        return fail("MP2 gradient: needs the exact-integral SCF (-qm.ri off); the RHF part of the MP2 gradient "
                    "with an RI-J/RI-JK reference is not implemented");
    const int n = m_nbf;
    const int n_occ = m_num_electrons / 2;
    if (m_mo.rows() != n || m_mo.cols() != n || m_energies.size() != n)
        return fail("MP2 gradient: no MO coefficients / orbital energies from the SCF");
    int nf = 0;
    if (opt.frozen_core)
        for (int Z : m_atoms) nf += qmint::frozenCoreOrbitals(Z);
    if (nf >= n_occ && n_occ > 0)
        return fail(fmt::format("MP2 gradient: frozen core ({} orbitals) leaves no active occupied orbital", nf));
    const int o = n_occ - nf, v = n - n_occ;
    res.n_frozen = nf;
    res.n_occ_active = o;
    res.n_virt = v;
    res.n_basis = n;
    if (o == 0 || v == 0) {  // no correlation energy: the RHF gradient is the MP2 gradient
        computeGradient();
        return !m_gradient_failed;
    }
    const double gamma_mb = std::pow((double)n, 4) * 8.0 / 1.0e6;
    if (gamma_mb > opt.max_memory_mb)
        return fail(fmt::format("MP2 gradient: the AO two-particle density ({:.1f} MB) exceeds -qm.mp2_max_memory_mb {:.1f}",
                                gamma_mb, opt.max_memory_mb));

    // ---- amplitudes and intermediates (CPU) ----
    qmint::ERITensor local;
    std::string src;
    const double* eri = mp2HostTensor(opt, res, local, src);
    if (!eri) return fail("MP2 gradient: no AO tensor (see the error above)");
    std::unique_ptr<double[]> GammaBuf;  // uninitialised; mp2GradientHost zeroes it in parallel
    try {
        GammaBuf.reset(new double[(size_t)n * n * n * n]);
    } catch (const std::bad_alloc&) {
        return fail(fmt::format("MP2 gradient: cannot allocate the AO two-particle density ({:.1f} MB)", gamma_mb));
    }
    double* const Gamma = GammaBuf.get();
    qmint::MP2GradIntermediates im;
    qmint::mp2GradientHost(eri, n, m_mo, m_energies, nf, n_occ, c_os, c_ss, m_threads, opt.batch_mb, Gamma, im, res);
    res.info = src + fmt::format("; CPU, {} threads", m_threads);
    if (!std::isfinite(res.e_os) || !std::isfinite(res.e_ss))
        return fail("MP2 gradient: non-finite correlation energy (degenerate HOMO/LUMO?)");
    // A tensor rebuilt for MP2 (integral-direct SCF) also serves the J/K builds of the Z-vector: the integral-direct
    // builder has a screening noise floor (absolute, per quartet) that stalled the iterations at ~1e-9 (naphthalene).
    // With a device J/K builder of this geometry (stored device tensor) the device path is used instead: exact and faster.
    const bool dev_jk = m_gpu_jk && m_gpu_jk->n() == n && !m_gpu_failed;
    const qmint::ERITensor* zeri = (local.empty() || dev_jk) ? nullptr : &local;

    // ---- MO Lagrangian ----
    const auto t_z0 = std::chrono::steady_clock::now();
    const Matrix& C = m_mo;
    const Vector& e = m_energies;
    const Matrix L1 = 4.0 * (C.transpose() * im.Y1);  // n x o: L1_{p i}
    const Matrix L2 = 4.0 * (C.transpose() * im.Y2);  // n x v: L2_{p a}
    Matrix P2 = Matrix::Zero(n, n);
    P2.block(nf, nf, o, o) = im.Poo;
    P2.block(n_occ, n_occ, v, v) = im.Pvv;
    Matrix w = Matrix::Zero(std::max(nf, 1), o);  // frozen-active rotation multipliers
    for (int I = 0; I < nf; ++I)
        for (int jj = 0; jj < o; ++jj) {
            const int j = nf + jj;
            w(I, jj) = L1(I, jj) / (e(j) - e(I));
            P2(I, j) = P2(j, I) = 0.5 * w(I, jj);
        }
    const Matrix G = mp2G4(P2, zeri);
    Matrix Lag(v, n_occ);  // L_ak
    for (int a = 0; a < v; ++a)
        for (int k = 0; k < n_occ; ++k)
            Lag(a, k) = (k >= nf ? L1(n_occ + a, k - nf) : 0.0) - L2(k, a) + G(n_occ + a, k);

    // ---- Z-vector: (e_a - e_k) z_ak + G4(Z)_ak = -L_ak, preconditioned conjugate gradient ----
    auto toSym = [&](const Matrix& z) {
        Matrix Zs = Matrix::Zero(n, n);
        for (int a = 0; a < v; ++a)
            for (int k = 0; k < n_occ; ++k) Zs(n_occ + a, k) = Zs(k, n_occ + a) = 0.5 * z(a, k);
        return Zs;
    };
    auto applyM = [&](const Matrix& z) {
        const Matrix Gz = mp2G4(toSym(z), zeri);
        Matrix r(v, n_occ);
        for (int a = 0; a < v; ++a)
            for (int k = 0; k < n_occ; ++k) r(a, k) = (e(n_occ + a) - e(k)) * z(a, k) + Gz(n_occ + a, k);
        return r;
    };
    Matrix diag(v, n_occ);
    for (int a = 0; a < v; ++a)
        for (int k = 0; k < n_occ; ++k) diag(a, k) = e(n_occ + a) - e(k);
    const Matrix rhs = -Lag;
    Matrix z = rhs.cwiseQuotient(diag);
    Matrix r = rhs - applyM(z);
    Matrix pz = r.cwiseQuotient(diag);
    double rz = r.cwiseProduct(pz).sum();
    const double ztol = 1.0e-10;
    const int zmax = 100;
    int it = 0;
    double rmax = r.cwiseAbs().maxCoeff();
    while (rmax > ztol && it < zmax) {
        const Matrix Ap = applyM(pz);
        const double alpha = rz / pz.cwiseProduct(Ap).sum();
        z += alpha * pz;
        r -= alpha * Ap;
        ++it;
        rmax = r.cwiseAbs().maxCoeff();
        CurcumaLogger::result(fmt::format("MP2 Z-vector iteration {:3d}  max|R| = {:.2e}", it, rmax));
        if (rmax <= ztol) break;
        const Matrix s = r.cwiseQuotient(diag);
        const double rz_new = r.cwiseProduct(s).sum();
        pz = s + (rz_new / rz) * pz;
        rz = rz_new;
    }
    res.zvector_iterations = it;
    res.zvector_residual = rmax;
    if (rmax > 1.0e-6)
        return fail(fmt::format("MP2 gradient: Z-vector equations not converged ({} iterations, max residual {:.2e})", it, rmax));
    if (rmax > ztol && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("MP2 gradient: Z-vector stopped at max residual {:.2e} after {} iterations", rmax, it));
    res.t_grad_zvector_ms = msSince(t_z0);

    // ---- relaxed density, energy-weighted density, 2PDM ----
    const auto t_d0 = std::chrono::steady_clock::now();
    P2 += toSym(z);
    const Matrix Gt = mp2G4(P2, zeri);
    Matrix Ws = Matrix::Zero(n, n);
    for (int k = 0; k < n_occ; ++k)
        for (int l = 0; l < n_occ; ++l) Ws(k, l) = -0.5 * Gt(k, l);
    for (int ii = 0; ii < o; ++ii)
        for (int jj = 0; jj < o; ++jj) {
            const int i = nf + ii, j = nf + jj;
            const double yij = L1(i, jj) + 2.0 * P2(i, j) * e(i);
            const double yji = L1(j, ii) + 2.0 * P2(j, i) * e(j);
            Ws(i, j) += -0.25 * (yij + yji);
        }
    for (int aa = 0; aa < v; ++aa)
        for (int bb = 0; bb < v; ++bb) {
            const int a = n_occ + aa, b = n_occ + bb;
            const double yab = L2(a, bb) + 2.0 * P2(a, b) * e(a);
            const double yba = L2(b, aa) + 2.0 * P2(b, a) * e(b);
            Ws(a, b) += -0.25 * (yab + yba);
        }
    for (int aa = 0; aa < v; ++aa)
        for (int k = 0; k < n_occ; ++k) {
            const double x = 0.5 * (-L2(k, aa) - z(aa, k) * e(k));
            Ws(n_occ + aa, k) += x;
            Ws(k, n_occ + aa) += x;
        }
    for (int I = 0; I < nf; ++I)
        for (int jj = 0; jj < o; ++jj) {
            const double x = -0.5 * w(I, jj) * e(nf + jj);
            Ws(I, nf + jj) += x;
            Ws(nf + jj, I) += x;
        }
    Matrix Whf = Matrix::Zero(n, n);
    for (int i = 0; i < n_occ; ++i) Whf += 2.0 * e(i) * (C.col(i) * C.col(i).transpose());
    const Matrix P2ao = C * P2 * C.transpose();
    const Matrix Wtot = Whf - C * Ws * C.transpose();
    Matrix Pc = m_density, P2c = P2ao, Wc = Wtot;
    if (m_Q.size() != 0) {
        Pc = m_Q * m_density * m_Q.transpose();
        P2c = m_Q * P2ao * m_Q.transpose();
        Wc = m_Q * Wtot * m_Q.transpose();
    }
    // Non-separable part D_ns = 4 Gamma^sym = (sum of the 8 index permutations of Gamma) / 2, in place.
    {
        double* Gm = Gamma;
        const size_t N = (size_t)n;
        auto at = [&](size_t a, size_t b, size_t c, size_t d) -> double& { return Gm[((a * N + b) * N + c) * N + d]; };
#ifdef _OPENMP
#pragma omp parallel for num_threads(std::max(1, m_threads)) schedule(dynamic, 1)
#endif
        for (int m = 0; m < n; ++m)
            for (int nn = 0; nn <= m; ++nn) {
                const size_t mn = (size_t)m * (m + 1) / 2 + nn;
                for (int l = 0; l < n; ++l)
                    for (int s = 0; s <= l; ++s) {
                        const size_t ls = (size_t)l * (l + 1) / 2 + s;
                        if (ls > mn) break;
                        double& v0 = at(m, nn, l, s); double& v1 = at(nn, m, l, s);
                        double& v2 = at(m, nn, s, l); double& v3 = at(nn, m, s, l);
                        double& v4 = at(l, s, m, nn); double& v5 = at(s, l, m, nn);
                        double& v6 = at(l, s, nn, m); double& v7 = at(s, l, nn, m);
                        const double sum = 0.5 * (v0 + v1 + v2 + v3 + v4 + v5 + v6 + v7);
                        v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = sum;
                    }
            }
    }
    res.t_grad_density_ms = msSince(t_d0);

    // ---- derivative integrals ----
    Matrix posBohr = Matrix::Zero(m_atomcount, 3);
    for (int i = 0; i < m_atomcount; ++i)
        for (int k = 0; k < 3; ++k) posBohr(i, k) = CurcumaUnit::Length::angstrom_to_bohr(m_geometry(i, k));
    const auto t_1e0 = std::chrono::steady_clock::now();
    const Matrix Ptc = Pc + P2c;
    Matrix g1;
    std::string dev1_info;
    bool on_device1 = false;
    const std::string be1 = gpuGradBackend(m_gpu_grad1e, "gpu_grad1e", "qm_grad1e", m_gpu_grad1e_failed);
    // the 2e part on the device is controlled by -qm.gpu_grad (as the RHF Q-G7 part), counted as qm_mp2_grad
    const std::string be2 = gpuGradBackend(m_gpu_grad, "gpu_grad", "qm_mp2_grad", m_gpu_mp2_grad_failed);
    qmint::EriGpuTables tables;
    std::string tab_why;
    bool tables_ok = false;
    if (!be1.empty() || !be2.empty())
        tables_ok = qmint::buildEriGradGpuTables(m_gto_basis, m_eri_screening, tables, tab_why);
    if (!be1.empty())
        on_device1 = gpuGradientOneElectron(tables_ok ? &tables : nullptr, tab_why, be1, posBohr, Ptc, Wc, g1, dev1_info);
    if (!on_device1) g1 = qmint::gradientOneElectron(m_gto_basis, m_atoms, posBohr, Ptc, Wc, m_threads);
    res.t_grad_deriv1e_ms = msSince(t_1e0);
    const auto t_2e0 = std::chrono::steady_clock::now();
    Matrix g2;
    bool on_device2 = false;
    std::string dev2_info;
    if (!be2.empty()) {
        // cartesian non-separable 2PDM for the device (the CPU path transforms per shell quartet instead)
        std::string why = tab_why;
        std::unique_ptr<qmint::QMGpuMP2Grad2e> dev;
        if (tables_ok) {
            qmint::QMGpuGradExtra extra;
            extra.natoms = m_atomcount;
            extra.sh_atom = tables.sh_atom.data();
            extra.nprimpair = (int)tables.pp_ea.size();
            extra.pp_ea = tables.pp_ea.data();
            extra.pp_eb = tables.pp_eb.data();
            dev = gpu_plugin::createQmMp2Grad2e(be2, tables.view(), extra, m_gpu_device, why);
        }
        if (dev) {
            const auto tc = std::chrono::steady_clock::now();
            std::unique_ptr<double[]> Gc;
            const double* Gcart = Gamma;
            if (m_Q.size() != 0) {
                qmint::mp2ActiveToCartesian(Gamma, m_Q, Gc, m_threads);
                Gcart = Gc.get();
            }
            const double t_cart = msSince(tc);
            Matrix g = Matrix::Zero(m_atomcount, 3);
            const Matrix Pr = Pc, P2r = P2c;  // row-major contiguous
            if (dev->compute(Pr.data(), P2r.data(), Gcart, g.data())) {
                g2 = g;
                on_device2 = true;
                dev2_info = fmt::format("{} ({}; 2PDM to cartesian {:.1f} ms)", dev->name(), dev->stats(), t_cart);
            } else {
                why = "kernel failed";
            }
        }
        if (!on_device2) {
            m_gpu_mp2_grad_failed = true;
            curcuma::reportGpuFallback("qm_mp2_grad", "device MP2 2e gradient not possible: " + why);
        }
    }
    if (!on_device2 || (std::getenv("CURCUMA_QM_GPU_GRAD_CHECK") && std::getenv("CURCUMA_QM_GPU_GRAD_CHECK")[0] == '1')) {
        const auto tcpu = std::chrono::steady_clock::now();
        const Matrix gc = qmint::gradientTwoElectronMP2(m_gto_basis, Pc, P2c, Gamma, m_Q.size() != 0 ? &m_Q : nullptr,
                                                        m_atomcount, m_threads);
        if (on_device2) {
            // CURCUMA_QM_GPU_GRAD_CHECK=1: one comparison line (parsed by test_cases/qm_mp2/check_mp2_gradient.py)
            const Eigen::RowVector3d sd = g2.colwise().sum();
            CurcumaLogger::result(fmt::format("QM GPU MP2 grad check: natoms={} max|dev-cpu|={:.3e} max|cpu|={:.3e} "
                                              "sum_dev=({:.3e},{:.3e},{:.3e}) dev_ms={:.1f} cpu_ms={:.1f}",
                                              m_atomcount, (g2 - gc).cwiseAbs().maxCoeff(), gc.cwiseAbs().maxCoeff(), sd(0),
                                              sd(1), sd(2), msSince(t_2e0) - msSince(tcpu), msSince(tcpu)));
        } else {
            g2 = gc;
        }
    }
    res.t_grad_deriv2e_ms = msSince(t_2e0);
    const Matrix gn = calculateCoreRepulsionGradient();
    m_gradient = g1 + g2 + gn;
    m_gradient_parts = { g1, g2, gn };
    if (on_device1) res.info += "; 1e derivative integrals on the device (" + dev1_info + ")";
    if (on_device2) res.info += "; 2e derivative integrals on the device (" + dev2_info + ")";
    if (const char* tr = std::getenv("CURCUMA_QM_SCF_TRACE"); tr && tr[0] == '1')
        fmt::print(stderr, "QM MP2 gradient trace: amplitudes_ms={:.1f} backtrans_ms={:.1f} zvector_ms={:.1f} ({} it) "
                           "density_ms={:.1f} deriv1e_ms={:.1f} ({}) deriv2e_ms={:.1f} ({})\n",
                   res.t_grad_amplitudes_ms, res.t_grad_backtrans_ms, res.t_grad_zvector_ms, it, res.t_grad_density_ms,
                   res.t_grad_deriv1e_ms, on_device1 ? "device" : "cpu", res.t_grad_deriv2e_ms, on_device2 ? "device" : "cpu");
    return true;
}
