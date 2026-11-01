/*
 * < CUDA Kohn-Sham exchange-correlation of the native QM engine: AO values, density, V_xc, E_xc, XC gradient >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, GPU KS-DFT G1/G3 of docs/QM_GPU_ROADMAP.md section 4). AI-generated, machine-tested
// (RTX 5090 only).
//
// Device version of qmxc::buildVxc, qmxc::buildXCGradient and qmgrid::gridWeightGradient (formulas: qm_xc.h,
// qm_grid.h). The point kernels (functionals, Becke polynomials, one AO component) are the inline functions of
// qm_xc_kernels.h, the same source the CPU compiles. Per streamed block of B grid points:
//
//   AO        phi, d phi, d2 phi (B x nbf, column-major, point index fastest)   aoKernel, one thread per (point, shell)
//   density   X = phi P (cublasDgemm), rho = sum_m phi X, grad rho = 2 sum_m d phi X        densityXCKernel
//   XC        eps_x, eps_c, v_rho, v_sigma per point (qmxc::kern::evaluatePoint)            densityXCKernel
//   V_xc      Y = (w v_rho / 2) phi + 2 w v_sigma (grad rho . grad phi),  W += phi^T Y (DGEMM), V = W + W^T
//   gradient  Z = w v_rho phi + sum_j v_j d_j phi,  ZP = Z P (GEMM),
//             T_k = d_k phi o ZP + (GGA) sum_j v_j d_k d_j phi o (phi P);  basis part -2 sum_q T_k per function,
//             point part -2 sum_m T_k per point (grid response: points move with their atom)
//
// The grid (points, weights) stays on the device for the lifetime of the object (one geometry); AO values exist
// only per block (memory O(B nbf)). E_xc, N and the point counts are summed from per-point arrays over the WHOLE grid
// in fixed 256-point chunks (independent of B); V_xc is accumulated block by block (fixed order for a fixed B). No
// floating-point atomics: repeated runs are bitwise identical for the same B and device (checked on the RTX 5090,
// docs/QM_GPU_ROADMAP.md section 4). The Becke weight derivative runs one thread block per point at a time with
// O(natoms^2) scratch per block in global memory; per-block partial gradients are summed in block order.
//
// AO screening (Oct 2026, -qm.xc_ao_threshold > 0): the host sorts the grid into compact blocks and lists the
// shells that reach each block (qmgrid::BlockScreen, qm_grid.cpp; the one radius formula is
// AOEvaluator::shellRadii). Here consecutive blocks are merged into streamed groups (<= kScreenGroupPoints or
// block_points points) with the union S of their shells: AO values only for S (compact columns; zeros for a shell
// that the point's own host block does not keep, so each point sees the CPU's functions), P_SS gathered,
// X = phi_S P_SS, W_SS = phi_S^T Y_S (DGEMM) scattered into W(S, S) group by group (one stream, no atomics), the
// gradient columns scattered into gfun(S). Without screening the kernels run as before (shell list nullptr).
//
// Weight derivative (Oct 2026): no natoms limit any more. The old kernel kept 8 natoms doubles in shared memory
// (48 KB: natoms <= 700) and S, dS, mu, "P without one factor" (4 natoms^2) in global scratch per block. Now all
// per-point vectors are in global scratch, mu is recomputed from the distances (the same expression, so the same
// bits), S and dS are stored as one triangle of h and h' (2 natoms^2 in total with EX); a fixed number of virtual
// blocks (4 per SM, independent of the free memory) is launched in waves that fit the memory, so the summation
// order and the result do not depend on the free memory.
//
// Meta-GGA (Oct 2026, WP11 stage 3; TPSS, SCAN): device version of buildVxcMeta / buildXCGradientMeta (qm_xc.cpp).
// The point functional is qmxc::kern::evaluatePointMGGA (dual numbers, qm_xc_kernels.h), the same source as the CPU.
// Per block one strided-batched GEMM forms [X, DP_x, DP_y, DP_z] = [phi, d_x phi, d_y phi, d_z phi] P (four
// consecutive planes), tau_q = 1/2 sum_k sum_m d_k phi_qm (DP_k)_qm, then
//   V_xc      W += phi^T Y + sum_k (d_k phi)^T diag(w v_tau / 4) (d_k phi)   (V = W + W^T doubles the tau part to
//             sum_k (d_k phi)^T diag(w v_tau / 2) (d_k phi), the term of buildVxcMeta; four DGEMMs, fixed order)
//   gradient  T_k += (w v_tau / 2) sum_j d_k d_j phi o DP_j                    (AO deriv 2, as for a GGA)
// The LDA/GGA kernels are untouched (separate kernels; the gradient kernels are templates, <false> = the old code).
//
// Double hybrids (Oct 31 2026, WP12 stage 3): P86 (PZ81 / VWN5 local part) and B95 need no code here, they come in
// through evaluatePoint / evaluatePointMGGA of qm_xc_kernels.h (GGA and meta-GGA paths). For the double-hybrid
// gradient the object (QMGpuXCSetup::fxc) also forms F_xc[D] and the P2-directional XC gradient with the second
// derivatives of qm_xc_fxc_kernels.h: separate kernels (densityPairKernel, respRowKernel, respColKernel), so the
// LDA/GGA/meta-GGA kernels and, without the flag, their block sizes are unchanged.

#include "../qm_gpu_xc_backend.h"
#include "../qm_xc_fxc_kernels.h"
#include "../qm_xc_kernels.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr int kMaxTerms = 8;
constexpr int kMaxPrim = 32;      // primitives per shell
constexpr int kMaxComp = 15;      // cartesian components per shell (L <= 4; 10 until Nov 2026)
constexpr int kChunk = 256;       // points per energy-reduction chunk (fixed, independent of B)
constexpr int kScreenGroupPoints = 2048;  // AO screening: points per streamed group (measured, QM_GPU_ROADMAP section 4)

struct DevTerms {
    int n = 0;
    qmxc::XCTerm t[kMaxTerms];
};

struct DevShells {
    const double *xyz, *amin, *exps, *coef, *T;
    const int *nprim, *prim_off, *ncomp, *comp_off, *coef_off, *lmn, *out_first, *out_width, *T_off;
};

/// AO values (and derivatives) of shell blockIdx.y at points [begin, begin + count): ao[(k nbf + mu) B + q].
/// Same arithmetic as AOEvaluator::evaluate (qm_grid.cpp): radial sums R0, R1, R2 in primitive order, components via
/// qmgrid::kern::aoComponent, active functions as sum over the nonzero T entries in component order.
// AO screening: shl[blockIdx.y] is the shell and col[blockIdx.y] its first (compact) column; nbf is then the number
// of kept functions. The group consists of nblk host blocks starting at bstart[k] (relative to begin); a point of
// block k gets the values of shell blockIdx.y only when mask[k nsh + blockIdx.y] is set, else zeros: every point
// sees exactly the functions the CPU uses for it. shl == nullptr: every shell, columns out_first (unscreened layout).
__global__ void aoKernel(DevShells s, const double* __restrict__ pts, int begin, int count, int deriv, int nbf, int B,
                         double* __restrict__ ao, const int* __restrict__ shl, const int* __restrict__ col,
                         const int* __restrict__ bstart, const unsigned char* __restrict__ mask, int nblk)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const int sh = shl ? shl[blockIdx.y] : blockIdx.y;
    const int nder = deriv >= 2 ? 10 : (deriv >= 1 ? 4 : 1);
    const int p = begin + q;
    const double x = pts[3 * p] - s.xyz[3 * sh], y = pts[3 * p + 1] - s.xyz[3 * sh + 1], z = pts[3 * p + 2] - s.xyz[3 * sh + 2];
    const double r2 = x * x + y * y + z * z;
    const int nc = s.ncomp[sh], first = shl ? col[blockIdx.y] : s.out_first[sh], width = s.out_width[sh], toff = s.T_off[sh];
    const size_t plane = (size_t)nbf * B;
    bool keep = true;
    if (mask) {
        int k = 0;
        while (k + 1 < nblk && bstart[k + 1] <= q) ++k;
        keep = mask[(size_t)k * gridDim.y + blockIdx.y] != 0;
    }
    if (!keep || s.amin[sh] * r2 > 80.0) {
        for (int k = 0; k < nder; ++k)
            for (int w = 0; w < width; ++w) ao[k * plane + (size_t)(first + w) * B + q] = 0.0;
        return;
    }
    const int np = s.nprim[sh], po = s.prim_off[sh];
    double e[kMaxPrim], ae[kMaxPrim], aae[kMaxPrim];
    for (int k = 0; k < np; ++k) {
        const double a = s.exps[po + k];
        e[k] = exp(-a * r2);
        ae[k] = a * e[k];
        aae[k] = a * ae[k];
    }
    double cart[kMaxComp][10];
    for (int c = 0; c < nc; ++c) {
        const int ci = s.comp_off[sh] + c;
        const double* cc = s.coef + s.coef_off[ci];
        double R0 = 0.0, R1 = 0.0, R2 = 0.0;
        for (int k = 0; k < np; ++k) {
            R0 += cc[k] * e[k];
            R1 += cc[k] * ae[k];
        }
        if (deriv >= 2)
            for (int k = 0; k < np; ++k) R2 += cc[k] * aae[k];
        qmgrid::kern::aoComponent(s.lmn[3 * ci], s.lmn[3 * ci + 1], s.lmn[3 * ci + 2], x, y, z, R0, R1, R2, deriv, cart[c]);
    }
    if (toff < 0) {
        for (int k = 0; k < nder; ++k)
            for (int c = 0; c < nc; ++c) ao[k * plane + (size_t)(first + c) * B + q] = cart[c][k];
        return;
    }
    const double* T = s.T + toff;
    for (int k = 0; k < nder; ++k)
        for (int w = 0; w < width; ++w) {
            double acc = 0.0;
            for (int c = 0; c < nc; ++c) {
                const double t = T[c * width + w];
                if (t != 0.0) acc += t * cart[c][k];
            }
            ao[k * plane + (size_t)(first + w) * B + q] = acc;
        }
}

/// Per point: rho, grad rho (from X = phi P), the functional, and the coefficient vectors of V_xc (mode 0) or of
/// the gradient (mode 1). Writes the whole-grid arrays exw, ecw, ne, used (energy sums) and, mode 1, edens.
///   mode 0: c0 = w v_rho / 2, cj = 2 w v_sigma d_j rho          (Y = c0 phi + sum_j cj d_j phi, buildVxc)
///   mode 1: c0 = w v_rho,     cj = 2 w v_sigma d_j rho if used   (Z, buildXCGradient)
__global__ void densityXCKernel(const double* __restrict__ ao, const double* __restrict__ X, int nbf, int B, int count,
                                int begin, const double* __restrict__ wgt, DevTerms terms, int gga, double rho_min,
                                int mode, double* __restrict__ c0, double* __restrict__ cx, double* __restrict__ cy,
                                double* __restrict__ cz, double* __restrict__ exw, double* __restrict__ ecw,
                                double* __restrict__ ne, double* __restrict__ used, double* __restrict__ edens)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const size_t plane = (size_t)nbf * B;
    double rho = 0.0, gx = 0.0, gy = 0.0, gz = 0.0;
    for (int m = 0; m < nbf; ++m) {
        const size_t i = (size_t)m * B + q;
        const double xv = X[i];
        rho += ao[i] * xv;
        if (gga) {
            gx += ao[plane + i] * xv;
            gy += ao[2 * plane + i] * xv;
            gz += ao[3 * plane + i] * xv;
        }
    }
    gx *= 2.0;
    gy *= 2.0;
    gz *= 2.0;
    const double sigma = gx * gx + gy * gy + gz * gz;
    double ex = 0.0, ec = 0.0, vr = 0.0, vs = 0.0;
    const bool ok = qmxc::kern::evaluatePoint(terms.t, terms.n, gga != 0, rho_min, rho, sigma, ex, ec, vr, vs);
    const int p = begin + q;
    const double w = wgt[p];
    const double wv = w * vr;
    ne[p] = w * rho;
    exw[p] = ok ? w * rho * ex : 0.0;
    ecw[p] = ok ? w * rho * ec : 0.0;
    used[p] = ok ? 1.0 : 0.0;
    const double w2s = 2.0 * w * vs;
    if (mode == 0) {
        c0[q] = 0.5 * wv;
        cx[q] = w2s * gx;
        cy[q] = w2s * gy;
        cz[q] = w2s * gz;
    } else {
        c0[q] = wv;
        cx[q] = ok && gga ? w2s * gx : 0.0;
        cy[q] = ok && gga ? w2s * gy : 0.0;
        cz[q] = ok && gga ? w2s * gz : 0.0;
        edens[p] = ok ? rho * (ex + ec) : 0.0;
    }
}

/// Y = c0 phi (+ cx dx + cy dy + cz dz for a GGA), B x nbf column-major (rows q < count used).
__global__ void combineKernel(const double* __restrict__ ao, int nbf, int B, int count, int gga, const double* __restrict__ c0,
                              const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
                              double* __restrict__ Y)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    const size_t plane = (size_t)nbf * B;
    if (k >= plane) return;
    const int q = (int)(k % B);
    if (q >= count) {
        Y[k] = 0.0;
        return;
    }
    double y = c0[q] * ao[k];
    if (gga) {
        y += cx[q] * ao[plane + k];
        y += cy[q] * ao[2 * plane + k];
        y += cz[q] * ao[3 * plane + k];
    }
    Y[k] = y;
}

/// Meta-GGA version of densityXCKernel (WP11 stage 3): X holds four planes [phi P, d_x phi P, d_y phi P, d_z phi P]
/// (plane = nbf B); tau = 1/2 sum_k sum_m d_k phi (d_k phi P), the functional is qmxc::kern::evaluatePointMGGA.
///   mode 0: c0 = w v_rho / 2, cj = 2 w v_sigma d_j rho, ct = w v_tau / 4   (W += phi^T Y + sum_k d_k^T ct d_k)
///   mode 1: c0 = w v_rho, cj = 2 w v_sigma d_j rho, ct = w v_tau / 2 if used, else zero (buildXCGradientMeta)
/// Claude Generated (Oct 2026).
__global__ void densityXCMetaKernel(const double* __restrict__ ao, const double* __restrict__ X, int nbf, int B, int count,
                                    int begin, const double* __restrict__ wgt, DevTerms terms, double rho_min, int mode,
                                    double* __restrict__ c0, double* __restrict__ cx, double* __restrict__ cy,
                                    double* __restrict__ cz, double* __restrict__ ct, double* __restrict__ exw,
                                    double* __restrict__ ecw, double* __restrict__ ne, double* __restrict__ used,
                                    double* __restrict__ edens)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const size_t plane = (size_t)nbf * B;
    double rho = 0.0, gx = 0.0, gy = 0.0, gz = 0.0, tau = 0.0;
    for (int m = 0; m < nbf; ++m) {
        const size_t i = (size_t)m * B + q;
        const double xv = X[i];
        const double dx = ao[plane + i], dy = ao[2 * plane + i], dz = ao[3 * plane + i];
        rho += ao[i] * xv;
        gx += dx * xv;
        gy += dy * xv;
        gz += dz * xv;
        tau += dx * X[plane + i] + dy * X[2 * plane + i] + dz * X[3 * plane + i];
    }
    gx *= 2.0;
    gy *= 2.0;
    gz *= 2.0;
    tau *= 0.5;
    const double sigma = gx * gx + gy * gy + gz * gz;
    double ex = 0.0, ec = 0.0, vr = 0.0, vs = 0.0, vt = 0.0;
    const bool ok = qmxc::kern::evaluatePointMGGA(terms.t, terms.n, rho_min, rho, sigma, tau, ex, ec, vr, vs, vt);
    const int p = begin + q;
    const double w = wgt[p];
    ne[p] = w * rho;
    exw[p] = ok ? w * rho * ex : 0.0;
    ecw[p] = ok ? w * rho * ec : 0.0;
    used[p] = ok ? 1.0 : 0.0;
    const double w2s = 2.0 * w * vs;
    if (mode == 0) {
        c0[q] = 0.5 * w * vr;
        cx[q] = w2s * gx;
        cy[q] = w2s * gy;
        cz[q] = w2s * gz;
        ct[q] = 0.25 * w * vt;
    } else {
        c0[q] = w * vr;
        cx[q] = ok ? w2s * gx : 0.0;
        cy[q] = ok ? w2s * gy : 0.0;
        cz[q] = ok ? w2s * gz : 0.0;
        ct[q] = ok ? 0.5 * w * vt : 0.0;
        edens[p] = ok ? rho * (ex + ec) : 0.0;
    }
}

/// Meta-GGA V_xc: Y_j = ct o d_j phi for j = x, y, z (planes 1..3 of ao into planes 1..3 of Y, B x nbf column-major
/// each; zero rows beyond count). Claude Generated (Oct 2026).
__global__ void combineTauKernel(const double* __restrict__ ao, int nbf, int B, int count, const double* __restrict__ ct,
                                 double* __restrict__ Y)
{
    const size_t plane = (size_t)nbf * B;
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= 3 * plane) return;
    const int q = (int)(k % B);
    Y[plane + k] = q < count ? ct[q] * ao[plane + k] : 0.0;
}

/// T_k(q, m) of buildXCGradient (k = 0, 1, 2); kMeta (meta-GGA, planes DP_x, DP_y, DP_z) also the tau term
/// wt sum_j d_k d_j phi (DP_j) of buildXCGradientMeta, added last as on the CPU. A template parameter, not a runtime
/// branch: with a runtime branch the LDA/GGA gradient changed in the last bits (1.7e-15 Eh/A, other FMA contraction).
template <bool kMeta>
__device__ inline void gradTerms(const double* ao, const double* ZP, const double* X, size_t plane, size_t i, int gga,
                                 double vx, double vy, double vz, double t[3], const double* DP, double wt)
{
    const double zp = ZP[i];
    t[0] = ao[plane + i] * zp;
    t[1] = ao[2 * plane + i] * zp;
    t[2] = ao[3 * plane + i] * zp;
    if (gga) {
        const double px = X[i];
        const double xx = ao[4 * plane + i], xy = ao[5 * plane + i], xz = ao[6 * plane + i];
        const double yy = ao[7 * plane + i], yz = ao[8 * plane + i], zz = ao[9 * plane + i];
        t[0] += (vx * xx + vy * xy + vz * xz) * px;
        t[1] += (vx * xy + vy * yy + vz * yz) * px;
        t[2] += (vx * xz + vy * yz + vz * zz) * px;
        if constexpr (kMeta) {
            const double p0 = DP[i], p1 = DP[plane + i], p2 = DP[2 * plane + i];
            t[0] += wt * (xx * p0 + xy * p1 + xz * p2);
            t[1] += wt * (xy * p0 + yy * p1 + yz * p2);
            t[2] += wt * (xz * p0 + yz * p1 + zz * p2);
        }
    }
}

/// Point part: spoint[(begin + q) 3 + k] = -2 sum_m T_k(q, m).
template <bool kMeta>
__global__ void gradRowKernel(const double* __restrict__ ao, const double* __restrict__ ZP, const double* __restrict__ X,
                              int nbf, int B, int count, int begin, int gga, const double* __restrict__ cx,
                              const double* __restrict__ cy, const double* __restrict__ cz, double* __restrict__ spoint,
                              const double* __restrict__ DP, const double* __restrict__ ct)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const size_t plane = (size_t)nbf * B;
    const double wt = kMeta ? ct[q] : 0.0;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int m = 0; m < nbf; ++m) {
        double t[3];
        gradTerms<kMeta>(ao, ZP, X, plane, (size_t)m * B + q, gga, cx[q], cy[q], cz[q], t, DP, wt);
        s[0] += t[0];
        s[1] += t[1];
        s[2] += t[2];
    }
    for (int k = 0; k < 3; ++k) spoint[3 * (size_t)(begin + q) + k] = -2.0 * s[k];
}

/// Basis part: gfun[3 m' + k] += sum_q T_k(q, m) (one block per function, fixed tree order); m' = idx[m] with AO
/// screening (compact column m of the group), else m.
template <bool kMeta>
__global__ void gradColKernel(const double* __restrict__ ao, const double* __restrict__ ZP, const double* __restrict__ X,
                              int nbf, int B, int count, int gga, const double* __restrict__ cx, const double* __restrict__ cy,
                              const double* __restrict__ cz, double* __restrict__ gfun, const int* __restrict__ idx,
                              const double* __restrict__ DP, const double* __restrict__ ct)
{
    __shared__ double sh[3][256];
    const int m = blockIdx.x;
    const size_t plane = (size_t)nbf * B;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int q = threadIdx.x; q < count; q += blockDim.x) {
        double t[3];
        gradTerms<kMeta>(ao, ZP, X, plane, (size_t)m * B + q, gga, cx[q], cy[q], cz[q], t, DP, kMeta ? ct[q] : 0.0);
        s[0] += t[0];
        s[1] += t[1];
        s[2] += t[2];
    }
    for (int k = 0; k < 3; ++k) sh[k][threadIdx.x] = s[k];
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o)
            for (int k = 0; k < 3; ++k) sh[k][threadIdx.x] += sh[k][threadIdx.x + o];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        const int mg = idx ? idx[m] : m;
        for (int k = 0; k < 3; ++k) gfun[3 * mg + k] += sh[k][0];
    }
}

/// AO screening: Ps(i, j) = P(idx[i], idx[j]) (ns x ns, column-major; P n x n symmetric).
__global__ void gatherKernel(const double* __restrict__ P, int n, const int* __restrict__ idx, int ns, double* __restrict__ Ps)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (size_t)ns * ns) return;
    const int i = (int)(k % ns), j = (int)(k / ns);
    Ps[k] = P[(size_t)idx[i] + (size_t)idx[j] * n];
}

/// AO screening: W(idx[i], idx[j]) += Ws(i, j); the indices of one group are distinct, groups run one after the
/// other on one stream: no atomics, fixed order.
__global__ void scatterAddKernel(const double* __restrict__ Ws, int ns, const int* __restrict__ idx, int n, double* __restrict__ W)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (size_t)ns * ns) return;
    const int i = (int)(k % ns), j = (int)(k / ns);
    W[(size_t)idx[i] + (size_t)idx[j] * n] += Ws[k];
}

/// Sum of `n` values per chunk of kChunk points (fixed tree order): part[c * narr + j] for arrays a0..a3.
__global__ void chunkSumKernel(const double* __restrict__ a0, const double* __restrict__ a1, const double* __restrict__ a2,
                               const double* __restrict__ a3, int n, double* __restrict__ part)
{
    __shared__ double sh[4][kChunk];
    const int p = blockIdx.x * kChunk + threadIdx.x;
    const double* a[4] = { a0, a1, a2, a3 };
    for (int j = 0; j < 4; ++j) sh[j][threadIdx.x] = p < n ? a[j][p] : 0.0;
    __syncthreads();
    for (int o = kChunk / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o)
            for (int j = 0; j < 4; ++j) sh[j][threadIdx.x] += sh[j][threadIdx.x + o];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        for (int j = 0; j < 4; ++j) part[4 * (size_t)blockIdx.x + j] = sh[j][0];
}

/// out[j] = sum over chunks of part[c * 4 + j] (one block of 256 threads, strided then tree).
__global__ void finalSumKernel(const double* __restrict__ part, int nchunk, int narr, double* __restrict__ out)
{
    __shared__ double sh[256];
    for (int j = 0; j < narr; ++j) {
        double s = 0.0;
        for (int c = threadIdx.x; c < nchunk; c += blockDim.x) s += part[(size_t)c * narr + j];
        sh[threadIdx.x] = s;
        __syncthreads();
        for (int o = blockDim.x / 2; o > 0; o >>= 1) {
            if (threadIdx.x < o) sh[threadIdx.x] += sh[threadIdx.x + o];
            __syncthreads();
        }
        if (threadIdx.x == 0) out[j] = sh[0];
        __syncthreads();
    }
}

/// Point part per atom: gpts[3 A + k] = - sum_{p of A} spoint[3 p + k] (one block per atom, fixed order).
__global__ void atomPointSumKernel(const double* __restrict__ spoint, const int* __restrict__ atom_offset, double* __restrict__ gpts)
{
    __shared__ double sh[3][256];
    const int A = blockIdx.x;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int p = atom_offset[A] + threadIdx.x; p < atom_offset[A + 1]; p += blockDim.x)
        for (int k = 0; k < 3; ++k) s[k] += spoint[3 * (size_t)p + k];
    for (int k = 0; k < 3; ++k) sh[k][threadIdx.x] = s[k];
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o)
            for (int k = 0; k < 3; ++k) sh[k][threadIdx.x] += sh[k][threadIdx.x + o];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        for (int k = 0; k < 3; ++k) gpts[3 * A + k] = -sh[k][0];
}

/// Becke weight derivative (qmgrid::gridWeightGradient, formulas in qm_grid.cpp section WP8): virtual block
/// v = v0 + blockIdx.x of G handles the points v, v + G, v + 2G, ... one at a time; cell functions h, h' per pair
/// (lower triangle, tri(i, j) = i (i - 1)/2 + j for j < i) and "P without one factor" EX (full, prefix/suffix
/// products as on the CPU) in this block's global scratch, the per-point vectors there too (no shared-memory limit
/// on natoms); d w_p / d R_A for A != owner by thread A, the owner from translational invariance; accumulated into
/// slice v of bgrad (natoms x 3), no atomics. S(i, j) = 1/2 (1 - h), S(j, i) = 1/2 (1 + h), dS = -h'/2 both ways,
/// mu(A, D) recomputed as the CPU stores it (j < i: 1/R_ij (d_i - d_j), mu(j, i) = -mu(i, j)): the same bits as the
/// previous kernel with stored matrices. Claude Generated (Oct 2026).
__global__ void weightGradKernel(const double* __restrict__ pts, const double* __restrict__ quadw, const int* __restrict__ patom,
                                 const double* __restrict__ f, int npts, int na, const double* __restrict__ axyz,
                                 const double* __restrict__ invd, const double* __restrict__ aadj, int stratmann, int G, int v0,
                                 size_t per_block, double* __restrict__ scratch, double* __restrict__ bgrad)
{
    __shared__ double zsum;
    const int v = v0 + blockIdx.x;
    if (v >= G) return;
    const size_t nn = (size_t)na * na, nt = (size_t)na * (na - 1) / 2;
    double* EX = scratch + (size_t)blockIdx.x * per_block;
    double* H = EX + nn;
    double* DH = H + nt;
    double* dist = DH + nt;
    double *ux = dist + na, *uy = ux + na, *uz = uy + na, *Pv = uz + na, *dwA = Pv + na;
    double* bg = bgrad + (size_t)v * na * 3;
    for (int p = v; p < npts; p += G) {
        if (f[p] == 0.0) continue;  // uniform over the block
        const int B = patom[p];
        for (int A = threadIdx.x; A < na; A += blockDim.x) {
            const double dx = pts[3 * p] - axyz[3 * A], dy = pts[3 * p + 1] - axyz[3 * A + 1], dz = pts[3 * p + 2] - axyz[3 * A + 2];
            dist[A] = sqrt(dx * dx + dy * dy + dz * dz);
            const double inv = dist[A] > 0.0 ? 1.0 / dist[A] : 0.0;
            ux[A] = dx * inv;
            uy[A] = dy * inv;
            uz[A] = dz * inv;
        }
        __syncthreads();
        for (size_t t = threadIdx.x; t < nt; t += blockDim.x) {
            // t -> (i, j), j < i: i = floor((1 + sqrt(1 + 8 t)) / 2), corrected for rounding
            int i = (int)((1.0 + sqrt(1.0 + 8.0 * (double)t)) * 0.5);
            while ((size_t)i * (i - 1) / 2 > t) --i;
            while ((size_t)(i + 1) * i / 2 <= t) ++i;
            const int j = (int)(t - (size_t)i * (i - 1) / 2);
            const size_t idx = (size_t)i * na + j;
            const double mu = invd[idx] * (dist[i] - dist[j]);
            const double nu = mu + aadj[idx] * (1.0 - mu * mu);
            H[t] = stratmann ? qmgrid::kern::stratmannIterate(nu) : qmgrid::kern::beckeIterate(nu);
            DH[t] = (stratmann ? qmgrid::kern::stratmannIterateDerivative(nu) : qmgrid::kern::beckeIterateDerivative(nu))
                * (1.0 - 2.0 * aadj[idx] * mu);
        }
        __syncthreads();
        for (int i = threadIdx.x; i < na; i += blockDim.x) {
            double* row = EX + (size_t)i * na;
            double pre = 1.0;
            for (int j = 0; j < na; ++j) {
                row[j] = pre;
                const double Sij = j == i ? 1.0
                    : (j < i ? 0.5 * (1.0 - H[(size_t)i * (i - 1) / 2 + j]) : 0.5 * (1.0 + H[(size_t)j * (j - 1) / 2 + i]));
                pre = pre * Sij;
            }
            Pv[i] = pre;
            double suf = 1.0;
            for (int j = na - 1; j >= 0; --j) {
                row[j] = (j == i) ? 0.0 : row[j] * suf;
                const double Sij = j == i ? 1.0
                    : (j < i ? 0.5 * (1.0 - H[(size_t)i * (i - 1) / 2 + j]) : 0.5 * (1.0 + H[(size_t)j * (j - 1) / 2 + i]));
                suf = suf * Sij;
            }
        }
        __syncthreads();
        if (threadIdx.x == 0) {
            double zs = 0.0;
            for (int i = 0; i < na; ++i) zs += Pv[i];
            zsum = zs;
        }
        __syncthreads();
        const double Zs = zsum;
        if (!(Zs > 0.0)) {
            __syncthreads();
            continue;
        }
        const double fw = f[p] * quadw[p];
        const double invZ = 1.0 / Zs, PB = Pv[B];
        for (int A = threadIdx.x; A < na; A += blockDim.x) {
            if (A == B) {
                dwA[3 * A] = dwA[3 * A + 1] = dwA[3 * A + 2] = 0.0;
                continue;
            }
            const double uA[3] = { ux[A], uy[A], uz[A] };
            double sa[3] = { 0.0, 0.0, 0.0 }, dPB[3] = { 0.0, 0.0, 0.0 };
            for (int D = 0; D < na; ++D) {
                if (D == A) continue;
                const size_t ad = (size_t)A * na + D, da = (size_t)D * na + A;
                const double iR = invd[ad];
                const double mu = A > D ? invd[ad] * (dist[A] - dist[D]) : -(invd[da] * (dist[D] - dist[A]));
                const double dS = -0.5 * (A > D ? DH[(size_t)A * (A - 1) / 2 + D] : DH[(size_t)D * (D - 1) / 2 + A]);
                const double cA = EX[ad] * dS, cD = -EX[da] * dS;
                for (int k = 0; k < 3; ++k) {
                    const double dmu = -uA[k] * iR - mu * (axyz[3 * A + k] - axyz[3 * D + k]) * iR * iR;
                    sa[k] += cA * dmu + cD * dmu;
                    if (D == B) dPB[k] += cD * dmu;
                }
            }
            for (int k = 0; k < 3; ++k) {
                const double dw = fw * (dPB[k] * invZ - PB * invZ * invZ * sa[k]);
                bg[3 * A + k] += dw;
                dwA[3 * A + k] = dw;
            }
        }
        __syncthreads();
        if (threadIdx.x < 3) {
            double sb = 0.0;
            for (int A = 0; A < na; ++A) sb += dwA[3 * A + threadIdx.x];
            bg[3 * B + threadIdx.x] -= sb;
        }
        __syncthreads();
    }
}

/// V = W + W^T (n x n), exactly symmetric.
__global__ void symmetrizeKernel(const double* __restrict__ W, double* __restrict__ V, int n)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (size_t)n * n) return;
    const size_t i = k % n, j = k / n;
    V[k] = W[k] + W[j + i * n];
}

// ---------------------------------------------------------------------------------------------------------------
// Double-hybrid gradient parts (WP12 stage 3, Claude Generated, Oct 31 2026): XC kernel response F_xc[D] of the
// Z-vector equations and the P2-directional derivative of the XC gradient, device versions of qmxc::buildFxc and
// qmxc::buildXCResponseGradient (qm_xc_response.cpp, formulas in qm_xc_fxc.h). The second derivatives come from
// qmxc::fxc::evaluatePointFxc (qm_xc_fxc_kernels.h, the host source). Two densities per block: Xa = phi P (the KS
// density) and Xb = phi D (or phi P2), each with the d_k phi planes for a meta-GGA. No FP atomics; the block sizes are
// those of the object (fixed for a fixed setup), so a repeated call gives the same bits.
// ---------------------------------------------------------------------------------------------------------------

/// Per point: rho, grad rho, tau of the density a (KS density P) and b (D or P2), first and second functional
/// derivatives (v, H) at a, the directional derivatives dv_x = sum_y H_xy y_b, and the coefficient vectors:
///   mode 0 (F_xc):     c0 = w dv_rho / 2, cj = 2 w (dv_sigma d_j rho_a + v_sigma d_j rho_b), ct = w dv_tau / 4
///                      (W += phi^T Y + sum_k d_k phi^T diag(ct) d_k phi, V = W + W^T: the X + X^T + tau term of buildFxc)
///   mode 1 (gradient): c0 = w v_rho, cj = 2 w v_sigma d_j rho_a, ct = w v_tau / 2 (Z and wt of the CPU),
///                      d0 = w dv_rho, dj = 2 w (dv_sigma d_j rho_a + v_sigma d_j rho_b), dt = w dv_tau / 2 (dZ, dwt),
///                      edens = v_rho rho_b + v_sigma sigma_b + v_tau tau_b, used = 1 (sigma_b = 2 grad rho_a . grad rho_b)
/// Points below the threshold get zeros (as the CPU skips them). Claude Generated (Oct 31 2026).
__global__ void densityPairKernel(const double* __restrict__ ao, const double* __restrict__ Xa, const double* __restrict__ Xb,
                                  int nbf, int B, int count, int begin, const double* __restrict__ wgt, DevTerms terms, int gga,
                                  int mgga, double rho_min, int mode, double* __restrict__ c0, double* __restrict__ cx,
                                  double* __restrict__ cy, double* __restrict__ cz, double* __restrict__ ct,
                                  double* __restrict__ d0, double* __restrict__ dx, double* __restrict__ dy,
                                  double* __restrict__ dz, double* __restrict__ dt, double* __restrict__ used,
                                  double* __restrict__ edens)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const size_t plane = (size_t)nbf * B;
    double ra = 0.0, rb = 0.0, ga[3] = { 0.0, 0.0, 0.0 }, gb[3] = { 0.0, 0.0, 0.0 }, ta = 0.0, tb = 0.0;
    for (int m = 0; m < nbf; ++m) {
        const size_t i = (size_t)m * B + q;
        const double phi = ao[i], xa = Xa[i], xb = Xb[i];
        ra += phi * xa;
        rb += phi * xb;
        if (gga) {
            for (int k = 0; k < 3; ++k) {
                const double d = ao[(k + 1) * plane + i];
                ga[k] += d * xa;
                gb[k] += d * xb;
                if (mgga) {
                    ta += d * Xa[(k + 1) * plane + i];
                    tb += d * Xb[(k + 1) * plane + i];
                }
            }
        }
    }
    for (int k = 0; k < 3; ++k) {
        ga[k] *= 2.0;
        gb[k] *= 2.0;
    }
    ta *= 0.5;
    tb *= 0.5;
    const double sigma = gga ? ga[0] * ga[0] + ga[1] * ga[1] + ga[2] * ga[2] : 0.0;
    double v[3], H[6];
    const bool ok = qmxc::fxc::evaluatePointFxc(terms.t, terms.n, rho_min, ra, sigma, ta, v, H);
    const int p = begin + q;
    if (!ok) {
        c0[q] = cx[q] = cy[q] = cz[q] = 0.0;
        if (mgga) ct[q] = 0.0;
        if (mode == 1) {
            d0[q] = dx[q] = dy[q] = dz[q] = 0.0;
            if (mgga) dt[q] = 0.0;
            used[p] = 0.0;
            edens[p] = 0.0;
        }
        return;
    }
    const double w = wgt[p];
    const double s1 = gga ? 2.0 * (ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2]) : 0.0;
    const double dvr = H[0] * rb + H[1] * s1 + H[2] * tb;
    const double dvs = H[1] * rb + H[3] * s1 + H[4] * tb;
    const double dvt = H[2] * rb + H[4] * s1 + H[5] * tb;
    double* cj[3] = { cx, cy, cz };
    if (mode == 0) {
        c0[q] = 0.5 * (w * dvr);
        for (int k = 0; k < 3; ++k) cj[k][q] = 2.0 * w * (dvs * ga[k] + v[1] * gb[k]);
        if (mgga) ct[q] = 0.25 * w * dvt;
        return;
    }
    double* dj[3] = { dx, dy, dz };
    c0[q] = w * v[0];
    d0[q] = w * dvr;
    for (int k = 0; k < 3; ++k) {
        cj[k][q] = 2.0 * w * v[1] * ga[k];
        dj[k][q] = 2.0 * w * (dvs * ga[k] + v[1] * gb[k]);
    }
    if (mgga) {
        ct[q] = 0.5 * w * v[2];
        dt[q] = 0.5 * w * dvt;
    }
    used[p] = 1.0;
    edens[p] = v[0] * rb + v[1] * s1 + v[2] * tb;
}

/// T_k(q, m) of buildXCResponseGradient: d_k phi o ZP + sum_j d_k d_j phi o (vj (phi P2) + dvj (phi P))
/// (+ meta-GGA: sum_j d_k d_j phi o (wt (d_j phi P2) + dwt (d_j phi P))), ZP = Z P2 + dZ P. Xa = phi P planes,
/// Xb = phi P2 planes (each [phi, d_x phi, d_y phi, d_z phi] X for a meta-GGA). Claude Generated (Oct 31 2026).
__device__ inline void respTerms(const double* ao, const double* ZP, const double* Xa, const double* Xb, size_t plane,
                                 size_t i, int gga, int mgga, const double vj[3], const double dvj[3], double wt, double dwt,
                                 double t[3])
{
    const double zp = ZP[i];
    t[0] = ao[plane + i] * zp;
    t[1] = ao[2 * plane + i] * zp;
    t[2] = ao[3 * plane + i] * zp;
    if (!gga) return;
    const double a = Xb[i], b = Xa[i];
    const double u0 = vj[0] * a + dvj[0] * b, u1 = vj[1] * a + dvj[1] * b, u2 = vj[2] * a + dvj[2] * b;
    const double xx = ao[4 * plane + i], xy = ao[5 * plane + i], xz = ao[6 * plane + i];
    const double yy = ao[7 * plane + i], yz = ao[8 * plane + i], zz = ao[9 * plane + i];
    t[0] += xx * u0 + xy * u1 + xz * u2;
    t[1] += xy * u0 + yy * u1 + yz * u2;
    t[2] += xz * u0 + yz * u1 + zz * u2;
    if (mgga) {
        const double s0 = wt * Xb[plane + i] + dwt * Xa[plane + i];
        const double s1 = wt * Xb[2 * plane + i] + dwt * Xa[2 * plane + i];
        const double s2 = wt * Xb[3 * plane + i] + dwt * Xa[3 * plane + i];
        t[0] += xx * s0 + xy * s1 + xz * s2;
        t[1] += xy * s0 + yy * s1 + yz * s2;
        t[2] += xz * s0 + yz * s1 + zz * s2;
    }
}

/// Response gradient, point part: spoint[(begin + q) 3 + k] = -2 sum_m T_k(q, m). Claude Generated (Oct 31 2026).
__global__ void respRowKernel(const double* __restrict__ ao, const double* __restrict__ ZP, const double* __restrict__ Xa,
                              const double* __restrict__ Xb, int nbf, int B, int count, int begin, int gga, int mgga,
                              const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
                              const double* __restrict__ dx, const double* __restrict__ dy, const double* __restrict__ dz,
                              const double* __restrict__ ct, const double* __restrict__ dt, double* __restrict__ spoint)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const size_t plane = (size_t)nbf * B;
    const double vj[3] = { cx[q], cy[q], cz[q] }, dvj[3] = { dx[q], dy[q], dz[q] };
    const double wt = mgga ? ct[q] : 0.0, dwt = mgga ? dt[q] : 0.0;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int m = 0; m < nbf; ++m) {
        double t[3];
        respTerms(ao, ZP, Xa, Xb, plane, (size_t)m * B + q, gga, mgga, vj, dvj, wt, dwt, t);
        s[0] += t[0];
        s[1] += t[1];
        s[2] += t[2];
    }
    for (int k = 0; k < 3; ++k) spoint[3 * (size_t)(begin + q) + k] = -2.0 * s[k];
}

/// Response gradient, basis part: gfun[3 m' + k] += sum_q T_k(q, m) (one block per function, fixed tree order).
/// Claude Generated (Oct 31 2026).
__global__ void respColKernel(const double* __restrict__ ao, const double* __restrict__ ZP, const double* __restrict__ Xa,
                              const double* __restrict__ Xb, int nbf, int B, int count, int gga, int mgga,
                              const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
                              const double* __restrict__ dx, const double* __restrict__ dy, const double* __restrict__ dz,
                              const double* __restrict__ ct, const double* __restrict__ dt, double* __restrict__ gfun,
                              const int* __restrict__ idx)
{
    __shared__ double sh[3][256];
    const int m = blockIdx.x;
    const size_t plane = (size_t)nbf * B;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int q = threadIdx.x; q < count; q += blockDim.x) {
        const double vj[3] = { cx[q], cy[q], cz[q] }, dvj[3] = { dx[q], dy[q], dz[q] };
        double t[3];
        respTerms(ao, ZP, Xa, Xb, plane, (size_t)m * B + q, gga, mgga, vj, dvj, mgga ? ct[q] : 0.0, mgga ? dt[q] : 0.0, t);
        s[0] += t[0];
        s[1] += t[1];
        s[2] += t[2];
    }
    for (int k = 0; k < 3; ++k) sh[k][threadIdx.x] = s[k];
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o)
            for (int k = 0; k < 3; ++k) sh[k][threadIdx.x] += sh[k][threadIdx.x + o];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        const int mg = idx ? idx[m] : m;
        for (int k = 0; k < 3; ++k) gfun[3 * mg + k] += sh[k][0];
    }
}

inline unsigned g1(size_t n, int t = 256) { return (unsigned)((n + t - 1) / t); }

class QmGpuXCCuda final : public qmint::QMGpuXC {
public:
    ~QmGpuXCCuda() override
    {
        for (void* p : m_owned) cudaFree(p);
        if (m_cublas) cublasDestroy(m_cublas);
    }

    bool init(const qmint::QMGpuXCSetup& s, int device, std::string& why)
    {
        m_nbf = s.nbf;
        m_npts = s.npts;
        m_na = s.natoms;
        m_gga = s.gga;
        m_rho_min = s.rho_min;
        m_nshell = s.nshell;
        m_stratmann = s.stratmann;
        m_mgga = s.mgga;
        m_fxc = s.fxc;
        if (m_nbf <= 0 || m_npts <= 0 || m_na <= 0) { why = "empty grid or basis"; return false; }
        if (m_mgga && !m_gga) { why = "meta-GGA setup without the GGA flag"; return false; }
        if (s.nterms <= 0 || s.nterms > kMaxTerms) { why = "functional with more than 8 terms"; return false; }
        if (s.max_ncomp > kMaxComp) { why = "shells above g (more than 15 cartesian components) not supported"; return false; }
        for (int i = 0; i < m_nshell; ++i)
            if (s.sh_nprim[i] > kMaxPrim) { why = "shell with more than 32 primitives"; return false; }
        m_terms.n = s.nterms;
        for (int t = 0; t < s.nterms; ++t) {
            m_terms.t[t].kernel = static_cast<qmxc::XCKernel>(s.term_kernel[t]);
            m_terms.t[t].weight = s.term_weight[t];
        }
        // basis function -> atom (host, for the basis part of the gradient)
        m_bf_atom.assign(m_nbf, 0);
        for (int sh = 0; sh < m_nshell; ++sh)
            for (int w = 0; w < s.sh_out_width[sh]; ++w) m_bf_atom[s.sh_out_first[sh] + w] = s.sh_atom[sh];
        if (cublasCreate(&m_cublas) != CUBLAS_STATUS_SUCCESS) { why = "cublasCreate failed"; return false; }
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        m_sms = prop.multiProcessorCount > 0 ? prop.multiProcessorCount : 16;

        const size_t np = (size_t)m_npts, nb = (size_t)m_nbf;
        m_nchunk = (int)((np + kChunk - 1) / kChunk);
        bool ok = up(&m_pts, s.points, 3 * np) && up(&m_w, s.weights, np) && up(&m_qw, s.quad_weights, np)
            && up(&m_patom, s.point_atom, np) && up(&m_aoff, s.atom_offset, (size_t)m_na + 1)
            && up(&m_axyz, s.atom_xyz, 3 * (size_t)m_na) && up(&m_invd, s.inv_dist, (size_t)m_na * m_na)
            && up(&m_aadj, s.a_adjust, (size_t)m_na * m_na)
            && up(const_cast<double**>(&m_sh.xyz), s.sh_xyz, 3 * (size_t)m_nshell)
            && up(const_cast<double**>(&m_sh.amin), s.sh_amin, (size_t)m_nshell)
            && up(const_cast<double**>(&m_sh.exps), s.exps, (size_t)std::max(1, s.nexps))
            && up(const_cast<double**>(&m_sh.coef), s.coef, (size_t)std::max(1, s.ncoef))
            && up(const_cast<double**>(&m_sh.T), s.T, (size_t)std::max(1, s.nT))
            && up(const_cast<int**>(&m_sh.nprim), s.sh_nprim, (size_t)m_nshell)
            && up(const_cast<int**>(&m_sh.prim_off), s.sh_prim_off, (size_t)m_nshell)
            && up(const_cast<int**>(&m_sh.ncomp), s.sh_ncomp, (size_t)m_nshell)
            && up(const_cast<int**>(&m_sh.comp_off), s.sh_comp_off, (size_t)m_nshell)
            && up(const_cast<int**>(&m_sh.coef_off), s.coef_off, (size_t)std::max(1, s.ncomp_total))
            && up(const_cast<int**>(&m_sh.lmn), s.lmn, (size_t)std::max(3, s.nlmn))
            && up(const_cast<int**>(&m_sh.out_first), s.sh_out_first, (size_t)m_nshell)
            && up(const_cast<int**>(&m_sh.out_width), s.sh_out_width, (size_t)m_nshell)
            && up(const_cast<int**>(&m_sh.T_off), s.sh_T_off, (size_t)m_nshell)
            && alloc(&m_exw, np) && alloc(&m_ecw, np) && alloc(&m_ne, np) && alloc(&m_used, np) && alloc(&m_edens, np)
            && alloc(&m_part, 4 * (size_t)m_nchunk) && alloc(&m_red, 8) && alloc(&m_P, nb * nb) && alloc(&m_V, nb * nb) && alloc(&m_W, nb * nb)
            && alloc(&m_gfun, 3 * nb) && alloc(&m_gpts, 3 * (size_t)m_na);
        if (!ok) { why = "device allocation / upload of the grid failed"; return false; }

        // AO screening: merge consecutive host blocks into groups of <= target points, union of their shells
        if (s.screen_nblocks > 0 && !buildGroups(s, why)) return false;

        // Block size from the memory budget: AO values for the gradient (deriv 2 for a GGA: 10 planes, else 4),
        // X, Y, ZP (3 planes; meta-GGA: 3 more each for d_k phi P and ct d_k phi), 5 point vectors (meta-GGA: 6);
        // whole-grid arrays and the n^2 matrices are already allocated.
        size_t free_b = 0, total_b = 0;
        cudaMemGetInfo(&free_b, &total_b);
        double budget = 0.6 * (double)free_b;
        if (s.max_device_mb > 0.0) budget = std::min(budget, s.max_device_mb * 1.0e6 - (double)m_bytes);
        m_nder = m_gga ? 10 : 4;
        const bool screened = !m_groups.empty();
        const size_t ncol = screened ? (size_t)std::max(1, m_sub_max) : nb;  // columns of the per-block planes
        const int nxp = m_mgga ? 4 : 1;  // planes of X: phi P (+ d_k phi P for a meta-GGA)
        // WP12 stage 3: with fxc the second density (nxp planes), the dZ plane and five more point vectors; without it
        // the expression (and with it the block size of LDA/GGA/meta-GGA) is the one of before
        const double per_point = m_fxc ? 8.0 * ((double)ncol * (m_nder + 1 + 2 * nxp + nxp + 1) + (m_mgga ? 9.0 : 8.0) + 5.0)
                                       : 8.0 * ((double)ncol * (m_nder + 1 + 2 * nxp) + (m_mgga ? 9.0 : 8.0));
        long long B = s.block_points > 0 ? s.block_points : (long long)(budget / per_point);
        B = std::min<long long>(B, ((long long)np + 127) / 128 * 128);
        if (s.block_points <= 0) B = std::min<long long>(B, 32768);
        B = B / 32 * 32;
        if (screened) B = (m_group_max + 31) / 32 * 32;  // the leading dimension holds the largest group
        if (B < 32 || ((s.block_points > 0 || screened) && (double)B * per_point > budget)) {
            char buf[200];
            std::snprintf(buf, sizeof buf, "not enough device memory for a block of grid points (budget %.1f MB, %.3f MB per 32 points)",
                          budget / 1.0e6, 32.0 * per_point / 1.0e6);
            why = buf;
            return false;
        }
        m_B = (int)B;
        const size_t plane = (size_t)m_B * ncol;
        ok = alloc(&m_ao, plane * m_nder) && alloc(&m_X, plane * nxp) && alloc(&m_Y, plane * nxp) && alloc(&m_ZP, plane)
            && alloc(&m_c0, (size_t)m_B) && alloc(&m_cx, (size_t)m_B) && alloc(&m_cy, (size_t)m_B) && alloc(&m_cz, (size_t)m_B)
            && alloc(&m_spoint, 3 * np) && (!m_mgga || alloc(&m_ct, (size_t)m_B));
        if (ok && m_fxc)  // WP12 stage 3: second density, dZ and the directional coefficients
            ok = alloc(&m_X2, plane * nxp) && alloc(&m_Y2, plane) && alloc(&m_P2, nb * nb) && alloc(&m_d0, (size_t)m_B)
                && alloc(&m_dx, (size_t)m_B) && alloc(&m_dy, (size_t)m_B) && alloc(&m_dz, (size_t)m_B)
                && (!m_mgga || alloc(&m_dt, (size_t)m_B));
        if (ok && screened) {
            const size_t ms = (size_t)m_sub_max * m_sub_max;
            ok = alloc(&m_Ps, std::max<size_t>(1, ms)) && alloc(&m_Ws, std::max<size_t>(1, ms))
                && (!m_fxc || alloc(&m_P2s, std::max<size_t>(1, ms)))
                && up(&m_g_shells, m_h_shells.data(), m_h_shells.size()) && up(&m_g_cols, m_h_cols.data(), m_h_cols.size())
                && up(&m_g_idx, m_h_idx.data(), m_h_idx.size()) && up(&m_g_bstart, m_h_bstart.data(), m_h_bstart.size())
                && up(&m_g_mask, m_h_mask.data(), m_h_mask.size());
        }
        if (!ok) { why = "device allocation of the grid blocks failed"; return false; }
        return cudaGetLastError() == cudaSuccess;
    }

    std::string name() const override { return m_name; }
    int nbf() const override { return m_nbf; }
    int blockPoints() const override { return m_B; }

    bool vxc(const double* P, double* V, qmint::QMGpuXCScalars& out) override
    {
        const size_t bytes = (size_t)m_nbf * m_nbf * sizeof(double);
        if (cudaMemcpy(m_P, P, bytes, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (!vxcDevice(m_P, m_V, out)) return false;
        return cudaMemcpy(V, m_V, bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
    }

    bool vxcDevice(const double* dP, double* dV, qmint::QMGpuXCScalars& out) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        const int n = m_nbf;
        const int deriv = m_gga ? 1 : 0;
        const double one = 1.0, zero = 0.0;
        if (!m_groups.empty()) {  // AO screening
            if (!clearPointArrays() || cudaMemsetAsync(m_W, 0, (size_t)n * n * sizeof(double), 0) != cudaSuccess) return false;
            for (const Group& g : m_groups) {
                if (g.ns == 0) continue;
                if (!groupDensity(dP, g, deriv, 0)) return false;
                combineKernel<<<g1((size_t)m_B * g.ns), 256>>>(m_ao, g.ns, m_B, g.count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
                if (m_mgga) {
                    if (!metaVxcProduct(g.ns, g.count, m_Ws, false)) return false;
                } else if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, g.ns, g.ns, g.count, &one, m_ao, m_B, m_Y, m_B, &zero,
                                       m_Ws, g.ns) != CUBLAS_STATUS_SUCCESS)
                    return false;
                scatterAddKernel<<<g1((size_t)g.ns * g.ns), 256>>>(m_Ws, g.ns, m_g_idx + g.idx_off, n, m_W);
                profile(3);
            }
        }
        for (int begin = 0; m_groups.empty() && begin < m_npts; begin += m_B) {
            const int count = std::min(m_B, m_npts - begin);
            if (!blockDensity(dP, begin, count, deriv, 0)) return false;
            combineKernel<<<g1((size_t)m_B * n), 256>>>(m_ao, n, m_B, count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
            // W += phi^T Y (cublasDsyr2k was ~40 ms per call for n = 24, k = 29000 on the RTX 5090; DGEMM is not)
            if (m_mgga) {
                if (!metaVxcProduct(n, count, m_W, begin != 0)) return false;
            } else if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, count, &one, m_ao, m_B, m_Y, m_B,
                                   begin == 0 ? &zero : &one, m_W, n) != CUBLAS_STATUS_SUCCESS)
                return false;
            profile(3);
        }
        symmetrizeKernel<<<g1((size_t)n * n), 256>>>(m_W, dV, n);  // V = W + W^T (CPU: V += X + X^T per block)
        if (!energySums(out)) return false;
        profile(4);
        ++m_nvxc;
        m_t_vxc += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return cudaGetLastError() == cudaSuccess;
    }

    bool gradient(const double* P, bool response, double* g_basis, double* g_points, double* g_weights,
                  long long& points_used) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        const int n = m_nbf;
        if (cudaMemcpy(m_P, P, (size_t)n * n * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        cudaMemsetAsync(m_gfun, 0, 3 * (size_t)n * sizeof(double), 0);
        const int deriv = m_gga ? 2 : 1;
        const double one = 1.0, zero = 0.0;
        if (!m_groups.empty()) {  // AO screening
            if (!clearPointArrays()) return false;
            for (const Group& g : m_groups) {
                if (g.ns == 0) continue;
                if (!groupDensity(m_P, g, deriv, 1)) return false;
                combineKernel<<<g1((size_t)m_B * g.ns), 256>>>(m_ao, g.ns, m_B, g.count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
                if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, g.count, g.ns, g.ns, &one, m_Y, m_B, m_Ps, g.ns, &zero, m_ZP,
                                m_B) != CUBLAS_STATUS_SUCCESS)
                    return false;
                const double* DPg = m_X + (size_t)g.ns * m_B;  // meta-GGA: planes d_k phi P
                if (m_mgga) {
                    gradRowKernel<true><<<g1(g.count, 128), 128>>>(m_ao, m_ZP, m_X, g.ns, m_B, g.count, g.begin, m_gga, m_cx, m_cy,
                                                                   m_cz, m_spoint, DPg, m_ct);
                    gradColKernel<true><<<g.ns, 256>>>(m_ao, m_ZP, m_X, g.ns, m_B, g.count, m_gga, m_cx, m_cy, m_cz, m_gfun,
                                                       m_g_idx + g.idx_off, DPg, m_ct);
                } else {
                    gradRowKernel<false><<<g1(g.count, 128), 128>>>(m_ao, m_ZP, m_X, g.ns, m_B, g.count, g.begin, m_gga, m_cx, m_cy,
                                                                    m_cz, m_spoint, nullptr, nullptr);
                    gradColKernel<false><<<g.ns, 256>>>(m_ao, m_ZP, m_X, g.ns, m_B, g.count, m_gga, m_cx, m_cy, m_cz, m_gfun,
                                                        m_g_idx + g.idx_off, nullptr, nullptr);
                }
            }
        }
        for (int begin = 0; m_groups.empty() && begin < m_npts; begin += m_B) {
            const int count = std::min(m_B, m_npts - begin);
            if (!blockDensity(m_P, begin, count, deriv, 1)) return false;
            // Z = w v_rho phi + sum_j v_j d_j phi (zero rows beyond count), ZP = Z P
            combineKernel<<<g1((size_t)m_B * n), 256>>>(m_ao, n, m_B, count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
            if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, n, n, &one, m_Y, m_B, m_P, n, &zero, m_ZP, m_B)
                != CUBLAS_STATUS_SUCCESS)
                return false;
            const double* DPb = m_X + (size_t)n * m_B;  // meta-GGA: planes d_k phi P
            if (m_mgga) {
                gradRowKernel<true><<<g1(count, 128), 128>>>(m_ao, m_ZP, m_X, n, m_B, count, begin, m_gga, m_cx, m_cy, m_cz,
                                                             m_spoint, DPb, m_ct);
                gradColKernel<true><<<n, 256>>>(m_ao, m_ZP, m_X, n, m_B, count, m_gga, m_cx, m_cy, m_cz, m_gfun, nullptr, DPb, m_ct);
            } else {
                gradRowKernel<false><<<g1(count, 128), 128>>>(m_ao, m_ZP, m_X, n, m_B, count, begin, m_gga, m_cx, m_cy, m_cz,
                                                              m_spoint, nullptr, nullptr);
                gradColKernel<false><<<n, 256>>>(m_ao, m_ZP, m_X, n, m_B, count, m_gga, m_cx, m_cy, m_cz, m_gfun, nullptr, nullptr,
                                                 nullptr);
            }
        }
        if (cudaGetLastError() != cudaSuccess) return false;
        qmint::QMGpuXCScalars sc;
        if (!energySums(sc)) return false;
        points_used = sc.points_used;
        std::vector<double> gfun(3 * (size_t)n);
        if (cudaMemcpy(gfun.data(), m_gfun, gfun.size() * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        for (int k = 0; k < 3 * m_na; ++k) g_basis[k] = g_points[k] = g_weights[k] = 0.0;
        for (int m = 0; m < n; ++m)
            for (int k = 0; k < 3; ++k) g_basis[3 * m_bf_atom[m] + k] += -2.0 * gfun[3 * (size_t)m + k];
        const auto t1 = std::chrono::steady_clock::now();
        if (response) {
            atomPointSumKernel<<<m_na, 256>>>(m_spoint, m_aoff, m_gpts);
            if (cudaMemcpy(g_points, m_gpts, 3 * (size_t)m_na * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess)
                return false;
            if (m_na >= 2 && !weightGradientDevice(g_weights)) return false;
        }
        const auto t2 = std::chrono::steady_clock::now();
        m_t_grad += std::chrono::duration<double, std::milli>(t1 - t0).count();
        m_t_wgrad += std::chrono::duration<double, std::milli>(t2 - t1).count();
        ++m_ngrad;
        return cudaGetLastError() == cudaSuccess;
    }

    std::string stats() const override
    {
        char buf[320];
        std::snprintf(buf, sizeof buf,
                      "%d grid points in blocks of %d, %d V_xc builds %.1f ms, %d XC gradients %.1f ms (weight derivative %.1f ms), "
                      "device memory %.1f MB",
                      m_npts, m_B, m_nvxc, m_t_vxc, m_ngrad, m_t_grad, m_t_wgrad, m_bytes / 1.0e6);
        std::string out = buf;
        if (m_fxc) {
            char fx[160];
            std::snprintf(fx, sizeof fx, "; %d F_xc builds %.1f ms, %d XC response gradients %.1f ms", m_nfxc, m_t_fxc, m_nresp,
                          m_t_resp);
            out += fx;
        }
        if (!m_groups.empty()) {
            char sc[200];
            std::snprintf(sc, sizeof sc, "; AO screening: %d groups (<= %d points), %.1f %% of the functions kept on average",
                          (int)m_groups.size(), m_group_max, 100.0 * m_kept_frac);
            out += sc;
        }
        if (m_wgrad_waves > 0) {
            char wg[160];
            std::snprintf(wg, sizeof wg, "; weight derivative: %d virtual blocks in %d wave(s), scratch %.1f MB", m_wgrad_G,
                          m_wgrad_waves, m_wgrad_scratch_mb);
            out += wg;
        }
        if (!m_prof_on) return out;
        char p[200];
        std::snprintf(p, sizeof p, "; V_xc stages (synchronised): AO %.1f ms, X = phi P %.1f ms, density+XC %.1f ms, "
                      "W = phi^T Y %.1f ms, sums %.1f ms", m_prof[0], m_prof[1], m_prof[2], m_prof[3], m_prof[4]);
        return out + p;
    }

    // WP12 stage 3 (Claude Generated, Oct 31 2026): F_xc[D] of the double-hybrid Z-vector, the device version of
    // qmxc::buildFxc. Blocks/groups as vxcDevice (AO deriv 1 for a GGA), two densities per block, then the V_xc
    // assembly of that path (combineKernel, phi^T Y, meta-GGA tau products, V = W + W^T).
    bool fxc(const double* P, const double* D, double* V) override
    {
        if (!m_fxc) return false;
        const auto t0 = std::chrono::steady_clock::now();
        const int n = m_nbf;
        const size_t bytes = (size_t)n * n * sizeof(double);
        if (cudaMemcpy(m_P, P, bytes, cudaMemcpyHostToDevice) != cudaSuccess
            || cudaMemcpy(m_P2, D, bytes, cudaMemcpyHostToDevice) != cudaSuccess)
            return false;
        const int deriv = m_gga ? 1 : 0;
        const double one = 1.0, zero = 0.0;
        if (!m_groups.empty()) {  // AO screening
            if (cudaMemsetAsync(m_W, 0, bytes, 0) != cudaSuccess) return false;
            for (const Group& g : m_groups) {
                if (g.ns == 0) continue;
                if (!groupPair(g, deriv, 0)) return false;
                combineKernel<<<g1((size_t)m_B * g.ns), 256>>>(m_ao, g.ns, m_B, g.count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
                if (m_mgga) {
                    if (!metaVxcProduct(g.ns, g.count, m_Ws, false)) return false;
                } else if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, g.ns, g.ns, g.count, &one, m_ao, m_B, m_Y, m_B, &zero,
                                       m_Ws, g.ns) != CUBLAS_STATUS_SUCCESS)
                    return false;
                scatterAddKernel<<<g1((size_t)g.ns * g.ns), 256>>>(m_Ws, g.ns, m_g_idx + g.idx_off, n, m_W);
            }
        }
        for (int begin = 0; m_groups.empty() && begin < m_npts; begin += m_B) {
            const int count = std::min(m_B, m_npts - begin);
            if (!blockPair(begin, count, deriv, 0)) return false;
            combineKernel<<<g1((size_t)m_B * n), 256>>>(m_ao, n, m_B, count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
            if (m_mgga) {
                if (!metaVxcProduct(n, count, m_W, begin != 0)) return false;
            } else if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, count, &one, m_ao, m_B, m_Y, m_B,
                                   begin == 0 ? &zero : &one, m_W, n) != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        symmetrizeKernel<<<g1((size_t)n * n), 256>>>(m_W, m_V, n);
        if (cudaGetLastError() != cudaSuccess || cudaMemcpy(V, m_V, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        ++m_nfxc;
        m_t_fxc += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    // WP12 stage 3 (Claude Generated, Oct 31 2026): the P2-directional XC gradient, the device version of
    // qmxc::buildXCResponseGradient plus the Becke weight part of its energy density (weightGradientDevice).
    bool responseGradient(const double* P, const double* P2, bool response, double* g_basis, double* g_points,
                          double* g_weights, long long& points_used) override
    {
        if (!m_fxc) return false;
        const auto t0 = std::chrono::steady_clock::now();
        const int n = m_nbf;
        const size_t bytes = (size_t)n * n * sizeof(double);
        if (cudaMemcpy(m_P, P, bytes, cudaMemcpyHostToDevice) != cudaSuccess
            || cudaMemcpy(m_P2, P2, bytes, cudaMemcpyHostToDevice) != cudaSuccess)
            return false;
        if (cudaMemsetAsync(m_gfun, 0, 3 * (size_t)n * sizeof(double), 0) != cudaSuccess || !clearPointArrays()) return false;
        const int deriv = m_gga ? 2 : 1;
        const double one = 1.0, zero = 0.0;
        // ZP = Z P2 + dZ P (Z in plane 0 of m_Y, dZ in m_Y2), then the T_k row and column sums
        auto contract = [&](int ns, int count, int begin, const double* Pa, const double* Pb, const int* idx) {
            combineKernel<<<g1((size_t)m_B * ns), 256>>>(m_ao, ns, m_B, count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
            combineKernel<<<g1((size_t)m_B * ns), 256>>>(m_ao, ns, m_B, count, m_gga, m_d0, m_dx, m_dy, m_dz, m_Y2);
            if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, ns, ns, &one, m_Y, m_B, Pb, ns, &zero, m_ZP, m_B)
                    != CUBLAS_STATUS_SUCCESS
                || cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, ns, ns, &one, m_Y2, m_B, Pa, ns, &one, m_ZP, m_B)
                    != CUBLAS_STATUS_SUCCESS)
                return false;
            respRowKernel<<<g1(count, 128), 128>>>(m_ao, m_ZP, m_X, m_X2, ns, m_B, count, begin, m_gga, m_mgga, m_cx, m_cy, m_cz,
                                                   m_dx, m_dy, m_dz, m_ct, m_dt, m_spoint);
            respColKernel<<<ns, 256>>>(m_ao, m_ZP, m_X, m_X2, ns, m_B, count, m_gga, m_mgga, m_cx, m_cy, m_cz, m_dx, m_dy, m_dz,
                                       m_ct, m_dt, m_gfun, idx);
            return cudaGetLastError() == cudaSuccess;
        };
        if (!m_groups.empty()) {  // AO screening
            for (const Group& g : m_groups) {
                if (g.ns == 0) continue;
                if (!groupPair(g, deriv, 1) || !contract(g.ns, g.count, g.begin, m_Ps, m_P2s, m_g_idx + g.idx_off)) return false;
            }
        }
        for (int begin = 0; m_groups.empty() && begin < m_npts; begin += m_B) {
            const int count = std::min(m_B, m_npts - begin);
            if (!blockPair(begin, count, deriv, 1) || !contract(n, count, begin, m_P, m_P2, nullptr)) return false;
        }
        qmint::QMGpuXCScalars sc;
        if (!energySums(sc)) return false;  // exw/ecw/ne are zero here; only the point count is read
        points_used = sc.points_used;
        std::vector<double> gfun(3 * (size_t)n);
        if (cudaMemcpy(gfun.data(), m_gfun, gfun.size() * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        for (int k = 0; k < 3 * m_na; ++k) g_basis[k] = g_points[k] = g_weights[k] = 0.0;
        for (int m = 0; m < n; ++m)
            for (int k = 0; k < 3; ++k) g_basis[3 * m_bf_atom[m] + k] += -2.0 * gfun[3 * (size_t)m + k];
        if (response) {
            atomPointSumKernel<<<m_na, 256>>>(m_spoint, m_aoff, m_gpts);
            if (cudaMemcpy(g_points, m_gpts, 3 * (size_t)m_na * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess)
                return false;
            if (m_na >= 2 && !weightGradientDevice(g_weights)) return false;
        }
        ++m_nresp;
        m_t_resp += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return cudaGetLastError() == cudaSuccess;
    }

    bool weightGradient(const double* f, double* g) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < 3 * m_na; ++k) g[k] = 0.0;
        if (cudaMemcpy(m_edens, f, (size_t)m_npts * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        const bool ok = m_na < 2 || weightGradientDevice(g);
        m_t_wgrad += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return ok && cudaGetLastError() == cudaSuccess;
    }

private:
    /// AO screening: one streamed group of consecutive host blocks, its kept shells (with compact first columns)
    /// and functions.
    struct Group {
        int begin = 0, count = 0, ns = 0, nsh = 0, sh_off = 0, idx_off = 0;
        int nblk = 0, bstart_off = 0, mask_off = 0;  // host blocks of the group, their starts, shell masks
    };

    /// Groups of <= target points from the host blocks (union of their shell lists, ascending). Claude Generated.
    bool buildGroups(const qmint::QMGpuXCSetup& s, std::string& why)
    {
        int target = s.block_points > 0 ? s.block_points : kScreenGroupPoints;
        for (int b = 0; b < s.screen_nblocks; ++b) target = std::max(target, s.screen_count[b]);  // one host block fits
        std::vector<char> mark(m_nshell, 0);
        std::vector<int> width(m_nshell), first(m_nshell);
        for (int i = 0; i < m_nshell; ++i) {
            width[i] = s.sh_out_width[i];
            first[i] = s.sh_out_first[i];
        }
        long long kept = 0, total = 0;
        int b = 0;
        while (b < s.screen_nblocks) {
            Group g;
            g.begin = s.screen_begin[b];
            int e = b;
            while (e < s.screen_nblocks && g.count + s.screen_count[e] <= target && s.screen_begin[e] == g.begin + g.count) {
                g.count += s.screen_count[e];
                for (int k = s.screen_shell_off[e]; k < s.screen_shell_off[e + 1]; ++k) mark[s.screen_shells[k]] = 1;
                ++e;
            }
            if (e == b) {
                why = "AO screening: a block of the host screen is larger than the device group";
                return false;
            }
            g.sh_off = (int)m_h_shells.size();
            g.idx_off = (int)m_h_idx.size();
            std::vector<int> gsh;
            for (int sh = 0; sh < m_nshell; ++sh) {
                if (!mark[sh]) continue;
                mark[sh] = 0;
                gsh.push_back(sh);
                m_h_shells.push_back(sh);
                m_h_cols.push_back(g.ns);
                for (int w = 0; w < width[sh]; ++w) m_h_idx.push_back(first[sh] + w);
                g.ns += width[sh];
                ++g.nsh;
            }
            // per host block of the group: which of the group's shells it keeps (the AO kernel writes zeros for the
            // others, so every point sees exactly the functions of its own block, as on the CPU)
            g.nblk = e - b;
            g.bstart_off = (int)m_h_bstart.size();
            g.mask_off = (int)m_h_mask.size();
            for (int k = b; k < e; ++k) {
                m_h_bstart.push_back(s.screen_begin[k] - g.begin);
                std::vector<char> in(m_nshell, 0);
                for (int t = s.screen_shell_off[k]; t < s.screen_shell_off[k + 1]; ++t) in[s.screen_shells[t]] = 1;
                for (int sh : gsh) m_h_mask.push_back((unsigned char)in[sh]);
            }
            m_h_bstart.push_back(g.count);
            m_group_max = std::max(m_group_max, g.count);
            m_sub_max = std::max(m_sub_max, g.ns);
            kept += (long long)g.ns * g.count;
            total += (long long)m_nbf * g.count;
            m_groups.push_back(g);
            b = e;
        }
        m_kept_frac = total > 0 ? (double)kept / (double)total : 1.0;
        if (m_h_shells.empty()) m_h_shells.push_back(0), m_h_cols.push_back(0);  // non-empty uploads
        if (m_h_idx.empty()) m_h_idx.push_back(0);
        if (m_h_mask.empty()) m_h_mask.push_back(0);
        return true;
    }

    /// Per-point arrays of the whole grid to zero (AO screening: groups without a shell write nothing).
    bool clearPointArrays()
    {
        const size_t b = (size_t)m_npts * sizeof(double);
        return cudaMemsetAsync(m_exw, 0, b, 0) == cudaSuccess && cudaMemsetAsync(m_ecw, 0, b, 0) == cudaSuccess
            && cudaMemsetAsync(m_ne, 0, b, 0) == cudaSuccess && cudaMemsetAsync(m_used, 0, b, 0) == cudaSuccess
            && cudaMemsetAsync(m_edens, 0, b, 0) == cudaSuccess && cudaMemsetAsync(m_spoint, 0, 3 * b, 0) == cudaSuccess;
    }

    /// AO screening: AO values of the kept shells of group g, P_SS (m_Ps), X = phi_S P_SS, per-point functional.
    bool groupDensity(const double* dP, const Group& g, int deriv, int mode)
    {
        const double one = 1.0, zero = 0.0;
        profile(-1);
        aoKernel<<<dim3(g1(g.count, 128), g.nsh), 128>>>(m_sh, m_pts, g.begin, g.count, deriv, g.ns, m_B, m_ao,
                                                         m_g_shells + g.sh_off, m_g_cols + g.sh_off, m_g_bstart + g.bstart_off,
                                                         m_g_mask + g.mask_off, g.nblk);
        gatherKernel<<<g1((size_t)g.ns * g.ns), 256>>>(dP, m_nbf, m_g_idx + g.idx_off, g.ns, m_Ps);
        if (cudaGetLastError() != cudaSuccess) return false;
        profile(0);
        if (m_mgga) return metaDensity(m_Ps, g.ns, g.begin, g.count, mode);
        if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, g.count, g.ns, g.ns, &one, m_ao, m_B, m_Ps, g.ns, &zero, m_X, m_B)
            != CUBLAS_STATUS_SUCCESS)
            return false;
        profile(1);
        densityXCKernel<<<g1(g.count, 128), 128>>>(m_ao, m_X, g.ns, m_B, g.count, g.begin, m_w, m_terms, m_gga, m_rho_min, mode,
                                                   m_c0, m_cx, m_cy, m_cz, m_exw, m_ecw, m_ne, m_used, m_edens);
        profile(2);
        return cudaGetLastError() == cudaSuccess;
    }

    /// AO values of the block, X = phi P, per-point functional and coefficient vectors (mode 0: V_xc, 1: gradient).
    bool blockDensity(const double* dP, int begin, int count, int deriv, int mode)
    {
        const int n = m_nbf;
        const double one = 1.0, zero = 0.0;
        profile(-1);
        aoKernel<<<dim3(g1(count, 128), m_nshell), 128>>>(m_sh, m_pts, begin, count, deriv, n, m_B, m_ao, nullptr, nullptr,
                                                          nullptr, nullptr, 0);
        if (cudaGetLastError() != cudaSuccess) return false;
        profile(0);
        if (m_mgga) return metaDensity(dP, n, begin, count, mode);
        if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, n, n, &one, m_ao, m_B, dP, n, &zero, m_X, m_B)
            != CUBLAS_STATUS_SUCCESS)
            return false;
        profile(1);
        densityXCKernel<<<g1(count, 128), 128>>>(m_ao, m_X, n, m_B, count, begin, m_w, m_terms, m_gga, m_rho_min, mode, m_c0,
                                                 m_cx, m_cy, m_cz, m_exw, m_ecw, m_ne, m_used, m_edens);
        profile(2);
        return cudaGetLastError() == cudaSuccess;
    }

    /// WP12 stage 3: Xa = phi Pa, Xb = phi Pb (meta-GGA: the four planes each, one strided-batched GEMM per density),
    /// then densityPairKernel. ns columns (n or the kept functions of a group). Claude Generated (Oct 31 2026).
    bool pairDensity(const double* Pa, const double* Pb, int ns, int begin, int count, int mode)
    {
        const double one = 1.0, zero = 0.0;
        const long long plane = (long long)ns * m_B;
        const int nxp = m_mgga ? 4 : 1;
        if (cublasDgemmStridedBatched(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, ns, ns, &one, m_ao, m_B, plane, Pa, ns, 0, &zero,
                                      m_X, m_B, plane, nxp) != CUBLAS_STATUS_SUCCESS
            || cublasDgemmStridedBatched(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, ns, ns, &one, m_ao, m_B, plane, Pb, ns, 0,
                                         &zero, m_X2, m_B, plane, nxp) != CUBLAS_STATUS_SUCCESS)
            return false;
        densityPairKernel<<<g1(count, 128), 128>>>(m_ao, m_X, m_X2, ns, m_B, count, begin, m_w, m_terms, m_gga, m_mgga, m_rho_min,
                                                   mode, m_c0, m_cx, m_cy, m_cz, m_ct, m_d0, m_dx, m_dy, m_dz, m_dt, m_used,
                                                   m_edens);
        return cudaGetLastError() == cudaSuccess;
    }

    /// WP12 stage 3: AO values of group g, P and the second density (m_P, m_P2) gathered to the kept functions, the
    /// pair densities. Claude Generated (Oct 31 2026).
    bool groupPair(const Group& g, int deriv, int mode)
    {
        aoKernel<<<dim3(g1(g.count, 128), g.nsh), 128>>>(m_sh, m_pts, g.begin, g.count, deriv, g.ns, m_B, m_ao,
                                                         m_g_shells + g.sh_off, m_g_cols + g.sh_off, m_g_bstart + g.bstart_off,
                                                         m_g_mask + g.mask_off, g.nblk);
        gatherKernel<<<g1((size_t)g.ns * g.ns), 256>>>(m_P, m_nbf, m_g_idx + g.idx_off, g.ns, m_Ps);
        gatherKernel<<<g1((size_t)g.ns * g.ns), 256>>>(m_P2, m_nbf, m_g_idx + g.idx_off, g.ns, m_P2s);
        if (cudaGetLastError() != cudaSuccess) return false;
        return pairDensity(m_Ps, m_P2s, g.ns, g.begin, g.count, mode);
    }

    /// WP12 stage 3: AO values of the block [begin, begin + count) (all functions) and the pair densities of m_P, m_P2.
    bool blockPair(int begin, int count, int deriv, int mode)
    {
        aoKernel<<<dim3(g1(count, 128), m_nshell), 128>>>(m_sh, m_pts, begin, count, deriv, m_nbf, m_B, m_ao, nullptr, nullptr,
                                                          nullptr, nullptr, 0);
        if (cudaGetLastError() != cudaSuccess) return false;
        return pairDensity(m_P, m_P2, m_nbf, begin, count, mode);
    }

    /// Meta-GGA (WP11 stage 3): [X, DP_x, DP_y, DP_z] = [phi, d_x phi, d_y phi, d_z phi] P (one strided-batched
    /// GEMM over the four consecutive AO planes, P (ns x ns) shared), then the per-point functional with tau.
    /// ns = columns of the block (n, or the kept functions of a screened group). Claude Generated (Oct 2026).
    bool metaDensity(const double* P, int ns, int begin, int count, int mode)
    {
        const double one = 1.0, zero = 0.0;
        const long long plane = (long long)ns * m_B;
        if (cublasDgemmStridedBatched(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, ns, ns, &one, m_ao, m_B, plane, P, ns, 0, &zero,
                                      m_X, m_B, plane, 4) != CUBLAS_STATUS_SUCCESS)
            return false;
        profile(1);
        densityXCMetaKernel<<<g1(count, 128), 128>>>(m_ao, m_X, ns, m_B, count, begin, m_w, m_terms, m_rho_min, mode, m_c0, m_cx,
                                                     m_cy, m_cz, m_ct, m_exw, m_ecw, m_ne, m_used, m_edens);
        profile(2);
        return cudaGetLastError() == cudaSuccess;
    }

    /// Meta-GGA V_xc of one block: Wb (ns x ns) (+)= phi^T Y + sum_k (d_k phi)^T diag(ct) (d_k phi), ct = w v_tau / 4
    /// (densityXCMetaKernel mode 0; Y in plane 0 of m_Y from combineKernel, the planes ct d_k phi from
    /// combineTauKernel): four DGEMMs accumulated in fixed order on one stream. A strided-batched GEMM into scratch
    /// plus a sum was 3x slower on the RTX 5090 (benzene 150x974: 9.2 s instead of 3.1 s for 13 builds).
    /// Claude Generated (Oct 2026).
    bool metaVxcProduct(int ns, int count, double* Wb, bool accumulate)
    {
        const double one = 1.0, zero = 0.0;
        const size_t plane = (size_t)ns * m_B;
        combineTauKernel<<<g1(3 * plane), 256>>>(m_ao, ns, m_B, count, m_ct, m_Y);
        for (int k = 0; k < 4; ++k)
            if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, ns, ns, count, &one, m_ao + k * plane, m_B, m_Y + k * plane, m_B,
                            (k == 0 && !accumulate) ? &zero : &one, Wb, ns) != CUBLAS_STATUS_SUCCESS)
                return false;
        return cudaGetLastError() == cudaSuccess;
    }

    /// CURCUMA_QM_GPU_XC_PROFILE=1: synchronise after each stage and accumulate its time (stage -1 = start mark);
    /// printed by stats(). Diagnostic only, the numbers do not change.
    void profile(int stage)
    {
        if (!m_prof_on) return;
        cudaDeviceSynchronize();
        const auto t = std::chrono::steady_clock::now();
        if (stage >= 0 && stage < 5) m_prof[stage] += std::chrono::duration<double, std::milli>(t - m_prof_t).count();
        m_prof_t = t;
    }

    /// E_x, E_c, N, points used over the whole grid (fixed chunks), synchronous download.
    bool energySums(qmint::QMGpuXCScalars& out)
    {
        chunkSumKernel<<<m_nchunk, kChunk>>>(m_exw, m_ecw, m_ne, m_used, m_npts, m_part);
        finalSumKernel<<<1, 256>>>(m_part, m_nchunk, 4, m_red);
        double r[4];
        if (cudaGetLastError() != cudaSuccess || cudaMemcpy(r, m_red, sizeof r, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        out.ex = r[0];
        out.ec = r[1];
        out.nelec = r[2];
        out.points_used = (long long)(r[3] + 0.5);
        return true;
    }

    /// Becke weight derivative contracted with m_edens (energy density of the last gradient(), or the host f of
    /// weightGradient()); adds into g. G = 4 virtual blocks per SM (fixed: the summation order does not depend on
    /// the free memory), launched in waves of as many blocks as fit 50 % of the free memory (2 natoms^2 + 8 natoms
    /// doubles each); fails only when not even one block fits. Claude Generated (Oct 2026: no natoms limit).
    bool weightGradientDevice(double* g)
    {
        const size_t nn = (size_t)m_na * m_na, nt = (size_t)m_na * (m_na - 1) / 2;
        const size_t per_block = nn + 2 * nt + 8 * (size_t)m_na;
        const int G = std::max(1, std::min(4 * m_sms, m_npts));
        size_t free_b = 0, total_b = 0;
        cudaMemGetInfo(&free_b, &total_b);
        const double bgrad_b = 8.0 * (double)G * m_na * 3;
        const double avail = 0.5 * (double)free_b - bgrad_b;
        const int Gw = (int)std::min<double>(G, std::floor(avail / (8.0 * (double)per_block)));
        if (Gw < 1) return false;
        double *scratch = nullptr, *bg = nullptr;
        if (cudaMalloc(&scratch, (size_t)Gw * per_block * sizeof(double)) != cudaSuccess) return false;
        if (cudaMalloc(&bg, (size_t)G * m_na * 3 * sizeof(double)) != cudaSuccess) {
            cudaFree(scratch);
            return false;
        }
        cudaMemset(bg, 0, (size_t)G * m_na * 3 * sizeof(double));
        int waves = 0;
        for (int v0 = 0; v0 < G; v0 += Gw, ++waves)
            weightGradKernel<<<std::min(Gw, G - v0), 128>>>(m_pts, m_qw, m_patom, m_edens, m_npts, m_na, m_axyz, m_invd, m_aadj,
                                                            m_stratmann, G, v0, per_block, scratch, bg);
        std::vector<double> h((size_t)G * m_na * 3);
        const bool ok = cudaGetLastError() == cudaSuccess
            && cudaMemcpy(h.data(), bg, h.size() * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
        cudaFree(scratch);
        cudaFree(bg);
        m_wgrad_G = G;
        m_wgrad_waves = waves;
        m_wgrad_scratch_mb = 8.0 * ((double)Gw * per_block + (double)G * m_na * 3) / 1.0e6;
        if (!ok) return false;
        for (int b = 0; b < G; ++b)
            for (int k = 0; k < 3 * m_na; ++k) g[k] += h[(size_t)b * m_na * 3 + k];
        return true;
    }

    template <typename T>
    bool alloc(T** p, size_t count)
    {
        if (cudaMalloc(p, count * sizeof(T)) != cudaSuccess) { *p = nullptr; return false; }
        m_owned.push_back(*p);
        m_bytes += count * sizeof(T);
        return true;
    }
    template <typename T>
    bool up(T** p, const T* host, size_t count)
    {
        if (!host) return alloc(p, count) && cudaMemset(*p, 0, count * sizeof(T)) == cudaSuccess;
        return alloc(p, count) && cudaMemcpy(*p, host, count * sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess;
    }

    cublasHandle_t m_cublas = nullptr;
    std::vector<void*> m_owned;
    std::vector<int> m_bf_atom;
    DevShells m_sh{};
    DevTerms m_terms;
    int m_nbf = 0, m_npts = 0, m_na = 0, m_gga = 0, m_nshell = 0, m_stratmann = 0, m_B = 0, m_nder = 4, m_nchunk = 0, m_sms = 16;
    double m_rho_min = 0.0;
    double *m_pts = nullptr, *m_w = nullptr, *m_qw = nullptr, *m_axyz = nullptr, *m_invd = nullptr, *m_aadj = nullptr;
    int *m_patom = nullptr, *m_aoff = nullptr;
    double *m_exw = nullptr, *m_ecw = nullptr, *m_ne = nullptr, *m_used = nullptr, *m_edens = nullptr, *m_part = nullptr,
           *m_red = nullptr, *m_P = nullptr, *m_V = nullptr, *m_W = nullptr, *m_gfun = nullptr, *m_gpts = nullptr, *m_ao = nullptr,
           *m_X = nullptr, *m_Y = nullptr, *m_ZP = nullptr, *m_c0 = nullptr, *m_cx = nullptr, *m_cy = nullptr, *m_cz = nullptr,
           *m_spoint = nullptr;
    int m_mgga = 0;          // meta-GGA (WP11 stage 3): m_X has 4 planes, m_ct the per-point v_tau coefficient
    double* m_ct = nullptr;
    // double-hybrid gradient (WP12 stage 3): second density planes, dZ, D or P2 (full and gathered), directional
    // coefficients
    int m_fxc = 0;
    double *m_X2 = nullptr, *m_Y2 = nullptr, *m_P2 = nullptr, *m_P2s = nullptr, *m_d0 = nullptr, *m_dx = nullptr,
           *m_dy = nullptr, *m_dz = nullptr, *m_dt = nullptr;
    int m_nfxc = 0, m_nresp = 0;
    double m_t_fxc = 0.0, m_t_resp = 0.0;
    size_t m_bytes = 0;
    // AO screening (groups on the host, shell/column/function lists on the device)
    std::vector<Group> m_groups;
    std::vector<int> m_h_shells, m_h_cols, m_h_idx, m_h_bstart;
    std::vector<unsigned char> m_h_mask;
    int *m_g_shells = nullptr, *m_g_cols = nullptr, *m_g_idx = nullptr, *m_g_bstart = nullptr;
    unsigned char* m_g_mask = nullptr;
    double *m_Ps = nullptr, *m_Ws = nullptr;
    int m_group_max = 0, m_sub_max = 0;
    double m_kept_frac = 1.0;
    int m_wgrad_G = 0, m_wgrad_waves = 0;
    double m_wgrad_scratch_mb = 0.0;
    int m_nvxc = 0, m_ngrad = 0;
    double m_t_vxc = 0.0, m_t_grad = 0.0, m_t_wgrad = 0.0;
    std::string m_name = "cuda";
    bool m_prof_on = std::getenv("CURCUMA_QM_GPU_XC_PROFILE") != nullptr;
    double m_prof[5] = { 0, 0, 0, 0, 0 };
    std::chrono::steady_clock::time_point m_prof_t = std::chrono::steady_clock::now();
};

void setWhyXC(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

}  // namespace

// Device XC object (QMGpuXC*) for one geometry, or nullptr with a reason. Claude Generated (Oct 2026, GPU KS-DFT).
// Renamed from _create_qm_xc with the changed vtable (weightGradient) and the appended screening fields, then from
// _create_qm_xc_screen with the appended meta-GGA flag (WP11 stage 3), then from _create_qm_xc_mgga with the appended
// fxc flag and the vtable entries fxc/responseGradient (WP12 stage 3): an old host never hands a short setup to it.
extern "C" void* curcuma_cuda_create_qm_xc_fxc(const qmint::QMGpuXCSetup* s, int device, char* why, int whylen)
{
    if (!s) { setWhyXC(why, whylen, "no setup"); return nullptr; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyXC(why, whylen, "cudaSetDevice failed"); return nullptr; }
    auto x = std::make_unique<QmGpuXCCuda>();
    std::string err;
    if (!x->init(*s, device, err)) {
        cudaGetLastError();
        setWhyXC(why, whylen, err);
        return nullptr;
    }
    return static_cast<qmint::QMGpuXC*>(x.release());
}
