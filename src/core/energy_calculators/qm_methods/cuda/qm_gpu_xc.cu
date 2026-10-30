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

#include "../qm_gpu_xc_backend.h"
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
constexpr int kMaxComp = 10;      // cartesian components per shell (L <= 3)
constexpr int kChunk = 256;       // points per energy-reduction chunk (fixed, independent of B)
constexpr int kMaxWeightAtoms = 700;  // shared memory of the weight-derivative kernel: 8 natoms doubles

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
__global__ void aoKernel(DevShells s, const double* __restrict__ pts, int begin, int count, int deriv, int nbf, int B,
                         double* __restrict__ ao)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const int sh = blockIdx.y;
    const int nder = deriv >= 2 ? 10 : (deriv >= 1 ? 4 : 1);
    const int p = begin + q;
    const double x = pts[3 * p] - s.xyz[3 * sh], y = pts[3 * p + 1] - s.xyz[3 * sh + 1], z = pts[3 * p + 2] - s.xyz[3 * sh + 2];
    const double r2 = x * x + y * y + z * z;
    const int nc = s.ncomp[sh], first = s.out_first[sh], width = s.out_width[sh], toff = s.T_off[sh];
    const size_t plane = (size_t)nbf * B;
    if (s.amin[sh] * r2 > 80.0) {
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

/// T_k(q, m) of buildXCGradient (k = 0, 1, 2).
__device__ inline void gradTerms(const double* ao, const double* ZP, const double* X, size_t plane, size_t i, int gga,
                                 double vx, double vy, double vz, double t[3])
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
    }
}

/// Point part: spoint[(begin + q) 3 + k] = -2 sum_m T_k(q, m).
__global__ void gradRowKernel(const double* __restrict__ ao, const double* __restrict__ ZP, const double* __restrict__ X,
                              int nbf, int B, int count, int begin, int gga, const double* __restrict__ cx,
                              const double* __restrict__ cy, const double* __restrict__ cz, double* __restrict__ spoint)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= count) return;
    const size_t plane = (size_t)nbf * B;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int m = 0; m < nbf; ++m) {
        double t[3];
        gradTerms(ao, ZP, X, plane, (size_t)m * B + q, gga, cx[q], cy[q], cz[q], t);
        s[0] += t[0];
        s[1] += t[1];
        s[2] += t[2];
    }
    for (int k = 0; k < 3; ++k) spoint[3 * (size_t)(begin + q) + k] = -2.0 * s[k];
}

/// Basis part: gfun[3 m + k] += sum_q T_k(q, m) (one block per function, fixed tree order).
__global__ void gradColKernel(const double* __restrict__ ao, const double* __restrict__ ZP, const double* __restrict__ X,
                              int nbf, int B, int count, int gga, const double* __restrict__ cx, const double* __restrict__ cy,
                              const double* __restrict__ cz, double* __restrict__ gfun)
{
    __shared__ double sh[3][256];
    const int m = blockIdx.x;
    const size_t plane = (size_t)nbf * B;
    double s[3] = { 0.0, 0.0, 0.0 };
    for (int q = threadIdx.x; q < count; q += blockDim.x) {
        double t[3];
        gradTerms(ao, ZP, X, plane, (size_t)m * B + q, gga, cx[q], cy[q], cz[q], t);
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
    if (threadIdx.x == 0)
        for (int k = 0; k < 3; ++k) gfun[3 * m + k] += sh[k][0];
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

/// Becke weight derivative (qmgrid::gridWeightGradient, formulas in qm_grid.cpp section WP8): one thread block per
/// point at a time (points blockIdx.x, + gridDim.x, ...), cell functions in per-block global scratch (S, dS, mu,
/// "P without one factor" from prefix/suffix products as on the CPU), d w_p / d R_A for A != owner by thread A, the
/// owner from translational invariance; accumulated into this block's slice of bgrad (natoms x 3), no atomics.
__global__ void weightGradKernel(const double* __restrict__ pts, const double* __restrict__ quadw, const int* __restrict__ patom,
                                 const double* __restrict__ f, int npts, int na, const double* __restrict__ axyz,
                                 const double* __restrict__ invd, const double* __restrict__ aadj, int stratmann,
                                 double* __restrict__ scratch, double* __restrict__ bgrad)
{
    extern __shared__ double sm[];
    double *dist = sm, *ux = sm + na, *uy = sm + 2 * na, *uz = sm + 3 * na, *Pv = sm + 4 * na, *dwA = sm + 5 * na;
    double* misc = sm + 8 * na;
    const size_t nn = (size_t)na * na;
    double* S = scratch + (size_t)blockIdx.x * 4 * nn;
    double* dS = S + nn;
    double* MU = dS + nn;
    double* EX = MU + nn;
    double* bg = bgrad + (size_t)blockIdx.x * na * 3;
    for (int p = blockIdx.x; p < npts; p += gridDim.x) {
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
        for (size_t idx = threadIdx.x; idx < nn; idx += blockDim.x) {
            const int i = (int)(idx / na), j = (int)(idx % na);
            if (i == j) {
                S[idx] = 1.0;
                dS[idx] = 0.0;
                MU[idx] = 0.0;
                continue;
            }
            if (j > i) continue;
            const double mu = invd[idx] * (dist[i] - dist[j]);
            const double nu = mu + aadj[idx] * (1.0 - mu * mu);
            const double h = stratmann ? qmgrid::kern::stratmannIterate(nu) : qmgrid::kern::beckeIterate(nu);
            const double dh = (stratmann ? qmgrid::kern::stratmannIterateDerivative(nu) : qmgrid::kern::beckeIterateDerivative(nu))
                * (1.0 - 2.0 * aadj[idx] * mu);
            const size_t ji = (size_t)j * na + i;
            S[idx] = 0.5 * (1.0 - h);
            S[ji] = 0.5 * (1.0 + h);
            dS[idx] = -0.5 * dh;
            dS[ji] = -0.5 * dh;
            MU[idx] = mu;
            MU[ji] = -mu;
        }
        __syncthreads();
        for (int i = threadIdx.x; i < na; i += blockDim.x) {
            double* row = EX + (size_t)i * na;
            const double* Si = S + (size_t)i * na;
            double pre = 1.0;
            for (int j = 0; j < na; ++j) {
                row[j] = pre;
                pre = pre * (j == i ? 1.0 : Si[j]);
            }
            Pv[i] = pre;
            double suf = 1.0;
            for (int j = na - 1; j >= 0; --j) {
                row[j] = (j == i) ? 0.0 : row[j] * suf;
                suf = suf * (j == i ? 1.0 : Si[j]);
            }
        }
        __syncthreads();
        if (threadIdx.x == 0) {
            double zs = 0.0;
            for (int i = 0; i < na; ++i) zs += Pv[i];
            misc[0] = zs;
        }
        __syncthreads();
        const double Zs = misc[0];
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
            double s[3] = { 0.0, 0.0, 0.0 }, dPB[3] = { 0.0, 0.0, 0.0 };
            for (int D = 0; D < na; ++D) {
                if (D == A) continue;
                const size_t ad = (size_t)A * na + D, da = (size_t)D * na + A;
                const double iR = invd[ad], mu = MU[ad];
                const double cA = EX[ad] * dS[ad], cD = -EX[da] * dS[da];
                for (int k = 0; k < 3; ++k) {
                    const double dmu = -uA[k] * iR - mu * (axyz[3 * A + k] - axyz[3 * D + k]) * iR * iR;
                    s[k] += cA * dmu + cD * dmu;
                    if (D == B) dPB[k] += cD * dmu;
                }
            }
            for (int k = 0; k < 3; ++k) {
                const double dw = fw * (dPB[k] * invZ - PB * invZ * invZ * s[k]);
                bg[3 * A + k] += dw;
                dwA[3 * A + k] = dw;
            }
        }
        __syncthreads();
        if (threadIdx.x < 3) {
            double s = 0.0;
            for (int A = 0; A < na; ++A) s += dwA[3 * A + threadIdx.x];
            bg[3 * B + threadIdx.x] -= s;
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
        if (m_nbf <= 0 || m_npts <= 0 || m_na <= 0) { why = "empty grid or basis"; return false; }
        if (s.nterms <= 0 || s.nterms > kMaxTerms) { why = "functional with more than 8 terms"; return false; }
        if (s.max_ncomp > kMaxComp) { why = "shells above f (more than 10 cartesian components) not supported"; return false; }
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

        // Block size from the memory budget: AO values for the gradient (deriv 2 for a GGA: 10 planes, else 4),
        // X, Y, ZP (3 planes), 5 point vectors; whole-grid arrays and the n^2 matrices are already allocated.
        size_t free_b = 0, total_b = 0;
        cudaMemGetInfo(&free_b, &total_b);
        double budget = 0.6 * (double)free_b;
        if (s.max_device_mb > 0.0) budget = std::min(budget, s.max_device_mb * 1.0e6 - (double)m_bytes);
        m_nder = m_gga ? 10 : 4;
        const double per_point = 8.0 * ((double)nb * (m_nder + 3) + 8.0);
        long long B = s.block_points > 0 ? s.block_points : (long long)(budget / per_point);
        B = std::min<long long>(B, ((long long)np + 127) / 128 * 128);
        if (s.block_points <= 0) B = std::min<long long>(B, 32768);
        B = B / 32 * 32;
        if (B < 32 || (s.block_points > 0 && (double)B * per_point > budget)) {
            char buf[200];
            std::snprintf(buf, sizeof buf, "not enough device memory for a block of grid points (budget %.1f MB, %.3f MB per 32 points)",
                          budget / 1.0e6, 32.0 * per_point / 1.0e6);
            why = buf;
            return false;
        }
        m_B = (int)B;
        const size_t plane = (size_t)m_B * nb;
        ok = alloc(&m_ao, plane * m_nder) && alloc(&m_X, plane) && alloc(&m_Y, plane) && alloc(&m_ZP, plane)
            && alloc(&m_c0, (size_t)m_B) && alloc(&m_cx, (size_t)m_B) && alloc(&m_cy, (size_t)m_B) && alloc(&m_cz, (size_t)m_B)
            && alloc(&m_spoint, 3 * np);
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
        for (int begin = 0; begin < m_npts; begin += m_B) {
            const int count = std::min(m_B, m_npts - begin);
            if (!blockDensity(dP, begin, count, deriv, 0)) return false;
            combineKernel<<<g1((size_t)m_B * n), 256>>>(m_ao, n, m_B, count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
            // W += phi^T Y (cublasDsyr2k was ~40 ms per call for n = 24, k = 29000 on the RTX 5090; DGEMM is not)
            if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, count, &one, m_ao, m_B, m_Y, m_B, begin == 0 ? &zero : &one,
                            m_W, n) != CUBLAS_STATUS_SUCCESS)
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
        for (int begin = 0; begin < m_npts; begin += m_B) {
            const int count = std::min(m_B, m_npts - begin);
            if (!blockDensity(m_P, begin, count, deriv, 1)) return false;
            // Z = w v_rho phi + sum_j v_j d_j phi (zero rows beyond count), ZP = Z P
            combineKernel<<<g1((size_t)m_B * n), 256>>>(m_ao, n, m_B, count, m_gga, m_c0, m_cx, m_cy, m_cz, m_Y);
            if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, n, n, &one, m_Y, m_B, m_P, n, &zero, m_ZP, m_B)
                != CUBLAS_STATUS_SUCCESS)
                return false;
            gradRowKernel<<<g1(count, 128), 128>>>(m_ao, m_ZP, m_X, n, m_B, count, begin, m_gga, m_cx, m_cy, m_cz, m_spoint);
            gradColKernel<<<n, 256>>>(m_ao, m_ZP, m_X, n, m_B, count, m_gga, m_cx, m_cy, m_cz, m_gfun);
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
            if (m_na >= 2 && !weightGradient(g_weights)) return false;
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
        if (!m_prof_on) return buf;
        char p[200];
        std::snprintf(p, sizeof p, "; V_xc stages (synchronised): AO %.1f ms, X = phi P %.1f ms, density+XC %.1f ms, "
                      "W = phi^T Y %.1f ms, sums %.1f ms", m_prof[0], m_prof[1], m_prof[2], m_prof[3], m_prof[4]);
        return std::string(buf) + p;
    }

private:
    /// AO values of the block, X = phi P, per-point functional and coefficient vectors (mode 0: V_xc, 1: gradient).
    bool blockDensity(const double* dP, int begin, int count, int deriv, int mode)
    {
        const int n = m_nbf;
        const double one = 1.0, zero = 0.0;
        profile(-1);
        aoKernel<<<dim3(g1(count, 128), m_nshell), 128>>>(m_sh, m_pts, begin, count, deriv, n, m_B, m_ao);
        if (cudaGetLastError() != cudaSuccess) return false;
        profile(0);
        if (cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, count, n, n, &one, m_ao, m_B, dP, n, &zero, m_X, m_B)
            != CUBLAS_STATUS_SUCCESS)
            return false;
        profile(1);
        densityXCKernel<<<g1(count, 128), 128>>>(m_ao, m_X, n, m_B, count, begin, m_w, m_terms, m_gga, m_rho_min, mode, m_c0,
                                                 m_cx, m_cy, m_cz, m_exw, m_ecw, m_ne, m_used, m_edens);
        profile(2);
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

    /// Becke weight derivative contracted with the energy density of the last gradient() (m_edens).
    bool weightGradient(double* g)
    {
        if (m_na > kMaxWeightAtoms) return false;
        const size_t nn = (size_t)m_na * m_na;
        size_t free_b = 0, total_b = 0;
        cudaMemGetInfo(&free_b, &total_b);
        const double per_block = 8.0 * (4.0 * nn + 3.0 * m_na);
        int G = std::min(4 * m_sms, (int)std::max(1.0, 0.5 * (double)free_b / per_block));
        G = std::max(1, std::min(G, m_npts));
        double *scratch = nullptr, *bg = nullptr;
        if (cudaMalloc(&scratch, (size_t)G * 4 * nn * sizeof(double)) != cudaSuccess) return false;
        if (cudaMalloc(&bg, (size_t)G * m_na * 3 * sizeof(double)) != cudaSuccess) {
            cudaFree(scratch);
            return false;
        }
        cudaMemset(bg, 0, (size_t)G * m_na * 3 * sizeof(double));
        const size_t shmem = (8 * (size_t)m_na + 2) * sizeof(double);
        weightGradKernel<<<G, 128, shmem>>>(m_pts, m_qw, m_patom, m_edens, m_npts, m_na, m_axyz, m_invd, m_aadj, m_stratmann,
                                            scratch, bg);
        std::vector<double> h((size_t)G * m_na * 3);
        const bool ok = cudaGetLastError() == cudaSuccess
            && cudaMemcpy(h.data(), bg, h.size() * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
        cudaFree(scratch);
        cudaFree(bg);
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
    size_t m_bytes = 0;
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
extern "C" void* curcuma_cuda_create_qm_xc(const qmint::QMGpuXCSetup* s, int device, char* why, int whylen)
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
