/*
 * <Device Ewald sum of the GFN-FF Gaussian charges (CUDA)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (Nov 2026, GPU_PERIODIC_GAPS A2, stages S3a/S3b). See ewald_gpu.h for what is computed;
 * the formulas are those of ewald.h/.cpp (B0/B2/N1), referenced per kernel below.
 */

#include "ewald_gpu.h"

#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {

constexpr unsigned FULL_MASK = 0xffffffffu;
constexpr int KGROUP = 8;           ///< k vectors per structure-factor group (same n_x, n_y, consecutive n_z)
constexpr int KPASS_ATOMS = 128;    ///< atoms per block of the per-atom reciprocal pass
constexpr double TWO_OVER_SQRT_PI = 1.1283791670955126;  ///< 2/sqrt(pi)
constexpr double SQRT_PI = 1.7724538509055159;

__device__ __forceinline__ double ew_mic(double d, double L, double inv_L)
{
    double s = d * inv_L;  // PBCUtils::minimumImageOrtho
    s -= rint(s);
    return L * s;
}

__device__ __forceinline__ float ew_micf(float d, float L, float inv_L)
{
    float s = d * inv_L;
    s -= rintf(s);
    return L * s;
}

/// EwaldSum pairKernel (ewald.cpp) for Gaussian charges: K = erfc(a r)/r - erfc(g r)/r (g r < 6), dK/dr
__device__ __forceinline__ double ew_pair_kernel(double r, double alpha, double alp_sum, double& dK)
{
    const double erfc_ar = erfc(alpha * r);
    const double inv_r = 1.0 / r;
    dK = -erfc_ar * inv_r * inv_r - (2.0 * alpha / SQRT_PI) * exp(-alpha * alpha * r * r) * inv_r;
    double K = erfc_ar * inv_r;
    const double g = 1.0 / sqrt(alp_sum);
    const double x = g * r;
    if (x < curcuma::ewald::GAUSSIAN_CORRECTION_XMAX) {
        const double erfc_g = erfc(x);
        K -= erfc_g * inv_r;
        dK += erfc_g * inv_r * inv_r + g * TWO_OVER_SQRT_PI * exp(-x * x) * inv_r;
    }
    return K;
}

struct CellArgs {
    double L[3], inv_L[3];
    float Lf[3], inv_Lf[3];
    double rc2;
    float rcf2;   ///< (r_cut + 0.05 Bohr)^2: FP32 pre-screen
};

/// Float copy of the positions wrapped into [0, L) (pre-screen only)
__global__ void k_ew_wrapf(int N, const double* __restrict__ pos, CellArgs c, float* __restrict__ pf)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= 3 * N) return;
    const int a = t % 3;
    double s = pos[t] * c.inv_L[a];
    s -= floor(s);
    pf[t] = static_cast<float>(c.L[a] * s);
}

/// FP32 candidate test (wrapped float positions, margin in rcf2). The exact FP64 test r^2 < r_cut^2 is made per entry
/// in k_ew_kernel: a candidate outside r_cut gets K = dK = 0, which is the same sum as leaving it out.
__device__ __forceinline__ bool ew_candidate(int j, const float* __restrict__ pf, const CellArgs& c, float fx, float fy, float fz)
{
    const float ex = ew_micf(pf[3 * j] - fx, c.Lf[0], c.inv_Lf[0]);
    const float ey = ew_micf(pf[3 * j + 1] - fy, c.Lf[1], c.inv_Lf[1]);
    const float ez = ew_micf(pf[3 * j + 2] - fz, c.Lf[2], c.inv_Lf[2]);
    return ex * ex + ey * ey + ez * ez < c.rcf2;
}

/// One warp per atom i: number of candidates j != i (count[N] = 0 for the scan)
__global__ void k_ew_count(int N, const float* __restrict__ pf, CellArgs c, int* __restrict__ count)
{
    const int i = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    const int lane = threadIdx.x & 31;
    if (i > N) return;
    if (i == N) { if (lane == 0) count[N] = 0; return; }
    const float fx = pf[3 * i], fy = pf[3 * i + 1], fz = pf[3 * i + 2];
    int n = 0;
    for (int j = lane; j < N; j += 32)
        if (j != i && ew_candidate(j, pf, c, fx, fy, fz)) ++n;
    for (int o = 16; o > 0; o >>= 1) n += __shfl_down_sync(FULL_MASK, n, o);
    if (lane == 0) count[i] = n;
}

/// One warp per atom i: candidates in ascending j, owner atom per entry
__global__ void k_ew_fill(int N, const float* __restrict__ pf, CellArgs c, const int* __restrict__ start, int* __restrict__ nbr,
    int* __restrict__ owner)
{
    const int i = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    const int lane = threadIdx.x & 31;
    if (i >= N) return;
    const float fx = pf[3 * i], fy = pf[3 * i + 1], fz = pf[3 * i + 2];
    int off = start[i];
    for (int base = 0; base < N; base += 32) {
        const int j = base + lane;
        const bool hit = (j < N && j != i) && ew_candidate(j, pf, c, fx, fy, fz);
        const unsigned m = __ballot_sync(FULL_MASK, hit);
        if (hit) {
            const int e = off + __popc(m & ((1u << lane) - 1u));
            nbr[e] = j;
            owner[e] = i;
        }
        off += __popc(m);
    }
}

/// One thread per entry: exact minimum-image test r^2 < r_cut^2 (FP64, PBCUtils::minimumImageOrtho) and the kernel
__global__ void k_ew_kernel(long long M, const double* __restrict__ pos, const double* __restrict__ alp, CellArgs c, double alpha,
    const int* __restrict__ owner, const int* __restrict__ nbr, double* __restrict__ K, double* __restrict__ dK)
{
    const long long e = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (e >= M) return;
    const int i = owner[e], j = nbr[e];
    const double dx = ew_mic(pos[3 * j] - pos[3 * i], c.L[0], c.inv_L[0]);
    const double dy = ew_mic(pos[3 * j + 1] - pos[3 * i + 1], c.L[1], c.inv_L[1]);
    const double dz = ew_mic(pos[3 * j + 2] - pos[3 * i + 2], c.L[2], c.inv_L[2]);
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (!(r2 < c.rc2)) {
        K[e] = 0.0;
        dK[e] = 0.0;
        return;
    }
    double d;
    K[e] = ew_pair_kernel(sqrt(r2), alpha, alp[i] + alp[j], d);
    dK[e] = d;
}

struct PhaseArgs {
    double* re[3];
    double* im[3];
    int nmax[3];
};

/// cos/sin(2 pi n s), s = wrapped fractional coordinate, n = 0..nmax_a; layout [n * N + i] (EwaldPhaseTable)
__global__ void k_ew_phase(int N, const double* __restrict__ pos, CellArgs c, PhaseArgs t)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= 3 * N) return;
    const int a = idx / N, i = idx % N;
    double s = pos[3 * i + a] * c.inv_L[a];
    s -= floor(s);
    for (int n = 0; n <= t.nmax[a]; ++n) {
        double sn, cs;
        sincospi(2.0 * n * s, &sn, &cs);
        t.re[a][static_cast<size_t>(n) * N + i] = cs;
        t.im[a][static_cast<size_t>(n) * N + i] = sn;
    }
}

struct KArgs {
    const int4* groups;   ///< {first k, count, n_x, n_y}
    const int* knz;       ///< n_z per k
    const double* kw;     ///< kx ky kz weight per k (weight doubled for -k, EwaldSum::kVectorTable)
    int ngroups;
    int nk;
};

/// e_i(k) = exp(i k.r_i) from the phase table: exy for (n_x, n_y), times ez(n_z) (EwaldSum::reciprocalRange signs)
__device__ __forceinline__ void ew_exy(const PhaseArgs& t, int N, int i, int nx, int ny, double& exr, double& exi)
{
    const int ax = nx >= 0 ? nx : -nx, ay = ny >= 0 ? ny : -ny;
    const double xr = t.re[0][static_cast<size_t>(ax) * N + i];
    const double b = (nx >= 0 ? 1.0 : -1.0) * t.im[0][static_cast<size_t>(ax) * N + i];
    const double yr = t.re[1][static_cast<size_t>(ay) * N + i];
    const double d = (ny >= 0 ? 1.0 : -1.0) * t.im[1][static_cast<size_t>(ay) * N + i];
    exr = xr * yr - b * d;
    exi = xr * d + b * yr;
}

/// One block (256 threads) per k group: S(k) = sum_i x_i e_i(k), fixed-order warp/block tree. Per atom the products
/// x exy are formed once per group, then each k costs four FMA (sign of n_z by a block-uniform branch). Also writes
/// Sw(k) = 2 w_k S(k), the factor of the per-atom passes.
__global__ void k_ew_sf(int N, KArgs k, PhaseArgs t, const double* __restrict__ x, double* __restrict__ S, double* __restrict__ Sw)
{
    const int4 g = k.groups[blockIdx.x];
    double ar[KGROUP], ai[KGROUP];
    int nz[KGROUP];
#pragma unroll
    for (int c = 0; c < KGROUP; ++c) {
        ar[c] = ai[c] = 0.0;
        nz[c] = (c < g.y) ? k.knz[g.x + c] : 0;
    }
    for (int i = threadIdx.x; i < N; i += blockDim.x) {
        double exr, exi;
        ew_exy(t, N, i, g.z, g.w, exr, exi);
        const double pr = x[i] * exr, pi = x[i] * exi;
#pragma unroll
        for (int c = 0; c < KGROUP; ++c) {
            if (c < g.y) {
                const int az = nz[c] >= 0 ? nz[c] : -nz[c];
                const double zr = t.re[2][static_cast<size_t>(az) * N + i];
                const double zi = t.im[2][static_cast<size_t>(az) * N + i];
                if (nz[c] >= 0) {   // e = exy (zr + i zi)
                    ar[c] += pr * zr - pi * zi;
                    ai[c] += pr * zi + pi * zr;
                } else {            // e = exy (zr - i zi)
                    ar[c] += pr * zr + pi * zi;
                    ai[c] += pi * zr - pr * zi;
                }
            }
        }
    }
    __shared__ double sh[8][2 * KGROUP];
    const int warp = threadIdx.x / 32, lane = threadIdx.x & 31;
#pragma unroll
    for (int c = 0; c < KGROUP; ++c) {
        double vr = ar[c], vi = ai[c];
        for (int o = 16; o > 0; o >>= 1) {
            vr += __shfl_down_sync(FULL_MASK, vr, o);
            vi += __shfl_down_sync(FULL_MASK, vi, o);
        }
        if (lane == 0) { sh[warp][2 * c] = vr; sh[warp][2 * c + 1] = vi; }
    }
    __syncthreads();
    const int nw = blockDim.x / 32;
    if (threadIdx.x < 2 * g.y) {
        double s = 0.0;
        for (int w = 0; w < nw; ++w) s += sh[w][threadIdx.x];
        const int kk = g.x + threadIdx.x / 2;
        S[2 * g.x + threadIdx.x] = s;  // [2 k] = Re, [2 k + 1] = Im
        Sw[2 * g.x + threadIdx.x] = 2.0 * k.kw[4 * kk + 3] * s;
    }
}

/// Per-atom reciprocal pass over the groups of one slice, partial results per slice (summed in slice order later).
///  potentials: phi_i += Re(conj(e_i) Sw) = exr U + exi V, U = sum_c (zr Swr + zi Swi), V = sum_c (zr Swi - zi Swr)
///  gradient:   g_i  += q_i Im(conj(e_i) Sw) k = q_i (exr P - exi O) k with P = zr Swi - zi Swr, O = zi Swi + zr Swr
///  (e = exy ez, ez = zr + i zi with zi carrying the sign of n_z). (k_x, k_y) are fixed per group, k_z per k.
template <bool GRAD>
__global__ void k_ew_kpass(int N, KArgs k, PhaseArgs t, const double* __restrict__ Sw, const double* __restrict__ q,
    int nslice, double* __restrict__ part)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int sl = blockIdx.y;
    if (i >= N) return;
    const int gb = static_cast<int>(static_cast<long long>(k.ngroups) * sl / nslice);
    const int ge = static_cast<int>(static_cast<long long>(k.ngroups) * (sl + 1) / nslice);
    double p0 = 0.0, p1 = 0.0, p2 = 0.0;
    for (int gi = gb; gi < ge; ++gi) {
        const int4 g = k.groups[gi];
        double exr, exi;
        ew_exy(t, N, i, g.z, g.w, exr, exi);
        double U = 0.0, V = 0.0, Tz = 0.0;
        for (int c = 0; c < g.y; ++c) {
            const int kk = g.x + c;
            const int nz = k.knz[kk];
            const int az = nz >= 0 ? nz : -nz;
            const double zr = t.re[2][static_cast<size_t>(az) * N + i];
            const double zi = (nz >= 0 ? 1.0 : -1.0) * t.im[2][static_cast<size_t>(az) * N + i];
            const double Sr = Sw[2 * kk], Si = Sw[2 * kk + 1];
            if (GRAD) {
                const double P = zr * Si - zi * Sr, O = zi * Si + zr * Sr;
                const double tk = exr * P - exi * O;
                U += tk;                    // sum over the group (k_x, k_y common)
                Tz += tk * k.kw[4 * kk + 2];
            } else {
                U += zr * Sr + zi * Si;
                V += zr * Si - zi * Sr;
            }
        }
        if (GRAD) {
            p0 += U * k.kw[4 * g.x];
            p1 += U * k.kw[4 * g.x + 1];
            p2 += Tz;
        } else {
            p0 += exr * U + exi * V;
        }
    }
    const size_t o = static_cast<size_t>(sl) * N + i;
    if (GRAD) {
        const double qi = q[i];
        part[3 * o] = qi * p0;
        part[3 * o + 1] = qi * p1;
        part[3 * o + 2] = qi * p2;
    } else {
        part[o] = p0;
    }
}

/// One warp per atom: y_i = sum_j K_ij x_j + sum_slices recip part (potentials, the pcg operator)
__global__ void k_ew_real_pot(int N, const int* __restrict__ start, const int* __restrict__ nbr, const double* __restrict__ K,
    const double* __restrict__ x, int nslice, const double* __restrict__ part, double* __restrict__ y)
{
    const int i = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    const int lane = threadIdx.x & 31;
    if (i >= N) return;
    double s = 0.0;
    for (int e = start[i] + lane; e < start[i + 1]; e += 32) s += K[e] * x[nbr[e]];
    for (int o = 16; o > 0; o >>= 1) s += __shfl_down_sync(FULL_MASK, s, o);
    if (lane == 0) {
        double r = 0.0;
        for (int sl = 0; sl < nslice; ++sl) r += part[static_cast<size_t>(sl) * N + i];
        y[i] = s + r;
    }
}

/// One warp per atom: real-space energy share 1/2 q_i sum_j K q_j, gradient -sum_j q_i q_j dK/r d_ij (d = r_j - r_i),
/// virial share 1/2 sum_j (-q_i q_j dK/r) d_a d_b (EwaldSum::realSpaceRangeCached / realSpaceVirialRange, each pair
/// seen from both atoms); the gradient adds the reciprocal slices.
template <bool GRAD, bool VIR>
__global__ void k_ew_real_eg(int N, const double* __restrict__ pos, CellArgs c, const int* __restrict__ start,
    const int* __restrict__ nbr, const double* __restrict__ K, const double* __restrict__ dK, const double* __restrict__ q,
    int nslice, const double* __restrict__ part, double* __restrict__ e_atom, double* __restrict__ gx,
    double* __restrict__ gy, double* __restrict__ gz, double* __restrict__ v_atom)
{
    const int i = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    const int lane = threadIdx.x & 31;
    if (i >= N) return;
    const double qi = q[i];
    const double pix = pos[3 * i], piy = pos[3 * i + 1], piz = pos[3 * i + 2];
    double e = 0.0, g0 = 0.0, g1 = 0.0, g2 = 0.0;
    double v[6] = { 0, 0, 0, 0, 0, 0 };  // xx yy zz xy xz yz
    for (int ee = start[i] + lane; ee < start[i + 1]; ee += 32) {
        const int j = nbr[ee];
        const double qj = q[j];
        e += K[ee] * qj;
        if (GRAD) {
            const double dx = ew_mic(pos[3 * j] - pix, c.L[0], c.inv_L[0]);
            const double dy = ew_mic(pos[3 * j + 1] - piy, c.L[1], c.inv_L[1]);
            const double dz = ew_mic(pos[3 * j + 2] - piz, c.L[2], c.inv_L[2]);
            const double r = sqrt(dx * dx + dy * dy + dz * dz);
            const double f = qi * qj * dK[ee] / r;
            g0 -= f * dx;
            g1 -= f * dy;
            g2 -= f * dz;
            if (VIR) {
                const double h = -0.5 * f;
                v[0] += h * dx * dx; v[1] += h * dy * dy; v[2] += h * dz * dz;
                v[3] += h * dx * dy; v[4] += h * dx * dz; v[5] += h * dy * dz;
            }
        }
    }
    for (int o = 16; o > 0; o >>= 1) {
        e += __shfl_down_sync(FULL_MASK, e, o);
        if (GRAD) {
            g0 += __shfl_down_sync(FULL_MASK, g0, o);
            g1 += __shfl_down_sync(FULL_MASK, g1, o);
            g2 += __shfl_down_sync(FULL_MASK, g2, o);
        }
        if (VIR)
            for (int a = 0; a < 6; ++a) v[a] += __shfl_down_sync(FULL_MASK, v[a], o);
    }
    if (lane == 0) {
        e_atom[i] = 0.5 * qi * e;
        if (GRAD) {
            double r0 = 0.0, r1 = 0.0, r2 = 0.0;
            for (int sl = 0; sl < nslice; ++sl) {
                const size_t o = 3 * (static_cast<size_t>(sl) * N + i);
                r0 += part[o];
                r1 += part[o + 1];
                r2 += part[o + 2];
            }
            gx[i] = g0 + r0;
            gy[i] = g1 + r1;
            gz[i] = g2 + r2;
        }
        if (VIR)
            for (int a = 0; a < 6; ++a) v_atom[6 * static_cast<size_t>(i) + a] = v[a];
    }
}

/// One block of 1024 threads: sum of `nval` values per item over n items in fixed order (strided + tree)
template <int NV>
__global__ void k_ew_sum_items(int n, const double* __restrict__ vals, double* __restrict__ out)
{
    double s[NV];
#pragma unroll
    for (int a = 0; a < NV; ++a) s[a] = 0.0;
    for (int it = threadIdx.x; it < n; it += blockDim.x)
#pragma unroll
        for (int a = 0; a < NV; ++a) s[a] += vals[static_cast<size_t>(it) * NV + a];
    __shared__ double sh[32][NV];
    const int warp = threadIdx.x / 32, lane = threadIdx.x & 31;
#pragma unroll
    for (int a = 0; a < NV; ++a) {
        double v = s[a];
        for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(FULL_MASK, v, o);
        if (lane == 0) sh[warp][a] = v;
    }
    __syncthreads();
    if (threadIdx.x < NV) {
        double t = 0.0;
        for (int w = 0; w < static_cast<int>(blockDim.x / 32); ++w) t += sh[w][threadIdx.x];
        out[threadIdx.x] = t;
    }
}

/// Per k: E_k = w |S|^2 and E_k [delta_ab - 2 k_a k_b (1/k^2 + 1/(4 a^2))] (EwaldSum::reciprocalVirialRange)
__global__ void k_ew_kitems(int nk, const double* __restrict__ kw, const double* __restrict__ S, double alpha,
    double* __restrict__ items)
{
    const int kk = blockIdx.x * blockDim.x + threadIdx.x;
    if (kk >= nk) return;
    const double kx = kw[4 * kk], ky = kw[4 * kk + 1], kz = kw[4 * kk + 2];
    const double Sr = S[2 * kk], Si = S[2 * kk + 1];
    const double ek = kw[4 * kk + 3] * (Sr * Sr + Si * Si);
    const double k2 = kx * kx + ky * ky + kz * kz;
    const double cc = 2.0 * (1.0 / k2 + 1.0 / (4.0 * alpha * alpha));
    double* o = items + 7 * static_cast<size_t>(kk);
    o[0] = ek;
    o[1] = ek * (1.0 - cc * kx * kx);
    o[2] = ek * (1.0 - cc * ky * ky);
    o[3] = ek * (1.0 - cc * kz * kz);
    o[4] = -ek * cc * kx * ky;
    o[5] = -ek * cc * kx * kz;
    o[6] = -ek * cc * ky * kz;
}

template <class T>
struct DevBuf {
    T* p = nullptr;
    size_t n = 0;
    cudaError_t ensure(size_t cnt)
    {
        if (cnt <= n && p) return cudaSuccess;
        if (p) cudaFree(p);
        p = nullptr;
        n = 0;
        const cudaError_t e = cudaMalloc(&p, std::max<size_t>(cnt, 1) * sizeof(T));
        if (e == cudaSuccess) n = cnt;
        return e;
    }
    ~DevBuf()
    {
        if (p) cudaFree(p);
    }
};

template <class T>
struct HostPinned {
    T* p = nullptr;
    size_t n = 0;
    cudaError_t ensure(size_t cnt)
    {
        if (cnt <= n && p) return cudaSuccess;
        if (p) cudaFreeHost(p);
        p = nullptr;
        n = 0;
        const cudaError_t e = cudaMallocHost(&p, std::max<size_t>(cnt, 1) * sizeof(T));
        if (e == cudaSuccess) n = cnt;
        return e;
    }
    ~HostPinned()
    {
        if (p) cudaFreeHost(p);
    }
};

} // namespace

struct GfnffEwaldGPU::Impl {
    cudaStream_t stream = nullptr;
    std::string err;
    std::string fail_inject;   ///< CURCUMA_GFNFF_EWALD_DEVICE_FAIL
    bool profile = false;      ///< CURCUMA_GFNFF_EWALD_PROFILE=1: synchronise after each phase and print EWALDPROF lines (stderr)
    std::chrono::steady_clock::time_point tp;
    std::string prof;
    void mark(const char* what)
    {
        if (!profile) return;
        cudaStreamSynchronize(stream);
        const auto now = std::chrono::steady_clock::now();
        char buf[64];
        std::snprintf(buf, sizeof(buf), " %s=%.3f", what, std::chrono::duration<double, std::milli>(now - tp).count());
        prof += buf;
        tp = now;
    }
    void begin()
    {
        if (!profile) return;
        prof.clear();
        tp = std::chrono::steady_clock::now();
    }
    void flush(const char* tag)
    {
        if (profile) std::fprintf(stderr, "EWALDPROF %s%s\n", tag, prof.c_str());
    }
    double last_prepare_ms = 0.0;

    // key of the prepared data
    bool valid = false;
    int N = 0;
    double L[3] = { 0, 0, 0 };
    double alpha = 0.0, r_cut = 0.0, k_max = 0.0;
    int nk = 0, ngroups = 0, nslice = 1;
    int nmax[3] = { 0, 0, 0 };
    Eigen::MatrixXd pos_key;
    Eigen::VectorXd alp_key;
    CellArgs cell {};

    DevBuf<double> pos, alp, K, dK, re[3], im[3], kw, S, Sw, x, y, part, e_atom, v_atom, g, kitems, sums;
    DevBuf<float> pf;
    DevBuf<int> count, start, nbr, owner, knz;
    DevBuf<int4> groups;
    DevBuf<unsigned char> scan_tmp;
    HostPinned<double> h_vec, h_out;
    long long npairs_full = 0;

    bool ck(cudaError_t e, const char* what)
    {
        if (e == cudaSuccess) return true;
        err = std::string("GfnffEwaldGPU: ") + what + ": " + cudaGetErrorString(e);
        valid = false;
        return false;
    }
    PhaseArgs phases()
    {
        PhaseArgs t;
        for (int a = 0; a < 3; ++a) {
            t.re[a] = re[a].p;
            t.im[a] = im[a].p;
            t.nmax[a] = nmax[a];
        }
        return t;
    }
    KArgs kargs() const
    {
        KArgs k;
        k.groups = groups.p;
        k.knz = knz.p;
        k.kw = kw.p;
        k.ngroups = ngroups;
        k.nk = nk;
        return k;
    }
    /// structure factors of the device vector v (N)
    void structureFactors(const double* v)
    {
        k_ew_sf<<<ngroups, 256, 0, stream>>>(N, kargs(), phases(), v, S.p, Sw.p);
    }
};

GfnffEwaldGPU::GfnffEwaldGPU()
    : m(new Impl)
{
    if (const char* f = std::getenv("CURCUMA_GFNFF_EWALD_DEVICE_FAIL"))
        m->fail_inject = f;
    if (const char* f = std::getenv("CURCUMA_GFNFF_EWALD_PROFILE"))
        m->profile = f[0] && std::strcmp(f, "0") != 0;
    if (cudaStreamCreateWithFlags(&m->stream, cudaStreamNonBlocking) != cudaSuccess)
        m->stream = nullptr;
}

GfnffEwaldGPU::~GfnffEwaldGPU()
{
    if (m && m->stream)
        cudaStreamDestroy(m->stream);
}

std::string GfnffEwaldGPU::lastError() const { return m->err; }
double GfnffEwaldGPU::lastPrepareMs() const { return m->last_prepare_ms; }

// Claude Generated (Nov 2026, A2/S3): pair list, kernel values and phase table of this geometry on the device
bool GfnffEwaldGPU::prepare(const curcuma::ewald::EwaldSum& ew, const Eigen::MatrixXd& positions, const Eigen::VectorXd& alp)
{
    Impl& I = *m;
    I.last_prepare_ms = 0.0;
    if (!I.stream) {
        I.err = "GfnffEwaldGPU: no CUDA stream";
        return false;
    }
    if (I.fail_inject == "prepare") {
        I.err = "GfnffEwaldGPU: fault injection (CURCUMA_GFNFF_EWALD_DEVICE_FAIL=prepare)";
        I.valid = false;
        return false;
    }
    const int N = static_cast<int>(positions.rows());
    if (positions.cols() != 3 || alp.size() != N || N < 1) {
        I.err = "GfnffEwaldGPU: positions must be N x 3 and alp of length N";
        return false;
    }
    const curcuma::ewald::EwaldParameters& par = ew.parameters();
    const PBCUtils::OrthoCell& cell = ew.cell();
    const bool same_setup = I.valid && I.N == N && I.alpha == par.alpha && I.r_cut == par.r_cut && I.k_max == par.k_max
        && I.L[0] == cell.L[0] && I.L[1] == cell.L[1] && I.L[2] == cell.L[2] && I.nk == ew.numberOfKVectors();
    if (same_setup) {
        bool same = true;
        for (Eigen::Index c = 0; c < 3 && same; ++c)
            for (Eigen::Index r = 0; r < N; ++r)
                if (I.pos_key(r, c) != positions(r, c)) { same = false; break; }
        for (Eigen::Index r = 0; r < N && same; ++r)
            if (I.alp_key(r) != alp(r)) same = false;
        if (same) {
            ++m_stats.key_hits;
            return true;
        }
    }
    const auto t0 = std::chrono::steady_clock::now();
    I.valid = false;
    I.begin();
    if (!same_setup) {
        // k tables: groups of <= KGROUP consecutive k with the same (n_x, n_y) and consecutive n_z (list order of EwaldSum)
        const std::vector<double> kt = ew.kVectorTable();
        const std::vector<int> kn = ew.kVectorIndices();
        const int nk = ew.numberOfKVectors();
        std::vector<int4> grp;
        std::vector<int> nz(nk);
        for (int kk = 0; kk < nk; ++kk) {
            nz[kk] = kn[3 * kk + 2];
            const bool extend = !grp.empty() && grp.back().y < KGROUP && grp.back().z == kn[3 * kk] && grp.back().w == kn[3 * kk + 1]
                && kn[3 * (kk - 1) + 2] + 1 == kn[3 * kk + 2];
            if (extend)
                ++grp.back().y;
            else
                grp.push_back(make_int4(kk, 1, kn[3 * kk], kn[3 * kk + 1]));
        }
        I.nk = nk;
        I.ngroups = static_cast<int>(grp.size());
        for (int a = 0; a < 3; ++a)
            I.nmax[a] = ew.nmax(a);
        if (!I.ck(I.kw.ensure(kt.size()), "alloc k table") || !I.ck(I.knz.ensure(nk), "alloc k table")
            || !I.ck(I.groups.ensure(grp.size()), "alloc k groups") || !I.ck(I.S.ensure(2 * static_cast<size_t>(nk)), "alloc S") || !I.ck(I.Sw.ensure(2 * static_cast<size_t>(nk)), "alloc S")
            || !I.ck(I.kitems.ensure(7 * static_cast<size_t>(nk)), "alloc k items"))
            return false;
        if (!I.ck(cudaMemcpyAsync(I.kw.p, kt.data(), kt.size() * sizeof(double), cudaMemcpyHostToDevice, I.stream), "upload k table")
            || !I.ck(cudaMemcpyAsync(I.knz.p, nz.data(), nz.size() * sizeof(int), cudaMemcpyHostToDevice, I.stream), "upload k table")
            || !I.ck(cudaMemcpyAsync(I.groups.p, grp.data(), grp.size() * sizeof(int4), cudaMemcpyHostToDevice, I.stream), "upload groups")
            || !I.ck(cudaStreamSynchronize(I.stream), "k table sync"))
            return false;
        I.N = N;
        I.alpha = par.alpha;
        I.r_cut = par.r_cut;
        I.k_max = par.k_max;
        for (int a = 0; a < 3; ++a) {
            I.L[a] = cell.L[a];
            I.cell.L[a] = cell.L[a];
            I.cell.inv_L[a] = cell.inv_L[a];
            I.cell.Lf[a] = static_cast<float>(cell.L[a]);
            I.cell.inv_Lf[a] = static_cast<float>(cell.inv_L[a]);
        }
        I.cell.rc2 = par.r_cut * par.r_cut;
        I.cell.rcf2 = static_cast<float>((par.r_cut + 0.05) * (par.r_cut + 0.05));
        // slices of the per-atom reciprocal pass: enough blocks for the device (fixed per N and k list: deterministic)
        const int atom_blocks = (N + KPASS_ATOMS - 1) / KPASS_ATOMS;
        I.nslice = std::max(1, std::min(32, std::min(I.ngroups, (512 + atom_blocks - 1) / atom_blocks)));
        if (!I.ck(I.pos.ensure(3 * static_cast<size_t>(N)), "alloc positions") || !I.ck(I.pf.ensure(3 * static_cast<size_t>(N)), "alloc positions")
            || !I.ck(I.alp.ensure(N), "alloc alp") || !I.ck(I.count.ensure(N + 1), "alloc counts")
            || !I.ck(I.start.ensure(N + 1), "alloc starts") || !I.ck(I.x.ensure(N), "alloc x") || !I.ck(I.y.ensure(N), "alloc y")
            || !I.ck(I.part.ensure(3 * static_cast<size_t>(I.nslice) * N), "alloc partials")
            || !I.ck(I.e_atom.ensure(N), "alloc atom energies") || !I.ck(I.v_atom.ensure(6 * static_cast<size_t>(N)), "alloc atom virials")
            || !I.ck(I.g.ensure(3 * static_cast<size_t>(N)), "alloc gradient") || !I.ck(I.sums.ensure(16), "alloc sums")
            || !I.ck(I.h_vec.ensure(3 * static_cast<size_t>(N)), "alloc pinned") || !I.ck(I.h_out.ensure(16), "alloc pinned"))
            return false;
        for (int a = 0; a < 3; ++a)
            if (!I.ck(I.re[a].ensure(static_cast<size_t>(I.nmax[a] + 1) * N), "alloc phase table")
                || !I.ck(I.im[a].ensure(static_cast<size_t>(I.nmax[a] + 1) * N), "alloc phase table"))
                return false;
    }
    // positions atom-major (x y z per atom), widths
    for (int i = 0; i < N; ++i)
        for (int a = 0; a < 3; ++a)
            I.h_vec.p[3 * i + a] = positions(i, a);
    if (!I.ck(cudaMemcpyAsync(I.pos.p, I.h_vec.p, 3 * static_cast<size_t>(N) * sizeof(double), cudaMemcpyHostToDevice, I.stream), "upload positions"))
        return false;
    std::vector<double> alp_h(alp.data(), alp.data() + N);
    if (!I.ck(cudaMemcpyAsync(I.alp.p, alp_h.data(), N * sizeof(double), cudaMemcpyHostToDevice, I.stream), "upload alp"))
        return false;
    I.mark(same_setup ? "upload" : "setup+upload");
    k_ew_wrapf<<<(3 * N + 255) / 256, 256, 0, I.stream>>>(N, I.pos.p, I.cell, I.pf.p);
    k_ew_count<<<(N + 1 + 7) / 8, 256, 0, I.stream>>>(N, I.pf.p, I.cell, I.count.p);
    size_t tmp_bytes = 0;
    if (!I.ck(cub::DeviceScan::ExclusiveSum(nullptr, tmp_bytes, I.count.p, I.start.p, N + 1, I.stream), "scan size")
        || !I.ck(I.scan_tmp.ensure(tmp_bytes), "alloc scan")
        || !I.ck(cub::DeviceScan::ExclusiveSum(I.scan_tmp.p, tmp_bytes, I.count.p, I.start.p, N + 1, I.stream), "scan"))
        return false;
    int total = 0;
    if (!I.ck(cudaMemcpyAsync(&total, I.start.p + N, sizeof(int), cudaMemcpyDeviceToHost, I.stream), "download pair count")
        || !I.ck(cudaStreamSynchronize(I.stream), "pair count sync"))
        return false;
    I.mark("count+scan");
    if (total < 0) {
        I.err = "GfnffEwaldGPU: pair count overflow";
        return false;
    }
    if (!I.ck(I.nbr.ensure(total), "alloc pair list") || !I.ck(I.owner.ensure(total), "alloc pair list")
        || !I.ck(I.K.ensure(total), "alloc kernel") || !I.ck(I.dK.ensure(total), "alloc kernel"))
        return false;
    k_ew_fill<<<(N + 7) / 8, 256, 0, I.stream>>>(N, I.pf.p, I.cell, I.start.p, I.nbr.p, I.owner.p);
    if (total > 0)
        k_ew_kernel<<<static_cast<unsigned>((total + 255) / 256), 256, 0, I.stream>>>(total, I.pos.p, I.alp.p, I.cell, I.alpha,
            I.owner.p, I.nbr.p, I.K.p, I.dK.p);
    I.mark("fill");
    k_ew_phase<<<(3 * N + 127) / 128, 128, 0, I.stream>>>(N, I.pos.p, I.cell, I.phases());
    if (!I.ck(cudaGetLastError(), "launch") || !I.ck(cudaStreamSynchronize(I.stream), "build sync"))
        return false;
    I.mark("phase");
    if (I.profile) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), " N=%d nk=%d groups=%d slices=%d entries=%d", N, I.nk, I.ngroups, I.nslice, total);
        I.prof += buf;
    }
    I.flush("prepare");
    I.npairs_full = total;
    I.pos_key = positions;
    I.alp_key = alp;
    I.valid = true;
    ++m_stats.rebuilds;
    m_stats.full_pairs = total;
    I.last_prepare_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

// Claude Generated (Nov 2026, A2/S3b): y = (A_real + A_recip) x, the operator of the periodic pcg
bool GfnffEwaldGPU::potentials(const Eigen::VectorXd& xin, Eigen::VectorXd& yout)
{
    Impl& I = *m;
    if (!I.valid) {
        I.err = "GfnffEwaldGPU: potentials() before a successful prepare()";
        return false;
    }
    if (I.fail_inject == "potentials") {
        I.err = "GfnffEwaldGPU: fault injection (CURCUMA_GFNFF_EWALD_DEVICE_FAIL=potentials)";
        return false;
    }
    if (xin.size() != I.N) {
        I.err = "GfnffEwaldGPU: vector length differs from the prepared geometry";
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const int N = I.N;
    std::memcpy(I.h_vec.p, xin.data(), N * sizeof(double));
    if (!I.ck(cudaMemcpyAsync(I.x.p, I.h_vec.p, N * sizeof(double), cudaMemcpyHostToDevice, I.stream), "upload x"))
        return false;
    I.begin();
    I.structureFactors(I.x.p);
    I.mark("sf");
    k_ew_kpass<false><<<dim3((N + KPASS_ATOMS - 1) / KPASS_ATOMS, I.nslice), KPASS_ATOMS, 0, I.stream>>>(
        N, I.kargs(), I.phases(), I.Sw.p, nullptr, I.nslice, I.part.p);
    I.mark("kpot");
    k_ew_real_pot<<<(N + 7) / 8, 256, 0, I.stream>>>(N, I.start.p, I.nbr.p, I.K.p, I.x.p, I.nslice, I.part.p, I.y.p);
    I.mark("real");
    if (!I.ck(cudaGetLastError(), "launch potentials")
        || !I.ck(cudaMemcpyAsync(I.h_vec.p, I.y.p, N * sizeof(double), cudaMemcpyDeviceToHost, I.stream), "download y")
        || !I.ck(cudaStreamSynchronize(I.stream), "potentials sync"))
        return false;
    I.mark("down");
    I.flush("potentials");
    yout.resize(N);
    std::memcpy(yout.data(), I.h_vec.p, N * sizeof(double));
    ++m_stats.potentials_calls;
    m_stats.potentials_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

// Claude Generated (Nov 2026, A2/S3a): pair energy, gradient and virial at fixed charges
bool GfnffEwaldGPU::evaluate(const Eigen::VectorXd& qin, curcuma::ewald::EwaldEnergyTerms& terms, Eigen::MatrixXd* grad,
    Eigen::Matrix3d* virial)
{
    Impl& I = *m;
    if (!I.valid) {
        I.err = "GfnffEwaldGPU: evaluate() before a successful prepare()";
        return false;
    }
    if (I.fail_inject == "evaluate") {
        I.err = "GfnffEwaldGPU: fault injection (CURCUMA_GFNFF_EWALD_DEVICE_FAIL=evaluate)";
        return false;
    }
    if (qin.size() != I.N) {
        I.err = "GfnffEwaldGPU: charge vector length differs from the prepared geometry";
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const int N = I.N;
    const bool want_g = grad != nullptr;
    const bool want_v = virial != nullptr && want_g;
    std::memcpy(I.h_vec.p, qin.data(), N * sizeof(double));
    if (!I.ck(cudaMemcpyAsync(I.x.p, I.h_vec.p, N * sizeof(double), cudaMemcpyHostToDevice, I.stream), "upload q"))
        return false;
    I.structureFactors(I.x.p);
    if (want_g)
        k_ew_kpass<true><<<dim3((N + KPASS_ATOMS - 1) / KPASS_ATOMS, I.nslice), KPASS_ATOMS, 0, I.stream>>>(
            N, I.kargs(), I.phases(), I.Sw.p, I.x.p, I.nslice, I.part.p);
    const dim3 rg((N + 7) / 8), rb(256);
    double* gx = I.g.p;
    double* gy = I.g.p + N;
    double* gz = I.g.p + 2 * static_cast<size_t>(N);
    if (want_v)
        k_ew_real_eg<true, true><<<rg, rb, 0, I.stream>>>(N, I.pos.p, I.cell, I.start.p, I.nbr.p, I.K.p, I.dK.p, I.x.p,
            I.nslice, I.part.p, I.e_atom.p, gx, gy, gz, I.v_atom.p);
    else if (want_g)
        k_ew_real_eg<true, false><<<rg, rb, 0, I.stream>>>(N, I.pos.p, I.cell, I.start.p, I.nbr.p, I.K.p, I.dK.p, I.x.p,
            I.nslice, I.part.p, I.e_atom.p, gx, gy, gz, I.v_atom.p);
    else
        k_ew_real_eg<false, false><<<rg, rb, 0, I.stream>>>(N, I.pos.p, I.cell, I.start.p, I.nbr.p, I.K.p, I.dK.p, I.x.p,
            I.nslice, I.part.p, I.e_atom.p, gx, gy, gz, I.v_atom.p);
    k_ew_kitems<<<(I.nk + 255) / 256, 256, 0, I.stream>>>(I.nk, I.kw.p, I.S.p, I.alpha, I.kitems.p);
    k_ew_sum_items<7><<<1, 1024, 0, I.stream>>>(I.nk, I.kitems.p, I.sums.p);
    k_ew_sum_items<1><<<1, 1024, 0, I.stream>>>(N, I.e_atom.p, I.sums.p + 7);
    if (want_v)
        k_ew_sum_items<6><<<1, 1024, 0, I.stream>>>(N, I.v_atom.p, I.sums.p + 8);
    if (!I.ck(cudaGetLastError(), "launch evaluate")
        || !I.ck(cudaMemcpyAsync(I.h_out.p, I.sums.p, 14 * sizeof(double), cudaMemcpyDeviceToHost, I.stream), "download sums"))
        return false;
    if (want_g && !I.ck(cudaMemcpyAsync(I.h_vec.p, I.g.p, 3 * static_cast<size_t>(N) * sizeof(double), cudaMemcpyDeviceToHost, I.stream),
                     "download gradient"))
        return false;
    if (!I.ck(cudaStreamSynchronize(I.stream), "evaluate sync"))
        return false;
    const double* o = I.h_out.p;
    terms.reciprocal = o[0];
    terms.real = o[7];
    if (want_g) {
        grad->resize(N, 3);
        std::memcpy(grad->data(), I.h_vec.p, 3 * static_cast<size_t>(N) * sizeof(double));  // column-major = gx, gy, gz
    }
    if (want_v) {
        Eigen::Matrix3d v;
        v(0, 0) = o[1] + o[8];
        v(1, 1) = o[2] + o[9];
        v(2, 2) = o[3] + o[10];
        v(0, 1) = v(1, 0) = o[4] + o[11];
        v(0, 2) = v(2, 0) = o[5] + o[12];
        v(1, 2) = v(2, 1) = o[6] + o[13];
        *virial = v;
    } else if (virial) {
        virial->setZero();
    }
    ++m_stats.evaluate_calls;
    m_stats.evaluate_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}
