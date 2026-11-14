/*
 * < CUDA evaluation and fitting of the Delta-ML model >
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
 * Claude Generated (July 2026, rewritten October 2026): GPU path of the Delta-ML work package.
 *
 * What runs here (October 2026): the whole correction. Coordinates and per-atom element
 * indices go in; energy, extrapolation indicator and gradient come out. In between the device
 * builds the neighbour lists, the ACSF descriptor (Behler G2/G4, identical definitions to
 * descriptor.cpp), the Gaussian kernel against the resident centres, the feature-space adjoint
 * W and its contraction with the descriptor derivatives into forces. Neither the descriptor
 * matrix nor the compact Jacobian ever crosses the bus.
 *
 * Why: measured on 58 and 90 atoms (docs/DELTA_ML_GPU_2026-10.md), the July version moved only
 * the kernel and spent 80-90 % of the step on the host building descriptor and Jacobian and
 * contracting them; the fit spent two thirds of its time building the design block on the host.
 *
 * Determinism: no floating-point atomics. Every sum runs in a fixed order, and the per-atom
 * gradient is gathered over centres in ascending index, the same order descriptor.cpp /
 * krr_model.cpp accumulate in. Results agree with the CPU to rounding (exp/pow/sqrt and FMA
 * contraction differ), not bit for bit.
 *
 * Fit: per training structure the design block (one energy row plus 3N gradient rows against
 * the M centres) is assembled on the device from the device-side compact Jacobian and folded
 * into the normal equations with cuBLAS DSYRK/DGEMV, as before. Only the labels go up.
 *
 * Storage convention of the design block: row-major (rows x M). Read with leading dimension M
 * it is the column-major Phi^T that DSYRK/DGEMV want, so nothing is transposed.
 */

#include "../ml_gpu_abi.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

constexpr int kThreads = 128; // multiple of 32 (warp ballot in the neighbour kernel)
/// Per neighbour: dx, dy, dz, r, fc_radial(r), dfc_radial(r), fc_angular(r), dfc_angular(r).
/// The cutoff values are computed once here instead of once per pair and triplet: on a card
/// with 1/64-rate FP64 the transcendental functions are what the step time is made of.
constexpr int kNbv = 8;

// =================================================================================
// Device helpers -- same formulas as ACSFDescriptor::cutoff / cutoffDerivative
// =================================================================================

__device__ inline double fcut(double r, double rc)
{
    if (r >= rc)
        return 0.0;
    return 0.5 * (cos(M_PI * r / rc) + 1.0);
}

__device__ inline double dfcut(double r, double rc)
{
    if (r >= rc)
        return 0.0;
    return -0.5 * M_PI / rc * sin(M_PI * r / rc);
}

/**
 * @brief base^zeta and base^(zeta-1). Integer zeta (the default grid 1, 4, 16) by repeated
 * squaring, otherwise pow(); agrees with std::pow of descriptor.cpp to a few ulp.
 */
__device__ inline void powPair(double base, double zeta, int zi, double& P, double& Pm1)
{
    if (zi > 0) {
        double r = 1.0, b = base;
        for (int e = zi - 1; e; e >>= 1) {
            if (e & 1)
                r *= b;
            b *= b;
        }
        Pm1 = r;
        P = r * base;
    } else {
        P = pow(base, zeta);
        Pm1 = pow(base, zeta - 1.0);
    }
}

/*! \brief Position of the pair (ta < va) in the row-wise enumeration of the na x na upper triangle. */
__device__ inline int tripletIndex(int ta, int va, int na)
{
    return ta * na - (ta * (ta + 1)) / 2 + (va - ta - 1);
}

/*! \brief Inverse of tripletIndex. */
__device__ inline void tripletPair(int k, int na, int& ta, int& va)
{
    ta = 0;
    int row = na - 1;
    while (k >= row) {
        k -= row;
        ++ta;
        --row;
    }
    va = ta + 1 + k;
}

/*! \brief Shared-memory bytes per staged triplet: descriptor (NA values + pair id), contraction (two 3-vectors). */
inline int descriptorTripletBytes(const CurcumaMLDeviceConfig& c)
{
    return c.n_zeta * c.n_lambda * static_cast<int>(sizeof(double)) + static_cast<int>(sizeof(int));
}
constexpr int kContractTripletBytes = 6 * sizeof(double);
/// Dynamic shared memory for the triplet stage; beyond it the kernels fall back to the
/// unstaged loops, which are slower but give the same sums in the same order.
constexpr int kStageBytes = 40 * 1024;

/*! \brief ACSFConfig::pairIndex: row-wise upper-triangle index including the diagonal. */
__device__ inline int pairIndex(int a, int b, int s)
{
    const int lo = a < b ? a : b;
    const int hi = a < b ? b : a;
    return lo * s - (lo * (lo - 1)) / 2 + (hi - lo);
}

/*! \brief Slot of atom @p a in the ascending neighbour list of length n, or -1. */
__device__ inline int findSlot(const int* list, int n, int a)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        const int mid = (lo + hi) >> 1;
        const int v = list[mid];
        if (v == a)
            return mid;
        if (v < a)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

/**
 * @brief Collect, in ascending slot order, the neighbours inside the angular cutoff.
 * Run by one thread; the list is short (atoms within rcut_angular) and the order matters.
 */
__device__ inline int angularList(int nnb, const double* nbv, double rc_a, int* ang)
{
    int na = 0;
    for (int t = 0; t < nnb; ++t)
        if (nbv[kNbv * t + 3] < rc_a)
            ang[na++] = t;
    return na;
}

// =================================================================================
// Kernels
// =================================================================================

/**
 * @brief Ordered neighbour list per centre atom: j ascending, j != i, both elements known,
 * 1e-8 <= R_ij < max(rcut_radial, rcut_angular) -- the filter of ACSFDescriptor::compute.
 * nbv holds (dx, dy, dz, r) with d = R_j - R_i.
 */
__global__ void kNeighbours(int n, int cap, double rc_r, double rc_a, const double* __restrict__ R,
    const int* __restrict__ sidx, int* __restrict__ nnb, int* __restrict__ nb,
    double* __restrict__ nbv, int* __restrict__ flag)
{
    __shared__ int warp_cnt[kThreads / 32];
    __shared__ int base;
    const int i = blockIdx.x;
    if (threadIdx.x == 0)
        base = 0;
    __syncthreads();
    if (sidx[i] < 0) {
        if (threadIdx.x == 0)
            nnb[i] = 0;
        return;
    }
    const double xi = R[3 * i], yi = R[3 * i + 1], zi = R[3 * i + 2];
    const int lane = threadIdx.x & 31;
    const int w = threadIdx.x >> 5;
    const double rc_max = rc_r > rc_a ? rc_r : rc_a;

    for (int j0 = 0; j0 < n; j0 += blockDim.x) {
        const int j = j0 + threadIdx.x;
        bool ok = false;
        double vx = 0, vy = 0, vz = 0, r = 0;
        if (j < n && j != i && sidx[j] >= 0) {
            vx = R[3 * j] - xi;
            vy = R[3 * j + 1] - yi;
            vz = R[3 * j + 2] - zi;
            r = sqrt(vx * vx + vy * vy + vz * vz);
            ok = (r < rc_max && r >= 1e-8);
        }
        const unsigned ball = __ballot_sync(0xffffffffu, ok);
        if (lane == 0)
            warp_cnt[w] = __popc(ball);
        __syncthreads();
        int off = base;
        for (int q = 0; q < w; ++q)
            off += warp_cnt[q];
        off += __popc(ball & ((1u << lane) - 1u));
        if (ok && off < cap) {
            const size_t s = static_cast<size_t>(i) * cap + off;
            nb[s] = j;
            double* o = nbv + kNbv * s;
            o[0] = vx;
            o[1] = vy;
            o[2] = vz;
            o[3] = r;
            o[4] = fcut(r, rc_r);
            o[5] = dfcut(r, rc_r);
            o[6] = fcut(r, rc_a);
            o[7] = dfcut(r, rc_a);
        }
        __syncthreads();
        if (threadIdx.x == 0) {
            int tot = 0;
            for (int q = 0; q < static_cast<int>(blockDim.x) / 32; ++q)
                tot += warp_cnt[q];
            base += tot;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        // Clamped so the following kernels never read past this row; the flag makes the
        // whole call fail before any result leaves the device.
        nnb[i] = base < cap ? base : cap;
        if (base > cap)
            atomicExch(flag, 1);
    }
}

/**
 * @brief Standardised ACSF descriptor, one block per centre atom.
 * Angular triplets are first evaluated once each into shared memory (one thread per triplet),
 * then one thread per feature adds them up; radial features loop over the neighbours. Sums run
 * in ascending neighbour / pair (t < v) order as in ACSFDescriptor::compute.
 * Xs = (x - mean) * inv_scale (DeltaMLModel::standardise).
 */
__global__ void kDescriptor(int F, int cap, int tripcap, CurcumaMLDeviceConfig c, const int* __restrict__ sidx,
    const int* __restrict__ nnb_all, const int* __restrict__ nb_all,
    const double* __restrict__ nbv_all, const double* __restrict__ mean,
    const double* __restrict__ inv_scale, double* __restrict__ Xs)
{
    __shared__ int ang[CURCUMA_ML_MAX_NEIGHBOURS];
    __shared__ int na;
    extern __shared__ double stage[]; // tripcap x NA values, then tripcap pair ids
    const int i = blockIdx.x;
    const bool known = sidx[i] >= 0;
    const int nnb = known ? nnb_all[i] : 0;
    const int* nb = nb_all + static_cast<size_t>(i) * cap;
    const double* nbv = nbv_all + static_cast<size_t>(i) * cap * kNbv;
    if (threadIdx.x == 0)
        na = angularList(nnb, nbv, c.rcut_angular, ang);
    __syncthreads();

    const int S = c.n_species;
    const int NR = c.n_radial;
    const int NA = c.n_zeta * c.n_lambda;

    // Triplet stage: every angular triplet (t < v) is evaluated once, by one thread, and its
    // NA terms are kept in shared memory; the feature threads below then only add them up in
    // ascending triplet order (= the loop order of descriptor.cpp).
    const int ntrip = na * (na - 1) / 2;
    const bool staged = known && ntrip <= tripcap;
    double* tval = stage;
    int* tpair = reinterpret_cast<int*>(stage + static_cast<size_t>(tripcap) * NA);
    if (staged) {
        const double rc = c.rcut_angular;
        for (int k = threadIdx.x; k < ntrip; k += blockDim.x) {
            int ta, va;
            tripletPair(k, na, ta, va);
            const double* vt = nbv + kNbv * ang[ta];
            const double* vv = nbv + kNbv * ang[va];
            const double a = vt[3], b = vv[3];
            const double jx = vv[0] - vt[0], jy = vv[1] - vt[1], jz = vv[2] - vt[2];
            const double cc = sqrt(jx * jx + jy * jy + jz * jz);
            if (cc >= rc || cc < 1e-8) {
                tpair[k] = -1;
                continue;
            }
            tpair[k] = pairIndex(sidx[nb[ang[ta]]], sidx[nb[ang[va]]], S);
            const double cos_theta = (vt[0] * vv[0] + vt[1] * vv[1] + vt[2] * vv[2]) / (a * b);
            const double Fc = vt[6] * vv[6] * fcut(cc, rc);
            const double E = exp(-c.angular_eta * (a * a + b * b + cc * cc));
            int q = 0;
            for (int z = 0; z < c.n_zeta; ++z) {
                for (int l = 0; l < c.n_lambda; ++l, ++q) {
                    double base = 1.0 + c.lambda[l] * cos_theta;
                    if (base < 0.0)
                        base = 0.0;
                    double P, Pm1;
                    powPair(base, c.zeta[z], c.zeta_int[z], P, Pm1);
                    tval[static_cast<size_t>(k) * NA + q] = c.zeta_pref[z] * P * E * Fc;
                }
            }
        }
    }
    __syncthreads();

    for (int f = threadIdx.x; f < F; f += blockDim.x) {
        double x = 0.0;
        if (known && f < S * NR) {
            const int s = f / NR;
            const int p = f % NR;
            for (int t = 0; t < nnb; ++t) {
                const double r = nbv[kNbv * t + 3];
                if (r >= c.rcut_radial || sidx[nb[t]] != s)
                    continue;
                const double dr = r - c.radial_rs[p];
                x += exp(-c.radial_eta * dr * dr) * nbv[kNbv * t + 4];
            }
        } else if (staged) {
            const int fa = f - S * NR;
            const int pi = fa / NA;
            const int q = fa % NA;
            for (int k = 0; k < ntrip; ++k)
                if (tpair[k] == pi)
                    x += tval[static_cast<size_t>(k) * NA + q];
        } else if (known) {
            const int fa = f - S * NR;
            const int pi = fa / NA;
            const int q = fa % NA;
            const int z = q / c.n_lambda;
            const int l = q % c.n_lambda;
            const double zeta = c.zeta[z];
            const int zi = c.zeta_int[z];
            const double pref = c.zeta_pref[z];
            const double lambda = c.lambda[l];
            const double rc = c.rcut_angular;
            for (int ta = 0; ta < na; ++ta) {
                const int t = ang[ta];
                const int st = sidx[nb[t]];
                const double* vt = nbv + kNbv * t;
                const double a = vt[3];
                for (int va = ta + 1; va < na; ++va) {
                    const int v = ang[va];
                    if (pairIndex(st, sidx[nb[v]], S) != pi)
                        continue;
                    const double* vv = nbv + kNbv * v;
                    const double b = vv[3];
                    const double jx = vv[0] - vt[0], jy = vv[1] - vt[1], jz = vv[2] - vt[2];
                    const double cc = sqrt(jx * jx + jy * jy + jz * jz);
                    if (cc >= rc || cc < 1e-8)
                        continue;
                    const double cos_theta = (vt[0] * vv[0] + vt[1] * vv[1] + vt[2] * vv[2]) / (a * b);
                    const double Fc = vt[6] * vv[6] * fcut(cc, rc);
                    const double E = exp(-c.angular_eta * (a * a + b * b + cc * cc));
                    double base = 1.0 + lambda * cos_theta;
                    if (base < 0.0)
                        base = 0.0;
                    double P, Pm1;
                    powPair(base, zeta, zi, P, Pm1);
                    x += pref * P * E * Fc;
                }
            }
        }
        Xs[static_cast<size_t>(i) * F + f] = (x - mean[f]) * inv_scale[f];
    }
}

/*! \brief k[i,m] = exp(-|x_i - z_m|^2 / 2 sigma^2) inside atom i's element block, else 0. */
__global__ void kKernel(int natoms, int M, int F, double inv_two_sigma2, const double* __restrict__ X,
    const double* __restrict__ Z, const int* __restrict__ start, const int* __restrict__ count,
    double* __restrict__ k)
{
    const int m = blockIdx.x * blockDim.x + threadIdx.x;
    const int i = blockIdx.y;
    if (i >= natoms || m >= M)
        return;
    const int s0 = start[i];
    const int c0 = count[i];
    if (c0 <= 0 || m < s0 || m >= s0 + c0) {
        k[static_cast<size_t>(i) * M + m] = 0.0;
        return;
    }
    double d2 = 0.0;
    for (int f = 0; f < F; ++f) {
        const double diff = X[static_cast<size_t>(i) * F + f] - Z[static_cast<size_t>(m) * F + f];
        d2 += diff * diff;
    }
    k[static_cast<size_t>(i) * M + m] = exp(-d2 * inv_two_sigma2);
}

/*! \brief Per atom: s_i = sum_m beta_m k_im and the largest k (the extrapolation indicator). */
__global__ void kReduce(int M, const double* __restrict__ k, const double* __restrict__ beta,
    double* __restrict__ s, double* __restrict__ kmax)
{
    __shared__ double sh_s[kThreads];
    __shared__ double sh_m[kThreads];
    const int i = blockIdx.x;
    double acc = 0.0;
    double best = 0.0;
    for (int m = threadIdx.x; m < M; m += blockDim.x) {
        const double kv = k[static_cast<size_t>(i) * M + m];
        acc += beta[m] * kv;
        best = fmax(best, kv);
    }
    sh_s[threadIdx.x] = acc;
    sh_m[threadIdx.x] = best;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            sh_s[threadIdx.x] += sh_s[threadIdx.x + stride];
            sh_m[threadIdx.x] = fmax(sh_m[threadIdx.x], sh_m[threadIdx.x + stride]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        s[i] = sh_s[0];
        kmax[i] = sh_m[0];
    }
}

/**
 * @brief Adjoint in raw descriptor space: W_if = (s_i x_if - sum_m beta_m k_im z_mf) * (-1/sigma^2) / scale_f.
 * Block (i, feature tile): threadIdx.x picks the feature, threadIdx.y one of kAdjSplit
 * interleaved parts of the centre block; the parts are added in fixed order.
 */
constexpr int kAdjSplit = 8;
__global__ void kAdjoint(int natoms, int M, int F, double neg_inv_sigma2, const double* __restrict__ X,
    const double* __restrict__ Z, const double* __restrict__ k, const double* __restrict__ beta,
    const double* __restrict__ s, const double* __restrict__ inv_scale, const int* __restrict__ start,
    const int* __restrict__ count, double* __restrict__ W)
{
    __shared__ double part[kAdjSplit][32];
    const int f = blockIdx.x * 32 + threadIdx.x;
    const int i = blockIdx.y;
    const int s0 = start[i];
    const int c0 = count[i];
    double dot = 0.0;
    if (f < F)
        for (int m = s0 + threadIdx.y; m < s0 + c0; m += kAdjSplit)
            dot += beta[m] * k[static_cast<size_t>(i) * M + m] * Z[static_cast<size_t>(m) * F + f];
    part[threadIdx.y][threadIdx.x] = dot;
    __syncthreads();
    if (threadIdx.y == 0 && f < F) {
        double sum = 0.0;
        for (int y = 0; y < kAdjSplit; ++y)
            sum += part[y][threadIdx.x];
        const double v = s[i] * X[static_cast<size_t>(i) * F + f] - sum;
        W[static_cast<size_t>(i) * F + f] = v * neg_inv_sigma2 * inv_scale[f];
    }
}

/**
 * @brief Derivative of one centre's descriptor with respect to neighbour slot t.
 *
 * Shared by the force contraction (Jrow == nullptr: returns sum_f w_f dX_f/dR_t in g) and
 * the Jacobian build (Jrow != nullptr: writes dX_f/dR_t,gamma into Jrow[gamma*F + f]).
 * Formulas are those of ACSFDescriptor::compute: slot t takes the role "j" of the triplets
 * (t, v > t) and the role "k" of the triplets (v < t, t).
 */
__device__ void slotDerivative(int t, int F, const CurcumaMLDeviceConfig& c, const int* sidx,
    const int* nb, const double* nbv, const int* ang, int na, const double* w, double* Jrow,
    double g[3], bool with_angular = true)
{
    g[0] = g[1] = g[2] = 0.0;
    const int S = c.n_species;
    const int NR = c.n_radial;
    const int NA = c.n_zeta * c.n_lambda;
    const double* vt = nbv + kNbv * t;
    const double r = vt[3];

    // --- radial (G2) ---
    if (r < c.rcut_radial) {
        const int s = sidx[nb[t]];
        const double fc = vt[4];
        const double dfc = vt[5];
        const double u[3] = { vt[0] / r, vt[1] / r, vt[2] / r };
        double D = 0.0;
        for (int p = 0; p < NR; ++p) {
            const double dr = r - c.radial_rs[p];
            const double e = exp(-c.radial_eta * dr * dr);
            const double dval = e * (-2.0 * c.radial_eta * dr) * fc + e * dfc;
            const int f = s * NR + p;
            if (Jrow) {
                for (int gg = 0; gg < 3; ++gg)
                    Jrow[gg * F + f] += dval * u[gg];
            } else {
                D += dval * w[f];
            }
        }
        for (int gg = 0; gg < 3; ++gg)
            g[gg] += D * u[gg];
    }

    // --- angular (G4) ---
    const double rc = c.rcut_angular;
    if (!with_angular || r >= rc)
        return;
    for (int va = 0; va < na; ++va) {
        const int v = ang[va];
        if (v == t)
            continue;
        const bool t_is_j = t < v;
        const int jj = t_is_j ? t : v;
        const int kk = t_is_j ? v : t;
        const double* pj = nbv + kNbv * jj;
        const double* pk = nbv + kNbv * kk;
        const double a = pj[3];
        const double b = pk[3];
        const double rjk[3] = { pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2] };
        const double cc = sqrt(rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]);
        if (cc >= rc || cc < 1e-8)
            continue;

        const double cos_theta = (pj[0] * pk[0] + pj[1] * pk[1] + pj[2] * pk[2]) / (a * b);
        const int pi = pairIndex(sidx[nb[jj]], sidx[nb[kk]], S);
        const double fa = pj[6], fb = pk[6];
        const double dfa = pj[7], dfb = pk[7];
        double sn, cs;
        sincos(M_PI * cc / rc, &sn, &cs); // cc < rc here, so this is fcut/dfcut of descriptor.cpp
        const double fcc = 0.5 * (cs + 1.0);
        const double dfcc = -0.5 * M_PI / rc * sn;
        const double Fc = fa * fb * fcc;
        const double E = exp(-c.angular_eta * (a * a + b * b + cc * cc));

        // Geometric vectors of this slot's role.
        double v_cos[3], v_r[3], v_c[3];
        for (int gg = 0; gg < 3; ++gg) {
            const double u_c = rjk[gg] / cc;
            if (t_is_j) {
                v_cos[gg] = pk[gg] / (a * b) - cos_theta * pj[gg] / (a * a); // dcos/dr_j
                v_r[gg] = pj[gg] / a;                                         // da/dr_j
                v_c[gg] = -u_c;                                               // dc/dr_j
            } else {
                v_cos[gg] = pj[gg] / (a * b) - cos_theta * pk[gg] / (b * b); // dcos/dr_k
                v_r[gg] = pk[gg] / b;                                         // db/dr_k
                v_c[gg] = u_c;                                                // dc/dr_k
            }
        }

        double A_cos = 0.0, A_r = 0.0, A_c = 0.0;
        int q = 0;
        for (int z = 0; z < c.n_zeta; ++z) {
            const double zeta = c.zeta[z];
            const int zi = c.zeta_int[z];
            const double pref = c.zeta_pref[z];
            for (int l = 0; l < c.n_lambda; ++l, ++q) {
                const double lambda = c.lambda[l];
                double base = 1.0 + lambda * cos_theta;
                if (base < 0.0)
                    base = 0.0;
                double P, Pm1;
                powPair(base, zeta, zi, P, Pm1);
                const double dP = zeta * lambda * Pm1;
                const double dt_dcos = pref * dP * E * Fc;
                const double dt_dr = t_is_j
                    ? pref * P * ((-2.0 * c.angular_eta * a) * E * Fc + E * dfa * fb * fcc)
                    : pref * P * ((-2.0 * c.angular_eta * b) * E * Fc + E * fa * dfb * fcc);
                const double dt_dc = pref * P * ((-2.0 * c.angular_eta * cc) * E * Fc + E * fa * fb * dfcc);
                const int f = S * NR + pi * NA + q;
                if (Jrow) {
                    for (int gg = 0; gg < 3; ++gg)
                        Jrow[gg * F + f] += dt_dcos * v_cos[gg] + dt_dr * v_r[gg] + dt_dc * v_c[gg];
                } else {
                    A_cos += w[f] * dt_dcos;
                    A_r += w[f] * dt_dr;
                    A_c += w[f] * dt_dc;
                }
            }
        }
        for (int gg = 0; gg < 3; ++gg)
            g[gg] += A_cos * v_cos[gg] + A_r * v_r[gg] + A_c * v_c[gg];
    }
}

/**
 * @brief Contracted derivative of one angular triplet (j = slot jj, k = slot kk, jj < kk) with
 * the adjoint w: gj = sum_f w_f dG4_f/dR_j, gk likewise (descriptor.cpp formulas).
 * Returns false for a j-k pair outside the cutoff (no contribution).
 */
__device__ bool tripletGradient(int jj, int kk, int F, const CurcumaMLDeviceConfig& c, const int* sidx,
    const int* nb, const double* nbv, const double* w, double gj[3], double gk[3])
{
    const double rc = c.rcut_angular;
    const double* pj = nbv + kNbv * jj;
    const double* pk = nbv + kNbv * kk;
    const double a = pj[3];
    const double b = pk[3];
    const double rjk[3] = { pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2] };
    const double cc = sqrt(rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]);
    if (cc >= rc || cc < 1e-8)
        return false;
    const int S = c.n_species;
    const int NA = c.n_zeta * c.n_lambda;
    const double cos_theta = (pj[0] * pk[0] + pj[1] * pk[1] + pj[2] * pk[2]) / (a * b);
    const int pi = pairIndex(sidx[nb[jj]], sidx[nb[kk]], S);
    const double fa = pj[6], fb = pk[6];
    const double dfa = pj[7], dfb = pk[7];
    double sn, cs;
    sincos(M_PI * cc / rc, &sn, &cs);
    const double fcc = 0.5 * (cs + 1.0);
    const double dfcc = -0.5 * M_PI / rc * sn;
    const double Fc = fa * fb * fcc;
    const double E = exp(-c.angular_eta * (a * a + b * b + cc * cc));

    double A_cos = 0.0, A_a = 0.0, A_b = 0.0, A_c = 0.0;
    const double* wf = w + S * c.n_radial + pi * NA;
    int q = 0;
    for (int z = 0; z < c.n_zeta; ++z) {
        const double zeta = c.zeta[z];
        const double pref = c.zeta_pref[z];
        for (int l = 0; l < c.n_lambda; ++l, ++q) {
            const double lambda = c.lambda[l];
            double base = 1.0 + lambda * cos_theta;
            if (base < 0.0)
                base = 0.0;
            double P, Pm1;
            powPair(base, zeta, c.zeta_int[z], P, Pm1);
            const double ww = wf[q];
            A_cos += ww * (pref * zeta * lambda * Pm1 * E * Fc);
            A_a += ww * (pref * P * ((-2.0 * c.angular_eta * a) * E * Fc + E * dfa * fb * fcc));
            A_b += ww * (pref * P * ((-2.0 * c.angular_eta * b) * E * Fc + E * fa * dfb * fcc));
            A_c += ww * (pref * P * ((-2.0 * c.angular_eta * cc) * E * Fc + E * fa * fb * dfcc));
        }
    }
    for (int g = 0; g < 3; ++g) {
        const double u_c = rjk[g] / cc;
        const double dcos_dj = pk[g] / (a * b) - cos_theta * pj[g] / (a * a);
        const double dcos_dk = pj[g] / (a * b) - cos_theta * pk[g] / (b * b);
        gj[g] = A_cos * dcos_dj + A_a * (pj[g] / a) - A_c * u_c;
        gk[g] = A_cos * dcos_dk + A_b * (pk[g] / b) + A_c * u_c;
    }
    return true;
}

/**
 * @brief Force contraction per centre: C[i][t+1] = sum_f W_if dX_if/dR_{nb_t}, and the centre
 * slot C[i][0] = -sum_t C[i][t+1] (translation invariance, as in descriptor.cpp).
 */
__global__ void kContract(int F, int cap, int tripcap, CurcumaMLDeviceConfig c, const int* __restrict__ sidx,
    const int* __restrict__ nnb_all, const int* __restrict__ nb_all, const double* __restrict__ nbv_all,
    const double* __restrict__ W, double* __restrict__ C)
{
    __shared__ int ang[CURCUMA_ML_MAX_NEIGHBOURS];
    __shared__ int angpos[CURCUMA_ML_MAX_NEIGHBOURS];
    __shared__ int na;
    extern __shared__ double stage[]; // tripcap x (gj, gk)
    const int i = blockIdx.x;
    const int nnb = sidx[i] >= 0 ? nnb_all[i] : 0;
    const int* nb = nb_all + static_cast<size_t>(i) * cap;
    const double* nbv = nbv_all + static_cast<size_t>(i) * cap * kNbv;
    double* Ci = C + static_cast<size_t>(i) * (cap + 1) * 3;
    if (threadIdx.x == 0) {
        na = angularList(nnb, nbv, c.rcut_angular, ang);
        for (int t = 0; t < nnb; ++t)
            angpos[t] = -1;
        for (int ta = 0; ta < na; ++ta)
            angpos[ang[ta]] = ta;
    }
    __syncthreads();

    const double* w = W + static_cast<size_t>(i) * F;
    const int ntrip = na * (na - 1) / 2;
    const bool staged = ntrip <= tripcap;
    if (staged) {
        // Each triplet once: its two contracted 3-vectors go to shared memory.
        for (int k = threadIdx.x; k < ntrip; k += blockDim.x) {
            int ta, va;
            tripletPair(k, na, ta, va);
            double* o = stage + 6 * static_cast<size_t>(k);
            if (!tripletGradient(ang[ta], ang[va], F, c, sidx, nb, nbv, w, o, o + 3))
                for (int e = 0; e < 6; ++e)
                    o[e] = 0.0;
        }
    }
    __syncthreads();

    for (int t = threadIdx.x; t < nnb; t += blockDim.x) {
        double g[3] = { 0.0, 0.0, 0.0 };
        if (staged) {
            // radial part only, then the staged triplets of this slot in ascending partner order
            slotDerivative(t, F, c, sidx, nb, nbv, ang, na, w, nullptr, g, false);
            const int ta = angpos[t];
            if (ta >= 0) {
                for (int va = 0; va < na; ++va) {
                    if (va == ta)
                        continue;
                    const double* o = (ta < va) ? stage + 6 * static_cast<size_t>(tripletIndex(ta, va, na))
                                                : stage + 6 * static_cast<size_t>(tripletIndex(va, ta, na)) + 3;
                    g[0] += o[0];
                    g[1] += o[1];
                    g[2] += o[2];
                }
            }
        } else {
            slotDerivative(t, F, c, sidx, nb, nbv, ang, na, w, nullptr, g);
        }
        for (int gg = 0; gg < 3; ++gg)
            Ci[3 * (t + 1) + gg] = g[gg];
    }
    __syncthreads();
    if (threadIdx.x < 3) {
        double acc = 0.0;
        for (int t = 0; t < nnb; ++t)
            acc += Ci[3 * (t + 1) + threadIdx.x];
        Ci[threadIdx.x] = -acc;
    }
}

/**
 * @brief Gather the gradient of atom a over all centres i whose descriptor depends on it, in
 * ascending i -- the accumulation order of the host code. The neighbour relation is symmetric,
 * so those centres are a itself plus a's own neighbour list. One block per atom: the threads
 * look the contributions up in parallel, three threads add them in list order.
 */
__global__ void kGather(int n, int cap, const int* __restrict__ sidx, const int* __restrict__ nnb_all,
    const int* __restrict__ nb_all, const double* __restrict__ C, double* __restrict__ grad)
{
    __shared__ double part[3 * (CURCUMA_ML_MAX_NEIGHBOURS + 1)];
    const int a = blockIdx.x;
    if (a >= n)
        return;
    const int nna = sidx[a] >= 0 ? nnb_all[a] : 0;
    const int* nba = nb_all + static_cast<size_t>(a) * cap;
    // entry e of the ordered list {a} U nb(a): position of a among its neighbours
    int self = 0;
    while (self < nna && nba[self] < a)
        ++self;
    for (int e = threadIdx.x; e <= nna; e += blockDim.x) {
        const double* ci = nullptr;
        if (e == self) {
            ci = C + static_cast<size_t>(a) * (cap + 1) * 3; // slot 0 of centre a
        } else {
            const int i = nba[e < self ? e : e - 1];
            const int s = findSlot(nb_all + static_cast<size_t>(i) * cap, nnb_all[i], a);
            if (s >= 0)
                ci = C + (static_cast<size_t>(i) * (cap + 1) + s + 1) * 3;
        }
        for (int g = 0; g < 3; ++g)
            part[3 * e + g] = (ci && sidx[a] >= 0) ? ci[g] : 0.0;
    }
    __syncthreads();
    if (threadIdx.x < 3) {
        double acc = 0.0;
        for (int e = 0; e <= nna; ++e)
            acc += part[3 * e + threadIdx.x];
        grad[3 * a + threadIdx.x] = acc;
    }
}

/*! \brief out[0] = sum_i s_i, out[1] = max_i (1 - kmax_i) with 1 for an unknown element. */
__global__ void kFinal(int n, const int* __restrict__ count, const double* __restrict__ s,
    const double* __restrict__ kmax, double* __restrict__ out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0)
        return;
    double e = 0.0, worst = 0.0;
    for (int i = 0; i < n; ++i) {
        e += s[i];
        worst = fmax(worst, count[i] > 0 ? 1.0 - kmax[i] : 1.0);
    }
    out[0] = e;
    out[1] = worst;
}

/**
 * @brief Compact Jacobian per centre on the device: J[i][slot][gamma][f], slot 0 the centre.
 * One thread per neighbour slot writes only its own row, so no atomics are needed.
 */
__global__ void kJacobian(int F, int cap, int capJ, CurcumaMLDeviceConfig c, const int* __restrict__ sidx,
    const int* __restrict__ nnb_all, const int* __restrict__ nb_all, const double* __restrict__ nbv_all,
    double* __restrict__ J)
{
    __shared__ int ang[CURCUMA_ML_MAX_NEIGHBOURS];
    __shared__ int na;
    const int i = blockIdx.x;
    const int nnb = sidx[i] >= 0 ? nnb_all[i] : 0;
    const int* nb = nb_all + static_cast<size_t>(i) * cap;
    const double* nbv = nbv_all + static_cast<size_t>(i) * cap * kNbv;
    double* Ji = J + static_cast<size_t>(i) * (capJ + 1) * 3 * F;
    if (threadIdx.x == 0)
        na = angularList(nnb, nbv, c.rcut_angular, ang);
    __syncthreads();
    for (int t = threadIdx.x; t < nnb; t += blockDim.x) {
        double g[3];
        slotDerivative(t, F, c, sidx, nb, nbv, ang, na, nullptr, Ji + static_cast<size_t>(t + 1) * 3 * F, g);
    }
    __syncthreads();
    for (int e = threadIdx.x; e < 3 * F; e += blockDim.x) {
        double acc = 0.0;
        for (int t = 0; t < nnb; ++t)
            acc += Ji[static_cast<size_t>(t + 1) * 3 * F + e];
        Ji[e] = -acc;
    }
}

/*! \brief Energy row of the design block: phi_m = sum_i k_im, i ascending. */
__global__ void kEnergyRow(int n, int M, const double* __restrict__ k, double* __restrict__ phi)
{
    const int m = blockIdx.x * blockDim.x + threadIdx.x;
    if (m >= M)
        return;
    double acc = 0.0;
    for (int i = 0; i < n; ++i)
        acc += k[static_cast<size_t>(i) * M + m];
    phi[m] = acc;
}

/**
 * @brief Gradient rows of the design block (krr_model.cpp, "Gradient rows"):
 *   PhiG[3a+g][m] = sum_{i ~ a, i in block(m)} sum_f J_i[a,g,f] * Dk_i[m,f],
 *   Dk_i[m,f] = -(k_im / sigma^2) (x_if - z_mf) / scale_f,
 * summed over centres i in ascending order (a itself and a's neighbours).
 */
__global__ void kDesign(int n, int M, int F, int cap, int capJ, double neg_inv_sigma2,
    const int* __restrict__ sidx, const int* __restrict__ nnb_all, const int* __restrict__ nb_all,
    const int* __restrict__ start, const int* __restrict__ count, const double* __restrict__ J,
    const double* __restrict__ Xs, const double* __restrict__ Zt, const double* __restrict__ inv_scale,
    const double* __restrict__ k, double* __restrict__ PhiG)
{
    const int m = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = blockIdx.y;
    if (m >= M || a >= n)
        return;
    double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0;
    if (sidx[a] >= 0) {
        const int nna = nnb_all[a];
        const int* nba = nb_all + static_cast<size_t>(a) * cap;
        bool self_done = false;
        for (int idx = 0; idx <= nna; ++idx) {
            int i, slot;
            if (!self_done && (idx == nna || nba[idx] > a)) {
                i = a;
                slot = 0;
                self_done = true;
                --idx;
            } else {
                if (idx == nna)
                    break;
                i = nba[idx];
                slot = -1;
            }
            if (m < start[i] || m >= start[i] + count[i])
                continue;
            if (slot < 0) {
                const int s = findSlot(nb_all + static_cast<size_t>(i) * cap, nnb_all[i], a);
                if (s < 0)
                    continue;
                slot = s + 1;
            }
            const double kv = k[static_cast<size_t>(i) * M + m];
            const double* Jr = J + (static_cast<size_t>(i) * (capJ + 1) + slot) * 3 * F;
            const double* xi = Xs + static_cast<size_t>(i) * F;
            double s0 = 0.0, s1 = 0.0, s2 = 0.0;
            for (int f = 0; f < F; ++f) {
                const double d = (xi[f] - Zt[static_cast<size_t>(f) * M + m]) * inv_scale[f];
                s0 += Jr[f] * d;
                s1 += Jr[F + f] * d;
                s2 += Jr[2 * F + f] * d;
            }
            acc0 += kv * s0;
            acc1 += kv * s1;
            acc2 += kv * s2;
        }
    }
    PhiG[static_cast<size_t>(3 * a) * M + m] = acc0 * neg_inv_sigma2;
    PhiG[static_cast<size_t>(3 * a + 1) * M + m] = acc1 * neg_inv_sigma2;
    PhiG[static_cast<size_t>(3 * a + 2) * M + m] = acc2 * neg_inv_sigma2;
}

// =================================================================================
// Context
// =================================================================================

struct MLModelContext {
    CurcumaMLDeviceConfig cfg;
    int M = 0;
    int F = 0;
    double sigma = 1.0;
    bool have_beta = false;

    // resident model
    double* d_Z = nullptr;  ///< M x F, standardised centres
    double* d_Zt = nullptr; ///< F x M, the same transposed (coalesced reads in kDesign)
    double* d_mean = nullptr;
    double* d_inv_scale = nullptr;
    double* d_beta = nullptr;

    // per structure, grown on demand
    int cap_atoms = 0;
    int cap_nb = 0;
    double* d_R = nullptr;
    int* d_sidx = nullptr;
    int* d_start = nullptr;
    int* d_count = nullptr;
    int* d_nnb = nullptr;
    int* d_nb = nullptr;
    double* d_nbv = nullptr;
    int* d_flag = nullptr;
    double* d_X = nullptr;
    double* d_k = nullptr;
    double* d_s = nullptr;
    double* d_kmax = nullptr;
    double* d_W = nullptr;
    double* d_C = nullptr;
    double* d_out = nullptr; ///< 2 + 3N: energy, indicator, gradient

    // element ids already on the device, to skip re-uploading them every step
    std::vector<int> h_sidx, h_start, h_count;

    // fit
    cublasHandle_t blas = nullptr;
    double* d_A = nullptr;
    double* d_b = nullptr;
    double* d_J = nullptr;
    size_t cap_J = 0;
    double* d_PhiG = nullptr;
    int cap_rows = 0;
    double* d_phi = nullptr;
    double* d_y = nullptr;
    std::vector<int> h_nnb;
};

template <typename T>
void freeDev(T*& p)
{
    if (p)
        cudaFree(p);
    p = nullptr;
}

template <typename T>
bool allocDev(T*& p, size_t count)
{
    freeDev(p);
    return cudaMalloc(&p, count * sizeof(T)) == cudaSuccess;
}

/*! \brief Grow the per-structure buffers for @p n atoms. */
bool ensureAtoms(MLModelContext* c, int n)
{
    const int cap = std::max(1, std::min(n - 1, CURCUMA_ML_MAX_NEIGHBOURS));
    if (n <= c->cap_atoms && cap <= c->cap_nb)
        return true;
    const int na = std::max(n, c->cap_atoms);
    const int nc = std::max(cap, c->cap_nb);
    const size_t N = static_cast<size_t>(na);
    const bool ok = allocDev(c->d_R, 3 * N) && allocDev(c->d_sidx, N) && allocDev(c->d_start, N)
        && allocDev(c->d_count, N) && allocDev(c->d_nnb, N) && allocDev(c->d_nb, N * nc)
        && allocDev(c->d_nbv, N * nc * kNbv) && allocDev(c->d_X, N * c->F)
        && allocDev(c->d_k, N * c->M) && allocDev(c->d_s, N) && allocDev(c->d_kmax, N)
        && allocDev(c->d_W, N * c->F) && allocDev(c->d_C, N * (nc + 1) * 3)
        && allocDev(c->d_out, 2 + 3 * N);
    c->h_sidx.clear(); // device copies are gone
    if (!ok) {
        c->cap_atoms = c->cap_nb = 0;
        return false;
    }
    c->cap_atoms = na;
    c->cap_nb = nc;
    return true;
}

/*! \brief Upload geometry (always) and element ids (only when they changed). */
bool upload(MLModelContext* c, int n, const double* R, const int* sidx, const int* start, const int* count)
{
    if (cudaMemcpy(c->d_R, R, 3 * static_cast<size_t>(n) * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        return false;
    const bool same = static_cast<int>(c->h_sidx.size()) == n
        && std::equal(sidx, sidx + n, c->h_sidx.begin()) && std::equal(start, start + n, c->h_start.begin())
        && std::equal(count, count + n, c->h_count.begin());
    if (same)
        return true;
    const size_t b = static_cast<size_t>(n) * sizeof(int);
    if (cudaMemcpy(c->d_sidx, sidx, b, cudaMemcpyHostToDevice) != cudaSuccess
        || cudaMemcpy(c->d_start, start, b, cudaMemcpyHostToDevice) != cudaSuccess
        || cudaMemcpy(c->d_count, count, b, cudaMemcpyHostToDevice) != cudaSuccess)
        return false;
    c->h_sidx.assign(sidx, sidx + n);
    c->h_start.assign(start, start + n);
    c->h_count.assign(count, count + n);
    return true;
}

/*! \brief Neighbour lists, standardised descriptor and kernel matrix for the uploaded structure. */
int descriptorAndKernel(MLModelContext* c, int n)
{
    const int cap = std::max(1, std::min(n - 1, CURCUMA_ML_MAX_NEIGHBOURS));
    if (cudaMemset(c->d_flag, 0, sizeof(int)) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    kNeighbours<<<n, kThreads>>>(n, cap, c->cfg.rcut_radial, c->cfg.rcut_angular, c->d_R, c->d_sidx, c->d_nnb, c->d_nb, c->d_nbv, c->d_flag);
    const int tripcap = kStageBytes / descriptorTripletBytes(c->cfg);
    kDescriptor<<<n, kThreads, kStageBytes>>>(c->F, cap, tripcap, c->cfg, c->d_sidx, c->d_nnb, c->d_nb, c->d_nbv, c->d_mean,
        c->d_inv_scale, c->d_X);
    dim3 gridK((c->M + kThreads - 1) / kThreads, n);
    kKernel<<<gridK, kThreads>>>(n, c->M, c->F, 1.0 / (2.0 * c->sigma * c->sigma), c->d_X, c->d_Z,
        c->d_start, c->d_count, c->d_k);
    if (cudaGetLastError() != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    return CURCUMA_ML_OK;
}

int checkOverflow(MLModelContext* c)
{
    int flag = 0;
    if (cudaMemcpy(&flag, c->d_flag, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    return flag ? CURCUMA_ML_ERR_NEIGHBOURS : CURCUMA_ML_OK;
}

} // namespace

extern "C" {

void curcuma_cuda_ml_model_destroy(void* handle);

/**
 * @brief Upload a model: descriptor grid, standardised centres (M x F, row major), feature
 * mean and inverse scale, kernel width. Weights come later (set_beta). nullptr on failure,
 * including "no CUDA device".
 */
void* curcuma_cuda_ml_model_create(const CurcumaMLDeviceConfig* cfg, int M, int F, const double* centers,
    const double* mean, const double* inv_scale, double sigma)
{
    if (!cfg || M <= 0 || F <= 0 || !centers || !mean || !inv_scale)
        return nullptr;
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || devices <= 0)
        return nullptr;

    MLModelContext* c = new MLModelContext();
    c->cfg = *cfg;
    c->M = M;
    c->F = F;
    c->sigma = sigma;

    std::vector<double> zt(static_cast<size_t>(M) * F);
    for (int m = 0; m < M; ++m)
        for (int f = 0; f < F; ++f)
            zt[static_cast<size_t>(f) * M + m] = centers[static_cast<size_t>(m) * F + f];

    const size_t zb = static_cast<size_t>(M) * F * sizeof(double);
    const bool ok = allocDev(c->d_Z, static_cast<size_t>(M) * F) && allocDev(c->d_Zt, static_cast<size_t>(M) * F)
        && allocDev(c->d_mean, F) && allocDev(c->d_inv_scale, F) && allocDev(c->d_beta, M)
        && allocDev(c->d_flag, 1)
        && cudaMemcpy(c->d_Z, centers, zb, cudaMemcpyHostToDevice) == cudaSuccess
        && cudaMemcpy(c->d_Zt, zt.data(), zb, cudaMemcpyHostToDevice) == cudaSuccess
        && cudaMemcpy(c->d_mean, mean, F * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess
        && cudaMemcpy(c->d_inv_scale, inv_scale, F * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
    if (!ok) {
        curcuma_cuda_ml_model_destroy(c);
        return nullptr;
    }
    return c;
}

/*! \brief Upload (or replace) the regression weights, length M. */
int curcuma_cuda_ml_model_set_beta(void* handle, const double* beta)
{
    MLModelContext* c = static_cast<MLModelContext*>(handle);
    if (!c || !beta)
        return CURCUMA_ML_ERR_ARGS;
    if (cudaMemcpy(c->d_beta, beta, static_cast<size_t>(c->M) * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    c->have_beta = true;
    return CURCUMA_ML_OK;
}

/**
 * @brief One prediction.
 * @param R      natoms x 3 coordinates [Angstrom], row major
 * @param sidx   per atom, index into the descriptor species list, -1 for an unknown element
 * @param start  per atom, first centre of its element block
 * @param count  per atom, size of that block (0 if the model has no centres for it)
 * @param out    2 + 3*natoms: kernel energy (no element baseline) [Eh], indicator, and, if
 *               want_grad, the gradient [Eh/Angstrom] row major
 */
int curcuma_cuda_ml_model_predict(void* handle, int natoms, const double* R, const int* sidx,
    const int* start, const int* count, int want_grad, double* out)
{
    MLModelContext* c = static_cast<MLModelContext*>(handle);
    if (!c || natoms <= 0 || !R || !sidx || !start || !count || !out)
        return CURCUMA_ML_ERR_ARGS;
    if (!c->have_beta)
        return CURCUMA_ML_ERR_NOT_READY;
    if (!ensureAtoms(c, natoms))
        return CURCUMA_ML_ERR_ALLOC;
    if (!upload(c, natoms, R, sidx, start, count))
        return CURCUMA_ML_ERR_CUDA;

    const int cap = std::max(1, std::min(natoms - 1, CURCUMA_ML_MAX_NEIGHBOURS));
    int rc = descriptorAndKernel(c, natoms);
    if (rc != CURCUMA_ML_OK)
        return rc;
    kReduce<<<natoms, kThreads>>>(c->M, c->d_k, c->d_beta, c->d_s, c->d_kmax);
    if (want_grad) {
        dim3 gridW((c->F + 31) / 32, natoms);
        kAdjoint<<<gridW, dim3(32, kAdjSplit)>>>(natoms, c->M, c->F, -1.0 / (c->sigma * c->sigma), c->d_X, c->d_Z, c->d_k,
            c->d_beta, c->d_s, c->d_inv_scale, c->d_start, c->d_count, c->d_W);
        kContract<<<natoms, kThreads, kStageBytes>>>(c->F, cap, kStageBytes / kContractTripletBytes, c->cfg, c->d_sidx, c->d_nnb, c->d_nb, c->d_nbv, c->d_W, c->d_C);
        kGather<<<natoms, 64>>>(natoms, cap, c->d_sidx, c->d_nnb, c->d_nb,
            c->d_C, c->d_out + 2);
    }
    kFinal<<<1, 1>>>(natoms, c->d_count, c->d_s, c->d_kmax, c->d_out);
    if (cudaGetLastError() != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;

    rc = checkOverflow(c);
    if (rc != CURCUMA_ML_OK)
        return rc;
    const size_t nout = want_grad ? 2 + 3 * static_cast<size_t>(natoms) : 2;
    if (cudaMemcpy(out, c->d_out, nout * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    return CURCUMA_ML_OK;
}

/*! \brief Allocate and zero the normal equations A (M x M, lower triangle) and b. */
int curcuma_cuda_ml_model_accum_begin(void* handle)
{
    MLModelContext* c = static_cast<MLModelContext*>(handle);
    if (!c)
        return CURCUMA_ML_ERR_ARGS;
    if (!c->blas && cublasCreate(&c->blas) != CUBLAS_STATUS_SUCCESS)
        return CURCUMA_ML_ERR_BLAS;
    const size_t M = static_cast<size_t>(c->M);
    if (!c->d_A && (!allocDev(c->d_A, M * M) || !allocDev(c->d_b, M) || !allocDev(c->d_phi, M)))
        return CURCUMA_ML_ERR_ALLOC;
    if (cudaMemset(c->d_A, 0, M * M * sizeof(double)) != cudaSuccess
        || cudaMemset(c->d_b, 0, M * sizeof(double)) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    return CURCUMA_ML_OK;
}

/**
 * @brief Fold one training structure into the normal equations:
 *   A += wE phi phi^T + wF PhiG^T PhiG,  b += wE phi y_energy + wF PhiG^T y_grad.
 * @param y_grad 3*natoms gradient labels (reference minus base), or nullptr for an energy-only fit
 */
int curcuma_cuda_ml_model_accum_add(void* handle, int natoms, const double* R, const int* sidx,
    const int* start, const int* count, double y_energy, const double* y_grad, double wE, double wF)
{
    MLModelContext* c = static_cast<MLModelContext*>(handle);
    if (!c || natoms <= 0 || !R || !sidx || !start || !count)
        return CURCUMA_ML_ERR_ARGS;
    if (!c->blas || !c->d_A)
        return CURCUMA_ML_ERR_NOT_READY;
    if (!ensureAtoms(c, natoms))
        return CURCUMA_ML_ERR_ALLOC;
    if (!upload(c, natoms, R, sidx, start, count))
        return CURCUMA_ML_ERR_CUDA;
    int rc = descriptorAndKernel(c, natoms);
    if (rc != CURCUMA_ML_OK)
        return rc;
    rc = checkOverflow(c);
    if (rc != CURCUMA_ML_OK)
        return rc;

    const int M = c->M;
    const double one = 1.0;

    // --- energy row ---
    kEnergyRow<<<(M + kThreads - 1) / kThreads, kThreads>>>(natoms, M, c->d_k, c->d_phi);
    const int rows = 3 * natoms;
    if (rows > c->cap_rows) {
        if (!allocDev(c->d_y, rows) || !allocDev(c->d_PhiG, static_cast<size_t>(rows) * M)) {
            c->cap_rows = 0;
            return CURCUMA_ML_ERR_ALLOC;
        }
        c->cap_rows = rows;
    }
    if (cudaMemcpy(c->d_y, &y_energy, sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    if (cublasDsyrk(c->blas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, M, 1, &wE, c->d_phi, M, &one, c->d_A, M)
            != CUBLAS_STATUS_SUCCESS
        || cublasDgemv(c->blas, CUBLAS_OP_N, M, 1, &wE, c->d_phi, M, c->d_y, 1, &one, c->d_b, 1)
            != CUBLAS_STATUS_SUCCESS)
        return CURCUMA_ML_ERR_BLAS;
    if (!y_grad)
        return cudaDeviceSynchronize() == cudaSuccess ? CURCUMA_ML_OK : CURCUMA_ML_ERR_CUDA;

    // --- compact Jacobian, sized by the largest neighbour list of this structure ---
    const int cap = std::max(1, std::min(natoms - 1, CURCUMA_ML_MAX_NEIGHBOURS));
    c->h_nnb.resize(natoms);
    if (cudaMemcpy(c->h_nnb.data(), c->d_nnb, natoms * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    const int capJ = std::max(1, *std::max_element(c->h_nnb.begin(), c->h_nnb.end()));
    const size_t jsize = static_cast<size_t>(natoms) * (capJ + 1) * 3 * c->F;
    if (jsize > c->cap_J) {
        if (!allocDev(c->d_J, jsize)) {
            c->cap_J = 0;
            return CURCUMA_ML_ERR_ALLOC;
        }
        c->cap_J = jsize;
    }
    if (cudaMemset(c->d_J, 0, jsize * sizeof(double)) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    kJacobian<<<natoms, kThreads>>>(c->F, cap, capJ, c->cfg, c->d_sidx, c->d_nnb, c->d_nb, c->d_nbv, c->d_J);

    // --- gradient rows ---
    dim3 gridD((M + kThreads - 1) / kThreads, natoms);
    kDesign<<<gridD, kThreads>>>(natoms, M, c->F, cap, capJ, -1.0 / (c->sigma * c->sigma), c->d_sidx, c->d_nnb,
        c->d_nb, c->d_start, c->d_count, c->d_J, c->d_X, c->d_Zt, c->d_inv_scale, c->d_k, c->d_PhiG);
    if (cudaGetLastError() != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    if (cudaMemcpy(c->d_y, y_grad, rows * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    if (cublasDsyrk(c->blas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, M, rows, &wF, c->d_PhiG, M, &one, c->d_A, M)
            != CUBLAS_STATUS_SUCCESS
        || cublasDgemv(c->blas, CUBLAS_OP_N, M, rows, &wF, c->d_PhiG, M, c->d_y, 1, &one, c->d_b, 1)
            != CUBLAS_STATUS_SUCCESS)
        return CURCUMA_ML_ERR_BLAS;
    return CURCUMA_ML_OK;
}

/*! \brief Copy A (lower triangle, column major) and b back. */
int curcuma_cuda_ml_model_accum_fetch(void* handle, double* A, double* b)
{
    MLModelContext* c = static_cast<MLModelContext*>(handle);
    if (!c || !A || !b)
        return CURCUMA_ML_ERR_ARGS;
    if (!c->d_A)
        return CURCUMA_ML_ERR_NOT_READY;
    if (cudaDeviceSynchronize() != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    const size_t M = static_cast<size_t>(c->M);
    if (cudaMemcpy(A, c->d_A, M * M * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
        || cudaMemcpy(b, c->d_b, M * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess)
        return CURCUMA_ML_ERR_CUDA;
    // The M x M system is not needed for prediction; give the memory back.
    freeDev(c->d_A);
    freeDev(c->d_b);
    freeDev(c->d_phi);
    freeDev(c->d_J);
    freeDev(c->d_PhiG);
    freeDev(c->d_y);
    c->cap_J = 0;
    c->cap_rows = 0;
    return CURCUMA_ML_OK;
}

void curcuma_cuda_ml_model_destroy(void* handle)
{
    MLModelContext* c = static_cast<MLModelContext*>(handle);
    if (!c)
        return;
    freeDev(c->d_Z);
    freeDev(c->d_Zt);
    freeDev(c->d_mean);
    freeDev(c->d_inv_scale);
    freeDev(c->d_beta);
    freeDev(c->d_R);
    freeDev(c->d_sidx);
    freeDev(c->d_start);
    freeDev(c->d_count);
    freeDev(c->d_nnb);
    freeDev(c->d_nb);
    freeDev(c->d_nbv);
    freeDev(c->d_flag);
    freeDev(c->d_X);
    freeDev(c->d_k);
    freeDev(c->d_s);
    freeDev(c->d_kmax);
    freeDev(c->d_W);
    freeDev(c->d_C);
    freeDev(c->d_out);
    freeDev(c->d_A);
    freeDev(c->d_b);
    freeDev(c->d_J);
    freeDev(c->d_PhiG);
    freeDev(c->d_phi);
    freeDev(c->d_y);
    if (c->blas)
        cublasDestroy(c->blas);
    delete c;
}

} // extern "C"
