/*
 * <Native QM: point kernels of the semi-local ECP integrals, shared by the CPU (qm_ecp.cpp) and the GPU
 *  (cuda/qm_gpu_ecp.cu)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ONE source of the ECP building blocks for both sides (the pattern of qm_xc_kernels.h): scaled modified spherical
 * Bessel functions, binomials, Gaussian moments, the analytic three-centre overlap of the local channel, the
 * Gauss-Chebyshev radial point and the monomial index of the angular tables. Formulas and literature: qm_ecp.h.
 * A change here changes both sides; ctest -L qm_ecp (CPU vs PySCF/ORCA) and -L qm_gpu_ecp (device vs CPU) cover it.
 *
 * Plain C math from <cmath> (global namespace), no Eigen, no std containers: host code (g++) and __host__
 * __device__ code (nvcc). Device math may differ from glibc by an ulp, so the GPU agrees with the CPU to rounding,
 * not bitwise.
 *
 * Claude Generated (Nov 2026): functions moved here from qm_ecp.cpp with the same arithmetic (the CPU ECP energies
 * are unchanged, checked by md5 of the qm_ecp dumps). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include <cmath>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define QM_ECP_HD __host__ __device__ inline
#else
#define QM_ECP_HD inline
#endif

namespace qmecp {
namespace kern {

constexpr double kPi = 3.14159265358979323846;
constexpr double kExpCut = 46.0;  ///< e^-46 ~ 1e-20: terms below are dropped (same cut on both sides)
constexpr int kMaxBesselL = 16;    ///< largest lambda the device tables use (n_max 5 + l_proj 5 + margin)

/// One cartesian primitive c (x-Ax)^p0 (y-Ay)^p1 (z-Az)^p2 exp(-alpha |r-A|^2) of a function on centre A.
struct Prim {
    double c, alpha;
    int p[3];
};

/// e^(-x) i_l(x), l = 0..lmax (modified spherical Bessel function of the first kind, scaled). Series for x <= 50
/// (positive terms, no cancellation), the exact finite form above:
///   i_l(x) = 1/(2x) [ e^x sum_k (-1)^k a_k (2x)^-k + (-1)^(l+1) e^-x sum_k a_k (2x)^-k ],  a_k = (l+k)!/(k!(l-k)!)
QM_ECP_HD void scaledBessel(double x, int lmax, double* out)
{
    if (x < 1.0e-14) {
        out[0] = 1.0;
        for (int l = 1; l <= lmax; ++l) out[l] = 0.0;
        return;
    }
    if (x <= 50.0) {
        const double ex = exp(-x), h = 0.5 * x * x;
        double xl = 1.0, df = 1.0;  // x^l, (2l+1)!!
        for (int l = 0; l <= lmax; ++l) {
            if (l > 0) {
                xl *= x;
                df *= (2 * l + 1);
            }
            double t = xl / df, s = t;
            for (int k = 0; k < 2000; ++k) {
                t *= h / ((k + 1.0) * (2.0 * l + 2.0 * k + 3.0));
                s += t;
                if (t < 1.0e-17 * s) break;
            }
            out[l] = ex * s;
        }
        return;
    }
    const double e2 = exp(-2.0 * x), inv = 1.0 / (2.0 * x);
    for (int l = 0; l <= lmax; ++l) {
        double s1 = 0.0, s2 = 0.0, a = 1.0, p = 1.0;  // a_k, (2x)^-k
        for (int k = 0; k <= l; ++k) {
            if (k > 0) {
                a *= (double)(l + k) * (l - k + 1) / k;
                p *= inv;
            }
            s1 += ((k & 1) ? -a : a) * p;
            s2 += a * p;
        }
        out[l] = inv * (s1 + ((l & 1) ? 1.0 : -1.0) * e2 * s2);
    }
}

QM_ECP_HD double binom(int n, int k)
{
    double b = 1.0;
    for (int i = 1; i <= k; ++i) b = b * (n - k + i) / i;
    return b;
}

/// int t^k exp(-g t^2) dt over R
QM_ECP_HD double gaussMoment(int k, double g)
{
    if (k & 1) return 0.0;
    double df = 1.0;
    for (int i = k - 1; i > 1; i -= 2) df *= i;
    return df / pow(2.0 * g, (double)(k / 2)) * sqrt(kPi / g);
}

/// <a| exp(-zeta |r-C|^2) |b> of two primitives (without their coefficients), analytic Gaussian product.
QM_ECP_HD double threeCentreOverlap(const Prim& a, const double* A, const Prim& b, const double* B, double zeta,
                                    const double* C)
{
    const double g = a.alpha + b.alpha + zeta;
    double ab2 = 0, ac2 = 0, bc2 = 0, Q[3];
    for (int d = 0; d < 3; ++d) {
        ab2 += (A[d] - B[d]) * (A[d] - B[d]);
        ac2 += (A[d] - C[d]) * (A[d] - C[d]);
        bc2 += (B[d] - C[d]) * (B[d] - C[d]);
        Q[d] = (a.alpha * A[d] + b.alpha * B[d] + zeta * C[d]) / g;
    }
    const double ex = (a.alpha * b.alpha * ab2 + a.alpha * zeta * ac2 + b.alpha * zeta * bc2) / g;
    if (ex > kExpCut) return 0.0;
    double v = exp(-ex);
    for (int d = 0; d < 3; ++d) {
        double s = 0.0;
        const double qa = Q[d] - A[d], qb = Q[d] - B[d];
        for (int i = 0; i <= a.p[d]; ++i)
            for (int j = 0; j <= b.p[d]; ++j) {
                if ((i + j) & 1) continue;
                s += binom(a.p[d], i) * binom(b.p[d], j) * pow(qa, (double)(a.p[d] - i)) *
                     pow(qb, (double)(b.p[d] - j)) * gaussMoment(i + j, g);
            }
        v *= s;
    }
    return v;
}

/// Point i (1-based) of the n-point Gauss-Chebyshev (2nd kind) radial grid on r = ln(2/(1-x))/ln 2; w includes
/// dr/dx (not r^2). Grids with n = 2^k - 1 are nested: point i of n equals point 2i of 2n+1 bit for bit.
QM_ECP_HD void radialPoint(int i, int n, double& r, double& w)
{
    const double ln2 = log(2.0);
    const double th = i * kPi / (n + 1), x = cos(th);
    r = log(2.0 / (1.0 - x)) / ln2;
    w = kPi / (n + 1) * sin(th) / ((1.0 - x) * ln2);
}

/// Index of the monomial x^i y^j z^k (t = i + j + k) in the order of qmecp's MonoIndex: totals ascending, within a
/// total i descending, then j descending. Independent of the largest total.
QM_ECP_HD int monoIndex(int i, int j, int k)
{
    const int t = i + j + k;
    return t * (t + 1) * (t + 2) / 6 + (t - i) * (t - i + 1) / 2 + (t - i - j);
}

}  // namespace kern
}  // namespace qmecp
