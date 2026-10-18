/*
 * <Grid-free integration of a mean-force profile dG/ds by a natural cubic spline.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 * Header-only so that NebMD::writeMeanForceProfile and the unit test share one implementation.
 * Also holds the absolute TST rate / Eyring free energy across s* (tstRate, zHalfNear).
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "src/core/units.h" // MD time unit, Planck constant, Hartree energy (TST rate)

namespace curcuma {
namespace fe_profile {

/*! \brief Integrated free-energy profile on a fine grid (Claude Generated 2026).
 *
 * All energies in the unit of dgds * s (Eh when dG/ds is in Eh per unit s). G is zero at
 * the first support point. Two barrier definitions are returned:
 *  - barrier_forward   = max G - G(first support point)  (the NebMD report convention)
 *  - barrier_mean_ends = max G - (G(first) + G(last)) / 2
 */
struct SplineProfile {
    std::vector<double> s;       ///< fine grid (ngrid points from first to last support point)
    std::vector<double> G;       ///< G on the fine grid
    std::vector<double> G_nodes; ///< G at the support points (same integration)
    double gmax = 0.0;           ///< maximum of G over the fine grid and the support points
    double s_max = 0.0;          ///< position of that maximum
    double g_first = 0.0, g_last = 0.0;
    double barrier_forward = 0.0;
    double barrier_mean_ends = 0.0;
    double reaction = 0.0;       ///< G(last) - G(first)
    bool spline = false;         ///< false: fallback to the trapezoid rule on the support points
    std::string warning;         ///< non-empty if the fallback was used
};

/*! \brief Decomposition of the conditional mean force at one support point (Claude Generated 2026).
 *
 * With a z wall W the sampled potential is E + W and the estimator (I) is
 *     dG/ds = <b.grad(E + W)> - kT <div b> = grad_e + metric + wall,
 *     grad_e = <b.grad E>, metric = -kT <div b>, wall = <b.grad W>.
 * Inputs are the weighted sums over the samples of one support point: sum_wf = sum w b.grad(E+W),
 * sum_wwall = sum w b.grad W, sum_wd = sum w div b, sum_w = sum w (w = Z^-1/2 constrained, 1 restrained).
 */
struct MeanForceSplit {
    double grad_e = 0.0; ///< <b.grad E>, E = physical potential only [Eh per unit s]
    double metric = 0.0; ///< -kT <div b>                             [Eh per unit s]
    double wall = 0.0;   ///< <b.grad W>, W = z wall                    [Eh per unit s]
};

inline MeanForceSplit splitMeanForce(double sum_wf, double sum_wwall, double sum_wd, double sum_w, double kT)
{
    MeanForceSplit r;
    r.grad_e = (sum_wf - sum_wwall) / sum_w;
    r.metric = -kT * sum_wd / sum_w;
    r.wall = sum_wwall / sum_w;
    return r;
}

/*! \brief Trapezoid integral of y over x at the support points (cumulative, starts at 0). */
inline std::vector<double> integrateTrapezoid(const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> G(x.size(), 0.0);
    for (size_t k = 1; k < x.size(); ++k)
        G[k] = G[k - 1] + 0.5 * (y[k] + y[k - 1]) * (x[k] - x[k - 1]);
    return G;
}

/*! \brief Second derivatives M_k of the natural cubic spline through (x_k, y_k), M_0 = M_n-1 = 0.
 *
 * Tridiagonal system h_{k-1} M_{k-1} + 2(h_{k-1}+h_k) M_k + h_k M_{k+1}
 *   = 6 [(y_{k+1}-y_k)/h_k - (y_k-y_{k-1})/h_{k-1}], solved by the Thomas algorithm.
 */
inline std::vector<double> naturalSplineM(const std::vector<double>& x, const std::vector<double>& y)
{
    const size_t n = x.size();
    std::vector<double> M(n, 0.0);
    if (n < 3)
        return M;
    const size_t m = n - 2; // unknowns M_1..M_{n-2}
    std::vector<double> a(m), b(m), c(m), d(m);
    for (size_t i = 0; i < m; ++i) {
        const size_t k = i + 1;
        const double h0 = x[k] - x[k - 1], h1 = x[k + 1] - x[k];
        a[i] = h0;
        b[i] = 2.0 * (h0 + h1);
        c[i] = h1;
        d[i] = 6.0 * ((y[k + 1] - y[k]) / h1 - (y[k] - y[k - 1]) / h0);
    }
    for (size_t i = 1; i < m; ++i) {
        const double w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    M[m] = d[m - 1] / b[m - 1];
    for (size_t i = m - 1; i-- > 0;)
        M[i + 1] = (d[i] - c[i] * M[i + 2]) / b[i];
    return M;
}

/*! \brief Integrate dG/ds through a natural cubic spline, exactly per interval.
 *
 * On [s_k, s_k+1] with h = s_k+1 - s_k and t = s - s_k the spline is
 *   y(t) = y_k + b t + M_k/2 t^2 + (M_k+1 - M_k)/(6h) t^3,  b = (y_k+1 - y_k)/h - h (2 M_k + M_k+1)/6,
 * and G(t) = G_k + y_k t + b t^2/2 + M_k t^3/6 + (M_k+1 - M_k) t^4/(24 h) is its exact integral.
 * G is evaluated on ngrid equidistant points between the first and last support point (and
 * at the support points); max G is taken over both. Removes the systematic underestimate of
 * the barrier when the maximum falls between two support points.
 * With fewer than 3 support points, or non-increasing s after sorting, the trapezoid values
 * at the support points are used instead and `warning` is set.
 *
 * \param s        support points (sorted internally)
 * \param dgds_eh  dG/ds at the support points (Eh per unit s)
 * \param ngrid    fine grid size (>= 2)
 */
inline SplineProfile integrateSpline(const std::vector<double>& s, const std::vector<double>& dgds_eh, int ngrid = 8001)
{
    SplineProfile r;
    const size_t n = std::min(s.size(), dgds_eh.size());
    if (n == 0) {
        r.warning = "no support points";
        return r;
    }
    std::vector<std::pair<double, double>> p(n);
    for (size_t k = 0; k < n; ++k)
        p[k] = { s[k], dgds_eh[k] };
    std::sort(p.begin(), p.end());
    std::vector<double> x(n), y(n);
    for (size_t k = 0; k < n; ++k) {
        x[k] = p[k].first;
        y[k] = p[k].second;
    }
    bool distinct = true;
    for (size_t k = 1; k < n; ++k)
        if (!(x[k] - x[k - 1] > 1e-14))
            distinct = false;

    auto finish = [&r]() {
        r.g_first = r.G_nodes.front();
        r.g_last = r.G_nodes.back();
        r.barrier_forward = r.gmax - r.g_first;
        r.barrier_mean_ends = r.gmax - 0.5 * (r.g_first + r.g_last);
        r.reaction = r.g_last - r.g_first;
    };

    if (n < 3 || !distinct) {
        r.warning = n < 3 ? "fewer than 3 support points: spline not defined, using the trapezoid values at the support points"
                          : "repeated s values: spline not defined, using the trapezoid values at the support points";
        r.G_nodes = integrateTrapezoid(x, y);
        r.s = x;
        r.G = r.G_nodes;
        size_t km = 0;
        for (size_t k = 1; k < n; ++k)
            if (r.G_nodes[k] > r.G_nodes[km])
                km = k;
        r.gmax = r.G_nodes[km];
        r.s_max = x[km];
        finish();
        return r;
    }

    r.spline = true;
    const std::vector<double> M = naturalSplineM(x, y);
    auto segIntegral = [&](size_t k, double t) {
        const double h = x[k + 1] - x[k];
        const double b = (y[k + 1] - y[k]) / h - h * (2.0 * M[k] + M[k + 1]) / 6.0;
        return y[k] * t + 0.5 * b * t * t + M[k] * t * t * t / 6.0 + (M[k + 1] - M[k]) * t * t * t * t / (24.0 * h);
    };
    r.G_nodes.assign(n, 0.0);
    for (size_t k = 0; k + 1 < n; ++k)
        r.G_nodes[k + 1] = r.G_nodes[k] + segIntegral(k, x[k + 1] - x[k]);

    ngrid = std::max(ngrid, 2);
    r.s.resize(ngrid);
    r.G.resize(ngrid);
    size_t k = 0;
    for (int i = 0; i < ngrid; ++i) {
        const double si = (i == ngrid - 1) ? x.back() : x.front() + (x.back() - x.front()) * i / (ngrid - 1.0);
        while (k + 2 < n && si > x[k + 1])
            ++k;
        r.s[i] = si;
        r.G[i] = r.G_nodes[k] + segIntegral(k, si - x[k]);
    }
    r.gmax = r.G_nodes[0];
    r.s_max = x[0];
    for (size_t j = 0; j < n; ++j)
        if (r.G_nodes[j] > r.gmax) {
            r.gmax = r.G_nodes[j];
            r.s_max = x[j];
        }
    for (int i = 0; i < ngrid; ++i)
        if (r.G[i] > r.gmax) {
            r.gmax = r.G[i];
            r.s_max = r.s[i];
        }
    finish();
    return r;
}

/*! \brief Barrier and reaction free energy referenced to the basin free energies (Claude Generated 2026).
 *
 *  s* = position of the highest INTERIOR local maximum of G (a point not at either end with G not
 *  below its two neighbours; searched over s in [search_lo, search_hi], default: all points).
 *  Interior, not global: a profile that extends beyond the basins can rise at its edges above the
 *  barrier, and a global maximum would then sit on that wall. The two basins are A = {s <= s*} and B = {s >= s*} (all points of the profile on
 *  that side, including s*), and their free energies are
 *
 *      G_A = -kT ln int_A exp(-G(s)/kT) ds ,   G_B = -kT ln int_B exp(-G(s)/kT) ds
 *
 *  (trapezoid rule over the given points; G shifted by its minimum before exponentiation).
 *  barrier_fwd = G(s*) - G_A, barrier_rev = G(s*) - G_B, reaction = G_B - G_A, hence
 *  barrier_fwd - barrier_rev = reaction identically.
 *
 *  Caveat (scale): ds carries the unit of the CV coordinate s, so G_A and G_B - and the basin
 *  barriers - shift by -kT ln(c) if s is rescaled by c. The absolute value therefore depends on
 *  the CV parametrisation (frame count, extension, lambda); only compare basin barriers computed
 *  with the same definition on the same CV. The reaction dG is affected only through the
 *  relative basin widths, not by a global rescaling.
 *
 *  G and kT must share one unit (Eh or kcal/mol). s must be ascending (sorted internally
 *  otherwise). valid = false (with `warning`) when there is no interior local maximum (monotonic
 *  profile) or fewer than 3 points are given. */
struct BasinBarrier {
    bool valid = false;
    double s_star = 0.0, g_star = 0.0;
    double G_A = 0.0, G_B = 0.0;
    double barrier_fwd = 0.0, barrier_rev = 0.0, reaction = 0.0;
    std::string warning;
};

inline BasinBarrier basinBarrier(const std::vector<double>& s_in, const std::vector<double>& G_in, double kT,
                                 double search_lo = -1e300, double search_hi = 1e300)
{
    BasinBarrier r;
    const size_t n = std::min(s_in.size(), G_in.size());
    if (n < 3 || !(kT > 0.0)) {
        r.warning = "basin barrier needs at least 3 points and kT > 0";
        return r;
    }
    std::vector<std::pair<double, double>> p(n);
    for (size_t k = 0; k < n; ++k)
        p[k] = { s_in[k], G_in[k] };
    std::sort(p.begin(), p.end());
    // s* = highest INTERIOR local maximum (G[k] >= both neighbours, not the first/last point).
    // A global maximum would land on a wall at the edge whenever the profile extends beyond the
    // basins and rises there above the barrier (measured on the unit-test double well).
    size_t kmax = n;
    double gmin = p[0].second;
    for (size_t k = 0; k < n; ++k)
        gmin = std::min(gmin, p[k].second);
    for (size_t k = 1; k + 1 < n; ++k) {
        if (p[k].first < search_lo || p[k].first > search_hi)
            continue;
        if (p[k].second < p[k - 1].second || p[k].second < p[k + 1].second)
            continue;
        if (kmax == n || p[k].second > p[kmax].second)
            kmax = k;
    }
    if (kmax == n) {
        r.warning = "no interior maximum of G (monotonic profile or maximum at the edge): no barrier between two basins";
        return r;
    }
    r.s_star = p[kmax].first;
    r.g_star = p[kmax].second;
    auto basinIntegral = [&](size_t k0, size_t k1) {
        double I = 0.0;
        for (size_t k = k0 + 1; k <= k1; ++k)
            I += 0.5 * (std::exp(-(p[k].second - gmin) / kT) + std::exp(-(p[k - 1].second - gmin) / kT))
                * (p[k].first - p[k - 1].first);
        return I;
    };
    const double IA = basinIntegral(0, kmax), IB = basinIntegral(kmax, n - 1);
    if (!(IA > 0.0) || !(IB > 0.0)) {
        r.warning = "empty basin integral";
        return r;
    }
    r.G_A = gmin - kT * std::log(IA);
    r.G_B = gmin - kT * std::log(IB);
    r.barrier_fwd = r.g_star - r.G_A;
    r.barrier_rev = r.g_star - r.G_B;
    r.reaction = r.G_B - r.G_A;
    r.valid = true;
    return r;
}

/*! \brief Does the profile reach the basin bottoms? (Claude Generated 2026)
 *
 *  Criterion on the slope at the outermost support points, relative to the largest slope:
 *  the left basin is NOT covered if dG/ds(first) > frac * max|dG/ds| (G still rising steeply
 *  towards the barrier, the bottom lies further left); the right basin is NOT covered if
 *  -dG/ds(last) > frac * max|dG/ds|. A slope of the opposite sign at the edge means the
 *  bottom has been passed and counts as covered. s ascending (sorted internally). */
struct BasinCoverage {
    bool left = true, right = true;
    double slope_left = 0.0, slope_right = 0.0, max_abs = 0.0;
};

inline BasinCoverage basinCoverage(const std::vector<double>& s, const std::vector<double>& dgds, double frac = 0.5)
{
    BasinCoverage c;
    const size_t n = std::min(s.size(), dgds.size());
    if (n == 0)
        return c;
    std::vector<std::pair<double, double>> p(n);
    for (size_t k = 0; k < n; ++k)
        p[k] = { s[k], dgds[k] };
    std::sort(p.begin(), p.end());
    for (const auto& q : p)
        c.max_abs = std::max(c.max_abs, std::abs(q.second));
    c.slope_left = p.front().second;
    c.slope_right = p.back().second;
    c.left = !(c.slope_left > frac * c.max_abs);
    c.right = !(-c.slope_right > frac * c.max_abs);
    return c;
}

/*! \brief dG/ds at the points of a sampled curve G(s) by finite differences (Claude Generated 2026).
 *  One-sided at the ends, central (non-uniform) inside; s ascending. Used for the coverage test of
 *  a binned reference profile, where no mean force is available. */
inline std::vector<double> finiteDifferenceSlope(const std::vector<double>& s, const std::vector<double>& G)
{
    const size_t n = std::min(s.size(), G.size());
    std::vector<double> d(n, 0.0);
    if (n < 2)
        return d;
    d[0] = (G[1] - G[0]) / (s[1] - s[0]);
    d[n - 1] = (G[n - 1] - G[n - 2]) / (s[n - 1] - s[n - 2]);
    for (size_t k = 1; k + 1 < n; ++k)
        d[k] = (G[k + 1] - G[k - 1]) / (s[k + 1] - s[k - 1]);
    return d;
}

/*! \brief Bins used by the zHalfNear fallback (Claude Generated 2026).
 *
 *  Left side = populated bins with s_b <= s*, right side = s_b > s*. On each side the bins are taken
 *  in order of increasing distance from s* until the side holds at least ceil(min_samples/2)
 *  samples (or the side is exhausted). Exposed so that a caller can evaluate an expensive per-bin
 *  quantity only for the bins that will be used. */
inline std::vector<size_t> zHalfFallbackBins(const std::vector<double>& s_b, const std::vector<long>& n_b,
                                             double s_star, long min_samples)
{
    const size_t n = std::min(s_b.size(), n_b.size());
    const long need = (min_samples + 1) / 2;
    std::vector<std::pair<double, size_t>> left, right;
    for (size_t b = 0; b < n; ++b) {
        if (n_b[b] <= 0)
            continue;
        if (s_b[b] <= s_star)
            left.push_back({ s_star - s_b[b], b });
        else
            right.push_back({ s_b[b] - s_star, b });
    }
    std::sort(left.begin(), left.end());
    std::sort(right.begin(), right.end());
    std::vector<size_t> sel;
    for (const auto* side : { &left, &right }) {
        long got = 0;
        for (const auto& q : *side) {
            if (got >= need && got > 0)
                break;
            sel.push_back(q.second);
            got += n_b[q.second];
        }
    }
    return sel;
}

/*! \brief Mean <Z^{1/2}> near s* from per-bin estimates (Claude Generated 2026).
 *
 *  Inputs per bin b: position s_b (mean s of the bin, or its midpoint), sample count n_b and the
 *  bin estimate v_b of <Z^{1/2}>_s (constrained ensemble: v_b = n_b / sum_i w_i with w = Z^{-1/2};
 *  free MD: v_b = plain mean of Z^{1/2}). Bins with n_b <= 0 are unpopulated.
 *
 *  Window rule: sample-weighted mean of v_b over all populated bins with |s_b - s*| <= half_width.
 *  If the window holds fewer than min_samples samples, `fallback` is set (with `warning`) and the
 *  value is the linear interpolation in s between the sample-weighted means of the two sides
 *  (zHalfFallbackBins: nearest populated bins per side until ceil(min_samples/2) samples), or the
 *  one side's mean if s* lies outside the populated range. */
struct ZHalfEstimate {
    bool valid = false;
    double value = 0.0;   ///< <Z^{1/2}> in 1/(amu^{1/2} A) per unit s
    long samples = 0;     ///< samples that entered the value
    int bins = 0;         ///< bins that entered the value
    bool fallback = false;
    std::string warning;
};

inline ZHalfEstimate zHalfNear(const std::vector<double>& s_b, const std::vector<long>& n_b,
                               const std::vector<double>& v_b, double s_star, double half_width,
                               long min_samples)
{
    ZHalfEstimate r;
    const size_t n = std::min(s_b.size(), std::min(n_b.size(), v_b.size()));
    double num = 0.0;
    long cnt = 0;
    int nb = 0;
    for (size_t b = 0; b < n; ++b) {
        if (n_b[b] <= 0 || !std::isfinite(v_b[b]) || std::abs(s_b[b] - s_star) > half_width)
            continue;
        num += static_cast<double>(n_b[b]) * v_b[b];
        cnt += n_b[b];
        ++nb;
    }
    if (cnt >= min_samples && cnt > 0) {
        r.valid = true;
        r.value = num / static_cast<double>(cnt);
        r.samples = cnt;
        r.bins = nb;
        return r;
    }
    // fallback: per side, the nearest populated bins (growing outward) until that side holds
    // min_samples/2 samples; interpolate linearly between the two side means.
    const std::vector<size_t> sel = zHalfFallbackBins(s_b, n_b, s_star, min_samples);
    double sl = 0.0, vl = 0.0, sr = 0.0, vr = 0.0;
    long nl = 0, nr = 0;
    int bl = 0, br = 0;
    for (size_t b : sel) {
        if (!std::isfinite(v_b[b]))
            continue;
        const double w = static_cast<double>(n_b[b]);
        if (s_b[b] <= s_star) {
            sl += w * s_b[b]; vl += w * v_b[b]; nl += n_b[b]; ++bl;
        } else {
            sr += w * s_b[b]; vr += w * v_b[b]; nr += n_b[b]; ++br;
        }
    }
    r.fallback = true;
    if (nl == 0 && nr == 0) {
        r.warning = "no populated bin: <Z^1/2> at s* not defined";
        return r;
    }
    if (nl > 0 && nr > 0) {
        sl /= nl; vl /= nl; sr /= nr; vr /= nr;
        const double t = (sr - sl) > 1e-15 ? (s_star - sl) / (sr - sl) : 0.5;
        r.value = (1.0 - t) * vl + t * vr;
    } else {
        r.value = nl > 0 ? vl / nl : vr / nr;
    }
    r.samples = nl + nr;
    r.bins = bl + br;
    r.valid = true;
    r.warning = "only " + std::to_string(cnt) + " samples within s* +- " + std::to_string(half_width)
        + " (minimum " + std::to_string(min_samples) + "): <Z^1/2> interpolated between the sample-weighted means of the"
        " nearest populated bins on each side (>= " + std::to_string((min_samples + 1) / 2) + " samples per side if available)";
    return r;
}

/*! \brief Absolute transition-state-theory rate across s* and its Eyring free energy (Claude Generated 2026).
 *
 *  Derivation (classical TST on a coordinate s(R); checked by hand, see test_tst_rate.cpp):
 *   - sdot = sum_a grad_a s . v_a with Maxwell-Boltzmann v_a ~ N(0, kT/m_a), independent of R, so for
 *     fixed R sdot is Gaussian with mean 0 and variance kT Z(R), Z = sum_a |grad_a s|^2 / m_a.
 *   - E|x| = sigma sqrt(2/pi) for x ~ N(0, sigma^2), hence <|sdot|>_s = sqrt(2 kT / pi) <Z^{1/2}>_s,
 *     averaged over the TRUE conditional ensemble at fixed s. The one-way flux <sdot theta(sdot)> is
 *     half of it.
 *   - k_A->B = <delta(s - s*) sdot theta(sdot)> / <theta(s* - s)> = (1/2) <|sdot|>_{s*} p_A(s*),
 *     p_A(s*) = exp(-G(s*)/kT) / int_{s<=s*} exp(-G/kT) ds, with G = -kT ln P(s) the free energy of
 *     the MARGINAL density of s (what the Z^{-1/2}-reweighted constrained mean force integrates to).
 *     B -> A uses p_B over s >= s*.
 *   - Constrained (blue-moon) sampling has density proportional to Z^{1/2} times the conditional one,
 *     which is why every conditional average is reweighted with w = Z^{-1/2}:
 *     <Z^{1/2}>_s = sum_i w_i Z_i^{1/2} / sum_i w_i = N / sum_i w_i (since w_i Z_i^{1/2} = 1).
 *   - Scale invariance: for s' = a s, p' = p / a and <Z'^{1/2}> = a <Z^{1/2}>, so k is unchanged
 *     (unlike the basin barrier, which shifts by kT ln a).
 *   - Eyring: dG_TST = kT ln( (kB T / h) / k_TST ).
 *  Units: kT in Eh, G in Eh, s dimensionless or any unit, zhalf in s / (amu^{1/2} A). Then
 *  sqrt(kT) zhalf is in s per MD time unit sqrt(amu A^2 / Eh) = MD_TIME_UNIT_FS fs; converted to ps.
 *  The basin integrals use the trapezoid rule over the given points (G shifted by its minimum), with
 *  s* inserted by linear interpolation of G if it is not a grid point. */
struct TstRate {
    bool valid = false;
    double s_star = 0.0, g_star = 0.0;   ///< s* and G(s*) (Eh)
    double p_fwd = 0.0, p_rev = 0.0;      ///< p_A(s*), p_B(s*) in 1/(unit s)
    double sdot_abs = 0.0;                ///< <|sdot|>(s*) in s per ps
    double k_fwd = 0.0, k_rev = 0.0;      ///< rates in 1/ps
    double eyring_prefactor = 0.0;        ///< kB T / h in 1/ps
    double dG_fwd = 0.0, dG_rev = 0.0;    ///< Eyring activation free energies (Eh)
    std::string warning;
};

/*! \brief kB T / h in 1/ps for kT in Eh (Claude Generated 2026). */
inline double eyringPrefactorPerPs(double kT_eh)
{
    return kT_eh * CurcumaUnit::Constants::HARTREE_TO_JOULE / CurcumaUnit::Constants::PLANCK_SI * 1e-12;
}

/*! \brief Eyring activation free energy dG = kT ln((kB T/h) / k) in Eh; k in 1/ps, kT in Eh (Claude Generated 2026). */
inline double eyringFreeEnergy(double k_per_ps, double kT_eh)
{
    return kT_eh * std::log(eyringPrefactorPerPs(kT_eh) / k_per_ps);
}

inline TstRate tstRate(const std::vector<double>& s_in, const std::vector<double>& G_in, double kT_eh,
                       double s_star, double zhalf)
{
    TstRate r;
    const size_t n = std::min(s_in.size(), G_in.size());
    if (n < 2 || !(kT_eh > 0.0) || !(zhalf > 0.0) || !std::isfinite(zhalf)) {
        r.warning = "TST rate needs at least 2 points, kT > 0 and <Z^1/2> > 0";
        return r;
    }
    std::vector<std::pair<double, double>> p(n);
    for (size_t k = 0; k < n; ++k)
        p[k] = { s_in[k], G_in[k] };
    std::sort(p.begin(), p.end());
    if (s_star <= p.front().first || s_star >= p.back().first) {
        r.warning = "s* is not inside the profile";
        return r;
    }
    double gmin = p[0].second;
    for (const auto& q : p)
        gmin = std::min(gmin, q.second);
    size_t j = 1;
    while (j < n - 1 && p[j].first < s_star)
        ++j; // p[j-1].first < s_star <= p[j].first
    const double h = p[j].first - p[j - 1].first;
    const double t = h > 0.0 ? (s_star - p[j - 1].first) / h : 0.0;
    const double gs = (1.0 - t) * p[j - 1].second + t * p[j].second;
    auto e = [&](double g) { return std::exp(-(g - gmin) / kT_eh); };
    double IA = 0.0, IB = 0.0;
    for (size_t k = 1; k < j; ++k)
        IA += 0.5 * (e(p[k].second) + e(p[k - 1].second)) * (p[k].first - p[k - 1].first);
    IA += 0.5 * (e(gs) + e(p[j - 1].second)) * (s_star - p[j - 1].first);
    IB += 0.5 * (e(p[j].second) + e(gs)) * (p[j].first - s_star);
    for (size_t k = j + 1; k < n; ++k)
        IB += 0.5 * (e(p[k].second) + e(p[k - 1].second)) * (p[k].first - p[k - 1].first);
    if (!(IA > 0.0) || !(IB > 0.0)) {
        r.warning = "empty basin integral";
        return r;
    }
    const double md_time_ps = CurcumaUnit::Constants::MD_TIME_UNIT_FS * 1e-3;
    r.s_star = s_star;
    r.g_star = gs;
    r.p_fwd = e(gs) / IA;
    r.p_rev = e(gs) / IB;
    r.sdot_abs = std::sqrt(2.0 * kT_eh / CurcumaUnit::Constants::PI) * zhalf / md_time_ps;
    r.k_fwd = 0.5 * r.sdot_abs * r.p_fwd;
    r.k_rev = 0.5 * r.sdot_abs * r.p_rev;
    r.eyring_prefactor = eyringPrefactorPerPs(kT_eh);
    r.dG_fwd = eyringFreeEnergy(r.k_fwd, kT_eh);
    r.dG_rev = eyringFreeEnergy(r.k_rev, kT_eh);
    r.valid = std::isfinite(r.k_fwd) && std::isfinite(r.k_rev) && r.k_fwd > 0.0 && r.k_rev > 0.0;
    if (!r.valid)
        r.warning = "non-finite rate";
    return r;
}

} // namespace fe_profile
} // namespace curcuma
