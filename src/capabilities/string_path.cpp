/*
 * <Pure path/string helpers for the finite-temperature string - implementation.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested - not TESTED/APPROVED).
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
 */

#include "src/capabilities/string_path.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace curcuma {
namespace string_path {

namespace {
/*! Weight of atom a (uniform 1 when w is empty). Claude Generated 2026. */
inline double weightOf(const std::vector<double>& w, int a)
{
    return w.empty() ? 1.0 : w[static_cast<size_t>(a)];
}

/*! Weighted squared norm sum_a w_a |d_a|^2. Claude Generated 2026. */
double weightedSquaredNorm(const Geometry& d, const std::vector<double>& w)
{
    if (w.empty())
        return d.squaredNorm();
    double s = 0.0;
    for (int a = 0; a < d.rows(); ++a)
        s += w[static_cast<size_t>(a)] * d.row(a).squaredNorm();
    return s;
}
} // namespace

Alignment kabsch(const Geometry& src, const Geometry& dst, const std::vector<double>& w)
{
    // Claude Generated 2026: weighted Kabsch (see header). Same SVD/reflection handling as
    // PathCV::rmsdTo, generalised to per-atom weights.
    Alignment A;
    const int n = static_cast<int>(src.rows());
    double wsum = 0.0;
    for (int a = 0; a < n; ++a) {
        const double wa = weightOf(w, a);
        wsum += wa;
        A.c_src += wa * Eigen::RowVector3d(src(a, 0), src(a, 1), src(a, 2));
        A.c_dst += wa * Eigen::RowVector3d(dst(a, 0), dst(a, 1), dst(a, 2));
    }
    if (wsum <= 0.0)
        return A;
    A.c_src /= wsum;
    A.c_dst /= wsum;

    Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
    for (int a = 0; a < n; ++a) {
        const double wa = weightOf(w, a);
        const Eigen::Vector3d p(src(a, 0) - A.c_src(0), src(a, 1) - A.c_src(1), src(a, 2) - A.c_src(2));
        const Eigen::Vector3d q(dst(a, 0) - A.c_dst(0), dst(a, 1) - A.c_dst(1), dst(a, 2) - A.c_dst(2));
        H += wa * p * q.transpose();
    }
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < 0)
        D(2, 2) = -1.0;
    A.rot = svd.matrixV() * D * svd.matrixU().transpose();
    return A;
}

Geometry applyAlignment(const Geometry& g, const Alignment& A)
{
    // Claude Generated 2026: x' = (x - c_src) rot^T + c_dst for every atom row.
    Geometry out(g.rows(), 3);
    const Eigen::Matrix3d RT = A.rot.transpose();
    for (int a = 0; a < g.rows(); ++a) {
        const Eigen::RowVector3d x(g(a, 0) - A.c_src(0), g(a, 1) - A.c_src(1), g(a, 2) - A.c_src(2));
        const Eigen::RowVector3d y = x * RT + A.c_dst;
        out(a, 0) = y(0);
        out(a, 1) = y(1);
        out(a, 2) = y(2);
    }
    return out;
}

double weightedRmsd(const Geometry& a, const Geometry& b, const std::vector<double>& w, bool best_fit)
{
    // Claude Generated 2026
    double wsum = 0.0;
    for (int i = 0; i < a.rows(); ++i)
        wsum += weightOf(w, i);
    if (wsum <= 0.0 || a.rows() != b.rows())
        return 0.0;
    const Geometry d = best_fit ? Geometry(applyAlignment(a, kabsch(a, b, w)) - b) : Geometry(a - b);
    return std::sqrt(weightedSquaredNorm(d, w) / wsum);
}

std::vector<Geometry> reparametrizeEqualArclength(const std::vector<Geometry>& path,
                                                  const std::vector<double>& w)
{
    // Claude Generated 2026: piecewise-linear equal-arclength re-sampling (string method).
    // The arithmetic for empty w is deliberately the one of the former
    // NebMD::reparametrizeString() (Frobenius norm, same segment search, same interpolation
    // formula), so that function can delegate here bit-identically (test_string_path (d)).
    const int n = static_cast<int>(path.size());
    std::vector<Geometry> out = path;
    if (n < 3)
        return out;

    std::vector<double> s(n, 0.0);
    for (int i = 1; i < n; ++i) {
        const Geometry d = path[i] - path[i - 1];
        s[i] = s[i - 1] + (w.empty() ? d.norm() : std::sqrt(weightedSquaredNorm(d, w)));
    }
    const double total = s[n - 1];
    if (total < 1e-9)
        return out;

    for (int i = 1; i < n - 1; ++i) {
        const double s_target = total * static_cast<double>(i) / static_cast<double>(n - 1);
        int seg = 0;
        while (seg < n - 2 && s[seg + 1] < s_target)
            ++seg;
        const double ds = s[seg + 1] - s[seg];
        const double t = (ds > 1e-12) ? (s_target - s[seg]) / ds : 0.0;
        out[i] = (1.0 - t) * path[seg] + t * path[seg + 1];
    }
    return out;
}

std::vector<Geometry> alignConsecutive(const std::vector<Geometry>& path, const std::vector<double>& w)
{
    // Claude Generated 2026
    std::vector<Geometry> out = path;
    for (size_t i = 1; i < out.size(); ++i)
        out[i] = applyAlignment(out[i], kabsch(out[i], out[i - 1], w));
    return out;
}

double maxFrameDisplacement(const std::vector<Geometry>& a, const std::vector<Geometry>& b,
                            const std::vector<double>& w, bool best_fit)
{
    // Claude Generated 2026
    double dmax = 0.0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i)
        dmax = std::max(dmax, weightedRmsd(a[i], b[i], w, best_fit));
    return dmax;
}

std::vector<Geometry> smoothPath(const std::vector<Geometry>& path, double kappa)
{
    // Claude Generated 2026: one explicit Laplacian pass on the interior frames.
    std::vector<Geometry> out = path;
    if (kappa <= 0.0 || path.size() < 3)
        return out;
    for (size_t i = 1; i + 1 < path.size(); ++i)
        out[i] = (1.0 - kappa) * path[i] + (0.5 * kappa) * (path[i - 1] + path[i + 1]);
    return out;
}

std::vector<Geometry> resamplePath(const std::vector<Geometry>& path, int nout,
                                   const std::vector<double>& w)
{
    // Claude Generated 2026: equal-arclength re-sampling onto nout frames (see header). Same
    // segment search and interpolation formula as reparametrizeEqualArclength.
    const int n = static_cast<int>(path.size());
    if (n < 2 || nout < 2)
        return path;
    std::vector<double> s(n, 0.0);
    for (int i = 1; i < n; ++i) {
        const Geometry d = path[i] - path[i - 1];
        s[i] = s[i - 1] + (w.empty() ? d.norm() : std::sqrt(weightedSquaredNorm(d, w)));
    }
    const double total = s[n - 1];
    if (total < 1e-9 && nout == n)
        return path; // degenerate path: same convention as reparametrizeEqualArclength
    std::vector<Geometry> out(nout);
    out.front() = path.front();
    out.back() = path.back();
    for (int i = 1; i < nout - 1; ++i) {
        if (total < 1e-9) {
            out[i] = path.front();
            continue;
        }
        const double s_target = total * static_cast<double>(i) / static_cast<double>(nout - 1);
        int seg = 0;
        while (seg < n - 2 && s[seg + 1] < s_target)
            ++seg;
        const double ds = s[seg + 1] - s[seg];
        const double t = (ds > 1e-12) ? (s_target - s[seg]) / ds : 0.0;
        out[i] = (1.0 - t) * path[seg] + t * path[seg + 1];
    }
    return out;
}

std::vector<Geometry> extendPath(const std::vector<Geometry>& path, int m, const std::vector<double>& w)
{
    // Claude Generated 2026: straight continuation of the two end segments (see header).
    const int n = static_cast<int>(path.size());
    if (m <= 0 || n < 2)
        return path;
    const Geometry& R0 = path.front();
    const Geometry& RN = path.back();
    const Geometry R1a = applyAlignment(path[1], kabsch(path[1], R0, w));
    const Geometry RN2a = applyAlignment(path[n - 2], kabsch(path[n - 2], RN, w));
    const Geometry d0 = R0 - R1a;  // direction away from the path at the start
    const Geometry dN = RN - RN2a; // ... at the end
    std::vector<Geometry> out;
    out.reserve(static_cast<size_t>(n + 2 * m));
    for (int j = m; j >= 1; --j)
        out.push_back(R0 + static_cast<double>(j) * d0);
    for (const auto& g : path)
        out.push_back(g);
    for (int j = 1; j <= m; ++j)
        out.push_back(RN + static_cast<double>(j) * dN);
    return out;
}

int extensionFrames(double length_fraction, int nframes)
{
    // Claude Generated 2026: m = ceil(f L / d) + 1 with d = L/(F-1)  ->  ceil(f (F-1)) + 1.
    // The +1 is the outermost extension frame, which remains a plateau of s without a window.
    if (!(length_fraction > 0.0) || nframes < 2)
        return 0;
    const double steps = length_fraction * static_cast<double>(nframes - 1); // f L / d
    return static_cast<int>(std::ceil(steps - 1e-9)) + 1;
}

StepLimit limitStep(const std::vector<Geometry>& old_path, const std::vector<Geometry>& proposed,
                    double max_step, const std::vector<double>& w, bool fix_endpoints)
{
    // Claude Generated 2026: trust radius on the whole mean-path update (see header).
    StepLimit res;
    res.path = proposed;
    const int n = static_cast<int>(std::min(old_path.size(), proposed.size()));
    if (n == 0)
        return res;
    auto fixed = [&](int i) { return fix_endpoints && (i == 0 || i == n - 1); };
    for (int i = 0; i < n; ++i)
        if (!fixed(i))
            res.delta_max = std::max(res.delta_max, weightedRmsd(proposed[i], old_path[i], w, true));
    if (!(max_step > 0.0) || res.delta_max <= max_step)
        return res; // nothing to limit: proposed path bit-identical
    res.scale = max_step / res.delta_max;
    for (int i = 0; i < n; ++i) {
        if (fixed(i)) {
            res.path[i] = old_path[i];
            continue;
        }
        const Geometry P = applyAlignment(proposed[i], kabsch(proposed[i], old_path[i], w));
        res.path[i] = old_path[i] + res.scale * (P - old_path[i]);
    }
    return res;
}

MeanPathUpdate::~MeanPathUpdate() = default;

MeanPathUpdate updateMeanPath(const std::vector<Geometry>& old_path, const std::vector<Geometry>& means,
                              const std::vector<std::vector<Geometry>>& block_means, double k,
                              const std::vector<double>& w, double alpha, double gain, double gain_z,
                              double ftol, double max_step, double smooth, bool fix_endpoints)
{
    // Claude Generated 2026: see header. Order: mean -> tau/perp -> block statistics/F/SE/z ->
    // gain -> mixing -> smoothing -> reparametrisation -> delta/delta_perp -> limitStep.
    MeanPathUpdate res;
    const int N = static_cast<int>(old_path.size());
    res.path = old_path;
    res.perp_rmsd.assign(N, 0.0);
    res.f_rms.assign(N, 0.0);
    res.f_se.assign(N, 0.0);
    res.f_sig.assign(N, 0.0);
    res.z.assign(N, 0.0);
    res.gain.assign(N, 1.0);
    if (N == 0)
        return res;
    const int nat = static_cast<int>(old_path[0].rows());
    double wsum = 0.0;
    for (int a = 0; a < nat; ++a)
        wsum += weightOf(w, a);
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<Geometry>& Mnew = res.path;
    std::vector<Geometry> tau_unit(N); // weighted unit tangent of the OLD path (rows 0 if undefined)
    std::vector<char> updated(N, 0);

    // perpendicular projection of a deviation d with the unit tangent t (weighted dot product)
    auto perpOf = [&](const Geometry& d, const Geometry& t) {
        if (t.rows() == 0)
            return d;
        double par = 0.0;
        for (int a = 0; a < d.rows(); ++a)
            par += weightOf(w, a) * d.row(a).dot(t.row(a));
        return Geometry(d - par * t);
    };

    for (int i = 0; i < N; ++i) {
        const bool endpoint = (i == 0 || i == N - 1);
        if (i < static_cast<int>(means.size()) && means[i].rows() == nat && !(endpoint && fix_endpoints)) {
            updated[i] = 1;
        }
        Geometry tau = old_path[std::min(i + 1, N - 1)] - old_path[std::max(i - 1, 0)];
        const double tn = std::sqrt(weightedSquaredNorm(tau, w));
        if (tn > 1e-12)
            tau_unit[i] = tau / tn;
        if (!updated[i])
            continue;
        ++res.n_updated;
        const Geometry& mean = means[i];
        const Geometry dev = mean - old_path[i];
        const Geometry perp = perpOf(dev, tau_unit[i]);
        res.perp_rmsd[i] = std::sqrt(weightedSquaredNorm(perp, w) / std::max(wsum, 1e-300));
        res.perp_max = std::max(res.perp_max, res.perp_rmsd[i]);

        // mean force F_a = k w_a perp_a, RMS per atom
        double f2 = 0.0;
        for (int a = 0; a < nat; ++a) {
            const double ka = k * weightOf(w, a);
            f2 += ka * ka * perp.row(a).squaredNorm();
        }
        res.f_rms[i] = std::sqrt(f2 / std::max(nat, 1));

        // block statistics: sample variance of every force component over the blocks
        const std::vector<Geometry>* blocks = (i < static_cast<int>(block_means.size())) ? &block_means[i] : nullptr;
        int nb = 0;
        if (blocks)
            for (const auto& b : *blocks)
                if (b.rows() == nat)
                    ++nb;
        if (nb >= 2) {
            Geometry fsum = Geometry::Zero(nat, 3), f2sum = Geometry::Zero(nat, 3);
            for (const auto& b : *blocks) {
                if (b.rows() != nat)
                    continue;
                const Geometry pb = perpOf(b - old_path[i], tau_unit[i]);
                for (int a = 0; a < nat; ++a) {
                    const double ka = k * weightOf(w, a);
                    for (int c = 0; c < 3; ++c) { // scalar: avoids a GCC -Warray-bounds false positive on RowVector3d
                        const double fb = ka * pb(a, c);
                        fsum(a, c) += fb;
                        f2sum(a, c) += fb * fb;
                    }
                }
            }
            double s2tot = 0.0;
            for (int a = 0; a < nat; ++a)
                for (int c = 0; c < 3; ++c) {
                    const double m = fsum(a, c) / nb;
                    s2tot += std::max(f2sum(a, c) - nb * m * m, 0.0) / (nb - 1);
                }
            res.f_se[i] = std::sqrt(s2tot / nb / std::max(nat, 1));
        } else {
            res.f_se[i] = inf;
        }
        res.f_sig[i] = std::sqrt(std::max(res.f_rms[i] * res.f_rms[i] - (std::isfinite(res.f_se[i]) ? res.f_se[i] * res.f_se[i] : 0.0), 0.0));
        res.z[i] = (std::isfinite(res.f_se[i]) && res.f_se[i] > 0.0) ? res.f_sig[i] / res.f_se[i] : 0.0;

        // gain on the perpendicular part (only when the force is significant)
        double G = 1.0;
        // Gate: significant (z >= gain_z) AND above the convergence level (f_sig >= ftol > 0).
        // z alone is t-distributed with n_b - 1 dof (heavy tails): with 5 blocks z >= 3 fired in
        // 3.7 % of zero-signal frame-updates and each false trigger moved the frame by max_step.
        if (gain > 1.0 && ftol > 0.0 && std::isfinite(res.f_se[i]) && res.z[i] >= gain_z && res.f_sig[i] >= ftol
            && res.perp_rmsd[i] > 0.0) {
            G = gain;
            if (max_step > 0.0 && alpha > 0.0)
                G = std::min(G, max_step / (alpha * res.perp_rmsd[i]));
            G = std::max(G, 1.0);
        }
        res.gain[i] = G;
        if (G == 1.0)
            Mnew[i] = (1.0 - alpha) * old_path[i] + alpha * mean; // exact former expression
        else
            Mnew[i] = old_path[i] + alpha * dev + (alpha * (G - 1.0)) * perp;

        res.f_rms_max = std::max(res.f_rms_max, res.f_rms[i]);
        res.f_se_max = std::max(res.f_se_max, res.f_se[i]);
        res.f_sig_max = std::max(res.f_sig_max, res.f_sig[i]);
        res.z_max = std::max(res.z_max, res.z[i]);
        res.gain_max = std::max(res.gain_max, G);
    }
    if (smooth > 0.0)
        Mnew = smoothPath(Mnew, smooth);
    Mnew = reparametrizeEqualArclength(Mnew, w);
    if (fix_endpoints) { // exact, independent of the helper's guarantee
        Mnew.front() = old_path.front();
        Mnew.back() = old_path.back();
    }

    // delta (best fit, all frames as in the former code) and its perpendicular part
    double d2 = 0.0;
    int nd = 0;
    for (int i = 0; i < N; ++i) {
        const double r = weightedRmsd(Mnew[i], old_path[i], w, true);
        res.delta_max = std::max(res.delta_max, r);
        if (i > 0 && i < N - 1) {
            d2 += r * r;
            ++nd;
        }
        if (fix_endpoints && (i == 0 || i == N - 1))
            continue;
        const Geometry P = applyAlignment(Mnew[i], kabsch(Mnew[i], old_path[i], w));
        const Geometry sp = perpOf(P - old_path[i], tau_unit[i]);
        res.delta_perp_max = std::max(res.delta_perp_max, std::sqrt(weightedSquaredNorm(sp, w) / std::max(wsum, 1e-300)));
    }
    res.delta_rms = std::sqrt(d2 / std::max(nd, 1));

    if (max_step > 0.0 && res.delta_max > max_step) {
        StepLimit lim = limitStep(old_path, Mnew, max_step, w, fix_endpoints);
        Mnew = std::move(lim.path);
        res.scale = lim.scale;
    }
    for (int i = 1; i < N; ++i)
        res.arclength += weightedRmsd(Mnew[i], Mnew[i - 1], w, false);
    return res;
}

} // namespace string_path
} // namespace curcuma
