/*
 * <Path collective variables (Branduardi s,z) - implementation.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested pending - not TESTED/APPROVED).
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

#include "src/capabilities/path_cv.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "src/core/curcuma_logger.h"
#include "src/core/parameter_registry.h"

PathCV::PathCV(const std::vector<Geometry>& path, const Molecule& templ, double lambda,
               const std::vector<int>& atoms, Metric metric, const std::vector<AtomPair>& pairs,
               PairTransform transform)
    : m_path(path)
    , m_atoms(atoms)
    , m_template(templ)
    , m_metric(metric)
    , m_pairs(pairs)
    , m_transform(transform)
{
    if (m_path.size() < 2)
        CurcumaLogger::error("PathCV: need at least two reference frames.");
    if (m_metric == Metric::Distances) {
        // Claude Generated 2026: distance metric - pair coordinates of every frame, built once.
        if (m_pairs.empty())
            throw std::invalid_argument("PathCV: distance metric without atom pairs");
        m_path_x.reserve(m_path.size());
        for (const auto& g : m_path)
            m_path_x.push_back(pairValues(g, m_pairs, m_transform));
        m_lambda = (lambda > 0.0) ? lambda : computeLambdaFromPairs(m_path, m_pairs, m_transform);
        return;
    }
    // Claude Generated 2026: the reference frames restricted to the selected atoms, built once.
    m_path_sel.reserve(m_path.size());
    for (const auto& g : m_path)
        m_path_sel.push_back(selectRows(g));
    m_lambda = (lambda > 0.0) ? lambda : computeLambdaFromPath(m_path, m_atoms);
}

Vector PathCV::pairValues(const Geometry& R, const std::vector<AtomPair>& pairs, PairTransform t)
{
    // Claude Generated 2026: x_p = |R_a - R_b| or 1/|R_a - R_b|.
    Vector x(static_cast<int>(pairs.size()));
    for (size_t p = 0; p < pairs.size(); ++p) {
        const double d = (R.row(pairs[p].first) - R.row(pairs[p].second)).norm();
        x(static_cast<int>(p)) = (t == PairTransform::Inverse) ? 1.0 / d : d;
    }
    return x;
}

double PathCV::pairSpacing2(const Geometry& A, const Geometry& B, const std::vector<AtomPair>& pairs,
                            PairTransform t)
{
    // Claude Generated 2026
    if (pairs.empty())
        return 0.0;
    return (pairValues(A, pairs, t) - pairValues(B, pairs, t)).squaredNorm() / static_cast<double>(pairs.size());
}

double PathCV::computeLambdaFromPairs(const std::vector<Geometry>& path, const std::vector<AtomPair>& pairs,
                                      PairTransform t)
{
    // Claude Generated 2026: same rule as computeLambdaFromPath, lambda = 1/<d_{k,k+1}^2>, with d
    // the distance-space spacing of consecutive frames (the metric of the CV itself).
    if (path.size() < 2)
        return 1.0;
    double sum = 0.0;
    for (size_t i = 1; i < path.size(); ++i)
        sum += pairSpacing2(path[i], path[i - 1], pairs, t);
    const double mean_d2 = sum / static_cast<double>(path.size() - 1);
    return (mean_d2 > 1e-12) ? 1.0 / mean_d2 : 1.0;
}

bool PathCV::parsePairList(const std::string& spec, int natoms, std::vector<AtomPair>& pairs, std::string& error)
{
    // Claude Generated 2026: "1-4,1-9,4-9" -> {(0,3),(0,8),(3,8)}.
    pairs.clear();
    error.clear();
    auto trim = [](std::string t) {
        const auto b = t.find_first_not_of(" \t");
        const auto e = t.find_last_not_of(" \t");
        return b == std::string::npos ? std::string() : t.substr(b, e - b + 1);
    };
    const std::string s = trim(spec);
    size_t pos = 0;
    while (pos <= s.size()) {
        const size_t comma = s.find(',', pos);
        const std::string tok = trim(s.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos));
        pos = (comma == std::string::npos) ? s.size() + 1 : comma + 1;
        const size_t dash = tok.find('-');
        int i = 0, j = 0;
        bool ok = (dash != std::string::npos) && tok.find('-', dash + 1) == std::string::npos;
        if (ok) {
            const std::string a = trim(tok.substr(0, dash)), b = trim(tok.substr(dash + 1));
            try {
                size_t ua = 0, ub = 0;
                i = std::stoi(a, &ua);
                j = std::stoi(b, &ub);
                ok = (ua == a.size()) && (ub == b.size()) && !a.empty() && !b.empty();
            } catch (const std::exception&) {
                ok = false;
            }
        }
        if (!ok) {
            error = "cv_pairs: cannot read '" + tok + "' (expected 'auto' or 1-based pairs i-j, e.g. 1-4,1-9,4-9)";
            pairs.clear();
            return false;
        }
        if (i == j) {
            error = "cv_pairs: pair '" + tok + "' joins an atom with itself";
            pairs.clear();
            return false;
        }
        if (i < 1 || j < 1 || i > natoms || j > natoms) {
            error = "cv_pairs: index in '" + tok + "' outside 1.." + std::to_string(natoms);
            pairs.clear();
            return false;
        }
        const AtomPair p(std::min(i, j) - 1, std::max(i, j) - 1);
        if (std::find(pairs.begin(), pairs.end(), p) != pairs.end()) {
            error = "cv_pairs: pair " + std::to_string(p.first + 1) + "-" + std::to_string(p.second + 1)
                + " is given twice";
            pairs.clear();
            return false;
        }
        pairs.push_back(p);
    }
    if (pairs.size() < 2) {
        error = "cv_pairs '" + s + "' gives " + std::to_string(pairs.size())
            + " pair(s); the distance-space path CV needs at least 2";
        pairs.clear();
        return false;
    }
    return true;
}

std::vector<PathCV::AtomPair> PathCV::selectPairsBySpan(const std::vector<Geometry>& frames,
                                                        const std::vector<int>& candidates, double threshold,
                                                        std::vector<double>* spans)
{
    // Claude Generated 2026: pairs whose plain distance changes by more than `threshold` along
    // the frames (span = max - min). Order: lexicographic in (first, second).
    std::vector<AtomPair> out;
    if (spans)
        spans->clear();
    std::vector<int> c = candidates;
    std::sort(c.begin(), c.end());
    for (size_t u = 0; u < c.size(); ++u)
        for (size_t v = u + 1; v < c.size(); ++v) {
            double dmin = std::numeric_limits<double>::max(), dmax = -dmin;
            for (const auto& g : frames) {
                const double d = (g.row(c[u]) - g.row(c[v])).norm();
                dmin = std::min(dmin, d);
                dmax = std::max(dmax, d);
            }
            if (!frames.empty() && dmax - dmin > threshold) {
                out.emplace_back(c[u], c[v]);
                if (spans)
                    spans->push_back(dmax - dmin);
            }
        }
    return out;
}

std::vector<Geometry> PathCV::resampleByPairs(const std::vector<Geometry>& path, int nout,
                                              const std::vector<AtomPair>& pairs, PairTransform t, int sub)
{
    // Claude Generated 2026: equal arclength in the distance metric. Same structure as
    // string_path::resamplePath (cumulative arclength, linear Cartesian interpolation), but the
    // arclength of segment k is the distance-space length of the straight Cartesian segment,
    // integrated over `sub` sub-steps: x_p(R) is not linear along R_k + t (R_k+1 - R_k), so a
    // single chord would misplace the frames inside strongly curved segments.
    const int n = static_cast<int>(path.size());
    if (n < 2 || nout < 2 || pairs.empty())
        return path;
    sub = std::max(sub, 1);
    const double np = static_cast<double>(pairs.size());
    // fine grid: cumulative length at every sub-step of every segment
    std::vector<double> sgrid;      // cumulative length at grid point
    std::vector<double> tgrid;      // global parameter k + t of grid point
    sgrid.push_back(0.0);
    tgrid.push_back(0.0);
    Vector xprev = pairValues(path[0], pairs, t);
    for (int k = 0; k + 1 < n; ++k)
        for (int q = 1; q <= sub; ++q) {
            const double tt = static_cast<double>(q) / sub;
            const Geometry Rq = (1.0 - tt) * path[k] + tt * path[k + 1];
            const Vector xq = pairValues(Rq, pairs, t);
            sgrid.push_back(sgrid.back() + std::sqrt((xq - xprev).squaredNorm() / np));
            tgrid.push_back(k + tt);
            xprev = xq;
        }
    const double total = sgrid.back();
    if (total < 1e-12 && nout == n)
        return path;
    std::vector<Geometry> out(nout);
    out.front() = path.front();
    out.back() = path.back();
    const int ng = static_cast<int>(sgrid.size());
    for (int i = 1; i < nout - 1; ++i) {
        if (total < 1e-12) {
            out[i] = path.front();
            continue;
        }
        const double target = total * static_cast<double>(i) / static_cast<double>(nout - 1);
        int g = 0;
        while (g < ng - 2 && sgrid[g + 1] < target)
            ++g;
        const double dsg = sgrid[g + 1] - sgrid[g];
        const double f = (dsg > 1e-15) ? (target - sgrid[g]) / dsg : 0.0;
        const double T = tgrid[g] + f * (tgrid[g + 1] - tgrid[g]);
        int k = std::min(static_cast<int>(std::floor(T)), n - 2);
        const double tt = T - k;
        out[i] = (1.0 - tt) * path[k] + tt * path[k + 1];
    }
    return out;
}

Geometry PathCV::selectRows(const Geometry& g) const
{
    // Claude Generated 2026
    if (m_atoms.empty())
        return g;
    Geometry out(static_cast<int>(m_atoms.size()), 3);
    for (size_t k = 0; k < m_atoms.size(); ++k)
        out.row(static_cast<int>(k)) = g.row(m_atoms[k]);
    return out;
}

double PathCV::computeLambdaFromPath(const std::vector<Geometry>& path, const std::vector<int>& atoms)
{
    // lambda ~ 1/<d^2> with d the RMSD spacing of neighbouring frames, so that the
    // Gaussians of adjacent frames overlap appreciably. Using the plain (non-best-fit)
    // spacing here is sufficient for the scale estimate. With an atom selection the mean
    // square runs over the selected atoms only (the metric of the CV itself).
    if (path.size() < 2)
        return 1.0;
    const int nat = static_cast<int>(path[0].rows());
    double sum = 0.0;
    int cnt = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        double d2 = 0.0;
        if (atoms.empty()) {
            d2 = (path[i] - path[i - 1]).squaredNorm() / std::max(nat, 1);
        } else {
            for (int a : atoms)
                d2 += (path[i].row(a) - path[i - 1].row(a)).squaredNorm();
            d2 /= static_cast<double>(atoms.size());
        }
        sum += d2;
        ++cnt;
    }
    const double mean_d2 = (cnt > 0) ? sum / cnt : 1.0;
    return (mean_d2 > 1e-12) ? 1.0 / mean_d2 : 1.0;
}

bool PathCV::parseAtomSelection(const std::string& spec, const Molecule& mol, std::vector<int>& atoms,
                                std::string& error)
{
    // Claude Generated 2026: "all" | "heavy" | 1-based list with ranges ("1-6,9").
    atoms.clear();
    error.clear();
    const int nat = static_cast<int>(mol.AtomCount());
    auto trim = [](std::string t) {
        const auto b = t.find_first_not_of(" \t");
        const auto e = t.find_last_not_of(" \t");
        return b == std::string::npos ? std::string() : t.substr(b, e - b + 1);
    };
    const std::string s = trim(spec);
    if (s.empty() || s == "all")
        return true;
    if (s == "heavy") {
        for (int a = 0; a < nat; ++a)
            if (mol.Atom(a).first > 1)
                atoms.push_back(a);
    } else {
        std::vector<char> used(static_cast<size_t>(std::max(nat, 0)), 0);
        size_t pos = 0;
        while (pos <= s.size()) {
            const size_t comma = s.find(',', pos);
            const std::string tok = trim(s.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos));
            pos = (comma == std::string::npos) ? s.size() + 1 : comma + 1;
            int lo = 0, hi = 0;
            const size_t dash = tok.find('-');
            try {
                size_t used_chars = 0;
                if (dash == std::string::npos) {
                    lo = hi = std::stoi(tok, &used_chars);
                    if (used_chars != tok.size())
                        throw std::invalid_argument(tok);
                } else {
                    const std::string a = trim(tok.substr(0, dash)), b = trim(tok.substr(dash + 1));
                    lo = std::stoi(a, &used_chars);
                    if (used_chars != a.size())
                        throw std::invalid_argument(tok);
                    hi = std::stoi(b, &used_chars);
                    if (used_chars != b.size())
                        throw std::invalid_argument(tok);
                }
            } catch (const std::exception&) {
                error = "cv_atoms: cannot read '" + tok + "' (expected all, heavy or a 1-based list like 1-6,9)";
                return false;
            }
            if (lo > hi) {
                error = "cv_atoms: range '" + tok + "' is descending";
                return false;
            }
            if (lo < 1 || hi > nat) {
                error = "cv_atoms: index in '" + tok + "' outside 1.." + std::to_string(nat);
                return false;
            }
            for (int k = lo; k <= hi; ++k) {
                if (used[static_cast<size_t>(k - 1)]) {
                    error = "cv_atoms: atom " + std::to_string(k) + " is selected twice";
                    return false;
                }
                used[static_cast<size_t>(k - 1)] = 1;
                atoms.push_back(k - 1);
            }
        }
        std::sort(atoms.begin(), atoms.end());
    }
    if (atoms.size() < 3) {
        error = "cv_atoms '" + s + "' selects " + std::to_string(atoms.size())
            + " atom(s); the path CV needs at least 3 (fewer cannot fix the rotation about their axis)";
        atoms.clear();
        return false;
    }
    return true;
}

std::vector<double> PathCV::selectionWeights(const std::vector<int>& atoms, int natoms, double other)
{
    // Claude Generated 2026
    if (atoms.empty())
        return {};
    std::vector<double> w(static_cast<size_t>(natoms), other);
    for (int a : atoms)
        w[static_cast<size_t>(a)] = 1.0;
    return w;
}

double PathCV::rmsdTo(int i, const Geometry& R, Geometry* grad) const
{
    // Claude Generated 2026: minimal RMSD of R onto frame i (Kabsch, exact) and its gradient
    // in the LABORATORY frame of R.
    //
    // Replaces the RMSDDriver call. That driver aligns the target onto the reference and
    // returns the gradient in the aligned (reference) frame, so the gradient was only right
    // for geometries already oriented like the reference. Measured: relative error of ds/dR
    // 4 % unrotated, 49 % after a 30 degree rotation, 99 % after 90 degrees (a
    // finite-difference check on a non-collinear molecule, see test_path_cv). Since NEB-MD
    // lets the images rotate freely (no COM/rotation removal), that error grew with time.
    //
    //   d^2 = 1/N sum_a |Rot (p_a - <p>) - (q_a - <q>)|^2 ,  Rot = V D U^T from SVD(P^T Q)
    //   d * dd/dp_a = Rot^T (Rot (p_a - <p>) - (q_a - <q>)) / N     (envelope theorem)
    //
    // The returned gradient keeps the sign convention of the old GradientTimesRMSD():
    // it is MINUS d*dd/dR, and the caller folds that sign into its chain rule.
    //
    // Atom selection A (Claude Generated 2026): p_a, q_a, the centroids and N run over A only
    // (N = |A|); the gradient rows of atoms outside A are exactly zero. Empty A = all atoms.
    const bool all = m_atoms.empty();
    Geometry Rsel;
    if (!all)
        Rsel = selectRows(R);
    const Geometry& Rs = all ? R : Rsel;
    const Geometry& Q = all ? m_path[i] : m_path_sel[i];
    const int N = static_cast<int>(Rs.rows());
    const Eigen::RowVector3d cp = Rs.colwise().mean();
    const Eigen::RowVector3d cq = Q.colwise().mean();
    const Eigen::MatrixXd P = Rs.rowwise() - cp;
    const Eigen::MatrixXd Qc = Q.rowwise() - cq;
    const Eigen::Matrix3d H = P.transpose() * Qc;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < 0)
        D(2, 2) = -1.0;
    const Eigen::Matrix3d Rot = svd.matrixV() * D * svd.matrixU().transpose();
    // rows: aligned_a = Rot * p_a  ->  P * Rot^T
    const Eigen::MatrixXd res = P * Rot.transpose() - Qc;
    const double d = std::sqrt(res.squaredNorm() / std::max(N, 1));
    if (grad) {
        if (all) {
            *grad = -(res * Rot) / static_cast<double>(std::max(N, 1)); // rows: Rot^T res_a
        } else {
            const Geometry gsel = -(res * Rot) / static_cast<double>(std::max(N, 1));
            *grad = Geometry::Zero(R.rows(), 3);
            for (size_t k = 0; k < m_atoms.size(); ++k)
                grad->row(m_atoms[k]) = gsel.row(static_cast<int>(k));
        }
    }
    return d;
}

double PathCV::distTo(int i, const Geometry& R, const Vector& x, Geometry* grad) const
{
    // Claude Generated 2026: distance-space d_i and its gradient.
    //
    //   d_i^2 = 1/P sum_p (x_p - x_p,i)^2
    //   d_i dd_i/dR = 1/P sum_p (x_p - x_p,i) dx_p/dR
    //   dx_p/dR_a = +u/r, dx_p/dR_b = -u/r       (x = r = |u|, u = R_a - R_b)
    //   dx_p/dR_a = -u/r^3, dx_p/dR_b = +u/r^3   (x = 1/r)
    //
    // Returned is MINUS d_i dd_i/dR, the sign convention of rmsdTo(), so evaluate() and
    // everything built on it (z wall, SHAKE/RATTLE, mean force) are metric-agnostic.
    const Vector& xi = m_path_x[static_cast<size_t>(i)];
    const double P = static_cast<double>(m_pairs.size());
    const Vector diff = x - xi;
    const double d = std::sqrt(diff.squaredNorm() / P);
    if (grad) {
        *grad = Geometry::Zero(R.rows(), 3);
        for (size_t p = 0; p < m_pairs.size(); ++p) {
            const int a = m_pairs[p].first, b = m_pairs[p].second;
            const Eigen::RowVector3d u = R.row(a) - R.row(b);
            const double r = u.norm();
            const Eigen::RowVector3d dxda = (m_transform == PairTransform::Inverse) ? Eigen::RowVector3d(-u / (r * r * r))
                                                                                     : Eigen::RowVector3d(u / r);
            const double c = -diff(static_cast<int>(p)) / P;
            grad->row(a) += c * dxda;
            grad->row(b) -= c * dxda;
        }
    }
    return d;
}

void PathCV::allDistances(const Geometry& R, std::vector<double>& d, std::vector<Geometry>* gd) const
{
    // Claude Generated 2026: d_i (and gradients) of every frame in the active metric.
    const int N = static_cast<int>(m_path.size());
    d.assign(N, 0.0);
    if (gd)
        gd->assign(N, Geometry());
    if (m_metric == Metric::Distances) {
        const Vector x = pairValues(R, m_pairs, m_transform);
        for (int i = 0; i < N; ++i)
            d[i] = distTo(i, R, x, gd ? &(*gd)[i] : nullptr);
        return;
    }
    for (int i = 0; i < N; ++i)
        d[i] = rmsdTo(i, R, gd ? &(*gd)[i] : nullptr);
}

int PathCV::nearestFrame(const Geometry& R, double* z_out) const
{
    // Claude Generated 2026: argmin_i d_i plus z = -1/lambda ln sum_i exp(-lambda d_i^2) from
    // the same d_i (same max-shift as evaluate(), so z agrees with evaluate() bit for bit).
    std::vector<double> d;
    allDistances(R, d, nullptr);
    const int N = static_cast<int>(d.size());
    int best = 0;
    double emax = -std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i) {
        if (d[i] < d[best])
            best = i;
        emax = std::max(emax, -m_lambda * d[i] * d[i]);
    }
    if (z_out) {
        double wsum = 0.0;
        for (int i = 0; i < N; ++i)
            wsum += std::exp(-m_lambda * d[i] * d[i] - emax);
        if (wsum < 1e-300)
            wsum = 1e-300;
        *z_out = -(1.0 / m_lambda) * (std::log(wsum) + emax);
    }
    return best;
}

void PathCV::evaluate(const Geometry& R, double* s_out, double* z_out, Geometry* ds_dR,
                      Geometry* dz_dR) const
{
    const int N = static_cast<int>(m_path.size());
    const int nat = static_cast<int>(R.rows());
    if (N < 2)
        return;

    std::vector<double> d(N, 0.0);
    std::vector<Geometry> gd(N);
    const bool need_grad = (ds_dR != nullptr) || (dz_dR != nullptr);
    // Claude Generated 2026: d_i in the active metric (rmsd: the former loop over rmsdTo).
    allDistances(R, d, need_grad ? &gd : nullptr);
    if (!need_grad)
        gd.assign(N, Geometry());

    // Softmax weights with the max subtracted for numerical stability.
    double emax = -std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i)
        emax = std::max(emax, -m_lambda * d[i] * d[i]);

    std::vector<double> w(N, 0.0);
    double wsum = 0.0;
    for (int i = 0; i < N; ++i) {
        w[i] = std::exp(-m_lambda * d[i] * d[i] - emax);
        wsum += w[i];
    }
    if (wsum < 1e-300)
        wsum = 1e-300;

    double s = 0.0;
    for (int i = 0; i < N; ++i)
        s += static_cast<double>(i) * w[i];
    s /= (wsum * (N - 1));

    if (s_out)
        *s_out = s;
    if (z_out)
        *z_out = -(1.0 / m_lambda) * (std::log(wsum) + emax);

    // dz/dR = sum_i (w_i/wsum) * 2 * d_i * dd_i/dR. The driver gradient is
    // -d(RMSD)/dR (see the sign note below), hence the leading minus here.
    if (dz_dR) {
        *dz_dR = Geometry::Zero(nat, 3);
        for (int i = 0; i < N; ++i) {
            if (gd[i].rows() != nat || gd[i].cols() != 3)
                continue;
            *dz_dR += (-(w[i] / wsum) * 2.0) * gd[i];
        }
    }

    if (!ds_dR)
        return;

    // ds/dR = sum_i [ i/(N-1) - s ] * (w_i/wsum) * (-2 lambda d_i) * dd_i/dR
    //
    // NOTE on the sign: RMSDDriver::Gradient() returns (ref - tar)/(RMSD*N), i.e.
    // MINUS d(RMSD)/dR_target (RMSD itself is sqrt(sum|delta|^2/N), see
    // RMSDFunctions::getRMSD). Folding that minus into the prefactor turns the
    // -2*lambda*d of the chain rule into +2*lambda*d. Getting this wrong shows up as
    // an exact factor -2 against finite differences, which is what the unit test
    // caught.
    *ds_dR = Geometry::Zero(nat, 3);
    for (int i = 0; i < N; ++i) {
        // The factor d_i that the chain rule produces is already contained in gd[i]
        // (= d_i * dd_i/dR, see rmsdTo), so it must NOT appear here a second time.
        // Keeping the two together is what removes the 1/RMSD singularity: a structure
        // sitting exactly on a reference frame is now an ordinary point.
        const double pref = (static_cast<double>(i) / (N - 1) - s) * (w[i] / wsum)
            * (2.0 * m_lambda);
        if (gd[i].rows() == nat && gd[i].cols() == 3)
            *ds_dR += pref * gd[i];
    }
}

double PathCV::zWallEnergy(const Geometry& R, double z_wall, double k, Geometry* grad,
                           double* z_out) const
{
    // Claude Generated 2026: flat-bottom wall on the distance-from-path z ("tube").
    //
    // Why: with s held by a constraint, the surface s = s0 is a (3N-7)-dimensional set that can
    // contain regions far from the reference path (retro-Diels-Alder: the dissociated ethylene
    // region at z ~ 1.4-3.8 A^2 next to the TS). Nothing in the window dynamics restrains z, so
    // the windows drift there. W confines them to z <= z_wall without touching the tube interior.
    //
    //   W(z)  = 0                          for z <= z_wall
    //   W(z)  = k/2 (z - z_wall)^2         for z >  z_wall
    //   dW/dR = W'(z) dz/dR = k (z - z_wall) dz/dR
    //
    // dz/dR comes from evaluate() in the laboratory frame (rmsdTo() returns lab-frame RMSD
    // gradients), so dW/dR can be added to dE/dR directly. W and dW/dR are continuous at
    // z = z_wall; the second derivative jumps there (harmless for velocity Verlet at the
    // stiffness chosen by zWallDefaultK).
    //
    // Note on the free energy: the sampled ensemble is exp(-(E+W)/kT), so the PMF obtained with
    // the wall is the PMF of the z-restricted system. It depends on z_wall.
    const int nat = static_cast<int>(R.rows());
    double s = 0.0, z = 0.0;
    Geometry ds, dz;
    evaluate(R, &s, &z, grad ? &ds : nullptr, grad ? &dz : nullptr);
    if (z_out)
        *z_out = z;
    const double excess = z - z_wall;
    if (excess <= 0.0 || k <= 0.0) {
        if (grad)
            *grad = Geometry::Zero(nat, 3);
        return 0.0;
    }
    if (grad)
        *grad = (k * excess) * dz;
    return 0.5 * k * excess * excess;
}

double PathCV::zWallDefaultK(double kT, double z_wall)
{
    // Claude Generated 2026: stiffness so that a thermal excursion is stopped close to the wall.
    // Require W = 3 kT at a penetration delta = 0.2 z_wall beyond the wall:
    //     k/2 delta^2 = 3 kT   ->   k = 6 kT / (0.2 z_wall)^2 = 150 kT / z_wall^2
    // The Boltzmann factor at that depth is exp(-3) = 0.05, so few samples reach z > 1.2 z_wall.
    // Example 300 K (kT = 9.50e-4 Eh), z_wall = 0.5 A^2: k = 0.57 Eh/A^4. Force at that depth:
    // |dW/dR| = k delta |dz/dR| with |dz/dR| ~ 2 d / sqrt(N) (whole 3N vector, d = RMSD ~
    // sqrt(z)); 16 atoms: 0.57 * 0.1 * 2 * 0.77 / 4 = 0.022 Eh/A, below typical bonded forces.
    if (z_wall <= 0.0)
        return 0.0;
    const double delta = 0.2 * z_wall;
    return 2.0 * 3.0 * kT / (delta * delta);
}

double PathCV::s(const Geometry& R) const
{
    double s_val = 0.0;
    evaluate(R, &s_val, nullptr, nullptr);
    return s_val;
}

double PathCV::z(const Geometry& R) const
{
    double z_val = 0.0;
    evaluate(R, nullptr, &z_val, nullptr);
    return z_val;
}
