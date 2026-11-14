/*
 * < Ablation: do persistent-homology barcodes carry energy information beyond ACSF? >
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
 * Claude Generated (Oct 2026). Energy-only experiment, not part of the curcuma CLI.
 *
 * Question: does a local Vietoris-Rips barcode embedding predict E_ref - E_base of a conformer
 * set better than ACSF alone? Everything is energy-only so that variants are comparable (a
 * barcode has no usable analytic gradient, see docs/DELTA_ML.md).
 *
 * Model (all variants identical, only the per-atom feature vectors differ):
 *   Delta(S) = sum_T alpha_T K(S,T) + c,  K(S,T) = sum_{i in S, j in T, Z_i = Z_j} exp(-|x_i-x_j|^2 / 2 sigma^2)
 * dense kernel ridge on structures (the dense limit of DeltaMLModel without forces), nested
 * cross-validation: hyper-parameters (sigma scale, lambda) are chosen inside each training fold.
 *
 * Variants: ACSF | PI0 (H0 death histogram) | PI1 (H1 persistence image) | PI (both) | RDF (pair
 * distance histogram of the same patch, same dimension as PI: the control that says whether
 * topology matters or only the distances in the patch) | ACSF+PI | ACSF+RDF | ACSF+PIperm (PI
 * rows shuffled across conformers: a null that keeps the feature statistics, destroys the link
 * to the energy).
 */
#include "src/capabilities/persistentdiagram.h"
#include "src/core/global.h"
#include "src/core/ml/descriptor.h"
#include "src/core/ml/ml_dataset.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <thread>

using namespace curcuma::ml;

namespace {

constexpr double EH_TO_KCAL = 627.5094740631;

struct Options {
    std::string dataset;
    double radius = 5.0;   ///< patch radius around each atom [Angstrom]
    int dim_max = 1;       ///< highest homology dimension (1 = H0 + H1)
    int folds = 5;
    int repeats = 3;
    unsigned seed = 7;
    int threads = 0;
    int max_structures = 0; ///< 0 = all (debug aid)
    std::string json_out;
};

// ---------------------------------------------------------------------------------------------
// Barcode / histogram features of one atom's patch
// ---------------------------------------------------------------------------------------------
constexpr int N_H0 = 16;           // death histogram bins
constexpr double H0_LO = 0.8, H0_HI = 5.0;
constexpr int N_PI_B = 8, N_PI_P = 8; // H1 image: birth x persistence
constexpr double PI_B_LO = 1.5, PI_B_HI = 5.0, PI_P_LO = 0.0, PI_P_HI = 1.5;
constexpr double PI_SIG_B = 0.30, PI_SIG_P = 0.20;
constexpr int N_RDF = N_H0 + N_PI_B * N_PI_P; // 80: same dimension as the PI feature
constexpr double RDF_LO = 0.8, RDF_HI = 6.0;

/// Triangular ("linear interpolation") histogram: the value is continuous in the input.
void addTriangular(double v, double lo, double hi, int nb, double w, double* bins)
{
    const double u = (v - lo) / (hi - lo) * nb - 0.5;
    const int i0 = static_cast<int>(std::floor(u));
    const double f = u - i0;
    if (i0 >= 0 && i0 < nb)
        bins[i0] += w * (1.0 - f);
    if (i0 + 1 >= 0 && i0 + 1 < nb)
        bins[i0 + 1] += w * f;
}

struct PerStructureFeatures {
    Matrix acsf, pi0, pi1, rdf; ///< natoms x dim
};

void patchFeatures(const Geometry& R, int centre, double radius, int dim_max, double* pi0, double* pi1, double* rdf)
{
    const int n = static_cast<int>(R.rows());
    std::vector<int> idx;
    for (int j = 0; j < n; ++j)
        if ((R.row(j) - R.row(centre)).norm() <= radius)
            idx.push_back(j);
    const int m = static_cast<int>(idx.size());
    if (m < 2)
        return;

    // compressed lower-triangular distance matrix, row i holds entries j < i
    std::vector<float> dist;
    dist.reserve(static_cast<size_t>(m) * (m - 1) / 2);
    for (int a = 1; a < m; ++a)
        for (int b = 0; b < a; ++b) {
            const double d = (R.row(idx[a]) - R.row(idx[b])).norm();
            dist.push_back(static_cast<float>(d));
            addTriangular(d, RDF_LO, RDF_HI, N_RDF, 1.0, rdf);
        }

    // A standalone tool has no ParameterRegistry defaults; an empty config leaves m_defaults null.
    PersistentDiagram pd(json{ { "ripser", { { "dimension", dim_max }, { "ratio", 1.0 } } } }, true);
    pd.setDimension(dim_max);
    pd.setDistanceMatrix(dist);
    const auto bars = pd.generateBarcodes();

    auto it0 = bars.find(0);
    if (it0 != bars.end())
        for (const auto& p : it0->second)
            if (std::isfinite(p.second))
                addTriangular(p.second, H0_LO, H0_HI, N_H0, 1.0, pi0);

    auto it1 = bars.find(1);
    if (it1 != bars.end())
        for (const auto& p : it1->second) {
            if (!std::isfinite(p.second))
                continue;
            const double b = p.first, pers = p.second - p.first;
            for (int ib = 0; ib < N_PI_B; ++ib) {
                const double cb = PI_B_LO + (PI_B_HI - PI_B_LO) * (ib + 0.5) / N_PI_B;
                for (int ip = 0; ip < N_PI_P; ++ip) {
                    const double cp = PI_P_LO + (PI_P_HI - PI_P_LO) * (ip + 0.5) / N_PI_P;
                    pi1[ib * N_PI_P + ip] += pers
                        * std::exp(-0.5 * ((cb - b) * (cb - b) / (PI_SIG_B * PI_SIG_B) + (cp - pers) * (cp - pers) / (PI_SIG_P * PI_SIG_P)));
                }
            }
        }
}

void parallelFor(int n, int threads, const std::function<void(int)>& f)
{
    std::atomic<int> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&]() {
            for (int i = next++; i < n; i = next++)
                f(i);
        });
    for (auto& th : pool)
        th.join();
}

// ---------------------------------------------------------------------------------------------
// Per-atom feature blocks, standardisation, kernels
// ---------------------------------------------------------------------------------------------
using Blocks = std::vector<Matrix>; ///< one natoms x d matrix per structure

/// Standardise per element and column over the whole set (unsupervised, same for every variant),
/// drop constant columns, and scale the block so that its mean squared row norm is 1.
Blocks standardise(const Blocks& in, const std::vector<int>& Z)
{
    const int ns = static_cast<int>(in.size()), na = static_cast<int>(Z.size());
    const int d = static_cast<int>(in[0].cols());
    Blocks out(ns, Matrix::Zero(na, d));
    std::set<int> elems(Z.begin(), Z.end());
    double sumsq = 0.0;
    long cnt = 0;
    for (int e : elems) {
        std::vector<int> rows;
        for (int i = 0; i < na; ++i)
            if (Z[i] == e)
                rows.push_back(i);
        Vector mean = Vector::Zero(d), var = Vector::Zero(d);
        for (int s = 0; s < ns; ++s)
            for (int i : rows)
                mean += in[s].row(i).transpose();
        mean /= static_cast<double>(ns * rows.size());
        for (int s = 0; s < ns; ++s)
            for (int i : rows)
                var += (in[s].row(i).transpose() - mean).cwiseAbs2();
        var /= static_cast<double>(ns * rows.size());
        for (int s = 0; s < ns; ++s)
            for (int i : rows)
                for (int c = 0; c < d; ++c)
                    if (var(c) > 1e-14) {
                        const double v = (in[s](i, c) - mean(c)) / std::sqrt(var(c));
                        out[s](i, c) = v;
                        sumsq += v * v;
                    }
        cnt += static_cast<long>(ns) * rows.size();
    }
    const double scale = 1.0 / std::sqrt(sumsq / std::max(1L, cnt));
    for (auto& m : out)
        m *= scale;
    return out;
}

Blocks concat(const Blocks& a, const Blocks& b)
{
    Blocks out(a.size());
    for (size_t s = 0; s < a.size(); ++s) {
        out[s].resize(a[s].rows(), a[s].cols() + b[s].cols());
        out[s] << a[s], b[s];
    }
    return out;
}

/// Median distance between random same-element atom pairs of different structures.
double medianDistance(const Blocks& X, const std::vector<int>& Z, unsigned seed)
{
    std::mt19937 rng(seed);
    const int ns = static_cast<int>(X.size()), na = static_cast<int>(Z.size());
    std::vector<double> d;
    for (int k = 0; k < 20000; ++k) {
        const int s = rng() % ns, t = rng() % ns, i = rng() % na;
        if (s == t)
            continue;
        d.push_back((X[s].row(i) - X[t].row(i)).norm()); // same atom index = same element and role
    }
    std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
    return d[d.size() / 2];
}

/// Kernels for several sigma values in one pass over the distance matrices.
std::vector<Matrix> sumKernels(const Blocks& X, const std::vector<int>& Z, const std::vector<double>& sigmas, int threads)
{
    const int ns = static_cast<int>(X.size());
    const int d = static_cast<int>(X[0].cols());
    std::vector<Matrix> K(sigmas.size(), Matrix::Zero(ns, ns));
    std::set<int> elems(Z.begin(), Z.end());
    for (int e : elems) {
        std::vector<int> rows;
        for (int i = 0; i < static_cast<int>(Z.size()); ++i)
            if (Z[i] == e)
                rows.push_back(i);
        const int ne = static_cast<int>(rows.size());
        Matrix E(ns * ne, d);
        for (int s = 0; s < ns; ++s)
            for (int r = 0; r < ne; ++r)
                E.row(s * ne + r) = X[s].row(rows[r]);
        Vector nrm = E.rowwise().squaredNorm();
        parallelFor(ns, threads, [&](int s) {
            Matrix G = E.middleRows(s * ne, ne) * E.transpose(); // ne x (ns*ne)
            for (size_t k = 0; k < sigmas.size(); ++k) {
                const double inv = -1.0 / (2.0 * sigmas[k] * sigmas[k]);
                for (int t = 0; t < ns; ++t) {
                    double acc = 0.0;
                    for (int j = 0; j < ne; ++j)
                        for (int i = 0; i < ne; ++i) {
                            const double d2 = nrm(s * ne + i) + nrm(t * ne + j) - 2.0 * G(i, t * ne + j);
                            acc += std::exp(inv * std::max(0.0, d2));
                        }
                    K[k](s, t) += acc;
                }
            }
        });
    }
    return K;
}

// ---------------------------------------------------------------------------------------------
// Kernel ridge with nested cross-validation
// ---------------------------------------------------------------------------------------------
double centeredMAE(const Vector& err)
{
    return (err.array() - err.mean()).abs().mean();
}

Vector krrPredict(const Matrix& K, const std::vector<int>& train, const std::vector<int>& test, const Vector& y, double lambda_rel)
{
    const int n = static_cast<int>(train.size());
    Matrix Ktt(n, n);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            Ktt(a, b) = K(train[a], train[b]);
    const double scale = Ktt.diagonal().mean();
    Vector yt(n);
    for (int a = 0; a < n; ++a)
        yt(a) = y(train[a]);
    const double ymean = yt.mean();
    Matrix A = Ktt;
    A.diagonal().array() += lambda_rel * scale;
    Vector rhs = yt;
    rhs.array() -= ymean;
    Vector alpha = A.ldlt().solve(rhs);
    Vector pred(test.size());
    for (size_t a = 0; a < test.size(); ++a) {
        double v = ymean;
        for (int b = 0; b < n; ++b)
            v += alpha(b) * K(test[a], train[b]);
        pred(a) = v;
    }
    return pred;
}

struct Choice {
    int sigma = 0;
    double lambda = 0.0;
};

Choice selectHyper(const std::vector<Matrix>& Ks, const std::vector<int>& train, const Vector& y, const std::vector<double>& lambdas, unsigned seed)
{
    const int inner = 4;
    std::vector<int> perm = train;
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);
    Choice best;
    double bestErr = 1e300;
    for (size_t k = 0; k < Ks.size(); ++k)
        for (double lam : lambdas) {
            double err = 0.0;
            for (int f = 0; f < inner; ++f) {
                std::vector<int> tr, te;
                for (size_t a = 0; a < perm.size(); ++a)
                    (static_cast<int>(a) % inner == f ? te : tr).push_back(perm[a]);
                Vector p = krrPredict(Ks[k], tr, te, y, lam);
                Vector e(te.size());
                for (size_t a = 0; a < te.size(); ++a)
                    e(a) = p(a) - y(te[a]);
                err += centeredMAE(e) * te.size();
            }
            if (err < bestErr) {
                bestErr = err;
                best = { static_cast<int>(k), lam };
            }
        }
    return best;
}

double spearman(const Vector& a, const Vector& b)
{
    auto ranks = [](const Vector& v) {
        const int n = static_cast<int>(v.size());
        std::vector<int> o(n);
        std::iota(o.begin(), o.end(), 0);
        std::sort(o.begin(), o.end(), [&](int i, int j) { return v(i) < v(j); });
        Vector r(n);
        for (int i = 0; i < n;) {
            int j = i;
            while (j + 1 < n && v(o[j + 1]) == v(o[i]))
                ++j;
            for (int k = i; k <= j; ++k)
                r(o[k]) = 0.5 * (i + j);
            i = j + 1;
        }
        return r;
    };
    Vector ra = ranks(a), rb = ranks(b);
    ra.array() -= ra.mean();
    rb.array() -= rb.mean();
    return ra.dot(rb) / std::sqrt(ra.squaredNorm() * rb.squaredNorm());
}

struct VariantResult {
    std::string name;
    Vector oof_err; ///< out-of-fold prediction error [kcal/mol], mean over repeats
    double mae = 0, rmse = 0, spearman = 0;
    std::map<int, int> sigma_votes;
    double lambda_median = 0;
};

VariantResult evaluate(const std::string& name, const std::vector<Matrix>& Ks, const Vector& y, const Vector& e_base, const Vector& e_ref, const Options& opt)
{
    const int ns = static_cast<int>(y.size());
    const std::vector<double> lambdas = { 1e-8, 1e-6, 1e-4, 1e-2, 1.0 };
    VariantResult r;
    r.name = name;
    Vector err_sum = Vector::Zero(ns), pred_sum = Vector::Zero(ns);
    std::vector<double> lams;
    for (int rep = 0; rep < opt.repeats; ++rep) {
        std::vector<int> perm(ns);
        std::iota(perm.begin(), perm.end(), 0);
        std::mt19937 rng(opt.seed + 101 * rep);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int f = 0; f < opt.folds; ++f) {
            std::vector<int> tr, te;
            for (int a = 0; a < ns; ++a)
                (a % opt.folds == f ? te : tr).push_back(perm[a]);
            Choice c = selectHyper(Ks, tr, y, lambdas, opt.seed + 13 * f + rep);
            r.sigma_votes[c.sigma]++;
            lams.push_back(c.lambda);
            Vector p = krrPredict(Ks[c.sigma], tr, te, y, c.lambda);
            for (size_t a = 0; a < te.size(); ++a) {
                err_sum(te[a]) += p(a) - y(te[a]);
                pred_sum(te[a]) += p(a);
            }
        }
    }
    r.oof_err = err_sum / opt.repeats;
    Vector pred = pred_sum / opt.repeats;
    r.mae = centeredMAE(r.oof_err);
    r.rmse = std::sqrt((r.oof_err.array() - r.oof_err.mean()).square().mean());
    r.spearman = spearman(e_base + pred, e_ref);
    std::sort(lams.begin(), lams.end());
    r.lambda_median = lams[lams.size() / 2];
    return r;
}

/// 95 % bootstrap interval of mean(|e - mean e|) (and of the paired difference to a reference).
std::pair<double, double> bootstrapCI(const Vector& e, const Vector* ref, unsigned seed)
{
    const int n = static_cast<int>(e.size());
    std::mt19937 rng(seed);
    std::vector<double> v;
    for (int b = 0; b < 2000; ++b) {
        Vector s(n), t(n);
        for (int i = 0; i < n; ++i) {
            const int k = rng() % n;
            s(i) = e(k);
            if (ref)
                t(i) = (*ref)(k);
        }
        v.push_back(centeredMAE(s) - (ref ? centeredMAE(t) : 0.0));
    }
    std::sort(v.begin(), v.end());
    return { v[static_cast<size_t>(0.025 * v.size())], v[static_cast<size_t>(0.975 * v.size())] };
}

} // namespace

int main(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (a == "--radius") opt.radius = std::stod(next());
        else if (a == "--dim") opt.dim_max = std::stoi(next());
        else if (a == "--folds") opt.folds = std::stoi(next());
        else if (a == "--repeats") opt.repeats = std::stoi(next());
        else if (a == "--seed") opt.seed = std::stoul(next());
        else if (a == "--threads") opt.threads = std::stoi(next());
        else if (a == "--max-structures") opt.max_structures = std::stoi(next());
        else if (a == "--json") opt.json_out = next();
        else opt.dataset = a;
    }
    if (opt.dataset.empty()) {
        std::cerr << "usage: ml_barcode_ablation <dataset.json from -activelearn> [--radius 5.0] [--dim 1] [--folds 5] [--repeats 3]\n"
                     "                           [--seed 7] [--threads N] [--max-structures N] [--json out.json]\n";
        return 2;
    }
    if (opt.threads <= 0)
        opt.threads = std::max(1u, std::thread::hardware_concurrency());

    MLDataset data;
    if (!data.load(opt.dataset)) {
        std::cerr << "cannot read dataset " << opt.dataset << "\n";
        return 1;
    }
    int ns = static_cast<int>(data.size());
    if (opt.max_structures > 0)
        ns = std::min(ns, opt.max_structures);
    const std::vector<int> Z = data[0].elements;
    for (int s = 0; s < ns; ++s)
        if (data[s].elements != Z) {
            std::cerr << "structure " << s << " differs in element order; the ablation needs one molecule\n";
            return 1;
        }
    std::cout << "structures " << ns << ", atoms " << Z.size() << ", patch radius " << opt.radius << " A, homology dim <= " << opt.dim_max
              << ", threads " << opt.threads << "\n";

    Vector y(ns), e_base(ns), e_ref(ns);
    for (int s = 0; s < ns; ++s) {
        e_base(s) = data[s].energy_base * EH_TO_KCAL;
        e_ref(s) = data[s].energy_ref * EH_TO_KCAL;
        y(s) = e_ref(s) - e_base(s);
    }

    // ----- features ------------------------------------------------------------------------
    auto t0 = std::chrono::steady_clock::now();
    ACSFDescriptor acsf(ACSFConfig::defaultFor(data.species()));
    Blocks B_acsf(ns), B_pi0(ns), B_pi1(ns), B_rdf(ns);
    const int na = static_cast<int>(Z.size());
    for (int s = 0; s < ns; ++s) {
        acsf.compute(Z, data[s].geometry, B_acsf[s]);
        B_pi0[s] = Matrix::Zero(na, N_H0);
        B_pi1[s] = Matrix::Zero(na, N_PI_B * N_PI_P);
        B_rdf[s] = Matrix::Zero(na, N_RDF);
    }
    parallelFor(ns * na, opt.threads, [&](int k) {
        const int s = k / na, i = k % na;
        patchFeatures(data[s].geometry, i, opt.radius, opt.dim_max, B_pi0[s].row(i).data(), B_pi1[s].row(i).data(), B_rdf[s].row(i).data());
    });
    std::cout << "features: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s ("
              << ns * na << " patches)\n";

    {
        // sanity: an empty barcode would make every PI variant a constant
        double h0 = 0, h1 = 0, rd = 0;
        for (int s = 0; s < ns; ++s)
            for (int i = 0; i < na; ++i) {
                h0 += B_pi0[s].row(i).sum();
                h1 += B_pi1[s].row(i).sum();
                rd += B_rdf[s].row(i).sum();
            }
        const double np = static_cast<double>(ns) * na;
        std::cout << "mean per patch: H0 death count " << h0 / np << ", H1 image mass " << h1 / np << ", pair count in RDF range " << rd / np << "\n";
    }

    // ----- variants ------------------------------------------------------------------------
    Blocks S_acsf = standardise(B_acsf, Z), S_pi0 = standardise(B_pi0, Z), S_pi1 = standardise(B_pi1, Z), S_rdf = standardise(B_rdf, Z);
    Blocks S_pi = standardise(concat(B_pi0, B_pi1), Z);
    Blocks S_pi_perm(ns);
    {
        std::vector<int> perm(ns);
        std::iota(perm.begin(), perm.end(), 0);
        std::mt19937 rng(opt.seed + 999);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int s = 0; s < ns; ++s)
            S_pi_perm[s] = S_pi[perm[s]];
    }
    // concatenation of two unit-scale blocks: each keeps weight 1/sqrt(2) so the sum is again unit scale
    auto both = [&](const Blocks& a, const Blocks& b) {
        Blocks c = concat(a, b);
        for (auto& m : c)
            m *= 1.0 / std::sqrt(2.0);
        return c;
    };
    // small admixture: ACSF keeps weight 1, the second block weight w. Equal weighting (above) dilutes
    // ACSF, which could make "ACSF + X" look worse than "ACSF" for a reason unrelated to X.
    auto mix = [&](const Blocks& a, const Blocks& b, double w) {
        Blocks c = concat(a, b);
        const int da = static_cast<int>(a[0].cols());
        for (auto& m : c)
            m.rightCols(m.cols() - da) *= w;
        return c;
    };
    struct Variant {
        std::string name;
        Blocks X;
    };
    std::vector<Variant> variants = {
        { "ACSF", S_acsf },
        { "PI0 (H0)", S_pi0 },
        { "PI1 (H1)", S_pi1 },
        { "PI (H0+H1)", S_pi },
        { "RDF control", S_rdf },
        { "ACSF+PI", both(S_acsf, S_pi) },
        { "ACSF+RDF", both(S_acsf, S_rdf) },
        { "ACSF+PI shuffled", both(S_acsf, S_pi_perm) },
        { "ACSF+0.5PI", mix(S_acsf, S_pi, 0.5) },
        { "ACSF+0.5RDF", mix(S_acsf, S_rdf, 0.5) },
        { "ACSF+0.5PI shuf", mix(S_acsf, S_pi_perm, 0.5) },
        { "ACSF+0.25PI", mix(S_acsf, S_pi, 0.25) },
        { "ACSF+0.25RDF", mix(S_acsf, S_rdf, 0.25) },
        { "ACSF+0.25PI shuf", mix(S_acsf, S_pi_perm, 0.25) },
    };

    std::vector<VariantResult> results;
    for (const auto& v : variants) {
        auto tk = std::chrono::steady_clock::now();
        const double med = medianDistance(v.X, Z, opt.seed);
        std::vector<double> sigmas = { 0.5 * med, 1.0 * med, 2.0 * med };
        std::vector<Matrix> Ks = sumKernels(v.X, Z, sigmas, opt.threads);
        results.push_back(evaluate(v.name, Ks, y, e_base, e_ref, opt));
        std::cout << "  " << std::left << std::setw(18) << v.name << " kernel+CV " << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - tk).count() << " s\n"
                  << std::flush;
    }

    // ----- report --------------------------------------------------------------------------
    const double base_mae = centeredMAE(y.array() - y.mean());
    std::cout << "\nenergy-only, " << opt.folds << "-fold nested CV x " << opt.repeats << " repeats, n = " << ns
              << " conformers of ONE molecule\n"
              << "target: E_ref - E_base [kcal/mol], errors centred (a constant offset cancels in relative energies)\n"
              << "no model (E_base only): centred MAE " << std::setprecision(3) << base_mae << ", Spearman(E_base, E_ref) "
              << spearman(e_base, e_ref) << "\n\n";
    std::cout << std::left << std::setw(18) << "variant" << std::right << std::setw(10) << "MAE" << std::setw(20) << "95% CI (MAE)" << std::setw(24)
              << "dMAE vs ACSF [95% CI]" << std::setw(9) << "RMSE" << std::setw(10) << "Spearman" << "\n";
    json jout = json::array();
    const VariantResult& ref = results[0];
    for (const auto& r : results) {
        auto ci = bootstrapCI(r.oof_err, nullptr, opt.seed);
        std::ostringstream c1, c2;
        c1 << "[" << std::setprecision(3) << ci.first << "," << ci.second << "]";
        if (&r != &ref) {
            auto d = bootstrapCI(r.oof_err, &ref.oof_err, opt.seed);
            c2 << std::showpos << std::setprecision(3) << (r.mae - ref.mae) << " [" << d.first << "," << d.second << "]";
        } else {
            c2 << "-";
        }
        std::cout << std::left << std::setw(18) << r.name << std::right << std::fixed << std::setprecision(3) << std::setw(10) << r.mae << std::setw(20)
                  << c1.str() << std::setw(24) << c2.str() << std::setw(9) << r.rmse << std::setw(10) << r.spearman << "\n";
        jout.push_back({ { "variant", r.name }, { "mae", r.mae }, { "rmse", r.rmse }, { "spearman", r.spearman }, { "ci_lo", ci.first }, { "ci_hi", ci.second } });
    }
    std::cout << "\nreading: PI must beat RDF control and ACSF+PI must beat ACSF outside the paired interval; ACSF+PI shuffled must NOT.\n";
    if (!opt.json_out.empty()) {
        std::ofstream(opt.json_out) << json{ { "n", ns }, { "radius", opt.radius }, { "dim_max", opt.dim_max }, { "results", jout } }.dump(2);
    }
    return 0;
}
