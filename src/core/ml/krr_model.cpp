/*
 * < Sparse Kernel-Ridge model for Delta-Machine-Learning >
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
 * Claude Generated (July 2026): regression layer of the Delta-ML work package.
 *
 * Literature:
 *   A. P. Bartok, G. Csanyi, Int. J. Quantum Chem. 115, 1051 (2015) -- sparse GP / GAP
 *   R. Ramakrishnan et al., J. Chem. Theory Comput. 11, 2087 (2015) -- delta-learning
 */

#include "krr_model.h"
#include "kernel.h"

#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <random>
#include <set>

namespace curcuma {
namespace ml {

namespace {

    /*! \brief Ordinal ranks of a value list, ties broken by index (adequate for floats). */
    std::vector<double> ranks(const std::vector<double>& v)
    {
        std::vector<std::size_t> order(v.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&v](std::size_t a, std::size_t b) { return v[a] < v[b]; });
        std::vector<double> r(v.size(), 0.0);
        for (std::size_t k = 0; k < order.size(); ++k)
            r[order[k]] = static_cast<double>(k);
        return r;
    }

    double pearson(const std::vector<double>& a, const std::vector<double>& b)
    {
        if (a.size() < 2 || a.size() != b.size())
            return 0.0;
        const double n = static_cast<double>(a.size());
        const double ma = std::accumulate(a.begin(), a.end(), 0.0) / n;
        const double mb = std::accumulate(b.begin(), b.end(), 0.0) / n;
        double num = 0.0, da = 0.0, db = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const double x = a[i] - ma;
            const double y = b[i] - mb;
            num += x * y;
            da += x * x;
            db += y * y;
        }
        if (da <= 0.0 || db <= 0.0)
            return 0.0;
        return num / std::sqrt(da * db);
    }

    double spearman(const std::vector<double>& a, const std::vector<double>& b)
    {
        return pearson(ranks(a), ranks(b));
    }

} // namespace

// =================================================================================
// Helpers
// =================================================================================

void DeltaMLModel::standardise(Matrix& X) const
{
    if (m_feature_mean.size() != X.cols())
        return;
    X.rowwise() -= m_feature_mean.transpose();
    X = X * m_feature_scale.cwiseInverse().asDiagonal();
}

double DeltaMLModel::baselineEnergy(const std::vector<int>& elements) const
{
    double e = 0.0;
    for (int z : elements) {
        const auto it = m_element_energy.find(z);
        if (it != m_element_energy.end())
            e += it->second;
    }
    return e;
}

Vector DeltaMLModel::kernelRow(const Vector& x, int element, int& start, int& count) const
{
    const auto it = m_element_block.find(element);
    if (it == m_element_block.end()) {
        start = 0;
        count = 0;
        return Vector();
    }
    start = it->second.first;
    count = it->second.second;

    Vector d2(count);
    for (int m = 0; m < count; ++m)
        d2(m) = (m_centers.row(start + m).transpose() - x).squaredNorm();
    return (d2 * (-1.0 / (2.0 * m_sigma * m_sigma))).array().exp();
}

void DeltaMLModel::setAnchor(const std::vector<int>& elements, const Geometry& geometry)
{
    m_anchor_elements = elements;
    m_anchor_geometry = geometry;
}

// =================================================================================
// Training
// =================================================================================

bool DeltaMLModel::train(
    const MLDataset& data, const DeltaMLTrainingOptions& opt, DeltaMLTrainingReport& report)
{
    const auto t0 = std::chrono::high_resolution_clock::now();
    report = DeltaMLTrainingReport{};
    m_trained = false;

    if (data.empty()) {
        CurcumaLogger::error("DeltaMLModel: cannot train on an empty dataset");
        return false;
    }

    // --- descriptor -------------------------------------------------------------
    m_acsf = ACSFConfig::defaultFor(data.species());
    std::string cfg_error;
    if (!m_acsf.validate(cfg_error)) {
        CurcumaLogger::error_fmt("DeltaMLModel: invalid descriptor config: {}", cfg_error);
        return false;
    }
    const ACSFDescriptor desc(m_acsf);
    const int F = m_acsf.featureCount();
    const int N = static_cast<int>(data.size());

    // --- pass 1: descriptors ----------------------------------------------------
    std::vector<Matrix> X(N);
    int n_env = 0;
    for (int s = 0; s < N; ++s) {
        desc.compute(data[s].elements, data[s].geometry, X[s], nullptr);
        n_env += data[s].AtomCount();
    }
    if (n_env == 0) {
        CurcumaLogger::error("DeltaMLModel: dataset contains no atoms");
        return false;
    }

    // --- feature standardisation ------------------------------------------------
    m_feature_mean = Vector::Zero(F);
    for (int s = 0; s < N; ++s)
        m_feature_mean += X[s].colwise().sum().transpose();
    m_feature_mean /= static_cast<double>(n_env);

    Vector variance = Vector::Zero(F);
    for (int s = 0; s < N; ++s) {
        Matrix centred = X[s];
        centred.rowwise() -= m_feature_mean.transpose();
        variance += centred.array().square().colwise().sum().matrix().transpose();
    }
    variance /= static_cast<double>(n_env);

    m_feature_scale = Vector::Ones(F);
    for (int f = 0; f < F; ++f) {
        // A constant feature carries no information. Leaving its scale at one keeps the
        // centred value at zero; dividing by a tiny standard deviation would instead
        // amplify pure rounding noise into a dominant coordinate.
        if (variance(f) > 1e-12)
            m_feature_scale(f) = std::sqrt(variance(f));
    }

    for (int s = 0; s < N; ++s)
        standardise(X[s]);

    // --- kernel centres: farthest point sampling per element --------------------
    std::map<int, std::vector<std::pair<int, int>>> by_element; // Z -> (structure, atom)
    for (int s = 0; s < N; ++s)
        for (int i = 0; i < data[s].AtomCount(); ++i)
            by_element[data[s].elements[i]].emplace_back(s, i);

    int m_target = (opt.n_centers > 0) ? opt.n_centers : n_env;
    if (m_target > n_env)
        m_target = n_env;

    // The normal-equation matrix is M x M and the fit cost grows with M^2, so an
    // accidental "use everything" on a large dataset must fail loudly, not swap to death.
    const int kCenterLimit = 12000;
    if (m_target > kCenterLimit) {
        CurcumaLogger::error_fmt(
            "DeltaMLModel: {} kernel centres requested ({} environments available). The "
            "normal equations would need {:.1f} GB. Set n_centers to at most {}.",
            m_target, n_env, static_cast<double>(m_target) * m_target * 8.0 / 1e9, kCenterLimit);
        return false;
    }

    std::vector<int> center_elements;
    std::vector<std::pair<int, int>> center_source;
    center_elements.reserve(m_target);
    center_source.reserve(m_target);

    for (const auto& entry : by_element) {
        const int z = entry.first;
        const auto& envs = entry.second;
        const int n_z = static_cast<int>(envs.size());

        int m_z = static_cast<int>(std::llround(
            static_cast<double>(m_target) * static_cast<double>(n_z) / static_cast<double>(n_env)));
        m_z = std::max(1, std::min(m_z, n_z));

        // Gather this element's environments into one matrix for the sampling.
        Matrix E(n_z, F);
        for (int k = 0; k < n_z; ++k)
            E.row(k) = X[envs[k].first].row(envs[k].second);

        // Deterministic start: the environment farthest from the element centroid, so the
        // selection does not depend on the order structures happened to arrive in.
        const Vector centroid = E.colwise().mean().transpose();
        int current = 0;
        double best = -1.0;
        for (int k = 0; k < n_z; ++k) {
            const double d = (E.row(k).transpose() - centroid).squaredNorm();
            if (d > best) {
                best = d;
                current = k;
            }
        }

        std::vector<char> taken(n_z, 0);
        Vector mind = Vector::Constant(n_z, std::numeric_limits<double>::max());
        for (int pick = 0; pick < m_z; ++pick) {
            taken[current] = 1;
            center_elements.push_back(z);
            center_source.push_back(envs[current]);

            const Vector c = E.row(current).transpose();
            for (int k = 0; k < n_z; ++k) {
                if (taken[k])
                    continue;
                const double d = (E.row(k).transpose() - c).squaredNorm();
                if (d < mind(k))
                    mind(k) = d;
            }
            // Next centre: the environment worst covered so far.
            double far = -1.0;
            int next = -1;
            for (int k = 0; k < n_z; ++k) {
                if (!taken[k] && mind(k) > far) {
                    far = mind(k);
                    next = k;
                }
            }
            if (next < 0)
                break;
            current = next;
        }
    }

    const int M = static_cast<int>(center_elements.size());
    m_centers = Matrix::Zero(M, F);
    for (int m = 0; m < M; ++m)
        m_centers.row(m) = X[center_source[m].first].row(center_source[m].second);
    m_center_elements = center_elements;

    m_element_block.clear();
    for (int m = 0; m < M; ++m) {
        const int z = m_center_elements[m];
        auto it = m_element_block.find(z);
        if (it == m_element_block.end())
            m_element_block[z] = { m, 1 };
        else
            it->second.second += 1;
    }

    // --- kernel width -----------------------------------------------------------
    if (opt.kernel_sigma > 0.0) {
        m_sigma = opt.kernel_sigma;
    } else {
        // Median heuristic: sigma^2 = median(|x-y|^2)/2, evaluated on a bounded random
        // sample of same-element centre pairs.
        std::mt19937 rng(opt.seed);
        std::vector<double> d2;
        d2.reserve(4096);
        if (M > 1) {
            std::uniform_int_distribution<int> pick(0, M - 1);
            for (int t = 0; t < 4096; ++t) {
                const int a = pick(rng);
                const int b = pick(rng);
                if (a == b || m_center_elements[a] != m_center_elements[b])
                    continue;
                d2.push_back((m_centers.row(a) - m_centers.row(b)).squaredNorm());
            }
        }
        if (d2.empty()) {
            m_sigma = 1.0;
        } else {
            std::nth_element(d2.begin(), d2.begin() + d2.size() / 2, d2.end());
            m_sigma = std::sqrt(std::max(1e-12, d2[d2.size() / 2] / 2.0));
        }
    }

    // --- element-count baseline -------------------------------------------------
    // Absorbs the constant offset between two methods with different energy zeros.
    // Gradient-free by construction, so forces are untouched.
    const std::vector<int> species = m_acsf.species;
    const int S = static_cast<int>(species.size());
    Matrix C = Matrix::Zero(N, S);
    Vector y_delta = Vector::Zero(N);
    for (int s = 0; s < N; ++s) {
        for (int z : data[s].elements) {
            const int si = m_acsf.speciesIndex(z);
            if (si >= 0)
                C(s, si) += 1.0;
        }
        y_delta(s) = data[s].delta();
    }
    // Rank-revealing: for a single molecule every row of C is identical, so the system is
    // rank one and only the intercept is determined. QR with column pivoting handles that
    // without producing garbage.
    const Vector c_elem = C.colPivHouseholderQr().solve(y_delta);
    m_element_energy.clear();
    for (int si = 0; si < S; ++si)
        m_element_energy[species[si]] = c_elem(si);
    const Vector y_res = y_delta - C * c_elem;

    // --- normal equations -------------------------------------------------------
    Matrix A = Matrix::Zero(M, M);
    Vector b = Vector::Zero(M);
    const double wE = 1.0 / (opt.sigma_energy * opt.sigma_energy);
    const double wF = 1.0 / (opt.sigma_force * opt.sigma_force);
    const double inv_sigma2 = 1.0 / (m_sigma * m_sigma);
    const Vector inv_scale = m_feature_scale.cwiseInverse();

    const bool want_forces = opt.use_forces && data.hasGradients();

    for (int s = 0; s < N; ++s) {
        const int natoms = data[s].AtomCount();

        // Energy row: phi(m) = sum_i k(x_i, z_m)
        Vector phi = Vector::Zero(M);
        for (int i = 0; i < natoms; ++i) {
            int start = 0, count = 0;
            const Vector xi = X[s].row(i).transpose();
            const Vector k = kernelRow(xi, data[s].elements[i], start, count);
            if (count > 0)
                phi.segment(start, count) += k;
        }
        A.selfadjointView<Eigen::Lower>().rankUpdate(phi, wE);
        b += wE * phi * y_res(s);
        report.n_energy_rows += 1;

        if (!want_forces)
            continue;

        // Gradient rows. The compact Jacobian has to be recomputed here because storing it
        // for every structure would dominate memory; the geometric loops are cheap.
        Matrix Xraw;
        std::vector<AtomJacobian> jac;
        desc.compute(data[s].elements, data[s].geometry, Xraw, &jac);

        Matrix PhiG = Matrix::Zero(3 * natoms, M);
        for (int i = 0; i < natoms; ++i) {
            const auto it = m_element_block.find(data[s].elements[i]);
            if (it == m_element_block.end())
                continue;
            const int start = it->second.first;
            const int count = it->second.second;

            const Vector xi = X[s].row(i).transpose();
            const Matrix Zb = m_centers.middleRows(start, count);

            Vector k(count);
            for (int m = 0; m < count; ++m)
                k(m) = (Zb.row(m).transpose() - xi).squaredNorm();
            k = (k * (-0.5 * inv_sigma2)).array().exp();

            // dk(x_i,z_m)/dx_i = -(k_m / sigma^2) (x_i - z_m), stacked over m.
            Matrix Dk = (k * xi.transpose() - k.asDiagonal() * Zb) * (-inv_sigma2);
            // Chain out of the standardised feature space into raw descriptor space.
            Dk = Dk * inv_scale.asDiagonal();

            // (3*n_partner x F) * (F x count) -> contribution of atom i to the force rows
            const Matrix contrib = jac[i].d * Dk.transpose();
            for (std::size_t p = 0; p < jac[i].partners.size(); ++p) {
                const int atom = jac[i].partners[p];
                for (int g = 0; g < 3; ++g)
                    PhiG.row(3 * atom + g).segment(start, count)
                        += contrib.row(3 * static_cast<int>(p) + g);
            }
        }

        const Matrix dG = data[s].deltaGradient();
        Vector yg(3 * natoms);
        for (int c2 = 0; c2 < natoms; ++c2)
            for (int g = 0; g < 3; ++g)
                yg(3 * c2 + g) = dG(c2, g);

        A.selfadjointView<Eigen::Lower>().rankUpdate(PhiG.transpose(), wF);
        b += wF * PhiG.transpose() * yg;
        report.n_force_rows += 3 * natoms;
    }

    // rankUpdate only touched the lower triangle; materialise the full symmetric matrix
    // through a temporary, because writing A from a view of A itself would alias.
    {
        const Matrix full = A.selfadjointView<Eigen::Lower>();
        A = full;
    }

    // --- regulariser and solve --------------------------------------------------
    // The RKHS penalty lambda * beta^T K_MM beta is what keeps the fit from interpolating
    // label noise into a rippled surface -- invisible in the energy error, fatal for MD.
    Matrix Kmm = Matrix::Zero(M, M);
    for (const auto& blk : m_element_block) {
        const int start = blk.second.first;
        const int count = blk.second.second;
        Kmm.block(start, start, count, count)
            = gaussianKernel(m_centers.middleRows(start, count), m_centers.middleRows(start, count), m_sigma);
    }
    A += opt.lambda * Kmm;

    double jitter = 0.0;
    bool solved = false;
    for (int attempt = 0; attempt < 6 && !solved; ++attempt) {
        Matrix Aj = A;
        if (jitter > 0.0)
            Aj.diagonal().array() += jitter;

        Eigen::LDLT<Matrix> ldlt(Aj);
        if (ldlt.info() == Eigen::Success) {
            m_beta = ldlt.solve(b);
            if (m_beta.allFinite()) {
                solved = true;
                break;
            }
        }
        jitter = (jitter == 0.0) ? 1e-12 : jitter * 100.0;
    }

    if (!solved) {
        CurcumaLogger::error("DeltaMLModel: linear solve failed even with jitter; "
                             "try a larger lambda or fewer centres");
        return false;
    }

    m_trained = true;

    report.n_structures = N;
    report.n_environments = n_env;
    report.n_centers = M;
    report.n_features = F;
    report.kernel_sigma = m_sigma;
    report.lambda_used = opt.lambda;
    report.jitter_used = jitter;
    report.solver_ok = true;
    report.train_errors = evaluate(data);

    const auto t1 = std::chrono::high_resolution_clock::now();
    report.train_seconds = std::chrono::duration<double>(t1 - t0).count();

    return true;
}

// =================================================================================
// Prediction
// =================================================================================

double DeltaMLModel::predictEnergy(const std::vector<int>& elements, const Geometry& geometry) const
{
    if (!m_trained)
        return 0.0;

    const ACSFDescriptor desc(m_acsf);
    Matrix X;
    desc.compute(elements, geometry, X, nullptr);
    standardise(X);

    double energy = baselineEnergy(elements);
    for (int i = 0; i < static_cast<int>(elements.size()); ++i) {
        int start = 0, count = 0;
        const Vector k = kernelRow(X.row(i).transpose(), elements[i], start, count);
        if (count > 0)
            energy += k.dot(m_beta.segment(start, count));
    }
    return energy;
}

double DeltaMLModel::predictEnergyAndGradient(
    const std::vector<int>& elements, const Geometry& geometry, Matrix& gradient) const
{
    const int natoms = static_cast<int>(elements.size());
    gradient = Matrix::Zero(natoms, 3);
    if (!m_trained)
        return 0.0;

    const ACSFDescriptor desc(m_acsf);
    Matrix X;
    std::vector<AtomJacobian> jac;
    desc.compute(elements, geometry, X, &jac);

    Matrix Xs = X;
    standardise(Xs);

    const double inv_sigma2 = 1.0 / (m_sigma * m_sigma);
    const Vector inv_scale = m_feature_scale.cwiseInverse();

    double energy = baselineEnergy(elements);

    for (int i = 0; i < natoms; ++i) {
        const auto it = m_element_block.find(elements[i]);
        if (it == m_element_block.end())
            continue;
        const int start = it->second.first;
        const int count = it->second.second;

        const Vector xi = Xs.row(i).transpose();
        const Matrix Zb = m_centers.middleRows(start, count);
        const Vector beta_b = m_beta.segment(start, count);

        Vector k(count);
        for (int m = 0; m < count; ++m)
            k(m) = (Zb.row(m).transpose() - xi).squaredNorm();
        k = (k * (-0.5 * inv_sigma2)).array().exp();

        const double s = k.dot(beta_b);
        energy += s;

        // Adjoint in standardised feature space, then chained to raw descriptor space.
        const Vector bk = beta_b.cwiseProduct(k);
        Vector w = (s * xi - Zb.transpose() * bk) * (-inv_sigma2);
        w = w.cwiseProduct(inv_scale);

        // Scatter the contraction with the compact Jacobian onto the atoms involved.
        const Vector contrib = jac[i].d * w;
        for (std::size_t p = 0; p < jac[i].partners.size(); ++p) {
            const int atom = jac[i].partners[p];
            for (int g = 0; g < 3; ++g)
                gradient(atom, g) += contrib(3 * static_cast<int>(p) + g);
        }
    }

    return energy;
}

double DeltaMLModel::uncertainty(const std::vector<int>& elements, const Geometry& geometry) const
{
    if (!m_trained)
        return 1.0;

    const ACSFDescriptor desc(m_acsf);
    Matrix X;
    desc.compute(elements, geometry, X, nullptr);
    standardise(X);

    double worst = 0.0;
    for (int i = 0; i < static_cast<int>(elements.size()); ++i) {
        int start = 0, count = 0;
        const Vector k = kernelRow(X.row(i).transpose(), elements[i], start, count);
        // An element the model never saw is maximally unfamiliar.
        const double closeness = (count > 0) ? k.maxCoeff() : 0.0;
        worst = std::max(worst, 1.0 - closeness);
    }
    return worst;
}

// =================================================================================
// Evaluation
// =================================================================================

DeltaMLErrors DeltaMLModel::evaluate(const MLDataset& data) const
{
    DeltaMLErrors e;
    e.n_structures = static_cast<int>(data.size());
    if (data.empty())
        return e;

    std::vector<double> pred, ref, base;
    pred.reserve(data.size());
    ref.reserve(data.size());
    base.reserve(data.size());

    const bool with_forces = data.hasGradients();
    double f_num = 0.0, f_base_num = 0.0;
    long long f_count = 0;

    for (std::size_t s = 0; s < data.size(); ++s) {
        const MLStructure& st = data[s];
        double corr = 0.0;
        if (with_forces) {
            Matrix g;
            corr = predictEnergyAndGradient(st.elements, st.geometry, g);
            const Matrix gp = st.gradient_base + g;
            f_num += (gp - st.gradient_ref).array().square().sum();
            f_base_num += (st.gradient_base - st.gradient_ref).array().square().sum();
            f_count += static_cast<long long>(st.AtomCount()) * 3;
        } else {
            corr = predictEnergy(st.elements, st.geometry);
        }
        pred.push_back(st.energy_base + corr);
        ref.push_back(st.energy_ref);
        base.push_back(st.energy_base);
    }

    const double n = static_cast<double>(pred.size());
    double sae = 0.0, sse = 0.0, mean_err = 0.0;
    for (std::size_t i = 0; i < pred.size(); ++i) {
        const double d = pred[i] - ref[i];
        sae += std::abs(d);
        sse += d * d;
        mean_err += d;
    }
    mean_err /= n;
    e.energy_mae = sae / n;
    e.energy_rmse = std::sqrt(sse / n);

    double sae_c = 0.0;
    for (std::size_t i = 0; i < pred.size(); ++i)
        sae_c += std::abs(pred[i] - ref[i] - mean_err);
    e.energy_mae_centered = sae_c / n;

    double base_mean = 0.0;
    for (std::size_t i = 0; i < base.size(); ++i)
        base_mean += base[i] - ref[i];
    base_mean /= n;
    double base_sae_c = 0.0;
    for (std::size_t i = 0; i < base.size(); ++i)
        base_sae_c += std::abs(base[i] - ref[i] - base_mean);
    e.base_energy_mae_centered = base_sae_c / n;

    e.spearman = spearman(pred, ref);
    e.base_spearman = spearman(base, ref);

    if (f_count > 0) {
        e.force_rmse = std::sqrt(f_num / static_cast<double>(f_count));
        e.base_force_rmse = std::sqrt(f_base_num / static_cast<double>(f_count));
    }

    return e;
}

// =================================================================================
// Persistence
// =================================================================================

json DeltaMLModel::toJson() const
{
    json j;
    j["format"] = "curcuma-deltaml-model";
    j["version"] = 1;
    j["base_method"] = m_base_method;
    j["reference_method"] = m_ref_method;
    j["descriptor"] = m_acsf.toJson();
    j["kernel"] = { { "type", "gaussian" }, { "sigma", m_sigma } };

    const int M = static_cast<int>(m_centers.rows());
    const int F = static_cast<int>(m_centers.cols());
    j["n_centers"] = M;
    j["n_features"] = F;

    std::vector<double> flat;
    flat.reserve(static_cast<std::size_t>(M) * F);
    for (int m = 0; m < M; ++m)
        for (int f = 0; f < F; ++f)
            flat.push_back(m_centers(m, f));
    j["centers"] = flat;
    j["center_elements"] = m_center_elements;

    std::vector<double> beta(m_beta.data(), m_beta.data() + m_beta.size());
    j["beta"] = beta;
    std::vector<double> mean(m_feature_mean.data(), m_feature_mean.data() + m_feature_mean.size());
    std::vector<double> scale(m_feature_scale.data(), m_feature_scale.data() + m_feature_scale.size());
    j["feature_mean"] = mean;
    j["feature_scale"] = scale;

    std::vector<int> be_z;
    std::vector<double> be_v;
    for (const auto& kv : m_element_energy) {
        be_z.push_back(kv.first);
        be_v.push_back(kv.second);
    }
    j["baseline_elements"] = be_z;
    j["baseline_energies"] = be_v;

    if (hasAnchor()) {
        std::vector<double> ag;
        ag.reserve(static_cast<std::size_t>(m_anchor_geometry.rows()) * 3);
        for (int i = 0; i < m_anchor_geometry.rows(); ++i)
            for (int g = 0; g < 3; ++g)
                ag.push_back(m_anchor_geometry(i, g));
        j["anchor_elements"] = m_anchor_elements;
        j["anchor_geometry"] = ag;
    }

    return j;
}

bool DeltaMLModel::fromJson(const json& j)
{
    m_trained = false;

    if (j.value("format", std::string()) != "curcuma-deltaml-model") {
        CurcumaLogger::error("DeltaMLModel: file is not a curcuma delta-ML model");
        return false;
    }

    m_base_method = j.value("base_method", std::string("gfnff"));
    m_ref_method = j.value("reference_method", std::string());
    m_acsf = ACSFConfig::fromJson(j.value("descriptor", json::object()));
    m_sigma = j.contains("kernel") ? j["kernel"].value("sigma", 1.0) : 1.0;

    const int M = j.value("n_centers", 0);
    const int F = j.value("n_features", 0);
    if (M <= 0 || F <= 0) {
        CurcumaLogger::error("DeltaMLModel: model file has no kernel centres");
        return false;
    }

    const auto flat = j.value("centers", std::vector<double>{});
    if (static_cast<int>(flat.size()) != M * F) {
        CurcumaLogger::error("DeltaMLModel: centre matrix size does not match n_centers x n_features");
        return false;
    }
    m_centers = Matrix::Zero(M, F);
    for (int m = 0, k = 0; m < M; ++m)
        for (int f = 0; f < F; ++f, ++k)
            m_centers(m, f) = flat[k];

    m_center_elements = j.value("center_elements", std::vector<int>{});
    if (static_cast<int>(m_center_elements.size()) != M) {
        CurcumaLogger::error("DeltaMLModel: centre element list does not match n_centers");
        return false;
    }
    m_element_block.clear();
    for (int m = 0; m < M; ++m) {
        const int z = m_center_elements[m];
        auto it = m_element_block.find(z);
        if (it == m_element_block.end())
            m_element_block[z] = { m, 1 };
        else
            it->second.second += 1;
    }

    const auto beta = j.value("beta", std::vector<double>{});
    if (static_cast<int>(beta.size()) != M) {
        CurcumaLogger::error("DeltaMLModel: weight vector does not match n_centers");
        return false;
    }
    m_beta = Eigen::Map<const Vector>(beta.data(), M);

    const auto mean = j.value("feature_mean", std::vector<double>{});
    const auto scale = j.value("feature_scale", std::vector<double>{});
    if (static_cast<int>(mean.size()) != F || static_cast<int>(scale.size()) != F) {
        CurcumaLogger::error("DeltaMLModel: feature standardisation does not match n_features");
        return false;
    }
    m_feature_mean = Eigen::Map<const Vector>(mean.data(), F);
    m_feature_scale = Eigen::Map<const Vector>(scale.data(), F);

    m_element_energy.clear();
    const auto be_z = j.value("baseline_elements", std::vector<int>{});
    const auto be_v = j.value("baseline_energies", std::vector<double>{});
    for (std::size_t i = 0; i < be_z.size() && i < be_v.size(); ++i)
        m_element_energy[be_z[i]] = be_v[i];

    m_anchor_elements.clear();
    m_anchor_geometry = Geometry();
    if (j.contains("anchor_elements") && j.contains("anchor_geometry")) {
        m_anchor_elements = j["anchor_elements"].get<std::vector<int>>();
        const auto ag = j["anchor_geometry"].get<std::vector<double>>();
        const int na = static_cast<int>(m_anchor_elements.size());
        if (static_cast<int>(ag.size()) == 3 * na) {
            m_anchor_geometry = Geometry::Zero(na, 3);
            for (int i = 0, k = 0; i < na; ++i)
                for (int g = 0; g < 3; ++g, ++k)
                    m_anchor_geometry(i, g) = ag[k];
        } else {
            m_anchor_elements.clear();
        }
    }

    if (m_acsf.featureCount() != F) {
        CurcumaLogger::error_fmt(
            "DeltaMLModel: descriptor yields {} features but the model stores {}",
            m_acsf.featureCount(), F);
        return false;
    }

    m_trained = true;
    return true;
}

bool DeltaMLModel::save(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        CurcumaLogger::error_fmt("DeltaMLModel: cannot write '{}'", filename);
        return false;
    }
    out << toJson().dump(1);
    return true;
}

bool DeltaMLModel::load(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        CurcumaLogger::error_fmt("DeltaMLModel: cannot read '{}'", filename);
        return false;
    }
    json j;
    try {
        in >> j;
    } catch (const std::exception& ex) {
        CurcumaLogger::error_fmt("DeltaMLModel: malformed JSON in '{}': {}", filename, ex.what());
        return false;
    }
    return fromJson(j);
}

} // namespace ml
} // namespace curcuma
