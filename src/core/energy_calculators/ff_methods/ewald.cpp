/*
 * < Ewald summation for point charges in an orthorhombic periodic cell >
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
 * Claude Generated (Nov 2026, B0): direct Ewald sum, formulas in ewald.h.
 */

#include "ewald.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

namespace curcuma {
namespace ewald {

namespace {
constexpr double kPi = 3.14159265358979323846;
using cplx = std::complex<double>;
} // namespace

// Claude Generated (Nov 2026): parameter choice, derivation in ewald.h
EwaldParameters chooseEwaldParameters(const PBCUtils::OrthoCell& cell, double tolerance, double r_cut)
{
    if (!cell.active)
        throw std::invalid_argument("chooseEwaldParameters: the cell is not active");
    if (!(tolerance > 0.0 && tolerance < 1.0))
        throw std::invalid_argument("chooseEwaldParameters: tolerance must lie in (0, 1)");
    EwaldParameters p;
    p.r_cut = (r_cut > 0.0) ? r_cut : 0.5 * cell.minLength();
    const double s = std::sqrt(-std::log(tolerance)); // exp(-s^2) = eps
    p.alpha = s / p.r_cut;                             // exp(-(alpha r_cut)^2) = eps
    p.k_max = 2.0 * p.alpha * s;                       // exp(-k_max^2/(4 alpha^2)) = eps
    return p;
}

// Claude Generated (Nov 2026)
EwaldSum::EwaldSum(const PBCUtils::OrthoCell& cell, double alpha, double r_cut, double k_max)
    : m_cell(cell)
{
    if (!cell.active)
        throw std::invalid_argument("EwaldSum: the cell is not active");
    if (!(alpha > 0.0) || !(k_max > 0.0) || !(r_cut > 0.0))
        throw std::invalid_argument("EwaldSum: alpha, r_cut and k_max must be positive");
    // Only the nearest image of each pair is summed in real space; beyond L_min/2 a second image
    // would be closer than r_cut (minimum-image convention, see pbc_utils.h).
    if (r_cut > 0.5 * cell.minLength() * (1.0 + 1e-12))
        throw std::invalid_argument("EwaldSum: r_cut = " + std::to_string(r_cut)
            + " Bohr exceeds half the shortest box length " + std::to_string(0.5 * cell.minLength()) + " Bohr");

    m_par.alpha = alpha;
    m_par.r_cut = r_cut;
    m_par.k_max = k_max;
    m_volume = cell.L[0] * cell.L[1] * cell.L[2];

    // Reciprocal lattice of a rectangular box: k = 2 pi (n_x/L_x, n_y/L_y, n_z/L_z).
    for (int a = 0; a < 3; ++a)
        m_nmax[a] = static_cast<int>(std::floor(k_max * cell.L[a] / (2.0 * kPi)));

    // Half space: |S(-k)|^2 = |S(k)|^2, so of each pair (n, -n) only the one with the first
    // non-zero index positive is kept and its weight doubled. The order (n_x, n_y outer, n_z inner)
    // lets the evaluation reuse the product of the x and y phase factors along n_z.
    const double k_max2 = k_max * k_max;
    for (int nx = 0; nx <= m_nmax[0]; ++nx) {
        for (int ny = -m_nmax[1]; ny <= m_nmax[1]; ++ny) {
            if (nx == 0 && ny < 0)
                continue;
            for (int nz = -m_nmax[2]; nz <= m_nmax[2]; ++nz) {
                if (nx == 0 && ny == 0 && nz <= 0)
                    continue; // drops k = 0 (tin-foil) and the -n_z partners
                KVector kv;
                kv.n[0] = nx;
                kv.n[1] = ny;
                kv.n[2] = nz;
                for (int a = 0; a < 3; ++a)
                    kv.k[a] = 2.0 * kPi * kv.n[a] / cell.L[a];
                const double k2 = kv.k[0] * kv.k[0] + kv.k[1] * kv.k[1] + kv.k[2] * kv.k[2];
                if (k2 > k_max2)
                    continue;
                // 2 x (1/(2V)) x 4 pi / k^2 x exp(-k^2/(4 alpha^2))
                kv.weight = (4.0 * kPi / m_volume) * std::exp(-k2 / (4.0 * alpha * alpha)) / k2;
                m_kvectors.push_back(kv);
            }
        }
    }
}

// Claude Generated (Nov 2026)
double EwaldSum::realSpaceKernel(double r, double alpha, double& dkernel_dr)
{
    // K(r) = erfc(a r)/r,  dK/dr = -erfc(a r)/r^2 - (2 a/sqrt(pi)) exp(-a^2 r^2)/r
    const double erfc_ar = std::erfc(alpha * r);
    const double inv_r = 1.0 / r;
    dkernel_dr = -erfc_ar * inv_r * inv_r - (2.0 * alpha / std::sqrt(kPi)) * std::exp(-alpha * alpha * r * r) * inv_r;
    return erfc_ar * inv_r;
}

// Claude Generated (Nov 2026, B2): g_min = 1/sqrt(2 max alp); erfc(g_min r) < erfc(XMAX) for r > XMAX / g_min
double gaussianCorrectionRadius(const Eigen::VectorXd& alp)
{
    double amax = 0.0;
    for (Eigen::Index i = 0; i < alp.size(); ++i)
        amax = std::max(amax, alp(i));
    return amax > 0.0 ? GAUSSIAN_CORRECTION_XMAX * std::sqrt(2.0 * amax) : 0.0;
}

// Claude Generated (Nov 2026): real-space pair list with linked cells (Allen & Tildesley, Sec. 5.3, from memory).
// P1 (Nov 2026): the search is split into T ranges of the (cx, cy) cell columns (all-pairs path: of the outer
// atom index); range t writes its own list and the lists are concatenated in range order, so the result is the
// serial list for every T.
namespace {
struct PairSearch {
    const Eigen::MatrixXd& positions;
    const PBCUtils::OrthoCell& cell;
    double rc2;
    int N;
    int nc[3];
    int m = 1;      ///< neighbour cells per direction (cell edge >= r_cut / m)
    bool linked;
    std::vector<std::vector<int>> members;

    PairSearch(const Eigen::MatrixXd& pos, const PBCUtils::OrthoCell& c, double r_cut)
        : positions(pos), cell(c), rc2(r_cut * r_cut), N(static_cast<int>(pos.rows()))
    {
        // Cells of edge >= r_cut/m: every partner within r_cut sits in the (2m+1)^3 surrounding cells. P1 (Nov 2026):
        // m = 2 when every axis has >= 5 such cells (candidate volume 125 (r_cut/2)^3 instead of 27 r_cut^3, and
        // more slabs for the threads), else m = 1 with >= 3 cells (B0), else all pairs (the neighbours would not be distinct).
        for (int mm : { 2, 1 }) {
            for (int a = 0; a < 3; ++a)
                nc[a] = static_cast<int>(std::floor(cell.L[a] * mm / r_cut));
            m = mm;
            if (nc[0] >= 2 * mm + 1 && nc[1] >= 2 * mm + 1 && nc[2] >= 2 * mm + 1)
                break;
        }
        linked = (nc[0] >= 2 * m + 1 && nc[1] >= 2 * m + 1 && nc[2] >= 2 * m + 1);
        if (!linked)
            return;
        // Bin the atoms by their wrapped fractional coordinate.
        members.resize(static_cast<size_t>(nc[0]) * nc[1] * nc[2]);
        int cidx[3];
        for (int i = 0; i < N; ++i) {
            for (int a = 0; a < 3; ++a) {
                double s = positions(i, a) * cell.inv_L[a];
                s -= std::floor(s); // [0, 1)
                cidx[a] = std::min(static_cast<int>(s * nc[a]), nc[a] - 1);
            }
            members[(cidx[0] * nc[1] + cidx[1]) * nc[2] + cidx[2]].push_back(i);
        }
    }
    double micDistance2(int i, int j) const
    {
        double dx = positions(j, 0) - positions(i, 0);
        double dy = positions(j, 1) - positions(i, 1);
        double dz = positions(j, 2) - positions(i, 2);
        PBCUtils::applyMinimumImageOrtho(cell, dx, dy, dz);
        return dx * dx + dy * dy + dz * dz;
    }
    /// outer index: (cx, cy) column of cells, or the atom index i on the all-pairs path
    int outerCount() const { return linked ? nc[0] * nc[1] : N; }
    /// pairs of the outer range [o_begin, o_end), appended in serial order
    void run(int o_begin, int o_end, std::vector<std::pair<int, int>>& pairs) const
    {
        if (!linked) {
            for (int i = o_begin; i < o_end; ++i)
                for (int j = i + 1; j < N; ++j)
                    if (micDistance2(i, j) < rc2)
                        pairs.emplace_back(i, j);
            return;
        }
        // Each unordered pair is seen from both of its cells; i < j keeps exactly one visit.
        for (int o = o_begin; o < o_end; ++o) {
            const int cx = o / nc[1], cy = o % nc[1];
            for (int cz = 0; cz < nc[2]; ++cz) {
                const int c = (cx * nc[1] + cy) * nc[2] + cz;
                for (int ox = -m; ox <= m; ++ox)
                    for (int oy = -m; oy <= m; ++oy)
                        for (int oz = -m; oz <= m; ++oz) {
                            const int nx = (cx + ox + nc[0]) % nc[0];
                            const int ny = (cy + oy + nc[1]) % nc[1];
                            const int nz = (cz + oz + nc[2]) % nc[2];
                            const int nb = (nx * nc[1] + ny) * nc[2] + nz;
                            for (int i : members[c])
                                for (int j : members[nb])
                                    if (i < j && micDistance2(i, j) < rc2)
                                        pairs.emplace_back(i, j);
                        }
            }
        }
    }
};

/// Real-space kernel of one pair: erfc(a r)/r, minus erfc(g r)/r for Gaussian charges with g r < XMAX, and dK/dr.
/// One compiled function for the per-call and the cached path, so both give the same bits (P1): not inlined, so the
/// compiler cannot contract (FMA) or vectorise it differently in the two loops.
#if defined(__GNUC__)
__attribute__((noinline))
#endif
double pairKernel(double r, double alpha, double alp_sum, bool gaussian, double& dK)
{
    static const double two_over_sqrt_pi = 2.0 / std::sqrt(kPi); // the same expression as before P1 (bitwise)
    double K = EwaldSum::realSpaceKernel(r, alpha, dK);
    if (gaussian) {
        // erf(g r)/r = 1/r - erfc(g r)/r: the cloud pair is the point pair minus this term
        const double g = 1.0 / std::sqrt(alp_sum);
        const double x = g * r;
        if (x < GAUSSIAN_CORRECTION_XMAX) {
            const double erfc_g = std::erfc(x);
            const double inv_r = 1.0 / r;
            K -= erfc_g * inv_r;
            dK += erfc_g * inv_r * inv_r + g * two_over_sqrt_pi * std::exp(-x * x) * inv_r;
        }
    }
    return K;
}
} // namespace

std::vector<std::pair<int, int>> EwaldSum::realSpacePairs(const Eigen::MatrixXd& positions) const
{
    const PairSearch search(positions, m_cell, m_par.r_cut);
    std::vector<std::pair<int, int>> pairs;
    search.run(0, search.outerCount(), pairs);
    return pairs;
}

// ---- 1. Real space: sum_{i<j} q_i q_j K_ij(r) over minimum-image pairs with r < r_cut,
// K_ij = erfc(a r)/r, minus erfc(g_ij r)/r for Gaussian charges (B2, g_ij r < XMAX).
// Claude Generated (Nov 2026, B0; range form B2)
void EwaldSum::realSpaceRange(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
    const std::vector<std::pair<int, int>>& pairs, int begin, int end, const Eigen::VectorXd* gaussian_alp,
    bool want_gradient, bool want_potentials, EwaldResult& acc) const
{
    const double alpha = m_par.alpha;
    for (int idx = begin; idx < end; ++idx) {
        const int i = pairs[idx].first, j = pairs[idx].second;
        double d[3] = { positions(j, 0) - positions(i, 0), positions(j, 1) - positions(i, 1),
            positions(j, 2) - positions(i, 2) };
        PBCUtils::applyMinimumImageOrtho(m_cell, d[0], d[1], d[2]);
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        double dK;
        const double K = pairKernel(r, alpha, gaussian_alp ? (*gaussian_alp)(i) + (*gaussian_alp)(j) : 0.0,
            gaussian_alp != nullptr, dK);
        acc.terms.real += charges(i) * charges(j) * K;
        if (want_potentials) {
            acc.potentials(i) += charges(j) * K;
            acc.potentials(j) += charges(i) * K;
        }
        if (want_gradient) {
            // d = r_j - r_i, so dr/dr_j = d/r and dr/dr_i = -d/r
            const double f = charges(i) * charges(j) * dK / r;
            for (int a = 0; a < 3; ++a) {
                acc.gradient(i, a) -= f * d[a];
                acc.gradient(j, a) += f * d[a];
            }
        }
    }
}

// ---- 2a. Phase factors exp(i 2 pi n x_a / L_a) per atom and axis by the recursion
// e^{i(n+1)t} = e^{int} e^{it}; negative n is the complex conjugate. Claude Generated (Nov 2026)
// P1 (Nov 2026): stored as separate cos/sin arrays, n-major; the recursion is the same complex product as before.
EwaldPhaseTable EwaldSum::phaseTable(const Eigen::MatrixXd& positions) const
{
    const int N = static_cast<int>(positions.rows());
    EwaldPhaseTable t;
    t.natoms = N;
    for (int a = 0; a < 3; ++a) {
        t.nmax[a] = m_nmax[a];
        const int nn = m_nmax[a] + 1;
        t.re[a].assign(static_cast<size_t>(N) * nn, 1.0);
        t.im[a].assign(static_cast<size_t>(N) * nn, 0.0);
        for (int i = 0; i < N; ++i) {
            const double ph = 2.0 * kPi * positions(i, a) * m_cell.inv_L[a];
            const double c1 = std::cos(ph), s1 = std::sin(ph);
            double c = 1.0, sn = 0.0;
            for (int n = 1; n < nn; ++n) {
                const double cn = c * c1 - sn * s1;
                const double snn = c * s1 + sn * c1;
                c = cn;
                sn = snn;
                t.re[a][static_cast<size_t>(n) * N + i] = c;
                t.im[a][static_cast<size_t>(n) * N + i] = sn;
            }
        }
    }
    return t;
}

// ---- 2b. Reciprocal space over the k vectors [k_begin, k_end). This loop (N_k x N complex products)
// is the cost centre; the caller spreads the k range over threads. Claude Generated (Nov 2026)
// P1 (Nov 2026): real arithmetic on the structure-of-arrays phase table, atom loops written for the compiler's
// vectoriser (omp simd; the two reductions of S(k) are reassociated by it, so the bits differ from the scalar
// complex loop of B0/B2 at the 1e-16 relative level). e_i = exp(i k.r_i) = exy_i * ez_i with exy refreshed only
// when (n_x, n_y) changes (the k list runs n_z innermost).
void EwaldSum::reciprocalRange(const EwaldPhaseTable& table, const Eigen::VectorXd& charges, int k_begin, int k_end,
    bool want_gradient, bool want_potentials, EwaldResult& acc) const
{
    const int N = table.natoms;
    if (k_begin >= k_end || N == 0)
        return;
    std::vector<double> exr(N), exi(N), er(N), ei(N);
    const double* q = charges.data();
    double* pot = want_potentials ? acc.potentials.data() : nullptr;
    double* gx = want_gradient ? acc.gradient.col(0).data() : nullptr;
    double* gy = want_gradient ? acc.gradient.col(1).data() : nullptr;
    double* gz = want_gradient ? acc.gradient.col(2).data() : nullptr;
    double* pexr = exr.data();
    double* pexi = exi.data();
    double* per = er.data();
    double* pei = ei.data();
    int last_nx = -1, last_ny = 0;
    double e_rec = 0.0;
    for (int kidx = k_begin; kidx < k_end; ++kidx) {
        const KVector& kv = m_kvectors[kidx];
        if (kv.n[0] != last_nx || kv.n[1] != last_ny) {
            const size_t ox = static_cast<size_t>(std::abs(kv.n[0])) * N, oy = static_cast<size_t>(std::abs(kv.n[1])) * N;
            const double* xr = table.re[0].data() + ox;
            const double* xi = table.im[0].data() + ox;
            const double* yr = table.re[1].data() + oy;
            const double* yi = table.im[1].data() + oy;
            const double sx = kv.n[0] >= 0 ? 1.0 : -1.0, sy = kv.n[1] >= 0 ? 1.0 : -1.0;
#pragma omp simd
            for (int i = 0; i < N; ++i) {
                const double b = sx * xi[i], d = sy * yi[i];
                pexr[i] = xr[i] * yr[i] - b * d;
                pexi[i] = xr[i] * d + b * yr[i];
            }
            last_nx = kv.n[0];
            last_ny = kv.n[1];
        }
        const size_t oz = static_cast<size_t>(std::abs(kv.n[2])) * N;
        const double* zr = table.re[2].data() + oz;
        const double* zi = table.im[2].data() + oz;
        const double sz = kv.n[2] >= 0 ? 1.0 : -1.0;
        double Sr = 0.0, Si = 0.0; // structure factor S(k) = sum_j q_j exp(i k.r_j)
#pragma omp simd reduction(+ : Sr, Si)
        for (int i = 0; i < N; ++i) {
            const double zim = sz * zi[i];
            const double a = pexr[i] * zr[i] - pexi[i] * zim;
            const double b = pexr[i] * zim + pexi[i] * zr[i];
            per[i] = a;
            pei[i] = b;
            Sr += q[i] * a;
            Si += q[i] * b;
        }
        e_rec += kv.weight * (Sr * Sr + Si * Si); // weight already holds the -k partner
        if (pot) {
            // d|S|^2/dq_i = 2 Re(conj(e_i) S) = 2 (e_r S_r + e_i S_i)
            const double w2 = 2.0 * kv.weight;
#pragma omp simd
            for (int i = 0; i < N; ++i)
                pot[i] += w2 * (per[i] * Sr + pei[i] * Si);
        }
        if (gx) {
            // d|S|^2/dr_i = 2 q_i k Im(conj(e_i) S) = 2 q_i k (e_r S_i - e_i S_r)
            const double w2 = 2.0 * kv.weight;
            const double kx = kv.k[0], ky = kv.k[1], kz = kv.k[2];
#pragma omp simd
            for (int i = 0; i < N; ++i) {
                const double g = w2 * q[i] * (per[i] * Si - pei[i] * Sr);
                gx[i] += g * kx;
                gy[i] += g * ky;
                gz[i] += g * kz;
            }
        }
    }
    acc.terms.reciprocal += e_rec;
}

// ---- P1 (Nov 2026): shared geometry data, Claude Generated
size_t EwaldGeometryCache::bytes() const
{
    size_t b = pairs.capacity() * sizeof(std::pair<int, int>) + (kernel.capacity() + dkernel.capacity()) * sizeof(double);
    for (int a = 0; a < 3; ++a)
        b += (table.re[a].capacity() + table.im[a].capacity()) * sizeof(double);
    return b;
}

bool EwaldSum::cacheValidFor(const EwaldGeometryCache& cache, const Eigen::MatrixXd& positions, const Eigen::VectorXd* gaussian_alp) const
{
    if (!cache.valid || cache.alpha != m_par.alpha || cache.r_cut != m_par.r_cut || cache.k_max != m_par.k_max)
        return false;
    if (cache.pos.rows() != positions.rows() || cache.pos.cols() != positions.cols())
        return false;
    if (gaussian_alp ? (cache.alp.size() != gaussian_alp->size()) : (cache.alp.size() != 0))
        return false;
    for (Eigen::Index c = 0; c < positions.cols(); ++c)
        for (Eigen::Index r = 0; r < positions.rows(); ++r)
            if (cache.pos(r, c) != positions(r, c))
                return false;
    if (gaussian_alp)
        for (Eigen::Index i = 0; i < gaussian_alp->size(); ++i)
            if (cache.alp(i) != (*gaussian_alp)(i))
                return false;
    return true;
}

bool EwaldSum::buildGeometryCache(const Eigen::MatrixXd& positions, const Eigen::VectorXd* gaussian_alp, EwaldGeometryCache& cache,
    size_t max_kernel_bytes, int T, const EwaldParallelFor& run) const
{
    if (cacheValidFor(cache, positions, gaussian_alp)) {
        ++cache.reuses;
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    if (T < 1 || !run)
        T = 1;
    auto parallel = [&](const std::function<void(int)>& f) {
        if (T > 1)
            run(T, f);
        else
            f(0);
    };
    cache.valid = false;
    cache.pos = positions;
    if (gaussian_alp)
        cache.alp = *gaussian_alp;
    else
        cache.alp.resize(0);
    cache.alpha = m_par.alpha;
    cache.r_cut = m_par.r_cut;
    cache.k_max = m_par.k_max;

    // Pair list: slabs in parallel, concatenated in slab order (= the serial list)
    const PairSearch search(cache.pos, m_cell, m_par.r_cut);
    const int outer = search.outerCount();
    std::vector<std::vector<std::pair<int, int>>> part(T);
    parallel([&](int t) {
        const int b = static_cast<int>(static_cast<long long>(outer) * t / T);
        const int e = static_cast<int>(static_cast<long long>(outer) * (t + 1) / T);
        search.run(b, e, part[t]);
    });
    size_t npairs = 0;
    for (const auto& p : part)
        npairs += p.size();
    cache.pairs.clear();
    cache.pairs.reserve(npairs);
    for (const auto& p : part)
        cache.pairs.insert(cache.pairs.end(), p.begin(), p.end());

    const auto t1 = std::chrono::steady_clock::now();
    cache.pairs_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    // Kernel values once per geometry, unless over the cap (then realSpaceRangeCached evaluates them per call)
    const size_t kernel_bytes = 2 * sizeof(double) * npairs;
    if (kernel_bytes <= max_kernel_bytes && npairs > 0) {
        cache.kernel.resize(npairs);
        cache.dkernel.resize(npairs);
        const double alpha = m_par.alpha;
        const Eigen::MatrixXd& P = cache.pos;
        const Eigen::VectorXd* alp = gaussian_alp ? &cache.alp : nullptr;
        parallel([&](int t) {
            const size_t b = npairs * static_cast<size_t>(t) / T, e = npairs * static_cast<size_t>(t + 1) / T;
            for (size_t idx = b; idx < e; ++idx) {
                const int i = cache.pairs[idx].first, j = cache.pairs[idx].second;
                double d[3] = { P(j, 0) - P(i, 0), P(j, 1) - P(i, 1), P(j, 2) - P(i, 2) };
                PBCUtils::applyMinimumImageOrtho(m_cell, d[0], d[1], d[2]);
                const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                double dK;
                cache.kernel[idx] = pairKernel(r, alpha, alp ? (*alp)(i) + (*alp)(j) : 0.0, alp != nullptr, dK);
                cache.dkernel[idx] = dK;
            }
        });
    } else {
        std::vector<double>().swap(cache.kernel);
        std::vector<double>().swap(cache.dkernel);
    }
    const auto t2 = std::chrono::steady_clock::now();
    cache.kernel_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    cache.table = phaseTable(cache.pos);
    cache.table_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t2).count();
    cache.valid = true;
    ++cache.builds;
    cache.build_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

void EwaldSum::realSpaceRangeCached(const EwaldGeometryCache& cache, const Eigen::VectorXd& charges, int begin, int end,
    bool want_gradient, bool want_potentials, EwaldResult& acc) const
{
    if (!cache.hasKernel()) {
        realSpaceRange(cache.pos, charges, cache.pairs, begin, end, cache.alp.size() ? &cache.alp : nullptr,
            want_gradient, want_potentials, acc);
        return;
    }
    const double* q = charges.data();
    const double* K = cache.kernel.data();
    const auto* pr = cache.pairs.data();
    double e_real = 0.0;
    if (!want_gradient) {
        // EEQ matrix-vector product and energy: a gather over the stored kernel, no transcendental function
        double* pot = want_potentials ? acc.potentials.data() : nullptr;
        for (int idx = begin; idx < end; ++idx) {
            const int i = pr[idx].first, j = pr[idx].second;
            e_real += q[i] * q[j] * K[idx];
            if (pot) {
                pot[i] += q[j] * K[idx];
                pot[j] += q[i] * K[idx];
            }
        }
        acc.terms.real += e_real;
        return;
    }
    const Eigen::MatrixXd& P = cache.pos;
    const double* dK = cache.dkernel.data();
    for (int idx = begin; idx < end; ++idx) {
        const int i = pr[idx].first, j = pr[idx].second;
        double d[3] = { P(j, 0) - P(i, 0), P(j, 1) - P(i, 1), P(j, 2) - P(i, 2) };
        PBCUtils::applyMinimumImageOrtho(m_cell, d[0], d[1], d[2]);
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        acc.terms.real += q[i] * q[j] * K[idx];
        if (want_potentials) {
            acc.potentials(i) += q[j] * K[idx];
            acc.potentials(j) += q[i] * K[idx];
        }
        // d = r_j - r_i, so dr/dr_j = d/r and dr/dr_i = -d/r
        const double f = q[i] * q[j] * dK[idx] / r;
        for (int a = 0; a < 3; ++a) {
            acc.gradient(i, a) -= f * d[a];
            acc.gradient(j, a) += f * d[a];
        }
    }
}

// ---- Atomic virial (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED N1)
// Strain eps_ab (symmetric) at fixed fractional coordinates: r -> (1 + eps) r for every pair vector, V -> V (1 + tr eps),
// and k -> (1 - eps) k (k = 2 pi n / L stays conjugate to the scaled lattice, so k.r is invariant). Xi_ab = -dE/d(eps_ab).
//  Real space: E_pair = q_i q_j K(r), dr/d(eps_ab) = d_a d_b / r, so Xi_ab -= q_i q_j K'(r) d_a d_b / r.
//  Reciprocal: E_k = w_k |S(k)|^2 with w_k = (4 pi/V) exp(-k^2/(4 a^2)) / k^2 and |S|^2 strain-independent;
//    d(k^2)/d(eps_ab) = -2 k_a k_b and dV/d(eps_ab) = V delta_ab give
//    d ln w_k / d(eps_ab) = -delta_ab + 2 k_a k_b (1/k^2 + 1/(4 a^2)),  Xi_ab += E_k [delta_ab - 2 k_a k_b (1/k^2 + 1/(4 a^2))].
//    Trace: Xi += E_k [1 - k^2/(2 a^2)] (the B0 note "E_k [1 - 2 k_a^2 (1/k^2 + 1/(4 a^2))]" per axis is this tensor's
//    diagonal). Alpha is held fixed (the energy does not depend on it beyond the tolerance).
//  Self: -(a/sqrt pi) sum q^2, no strain dependence. Background -pi Q^2/(2 a^2 V) and Gaussian background
//    (pi Q/V) sum q alp: both ~ 1/V, so Xi_ab = E delta_ab (backgroundVirial).
void EwaldSum::realSpaceVirialRange(const EwaldGeometryCache& cache, const Eigen::VectorXd& charges, int begin, int end,
    Eigen::Matrix3d& virial) const
{
    const Eigen::MatrixXd& P = cache.pos;
    const bool stored = cache.hasKernel();
    const Eigen::VectorXd* alp = cache.alp.size() ? &cache.alp : nullptr;
    const double alpha = m_par.alpha;
    Eigen::Matrix3d xi = Eigen::Matrix3d::Zero();
    for (int idx = begin; idx < end; ++idx) {
        const int i = cache.pairs[idx].first, j = cache.pairs[idx].second;
        double d[3] = { P(j, 0) - P(i, 0), P(j, 1) - P(i, 1), P(j, 2) - P(i, 2) };
        PBCUtils::applyMinimumImageOrtho(m_cell, d[0], d[1], d[2]);
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        double dK;
        if (stored) {
            dK = cache.dkernel[idx];
        } else {
            pairKernel(r, alpha, alp ? (*alp)(i) + (*alp)(j) : 0.0, alp != nullptr, dK);
        }
        const double f = -charges(i) * charges(j) * dK / r;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                xi(a, b) += f * d[a] * d[b];
    }
    virial += xi;
}

void EwaldSum::reciprocalVirialRange(const EwaldPhaseTable& table, const Eigen::VectorXd& charges, int k_begin, int k_end,
    Eigen::Matrix3d& virial) const
{
    const int N = table.natoms;
    if (k_begin >= k_end || N == 0)
        return;
    std::vector<double> exr(N), exi(N);
    const double* q = charges.data();
    int last_nx = -1, last_ny = 0;
    const double inv_4a2 = 1.0 / (4.0 * m_par.alpha * m_par.alpha);
    Eigen::Matrix3d xi = Eigen::Matrix3d::Zero();
    for (int kidx = k_begin; kidx < k_end; ++kidx) {
        const KVector& kv = m_kvectors[kidx];
        if (kv.n[0] != last_nx || kv.n[1] != last_ny) {
            const size_t ox = static_cast<size_t>(std::abs(kv.n[0])) * N, oy = static_cast<size_t>(std::abs(kv.n[1])) * N;
            const double sx = kv.n[0] >= 0 ? 1.0 : -1.0, sy = kv.n[1] >= 0 ? 1.0 : -1.0;
            for (int i = 0; i < N; ++i) {
                const double xr = table.re[0][ox + i], b = sx * table.im[0][ox + i];
                const double yr = table.re[1][oy + i], dd = sy * table.im[1][oy + i];
                exr[i] = xr * yr - b * dd;
                exi[i] = xr * dd + b * yr;
            }
            last_nx = kv.n[0];
            last_ny = kv.n[1];
        }
        const size_t oz = static_cast<size_t>(std::abs(kv.n[2])) * N;
        const double sz = kv.n[2] >= 0 ? 1.0 : -1.0;
        double Sr = 0.0, Si = 0.0;
        for (int i = 0; i < N; ++i) {
            const double zr = table.re[2][oz + i], zim = sz * table.im[2][oz + i];
            Sr += q[i] * (exr[i] * zr - exi[i] * zim);
            Si += q[i] * (exr[i] * zim + exi[i] * zr);
        }
        const double ek = kv.weight * (Sr * Sr + Si * Si);
        const double k2 = kv.k[0] * kv.k[0] + kv.k[1] * kv.k[1] + kv.k[2] * kv.k[2];
        const double c = 2.0 * (1.0 / k2 + inv_4a2);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                xi(a, b) += ek * ((a == b ? 1.0 : 0.0) - c * kv.k[a] * kv.k[b]);
    }
    virial += xi;
}

// ---- 3. Self energy: each Gaussian screening cloud interacts with its own point charge in step 2.
// ---- 4. Neutralising background for a charged cell (Q != 0); B2: plus the Gaussian background
//         (pi Q/V) sum_i q_i alp_i (derivation in ewald.h). Claude Generated (Nov 2026)
void EwaldSum::selfAndBackground(const Eigen::VectorXd& charges, const Eigen::VectorXd* gaussian_alp,
    bool want_potentials, EwaldResult& acc) const
{
    const int N = static_cast<int>(charges.size());
    const double alpha = m_par.alpha;
    const double self_pref = alpha / std::sqrt(kPi);
    const double Q = charges.sum();
    const double bg_pref = kPi / (alpha * alpha * m_volume);
    acc.terms.self += -self_pref * charges.squaredNorm();
    acc.terms.background += -0.5 * bg_pref * Q * Q;
    double q_alp = 0.0; // sum_i q_i alp_i
    if (gaussian_alp) {
        q_alp = charges.dot(*gaussian_alp);
        acc.terms.gaussian_background += (kPi / m_volume) * Q * q_alp;
    }
    if (want_potentials) {
        for (int i = 0; i < N; ++i) {
            acc.potentials(i) += -2.0 * self_pref * charges(i) - bg_pref * Q;
            if (gaussian_alp) // d/dq_i of (pi/V) Q sum_j q_j alp_j
                acc.potentials(i) += (kPi / m_volume) * (q_alp + Q * (*gaussian_alp)(i));
        }
    }
}

// Claude Generated (Nov 2026): one pass over pairs and k vectors, composed of the range functions
EwaldResult EwaldSum::evaluate(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
    bool want_gradient, bool want_potentials, const Eigen::VectorXd* gaussian_alp) const
{
    const int N = static_cast<int>(positions.rows());
    if (positions.cols() != 3 || charges.size() != N)
        throw std::invalid_argument("EwaldSum: positions must be N x 3 and charges of length N");
    if (gaussian_alp && gaussian_alp->size() != N)
        throw std::invalid_argument("EwaldSum: gaussian_alp must have length N");

    EwaldResult res;
    if (want_gradient)
        res.gradient = Eigen::MatrixXd::Zero(N, 3);
    if (want_potentials)
        res.potentials = Eigen::VectorXd::Zero(N);

    const std::vector<std::pair<int, int>> pairs = realSpacePairs(positions);
    realSpaceRange(positions, charges, pairs, 0, static_cast<int>(pairs.size()), gaussian_alp,
        want_gradient, want_potentials, res);
    const EwaldPhaseTable table = phaseTable(positions);
    reciprocalRange(table, charges, 0, numberOfKVectors(), want_gradient, want_potentials, res);
    selfAndBackground(charges, gaussian_alp, want_potentials, res);
    return res;
}

double EwaldSum::energy(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
    const Eigen::VectorXd* gaussian_alp) const
{
    return evaluate(positions, charges, false, false, gaussian_alp).terms.total();
}

Eigen::MatrixXd EwaldSum::gradient(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
    const Eigen::VectorXd* gaussian_alp) const
{
    return evaluate(positions, charges, true, false, gaussian_alp).gradient;
}

Eigen::VectorXd EwaldSum::potentials(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
    const Eigen::VectorXd* gaussian_alp) const
{
    return evaluate(positions, charges, false, true, gaussian_alp).potentials;
}

// Claude Generated (Nov 2026): dense Ewald matrix, formula in ewald.h
Eigen::MatrixXd EwaldSum::matrix(const Eigen::MatrixXd& positions, const Eigen::VectorXd* gaussian_alp) const
{
    const int N = static_cast<int>(positions.rows());
    const double alpha = m_par.alpha;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);

    for (const auto& pr : realSpacePairs(positions)) {
        const int i = pr.first, j = pr.second;
        double dx = positions(j, 0) - positions(i, 0), dy = positions(j, 1) - positions(i, 1),
               dz = positions(j, 2) - positions(i, 2);
        PBCUtils::applyMinimumImageOrtho(m_cell, dx, dy, dz);
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        double dK;
        double K = realSpaceKernel(r, alpha, dK);
        if (gaussian_alp) {
            const double x = r / std::sqrt((*gaussian_alp)(i) + (*gaussian_alp)(j));
            if (x < GAUSSIAN_CORRECTION_XMAX)
                K -= std::erfc(x) / r;
        }
        A(i, j) += K;
        A(j, i) += K;
    }
    if (gaussian_alp) // (pi Q/V) sum_i q_i alp_i = 1/2 sum_ij q_i q_j (pi/V)(alp_i + alp_j)
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                A(i, j) += (kPi / m_volume) * ((*gaussian_alp)(i) + (*gaussian_alp)(j));

    // 2 w cos(k.(r_i - r_j)) = 2 w (cos_i cos_j + sin_i sin_j), the diagonal included. B3 (Nov 2026): the k
    // vectors are taken in blocks, C(i, k) = sqrt(2 w_k) cos(k.r_i), S likewise, and the block enters as the
    // symmetric rank update A += C C^T + S S^T (one GEMM-class product per block instead of N_k rank-2
    // updates); phases from the recursion table, not cos/sin per (k, atom). Lower triangle only, mirrored at the end.
    const EwaldPhaseTable tab = phaseTable(positions);
    auto phase = [&](int a, int i, int n) {
        const size_t o = static_cast<size_t>(std::abs(n)) * N + i;
        const cplx v(tab.re[a][o], tab.im[a][o]);
        return n >= 0 ? v : std::conj(v);
    };
    const int Nk = numberOfKVectors();
    const int KB = std::max(1, std::min(Nk, static_cast<int>((16u << 20) / (2 * sizeof(double) * std::max(1, N)))));
    Eigen::MatrixXd C(N, KB), S(N, KB);
    for (int k0 = 0; k0 < Nk; k0 += KB) {
        const int nb = std::min(KB, Nk - k0);
        for (int kk = 0; kk < nb; ++kk) {
            const KVector& kv = m_kvectors[k0 + kk];
            const double sw = std::sqrt(2.0 * kv.weight);
            for (int i = 0; i < N; ++i) {
                const cplx e = phase(0, i, kv.n[0]) * phase(1, i, kv.n[1]) * phase(2, i, kv.n[2]);
                C(i, kk) = sw * e.real();
                S(i, kk) = sw * e.imag();
            }
        }
        A.selfadjointView<Eigen::Lower>().rankUpdate(C.leftCols(nb));
        A.selfadjointView<Eigen::Lower>().rankUpdate(S.leftCols(nb));
    }
    A.triangularView<Eigen::StrictlyUpper>() = A.transpose();

    A.diagonal().array() -= 2.0 * alpha / std::sqrt(kPi);
    A.array() -= kPi / (alpha * alpha * m_volume);
    return A;
}

// Claude Generated (Nov 2026, B3): the diagonal of the reciprocal part is the same for every atom
double EwaldSum::reciprocalDiagonal() const
{
    double d = 0.0;
    for (const KVector& kv : m_kvectors)
        d += 2.0 * kv.weight;
    return d;
}

// Claude Generated (Nov 2026, B3): A_ii of matrix() without building it (Jacobi preconditioner)
Eigen::VectorXd EwaldSum::diagonal(int natoms, const Eigen::VectorXd* gaussian_alp) const
{
    const double alpha = m_par.alpha;
    const double base = reciprocalDiagonal() - 2.0 * alpha / std::sqrt(kPi) - kPi / (alpha * alpha * m_volume);
    Eigen::VectorXd d = Eigen::VectorXd::Constant(natoms, base);
    if (gaussian_alp)
        d += (2.0 * kPi / m_volume) * (*gaussian_alp);
    return d;
}

} // namespace ewald
} // namespace curcuma
