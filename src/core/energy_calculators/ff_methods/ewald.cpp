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
    // B2 hook: a Gaussian-charge correction -erfc(g_ij r)/r (and its derivative) belongs here.
    const double erfc_ar = std::erfc(alpha * r);
    const double inv_r = 1.0 / r;
    dkernel_dr = -erfc_ar * inv_r * inv_r - (2.0 * alpha / std::sqrt(kPi)) * std::exp(-alpha * alpha * r * r) * inv_r;
    return erfc_ar * inv_r;
}

// Claude Generated (Nov 2026): real-space pair list with linked cells (Allen & Tildesley, Sec. 5.3, from memory)
std::vector<std::pair<int, int>> EwaldSum::realSpacePairs(const Eigen::MatrixXd& positions) const
{
    const int N = static_cast<int>(positions.rows());
    const double rc2 = m_par.r_cut * m_par.r_cut;
    std::vector<std::pair<int, int>> pairs;

    auto micDistance2 = [&](int i, int j) {
        double dx = positions(j, 0) - positions(i, 0);
        double dy = positions(j, 1) - positions(i, 1);
        double dz = positions(j, 2) - positions(i, 2);
        PBCUtils::applyMinimumImageOrtho(m_cell, dx, dy, dz);
        return dx * dx + dy * dy + dz * dz;
    };

    // Cells of edge >= r_cut: every partner within r_cut sits in the 27 surrounding cells. With
    // fewer than 3 cells along an axis the 27 neighbours are not distinct, so fall back to all pairs.
    int nc[3];
    for (int a = 0; a < 3; ++a)
        nc[a] = static_cast<int>(std::floor(m_cell.L[a] / m_par.r_cut));
    if (nc[0] < 3 || nc[1] < 3 || nc[2] < 3) {
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                if (micDistance2(i, j) < rc2)
                    pairs.emplace_back(i, j);
        return pairs;
    }

    // Bin the atoms by their wrapped fractional coordinate.
    const int ncell = nc[0] * nc[1] * nc[2];
    std::vector<std::vector<int>> members(ncell);
    std::vector<int> cell_of(N);
    int cidx[3];
    for (int i = 0; i < N; ++i) {
        for (int a = 0; a < 3; ++a) {
            double s = positions(i, a) * m_cell.inv_L[a];
            s -= std::floor(s); // [0, 1)
            cidx[a] = std::min(static_cast<int>(s * nc[a]), nc[a] - 1);
        }
        cell_of[i] = (cidx[0] * nc[1] + cidx[1]) * nc[2] + cidx[2];
        members[cell_of[i]].push_back(i);
    }

    // Each unordered pair is seen from both of its cells; i < j keeps exactly one visit.
    for (int cx = 0; cx < nc[0]; ++cx)
        for (int cy = 0; cy < nc[1]; ++cy)
            for (int cz = 0; cz < nc[2]; ++cz) {
                const int c = (cx * nc[1] + cy) * nc[2] + cz;
                for (int ox = -1; ox <= 1; ++ox)
                    for (int oy = -1; oy <= 1; ++oy)
                        for (int oz = -1; oz <= 1; ++oz) {
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
    return pairs;
}

// Claude Generated (Nov 2026): one pass over pairs and k vectors
EwaldResult EwaldSum::evaluate(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
    bool want_gradient, bool want_potentials) const
{
    const int N = static_cast<int>(positions.rows());
    if (positions.cols() != 3 || charges.size() != N)
        throw std::invalid_argument("EwaldSum: positions must be N x 3 and charges of length N");
    const double alpha = m_par.alpha;

    EwaldResult res;
    if (want_gradient)
        res.gradient = Eigen::MatrixXd::Zero(N, 3);
    if (want_potentials)
        res.potentials = Eigen::VectorXd::Zero(N);

    // ---- 1. Real space: sum_{i<j} q_i q_j erfc(a r)/r over minimum-image pairs with r < r_cut.
    for (const auto& pr : realSpacePairs(positions)) {
        const int i = pr.first, j = pr.second;
        double d[3] = { positions(j, 0) - positions(i, 0), positions(j, 1) - positions(i, 1),
            positions(j, 2) - positions(i, 2) };
        PBCUtils::applyMinimumImageOrtho(m_cell, d[0], d[1], d[2]);
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        double dK;
        const double K = realSpaceKernel(r, alpha, dK);
        res.terms.real += charges(i) * charges(j) * K;
        if (want_potentials) {
            res.potentials(i) += charges(j) * K;
            res.potentials(j) += charges(i) * K;
        }
        if (want_gradient) {
            // d = r_j - r_i, so dr/dr_j = d/r and dr/dr_i = -d/r
            const double f = charges(i) * charges(j) * dK / r;
            for (int a = 0; a < 3; ++a) {
                res.gradient(i, a) -= f * d[a];
                res.gradient(j, a) += f * d[a];
            }
        }
    }

    // ---- 2. Reciprocal space. Phase factors exp(i 2 pi n x_a / L_a) per atom and axis by the
    // recursion e^{i(n+1)t} = e^{int} e^{it}; negative n is the complex conjugate.
    std::vector<cplx> table[3];
    for (int a = 0; a < 3; ++a) {
        const int nn = m_nmax[a] + 1;
        table[a].assign(static_cast<size_t>(N) * nn, cplx(1.0, 0.0));
        for (int i = 0; i < N; ++i) {
            const double t = 2.0 * kPi * positions(i, a) * m_cell.inv_L[a];
            const cplx e1(std::cos(t), std::sin(t));
            cplx* row = &table[a][static_cast<size_t>(i) * nn];
            for (int n = 1; n < nn; ++n)
                row[n] = row[n - 1] * e1;
        }
    }
    auto phase = [&](int a, int i, int n) {
        const cplx v = table[a][static_cast<size_t>(i) * (m_nmax[a] + 1) + std::abs(n)];
        return n >= 0 ? v : std::conj(v);
    };

    // This k loop is the cost centre (N_k x N complex products) and the place to parallelise.
    std::vector<cplx> exy(N), eikr(N);
    int last_nx = -1, last_ny = 0;
    for (const KVector& kv : m_kvectors) {
        if (kv.n[0] != last_nx || kv.n[1] != last_ny) {
            for (int i = 0; i < N; ++i)
                exy[i] = phase(0, i, kv.n[0]) * phase(1, i, kv.n[1]);
            last_nx = kv.n[0];
            last_ny = kv.n[1];
        }
        cplx S(0.0, 0.0); // structure factor S(k) = sum_j q_j exp(i k.r_j)
        for (int i = 0; i < N; ++i) {
            eikr[i] = exy[i] * phase(2, i, kv.n[2]);
            S += charges(i) * eikr[i];
        }
        res.terms.reciprocal += kv.weight * std::norm(S); // weight already holds the -k partner
        if (!want_potentials && !want_gradient)
            continue;
        for (int i = 0; i < N; ++i) {
            // conj(e_i) S = [cos_i C + sin_i S_s] + i [cos_i S_s - sin_i C]
            const cplx z = std::conj(eikr[i]) * S;
            if (want_potentials)
                res.potentials(i) += 2.0 * kv.weight * z.real(); // d|S|^2/dq_i = 2 Re(conj(e_i) S)
            if (want_gradient) {
                const double g = 2.0 * kv.weight * charges(i) * z.imag(); // d|S|^2/dr_i = 2 q_i k Im(conj(e_i) S)
                for (int a = 0; a < 3; ++a)
                    res.gradient(i, a) += g * kv.k[a];
            }
        }
    }

    // ---- 3. Self energy: each Gaussian screening cloud interacts with its own point charge in step 2.
    const double self_pref = alpha / std::sqrt(kPi);
    // ---- 4. Neutralising background for a charged cell (Q != 0).
    const double Q = charges.sum();
    const double bg_pref = kPi / (alpha * alpha * m_volume);
    res.terms.self = -self_pref * charges.squaredNorm();
    res.terms.background = -0.5 * bg_pref * Q * Q;
    if (want_potentials) {
        for (int i = 0; i < N; ++i)
            res.potentials(i) += -2.0 * self_pref * charges(i) - bg_pref * Q;
    }
    return res;
}

double EwaldSum::energy(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges) const
{
    return evaluate(positions, charges, false, false).terms.total();
}

Eigen::MatrixXd EwaldSum::gradient(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges) const
{
    return evaluate(positions, charges, true, false).gradient;
}

Eigen::VectorXd EwaldSum::potentials(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges) const
{
    return evaluate(positions, charges, false, true).potentials;
}

// Claude Generated (Nov 2026): dense Ewald matrix, formula in ewald.h
Eigen::MatrixXd EwaldSum::matrix(const Eigen::MatrixXd& positions) const
{
    const int N = static_cast<int>(positions.rows());
    const double alpha = m_par.alpha;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);

    for (const auto& pr : realSpacePairs(positions)) {
        const int i = pr.first, j = pr.second;
        double dx = positions(j, 0) - positions(i, 0), dy = positions(j, 1) - positions(i, 1),
               dz = positions(j, 2) - positions(i, 2);
        PBCUtils::applyMinimumImageOrtho(m_cell, dx, dy, dz);
        double dK;
        const double K = realSpaceKernel(std::sqrt(dx * dx + dy * dy + dz * dz), alpha, dK);
        A(i, j) += K;
        A(j, i) += K;
    }

    // 2 w cos(k.(r_i - r_j)) = 2 w (cos_i cos_j + sin_i sin_j), the diagonal included
    Eigen::VectorXd c(N), s(N);
    for (const KVector& kv : m_kvectors) {
        for (int i = 0; i < N; ++i) {
            const double t = kv.k[0] * positions(i, 0) + kv.k[1] * positions(i, 1) + kv.k[2] * positions(i, 2);
            c(i) = std::cos(t);
            s(i) = std::sin(t);
        }
        A.noalias() += (2.0 * kv.weight) * (c * c.transpose() + s * s.transpose());
    }

    A.diagonal().array() -= 2.0 * alpha / std::sqrt(kPi);
    A.array() -= kPi / (alpha * alpha * m_volume);
    return A;
}

} // namespace ewald
} // namespace curcuma
