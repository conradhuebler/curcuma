/*
 * <Native QM: XC kernel response F_xc[D] and the P2-directional XC gradient of the double hybrids (WP12 stage 2)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas: qm_xc_fxc.h (buildFxc, buildXCResponseGradient) and docs/DFT_ROADMAP/WP12_doppelhybride.md, stage 2.
 * CPU only. Blocks as the meta-GGA paths of qm_xc.cpp (AO screen blocks with the kept functions, else 128-point
 * blocks with all functions); a fixed number of chunks (32) with own accumulators added in chunk order, so the
 * result does not depend on the thread count.
 *
 * Claude Generated (Oct 2026, WP12 stage 2). AI-generated, machine-tested; human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_xc_fxc.h"

#include <algorithm>
#include <chrono>

namespace qmxc {

namespace {

/// P restricted to the kept functions (as gatherP in qm_xc.cpp).
Eigen::MatrixXd gatherSub(const Eigen::MatrixXd& P, const std::vector<int>& cols)
{
    const int ns = static_cast<int>(cols.size());
    Eigen::MatrixXd S(ns, ns);
    for (int j = 0; j < ns; ++j)
        for (int i = 0; i < ns; ++i) S(i, j) = P(cols[i], cols[j]);
    return S;
}

/// Density quantities of one symmetric AO matrix on a block: rho, grad rho (3 columns), tau (if mgga), phi X, d_j phi X.
struct BlockDensity {
    Eigen::VectorXd rho, tau;
    Eigen::MatrixXd grad;   // count x 3
    Eigen::MatrixXd Xphi;   // phi X_S (count x ns)
    Eigen::MatrixXd dX[3];  // d_j phi X_S
};

void blockDensity(const Eigen::MatrixXd& phi, const Eigen::MatrixXd* d, bool gga, bool mgga, const Eigen::MatrixXd& Xs,
                  BlockDensity& out)
{
    const int nb = static_cast<int>(phi.rows());
    out.Xphi.noalias() = phi * Xs;
    out.rho.resize(nb);
    for (int q = 0; q < nb; ++q) out.rho(q) = phi.row(q).dot(out.Xphi.row(q));
    out.grad = Eigen::MatrixXd::Zero(nb, 3);
    if (gga)
        for (int q = 0; q < nb; ++q)
            for (int k = 0; k < 3; ++k) out.grad(q, k) = 2.0 * d[k].row(q).dot(out.Xphi.row(q));
    out.tau = Eigen::VectorXd::Zero(nb);
    if (mgga)
        for (int k = 0; k < 3; ++k) {
            out.dX[k].noalias() = d[k] * Xs;
            for (int q = 0; q < nb; ++q) out.tau(q) += 0.5 * d[k].row(q).dot(out.dX[k].row(q));
        }
}

struct BlockLoop {
    int nblocks = 0, block = 128, npts = 0, n = 0;
    bool screened = false;
};

/// AO values of block b (deriv as requested) and the active-function list `cols`.
void evaluateBlock(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const qmgrid::BlockScreen* screen,
                   const BlockLoop& L, int b, int deriv, qmgrid::AOBlock& blk, std::vector<int>& cols, int& begin,
                   int& count, bool& empty)
{
    empty = false;
    if (L.screened) {
        begin = screen->begin[b];
        count = screen->count[b];
        const int off = screen->shell_off[b], nsh = screen->shell_off[b + 1] - off;
        if (nsh == 0) {
            empty = true;
            return;
        }
        ao.evaluateShells(grid.points, begin, count, deriv, screen->shells.data() + off, nsh, blk, cols);
    } else {
        begin = b * L.block;
        count = std::min(L.block, L.npts - begin);
        ao.evaluate(grid.points, begin, count, deriv, blk);
        cols.resize(L.n);
        for (int m = 0; m < L.n; ++m) cols[m] = m;
    }
}

double rhoMin(const XCFunctional& f, double rho_threshold)
{
    return kern::effectiveThreshold(f.needsGradient(), rho_threshold);
}

}  // namespace

// Claude Generated (Oct 2026, WP12 stage 2): formulas in qm_xc_fxc.h.
bool buildFxc(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P, const Matrix& D,
              const XCFunctional& f, double rho_threshold, int threads, Matrix& V, const qmgrid::BlockScreen* screen)
{
    if (!fxc::functionalCovered(f)) return false;
    const int n = ao.nbf();
    BlockLoop L;
    L.n = n;
    L.npts = grid.size();
    L.screened = screen && screen->active();
    L.nblocks = L.screened ? screen->nblocks() : (L.npts + L.block - 1) / L.block;
    const bool gga = f.needsGradient(), mgga = f.needsTau();
    const double rho_min = rhoMin(f, rho_threshold);
    const int nchunks = std::max(1, std::min(32, L.nblocks));
    const Eigen::MatrixXd Pd(P), Dd(D);
    std::vector<Eigen::MatrixXd> Vc(nchunks);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)L.nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)L.nblocks * (c + 1) / nchunks);
        Eigen::MatrixXd Vl = Eigen::MatrixXd::Zero(n, n);
        qmgrid::AOBlock blk;
        std::vector<int> cols;
        BlockDensity d0, d1;
        for (int b = b0; b < b1; ++b) {
            int begin = 0, count = 0;
            bool empty = false;
            evaluateBlock(grid, ao, screen, L, b, gga ? 1 : 0, blk, cols, begin, count, empty);
            if (empty) continue;
            const Eigen::MatrixXd phi(blk.phi);
            Eigen::MatrixXd dd[3];
            if (gga) {
                dd[0] = blk.dx;
                dd[1] = blk.dy;
                dd[2] = blk.dz;
            }
            blockDensity(phi, dd, gga, mgga, gatherSub(Pd, cols), d0);
            blockDensity(phi, dd, gga, mgga, gatherSub(Dd, cols), d1);
            Eigen::VectorXd wv = Eigen::VectorXd::Zero(count), wt = Eigen::VectorXd::Zero(count);
            Eigen::MatrixXd gv = Eigen::MatrixXd::Zero(count, 3);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double sigma = gga ? d0.grad.row(q).squaredNorm() : 0.0;
                double v[3], H[6];
                if (!fxc::evaluatePointFxc(f.terms.data(), static_cast<int>(f.terms.size()), rho_min, d0.rho(q), sigma,
                                           d0.tau(q), v, H))
                    continue;
                any = true;
                const double r1 = d1.rho(q), s1 = gga ? 2.0 * d0.grad.row(q).dot(d1.grad.row(q)) : 0.0, t1 = d1.tau(q);
                const double dvr = H[0] * r1 + H[1] * s1 + H[2] * t1;
                const double dvs = H[1] * r1 + H[3] * s1 + H[4] * t1;
                const double dvt = H[2] * r1 + H[4] * s1 + H[5] * t1;
                const double w = grid.weights(begin + q);
                wv(q) = w * dvr;
                wt(q) = 0.5 * w * dvt;
                for (int k = 0; k < 3; ++k) gv(q, k) = 2.0 * w * (dvs * d0.grad(q, k) + v[1] * d1.grad(q, k));
            }
            if (!any) continue;
            Eigen::MatrixXd Y = (0.5 * wv).asDiagonal() * phi;
            if (gga)
                for (int k = 0; k < 3; ++k) Y.noalias() += gv.col(k).asDiagonal() * dd[k];
            Eigen::MatrixXd X = phi.transpose() * Y;
            X += X.transpose().eval();
            if (mgga)
                for (int k = 0; k < 3; ++k) X.noalias() += dd[k].transpose() * (wt.asDiagonal() * dd[k]);
            const int ns = static_cast<int>(cols.size());
            for (int j = 0; j < ns; ++j)
                for (int i = 0; i < ns; ++i) Vl(cols[i], cols[j]) += X(i, j);
        }
        Vc[c] = std::move(Vl);
    }
    Eigen::MatrixXd Vs = Eigen::MatrixXd::Zero(n, n);
    for (int c = 0; c < nchunks; ++c) Vs += Vc[c];
    V = Matrix(0.5 * (Vs + Vs.transpose()));
    return true;
}

// Claude Generated (Oct 2026, WP12 stage 2): formulas in qm_xc_fxc.h.
bool buildXCResponseGradient(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                             const Matrix& P2, const XCFunctional& f, double rho_threshold, int threads,
                             const std::vector<int>& bf_atom, int natoms, XCGradient& out,
                             const qmgrid::BlockScreen* screen)
{
    if (!fxc::functionalCovered(f)) return false;
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    BlockLoop L;
    L.n = n;
    L.npts = grid.size();
    L.screened = screen && screen->active();
    L.nblocks = L.screened ? screen->nblocks() : (L.npts + L.block - 1) / L.block;
    const bool gga = f.needsGradient(), mgga = f.needsTau();
    const double rho_min = rhoMin(f, rho_threshold);
    const int nchunks = std::max(1, std::min(32, L.nblocks));
    const Eigen::MatrixXd Pd(P), P2d(P2);
    std::vector<Matrix> gb(nchunks);
    std::vector<long> used_c(nchunks, 0);
    out.edens = Vector::Zero(L.npts);
    Matrix spoint = Matrix::Zero(L.npts, 3);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)L.nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)L.nblocks * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        long used = 0;
        qmgrid::AOBlock blk;
        std::vector<int> cols;
        BlockDensity d0, d2;
        for (int b = b0; b < b1; ++b) {
            int begin = 0, count = 0;
            bool empty = false;
            evaluateBlock(grid, ao, screen, L, b, gga ? 2 : 1, blk, cols, begin, count, empty);
            if (empty) continue;
            const Eigen::MatrixXd phi(blk.phi);
            const Eigen::MatrixXd dd[3] = { Eigen::MatrixXd(blk.dx), Eigen::MatrixXd(blk.dy), Eigen::MatrixXd(blk.dz) };
            const Eigen::MatrixXd Ps = gatherSub(Pd, cols), P2s = gatherSub(P2d, cols);
            blockDensity(phi, dd, gga, mgga, Ps, d0);
            blockDensity(phi, dd, gga, mgga, P2s, d2);
            Eigen::VectorXd wv = Eigen::VectorXd::Zero(count), dwv = Eigen::VectorXd::Zero(count);
            Eigen::VectorXd wt = Eigen::VectorXd::Zero(count), dwt = Eigen::VectorXd::Zero(count);
            Eigen::MatrixXd vj = Eigen::MatrixXd::Zero(count, 3), dvj = Eigen::MatrixXd::Zero(count, 3);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double sigma = gga ? d0.grad.row(q).squaredNorm() : 0.0;
                double v[3], H[6];
                if (!fxc::evaluatePointFxc(f.terms.data(), static_cast<int>(f.terms.size()), rho_min, d0.rho(q), sigma,
                                           d0.tau(q), v, H))
                    continue;
                any = true;
                ++used;
                const double r2 = d2.rho(q), s2 = gga ? 2.0 * d0.grad.row(q).dot(d2.grad.row(q)) : 0.0, t2 = d2.tau(q);
                const double dvr = H[0] * r2 + H[1] * s2 + H[2] * t2;
                const double dvs = H[1] * r2 + H[3] * s2 + H[4] * t2;
                const double dvt = H[2] * r2 + H[4] * s2 + H[5] * t2;
                out.edens(begin + q) = v[0] * r2 + v[1] * s2 + v[2] * t2;
                const double w = grid.weights(begin + q);
                wv(q) = w * v[0];
                dwv(q) = w * dvr;
                wt(q) = 0.5 * w * v[2];
                dwt(q) = 0.5 * w * dvt;
                for (int k = 0; k < 3; ++k) {
                    vj(q, k) = 2.0 * w * v[1] * d0.grad(q, k);
                    dvj(q, k) = 2.0 * w * (dvs * d0.grad(q, k) + v[1] * d2.grad(q, k));
                }
            }
            if (!any) continue;
            Eigen::MatrixXd Zm = wv.asDiagonal() * phi, dZ = dwv.asDiagonal() * phi;
            if (gga)
                for (int k = 0; k < 3; ++k) {
                    Zm.noalias() += vj.col(k).asDiagonal() * dd[k];
                    dZ.noalias() += dvj.col(k).asDiagonal() * dd[k];
                }
            const Eigen::MatrixXd ZP = Zm * P2s + dZ * Ps;
            Eigen::MatrixXd T[3] = { dd[0].cwiseProduct(ZP), dd[1].cwiseProduct(ZP), dd[2].cwiseProduct(ZP) };
            if (gga) {
                const Eigen::MatrixXd xx(blk.dxx), xy(blk.dxy), xz(blk.dxz), yy(blk.dyy), yz(blk.dyz), zz(blk.dzz);
                // second AO derivatives d_k d_j phi: row k = {xx, xy, xz}, {xy, yy, yz}, {xz, yz, zz}
                const Eigen::MatrixXd* D2[3][3] = { { &xx, &xy, &xz }, { &xy, &yy, &yz }, { &xz, &yz, &zz } };
                for (int k = 0; k < 3; ++k)
                    for (int j = 0; j < 3; ++j) {
                        const Eigen::MatrixXd& dkj = *D2[k][j];
                        T[k].noalias() += (vj.col(j).asDiagonal() * dkj).cwiseProduct(d2.Xphi);
                        T[k].noalias() += (dvj.col(j).asDiagonal() * dkj).cwiseProduct(d0.Xphi);
                        if (mgga) {
                            T[k].noalias() += (wt.asDiagonal() * dkj).cwiseProduct(d2.dX[j]);
                            T[k].noalias() += (dwt.asDiagonal() * dkj).cwiseProduct(d0.dX[j]);
                        }
                    }
            }
            const int ns = static_cast<int>(cols.size());
            for (int k = 0; k < 3; ++k) {
                const Eigen::VectorXd colsum = T[k].colwise().sum().transpose();
                for (int m = 0; m < ns; ++m) g(bf_atom[cols[m]], k) += -2.0 * colsum(m);
                const Eigen::VectorXd rowsum = T[k].rowwise().sum();
                for (int q = 0; q < count; ++q) spoint(begin + q, k) = -2.0 * rowsum(q);
            }
        }
        gb[c] = std::move(g);
        used_c[c] = used;
    }
    out.basis = Matrix::Zero(natoms, 3);
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {
        out.basis += gb[c];
        out.points_used += used_c[c];
    }
    out.points = Matrix::Zero(natoms, 3);
    for (int p = 0; p < L.npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

}  // namespace qmxc
