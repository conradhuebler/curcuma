/*
 * <Native QM: spin-polarised exchange-correlation for open-shell Kohn-Sham (UKS), stage 5b>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas and literature: qm_xc_spin.h. The point kernels (dual numbers, spin-polarised correlation, exchange by
 * spin scaling) are in qm_xc_spin_kernels.h since stage 6a (Nov 2026), shared with the device; this file keeps the
 * grid loops (V_xc, XC gradient) of the CPU path.
 *
 * Claude Generated (Nov 2026, stage 5b). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_xc_spin.h"
#include "qm_xc_spin_kernels.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace qmxc {

XCResultSpin::XCResultSpin() = default;
XCResultSpin::XCResultSpin(const XCResultSpin&) = default;
XCResultSpin& XCResultSpin::operator=(const XCResultSpin&) = default;
XCResultSpin::~XCResultSpin() = default;

bool supportsSpinPolarised(const XCFunctional& f, std::string& why)
{
    for (const XCTerm& t : f.terms) {
        if (t.weight == 0.0) continue;
        switch (t.kernel) {
            case XCKernel::SlaterX:
            case XCKernel::PBE_X:
            case XCKernel::B88_X:
            case XCKernel::VWN5_C:
            case XCKernel::VWNRPA_C:
            case XCKernel::PW92_C:
            case XCKernel::PBE_C:
            case XCKernel::LYP_C:
                break;
            default:
                why = "no spin-polarised form of the kernel " + kernelName(t.kernel)
                    + " (open-shell Kohn-Sham: Slater, B88, PBE-X, VWN5, VWN-RPA, PW92, PBE-C, LYP, i.e. lda, pbe, b3lyp)";
                return false;
        }
    }
    return true;
}

bool evaluatePointSpin(const XCFunctional& f, bool gga, double rho_min, double ra, double rb, double saa, double sab,
                       double sbb, double& ex, double& ec, double v[5])
{
    // stage 6a (Nov 2026): the body lives in qm_xc_spin_kernels.h, shared with the device (cuda/qm_gpu_xc.cu)
    return kern::evaluatePointSpinTerms(f.terms.data(), static_cast<int>(f.terms.size()), gga, rho_min, ra, rb, saa, sab,
                                        sbb, ex, ec, v);
}

void buildVxcSpin(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& Pa, const Matrix& Pb,
                  const XCFunctional& f, double rho_threshold, int threads, XCResultSpin& out, int block)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = kern::effectiveThreshold(gga, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    // fixed chunk count (basis size only, two matrices per chunk): thread-count independent result
    const double mat_mb = 16.0 * double(n) * n / 1.0e6;
    const int nchunks = std::max(1, std::min({ 32, nblocks, static_cast<int>(128.0 / std::max(mat_mb, 1e-3)) }));
    std::vector<Eigen::MatrixXd> Vac(nchunks), Vbc(nchunks);
    std::vector<double> ex_c(nchunks, 0.0), ec_c(nchunks, 0.0), ne_c(nchunks, 0.0);
    std::vector<long> used_c(nchunks, 0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Eigen::MatrixXd Va = Eigen::MatrixXd::Zero(n, n), Vb = Eigen::MatrixXd::Zero(n, n);
        double ex_sum = 0.0, ec_sum = 0.0, ne_sum = 0.0;
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix ga, gb;
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 1 : 0, blk);
            Vector ra, rb;
            qmgrid::densityOnBlock(blk, Pa, ra, gga ? &ga : nullptr);
            qmgrid::densityOnBlock(blk, Pb, rb, gga ? &gb : nullptr);
            Eigen::VectorXd wva = Eigen::VectorXd::Zero(count), wvb = Eigen::VectorXd::Zero(count);
            Eigen::MatrixXd ua = Eigen::MatrixXd::Zero(count, 3), ub = Eigen::MatrixXd::Zero(count, 3);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                ne_sum += w * (ra(q) + rb(q));
                double saa = 0.0, sab = 0.0, sbb = 0.0;
                if (gga) {
                    saa = ga.row(q).squaredNorm();
                    sbb = gb.row(q).squaredNorm();
                    sab = ga.row(q).dot(gb.row(q));
                }
                double ex = 0.0, ec = 0.0, v[5];
                if (!evaluatePointSpin(f, gga, rho_min, ra(q), rb(q), saa, sab, sbb, ex, ec, v)) continue;
                ex_sum += w * ex;
                ec_sum += w * ec;
                ++used;
                any = true;
                wva(q) = w * v[0];
                wvb(q) = w * v[1];
                if (gga)
                    for (int k = 0; k < 3; ++k) {
                        ua(q, k) = w * (2.0 * v[2] * ga(q, k) + v[3] * gb(q, k));
                        ub(q, k) = w * (2.0 * v[4] * gb(q, k) + v[3] * ga(q, k));
                    }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi);
            if (!gga) {
                Va.noalias() += phi.transpose() * (wva.asDiagonal() * phi);
                Vb.noalias() += phi.transpose() * (wvb.asDiagonal() * phi);
                continue;
            }
            const Eigen::MatrixXd dx(blk.dx), dy(blk.dy), dz(blk.dz);
            // Y_s = (w v_rho_s / 2) phi + u_s . grad phi;  V_s += phi^T Y_s + Y_s^T phi
            Eigen::MatrixXd Ya = (0.5 * wva).asDiagonal() * phi;
            Ya.noalias() += ua.col(0).asDiagonal() * dx;
            Ya.noalias() += ua.col(1).asDiagonal() * dy;
            Ya.noalias() += ua.col(2).asDiagonal() * dz;
            Eigen::MatrixXd Yb = (0.5 * wvb).asDiagonal() * phi;
            Yb.noalias() += ub.col(0).asDiagonal() * dx;
            Yb.noalias() += ub.col(1).asDiagonal() * dy;
            Yb.noalias() += ub.col(2).asDiagonal() * dz;
            const Eigen::MatrixXd Xa = phi.transpose() * Ya, Xb = phi.transpose() * Yb;
            Va.noalias() += Xa;
            Va.noalias() += Xa.transpose();
            Vb.noalias() += Xb;
            Vb.noalias() += Xb.transpose();
        }
        Vac[c] = std::move(Va);
        Vbc[c] = std::move(Vb);
        ex_c[c] = ex_sum;
        ec_c[c] = ec_sum;
        ne_c[c] = ne_sum;
        used_c[c] = used;
    }
    Eigen::MatrixXd Va = Eigen::MatrixXd::Zero(n, n), Vb = Eigen::MatrixXd::Zero(n, n);
    out.ex = out.ec = out.nelec = 0.0;
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {  // chunk order: fixed summation order
        Va += Vac[c];
        Vb += Vbc[c];
        out.ex += ex_c[c];
        out.ec += ec_c[c];
        out.nelec += ne_c[c];
        out.points_used += used_c[c];
    }
    out.exc = out.ex + out.ec;
    out.Va = Matrix(0.5 * (Va + Va.transpose()));
    out.Vb = Matrix(0.5 * (Vb + Vb.transpose()));
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

void buildXCGradientSpin(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& Pa,
                         const Matrix& Pb, const XCFunctional& f, double rho_threshold, int threads,
                         const std::vector<int>& bf_atom, int natoms, XCGradient& out, int block)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = kern::effectiveThreshold(gga, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    const int nchunks = std::max(1, std::min(32, nblocks));
    const Eigen::MatrixXd Pad(Pa), Pbd(Pb);
    std::vector<Matrix> gbc(nchunks);
    std::vector<long> used_c(nchunks, 0);
    out.edens = Vector::Zero(npts);
    Matrix spoint = Matrix::Zero(npts, 3);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix ga, gb;
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 2 : 1, blk);
            Vector ra, rb;
            qmgrid::densityOnBlock(blk, Pa, ra, gga ? &ga : nullptr);
            qmgrid::densityOnBlock(blk, Pb, rb, gga ? &gb : nullptr);
            Eigen::VectorXd wv[2] = { Eigen::VectorXd::Zero(count), Eigen::VectorXd::Zero(count) };
            Eigen::MatrixXd u[2] = { Eigen::MatrixXd::Zero(count, 3), Eigen::MatrixXd::Zero(count, 3) };
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                double saa = 0.0, sab = 0.0, sbb = 0.0;
                if (gga) {
                    saa = ga.row(q).squaredNorm();
                    sbb = gb.row(q).squaredNorm();
                    sab = ga.row(q).dot(gb.row(q));
                }
                double ex = 0.0, ec = 0.0, v[5];
                if (!evaluatePointSpin(f, gga, rho_min, ra(q), rb(q), saa, sab, sbb, ex, ec, v)) continue;
                out.edens(begin + q) = ex + ec;
                ++used;
                any = true;
                wv[0](q) = w * v[0];
                wv[1](q) = w * v[1];
                if (gga)
                    for (int k = 0; k < 3; ++k) {
                        u[0](q, k) = w * (2.0 * v[2] * ga(q, k) + v[3] * gb(q, k));
                        u[1](q, k) = w * (2.0 * v[4] * gb(q, k) + v[3] * ga(q, k));
                    }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi), dx(blk.dx), dy(blk.dy), dz(blk.dz);
            Eigen::MatrixXd xx, xy, xz, yy, yz, zz;
            if (gga) {
                xx = blk.dxx;
                xy = blk.dxy;
                xz = blk.dxz;
                yy = blk.dyy;
                yz = blk.dyz;
                zz = blk.dzz;
            }
            Eigen::MatrixXd rowsum_all = Eigen::MatrixXd::Zero(count, 3);
            for (int s = 0; s < 2; ++s) {
                const Eigen::MatrixXd& Ps = s == 0 ? Pad : Pbd;
                Eigen::MatrixXd Zm = wv[s].asDiagonal() * phi;
                if (gga) {
                    Zm.noalias() += u[s].col(0).asDiagonal() * dx;
                    Zm.noalias() += u[s].col(1).asDiagonal() * dy;
                    Zm.noalias() += u[s].col(2).asDiagonal() * dz;
                }
                const Eigen::MatrixXd ZP = Zm * Ps;
                Eigen::MatrixXd T[3] = { dx.cwiseProduct(ZP), dy.cwiseProduct(ZP), dz.cwiseProduct(ZP) };
                if (gga) {
                    const Eigen::MatrixXd PhiP = phi * Ps;
                    const auto ux = u[s].col(0).asDiagonal(), uy = u[s].col(1).asDiagonal(), uz = u[s].col(2).asDiagonal();
                    T[0].noalias() += ((ux * xx) + (uy * xy) + (uz * xz)).cwiseProduct(PhiP);
                    T[1].noalias() += ((ux * xy) + (uy * yy) + (uz * yz)).cwiseProduct(PhiP);
                    T[2].noalias() += ((ux * xz) + (uy * yz) + (uz * zz)).cwiseProduct(PhiP);
                }
                for (int k = 0; k < 3; ++k) {
                    const Eigen::VectorXd colsum = T[k].colwise().sum().transpose();
                    for (int m = 0; m < n; ++m) g(bf_atom[m], k) += -2.0 * colsum(m);
                    rowsum_all.col(k) += T[k].rowwise().sum();
                }
            }
            for (int q = 0; q < count; ++q)
                for (int k = 0; k < 3; ++k) spoint(begin + q, k) = -2.0 * rowsum_all(q, k);
        }
        gbc[c] = std::move(g);
        used_c[c] = used;
    }
    out.basis = Matrix::Zero(natoms, 3);
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {
        out.basis += gbc[c];
        out.points_used += used_c[c];
    }
    out.points = Matrix::Zero(natoms, 3);
    for (int p = 0; p < npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace qmxc
