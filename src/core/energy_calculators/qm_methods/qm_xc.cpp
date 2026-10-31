/*
 * <Native QM: exchange-correlation functionals of closed-shell Kohn-Sham DFT (WP5: LDA, WP6: PBE, WP7: B3LYP,
 *  WP11: meta-GGA TPSS, SCAN)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas and literature: qm_xc.h. Checked point by point against libxc through PySCF
 * (pyscf.dft.libxc.eval_xc 'LDA_X', 'LDA_C_VWN', 'LDA_C_VWN_RPA'; ctest qm_lda_kernel; 'GGA_X_PBE', 'GGA_C_PBE',
 * 'LDA_C_PW_MOD'; ctest qm_pbe_kernel; 'GGA_X_B88', 'GGA_C_LYP'; ctest qm_b3lyp_kernel; 'MGGA_X_TPSS', 'MGGA_C_TPSS',
 * 'MGGA_X_SCAN', 'MGGA_C_SCAN'; ctest qm_mgga_kernel).
 *
 * xcDFT (TCCM winter school 2019) has a Slater exchange (lda_exchange_*.f90) and a "correlation" whose energy
 * formula Ec = -a sum rho w / (1 + b rho^(-1/3) w) puts the grid weight into the denominator; that is not a
 * functional of the density. It is not ported; VWN replaces it.
 *
 * Claude Generated (Oct 2026, WP5; WP6 PBE; WP7 B88, LYP, B3LYP; WP11 meta-GGA). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_xc.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

namespace qmxc {

namespace {
/// Density below which a point is left out of E_xc and V_xc for functional f (shared rule, qm_xc_kernels.h).
double effectiveThreshold(const XCFunctional& f, double rho_threshold)
{
    return kern::effectiveThreshold(f.needsGradient(), rho_threshold);
}
}  // namespace

bool parseVWNVariant(const std::string& name, VWNVariant& v)
{
    std::string s = name;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "vwn5" || s == "vwn") {
        v = VWNVariant::VWN5;
        return true;
    }
    if (s == "vwn_rpa" || s == "vwn-rpa" || s == "rpa" || s == "vwn3") {
        v = VWNVariant::VWN_RPA;
        return true;
    }
    return false;
}

std::string vwnName(VWNVariant v) { return v == VWNVariant::VWN5 ? "VWN5" : "VWN-RPA"; }

bool isExchangeKernel(XCKernel k) { return kern::isExchange(k); }

bool isGGAKernel(XCKernel k) { return kern::isGGA(k); }

bool isMGGAKernel(XCKernel k) { return kern::isMGGA(k); }

std::string kernelName(XCKernel k)
{
    switch (k) {
        case XCKernel::SlaterX:  return "Slater";
        case XCKernel::VWN5_C:   return "VWN5";
        case XCKernel::VWNRPA_C: return "VWN-RPA";
        case XCKernel::PW92_C:   return "PW92";
        case XCKernel::PBE_X:    return "PBE-X";
        case XCKernel::PBE_C:    return "PBE-C";
        case XCKernel::B88_X:    return "B88";
        case XCKernel::LYP_C:    return "LYP";
        case XCKernel::TPSS_X:   return "TPSS-X";
        case XCKernel::TPSS_C:   return "TPSS-C";
        case XCKernel::SCAN_X:   return "SCAN-X";
        case XCKernel::SCAN_C:   return "SCAN-C";
        case XCKernel::P86_C:    return "P86(PZ81)";
        case XCKernel::P86VWN_C: return "P86(VWN5)";
        case XCKernel::B95_C:    return "B95";
    }
    return "?";
}

XCFunctional::XCFunctional() = default;
XCFunctional::XCFunctional(const XCFunctional&) = default;
XCFunctional& XCFunctional::operator=(const XCFunctional&) = default;
XCFunctional::~XCFunctional() = default;

bool XCFunctional::needsGradient() const
{
    for (const XCTerm& t : terms)
        if ((isGGAKernel(t.kernel) || isMGGAKernel(t.kernel)) && t.weight != 0.0) return true;
    return false;
}

bool XCFunctional::needsTau() const
{
    for (const XCTerm& t : terms)
        if (isMGGAKernel(t.kernel) && t.weight != 0.0) return true;
    return false;
}

std::string XCFunctional::label() const
{
    if (!hasXC()) return "HF (no density functional)";
    if (!description.empty()) return description;
    std::string s;
    for (const XCTerm& t : terms) s += (s.empty() ? "" : " + ") + fmt::format("{:g} {}", t.weight, kernelName(t.kernel));
    return s;
}

XCFunctional ldaFunctional(VWNVariant v)
{
    XCFunctional f;
    f.name = "lda";
    f.exact_exchange = 0.0;
    f.terms = { { XCKernel::SlaterX, 1.0 }, { v == VWNVariant::VWN5 ? XCKernel::VWN5_C : XCKernel::VWNRPA_C, 1.0 } };
    f.description = std::string("SVWN") + (v == VWNVariant::VWN5 ? "5" : "-RPA") + " (Slater exchange + " + vwnName(v)
        + " correlation)";
    return f;
}

XCFunctional pbeFunctional()
{
    XCFunctional f;
    f.name = "pbe";
    f.exact_exchange = 0.0;
    f.terms = { { XCKernel::PBE_X, 1.0 }, { XCKernel::PBE_C, 1.0 } };
    f.description = "PBE-GGA (PBE exchange + PBE correlation, PW92 LDA part)";
    return f;
}

// WP7: B3LYP (Stephens et al. 1994) = 0.08 Slater + 0.72 B88 + 0.20 HF + 0.19 VWN + 0.81 LYP. The VWN variant is a
// choice of the caller (QMEngine: -qm.vwn_variant, default vwn5 as ORCA's B3LYP); it is always named in the label.
XCFunctional b3lypFunctional(VWNVariant v)
{
    XCFunctional f;
    f.name = "b3lyp";
    f.exact_exchange = 0.20;
    f.terms = { { XCKernel::SlaterX, 0.08 }, { XCKernel::B88_X, 0.72 },
                { v == VWNVariant::VWN5 ? XCKernel::VWN5_C : XCKernel::VWNRPA_C, 0.19 }, { XCKernel::LYP_C, 0.81 } };
    f.description = std::string("B3LYP-") + vwnName(v) + " (0.08 Slater + 0.72 B88 + 0.20 HF exchange, 0.19 "
        + vwnName(v) + " + 0.81 LYP correlation)";
    return f;
}

// WP11 (Claude Generated, Oct 2026): meta-GGA functionals as table rows (name, label, a_x, terms). A global
// meta-GGA hybrid would be one more row (e.g. TPSSh: { "tpssh", ..., 0.10, { { TPSS_X, 0.90 }, { TPSS_C, 1.0 } } },
// Staroverov et al. 2003), r2SCAN needs its own kernels first; neither is implemented or validated.
namespace {
struct MetaGGARow {
    const char* name;
    const char* description;
    double exact_exchange;
    XCTerm terms[2];
};
const MetaGGARow kMetaGGATable[] = {
    { "tpss", "TPSS meta-GGA (TPSS exchange + TPSS correlation, Tao et al. 2003)", 0.0,
      { { XCKernel::TPSS_X, 1.0 }, { XCKernel::TPSS_C, 1.0 } } },
    { "scan", "SCAN meta-GGA (SCAN exchange + SCAN correlation, Sun et al. 2015)", 0.0,
      { { XCKernel::SCAN_X, 1.0 }, { XCKernel::SCAN_C, 1.0 } } },
};
}  // namespace

bool metaGGAFunctional(const std::string& name, XCFunctional& f)
{
    for (const MetaGGARow& row : kMetaGGATable) {
        if (name != row.name) continue;
        f = XCFunctional();
        f.name = row.name;
        f.description = row.description;
        f.exact_exchange = row.exact_exchange;
        f.terms.assign(std::begin(row.terms), std::end(row.terms));
        return true;
    }
    return false;
}

// The point kernels live in qm_xc_kernels.h (one source for CPU and GPU, Oct 2026); formulas and derivations there.
void slaterExchange(double rho, double& eps, double& v) { kern::slaterExchange(rho, eps, v); }

void vwnCorrelation(double rho, VWNVariant variant, double& eps, double& v)
{
    kern::vwnCorrelation(rho, variant == VWNVariant::VWN_RPA, eps, v);
}

void pw92Correlation(double rho, double& eps, double& v) { kern::pw92Correlation(rho, eps, v); }

void pbeExchange(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    kern::pbeExchange(rho, sigma, eps, vrho, vsigma);
}

void pbeCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    kern::pbeCorrelation(rho, sigma, eps, vrho, vsigma);
}

void b88Exchange(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    kern::b88Exchange(rho, sigma, eps, vrho, vsigma);
}

void lypCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    kern::lypCorrelation(rho, sigma, eps, vrho, vsigma);
}

void evaluateKernel(XCKernel k, double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    if (!kern::evaluateKernel(k, rho, sigma, eps, vrho, vsigma)) eps = vrho = std::numeric_limits<double>::quiet_NaN();
}

void evaluateKernelMGGA(XCKernel k, double rho, double sigma, double tau, double& eps, double& vrho, double& vsigma,
                        double& vtau)
{
    if (!kern::evaluateKernelMGGA(k, rho, sigma, tau, eps, vrho, vsigma, vtau))
        eps = vrho = vsigma = vtau = std::numeric_limits<double>::quiet_NaN();
}

void evaluateXCMeta(const XCFunctional& f, const Vector& rho, const Vector& sigma, const Vector& tau,
                    double rho_threshold, Vector& ex, Vector& ec, Vector& vrho, Vector& vsigma, Vector& vtau)
{
    const int n = static_cast<int>(rho.size());
    const double rho_min = effectiveThreshold(f, rho_threshold);
    ex = Vector::Zero(n);
    ec = Vector::Zero(n);
    vrho = Vector::Zero(n);
    vsigma = Vector::Zero(n);
    vtau = Vector::Zero(n);
    for (int p = 0; p < n; ++p)
        kern::evaluatePointMGGA(f.terms.data(), static_cast<int>(f.terms.size()), rho_min, rho(p), sigma(p), tau(p),
                                ex(p), ec(p), vrho(p), vsigma(p), vtau(p));
}

void evaluateXC(const XCFunctional& f, const Vector& rho, const Vector& sigma, double rho_threshold, Vector& ex,
                Vector& ec, Vector& vrho, Vector& vsigma)
{
    const int n = static_cast<int>(rho.size());
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    ex = Vector::Zero(n);
    ec = Vector::Zero(n);
    vrho = Vector::Zero(n);
    vsigma = Vector::Zero(n);
    for (int p = 0; p < n; ++p)
        kern::evaluatePoint(f.terms.data(), static_cast<int>(f.terms.size()), gga, rho_min, rho(p), gga ? sigma(p) : 0.0,
                            ex(p), ec(p), vrho(p), vsigma(p));
}

namespace {

/// P restricted to the kept functions: P_SS(i, j) = P(cols[i], cols[j]).
Eigen::MatrixXd gatherP(const Eigen::MatrixXd& P, const std::vector<int>& cols)
{
    const int ns = static_cast<int>(cols.size());
    Eigen::MatrixXd S(ns, ns);
    for (int j = 0; j < ns; ++j)
        for (int i = 0; i < ns; ++i) S(i, j) = P(cols[i], cols[j]);
    return S;
}

/// rho and (GGA) grad rho on a block from the kept functions: rho_q = sum_ij phi_qi P_ij phi_qj.
void densitySubset(const Eigen::MatrixXd& phi, const Eigen::MatrixXd* dphi, const Eigen::MatrixXd& Ps, Vector& rho,
                   Matrix* grad, Eigen::MatrixXd& Pphi)
{
    const int nb = static_cast<int>(phi.rows());
    Pphi.noalias() = phi * Ps;
    rho.resize(nb);
    for (int q = 0; q < nb; ++q) rho(q) = phi.row(q).dot(Pphi.row(q));
    if (grad) {
        grad->resize(nb, 3);
        for (int q = 0; q < nb; ++q)
            for (int k = 0; k < 3; ++k) (*grad)(q, k) = 2.0 * dphi[k].row(q).dot(Pphi.row(q));
    }
}

/// buildVxc on the blocks of an active AO screen (qmgrid::BlockScreen): per block only the kept functions S,
/// X_SS = phi_S^T Y_S scattered into the chunk matrix over U_c (union of the chunk's functions), the chunk matrices
/// added into V in chunk order. Same formulas as the unscreened path.
/// Claude Generated (Oct 2026, AO screening).
void buildVxcScreened(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                      const XCFunctional& f, double rho_threshold, int threads, XCResult& out,
                      const qmgrid::BlockScreen& sc)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    const int nblocks = sc.nblocks();
    // Chunks of consecutive (spatially compact) blocks; every chunk accumulates V only over the union U_c of the
    // functions of its blocks (|U_c| x |U_c|, not n x n), so the number of chunks need not shrink with n. Start with
    // 64 chunks and halve while the chunk matrices need more than 1 GB: the count depends only on grid and basis,
    // never on `threads` (bitwise thread-independent).
    std::vector<char> fmark(n, 0);
    auto unionOf = [&](int b0, int b1, std::vector<int>& U) {
        U.clear();
        for (int b = b0; b < b1; ++b)
            for (int k = sc.shell_off[b]; k < sc.shell_off[b + 1]; ++k) {
                const int sh = sc.shells[k];
                for (int w = 0; w < ao.shellFunctions(sh); ++w) fmark[ao.shellFirstFunction(sh) + w] = 1;
            }
        for (int m = 0; m < n; ++m)
            if (fmark[m]) {
                U.push_back(m);
                fmark[m] = 0;
            }
    };
    int nchunks = std::max(1, std::min(64, nblocks));
    std::vector<std::vector<int>> Uc;
    while (true) {
        Uc.assign(nchunks, {});
        double mb = 0.0;
        for (int c = 0; c < nchunks; ++c) {
            unionOf(static_cast<int>((long long)nblocks * c / nchunks), static_cast<int>((long long)nblocks * (c + 1) / nchunks), Uc[c]);
            mb += 8.0 * double(Uc[c].size()) * double(Uc[c].size()) / 1.0e6;
        }
        if (mb <= 1024.0 || nchunks == 1) break;
        nchunks = std::max(1, nchunks / 2);
    }
    const Eigen::MatrixXd Pd(P);
    std::vector<Eigen::MatrixXd> Vc(nchunks);
    std::vector<double> exc_c(nchunks, 0.0), ec_c(nchunks, 0.0), ne_c(nchunks, 0.0);
    std::vector<long> used_c(nchunks, 0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        const std::vector<int>& U = Uc[c];
        std::vector<int> pos(n, -1);  // active function -> row/column of the chunk matrix
        for (size_t i = 0; i < U.size(); ++i) pos[U[i]] = static_cast<int>(i);
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(U.size(), U.size());
        double ex_sum = 0.0, ec_sum = 0.0, ne_sum = 0.0;
        long used = 0;
        qmgrid::AOBlock blk;
        std::vector<int> cols, lc;
        Matrix grad;
        Eigen::MatrixXd Pphi;
        for (int b = b0; b < b1; ++b) {
            const int begin = sc.begin[b], count = sc.count[b];
            const int off = sc.shell_off[b], nsh = sc.shell_off[b + 1] - off;
            if (nsh == 0) continue;  // no function reaches the block: rho = 0 there
            ao.evaluateShells(grid.points, begin, count, gga ? 1 : 0, sc.shells.data() + off, nsh, blk, cols);
            const Eigen::MatrixXd phi(blk.phi);  // count x ns, column-major
            Eigen::MatrixXd d[3];
            if (gga) {
                d[0] = blk.dx;
                d[1] = blk.dy;
                d[2] = blk.dz;
            }
            const Eigen::MatrixXd Ps = gatherP(Pd, cols);
            Vector rho, sigma;
            densitySubset(phi, d, Ps, rho, gga ? &grad : nullptr, Pphi);
            if (gga) sigma = grad.rowwise().squaredNorm();
            Vector ex, ec, vrho, vsigma;
            evaluateXC(f, rho, sigma, rho_threshold, ex, ec, vrho, vsigma);
            Eigen::VectorXd wv(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                ne_sum += w * rho(q);
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    ex_sum += w * rho(q) * ex(q);
                    ec_sum += w * rho(q) * ec(q);
                    ++used;
                    any = true;
                }
            }
            if (!any) continue;
            const int ns = static_cast<int>(cols.size());
            lc.resize(ns);
            for (int i = 0; i < ns; ++i) lc[i] = pos[cols[i]];
            if (!gga) {
                const Eigen::MatrixXd wphi = wv.asDiagonal() * phi;
                const Eigen::MatrixXd X = phi.transpose() * wphi;
                for (int j = 0; j < ns; ++j)
                    for (int i = 0; i < ns; ++i) V(lc[i], lc[j]) += X(i, j);
                continue;
            }
            Eigen::VectorXd gx(count), gy(count), gz(count);
            for (int q = 0; q < count; ++q) {
                const double w2s = 2.0 * grid.weights(begin + q) * vsigma(q);
                gx(q) = w2s * grad(q, 0);
                gy(q) = w2s * grad(q, 1);
                gz(q) = w2s * grad(q, 2);
            }
            Eigen::MatrixXd Y = (0.5 * wv).asDiagonal() * phi;
            Y.noalias() += gx.asDiagonal() * d[0];
            Y.noalias() += gy.asDiagonal() * d[1];
            Y.noalias() += gz.asDiagonal() * d[2];
            const Eigen::MatrixXd X = phi.transpose() * Y;
            for (int j = 0; j < ns; ++j)
                for (int i = 0; i < ns; ++i) {
                    double& v = V(lc[i], lc[j]);
                    v += X(i, j);
                    v += X(j, i);
                }
        }
        Vc[c] = std::move(V);
        exc_c[c] = ex_sum;
        ec_c[c] = ec_sum;
        ne_c[c] = ne_sum;
        used_c[c] = used;
    }
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
    out.ex = out.ec = out.nelec = 0.0;
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {  // chunk order: fixed summation order
        const std::vector<int>& U = Uc[c];
        const int nu = static_cast<int>(U.size());
        for (int j = 0; j < nu; ++j)
            for (int i = 0; i < nu; ++i) V(U[i], U[j]) += Vc[c](i, j);
        out.ex += exc_c[c];
        out.ec += ec_c[c];
        out.nelec += ne_c[c];
        out.points_used += used_c[c];
    }
    out.exc = out.ex + out.ec;
    out.V = Matrix(0.5 * (V + V.transpose()));
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

/// buildXCGradient on the blocks of an active AO screen (kept functions only, as buildVxcScreened).
/// Claude Generated (Oct 2026, AO screening).
void buildXCGradientScreened(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                             const XCFunctional& f, double rho_threshold, int threads, const std::vector<int>& bf_atom,
                             int natoms, XCGradient& out, const qmgrid::BlockScreen& sc)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    const int nblocks = sc.nblocks();
    const int nchunks = std::max(1, std::min(32, nblocks));
    const Eigen::MatrixXd Pd(P);
    std::vector<Matrix> gb(nchunks);
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
        std::vector<int> cols;
        Matrix grad;
        Eigen::MatrixXd PhiP;
        for (int b = b0; b < b1; ++b) {
            const int begin = sc.begin[b], count = sc.count[b];
            const int off = sc.shell_off[b], nsh = sc.shell_off[b + 1] - off;
            if (nsh == 0) continue;
            ao.evaluateShells(grid.points, begin, count, gga ? 2 : 1, sc.shells.data() + off, nsh, blk, cols);
            const Eigen::MatrixXd phi(blk.phi), dx(blk.dx), dy(blk.dy), dz(blk.dz);
            const Eigen::MatrixXd d[3] = { dx, dy, dz };
            const Eigen::MatrixXd Ps = gatherP(Pd, cols);
            Vector rho, sigma;
            densitySubset(phi, d, Ps, rho, gga ? &grad : nullptr, PhiP);  // PhiP = phi P_SS
            if (gga) sigma = grad.rowwise().squaredNorm();
            Vector ex, ec, vrho, vsigma;
            evaluateXC(f, rho, sigma, rho_threshold, ex, ec, vrho, vsigma);
            Eigen::VectorXd wv(count), vx(count), vy(count), vz(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                vx(q) = vy(q) = vz(q) = 0.0;
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    out.edens(begin + q) = rho(q) * (ex(q) + ec(q));
                    ++used;
                    any = true;
                    if (gga) {
                        const double w2s = 2.0 * w * vsigma(q);
                        vx(q) = w2s * grad(q, 0);
                        vy(q) = w2s * grad(q, 1);
                        vz(q) = w2s * grad(q, 2);
                    }
                }
            }
            if (!any) continue;
            Eigen::MatrixXd Zm = wv.asDiagonal() * phi;
            if (gga) {
                Zm.noalias() += vx.asDiagonal() * dx;
                Zm.noalias() += vy.asDiagonal() * dy;
                Zm.noalias() += vz.asDiagonal() * dz;
            }
            const Eigen::MatrixXd ZP = Zm * Ps;
            Eigen::MatrixXd T[3] = { dx.cwiseProduct(ZP), dy.cwiseProduct(ZP), dz.cwiseProduct(ZP) };
            if (gga) {
                const Eigen::MatrixXd xx(blk.dxx), xy(blk.dxy), xz(blk.dxz), yy(blk.dyy), yz(blk.dyz), zz(blk.dzz);
                T[0].noalias() += ((vx.asDiagonal() * xx) + (vy.asDiagonal() * xy) + (vz.asDiagonal() * xz)).cwiseProduct(PhiP);
                T[1].noalias() += ((vx.asDiagonal() * xy) + (vy.asDiagonal() * yy) + (vz.asDiagonal() * yz)).cwiseProduct(PhiP);
                T[2].noalias() += ((vx.asDiagonal() * xz) + (vy.asDiagonal() * yz) + (vz.asDiagonal() * zz)).cwiseProduct(PhiP);
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
    for (int p = 0; p < npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

/**
 * buildVxc for a meta-GGA functional (WP11 stage 1). Blocks: those of an active AO screen (kept functions only), else
 * uniform `block`-point blocks with all functions. Per block: rho, grad rho, tau = 1/2 sum_k (d_k phi P d_k phi^T)_qq,
 *   V += phi^T Y + Y^T phi + sum_k (d_k phi)^T diag(w v_tau / 2) (d_k phi),
 *   Y = (w v_rho / 2) phi + 2 w v_sigma (grad rho . grad phi)       (the GGA part as in buildVxc).
 * Chunks of consecutive blocks with their own matrix over the union of their functions (64 chunks, halved above 1 GB
 * of chunk matrices; never dependent on `threads`), added in chunk order: 1 vs N threads bitwise.
 * Claude Generated (Oct 2026, WP11).
 */
void buildVxcMeta(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                  const XCFunctional& f, double rho_threshold, int threads, XCResult& out, int block,
                  const qmgrid::BlockScreen* screen)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    const bool screened = screen && screen->active();
    if (block < 1) block = 128;
    const int nblocks = screened ? screen->nblocks() : (npts + block - 1) / block;
    std::vector<char> fmark(n, 0);
    auto unionOf = [&](int b0, int b1, std::vector<int>& U) {
        U.clear();
        if (!screened) {
            for (int m = 0; m < n; ++m) U.push_back(m);
            return;
        }
        for (int b = b0; b < b1; ++b)
            for (int k = screen->shell_off[b]; k < screen->shell_off[b + 1]; ++k) {
                const int sh = screen->shells[k];
                for (int w = 0; w < ao.shellFunctions(sh); ++w) fmark[ao.shellFirstFunction(sh) + w] = 1;
            }
        for (int m = 0; m < n; ++m)
            if (fmark[m]) {
                U.push_back(m);
                fmark[m] = 0;
            }
    };
    int nchunks = std::max(1, std::min(64, nblocks));
    std::vector<std::vector<int>> Uc;
    while (true) {
        Uc.assign(nchunks, {});
        double mb = 0.0;
        for (int c = 0; c < nchunks; ++c) {
            unionOf(static_cast<int>((long long)nblocks * c / nchunks), static_cast<int>((long long)nblocks * (c + 1) / nchunks), Uc[c]);
            mb += 8.0 * double(Uc[c].size()) * double(Uc[c].size()) / 1.0e6;
        }
        if (mb <= 1024.0 || nchunks == 1) break;
        nchunks = std::max(1, nchunks / 2);
    }
    const Eigen::MatrixXd Pd(P);
    std::vector<Eigen::MatrixXd> Vc(nchunks);
    std::vector<double> exc_c(nchunks, 0.0), ec_c(nchunks, 0.0), ne_c(nchunks, 0.0);
    std::vector<long> used_c(nchunks, 0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        const std::vector<int>& U = Uc[c];
        std::vector<int> pos(n, -1);
        for (size_t i = 0; i < U.size(); ++i) pos[U[i]] = static_cast<int>(i);
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(U.size(), U.size());
        double ex_sum = 0.0, ec_sum = 0.0, ne_sum = 0.0;
        long used = 0;
        qmgrid::AOBlock blk;
        std::vector<int> cols, lc;
        Matrix grad;
        Eigen::MatrixXd Pphi;
        for (int b = b0; b < b1; ++b) {
            int begin = 0, count = 0;
            if (screened) {
                begin = screen->begin[b];
                count = screen->count[b];
                const int off = screen->shell_off[b], nsh = screen->shell_off[b + 1] - off;
                if (nsh == 0) continue;  // no function reaches the block: rho = 0 there
                ao.evaluateShells(grid.points, begin, count, 1, screen->shells.data() + off, nsh, blk, cols);
            } else {
                begin = b * block;
                count = std::min(block, npts - begin);
                ao.evaluate(grid.points, begin, count, 1, blk);
                cols.resize(n);
                for (int m = 0; m < n; ++m) cols[m] = m;
            }
            const Eigen::MatrixXd phi(blk.phi);  // count x ns, column-major
            const Eigen::MatrixXd d[3] = { Eigen::MatrixXd(blk.dx), Eigen::MatrixXd(blk.dy), Eigen::MatrixXd(blk.dz) };
            const Eigen::MatrixXd Ps = gatherP(Pd, cols);
            Vector rho, sigma, tau = Vector::Zero(count);
            densitySubset(phi, d, Ps, rho, &grad, Pphi);
            sigma = grad.rowwise().squaredNorm();
            // tau_q = 1/2 sum_k sum_mn P_mn d_k phi_qm d_k phi_qn
            for (int k = 0; k < 3; ++k) {
                const Eigen::MatrixXd dP = d[k] * Ps;
                for (int q = 0; q < count; ++q) tau(q) += 0.5 * d[k].row(q).dot(dP.row(q));
            }
            Vector ex, ec, vrho, vsigma, vtau;
            evaluateXCMeta(f, rho, sigma, tau, rho_threshold, ex, ec, vrho, vsigma, vtau);
            Eigen::VectorXd wv(count), wt(count), gx(count), gy(count), gz(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                wt(q) = 0.5 * w * vtau(q);
                const double w2s = 2.0 * w * vsigma(q);
                gx(q) = w2s * grad(q, 0);
                gy(q) = w2s * grad(q, 1);
                gz(q) = w2s * grad(q, 2);
                ne_sum += w * rho(q);
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    ex_sum += w * rho(q) * ex(q);
                    ec_sum += w * rho(q) * ec(q);
                    ++used;
                    any = true;
                }
            }
            if (!any) continue;
            Eigen::MatrixXd Y = (0.5 * wv).asDiagonal() * phi;
            Y.noalias() += gx.asDiagonal() * d[0];
            Y.noalias() += gy.asDiagonal() * d[1];
            Y.noalias() += gz.asDiagonal() * d[2];
            Eigen::MatrixXd X = phi.transpose() * Y;
            X += X.transpose().eval();
            for (int k = 0; k < 3; ++k) X.noalias() += d[k].transpose() * (wt.asDiagonal() * d[k]);
            const int ns = static_cast<int>(cols.size());
            lc.resize(ns);
            for (int i = 0; i < ns; ++i) lc[i] = pos[cols[i]];
            for (int j = 0; j < ns; ++j)
                for (int i = 0; i < ns; ++i) V(lc[i], lc[j]) += X(i, j);
        }
        Vc[c] = std::move(V);
        exc_c[c] = ex_sum;
        ec_c[c] = ec_sum;
        ne_c[c] = ne_sum;
        used_c[c] = used;
    }
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
    out.ex = out.ec = out.nelec = 0.0;
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {  // chunk order: fixed summation order
        const std::vector<int>& U = Uc[c];
        const int nu = static_cast<int>(U.size());
        for (int j = 0; j < nu; ++j)
            for (int i = 0; i < nu; ++i) V(U[i], U[j]) += Vc[c](i, j);
        out.ex += exc_c[c];
        out.ec += ec_c[c];
        out.nelec += ne_c[c];
        out.points_used += used_c[c];
    }
    out.exc = out.ex + out.ec;
    out.V = Matrix(0.5 * (V + V.transpose()));
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

/**
 * buildXCGradient for a meta-GGA functional (WP11 stage 2). Blocks as buildVxcMeta (screen blocks with the kept
 * functions, else uniform `block`-point blocks with all functions); AO derivatives up to second order.
 *
 * E_xc = sum_g w_g e(rho, sigma, tau); with phi_m on atom A, d phi_m / d A_k = - d_k phi_m:
 *   d rho / d A_k     = -2 sum_{m on A} d_k phi_m (P phi)_m
 *   d sigma / d A_k   = 2 grad rho . d grad rho / d A_k   (GGA part, unchanged from buildXCGradient)
 *   d tau / d A_k     = - sum_{m on A} sum_j d_k d_j phi_m (P d_j phi)_m        (tau = 1/2 sum_j sum_mn P d_j phi d_j phi)
 * so per point and function (T_k as in buildXCGradient, dE/dA_k = -2 sum_g sum_{m on A} T_k):
 *   T_k = d_k phi (Z P) + sum_j v_j d_k d_j phi (P phi) + sum_j (w v_tau / 2) d_k d_j phi (P d_j phi),
 *   Z = w v_rho phi + sum_j v_j d_j phi,  v_j = 2 w v_sigma d_j rho.
 * Only second AO derivatives are needed (tau holds first derivatives, its nuclear derivative one more).
 * The point part (grid response) uses the same translation identity as the GGA path; the Becke weight derivative
 * needs only e_p = rho eps_xc (`edens`), which is the meta-GGA integrand here.
 * Derivatives are taken at libxc's clamped point (tau floor, sigma <= 8 rho tau), as in V_xc: where a clamp is active
 * the gradient is not the derivative of the clamped energy (for P positive semidefinite tau >= tau_W holds exactly,
 * the clamp only acts at rounding level). Fixed chunk count (32): 1 vs N threads bitwise.
 * Claude Generated (Oct 2026, WP11 stage 2).
 */
void buildXCGradientMeta(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                         const XCFunctional& f, double rho_threshold, int threads, const std::vector<int>& bf_atom,
                         int natoms, XCGradient& out, int block, const qmgrid::BlockScreen* screen)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    const bool screened = screen && screen->active();
    if (block < 1) block = 128;
    const int nblocks = screened ? screen->nblocks() : (npts + block - 1) / block;
    const int nchunks = std::max(1, std::min(32, nblocks));
    const Eigen::MatrixXd Pd(P);
    std::vector<Matrix> gb(nchunks);
    std::vector<long> used_c(nchunks, 0);
    out.edens = Vector::Zero(npts);
    Matrix spoint = Matrix::Zero(npts, 3);  // per point: sum over all functions of the basis part

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        long used = 0;
        qmgrid::AOBlock blk;
        std::vector<int> cols;
        Matrix grad;
        Eigen::MatrixXd PhiP;
        for (int b = b0; b < b1; ++b) {
            int begin = 0, count = 0;
            if (screened) {
                begin = screen->begin[b];
                count = screen->count[b];
                const int off = screen->shell_off[b], nsh = screen->shell_off[b + 1] - off;
                if (nsh == 0) continue;  // no function reaches the block: rho = 0 there
                ao.evaluateShells(grid.points, begin, count, 2, screen->shells.data() + off, nsh, blk, cols);
            } else {
                begin = b * block;
                count = std::min(block, npts - begin);
                ao.evaluate(grid.points, begin, count, 2, blk);
                cols.resize(n);
                for (int m = 0; m < n; ++m) cols[m] = m;
            }
            const Eigen::MatrixXd phi(blk.phi);
            const Eigen::MatrixXd d[3] = { Eigen::MatrixXd(blk.dx), Eigen::MatrixXd(blk.dy), Eigen::MatrixXd(blk.dz) };
            const Eigen::MatrixXd Ps = gatherP(Pd, cols);
            Vector rho, sigma, tau = Vector::Zero(count);
            densitySubset(phi, d, Ps, rho, &grad, PhiP);  // PhiP = phi P_SS
            sigma = grad.rowwise().squaredNorm();
            Eigen::MatrixXd dP[3];  // (d_j phi) P_SS
            for (int k = 0; k < 3; ++k) {
                dP[k].noalias() = d[k] * Ps;
                for (int q = 0; q < count; ++q) tau(q) += 0.5 * d[k].row(q).dot(dP[k].row(q));
            }
            Vector ex, ec, vrho, vsigma, vtau;
            evaluateXCMeta(f, rho, sigma, tau, rho_threshold, ex, ec, vrho, vsigma, vtau);
            Eigen::VectorXd wv(count), vx(count), vy(count), vz(count), wt(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                vx(q) = vy(q) = vz(q) = wt(q) = 0.0;
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    out.edens(begin + q) = rho(q) * (ex(q) + ec(q));
                    ++used;
                    any = true;
                    const double w2s = 2.0 * w * vsigma(q);
                    vx(q) = w2s * grad(q, 0);
                    vy(q) = w2s * grad(q, 1);
                    vz(q) = w2s * grad(q, 2);
                    wt(q) = 0.5 * w * vtau(q);
                }
            }
            if (!any) continue;
            Eigen::MatrixXd Zm = wv.asDiagonal() * phi;
            Zm.noalias() += vx.asDiagonal() * d[0];
            Zm.noalias() += vy.asDiagonal() * d[1];
            Zm.noalias() += vz.asDiagonal() * d[2];
            const Eigen::MatrixXd ZP = Zm * Ps;
            Eigen::MatrixXd T[3] = { d[0].cwiseProduct(ZP), d[1].cwiseProduct(ZP), d[2].cwiseProduct(ZP) };
            const Eigen::MatrixXd xx(blk.dxx), xy(blk.dxy), xz(blk.dxz), yy(blk.dyy), yz(blk.dyz), zz(blk.dzz);
            // GGA: sum_j v_j d_k d_j phi_m (P phi)_m
            T[0].noalias() += ((vx.asDiagonal() * xx) + (vy.asDiagonal() * xy) + (vz.asDiagonal() * xz)).cwiseProduct(PhiP);
            T[1].noalias() += ((vx.asDiagonal() * xy) + (vy.asDiagonal() * yy) + (vz.asDiagonal() * yz)).cwiseProduct(PhiP);
            T[2].noalias() += ((vx.asDiagonal() * xz) + (vy.asDiagonal() * yz) + (vz.asDiagonal() * zz)).cwiseProduct(PhiP);
            // tau: sum_j (w v_tau / 2) d_k d_j phi_m (P d_j phi)_m
            T[0].noalias() += (wt.asDiagonal() * (xx.cwiseProduct(dP[0]) + xy.cwiseProduct(dP[1]) + xz.cwiseProduct(dP[2])));
            T[1].noalias() += (wt.asDiagonal() * (xy.cwiseProduct(dP[0]) + yy.cwiseProduct(dP[1]) + yz.cwiseProduct(dP[2])));
            T[2].noalias() += (wt.asDiagonal() * (xz.cwiseProduct(dP[0]) + yz.cwiseProduct(dP[1]) + zz.cwiseProduct(dP[2])));
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
    for (int p = 0; p < npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

void buildVxc(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
              const XCFunctional& f, double rho_threshold, int threads, XCResult& out, int block,
              const qmgrid::BlockScreen* screen)
{
    if (f.needsTau()) {  // WP11 meta-GGA: own function, the LDA/GGA paths below are unchanged
        buildVxcMeta(grid, ao, P, f, rho_threshold, threads, out, block, screen);
        return;
    }
    if (screen && screen->active()) {  // AO screening (Oct 2026); the unscreened path below is unchanged
        buildVxcScreened(grid, ao, P, f, rho_threshold, threads, out, *screen);
        return;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    // Fixed chunk count from the basis size only (memory budget ~128 MB for the per-chunk matrices), so the
    // summation order, and with it every bit of the result, is the same for any thread count.
    const double mat_mb = 8.0 * double(n) * n / 1.0e6;
    const int nchunks = std::max(1, std::min({ 32, nblocks, static_cast<int>(128.0 / std::max(mat_mb, 1e-3)) }));

    std::vector<Eigen::MatrixXd> Vc(nchunks);
    std::vector<double> exc_c(nchunks, 0.0), ec_c(nchunks, 0.0), ne_c(nchunks, 0.0);
    std::vector<long> used_c(nchunks, 0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
        double ex_sum = 0.0, ec_sum = 0.0, ne_sum = 0.0;
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix grad;  // grad rho on the block (GGA only)
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 1 : 0, blk);
            Vector rho, sigma;
            qmgrid::densityOnBlock(blk, P, rho, gga ? &grad : nullptr);
            if (gga) sigma = grad.rowwise().squaredNorm();
            Vector ex, ec, vrho, vsigma;
            evaluateXC(f, rho, sigma, rho_threshold, ex, ec, vrho, vsigma);
            // weighted potential per point; points below the density threshold carry vrho = 0
            Eigen::VectorXd wv(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                ne_sum += w * rho(q);
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    ex_sum += w * rho(q) * ex(q);
                    ec_sum += w * rho(q) * ec(q);
                    ++used;
                    any = true;
                }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi);  // count x n, column-major
            if (!gga) {
                // LDA: V += phi^T diag(w v) phi  (one GEMM per block)
                const Eigen::MatrixXd wphi = wv.asDiagonal() * phi;
                V.noalias() += phi.transpose() * wphi;
                continue;
            }
            // GGA (WP6): Y = (w v_rho / 2) phi + 2 w v_sigma (grad rho . grad phi), V += phi^T Y + Y^T phi
            //            = sum_g w [v_rho phi_m phi_n + 2 v_sigma grad rho . (phi_m grad phi_n + grad phi_m phi_n)]
            Eigen::VectorXd gx(count), gy(count), gz(count);
            for (int q = 0; q < count; ++q) {
                const double w2s = 2.0 * grid.weights(begin + q) * vsigma(q);
                gx(q) = w2s * grad(q, 0);
                gy(q) = w2s * grad(q, 1);
                gz(q) = w2s * grad(q, 2);
            }
            Eigen::MatrixXd Y = (0.5 * wv).asDiagonal() * phi;
            Y.noalias() += gx.asDiagonal() * Eigen::MatrixXd(blk.dx);
            Y.noalias() += gy.asDiagonal() * Eigen::MatrixXd(blk.dy);
            Y.noalias() += gz.asDiagonal() * Eigen::MatrixXd(blk.dz);
            const Eigen::MatrixXd X = phi.transpose() * Y;
            V.noalias() += X;
            V.noalias() += X.transpose();
        }
        Vc[c] = std::move(V);
        exc_c[c] = ex_sum;
        ec_c[c] = ec_sum;
        ne_c[c] = ne_sum;
        used_c[c] = used;
    }

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
    out.ex = out.ec = out.nelec = 0.0;
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {
        V += Vc[c];
        out.ex += exc_c[c];
        out.ec += ec_c[c];
        out.nelec += ne_c[c];
        out.points_used += used_c[c];
    }
    out.exc = out.ex + out.ec;
    out.V = Matrix(0.5 * (V + V.transpose()));
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// WP8 (Claude Generated, Oct 2026): nuclear gradient of E_xc, derivation in qm_xc.h (buildXCGradient).
void buildXCGradient(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                     const XCFunctional& f, double rho_threshold, int threads, const std::vector<int>& bf_atom,
                     int natoms, XCGradient& out, int block, const qmgrid::BlockScreen* screen)
{
    if (f.needsTau()) {  // WP11 stage 2: meta-GGA in its own function, the LDA/GGA paths below are unchanged
        buildXCGradientMeta(grid, ao, P, f, rho_threshold, threads, bf_atom, natoms, out, block, screen);
        return;
    }
    if (screen && screen->active()) {  // AO screening (Oct 2026); the unscreened path below is unchanged
        buildXCGradientScreened(grid, ao, P, f, rho_threshold, threads, bf_atom, natoms, out, *screen);
        return;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    const int nchunks = std::max(1, std::min(32, nblocks));
    const Eigen::MatrixXd Pd(P);
    std::vector<Matrix> gb(nchunks);
    std::vector<long> used_c(nchunks, 0);
    out.edens = Vector::Zero(npts);
    Matrix spoint = Matrix::Zero(npts, 3);  // per point: sum over all functions of the basis part

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix grad;
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 2 : 1, blk);
            Vector rho, sigma;
            qmgrid::densityOnBlock(blk, P, rho, gga ? &grad : nullptr);
            if (gga) sigma = grad.rowwise().squaredNorm();
            Vector ex, ec, vrho, vsigma;
            evaluateXC(f, rho, sigma, rho_threshold, ex, ec, vrho, vsigma);
            Eigen::VectorXd wv(count), vx(count), vy(count), vz(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                vx(q) = vy(q) = vz(q) = 0.0;
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    out.edens(begin + q) = rho(q) * (ex(q) + ec(q));
                    ++used;
                    any = true;
                    if (gga) {
                        const double w2s = 2.0 * w * vsigma(q);
                        vx(q) = w2s * grad(q, 0);
                        vy(q) = w2s * grad(q, 1);
                        vz(q) = w2s * grad(q, 2);
                    }
                }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi), dx(blk.dx), dy(blk.dy), dz(blk.dz);
            Eigen::MatrixXd Zm = wv.asDiagonal() * phi;
            if (gga) {
                Zm.noalias() += vx.asDiagonal() * dx;
                Zm.noalias() += vy.asDiagonal() * dy;
                Zm.noalias() += vz.asDiagonal() * dz;
            }
            const Eigen::MatrixXd ZP = Zm * Pd;  // (Z P)_pm, P symmetric
            Eigen::MatrixXd T[3] = { dx.cwiseProduct(ZP), dy.cwiseProduct(ZP), dz.cwiseProduct(ZP) };
            if (gga) {
                const Eigen::MatrixXd PhiP = phi * Pd;
                const Eigen::MatrixXd xx(blk.dxx), xy(blk.dxy), xz(blk.dxz), yy(blk.dyy), yz(blk.dyz), zz(blk.dzz);
                // sum_j v_j d_k d_j phi_m (P phi)_m
                T[0].noalias() += ((vx.asDiagonal() * xx) + (vy.asDiagonal() * xy) + (vz.asDiagonal() * xz)).cwiseProduct(PhiP);
                T[1].noalias() += ((vx.asDiagonal() * xy) + (vy.asDiagonal() * yy) + (vz.asDiagonal() * yz)).cwiseProduct(PhiP);
                T[2].noalias() += ((vx.asDiagonal() * xz) + (vy.asDiagonal() * yz) + (vz.asDiagonal() * zz)).cwiseProduct(PhiP);
            }
            for (int k = 0; k < 3; ++k) {
                const Eigen::VectorXd colsum = T[k].colwise().sum().transpose();  // per function
                for (int m = 0; m < n; ++m) g(bf_atom[m], k) += -2.0 * colsum(m);
                const Eigen::VectorXd rowsum = T[k].rowwise().sum();               // per point
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
    // points of atom B move with B: d e(r_p)/d R_B (point part) = d_r e = - (sum of all basis parts of the point)
    out.points = Matrix::Zero(natoms, 3);
    for (int p = 0; p < npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace qmxc
