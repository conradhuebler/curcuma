/*
 * Unit tests for the stand-alone Ewald sum (work package B0)
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Nov 2026, B0). Modes:
 *   test_ewald          correctness tests (a)-(g), see docs/WP-PERIODIC-NONBONDED.md, section B0
 *   test_ewald perf     timing scan (h), N = 1000, 3000, 6000 at water atom density
 * Units: Bohr, e, Eh. Exit status 0 = all checks passed.
 *
 * Independent references (not Ewald): (a) Evjen summation for NaCl, (g) neutral-cube sum for the
 * simple cubic lattice in a uniform background. Test (c) compares against a second, deliberately
 * naive Ewald implementation in long double: that checks self-consistency, not the physics.
 */

#include "src/core/energy_calculators/ff_methods/ewald.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using curcuma::ewald::EwaldParameters;
using curcuma::ewald::EwaldSum;
using curcuma::ewald::chooseEwaldParameters;

namespace {

int g_failures = 0;

void check(bool ok, const char* what, double value, double limit)
{
    std::printf("  [%s] %-58s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what, value, limit);
    if (!ok)
        ++g_failures;
}

PBCUtils::OrthoCell makeCell(double lx, double ly, double lz)
{
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m(0, 0) = lx;
    m(1, 1) = ly;
    m(2, 2) = lz;
    std::string err;
    PBCUtils::OrthoCell c = PBCUtils::makeOrthoCell(m, err);
    if (!c.active) {
        std::printf("makeOrthoCell failed: %s\n", err.c_str());
        std::exit(2);
    }
    return c;
}

double relDiff(double a, double b) { return std::abs(a - b) / std::max(std::abs(a), std::abs(b)); }

double relMatDiff(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b)
{
    return (a - b).cwiseAbs().maxCoeff() / std::max(a.cwiseAbs().maxCoeff(), b.cwiseAbs().maxCoeff());
}

/// Random positions in the box with a minimum-image distance of at least dmin, neutral charges in [-1, 1].
void randomSystem(int N, const PBCUtils::OrthoCell& cell, double dmin, unsigned seed,
    Eigen::MatrixXd& pos, Eigen::VectorXd& q)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    pos.resize(N, 3);
    q.resize(N);
    int placed = 0;
    while (placed < N) {
        double p[3] = { u(rng) * cell.L[0], u(rng) * cell.L[1], u(rng) * cell.L[2] };
        bool ok = true;
        for (int j = 0; j < placed && ok; ++j) {
            double d[3] = { p[0] - pos(j, 0), p[1] - pos(j, 1), p[2] - pos(j, 2) };
            PBCUtils::applyMinimumImageOrtho(cell, d[0], d[1], d[2]);
            ok = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) >= dmin * dmin;
        }
        if (!ok)
            continue;
        for (int a = 0; a < 3; ++a)
            pos(placed, a) = p[a];
        q(placed) = 2.0 * u(rng) - 1.0;
        ++placed;
    }
    q.array() -= q.mean();
}

// ---------------------------------------------------------------------------------------------
// (a) NaCl: Ewald and Evjen
// ---------------------------------------------------------------------------------------------

/// Evjen (1932) sum for the rock-salt lattice: neutral cube -n..n around one ion, face/edge/corner
/// ions weighted 1/2, 1/4, 1/8. Returns the Madelung constant M (positive), unit spacing r0 = 1.
double evjenNaCl(int n)
{
    long double sum = 0.0L; // 4.2 million terms at n = 80: a double accumulator loses ~1e-9
    for (int i = -n; i <= n; ++i)
        for (int j = -n; j <= n; ++j)
            for (int k = -n; k <= n; ++k) {
                if (i == 0 && j == 0 && k == 0)
                    continue;
                double w = 1.0;
                if (std::abs(i) == n)
                    w *= 0.5;
                if (std::abs(j) == n)
                    w *= 0.5;
                if (std::abs(k) == n)
                    w *= 0.5;
                const double sign = ((i + j + k) % 2 == 0) ? 1.0 : -1.0;
                sum += sign * w / std::sqrt(double(i * i + j * j + k * k));
            }
    return -(double)sum;
}

void testNaCl()
{
    std::printf("(a) NaCl rock salt, conventional cubic cell, 8 ions\n");
    const double r0 = 5.0; // nearest-neighbour distance, Bohr (any value; E scales as 1/r0)
    const double L = 2.0 * r0;
    PBCUtils::OrthoCell cell = makeCell(L, L, L);
    Eigen::MatrixXd pos(8, 3);
    Eigen::VectorXd q(8);
    const double fcc[4][3] = { { 0, 0, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 0.5, 0.5, 0 } };
    for (int s = 0; s < 4; ++s) {
        for (int a = 0; a < 3; ++a) {
            pos(s, a) = fcc[s][a] * L;
            pos(s + 4, a) = fcc[s][a] * L + (a == 0 ? 0.5 * L : 0.0);
        }
        q(s) = 1.0;
        q(s + 4) = -1.0;
    }

    const double m40 = evjenNaCl(40), m80 = evjenNaCl(80);
    const double m_evjen = m80 + (m80 - m40) / 15.0; // error ~ n^-4 (measured: x16 per doubling)
    std::printf("  Evjen n=40 %.10f  n=80 %.10f  Richardson %.10f\n", m40, m80, m_evjen);

    const double M_lit6 = 1.747565;          // 6-digit literature value
    const double M_lit = 1.747564594633;     // longer value quoted in the web search result (Borwein et al.)
    double m_ew[3];
    const double rcuts[3] = { 0.5 * L, 0.4 * L, 0.3 * L };
    for (int t = 0; t < 3; ++t) {
        EwaldSum ew(cell, chooseEwaldParameters(cell, 1e-14, rcuts[t]));
        m_ew[t] = -ew.energy(pos, q) / 4.0 * r0; // E per ion pair = -M q^2 / r0
        std::printf("  Ewald r_cut %.2f alpha %.4f N_k %d: M = %.12f\n", rcuts[t], ew.parameters().alpha,
            ew.numberOfKVectors(), m_ew[t]);
    }
    check(relDiff(m_ew[0], m_evjen) < 1e-7, "Ewald vs Evjen (independent) rel", relDiff(m_ew[0], m_evjen), 1e-7);
    check(relDiff(m_ew[0], M_lit) < 1e-7, "Ewald vs 1.747564594633 rel", relDiff(m_ew[0], M_lit), 1e-7);
    check(std::abs(m_ew[0] - M_lit6) <= 5e-7, "Ewald vs 1.747565 (6-digit rounding) abs", std::abs(m_ew[0] - M_lit6), 5e-7);
    const double spread = std::max(relDiff(m_ew[0], m_ew[1]), relDiff(m_ew[0], m_ew[2]));
    check(spread < 1e-12, "NaCl alpha spread rel", spread, 1e-12);
}

// ---------------------------------------------------------------------------------------------
// (b) alpha independence
// ---------------------------------------------------------------------------------------------

void testAlphaIndependence(const Eigen::MatrixXd& pos, const Eigen::VectorXd& q, const PBCUtils::OrthoCell& cell)
{
    std::printf("(b) alpha independence, N=%d neutral, L=%.0f Bohr\n", int(q.size()), cell.L[0]);
    const double rcuts[3] = { 15.0, 12.0, 9.0 }; // 9 Bohr: 3 cells per axis, linked-cell path
    curcuma::ewald::EwaldResult r[3];
    for (int t = 0; t < 3; ++t) {
        EwaldSum ew(cell, chooseEwaldParameters(cell, 1e-12, rcuts[t]));
        r[t] = ew.evaluate(pos, q, true, true);
        std::printf("  r_cut %5.1f alpha %.4f k_max %.3f N_k %6d  E = %.14f (real %.6f recip %.6f self %.6f)\n",
            rcuts[t], ew.parameters().alpha, ew.parameters().k_max, ew.numberOfKVectors(), r[t].terms.total(),
            r[t].terms.real, r[t].terms.reciprocal, r[t].terms.self);
    }
    double dE = 0, dG = 0, dP = 0;
    for (int t = 1; t < 3; ++t) {
        dE = std::max(dE, relDiff(r[0].terms.total(), r[t].terms.total()));
        dG = std::max(dG, relMatDiff(r[0].gradient, r[t].gradient));
        dP = std::max(dP, relMatDiff(r[0].potentials, r[t].potentials));
    }
    check(dE < 1e-8, "energy spread over 3 alphas rel", dE, 1e-8);
    check(dG < 1e-8, "gradient spread (max abs / max |g|)", dG, 1e-8);
    check(dP < 1e-8, "potential spread (max abs / max |phi|)", dP, 1e-8);
}

// ---------------------------------------------------------------------------------------------
// (c) naive long-double Ewald with very large parameters (self-consistency reference)
// ---------------------------------------------------------------------------------------------

void naiveEwaldLD(const Eigen::MatrixXd& pos, const Eigen::VectorXd& q, const PBCUtils::OrthoCell& cell,
    long double alpha, int nimg, int nk, long double& E, std::vector<long double>& phi)
{
    const int N = int(q.size());
    const long double pi = 3.14159265358979323846264338327950288L;
    const long double V = (long double)cell.L[0] * cell.L[1] * cell.L[2];
    phi.assign(N, 0.0L);
    // real space over all images in a cube, every ordered pair (i, j, n) except i = j at n = 0
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int a = -nimg; a <= nimg; ++a)
                for (int b = -nimg; b <= nimg; ++b)
                    for (int c = -nimg; c <= nimg; ++c) {
                        if (i == j && a == 0 && b == 0 && c == 0)
                            continue;
                        long double d0 = (long double)pos(j, 0) - pos(i, 0) + a * (long double)cell.L[0];
                        long double d1 = (long double)pos(j, 1) - pos(i, 1) + b * (long double)cell.L[1];
                        long double d2 = (long double)pos(j, 2) - pos(i, 2) + c * (long double)cell.L[2];
                        long double r = std::sqrt(d0 * d0 + d1 * d1 + d2 * d2);
                        phi[i] += q(j) * std::erfc(alpha * r) / r;
                    }
    // reciprocal space over the full cube of n (no half-space, no recursion)
    for (int a = -nk; a <= nk; ++a)
        for (int b = -nk; b <= nk; ++b)
            for (int c = -nk; c <= nk; ++c) {
                if (a == 0 && b == 0 && c == 0)
                    continue;
                long double k[3] = { 2 * pi * a / cell.L[0], 2 * pi * b / cell.L[1], 2 * pi * c / cell.L[2] };
                long double k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                long double f = 4 * pi / V * std::exp(-k2 / (4 * alpha * alpha)) / k2;
                long double C = 0, S = 0;
                std::vector<long double> cs(N), sn(N);
                for (int i = 0; i < N; ++i) {
                    long double t = k[0] * pos(i, 0) + k[1] * pos(i, 1) + k[2] * pos(i, 2);
                    cs[i] = std::cos(t);
                    sn[i] = std::sin(t);
                    C += q(i) * cs[i];
                    S += q(i) * sn[i];
                }
                for (int i = 0; i < N; ++i)
                    phi[i] += f * (cs[i] * C + sn[i] * S);
            }
    long double Q = 0;
    for (int i = 0; i < N; ++i)
        Q += q(i);
    E = 0;
    for (int i = 0; i < N; ++i) {
        phi[i] += -2 * alpha / std::sqrt(pi) * q(i) - pi * Q / (alpha * alpha * V);
        E += 0.5L * q(i) * phi[i];
    }
}

void testLongDoubleReference()
{
    std::printf("(c) long-double naive Ewald reference (self-consistency, not independent)\n");
    PBCUtils::OrthoCell cell = makeCell(9.0, 10.0, 11.0);
    Eigen::MatrixXd pos;
    Eigen::VectorXd q;
    randomSystem(12, cell, 1.5, 7, pos, q);
    // alpha 0.7: real images |n| <= 4 leave erfc(0.7 * 31) ~ 1e-210; k up to 24: exp(-k^2/4a^2) < 1e-40
    long double Eref;
    std::vector<long double> phiref;
    naiveEwaldLD(pos, q, cell, 0.7L, 4, 24, Eref, phiref);
    std::printf("  reference E = %.16Lf\n", Eref);
    double last = 1.0;
    bool monotone = true;
    for (double eps : { 1e-6, 1e-8, 1e-10, 1e-12, 1e-14 }) {
        EwaldSum ew(cell, chooseEwaldParameters(cell, eps));
        auto r = ew.evaluate(pos, q, false, true);
        double dphi = 0, maxphi = 0;
        for (int i = 0; i < q.size(); ++i) {
            dphi = std::max(dphi, std::abs(r.potentials(i) - (double)phiref[i]));
            maxphi = std::max(maxphi, std::abs((double)phiref[i]));
        }
        const double dE = std::abs(r.terms.total() - (double)Eref) / std::abs((double)Eref);
        std::printf("  eps %.0e: |dE|/|E| = %.3e  max|dphi|/max|phi| = %.3e\n", eps, dE, dphi / maxphi);
        if (eps > 1e-12 && dE > last)
            monotone = false;
        last = dE;
        if (eps == 1e-12)
            check(dE < 1e-11 && dphi / maxphi < 1e-11, "eps 1e-12 vs long-double reference", std::max(dE, dphi / maxphi), 1e-11);
        if (eps == 1e-8)
            check(dE < 1e-7, "eps 1e-8 vs long-double reference (energy)", dE, 1e-7);
    }
    check(monotone, "error decreases from eps 1e-6 to 1e-12", monotone ? 0.0 : 1.0, 0.5);
}

// ---------------------------------------------------------------------------------------------
// (d) finite differences, (f) dense matrix
// ---------------------------------------------------------------------------------------------

void testFiniteDifferencesAndMatrix()
{
    std::printf("(d) finite differences and (f) dense matrix, N=16, cell 12 x 13 x 14 Bohr\n");
    PBCUtils::OrthoCell cell = makeCell(12.0, 13.0, 14.0);
    Eigen::MatrixXd pos;
    Eigen::VectorXd q;
    randomSystem(16, cell, 1.5, 11, pos, q);
    q(0) += 0.3; // charged cell: the background term enters every derivative
    EwaldSum ew(cell, chooseEwaldParameters(cell, 1e-13));
    auto r = ew.evaluate(pos, q, true, true);
    const double gmax = r.gradient.cwiseAbs().maxCoeff();

    double err[2];
    const double hs[2] = { 1e-4, 2e-4 };
    for (int t = 0; t < 2; ++t) {
        err[t] = 0;
        for (int i = 0; i < q.size(); ++i)
            for (int a = 0; a < 3; ++a) {
                Eigen::MatrixXd p = pos;
                p(i, a) += hs[t];
                const double ep = ew.energy(p, q);
                p(i, a) -= 2 * hs[t];
                const double em = ew.energy(p, q);
                err[t] = std::max(err[t], std::abs((ep - em) / (2 * hs[t]) - r.gradient(i, a)));
            }
        std::printf("  h = %.0e Bohr: max |FD - analytic| = %.3e Eh/Bohr (max |g| %.3e)\n", hs[t], err[t], gmax);
    }
    std::printf("  ratio err(2h)/err(h) = %.2f (4 for a quadratic truncation error)\n", err[1] / err[0]);
    check(err[0] / gmax < 1e-7, "gradient FD h=1e-4 (rel to max |g|)", err[0] / gmax, 1e-7);
    check(err[1] / err[0] > 3.0 && err[1] / err[0] < 5.0, "FD error quadratic in h (ratio in 3..5)", err[1] / err[0], 5.0);

    double perr = 0;
    const double hq = 1e-3; // E is quadratic in q: the central difference is exact up to rounding
    for (int i = 0; i < q.size(); ++i) {
        Eigen::VectorXd qq = q;
        qq(i) += hq;
        const double ep = ew.energy(pos, qq);
        qq(i) -= 2 * hq;
        const double em = ew.energy(pos, qq);
        perr = std::max(perr, std::abs((ep - em) / (2 * hq) - r.potentials(i)));
    }
    const double pmax = r.potentials.cwiseAbs().maxCoeff();
    check(perr / pmax < 1e-9, "potential vs FD dE/dq (rel to max |phi|)", perr / pmax, 1e-9);
    const double ehalf = 0.5 * q.dot(r.potentials);
    check(relDiff(ehalf, r.terms.total()) < 1e-13, "E == 1/2 sum q phi rel", relDiff(ehalf, r.terms.total()), 1e-13);

    Eigen::MatrixXd A = ew.matrix(pos);
    Eigen::VectorXd Aq = A * q;
    check(relMatDiff(Aq, r.potentials) < 1e-12, "dense A q == phi", relMatDiff(Aq, r.potentials), 1e-12);
    check((A - A.transpose()).cwiseAbs().maxCoeff() < 1e-14, "A symmetric (max abs)", (A - A.transpose()).cwiseAbs().maxCoeff(), 1e-14);
}

// ---------------------------------------------------------------------------------------------
// (e) translation invariance, (f) supercell
// ---------------------------------------------------------------------------------------------

void testTranslationAndSupercell(const Eigen::MatrixXd& pos, const Eigen::VectorXd& q, const PBCUtils::OrthoCell& cell)
{
    std::printf("(e) translation invariance and (f) 2x1x1 supercell, N=%d\n", int(q.size()));
    const EwaldParameters par = chooseEwaldParameters(cell, 1e-10, 12.0);
    EwaldSum ew(cell, par);
    auto r0 = ew.evaluate(pos, q, true, true);

    double worst = 0;
    const double shifts[3][3] = { { 3.7, -11.2, 0.4 }, { cell.L[0], 0, 0 }, { -cell.L[0], 2 * cell.L[1], cell.L[2] } };
    for (const auto& t : shifts) {
        Eigen::MatrixXd p = pos;
        for (int i = 0; i < p.rows(); ++i)
            for (int a = 0; a < 3; ++a)
                p(i, a) += t[a];
        auto r = ew.evaluate(p, q, true, false);
        worst = std::max(worst, std::max(relDiff(r.terms.total(), r0.terms.total()), relMatDiff(r.gradient, r0.gradient)));
    }
    check(worst < 1e-12, "translate all atoms (E rel, grad rel to max)", worst, 1e-12);

    worst = 0;
    for (int i : { 0, 5, 17 }) {
        Eigen::MatrixXd p = pos;
        p(i, 0) += cell.L[0];
        p(i, 2) -= cell.L[2];
        auto r = ew.evaluate(p, q, true, false);
        worst = std::max(worst, std::max(relDiff(r.terms.total(), r0.terms.total()), relMatDiff(r.gradient, r0.gradient)));
    }
    check(worst < 1e-12, "shift single atoms by L (E rel, grad rel to max)", worst, 1e-12);

    // Supercell: two copies along x, same alpha, r_cut and k_max (r_cut stays below L_y/2, L_z/2).
    const int N = int(q.size());
    PBCUtils::OrthoCell super = makeCell(2 * cell.L[0], cell.L[1], cell.L[2]);
    Eigen::MatrixXd p2(2 * N, 3);
    Eigen::VectorXd q2(2 * N);
    for (int i = 0; i < N; ++i) {
        p2.row(i) = pos.row(i);
        p2.row(i + N) = pos.row(i);
        p2(i + N, 0) += cell.L[0];
        q2(i) = q2(i + N) = q(i);
    }
    EwaldSum ew2(super, par);
    auto r2 = ew2.evaluate(p2, q2, true, true);
    const double dE = relDiff(r2.terms.total() / (2 * N), r0.terms.total() / N);
    std::printf("  E/N cell %.14f  supercell %.14f\n", r0.terms.total() / N, r2.terms.total() / (2 * N));
    check(dE < 1e-10, "supercell energy per atom rel", dE, 1e-10);
    check(relMatDiff(r2.gradient.topRows(N), r0.gradient) < 1e-10, "supercell gradient rel", relMatDiff(r2.gradient.topRows(N), r0.gradient), 1e-10);
    check(relMatDiff(r2.potentials.head(N), r0.potentials) < 1e-10, "supercell potentials rel", relMatDiff(r2.potentials.head(N), r0.potentials), 1e-10);
}

// ---------------------------------------------------------------------------------------------
// (g) single charge in a neutralising background (simple cubic Wigner lattice)
// ---------------------------------------------------------------------------------------------

/*! Independent value of xi from cubic cells: each cell = point charge at its centre plus uniform
 * background -1/L^3 filling the cube. Summing complete cells over -m..m in every direction gives
 *     phi_cubes(0) L = sum_{n != 0} 1/|n| - (2m+1)^2 C,  C = int_{unit cube} d^3r / r = 3 ln(2 + sqrt 3) - pi/2
 * (closed form; checked against 200-point Gauss-Legendre quadrature to 1e-15 while writing this test).
 * The cell has zero dipole and cubic symmetry, so the cell sum converges absolutely, but its mean
 * potential differs from Ewald's (k = 0 dropped) by -(2 pi/(3V)) int_cell rho r^2 = +pi/(6L)
 * (second moment of the background in the cube: -L^2/4). Hence xi = -(phi_cubes L - pi/6),
 * and E = 1/2 q phi_Ewald = -xi q^2/(2L). Error ~ m^-2 (measured), removed by Richardson.
 */
double wignerCubeSum(int m)
{
    const double C = 3.0 * std::log(2.0 + std::sqrt(3.0)) - 0.5 * M_PI;
    long double s = 0.0L; // cancellation against (2m+1)^2 C ~ 6e4: accumulate in long double
    for (int i = -m; i <= m; ++i)
        for (int j = -m; j <= m; ++j)
            for (int k = -m; k <= m; ++k)
                if (i || j || k)
                    s += 1.0L / std::sqrt((long double)(i * i + j * j + k * k));
    const double phiL = (double)s - double(2 * m + 1) * (2 * m + 1) * C;
    return -(phiL - M_PI / 6.0);
}

void testChargedCell()
{
    std::printf("(g) single charge with neutralising background, cubic L = 10 Bohr\n");
    const double x40 = wignerCubeSum(40), x80 = wignerCubeSum(80);
    const double xi_cube = x80 + (x80 - x40) / 3.0;
    std::printf("  cube sum m=40 %.10f  m=80 %.10f  Richardson %.10f\n", x40, x80, xi_cube);
    const double L = 10.0;
    PBCUtils::OrthoCell cell = makeCell(L, L, L);
    Eigen::MatrixXd pos(1, 3);
    pos << 1.3, 7.2, 4.4;
    double xi[3];
    for (int t = 0; t < 3; ++t) {
        const double rc[3] = { 5.0, 4.0, 3.0 };
        EwaldSum ew(cell, chooseEwaldParameters(cell, 1e-14, rc[t]));
        Eigen::VectorXd q(1);
        q(0) = -2.0; // q^2 scaling
        xi[t] = -2.0 * L * ew.energy(pos, q) / (q(0) * q(0));
        std::printf("  alpha %.4f: xi = %.12f\n", ew.parameters().alpha, xi[t]);
    }
    const double spread = std::max(relDiff(xi[0], xi[1]), relDiff(xi[0], xi[2]));
    // 1e-11: rounding of ~8e4 k terms at the smallest r_cut (measured 1.1e-12), required is 1e-8
    check(spread < 1e-11, "xi alpha spread rel", spread, 1e-11);
    check(relDiff(xi[0], xi_cube) < 2e-7, "Ewald xi vs independent cube sum rel", relDiff(xi[0], xi_cube), 2e-7);
    check(std::abs(xi[0] - 2.837297) <= 5e-7, "Ewald xi vs 2.837297 (memory value) abs", std::abs(xi[0] - 2.837297), 5e-7);
}

// ---------------------------------------------------------------------------------------------
// (h) timing
// ---------------------------------------------------------------------------------------------

void performanceScan()
{
    std::printf("(h) timing, tolerance 1e-6, one thread, energy + gradient + potentials\n");
    // water: 0.0334 molecules / A^3 = 0.1003 atoms / A^3; 1 A = 1.8897261 Bohr
    const double rho = 0.1003 / std::pow(1.8897261, 3);
    std::printf("  %6s %7s %7s %8s %8s %10s\n", "N", "L", "r_cut", "alpha", "N_k", "time/s");
    double best_t[3];
    const int Ns[3] = { 1000, 3000, 6000 };
    for (int s = 0; s < 3; ++s) {
        const int N = Ns[s];
        const double L = std::cbrt(N / rho);
        PBCUtils::OrthoCell cell = makeCell(L, L, L);
        Eigen::MatrixXd pos;
        Eigen::VectorXd q;
        randomSystem(N, cell, 0.0, 100 + s, pos, q);
        best_t[s] = 1e30;
        double best_rc = 0;
        for (int t = 0; t < 8; ++t) {
            const double rc = 8.0 + t * (0.5 * L - 8.0) / 7.0;
            EwaldSum ew(cell, chooseEwaldParameters(cell, 1e-6, rc));
            auto t0 = std::chrono::steady_clock::now();
            auto r = ew.evaluate(pos, q, true, true);
            const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::printf("  %6d %7.2f %7.2f %8.4f %8d %10.3f\n", N, L, rc, ew.parameters().alpha, ew.numberOfKVectors(), dt);
            if (dt < best_t[s]) {
                best_t[s] = dt;
                best_rc = rc;
            }
            (void)r;
        }
        if (s == 0) { // measured accuracy of the 1e-6 setting at the fastest r_cut
            EwaldSum fast(cell, chooseEwaldParameters(cell, 1e-6, best_rc));
            EwaldSum tight(cell, chooseEwaldParameters(cell, 1e-12));
            auto a = fast.evaluate(pos, q, true, false), b = tight.evaluate(pos, q, true, false);
            std::printf("  N=1000 accuracy at tolerance 1e-6, r_cut %.2f: E rel %.2e, gradient rel-to-max %.2e\n",
                best_rc, relDiff(a.terms.total(), b.terms.total()), relMatDiff(a.gradient, b.gradient));
        }
        std::printf("  N=%d fastest: r_cut %.2f Bohr, %.3f s\n", N, best_rc, best_t[s]);
    }
    std::printf("  scaling exponent 1000->6000: %.2f, 3000->6000: %.2f\n",
        std::log(best_t[2] / best_t[0]) / std::log(6.0), std::log(best_t[2] / best_t[1]) / std::log(2.0));
}

} // namespace

int main(int argc, char** argv)
{
    if (argc > 1 && std::strcmp(argv[1], "perf") == 0) {
        performanceScan();
        return 0;
    }
    testNaCl();

    PBCUtils::OrthoCell cell30 = makeCell(30.0, 30.0, 30.0);
    Eigen::MatrixXd pos;
    Eigen::VectorXd q;
    randomSystem(64, cell30, 1.5, 3, pos, q);
    testAlphaIndependence(pos, q, cell30);
    testLongDoubleReference();
    testFiniteDifferencesAndMatrix();
    testTranslationAndSupercell(pos, q, cell30);
    testChargedCell();

    std::printf("test_ewald: %d failure(s)\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
