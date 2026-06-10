/*
 * test_purification — validate the 0 K density-matrix purification kernel
 * (curcuma::eigsolver::purifyDensity, the GEMM-only `-eigensolver purify` density path)
 * against the occupied-subspace projector from Eigen's SelfAdjointEigenSolver.
 *
 * Checks, on gapped symmetric matrices of various sizes and occupations:
 *   - idempotency  P² ≈ P
 *   - electron count  Tr(P) = nocc
 *   - P matches the reference projector onto the nocc lowest eigenstates
 *   - energy-weighted density  W = P·A·P  matches  Σ_{i<nocc} ε_i c_i c_iᵀ
 *
 * Exit code 0 = pass, 1 = fail (ctest).
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. Claude Generated. GPL-3.0.
 */

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

#include "src/core/energy_calculators/qm_methods/native_eigensolver.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

static std::mt19937 rng(2026);

static MatrixXd randOrth(int n)
{
    MatrixXd M(n, n);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) M(i, j) = nd(rng);
    return Eigen::HouseholderQR<MatrixXd>(M).householderQ();
}

// Symmetric matrix with a clear gap between eigenstate nocc-1 (occupied) and nocc (virtual).
static MatrixXd gappedMatrix(int n, int nocc, double gap)
{
    VectorXd lam(n);
    std::uniform_real_distribution<double> lo(-2.0, -0.5), hi(0.5, 2.0);
    for (int i = 0; i < nocc; ++i) lam(i) = lo(rng);            // occupied band (low)
    for (int i = nocc; i < n; ++i) lam(i) = hi(rng) + gap;      // virtual band (high, separated)
    std::sort(lam.data(), lam.data() + n);
    const MatrixXd Q = randOrth(n);
    return Q * lam.asDiagonal() * Q.transpose();
}

static bool check(const char* name, int n, int nocc, double gap)
{
    const MatrixXd A = gappedMatrix(n, nocc, gap);
    MatrixXd P, W;
    if (!curcuma::eigsolver::purifyDensity(A, nocc, P, W)) {
        std::printf("%-16s n=%4d nocc=%4d  FAILED (no convergence)\n", name, n, nocc);
        return false;
    }
    // Reference: projector + energy-weighted density from the nocc lowest eigenstates.
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    const MatrixXd Cocc = es.eigenvectors().leftCols(nocc);
    const MatrixXd Pref = Cocc * Cocc.transpose();
    const VectorXd eocc = es.eigenvalues().head(nocc);
    const MatrixXd Wref = Cocc * eocc.asDiagonal() * Cocc.transpose();

    const double idem  = (P * P - P).cwiseAbs().maxCoeff();
    const double trerr = std::fabs(P.trace() - double(nocc));
    const double dP    = (P - Pref).cwiseAbs().maxCoeff();
    const double dW    = (W - Wref).cwiseAbs().maxCoeff();

    const bool ok = (idem < 1e-8) && (trerr < 1e-8) && (dP < 1e-7) && (dW < 1e-7);
    std::printf("%-16s n=%4d nocc=%4d  idem=%.2e dTr=%.2e dP=%.2e dW=%.2e  %s\n",
                name, n, nocc, idem, trerr, dP, dW, ok ? "PASS" : "FAIL");
    return ok;
}

int main()
{
    bool ok = true;
    ok &= check("small",      20,   8, 0.5);
    ok &= check("mid",       137,  60, 0.3);
    ok &= check("half",      200, 100, 0.4);
    ok &= check("sparse-occ", 300,  40, 0.5);
    ok &= check("dense-occ",  300, 250, 0.5);
    ok &= check("gfn-like",   558, 280, 0.2);

    std::printf("\nOverall: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
