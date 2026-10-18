/*
 * test_native_eigensolver — validate the self-contained symmetric eigensolver
 * (curcuma::eigsolver::solveSymmetric, the `-eigensolver native` backend = Householder
 * tridiagonalization + Cuppen divide-and-conquer) against Eigen's SelfAdjointEigenSolver.
 *
 * Covers the cases the divide-and-conquer must get right:
 *   - random well-separated spectra (deep recursion, secular accuracy)
 *   - exact eigenvalue multiplicities and identity-like spectra (deflation)
 *   - tight clusters (deflation + Löwner/Gu-Eisenstat orthogonality)
 * Checks eigenvalues, reconstruction V·Λ·Vᵀ ≈ A, and orthonormality Vᵀ·V ≈ I.
 *
 * Exit code 0 = pass, 1 = fail (ctest).
 * Copyright (C) 2019 - 2026 Conrad Hübler. Claude Generated. GPL-3.0.
 */

#include <Eigen/Dense>
#include <cstdio>
#include <random>

#include "src/core/energy_calculators/qm_methods/native_eigensolver.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

static std::mt19937 rng(12345);

static MatrixXd randOrth(int n)
{
    MatrixXd M(n, n);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) M(i, j) = nd(rng);
    return Eigen::HouseholderQR<MatrixXd>(M).householderQ();
}

static MatrixXd fromSpectrum(const VectorXd& lam)
{
    const MatrixXd Q = randOrth(static_cast<int>(lam.size()));
    return Q * lam.asDiagonal() * Q.transpose();
}

static bool check(const char* name, const MatrixXd& A)
{
    const int n = static_cast<int>(A.rows());
    VectorXd ev;
    MatrixXd vc;
    if (!curcuma::eigsolver::solveSymmetric(A, ev, vc)) {
        std::printf("%-20s n=%4d  FAILED (solve)\n", name, n);
        return false;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    const double scale = std::max(1.0, A.cwiseAbs().maxCoeff());
    const double de    = (ev - es.eigenvalues()).cwiseAbs().maxCoeff() / scale;
    const double recon = (vc * ev.asDiagonal() * vc.transpose() - A).cwiseAbs().maxCoeff() / scale;
    const double orth  = (vc.transpose() * vc - MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff();

    const bool ok = (de < 1e-10) && (recon < 1e-11) && (orth < 1e-10);
    std::printf("%-20s n=%4d  d(eval)=%.2e recon=%.2e orth=%.2e  %s\n",
                name, n, de, recon, orth, ok ? "PASS" : "FAIL");
    return ok;
}

int main()
{
    bool ok = true;
    std::normal_distribution<double> nd(0.0, 1.0);

    for (int n : {1, 2, 3, 5, 13, 40, 137, 300, 558}) {
        MatrixXd M(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) M(i, j) = nd(rng);
        ok &= check("random", 0.5 * (M + M.transpose()));
    }
    { VectorXd l(8);  l << 1, 1, 1, 2, 2, 5, 5, 9;                 ok &= check("degenerate-mult", fromSpectrum(l)); }
    { VectorXd l(50); for (int i = 0; i < 50; ++i) l(i) = (i < 25) ? 3.0 : 7.0;
                                                                   ok &= check("two-clusters", fromSpectrum(l)); }
    { VectorXd l = VectorXd::Ones(60);                             ok &= check("identity-like", fromSpectrum(l)); }
    { VectorXd l(40); for (int i = 0; i < 40; ++i) l(i) = (i < 20) ? 1.0 + i * 1e-11 : 5.0 + (i - 20) * 1e-10;
                                                                   ok &= check("tight-clusters", fromSpectrum(l)); }
    { VectorXd l(200); for (int i = 0; i < 200; ++i) l(i) = double(i / 20);
                                                                   ok &= check("blocky-200", fromSpectrum(l)); }

    std::printf("\nOverall: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
