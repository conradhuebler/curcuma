/*
 * test_native_eigensolver — validate the self-contained symmetric eigensolver
 * (curcuma::eigsolver::solveSymmetric, the `-eigensolver native` backend) against
 * Eigen's SelfAdjointEigenSolver on random symmetric matrices.
 *
 * Checks, over a range of sizes incl. the GFN minimal-basis range:
 *   - eigenvalues match (ascending) to a tight tolerance
 *   - reconstruction  V·diag(λ)·Vᵀ ≈ A
 *   - orthonormality  Vᵀ·V ≈ I
 *
 * Exit code 0 = pass, 1 = fail (ctest).
 * Copyright (C) 2019 - 2026 Conrad Hübler. Claude Generated. GPL-3.0.
 */

#include <Eigen/Dense>
#include <cstdio>
#include <random>

#include "src/core/energy_calculators/qm_methods/native_eigensolver.h"

int main()
{
    std::mt19937 rng(12345);
    std::normal_distribution<double> nd(0.0, 1.0);

    // Scaled relative tolerances (eigenvalue error grows ~n·eps·||A||).
    const double tol_eval  = 1e-8;
    const double tol_recon = 1e-9;
    const double tol_orth  = 1e-10;

    bool ok = true;
    std::printf("%6s  %12s  %12s  %12s  %s\n", "n", "d(eval)", "recon", "orth", "result");
    for (int n : {1, 2, 3, 5, 13, 40, 137, 300, 558}) {
        Eigen::MatrixXd M(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) M(i, j) = nd(rng);
        const Eigen::MatrixXd A = 0.5 * (M + M.transpose());

        Eigen::VectorXd ev;
        Eigen::MatrixXd vc;
        if (!curcuma::eigsolver::solveSymmetric(A, ev, vc)) {
            std::printf("n=%d: solveSymmetric FAILED\n", n);
            return 1;
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

        const double de    = (ev - es.eigenvalues()).cwiseAbs().maxCoeff();
        const double recon = (vc * ev.asDiagonal() * vc.transpose() - A).cwiseAbs().maxCoeff();
        const double orth  = (vc.transpose() * vc
                              - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff();

        const bool pass = (de < tol_eval) && (recon < tol_recon) && (orth < tol_orth);
        ok = ok && pass;
        std::printf("%6d  %12.3e  %12.3e  %12.3e  %s\n", n, de, recon, orth,
                    pass ? "PASS" : "FAIL");
    }
    std::printf("\nOverall: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
