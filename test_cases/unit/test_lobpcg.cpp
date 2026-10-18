/*
 * test_lobpcg — validate the seeded block LOBPCG kernel
 * (curcuma::eigsolver::lobpcgLowest, the iterative `-eigensolver lobpcg` path) against the
 * lowest-k eigenpairs from Eigen's SelfAdjointEigenSolver.
 *
 * Covers what the LOBPCG must get right:
 *   - eigenvalues of the k lowest states (cold random seed)
 *   - the invariant subspace ‖X·Xᵀ − Σ_{i<k} v_i v_iᵀ‖ (degeneracy/sign-robust)
 *   - warm start: seeding with the previous solution on a perturbed matrix converges fast + correct
 *
 * LOBPCG is validated here in its proper regime (k ≪ n); at k ≈ n/2 it may not converge within the
 * iteration cap (the documented dense net-loss), in which case the SCF wiring falls back to dsyevd.
 *
 * Exit code 0 = pass, 1 = fail (ctest).
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. Claude Generated. GPL-3.0.
 */

#include <Eigen/Dense>
#include <cstdio>
#include <random>

#include "src/core/energy_calculators/qm_methods/native_eigensolver.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

static std::mt19937 rng(7);

static MatrixXd randSym(int n)
{
    MatrixXd M(n, n);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) M(i, j) = nd(rng);
    return 0.5 * (M + M.transpose());
}

// Subspace distance between two orthonormal bases (columns): ‖Pa − Pb‖∞ of their projectors.
static double subspaceDist(const MatrixXd& Xa, const MatrixXd& Xb)
{
    return (Xa * Xa.transpose() - Xb * Xb.transpose()).cwiseAbs().maxCoeff();
}

// Block size = k wanted + guard vectors (the gapless boundary state converges slowly, so the
// guards absorb it and let the lowest k converge — exactly how the SCF requests nocc + buffer).
static int guardFor(int k) { return std::max(6, k / 2); }

static bool checkCold(const char* name, int n, int k)
{
    const MatrixXd A = randSym(n);
    const int kb = k + guardFor(k);
    MatrixXd X;                       // empty → cold random seed
    VectorXd ev;
    if (!curcuma::eigsolver::lobpcgLowest(A, kb, X, ev, 200, 1e-9, k)) {
        std::printf("%-14s n=%4d k=%3d  FAILED (no convergence)\n", name, n, k);
        return false;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    const double scale = std::max(1.0, A.cwiseAbs().maxCoeff());
    const double de = (ev.head(k) - es.eigenvalues().head(k)).cwiseAbs().maxCoeff() / scale;
    const double ds = subspaceDist(X.leftCols(k), es.eigenvectors().leftCols(k));
    const bool ok = (de < 1e-7) && (ds < 1e-6);
    std::printf("%-14s n=%4d k=%3d  d(eval)=%.2e d(subspace)=%.2e  %s\n",
                name, n, k, de, ds, ok ? "PASS" : "FAIL");
    return ok;
}

// Warm start: solve A, perturb to A+εΔ, re-solve SEEDED with the previous vectors; check it
// converges (few iters) to the new lowest-k subspace — the SCF recycling use case.
static bool checkWarm(const char* name, int n, int k)
{
    const MatrixXd A0 = randSym(n);
    const int kb = k + guardFor(k);
    MatrixXd X;
    VectorXd ev;
    if (!curcuma::eigsolver::lobpcgLowest(A0, kb, X, ev, 200, 1e-9, k)) {
        std::printf("%-14s n=%4d k=%3d  FAILED (cold)\n", name, n, k);
        return false;
    }
    const MatrixXd A1 = A0 + 1e-3 * randSym(n);          // small perturbation (slow SCF drift)
    if (!curcuma::eigsolver::lobpcgLowest(A1, kb, X, ev, 30, 1e-9, k)) {  // tight cap: seed must help
        std::printf("%-14s n=%4d k=%3d  FAILED (warm, >30 iters)\n", name, n, k);
        return false;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A1);
    const double scale = std::max(1.0, A1.cwiseAbs().maxCoeff());
    const double de = (ev.head(k) - es.eigenvalues().head(k)).cwiseAbs().maxCoeff() / scale;
    const double ds = subspaceDist(X.leftCols(k), es.eigenvectors().leftCols(k));
    const bool ok = (de < 1e-7) && (ds < 1e-6);
    std::printf("%-14s n=%4d k=%3d  d(eval)=%.2e d(subspace)=%.2e  %s\n",
                name, n, k, de, ds, ok ? "PASS" : "FAIL");
    return ok;
}

int main()
{
    bool ok = true;
    ok &= checkCold("cold-small",   60,   6);
    ok &= checkCold("cold-mid",    200,  16);
    ok &= checkCold("cold-large",  400,  20);
    ok &= checkWarm("warm-mid",    200,  16);
    ok &= checkWarm("warm-large",  400,  24);

    std::printf("\nOverall: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
