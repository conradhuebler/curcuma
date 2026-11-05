/*
 * <Native QM: SCF building blocks shared by the closed-shell (qm_scf.cpp) and the open-shell (qm_uscf.cpp) loop>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Moved unchanged out of qm_scf.cpp in Nov 2026 (stage 5b, UHF/UKS) so that both SCF loops use the same DIIS/ADIIS
 * subspace logic, convergence test and threaded GEMM. Internal header: included only by qm_scf.cpp and qm_uscf.cpp;
 * everything sits in an unnamed namespace (one copy per translation unit, as before the move).
 *
 * Claude Generated (Sep/Oct 2026, moved Nov 2026).
 *
 * This program is free software under GPL-3.0
 */
#pragma once

#include "src/core/global.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace {

/// C = op(A) * B with the columns of B (and C) split into one block per OpenMP thread, each
/// block one sequential BLAS call through Eigen (the release build links the sequential MKL,
/// so a plain product runs on one core). Same arithmetic per column; blocks only change
/// which call computes a column. Used for the n^3 products of the SCF loop (DIIS error,
/// X^T F X, density), which took seconds per iteration at ~2000 basis functions.
/// Claude Generated (Oct 2026).
Eigen::MatrixXd gemmThreaded(const Eigen::MatrixXd& A, bool transA, const Eigen::MatrixXd& B, int threads)
{
    const int m = static_cast<int>(transA ? A.cols() : A.rows());
    const int n = static_cast<int>(B.cols());
    Eigen::MatrixXd C(m, n);
    const int nb = std::max(1, std::min(threads, n / 64));
    if (nb <= 1) {
        if (transA)
            C.noalias() = A.transpose() * B;
        else
            C.noalias() = A * B;
        return C;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nb)
#endif
    for (int t = 0; t < nb; ++t) {
        const int c0 = static_cast<int>((long)n * t / nb), c1 = static_cast<int>((long)n * (t + 1) / nb);
        if (c1 <= c0) continue;
        if (transA)
            C.middleCols(c0, c1 - c0).noalias() = A.transpose() * B.middleCols(c0, c1 - c0);
        else
            C.middleCols(c0, c1 - c0).noalias() = A * B.middleCols(c0, c1 - c0);
    }
    return C;
}
// =================================================================================
// SCF convergence accelerator for the QM engine: Pulay DIIS and ADIIS
// (Claude Generated, Sep 2026). Self-contained here on purpose -- the shared
// DIISAccelerator (diis_accelerator.h) also drives the native xTB SCF, whose
// energies are gated at 1e-8 against tblite, and stays untouched.
//
//  DIIS  (P. Pulay, Chem. Phys. Lett. 73, 393 (1980); J. Comput. Chem. 3, 556 (1982)):
//        F = sum c_i F_i with sum c_i = 1 minimising |sum c_i e_i|^2, where
//        e_i = X^T (F_i P_i S - S P_i F_i) X is the commutator in the ORTHONORMAL
//        basis (X = S^-1/2) -- the AO-basis commutator is not a proper residual
//        when S is far from the identity, which is the Pulay-1982 recommendation.
//  ADIIS (X. Hu, W. Yang, J. Chem. Phys. 132, 054109 (2010)): far from convergence
//        the commutator is a poor guide; instead minimise the quadratic (ARH) model
//        of the energy over convex combinations of the stored densities,
//          f(c) = sum c_i Tr[(P_i - P_n) F_n] + 1/2 sum c_i c_j Tr[(P_i - P_n)(F_j - F_n)],
//        c_i >= 0, sum c_i = 1  (n = newest entry; spin-summed P, so dE = Tr(F dP)).
//  Blend (R. Garza, G. Scuseria, J. Chem. Phys. 137, 054110 (2012)): ADIIS while
//        max|e| > 1e-1, DIIS below 1e-4, linear mix of the two Fock matrices between.
// =================================================================================
// Claude Generated (Oct 2026, Q-G4): the accelerator is split into the subspace logic (which
// entries exist, which coefficients to use) and the storage of F, P, e. The host loop keeps the
// matrices in HostDiisStore; the device-resident loop keeps them on the device (QMGpuSCF) and
// only ships the m x m Gram/ADIIS matrices here. Both loops therefore solve the SAME small
// problems with the SAME code; only the O(m^2 n^2) inner products are computed elsewhere.
class QMScfAccelerator {
public:
    /// History depth max_size is clamped to 2..64 (one bit per slot in m_full_mask).
    QMScfAccelerator(int max_size, bool adiis) : m_max(std::clamp(max_size, 2, 64)), m_adiis(adiis) {}

    /// Number of storage slots the history needs.
    int capacity() const { return m_max; }

    /// Reserve the storage slot for a new (F, P, e) entry; the oldest entry is dropped when
    /// the history is full and its slot reused, otherwise the lowest free slot is taken.
    /// full_build: F came from a full J/K build (see dropIncrementalEntries()). Returns the slot.
    int push(bool full_build = true)
    {
        int slot = 0;
        if ((int)m_slots.size() >= m_max) {
            slot = m_slots.front();
            m_slots.erase(m_slots.begin());
        } else {
            while (std::find(m_slots.begin(), m_slots.end(), slot) != m_slots.end()) ++slot;
        }
        m_slots.push_back(slot);
        const std::uint64_t bit = std::uint64_t(1) << slot;
        m_full_mask = full_build ? (m_full_mask | bit) : (m_full_mask & ~bit);
        return slot;
    }
    /// Claude Generated (Oct 2026): drop every entry whose Fock matrix came from an incremental
    /// J/K build. Called when DirectFockCadence switches to full builds: those Fock matrices carry
    /// the incremental screening noise, and keeping them in the Pulay subspace costs iterations
    /// (measured: 1-6 fewer iterations without them, docs/QM_GPU_ROADMAP.md section 2).
    void dropIncrementalEntries()
    {
        std::vector<int> kept;
        for (int slot : m_slots)
            if (m_full_mask & (std::uint64_t(1) << slot)) kept.push_back(slot);
        m_slots = kept;
    }
    void setLastError(double e) { m_last_err = e; }
    /// Slots of the stored entries, oldest first.
    const std::vector<int>& slots() const { return m_slots; }
    /// Which subspace data coefficients() will read for the current entries/error.
    bool needsDiis() const { return m_slots.size() >= 2 && (!m_adiis || m_last_err <= 1e-1); }
    bool needsAdiis() const { return m_slots.size() >= 2 && m_adiis && m_last_err >= 1e-4; }

    /// Coefficients c over slots() (oldest first) of the Fock matrix to diagonalise next,
    /// F = sum c_i F_i. G = error Gram matrix <e_i, e_j> (needsDiis()), a/M = ADIIS model
    /// (needsAdiis(), see HostDiisStore::adiisData()).
    Eigen::VectorXd coefficients(const Eigen::MatrixXd& G, const Eigen::VectorXd& a, const Eigen::MatrixXd& M) const
    {
        const int m = static_cast<int>(m_slots.size());
        if (m < 2) return unitNewest(m);
        if (!m_adiis || m_last_err < 1e-4) return diis(G);
        if (m_last_err > 1e-1) return adiis(a, M);
        const double w = (std::log10(m_last_err) + 4.0) / 3.0;  // 0 at 1e-4 .. 1 at 1e-1
        return w * adiis(a, M) + (1.0 - w) * diis(G);
    }

private:
    int m_max;
    bool m_adiis;
    double m_last_err = 0.0;
    std::vector<int> m_slots;
    std::uint64_t m_full_mask = 0;  // bit s set: F in slot s came from a full J/K build (slots < 64)

    static Eigen::VectorXd unitNewest(int m)
    {
        Eigen::VectorXd c = Eigen::VectorXd::Zero(std::max(1, m));
        c(c.size() - 1) = 1.0;
        return c;
    }

    static Eigen::VectorXd diis(const Eigen::MatrixXd& G)
    {
        // Solve the bordered Pulay system; if it is singular (nearly dependent error
        // vectors), drop the oldest entries until it is not.
        const int m = static_cast<int>(G.rows());
        for (int first = 0; first + 2 <= m; ++first) {
            const int n = m - first;
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n + 1, n + 1);
            B.topLeftCorner(n, n) = G.bottomRightCorner(n, n);
            const double scale = B.topLeftCorner(n, n).diagonal().maxCoeff();
            if (scale <= 0.0) return unitNewest(m);
            B.topLeftCorner(n, n) /= scale;  // conditioning only; the solution is scale-free
            for (int i = 0; i < n; ++i) B(i, n) = B(n, i) = -1.0;
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
            rhs(n) = -1.0;
            Eigen::FullPivLU<Eigen::MatrixXd> lu(B);
            if (lu.rank() < n + 1 || lu.rcond() < 1e-14) continue;
            const Eigen::VectorXd c = lu.solve(rhs).head(n);
            if (!c.allFinite()) continue;
            Eigen::VectorXd full = Eigen::VectorXd::Zero(m);
            full.tail(n) = c;
            return full;
        }
        return unitNewest(m);
    }

    static Eigen::VectorXd adiis(const Eigen::VectorXd& a, Eigen::MatrixXd M)
    {
        const int n = static_cast<int>(a.size());
        M = 0.5 * (M + M.transpose()).eval();
        // c_i = t_i^2 / sum t^2 keeps c on the simplex; gradient descent with
        // backtracking on t (n <= subspace size, so this is cheap).
        auto fval = [&](const Eigen::VectorXd& c) { return a.dot(c) + 0.5 * c.dot(M * c); };
        auto toC = [](const Eigen::VectorXd& t) { return Eigen::VectorXd(t.cwiseAbs2() / t.squaredNorm()); };
        Eigen::VectorXd t = Eigen::VectorXd::Constant(n, 1.0);
        t(n - 1) = 2.0;  // lean towards the newest density
        Eigen::VectorXd c = toC(t);
        double f = fval(c), step = 1.0;
        for (int it = 0; it < 500; ++it) {
            const Eigen::VectorXd g = a + M * c;
            const double cg = c.dot(g);
            Eigen::VectorXd gt(n);
            const double S = t.squaredNorm();
            for (int k = 0; k < n; ++k) gt(k) = 2.0 * t(k) / S * (g(k) - cg);
            if (gt.norm() < 1e-12) break;
            bool accepted = false;
            for (int ls = 0; ls < 40; ++ls) {
                const Eigen::VectorXd tn = t - step * gt;
                const Eigen::VectorXd cn = toC(tn);
                const double fn = fval(cn);
                if (fn < f - 1e-4 * step * gt.squaredNorm()) {
                    t = tn; c = cn; f = fn; step *= 2.0; accepted = true;
                    break;
                }
                step *= 0.5;
            }
            if (!accepted) break;
        }
        return c;
    }
};

// Claude Generated (Oct 25 2026): SCF convergence test shared by the host and the device loop
// (-qm.scf_convergence auto|density). Per iteration n the loop feeds
//   E_n  = 1/2 Tr(P_n (H + F(P_n)))     energy of the density whose Fock matrix was just built,
//   g_n  = ||X^T (F_n P_n S - S P_n F_n) X||_F   orbital gradient of that density,
//   dP_n = ||P_(n+1) - P_n||_F           density change of the step just taken.
// auto: converged when |E_n - E_(n-1)| < thr_E, g_n < thr_G and dP_n < thr_dP (all three, as ORCA's
// TolE/TolErr/TolRMSP family and PySCF's conv_tol/conv_tol_grad); the loop then finishes with P_(n+1),
// one step beyond the tested density. The HF energy is stationary, its error is second order in g;
// gradients and correlation energies are first order, hence the tight thr_G for them (QMEngine
// m_scf_tight_now). density: dP_n < thr_dP only (the rule before Oct 25 2026).
// Rationale and measurements: docs/QM_SCF_CONVERGENCE_2026-10.md.
struct ScfConvergence {
    bool density_only = false;
    double thr_dP = 1e-6, thr_E = 1e-9, thr_G = 1e-5;
    bool have_prev = false;
    double e_prev = 0.0;
    double dE = std::numeric_limits<double>::infinity();
    double grad = std::numeric_limits<double>::infinity();

    /// Feed E_n and g_n of this iteration's Fock build.
    void observe(double e_n, double g_n)
    {
        dE = have_prev ? e_n - e_prev : std::numeric_limits<double>::infinity();
        e_prev = e_n;
        have_prev = true;
        grad = g_n;
    }
    /// The energy of the previous build is not comparable (precision or screening phase change).
    void resetEnergy() { have_prev = false; }
    bool converged(double delta_P) const
    {
        if (!(delta_P < thr_dP)) return false;
        if (density_only) return true;
        return std::abs(dE) < thr_E && grad < thr_G;
    }
};

/**
 * @brief Level shift of the native HF/KS SCF (-qm.scf_level_shift b > 0; V. R. Saunders, I. H. Hillier,
 * Int. J. Quantum Chem. 7, 699 (1973)). Claude Generated (Nov 2026).
 *
 * The matrix diagonalised in one iteration is F' = F + b (S - w S P S) with P the density F was built from (w = 1/2
 * for the spin-summed closed-shell P, 1 for a spin density). In the Lowdin basis S - w S P S is the projector onto the
 * virtual space of P, so the virtual orbitals of the current density move up by b and the occupied-virtual rotation of
 * a step shrinks from about g / (e_a - e_i) to g / (e_a - e_i + b): a damped, more robust SCF. At a fixed point P
 * commutes with F' and therefore with F, so the fixed point is the unshifted SCF solution. The DIIS error, the energy
 * and every convergence test use the unshifted F; the shift is removed when the orbital gradient falls below `off`
 * (> 0). When a shifted iteration meets every convergence test, the loop diagonalises the unshifted (extrapolated) F
 * of that iteration once more and returns those orbitals, orbital energies and that density, i.e. the same rule as the
 * run without shift (a further unshifted loop was tried first: plain undamped Roothaan steps can leave a converged
 * point again, ozone/def2-SVP did not reconverge).
 */
struct ScfLevelShift {
    double b = 0.0, off = 0.0;
    bool on = false;
    ScfLevelShift(double shift, double off_grad) : b(shift), off(off_grad), on(shift > 0.0) {}
    /// After the convergence measures of this Fock build: early removal below the gradient threshold.
    void observe(double grad)
    {
        if (on && off > 0.0 && grad < off) on = false;
    }
    /// F + b (S - w S P S).
    Matrix apply(const Matrix& F, const Matrix& S, const Matrix& P, double w, int threads) const
    {
        const Eigen::MatrixXd Sc(S), Pc(P);
        const Eigen::MatrixXd SPS = gemmThreaded(gemmThreaded(Sc, false, Pc, threads), false, Sc, threads);
        return F + Matrix(b * (Sc - w * SPS));
    }
};

ScfConvergence makeConvergence(bool density_only, double thr_dP, double thr_E, double thr_G)
{
    ScfConvergence c;
    c.density_only = density_only;
    c.thr_dP = thr_dP;
    c.thr_E = thr_E;
    c.thr_G = thr_G;
    return c;
}

}  // namespace
