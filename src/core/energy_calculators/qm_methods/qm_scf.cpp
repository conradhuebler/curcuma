/*
 * <Native KS-DFT -- WP3 closed-shell HF-SCF>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Closed-shell Roothaan-Hall SCF for the native DFT engine. This is the
 * "hard ERI gate" of the roadmap: it exercises the WP1 1e matrices and the
 * WP2 4-centre ERI end-to-end and, for the first time, makes the MO spectrum
 * + total energy vs ORCA a valid comparison (ORCA orbital energies are the SCF
 * Fock spectrum H + 2J - K, not the 1e Hcore that WP1 had alone).
 *
 * Convention (Szabo-Ostlund / PySCF, spin-summed density):
 *   P_munu  = 2 * sum_{i=1}^{n_occ} C_mui C_nui           (factor 2 = spin)
 *   J_munu = sum_{lam,sig} P_lamsig (mu nu | lam sig)     (Coulomb)
 *   K_munu = sum_{lam,sig} P_lamsig (mu lam | nu sig)     (exchange)
 *   F      = H + J - 1/2 K                                 (RHF Fock)
 *   E_elec = 1/2 * Tr(P (H + F))                           (total electronic)
 *   components: ET = Tr(P T), EV = Tr(P V), EJ = 1/2 Tr(P J),
 *               Ex = -1/4 Tr(P K);  E_elec = (ET + EV) + EJ + Ex.
 * (The roadmap note used the equivalent spatial-density form F = H + 2J - K;
 * both give the same Fock eigenvalues and the same total energy. The form here
 * is the standard spin-summed one and matches buildCoulomb/buildExchange with a
 * factor-2 density directly.)
 *
 * SCF driver mirrors the NDDO Fock-DIIS pattern (nddo.cpp runSCF): build Fock,
 * DIIS-extrapolate once the history is long enough, diagonalize, rebuild the
 * density, damp, check ||dP||. The generalized eigenproblem F C = S C eps is
 * reduced to the standard one via the one-time Lowdin orthonormalizer
 * X = S^{-1/2} (no generalized eigensolver exists in curcuma::eigsolver).
 *
 * Literature:
 *   - P. Pulay, Chem. Phys. Lett. 73, 393 (1980); J. Comput. Chem. 3, 556 (1982).
 *   - T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure
 *     Theory (Wiley, 2000), ch. 10 (SCF) + ch. 9.9 (ERI).
 *   - A. Szabo, N. S. Ostlund, Modern Quantum Chemistry (Dover, 1996), ch. 3.
 *   - xcDFT RKS.f90 (TCCM winter school 2019: DFT) -- structural SCF template.
 *
 * Claude Generated: WP3 HF-SCF
 *
 * This program is free software under GPL-3.0
 */

#include "qm_citations.h"
#include "qm_engine.h"
#include "qm_ri.h"  // RIJK::coulomb (WP5 J-only build)

#include <cstdlib>
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "src/core/charge_extrapolation.h"

#include "native_eigensolver.h"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>

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
// Symmetric Loewdin orthonormalizer X = S^{-1/2} via S = U diag(w) U^T,
// X = U diag(1/sqrt(w)) U^T. Returns true on success; leaves X empty on a
// singular / indefinite S (caller reports the error).
// Claude Generated: Lowdin S^{-1/2} for the generalized SCF eigenproblem.
bool lowdinOrthonormalizer(const Eigen::MatrixXd& S, Eigen::MatrixXd& X)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    if (es.info() != Eigen::Success) return false;
    const Eigen::VectorXd w = es.eigenvalues();
    const Eigen::MatrixXd U = es.eigenvectors();
    // Reject a singular / indefinite overlap (linearly dependent basis).
    for (int i = 0; i < w.size(); ++i)
        if (w(i) < 1.0e-10) return false;
    Eigen::VectorXd inv_sqrt(w.size());
    for (int i = 0; i < w.size(); ++i) inv_sqrt(i) = 1.0 / std::sqrt(w(i));
    X = U * inv_sqrt.asDiagonal() * U.transpose();
    return true;
}
}  // namespace


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
namespace {

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

// Integral-direct SCF: full J/K rebuild after this many incremental builds (host and device loop).
constexpr int kDirectRebuild = 8;

// Claude Generated (Oct 2026): Fock-build cadence of the integral-direct SCF, shared by the host
// and the device loop. Incremental builds on dP = P - P_ref are screened with the same absolute
// threshold as a full build (Q_AB Q_CD max|dP| < eri_screening), so each one drops terms of up to
// eri_screening per quartet whatever the size of dP. That error does not shrink with dP: the
// accumulated Fock matrix carries a path-dependent noise of ~1e-10 Eh (naphthalene/def2-SVP at
// eri_screening 1e-12), which shows up as a floor of ||dP||_F ~ 1e-8 (benzene) to 3e-7 (anthracene,
// pyrene); below it the SCF cannot converge (docs/QM_GPU_ROADMAP.md section 2, "Incremental-build floor").
// A full build is a fixed function of P (its screen only sees P), so full builds converge like the
// stored tensor. Rule: once ||dP|| (below kIncrementalStallMaxDP) has not fallen below
// kIncrementalProgress times its smallest earlier value for kIncrementalStallCount iterations in a
// row, the rest of this SCF uses full builds. A run that converges before it stalls (every default
// run measured) is unchanged; a false trigger costs time, never accuracy. The constants were chosen
// on recorded ||dP|| traces at scf_threshold 1e-10: no trigger in the stored-tensor and full-build
// traces of benzene, naphthalene, anthracene; the pyrene full-build trace (slow tail near 2e-10)
// would trigger 3 iterations before convergence; on incremental traces it fires 1-4 iterations
// after the floor is reached. At the switch the DIIS history drops its entries from incremental
// builds (QMScfAccelerator::dropIncrementalEntries()).
constexpr int kIncrementalStallCount = 2;
constexpr double kIncrementalProgress = 0.5;
constexpr double kIncrementalStallMaxDP = 1.0e-4;

struct DirectFockCadence {
    bool direct = false;
    bool have_ref = false;          // a full build exists (J/K of P_ref accumulated)
    int since_full = 0;             // incremental builds since the last full one
    bool incremental_off = false;   // stall seen: full builds for the rest of this SCF
    double best_dP = std::numeric_limits<double>::infinity();
    int stall = 0;
    // Claude Generated (Oct 2026): early phase with a looser screening threshold (-qm.scf_early_screening).
    // It ends once ||dP|| < -qm.scf_early_switch; the next build is then a full one with eri_screening, so
    // nothing accumulated under the loose threshold survives into the converged Fock matrix.
    bool early = false;
    // Claude Generated (Oct 2026, Q-G5): FP32 phase of the device direct builder (-qm.gpu_eri_precision).
    // It ends once ||dP|| < -qm.gpu_fp32_switch or after kFP32NoProgress iterations without a new
    // smallest ||dP|| (FP32 noise floor); the next build is then a full FP64 one, so no FP32 J/K
    // contribution survives, and convergence is only accepted with fp32 == false.
    bool fp32 = false;
    int fp32_builds = 0;
    double fp32_best = std::numeric_limits<double>::infinity();
    int fp32_noprog = 0;
    static constexpr int kFP32NoProgress = 3;

    explicit DirectFockCadence(bool d, bool early_phase = false, bool fp32_phase = false)
        : direct(d), early(d && early_phase), fp32(d && fp32_phase) {}
    /// Feed ||P_new - P||_F; true when the FP32 phase just ended.
    bool endFP32(double delta_P, double switch_dP)
    {
        if (!fp32) return false;
        if (delta_P < fp32_best) {
            fp32_best = delta_P;
            fp32_noprog = 0;
        } else {
            ++fp32_noprog;
        }
        if (delta_P >= switch_dP && fp32_noprog < kFP32NoProgress) return false;
        fp32 = false;
        have_ref = false;  // next build: full, FP64
        return true;
    }
    /// Feed ||P_new - P||_F; true when the early (loose-screening) phase just ended.
    bool endEarly(double delta_P, double switch_dP)
    {
        if (!early || delta_P >= switch_dP) return false;
        early = false;
        have_ref = false;  // next build: full, with the final threshold
        return true;
    }
    bool nextFull() const { return !direct || !have_ref || incremental_off || since_full >= kDirectRebuild; }
    void built(bool full)
    {
        have_ref = true;
        since_full = full ? 0 : since_full + 1;
    }
    /// Feed ||P_new - P||_F of this iteration; true when incremental builds were just switched off.
    bool observe(double delta_P)
    {
        if (!direct || incremental_off) return false;
        const bool progress = delta_P < kIncrementalProgress * best_dP;
        best_dP = std::min(best_dP, delta_P);
        if (progress) {
            stall = 0;
            return false;
        }
        if (delta_P < kIncrementalStallMaxDP && ++stall >= kIncrementalStallCount) {
            incremental_off = true;
            return true;
        }
        return false;
    }
};

// One line at verbosity >= 2 when the incremental scheme is switched off (diagnostic only).
void reportIncrementalStall(int iter, double best_dP, double screening)
{
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "QM: integral-direct SCF: incremental J/K builds stalled near ||dP|| = {:.2e} (screening noise of "
            "eri_screening {:.1e}); full J/K builds from iteration {} on",
            best_dP, screening, iter + 2));
}

// One line at verbosity >= 2 when the FP32 phase ends (diagnostic only). Claude Generated (Oct 2026).
void reportFP32End(int iter, double delta_P, int builds)
{
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "QM: integral-direct SCF: ||dP|| = {:.2e}, FP32 phase ends after {} FP32 J/K builds; full FP64 J/K build "
            "from iteration {} on",
            delta_P, builds, iter + 2));
}

// One line at verbosity >= 2 when the early (loose-screening) phase ends (diagnostic only).
void reportEarlyEnd(int iter, double delta_P)
{
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "QM: integral-direct SCF: ||dP|| = {:.2e}, early screening phase ends; full J/K build with eri_screening "
            "from iteration {} on",
            delta_P, iter + 2));
}

// Host storage of the DIIS/ADIIS history (the device loop keeps the same data in QMGpuSCF).
struct HostDiisStore {
    std::vector<Matrix> F, P, e;  // indexed by slot

    explicit HostDiisStore(int slots) : F(slots), P(slots), e(slots) {}
    ~HostDiisStore();  // out of line: -Winline in the grown runSCF (Oct 2026)

    /// Store F (built from P) in `slot`, form its orthonormal-basis error and return max|e|.
    double push(int slot, const Matrix& Fi, const Matrix& Pi, const Matrix& S, const Matrix& X, int threads)
    {
        e[slot] = orbitalGradient(Fi, Pi, S, X, threads);
        F[slot] = Fi;
        P[slot] = Pi;
        return e[slot].cwiseAbs().maxCoeff();
    }
    /// Orbital gradient (= DIIS error) of the closed-shell SCF in the Lowdin-orthonormal basis,
    /// e = X^T (F P S - S P F) X; zero exactly at a stationary density (Pulay 1982).
    static Matrix orbitalGradient(const Matrix& Fi, const Matrix& Pi, const Matrix& S, const Matrix& X, int threads)
    {
        const Eigen::MatrixXd Fc(Fi), Pc(Pi), Sc(S), Xc(X);
        const Eigen::MatrixXd FPS = gemmThreaded(gemmThreaded(Fc, false, Pc, threads), false, Sc, threads);
        const Eigen::MatrixXd comm = FPS - FPS.transpose();
        return Matrix(gemmThreaded(gemmThreaded(Xc, true, comm, threads), false, Xc, threads));
    }
    Eigen::MatrixXd gram(const std::vector<int>& sl) const
    {
        const int m = static_cast<int>(sl.size());
        Eigen::MatrixXd G(m, m);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j <= i; ++j) G(i, j) = G(j, i) = e[sl[i]].cwiseProduct(e[sl[j]]).sum();
        return G;
    }
    /// ADIIS (ARH) model, n = newest: a_i = Tr[(P_i - P_n) F_n], M_ij = Tr[(P_i - P_n)(F_j - F_n)].
    void adiisData(const std::vector<int>& sl, Eigen::VectorXd& a, Eigen::MatrixXd& M) const
    {
        const int m = static_cast<int>(sl.size());
        const Matrix& Pn = P[sl.back()];
        const Matrix& Fn = F[sl.back()];
        a.resize(m);
        M.resize(m, m);
        for (int i = 0; i < m; ++i) {
            const Matrix dP = P[sl[i]] - Pn;
            a(i) = dP.cwiseProduct(Fn).sum();
            for (int j = 0; j < m; ++j) M(i, j) = dP.cwiseProduct(F[sl[j]] - Fn).sum();
        }
    }
    Matrix combine(const std::vector<int>& sl, const Eigen::VectorXd& c) const
    {
        Matrix out = Matrix::Zero(F[sl[0]].rows(), F[sl[0]].cols());
        for (int i = 0; i < c.size(); ++i) out += c(i) * F[sl[i]];
        return out;
    }
};

HostDiisStore::~HostDiisStore() = default;

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

ScfConvergence makeConvergence(bool density_only, double thr_dP, double thr_E, double thr_G)
{
    ScfConvergence c;
    c.density_only = density_only;
    c.thr_dP = thr_dP;
    c.thr_E = thr_E;
    c.thr_G = thr_G;
    return c;
}

// Claude Generated (Oct 27 2026): verified stop of the RHF loop (-qm.scf_verified_stop, host and device loop).
// The rule before: converged at iteration n when |dE_n|, g_n (orbital gradient of P_n) and dP_n = ||P_(n+1) - P_n||
// are small; the loop then returned P_(n+1) after one more full J/K build. Here the test is applied to the density
// that is returned: at the top of iteration n, right after a FULL (exact) build of P_n, |dE_n| < thr_E, g_n < thr_G
// and dP_(n-1) = ||P_n - P_(n-1)|| < thr_dP mean P_n, F(P_n) and the orbitals C_(n-1) that P_n was built from are the
// result (the same structure as before: the orbitals come from the previous, DIIS-extrapolated Fock matrix). This
// saves the final build, and the returned density meets the gradient threshold itself. Integral-direct builds are
// mostly incremental; when g_(n-1)^2 / g_(n-2) (geometric decay) and dP_(n-1) predict convergence, build n is made
// a full one. A wrong prediction costs one full instead of one incremental build. Needs undamped updates
// (P_n = 2 C_(n-1),occ C_(n-1),occ^T), i.e. diis/adiis, and scf_convergence auto.
// Energy-only SCFs only: the old rule returned P_(n+1), about one DIIS step (factor ~0.25 in g) better than the
// tested P_n; the HF energy is second order in g (no visible change), but gradients and correlation energies are
// first order (methanol HF gradient 5e-8 instead of 5e-9 Eh/A from the tight reference), so tight SCFs
// (m_scf_tight_now) keep the old rule.
struct ScfVerifiedStop {
    bool enabled = false;
    bool force_full = false;           // the next build is a full one
    bool have_prev = false;            // C_prev / eps_prev / dP_prev belong to the density of the next build
    double dP_prev = std::numeric_limits<double>::infinity();
    double g_older = std::numeric_limits<double>::infinity();  // orbital gradient one iteration before conv.grad
    Matrix C_prev;
    Vector eps_prev;

    explicit ScfVerifiedStop(bool on) : enabled(on) {}
    ~ScfVerifiedStop();  // out of line: -Winline (as HostDiisStore)
    /// Top of an iteration, after conv.observe(): true when this build's density is the converged result.
    bool accept(const ScfConvergence& conv, bool exact_full_build) const
    {
        return enabled && have_prev && exact_full_build && conv.converged(dP_prev);
    }
    /// End of an iteration that did not converge: remember dP and decide whether the next build must be full.
    void afterStep(const ScfConvergence& conv, double delta_P, bool may_force)
    {
        if (!enabled) return;
        have_prev = true;
        dP_prev = delta_P;
        const double g = conv.grad;
        const double g_pred = (std::isfinite(g_older) && g < g_older && g_older > 0.0) ? g * g / g_older : g;
        force_full = may_force && delta_P < conv.thr_dP && g_pred < conv.thr_G;
        g_older = g;
    }
};

ScfVerifiedStop::~ScfVerifiedStop() = default;

}  // namespace

// =================================================================================
// SCF driver
// =================================================================================

bool QMEngine::runSCF()
{
    // Claude Generated: WP3 closed-shell RHF SCF (DIIS / plain damping).
    const int n = m_nbf;
    if (n <= 0) {
        CurcumaLogger::error("QM SCF: no basis functions (InitialiseMolecule failed?)");
        return false;
    }
    // Closed-shell only: even electron count.
    if (m_num_electrons < 0 || (m_num_electrons % 2) != 0) {
        CurcumaLogger::error(fmt::format(
            "QM SCF: closed-shell HF needs an even electron count, got {}",
            m_num_electrons));
        return false;
    }
    const int n_occ = m_num_electrons / 2;
    if (n_occ > n) {
        CurcumaLogger::error(fmt::format(
            "QM SCF: {} electrons need {} occupied orbitals but basis has only {}",
            m_num_electrons, n_occ, n));
        return false;
    }

    // One-time setup: Lowdin orthonormalizer + active-basis ERI.
    buildOrthonormalizer();
    if (m_X.size() == 0) {
        CurcumaLogger::error("QM SCF: S^{-1/2} construction failed (linear basis?)");
        return false;
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        Eigen::MatrixXd Scol(m_S), Hcol(m_H);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esS(Scol);
        const Eigen::VectorXd sw = esS.eigenvalues();
        std::string svals = "S eigvals:";
        for (int i = 0; i < sw.size(); ++i) svals += fmt::format(" {:.4e}", sw(i));
        CurcumaLogger::info(svals);
        // Cross-check: direct generalized eigensolver for the bare H (iter 0 Fock).
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Hcol, Scol);
        const Eigen::VectorXd hv = ges.eigenvalues();
        std::string hvals = "H gen-eigvals:";
        for (int i = 0; i < std::min<int>(5, hv.size()); ++i)
            hvals += fmt::format(" {:.4f}", hv(i));
        CurcumaLogger::info(hvals);
    }
    // Two-electron part: stored tensor (built once here) or integral-direct. WP10: RI-JK needs
    // neither (no 4-index tensor, full builds only); RI-J keeps the exact path for the exchange.
    const bool ri_jk = riActive() && m_ri_mode == "jk";
    if (riActive() && !ensureRI()) return false;
    // WP5: Kohn-Sham functional -> grid of this geometry (built once, cached until the geometry changes)
    if (isKohnSham()) {
        if (!ensureXCGrid()) return false;
        qmcite::kohnSham(m_xc, m_grid_options.radial == "treutler");
    }
    const bool direct = !ri_jk && useDirectSCF();
    if (direct) qmcite::directScf();
    if (m_scf_mode == "adiis") CitationRegistry::cite("adiis");
    if (!ri_jk && !direct && !ensureGpuJK()) {  // KS: J only, from the same tensor
        const qmint::ERITensor& eri = eriActive();  // builds (and caches) on first call
        if (eri.n() != n) {
            CurcumaLogger::error(fmt::format(
                "QM SCF: active ERI dimension {} != basis dimension {}", eri.n(), n));
            return false;
        }
    }
    m_direct_quartets = 0;

    const bool use_diis = (m_scf_mode != "plain");   // diis | adiis (other names are rejected by the constructor)
    // Fraction of the OLD density kept per step: plain mode 1 - scf_damping (default 0.5, the value before Oct 27 2026)
    const double damping = use_diis ? 0.0 : 1.0 - m_scf_damping;

    // Initial guess (m_scf_guess): SAD (default) or the bare core Hamiltonian
    // (zero density -> first Fock = H).
    m_density = buildInitialGuess();
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::param("SCF guess",
                             fmt::format("{} (Tr(PS) = {:.6f})", m_last_warm_start ? "warm start" : m_scf_guess,
                                         m_density.cwiseProduct(m_S).sum()));

    // Q-G4: device-resident loop (same iterations on device matrices); -1 = not usable or a
    // device error, then the host loop below runs from the same start density.
    // KS (WP5): the device loop has no V_xc, the host loop runs (J/K may still come from the device).
    if (!riActive() && !isKohnSham() && gpuRequested() && m_gpu_scf != "off" && m_gpu_scf != "false" && !m_gpu_scf_failed) {
        const Matrix guess = m_density;
        const int r = runSCFDevice(n_occ, direct, use_diis, damping);
        if (r >= 0) return r == 1;
        m_density = guess;
        m_direct_quartets = 0;
    }

    QMScfAccelerator diis(m_diis_subspace, m_scf_mode == "adiis");
    HostDiisStore store(diis.capacity());
    ScfConvergence conv = makeConvergence(m_scf_convergence == "density", m_scf_threshold, m_scf_energy_threshold,
                                          m_scf_tight_now ? m_scf_gradient_threshold_tight : m_scf_gradient_threshold);
    const auto t_scf0 = std::chrono::steady_clock::now();
    double t_fock_ms = 0.0, t_eig_ms = 0.0;

    m_scf_converged = false;
    m_scf_iterations = 0;

    // Incremental Fock build (integral-direct only): J and K are linear in P, so
    // J(P_n) = J(P_ref) + J(P_n - P_ref). The density change shrinks as the SCF
    // converges and the density-weighted screening then skips most quartets.
    // A full build every kDirectRebuild iterations bounds the accumulated
    // screening error; the final energy always comes from a full build.
    // When ||dP|| stalls at the incremental screening noise, full builds only
    // (DirectFockCadence).
    Matrix J_acc, K_acc, P_ref;
    DirectFockCadence cadence(direct, useEarlyScreening(), useGpuFP32Phase(direct));
    // Verified stop (Claude Generated, Oct 27 2026; -qm.scf_verified_stop): see ScfVerifiedStop.
    ScfVerifiedStop vstop(m_scf_verified_stop && use_diis && !conv.density_only && !m_scf_tight_now && m_ri_mp_stage == 0);
    auto fockFromDensity = [&](const Matrix& P) -> Matrix {
        const bool full = cadence.nextFull() || vstop.force_full;
        vstop.force_full = false;
        m_direct_thr_now = cadence.early ? m_scf_early_screening : -1.0;
        m_direct_fp32_now = cadence.fp32;
        if (cadence.fp32) ++cadence.fp32_builds;
        if (full) {
            buildJKForSCF(P, J_acc, K_acc);
        } else {
            Matrix dJ, dK;
            buildJKForSCF(P - P_ref, dJ, dK);
            J_acc += dJ;
            K_acc += dK;
        }
        cadence.built(full);
        P_ref = P;
        // HF: H + J - K/2; KS: + V_xc[P] - (a_x/2) K (V_xc is not linear in P: always from the full density)
        return composeFock(P, J_acc, K_acc);
    };

    for (int iter = 0; iter < m_scf_max_iter; ++iter) {
        m_scf_iterations = iter + 1;

        // Fock from the current density (F = H + J - 1/2 K, spin-summed convention).
        const auto t_f0 = std::chrono::steady_clock::now();
        Matrix fock = fockFromDensity(m_density);
        t_fock_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_f0).count();

        if (CurcumaLogger::get_verbosity() >= 3 && iter < 2) {
            // Inspect J/K in the MO basis of the *previous* orbitals (m_prev_C).
            Matrix Jm, Km;
            buildJK(m_density, Jm, Km);
            CurcumaLogger::info(fmt::format(
                "DBG Fock diag[0..3]={:.4} {:.4} {:.4} {:.4}  J_00(AO)={:.4} K_00(AO)={:.4} Tr(PJ)={:.6} Tr(PK)={:.6}",
                fock(0, 0), n > 1 ? fock(1, 1) : 0, n > 2 ? fock(2, 2) : 0, n > 3 ? fock(3, 3) : 0,
                Jm(0, 0), Km(0, 0),
                m_density.cwiseProduct(Jm).sum(), m_density.cwiseProduct(Km).sum()));
        }

        // DIIS: push the commutator error e = F P S - S P F every iter once the
        // history can support it; use the extrapolated Fock for the diagonalization.
        Matrix fock_use = fock;
        double grad_norm = 0.0;  // ||e||_F of this Fock build (convergence test)
        if (use_diis && iter >= m_diis_start) {
            // FP32 RI builds (mixed precision) and FP32 direct builds (Q-G5) count as non-full: dropped
            // when the SCF switches to FP64
            const int slot = diis.push(cadence.since_full == 0 && !riLowPrecision() && !cadence.fp32);
            diis.setLastError(store.push(slot, fock, m_density, m_S, m_X, m_threads));
            grad_norm = store.e[slot].norm();
            Eigen::MatrixXd G, M;
            Eigen::VectorXd a;
            if (diis.needsDiis()) G = store.gram(diis.slots());
            if (diis.needsAdiis()) store.adiisData(diis.slots(), a, M);
            fock_use = store.combine(diis.slots(), diis.coefficients(G, a, M));
        } else if (!conv.density_only) {
            grad_norm = HostDiisStore::orbitalGradient(fock, m_density, m_S, m_X, m_threads).norm();
        }
        // energy of the density this Fock matrix was built from (variational, unlike Tr(P_new (H + F)))
        conv.observe(scfEnergy(m_density, fock), grad_norm);

        // Verified stop: this full build of P_n meets every condition -> P_n, F(P_n) and the orbitals C_(n-1)
        // that P_n was built from are the result; no further diagonalisation and no final build.
        if (vstop.accept(conv, cadence.since_full == 0 && !cadence.early && !cadence.fp32 && !riLowPrecision())) {
            m_scf_last_dE = conv.dE;
            m_scf_last_grad = conv.grad;
            m_direct_thr_now = -1.0;
            m_direct_fp32_now = false;
            finishConvergedSCF(vstop.C_prev, vstop.eps_prev, n_occ, iter + 1, vstop.dP_prev, direct,
                               std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scf0).count(),
                               t_fock_ms, t_eig_ms, std::string(), &J_acc, &K_acc);
            return true;
        }

        // Solve F C = S C eps (Lowdin-reduced standard problem).
        Matrix C;
        Vector eps;
        const auto t_e0 = std::chrono::steady_clock::now();
        if (!solveFock(fock_use, C, eps)) {
            CurcumaLogger::error(fmt::format(
                "QM SCF: diagonalization failed at iteration {}", iter + 1));
            return false;
        }
        t_eig_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_e0).count();

        // New density from the occupied MOs.
        Matrix density_new = buildDensity(C, n_occ);
        if (riActive() && m_ri_mode == "jk") {
            // RI-JK: keep the orbital factor of this density for the exchange build (QMEngine::riOccFactor)
            m_ri_occ_factor = std::sqrt(2.0) * Eigen::MatrixXd(C.leftCols(n_occ));
            m_ri_occ_density = density_new;
        }

        // Convergence: Frobenius norm of the density change.
        const double delta_P = (density_new - m_density).norm();

        if (CurcumaLogger::get_verbosity() >= 3 && iter < 3) {
            const double trPS = (m_density.cwiseProduct(m_S).sum());
            const double trPSnew = (density_new.cwiseProduct(m_S).sum());
            const double cSc = (C.col(0).transpose() * m_S * C.col(0)).value();
            CurcumaLogger::info(fmt::format(
                "DBG iter {} Tr(PS)={} Tr(PS_new)={} c0Sc={} eps0={:.6} eps1={:.6} ||P||={:.3}",
                iter + 1, trPS, trPSnew, cSc, eps.size() > 0 ? eps(0) : 0.0,
                eps.size() > 1 ? eps(1) : 0.0, density_new.norm()));
        }

        if (CurcumaLogger::get_verbosity() >= 1)  // one line per iteration from level 1 on (was level 3)
            CurcumaLogger::result(fmt::format("SCF iter {:3d}: E_elec = {:.10f} Eh, dE = {:.2e}, |g| = {:.2e}, dP = {:.3e}",
                                              iter + 1, conv.e_prev, conv.have_prev && iter > 0 ? conv.dE : 0.0,
                                              conv.grad, delta_P));

        // Never converged on an FP32 RI build (mixed precision): the FP32 noise floor is not the FP64 fixed point
        if (conv.converged(delta_P) && !cadence.early && !cadence.fp32 && !riLowPrecision()) {
            m_scf_last_dE = conv.dE;
            m_scf_last_grad = conv.grad;
            riRequestExactBuild();  // mixed-precision RI: the final J/K build is exact FP64
            m_density = density_new;
            m_direct_thr_now = -1.0;
            m_direct_fp32_now = false;
            finishConvergedSCF(C, eps, n_occ, iter + 1, delta_P, direct,
                               std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scf0).count(),
                               t_fock_ms, t_eig_ms, std::string());
            return true;
        }
        if (cadence.endEarly(delta_P, m_scf_early_switch)) {
            reportEarlyEnd(iter, delta_P);
            conv.resetEnergy();  // next energy from a build with another screening threshold
        }
        if (cadence.endFP32(delta_P, gpuFP32Switch())) {
            reportFP32End(iter, delta_P, cadence.fp32_builds);
            diis.dropIncrementalEntries();  // Fock matrices of FP32 builds leave the subspace
            conv.resetEnergy();
        }
        if (riEndLowPrecision(delta_P)) {  // FP32 -> exact FP64 RI builds (host loop only: RI never runs the device loop)
            diis.dropIncrementalEntries();
            conv.resetEnergy();
        }
        if (cadence.observe(delta_P)) {
            reportIncrementalStall(iter, cadence.best_dP, m_eri_screening);
            diis.dropIncrementalEntries();  // keep only Fock matrices from full builds
        }
        vstop.afterStep(conv, delta_P, direct && !cadence.early && !cadence.fp32 && !riLowPrecision());
        vstop.C_prev = std::move(C);
        vstop.eps_prev = std::move(eps);

        // Density update: damped (plain) or direct (diis -- the extrapolation
        // already smoothed the Fock). Reuses NDDO's update form.
        m_density = damping * m_density + (1.0 - damping) * density_new;
    }

    // Ran out of iterations: keep the last density/Fock for diagnostics.
    m_direct_thr_now = -1.0;
    m_direct_fp32_now = false;
    m_fock = buildFock(m_density);
    m_e_elec = scfEnergy(m_density, m_fock);
    return false;
}

// Claude Generated (Oct 2026, Q-G4): common tail of the host and the device loop. The final
// Fock matrix and the energy components always come from ONE full J/K build of the converged
// density through buildJK() (never incremental), identical for both loops.
void QMEngine::finishConvergedSCF(const Matrix& C, const Vector& eps, int n_occ, int iterations, double delta_P,
                                  bool direct, double loop_ms, double fock_ms, double eig_ms,
                                  const std::string& device_stats, const Matrix* J_in, const Matrix* K_in)
{
    Matrix J, K;
    if (J_in && K_in && J_in->rows() == m_nbf && K_in->rows() == m_nbf) {
        J = *J_in;  // verified stop: the loop's last (full, exact) build of m_density
        K = *K_in;
    } else {
        buildJKForSCF(m_density, J, K);
    }
    m_fock = composeFock(m_density, J, K);
    m_scf_converged = true;
    m_scf_tight_done = m_scf_tight_now;

    // Final energy + components from the converged density and Fock.
    m_et = m_density.cwiseProduct(m_T).sum();
    m_ev = m_density.cwiseProduct(m_V).sum();
    m_ej = 0.5 * m_density.cwiseProduct(J).sum();
    m_ex = isKohnSham() ? -0.25 * m_xc.exact_exchange * m_density.cwiseProduct(K).sum()
                        : -0.25 * m_density.cwiseProduct(K).sum();
    m_e_elec = scfEnergy(m_density, m_fock);
    if (isKohnSham()) {
        const qmxc::XCResult& xc = xcFor(m_density);
        m_exc = xc.exc;
        m_exc_x = xc.ex;
        m_exc_c = xc.ec;
        m_exc_nelec = xc.nelec;
    }

    // Store the MO spectrum for wrapper compatibility (QMDriver slots).
    m_mo = C;
    // Keep the occupied orbitals for the next geometry's warm start.
    m_warm_C = C.leftCols(n_occ);
    m_warm_atoms = m_atoms;
    m_energies = eps;
    // -qm.scf_extrapolation (Oct 27 2026): converged density per geometry, newest first; a second SCF of the same
    // geometry (tightening) replaces its entry. Trimmed to what the predictor uses (as the xtb history).
    if (m_scf_extrapolation != "none") {
        if (!m_dens_hist.empty() && m_dens_hist.front().rows() != m_nbf) m_dens_hist.clear();
        if (!m_dens_hist.empty() && m_dens_hist_serial == m_geom_serial)
            m_dens_hist.front() = m_density;
        else
            m_dens_hist.push_front(m_density);
        m_dens_hist_serial = m_geom_serial;
        const size_t maxlen = m_scf_extrapolation == "aspc" ? static_cast<size_t>(std::max(m_scf_extrap_order + 2, 2))
                                                            : static_cast<size_t>(std::min(2 * m_scf_extrap_order + 2, 10));
        while (m_dens_hist.size() > maxlen) m_dens_hist.pop_back();
    }

    // Diagnostic (Q-G4): -opt and batch -sp silence the method, so the loop choice and warm start
    // are not visible there; CURCUMA_QM_SCF_TRACE=1 prints one ASCII line per converged SCF.
    if (const char* tr = std::getenv("CURCUMA_QM_SCF_TRACE"); tr && tr[0] == '1')
        fmt::print(stderr, "QM SCF trace: loop={} n={} iterations={} start={} E_elec={:.10f}\n",
                   device_stats.empty() ? "host" : "device", m_nbf, iterations,
                   m_last_warm_start ? "warm" : m_scf_guess, m_e_elec);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format(
            "QM {} SCF converged in {} iterations (dP = {:.3e}{}, start: {}{})",
            isKohnSham() ? "KS-" + m_xc.label().substr(0, m_xc.label().find(' ')) : std::string("HF"), iterations, delta_P,
            m_scf_convergence == "density"
                ? std::string()
                : fmt::format(", dE = {:.1e}, |g| = {:.1e} < {:.0e}", m_scf_last_dE, m_scf_last_grad,
                              m_scf_tight_now ? m_scf_gradient_threshold_tight : m_scf_gradient_threshold),
            m_scf_continue ? "tightened from the converged density" : (m_last_warm_start ? "previous geometry" : m_scf_guess),
            device_stats.empty() ? std::string() : ", device-resident loop"));
        CurcumaLogger::info(fmt::format("QM: SCF loop {:.1f} ms, of which Fock (J/K) builds {:.1f} ms, eigensolve{} {:.1f} ms{}",
            loop_ms, fock_ms, device_stats.empty() ? "" : " + density", eig_ms,
            direct ? fmt::format(" (integral-direct{}, {} shell quartets computed)",
                                 m_gpu_djk ? " on " + m_gpu_djk->name() : std::string(), m_direct_quartets)
                   : std::string()));
        if (!device_stats.empty()) CurcumaLogger::info("QM: device SCF: " + device_stats);
        if (m_gpu_djk) CurcumaLogger::info("QM: device direct J/K: " + m_gpu_djk->stats());
        if (isKohnSham())
            CurcumaLogger::info(fmt::format("QM: V_xc builds {} ({:.1f} ms; grid {} points, {} above rho {:.0e} in the "
                                            "last build; {} threads)", m_xc_builds, m_xc_ms, m_xc_grid.size(),
                                            m_xc_last.points_used, m_xc_rho_threshold, m_threads));
    }
}

// Claude Generated (Oct 2026, Q-G4): device J/K builder for the resident loop -- the same
// object buildJK() would use for this SCF (stored device tensor, or device integral-direct).
qmint::QMGpuJK* QMEngine::gpuScfJK(bool direct) const
{
    if (direct) return (gpuDirectAllowed() && ensureGpuDirect()) ? m_gpu_djk.get() : nullptr;
    return ensureGpuJK() ? m_gpu_jk.get() : nullptr;
}

// Claude Generated (Oct 2026, Q-G4 of docs/QM_GPU_ROADMAP.md): device-resident SCF loop.
// Step for step the host loop of runSCF(): same Fock cadence (full / incremental on dP with a
// full rebuild every kDirectRebuild builds), same DIIS/ADIIS subspace logic (QMScfAccelerator,
// only the inner products come from the device), same convergence test ||P_new - P||_F <
// scf_threshold, same damping, same final full J/K build (finishConvergedSCF). P, F, the DIIS
// history and C stay on the device; per iteration the host gets max|e|, the m x m subspace
// matrices, ||dP|| and E, and sends m coefficients. Any device error: counted fallback
// ("qm_scf"), return -1, and the caller restarts the host loop from the same guess.
int QMEngine::runSCFDevice(int n_occ, bool direct, bool use_diis, double damping)
{
    const int n = m_nbf;
    const bool forced = (m_gpu_scf == "on" || m_gpu_scf == "true");
    if (!forced && m_gpu_scf != "auto") {
        static bool warned = false;
        if (!warned && CurcumaLogger::get_verbosity() >= 1) {
            warned = true;
            CurcumaLogger::warn(fmt::format("QM: unknown -qm.gpu_scf '{}' (use auto|on|off); using auto", m_gpu_scf));
        }
    }
    if (!forced && n < m_gpu_scf_min_nbf) return -1;  // auto: host loop is faster for small bases
    qmint::QMGpuJK* jk = gpuScfJK(direct);
    if (!jk) {
        if (forced) curcuma::reportGpuFallback("qm_scf", "J/K not on the device; host SCF loop");
        return -1;
    }
    auto fail = [&](const std::string& what) {
        m_gpu_scf_failed = true;
        m_gpu_scf_state.reset();
        curcuma::reportGpuFallback("qm_scf", "device SCF loop: " + what + " failed; host SCF loop");
        return -1;
    };
    const auto t_setup0 = std::chrono::steady_clock::now();
    if (!m_gpu_scf_state || m_gpu_scf_state->n() != n) {
        m_gpu_scf_state.reset();
        const std::string backend = m_gpu;
        std::string why;
        m_gpu_scf_state = gpu_plugin::createQmScf(backend, n, m_gpu_device, why);
        if (!m_gpu_scf_state) return fail("state creation (" + why + ")");
    }
    qmint::QMGpuSCF& dev = *m_gpu_scf_state;
    QMScfAccelerator diis(m_diis_subspace, m_scf_mode == "adiis");
    if (!dev.setup(m_H.data(), m_S.data(), m_X.data(), n_occ, diis.capacity())) return fail("setup");
    if (!dev.setDensity(m_density.data())) return fail("density upload");
    const double setup_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_setup0).count();

    const auto t_scf0 = std::chrono::steady_clock::now();
    double t_fock_ms = 0.0, t_eig_ms = 0.0;
    m_scf_converged = false;
    m_scf_iterations = 0;
    DirectFockCadence cadence(direct, useEarlyScreening(), useGpuFP32Phase(direct));
    ScfConvergence conv = makeConvergence(m_scf_convergence == "density", m_scf_threshold, m_scf_energy_threshold,
                                          m_scf_tight_now ? m_scf_gradient_threshold_tight : m_scf_gradient_threshold);
    // screening of the device direct builder (early phase / eri_screening, primitive mode) and the
    // integral precision (Q-G5: FP32 phase / FP64)
    auto setJKScreening = [&](bool early) {
        if (direct && !jk->setFP32(cadence.fp32)) cadence.fp32 = false;  // builder without FP32: FP64 throughout
        if (direct && !jk->setScreening(early ? m_scf_early_screening : m_eri_screening, m_eri_prim_density)) {
            cadence.early = false;  // builder without the hook: final threshold throughout
            return false;
        }
        return true;
    };

    ScfVerifiedStop vstop(m_scf_verified_stop && use_diis && !conv.density_only && !m_scf_tight_now);  // as the host loop (Oct 27 2026)

    for (int iter = 0; iter < m_scf_max_iter; ++iter) {
        m_scf_iterations = iter + 1;

        // Fock (host: fockFromDensity)
        const bool full = cadence.nextFull() || vstop.force_full;
        vstop.force_full = false;
        setJKScreening(cadence.early);
        const bool fp32_build = cadence.fp32;
        if (fp32_build) ++cadence.fp32_builds;
        const auto t_f0 = std::chrono::steady_clock::now();
        if (!dev.buildFock(jk, full)) return fail("Fock build");
        cadence.built(full);
        if (direct) m_direct_quartets += jk->lastQuartets();
        const double fock_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_f0).count();
        t_fock_ms += fock_ms;

        // DIIS / ADIIS: history on the device, subspace solve here
        const auto t_d0 = std::chrono::steady_clock::now();
        int diis_slot = -1;  // slot holding this build's error vector (convergence test), -1 = none
        if (use_diis && iter >= m_diis_start) {
            const int slot = diis.push(cadence.since_full == 0 && !fp32_build);  // FP32 entries leave at the switch
            diis_slot = slot;
            double err = 0.0;
            if (!dev.diisPush(slot, err)) return fail("DIIS error vector");
            diis.setLastError(err);
            const std::vector<int>& sl = diis.slots();
            const int m = static_cast<int>(sl.size());
            Eigen::MatrixXd G, M;
            Eigen::VectorXd a;
            if (diis.needsDiis()) {
                std::vector<double> B((size_t)m * m);
                if (!dev.diisErrorGram(sl.data(), m, B.data())) return fail("DIIS Gram matrix");
                G = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(B.data(), m, m);
            }
            if (diis.needsAdiis()) {
                std::vector<double> av((size_t)m), Mv((size_t)m * m);
                if (!dev.adiisData(sl.data(), m, av.data(), Mv.data())) return fail("ADIIS model");
                a = Eigen::Map<const Eigen::VectorXd>(av.data(), m);
                M = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(Mv.data(), m, m);
            }
            const Eigen::VectorXd c = diis.coefficients(G, a, M);
            if (!dev.selectFock(sl.data(), c.data(), m)) return fail("DIIS extrapolation");
        } else if (!dev.selectFock(nullptr, nullptr, 0)) {
            return fail("Fock selection");
        }
        // convergence scalars of this build (host loop: conv.observe), Claude Generated (Oct 25 2026)
        if (!conv.density_only) {
            double e_n = 0.0, g_n = 0.0;
            if (!dev.scfMeasures(diis_slot, e_n, g_n)) return fail("SCF energy / orbital gradient");
            conv.observe(e_n, g_n);
        }
        // Verified stop (host loop: the same test): the current device P, its J/K and the orbitals of the last
        // eigensolve are the result.
        if (vstop.accept(conv, cadence.since_full == 0 && !cadence.early && !cadence.fp32)) {
            m_scf_last_dE = conv.dE;
            m_scf_last_grad = conv.grad;
            cadence.fp32 = false;
            setJKScreening(false);
            Matrix P(n, n), C(n, n), J(n, n), K(n, n);
            Vector eps(n);
            if (!dev.downloadDensity(false, P.data()) || !dev.downloadOrbitals(C.data(), eps.data()))
                return fail("result download");
            const bool have_jk = dev.downloadJK(J.data(), K.data());  // false: finishConvergedSCF builds J/K
            m_density = P;
            finishConvergedSCF(C, eps, n_occ, iter + 1, vstop.dP_prev, direct,
                               std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scf0).count(),
                               t_fock_ms, t_eig_ms,
                               fmt::format("{} on {} (state setup {:.1f} ms)", dev.stats(), dev.name(), setup_ms),
                               have_jk ? &J : nullptr, have_jk ? &K : nullptr);
            return 1;
        }

        const double diis_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_d0).count();

        // eigensolve, density, convergence scalars
        double delta_P = 0.0, e_elec = 0.0;
        const auto t_e0 = std::chrono::steady_clock::now();
        if (!dev.solveAndDensity(delta_P, e_elec)) return fail("eigensolve / density");
        const double eig_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_e0).count();
        t_eig_ms += eig_ms;

        // printed energy: E(P_n) of this build (auto), as the host loop; density mode keeps Tr(P_new (H + F)) / 2
        const double e_print = conv.density_only ? e_elec : conv.e_prev;
        const std::string conv_txt = conv.density_only ? std::string()
                                                       : fmt::format(", dE = {:.2e}, |g| = {:.2e}", iter > 0 ? conv.dE : 0.0, conv.grad);
        if (CurcumaLogger::get_verbosity() == 1 || CurcumaLogger::get_verbosity() == 2)
            CurcumaLogger::result(fmt::format("SCF iter {:3d}: E_elec = {:.10f} Eh{}, dP = {:.3e} (device, Fock {:.0f} ms {}{})",
                                              iter + 1, e_print, conv_txt, delta_P, fock_ms, full ? "full" : "incremental",
                                              fp32_build ? ", FP32" : ""));
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::param(fmt::format("SCF iter {}", iter + 1),
                                 fmt::format("E_elec = {:.10f} Eh{}, dP = {:.3e} (device; Fock {:.2f} ms ({}{}), DIIS {:.2f} ms, "
                                             "eigensolve + density {:.2f} ms)",
                                             e_print, conv_txt, delta_P, fock_ms, full ? "full" : "incremental",
                                             direct ? fmt::format(", {} quartets", jk->lastQuartets()) : std::string(),
                                             diis_ms, eig_ms));

        if (conv.converged(delta_P) && !cadence.early && !cadence.fp32) {
            m_scf_last_dE = conv.dE;
            m_scf_last_grad = conv.grad;
            cadence.fp32 = false;
            setJKScreening(false);
            Matrix P(n, n), C(n, n);
            Vector eps(n);
            if (!dev.downloadDensity(true, P.data()) || !dev.downloadOrbitals(C.data(), eps.data()))
                return fail("result download");
            m_density = P;
            finishConvergedSCF(C, eps, n_occ, iter + 1, delta_P, direct,
                               std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scf0).count(),
                               t_fock_ms, t_eig_ms,
                               fmt::format("{} on {} (state setup {:.1f} ms)", dev.stats(), dev.name(), setup_ms));
            return 1;
        }
        if (cadence.endEarly(delta_P, m_scf_early_switch)) {
            reportEarlyEnd(iter, delta_P);
            conv.resetEnergy();
        }
        if (cadence.endFP32(delta_P, gpuFP32Switch())) {
            reportFP32End(iter, delta_P, cadence.fp32_builds);
            diis.dropIncrementalEntries();  // Fock matrices of FP32 builds leave the subspace
            conv.resetEnergy();
        }
        if (cadence.observe(delta_P)) {
            reportIncrementalStall(iter, cadence.best_dP, m_eri_screening);
            diis.dropIncrementalEntries();  // keep only Fock matrices from full builds
        }
        vstop.afterStep(conv, delta_P, direct && !cadence.early && !cadence.fp32);
        if (!dev.acceptDensity(damping)) return fail("density update");
    }

    // Ran out of iterations: as the host loop, last density/Fock for diagnostics.
    cadence.fp32 = false;
    setJKScreening(false);
    m_direct_thr_now = -1.0;
    Matrix P(n, n);
    if (!dev.downloadDensity(false, P.data())) return fail("density download");
    m_density = P;
    m_fock = buildFock(m_density);
    m_e_elec = 0.5 * m_density.cwiseProduct(m_H + m_fock).sum();
    return 0;
}

// =================================================================================
// SCF initial guess
// =================================================================================

// Atom index of each active basis function. The cartesian -> spherical transform is
// block-diagonal per shell (Q's column j has nonzeros only in its own shell's
// cartesian rows), so every active AO belongs to exactly one atom.
std::vector<int> QMEngine::activeAtomIndex() const
{
    const int n = m_nbf;
    std::vector<int> ao_atom(n, -1);
    const int ncart = static_cast<int>(m_gto_basis.size());
    if (m_Q.size() == 0) {
        for (int i = 0; i < n && i < ncart; ++i) ao_atom[i] = m_gto_basis[i].atom;
        return ao_atom;
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < ncart; ++i) {
            if (std::abs(m_Q(i, j)) > 1.0e-10) { ao_atom[j] = m_gto_basis[i].atom; break; }
        }
    }
    return ao_atom;
}

// Claude Generated: SAD (superposition of atomic densities) initial guess.
// For every atom, diagonalize that atom's block of the core Hamiltonian in its own
// atomic basis and fill the atom's electrons into the resulting atomic orbitals
// (aufbau; the last orbital is fractionally occupied when Z is odd). Summing the
// atomic densities gives Tr(P S) = N and a start far closer to the physical
// solution than the bare core guess -- which is what keeps the HF SCF from settling
// on a secondary solution (BH is the case in the validation set).
Matrix QMEngine::buildAtomicGuess() const
{
    const int n = m_nbf;
    const std::vector<int> ao_atom = activeAtomIndex();
    Matrix P = Matrix::Zero(n, n);
    double occupied = 0.0;

    const Eigen::MatrixXd Hcol(m_H);   // RowMajor -> ColMajor for Eigen
    const Eigen::MatrixXd Scol(m_S);

    for (int a = 0; a < m_atomcount; ++a) {
        std::vector<int> idx;
        for (int i = 0; i < n; ++i)
            if (ao_atom[i] == a) idx.push_back(i);
        const int m = static_cast<int>(idx.size());
        if (m == 0) continue;

        Eigen::MatrixXd Hs(m, m), Ss(m, m);
        for (int p = 0; p < m; ++p)
            for (int q = 0; q < m; ++q) {
                Hs(p, q) = Hcol(idx[p], idx[q]);
                Ss(p, q) = Scol(idx[p], idx[q]);
            }
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Hs, Ss);
        if (ges.info() != Eigen::Success) continue;   // skip: keep the atom's block empty
        const Eigen::MatrixXd C = ges.eigenvectors();

        // Aufbau over this atom's orbitals with fractional occupation of the last.
        double z = static_cast<double>(m_atoms[a]);
        for (int k = 0; k < m && z > 0.0; ++k) {
            const double f = std::min(2.0, z);
            z -= f;
            occupied += f;
            for (int p = 0; p < m; ++p)
                for (int q = 0; q < m; ++q)
                    P(idx[p], idx[q]) += f * C(p, k) * C(q, k);
        }
    }

    // Charged systems: scale so the guess carries exactly the right electron count.
    const double nelec = static_cast<double>(m_num_electrons);
    if (occupied > 1.0e-12 && std::abs(occupied - nelec) > 1.0e-12)
        P *= nelec / occupied;
    return P;
}

// Claude Generated (Oct 2026): "atom" guess -- superposition of atomic RHF densities.
// The "sad" guess above fills each atom's electrons into the eigenvectors of its block of the
// MOLECULAR core Hamiltonian (no electron repulsion, attraction of all nuclei), which for long
// chains starts the SCF far from the solution (C40H82: E_elec(iter 1) 9 % above the converged value,
// ||dP|| above 50 for the first 9 iterations). Here every element gets one small SCF of the free atom
// in its own basis (core Hamiltonian of its own nucleus only, exact J/K of the atomic basis):
//   F = H_atom + J(P) - 1/2 K(P),  F C = S C eps,  P = sum_k n_k C_k C_k^T,
// with spherically averaged fractional occupations: orbitals with the same eigenvalue (within
// 1e-5 Eh) form one group and share the electrons left for that group equally (e.g. C: 1s^2 2s^2
// 2p^(2/3,2/3,2/3)), which keeps the density spherical. This is the usual SAD of molecular codes
// (e.g. Van Lenthe et al., J. Comput. Chem. 27, 926 (2006); PySCF init_guess "atom"). Only the start
// density changes; the SCF fixed point is the same.
Matrix QMEngine::buildAtomicHFGuess() const
{
    const int n = m_nbf;
    const std::vector<int> ao_atom = activeAtomIndex();
    const int ncart = static_cast<int>(m_gto_basis.size());
    Matrix P = Matrix::Zero(n, n);
    double occupied = 0.0;
    std::map<int, Eigen::MatrixXd> cache;  // element -> atomic density in the atom's active basis

    for (int a = 0; a < m_atomcount; ++a) {
        std::vector<int> idx, cidx;
        for (int i = 0; i < n; ++i)
            if (ao_atom[i] == a) idx.push_back(i);
        for (int i = 0; i < ncart; ++i)
            if (m_gto_basis[i].atom == a) cidx.push_back(i);
        const int m = static_cast<int>(idx.size());
        if (m == 0) continue;
        const int Z = m_atoms[a];
        auto it = cache.find(Z);
        if (it == cache.end() || it->second.rows() != m) {
            std::vector<GTO::Orbital> sub;
            for (int i : cidx) sub.push_back(m_gto_basis[i]);
            Matrix pos(1, 3);
            for (int k = 0; k < 3; ++k) pos(0, k) = CurcumaUnit::Length::angstrom_to_bohr(m_geometry(a, k));
            const std::vector<int> zv{ Z };
            Matrix Sa = qmint::buildOverlap(sub);
            Matrix Ha = qmint::buildKinetic(sub) + qmint::buildNuclearAttraction(sub, zv, pos);
            Matrix Qa;
            if (m_Q.size() != 0) {
                Qa = Matrix::Zero((int)cidx.size(), m);
                for (int r = 0; r < (int)cidx.size(); ++r)
                    for (int c = 0; c < m; ++c) Qa(r, c) = m_Q(cidx[r], idx[c]);
                Sa = qmint::applySphericalTransform(Sa, Qa);
                Ha = qmint::applySphericalTransform(Ha, Qa);
            }
            if (Sa.rows() != m) return buildAtomicGuess();  // basis bookkeeping mismatch: old guess
            const qmint::ERITensor eri = qmint::buildERI(sub, 1, 0.0, m_Q.size() != 0 ? &Qa : nullptr);
            const Eigen::MatrixXd Scol(Sa);
            Eigen::MatrixXd Pa = Eigen::MatrixXd::Zero(m, m);
            Eigen::MatrixXd Pold = Pa;
            for (int it2 = 0; it2 < 200; ++it2) {
                const Matrix Pr(Pa);
                const Matrix Fa = Ha + qmint::buildCoulomb(eri, Pr, 1) - 0.5 * qmint::buildExchange(eri, Pr, 1);
                Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Eigen::MatrixXd(Fa), Scol);
                if (ges.info() != Eigen::Success) break;
                const Eigen::VectorXd e = ges.eigenvalues();
                const Eigen::MatrixXd C = ges.eigenvectors();
                Eigen::MatrixXd Pn = Eigen::MatrixXd::Zero(m, m);
                double left = static_cast<double>(Z);
                for (int k = 0; k < m && left > 1.0e-12;) {
                    int g = 1;
                    while (k + g < m && std::abs(e(k + g) - e(k)) < 1.0e-5) ++g;
                    const double per = std::min(2.0, left / g);
                    for (int q = k; q < k + g; ++q) Pn += per * C.col(q) * C.col(q).transpose();
                    left -= per * g;
                    k += g;
                }
                // damped update (fractional occupations can oscillate between groups)
                Pa = it2 == 0 ? Pn : 0.5 * (Pn + Pa);
                if (it2 > 0 && (Pa - Pold).norm() < 1.0e-9) break;
                Pold = Pa;
            }
            cache[Z] = Pa;
            it = cache.find(Z);
        }
        const Eigen::MatrixXd& Pa = it->second;
        for (int p = 0; p < m; ++p)
            for (int q = 0; q < m; ++q) P(idx[p], idx[q]) = Pa(p, q);
        occupied += static_cast<double>(Z);
    }
    const double nelec = static_cast<double>(m_num_electrons);
    if (occupied > 1.0e-12 && std::abs(occupied - nelec) > 1.0e-12)
        P *= nelec / occupied;
    return P;
}

// Claude Generated (Sep 2026): warm start for the next geometry of the same molecule.
// The basis functions moved with the atoms, so the old occupied orbitals C are no
// longer orthonormal in the new metric S. Symmetric (Loewdin) re-orthonormalisation
//   C' = C (C^T S C)^{-1/2},   P = 2 C' C'^T
// restores Tr(P S) = N and idempotency (P S P = 2 P) exactly, and is the closest
// orthonormal set to C -- a much better start than superposed atomic densities once
// the steps are small (optimisation, MD).
Matrix QMEngine::buildWarmStartGuess() const
{
    const int n_occ = m_num_electrons / 2;
    if (!m_scf_warm_start || m_warm_C.size() == 0 || m_warm_atoms != m_atoms
        || m_warm_C.rows() != m_nbf || m_warm_C.cols() != n_occ || n_occ == 0)
        return Matrix();
    const Eigen::MatrixXd C(m_warm_C);
    const Eigen::MatrixXd M = C.transpose() * Eigen::MatrixXd(m_S) * C;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() < 1.0e-8)
        return Matrix();  // near-singular: the geometry step was too large, fall back
    const Eigen::MatrixXd Minvsqrt = es.eigenvectors()
        * es.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() * es.eigenvectors().transpose();
    const Eigen::MatrixXd Cp = C * Minvsqrt;
    return Matrix(2.0 * Cp * Cp.transpose());
}

// Claude Generated (Oct 27 2026): -qm.scf_extrapolation aspc|gauss. P_pred = sum_j w_j P_j over the converged
// densities of the last geometries (weights of src/core/charge_extrapolation.h, sum w_j = 1, as the xtb SCC
// extrapolation), then projected onto an idempotent closed-shell density in the CURRENT metric: natural orbitals of
// M = (S X)^T P_pred (S X) (X = S^-1/2, so S X = S^1/2), the n_occ largest occupations kept, C = X V_occ,
// P = 2 C C^T. Empty (warm start instead) when the history is too short, the molecule or basis changed, or the
// predicted occupations have no clear gap (lambda_occ - lambda_virt < 1; exact: 2 and 0).
Matrix QMEngine::buildExtrapolatedGuess() const
{
    const int n_occ = m_num_electrons / 2;
    if (m_scf_extrapolation == "none" || !m_scf_warm_start || m_dens_hist.size() < 2 || m_warm_atoms != m_atoms
        || n_occ == 0 || m_X.size() == 0 || m_dens_hist_serial == m_geom_serial)
        return Matrix();
    for (const Matrix& P : m_dens_hist)
        if (P.rows() != m_nbf || P.cols() != m_nbf) return Matrix();
    const Eigen::VectorXd w = curcuma::extrapolation::weights(m_scf_extrapolation, m_scf_extrap_order,
                                                              static_cast<int>(m_dens_hist.size()));
    if (w.size() < 2) return Matrix();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(m_nbf, m_nbf);
    for (int j = 0; j < w.size(); ++j) P.noalias() += w(j) * Eigen::MatrixXd(m_dens_hist[j]);
    if (!P.allFinite()) return Matrix();
    const Eigen::MatrixXd Xc(m_X), Sc(m_S);
    const Eigen::MatrixXd SX = gemmThreaded(Sc, false, Xc, m_threads);
    const Eigen::MatrixXd PSX = gemmThreaded(P, false, SX, m_threads);
    Eigen::MatrixXd M = gemmThreaded(SX, true, PSX, m_threads);
    M = 0.5 * (M + M.transpose());
    Eigen::VectorXd lam;
    Eigen::MatrixXd V;
    if (!curcuma::eigsolver::solveSymmetric(M, lam, V, m_threads)) return Matrix();
    const int m = static_cast<int>(lam.size());
    if (m <= n_occ) return Matrix();
    const double lam_occ = lam(m - n_occ), lam_virt = lam(m - n_occ - 1);  // ascending order
    if (!(lam_occ - lam_virt > 1.0)) return Matrix();
    const Eigen::MatrixXd Vo = V.rightCols(n_occ);
    const Eigen::MatrixXd CoT = gemmThreaded(Vo, true, Xc.transpose(), m_threads);  // (X V_occ)^T
    const Eigen::MatrixXd Ph = gemmThreaded(CoT, true, CoT, m_threads);
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::param("SCF guess", fmt::format("{} extrapolation over {} densities (natural occupations: "
                                                      "lowest occupied {:.4f}, highest virtual {:.4f})",
                                                      m_scf_extrapolation, static_cast<int>(w.size()), lam_occ, lam_virt));
    CurcumaLogger::citation(m_scf_extrapolation == "aspc" ? "aspc" : "density_extrapolation");
    return Matrix(Ph + Ph.transpose());  // = 2 C_occ C_occ^T
}

Matrix QMEngine::buildInitialGuess() const
{
    m_last_warm_start = false;
    if (m_scf_continue && m_density.rows() == m_nbf && m_density.cols() == m_nbf) {
        m_last_warm_start = true;  // tightening a converged SCF of this geometry: its density (Oct 25 2026)
        return m_density;
    }
    const Matrix predicted = buildExtrapolatedGuess();  // -qm.scf_extrapolation aspc|gauss (Oct 27 2026)
    if (predicted.size() != 0) {
        m_last_warm_start = true;
        return predicted;
    }
    const Matrix warm = buildWarmStartGuess();
    if (warm.size() != 0) {
        m_last_warm_start = true;
        return warm;
    }
    if (m_scf_guess == "h0")
        return Matrix::Zero(m_nbf, m_nbf);   // bare core: first Fock = H
    if (m_scf_guess == "atom")
        return buildAtomicHFGuess();
    if (m_scf_guess != "sad" && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format(
            "QM: unknown -qm.scf_guess '{}' (use sad|atom|h0); using sad", m_scf_guess));
    return buildAtomicGuess();
}

// =================================================================================
// Fock build: F = H + J - 1/2 K  (closed-shell RHF, spin-summed density)
// =================================================================================

Matrix QMEngine::buildFock(const Matrix& P) const
{
    Matrix J, K;
    buildJKForSCF(P, J, K);
    return composeFock(P, J, K);
}

// =================================================================================
// WP5 Kohn-Sham pieces (Claude Generated, Oct 2026): F = H + J - (a_x/2) K + V_xc, E with E_xc
// =================================================================================

void QMEngine::invalidateXC() const
{
    m_xc_grid_ready = false;
    m_xc_grid = qmgrid::MolecularGrid();
    m_xc_ao = qmgrid::AOEvaluator();
    m_xc_last_valid = false;
    m_xc_last_P = Matrix();
    m_xc_failed = false;
}

bool QMEngine::ensureXCGrid() const
{
    if (m_xc_grid_ready) return true;
    std::string why;
    if (m_xc_frozen_grid) {
        m_xc_grid = *m_xc_frozen_grid;  // WP8 validation hook: grid fixed in space
    } else if (!buildGrid(m_xc_grid, why)) {
        m_xc_failed = true;
        CurcumaLogger::error("QM: DFT grid could not be built: " + why);
        return false;
    }
    m_xc_ao = qmgrid::AOEvaluator(m_gto_basis, m_Q);  // AO values in the active basis (basis centres of this geometry)
    m_xc_grid_ready = true;
    m_xc_last_valid = false;
    m_xc_builds = 0;
    m_xc_ms = 0.0;
    return true;
}

// WP8 (Claude Generated, Oct 2026): validation hook, grid fixed in space (see qm_engine.h).
void QMEngine::freezeXCGrid(const qmgrid::MolecularGrid& grid)
{
    m_xc_frozen_grid = std::make_shared<const qmgrid::MolecularGrid>(grid);
    invalidateXC();
}

// WP8 (Claude Generated, Oct 2026): XC part of the Kohn-Sham nuclear gradient of the converged density (Eh/Bohr).
//   g_basis    = the basis-function part (qmxc::buildXCGradient), the whole XC gradient of a grid fixed in space;
//   g_response = with -qm.grid_response: point motion (qmxc XCGradient::points) + Becke weight derivative
//                (qmgrid::gridWeightGradient), so that g_basis + g_response is the derivative of the energy whose
//                grid moves with the atoms (the energy curcuma reports); zero otherwise.
bool QMEngine::xcGradient(Matrix& g_basis, Matrix& g_response, std::string& why) const
{
    g_response = Matrix::Zero(m_atomcount, 3);
    if (!ensureXCGrid()) {
        why = "no DFT grid";
        g_basis = Matrix::Zero(m_atomcount, 3);
        return false;
    }
    std::vector<int> bf_atom(m_nbf, 0);
    if (m_Q.size() != 0) {
        for (int a = 0; a < m_nbf; ++a)
            for (int c = 0; c < m_Q.rows(); ++c)
                if (m_Q(c, a) != 0.0) {
                    bf_atom[a] = m_gto_basis[c].atom;
                    break;
                }
    } else {
        for (int a = 0; a < m_nbf; ++a) bf_atom[a] = m_gto_basis[a].atom;
    }
    qmxc::XCGradient xg;
    qmxc::buildXCGradient(m_xc_grid, m_xc_ao, m_density, m_xc, m_xc_rho_threshold, m_threads, bf_atom, m_atomcount, xg);
    g_basis = xg.basis;
    double wms = 0.0;
    if (m_grid_response) {
        if (m_xc_frozen_grid) {
            why = "grid response with a frozen grid (validation hook) is undefined; use -qm.grid_response false";
            return false;
        }
        const auto t0 = std::chrono::steady_clock::now();
        Matrix gw = Matrix::Zero(m_atomcount, 3);
        if (!qmgrid::gridWeightGradient(m_atoms, atomPositionsBohr(), m_grid_options, m_xc_grid, xg.edens, gw, why,
                                        m_threads))
            return false;
        g_response = xg.points + gw;
        wms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("QM: XC gradient {:.1f} ms ({} points above threshold), weight derivative "
                                        "{:.1f} ms; |g_xc basis| {:.3e}, |g_response| {:.3e} Eh/Bohr",
                                        xg.ms, xg.points_used, wms, g_basis.norm(), g_response.norm()));
    return true;
}

const qmxc::XCResult& QMEngine::xcFor(const Matrix& P) const
{
    if (m_xc_last_valid && m_xc_last_P.rows() == P.rows() && m_xc_last_P.cols() == P.cols()
        && (m_xc_last_P.array() == P.array()).all())
        return m_xc_last;
    if (!ensureXCGrid()) {
        // never a silent zero V_xc: NaN makes the SCF fail loudly
        m_xc_last = qmxc::XCResult();
        m_xc_last.V = Matrix::Constant(m_nbf, m_nbf, std::numeric_limits<double>::quiet_NaN());
        m_xc_last.exc = std::numeric_limits<double>::quiet_NaN();
        m_xc_last_valid = false;
        return m_xc_last;
    }
    qmxc::buildVxc(m_xc_grid, m_xc_ao, P, m_xc, m_xc_rho_threshold, m_threads, m_xc_last);
    m_xc_last_P = P;
    m_xc_last_valid = true;
    ++m_xc_builds;
    m_xc_ms += m_xc_last.ms;
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("QM: V_xc build {}: E_xc = {:.10f} Eh (E_x {:.10f}, E_c {:.10f}), N_grid = {:.8f}, "
                                        "{:.1f} ms", m_xc_builds, m_xc_last.exc, m_xc_last.ex, m_xc_last.ec,
                                        m_xc_last.nelec, m_xc_last.ms));
    return m_xc_last;
}

Matrix QMEngine::composeFock(const Matrix& P, const Matrix& J, const Matrix& K) const
{
    if (!isKohnSham()) return m_H + J - 0.5 * K;  // HF: unchanged expression
    Matrix F = m_H + J + xcFor(P).V;
    if (m_xc.exact_exchange != 0.0) F -= (0.5 * m_xc.exact_exchange) * K;  // hybrids (WP7)
    return F;
}

double QMEngine::scfEnergy(const Matrix& P, const Matrix& F) const
{
    const double e = 0.5 * P.cwiseProduct(m_H + F).sum();
    if (!isKohnSham()) return e;
    // 1/2 Tr(P(H+F)) counts 1/2 Tr(P V_xc) instead of E_xc: E = Tr(Ph) + 1/2 Tr(PJ) - a_x/4 Tr(PK) + E_xc
    const qmxc::XCResult& xc = xcFor(P);
    return e - 0.5 * P.cwiseProduct(xc.V).sum() + xc.exc;
}

void QMEngine::buildJKForSCF(const Matrix& P, Matrix& J, Matrix& K) const
{
    if (!isKohnSham() || m_xc.exact_exchange != 0.0) {
        buildJK(P, J, K);
        return;
    }
    buildCoulombOnly(P, J);
    K = Matrix::Zero(m_nbf, m_nbf);
}

// J(P) alone (pure density functionals). Stored tensor: buildCoulomb only (the exchange contraction is skipped);
// RI (-qm.ri j|jk): J from the fit, no exchange; device stored tensor: the J/K kernel, K discarded; integral-direct
// (CPU DirectJK or the device builder): J/K digest of the same quartets, K discarded -- a J-only direct kernel with
// its own (looser) screen does not exist yet. RI mixed precision keeps the full buildJK (its stage machine counts
// exact builds there).
void QMEngine::buildCoulombOnly(const Matrix& P, Matrix& J) const
{
    if (riActive()) {
        if (m_ri_mp_stage != 0) {
            Matrix K;
            buildJK(P, J, K);
            return;
        }
        if (!ensureRI()) {
            J = Matrix::Constant(m_nbf, m_nbf, std::numeric_limits<double>::quiet_NaN());
            return;
        }
        const auto t0 = std::chrono::steady_clock::now();
        if (!gpuRIBuild(P, J, nullptr)) J = m_ri->coulomb(P);
        m_ri_build_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        ++m_ri_builds;
        return;
    }
    if (!useDirectSCF()) {
        Matrix K;
        if (ensureGpuJK() && gpuBuildJK(P, J, K)) return;
        if (!useDirectSCF()) {
            J = qmint::buildCoulomb(eriActive(), P, m_threads);
            return;
        }
    }
    Matrix K;
    buildJKExact(P, J, K);
}

// Claude Generated (Oct 2026): the early loose-screening phase applies to a cold-started
// integral-direct SCF only (a warm start begins below any sensible switch threshold).
bool QMEngine::useEarlyScreening() const
{
    return m_scf_early_screening > 0.0 && m_scf_early_screening > m_eri_screening && !m_last_warm_start;
}

// Claude Generated (Oct 2026, Q-G5): -qm.gpu_eri_precision fp64|fp32|mixed|auto. True when this SCF
// starts with FP32 device direct builds: cold start, exact J/K integral-direct on the device, a
// builder with FP32 kernels, and for auto a device whose FP64 rate is <= 1/16 of FP32.
bool QMEngine::useGpuFP32Phase(bool direct) const
{
    std::string mode = m_gpu_eri_precision;
    std::transform(mode.begin(), mode.end(), mode.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    if (mode != "fp64" && mode != "fp32" && mode != "mixed" && mode != "auto") {
        static bool warned = false;
        if (!warned && CurcumaLogger::get_verbosity() >= 1) {
            warned = true;
            CurcumaLogger::warn(fmt::format("QM: unknown -qm.gpu_eri_precision '{}' (use fp64|fp32|mixed|auto); using auto",
                                            m_gpu_eri_precision));
        }
        mode = "auto";
    }
    if (mode == "fp64" || !direct || m_last_warm_start || riActive()) return false;
    const bool explicit_fp32 = mode != "auto";
    auto unavailable = [&](const std::string& why) {
        static bool warned = false;
        if (explicit_fp32 && !warned && CurcumaLogger::get_verbosity() >= 1) {
            warned = true;
            CurcumaLogger::warn("QM: -qm.gpu_eri_precision " + mode + " not available (" + why + "); FP64 J/K builds");
        }
        return false;
    };
    if (!gpuDirectAllowed() || !ensureGpuDirect()) return unavailable("needs the device integral-direct J/K");
    // auto: FP64-poor device and not a small basis (benzene, 114 functions, was 16 % slower with the
    // FP32 phase: one extra iteration outweighs the cheaper builds; naphthalene, 180, 10 % faster)
    constexpr int kFP32AutoMinNbf = 150;
    if (mode == "auto" && (m_gpu_djk->fp32Ratio() < 16 || m_nbf < kFP32AutoMinNbf)) return false;
    if (!m_gpu_djk->setFP32(true)) return unavailable("no FP32 kernels for this basis on " + m_gpu_djk->name());
    m_gpu_djk->setFP32(false);
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format("QM: FP32 integrals until ||dP|| < {:.1e} (-qm.gpu_eri_precision {}, FP32:FP64 rate {} on {})",
                                        gpuFP32Switch(), mode, m_gpu_djk->fp32Ratio(), m_gpu_djk->name()));
    return true;
}

double QMEngine::gpuFP32Switch() const
{
    std::string mode = m_gpu_eri_precision;
    std::transform(mode.begin(), mode.end(), mode.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return mode == "fp32" ? std::min(m_gpu_fp32_switch, 1.0e-5) : m_gpu_fp32_switch;
}

// Claude Generated (Sep 2026): -qm.scf_direct on|off|auto. auto keeps the stored
// tensor (fastest J/K) as long as it fits in -qm.eri_max_memory_mb.
bool QMEngine::useDirectSCF() const
{
    // Q-G3: -qm.gpu_direct on forces integral-direct J/K on the device for every -gpu run.
    if (m_gpu_direct == "on" && gpuRequested() && !m_gpu_direct_failed) return true;
    if (m_scf_direct == "on" || m_scf_direct == "true") return true;
    if (m_scf_direct == "off" || m_scf_direct == "false") return false;
    if (m_scf_direct != "auto" && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.scf_direct '{}' (use auto|on|off); using auto", m_scf_direct));
    if (hostTensorAllowed()) return false;
    // Q-G2: a tensor built on the device needs no host memory, so auto keeps the stored
    // (device) tensor as long as it fits into device memory.
    if (gpuRequested() && m_gpu_eri && ensureGpuJK()) return false;
    return true;
}

// Claude Generated (Oct 2026, Q-G1/Q-G2): J/K on the device from a device-resident tensor.
// Q-G2 (default, -qm.gpu_eri true): the tensor is computed on the device from the flat
// shell-pair tables (qmint::buildEriGpuTables) and never exists on the host.
// Q-G1 (fallback): CPU tensor, uploaded once per geometry (only if the host may hold it).
// m_gpu_jk is reset with the ERI caches. Every failure falls back to the CPU and is
// counted (gpu_fallback): "qm_eri" for the device build, "qm_jk" for the J/K path.
bool QMEngine::gpuRequested() const
{
    return !(m_gpu == "none" || m_gpu == "cpu" || m_gpu.empty() || m_gpu_failed);
}

bool QMEngine::hostTensorAllowed() const
{
    if (m_scf_direct == "off" || m_scf_direct == "false") return true;
    const double n = static_cast<double>(m_nbf);
    return n * n * n * n * 8.0 / 1.0e6 <= m_eri_max_memory_mb;
}

bool QMEngine::ensureGpuJK() const
{
    if (!gpuRequested()) return false;
    if (m_gpu_jk) return true;
    const std::string backend = m_gpu;
    if (backend == "none" || !gpu_plugin::available(backend)) {
        m_gpu_failed = true;
        curcuma::reportGpuFallback("qm_jk", "no GPU plugin for -gpu " + m_gpu);
        return false;
    }
    if (m_gpu_eri && !m_gpu_eri_failed) {
        const auto t0 = std::chrono::steady_clock::now();
        qmint::EriGpuTables tables;
        std::string why;
        bool ok = true;
        const double tensor_mb = std::pow((double)m_nbf, 4) * 8.0 / 1.0e6;
        // Fit check BEFORE the host tables (Claude Generated, Oct 2026): the tables of the stored
        // device tensor include the O(npair^2) quartet list (~8 s at 970 and ~160 s at 1930 basis
        // functions) and were built even when the n^4 tensor was larger than the whole card. The
        // card's total memory is a static device property (no context); the plugin's own margin
        // (3 n^2 matrices + 256 MB) is added. Below the total, the free memory is still checked by
        // the plugin after the tables exist, as before.
        const json info = gpu_plugin::deviceInfo(backend, m_gpu_device >= 0 ? m_gpu_device : 0);
        const double card_mb = info.contains("memory_total_bytes") && info["memory_total_bytes"].is_number()
            ? info["memory_total_bytes"].get<double>() / 1.0e6
            : 0.0;
        const double need_mb = tensor_mb + 3.0 * (double)m_nbf * m_nbf * 8.0 / 1.0e6 + 256.0 * 1.048576;
        if (m_gpu_max_memory_mb > 0.0 && tensor_mb > m_gpu_max_memory_mb) {
            // -qm.gpu_max_memory_mb: decided on the host, no tables built
            ok = false;
            why = fmt::format("ERI tensor {:.1f} MB above -qm.gpu_max_memory_mb {:.1f}", tensor_mb, m_gpu_max_memory_mb);
        } else if (card_mb > 0.0 && need_mb > card_mb) {
            ok = false;
            why = fmt::format("ERI tensor {:.1f} MB does not fit the device memory {:.1f} MB (decided before building "
                              "the host tables)", tensor_mb, card_mb);
        } else {
            ok = qmint::buildEriGpuTables(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr, tables, why);
            tables.device_mem_limit_mb = m_gpu_max_memory_mb;
        }
        const auto t1 = std::chrono::steady_clock::now();
        if (ok) m_gpu_jk = gpu_plugin::createQmEriJK(backend, tables.view(), m_gpu_device, why);
        if (m_gpu_jk) {
            if (CurcumaLogger::get_verbosity() >= 2) {
                const auto t2 = std::chrono::steady_clock::now();
                const int n = m_gpu_jk->n();
                CurcumaLogger::info(fmt::format(
                    "QM: built 4-centre ERI on {} ({} basis functions, {:.1f} MB on the device) in {:.1f} ms "
                    "(host tables {:.1f} ms; device: {}); J/K on {}",
                    m_gpu_jk->name(), n, (double)n * n * n * n * 8.0 / 1.0e6,
                    std::chrono::duration<double, std::milli>(t2 - t0).count(),
                    std::chrono::duration<double, std::milli>(t1 - t0).count(), m_gpu_jk->stats(), m_gpu_jk->name()));
            }
            // Diagnostic (test_cases/qm_gpu): element-wise device tensor vs the CPU buildERI.
            if (const char* chk = std::getenv("CURCUMA_QM_GPU_ERI_CHECK"); chk && chk[0] == '1')
                checkGpuEri();
            return true;
        }
        m_gpu_eri_failed = true;
        if (!hostTensorAllowed() && gpuDirectAllowed()) {
            // Q-G3: not a CPU fallback -- the SCF goes integral-direct on the device next.
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info("QM: no stored device ERI tensor (" + why + "); J/K integral-direct on the device");
        } else {
            curcuma::reportGpuFallback("qm_eri", "device ERI build not possible: " + why);
        }
    }
    // Q-G1: CPU tensor + upload. Not when the host may not hold the tensor (direct SCF then).
    if (!hostTensorAllowed()) return false;
    const qmint::ERITensor& eri = eriActive();
    const auto t0 = std::chrono::steady_clock::now();
    m_gpu_jk = gpu_plugin::createQmJK(backend, eri.data(), eri.n(), m_gpu_device);
    if (!m_gpu_jk) {
        m_gpu_failed = true;
        curcuma::reportGpuFallback("qm_jk", "device upload of the ERI tensor failed (" + backend + ")");
        return false;
    }
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format("QM: J/K on {} ({} basis functions, ERI uploaded in {:.1f} ms)",
            m_gpu_jk->name(), eri.n(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count()));
    return true;
}

// Claude Generated (Oct 2026, Q-G2): CURCUMA_QM_GPU_ERI_CHECK=1 downloads the device tensor and
// compares it element by element with the CPU tensor (same screening, same 5d/6d basis).
// Prints one ASCII line that test_cases/qm_gpu/check_qm_gpu_eri.py parses. Diagnostic only:
// it costs a full CPU ERI build and n^4 doubles of host memory.
void QMEngine::checkGpuEri() const
{
    if (!m_gpu_jk) return;
    const int n = m_gpu_jk->n();
    std::vector<double> dev((size_t)n * n * n * n);
    if (!m_gpu_jk->downloadERI(dev.data())) {
        CurcumaLogger::error("QM GPU ERI check: download of the device tensor failed");
        return;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const qmint::ERITensor cpu = qmint::buildERI(m_gto_basis, m_threads, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
    const double cpu_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    double max_abs = 0.0, max_ref = 0.0;
    size_t nbad = 0, zero_mismatch = 0;
    const double* c = cpu.data();
    for (size_t k = 0; k < dev.size(); ++k) {
        const double d = std::abs(dev[k] - c[k]);
        max_abs = std::max(max_abs, d);
        max_ref = std::max(max_ref, std::abs(c[k]));
        if (d > 1e-12) ++nbad;
        if ((dev[k] == 0.0) != (c[k] == 0.0)) ++zero_mismatch;
    }
    CurcumaLogger::result(fmt::format(
        "QM GPU ERI check: n={} elements={} max|dev-cpu|={:.3e} max|cpu|={:.3e} above_1e-12={} zero_pattern_mismatch={} cpu_build_ms={:.1f}",
        n, dev.size(), max_abs, max_ref, nbad, zero_mismatch, cpu_ms));
}

bool QMEngine::gpuBuildJK(const Matrix& P, Matrix& J, Matrix& K) const
{
    if (!m_gpu_jk) return false;
    const int n = m_gpu_jk->n();
    const Matrix Pc = P;
    J = Matrix::Zero(n, n);
    K = Matrix::Zero(n, n);
    if (!m_gpu_jk->build(Pc.data(), J.data(), K.data())) {
        m_gpu_failed = true;
        m_gpu_jk.reset();
        curcuma::reportGpuFallback("qm_jk", "device J/K build failed");
        return false;
    }
    return true;
}

// Claude Generated (Oct 2026, Q-G3): integral-direct J/K on the device. Used whenever the SCF
// runs integral-direct (useDirectSCF()) with -gpu and -qm.gpu_direct auto|on. The device
// builder recomputes the screened shell quartets per call exactly like qmint::DirectJK (same
// tables, same density-weighted screen), so the incremental-on-dP scheme of runSCF() applies
// unchanged. Any failure: CPU DirectJK, counted as "qm_direct".
bool QMEngine::gpuDirectAllowed() const
{
    static bool warned = false;
    if (!warned && m_gpu_direct != "auto" && m_gpu_direct != "on" && m_gpu_direct != "off" && CurcumaLogger::get_verbosity() >= 1) {
        warned = true;
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.gpu_direct '{}' (use auto|on|off); using auto", m_gpu_direct));
    }
    return gpuRequested() && m_gpu_direct != "off" && !m_gpu_direct_failed;
}

bool QMEngine::ensureGpuDirect() const
{
    if (m_gpu_djk) return true;
    if (!gpuDirectAllowed()) return false;
    const std::string backend = m_gpu;
    if (backend == "none" || !gpu_plugin::available(backend)) {
        m_gpu_direct_failed = true;
        curcuma::reportGpuFallback("qm_direct", "no GPU plugin for -gpu " + m_gpu);
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    qmint::EriGpuTables tables;
    std::string why;
    const bool ok = qmint::buildEriGpuTables(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr,
                                             tables, why, /*list_quartets=*/false);
    const auto t1 = std::chrono::steady_clock::now();
    if (ok) m_gpu_djk = gpu_plugin::createQmDirectJK(backend, tables.view(), m_gpu_device, why);
    if (!m_gpu_djk) {
        m_gpu_direct_failed = true;
        curcuma::reportGpuFallback("qm_direct", "device integral-direct J/K not possible: " + why);
        return false;
    }
    if (CurcumaLogger::get_verbosity() >= 2) {
        const int n = m_gpu_djk->n();
        CurcumaLogger::info(fmt::format(
            "QM: integral-direct J/K on {} ({} basis functions; stored tensor would need {:.1f} MB; setup {:.1f} ms, "
            "host tables {:.1f} ms)",
            m_gpu_djk->name(), n, std::pow((double)n, 4) * 8.0 / 1.0e6,
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(),
            std::chrono::duration<double, std::milli>(t1 - t0).count()));
    }
    return true;
}

bool QMEngine::gpuDirectBuild(const Matrix& P, Matrix& J, Matrix& K) const
{
    if (!m_gpu_djk) return false;
    const int n = m_gpu_djk->n();
    const Matrix Pc = P;
    J = Matrix::Zero(n, n);
    K = Matrix::Zero(n, n);
    if (!m_gpu_djk->build(Pc.data(), J.data(), K.data())) {
        m_gpu_direct_failed = true;
        m_gpu_djk.reset();
        curcuma::reportGpuFallback("qm_direct", "device integral-direct J/K build failed");
        return false;
    }
    m_direct_quartets += m_gpu_djk->lastQuartets();
    if (!m_gpu_jk_checked) {
        m_gpu_jk_checked = true;
        if (const char* chk = std::getenv("CURCUMA_QM_GPU_JK_CHECK"); chk && chk[0] == '1')
            checkGpuDirectJK(Pc, J, K);
    }
    return true;
}

// Claude Generated (Oct 2026, Q-G3): CURCUMA_QM_GPU_JK_CHECK=1 compares J and K of the FIRST
// device integral-direct build (a full density) with (a) the CPU DirectJK, (b) the stored
// device tensor (Q-G2) if it fits on the card, (c) the stored CPU tensor if the host may hold
// it. One ASCII line, parsed by test_cases/qm_gpu/check_qm_gpu_direct.py. Diagnostic only.
void QMEngine::checkGpuDirectJK(const Matrix& P, const Matrix& J, const Matrix& K) const
{
    auto maxdiff = [](const Matrix& a, const Matrix& b) {
        return (a.rows() == b.rows() && a.cols() == b.cols()) ? (a - b).cwiseAbs().maxCoeff() : -1.0;
    };
    const long dev_q = (long)m_gpu_djk->lastQuartets();
    // (a) CPU integral-direct, same tables and screen
    const auto t0 = std::chrono::steady_clock::now();
    const qmint::DirectJK cpu(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
    Matrix Jc, Kc;
    const long cpu_q = cpu.build(P, Jc, Kc, m_threads, m_direct_thr_now, m_eri_prim_density);
    const double cpu_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::string line = fmt::format(
        "QM GPU JK check: n={} quartets_dev={} quartets_cpu={} max|J|={:.3e} max|K|={:.3e} "
        "dJ_cpu_direct={:.3e} dK_cpu_direct={:.3e} cpu_direct_ms={:.1f}",
        J.rows(), dev_q, cpu_q, J.cwiseAbs().maxCoeff(), K.cwiseAbs().maxCoeff(), maxdiff(J, Jc), maxdiff(K, Kc), cpu_ms);
    // (b) stored device tensor
    const std::string backend = m_gpu;
    qmint::EriGpuTables tables;
    std::string why;
    std::unique_ptr<qmint::QMGpuJK> stored;
    if (qmint::buildEriGpuTables(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr, tables, why))
        stored = gpu_plugin::createQmEriJK(backend, tables.view(), m_gpu_device, why);
    if (stored) {
        const int n = stored->n();
        Matrix Js = Matrix::Zero(n, n), Ks = Matrix::Zero(n, n);
        if (stored->build(P.data(), Js.data(), Ks.data()))
            line += fmt::format(" dJ_dev_stored={:.3e} dK_dev_stored={:.3e}", maxdiff(J, Js), maxdiff(K, Ks));
        else
            line += " dJ_dev_stored=failed";
    } else {
        line += " dJ_dev_stored=n/a(" + why + ")";
    }
    stored.reset();
    // (c) stored CPU tensor
    if (hostTensorAllowed()) {
        const qmint::ERITensor eri = qmint::buildERI(m_gto_basis, m_threads, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
        const Matrix Jh = qmint::buildCoulomb(eri, P, m_threads), Kh = qmint::buildExchange(eri, P, m_threads);
        line += fmt::format(" dJ_cpu_stored={:.3e} dK_cpu_stored={:.3e}", maxdiff(J, Jh), maxdiff(K, Kh));
    } else {
        line += " dJ_cpu_stored=n/a(host limit)";
    }
    CurcumaLogger::result(line);
}

// J(P), K(P) in the active basis: from the stored tensor (device or host), or
// integral-direct (device first with -gpu, Q-G3; shell-pair tables built on the first call
// per geometry).
void QMEngine::buildJKExact(const Matrix& P, Matrix& J, Matrix& K) const
{
    if (!useDirectSCF()) {
        if (ensureGpuJK() && gpuBuildJK(P, J, K)) return;
        // re-evaluated: after a device failure the host may not be allowed to hold the tensor
        if (!useDirectSCF()) {
            const qmint::ERITensor& eri = eriActive();
            J = qmint::buildCoulomb(eri, P, m_threads);
            K = qmint::buildExchange(eri, P, m_threads);
            return;
        }
    }
    if (gpuDirectAllowed() && ensureGpuDirect()) {
        m_gpu_djk->setScreening(m_direct_thr_now >= 0.0 ? m_direct_thr_now : m_eri_screening, m_eri_prim_density);
        m_gpu_djk->setFP32(m_direct_fp32_now);  // Q-G5: FP32 only inside the FP32 phase of the host loop
        if (gpuDirectBuild(P, J, K)) return;
    }
    if (!m_direct.ready()) {
        m_direct = qmint::DirectJK(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "QM: integral-direct J/K ({} basis functions; stored tensor would need {:.1f} MB)",
                m_direct.n(), std::pow((double)m_direct.n(), 4) * 8.0 / 1.0e6));
    }
    m_direct_quartets += m_direct.build(P, J, K, m_threads, m_direct_thr_now, m_eri_prim_density);
}

// =================================================================================
// Generalized eigenproblem via the Lowdin reduce: F C = S C eps
//   F' = X^T F X  (symmetric),  F' C' = C' eps,  C = X C'.
// =================================================================================

bool QMEngine::solveFock(const Matrix& F, Matrix& C, Vector& eps) const
{
    const Eigen::MatrixXd Fcol(F);          // RowMajor -> ColMajor copy
    const Eigen::MatrixXd Xc(m_X);
    const Eigen::MatrixXd Fp = gemmThreaded(gemmThreaded(Xc, true, Fcol, m_threads), false, Xc, m_threads);
    Eigen::VectorXd ev;
    Eigen::MatrixXd evecs;
    if (!curcuma::eigsolver::solveSymmetric(Fp, ev, evecs, m_threads))
        return false;
    eps = ev;                                // ascending orbital energies
    C = gemmThreaded(Xc, false, evecs, m_threads);  // back to the AO basis
    return true;
}

// =================================================================================
// Closed-shell density: P = 2 * sum_{i<n_occ} C_i C_i^T
// =================================================================================

Matrix QMEngine::buildDensity(const Matrix& C, int n_occ) const
{
    // P = 2 C_occ C_occ^T as one (threaded) matrix product; the former sum of n_occ rank-1
    // updates streamed the whole n x n matrix n_occ times (~1 s per call at 2000 functions).
    const Eigen::MatrixXd Co(C.leftCols(n_occ));
    const Eigen::MatrixXd CoT = Co.transpose();
    const Eigen::MatrixXd P = gemmThreaded(CoT, true, CoT, m_threads);
    return P + P.transpose();  // = 2 C_occ C_occ^T, exactly symmetric (blocks may round differently)
}

// =================================================================================
// Lowdin orthonormalizer X = S^{-1/2}
// =================================================================================

void QMEngine::buildOrthonormalizer()
{
    const Eigen::MatrixXd Scol(m_S);         // RowMajor -> ColMajor copy
    Eigen::MatrixXd X;
    if (!lowdinOrthonormalizer(Scol, X)) {
        CurcumaLogger::error("QM SCF: overlap matrix is singular (basis linearly dependent)");
        m_X = Matrix();
        return;
    }
    m_X = X;                                  // store as RowMajor
}

// =================================================================================
// Active-basis ERI (cartesian, or the spherical-5d transform of it)
// =================================================================================

const qmint::ERITensor& QMEngine::eriActive() const
{
    // No d shells / cartesian_d: the cartesian tensor IS the active one -- return
    // it directly instead of keeping a second n^4 copy.
    if (m_Q.size() == 0)
        return cartesianERI();
    if (!m_eri_active_ready) {
        // Spherical: transform every shell quartet as it is computed, so the
        // cartesian tensor (1.7 GB for benzene/def2-SVP) is never formed.
        const auto t0 = std::chrono::steady_clock::now();
        m_eri_active = qmint::buildERI(m_gto_basis, m_threads, m_eri_screening, &m_Q);
        m_eri_active_ready = true;
        if (CurcumaLogger::get_verbosity() >= 2) {
            const int n = m_eri_active.n();
            CurcumaLogger::info(fmt::format(
                "QM: built 4-centre ERI (spherical, {} basis functions, {:.1f} MB) in {:.1f} ms ({} threads)",
                n, (double)n * n * n * n * 8.0 / 1.0e6,
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(),
                m_threads));
        }
    }
    return m_eri_active;
}

// =================================================================================
// Energy components {ET, EV, EJ, Ex, E_elec} for the HF run (Hartree)
// =================================================================================

std::vector<double> QMEngine::energyComponents() const
{
    return { m_et, m_ev, m_ej, m_ex, m_e_elec };
}