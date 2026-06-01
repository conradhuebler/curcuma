/*
 * <xTB SCF driver — Fock builder, eigensolver, population update>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Ported from test_xtb_scf_snapshot.cpp (tblite-validated kernels).
 *
 * Fock matrix assembly:
 *   F = H0 - 0.5 * S * diag(v_ao)  (diagonal expansion of shell potential)
 *
 * For GFN2, multipole contributions are added:
 *   vdp, vqp → F += -0.5 * (vdp contributions via dp_int + vqp contributions via qp_int)
 *
 * Generalized eigenproblem: F C = ε S C  (closed-shell)
 * Mulliken populations: P = 2 * C_occ * C_occ^T
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "native_eigensolver.h"
#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>

#include <chrono>
#include <cmath>
#include <vector>

// When MKL/BLAS is linked (EIGEN_USE_BLAS is set by CMake via USE_BLAS/USE_MKL,
// see curcuma_eigen_config.h) we call LAPACK's divide-and-conquer symmetric
// eigensolver dsyevd directly. It is several times faster than Eigen's built-in
// QR path for full eigenvector computation, which dominates the SCF on large
// systems. No LAPACKE header is needed — we declare the Fortran symbol. Without
// BLAS, solveEigenStandard() falls back to Eigen's SelfAdjointEigenSolver.
#if defined(EIGEN_USE_BLAS) || defined(USE_BLAS) || defined(USE_MKL)
#  define CURCUMA_XTB_HAVE_LAPACK_SYEVD 1
extern "C" void dsyevd_(const char* jobz, const char* uplo, const int* n,
                        double* a, const int* lda, double* w,
                        double* work, const int* lwork,
                        int* iwork, const int* liwork, int* info);
// Cholesky-based reduction of the generalized problem F C = S C ε. The overlap
// S is constant over the SCF, so we factor it once (Eigen LLT → L, BLAS-backed)
// and per iteration reduce F to standard form with LAPACK dsygst (triangular,
// ~3x fewer flops than the former dense S^{-1/2} double-GEMM), then back-transform
// eigenvectors with one triangular solve (Eigen's triangularView, BLAS dtrsm).
// dsygst is LAPACK and not declared by the BLAS backend; dtrsm/dpotrf WOULD
// collide with the BLAS-backend prototypes, so we route those through Eigen.
extern "C" void dsygst_(const int* itype, const char* uplo, const int* n,
                        double* a, const int* lda, const double* b,
                        const int* ldb, int* info);
#endif

// The whole native xTB Calculation() runs MKL-serial via curcuma::xtb::
// MklSerialScope (declared in xtb_native.h) — MKL's per-call thread team is
// pure overhead for the serial SCF up to at least 231 atoms (measured). The
// canonical C symbol MKL_Set_Num_Threads_Local lives there; nothing extra is
// needed here.

namespace curcuma::xtb {

namespace {
// Solve the standard symmetric eigenproblem A·C = C·ε for a symmetric matrix A.
// On success fills `evals` (ascending) and `evecs` (columns = eigenvectors).
// Uses LAPACK dsyevd (divide-and-conquer) when available, else Eigen.
// `A` is taken by value as a column-major buffer LAPACK can overwrite.
[[maybe_unused]] bool solveStandardSymmetric(Eigen::MatrixXd A,
                            Eigen::VectorXd& evals,
                            Eigen::MatrixXd& evecs)
{
    const int n = static_cast<int>(A.rows());
    if (n == 0) return false;
    evals.resize(n);
#ifdef CURCUMA_XTB_HAVE_LAPACK_SYEVD
    // A is Eigen::MatrixXd → column-major contiguous, exactly what LAPACK wants.
    // For a symmetric matrix the triangle choice is immaterial (full matrix set).
    const char jobz = 'V', uplo = 'U';
    int info = 0;
    int lwork = 1 + 6 * n + 2 * n * n;
    int liwork = 3 + 5 * n;
    std::vector<double> work(static_cast<size_t>(lwork));
    std::vector<int>    iwork(static_cast<size_t>(liwork));
    dsyevd_(&jobz, &uplo, &n, A.data(), &n, evals.data(),
            work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info != 0) return false;
    evecs = std::move(A);  // dsyevd overwrote A with eigenvectors (column-major)
    return true;
#else
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    if (es.info() != Eigen::Success) return false;
    evals = es.eigenvalues();
    evecs = es.eigenvectors();
    return true;
#endif
}
} // namespace

/* ------------------------------------------------------------------ *
 *  Expand shell + atom potential to AO-resolved v_ao:                *
 *    v_ao(mu) = v_sh[shell(mu)] + v_at[atom(mu)]                     *
 * ------------------------------------------------------------------ */
static void expand_potential(const BasisMap& basis,
                             const Potential& pot,
                             Eigen::VectorXd& v_ao)
{
    v_ao.resize(basis.nao);
    for (int ao = 0; ao < basis.nao; ++ao) {
        const int sh = basis.ao2sh[ao];
        const int at = basis.ao2at[ao];
        v_ao(ao) = pot.v_sh(sh) + pot.v_at(at);
    }
}

/* ------------------------------------------------------------------ *
 *  Build Fock matrix from H0, S, and the current potential.          *
 *                                                                    *
 *  Isotropic part (tblite scf/potential.f90):                        *
 *    F_μν = H0_μν - 0.5 * S_μν * (v_ao(μ) + v_ao(ν))               *
 *                                                                    *
 *  GFN2 multipole (tblite add_vmp_to_h1):                            *
 *    F_μν -= 0.5 * [dp_int[:,μ,ν]·vdp(atom_ν) + dp_int[:,ν,μ]·vdp(atom_μ)
 *                   + qp_int[:,μ,ν]·vqp(atom_ν) + qp_int[:,ν,μ]·vqp(atom_μ)]
 * ------------------------------------------------------------------ */
Matrix XTB::buildFock(const Matrix& H0,
                       const Matrix& S,
                       const Potential& pot) const
{
    const int nao = m_basis.nao;
    Eigen::VectorXd v_ao;
    expand_potential(m_basis, pot, v_ao);

    Matrix F = H0;
    const bool gfn2_mp = (m_method == MethodType::GFN2 && m_mp_initialized);

    // Parallel over the row AO mu (Claude Generated): each mu writes only row mu of
    // F, so the stripes are independent (disjoint rows, bit-identical to serial).
    // Combines the isotropic shift and the GFN2 multipole Fock contribution
    // (tblite add_vmp_to_h1) into one pass to dispatch the pool only once per build.
    const int n_threads = effectiveIntraThreads(nao);
    parallelStripes(n_threads, [&](int tid, int nth) {
    for (int mu = tid; mu < nao; mu += nth) {
        const int iat = gfn2_mp ? m_basis.ao2at[mu] : 0;
        for (int nu = 0; nu < nao; ++nu) {
            double f = -0.5 * S(mu, nu) * (v_ao(mu) + v_ao(nu));
            if (gfn2_mp) {
                const int jat = m_basis.ao2at[nu];
                double dd = 0.0;
                for (int k = 0; k < 3; ++k)
                    dd += m_dp_int[k](mu, nu) * pot.v_dp(k, jat)
                        + m_dp_int[k](nu, mu) * pot.v_dp(k, iat);
                double qq = 0.0;
                for (int k = 0; k < 6; ++k)
                    qq += m_qp_int[k](mu, nu) * pot.v_qp(k, jat)
                        + m_qp_int[k](nu, mu) * pot.v_qp(k, iat);
                f -= 0.5 * (dd + qq);
            }
            F(mu, nu) += f;
        }
    }
    });  // parallelStripes over mu
    return F;
}

/* ------------------------------------------------------------------ *
 *  Solve generalized eigenvalue problem: F C = ε S C                 *
 *  using Eigen's GeneralizedSelfAdjointEigenSolver.                  *
 *                                                                    *
 *  Two occupation schemes controlled by m_electronic_temp:           *
 *    temp == 0  → integer closed-shell (2/0 per orbital)             *
 *    temp > 0   → Fermi-Dirac smearing (fractional occupation)       *
 *                                                                    *
 *  Default: 300 K (matches TBLite default).                          *
 * ------------------------------------------------------------------ */
void XTB::buildOrthonormalizer()
{
    const int nao = m_basis.nao;
    if (nao == 0 || m_S.rows() != nao) { m_X.resize(0, 0); return; }

#ifdef CURCUMA_XTB_HAVE_LAPACK_SYEVD
    // Cache the lower Cholesky factor L of the constant overlap, S = L·Lᵀ
    // (Eigen LLT, BLAS-backed, ~nao³/3 flops — far cheaper than the former full
    // eigendecomposition of S). solveEigen() reduces each F to standard form with
    // dsygst and back-transforms eigenvectors with a triangular solve; m_X holds L.
    Eigen::LLT<Matrix> llt(m_S);
    if (llt.info() != Eigen::Success) {
        // S not positive definite (near-linear-dependent basis): fall back to the
        // robust generalized solver by leaving m_X empty.
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::warn("Overlap Cholesky failed; using generalized eigensolver");
        m_X.resize(0, 0);
        return;
    }
    m_X = llt.matrixL();             // dense lower-triangular L (upper part zero)
#else
    // No LAPACK: X = S^{-1/2} (Löwdin) via one eigendecomposition of the overlap.
    Eigen::SelfAdjointEigenSolver<Matrix> es(m_S);
    if (es.info() != Eigen::Success) { m_X.resize(0, 0); return; }
    const double smin = es.eigenvalues().minCoeff();
    if (smin < 1.0e-8) { m_X.resize(0, 0); return; }
    const Eigen::VectorXd inv_sqrt = es.eigenvalues().cwiseSqrt().cwiseInverse();
    const Matrix Us = es.eigenvectors() * inv_sqrt.asDiagonal();
    m_X.noalias() = Us * es.eigenvectors().transpose();
#endif
}

bool XTB::solveEigen(const Matrix& F, const Matrix& S)
{
    const int nao = m_basis.nao;

    // Eigenvalues + eigenvectors of F in the S metric. With the cached
    // orthonormalizer m_X = S^{-1/2} we reduce the generalized problem
    //   F C = S C ε   to the standard one   (X F X) C~ = C~ ε,   C = X C~,
    // which avoids the per-iteration Cholesky/transform of the constant S and
    // lets the fast divide-and-conquer dsyevd do the heavy lifting. The result
    // satisfies Cᵀ S C = I, identical to the generalized solver, so the gradient
    // and CPSCF paths are unaffected.
    if (m_X.rows() == nao && m_X.cols() == nao) {
        // Sub-phase timing for the verbosity>=3 SCF profile (a few steady_clock
        // reads per iteration — negligible vs the ms-scale GEMM/eigensolve work).
        using clk = std::chrono::steady_clock;
        auto ms = [](clk::time_point a, clk::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a).count();
        };
#ifdef CURCUMA_XTB_HAVE_LAPACK_SYEVD
        // m_X holds the lower Cholesky factor L of S (S = L·Lᵀ). Reduce F to the
        // standard problem  Ftil = L⁻¹·F·L⁻ᵀ  (dsygst, itype=1, lower), solve it
        // (dsyevd), then back-transform the eigenvectors  C = L⁻ᵀ·C~  by solving
        // Lᵀ·C = C~ (dtrsm). All three steps use the triangular L, so this costs
        // ~2·nao³ versus the ~6·nao³ of the former dense S^{-1/2} double-GEMM.
        const int n = nao;
        const char uplo = 'L';
        int info = 0;
        const bool use_native = (m_eigensolver == "native" || m_eigensolver == "dnc");
        const auto te0 = clk::now();
        Eigen::MatrixXd A(n, n);                 // holds the standard-form matrix, then its eigenvectors
        Eigen::VectorXd eps(n);
        if (!use_native) {
            // MKL path: reduce F to standard form (dsygst, in-place lower triangle) then
            // solve with the divide-and-conquer dsyevd.
            A = F;
            int itype = 1;
            dsygst_(&itype, &uplo, &n, A.data(), &n, m_X.data(), &n, &info);
            if (info != 0) return false;
            const auto te1 = clk::now();
            int lwork = 1 + 6 * n + 2 * n * n, liwork = 3 + 5 * n;
            std::vector<double> work(static_cast<size_t>(lwork));
            std::vector<int>    iwork(static_cast<size_t>(liwork));
            const char jobz = 'V';
            dsyevd_(&jobz, &uplo, &n, A.data(), &n, eps.data(),
                    work.data(), &lwork, iwork.data(), &liwork, &info);
            if (info != 0) return false;
            const auto te2 = clk::now();
            m_t_xfx  += ms(te0, te1);            // reduce (dsygst)
            m_t_diag += ms(te1, te2);            // dsyevd
        } else {
            // Native path (no LAPACK eigensolve): reduce A = L⁻¹·F·L⁻ᵀ with Eigen
            // triangular solves (BLAS dtrsm), then diagonalise with our own self-contained
            // solver. Y = L⁻¹·F; A = L⁻¹·Yᵀ = L⁻¹·F·L⁻ᵀ (F symmetric). Claude Generated.
            const Eigen::MatrixXd Y = m_X.triangularView<Eigen::Lower>().solve(F);
            A = m_X.triangularView<Eigen::Lower>().solve(Y.transpose());
            const auto te1 = clk::now();
            Eigen::MatrixXd Cstd;
            // Parallelise the D&C's independent subtrees on the gated intra-thread budget
            // (serial under molecule-level parallelism / for small bases).
            if (!curcuma::eigsolver::solveSymmetric(A, eps, Cstd, effectiveIntraThreads(n)))
                return false;
            A = std::move(Cstd);                 // standard-problem eigenvectors
            const auto te2 = clk::now();
            m_t_xfx  += ms(te0, te1);            // reduce (triangular solves)
            m_t_diag += ms(te1, te2);            // native eigensolve
        }
        // Back-transform C = L⁻ᵀ·C~  ⇔ solve Lᵀ·C = C~ (Eigen triangular solve → BLAS dtrsm).
        const auto te2b = clk::now();
        m_X.triangularView<Eigen::Lower>().transpose().solveInPlace(A);
        const auto te3 = clk::now();
        m_wfn.eps = std::move(eps);
        m_wfn.C   = std::move(A);               // A now holds the generalized eigenvectors
        m_t_back += ms(te2b, te3);              // back-transform (dtrsm)
#else
        const auto te0 = clk::now();
        Eigen::MatrixXd Ftil = m_X * F * m_X;   // X = S^{-1/2}: Xᵀ F X (dense GEMMs)
        const auto te1 = clk::now();
        Eigen::VectorXd eps;
        Eigen::MatrixXd Ctil;
        if (!solveStandardSymmetric(std::move(Ftil), eps, Ctil)) return false;
        const auto te2 = clk::now();
        m_wfn.eps = eps;
        m_wfn.C.noalias() = m_X * Ctil;         // back-transform C = X C~
        const auto te3 = clk::now();
        m_t_xfx += ms(te0, te1);
        m_t_diag += ms(te1, te2);
        m_t_back += ms(te2, te3);
#endif
    } else {
        // Fallback: generalized solver (m_X unavailable / near-singular overlap).
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, S);
        if (solver.info() != Eigen::Success) return false;
        m_wfn.eps = solver.eigenvalues();
        m_wfn.C   = solver.eigenvectors();
    }

    const auto t_dens0 = std::chrono::steady_clock::now();
    if (m_electronic_temp > 0.0) {
        // Fermi-Dirac smearing: bisect for Fermi level, build fractional-occupation density
        const double kT = m_electronic_temp * 3.166808e-6;  // K → Hartree
        const double n_elec = m_wfn.nocc;

        double mu_lo = m_wfn.eps.minCoeff() - 1.0;
        double mu_hi = m_wfn.eps.maxCoeff() + 1.0;
        for (int bisect = 0; bisect < 100; ++bisect) {
            const double mu = 0.5 * (mu_lo + mu_hi);
            double n_sum = 0.0;
            for (int i = 0; i < nao; ++i) {
                const double x = (m_wfn.eps(i) - mu) / kT;
                n_sum += 2.0 / (1.0 + std::exp(std::min(x, 500.0)));
            }
            if (n_sum > n_elec) mu_hi = mu;
            else                mu_lo = mu;
            if (mu_hi - mu_lo < 1e-14) break;
        }
        const double mu_f = 0.5 * (mu_lo + mu_hi);

        // Fractional occupations, then P = C·diag(occ)·Cᵀ as a single weighted GEMM
        // (replaces the per-orbital rank-1 accumulation). Eigenvalues are ascending,
        // so occupations are descending: only the leading `ncol` columns carry a
        // non-negligible weight. Restricting the GEMM to leftCols(ncol) cuts it from
        // nao³ to nao²·ncol (ncol ≈ nocc); the dropped columns weigh < 1e-12.
        Eigen::VectorXd occ(nao);
        int ncol = 0;
        for (int i = 0; i < nao; ++i) {
            const double x = (m_wfn.eps(i) - mu_f) / kT;
            occ(i) = 2.0 / (1.0 + std::exp(std::min(x, 500.0)));
            if (occ(i) > 1.0e-12) ncol = i + 1;
        }
        const auto Cocc = m_wfn.C.leftCols(ncol);
        const Matrix Cw = Cocc * occ.head(ncol).asDiagonal();
        m_wfn.P.noalias() = Cw * Cocc.transpose();
    } else {
        // Integer occupation: closed-shell, 2 electrons per occupied orbital
        const int nocc = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
        if (nocc < 0 || nocc > nao) return false;
        m_wfn.P.noalias() =
            2.0 * m_wfn.C.leftCols(nocc) * m_wfn.C.leftCols(nocc).transpose();
    }
    m_t_dens += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t_dens0).count();

    return true;
}

/* ------------------------------------------------------------------ *
 *  Saunders–Hillier level shift (Claude Generated).                  *
 *                                                                    *
 *  Raises the virtual-orbital energies by `shift` to damp the        *
 *  density's response to changes in the Fock matrix, which keeps a   *
 *  charge-sloshing iteration inside the convergence basin:           *
 *                                                                    *
 *     F' = F + shift·(S − ½·S·P·S)                                   *
 *                                                                    *
 *  The bracket is the virtual-space projector in the S metric:       *
 *  with the closed-shell P = 2·C_occ·C_occᵀ one has ½·S·P·S =        *
 *  S·C_occ·C_occᵀ·S, and S − that = S·C_virt·C_virtᵀ·S. Hence F'     *
 *  leaves the occupied orbitals (and the density) invariant at       *
 *  self-consistency and only shifts the virtuals up by `shift`.      *
 *  Reference: V. R. Saunders, I. H. Hillier, Int. J. Quantum Chem.   *
 *  7, 699 (1973).                                                    *
 * ------------------------------------------------------------------ */
Matrix XTB::applyLevelShift(const Matrix& F, const Matrix& S,
                            const Matrix& P, double shift) const
{
    if (shift == 0.0 || P.rows() != S.rows() || P.cols() != S.cols())
        return F;
    // S·P·S is symmetric; 0.5 because the closed-shell P already carries the
    // factor 2 (½·P = C_occ·C_occᵀ).
    return F + shift * (S - 0.5 * (S * P * S));
}

/* ------------------------------------------------------------------ *
 *  Update Mulliken populations from density and overlap.             *
 *    n_sh(s) = sum_{μ∈s, ν} P_μν * S_νμ                             *
 *    q_sh(s) = n0_sh(s) - n_sh(s)                                    *
 *    n_at(i) = sum_{s∈i} n_sh(s)                                     *
 * ------------------------------------------------------------------ */
void XTB::updatePopulations(const Matrix& S)
{
    const int nao = m_basis.nao;
    const int nsh = m_basis.nsh;

    // Shell populations via Mulliken: n_sh = sum_{μ∈sh, ν} P_μν * S_νμ
    // (Serial: ~1.3 ms/it — too small to amortise a per-iteration pool dispatch.)
    Vector n_sh = Vector::Zero(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int ao_start = m_basis.iao_sh[s];
        const int ao_nao   = m_basis.nao_sh[s];
        for (int mu = ao_start; mu < ao_start + ao_nao; ++mu) {
            for (int nu = 0; nu < nao; ++nu) {
                n_sh(s) += m_wfn.P(mu, nu) * S(nu, mu);
            }
        }
    }

    // Compute shell charges relative to reference
    m_wfn.q_sh = m_wfn.n0_sh - n_sh;

    // Accumulate atomic populations and charges
    m_wfn.n_at.setZero(m_atomcount);
    m_wfn.q_at.setZero(m_atomcount);
    for (int s = 0; s < nsh; ++s) {
        const int iat = m_basis.sh2at[s];
        m_wfn.n_at(iat) += n_sh(s);
    }
    for (int i = 0; i < m_atomcount; ++i) {
        m_wfn.q_at(i) = static_cast<double>(m_atoms[i]) - m_wfn.n_at(i);
        // valence-only charge (GFN convention): Z_valence - n_at
        // Actually, GFN uses q = n0_at - n_at (where n0 includes only valence)
        m_wfn.q_at(i) = m_wfn.n0_at(i) - m_wfn.n_at(i);
    }

    // GFN2 atomic multipoles via Mulliken on tblite-convention integrals
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        m_wfn.dp_at.setZero(3, m_atomcount);
        m_wfn.qp_at.setZero(6, m_atomcount);
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            double d0 = 0, d1 = 0, d2 = 0;
            double q0 = 0, q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0;
            for (int nu = 0; nu < nao; ++nu) {
                const double Pji = m_wfn.P(nu, mu);
                d0 += Pji * m_dp_int[0](nu, mu);
                d1 += Pji * m_dp_int[1](nu, mu);
                d2 += Pji * m_dp_int[2](nu, mu);
                q0 += Pji * m_qp_int[0](nu, mu);
                q1 += Pji * m_qp_int[1](nu, mu);
                q2 += Pji * m_qp_int[2](nu, mu);
                q3 += Pji * m_qp_int[3](nu, mu);
                q4 += Pji * m_qp_int[4](nu, mu);
                q5 += Pji * m_qp_int[5](nu, mu);
            }
            // Sign: multipoles are valence deviations (like q_sh)
            m_wfn.dp_at(0, iat) -= d0;
            m_wfn.dp_at(1, iat) -= d1;
            m_wfn.dp_at(2, iat) -= d2;
            m_wfn.qp_at(0, iat) -= q0;
            m_wfn.qp_at(1, iat) -= q1;
            m_wfn.qp_at(2, iat) -= q2;
            m_wfn.qp_at(3, iat) -= q3;
            m_wfn.qp_at(4, iat) -= q4;
            m_wfn.qp_at(5, iat) -= q5;
        }
    }
}

/* ------------------------------------------------------------------ *
 *  SCF convergence check.                                            *
 *  Returns true when max |Δq| < threshold and |ΔE| < threshold.      *
 * ------------------------------------------------------------------ */
static bool checkConvergence(const Vector& q_sh_old,
                             const Vector& q_sh_new,
                             double e_old, double e_new,
                             double threshold)
{
    const double dq = (q_sh_new - q_sh_old).cwiseAbs().maxCoeff();
    const double de = std::fabs(e_new - e_old);
    return (dq < threshold && de < threshold);
}

/* ------------------------------------------------------------------ *
 *  buildFockFromPotential()                                          *
 *                                                                    *
 *  H0-free delta-Fock from a potential perturbation δpot:           *
 *    δF_μν = -0.5·S_μν·(δv_ao(μ) + δv_ao(ν))                        *
 *  + GFN2 multipole terms if δpot.v_dp/v_qp are populated.          *
 *                                                                    *
 *  Used in the CPSCF orbital-Hessian to compute (K·z)_ai:           *
 *    δF = buildFockFromPotential(δpot)  → (Cᵀ δF C)_ai             *
 * ------------------------------------------------------------------ */
Matrix XTB::buildFockFromPotential(const Potential& dpot) const
{
    const int nao = m_basis.nao;
    const int nat = m_atomcount;
    Eigen::VectorXd v_ao;
    expand_potential(m_basis, dpot, v_ao);

    Matrix dF = Matrix::Zero(nao, nao);
    for (int mu = 0; mu < nao; ++mu)
        for (int nu = 0; nu < nao; ++nu)
            dF(mu, nu) -= 0.5 * m_S(mu, nu) * (v_ao(mu) + v_ao(nu));

    if (m_method == MethodType::GFN2 && m_mp_initialized
        && dpot.v_dp.cols() == nat && dpot.v_qp.cols() == nat) {
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            for (int nu = 0; nu < nao; ++nu) {
                const int jat = m_basis.ao2at[nu];
                double dd = 0.0;
                for (int k = 0; k < 3; ++k)
                    dd += m_dp_int[k](mu, nu) * dpot.v_dp(k, jat)
                        + m_dp_int[k](nu, mu) * dpot.v_dp(k, iat);
                double qq = 0.0;
                for (int k = 0; k < 6; ++k)
                    qq += m_qp_int[k](mu, nu) * dpot.v_qp(k, jat)
                        + m_qp_int[k](nu, mu) * dpot.v_qp(k, iat);
                dF(mu, nu) -= 0.5 * (dd + qq);
            }
        }
    }
    return dF;
}

/* ------------------------------------------------------------------ *
 *  mullikenFromDensity()                                             *
 *                                                                    *
 *  Compute Mulliken charges and multipoles from a density            *
 *  perturbation δP without modifying any member state.              *
 *                                                                    *
 *    δn_sh(s) = Σ_{μ∈s,ν} δP_μν · S_νμ                              *
 *    δq_sh    = -δn_sh     (n0 constant)                             *
 *    δq_at    = -Σ_{s∈i} δn_sh(s)                                   *
 *    δdp_at, δqp_at  via tblite-convention integrals (GFN2 only)    *
 * ------------------------------------------------------------------ */
void XTB::mullikenFromDensity(const Matrix& dP,
                               Vector& dq_sh, Vector& dq_at,
                               Eigen::MatrixXd& ddp_at,
                               Eigen::MatrixXd& dqp_at) const
{
    const int nao = m_basis.nao;
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;

    Vector dn_sh = Vector::Zero(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int ao_start = m_basis.iao_sh[s];
        const int ao_nao   = m_basis.nao_sh[s];
        for (int mu = ao_start; mu < ao_start + ao_nao; ++mu)
            for (int nu = 0; nu < nao; ++nu)
                dn_sh(s) += dP(mu, nu) * m_S(nu, mu);
    }
    dq_sh = -dn_sh;

    dq_at = Vector::Zero(nat);
    for (int s = 0; s < nsh; ++s)
        dq_at(m_basis.sh2at[s]) -= dn_sh(s);

    ddp_at = Eigen::MatrixXd::Zero(3, nat);
    dqp_at = Eigen::MatrixXd::Zero(6, nat);
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            for (int nu = 0; nu < nao; ++nu) {
                const double dPji = dP(nu, mu);
                for (int k = 0; k < 3; ++k)
                    ddp_at(k, iat) -= dPji * m_dp_int[k](nu, mu);
                for (int k = 0; k < 6; ++k)
                    dqp_at(k, iat) -= dPji * m_qp_int[k](nu, mu);
            }
        }
    }
}

} // namespace curcuma::xtb
