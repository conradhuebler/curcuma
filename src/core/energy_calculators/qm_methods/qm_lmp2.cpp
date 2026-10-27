/*
 * <Local PNO-based RI-MP2 (closed shell, energy only), WP15 -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Method, equations and literature: qm_lmp2.h. Python counterpart of every step:
 * scripts/local_mp2/lmp2_common.py (functions named in the comments below).
 *
 * Layout of this file (top to bottom = order of QMEngine::computeLMP2):
 *   1. dipole integrals <mu|r|nu> (Boys localisation, orbital centroids, dipole estimate)
 *   2. localisation: Cholesky start orbitals, Jacobi sweeps for Pipek-Mezey / Boys
 *   3. Boughton-Pulay domains, PAOs, pseudo-canonical pair spaces
 *   4. QMEngine::computeLMP2: half-transformed RI factors B_{i,r}^P, PNOs, distant pairs, the LMP2 iterations
 *
 * Memory: B_{i,r}^P only for the PAOs r on the atoms of the extended domain of i (the union of [j] over all
 * kept partners j); per kept pair its PNO coefficients, K and amplitudes; no o^2 v^2 tensor.
 *
 * Claude Generated (Oct 2026, WP15). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_lmp2.h"
#include "src/core/gpu_fallback.h"

#include "qm_engine.h"
#include "qm_mp2.h"  // frozenCoreOrbitals
#include "qm_ri.h"

#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>

qmint::LMP2Options::LMP2Options() = default;
qmint::LMP2Options::LMP2Options(const LMP2Options&) = default;
qmint::LMP2Options::LMP2Options(LMP2Options&&) noexcept = default;
qmint::LMP2Options& qmint::LMP2Options::operator=(const LMP2Options&) = default;
qmint::LMP2Options& qmint::LMP2Options::operator=(LMP2Options&&) noexcept = default;
qmint::LMP2Options::~LMP2Options() = default;
qmint::LMP2Result::LMP2Result() = default;
qmint::LMP2Result::LMP2Result(const LMP2Result&) = default;
qmint::LMP2Result::LMP2Result(LMP2Result&&) noexcept = default;
qmint::LMP2Result& qmint::LMP2Result::operator=(const LMP2Result&) = default;
qmint::LMP2Result& qmint::LMP2Result::operator=(LMP2Result&&) noexcept = default;
qmint::LMP2Result::~LMP2Result() = default;

namespace {

using MatX = Eigen::MatrixXd;  // column-major, local work matrices
using VecX = Eigen::VectorXd;

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

std::string auxLabel(const std::string& path)
{
    std::string s = path;
    const size_t slash = s.find_last_of('/');
    if (slash != std::string::npos) s = s.substr(slash + 1);
    const size_t dot = s.rfind(".turbomole");
    if (dot != std::string::npos) s = s.substr(0, dot);
    return s;
}

// =====================================================================================================
// 1. Dipole integrals
// =====================================================================================================

/// 1D Obara-Saika overlap table s[a][b] (without the Gaussian prefactor) for powers a <= la, b <= lb:
///   s[a+1][b] = PA s[a][b] + (a s[a-1][b] + b s[a][b-1]) / (2g),  s[a][b+1] = PB s[a][b] + (...) / (2g).
/// Claude Generated (Oct 2026).
void overlapTable1D(int la, int lb, double PA, double PB, double g, double s[6][7])
{
    const double h = 0.5 / g;
    for (int a = 0; a <= la; ++a)
        for (int b = 0; b <= lb; ++b) {
            double v;
            if (a == 0 && b == 0)
                v = 1.0;
            else if (b == 0)
                v = PA * s[a - 1][0] + (a > 1 ? (a - 1) * h * s[a - 2][0] : 0.0);
            else
                v = PB * s[a][b - 1] + (a > 0 ? a * h * s[a - 1][b - 1] : 0.0) + (b > 1 ? (b - 1) * h * s[a][b - 2] : 0.0);
            s[a][b] = v;
        }
}

/// @brief Overlap S and dipole matrices <mu|x|nu>, <mu|y|nu>, <mu|z|nu> (origin 0, Bohr) over the flat cartesian
/// basis (coefficients include the primitive normalisation, as in qm_integrals.cpp). Uses
/// <a|x|b> = <a|b + 1_x> + B_x <a|b> (the power of b raised by one about its own centre).
/// S is returned for the self-check against the engine's overlap. Claude Generated (Oct 2026).
void cartesianDipoleIntegrals(const std::vector<GTO::Orbital>& basis, Matrix& S, Matrix& X, Matrix& Y, Matrix& Z)
{
    const int n = (int)basis.size();
    S = Matrix::Zero(n, n);
    X = Matrix::Zero(n, n);
    Y = Matrix::Zero(n, n);
    Z = Matrix::Zero(n, n);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        const GTO::Orbital& a = basis[i];
        int l1, m1, n1;
        GTO::orbitalTypeToComponents(a.type, l1, m1, n1);
        for (int j = i; j < n; ++j) {
            const GTO::Orbital& b = basis[j];
            int l2, m2, n2;
            GTO::orbitalTypeToComponents(b.type, l2, m2, n2);
            double ss = 0, sx = 0, sy = 0, sz = 0;
            for (size_t ia = 0; ia < a.exponents.size(); ++ia)
                for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
                    const double al = a.exponents[ia], be = b.exponents[ib], g = al + be;
                    const double Px = (al * a.x + be * b.x) / g, Py = (al * a.y + be * b.y) / g, Pz = (al * a.z + be * b.z) / g;
                    const double R2 = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
                    const double K = std::exp(-al * be / g * R2);
                    if (K < 1e-300) continue;
                    const double pref = a.coefficients[ia] * b.coefficients[ib] * std::pow(M_PI / g, 1.5) * K;
                    double tx[6][7], ty[6][7], tz[6][7];
                    overlapTable1D(l1, l2 + 1, Px - a.x, Px - b.x, g, tx);
                    overlapTable1D(m1, m2 + 1, Py - a.y, Py - b.y, g, ty);
                    overlapTable1D(n1, n2 + 1, Pz - a.z, Pz - b.z, g, tz);
                    const double ox = tx[l1][l2], oy = ty[m1][m2], oz = tz[n1][n2];
                    ss += pref * ox * oy * oz;
                    sx += pref * (tx[l1][l2 + 1] + b.x * ox) * oy * oz;
                    sy += pref * ox * (ty[m1][m2 + 1] + b.y * oy) * oz;
                    sz += pref * ox * oy * (tz[n1][n2 + 1] + b.z * oz);
                }
            S(i, j) = S(j, i) = ss;
            X(i, j) = X(j, i) = sx;
            Y(i, j) = Y(j, i) = sy;
            Z(i, j) = Z(j, i) = sz;
        }
    }
}

MatX toEigen(const Matrix& M)
{
    return MatX(M);
}

// =====================================================================================================
// 2. Localisation
// =====================================================================================================

/// @brief Cholesky orbitals as start of the localisation (Aquilante, Pedersen, Koch, J. Chem. Phys. 125, 174101
/// (2006)): pivoted Cholesky of the AO density D = C C^T (rank o). The factor columns span the occupied space and
/// are S-orthonormal (L = C U with U orthogonal); a Loewdin step removes the rounding. Claude Generated (Oct 2026).
bool choleskyOrbitals(const MatX& Cval, const MatX& S, MatX& L)
{
    const int n = (int)Cval.rows(), o = (int)Cval.cols();
    const MatX D = Cval * Cval.transpose();
    VecX diag = D.diagonal();
    std::vector<char> used(n, 0);
    L = MatX::Zero(n, o);
    for (int k = 0; k < o; ++k) {
        int p = -1;
        double best = -1.0;
        for (int r = 0; r < n; ++r)
            if (!used[r] && diag(r) > best) {
                best = diag(r);
                p = r;
            }
        if (p < 0 || best < 1e-12) return false;
        used[p] = 1;
        VecX col = D.col(p);
        if (k > 0) col -= L.leftCols(k) * L.row(p).leftCols(k).transpose();
        col /= std::sqrt(best);
        L.col(k) = col;
        diag -= col.cwiseAbs2();
    }
    const MatX M = L.transpose() * S * L;
    Eigen::SelfAdjointEigenSolver<MatX> es(M);
    if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() < 1e-8) return false;
    L = L * (es.eigenvectors() * es.eigenvalues().cwiseInverse().cwiseSqrt().asDiagonal() * es.eigenvectors().transpose());
    return true;
}

/// @brief Jacobi sweeps maximising F = sum_M sum_s (M_ss)^2 over orthogonal rotations of the orbitals.
/// PM: M = Mulliken population matrices Q^A_st = 1/2 sum_{mu in A} (C_mu,s (SC)_mu,t + C_mu,t (SC)_mu,s), one per
/// atom; Boys: M = <s|x|t>, <s|y|t>, <s|z|t>. For a pair (s,t) and s' = c s + sn t, t' = -sn s + c t:
///   dF = A (1 - cos 4g) + B sin 4g,  A = sum_M [M_st^2 - (M_ss - M_tt)^2 / 4],  B = sum_M M_st (M_ss - M_tt),
/// maximal at (cos 4g, sin 4g) = (-A, B) / sqrt(A^2 + B^2) (Pipek, Mezey 1989, eq. 22-24; Edmiston-Ruedenberg form).
/// Disjoint pairs are rotated together (round-robin schedule, o - 1 rounds per sweep), so the result does not
/// depend on the thread count. Returns the number of sweeps (negative = not converged).
/// Claude Generated (Oct 2026).
int jacobiLocalize(MatX& L, std::vector<MatX>& M, int threads, double tol, int max_sweeps, double& functional)
{
    const int o = (int)L.cols();
    const int nm = (int)M.size();
    auto F = [&]() {
        double f = 0.0;
        for (const MatX& m : M) f += m.diagonal().squaredNorm();
        return f;
    };
    functional = F();
    if (o < 2) return 0;
    const int np = o + (o % 2);  // players incl. one dummy for odd o
    std::vector<int> ring(np);
    std::iota(ring.begin(), ring.end(), 0);
    std::vector<int> ps(np / 2), pt(np / 2);
    std::vector<double> cs(np / 2), sn(np / 2), gain(np / 2);
    for (int sweep = 1; sweep <= max_sweeps; ++sweep) {
        const double f0 = functional;
        for (int round = 0; round < np - 1; ++round) {
            // circle method: ring[0] fixed, pairs (ring[k], ring[np-1-k])
            for (int k = 0; k < np / 2; ++k) {
                ps[k] = std::min(ring[k], ring[np - 1 - k]);
                pt[k] = std::max(ring[k], ring[np - 1 - k]);
            }
#pragma omp parallel for schedule(static) num_threads(threads)
            for (int k = 0; k < np / 2; ++k) {
                cs[k] = 1.0;
                sn[k] = 0.0;
                gain[k] = 0.0;
                const int s = ps[k], t = pt[k];
                if (t >= o) continue;  // dummy partner
                double A = 0.0, B = 0.0;
                for (int m = 0; m < nm; ++m) {
                    const double mst = M[m](s, t), d = M[m](s, s) - M[m](t, t);
                    A += mst * mst - 0.25 * d * d;
                    B += mst * d;
                }
                const double r = std::sqrt(A * A + B * B);
                if (r < 1e-14) continue;
                const double g = 0.25 * std::atan2(B, -A);
                cs[k] = std::cos(g);
                sn[k] = std::sin(g);
                gain[k] = A + r;
            }
            // rotate every matrix (rows, then columns) and the orbital coefficients
#pragma omp parallel for schedule(dynamic) num_threads(threads)
            for (int m = 0; m < nm; ++m) {
                MatX& Mm = M[m];
                for (int k = 0; k < np / 2; ++k) {
                    if (sn[k] == 0.0) continue;
                    const int s = ps[k], t = pt[k];
                    const double c = cs[k], si = sn[k];
                    for (int q = 0; q < o; ++q) {
                        const double a = Mm(s, q), b = Mm(t, q);
                        Mm(s, q) = c * a + si * b;
                        Mm(t, q) = -si * a + c * b;
                    }
                }
                for (int k = 0; k < np / 2; ++k) {
                    if (sn[k] == 0.0) continue;
                    const int s = ps[k], t = pt[k];
                    const double c = cs[k], si = sn[k];
                    for (int q = 0; q < o; ++q) {
                        const double a = Mm(q, s), b = Mm(q, t);
                        Mm(q, s) = c * a + si * b;
                        Mm(q, t) = -si * a + c * b;
                    }
                }
            }
#pragma omp parallel for schedule(static) num_threads(threads)
            for (int k = 0; k < np / 2; ++k) {
                if (sn[k] == 0.0) continue;
                const int s = ps[k], t = pt[k];
                const VecX a = L.col(s), b = L.col(t);
                L.col(s) = cs[k] * a + sn[k] * b;
                L.col(t) = -sn[k] * a + cs[k] * b;
            }
            // advance the ring (all but element 0 move one position)
            const int last = ring[np - 1];
            for (int k = np - 1; k > 1; --k) ring[k] = ring[k - 1];
            ring[1] = last;
        }
        functional = F();
        if (std::abs(functional - f0) < tol * std::max(1.0, std::abs(functional))) return sweep;
    }
    return -max_sweeps;
}

constexpr double kLocTol = 1e-11;   // relative change of the functional per sweep (both routes)
constexpr int kLocMaxSweeps = 500;

/// @brief C = A B with the columns of B split over OpenMP threads; each block is one sequential BLAS call (the
/// linked MKL is sequential, see qmint::gemmABt in qm_mp2.h). Claude Generated (Oct 2026).
MatX productParallel(const MatX& A, const MatX& B, int threads)
{
    MatX C(A.rows(), B.cols());
    const int nc = (int)B.cols();
    const int nb = std::max(1, std::min(threads, nc));
    const int chunk = (nc + nb - 1) / std::max(1, nb);
#pragma omp parallel for schedule(static) num_threads(nb)
    for (int b = 0; b < nb; ++b) {
        const int c0 = b * chunk, w = std::min(chunk, nc - c0);
        if (w > 0) C.middleCols(c0, w).noalias() = A * B.middleCols(c0, w);
    }
    return C;
}

/// @brief C = A^T B, columns of B split over threads (as productParallel). Claude Generated (Oct 2026).
MatX productParallelAtB(const MatX& A, const MatX& B, int threads)
{
    MatX C(A.cols(), B.cols());
    const int nc = (int)B.cols();
    const int nb = std::max(1, std::min(threads, nc));
    const int chunk = (nc + nb - 1) / std::max(1, nb);
#pragma omp parallel for schedule(static) num_threads(nb)
    for (int b = 0; b < nb; ++b) {
        const int c0 = b * chunk, w = std::min(chunk, nc - c0);
        if (w > 0) C.middleCols(c0, w).noalias() = A.transpose() * B.middleCols(c0, w);
    }
    return C;
}

/// The circle-method round-robin schedule of jacobiLocalize (same order), used by pmLocalizeFast: round r of a sweep pairs (ring[k],
/// ring[np-1-k]); np = o rounded up to even (one dummy player). Claude Generated (Oct 2026).
struct RoundRobin {
    int np;
    std::vector<int> ring;
    explicit RoundRobin(int o) : np(o + (o % 2)), ring(np) { std::iota(ring.begin(), ring.end(), 0); }
    void pairs(std::vector<int>& ps, std::vector<int>& pt) const
    {
        for (int k = 0; k < np / 2; ++k) {
            ps[k] = std::min(ring[k], ring[np - 1 - k]);
            pt[k] = std::max(ring[k], ring[np - 1 - k]);
        }
    }
    void advance()
    {
        const int last = ring[np - 1];
        for (int k = np - 1; k > 1; --k) ring[k] = ring[k - 1];
        ring[1] = last;
    }
};

/// Counters of pmLocalizeFast (over all sweeps). Claude Generated (Oct 2026).
struct PMSweepStats {
    long pairs_total = 0, pairs_evaluated = 0, pairs_rotated = 0;  ///< evaluated = Q^A_st formed
    double mean_support = 0.0;  ///< atoms per orbital support at the end
};

/// @brief Pipek-Mezey Jacobi sweeps without the natm population matrices (-qm.lmp2_localizer_algo fast).
/// Same objective F = sum_A sum_s (Q^A_ss)^2, pair formula, round-robin schedule and convergence test as
/// jacobiLocalize, but per pair (s,t) Q^A_st is formed from the columns s, t of L and SL = S L,
///   Q^A_st = 1/2 sum_{mu in A} (L_mu,s (SL)_mu,t + L_mu,t (SL)_mu,s),
/// and the populations P_As = Q^A_ss with sig_s = sum_A P_As^2 are refreshed for the two rotated orbitals only.
/// Cost per sweep O(o^2 n) unscreened instead of O(natm o^3) (rotating natm o x o matrices), parallel over the
/// disjoint pairs of a round (thread-count independent).
/// Screening (screen > 0), w_As = max(||L_{A,s}||, ||(SL)_{A,s}||):
///  (1) atom A is in the support of s if w_As > screen; a pair is evaluated over the atoms common to both supports and
///      skipped if there is none: outside them |Q^A_st| <= w_As w_At <= ~screen (Cauchy-Schwarz), and the
///      (P_As - P_At)^2 terms outside are taken exactly from sig_s, sig_t;
///  (2) before Q^A_st is formed, the rotation angle is bounded from the stored w and P alone:
///      |B| <= Bmax = sum_A w_As w_At (|P_As| + |P_At|),  A <= Aub = -(sig_s + sig_t)/4 + sum_A (w_As w_At)^2
///      + 1/2 sum_A |P_As P_At|, and |gamma| <= Bmax / (4 |Aub|) if Aub < 0; the pair is skipped if this bound is
///      <= screen (a rotation by at most screen radians).
/// The O(n) work (Q^A_st, rotation, refresh) then falls on the pairs that rotate; the support walk is O(o^2 k).
/// screen = 0 keeps every pair over all atoms (the classic arithmetic in O(o^2 n)). F is always the exact sum over
/// all atoms.
/// Returns the number of sweeps (negative = not converged). Claude Generated (Oct 2026).
int pmLocalizeFast(MatX& L, const MatX& S, const std::vector<std::vector<int>>& atom_aos, int threads, double tol,
                   int max_sweeps, double screen, double& functional, PMSweepStats& st)
{
    const int n = (int)L.rows(), o = (int)L.cols(), natm = (int)atom_aos.size();
    MatX SL = productParallel(S, L, threads);
    MatX pop(natm, o), W(natm, o);  // populations P_As and w_As = max(||L_{A,s}||, ||(SL)_{A,s}||)
    std::vector<double> sig(o, 0.0);
    std::vector<std::vector<int>> supp(o);
    const double screen2 = screen * screen;
    auto refresh = [&](int s) {
        double sg = 0.0;
        supp[s].clear();
        const double* l = L.col(s).data();
        const double* q = SL.col(s).data();
        for (int A = 0; A < natm; ++A) {
            double p = 0.0, a2 = 0.0, b2 = 0.0;
            for (int mu : atom_aos[A]) {
                p += l[mu] * q[mu];
                a2 += l[mu] * l[mu];
                b2 += q[mu] * q[mu];
            }
            pop(A, s) = p;
            W(A, s) = std::sqrt(std::max(a2, b2));
            sg += p * p;
            if (!(screen > 0.0) || std::max(a2, b2) > screen2) supp[s].push_back(A);
        }
        sig[s] = sg;
    };
    auto total = [&]() {
        double f = 0.0;
        for (int s = 0; s < o; ++s) f += sig[s];
        return f;
    };
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int s = 0; s < o; ++s) refresh(s);
    functional = total();
    auto finish = [&](int ret) {
        double sup = 0.0;
        for (int s = 0; s < o; ++s) sup += supp[s].size();
        st.mean_support = o ? sup / o : 0.0;
        return ret;
    };
    if (o < 2) return finish(0);
    RoundRobin rr(o);
    const int npair = rr.np / 2;
    std::vector<int> ps(npair), pt(npair);
    for (int sweep = 1; sweep <= max_sweeps; ++sweep) {
        const double f0 = functional;
        long n_eval = 0, n_rot = 0;
        for (int round = 0; round < rr.np - 1; ++round) {
            rr.pairs(ps, pt);
#pragma omp parallel for schedule(dynamic, 4) num_threads(threads) reduction(+ : n_eval, n_rot)
            for (int k = 0; k < npair; ++k) {
                const int s = ps[k], t = pt[k];
                if (t >= o) continue;
                const double* ls = L.col(s).data();
                const double* lt = L.col(t).data();
                const double* qs = SL.col(s).data();
                const double* qt = SL.col(t).data();
                const std::vector<int>& us = supp[s];
                const std::vector<int>& ut = supp[t];
                std::vector<int> common;  // atoms in both supports
                common.reserve(std::min(us.size(), ut.size()));
                std::set_intersection(us.begin(), us.end(), ut.begin(), ut.end(), std::back_inserter(common));
                if (common.empty()) continue;
                if (screen > 0.0) {  // (2) bound on the rotation angle from w and P
                    double bmax = 0.0, s2 = 0.0, cc = 0.0;
                    for (int At : common) {
                        const double ww = W(At, s) * W(At, t), pa = pop(At, s), pb = pop(At, t);
                        bmax += ww * (std::abs(pa) + std::abs(pb));
                        s2 += ww * ww;
                        cc += std::abs(pa * pb);
                    }
                    const double aub = -0.25 * (sig[s] + sig[t]) + s2 + 0.5 * cc;
                    if (aub < 0.0 && bmax <= 4.0 * std::abs(aub) * screen) continue;
                }
                double A = 0.0, B = 0.0, ps2 = 0.0, pt2 = 0.0;
                for (int At : common) {
                    double qst = 0.0;
                    for (int mu : atom_aos[At]) qst += ls[mu] * qt[mu] + lt[mu] * qs[mu];
                    qst *= 0.5;
                    const double pa = pop(At, s), pb = pop(At, t), d = pa - pb;
                    A += qst * qst - 0.25 * d * d;
                    B += qst * d;
                    ps2 += pa * pa;
                    pt2 += pb * pb;
                }
                ++n_eval;
                A -= 0.25 * ((sig[s] - ps2) + (sig[t] - pt2));  // (P_As - P_At)^2 on the remaining atoms
                const double r = std::sqrt(A * A + B * B);
                if (r < 1e-14) continue;
                const double g = 0.25 * std::atan2(B, -A);
                const double c = std::cos(g), si = std::sin(g);
                if (si == 0.0) continue;
                ++n_rot;
                double* Ls = L.col(s).data();
                double* Lt = L.col(t).data();
                double* Qs = SL.col(s).data();
                double* Qt = SL.col(t).data();
                for (int mu = 0; mu < n; ++mu) {
                    const double a1 = Ls[mu], b1 = Lt[mu];
                    Ls[mu] = c * a1 + si * b1;
                    Lt[mu] = -si * a1 + c * b1;
                    const double a2 = Qs[mu], b2 = Qt[mu];
                    Qs[mu] = c * a2 + si * b2;
                    Qt[mu] = -si * a2 + c * b2;
                }
                refresh(s);
                refresh(t);
            }
            rr.advance();
        }
        st.pairs_total += (long)o * (o - 1) / 2;
        st.pairs_evaluated += n_eval;
        st.pairs_rotated += n_rot;
        functional = total();
        if (std::abs(functional - f0) < tol * std::max(1.0, std::abs(functional))) return finish(sweep);
    }
    return finish(-max_sweeps);
}

// =====================================================================================================
// 3. Domains and pair spaces
// =====================================================================================================

/// @brief Boughton-Pulay domain of the orbital c (lmp2_common.bp_domain): atoms by decreasing Mulliken gross
/// population until 1 - <c|P_A|c> < thr, P_A = projector onto the AOs of the selected atoms (least-squares fit).
/// Claude Generated (Oct 2026).
std::vector<int> bpDomain(const VecX& c, const MatX& S, const std::vector<std::vector<int>>& atom_aos, double thr)
{
    const VecX Sc = S * c;
    const int natm = (int)atom_aos.size();
    std::vector<double> pop(natm, 0.0);
    for (int A = 0; A < natm; ++A)
        for (int mu : atom_aos[A]) pop[A] += c(mu) * Sc(mu);
    std::vector<int> order(natm);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return pop[a] > pop[b]; });
    std::vector<int> sel, idx;
    for (int A : order) {
        sel.push_back(A);
        idx.insert(idx.end(), atom_aos[A].begin(), atom_aos[A].end());
        const int m = (int)idx.size();
        MatX Ss(m, m);
        VecX b(m);
        for (int r = 0; r < m; ++r) {
            b(r) = Sc(idx[r]);
            for (int q = 0; q < m; ++q) Ss(r, q) = S(idx[r], idx[q]);
        }
        const double f = 1.0 - b.dot(Ss.ldlt().solve(b));
        if (f < thr) break;
    }
    std::sort(sel.begin(), sel.end());
    return sel;
}

/// Orthonormal pseudo-canonical virtual space spanned by the PAOs `idx` (lmp2_common.pair_space):
/// Xc = coefficients over these PAOs (m_ao x m), eps = pseudo-canonical energies. Claude Generated (Oct 2026).
struct PaoSpace {
    MatX Xc;
    VecX eps;
};
constexpr double kPaoLindep = 1.0e-6;  // eigenvalue cut of the PAO overlap (study: "1e-6 redundancy cut")

bool paoSpace(const std::vector<int>& idx, const MatX& Spao, const MatX& Fpao, PaoSpace& out)
{
    const int m = (int)idx.size();
    MatX Ss(m, m), Fs(m, m);
    for (int q = 0; q < m; ++q)
        for (int r = 0; r < m; ++r) {
            Ss(r, q) = Spao(idx[r], idx[q]);
            Fs(r, q) = Fpao(idx[r], idx[q]);
        }
    Eigen::SelfAdjointEigenSolver<MatX> es(Ss);
    if (es.info() != Eigen::Success) return false;
    int keep = 0;
    for (int r = 0; r < m; ++r)
        if (es.eigenvalues()(r) > kPaoLindep) ++keep;
    MatX X(m, keep);
    for (int r = 0, c = 0; r < m; ++r)
        if (es.eigenvalues()(r) > kPaoLindep) X.col(c++) = es.eigenvectors().col(r) / std::sqrt(es.eigenvalues()(r));
    const MatX Fx = X.transpose() * Fs * X;
    Eigen::SelfAdjointEigenSolver<MatX> ef(Fx);
    if (ef.info() != Eigen::Success) return false;
    out.Xc = X * ef.eigenvectors();
    out.eps = ef.eigenvalues();
    return true;
}

/// Closed-shell pair energy (2 - delta_ij) sum_ab T_ab (2 K_ab - K_ba).
double pairEnergy(const MatX& T, const MatX& K, bool diag)
{
    if (T.size() == 0) return 0.0;
    return (diag ? 1.0 : 2.0) * (T.array() * (2.0 * K - K.transpose()).array()).sum();
}

/// Everything one kept pair carries through the iterations.
struct PairData {
    int i = 0, j = 0;
    std::vector<int> idx;  ///< PAO (= AO) indices of [i] u [j]
    MatX c;                ///< PNO coefficients over the PAOs idx (m_ao x npno)
    VecX eps;              ///< pseudo-canonical PNO energies
    MatX K;                ///< (ia|jb) in the PNO basis
    MatX den;              ///< e_a + e_b - F_ii - F_jj
    double e_sc_base = 0.0, e_sc_pno = 0.0;
    int m_pao = 0;         ///< pair (PAO) space dimension before the PNO cut
};

/// One coupling term of the residual of pair p: -f S_pq X_q S_pq^T with X_q = T_q or T_q^T.
struct Coupling {
    int q = 0;
    int slot = 0;  ///< index of S_pq in the pair's overlap list
    double f = 0.0;
    bool transpose = false;
};

/// Pulay DIIS on flat amplitude vectors (error = preconditioned residual). Claude Generated (Oct 2026).
class Diis {
public:
    explicit Diis(int nmax) : m_nmax(nmax) {}
    VecX update(const VecX& x, const VecX& e)
    {
        if (m_nmax < 2) return x;
        m_x.push_back(x);
        m_e.push_back(e);
        if ((int)m_x.size() > m_nmax) {
            m_x.erase(m_x.begin());
            m_e.erase(m_e.begin());
        }
        const int n = (int)m_x.size();
        if (n < 2) return x;
        MatX B = MatX::Constant(n + 1, n + 1, -1.0);
        B(n, n) = 0.0;
        for (int a = 0; a < n; ++a)
            for (int b = a; b < n; ++b) B(a, b) = B(b, a) = m_e[a].dot(m_e[b]);
        // scale the error block to O(1): near convergence its entries are ~1e-20 next to the -1 border,
        // and the unscaled system stalls the iteration at max|R| ~ 1e-9 (measured, water dimer)
        const double scale = B.topLeftCorner(n, n).diagonal().maxCoeff();
        if (scale > 0.0) B.topLeftCorner(n, n) /= scale;
        VecX rhs = VecX::Zero(n + 1);
        rhs(n) = -1.0;
        const VecX c = B.colPivHouseholderQr().solve(rhs);
        if (!c.allFinite()) return x;
        VecX out = VecX::Zero(x.size());
        for (int a = 0; a < n; ++a) out += c(a) * m_x[a];
        return out;
    }

private:
    int m_nmax;
    std::vector<VecX> m_x, m_e;
};

// =====================================================================================================
// 3b. Local fitting with integral-direct 3-index integrals (Claude Generated, Oct 2026)
// =====================================================================================================
//
// Instead of the stored global tensor, every orbital i keeps the RAW half-transformed integrals
//   X_i(P, r) = (P | i r) = sum_mn L_mi Cpao_nr (P|mn),   P in G_i,  r in ext(i),
// where G_i = union of the fit domains F_j of all kept partners j of i (incl. i itself) and ext(i) are the PAOs of
// the extended domain (as for the global B factors). They are computed aux shell by aux shell without storing
// (P|mn): for aux shell c on atom K, every orbital i with K in G_i gets its rows of X_i from
//   T_p(m, n) = (p|mn)  ->  Y = L^T T  ->  X_i(p, :) = Y_i(p, :) Cpao(:, ext(i)),
// with m restricted to the AO support of the orbitals (|L_mi| > ao_cut) and n to the support of their PAOs
// (|Cpao_nr| > ao_cut); dropped coefficients are set to zero, so ao_cut = 0 is the exact half transformation.
// Then the pair integrals are fitted locally:
//   pair    : (ia|jb) = sum_{P,Q in F_ij} (ia|P) [J_F^-1]_PQ (Q|jb),  F_ij = F_i u F_j, J_F = metric block
//             (symmetric fit within one domain; with F = all functions this is the global fit),
//   orbital : robust fit (Dunlap) with the orbital domains, d_ia = J_Fi^-1 (P|ia)_{P in F_i}:
//             (ia|jb) = d_ia^T (P|jb)_{F_i} + (ia|Q)_{F_j} d_jb - d_ia^T J_{Fi,Fj} d_jb
//             (error second order in the fitting errors of ia and jb; F = all gives the global fit again).
// B. I. Dunlap, J. Mol. Struct. THEOCHEM 529, 37 (2000) (robust fitting); H.-J. Werner, F. R. Manby, P. J. Knowles,
// J. Chem. Phys. 118, 8149 (2003) (local fit domains in LMP2); both recalled, not looked up.

/// Factor of one fit domain: B = F X fitted (cholesky: F = L^-1 of the Coulomb-normalised metric block; eigen:
/// F = w^-1/2 U^T with eigenvalues below the cutoff dropped, as RIJK). Rows of X must already carry the
/// normalisation dn (X_n = dn X). Claude Generated (Oct 2026).
struct FitFactor {
    bool chol = true;
    MatX L;  ///< cholesky: lower factor
    MatX U;  ///< eigen: kept eigenvectors
    VecX w;  ///< eigen: kept eigenvalues
    /// X := F X (rows = domain size before, F rows after)
    void fit(MatX& X) const
    {
        if (chol)
            L.triangularView<Eigen::Lower>().solveInPlace(X);
        else
            X = (w.cwiseSqrt().cwiseInverse().asDiagonal() * (U.transpose() * X)).eval();
    }
    /// X := Jn^-1 X (normalised metric of the domain)
    void solve(MatX& X) const
    {
        if (chol) {
            L.triangularView<Eigen::Lower>().solveInPlace(X);
            L.transpose().triangularView<Eigen::Upper>().solveInPlace(X);
        } else {
            X = (U * (w.cwiseInverse().asDiagonal() * (U.transpose() * X))).eval();
        }
    }
};

/// Normalised metric block Jn(F) = dn J dn over the global aux indices fidx, factorised. Cholesky unless the
/// engine asks for eigen or the Cholesky fails (then eigen with the cutoff, as RIJK::setup).
bool makeFitFactor(const Matrix& J, const VecX& dn, const std::vector<int>& fidx, bool want_eigen, double cutoff, FitFactor& f)
{
    const int m = (int)fidx.size();
    MatX Jn(m, m);
    for (int q = 0; q < m; ++q)
        for (int r = 0; r < m; ++r) Jn(r, q) = J(fidx[r], fidx[q]) * dn(fidx[r]) * dn(fidx[q]);
    if (!want_eigen) {
        Eigen::LLT<MatX> llt(Jn);
        if (llt.info() == Eigen::Success) {
            f.chol = true;
            f.L = llt.matrixL();
            return true;
        }
    }
    Eigen::SelfAdjointEigenSolver<MatX> es(Jn);
    if (es.info() != Eigen::Success) return false;
    int keep = 0;
    for (int k = 0; k < m; ++k)
        if (es.eigenvalues()(k) > cutoff) ++keep;
    f.chol = false;
    f.U.resize(m, keep);
    f.w.resize(keep);
    for (int k = 0, c = 0; k < m; ++k)
        if (es.eigenvalues()(k) > cutoff) {
            f.U.col(c) = es.eigenvectors().col(k);
            f.w(c++) = es.eigenvalues()(k);
        }
    return true;
}

/// @brief Integral-direct raw factors X_i(P, r) = (P|i r) for P in G_i (atoms g_atoms[i], rows atom by atom in
/// ascending atom order, x_atom_row[i][A] = first row of atom A or -1) and r in ext_ao[i] (columns in that order).
/// See the block comment above. Returns the number of computed (aux shell | shell pair) blocks.
/// Claude Generated (Oct 2026).
long localThreeIndex(const qmri::ThreeIndexDirect& eri, const MatX& L, const MatX& Cpao, double ao_cut,
                     const std::vector<std::vector<int>>& ext_ao, const std::vector<std::vector<int>>& g_atoms,
                     const std::vector<int>& aux_atom_first, const std::vector<std::vector<int>>& x_atom_row, int natm,
                     int threads, std::vector<Matrix>& X, double t_split[3])
{
    const int n = (int)L.rows(), o = (int)L.cols(), ns = eri.nShells();
    std::vector<int> ao_shell(n, -1);
    for (int A = 0; A < ns; ++A)
        for (int k = 0; k < eri.shellWidth(A); ++k) ao_shell[eri.shellFirst(A) + k] = A;
    // truncated coefficients (exact zeros for |C| <= ao_cut)
    MatX Lt = L, Ct = Cpao;
    for (Eigen::Index k = 0; k < Lt.size(); ++k)
        if (std::abs(Lt.data()[k]) <= ao_cut) Lt.data()[k] = 0.0;
    for (Eigen::Index k = 0; k < Ct.size(); ++k)
        if (std::abs(Ct.data()[k]) <= ao_cut) Ct.data()[k] = 0.0;
    // shell supports: of each orbital (m side) and of the PAOs of its extended domain (n side)
    std::vector<std::vector<int>> lsh(o), nsh(o), nsao(o);
    std::vector<MatX> CP(o);  // CP[i] = Ct(nsao[i], ext_ao[i])
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < o; ++i) {
        std::vector<char> ml(ns, 0), mn(ns, 0);
        for (int mu = 0; mu < n; ++mu)
            if (Lt(mu, i) != 0.0) ml[ao_shell[mu]] = 1;
        for (int r : ext_ao[i])
            for (int nu = 0; nu < n; ++nu)
                if (Ct(nu, r) != 0.0) mn[ao_shell[nu]] = 1;
        for (int A = 0; A < ns; ++A) {
            if (ml[A]) lsh[i].push_back(A);
            if (mn[A]) {
                nsh[i].push_back(A);
                for (int k = 0; k < eri.shellWidth(A); ++k) nsao[i].push_back(eri.shellFirst(A) + k);
            }
        }
        const auto& e = ext_ao[i];
        CP[i].resize(nsao[i].size(), e.size());
        for (size_t c = 0; c < e.size(); ++c)
            for (size_t r = 0; r < nsao[i].size(); ++r) CP[i](r, c) = Ct(nsao[i][r], e[c]);
    }
    // orbitals whose G_i contains atom K
    std::vector<std::vector<int>> orbs_of_atom(natm);
    for (int i = 0; i < o; ++i)
        for (int A : g_atoms[i]) orbs_of_atom[A].push_back(i);
    const int nc = eri.nAuxShells();
    long nblocks = 0;
    double t_int = 0.0, t_y = 0.0, t_x = 0.0;
#pragma omp parallel num_threads(threads) reduction(+ : nblocks, t_int, t_y, t_x)
    {
        std::vector<char> inM(ns, 0), inU(ns, 0);
        std::vector<int> mrow(ns, -1), ucol(ns, -1);
        std::vector<double> blk;
        std::vector<std::pair<int, int>> mirror;
#pragma omp for schedule(dynamic, 1)
        for (int c = 0; c < nc; ++c) {
            const int K = eri.auxAtom(c);
            const std::vector<int>& I = orbs_of_atom[K];
            if (I.empty()) continue;
            const int wc = eri.auxWidth(c);
            std::vector<int> Msh, Ush;
            for (int i : I) {
                for (int A : lsh[i]) inM[A] = 1;
                for (int B : nsh[i]) inU[B] = 1;
            }
            int nM = 0, nU = 0;
            for (int A = 0; A < ns; ++A) {
                if (inM[A]) {
                    Msh.push_back(A);
                    mrow[A] = nM;
                    nM += eri.shellWidth(A);
                }
                if (inU[A]) {
                    Ush.push_back(A);
                    ucol[A] = nU;
                    nU += eri.shellWidth(A);
                }
            }
            // T(m, p * nU + u) = (p|m u), m in the orbital supports, u in the PAO supports
            const auto t0 = std::chrono::steady_clock::now();
            MatX T = MatX::Zero(nM, (Eigen::Index)wc * nU);
            mirror.clear();
            for (int A : Msh) {
                const int wA = eri.shellWidth(A);
                for (int B : eri.neighbours(A)) {
                    if (!inU[B]) continue;
                    if (A > B && inM[B] && inU[A]) {  // (B, A) is computed as well: copy it below
                        mirror.push_back({ A, B });
                        continue;
                    }
                    if (!eri.block(c, A, B, blk)) continue;
                    ++nblocks;
                    const int wB = eri.shellWidth(B);
                    for (int a = 0; a < wA; ++a)
                        for (int b = 0; b < wB; ++b) {
                            const double* v = &blk[((size_t)a * wB + b) * wc];
                            for (int p = 0; p < wc; ++p) T(mrow[A] + a, (Eigen::Index)p * nU + ucol[B] + b) = v[p];
                        }
                }
            }
            for (const auto& [A, B] : mirror) {
                const int wA = eri.shellWidth(A), wB = eri.shellWidth(B);
                for (int p = 0; p < wc; ++p)
                    for (int a = 0; a < wA; ++a)
                        for (int b = 0; b < wB; ++b)
                            T(mrow[A] + a, (Eigen::Index)p * nU + ucol[B] + b) = T(mrow[B] + b, (Eigen::Index)p * nU + ucol[A] + a);
            }
            t_int += msSince(t0);
            // Y = L_M^T T: (|I| x wc nU)
            const auto t1 = std::chrono::steady_clock::now();
            MatX LM(nM, (Eigen::Index)I.size());
            for (int A : Msh)
                for (int a = 0; a < eri.shellWidth(A); ++a)
                    for (size_t ii = 0; ii < I.size(); ++ii) LM(mrow[A] + a, ii) = Lt(eri.shellFirst(A) + a, I[ii]);
            const MatX Y = LM.transpose() * T;
            t_y += msSince(t1);
            const auto t2 = std::chrono::steady_clock::now();
            for (size_t ii = 0; ii < I.size(); ++ii) {
                const int i = I[ii];
                const int row0 = x_atom_row[i][K] + (eri.auxFirst(c) - aux_atom_first[K]);
                const auto& nl = nsao[i];
                MatX Yi(wc, (Eigen::Index)nl.size());
                for (size_t k = 0; k < nl.size(); ++k) {
                    const int mu = nl[k], A = ao_shell[mu];
                    const Eigen::Index u = ucol[A] + (mu - eri.shellFirst(A));
                    for (int p = 0; p < wc; ++p) Yi(p, k) = Y(ii, (Eigen::Index)p * nU + u);
                }
                X[i].middleRows(row0, wc).noalias() = Yi * CP[i];
            }
            t_x += msSince(t2);
            for (int A : Msh) {
                inM[A] = 0;
                mrow[A] = -1;
            }
            for (int A : Ush) {
                inU[A] = 0;
                ucol[A] = -1;
            }
        }
    }
    t_split[0] = t_int;
    t_split[1] = t_y;
    t_split[2] = t_x;
    return nblocks;
}

}  // namespace

// =====================================================================================================
// 4. QMEngine::computeLMP2
// =====================================================================================================

// Claude Generated (Oct 2026, WP15). Needs a converged closed-shell RHF (Calculation() first).
bool QMEngine::computeLMP2(const qmint::LMP2Options& opt, qmint::LMP2Result& res)
{
    res = qmint::LMP2Result();
    if (m_functional != QMFunctional::HF) {
        CurcumaLogger::error("LMP2: needs the HF reference (QMFunctional::HF)");
        return false;
    }
    if (!m_scf_converged) {
        CurcumaLogger::error("LMP2: the RHF SCF is not converged -- no LMP2 energy");
        return false;
    }
    if (opt.localizer != "pm" && opt.localizer != "boys") {
        CurcumaLogger::error("LMP2: unknown -qm.lmp2_localizer '" + opt.localizer + "' (use pm | boys)");
        return false;
    }
    if (opt.localizer_algo != "fast" && opt.localizer_algo != "classic") {
        CurcumaLogger::error("LMP2: unknown -qm.lmp2_localizer_algo '" + opt.localizer_algo + "' (use fast | classic)");
        return false;
    }
    if (!(opt.localizer_screen >= 0.0)) {
        CurcumaLogger::error("LMP2: -qm.lmp2_localizer_screen must be >= 0");
        return false;
    }
    const bool fast_loc = opt.localizer_algo == "fast";
    // Claude Generated (Oct 27 2026): as CCSD, the CPU-only part of a -gpu run is a counted fallback (-gpu_strict
    // stops here); the RHF before it may have run on the device.
    if (gpuRequested())
        curcuma::reportGpuFallback("qm_lmp2", "LMP2 has no device implementation; localisation, PNOs and pair energies on the CPU");
    if (opt.distant_pairs != "dipole" && opt.distant_pairs != "drop") {
        CurcumaLogger::error("LMP2: unknown -qm.lmp2_distant_pairs '" + opt.distant_pairs + "' (use dipole | drop)");
        return false;
    }
    if (opt.fit != "global" && opt.fit != "local") {
        CurcumaLogger::error("LMP2: unknown -qm.lmp2_fit '" + opt.fit + "' (use global | local)");
        return false;
    }
    if (opt.fit_domain != "pair" && opt.fit_domain != "orbital") {
        CurcumaLogger::error("LMP2: unknown -qm.lmp2_fit_domain '" + opt.fit_domain + "' (use pair | orbital)");
        return false;
    }
    if (opt.fit == "local" && (!(opt.fit_radius > 0.0) || !(opt.ao_cut >= 0.0))) {
        CurcumaLogger::error("LMP2: -qm.lmp2_fit_radius must be > 0 and -qm.lmp2_ao_cut >= 0");
        return false;
    }
    const bool local_fit = opt.fit == "local";
    const bool robust_fit = local_fit && opt.fit_domain == "orbital";
    const int n = m_nbf, natm = m_atomcount;
    const int n_occ = m_num_electrons / 2;
    int n_frozen = 0;
    if (opt.frozen_core)
        for (int Z : m_atoms) n_frozen += qmint::frozenCoreOrbitals(Z);
    if (n_frozen >= n_occ && n_occ > 0) {
        CurcumaLogger::error(fmt::format("LMP2: frozen core ({} orbitals) leaves no active occupied orbital ({} occupied)",
                                         n_frozen, n_occ));
        return false;
    }
    if (m_mo.rows() != n || m_mo.cols() != n || m_energies.size() != n) {
        CurcumaLogger::error("LMP2: no MO coefficients / orbital energies from the SCF");
        return false;
    }
    const int o = n_occ - n_frozen, v = n - n_occ;
    const int threads = std::max(1, m_threads);
    res.n_frozen = n_frozen;
    res.n_occ_active = o;
    res.n_virt = v;
    res.n_basis = n;
    res.n_atoms = natm;
    if (riActive() && m_ri) res.scf_ri = fmt::format("RI-{} {}", m_ri_mode == "j" ? "J" : "JK", auxLabel(m_ri_file));
    std::string file;
    try {
        file = qmri::resolveAuxBasisFile(opt.ri_aux, "jk");
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("LMP2: ") + e.what());
        return false;
    }
    res.aux_name = auxLabel(file);
    if (o == 0 || v == 0) return true;

    // ---- common matrices (column-major copies) ----
    const MatX C = toEigen(m_mo), S = toEigen(m_S);
    const VecX eps = m_energies;
    const MatX Cocc = C.leftCols(n_occ), Cval = C.middleCols(n_frozen, o);
    const MatX SC = S * C;
    const MatX Fao = SC * eps.asDiagonal() * SC.transpose();  // F consistent with the converged orbitals
    const std::vector<int> ao_atom = activeAtomIndex();
    std::vector<std::vector<int>> atom_aos(natm);
    for (int mu = 0; mu < n; ++mu) {
        if (ao_atom[mu] < 0 || ao_atom[mu] >= natm) {
            CurcumaLogger::error("LMP2: basis function without atom");
            return false;
        }
        atom_aos[ao_atom[mu]].push_back(mu);
    }
    MatX xyz(natm, 3);
    for (int a = 0; a < natm; ++a)
        for (int k = 0; k < 3; ++k) xyz(a, k) = CurcumaUnit::Length::angstrom_to_bohr(m_geometry(a, k));

    // ---- 1. localisation (lmp2_common.localise) ----
    // Sub-timings go to res.localization_info (printed at verbosity 2). Claude Generated (Oct 2026).
    const auto t_loc = std::chrono::steady_clock::now();
    auto t_sub = t_loc;
    std::string loc_info;
    auto lap = [&](const char* name) {
        loc_info += fmt::format("{}{} {:.1f}", loc_info.empty() ? "" : ", ", name, msSince(t_sub));
        t_sub = std::chrono::steady_clock::now();
    };
    std::vector<MatX> Dip(3);  // active-basis dipole matrices (Bohr)
    {
        Matrix Sc, Xc, Yc, Zc;
        cartesianDipoleIntegrals(m_gto_basis, Sc, Xc, Yc, Zc);
        // cartesian -> spherical with the sparse transform (Q^T M Q as dense GEMMs was ~n^3, 0.2 s at n = 970)
        Eigen::SparseMatrix<double> Qs;
        const bool sph = m_Q.size() > 0;
        if (sph) Qs = MatX(m_Q).sparseView();
        auto toActive = [&](const Matrix& Mc) -> MatX {
            if (!sph) return toEigen(Mc);
            const MatX MQ = toEigen(Mc) * Qs;
            return MatX(Qs.transpose() * MQ);
        };
        const double dS = (toActive(Sc) - toEigen(m_S)).cwiseAbs().maxCoeff();
        if (dS > 1e-8) {  // self-check of the integral code against the engine's overlap
            CurcumaLogger::error(fmt::format("LMP2: dipole integral self-check failed (overlap differs by {:.1e})", dS));
            return false;
        }
        Dip[0] = toActive(Xc);
        Dip[1] = toActive(Yc);
        Dip[2] = toActive(Zc);
    }
    lap("dipole integrals");
    MatX L;
    if (!choleskyOrbitals(Cval, S, L)) {
        CurcumaLogger::error("LMP2: Cholesky start orbitals failed (occupied space rank-deficient?)");
        return false;
    }
    lap("Cholesky start");
    {
        double f = 0.0;
        int sweeps = 0;
        if (fast_loc && opt.localizer == "pm") {
            PMSweepStats st;
            sweeps = pmLocalizeFast(L, S, atom_aos, threads, kLocTol, kLocMaxSweeps, opt.localizer_screen, f, st);
            lap("PM sweeps");
            loc_info += fmt::format(" (pairs evaluated {:.1f} %, rotated {:.1f} %, mean support {:.1f} of {} atoms)",
                                    st.pairs_total ? 100.0 * st.pairs_evaluated / st.pairs_total : 0.0,
                                    st.pairs_total ? 100.0 * st.pairs_rotated / st.pairs_total : 0.0, st.mean_support, natm);
        } else {
            std::vector<MatX> M;
            if (opt.localizer == "boys") {
                for (int c = 0; c < 3; ++c) {
                    const MatX DL = productParallel(Dip[c], L, threads);
                    M.push_back(L.transpose() * DL);
                }
            } else {
                const MatX SL = S * L;
                M.resize(natm);
#pragma omp parallel for schedule(dynamic) num_threads(threads)
                for (int A = 0; A < natm; ++A) {
                    MatX q = MatX::Zero(o, o);
                    for (int mu : atom_aos[A]) q += L.row(mu).transpose() * SL.row(mu);
                    M[A] = 0.5 * (q + q.transpose());
                }
            }
            lap("matrices");
            // Boys (3 matrices, O(o^3) per sweep) keeps the matrix route in both algorithms: a pair-parallel variant
            // was slower up to pentacene (parallel-region overhead per round)
            sweeps = jacobiLocalize(L, M, threads, kLocTol, kLocMaxSweeps, f);
            lap("sweeps");
        }
        res.localization_sweeps = sweeps;
        res.localization_functional = f;
        if (sweeps < 0 && CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::warn(fmt::format("LMP2: localisation not converged after {} sweeps (functional {:.10f}); the "
                                            "energy with truncation depends slightly on the orbitals",
                                            -sweeps, f));
    }
    // deterministic order (centroid x, y, z rounded to 1e-6 Bohr) and sign (largest coefficient positive);
    // DipL[c] = Dip[c] L is kept (columns reordered with L) for the distant-pair estimate below
    std::vector<MatX> DipL(3);
    for (int c = 0; c < 3; ++c) DipL[c] = productParallel(Dip[c], L, threads);
    MatX cen(o, 3);
    for (int i = 0; i < o; ++i)
        for (int c = 0; c < 3; ++c) cen(i, c) = L.col(i).dot(DipL[c].col(i));
    {
        std::vector<int> order(o);
        std::iota(order.begin(), order.end(), 0);
        auto key = [&](int i, int c) { return std::llround(cen(i, c) * 1e6); };
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            for (int c = 0; c < 3; ++c)
                if (key(a, c) != key(b, c)) return key(a, c) < key(b, c);
            return false;
        });
        MatX L2(n, o), cen2(o, 3);
        std::vector<MatX> DipL2(3, MatX(n, o));
        for (int k = 0; k < o; ++k) {
            L2.col(k) = L.col(order[k]);
            cen2.row(k) = cen.row(order[k]);
            for (int c = 0; c < 3; ++c) DipL2[c].col(k) = DipL[c].col(order[k]);
            Eigen::Index imax;
            L2.col(k).cwiseAbs().maxCoeff(&imax);
            if (L2(imax, k) < 0) {
                L2.col(k) *= -1.0;
                for (int c = 0; c < 3; ++c) DipL2[c].col(k) *= -1.0;
            }
        }
        L.swap(L2);
        cen.swap(cen2);
        DipL.swap(DipL2);
    }
    lap("centroids");
    const MatX U = Cval.transpose() * productParallel(S, L, threads);  // L = Cval U
    const MatX Focc = U.transpose() * eps.segment(n_frozen, o).asDiagonal() * U;
    lap("U");
    res.t_localize_ms = msSince(t_loc);
    res.localization_info = loc_info;

    // ---- 2. domains (bp_domain + extend_domain) and pairs ----
    const auto t_pairs = std::chrono::steady_clock::now();
    std::vector<std::vector<int>> dom(o), bpdom(o);
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < o; ++i) {
        std::vector<int> bp = bpDomain(L.col(i), S, atom_aos, opt.bp_threshold);
        bpdom[i] = bp;
        std::vector<int> ext;
        for (int A = 0; A < natm; ++A) {
            bool in = false;
            for (int B : bp)
                if ((xyz.row(A) - xyz.row(B)).norm() <= opt.domain_radius) {
                    in = true;
                    break;
                }
            if (in || std::binary_search(bp.begin(), bp.end(), A)) ext.push_back(A);
        }
        dom[i] = ext;
    }
    double dom_sum = 0.0;
    for (const auto& d : dom) dom_sum += d.size();
    res.mean_domain_atoms = dom_sum / o;

    std::vector<int> pair_index((size_t)o * o, -1);  // kept pairs, both orders
    std::vector<std::pair<int, int>> kept, distant;
    for (int i = 0; i < o; ++i)
        for (int j = i; j < o; ++j) {
            const double r = (cen.row(i) - cen.row(j)).norm();
            if (r <= opt.pair_cut) {
                pair_index[(size_t)i * o + j] = pair_index[(size_t)j * o + i] = (int)kept.size();
                kept.push_back({i, j});
            } else {
                distant.push_back({i, j});
            }
        }
    const int npairs = (int)kept.size();
    res.n_pairs_total = (long)o * (o + 1) / 2;
    res.n_pairs_kept = npairs;
    res.n_pairs_distant = (long)distant.size();

    // PAOs (lmp2_common.paos): columns of (1 - D_occ S), D_occ incl. the core, normalised; PAO overlap and Fock
    const MatX Praw = MatX::Identity(n, n) - Cocc * (Cocc.transpose() * S);
    MatX Spao = Praw.transpose() * S * Praw;
    VecX inv_norm(n);
    for (int r = 0; r < n; ++r) inv_norm(r) = Spao(r, r) > 1e-14 ? 1.0 / std::sqrt(Spao(r, r)) : 0.0;
    const MatX Cpao = Praw * inv_norm.asDiagonal();
    Spao = inv_norm.asDiagonal() * Spao * inv_norm.asDiagonal();
    const MatX Fpao = Cpao.transpose() * Fao * Cpao;

    // extended domain of i: union of [j] over all kept partners j (incl. i); B_{i,r} is kept for its PAOs
    std::vector<std::vector<int>> ext_ao(o);
    std::vector<std::vector<int>> ext_pos(o);  // AO -> column in B[i], -1 if absent
    for (int i = 0; i < o; ++i) {
        std::vector<char> atom_in(natm, 0);
        for (int j = 0; j < o; ++j)
            if (pair_index[(size_t)i * o + j] >= 0)
                for (int A : dom[j]) atom_in[A] = 1;
        ext_pos[i].assign(n, -1);
        for (int mu = 0; mu < n; ++mu)
            if (atom_in[ao_atom[mu]]) {
                ext_pos[i][mu] = (int)ext_ao[i].size();
                ext_ao[i].push_back(mu);
            }
    }
    std::vector<size_t> ext_off(o + 1, 0);
    for (int i = 0; i < o; ++i) ext_off[i + 1] = ext_off[i] + ext_ao[i].size();
    const size_t ncol = ext_off[o];

    const double t_domains_ms = msSince(t_pairs);  // domains, pairs, PAOs

    // ---- 3. 3-index factors ----
    // global: fitted B_{i,r}^P from the stored RI tensor (QMEngine::correlationRI);
    // local : raw X_i(P, r) = (P|i r), P in G_i, integral-direct (localThreeIndex), fitted per pair in step 4.
    std::string ri_src, why;
    double b_mb = 0.0;
    int naux = 0;
    std::vector<MatX> Bf;                       // global: Bf[i] naux x |ext_ao[i]|, column r = B_{i, ext_ao[i][r]}^P
    std::vector<Matrix> Xl;                     // local: X_i, rows G_i (atom by atom), columns ext_ao[i]
    std::vector<std::vector<int>> f_atoms(o);   // local: fit domain F_i (atoms)
    std::vector<std::vector<int>> x_atom_row;   // local: first row of atom A in X_i (-1: absent)
    std::vector<std::vector<int>> g_atoms(o);   // local: atoms of the rows of X_i (G_i = union of partner F_j)
    std::vector<int> aux_atom_first(natm + 1, 0);
    Matrix Jraw;  // local: raw Coulomb metric (naux x naux)
    VecX dn;      // local: 1/sqrt((P|P))
    if (!local_fit) {
        double setup_ms = 0.0;
        std::shared_ptr<qmri::RIJK> ri = correlationRI(file, ri_src, setup_ms, why);
        if (!ri) {
            CurcumaLogger::error("LMP2: density fitting setup with " + file + " not possible: " + why);
            return false;
        }
        res.t_ri_setup_ms = setup_ms;
        const int nrows = (int)ri->tensor().rows();
        b_mb = (double)ncol * nrows * 8.0 / 1.0e6;
        if (2.0 * b_mb > opt.max_memory_mb) {
            CurcumaLogger::error(fmt::format("LMP2: the factors B_(i,r)^P need 2 x {:.1f} MB ({} columns over the extended domains x {} "
                                             "auxiliary functions, plus one copy during the fit), above -qm.lmp2_max_memory_mb {:g}",
                                             b_mb, ncol, nrows, opt.max_memory_mb));
            return false;
        }
        const auto t_tr = std::chrono::steady_clock::now();
        Bf.resize(o);
        {
            Matrix Xr(nrows, (Eigen::Index)ncol);  // (Q|i r), row Q
            const MatX Lt = L.transpose();
#pragma omp parallel num_threads(threads)
            {
                MatX TQ(n, n), Y(o, n), Zq(o, n);
#pragma omp for schedule(dynamic, 1)
                for (int Q = 0; Q < nrows; ++Q) {
                    ri->unpackRow(Q, TQ.data());  // symmetric: row- and column-major identical
                    Y.noalias() = Lt * TQ;
                    Zq.noalias() = Y * Cpao;  // (Q| i r) with PAO r
                    double* row = Xr.data() + (size_t)Q * ncol;
                    for (int i = 0; i < o; ++i) {
                        double* dst = row + ext_off[i];
                        const auto& e = ext_ao[i];
                        for (size_t r = 0; r < e.size(); ++r) dst[r] = Zq(i, e[r]);
                    }
                }
            }
            ri->fitRows(Xr);  // metric factor (no-op when the SCF tensor already is B)
            ri.reset();       // an own RI setup is released here (the SCF one stays with the engine)
            naux = (int)Xr.rows();
#pragma omp parallel for schedule(dynamic) num_threads(threads)
            for (int i = 0; i < o; ++i) Bf[i] = Xr.block(0, (Eigen::Index)ext_off[i], naux, (Eigen::Index)ext_ao[i].size());
        }
        res.n_aux = naux;
        res.t_transform_ms = msSince(t_tr);
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("LMP2: {}; B factors {:.0f} MB ({} columns)", ri_src, b_mb, ncol));
    } else {
        // local fit: auxiliary shells, integral-direct engine, metric (no 3-index tensor anywhere)
        const auto t_setup = std::chrono::steady_clock::now();
        qmri::ThreeIndexDirect eri;
        try {
            const qmri::AuxBasisMap aux = qmri::parseTurbomoleBasis(file);
            Matrix pos(natm, 3);
            for (int a = 0; a < natm; ++a)
                for (int k = 0; k < 3; ++k) pos(a, k) = xyz(a, k);
            const std::vector<qmri::Shell> sh = qmri::auxiliaryShells(aux, m_atoms, pos);
            if (!eri.setup(m_gto_basis, m_Q.size() != 0 ? &m_Q : nullptr, sh, !m_cartesian_d, m_eri_screening, threads, why)) {
                CurcumaLogger::error("LMP2: integral-direct 3-index setup failed: " + why);
                return false;
            }
            Jraw = qmri::coulombMetric(sh, !m_cartesian_d, threads);
        } catch (const std::exception& e) {
            CurcumaLogger::error(std::string("LMP2: local fit setup with ") + file + ": " + e.what());
            return false;
        }
        if (eri.n() != n) {
            CurcumaLogger::error(fmt::format("LMP2: integral-direct basis has {} functions, the engine {}", eri.n(), n));
            return false;
        }
        naux = eri.naux();
        res.n_aux = naux;
        dn.resize(naux);
        for (int P = 0; P < naux; ++P) {
            if (!(Jraw(P, P) > 0.0)) {
                CurcumaLogger::error(fmt::format("LMP2: auxiliary function {} has (P|P) = {:.3e}", P, Jraw(P, P)));
                return false;
            }
            dn(P) = 1.0 / std::sqrt(Jraw(P, P));
        }
        // auxiliary functions of one atom are contiguous (auxiliaryShells): first function per atom
        std::vector<int> cnt(natm, 0);
        for (int c = 0; c < eri.nAuxShells(); ++c) cnt[eri.auxAtom(c)] += eri.auxWidth(c);
        for (int A = 0; A < natm; ++A) aux_atom_first[A + 1] = aux_atom_first[A] + cnt[A];
        for (int c = 0; c < eri.nAuxShells(); ++c) {
            const int A = eri.auxAtom(c);
            if (eri.auxFirst(c) < aux_atom_first[A] || eri.auxFirst(c) + eri.auxWidth(c) > aux_atom_first[A + 1]) {
                CurcumaLogger::error("LMP2: auxiliary functions of an atom are not contiguous");
                return false;
            }
        }
        // fit domains F_i: atoms within fit_radius of a Boughton-Pulay atom of i; G_i: union of F_j over kept partners
        for (int i = 0; i < o; ++i)
            for (int A = 0; A < natm; ++A)
                for (int B : bpdom[i])
                    if (A == B || (xyz.row(A) - xyz.row(B)).norm() <= opt.fit_radius) {
                        f_atoms[i].push_back(A);
                        break;
                    }
        x_atom_row.assign(o, std::vector<int>(natm, -1));
        double x_mb = 0.0, cp_mb = 0.0, d_mb = 0.0, rows_sum = 0.0, f_sum = 0.0;
        std::vector<int> x_rows(o, 0);
        for (int i = 0; i < o; ++i) {
            std::vector<char> in(natm, 0);
            for (int j = 0; j < o; ++j)
                if (pair_index[(size_t)i * o + j] >= 0)
                    for (int A : f_atoms[j]) in[A] = 1;
            for (int A = 0; A < natm; ++A)
                if (in[A]) {
                    g_atoms[i].push_back(A);
                    x_atom_row[i][A] = x_rows[i];
                    x_rows[i] += cnt[A];
                }
            int fi = 0;
            for (int A : f_atoms[i]) fi += cnt[A];
            x_mb += (double)x_rows[i] * ext_ao[i].size() * 8.0 / 1e6;
            cp_mb += (double)n * ext_ao[i].size() * 8.0 / 1e6;  // upper bound of the PAO blocks (support <= n)
            if (robust_fit) d_mb += (double)(fi + x_rows[i]) * ext_ao[i].size() * 8.0 / 1e6;  // d_i and Z_i
            rows_sum += x_rows[i];
            f_sum += fi;
        }
        res.mean_x_rows = rows_sum / o;
        if (robust_fit) res.mean_fit_aux = f_sum / o;
        const double metric_mb = (double)naux * naux * 8.0 / 1e6;
        res.factors_mb = x_mb + d_mb + metric_mb;
        if (x_mb + cp_mb + d_mb + metric_mb > opt.max_memory_mb) {
            CurcumaLogger::error(fmt::format("LMP2 (local fit): the 3-index factors (P|i r) need {:.1f} MB, the PAO blocks up to {:.1f} MB, "
                                             "the fit coefficients {:.1f} MB and the metric {:.1f} MB, above -qm.lmp2_max_memory_mb {:g}",
                                             x_mb, cp_mb, d_mb, metric_mb, opt.max_memory_mb));
            return false;
        }
        res.t_ri_setup_ms = msSince(t_setup);
        const auto t_tr = std::chrono::steady_clock::now();
        Xl.resize(o);
        for (int i = 0; i < o; ++i) Xl[i] = Matrix::Zero(x_rows[i], (Eigen::Index)ext_ao[i].size());
        double tsplit[3] = { 0.0, 0.0, 0.0 };
        res.n_int_blocks =
            localThreeIndex(eri, L, Cpao, opt.ao_cut, ext_ao, g_atoms, aux_atom_first, x_atom_row, natm, threads, Xl, tsplit);
        res.t_transform_ms = msSince(t_tr);
        ri_src = fmt::format("local fit ({} domains, radius {:g} Bohr, ao_cut {:g}), integral-direct 3-index integrals ({} shell pairs, "
                             "{} blocks; thread-summed ms: integrals {:.0f}, L transformation {:.0f}, PAO transformation {:.0f}), X rows mean {:.0f} of {}",
                             opt.fit_domain, opt.fit_radius, opt.ao_cut, eri.nPairs(), res.n_int_blocks, tsplit[0], tsplit[1],
                             tsplit[2], res.mean_x_rows, naux);
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("LMP2: {}; factors {:.0f} MB, metric {:.0f} MB", ri_src, x_mb + d_mb, metric_mb));
    }
    res.fit_label = local_fit ? "local-" + opt.fit_domain : "global";

    // ---- 4. pair spaces, K in the PAO space, PNOs (lmp2_common.pair_space, pair_ints, make_pnos) ----
    std::vector<PairData> P(npairs);
    std::vector<char> pair_ok(npairs, 1);
    // pair space of a kept pair (PAOs of [i] u [j], pseudo-canonical)
    auto pairSpace = [&](int p, PaoSpace& sp) {
        const int i = kept[p].first, j = kept[p].second;
        PairData& pd = P[p];
        pd.i = i;
        pd.j = j;
        std::vector<int> atoms = dom[i];
        atoms.insert(atoms.end(), dom[j].begin(), dom[j].end());
        std::sort(atoms.begin(), atoms.end());
        atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
        for (int A : atoms) pd.idx.insert(pd.idx.end(), atom_aos[A].begin(), atom_aos[A].end());
        std::sort(pd.idx.begin(), pd.idx.end());
        if (!paoSpace(pd.idx, Spao, Fpao, sp)) return false;
        pd.m_pao = (int)sp.eps.size();
        return true;
    };
    // semicanonical amplitudes, pair density, PNOs from K = (ia|jb) in the pseudo-canonical pair space
    auto finishPair = [&](int p, const PaoSpace& sp, const MatX& K) {
        PairData& pd = P[p];
        const int i = pd.i, j = pd.j, m = (int)sp.eps.size();
        const double Fii = Focc(i, i), Fjj = Focc(j, j);
        MatX den(m, m);
        for (int b = 0; b < m; ++b)
            for (int a = 0; a < m; ++a) den(a, b) = sp.eps(a) + sp.eps(b) - Fii - Fjj;
        const MatX T = -K.cwiseQuotient(den);
        pd.e_sc_base = pairEnergy(T, K, i == j);
        // pair density D = T~^T T + T~ T^T, T~ = (4T - 2T^T) / (1 + delta_ij)
        const MatX Tt = (4.0 * T - 2.0 * T.transpose()) / (i == j ? 2.0 : 1.0);
        MatX D = Tt.transpose() * T + Tt * T.transpose();
        D = 0.5 * (D + D.transpose()).eval();
        Eigen::SelfAdjointEigenSolver<MatX> ed(D);
        std::vector<int> keep;
        for (int r = 0; r < m; ++r)
            if (opt.pno_cut <= 0.0 || ed.eigenvalues()(r) >= opt.pno_cut) keep.push_back(r);
        const int npno = (int)keep.size();
        MatX Qm(m, npno);
        for (int r = 0; r < npno; ++r) Qm.col(r) = ed.eigenvectors().col(keep[r]);
        MatX QW = Qm;
        if (npno > 0) {
            Eigen::SelfAdjointEigenSolver<MatX> ep(Qm.transpose() * sp.eps.asDiagonal() * Qm);
            QW = Qm * ep.eigenvectors();
            pd.eps = ep.eigenvalues();
        }
        pd.K = QW.transpose() * K * QW;
        pd.den.resize(npno, npno);
        for (int b = 0; b < npno; ++b)
            for (int a = 0; a < npno; ++a) pd.den(a, b) = pd.eps(a) + pd.eps(b) - Fii - Fjj;
        pd.e_sc_pno = npno ? pairEnergy(-pd.K.cwiseQuotient(pd.den), pd.K, i == j) : 0.0;
        pd.c = sp.Xc * QW;
    };
    // local fit: global aux indices of an atom set, and rows `atoms` of X_k over the pair space (scaled by dn if asked)
    auto auxIndices = [&](const std::vector<int>& atoms) {
        std::vector<int> f;
        for (int A : atoms)
            for (int q = aux_atom_first[A]; q < aux_atom_first[A + 1]; ++q) f.push_back(q);
        return f;
    };
    auto localRows = [&](const std::vector<Matrix>& Xs, int k, const std::vector<int>& atoms, const PairData& pd, const PaoSpace& sp,
                         bool scale) {
        const int mao = (int)pd.idx.size();
        int f = 0;
        for (int A : atoms) f += aux_atom_first[A + 1] - aux_atom_first[A];
        std::vector<int> cols(mao);
        for (int r = 0; r < mao; ++r) cols[r] = ext_pos[k][pd.idx[r]];
        MatX G(f, mao);
        int row = 0;
        for (int A : atoms) {
            const int x0 = x_atom_row[k][A], na = aux_atom_first[A + 1] - aux_atom_first[A];
            for (int t = 0; t < na; ++t, ++row) {
                const double* src = Xs[k].data() + (size_t)(x0 + t) * Xs[k].cols();
                const double s = scale ? dn(aux_atom_first[A] + t) : 1.0;
                for (int r = 0; r < mao; ++r) G(row, r) = s * src[cols[r]];
            }
        }
        return MatX(G * sp.Xc);
    };
    const bool want_eigen = m_ri_metric == "eigen";
    std::vector<char> fit_ok(std::max(1, npairs), 1);
    const auto t_pairstep = std::chrono::steady_clock::now();
    if (!local_fit) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (int p = 0; p < npairs; ++p) {
            PaoSpace sp;
            if (!pairSpace(p, sp)) {
                pair_ok[p] = 0;
                continue;
            }
            const PairData& pd = P[p];
            const int i = pd.i, j = pd.j, mao = (int)pd.idx.size();
            // B_ia^P = sum_r Xc_ra B_{i,r}^P over the pair's PAOs
            auto half = [&](int k) {
                MatX G(naux, mao);
                for (int r = 0; r < mao; ++r) G.col(r) = Bf[k].col(ext_pos[k][pd.idx[r]]);
                return MatX(G * sp.Xc);
            };
            const MatX Bi = half(i);
            const MatX K = (i == j) ? MatX(Bi.transpose() * Bi) : MatX(Bi.transpose() * half(j));
            finishPair(p, sp, K);
        }
    } else if (!robust_fit) {
        // pair fit domains F_ij = F_i u F_j; pairs with the same domain share one metric factor
        std::map<std::vector<int>, std::vector<int>> by_domain;
        for (int p = 0; p < npairs; ++p) {
            const int i = kept[p].first, j = kept[p].second;
            std::vector<int> key = f_atoms[i];
            key.insert(key.end(), f_atoms[j].begin(), f_atoms[j].end());
            std::sort(key.begin(), key.end());
            key.erase(std::unique(key.begin(), key.end()), key.end());
            by_domain[key].push_back(p);
        }
        std::vector<std::pair<std::vector<int>, std::vector<int>>> groups(by_domain.begin(), by_domain.end());
        std::map<std::vector<int>, std::vector<int>>().swap(by_domain);
        res.n_fit_factors = (long)groups.size();
        double fsum = 0.0, fit_ms = 0.0;
        // chunks of groups: factors of one chunk in parallel, then all pairs of the chunk in parallel (memory: at
        // most `chunk` factors at a time; a single large domain still threads over its pairs)
        const int ngroups = (int)groups.size(), chunk = std::max(1, 4 * threads);
        for (int g0 = 0; g0 < ngroups; g0 += chunk) {
            const int g1 = std::min(ngroups, g0 + chunk);
            std::vector<FitFactor> ff(g1 - g0);
            std::vector<std::vector<int>> fidx(g1 - g0);
            std::vector<char> gok(g1 - g0, 1);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+ : fit_ms)
            for (int g = g0; g < g1; ++g) {
                const auto t0 = std::chrono::steady_clock::now();
                fidx[g - g0] = auxIndices(groups[g].first);
                gok[g - g0] = makeFitFactor(Jraw, dn, fidx[g - g0], want_eigen, m_ri_metric_cutoff, ff[g - g0]) ? 1 : 0;
                fit_ms += msSince(t0);
            }
            std::vector<std::pair<int, int>> work;  // (pair, group - g0)
            for (int g = g0; g < g1; ++g)
                for (int p : groups[g].second) work.push_back({ p, g - g0 });
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+ : fsum, fit_ms)
            for (int w = 0; w < (int)work.size(); ++w) {
                const int p = work[w].first, gl = work[w].second;
                const std::vector<int>& atoms = groups[g0 + gl].first;
                fsum += (double)fidx[gl].size();
                if (!gok[gl]) {
                    fit_ok[p] = 0;
                    continue;
                }
                PaoSpace sp;
                if (!pairSpace(p, sp)) {
                    pair_ok[p] = 0;
                    continue;
                }
                const PairData& pd = P[p];
                const auto t1 = std::chrono::steady_clock::now();
                MatX Bi = localRows(Xl, pd.i, atoms, pd, sp, true);
                ff[gl].fit(Bi);
                MatX K;
                if (pd.i == pd.j) {
                    K = Bi.transpose() * Bi;
                } else {
                    MatX Bj = localRows(Xl, pd.j, atoms, pd, sp, true);
                    ff[gl].fit(Bj);
                    K = Bi.transpose() * Bj;
                }
                fit_ms += msSince(t1);
                finishPair(p, sp, K);
            }
        }
        res.mean_fit_aux = npairs ? fsum / npairs : 0.0;
        res.t_fit_ms = fit_ms;  // summed over threads
    } else {
        // robust fit with orbital domains, once per orbital over all columns of X_i:
        //   d_i = J_Fi^-1 X_i(F_i),   Z_i = X_i - J(G_i, F_i) d_i   (rows G_i, the residual of the fit of (P|i r))
        // so that per pair (ia|jb) = d_i^T Z_j(F_i) + X_i(F_j)^T d_j, which is the robust expression
        //   d_i^T X_j(F_i) + X_i(F_j)^T d_j - d_i^T J(F_i, F_j) d_j   (F_i is part of G_j for every partner).
        std::vector<MatX> Dl(o);
        std::vector<Matrix> Zl(o);
        std::vector<std::vector<int>> fidx(o);
        std::vector<char> ook(o, 1);
        res.n_fit_factors = o;
        double fit_ms = 0.0;
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+ : fit_ms)
        for (int i = 0; i < o; ++i) {
            const auto t0 = std::chrono::steady_clock::now();
            fidx[i] = auxIndices(f_atoms[i]);
            FitFactor ff;
            if (!makeFitFactor(Jraw, dn, fidx[i], want_eigen, m_ri_metric_cutoff, ff)) {
                ook[i] = 0;
                continue;
            }
            const Eigen::Index nc = Xl[i].cols();
            MatX Y((Eigen::Index)fidx[i].size(), nc);
            int row = 0;
            for (int A : f_atoms[i]) {
                const int x0 = x_atom_row[i][A], na = aux_atom_first[A + 1] - aux_atom_first[A];
                for (int t = 0; t < na; ++t, ++row) Y.row(row) = dn(aux_atom_first[A] + t) * Xl[i].row(x0 + t);
            }
            ff.solve(Y);  // Jn^-1 dn X
            for (Eigen::Index r = 0; r < Y.rows(); ++r) Y.row(r) *= dn(fidx[i][r]);  // d = dn Jn^-1 dn X = J^-1 X
            // Z_i = X_i - J(G_i, F_i) d_i
            std::vector<int> gidx = auxIndices(g_atoms[i]);
            MatX Jg((Eigen::Index)gidx.size(), (Eigen::Index)fidx[i].size());
            for (size_t q = 0; q < fidx[i].size(); ++q)
                for (size_t r = 0; r < gidx.size(); ++r) Jg(r, q) = Jraw(gidx[r], fidx[i][q]);
            Zl[i] = Xl[i];
            Zl[i].noalias() -= Jg * Y;
            Dl[i].swap(Y);
            fit_ms += msSince(t0);
        }
        for (int i = 0; i < o; ++i)
            if (!ook[i]) {
                CurcumaLogger::error("LMP2: factorisation of an orbital fit domain metric failed");
                return false;
            }
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+ : fit_ms)
        for (int p = 0; p < npairs; ++p) {
            PaoSpace sp;
            if (!pairSpace(p, sp)) {
                pair_ok[p] = 0;
                continue;
            }
            const PairData& pd = P[p];
            const int i = pd.i, j = pd.j, mao = (int)pd.idx.size();
            const auto t0 = std::chrono::steady_clock::now();
            auto coef = [&](int k) {  // d_k over the pair space: |F_k| x m
                MatX G(Dl[k].rows(), mao);
                for (int r = 0; r < mao; ++r) G.col(r) = Dl[k].col(ext_pos[k][pd.idx[r]]);
                return MatX(G * sp.Xc);
            };
            const MatX Di = coef(i);
            const MatX Zj_Fi = localRows(Zl, j, f_atoms[i], pd, sp, false);  // (P|jb) - (P|fitted jb), P in F_i
            const MatX Xi_Fj = localRows(Xl, i, f_atoms[j], pd, sp, false);  // (ia|Q), Q in F_j
            MatX K = Di.transpose() * Zj_Fi;
            if (i == j)
                K += Xi_Fj.transpose() * Di;
            else
                K += Xi_Fj.transpose() * coef(j);
            if (i == j) K = (0.5 * (K + K.transpose())).eval();
            fit_ms += msSince(t0);
            finishPair(p, sp, K);
        }
        res.t_fit_ms = fit_ms;  // summed over threads
    }
    for (int p = 0; p < npairs; ++p)
        if (!pair_ok[p] || !fit_ok[p]) {
            CurcumaLogger::error(!pair_ok[p] ? "LMP2: eigensolver failed in a pair space"
                                             : "LMP2: factorisation of a pair fit domain metric failed");
            return false;
        }
    std::vector<MatX>().swap(Bf);
    std::vector<Matrix>().swap(Xl);
    Matrix().swap(Jraw);
    long amp = 0, pno_sum = 0, pao_sum = 0;
    double pno_mb = 0.0;
    for (const auto& pd : P) {
        const long np = pd.eps.size();
        amp += np * np;
        pno_sum += np;
        pao_sum += pd.m_pao;
        res.max_pno = std::max(res.max_pno, (int)np);
        pno_mb += (pd.c.size() + 3.0 * np * np) * 8.0 / 1e6;  // c, K, den, T
        if (opt.pno_correction) res.e_pno_corr += pd.e_sc_base - pd.e_sc_pno;
    }
    res.mean_pno = npairs ? (double)pno_sum / npairs : 0.0;
    res.mean_pao = npairs ? (double)pao_sum / npairs : 0.0;
    res.amp_fraction = (double)amp / ((double)o * (o + 1) / 2 * (double)v * v);
    // pair step = domains/PAOs + pair spaces/PNOs, without the RI setup and the B transformation in between
    res.t_pairs_ms = t_domains_ms + msSince(t_pairstep);

    // ---- 5. distant pairs: dipole estimate over the PAO spaces of [i] and [j] ----
    const auto t_dist = std::chrono::steady_clock::now();
    if (opt.distant_pairs == "dipole" && !distant.empty()) {
        std::vector<char> need(o, 0);
        for (const auto& pr : distant) need[pr.first] = need[pr.second] = 1;
        std::vector<VecX> oe(o);
        std::vector<MatX> mu(o);  // 3 x m_i: <i|r|a> over the pseudo-canonical PAO space of [i]
        std::vector<MatX> W(3);   // PAO^T r L: n x o
        for (int c = 0; c < 3; ++c) W[c] = productParallelAtB(Cpao, DipL[c], threads);
        std::vector<char> ok(o, 1);
#pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (int i = 0; i < o; ++i) {
            if (!need[i]) continue;
            std::vector<int> idx;
            for (int A : dom[i]) idx.insert(idx.end(), atom_aos[A].begin(), atom_aos[A].end());
            std::sort(idx.begin(), idx.end());
            PaoSpace sp;
            if (!paoSpace(idx, Spao, Fpao, sp)) {
                ok[i] = 0;
                continue;
            }
            oe[i] = sp.eps;
            mu[i].resize(3, sp.eps.size());
            for (int c = 0; c < 3; ++c) {
                VecX w(idx.size());
                for (size_t r = 0; r < idx.size(); ++r) w(r) = W[c](idx[r], i);
                mu[i].row(c) = (sp.Xc.transpose() * w).transpose();
            }
        }
        for (int i = 0; i < o; ++i)
            if (!ok[i]) {
                CurcumaLogger::error("LMP2: eigensolver failed in an orbital domain");
                return false;
            }
        std::vector<double> ed(distant.size(), 0.0);
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
        for (size_t d = 0; d < distant.size(); ++d) {
            const int i = distant[d].first, j = distant[d].second;
            const Eigen::Vector3d Rv = (cen.row(j) - cen.row(i)).transpose();
            const double R = Rv.norm();
            const Eigen::Vector3d nv = Rv / R;
            const Eigen::Matrix3d M = (Eigen::Matrix3d::Identity() - 3.0 * nv * nv.transpose()) / (R * R * R);
            const MatX G = mu[i].transpose() * (M * mu[j]);  // (ia|jb) in the dipole approximation
            const double Fs = Focc(i, i) + Focc(j, j);
            double e = 0.0;
            for (int b = 0; b < G.cols(); ++b)
                for (int a = 0; a < G.rows(); ++a) e += G(a, b) * G(a, b) / (oe[i](a) + oe[j](b) - Fs);
            ed[d] = -4.0 * e;
        }
        for (double e : ed) res.e_distant += e;  // fixed order
    }
    res.t_distant_ms = msSince(t_dist);

    // ---- 6. couplings and PNO overlaps S_pq = c_p^T Spao(idx_p, idx_q) c_q ----
    const auto t_ov = std::chrono::steady_clock::now();
    std::vector<std::vector<Coupling>> cpl(npairs);
    std::vector<std::vector<int>> slot_q(npairs);  // q of each overlap slot
    long ncpl = 0;
    double ov_mb = 0.0;
    for (int p = 0; p < npairs; ++p) {
        const int i = P[p].i, j = P[p].j;
        std::map<int, int> slot;
        auto add = [&](int q, double f, bool tr) {
            if (q < 0) return;
            auto it = slot.find(q);
            int s;
            if (it == slot.end()) {
                s = (int)slot_q[p].size();
                slot[q] = s;
                slot_q[p].push_back(q);
                ov_mb += (double)P[p].eps.size() * P[q].eps.size() * 8.0 / 1e6;
            } else {
                s = it->second;
            }
            cpl[p].push_back({q, s, f, tr});
        };
        for (int k = 0; k < o; ++k) {
            if (k != i && std::abs(Focc(i, k)) > opt.fock_cut) add(pair_index[(size_t)k * o + j], Focc(i, k), k > j);
            if (k != j && std::abs(Focc(k, j)) > opt.fock_cut) add(pair_index[(size_t)i * o + k], Focc(k, j), i > k);
        }
        ncpl += (long)slot_q[p].size();
    }
    res.n_couplings = ncpl;
    res.memory_mb = (local_fit ? res.factors_mb : b_mb) + pno_mb + ov_mb;
    if (pno_mb + ov_mb > opt.max_memory_mb) {
        CurcumaLogger::error(fmt::format("LMP2: PNO data ({:.0f} MB) and pair overlaps ({:.0f} MB, {} couplings) exceed "
                                         "-qm.lmp2_max_memory_mb {:g}",
                                         pno_mb, ov_mb, ncpl, opt.max_memory_mb));
        return false;
    }
    std::vector<std::vector<MatX>> Sov(npairs);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int p = 0; p < npairs; ++p) {
        const PairData& pp = P[p];
        const int mao = (int)pp.idx.size();
        // Pp = Spao(:, idx_p) c_p (n x npno_p), then S_pq = Pp(idx_q, :)^T c_q
        MatX G(n, mao);
        for (int r = 0; r < mao; ++r) G.col(r) = Spao.col(pp.idx[r]);
        const MatX Pp = G * pp.c;
        Sov[p].resize(slot_q[p].size());
        for (size_t s = 0; s < slot_q[p].size(); ++s) {
            const PairData& pq = P[slot_q[p][s]];
            MatX Pr(pq.idx.size(), Pp.cols());
            for (size_t r = 0; r < pq.idx.size(); ++r) Pr.row(r) = Pp.row(pq.idx[r]);
            Sov[p][s] = Pr.transpose() * pq.c;
        }
    }
    res.t_overlap_ms = msSince(t_ov);

    // ---- 7. LMP2 iterations (lmp2_common.solve_lmp2 in the PNO spaces) ----
    const auto t_it = std::chrono::steady_clock::now();
    std::vector<size_t> off(npairs + 1, 0);
    for (int p = 0; p < npairs; ++p) off[p + 1] = off[p] + (size_t)P[p].K.size();
    const size_t ntot = off[npairs];
    VecX t(ntot), R(ntot);
    for (int p = 0; p < npairs; ++p)
        Eigen::Map<MatX>(t.data() + off[p], P[p].K.rows(), P[p].K.cols()) = -P[p].K.cwiseQuotient(P[p].den);
    auto energy = [&](const VecX& x) {
        std::vector<double> e(npairs);
#pragma omp parallel for schedule(static) num_threads(threads)
        for (int p = 0; p < npairs; ++p) {
            const Eigen::Map<const MatX> T(x.data() + off[p], P[p].K.rows(), P[p].K.cols());
            e[p] = pairEnergy(T, P[p].K, P[p].i == P[p].j);
        }
        double s = 0.0;
        for (double v_ : e) s += v_;  // fixed order
        return s;
    };
    double e_old = energy(t);
    Diis diis(opt.diis);
    bool converged = false;
    int it = 0;
    std::vector<double> rmax_p(npairs);
    for (it = 1; it <= opt.max_iter; ++it) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (int p = 0; p < npairs; ++p) {
            const PairData& pd = P[p];
            const Eigen::Index np = pd.K.rows();
            const Eigen::Map<const MatX> Tp(t.data() + off[p], np, np);
            Eigen::Map<MatX> Rp(R.data() + off[p], np, np);
            Rp = pd.K + pd.den.cwiseProduct(Tp);
            for (const Coupling& c : cpl[p]) {
                const Eigen::Index nq = P[c.q].K.rows();
                if (np == 0 || nq == 0) continue;
                const Eigen::Map<const MatX> Tq(t.data() + off[c.q], nq, nq);
                const MatX& Sq = Sov[p][c.slot];
                const MatX X = c.transpose ? MatX(Sq * Tq.transpose()) : MatX(Sq * Tq);
                Rp.noalias() -= c.f * (X * Sq.transpose());
            }
            rmax_p[p] = np ? Rp.cwiseAbs().maxCoeff() : 0.0;
        }
        double rmax = 0.0;
        for (double r : rmax_p) rmax = std::max(rmax, r);
        VecX step(ntot);
        for (int p = 0; p < npairs; ++p) {
            const Eigen::Index np = P[p].K.rows();
            Eigen::Map<MatX>(step.data() + off[p], np, np) =
                Eigen::Map<const MatX>(R.data() + off[p], np, np).cwiseQuotient(P[p].den);
        }
        t = diis.update(VecX(t - step), step);
        const double e = energy(t);
        if (!std::isfinite(e)) {
            CurcumaLogger::error("LMP2: non-finite energy in the iterations");
            return false;
        }
        res.last_de = e - e_old;
        res.max_residual = rmax;
        if (CurcumaLogger::get_verbosity() >= 1)  // iterative progress from level 1 on (was 3)
            CurcumaLogger::result(fmt::format("LMP2 it {:3d}  E = {:.12f}  dE = {:.2e}  max|R| = {:.2e}", it, e, e - e_old, rmax));
        e_old = e;
        if (std::abs(res.last_de) < opt.e_threshold && rmax < opt.r_threshold) {
            converged = true;
            break;
        }
    }
    res.iterations = std::min(it, opt.max_iter);
    res.e_iter = e_old;
    res.t_iter_ms = msSince(t_it);
    if (!converged) {
        CurcumaLogger::error(fmt::format("LMP2: amplitude equations not converged in {} iterations (last dE {:.1e}, max|R| {:.1e})",
                                         opt.max_iter, res.last_de, res.max_residual));
        return false;
    }
    res.info = fmt::format("{}; {} threads", ri_src, threads);
    return true;
}
