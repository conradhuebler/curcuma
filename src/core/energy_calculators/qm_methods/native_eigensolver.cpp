/*
 * <Native symmetric eigensolver — implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated. GPL-3.0.
 *
 * Householder tridiagonalization (tred2) + implicit-shift QL with eigenvector
 * accumulation (tql2). Classic, dependency-free symmetric eigensolver (EISPACK /
 * Numerical Recipes lineage), 0-indexed and expressed on Eigen storage. Used as the
 * `-eigensolver native` alternative to MKL dsyevd; see native_eigensolver.h. The
 * tridiagonal QL step is the drop-in point for a future Cuppen divide-and-conquer
 * solver (same d/e in, eigenpairs out).
 */

#include "native_eigensolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace curcuma::eigsolver {

namespace {

// sqrt(a^2 + b^2) without destructive under/overflow (Numerical Recipes pythag).
inline double pythag(double a, double b)
{
    const double aa = std::fabs(a), ab = std::fabs(b);
    if (aa > ab) { const double r = ab / aa; return aa * std::sqrt(1.0 + r * r); }
    if (ab == 0.0) return 0.0;
    const double r = aa / ab; return ab * std::sqrt(1.0 + r * r);
}

inline double sign(double a, double b) { return (b >= 0.0) ? std::fabs(a) : -std::fabs(a); }

/* ------------------------------------------------------------------ *
 *  tred2 — Householder reduction of a real symmetric matrix to       *
 *  tridiagonal form. On entry `z` holds the symmetric matrix; on     *
 *  exit `z` is the orthogonal transform Q (A = Q·T·Qᵀ), `d` the       *
 *  diagonal of T, `e` its sub-diagonal (e[0]=0).                     *
 * ------------------------------------------------------------------ */
void tred2(Eigen::MatrixXd& z, Eigen::VectorXd& d, Eigen::VectorXd& e)
{
    const int n = static_cast<int>(z.rows());
    for (int i = n - 1; i >= 1; --i) {
        const int l = i - 1;
        double h = 0.0, scale = 0.0;
        if (l > 0) {
            for (int k = 0; k <= l; ++k) scale += std::fabs(z(i, k));
            if (scale == 0.0) {
                e(i) = z(i, l);
            } else {
                for (int k = 0; k <= l; ++k) { z(i, k) /= scale; h += z(i, k) * z(i, k); }
                double f = z(i, l);
                const double g = (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
                e(i) = scale * g;
                h -= f * g;
                z(i, l) = f - g;
                f = 0.0;
                for (int j = 0; j <= l; ++j) {
                    z(j, i) = z(i, j) / h;
                    double gg = 0.0;
                    for (int k = 0; k <= j; ++k)     gg += z(j, k) * z(i, k);
                    for (int k = j + 1; k <= l; ++k) gg += z(k, j) * z(i, k);
                    e(j) = gg / h;
                    f += e(j) * z(i, j);
                }
                const double hh = f / (h + h);
                for (int j = 0; j <= l; ++j) {
                    f = z(i, j);
                    const double gg = e(j) - hh * f;
                    e(j) = gg;
                    for (int k = 0; k <= j; ++k)
                        z(j, k) -= (f * e(k) + gg * z(i, k));
                }
            }
        } else {
            e(i) = z(i, l);
        }
        d(i) = h;
    }
    d(0) = 0.0;
    e(0) = 0.0;
    // Accumulate the transformation matrices.
    for (int i = 0; i < n; ++i) {
        const int l = i - 1;
        if (d(i) != 0.0) {
            for (int j = 0; j <= l; ++j) {
                double g = 0.0;
                for (int k = 0; k <= l; ++k) g += z(i, k) * z(k, j);
                for (int k = 0; k <= l; ++k) z(k, j) -= g * z(k, i);
            }
        }
        d(i) = z(i, i);
        z(i, i) = 1.0;
        for (int j = 0; j <= l; ++j) { z(j, i) = 0.0; z(i, j) = 0.0; }
    }
}

/* ------------------------------------------------------------------ *
 *  tql2 — eigenvalues/vectors of a symmetric tridiagonal matrix by   *
 *  the QL algorithm with implicit shifts. `d` holds the diagonal     *
 *  (→ eigenvalues), `e` the sub-diagonal, `z` the accumulated        *
 *  transform from tred2 (→ eigenvectors, columns). Returns false if  *
 *  an eigenvalue needs > 50 iterations to isolate.                   *
 * ------------------------------------------------------------------ */
bool tql2(Eigen::MatrixXd& z, Eigen::VectorXd& d, Eigen::VectorXd& e)
{
    const int n = static_cast<int>(z.rows());
    for (int i = 1; i < n; ++i) e(i - 1) = e(i);
    e(n - 1) = 0.0;

    for (int l = 0; l < n; ++l) {
        int iter = 0;
        int m = 0;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::fabs(d(m)) + std::fabs(d(m + 1));
                if (std::fabs(e(m)) <= std::numeric_limits<double>::epsilon() * dd) break;
            }
            if (m != l) {
                if (iter++ == 50) return false;
                double g = (d(l + 1) - d(l)) / (2.0 * e(l));
                double r = pythag(g, 1.0);
                g = d(m) - d(l) + e(l) / (g + sign(r, g));
                double s = 1.0, c = 1.0, p = 0.0;
                int i = m - 1;
                for (; i >= l; --i) {
                    double f = s * e(i);
                    const double b = c * e(i);
                    r = pythag(f, g);
                    e(i + 1) = r;
                    if (r == 0.0) { d(i + 1) -= p; e(m) = 0.0; break; }
                    s = f / r;
                    c = g / r;
                    g = d(i + 1) - p;
                    r = (d(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d(i + 1) = g + p;
                    g = c * r - b;
                    for (int k = 0; k < n; ++k) {
                        f = z(k, i + 1);
                        z(k, i + 1) = s * z(k, i) + c * f;
                        z(k, i)     = c * z(k, i) - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d(l) -= p;
                e(l) = g;
                e(m) = 0.0;
            }
        } while (m != l);
    }
    return true;
}

/* ================================================================== *
 *  Cuppen divide-and-conquer for the symmetric tridiagonal problem.   *
 *  (Cuppen 1981; Dongarra–Sorensen 1987; Gu–Eisenstat 1994.)          *
 * ================================================================== */

// Ascending sort permutation of a vector.
std::vector<int> argsortAsc(const Eigen::VectorXd& v)
{
    std::vector<int> p(static_cast<size_t>(v.size()));
    for (size_t i = 0; i < p.size(); ++i) p[i] = static_cast<int>(i);
    std::sort(p.begin(), p.end(), [&](int a, int b) { return v(a) < v(b); });
    return p;
}

// Base case: full tridiagonal (diag, off) via QL with eigenvector accumulation.
// `off(k)` connects rows (k, k+1). Returns ascending eigenvalues + eigenvectors (columns).
bool solveTriQL(const Eigen::VectorXd& diag, const Eigen::VectorXd& off,
                Eigen::VectorXd& eval, Eigen::MatrixXd& evec)
{
    const int n = static_cast<int>(diag.size());
    evec = Eigen::MatrixXd::Identity(n, n);
    Eigen::VectorXd d = diag;
    Eigen::VectorXd e(n);
    e(0) = 0.0;
    for (int k = 0; k < n - 1; ++k) e(k + 1) = off(k);
    if (!tql2(evec, d, e)) return false;
    const std::vector<int> idx = argsortAsc(d);
    eval.resize(n);
    Eigen::MatrixXd V(n, n);
    for (int j = 0; j < n; ++j) { eval(j) = d(idx[j]); V.col(j) = evec.col(idx[j]); }
    evec = std::move(V);
    return true;
}

// Eigenproblem of diag(D) + rho·z·zᵀ (rho>0, D strictly ascending, z all non-negligible).
// Roots of the secular equation 1 + rho·Σ z_i²/(D_i−λ)=0 lie one per gap (D_j, D_{j+1})
// plus one in (D_{k-1}, D_{k-1}+rho·‖z‖²); found by bisection. Eigenvectors via the
// Gu–Eisenstat (Löwner) reconstructed weights for numerical orthogonality.
void secularSolve(const Eigen::VectorXd& D, const Eigen::VectorXd& z, double rho,
                  Eigen::VectorXd& lam, Eigen::MatrixXd& W)
{
    const int k = static_cast<int>(D.size());
    lam.resize(k);
    const double znorm2 = z.squaredNorm();
    const double eps = std::numeric_limits<double>::epsilon();

    // dml(i,j) = D_i − λ_j, computed WITHOUT cancellation: anchor each root at its nearer
    // pole (dlaed4 trick) and solve for the shift τ, so D_i − λ = (D_i − anchor) − τ stays
    // accurate even when the root sits next to a pole. Direct bisection on λ loses the
    // relative precision of the small denominators and corrupts the eigenvectors.
    Eigen::MatrixXd dml(k, k);
    for (int j = 0; j < k; ++j) {
        const double left  = D(j);
        const double right = (j < k - 1) ? D(j + 1) : (D(k - 1) + rho * znorm2);
        // f increasing on (left,right) from −∞ to +∞; pick the nearer pole as origin.
        const double mid = 0.5 * (left + right);
        double fmid = 1.0;
        for (int i = 0; i < k; ++i) fmid += rho * z(i) * z(i) / (D(i) - mid);
        const double origin = (j < k - 1 && fmid < 0.0) ? D(j + 1) : D(j);

        // Solve h(τ)=1+ρ·Σ z_i²/((D_i−origin)−τ)=0 for τ = λ−origin in (left−origin, right−origin).
        // Bisect τ to FULL machine precision (relative to τ itself, not |origin|): the shift τ
        // is exactly the small distance to the nearer pole, so its relative accuracy is what
        // the eigenvector denominators dml need. Natural fp convergence (tau pinned between
        // adjacent floats) terminates the loop. Claude Generated.
        double a = left - origin, b = right - origin;
        double tau = 0.5 * (a + b);
        for (int it = 0; it < 200; ++it) {
            const double t = 0.5 * (a + b);
            if (t <= a || t >= b) break;   // fp resolution reached
            double h = 1.0;
            for (int i = 0; i < k; ++i) h += rho * z(i) * z(i) / ((D(i) - origin) - t);
            if (h < 0.0) a = t; else b = t;
            tau = t;
        }
        (void)eps;
        lam(j) = origin + tau;
        for (int i = 0; i < k; ++i) dml(i, j) = (D(i) - origin) - tau;  // D_i − λ_j (accurate)
    }

    // Löwner-reconstructed weights: z̃_i² = ∏_j (λ_j−D_i) / ∏_{j≠i}(D_j−D_i), each ratio O(1).
    // Gives eigenvectors that are orthogonal even for tight (post-deflation) gaps. Uses the
    // accurate dml (λ_j−D_i = −dml(i,j)).
    Eigen::VectorXd zt(k);
    for (int i = 0; i < k; ++i) {
        double prod = -dml(i, i);                              // λ_i − D_i
        for (int j = 0; j < k; ++j)
            if (j != i) prod *= (-dml(i, j)) / (D(j) - D(i));  // (λ_j−D_i)/(D_j−D_i)
        zt(i) = std::sqrt(std::max(prod, 0.0)) * (z(i) >= 0.0 ? 1.0 : -1.0);
    }
    W.resize(k, k);
    for (int j = 0; j < k; ++j) {
        for (int i = 0; i < k; ++i) W(i, j) = zt(i) / dml(i, j);  // z̃_i/(D_i−λ_j)
        const double nrm = W.col(j).norm();
        if (nrm > 0.0) W.col(j) /= nrm;
    }
}

// Eigenproblem of diag(D_in) + rho_in·z_in·z_inᵀ with deflation. Returns eigenvalues
// `eval` and eigenvectors `W` (columns), eval(j) ↔ W.col(j), in the input basis/order.
void rank1Eigen(const Eigen::VectorXd& D_in, const Eigen::VectorXd& z_in, double rho_in,
                Eigen::VectorXd& eval, Eigen::MatrixXd& W)
{
    const int n = static_cast<int>(D_in.size());
    eval.resize(n);
    W = Eigen::MatrixXd::Zero(n, n);

    const double zn = z_in.norm();
    if (zn < std::numeric_limits<double>::min()) {   // no coupling
        eval = D_in;
        W.setIdentity();
        return;
    }
    // Normalise z (‖z‖=1) and fold the sign of rho (solve for −D if rho<0, negate at end).
    Eigen::VectorXd z = z_in / zn;
    double rho = rho_in * zn * zn;
    const bool flip = (rho < 0.0);
    Eigen::VectorXd D = flip ? (-D_in).eval() : D_in;
    if (flip) rho = -rho;

    // Sort D ascending (work in the sorted basis; row i ↔ original index perm[i]).
    const std::vector<int> perm = argsortAsc(D);
    Eigen::VectorXd Ds(n), zs(n);
    for (int i = 0; i < n; ++i) { Ds(i) = D(perm[i]); zs(i) = z(perm[i]); }

    const double tol = 8.0 * std::numeric_limits<double>::epsilon()
                       * (Ds.cwiseAbs().maxCoeff() + std::fabs(rho));

    // Deflation. Givens (plane p,i; c,s) recorded for degenerate-D merges; `deflated`
    // marks indices whose eigenpair is fixed (a unit vector e_i in the rotated basis).
    struct Givens { int p, i; double c, s; };
    std::vector<Givens> givens;
    std::vector<char> deflated(n, 0);
    int jprev = -1;
    for (int i = 0; i < n; ++i) {
        if (jprev < 0) { jprev = i; continue; }
        if (std::fabs(Ds(i) - Ds(jprev)) <= tol) {
            // Degenerate diagonal: rotate to pour zs(i) into zs(jprev), deflating i.
            const double r = pythag(zs(jprev), zs(i));
            const double c = (r > 0.0) ? zs(jprev) / r : 1.0;
            const double s = (r > 0.0) ? zs(i) / r : 0.0;
            givens.push_back({jprev, i, c, s});
            zs(jprev) = r;
            zs(i) = 0.0;
            deflated[i] = 1;
        } else {
            jprev = i;
        }
    }
    for (int i = 0; i < n; ++i)
        if (!deflated[i] && std::fabs(zs(i)) <= tol) deflated[i] = 1;  // negligible weight

    // Secular set (non-deflated).
    std::vector<int> sec;
    for (int i = 0; i < n; ++i) if (!deflated[i]) sec.push_back(i);

    Eigen::VectorXd lam_sec;
    Eigen::MatrixXd W_sec;
    if (!sec.empty()) {
        const int ksz = static_cast<int>(sec.size());
        Eigen::VectorXd delta(ksz), zeta(ksz);
        for (int t = 0; t < ksz; ++t) { delta(t) = Ds(sec[t]); zeta(t) = zs(sec[t]); }
        secularSolve(delta, zeta, rho, lam_sec, W_sec);
    }

    // Assemble eigenpairs in the rotated sorted basis (M'): deflated → e_i, secular → W_sec.
    Eigen::MatrixXd Wp = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd evalp(n);
    int slot = 0;
    for (int i = 0; i < n; ++i) {
        if (deflated[i]) { Wp(i, slot) = 1.0; evalp(slot) = Ds(i); ++slot; }
    }
    for (int t = 0; t < static_cast<int>(sec.size()); ++t) {
        for (int u = 0; u < static_cast<int>(sec.size()); ++u)
            Wp(sec[u], slot) = W_sec(u, t);
        evalp(slot) = lam_sec(t);
        ++slot;
    }

    // Undo the Givens rotations (apply Gᵀ in reverse) → eigenvectors in the sorted basis.
    for (int g = static_cast<int>(givens.size()) - 1; g >= 0; --g) {
        const Givens& G = givens[g];
        for (int col = 0; col < n; ++col) {
            const double a = Wp(G.p, col), b = Wp(G.i, col);
            Wp(G.p, col) = G.c * a - G.s * b;
            Wp(G.i, col) = G.s * a + G.c * b;
        }
    }

    // Scatter back to the original basis (row perm[i] ← sorted row i) and unflip eigenvalues.
    for (int i = 0; i < n; ++i) W.row(perm[i]) = Wp.row(i);
    eval = flip ? (-evalp).eval() : evalp;
}

// Recursive Cuppen divide-and-conquer for the tridiagonal (diag d, off-diagonal e).
bool tridiagDC(const Eigen::VectorXd& d, const Eigen::VectorXd& e,
               Eigen::VectorXd& eval, Eigen::MatrixXd& evec)
{
    const int n = static_cast<int>(d.size());
    if (n == 1) { eval = d; evec = Eigen::MatrixXd::Constant(1, 1, 1.0); return true; }
    constexpr int kCutoff = 32;
    if (n <= kCutoff) return solveTriQL(d, e, eval, evec);

    const int m = n / 2;
    const double rho = e(m - 1);                 // tearing off-diagonal

    Eigen::VectorXd d1 = d.head(m);              d1(m - 1) -= rho;
    Eigen::VectorXd d2 = d.tail(n - m);          d2(0)     -= rho;
    Eigen::VectorXd e1 = e.head(m - 1);
    Eigen::VectorXd e2 = e.tail(n - m - 1);

    Eigen::VectorXd ev1, ev2;
    Eigen::MatrixXd V1, V2;
    if (!tridiagDC(d1, e1, ev1, V1)) return false;
    if (!tridiagDC(d2, e2, ev2, V2)) return false;

    // Merge: M = diag([ev1; ev2]) + rho·z·zᵀ, z = [last row of V1; first row of V2].
    Eigen::VectorXd D(n);
    D.head(m) = ev1;  D.tail(n - m) = ev2;
    Eigen::VectorXd z(n);
    z.head(m)     = V1.row(m - 1).transpose();
    z.tail(n - m) = V2.row(0).transpose();

    Eigen::VectorXd lam;
    Eigen::MatrixXd Wm;
    rank1Eigen(D, z, rho, lam, Wm);

    eval = lam;
    evec.resize(n, n);                            // evec = blockdiag(V1,V2)·Wm
    evec.topRows(m)        = V1 * Wm.topRows(m);
    evec.bottomRows(n - m) = V2 * Wm.bottomRows(n - m);
    return true;
}

} // namespace

bool solveSymmetric(const Eigen::MatrixXd& A,
                    Eigen::VectorXd& evals,
                    Eigen::MatrixXd& evecs)
{
    const int n = static_cast<int>(A.rows());
    if (n == 0 || A.cols() != n) return false;

    // 1. Householder tridiagonalization: A = Q·T·Qᵀ. Symmetrise a copy first to guard
    //    against tiny input asymmetries (tred2 works in place on this copy → Q).
    Eigen::MatrixXd Q = 0.5 * (A + A.transpose());
    Eigen::VectorXd d(n), efull(n);
    tred2(Q, d, efull);          // d = diagonal of T; efull(i) = subdiag(i-1,i), efull(0)=0

    if (n == 1) { evals = d; evecs = Q; return true; }

    // 2. Tridiagonal eigenproblem T = V·Λ·Vᵀ by Cuppen divide-and-conquer.
    Eigen::VectorXd eoff(n - 1);
    for (int k = 0; k < n - 1; ++k) eoff(k) = efull(k + 1);   // off(k) connects (k,k+1)
    Eigen::VectorXd ev;
    Eigen::MatrixXd V;
    if (!tridiagDC(d, eoff, ev, V)) return false;

    // 3. Back-transform eigenvectors C = Q·V and sort eigenpairs ascending.
    const Eigen::MatrixXd QV = Q * V;
    const std::vector<int> idx = argsortAsc(ev);
    evals.resize(n);
    evecs.resize(n, n);
    for (int j = 0; j < n; ++j) { evals(j) = ev(idx[j]); evecs.col(j) = QV.col(idx[j]); }
    return true;
}

} // namespace curcuma::eigsolver
