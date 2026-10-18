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

#include "qm_engine.h"
#include "qm_mp2.h"  // frozenCoreOrbitals
#include "qm_ri.h"

#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>

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
    if (opt.distant_pairs != "dipole" && opt.distant_pairs != "drop") {
        CurcumaLogger::error("LMP2: unknown -qm.lmp2_distant_pairs '" + opt.distant_pairs + "' (use dipole | drop)");
        return false;
    }
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
    const auto t_loc = std::chrono::steady_clock::now();
    std::vector<MatX> Dip(3);  // active-basis dipole matrices (Bohr)
    {
        Matrix Sc, Xc, Yc, Zc;
        cartesianDipoleIntegrals(m_gto_basis, Sc, Xc, Yc, Zc);
        const Matrix Sa = qmint::applySphericalTransform(Sc, m_Q);
        const double dS = (Sa - m_S).cwiseAbs().maxCoeff();
        if (dS > 1e-8) {  // self-check of the integral code against the engine's overlap
            CurcumaLogger::error(fmt::format("LMP2: dipole integral self-check failed (overlap differs by {:.1e})", dS));
            return false;
        }
        Dip[0] = toEigen(qmint::applySphericalTransform(Xc, m_Q));
        Dip[1] = toEigen(qmint::applySphericalTransform(Yc, m_Q));
        Dip[2] = toEigen(qmint::applySphericalTransform(Zc, m_Q));
    }
    MatX L;
    if (!choleskyOrbitals(Cval, S, L)) {
        CurcumaLogger::error("LMP2: Cholesky start orbitals failed (occupied space rank-deficient?)");
        return false;
    }
    {
        std::vector<MatX> M;
        if (opt.localizer == "boys") {
            for (int c = 0; c < 3; ++c) M.push_back(L.transpose() * Dip[c] * L);
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
        double f = 0.0;
        const int sweeps = jacobiLocalize(L, M, threads, 1e-11, 500, f);
        res.localization_sweeps = sweeps;
        if (sweeps < 0 && CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::warn(fmt::format("LMP2: localisation not converged after {} sweeps (functional {:.10f}); the "
                                            "energy with truncation depends slightly on the orbitals",
                                            -sweeps, f));
    }
    // deterministic order (centroid x, y, z rounded to 1e-6 Bohr) and sign (largest coefficient positive)
    MatX cen(o, 3);
    for (int i = 0; i < o; ++i)
        for (int c = 0; c < 3; ++c) cen(i, c) = L.col(i).dot(Dip[c] * L.col(i));
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
        for (int k = 0; k < o; ++k) {
            L2.col(k) = L.col(order[k]);
            cen2.row(k) = cen.row(order[k]);
            Eigen::Index imax;
            L2.col(k).cwiseAbs().maxCoeff(&imax);
            if (L2(imax, k) < 0) L2.col(k) *= -1.0;
        }
        L.swap(L2);
        cen.swap(cen2);
    }
    const MatX U = Cval.transpose() * S * L;  // L = Cval U
    const MatX Focc = U.transpose() * eps.segment(n_frozen, o).asDiagonal() * U;
    res.t_localize_ms = msSince(t_loc);

    // ---- 2. domains (bp_domain + extend_domain) and pairs ----
    const auto t_pairs = std::chrono::steady_clock::now();
    std::vector<std::vector<int>> dom(o);
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < o; ++i) {
        std::vector<int> bp = bpDomain(L.col(i), S, atom_aos, opt.bp_threshold);
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

    // ---- 3. RI tensor and the half-transformed factors B_{i,r}^P ----
    std::string ri_src, why;
    double setup_ms = 0.0;
    std::shared_ptr<qmri::RIJK> ri = correlationRI(file, ri_src, setup_ms, why);
    if (!ri) {
        CurcumaLogger::error("LMP2: density fitting setup with " + file + " not possible: " + why);
        return false;
    }
    res.t_ri_setup_ms = setup_ms;
    const int nrows = (int)ri->tensor().rows();
    const double b_mb = (double)ncol * nrows * 8.0 / 1.0e6;
    if (2.0 * b_mb > opt.max_memory_mb) {
        CurcumaLogger::error(fmt::format("LMP2: the factors B_(i,r)^P need 2 x {:.1f} MB ({} columns over the extended domains x {} "
                                         "auxiliary functions, plus one copy during the fit), above -qm.lmp2_max_memory_mb {:g}",
                                         b_mb, ncol, nrows, opt.max_memory_mb));
        return false;
    }
    const auto t_tr = std::chrono::steady_clock::now();
    std::vector<MatX> Bf(o);  // Bf[i]: naux x |ext_ao[i]|, column r = B_{i, ext_ao[i][r]}^P
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
        const int naux = (int)Xr.rows();
        res.n_aux = naux;
#pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (int i = 0; i < o; ++i) Bf[i] = Xr.block(0, (Eigen::Index)ext_off[i], naux, (Eigen::Index)ext_ao[i].size());
    }
    res.t_transform_ms = msSince(t_tr);
    const int naux = res.n_aux;
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("LMP2: {}; B factors {:.0f} MB ({} columns)", ri_src, b_mb, ncol));

    // ---- 4. pair spaces, K in the PAO space, PNOs (lmp2_common.pair_space, pair_ints, make_pnos) ----
    std::vector<PairData> P(npairs);
    std::vector<char> pair_ok(npairs, 1);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int p = 0; p < npairs; ++p) {
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
        PaoSpace sp;
        if (!paoSpace(pd.idx, Spao, Fpao, sp)) {
            pair_ok[p] = 0;
            continue;
        }
        const int mao = (int)pd.idx.size(), m = (int)sp.eps.size();
        pd.m_pao = m;
        // B_ia^P = sum_r Xc_ra B_{i,r}^P over the pair's PAOs
        auto half = [&](int k) {
            MatX G(naux, mao);
            for (int r = 0; r < mao; ++r) G.col(r) = Bf[k].col(ext_pos[k][pd.idx[r]]);
            return MatX(G * sp.Xc);
        };
        const MatX Bi = half(i);
        const MatX K = (i == j) ? MatX(Bi.transpose() * Bi) : MatX(Bi.transpose() * half(j));
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
    }
    for (int p = 0; p < npairs; ++p)
        if (!pair_ok[p]) {
            CurcumaLogger::error("LMP2: eigensolver failed in a pair space");
            return false;
        }
    std::vector<MatX>().swap(Bf);
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
    res.t_pairs_ms = msSince(t_pairs) - res.t_transform_ms - res.t_ri_setup_ms;

    // ---- 5. distant pairs: dipole estimate over the PAO spaces of [i] and [j] ----
    const auto t_dist = std::chrono::steady_clock::now();
    if (opt.distant_pairs == "dipole" && !distant.empty()) {
        std::vector<char> need(o, 0);
        for (const auto& pr : distant) need[pr.first] = need[pr.second] = 1;
        std::vector<VecX> oe(o);
        std::vector<MatX> mu(o);  // 3 x m_i: <i|r|a> over the pseudo-canonical PAO space of [i]
        std::vector<MatX> W(3);   // PAO^T r L: n x o
        for (int c = 0; c < 3; ++c) W[c] = Cpao.transpose() * (Dip[c] * L);
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
    res.memory_mb = b_mb + pno_mb + ov_mb;
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
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("LMP2 it {:3d}  E = {:.12f}  dE = {:.2e}  max|R| = {:.2e}", it, e, e - e_old, rmax));
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
