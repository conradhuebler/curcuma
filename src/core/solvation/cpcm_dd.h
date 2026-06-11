/*
 * <Domain-decomposition ddCOSMO core for the native CPCM solvation model>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (June 2026): one-to-one C++ port of the domain-decomposition
 * ddCOSMO solver from tblite (external/tblite/src/tblite/solvation/cpcm_dd.f90),
 * which underlies tblite's CPCM continuum-solvation model. The math (Lipparini
 * et al., domain-decomposition COSMO) is unchanged; only the host language and
 * the 1-based -> 0-based, column-major bookkeeping differ.
 *
 * NOTE: all 2D state uses Eigen's default COLUMN-major layout so that a column
 * (e.g. sigma.col(iat), basis.col(ig)) is contiguous, matching the Fortran
 * arrays self%sigma(:, iat) etc. (Curcuma's global `Matrix` is row-major and
 * is therefore NOT used here.)
 */

#pragma once

#include "src/core/solvation/lebedev_grid.h"

#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace Curcuma {
namespace Solvation {

/// ddCOSMO control parameters (tblite domain_decomposition_input).
struct DomainDecompositionInput {
    int lmax = 6;        ///< max angular momentum of the spherical-harmonic basis
    double conv = 1e-8;  ///< iterative-solver threshold
    double eta = 2.0;    ///< switch-function regularisation width
};

/**
 * @brief Domain-decomposition ddCOSMO calculator (tblite domain_decomposition).
 *
 * @ref init builds the geometry-independent basis (spherical harmonics on the
 * Lebedev grid, normalisation factors). @ref update rebuilds the cavity for a
 * new geometry. @ref solveCosmoDirect / @ref solveCosmoAdjoint solve the forward
 * and adjoint ddCOSMO linear systems; @ref getDeriv assembles the gradient.
 */
class DomainDecomposition {
public:
    using Mat = Eigen::MatrixXd; // column-major

    /// Initialise basis + factors. `grid` carries unit-sphere points and raw
    /// Lebedev weights (sum = 1); they are scaled by 4*pi internally.
    void init(const DomainDecompositionInput& input,
              const std::vector<double>& rvdw,
              const std::vector<LebedevPoint>& grid);

    /// Rebuild the cavity (neighbour list, characteristic function, points)
    /// for the geometry `xyz` (shape 3 x nat, atomic units / Bohr).
    void update(const Mat& xyz);

    // --- ddCOSMO matrix-vector products (x, y shape nylm x nat) ---
    void lx(const Mat& x, Mat& y) const;      ///< y = L x (off-diagonal blocks, negated)
    void lstarx(const Mat& x, Mat& y) const;  ///< y = L* x (adjoint off-diagonal, negated)
    void ldm1x(const Mat& x, Mat& y) const;    ///< y = D^-1 x (inverse diagonal block)
    double hnorm(const Mat& x) const;          ///< H^-1/2 rms norm of the increment

    // --- linear solvers ---
    /// Solve L sigma = G. If `cart`, the RHS is assembled from the cavity
    /// potential `phi` (length ncav); otherwise `glm` is used directly.
    void solveCosmoDirect(bool cart, const std::vector<double>& phi,
                          const Mat& glm, Mat& sigma, bool restart) const;
    /// Solve L* s = psi (adjoint). `accuracy` < 0 uses the default conv.
    void solveCosmoAdjoint(const Mat& psi, Mat& sigma, bool restart,
                           double accuracy = -1.0) const;

    // --- gradient ---
    void getDeriv(double keps, const std::vector<double>& phi,
                  const Mat& sigma, const Mat& s, Mat& gradient) const;
    void getZeta(double keps, const Mat& s, std::vector<double>& zeta) const;

    // --- sizes (public so the CPCM wrapper can allocate) ---
    int nat = 0;     ///< number of spheres/atoms
    int ngrid = 0;   ///< Lebedev points per sphere
    int nylm = 0;    ///< number of spherical-harmonic basis functions = (lmax+1)^2
    int ncav = 0;    ///< number of active cavity points
    int lmax = 6;
    double conv = 1e-8;
    double eta = 2.0;

    const Mat& cavityPoints() const { return m_ccav; } ///< 3 x ncav active cavity points
    const Mat& atomCoords() const { return m_xyz; }    ///< 3 x nat atom coordinates (Bohr)
    /// Owning atom (0-based) of each active cavity point, in ccav order.
    const std::vector<int>& cavityAtom() const { return m_cavAtom; }
    const Eigen::Map<const Eigen::VectorXd> faclVec() const {
        return Eigen::Map<const Eigen::VectorXd>(m_facl.data(), nylm);
    }

private:
    static constexpr double se = -1.0;   ///< interior shift for the switch function
    static constexpr int nngmax = 100;   ///< max neighbours per sphere (CSR capacity)
    static constexpr int ndiis = 25;     ///< DIIS history length

    // geometry-independent state
    std::vector<double> m_rvdw;  // [nat]
    std::vector<double> m_w;     // [ngrid] Lebedev weights * 4pi
    Mat m_grid;                  // 3 x ngrid (unit sphere)
    Mat m_basis;                 // nylm x ngrid (Y_lm at grid points)
    std::vector<double> m_fact;  // [2*lmax+1] factorials k!
    std::vector<double> m_facl;  // [nylm] (2l+1)/(4pi)
    std::vector<double> m_facs;  // [nylm] real-harmonic normalisation
    bool m_grad = true;

    // geometry-dependent state
    Mat m_xyz;                   // 3 x nat
    std::vector<int> m_inl;      // [nat+1] CSR offsets (0-based)
    std::vector<int> m_nl;       // [nat*nngmax] neighbour atom indices (0-based)
    Mat m_fi;                    // ngrid x nat
    Mat m_ui;                    // ngrid x nat
    std::vector<double> m_zi;    // 3*ngrid*nat (switch-region gradient direction)
    Mat m_ccav;                  // 3 x ncav
    std::vector<int> m_cavAtom;  // [ncav] owning atom of each cavity point

    // cavity builders
    void mknnl();
    void mkfiui();
    void mkccav();

    // spherical-harmonic helpers (raw pointers, length nylm / lmax+1)
    void ylmbas(const double x[3], double* basloc, double* vplm,
                double* vcos, double* vsin) const;
    void dbasis(const double x[3], double* basloc, double* dbsloc /*3*nylm*/,
                double* vplm, double* vcos, double* vsin) const;
    void polleg(double x, double y, double* plm) const;
    void trgev(double x, double y, double* cx, double* sx) const;
    double intmlp(double t, const double* sigma_col, const double* basloc) const;

    // RHS / norm helpers
    void wghpot(const std::vector<double>& phi, Mat& g) const;
    void intrhs(int iat, const double* x_grid, double* xlm) const;
    void adjrhs(int iat, const Mat& xi, double* vlm, double* basloc,
                double* vplm, double* vcos, double* vsin) const;
    void calcv(bool first, int iat, double* pot, const Mat& sigma,
               double* basloc, double* vplm, double* vcos, double* vsin) const;
    void hsnorm(const double* u, double& unorm) const;

    // gradient kernels
    void fdoka(int iat, const Mat& sigma, const double* xi_col, double* basloc,
               double* dbsloc, double* vplm, double* vcos, double* vsin, double fx[3]) const;
    void fdokb(int iat, const Mat& sigma, const Mat& xi, double* basloc,
               double* dbsloc, double* vplm, double* vcos, double* vsin, double fx[3]) const;
    void fdoga(int iat, const Mat& xi, const Mat& phi, double fx[3]) const;

    // Jacobi/DIIS driver; `matvec` is lx (direct) or lstarx (adjoint)
    void jacobiDiis(int n, double tol, const Mat& rhs, Mat& x, int& n_iter, bool& ok,
                    const std::function<void(const Mat&, Mat&)>& matvec) const;

    // index of the real spherical harmonic (l, m), 0-based into nylm arrays
    static int idx(int l, int m) { return l * l + l + m; }
    // switch function and its derivative (5th-degree polynomial)
    static double fsw(double t, double s, double eta);
    static double dfsw(double t, double s, double eta);
};

} // namespace Solvation
} // namespace Curcuma
