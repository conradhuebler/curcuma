/*
 * < Ewald summation for point charges in an orthorhombic periodic cell >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Nov 2026, work package B0 of the periodic box plan): a stand-alone
 * direct Ewald sum for point charges. Not connected to GFN-FF or the EEQ yet (B2/B3).
 * Status: AI-generated, machine-tested (test_cases/unit/test_ewald.cpp), no human testing.
 *
 * ---------------------------------------------------------------------------------------
 * What is computed (atomic units: positions in Bohr, charges in e, energy in Eh, no 4 pi eps0)
 * ---------------------------------------------------------------------------------------
 * Periodic Coulomb energy of N point charges q_i at r_i in a box L_x x L_y x L_z (V = volume):
 *
 *   E = 1/2 sum_n' sum_i sum_j q_i q_j / |r_ij + n L|      (n = image vector, ' = no i=j at n=0)
 *
 * The sum is only conditionally convergent. Ewald (1921) splits 1/r = erfc(a r)/r + erf(a r)/r
 * (a = alpha, the splitting parameter, 1/Bohr) and sums the short-range part in real space and
 * the smooth part in reciprocal space (Allen & Tildesley, "Computer Simulation of Liquids",
 * 2nd ed. 2017, Sec. 6.2; Frenkel & Smit, "Understanding Molecular Simulation", 2nd ed.,
 * Sec. 12.1 -- section numbers from memory, unverified; the formulas below are re-derived
 * and checked by the unit tests, not taken from a page):
 *
 *   E_real  = sum_{i<j} q_i q_j erfc(a r_ij) / r_ij                    (minimum image, r_ij < r_cut)
 *   E_recip = (1/(2V)) sum_{k != 0} (4 pi / k^2) exp(-k^2/(4 a^2)) |S(k)|^2,
 *             S(k) = sum_j q_j exp(i k.r_j),  k = 2 pi (n_x/L_x, n_y/L_y, n_z/L_z)
 *   E_self  = -(a / sqrt(pi)) sum_i q_i^2
 *   E_bg    = -pi Q^2 / (2 a^2 V),  Q = sum_i q_i  (uniform neutralising background; 0 if Q = 0)
 *
 * Boundary condition: "tin-foil" (the infinite lattice is embedded in a conductor, eps' = inf).
 * The surface term 2 pi |sum_i q_i r_i|^2 / ((2 eps' + 1) V) of de Leeuw, Perram & Smith (1980)
 * is therefore zero and NOT added; the k = 0 term is dropped (from memory, unverified citation).
 *
 * Real space under the minimum image only: requires r_cut <= L_min/2. Contributions of further
 * images (including an atom with its own images at distance >= L) are neglected; with the
 * parameter choice below a r_cut = sqrt(-ln eps), so they are below eps.
 *
 * Potentials phi_i = dE/dq_i. E is a quadratic form in q, E = 1/2 q^T A q, so phi = A q and
 * E = 1/2 sum_i q_i phi_i. A is the Ewald matrix the periodic EEQ needs (B3):
 *   A_ij = erfc(a r_ij)/r_ij [i != j, MIC, r < r_cut]
 *        + (4 pi / V) sum_{k != 0} exp(-k^2/(4a^2))/k^2 cos(k.(r_i - r_j))
 *        - (2 a / sqrt(pi)) delta_ij - pi / (a^2 V)
 *
 * ---------------------------------------------------------------------------------------
 * B2 (Nov 2026): Gaussian charges (GFN-FF EEQ), optional per-atom width alp_i (Bohr^2)
 * ---------------------------------------------------------------------------------------
 * GFN-FF interacts its EEQ charges as Gaussian clouds, pair term q_i q_j erf(g_ij r)/r with
 * g_ij = 1/sqrt(alp_i + alp_j) (alp_i = alpeeq_i = (alpha_i + ff q_i)^2, gfnff_ini.f90:718-725), and keeps
 * the Gaussian self-energy q_i^2/sqrt(2 pi alp_i) in its EEQ self term. With erf = 1 - erfc the periodic
 * pair energy of the clouds, for the same tin-foil / k = 0-dropped convention as the point charges, is
 *
 *   E_G = E_pt(real + reciprocal + self + background)
 *         - sum_{(i,j,n) != (i,i,0)} 1/2 q_i q_j erfc(g_ij |r_ij + n|)/|r_ij + n|        (short range, MIC)
 *         + (pi Q / V) sum_i q_i alp_i                                                   (Q != 0 only)
 *
 * Derivation (Poisson summation of the absolutely convergent correction h_ij(r) = -erfc(g_ij r)/r with
 * h_ij(k = 0) = -pi (alp_i + alp_j)): the k != 0 part of the correction is exactly the difference
 * between the Gaussian and the point-charge reciprocal sums, the k = 0 part, -(pi Q/V) sum_i q_i alp_i,
 * is NOT in the convention "k = 0 dropped" and has to be added back; it is the interaction of the
 * neutralising background with the difference between a cloud and a point (space average of
 * -q erfc(r/sqrt alp)/r is -pi q alp). For Q = 0 it vanishes, and E_G -> sum_{i<j} q_i q_j erf(g r)/r
 * for L -> infinity (the point-Ewald sum tends to sum_{i<j} q_i q_j / r, the correction to
 * -sum_{i<j} q_i q_j erfc(g r)/r; the (i,i,0) term of the correction is the Gaussian minus the point
 * self-energy, which both sides exclude, and the images (i,i,n != 0) carry erfc(g L) < 1e-16).
 * The correction is applied only for g_ij r < GAUSSIAN_CORRECTION_XMAX (erfc(6) = 2e-17, below one
 * ulp of the 1/r it corrects); the real-space loop radius r_cut must therefore reach
 * gaussianCorrectionRadius(alp) = XMAX / g_min, g_min = 1/sqrt(2 max_i alp_i), and stay below L_min/2.
 * Derivation and tests: docs/WP-PERIODIC-NONBONDED.md, section B2.
 *
 * Not implemented: PME/particle-mesh (cost here is O(N^1.5) at best), triclinic cells, the
 * virial/pressure tensor, the dipole (vacuum) boundary term, slab correction. Threading is done by the
 * caller over the range functions (FFWorkspace, per-partition accumulators, deterministic order).
 * A virial formula for diagonal strain (derived, untested) is in docs/WP-PERIODIC-NONBONDED.md, section B0.
 */

#pragma once

#include "src/tools/pbc_utils.h"

#include <Eigen/Dense>

#include <complex>
#include <utility>
#include <vector>

namespace curcuma {
namespace ewald {

/// B2: the Gaussian correction -erfc(g r)/r is applied for g r < XMAX; erfc(6) = 2.2e-17 (Claude Generated, Nov 2026)
constexpr double GAUSSIAN_CORRECTION_XMAX = 6.0;

/*! \brief Radius beyond which every Gaussian correction is below XMAX: XMAX * sqrt(2 max_i alp_i) - Claude Generated (Nov 2026, B2) */
double gaussianCorrectionRadius(const Eigen::VectorXd& alp);

/*! \brief Per-atom, per-axis phase factors exp(i 2 pi n x_a / L_a), n = 0..nmax - Claude Generated (Nov 2026, B2)
 * Built once per geometry (phaseTable), read by reciprocalRange() from any thread.
 */
struct EwaldPhaseTable {
    int natoms = 0;
    int nmax[3] = { 0, 0, 0 };
    std::vector<std::complex<double>> table[3]; ///< [atom * (nmax+1) + n]
};

/*! \brief Ewald parameters (atomic units) - Claude Generated (Nov 2026) */
struct EwaldParameters {
    double alpha = 0.0; ///< splitting parameter a (1/Bohr)
    double r_cut = 0.0; ///< real-space cutoff (Bohr), <= L_min/2
    double k_max = 0.0; ///< reciprocal-space cutoff |k| <= k_max (1/Bohr)
};

/*! \brief Choose alpha and k_max for a target tolerance eps at a given real-space cutoff - Claude Generated (Nov 2026)
 * \param cell      active orthorhombic cell (Bohr)
 * \param tolerance eps, e.g. 1e-8; the largest neglected term relative to its prefactor
 * \param r_cut     real-space cutoff (Bohr); <= 0 selects L_min/2
 *
 * Derivation (per-term criterion, not a rigorous bound on the total error):
 *  - Real space: the first neglected pair term is q_i q_j erfc(a r_cut)/r_cut. For x > 1,
 *    erfc(x) < exp(-x^2)/(x sqrt(pi)) < exp(-x^2). Requiring exp(-(a r_cut)^2) = eps gives
 *        a = sqrt(-ln eps) / r_cut.
 *  - Reciprocal space: the k-terms carry the Gaussian factor exp(-k^2/(4 a^2)) (times
 *    4 pi/k^2 |S|^2/V, which decays as well). Requiring exp(-k_max^2/(4 a^2)) = eps gives
 *        k_max = 2 a sqrt(-ln eps).
 *  Refined error estimates exist (Kolafa & Perram 1992, from memory, unverified); the unit
 *  tests measure the actual error against tighter settings instead.
 */
EwaldParameters chooseEwaldParameters(const PBCUtils::OrthoCell& cell, double tolerance, double r_cut = -1.0);

/*! \brief Energy split into the four Ewald terms (Eh) - Claude Generated (Nov 2026) */
struct EwaldEnergyTerms {
    double real = 0.0;
    double reciprocal = 0.0;
    double self = 0.0;
    double background = 0.0;
    /// B2: (pi Q/V) sum_i q_i alp_i for Gaussian charges in a charged cell; 0 for point charges or Q = 0
    double gaussian_background = 0.0;
    double total() const { return real + reciprocal + self + background + gaussian_background; }
    EwaldEnergyTerms& operator+=(const EwaldEnergyTerms& o)
    {
        real += o.real;
        reciprocal += o.reciprocal;
        self += o.self;
        background += o.background;
        gaussian_background += o.gaussian_background;
        return *this;
    }
};

/*! \brief Everything one evaluation can return - Claude Generated (Nov 2026) */
struct EwaldResult {
    EwaldEnergyTerms terms;
    Eigen::MatrixXd gradient;  ///< N x 3, dE/dr_i (Eh/Bohr); empty unless requested
    Eigen::VectorXd potentials; ///< N, phi_i = dE/dq_i (Eh/e); empty unless requested
};

/*! \brief Direct Ewald sum for point charges in an orthorhombic cell - Claude Generated (Nov 2026)
 *
 * Positions: N x 3 matrix (Bohr, rows are atoms, need not be wrapped into the box).
 * Charges: N-vector (e). Every call is independent (no cached geometry).
 */
class EwaldSum {
public:
    /*! \brief Set up the cell, alpha, r_cut and the list of k vectors
     * \throws std::invalid_argument for an inactive cell, alpha <= 0, r_cut > L_min/2 or k_max <= 0
     */
    EwaldSum(const PBCUtils::OrthoCell& cell, double alpha, double r_cut, double k_max);
    EwaldSum(const PBCUtils::OrthoCell& cell, const EwaldParameters& p)
        : EwaldSum(cell, p.alpha, p.r_cut, p.k_max)
    {
    }

    /// Total energy E (Eh). gaussian_alp (B2): per-atom alp_i (Bohr^2) of Gaussian charges, nullptr = points.
    double energy(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
        const Eigen::VectorXd* gaussian_alp = nullptr) const;
    /// Gradient dE/dr_i (N x 3, Eh/Bohr; the force is its negative, as everywhere in curcuma).
    Eigen::MatrixXd gradient(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
        const Eigen::VectorXd* gaussian_alp = nullptr) const;
    /// Potentials phi_i = dE/dq_i = (A q)_i without building A (what an iterative EEQ solver needs).
    Eigen::VectorXd potentials(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
        const Eigen::VectorXd* gaussian_alp = nullptr) const;
    /// Energy terms plus, on request, gradient and potentials in one pass over pairs and k vectors.
    EwaldResult evaluate(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
        bool want_gradient, bool want_potentials, const Eigen::VectorXd* gaussian_alp = nullptr) const;
    /// Dense Ewald matrix A (N x N, symmetric), A q = phi. O(N^2 N_k): small test systems only.
    /// With gaussian_alp: A_ij -= erfc(g_ij r)/r (pairs, g r < XMAX) and A_ij += (pi/V)(alp_i + alp_j) (all i, j).
    Eigen::MatrixXd matrix(const Eigen::MatrixXd& positions, const Eigen::VectorXd* gaussian_alp = nullptr) const;

    /*! \brief Real-space pair kernel erfc(a r)/r and its r-derivative (point charges) */
    static double realSpaceKernel(double r, double alpha, double& dkernel_dr);

    // ---- Range interface (B2, Claude Generated Nov 2026): the same four pieces evaluate() is made of,
    // callable per index range from several threads. Every function ADDS into `acc` (terms, and
    // gradient / potentials when they have size N x 3 / N); the caller allocates and reduces.

    /// Pairs (i < j) with minimum-image distance < r_cut; linked cells when every axis has >= 3 cells.
    std::vector<std::pair<int, int>> realSpacePairs(const Eigen::MatrixXd& positions) const;
    /// Phase factors per atom and axis for reciprocalRange().
    EwaldPhaseTable phaseTable(const Eigen::MatrixXd& positions) const;
    /// Real-space terms of pairs[begin, end): point kernel, minus the Gaussian correction when gaussian_alp is given.
    void realSpaceRange(const Eigen::MatrixXd& positions, const Eigen::VectorXd& charges,
        const std::vector<std::pair<int, int>>& pairs, int begin, int end, const Eigen::VectorXd* gaussian_alp,
        bool want_gradient, bool want_potentials, EwaldResult& acc) const;
    /// Reciprocal terms of the k vectors [k_begin, k_end) (indices into the stored half-space list).
    void reciprocalRange(const EwaldPhaseTable& table, const Eigen::VectorXd& charges, int k_begin, int k_end,
        bool want_gradient, bool want_potentials, EwaldResult& acc) const;
    /// Self, neutralising background and (with gaussian_alp) the Gaussian background term; O(N), one call.
    void selfAndBackground(const Eigen::VectorXd& charges, const Eigen::VectorXd* gaussian_alp,
        bool want_potentials, EwaldResult& acc) const;

    const EwaldParameters& parameters() const { return m_par; }
    int numberOfKVectors() const { return static_cast<int>(m_kvectors.size()); }
    double volume() const { return m_volume; }

private:
    /// One reciprocal vector of the half space (n and -n give the same term; only one is stored).
    struct KVector {
        int n[3];      ///< integer indices
        double k[3];   ///< 2 pi n_a / L_a (1/Bohr)
        double weight; ///< (4 pi / V) exp(-k^2/(4 a^2)) / k^2; already doubled for the -k partner
    };

    PBCUtils::OrthoCell m_cell;
    EwaldParameters m_par;
    double m_volume = 0.0;
    int m_nmax[3] = { 0, 0, 0 };
    std::vector<KVector> m_kvectors;
};

} // namespace ewald
} // namespace curcuma
