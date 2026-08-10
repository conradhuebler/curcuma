/*
 * <QMDFF force-constant parametrisation from a QM Hessian>
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
 * =============================================================================
 * WHAT THIS IS
 * =============================================================================
 *
 * QMDFF derives its force constants from a quantum-mechanical Hessian
 *   S. Grimme, J. Chem. Theory Comput. 2014, 10, 4497-4514 (DOI 10.1021/ct500573f)
 * This file implements that step for curcuma.
 *
 * The central observation is that the QMDFF energy is EXACTLY LINEAR in every
 * constant we want to determine — the bond force constant, the angle force
 * constant, and the torsion / out-of-plane barrier scale:
 *
 *     H_FF(x0) = sum_p  k_p * U_p(x0)
 *
 * where U_p is the Hessian of term p evaluated with its own constant set to 1.
 * Fitting {k_p} to a reference Hessian is therefore a single LINEAR least-squares
 * problem, not an iterative optimisation — no force-field Hessian is ever needed
 * during the fit. (The previous implementation ran Levenberg-Marquardt with a
 * numerically differentiated full force-field Hessian per trial vector, which cost
 * ~1e6 force-field gradient calls per step for a 50-atom molecule.)
 *
 * Two stages:
 *   1. a per-term starting guess — Seminario projection for bonds and angles,
 *      the Rayleigh quotient <U_p,H>/<U_p,U_p> for torsions and out-of-plane
 *      terms (Seminario's method has no defensible analogue there, see below);
 *   2. a global least-squares refinement that accounts for the coupling between
 *      overlapping terms, regularised toward the stage-1 guess and optionally
 *      constrained to non-negative constants.
 *
 * References for stage 1:
 *   J. M. Seminario, Int. J. Quantum Chem. 60 (1996) 1271-1277.
 *   A. E. A. Allen, M. C. Payne, D. J. Cole, J. Chem. Theory Comput. 14 (2018) 274-281
 *     ("Modified Seminario Method" — symmetrise the interaction block).
 *
 * =============================================================================
 * UNITS
 * =============================================================================
 * Everything inside this class is ATOMIC UNITS: Bohr, Hartree, Hartree/Bohr^2.
 * Geometry is accepted in Angstrom (curcuma's UFF/QMDFF convention) and converted
 * internally. The reference Hessian is accepted in Hartree/Angstrom^2 (what
 * Hessian::getRawHessian() returns) and converted by hessianAngstromToBohr(),
 * which is the ONE conversion point in the whole file.
 *
 * Claude Generated (August 2026)
 */

#pragma once

#include "qmdff_terms.h"
#include "src/core/global.h"

#include <Eigen/Dense>

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "json.hpp"
using json = nlohmann::json;

namespace curcuma::qmdff {

/// Which QMDFF term list a fitted constant belongs to.
/**
 * @brief Which force field's functional forms the fitted terms use.
 *
 * Both are linear in exactly one constant per term, so the same linear Hessian fit
 * applies. They differ in the shapes and in which JSON keys carry the term lists:
 *   QMDFF  — bonds / angles / qmdff_torsions   (kernels in qmdff_terms.h)
 *   GFN_FF — bonds / angles / dihedrals / inversions (kernels in gfnff_unit_terms.h)
 *
 * GFN-FF is the better starting point for conformer work: its non-covalent block
 * (hydrogen bonds, D4 dispersion, EEQ electrostatics) carries ~80% of the conformational
 * energy spread and is measurably better than QMDFF's, while the Hessian fit improves
 * the bonded constants it is actually able to constrain. See docs/QMDFF_CONFORMER_USECASE.md.
 */
enum class TermSource : int { QMDFF = 0,
    GFN_FF = 1 };

enum class TermKind : int { Bond = 0,
    Angle = 1,
    Torsion = 2,
    OutOfPlane = 3,
    /// GFN-FF keeps its extra sp3-sp3 gauche torsions in a separate list; they behave
    /// exactly like Torsion but must be written back to `extra_dihedrals`.
    ExtraTorsion = 4 };

const char* termKindName(TermKind kind);

/**
 * @brief Hessian of ONE bonded term with its force constant / scale set to 1.
 *
 * Sparse by construction: a term touches at most four atoms, so only a 12x12
 * block is ever non-zero. Fixed storage keeps the build allocation-free
 * (1.2 kB per term, ~3.5 MB for 3000 terms).
 */
struct TermHessian {
    TermKind kind = TermKind::Bond;
    int list_index = 0; ///< index into parameters["bonds" | "angles" | "qmdff_torsions"]
    int natoms = 2;     ///< 2 (bond), 3 (angle) or 4 (torsion / out-of-plane)
    int group = -1;     ///< fit-parameter index; equals the term index unless terms are tied
    std::array<int, 4> atom{ { -1, -1, -1, -1 } };
    Eigen::Matrix<double, 12, 12> block = Eigen::Matrix<double, 12, 12>::Zero();
    /// Energy and gradient of the term at unit constant, at the same geometry as `block`.
    /// The energy is linear in the constant just like the Hessian is, which is what lets
    /// energy and gradient residuals join the same linear system.
    double energy = 0.0;
    Eigen::Matrix<double, 12, 1> grad = Eigen::Matrix<double, 12, 1>::Zero();
    double norm2 = 0.0; ///< <U,U>, cached

    // -------------------------------------------------------------------------
    // GFN-FF only: the dynamic-r0 coordination-number response
    // -------------------------------------------------------------------------
    /**
     * GFN-FF's bond reference distance follows the coordination numbers,
     *     r0 = (r0_base_i + cnfak_i cn_i + r0_base_j + cnfak_j cn_j + rabshift) ff,
     * and the force field carries that through to the gradient as a chain-rule term
     * dE/dcn * dcn/dx. CN reaches every neighbour, so the exact unit Hessian of a bond
     * is NOT confined to the two bonded atoms. It is, however, still exactly linear in
     * the force constant and of rank two:
     *
     *     U = block(local)  +  extra_coeff * (u w^T + w u^T)
     *
     * with u = dr/dx (nonzero on i and j only) and w = ff (cnfak_i dcn_i/dx +
     * cnfak_j dcn_j/dx) frozen at the reference geometry. Storing the factorisation
     * instead of a dense 3N x 3N matrix keeps the assembly affordable.
     *
     * Both vectors are full 3N and empty when the correction does not apply (every
     * non-bond term, and every term when the source is QMDFF).
     */
    double extra_coeff = 0.0;
    Eigen::VectorXd extra_u; ///< dr/dx, full 3N
    Eigen::VectorXd extra_w; ///< ff * sum_a cnfak_a dcn_a/dx, full 3N, frozen
    double extra_grad_coeff = 0.0; ///< unit gradient picks up extra_grad_coeff * w

    bool hasExtra() const { return extra_u.size() > 0 && extra_w.size() > 0; }

    double at(int a, int al, int b, int be) const { return block(3 * a + al, 3 * b + be); }

    /// H += scale * U, H being the full 3N x 3N matrix.
    void addTo(Matrix& H, double scale) const
    {
        for (int a = 0; a < natoms; ++a)
            for (int al = 0; al < 3; ++al)
                for (int b = 0; b < natoms; ++b)
                    for (int be = 0; be < 3; ++be)
                        H(3 * atom[a] + al, 3 * atom[b] + be) += scale * at(a, al, b, be);
        if (hasExtra() && extra_coeff != 0.0)
            H.noalias() += (scale * extra_coeff)
                * (extra_u * extra_w.transpose() + extra_w * extra_u.transpose());
    }

    /// g += scale * (unit gradient), g being the full 3N vector.
    void addGradTo(Vector& g, double scale) const
    {
        for (int a = 0; a < natoms; ++a)
            for (int c = 0; c < 3; ++c)
                g(3 * atom[a] + c) += scale * grad(3 * a + c);
        if (hasExtra() && extra_grad_coeff != 0.0)
            g.noalias() += (scale * extra_grad_coeff) * extra_w;
    }

    /**
     * @brief <U, H> under the inner product sum_uv M_u M_v U_uv H_uv, M = weight^2.
     * @param mass_weight per-coordinate weight SQUARED (pass an empty vector for M = 1).
     * @note H must be symmetric, which every Hessian here is.
     */
    double dotHessian(const Matrix& H, const Vector& mass_weight) const
    {
        const bool w_on = mass_weight.size() > 0;
        double sum = 0.0;
        for (int a = 0; a < natoms; ++a)
            for (int al = 0; al < 3; ++al) {
                const int mu = 3 * atom[a] + al;
                for (int b = 0; b < natoms; ++b)
                    for (int be = 0; be < 3; ++be) {
                        const int nu = 3 * atom[b] + be;
                        const double m = w_on ? mass_weight(mu) * mass_weight(nu) : 1.0;
                        sum += m * at(a, al, b, be) * H(mu, nu);
                    }
            }
        if (hasExtra() && extra_coeff != 0.0) {
            const Vector mu_v = w_on ? Vector(mass_weight.cwiseProduct(extra_u)) : extra_u;
            const Vector mw_v = w_on ? Vector(mass_weight.cwiseProduct(extra_w)) : extra_w;
            sum += 2.0 * extra_coeff * mu_v.dot(H * mw_v);
        }
        return sum;
    }

    /// <U_p, U_q> under the same inner product as dotHessian().
    double dotTerm(const TermHessian& o, const Vector& mass_weight) const
    {
        const bool w_on = mass_weight.size() > 0;
        auto dotM = [&](const Vector& x, const Vector& y) {
            return w_on ? x.cwiseProduct(mass_weight).dot(y) : x.dot(y);
        };

        // local x local — only the atoms the two terms share contribute
        double sum = 0.0;
        for (int a = 0; a < natoms; ++a) {
            int qa = -1;
            for (int x = 0; x < o.natoms; ++x)
                if (o.atom[x] == atom[a]) { qa = x; break; }
            if (qa < 0) continue;
            for (int b = 0; b < natoms; ++b) {
                int qb = -1;
                for (int x = 0; x < o.natoms; ++x)
                    if (o.atom[x] == atom[b]) { qb = x; break; }
                if (qb < 0) continue;
                for (int al = 0; al < 3; ++al)
                    for (int be = 0; be < 3; ++be) {
                        const double m = w_on
                            ? mass_weight(3 * atom[a] + al) * mass_weight(3 * atom[b] + be)
                            : 1.0;
                        sum += m * at(a, al, b, be) * o.at(qa, al, qb, be);
                    }
            }
        }

        // local x rank-2 (both directions) and rank-2 x rank-2
        const bool e_p = hasExtra() && extra_coeff != 0.0;
        const bool e_q = o.hasExtra() && o.extra_coeff != 0.0;
        if (e_q)
            sum += 2.0 * o.extra_coeff * localTimes(o.extra_u, o.extra_w, mass_weight);
        if (e_p)
            sum += 2.0 * extra_coeff * o.localTimes(extra_u, extra_w, mass_weight);
        if (e_p && e_q)
            sum += 2.0 * extra_coeff * o.extra_coeff
                * (dotM(extra_u, o.extra_u) * dotM(extra_w, o.extra_w)
                    + dotM(extra_u, o.extra_w) * dotM(extra_w, o.extra_u));
        return sum;
    }

private:
    /// sum_uv M_u M_v block_uv x_u y_v, over this term's own atoms only.
    double localTimes(const Vector& x, const Vector& y, const Vector& mass_weight) const
    {
        const bool w_on = mass_weight.size() > 0;
        double sum = 0.0;
        for (int a = 0; a < natoms; ++a)
            for (int al = 0; al < 3; ++al) {
                const int mu = 3 * atom[a] + al;
                for (int b = 0; b < natoms; ++b)
                    for (int be = 0; be < 3; ++be) {
                        const int nu = 3 * atom[b] + be;
                        const double m = w_on ? mass_weight(mu) * mass_weight(nu) : 1.0;
                        sum += m * at(a, al, b, be) * x(mu) * y(nu);
                    }
            }
        return sum;
    }
};

/**
 * @brief One fitting point: a geometry with reference data attached.
 *
 * A Hessian constrains the CURVATURE at one structure. Conformer ranking is about
 * differences in E between distant structures, which curvature alone cannot determine
 * — measured on a 107-atom peptide, a single-basin Hessian fit ranked conformers worse
 * than plain GFN-FF (docs/QMDFF_CONFORMER_USECASE.md). Supplying several basins, and
 * per basin also the reference ENERGY and GRADIENT, is what makes relative basin depths
 * enter the fit. All three residual types are linear in the same constants, so they
 * stack into one least-squares system.
 *
 * Anything left empty is simply not used as a constraint.
 */
struct BasinData {
    Matrix geometry_angstrom;   ///< N x 3
    Matrix hessian_angstrom;    ///< 3N x 3N, Hartree/Angstrom^2, RAW (see Hessian::getRawHessian)
    Matrix gradient_angstrom;   ///< N x 3, Hartree/Angstrom; empty = no gradient constraint
    double energy = 0.0;        ///< reference total energy, Hartree
    bool has_energy = false;
    double weight = 1.0;        ///< relative weight of this basin
    /// GFN-FF only: this basin's own coordination numbers and their frozen derivative.
    /// Both are geometry dependent, so a multi-basin fit that leaves them empty falls
    /// back to basin 0's values and its bond unit Hessians are correspondingly off.
    Vector cn;
    Matrix dcn;
};

/// Knobs of the fit. Defaults are the recommended settings; see docs/QMDFF.md.
struct FitOptions {
    TermSource source = TermSource::QMDFF; ///< which functional forms to fit
    double fd_step_bohr = 1.0e-4; ///< central-difference step for the unit Hessians
    bool mass_weighted = false;   ///< false = plain Cartesian Frobenius inner product

    /**
     * @brief Project translations and rotations out of the target Hessian.
     *
     * OFF by default, and that is a correctness requirement, not a preference. A term's
     * Hessian annihilates infinitesimal rotations only where that term's gradient
     * vanishes: H.R = -dg/dtheta. Bonds and angles sit at their reference by construction
     * and are therefore invariant, but a QMDFF TORSION generally does not, so its unit
     * Hessian carries a genuine rotational component. Projecting the target while the
     * model keeps that component makes the least-squares problem inconsistent — measured
     * on CH3OCH3 it turned an exactly recoverable synthetic system into a 3.5e-3 residual
     * with force constants off by 22%. Projecting the unit Hessians too would fix the
     * inconsistency but makes every U_p dense (3N x 3N) and destroys the O(N) assembly.
     *
     * The information the projection was meant to provide is reported instead, as
     * FitReport::tr_content.
     */
    bool project_tr = false;

    double lambda = -1.0;         ///< Tikhonov strength; < 0 selects it automatically
    std::vector<double> lambda_ladder{ 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0 };
    double lambda_tol = 0.05;     ///< accept the smallest lambda within 5% of the best residual

    /**
     * @brief Separate, much stronger Tikhonov prior for torsion and out-of-plane scales.
     *
     * A torsion constant is a BARRIER HEIGHT, but a Hessian only sees the curvature at one
     * geometry, so the data barely constrains it — and what little signal there is has
     * usually already been absorbed by the overlapping bond and angle terms. Left to the
     * unconstrained solve, torsion scales are driven to their lower bound, which silently
     * produces a force field with free internal rotation.
     *
     * lambda = 1 makes the prior (stay at the incoming, rule-derived barrier) exactly as
     * strong as the data, so a torsion only moves when the Hessian genuinely demands it.
     * Set to 0 to fit torsions on the same footing as bonds and angles.
     */
    double lambda_torsion = 1.0;

    bool nonnegative = true;      ///< constrain the constants to their lower bounds
    double bond_floor = 1.0e-4;   ///< Hartree
    double angle_floor = 1.0e-5;  ///< Hartree
    double torsion_floor = 0.0;
    double oop_floor = 0.0;

    /// Relative weight of each residual block. Each block is normalised by its own norm
    /// first, so these are dimensionless and comparable despite the units differing
    /// (Hartree vs Hartree/Bohr vs Hartree/Bohr^2).
    double weight_hessian = 1.0;
    double weight_energy = 1.0; ///< the block that can constrain conformer ranking
    /**
     * @brief Weight of the gradient block. DEFAULT 0, and that is a correctness default.
     *
     * r0/theta0 come from the structure the parameter set was generated on. At any OTHER
     * basin the bonded terms are displaced from their reference, so the force field has a
     * large gradient there, while the reference gradient at an optimised conformer is
     * ~0. Asking the fit to reproduce that with non-negative constants has exactly one
     * solution: switch the terms off. Measured on the 107-atom peptide with 4 basins, the
     * gradient block drove the relative residual from 0.18 to 1.01 (R^2 = -0.03) and
     * pinned 120 of 364 constants at their floor.
     *
     * Only enable it if every basin shares the reference internal coordinates.
     */
    double weight_gradient = 0.0;

    bool fit_torsions = true;
    bool fit_oop = true;
    /// Tie all torsions about the same central bond to one shared scale. Torsions about a
    /// single bond have near-proportional unit Hessians, so their individual values are not
    /// identifiable; only their sum is. Recommended, and how UFF and QMDFF both think.
    bool tie_torsions_by_central_bond = true;

    int verbosity = 1;
};

/// Everything the caller needs to judge whether the fit is trustworthy.
struct FitReport {
    int n_terms = 0, n_params = 0;
    int n_bonds = 0, n_angles = 0, n_torsions = 0, n_oop = 0;
    double lambda_used = 0.0;
    double rel_residual_initial = 0.0; ///< with the stage-1 guess
    double rel_residual_fitted = 0.0;
    double r_squared = 0.0;
    int n_at_floor = 0;  ///< constants pushed to their lower bound by the constraint
    int n_frozen = 0;    ///< terms excluded because their unit Hessian is numerically zero
    double target_norm = 0.0;
    /// Fraction of the target Hessian living in the translation/rotation subspace. Nonzero
    /// means the geometry is not a stationary point of the reference method; large values
    /// mean a correspondingly large part of the target is not representable by any
    /// internal-coordinate force field.
    double tr_content = 0.0;
    int n_basins = 1;
    double rel_residual_energy = 0.0;   ///< relative residual of the energy block alone
    double rel_residual_gradient = 0.0;
    Vector k_initial, k_fitted;
    Matrix hessian_target;    ///< projected, Hartree/Bohr^2
    Matrix hessian_predicted; ///< sum_p k_p U_p + H_nonbonded, Hartree/Bohr^2
    std::vector<std::pair<double, double>> lambda_curve; ///< (lambda, relative residual)
    std::vector<std::string> warnings;

    json toJson() const;
};

/**
 * @brief Determines QMDFF force constants by fitting a reference Hessian.
 *
 * Deliberately free of any dependency on `capabilities/hessian.h`,
 * `EnergyCalculator` or the thread pool: it takes a matrix in and returns a
 * parameter set, which keeps the round-trip test fast and the layering clean
 * (core must not depend on capabilities).
 */
class QMDFFParametrisation {
public:
    /**
     * @param atoms             Element numbers, size N
     * @param geometry_angstrom N x 3 reference geometry (Angstrom)
     * @param parameters        A QMDFF parameter set as produced by
     *                          ForceFieldGenerator::getParameter() — needs "bonds",
     *                          "angles" and optionally "qmdff_torsions", "qmdff_ncis",
     *                          "qmdff_hbonds"
     */
    QMDFFParametrisation(const std::vector<int>& atoms,
        const Matrix& geometry_angstrom,
        const json& parameters,
        const FitOptions& options = FitOptions());

    /**
     * @brief Coordination numbers, required for GFN-FF's DYNAMIC bond r0.
     *
     * FFWorkspace computes r0 = (r0_base_i + cnfak_i cn_i + r0_base_j + cnfak_j cn_j
     * + rabshift) * ff whenever CN is available, and falls back to Bond::r0_ij otherwise.
     * The fit must use the same r0 or every bond unit Hessian is evaluated at the wrong
     * reference — set this to whatever the force field itself is using.
     *
     * REBUILDS the term data. The constructor already built it, at that point without any
     * coordination numbers, so every GFN-FF bond silently used the static r0_ij fallback.
     * That is a ~5e-3 Bohr error in r0, invisible in the longitudinal curvature but a 2%
     * error in the transverse one (which is proportional to r - r0), and it was the whole
     * residual of the GFN-FF fit verification.
     */
    void setCoordinationNumbers(const Vector& cn);

    /**
     * @brief The FROZEN coordination-number derivative, dcn(a, 3*i+c) = d cn_a / d x_(i,c).
     *
     * In Bohr^-1, at the reference geometry. Supplying it makes the GFN-FF bond unit
     * Hessians EXACT: the force field propagates the dynamic r0 through
     * dE/dcn * dcn/dx, a contribution that reaches beyond the two bonded atoms and that
     * a purely local 12x12 block cannot represent. Without it the bond kernel is ~1.5%
     * off (measured on CH3OCH3), with it the model reproduces the production Hessian.
     *
     * Like setCoordinationNumbers(), this REBUILDS the term data.
     */
    void setCNDerivatives(const Matrix& dcn);

    /**
     * @brief Fit the force constants to a reference Hessian.
     * @param hessian_qm_angstrom Raw Cartesian Hessian in Hartree/Angstrom^2
     *                            (Hessian::getRawHessian(), NOT getHessian()).
     * @param report              Optional diagnostics.
     * @return The input parameter set with the fitted constants substituted.
     */
    json fit(const Matrix& hessian_qm_angstrom, FitReport* report = nullptr);

    /**
     * @brief Fit one parameter set to several basins at once.
     *
     * The model is linear in the constants at EVERY geometry, so k basins simply add
     * their normal equations: A = sum_b A(x_b), b = sum_b b(x_b). One parameter set for
     * the whole ensemble — no basin-assignment problem and no discontinuity between
     * basins, unlike refitting separately per basin.
     *
     * The topology (which bonds/angles/torsions exist, and r0/theta0) comes from the
     * parameter set handed to the constructor, i.e. from ITS geometry; the basins only
     * supply reference data. All basins must therefore be the same molecule in the same
     * atom order and the same protomer.
     */
    json fitMultiBasin(const std::vector<BasinData>& basins, FitReport* report = nullptr);

    /// The single Angstrom^-2 -> Bohr^-2 conversion point of this module.
    static Matrix hessianAngstromToBohr(const Matrix& h);

    /**
     * @brief Supply the Hessian/energy/gradient of everything that is NOT fitted.
     *
     * For QMDFF the remainder is built internally from the qmdff_ncis / qmdff_hbonds
     * lists. GFN-FF's non-covalent block (D4 dispersion, EEQ electrostatics, Born-Mayer
     * repulsion, HB/XB, BATM/ATM) has no such JSON representation, so the caller supplies
     * it instead — most cheaply BY DIFFERENCE, exploiting the same linearity the fit rests
     * on:  H_nb = H_FF(k_current) - sum_p k_current_p U_p.
     *
     * That costs one force-field Hessian and covers every non-fitted term automatically.
     * It does not hide a wrong unit Hessian: the post-fit verification compares against a
     * force-field Hessian at the NEW constants, so any error in U_p survives as
     * (k_new - k_old)(U_p - U_p_true) and shows up in the residual.
     *
     * All arguments in atomic units (Bohr, Hartree).
     */
    void setExternalNonbonded(const Matrix& hessian, double energy, const Matrix& gradient)
    {
        m_h_nonbonded = hessian;
        m_e_nonbonded = energy;
        m_g_nonbonded = gradient;
    }

    /**
     * @brief Hessian of the fitted terms at the constants the parameter set came in with.
     *
     * `H_nonbonded = H_forcefield - currentBondedHessian()` is the cheapest way to obtain
     * the non-fitted remainder for a force field (like GFN-FF) whose non-covalent block
     * has no JSON representation. Kept inside the engine because it must use exactly the
     * same term ordering and sign convention as the fit — duplicating the mapping in the
     * caller silently produced a sign error on GFN-FF bonds.
     *
     * @return 3N x 3N in Hartree/Bohr^2
     */
    Matrix currentBondedHessian() const;

    /// The incoming constants in the INTERNAL (sign-flipped) convention, matching U_p.
    Vector currentConstantsInternalPublic() const { return currentConstantsInternal(); }

    /// The incoming constants, in the force field's own sign convention (one per parameter).
    Vector currentConstants() const;

    /// Access the unit Hessians (valid after the constructor). Exposed for testing.
    const std::vector<TermHessian>& termHessians() const { return m_terms; }

    /// Hessian of all non-fitted terms (NCI + HB/XB), Hartree/Bohr^2. Valid after ctor.
    const Matrix& nonbondedHessian() const { return m_h_nonbonded; }

private:
    // --- construction -------------------------------------------------------
    /// Build unit-constant Hessian, energy and gradient of every fitted term at `geom_bohr`.
    void buildTermData(const Matrix& geom_bohr, std::vector<TermHessian>& terms) const;
    /// Hessian, energy and gradient of everything that is NOT fitted (NCI + HB/XB).
    void buildNonbonded(const Matrix& geom_bohr, Matrix& hessian, double& energy,
        Matrix& gradient) const;
    void assignParameterGroups();

    // --- stages -------------------------------------------------------------
    Vector currentConstantsInternal() const;
    Vector seminarioGuess(const Matrix& h_target) const;
    /// Shared tail of fit() and fitMultiBasin(): stage-1 guess, lambda selection,
    /// constrained solve and write-back. `A`, `b` and `target_norm2` are whatever the
    /// caller assembled — one basin's Hessian, or several basins' Hessian + energy +
    /// gradient blocks.
    json solveAndWriteBack(const Matrix& A, const Vector& b, double target_norm2,
        const Matrix& h_target_for_guess, FitReport& rep) const;
    /// @param lambda per-parameter Tikhonov strength (bonds/angles vs torsions differ)
    Vector solveGlobal(const Matrix& A, const Vector& b, const Vector& k0,
        const Vector& lambda, const Vector& floors) const;

    /// Bounded Lawson-Hanson active set: minimise ||Ak - b|| subject to k >= floors.
    static Vector boundedNNLS(const Matrix& A, const Vector& b, const Vector& floors,
        int max_iter);

    Matrix projectorTR() const;

    // --- data ---------------------------------------------------------------
    std::vector<int> m_atoms;
    Matrix m_geometry_bohr; ///< N x 3, Bohr
    json m_parameters;
    FitOptions m_options;
    int m_natoms = 0;

    std::vector<TermHessian> m_terms;
    int m_nparams = 0;
    Vector m_floors;
    std::vector<char> m_is_torsion_group; ///< per parameter: torsion or out-of-plane?
    Matrix m_h_nonbonded; ///< 3N x 3N, Hartree/Bohr^2
    double m_e_nonbonded = 0.0;
    Matrix m_g_nonbonded; ///< N x 3, Hartree/Bohr
    Vector m_cn;          ///< coordination numbers for GFN-FF's dynamic bond r0 (may be empty)
    Matrix m_dcn;         ///< dcn(a, 3*i+c) = d cn_a / d x_(i,c), Bohr^-1, frozen (may be empty)
    Vector m_weight;      ///< 3N diagonal of the inner product (all ones unless mass-weighted)
};

} // namespace curcuma::qmdff
