/*
 * <Fragment-wise SCF preconditioner for the native xTB Broyden mixer>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Theory:
 *   The native GFN SCF iterates the fixed-point map x -> G(x) on the SCC vector
 *       x = [ q_sh (nsh) | dp (3*nat, index nsh + k*nat + atom) | qp (6*nat, index nsh + 3*nat + k*nat + atom) ]
 *   (GFN2; GFN1 mixes q_sh only) with the residual F(x) = G(x) - x and its Jacobian
 *       J = dF/dx = G'(x) - I.
 *   A Newton step would be dx = -J^-1 F. The modified Broyden mixer (Johnson 1988; Eyert 1996,
 *   see broyden_mixer.h) builds an inverse-Jacobian estimate from the iteration history on top
 *   of a seed, which is the scalar step alpha*F (alpha = scf_damping). This module replaces the
 *   seed by M*F with M ~ -J^-1 (BroydenMixer::setStepHook): the closer M is to -J^-1, the less
 *   the secant updates have to learn and the fewer SCF iterations are needed. The fixed point
 *   G(x*) = x* is unchanged, because the convergence test still acts on F itself; only the path
 *   to it changes.
 *
 *   'frag' builds M from the molecules (covalent fragments) of the system: for every fragment
 *   the Jacobian J_f of the isolated, neutral fragment is computed by central finite differences
 *   of its own SCC map at its own converged state, and B_f = -J_f^-1 is embedded at the global
 *   indices of that fragment's shells and atoms. Coupling between fragments is not represented
 *   in M; Broyden learns it from the history.
 *
 *   Form (additive Schwarz):  M = sum_i R_i^T D_i B_i R_i  +  C
 *     R_i  restriction to the index set of block i (gather), R_i^T the prolongation (scatter-add)
 *     B_i  local block (here -J_f^-1 of one fragment)
 *     D_i  diagonal weights of block i (here 1: the blocks do not overlap)
 *     C    diagonal fallback alpha on every index that no block covers.
 *   This is the one-level part of a Schwarz decomposition (no coarse level).
 *   Large fragments (scf_precond_segment_atoms > 0, Oct 2026): restricted additive Schwarz. The
 *   fragment is cut into segments (segmentFragment); each segment's Jacobian is computed on an
 *   overlapping, hydrogen-capped local cluster, B_i = -J_loc^-1 is restricted to the cluster's real
 *   atoms and only the rows of the segment's own atoms are kept (D_i = 1 there, 0 elsewhere; the
 *   zero rows are stored, so the weights stay empty and the device copy is unchanged).
 *   Cost of apply(): O(sum_i dim_i^2); no dense dim x dim matrix is formed.
 *
 * References:
 *   D. D. Johnson, Phys. Rev. B 38, 12807 (1988) - modified Broyden method.
 *   V. Eyert, J. Comput. Phys. 124, 271 (1996) - Broyden methods, inverse-Jacobian seed.
 *   B. Smith, P. Bjorstad, W. Gropp, Domain Decomposition (Cambridge 1996) - additive Schwarz.
 *   X.-C. Cai, M. Sarkis, SIAM J. Sci. Comput. 21, 792 (1999) - restricted additive Schwarz.
 *
 * Claude Generated (Oct 2026). AI-generated, machine-tested only; human production testing
 * pending. See docs/SCF_PRECOND.md.
 */

#pragma once

#include "src/core/global.h"

#include <Eigen/Dense>
#include <functional>
#include <string>
#include <vector>

namespace curcuma::xtb {

/// Dense row-major matrix type of the preconditioner (same as the project Matrix).
using PrecondMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/**
 * @brief Residual-based fallback from a Broyden preconditioner M to the scalar step
 *        (Claude Generated, Oct 2026).
 *
 * Watches f_n = ||F_n||_2 of the full SCC residual F = G(x_n) - x_n, which the convergence test
 * computes anyway; no extra evaluation of the map G. Rules (spec string, -scf_precond_guard or
 * CURCUMA_SCF_PRECOND_GUARD):
 *   g1:<tau>[:keep][:best]          growth: fires at n >= 2 if f_n > tau * min_{k<n} f_k
 *   g2:<n0>[:<ratio>][:keep][:best] stagnation: fires at n == n0 if f_n / f_0 >= ratio; default ratio
 *                                   is the worst scalar (no M) reduction of the prototype test set,
 *                                   rounded up: 0.15 for n0 = 4, 0.065 for n0 = 6 (24 scalar SCFs,
 *                                   6 systems x 4 geometries)
 * "keep" keeps the Broyden history (BroydenMixer::dropStepHook), otherwise it is reset.
 * "best" restarts from the input x_k with the smallest ||F_k|| seen so far (its residual F_k is
 * stored, so the scalar step x_k + alpha*F_k needs no evaluation); without it the scalar step
 * starts from the current x_n.
 * Empty or "off": inactive.
 */
struct ScfPrecondGuard {
    enum class Rule { Off, Growth, Stagnation } rule = Rule::Off;
    double tau = 0.0;     ///< G1 growth factor
    int    n0 = 0;        ///< G2 check iteration
    double ratio = 0.0;   ///< G2 required reduction f_n0 / f_0
    bool   keep = false;  ///< keep the Broyden history on fallback
    bool   best = false;  ///< restart from the best input seen so far
    Vector x_best, F_best;  ///< input and residual with the smallest ||F|| so far (only with best)
    std::string name;
    double f0 = 0.0, fmin = 0.0;

    // Special members defined out of line (xtb_precond.cpp) so that -Winline stays quiet.
    ScfPrecondGuard();
    ~ScfPrecondGuard();
    ScfPrecondGuard(const ScfPrecondGuard&);
    ScfPrecondGuard(ScfPrecondGuard&&) noexcept;
    ScfPrecondGuard& operator=(const ScfPrecondGuard&);
    ScfPrecondGuard& operator=(ScfPrecondGuard&&) noexcept;

    /// @brief Parse a guard spec; unknown text warns (naming @p origin, only if @p warn) and leaves
    ///        the guard off.
    static ScfPrecondGuard parse(const std::string& spec, const std::string& origin, bool warn = true);

    /// @brief Parse CURCUMA_SCF_PRECOND_GUARD (unset or "off": guard off).
    static ScfPrecondGuard fromEnv();

    /// @brief Record f_n = ||F_n|| of iteration n (input x, residual F); true if the rule says:
    ///        drop M for the rest of this SCF.
    bool fires(int n, const Vector& x, const Vector& F);

    /// @brief The rule of fires() on the norm alone (Claude Generated, Oct 2026): records
    ///        f_n = @p f, sets @p new_best (if non-null) when x_n/F_n of this iteration are the
    ///        new best pair that "best" has to keep, and returns true if M is to be dropped. The
    ///        caller stores x_n/F_n itself (the device Broyden keeps them as device copies).
    bool firesNorm(int n, double f, bool* new_best);
};

/**
 * @brief Additive-Schwarz block preconditioner M = sum_i R_i^T D_i B_i R_i + C
 *        (Claude Generated, Oct 2026). See the file header for the mathematics.
 */
class SchwarzPreconditioner {
public:
    /// One block: global indices (R_i), local matrix B_i, optional weights D_i (empty = 1).
    struct Block {
        std::vector<int> index;
        PrecondMatrix    B;
        Vector           weight;
        Block();
        ~Block();
        Block(Block&&) noexcept;
        Block& operator=(Block&&) noexcept;
    };

    /// @param dim   length of the SCC vector
    /// @param alpha diagonal fallback C on every index no block covers (scf_damping)
    SchwarzPreconditioner(int dim, double alpha);

    /// @brief Add a block. Indices must lie in [0, dim); B must be index.size() square.
    /// @param covers local positions k (index[k]) on which the fallback C is removed; empty = all.
    ///        A restricted-additive-Schwarz block (Claude Generated, Oct 2026) covers only the rows
    ///        of its own atoms; its other rows are zero and its other indices keep C unless another
    ///        block owns them.
    /// @return false (and nothing added) on inconsistent input.
    bool addBlock(Block block, const std::vector<char>& covers = {});

    /// @brief M*v, O(sum_i dim_i^2).
    Vector apply(const Vector& v) const;

    int dim() const { return m_dim; }
    double alpha() const { return m_alpha; }
    std::size_t blockCount() const { return m_blocks.size(); }
    /// Number of doubles stored in the blocks and weights (memory = 8 bytes each).
    std::size_t storedDoubles() const;

    /**
     * @brief Flat copy of M for a device upload (Claude Generated, Oct 2026).
     *
     * Block b covers idx[idx_off[b] .. idx_off[b+1]) and its row-major matrix starts at
     * B[b_off[b]]; fallback is C (alpha on uncovered indices, 0 on covered ones); max_dim is the
     * largest block. Returns false (outputs undefined) if a block carries weights D_i, which the
     * device kernel does not apply, or if an offset does not fit an int.
     */
    bool flatten(std::vector<int>& idx_off, std::vector<int>& idx, std::vector<int>& b_off,
                 std::vector<double>& B, std::vector<double>& fallback, int& max_dim) const;

private:
    int m_dim;
    double m_alpha;
    std::vector<Block> m_blocks;
    Vector m_fallback;   ///< alpha on uncovered indices, 0 on covered ones
};

/**
 * @brief Central-difference Jacobian of the SCC residual (Claude Generated, Oct 2026).
 *
 * J[i,k] = (G_i(x*+h e_k) - G_i(x*-h e_k)) / 2h - delta_ik. When @p J_fwd is non-null it also
 * receives the forward difference (G_i(x*+h e_k) - g_star_i) / h - delta_ik from the same
 * evaluations. 2*dim evaluations of @p evalG. Row-major dim x dim, row = output component.
 */
void numericalScfJacobian(const Vector& x_star, const Vector& g_star,
                          const std::function<Vector(const Vector&)>& evalG, double h,
                          std::vector<double>& J, std::vector<double>* J_fwd);

/**
 * @brief Covalent fragments (molecules) by a bond graph (Claude Generated, Oct 2026).
 *
 * Atoms i, j are bonded if r_ij < 1.25 * (r_cov(i) + r_cov(j)), with the single-bond covalent
 * radii of Pyykko and Atsumi (Chem. Eur. J. 15 (2009) 186; Elements::CovalentRadius). Fragments
 * are the connected components, each sorted by atom index, ordered by their smallest atom index.
 * O(nat^2).
 *
 * @param atoms     atomic numbers
 * @param geometry  nat x 3 coordinates in Angstrom
 */
std::vector<std::vector<int>> detectCovalentFragments(const std::vector<int>& atoms,
                                                      const Matrix& geometry);

/**
 * @brief One segment of a large fragment for the restricted additive Schwarz (RAS) build
 *        (Claude Generated, Oct 2026).
 *
 * own:    atoms whose rows of M this segment provides (disjoint between segments of a fragment)
 * domain: own atoms plus the neighbourhood the local cluster is computed for, ascending indices
 * caps:   ncap x 3 positions (Angstrom) of the hydrogen atoms that saturate the cut bonds
 */
struct ScfPrecondSegment {
    std::vector<int> own;
    std::vector<int> domain;
    Matrix caps;
    ScfPrecondSegment();
    ~ScfPrecondSegment();
    ScfPrecondSegment(const ScfPrecondSegment&);
    ScfPrecondSegment(ScfPrecondSegment&&) noexcept;
    ScfPrecondSegment& operator=(const ScfPrecondSegment&);
    ScfPrecondSegment& operator=(ScfPrecondSegment&&) noexcept;
};

/**
 * @brief Cut one covalent fragment into segments with overlapping domains for a restricted
 *        additive Schwarz preconditioner (Claude Generated, Oct 2026).
 *
 * Restricted additive Schwarz (Cai and Sarkis, SIAM J. Sci. Comput. 21 (1999) 792; Smith, Bjorstad,
 * Gropp, Domain Decomposition, 1996): every subdomain i is solved on an overlapping domain, but
 * only the rows of its own, non-overlapping part enter M = sum_i R_i^T D_i B_i R_i (D_i = 1 on own
 * components, 0 elsewhere), so no component is counted twice.
 *
 * Construction: only single bonds outside rings are cut. A heavy-heavy bond (bond criterion of
 * detectCovalentFragments) is cuttable if it is a bridge of the heavy-atom bond graph and
 * r_ij / (r_cov,i + r_cov,j) >= 0.95 (Pyykko-Atsumi radii; aromatic ~0.92, C=O ~0.88, amide C-N
 * ~0.94, single ~1.0). Heavy atoms (Z > 1; all atoms if the fragment has none) joined by uncuttable
 * bonds form units that are never split. Units in breadth-first order (start at the smallest index,
 * neighbours ascending) are collected into a segment until it holds >= @p segment_heavy_atoms heavy
 * atoms; each H joins its nearest heavy atom. Domain: every unit with a heavy atom within @p radius
 * (Angstrom) of an own heavy atom, with the H. Each cut heavy-heavy bond h-k (k outside) gets a
 * hydrogen cap on the axis h->k at r_cov(h) + r_cov(H). O(nf^2) for nf fragment atoms.
 *
 * @param atoms     atomic numbers of the whole system
 * @param geometry  nat x 3 coordinates in Angstrom
 * @param frag      atom indices of the fragment (ascending)
 */
std::vector<ScfPrecondSegment> segmentFragment(const std::vector<int>& atoms, const Matrix& geometry,
                                               const std::vector<int>& frag, int segment_heavy_atoms,
                                               double radius);

} // namespace curcuma::xtb
