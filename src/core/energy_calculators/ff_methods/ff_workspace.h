/*
 * <FFWorkspace - Unified Force Field Workspace for Curcuma>
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
 * Claude Generated (March 2026): Unified workspace replacing ForceField+ForceFieldThread
 * for GFN-FF calculations. Single shared workspace with per-partition accumulators.
 *
 * Architecture:
 *   T=1: Direct function calls, zero pool overhead
 *   T>1: pool->enqueue() per partition, barrier between phases, reduce at end
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 */

#pragma once

#include "src/core/global.h"
#include <limits>
#include "gfnff_parameters.h"
#include "ff_terms.h"  // Bond, Angle, Dihedral, Inversion, vdW, EQ, CNDerivStore, GeoGradMatrix
#include "gfnff_param_tables.h"  // Claude Generated (Sep 2026): runtime gen scalars
#include "rev_bond_order.h"      // Claude Generated (Sep 2026): rev-gfnff stage 1 switching functions

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <Eigen/Dense>
#include <array>
#include <functional>
#include <Eigen/Sparse>

#include <functional>
#include <future>
#include <memory>
#include <vector>

/**
 * @brief Energy components for all force field terms
 *
 * Claude Generated (March 2026): Replaces scattered double members
 * in ForceFieldThread. Provides reset() and operator+= for reduction.
 */
struct FFEnergyComponents {
    double bond = 0, angle = 0, dihedral = 0, inversion = 0;
    double dispersion = 0;
    double vdw = 0, rep = 0;                    // UFF/QMDFF LJ non-bonded (attractive/repulsive)
    double bonded_rep = 0, nonbonded_rep = 0;   // GFN-FF exponential repulsion
    double coulomb = 0, hbond = 0, xbond = 0;
    double atm = 0, batm = 0, stors = 0;
    double over_coord = 0;                      // rev-gfnff stage 1: over-coordination penalty (Sep 2026)
    double sqe_hardness = 0;                    // rev-gfnff stage 2: bond hardness 1/2 kappa_ij p_ij^2 (Sep 2026)
    // Claude Generated (May 2026, HB-investigation): per-case HB diagnostic split.
    // Sum of cases equals hbond. Counts compare against Fortran nhb1/nhb2.
    double hbond_case1 = 0, hbond_case2 = 0, hbond_case3 = 0, hbond_case4 = 0;
    int hbond_case1_count = 0, hbond_case2_count = 0, hbond_case3_count = 0, hbond_case4_count = 0;

    void reset() {
        bond = angle = dihedral = inversion = 0;
        dispersion = 0;
        vdw = rep = 0;
        bonded_rep = nonbonded_rep = 0;
        coulomb = hbond = xbond = 0;
        atm = batm = stors = 0;
        over_coord = 0;
        sqe_hardness = 0;
        hbond_case1 = hbond_case2 = hbond_case3 = hbond_case4 = 0;
        hbond_case1_count = hbond_case2_count = hbond_case3_count = hbond_case4_count = 0;
    }

    FFEnergyComponents& operator+=(const FFEnergyComponents& o) {
        bond += o.bond; angle += o.angle; dihedral += o.dihedral;
        inversion += o.inversion; dispersion += o.dispersion;
        vdw += o.vdw; rep += o.rep;
        bonded_rep += o.bonded_rep; nonbonded_rep += o.nonbonded_rep;
        coulomb += o.coulomb; hbond += o.hbond; xbond += o.xbond;
        atm += o.atm; batm += o.batm; stors += o.stors;
        over_coord += o.over_coord;
        sqe_hardness += o.sqe_hardness;
        // Claude Generated (May 2026, HB-investigation): per-case reduction
        hbond_case1 += o.hbond_case1; hbond_case2 += o.hbond_case2;
        hbond_case3 += o.hbond_case3; hbond_case4 += o.hbond_case4;
        hbond_case1_count += o.hbond_case1_count;
        hbond_case2_count += o.hbond_case2_count;
        hbond_case3_count += o.hbond_case3_count;
        hbond_case4_count += o.hbond_case4_count;
        return *this;
    }

    double total() const {
        return bond + angle + dihedral + inversion + dispersion +
               vdw + rep +
               bonded_rep + nonbonded_rep + coulomb + hbond + xbond +
               atm + batm + stors + over_coord + sqe_hardness;
    }
};

/**
 * @brief Per-term timing for FF energy calculation (CPU-sum across threads)
 *
 * Claude Generated (May 2026): Tracks ms spent in each calc* function.
 * Each accumulator owns one; reduce() sums them for CPU-sum timing display
 * in GFNFFEnergyReport.
 */
struct FFTermTimings {
    double bonds = -1, angles = -1, dihedrals = -1, inversions = -1, stors = -1;
    double dispersion = -1, bonded_rep = -1, nonbonded_rep = -1;
    double coulomb = -1, hbond = -1, xbond = -1, atm = -1, batm = -1;

    void reset() {
        bonds = angles = dihedrals = inversions = stors = -1;
        dispersion = bonded_rep = nonbonded_rep = -1;
        coulomb = hbond = xbond = atm = batm = -1;
    }
    // Sum across partitions: -1 means "not measured" — only sum positive contributions.
    static void merge(double& dst, double src) {
        if (src < 0) return;
        if (dst < 0) dst = 0;
        dst += src;
    }
    FFTermTimings& operator+=(const FFTermTimings& o) {
        merge(bonds, o.bonds); merge(angles, o.angles);
        merge(dihedrals, o.dihedrals); merge(inversions, o.inversions);
        merge(stors, o.stors); merge(dispersion, o.dispersion);
        merge(bonded_rep, o.bonded_rep); merge(nonbonded_rep, o.nonbonded_rep);
        merge(coulomb, o.coulomb); merge(hbond, o.hbond); merge(xbond, o.xbond);
        merge(atm, o.atm); merge(batm, o.batm);
        return *this;
    }
};

/**
 * @brief Per-partition accumulator for gradient and energy
 *
 * Claude Generated (March 2026): Each partition (thread) writes to its own
 * accumulator. After all partitions complete, accumulators are reduced.
 * For T=1, acc[0] IS the result (no reduce needed).
 */
struct FFAccumulator {
    GeoGradMatrix gradient;          ///< N×3 gradient accumulator (WP-G: RowMajor)
    /// rev-gfnff stage 1b: acc += scale * o for everything a bonded kernel writes (energies, gradient, dEdcn, components)
    Vector dEdcn;             ///< N dE/dCN for chain-rule
    Vector dEdcn_bond;        ///< N bond-only dE/dCN for per-component attribution
    Vector dEdshare;          ///< rev-gfnff 3a(ii): N dE/d(sum_i w_ik) of the valence share
    FFEnergyComponents energy;
    FFTermTimings timings;    ///< Claude Generated (May 2026): per-term ms (this partition only)

    // Optional per-component gradients (only allocated if store_components=true)
    GeoGradMatrix grad_bond, grad_angle, grad_torsion, grad_repulsion;
    GeoGradMatrix grad_coulomb, grad_dispersion, grad_hb, grad_xb, grad_batm, grad_atm;

    bool has_components = false;

    void reset(int natoms, bool do_gradient, bool store_components) {
        energy.reset();
        timings.reset();
        if (do_gradient) {
            if (gradient.rows() != natoms || gradient.cols() != 3)
                gradient.resize(natoms, 3);
            gradient.setZero();
            if (dEdcn.size() != natoms) dEdcn.resize(natoms);
            dEdcn.setZero();
            if (dEdcn_bond.size() != natoms) dEdcn_bond.resize(natoms);
            dEdcn_bond.setZero();
            if (dEdshare.size() != natoms) dEdshare.resize(natoms);
            dEdshare.setZero();
        }
        has_components = store_components;
        if (store_components && do_gradient) {
            auto resetMat = [&](GeoGradMatrix& m) {  // WP-G: lambda matches RowMajor member type
                if (m.rows() != natoms || m.cols() != 3) m.resize(natoms, 3);
                m.setZero();
            };
            resetMat(grad_bond); resetMat(grad_angle); resetMat(grad_torsion);
            resetMat(grad_repulsion); resetMat(grad_coulomb); resetMat(grad_dispersion);
            resetMat(grad_hb); resetMat(grad_xb); resetMat(grad_batm); resetMat(grad_atm);
        }
    }
};

/**
 * @brief Index ranges for one partition into the master interaction lists
 *
 * Claude Generated (March 2026): Linear ranges (begin/end) for simple lists,
 * index vectors for three-body terms that need modulo-based distribution.
 */
/**
 * @brief rev-gfnff stage 1 settings of the workspace (Claude Generated, Sep 2026)
 *
 * Off by default: every kernel then evaluates exactly the GFN-FF expression. With
 * `enabled` the bonded terms are multiplied by the continuous bond order of their bonds,
 * the bonded/non-bonded repulsion of a pair is blended with the same weight, and the
 * over-coordination energy is added (docs/REV_GFNFF_ROADMAP.md, WP3; rev_bond_order.h).
 */
/**
 * @brief rev-gfnff stage 2: one split-charge pair of the current topology corner
 *
 * Claude Generated (Sep 2026), docs/REV_GFNFF_STAGE2.md. The charges themselves are solved
 * by EEQSolver::calculateSplitCharges(); the workspace only carries the resulting split
 * charge p of every pair so that the hardness energy and its r-derivative can be evaluated
 * on the same geometry as every other term:
 *
 *     E_sqe = sum_(ij) 1/2 kappa0_ij / b_ij(r) p_ij^2
 *     dE/dr = -1/2 p_ij^2 kappa0_ij / b_ij^2 db_ij/dr        (E is variational in p)
 *
 * b is the bond order of the E_over switch (RevSettings::R2 / bo2_width), so the pair term
 * grows without bound as the bond breaks — which is exactly what pins the charge back onto
 * the separating fragment.
 */
struct SqePairData {
    int i = -1, j = -1;    ///< the pair; p flows i -> j
    double p = 0.0;        ///< split charge (e), from the SQE solve of this corner
    double kappa0 = 0.0;   ///< 1/2 (kappa_Z(i) + kappa_Z(j)) in Eh
    /// rev-gfnff P3 (Claude Generated, Sep 23, 2026): flat excess-electron hardness x_ij kappa_x
    /// (Eh), added to kappa(b) with no b-derivative - the same number EEQSolver::SqePair carries.
    double kappa_x = 0.0;
    /// rev-gfnff P3 alternative "frac" (Claude Generated, Sep 23, 2026): fractional-charge
    /// correction factor c (see EEQSolver::SqePair::frac_c); E_x = 1/2 c K_ij(r) q_i q_j.
    double frac_c = 0.0;
    /// rev-gfnff P3 alternative "harris" (Claude Generated, Sep 24, 2026,
    /// _log/P2P3_HARRIS_STATUS.md): topological excess-electron count x_ij of the pair (0 = inert).
    /// Adds E = x_ij g(r_ij) (RevHarrisTable::harrisG) - never enters the charge solve.
    double harris_x = 0.0;
    /// rev_excess_bond_extend: the pair's b is revOrder(r / b_scale) (see EEQSolver::SqePair::b_scale)
    double b_scale = 1.0;
};

struct RevSettings {
    bool enabled = false;
    bool bond_weight = true;      ///< bond well  x  b_ij
    bool term_weights = true;     ///< angle/torsion/inversion damping  x  product of bond weights
    bool blend_repulsion = true;  ///< E_rep = b E_bonded + (1-b) E_nonbonded on every pair that carries both sets
    bool over_coord = true;       ///< E_over,i = p_Z sp(sum_j b_ij BO_ij - Val_Z)^2
    /// rev-gfnff stage 3a(ii) (Sep 2026): valence share of the bond well,
    ///     E_ij = -k_b e^{-a dr^2} * w_ij * c_ij,   c_ij = 1/2 (f_i + f_j),
    ///     f_i = clip((Val_i - sum_{k != j} w_ik) / w_ij, 0, 1)   (C1 soft clip)
    /// with S = the term weight w of the pair (the same switch the well is multiplied with) and
    ///     Val_i = Val_Z(i) + softplus(sum_k shareClip(b_ik) - Val_Z(i))
    /// the HYPERVALENT-CORRECT effective valence: the nominal sigma valence plus a smooth,
    /// saturating term in the settled-partner count (b = the tight bond order). It equals the
    /// nominal valence bit-for-bit for every atom with at most Val_Z partners (all ordinary
    /// chemistry, saturated or not), so an ammonium, hydronium or perchlorate keeps c = 1 on its
    /// genuine bonds, while a migrating H that hangs two partial bonds off the same valence still
    /// shares them. Removes the missing valence conservation of pairwise wells: two wells that
    /// share one valence sum to one well's worth instead of two.
    bool valence_share = true;
    /// rev-gfnff stage 3a(ii) (Claude Generated, Sep 14, 2026): the SMOOTH 1,3 proxy of the share
    /// - a pair does not claim valence for the bond order that leaks onto it from a SETTLED
    /// shared partner (g_p = shareClip(1 - sum_k sigma_ik sigma_jk), sigma = shareSettled(b)),
    /// which is what separates a 1,3 contact (the F...F pairs of a compressed BF4-) from the
    /// migrating pair of an exchange transition state without a discrete topology test. See the
    /// member documentation.
    /// DEFAULT OFF (measured, test_cases/revgfnff/_log/PROXY_STATUS.md): it fixes the BF4- probe
    /// exactly (0.0000 vs +569.70 kcal/mol), is FD-exact and preserves every bit-identity, but it
    /// gives a 1,3 contact pair the FULL well of its pair (c = 1) where the plain share suppresses
    /// it, and in hot react MD those wells appear and vanish -> max |dE_jump| 3331 vs 471 kJ/mol
    /// and 24 vs 3 events >= 50 kJ/mol on the 22-cell grid. So it is an opt-in experiment, not the
    /// default; the design decision it feeds is recorded in the vault note of
    /// docs/REV_GFNFF_ROADMAP.md.
    bool share_onethree = false;
    /// rev-gfnff stage 3a(ii) (Claude Generated, Sep 15, 2026): hydrogen keeps its NOMINAL valence
    /// 1 in the share budget. The softplus budget Val_i = Val_Z + G(settled_i - Val_Z) is right for
    /// a hypervalent centre (an ammonium N really does carry four bonds) but wrong for hydrogen: an
    /// H between two partners is a 3c-2e bridge with ONE valence, and letting its budget grow to 2
    /// hands both of its partial wells a full share the moment the second partner's tight bond
    /// order crosses the settled window - a step change of the bond energy with no topology event
    /// (measured: c2h6/T2000 frame 16, Val(H) 1.06 -> 1.93 and the bond term -1.05 -> -1.52 Eh
    /// inside one 0.25 fs step, which then drives the pair to r = 0.56 a0 and 62 000 K). With the
    /// flag on, Val_H = Val_Z(H) = 1 exactly and its derivative channel is 0, so a bridging H
    /// shares its one valence between its two wells. Every other element is untouched.
    /// DEFAULT ON since Sep 18, 2026 (operator decision after FABLE_REVIEW_2 A.2); false
    /// reproduces the pre-Sep-18 behaviour.
    bool budget_fix_h = true;
    /// rev-gfnff stage 3a(ii) (Claude Generated, Sep 18, 2026): the "conserving" share of
    /// FABLE_REVIEW_2 A.5, selected by -gfnff.rev_share_form conserving. Two changes at once,
    /// because neither works without the other (A.5 measured both halves separately):
    ///   (1) f_i = min(1, Val_i / S_i) per ATOM with S_i = sum_k w_ik g_ik, c_ij = f_i f_j.
    ///       The delivered rule is a LEFT-OVER rule, f_i = clip((Val_i - sum_{k != j} w_ik)/w_ij):
    ///       with w ~ 1 on every partner, one partner too many makes every pair of that atom see
    ///       "nothing left", so a 5-coordinate carbon hands out 0 of its 4 valences instead of
    ///       4/5 each. The conserving form satisfies sum_j f_i w_ij = min(Val_i, S_i) exactly.
    ///       The PRODUCT, not the mean: with the mean the rkt06 exchange TS pair gets c = 0.75
    ///       and the path breaks (A.5).
    ///   (2) the excess budget is granted by CHARGE, not by element:
    ///       Val_i = Val_Z + min(G(S_i - Val_Z), X_i), X_i = 0 for H and F, 1 for group 13,
    ///       6 - Val_Z for period >= 3 groups 15-17, and clip(Q_i) otherwise, with Q_i the
    ///       topological (phase-1 EEQ) charge of atom i plus that of its H partners in this
    ///       corner - a per-corner constant, so no chain rule runs through it and a change is
    ///       carried by the existing s-blend. What separates NH4+ (four real bonds) from
    ///       NH3 + H (no bond) is neither the element nor the geometry but the electron count,
    ///       and the only electron-count information a force field has is the charge.
    /// DEFAULT ON since Sep 19, 2026 (operator decision), together with share_donor_rule below -
    /// which is what made the flip possible: the dative/ylide regression that held it back
    /// (+73 to +110 kcal/mol on H3N-BH3, amine oxides, N-ylides) is 0.00 with the donor rule, and
    /// what the mode buys is the artificial radical adducts (class-C dev min -87..-107 -> -1.5..
    /// +0.0 kcal/mol) and the grid's per-step continuity (487 -> 72 events >= 50 kJ/mol over the
    /// 20 react-MD cells). `-gfnff.rev_share_form delivered` reproduces the old behaviour.
    bool share_conserving = true;
    /// rev-gfnff 3a(ii) "conserving": half-width of the C1 smooth min in Val/S units. The min is
    /// exactly 1 above 1 and exactly Val/S below 1 - a, so this only smooths the corner between.
    double share_min_width = 0.1;
    /// rev-gfnff stage 3a(ii) (Claude Generated, Sep 19, 2026): the DONOR RULE of the conserving
    /// share, FABLE_REVIEW_2 A.5's open item and WORK_STATUS 3.8(1). The charge rule above cannot
    /// see a DATIVE bond: a dative bond puts a whole valence into the acceptor's empty orbital,
    /// but the donor's EEQ charge is ~+0.2, not +1, so X_i stays ~0.2 and all four wells of an
    /// amine borane's nitrogen are scaled by ~0.8 - measured as +73 to +110 kcal/mol on H3N-BH3,
    /// amine oxides and N-ylides, where the delivered rule is inert. The rule: atom i is granted
    /// X_i >= 1 if, IN THIS CORNER, it has a partner j that is either (a) a group-13 element
    /// (B, Al, ... - an empty p orbital, which no bond count can show) or (b) carries fewer
    /// partners than its own nominal sigma valence, i.e. a free coordination site (the amine
    /// oxide's one-coordinate O, the ylide's three-coordinate C). It is granted ON TOP of the
    /// charge rule (X_i = max(clip(Q_i), 1)), never below it, and only in the charge branch -
    /// group 13 and the period >= 3 expansion already have a larger cap. Both tests read the
    /// corner's own bond list, so this is a per-corner constant exactly like Q_i: no chain rule,
    /// and a change is carried by the existing s-blend. DEFAULT ON, and only read when
    /// share_conserving is on; -gfnff.rev_share_donor_rule false is the ablation arm.
    bool share_donor_rule = true;
    /// rev-gfnff stage 3a(iii) (Claude Generated, Sep 18, 2026): the BOND WELL FORM.
    /// 0 = gauss (delivered, bit-identical), 1 = MG, 2 = erf-Morse. Both new forms are
    ///     E = -D (2y - y^2)
    /// with y = exp(-(a x + beta x^2)) (MG) or y = erfc((x - u)/sigma)/erfc(-u/sigma)
    /// (erf-Morse), x = r - r0. Both have E(r0) = -D and E'(r0) = 0 identically, and both are
    /// CURVATURE-PINNED: a (MG) resp. u (erf-Morse, by bisection) is chosen so that
    /// E''(r0) = 2 alpha |k_b|, the delivered Gaussian's own force constant. So r_min and the
    /// force constant are reproduced by construction and only the DEPTH (s = D/|k_b|) and the
    /// TAIL are fitted, per element pair, from the class-A reference scans (rev_well_table.h).
    /// The term weight w is NOT applied to these wells: they decay by themselves (the fit
    /// measures |E_pair| at the last grid point as median 0.000, max 0.69 kcal/mol), so
    /// multiplying by w would truncate the tail the fit just put there.
    /// 0 = gauss, 1 = mg, 2 = erfmorse, 3 = mg2, 4 = mg3 (the bond-order-resolved MG well).
    /// DEFAULT 4 (MG3) since Sep 22, 2026 (operator decision; 1 = MG from Sep 19 to Sep 22).
    /// This in-class initializer is inert: GFNFF::setupRevSettings() assigns well_form
    /// unconditionally in the constructor, and the whole rev path is gated on `enabled`
    /// (false here), so plain -method gfnff cannot see this - verified, not assumed.
    int well_form = 4;
    /// rev-gfnff 3b diagnostic: when >= 0, every bond's continuous order is replaced by this
    /// value before the order-resolved well table is read (PARAM rev_well_order_override).
    double well_order_override = -1.0;
    /// rev-gfnff stage 3a(ii) (Sep 2026): "an H is never sp" - an sp hydrogen is not treated as
    /// a bridging atom, so its bond keeps the full strength instead of the reference's 0.30
    /// scaling. See the comment at the rule in gfnff_method.cpp.
    /// SUPERSEDED by h_scope/h_scope_h1 below (FABLE_BOND_STATE_2.md sec 7.6, "Q5"): this flag
    /// stays as its own, narrower, independent mechanism (see the note at h_scope) rather than
    /// being removed or aliased, since it is the only handle left when h_scope_h1 is off.
    bool h_not_sp = true;
    /// rev-gfnff Q5 (Claude Generated, Sep 2026; FABLE_BOND_STATE_2.md sec 7.6, "H-scope"):
    /// master switch for the hydrogen-perception rule set. A hydrogen bridging two partners
    /// (an X-H-Y 3c-4e/3c-2e bond, a geminal migrating H) is today given hyb = 1 ("sp") by
    /// determineHybridizationFortran, which lets it trigger four rules meant for a genuine sp
    /// centre: the bsmat[1][*] bond-strength column (1.3234x instead of the terminal-H
    /// bsmat[hyb_X][0]), 3-ring membership through the H-H/X-H edges of a bridge triangle
    /// (ringf 1.18 plus, for a bridging carbon, the sibling AND bridging C-H fxh 1.05), an
    /// angle centred on the bridging H with theta0 = 180 deg (a bent bridge is then penalised
    /// towards linear), and counting as an sp/sp2 "picon" neighbour for pi-conjugation. All
    /// four are geometry-free per-corner constants (element + partner count + ring/pi
    /// membership of the CORNER's own bond list), so a change is carried by the existing
    /// s-blend exactly like every other corner quantity - no new derivative. Measured
    /// (test_cases/revgfnff/_log/H_SCOPE_IMPL_STATUS.md): FHF- De -120.6 -> -76.9 kcal/mol
    /// (known ~-45; plain gfnff -74.9), CH5+ probe +70.9..+87.8 kcal/mol, rkt06 (14-point
    /// collinear H+H2 path) exactly 0.00 change at every point (H-H bond strength is a pure
    /// Z==1&&Z==1 check, independent of hybridization, and the bridging angle is exactly 180
    /// deg by symmetry there). DEFAULT OFF: bit-identical over GMTKN55+MOR41+S30L-CI (2647
    /// structures) when off; touches exactly the 80 structures with a genuinely 2+-coordinate
    /// hydrogen (FABLE_BOND_STATE_2.md sec 7.6/7.3) when on, none other. Sub-switches
    /// h_scope_h1/h2/r1 below are ablation arms, read only when this is true; a hydrogen never
    /// counting as a picon neighbour (the design's "P1") is a structural CONSEQUENCE of h1 (a
    /// hydrogen with hyb forced to 0 can never satisfy the hyb==1||2 picon test) and has no
    /// separate switch - there is no independent code path to gate.
    bool h_scope = false;
    /// rev-gfnff Q5 "H1": every atom with Z == 1 gets hyb = 0, whatever its partner count (the
    /// grp == 1 branch in determineHybridizationFortran stays reachable for Li/Na/K - this
    /// tests the ELEMENT, not the periodic group). Subsumes h_not_sp for real hydrogen (once
    /// hyb(H) is forced to 0 it can never equal 1, so the is_bridge/0.30-scaling block that
    /// h_not_sp guards becomes structurally unreachable for Z==1, independent of h_not_sp's own
    /// value). DEFAULT true (only read when h_scope is true).
    bool h_scope_h1 = true;
    /// rev-gfnff Q5 "H2": no GFN-FF angle term is ever centred on a Z==1 atom. Required
    /// together with h1 - h1 alone (hyb(H) = 0) gives a spurious tetrahedral theta0 = 109.5 deg
    /// at the bridging H, which was tried and rejected once already (see the note at
    /// determineHybridizationFortran's return). DEFAULT true (only read when h_scope is true).
    bool h_scope_h2 = true;
    /// rev-gfnff Q5 "R1": ring enumeration excludes every Z==1 atom entirely - no ringf and no
    /// 3-ring fxh correction reaches ANY bond of a hydrogen-bridged ring, including the
    /// bridging bonds themselves (FABLE_BOND_STATE_2.md sec 7.6 found this must apply to the
    /// bridging C-H too, not just its siblings: fxh is keyed on the CARBON's ring membership,
    /// so excluding H from ring perception clears it for every C-H of that carbon at once).
    /// DEFAULT true (only read when h_scope is true).
    bool h_scope_r1 = true;
    bool blend = true;            ///< stage 1b: dual-topology blending of the bonded terms over a transition
    double bo_center = 2.0;       ///< term WEIGHT switch: R = f_b (rcov_i + rcov_j) fat_i fat_j (wide: the well decays by itself)
    double bo_width = -7.5;       ///< k of the weight switch (negative: w -> 1 inside R)
    double w_join = 0.05;         ///< the react scan's join weight (rev_bo_form): the TERM weight is (w - w_join)/(1 - w_join) clamped, so a pair joins and leaves the lists at exactly zero weight (Sep 12, 2026: the join used to cost 1-5 kJ/mol, the dominant NVE drift)
    double bo2_center = 1.4;      ///< BOND ORDER switch for E_over (Sep 12, 2026: softened from 1.3/-16, which sat inside the thermal amplitude of X-H bonds and made E_over a stiff wall at 0.5 fs)
    double bo2_width = -6.0;      ///< k of the bond-order switch (softened Sep 12, 2026 for 0.5 fs stability)
    double bo3_center = 1.6;      ///< stage 1b TRANSITION coordinate switch: the neighbour re-parametrisation of a topology change blends in over rev_tr_begin..rev_tr_end of this order, i.e. between 1.63x (1,3 pairs read 0.02) and 1.31x (equilibrium bonds read > 0.95)
    double bo3_width = -8.0;      ///< k of the transition switch (-5 was tried for 0.5 fs and was LESS stable: the wider window lets more transitions overlap)
    double bo4_center = 1.7;      ///< REPULSION BLEND switch: 0.9993 at 1.38x (the turning point of a hot X-H bond; at 1.5x/-10 the 4.5 % of the huge non-bonded repulsion there blew up 0.5 fs MD), 0.5 at 1.7x where a bond really breaks, 0.02 at 1.9x; 1,3 and 1,4 pairs are excluded topologically
    double bo4_width = -12.0;     ///< k of the repulsion blend switch
    double bo5_center = 1.3;      ///< REPULSION BLEND switch of NON-BONDED-list pairs: 0.5 at 1.3x, 1e-4 at 1.63x (1,3 / 1,4 distances), 1e-9 at H-bond distance - a pair that really approaches has joined the bonded list (w > 0.05 at 2.3x) long before, where bo4 takes over continuously (both switches read ~0 there)
    double bo5_width = -12.0;     ///< k of the non-bonded repulsion blend switch
    double over_k = 10.0;         ///< softplus steepness
    double over_shift = 0.5;      ///< penalty argument is (bo_sum - valence - over_shift): no penalty for a saturated atom
    std::vector<double> rcov;     ///< per atom, Bohr (GFNFF covalent radius, as the react scan uses it)
    std::vector<double> fat;      ///< per atom, element scaling of the react scan threshold
    std::vector<double> over_p;   ///< per atom, penalty prefactor p_Z (Eh)
    std::vector<double> valence;  ///< per atom, nominal sigma valence Val_Z

    /// switching radius of the term weight of a pair
    double R(int i, int j) const { return bo_center * (rcov[i] + rcov[j]) * fat[i] * fat[j]; }
    /// switching radius of the bond order of a pair
    double R2(int i, int j) const { return bo2_center * (rcov[i] + rcov[j]) * fat[i] * fat[j]; }
    /// switching radius of the transition coordinate of a pair (stage 1b)
    double R3(int i, int j) const { return bo3_center * (rcov[i] + rcov[j]) * fat[i] * fat[j]; }
    /// switching radius of the repulsion blend of a pair
    double R4(int i, int j) const { return bo4_center * (rcov[i] + rcov[j]) * fat[i] * fat[j]; }
    double R5(int i, int j) const { return bo5_center * (rcov[i] + rcov[j]) * fat[i] * fat[j]; }
};

/**
 * @brief rev-gfnff stage 1b: one topology transition in flight (Claude Generated, Sep 2026)
 *
 * While the term weight w of the pair (i, j) crosses the window [w_a, w_b] the workspace
 * carries the bonded lists of BOTH topologies: the primary (new, after the event) and the
 * alternative (old, before the event), and evaluates
 *     E_bonded = (1 - s) E_alt + s E_primary,   s = clamp((w - w_a) / (w_b - w_a), 0, 1).
 * Forming: w_a = 0.05 -> w_b = 0.5 (the new bond and the re-selected neighbour parameters
 * grow in as the pair approaches); breaking: w_a = 0.5 -> w_b = 0.02 (the old topology fades
 * out while the bond stretches). The list swap therefore costs nothing at either end.
 */
struct RevTransition {
    bool active = false;
    bool forming = true;            ///< formation (s grows as the pair closes) or break (s grows as it opens)
    bool tight = false;             ///< coordinate is the bond-order switch (1,3 ring closure) instead of the transition switch
    bool well_blend = false;        ///< a FORMING pair: its own well lives only in the new corners, so the blend carries it in over s (rev_form_switch = order). false = the well is copied into every corner, which is only energy-neutral when the join sits where the term weight is ~0 (rev_form_switch = weight).
    int i = -1, j = -1;
    double w_a = 0.05, w_b = 0.5;   ///< window of the transition COORDINATE c (bo3 switch): s = (c - w_a) / (w_b - w_a), clamped to [0, 1]
    // per step (set by FFWorkspace::updateBlend)
    double w = 0.0, dwdr = 0.0, r = 0.0, s = 0.0, dsdw = 0.0; ///< w = term weight of the pair (reporting only); dwdr = dc/dr (gradient of s)
    double c = 0.0;                 ///< transition coordinate of the pair (bo3 order)
};

struct PartitionRanges {
    std::pair<int,int> bonds = {0,0};
    std::pair<int,int> angles = {0,0};
    std::pair<int,int> dihedrals = {0,0};
    std::pair<int,int> extra_dihedrals = {0,0};
    std::pair<int,int> inversions = {0,0};
    std::pair<int,int> storsions = {0,0};
    std::pair<int,int> dispersions = {0,0};
    std::pair<int,int> d4_dispersions = {0,0};
    std::pair<int,int> bonded_reps = {0,0};
    std::pair<int,int> nonbonded_reps = {0,0};
    std::pair<int,int> coulombs = {0,0};
    /// Implicit Coulomb (no stored pair list): the ATOM range [first, second) this partition
    /// owns as the outer index i, balanced by pair count. Claude Generated (Sep 2026).
    std::pair<int,int> coulomb_atoms = {0,0};
    std::pair<int,int> hbonds = {0,0};
    std::pair<int,int> xbonds = {0,0};
    std::pair<int,int> atm_triples = {0,0};
    std::pair<int,int> batm_triples = {0,0};
    std::pair<int,int> vdws = {0,0};            // UFF/QMDFF LJ pairs
};

/**
 * @brief Unified force field workspace — shared state + partitioned accumulators
 *
 * Claude Generated (March 2026): Replaces ForceField+ForceFieldThread for GFN-FF.
 *
 * Key differences from ForceFieldThread:
 *   - Interaction lists are shared (ranges on master), not copied per thread
 *   - Only accumulators are per-partition
 *   - T=1 path has zero pool overhead
 *   - postProcess() handles Coulomb TERM 2+3 and dEdcn chain-rule (sequential)
 *
 * Usage:
 *   FFWorkspace ws(num_threads);
 *   ws.setInteractionLists(std::move(params));
 *   ws.setAtomTypes(atoms);
 *   ws.partition();
 *   // Per step:
 *   ws.setGeometry(geom);
 *   ws.setEEQCharges(q);
 *   ws.setD3CN(cn);
 *   double E = ws.calculate(gradient);
 */
class FFWorkspace {
public:
    explicit FFWorkspace(int num_threads = 1);

    // === Init (once after parameter generation) ===

    /// Move interaction lists from GFNFFParameterSet
    void setInteractionLists(GFNFFParameterSet&& params);
    /// Claude Generated (Sep 2026, rev-gfnff): the runtime parameter tables the kernels read their gen scalars from
    void setTables(std::shared_ptr<const GFNFFTables> tables) { m_tables = std::move(tables); }
    const GFNFFTables& tables() const { return *m_tables; }

    /// Set atom types (element numbers)
    void setAtomTypes(const std::vector<int>& atoms);
    /// Claude Generated (Sep 2026, rev-gfnff stage 1): switch the continuous-bond-order kernels on
    void setRev(const RevSettings& rev) { m_rev = rev; }
    const RevSettings& rev() const { return m_rev; }
    /// per-atom bond-order sums of the last calculate() (rev mode; empty otherwise)
    const Vector& revBondOrderSum() const { return m_rev_bo_sum; }
    /// rev-gfnff stage 1b: install `params` as the new primary lists and keep the CURRENT bonded
    /// lists as the alternative topology, blended on the weight of pair (i, j) over [w_a, w_b]
    // ---- rev-gfnff stage 1b (Sep 2026): multi-transition blending over 2^k topology corners.
    // Every transition t in flight has a coordinate s_t in [0, 1]; corner mask b holds the force
    // field of the topology "base + events with bit set"; E = sum_b W_b E_b with
    // W_b = prod_t (b_t ? s_t : 1 - s_t). Each corner is a complete FFWorkspace state (bonded and
    // non-bonded lists, charges, e0) swapped into the evaluation slot in O(1), so no kernel knows
    // about the blend. The slot holds the all-ones corner (the expected topology) between calls.
    struct TopologyState {
        std::vector<Bond> bonds;
        std::vector<Angle> angles;
        std::vector<Dihedral> dihedrals, extra_dihedrals;
        std::vector<Inversion> inversions;
        std::vector<GFNFFSTorsion> storsions;
        std::vector<GFNFFDispersion> dispersions, d4_dispersions;
        std::vector<GFNFFRepulsion> bonded_reps, nonbonded_reps;
        std::vector<GFNFFCoulomb> coulombs;
        std::vector<GFNFFHydrogenBond> hbonds;
        std::vector<GFNFFHalogenBond> xbonds;
        std::vector<ATMTriple> atm_triples;
        std::vector<GFNFFBatmTriple> batm_triples;
        std::vector<BondHBEntry> bond_hb_data;
        std::vector<HBGradEntry> hb_grad_entries;
        std::vector<int> hb_grad_offsets, hb_grad_list;
        std::set<std::pair<int, int>> bonded_pairs;
        std::vector<PartitionRanges> partitions;
        Vector eeq_charges, topology_charges, coul_chi_base, coul_gam, coul_alp, coul_cnf, coul_chi_static;
        std::vector<SqePairData> sqe_pairs; ///< rev-gfnff stage 2: this corner's split-charge pairs
        double e0 = 0.0;
    };
    /// start a transition: params[m] is the parameter set of corner (m | 1 << k) for every existing corner m
    void beginTransition(RevTransition tr, std::vector<GFNFFParameterSet>&& params);
    /// end transition t: keep the corners with its bit set (completed) or cleared (reverted)
    void endTransition(int t, bool keep_new);
    const std::vector<RevTransition>& transitions() const { return m_transitions; }
    bool transitionActive() const { return !m_transitions.empty(); }
    const std::vector<Bond>& bonds() const { return m_bonds; } ///< bond list of the slot corner (fading wells are copied from here)
    /// called with the corner mask after that corner was swapped into the slot and before it is evaluated (per-corner EEQ charges)
    void setCornerPrepare(std::function<void(int)> f) { m_corner_prepare = std::move(f); }
    /// rev-gfnff stage 2 (Sep 2026): install the split-charge pairs of the corner currently in the
    /// slot. Called per step from the same callback that sets the corner's EEQ charges, so the data
    /// travels with the corner through swapState().
    void setSqePairs(std::vector<SqePairData> pairs) { m_sqe_pairs = std::move(pairs); }
    /// rev-gfnff stage 2: pairs with b below this are treated as rigid (no hardness term)
    void setSqeBmin(double b) { m_sqe_bmin = b; }
    /// rev-gfnff stage 2 "B2": functional form of kappa(b) (EEQSolver::SqeKappaForm as an int)
    /// and the exponent of the Power form. MUST be the same values the SQE solve used, otherwise
    /// p is not stationary for the kappa this kernel differentiates and the gradient is wrong.
    void setSqeKappaForm(int form, double exponent) { m_sqe_kappa_form = form; m_sqe_kappa_exponent = exponent; }
    const std::vector<SqePairData>& sqePairs() const { return m_sqe_pairs; }
    void updateTransitions();

    /// Create partition ranges and allocate accumulators
    void partition();

    // === Per-step state updates ===

    /// Set current geometry (Bohr coordinates)
    void setGeometry(const Matrix& geom) { m_geometry = geom; }

    /// Set Phase-2 EEQ charges (geometry-dependent)
    void setEEQCharges(const Vector& q) { m_eeq_charges = q; }

    /// Set Phase-1 topology charges (fixed)
    void setTopologyCharges(const Vector& q) { m_topology_charges = q; }

    /// Set D3 coordination numbers (for dynamic r0)
    void setD3CN(const Vector& cn) { m_d3_cn = cn; }

    /// Set the CN read by the Coulomb self-energy's chi(CN) = chi_base + cnf*sqrt(CN) term.
    /// Claude Generated (Sep 2026): must be refreshed on EVERY evaluation, energy-only calls
    /// included - before, only setCNDerivatives() (gradient calls) wrote m_cn, so an energy-only
    /// call on a reused calculator evaluated the EN term with the CN of the last gradient
    /// geometry (see test_cases/revgfnff/_log/STALE_CN_STATUS.md).
    void setCN(const Vector& cn) { m_cn = cn; }

    /// Visit every D4 pair list - the installed one and each stored rev-gfnff corner's - so
    /// GFNFF can refresh the per-pair C6(CN) every call (stale-CN package, Sep 2026). Only C6
    /// may be changed; the lists themselves stay fixed. Returns the number of lists visited.
    int forEachD4PairList(const std::function<void(std::vector<GFNFFDispersion>&)>& f)
    {
        int n = 0;
        if (!m_d4_dispersions.empty()) { f(m_d4_dispersions); ++n; }
        for (auto& c : m_corners)
            if (!c.d4_dispersions.empty()) { f(c.d4_dispersions); ++n; }
        return n;
    }


    /// Set CN, CNF, and CN derivatives (gradient only)
    /// Claude Generated (WP4, May 2026): dcn now CNDerivStore (pair-list) instead of std::vector<SpMatrix>
    void setCNDerivatives(const Vector& cn, const Vector& cnf,
                          const CNDerivStore& dcn);

    /// Set dc6dcn pointer for D4 dispersion CN gradient
    void setDC6DCNPtr(const Matrix* ptr) { m_dc6dcn_ptr = ptr; }

    /// Set baseline energy (e0 from parameter set)
    void setE0(double e0) { m_e0 = e0; }

    // === Main calculation ===

    /// Calculate total energy (and gradient if requested)
    double calculate(bool gradient);

    // === Results ===

    const GeoGradMatrix& gradient() const { return m_result_gradient; }
    const FFEnergyComponents& energyComponents() const { return m_result_energy; }
    const FFTermTimings& termTimings() const { return m_result_timings; }
    int threadCount() const { return m_num_threads; }
    const Vector& dEdcnTotal() const { return m_dEdcn_total; }
    const Vector& dEdcnBondTotal() const { return m_dEdcn_bond_total; }

    /// Gradient before CN chain-rule (diagnostic, valid after calculate with gradient)
    const GeoGradMatrix& gradientBeforeCN() const { return m_grad_before_cn; }

    // Per-component gradient getters (only valid if store_components=true)
    const GeoGradMatrix& gradientBond() const { return m_result_grad_bond; }
    const GeoGradMatrix& gradientAngle() const { return m_result_grad_angle; }
    const GeoGradMatrix& gradientTorsion() const { return m_result_grad_torsion; }
    const GeoGradMatrix& gradientRepulsion() const { return m_result_grad_repulsion; }
    const GeoGradMatrix& gradientCoulomb() const { return m_result_grad_coulomb; }
    const GeoGradMatrix& gradientDispersion() const { return m_result_grad_dispersion; }
    const GeoGradMatrix& gradientHB() const { return m_result_grad_hb; }
    const GeoGradMatrix& gradientXB() const { return m_result_grad_xb; }
    const GeoGradMatrix& gradientBATM() const { return m_result_grad_batm; }
    const GeoGradMatrix& gradientATM() const { return m_result_grad_atm; }

    // === Configuration ===

    void setStoreGradientComponents(bool v) { m_store_components = v; }
    void setPool(CxxThreadPool* pool) { m_pool = pool; }
    int numThreads() const { return m_num_threads; }

    // Term enable flags (match ForceField flags)
    void setDispersionEnabled(bool v) { m_dispersion_enabled = v; }
    void setHBondEnabled(bool v) { m_hbond_enabled = v; }
    void setRepulsionEnabled(bool v) { m_repulsion_enabled = v; }
    void setCoulombEnabled(bool v) { m_coulomb_enabled = v; }

    // Coulomb self-energy parameters (extracted from pairs at init)
    void setCoulombSelfEnergyParams(const Vector& chi_base, const Vector& gam,
                                     const Vector& alp, const Vector& cnf,
                                     const Vector& chi_static);

    // Bond-HB data for coordination number calculation
    void setBondHBData(const std::vector<BondHBEntry>& data) { m_bond_hb_data = data; }

    // Dynamic HB/XB list updates (for MD simulations)
    /**
     * @brief Full interaction-list swap + re-partition for a live workspace.
     *
     * Claude Generated (Aug 2026): used by the GFN-FF react topology mode when the
     * bond topology changed. Atom types, thread pool and partition count survive;
     * all master lists (incl. the bonded/non-bonded repulsion partition and the
     * Coulomb self-energy params) are replaced by the freshly generated set.
     */
    void rebuildInteractionLists(GFNFFParameterSet&& params);

    void updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds);
    void updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds);

    /**
     * @brief Replace the bonded/non-bonded repulsion pair lists and re-partition.
     *
     * Claude Generated (Sep 2026): the non-bonded repulsion list is built once from a
     * hard 20 Bohr distance cutoff (GFNFF::generateRepulsionPairsNative()) and, unlike
     * HB/XB, was never refreshed during MD — a pair that starts beyond the cutoff and
     * diffuses inside it is never added, so it can pass through the geometric wall with
     * zero repulsive force (see GFNFF::updateNonbondedRepulsionIfNeeded()). Mirrors
     * updateHBonds()/updateXBonds(): cheap vector swap + range re-partition, no other
     * workspace state touched.
     */
    void updateRepulsion(const std::vector<GFNFFRepulsion>& bonded_reps,
                          const std::vector<GFNFFRepulsion>& nonbonded_reps);

    /**
     * @brief Replace the D4 dispersion pair list and re-partition.
     *
     * Claude Generated (Sep 2026): counterpart to updateRepulsion() for the D4 pair list
     * (see GFNFF::updateDispersionPairsIfNeeded()). Takes the list by rvalue: at 7320 atoms
     * it holds ~9.5 M pairs (~0.7 GB), so a copy per rebuild is not affordable.
     */
    void updateD4Dispersions(std::vector<GFNFFDispersion>&& pairs);

    /// Current D4 pair list (read-only; GPU re-upload source after a rebuild)
    const std::vector<GFNFFDispersion>& d4Dispersions() const { return m_d4_dispersions; }

    /**
     * @brief Replace the explicit Coulomb pair list and re-partition.
     *
     * Claude Generated (Sep 2026): only used when an explicit, distance-truncated Coulomb list
     * exists (eeq_distance_cutoff > 0, or coulomb_implicit=false). The implicit path enumerates
     * every pair every step and has no list to go stale. See GFNFF::updateCoulombPairsIfNeeded().
     */
    void updateCoulombPairs(std::vector<GFNFFCoulomb>&& pairs);

    /// Current explicit Coulomb pair list (read-only; GPU re-upload source)
    const std::vector<GFNFFCoulomb>& coulombPairs() const { return m_coulombs; }

    // Master bond list accessor bonds(): declared once above (rev-gfnff stage 1b, slot corner).
    // feature/multi-gpu added an identical one for the bond-HB cross reference; merged Sep 25, 2026.

    /**
     * @brief Replace the bond-HB cross reference after an HB re-detection.
     *
     * Claude Generated (Sep 2026): the HB-modified bond term (egbond_hb, bond.nr_hb >= 1) and
     * the H-bond coordination number (computeHBCoordinationNumbers) read m_bond_hb_data and
     * bond.nr_hb. Both were filled once at setup and never refreshed on the CPU after
     * updateHBXBIfNeeded() replaced the HB list — the GPU path already did this
     * (GFNFF::rebuildBondHBData() -> updateBondHBMetadata()), so the two engines diverged
     * after the first HB re-detection. Mirrors the GPU update.
     *
     * @param nr_hb  One entry per bond in m_bonds order (number of N/O acceptors of its A-H)
     * @param data   The (A, H, B-atoms) entries for all bonds with nr_hb >= 1
     */
    void updateBondHBData(const std::vector<int>& nr_hb, std::vector<BondHBEntry>&& data);

    // Access master interaction list sizes (for diagnostics)
    int bondCount() const { return static_cast<int>(m_bonds.size()); }
    int dispersionPairCount() const { return static_cast<int>(m_dispersions.size() + m_d4_dispersions.size()); }

    int getHBondCount() const { return static_cast<int>(m_hbonds.size()); }
    int getXBondCount() const { return static_cast<int>(m_xbonds.size()); }

    /// rev-gfnff pair-validity gate (Claude Generated, Sep 2026; FABLE_BOND_STATE_2.md sec
    /// 2.1-rev, cap_i / X_i): the FULL per-atom budget cap of the "conserving" valence share,
    /// factored out of prepareConservingShare()'s per-atom loop (ff_workspace_gfnff.cpp) so the
    /// pair-validity gate (GFNFF::findInvalidPairValidityPairs, gfnff_pair_validity.cpp) reads
    /// the IDENTICAL cap formula instead of a second, drifting implementation - the gate runs
    /// inside GFNFF, before any FFWorkspace instance exists for the corner under test, so it
    /// cannot call prepareConservingShare() itself (that needs a populated m_bonds/m_atom_types/
    /// m_topology_charges/m_rev_share_sum). This static function has no such dependency: it is a
    /// pure function of its five scalar arguments plus GFNFFParameters::periodic_group (defined
    /// in ff_workspace_gfnff.cpp, which already includes gfnff_par.h).
    /// @param Z                atomic number
    /// @param qgroup_i         topological (Phase-1 EEQ) charge of the atom plus its H partners
    /// @param donor_i          true if the atom donates into a deficient/group-13 partner's
    ///                         orbital in this corner (the donor rule, rev_share_donor_rule)
    /// @param valz             the atom's nominal sigma valence Val_Z(i) (GFNFF::revValence)
    /// @param fix_h            exclude hydrogen from hypervalent growth (rev_budget_fix_h)
    /// @param delivered_growth [out] true for a d-block metal (FABLE_REVIEW_2 A.5 does not
    ///                         budget metals; the caller reads is_metal instead of a finite cap)
    static double shareCapForAtom(int Z, double qgroup_i, bool donor_i, double valz, bool fix_h,
                                  bool& delivered_growth);

private:
    int m_natoms = 0;
    int m_num_threads = 1;
    CxxThreadPool* m_pool = nullptr;

    // === Shared state (read-only per step) ===
    GeoGradMatrix m_geometry;  // WP-G: RowMajor for contiguous row(i) reads in inner loops
    std::vector<int> m_atom_types;
    std::shared_ptr<const GFNFFTables> m_tables = GFNFFTables::defaults(); ///< runtime tables (Sep 2026)
    RevSettings m_rev;         ///< rev-gfnff stage 1 (Sep 2026)
    Vector m_rev_bo_sum;       ///< per-atom sum_j b_ij BO_ij of the last step (rev mode)
    std::vector<SqePairData> m_sqe_pairs; ///< rev-gfnff stage 2: split-charge pairs of the slot corner
    double m_sqe_bmin = 1e-3;  ///< rev-gfnff stage 2: bond-order floor of the hardness term
    int m_sqe_kappa_form = 0;         ///< rev-gfnff stage 2 "B2": EEQSolver::SqeKappaForm as an int
    double m_sqe_kappa_exponent = 3.0;///< rev-gfnff stage 2 "B2": n of the Power form
    std::vector<RevTransition> m_transitions;      ///< stage 1b: transitions in flight (bit t of a corner mask)
    std::vector<TopologyState> m_corners;          ///< corner states by mask; the slot's own entry is a placeholder
    int m_slot_mask = 0;                           ///< which corner the member slot currently holds
    std::function<void(int)> m_corner_prepare;
    std::vector<double> m_corner_energy;
    std::vector<GeoGradMatrix> m_corner_gradient;
    std::vector<FFEnergyComponents> m_corner_components;
    void swapState(TopologyState& st);             ///< O(1) exchange of the slot with a stored corner
    double calculateSingle(bool gradient);         ///< one topology (the pre-stage-1b calculate())
    Vector m_eeq_charges, m_topology_charges, m_d3_cn;
    /// rev-gfnff stage 3a(i): GFN-FF CN radii in Bohr (GFNFFParameters::gfnff_cn_rcov_bohr), set once
    std::vector<double> m_rev_cn_rcov;
    Vector m_cn, m_cnf;
    CNDerivStore m_dcn;  // Claude Generated (WP4, May 2026): pair-list replaces std::vector<SpMatrix>
    const Matrix* m_dc6dcn_ptr = nullptr;
    double m_e0 = 0.0;

    // Coulomb self-energy parameters (O(N), extracted at init)
    Vector m_coul_chi_base, m_coul_gam, m_coul_alp, m_coul_cnf, m_coul_chi_static;
    /// Claude Generated (Sep 2026): evaluate the N^2/2 Coulomb pairs on the fly from the per-atom
    /// data (charges + alpeeq) instead of reading a stored pair list. Set from the parameter set
    /// (GFN-FF only); see calcCoulomb().
    bool   m_coulomb_implicit = false;
    double m_coulomb_implicit_rcut = std::numeric_limits<double>::infinity();  // no cutoff, as the reference (Sep 2026)

    // Term-enable flags
    bool m_dispersion_enabled = true;
    bool m_hbond_enabled = true;
    bool m_repulsion_enabled = true;
    bool m_coulomb_enabled = true;

    // === Method type (set by setInteractionLists) ===
    FFMethodType m_method_type = FFMethodType::GFN_FF;
    double m_au = 1.0;  ///< Distance unit factor: 1.0 (GFN-FF/Bohr), 1.889726125 (UFF/QMDFF Å→Bohr)

    // === Master interaction lists (owned, moved from ForceFieldParameterSet) ===
    std::vector<Bond> m_bonds;
    std::vector<Angle> m_angles;
    std::vector<Dihedral> m_dihedrals, m_extra_dihedrals;
    std::vector<Inversion> m_inversions;
    std::vector<GFNFFSTorsion> m_storsions;
    std::vector<GFNFFDispersion> m_dispersions, m_d4_dispersions;
    std::vector<GFNFFRepulsion> m_bonded_reps, m_nonbonded_reps;
    std::vector<GFNFFCoulomb> m_coulombs;
    std::vector<GFNFFHydrogenBond> m_hbonds;
    std::vector<GFNFFHalogenBond> m_xbonds;
    std::vector<ATMTriple> m_atm_triples;
    std::vector<GFNFFBatmTriple> m_batm_triples;
    std::vector<BondHBEntry> m_bond_hb_data;
    std::vector<HBGradEntry> m_hb_grad_entries;
    // CSR index of m_hb_grad_entries by H atom, insertion order preserved (B3, Sep 2026):
    // calcBonds() visits only the entries of its own H instead of scanning all of them.
    std::vector<int> m_hb_grad_offsets, m_hb_grad_list;
    std::vector<vdW> m_vdws;                    ///< UFF/QMDFF LJ non-bonded pairs

    // Cached bonded pairs for fast repulsion lookup
    std::set<std::pair<int,int>> m_bonded_pairs;

    // === Partitions + Accumulators ===
    std::vector<PartitionRanges> m_partitions;
    std::vector<FFAccumulator> m_accumulators;

    // === Result storage ===
    GeoGradMatrix m_result_gradient;  // WP-G: RowMajor
    FFEnergyComponents m_result_energy;
    FFTermTimings m_result_timings;  // Claude Generated (May 2026): aggregated per-term timing
    Vector m_dEdcn_total, m_dEdcn_bond_total;
    /// rev-gfnff 3a(ii): sum_k w_ik per atom - the CLAIM of i's other partners - plus the reduced
    /// dE/d(sum_i w_ik) coefficient of the Lambda chain rule
    Vector m_rev_share_sum;
    /// rev-gfnff 3a(ii): the reduced dE/d(sum_i w_ik) coefficient of the Lambda chain rule
    Vector m_dEdshare_total;
    /// rev-gfnff 3a(ii): the effective valence of the corner and its derivative w.r.t. the
    /// settled count. val_i = Val_Z + G(settled_i - Val_Z) (see shareExcess); dval_i = dG/dx.
    /// Both are consumed by calcBonds, which carries their whole chain rule locally.
    Vector m_rev_share_val, m_rev_share_dval;
    /// rev-gfnff 3a(ii) "conserving" share (Claude Generated, Sep 18, 2026, FABLE_REVIEW_2 A.5).
    /// f_i = smoothmin(1, Val_i / S_i) is a PER-ATOM quantity there (the delivered rule's f is
    /// per pair), so it is computed once per corner next to the sums and read by calcBonds.
    /// m_rev_share_dfdS(i) = d f_i / d S_i, which already contains d Val_i / d S_i - the excess
    /// budget of this mode is built from the same wide sum S_i, not from the settled count, so
    /// the whole chain rule rides the existing Lambda pass over the term weights and
    /// m_rev_share_dval is exactly 0 in this mode.
    Vector m_rev_share_f, m_rev_share_dfdS;
    /// rev-gfnff 3a(ii) "conserving": the per-atom excess budget cap X_i (element/charge rule).
    /// Diagnostic only after prepareValenceShare; kept for the CURCUMA_SHAREDUMP table.
    Vector m_rev_share_cap;
    /// rev-gfnff 3a(iii): per BOND (m_bonds order) the well-form parameters of this corner.
    /// D > 0 = the depth, p1 = a (MG) or u (erf-Morse), p2 = beta (MG) or sigma (erf-Morse),
    /// all in ATOMIC units. form = 0 means "no table entry, use the Gaussian" for that bond.
    /// rev-gfnff 3a(iii): one bond's well parameters. form 1 = MG (p1 = a, p2 = beta),
    /// 2 = erf-Morse (p1 = u, p2 = sigma), 3 = the free-curvature MG of stage 3a(iii) step 2 /
    /// stage 3b (same p1/p2 meaning, plus dr0). dr0 is the fitted r0 OFFSET in Bohr, added to the
    /// model's own (dynamic) r0 before the well coordinate x = r - r0 - dr0 is formed; it is 0
    /// for forms 1 and 2, which inherit r0 unchanged.
    struct RevWellPar { double D = 0.0, p1 = 0.0, p2 = 0.0; int form = 0; double dr0 = 0.0;
        /// rev-gfnff P3 (Sep 23, 2026): weight of the UNCAPPED inner branch (0 = the y = 2 cap as
        /// always; 1 = a half-order row's full Morse-type wall). RevWellTableV2::halfOrderWeight.
        double uncap = 0.0; };
    std::vector<RevWellPar> m_rev_well;
    /// rev-gfnff 3a(iii): the (fc, exponent, z_i, z_j) the current m_rev_well was built from.
    /// prepareWellForms() runs on every energy call (next to the share's own pass), but its
    /// inputs are per-bond CONSTANTS, so recomputing the erf-Morse bisection every step is pure
    /// cost - measured at 1.4x the whole react-MD wall time before this cache.
    /// (fc, exponent, z_i, z_j, rev_order) - rev_order is in the stamp because the
    /// bond-order-resolved table (form 3 / 'mg3') keys on it, and a rebuild can change it while
    /// leaving everything else alone. Claude Generated (Sep 20, 2026).
    /// + rev_pi_excess (Sep 25, 2026, P3 pi* prototype): it selects the pi-excess row.
    std::vector<std::array<double, 6>> m_rev_well_stamp;
    int m_rev_well_stamp_form = -1;
    /// rev-gfnff 3a(ii) (Claude Generated, Sep 14, 2026): the SMOOTH 1,3 proxy - the genuineness
    /// g_p of every pair of the corner's bond list, in m_bonds order. A compact polyhedron's
    /// perception carries 1,3 contacts as bonds - the six F...F contacts of a tetrahedral BF4-
    /// (1.867 A, tight bond order 0.4739) are the clearest case - and the wide term weight reads
    /// ~1 for them, so without a mask they claim the full valence of both ends and halve every
    /// genuine bond (BF4- at B-F = 1.143 A: +569.7 kcal/mol). The discrete test "the two ends
    /// share a BONDED neighbour" separates them from the migrating pair of an exchange transition
    /// state (tight bond order 0.4985, every continuous quantity within 5 % - see
    /// test_cases/revgfnff/_log/CIJ_STATUS.md), but a topology test is a switch, so this is its
    /// continuous stand-in: the BOND-ORDER LEAK
    ///     t_p = sum_{k != i,j} sigma_ik sigma_jk,   sigma = shareSettled(b) in [0, 1],
    ///     g_p = shareClip(1 - t_p)                  (1 = genuine bond, 0 = pure 1,3 contact),
    /// i.e. how much the bond order of i's and j's SHARED settled partner leaks onto this pair.
    /// sigma is the settled weight of the tight switch, so a 1,3 contact of two ends that are only
    /// contacts themselves (sigma = 0, e.g. B in BF4-) does not create a leak: the six F...F pairs
    /// read t = sigma_BF^2 = 1 exactly while the four B-F bonds read t = 0 exactly. Nothing is
    /// counted and no threshold is applied - t is a smooth function of the existing bond orders,
    /// C1 in every coordinate (shareClip is C1 and exactly 0/1 outside [0,1]).
    /// Consumed as: the claim of a pair on its ends' valence is w_p g_p (so a 1,3 contact consumes
    /// nothing), and its own well is multiplied by c_p = 1 - g_p (1 - (f_i + f_j)/2), which reads
    /// exactly 1 for a 1,3 contact and reduces to the plain share for a genuine bond. A pair with
    /// g = 0 keeps its own term, exactly as the over-coordination term keeps its 1,3 repulsion.
    std::vector<double> m_rev_share_g;      ///< g_p (per pair, m_bonds order)
    std::vector<double> m_rev_share_dg;     ///< dE/dg_p, accumulated by calcBonds
    std::vector<double> m_rev_share_gclip;  ///< dg_p/d t_p = -shareClipD(1 - t_p)
    std::vector<double> m_rev_share_sig;    ///< sigma_p = shareSettled(b_p)
    std::vector<double> m_rev_share_dsig;   ///< d sigma_p / d r_p (the three-body chain rule)
    std::vector<int> m_rev_share_stamp;     ///< scratch for the common-neighbour lookup
    std::vector<std::vector<int>> m_rev_adj; ///< per atom, the bond indices of the corner (m_bonds order)
    /// the other end of bond q as seen from atom a (rev-gfnff 3a(ii) 1,3 proxy)
    int otherEndOf(int q, int a) const { return (m_bonds[q].i == a) ? m_bonds[q].j : m_bonds[q].i; }
    GeoGradMatrix m_grad_before_cn;  ///< Gradient snapshot before CN chain-rule (diagnostic) — WP-G: RowMajor
    bool m_store_components = false;
    bool m_do_gradient = false;

    // Per-component result gradients
    // WP-G: RowMajor result gradient buffers
    GeoGradMatrix m_result_grad_bond, m_result_grad_angle, m_result_grad_torsion;
    GeoGradMatrix m_result_grad_repulsion, m_result_grad_coulomb, m_result_grad_dispersion;
    GeoGradMatrix m_result_grad_hb, m_result_grad_xb, m_result_grad_batm, m_result_grad_atm;

    // === Core execution ===
    void executeGFNFF(int partition);
    void executeUFF(int partition);    ///< Claude Generated (March 2026): UFF energy/gradient
    void executeCG(int partition);     ///< Sep 2026: coarse-grained LJ spheres/ellipsoids (ff_workspace_cg.cpp)
    void calcCGPairs(int partition);
    void executeQMDFF(int partition);  ///< Claude Generated (March 2026): QMDFF energy/gradient
    void postProcess(bool gradient);
    void reduce();

    // === GFN-FF energy term calculators (ported from ForceFieldThread) ===
    // The six bonded kernels take the accumulator, the list and the index range explicitly so
    // that the same code evaluates the primary and the alternative topology of a blend
    // (rev-gfnff stage 1b, Sep 2026).
    void calcBonds(FFAccumulator& acc, const std::vector<Bond>& list, std::pair<int, int> range);
    void calcAngles(FFAccumulator& acc, const std::vector<Angle>& list, std::pair<int, int> range);
    void calcDihedrals(FFAccumulator& acc, const std::vector<Dihedral>& list, std::pair<int, int> range);
    void calcExtraTorsions(FFAccumulator& acc, const std::vector<Dihedral>& list, std::pair<int, int> range);
    void calcInversions(FFAccumulator& acc, const std::vector<Inversion>& list, std::pair<int, int> range);
    void calcSTorsions(FFAccumulator& acc, const std::vector<GFNFFSTorsion>& list, std::pair<int, int> range);
    /// all six bonded kernels of one topology into `acc`
    void runBonded(FFAccumulator& acc, const PartitionRanges& pr, FFTermTimings* timings);
    /// rev-gfnff: continuous bond order of a pair at distance r (Bohr); dw = db/dr
    double revWeight(int i, int j, double r, double* dw) const {
        double d = 0.0;
        const double w = RevGFNFF::bondOrder(r, m_rev.R(i, j), m_rev.bo_width, &d);
        const double denom = 1.0 - m_rev.w_join;
        if (w <= m_rev.w_join || denom <= 0.0) {
            if (dw) *dw = 0.0;
            return 0.0;
        }
        if (dw) *dw = d / denom;
        return (w - m_rev.w_join) / denom;
    }
    /// the raw (unshifted) term-weight switch, as the react scan sees it
    double revWeightRaw(int i, int j, double r, double* dw) const {
        return RevGFNFF::bondOrder(r, m_rev.R(i, j), m_rev.bo_width, dw);
    }
    /// rev-gfnff: bond order of a pair (tight switch: E_over sums and the repulsion blend)
    double revOrder(int i, int j, double r, double* db) const {
        return RevGFNFF::bondOrder(r, m_rev.R2(i, j), m_rev.bo2_width, db);
    }
    /// rev-gfnff: repulsion blend weight of a pair (see RevSettings::bo4_center)
    double revBlend(int i, int j, double r, double* db) const {
        return RevGFNFF::bondOrder(r, m_rev.R4(i, j), m_rev.bo4_width, db);
    }
    /// rev-gfnff: repulsion blend weight of a NON-BONDED-list pair (tight switch, see RevSettings::bo5_center)
    double revBlendNB(int i, int j, double r, double* db) const {
        return RevGFNFF::bondOrder(r, m_rev.R5(i, j), m_rev.bo5_width, db);
    }
    /// rev-gfnff stage 1b: transition coordinate of a pair (medium switch, see RevSettings::bo3_center)
    double revCoord(int i, int j, double r, double* dc) const {
        return RevGFNFF::bondOrder(r, m_rev.R3(i, j), m_rev.bo3_width, dc);
    }
    /// rev-gfnff 3a(ii): C1 soft clip on [0, 1] — exactly 0 below and exactly 1 above, so an
    /// equilibrium bond (f > 1) is multiplied by exactly 1.0 and the term stays bit-identical.
    /// No free width: the smoothstep runs over the whole unit interval of f.
    static double shareClip(double x)
    {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;
        return x * x * (3.0 - 2.0 * x);
    }
    /// d shareClip / dx (0 at both ends, so the force is continuous)
    static double shareClipD(double x)
    {
        if (x <= 0.0 || x >= 1.0) return 0.0;
        return 6.0 * x * (1.0 - x);
    }
    /// rev-gfnff 3a(ii) (Claude Generated, Sep 2026): smooth excess valence - the softplus
    ///     G(x) = ln(1 + e^(beta x)) / beta,   0 for x <= 0, x for x >> 1/beta.
    /// The effective valence of an atom is Val_Z + G(sum_k sigma(b_ik) - Val_Z) with sigma the
    /// settled-partner weight (shareClip of the TIGHT bond order b, the same switch the
    /// over-coordination term sums). While an atom carries at most Val_Z partners the sum stays
    /// at or below Val_Z, G is 0 and the nominal valence is untouched - which is every ordinary
    /// equilibrium, including a bond stretched inside a molecule (fewer bonds, not more). G only
    /// grows once an atom has MORE partners than its nominal valence, which is also the only case
    /// in which the share can bite at all (f_i >= 1 follows from n_i <= Val_Z). It is what makes
    /// the valence hypervalent-correct: the fourth bond of an ammonium, the third of a hydronium
    /// and the fourth of a perchlorate are settled partners and are credited instead of halving
    /// the genuine bonds. beta is the one new GLOBAL constant of stage 3a(ii) - it sets how
    /// sharply the settled count crosses the nominal valence - and no element data is added.
    static constexpr double kShareExcessBeta = 50.0;
    static double shareExcess(double x)
    {
        const double y = kShareExcessBeta * x;
        if (y > 700.0) return x;                       // exp() overflow guard: G -> x
        return std::log1p(std::exp(y)) / kShareExcessBeta;
    }
    /// d G / dx = sigmoid(beta x) (0 and 1 at the ends, no kink anywhere)
    static double shareExcessD(double x)
    {
        const double y = kShareExcessBeta * x;
        if (y > 700.0) return 1.0;
        if (y < -700.0) return 0.0;
        return 1.0 / (1.0 + std::exp(-y));
    }
    /// rev-gfnff 3a(ii): the settled weight of ONE partner - the same unit-interval smoothstep,
    /// applied to the RESCALED tight bond order 2b - 1. "Settled" therefore means "past the middle
    /// of this pair's switching radius" (b > 1/2, i.e. r < R2), and the transition is a wide,
    /// C1-continuous ramp (b runs 0.05..0.95 over 0.39 R2), not a threshold. Rescaling matters:
    /// shareClip(b) itself credits a HALF-formed bond with half a valence, so at a 5-coordinate
    /// exchange carbon (four full bonds + one half-formed) the settled count already exceeds the
    /// nominal valence and the effective valence rises with it, cancelling the share exactly where
    /// it is needed. With 2b - 1 the four full bonds still read 0.9997 each while the half-formed
    /// one reads 0, so N_i stays below Val_Z and the share keeps its full strength. No new
    /// parameter: same smoothstep, rescaled argument.
    static double shareSettled(double b_order) { return shareClip(2.0 * b_order - 1.0); }
    /// d shareSettled / db (chain rule of the 2b - 1 rescaling)
    static double shareSettledD(double b_order) { return 2.0 * shareClipD(2.0 * b_order - 1.0); }
    /// rev-gfnff 3a(ii) "conserving" share (Claude Generated, Sep 18, 2026): a C1 smooth min(1, x)
    /// that is EXACTLY 1 for x >= 1 and EXACTLY x for x <= 1 - a, with a cubic joining the two.
    /// The exactness on both sides is the point, not a nicety: an equilibrium atom has
    /// Val_i >= S_i, i.e. x >= 1, so its share factor must come out as the literal 1.0 for the
    /// term to stay bit-identical; and an over-claimed atom must get the literal Val_i/S_i for
    /// sum_j f_i w_ij = min(Val_i, S_i) - the valence conservation the mode exists for - to hold.
    /// A softplus-based min satisfies neither (it is off by ln(2)/beta at x = 1).
    /// h(u) = -u^3/a^2 - 2u^2/a on [-a, 0] is the unique cubic with h(-a) = -a, h'(-a) = 1,
    /// h(0) = h'(0) = 0; it is monotone there (h' = -(u/a)(3u/a + 4) > 0).
    static double shareMinOne(double x, double a)
    {
        const double u = x - 1.0;
        if (u >= 0.0) return 1.0;
        if (u <= -a) return x;
        return 1.0 - u * u * u / (a * a) - 2.0 * u * u / a;
    }
    /// d shareMinOne / dx
    static double shareMinOneD(double x, double a)
    {
        const double u = x - 1.0;
        if (u >= 0.0) return 0.0;
        if (u <= -a) return 1.0;
        return -3.0 * u * u / (a * a) - 4.0 * u / a;
    }
    /// rev-gfnff 3a(ii): per-atom sum_k m_ik w_ik of the corner's bond list (main thread, per step)
    void prepareValenceShare();
    /// rev-gfnff 3a(ii) "conserving": per-atom budget cap, effective valence, f_i and df_i/dS_i
    void prepareConservingShare(bool fix_h);
    /// rev-gfnff 3a(iii): per-bond well-form parameters of the corner (D, a/u, tail), main thread
    void prepareWellForms();
    /// rev-gfnff 3a(iii): the erf-Morse offset u from 2 D h(-u/sigma)^2/sigma^2 = K (bisection)
    static double wellErfMorseU(double D, double sigma, double K);
    /// rev-gfnff 3a(ii): the chain rule of that sum and of the effective valence (after the partitions)
    void applyValenceShareGradient();
    /// rev-gfnff: over-coordination energy + gradient, main thread, after the partitions (Sep 2026)
    void calcOverCoordination(bool gradient);
    /// rev-gfnff stage 2: bond-hardness energy + gradient of the split charges (Sep 2026)
    void calcSqeHardness(bool gradient);
    void calcDispersion(int p);
    void calcD4Dispersion(int p);
    void calcBondedRepulsion(int p);
    void calcNonbondedRepulsion(int p);
    void calcCoulomb(int p);
    void calcHydrogenBonds(int p);
    void calcHalogenBonds(int p);
    void calcATM(int p);
    void calcATMGradient(int p);
    void calcBATM(int p);
    void computeHBCoordinationNumbers(int p);

    // === UFF/QMDFF energy term calculators (Claude Generated March 2026) ===
    void calcUFFBonds(int p);
    void calcUFFAngles(int p);
    void calcUFFDihedrals(int p);
    void calcUFFInversions(int p);
    void calcUFFvdW(int p);
    void calcQMDFFBonds(int p);
    void calcQMDFFAngles(int p);

    // === Helpers ===

    /// Linear partition range for T threads
    static std::pair<int,int> linearRange(int total, int t, int T) {
        int begin = t * total / T;
        int end = (t + 1) * total / T;
        return {begin, end};
    }
};
