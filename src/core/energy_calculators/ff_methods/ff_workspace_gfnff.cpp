/*
 * <FFWorkspace GFN-FF Energy Term Calculators>
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
 * Claude Generated (March 2026): GFN-FF energy term calculators ported from
 * ForceFieldThread to FFWorkspace. Physics formulas are identical — only the
 * data access pattern changes (ranges on shared master lists, per-partition
 * accumulators instead of per-thread copies).
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 * Reference: Fortran gfnff_engrad.F90 (energy and gradient formulas)
 */

#include "ff_workspace.h"
#include "cn_calculator.h"
#include "gfnff_par.h"
#include "rev_well_table.h"   // rev-gfnff 3a(iii): AI-fitted per-element-pair well parameters
#include "forcefieldfunctions.h"
#include "gfnff_geometry.h"
#include "src/core/units.h"
#include "src/core/curcuma_logger.h"
#include "src/core/math_compat.h"

#include <fmt/core.h>
#include <cstdlib>
#include <fmt/format.h>

#include <cmath>
#include <unordered_map>

// ============================================================================
// HB Coordination Numbers (must run before bonds, only on partition 0)
// Reference: Fortran gfnff_engrad.F90:361, gfnff_data_types.f90:88
// ============================================================================

void FFWorkspace::computeHBCoordinationNumbers(int p)
{
    // Only partition 0 runs this — it modifies the shared m_bonds array,
    // but each bond.hb_cn_H is written by only one entry, so no race.
    if (p != 0) return;
    if (m_bond_hb_data.empty()) return;

    static const std::vector<double>& rcov_base = GFNFFParameters::covalent_rad_d3;
    constexpr double rcov_43 = 4.0 / 3.0;
    constexpr double kn = 27.5;
    constexpr double rcov_scal = 1.78;
    constexpr double thr = 900.0;

    std::unordered_map<int, double> hb_cn_map;
    m_hb_grad_entries.clear();
    constexpr double inv_sqrt_pi = 0.5641895835477563;

    for (const auto& entry : m_bond_hb_data) {
        int H = entry.H;
        int ati = m_atom_types[H];

        for (int B : entry.B_atoms) {
            int atj = m_atom_types[B];
            double dx = m_geometry(B, 0) - m_geometry(H, 0);
            double dy = m_geometry(B, 1) - m_geometry(H, 1);
            double dz = m_geometry(B, 2) - m_geometry(H, 2);
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > thr) continue;
            double r = std::sqrt(r2);

            double rcovij = rcov_scal * rcov_43 * (rcov_base[ati - 1] + rcov_base[atj - 1]);
            double arg = -kn * (r - rcovij) / rcovij;
            double tmp = 0.5 * (1.0 + curcuma_erf(arg));
            hb_cn_map[H] += tmp;

            double dCN_dr = inv_sqrt_pi * (-kn / rcovij) * std::exp(-arg * arg) / r;
            Eigen::Vector3d r_HB(dx, dy, dz);
            m_hb_grad_entries.push_back({H, B, -dCN_dr * r_HB, dCN_dr * r_HB});
        }
    }

    for (auto& bond : m_bonds) {
        if (bond.nr_hb < 1) continue;
        int H = -1;
        if (m_atom_types[bond.i] == 1) H = bond.i;
        else if (m_atom_types[bond.j] == 1) H = bond.j;
        if (H >= 0) {
            auto it = hb_cn_map.find(H);
            bond.hb_cn_H = (it != hb_cn_map.end()) ? it->second : 0.0;
        }
    }

    // CSR index by H atom. The per-H entry order equals the order the former linear
    // scan visited them, so the gradient accumulation order (and result) is unchanged.
    const int nat = static_cast<int>(m_atom_types.size());
    m_hb_grad_offsets.assign(nat + 1, 0);
    for (const auto& e : m_hb_grad_entries) ++m_hb_grad_offsets[e.H_atom + 1];
    for (int a = 0; a < nat; ++a) m_hb_grad_offsets[a + 1] += m_hb_grad_offsets[a];
    m_hb_grad_list.resize(m_hb_grad_entries.size());
    std::vector<int> fill(m_hb_grad_offsets.begin(), m_hb_grad_offsets.end() - 1);
    for (int k = 0; k < static_cast<int>(m_hb_grad_entries.size()); ++k)
        m_hb_grad_list[fill[m_hb_grad_entries[k].H_atom]++] = k;
}

// ============================================================================
// Bond Stretching (GFN-FF exponential potential)
// Reference: Fortran gfnff_engrad.F90:675-721
// ============================================================================

void FFWorkspace::calcBonds(FFAccumulator& acc, const std::vector<Bond>& list, std::pair<int, int> range)
{
    auto [begin, end] = range;
    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    bool use_dynamic_r0 = (m_d3_cn.size() > 0);
    // Claude Generated (Sep 14, 2026): the per-pair well-depth dump, same env gate as the share
    // table of prepareValenceShare(). Read once per process.
    static const bool s_share_dump = [] {
        const char* d = std::getenv("CURCUMA_SHAREDUMP");
        return d && d[0] == '1';
    }();
    // rev-gfnff stage 3a(i) (Claude Generated, Sep 2026): take the pair's OWN erf-CN
    // contribution out of the CN that builds its r0 (see the block below). rev-only: with
    // rev disabled the r0 is bit-identical to before.
    const bool rev_cn_pair = m_rev.enabled && use_dynamic_r0
                             && m_rev_cn_rcov.size() == static_cast<size_t>(m_natoms);
    // rev-gfnff stage 3a(ii) (Claude Generated, Sep 2026): the valence share factor c_ij. It
    // needs the per-atom bond-order sums of THIS corner, built by prepareValenceShare() on the
    // main thread (see the block in calcBonds below).
    const bool rev_share = m_rev.enabled && m_rev.valence_share
                           && m_rev_share_sum.size() == m_natoms
                           && m_rev_share_val.size() == m_natoms
                           && static_cast<int>(m_rev.valence.size()) == m_natoms;
    // rev-gfnff 3a(ii) (Sep 18, 2026): the conserving share reads the PER-ATOM f built by
    // prepareConservingShare; if that vector is missing the delivered branch is taken.
    const bool rev_share_conserving = rev_share && m_rev.share_conserving
                                      && m_rev_share_f.size() == m_natoms
                                      && m_rev_share_dfdS.size() == m_natoms;
    // rev-gfnff stage 3a(iii) (Sep 18, 2026): the per-bond well-form parameters of this corner,
    // built by prepareWellForms() on the main thread (the erf-Morse bisection must not run per
    // energy call). Bonds whose element pair has no class-A entry keep form = 0 = the Gaussian.
    const bool rev_well_on = m_rev.enabled && m_rev.well_form != 0
                             && m_rev_well.size() == list.size();

    for (int idx = begin; idx < end; ++idx) {
        const auto& bond = list[idx];

        Eigen::VectorXd vi = m_geometry.row(bond.i);
        Eigen::VectorXd vj = m_geometry.row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(vi, vj, derivate, m_do_gradient);

        double dcpair_dr = 0.0;   // d cn_ij / dr, needed by the gradient below
        double r0_ij;
        if (use_dynamic_r0 && bond.z_i > 0 && bond.z_j > 0 &&
            bond.i < m_d3_cn.size() && bond.j < m_d3_cn.size()) {
            double cn_i = m_d3_cn(bond.i);
            double cn_j = m_d3_cn(bond.j);
            if (rev_cn_pair) {
                // rev-gfnff stage 3a(i): r0 = (r0_base + cnfak*CN)*ff is built from
                //     CN_i' = CN_i - cn_ij(r) + 1,
                // i.e. the partner counts as PRESENT even while the pair stretches, where
                // cn_ij is the pair's own contribution to the raw CN
                //     cn_ij(r) = 0.5*(1 + erf(kn*(r - R)/R)),  R = rcov_i + rcov_j  (Bohr)
                // (CNCalculator::pairCNContribution — the same expression the CN itself is built
                // from). Without it, atom i loses j's count as the bond opens, r0 shrinks and the
                // well retreats from the departing atom: a positive feedback with no physical
                // counterpart (a real C-H bond SHORTENS by ~0.012 A when stretched), worth
                // +9..+22 kcal/mol on every X-H bond at 1.3-1.4 r_eq. At r_eq cn_ij is 0.97-0.99,
                // so the correction is +0.01..0.03 in CN and r0 barely moves. This is the
                // continuous form of the frozen-CN (gfnff-fast) reference: it removes only the
                // pair's own fading, not the neighbour re-parametrisation a real topology change
                // produces. 0 new parameters; smooth in r (no branch).
                const double rcov_sum = m_rev_cn_rcov[bond.i] + m_rev_cn_rcov[bond.j];
                if (rcov_sum > 0.0) {
                    const double cpair = CNCalculator::pairCNContribution(rij, rcov_sum);
                    cn_i += 1.0 - cpair;
                    cn_j += 1.0 - cpair;
                    if (m_do_gradient)
                        dcpair_dr = CNCalculator::pairCNContributionDerivative(rij, rcov_sum);
                }
            }
            double ra = bond.r0_base_i + bond.cnfak_i * cn_i;
            double rb = bond.r0_base_j + bond.cnfak_j * cn_j;
            r0_ij = (ra + rb + bond.rabshift) * bond.ff;
        } else {
            r0_ij = bond.r0_ij;
        }

        double dr = rij - r0_ij;
        double alpha_orig = bond.exponent;
        double alpha = alpha_orig;
        double k_b = bond.fc;

        // egbond_hb: Modified exponent for HB X-H bonds
        // Reference: Fortran gfnff_engrad.F90:957-958
        if (bond.nr_hb >= 1) {
            constexpr double VBOND_SCALE = 0.9;
            double t1 = 1.0 - VBOND_SCALE;
            alpha = (-t1 * bond.hb_cn_H + 1.0) * alpha_orig;
        }

        double exp_term = std::exp(-alpha * dr * dr);
        double energy = k_b * exp_term;
        // rev-gfnff stage 1 (Sep 2026): the well is switched off smoothly by the continuous
        // bond order instead of being dropped when the bond leaves the topology.
        double w = 1.0, dwdr = 0.0;
        if (m_rev.enabled && m_rev.bond_weight)
            w = revWeight(bond.i, bond.j, rij, &dwdr);
        // rev-gfnff stage 3a(iii) (Claude Generated, Sep 18, 2026): the two alternative well
        // forms. `well` is the pair's well BEFORE the term weight and the valence share;
        // `dwell_dx` its derivative in x = r - r0. The new forms are NOT multiplied by w (they
        // decay by themselves - see RevSettings::well_form), the delivered Gaussian is.
        double well = k_b * exp_term;
        double dwell_dx = -2.0 * alpha * dr * well;
        bool new_form = false;
        if (rev_well_on && idx < static_cast<int>(m_rev_well.size()) && m_rev_well[idx].form != 0) {
            const RevWellPar& wp = m_rev_well[idx];
            double y, dy_dx;
            if (wp.form == 1) {
                // MG: y = exp(-(a x + beta x^2))
                const double phi = wp.p1 * dr + wp.p2 * dr * dr;
                const double yr = std::exp(-std::min(std::max(phi, -50.0), 200.0));
                y = yr;
                dy_dx = -(wp.p1 + 2.0 * wp.p2 * dr) * yr;
            } else {
                // erf-Morse: y = erfc((x - u)/sigma) / erfc(-u/sigma)
                const double n0 = std::erfc(-wp.p1 / wp.p2);
                const double z = (dr - wp.p1) / wp.p2;
                y = std::erfc(z) / n0;
                dy_dx = -2.0 / (1.7724538509055159 * wp.p2 * n0) * std::exp(-z * z);
            }
            // The inner side is capped at y = 2 (C1, exact below y = 1.6, so the fitted region
            // and the minimum itself are untouched): E = -D(2y - y^2) = D((y-1)^2 - 1) is then
            // bounded in [-D, 0] exactly as the Gaussian is in [k_b, 0]. The repulsive wall stays
            // the repulsion term's job, as it is for the Gaussian.
            constexpr double kYCap = 2.0, kYCapWidth = 0.2;
            const double t = y / kYCap;
            const double yc = kYCap * shareMinOne(t, kYCapWidth);
            const double dyc_dy = shareMinOneD(t, kYCapWidth);
            well = -wp.D * (2.0 * yc - yc * yc);
            dwell_dx = -wp.D * (2.0 - 2.0 * yc) * dyc_dy * dy_dx;
            new_form = true;
        }
        energy = new_form ? well : well * w;
        // rev-gfnff stage 3a(ii) (Claude Generated, Sep 2026): valence share.
        //     E = -k_b e^{-a dr^2} w c,     c = 1/2 (f_i + f_j),
        //     f_i = clip((Val_i - sum_{k != j} w_ik) / w_ij, 0, 1),
        // w_ij is the pair's OWN term weight, so f_i reads how much of atom i's valence the OTHER
        // partners already claim. A lone bond always has f > 1 (sum = 0) and a saturated
        // equilibrium atom has sum_{k != j} w = (n_i - 1) w < Val_i - w at w <= 1, so both are
        // clipped to exactly 1 and the equilibrium term stays bit-identical. It bites where two
        // partners share one valence (the exchange transition state), where the two wells then add
        // to one well's worth instead of two. Zero element data: the clip width is the unit
        // interval itself (shareClip) and Val_i is the over-coordination term's own valence table,
        // made hypervalent-correct by the settled-partner count (see prepareValenceShare).
        // What the share factor multiplies: the delivered Gaussian TIMES the term weight, or the
        // new well form on its own. Every dE/d(share input) below is this number times a purely
        // share-side derivative, so the two well forms share one code path from here on.
        const double base = energy;
        double cshare = 1.0, dcdw = 0.0;
        if (rev_share && w > 1e-12 && idx < static_cast<int>(m_rev_share_g.size())) {
            const double inv = 1.0 / w;
            // the "K" of the share formulas below is base/w: the share's own derivatives carry a
            // 1/w (f_i has w in its denominator), which is a property of the share, not the well.
            // For the delivered Gaussian that is the well itself - written that way rather than as
            // base*inv, because (K*w)*(1/w) is NOT bitwise K and the react MD amplifies one ulp
            // (Known Issue #33): measured, the round trip alone moved the 20-cell grid from 1186
            // to 1263 rebuilds.
            const double Kshare = new_form ? base * inv : well;
            // rev-gfnff stage 3a(ii) Sep 14, 2026: g = the genuineness of THIS pair (the smooth
            // 1,3 proxy, m_bonds order). It enters twice - the pair claims only w g of valence on
            // each end (so it is inside the sums, built by prepareValenceShare), and its own well
            // is multiplied by c = 1 - g (1 - (f_i + f_j)/2): g = 1 gives the plain share, g = 0
            // (a 1,3 contact) leaves the well exactly as it is.
            const double g = m_rev_share_g[idx];
            if (rev_share_conserving) {
                // rev-gfnff 3a(ii) "conserving" (Claude Generated, Sep 18, 2026, FABLE_REVIEW_2
                // A.5): f is a PER-ATOM factor, f_i = min(1, Val_i/S_i), and the pair takes the
                // PRODUCT. Both f and df/dS come from prepareConservingShare; the pair's own w
                // enters only through S_i, which the Lambda pass of applyValenceShareGradient
                // already carries, so there is no dc/dw term here and no settled-count channel.
                const double fi = m_rev_share_f(bond.i), fj = m_rev_share_f(bond.j);
                cshare = 1.0 - g * (1.0 - fi * fj);
                energy *= cshare;
                if (m_do_gradient) {
                    // dE/dS_i = (dE/df_i) (df_i/dS_i) = base g f_j * dfdS_i - stored directly as
                    // dE/d sum_i, which is what the Lambda pass multiplies by g_p dw/dr.
                    acc.dEdshare(bond.i) += base * g * fj * m_rev_share_dfdS(bond.i);
                    acc.dEdshare(bond.j) += base * g * fi * m_rev_share_dfdS(bond.j);
                    // dE/dg at fixed sums: c = 1 - g (1 - f_i f_j)
                    m_rev_share_dg[idx] += base * (-(1.0 - fi * fj));
                }
            } else {
            const double gi = (m_rev_share_val(bond.i) - m_rev_share_sum(bond.i) + w * g) * inv;
            const double gj = (m_rev_share_val(bond.j) - m_rev_share_sum(bond.j) + w * g) * inv;
            const double fi = shareClip(gi), fj = shareClip(gj);
            cshare = 1.0 - g * (1.0 - 0.5 * (fi + fj));
            energy *= cshare;
            if (m_do_gradient) {
                const double di = shareClipD(gi), dj = shareClipD(gj);
                // Claude Generated (Sep 13, 2026): d c / d w at fixed SUMS - the sum's own
                // w-dependence is supplied by the Lambda pass of applyValenceShareGradient.
                // With u_i = (Val_i - sum_i + w g)/w = (Val_i - sum_i)/w + g at constant sum_i,
                // du_i/dw = -(u_i - g)/w, and c = 1 - g (1 - (f_i + f_j)/2) carries the g in front:
                //     dc/dw = -(g/2) [ f_i'(u_i - g) + f_j'(u_j - g) ] / w.
                // The "sums free" form -g_i/w (the earlier code, whose comment already claimed
                // "fixed sums") double-counts Lambda: at the rkt06 exchange TS the analytic
                // gradient was then 0.4 % short of the FD of the share energy (1.4e-04 vs
                // 4.7e-06 Eh/Angstrom after the fix).
                dcdw = -0.5 * g * (di * (gi - g) + dj * (gj - g)) * inv;
                // dE_pair/d sum_i = -Kshare g f_i' / 2 (the pair's own S does not enter its own
                // sum, so no extra own-pair term here)
                acc.dEdshare(bond.i) += -0.5 * Kshare * g * di;
                acc.dEdshare(bond.j) += -0.5 * Kshare * g * dj;
                // dE_pair/d g = K w d c/d g, with f_i = clip((Val_i - sum_i)/w + g) so that
                // d f_i/d g = f_i' = shareClipD(g_i) = di:
                //     d c/d g = -(1 - (f_i + f_j)/2) + (g/2)(di + dj).
                // The g-dependence of the SUMS is the second half of the same derivative and is
                // applied in applyValenceShareGradient (it spans three atoms).
                m_rev_share_dg[idx] += base * (-(1.0 - 0.5 * (fi + fj)) + 0.5 * g * (di + dj));
            }
            }
        }
        acc.energy.bond += energy;
        // Claude Generated (Sep 14, 2026): CURCUMA_SHAREDUMP=1 also prints the per-pair WELL DEPTH
        // D_p = -k_b e^{-a dr^2} w (positive) next to the share factor it is multiplied by. The
        // prepareValenceShare() dump above cannot print it: the dynamic r0 (and hence the Gaussian)
        // is formed here, not there. D_p and the per-atom budget Val are the two inputs the
        // offline QP prototype of FABLE_BOND_STATE 3.2 needs. Env-gated, zero cost when unset;
        // with -threads 1 the order is the bond-list order.
        if (s_share_dump)
            CurcumaLogger::result(fmt::format(
                "shareD {:3d} {:3d}-{:3d} r {:9.5f} D {:14.8f} w {:9.6f} c {:9.6f} E {:14.8f}",
                idx, bond.i + 1, bond.j + 1, rij, -base, w, cshare, energy));

        if (m_do_gradient) {
            // d(w E_gauss)/dr = w dE_gauss/dr + E_gauss dw/dr, E_gauss = energy / w.
            // Only the first part depends on r0 and therefore on CN: dEdr stays the
            // Gaussian derivative for the chain rule below, the weight part goes to the
            // Cartesian gradient directly.
            // dE/dr at fixed w and c. Both well forms depend on r only through x = r - r0, so
            // this is also -dE/dr0 and the CN chain rule below is unchanged.
            // Same ulp caveat as Kshare above: for the Gaussian the delivered association
            // (-2 alpha dr) * energy is kept verbatim, which is bitwise what it always was.
            double dEdr = new_form ? (dwell_dx * cshare) : (-2.0 * alpha * dr * energy);
            // dE/dr through the term weight. Delivered: E = K w c, so dE/dw = K (c + w dc/dw).
            // New forms: E = well c and the weight is NOT on the well, so only the share's own
            // dc/dw survives, dE/dw = base dc/dw.
            const double dEdw = new_form ? (dcdw != 0.0 ? base * dcdw : 0.0)
                                         : well * (cshare + (dcdw != 0.0 ? w * dcdw : 0.0));
            double dEdr_cart = dEdr + ((w != 1.0 || dwdr != 0.0) ? dEdw * dwdr : 0.0);
            acc.gradient.row(bond.i) += dEdr_cart * derivate.row(0);
            acc.gradient.row(bond.j) += dEdr_cart * derivate.row(1);

            // HB alpha-modulation chain-rule gradient. Not applicable to the new well forms:
            // they do not use the modulated alpha (see prepareWellForms), so there is nothing to
            // differentiate - and applying the Gaussian's dE/d hb_cn_H to them would be wrong.
            if (bond.nr_hb >= 1 && !new_form) {
                constexpr double t1 = 0.1;
                double zz = t1 * alpha_orig * dr * dr * energy;
                int H = (m_atom_types[bond.i] == 1) ? bond.i : bond.j;
                if (H + 1 < static_cast<int>(m_hb_grad_offsets.size())) {
                    for (int k = m_hb_grad_offsets[H]; k < m_hb_grad_offsets[H + 1]; ++k) {
                        const auto& hbg = m_hb_grad_entries[m_hb_grad_list[k]];
                        acc.gradient.row(H) += zz * hbg.dCN_dH.transpose();
                        acc.gradient.row(hbg.B_atom) += zz * hbg.dCN_dB.transpose();
                    }
                }
            }

            // Bond dr0/dCN chain-rule contribution
            if (use_dynamic_r0 && bond.z_i > 0 && bond.z_j > 0) {
                double yy = -dEdr;
                acc.dEdcn(bond.i) += yy * bond.ff * bond.cnfak_i;
                acc.dEdcn(bond.j) += yy * bond.ff * bond.cnfak_j;
                acc.dEdcn_bond(bond.i) += yy * bond.ff * bond.cnfak_i;
                acc.dEdcn_bond(bond.j) += yy * bond.ff * bond.cnfak_j;

                // rev-gfnff stage 3a(i): with CN_i' = CN_i - cn_ij(r) + 1 the CN chain rule above
                // still supplies dE/dCN * dCN_i/dx (which contains the pair's own contribution),
                // so only its pair part has to be removed again: dCN_i'/dx = dCN_i/dx - d cn_ij/dx,
                // with d cn_ij/dx_i = +(dcn/dr) * d r/dx_i. The same dE/dCN coefficient is used,
                // hence dE/dx_i = -(dE/dCN_i + dE/dCN_j) * (dcn/dr) * dr/dx_i.
                if (rev_cn_pair && dcpair_dr != 0.0) {
                    const double coef = -(yy * bond.ff * (bond.cnfak_i + bond.cnfak_j)) * dcpair_dr;
                    acc.gradient.row(bond.i) += coef * derivate.row(0);
                    acc.gradient.row(bond.j) += coef * derivate.row(1);
                }
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_bond += (acc.gradient - grad_before);
}

// ============================================================================
// Angle Bending (GFN-FF cosine + distance damping)
// Reference: Fortran gfnff_engrad.F90:857-916
// ============================================================================

void FFWorkspace::calcAngles(FFAccumulator& acc, const std::vector<Angle>& list, std::pair<int, int> range)
{
    auto [begin, end] = range;
    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double pi = 3.14159265358979323846;
    const double linear_threshold = 1.0e-6;
    const double atcuta = tables().gen.atcuta; // runtime table (Sep 2026)
    constexpr double rcov_scale_angle = 4.0 / 3.0;

    auto get_rcov_bohr = [&](int atomic_number) -> double {
        if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_angle;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_angle;
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& angle = list[idx];
        auto i = m_geometry.row(angle.i);
        auto j = m_geometry.row(angle.j);
        auto k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_do_gradient);
        costheta = std::max(-1.0, std::min(1.0, costheta));

        double theta = std::acos(costheta);
        double theta0 = angle.theta0_ijk;
        double k_ijk = angle.fc;
        double energy, dedtheta;

        if (std::abs(pi - theta0) < linear_threshold) {
            double dtheta = theta - theta0;
            energy = k_ijk * dtheta * dtheta;
            dedtheta = 2.0 * k_ijk * dtheta;
        } else {
            double costheta0 = std::cos(theta0);
            double dcostheta = costheta - costheta0;
            energy = k_ijk * dcostheta * dcostheta;
            double sintheta = std::sin(theta);
            dedtheta = 2.0 * k_ijk * sintheta * (costheta0 - costheta);
        }

        double r_ij_sq = (i - j).squaredNorm();
        double r_jk_sq = (k - j).squaredNorm();

        double rcov_i = get_rcov_bohr(m_atom_types[angle.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[angle.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[angle.k]);

        double rcut_ij_sq = atcuta * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk_sq = atcuta * (rcov_j + rcov_k) * (rcov_j + rcov_k);

        double rr_ij = (r_ij_sq / rcut_ij_sq); rr_ij = rr_ij * rr_ij;
        double rr_jk = (r_jk_sq / rcut_jk_sq); rr_jk = rr_jk * rr_jk;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);

        double damp2ij = (r_ij_sq > 1e-8) ? -2.0 * 2.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
        double damp2jk = (r_jk_sq > 1e-8) ? -2.0 * 2.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;

        // rev-gfnff stage 1 (Sep 2026): each bond of the angle contributes its continuous
        // bond order as an extra factor; damp2 stays "(d factor/dr)/r" so the gradient code
        // below is unchanged.
        if (m_rev.enabled && m_rev.term_weights) {
            double r_ij = std::sqrt(r_ij_sq), r_jk = std::sqrt(r_jk_sq);
            double dw_ij = 0.0, dw_jk = 0.0;
            double w_ij = revWeight(angle.i, angle.j, r_ij, &dw_ij);
            double w_jk = revWeight(angle.j, angle.k, r_jk, &dw_jk);
            damp2ij = damp2ij * w_ij + damp_ij * dw_ij / std::max(r_ij, 1e-8);
            damp2jk = damp2jk * w_jk + damp_jk * dw_jk / std::max(r_jk, 1e-8);
            damp_ij *= w_ij;
            damp_jk *= w_jk;
        }
        double damp = damp_ij * damp_jk;

        acc.energy.angle += energy * damp;

        if (m_do_gradient) {
            // Claude Generated (Mar 2026): Use intermediate Vector to force column orientation
            // Mixing derivate.row() (1×3 row) with VectorXd (3×1 column) in += is undefined
            // in Eigen. Force all to column vectors first, matching ForceFieldThread pattern.
            Vector vab = (i - j).transpose();
            Vector vcb = (k - j).transpose();

            Vector term1 = energy * damp2ij * damp_jk * vab;
            Vector term2 = energy * damp2jk * damp_ij * vcb;

            Vector grad_i = dedtheta * damp * Vector(derivate.row(0));
            Vector grad_j = dedtheta * damp * Vector(derivate.row(1));
            Vector grad_k = dedtheta * damp * Vector(derivate.row(2));

            acc.gradient.row(angle.i) += (grad_i + term1).transpose();
            acc.gradient.row(angle.j) += (grad_j - term1 - term2).transpose();
            acc.gradient.row(angle.k) += (grad_k + term2).transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_angle += (acc.gradient - grad_before);
}

// ============================================================================
// Dihedral Torsion (GFN-FF with distance damping)
// Reference: Fortran gfnff_engrad.F90:1041-1122
// ============================================================================

void FFWorkspace::calcDihedrals(FFAccumulator& acc, const std::vector<Dihedral>& list, std::pair<int, int> range)
{
    auto [begin, end] = range;
    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double atcutt = tables().gen.atcutt;         // runtime tables (Sep 2026)
    const double atcutt_nci = tables().gen.atcutt_nci;
    constexpr double rcov_scale_tors = 4.0 / 3.0;

    auto get_rcov_bohr = [&](int atomic_number) -> double {
        if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_tors;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_tors;
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& dih = list[idx];

        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            m_geometry.row(dih.i).transpose(),
            m_geometry.row(dih.j).transpose(),
            m_geometry.row(dih.k).transpose(),
            m_geometry.row(dih.l).transpose(),
            derivate, m_do_gradient);

        double V = dih.V;
        double n = dih.n;
        double phi0 = dih.phi0;
        // Primary torsion: E = V * (1 + cos(n*(phi - phi0) + pi)) * damp
        // Reference: gfnff_engrad.F90:1268 — c1 = n*(phi-phi0) + pi
        double c1 = n * (phi - phi0) + M_PI;
        double energy = V * (1.0 + std::cos(c1));

        // Three-bond distance damping
        Eigen::Vector3d ri = m_geometry.row(dih.i).transpose();
        Eigen::Vector3d rj = m_geometry.row(dih.j).transpose();
        Eigen::Vector3d rk = m_geometry.row(dih.k).transpose();
        Eigen::Vector3d rl = m_geometry.row(dih.l).transpose();

        double r_ij_sq = (ri - rj).squaredNorm();
        double r_jk_sq = (rj - rk).squaredNorm();
        double r_kl_sq = (rk - rl).squaredNorm();

        double atcut = dih.is_nci ? atcutt_nci : atcutt;

        double rcov_i = get_rcov_bohr(m_atom_types[dih.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[dih.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[dih.k]);
        double rcov_l = get_rcov_bohr(m_atom_types[dih.l]);

        double rcut_ij = atcut * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk = atcut * (rcov_j + rcov_k) * (rcov_j + rcov_k);
        double rcut_kl = atcut * (rcov_k + rcov_l) * (rcov_k + rcov_l);

        double rr_ij = (r_ij_sq / rcut_ij); rr_ij *= rr_ij;
        double rr_jk = (r_jk_sq / rcut_jk); rr_jk *= rr_jk;
        double rr_kl = (r_kl_sq / rcut_kl); rr_kl *= rr_kl;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);
        double damp_kl = 1.0 / (1.0 + rr_kl);
        double damp2ij = (r_ij_sq > 1e-8) ? -4.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
        double damp2jk = (r_jk_sq > 1e-8) ? -4.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;
        double damp2kl = (r_kl_sq > 1e-8) ? -4.0 * rr_kl / (r_kl_sq * (1.0 + rr_kl) * (1.0 + rr_kl)) : 0.0;
        int rev_extra_i = -1, rev_extra_l = -1;
        double rev_d1 = 0.0, rev_d3 = 0.0, rev_r1 = 0.0, rev_r3 = 0.0, rev_w1 = 1.0, rev_w3 = 1.0, rev_damp0_ij = 0.0, rev_damp0_kl = 0.0;
        if (m_rev.enabled && m_rev.term_weights) { // rev-gfnff stage 1 (Sep 2026), see calcAngles
            double r1 = std::sqrt(r_ij_sq), r2 = std::sqrt(r_jk_sq), r3 = std::sqrt(r_kl_sq);
            double d1 = 0.0, d2 = 0.0, d3 = 0.0;
            // The stored quartet is NOT the bonded chain: the central bond is j-k, but the outer
            // atoms hang crosswise (i on k, l on j; the GFN-FF damping deliberately uses the
            // 1,3 distances i-j and k-l). The rev weight of an outer atom must be the weight of
            // its BOND, so it is taken from the bonded partner (Sep 12, 2026: with the wrong pairs
            // the weights were 0.5-0.96 at equilibrium and halved caffeine's torsion term).
            const int pi = m_bonded_pairs.count({ dih.i, dih.k }) ? dih.k : dih.j;
            const int pl = m_bonded_pairs.count({ dih.l, dih.j }) ? dih.j : dih.k;
            const double r1b = (m_geometry.row(dih.i) - m_geometry.row(pi)).norm();
            const double r3b = (m_geometry.row(dih.l) - m_geometry.row(pl)).norm();
            double w1 = revWeight(dih.i, pi, r1b, &d1), w2 = revWeight(dih.j, dih.k, r2, &d2), w3 = revWeight(dih.l, pl, r3b, &d3);
            (void)r1; (void)r3;
            // d w1 / d r acts on the bond i-pi, not on i-j: fold it into the Cartesian gradient of
            // that pair directly (below) and keep damp2 for the j-k chain-rule only
            damp2jk = damp2jk * w2 + damp_jk * d2 / std::max(r2, 1e-8);
            damp2ij *= w1; damp2kl *= w3;
            rev_extra_i = pi; rev_extra_l = pl; rev_d1 = d1; rev_d3 = d3; rev_r1 = r1b; rev_r3 = r3b; rev_w1 = w1; rev_w3 = w3;
            rev_damp0_ij = damp_ij; rev_damp0_kl = damp_kl;
            damp_ij *= w1; damp_jk *= w2; damp_kl *= w3;
        }
        double damp = damp_ij * damp_jk * damp_kl;

        acc.energy.dihedral += energy * damp;

        if (m_do_gradient) {
            double dEdphi = -V * n * std::sin(c1) * damp;

            // derivate.row() is 1×3: no .transpose() needed for row() +=
            acc.gradient.row(dih.i) += dEdphi * derivate.row(0);
            acc.gradient.row(dih.j) += dEdphi * derivate.row(1);
            acc.gradient.row(dih.k) += dEdphi * derivate.row(2);
            acc.gradient.row(dih.l) += dEdphi * derivate.row(3);

            // Damping gradient terms (damp2 = (d factor/dr)/r, bond weights included above)

            Eigen::Vector3d vij = ri - rj;
            Eigen::Vector3d vjk = rj - rk;
            Eigen::Vector3d vkl = rk - rl;

            Eigen::Vector3d t1 = energy * damp2ij * damp_jk * damp_kl * vij;
            Eigen::Vector3d t2 = energy * damp_ij * damp2jk * damp_kl * vjk;
            Eigen::Vector3d t3 = energy * damp_ij * damp_jk * damp2kl * vkl;

            // .transpose() converts Vector3d (3×1) to RowVector (1×3) for row() +=
            acc.gradient.row(dih.i) += t1.transpose();
            acc.gradient.row(dih.j) += (-t1 + t2).transpose();
            acc.gradient.row(dih.k) += (-t2 + t3).transpose();
            acc.gradient.row(dih.l) += (-t3).transpose();
            if (rev_extra_i >= 0) { // d w1 / d r on the bond i-pi and d w3 / d r on l-pl (rev-gfnff, Sep 12, 2026)
                if (rev_r1 > 1e-8) {
                    const double f1 = energy * rev_damp0_ij * damp_jk * damp_kl * rev_d1 / rev_r1;
                    Eigen::Vector3d g1 = f1 * (m_geometry.row(dih.i) - m_geometry.row(rev_extra_i)).transpose();
                    acc.gradient.row(dih.i) += g1.transpose();
                    acc.gradient.row(rev_extra_i) -= g1.transpose();
                }
                if (rev_r3 > 1e-8) {
                    const double f3 = energy * damp_ij * damp_jk * rev_damp0_kl * rev_d3 / rev_r3;
                    Eigen::Vector3d g3 = f3 * (m_geometry.row(dih.l) - m_geometry.row(rev_extra_l)).transpose();
                    acc.gradient.row(dih.l) += g3.transpose();
                    acc.gradient.row(rev_extra_l) -= g3.transpose();
                }
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Extra Torsions (sp3-sp3 gauche, same formula)
// ============================================================================

void FFWorkspace::calcExtraTorsions(FFAccumulator& acc, const std::vector<Dihedral>& list, std::pair<int, int> range)
{
    auto [begin, end] = range;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double atcutt = tables().gen.atcutt; // runtime table (Sep 2026)
    constexpr double rcov_scale_tors = 4.0 / 3.0;

    auto get_rcov_bohr = [&](int atomic_number) -> double {
        if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_tors;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_tors;
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& dih = list[idx];

        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            m_geometry.row(dih.i).transpose(),
            m_geometry.row(dih.j).transpose(),
            m_geometry.row(dih.k).transpose(),
            m_geometry.row(dih.l).transpose(),
            derivate, m_do_gradient);

        double V = dih.V;
        double n = dih.n;
        double phi0 = dih.phi0;
        // Extra torsion: same formula as primary — c1 = n*(phi-phi0) + π
        // Reference: gfnff_engrad.F90:1268 — ALL torsions with tlist(5)>0 use +π
        double c1 = n * (phi - phi0) + M_PI;
        double energy = V * (1.0 + std::cos(c1));

        Eigen::Vector3d ri = m_geometry.row(dih.i).transpose();
        Eigen::Vector3d rj = m_geometry.row(dih.j).transpose();
        Eigen::Vector3d rk = m_geometry.row(dih.k).transpose();
        Eigen::Vector3d rl = m_geometry.row(dih.l).transpose();

        double r_ij_sq = (ri - rj).squaredNorm();
        double r_jk_sq = (rj - rk).squaredNorm();
        double r_kl_sq = (rk - rl).squaredNorm();

        double rcov_i = get_rcov_bohr(m_atom_types[dih.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[dih.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[dih.k]);
        double rcov_l = get_rcov_bohr(m_atom_types[dih.l]);

        double rcut_ij = atcutt * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk = atcutt * (rcov_j + rcov_k) * (rcov_j + rcov_k);
        double rcut_kl = atcutt * (rcov_k + rcov_l) * (rcov_k + rcov_l);

        double rr_ij = (r_ij_sq / rcut_ij); rr_ij *= rr_ij;
        double rr_jk = (r_jk_sq / rcut_jk); rr_jk *= rr_jk;
        double rr_kl = (r_kl_sq / rcut_kl); rr_kl *= rr_kl;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);
        double damp_kl = 1.0 / (1.0 + rr_kl);
        double damp2ij = (r_ij_sq > 1e-8) ? -4.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
        double damp2jk = (r_jk_sq > 1e-8) ? -4.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;
        double damp2kl = (r_kl_sq > 1e-8) ? -4.0 * rr_kl / (r_kl_sq * (1.0 + rr_kl) * (1.0 + rr_kl)) : 0.0;
        int rev_extra_i = -1, rev_extra_l = -1;
        double rev_d1 = 0.0, rev_d3 = 0.0, rev_r1 = 0.0, rev_r3 = 0.0, rev_w1 = 1.0, rev_w3 = 1.0, rev_damp0_ij = 0.0, rev_damp0_kl = 0.0;
        if (m_rev.enabled && m_rev.term_weights) { // rev-gfnff stage 1 (Sep 2026), see calcAngles
            double r1 = std::sqrt(r_ij_sq), r2 = std::sqrt(r_jk_sq), r3 = std::sqrt(r_kl_sq);
            double d1 = 0.0, d2 = 0.0, d3 = 0.0;
            // The stored quartet is NOT the bonded chain: the central bond is j-k, but the outer
            // atoms hang crosswise (i on k, l on j; the GFN-FF damping deliberately uses the
            // 1,3 distances i-j and k-l). The rev weight of an outer atom must be the weight of
            // its BOND, so it is taken from the bonded partner (Sep 12, 2026: with the wrong pairs
            // the weights were 0.5-0.96 at equilibrium and halved caffeine's torsion term).
            const int pi = m_bonded_pairs.count({ dih.i, dih.k }) ? dih.k : dih.j;
            const int pl = m_bonded_pairs.count({ dih.l, dih.j }) ? dih.j : dih.k;
            const double r1b = (m_geometry.row(dih.i) - m_geometry.row(pi)).norm();
            const double r3b = (m_geometry.row(dih.l) - m_geometry.row(pl)).norm();
            double w1 = revWeight(dih.i, pi, r1b, &d1), w2 = revWeight(dih.j, dih.k, r2, &d2), w3 = revWeight(dih.l, pl, r3b, &d3);
            (void)r1; (void)r3;
            // d w1 / d r acts on the bond i-pi, not on i-j: fold it into the Cartesian gradient of
            // that pair directly (below) and keep damp2 for the j-k chain-rule only
            damp2jk = damp2jk * w2 + damp_jk * d2 / std::max(r2, 1e-8);
            damp2ij *= w1; damp2kl *= w3;
            rev_extra_i = pi; rev_extra_l = pl; rev_d1 = d1; rev_d3 = d3; rev_r1 = r1b; rev_r3 = r3b; rev_w1 = w1; rev_w3 = w3;
            rev_damp0_ij = damp_ij; rev_damp0_kl = damp_kl;
            damp_ij *= w1; damp_jk *= w2; damp_kl *= w3;
        }
        double damp = damp_ij * damp_jk * damp_kl;

        acc.energy.dihedral += energy * damp;

        if (m_do_gradient) {
            double dEdphi = -V * n * std::sin(c1) * damp;
            // derivate.row() is 1×3: no .transpose() needed for row() +=
            acc.gradient.row(dih.i) += dEdphi * derivate.row(0);
            acc.gradient.row(dih.j) += dEdphi * derivate.row(1);
            acc.gradient.row(dih.k) += dEdphi * derivate.row(2);
            acc.gradient.row(dih.l) += dEdphi * derivate.row(3);


            Eigen::Vector3d t1 = energy * damp2ij * damp_jk * damp_kl * (ri - rj);
            Eigen::Vector3d t2 = energy * damp_ij * damp2jk * damp_kl * (rj - rk);
            Eigen::Vector3d t3 = energy * damp_ij * damp_jk * damp2kl * (rk - rl);

            // .transpose() converts Vector3d (3×1) to RowVector (1×3) for row() +=
            acc.gradient.row(dih.i) += t1.transpose();
            acc.gradient.row(dih.j) += (-t1 + t2).transpose();
            acc.gradient.row(dih.k) += (-t2 + t3).transpose();
            acc.gradient.row(dih.l) += (-t3).transpose();
            if (rev_extra_i >= 0) { // d w1 / d r on the bond i-pi and d w3 / d r on l-pl (rev-gfnff, Sep 12, 2026)
                if (rev_r1 > 1e-8) {
                    const double f1 = energy * rev_damp0_ij * damp_jk * damp_kl * rev_d1 / rev_r1;
                    Eigen::Vector3d g1 = f1 * (m_geometry.row(dih.i) - m_geometry.row(rev_extra_i)).transpose();
                    acc.gradient.row(dih.i) += g1.transpose();
                    acc.gradient.row(rev_extra_i) -= g1.transpose();
                }
                if (rev_r3 > 1e-8) {
                    const double f3 = energy * damp_ij * damp_jk * rev_damp0_kl * rev_d3 / rev_r3;
                    Eigen::Vector3d g3 = f3 * (m_geometry.row(dih.l) - m_geometry.row(rev_extra_l)).transpose();
                    acc.gradient.row(dih.l) += g3.transpose();
                    acc.gradient.row(rev_extra_l) -= g3.transpose();
                }
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Inversions (out-of-plane bending)
// Reference: ForceFieldThread::CalculateGFNFFInversionContribution
// Uses domegadr analytical derivatives from gfnff_inversions.cpp
// ============================================================================

void FFWorkspace::calcInversions(FFAccumulator& acc, const std::vector<Inversion>& list, std::pair<int, int> range)
{
    // Claude Generated (Mar 2026): Complete inversion with gradient
    // Ported from ForceFieldThread::CalculateGFNFFInversionContribution
    // Reference: gfnff_engrad.F90:1355-1387
    auto [begin, end] = range;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    constexpr double rcov_scale = 4.0 / 3.0;

    auto get_rcov = [&](int at) -> double {
        if (at >= 1 && at <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[at - 1] * rcov_scale;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale;
    };

    const double atcutt_inv = tables().gen.atcutt; // runtime table (Sep 2026), captured by the lambdas
    auto calc_damp = [atcutt_inv](double rsq, double rcov_a, double rcov_b) -> double {
        double rcut = atcutt_inv * (rcov_a + rcov_b) * (rcov_a + rcov_b);
        double rr = (rsq / rcut) * (rsq / rcut);
        return 1.0 / (1.0 + rr);
    };

    auto calc_ddamp = [atcutt_inv](double r2_val, double rcov_a, double rcov_b) -> double {
        if (r2_val < 1e-8) return 0.0;
        double rcut_val = atcutt_inv * (rcov_a + rcov_b) * (rcov_a + rcov_b);
        double rr_val = (r2_val / rcut_val) * (r2_val / rcut_val);
        double one_plus_rr = 1.0 + rr_val;
        return -4.0 * rr_val / (r2_val * one_plus_rr * one_plus_rr);
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& inv = list[idx];

        // Atom layout: i=center, j=nb1, k=nb2, l=nb3
        Eigen::Vector3d r_center = m_geometry.row(inv.i).transpose();
        Eigen::Vector3d r_nb1 = m_geometry.row(inv.j).transpose();
        Eigen::Vector3d r_nb2 = m_geometry.row(inv.k).transpose();
        Eigen::Vector3d r_nb3 = m_geometry.row(inv.l).transpose();

        // Out-of-plane angle via calculateOutOfPlaneAngle
        Matrix derivate;
        double omega = GFNFF_Geometry::calculateOutOfPlaneAngle(
            r_center, r_nb1, r_nb2, r_nb3, derivate, m_do_gradient);

        // Damping: nb1 as hub (NOT star topology)
        // Reference: gfnff_engrad.F90:1356-1365
        double rcov_c = get_rcov(m_atom_types[inv.i]);
        double rcov_1 = get_rcov(m_atom_types[inv.j]);
        double rcov_2 = get_rcov(m_atom_types[inv.k]);
        double rcov_3 = get_rcov(m_atom_types[inv.l]);

        double rij_sq = (r_nb1 - r_center).squaredNorm();  // nb1-center: bond
        double rjk_sq = (r_nb1 - r_nb2).squaredNorm();     // nb1-nb2: 1-3 distance
        double rjl_sq = (r_nb1 - r_nb3).squaredNorm();     // nb1-nb3: 1-3 distance

        double damp_ij = calc_damp(rij_sq, rcov_c, rcov_1);
        double damp_jk = calc_damp(rjk_sq, rcov_2, rcov_1);
        double damp_jl = calc_damp(rjl_sq, rcov_1, rcov_3);
        // rev-gfnff stage 1 (Sep 2026): the three centre-neighbour bonds carry the weights
        double w_inv = 1.0, dw_c1 = 0.0, dw_c2 = 0.0, dw_c3 = 0.0, r_c1 = 0.0, r_c2 = 0.0, r_c3 = 0.0;
        if (m_rev.enabled && m_rev.term_weights) {
            r_c1 = std::sqrt(rij_sq);
            r_c2 = (r_nb2 - r_center).norm();
            r_c3 = (r_nb3 - r_center).norm();
            const double w1 = revWeight(inv.i, inv.j, r_c1, &dw_c1);
            const double w2 = revWeight(inv.i, inv.k, r_c2, &dw_c2);
            const double w3 = revWeight(inv.i, inv.l, r_c3, &dw_c3);
            w_inv = w1 * w2 * w3;
            // d(w_inv)/dr_ck expressed through the other two weights
            dw_c1 *= w2 * w3; dw_c2 *= w1 * w3; dw_c3 *= w1 * w2;
        }
        double damp = damp_ij * damp_jk * damp_jl * w_inv;

        // Energy and dE/domega
        double V = inv.fc;
        double et = 0.0;
        double dEdomega = 0.0;

        if (inv.potential_type == 0) {
            // Planar sp2: E = V*(1 - cos(omega)) * damp
            et = V * (1.0 - std::cos(omega));
            dEdomega = V * std::sin(omega) * damp;
        } else {
            // Saturated N: E = V*(cos(omega) - cos(omega0))^2 * damp
            double diff = std::cos(omega) - std::cos(inv.omega0);
            et = V * diff * diff;
            dEdomega = -2.0 * V * std::sin(omega) * diff * damp;
        }

        acc.energy.inversion += et * damp;

        if (m_do_gradient) {
            // PART 1: Omega derivative
            // Note: derivate.row() is 1×3; no transpose needed for row() += row()
            acc.gradient.row(inv.i) += dEdomega * derivate.row(0);
            acc.gradient.row(inv.j) += dEdomega * derivate.row(1);
            acc.gradient.row(inv.k) += dEdomega * derivate.row(2);
            acc.gradient.row(inv.l) += dEdomega * derivate.row(3);

            // PART 2: Damping gradient (nb1 as hub)
            // Reference: gfnff_engrad.F90:1379-1385
            double ddamp_ij = calc_ddamp(rij_sq, rcov_c, rcov_1);
            double ddamp_jk = calc_ddamp(rjk_sq, rcov_2, rcov_1);
            double ddamp_jl = calc_ddamp(rjl_sq, rcov_1, rcov_3);

            Eigen::Vector3d vab = r_nb1 - r_center;  // j - i
            Eigen::Vector3d vcb = r_nb1 - r_nb2;     // j - k
            Eigen::Vector3d vdc = r_nb1 - r_nb3;     // j - l

            Eigen::Vector3d term1 = (et * ddamp_ij * damp_jk * damp_jl) * vab;
            Eigen::Vector3d term2 = (et * ddamp_jk * damp_ij * damp_jl) * vcb;
            Eigen::Vector3d term3 = (et * ddamp_jl * damp_ij * damp_jk) * vdc;

            // center: -term1, nb1(hub): +term1+term2+term3, nb2: -term2, nb3: -term3
            // Note: .transpose() converts Vector3d (3×1) to RowVector (1×3) for row() +=
            // (the damping product below excludes w_inv, whose gradient follows separately)
            const double dpl = damp_ij * damp_jk * damp_jl;
            acc.gradient.row(inv.i) -= (w_inv * term1).transpose();
            acc.gradient.row(inv.j) += (w_inv * (term1 + term2 + term3)).transpose();
            acc.gradient.row(inv.k) -= (w_inv * term2).transpose();
            acc.gradient.row(inv.l) -= (w_inv * term3).transpose();
            if (m_rev.enabled && m_rev.term_weights) {
                // dE/dr_ck of the bond-order weights: et * dpl * d(w_inv)/dr_ck along the centre-neighbour vector
                auto add_pair = [&](int a, int b, const Eigen::Vector3d& vab, double r, double dw) {
                    if (r < 1e-8) return;
                    Eigen::Vector3d g = (et * dpl * dw / r) * vab;
                    acc.gradient.row(a) -= g.transpose();
                    acc.gradient.row(b) += g.transpose();
                };
                add_pair(inv.i, inv.j, r_nb1 - r_center, r_c1, dw_c1);
                add_pair(inv.i, inv.k, r_nb2 - r_center, r_c2, dw_c2);
                add_pair(inv.i, inv.l, r_nb3 - r_center, r_c3, dw_c3);
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Triple Bond Torsions (sTors_eg)
// Reference: Fortran gfnff_engrad.F90:3454
// ============================================================================

void FFWorkspace::calcSTorsions(FFAccumulator& acc, const std::vector<GFNFFSTorsion>& list, std::pair<int, int> range)
{
    auto [begin, end] = range;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& stor = list[idx];
        Matrix derivate;
        // GFN-FF uses Bohr coordinates internally (m_au = 1.0)
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            m_geometry.row(stor.i).transpose(),
            m_geometry.row(stor.j).transpose(),
            m_geometry.row(stor.k).transpose(),
            m_geometry.row(stor.l).transpose(),
            derivate, m_do_gradient);

        double erefhalf = stor.erefhalf;
        double energy = -erefhalf * std::cos(2.0 * phi) + erefhalf;
        acc.energy.stors += energy;

        if (m_do_gradient) {
            double dEdphi = 2.0 * erefhalf * std::sin(2.0 * phi);
            // derivate.row() is 1×3: no .transpose() needed for row() +=
            acc.gradient.row(stor.i) += dEdphi * derivate.row(0);
            acc.gradient.row(stor.j) += dEdphi * derivate.row(1);
            acc.gradient.row(stor.k) += dEdphi * derivate.row(2);
            acc.gradient.row(stor.l) += dEdphi * derivate.row(3);
        }
    }

    // Track sTorsion gradient in torsion component (matches ForceFieldThread)
    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Dispersion (GFN-FF modified BJ damping)
// Reference: gfnff_gdisp0.f90:365-377
// ============================================================================

void FFWorkspace::calcDispersion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].dispersions;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& disp = m_dispersions[idx];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();

        if (rij > disp.r_cut || rij < 1e-8) continue;

        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;
        double r0_6 = disp.r0_squared * disp.r0_squared * disp.r0_squared;
        double t6 = 1.0 / (r6 + r0_6);
        double r8 = r6 * r2;
        double r0_8 = r0_6 * disp.r0_squared;
        double t8 = 1.0 / (r8 + r0_8);

        double disp_sum = t6 + 2.0 * disp.r4r2ij * t8;
        double energy = -disp.C6 * disp_sum * disp.zetac6;
        acc.energy.dispersion += energy;

        if (m_do_gradient) {
            double d6 = -6.0 * r2 * r2 * t6 * t6;
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
            double ddisp_dr2 = d6 + 2.0 * disp.r4r2ij * d8;
            double dEdr = -disp.C6 * disp.zetac6 * ddisp_dr2 * rij;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(disp.i) += grad.transpose();
            acc.gradient.row(disp.j) -= grad.transpose();

            // dc6dcn chain-rule
            if (m_dc6dcn_ptr && m_dc6dcn_ptr->size() > 0 &&
                disp.i < m_dc6dcn_ptr->rows() && disp.j < m_dc6dcn_ptr->cols()) {
                double disp_value = disp_sum * disp.zetac6;
                acc.dEdcn(disp.i) -= (*m_dc6dcn_ptr)(disp.i, disp.j) * disp_value;
                acc.dEdcn(disp.j) -= (*m_dc6dcn_ptr)(disp.j, disp.i) * disp_value;
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_dispersion += (acc.gradient - grad_before);
}

// ============================================================================
// D4 Dispersion (same formula, separate parameter list)
// Reference: gfnff_gdisp0.f90:365-377
// ============================================================================

void FFWorkspace::calcD4Dispersion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].d4_dispersions;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& disp = m_d4_dispersions[idx];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();

        if (rij > disp.r_cut || rij < 1e-10) continue;

        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;
        double r0_6 = disp.r0_squared * disp.r0_squared * disp.r0_squared;
        double t6 = 1.0 / (r6 + r0_6);
        double r8 = r6 * r2;
        double r0_8 = r0_6 * disp.r0_squared;
        double t8 = 1.0 / (r8 + r0_8);

        double disp_sum = t6 + 2.0 * disp.r4r2ij * t8;
        double pair_energy = -disp.C6 * disp_sum * disp.zetac6;
        acc.energy.dispersion += pair_energy;

        if (m_do_gradient) {
            double d6 = -6.0 * r2 * r2 * t6 * t6;
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
            double ddisp_dr2 = d6 + 2.0 * disp.r4r2ij * d8;
            double dEdr = -disp.C6 * disp.zetac6 * ddisp_dr2 * rij;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(disp.i) += grad.transpose();
            acc.gradient.row(disp.j) -= grad.transpose();

            if (m_dc6dcn_ptr && m_dc6dcn_ptr->size() > 0 &&
                disp.i < m_dc6dcn_ptr->rows() && disp.j < m_dc6dcn_ptr->cols()) {
                double disp_value = disp_sum * disp.zetac6;
                acc.dEdcn(disp.i) -= (*m_dc6dcn_ptr)(disp.i, disp.j) * disp_value;
                acc.dEdcn(disp.j) -= (*m_dc6dcn_ptr)(disp.j, disp.i) * disp_value;
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_dispersion += (acc.gradient - grad_before);
}

// ============================================================================
// Bonded Repulsion — E = repab * exp(-α*r^1.5) / r
// Reference: Fortran gfnff_engrad.F90:467-495
// ============================================================================

void FFWorkspace::calcBondedRepulsion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].bonded_reps;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& rep = m_bonded_reps[idx];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();
        if (rij > rep.r_cut || rij < 1e-8) continue;

        double r_1_5 = rij * std::sqrt(rij);
        double base_energy, dEdr;
        if (rep.blend && m_rev.enabled && m_rev.blend_repulsion) {
            // rev-gfnff stage 1 (Sep 2026): b E_bonded + (1-b) E_nonbonded with the continuous
            // bond order b of the pair; the list membership no longer decides the parameter set.
            double dw = 0.0;
            const double w = revBlend(rep.i, rep.j, rij, &dw); // bonded list: crossover at 1.7x/-12, 1 - 1e-9 at 1.1x, 0.9993 at the turning point of a hot X-H bond (1.38x)
            const double eb = rep.repab_b * std::exp(-rep.alpha_b * r_1_5) / rij;
            const double en = rep.repab_n * std::exp(-rep.alpha_n * r_1_5) / rij;
            base_energy = w * eb + (1.0 - w) * en;
            const double deb = -eb / rij - 1.5 * rep.alpha_b * std::sqrt(rij) * eb;
            const double den = -en / rij - 1.5 * rep.alpha_n * std::sqrt(rij) * en;
            dEdr = w * deb + (1.0 - w) * den + (eb - en) * dw;
        } else {
            double exp_term = std::exp(-rep.alpha * r_1_5);
            base_energy = rep.repab * exp_term / rij;
            dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
        }
        acc.energy.bonded_rep += base_energy;

        if (m_do_gradient) {
            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(rep.i) += grad.transpose();
            acc.gradient.row(rep.j) -= grad.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_repulsion += (acc.gradient - grad_before);
}

// ============================================================================
// Non-bonded Repulsion — same formula, different parameter set
// Reference: Fortran gfnff_engrad.F90:255-276
// ============================================================================

void FFWorkspace::calcNonbondedRepulsion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].nonbonded_reps;
    if (begin == end) return;
    // Non-bonded repulsion goes to total gradient only, NOT to g_rep component

    for (int idx = begin; idx < end; ++idx) {
        const auto& rep = m_nonbonded_reps[idx];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();
        if (rij > rep.r_cut || rij < 1e-8) continue;

        double r_1_5 = rij * std::sqrt(rij);
        double base_energy, dEdr;
        if (rep.blend && m_rev.enabled && m_rev.blend_repulsion) {
            // rev-gfnff stage 1 (Sep 2026): b E_bonded + (1-b) E_nonbonded with the continuous
            // bond order b of the pair; the list membership no longer decides the parameter set.
            double dw = 0.0;
            const double w = revBlendNB(rep.i, rep.j, rij, &dw); // non-bonded list: tight switch (1.3x/-12), ~0 at 1,4 / H-bond distances
            const double eb = rep.repab_b * std::exp(-rep.alpha_b * r_1_5) / rij;
            const double en = rep.repab_n * std::exp(-rep.alpha_n * r_1_5) / rij;
            base_energy = w * eb + (1.0 - w) * en;
            const double deb = -eb / rij - 1.5 * rep.alpha_b * std::sqrt(rij) * eb;
            const double den = -en / rij - 1.5 * rep.alpha_n * std::sqrt(rij) * en;
            dEdr = w * deb + (1.0 - w) * den + (eb - en) * dw;
        } else {
            double exp_term = std::exp(-rep.alpha * r_1_5);
            base_energy = rep.repab * exp_term / rij;
            dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
        }
        acc.energy.nonbonded_rep += base_energy;

        if (m_do_gradient) {
            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(rep.i) += grad.transpose();
            acc.gradient.row(rep.j) -= grad.transpose();
        }
    }
}

// ============================================================================
// Coulomb (TERM 1: pairwise electrostatics)
// TERM 2+3 (self-energy) handled in postProcess()
// Reference: Fortran gfnff_engrad.F90:1378-1389
// ============================================================================

void FFWorkspace::calcCoulomb(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].coulombs;
    const auto [atom_begin, atom_end] = m_partitions[p].coulomb_atoms;
    if (begin == end && atom_begin == atom_end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    // Claude Generated (Sep 2026): implicit pairs. Everything a pair needs is per-atom - the EEQ
    // charge and alpeeq (gamma_ij = 1/sqrt(alp_i+alp_j)) - so the N^2/2 list does not have to
    // exist. Storing it costs 128 bytes per pair: 3.4 GB and ~0.5 s of pure write bandwidth at
    // 7320 atoms, which no amount of threading removes (measured). Same formula and the same
    // i<j pair set as the stored path below; the CUDA k_coulomb_implicit kernel is the device
    // twin of this loop.
    if (atom_begin < atom_end) {
        const double sqrt_pi = 1.772453850905516;
        const bool have_q = (m_eeq_charges.size() == m_natoms);
        const bool have_alp = (m_coul_alp.size() == m_natoms);
        if (have_q && have_alp) {
            for (int i = atom_begin; i < atom_end; ++i) {
                const double qi = m_eeq_charges(i);
                const double alp_i = m_coul_alp(i);
                if (std::isnan(qi) || alp_i <= 0.0) continue;
                const Eigen::Vector3d ri = m_geometry.row(i);
                for (int j = i + 1; j < m_natoms; ++j) {
                    const double qj = m_eeq_charges(j);
                    const double alp_j = m_coul_alp(j);
                    if (std::isnan(qj) || alp_j <= 0.0) continue;
                    const Eigen::Vector3d rij_vec = ri - m_geometry.row(j).transpose();
                    const double rij = rij_vec.norm();
                    if (rij > m_coulomb_implicit_rcut || rij < 1e-10) continue;

                    const double gamma_ij = 1.0 / std::sqrt(alp_i + alp_j);
                    const double gamma_r = gamma_ij * rij;
                    const double erf_term = curcuma_erf(gamma_r);
                    acc.energy.coulomb += qi * qj * erf_term / rij;

                    if (m_do_gradient) {
                        const double exp_term = std::exp(-gamma_r * gamma_r);
                        const double derf_dr = gamma_ij * exp_term * (2.0 / sqrt_pi);
                        const double dEdr_pair = qi * qj * (derf_dr / rij - erf_term / (rij * rij));
                        const Eigen::Vector3d grad = dEdr_pair * rij_vec / rij;
                        acc.gradient.row(i) += grad.transpose();
                        acc.gradient.row(j) -= grad.transpose();
                    }
                }
            }
        }
        if (acc.has_components && m_do_gradient)
            acc.grad_coulomb += (acc.gradient - grad_before);
        if (begin == end) return;
        if (acc.has_components && m_do_gradient) grad_before = acc.gradient;
    }

    for (int idx = begin; idx < end; ++idx) {
        const auto& coul = m_coulombs[idx];

        Eigen::Vector3d ri = m_geometry.row(coul.i);
        Eigen::Vector3d rj = m_geometry.row(coul.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();
        if (rij > coul.r_cut || rij < 1e-10) continue;

        // Dynamic EEQ charges (fall back to static if unavailable/NaN)
        double qi = coul.q_i, qj = coul.q_j;
        if (m_eeq_charges.size() > 0) {
            double qi_dyn = m_eeq_charges(coul.i);
            double qj_dyn = m_eeq_charges(coul.j);
            if (!std::isnan(qi_dyn) && !std::isnan(qj_dyn)) {
                qi = qi_dyn;
                qj = qj_dyn;
            }
        }

        double gamma_r = coul.gamma_ij * rij;
        double erf_term = curcuma_erf(gamma_r);
        double energy_pair = qi * qj * erf_term / rij;
        acc.energy.coulomb += energy_pair;

        if (m_do_gradient) {
            const double sqrt_pi = 1.772453850905516;
            double exp_term = std::exp(-gamma_r * gamma_r);
            double derf_dr = coul.gamma_ij * exp_term * (2.0 / sqrt_pi);
            double dEdr_pair = qi * qj * (derf_dr / rij - erf_term / (rij * rij));

            Eigen::Vector3d grad = dEdr_pair * rij_vec / rij;
            acc.gradient.row(coul.i) += grad.transpose();
            acc.gradient.row(coul.j) -= grad.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_coulomb += (acc.gradient - grad_before);
}

// ============================================================================
// HB/XB Damping Helper Functions (file-local)
// ============================================================================

namespace {

inline double ws_damping_out_of_line(double r_AH, double r_HB, double r_AB, double radab, double bacut)
{
    double ratio = (r_AH + r_HB) / r_AB;
    double exponent = (bacut / radab) * (ratio - 1.0);
    if (exponent > 15.0) return 0.0;
    return 2.0 / (1.0 + std::exp(exponent));
}

inline double ws_damping_short_range(double r, double r_vdw, double scut, double alp)
{
    double ratio = scut * r_vdw / (r * r);
    return 1.0 / (1.0 + std::pow(ratio, alp));
}

inline double ws_damping_long_range(double r, double longcut, double alp)
{
    return 1.0 / (1.0 + std::pow(r * r / longcut, alp));
}

inline double ws_charge_scaling(double q, double st, double sf)
{
    double exp_term = std::exp(st * q);
    return exp_term / (exp_term + sf);
}

} // anonymous namespace

// ============================================================================
// Hydrogen Bonds (three-body A-H...B)
// Reference: gfnff_engrad.F90 - abhgfnff_eg1, abhgfnff_eg2new, abhgfnff_eg3
// ============================================================================

void FFWorkspace::calcHydrogenBonds(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].hbonds;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& hb = m_hbonds[idx];

        Eigen::Vector3d pos_A = m_geometry.row(hb.i).transpose();
        Eigen::Vector3d pos_H = m_geometry.row(hb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(hb.k).transpose();

        double r_AH = (pos_H - pos_A).norm();
        double r_HB = (pos_B - pos_H).norm();
        double r_AB = (pos_B - pos_A).norm();

        if (r_AB > hb.r_cut) continue;

        double r_AH_4 = r_AH * r_AH * r_AH * r_AH;
        double r_HB_4 = r_HB * r_HB * r_HB * r_HB;
        double denom_DA = 1.0 / (r_AH_4 + r_HB_4);

        // HBond uses Phase 1 topology charges from struct (same as ForceFieldThread)
        // NOT Phase 2 EEQ charges — Fortran gfnff_engrad uses nhb1/nhb2 list charges
        double Q_H = ws_charge_scaling(hb.q_H, HB_ST, HB_SF);
        double Q_A = ws_charge_scaling(-hb.q_A, HB_ST, HB_SF);
        double Q_B = ws_charge_scaling(-hb.q_B, HB_ST, HB_SF);

        double bas = (Q_A * hb.basicity_A * r_AH_4 + Q_B * hb.basicity_B * r_HB_4) * denom_DA;
        double aci = (hb.acidity_B * r_AH_4 + hb.acidity_A * r_HB_4) * denom_DA;

        int elem_A = m_atom_types[hb.i];
        int elem_B = m_atom_types[hb.k];
        double r_vdw_AB = covalent_radii[elem_A - 1] + covalent_radii[elem_B - 1];

        double damp_short = ws_damping_short_range(r_AB, r_vdw_AB, HB_SCUT, HB_ALP);
        double damp_long = ws_damping_long_range(r_AB, HB_LONGCUT, HB_ALP);
        double damp_env = damp_short * damp_long;

        double damp_outl = ws_damping_out_of_line(r_AH, r_HB, r_AB, r_vdw_AB, HB_BACUT);

        double outl_nb_tot = 1.0;
        if (hb.case_type >= 2) {
            double hbnbcut_save = (elem_B == 7 && hb.neighbors_B.size() == 1) ? 2.0 : HB_NBCUT;
            for (int nb : hb.neighbors_B) {
                Eigen::Vector3d pos_nb = m_geometry.row(nb).transpose();
                double r_Anb = (pos_nb - pos_A).norm();
                double r_Bnb = (pos_nb - pos_B).norm();
                double expo_nb = (hbnbcut_save / r_vdw_AB) * ((r_Anb + r_Bnb) / r_AB - 1.0);
                outl_nb_tot *= (2.0 / (1.0 + std::exp(-expo_nb)) - 1.0);
            }
        }

        double rdamp;
        if (hb.case_type >= 2) {
            rdamp = damp_env * (1.8 / (r_HB * r_HB * r_HB) - 0.8 / (r_AB * r_AB * r_AB));
        } else {
            rdamp = damp_env / (r_AB * r_AB * r_AB);
        }

        // Case 4: Virtual LP (with gradient variables)
        // Claude Generated (Mar 2026): Port from ForceFieldThread with LP gradient support
        double outl_lp = 1.0;
        Eigen::Vector3d lp_pos = pos_B;
        double lp_dist = 0.0;
        int nbb_lp = 0;
        Eigen::Vector3d lp_vector = Eigen::Vector3d::Zero();
        if (hb.case_type == 4) {
            int z_B = m_atom_types[hb.k];
            double repz_B = (z_B >= 1 && z_B <= static_cast<int>(tables().repz.size()))
                          ? tables().repz[z_B - 1] : 1.0;
            lp_dist = 0.50 - 0.018 * repz_B;
            static constexpr double HBLPCUT = 56.0;

            nbb_lp = static_cast<int>(hb.neighbors_B.size());
            for (int nb_idx : hb.neighbors_B) {
                lp_vector += (m_geometry.row(nb_idx).transpose() - pos_B);
            }
            double vnorm = lp_vector.norm();
            if (vnorm > 1e-10 && nbb_lp > 0) {
                lp_pos = pos_B + (-lp_dist) * (lp_vector / vnorm);
                double ralp = (pos_A - lp_pos).norm();
                double expo_lp_val = (HBLPCUT / r_vdw_AB) * ((ralp + lp_dist + 1e-12) / r_AB - 1.0);
                outl_lp = 2.0 / (1.0 + std::exp(expo_lp_val));
            } else {
                nbb_lp = 0;
                outl_lp = 1.0;
            }
        }

        double qhoutl = Q_H * damp_outl * outl_nb_tot * outl_lp;

        // Case 3: eangl + etors (with gradient data)
        // Claude Generated (Mar 2026): Full port from ForceFieldThread with gradient
        // Reference: gfnff_engrad.F90:2574-2583 (egbend_nci_mul), 2529-2557 (egtors_nci_mul)
        double eangl = 1.0, etors = 1.0;
        struct TorsGrad {
            double energy;
            Eigen::Matrix<double, 4, 3> grad;
            int D_idx;
        };
        std::vector<TorsGrad> tors_data;
        Eigen::Matrix<double, 3, 3> gangl_3body = Eigen::Matrix<double, 3, 3>::Zero();
        int jj_idx = -1, kk_idx = -1, ll_idx = -1;

        if (hb.case_type == 3 && hb.acceptor_parent_index != -1) {
            int B_idx = hb.k;
            int C_idx = hb.acceptor_parent_index;
            int H_idx = hb.j;
            jj_idx = B_idx;
            kk_idx = C_idx;
            ll_idx = H_idx;

            // eangl: angle bending H...B=C (matching ForceFieldThread exactly)
            {
                double c0 = 120.0 * M_PI / 180.0;
                double fc_bend = 1.0 - BEND_HB;
                double kijk = fc_bend / ((std::cos(0.0) - std::cos(c0)) * (std::cos(0.0) - std::cos(c0)));

                Eigen::Vector3d va = m_geometry.row(C_idx).transpose();
                Eigen::Vector3d vb = m_geometry.row(B_idx).transpose();
                Eigen::Vector3d vc = m_geometry.row(H_idx).transpose();

                Eigen::Vector3d vab = va - vb;
                Eigen::Vector3d vcb = vc - vb;
                double rab2_a = vab.squaredNorm();
                double rcb2_a = vcb.squaredNorm();
                Eigen::Vector3d vp = vcb.cross(vab);
                double rp = vp.norm() + 1e-14;

                double cosa = vab.dot(vcb) / (std::sqrt(rab2_a) * std::sqrt(rcb2_a) + 1e-14);
                cosa = std::clamp(cosa, -1.0, 1.0);
                double theta_a = std::acos(cosa);

                double ea, deddt;
                if (M_PI - c0 < 1e-6) {
                    double dt = theta_a - c0;
                    ea = kijk * dt * dt;
                    deddt = 2.0 * kijk * dt;
                } else {
                    ea = kijk * (cosa - std::cos(c0)) * (cosa - std::cos(c0));
                    deddt = 2.0 * kijk * std::sin(theta_a) * (std::cos(c0) - cosa);
                }
                eangl = 1.0 - ea;

                if (m_do_gradient) {
                    Eigen::Vector3d deda_v = vab.cross(vp) * (-deddt / (rab2_a * rp));
                    Eigen::Vector3d dedc_v = vcb.cross(vp) * (deddt / (rcb2_a * rp));
                    Eigen::Vector3d dedb_v = deda_v + dedc_v;
                    gangl_3body.row(0) = dedb_v.transpose();
                    gangl_3body.row(1) = -deda_v.transpose();
                    gangl_3body.row(2) = -dedc_v.transpose();
                }
            }

            // etors: product of torsion terms D-B-C-H (with gradient)
            for (int D_idx : hb.neighbors_C) {
                if (D_idx == B_idx) continue;
                Matrix dihedral_grad;
                double phi = GFNFF_Geometry::calculateDihedralAngle(
                    m_geometry.row(D_idx).transpose(),
                    m_geometry.row(B_idx).transpose(),
                    m_geometry.row(C_idx).transpose(),
                    m_geometry.row(H_idx).transpose(),
                    dihedral_grad, m_do_gradient);

                double tshift = TORS_HB;
                double fc_tors = (1.0 - tshift) / 2.0;
                double phi0_tors = M_PI / 2.0;
                int rn = 2;
                double dphi1 = phi - phi0_tors;
                double c1 = rn * dphi1 + M_PI;
                double et = (1.0 + std::cos(c1)) * fc_tors + tshift;
                double dij = -rn * std::sin(c1) * fc_tors;

                TorsGrad tg;
                tg.energy = et;
                tg.D_idx = D_idx;
                tg.grad = Eigen::Matrix<double, 4, 3>::Zero();
                if (m_do_gradient && dihedral_grad.rows() == 4) {
                    for (int a = 0; a < 4; ++a) {
                        tg.grad.row(a) = dij * dihedral_grad.row(a);
                    }
                }
                tors_data.push_back(tg);
            }

            etors = 1.0;
            for (const auto& tg : tors_data) {
                etors *= tg.energy;
            }
        }

        // Energy: case-specific formula (matching ForceFieldThread)
        double global_scale = 1.0;
        if (hb.case_type == 2 || hb.case_type == 4) global_scale = XHACI_GLOBABH;
        else if (hb.case_type == 3) global_scale = XHACI_COH;

        double E_HB;
        if (hb.case_type >= 2) {
            double const_val = hb.acidity_A * hb.basicity_B * Q_A * Q_B * global_scale;
            E_HB = -rdamp * qhoutl * const_val * eangl * etors;
        } else {
            E_HB = -bas * aci * rdamp * qhoutl;
        }
        acc.energy.hbond += E_HB;

        static const bool hb_dump = (std::getenv("CURCUMA_HB_DUMP") != nullptr);  // once, not per triple
        if (hb_dump && std::abs(E_HB) > 1e-13) {
            fmt::print("HBTRIP A={:3d} H={:3d} B={:3d} case={} E={:.12f}\n",
                       hb.i+1, hb.j+1, hb.k+1, hb.case_type, E_HB);
        }

        // Claude Generated (May 2026, HB-investigation): per-case split for Fortran comparison
        switch (hb.case_type) {
            case 1: acc.energy.hbond_case1 += E_HB; ++acc.energy.hbond_case1_count; break;
            case 2: acc.energy.hbond_case2 += E_HB; ++acc.energy.hbond_case2_count; break;
            case 3: acc.energy.hbond_case3 += E_HB; ++acc.energy.hbond_case3_count; break;
            case 4: acc.energy.hbond_case4 += E_HB; ++acc.energy.hbond_case4_count; break;
            default: break;
        }

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Claude Generated (Mar 2026): Complete port from ForceFieldThread
        // Reference: gfnff_engrad.F90 - abhgfnff_eg1/eg2new/eg2_rnr/eg3
        if (m_do_gradient) {
            // Fortran convention distance vectors (all in Bohr)
            Eigen::Vector3d drab = pos_A - pos_B;  // A - B
            Eigen::Vector3d drah = pos_A - pos_H;  // A - H
            Eigen::Vector3d drbh = pos_B - pos_H;  // B - H

            double rab2 = r_AB * r_AB;
            double rbh2 = r_HB * r_HB;
            double rah2 = r_AH * r_AH;
            double rahprbh = r_AH + r_HB + 1e-12;

            // Damping derivative intermediates
            double ratio1 = std::pow(rab2 / HB_LONGCUT, HB_ALP);
            double shortcut = HB_SCUT * r_vdw_AB;
            double ratio3 = std::pow(shortcut / rab2, HB_ALP);
            double ddamp = (-2.0 * HB_ALP * ratio1 / (1.0 + ratio1))
                         + ( 2.0 * HB_ALP * ratio3 / (1.0 + ratio3));

            // Out-of-line intermediates for gradient
            double expo = (HB_BACUT / r_vdw_AB) * (rahprbh / r_AB - 1.0);
            if (expo > 15.0) continue;  // Fortran early return (gfnff_engrad.F90:1781)
            double ratio2 = std::exp(expo);

            Eigen::Vector3d ga = Eigen::Vector3d::Zero();
            Eigen::Vector3d gb = Eigen::Vector3d::Zero();
            Eigen::Vector3d gh = Eigen::Vector3d::Zero();

            if (hb.case_type >= 2) {
                // ===== Case 2/3/4: abhgfnff_eg2new gradient =====
                const double p_bh = 1.0 + tables().gen.hbabmix;   // 1 + hbabmix (runtime table)
                const double p_ab = -tables().gen.hbabmix;  // -hbabmix
                double rbhdamp = damp_env * p_bh / (rbh2 * r_HB);
                double rabdamp = damp_env * p_ab / (rab2 * r_AB);

                double const_val = hb.acidity_A * hb.basicity_B * Q_A * Q_B * global_scale;
                double dterm  = -qhoutl * eangl * etors * const_val;
                double aterm  = -rdamp * Q_H * outl_nb_tot * outl_lp * eangl * etors * const_val;
                double nbterm = -rdamp * Q_H * damp_outl * outl_lp * eangl * etors * const_val;

                // Damping part: rab
                double gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2;
                gi *= dterm;
                Eigen::Vector3d dg = gi * drab;
                ga = dg;
                gb = -dg;

                // Damping part: rbh
                gi = -3.0 * rbhdamp / rbh2;
                gi *= dterm;
                dg = gi * drbh;
                gb += dg;
                gh = -dg;

                // Out-of-line: rab
                double tmp1 = -2.0 * aterm * ratio2 * expo
                            / ((1.0 + ratio2) * (1.0 + ratio2))
                            / (rahprbh - r_AB);
                gi = -tmp1 * rahprbh / rab2;
                dg = gi * drab;
                ga += dg;
                gb -= dg;

                // Out-of-line: rah, rbh
                gi = tmp1 / r_AH;
                Eigen::Vector3d dga_outl = gi * drah;
                ga += dga_outl;
                gi = tmp1 / r_HB;
                Eigen::Vector3d dgb_outl = gi * drbh;
                gb += dgb_outl;
                gh += -dga_outl - dgb_outl;

                // Neighbor out-of-line gradient (Case >= 2)
                double hbnbcut_g = (elem_B == 7 && hb.neighbors_B.size() == 1) ? 2.0 : HB_NBCUT;
                for (size_t nb_i = 0; nb_i < hb.neighbors_B.size(); ++nb_i) {
                    int nb = hb.neighbors_B[nb_i];
                    Eigen::Vector3d pos_nb = m_geometry.row(nb).transpose();
                    Eigen::Vector3d dranb = pos_A - pos_nb;
                    Eigen::Vector3d drbnb = pos_B - pos_nb;
                    double ranb = dranb.norm();
                    double rbnb = drbnb.norm();
                    double ranbprbnb = ranb + rbnb + 1e-12;

                    double expo_nb_i = (hbnbcut_g / r_vdw_AB) * (ranbprbnb / r_AB - 1.0);
                    double ratio2_nb_i = std::exp(-expo_nb_i);
                    double outl_nb_i = 2.0 / (1.0 + ratio2_nb_i) - 1.0;

                    double outl_nb_others = 1.0;
                    if (std::abs(outl_nb_i) > 1e-12) {
                        outl_nb_others = outl_nb_tot / outl_nb_i;
                    }

                    double tmp2 = 2.0 * nbterm * outl_nb_others * ratio2_nb_i * expo_nb_i
                                / ((1.0 + ratio2_nb_i) * (1.0 + ratio2_nb_i))
                                / (ranbprbnb - r_AB);

                    double gi_nb = -tmp2 * ranbprbnb / rab2;
                    dg = gi_nb * drab;
                    ga += dg;
                    gb -= dg;

                    gi_nb = tmp2 / ranb;
                    Eigen::Vector3d dga_nb = gi_nb * dranb;
                    ga += dga_nb;
                    gi_nb = tmp2 / rbnb;
                    Eigen::Vector3d dgb_nb = gi_nb * drbnb;
                    gb += dgb_nb;
                    acc.gradient.row(nb) += (-dga_nb - dgb_nb).transpose();
                }

                // Case 4: LP out-of-line gradient
                if (hb.case_type == 4 && nbb_lp > 0) {
                    double lpterm = -rdamp * Q_H * damp_outl * outl_nb_tot * const_val;
                    double ralp = (pos_A - lp_pos).norm();
                    double rblp = lp_dist;
                    double ralpprblp = ralp + rblp + 1e-12;
                    static constexpr double HBLPCUT = 56.0;
                    double expo_lp_val = (HBLPCUT / r_vdw_AB) * (ralpprblp / r_AB - 1.0);
                    double ratio2_lp_val = std::exp(expo_lp_val);

                    // LP out-of-line: rab
                    double tmp3 = -2.0 * lpterm * ratio2_lp_val * expo_lp_val
                                / ((1.0 + ratio2_lp_val) * (1.0 + ratio2_lp_val))
                                / (ralpprblp - r_AB);
                    double gi_lp = -tmp3 * ralpprblp / rab2;
                    Eigen::Vector3d dg_lp = gi_lp * drab;
                    ga += dg_lp;
                    gb -= dg_lp;

                    // LP out-of-line: ralp
                    Eigen::Vector3d dralp = pos_A - lp_pos;
                    gi_lp = tmp3 / (ralp + 1e-12);
                    Eigen::Vector3d dga_lp = gi_lp * dralp;
                    ga += dga_lp;

                    // Fortran: gb -= dga (uses dga, not dgb)
                    gb -= dga_lp;
                    Eigen::Vector3d glp = -dga_lp;

                    // LP neighbor chain rule
                    double vnorm = lp_vector.norm();
                    if (vnorm > 1e-10) {
                        Eigen::Matrix3d gii = Eigen::Matrix3d::Zero();
                        for (int col = 0; col < 3; ++col) {
                            Eigen::Vector3d unit_vec = Eigen::Vector3d::Zero();
                            unit_vec(col) = -1.0;
                            gii.col(col) = -lp_dist * static_cast<double>(nbb_lp)
                                         * (unit_vec / vnorm + lp_vector * lp_vector(col) / std::pow(vnorm, 3.0));
                        }
                        Eigen::Vector3d gnb_lp = gii * glp;
                        gb += gnb_lp;
                        Eigen::Vector3d gnb_lp_share = gnb_lp / static_cast<double>(nbb_lp);
                        for (int nb_idx : hb.neighbors_B) {
                            acc.gradient.row(nb_idx) -= gnb_lp_share.transpose();
                        }
                    }
                }

                // Case 3: angle bending and torsion gradient contributions
                if (hb.case_type == 3 && jj_idx >= 0) {
                    double bterm_c3 = -rdamp * qhoutl * etors * const_val;
                    double tterm_c3 = -rdamp * qhoutl * eangl * const_val;

                    acc.gradient.row(jj_idx) += bterm_c3 * gangl_3body.row(0);
                    acc.gradient.row(kk_idx) += bterm_c3 * gangl_3body.row(1);
                    acc.gradient.row(ll_idx) += bterm_c3 * gangl_3body.row(2);

                    for (size_t k = 0; k < tors_data.size(); ++k) {
                        double factor_k = (std::abs(tors_data[k].energy) > 1e-12)
                                        ? etors / tors_data[k].energy : 0.0;
                        double t_k = factor_k * tterm_c3;
                        acc.gradient.row(tors_data[k].D_idx) += t_k * tors_data[k].grad.row(0);
                        acc.gradient.row(jj_idx)             += t_k * tors_data[k].grad.row(1);
                        acc.gradient.row(kk_idx)             += t_k * tors_data[k].grad.row(2);
                        acc.gradient.row(ll_idx)             += t_k * tors_data[k].grad.row(3);
                    }
                }
            } else {
                // ===== Case 1: abhgfnff_eg1 gradient =====
                double caa = Q_A * hb.basicity_A;
                double cbb = Q_B * hb.basicity_B;

                double rterm = -aci * rdamp * qhoutl;
                double dterm = -aci * bas * qhoutl;
                double sterm = -rdamp * bas * qhoutl;
                double aterm = -aci * bas * rdamp * Q_H;

                double denom_val = 1.0 / (r_AH_4 + r_HB_4);
                double tmp = denom_val * denom_val * 4.0;
                double dd24a = rah2 * r_HB_4 * tmp;
                double dd24b = rbh2 * r_AH_4 * tmp;

                // Donor-acceptor: bas
                double gi = (caa - cbb) * dd24a * rterm;
                ga = gi * drah;
                gi = (cbb - caa) * dd24b * rterm;
                gb = gi * drbh;
                gh = -ga - gb;

                // Donor-acceptor: aci
                gi = (hb.acidity_B - hb.acidity_A) * dd24a;
                Eigen::Vector3d dga_aci = gi * drah * sterm;
                ga += dga_aci;
                gi = (hb.acidity_A - hb.acidity_B) * dd24b;
                Eigen::Vector3d dgb_aci = gi * drbh * sterm;
                gb += dgb_aci;
                gh += -dga_aci - dgb_aci;

                // Damping: rab
                gi = rdamp * (ddamp - 3.0) / rab2;
                Eigen::Vector3d dg = gi * drab * dterm;
                ga += dg;
                gb -= dg;

                // Out-of-line: rab
                gi = aterm * 2.0 * ratio2 * expo * rahprbh
                   / ((1.0 + ratio2) * (1.0 + ratio2))
                   / (rahprbh - r_AB) / rab2;
                dg = gi * drab;
                ga += dg;
                gb -= dg;

                // Out-of-line: rah, rbh
                double tmp_outl = -2.0 * aterm * ratio2 * expo
                                / ((1.0 + ratio2) * (1.0 + ratio2))
                                / (rahprbh - r_AB);
                Eigen::Vector3d dga_outl = drah * tmp_outl / r_AH;
                ga += dga_outl;
                Eigen::Vector3d dgb_outl = drbh * tmp_outl / r_HB;
                gb += dgb_outl;
                gh += -dga_outl - dgb_outl;
            }

            // Accumulate: A=hb.i, B=hb.k, H=hb.j
            acc.gradient.row(hb.i) += ga.transpose();
            acc.gradient.row(hb.k) += gb.transpose();
            acc.gradient.row(hb.j) += gh.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_hb += (acc.gradient - grad_before);
}

// ============================================================================
// Halogen Bonds (three-body A-X...B)
// Reference: gfnff_engrad.F90 - rbxgfnff_eg
// ============================================================================

void FFWorkspace::calcHalogenBonds(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].xbonds;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& xb = m_xbonds[idx];

        Eigen::Vector3d pos_A = m_geometry.row(xb.i).transpose();
        Eigen::Vector3d pos_X = m_geometry.row(xb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(xb.k).transpose();

        Eigen::Vector3d r_AX_vec = pos_X - pos_A;
        Eigen::Vector3d r_XB_vec = pos_B - pos_X;
        Eigen::Vector3d r_AB_vec = pos_B - pos_A;

        double r_AX = r_AX_vec.norm();
        double r_XB = r_XB_vec.norm();
        double r_AB = r_AB_vec.norm();

        if (r_XB > xb.r_cut) continue;

        // Claude Generated (Mar 2026): Fixed to match ForceFieldThread exactly
        int elem_A = m_atom_types[xb.i];
        int elem_B = m_atom_types[xb.k];
        double r_vdw_AB = covalent_radii[elem_A - 1] + covalent_radii[elem_B - 1];

        double damp_short = ws_damping_short_range(r_XB, r_vdw_AB, XB_SCUT, HB_ALP);
        double damp_long = ws_damping_long_range(r_XB, HB_LONGCUT_XB, HB_ALP);

        // XB out-of-line: uses xbacut directly, not divided by radab
        double ratio_outl = (r_AX + r_XB) / r_AB;
        double expo_outl = XB_BACUT * (ratio_outl - 1.0);
        if (expo_outl > 15.0) continue;  // Early exit matching Fortran
        double damp_outl = 2.0 / (1.0 + std::exp(expo_outl));

        double Q_X = ws_charge_scaling(xb.q_X, XB_ST, XB_SF);
        double Q_B = ws_charge_scaling(-xb.q_B, XB_ST, XB_SF);

        // Energy: R_damp = damp_short * damp_long * damp_outl / r_XB^3
        // (Fixed: was r_AB^3, now r_XB^3 matching ForceFieldThread/Fortran)
        double R_damp = damp_short * damp_long * damp_outl / (r_XB * r_XB * r_XB);
        double E_XB = -R_damp * Q_B * xb.acidity_X * Q_X;
        acc.energy.xbond += E_XB;

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Claude Generated (Mar 2026): Port from ForceFieldThread
        // Reference: gfnff_engrad.F90 - rbxgfnff_eg()
        if (m_do_gradient) {
            // Short-range damping derivative w.r.t. r_XB
            double ratio_short = XB_SCUT * r_vdw_AB / (r_XB * r_XB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term
                                  / (r_XB * (1.0 + damp_short_term));

            // Long-range damping derivative w.r.t. r_XB
            double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_XB * damp_long * damp_long_term
                                 / (HB_LONGCUT_XB * (1.0 + damp_long_term));

            // Out-of-line damping derivatives
            double exp_term = std::exp(expo_outl);
            double denom_outl = 1.0 + exp_term;
            double ddamp_outl_drAX = -2.0 * exp_term * XB_BACUT
                                   / (r_AB * denom_outl * denom_outl);
            double ddamp_outl_drXB = ddamp_outl_drAX;
            double ddamp_outl_drAB = 2.0 * exp_term * XB_BACUT * (r_AX + r_XB)
                                   / (r_AB * r_AB * denom_outl * denom_outl);

            // R_damp chain rule derivatives
            double dRdamp_drXB = (ddamp_short_dr * damp_long * damp_outl
                                + damp_short * ddamp_long_dr * damp_outl
                                + damp_short * damp_long * ddamp_outl_drXB) / (r_XB * r_XB * r_XB)
                               - 3.0 * R_damp / r_XB;
            double dRdamp_drAX = damp_short * damp_long * ddamp_outl_drAX / (r_XB * r_XB * r_XB);
            double dRdamp_drAB = damp_short * damp_long * ddamp_outl_drAB / (r_XB * r_XB * r_XB);

            double E_prefactor = -Q_B * xb.acidity_X * Q_X;
            double dE_drXB = E_prefactor * dRdamp_drXB;
            double dE_drAX = E_prefactor * dRdamp_drAX;
            double dE_drAB = E_prefactor * dRdamp_drAB;

            Eigen::Vector3d grad_rAX_unit = r_AX_vec / (r_AX + 1e-14);
            Eigen::Vector3d grad_rXB_unit = r_XB_vec / (r_XB + 1e-14);
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / (r_AB + 1e-14);

            // A: dr_AX/dA = -unit_AX, dr_AB/dA = -unit_AB
            Eigen::Vector3d grad_A = -dE_drAX * grad_rAX_unit - dE_drAB * grad_rAB_unit;
            // X: dr_AX/dX = +unit_AX, dr_XB/dX = -unit_XB
            Eigen::Vector3d grad_X = dE_drAX * grad_rAX_unit - dE_drXB * grad_rXB_unit;
            // B: dr_XB/dB = +unit_XB, dr_AB/dB = +unit_AB
            Eigen::Vector3d grad_B = dE_drXB * grad_rXB_unit + dE_drAB * grad_rAB_unit;

            acc.gradient.row(xb.i) += grad_A.transpose();
            acc.gradient.row(xb.j) += grad_X.transpose();
            acc.gradient.row(xb.k) += grad_B.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_xb += (acc.gradient - grad_before);
}

// ============================================================================
// ATM Three-Body Dispersion (Axilrod-Teller-Muto)
// Reference: external/cpp-d4/src/damping/atm.cpp:70-138
// ============================================================================

void FFWorkspace::calcATM(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].atm_triples;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& triple = m_atm_triples[idx];

        Eigen::Vector3d pos_i = m_geometry.row(triple.i).transpose();
        Eigen::Vector3d pos_j = m_geometry.row(triple.j).transpose();
        Eigen::Vector3d pos_k = m_geometry.row(triple.k).transpose();

        double rij = (pos_i - pos_j).norm();
        double rik = (pos_i - pos_k).norm();
        double rjk = (pos_j - pos_k).norm();

        double r2ij = rij * rij, r2ik = rik * rik, r2jk = rjk * rjk;

        double c9 = triple.s9 * std::sqrt(std::fabs(triple.C6_ij * triple.C6_ik * triple.C6_jk));

        int zi = m_atom_types[triple.i], zj = m_atom_types[triple.j], zk = m_atom_types[triple.k];
        // Claude Generated (May 2026, GPU/CPU 8.9 µEh fix): GFN-FF covalent radii
        // (covalent_rad_d3). Earlier rcov_bohr (= r0_gfnff) was a GFN-FF-specific bond-r0
        // table, not the covalent radii ATM expects. Drove the polymer ATM mismatch.
        // (Jul 2026) The GPU uploads its rcov from this same covalent_rad_d3 array now.
        double r_cov_i = (zi > 0 && zi <= static_cast<int>(covalent_rad_d3.size())) ? covalent_rad_d3[zi - 1] : 1.0;
        double r_cov_j = (zj > 0 && zj <= static_cast<int>(covalent_rad_d3.size())) ? covalent_rad_d3[zj - 1] : 1.0;
        double r_cov_k = (zk > 0 && zk <= static_cast<int>(covalent_rad_d3.size())) ? covalent_rad_d3[zk - 1] : 1.0;

        double r0ij = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_j) + triple.a2;
        double r0ik = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_k) + triple.a2;
        double r0jk = triple.a1 * std::sqrt(3.0 * r_cov_j * r_cov_k) + triple.a2;

        double rijk = rij * rik * rjk;
        double r2ijk = r2ij * r2ik * r2jk;
        double r3ijk = rijk * r2ijk;

        double fdmp = 1.0 / (1.0 + 6.0 * std::pow(r0ij * r0ik * r0jk / rijk, triple.alp / 3.0));

        double A = r2ij + r2jk - r2ik;
        double B = r2ij + r2ik - r2jk;
        double C = r2ik + r2jk - r2ij;
        double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

        acc.energy.atm += ang * fdmp * c9 / 3.0 * triple.triple_scale;
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_atm += (acc.gradient - grad_before);
}

// ============================================================================
// ATM Gradient (analytical)
// Reference: external/cpp-d4/src/damping/atm.cpp:141-289
// ============================================================================

void FFWorkspace::calcATMGradient(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].atm_triples;
    if (begin == end) return;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& triple = m_atm_triples[idx];

        Eigen::Vector3d pos_i = m_geometry.row(triple.i).transpose();
        Eigen::Vector3d pos_j = m_geometry.row(triple.j).transpose();
        Eigen::Vector3d pos_k = m_geometry.row(triple.k).transpose();

        Eigen::Vector3d rij_vec = pos_j - pos_i;
        Eigen::Vector3d rik_vec = pos_k - pos_i;
        Eigen::Vector3d rjk_vec = pos_k - pos_j;

        double rij = rij_vec.norm(), rik = rik_vec.norm(), rjk = rjk_vec.norm();
        double r2ij = rij * rij, r2ik = rik * rik, r2jk = rjk * rjk;

        // Negative c9 for gradient (dispersion is attractive)
        double c9 = -triple.s9 * std::sqrt(std::fabs(triple.C6_ij * triple.C6_ik * triple.C6_jk));

        int zi = m_atom_types[triple.i], zj = m_atom_types[triple.j], zk = m_atom_types[triple.k];
        // Claude Generated (May 2026): D3 covalent radii (matches GPU + ATM theory).
        double r_cov_i = (zi > 0 && zi <= static_cast<int>(covalent_rad_d3.size())) ? covalent_rad_d3[zi - 1] : 1.0;
        double r_cov_j = (zj > 0 && zj <= static_cast<int>(covalent_rad_d3.size())) ? covalent_rad_d3[zj - 1] : 1.0;
        double r_cov_k = (zk > 0 && zk <= static_cast<int>(covalent_rad_d3.size())) ? covalent_rad_d3[zk - 1] : 1.0;

        double r0ij = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_j) + triple.a2;
        double r0ik = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_k) + triple.a2;
        double r0jk = triple.a1 * std::sqrt(3.0 * r_cov_j * r_cov_k) + triple.a2;
        double r0ijk = r0ij * r0ik * r0jk;

        double rijk = rij * rik * rjk;
        double r2ijk = r2ij * r2ik * r2jk;
        double r3ijk = rijk * r2ijk;
        double r5ijk = r2ijk * r3ijk;

        double tmp = std::pow(r0ijk / rijk, triple.alp / 3.0);
        double fdmp = 1.0 / (1.0 + 6.0 * tmp);
        double dfdmp = -2.0 * triple.alp * tmp * fdmp * fdmp;

        double A = r2ij + r2jk - r2ik;
        double B = r2ij + r2ik - r2jk;
        double C = r2ik + r2jk - r2ij;
        double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

        double dang_ij = -0.375 * (std::pow(r2ij, 3) + std::pow(r2ij, 2) * (r2jk + r2ik)
                          + r2ij * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * std::pow(r2ik, 2))
                          - 5.0 * std::pow(r2jk - r2ik, 2) * (r2jk + r2ik)) / r5ijk;

        double dang_ik = -0.375 * (std::pow(r2ik, 3) + std::pow(r2ik, 2) * (r2jk + r2ij)
                          + r2ik * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * std::pow(r2ij, 2))
                          - 5.0 * std::pow(r2jk - r2ij, 2) * (r2jk + r2ij)) / r5ijk;

        double dang_jk = -0.375 * (std::pow(r2jk, 3) + std::pow(r2jk, 2) * (r2ik + r2ij)
                          + r2jk * (3.0 * std::pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * std::pow(r2ij, 2))
                          - 5.0 * std::pow(r2ik - r2ij, 2) * (r2ik + r2ij)) / r5ijk;

        double prefactor = c9 * triple.triple_scale / 3.0;
        Eigen::Vector3d dgij = prefactor * (-dang_ij * fdmp + ang * dfdmp) / r2ij * rij_vec;
        Eigen::Vector3d dgik = prefactor * (-dang_ik * fdmp + ang * dfdmp) / r2ik * rik_vec;
        Eigen::Vector3d dgjk = prefactor * (-dang_jk * fdmp + ang * dfdmp) / r2jk * rjk_vec;

        acc.gradient.row(triple.i) += -(dgij + dgik);
        acc.gradient.row(triple.j) += (dgij - dgjk);
        acc.gradient.row(triple.k) += (dgik + dgjk);
    }
}

// ============================================================================
// Bonded ATM (BATM) for 1,4-pairs
// Reference: gfnff_engrad.F90:3267-3334 (batmgfnff_eg)
// ============================================================================

void FFWorkspace::calcBATM(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].batm_triples;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double fqq = 3.0;

    for (int idx = begin; idx < end; ++idx) {
        const auto& batm = m_batm_triples[idx];

        Eigen::Vector3d i_pos = m_geometry.row(batm.i).transpose();
        Eigen::Vector3d j_pos = m_geometry.row(batm.j).transpose();
        Eigen::Vector3d k_pos = m_geometry.row(batm.k).transpose();

        Eigen::Vector3d rij_vec = j_pos - i_pos;
        Eigen::Vector3d rik_vec = k_pos - i_pos;
        Eigen::Vector3d rjk_vec = k_pos - j_pos;

        double r2ij = rij_vec.squaredNorm();
        double r2jk = rjk_vec.squaredNorm();
        double r2ik = rik_vec.squaredNorm();

        double rij = std::sqrt(r2ij), rjk = std::sqrt(r2jk), rik = std::sqrt(r2ik);

        double rijk3 = r2ij * r2jk * r2ik;
        double mijk = -r2ij + r2jk + r2ik;
        double imjk = r2ij - r2jk + r2ik;
        double ijmk = r2ij + r2jk - r2ik;

        double ang = 0.375 * ijmk * imjk * mijk / rijk3;
        double rav3 = std::pow(rijk3, 1.5);
        double angr9 = (ang + 1.0) / rav3;

        // Phase-1 topology charges (fixed)
        double fi = std::min(std::max(1.0 - fqq * m_topology_charges(batm.i), -4.0), 4.0);
        double fj = std::min(std::max(1.0 - fqq * m_topology_charges(batm.j), -4.0), 4.0);
        double fk = std::min(std::max(1.0 - fqq * m_topology_charges(batm.k), -4.0), 4.0);

        double c9 = fi * fj * fk * batm.zb3atm_i * batm.zb3atm_j * batm.zb3atm_k;
        double energy = c9 * angr9;
        acc.energy.batm += energy;

        if (m_do_gradient) {
            double dang_ij = -0.375 * (std::pow(r2ij, 3) + std::pow(r2ij, 2) * (r2jk + r2ik)
                                  + r2ij * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * std::pow(r2ik, 2))
                                  - 5.0 * std::pow(r2jk - r2ik, 2) * (r2jk + r2ik))
                                  / (rij * rijk3 * rav3);

            double dang_jk = -0.375 * (std::pow(r2jk, 3) + std::pow(r2jk, 2) * (r2ik + r2ij)
                                  + r2jk * (3.0 * std::pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * std::pow(r2ij, 2))
                                  - 5.0 * std::pow(r2ik - r2ij, 2) * (r2ik + r2ij))
                                  / (rjk * rijk3 * rav3);

            double dang_ik = -0.375 * (std::pow(r2ik, 3) + std::pow(r2ik, 2) * (r2jk + r2ij)
                                  + r2ik * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * std::pow(r2ij, 2))
                                  - 5.0 * std::pow(r2jk - r2ij, 2) * (r2jk + r2ij))
                                  / (rik * rijk3 * rav3);

            Eigen::Vector3d dgij = -dang_ij * c9 * (rij_vec / rij);
            Eigen::Vector3d dgjk = -dang_jk * c9 * (rjk_vec / rjk);
            Eigen::Vector3d dgik = -dang_ik * c9 * (rik_vec / rik);

            acc.gradient.row(batm.j) += (-dgij + dgjk);
            acc.gradient.row(batm.k) += (-dgik - dgjk);
            acc.gradient.row(batm.i) += (dgij + dgik);
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_batm += (acc.gradient - grad_before);
}

// ============================================================================
// rev-gfnff stage 3a(ii) (Claude Generated, Sep 2026): valence share of the bond well
//   E_ij = -k_b e^{-a dr^2} w_ij c_ij,  c_ij = 1/2 (f_i + f_j),
//   f_i  = clip((Val_i - sum_{k != j} w_ik) / w_ij, 0, 1),
//   Val_i = Val_Z(i) + softplus(sum_k shareClip(b_ik) - Val_Z(i))      (b = the tight bond order)
// Pairwise wells carry no valence conservation: at an exchange transition state the well of the
// breaking and the well of the forming bond are both evaluated at full depth, so the pair sum is
// two wells where one bond's worth of valence is available. c_ij makes the two wells share it.
//
// The per-atom sums are built HERE, once per topology corner, from the corner's own bond list -
// so the stage-1b corner blend carries a change of the sums in over s exactly like every other
// neighbour re-parameterisation of a topology swap (a pair joins the list at s = 0, where its
// corner has weight 0). Val_Z is the per-element sigma valence the over-coordination term already
// uses (revValence(Z), no second table); the clip is C1 and exactly 0/1 outside [0,1].
//
// Val_i is the HYPERVALENT-CORRECT reading of that table: the nominal valence plus a saturating
// term in the settled-partner count. Since sum_k shareClip(b_ik) <= n_i (the bond count) and
// shareClip <= 1, Val_i equals Val_Z(i) EXACTLY whenever n_i <= Val_Z(i) - so a saturated
// equilibrium, an unsaturated atom and a stretched bond all keep the nominal valence and the
// term stays bit-identical - and it grows only for an atom with MORE partners than its nominal
// valence. An ammonium's fourth bond, a hydronium's third and a perchlorate's fourth are settled
// partners and are credited (c = 1 on every genuine bond); a migrating H that hangs two PARTIAL
// bonds off one valence has a settled count at or below 1, so it still shares them.
// ============================================================================

void FFWorkspace::prepareValenceShare()
{
    const int N = m_natoms;
    const int nb = static_cast<int>(m_bonds.size());
    if (N == 0 || static_cast<int>(m_rev.valence.size()) != N) {
        m_rev_share_sum.resize(0);
        m_rev_share_val.resize(0);
        m_rev_share_dval.resize(0);
        m_rev_share_g.clear();
        m_rev_share_dg.clear();
        m_rev_share_gclip.clear();
        m_rev_share_sig.clear();
        m_rev_share_dsig.clear();
        return;
    }
    if (m_rev_share_sum.size() != N) {
        m_rev_share_sum.resize(N);
        m_rev_share_val.resize(N);
        m_rev_share_dval.resize(N);
    }
    m_rev_share_sum.setZero();
    m_rev_share_val.setZero();
    m_rev_share_dval.setZero();
    m_rev_share_g.assign(nb, 1.0);
    m_rev_share_dg.assign(nb, 0.0);
    m_rev_share_gclip.assign(nb, 0.0);
    m_rev_share_sig.assign(nb, 0.0);
    m_rev_share_dsig.assign(nb, 0.0);
    // rev-gfnff stage 3a(ii) (Claude Generated, Sep 14, 2026): the smooth 1,3 proxy needs, per
    // pair, the settled weights of the pairs that share a partner with it, so the corner's bond
    // list is indexed by atom once. `adj` holds BOND INDICES (m_bonds order); the pair's own
    // weight w and settled weight sigma are computed in the same pass.
    m_rev_adj.assign(N, {});
    std::vector<double> wp(nb, 0.0);
    for (int p = 0; p < nb; ++p) {
        const Bond& b = m_bonds[p];
        if (b.i < 0 || b.j < 0 || b.i >= N || b.j >= N)
            continue;
        const Eigen::Vector3d d = m_geometry.row(b.i) - m_geometry.row(b.j);
        const double r = d.norm();
        if (r < 1e-8)
            continue;
        m_rev_adj[b.i].push_back(p);
        m_rev_adj[b.j].push_back(p);
        double dbdr = 0.0;
        const double bo = revOrder(b.i, b.j, r, &dbdr);
        m_rev_share_sig[p] = shareSettled(bo);
        m_rev_share_dsig[p] = shareSettledD(bo) * dbdr;
        wp[p] = revWeight(b.i, b.j, r, nullptr);
    }
    // Pass 2: the genuineness g_p of every pair - g = shareClip(1 - t) with
    //     t_p = sum_{k != i,j} sigma_ik sigma_jk
    // the bond-order leak of a SETTLED shared partner (k is a neighbour of BOTH ends). The two
    // ends themselves are never counted: a bond index of i's list whose other end is j cannot be
    // matched by an entry of j's list (whose other ends are != j). With the proxy off, g = 1 and
    // t = 0 for every pair, i.e. exactly the unmasked share.
    const bool onethree = m_rev.share_onethree;
    if (onethree) {
        if (static_cast<int>(m_rev_share_stamp.size()) != N)
            m_rev_share_stamp.assign(N, -1);
        else
            std::fill(m_rev_share_stamp.begin(), m_rev_share_stamp.end(), -1);
        for (int p = 0; p < nb; ++p) {
            const Bond& b = m_bonds[p];
            if (b.i < 0 || b.j < 0 || b.i >= N || b.j >= N)
                continue;
            // stamp[k] = the BOND INDEX of the pair (i, k) while k is a neighbour of i, -1
            // otherwise. The entries are cleared again after the pair, so the array is all -1
            // between pairs (the invariant the membership test relies on).
            for (int qa : m_rev_adj[b.i])
                m_rev_share_stamp[otherEndOf(qa, b.i)] = qa;
            double t = 0.0;
            for (int qb : m_rev_adj[b.j]) {
                const int k = otherEndOf(qb, b.j);
                const int qa = m_rev_share_stamp[k];
                if (qa >= 0)
                    t += m_rev_share_sig[qa] * m_rev_share_sig[qb];
            }
            for (int qa : m_rev_adj[b.i])
                m_rev_share_stamp[otherEndOf(qa, b.i)] = -1;
            const double leak = 1.0 - t;
            m_rev_share_g[p] = shareClip(leak);
            m_rev_share_gclip[p] = -shareClipD(leak);
        }
    }
    // Pass 3: the per-atom sums. sum_i = sum_k w_ik g_ik - the CLAIM of i's partners on its
    // valence, which a 1,3 contact (g = 0) does not make. The settled count that drives the
    // effective valence below stays over ALL partners: sigma is exactly 0 for a 1,3 contact of
    // two ends that are not themselves settled (the F...F pairs of a compressed BF4-), and where
    // it is not, the partner genuinely carries bond order toward both ends.
    for (int p = 0; p < nb; ++p) {
        const Bond& b = m_bonds[p];
        if (b.i < 0 || b.j < 0 || b.i >= N || b.j >= N)
            continue;
        const double claim = wp[p] * m_rev_share_g[p];
        m_rev_share_sum(b.i) += claim;
        m_rev_share_sum(b.j) += claim;
        m_rev_share_val(b.i) += m_rev_share_sig[p];
        m_rev_share_val(b.j) += m_rev_share_sig[p];
    }
    // Pass 4: the effective valence. Val_Z + G(settled - Val_Z) is exactly the nominal valence
    // for every atom with at most Val_Z partners (a saturated atom, a stretched bond, a radical)
    // and only grows for an atom that carries MORE partners than its nominal valence - the only
    // case in which the share can bite at all.
    // Claude Generated (Sep 15, 2026): -gfnff.rev_budget_fix_h true exempts HYDROGEN from that
    // growth. The softplus is right for a hypervalent centre but wrong for an H: a hydrogen has one
    // valence, and a bridging H (a just-formed H2 still bonded to its carbon) is a 3c-2e bond whose
    // two partial wells must SHARE that one valence. Without the exemption its budget reaches 2 as
    // soon as the second partner's tight bond order crosses the settled window, and BOTH wells jump
    // from half share to full share inside one step with no topology event - the measured origin of
    // the hot react-MD blow-ups (see the RevSettings member documentation). The derivative channel
    // is zeroed with it: m_rev_share_dval is read only as the dVal_i/d(settled) factor of
    // applyValenceShareGradient, so a constant budget must contribute exactly nothing there.
    const bool fix_h = m_rev.budget_fix_h && static_cast<int>(m_atom_types.size()) == N;
    // Claude Generated (Sep 18, 2026): the "conserving" share of FABLE_REVIEW_2 A.5 replaces this
    // whole pass (see RevSettings::share_conserving). Its budget is built from the SAME wide sum
    // S_i that the share divides by, not from the settled count - that is what makes NH4+ / H3O+
    // come out at exactly Val = S, i.e. f = 1 and a bit-identical energy, instead of the +0.17 /
    // +0.05 kcal/mol the review measured offline with the settled-count budget.
    const bool conserving = m_rev.share_conserving;
    if (conserving) {
        prepareConservingShare(fix_h);
    } else {
        m_rev_share_f.resize(0);
        m_rev_share_dfdS.resize(0);
        m_rev_share_cap.resize(0);
        for (int i = 0; i < N; ++i) {
            if (fix_h && m_atom_types[i] == 1) {
                m_rev_share_dval(i) = 0.0;
                m_rev_share_val(i) = m_rev.valence[i];
                continue;
            }
            const double x = m_rev_share_val(i) - m_rev.valence[i];
            m_rev_share_dval(i) = shareExcessD(x);
            m_rev_share_val(i) = m_rev.valence[i] + shareExcess(x);
        }
    }
    // Claude Generated (Sep 14, 2026): CURCUMA_SHAREDUMP=1 prints the per-pair share table of the
    // corner that is being evaluated - the only place the 1,3 proxy and the share's argument are
    // visible (same convention as CURCUMA_NBDIAG / CURCUMA_BONDDUMP / CURCUMA_REVDUMP: zero cost
    // when unset). Columns: pair, r, w, tight b, g, sum/Val at both ends, the share arguments and
    // the resulting c.
    if (const char* d = std::getenv("CURCUMA_SHAREDUMP"); d && d[0] == '1') {
        CurcumaLogger::result(fmt::format(
            "share dump: corner with {} bonds, onethree {}", nb, onethree ? 1 : 0));
        for (int p = 0; p < nb; ++p) {
            const Bond& b = m_bonds[p];
            if (b.i < 0 || b.j < 0 || b.i >= N || b.j >= N)
                continue;
            const double r = (m_geometry.row(b.i) - m_geometry.row(b.j)).norm();
            const double w = revWeight(b.i, b.j, r, nullptr);
            const double bo = revOrder(b.i, b.j, r, nullptr);
            double c1, c2, fi, fj, cshare;
            if (conserving) {
                // u = Val/S (the argument of the smooth min), f = the per-atom share factor
                c1 = m_rev_share_sum(b.i) > 1e-12 ? m_rev_share_val(b.i) / m_rev_share_sum(b.i) : 1.0;
                c2 = m_rev_share_sum(b.j) > 1e-12 ? m_rev_share_val(b.j) / m_rev_share_sum(b.j) : 1.0;
                fi = m_rev_share_f(b.i);
                fj = m_rev_share_f(b.j);
                cshare = 1.0 - m_rev_share_g[p] * (1.0 - fi * fj);
            } else {
                c1 = (m_rev_share_val(b.i) - m_rev_share_sum(b.i) + w * m_rev_share_g[p])
                     / (w > 1e-12 ? w : 1e-12);
                c2 = (m_rev_share_val(b.j) - m_rev_share_sum(b.j) + w * m_rev_share_g[p])
                     / (w > 1e-12 ? w : 1e-12);
                fi = shareClip(c1);
                fj = shareClip(c2);
                cshare = 1.0 - m_rev_share_g[p] * (1.0 - 0.5 * (fi + fj));
            }
            CurcumaLogger::result(fmt::format(
                "share {:3d} {:3d}-{:3d} r {:8.4f} w {:7.4f} b {:7.4f} sig {:7.4f} g {:7.4f} "
                "sum {:7.4f}/{:7.4f} Val {:7.4f}/{:7.4f} u {:7.4f}/{:7.4f} f {:7.4f}/{:7.4f} c {:7.4f}",
                p, b.i + 1, b.j + 1, r, w, bo,
                m_rev_share_sig[p], m_rev_share_g[p],
                m_rev_share_sum(b.i), m_rev_share_sum(b.j),
                m_rev_share_val(b.i), m_rev_share_val(b.j), c1, c2, fi, fj, cshare));
        }
        if (conserving && m_rev_share_cap.size() == N) {
            for (int i = 0; i < N; ++i)
                CurcumaLogger::result(fmt::format(
                    "shareA {:3d} Z {:3d} S {:9.5f} ValZ {:6.3f} cap {:7.4f} Val {:9.5f} "
                    "f {:9.6f} dfdS {:12.6e}",
                    i + 1, (static_cast<int>(m_atom_types.size()) == N ? m_atom_types[i] : 0),
                    m_rev_share_sum(i), m_rev.valence[i], m_rev_share_cap(i),
                    m_rev_share_val(i), m_rev_share_f(i), m_rev_share_dfdS(i)));
        }
    }
}

// ============================================================================
// rev-gfnff stage 3a(ii), the "conserving" share (Claude Generated, Sep 18, 2026)
// FABLE_REVIEW_2 A.5. Per atom:
//     S_i    = sum_k w_ik g_ik                       (built by Pass 3 above)
//     X_i    = the excess-budget cap, element/charge rule (capForAtom below)
//     Val_i  = Val_Z + min(G(S_i - Val_Z), X_i)      G = shareExcess, min = shareMinOne
//     f_i    = min(1, Val_i / S_i)                   shareMinOne again
//     c_ij   = 1 - g_ij (1 - f_i f_j)                = f_i f_j for an ordinary pair (g = 1)
// and d f_i / d S_i is stored, because EVERY geometry dependence of this mode runs through the
// term weights w that build S_i - the budget included. That is why m_rev_share_dval is 0 here:
// the settled-count channel of the delivered rule has no counterpart, and the existing Lambda
// pass over dw/dr carries the whole chain rule (applyValenceShareGradient).
// ============================================================================
// ============================================================================
// rev-gfnff stage 3a(iii): the bond-well forms (Claude Generated, Sep 18, 2026)
// FABLE_REVIEW_2 B. Both new forms are E = -D (2y - y^2) with E(r0) = -D and E'(r0) = 0
// identically, and both are curvature-pinned to the delivered Gaussian's 2 alpha |k_b|, so a
// (MG) resp. u (erf-Morse) follows from the depth. Only the depth scale s = D/|k_b| and the tail
// are fitted, per element pair (rev_well_table.h, from scripts/revgfnff_wellfit.py).
//
// The bisection for u is the erf-Morse form's extra setup cost (the MG form has a closed form)
// and is the reason this runs once per corner, on the main thread, instead of per energy call.
// ============================================================================
double FFWorkspace::wellErfMorseU(double D, double sigma, double K)
{
    // g(u) = (2/(sigma sqrt(pi))) e^{-z^2}/erfc(z) at z = -u/sigma is |y'(0)| and decreases
    // monotonically in u; the pinned curvature 2 D g^2 = K gives g = sqrt(K/(2D)).
    const double tgt = std::sqrt(K / (2.0 * D));
    auto g = [&](double u) {
        const double z = -u / sigma;
        if (z > 25.0)
            return 2.0 * z / sigma;                     // erfc asymptotics, avoids 0/0
        return 2.0 / (sigma * 1.7724538509055159) * std::exp(-z * z) / std::erfc(z);
    };
    double lo = -40.0 * sigma, hi = 8.0 * sigma;
    if (g(lo) < tgt || g(hi) > tgt)
        return std::numeric_limits<double>::quiet_NaN();
    for (int it = 0; it < 80; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (g(mid) > tgt)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

void FFWorkspace::prepareWellForms()
{
    const int nb = static_cast<int>(m_bonds.size());
    if (!m_rev.enabled || m_rev.well_form == 0 || nb == 0) {
        m_rev_well.clear();
        m_rev_well_stamp.clear();
        m_rev_well_stamp_form = -1;
        return;
    }
    // Cache: the inputs are per-bond constants (fc, exponent and the two element numbers), so
    // the whole pass - and in particular the erf-Morse bisection - only has to run when the bond
    // list or its parameters change, not on every energy call.
    bool same = (m_rev_well_stamp_form == m_rev.well_form)
                && m_rev_well_stamp.size() == static_cast<size_t>(nb)
                && m_rev_well.size() == static_cast<size_t>(nb);
    if (same) {
        for (int p = 0; p < nb; ++p) {
            const Bond& b = m_bonds[p];
            const auto& st = m_rev_well_stamp[p];
            if (st[0] != b.fc || st[1] != b.exponent
                || st[2] != static_cast<double>(b.z_i) || st[3] != static_cast<double>(b.z_j)) {
                same = false;
                break;
            }
        }
    }
    if (same)
        return;
    m_rev_well_stamp.resize(nb);
    m_rev_well_stamp_form = m_rev.well_form;
    for (int p = 0; p < nb; ++p) {
        const Bond& b = m_bonds[p];
        m_rev_well_stamp[p] = { b.fc, b.exponent, static_cast<double>(b.z_i), static_cast<double>(b.z_j) };
    }
    m_rev_well.assign(nb, RevWellPar{});
    // The table is in Angstrom units (that is how the class-A fit reports it and how the header
    // reads); the workspace works in Bohr. beta [1/A^2] -> [1/Bohr^2] multiplies by (A/Bohr)^2
    // and sigma [A] -> [Bohr] divides by it.
    constexpr double kBohrPerAng = 1.8897261246257702;
    std::vector<std::pair<int, int>> missing;
    for (int p = 0; p < nb; ++p) {
        const Bond& b = m_bonds[p];
        if (b.z_i <= 0 || b.z_j <= 0)
            continue;
        const RevWellTable::Entry* e = RevWellTable::find(b.z_i, b.z_j);
        if (!e) {
            const int z1 = std::min(b.z_i, b.z_j), z2 = std::max(b.z_i, b.z_j);
            if (std::find(missing.begin(), missing.end(), std::make_pair(z1, z2)) == missing.end())
                missing.emplace_back(z1, z2);
            continue;                                   // form stays 0 -> Gaussian for this bond
        }
        const double kb = std::abs(b.fc);
        if (kb <= 0.0)
            continue;
        // The HB alpha modulation (egbond_hb) is DELIBERATELY not applied to the new forms: it
        // would enter through a resp. u and its chain rule (dE/d hb_cn_H) has no counterpart
        // there. Using alpha_orig keeps the well and its gradient exactly consistent; the cost is
        // that a hydrogen-bond donor's X-H bond is not softened in these forms. Documented, not
        // measured against a reference (the class-A set has no hydrogen bond).
        const double alpha = b.exponent;
        const double K = 2.0 * alpha * kb;
        if (m_rev.well_form == 1) {
            const double D = e->mg_s * kb;
            m_rev_well[p] = RevWellPar{ D, std::sqrt(alpha / e->mg_s),
                                        e->mg_beta / (kBohrPerAng * kBohrPerAng), 1 };
        } else {
            const double D = e->em_s * kb;
            const double sigma = e->em_sigma * kBohrPerAng;
            const double u = wellErfMorseU(D, sigma, K);
            if (std::isfinite(u))
                m_rev_well[p] = RevWellPar{ D, u, sigma, 2 };
        }
    }
    if (!missing.empty() && CurcumaLogger::get_verbosity() >= 2) {
        std::string s;
        for (const auto& m : missing)
            s += fmt::format("{}{}-{}", s.empty() ? "" : ", ", m.first, m.second);
        CurcumaLogger::warn(fmt::format(
            "rev-gfnff well form: no class-A parameters for element pair(s) {} - those bonds keep "
            "the delivered Gaussian", s));
    }
}

void FFWorkspace::prepareConservingShare(bool fix_h)
{
    const int N = m_natoms;
    if (m_rev_share_f.size() != N) {
        m_rev_share_f.resize(N);
        m_rev_share_dfdS.resize(N);
        m_rev_share_cap.resize(N);
    }
    const bool have_types = static_cast<int>(m_atom_types.size()) == N;
    const bool have_q = m_topology_charges.size() == N;
    // Q_i = the topological (phase-1 EEQ) charge of atom i plus that of its H partners IN THIS
    // CORNER. It is what separates NH4+ from NH3 + H: the +1 of an ammonium sits on the N-H4
    // group as a whole (the N itself is negative), and a neutral CH4 + H group carries 0 however
    // close the radical is. Built from the corner's own bond list, so it is a constant inside a
    // corner and a change is carried by the s-blend, exactly like every other corner quantity.
    std::vector<double> qgroup(N, 0.0);
    if (have_q) {
        for (int i = 0; i < N; ++i)
            qgroup[i] = m_topology_charges(i);
        if (have_types) {
            for (int p = 0; p < static_cast<int>(m_bonds.size()); ++p) {
                const Bond& b = m_bonds[p];
                if (b.i < 0 || b.j < 0 || b.i >= N || b.j >= N)
                    continue;
                if (m_atom_types[b.j] == 1)
                    qgroup[b.i] += m_topology_charges(b.j);
                if (m_atom_types[b.i] == 1)
                    qgroup[b.j] += m_topology_charges(b.i);
            }
        }
    }
    const double a = m_rev.share_min_width > 1e-6 ? m_rev.share_min_width : 1e-6;
    for (int i = 0; i < N; ++i) {
        const double S = m_rev_share_sum(i);
        const double valz = m_rev.valence[i];
        const int Z = have_types ? m_atom_types[i] : 0;
        // --- the cap X_i ---------------------------------------------------------------
        double cap = 0.0;
        bool delivered_growth = false;   // transition metals, see below
        if (Z == 1 || Z == 9) {
            cap = 0.0;                   // H is never hypervalent; F is never hypervalent
        } else {
            const int grp = (Z >= 1 && Z <= 86) ? GFNFFParameters::periodic_group[Z - 1] : 0;
            const int period = Z <= 2 ? 1 : Z <= 10 ? 2 : Z <= 18 ? 3 : Z <= 36 ? 4 : Z <= 54 ? 5 : 6;
            if (grp < 0) {
                // A d-block element: periodic_group is negative for it, and FABLE_REVIEW_2 A.5
                // does not cover metals. Their coordination numbers routinely exceed any sigma
                // valence, so the charge rule would scale EVERY metal-ligand well by Val_Z/CN.
                // Deliberate, documented carve-out: a transition metal keeps the delivered
                // growth (cap = the softplus itself), i.e. this mode does not touch it. NOT
                // measured - no metal is in any rev-gfnff reference set.
                delivered_growth = true;
            } else if (grp == 3) {
                // NB: GFNFFParameters::periodic_group uses MAIN-GROUP numbering 1-8, not IUPAC
                // 1-18 - so "group 13" of FABLE_REVIEW_2 A.5 (B, Al, Ga, In, Tl) is 3 here and
                // "groups 15-17" (N/P/As.., O/S/Se.., F/Cl/Br..) are 5-7. Measured the hard way:
                // with the IUPAC numbers no rule ever fired and ClO4- fell to the charge rule
                // (cap 0.65 instead of 5), costing +164 kcal/mol.
                cap = 1.0;               // the empty orbital: BF4-, BH4-, AlCl4-, H3N-BH3
            } else if (period >= 3 && grp >= 5 && grp <= 7) {
                cap = 6.0 - valz;        // the octet expansion the valence table already grants P/S
            } else {
                cap = shareClip(qgroup[i]);   // C, N, O and the rest: granted by charge
            }
        }
        if (fix_h && Z == 1)
            cap = 0.0;
        m_rev_share_cap(i) = delivered_growth ? 99.0 : cap;
        // --- the budget ----------------------------------------------------------------
        const double exc_arg = S - valz;
        const double G = shareExcess(exc_arg);
        const double dG = shareExcessD(exc_arg);          // d G / d S
        double val, dval_dS;
        if (delivered_growth) {
            val = valz + G;
            dval_dS = dG;
        } else if (cap <= 1e-12) {
            val = valz;
            dval_dS = 0.0;
        } else {
            // min(G, cap) = cap * shareMinOne(G/cap): exactly G below cap (1 - a), exactly cap
            // above it, C1 in between - and, crucially, EXACTLY cap once the softplus has
            // saturated, which is what makes an ammonium's Val equal its S to the last bit.
            const double t = G / cap;
            val = valz + cap * shareMinOne(t, a);
            dval_dS = shareMinOneD(t, a) * dG;
        }
        m_rev_share_val(i) = val;
        m_rev_share_dval(i) = 0.0;       // the settled-count channel does not exist in this mode
        // --- the share factor ----------------------------------------------------------
        if (S <= 1e-12) {
            m_rev_share_f(i) = 1.0;
            m_rev_share_dfdS(i) = 0.0;
            continue;
        }
        const double x = val / S;
        m_rev_share_f(i) = shareMinOne(x, a);
        // d f/d S = f'(x) * (val' S - val)/S^2
        m_rev_share_dfdS(i) = shareMinOneD(x, a) * (dval_dS * S - val) / (S * S);
    }
}

void FFWorkspace::applyValenceShareGradient()
{
    if (m_dEdshare_total.size() != m_natoms || m_natoms == 0)
        return;
    if (m_rev_share_val.size() != m_natoms || m_rev_share_dval.size() != m_natoms)
        return;
    // Lambda_i = sum over i's pairs of dE_pair/d sum_i. Every sum depends on the geometry only
    // through the term weights of the atom's own pairs, so the second pass is the same dw/dr
    // switch over the same bond list - d sum_i/dx = sum_k (dw_ik/dr) d r_ik/dx.
    // The effective valence adds a second chain: dE/dVal_i = -Lambda_i (the pair carries
    // +1/2 K g f_i' for Val_i and -1/2 K g f_i' for sum_i, i.e. the two are exact negatives), and
    // dVal_i/dr = dG/dx * shareSettledD(b_ik) * db_ik/dr, so the same bond list carries it too.
    const bool share_masked = m_rev_share_g.size() == m_bonds.size();
    for (int p = 0; p < static_cast<int>(m_bonds.size()); ++p) {
        const Bond& b = m_bonds[p];
        if (b.i < 0 || b.j < 0 || b.i >= m_natoms || b.j >= m_natoms)
            continue;
        const Eigen::Vector3d d = m_geometry.row(b.i) - m_geometry.row(b.j);
        const double r = d.norm();
        if (r < 1e-8)
            continue;
        double dwdr = 0.0;
        revWeight(b.i, b.j, r, &dwdr);
        double dbdr = 0.0;
        const double bo = revOrder(b.i, b.j, r, &dbdr);
        // dE/d sum: the term weight of this pair feeds both atoms' sums, and the pair may only be
        // CLAIMED in proportion to its genuineness, so the whole channel carries g_p (g = 1 for an
        // ordinary bond, and the pair's own g enters its own f through the credit-back term, which
        // calcBonds' dcdw carries). Verified against the FD of the energy with g = 1 and g < 1:
        // without the factor the share's sum channel is short by (1 - g_p) times the atoms' claims.
        double dEdr_sum = 0.0;
        if (dwdr != 0.0) {
            const double coeff = share_masked ? m_rev_share_g[p]
                                                    * (m_dEdshare_total(b.i) + m_dEdshare_total(b.j))
                                              : (m_dEdshare_total(b.i) + m_dEdshare_total(b.j));
            dEdr_sum = coeff * dwdr;
        }
        // dE/d Val: only the SETTLED part of the pair changes Val (shareSettledD is 0 outside the
        // unit interval of b, so a pair that is fully or barely bonded contributes nothing).
        double dEdr_val = 0.0;
        if (dbdr != 0.0) {
            const double scd = shareSettledD(bo);
            if (scd != 0.0)
                dEdr_val = -(m_dEdshare_total(b.i) * m_rev_share_dval(b.i)
                             + m_dEdshare_total(b.j) * m_rev_share_dval(b.j)) * scd * dbdr;
        }
        const double dEdr = dEdr_sum + dEdr_val;
        if (dEdr == 0.0)
            continue;
        const double f = dEdr / r;
        const Eigen::Vector3d g = f * d;
        m_result_gradient.row(b.i) += g.transpose();
        m_result_gradient.row(b.j) -= g.transpose();
    }
    // rev-gfnff stage 3a(ii) (Claude Generated, Sep 14, 2026): the THREE-BODY chain rule of the
    // smooth 1,3 proxy. g_p = shareClip(1 - t_p) with t_p = sum_k sigma_ik sigma_jk, so
    //     dg_p/dx = (dg_p/dt_p) * sum_k [ sigma_jk (d sigma_ik/dr_ik) dr_ik/dx
    //                                     + sigma_ik (d sigma_jk/dr_jk) dr_jk/dx ],
    // i.e. every shared settled partner k of a pair pulls on the two bonds i-k and j-k. The
    // coefficient of g_p in the energy is the pair's own dE/dg (accumulated by calcBonds) PLUS the
    // second half of the same derivative - g_p also feeds the SUMS of its two atoms, which the
    // sum chain above carries as dE/d sum_i = dEdshare(i) - so
    //     Lambda_p = dE/dg_p + w_p (dEdshare(i) + dEdshare(j)).
    // Same place and the same per-corner bond list as the sum chain, so the corner blend carries
    // it identically.
    {
        bool any = false;
        for (const double v : m_rev_share_dg)
            if (v != 0.0) { any = true; break; }
        if (any && m_rev_share_gclip.size() == m_bonds.size()
            && m_rev_share_stamp.size() == static_cast<size_t>(m_natoms)) {
            std::fill(m_rev_share_stamp.begin(), m_rev_share_stamp.end(), -1);
            auto addRadial = [&](int a, int c, double dEdr) {
                const Eigen::Vector3d d = m_geometry.row(a) - m_geometry.row(c);
                const double r = d.norm();
                if (r < 1e-8 || dEdr == 0.0)
                    return;
                const double f = dEdr / r;
                m_result_gradient.row(a) += (f * d).transpose();
                m_result_gradient.row(c) -= (f * d).transpose();
            };
            for (int p = 0; p < static_cast<int>(m_bonds.size()); ++p) {
                const Bond& b = m_bonds[p];
                if (b.i < 0 || b.j < 0 || b.i >= m_natoms || b.j >= m_natoms)
                    continue;
                const Eigen::Vector3d d = m_geometry.row(b.i) - m_geometry.row(b.j);
                const double r = d.norm();
                if (r < 1e-8)
                    continue;
                const double wp = revWeight(b.i, b.j, r, nullptr);
                const double lam = m_rev_share_dg[p]
                                   + wp * (m_dEdshare_total(b.i) + m_dEdshare_total(b.j));
                const double cd = m_rev_share_gclip[p] * lam;
                if (cd == 0.0)
                    continue;
                for (int qa : m_rev_adj[b.i])
                    m_rev_share_stamp[otherEndOf(qa, b.i)] = qa;
                for (int qb : m_rev_adj[b.j]) {
                    const int k = otherEndOf(qb, b.j);
                    const int qa = m_rev_share_stamp[k];
                    if (qa < 0)
                        continue;
                    // d sigma_ik / dr_ik dr_ik/dx  and  d sigma_jk / dr_jk dr_jk/dx
                    addRadial(b.i, k, cd * m_rev_share_sig[qb] * m_rev_share_dsig[qa]);
                    addRadial(b.j, k, cd * m_rev_share_sig[qa] * m_rev_share_dsig[qb]);
                }
                for (int qa : m_rev_adj[b.i])
                    m_rev_share_stamp[otherEndOf(qa, b.i)] = -1;
            }
        }
    }
}

// ============================================================================
// rev-gfnff stage 1: over-coordination energy (Claude Generated, Sep 2026)
//   bo_sum_i = sum_j b_ij BO_ij over every pair that carries the blended repulsion
//   E_over  = sum_i p_i sp(bo_sum_i - Val_i)^2,  sp = softplus (rev_bond_order.h)
// Two passes over the pairs: the sums first, then the gradient through d b_ij / dr.
// Runs on the main thread after the partitions (needs the complete sums).
// ============================================================================

void FFWorkspace::calcOverCoordination(bool gradient)
{
    const int N = m_natoms;
    if (N == 0 || static_cast<int>(m_rev.over_p.size()) != N || static_cast<int>(m_rev.valence.size()) != N)
        return;
    m_rev_bo_sum = Vector::Zero(N);
    auto pair_r = [&](int i, int j, Eigen::Vector3d* vec) {
        Eigen::Vector3d ri = m_geometry.row(i), rj = m_geometry.row(j);
        Eigen::Vector3d d = ri - rj;
        if (vec) *vec = d;
        return d.norm();
    };
    auto accumulate = [&](const std::vector<GFNFFRepulsion>& list) {
        for (const auto& rep : list) {
            if (!rep.blend) continue;
            const double r = pair_r(rep.i, rep.j, nullptr);
            if (r < 1e-8) continue;
            const double b = revOrder(rep.i, rep.j, r, nullptr) * rep.bo_mult;
            m_rev_bo_sum(rep.i) += b;
            m_rev_bo_sum(rep.j) += b;
        }
    };
    accumulate(m_bonded_reps);
    accumulate(m_nonbonded_reps);

    double e_over = 0.0;
    Vector dEdsum = Vector::Zero(N); // dE_over / d bo_sum_i
    for (int i = 0; i < N; ++i) {
        double dsp = 0.0;
        const double sp = RevGFNFF::softplus(m_rev_bo_sum(i) - m_rev.valence[i] - m_rev.over_shift, m_rev.over_k, &dsp);
        e_over += m_rev.over_p[i] * sp * sp;
        dEdsum(i) = 2.0 * m_rev.over_p[i] * sp * dsp;
    }
    m_result_energy.over_coord += e_over;
    if (!gradient)
        return;
    auto apply = [&](const std::vector<GFNFFRepulsion>& list) {
        for (const auto& rep : list) {
            if (!rep.blend) continue;
            Eigen::Vector3d vec;
            const double r = pair_r(rep.i, rep.j, &vec);
            if (r < 1e-8) continue;
            double dbdr = 0.0;
            revOrder(rep.i, rep.j, r, &dbdr);
            const double f = (dEdsum(rep.i) + dEdsum(rep.j)) * rep.bo_mult * dbdr / r;
            {
                Eigen::Vector3d g = f * vec;
                m_result_gradient.row(rep.i) += g.transpose();
                m_result_gradient.row(rep.j) -= g.transpose();
            }
        }
    };
    apply(m_bonded_reps);
    apply(m_nonbonded_reps);
}

// ============================================================================
// rev-gfnff stage 2: bond hardness of the split-charge model (Claude Generated, Sep 2026)
//   E_sqe   = sum_(ij) 1/2 kappa_ij p_ij^2,   kappa_ij(b) = kappa0_ij / b_ij(r)
//   dE/dr   = 1/2 p_ij^2 dkappa/db db/dr = -1/2 p_ij^2 kappa0_ij / b_ij^2 db_ij/dr
// E is variational in p (the SQE system is solved to dE/dp = 0 before the kernel runs), so
// the only explicit r-dependence left is b_ij — the p_ij are held fixed here. b is the same
// bond order the over-coordination term uses (RevSettings::R2 / bo2_width).
// Runs on the main thread after the partitions, next to calcOverCoordination().
// docs/REV_GFNFF_STAGE2.md
// ============================================================================

void FFWorkspace::calcSqeHardness(bool gradient)
{
    if (m_sqe_pairs.empty() || m_natoms == 0)
        return;
    const double bmin = std::max(m_sqe_bmin, 1e-12);
    double e_sqe = 0.0;
    for (const SqePairData& sp : m_sqe_pairs) {
        if (sp.i < 0 || sp.j < 0 || sp.i >= m_natoms || sp.j >= m_natoms)
            continue;
        if (sp.kappa0 == 0.0 || sp.p == 0.0)
            continue;
        Eigen::Vector3d ri = m_geometry.row(sp.i), rj = m_geometry.row(sp.j);
        const Eigen::Vector3d d = ri - rj;
        const double r = d.norm();
        if (r < 1e-8)
            continue;
        double dbdr = 0.0;
        const double b_raw = revOrder(sp.i, sp.j, r, gradient ? &dbdr : nullptr);
        // Below the floor the pair is rigid: the solver has already forced p = 0 there, and
        // clamping keeps kappa (and its derivative) finite if a pair drifts out mid-step.
        const bool clamped = (b_raw <= bmin);
        const double b = clamped ? bmin : b_raw;
        const double pp = sp.p * sp.p;
        e_sqe += 0.5 * sp.kappa0 * pp / b;
        if (!gradient || clamped)
            continue;
        const double dEdr = -0.5 * pp * sp.kappa0 / (b * b) * dbdr;
        const Eigen::Vector3d g = (dEdr / r) * d;
        m_result_gradient.row(sp.i) += g.transpose();
        m_result_gradient.row(sp.j) -= g.transpose();
    }
    m_result_energy.sqe_hardness += e_sqe;
}
