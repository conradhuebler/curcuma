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
#include "src/tools/pbc_utils.h"  // PBCUtils::OrthoCell (Nov 2026, A1)
#include "ewald.h"                // curcuma::ewald::EwaldSum, periodic GFN-FF Coulomb (Nov 2026, B2)

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <Eigen/Dense>
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
    // Claude Generated (Nov 2026, A3c): homogeneous long-range D4 tail, periodic cell only (0 otherwise);
    // set in postProcess() from setDispersionTail(), not by a kernel
    double disp_tail = 0;
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
        disp_tail = 0;
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
        disp_tail += o.disp_tail;
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
               atm + batm + stors + disp_tail;  // disp_tail last: + 0.0 leaves the sum unchanged without a cell
    }
};

/**
 * @brief Atomic virial tensor of the GFN-FF energy per term (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED N1)
 *
 * Xi_ab = -dE/d(eps_ab) under an affine strain of the cell AND all atoms at fixed fractional coordinates
 * (atomic virial: molecules are deformed with the cell). Each kernel contracts the positions it actually used
 * (its local minimum-image frame, see placeMinimumImage) with the gradient rows it wrote: Xi_ab -= x_a G_b, summed
 * over the atoms of the term; the sum of a term's gradient rows is zero, so the frame origin drops out. Pair terms
 * give Xi_ab -= d_a G_{j,b} with d = x_j - x_i (minimum image). The parts that are not a kernel sum (CN chain rule,
 * Ewald background terms, D4 tail) are added in postProcess(). Units: Eh. Only filled when setWantVirial(true);
 * the default path does not read or write it (kernels are instantiated without the virial code).
 *
 * V-LRC (Nov 2026): `tail_lrc` is NOT part of -dE/d(eps) of the model energy. It is the classical long-range
 * pressure correction of the pairs that cross the sharp D4 evaluation radius r_c under a strain (statistical
 * estimator for g(r) = 1 beyond r_c, derivation at FFWorkspace::postProcessVirial). total() stays the model
 * derivative (finite-difference gate pbc_virial); pressureTotal() = total() + tail_lrc is what the reported
 * pressure uses (tail_lrc is zero with pressure_lrc false).
 */
struct FFVirialComponents {
    Eigen::Matrix3d bond, angle, dihedral, inversion, stors;
    Eigen::Matrix3d dispersion;   ///< D4 (and D3) pair kernel
    Eigen::Matrix3d disp_tail;    ///< homogeneous tail: E_tail * I (fixed r_c, E_tail ~ 1/V)
    Eigen::Matrix3d tail_lrc;     ///< V-LRC: (3 E6_tail + 5 E8_tail)/3 * I, crossing-pair correction (not a model derivative)
    Eigen::Matrix3d repulsion;    ///< bonded + non-bonded repulsion
    Eigen::Matrix3d coulomb;      ///< Coulomb pair term (plain, or Ewald real + reciprocal)
    Eigen::Matrix3d coulomb_bg;   ///< Ewald neutralising and Gaussian background: (E_bg + E_gbg) * I
    Eigen::Matrix3d hbond, xbond, batm, atm;
    Eigen::Matrix3d cn_bond, cn_dispersion, cn_coulomb; ///< CN chain rule, attributed as the gradient components are

    void reset()
    {
        for (Eigen::Matrix3d* m : { &bond, &angle, &dihedral, &inversion, &stors, &dispersion, &disp_tail, &tail_lrc,
                 &repulsion, &coulomb, &coulomb_bg, &hbond, &xbond, &batm, &atm, &cn_bond, &cn_dispersion, &cn_coulomb })
            m->setZero();
    }
    FFVirialComponents& operator+=(const FFVirialComponents& o)
    {
        bond += o.bond; angle += o.angle; dihedral += o.dihedral; inversion += o.inversion; stors += o.stors;
        dispersion += o.dispersion; disp_tail += o.disp_tail; tail_lrc += o.tail_lrc; repulsion += o.repulsion;
        coulomb += o.coulomb; coulomb_bg += o.coulomb_bg; hbond += o.hbond; xbond += o.xbond;
        batm += o.batm; atm += o.atm; cn_bond += o.cn_bond; cn_dispersion += o.cn_dispersion; cn_coulomb += o.cn_coulomb;
        return *this;
    }
    /// Per-energy-term totals (kernel part plus its CN chain-rule part), matching the FFEnergyComponents fields
    Eigen::Matrix3d bondTotal() const { return bond + cn_bond; }
    Eigen::Matrix3d dispersionTotal() const { return dispersion + disp_tail + cn_dispersion; }
    Eigen::Matrix3d coulombTotal() const { return coulomb + coulomb_bg + cn_coulomb; }
    /// -dE/d(eps) of the model energy (sum of the terms above; tail_lrc NOT included)
    Eigen::Matrix3d total() const
    {
        return bondTotal() + angle + dihedral + inversion + stors + dispersionTotal() + repulsion + coulombTotal()
            + hbond + xbond + batm + atm;
    }
    /// V-LRC: the virial of the reported pressure, total() + tail_lrc (bitwise total() when tail_lrc is zero)
    Eigen::Matrix3d pressureTotal() const { return total() + tail_lrc; }
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
    Vector dEdcn;             ///< N dE/dCN for chain-rule
    Vector dEdcn_bond;        ///< N bond-only dE/dCN for per-component attribution
    FFEnergyComponents energy;
    FFTermTimings timings;    ///< Claude Generated (May 2026): per-term ms (this partition only)

    // Optional per-component gradients (only allocated if store_components=true)
    GeoGradMatrix grad_bond, grad_angle, grad_torsion, grad_repulsion;
    GeoGradMatrix grad_coulomb, grad_dispersion, grad_hb, grad_xb, grad_batm, grad_atm;

    bool has_components = false;

    /// Claude Generated (Nov 2026, B2): periodic Coulomb (Ewald) of this partition. Real and reciprocal
    /// energies for the term report; coulomb_potentials (dE_pair/dq_i, size N) only when requested.
    double ewald_real = 0.0, ewald_recip = 0.0;
    double ewald_real_ms = 0.0, ewald_recip_ms = 0.0;  ///< P1 (Nov 2026): CPU time of the two Ewald parts in this partition
    Vector coulomb_potentials;
    /// N1 (Nov 2026): atomic virial of this partition; written only by the VIRIAL kernel instantiations
    FFVirialComponents virial;

    void reset(int natoms, bool do_gradient, bool store_components, bool want_potentials = false, bool want_virial = false) {
        energy.reset();
        timings.reset();
        ewald_real = ewald_recip = 0.0;
        ewald_real_ms = ewald_recip_ms = 0.0;
        if (want_virial)
            virial.reset();
        if (want_potentials) {
            if (coulomb_potentials.size() != natoms) coulomb_potentials.resize(natoms);
            coulomb_potentials.setZero();
        } else {
            coulomb_potentials.resize(0);
        }
        if (do_gradient) {
            if (gradient.rows() != natoms || gradient.cols() != 3)
                gradient.resize(natoms, 3);
            gradient.setZero();
            if (dEdcn.size() != natoms) dEdcn.resize(natoms);
            dEdcn.setZero();
            if (dEdcn_bond.size() != natoms) dEdcn_bond.resize(natoms);
            dEdcn_bond.setZero();
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

    /// Set atom types (element numbers)
    void setAtomTypes(const std::vector<int>& atoms);

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
    /// Claude Generated (Sep 2026, Known Issue #32): must be refreshed on EVERY evaluation,
    /// energy-only calls included - before, only setCNDerivatives() (gradient calls) wrote
    /// m_cn, so an energy-only call on a reused calculator evaluated the EN term with the CN
    /// of the last gradient geometry.
    void setCN(const Vector& cn) { m_cn = cn; }

    /// Set CN, CNF, and CN derivatives (gradient only)
    /// Claude Generated (WP4, May 2026): dcn now CNDerivStore (pair-list) instead of std::vector<SpMatrix>
    void setCNDerivatives(const Vector& cn, const Vector& cnf,
                          const CNDerivStore& dcn);

    /// Set the per-pair dC6/dCN source for the D4 dispersion CN gradient (Claude Generated,
    /// Oct 2026: per-pair arrays of the D4 generator instead of a dense N x N matrix)
    void setDC6DCNPtr(const D4DC6DCNPairs* ptr) { m_dc6dcn_ptr = ptr; }

    /// Fused D4 pair C6 (Claude Generated, Oct 2026): with a view set, calcD4Dispersion() computes
    /// each pair's C6 and dC6/dCN from per-atom data (d4_pair_c6.h) and ignores
    /// GFNFFDispersion::C6 and the setDC6DCNPtr() arrays' values (the pointer still gates the
    /// dC6/dCN term). nullptr = stored per-pair values. Set by GFNFF::prepareCNAndEEQ().
    void setD4FusedC6View(const D4FusedC6View* view) { m_d4_fused_c6 = view; }

    /// Homogeneous long-range D4 tail (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED A3c): energy and
    /// dE_tail/dCN computed by GFNFF (D4ParameterGenerator::homogeneousDispersionTail) before calculate().
    /// postProcess() stores the energy as FFEnergyComponents::disp_tail and adds dEdcn to the CN chain
    /// rule (gradient calls only, dEdcn of size natoms). clearDispersionTail(): no tail (no cell, switched
    /// off); every following calculate() is then bitwise what it was without this term.
    /// V-LRC (Nov 2026): e6/e8 are the C6 and C8 parts of `energy` (E6_tail + E8_tail = energy up to rounding);
    /// `pressure_lrc` true fills FFVirialComponents::tail_lrc = (3 e6 + 5 e8)/3 * I in postProcessVirial().
    void setDispersionTail(double energy, const Vector& dEdcn, double e6 = 0.0, double e8 = 0.0, bool pressure_lrc = false)
    {
        m_disp_tail_active = true;
        m_disp_tail_energy = energy;
        m_disp_tail_dEdcn = dEdcn;
        m_disp_tail_e6 = e6;
        m_disp_tail_e8 = e8;
        m_disp_tail_lrc = pressure_lrc;
        if (m_dispersion_scale != 1.0) {  // dispersion_scale diagnostic (Nov 2026): the 1.0 path above is the former code
            m_disp_tail_energy *= m_dispersion_scale;
            m_disp_tail_dEdcn *= m_dispersion_scale;
            m_disp_tail_e6 *= m_dispersion_scale;
            m_disp_tail_e8 *= m_dispersion_scale;
        }
    }
    void clearDispersionTail()
    {
        m_disp_tail_active = false;
        m_disp_tail_energy = 0.0;
        m_disp_tail_dEdcn.resize(0);
        m_disp_tail_e6 = m_disp_tail_e8 = 0.0;
        m_disp_tail_lrc = false;
    }
    /// V-LRC: the C6 and C8 parts of the last tail energy (0 without a tail), for tests and reports
    double dispersionTailE6() const { return m_disp_tail_active ? m_disp_tail_e6 : 0.0; }
    double dispersionTailE8() const { return m_disp_tail_active ? m_disp_tail_e8 : 0.0; }

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
    /// Claude Generated (Oct 2026): exact far-field shortcut of the Coulomb pair loop (see calcCoulomb).
    void setFarfieldShortcut(bool v) { m_farfield_shortcut = v; }
    /// Claude Generated (Oct 2026): TERM 1b chain-rule term set to 0 at CN_i <= 0 (gfnff.cn_zero_guard, see postProcess).
    void setCnZeroGuard(bool v) { m_cn_zero_guard = v; }
    /// Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED A1): orthorhombic cell in Bohr. Stored only;
    /// no energy term applies the minimum image yet (A2). Inactive by default = non-periodic path.
    void setPeriodicCell(const PBCUtils::OrthoCell& cell_bohr) { m_cell_bohr = cell_bohr; }
    const PBCUtils::OrthoCell& periodicCell() const { return m_cell_bohr; }
    bool hasPeriodicCell() const { return m_cell_bohr.active; }
    int numThreads() const { return m_num_threads; }

    // === Periodic Coulomb by Ewald summation (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED B2) ===
    /// Replace the GFN-FF Coulomb pair term by the Ewald sum of the Gaussian charges (ewald.h, B2):
    /// calcCoulombEwald() instead of calcCoulomb(); the self, background and Gaussian-background terms
    /// are added in postProcess(). Needs the periodic cell (setPeriodicCell) and the per-atom alp
    /// (setInteractionLists / setCoulombSelfEnergyParams); throws without a cell. The charges stay the
    /// ones last set by setEEQCharges(); the GFN-FF self-energy (TERM 2+3) is unchanged.
    void setEwald(const curcuma::ewald::EwaldParameters& p);
    void clearEwald() { m_ewald.reset(); m_coulomb_ewald = false; }
    /// P1 (Nov 2026): geometry data shared with the EEQ solver (built by whoever sees a geometry first) and the cap
    /// of its stored kernel arrays in bytes (above it the kernel is evaluated per call, as before P1).
    void setEwaldGeometryCache(std::shared_ptr<curcuma::ewald::EwaldGeometryCache> cache, size_t kernel_cap_bytes)
    {
        m_ewald_cache = std::move(cache);
        m_ewald_kernel_cap_bytes = kernel_cap_bytes;
    }
    const curcuma::ewald::EwaldGeometryCache* ewaldGeometryCache() const { return m_ewald_cache.get(); }
    bool ewaldActive() const { return m_coulomb_ewald && m_ewald != nullptr; }
    const curcuma::ewald::EwaldSum* ewald() const { return m_ewald.get(); }
    /// Also accumulate dE_pair/dq_i (periodic Coulomb pair energy only, for the periodic EEQ, B3).
    void setWantCoulombPotentials(bool v) { m_want_coulomb_potentials = v; }
    /// phi_i = dE_pair/dq_i after calculate() with setWantCoulombPotentials(true); empty otherwise.
    const Vector& coulombPotentials() const { return m_result_coulomb_potentials; }
    /// Periodic Coulomb pair energy (real + reciprocal + self + background + Gaussian background) of the last calculate().
    double ewaldPairEnergy() const { return m_ewald_terms.total(); }
    const curcuma::ewald::EwaldEnergyTerms& ewaldTerms() const { return m_ewald_terms; }
    /// Milliseconds of the per-step Ewald setup (pair list, kernel values, phase table) in the last calculate(); ~0 when
    /// the EEQ solver already built them for this geometry (P1 shared cache).
    double ewaldSetupMs() const { return m_ewald_setup_ms; }
    /// P1 (Nov 2026): CPU-sum milliseconds of the Ewald real-space and reciprocal parts in the last calculate().
    double ewaldRealMs() const { return m_ewald_real_ms; }
    double ewaldRecipMs() const { return m_ewald_recip_ms; }
    /*! \brief Claude Generated (Nov 2026, S2): the periodic Coulomb pair energy alone (real + reciprocal + self +
     * background + Gaussian background) and, on request, its gradient at fixed charges, with no other term.
     * The GPU hybrid (local terms on the device, EEQ and Ewald on the host) adds this to the device result.
     * Same code as calculate(): the shared geometry cache (P1), calcCoulombEwald() per partition over the pool,
     * reduced in partition order (reproducible at a fixed thread count), selfAndBackground() of postProcess().
     * The GFN-FF self term (TERM 2+3) and TERM 1b are NOT included (they stay on the device, k_coulomb_postprocess).
     * \param geometry_bohr N x 3 positions (Bohr); \param charges EEQ charges (e); \param gradient also dE/dr
     * \param gradient_out receives N x 3 Eh/Bohr when gradient is set (resized); \param terms_out the split (may be null)
     * \return E_pair (Eh). Throws without setEwald() or without per-atom alp. */
    double evaluateEwaldCoulomb(const Matrix& geometry_bohr, const Vector& charges, bool gradient,
                                Matrix* gradient_out, curcuma::ewald::EwaldEnergyTerms* terms_out);
    // === Atomic virial (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED N1) ===
    /// Request the virial tensor Xi_ab = -dE/d(eps_ab) (Eh) with every gradient evaluation. Off (default): every
    /// kernel runs its non-virial instantiation and nothing here costs anything. Meaningful only for GFN-FF with a
    /// periodic cell (the minimum-image frames define the atomic virial); without a cell the sum is the plain
    /// -sum_a x_a . G_a of the cluster, which is reported but not used.
    void setWantVirial(bool v) { m_want_virial = v; }
    bool wantVirial() const { return m_want_virial; }
    /// True after calculate(true) with the virial requested; evaluateEwaldCoulomb() (GPU hybrid) clears it.
    bool virialValid() const { return m_virial_valid; }
    /// Xi (3x3, Eh) of the last calculate(true) for the reported pressure: the model derivative plus the V-LRC
    /// crossing-pair correction (FFVirialComponents::pressureTotal; equal to modelVirial() with pressure_lrc false).
    Eigen::Matrix3d virial() const { return m_result_virial.pressureTotal(); }
    /// -dE/d(eps) of the model energy alone (FFVirialComponents::total), the quantity the FD gate checks
    Eigen::Matrix3d modelVirial() const { return m_result_virial.total(); }
    const FFVirialComponents& virialComponents() const { return m_result_virial; }

    /// Per-atom Coulomb self-energy inputs (setCoulombSelfEnergyParams order): tests copy them between instances.
    const Vector& coulombChiBase() const { return m_coul_chi_base; }
    const Vector& coulombGam() const { return m_coul_gam; }
    const Vector& coulombAlp() const { return m_coul_alp; }
    const Vector& coulombCnf() const { return m_coul_cnf; }
    const Vector& coulombChiStatic() const { return m_coul_chi_static; }

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
    /// Claude Generated (Nov 2026, A4a): current non-bonded repulsion pair list (read-only; tests)
    const std::vector<GFNFFRepulsion>& nonbondedRepulsions() const { return m_nonbonded_reps; }

    /**
     * @brief Mutable access to the D4 pair list, for the per-step C6 refresh only.
     *
     * Claude Generated (Sep 2026): GFNFF::refreshDispersionC6() rewrites each pair's C6 from
     * the current-step Gaussian weights. The pair set itself must not be changed through this
     * accessor (the partition ranges would go stale) — use updateD4Dispersions() for that.
     */
    std::vector<GFNFFDispersion>& d4DispersionsForC6Refresh() { return m_d4_dispersions; }

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

    /// Master bond list (read-only; used to re-derive the bond-HB cross reference)
    const std::vector<Bond>& bonds() const { return m_bonds; }
    /// Claude Generated (Nov 2026, A4c): current bonded-term lists (read-only; topology tests)
    const std::vector<Angle>& angles() const { return m_angles; }
    const std::vector<Dihedral>& dihedrals() const { return m_dihedrals; }
    const std::vector<Inversion>& inversions() const { return m_inversions; }

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

    /// Visit the installed D4 pair list so GFNFF can refresh the per-pair C6(CN) every call
    /// (Claude Generated, Sep 2026, Known Issue #32). Only C6 may be changed; the list itself
    /// stays fixed. Returns the number of lists visited (1 if a D4 pair list exists, else 0).
    int forEachD4PairList(const std::function<void(std::vector<GFNFFDispersion>&)>& f)
    {
        int n = 0;
        if (!m_d4_dispersions.empty()) {
            f(m_d4_dispersions);
            ++n;
        }
        return n;
    }
    int getHBondCount() const { return static_cast<int>(m_hbonds.size()); }

    /// Record the energy of every hydrogen-bond triple at each evaluation (CPU path only).
    /// Off by default: it costs one double per triple and a write per triple and step.
    /// Claude Generated (Oct 2026).
    void setRecordHBondEnergies(bool on);
    bool recordsHBondEnergies() const { return m_record_hb_energies; }

    /// Diagnostic factor on the hydrogen-bond term (PARAM gfnff.hb_scale, Claude Generated Nov 2026).
    /// Scales E_HB (all four cases and the per-case split), its gradient rows, its virial component and the
    /// recorded per-triple energies of calcHydrogenBonds() by s. Not scaled: the HB detection (same triples,
    /// unscaled pair-energy filter), the HB-CN modulation of the X-H bond strength (bond term), halogen bonds.
    /// 1.0 (default, the reference) runs the unchanged kernel; any other value goes through
    /// calcHydrogenBondsScaled() and deviates from the reference parametrisation.
    void setHBScale(double s) { m_hb_scale = s; }
    double hbScale() const { return m_hb_scale; }

    /// Diagnostic factors on three further GFN-FF terms (PARAMs gfnff.repulsion_scale, dispersion_scale,
    /// coulomb_scale; Claude Generated Nov 2026, "term-scale" switches, docs/REV_GFNFF_TODO.md #19). Same scheme as
    /// hb_scale: 1.0 (default, the reference) runs the unchanged kernels; any other value runs the unchanged kernel
    /// and multiplies its contribution to the partition accumulator by s (calcTermScaled).
    ///  - repulsion: bonded and non-bonded repulsion (energies, gradient, grad_repulsion, virial.repulsion).
    ///  - dispersion: D4 (and D3) pair term with its dE/dCN, plus the homogeneous tail (energy, dE/dCN, E6/E8 and so
    ///    virial.disp_tail and virial.tail_lrc, scaled in setDispersionTail). The ATM three-body term and BATM are not scaled.
    ///  - coulomb: the whole electrostatic component (pair or Ewald real/reciprocal, GFN-FF self term -q chi + q^2/2
    ///    (gam + ...), Ewald self/background/Gaussian background, TERM 1b chi(CN) chain rule, their virial parts,
    ///    the reported Ewald split and the pair potentials). The EEQ charges are NOT changed: s E_EEQ(q) has the same
    ///    constrained minimum, so the fixed-charge gradient times s is still the full gradient.
    /// Call before the first calculate() (the tail is scaled when GFNFF sets it).
    void setRepulsionScale(double s) { m_repulsion_scale = s; }
    void setDispersionScale(double s) { m_dispersion_scale = s; }
    void setCoulombScale(double s) { m_coulomb_scale = s; }
    double repulsionScale() const { return m_repulsion_scale; }
    double dispersionScale() const { return m_dispersion_scale; }
    double coulombScale() const { return m_coulomb_scale; }

    /// Energy of each triple of the list from the last evaluation, in list order. Empty unless
    /// recording was on before that evaluation and the list has not changed since.
    std::vector<GFNFFHBondEnergy> hydrogenBondEnergies() const;
    int getXBondCount() const { return static_cast<int>(m_xbonds.size()); }

private:
    int m_natoms = 0;
    int m_num_threads = 1;
    CxxThreadPool* m_pool = nullptr;

    // === Shared state (read-only per step) ===
    GeoGradMatrix m_geometry;  // WP-G: RowMajor for contiguous row(i) reads in inner loops
    std::vector<int> m_atom_types;
    Vector m_eeq_charges, m_topology_charges, m_d3_cn;
    Vector m_cn, m_cnf;
    CNDerivStore m_dcn;  // Claude Generated (WP4, May 2026): pair-list replaces std::vector<SpMatrix>
    const D4DC6DCNPairs* m_dc6dcn_ptr = nullptr;
    const D4FusedC6View* m_d4_fused_c6 = nullptr;  // Claude Generated (Oct 2026), see setD4FusedC6View()
    bool m_disp_tail_active = false;  // Claude Generated (Nov 2026, A3c), see setDispersionTail()
    double m_disp_tail_energy = 0.0;
    Vector m_disp_tail_dEdcn;
    double m_disp_tail_e6 = 0.0;    // V-LRC (Nov 2026): C6 and C8 parts of m_disp_tail_energy
    double m_disp_tail_e8 = 0.0;
    bool m_disp_tail_lrc = false;   // V-LRC: PARAM gfnff.pressure_lrc, fills virial.tail_lrc
    double m_e0 = 0.0;

    // Coulomb self-energy parameters (O(N), extracted at init)
    Vector m_coul_chi_base, m_coul_gam, m_coul_alp, m_coul_cnf, m_coul_chi_static;
    /// Claude Generated (Sep 2026): evaluate the N^2/2 Coulomb pairs on the fly from the per-atom
    /// data (charges + alpeeq) instead of reading a stored pair list. Set from the parameter set
    /// (GFN-FF only); see calcCoulomb().
    bool   m_coulomb_implicit = false;
    double m_coulomb_implicit_rcut = std::numeric_limits<double>::infinity();  // no cutoff, as the reference (Sep 2026)
    /// Claude Generated (Oct 2026): skip erf/exp for gamma_ij*r >= 6.5 (erf = 1 exactly in FP64), see calcCoulomb.
    bool   m_farfield_shortcut = true;
    /// Claude Generated (Oct 2026): opt-in, default false = reference formula (gfnff.cn_zero_guard).
    bool   m_cn_zero_guard = false;
    /// Claude Generated (Nov 2026, A1): periodic orthorhombic cell in Bohr (inactive = no PBC).
    PBCUtils::OrthoCell m_cell_bohr;
    /// Claude Generated (Nov 2026, B2): Ewald sum of the GFN-FF Coulomb term, see setEwald().
    std::unique_ptr<curcuma::ewald::EwaldSum> m_ewald;
    bool m_coulomb_ewald = false;
    bool m_want_coulomb_potentials = false;
    /// P1 (Nov 2026): pair list, real-space kernel values and phase table of the current geometry, shared with the
    /// periodic EEQ solver when GFNFF hands both the same object (setEwaldGeometryCache); else private.
    std::shared_ptr<curcuma::ewald::EwaldGeometryCache> m_ewald_cache;
    size_t m_ewald_kernel_cap_bytes = size_t(1024) << 20; ///< P1: memory cap of the stored kernel (ewald_kernel_cache_mb)
    curcuma::ewald::EwaldEnergyTerms m_ewald_terms;   ///< result split of the last calculate()
    Vector m_result_coulomb_potentials;
    double m_ewald_setup_ms = 0.0;
    double m_ewald_real_ms = 0.0, m_ewald_recip_ms = 0.0;  ///< P1 profile, CPU sums over partitions
    /// N1 (Nov 2026): atomic virial requested / result of the last calculate(true) / result is from a full calculate()
    bool m_want_virial = false;
    bool m_virial_valid = false;
    FFVirialComponents m_result_virial;

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
    bool m_record_hb_energies = false;
    double m_hb_scale = 1.0;                  ///< hb_scale diagnostic factor (setHBScale); 1.0 = reference, plugin ABI 5
    double m_repulsion_scale = 1.0;           ///< term-scale diagnostic factors (setRepulsionScale, ...); 1.0 = reference, plugin ABI 7
    double m_dispersion_scale = 1.0;
    double m_coulomb_scale = 1.0;
    std::vector<double> m_hb_triple_energy;   ///< per triple of m_hbonds, Eh; filled by calcHydrogenBonds when recording
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
    // N1 (Nov 2026): each kernel is a template on VIRIAL; the false instantiation is the former kernel statement for
    // statement (if constexpr blocks discarded), the true one also contracts frame positions with gradient rows.
    // The non-template name dispatches on m_want_virial outside the loops.
    void calcBonds(int p);
    void calcAngles(int p);
    void calcDihedrals(int p);
    void calcExtraTorsions(int p);
    void calcInversions(int p);
    void calcSTorsions(int p);
    void calcDispersion(int p);
    void calcD4Dispersion(int p);
    void calcBondedRepulsion(int p);
    void calcNonbondedRepulsion(int p);
    void calcCoulomb(int p);
    /// Claude Generated (Nov 2026, B2): periodic Coulomb, real-space pairs and k vectors of partition p.
    void calcCoulombEwald(int p);
    void calcHydrogenBonds(int p);
    void calcHalogenBonds(int p);
    void calcATM(int p);
    void calcATMGradient(int p);
    void calcBATM(int p);
    void computeHBCoordinationNumbers(int p);
    template <bool VIRIAL> void calcBondsImpl(int p);
    template <bool VIRIAL> void calcAnglesImpl(int p);
    template <bool VIRIAL> void calcDihedralsImpl(int p);
    template <bool VIRIAL> void calcExtraTorsionsImpl(int p);
    template <bool VIRIAL> void calcInversionsImpl(int p);
    template <bool VIRIAL> void calcSTorsionsImpl(int p);
    template <bool VIRIAL> void calcDispersionImpl(int p);
    template <bool VIRIAL> void calcD4DispersionImpl(int p);
    template <bool VIRIAL> void calcBondedRepulsionImpl(int p);
    template <bool VIRIAL> void calcNonbondedRepulsionImpl(int p);
    template <bool VIRIAL> void calcCoulombImpl(int p);
    template <bool VIRIAL> void calcCoulombEwaldImpl(int p);
    template <bool VIRIAL> void calcHydrogenBondsImpl(int p);
    void calcHydrogenBondsScaled(int p);  ///< hb_scale != 1: the unchanged kernel, then its contributions times s
    enum class ScaledTerm { Repulsion, Dispersion, Coulomb };
    /// term-scale != 1 (Nov 2026): runs `kernel` (the unchanged calc*Impl), then multiplies what it added to the
    /// accumulator of partition p (fields of `term`) by s
    void calcTermScaled(int p, ScaledTerm term, double s, const std::function<void()>& kernel);
    template <bool VIRIAL> void calcHalogenBondsImpl(int p);
    template <bool VIRIAL> void calcATMGradientImpl(int p);
    template <bool VIRIAL> void calcBATMImpl(int p);
    /// N1: CN chain rule, Ewald background and D4 tail parts of the virial (postProcess, gradient calls only)
    void postProcessVirial(bool has_term1b, const Vector& qtmp, bool tail_cn);

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
