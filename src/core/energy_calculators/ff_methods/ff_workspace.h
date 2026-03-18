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
#include "gfnff_parameters.h"
#include "forcefieldthread.h"  // For Bond, Angle, Dihedral, Inversion struct definitions

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
    double bonded_rep = 0, nonbonded_rep = 0;
    double coulomb = 0, hbond = 0, xbond = 0;
    double atm = 0, batm = 0, stors = 0;

    void reset() {
        bond = angle = dihedral = inversion = 0;
        dispersion = 0;
        bonded_rep = nonbonded_rep = 0;
        coulomb = hbond = xbond = 0;
        atm = batm = stors = 0;
    }

    FFEnergyComponents& operator+=(const FFEnergyComponents& o) {
        bond += o.bond; angle += o.angle; dihedral += o.dihedral;
        inversion += o.inversion; dispersion += o.dispersion;
        bonded_rep += o.bonded_rep; nonbonded_rep += o.nonbonded_rep;
        coulomb += o.coulomb; hbond += o.hbond; xbond += o.xbond;
        atm += o.atm; batm += o.batm; stors += o.stors;
        return *this;
    }

    double total() const {
        return bond + angle + dihedral + inversion + dispersion +
               bonded_rep + nonbonded_rep + coulomb + hbond + xbond +
               atm + batm + stors;
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
    Matrix gradient;          ///< N×3 gradient accumulator
    Vector dEdcn;             ///< N dE/dCN for chain-rule
    Vector dEdcn_bond;        ///< N bond-only dE/dCN for per-component attribution
    FFEnergyComponents energy;

    // Optional per-component gradients (only allocated if store_components=true)
    Matrix grad_bond, grad_angle, grad_torsion, grad_repulsion;
    Matrix grad_coulomb, grad_dispersion, grad_hb, grad_xb, grad_batm, grad_atm;

    bool has_components = false;

    void reset(int natoms, bool do_gradient, bool store_components) {
        energy.reset();
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
            auto resetMat = [&](Matrix& m) {
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
    std::pair<int,int> hbonds = {0,0};
    std::pair<int,int> xbonds = {0,0};
    std::pair<int,int> atm_triples = {0,0};
    std::pair<int,int> batm_triples = {0,0};
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

    /// Set CN, CNF, and sparse dcn derivatives (gradient only)
    void setCNDerivatives(const Vector& cn, const Vector& cnf,
                          const std::vector<SpMatrix>& dcn);

    /// Set dc6dcn pointer for D4 dispersion CN gradient
    void setDC6DCNPtr(const Matrix* ptr) { m_dc6dcn_ptr = ptr; }

    /// Set baseline energy (e0 from parameter set)
    void setE0(double e0) { m_e0 = e0; }

    // === Main calculation ===

    /// Calculate total energy (and gradient if requested)
    double calculate(bool gradient);

    // === Results ===

    const Matrix& gradient() const { return m_result_gradient; }
    const FFEnergyComponents& energyComponents() const { return m_result_energy; }
    const Vector& dEdcnTotal() const { return m_dEdcn_total; }
    const Vector& dEdcnBondTotal() const { return m_dEdcn_bond_total; }

    // Per-component gradient getters (only valid if store_components=true)
    const Matrix& gradientBond() const { return m_result_grad_bond; }
    const Matrix& gradientAngle() const { return m_result_grad_angle; }
    const Matrix& gradientTorsion() const { return m_result_grad_torsion; }
    const Matrix& gradientRepulsion() const { return m_result_grad_repulsion; }
    const Matrix& gradientCoulomb() const { return m_result_grad_coulomb; }
    const Matrix& gradientDispersion() const { return m_result_grad_dispersion; }
    const Matrix& gradientHB() const { return m_result_grad_hb; }
    const Matrix& gradientXB() const { return m_result_grad_xb; }
    const Matrix& gradientBATM() const { return m_result_grad_batm; }
    const Matrix& gradientATM() const { return m_result_grad_atm; }

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
    void updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds);
    void updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds);

    // Access master interaction list sizes (for diagnostics)
    int bondCount() const { return static_cast<int>(m_bonds.size()); }
    int dispersionPairCount() const { return static_cast<int>(m_dispersions.size() + m_d4_dispersions.size()); }
    int getHBondCount() const { return static_cast<int>(m_hbonds.size()); }
    int getXBondCount() const { return static_cast<int>(m_xbonds.size()); }

private:
    int m_natoms = 0;
    int m_num_threads = 1;
    CxxThreadPool* m_pool = nullptr;

    // === Shared state (read-only per step) ===
    Matrix m_geometry;
    std::vector<int> m_atom_types;
    Vector m_eeq_charges, m_topology_charges, m_d3_cn;
    Vector m_cn, m_cnf;
    std::vector<SpMatrix> m_dcn;
    const Matrix* m_dc6dcn_ptr = nullptr;
    double m_e0 = 0.0;

    // Coulomb self-energy parameters (O(N), extracted at init)
    Vector m_coul_chi_base, m_coul_gam, m_coul_alp, m_coul_cnf, m_coul_chi_static;

    // Term-enable flags
    bool m_dispersion_enabled = true;
    bool m_hbond_enabled = true;
    bool m_repulsion_enabled = true;
    bool m_coulomb_enabled = true;

    // === Master interaction lists (owned, moved from GFNFFParameterSet) ===
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

    // Cached bonded pairs for fast repulsion lookup
    std::set<std::pair<int,int>> m_bonded_pairs;

    // === Partitions + Accumulators ===
    std::vector<PartitionRanges> m_partitions;
    std::vector<FFAccumulator> m_accumulators;

    // === Result storage ===
    Matrix m_result_gradient;
    FFEnergyComponents m_result_energy;
    Vector m_dEdcn_total, m_dEdcn_bond_total;
    bool m_store_components = false;
    bool m_do_gradient = false;

    // Per-component result gradients
    Matrix m_result_grad_bond, m_result_grad_angle, m_result_grad_torsion;
    Matrix m_result_grad_repulsion, m_result_grad_coulomb, m_result_grad_dispersion;
    Matrix m_result_grad_hb, m_result_grad_xb, m_result_grad_batm, m_result_grad_atm;

    // === Core execution ===
    void executeGFNFF(int partition);
    void postProcess(bool gradient);
    void reduce();

    // === GFN-FF energy term calculators (ported from ForceFieldThread) ===
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
    void calcHydrogenBonds(int p);
    void calcHalogenBonds(int p);
    void calcATM(int p);
    void calcATMGradient(int p);
    void calcBATM(int p);
    void computeHBCoordinationNumbers(int p);

    // === Helpers ===

    /// Linear partition range for T threads
    static std::pair<int,int> linearRange(int total, int t, int T) {
        int begin = t * total / T;
        int end = (t + 1) * total / T;
        return {begin, end};
    }
};
