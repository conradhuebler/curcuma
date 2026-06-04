/*
 * D4 single-shot EEQ charge model — charges + analytical dq/dx response
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License v3 (or later).
 *
 * Claude Generated 2026 — AP ∂q/∂x (Phase 2).
 *
 * Canonical single-shot electronegativity-equilibration (EEQ) charge model for
 * the D4 dispersion zeta scaling, in the spirit of Caldeweyher et al. 2019
 * (dftd4). Unlike the GFN-FF two-phase solver (topological Phase 1 + environment
 * corrections), this is ONE smooth linear system, so its derivative w.r.t.
 * nuclear coordinates has a clean closed form — required for the D4 q-response
 * gradient term ∂E_D4/∂q · ∂q/∂x.
 *
 * Linear system (augmented with charge conservation):
 *     A_ii = γ_i + sqrt(2/π)/sqrt(α_i)
 *     A_ij = erf(γ_ij·r_ij)/r_ij,   γ_ij = 1/sqrt(α_i + α_j)
 *     b_i  = -χ_i + κ_i·sqrt(CN_i)
 *     [ A  1 ] [ q ]   [ b ]
 *     [ 1  0 ] [ λ ] = [ Q ]
 * Element parameters (χ, γ, α, κ) are the angewChem2020 EEQ set (gfnff_par.h);
 * α is stored squared, matching EEQSolver::getParameters. CN is the GFN-FF
 * log-compressed erf coordination number (CNCalculator::calculateGFNFFCN), the
 * same definition used for the D4 C6 weighting.
 *
 * Reference: E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019).
 */

#pragma once

#include "src/core/global.h"

#include <Eigen/Dense>
#include <vector>

namespace curcuma::dispersion {

class D4ChargeModel {
public:
    D4ChargeModel() = default;

    /**
     * @brief Solve the single-shot EEQ system and cache state for the gradient.
     *
     * @param atoms        Atomic numbers (1-based), size N
     * @param geom_bohr    Coordinates in Bohr, N×3
     * @param total_charge Total molecular charge Q
     * @return Atomic partial charges q (size N)
     *
     * Caches the LU factorisation of the augmented matrix plus the per-atom CN,
     * raw CN and log-compression factors so addChargeResponseGradient() can run
     * without re-solving or recomputing the CN.
     */
    Vector computeCharges(const std::vector<int>& atoms,
                          const Matrix& geom_bohr,
                          double total_charge);

    /**
     * @brief Accumulate the D4 charge-response gradient into grad_out.
     *
     * grad_out(m,:) += Σ_A dEdq(A) · ∂q_A/∂R_m   (Hartree/Bohr).
     *
     * Uses the adjoint (Z-vector) trick: solve M·z = [dEdq; 0] once, then
     * contract z against the closed-form geometry derivatives of A and b.
     * Must be called after computeCharges() for the same geometry.
     *
     * @param dEdq     Per-atom dE_D4/dq (size N), from D4Evaluator
     * @param grad_out N×3 Cartesian gradient accumulator (Hartree/Bohr)
     */
    void addChargeResponseGradient(const Vector& dEdq, Matrix& grad_out) const;

    const Vector& charges() const { return m_q; }
    bool valid() const { return m_n > 0; }

    /**
     * @brief Resolve the per-atom EEQ parameters (angewChem2020 set).
     *
     * Single source of truth for the element-table lookups used by the EEQ
     * system: χ_i, γ_i, α_i² (stored squared, matching getParameters), κ_i (CN
     * scaling) and the 4/3-scaled covalent radius in Bohr. computeCharges() uses
     * it; the GPU port (Stage 5) calls it to feed the device solver the identical
     * parameters. Out-of-range Z (not 1..86) → χ=1, γ=0, α²=1, κ=0, rcov=0.
     * Each output vector is resized to atoms.size(). Claude Generated 2026.
     */
    static void resolveParams(const std::vector<int>& atoms,
                              std::vector<double>& chi,
                              std::vector<double>& gam,
                              std::vector<double>& alpha_sq,
                              std::vector<double>& cnf,
                              std::vector<double>& rcov_bohr);

private:
    // Cached state from the last computeCharges() call.
    int m_n = 0;
    std::vector<int> m_atoms;
    Matrix m_geom;                       ///< N×3, Bohr
    Vector m_q;                          ///< charges (N)
    Eigen::PartialPivLU<Matrix> m_lu;    ///< factorisation of augmented M (N+1)

    // Per-atom EEQ parameters (resolved once per geometry).
    std::vector<double> m_alp;           ///< α_i (squared, as in getParameters)
    std::vector<double> m_cnf;           ///< κ_i (CN scaling)

    // CN data needed for the b-term derivative.
    std::vector<double> m_cn;            ///< log-compressed CN_i
    std::vector<double> m_cn_raw;        ///< raw erf-sum CN_i (pre log compression)
    std::vector<double> m_rcov_bohr;     ///< scaled covalent radii (4/3 · rcov · Å→Bohr)
};

}  // namespace curcuma::dispersion
