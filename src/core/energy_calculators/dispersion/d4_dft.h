/*
 * <DFT-D4 for the native DFT methods (WP12 stage 4): functional parameters + the native D4 pipeline>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * D4 with EEQ charges and the ATM three-body term as dftd4 evaluates it for a density functional ("bj-eeq-atm",
 * E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)), built from the native D4 code of GFN2:
 *
 *   1. q   = EEQ(2019) charges + dq/dR          D4ChargeModel, Model::Eeq2019 (multicharge new_eeq2019_model)
 *   2. CN  = EN-weighted covalent CN, 30 Bohr     computeD4CovalentCN / addD4CovalentCNGradient (d4_ncoord)
 *   3. C6(CN, q) per reference state             D4ParameterGenerator, ChargeVariant::Eeq (clsq/clsh), covalent
 *                                                 CN on, ga 3, gc 2, wf 6 (the code native GFN2 uses)
 *   4. E2  = -sum_ij C6 (s6 t6 + s8 r4r2_ij t8),  D4Evaluator::computeEnergyAndGradient (per-reference charge)
 *      t_n = 1/(r^n + (a1 sqrt(r4r2_ij) + a2)^n)
 *   5. E3  = ATM, C6 at q = 0, s9, alpha 16,     D4Evaluator::computeATM, 40 Bohr (dftd4 disp3 default)
 *   6. gradient = direct + dE/dCN dCN/dR + dE/dq dq/dR (adjoint of the EEQ system)
 *
 * Only the parameters (s6, s8, a1, a2, s9) depend on the functional; d4DftParameters() is their table with the
 * source of each set. Native GFN2 (xtb_native.cpp calcDispersionEnergy) runs steps 2-6 itself with its SCF
 * Mulliken charges and the GFN2 reference charges; it does not call this class.
 *
 * Claude Generated (Oct 2026, WP12 stage 4, docs/DFT_ROADMAP/WP12_doppelhybride.md). AI-generated, machine-tested;
 * human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "d4_charge_model.h"
#include "src/core/global.h"

#include <memory>
#include <string>
#include <vector>

class D4ParameterGenerator;

namespace curcuma::dispersion {

class D4Evaluator;

/// One D4(EEQ)-ATM rational-damping parameter set (a2 in Bohr) with its source. Claude Generated (WP12 stage 4).
struct D4DftParams {
    const char* functional;  ///< curcuma functional key (lower case): pbe, b3lyp, tpss, scan, b2plyp, ...
    double s6, s8, a1, a2, s9;
    const char* source;
};

/// The table (d4_dft.cpp). Claude Generated (Oct 2026, WP12 stage 4).
const std::vector<D4DftParams>& d4DftParameters();
/// Row by functional key; nullptr if no D4 parameters are tabulated. Claude Generated (WP12 stage 4).
const D4DftParams* findD4DftParameters(const std::string& functional);

/**
 * @brief Energy and gradient of DFT-D4 (EEQ charges, ATM) for one parameter set.
 *
 * Reusable across geometries and molecules (the generator rebuilds its C6 references when a new element shows
 * up). Thread safety: one instance per thread.
 */
class D4DftDispersion {
public:
    explicit D4DftDispersion(const D4DftParams& params);
    ~D4DftDispersion();

    /**
     * @brief E_D4 (Eh) and, if gradient, dE_D4/dR (Eh/Bohr, N x 3) at geom_bohr.
     * @return false with error() set for elements without EEQ(2019) parameters (Z > 103) or a failed setup.
     */
    bool compute(const std::vector<int>& atoms, const Matrix& geom_bohr, double total_charge, bool gradient,
                 double& energy, Matrix& gradient_bohr);

    const D4DftParams& params() const { return m_params; }
    const std::string& error() const { return m_error; }
    /// Diagnostics of the last compute(): EEQ charges, D4 CN, two-body and ATM energies.
    const Vector& charges() const { return m_q; }
    const std::vector<double>& coordinationNumbers() const { return m_cn; }
    double twoBodyEnergy() const { return m_e2; }
    double threeBodyEnergy() const { return m_e3; }

    /// dftd4 real-space cutoffs (cutoff.f90 defaults; EEQ CN 25 Bohr is fixed in D4ChargeModel).
    static constexpr double kCnCutoffBohr = 30.0;
    static constexpr double kAtmCutoffBohr = 40.0;

private:
    D4DftParams m_params;
    std::unique_ptr<D4ParameterGenerator> m_generator;
    std::unique_ptr<D4Evaluator> m_evaluator;
    D4ChargeModel m_eeq;
    Vector m_q;
    std::vector<double> m_cn;
    double m_e2 = 0.0, m_e3 = 0.0;
    std::string m_error;
};

}  // namespace curcuma::dispersion
