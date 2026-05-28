/*
 * D4 Dispersion Evaluator — energy + analytical gradient
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License v3 (or later).
 *
 * Claude Generated 2026 — extracted from GFN-FF's ForceFieldThread D4 path,
 * generalised to serve both GFN-FF and native GFN2.
 *
 * Pedagogical note on damping formulas:
 *   The two damping variants "modified BJ" (GFN-FF Fortran) and standard
 *   Becke-Johnson D4 (Caldeweyher 2019) are mathematically the same expression
 *
 *        E_pair = -ζc6 · C6 · ( s6·t6  +  s8·r4r2_ij·t8 )
 *
 *   with t_n = 1/(r^n + R0^n) and R0 = a1·sqrt(r4r2_ij) + a2.
 *   They differ only in the scaling values: GFN-FF hard-codes s6=1, s8=2;
 *   GFN2 uses s6=1, s8=2.7 (Bannwarth 2019). The enum DampingFormula keeps
 *   the call site declarative and gives a hook for future variants that may
 *   need a structurally different kernel.
 *
 * Reference: gfnff_gdisp0.f90:365-377 (Fortran), Caldeweyher et al. J. Chem.
 *            Phys. 150, 154122 (2019), Bannwarth et al. JCTC 15, 1652 (2019).
 */

#pragma once

#include "d4param_generator.h"
#include "src/core/energy_calculators/ff_methods/gfnff_parameters.h"  // GFNFFDispersion
#include "src/core/global.h"

#include <Eigen/Dense>
#include <vector>

namespace curcuma::dispersion {

enum class DampingFormula {
    ModifiedBJ_GFNFF,  ///< GFN-FF default (s6=1, s8=2)
    StandardBJ_D4      ///< Standard D4 BJ damping (Caldeweyher 2019)
};

struct D4Params {
    double s6 = 0.0;
    double s8 = 0.0;
    double s9 = 0.0;     ///< Three-body ATM scaling (s9; GFN2=5.0). 0 disables ATM.
    double a1 = 0.0;
    double a2 = 0.0;     ///< In Bohr
    double alpha = 16.0; ///< ATM zero-damping exponent (dftd4 alp; GFN2=16.0)
    DampingFormula damping = DampingFormula::ModifiedBJ_GFNFF;
    /// AP6b: native GFN2 uses the exact dftd4 per-reference charge weighting
    /// (D4ParameterGenerator::weightedC6Gfn2). GFN-FF leaves this false and keeps
    /// its CN-only + single-zeta-prefactor pair data (pair.C6·pair.zetac6).
    bool per_reference_charge = false;
};

/**
 * @brief Stateless D4 dispersion energy + gradient kernel.
 *
 * The evaluator borrows a non-owning pointer to a D4ParameterGenerator that
 * carries the reference data (C6 tables, r4r2 polarizabilities, dc6/dCN).
 * It provides two entry points:
 *
 *   1. pairEnergyAndGradient(): per-pair primitive used by GFN-FF's
 *      ForceFieldThread, which already maintains a precomputed pair list.
 *
 *   2. computeEnergyAndGradient(): whole-molecule loop used by native GFN2,
 *      which has no precomputed pair structures.
 *
 * Thread safety: the evaluator itself is stateless; concurrent calls are
 * safe as long as the backing D4ParameterGenerator's read-only data is
 * stable for the geometry under evaluation.
 *
 * CUDA hook: launchGpuKernel() is a documented extension point — empty in
 * this AP. GFN-FF's existing GPU pipeline does not go through this class;
 * future GFN2-CUDA work can plug in a damping-aware kernel here.
 */
class D4Evaluator {
public:
    D4Evaluator(D4ParameterGenerator* data, const D4Params& params);
    virtual ~D4Evaluator() = default;

    // ---------- per-pair (ForceFieldThread style) ----------
    //
    // Inputs:
    //   pair        precomputed C6, r4r2ij, r0_squared, zetac6, r_cut
    //   rij_bohr    vector r_i - r_j in Bohr
    //
    // Outputs:
    //   E             pair energy (Hartree)
    //   dE_dr_vec     dE/dr_i in Hartree/Bohr (negate for atom j)
    //   dE_dCN_i,j    chain-rule accumulators -> caller distributes via dCN/dx
    //   disp_sum_out  s6·t6 + s8·r4r2ij·t8 (so the caller can reuse it for
    //                 the CN chain rule without recomputing)
    //
    // Returns false if the pair is beyond cutoff or degenerate (E=0, grad=0).
    bool pairEnergyAndGradient(const GFNFFDispersion& pair,
                               const Eigen::Vector3d& rij_bohr,
                               bool with_gradient,
                               double& E,
                               Eigen::Vector3d& dE_dr_vec,
                               double& dE_dCN_i,
                               double& dE_dCN_j,
                               double& disp_sum_out) const;

    // ---------- whole-molecule (GFN2 style) ----------
    //
    // Loops i<j over atoms, sums pair energies, accumulates gradients +
    // dE/dCN vector. The caller multiplies dEdCN_out by dCN/dx_k (computed
    // elsewhere) to obtain the CN-chain-rule contribution to the Cartesian
    // gradient.
    //
    // dEdq_out holds the per-atom charge derivative dE_D4/dq_A. It is the
    // first half of the charge-response chain rule; the caller multiplies it
    // by dq_A/dx (from EEQ response or CPSCF) to obtain the Cartesian
    // contribution. Only the zetac6 factor is charge-dependent (the CN-only
    // C6 weighting is not), so this term is exact for the GFN-FF/GFN2 model.
    // Pass with_dEdq=false to skip it (e.g. static-prefactor mode).
    //
    // Side effect: refreshes the underlying generator's CN-dependent state
    // (dc6dcn matrix) via D4ParameterGenerator::updateCNValuesForGradient()
    // if with_gradient is true.
    //
    // Returns total D4 dispersion energy (Hartree).
    double computeEnergyAndGradient(const std::vector<int>& atoms,
                                    const Matrix& geometry_bohr,
                                    bool with_gradient,
                                    Matrix& gradient_out,
                                    Vector& dEdCN_out,
                                    Vector& dEdq_out,
                                    bool with_dEdq = false);

    // ---------- ATM three-body term (GFN2) ----------
    //
    // Axilrod-Teller-Muto triple-dipole dispersion with Chai--Head-Gordon zero
    // damping + BJ critical radii. Exact port of dftd4 get_atm_dispersion
    // (damping/atm.f90). Returns the ATM energy (Hartree) and, if with_gradient,
    // ACCUMULATES the Cartesian gradient (Hartree/Bohr) into gradient_out and the
    // CN chain-rule term into dEdCN_out (both are += , not overwritten — call
    // after computeEnergyAndGradient so the 2-body contributions stay).
    //
    // dftd4 evaluates the ATM at the q=0 (CN-only) reference C6, so the term is
    // charge-INDEPENDENT — there is no q-response (matches tblite
    // get_dispersion_nonsc). No-op when m_params.s9 == 0.
    //
    // cutoff_bohr: real-space cutoff for the 3-body sum (tblite GFN2 uses 25.0).
    double computeATM(const std::vector<int>& atoms,
                      const Matrix& geometry_bohr,
                      bool with_gradient,
                      Matrix& gradient_out,
                      Vector& dEdCN_out,
                      double cutoff_bohr = 25.0);

    // ---------- CUDA hook (not used in this AP) ----------
    //
    // Place-holder for a future GFN2-CUDA path. GFN-FF's existing kernel
    // (gfnff_kernels.cu::k_dispersion) does NOT route through here — it
    // consumes the precomputed pair data directly. Adding a GFN2 GPU path
    // later will likely add a new kernel variant selected by m_params.damping.
    virtual void launchGpuKernel(/* TODO: device buffers */) const {}

    const D4Params& params() const { return m_params; }

private:
    // The unified single-pair formula (see header comment).
    // Returns disp_sum = s6·t6 + s8·r4r2ij·t8 so callers can reuse it.
    inline double evalDispSum(double t6, double t8, double r4r2ij) const;

    D4ParameterGenerator* m_data;  // non-owning
    D4Params m_params;
};

}  // namespace curcuma::dispersion
