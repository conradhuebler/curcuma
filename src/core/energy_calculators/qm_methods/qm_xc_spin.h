/*
 * <Native QM: spin-polarised exchange-correlation for open-shell Kohn-Sham (UKS), stage 5b>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Unrestricted Kohn-Sham (J. A. Pople, P. M. W. Gill, B. G. Johnson, Chem. Phys. Lett. 199, 557 (1992)) needs the
 * functional of the two spin densities rho_a, rho_b and, for a GGA, of the three invariants
 *     sigma_aa = |grad rho_a|^2,  sigma_ab = grad rho_a . grad rho_b,  sigma_bb = |grad rho_b|^2
 * (libxc's spin-polarised variables). Per spin s the Kohn-Sham matrix is
 *     F_s = h + J[P_a + P_b] - a_x K[P_s] + V_xc,s,
 *     V_xc,s,mn = sum_g w_g [ v_rho_s phi_m phi_n + (2 v_sigma_ss grad rho_s + v_sigma_ab grad rho_s') . grad(phi_m phi_n) ].
 *
 * Exchange: the spin-scaling relation (G. L. Oliver, J. P. Perdew, Phys. Rev. A 20, 397 (1979))
 *     E_x[rho_a, rho_b] = 1/2 E_x[2 rho_a] + 1/2 E_x[2 rho_b]   (sigma -> 4 sigma_ss),
 * so every closed-shell exchange kernel of qm_xc_kernels.h (Slater, PBE-X, B88) is reused unchanged; per spin
 *     e = rho_s eps(2 rho_s, 4 sigma_ss),  v_rho_s = v_rho(2 rho_s, 4 sigma_ss),  v_sigma_ss = 2 v_sigma(...).
 * Correlation has no such relation; the spin-polarised forms (formulas in qm_xc_spin.cpp, checked point by point
 * against libxc 7.0.0 with spin = 1 in ctest qm_uks_kernel) are written once with 5-variable dual numbers:
 *   - VWN5   (libxc LDA_C_VWN):     eps = e_P + alpha f(z)(1 - z^4)/f''(0) + (e_F - e_P) f(z) z^4  (VWN 1980 eq. 4.4)
 *   - VWN-RPA(libxc LDA_C_VWN_RPA): eps = e_P (1 - f(z)) + e_F f(z), RPA fits (libxc's VWN_RPA interpolation)
 *   - PW92   (libxc LDA_C_PW_MOD):  Perdew-Wang 1992 eq. 8 with the three G fits
 *   - PBE-C  (libxc GGA_C_PBE):     PBE 1996 eqs. 3, 7, 8 with phi(z) = ((1+z)^(2/3) + (1-z)^(2/3))/2
 *   - LYP    (libxc GGA_C_LYP):     Miehlich, Savin, Stoll, Preuss, Chem. Phys. Lett. 157, 200 (1989), eq. 2
 * Not spin-polarised (open shell rejected with a message): the meta-GGAs TPSS/SCAN/B95 and P86 (double hybrids are
 * closed-shell only anyway).
 *
 * No AO screening (every block uses every function; -qm.xc_ao_threshold is ignored by the open-shell path). The point
 * kernels live in qm_xc_spin_kernels.h (host and device since stage 6a; the CUDA path is QMGpuXC::vxcSpin/gradientSpin).
 * Claude Generated (Nov 2026, stage 5b). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "qm_xc.h"

#include <string>
#include <vector>

namespace qmxc {

/// True when every term of f has a spin-polarised form here; else false and `why` names the first kernel without one.
bool supportsSpinPolarised(const XCFunctional& f, std::string& why);

/**
 * @brief Spin-polarised functional at one grid point (all terms of f, weights applied).
 *
 * In: rho_a, rho_b >= 0, sigma_aa, sigma_ab, sigma_bb (ignored unless gga). The point is skipped (false, all outputs
 * zero) when rho_a + rho_b < rho_min; an exchange spin channel is skipped when 2 rho_s < rho_min.
 * Out: ex, ec = exchange and correlation ENERGY DENSITIES (per volume, rho eps); v[0..4] = d e / d (rho_a, rho_b,
 * sigma_aa, sigma_ab, sigma_bb). Claude Generated (Nov 2026).
 */
bool evaluatePointSpin(const XCFunctional& f, bool gga, double rho_min, double ra, double rb, double saa, double sab,
                       double sbb, double& ex, double& ec, double v[5]);

/// Result of one spin-polarised V_xc build.
struct XCResultSpin {
    Matrix Va, Vb;          ///< V_xc,alpha and V_xc,beta (active AO basis, symmetric)
    double exc = 0.0, ex = 0.0, ec = 0.0;
    double nelec = 0.0;     ///< sum_g w_g (rho_a + rho_b)
    long points_used = 0;
    double ms = 0.0;
    // out of line: the implicit inline ones trip -Winline (as XCFunctional)
    XCResultSpin();
    XCResultSpin(const XCResultSpin&);
    XCResultSpin& operator=(const XCResultSpin&);
    ~XCResultSpin();
};

/// V_xc,a, V_xc,b and E_xc of the spin densities Pa, Pb on the grid (blocks of `block` points, fixed chunking as
/// buildVxc: the result does not depend on `threads`). Claude Generated (Nov 2026).
void buildVxcSpin(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& Pa, const Matrix& Pb,
                  const XCFunctional& f, double rho_threshold, int threads, XCResultSpin& out, int block = 128);

/**
 * @brief Basis-function and point parts of dE_xc/dR for the spin densities Pa, Pb (as buildXCGradient per spin):
 *   dE_xc/dA_k = -2 sum_s sum_g sum_{m on A} [ d_k phi_m (P_s Z_s)_m + sum_j u^s_j d_k d_j phi_m (P_s phi)_m ],
 *   Z_s = w v_rho_s phi + sum_j u^s_j d_j phi,  u^s_j = w (2 v_sigma_ss d_j rho_s + v_sigma_ab d_j rho_s').
 * edens = e_x + e_c per point for the Becke weight derivative. Claude Generated (Nov 2026).
 */
void buildXCGradientSpin(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& Pa,
                         const Matrix& Pb, const XCFunctional& f, double rho_threshold, int threads,
                         const std::vector<int>& bf_atom, int natoms, XCGradient& out, int block = 128);

}  // namespace qmxc
