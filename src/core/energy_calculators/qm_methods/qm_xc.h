/*
 * <Native QM: exchange-correlation functionals of closed-shell Kohn-Sham DFT (WP5: LDA, WP6: PBE, WP7: B3LYP,
 *  WP11: meta-GGA TPSS and SCAN)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Kohn-Sham DFT (W. Kohn, L. J. Sham, Phys. Rev. 140, A1133 (1965)) replaces the Hartree-Fock exchange by an
 * exchange-correlation energy E_xc[rho] = integral rho(r) eps_xc(rho(r)) d^3r. For a closed-shell density
 * (rho_alpha = rho_beta = rho/2) the Kohn-Sham matrix is
 *
 *     F = h + J[P] + V_xc - (a_x / 2) K[P],      V_xc,mn = integral v_xc(r) phi_m(r) phi_n(r) d^3r,
 *     v_xc = d(rho eps_xc) / d rho   (LDA; a GGA adds the sigma term of buildVxc),
 *     E    = Tr(P h) + 1/2 Tr(P J) - (a_x / 4) Tr(P K) + E_xc + E_nn,
 *
 * with a_x = 0 for a pure density functional (LDA, WP5; GGA, WP6) and a_x > 0 for a hybrid (B3LYP, WP7). The
 * integrals over r are sums over the molecular grid of qm_grid.h (weights w_g, Becke partition).
 *
 * WP5 (this file): local density approximation, spin-restricted (zeta = 0):
 *   - exchange: Dirac/Slater, eps_x = -(3/4) (3/pi)^(1/3) rho^(1/3), v_x = (4/3) eps_x
 *     (P. A. M. Dirac, Proc. Cambridge Philos. Soc. 26, 376 (1930); J. C. Slater, Phys. Rev. 81, 385 (1951)).
 *   - correlation: Vosko-Wilk-Nusair Pade fit of the paramagnetic correlation energy (S. H. Vosko, L. Wilk,
 *     M. Nusair, Can. J. Phys. 58, 1200 (1980)), x = sqrt(r_s), r_s = (3 / (4 pi rho))^(1/3):
 *       eps_c = A [ ln(x^2/X(x)) + 2b/Q atan(Q/(2x+b))
 *                   - b x0/X(x0) ( ln((x-x0)^2/X(x)) + 2(b+2 x0)/Q atan(Q/(2x+b)) ) ],
 *       X(x) = x^2 + b x + c,  Q = sqrt(4c - b^2),  v_c = eps_c - (x/6) d eps_c/dx.
 *     Two parameter sets, both A = 0.0310907 Eh (= (1 - ln 2)/pi^2, the high-density limit eps_c -> A ln r_s):
 *       vwn5    : fit to the Ceperley-Alder Monte-Carlo data (VWN's "functional V"),
 *                 x0 = -0.10498, b = 3.72744, c = 12.9352     (libxc LDA_C_VWN; ORCA "VWN5", LSD, B3LYP)
 *       vwn_rpa : fit to the RPA correlation energy (VWN's "functional III" in the closed-shell case),
 *                 x0 = -0.409286, b = 13.0720, c = 42.7198    (libxc LDA_C_VWN_RPA; PySCF "VWN3"; the correlation of
 *                 Gaussian's and PySCF's/libxc's B3LYP, ORCA "B3LYP/G" / "VWN3")
 *     For zeta = 0 only the paramagnetic fit enters; the spin interpolation of the variants is not used.
 *     NOTE the name trap: libxc's LDA_C_VWN_1 ... LDA_C_VWN_4 all equal VWN5 for a closed shell; only
 *     LDA_C_VWN_RPA (id 8) is the RPA fit. PySCF "LDA" alone is Slater exchange WITHOUT correlation.
 *
 * WP6: generalised gradient approximation PBE (J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865
 * (1996); equation numbers below refer to that paper), closed shell (zeta = 0, phi = 1), as libxc GGA_X_PBE +
 * GGA_C_PBE (= PySCF/ORCA "PBE"). A GGA energy density depends on rho and sigma = |grad rho|^2 (total density):
 *   - exchange (eqs. 10, 14): eps_x = eps_x^LDA(rho) F_x(s), F_x = 1 + kappa - kappa / (1 + mu s^2 / kappa),
 *       s = |grad rho| / (2 k_F rho), k_F = (3 pi^2 rho)^(1/3), so s^2 = sigma / (4 (3 pi^2)^(2/3) rho^(8/3)),
 *       kappa = 0.804, mu = beta pi^2 / 3 = 0.2195149727645171.
 *   - correlation (eqs. 3, 7, 8): eps_c = eps_c^PW92(r_s) + H(r_s, t),
 *       H = gamma ln(1 + (beta/gamma) t^2 (1 + A t^2) / (1 + A t^2 + A^2 t^4)),
 *       A = (beta/gamma) / (exp(-eps_c^PW92 / gamma) - 1),  t = |grad rho| / (2 k_s rho),  k_s = sqrt(4 k_F / pi),
 *       so t^2 = sigma pi / (16 k_F rho^2); beta = 0.06672455060314922, gamma = (1 - ln 2) / pi^2.
 *     eps_c^PW92 (J. P. Perdew, Y. Wang, Phys. Rev. B 45, 13244 (1992), eq. 10, paramagnetic set), NOT VWN:
 *       G(r_s) = -2 A0 (1 + alpha1 r_s) ln(1 + 1 / (2 A0 (beta1 r_s^(1/2) + beta2 r_s + beta3 r_s^(3/2) + beta4 r_s^2))),
 *       with libxc's LDA_C_PW_MOD constants A0 = 0.0310907 (the paper rounds to 0.031091), alpha1 = 0.21370,
 *       beta1..4 = 7.5957, 3.5876, 1.6382, 0.49294 (libxc's GGA_C_PBE uses exactly these).
 *   v_rho = d(rho eps)/d rho and v_sigma = d(rho eps)/d sigma are analytic (derivations in qm_xc.cpp).
 *   sigma enters only as sigma (never sqrt(sigma)), so sigma = 0 needs no special case. Numerical thresholds of
 *   libxc (rho < 1e-15 exchange / 1e-12 correlation dropped, sigma raised to 1e-40 / 1e-32) are NOT copied into the
 *   kernels; the SCF drops rho < -qm.xc_density_threshold (default 1e-12) for all terms, GGA kernels in addition
 *   rho < 1e-20 (the formulas give NaN near 1e-27) (WP6 doc).
 *
 * WP7: hybrid B3LYP (P. J. Stephens, F. J. Devlin, C. F. Chabalowski, M. J. Frisch, J. Phys. Chem. 98, 11623
 * (1994); three-parameter form of A. D. Becke, J. Chem. Phys. 98, 5648 (1993)):
 *     E_xc = 0.08 E_x^Slater + 0.72 E_x^B88 + 0.20 E_x^HF + 0.19 E_c^VWN + 0.81 E_c^LYP
 *   (a0 = 0.20, ax = 0.72, ac = 0.81; Slater weight 1 - a0 - ax, VWN weight 1 - ac), libxc HYB_GGA_XC_B3LYP (VWN-RPA)
 *   and HYB_GGA_XC_B3LYP5 (VWN5). Which VWN is "B3LYP" depends on the program (WP5 doc table): Gaussian, libxc and
 *   PySCF "B3LYP" use the RPA fit, ORCA "B3LYP" and PySCF "B3LYP5" use VWN5, ORCA "B3LYP/G" the RPA fit.
 *   - B88 exchange (A. D. Becke, Phys. Rev. A 38, 3098 (1988), eq. 8), per spin s with x_s = |grad rho_s| / rho_s^(4/3):
 *       e_x = sum_s [ -C_x rho_s^(4/3) - beta rho_s^(4/3) x_s^2 / (1 + 6 beta x_s asinh x_s) ],
 *       C_x = (3/2) (3/(4 pi))^(1/3), beta = 0.0042; closed shell rho_s = rho/2, x_s^2 = 2^(2/3) sigma / rho^(8/3).
 *   - LYP correlation (C. Lee, W. Yang, R. G. Parr, Phys. Rev. B 37, 785 (1988)) in the form without the Laplacian
 *     of rho obtained by partial integration (B. Miehlich, A. Savin, H. Stoll, H. Preuss, Chem. Phys. Lett. 157, 200
 *     (1989), eq. 2; the form libxc's GGA_C_LYP evaluates). For rho_a = rho_b = rho/2, |grad rho_s|^2 = sigma/4 it
 *     reduces to (derivation in qm_xc.cpp)
 *       e_c = -a rho / (1 + d r) - a b exp(-c r) / (1 + d r) [ C_F rho - sigma r^5 (1/24 + 7 delta / 72) ],
 *       r = rho^(-1/3), delta = c r + d r / (1 + d r), C_F = (3/10) (3 pi^2)^(2/3),
 *       a = 0.04918, b = 0.132, c = 0.2533, d = 0.349.
 *   The exact-exchange part is the K build of the SCF: F gets -(a_x/2) K, E gets -(a_x/4) Tr(PK) (a_x = 0.2).
 *   Further global hybrids (e.g. PBE0 = 0.75 PBE-X + PBE-C + 0.25 HF) are one more table function like
 *   b3lypFunctional(); not added (not validated).
 *
 * WP11 stage 1 (Oct 2026): meta-GGAs, rung 3, which also depend on the kinetic energy density
 *   tau = 1/2 sum_i^occ |grad psi_i|^2 (both spins) = 1/2 sum_mn P_mn grad phi_m . grad phi_n
 * (the convention of libxc's spin-unpolarised tau and of PySCF's eval_rho):
 *   - TPSS (J. Tao, J. P. Perdew, V. N. Staroverov, G. E. Scuseria, Phys. Rev. Lett. 91, 146401 (2003)), libxc
 *     MGGA_X_TPSS + MGGA_C_TPSS (= PySCF 'TPSS', ORCA 'TPSS');
 *   - SCAN (J. Sun, A. Ruzsinszky, J. P. Perdew, Phys. Rev. Lett. 115, 036402 (2015)), libxc MGGA_X_SCAN +
 *     MGGA_C_SCAN (= PySCF 'SCAN').
 *   Formulas, equation numbers, constants, libxc thresholds and the derivative scheme (dual numbers): the meta-GGA
 *   block of qm_xc_kernels.h. V_xc gains (1/2) v_tau grad phi_m . grad phi_n (buildVxc). Energies and (stage 2)
 *   the analytic gradient (buildXCGradient); stage 3: the same quantities on CUDA (cuda/qm_gpu_xc.cu, same point
 *   kernels; docs/DFT_ROADMAP/WP11_meta_gga.md).
 *
 * Claude Generated (Oct 2026, WP5; WP6 GGA; WP7 B3LYP; WP11 meta-GGA). AI-generated, machine-tested only
 * (docs/DFT_ROADMAP/WP5_lda_vwn5.md, WP6_pbe_gga.md, WP7_b3lyp_hybrid.md, WP11_meta_gga.md).
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (row-major), Vector

#include "qm_grid.h"
#include "qm_xc_kernels.h"  // point kernels shared with cuda/qm_gpu_xc.cu (one source for CPU and GPU)

#include <string>
#include <vector>

namespace qmxc {

/// Parameter set of the VWN correlation (see the file header).
enum class VWNVariant {
    VWN5,   ///< Monte-Carlo fit (VWN functional V); libxc LDA_C_VWN, ORCA VWN5
    VWN_RPA ///< RPA fit (VWN functional III for zeta = 0); libxc LDA_C_VWN_RPA, PySCF VWN3, ORCA VWN3
};

/// Parse "vwn5" | "vwn_rpa" (aliases "vwn3", "rpa", "vwn-rpa"); false for anything else. Claude Generated (WP5).
bool parseVWNVariant(const std::string& name, VWNVariant& v);
/// "VWN5" | "VWN-RPA" (output label).
std::string vwnName(VWNVariant v);

// XCKernel and XCTerm are defined in qm_xc_kernels.h (shared with the GPU kernels, Oct 2026).

/// True for an exchange kernel (its energy is reported as E_x), false for correlation.
bool isExchangeKernel(XCKernel k);
/// True when the kernel depends on sigma = |grad rho|^2 (GGA).
bool isGGAKernel(XCKernel k);
/// True when the kernel depends on tau (meta-GGA, WP11); it needs sigma as well.
bool isMGGAKernel(XCKernel k);
/// Short kernel name for the output ("Slater", "VWN5", "PBE-X", ...).
std::string kernelName(XCKernel k);

/**
 * @brief Closed-shell exchange-correlation functional: what the SCF needs to know.
 *
 * E_xc = sum_terms weight * integral rho eps_kernel(rho, sigma) plus a_x (exact_exchange) times the HF exchange.
 * WP5: lda (Slater + VWN), WP6: pbe (PBE_X + PBE_C), WP7: b3lyp (hybrid, exact_exchange = 0.2).
 */
struct XCFunctional {
    std::string name = "none";        ///< method name (lda, pbe, ...)
    std::string description;          ///< output label, e.g. "SVWN5 (Slater exchange + VWN5 correlation)"
    double exact_exchange = 1.0;      ///< a_x: fraction of HF exchange in F and E (1 = HF, 0 = pure DFT, 0.2 = B3LYP)
    std::vector<XCTerm> terms;        ///< empty = no density functional (HF)
    /// True when there is any density-functional part (else the SCF is plain HF).
    bool hasXC() const { return !terms.empty(); }
    /// GGA: the kernels need sigma = |grad rho|^2 (AO derivatives on the grid) and return v_sigma (also true for a
    /// meta-GGA).
    bool needsGradient() const;
    /// Meta-GGA (WP11): the kernels also need tau and return v_tau; V_xc and E_xc then come from the meta-GGA branch
    /// of buildVxc (CPU only, no analytic gradient, no device kernel yet).
    bool needsTau() const;
    /// Human-readable label for the output.
    std::string label() const;
    // Special members out of line: the implicit inline ones trip -Winline (as qmgrid::GridOptions).
    XCFunctional();
    XCFunctional(const XCFunctional&);
    XCFunctional& operator=(const XCFunctional&);
    ~XCFunctional();
};

/// The functional of `-method lda` with the given VWN variant (Slater + VWN, a_x = 0). Claude Generated (WP5).
XCFunctional ldaFunctional(VWNVariant v);
/// The functional of `-method pbe`: PBE exchange + PBE correlation (PW92 LDA part), a_x = 0. Claude Generated (WP6).
XCFunctional pbeFunctional();
/// The functional of `-method b3lyp`: 0.08 Slater + 0.72 B88 + 0.19 VWN + 0.81 LYP, a_x = 0.20. The VWN variant
/// (VWN5: ORCA B3LYP, PySCF B3LYP5; VWN-RPA: Gaussian/libxc/PySCF B3LYP, ORCA B3LYP/G) is named in the label.
/// Claude Generated (WP7).
XCFunctional b3lypFunctional(VWNVariant v);
/**
 * @brief The meta-GGA functionals by method name (WP11): "tpss" (TPSS_X + TPSS_C), "scan" (SCAN_X + SCAN_C), a_x = 0.
 *
 * One table (qm_xc.cpp, kMetaGGATable): name, label, exact-exchange fraction, terms. A global meta-GGA hybrid is one
 * more row once its kernels exist (TPSSh: 0.9 TPSS_X + TPSS_C with a_x = 0.1; r2SCAN needs its own kernels); neither
 * is implemented or validated. False for an unknown name. Claude Generated (WP11).
 */
bool metaGGAFunctional(const std::string& name, XCFunctional& f);

// ---------------------------------------------------------------------------------------------------------------
// Point kernels (closed shell, rho = total density > 0). eps = energy per electron (Eh), v = d(rho eps)/d rho,
// vsigma = d(rho eps)/d sigma with sigma = |grad rho|^2 (total density, as libxc's spin-unpolarised sigma).
// ---------------------------------------------------------------------------------------------------------------

/// Dirac/Slater exchange: eps_x = -(3/4) (3 rho / pi)^(1/3), v_x = -(3 rho / pi)^(1/3). Claude Generated (WP5).
void slaterExchange(double rho, double& eps, double& v);
/// VWN paramagnetic correlation (variant per file header). Claude Generated (WP5).
void vwnCorrelation(double rho, VWNVariant variant, double& eps, double& v);
/// Perdew-Wang 1992 paramagnetic correlation (libxc LDA_C_PW_MOD constants). Claude Generated (WP6).
void pw92Correlation(double rho, double& eps, double& v);
/// PBE exchange (PBE 1996 eqs. 10, 14). Claude Generated (WP6).
void pbeExchange(double rho, double sigma, double& eps, double& vrho, double& vsigma);
/// PBE correlation (PBE 1996 eqs. 3, 7, 8 with phi = 1). Claude Generated (WP6).
void pbeCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma);
/// Becke 1988 exchange (Becke 1988 eq. 8, closed shell). Claude Generated (WP7).
void b88Exchange(double rho, double sigma, double& eps, double& vrho, double& vsigma);
/// Lee-Yang-Parr correlation, closed shell, Miehlich et al. 1989 eq. 2 (no Laplacian). Claude Generated (WP7).
void lypCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma);
/// Any kernel by enum (sigma ignored by LDA kernels, vsigma = 0 for them). Claude Generated (WP6).
void evaluateKernel(XCKernel k, double rho, double sigma, double& eps, double& vrho, double& vsigma);
/// A meta-GGA kernel (TPSS_X, TPSS_C, SCAN_X, SCAN_C) at (rho, sigma, tau): eps, v_rho, v_sigma, v_tau, with libxc's
/// tau floor and sigma <= 8 rho tau (qm_xc_kernels.h). NaN for a kernel that is not a meta-GGA. Claude Generated (WP11).
void evaluateKernelMGGA(XCKernel k, double rho, double sigma, double tau, double& eps, double& vrho, double& vsigma,
                        double& vtau);

/**
 * @brief Functional values on a block of points.
 *
 * rho (npts), sigma (npts, only read when f.needsGradient(); may be empty for an LDA); points with
 * rho < rho_threshold get exc = vrho = vsigma = 0 (they are left out of E_xc and V_xc).
 * Outputs per point: ex, ec (energy per electron of the exchange and correlation terms, weights applied),
 * vrho = sum weight * v, vsigma = sum weight * vsigma (zero for an LDA). Claude Generated (WP5, sigma WP6).
 */
void evaluateXC(const XCFunctional& f, const Vector& rho, const Vector& sigma, double rho_threshold, Vector& ex,
                Vector& ec, Vector& vrho, Vector& vsigma);
/// evaluateXC for a meta-GGA functional (f.needsTau()): also tau in and v_tau out. Claude Generated (WP11).
void evaluateXCMeta(const XCFunctional& f, const Vector& rho, const Vector& sigma, const Vector& tau,
                    double rho_threshold, Vector& ex, Vector& ec, Vector& vrho, Vector& vsigma, Vector& vtau);

/// Result of one V_xc build.
struct XCResult {
    Matrix V;               ///< V_xc in the active AO basis (nbf x nbf, symmetric)
    double exc = 0.0;       ///< E_xc = sum_g w_g rho_g (eps_x + eps_c)
    double ex = 0.0;        ///< exchange part of E_xc
    double ec = 0.0;        ///< correlation part of E_xc
    double nelec = 0.0;     ///< sum_g w_g rho_g (electron count on the grid, diagnostic)
    long points_used = 0;   ///< grid points with rho >= threshold
    double ms = 0.0;        ///< wall time of the build
};

/**
 * @brief V_xc and E_xc of the density P on the molecular grid.
 *
 *   LDA: V_mn = sum_g w_g v_rho phi_m phi_n,
 *   GGA: V_mn = sum_g w_g [ v_rho phi_m phi_n + 2 v_sigma grad rho . grad(phi_m phi_n) ],
 *        (dE_xc/dP_mn with sigma = grad rho . grad rho and d grad rho / dP_mn = grad(phi_m phi_n)),
 *        built as X = phi^T [ (w v_rho / 2) phi + 2 w v_sigma (grad rho . grad phi) ], V = X + X^T,
 *   E_xc = sum_g w_g rho_g eps_xc(rho_g, sigma_g),
 *   meta-GGA (WP11, f.needsTau()): V_mn = sum_g w_g [ v_rho phi_m phi_n + 2 v_sigma grad rho . grad(phi_m phi_n)
 *        + (1/2) v_tau grad phi_m . grad phi_n ]   (d tau / d P_mn = (1/2) grad phi_m . grad phi_n), the tau part
 *        built as sum_k (d_k phi)^T diag(w v_tau / 2) (d_k phi); own function (buildVxcMeta in qm_xc.cpp), screened
 *        and unscreened, chunking as the screened path; the LDA/GGA code paths are not touched,
 *
 * in blocks of `block` points (AOEvaluator values, rho via qmgrid::densityOnBlock). The blocks are split into a
 * fixed number of contiguous chunks that depends only on nbf (not on `threads`); every chunk accumulates its own
 * matrix and the chunks are added in order, so the result does not depend on the thread count.
 * AO screening (Oct 2026): with an active `screen` (qmgrid::BlockScreen, -qm.xc_ao_threshold > 0) the blocks are
 * those of the screen (`block` ignored) and every block works on the kept functions S only: phi_S, P_SS, the small
 * matrix X_SS = phi_S^T Y_S is scattered into the chunk matrix over the union of the chunk's functions (64 chunks, halved above 1 GB of chunk matrices; independent of `threads`). Without a screen
 * (or threshold 0) the code path is the unscreened one, bitwise as before.
 * Claude Generated (WP5; screening Oct 2026).
 */
void buildVxc(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
              const XCFunctional& f, double rho_threshold, int threads, XCResult& out, int block = 128,
              const qmgrid::BlockScreen* screen = nullptr);

/// Nuclear gradient of E_xc on the molecular grid (WP8), Eh/Bohr, natoms x 3 each.
struct XCGradient {
    Matrix basis;          ///< basis-function (AO centre) part: the gradient without grid response
    Matrix points;         ///< grid response, part 1: the points of atom B move with B (sum_p w_p grad_r e(r_p) on B)
    Vector edens;          ///< energy density e_p = rho_p eps_xc(r_p) per point (for qmgrid::gridWeightGradient)
    long points_used = 0;  ///< points with rho >= threshold
    double ms = 0.0;       ///< wall time
};

/**
 * @brief Basis-function and point parts of dE_xc/dR_A for the closed-shell density P (active AO basis).
 *
 * With phi_m on centre A, d phi_m / d A_k = - d_k phi_m, and the variational KS density (only explicit
 * derivatives; the orbital response is the -W dS term of the one-electron gradient):
 *   LDA: dE_xc/dA_k = -2 sum_g w_g v_rho sum_{m on A} d_k phi_m (P phi)_m,
 *   GGA: dE_xc/dA_k = -2 sum_g w_g sum_{m on A} [ d_k phi_m (P Z)_m + sum_j v_j d_k d_j phi_m (P phi)_m ],
 *        Z = v_rho phi + sum_j v_j d_j phi,  v_j = 2 v_sigma d_j rho
 * (Johnson, Gill, Pople, J. Chem. Phys. 98, 5612 (1993), eqs. 9-13 for the closed shell; second AO derivatives
 * from AOEvaluator deriv 2). `points`: for a point of atom B the translation identity d_r e = - sum_C (d e / d R_C)
 * moves - sum of all basis parts of that point onto B. bf_atom[m] = atom of active function m.
 * Fixed chunking from nbf as buildVxc: independent of `threads`. With an active `screen` the blocks and kept
 * functions of the screen (as buildVxc).
 *   meta-GGA (WP11 stage 2, f.needsTau()): with d tau / d A_k = - sum_{m on A} sum_j d_k d_j phi_m (P d_j phi)_m the
 *   bracket gains sum_j (v_tau / 2) d_k d_j phi_m (P d_j phi)_m (still AO deriv 2); own function (buildXCGradientMeta
 *   in qm_xc.cpp), the LDA/GGA paths are unchanged. `edens` is rho eps_xc(rho, sigma, tau), so the Becke weight
 *   derivative (qmgrid::gridWeightGradient) and the point part need no meta-GGA term of their own.
 * Claude Generated (Oct 2026, WP8; screening Oct 2026; meta-GGA WP11 stage 2).
 */
void buildXCGradient(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                     const XCFunctional& f, double rho_threshold, int threads, const std::vector<int>& bf_atom,
                     int natoms, XCGradient& out, int block = 128, const qmgrid::BlockScreen* screen = nullptr);

}  // namespace qmxc
