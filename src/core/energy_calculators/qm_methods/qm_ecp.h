/*
 * <Native QM: semi-local effective core potentials (def2-ECP)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only (vs PySCF 2.14 / ORCA 6.1 on small molecules,
 * see docs/NATIVE_QM_IMPLEMENTATION.md "Elements beyond Ne and ECPs").
 *
 * Semi-local (Pople/Hay-Wadt/Stuttgart form) pseudopotential of a centre C with local channel L:
 *
 *   U_C(r) = U_L(r) + sum_{l=0}^{L-1} sum_m |lm> [U_l(r) - U_L(r)] <lm|,      U_x(r) = sum_k c_k r^(n_k-2) e^(-z_k r^2)
 *
 * (r measured from C, |lm> real spherical harmonics about C). Turbomole '$ecp' text as built in from the BSE
 * (qm_basis_data_def2_ecp.h): "ncore = N lmax = L", block "<L letter>" = U_L, blocks "<l>-<L>" = U_l - U_L, each line
 * "c n z". def2-ECP has only n = 2 (r^0), the code accepts only n = 2 and rejects others.
 *
 * Integrals over contracted cartesian Gaussians phi_a (QMEngine::m_gto_basis, coefficients carry the primitive and
 * contraction normalisation):
 *   local  <a|U_L|b>: three-centre Gaussian overlap, analytic (Gaussian product of a, b and e^(-z r_C^2)).
 *   semi-local <a|P_l dU_l P_l|b> = int_0^inf r^2 dU_l(r) sum_m F_a,lm(r) F_b,lm(r) dr with the angular projections
 *     F_a,lm(r) = int dOmega phi_a(C + r Omega) Y_lm(Omega). With D = A - C and the expansion
 *     e^(2 alpha D r cos g) = sum_lambda (2 lambda + 1) i_lambda(2 alpha D r) P_lambda(cos g):
 *       F_a,lm(r) = sum_{ijk} B_ijk r^(i+j+k) sum_lambda e^(-alpha (r - D)^2) [e^(-x) i_lambda(x)]_(x = 2 alpha D r) T_ijk,lambda,lm
 *       T_ijk,lambda,lm = (2 lambda + 1) int dOmega x^i y^j z^k P_lambda(Omega . D/|D|) Y_lm(Omega)  (unit sphere)
 *     B_ijk from the binomial expansion of (x - A_x)^a ... about C (McMurchie & Davidson, J. Comput. Phys. 44, 289
 *     (1981); Flores-Moreno et al., J. Comput. Chem. 27, 1009 (2006), half-numerical scheme). T is a polynomial
 *     integral, done EXACTLY with a Lebedev grid of sufficient degree; the radial integral is a Gauss-Chebyshev
 *     (2nd kind) quadrature on r = ln(2/(1-x))/ln 2 whose order (2^k - 1 points) is doubled until the matrix changes
 *     by less than max(`radial_tol` = 1e-13 Eh, 1e-14 max|V|); CURCUMA_QM_ECP_DEBUG=1 prints the change per order.
 *   gradient: derivative functions d phi/dA = 2 alpha phi(a+1) - a phi(a-1) through the same code; the ECP centre
 *     by translational invariance, dV/dC = -(dV/dA + dV/dB); the quadrature grid moves with C, so the gradient is
 *     the exact derivative of the quadrature energy at a fixed order.
 */

#pragma once

#include "src/core/global.h"
#include "GTOIntegrals.hpp"
#include "qm_gpu_ecp_backend.h"  // device seam (CUDA-free interface), Nov 2026

#include <string>
#include <vector>

namespace qmecp {

/// One term c * r^(n-2) * exp(-zeta r^2) of an ECP channel.
struct EcpTerm {
    double coef = 0.0;
    int n = 2;
    double zeta = 0.0;
};

/// ECP of one element: core electrons it replaces, local channel L (= "lmax" of the Turbomole text) and the
/// semi-local differences U_l - U_L for l = 0..L-1.
struct ElementECP {
    int ncore = 0;
    int lmax = 0;
    std::vector<EcpTerm> local;
    std::vector<std::vector<EcpTerm>> semilocal;  ///< [l], l < lmax
    bool empty() const { return ncore == 0 && local.empty() && semilocal.empty(); }
};

/// Parse one Turbomole '$ecp' element block (e.g. qmbasis::embeddedElementBlock("def2-ECP", Z)). False with `why`
/// on any format error or an r power other than n = 2.
bool parseTurbomoleECP(const std::string& block, ElementECP& out, std::string& why);

/// An atom of the molecule carrying an ECP (position in Bohr).
struct EcpCentre {
    int atom = -1;
    double x = 0, y = 0, z = 0;
    const ElementECP* ecp = nullptr;
};

struct EcpOptions {
    double radial_tol = 1.0e-13;  ///< max |dV| between two quadrature orders to accept the finer one (Eh); the
                                  ///< rounding floor floor_rel max|V| is accepted too
    double floor_rel = 1.0e-14;   ///< rounding floor relative to max|V| (Nov 2026: was 2e-15, inside the GEMM noise)
    int min_level = 7;            ///< first order 2^level - 1 radial points
    int max_level = 13;           ///< last order tried (8191 points); not converged there is an error
    int threads = 1;
};

/// Per-centre radial order the energy build accepted (reused by the gradient so energy and gradient share one
/// quadrature). Filled by buildECP.
struct EcpLevels {
    std::vector<int> level;  ///< per centre
};

/// V^ECP (ncart x ncart, cartesian basis of `basis`), local + semi-local, summed over all centres. False with
/// `why` if a quadrature did not converge.
bool buildECP(const std::vector<GTO::Orbital>& basis, const std::vector<EcpCentre>& centres, const EcpOptions& opt,
              Matrix& V, EcpLevels& levels, std::string& why);

/// dE/dR (natoms x 3, Eh/Bohr) of E = sum_ab P_ab V^ECP_ab for a symmetric cartesian density P (the density that
/// multiplies the one-electron Hamiltonian, e.g. the relaxed density of MP2). `levels` from buildECP of the
/// same geometry (empty: converge again).
Matrix gradientECP(const std::vector<GTO::Orbital>& basis, const std::vector<EcpCentre>& centres, const Matrix& P,
                   int natoms, const EcpOptions& opt, const EcpLevels& levels);

// ---------------------------------------------------------------------------------------------------------------
// Device path (Claude Generated, Nov 2026; seam qm_gpu_ecp_backend.h, kernels cuda/qm_gpu_ecp.cu). The host keeps
// the order loop, its acceptance test and the angular tables; the device object computes the local matrix, the
// semi-local matrix at one order and the gradient pieces. Same quantities as buildECP / gradientECP.
// ---------------------------------------------------------------------------------------------------------------

/// Host arrays behind a QMGpuEcpSetup: the functions of `basis` (and with `gradient` their centre derivatives).
struct EcpDeviceFunctions {
    std::vector<int> fn_atom, fn_term, term_p, term_slot, slot_atom;
    std::vector<double> fn_xyz, term_c, term_alpha, slot_alpha;
    int nbasis = 0, nfn = 0, natoms = 0;
    QMGpuEcpSetup view() const;
};
EcpDeviceFunctions deviceFunctions(const std::vector<GTO::Orbital>& basis, int natoms, bool gradient);

/// What the device path did (for the log line and the tests).
struct EcpDeviceStats {
    int centres = 0;
    int cpu_order_centres = 0;   ///< centres whose order test was too close to the threshold: CPU loop used
    double host_tables_ms = 0.0; ///< angular tables on the host
    double max_order_margin = 0.0;
};

/// buildECP with the heavy parts on `dev` (created from deviceFunctions(basis, natoms, false)). The order of each
/// centre is accepted with the CPU test on the device matrices; when the change lies within `order_margin` of the
/// threshold, that centre's semi-local loop runs on the CPU (so the order is the CPU's). False with `why` on a device
/// error (the caller falls back to buildECP) or a quadrature that did not converge.
bool buildECPDevice(const std::vector<GTO::Orbital>& basis, const std::vector<EcpCentre>& centres, const EcpOptions& opt,
                    QMGpuECP& dev, Matrix& V, EcpLevels& levels, EcpDeviceStats& stats, std::string& why);

/// gradientECP with the heavy parts on `dev` (created from deviceFunctions(basis, natoms, true)); the orders come
/// from `levels` (missing ones are converged on the CPU first, as gradientECP does). False on a device error.
bool gradientECPDevice(const std::vector<GTO::Orbital>& basis, const std::vector<EcpCentre>& centres, const Matrix& P,
                       int natoms, const EcpOptions& opt, const EcpLevels& levels, QMGpuECP& dev, Matrix& g,
                       EcpDeviceStats& stats, std::string& why);

}  // namespace qmecp
