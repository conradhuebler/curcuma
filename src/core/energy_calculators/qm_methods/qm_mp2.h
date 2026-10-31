/*
 * <Canonical closed-shell MP2 on top of the native RHF (WP11)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Second-order Moller-Plesset perturbation theory for a closed-shell RHF reference
 * (C. Moller, M. S. Plesset, Phys. Rev. 46, 618 (1934); spin-adapted closed-shell form
 * e.g. Szabo/Ostlund, Modern Quantum Chemistry, eq. 6.74):
 *
 *   E_MP2 = sum_{ij in occ} sum_{ab in vir} (ia|jb) [2 (ia|jb) - (ib|ja)] / (e_i + e_j - e_a - e_b)
 *
 * split into the opposite-spin and same-spin pair contributions
 *
 *   E_OS = sum (ia|jb)^2 / D_ijab,   E_SS = sum (ia|jb) [(ia|jb) - (ib|ja)] / D_ijab,
 *
 * with the spin-component-scaled variants SCS-MP2 (S. Grimme, J. Chem. Phys. 118, 9095
 * (2003): 1.2 E_OS + 1/3 E_SS) and SOS-MP2 (Y. Jung et al., J. Chem. Phys. 121, 9793
 * (2004): 1.3 E_OS).
 *
 * The MO integrals (ia|jb) come from the AO tensor (mu nu|lam sig) by four one-index
 * transformations (O(n^5), "canonical" MP2, no resolution of the identity), done in
 * batches over the occupied index i:
 *   1. (i nu|lam sig) = sum_mu C_mu,i (mu nu|lam sig)          one GEMM per i batch
 *   2. (i a |lam sig) = sum_nu C_nu,a (i nu|lam sig)           per i
 *   3. (i a |j   sig) = sum_lam C_lam,j (i a|lam sig)          per (i,a)
 *   4. (i a |j   b  ) = sum_sig (i a|j sig) C_sig,b             per (i,a)
 * Both (ia|jb) and (ib|ja) belong to the same i, so the energy of an i batch needs no
 * other batch: extra memory O(batch * n^3 + threads * v n^2), never the o^2 v^2 tensor.
 *
 * The parameters below extend the `qm` module (a second PARAM block of the same module, as
 * `xtb` has two); they are read by MP2Method (mp2_method.cpp), not by the QMEngine ctor.
 *
 * Claude Generated (Oct 2026, WP11 of docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md).
 * AI-generated, machine-tested; see the WP11 section of that file for what was compared.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"
#include "src/core/parameter_macros.h"

#include "GTOIntegrals.hpp"  // GTO::Orbital (gradientTwoElectronMP2)

#include <memory>
#include <string>
#include <vector>

// Claude Generated (Oct 2026, WP11): MP2 settings of the qm module.
BEGIN_PARAMETER_DEFINITION(qm)
    // Claude Generated (Oct 31 2026): "auto" marker (as the scf module) so that double hybrids can default to true while
    // mp2/ri-mp2/scs-mp2/ccsd/lmp2 keep false; -export_run writes "auto" and a replay keeps the method default.
    PARAM(frozen_core, String, "auto",
          "Correlated methods: keep the chemical core orbitals out of the correlation treatment (H-He none, Li-Ne 1s, Na-Ar 1s2s2p, K-Kr [Ne]3s3p; per atom, independent of the charge). The lowest n_core orbitals of the reference are frozen. auto (default) = the method's own default: false for mp2, scs-mp2, sos-mp2, ri-mp2, scs-ri-mp2, sos-ri-mp2, ccsd, lmp2; true for the double hybrids (b2plyp, dsd-*, ..., since Oct 31 2026, as ORCA). true | false set it explicitly for every method.",
          "Correlation", {})
    PARAM(mp2_max_memory_mb, Double, 32000.0,
          "MP2: largest AO tensor (n^4 * 8 bytes, in MB) the MP2 step may build on the host when the SCF did not store one (integral-direct SCF). Above it MP2 stops with an error (no integral-direct MP2 yet).",
          "Correlation", {})
    PARAM(mp2_batch_mb, Double, 2000.0,
          "MP2: memory for one batch of first-quarter-transformed integrals (i nu|lam sig) (batch * n^3 * 8 bytes, in MB); at least one occupied orbital per batch.",
          "Correlation", {})
    PARAM(gpu_mp2, String, "auto",
          "With -gpu cuda|auto, MP2 transformation and energy sum on the device (cuBLAS DGEMM, FP64): auto | on (device whenever -gpu is active) | off (CPU). Reuses the stored device ERI tensor of the SCF, else uploads the host tensor. Any device error falls back to the CPU and is counted as qm_mp2.",
          "Correlation", {})
    // Claude Generated (Oct 2026): RI-MP2 (methods ri-mp2, scs-ri-mp2, sos-ri-mp2), implementation qm_rimp2.cpp
    PARAM(ri_mp2_aux, String, "def2-universal-JKFIT",
          "RI-MP2: auxiliary basis for the correlation fit (ia|jb) ~ sum_P B_ia^P B_jb^P: def2-universal-JKFIT (default) | def2-universal-JFIT (def2/JK, def2/J accepted) | path to a Turbomole-format file. No def2-SVP/C file ships with curcuma, so RI-MP2 here uses a JK/J fitting set as correlation fitting basis: a method variant, not identical to RI-MP2/def2-SVP/C of other programs. Independent of -qm.ri (the SCF stays exact unless -qm.ri j|jk is given).",
          "Correlation", {})
END_PARAMETER_DEFINITION

namespace qmint {

/// Settings of one MP2 evaluation (filled by MP2Method from the qm module).
struct MP2Options {
    bool frozen_core = false;
    double max_memory_mb = 32000.0;  ///< host n^4 tensor limit when the SCF stored none
    double batch_mb = 2000.0;        ///< first-quarter batch memory
    std::string gpu_mp2 = "auto";    ///< auto | on | off (canonical and RI path)
    bool ri = false;                 ///< RI-MP2 instead of canonical MP2 (QMEngine::computeRIMP2)
    std::string ri_aux = "def2-universal-JKFIT";  ///< RI-MP2 auxiliary basis (-qm.ri_mp2_aux)
};

/// Result of one MP2 evaluation (Hartree, milliseconds).
struct MP2Result {
    double e_os = 0.0;          ///< opposite-spin pair energy
    double e_ss = 0.0;          ///< same-spin pair energy
    double e_corr() const { return e_os + e_ss; }
    double e_scs() const { return 1.2 * e_os + e_ss / 3.0; }  ///< Grimme 2003
    double e_sos() const { return 1.3 * e_os; }               ///< Jung et al. 2004
    int n_frozen = 0, n_occ_active = 0, n_virt = 0, n_basis = 0;
    double t_integrals_ms = 0.0;  ///< building / downloading / uploading the AO tensor
    double t_transform_ms = 0.0;  ///< four index transformations
    double t_energy_ms = 0.0;     ///< pair-energy sum
    bool on_device = false;
    std::string info;             ///< where the tensor came from, batches, device stats
    // RI-MP2 only (Claude Generated, Oct 2026)
    bool ri = false;
    int n_aux = 0;                ///< fitted auxiliary functions (rows of B)
    std::string aux_name;         ///< e.g. def2-universal-JKFIT
    double t_ri_setup_ms = 0.0;   ///< 3-index integrals (mn|P) + metric (0 when the SCF's RI tensor is reused)
    std::string scf_ri;           ///< "" (exact SCF) or e.g. "RI-JK def2-universal-JKFIT"
    // Analytic gradient only (Claude Generated, Oct 2026; QMEngine::computeMP2Gradient, qm_mp2_grad.cpp)
    bool gradient = false;
    double t_grad_amplitudes_ms = 0.0;  ///< (ia|jb), amplitudes, P_ij/P_ab, (i nu|jb)-contraction
    double t_grad_backtrans_ms = 0.0;   ///< back-transformation to the AO 2PDM and the (kappa nu|lam sig) R contraction
    double t_grad_zvector_ms = 0.0;     ///< Lagrangian J/K build + Z-vector (CPHF) iterations
    double t_grad_density_ms = 0.0;     ///< relaxed density, energy-weighted density, 2PDM symmetrisation
    double t_grad_deriv2e_ms = 0.0;     ///< derivative ERIs contracted with the full 2PDM
    double t_grad_deriv1e_ms = 0.0;     ///< derivative S/T/V contracted with P and W
    int zvector_iterations = 0;
    double zvector_residual = 0.0;
};

/// Number of frozen core orbitals for an atom of nuclear charge Z (chemical core of the
/// preceding noble gas: 0 for H-He, 1 for Li-Ne, 5 for Na-Ar, 9 for K-Kr, 18 for Rb-Xe).
/// Claude Generated (Oct 2026, WP11).
int frozenCoreOrbitals(int Z);

/// @brief Canonical closed-shell MP2 pair energies from the AO tensor on the host.
/// @param eri    active-basis AO tensor, n^4 doubles, layout of ERITensor ((mu nu|lam sig) at
///               ((mu*n+nu)*n+lam)*n+sig)
/// @param C      MO coefficients (n x n, columns = MOs, ascending orbital energy)
/// @param eps    orbital energies (n)
/// @param n_frozen, n_occ   frozen and total doubly occupied orbitals (frozen ones are the lowest)
/// @param threads  OpenMP threads (BLAS for the batch GEMM, OpenMP over i otherwise)
/// @param batch_mb memory for one (i nu|lam sig) batch
/// Fills e_os, e_ss and the transformation / energy timings of `res`.
/// Claude Generated (Oct 2026, WP11).
void mp2EnergyHost(const double* eri, int n, const Matrix& C, const Vector& eps, int n_frozen, int n_occ,
                   int threads, double batch_mb, MP2Result& res);

/// @brief C (m x n) = A (m x k) * B(n x k)^T, all row-major and contiguous (lda = ldb = k, ldc = n).
/// Own register-blocked kernel (4 x 4 dot products, single thread): the build links the system BLAS
/// through Eigen, which on the development machine is single-threaded MKL (sequential layer,
/// ~68 GFLOPS per thread in a standalone Eigen test, see WP10 doc "BLAS (Korrektur)"); this kernel runs at ~29 GFLOPS
/// per thread and wins only through threading. Callers thread over independent calls.
/// Claude Generated (Oct 2026, RI-MP2).
void gemmABt(const double* A, int m, const double* B, int n, int k, double* C);

/// @brief RI-MP2 pair energies from the fitted 3-index tensor on the host.
/// @param B   fitted (ia|P) factors, layout [i][a][P] (o*v rows of naux), so that
///            (ia|jb) = sum_P B[i][a][P] B[j][b][P]
/// @param eo, ev  active occupied (o) and virtual (v) orbital energies
/// Same E_OS / E_SS expression as mp2EnergyHost; loops over pairs i >= j (the (j,i) pair has the
/// same energy), one v x v GEMM per pair, pair energies summed in fixed order afterwards (result
/// independent of the thread count). Fills e_os, e_ss, t_energy_ms.
/// Claude Generated (Oct 2026, RI-MP2).
void riMp2EnergyHost(const double* B, int o, int v, int naux, const double* eo, const double* ev, int threads,
                     MP2Result& res);

/// Intermediates of the closed-shell MP2 gradient that the host routine hands back (Claude Generated,
/// Oct 2026). o = active occupied, v = virtual, n = basis functions. With
///   t_ijab = (ia|jb) / (e_i + e_j - e_a - e_b),  T~_ijab = (c_os + c_ss) t_ijab - c_ss t_ijba
/// (MP2: c_os = c_ss = 1, T~ = 2t - t^ba):
struct MP2GradIntermediates {
    Matrix Poo;  ///< o x o:  P_ij = -2 sum_kab T~_ikab t_jkab   (unrelaxed occupied block)
    Matrix Pvv;  ///< v x v:  P_ab =  2 sum_ijc T~_ijac t_ijbc   (unrelaxed virtual block)
    Matrix Y1;   ///< n x o:  Y1_{kappa i} = sum_{nu lam sig} (kappa nu|lam sig) R_i^{nu lam sig}
    Matrix Y2;   ///< n x v:  Y2_{nu a}    = sum_{ijb} (i nu|jb) T~_ijab
};

/// @brief Canonical MP2 energy plus the gradient intermediates on the host (Claude Generated, Oct 2026).
/// @param eri    active-basis AO tensor (as mp2EnergyHost)
/// @param c_os, c_ss  spin-component scaling of the energy whose gradient is wanted (1, 1 for MP2)
/// @param Gamma  output, n^4 doubles: back-transformed amplitudes Gamma_{mu nu lam sig} =
///               sum_{ijab} T~_ijab C_mu,i C_nu,a C_lam,j C_sig,b (NOT symmetrised)
/// Fills res.e_os/e_ss, out, and res.t_grad_amplitudes_ms / t_grad_backtrans_ms. Results do not depend on
/// the thread count (fixed summation orders).
void mp2GradientHost(const double* eri, int n, const Matrix& C, const Vector& eps, int n_frozen, int n_occ,
                     double c_os, double c_ss, int threads, double batch_mb, double* Gamma, MP2GradIntermediates& out,
                     MP2Result& res);

/// @brief Gamma_cart = (Q x Q x Q x Q) Gamma_active for an n_act^4 tensor (Q: cartesian x active, block-diagonal
/// per shell); out gets ncart^4 doubles. Used for the device MP2 2e gradient. Claude Generated (Oct 2026).
void mp2ActiveToCartesian(const double* in, const Matrix& Q, std::unique_ptr<double[]>& out, int threads);

/// @brief Two-electron part 1/2 sum D (mn|ls)^A (Eh/Bohr) of the MP2 gradient, D = RHF part of P (as
/// gradientTwoElectron) + the separable P2 x P part + Dns. P, P2: CARTESIAN densities; Dns: non-separable
/// 2PDM in the ACTIVE basis (n_act^4, full 8-fold symmetry), taken to the cartesian basis with Q per shell
/// quartet (Q = nullptr or empty: active = cartesian). Implemented in qm_integrals.cpp (declared here so that
/// qm_integrals.hpp, included by the CUDA units, stays untouched). Claude Generated (Oct 2026).
Matrix gradientTwoElectronMP2(const std::vector<GTO::Orbital>& basis, const Matrix& P, const Matrix& P2,
                              const double* Dns, const Matrix* Q, int natoms, int threads = 1);

}  // namespace qmint
