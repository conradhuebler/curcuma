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

#include <string>
#include <vector>

// Claude Generated (Oct 2026, WP11): MP2 settings of the qm module.
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(frozen_core, Bool, false,
          "MP2: keep the chemical core orbitals out of the correlation treatment (H-He none, Li-Ne 1s, Na-Ar 1s2s2p, K-Kr [Ne]3s3p; per atom, independent of the charge). The lowest n_core RHF orbitals are frozen.",
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
END_PARAMETER_DEFINITION

namespace qmint {

/// Settings of one MP2 evaluation (filled by MP2Method from the qm module).
struct MP2Options {
    bool frozen_core = false;
    double max_memory_mb = 32000.0;  ///< host n^4 tensor limit when the SCF stored none
    double batch_mb = 2000.0;        ///< first-quarter batch memory
    std::string gpu_mp2 = "auto";    ///< auto | on | off
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

}  // namespace qmint
