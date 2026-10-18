/*
 * <Canonical closed-shell CCSD on top of the native RHF (WP13)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Coupled-cluster singles and doubles for a closed-shell RHF reference, spin-adapted, canonical
 * orbitals, exact (non-RI) MO integrals. References:
 *   - G. D. Purvis, R. J. Bartlett, J. Chem. Phys. 76, 1910 (1982)  (CCSD)
 *   - S. Hirata, P.-D. Fan, A. A. Auer, M. Nooijen, P. Piecuch, J. Chem. Phys. 120, 2581 (2004),
 *     Eqs. 35-45: the spin-adapted closed-shell amplitude equations used here, in the form of
 *     PySCF's pyscf/cc/rintermediates.py and rccsd.update_amps (chemists' notation integrals).
 *   - T. J. Lee, P. R. Taylor, Int. J. Quantum Chem. Symp. 23, 199 (1989)  (T1 diagnostic)
 *
 * The C++ code is a port of the validated Python prototype scripts/local_ccsd/ccsd_canonical.py
 * (docs/LOCAL_CCSD_PROTOTYPE_2026-10.md); every contraction there appears here with the same
 * index string, so the two can be read side by side.
 *
 * Unknowns: t1[i][a], t2[i][j][a][b] (i,j active occupied, a,b virtual), t2_ijab = t2_jiba.
 * Energy (closed shell, tau = t2 + t1 t1):
 *
 *   E_corr = 2 sum_ia f_ia t_ia + sum_ijab tau_ijab [2 (ia|jb) - (ib|ja)]
 *
 * Amplitude equations as residuals R = N - e t (e_ia = f_ii - f_aa, e_ijab = e_ia + e_jb; N is the
 * numerator of the textbook Jacobi update t = N / e). R = 0 at the solution. Solver: Jacobi step
 * t <- t + R / e with Pulay DIIS on (t1, t2) (P. Pulay, Chem. Phys. Lett. 73, 393 (1980);
 * G. E. Scuseria, T. J. Lee, H. F. Schaefer, Chem. Phys. Lett. 130, 236 (1986) for CC), start t1 = 0,
 * t2 = MP2 amplitudes, so the iteration-0 energy is the MP2 energy.
 *
 * Cost: memory O(v^4) for (ac|bd) (stored as one v^2 x v^2 matrix) plus the full MO tensor during
 * the transformation; time O(n^5) for the transformation and O(o^2 v^4) per iteration (the
 * "particle-particle ladder" GEMM). Practical up to roughly 150-200 basis functions.
 *
 * The parameters below extend the `qm` module (a third PARAM block of the module, as qm_mp2.h);
 * they are read by CCSDMethod (ccsd_method.cpp). -qm.frozen_core (qm_mp2.h) applies as for MP2.
 *
 * Claude Generated (Oct 2026, WP13 of docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md).
 * AI-generated, machine-tested; see the WP13 section of that file for what was compared.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>

// Claude Generated (Oct 2026, WP13): CCSD settings of the qm module.
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(cc_max_iter, Int, 100,
          "CCSD: maximum number of amplitude iterations. An unconverged CCSD is an error (no energy).",
          "Correlation", {})
    PARAM(cc_threshold, Double, 1.0e-10,
          "CCSD: convergence threshold on the change of the correlation energy between two iterations (Eh). Convergence also needs max|R| below -qm.cc_residual_threshold.",
          "Correlation", {})
    PARAM(cc_residual_threshold, Double, 1.0e-8,
          "CCSD: convergence threshold on the largest element of the amplitude residuals R1, R2 (Eh).",
          "Correlation", {})
    PARAM(cc_diis, Int, 8,
          "CCSD: number of amplitude vectors kept for Pulay DIIS (0 or 1 = plain Jacobi iteration without extrapolation).",
          "Correlation", {})
    PARAM(cc_max_memory_mb, Double, 32000.0,
          "CCSD: memory limit (MB) for the MO integrals (incl. the v^4 block), the transformation buffers and the amplitudes. Above it CCSD stops with an error that names the needed size. The AO tensor itself is governed by -qm.mp2_max_memory_mb.",
          "Correlation", {})
END_PARAMETER_DEFINITION

namespace qmint {

/// Settings of one CCSD evaluation (filled by CCSDMethod from the qm module).
struct CCSDOptions {
    bool frozen_core = false;
    int max_iter = 100;
    double e_threshold = 1.0e-10;   ///< |E(it) - E(it-1)|
    double r_threshold = 1.0e-8;    ///< max |R1|, |R2|
    int diis = 8;                   ///< DIIS vectors (<= 1: no DIIS)
    double max_memory_mb = 32000.0; ///< MO integrals + buffers + amplitudes
    double ao_max_memory_mb = 32000.0;  ///< AO tensor rebuild limit (-qm.mp2_max_memory_mb)
};

/// Result of one CCSD evaluation (Hartree, milliseconds).
struct CCSDResult {
    double e_corr = 0.0;        ///< CCSD correlation energy
    double e_mp2 = 0.0;         ///< iteration-0 energy (t1 = 0, t2 = MP2) = MP2 correlation energy
    double t1_diagnostic = 0.0; ///< sqrt(sum t_ia^2 / (2 o)) (Lee/Taylor; PySCF's convention)
    double max_residual = 0.0;  ///< max |R| of the last iteration
    double last_de = 0.0;       ///< energy change of the last iteration
    int iterations = 0;
    bool converged = false;
    int n_frozen = 0, n_occ_active = 0, n_virt = 0, n_basis = 0;
    double memory_mb = 0.0;     ///< estimated peak memory of the CCSD part (without the AO tensor)
    double t_integrals_ms = 0.0;  ///< AO tensor (rebuild or download)
    double t_transform_ms = 0.0;  ///< AO -> MO transformation and block extraction
    double t_iterations_ms = 0.0; ///< amplitude iterations
    std::string info;
};

/// A dense row-major tensor of rank 0-4 (rank 0: one number). Only a container: the
/// contractions are explicit calls of ccsdContract (qm_ccsd.cpp).
struct CCTensor {
    std::vector<int> dims;
    std::vector<double> data;
    CCTensor() = default;
    explicit CCTensor(std::vector<int> d) : dims(std::move(d))
    {
        size_t n = 1;
        for (int x : dims) n *= (size_t)x;
        data.assign(n, 0.0);
    }
    ~CCTensor();  // out of line (qm_ccsd.cpp): avoids -Winline noise
    CCTensor(const CCTensor&) = default;
    CCTensor(CCTensor&&) = default;
    CCTensor& operator=(const CCTensor&) = default;
    CCTensor& operator=(CCTensor&&) = default;
    size_t size() const { return data.size(); }
};

/// MO integral blocks over the active occupied (o) and virtual (v) orbitals, chemists' notation:
///   oooo[i][j][k][l] = (ij|kl), ovoo[i][a][j][k] = (ia|jk), ovov[i][a][j][b] = (ia|jb),
///   oovv[i][j][a][b] = (ij|ab), ovvo[i][a][b][j] = (ia|bj), ovvv[i][a][b][c] = (ia|bc),
///   vvvv: (v^2 x v^2) matrix V[(a,b)][(c,d)] = (ac|bd) (physicists' order, for the ladder GEMM).
/// foo, fvv, fov: Fock matrix blocks in the MO basis (canonical: diagonal, f_ov = 0).
struct CCSDIntegrals {
    int o = 0, v = 0;
    CCTensor oooo, ovoo, ovov, oovv, ovvo, ovvv;
    std::vector<double> vvvv;
    CCTensor foo, fvv, fov;
    CCSDIntegrals() = default;
    ~CCSDIntegrals();  // out of line (qm_ccsd.cpp)
};

/// @brief Builds the CCSD MO integral blocks from the active-basis AO tensor.
/// Four one-index transformations of the full tensor to the active MO space (m = n - n_frozen
/// orbitals), each one GEMM Out[q r s][p'] = sum_p In[p][q r s] C[p][p'] (threaded over rows),
/// then the blocks are copied out. Canonical Fock blocks from the orbital energies.
/// Claude Generated (Oct 2026, WP13).
void ccsdIntegrals(const double* eri, int n, const Matrix& C, const Vector& eps, int n_frozen, int n_occ,
                   int threads, CCSDIntegrals& I);

/// @brief Peak memory (MB) of ccsdIntegrals + ccsdSolve for the given sizes (estimate, used for the
/// -qm.cc_max_memory_mb check). Claude Generated (Oct 2026, WP13).
double ccsdMemoryMB(int n, int n_frozen, int n_occ, int diis);

/// @brief Closed-shell CCSD energy E = 2 f_ia t_ia + tau_ijab [2 (ia|jb) - (ib|ja)].
double ccsdEnergy(const CCTensor& t1, const CCTensor& t2, const CCSDIntegrals& I, int threads);

/// @brief CCSD residuals R1 = N1 - e_ia t1, R2 = N2 - e_ijab t2 (Hirata et al. 2004, Eqs. 35-45).
void ccsdResidual(const CCTensor& t1, const CCTensor& t2, const CCSDIntegrals& I, int threads, CCTensor& R1,
                  CCTensor& R2);

/// @brief Solves the CCSD amplitude equations (Jacobi + DIIS), fills e_corr, e_mp2, t1_diagnostic,
/// iterations, converged, max_residual, last_de, t_iterations_ms. verbosity >= 2 prints each iteration.
void ccsdSolve(const CCSDIntegrals& I, const CCSDOptions& opt, int threads, CCSDResult& res);

}  // namespace qmint
