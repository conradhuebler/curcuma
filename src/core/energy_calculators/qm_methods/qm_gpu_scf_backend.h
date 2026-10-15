/*
 * < Device seam of the native QM engine: resident SCF loop (Q-G4) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, Q-G4 of docs/QM_GPU_ROADMAP.md). Separate from qm_gpu_backend.h on
// purpose: the heavy ERI kernel TUs (cuda/qm_gpu_eri.cu, cuda/qm_gpu_direct.cu) include that header
// and take ~9 min each to compile; this one is only seen by cuda/qm_gpu_scf.cu and the core.

#include "qm_gpu_backend.h"

#include <string>

namespace qmint {

/// Q-G4 (Oct 2026): device-resident closed-shell RHF SCF state. The SCF driver itself stays in
/// QMEngine::runSCF (qm_scf.cpp); this object only owns the device matrices and performs the
/// steps of one iteration on them, so P, F, the DIIS history and the orbitals never leave the
/// device inside the loop. Per iteration the host receives O(1) scalars (max|e|, ||dP||, E)
/// plus the O(m^2) DIIS/ADIIS subspace matrices, and sends back the m coefficients.
/// Host matrices are n x n row-major (all symmetric except X, which is treated as general).
/// Every call returns false on a device error; the caller then falls back to the host loop.
class QMGpuSCF {
public:
    virtual ~QMGpuSCF() = default;
    virtual std::string name() const = 0;
    virtual int n() const = 0;
    /// Per SCF: core Hamiltonian H, overlap S, orthonormaliser X (F C = S C eps is solved as
    /// eig(X^T F X), C = X C'), number of doubly occupied orbitals, number of DIIS slots.
    virtual bool setup(const double* H, const double* S, const double* X, int n_occ, int nslots) = 0;
    /// Current density P (the start guess).
    virtual bool setDensity(const double* P) = 0;
    /// F = H + J - 1/2 K from the current density with the device J/K builder `jk`:
    /// full: J, K := JK(P); else J += JK(P - P_ref), K += ... (incremental, linear in P);
    /// then P_ref := P. Same arithmetic as the host loop's fockFromDensity().
    virtual bool buildFock(QMGpuJK* jk, bool full) = 0;
    /// Copy F and P into DIIS slot `slot` and form there e = X^T (F P S - S P F) X; max|e| out.
    virtual bool diisPush(int slot, double& max_err) = 0;
    /// B_ij = sum_kl e_i(kl) e_j(kl) for the given slots (m x m row-major).
    virtual bool diisErrorGram(const int* slots, int m, double* B) = 0;
    /// ADIIS model with n = slots[m-1] the newest entry:
    ///   a_i = sum (P_i - P_n) F_n,  M_ij = sum (P_i - P_n)(F_j - F_n)  (M not symmetrised).
    virtual bool adiisData(const int* slots, int m, double* a, double* M) = 0;
    /// Fock matrix to diagonalise: sum_i c_i F_{slots[i]}, or the plain F when m == 0.
    virtual bool selectFock(const int* slots, const double* c, int m) = 0;
    /// Diagonalise the selected Fock matrix, P_new = 2 C_occ C_occ^T; returns ||P_new - P||_F
    /// and E_elec = 1/2 Tr(P_new (H + F)) with the (non-extrapolated) F of buildFock().
    virtual bool solveAndDensity(double& delta_P, double& e_elec) = 0;
    /// P := damping P + (1 - damping) P_new.
    virtual bool acceptDensity(double damping) = 0;
    /// P_new (newest = true) or the current P, n x n row-major.
    virtual bool downloadDensity(bool newest, double* P) const = 0;
    /// Orbitals of the last solveAndDensity(): C (n x n, row-major, columns = MOs), eps (n).
    virtual bool downloadOrbitals(double* C, double* eps) const = 0;
    /// Cumulative device timings since setup() (one ASCII line).
    virtual std::string stats() const { return std::string(); }
};

} // namespace qmint
