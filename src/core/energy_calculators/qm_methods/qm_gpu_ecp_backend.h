/*
 * < Device seam of the native QM engine: ECP matrix and ECP gradient (semi-local def2-ECP) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Nov 2026, ECP on the device, docs/QM_GPU_ROADMAP.md section 4.6). Own header (pattern of
// qm_gpu_grad1e_backend.h), seen only by cuda/qm_gpu_ecp.cu and qm_ecp.cpp. No CUDA types: the core drives the
// device through this interface, the plugin implements it.
//
// Division of work (qm_ecp.cpp "device path"): the host keeps the radial-order loop and its acceptance test,
// the angular tables T_ijk,lambda,lm (Lebedev, per atom and ECP centre) and every decision; the device does the
// heavy parts at one quadrature order: the local three-centre overlaps, the radial factors, the angular projections
// F_f,lm(r_k), the weighted GEMM V = sum F diag(w r^2 dU_l) F^T and, for the gradient, X = P FW and the row dots
// with the derivative-function projections. All arrays are host memory, matrices row-major.
//
// Layout check: both structs carry their size in `bytes` (set by the host to sizeof), the plugin refuses a
// mismatch, and the entry name carries the version (curcuma_<backend>_create_qm_ecp_v1).

#include <string>

namespace qmecp {

/// Functions of one device object: `nbasis` basis functions (cartesian, the order of QMEngine::m_gto_basis), for a
/// gradient object followed by the three centre-derivative functions of each, index nbasis + 3 i + d. Each function
/// is a list of cartesian primitives (terms) on one atom.
struct QMGpuEcpSetup {
    int bytes = (int)sizeof(QMGpuEcpSetup);
    int nbasis = 0;
    int nfn = 0;                       ///< nbasis (matrix object) or 4 nbasis (gradient object)
    const int* fn_atom = nullptr;      ///< per function: atom index
    const int* fn_term = nullptr;      ///< nfn + 1 offsets into the term arrays
    int nterm = 0;
    const double* term_c = nullptr;      ///< coefficient (primitive and contraction normalisation folded in)
    const double* term_alpha = nullptr;  ///< exponent
    const int* term_p = nullptr;         ///< 3 per term: cartesian powers
    const int* term_slot = nullptr;      ///< per term: radial slot = unique (atom, exponent)
    int nslot = 0;
    const int* slot_atom = nullptr;
    const double* slot_alpha = nullptr;
    int natoms = 0;                      ///< atoms of the molecule (functions and centres index them)
    const double* fn_xyz = nullptr;      ///< per function: centre (Bohr), 3 per function
};

/// One ECP centre as the device needs it (setCentre). Angular tables per atom from the host (qm_ecp.cpp
/// angularTable with the atom's largest power, the same tables as the CPU path): T[(mono * (lam_max+1) + lam) *
/// lproj^2 + lm], lam_max = ang_nmax[a] + lproj, mono = kern::monoIndex.
struct QMGpuEcpCentre {
    int bytes = (int)sizeof(QMGpuEcpCentre);
    double C[3] = { 0, 0, 0 };           ///< centre (Bohr)
    int nlocal = 0;
    const double* local_coef = nullptr;  ///< local channel U_L: c_k exp(-zeta_k r^2)
    const double* local_zeta = nullptr;
    int lproj = 0;                       ///< semi-local channels l < lproj (= ElementECP::lmax), <= 5 on the device
    const int* semi_off = nullptr;       ///< lproj + 1 offsets into semi_coef/semi_zeta
    const double* semi_coef = nullptr;
    const double* semi_zeta = nullptr;
    double rcut = 0.0;                   ///< radius beyond which every semi-local channel is below 1e-20
    const double* atom_D = nullptr;      ///< 3 per atom: D = A - C (Bohr), as the CPU computes it
    const double* atom_d = nullptr;      ///< per atom: |D|
    const long long* ang_off = nullptr;  ///< per atom: offset into ang_T, -1 = no functions on the atom
    const int* ang_nmax = nullptr;       ///< per atom: largest cartesian power of its functions
    const double* ang_T = nullptr;
    long long ang_size = 0;
};

class QMGpuECP {
public:
    virtual ~QMGpuECP() = default;
    virtual std::string name() const = 0;
    /// Upload one centre (tables, channels); the calls below refer to it. False on a device error.
    virtual bool setCentre(const QMGpuEcpCentre& c) = 0;
    /// Local part <a|U_L|b> of the current centre, nbasis x nbasis, overwritten.
    virtual bool local(double* V) = 0;
    /// Semi-local part of the current centre at the npts-point radial order, nbasis x nbasis, overwritten.
    virtual bool semilocal(int npts, double* V) = 0;
    /// Gradient pieces of the current centre (gradient object only): G(a, d) = sum_b P_ab <d_d a| U_C |b>,
    /// local + semi-local at npts points; P nbasis x nbasis, G nbasis x 3, overwritten.
    virtual bool gradient(int npts, const double* P, double* G) = 0;
    /// One ASCII line of device timings, empty if none.
    virtual std::string stats() const { return std::string(); }
};

}  // namespace qmecp
