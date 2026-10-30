/*
 * < Device seam of the native QM engine (ERI build, J/K builds) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include <cstddef>
#include <string>

// Claude Generated (Oct 2026, Q-G1..Q-G4 of docs/QM_GPU_ROADMAP.md): backend-free interface the
// QMEngine drives a GPU through. The concrete classes live in a runtime plugin
// (libcurcuma_cuda.so, entries curcuma_cuda_create_qm_jk / curcuma_cuda_create_qm_eri_jk /
// curcuma_cuda_create_qm_direct_jk / curcuma_cuda_create_qm_scf_ks), so the core carries no CUDA symbol.
namespace qmint {

/// Coulomb/exchange builder on a device: either holding the stored ERI tensor (Q-G1/Q-G2) or
/// integral-direct (Q-G3: every screened shell quartet recomputed per build and contracted
/// with P straight away, as qmint::DirectJK; O(n^2) device memory).
/// Same definitions as buildCoulomb/buildExchange in qm_integrals.hpp:
///   J_mn = sum_ls P_ls (mn|ls),  K_mn = sum_ls P_ls (ml|ns)
/// P, J, K are n x n row-major, contiguous.
class QMGpuJK {
public:
    virtual ~QMGpuJK() = default;
    /// Returns false on a device error; the caller then falls back to the CPU.
    virtual bool build(const double* P, double* J, double* K) = 0;
    virtual std::string name() const = 0;
    virtual int n() const = 0;
    /// Copy the device tensor (n^4 doubles, layout of ERITensor) to the host. Used by the
    /// CPU-vs-GPU tensor test; false if not supported or on a device error. (Q-G2)
    virtual bool downloadERI(double* /*host*/) const { return false; }
    /// One line of device-side timings (context, uploads, kernels), empty if none. (Q-G2)
    virtual std::string stats() const { return std::string(); }
    /// Integral-direct builders: shell quartets computed in the last build() (Q-G3), else 0.
    virtual long long lastQuartets() const { return 0; }
    /// Q-G4 (resident SCF): same as build(), but P, J, K are DEVICE pointers (n x n, on the
    /// device this builder runs on); J and K come back symmetrised exactly as from build().
    /// Synchronous: the result is complete on return. False = unsupported or device error.
    virtual bool buildDevice(const double* /*dP*/, double* /*dJ*/, double* /*dK*/) { return false; }
    /// WP11 (MP2 on the device, Claude Generated Oct 2026): device pointer of the stored tensor
    /// (n^4 doubles, ERITensor layout, on the device this builder runs on), nullptr if the
    /// builder holds none (integral-direct) or the backend does not expose it. Read-only.
    virtual const double* deviceERI() const { return nullptr; }
    /// Integral-direct builders (Claude Generated, Oct 2026): density-weighted screening threshold
    /// of the following builds (primitive threshold thr * 1e-3) and whether the primitive screen is
    /// relative to the quartet's density factor (as qmint::DirectJK::build). False = not supported.
    virtual bool setScreening(double /*thr*/, bool /*prim_density*/) { return false; }
    /// Integral-direct builders (Q-G5, Claude Generated Oct 2026): FP32 integrals for the following
    /// builds (early SCF phase of -qm.gpu_eri_precision mixed); FP64 accumulation of J/K. False =
    /// not supported (the builder stays FP64). The SCF accepts convergence only on FP64 builds.
    virtual bool setFP32(bool /*on*/) { return false; }
    /// FP32:FP64 throughput ratio of the device (e.g. 64 on a GeForce RTX 5090, 2 on an HPC card);
    /// 0 = unknown. Used by -qm.gpu_eri_precision auto.
    virtual int fp32Ratio() const { return 0; }
};

/// Flat, pointer-only description of everything the device ERI kernel needs (Q-G2).
/// Filled on the host by qmint::EriGpuTables::view() (qm_integrals.cpp) from exactly the
/// same shell-pair tables, Schwarz factors and primitive ordering the CPU buildERI uses,
/// so the screening decisions are identical. All arrays are host memory; the plugin
/// copies them. Indices refer to:
///   shells       s = 0..nshell-1 (same grouping as the CPU kernel)
///   shell pairs  i = 0..npair-1, A <= B, CPU order
///   prim pairs   k = pr_prim_off[i] .. + pr_nprim[i]-1, sorted by Schwarz factor (largest first)
///   quartets     canonical (i, j), j >= i, Schwarz-screened, grouped by angular-momentum class
struct QMGpuEriTables {
    int n_out = 0;           ///< active basis dimension; the tensor has n_out^4 entries
    int nshell = 0, npair = 0, nprimpair = 0;
    double prim_thr = 0.0;   ///< primitive-quartet threshold (screening * 1e-3, as on the CPU)
    double pref = 0.0;       ///< 2 pi^(5/2), computed on the host as in the CPU kernel
    // per shell
    const int* sh_L = nullptr;
    const int* sh_ncomp = nullptr;
    const int* sh_nprim = nullptr;
    const int* sh_lmn_off = nullptr;   ///< into lmn (3 ints per component)
    const int* sh_coef_off = nullptr;  ///< into coef, layout [component][primitive]
    const int* sh_out_first = nullptr; ///< first index in the active basis
    const int* sh_out_width = nullptr; ///< number of active functions
    const int* sh_T_off = nullptr;     ///< into T (ncomp x width, row-major); -1 = identity
    const int* lmn = nullptr; int n_lmn = 0;
    const double* coef = nullptr; int n_coef = 0;
    const double* T = nullptr; int n_T = 0;
    // per shell pair
    const int* pr_A = nullptr;
    const int* pr_B = nullptr;
    const int* pr_prim_off = nullptr;
    const int* pr_nprim = nullptr;
    // per primitive pair
    const int* pp_ia = nullptr;
    const int* pp_ib = nullptr;
    const double* pp_p = nullptr;      ///< composite exponent
    const double* pp_P = nullptr;      ///< product centre, 3 per pair
    const double* pp_schwarz = nullptr;
    const long long* pp_E_off = nullptr; ///< into E: Ex, Ey, Ez consecutive, each (LA+LB+1)(LA+1)(LB+1)
    const double* E = nullptr; long long n_E = 0;
    // canonical quartets after shell-level screening, grouped by class
    const int* q_bra = nullptr;        ///< shell-pair index i
    const int* q_ket = nullptr;        ///< shell-pair index j >= i
    long long nquartet = 0;
    int nclass = 0;
    const long long* class_off = nullptr; ///< nclass+1 offsets into q_bra/q_ket
    const int* class_L = nullptr;      ///< 4 per class: LA LB LC LD
    // Boys table of qm_integrals.cpp (grid + Taylor), copied verbatim
    const double* boys_f = nullptr;
    int boys_npts = 0, boys_ntab = 0, boys_ntaylor = 0;
    double boys_dT = 0.0, boys_Tmax = 0.0;
    // Q-G3 (Oct 2026)
    double device_mem_limit_mb = 0.0;  ///< > 0: treat the card as having at most this much free memory
    double screening = 0.0;            ///< Schwarz / density-weighted threshold (-qm.eri_screening)
    const double* pr_schwarz = nullptr;  ///< per shell pair: sqrt(max |(ab|ab)|) (ShellPair::schwarz)
    /// Integral-direct only (built with list_quartets = false; q_bra/q_ket then empty):
    /// shell pairs grouped by pair class k = LA*3 + LB (9 classes), CPU order inside a class.
    const int* pc_off = nullptr;       ///< 10 offsets into pc_pairs
    const int* pc_pairs = nullptr;
};

} // namespace qmint
