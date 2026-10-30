/*
 * < Device seam of the native QM engine: Kohn-Sham exchange-correlation on the GPU >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, GPU KS-DFT G1-G3 of docs/QM_GPU_ROADMAP.md section 4). Light header (no Eigen): seen by
// cuda/qm_gpu_xc.cu, cuda/qm_gpu_scf.cu (via qm_gpu_scf_backend.h) and the core. The point kernels are those of
// qm_xc_kernels.h (one source for CPU and GPU).

#include <string>

namespace qmint {

/// Everything the device XC object needs, host arrays (the plugin copies them). Built by QMEngine::ensureGpuXC from
/// the molecular grid (qm_grid.h), AOEvaluator::flatten and the functional (qm_xc.h).
struct QMGpuXCSetup {
    int nbf = 0;                      ///< active basis dimension (V_xc, P)
    int npts = 0;                     ///< grid points
    int natoms = 0;
    // grid (Bohr): points npts x 3 row-major, Becke-partitioned weights, quadrature weights before the partition,
    // owning atom per point (points of an atom contiguous, atom_offset natoms + 1)
    const double* points = nullptr;
    const double* weights = nullptr;
    const double* quad_weights = nullptr;
    const int* point_atom = nullptr;
    const int* atom_offset = nullptr;
    // AO shells (qmgrid::FlatShells)
    int nshell = 0, max_ncomp = 0, nexps = 0, ncoef = 0, nlmn = 0, nT = 0, ncomp_total = 0;
    const double* sh_xyz = nullptr;   ///< 3 per shell
    const double* sh_amin = nullptr;
    const int* sh_nprim = nullptr;
    const int* sh_prim_off = nullptr;
    const int* sh_ncomp = nullptr;
    const int* sh_comp_off = nullptr;
    const int* sh_out_first = nullptr;
    const int* sh_out_width = nullptr;
    const int* sh_T_off = nullptr;
    const int* sh_atom = nullptr;     ///< atom of the shell (basis part of the gradient)
    const double* exps = nullptr;
    const double* coef = nullptr;
    const int* coef_off = nullptr;    ///< per component
    const int* lmn = nullptr;         ///< 3 per component
    const double* T = nullptr;
    // functional: terms (qmxc::XCKernel as int, weights), GGA flag, density threshold (already the effective one)
    int nterms = 0;
    const int* term_kernel = nullptr;
    const double* term_weight = nullptr;
    int gga = 0;
    double rho_min = 0.0;
    // grid response (Becke weight derivative): atom positions (Bohr, natoms x 3), 1/R_ij and a_ij (natoms x natoms)
    const double* atom_xyz = nullptr;
    const double* inv_dist = nullptr;
    const double* a_adjust = nullptr;
    int stratmann = 0;
    // memory: > 0 caps the device memory the object may use (MB, -qm.gpu_xc_max_memory_mb); block_points > 0 fixes
    // the number of grid points per streamed block (testing), else chosen from the memory budget
    double max_device_mb = 0.0;
    int block_points = 0;
    // AO screening (Claude Generated, Oct 2026, -qm.xc_ao_threshold; appended so that the layout above is unchanged):
    // the blocks of qmgrid::BlockScreen on the (sorted) grid above, block b = points [screen_begin[b], + screen_count[b])
    // with the kept shells screen_shells[screen_shell_off[b] .. screen_shell_off[b+1]) (ascending). The device merges
    // consecutive blocks into streamed groups of at most block_points points (0: kScreenGroupPoints of
    // cuda/qm_gpu_xc.cu) and works on the union of their shells. screen_nblocks = 0: no screening.
    int screen_nblocks = 0;
    const int* screen_begin = nullptr;
    const int* screen_count = nullptr;
    const int* screen_shell_off = nullptr;
    const int* screen_shells = nullptr;
};

/// E_xc parts of one V_xc build (sums over the grid, as qmxc::XCResult).
struct QMGpuXCScalars {
    double ex = 0.0, ec = 0.0, nelec = 0.0;
    long long points_used = 0;
};

/// Closed-shell Kohn-Sham exchange-correlation on a device: AO values (deriv 0/1/2, spherical transform), density
/// and density gradient, the functional, V_xc and E_xc, and the nuclear XC gradient incl. the grid response.
/// Same quantities as qmxc::buildVxc / qmxc::buildXCGradient / qmgrid::gridWeightGradient. The grid stays resident;
/// AO values exist only per block of points (streamed), so memory is O(npts + block * nbf + nbf^2).
class QMGpuXC {
public:
    virtual ~QMGpuXC() = default;
    virtual std::string name() const = 0;
    virtual int nbf() const = 0;
    /// V_xc (n x n, symmetric; row-major = column-major) and the scalars for the host density P.
    virtual bool vxc(const double* P, double* V, QMGpuXCScalars& out) = 0;
    /// Same with DEVICE pointers (P, V: n x n on this object's device; used by the resident SCF loop). Synchronous.
    virtual bool vxcDevice(const double* dP, double* dV, QMGpuXCScalars& out) = 0;
    /// Nuclear gradient (Eh/Bohr, natoms x 3 row-major each, overwritten): basis-function part, point part of the
    /// grid response (points move with their atom), Becke weight derivative (only when response is true, else zero).
    virtual bool gradient(const double* P, bool response, double* g_basis, double* g_points, double* g_weights,
                          long long& points_used) = 0;
    /// Points per streamed block chosen from the memory budget.
    virtual int blockPoints() const = 0;
    /// Becke weight derivative alone (grid response part 2): g (natoms x 3, row-major, overwritten) =
    /// sum_p f_p d w_p / d R_A for the HOST array f (npts), as qmgrid::gridWeightGradient. Any number of atoms
    /// (O(natoms^2) device scratch per thread block, streamed in waves). Claude Generated (Oct 2026).
    virtual bool weightGradient(const double* f, double* g) = 0;
    /// One ASCII line of cumulative device timings.
    virtual std::string stats() const { return std::string(); }
};

}  // namespace qmint
