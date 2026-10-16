/*
 * < Device seam of the native QM engine: 3-index part of the RI gradient >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, RI gradient). Own header (pattern of qm_gpu_ri_backend.h) so only the new
// cuda/qm_gpu_ri_grad.cu and the core see it; the heavy ERI kernel TUs are not recompiled.
//
// The device evaluates   g_X = sum_{(A,B) canonical} w_AB sum_C sum_{a,b,c} W^{ABC}_{abc} d/dX (ab|c)
// for the derivative on the centres of A and B, the centre of C getting -(dA + dB) (translational
// invariance), i.e. the quantity of qmri::RIJK::threeIndexGradient. The host (qmri::RIJK) prepares
// everything that is not an integral: shells, Hermite tables of the primitive pairs (bra pairs raised
// by one on both shells, auxiliary kets (C, unit s) not raised), the screened task list (pair, aux shell)
// and the cartesian weight blocks W (3-index density back-transformed per block). FP64 only.

namespace qmint {

struct QMGpuRIGradInput {
    int natoms = 0;
    // Shells: orbital shells 0..aux_base-1, auxiliary shells aux_base..nshell-1.
    int nshell = 0, aux_base = 0;
    const int* sh_L = nullptr;
    const int* sh_atom = nullptr;
    const int* sh_ncomp = nullptr;
    const int* sh_lmn_off = nullptr;   ///< into lmn (3 ints per component)
    const int* lmn = nullptr;
    const int* sh_nprim = nullptr;
    const int* sh_coef_off = nullptr;  ///< into coef: [component][primitive], normalisation included
    const double* coef = nullptr;
    const int* sh_exp_off = nullptr;   ///< into exps
    const double* exps = nullptr;
    int nlmn = 0, ncoef = 0, nexp = 0;
    // Bra pairs (canonical orbital shell pairs A <= B), Hermite tables raised by one on both shells.
    int npairs = 0;
    const int* pair_A = nullptr;
    const int* pair_B = nullptr;
    const int* pair_pp_off = nullptr;  ///< first primitive pair
    const int* pair_npp = nullptr;
    // Ket primitive pairs of every auxiliary shell (local index c = shell - aux_base)
    const int* aux_pp_off = nullptr;
    const int* aux_npp = nullptr;
    // Primitive pairs: ia, ib, (p, Px, Py, Pz), offset of the Ex table (Ey, Ez follow at +sz, +2sz)
    int npp = 0;
    const int* pp_ia = nullptr;
    const int* pp_ib = nullptr;
    const double* pp_geom = nullptr;   ///< 4 per primitive pair
    const long long* pp_E_off = nullptr;
    long long nE = 0;
    const double* E = nullptr;
    // Tasks per bra pair k: task_start[k] .. task_start[k+1]-1 (CSR), auxiliary shell (local) and weights
    const long long* task_start = nullptr;  ///< npairs + 1
    long long ntask = 0;
    const int* task_aux = nullptr;
    const long long* task_woff = nullptr;
    long long nW = 0;
    const double* W = nullptr;              ///< per task [a][b][c] cartesian (nA x nB x nC)
};

struct QMGpuRIGradOutput {
    double ms_upload = 0.0, ms_kernel = 0.0;
};

} // namespace qmint
