/*
 * < CUDA ERI kernels for quartet classes with f and g shells (stage 6b) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Nov 2026, stage 6b of docs/QM_GPU_ROADMAP.md). Seam between the two ERI users
// (qm_gpu_eri.cu: stored tensor, Q-G2; qm_gpu_direct.cu: integral-direct J/K, Q-G3) and the generic
// McMurchie-Davidson kernel of qm_gpu_eri_generic.cu for every quartet class that contains an f or g
// shell (L = 3, 4). The compile-time kernels of qm_gpu_eri_kernels.cuh / qm_gpu_direct_md_impl.cuh
// stay the only path for s, p, d classes (their code and their bits are unchanged); a basis without
// f/g shells never reaches this file.

#include "qm_gpu_eri_kernels.cuh"

#include <string>
#include <vector>

namespace qmgpu {

/// Device tables of the generic kernel. Per primitive pair k of shell pair (A, B):
///   F[k][a'][b'][h] = sum_{a,b} T_A(a,a') T_B(b,b') c_a c_b E^x_t E^y_u E^z_v,   h = (t,u,v), t+u+v <= LA+LB,
/// in the ACTIVE basis (cartesian -> spherical transform and contraction folded in, as the F tables of
/// qm_gpu_direct_md.cu), h in DEGREE order (all t+u+v = 0, then 1, ...; inside a degree by t, then u).
/// scratch: one slot of `slot` doubles per resident block (ket Hermite sums and the quartet block).
struct GenTables {
    const double* F = nullptr;
    const long long* F_off = nullptr;  ///< per primitive pair
    const int* sh_L = nullptr;         ///< per shell
    const int* width = nullptr;        ///< per shell: active width (sh_out_width)
    double* scratch = nullptr;
    long long slot = 0;                ///< doubles per block slot
    int nslots = 0;                    ///< = maximum grid of one launch
    size_t bytes = 0;                  ///< device memory of these tables
};

/// Highest shell angular momentum of the tables (the host builds pair classes with base 3 for L <= 2 and
/// base 5 for L <= 4, qmint::buildEriGpuTables).
inline int maxShellL(const qmint::QMGpuEriTables* t)
{
    int m = 0;
    for (int s = 0; s < t->nshell; ++s) m = t->sh_L[s] > m ? t->sh_L[s] : m;
    return m;
}
inline int pairClassBase(const qmint::QMGpuEriTables* t) { return maxShellL(t) <= 2 ? 3 : 5; }

/// Build the F tables and the scratch slots on the device (d: the uploaded common tables of t).
/// false + reason on a device error or out of memory.
bool genInit(const qmint::QMGpuEriTables* t, const DevTables& d, GenTables& g, std::vector<void*>& owned,
             std::string& why);

/// Launch class (LA LB | LC LD) on the quartets of src: direct = Q-G3 digest into d.J/d.K (and d.K2), else
/// write the stored tensor out (set8, Q-G2). All launches must go to ONE stream (they share the scratch).
bool genLaunch(int LA, int LB, int LC, int LD, const DevTables& d, const GenTables& g, const QSrc& src, long long nq,
               bool direct, double* out, cudaStream_t st);

} // namespace qmgpu
