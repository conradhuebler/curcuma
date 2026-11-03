/*
 * < CUDA integral-direct J/K: thread-per-quartet kernels for the large classes (Q-G3) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 20 2026). Seam between qm_gpu_direct.cu (builder, enumeration, launches)
// and qm_gpu_direct_md.cu (kernels for the classes with more than 27 cartesian components, which
// the thread kernel of qm_gpu_eri_kernels.cuh does not cover). Separate TU: the template instances
// compile in parallel with the other kernel files. See the header of qm_gpu_direct_md.cu.

#include "qm_gpu_eri_kernels.cuh"

#include <string>
#include <vector>

namespace qmgpu {

/// Per primitive pair k: the full Hermite expansion of every active component pair,
///   F[k][a'][b'][h] = sum_{a,b} T_A(a,a') T_B(b,b') c_a c_b E^x_t E^y_u E^z_v,  h = (t,u,v), t+u+v <= LA+LB,
/// contraction coefficients and the 6d->5d transform folded in (h in the simplex order of sidx).
struct MdTables {
    const double* F = nullptr;
    const long long* F_off = nullptr;  ///< per primitive pair
    // FP32 copies of F, the density and the Boys table (read by the RT = float instances, Q-G5)
    const float* F32 = nullptr;
    const float* P32 = nullptr;
    const float* boys_f32 = nullptr;
    long long nF = 0;     ///< number of F entries (set by mdInit)
    long long nBoys = 0;  ///< number of Boys table entries (set by mdInit)
};

/// Build the F tables on the device from the uploaded tables d (needs t for the pair -> primitive
/// map). false + reason if the basis is not supported (cartesian d, unexpected component order or
/// transform pattern) or on a device error; the builder then keeps its old kernels.
/// dropped (stage 6b): largest |T| entry outside the 5d pattern that was accepted as zero (<= 1e-14).
bool mdInit(const qmint::QMGpuEriTables* t, const DevTables& d, MdTables& m, std::vector<void*>& owned,
            std::string& why, size_t& bytes, double* dropped = nullptr);

/// Launch the thread-per-quartet kernel of class (LA LB | LC LD) on the quartets of src; false if
/// this class has no instance (the caller uses its old kernel).
bool mdLaunch(int LA, int LB, int LC, int LD, const DevTables& d, const MdTables& m, const QSrc& src, long long nq,
              cudaStream_t st);

/// FP32 phase (qm_gpu_direct_fp32.cu, Claude Generated Oct 20 2026): FP32 copies of F and of the
/// Boys table (once), a device buffer for the FP32 density; false + reason on a device error.
bool md32Init(const DevTables& d, MdTables& m, std::vector<void*>& owned, std::string& why, size_t& bytes);

/// FP32 copy of the density (n x n) into m.P32 on stream st.
void md32Density(const double* P, const MdTables& m, size_t n2, cudaStream_t st);

/// FP32 instance of mdLaunch (same classes, same screening; FP64 atomics into Jt, Kt).
bool mdLaunch32(int LA, int LB, int LC, int LD, const DevTables& d, const MdTables& m, const QSrc& src, long long nq,
                cudaStream_t st);

} // namespace qmgpu
