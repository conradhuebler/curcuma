/*
 * <Embedded SPIR-V for the Vulkan FP64 two-sided Jacobi eigensolver>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (2026-06): SPIR-V byte arrays compiled from the GLSL compute
 * shaders in this directory (angles/col/row/vec). Committed so the build needs no
 * glslc at compile time. Regenerate with ./compile_shaders.sh after editing a .comp.
 */

#pragma once

#ifdef USE_VULKAN

#include <cstdint>

namespace curcuma {
namespace vk {
namespace shaders {

static const uint32_t angles_spv[] =
#include "angles.spv.inc"
;
static const uint32_t col_spv[] =
#include "col.spv.inc"
;
static const uint32_t row_spv[] =
#include "row.spv.inc"
;
static const uint32_t vec_spv[] =
#include "vec.spv.inc"
;
static const uint32_t gemm_spv[] =
#include "gemm.spv.inc"
;
static const uint32_t scale_cols_spv[] =
#include "scale_cols.spv.inc"
;
static const uint32_t fock_spv[] =
#include "fock.spv.inc"
;
static const uint32_t popband_spv[] =
#include "popband.spv.inc"
;
// Stage 3: on-device integral build (CN / self-energy / overlap+H0 / Coulomb gamma).
static const uint32_t cn_spv[] =
#include "cn.spv.inc"
;
static const uint32_t self_energy_spv[] =
#include "self_energy.spv.inc"
;
static const uint32_t overlap_h0_spv[] =
#include "overlap_h0.spv.inc"
;
static const uint32_t gamma_spv[] =
#include "gamma.spv.inc"
;
// Stage 4: on-device GFN1 nuclear gradient (repulsion / Coulomb / H0-Pulay+CN).
// Per-atom GATHER kernels — no FP64 atomics (Vulkan has none). Claude Generated (V-AP1).
static const uint32_t grad_rep_spv[] =
#include "grad_rep.spv.inc"
;
static const uint32_t grad_coulomb_spv[] =
#include "grad_coulomb.spv.inc"
;
static const uint32_t grad_pulay_spv[] =
#include "grad_pulay.spv.inc"
;
// Stage 3m: GFN2 AO multipole integrals (dp_int/qp_int) on device. Claude Generated (V-AP2).
static const uint32_t multipole_ints_spv[] =
#include "multipole_ints.spv.inc"
;
// Stage 2b: GFN2 device-resident multipole SCF — anisotropic Fock term + atomic moments.
static const uint32_t fock_multipole_spv[] =
#include "fock_multipole.spv.inc"
;
static const uint32_t multipole_moments_spv[] =
#include "multipole_moments.spv.inc"
;
// X-AP3: FP32 two-sided Jacobi (opt-in -scf_mixed_precision; FP32 is ~16× FP64 on an iGPU).
static const uint32_t angles_f32_spv[] =
#include "angles_f32.spv.inc"
;
static const uint32_t col_f32_spv[] =
#include "col_f32.spv.inc"
;
static const uint32_t row_f32_spv[] =
#include "row_f32.spv.inc"
;
static const uint32_t vec_f32_spv[] =
#include "vec_f32.spv.inc"
;
// GPU Householder tridiagonalization eigensolver (replaces the cyclic Jacobi):
// matvec + symmetric rank-2 trailing-block update + reflector back-transform.
static const uint32_t tri_matvec_spv[] =
#include "tri_matvec.spv.inc"
;
static const uint32_t tri_rank2_spv[] =
#include "tri_rank2.spv.inc"
;
static const uint32_t tri_applyl_spv[] =
#include "tri_applyl.spv.inc"
;
// EIG-1: fully-GPU reduction — per-column reflector build (tri_house) + w-vector build
// (tri_kw), so the host stays out of the per-column loop (no per-step fence sync).
static const uint32_t tri_house_spv[] =
#include "tri_house.spv.inc"
;
static const uint32_t tri_kw_spv[] =
#include "tri_kw.spv.inc"
;
// EIG-4: FP32 mixed-precision variants of the tridiagonalization eigensolver kernels
// (used far from SCF convergence; FP64 corrector near it keeps the converged energy FP64).
static const uint32_t tri_house_f32_spv[] =
#include "tri_house_f32.spv.inc"
;
static const uint32_t tri_matvec_f32_spv[] =
#include "tri_matvec_f32.spv.inc"
;
static const uint32_t tri_kw_f32_spv[] =
#include "tri_kw_f32.spv.inc"
;
static const uint32_t tri_rank2_f32_spv[] =
#include "tri_rank2_f32.spv.inc"
;
static const uint32_t tri_applyl_f32_spv[] =
#include "tri_applyl_f32.spv.inc"
;
// EIG-2B: WY-blocked eigenvector back-transform — scatter the compact panel reflectors
// to full layout (tri_vfull), build the b×b compact-WY factor T (wy_buildt), and apply
// the panel as three rectangular FP64 GEMMs (gemm_g) instead of n−3 per-reflector passes.
static const uint32_t gemm_g_spv[] =
#include "gemm_g.spv.inc"
;
static const uint32_t tri_vfull_spv[] =
#include "tri_vfull.spv.inc"
;
static const uint32_t wy_buildt_spv[] =
#include "wy_buildt.spv.inc"
;

} // namespace shaders
} // namespace vk
} // namespace curcuma

#endif // USE_VULKAN
