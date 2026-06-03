/*
 * Shared helper for the Stage-3 GPU integral component tests: pack the flattened
 * GpuBasisFlat / GpuH0Flat (from XTB::exportGpuBasis) into the device-facing
 * XtbGpuBasisData pointer bundle.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 3).
 */

#pragma once

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"

inline curcuma::xtb::gpu::XtbGpuContext::XtbGpuBasisData
makeGpuBasisData(const curcuma::xtb::GpuBasisFlat& bf, const curcuma::xtb::GpuH0Flat& hf)
{
    curcuma::xtb::gpu::XtbGpuContext::XtbGpuBasisData bd;
    bd.nat         = bf.nat;
    bd.nsh         = bf.nsh;
    bd.nao         = bf.nao;
    bd.is_gfn2     = bf.is_gfn2;
    bd.z           = bf.z.data();
    bd.sh2at       = bf.sh2at.data();
    bd.ang_sh      = bf.ang_sh.data();
    bd.iao_sh      = bf.iao_sh.data();
    bd.nao_sh      = bf.nao_sh.data();
    bd.sh_nprim    = bf.sh_nprim.data();
    bd.sh_prim_off = bf.sh_prim_off.data();
    bd.nprim_total = static_cast<int>(bf.prim_alpha.size());
    bd.prim_alpha  = bf.prim_alpha.data();
    bd.prim_coeff  = bf.prim_coeff.data();
    bd.sh_zeta     = bf.sh_zeta.data();
    bd.selfenergy  = hf.selfenergy.data();
    bd.kcn         = hf.kcn.data();
    bd.shpoly      = hf.shpoly.data();
    bd.valence     = bf.valence.empty() ? nullptr : bf.valence.data();
    bd.shell_hardness = bf.shell_hardness.empty() ? nullptr : bf.shell_hardness.data();
    bd.ao2at       = bf.ao2at.empty() ? nullptr : bf.ao2at.data();
    bd.ao2sh       = bf.ao2sh.empty() ? nullptr : bf.ao2sh.data();
    return bd;
}
