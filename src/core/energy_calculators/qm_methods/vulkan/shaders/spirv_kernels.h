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

} // namespace shaders
} // namespace vk
} // namespace curcuma

#endif // USE_VULKAN
