/*
 * <GFN-FF GPU SoA - ROCm include shim>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): the SoA layout is backend-neutral (runtime API spelled
 * through ff_methods/gpu_rt.h), so there is only ONE definition, in cuda/gfnff_soa.h.
 * The ROCm-only gather-kernel CSR members live there behind #ifdef USE_ROCM.
 */

#pragma once

#include "../cuda/gfnff_soa.h"
