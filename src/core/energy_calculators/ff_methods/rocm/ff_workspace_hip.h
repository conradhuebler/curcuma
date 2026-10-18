/*
 * <FFWorkspaceHip - ROCm include shim for the shared FFWorkspaceGPU declaration>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): cuda/ff_workspace_gpu.h is backend-neutral (Pimpl: no
 * CUDA/HIP type appears in the interface), so the ROCm build uses that one declaration
 * instead of a hipified copy that has to be kept in sync by hand.
 *
 * Why the rename instead of one shared class name: the HIP port was generated with
 * hipify, so rocm/gfnff_rocm.hip defines its members as FFWorkspaceHip::... (7k lines,
 * only compilable with a ROCm SDK).  Keeping the name also keeps the two plugins'
 * exported symbols distinct - gpu_plugin.cpp dlopen's libcurcuma_cuda.so and
 * libcurcuma_rocm.so with RTLD_GLOBAL, so identical strong symbols in both would bind
 * across plugin boundaries.  The two object-like macros below give the shared
 * declaration its ROCm names; nothing else in this build spells "FFWorkspaceGPU".
 */

#pragma once

#define FFWorkspaceGPU     FFWorkspaceHip
#define FFWorkspaceGPUImpl FFWorkspaceHipImpl

#include "../cuda/ff_workspace_gpu.h"
