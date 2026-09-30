/*
 * <EEQSolverHip - ROCm include shim for the shared EEQSolverGPU declaration>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): cuda/eeq_solver_gpu.h is backend-neutral (Pimpl: neither
 * cuSOLVER nor rocSOLVER appears in the interface), so the ROCm build uses that one
 * declaration instead of a hipified copy that has to be kept in sync by hand.
 *
 * Why the rename instead of one shared class name: the HIP port was generated with
 * hipify, so rocm/eeq_solver_hip.hiph (included by gfnff_rocm.hip) defines its members
 * as EEQSolverHip::... and only compiles with a ROCm SDK.  Keeping the name also keeps
 * the two plugins' exported symbols distinct - gpu_plugin.cpp dlopen's
 * libcurcuma_cuda.so and libcurcuma_rocm.so with RTLD_GLOBAL, so identical strong
 * symbols in both would bind across plugin boundaries.  The two object-like macros
 * below give the shared declaration its ROCm names; nothing else in this build spells
 * "EEQSolverGPU".
 */

#pragma once

#define EEQSolverGPU     EEQSolverHip
#define EEQSolverGPUImpl EEQSolverHipImpl

#include "../cuda/eeq_solver_gpu.h"
