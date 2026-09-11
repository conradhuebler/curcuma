/*
 * <gpu_rt.h - CUDA / HIP runtime-API compatibility layer>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): one neutral spelling for the handful of runtime-API
 * calls that the SHARED GFN-FF GPU headers need, so those headers exist exactly once
 * instead of twice (cuda/*.h + a hipified rocm/*_hip.h copy).
 *
 * The CUDA and HIP runtime APIs are name-for-name identical for everything used here
 * (allocate / copy / memset / stream handle / error string), so the mapping is a pure
 * rename:
 *
 *      gpuMalloc      -> cudaMalloc      / hipMalloc
 *      gpuFree        -> cudaFree        / hipFree
 *      gpuMemcpy      -> cudaMemcpy      / hipMemcpy
 *      gpuMemcpyAsync -> cudaMemcpyAsync / hipMemcpyAsync
 *      gpuMemsetAsync -> cudaMemsetAsync / hipMemsetAsync
 *      gpuStream_t    -> cudaStream_t    / hipStream_t
 *      gpuError_t     -> cudaError_t     / hipError_t
 *      ...
 *
 * Backend selection: the HIP toolchain always defines __HIP_PLATFORM_AMD__ (hipcc for
 * the device TU, and the g++ compile of the ROCm plugin sources — see CMakeLists), so
 * that macro alone decides.  USE_CUDA / USE_ROCM are only used for the "is any GPU
 * backend compiled at all" guard; they are plugin-private definitions and never both
 * present in the same translation unit.
 *
 * This header is HOST-side only (thin inline wrappers, no __device__ code); device
 * kernels keep using the native spelling of their own toolchain.
 */

#pragma once

#if defined(USE_CUDA) || defined(USE_ROCM)

#include <cstddef>

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIP__)

#include <hip/hip_runtime.h>

/// Literal prefix used in runtime error messages ("hipMalloc failed: ...").
#define GPU_RT_PREFIX "hip"

using gpuError_t  = hipError_t;
using gpuStream_t = hipStream_t;

constexpr gpuError_t gpuSuccess = hipSuccess;

constexpr hipMemcpyKind gpuMemcpyHostToDevice   = hipMemcpyHostToDevice;
constexpr hipMemcpyKind gpuMemcpyDeviceToHost   = hipMemcpyDeviceToHost;
constexpr hipMemcpyKind gpuMemcpyDeviceToDevice = hipMemcpyDeviceToDevice;

using gpuMemcpyKind = hipMemcpyKind;

inline gpuError_t gpuMalloc(void** ptr, std::size_t bytes) { return hipMalloc(ptr, bytes); }
inline gpuError_t gpuFree(void* ptr) { return hipFree(ptr); }
inline gpuError_t gpuMemcpy(void* dst, const void* src, std::size_t bytes, gpuMemcpyKind kind)
{
    return hipMemcpy(dst, src, bytes, kind);
}
inline gpuError_t gpuMemcpyAsync(void* dst, const void* src, std::size_t bytes,
                                 gpuMemcpyKind kind, gpuStream_t stream)
{
    return hipMemcpyAsync(dst, src, bytes, kind, stream);
}
inline gpuError_t gpuMemsetAsync(void* ptr, int value, std::size_t bytes, gpuStream_t stream)
{
    return hipMemsetAsync(ptr, value, bytes, stream);
}
inline const char* gpuGetErrorString(gpuError_t err) { return hipGetErrorString(err); }

#else // CUDA

#include <cuda_runtime.h>

/// Literal prefix used in runtime error messages ("cudaMalloc failed: ...").
#define GPU_RT_PREFIX "cuda"

using gpuError_t  = cudaError_t;
using gpuStream_t = cudaStream_t;

constexpr gpuError_t gpuSuccess = cudaSuccess;

constexpr cudaMemcpyKind gpuMemcpyHostToDevice   = cudaMemcpyHostToDevice;
constexpr cudaMemcpyKind gpuMemcpyDeviceToHost   = cudaMemcpyDeviceToHost;
constexpr cudaMemcpyKind gpuMemcpyDeviceToDevice = cudaMemcpyDeviceToDevice;

using gpuMemcpyKind = cudaMemcpyKind;

inline gpuError_t gpuMalloc(void** ptr, std::size_t bytes) { return cudaMalloc(ptr, bytes); }
inline gpuError_t gpuFree(void* ptr) { return cudaFree(ptr); }
inline gpuError_t gpuMemcpy(void* dst, const void* src, std::size_t bytes, gpuMemcpyKind kind)
{
    return cudaMemcpy(dst, src, bytes, kind);
}
inline gpuError_t gpuMemcpyAsync(void* dst, const void* src, std::size_t bytes,
                                 gpuMemcpyKind kind, gpuStream_t stream)
{
    return cudaMemcpyAsync(dst, src, bytes, kind, stream);
}
inline gpuError_t gpuMemsetAsync(void* ptr, int value, std::size_t bytes, gpuStream_t stream)
{
    return cudaMemsetAsync(ptr, value, bytes, stream);
}
inline const char* gpuGetErrorString(gpuError_t err) { return cudaGetErrorString(err); }

#endif // backend

#endif // USE_CUDA || USE_ROCM
