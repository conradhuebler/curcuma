/*
 * <GFNFFGPUComputationalMethod — CUDA instantiation of the shared GFN-FF GPU wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for gfnff GPU path.
 * Claude Generated (Sep 2026): the wrapper logic lives in the shared template
 * GFNFFGpuMethodImpl<Backend> (gfnff_gpu_method_impl.h).  This translation unit — the
 * only one compiled by nvcc — provides the CUDA traits' one runtime call and the
 * explicit instantiation, so all wrapper code is emitted here exactly as before.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu cuda
 *   ./curcuma -sp mol.xyz -method gfnff -gpu auto   # GPU if available
 */

#ifdef USE_CUDA

#include "gfnff_gpu_method.h"

#include <cuda_runtime.h>

// ---------------------------------------------------------------------------
// CUDA backend traits: the single runtime-API call the wrapper needs
// ---------------------------------------------------------------------------

void GFNFFCudaBackend::downloadDoubles(double* host, const double* device, int n)
{
    cudaMemcpy(host, device, n * sizeof(double), cudaMemcpyDeviceToHost);
}

int GFNFFCudaBackend::deviceCount()
{
    int n = 0;
    return cudaGetDeviceCount(&n) == cudaSuccess ? n : 0;
}

bool GFNFFCudaBackend::setDevice(int device)
{
    int cur = -1;
    if (cudaGetDevice(&cur) == cudaSuccess && cur == device) return true;
    return cudaSetDevice(device) == cudaSuccess;
}

// ---------------------------------------------------------------------------
// Explicit instantiation of the shared wrapper for the CUDA backend
// ---------------------------------------------------------------------------

template class GFNFFGpuMethodImpl<GFNFFCudaBackend>;

#endif // USE_CUDA
