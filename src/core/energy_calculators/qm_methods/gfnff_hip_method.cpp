/*
 * <GFNFFHipComputationalMethod — ROCm instantiation of the shared GFN-FF GPU wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for gfnff GPU path.
 * Claude Generated (Sep 2026): the wrapper logic lives in the shared template
 * GFNFFGpuMethodImpl<Backend> (gfnff_gpu_method_impl.h).  This translation unit — the
 * g++-compiled host side of the ROCm plugin — provides the HIP traits' one runtime call
 * and the explicit instantiation, so all wrapper code is emitted here as before.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu rocm
 *   ./curcuma -sp mol.xyz -method gfnff -gpu auto   # GPU if available
 */

#ifdef USE_ROCM

// The project headers come first (as in the CUDA sibling): <hip/hip_runtime.h> defines
// __noinline__ as a macro, which corrupts the [[__gnu__::__noinline__]] attribute in
// libstdc++'s <format> (pulled in via global.h -> <chrono>).
#include "gfnff_hip_method.h"

#include <hip/hip_runtime.h>  // host-side hipMemcpy for the GPU-Schur charge download

// ---------------------------------------------------------------------------
// ROCm backend traits: the single runtime-API call the wrapper needs
// ---------------------------------------------------------------------------

void GFNFFRocmBackend::downloadDoubles(double* host, const double* device, int n)
{
    hipMemcpy(host, device, n * sizeof(double), hipMemcpyDeviceToHost);
}

// ---------------------------------------------------------------------------
// Explicit instantiation of the shared wrapper for the ROCm backend
// ---------------------------------------------------------------------------

template class GFNFFGpuMethodImpl<GFNFFRocmBackend>;

#endif // USE_ROCM
