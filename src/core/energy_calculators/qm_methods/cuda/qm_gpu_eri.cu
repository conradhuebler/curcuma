/*
 * < CUDA 4-centre ERI build of the native QM engine (Q-G2) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, Q-G2 of docs/QM_GPU_ROADMAP.md).
//
// Builds the stored ERI tensor (mu nu | lam sig) of the native QM engine directly in device
// memory and hands it to the Q-G1 J/K builder (QmGpuJkCuda), so the tensor never crosses
// the PCIe bus. Reference: qmint::buildERI in qm_integrals.cpp; the host flattens the very
// same shell-pair tables (qmint::buildEriGpuTables), so shells, primitive-pair order,
// Hermite coefficients and both Schwarz screens are identical to the CPU.
// Since Q-G3 (Oct 2026) the kernels themselves live in qm_gpu_eri_kernels.cuh (shared with
// the integral-direct path qm_gpu_direct.cu); this file instantiates them with DIRECT = false.
//
// McMurchie-Davidson (Helgaker, Jorgensen, Olsen, ch. 9.9; McMurchie, Davidson,
// J. Comput. Phys. 26, 218 (1978)):
//   (ab|cd) = sum_prim c_a c_b c_c c_d 2 pi^(5/2) / (p q sqrt(p+q))
//             sum_{tuv} E^ab_t E^ab_u E^ab_v sum_{tau nu phi} (-1)^(tau+nu+phi)
//             E^cd_tau E^cd_nu E^cd_phi R^0_{t+tau, u+nu, v+phi}(rho, P-Q)
//
// Work decomposition (deliberately simple, correctness first):
//   * one thread block per canonical, Schwarz-surviving shell quartet; one launch per
//     angular-momentum class (LA,LB,LC,LD), so the dynamic shared memory is sized per class;
//   * per primitive quartet: thread 0 evaluates the Boys function (same grid + Taylor table
//     as the CPU), the block builds R^0 on the simplex t+u+v <= L level by level in shared
//     memory, then every thread accumulates its own cartesian component quartets in FP64
//     (each component owned by exactly one thread: no atomics, deterministic);
//   * the cartesian -> spherical 6d->5d transform is applied per quartet in shared memory,
//     one axis at a time (as transformAxis on the CPU);
//   * every tensor element is written by exactly one component of exactly one quartet
//     (canonical representative of its 8-fold orbit), then copied to all 8 permutations.
// Limits: s, p, d shells (L <= 2) in these kernels; classes with an f or g shell (stage 6b, Nov 2026) go to
// the generic kernel of qm_gpu_eri_generic.cu (one stream, shared scratch); FP64 throughout.

#include "qm_gpu_jk_cuda.h"
#include "qm_gpu_eri_kernels.cuh"  // kernels shared with the Q-G3 integral-direct path
#include "qm_gpu_eri_generic.cuh"  // stage 6b: classes with f/g shells

#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>


// Build the ERI tensor on the device and return a J/K builder that owns it (QMGpuJK*), or
// nullptr with a reason in `why` (unsupported basis, out of device memory, CUDA error).
extern "C" void* curcuma_cuda_create_qm_eri_jk_spin(const qmint::QMGpuEriTables* t, int device, char* why, int whylen)
{
    using clk = std::chrono::steady_clock;
    auto ms = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const auto t0 = clk::now();
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return nullptr; }
    cudaFree(nullptr);  // context creation, timed separately
    {
        // lazy module loading (CUDA_MODULE_LOADING=LAZY): load both kernels now, so the
        // kernel timing below is kernel work only
        cudaFuncAttributes a{};
        cudaFuncGetAttributes(&a, eriClassKernel<false>);
        qmgpu::warmJkKernel();
    }
    const auto t1 = clk::now();

    const size_t n = (size_t)t->n_out;
    const size_t tensor_bytes = n * n * n * n * sizeof(double);
    const size_t table_bytes = (size_t)t->n_E * 8 + (size_t)t->nquartet * 8 + (size_t)t->nprimpair * 64 + (1u << 20);
    size_t free_b = 0, total_b = 0;
    cudaMemGetInfo(&free_b, &total_b);
    // -qm.gpu_max_memory_mb (Q-G3): pretend the card is smaller, to force/test the direct path
    if (t->device_mem_limit_mb > 0.0 && (double)free_b > t->device_mem_limit_mb * 1.0e6)
        free_b = (size_t)(t->device_mem_limit_mb * 1.0e6);
    const size_t margin = 256ull << 20;
    if (tensor_bytes + table_bytes + 3 * n * n * 8 + margin > free_b) {
        char buf[256];
        std::snprintf(buf, sizeof buf, "ERI tensor %.1f MB does not fit into free device memory %.1f MB",
                      tensor_bytes / 1.0e6, free_b / 1.0e6);
        setWhy(why, whylen, buf);
        return nullptr;
    }

    std::vector<void*> owned;
    auto cleanup = [&owned]() { for (void* p : owned) cudaFree(p); owned.clear(); };
    DevTables d{};
    std::string err;
    bool ok = uploadCommonTables(t, d, owned, err, 4)
        && upload(const_cast<int**>(&d.q_bra), t->q_bra, (size_t)t->nquartet, owned)
        && upload(const_cast<int**>(&d.q_ket), t->q_ket, (size_t)t->nquartet, owned);
    // stage 6b (Claude Generated, Nov 2026): f/g shells -> F tables and scratch of the generic kernel
    qmgpu::GenTables gen{};
    const bool has_fg = qmgpu::maxShellL(t) > 2;
    if (ok && has_fg) ok = qmgpu::genInit(t, d, gen, owned, err);
    if (!ok) { cleanup(); setWhy(why, whylen, err.empty() ? "upload of the ERI tables failed" : err); return nullptr; }

    double* dE = nullptr;
    if (cudaMalloc(&dE, tensor_bytes) != cudaSuccess) {
        cleanup();
        setWhy(why, whylen, "cudaMalloc of the ERI tensor failed");
        return nullptr;
    }
    cudaMemset(dE, 0, tensor_bytes);
    cudaDeviceSynchronize();
    const auto t2 = clk::now();

    cudaEvent_t e0, e1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    // CURCUMA_QM_GPU_ERI_PROFILE=1: per-class kernel times on stderr (serialises the launches)
    const char* penv = std::getenv("CURCUMA_QM_GPU_ERI_PROFILE");
    const bool prof = penv && penv[0] == '1';
    cudaEvent_t p0, p1;
    cudaEventCreate(&p0);
    cudaEventCreate(&p1);
    bool spherical = false;
    const bool order_ok = threadKernelUsable(t, spherical);
    cudaEventRecord(e0);
    cudaEventRecord(p0);
    // Classes are independent (disjoint tensor elements), so they run on several streams:
    // for small molecules the per-class latency, not the work, sets the time.
    constexpr int kStreams = 8;
    cudaStream_t streams[kStreams];
    for (int k = 0; k < kStreams; ++k) cudaStreamCreateWithFlags(&streams[k], cudaStreamNonBlocking);
    const char* tenv = std::getenv("CURCUMA_QM_GPU_ERI_BLOCK_ONLY");  // diagnostic: block kernel everywhere
    const bool thread_ok = order_ok && !(tenv && tenv[0] == '1');
    int nthread_classes = 0, ngen_classes = 0;
    bool gen_failed = false;
    for (int c = 0; c < t->nclass; ++c) {
        const long long q0 = t->class_off[c], nq = t->class_off[c + 1] - q0;
        if (nq <= 0) continue;
        const int LA = t->class_L[4 * c], LB = t->class_L[4 * c + 1], LC = t->class_L[4 * c + 2], LD = t->class_L[4 * c + 3];
        const int ncart = ncompOf(LA) * ncompOf(LB) * ncompOf(LC) * ncompOf(LD);
        const bool fg = LA > 2 || LB > 2 || LC > 2 || LD > 2;
        // f/g classes share the generic kernel's scratch: always stream 0
        cudaStream_t st = (prof || fg) ? streams[0] : streams[c % kStreams];
        const QSrc src{ q0, nullptr, nullptr, nullptr, nullptr, 0 };
        bool threaded = false;
        if (fg) {
            if (!qmgpu::genLaunch(LA, LB, LC, LD, d, gen, src, nq, false, dE, st)) { gen_failed = true; break; }
            ++ngen_classes;
        } else {
            threaded = thread_ok && launchThreadClass<false>(LA, LB, LC, LD, d, src, nq, dE, st, spherical);
            if (threaded)
                ++nthread_classes;
            else
                launchBlockClass<false>(LA, LB, LC, LD, d, src, nq, dE, st);
        }
        if (cudaGetLastError() != cudaSuccess) break;
        if (prof) {
            cudaEventRecord(p1, st);
            cudaEventSynchronize(p1);
            float pm = 0.0f;
            cudaEventElapsedTime(&pm, p0, p1);
            std::fprintf(stderr, "qm_gpu_eri class %d%d%d%d quartets %lld ncart %d %s: %.3f ms\n", LA, LB, LC, LD,
                         nq, ncart, fg ? "generic" : threaded ? "thread" : "block", (double)pm);
            cudaEventRecord(p0, st);
        }
    }
    for (int k = 0; k < kStreams; ++k) {
        cudaStreamSynchronize(streams[k]);
        cudaStreamDestroy(streams[k]);
    }
    cudaEventRecord(e1);
    cudaError_t kerr = cudaEventSynchronize(e1) == cudaSuccess ? cudaGetLastError() : cudaErrorUnknown;
    if (kerr == cudaSuccess && gen_failed) kerr = cudaErrorLaunchFailure;
    float kms = 0.0f;
    cudaEventElapsedTime(&kms, e0, e1);
    cudaEventDestroy(e0);
    cudaEventDestroy(e1);
    cudaEventDestroy(p0);
    cudaEventDestroy(p1);
    cleanup();
    if (kerr != cudaSuccess) {
        cudaFree(dE);
        setWhy(why, whylen, std::string("ERI kernel failed: ") + cudaGetErrorString(kerr));
        return nullptr;
    }
    const auto t3 = clk::now();

    auto* jk = new qmgpu::QmGpuJkCuda(dE, (int)n, device, true);
    if (!jk->ok()) {
        delete jk;
        setWhy(why, whylen, "allocation of the J/K work buffers failed");
        return nullptr;
    }
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "context+module load %.1f ms, tables+tensor alloc %.1f ms, ERI kernels %.1f ms (%lld quartets, %d classes, %d thread-per-quartet, %d f/g generic), total %.1f ms",
                  ms(t0, t1), ms(t1, t2), (double)kms, t->nquartet, t->nclass, nthread_classes, ngen_classes, ms(t0, t3));
    jk->setStats(buf);
    return static_cast<qmint::QMGpuJK*>(jk);
}
