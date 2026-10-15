/*
 * < CUDA integral-direct J/K of the native QM engine (Q-G3) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, Q-G3 of docs/QM_GPU_ROADMAP.md).
//
// Integral-direct Coulomb/exchange build on the device: no stored ERI tensor, device memory
// O(n^2) (P, Pmax, J, K) plus the O(npair) shell-pair tables. Device counterpart of
// qmint::DirectJK::build (qm_integrals.cpp; Almloef, Faegri, Korsell, J. Comput. Chem. 3, 385
// (1982); Haeser, Ahlrichs, J. Comput. Chem. 10, 104 (1989)):
//   1. Pmax(A,B) = max |P| over the shell block (A,B)                      (kernel pmaxKernel)
//   2. every canonical shell quartet (i <= j) with
//        Q_i Q_j max(Pmax over AB, CD, AC, AD, BC, BD) >= thr
//      is computed with the Q-G2 kernels (qm_gpu_eri_kernels.cuh, DIRECT = true) and digested
//      at once into the unsymmetrised Jt, Kt (FP64 atomicAdd, weight deg/8);
//   3. J = (Jt + Jt^T)/2, K = (Kt + Kt^T)/2 on the host.
// Quartets are enumerated per pair class: the launch for quartet class (X, Y) runs over the
// rectangle (pairs of class X) x (pairs of class Y), the non-canonical half returns at once,
// so no O(npair^2) quartet list exists anywhere. The incremental-on-dP scheme and the
// full-rebuild cadence are the host SCF's (qm_scf.cpp), unchanged: this builder is linear in P.
// Atomic accumulation makes the last bits of J/K depend on the thread schedule (not
// run-to-run bit-identical).

#include "qm_gpu_eri_kernels.cuh"

#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

/// Pmax(A,B) = Pmax(B,A) = max |P| over the shell block, for A <= B (as DirectJK).
__global__ void pmaxKernel(const double* __restrict__ P, const int* __restrict__ first, const int* __restrict__ width,
                           int ns, int n, double* __restrict__ Pmax)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= ns * ns) return;
    const int A = k / ns, B = k % ns;
    if (A > B) return;
    double m = 0.0;
    for (int a = 0; a < width[A]; ++a)
        for (int b = 0; b < width[B]; ++b) m = fmax(m, fabs(P[(size_t)(first[A] + a) * n + first[B] + b]));
    Pmax[A * ns + B] = m;
    Pmax[B * ns + A] = m;
}

class QmGpuDirectJkCuda final : public qmint::QMGpuJK {
public:
    QmGpuDirectJkCuda() = default;
    ~QmGpuDirectJkCuda() override
    {
        for (void* p : m_owned) cudaFree(p);
        for (auto& s : m_streams)
            if (s) cudaStreamDestroy(s);
        if (m_e0) cudaEventDestroy(m_e0);
        if (m_e1) cudaEventDestroy(m_e1);
    }

    bool init(const qmint::QMGpuEriTables* t, int device, std::string& why)
    {
        if (!t->pc_off || !t->pr_schwarz) { why = "tables without pair classes (built with list_quartets)"; return false; }
        if (!uploadCommonTables(t, m_d, m_owned, why)) return false;
        m_n = t->n_out;
        m_ns = t->nshell;
        m_d.thr = t->screening;
        m_spherical = false;
        const char* tenv = std::getenv("CURCUMA_QM_GPU_ERI_BLOCK_ONLY");  // diagnostic: block kernel everywhere
        m_thread_ok = threadKernelUsable(t, m_spherical) && !(tenv && tenv[0] == '1');
        m_pc_off.assign(t->pc_off, t->pc_off + 10);
        const size_t n2 = (size_t)m_n * m_n;
        bool ok = upload(const_cast<double**>(&m_d.pr_schwarz), t->pr_schwarz, (size_t)t->npair, m_owned)
            && upload(&m_pc_pairs, t->pc_pairs, (size_t)t->npair, m_owned)
            && upload(&m_width, t->sh_out_width, (size_t)t->nshell, m_owned)
            && alloc(&m_P, n2) && alloc(&m_J, n2) && alloc(&m_K, n2) && alloc(&m_Pmax, (size_t)m_ns * m_ns)
            && alloc(&m_counter, 1);
        if (!ok) { why = "device allocation / upload for the direct J/K failed"; return false; }
        m_d.P = m_P;
        m_d.Pmax = m_Pmax;
        m_d.J = m_J;
        m_d.K = m_K;
        m_d.counter = m_counter;
        m_first = const_cast<int*>(m_d.sh_out_first);
        for (auto& s : m_streams)
            if (cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking) != cudaSuccess) { why = "stream creation failed"; return false; }
        cudaEventCreate(&m_e0);
        cudaEventCreate(&m_e1);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        m_bytes = 0;
        return cudaGetLastError() == cudaSuccess;
    }

    bool build(const double* P, double* J, double* K) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n2 = (size_t)m_n * m_n;
        if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        cudaMemset(m_J, 0, n2 * sizeof(double));
        cudaMemset(m_K, 0, n2 * sizeof(double));
        cudaMemset(m_counter, 0, sizeof(unsigned long long));
        cudaEventRecord(m_e0);
        const int nsn = m_ns * m_ns;
        pmaxKernel<<<(nsn + 127) / 128, 128>>>(m_P, m_first, m_width, m_ns, m_n, m_Pmax);
        cudaEventRecord(m_e1);  // later launches on the streams wait for Pmax and the memsets
        for (auto& s : m_streams) cudaStreamWaitEvent(s, m_e1, 0);

        int launch = 0;
        for (int X = 0; X < 9; ++X) {
            const long long nX = m_pc_off[X + 1] - m_pc_off[X];
            if (nX == 0) continue;
            for (int Y = 0; Y < 9; ++Y) {
                const long long nY = m_pc_off[Y + 1] - m_pc_off[Y];
                if (nY == 0) continue;
                const int LA = X / 3, LB = X % 3, LC = Y / 3, LD = Y % 3;
                const QSrc src{ 0, m_pc_pairs + m_pc_off[X], m_pc_pairs + m_pc_off[Y], nY };
                const long long nq = nX * nY;
                cudaStream_t st = m_streams[launch++ % kStreams];
                const bool threaded = m_thread_ok && launchThreadClass<true>(LA, LB, LC, LD, m_d, src, nq, nullptr, st, m_spherical);
                if (!threaded) launchBlockClass<true>(LA, LB, LC, LD, m_d, src, nq, nullptr, st);
                if (cudaGetLastError() != cudaSuccess) return false;
            }
        }
        for (auto& s : m_streams)
            if (cudaStreamSynchronize(s) != cudaSuccess) return false;
        cudaEventRecord(m_e1);
        if (cudaEventSynchronize(m_e1) != cudaSuccess || cudaGetLastError() != cudaSuccess) return false;
        float kms = 0.0f;
        cudaEventElapsedTime(&kms, m_e0, m_e1);
        unsigned long long cnt = 0;
        if (cudaMemcpy(J, m_J, n2 * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
            || cudaMemcpy(K, m_K, n2 * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
            || cudaMemcpy(&cnt, m_counter, sizeof cnt, cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        // J = (Jt + Jt^T)/2, K = (Kt + Kt^T)/2, as DirectJK::build
        const int n = m_n;
        for (int a = 0; a < n; ++a)
            for (int b = a; b < n; ++b) {
                const double js = 0.5 * (J[(size_t)a * n + b] + J[(size_t)b * n + a]);
                const double ks = 0.5 * (K[(size_t)a * n + b] + K[(size_t)b * n + a]);
                J[(size_t)a * n + b] = J[(size_t)b * n + a] = js;
                K[(size_t)a * n + b] = K[(size_t)b * n + a] = ks;
            }
        m_last_quartets = (long long)cnt;
        ++m_builds;
        m_kernel_ms += kms;
        m_total_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    std::string name() const override { return m_name; }
    int n() const override { return m_n; }
    long long lastQuartets() const override { return m_last_quartets; }
    std::string stats() const override
    {
        char buf[256];
        std::snprintf(buf, sizeof buf, "%d direct J/K builds, kernels %.1f ms, incl. transfers %.1f ms, device memory %.1f MB",
                      m_builds, m_kernel_ms, m_total_ms, m_bytes / 1.0e6);
        return buf;
    }

private:
    template <typename T>
    bool alloc(T** p, size_t count)
    {
        if (cudaMalloc(p, count * sizeof(T)) != cudaSuccess) { *p = nullptr; return false; }
        m_owned.push_back(*p);
        return true;
    }

    static constexpr int kStreams = 8;
    DevTables m_d{};
    std::vector<void*> m_owned;
    std::vector<long long> m_pc_off;
    int* m_pc_pairs = nullptr;
    int* m_width = nullptr;
    int* m_first = nullptr;
    double *m_P = nullptr, *m_J = nullptr, *m_K = nullptr, *m_Pmax = nullptr;
    unsigned long long* m_counter = nullptr;
    cudaStream_t m_streams[kStreams] = {};
    cudaEvent_t m_e0 = nullptr, m_e1 = nullptr;
    int m_n = 0, m_ns = 0, m_builds = 0;
    bool m_spherical = false, m_thread_ok = false;
    long long m_last_quartets = 0;
    double m_kernel_ms = 0.0, m_total_ms = 0.0;
    size_t m_bytes = 0;
    std::string m_name = "cuda";

public:
    void setBytes(size_t b) { m_bytes = b; }
};

} // namespace

// Integral-direct J/K builder (QMGpuJK*), or nullptr with a reason in `why`. Device memory:
// the shell/pair tables and 4 n^2 + nshell^2 doubles; no n^4 tensor. Claude Generated (Q-G3).
extern "C" void* curcuma_cuda_create_qm_direct_jk(const qmint::QMGpuEriTables* t, int device, char* why, int whylen)
{
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return nullptr; }
    cudaFree(nullptr);  // context creation
    size_t free0 = 0, total = 0;
    cudaMemGetInfo(&free0, &total);
    auto b = std::make_unique<QmGpuDirectJkCuda>();
    std::string err;
    if (!b->init(t, device, err)) {
        cudaGetLastError();  // clear a sticky-free error state for the CPU fallback
        setWhy(why, whylen, err);
        return nullptr;
    }
    size_t free1 = 0;
    cudaMemGetInfo(&free1, &total);
    b->setBytes(free0 > free1 ? free0 - free1 : 0);
    return static_cast<qmint::QMGpuJK*>(b.release());
}
