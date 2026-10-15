/*
 * < CUDA J/K builder of the native QM engine: plugin-internal declaration >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, Q-G2): shared between qm_gpu_jk.cu (Q-G1, tensor uploaded
// from the host) and qm_gpu_eri.cu (Q-G2, tensor built on the device). Only compiled into
// libcurcuma_cuda.so.

#include "../qm_gpu_backend.h"

#include <string>

namespace qmgpu {

/// Force the (lazy) module load of the J/K kernel, so it is not timed as kernel work.
void warmJkKernel();

class QmGpuJkCuda final : public qmint::QMGpuJK {
public:
    /// Q-G1: allocate and upload the host tensor eri (n^4 doubles).
    QmGpuJkCuda(const double* eri, int n, int device);
    /// Q-G2: adopt a device tensor d_eri (n^4 doubles, cudaMalloc'ed); freed by this object.
    QmGpuJkCuda(double* d_eri, int n, int device, bool adopt);
    ~QmGpuJkCuda() override;
    bool ok() const { return m_E != nullptr; }

    bool build(const double* P, double* J, double* K) override;
    /// Q-G4: P, J, K on the device (resident SCF loop); synchronous.
    bool buildDevice(const double* dP, double* dJ, double* dK) override;
    std::string name() const override { return m_name; }
    int n() const override { return m_n; }
    bool downloadERI(double* host) const override;
    std::string stats() const override { return m_stats; }
    void setStats(const std::string& s) { m_stats = s; }

private:
    bool allocWork(int device);
    void release();
    int m_n;
    double *m_E = nullptr, *m_P = nullptr, *m_J = nullptr, *m_K = nullptr;
    std::string m_name = "cuda";
    std::string m_stats;
};

} // namespace qmgpu
