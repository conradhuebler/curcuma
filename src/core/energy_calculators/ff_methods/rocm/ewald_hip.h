/*
 * <Device Ewald sum of the GFN-FF Gaussian charges (ROCm/HIP)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (Nov 24, 2026, ROCm parity B1, docs/ROCM_PARITY_B1.md). 🤖 AI-generated, ⚙️ machine-tested on one
 * RX 9060 XT (gfx1200); human production testing pending.
 *
 * GfnffEwaldHip implements curcuma::ewald::EwaldDeviceOps (ewald.h) for the ROCm periodic hybrid of GFN-FF. It is the
 * HIP port of GfnffEwaldGPU (cuda/ewald_gpu.h/.cu, S3 of docs/WP-PERIODIC-NONBONDED.md): same pair list (full list per
 * atom, ascending j, FP32 pre-screen with a 0.05 Bohr margin, exact FP64 minimum-image test per entry), same kernel
 * K_ij = erfc(a r)/r - erfc(g_ij r)/r with dK/dr, same phase table (sincospi), same structure-factor groups of 8 k
 * vectors and per-atom slices, same fixed-order sums (no atomics). The "warp" of the CUDA kernels is a group of 32
 * lanes here (shuffles with width 32, ballot shifted to the lane group), so the kernels also run on wave64 hardware.
 * Self, background and Gaussian background stay on the host (EwaldSum::selfAndBackground).
 * Fault injection: CURCUMA_GFNFF_EWALD_DEVICE_FAIL=prepare|potentials|evaluate; profile CURCUMA_GFNFF_EWALD_PROFILE=1.
 *
 * No HIP type appears here: the g++-compiled wrapper (gfnff_hip_method.cpp) includes this header; the definitions are
 * in rocm/ewald_hip.hip (own hipcc translation unit).
 */

#pragma once

#include "../ewald.h"

#include <memory>
#include <string>

class GfnffEwaldHip : public curcuma::ewald::EwaldDeviceOps {
public:
    GfnffEwaldHip();
    ~GfnffEwaldHip() override;
    GfnffEwaldHip(const GfnffEwaldHip&) = delete;
    GfnffEwaldHip& operator=(const GfnffEwaldHip&) = delete;

    bool prepare(const curcuma::ewald::EwaldSum& ew, const Eigen::MatrixXd& positions, const Eigen::VectorXd& alp) override;
    bool potentials(const Eigen::VectorXd& x, Eigen::VectorXd& y) override;
    bool evaluate(const Eigen::VectorXd& q, curcuma::ewald::EwaldEnergyTerms& terms, Eigen::MatrixXd* grad,
        Eigen::Matrix3d* virial) override;
    std::string lastError() const override;
    double lastPrepareMs() const override;

    /// Profile counters (step diagnostics, tests); same fields as GfnffEwaldGPU::Stats
    struct Stats {
        long rebuilds = 0;          ///< prepare() calls that rebuilt the device data
        long key_hits = 0;          ///< prepare() calls that found the data current
        long potentials_calls = 0;
        long evaluate_calls = 0;
        long long full_pairs = 0;   ///< entries of the full neighbour list (2 x pairs) of the last build
        double potentials_ms = 0.0; ///< accumulated wall time of potentials()
        double evaluate_ms = 0.0;   ///< accumulated wall time of evaluate()
    };
    const Stats& stats() const { return m_stats; }

    /// True when the constructor got a HIP stream (a device is present)
    bool usable() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m;
    Stats m_stats;
};
