/*
 * <Device Ewald sum of the GFN-FF Gaussian charges (CUDA)>
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
 * Claude Generated (Nov 2026, GPU_PERIODIC_GAPS A2, stages S3a/S3b of docs/WP-PERIODIC-NONBONDED.md).
 * Status: AI-generated, machine-tested (ctest gfnff_gpu_pbc_ewald), no human testing.
 *
 * GfnffEwaldGPU implements curcuma::ewald::EwaldDeviceOps (ewald.h) for the CUDA periodic hybrid of GFN-FF:
 * the real-space and reciprocal parts of the Ewald sum of Gaussian charges, the same quantities as
 * EwaldSum::realSpaceRangeCached / reciprocalRange / realSpaceVirialRange / reciprocalVirialRange on the CPU.
 *
 *  - Pair list: per atom i the full list of candidates j != i in ascending j (one warp per atom, FP32 test with a
 *    0.05 Bohr margin on r_cut), then one thread per entry: the exact FP64 minimum-image test r^2 < r_cut^2 and the
 *    kernel K_ij = erfc(a r)/r - erfc(g_ij r)/r (g_ij r < 6) with dK/dr (EwaldSum pairKernel); entries beyond r_cut
 *    keep K = dK = 0, which sums to the same as the CPU list without them.
 *  - Phase table: cos/sin(2 pi n s_a), s_a = fractional coordinate wrapped to [0, 1), n = 0..nmax_a (sincospi);
 *    the CPU builds the same numbers by a recursion, so the two agree to rounding, not bitwise.
 *  - Reciprocal: structure factors S(k) for groups of up to 8 k vectors with the same (n_x, n_y) (one block per group,
 *    fixed-order block tree), then one thread per atom and slice of groups (potentials or gradient), the slices summed
 *    in order. The algebra is arranged so that each (k, atom) costs about four FP64 multiply-adds (FP64 is 1/64 of
 *    FP32 on the RTX 5080, and these passes are FP64-bound).
 *  - Every sum has a fixed order (warp shuffles, block trees, one-block final reductions): run to run bitwise.
 * Self, background and Gaussian background stay on the host (EwaldSum::selfAndBackground).
 * Fault injection for the fallback tests: CURCUMA_GFNFF_EWALD_DEVICE_FAIL=prepare|potentials|evaluate.
 */

#pragma once

#include "../ewald.h"

#include <memory>
#include <string>

class GfnffEwaldGPU : public curcuma::ewald::EwaldDeviceOps {
public:
    GfnffEwaldGPU();
    ~GfnffEwaldGPU() override;
    GfnffEwaldGPU(const GfnffEwaldGPU&) = delete;
    GfnffEwaldGPU& operator=(const GfnffEwaldGPU&) = delete;

    bool prepare(const curcuma::ewald::EwaldSum& ew, const Eigen::MatrixXd& positions, const Eigen::VectorXd& alp) override;
    bool potentials(const Eigen::VectorXd& x, Eigen::VectorXd& y) override;
    bool evaluate(const Eigen::VectorXd& q, curcuma::ewald::EwaldEnergyTerms& terms, Eigen::MatrixXd* grad,
        Eigen::Matrix3d* virial) override;
    std::string lastError() const override;
    double lastPrepareMs() const override;

    /// Profile counters (step diagnostics, tests)
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

private:
    struct Impl;
    std::unique_ptr<Impl> m;
    Stats m_stats;
};
