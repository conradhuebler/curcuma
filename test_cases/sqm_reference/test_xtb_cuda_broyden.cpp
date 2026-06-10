/*
 * Stage 6 (S6.4) component test: device modified-Broyden mixer vs the host
 * BroydenMixer (Johnson 1988).
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 6).
 *
 * Standalone numeric kernel test — no molecule. Drives a deterministic contractive
 * fixed-point map G(x) = x* + J·(x − x*) with ‖J‖_inf = 0.6 (so the iteration
 * converges and the residual genuinely changes direction across steps, exercising
 * the Broyden history). At each step BOTH mixers are fed the SAME (vin, vout=G(vin)):
 * the host BroydenMixer::update (Eigen .inverse() of the w0²I+a Gram system) and the
 * device XtbGpuContext::broydenBegin/broydenUpdate (cuBLAS Gram + a single-block
 * M×M regularised solve, FIFO history ring). The next input is advanced by the host
 * result so the inputs stay identical, and the configs run past max_hist to exercise
 * the history eviction.
 *
 * The device keeps history in a ring buffer, the host in an erase-front FIFO — both
 * retain the same most-recent-max_hist SET, and vnext = vin + α·F − Σ_k γ_k u_k is
 * permutation-invariant over that set, so they agree to round-off (gate 1e-10) even
 * after eviction.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_broyden [--quiet] [--tol V]
 */

#include "src/core/energy_calculators/qm_methods/broyden_mixer.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using curcuma::xtb::gpu::XtbGpuContext;

namespace {

// Deterministic contractive map G(x) = x* + J·(x − x*), ‖J‖_inf scaled to 0.6.
struct ContractiveMap {
    Eigen::MatrixXd J;
    Vector xstar;
    explicit ContractiveMap(int N)
    {
        J.resize(N, N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                J(i, j) = std::sin(1.0 + 1.7 * i - 0.9 * j);   // deterministic in [−1,1]
        double maxrow = 0.0;
        for (int i = 0; i < N; ++i) {
            double row = 0.0;
            for (int j = 0; j < N; ++j) row += std::fabs(J(i, j));
            maxrow = std::max(maxrow, row);
        }
        if (maxrow > 0.0) J *= 0.6 / maxrow;                   // ‖J‖_inf = 0.6 < 1
        xstar.resize(N);
        for (int i = 0; i < N; ++i) xstar(i) = 0.2 * std::cos(0.5 * i + 0.3);
    }
    Vector operator()(const Vector& x) const { return xstar + J * (x - xstar); }
    Vector start() const
    {
        Vector x0(xstar.size());
        for (int i = 0; i < xstar.size(); ++i) x0(i) = xstar(i) + 0.7 * std::sin(0.9 * i + 0.1);
        return x0;
    }
};

// Run one (N, max_hist) config; return the max |host − device| over all steps,
// or a negative value on a device failure.
double runConfig(XtbGpuContext& ctx, int N, int max_hist, int nsteps, double alpha, double w0)
{
    BroydenMixer host(alpha, max_hist, w0);
    if (!ctx.broydenBegin(N, alpha, max_hist, w0)) return -1.0;

    ContractiveMap G(N);
    Vector vin = G.start();
    std::vector<double> vnext_gpu(N, 0.0);
    double dmax = 0.0;
    for (int step = 0; step < nsteps; ++step) {
        const Vector vout = G(vin);
        const Vector host_next = host.update(vin, vout);
        if (!ctx.broydenUpdate(N, vin.data(), vout.data(), vnext_gpu.data())) return -1.0;
        for (int i = 0; i < N; ++i)
            dmax = std::max(dmax, std::fabs(host_next(i) - vnext_gpu[i]));
        vin = host_next;   // feed both mixers the same next input
    }
    return dmax;
}

} // namespace

int main(int argc, char* argv[])
{
    CurcumaLogger::set_verbosity(0);

    double tol = 1.0e-10;
    bool quiet = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--quiet") quiet = true;
        else if (a == "--tol" && i + 1 < argc) tol = std::stod(argv[++i]);
    }

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_broyden: no usable CUDA device\n";
        return 0;
    }

    // (N, max_hist): below / at the cap, and a case that overflows the ring.
    const struct { int N, max_hist; } configs[] = {{12, 6}, {20, 20}, {31, 8}};
    const int nsteps = 14;
    const double alpha = 0.25, w0 = 0.01;

    int failures = 0;
    for (const auto& c : configs) {
        const double dmax = runConfig(ctx, c.N, c.max_hist, nsteps, alpha, w0);
        if (dmax < 0.0) {
            std::cerr << "broyden config N=" << c.N << " max_hist=" << c.max_hist
                      << ": device update failed\n";
            ++failures; continue;
        }
        const bool ok = (dmax <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s broyden N=%-3d max_hist=%-2d steps=%d: max|d_vnext|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", c.N, c.max_hist, nsteps, dmax, tol);
    }
    return failures == 0 ? 0 : 1;
}
