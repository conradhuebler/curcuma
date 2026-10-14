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
 * Preconditioned configs (Claude Generated, Oct 2026, -scf_precond frag on the device): the
 * host mixer gets the step hook M*v of a SchwarzPreconditioner (non-overlapping random blocks,
 * some indices left on the scalar alpha), the device the flattened M
 * (broydenSetPreconditioner). Guard configs drive the host as XTB::Calculation does
 * (ScfPrecondGuard::fires, BroydenMixer::dropStepHook, restart from x_best/F_best) and the device
 * through the callback on a second guard copy (firesNorm); rule g2:3:0 always fires at step 3,
 * so keep / reset / best are all exercised.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_broyden [--quiet] [--tol V]
 */

#include "src/core/energy_calculators/qm_methods/broyden_mixer.h"
#include "src/core/energy_calculators/qm_methods/xtb_precond.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
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

// Claude Generated (Oct 2026): block-diagonal M over a deterministic permutation of the indices;
// the last few indices stay uncovered (scalar alpha). Blocks are diagonally dominant around
// alpha so the preconditioned iteration still converges.
curcuma::xtb::SchwarzPreconditioner makePreconditioner(int N, double alpha)
{
    curcuma::xtb::SchwarzPreconditioner P(N, alpha);
    std::vector<int> perm(N);
    for (int i = 0; i < N; ++i) perm[i] = (7 * i + 3) % N;   // N not divisible by 7 below
    const int sizes[] = { 5, 1, 7, 4 };
    int pos = 0;
    for (int b = 0; pos < N - 3; ++b) {
        const int n = std::min(sizes[b % 4], N - 3 - pos);
        curcuma::xtb::SchwarzPreconditioner::Block blk;
        blk.B.resize(n, n);
        for (int r = 0; r < n; ++r) {
            blk.index.push_back(perm[pos + r]);
            for (int c = 0; c < n; ++c)
                blk.B(r, c) = (r == c ? 0.6 : 0.0) + 0.08 * std::sin(1.3 + 0.7 * r - 1.1 * c + b);
        }
        P.addBlock(std::move(blk));
        pos += n;
    }
    return P;
}

// Preconditioned run, optionally with a guard spec; host = reference semantics of XTB::Calculation.
double runPrecondConfig(XtbGpuContext& ctx, int N, int max_hist, int nsteps, double alpha, double w0,
                        const std::string& guard_spec, unsigned long long tag, int& fired_dev)
{
    using curcuma::xtb::ScfPrecondGuard;
    const auto P = makePreconditioner(N, alpha);
    BroydenMixer host(alpha, max_hist, w0);
    host.setStepHook([&P](const Vector& v) -> Vector { return P.apply(v); });
    ScfPrecondGuard guard_host = ScfPrecondGuard::parse(guard_spec, "test", true);
    auto guard_dev = std::make_shared<ScfPrecondGuard>(guard_host);

    if (!ctx.broydenBegin(N, alpha, max_hist, w0)) return -1.0;
    std::vector<int> ioff, idx, boff;
    std::vector<double> B, fb;
    int maxdim = 0;
    if (!P.flatten(ioff, idx, boff, B, fb, maxdim)) return -1.0;
    XtbGpuContext::BroydenPrecond pc;
    pc.N = N; pc.nblocks = static_cast<int>(ioff.size()) - 1; pc.max_dim = maxdim;
    pc.idx_off = ioff.data(); pc.idx = idx.data(); pc.b_off = boff.data();
    pc.B = B.data(); pc.fallback = fb.data();
    XtbGpuContext::BroydenGuard cb;
    if (guard_host.rule != ScfPrecondGuard::Rule::Off)
        cb = [guard_dev](int n, double f, bool* store) { return guard_dev->firesNorm(n, f, store); };
    if (!ctx.broydenSetPreconditioner(&pc, tag, cb, guard_host.keep, guard_host.best)) return -1.0;

    ContractiveMap G(N);
    Vector vin = G.start();
    std::vector<double> vnext_gpu(N, 0.0);
    double dmax = 0.0;
    for (int step = 0; step < nsteps; ++step) {
        const Vector vout = G(vin);
        bool restart_best = false;
        if (guard_host.rule != ScfPrecondGuard::Rule::Off && host.hasStepHook()
            && guard_host.fires(step, vin, vout - vin)) {
            host.dropStepHook(guard_host.keep);
            restart_best = guard_host.best;
        }
        const Vector host_next = restart_best
            ? host.update(guard_host.x_best, guard_host.x_best + guard_host.F_best)
            : host.update(vin, vout);
        if (!ctx.broydenUpdate(N, vin.data(), vout.data(), vnext_gpu.data())) return -1.0;
        for (int i = 0; i < N; ++i)
            dmax = std::max(dmax, std::fabs(host_next(i) - vnext_gpu[i]));
        vin = host_next;
    }
    fired_dev = ctx.broydenGuardFiredAt();
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

    // Claude Generated (Oct 2026): preconditioned device Broyden (-scf_precond frag) vs the host
    // mixer with the step hook, without and with the guard (fires at step 3 by construction).
    const struct { int N, max_hist; const char* guard; int fired; } pconfigs[] = {
        { 31, 8, "off", -1 }, { 40, 20, "off", -1 }, { 31, 6, "g2:3:0", 3 },
        { 31, 6, "g2:3:0:keep", 3 }, { 31, 6, "g2:3:0:best", 3 }, { 40, 8, "g2:3:0:keep:best", 3 },
        { 31, 8, "g2:4", -2 } };   // -2: firing not prescribed (contractive map), only agreement
    unsigned long long tag = 1;
    for (const auto& c : pconfigs) {
        int fired = -1;
        const double dmax = runPrecondConfig(ctx, c.N, c.max_hist, nsteps, alpha, w0, c.guard, tag++, fired);
        const bool fired_ok = (c.fired == -2) || (fired == c.fired);
        const bool ok = dmax >= 0.0 && dmax <= tol && fired_ok;
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s broyden+M N=%-3d max_hist=%-2d guard=%-16s steps=%d: max|d_vnext|=%.3e tol=%.1e fired=%d\n",
                        ok ? "PASS" : "FAIL", c.N, c.max_hist, c.guard, nsteps, dmax, tol, fired);
    }
    // Reuse of the resident upload (pc = nullptr, same tag) after a new broydenBegin; a wrong tag fails.
    {
        const bool begin_ok = ctx.broydenBegin(31, alpha, 8, w0);
        const unsigned long long t = ctx.broydenPreconditionerTag();
        const bool reuse_ok = begin_ok && ctx.broydenSetPreconditioner(nullptr, t, nullptr, false, false)
                              && ctx.broydenPreconditionerActive();
        const bool wrong_ok = ctx.broydenBegin(31, alpha, 8, w0)
                              && !ctx.broydenSetPreconditioner(nullptr, t + 1000, nullptr, false, false)
                              && !ctx.broydenPreconditionerActive();
        if (!reuse_ok || !wrong_ok) ++failures;
        if (!quiet || !reuse_ok || !wrong_ok)
            std::printf("%s broyden+M resident reuse by tag (reuse %d, wrong tag rejected %d)\n",
                        (reuse_ok && wrong_ok) ? "PASS" : "FAIL", reuse_ok, wrong_ok);
    }
    return failures == 0 ? 0 : 1;
}
