/*
 * Multi-step EEQ PCG warm-start extrapolation validation (GFN-FF).
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Drives GFN-FF over a smooth synthetic trajectory of a multi-fragment system
 * (so the iterative PCG EEQ solver is used) and compares the default 1-step PCG
 * warm-start (eeq_extrapolation=none) with the opt-in multi-step schemes
 * (aspc / gauss). Because PCG converges to the same charges regardless of the
 * initial guess, the energies must be identical to the 'none' run — the
 * extrapolation only changes the PCG iteration count, never the result.
 *
 * Asserts, for each scheme:
 *   (A) Correctness — every-frame energy equals the 'none' run within a tight gate.
 *   (B) Activation  — the extrapolation path actually fired (count > 0 for
 *       aspc/gauss, == 0 for none), via GFNFFComputationalMethod::eeqPcgExtrapolationCount().
 *
 *   test_eeq_extrapolation <molecule.xyz> [nsteps]
 *
 * Claude Generated (EEQ extrapolation WP, June 2026). GPL-3.0.
 */

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/global.h"
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"

#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

struct Result {
    std::vector<double> energies;
    long extrap_count = -1;
    bool ok = false;
};

// Run GFN-FF (forced PCG EEQ) over a smooth synthetic trajectory with the given
// extrapolation scheme and record the per-frame energy + extrapolation count.
Result run(const Mol& base, int nsteps, const std::string& mode)
{
    Result r;
    json config = {
        { "verbosity", 0 },
        { "threads", 1 },
        { "gfnff", json::object() },
        { "eeq_solver", {
            { "solve_method", "pcg" },
            { "eeq_extrapolation", mode },
            { "eeq_extrapolation_order", 3 }
        }}
    };

    EnergyCalculator ec("gfnff", config);

    Mol mol = base;
    const int nat = mol.m_number_atoms;
    const Geometry geom0 = mol.m_geometry;

    ec.setMolecule(mol);

    for (int k = 0; k <= nsteps; ++k) {
        Matrix geom = geom0;
        // Small, smooth, deterministic per-atom displacement (Angstrom). Same for
        // every scheme, so the geometry sequence is identical.
        if (k > 0) {
            for (int i = 0; i < nat; ++i)
                for (int d = 0; d < 3; ++d)
                    geom(i, d) = geom0(i, d) + 0.01 * std::sin(0.30 * k + 0.5 * i + d);
            ec.updateGeometry(geom);
        }
        const double e = ec.CalculateEnergy(false);
        r.energies.push_back(e);
        if (!std::isfinite(e)) return r;
    }

    auto* m = dynamic_cast<GFNFFComputationalMethod*>(ec.Interface());
    r.extrap_count = m ? m->eeqPcgExtrapolationCount() : -1;
    r.ok = true;
    return r;
}

bool check(const Mol& base, int nsteps, const std::string& mode, const Result& none)
{
    const double E_GATE = 1.0e-6; // PCG converges to the same charges -> same energy
    Result run_m = run(base, nsteps, mode);
    if (!run_m.ok || run_m.energies.size() != none.energies.size()) {
        std::printf("[%s] FAIL: run did not finish / size mismatch\n", mode.c_str());
        return false;
    }

    double max_dE = 0.0;
    for (size_t k = 0; k < none.energies.size(); ++k)
        max_dE = std::max(max_dE, std::fabs(run_m.energies[k] - none.energies[k]));

    std::printf("[%s] max|dE vs none| over %zu frames = %.2e ; extrapolation applied %ld times\n",
                mode.c_str(), none.energies.size(), max_dE, run_m.extrap_count);

    bool ok = true;
    if (max_dE > E_GATE) {
        std::printf("[%s] FAIL: energy differs from none by %.2e > %.0e\n", mode.c_str(), max_dE, E_GATE);
        ok = false;
    }
    if (run_m.extrap_count <= 0) {
        std::printf("[%s] FAIL: extrapolation never fired (count=%ld) -- PCG path not exercised?\n",
                    mode.c_str(), run_m.extrap_count);
        ok = false;
    }
    if (ok)
        std::printf("[%s] PASS\n", mode.c_str());
    return ok;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <molecule.xyz> [nsteps]\n", argv[0]);
        return 2;
    }
    const std::string xyz = argv[1];
    const int nsteps = (argc > 2) ? std::atoi(argv[2]) : 12;

    curcuma::Molecule molecule(xyz);
    Mol base = molecule.getMolInfo();
    if (base.m_number_atoms <= 0) {
        std::fprintf(stderr, "failed to load %s\n", xyz.c_str());
        return 2;
    }

    Result none = run(base, nsteps, "none");
    if (!none.ok) {
        std::printf("[none] FAIL: baseline run did not finish\n");
        std::printf("FAILURES\n");
        return 1;
    }
    std::printf("[none] %zu frames, final E=%.10f, extrapolation applied %ld times (expect 0)\n",
                none.energies.size(), none.energies.back(), none.extrap_count);

    bool ok = (none.extrap_count == 0);
    if (!ok)
        std::printf("[none] FAIL: extrapolation fired with mode=none (count=%ld)\n", none.extrap_count);

    ok &= check(base, nsteps, "aspc", none);
    ok &= check(base, nsteps, "gauss", none);

    std::printf("%s\n", ok ? "ALL PASS" : "FAILURES");
    return ok ? 0 : 1;
}
