/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * GFN-FF energy-only call at a NEW geometry vs a gradient call at that geometry.
 *
 * Claude Generated (Sep 2026, MULTI_GPU_GAPS F-11). An energy-only call (optimizer line search,
 * energy-only scan) must give the same energy as a gradient call at the same geometry. The GPU
 * path refreshed the D4 C6 only in gradient calls, so an energy-only call after a geometry
 * change used the C6 of the previous geometry. The CPU path had the same class of defect in the
 * Coulomb self-energy: the workspace CN behind chi = chi_base + cnf*sqrt(CN) was set only by
 * gradient calls (triose: 4.6e-3 Eh). Both backends are covered ("none" = CPU).
 *
 * Both calculators get the identical history (two gradient calls at A, then geometry B), so the
 * topology perceived at A is kept in both; a FRESH calculator at B is not a valid reference,
 * because it perceives the topology at B (4.6e-3 Eh apart on triose, CPU and GPU alike).
 *
 * Usage: test_gfnff_energy_only_geometry <molecule.xyz> [gpu backend, default "cuda"] [verbosity]
 * Exit 0 when |E(energy-only at B) - E(gradient call at B)| <= 1e-8 Eh.
 */
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    const int verbosity = (argc > 3) ? std::atoi(argv[3]) : 0;
    CurcumaLogger::set_verbosity(verbosity);
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s molecule.xyz [cuda|none]\n", argv[0]);
        return 2;
    }
    const std::string gpu = (argc > 2) ? argv[2] : "cuda";
    json config = { { "verbosity", verbosity }, { "threads", 1 }, { "gfnff", json::object() }, { "gpu", gpu } };

    curcuma::Molecule molecule(argv[1]);
    Mol mol_a = molecule.getMolInfo();

    // Geometry B: every atom displaced by up to 0.08 Angstrom (fixed seed), enough to move the
    // coordination numbers and hence the D4 C6.
    Mol mol_b = mol_a;
    std::mt19937 rng(4711);
    std::uniform_real_distribution<double> d(-0.08, 0.08);
    for (int i = 0; i < mol_b.m_geometry.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            mol_b.m_geometry(i, k) += d(rng);

    // Same history for both: two gradient calls at A (as in an optimizer; the second one runs
    // the device EEQ path), then geometry B.
    auto at_b = [&](bool gradient_at_b) {
        EnergyCalculator calc("gfnff", config);
        calc.setMolecule(mol_a);
        calc.CalculateEnergy(true);
        calc.CalculateEnergy(true);
        calc.updateGeometry(mol_b.m_geometry);
        return calc.CalculateEnergy(gradient_at_b);
    };
    const double e_energy_only = at_b(false);
    const double e_gradient    = at_b(true);

    const double diff = e_energy_only - e_gradient;
    std::printf("backend=%s  E(energy-only at B)=%.12f  E(gradient call at B)=%.12f  diff=%.3e Eh\n",
                gpu.c_str(), e_energy_only, e_gradient, diff);

    // Claude Generated (Oct 2026): frozen coordination numbers (-gfnff.static_cn). The ROCm fused
    // dispersion gather forms C6 per edge in gradient calls and leaves the per-pair C6 array stale;
    // an energy-only call with frozen CN does not refresh the Gaussian weights and must refill that
    // array from the stored per-atom data. Setup happens at A, the first (gradient) call at B, so
    // the stale per-pair C6 belongs to A and the frozen CN to B; the energy-only call at B must equal
    // the gradient call at B.
    json config_static = config;
    config_static["gfnff"]["static_cn"] = true;
    EnergyCalculator calc_static("gfnff", config_static);
    calc_static.setMolecule(mol_a);
    calc_static.updateGeometry(mol_b.m_geometry);
    const double e_grad_b = calc_static.CalculateEnergy(true);
    const double e_eo_same = calc_static.CalculateEnergy(false);
    const double diff_same = e_eo_same - e_grad_b;
    std::printf("backend=%s  static_cn: E(gradient call at B)=%.12f  E(energy-only at B again)=%.12f  diff=%.3e Eh\n",
                gpu.c_str(), e_grad_b, e_eo_same, diff_same);
    return (std::abs(diff) <= 1e-8 && std::abs(diff_same) <= 1e-8) ? 0 : 1;
}
