/*
 * test_gfnff_frag_charge_history.cpp - frag_charge_model ensemble: topology-history independence
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 27, 2026) - AI-generated, machine-tested; human production testing pending.
 *
 * Check 5 ("HISTORY") of cli_gfnff_05_frag_charge_ensemble, moved here because it needs ONE
 * calculator instance reused across two geometries, which master's -sp CLI cannot express
 * (the reactff2-llm original used its -batch_reuse_topology mode). See Known Issue #31 and
 * docs/FRAG_CHARGE_MODEL.md.
 *
 * Cl2- is set up at 2.50 A (one fragment: the pass-1 split sits at ~2.64 A) and then moved to
 * 2.70 A (two fragments) WITHOUT a topology rebuild (displacement 0.38 Bohr < the 0.5 Bohr
 * rebuild threshold). The ensemble model re-derives the base fragments from the CURRENT geometry,
 * so the kept-topology energy must equal a fresh evaluation at 2.70 A:
 *   ensemble (s_max 1.1): |E_kept - E_fresh| < 0.01 kcal/mol
 *   reference (liveness): |E_kept - E_fresh| > 1 kcal/mol, which proves the topology was kept
 *
 * All evaluations are GRADIENT calls, as in the original (-gradient true). On master a reused
 * calculator's ENERGY-ONLY call keeps the CN of the previous gradient call (stale-CN bug, fixed on
 * reactff2-llm as its Known Issue #35 and deliberately not part of this port): with the topology
 * unchanged, Cl2- built at 2.40 A and evaluated at 2.55 A is then off by +0.98 kcal/mol
 * energy-only vs +0.0013 kcal/mol with gradient calls, identically before and after this port.
 * The remaining ~1e-3 kcal/mol is the stale D4 C6 of the same issue; it is why the history
 * margin here is 0.0027 vs 0.01 kcal/mol instead of ~0.
 *
 * Usage: test_gfnff_frag_charge_history <path to cl2m_ea25.xyz>  (only the elements are used;
 * the Cl-Cl distance is the scan coordinate of the check, set below).
 */
#include <cmath>
#include <iomanip>
#include <iostream>

#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"

#include "json.hpp"
using json = nlohmann::json;

namespace {
constexpr double kKcal = 627.509474;

Matrix cl2Geometry(double r, double dx, double dy)
{
    Matrix g = Matrix::Zero(2, 3);
    g(1, 0) = dx;
    g(1, 1) = dy;
    g(1, 2) = r;
    return g;
}

json config(const std::string& model)
{
    return { { "verbosity", 0 }, { "threads", 1 },
        { "gfnff", { { "frag_charge_model", model }, { "frag_charge_s_max", 1.1 }, { "cache_topology", false } } } };
}

/// E at 2.70 A with the topology built at 2.50 A (kept) and freshly built at 2.70 A
bool keptAndFresh(const curcuma::Molecule& base, const std::string& model, double& e_kept, double& e_fresh)
{
    const Matrix g_build = cl2Geometry(2.50, 0.0, 0.0);
    const Matrix g_eval = cl2Geometry(2.70, 0.02, -0.01);

    Mol mol = base.getMolInfo();
    mol.m_geometry = g_build;
    EnergyCalculator kept("gfnff", config(model));
    kept.setMolecule(mol);
    kept.CalculateEnergy(true);
    kept.updateGeometry(g_eval);
    e_kept = kept.CalculateEnergy(true);

    mol.m_geometry = g_eval;
    EnergyCalculator fresh("gfnff", config(model));
    fresh.setMolecule(mol);
    e_fresh = fresh.CalculateEnergy(true);
    return !kept.Error() && !fresh.Error() && std::isfinite(e_kept) && std::isfinite(e_fresh);
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: test_gfnff_frag_charge_history <cl2m_ea25.xyz>\n";
        return 2;
    }
    curcuma::Molecule base(argv[1]);
    if (base.AtomCount() != 2) {
        std::cerr << "expected a two-atom Cl2 file, got " << base.AtomCount() << " atoms\n";
        return 2;
    }
    base.setCharge(-1);
    CurcumaLogger::set_verbosity(0);

    int fails = 0;
    std::cout << std::fixed << std::setprecision(4);

    double ek = 0, ef = 0;
    const bool ok_ens = keptAndFresh(base, "ensemble", ek, ef);
    const double d_ens = (ek - ef) * kKcal;
    const bool pass_ens = ok_ens && std::abs(d_ens) < 0.01;
    std::cout << (pass_ens ? "PASS" : "FAIL") << " history ensemble: kept - fresh = " << d_ens << " kcal/mol\n";
    fails += pass_ens ? 0 : 1;

    const bool ok_ref = keptAndFresh(base, "reference", ek, ef);
    const double d_ref = (ek - ef) * kKcal;
    const bool pass_ref = ok_ref && std::abs(d_ref) > 1.0;
    std::cout << (pass_ref ? "PASS" : "FAIL") << " history reference (liveness, topology really kept): kept - fresh = "
              << d_ref << " kcal/mol\n";
    fails += pass_ref ? 0 : 1;

    return fails == 0 ? 0 : 1;
}
