/*
 * Regression test: a native QM method must follow updateGeometry().
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Until Sep 2026 the QM engine (then `DFT`) did not override UpdateMolecule(),
 * so after updateGeometry() it kept the first geometry's integrals and SCF and
 * returned the first energy again -- invisible to single-point tests, fatal for
 * any scan/opt/MD. This test computes geometry 2 after geometry 1 on one method
 * object and compares with a fresh object, for `hf` and `hf-3c`.
 *
 * Also checks (Sep 2026): the SCF warm start (previous geometry's orbitals) must
 * reach the same energy as a cold SAD start in fewer iterations, and setMolecule()
 * with a DIFFERENT molecule on the same object must not reuse the old integrals.
 *
 * Claude Generated (Sep 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/hf3c_method.h"
#include "src/core/energy_calculators/qm_methods/qm_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"

#include <cmath>
#include <cstdio>
#include <memory>

static Mol makeBH(double r)
{
    Mol m;
    m.m_number_atoms = 2;
    m.m_atoms = { 5, 1 };
    m.m_geometry = Geometry::Zero(2, 3);
    m.m_geometry(1, 2) = r;
    m.m_charge = 0;
    m.m_spin = 0;
    return m;
}

template <typename Make>
static bool check(const char* name, Make make)
{
    const Mol a = makeBH(1.232), b = makeBH(1.400);
    auto moving = make();
    moving->setMolecule(a);
    const double e_a = moving->calculateEnergy(false);
    moving->updateGeometry(b.m_geometry);
    const double e_b_moved = moving->calculateEnergy(false);

    auto fresh = make();
    fresh->setMolecule(b);
    const double e_b_fresh = fresh->calculateEnergy(false);

    const double d = e_b_moved - e_b_fresh;
    const bool ok = std::abs(d) < 1e-9 && std::abs(e_a - e_b_fresh) > 1e-4;
    std::printf("%-6s E(a)=%.10f  E(b) after update=%.10f  fresh=%.10f  diff=%.2e  %s\n",
                name, e_a, e_b_moved, e_b_fresh, d, ok ? "ok" : "FAIL");
    return ok;
}

static Mol makeH2O()
{
    Mol m;
    m.m_number_atoms = 3;
    m.m_atoms = { 8, 1, 1 };
    m.m_geometry = Geometry::Zero(3, 3);
    m.m_geometry << 0.0, 0.0, 0.0, 0.758602, 0.0, 0.504284, -0.758602, 0.0, 0.504284;
    m.m_charge = 0;
    m.m_spin = 0;
    return m;
}

// Warm start: step a distorted water along a small displacement with one engine
// (warm) and compare with a cold engine at every step. Same energy, fewer iterations.
static bool checkWarmStart()
{
    json cfg = { { "basis", "def2-SVP" }, { "scf_threshold", 1e-9 } };
    json cold_cfg = cfg;
    cold_cfg["scf_warm_start"] = false;
    Mol m = makeH2O();
    QMEngine warm(QMFunctional::HF, cfg), cold(QMFunctional::HF, cold_cfg);
    warm.QMInterface::InitialiseMolecule(m);
    cold.QMInterface::InitialiseMolecule(m);
    warm.Calculation(false);
    cold.Calculation(false);
    int it_warm = 0, it_cold = 0, n_warm_started = 0;
    double dmax = 0.0;
    for (int step = 1; step <= 5; ++step) {
        m.m_geometry(1, 0) += 0.01;  // stretch one O-H by 0.01 A per step
        m.m_geometry(2, 2) -= 0.005;
        warm.UpdateMolecule(Matrix(m.m_geometry));
        cold.UpdateMolecule(Matrix(m.m_geometry));
        const double ew = warm.Calculation(false), ec = cold.Calculation(false);
        dmax = std::max(dmax, std::abs(ew - ec));
        it_warm += warm.scfIterations();
        it_cold += cold.scfIterations();
        n_warm_started += warm.lastScfWarmStarted() ? 1 : 0;
    }
    const bool ok = dmax < 1e-9 && it_warm < it_cold && n_warm_started == 5 && !cold.lastScfWarmStarted();
    std::printf("warm start: 5 steps, SCF iterations warm %d vs cold %d, max|dE| %.2e, warm-started %d/5  %s\n",
                it_warm, it_cold, dmax, n_warm_started, ok ? "ok" : "FAIL");
    return ok;
}

// setMolecule() with a different molecule on the same object must equal a fresh object.
static bool checkNewMolecule()
{
    const json cfg = { { "qm", { { "basis", "MINIX" }, { "scf_threshold", 1e-9 } } } };
    auto reused = std::make_unique<HF3CMethod>(cfg);
    reused->setMolecule(makeBH(1.232));
    reused->calculateEnergy(false);
    reused->setMolecule(makeH2O());
    const double e_reused = reused->calculateEnergy(false);
    auto fresh = std::make_unique<HF3CMethod>(cfg);
    fresh->setMolecule(makeH2O());
    const double e_fresh = fresh->calculateEnergy(false);
    const bool ok = std::abs(e_reused - e_fresh) < 1e-9;
    std::printf("new molecule on a reused object: %.10f vs fresh %.10f  %s\n", e_reused, e_fresh, ok ? "ok" : "FAIL");
    return ok;
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const json cfg = { { "qm", { { "basis", "MINIX" }, { "scf_threshold", 1e-9 } } } };
    bool ok = true;
    ok &= check("hf", [&] { return std::make_unique<QMMethod>(QMFunctional::HF, cfg); });
    ok &= check("hf-3c", [&] { return std::make_unique<HF3CMethod>(cfg); });
    ok &= checkWarmStart();
    ok &= checkNewMolecule();
    return ok ? 0 : 1;
}
