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
 * Claude Generated (Sep 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/hf3c_method.h"
#include "src/core/energy_calculators/qm_methods/qm_method.h"
#include "src/core/curcuma_logger.h"

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

int main()
{
    CurcumaLogger::set_verbosity(0);
    const json cfg = { { "qm", { { "basis", "MINIX" }, { "scf_threshold", 1e-9 } } } };
    bool ok = true;
    ok &= check("hf", [&] { return std::make_unique<QMMethod>(QMFunctional::HF, cfg); });
    ok &= check("hf-3c", [&] { return std::make_unique<HF3CMethod>(cfg); });
    return ok ? 0 : 1;
}
