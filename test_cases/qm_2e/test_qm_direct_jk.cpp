/*
 * Integral-direct J/K and SCF vs the stored-tensor path.
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * (1) qmint::DirectJK::build() against buildERI + buildCoulomb/buildExchange with
 *     screening off, for a random symmetric density, spherical 5d and cartesian 6d,
 *     1 and 4 threads -- must agree to rounding.
 * (2) The HF SCF with -qm.scf_direct on (incremental Fock builds, density-weighted
 *     screening at the default 1e-12) must reach the stored-tensor energy.
 * (3) The same along a geometry sequence (UpdateMolecule + warm start), which
 *     checks that the direct shell-pair tables are rebuilt per geometry.
 *
 * Claude Generated (Sep 2026). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/energy_calculators/qm_methods/qm_integrals.hpp"
#include "core/test_molecule_registry.h"

#include <cmath>
#include <cstdio>
#include <random>

// Molecules come from the shared test registry (test_cases/CLAUDE.md), in Angstrom.
static Mol registryMol(const char* name)
{
    return TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
}
static Mol water() { return registryMol("H2O"); }
static Mol methanol() { return registryMol("CH3OH"); }

// (1) kernel check
static bool checkKernel(const char* name, const Mol& mol, bool cartesian)
{
    json cfg = { { "basis", "def2-SVP" }, { "cartesian_d", cartesian } };
    QMEngine eng(QMFunctional::HF, cfg);
    eng.QMInterface::InitialiseMolecule(mol);
    const auto& basis = eng.gtoBasis();
    Matrix Q;
    if (!cartesian) Q = qmint::buildSphericalTransform(basis, qmint::buildOverlap(basis));
    const Matrix* Qp = Q.size() ? &Q : nullptr;

    const qmint::ERITensor eri = qmint::buildERI(basis, 1, 0.0, Qp);
    const int n = eri.n();
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(-0.5, 0.5);
    Matrix P(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) P(i, j) = P(j, i) = u(rng);
    const Matrix J0 = qmint::buildCoulomb(eri, P), K0 = qmint::buildExchange(eri, P);

    bool ok = true;
    for (int threads : { 1, 4 }) {
        const qmint::DirectJK direct(basis, 0.0, Qp);
        Matrix J, K;
        direct.build(P, J, K, threads);
        const double dj = (J - J0).cwiseAbs().maxCoeff(), dk = (K - K0).cwiseAbs().maxCoeff();
        const bool good = direct.n() == n && dj < 1e-11 && dk < 1e-11;
        std::printf("kernel %-14s %s n=%3d threads=%d  max|dJ|=%.2e max|dK|=%.2e  %s\n", name,
                    cartesian ? "6d" : "5d", n, threads, dj, dk, good ? "ok" : "FAIL");
        ok &= good;
    }
    return ok;
}

// (2) SCF energy, stored vs direct
static bool checkSCF(const char* name, const Mol& mol)
{
    json base = { { "basis", "def2-SVP" }, { "scf_threshold", 1e-9 }, { "threads", 4 } };
    json stored = base, direct = base;
    stored["scf_direct"] = "off";
    direct["scf_direct"] = "on";
    QMEngine es(QMFunctional::HF, stored), ed(QMFunctional::HF, direct);
    es.QMInterface::InitialiseMolecule(mol);
    ed.QMInterface::InitialiseMolecule(mol);
    const double e_s = es.Calculation(false), e_d = ed.Calculation(false);
    const bool ok = es.scfConverged() && ed.scfConverged() && std::abs(e_s - e_d) < 1e-9;
    std::printf("scf    %-14s stored %.10f (%d it)  direct %.10f (%d it)  diff %.2e  %s\n", name, e_s,
                es.scfIterations(), e_d, ed.scfIterations(), e_d - e_s, ok ? "ok" : "FAIL");
    return ok;
}

// (3) geometry sequence with warm start, direct vs stored
static bool checkGeometrySequence()
{
    json base = { { "basis", "def2-SVP" }, { "scf_threshold", 1e-9 } };
    json stored = base, direct = base;
    stored["scf_direct"] = "off";
    direct["scf_direct"] = "on";
    Mol m = water();
    QMEngine es(QMFunctional::HF, stored), ed(QMFunctional::HF, direct);
    es.QMInterface::InitialiseMolecule(m);
    ed.QMInterface::InitialiseMolecule(m);
    double dmax = std::abs(es.Calculation(false) - ed.Calculation(false));
    double e_first = 0.0, e_last = 0.0;
    for (int step = 1; step <= 4; ++step) {
        m.m_geometry(1, 0) += 0.02;
        es.UpdateMolecule(Matrix(m.m_geometry));
        ed.UpdateMolecule(Matrix(m.m_geometry));
        const double a = es.Calculation(false), b = ed.Calculation(false);
        dmax = std::max(dmax, std::abs(a - b));
        if (step == 1) e_first = b;
        e_last = b;
    }
    const bool ok = dmax < 1e-9 && std::abs(e_last - e_first) > 1e-4;
    std::printf("geometry sequence (4 steps, warm start): max|E_direct - E_stored| %.2e  %s\n", dmax, ok ? "ok" : "FAIL");
    return ok;
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    bool ok = true;
    ok &= checkKernel("water", water(), false);
    ok &= checkKernel("water", water(), true);
    ok &= checkKernel("methanol", methanol(), false);
    ok &= checkSCF("water", water());
    ok &= checkSCF("methanol", methanol());
    ok &= checkGeometrySequence();
    return ok ? 0 : 1;
}
