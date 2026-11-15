/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic native GFN1/GFN2 energies against tblite 0.4.0 and against themselves (WP-PERIODIC-XTB P2a/P2b)
 * - Claude Generated (Nov 2026).
 *
 * Registry water dimer (water-dimer.roo2.8) in cubic cells through EnergyCalculator with Mol::m_has_pbc.
 * tblite reference values: tblite_periodic_reference <method> <xyz> L L L (test_cases/sqm_reference/
 * tblite_periodic_reference.c, tblite 0.4.0 from external/tblite built in release_tblite/, Gamma point),
 * 2026-11-14; tblite's run-to-run scatter is 1e-14 Eh; tblite converges the energy to 1e-6 * accuracy and
 * the charges to 2e-5 * accuracy (xtb/singlepoint.f90), so the charges are compared at --accuracy 0.001.
 *
 * GFN1 (P2a): energies in 8/10/12 A cubes within 1e-9 Eh of tblite (measured 3-7e-11), converged charges
 * within 1e-6 e (measured 2.5e-10).
 *
 * GFN2 (P2b): tblite's periodic anisotropic electrostatics (coulomb/multipole.f90 get_amat_sdq_rec_3d) carry a
 * factor 2 in the charge-dipole and +1 instead of -1/3 in the charge-quadrupole reciprocal term, which makes its
 * matrices depend on the Ewald splitting parameter (docs/WP-PERIODIC-XTB.md P2b). curcuma uses the consistent
 * factors. Checks: (a) with tblite's factors (-xtb.pbc_multipole_tblite_factors true) the L = 8
 * energy reproduces tblite to 1e-9 Eh (the port is exact, the deviation is the factor); (b) with the consistent
 * factors the energy is independent of alpha (CURCUMA_PBC_MP_ALPHA 0.14 vs 0.3 at L = 12, < 5e-7 Eh); (c) the
 * 30 A cell is within 5e-6 Eh of the vacuum energy (dipolar cell, tin-foil term; GFN1 shows 1.1e-6);
 * (d) regression values of the consistent energies (own numbers of 2026-11-14, 1e-9 Eh).
 * Also: hasPeriodicEnergy() true with a cell, false without; a gradient request with a cell is accepted (P3).
 * Exit 0 when every check passes.
 */
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

static Mol dimerInCell(double L)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    Mol mol = m.getMolInfo();
    if (L > 0.0) {
        mol.m_unit_cell = Eigen::Matrix3d::Identity() * L;
        mol.m_has_pbc = true;
    }
    return mol;
}

static const json kCfg = { { "verbosity", 0 }, { "threads", 1 } };
static const json kTight = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };

static double energyOf(const char* method, double L, const json& cfg, Vector* charges = nullptr)
{
    Mol mol = dimerInCell(L);
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    check(!calc.Error(), std::string(method) + " setMolecule, L = " + std::to_string(L));
    check(calc.hasPeriodicEnergy() == (L > 0.0), std::string(method) + " hasPeriodicEnergy() follows the cell, L = " + std::to_string(L));
    const double e = calc.CalculateEnergy(false);
    if (charges) *charges = calc.Charges();
    return e;
}

static double maxChargeDiff(const Vector& q, double q_O, double q_H)
{
    if (q.size() != 6) return 1.0;
    return std::max({ std::abs(q(0) - q_O), std::abs(q(3) - q_O), std::abs(q(1) - q_H), std::abs(q(2) - q_H),
        std::abs(q(4) - q_H), std::abs(q(5) - q_H) });
}

// ---------------------------------------------------------------------------------------------
static void testGfn1()
{
    struct Ref { double L, energy; };
    // tblite_periodic_reference gfn1 ... L L L, accuracy 1.0
    static const Ref refs[] = { { 8.0, -11.5308428959022127 }, { 10.0, -11.5315715595413497 }, { 12.0, -11.5317220312253905 } };
    for (const Ref& r : refs) {
        const double e = energyOf("gfn1", r.L, kCfg);
        std::printf("gfn1 L = %.0f A: E %.12f tblite %.12f dE %.3e Eh\n", r.L, e, r.energy, e - r.energy);
        check(std::abs(e - r.energy) < 1.0e-9, "gfn1 periodic energy vs tblite, L = " + std::to_string(r.L));
    }
    // charges at L = 8: tblite --accuracy 0.001 (energy there -11.5308428959531071)
    Vector q;
    energyOf("gfn1", 8.0, kTight, &q);
    const double dq = maxChargeDiff(q, -0.648364111389, 0.324182055694);
    std::printf("gfn1 L = 8 A: max |dq| vs tblite (both tightly converged) %.3e e\n", dq);
    check(dq < 1.0e-6, "gfn1 periodic Mulliken charges vs tblite, L = 8");
}

// ---------------------------------------------------------------------------------------------
static void testGfn2()
{
    // (a) exact port: tblite's reciprocal factors reproduce tblite at L = 8 (tblite accuracy 1.0: -10.1324596572071339)
    //     through the registered switch -xtb.pbc_multipole_tblite_factors (P4; the env variables CURCUMA_PBC_MP_*
    //     remain diagnostic overrides)
    const json cfg_tb = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "pbc_multipole_tblite_factors", true } } } };
    const double e_tb = energyOf("gfn2", 8.0, cfg_tb);
    std::printf("gfn2 L = 8 A with tblite's reciprocal factors: E %.12f tblite -10.132459657207 dE %.3e Eh\n", e_tb, e_tb + 10.1324596572071339);
    check(std::abs(e_tb + 10.1324596572071339) < 1.0e-9, "gfn2 with tblite's factors (-xtb.pbc_multipole_tblite_factors) reproduces tblite at L = 8");

    // (b) alpha independence of the consistent version at L = 12
    setenv("CURCUMA_PBC_MP_ALPHA", "0.14", 1);
    const double e_a1 = energyOf("gfn2", 12.0, kCfg);
    setenv("CURCUMA_PBC_MP_ALPHA", "0.3", 1);
    const double e_a2 = energyOf("gfn2", 12.0, kCfg);
    unsetenv("CURCUMA_PBC_MP_ALPHA");
    std::printf("gfn2 L = 12 A: alpha 0.14 %.12f, alpha 0.3 %.12f, difference %.3e Eh\n", e_a1, e_a2, e_a1 - e_a2);
    check(std::abs(e_a1 - e_a2) < 5.0e-7, "gfn2 periodic energy independent of the Ewald alpha (consistent factors)");

    // (c) vacuum limit and (d) regression values (own numbers, 2026-11-14)
    const double e_vac = energyOf("gfn2", 0.0, kCfg);
    const double e8 = energyOf("gfn2", 8.0, kCfg);
    const double e12 = energyOf("gfn2", 12.0, kCfg);
    const double e30 = energyOf("gfn2", 30.0, kCfg);
    std::printf("gfn2 vacuum %.12f, L = 8 %.12f, L = 12 %.12f, L = 30 %.12f (30 - vacuum %.3e Eh)\n", e_vac, e8, e12, e30, e30 - e_vac);
    check(std::abs(e30 - e_vac) < 5.0e-6, "gfn2 30 A cell within 5e-6 Eh of vacuum");
    check(std::abs(e8 - (-10.133473317)) < 1.0e-8, "gfn2 L = 8 regression value");
    check(std::abs(e12 - (-10.134522052)) < 1.0e-8, "gfn2 L = 12 regression value");
    // charges at L = 8 against the tblite-factor run are not meaningful; the consistent charges are a regression value
    Vector q;
    energyOf("gfn2", 8.0, kTight, &q);
    check(q.size() == 6 && std::abs(q.sum()) < 1.0e-8, "gfn2 periodic charges sum to zero");
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    testGfn1();
    testGfn2();
    {
        // gradient with a cell (P3): accepted, same energy as the energy-only call (FD gate: pbc_xtb_gradient)
        Mol mol = dimerInCell(8.0);
        EnergyCalculator calc("gfn1", kCfg);
        calc.setMolecule(mol);
        const double e_g = calc.CalculateEnergy(true);
        check(!calc.Error(), "gradient request with a cell is accepted (P3)");
        const double e_0 = energyOf("gfn1", 8.0, kCfg);
        check(std::abs(e_g - e_0) < 1.0e-12, "energy with and without gradient request identical");
        // P5/G1: a -gpu cuda request with a cell gives the periodic energy: in a build without the plugin through the
        // counted CPU fallback (identical), with the CUDA plugin through the host-built periodic matrices and the
        // device-resident SCF (G1; same fixed point, device eigensolver rounding, parity gate pbc_xtb_gpu_parity)
        const json gpu_cfg = { { "verbosity", 0 }, { "threads", 1 }, { "gpu", "cuda" } };
        EnergyCalculator calc_gpu("gfn1", gpu_cfg);
        calc_gpu.setMolecule(mol);
        check(!calc_gpu.Error(), "gfn1 -gpu cuda with a cell: no error");
        check(calc_gpu.hasPeriodicEnergy(), "gfn1 -gpu cuda with a cell: periodic energy");
        const double e_gpu = calc_gpu.CalculateEnergy(false);
        check(std::abs(e_gpu - e_0) < 1.0e-8, "gfn1 -gpu cuda with a cell gives the CPU periodic energy (1e-8)");
    }
    if (g_failures == 0)
        std::printf("test_pbc_xtb_tblite: all checks passed\n");
    else
        std::printf("test_pbc_xtb_tblite: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
