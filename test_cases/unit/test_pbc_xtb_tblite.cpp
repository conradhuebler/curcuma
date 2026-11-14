/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic native GFN1 against tblite 0.4.0 (WP-PERIODIC-XTB P2a) - Claude Generated (Nov 2026).
 *
 * Registry water dimer (water-dimer.roo2.8) in cubic cells of 8, 10 and 12 A through EnergyCalculator with
 * Mol::m_has_pbc. Reference values: tblite_periodic_reference gfn1 <xyz> L L L (test_cases/sqm_reference/
 * tblite_periodic_reference.c, tblite 0.4.0 from external/tblite built in release_tblite/, accuracy 1.0,
 * Gamma point), 2026-11-14; tblite's run-to-run scatter is 1e-14 Eh. Tolerance 1e-9 Eh (measured 3-7e-11,
 * docs/WP-PERIODIC-XTB.md P2a) and 1e-6 e for the Mulliken charges (SCF threshold of both codes).
 * Also: hasPeriodicEnergy() true with a cell, false without; gfn2 with a cell and a gradient request with a cell are
 * refused (P2b / P3 pending). Exit 0 when every check passes.
 */
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
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

struct Ref {
    double L;
    double energy;
    double q_O;
    double q_H;
};
// tblite_periodic_reference gfn1 water-dimer.roo2.8.xyz L L L (2026-11-14), accuracy 1.0 energies (tblite converges
// the charges only to pconv = 2e-5 * accuracy, the energy to 1e-6 * accuracy; at --accuracy 0.001 the L = 8 energy is
// -11.5308428959531071, 5e-11 away). The charges are the --accuracy 0.001 values (pconv 2e-8).
static const Ref kRefs[] = {
    { 8.0, -11.5308428959022127, -0.648364111389, 0.324182055694 },
    { 10.0, -11.5315715595413497, 0.0, 0.0 },
    { 12.0, -11.5317220312253905, 0.0, 0.0 },
};

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

int main()
{
    CurcumaLogger::set_verbosity(0);
    const json cfg = { { "verbosity", 0 }, { "threads", 1 } };

    for (const Ref& r : kRefs) {
        Mol mol = dimerInCell(r.L);
        EnergyCalculator calc("gfn1", cfg);
        calc.setMolecule(mol);
        check(!calc.Error(), "gfn1 setMolecule with a cell L = " + std::to_string(r.L));
        check(calc.hasPeriodicEnergy(), "hasPeriodicEnergy() true with a cell");
        const double e = calc.CalculateEnergy(false);
        const double dE = e - r.energy;
        std::printf("gfn1 L = %.0f A: E %.12f tblite %.12f dE %.3e Eh\n", r.L, e, r.energy, dE);
        check(std::abs(dE) < 1.0e-9, "gfn1 periodic energy vs tblite, L = " + std::to_string(r.L));
        if (r.q_O != 0.0) {
            // charges: at the default SCF threshold (1e-5 in max |dq_shell|) they are 1e-5 from tblite; a tight
            // threshold shows the converged charges (tblite accuracy 1 converges them to ~1e-6)
            const json tight = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
            EnergyCalculator calc_t("gfn1", tight);
            calc_t.setMolecule(mol);
            calc_t.CalculateEnergy(false);
            const Vector q = calc_t.Charges();
            check(q.size() == 6, "charges returned");
            if (q.size() == 6) {
                const double dq = std::max({ std::abs(q(0) - r.q_O), std::abs(q(3) - r.q_O), std::abs(q(1) - r.q_H),
                    std::abs(q(2) - r.q_H), std::abs(q(4) - r.q_H), std::abs(q(5) - r.q_H) });
                std::printf("gfn1 L = %.0f A: max |dq| vs tblite %.3e e\n", r.L, dq);
                check(dq < 1.0e-6, "gfn1 periodic Mulliken charges vs tblite (both tightly converged), L = 8");
            }
        }
    }

    {
        Mol mol = dimerInCell(0.0);
        EnergyCalculator calc("gfn1", cfg);
        calc.setMolecule(mol);
        check(!calc.Error(), "gfn1 setMolecule without a cell");
        check(!calc.hasPeriodicEnergy(), "hasPeriodicEnergy() false without a cell");
    }
    {
        // gfn2 with a cell: refused until P2b (error, no silent vacuum energy)
        Mol mol = dimerInCell(8.0);
        EnergyCalculator calc("gfn2", cfg);
        bool refused = false;
        try {
            calc.setMolecule(mol);
            refused = calc.Error();
        } catch (const std::exception&) {
            refused = true;
        }
        check(refused, "gfn2 with a cell is refused (P2b pending)");
    }
    {
        // gradient with a cell: refused until P3
        Mol mol = dimerInCell(8.0);
        EnergyCalculator calc("gfn1", cfg);
        calc.setMolecule(mol);
        check(!calc.Error(), "gfn1 setMolecule with a cell (gradient probe)");
        bool threw = false;
        try {
            calc.CalculateEnergy(true);
        } catch (const std::exception&) {
            threw = true;
        }
        check(threw || calc.Error(), "gradient request with a cell is refused (P3 pending)");
    }

    if (g_failures == 0)
        std::printf("test_pbc_xtb_tblite: all checks passed\n");
    else
        std::printf("test_pbc_xtb_tblite: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
