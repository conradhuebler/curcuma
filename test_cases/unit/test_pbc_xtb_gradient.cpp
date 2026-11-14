/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic native GFN1/GFN2 gradient: finite-difference gate of the assembled surface (WP-PERIODIC-XTB P3)
 * - Claude Generated (Nov 2026).
 *
 * Registry water dimer in an 8 A cube (strong image interactions, every term active) and in a 12 A cube,
 * through EnergyCalculator (gradient in Eh/Angstrom) with a tight SCF (scf_threshold 1e-9, so that the
 * finite differences of the energy and the analytic gradient refer to the same converged state). Central
 * differences with h = 1e-4 A in Angstrom, as the root CLAUDE.md gate demands. Tolerance 2e-6 Eh/A (the
 * vacuum xtb_gradient tests use 5e-4; measured values are printed). Additionally the periodic gradient must
 * be translation invariant (forces sum to zero) and the GFN1 forces are compared with tblite 0.4.0
 * (tblite_periodic_reference gfn1 ... 8 8 8 --gradient, 2026-11-14, Eh/Bohr) to 1e-5 Eh/Bohr (tblite's
 * gradient is at its accuracy-1 convergence; its periodic ATM gradient issue concerns GFN2 only).
 * Exit 0 when every check passes.
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

static Mol dimerInCell(double L)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    Mol mol = m.getMolInfo();
    mol.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    mol.m_has_pbc = true;
    return mol;
}

// tblite_periodic_reference gfn1 water-dimer.roo2.8.xyz 8 8 8 --gradient (accuracy 1.0), Eh/Bohr
static const double kTbliteGfn1L8[6][3] = {
    { -7.52098928038721e-17, -1.94593934700634e-16, 1.62394217157041e-02 },
    { -3.67990189206180e-17, 2.52744751799374e-04, -6.40680079489740e-03 },
    { 1.10910556203349e-16, -2.52744751799172e-04, -6.40680079489705e-03 },
    { 1.46521375646534e-17, 1.39862080250630e-17, -1.62394217157029e-02 },
    { -2.74464165895013e-18, 2.52744751798836e-04, 6.40680079489671e-03 },
    { -1.08091449918214e-17, -2.52744751798833e-04, 6.40680079489654e-03 },
};

static void fdGate(const char* method, double L)
{
    const double BOHR = 0.529177210903;
    const json cfg = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    Mol mol = dimerInCell(L);
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    check(!calc.Error(), std::string(method) + " setMolecule with a cell");
    calc.CalculateEnergy(true);
    check(!calc.Error(), std::string(method) + " energy + gradient with a cell");
    const Matrix g = calc.Gradient();   // Eh/A
    const int nat = mol.m_number_atoms;
    check(g.rows() == nat && g.cols() == 3, "gradient shape");
    if (g.rows() != nat) return;

    const double h = 1.0e-4;
    double max_err = 0.0, max_g = 0.0;
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) {
            Matrix gp = mol.m_geometry, gm = mol.m_geometry;
            gp(i, k) += h;
            gm(i, k) -= h;
            calc.updateGeometry(gp);
            const double ep = calc.CalculateEnergy(false);
            calc.updateGeometry(gm);
            const double em = calc.CalculateEnergy(false);
            const double fd = (ep - em) / (2.0 * h);
            max_err = std::max(max_err, std::abs(fd - g(i, k)));
            max_g = std::max(max_g, std::abs(g(i, k)));
        }
    calc.updateGeometry(mol.m_geometry);
    const Eigen::RowVector3d fsum = g.colwise().sum();
    std::printf("%s L = %.0f A: analytic vs central FD (h 1e-4 A): max |diff| %.3e Eh/A, max |g| %.3e, |sum g| %.3e\n",
                method, L, max_err, max_g, fsum.norm());
    check(max_err < 2.0e-6, std::string(method) + " periodic gradient vs FD < 2e-6 Eh/A, L = " + std::to_string(L));
    check(fsum.norm() < 1.0e-7, std::string(method) + " periodic forces sum to zero, L = " + std::to_string(L));

    if (std::string(method) == "gfn1" && L == 8.0) {
        double dmax = 0.0;
        for (int i = 0; i < nat; ++i)
            for (int k = 0; k < 3; ++k) dmax = std::max(dmax, std::abs(g(i, k) * BOHR - kTbliteGfn1L8[i][k]));
        std::printf("gfn1 L = 8 A: max |g - g_tblite| %.3e Eh/Bohr\n", dmax);
        check(dmax < 1.0e-5, "gfn1 periodic gradient vs tblite, L = 8");
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    fdGate("gfn1", 8.0);
    fdGate("gfn2", 8.0);
    fdGate("gfn1", 12.0);
    fdGate("gfn2", 12.0);
    if (g_failures == 0)
        std::printf("test_pbc_xtb_gradient: all checks passed\n");
    else
        std::printf("test_pbc_xtb_gradient: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
