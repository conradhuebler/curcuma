/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Atomic virial of the periodic native GFN1/GFN2: finite-difference gate (WP-PERIODIC-XTB P4) - Claude Generated
 * (Nov 2026).
 *
 * Xi_ab = -dE/d(eps_ab) at fixed fractional coordinates. Gate: central differences of the energy under a diagonal
 * strain (orthorhombic cells only: eps_xx, eps_yy, eps_zz separately and the isotropic strain), positions and cell
 * scaled together, h = 1e-5, tight SCF (scf_threshold 1e-9), registry water dimer in an 8 A cube (gfn1, gfn2) and a
 * 12 A cube (gfn2). Target <= 1e-7 Eh per component (coordinator); the measured values are printed. The trace of the
 * virial must equal the isotropic derivative. For gfn1 at 8 A the virial is compared with tblite 0.4.0's
 * (tblite_periodic_reference gfn1 ... 8 8 8 --gradient, "virial_eh", 2026-11-14; tblite reports sigma = dE/d(eps)
 * at its accuracy-1 convergence, so Xi = -sigma is compared to 1e-5 Eh; gfn2 is not compared, see P2b).
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
    mol.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    mol.m_has_pbc = true;
    return mol;
}

// tblite_periodic_reference gfn1 water-dimer.roo2.8.xyz 8 8 8 --gradient, virial_eh (Eh), 2026-11-14
static const double kTbliteSigmaGfn1L8[3][3] = {
    { -5.66446870854992e-04, -9.98712910954569e-17, -4.72829785990769e-17 },
    { -9.98713346669764e-17, 8.74301250785553e-04, -2.50152546246718e-16 },
    { -4.72827708604510e-17, -2.50152546246718e-16, 1.28454378124731e-02 },
};

static double strainedEnergy(const char* method, const Mol& mol0, const Eigen::Matrix3d& eps)
{
    const json cfg = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    Mol mol = mol0;
    const Eigen::Matrix3d F = Eigen::Matrix3d::Identity() + eps;
    mol.m_unit_cell = F * mol0.m_unit_cell;   // diagonal strains keep the cell orthorhombic
    for (int i = 0; i < mol.m_number_atoms; ++i)
        mol.m_geometry.row(i) = (F * mol0.m_geometry.row(i).transpose()).transpose();
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    const double e = calc.CalculateEnergy(false);
    check(!calc.Error(), std::string(method) + " strained energy");
    return e;
}

static void gate(const char* method, double L, bool compare_tblite)
{
    const json cfg = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    const Mol mol = dimerInCell(L);
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    calc.setWantVirial(true);
    calc.CalculateEnergy(true);
    check(!calc.Error() && calc.virialValid(), std::string(method) + " virial computed with the gradient");
    const Eigen::Matrix3d Xi = calc.virial();

    const double h = 1.0e-5;
    double max_err = 0.0;
    for (int a = 0; a < 3; ++a) {
        Eigen::Matrix3d ep = Eigen::Matrix3d::Zero(), em = Eigen::Matrix3d::Zero();
        ep(a, a) = h;
        em(a, a) = -h;
        const double fd = -(strainedEnergy(method, mol, ep) - strainedEnergy(method, mol, em)) / (2.0 * h);
        max_err = std::max(max_err, std::abs(fd - Xi(a, a)));
        std::printf("%s L = %.0f A: Xi_%c%c analytic %.9e  FD %.9e  diff %.3e Eh\n", method, L, "xyz"[a], "xyz"[a], Xi(a, a), fd, fd - Xi(a, a));
    }
    {
        const Eigen::Matrix3d ep = Eigen::Matrix3d::Identity() * h, em = Eigen::Matrix3d::Identity() * (-h);
        const double fd = -(strainedEnergy(method, mol, ep) - strainedEnergy(method, mol, em)) / (2.0 * h);
        max_err = std::max(max_err, std::abs(fd - Xi.trace()));
        std::printf("%s L = %.0f A: tr Xi analytic %.9e  FD %.9e  diff %.3e Eh\n", method, L, Xi.trace(), fd, fd - Xi.trace());
    }
    check(max_err < 1.0e-7, std::string(method) + " virial vs cell-scaling FD <= 1e-7 Eh, L = " + std::to_string(L));
    check((Xi - Xi.transpose()).norm() < 1e-12, std::string(method) + " virial symmetric");
    if (compare_tblite) {
        double d = 0.0;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) d = std::max(d, std::abs(Xi(a, b) + kTbliteSigmaGfn1L8[a][b]));
        std::printf("%s L = %.0f A: max |Xi + sigma_tblite| %.3e Eh (tblite sigma = dE/d eps)\n", method, L, d);
        check(d < 1.0e-5, "gfn1 virial vs tblite sigma, L = 8");
    }
}

// ---------------------------------------------------------------------------------------------
// Multipole electrostatics alone (XTB probes, fixed moments, fixed damping radii, fixed Ewald alpha through the
// environment override): analytic virial vs FD of the multipole energy under diagonal strains, per moment class.
#include "src/core/energy_calculators/qm_methods/xtb_native.h"
static void multipoleIsolation()
{
    using curcuma::xtb::MethodType;
    const double L = 8.0, h = 1.0e-5;
    const Mol mol = dimerInCell(L);
    curcuma::xtb::XTB x(MethodType::GFN2);
    x.InitialiseMolecule(mol);
    std::string err;
    x.setPeriodicCell(mol.m_unit_cell, err);
    setenv("CURCUMA_PBC_MP_ALPHA", "0.2", 1);
    x.Calculation(false);
    const Vector q0 = x.getCharges();
    const Matrix dp0 = x.dipoleMoments(), qp0 = x.quadrupoleMoments();
    const int nat = mol.m_number_atoms;
    struct Case { const char* name; bool use_q, use_dp, use_qp; };
    const Case cases[] = { { "q+dp (SD, DD)", true, true, false }, { "dp only (DD)", false, true, false },
        { "q+qp (SQ)", true, false, true }, { "all", true, true, true } };
    for (const Case& c : cases) {
        const Vector q = c.use_q ? q0 : Vector::Zero(nat);
        const Matrix dp = c.use_dp ? dp0 : Matrix::Zero(3, nat);
        const Matrix qp = c.use_qp ? qp0 : Matrix::Zero(6, nat);
        x.multipoleEnergyProbe(mol.m_unit_cell, mol.m_geometry, q, dp, qp);
        const Eigen::Matrix3d Xi = x.multipoleVirialProbe(q, dp, qp);
        for (int a = 0; a < 3; ++a) {
            double e[2];
            for (int s = 0; s < 2; ++s) {
                const double f = 1.0 + (s == 0 ? h : -h);
                Matrix g = mol.m_geometry;
                g.col(a) *= f;
                Eigen::Matrix3d cell = mol.m_unit_cell;
                cell(a, a) *= f;
                e[s] = x.multipoleEnergyProbe(cell, g, q, dp, qp);
            }
            const double fd = -(e[0] - e[1]) / (2.0 * h);
            std::printf("mp-iso %-14s Xi_%c%c analytic %.9e FD %.9e diff %.3e\n", c.name, "xyz"[a], "xyz"[a], Xi(a, a), fd, fd - Xi(a, a));
        }
        x.multipoleEnergyProbe(mol.m_unit_cell, mol.m_geometry, q, dp, qp);
    }
    unsetenv("CURCUMA_PBC_MP_ALPHA");
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    if (std::getenv("PBC_VIRIAL_MP_ISO")) { multipoleIsolation(); return 0; }
    gate("gfn1", 8.0, true);
    gate("gfn2", 8.0, false);
    gate("gfn2", 12.0, false);
    if (g_failures == 0)
        std::printf("test_pbc_xtb_virial: all checks passed\n");
    else
        std::printf("test_pbc_xtb_virial: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
