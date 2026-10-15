/**
 * @file test_gfnff_stale_cn.cpp
 * @brief Regression guard (Known Issue #32): energy-only GFN-FF calls on a REUSED calculator
 *        must use the CN and D4 C6 of their own geometry.
 *
 * Before the fix two CN-dependent quantities were refreshed on gradient calls only:
 *   Fix A - the CN read by the Coulomb self-energy chi(CN) = chi_base + cnf*sqrt(CN)
 *   Fix B - the per-pair D4 C6(CN), frozen at the topology-build geometry
 * so an energy-only Calculation(false) on a reused instance (NumGrad, line searches, batch
 * reuse, energy-based Hessians) evaluated them at another geometry. A per-structure
 * single-point benchmark never reuses an instance and cannot see this.
 *
 * Checks (the three of reactff2-llm's CLI test gfnff/06_stale_cn_energy_only, expressed on
 * the GFNFF class directly because this branch has no -batch mode):
 *   1. FD      - Cl2- (q=-1) at 2.73 A: analytic gradient vs GFNFF::NumGrad(), whose
 *                energies are energy-only calls on the same instance.        < 1e-4 Eh/A
 *   2. E vs G  - reused instance, energy-only vs gradient call at the same geometry:
 *                Cl2- built at 2.73 A, evaluated at 2.65 A; triose, 3 frames perturbed
 *                by up to 0.05 A, each reached WITHOUT a prior gradient call. < 1e-10 Eh
 *   3. COULOMB - reused energy-only Coulomb term at 2.65 A == a fresh instance's. < 1e-9 Eh
 *   4. FRESH   - reused energy-only energy == a fresh instance's at the same geometry:
 *                Cl2- total at 2.65 A (< 1e-9 Eh); triose D4 dispersion term over the 3
 *                frames (< 1e-10 Eh). Checks 1-3 are satisfied by fix A alone (FD residual
 *                then 2.4e-5, and the gradient call does not refresh C6 either); these fail
 *                when fix B (the D4 C6 refresh) is missing. The triose TOTAL is not compared:
 *                a reused instance only re-detects its H-bond list on an RMSD gate, which
 *                leaves a separate ~3e-7 Eh H-bond difference on one frame (not this bug).
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (Sep 2026) - AI-generated, machine-tested; human production testing pending
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>

#include <Eigen/Dense>

#include "test_molecule_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"

namespace {

constexpr double kBohrToAng = 0.529177210903;
int g_fails = 0;

void check(const std::string& name, bool ok, double value, const std::string& what)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.2e", value);
    std::cout << (ok ? "PASS " : "FAIL ") << name << ": " << what << " = " << buf << std::endl;
    if (!ok)
        ++g_fails;
}

json config()
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = 1;
    c["cache_topology"] = false;  // never pick up a .topo.json written for another geometry
    return c;
}

/// Cl2- from the registry (2.73 A), second atom moved along the bond to length r (Angstrom).
Mol cl2_anion(double r)
{
    Mol m = TestMolecules::TestMoleculeRegistry::createMolecule("Cl2", false).getMolInfo();
    const Eigen::Vector3d a = m.m_geometry.row(0).transpose();
    const Eigen::Vector3d b = m.m_geometry.row(1).transpose();
    m.m_geometry.row(1) = (a + (b - a).normalized() * r).transpose();
    m.m_charge = -1;
    return m;
}

} // namespace

int main()
{
    CurcumaLogger::set_verbosity(0);
    std::cout << "\n=== GFN-FF energy-only calls on a reused calculator (stale CN / C6) ===\n"
              << std::endl;

    // 1. FD on a reused calculator: NumGrad() differentiates energy-only calls on this instance
    {
        GFNFF ff(config());
        if (!ff.InitialiseMolecule(cl2_anion(2.73))) {
            std::cout << "FAIL InitialiseMolecule(Cl2-)" << std::endl;
            return 1;
        }
        ff.Calculation(true);
        const Matrix g_an = ff.Gradient();     // Eh/Bohr
        const Matrix g_fd = ff.NumGrad(1e-4);  // Eh/Bohr
        const double dev = (g_an - g_fd).cwiseAbs().maxCoeff() / kBohrToAng;
        check("FD Cl2- 2.73 A, reused energy-only", dev < 1e-4, dev, "max |g_an - g_fd| [Eh/A]");
    }

    // 2./3. energy-only vs gradient call on a reused calculator, and Coulomb vs a fresh one
    {
        GFNFF ff(config());
        ff.InitialiseMolecule(cl2_anion(2.73));
        ff.Calculation(false);
        ff.UpdateMolecule(cl2_anion(2.65).m_geometry);
        const double e_only = ff.Calculation(false);
        const double coul_reused = ff.CoulombEnergy();
        const double e_grad = ff.Calculation(true);

        GFNFF fresh(config());
        fresh.InitialiseMolecule(cl2_anion(2.65));
        const double e_fresh = fresh.Calculation(false);
        const double coul_fresh = fresh.CoulombEnergy();

        const double de = std::abs(e_only - e_grad);
        const double dc = std::abs(coul_reused - coul_fresh);
        check("E vs G Cl2- 2.73 -> 2.65 A", de < 1e-10, de, "|E_energy-only - E_gradient| [Eh]");
        check("Coulomb reused vs fresh Cl2- 2.65 A", dc < 1e-9, dc, "|dCoulomb| [Eh]");
        const double dt = std::abs(e_only - e_fresh);
        check("Total E reused vs fresh Cl2- 2.65 A", dt < 1e-9, dt, "|E_reused - E_fresh| [Eh]");
    }

    // 2b. same on a larger polar molecule, several frames, gradient call only AFTER each E call
    {
        const Mol mol = TestMolecules::TestMoleculeRegistry::createMolecule("triose", false).getMolInfo();
        GFNFF ff(config());
        ff.InitialiseMolecule(mol);
        ff.Calculation(true);
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(-0.05, 0.05);
        double d = 0.0, d_disp = 0.0;
        for (int f = 0; f < 3; ++f) {
            Geometry g = mol.m_geometry;
            for (int i = 0; i < g.rows(); ++i)
                for (int k = 0; k < 3; ++k)
                    g(i, k) += u(rng);
            ff.UpdateMolecule(g);
            const double e_only = ff.Calculation(false);
            const double disp_reused = ff.DispersionEnergy();  // energy-only call's D4 term
            const double e_grad = ff.Calculation(true);
            d = std::max(d, std::abs(e_only - e_grad));

            Mol frame = mol;
            frame.m_geometry = g;
            GFNFF fresh(config());
            fresh.InitialiseMolecule(frame);
            fresh.Calculation(false);
            d_disp = std::max(d_disp, std::abs(disp_reused - fresh.DispersionEnergy()));
        }
        check("E vs G triose (3 perturbed frames)", d < 1e-10, d, "max |E_energy-only - E_gradient| [Eh]");
        check("D4 dispersion reused vs fresh triose (3 frames)", d_disp < 1e-10, d_disp, "max |E_disp,reused - E_disp,fresh| [Eh]");
    }

    std::cout << (g_fails ? "FAILURES PRESENT" : "ALL PASS") << "\n" << std::endl;
    return g_fails ? 1 : 0;
}
