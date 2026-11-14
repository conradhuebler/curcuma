/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic classical terms of the native GFN1/GFN2 (WP-PERIODIC-XTB P1) - Claude Generated (Nov 2026).
 *
 * The XTB engine is used directly (no SCF): coordination number of the self-energies, repulsion,
 * halogen bond (GFN1, zero for water) and dispersion (GFN1 D3, GFN2 D4 incl. ATM) at fixed charges,
 * on the registry water dimer.
 *   (a) vacuum limit: a 60 A cube equals the cell-free value (every image beyond the cutoffs);
 *   (b) supercell: the dimer in an 8 A cube against the same dimer replicated into a 16 x 8 x 8 cell,
 *       energy per cell equal, CN per atom equal;
 *   (c) translation of one molecule by a lattice vector and an atom-wise wrap leave everything unchanged;
 *   (d) brute force: repulsion of the small cell against a direct 13^3-image sum written here;
 *   (e) dispersion gradient (D3 / D4 with ATM and CN chain rule) against central differences of the
 *       periodic dispersion energy at fixed charges (h = 1e-4 A).
 * Exit 0 when every check passes; each failure prints one line.
 */
#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/parameters/gfn1_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn2_params.hpp"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "test_molecule_registry.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using curcuma::xtb::MethodType;
using curcuma::xtb::XTB;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}
static void checkClose(double a, double b, double tol, const std::string& what)
{
    const double d = std::abs(a - b);
    if (!(d <= tol)) {
        ++g_failures;
        std::printf("FAIL: %s: %.15e vs %.15e (diff %.3e > %.1e)\n", what.c_str(), a, b, d, tol);
    } else {
        std::printf("ok   %s: diff %.3e\n", what.c_str(), d);
    }
}

static Mol waterDimer()
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    return m.getMolInfo();
}

static std::unique_ptr<XTB> makeEngine(MethodType t, const Mol& mol, double L_x, double L_y, double L_z)
{
    auto x = std::make_unique<XTB>(t);
    x->InitialiseMolecule(mol);
    if (L_x > 0.0) {
        Eigen::Matrix3d cell = Eigen::Matrix3d::Zero();
        cell(0, 0) = L_x; cell(1, 1) = L_y; cell(2, 2) = L_z;
        std::string err;
        check(x->setPeriodicCell(cell, err), "setPeriodicCell: " + err);
    }
    return x;
}

static Vector testCharges(const Mol& mol)
{
    Vector q(mol.m_number_atoms);
    for (int i = 0; i < mol.m_number_atoms; ++i) q(i) = (mol.m_atoms[i] == 8) ? -0.6 : 0.3;
    return q;
}

struct Terms {
    double rep = 0.0, xb = 0.0, disp = 0.0;
    Vector cn;
    Matrix dgrad;
};
static Terms terms(XTB& x, const Vector& q, bool grad)
{
    Terms t;
    t.rep = x.repulsionEnergyNow();
    t.xb = x.halogenBondEnergyNow();
    t.cn = x.coordinationNumbersNow();
    t.disp = x.dispersionEnergyNow(q, grad);
    if (grad) t.dgrad = x.dispersionGradientNow();
    return t;
}

static const char* mname(MethodType t) { return t == MethodType::GFN1 ? "gfn1" : "gfn2"; }

// ---------------------------------------------------------------------------------------------
static void testVacuumLimit(MethodType t)
{
    const Mol mol = waterDimer();
    const Vector q = testCharges(mol);
    auto free = makeEngine(t, mol, 0, 0, 0);
    auto big = makeEngine(t, mol, 60.0, 60.0, 60.0);
    const Terms a = terms(*free, q, true), b = terms(*big, q, true);
    const std::string p = std::string(mname(t)) + " vacuum limit (60 A cube vs no cell) ";
    checkClose(a.rep, b.rep, 1e-12, p + "repulsion");
    checkClose(a.xb, b.xb, 1e-14, p + "halogen bond");
    checkClose(a.disp, b.disp, 1e-11, p + "dispersion");
    checkClose((a.cn - b.cn).cwiseAbs().maxCoeff(), 0.0, 1e-12, p + "CN");
    checkClose((a.dgrad - b.dgrad).cwiseAbs().maxCoeff(), 0.0, 1e-11, p + "dispersion gradient");
}

// ---------------------------------------------------------------------------------------------
static void testSupercell(MethodType t)
{
    const Mol mol = waterDimer();
    const double L = 8.0;
    Mol super = mol;
    super.m_number_atoms = 2 * mol.m_number_atoms;
    super.m_geometry.resize(super.m_number_atoms, 3);
    super.m_atoms.resize(super.m_number_atoms);
    for (int i = 0; i < mol.m_number_atoms; ++i) {
        super.m_geometry.row(i) = mol.m_geometry.row(i);
        super.m_geometry.row(i + mol.m_number_atoms) = mol.m_geometry.row(i) + Eigen::RowVector3d(L, 0.0, 0.0);
        super.m_atoms[i] = mol.m_atoms[i];
        super.m_atoms[i + mol.m_number_atoms] = mol.m_atoms[i];
    }
    Vector q1 = testCharges(mol), q2(super.m_number_atoms);
    q2 << q1, q1;
    auto one = makeEngine(t, mol, L, L, L);
    auto two = makeEngine(t, super, 2 * L, L, L);
    const Terms a = terms(*one, q1, true), b = terms(*two, q2, true);
    const std::string p = std::string(mname(t)) + " 2x1x1 supercell (8 A) ";
    checkClose(2.0 * a.rep, b.rep, 1e-12 * std::abs(b.rep), p + "repulsion per cell");
    checkClose(2.0 * a.disp, b.disp, 1e-11 * std::abs(b.disp), p + "dispersion per cell");
    double dcn = 0.0, dg = 0.0;
    for (int i = 0; i < mol.m_number_atoms; ++i) {
        dcn = std::max(dcn, std::abs(a.cn(i) - b.cn(i)));
        dcn = std::max(dcn, std::abs(a.cn(i) - b.cn(i + mol.m_number_atoms)));
        for (int k = 0; k < 3; ++k) {
            dg = std::max(dg, std::abs(a.dgrad(i, k) - b.dgrad(i, k)));
            dg = std::max(dg, std::abs(a.dgrad(i, k) - b.dgrad(i + mol.m_number_atoms, k)));
        }
    }
    checkClose(dcn, 0.0, 1e-12, p + "CN per atom");
    checkClose(dg, 0.0, 1e-11, p + "dispersion gradient per atom");
    // images matter in this cell: the periodic values differ from the cell-free ones
    auto free = makeEngine(t, mol, 0, 0, 0);
    const Terms f = terms(*free, q1, false);
    // (GFN1 repulsion decays as exp(-a r^1.5): below 1e-12 at 8 A; GFN2 with r^1 shows 2e-8)
    check(std::abs(a.disp - f.disp) > 1e-8, p + "dispersion sees images");
    check(std::abs(a.cn(0) - f.cn(0)) > 1e-7, p + "CN sees images");
    std::printf("info %s 8 A cube: E_rep %.12f (free %.12f), E_disp %.12f (free %.12f), CN(O) %.8f (free %.8f)\n",
                mname(t), a.rep, f.rep, a.disp, f.disp, a.cn(0), f.cn(0));
}

// ---------------------------------------------------------------------------------------------
static void testTranslationAndWrap(MethodType t)
{
    const Mol mol = waterDimer();
    const Vector q = testCharges(mol);
    const double L = 8.0;
    auto ref = makeEngine(t, mol, L, L, L);
    const Terms a = terms(*ref, q, true);
    // second molecule (atoms 3-5) shifted by (L, -L, 2L)
    Mol shifted = mol;
    for (int i = 3; i < 6; ++i) shifted.m_geometry.row(i) += Eigen::RowVector3d(L, -L, 2 * L);
    auto s = makeEngine(t, shifted, L, L, L);
    const Terms b = terms(*s, q, true);
    // atom-wise wrap into [0, L) after moving the whole dimer by (0.37 L, 0.61 L, 0.83 L)
    Mol torn = mol;
    for (int i = 0; i < 6; ++i) {
        torn.m_geometry.row(i) += Eigen::RowVector3d(0.37 * L, 0.61 * L, 0.83 * L);
        for (int k = 0; k < 3; ++k) torn.m_geometry(i, k) -= L * std::floor(torn.m_geometry(i, k) / L);
    }
    auto w = makeEngine(t, torn, L, L, L);
    const Terms c = terms(*w, q, true);
    const std::string p = std::string(mname(t)) + " 8 A cube ";
    checkClose(a.rep, b.rep, 1e-12, p + "molecule shifted by a lattice vector: repulsion");
    checkClose(a.disp, b.disp, 1e-11, p + "molecule shifted by a lattice vector: dispersion");
    checkClose((a.cn - b.cn).cwiseAbs().maxCoeff(), 0.0, 1e-12, p + "shift: CN");
    checkClose((a.dgrad - b.dgrad).cwiseAbs().maxCoeff(), 0.0, 1e-11, p + "shift: dispersion gradient");
    checkClose(a.rep, c.rep, 1e-11, p + "atom-wise wrap: repulsion");
    checkClose(a.disp, c.disp, 1e-10, p + "atom-wise wrap: dispersion");
    checkClose((a.cn - c.cn).cwiseAbs().maxCoeff(), 0.0, 1e-11, p + "wrap: CN");
    checkClose((a.dgrad - c.dgrad).cwiseAbs().maxCoeff(), 0.0, 1e-10, p + "wrap: dispersion gradient");
}

// ---------------------------------------------------------------------------------------------
// Brute-force repulsion: every image n in [-6, 6]^3 of an 8 A cube within 25 Bohr (6 * 15.1 Bohr > 25 + 13).
static double bruteRepulsion(MethodType t, const Mol& mol, double L_A)
{
    const double au = 1.0 / 0.529177210903;
    const double L = L_A * au, cut2 = 25.0 * 25.0;
    const int nat = mol.m_number_atoms;
    const bool g1 = (t == MethodType::GFN1);
    using namespace curcuma::xtb;
    double e = 0.0;
    for (int i = 0; i < nat; ++i)
        for (int j = 0; j < nat; ++j) {
            const int zi = mol.m_atoms[i], zj = mol.m_atoms[j];
            const double alfi = g1 ? gfn1_params::rep_alpha[zi - 1] : gfn2_params::rep_alpha[zi - 1];
            const double alfj = g1 ? gfn1_params::rep_alpha[zj - 1] : gfn2_params::rep_alpha[zj - 1];
            const double zi_eff = g1 ? gfn1_params::rep_zeff[zi - 1] : gfn2_params::rep_zeff[zi - 1];
            const double zj_eff = g1 ? gfn1_params::rep_zeff[zj - 1] : gfn2_params::rep_zeff[zj - 1];
            double kexp = g1 ? gfn1_params::rep_kexp : gfn2_params::rep_kexp;
            if (!g1 && zi <= 2 && zj <= 2) kexp = gfn2_params::rep_kexp_light;
            const double rexp = g1 ? gfn1_params::rep_rexp : gfn2_params::rep_rexp;
            const double ap = std::sqrt(alfi * alfj);
            for (int nx = -6; nx <= 6; ++nx)
                for (int ny = -6; ny <= 6; ++ny)
                    for (int nz = -6; nz <= 6; ++nz) {
                        if (i == j && nx == 0 && ny == 0 && nz == 0) continue;
                        const double dx = (mol.m_geometry(i, 0) - mol.m_geometry(j, 0)) * au - nx * L;
                        const double dy = (mol.m_geometry(i, 1) - mol.m_geometry(j, 1)) * au - ny * L;
                        const double dz = (mol.m_geometry(i, 2) - mol.m_geometry(j, 2)) * au - nz * L;
                        const double r2 = dx * dx + dy * dy + dz * dz;
                        if (r2 >= cut2) continue;
                        const double r = std::sqrt(r2);
                        e += 0.5 * zi_eff * zj_eff / std::pow(r, rexp) * std::exp(-ap * std::pow(r, kexp));
                    }
        }
    return e;
}

static void testBruteForce(MethodType t)
{
    const Mol mol = waterDimer();
    auto x = makeEngine(t, mol, 8.0, 8.0, 8.0);
    const double e = x->repulsionEnergyNow();
    const double b = bruteRepulsion(t, mol, 8.0);
    checkClose(e, b, 1e-13 * std::max(1.0, std::abs(b)), std::string(mname(t)) + " repulsion vs 13^3-image brute force (8 A cube)");
}

// ---------------------------------------------------------------------------------------------
static void testDispersionFD(MethodType t, double L)
{
    const Mol mol = waterDimer();
    const Vector q = testCharges(mol);
    auto x = makeEngine(t, mol, L, L, L);
    x->dispersionEnergyNow(q, true);
    const Matrix g = x->dispersionGradientNow();   // Eh/Bohr
    const double h = 1.0e-4;                        // Angstrom
    const double au = 1.0 / 0.529177210903;
    double max_err = 0.0, max_g = 0.0;
    for (int i = 0; i < mol.m_number_atoms; ++i)
        for (int k = 0; k < 3; ++k) {
            Matrix gp = mol.m_geometry, gm = mol.m_geometry;
            gp(i, k) += h;
            gm(i, k) -= h;
            x->UpdateMolecule(gp);
            const double ep = x->dispersionEnergyNow(q, false);
            x->UpdateMolecule(gm);
            const double em = x->dispersionEnergyNow(q, false);
            const double fd = (ep - em) / (2.0 * h) / au;   // Eh/Bohr
            max_err = std::max(max_err, std::abs(fd - g(i, k)));
            max_g = std::max(max_g, std::abs(g(i, k)));
        }
    x->UpdateMolecule(mol.m_geometry);
    std::printf("info %s dispersion gradient vs FD (h 1e-4 A, L %.0f A): max |diff| %.3e Eh/Bohr, max |g| %.3e\n",
                mname(t), L, max_err, max_g);
    check(max_err < 2e-8, std::string(mname(t)) + " dispersion gradient vs FD < 2e-8 Eh/Bohr, L = " + std::to_string(L));
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    for (MethodType t : { MethodType::GFN1, MethodType::GFN2 }) {
        testVacuumLimit(t);
        testSupercell(t);
        testTranslationAndWrap(t);
        testBruteForce(t);
        testDispersionFD(t, 0.0);
        testDispersionFD(t, 8.0);
    }
    if (g_failures == 0)
        std::printf("test_pbc_xtb_classical: all checks passed\n");
    else
        std::printf("test_pbc_xtb_classical: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
