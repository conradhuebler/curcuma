/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF Coulomb pair term by Ewald summation (WP-PERIODIC-NONBONDED B2) - Claude Generated (Nov 2026).
 *
 * The GFN-FF pair term sum_{i<j} q_i q_j erf(g_ij r)/r (Gaussian charges, g_ij = 1/sqrt(alp_i + alp_j)) summed
 * over all periodic images: point-charge Ewald sum (B0) minus the real-space correction erfc(g_ij r)/r per
 * pair, plus the Gaussian background (pi Q/V) sum q_i alp_i for a charged cell (ewald.h). Charges stay the
 * non-periodic EEQ ones (B3 pending), so every comparison below is at FIXED charges. Switch: coulomb_ewald.
 *
 * (a)  alpha independence with the real GFN-FF alp (water box, 20 H2O, L = 44 Bohr, 3 alphas, tol 1e-10):
 *      Coulomb energy, total gradient, potentials dE_pair/dq.
 * (b)  large-cell limit: H2O dimer (dipole) and benzene (quadrupole) in L = 150 / 200 / 300 Bohr against the
 *      non-periodic Coulomb energy of the plain path (same charges); difference and its scaling with L.
 * (c)  fixed charges: translation, atom-wise wrap (torn molecules), 2x1x1 supercell: Coulomb energy per atom
 *      and Coulomb gradient equal to the base cell (1e-12).
 * (d)  central finite differences at fixed charges (h = 1e-4 and 2e-4 Angstrom, dispersion_tail off), dE/dq by
 *      finite differences against the potentials.
 * (e)  charged cell (H3O+ 2 H2O plus 18 H2O, Q = +1, L = 44 Bohr): alpha independence with background.
 * (e2) EwaldSum with Gaussian alp against an independent pure reciprocal-space sum of Gaussian clouds
 *      (form factor exp(-k^2 alp/4), long double, charged cell): energy and potentials. This is the test that
 *      catches a wrong or missing Gaussian background term, which alpha independence cannot see.
 * (f)  point-charge limit alp -> 0 equals the B0 point-charge Ewald sum.
 * (g0) no cell: coulomb_ewald true is bitwise what coulomb_ewald false gives (the CLI check (g) covers 7 structures).
 *
 * Mode `--perf <xyz> <L_bohr>` (not run by ctest): Coulomb term time, Ewald in the cell vs the plain O(N^2) path.
 * Water geometry: registry "H2O_dimer", moved rigidly only. Exit 0 when every check passes.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/energy_calculators/ff_methods/ewald.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;
using curcuma::ewald::EwaldSum;

static int g_failures = 0;
static void check(bool ok, const char* what, double value, double limit)
{
    std::printf("  [%s] %-64s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what, value, limit);
    if (!ok)
        ++g_failures;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
static const double PI = 3.14159265358979323846;

static double relDiff(double a, double b) { return std::abs(a - b) / std::max(1e-300, std::max(std::abs(a), std::abs(b))); }
static double relMatDiff(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b)
{
    return (a - b).cwiseAbs().maxCoeff() / std::max(1e-300, std::max(a.cwiseAbs().maxCoeff(), b.cwiseAbs().maxCoeff()));
}

static Mol registryMol(const char* name, int charge = 0)
{
    Mol m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
    m.m_charge = charge;
    return m;
}

static Mol withCell(Mol m, double lx, double ly, double lz)
{
    m.m_has_pbc = true;
    Eigen::Matrix3d c = Eigen::Matrix3d::Zero();
    c(0, 0) = lx;
    c(1, 1) = ly;
    c(2, 2) = lz;
    m.m_unit_cell = c;
    return m;
}

/// Rigid random packing of copies of `units` into a cube of L (Angstrom), every intermolecular atom pair at
/// a minimum-image distance >= min_any; units[0] is placed first. Not a liquid.
static Mol randomPack(const std::vector<Mol>& units, double L_ang, unsigned seed, double min_any)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Mol m = units[0];
    m.m_atoms.clear();
    m.m_geometry = Geometry(0, 3);
    int charge = 0;
    auto mic = [&](Eigen::Vector3d v) {
        for (int k = 0; k < 3; ++k)
            v[k] -= L_ang * std::rint(v[k] / L_ang);
        return v;
    };
    for (const Mol& unit : units) {
        const Eigen::Vector3d centre = unit.m_geometry.colwise().mean().transpose();
        for (long tries = 0; tries < 200000; ++tries) {
            Eigen::Vector3d p(u(rng) * L_ang, u(rng) * L_ang, u(rng) * L_ang);
            Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
            qr.normalize();
            const Eigen::Matrix3d R = qr.toRotationMatrix();
            Geometry x(unit.m_number_atoms, 3);
            for (int i = 0; i < unit.m_number_atoms; ++i)
                x.row(i) = (p + R * (unit.m_geometry.row(i).transpose() - centre)).transpose();
            bool ok = true;
            for (int a = 0; a < x.rows() && ok; ++a)
                for (int b = 0; b < m.m_geometry.rows(); ++b)
                    if (mic(x.row(a).transpose() - m.m_geometry.row(b).transpose()).norm() < min_any) {
                        ok = false;
                        break;
                    }
            if (!ok)
                continue;
            const int n0 = static_cast<int>(m.m_geometry.rows());
            m.m_geometry.conservativeResize(n0 + x.rows(), 3);
            m.m_geometry.bottomRows(x.rows()) = x;
            for (int i = 0; i < unit.m_number_atoms; ++i)
                m.m_atoms.push_back(unit.m_atoms[i]);
            charge += unit.m_charge;
            break;
        }
    }
    m.m_number_atoms = static_cast<int>(m.m_geometry.rows());
    m.m_charge = charge;
    m.m_has_pbc = false;
    return m;
}

static Mol waterBox(int nwat, double L_bohr, unsigned seed, const std::vector<Mol>& first = {}, double min_any = 2.6)
{
    const Mol dimer = registryMol("H2O_dimer");
    Mol w = dimer;
    w.m_geometry = dimer.m_geometry.topRows(3);
    w.m_atoms.assign(dimer.m_atoms.begin(), dimer.m_atoms.begin() + 3);
    w.m_number_atoms = 3;
    std::vector<Mol> units = first;
    for (int i = 0; i < nwat; ++i)
        units.push_back(w);
    const double L_ang = L_bohr / A2B;
    return withCell(randomPack(units, L_ang, seed, min_any), L_ang, L_ang, L_ang);
}

struct Run {
    std::unique_ptr<GFNFF> ff;
    bool ok = false;
    double e_tot = 0.0, e_coul = 0.0;
    Matrix grad;    // Eh/Bohr
    Vector q, phi;  // EEQ charges, dE_pair/dq (Ewald only)
};

static json config(const json& extra, int threads = 1)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    return c;
}

static Run runGFNFF(const Mol& mol, const json& extra, bool gradient, bool potentials = false, int threads = 1)
{
    Run r;
    try {
        r.ff = std::make_unique<GFNFF>(config(extra, threads));
        if (!r.ff->InitialiseMolecule(mol))
            return r;
        if (r.ff->getWorkspace())
            r.ff->getWorkspace()->setStoreGradientComponents(true);
        r.ff->setWantCoulombPotentials(potentials);
        r.e_tot = r.ff->Calculation(gradient);
        r.e_coul = r.ff->CoulombEnergy();
        if (gradient)
            r.grad = r.ff->Gradient();
        r.q = r.ff->Charges();
        if (potentials)
            r.phi = r.ff->coulombPotentials();
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("GFNFF threw: %s\n", e.what());
    }
    return r;
}

static json ewaldFlags(double alpha, double r_cut, double tol = 1e-10)
{
    json e = json::object();
    e["coulomb_ewald"] = true;
    e["ewald_tolerance"] = tol;
    if (alpha > 0.0)
        e["ewald_alpha"] = alpha;
    if (r_cut > 0.0)
        e["ewald_r_cut"] = r_cut;
    return e;
}

// ---------------------------------------------------------------------------------------------
// (a) alpha independence, neutral water box
static void testAlphaIndependence(const Mol& box, const char* label, const json& base)
{
    std::printf("(a/e) alpha independence: %s, %d atoms, charge %d\n", label, box.m_number_atoms, box.m_charge);
    const double alphas[3] = { 0.25, 0.32, 0.42 }; // r_cut 19 Bohr: erfc(0.25*19) = 1.8e-11
    Run r[3];
    for (int t = 0; t < 3; ++t) {
        json e = ewaldFlags(alphas[t], 19.0);
        for (auto it = base.begin(); it != base.end(); ++it)
            e[it.key()] = it.value();
        r[t] = runGFNFF(box, e, true, true);
        if (!r[t].ok) {
            check(false, "run", 1, 0);
            return;
        }
        const auto& terms = r[t].ff->getWorkspace()->ewaldTerms();
        std::printf("  alpha %.2f: E_coul %.12f  (real %.6f recip %.6f self %.6f bg %.6f gauss_bg %.6f) N_k %d\n", alphas[t], r[t].e_coul,
            terms.real, terms.reciprocal, terms.self, terms.background, terms.gaussian_background,
            r[t].ff->getWorkspace()->ewald()->numberOfKVectors());
    }
    double dE = 0, dG = 0, dP = 0, dQ = 0;
    for (int t = 1; t < 3; ++t) {
        dE = std::max(dE, relDiff(r[0].e_coul, r[t].e_coul));
        dG = std::max(dG, relMatDiff(r[0].grad, r[t].grad));
        dP = std::max(dP, relMatDiff(r[0].phi, r[t].phi));
        dQ = std::max(dQ, (r[0].q - r[t].q).cwiseAbs().maxCoeff());
    }
    check(dQ == 0.0, "EEQ charges identical across alphas (non-periodic EEQ)", dQ, 1e-300);
    check(dE < 1e-9, "Coulomb energy spread over 3 alphas rel", dE, 1e-9);
    check(dG < 1e-8, "gradient spread (max abs / max |g|)", dG, 1e-8);
    check(dP < 1e-8, "potential spread (max abs / max |phi|)", dP, 1e-8);
}

// ---------------------------------------------------------------------------------------------
// (b) large-cell limit
static void testLargeCellLimit(const char* name, int charge, const char* label)
{
    std::printf("(b) large-cell limit: %s\n", label);
    const Mol mol = registryMol(name, charge);
    json plain = json::object();
    const Run ref = runGFNFF(mol, plain, false);
    if (!ref.ok) {
        check(false, "reference run", 1, 0);
        return;
    }
    Eigen::Vector3d dip = Eigen::Vector3d::Zero();
    for (int i = 0; i < mol.m_number_atoms; ++i)
        dip += ref.q(i) * mol.m_geometry.row(i).transpose() * A2B;
    std::printf("  non-periodic E_coul %.12f Eh, |dipole| %.4f e Bohr, |q| %.4f, sum q %.2e\n", ref.e_coul, dip.norm(), ref.q.norm(), ref.q.sum());
    const double Ls[3] = { 150.0, 200.0, 300.0 };
    double diff[3];
    for (int t = 0; t < 3; ++t) {
        const Mol cell = withCell(mol, Ls[t] / A2B, Ls[t] / A2B, Ls[t] / A2B);
        const Run r = runGFNFF(cell, ewaldFlags(0.0, 0.0, 1e-12), false);
        if (!r.ok) {
            check(false, "periodic run", 1, 0);
            return;
        }
        diff[t] = r.e_coul - ref.e_coul;
        const double dq = (r.q - ref.q).cwiseAbs().maxCoeff();
        std::printf("  L %5.0f Bohr: E_coul(pbc) - E_coul %.3e Eh  (charges equal to %.1e; r_cut %.2f alpha %.4f N_k %d)\n", Ls[t], diff[t], dq,
            r.ff->getWorkspace()->ewald()->parameters().r_cut, r.ff->getWorkspace()->ewald()->parameters().alpha,
            r.ff->getWorkspace()->ewald()->numberOfKVectors());
    }
    const double expo = std::log(std::abs(diff[0]) / std::abs(diff[2])) / std::log(Ls[2] / Ls[0]);
    std::printf("  scaling exponent of |dE| with L: %.2f (3 for a dipole lattice, 5 for a quadrupole lattice)\n", expo);
    check(std::abs(diff[2]) < 2e-6, "|E_coul(pbc, L=300) - E_coul| < 2e-6 Eh", std::abs(diff[2]), 2e-6);
    check(std::abs(diff[2]) < std::abs(diff[0]), "|dE| decreases with L", std::abs(diff[2]) / std::abs(diff[0]), 1.0);
}

// ---------------------------------------------------------------------------------------------
// (c) fixed charges: translation, wrap, supercell
struct Fixed {
    double e_coul = 0.0;
    Matrix g_coul;
};
/// Per-atom Coulomb inputs of the base instance (charges and the Phase-1-dependent alp, gam, chi): the non-periodic
/// Phase-1 EEQ gives a torn or replicated arrangement slightly different topology charges (A4c), so a fixed-charge
/// comparison has to fix these too. `copies` replicates them for a supercell.
struct FixedInputs {
    Vector q, chi_base, gam, alp, cnf, chi_static;
    FixedInputs replicate(int copies) const
    {
        FixedInputs r;
        auto rep = [&](const Vector& v) {
            Vector out(v.size() * copies);
            for (int c = 0; c < copies; ++c)
                out.segment(c * v.size(), v.size()) = v;
            return out;
        };
        r.q = rep(q);
        r.chi_base = rep(chi_base);
        r.gam = rep(gam);
        r.alp = rep(alp);
        r.cnf = rep(cnf);
        r.chi_static = rep(chi_static);
        return r;
    }
};
static Fixed fixedChargeCoulomb(const Mol& mol, const json& flags, const FixedInputs& in)
{
    Fixed f;
    Run r = runGFNFF(mol, flags, true);
    if (!r.ok) {
        check(false, "fixed-charge run", 1, 0);
        return f;
    }
    FFWorkspace* ws = r.ff->getWorkspace();
    ws->setCoulombSelfEnergyParams(in.chi_base, in.gam, in.alp, in.cnf, in.chi_static);
    ws->setEEQCharges(in.q);
    ws->calculate(true);
    f.e_coul = ws->energyComponents().coulomb;
    f.g_coul = ws->gradientCoulomb();
    return f;
}

static void testFixedChargeInvariance(const Mol& box)
{
    std::printf("(c) fixed charges: translation, atom-wise wrap, 2x1x1 supercell (%d atoms)\n", box.m_number_atoms);
    const json flags = ewaldFlags(0.0, 0.0, 1e-10);
    const Run base = runGFNFF(box, flags, true);
    if (!base.ok) {
        check(false, "base run", 1, 0);
        return;
    }
    FixedInputs in;
    {
        const FFWorkspace* ws = base.ff->getWorkspace();
        in.q = base.q;
        in.chi_base = ws->coulombChiBase();
        in.gam = ws->coulombGam();
        in.alp = ws->coulombAlp();
        in.cnf = ws->coulombCnf();
        in.chi_static = ws->coulombChiStatic();
    }
    const Fixed f0 = fixedChargeCoulomb(box, flags, in);
    const int N = box.m_number_atoms;
    const double L = box.m_unit_cell(0, 0);

    Mol shifted = box;
    for (int i = 0; i < N; ++i) {
        shifted.m_geometry(i, 0) += 0.37 * L;
        shifted.m_geometry(i, 1) -= 1.21 * L;
        shifted.m_geometry(i, 2) += 0.05 * L;
    }
    const Fixed f1 = fixedChargeCoulomb(shifted, flags, in);
    check(relDiff(f0.e_coul, f1.e_coul) < 1e-12, "translation: Coulomb energy rel", relDiff(f0.e_coul, f1.e_coul), 1e-12);
    check(relMatDiff(f0.g_coul, f1.g_coul) < 1e-11, "translation: Coulomb gradient rel", relMatDiff(f0.g_coul, f1.g_coul), 1e-11);

    Mol wrapped = shifted;
    for (int i = 0; i < N; ++i)
        for (int a = 0; a < 3; ++a)
            wrapped.m_geometry(i, a) -= L * std::floor(wrapped.m_geometry(i, a) / L);
    const Fixed f2 = fixedChargeCoulomb(wrapped, flags, in);
    check(relDiff(f0.e_coul, f2.e_coul) < 1e-12, "atom-wise wrap: Coulomb energy rel", relDiff(f0.e_coul, f2.e_coul), 1e-12);
    check(relMatDiff(f0.g_coul, f2.g_coul) < 1e-11, "atom-wise wrap: Coulomb gradient rel", relMatDiff(f0.g_coul, f2.g_coul), 1e-11);

    Mol super = box;
    super.m_geometry = Geometry(2 * N, 3);
    super.m_atoms.clear();
    for (int c = 0; c < 2; ++c)
        for (int i = 0; i < N; ++i) {
            super.m_geometry.row(c * N + i) = box.m_geometry.row(i);
            super.m_geometry(c * N + i, 0) += c * L;
            super.m_atoms.push_back(box.m_atoms[i]);
        }
    super.m_number_atoms = 2 * N;
    super.m_charge = 2 * box.m_charge;
    super.m_unit_cell(0, 0) = 2 * L;
    const Fixed f3 = fixedChargeCoulomb(super, flags, in.replicate(2));
    const double dE = relDiff(f0.e_coul / N, f3.e_coul / (2 * N));
    std::printf("  E_coul/N cell %.14f  supercell %.14f\n", f0.e_coul / N, f3.e_coul / (2 * N));
    check(dE < 1e-12, "2x1x1 supercell: Coulomb energy per atom rel", dE, 1e-12);
    check(relMatDiff(f0.g_coul, f3.g_coul.topRows(N)) < 1e-11, "2x1x1 supercell: Coulomb gradient rel (first copy)",
        relMatDiff(f0.g_coul, f3.g_coul.topRows(N)), 1e-11);
    check(relMatDiff(f0.g_coul, f3.g_coul.bottomRows(N)) < 1e-11, "2x1x1 supercell: Coulomb gradient rel (second copy)",
        relMatDiff(f0.g_coul, f3.g_coul.bottomRows(N)), 1e-11);

    // thread count: 4 threads bitwise equal to 1 thread is NOT required (the reduction order differs);
    // the partitions are fixed, so the result at 4 threads is reproducible and equal to 1e-13.
    Run t4 = runGFNFF(box, flags, true, true, 4);
    Run t1 = runGFNFF(box, flags, true, true, 1);
    check(relDiff(t1.e_coul, t4.e_coul) < 1e-13, "4 threads vs 1 thread: Coulomb energy rel", relDiff(t1.e_coul, t4.e_coul), 1e-13);
    check(relMatDiff(t1.grad, t4.grad) < 1e-12, "4 threads vs 1 thread: gradient rel", relMatDiff(t1.grad, t4.grad), 1e-12);
    check(relMatDiff(t1.phi, t4.phi) < 1e-12, "4 threads vs 1 thread: potentials rel", relMatDiff(t1.phi, t4.phi), 1e-12);
    Run t4b = runGFNFF(box, flags, true, true, 4);
    check(t4.e_coul == t4b.e_coul && (t4.grad - t4b.grad).cwiseAbs().maxCoeff() == 0.0, "4 threads: run-to-run bitwise", 0.0, 1e-300);
}

// ---------------------------------------------------------------------------------------------
// (d) finite differences at fixed charges
/// Central FD of the workspace Coulomb component at fixed charges (CN recomputed per displacement, as
/// GFNFF::NumGradFixedCharges does for the total); returns max |FD - gradientCoulomb| and its atom.
static double coulombFD(GFNFF& ff, const Mol& mol, double h_bohr, int& worst_atom)
{
    FFWorkspace* ws = ff.getWorkspace();
    Eigen::MatrixXd geom = mol.m_geometry * A2B;
    const PBCUtils::OrthoCell* cell = ws->hasPeriodicCell() ? &ws->periodicCell() : nullptr;
    auto energy_at = [&]() {
        auto cn = CNCalculator::calculateGFNFFCN(mol.m_atoms, geom, 1600.0, -7.5, 4.4, nullptr, cell);
        ws->setGeometry(geom);
        const Vector cnv = Vector::Map(cn.data(), cn.size()).eval();
        ws->setD3CN(cnv);
        ws->setCN(cnv);
        ws->calculate(false);
        return ws->energyComponents().coulomb;
    };
    energy_at();
    ws->calculate(true);
    const Matrix g = ws->gradientCoulomb();
    double worst = 0.0;
    worst_atom = -1;
    for (int i = 0; i < mol.m_number_atoms; ++i)
        for (int a = 0; a < 3; ++a) {
            geom(i, a) += h_bohr;
            const double ep = energy_at();
            geom(i, a) -= 2 * h_bohr;
            const double em = energy_at();
            geom(i, a) += h_bohr;
            const double d = std::abs((ep - em) / (2 * h_bohr) - g(i, a));
            if (d > worst) {
                worst = d;
                worst_atom = i;
            }
        }
    energy_at();
    ws->calculate(true);
    return worst;
}

static void testFiniteDifferences(const Mol& box, json flags, const char* label, bool with_checks)
{
    std::printf("(d) finite differences at fixed charges: %s (%d atoms, dispersion_tail off)\n", label, box.m_number_atoms);
    flags["dispersion_tail"] = false; // NumGradFixedCharges keeps the last tail; its dE/dCN would not be in the FD
    Run r = runGFNFF(box, flags, true, true);
    if (!r.ok) {
        check(false, "run", 1, 0);
        return;
    }
    const double gmax = r.grad.cwiseAbs().maxCoeff();
    const double hs[2] = { 1e-4 * A2B, 2e-4 * A2B };
    double err[2], cerr[2];
    for (int t = 0; t < 2; ++t) {
        const Matrix fd = r.ff->NumGradFixedCharges(hs[t]);
        Eigen::Index wi, wj;
        err[t] = (fd - r.grad).cwiseAbs().maxCoeff(&wi, &wj);
        int catom = -1;
        cerr[t] = coulombFD(*r.ff, box, hs[t], catom);
        std::printf("  h = %.0e Angstrom: total max |FD - analytic| = %.3e Eh/Bohr at atom %d (max |g| %.3e); Coulomb only %.3e at atom %d\n",
            hs[t] / A2B, err[t], int(wi), gmax, cerr[t], catom);
    }
    std::printf("  ratio err(2h)/err(h) = %.2f (total), %.2f (Coulomb)\n", err[1] / err[0], cerr[1] / cerr[0]);
    if (!with_checks)
        return;
    // The total shows an h-independent 1.5e-6 Eh/Bohr discrepancy that is also there without a cell and without
    // Ewald (`--fd-diagnose`, not in the Coulomb component): pre-existing, outside B2, reported in TODO.md.
    // B2 is judged on the Coulomb component: truncation error quadratic in h.
    check(cerr[0] / gmax < 5e-7, "Coulomb gradient vs FD of the Coulomb component, h = 1e-4 A (rel to max |g|)", cerr[0] / gmax, 5e-7);
    check(cerr[1] / cerr[0] > 3.0 && cerr[1] / cerr[0] < 5.0, "Coulomb FD error quadratic in h (ratio in 3..5)", cerr[1] / cerr[0], 5.0);
    check(err[0] / gmax < 1e-3, "total gradient vs NumGradFixedCharges (loose, pre-existing non-Coulomb discrepancy)", err[0] / gmax, 1e-3);

    // dE_pair/dq against the potentials: E_pair is quadratic in q, the central difference is exact up to rounding
    FFWorkspace* ws = r.ff->getWorkspace();
    ws->setGeometry(box.m_geometry * A2B);
    ws->calculate(true);
    const Vector q0 = r.q;
    const double hq = 1e-3;
    double perr = 0.0;
    for (int i = 0; i < q0.size(); ++i) {
        Vector qq = q0;
        qq(i) += hq;
        ws->setEEQCharges(qq);
        ws->calculate(false);
        const double ep = ws->ewaldPairEnergy();
        qq(i) -= 2 * hq;
        ws->setEEQCharges(qq);
        ws->calculate(false);
        const double em = ws->ewaldPairEnergy();
        perr = std::max(perr, std::abs((ep - em) / (2 * hq) - r.phi(i)));
    }
    ws->setEEQCharges(q0);
    ws->calculate(false);
    const double pmax = r.phi.cwiseAbs().maxCoeff();
    check(perr / pmax < 1e-9, "potentials vs FD dE_pair/dq (rel to max |phi|)", perr / pmax, 1e-9);
    const double ehalf = 0.5 * q0.dot(r.phi);
    check(relDiff(ehalf, ws->ewaldPairEnergy()) < 1e-12, "E_pair == 1/2 sum q phi rel", relDiff(ehalf, ws->ewaldPairEnergy()), 1e-12);
}

// ---------------------------------------------------------------------------------------------
// (e2) independent reciprocal-space reference for Gaussian clouds (charged cell), (f) point-charge limit
static void gaussianKSpaceReference(const Eigen::MatrixXd& pos, const Eigen::VectorXd& q, const Eigen::VectorXd& alp,
    const PBCUtils::OrthoCell& cell, int nk, long double& E, std::vector<long double>& phi)
{
    // E = (1/2V) sum_{k != 0} (4 pi/k^2) |sum_i q_i exp(-k^2 alp_i/4) e^{i k r_i}|^2 - sum_i q_i^2 / sqrt(2 pi alp_i)
    // Absolutely convergent for Gaussian clouds; the k = 0 term is dropped (same convention as the Ewald sum).
    const int N = int(q.size());
    const long double pi = 3.14159265358979323846264338327950288L;
    const long double V = (long double)cell.L[0] * cell.L[1] * cell.L[2];
    phi.assign(N, 0.0L);
    std::vector<long double> cs(N), sn(N), f(N);
    for (int a = -nk; a <= nk; ++a)
        for (int b = -nk; b <= nk; ++b)
            for (int c = -nk; c <= nk; ++c) {
                if (a == 0 && b == 0 && c == 0)
                    continue;
                long double k[3] = { 2 * pi * a / cell.L[0], 2 * pi * b / cell.L[1], 2 * pi * c / cell.L[2] };
                long double k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                long double w = 4 * pi / (V * k2);
                long double C = 0, S = 0;
                for (int i = 0; i < N; ++i) {
                    long double t = k[0] * pos(i, 0) + k[1] * pos(i, 1) + k[2] * pos(i, 2);
                    f[i] = std::exp(-k2 * (long double)alp(i) / 4);
                    cs[i] = f[i] * std::cos(t);
                    sn[i] = f[i] * std::sin(t);
                    C += q(i) * cs[i];
                    S += q(i) * sn[i];
                }
                for (int i = 0; i < N; ++i)
                    phi[i] += w * (cs[i] * C + sn[i] * S);
            }
    E = 0;
    for (int i = 0; i < N; ++i) {
        phi[i] -= 2 * q(i) / std::sqrt(2 * pi * (long double)alp(i));
        E += 0.5L * q(i) * phi[i];
    }
}

static void testGaussianReference()
{
    std::printf("(e2) Gaussian clouds, charged cell 12 x 13 x 14 Bohr: Ewald + correction vs pure k-space sum\n");
    Eigen::Matrix3d cm = Eigen::Matrix3d::Zero();
    cm(0, 0) = 12.0;
    cm(1, 1) = 13.0;
    cm(2, 2) = 14.0;
    std::string err;
    const PBCUtils::OrthoCell cell = PBCUtils::makeOrthoCell(cm, err);
    const int N = 10;
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Eigen::MatrixXd pos(N, 3);
    Eigen::VectorXd q(N), alp(N);
    for (int i = 0; i < N; ++i) {
        for (int a = 0; a < 3; ++a)
            pos(i, a) = u(rng) * cell.L[a];
        q(i) = 2.0 * u(rng) - 1.0;
        alp(i) = 0.34 + 3.9 * u(rng); // GFN-FF alpeeq range: H 0.34 ... Cl 4.2 Bohr^2
    }
    q(0) += 0.7; // Q != 0: the background and the Gaussian background enter
    std::printf("  Q = %.4f e, alp in [%.2f, %.2f], r_g = %.2f Bohr\n", q.sum(), alp.minCoeff(), alp.maxCoeff(),
        curcuma::ewald::gaussianCorrectionRadius(alp));
    // The real-space loop must reach r_g and stay below L_min/2 = 6 Bohr; with the GFN-FF widths r_g = 17 Bohr,
    // i.e. in this small test cell the clouds would overlap their own images (that is what the L/2 check in
    // GFNFF refuses). The formula test therefore uses narrower clouds, alp in [0.08, 0.38] Bohr^2: r_g = 5.2.
    alp = (0.05 + 0.33 * alp.array() / 4.24).matrix();
    std::printf("  used alp in [%.3f, %.3f], r_g = %.2f Bohr (< L_min/2 = 6)\n", alp.minCoeff(), alp.maxCoeff(),
        curcuma::ewald::gaussianCorrectionRadius(alp));
    // nk = 60: k up to 2 pi 60/14 = 27 1/Bohr, the smallest cloud pair carries exp(-k^2 alp_min/2) = exp(-29) there
    long double Eref;
    std::vector<long double> phiref;
    gaussianKSpaceReference(pos, q, alp, cell, 60, Eref, phiref);
    std::printf("  reference E = %.14Lf\n", Eref);
    double spread = 0.0, dE = 0.0, dphi = 0.0, gauss_bg = 0.0;
    double E_first = 0.0;
    for (int t = 0; t < 3; ++t) {
        const double rc[3] = { 5.9, 5.5, 5.2 };
        EwaldSum ew(cell, curcuma::ewald::chooseEwaldParameters(cell, 1e-13, rc[t]));
        const auto r = ew.evaluate(pos, q, false, true, &alp);
        if (t == 0) {
            E_first = r.terms.total();
            gauss_bg = r.terms.gaussian_background;
            dE = std::abs(r.terms.total() - (double)Eref) / std::abs((double)Eref);
            double maxphi = 0.0;
            for (int i = 0; i < N; ++i) {
                dphi = std::max(dphi, std::abs(r.potentials(i) - (double)phiref[i]));
                maxphi = std::max(maxphi, std::abs((double)phiref[i]));
            }
            dphi /= maxphi;
            // the same without the Gaussian background: the k-space reference must then disagree
            const double no_bg = r.terms.total() - r.terms.gaussian_background;
            std::printf("  Ewald E = %.14f (gaussian background %.6e); without it rel error %.3e\n", r.terms.total(),
                r.terms.gaussian_background, std::abs(no_bg - (double)Eref) / std::abs((double)Eref));
            check(std::abs(no_bg - (double)Eref) / std::abs((double)Eref) > 1e-6, "k-space reference discriminates the Gaussian background",
                std::abs(no_bg - (double)Eref) / std::abs((double)Eref), 1e-6);
        }
        spread = std::max(spread, relDiff(E_first, r.terms.total()));
    }
    (void)gauss_bg;
    check(dE < 1e-11, "energy vs pure k-space Gaussian sum rel", dE, 1e-11);
    check(dphi < 1e-11, "potentials vs pure k-space Gaussian sum (rel to max |phi|)", dphi, 1e-11);
    check(spread < 1e-12, "alpha spread rel (3 r_cut)", spread, 1e-12);

    // gradient of the Gaussian form against finite differences (charged cell)
    {
        EwaldSum ew(cell, curcuma::ewald::chooseEwaldParameters(cell, 1e-13, 5.9));
        const auto r = ew.evaluate(pos, q, true, false, &alp);
        double fderr = 0.0;
        const double h = 1e-4;
        for (int i = 0; i < N; ++i)
            for (int a = 0; a < 3; ++a) {
                Eigen::MatrixXd p = pos;
                p(i, a) += h;
                const double ep = ew.energy(p, q, &alp);
                p(i, a) -= 2 * h;
                const double em = ew.energy(p, q, &alp);
                fderr = std::max(fderr, std::abs((ep - em) / (2 * h) - r.gradient(i, a)));
            }
        check(fderr / r.gradient.cwiseAbs().maxCoeff() < 1e-7, "Gaussian gradient vs FD (rel to max |g|)",
            fderr / r.gradient.cwiseAbs().maxCoeff(), 1e-7);
        const Eigen::MatrixXd A = ew.matrix(pos, &alp);
        const auto rp = ew.evaluate(pos, q, false, true, &alp);
        check(relMatDiff(A * q, rp.potentials) < 1e-12, "dense Gaussian A q == phi", relMatDiff(A * q, rp.potentials), 1e-12);
    }

    std::printf("(f) point-charge limit alp = 1e-12\n");
    {
        EwaldSum ew(cell, curcuma::ewald::chooseEwaldParameters(cell, 1e-13, 5.9));
        Eigen::VectorXd tiny = Eigen::VectorXd::Constant(N, 1e-12);
        const auto a = ew.evaluate(pos, q, true, true, &tiny);
        const auto b = ew.evaluate(pos, q, true, true);
        check(relDiff(a.terms.total(), b.terms.total()) < 1e-12, "energy == point-charge Ewald rel", relDiff(a.terms.total(), b.terms.total()), 1e-12);
        check(relMatDiff(a.gradient, b.gradient) < 1e-12, "gradient == point-charge Ewald rel", relMatDiff(a.gradient, b.gradient), 1e-12);
        check(relMatDiff(a.potentials, b.potentials) < 1e-12, "potentials == point-charge Ewald rel", relMatDiff(a.potentials, b.potentials), 1e-12);
    }
}

// ---------------------------------------------------------------------------------------------
// (g0) no cell: the switch has no effect
static void testNoCell()
{
    std::printf("(g0) no cell: coulomb_ewald true == false bitwise\n");
    const Mol mol = registryMol("H2O_dimer");
    json on = json::object();
    on["coulomb_ewald"] = true;
    const Run a = runGFNFF(mol, on, true), b = runGFNFF(mol, json::object(), true);
    check(a.ok && b.ok && a.e_tot == b.e_tot && a.e_coul == b.e_coul && (a.grad - b.grad).cwiseAbs().maxCoeff() == 0.0,
        "energy and gradient bitwise", std::abs(a.e_tot - b.e_tot), 1e-300);
    check(a.ok && a.ff->getWorkspace() && !a.ff->getWorkspace()->ewaldActive(), "Ewald inactive without a cell", 0.0, 1e-300);
}

// ---------------------------------------------------------------------------------------------
// perf mode: Coulomb term CPU-sum (termTimings, summed over partitions) and the wall time of three
// repeated workspace calculate(true) calls (median), Ewald in the cell vs the plain O(N^2) path.
static void perfOne(const Mol& mol, const json& flags, int threads, const char* label)
{
    Run r = runGFNFF(mol, flags, true, false, threads);
    if (!r.ok) {
        std::printf("  %s: failed\n", label);
        return;
    }
    FFWorkspace* ws = r.ff->getWorkspace();
    double wall[3], coul[3];
    for (int k = 0; k < 3; ++k) {
        auto t0 = std::chrono::steady_clock::now();
        ws->calculate(true);
        wall[k] = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        coul[k] = ws->termTimings().coulomb + ws->ewaldSetupMs();
    }
    std::sort(wall, wall + 3);
    std::sort(coul, coul + 3);
    if (ws->ewaldActive())
        std::printf("  %-34s %d thr: r_cut %5.2f alpha %.4f N_k %6d pairs %8zu | Coulomb CPU-sum %8.2f ms | workspace wall %8.2f ms | E_coul %.8f\n",
            label, threads, ws->ewald()->parameters().r_cut, ws->ewald()->parameters().alpha, ws->ewald()->numberOfKVectors(),
            ws->ewald()->realSpacePairs(mol.m_geometry * A2B).size(), coul[1], wall[1], r.e_coul);
    else
        std::printf("  %-34s %d thr: %58s | Coulomb CPU-sum %8.2f ms | workspace wall %8.2f ms | E_coul %.8f\n",
            label, threads, "plain O(N^2), all pairs", coul[1], wall[1], r.e_coul);
}

static void perf(const char* xyz, double L_bohr)
{
    curcuma::Molecule m(xyz);
    Mol mol = m.getMolInfo();
    std::printf("perf: %s, %d atoms, charge %d, cubic cell L = %.1f Bohr\n", xyz, mol.m_number_atoms, mol.m_charge, L_bohr);
    json plain = json::object();
    const Mol cell = withCell(mol, L_bohr / A2B, L_bohr / A2B, L_bohr / A2B);
    for (int threads : { 1, 4 }) {
        perfOne(mol, plain, threads, "no cell");
        for (double rc : { 0.0, 20.0, 30.0, 40.0 }) {
            char label[64];
            std::snprintf(label, sizeof(label), "cell, Ewald r_cut %s", rc > 0 ? std::to_string(int(rc)).c_str() : "auto");
            perfOne(cell, ewaldFlags(0.0, rc, 1e-8), threads, label);
        }
    }
}

static void perfBox(int nwat, double L_bohr)
{
    const Mol box = waterBox(nwat, L_bohr, 20261107, {}, 1.7); // 1.7 A: no intermolecular bonds, packable at water density
    std::printf("perf: random water box, %d atoms in L = %.1f Bohr (density %.4f atoms/A^3)\n", box.m_number_atoms, L_bohr,
        box.m_number_atoms / std::pow(L_bohr / A2B, 3));
    json plain = json::object();
    Mol nocell = box;
    nocell.m_has_pbc = false;
    for (int threads : { 1, 4 }) {
        perfOne(nocell, plain, threads, "no cell");
        for (double rc : { 0.0, 15.0, 20.0, 25.0 }) {
            char label[64];
            std::snprintf(label, sizeof(label), "cell, Ewald r_cut %s", rc > 0 ? std::to_string(int(rc)).c_str() : "auto");
            perfOne(box, ewaldFlags(0.0, rc, 1e-8), threads, label);
        }
    }
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc > 3 && std::strcmp(argv[1], "--perf") == 0) {
        perf(argv[2], std::atof(argv[3]));
        return 0;
    }
    if (argc > 3 && std::strcmp(argv[1], "--perf-box") == 0) {
        perfBox(std::atoi(argv[2]), std::atof(argv[3]));
        return 0;
    }
    if (argc > 1 && std::strcmp(argv[1], "--fd-diagnose") == 0) {
        // FD check of the same box without Ewald (cell, plain Coulomb) and without a cell: locates a
        // discrepancy that is not B2's
        const Mol box = waterBox(20, 44.0, 20261105);
        json plain = json::object();
        testFiniteDifferences(box, plain, "cell, plain Coulomb (coulomb_ewald false)", false);
        Mol nocell = box;
        nocell.m_has_pbc = false;
        testFiniteDifferences(nocell, plain, "no cell", false);
        testFiniteDifferences(box, ewaldFlags(0.0, 0.0, 1e-12), "cell + Ewald", false);
        return 0;
    }
    testNoCell();
    testGaussianReference();

    const Mol box = waterBox(20, 44.0, 20261105);
    testAlphaIndependence(box, "20 H2O, L = 44 Bohr", json::object());
    testFixedChargeInvariance(box);
    testFiniteDifferences(box, ewaldFlags(0.0, 0.0, 1e-12), "cell + Ewald", true);
    testLargeCellLimit("H2O", 0, "water (dipole)");
    testLargeCellLimit("H2O_dimer", 0, "water dimer (registry geometry: no net dipole)");
    testLargeCellLimit("C6H6", 0, "benzene (quadrupole)");
    const Mol charged = waterBox(18, 44.0, 20261106, { registryMol("hydronium_water2", 1) });
    testAlphaIndependence(charged, "H3O+ 2 H2O + 18 H2O, L = 44 Bohr (Q = +1)", json::object());

    std::printf("test_pbc_coulomb_ewald: %d failure(s)\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
