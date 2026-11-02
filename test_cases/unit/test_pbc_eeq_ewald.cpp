/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic EEQ with the Ewald matrix (WP-PERIODIC-NONBONDED B3) - Claude Generated (Nov 2026).
 *
 * With a cell the Phase-2 EEQ charges minimise E(q) = sum_i [-chi_i q_i + 1/2 (gam_i + sqrt(2/pi)/sqrt(alp_i)) q_i^2]
 * + 1/2 q^T A_ew q with the Gaussian-charge Ewald matrix A_ew (ewald.h, B2), and the workspace evaluates the same
 * functional (coulomb_ewald, default on with a cell). Checks:
 *
 * (a)  dense cholesky / lu == matrix-free pcg (same charges to 1e-9, same energy)
 * (b)  alpha_E independence of charges, energy and total gradient (neutral box and a charged cell Q = +1)
 * (c)  large-cell limit L = 200 / 300 Bohr: charges and energy equal to the non-periodic path (dipole image term
 *      -2 pi M^2/(3V) reported)
 * (d)  translation, atom-wise wrap (torn molecules) and 2x1x1 supercell: charges per atom, energy, gradient
 * (e)  THE GATE: analytic total gradient against central finite differences of the total energy with the charges
 *      re-solved at every displaced geometry (h = 1e-4 / 2e-4 Angstrom, Richardson), neutral and Q = +1, dense
 *      (exact) and pcg (eeq_ppcg_tol 1e-14) solver; limit 1e-7 Eh/Angstrom, 1e-8 after Richardson
 * (f)  MIC-only EEQ (minimum-image erf(g r)/r matrix, built here) against the Ewald EEQ on a 50-water box:
 *      largest charge deviation and the EEQ energy difference (the measurement planned as B1)
 * (g0) no cell: coulomb_ewald true (default) is bitwise what coulomb_ewald false gives
 *
 * Mode `--perf-box <nwat> <L_bohr> [threads]` (not run by ctest): Calculation(true) wall time, EEQ solve time and
 * pcg iterations, Coulomb term time, with and without the cell. Water geometry: registry "H2O_dimer", moved rigidly.
 * Exit 0 when every check passes.
 */
#include "src/core/curcuma_logger.h"
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

static int g_failures = 0;
static void check(bool ok, const char* what, double value, double limit)
{
    std::printf("  [%s] %-70s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what, value, limit);
    if (!ok)
        ++g_failures;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
static const double PI = 3.14159265358979323846;

static double relDiff(double a, double b) { return std::abs(a - b) / std::max(1e-300, std::max(std::abs(a), std::abs(b))); }
static double maxAbsDiff(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) { return (a - b).cwiseAbs().maxCoeff(); }

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

/// Rigid random packing of copies of `units` into a cube of L (Angstrom), every intermolecular atom pair at a
/// minimum-image distance >= min_any; units[0] is placed first. Not a liquid.
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

static Mol waterUnit()
{
    const Mol dimer = registryMol("H2O_dimer");
    Mol w = dimer;
    w.m_geometry = dimer.m_geometry.topRows(3);
    w.m_atoms.assign(dimer.m_atoms.begin(), dimer.m_atoms.begin() + 3);
    w.m_number_atoms = 3;
    return w;
}

static Mol waterBox(int nwat, double L_bohr, unsigned seed, const std::vector<Mol>& first = {}, double min_any = 2.6)
{
    std::vector<Mol> units = first;
    for (int i = 0; i < nwat; ++i)
        units.push_back(waterUnit());
    const double L_ang = L_bohr / A2B;
    return withCell(randomPack(units, L_ang, seed, min_any), L_ang, L_ang, L_ang);
}

struct Run {
    std::unique_ptr<GFNFF> ff;
    bool ok = false;
    double e_tot = 0.0, e_coul = 0.0;
    Matrix grad; // Eh/Bohr
    Vector q;
    int iters = 0;
    std::string path;
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

static Run runGFNFF(const Mol& mol, const json& extra, bool gradient, int threads = 1)
{
    Run r;
    try {
        r.ff = std::make_unique<GFNFF>(config(extra, threads));
        if (!r.ff->InitialiseMolecule(mol))
            return r;
        r.e_tot = r.ff->Calculation(gradient);
        r.e_coul = r.ff->CoulombEnergy();
        if (gradient)
            r.grad = r.ff->Gradient();
        r.q = r.ff->Charges();
        if (r.ff->eeqSolver()) {
            r.iters = r.ff->eeqSolver()->lastPeriodicIterations();
            r.path = r.ff->eeqSolver()->lastPeriodicPath();
        }
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("GFNFF threw: %s\n", e.what());
    }
    return r;
}

static json flags(const char* solver = nullptr, double alpha = 0.0, double r_cut = 0.0, double tol = 1e-10)
{
    json e = json::object();
    e["ewald_tolerance"] = tol;
    if (solver)
        e["solve_method"] = solver;
    if (alpha > 0.0)
        e["ewald_alpha"] = alpha;
    if (r_cut > 0.0)
        e["ewald_r_cut"] = r_cut;
    return e;
}

// ---------------------------------------------------------------------------------------------
// (a) dense == matrix-free
static void testDenseVsPcg(const Mol& box, const char* label)
{
    std::printf("(a) dense cholesky / lu vs matrix-free pcg: %s (%d atoms, charge %d)\n", label, box.m_number_atoms, box.m_charge);
    const Run chol = runGFNFF(box, flags("cholesky"), true);
    const Run lu = runGFNFF(box, flags("lu"), true);
    json pcg = flags("pcg");
    pcg["eeq_ppcg_tol"] = 1e-13;
    const Run cg = runGFNFF(box, pcg, true);
    if (!chol.ok || !lu.ok || !cg.ok) {
        check(false, "runs", 1, 0);
        return;
    }
    std::printf("  paths: '%s' / '%s' / '%s' (%d iterations); E_tot %.12f / %.12f / %.12f\n", chol.path.c_str(), lu.path.c_str(),
        cg.path.c_str(), cg.iters, chol.e_tot, lu.e_tot, cg.e_tot);
    check(chol.path == "dense cholesky" && cg.path == "matrix-free pcg", "paths as requested", 0.0, 1e-300);
    check(maxAbsDiff(chol.q, lu.q) < 1e-11, "charges cholesky vs lu (e)", maxAbsDiff(chol.q, lu.q), 1e-11);
    check(maxAbsDiff(chol.q, cg.q) < 1e-9, "charges cholesky vs pcg (e)", maxAbsDiff(chol.q, cg.q), 1e-9);
    check(std::abs(chol.e_tot - cg.e_tot) < 1e-10, "energy cholesky vs pcg (Eh)", std::abs(chol.e_tot - cg.e_tot), 1e-10);
    check(maxAbsDiff(chol.grad, cg.grad) < 1e-9, "gradient cholesky vs pcg (Eh/Bohr)", maxAbsDiff(chol.grad, cg.grad), 1e-9);
    json ppcg = flags("ppcg");
    const Run pp = runGFNFF(box, ppcg, false);
    check(pp.ok && pp.path == "matrix-free pcg", "ppcg redirected to the matrix-free pcg", 0.0, 1e-300);
    bool refused = false;
    try {
        GFNFF ff(config(flags("batched"), 1));
        refused = !ff.InitialiseMolecule(box);
        if (!refused)
            ff.Calculation(false);
    } catch (const std::exception&) {
        refused = true;
    }
    check(refused, "batched refused with a cell", 0.0, 1e-300);
    // auto: dense below the threshold, pcg above
    json a1 = flags("auto");
    a1["eeq_ewald_dense_max_atoms"] = box.m_number_atoms;
    json a2 = flags("auto");
    a2["eeq_ewald_dense_max_atoms"] = box.m_number_atoms - 1;
    const Run r1 = runGFNFF(box, a1, false), r2 = runGFNFF(box, a2, false);
    check(r1.ok && r2.ok && r1.path == "dense cholesky" && r2.path == "matrix-free pcg", "auto: dense up to eeq_ewald_dense_max_atoms, pcg above", 0.0, 1e-300);
    // four threads: pcg and dense reproducible and equal to one thread (not bitwise: reduction order)
    const Run c4 = runGFNFF(box, pcg, true, 4), c4b = runGFNFF(box, pcg, true, 4);
    check(maxAbsDiff(cg.q, c4.q) < 1e-11, "pcg 4 threads vs 1 thread: charges (e)", maxAbsDiff(cg.q, c4.q), 1e-11);
    check(c4.e_tot == c4b.e_tot && maxAbsDiff(c4.q, c4b.q) == 0.0, "pcg 4 threads: run-to-run bitwise", 0.0, 1e-300);
}

// ---------------------------------------------------------------------------------------------
// (b) alpha independence
static void testAlphaIndependence(const Mol& box, const char* label)
{
    std::printf("(b) alpha independence: %s, %d atoms, charge %d\n", label, box.m_number_atoms, box.m_charge);
    const double alphas[3] = { 0.25, 0.32, 0.42 }; // r_cut 19 Bohr: erfc(0.25*19) = 1.8e-11
    Run r[3];
    for (int t = 0; t < 3; ++t) {
        r[t] = runGFNFF(box, flags("cholesky", alphas[t], 19.0), true);
        if (!r[t].ok) {
            check(false, "run", 1, 0);
            return;
        }
        std::printf("  alpha %.2f: E_tot %.12f  E_coul %.12f  N_k %d\n", alphas[t], r[t].e_tot, r[t].e_coul,
            r[t].ff->getWorkspace()->ewald()->numberOfKVectors());
    }
    double dE = 0, dG = 0, dQ = 0;
    const double gmax = r[0].grad.cwiseAbs().maxCoeff();
    for (int t = 1; t < 3; ++t) {
        dE = std::max(dE, relDiff(r[0].e_tot, r[t].e_tot));
        dG = std::max(dG, maxAbsDiff(r[0].grad, r[t].grad) / gmax);
        dQ = std::max(dQ, maxAbsDiff(r[0].q, r[t].q));
    }
    check(dQ < 1e-9, "charge spread over 3 alphas (e)", dQ, 1e-9);
    check(dE < 1e-10, "total energy spread rel", dE, 1e-10);
    check(dG < 1e-8, "total gradient spread (max abs / max |g|)", dG, 1e-8);
}

// ---------------------------------------------------------------------------------------------
// (c) large-cell limit
static void testLargeCellLimit(const Mol& mol, const char* label)
{
    std::printf("(c) large-cell limit: %s (%d atoms)\n", label, mol.m_number_atoms);
    json plain = json::object();
    const Run ref = runGFNFF(mol, plain, true);
    if (!ref.ok) {
        check(false, "reference run", 1, 0);
        return;
    }
    Eigen::Vector3d dip = Eigen::Vector3d::Zero();
    for (int i = 0; i < mol.m_number_atoms; ++i)
        dip += ref.q(i) * mol.m_geometry.row(i).transpose() * A2B;
    std::printf("  non-periodic E_tot %.12f Eh, |dipole| %.4f e Bohr\n", ref.e_tot, dip.norm());
    const double Ls[2] = { 200.0, 300.0 };
    double dq[2], de[2], dg[2];
    for (int t = 0; t < 2; ++t) {
        const Mol cell = withCell(mol, Ls[t] / A2B, Ls[t] / A2B, Ls[t] / A2B);
        json f = flags("cholesky", 0.0, 0.0, 1e-12);
        f["dispersion_tail"] = false; // the pair term alone equals the cluster (the tail is -1e-10 Eh here)
        const Run r = runGFNFF(cell, f, true);
        if (!r.ok) {
            check(false, "periodic run", 1, 0);
            return;
        }
        dq[t] = maxAbsDiff(r.q, ref.q);
        de[t] = r.e_tot - ref.e_tot;
        dg[t] = maxAbsDiff(r.grad, ref.grad);
        const double V = std::pow(Ls[t], 3);
        std::printf("  L %4.0f Bohr: max |dq| %.2e e, E(pbc) - E %.3e Eh (dipole term -2 pi M^2/(3V) = %.3e), max |dg| %.2e Eh/Bohr\n",
            Ls[t], dq[t], de[t], -2.0 * PI * dip.squaredNorm() / (3.0 * V), dg[t]);
    }
    // A cluster with a dipole M feels its images under the tin-foil boundary condition: -2 pi M^2/(3V) in the energy
    // and the corresponding 1/L^3 shifts in charges and gradient (physical, not an error of the sum). The checks
    // remove the dipole term from the energy and require the 1/L^3 scaling of the rest.
    const double dip_term[2] = { -2.0 * PI * dip.squaredNorm() / (3.0 * std::pow(Ls[0], 3)), -2.0 * PI * dip.squaredNorm() / (3.0 * std::pow(Ls[1], 3)) };
    check(std::abs(de[1] - dip_term[1]) < 1e-8, "|E(pbc, L = 300) - E - dipole term| (Eh)", std::abs(de[1] - dip_term[1]), 1e-8);
    check(dq[1] < 2e-6, "charges at L = 300 equal the non-periodic ones (e)", dq[1], 2e-6);
    check(dg[1] < 2e-7, "gradient at L = 300 equal the non-periodic one (Eh/Bohr)", dg[1], 2e-7);
    const double expo_q = std::log(dq[0] / dq[1]) / std::log(Ls[1] / Ls[0]);
    std::printf("  scaling exponent of max |dq| with L: %.2f (3 for the dipole image term)\n", expo_q);
    check(dq[1] < dq[0] && std::abs(de[1]) < std::abs(de[0]) && (dq[1] < 1e-8 || expo_q > 2.5), "deviations decrease with L (1/L^3 for a dipole)", expo_q, 2.5);
}

// ---------------------------------------------------------------------------------------------
// (d) invariances with the full periodic EEQ
static void testInvariances(const Mol& box, const char* label)
{
    std::printf("(d) translation, atom-wise wrap, 2x1x1 supercell: %s (%d atoms, charge %d)\n", label, box.m_number_atoms, box.m_charge);
    // eps 1e-12: the supercell has other k vectors and pairs than the cell, so the two sums differ by the
    // Ewald truncation error (3e-10 relative at eps 1e-10 on the charged box), not by the EEQ.
    json f = flags("cholesky", 0.0, 0.0, 1e-12);
    f["hb_min_pair_energy_eh"] = 0.0; // the |E_HB| < 1e-9 candidate filter is the only term that differs in a supercell (4e-11 Eh)
    const Run base = runGFNFF(box, f, true);
    if (!base.ok) {
        check(false, "base run", 1, 0);
        return;
    }
    const int N = box.m_number_atoms;
    // With the filter on, the charged box differed by 1.8e-9 Eh (one dropped HB candidate); with it off both boxes,
    // neutral and Q = +1 (ensemble charge model, two hydroniums in the supercell), agree to 1e-15 (--supercell-terms).
    const double e_lim = 1e-12, g_lim = 1e-12;
    const double L = box.m_unit_cell(0, 0);
    Mol shifted = box;
    for (int i = 0; i < N; ++i) {
        shifted.m_geometry(i, 0) += 0.37 * L;
        shifted.m_geometry(i, 1) -= 1.21 * L;
        shifted.m_geometry(i, 2) += 0.05 * L;
    }
    const Run s = runGFNFF(shifted, f, true);
    check(maxAbsDiff(base.q, s.q) < 1e-11, "translation: charges (e)", maxAbsDiff(base.q, s.q), 1e-11);
    check(relDiff(base.e_tot, s.e_tot) < 1e-12, "translation: energy rel", relDiff(base.e_tot, s.e_tot), 1e-12);
    check(maxAbsDiff(base.grad, s.grad) < 1e-11, "translation: gradient (Eh/Bohr)", maxAbsDiff(base.grad, s.grad), 1e-11);

    Mol wrapped = shifted;
    int torn = 0;
    for (int i = 0; i < N; ++i)
        for (int a = 0; a < 3; ++a) {
            const double w = wrapped.m_geometry(i, a) - L * std::floor(wrapped.m_geometry(i, a) / L);
            if (i % 3 != 0 && std::abs(w - (wrapped.m_geometry(i - i % 3, a) - L * std::floor(wrapped.m_geometry(i - i % 3, a) / L))) > L / 2)
                ++torn;
            wrapped.m_geometry(i, a) = w;
        }
    const Run w = runGFNFF(wrapped, f, true);
    std::printf("  atom-wise wrap: %d atom coordinates torn from their oxygen\n", torn);
    check(maxAbsDiff(base.q, w.q) < 1e-11, "atom-wise wrap: charges (e)", maxAbsDiff(base.q, w.q), 1e-11);
    check(relDiff(base.e_tot, w.e_tot) < 1e-12, "atom-wise wrap: energy rel", relDiff(base.e_tot, w.e_tot), 1e-12);
    check(maxAbsDiff(base.grad, w.grad) < 1e-11, "atom-wise wrap: gradient (Eh/Bohr)", maxAbsDiff(base.grad, w.grad), 1e-11);

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
    const Run sc = runGFNFF(super, f, true);
    const double dq = std::max(maxAbsDiff(base.q, sc.q.head(N)), maxAbsDiff(base.q, sc.q.tail(N)));
    std::printf("  E_tot cell %.12f  supercell/2 %.12f\n", base.e_tot, sc.e_tot / 2);
    check(dq < 1e-10, "2x1x1 supercell: charges per atom (e)", dq, 1e-10);
    check(relDiff(base.e_tot, sc.e_tot / 2) < e_lim, "2x1x1 supercell: energy = 2 x cell rel", relDiff(base.e_tot, sc.e_tot / 2), e_lim);
    check(std::max(maxAbsDiff(base.grad, sc.grad.topRows(N)), maxAbsDiff(base.grad, sc.grad.bottomRows(N))) < g_lim,
        "2x1x1 supercell: gradient per copy (Eh/Bohr)", std::max(maxAbsDiff(base.grad, sc.grad.topRows(N)), maxAbsDiff(base.grad, sc.grad.bottomRows(N))), g_lim);
}

// ---------------------------------------------------------------------------------------------
// (e) the gate: analytic total gradient vs central FD of the total energy with re-solved charges
static Matrix fdTotal(GFNFF& ff, const Mol& mol, double h_ang)
{
    Geometry geom = mol.m_geometry;
    Matrix fd = Matrix::Zero(mol.m_number_atoms, 3);
    auto energy_at = [&]() {
        ff.UpdateMolecule(geom);
        return ff.Calculation(false);
    };
    for (int i = 0; i < mol.m_number_atoms; ++i)
        for (int a = 0; a < 3; ++a) {
            geom(i, a) += h_ang;
            const double ep = energy_at();
            geom(i, a) -= 2 * h_ang;
            const double em = energy_at();
            geom(i, a) += h_ang;
            fd(i, a) = (ep - em) / (2 * h_ang); // Eh/Angstrom
        }
    energy_at();
    return fd;
}

static void testFiniteDifferences(const Mol& box, const json& solver_flags, const char* label)
{
    std::printf("(e) total gradient vs FD of the total energy (charges re-solved): %s (%d atoms, charge %d)\n", label, box.m_number_atoms, box.m_charge);
    json f = solver_flags;
    f["hb_update_force_every"] = 1;   // H-bond lists follow every displaced geometry
    f["hb_min_pair_energy_eh"] = 0.0; // no candidate pops in or out of the list by its energy filter
    Run r = runGFNFF(box, f, true);
    if (!r.ok) {
        check(false, "run", 1, 0);
        return;
    }
    const Matrix g = r.grad * A2B; // Eh/Bohr -> Eh/Angstrom
    const double gmax = g.cwiseAbs().maxCoeff();
    const double L = box.m_unit_cell(0, 0);
    double min_face = L;
    for (int i = 0; i < box.m_number_atoms; ++i)
        for (int a = 0; a < 3; ++a) {
            const double x = box.m_geometry(i, a) - L * std::floor(box.m_geometry(i, a) / L);
            min_face = std::min(min_face, std::min(x, L - x));
        }
    std::printf("  path '%s', max |g| %.3e Eh/A, closest atom to a face %.3f A\n", r.path.c_str(), gmax, min_face);
    const double hs[2] = { 1e-4, 2e-4 };
    Matrix fd[2];
    double err[2];
    Eigen::Index wi[2], wj[2];
    for (int t = 0; t < 2; ++t) {
        fd[t] = fdTotal(*r.ff, box, hs[t]);
        err[t] = (fd[t] - g).cwiseAbs().maxCoeff(&wi[t], &wj[t]);
        std::printf("  h = %.0e A: max |FD - analytic| = %.3e Eh/A at atom %d (xyz %d)\n", hs[t], err[t], int(wi[t]), int(wj[t]));
    }
    const Matrix rich = (4.0 * fd[0] - fd[1]) / 3.0;
    const double err_r = (rich - g).cwiseAbs().maxCoeff();
    std::printf("  ratio err(2h)/err(h) = %.2f; Richardson max |FD - analytic| = %.3e Eh/A\n", err[1] / err[0], err_r);
    check(err[0] < 1e-7, "max |FD - analytic| at h = 1e-4 A (Eh/A)", err[0], 1e-7);
    check(err[1] / err[0] > 3.0 && err[1] / err[0] < 5.0, "quadratic in h (ratio 3..5)", err[1] / err[0], 5.0);
    check(err_r < 1e-8, "Richardson-extrapolated max |FD - analytic| (Eh/A)", err_r, 1e-8);
}

// ---------------------------------------------------------------------------------------------
// (f) MIC-only EEQ against the Ewald EEQ (the measurement planned as B1)
static void testMicOnlyEEQ(const Mol& box, const char* label)
{
    std::printf("(f) MIC-only EEQ vs Ewald EEQ: %s (%d atoms, %.4f atoms/A^3)\n", label, box.m_number_atoms,
        box.m_number_atoms / std::pow(box.m_unit_cell(0, 0), 3));
    const Run r = runGFNFF(box, flags("cholesky", 0.0, 0.0, 1e-12), false);
    if (!r.ok) {
        check(false, "run", 1, 0);
        return;
    }
    const FFWorkspace* ws = r.ff->getWorkspace();
    const int N = box.m_number_atoms;
    const Vector& alp = ws->coulombAlp();
    const Vector& gam = ws->coulombGam();
    const Vector& chi_base = ws->coulombChiBase();
    const Vector& cnf = ws->coulombCnf();
    const Vector& cn = r.ff->getLastCN();
    Vector chi(N);
    for (int i = 0; i < N; ++i)
        chi(i) = chi_base(i) + cnf(i) * std::sqrt(std::max(cn(i), 0.0));
    const Eigen::MatrixXd pos = box.m_geometry * A2B;
    const PBCUtils::OrthoCell& cell = ws->periodicCell();
    const double TSQRT2PI = 0.797884560802866;
    // Both matrices with the same local diagonal and the instance's own fragment constraints (a dense packing can
    // merge waters into one fragment: 970 waters at 1.7 A gave 939 fragments, so "one row per water" would be wrong).
    const auto& topo = r.ff->getTopologyInfo();
    const int nfrag = topo.nfrag;
    std::printf("  %d fragments (instance topology)\n", nfrag);
    auto solveAug = [&](const Eigen::MatrixXd& Ann) {
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N + nfrag, N + nfrag);
        A.topLeftCorner(N, N) = Ann;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(N + nfrag);
        b.head(N) = chi;
        for (int i = 0; i < N; ++i) {
            A(N + topo.fraglist[i] - 1, i) = 1.0;
            A(i, N + topo.fraglist[i] - 1) = 1.0;
        }
        for (int f = 0; f < nfrag; ++f)
            b(N + f) = f < static_cast<int>(topo.qfrag.size()) ? topo.qfrag[f] : 0.0;
        const Eigen::VectorXd sol = A.partialPivLu().solve(b);
        return Vector(sol.head(N));
    };
    Eigen::MatrixXd A_ew = ws->ewald()->matrix(pos, &alp);
    Eigen::MatrixXd A_mic = Eigen::MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        A_ew(i, i) += gam(i) + TSQRT2PI / std::sqrt(alp(i));
        A_mic(i, i) = gam(i) + TSQRT2PI / std::sqrt(alp(i));
        for (int j = 0; j < i; ++j) {
            double dx = pos(j, 0) - pos(i, 0), dy = pos(j, 1) - pos(i, 1), dz = pos(j, 2) - pos(i, 2);
            PBCUtils::applyMinimumImageOrtho(cell, dx, dy, dz);
            const double rr = std::sqrt(dx * dx + dy * dy + dz * dz);
            A_mic(i, j) = A_mic(j, i) = std::erf(rr / std::sqrt(alp(i) + alp(j))) / rr;
        }
    }
    const Vector q_ew = solveAug(A_ew), q_mic = solveAug(A_mic);
    const double e_ew = 0.5 * q_ew.dot(A_ew * q_ew) - chi.dot(q_ew);
    const double e_mic = 0.5 * q_mic.dot(A_mic * q_mic) - chi.dot(q_mic);
    const double e_mic_in_ew = 0.5 * q_mic.dot(A_ew * q_mic) - chi.dot(q_mic);
    const double consistency = maxAbsDiff(q_ew, r.q);
    std::printf("  test-built dense Ewald EEQ vs instance charges: %.2e e (checks the inputs used here)\n", consistency);
    std::printf("  max |q_mic - q_ew| = %.3e e (max |q| %.3f), mean |dq| %.3e e\n", maxAbsDiff(q_mic, q_ew), q_ew.cwiseAbs().maxCoeff(),
        (q_mic - q_ew).cwiseAbs().mean());
    std::printf("  E_EEQ(Ewald, q_ew) %.8f  E_EEQ(MIC, q_mic) %.8f  difference %.3e Eh; E_EEQ(Ewald, q_mic) - E_EEQ(Ewald, q_ew) %.3e Eh\n",
        e_ew, e_mic, e_mic - e_ew, e_mic_in_ew - e_ew);
    // The workspace Coulomb component equals E_EEQ(q) exactly (pair + self + chi terms), every other term is charge
    // independent (Phase-1 charges), so E_tot(MIC) - E_tot(Ewald) = E_EEQ(MIC, q_mic) - E_EEQ(Ewald, q_ew).
    std::printf("  E_tot(Ewald) %.8f Eh, E_coul component %.8f (= E_EEQ(Ewald, q_ew) to %.1e); E_tot(MIC) - E_tot(Ewald) = %.3e Eh\n",
        r.e_tot, r.e_coul, std::abs(r.e_coul - e_ew), e_mic - e_ew);
    check(consistency < 1e-8, "test-built Ewald EEQ reproduces the instance charges (e)", consistency, 1e-8);
    check(e_mic_in_ew >= e_ew - 1e-12, "Ewald charges minimise the Ewald functional (MIC charges lie above)", e_mic_in_ew - e_ew, 1e-12);
}

// ---------------------------------------------------------------------------------------------
// (g0) no cell: default (true) and false bitwise equal
static void testNoCell()
{
    std::printf("(g0) no cell: coulomb_ewald default == false bitwise\n");
    const Mol mol = registryMol("H2O_dimer");
    json off = json::object();
    off["coulomb_ewald"] = false;
    const Run a = runGFNFF(mol, json::object(), true), b = runGFNFF(mol, off, true);
    check(a.ok && b.ok && a.e_tot == b.e_tot && a.e_coul == b.e_coul && maxAbsDiff(a.grad, b.grad) == 0.0 && maxAbsDiff(a.q, b.q) == 0.0,
        "energy, charges and gradient bitwise", std::abs(a.e_tot - b.e_tot), 1e-300);
    check(a.ok && !a.ff->periodicEEQActive() && !a.ff->getWorkspace()->ewaldActive(), "Ewald inactive without a cell", 0.0, 1e-300);
}

// ---------------------------------------------------------------------------------------------
// perf mode
static void perfOne(const Mol& mol, const json& f, int threads, const char* label)
{
    Run r = runGFNFF(mol, f, true, threads);
    if (!r.ok) {
        std::printf("  %s: failed\n", label);
        return;
    }
    FFWorkspace* ws = r.ff->getWorkspace();
    double wall[3], coul[3], eeq[3];
    int iters = 0;
    Geometry geom = mol.m_geometry;
    for (int k = 0; k < 3; ++k) {
        geom(0, 0) += 1e-3; // a new geometry: no Phase-2 skip
        r.ff->UpdateMolecule(geom);
        auto t0 = std::chrono::steady_clock::now();
        r.ff->Calculation(true);
        wall[k] = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        coul[k] = ws->termTimings().coulomb + ws->ewaldSetupMs();
        eeq[k] = r.ff->eeqSolver() ? r.ff->eeqSolver()->lastPeriodicSolveMs() : 0.0;
        iters = r.ff->eeqSolver() ? r.ff->eeqSolver()->lastPeriodicIterations() : 0;
    }
    std::sort(wall, wall + 3);
    std::sort(coul, coul + 3);
    std::sort(eeq, eeq + 3);
    if (ws->ewaldActive())
        std::printf("  %-28s %d thr: r_cut %5.2f N_k %6d | EEQ %s %4d it %8.1f ms | Coulomb CPU-sum %8.1f ms | Calculation(true) wall %8.1f ms | E %.8f\n",
            label, threads, ws->ewald()->parameters().r_cut, ws->ewald()->numberOfKVectors(), r.path.c_str(), iters, eeq[1], coul[1], wall[1], r.e_tot);
    else
        std::printf("  %-28s %d thr: %22s | EEQ (non-periodic, not timed here) | Coulomb CPU-sum %8.1f ms | Calculation(true) wall %8.1f ms | E %.8f\n",
            label, threads, "plain", coul[1], wall[1], r.e_tot);
}

static void perfBox(int nwat, double L_bohr, int threads)
{
    const Mol box = waterBox(nwat, L_bohr, 20261107, {}, 1.7);
    std::printf("perf: random water box, %d atoms in L = %.1f Bohr (density %.4f atoms/A^3), %d threads\n", box.m_number_atoms, L_bohr,
        box.m_number_atoms / std::pow(L_bohr / A2B, 3), threads);
    Mol nocell = box;
    nocell.m_has_pbc = false;
    perfOne(nocell, json::object(), threads, "no cell");
    json pcg = flags("pcg", 0.0, 0.0, 1e-8);
    perfOne(box, pcg, threads, "cell, pcg, r_cut auto");
    json pcg25 = flags("pcg", 0.0, 25.0, 1e-8);
    perfOne(box, pcg25, threads, "cell, pcg, r_cut 25");
    if (box.m_number_atoms <= 1500)
        perfOne(box, flags("cholesky", 0.0, 0.0, 1e-8), threads, "cell, dense cholesky");
}

static void supercellTerms(const Mol& box, const json& f, const char* label)
{
    const int N = box.m_number_atoms;
    const double L = box.m_unit_cell(0, 0);
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
    const Run a = runGFNFF(box, f, true), b = runGFNFF(super, f, true);
    const auto& ea = a.ff->getWorkspace()->energyComponents();
    const auto& eb = b.ff->getWorkspace()->energyComponents();
    std::printf("%s: term  cell  supercell/2  difference\n", label);
    auto line = [&](const char* n, double x, double y) { std::printf("  %-14s %18.12f %18.12f %12.3e\n", n, x, y / 2, x - y / 2); };
    line("total", a.e_tot, b.e_tot);
    line("bond", ea.bond, eb.bond);
    line("angle", ea.angle, eb.angle);
    line("dihedral", ea.dihedral, eb.dihedral);
    line("inversion", ea.inversion, eb.inversion);
    line("bonded_rep", ea.bonded_rep, eb.bonded_rep);
    line("nonbonded_rep", ea.nonbonded_rep, eb.nonbonded_rep);
    line("dispersion", ea.dispersion, eb.dispersion);
    line("disp_tail", ea.disp_tail, eb.disp_tail);
    line("coulomb", ea.coulomb, eb.coulomb);
    line("hbond", ea.hbond, eb.hbond);
    line("xbond", ea.xbond, eb.xbond);
    line("batm", ea.batm, eb.batm);
    const auto& ta = a.ff->getWorkspace()->ewaldTerms();
    const auto& tb = b.ff->getWorkspace()->ewaldTerms();
    line("ewald real", ta.real, tb.real);
    line("ewald recip", ta.reciprocal, tb.reciprocal);
    line("ewald self", ta.self, tb.self);
    line("ewald bg", ta.background, tb.background);
    line("ewald gauss_bg", ta.gaussian_background, tb.gaussian_background);
    std::printf("  charges per atom max diff %.2e e; sum q cell %.12f supercell %.12f\n",
        std::max(maxAbsDiff(a.q, b.q.head(N)), maxAbsDiff(a.q, b.q.tail(N))), a.q.sum(), b.q.sum());
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc > 3 && std::strcmp(argv[1], "--mic-vs-ewald") == 0) {
        // MIC-only EEQ against the Ewald EEQ at a chosen size and density (random rigid packing, NOT an equilibrated liquid)
        char label[96];
        std::snprintf(label, sizeof(label), "%s H2O, L = %s Bohr, random packing", argv[2], argv[3]);
        testMicOnlyEEQ(waterBox(std::atoi(argv[2]), std::atof(argv[3]), 20261107, {}, 1.7), label);
        return g_failures == 0 ? 0 : 1;
    }
    if (argc > 1 && std::strcmp(argv[1], "--supercell-terms") == 0) {
        const Mol charged = waterBox(18, 44.0, 20261106, { registryMol("hydronium_water2", 1) });
        json f = flags("cholesky", 0.0, 0.0, 1e-12);
        supercellTerms(charged, f, "Q = +1, frag_charge_model ensemble (default)");
        f["frag_charge_model"] = "reference";
        supercellTerms(charged, f, "Q = +1, frag_charge_model reference");
        const Mol box = waterBox(20, 44.0, 20261105);
        supercellTerms(box, flags("cholesky", 0.0, 0.0, 1e-12), "neutral 20 H2O");
        return 0;
    }
    if (argc > 3 && std::strcmp(argv[1], "--perf-box") == 0) {
        perfBox(std::atoi(argv[2]), std::atof(argv[3]), argc > 4 ? std::atoi(argv[4]) : 4);
        return 0;
    }
    testNoCell();
    const Mol box = waterBox(20, 44.0, 20261105);
    const Mol charged = waterBox(18, 44.0, 20261106, { registryMol("hydronium_water2", 1) });
    testDenseVsPcg(box, "20 H2O, L = 44 Bohr");
    testAlphaIndependence(box, "20 H2O, L = 44 Bohr");
    testAlphaIndependence(charged, "H3O+ 2 H2O + 18 H2O, L = 44 Bohr (Q = +1)");
    testLargeCellLimit(registryMol("H2O_dimer"), "water dimer");
    testLargeCellLimit(registryMol("C6H6"), "benzene");
    {
        Mol cluster = waterBox(8, 16.0, 20261108, {}, 2.4);
        cluster.m_has_pbc = false;
        testLargeCellLimit(cluster, "8-water cluster (random packing in 8.5 A)");
    }
    testInvariances(box, "20 H2O");
    testInvariances(charged, "H3O+ 2 H2O + 18 H2O (Q = +1)");
    testFiniteDifferences(box, flags("cholesky", 0.0, 0.0, 1e-12), "20 H2O, dense cholesky");
    {
        json pcg = flags("pcg", 0.0, 0.0, 1e-12);
        pcg["eeq_ppcg_tol"] = 1e-14;
        testFiniteDifferences(box, pcg, "20 H2O, matrix-free pcg (tol 1e-14)");
    }
    testFiniteDifferences(charged, flags("cholesky", 0.0, 0.0, 1e-12), "H3O+ 2 H2O + 18 H2O (Q = +1), dense cholesky");
    testMicOnlyEEQ(waterBox(50, 44.0, 20261109, {}, 2.4), "50 H2O, L = 44 Bohr");

    std::printf("test_pbc_eeq_ewald: %d failure(s)\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
