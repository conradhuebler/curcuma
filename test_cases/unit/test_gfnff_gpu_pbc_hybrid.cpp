/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * S2 of the periodic GPU port (Claude Generated, Nov 2026): the CUDA hybrid with an orthorhombic cell and
 * coulomb_ewald (default): local terms with the minimum image on the device (S1), periodic EEQ (Ewald matrix, P1 pcg)
 * and the Ewald Coulomb pair term on the host, against the CPU periodic path.
 *
 * Boxes (as in pbc_eeq_ewald, rigid random packings, not liquids): 20 and 50 waters in L = 44 Bohr, H3O+(H2O)2 +
 * 18 H2O (Q = +1, L = 44 Bohr, frag_charge_model reference so that the GPU wrapper keeps the device), 970 waters in
 * L = 58 Bohr (2910 atoms, pcg above eeq_ewald_dense_max_atoms).
 *
 * Checks (thresholds fixed beforehand):
 *  (i)   GPU hybrid vs CPU: charges |dq| <= 1e-12 e; at the default tolerance |dE| <= 1e-9 Eh per system; with
 *        eeq_ewald_tol 1e-14 on both sides |dE| <= 1e-11 Eh and gradient max |dG| <= 1e-12 Eh/A.
 *  (ii)  FD gate on the hybrid: analytic total gradient vs central FD of the total energy with re-solved charges,
 *        h = 1e-4 / 2e-4 A (HB lists every call, candidate filter off), Richardson <= 1e-8 Eh/A, neutral and Q = +1.
 *  (iii) invariances of the hybrid: translation by (L,-L,2L), atom-wise wrap (torn), 2x1x1 supercell: charges per atom
 *        <= 1e-12 e, energy <= 1e-11 Eh (supercell: 2 x cell, relative 1e-12), gradient per atom <= 1e-12 Eh/A.
 *  (iv)  two hybrid instances (1 thread, and 4 threads): energy and gradient bit for bit identical.
 * Exit 77 (ctest SKIP_RETURN_CODE) when built without CUDA.
 */
#ifdef USE_CUDA
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/qm_methods/gfnff_gpu_method.h"
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}
static std::string sci(double v)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.2e", v);
    return b;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;

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

/// Rigid random packing of copies of `units` into a cube of L (Angstrom), minimum-image distance >= min_any between
/// atoms of different molecules; units[0] first. The same construction as test_pbc_eeq_ewald.cpp (same seeds give the
/// same boxes).
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

static std::vector<std::vector<int>> fragmentsOf(const Mol& m)
{
    curcuma::Molecule mol(m);
    return mol.GetFragments();
}

static json config(int threads, const json& extra)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["gpu_periodic"] = true;          // the S2 opt-in; coulomb_ewald stays at its default (true)
    c["hb_update_force_every"] = 1;    // HB/XB lists from the geometry of every call
    c["hb_min_pair_energy_eh"] = 0.0;  // no candidate pops in or out by its energy filter (supercell, FD)
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    return c;
}

struct Result {
    double e = 0.0;
    double e_ewald = 0.0;   // GPU: host Ewald pair part; CPU: workspace Ewald pair energy
    Matrix grad;            // Eh/A
    Vector q;
};

/// One CPU or GPU-hybrid instance through the ComputationalMethod wrappers (same gradient contract Eh/A).
struct Instance {
    std::unique_ptr<GFNFFComputationalMethod> cpu;
    std::unique_ptr<GFNFFGPUComputationalMethod> gpu;
    bool ok = false;
    Instance(bool on_gpu, const Mol& mol, const json& c)
    {
        try {
            if (on_gpu) {
                gpu = std::make_unique<GFNFFGPUComputationalMethod>("gfnff", c);
                ok = gpu->setMolecule(mol) && !gpu->hasError() && !gpu->requestsCpuFallback() && gpu->periodicHybridActive()
                    && gpu->hasPeriodicEnergy();
                if (!ok) std::printf("GPU init: %s%s hybrid=%d\n", gpu->getErrorMessage().c_str(),
                                     gpu->requestsCpuFallback() ? " (CPU fallback requested)" : "", gpu->periodicHybridActive() ? 1 : 0);
            } else {
                cpu = std::make_unique<GFNFFComputationalMethod>("gfnff", c);
                ok = cpu->setMolecule(mol) && !cpu->hasError() && cpu->hasPeriodicEnergy();
            }
        } catch (const std::exception& e) {
            std::printf("init threw: %s\n", e.what());
            ok = false;
        }
    }
    Result evaluate(const Geometry& g_ang, bool gradient)
    {
        Result r;
        if (gpu) {
            gpu->updateGeometry(g_ang);
            r.e = gpu->calculateEnergy(gradient);
            if (gradient) r.grad = gpu->getGradient();
            r.q = gpu->getCharges();
            r.e_ewald = gpu->hostEwaldEnergy();
        } else {
            cpu->updateGeometry(g_ang);
            r.e = cpu->calculateEnergy(gradient);
            if (gradient) r.grad = cpu->getGradient();
            r.q = cpu->getCharges();
            r.e_ewald = cpu->getGFNFF()->getWorkspace()->ewaldPairEnergy();
        }
        return r;
    }
};

static double maxAbs(const Matrix& a, const Matrix& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols())
        return 1e300;
    return a.rows() ? (a - b).cwiseAbs().maxCoeff() : 0.0;
}
static double maxAbsV(const Vector& a, const Vector& b)
{
    if (a.size() != b.size())
        return 1e300;
    return a.size() ? (a - b).cwiseAbs().maxCoeff() : 0.0;
}

static double g_worst_q = 0.0, g_worst_e_def = 0.0, g_worst_e_tight = 0.0, g_worst_g_tight = 0.0;

/// Every atom moved by a pseudo-random 0.005 A: at the set-up geometry both paths reuse the charges of the topology
/// stage (CPU CN); the moved geometry forces the per-step solve, on the GPU with the device CN (the seam under test).
static Geometry moved(const Geometry& g0, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-0.005, 0.005);
    Geometry g = g0;
    for (int i = 0; i < g.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            g(i, k) += u(rng);
    return g;
}

// (i) GPU hybrid vs CPU periodic, default tolerance and eeq_ewald_tol 1e-14; at the input geometry (charges of the
// set-up solve reused on both sides) and at a moved geometry (per-step solve, device CN on the GPU, aspc warm start)
static void testVsCpu(const char* label, const Mol& mol, const json& extra, int threads)
{
    const Geometry g_moved = moved(mol.m_geometry, 20261104);
    // default tolerance
    {
        Instance gpu(true, mol, config(threads, extra));
        Instance cpu(false, mol, config(threads, extra));
        check(gpu.ok && cpu.ok, std::string(label) + ": set-up (GPU hybrid, CPU)");
        if (gpu.ok && cpu.ok) {
            const Geometry* geos[2] = { &mol.m_geometry, &g_moved };
            const char* what[2] = { "input geometry", "moved geometry" };
            for (int s = 0; s < 2; ++s) {
                const Result g = gpu.evaluate(*geos[s], true), c = cpu.evaluate(*geos[s], true);
                const double dq = maxAbsV(g.q, c.q), de = std::abs(g.e - c.e), dg = maxAbs(g.grad, c.grad);
                g_worst_q = std::max(g_worst_q, dq);
                g_worst_e_def = std::max(g_worst_e_def, de);
                check(dq <= 1e-12, std::string(label) + " default tol, " + what[s] + ": charges GPU vs CPU (" + sci(dq) + " e)");
                check(de <= 1e-9, std::string(label) + " default tol, " + what[s] + ": energy GPU vs CPU (" + sci(de) + " Eh)");
                std::printf("%-34s default tol, %-15s E %.12f / %.12f (dE %.2e), E_ewald %.10f / %.10f, dq %.2e, dG %.2e Eh/A\n", label,
                            what[s], g.e, c.e, de, g.e_ewald, c.e_ewald, dq, dg);
            }
        }
    }
    // tight tolerance on both sides (pcg: 1e-14; the dense path below eeq_ewald_dense_max_atoms is exact anyway)
    {
        json e2 = extra;
        e2["eeq_ewald_tol"] = 1e-14;
        Instance gpu(true, mol, config(threads, e2));
        Instance cpu(false, mol, config(threads, e2));
        check(gpu.ok && cpu.ok, std::string(label) + ": set-up with eeq_ewald_tol 1e-14");
        if (gpu.ok && cpu.ok) {
            const Geometry* geos[2] = { &mol.m_geometry, &g_moved };
            const char* what[2] = { "input geometry", "moved geometry" };
            for (int s = 0; s < 2; ++s) {
                const Result g = gpu.evaluate(*geos[s], true), c = cpu.evaluate(*geos[s], true);
                const double dq = maxAbsV(g.q, c.q), de = std::abs(g.e - c.e), dg = maxAbs(g.grad, c.grad);
                g_worst_q = std::max(g_worst_q, dq);
                g_worst_e_tight = std::max(g_worst_e_tight, de);
                g_worst_g_tight = std::max(g_worst_g_tight, dg);
                check(dq <= 1e-12, std::string(label) + " tol 1e-14, " + what[s] + ": charges GPU vs CPU (" + sci(dq) + " e)");
                check(de <= 1e-11, std::string(label) + " tol 1e-14, " + what[s] + ": energy GPU vs CPU (" + sci(de) + " Eh)");
                check(dg <= 1e-12, std::string(label) + " tol 1e-14, " + what[s] + ": gradient GPU vs CPU (" + sci(dg) + " Eh/A)");
                check(std::abs(g.e_ewald - c.e_ewald) <= 1e-11, std::string(label) + " tol 1e-14, " + what[s] + ": Ewald pair part GPU vs CPU ("
                                                                     + sci(std::abs(g.e_ewald - c.e_ewald)) + " Eh)");
                std::printf("%-34s tol 1e-14,   %-15s dE %.2e Eh, dq %.2e e, dG %.2e Eh/A, max |g| %.2e, dE_ewald %.2e\n", label,
                            what[s], de, dq, dg, g.grad.cwiseAbs().maxCoeff(), std::abs(g.e_ewald - c.e_ewald));
            }
        }
    }
}

// (ii) FD gate of the hybrid total energy with re-solved charges
static void testFiniteDifference(const char* label, const Mol& mol, const json& extra, int nfd)
{
    json e2 = extra;
    e2["eeq_ewald_tol"] = 1e-14;
    Instance gpu(true, mol, config(1, e2));
    check(gpu.ok, std::string(label) + ": FD set-up");
    if (!gpu.ok)
        return;
    const Result r = gpu.evaluate(mol.m_geometry, true);
    const int n = mol.m_number_atoms;
    // atoms nearest a face first (the minimum-image seams), then the rest
    const double L = mol.m_unit_cell(0, 0);
    std::vector<std::pair<double, int>> dist;
    for (int i = 0; i < n; ++i) {
        double dmin = 1e300;
        for (int k = 0; k < 3; ++k) {
            const double x = mol.m_geometry(i, k) - L * std::floor(mol.m_geometry(i, k) / L);
            dmin = std::min({ dmin, x, L - x });
        }
        dist.push_back({ dmin, i });
    }
    std::sort(dist.begin(), dist.end());
    const double hs[2] = { 2e-4, 1e-4 };
    double err[2] = { 0, 0 }, err_r = 0.0, gmax = 0.0;
    for (int s = 0; s < nfd && s < n; ++s) {
        const int a = dist[s].second;
        for (int k = 0; k < 3; ++k) {
            double fd[2];
            for (int hi = 0; hi < 2; ++hi) {
                Geometry gp = mol.m_geometry, gm = mol.m_geometry;
                gp(a, k) += hs[hi];
                gm(a, k) -= hs[hi];
                const double ep = gpu.evaluate(gp, false).e, em = gpu.evaluate(gm, false).e;
                fd[hi] = (ep - em) / (2.0 * hs[hi]);
                err[hi] = std::max(err[hi], std::abs(fd[hi] - r.grad(a, k)));
            }
            err_r = std::max(err_r, std::abs((4.0 * fd[1] - fd[0]) / 3.0 - r.grad(a, k)));
            gmax = std::max(gmax, std::abs(r.grad(a, k)));
        }
    }
    check(err[1] <= 1e-7, std::string(label) + " FD: analytic vs central difference at h 1e-4 A (" + sci(err[1]) + " Eh/A)");
    check(err_r <= 1e-8, std::string(label) + " FD: Richardson-extrapolated difference (" + sci(err_r) + " Eh/A)");
    check(err[0] > 2.5 * err[1] || err[0] < 1e-10, std::string(label) + " FD: error quadratic in h (" + sci(err[0]) + " -> " + sci(err[1]) + ")");
    std::printf("%-34s FD (%d atoms nearest a face): h 2e-4 %.2e, h 1e-4 %.2e, Richardson %.2e Eh/A, max |g| %.2e\n", label,
                std::min(nfd, n), err[0], err[1], err_r, gmax);
}

// (iii) invariances of the hybrid
static void testInvariance(const char* label, const Mol& mol, const json& extra, bool do_supercell)
{
    json e2 = extra;
    e2["eeq_ewald_tol"] = 1e-14;
    const int n = mol.m_number_atoms;
    const double L = mol.m_unit_cell(0, 0);
    Instance ref(true, mol, config(1, e2));
    check(ref.ok, std::string(label) + ": invariance set-up");
    if (!ref.ok)
        return;
    const Result r0 = ref.evaluate(mol.m_geometry, true);

    auto compare = [&](const Result& r, const std::string& what, double scale_e) {
        const double dq = maxAbsV(r.q, r0.q), de = std::abs(r.e - scale_e * r0.e), dg = maxAbs(r.grad, r0.grad);
        const double tol_e = std::max(1e-11, 1e-12 * std::abs(r.e));
        check(dq <= 1e-12 && de <= tol_e && dg <= 1e-12, std::string(label) + " " + what + ": charges " + sci(dq) + " e, energy "
                                                              + sci(de) + " Eh, gradient " + sci(dg) + " Eh/A");
        std::printf("%-34s %-28s dq %.2e e, dE %.2e Eh, dG %.2e Eh/A\n", label, what.c_str(), dq, de, dg);
    };
    // translation by (L,-L,2L)
    {
        Geometry g = mol.m_geometry;
        for (int i = 0; i < n; ++i) {
            g(i, 0) += L;
            g(i, 1) -= L;
            g(i, 2) += 2.0 * L;
        }
        compare(ref.evaluate(g, true), "translated (L,-L,2L)", 1.0);
    }
    // atom-wise wrap (torn molecules), on the same instance and on one initialised torn (A4c on the host)
    {
        Geometry g = mol.m_geometry;
        int torn = 0;
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) {
                const double w = g(i, k) - L * std::floor(g(i, k) / L);
                if (w != g(i, k)) ++torn;
                g(i, k) = w;
            }
        compare(ref.evaluate(g, true), "torn (" + std::to_string(torn) + " coordinates)", 1.0);
        Mol mt = mol;
        mt.m_geometry = g;
        Instance it(true, mt, config(1, e2));
        check(it.ok, std::string(label) + ": torn input set up directly");
        if (it.ok)
            compare(it.evaluate(g, true), "torn, initialised torn", 1.0);
    }
    // 2x1x1 supercell: charges per atom, energy doubled, gradient per atom
    if (do_supercell) {
        Mol m2 = mol;
        m2.m_unit_cell = Eigen::Vector3d(2.0 * L, L, L).asDiagonal();
        m2.m_geometry = Geometry(2 * n, 3);
        m2.m_geometry.topRows(n) = mol.m_geometry;
        m2.m_geometry.bottomRows(n) = mol.m_geometry;
        for (int i = 0; i < n; ++i)
            m2.m_geometry(n + i, 0) += L;
        for (int i = 0; i < n; ++i)
            m2.m_atoms.push_back(mol.m_atoms[i]);
        m2.m_number_atoms = 2 * n;
        m2.m_charge = 2 * mol.m_charge;
        Instance i2(true, m2, config(1, e2));
        check(i2.ok, std::string(label) + ": supercell set-up");
        if (i2.ok) {
            const Result r2 = i2.evaluate(m2.m_geometry, true);
            Result top, bot;
            top.e = bot.e = r2.e;
            top.q = r2.q.head(n);
            bot.q = r2.q.tail(n);
            top.grad = r2.grad.topRows(n);
            bot.grad = r2.grad.bottomRows(n);
            compare(top, "2x1x1 supercell (first copy)", 2.0);
            compare(bot, "2x1x1 supercell (second copy)", 2.0);
        }
    }
}

// (iv) two hybrid instances bit for bit
static void testBitwise(const char* label, const Mol& mol, const json& extra, int threads)
{
    Instance a(true, mol, config(threads, extra)), b(true, mol, config(threads, extra));
    check(a.ok && b.ok, std::string(label) + ": bitwise set-up");
    if (!a.ok || !b.ok)
        return;
    // two steps each (the second from the aspc warm start of the first), compared step by step
    Geometry g1 = mol.m_geometry;
    g1(0, 0) += 0.01;
    const Geometry* steps[2] = { &mol.m_geometry, &g1 };
    for (const Geometry* g : steps) {
        const Result ra = a.evaluate(*g, true), rb = b.evaluate(*g, true);
        const bool same_e = std::memcmp(&ra.e, &rb.e, sizeof(double)) == 0;
        const bool same_g = ra.grad.rows() == rb.grad.rows()
            && std::memcmp(ra.grad.data(), rb.grad.data(), sizeof(double) * ra.grad.size()) == 0;
        const bool same_q = ra.q.size() == rb.q.size() && std::memcmp(ra.q.data(), rb.q.data(), sizeof(double) * ra.q.size()) == 0;
        check(same_e && same_g && same_q, std::string(label) + " (" + std::to_string(threads) + " threads): two hybrid instances bit-identical (energy "
                                              + (same_e ? "yes" : "NO") + ", gradient " + (same_g ? "yes" : "NO") + ", charges " + (same_q ? "yes" : "NO") + ")");
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const Mol w20 = waterBox(20, 44.0, 20261105);
    const Mol w50 = waterBox(50, 44.0, 20261105);
    const Mol q1 = waterBox(18, 44.0, 20261106, { registryMol("hydronium_water2", 1) });
    const json ref_model = json{ { "frag_charge_model", "reference" } };

    testVsCpu("50 waters L44 (dense)", w50, json::object(), 1);
    testVsCpu("Q = +1 (H3O+ 2 H2O + 18 H2O)", q1, ref_model, 1);
    testVsCpu("970 waters L58 (pcg)", waterBox(970, 58.0, 20261105), json::object(), 4);

    testFiniteDifference("20 waters L44", w20, json::object(), 6);
    testFiniteDifference("Q = +1", q1, ref_model, 6);

    testInvariance("20 waters L44", w20, json::object(), true);
    testInvariance("Q = +1", q1, ref_model, false);   // the reference charge rule puts the net charge on the first fragment

    testBitwise("50 waters L44", w50, json::object(), 1);
    testBitwise("50 waters L44", w50, json::object(), 4);

    std::printf("GPU hybrid vs CPU over every comparison: charges %.2e e, energy (default tol) %.2e Eh, energy (tol 1e-14) %.2e Eh, "
                "gradient (tol 1e-14) %.2e Eh/A\n", g_worst_q, g_worst_e_def, g_worst_e_tight, g_worst_g_tight);
    std::printf(g_failures == 0 ? "gfnff_gpu_pbc_hybrid: all checks passed\n" : "gfnff_gpu_pbc_hybrid: %d failure(s)\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
#else
#include <cstdio>
int main()
{
    std::printf("gfnff_gpu_pbc_hybrid: built without CUDA, skipped\n");
    return 77;
}
#endif
