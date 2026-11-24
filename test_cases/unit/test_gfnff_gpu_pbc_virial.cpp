/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Atomic virial of the GFN-FF CUDA periodic hybrid (docs/GPU_PERIODIC_GAPS.md A1) - Claude Generated (Nov 2026).
 *
 * The hybrid (-gfnff.gpu_periodic true, coulomb_ewald) evaluates the local terms on the device and the EEQ and the
 * Ewald Coulomb pair term on the host; the virial Xi_ab = -dE/d(eps_ab) (fixed fractional coordinates, the N1
 * definition) comes from the device kernels (VIR instantiations and *_virial pair kernels) plus the host Ewald part.
 *
 * Checks (thresholds fixed before the first run):
 *  (a) GPU vs CPU periodic path, per term (bond incl. HB-CN, angle, dihedral, inversion, sTors, D4 pair, repulsion,
 *      HB, XB, BATM, Coulomb Ewald pair, background, D4 tail, V-LRC, CN chain rule as one sum), the model total and
 *      the pressure virial: <= 1e-7 Eh per tensor component.
 *  (b) FD gate of the hybrid surface: -dE/d(eps) by central differences of the GPU energy (fresh instance per strained
 *      geometry, charges re-solved) for an isotropic strain and the three axial strains eps_xx, eps_yy, eps_zz,
 *      eps 1e-5 / 2e-5, Richardson: |Xi_an - Xi_FD| <= 1e-7 Eh (trace resp. diagonal component).
 *  (c) reproducibility: a second GPU instance gives every tensor bit for bit; the gradient and energy with the virial
 *      requested are bitwise those without.
 *  (d) cell update on the GPU instance (GFNFF::updateCell, scale 1.001) == fresh GPU instance on the scaled box
 *      (energy 1e-10 Eh, gradient 1e-10 Eh/A, virial 1e-9 Eh per component).
 *  (e) NPT smoke test: SimpleMD 20 steps, Berendsen every 2 steps, on the CPU and with -gpu cuda: both runs finish;
 *      cell edge and pressure reported, cell edges equal to 1e-6 A.
 * Extra structures: `test_gfnff_gpu_pbc_virial --xyz FILE [FILE...]` runs (a) on XYZ files with a unitcell comment.
 * Exit 77 (ctest SKIP_RETURN_CODE) when built without CUDA.
 *
 * B2 (Nov 24, 2026, Claude Generated, ROCm parity, docs/ROCM_PARITY_B2.md): the same source builds the ROCm variant
 * (target test_gfnff_rocm_pbc_virial, ctest gfnff_rocm_pbc_virial, compiled with GFNFF_PBC_VIRIAL_ROCM): the HIP method
 * class GFNFFHipComputationalMethod and `-gpu rocm` in (e), same boxes, checks and thresholds. Exit 77 there when no
 * HIP device answers.
 */
#if defined(GFNFF_PBC_VIRIAL_ROCM) || defined(USE_CUDA)
#include "src/capabilities/simplemd.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#ifdef GFNFF_PBC_VIRIAL_ROCM
#include "src/core/energy_calculators/qm_methods/gfnff_hip_method.h"
using GFNFFGPUComputationalMethod = GFNFFHipComputationalMethod;   // B2: the HIP method class under the CUDA name
static const char* kGpuName = "rocm";
#else
#include "src/core/energy_calculators/qm_methods/gfnff_gpu_method.h"
static const char* kGpuName = "cuda";
#endif
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"
#include "src/core/fileiterator.h"
#include "src/core/gpu_fallback.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;

#include "generated/parameter_registry.h"
static struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int g_failures = 0;
static void check(bool ok, const std::string& what, double value, double limit)
{
    std::printf("  [%s] %-84s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), value, limit);
    if (!ok)
        ++g_failures;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
static const double kL = 44.0;  // Bohr

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

/// Rigid random packing (the construction of test_pbc_virial.cpp / test_pbc_eeq_ewald.cpp; same seeds, same boxes)
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

static Mol centredInBox(Mol m, double L_bohr)
{
    const double L_ang = L_bohr / A2B;
    const Eigen::RowVector3d c = m.m_geometry.colwise().mean();
    for (int i = 0; i < m.m_number_atoms; ++i)
        m.m_geometry.row(i) += Eigen::RowVector3d::Constant(0.5 * L_ang) - c;
    return withCell(m, L_ang, L_ang, L_ang);
}

/// Acetic acid dimer plus F-I...pyridine 6.5 A away (as in test_pbc_virial.cpp): XB, inversions, dihedrals, BATM
static Mol compositeMol()
{
    Mol a = registryMol("acetic_acid_dimer");
    const Mol b = registryMol("fi_pyridine");
    const Eigen::RowVector3d ca = a.m_geometry.colwise().mean(), cb = b.m_geometry.colwise().mean();
    const Eigen::RowVector3d off = ca - cb + Eigen::RowVector3d(0.0, 0.0, 6.5);
    const int na = a.m_number_atoms, nb = b.m_number_atoms;
    Geometry g(na + nb, 3);
    g.topRows(na) = a.m_geometry;
    for (int i = 0; i < nb; ++i)
        g.row(na + i) = b.m_geometry.row(i) + off;
    a.m_geometry = g;
    for (int i = 0; i < nb; ++i)
        a.m_atoms.push_back(b.m_atoms[i]);
    a.m_number_atoms = na + nb;
    return a;
}

/// Settings of test_pbc_virial (fixed Ewald parameters, absolute D4 radius) plus the S2 opt-in. `fd` = the FD
/// protocol (Ewald tolerance 1e-14, fixed alpha and r_cut); without it the defaults of a production run.
static json config(int threads, bool fd)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["gpu_periodic"] = true;
    c["hb_update_force_every"] = 1;
    c["hb_min_pair_energy_eh"] = 0.0;
    c["eeq_ewald_tol"] = 1e-14;
    c["frag_charge_model"] = "reference";
    if (fd) {
        c["ewald_tolerance"] = 1e-14;
        c["ewald_alpha"] = 0.32;
        c["ewald_r_cut"] = 19.0;
        c["dispersion_cutoff_bohr"] = 20.0;
    }
    return c;
}

enum Term { T_BOND, T_ANGLE, T_DIHED, T_INV, T_STORS, T_DISP, T_REP, T_HB, T_XB, T_BATM, T_COUL, T_COULBG, T_TAIL,
    T_LRC, T_HBX, T_XBX, T_CN, T_TOTAL, T_PRESSURE, NT };
static const char* kName[NT] = { "bond (+HB-CN)", "angle", "dihedral", "inversion", "stors", "D4 pair", "repulsion", "hbond",
    "xbond", "batm", "Coulomb Ewald pair", "Ewald background", "D4 tail", "V-LRC", "HB crossing (HB-X)", "XB crossing (HB-X)", "CN chain rule (sum)", "TOTAL (model)",
    "pressure virial" };

struct Eval {
    bool ok = false;
    double e = 0.0;
    Matrix grad;
    Eigen::Matrix3d xi[NT];
};

static void fill(const FFVirialComponents& v, Eval& r)
{
    r.xi[T_BOND] = v.bond;
    r.xi[T_ANGLE] = v.angle;
    r.xi[T_DIHED] = v.dihedral;
    r.xi[T_INV] = v.inversion;
    r.xi[T_STORS] = v.stors;
    r.xi[T_DISP] = v.dispersion;
    r.xi[T_REP] = v.repulsion;
    r.xi[T_HB] = v.hbond;
    r.xi[T_XB] = v.xbond;
    r.xi[T_BATM] = v.batm;
    r.xi[T_COUL] = v.coulomb;
    r.xi[T_COULBG] = v.coulomb_bg;
    r.xi[T_TAIL] = v.disp_tail;
    r.xi[T_LRC] = v.tail_lrc;
    r.xi[T_HBX] = v.hb_crossing;   // HB-X (Nov 2026)
    r.xi[T_XBX] = v.xb_crossing;
    r.xi[T_CN] = v.cn_bond + v.cn_dispersion + v.cn_coulomb;
    r.xi[T_TOTAL] = v.total();
    r.xi[T_PRESSURE] = v.pressureTotal();
}

/// One fresh instance (GPU hybrid or CPU) through the ComputationalMethod wrappers; gradient in Eh/A.
static Eval evaluate(bool gpu, const Mol& m, const json& cfg, bool virial)
{
    Eval r;
    for (auto& x : r.xi)
        x.setZero();
    try {
        std::unique_ptr<ComputationalMethod> meth;
        if (gpu) {
            auto g = std::make_unique<GFNFFGPUComputationalMethod>("gfnff", cfg);
            if (!g->setMolecule(m) || g->hasError() || !g->periodicHybridActive()) {
                std::printf("  GPU set-up failed: %s\n", g->getErrorMessage().c_str());
                return r;
            }
            meth = std::move(g);
        } else {
            auto c = std::make_unique<GFNFFComputationalMethod>("gfnff", cfg);
            if (!c->setMolecule(m) || c->hasError())
                return r;
            meth = std::move(c);
        }
        GFNFF* ff = meth->gfnffInstance();
        if (!ff)
            return r;
        if (virial)
            ff->setWantVirial(true);
        r.e = meth->calculateEnergy(true);
        r.grad = meth->getGradient();
        if (virial) {
            if (!ff->virialValid() || !ff->virialComponents()) {
                std::printf("  virial not valid (%s)\n", gpu ? "GPU" : "CPU");
                return r;
            }
            fill(*ff->virialComponents(), r);
        }
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("  evaluate threw: %s\n", e.what());
    }
    return r;
}

static double g_worst_term = 0.0;

// (a) GPU vs CPU per term
static void testVsCpu(const char* label, const Mol& m, int threads, bool fd_protocol)
{
    std::printf("(a) GPU vs CPU: %s (%d atoms, L %.2f Bohr, %d threads)\n", label, m.m_number_atoms, m.m_unit_cell(0, 0) * A2B, threads);
    const json cfg = config(threads, fd_protocol);
    const auto t0 = std::chrono::steady_clock::now();
    const Eval g = evaluate(true, m, cfg, true);
    const auto t1 = std::chrono::steady_clock::now();
    const Eval c = evaluate(false, m, cfg, true);
    const auto t2 = std::chrono::steady_clock::now();
    check(g.ok && c.ok, std::string(label) + ": both virials computed", g.ok && c.ok ? 1 : 0, 1);
    if (!g.ok || !c.ok)
        return;
    std::printf("  set-up + one gradient call: GPU %.0f ms, CPU %.0f ms; dE %.2e Eh, max dG %.2e Eh/A\n",
        std::chrono::duration<double, std::milli>(t1 - t0).count(), std::chrono::duration<double, std::milli>(t2 - t1).count(),
        g.e - c.e, (g.grad - c.grad).cwiseAbs().maxCoeff());
    for (int t = 0; t < NT; ++t) {
        const double d = (g.xi[t] - c.xi[t]).cwiseAbs().maxCoeff();
        if (t < T_TOTAL)
            g_worst_term = std::max(g_worst_term, d);
        std::printf("    %-24s tr GPU %+.10e  CPU %+.10e  max|d| %.2e\n", kName[t], g.xi[t].trace(), c.xi[t].trace(), d);
        check(d <= 1e-7, std::string(label) + " " + kName[t] + ": max |Xi_GPU - Xi_CPU| per component (Eh)", d, 1e-7);
        if (d > 1e-7)
            for (int a = 0; a < 3; ++a)
                std::printf("      row %d GPU %+.6e %+.6e %+.6e   CPU %+.6e %+.6e %+.6e\n", a, g.xi[t](a, 0), g.xi[t](a, 1), g.xi[t](a, 2),
                    c.xi[t](a, 0), c.xi[t](a, 1), c.xi[t](a, 2));
    }
    const double V = m.m_unit_cell.determinant() * std::pow(A2B, 3);
    std::printf("  virial pressure GPU %.4f bar, CPU %.4f bar\n", CurcumaUnit::Pressure::au_to_bar(g.xi[T_PRESSURE].trace() / (3.0 * V)),
        CurcumaUnit::Pressure::au_to_bar(c.xi[T_PRESSURE].trace() / (3.0 * V)));
}

static Mol strained(const Mol& m, int axis, double eps)
{
    // axis -1: isotropic; else only that Cartesian axis (x_a and L_a scaled)
    Mol s = m;
    for (int k = 0; k < 3; ++k)
        if (axis < 0 || axis == k) {
            s.m_geometry.col(k) *= (1.0 + eps);
            s.m_unit_cell(k, k) *= (1.0 + eps);
        }
    return s;
}

// (b) FD gate of the hybrid energy under cell scaling
static void testFD(const char* label, const Mol& m)
{
    std::printf("(b) FD gate of the GPU hybrid: %s (%d atoms)\n", label, m.m_number_atoms);
    const json cfg = config(1, true);
    const Eval an = evaluate(true, m, cfg, true);
    check(an.ok, std::string(label) + ": analytic GPU virial", an.ok ? 1 : 0, 1);
    if (!an.ok)
        return;
    const double eps[4] = { -2e-5, -1e-5, 1e-5, 2e-5 };
    for (int axis = -1; axis < 3; ++axis) {
        double e[4];
        bool ok = true;
        for (int s = 0; s < 4; ++s) {
            const Eval es = evaluate(true, strained(m, axis, eps[s]), cfg, false);
            ok = ok && es.ok;
            e[s] = es.e;
        }
        check(ok, std::string(label) + ": strained GPU instances", ok ? 1 : 0, 1);
        if (!ok)
            return;
        const double fd1 = -(e[2] - e[1]) / 2e-5, fd2 = -(e[3] - e[0]) / 4e-5;
        const double rich = (4.0 * fd1 - fd2) / 3.0;
        const double xan = axis < 0 ? an.xi[T_TOTAL].trace() : an.xi[T_TOTAL](axis, axis);
        const char* nm[4] = { "isotropic (trace)", "eps_xx", "eps_yy", "eps_zz" };
        std::printf("    %-18s Xi_an %+.10e  FD(1e-5) %+.10e  Richardson %+.10e  |d| %.2e / %.2e\n", nm[axis + 1], xan, fd1, rich,
            std::abs(xan - fd1), std::abs(xan - rich));
        check(std::abs(xan - rich) <= 1e-7, std::string(label) + " FD " + nm[axis + 1] + ": |Xi_an - Xi_FD(Richardson)| (Eh)",
            std::abs(xan - rich), 1e-7);
    }
}

static bool bitEqual(const Eigen::Matrix3d& a, const Eigen::Matrix3d& b) { return std::memcmp(a.data(), b.data(), 9 * sizeof(double)) == 0; }

// (c) reproducibility and no effect of the request on energy and gradient
static void testBitwise(const char* label, const Mol& m)
{
    std::printf("(c) reproducibility: %s\n", label);
    const json cfg = config(4, false);
    const Eval a = evaluate(true, m, cfg, true), b = evaluate(true, m, cfg, true), p = evaluate(true, m, cfg, false);
    check(a.ok && b.ok && p.ok, std::string(label) + ": three GPU instances", 1, 1);
    if (!a.ok || !b.ok || !p.ok)
        return;
    bool same = a.e == b.e && a.grad.size() == b.grad.size()
        && std::memcmp(a.grad.data(), b.grad.data(), a.grad.size() * sizeof(double)) == 0;
    for (int t = 0; t < NT; ++t)
        same = same && bitEqual(a.xi[t], b.xi[t]);
    check(same, std::string(label) + ": second GPU instance: energy, gradient, every tensor bit for bit", same ? 0 : 1, 0);
    const bool plain = a.e == p.e && a.grad.size() == p.grad.size()
        && std::memcmp(a.grad.data(), p.grad.data(), a.grad.size() * sizeof(double)) == 0;
    check(plain, std::string(label) + ": energy and gradient with the virial == without (bitwise)",
        plain ? 0.0 : std::max(std::abs(a.e - p.e), (a.grad - p.grad).cwiseAbs().maxCoeff()), 0);
}

// (d) the GPU instance follows GFNFF::updateCell
static void testCellUpdate(const char* label, const Mol& m)
{
    std::printf("(d) cell update on the GPU instance: %s\n", label);
    const json cfg = config(1, true);
    const double mu = 1.001;
    Mol ms = m;
    ms.m_geometry *= mu;
    ms.m_unit_cell *= mu;
    try {
        GFNFFGPUComputationalMethod g("gfnff", cfg);
        if (!g.setMolecule(m) || !g.periodicHybridActive()) {
            check(false, "set-up", 0, 1);
            return;
        }
        g.gfnffInstance()->setWantVirial(true);
        g.calculateEnergy(true);
        g.gfnffInstance()->updateCell(ms.m_unit_cell);
        g.updateGeometry(ms.m_geometry);
        const double e = g.calculateEnergy(true);
        const Matrix grad = g.getGradient();
        Eval upd;
        upd.ok = g.gfnffInstance()->virialValid();
        if (upd.ok)
            fill(*g.gfnffInstance()->virialComponents(), upd);
        const Eval fresh = evaluate(true, ms, cfg, true);
        check(upd.ok && fresh.ok, std::string(label) + ": virial after the update and of the fresh instance", 1, 1);
        if (!upd.ok || !fresh.ok)
            return;
        const double dv = (upd.xi[T_PRESSURE] - fresh.xi[T_PRESSURE]).cwiseAbs().maxCoeff();
        std::printf("    E %.12f vs fresh %.12f; tr Xi %.10e vs %.10e\n", e, fresh.e, upd.xi[T_PRESSURE].trace(), fresh.xi[T_PRESSURE].trace());
        check(std::abs(e - fresh.e) <= 1e-10, std::string(label) + ": energy updated == fresh (Eh)", std::abs(e - fresh.e), 1e-10);
        check((grad - fresh.grad).cwiseAbs().maxCoeff() <= 1e-10, std::string(label) + ": gradient updated == fresh (Eh/A)",
            (grad - fresh.grad).cwiseAbs().maxCoeff(), 1e-10);
        check(dv <= 1e-9, std::string(label) + ": virial updated == fresh (Eh, per component)", dv, 1e-9);
    } catch (const std::exception& ex) {
        std::printf("  threw: %s\n", ex.what());
        check(false, "cell update ran", 0, 1);
    }
}

// (e) NPT smoke test through SimpleMD
static json mdController(bool gpu)
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = "gfnff";
    md["no_restart"] = true;
    md["thermostat"] = "csvr";
    md["coupling"] = 100.0;
    md["temperature"] = 300.0;
    md["seed"] = 7;
    md["dt"] = 0.5;
    md["max_time"] = 10.0;   // 20 steps
    md["no_center"] = true;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    md["barostat"] = "berendsen";
    md["pressure"] = 1.01325;
    md["tau_p"] = 100.0;
    md["compressibility"] = 4.5e-5;
    md["pcoupl_interval"] = 2;
    controller["gfnff"]["frag_charge_model"] = "reference";
    controller["gfnff"]["cache_topology"] = false;
    controller["gfnff"]["gpu_periodic"] = true;
    if (gpu) {
        md["gpu"] = kGpuName;
        controller["gpu"] = kGpuName;
    }
    return controller;
}

static void testNPT(const Mol& m)
{
    std::printf("(e) NPT smoke test, Berendsen, 20 steps: CPU vs -gpu %s (%d atoms)\n", kGpuName, m.m_number_atoms);
    double L[2] = { 0, 0 }, P[2] = { 0, 0 };
    long ev[2] = { 0, 0 };
    bool ok[2] = { false, false };
    for (int k = 0; k < 2; ++k) {
        try {
            int fb0 = 0;
            for (const auto& c : curcuma::gpuFallbackCounts()) fb0 += c.second;
            json controller = mdController(k == 1);
            SimpleMD sim(controller, true);
            sim.overrideBasename(k == 1 ? "gpu_pbc_virial_npt_gpu" : "gpu_pbc_virial_npt_cpu");
            sim.setMolecule(curcuma::Molecule(m));
            if (!sim.Initialise())
                continue;
            sim.prepareRun();
            int steps = 0;
            while (sim.step())
                ++steps;
            L[k] = sim.unitCell()(0, 0);
            P[k] = sim.pressureBar();
            ev[k] = sim.barostatEvents();
            ok[k] = steps >= 19 && (sim.stopReason() == SimpleMD::StopReason::Running || sim.stopReason() == SimpleMD::StopReason::MaxTime);
            int fb1 = 0;
            for (const auto& c : curcuma::gpuFallbackCounts()) fb1 += c.second;
            if (k == 1)
                check(fb1 == fb0, "(e) the GPU run counted no GPU fallback (it ran on the device)", fb1 - fb0, 0);
            std::printf("    %s: %d steps, stop %d, events %ld, L %.9f A, P %.4f bar\n", k ? "GPU" : "CPU", steps,
                static_cast<int>(sim.stopReason()), ev[k], L[k], P[k]);
        } catch (const std::exception& e) {
            std::printf("    %s threw: %s\n", k ? "GPU" : "CPU", e.what());
        }
    }
    check(ok[0] && ok[1], "(e) both NPT runs finished", ok[0] && ok[1] ? 1 : 0, 1);
    check(ev[0] == ev[1] && ev[1] > 0, "(e) same number of barostat couplings", static_cast<double>(ev[1]), static_cast<double>(ev[0]));
    check(std::abs(L[0] - L[1]) <= 1e-6, "(e) cell edge CPU vs GPU after 20 steps (A)", std::abs(L[0] - L[1]), 1e-6);
    std::printf("    pressure CPU %.4f GPU %.4f bar (difference %.2e bar)\n", P[0], P[1], P[1] - P[0]);
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc > 2 && std::strcmp(argv[1], "--xyz") == 0) {
        for (int a = 2; a < argc; ++a) {
            FileIterator it(argv[a], true);
            curcuma::Molecule mol = it.Next();
            Mol m = mol.getMolInfo();
            if (!m.m_has_pbc) {
                std::printf("%s: no unit cell in the comment line\n", argv[a]);
                ++g_failures;
                continue;
            }
            testVsCpu(argv[a], m, 4, false);
        }
        std::printf("worst per-term component difference: %.3e Eh\n", g_worst_term);
        std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
        return g_failures == 0 ? 0 : 1;
    }
#ifdef GFNFF_PBC_VIRIAL_ROCM
    {   // B2: skip (77) when no HIP device answers
        GFNFFGPUComputationalMethod p("gfnff", config(1, false));
        bool ok = false;
        try {
            ok = p.setMolecule(centredInBox(registryMol("H2O"), kL)) && !p.hasError() && p.periodicHybridActive();
        } catch (const std::exception&) {
            ok = false;
        }
        if (!ok) {
            std::printf("gfnff_rocm_pbc_virial: the ROCm periodic hybrid could not be set up (no HIP device?), skipped\n");
            return 77;
        }
    }
#endif
    const char* only = std::getenv("GPU_VIRIAL_ONLY");   // diagnostic: a/b/c/d/e
    auto run = [&](char c) { return !only || std::strchr(only, c); };
    const Mol box20 = waterBox(20, kL, 20261105);
    const Mol box60 = waterBox(60, kL, 20261108);
    const Mol charged = waterBox(18, kL, 20261106, { registryMol("hydronium_water2", 1) });
    const Mol caffeine = centredInBox(registryMol("caffeine"), kL);
    const Mol composite = centredInBox(compositeMol(), kL);
    if (run('a')) testVsCpu("20 H2O", box20, 1, true);
    if (run('a')) testVsCpu("60 H2O", box60, 4, false);
    if (run('a')) testVsCpu("H3O+(H2O)2 + 18 H2O, Q = +1", charged, 1, true);
    if (run('a')) testVsCpu("caffeine", caffeine, 1, true);
    if (run('a')) testVsCpu("acetic dimer + F-I...pyridine", composite, 1, true);
    if (run('b')) testFD("20 H2O", box20);
    if (run('b')) testFD("acetic dimer + F-I...pyridine", composite);
    if (run('c')) testBitwise("60 H2O", box60);
    if (run('d')) testCellUpdate("20 H2O", box20);
    if (run('e')) testNPT(box20);
    std::printf("worst per-term component difference: %.3e Eh\n", g_worst_term);
    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
    return g_failures == 0 ? 0 : 1;
}
#else
#include <cstdio>
int main()
{
    std::printf("built without CUDA: skipped\n");
    return 77;
}
#endif
