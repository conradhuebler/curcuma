/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * S1 of the periodic GPU port (Claude Generated, Nov 2026): the local GFN-FF terms with an orthorhombic cell on the
 * CUDA device against the CPU, term by term, without Coulomb/EEQ (coulomb false, coulomb_ewald false, gpu_periodic true).
 *
 * Structures (registry, 44 Bohr cube): water dimer, caffeine, hydronium-water2 (+1, frag_charge_model reference),
 * acetic acid dimer, acetic acid dimer + F-I...pyridine composite, 40 and 60 random waters (intact molecules).
 * States: cell centre, shifted by (L,-L,2L), fragments by different lattice vectors, centred on the corner, torn per
 * atom (the topology is built by the host from the reassembled geometry, A4c), 2x1x1 supercell.
 *
 * Checks (thresholds fixed beforehand):
 *  (a) GPU vs CPU per state: every term energy |dE| <= 1e-11 Eh, gradient (all terms without Coulomb, incl. the CN
 *      chain rule) max |dG| <= 1e-12 Eh/A; the same with dispersion / repulsion switched off (centre, torn), which
 *      isolates the D4 (+ tail) and repulsion gradients as differences (the hbond switch is asymmetric, see below).
 *  (b) GPU invariance: every state equals the centred one, terms <= 1e-11 Eh, gradient <= 1e-12 Eh/A, same HB count;
 *      BATM has no minimum image on either side (A3b), so the batm term and, where BATM is present, the gradient
 *      are left out of the invariance check (GPU == CPU still holds for them).
 *  (c) 2x1x1 supercell on the GPU (intact molecules, neutral systems): every term twice the cell (<= 1e-11 Eh per
 *      term, relative 1e-12 for large terms), per-atom gradient equal.
 *  (d) central finite difference (Angstrom, h = 2e-4 and 1e-4) of the GPU energy for the atoms nearest a face of the
 *      torn caffeine, composite and 40-water box (lists fixed): error <= 5e-8 Eh/A at h = 1e-4, quadratic in h,
 *      Richardson-extrapolated difference <= 2e-9 Eh/A.
 *  (e) Verlet lists on the device: wrapping one molecule by +L rebuilds no D4 / repulsion list and keeps the energy.
 *  (f) two device instances on the same geometry: energy and gradient bit for bit identical.
 * The 60-water box also runs with nb_cell_list_min_atoms 0 (device repulsion list build, cell-list HB branch).
 * Exit 77 (ctest SKIP_RETURN_CODE) when built without CUDA.
 *
 * A2 (Nov 2026, Claude Generated): the same source builds the ROCm variant (target test_gfnff_rocm_pbc_terms, ctest
 * gfnff_rocm_pbc_terms, compiled with GFNFF_PBC_TERMS_ROCM): the HIP method class GFNFFHipComputationalMethod in place
 * of the CUDA one, same structures, states and thresholds. Exit 77 there when no HIP device answers.
 */
#if defined(GFNFF_PBC_TERMS_ROCM) || defined(USE_CUDA)
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#ifdef GFNFF_PBC_TERMS_ROCM
#include "src/core/energy_calculators/qm_methods/gfnff_hip_method.h"
using GpuGfnffMethod = GFNFFHipComputationalMethod;
static const char* kTestName = "gfnff_rocm_pbc_terms";
#else
#include "src/core/energy_calculators/qm_methods/gfnff_gpu_method.h"
using GpuGfnffMethod = GFNFFGPUComputationalMethod;
static const char* kTestName = "gfnff_gpu_pbc_terms";
#endif
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
static const double kCellBohr = 44.0;
static const double kTolE = 1e-11;    // Eh per term
static const double kTolG = 1e-12;    // Eh/A

static Mol registryMol(const char* name, int charge = 0)
{
    Mol m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
    m.m_charge = charge;
    return m;
}

static Mol withCell(Mol mol, const Eigen::Vector3d& L_ang)
{
    mol.m_has_pbc = L_ang.minCoeff() > 0.0;
    if (mol.m_has_pbc)
        mol.m_unit_cell = L_ang.asDiagonal();
    return mol;
}

static Geometry centred(const Geometry& g0, const Eigen::Vector3d& at)
{
    Geometry g = g0;
    Eigen::RowVector3d c = g.colwise().mean();
    for (int i = 0; i < g.rows(); ++i)
        g.row(i) += at.transpose() - c;
    return g;
}

/// Acetic acid dimer plus F-I...pyridine 6.5 A away (as in pbc_hbond).
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

/// Random water box (intact molecules, minimum-image O-O >= 2.9 A), as in pbc_hbond / pbc_bonded_repulsion.
static Mol waterBox(int nwat, double L_ang, unsigned seed)
{
    const Mol d = registryMol("H2O_dimer");
    Eigen::MatrixXd w(3, 3);
    for (int i = 0; i < 3; ++i)
        w.row(i) = d.m_geometry.row(i) - d.m_geometry.row(0);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<Eigen::Vector3d> O;
    Mol m = d;
    m.m_atoms.clear();
    m.m_geometry = Geometry(3 * nwat, 3);
    int placed = 0, tries = 0;
    while (placed < nwat && tries < 200000) {
        ++tries;
        Eigen::Vector3d p(u(rng) * L_ang, u(rng) * L_ang, u(rng) * L_ang);
        bool ok = true;
        for (const auto& q : O) {
            Eigen::Vector3d dv = p - q;
            for (int k = 0; k < 3; ++k)
                dv[k] -= L_ang * std::rint(dv[k] / L_ang);
            if (dv.norm() < 2.9) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
        qr.normalize();
        const Eigen::Matrix3d R = qr.toRotationMatrix();
        for (int i = 0; i < 3; ++i) {
            m.m_geometry.row(3 * placed + i) = (p + R * w.row(i).transpose()).transpose();
            m.m_atoms.push_back(d.m_atoms[i]);
        }
        O.push_back(p);
        ++placed;
    }
    m.m_number_atoms = 3 * placed;
    m.m_geometry.conservativeResize(3 * placed, 3);
    check(placed == nwat, "water box placed all molecules");
    return m;
}

static std::vector<std::vector<int>> fragmentsOf(const Mol& m)
{
    curcuma::Molecule mol(m);
    return mol.GetFragments();
}

static json config(const json& extra)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = 1;
    c["cache_topology"] = false;
    c["coulomb"] = false;          // the Coulomb term and the EEQ enter no checked number (S3/S4)
    c["coulomb_ewald"] = false;    // the only setting the S1 guard releases for the device
    c["gpu_periodic"] = true;
    c["hb_update_force_every"] = 1;  // HB/XB lists from the geometry of every call
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    // the GPU wrapper reads the term switches from the "gfnff" scope only (the host GFNFF from the flat keys)
    json scope = json::object();
    for (const char* key : { "coulomb", "hbond", "dispersion", "repulsion" })
        if (c.contains(key)) scope[key] = c[key];
    c["gfnff"] = scope;
    return c;
}

// --------------------------------------------------------------------------------------------------------------
struct Terms {
    double bond = 0, angle = 0, dihedral = 0, inversion = 0, stors = 0, brep = 0, nbrep = 0, disp = 0, tail = 0,
           hb = 0, xb = 0, batm = 0, total = 0;
    int n_hb = 0;
    Matrix grad;   // Eh/A
    double operator[](int k) const
    {
        const double v[] = { bond, angle, dihedral, inversion, stors, brep, nbrep, disp, tail, hb, xb, batm };
        return v[k];
    }
};
static const char* kTermNames[] = { "bond", "angle", "dihedral", "inversion", "stors", "bonded_rep", "nonbonded_rep",
                                    "dispersion", "disp_tail", "hbond", "xbond", "batm" };
static constexpr int kNTerms = 12;

static Terms fromComponents(const FFEnergyComponents& c)
{
    Terms t;
    t.bond = c.bond; t.angle = c.angle; t.dihedral = c.dihedral; t.inversion = c.inversion; t.stors = c.stors;
    t.brep = c.bonded_rep; t.nbrep = c.nonbonded_rep; t.disp = c.dispersion; t.tail = c.disp_tail;
    t.hb = c.hbond; t.xb = c.xbond; t.batm = c.batm;
    return t;
}

/// One CPU or GPU instance (the ComputationalMethod wrappers, so the gradient contract Eh/A is the same).
struct Instance {
    std::unique_ptr<GFNFFComputationalMethod> cpu;
    std::unique_ptr<GpuGfnffMethod> gpu;
    bool ok = false;
    Instance(bool on_gpu, const Mol& mol, const json& c)
    {
        try {
            if (on_gpu) {
                gpu = std::make_unique<GpuGfnffMethod>("gfnff", c);
                ok = gpu->setMolecule(mol) && !gpu->hasError() && !gpu->requestsCpuFallback();
                if (!ok) std::printf("GPU init: %s%s\n", gpu->getErrorMessage().c_str(),
                                     gpu->requestsCpuFallback() ? " (CPU fallback requested)" : "");
            } else {
                cpu = std::make_unique<GFNFFComputationalMethod>("gfnff", c);
                ok = cpu->setMolecule(mol) && !cpu->hasError();
            }
        } catch (const std::exception& e) {
            std::printf("init threw: %s\n", e.what());
            ok = false;
        }
    }
    GFNFF* ff() const { return gpu ? gpu->getGFNFF() : cpu->getGFNFF(); }
    Terms evaluate(const Geometry& g_ang, bool gradient)
    {
        Terms t;
        if (gpu) {
            gpu->updateGeometry(g_ang);
            const double e = gpu->calculateEnergy(gradient);
            t = fromComponents(gpu->getGPUWorkspace()->energyComponents());
            t.total = e;
            if (gradient) t.grad = gpu->getGradient();
            t.n_hb = gpu->getGPUWorkspace()->hbondCount();
        } else {
            cpu->updateGeometry(g_ang);
            const double e = cpu->calculateEnergy(gradient);
            t = fromComponents(cpu->getGFNFF()->getWorkspace()->energyComponents());
            t.total = e;
            if (gradient) t.grad = cpu->getGradient();
            t.n_hb = static_cast<int>(cpu->getGFNFF()->getLastHBonds().size());
        }
        return t;
    }
};

static double maxAbs(const Matrix& a, const Matrix& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols())
        return 1e300;
    return a.rows() ? (a - b).cwiseAbs().maxCoeff() : 0.0;
}

/// max per-term |a - b| (and its name) and the gradient difference
struct Diff {
    double e = 0, g = 0;
    int worst = 0;
};
static Diff diffTerms(const Terms& a, const Terms& b, bool with_grad)
{
    Diff d;
    for (int k = 0; k < kNTerms; ++k) {
        const double x = std::abs(a[k] - b[k]);
        if (x > d.e) { d.e = x; d.worst = k; }
    }
    if (with_grad) d.g = maxAbs(a.grad, b.grad);
    return d;
}

struct Variants {
    Geometry mid, shift, frag, corner, torn;
};
static Variants makeVariants(const Mol& base)
{
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    Variants v;
    v.mid = centred(base.m_geometry, 0.5 * Lv);
    v.corner = centred(base.m_geometry, Eigen::Vector3d::Zero());
    v.shift = v.mid;
    for (int i = 0; i < n; ++i) {
        v.shift(i, 0) += L;
        v.shift(i, 1) -= L;
        v.shift(i, 2) += 2.0 * L;
    }
    v.frag = v.mid;
    const auto frags = fragmentsOf(base);
    const int shifts[4][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { -1, 1, 0 }, { 0, -1, 2 } };
    for (size_t f = 0; f < frags.size(); ++f)
        for (int a : frags[f])
            for (int k = 0; k < 3; ++k)
                v.frag(a, k) += shifts[f % 4][k] * L;
    v.torn = v.corner;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            v.torn(i, k) -= L * std::floor(v.torn(i, k) / L);
    return v;
}

static double g_worst_e = 0.0, g_worst_g = 0.0;   // over every GPU-vs-CPU comparison
static double g_worst_term_e[kNTerms] = { 0 };

static void compareOne(const std::string& label, const Terms& gpu, const Terms& cpu, double tol_e, double tol_g)
{
    const Diff d = diffTerms(gpu, cpu, true);
    for (int k = 0; k < kNTerms; ++k)
        g_worst_term_e[k] = std::max(g_worst_term_e[k], std::abs(gpu[k] - cpu[k]));
    g_worst_e = std::max(g_worst_e, d.e);
    g_worst_g = std::max(g_worst_g, d.g);
    check(d.e <= tol_e, label + ": GPU vs CPU terms (worst " + kTermNames[d.worst] + " " + sci(d.e) + " Eh)");
    check(d.g <= tol_g, label + ": GPU vs CPU gradient (" + sci(d.g) + " Eh/A)");
    check(gpu.n_hb == cpu.n_hb, label + ": HB count GPU " + std::to_string(gpu.n_hb) + " vs CPU " + std::to_string(cpu.n_hb));
    std::printf("%-46s terms %.2e (%s)  grad %.2e Eh/A  HB %d/%d\n", label.c_str(), d.e, kTermNames[d.worst], d.g,
                gpu.n_hb, cpu.n_hb);
}

// (a) + (b): GPU vs CPU per state, GPU invariance between states
static void testStructure(const char* label, const Mol& base, const json& extra)
{
    const double L = kCellBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    const Variants v = makeVariants(base);
    Mol m_mid = withCell(base, Lv);
    m_mid.m_geometry = v.mid;

    Instance gpu(true, m_mid, config(extra));
    Instance cpu(false, m_mid, config(extra));
    check(gpu.ok && cpu.ok, std::string(label) + ": set-up in the cell centre (GPU, CPU)");
    if (!gpu.ok || !cpu.ok)
        return;
    const Terms g_ref = gpu.evaluate(v.mid, true);
    const Terms c_ref = cpu.evaluate(v.mid, true);
    compareOne(std::string(label) + " centre", g_ref, c_ref, kTolE, kTolG);
    std::printf("%-46s E_tail GPU %.12e CPU %.12e, E_HB %.10e, E_disp %.10e, total %.12f\n", label, g_ref.tail,
                c_ref.tail, g_ref.hb, g_ref.disp, g_ref.total);
    check(std::abs(g_ref.tail) > 1e-9 || std::abs(c_ref.tail) < 1e-300,
          std::string(label) + ": dispersion tail present on the device (" + sci(g_ref.tail) + ")");

    struct V { const char* what; const Geometry* g; };
    const V variants[] = { { "shifted by (L,-L,2L)", &v.shift }, { "fragments by lattice vectors", &v.frag },
                           { "centred on the corner", &v.corner }, { "torn (per-atom wrap)", &v.torn } };
    // BATM (bonded ATM, 1,4 pairs) has no minimum image on the CPU nor on the GPU (A3b: parity, wrong for torn
    // molecules on both sides): GPU == CPU holds for every state (above), but a torn geometry is not invariant in
    // BATM. The invariance check therefore skips the batm term and, where BATM is non-zero, the gradient.
    const bool has_batm = std::abs(g_ref.batm) > 0.0;
    auto invariance = [&](const Terms& g, const std::string& what) {
        Diff inv = diffTerms(g, g_ref, true);
        double e = 0.0;
        int worst = 0;
        for (int k = 0; k < kNTerms; ++k) {
            if (std::string(kTermNames[k]) == "batm") continue;
            const double x = std::abs(g[k] - g_ref[k]);
            if (x > e) { e = x; worst = k; }
        }
        const bool grad_ok = has_batm || inv.g <= kTolG;
        check(e <= kTolE && grad_ok && g.n_hb == g_ref.n_hb,
              std::string(label) + " " + what + ": GPU == GPU centred (terms " + sci(e) + " " + kTermNames[worst] + ", grad "
                  + (has_batm ? std::string("skipped, BATM ") + sci(inv.g) : sci(inv.g)) + ", HB " + std::to_string(g.n_hb) + "/"
                  + std::to_string(g_ref.n_hb) + ")");
    };
    for (const auto& vv : variants) {
        const Terms g = gpu.evaluate(*vv.g, true);
        const Terms c = cpu.evaluate(*vv.g, true);
        compareOne(std::string(label) + " " + vv.what, g, c, kTolE, kTolG);
        invariance(g, vv.what);
    }

    // torn input initialised directly (A4c on the host): terms as the centred instance
    {
        Mol m_torn = withCell(base, Lv);
        m_torn.m_geometry = v.torn;
        Instance gt(true, m_torn, config(extra));
        check(gt.ok, std::string(label) + ": torn input set up directly on the GPU");
        if (gt.ok)
            invariance(gt.evaluate(v.torn, true), "torn, initialised torn");
    }

    // (a) with single terms off: isolates the D4 (+ tail) and repulsion gradients as differences. Not the hbond
    // switch: with `hbond false` the CPU drops the HB detection and with it the HB-CN modulation of the X-H bonds,
    // the GPU keeps that modulation (pre-existing asymmetry of a diagnostic switch, TODO.md); the HB gradient is
    // covered by the full comparison of every state above.
    const char* offs[] = { "dispersion", "repulsion" };
    for (const char* off : offs) {
        json e2 = extra;
        e2[off] = false;
        Instance go(true, m_mid, config(e2));
        Instance co(false, m_mid, config(e2));
        check(go.ok && co.ok, std::string(label) + ": set-up with " + off + " off");
        if (!go.ok || !co.ok)
            continue;
        for (const V& vv : { V{ "centre", &v.mid }, V{ "torn", &v.torn } }) {
            const Terms g = go.evaluate(*vv.g, true);
            const Terms c = co.evaluate(*vv.g, true);
            compareOne(std::string(label) + " " + vv.what + " (" + off + " off)", g, c, kTolE, kTolG);
        }
        // the isolated term gradient (full minus off) GPU vs CPU
        const Terms g = go.evaluate(v.torn, true), c = co.evaluate(v.torn, true);
        const Terms gf = gpu.evaluate(v.torn, true), cf = cpu.evaluate(v.torn, true);
        const Matrix dg = gf.grad - g.grad, dc = cf.grad - c.grad;
        const double dd = maxAbs(dg, dc);
        g_worst_g = std::max(g_worst_g, dd);
        check(dd <= 2.0 * kTolG, std::string(label) + ": isolated " + off + " gradient GPU vs CPU (" + sci(dd) + " Eh/A, max |g| "
                                     + sci(dg.cwiseAbs().maxCoeff()) + ")");
        std::printf("%-46s isolated %-10s gradient GPU vs CPU %.2e (max |g| %.2e Eh/A)\n", label, off, dd, dg.cwiseAbs().maxCoeff());
    }
}

// (c) 2x1x1 supercell on the GPU
static void testSupercell(const char* label, const Mol& base, const json& extra)
{
    // intact molecules (cell centre): BATM has no minimum image, so a torn copy would differ from its image
    const double L = kCellBohr / A2B;
    const Variants v = makeVariants(base);
    Mol m1 = withCell(base, Eigen::Vector3d(L, L, L));
    m1.m_geometry = v.mid;
    Mol m2 = m1;
    m2.m_has_pbc = true;
    m2.m_unit_cell = Eigen::Vector3d(2.0 * L, L, L).asDiagonal();
    const int n = m1.m_number_atoms;
    m2.m_geometry = Geometry(2 * n, 3);
    m2.m_geometry.topRows(n) = m1.m_geometry;
    m2.m_geometry.bottomRows(n) = m1.m_geometry;
    for (int i = 0; i < n; ++i)
        m2.m_geometry(n + i, 0) += L;
    for (int i = 0; i < n; ++i)
        m2.m_atoms.push_back(m1.m_atoms[i]);
    m2.m_number_atoms = 2 * n;
    m2.m_charge = 2 * m1.m_charge;
    Instance g1(true, m1, config(extra)), g2(true, m2, config(extra));
    check(g1.ok && g2.ok, std::string(label) + ": supercell set-up");
    if (!g1.ok || !g2.ok)
        return;
    const Terms t1 = g1.evaluate(m1.m_geometry, true), t2 = g2.evaluate(m2.m_geometry, true);
    double worst = 0.0, worst_d = 0.0;
    int wk = 0;
    for (int k = 0; k < kNTerms; ++k) {
        const double d = std::abs(t2[k] - 2.0 * t1[k]);
        const double tol = std::max(kTolE, 1e-12 * std::abs(t2[k]));
        if (d / tol > worst) { worst = d / tol; worst_d = d; wk = k; }
        check(d <= tol, std::string(label) + " supercell: " + kTermNames[k] + " = 2 x cell (" + sci(d) + " Eh)");
    }
    const Matrix top = t2.grad.topRows(n), bot = t2.grad.bottomRows(n);
    const double dg = std::max(maxAbs(top, t1.grad), maxAbs(bot, t1.grad));
    check(dg <= kTolG, std::string(label) + " supercell: per-atom gradient equals the cell (" + sci(dg) + " Eh/A)");
    check(t2.n_hb == 2 * t1.n_hb, std::string(label) + " supercell: HB count doubled");
    std::printf("%-46s supercell worst term %s (%.2e Eh), grad %.2e Eh/A, HB %d = 2 x %d\n", label, kTermNames[wk],
                worst_d, dg, t2.n_hb, t1.n_hb);
}

// (d) finite difference of the GPU energy on a torn geometry (lists fixed)
static void testFiniteDifference(const char* label, const Mol& base, const json& extra, int nfd)
{
    const double L = kCellBohr / A2B;
    const Variants v = makeVariants(base);
    Mol m = withCell(base, Eigen::Vector3d(L, L, L));
    m.m_geometry = v.torn;
    json e2 = extra;
    e2["hb_update_force_every"] = 0;   // fixed lists between the displaced points
    Instance gpu(true, m, config(e2));
    check(gpu.ok, std::string(label) + ": FD set-up");
    if (!gpu.ok)
        return;
    const Terms r = gpu.evaluate(v.torn, true);
    const int n = base.m_number_atoms;
    std::vector<std::pair<double, int>> dist;
    for (int i = 0; i < n; ++i) {
        double dmin = 1e300;
        for (int k = 0; k < 3; ++k)
            dmin = std::min({ dmin, v.torn(i, k), L - v.torn(i, k) });
        dist.push_back({ dmin, i });
    }
    std::sort(dist.begin(), dist.end());
    // central differences at h = 2e-4 and 1e-4 A; the truncation error is O(h^2) (ratio ~4), the Richardson
    // extrapolation (4 f(h) - f(2h))/3 removes it and is the number compared with the analytic gradient
    const double hs[2] = { 2e-4, 1e-4 };
    double err[2] = { 0, 0 }, err_r = 0.0, gmax = 0.0;
    for (int s = 0; s < nfd && s < n; ++s) {
        const int a = dist[s].second;
        for (int k = 0; k < 3; ++k) {
            double fd[2];
            for (int hi = 0; hi < 2; ++hi) {
                Geometry gp = v.torn, gm = v.torn;
                gp(a, k) += hs[hi];
                gm(a, k) -= hs[hi];
                const double ep = gpu.evaluate(gp, false).total, em = gpu.evaluate(gm, false).total;
                fd[hi] = (ep - em) / (2.0 * hs[hi]);
                err[hi] = std::max(err[hi], std::abs(fd[hi] - r.grad(a, k)));
            }
            err_r = std::max(err_r, std::abs((4.0 * fd[1] - fd[0]) / 3.0 - r.grad(a, k)));
            gmax = std::max(gmax, std::abs(r.grad(a, k)));
        }
    }
    check(err[1] <= 5e-8, std::string(label) + " FD: analytic vs central difference at h 1e-4 A (" + sci(err[1]) + " Eh/A)");
    check(err_r <= 2e-9, std::string(label) + " FD: Richardson-extrapolated difference (" + sci(err_r) + " Eh/A)");
    check(err[0] > 2.5 * err[1] || err[0] < 1e-10, std::string(label) + " FD: error quadratic in h (" + sci(err[0]) + " -> " + sci(err[1]) + ")");
    std::printf("%-46s FD (torn, %d atoms nearest a face): h 2e-4 %.2e, h 1e-4 %.2e, Richardson %.2e Eh/A, max |g| %.2e\n", label,
                std::min(nfd, n), err[0], err[1], err_r, gmax);
}

// (e) Verlet lists: a wrap by +L is no motion
static void testWrap(const char* label, const Mol& base, const json& extra)
{
    const double L = kCellBohr / A2B;
    const Variants v = makeVariants(base);
    Mol m = withCell(base, Eigen::Vector3d(L, L, L));
    m.m_geometry = v.mid;
    json e2 = extra;
    e2["nb_cell_list_min_atoms"] = 0;   // distance-filtered repulsion list (skin rule, device rebuild)
    Instance gpu(true, m, config(e2));
    check(gpu.ok, std::string(label) + ": wrap set-up");
    if (!gpu.ok)
        return;
    const Terms r0 = gpu.evaluate(v.mid, true);
    const long d0 = gpu.ff()->dispersionRebuildCount(), p0 = gpu.ff()->repulsionRebuildCount();
    const auto frags = fragmentsOf(base);
    Geometry g = v.mid;
    for (int a : frags[0])
        g(a, 0) += L;
    const Terms r1 = gpu.evaluate(g, true);
    const long d1 = gpu.ff()->dispersionRebuildCount(), p1 = gpu.ff()->repulsionRebuildCount();
    const Diff d = diffTerms(r1, r0, true);
    check(d1 == d0 && p1 == p0, std::string(label) + ": wrap by +L rebuilds no D4 / repulsion list (" + std::to_string(d1 - d0) + ", "
                                    + std::to_string(p1 - p0) + ")");
    check(d.e <= kTolE && d.g <= kTolG, std::string(label) + ": wrap by +L keeps terms and gradient (" + sci(d.e) + ", " + sci(d.g) + ")");
    // a real move rebuilds (control)
    Geometry g2 = v.mid;
    for (int a : frags[0])
        g2(a, 0) += 1.5;   // 1.5 A move of the first molecule
    gpu.evaluate(g2, true);
    const long p2 = gpu.ff()->repulsionRebuildCount();
    std::printf("%-46s wrap: D4 rebuilds %ld -> %ld, repulsion %ld -> %ld (1.5 A move: %ld), |dE| %.2e, |dG| %.2e\n", label,
                d0, d1, p0, p1, p2, d.e, d.g);
}

// (f) two device instances: bit for bit
static void testBitwise(const char* label, const Mol& base, const json& extra)
{
    const double L = kCellBohr / A2B;
    const Variants v = makeVariants(base);
    Mol m = withCell(base, Eigen::Vector3d(L, L, L));
    m.m_geometry = v.torn;
    Instance a(true, m, config(extra)), b(true, m, config(extra));
    check(a.ok && b.ok, std::string(label) + ": bitwise set-up");
    if (!a.ok || !b.ok)
        return;
    const Terms ta = a.evaluate(v.torn, true), tb = b.evaluate(v.torn, true);
    const bool same_e = std::memcmp(&ta.total, &tb.total, sizeof(double)) == 0;
    const bool same_g = ta.grad.rows() == tb.grad.rows()
        && std::memcmp(ta.grad.data(), tb.grad.data(), sizeof(double) * ta.grad.size()) == 0;
    check(same_e && same_g, std::string(label) + ": two device instances bit-identical (energy " + (same_e ? "yes" : "NO")
                                + ", gradient " + (same_g ? "yes" : "NO") + ")");
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const double L_ang = kCellBohr / A2B;
#ifdef GFNFF_PBC_TERMS_ROCM
    // A2: no HIP device (or no usable plugin state) -> skipped, as the CUDA build without CUDA; the probe has no cell
    {
        Mol probe = registryMol("H2O_dimer");
        Instance p(true, probe, config(json::object()));
        if (!p.ok) {
            std::printf("%s: the ROCm device workspace could not be set up (no HIP device?), skipped\n", kTestName);
            return 77;
        }
    }
#endif

    struct S { const char* label; Mol mol; json extra; };
    std::vector<S> set;
    set.push_back({ "water dimer", registryMol("H2O_dimer"), json::object() });
    set.push_back({ "caffeine", registryMol("caffeine"), json::object() });
    set.push_back({ "hydronium-water2 (+1)", registryMol("hydronium_water2", 1), json{ { "frag_charge_model", "reference" } } });
    set.push_back({ "acetic acid dimer", registryMol("acetic_acid_dimer"), json::object() });
    set.push_back({ "composite (acetic dimer + FI-pyridine)", compositeMol(), json::object() });
    set.push_back({ "40 waters", waterBox(40, L_ang, 7), json::object() });
    set.push_back({ "60 waters", waterBox(60, L_ang, 11), json::object() });
    set.push_back({ "60 waters (cell lists)", waterBox(60, L_ang, 11), json{ { "nb_cell_list_min_atoms", 0 } } });

    for (const auto& s : set)
        testStructure(s.label, s.mol, s.extra);
    for (const auto& s : set)
        if (s.mol.m_charge == 0)   // frag_charge_model reference puts the whole net charge on the first fragment: the two copies differ
            testSupercell(s.label, s.mol, s.extra);
    testFiniteDifference("caffeine", registryMol("caffeine"), json::object(), 8);
    testFiniteDifference("composite", compositeMol(), json::object(), 8);
    testFiniteDifference("40 waters", waterBox(40, L_ang, 7), json::object(), 4);
    testWrap("60 waters", waterBox(60, L_ang, 11), json::object());
    testBitwise("60 waters", waterBox(60, L_ang, 11), json::object());
    testBitwise("composite", compositeMol(), json::object());

    std::printf("GPU vs CPU over every comparison: worst term %.2e Eh, worst gradient %.2e Eh/A\n", g_worst_e, g_worst_g);
    for (int k = 0; k < kNTerms; ++k)
        std::printf("  %-14s max |dE| %.2e Eh\n", kTermNames[k], g_worst_term_e[k]);
    if (g_failures == 0)
        std::printf("%s: all checks passed\n", kTestName);
    else
        std::printf("%s: %d failure(s)\n", kTestName, g_failures);
    return g_failures == 0 ? 0 : 1;
}
#else
#include <cstdio>
int main()
{
    std::printf("gfnff_gpu_pbc_terms: built without CUDA, skipped\n");
    return 77;
}
#endif
