/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Device Ewald of the GFN-FF CUDA periodic hybrid (docs/GPU_PERIODIC_GAPS.md A2, stages S3/S4 of
 * docs/WP-PERIODIC-NONBONDED.md) - Claude Generated (Nov 2026).
 *
 * With -gfnff.gpu_ewald_device true (default of the hybrid) the Ewald pair sum (energy, gradient, virial at fixed charges)
 * and the operator of the periodic pcg EEQ run on the device (GfnffEwaldGPU); false keeps the S2 host Ewald.
 *
 * Checks (thresholds fixed before the first run):
 *  (a) device Ewald hybrid vs CPU periodic path and vs the host-Ewald hybrid, at the input geometry and after a move of
 *      every atom by up to 0.005 A (per-step solve): charges <= 1e-9 e, energy <= 1e-9 Eh, gradient <= 1e-9 Eh/A, model
 *      virial <= 1e-7 Eh per component; the device was used (operator calls and Ewald evaluations counted), no GPU
 *      fallback counted.
 *  (b) FD gate of the device-Ewald hybrid: total gradient vs central differences of the total energy (charges re-solved),
 *      h 2e-4 / 1e-4 A, Richardson, 6 atoms nearest a face: <= 1e-8 Eh/A; virial vs FD under isotropic and axial strain
 *      (fixed fractional coordinates, eps 1e-5 / 2e-5, Richardson): <= 1e-7 Eh.
 *  (c) a second GPU instance (two steps each): energy, gradient, charges and virial bit for bit.
 *  (d) fault injection CURCUMA_GFNFF_EWALD_DEVICE_FAIL=prepare|potentials|evaluate: the run continues on the host path
 *      (energy equal to the host-Ewald hybrid within 1e-10 Eh) and a GPU fallback is counted.
 * Extra structures: `test_gfnff_gpu_pbc_ewald --xyz FILE [FILE...]` runs (a) on XYZ files with a unitcell comment
 * (production settings, 4 threads) and prints the time of one gradient step on each path.
 * Exit 77 (ctest SKIP_RETURN_CODE) when built without CUDA.
 *
 * B1 (Nov 24, 2026, Claude Generated, ROCm parity, docs/ROCM_PARITY_B1.md): the same source builds the ROCm variant
 * (target test_gfnff_rocm_pbc_ewald, ctest gfnff_rocm_pbc_ewald, GFNFF_PBC_EWALD_ROCM): GFNFFHipComputationalMethod and
 * GfnffEwaldHip. The ROCm hybrid has no device virial yet (item B2 of ROCM_PARITY_PLAN_2026-11), so there the model
 * virial checks of (a)-(c) are not made (kHybridVirial false); instead (e) compares the device Ewald object itself
 * with the CPU EwaldSum range functions: real/reciprocal energy, gradient, potentials and the pair virial
 * (real + reciprocal) of the device evaluate(), limits 1e-10 Eh, 1e-10 Eh/Bohr, 1e-10 Eh/e, 1e-9 Eh.
 * Exit 77 there when no HIP device answers.
 */
#if defined(GFNFF_PBC_EWALD_ROCM) || defined(USE_CUDA)
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#ifdef GFNFF_PBC_EWALD_ROCM
#include "src/core/energy_calculators/ff_methods/rocm/ewald_hip.h"
#include "src/core/energy_calculators/qm_methods/gfnff_hip_method.h"
using GpuGfnffMethod = GFNFFHipComputationalMethod;
using DeviceEwald = GfnffEwaldHip;
static constexpr bool kHybridVirial = false;   // B2 of the ROCm parity plan
#else
#include "src/core/energy_calculators/ff_methods/cuda/ewald_gpu.h"
#include "src/core/energy_calculators/qm_methods/gfnff_gpu_method.h"
using GpuGfnffMethod = GFNFFGPUComputationalMethod;
using DeviceEwald = GfnffEwaldGPU;
static constexpr bool kHybridVirial = true;
#endif
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"
#include "src/core/fileiterator.h"
#include "src/core/gpu_fallback.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "src/core/elements.h"
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
    std::printf("  [%s] %-88s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), value, limit);
    if (!ok)
        ++g_failures;
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

/// Rigid random packing (the construction of test_gfnff_gpu_pbc_hybrid.cpp; same seeds, same boxes)
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

static json config(int threads, bool fd, int dev /* 1 device Ewald, 0 host Ewald */)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["gpu_periodic"] = true;
    c["gpu_ewald_device"] = dev == 1;
    c["frag_charge_model"] = "reference";
    if (fd) {
        c["hb_update_force_every"] = 1;
        c["hb_min_pair_energy_eh"] = 0.0;
        c["eeq_ewald_tol"] = 1e-14;
        c["ewald_tolerance"] = 1e-14;
        c["ewald_alpha"] = 0.32;
        c["ewald_r_cut"] = 19.0;
        c["dispersion_cutoff_bohr"] = 20.0;
    }
    return c;
}

enum class Path { CPU, GPU_HOST_EWALD, GPU_DEVICE_EWALD };
static const char* pathName(Path p) { return p == Path::CPU ? "CPU" : (p == Path::GPU_HOST_EWALD ? "GPU host Ewald" : "GPU device Ewald"); }

struct Result {
    bool ok = false;
    double e = 0.0;
    Matrix grad;
    Vector q;
    Eigen::Matrix3d vir = Eigen::Matrix3d::Zero();
    double ms = 0.0;   ///< wall time of the gradient call
};

/// One instance of a path; evaluate() moves it to a geometry and returns E, gradient (Eh/A), charges, model virial.
struct Instance {
    std::unique_ptr<ComputationalMethod> meth;
    GpuGfnffMethod* gpu = nullptr;
    bool ok = false;
    Instance(Path p, const Mol& m, const json& cfg)
    {
        try {
            if (p == Path::CPU) {
                auto c = std::make_unique<GFNFFComputationalMethod>("gfnff", cfg);
                ok = c->setMolecule(m) && !c->hasError();
                meth = std::move(c);
            } else {
                auto g = std::make_unique<GpuGfnffMethod>("gfnff", cfg);
                ok = g->setMolecule(m) && !g->hasError() && g->periodicHybridActive();
                if (ok && (p == Path::GPU_DEVICE_EWALD) != g->deviceEwaldActive()) {
                    std::printf("  device Ewald active = %d, expected %d\n", g->deviceEwaldActive(), p == Path::GPU_DEVICE_EWALD);
                    ok = false;
                }
                gpu = g.get();
                meth = std::move(g);
            }
        } catch (const std::exception& e) {
            std::printf("  set-up threw: %s\n", e.what());
            ok = false;
        }
    }
    Result evaluate(const Geometry& g_ang, bool virial, const Eigen::Matrix3d* cell = nullptr)
    {
        Result r;
        if (!ok)
            return r;
        try {
            GFNFF* ff = meth->gfnffInstance();
            if (cell)
                ff->updateCell(*cell);
            meth->updateGeometry(g_ang);
            ff->setWantVirial(virial);
            const auto t0 = std::chrono::steady_clock::now();
            r.e = meth->calculateEnergy(true);
            r.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            r.grad = meth->getGradient();
            r.q = meth->getCharges();
            if (virial) {
                if (!ff->virialValid() || !ff->virialComponents())
                    return r;
                r.vir = ff->virialComponents()->total();
            }
            r.ok = true;
        } catch (const std::exception& e) {
            std::printf("  evaluate threw: %s\n", e.what());
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
static int fallbackCount()
{
    int n = 0;
    for (const auto& c : curcuma::gpuFallbackCounts())
        n += c.second;
    return n;
}

static Geometry moved(const Mol& m, unsigned seed, double amp_ang)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-amp_ang, amp_ang);
    Geometry g = m.m_geometry;
    for (int i = 0; i < g.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            g(i, k) += u(rng);
    return g;
}

// (a) device Ewald vs CPU and vs the host-Ewald hybrid
/// The ionic box needed loose CPU limits until A16 (Nov 23, 2026): the CUDA CN ignored cn_cutoff_bohr, so Na-Cl pairs at
/// 10-11 Bohr entered C6 and the EEQ on the GPU only. Since the fix the device CN uses the CPU CN range and every box is
/// held to the strict limits (GPU_PERIODIC_GAPS.md A16).
static void testParity(const char* label, const Mol& m, int threads, bool fd)
{
    std::printf("(a) %s (%d atoms, L %.2f Bohr, %d threads)\n", label, m.m_number_atoms, m.m_unit_cell(0, 0) * A2B, threads);
    const int fb0 = fallbackCount();
    Instance dev(Path::GPU_DEVICE_EWALD, m, config(threads, fd, 1));
    Instance host(Path::GPU_HOST_EWALD, m, config(threads, fd, 0));
    Instance cpu(Path::CPU, m, config(threads, fd, 0));
    check(dev.ok && host.ok && cpu.ok, std::string(label) + ": three instances", dev.ok && host.ok && cpu.ok, 1);
    if (!dev.ok || !host.ok || !cpu.ok)
        return;
    const Geometry geo[2] = { m.m_geometry, moved(m, 777, 0.005) };
    const char* gname[2] = { "input geometry", "moved geometry (per-step solve)" };
    for (int s = 0; s < 2; ++s) {
        const Result d = dev.evaluate(geo[s], kHybridVirial), h = host.evaluate(geo[s], kHybridVirial), c = cpu.evaluate(geo[s], kHybridVirial);
        check(d.ok && h.ok && c.ok, std::string(label) + " " + gname[s] + ": evaluated", 1, 1);
        if (!d.ok || !h.ok || !c.ok)
            continue;
        std::printf("    %s: E dev %.12f host %.12f CPU %.12f; ms dev %.1f host %.1f CPU %.1f\n", gname[s], d.e, h.e, c.e, d.ms, h.ms, c.ms);
        for (int ref = 0; ref < 2; ++ref) {
            const Result& r = ref == 0 ? c : h;
            const std::string tag = std::string(label) + " " + gname[s] + " dev vs " + (ref == 0 ? "CPU" : "host-Ewald GPU");
            const double le = 1e-9, lg = 1e-9, lv = 1e-7;
            check(maxAbsV(d.q, r.q) <= 1e-9, tag + ": charges (e)", maxAbsV(d.q, r.q), 1e-9);
            check(std::abs(d.e - r.e) <= le, tag + ": energy (Eh)", std::abs(d.e - r.e), le);
            check(maxAbs(d.grad, r.grad) <= lg, tag + ": gradient (Eh/A)", maxAbs(d.grad, r.grad), lg);
            if (kHybridVirial) {
                const double dv = (d.vir - r.vir).cwiseAbs().maxCoeff();
                check(dv <= lv, tag + ": model virial per component (Eh)", dv, lv);
            }
        }
    }
    const auto* ew = dynamic_cast<const DeviceEwald*>(dev.gpu->deviceEwald());
    check(ew != nullptr, std::string(label) + ": device Ewald object present", ew != nullptr, 1);
    if (ew) {
        const auto& st = ew->stats();
        std::printf("    device stats: rebuilds %ld, key hits %ld, operator calls %ld, evaluations %ld, full pairs %lld\n", st.rebuilds,
            st.key_hits, st.potentials_calls, st.evaluate_calls, st.full_pairs);
        check(st.evaluate_calls >= 2, std::string(label) + ": Ewald energy evaluated on the device (calls)", static_cast<double>(st.evaluate_calls), 2);
        if (m.m_number_atoms > 1000)
            check(st.potentials_calls > 0, std::string(label) + ": pcg operator ran on the device (calls)", static_cast<double>(st.potentials_calls), 1);
    }
    check(fallbackCount() == fb0, std::string(label) + ": no GPU fallback counted", fallbackCount() - fb0, 0);
}

// (b) FD gates of the device-Ewald hybrid
static void testFD(const char* label, const Mol& m, bool force_pcg = false)
{
    std::printf("(b) FD gate, device Ewald: %s (%d atoms)%s\n", label, m.m_number_atoms, force_pcg ? ", pcg forced" : "");
    json cfg = config(1, true, 1);
    if (force_pcg)
        cfg["eeq_ewald_dense_max_atoms"] = 0;   // the matrix-free pcg with the device operator below 1000 atoms
    Instance g(Path::GPU_DEVICE_EWALD, m, cfg);
    check(g.ok, std::string(label) + ": instance", g.ok, 1);
    if (!g.ok)
        return;
    const Result an = g.evaluate(m.m_geometry, kHybridVirial);
    // 6 atoms nearest a face
    std::vector<std::pair<double, int>> face;
    const double L = m.m_unit_cell(0, 0);
    for (int i = 0; i < m.m_number_atoms; ++i) {
        double dmin = 1e9;
        for (int k = 0; k < 3; ++k) {
            double s = m.m_geometry(i, k) / L;
            s -= std::floor(s);
            dmin = std::min(dmin, std::min(s, 1.0 - s) * L);
        }
        face.push_back({ dmin, i });
    }
    std::sort(face.begin(), face.end());
    double worst = 0.0;
    for (int n = 0; n < 6; ++n) {
        const int i = face[n].second;
        for (int k = 0; k < 3; ++k) {
            double e[4];
            const double hs[4] = { -2e-4, -1e-4, 1e-4, 2e-4 };
            for (int s = 0; s < 4; ++s) {
                Geometry x = m.m_geometry;
                x(i, k) += hs[s];
                e[s] = g.evaluate(x, false).e;
            }
            const double fd1 = (e[2] - e[1]) / 2e-4, fd2 = (e[3] - e[0]) / 4e-4;
            const double rich = (4.0 * fd1 - fd2) / 3.0;
            worst = std::max(worst, std::abs(an.grad(i, k) - rich));
        }
    }
    check(worst <= 1e-8, std::string(label) + ": gradient vs Richardson FD, 6 atoms near a face (Eh/A)", worst, 1e-8);
    if (force_pcg) {
        const auto* ew = dynamic_cast<const DeviceEwald*>(g.gpu->deviceEwald());
        const long calls = ew ? ew->stats().potentials_calls : 0;
        check(calls > 0, std::string(label) + ": the FD energies used the device pcg operator (calls)", static_cast<double>(calls), 1);
    }
    if (!kHybridVirial)
        return;   // ROCm: no device virial of the hybrid (B2)
    // virial under strain (fresh instance per strained cell)
    const double eps[4] = { -2e-5, -1e-5, 1e-5, 2e-5 };
    for (int axis = -1; axis < 3; ++axis) {
        double e[4];
        bool ok = true;
        for (int s = 0; s < 4; ++s) {
            Mol st = m;
            for (int k = 0; k < 3; ++k)
                if (axis < 0 || axis == k) {
                    st.m_geometry.col(k) *= (1.0 + eps[s]);
                    st.m_unit_cell(k, k) *= (1.0 + eps[s]);
                }
            Instance gs(Path::GPU_DEVICE_EWALD, st, cfg);
            const Result r = gs.evaluate(st.m_geometry, false);
            ok = ok && r.ok;
            e[s] = r.e;
        }
        if (!ok) {
            check(false, std::string(label) + ": strained instances", 0, 1);
            return;
        }
        const double fd1 = -(e[2] - e[1]) / 2e-5, fd2 = -(e[3] - e[0]) / 4e-5;
        const double rich = (4.0 * fd1 - fd2) / 3.0;
        const double xan = axis < 0 ? an.vir.trace() : an.vir(axis, axis);
        const char* nm[4] = { "isotropic (trace)", "eps_xx", "eps_yy", "eps_zz" };
        check(std::abs(xan - rich) <= 1e-7, std::string(label) + " virial FD " + nm[axis + 1] + " (Eh)", std::abs(xan - rich), 1e-7);
    }
}

// (c) bitwise repeat
static void testBitwise(const char* label, const Mol& m, int threads)
{
    std::printf("(c) second GPU instance, device Ewald: %s\n", label);
    const json cfg = config(threads, false, 1);
    Instance a(Path::GPU_DEVICE_EWALD, m, cfg), b(Path::GPU_DEVICE_EWALD, m, cfg);
    if (!a.ok || !b.ok) {
        check(false, std::string(label) + ": instances", 0, 1);
        return;
    }
    bool same = true;
    for (int s = 0; s < 2; ++s) {
        const Geometry x = moved(m, 900 + s, 0.003);
        const Result ra = a.evaluate(x, kHybridVirial), rb = b.evaluate(x, kHybridVirial);
        same = same && ra.ok && rb.ok && ra.e == rb.e && ra.grad.size() == rb.grad.size()
            && std::memcmp(ra.grad.data(), rb.grad.data(), ra.grad.size() * sizeof(double)) == 0 && ra.q.size() == rb.q.size()
            && std::memcmp(ra.q.data(), rb.q.data(), ra.q.size() * sizeof(double)) == 0
            && std::memcmp(ra.vir.data(), rb.vir.data(), 9 * sizeof(double)) == 0;
    }
    check(same, std::string(label) + ": energy, gradient, charges, virial bit for bit over two steps", same ? 0 : 1, 0);
}

// (d) fault injection
static void testFallback(const Mol& m)
{
    std::printf("(d) fault injection (%d atoms)\n", m.m_number_atoms);
    const json cfg_host = config(4, false, 0);
    Instance host(Path::GPU_HOST_EWALD, m, cfg_host);
    const Geometry x = moved(m, 4242, 0.004);
    const Result rh = host.evaluate(x, kHybridVirial);
    for (const char* mode : { "prepare", "potentials", "evaluate" }) {
        setenv("CURCUMA_GFNFF_EWALD_DEVICE_FAIL", mode, 1);
        const int fb0 = fallbackCount();
        Instance dev(Path::GPU_DEVICE_EWALD, m, config(4, false, 1));
        const Result rd = dev.evaluate(x, kHybridVirial);
        unsetenv("CURCUMA_GFNFF_EWALD_DEVICE_FAIL");
        const int fb = fallbackCount() - fb0;
        check(rd.ok && rh.ok && std::abs(rd.e - rh.e) <= 1e-10, std::string("(d) ") + mode + ": energy equals the host-Ewald hybrid (Eh)",
            rd.ok && rh.ok ? std::abs(rd.e - rh.e) : 1e300, 1e-10);
        check(fb > 0, std::string("(d) ") + mode + ": GPU fallback counted", fb, 1);
    }
}

// (e) the device Ewald object alone against the CPU EwaldSum range functions (the CPU workspace's sum, widths and
// positions): real/reciprocal energy, gradient, potentials and the real + reciprocal pair virial. Claude Generated (B1).
[[maybe_unused]] static void testDeviceOps(const char* label, const Mol& m)
{
    std::printf("(e) device Ewald object vs CPU EwaldSum: %s (%d atoms)\n", label, m.m_number_atoms);
    Instance cpu(Path::CPU, m, config(1, false, 0));
    FFWorkspace* ws = cpu.ok ? cpu.meth->gfnffInstance()->getWorkspace() : nullptr;
    check(ws && ws->ewald(), std::string(label) + ": CPU Ewald set-up", ws && ws->ewald(), 1);
    if (!ws || !ws->ewald())
        return;
    const curcuma::ewald::EwaldSum& ew = *ws->ewald();
    const Eigen::MatrixXd pos = m.m_geometry * A2B;
    const Eigen::VectorXd alp = ws->coulombAlp();
    const int N = m.m_number_atoms;
    std::mt19937 rng(31);
    std::uniform_real_distribution<double> u(-0.6, 0.6);
    Eigen::VectorXd q(N);
    for (int i = 0; i < N; ++i) q(i) = u(rng);
    q.array() -= q.mean();
    // CPU reference: stored-kernel cache, then the range functions over everything
    curcuma::ewald::EwaldGeometryCache cache;
    ew.buildGeometryCache(pos, &alp, cache, size_t(1) << 32);
    const curcuma::ewald::EwaldPhaseTable table = ew.phaseTable(pos);
    curcuma::ewald::EwaldResult acc;
    acc.gradient = Eigen::MatrixXd::Zero(N, 3);
    acc.potentials = Eigen::VectorXd::Zero(N);
    ew.realSpaceRangeCached(cache, q, 0, static_cast<int>(cache.pairs.size()), true, true, acc);
    ew.reciprocalRange(table, q, 0, ew.numberOfKVectors(), true, true, acc);
    Eigen::Matrix3d vir = Eigen::Matrix3d::Zero();
    ew.realSpaceVirialRange(cache, q, 0, static_cast<int>(cache.pairs.size()), vir);
    ew.reciprocalVirialRange(table, q, 0, ew.numberOfKVectors(), vir);
    // device
    DeviceEwald dev;
    Eigen::VectorXd y;
    Eigen::MatrixXd g;
    Eigen::Matrix3d v;
    curcuma::ewald::EwaldEnergyTerms t;
    const bool ok = dev.prepare(ew, pos, alp) && dev.potentials(q, y) && dev.evaluate(q, t, &g, &v);
    check(ok, std::string(label) + ": device prepare/potentials/evaluate (" + dev.lastError() + ")", ok, 1);
    if (!ok)
        return;
    check(std::abs(t.real - acc.terms.real) <= 1e-10, std::string(label) + ": real-space energy (Eh)", std::abs(t.real - acc.terms.real), 1e-10);
    check(std::abs(t.reciprocal - acc.terms.reciprocal) <= 1e-10, std::string(label) + ": reciprocal energy (Eh)",
        std::abs(t.reciprocal - acc.terms.reciprocal), 1e-10);
    check(maxAbs(g, acc.gradient) <= 1e-10, std::string(label) + ": pair gradient (Eh/Bohr)", maxAbs(g, acc.gradient), 1e-10);
    check(maxAbsV(y, acc.potentials) <= 1e-10, std::string(label) + ": potentials A q (Eh/e)", maxAbsV(y, acc.potentials), 1e-10);
    const double dv = (v - vir).cwiseAbs().maxCoeff();
    check(dv <= 1e-9, std::string(label) + ": real + reciprocal pair virial per component (Eh)", dv, 1e-9);
    // second evaluate on the same key: bit for bit
    Eigen::MatrixXd g2;
    Eigen::Matrix3d v2;
    curcuma::ewald::EwaldEnergyTerms t2;
    const bool ok2 = dev.prepare(ew, pos, alp) && dev.evaluate(q, t2, &g2, &v2);
    const bool same = ok2 && t2.real == t.real && t2.reciprocal == t.reciprocal
        && std::memcmp(g2.data(), g.data(), g.size() * sizeof(double)) == 0 && std::memcmp(v2.data(), v.data(), 9 * sizeof(double)) == 0;
    check(same, std::string(label) + ": repeated evaluate bit for bit (key hit)", same ? 0 : 1, 0);
    check(dev.stats().key_hits >= 1, std::string(label) + ": second prepare was a key hit", static_cast<double>(dev.stats().key_hits), 1);
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
#ifdef GFNFF_PBC_EWALD_ROCM
    {   // B1: no HIP device -> skipped (77), as the CUDA build without CUDA
        DeviceEwald probe;
        if (!probe.usable()) {
            std::printf("gfnff_rocm_pbc_ewald: no HIP stream (no HIP device?), skipped\n");
            return 77;
        }
    }
#endif
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
            testParity(argv[a], m, 4, false);
        }
        std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
        return g_failures == 0 ? 0 : 1;
    }
    if (argc > 2 && std::strcmp(argv[1], "--bench") == 0) {   // diagnostic: device Ewald alone (no other GPU work), ms per call
        FileIterator it(argv[2], true);
        Mol m = it.Next().getMolInfo();
        Instance cpu(Path::CPU, m, config(4, false, 0));
        FFWorkspace* ws = cpu.ok ? cpu.meth->gfnffInstance()->getWorkspace() : nullptr;
        if (!ws || !ws->ewald())
            return 1;
        const curcuma::ewald::EwaldSum& ew = *ws->ewald();
        const Eigen::MatrixXd pos = m.m_geometry * A2B;
        const Eigen::VectorXd alp = ws->coulombAlp();
        DeviceEwald dev;
        Eigen::VectorXd x = Eigen::VectorXd::Random(m.m_number_atoms), y;
        Eigen::MatrixXd g;
        Eigen::Matrix3d v;
        curcuma::ewald::EwaldEnergyTerms t;
        for (int rep = 0; rep < 5; ++rep) {
            Eigen::MatrixXd p = pos;
            p(0, 0) += 1e-6 * (rep + 1);   // new key: full rebuild
            auto t0 = std::chrono::steady_clock::now();
            const bool ok1 = dev.prepare(ew, p, alp);
            auto t1 = std::chrono::steady_clock::now();
            const bool ok2 = dev.potentials(x, y);
            auto t2 = std::chrono::steady_clock::now();
            const bool ok3 = dev.evaluate(x, t, &g, &v);
            auto t3 = std::chrono::steady_clock::now();
            std::printf("bench N=%d prepare %.3f ms, potentials %.3f ms, evaluate(grad, virial) %.3f ms %s\n", m.m_number_atoms,
                std::chrono::duration<double, std::milli>(t1 - t0).count(), std::chrono::duration<double, std::milli>(t2 - t1).count(),
                std::chrono::duration<double, std::milli>(t3 - t2).count(), ok1 && ok2 && ok3 ? "" : "FAILED");
        }
        return 0;
    }
    if (argc > 2 && std::strcmp(argv[1], "--write-salt") == 0) {   // diagnostic: the ionic box as XYZ with a unitcell line
        const Mol salt = waterBox(380, 48.0, 20261123, { registryMol("na2cl2_ions") }, 2.4);
        FILE* f = std::fopen(argv[2], "w");
        if (!f)
            return 1;
        const double L = salt.m_unit_cell(0, 0);
        std::fprintf(f, "%d\nunitcell %.12f %.12f %.12f 90.0 90.0 90.0\n", salt.m_number_atoms, L, L, L);
        for (int i = 0; i < salt.m_number_atoms; ++i)
            std::fprintf(f, "%s %.10f %.10f %.10f\n", Elements::ElementAbbr[salt.m_atoms[i]].c_str(), salt.m_geometry(i, 0),
                salt.m_geometry(i, 1), salt.m_geometry(i, 2));
        std::fclose(f);
        return 0;
    }
    const char* only = std::getenv("GPU_EWALD_ONLY");   // diagnostic: a/b/c/d/e
    auto run = [&](char c) { return !only || std::strchr(only, c); };
    const Mol box20 = waterBox(20, 44.0, 20261105);
    const Mol box60 = waterBox(60, 44.0, 20261108);
    const Mol charged = waterBox(18, 44.0, 20261106, { registryMol("hydronium_water2", 1) });
    if (run('a')) testParity("20 H2O (FD protocol)", box20, 1, true);
    if (run('a')) testParity("60 H2O", box60, 4, false);
    if (run('a')) testParity("H3O+(H2O)2 + 18 H2O, Q = +1", charged, 1, true);
    if (run('a') || run('c') || run('d')) {
        const Mol box970 = waterBox(970, 58.0, 20261105);
        if (run('a')) testParity("970 H2O L 58 (pcg)", box970, 4, false);
        if (run('c')) testBitwise("970 H2O L 58 (pcg)", box970, 4);
        if (run('d')) testFallback(box970);
    }
    if (run('a')) {
        // ionic mixture: Na2Cl2 (four ions, neutral in total) in 380 waters, L 48 Bohr (0.7 g/cm3), 1144 atoms (matrix-free pcg)
        const Mol salt = waterBox(380, 48.0, 20261123, { registryMol("na2cl2_ions") }, 2.4);
        testParity("Na2Cl2 + 380 H2O L 48 (pcg)", salt, 4, false);
    }
#ifdef GFNFF_PBC_EWALD_ROCM   // (e) is the ROCm stand-in for the hybrid virial checks (not run on CUDA: not verified there)
    if (run('e')) testDeviceOps("20 H2O", box20);
    if (run('e')) testDeviceOps("H3O+(H2O)2 + 18 H2O, Q = +1", charged);
    if (run('e')) testDeviceOps("970 H2O L 58", waterBox(970, 58.0, 20261105));
#endif
    if (run('b')) testFD("20 H2O", box20);
    if (run('b')) testFD("H3O+(H2O)2 + 18 H2O, Q = +1", charged);
    if (run('b')) testFD("20 H2O", box20, true);
    if (run('b')) testFD("H3O+(H2O)2 + 18 H2O, Q = +1", charged, true);
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
