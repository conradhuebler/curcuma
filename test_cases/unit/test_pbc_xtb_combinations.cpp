/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic native GFN1/GFN2: combinations with a unit cell that were never run together - Claude Generated
 * (Nov 23, 2026, docs/GPU_PERIODIC_GAPS.md A6, A9).
 *
 * Argument "refusals" (ctest pbc_xtb_cell_refusals), "cuda" (ctest pbc_xtb_gpu_combinations) or "rocm" (ctest
 * pbc_xtb_gpu_combinations_rocm, Nov 24, 2026: part 2 on the ROCm plugin; the multi-GPU split and GFN-FF checks are CUDA only).
 * Part 1 (any build, no GPU needed): refusals with a cell, each must set an error instead of returning an energy:
 *   spin polarisation (-xtb.spin_polarization true), implicit solvation (-xtb.solvent water), CPCM-X
 *   (-xtb.solvent_model cpcmx), on the CPU and with -gpu cuda; -gpu vulkan with a cell: either a counted CPU fallback
 *   with the CPU periodic energy or an error, never an uncounted result; -gpu rocm with a cell (device path since
 *   Nov 24, 2026, ROCm parity A3): error, counted CPU fallback, or the CPU energy within 1e-8 Eh;
 *   invalid pbc_* values are an error.
 * Part 2 (CUDA plugin and device, else exit code 77 = ctest SKIP):
 *   open shell (OH radical, spin 1, registry "OH") and a charged cell (registry "hydronium_water2", charge +1) in an
 *   8/10 A cube on the host-driven CUDA loop against the CPU: |dE| <= 1e-8 Eh, max |dG| <= 1e-8 Eh/A,
 *   max |dXi| <= 1e-7 Eh, no stage fallback counted; an explicit multi-GPU split request with a cell is refused.
 * scf_threshold 1e-9. Exit 0 when every check passes.
 */
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"
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

static int fallbackTotal()
{
    int n = 0;
    for (const auto& kv : curcuma::gpuFallbackCounts()) n += kv.second;
    return n;
}
static int stageTotal()
{
    int n = 0;
    for (const auto& kv : curcuma::gpuFallbackCounts())
        if (kv.first.rfind("xTB periodic device stage ", 0) == 0) n += kv.second;
    return n;
}

static Mol inCell(const char* name, double L, int charge, double spin)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
    Mol mol = m.getMolInfo();
    mol.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    mol.m_has_pbc = true;
    mol.m_charge = charge;
    mol.m_spin = spin;
    return mol;
}

struct Result {
    double e = 0.0;
    Matrix g;
    Eigen::Matrix3d xi = Eigen::Matrix3d::Zero();
    bool ok = false;
};

static json baseCfg(const char* gpu)
{
    json cfg = { { "verbosity", 0 }, { "threads", 1 }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    if (gpu) cfg["gpu"] = gpu;
    return cfg;
}

static Result run(const char* method, const Mol& mol, const json& cfg)
{
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    Result r;
    if (calc.Error()) return r;
    calc.setWantVirial(true);
    r.e = calc.CalculateEnergy(true);
    r.ok = !calc.Error() && calc.virialValid();
    r.g = calc.Gradient();
    r.xi = calc.virial();
    return r;
}

/// true when setMolecule (or the first energy) of this configuration reports an error
static bool refused(const char* method, const Mol& mol, const json& cfg)
{
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    if (calc.Error()) return true;
    calc.CalculateEnergy(false);
    return calc.Error();
}

/// the refusal is caused by the cell: the same configuration runs without it
static bool refusedForCell(const char* method, Mol mol, const json& cfg)
{
    const bool with_cell = refused(method, mol, cfg);
    mol.m_has_pbc = false;
    mol.m_unit_cell = Eigen::Matrix3d::Zero();
    const bool without = refused(method, mol, cfg);
    if (without) std::printf("  note: %s refused also without a cell (configuration error, not the cell)\n", method);
    return with_cell && !without;
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    const std::string part = argc > 1 ? argv[1] : "refusals";
    const Mol dimer = inCell("H2O_dimer", 8.0, 0, 0.0);
    if (part != "refusals" && part != "cuda" && part != "rocm") {
        std::printf("usage: test_pbc_xtb_combinations refusals|cuda|rocm\n");
        return 2;
    }

    // ---- Part 1: refusals --------------------------------------------------------------------------------
    if (part == "refusals") {
    for (const char* method : { "gfn1", "gfn2" }) {
        for (const char* gpu : { static_cast<const char*>(nullptr), "cuda" }) {
            const std::string tag = std::string(method) + (gpu ? " -gpu cuda" : " CPU");
            json sp = baseCfg(gpu);
            sp["xtb"]["spin_polarization"] = true;
            check(refusedForCell(method, inCell("OH", 8.0, 0, 1.0), sp), tag + ": spin polarisation with a cell is refused");
            json solv = baseCfg(gpu);
            solv["xtb"]["solvent"] = "water";
            check(refusedForCell(method, dimer, solv), tag + ": implicit solvation with a cell is refused");
            json cpx = baseCfg(gpu);
            cpx["xtb"]["solvent_model"] = "cpcmx";
            cpx["xtb"]["solvent"] = "water";
            if (const char* db = std::getenv("CURCUMA_TEST_CPCMX_DB")) cpx["xtb"]["cpcmx_db"] = db;   // CPCM-X 1.1.0 water.cosmo
            // CPCM-X with -gpu is refused with and without a cell ("runs on the CPU only"), so only the CPU run tests the cell
            check(gpu ? refused(method, dimer, cpx) : refusedForCell(method, dimer, cpx), tag + ": CPCM-X with a cell is refused");
        }
    }
    std::printf("part 1: spin polarisation, solvation, CPCM-X with a cell refused (gfn1/gfn2 x CPU/-gpu cuda)\n");

    // P4 (Nov 23, 2026): an invalid pbc_* value is an error, no longer replaced silently by the default
    for (const auto& kv : { std::pair<const char*, double>{ "pbc_ewald_tolerance", 2.0 }, { "pbc_coulomb_switch_bohr", 0.5 },
                            { "pbc_integral_cutoff_bohr", -1.0 }, { "pbc_ewald_r_cut_bohr", -3.0 } }) {
        json bad = baseCfg(nullptr);
        bad["xtb"][kv.first] = kv.second;
        check(refused("gfn2", dimer, bad), std::string("invalid ") + kv.first + " = " + std::to_string(kv.second) + " is an error");
    }

    // ROCm / Vulkan with a cell: counted CPU fallback with the CPU energy, or an error
    const Result cpu_dimer = run("gfn1", dimer, baseCfg(nullptr));
    check(cpu_dimer.ok, "gfn1 CPU periodic reference run");
    for (const char* gpu : { "rocm", "vulkan" }) {
        const int fb0 = fallbackTotal();
        EnergyCalculator calc("gfn1", baseCfg(gpu));
        calc.setMolecule(dimer);
        double e = 0.0;
        bool err = calc.Error();
        if (!err) {
            e = calc.CalculateEnergy(false);
            err = calc.Error();
        }
        const int counted = fallbackTotal() - fb0;
        std::printf("-gpu %s with a cell: %s, %d fallback(s) counted, dE vs CPU %.2e Eh\n", gpu,
                    err ? "error" : "energy", counted, err ? 0.0 : e - cpu_dimer.e);
        // ROCm (A3, Nov 24, 2026): the device path runs with a cell; its energy must be the CPU one to 1e-8 Eh
        const bool rocm_device = std::string(gpu) == "rocm" && !err && std::abs(e - cpu_dimer.e) <= 1.0e-8;
        check(err || rocm_device || (counted >= 1 && std::abs(e - cpu_dimer.e) <= 1.0e-10),
              std::string("-gpu ") + gpu + " with a cell: honest error, counted CPU fallback with the CPU energy"
                  + (std::string(gpu) == "rocm" ? " or the ROCm device energy" : ""));
    }
    std::printf("test_pbc_xtb_combinations refusals: %s\n", g_failures == 0 ? "all checks passed" : "FAILURES");
    return g_failures == 0 ? 0 : 1;
    }

    // ---- Part 2: CUDA (or ROCm) parity of open shell and charged cells ----------------------------------------
    const char* backend = part.c_str();   // "cuda" or "rocm"
    const bool is_cuda = part == "cuda";
    bool any_gpu = false;
    {
        const int fb0 = fallbackTotal();
        const Result probe = run("gfn1", dimer, baseCfg(backend));
        any_gpu = probe.ok && fallbackTotal() == fb0;
    }
    if (any_gpu) {
        struct Case { const char* name; double L; int charge; double spin; };
        for (const char* method : { "gfn1", "gfn2" }) {
            for (const Case& c : { Case{ "OH", 8.0, 0, 1.0 }, Case{ "hydronium_water2", 10.0, 1, 0.0 } }) {
                const Mol mol = inCell(c.name, c.L, c.charge, c.spin);
                const Result cpu = run(method, mol, baseCfg(nullptr));
                const int st0 = stageTotal();
                const Result gpu = run(method, mol, baseCfg(backend));
                const Result gpu2 = run(method, mol, baseCfg(backend));
                const std::string tag = std::string(method) + " " + c.name + " (charge " + std::to_string(c.charge)
                    + ", spin " + std::to_string(static_cast<int>(c.spin)) + ")";
                check(cpu.ok && gpu.ok && gpu2.ok, tag + ": CPU and " + part + " runs succeed");
                if (!(cpu.ok && gpu.ok && gpu2.ok)) continue;
                const double de = std::abs(gpu.e - cpu.e);
                const double dg = (gpu.g - cpu.g).cwiseAbs().maxCoeff();
                const double dx = (gpu.xi - cpu.xi).cwiseAbs().maxCoeff();
                const bool same = gpu2.e == gpu.e && (gpu2.g - gpu.g).cwiseAbs().maxCoeff() == 0.0
                    && (gpu2.xi - gpu.xi).cwiseAbs().maxCoeff() == 0.0;
                std::printf("%s L = %.0f: E cpu %.12f %s %.12f |dE| %.2e Eh, max|dG| %.2e Eh/A, max|dXi| %.2e Eh, "
                            "stage fallbacks %d, second run %s\n", tag.c_str(), c.L, cpu.e, backend, gpu.e, de, dg, dx,
                            stageTotal() - st0, same ? "bitwise" : "DIFFERS");
                check(de <= 1.0e-8, tag + ": " + part + " vs CPU energy <= 1e-8 Eh");
                check(dg <= 1.0e-8, tag + ": " + part + " vs CPU gradient <= 1e-8 Eh/A");
                check(dx <= 1.0e-7, tag + ": " + part + " vs CPU virial <= 1e-7 Eh");
                check(stageTotal() == st0, tag + ": no periodic stage fallback");
                check(same, tag + ": second " + part + " run bitwise");
            }
        }
        if (is_cuda) {
        // explicit multi-GPU split with a cell: refused (never run with a cell); without a cell still accepted
        json split = baseCfg("cuda");
        split["gpu_eigensolver_devices"] = "0,1";
        check(refused("gfn2", dimer, split), "explicit -gpu_eigensolver_devices with a cell is refused");
        json dsplit = baseCfg("cuda");
        dsplit["gpu_density_devices"] = "0,1";
        check(refused("gfn2", dimer, dsplit), "explicit -gpu_density_devices with a cell is refused");
        Mol vac = dimer;
        vac.m_has_pbc = false;
        vac.m_unit_cell = Eigen::Matrix3d::Zero();
        check(!refused("gfn2", vac, split), "explicit -gpu_eigensolver_devices without a cell still runs");
        // A11: GFN-FF term-scale factors != 1 under -gpu cuda run on the CPU (counted fallback in createNativeGfnff; needs
        // the plugin, without it -gpu cuda resolves to the CPU before that check); the energy equals the CPU run
        {
            Mol vac = dimer;
            vac.m_has_pbc = false;
            vac.m_unit_cell = Eigen::Matrix3d::Zero();
            for (const char* key : { "hb_scale", "repulsion_scale", "dispersion_scale", "coulomb_scale" }) {
                json cpu_cfg = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", { { key, 1.5 } } } };
                json gpu_cfg = cpu_cfg;
                gpu_cfg["gpu"] = "cuda";
                EnergyCalculator c_cpu("gfnff", cpu_cfg);
                c_cpu.setMolecule(vac);
                const double e_cpu = c_cpu.CalculateEnergy(false);
                int scale_fb0 = 0;
                for (const auto& kv : curcuma::gpuFallbackCounts())
                    if (kv.first.find("term-scale") != std::string::npos) scale_fb0 += kv.second;
                EnergyCalculator c_gpu("gfnff", gpu_cfg);
                c_gpu.setMolecule(vac);
                const double e_gpu = c_gpu.CalculateEnergy(false);
                int scale_fb = -scale_fb0;
                for (const auto& kv : curcuma::gpuFallbackCounts())
                    if (kv.first.find("term-scale") != std::string::npos) scale_fb += kv.second;
                std::printf("gfnff -gfnff.%s 1.5 -gpu cuda: %d term-scale fallback(s), E %.12f vs CPU %.12f\n", key, scale_fb, e_gpu, e_cpu);
                check(scale_fb == 1, std::string("gfnff ") + key + " != 1 under -gpu cuda: one counted CPU fallback");
                check(!c_gpu.Error() && e_gpu == e_cpu, std::string("gfnff ") + key + " != 1 under -gpu cuda: CPU energy (bitwise)");
            }
        }
        }   // is_cuda
    } else {
        std::printf("part 2: SKIP (no %s plugin or device)\n", backend);
    }

    if (g_failures == 0)
        std::printf("test_pbc_xtb_combinations %s: %s\n", backend, any_gpu ? "all checks passed" : "SKIP");
    else
        std::printf("test_pbc_xtb_combinations %s: %d failures\n", backend, g_failures);
    if (g_failures != 0) return 1;
    return any_gpu ? 0 : 77;
}
