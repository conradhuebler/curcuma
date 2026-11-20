/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * CPU / CUDA parity of the periodic native GFN1/GFN2 (WP-PERIODIC-XTB G1) - Claude Generated (Nov 2026).
 *
 * G1: with a unit cell the host builds the periodic S/H0/L, gamma and multipole integrals and uploads them; the
 * CUDA device runs the resident Fock/eigensolve/density/moments; potential, dispersion, gradient and virial stay
 * on the host. Gate: energy, gradient and virial of `-gpu cuda` against the CPU on the registry water dimer in an
 * 8 A and a 12 A cube, the water trimer in a 6.5 A cube and the 216 H2O box (G2-G6), gfn1 and gfn2, scf_threshold 1e-9: |dE| <= 1e-8 Eh,
 * max |dG| <= 1e-8 Eh/A, max |dXi| <= 1e-7 Eh (coordinator targets); a second CUDA run from a fresh calculator must be bitwise identical (energy and gradient,
 * docs/GPU_REPRODUCIBILITY.md). The GPU path counts as used when the run added no GPU fallback
 * (curcuma::gpuFallbackCounts()); without the plugin or a device the first GPU request falls back, the test prints
 * SKIP and exits 0, so it never false-fails off-GPU. Exit 0 when every check passes.
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

static Mol inCell(const char* name, double L)
{
    curcuma::Molecule m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
    Mol mol = m.getMolInfo();
    mol.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    mol.m_has_pbc = true;
    return mol;
}

struct Result {
    double e = 0.0;
    Matrix g;
    Eigen::Matrix3d xi = Eigen::Matrix3d::Zero();
    bool ok = false;
};

static Result run(const char* method, const Mol& mol, bool gpu, int threads)
{
    json cfg = { { "verbosity", 0 }, { "threads", threads }, { "xtb", { { "scf_threshold", 1.0e-9 } } } };
    if (gpu) cfg["gpu"] = "cuda";
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    calc.setWantVirial(true);
    Result r;
    r.e = calc.CalculateEnergy(true);
    r.ok = !calc.Error() && calc.virialValid();
    r.g = calc.Gradient();
    r.xi = calc.virial();
    return r;
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    bool any_gpu = false;
    for (const char* method : { "gfn1", "gfn2" }) {
        // water dimer in 8 and 12 A, water trimer ("H2O") in a dense 6.5 A cube (G2-G6: many self images and
        // image triples in the ATM, multipole and Ewald paths); 216 H2O box (G2: many k vectors and a real-space
        // cutoff below L/2, the size where the device kernels are blocked; 4 host threads)
        struct Cell { const char* name; double L; int threads; };
        for (const Cell& c : { Cell{ "H2O_dimer", 8.0, 1 }, Cell{ "H2O_dimer", 12.0, 1 }, Cell{ "H2O", 6.5, 1 },
                               Cell{ "water216_box", 18.625793445945, 4 } }) {
            const double L = c.L;
            const Mol mol = inCell(c.name, L);
            const Result cpu = run(method, mol, false, c.threads);
            check(cpu.ok, std::string(method) + " CPU periodic energy/gradient/virial");
            const int fb0 = fallbackTotal();
            const Result gpu = run(method, mol, true, c.threads);
            if (fallbackTotal() != fb0) {
                std::printf("SKIP %s %s L = %.1f A: the -gpu cuda request fell back to the CPU (no plugin or device)\n", method, c.name, L);
                continue;
            }
            any_gpu = true;
            check(gpu.ok, std::string(method) + " CUDA periodic energy/gradient/virial");
            const double de = std::abs(gpu.e - cpu.e);
            const double dg = (gpu.g - cpu.g).cwiseAbs().maxCoeff();
            const double dx = (gpu.xi - cpu.xi).cwiseAbs().maxCoeff();
            std::printf("%s %s L = %.1f A: E cpu %.12f cuda %.12f |dE| %.3e Eh, max |dG| %.3e Eh/A, max |dXi| %.3e Eh\n",
                        method, c.name, L, cpu.e, gpu.e, de, dg, dx);
            check(de <= 1.0e-8, std::string(method) + " CUDA vs CPU energy <= 1e-8 Eh, L = " + std::to_string(L));
            check(dg <= 1.0e-8, std::string(method) + " CUDA vs CPU gradient <= 1e-8 Eh/A, L = " + std::to_string(L));
            check(dx <= 1.0e-7, std::string(method) + " CUDA vs CPU virial <= 1e-7 Eh, L = " + std::to_string(L));
            // run-to-run reproducibility of the CUDA path (fresh calculator)
            const Result gpu2 = run(method, mol, true, c.threads);
            const bool same = gpu2.ok && gpu2.e == gpu.e && (gpu2.g - gpu.g).cwiseAbs().maxCoeff() == 0.0
                && (gpu2.xi - gpu.xi).cwiseAbs().maxCoeff() == 0.0;
            std::printf("%s %s L = %.1f A: second CUDA run %s (dE %.1e)\n", method, c.name, L, same ? "bitwise identical" : "DIFFERS", gpu2.e - gpu.e);
            check(same, std::string(method) + " CUDA run-to-run bitwise, L = " + std::to_string(L));
        }
    }
    if (!any_gpu) {
        std::printf("test_pbc_xtb_gpu_parity: SKIP (no CUDA path ran)\n");
        return 0;
    }
    if (g_failures == 0)
        std::printf("test_pbc_xtb_gpu_parity: all checks passed\n");
    else
        std::printf("test_pbc_xtb_gpu_parity: %d failures\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
