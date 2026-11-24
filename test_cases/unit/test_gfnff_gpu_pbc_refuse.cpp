/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * S0 defence in depth (Claude Generated, Nov 2026): the CUDA GFN-FF wrapper itself refuses a periodic cell.
 *
 * EnergyCalculator and the delta-ML base swap a GPU GFN-FF with a cell for the CPU method before setMolecule()
 * (MethodFactory::periodicCellNeedsCpu, ctest pbc_gpu_guard). This test bypasses that guard by constructing
 * GFNFFGPUComputationalMethod directly and checks:
 *  (a) hasPeriodicEnergy() is false before and after a set-up without a cell
 *  (b) setMolecule() with a cell returns false and sets the error flag with a message
 *  (c) the same object type without a cell sets up and returns a finite energy
 * Needs a CUDA device and the plugin (only built with USE_CUDA).
 * B1 (Nov 24, 2026, Claude Generated, ROCm parity): the same source builds test_gfnff_rocm_pbc_refuse
 * (GFNFF_PBC_REFUSE_ROCM, HIP method class, same checks).
 */
#if defined(GFNFF_PBC_REFUSE_ROCM) || defined(USE_CUDA)
#include "src/core/curcuma_logger.h"
#ifdef GFNFF_PBC_REFUSE_ROCM
#include "src/core/energy_calculators/qm_methods/gfnff_hip_method.h"
using GFNFFGpuUnderTest = GFNFFHipComputationalMethod;
static const char* kTestName = "gfnff_rocm_pbc_refuse";
#else
#include "src/core/energy_calculators/qm_methods/gfnff_gpu_method.h"
using GFNFFGpuUnderTest = GFNFFGPUComputationalMethod;
static const char* kTestName = "gfnff_gpu_pbc_refuse";
#endif
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <string>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++g_failures;
}

int main()
{
    CurcumaLogger::set_verbosity(1);
    Mol plain = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
    Mol periodic = plain;
    periodic.m_has_pbc = true;
    periodic.m_unit_cell = Eigen::Matrix3d::Identity() * 200.0;

    GFNFFGpuUnderTest refused("gfnff", json::object());
    check(!refused.hasPeriodicEnergy(), "(a) hasPeriodicEnergy() false before set-up");
    const bool ok = refused.setMolecule(periodic);
    check(!ok, "(b) setMolecule() with a cell returns false");
    check(refused.hasError() && refused.getErrorMessage().find("periodic cell") != std::string::npos,
        "(b) error flag and message: " + refused.getErrorMessage());

    GFNFFGpuUnderTest normal("gfnff", json::object());
    check(normal.setMolecule(plain), "(c) setMolecule() without a cell succeeds");
    const double e = normal.calculateEnergy(false);
    check(std::isfinite(e) && e < 0.0, "(c) finite energy without a cell: " + std::to_string(e));
    check(!normal.hasPeriodicEnergy(), "(a) hasPeriodicEnergy() false after set-up");

    // S1/S2 (Nov 2026): the opt-in. gpu_periodic alone (coulomb_ewald default true) is the S2 hybrid: accepted, host
    // periodic EEQ + Ewald, hasPeriodicEnergy() true (S1 refused it with "S3/S4"); gpu_periodic with coulomb_ewald
    // false sets up, computes and keeps hasPeriodicEnergy() false.
    json opt = json::object();
    opt["gpu_periodic"] = true;
    GFNFFGpuUnderTest ewald("gfnff", opt);
    check(ewald.setMolecule(periodic), "(d) gpu_periodic with coulomb_ewald (default): setMolecule() accepted (S2 hybrid)");
    check(!ewald.hasError() && ewald.periodicHybridActive() && ewald.hasPeriodicEnergy(),
        "(d) hybrid active, hasPeriodicEnergy() true: " + ewald.getErrorMessage());
    const double e_hyb = ewald.calculateEnergy(true);
    check(std::isfinite(e_hyb) && e_hyb < 0.0 && ewald.hostEwaldEnergy() != 0.0,
        "(d) finite hybrid energy with a host Ewald part: " + std::to_string(e_hyb) + " (Ewald " + std::to_string(ewald.hostEwaldEnergy()) + ")");
    opt["coulomb_ewald"] = false;
    GFNFFGpuUnderTest released("gfnff", opt);
    check(released.setMolecule(periodic), "(e) gpu_periodic with coulomb_ewald false: setMolecule() accepted");
    const double e_cell = released.calculateEnergy(true);
    check(std::isfinite(e_cell) && e_cell < 0.0, "(e) finite energy with a cell: " + std::to_string(e_cell));
    check(!released.hasPeriodicEnergy(), "(e) hasPeriodicEnergy() stays false (S1: no Ewald on the device)");
    // the 200 A cube holds the whole dimer within L/2: the local terms equal the run without a cell to rounding
    GFNFFGpuUnderTest plain_ref("gfnff", json{ { "coulomb_ewald", false } });
    check(plain_ref.setMolecule(plain), "(e) reference without a cell");
    const double e_plain = plain_ref.calculateEnergy(true);
    check(std::abs(e_cell - e_plain) < 1e-9, "(e) 200 A cube == no cell to 1e-9 Eh (|dE| " + std::to_string(std::abs(e_cell - e_plain)) + ")");

    if (g_failures == 0)
        std::printf("%s: all checks passed\n", kTestName);
    else
        std::printf("%s: %d failure(s)\n", kTestName, g_failures);
    return g_failures == 0 ? 0 : 1;
}
#else
int main() { return 0; }
#endif
