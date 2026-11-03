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
 */
#ifdef USE_CUDA
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/gfnff_gpu_method.h"
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

    GFNFFGPUComputationalMethod refused("gfnff", json::object());
    check(!refused.hasPeriodicEnergy(), "(a) hasPeriodicEnergy() false before set-up");
    const bool ok = refused.setMolecule(periodic);
    check(!ok, "(b) setMolecule() with a cell returns false");
    check(refused.hasError() && refused.getErrorMessage().find("periodic cell") != std::string::npos,
        "(b) error flag and message: " + refused.getErrorMessage());

    GFNFFGPUComputationalMethod normal("gfnff", json::object());
    check(normal.setMolecule(plain), "(c) setMolecule() without a cell succeeds");
    const double e = normal.calculateEnergy(false);
    check(std::isfinite(e) && e < 0.0, "(c) finite energy without a cell: " + std::to_string(e));
    check(!normal.hasPeriodicEnergy(), "(a) hasPeriodicEnergy() false after set-up");

    std::printf(g_failures == 0 ? "gfnff_gpu_pbc_refuse: all checks passed\n" : "gfnff_gpu_pbc_refuse: %d failure(s)\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
#else
int main() { return 0; }
#endif
