/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * GPU guard for periodic cells (WP-PERIODIC-NONBONDED A3b) - Claude Generated (Nov 2026).
 *
 * The GPU GFN-FF paths have no minimum image, so MethodFactory::periodicCellNeedsCpu() must send a
 * GPU GFN-FF method with a cell to the CPU (EnergyCalculator::setMolecule(), deltaml base). This test
 * runs WITHOUT a GPU: a stand-in for the plugin wrapper (ComputationalMethod with a GFNFF object that
 * is not GFNFFComputationalMethod) takes the place of GFNFFGpuMethodImpl<Backend>.
 *
 * (a) stand-in, no cell: false, no fallback counted
 * (b) stand-in, cell: true, one fallback counted under the "periodic cell" category
 * (c) CPU GFNFFComputationalMethod with a cell: false (CPU path is the periodic one)
 * (d) method without a GFNFF object (EHT stand-in) with a cell: false
 * (e) -gpu_strict: the stand-in with a cell ends a forked child with exit code 3
 * (f) EnergyCalculator gfnff -gpu cuda, plugin absent here, with a cell: runs on the CPU, no
 *     "periodic cell" fallback (only the plugin-absent one), energy equal to the run without -gpu
 *
 * Not covered: a real CUDA/ROCm wrapper (needs the plugin and a device).
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/method_factory.h"
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"
#include "src/core/energycalculator.h"
#include "src/core/gpu_fallback.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

/// Stand-in for a GPU method: answers gfnffInstance() like the plugin wrapper, or nullptr (other methods).
class FakeMethod : public ComputationalMethod {
public:
    explicit FakeMethod(bool with_gfnff)
    {
        if (with_gfnff) m_gfnff = std::make_unique<GFNFF>(json::object());
    }
    GFNFF* gfnffInstance() const override { return m_gfnff.get(); }
    bool setMolecule(const Mol&) override { return true; }
    bool updateGeometry(const Matrix&) override { return true; }
    double calculateEnergy(bool) override { return 0.0; }
    Matrix getGradient() const override { return Matrix(); }
    Vector getCharges() const override { return Vector(); }
    Vector getBondOrders() const override { return Vector(); }
    Position getDipole() const override { return Position::Zero(); }
    bool hasGradient() const override { return false; }
    std::string getMethodName() const override { return "fake"; }
    bool isThreadSafe() const override { return false; }
    void setThreadCount(int) override {}
    void setParameters(const json&) override {}
    json getParameters() const override { return json::object(); }
    bool hasError() const override { return false; }
    json getEnergyDecomposition() const override { return json::object(); }

private:
    std::unique_ptr<GFNFF> m_gfnff;
};

static int periodicCount()
{
    int n = 0;
    for (const auto& c : curcuma::gpuFallbackCounts())
        if (c.first.find("periodic cell") != std::string::npos) n += c.second;
    return n;
}

/// Water dimer from the structure library (registry "H2O_dimer", Angstrom) with an optional 200 A cube.
static Mol waterDimer(bool cell)
{
    Mol m = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
    if (cell) {
        m.m_has_pbc = true;
        m.m_unit_cell = Eigen::Matrix3d::Identity() * 200.0;
    }
    return m;
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const Mol plain = waterDimer(false), periodic = waterDimer(true);

    FakeMethod gpu_like(true), other(false);
    check(MethodFactory::isGpuGfnffMethod(gpu_like), "(a) stand-in counts as GPU GFN-FF");
    check(!MethodFactory::periodicCellNeedsCpu(plain, gpu_like), "(a) no cell: no fallback");
    check(periodicCount() == 0, "(a) no cell: nothing counted");

    check(MethodFactory::periodicCellNeedsCpu(periodic, gpu_like), "(b) cell: fallback requested");
    check(periodicCount() == 1, "(b) cell: one fallback counted");

    GFNFFComputationalMethod cpu("gfnff", json::object());
    check(!MethodFactory::isGpuGfnffMethod(cpu), "(c) CPU GFN-FF is not a GPU method");
    check(!MethodFactory::periodicCellNeedsCpu(periodic, cpu), "(c) CPU GFN-FF with cell: no fallback");
    check(!MethodFactory::periodicCellNeedsCpu(periodic, other), "(d) method without GFNFF: no fallback");
    check(periodicCount() == 1, "(c,d) still one fallback");

    // (e) strict mode in a forked child (reportGpuFallback exits the process)
    std::fflush(stdout);
    const pid_t pid = fork();
    if (pid == 0) {
        curcuma::setGpuStrict(true);
        MethodFactory::periodicCellNeedsCpu(periodic, gpu_like);
        _exit(0);   // not reached under -gpu_strict
    }
    int status = 0;
    waitpid(pid, &status, 0);
    check(WIFEXITED(status) && WEXITSTATUS(status) == curcuma::kGpuStrictExitCode,
        "(e) -gpu_strict: child exits with code 3 (status " + std::to_string(status) + ")");

    // (f) EnergyCalculator: plugin absent -> CPU method from the factory, the guard stays silent
    json ctrl = { { "method", "gfnff" }, { "threads", 1 }, { "verbosity", 0 } };
    EnergyCalculator ref("gfnff", ctrl);
    ref.setMolecule(periodic);
    const double e_ref = ref.CalculateEnergy(false);
    json ctrl_gpu = ctrl;
    ctrl_gpu["gpu"] = "cuda";
    EnergyCalculator gpu_req("gfnff", ctrl_gpu);
    gpu_req.setMolecule(periodic);
    const double e_gpu = gpu_req.CalculateEnergy(false);
    check(!ref.Error() && !gpu_req.Error(), "(f) both calculators ran");
    check(std::isfinite(e_ref) && e_ref == e_gpu, "(f) -gpu cuda without plugin, cell: CPU energy");
    check(periodicCount() == 1, "(f) no periodic-cell fallback when the factory already chose the CPU");

    if (g_failures == 0)
        std::printf("pbc_gpu_guard: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
