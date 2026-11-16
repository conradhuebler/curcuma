/*
 * Optimizer default layers (docs/debt_2026-11/PARAMS.md H5). Claude Generated (Nov 2026).
 *
 * 1. Every key that both the registry module "opt" and OptimizerInterfaceJson define has the same
 *    default. The registry fills the config of -opt, OptimizerInterfaceJson the config of the direct
 *    callers (ConfSearch, ConfGen, SimpleMD, Docking, PolymerBuild, RMSDTraj); a difference gave the
 *    two paths different optimisers. Keys that mean something else in the two layers are listed
 *    with the reason.
 * 2. A direct call with an EMPTY config optimises ethane (GFN-FF) to convergence. With
 *    lbfgs_max_line_search 2 (the former LBFGSppOptimizer default) the line search fails at the first
 *    step, no step is taken and the optimisation ends not converged; measured: 2 and 3 fail, 4 and
 *    more converge, 5 and more give the same steps and energy as 20.
 *
 * Usage: optimizer_default_layers_test ethane.xyz
 */
#include "generated/parameter_registry.h"
#include "src/capabilities/optimizer_factory.h"
#include "src/capabilities/optimizer_interface.h"
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/parameter_registry.h"
#include "src/tools/formats.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <set>

static int g_failed = 0;
#define CHECK(cond, text)                                                  \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::cerr << "FAIL: " << text << " (" #cond ")\n";             \
            ++g_failed;                                                    \
        }                                                                  \
    } while (0)

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: optimizer_default_layers_test ethane.xyz\n";
        return 2;
    }
    initialize_generated_registry();

    // 1. the two layers
    const nlohmann::json registry = ParameterRegistry::getInstance().getDefaultJson("opt");
    // method: the optimiser in OptimizerInterfaceJson, the energy method in the registry.
    // verbose, print_output, optimizer_type: only in OptimizerInterfaceJson / different meaning.
    const std::set<std::string> different_meaning = { "method", "verbose", "verbosity", "threads" };
    int compared = 0;
    for (const auto& [key, value] : Optimization::OptimizerInterfaceJson.items()) {
        if (different_meaning.count(key) || !registry.contains(key))
            continue;
        ++compared;
        const bool same = (value.is_number() && registry[key].is_number())
            ? std::abs(value.get<double>() - registry[key].get<double>()) <= 1e-12 * (1.0 + std::abs(value.get<double>()))
            : value == registry[key];
        CHECK(same, "default of '" + key + "' differs: OptimizerInterfaceJson " + value.dump() + " registry " + registry[key].dump());
    }
    CHECK(compared >= 15, "compared at least 15 shared keys (got " + std::to_string(compared) + ")");

    // 2. direct call with an empty config
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "curcuma_default_layers";
    fs::create_directories(dir);
    const fs::path file = dir / "ethane.xyz";
    fs::copy_file(argv[1], file, fs::copy_options::overwrite_existing);
    Molecule mol = Files::LoadFile(file.string());
    EnergyCalculator calc("gfnff", nlohmann::json{ { "geometry_file", file.string() } });
    nlohmann::json config = { { "write_trajectory", false } };
    auto result = Optimization::OptimizationDispatcher::optimizeStructure(&mol, Optimization::OptimizerType::LBFGSPP, &calc, config);
    CHECK(result.success, "direct call with default LBFGSpp config converges for ethane");
    CHECK(result.iterations_performed >= 3, "the optimiser took steps (iterations " + std::to_string(result.iterations_performed) + ")");
    fs::remove_all(dir);

    std::cout << (g_failed ? "FAILED" : "all checks passed") << " (" << compared << " shared keys compared)\n";
    return g_failed ? 1 : 0;
}
