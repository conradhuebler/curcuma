/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026) - per-triple GFN-FF hydrogen-bond energies on a GPU backend against the CPU.
 *
 * Usage: test_hbond_energies_gpu <cuda|rocm>
 *
 * Without a loadable plugin or a visible device it prints a SKIP line and exits 0 (nothing compared).
 * With a device, for each registry molecule (H2O_dimer, monosaccharide), through the same calls a CPU
 * caller uses (EnergyCalculator -> Interface()->gfnffInstance() -> setRecordHBondEnergies /
 * hydrogenBondEnergies):
 *  - nothing is returned before recording is switched on, nor after switching on but before an evaluation,
 *  - every entry carries its own atoms and case type, each triple appears once, and the set of GPU triples
 *    equals the CPU set with each triple's energy within 1e-8 Eh (matched by (donor, hydrogen, acceptor,
 *    case type; for case 1, where the H is bonded to neither atom, donor and acceptor are taken as one
 *    unordered pair), not by position),
 *  - host detection (-gfnff.gpu_hb_detect_on_device false): additionally one entry per triple of
 *    getLastHBonds(), in that order, and the GPU list equals the CPU list triple by triple. With device
 *    detection (the CUDA default since Oct 14, 2026) the device builds the same set in another order and the
 *    host list of the wrapper is not the list the device evaluates, so these two order checks do not apply,
 *  - the GPU sum equals the CPU HBond term within 1e-6 Eh,
 *  - recording leaves the GPU energy unchanged (1e-10 Eh), a gradient call records the same values,
 *  - switching recording off empties the result again.
 *
 * Tolerances: k_hbonds reads no EEQ or CN quantity (the triple charges are frozen in the list), so CPU
 * and GPU evaluate the same FP64 formula and differ by rounding (exp/sqrt implementations, FMA
 * contraction) and the 1e-12 in the GPU out-of-line ratio; 1e-8 Eh per triple is far above that and
 * matches the "8 decimals" CPU/ROCm agreement of whole GFN-FF energies (ff_methods/CLAUDE.md). The sum
 * uses the energy tolerance of test_gfnff_gpu.cpp (1e-6 Eh). 1e-10 Eh for on/off allows a different
 * order of the atomic block sums; GPU results are not run-to-run reproducible in general (CLAUDE.md #36).
 */

#include "src/core/energycalculator.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/gpu_fallback.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"
#include "test_molecule_registry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int g_failed = 0;

static void check(bool ok, const std::string& what)
{
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok)
        ++g_failed;
}

static double sum(const std::vector<GFNFFHBondEnergy>& v)
{
    double s = 0.0;
    for (const auto& t : v)
        s += t.energy;
    return s;
}

static bool sameAtoms(const std::vector<GFNFFHBondEnergy>& e, const std::vector<GFNFFHydrogenBond>& list)
{
    if (e.size() != list.size())
        return false;
    for (std::size_t i = 0; i < list.size(); ++i)
        if (e[i].donor != list[i].i || e[i].hydrogen != list[i].j || e[i].acceptor != list[i].k
            || e[i].case_type != list[i].case_type)
            return false;
    return true;
}

// Claude Generated (Nov 2026): energies keyed by triple identity; duplicate is set when a triple appears twice.
using HBKey = std::array<int, 4>;
static std::map<HBKey, double> keyed(const std::vector<GFNFFHBondEnergy>& v, bool& duplicate)
{
    std::map<HBKey, double> m;
    for (const auto& t : v)
        // Case 1 (an H bonded to neither atom, gfnff.h "Case 1 has no such coupling") treats A and B alike: the host
        // lists the higher atom index first, the device detection the lower one; one canonical order for both.
        if (!m.emplace(t.case_type == 1 ? HBKey { std::max(t.donor, t.acceptor), t.hydrogen, std::min(t.donor, t.acceptor), 1 }
                                        : HBKey { t.donor, t.hydrogen, t.acceptor, t.case_type }, t.energy).second)
            duplicate = true;
    return m;
}

static void compareMolecule(const std::string& name, const std::string& backend, bool device_detect)
{
    std::printf("%s on %s, GPU H-bond detection on the %s\n", name.c_str(), backend.c_str(), device_detect ? "device" : "host");
    const curcuma::Molecule molecule = TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
    const Mol mol = molecule.getMolInfo();

    // CPU reference
    json cpu_config = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", json::object() } };
    EnergyCalculator cpu("gfnff", cpu_config);
    cpu.setMolecule(mol);
    GFNFF* cpu_ff = cpu.Interface() ? cpu.Interface()->gfnffInstance() : nullptr;
    check(cpu_ff != nullptr, "CPU GFN-FF instance");
    if (!cpu_ff)
        return;
    cpu_ff->setRecordHBondEnergies(true);
    cpu.CalculateEnergy(false);
    const std::vector<GFNFFHBondEnergy> cpu_e = cpu_ff->hydrogenBondEnergies();
    const std::vector<GFNFFHydrogenBond>& cpu_list = cpu_ff->getLastHBonds();
    const double cpu_term = cpu.getEnergyDecomposition().value("HBond", 0.0);

    // GPU, through the same calls
    const auto fallbacks_before = curcuma::gpuFallbackCounts();
    json gpu_config = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", { { "gpu_hb_detect_on_device", device_detect } } }, { "gpu", backend } };
    EnergyCalculator gpu("gfnff", gpu_config);
    gpu.setMolecule(mol);
    check(curcuma::gpuFallbackCounts() == fallbacks_before, "the " + backend + " method runs (no GPU fallback)");
    GFNFF* gpu_ff = gpu.Interface() ? gpu.Interface()->gfnffInstance() : nullptr;
    check(gpu_ff != nullptr, "GPU GFN-FF instance");
    if (!gpu_ff)
        return;

    const double e_off = gpu.CalculateEnergy(false);
    check(gpu_ff->hydrogenBondEnergies().empty(), "nothing is returned unless recording is on");
    gpu_ff->setRecordHBondEnergies(true);
    check(gpu_ff->hydrogenBondEnergies().empty(), "switched on, nothing before the next evaluation");
    const double e_on = gpu.CalculateEnergy(false);
    const std::vector<GFNFFHBondEnergy> gpu_e = gpu_ff->hydrogenBondEnergies();
    const std::vector<GFNFFHydrogenBond>& gpu_list = gpu_ff->getLastHBonds();

    bool dup_gpu = false, dup_cpu = false;
    const std::map<HBKey, double> gpu_m = keyed(gpu_e, dup_gpu);
    const std::map<HBKey, double> cpu_m = keyed(cpu_e, dup_cpu);
    check(!gpu_e.empty() && !dup_gpu && !dup_cpu, "GPU entries: " + std::to_string(gpu_e.size()) + ", each triple once");
    bool same_set = gpu_m.size() == cpu_m.size();
    for (const auto& kv : gpu_m)
        if (!cpu_m.count(kv.first))
            same_set = false;
    check(same_set, "GPU and CPU have the same set of triples (" + std::to_string(gpu_m.size()) + " against "
                        + std::to_string(cpu_m.size()) + ")");
    if (!same_set) {
        // failure aid: the first triples only one side has
        int shown = 0;
        for (const auto& kv : gpu_m)
            if (!cpu_m.count(kv.first) && shown++ < 6)
                std::printf("      only GPU: donor %d hydrogen %d acceptor %d case %d  E = %.3e Eh\n", kv.first[0],
                    kv.first[1], kv.first[2], kv.first[3], kv.second);
        shown = 0;
        for (const auto& kv : cpu_m)
            if (!gpu_m.count(kv.first) && shown++ < 6)
                std::printf("      only CPU: donor %d hydrogen %d acceptor %d case %d  E = %.3e Eh\n", kv.first[0],
                    kv.first[1], kv.first[2], kv.first[3], kv.second);
    }
    if (!device_detect) {
        // host detection: the entries follow getLastHBonds() exactly, as on the CPU
        check(!gpu_e.empty() && gpu_e.size() == gpu_list.size(),
            "one entry per triple of getLastHBonds() (" + std::to_string(gpu_e.size()) + ")");
        check(sameAtoms(gpu_e, gpu_list), "entries in list order with the list's atoms and case type");
        check(gpu_e.size() == cpu_e.size() && sameAtoms(gpu_e, cpu_list), "GPU and CPU have the same triples, same order");
    }

    double max_diff = 0.0;
    if (same_set)
        for (const auto& kv : gpu_m)
            max_diff = std::max(max_diff, std::fabs(kv.second - cpu_m.at(kv.first)));
    char buf[160];
    std::snprintf(buf, sizeof(buf), "per-triple energy GPU vs CPU, max |dE| = %.3e Eh (tol 1e-8)", max_diff);
    check(same_set && max_diff < 1e-8, buf);
    std::snprintf(buf, sizeof(buf), "GPU sum %.10f vs CPU HBond term %.10f Eh (tol 1e-6)", sum(gpu_e), cpu_term);
    check(std::fabs(sum(gpu_e) - cpu_term) < 1e-6, buf);
    std::snprintf(buf, sizeof(buf), "recording leaves the GPU energy unchanged (%.3e Eh)", e_on - e_off);
    check(std::fabs(e_on - e_off) < 1e-10, buf);

    gpu.CalculateEnergy(true);
    const std::vector<GFNFFHBondEnergy> gpu_e_grad = gpu_ff->hydrogenBondEnergies();
    double max_grad_diff = gpu_e_grad.size() == gpu_e.size() ? 0.0 : 1.0;
    for (std::size_t i = 0; i < gpu_e_grad.size() && i < gpu_e.size(); ++i)
        max_grad_diff = std::max(max_grad_diff, std::fabs(gpu_e_grad[i].energy - gpu_e[i].energy));
    check(max_grad_diff < 1e-12, "a gradient call records the same values");

    gpu_ff->setRecordHBondEnergies(false);
    check(gpu_ff->hydrogenBondEnergies().empty(), "switched off, nothing is returned");
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s cuda|rocm\n", argv[0]);
        return 2;
    }
    const std::string backend = argv[1];
    std::printf("hbond energies gpu (%s)\n", backend.c_str());
    if (!gpu_plugin::available(backend)) {
        std::printf("SKIP: the %s plugin cannot be loaded here; nothing compared\n", backend.c_str());
        return 0;
    }
    if (gpu_plugin::deviceCount(backend) <= 0) {
        std::printf("SKIP: the %s plugin sees no device; nothing compared\n", backend.c_str());
        return 0;
    }

    for (const bool device_detect : { true, false }) {
        compareMolecule("H2O_dimer", backend, device_detect);
        compareMolecule("monosaccharide", backend, device_detect);
    }

    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
