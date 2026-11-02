/*
 * Full-precision energy dumper for the embedded-basis tests (qm_basis_*).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Runs one native QM method through EnergyCalculator (the same dispatch as `curcuma -sp`) and prints the total
 * energy with 17 significant digits, so that check_qm_basis_embed.py can compare a run with the basis read from a
 * file against a run with the embedded copy bit for bit (the CLI prints only 8 decimals).
 *
 *   dump_qm_energy <xyz> <method> [--qm KEY VALUE]... [--verbosity N] [--gradient] [--threads N]
 *
 * --gradient (Nov 2026, used by test_cases/qm_fg): also prints the analytic gradient, one line per atom,
 * "grad <atom> gx gy gz" in Eh/Bohr (the EnergyCalculator contract is Eh/Angstrom; converted here with units.h).
 * --charge Q / --spin S (Nov 2026, stage 5b, test_cases/qm_uhf): molecular charge and unpaired electrons; an
 * unrestricted run also prints "s2 <S^2>".
 * --gpu NAME (Nov 2026, stage 6a, test_cases/qm_uhf GPU tests): top-level "gpu" of the config (cuda, ...); every
 * counted GPU fallback of the run is printed at the end as "fallback <category> <count>".
 *
 * VALUE is parsed as JSON when possible (true, 1e-10, ...), else taken as a string.
 *
 * Claude Generated (Nov 2026). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/gpu_fallback.h"
#include "src/core/molecule.h"
#include "src/core/units.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: dump_qm_energy <xyz> <method> [--qm KEY VALUE]... [--verbosity N]\n";
        return 1;
    }
    const std::string xyz = argv[1], method = argv[2];
    int verbosity = 0, threads = 1;
    bool gradient = false;
    int charge = 0, spin = 0;
    bool set_charge = false, set_spin = false;
    std::string gpu;
    json qm = { { "scf_threshold", 1e-10 }, { "scf_max_iterations", 200 } };
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--qm" && i + 2 < argc) {
            const std::string key = argv[++i], val = argv[++i];
            const json parsed = json::parse(val, nullptr, false);
            qm[key] = parsed.is_discarded() ? json(val) : parsed;
        } else if (a == "--verbosity" && i + 1 < argc) {
            verbosity = std::stoi(argv[++i]);
        } else if (a == "--threads" && i + 1 < argc) {
            threads = std::stoi(argv[++i]);
        } else if (a == "--charge" && i + 1 < argc) {
            charge = std::stoi(argv[++i]);
            set_charge = true;
        } else if (a == "--spin" && i + 1 < argc) {
            spin = std::stoi(argv[++i]);
            set_spin = true;
        } else if (a == "--gpu" && i + 1 < argc) {
            gpu = argv[++i];
        } else if (a == "--gradient") {
            gradient = true;
        }
    }
    CurcumaLogger::set_verbosity(verbosity);
    json config = { { "verbosity", verbosity }, { "threads", threads }, { "qm", qm } };
    if (!gpu.empty()) config["gpu"] = gpu;
    curcuma::Molecule mol(xyz);
    if (set_charge) mol.setCharge(charge);  // only when given: the xyz comment may carry them
    if (set_spin) mol.setSpin(spin);
    EnergyCalculator calc(method, config);
    calc.setMolecule(mol.getMolInfo());
    const double e = calc.CalculateEnergy(gradient);
    if (calc.Error() || !std::isfinite(e)) {
        std::cerr << "dump_qm_energy: calculation failed\n";
        return 2;
    }
    std::fflush(stdout);
    std::printf("\nenergy %.17g\n", e);  // own line: the logger may leave a colour code on the current one
    const json dec = calc.getEnergyDecomposition();
    if (dec.contains("s2") && dec["s2"].is_number()) std::printf("s2 %.17g\n", dec["s2"].get<double>());
    if (gradient) {
        const Matrix& g = calc.Gradient();  // Eh/Angstrom
        const double f = 1.0 / CurcumaUnit::Length::angstrom_to_bohr(1.0);  // -> Eh/Bohr
        for (int i = 0; i < g.rows(); ++i)
            std::printf("grad %d %.17g %.17g %.17g\n", i, g(i, 0) * f, g(i, 1) * f, g(i, 2) * f);
    }
    for (const auto& c : curcuma::gpuFallbackCounts()) std::printf("fallback %s %d\n", c.first.c_str(), c.second);
    return 0;
}
