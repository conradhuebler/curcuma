/*
 * test_gfnff_grad_traj.cpp — per-frame gradient comparison on a fixed trajectory
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (June 2026): the clean force-backend diagnostic.
 *
 * Motivation: MD "Exchange with heat bath" is NOT a valid cross-backend force
 * test — it conflates force differences with integrator/thermostat/velocity-init
 * differences and chaotic amplification (see docs/GPU_GFNNF_DISCREPANCIES.md). The
 * correct metric is to evaluate the gradient of EACH backend at EACH frame of ONE
 * FIXED trajectory and report max||g_a(x_k) - g_b(x_k)|| per frame. No integrator,
 * no chaos — a pure force-backend comparison.
 *
 * Usage:
 *   test_gfnff_grad_traj <trajectory.xyz> [method=gfnff] [gpu_backend=cuda]
 * Compares backend A = "<method> (CPU)" vs backend B = "<method> -gpu <backend>".
 * When CUDA is not compiled in, the gpu config falls back to CPU and the diff is 0.
 */
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "json.hpp"
using json = nlohmann::json;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <trajectory.xyz> [method=gfnff] [gpu_backend=cuda]\n";
        return 1;
    }
    const std::string traj = argv[1];
    const std::string method = (argc > 2) ? argv[2] : "gfnff";
    const std::string gpu = (argc > 3) ? argv[3] : "cuda";

    // Backend A: CPU. Backend B: same method on the GPU backend.
    json cfg_a = { { "verbosity", 0 }, { "threads", 1 }, { method, json::object() } };
    json cfg_b = cfg_a;
    cfg_b["gpu"] = gpu;

    FileIterator it;
    it.setFile(traj);

    int frame = 0;
    double global_max = 0.0;
    int global_max_frame = -1;
    double sum_max = 0.0;

    std::cout << "Per-frame gradient comparison: " << method << " (CPU) vs " << method
              << " -gpu " << gpu << "\n";
    std::cout << "trajectory: " << traj << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << std::setw(6) << "frame" << std::setw(20) << "max|g_a-g_b|"
              << std::setw(18) << "||dg||/||g_a||" << std::setw(16) << "dE [Eh]" << "\n";

    while (!it.AtEnd()) {
        Molecule mol = it.Next();
        if (mol.AtomCount() == 0)
            break;

        EnergyCalculator a(method, cfg_a);
        a.setMolecule(mol.getMolInfo());
        const double ea = a.CalculateEnergy(true);
        Matrix ga = a.Gradient();

        EnergyCalculator b(method, cfg_b);
        b.setMolecule(mol.getMolInfo());
        const double eb = b.CalculateEnergy(true);
        Matrix gb = b.Gradient();

        const double maxd = (ga - gb).cwiseAbs().maxCoeff();
        const double rel = (ga - gb).norm() / (ga.norm() + 1e-30);
        if (maxd > global_max) {
            global_max = maxd;
            global_max_frame = frame;
        }
        sum_max += maxd;

        std::cout << std::setw(6) << frame << std::setw(20) << maxd << std::setw(18) << rel
                  << std::setw(16) << (eb - ea) << "\n";
        frame++;
    }

    std::cout << "------------------------------------------------------------\n";
    if (frame > 0) {
        std::cout << "frames = " << frame << "   mean max|g_a-g_b| = " << (sum_max / frame)
                  << "   global max = " << global_max << " (frame " << global_max_frame << ")\n";
    } else {
        std::cout << "No frames read from " << traj << "\n";
        return 1;
    }
    return 0;
}
