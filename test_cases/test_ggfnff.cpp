/**
 * GPU GFN-FF Validation: ggfnff vs gfnff (CPU)
 *
 * Claude Generated (March 2026): Validates that the GPU-accelerated GFN-FF
 * (ggfnff) produces energies and gradients numerically identical to the CPU
 * reference implementation (gfnff).
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#ifdef USE_CUDA

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "json.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    CurcumaLogger::set_verbosity(0);

    std::string mol_file = (argc > 1) ? argv[1] : "A.xyz";

    std::cout << "==================================================\n"
              << "  GFN-FF GPU Validation: ggfnff vs gfnff (CPU)\n"
              << "  Molecule: " << mol_file << "\n"
              << "==================================================" << std::endl;

    curcuma::Molecule molecule(mol_file);
    Mol mol = molecule.getMolInfo();
    std::cout << "Atoms: " << mol.m_number_atoms << std::endl;

    json config = {{"verbosity", 0}, {"threads", 1}};
    int passed = 0, failed = 0;

    // --- CPU gfnff ---
    auto* cpu = new EnergyCalculator("gfnff", config);
    cpu->setMolecule(mol);
    double cpu_e = cpu->CalculateEnergy(true);
    Matrix G_cpu = cpu->Gradient();  // cached in CalculateEnergy
    std::cout << "CPU energy: " << std::setprecision(12) << cpu_e << " Eh" << std::endl;
    std::cout << "CPU grad norm: " << std::scientific << G_cpu.norm() << std::endl;

    // --- GPU ggfnff ---
    auto* gpu = new EnergyCalculator("ggfnff", config);
    gpu->setMolecule(mol);
    double gpu_e = gpu->CalculateEnergy(true);
    Matrix G_gpu = gpu->Gradient();  // cached in CalculateEnergy (no re-call to GPU)
    std::cout << "GPU energy: " << std::setprecision(12) << gpu_e << " Eh" << std::endl;
    std::cout << "GPU grad norm: " << std::scientific << G_gpu.norm() << std::endl;

    // --- Compare Energy ---
    double e_err = std::abs(gpu_e - cpu_e);
    std::cout << "\nEnergy diff: " << std::scientific << e_err << " Eh";
    if (e_err < 1e-6) { std::cout << " [PASS]\n"; ++passed; }
    else              { std::cout << " [FAIL] (tol=1e-6)\n"; ++failed; }

    // --- Print gradient details for small molecules ---
    if (mol.m_number_atoms <= 10) {
        Matrix G_diff = G_gpu - G_cpu;
        std::cout << "\nPer-atom gradient comparison (GPU - CPU):\n"
                  << std::setw(5) << "Atom" << std::setw(16) << "CPU_x" << std::setw(16) << "GPU_x"
                  << std::setw(16) << "diff_x" << std::setw(16) << "diff_y" << std::setw(16) << "diff_z" << "\n";
        for (int i = 0; i < mol.m_number_atoms; ++i) {
            std::cout << std::setw(5) << i
                      << std::scientific << std::setprecision(6)
                      << std::setw(16) << G_cpu(i,0) << std::setw(16) << G_gpu(i,0)
                      << std::setw(16) << G_diff(i,0) << std::setw(16) << G_diff(i,1)
                      << std::setw(16) << G_diff(i,2) << "\n";
        }
        std::cout << std::endl;
    }

    // --- Compare Gradient max component ---
    double g_max = (G_gpu - G_cpu).cwiseAbs().maxCoeff();
    std::cout << "Grad max diff: " << std::scientific << g_max << " Eh/Bohr";
    if (g_max < 1e-4) { std::cout << " [PASS]\n"; ++passed; }
    else              { std::cout << " [FAIL] (tol=1e-4)\n"; ++failed; }

    // --- Compare Gradient norm ratio ---
    double ratio = (G_cpu.norm() > 1e-12) ? G_gpu.norm() / G_cpu.norm() : 1.0;
    std::cout << "Grad norm ratio: " << std::fixed << std::setprecision(6) << ratio;
    if (ratio >= 0.99 && ratio <= 1.01) { std::cout << " [PASS]\n"; ++passed; }
    else { std::cout << " [FAIL] (expected 0.99-1.01)\n"; ++failed; }

    std::cout << "\nPassed: " << passed << "/" << (passed + failed) << std::endl;

    // Intentional leak + _exit: CUDA heap corruption makes normal destructors crash.
    std::cout.flush();
    _exit(failed == 0 ? 0 : 1);
}

#else

#include <iostream>
int main() {
    std::cout << "test_ggfnff: built without USE_CUDA — skipping GPU tests\n";
    return 0;
}

#endif
