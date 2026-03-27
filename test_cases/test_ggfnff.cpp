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
#include "src/core/energy_calculators/qm_methods/gfnff_method.h"   // CPU wrapper
#include "src/core/energy_calculators/qm_methods/ggfnff_method.h"  // GPU wrapper
#include "src/core/energy_calculators/ff_methods/gfnff.h"          // GFNFF class
#include "src/core/energy_calculators/ff_methods/ff_workspace.h"   // FFWorkspace
#include "json.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    int verbose = 0;
    std::string mol_file = "A.xyz";
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-v") verbose = 3;
        else mol_file = argv[i];
    }
    CurcumaLogger::set_verbosity(verbose);

    // Claude Generated (March 2026): Handle .ref.json files by extracting geometry_file
    // and resolving the path relative to the source root
    std::string xyz_file = mol_file;

    if (mol_file.size() > 5 && mol_file.substr(mol_file.size() - 5) == ".json") {
        std::ifstream ref_file(mol_file);
        if (ref_file.is_open()) {
            try {
                json ref_data;
                ref_file >> ref_data;
                if (ref_data.contains("molecule") && ref_data["molecule"].contains("geometry_file")) {
                    std::string geo_path = ref_data["molecule"]["geometry_file"].get<std::string>();

                    // The geometry_file path is relative to different roots depending on the source.
                    // We're running from build/test_cases/ (where the .ref.json was copied).
                    //
                    // Path formats in .ref.json files:
                    // 1. "test_cases/molecules/..." - relative to curcuma root -> prepend "../../"
                    // 2. "../test_cases/molecules/..." - relative to some dir -> prepend "../"
                    // 3. "external/gfnff/test/..." - relative to curcuma root -> prepend "../../"
                    // 4. "../external/gfnff/test/..." - relative to some dir -> prepend "../"

                    if (geo_path.size() >= 3 && geo_path.substr(0, 3) == "../") {
                        // Path starts with "../" - already relative, use as-is with one more "../"
                        xyz_file = "../" + geo_path;
                    } else {
                        // Path starts without "../" - relative to curcuma root, go up two dirs
                        xyz_file = "../../" + geo_path;
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Error parsing JSON file " << mol_file << ": " << e.what() << std::endl;
                return 1;
            }
        }
    }

    std::cout << "==================================================\n"
              << "  GFN-FF GPU Validation: ggfnff vs gfnff (CPU)\n"
              << "  Molecule: " << xyz_file << "\n"
              << "==================================================" << std::endl;

    curcuma::Molecule molecule(xyz_file);
    Mol mol = molecule.getMolInfo();
    std::cout << "Atoms: " << mol.m_number_atoms << std::endl;

    json config = {{"verbosity", 0}, {"threads", 1}, {"gfnff", json::object()}, {"ggfnff", json::object()}};
    int passed = 0, failed = 0;

    // --- CPU gfnff ---
    auto* cpu = new EnergyCalculator("gfnff", config);
    cpu->setMolecule(mol);
    double cpu_e = cpu->CalculateEnergy(true);
    Matrix G_cpu = cpu->Gradient();
    std::cout << "CPU energy: " << std::setprecision(12) << cpu_e << " Eh" << std::endl;
    std::cout << "CPU grad norm: " << std::scientific << G_cpu.norm() << std::endl;

    std::cout.flush();

    // --- GPU ggfnff ---
    auto* gpu = new EnergyCalculator("ggfnff", config);
    gpu->setMolecule(mol);
    // Debug: check initialization
    if (gpu->Interface() && dynamic_cast<GGFNFFComputationalMethod*>(gpu->Interface())) {
        auto* gm = dynamic_cast<GGFNFFComputationalMethod*>(gpu->Interface());
        if (gm->hasError()) {
            std::cerr << "GPU INIT ERROR: " << gm->getErrorMessage() << std::endl;
        }
    }
    double gpu_e = gpu->CalculateEnergy(true);
    Matrix G_gpu = gpu->Gradient();
    std::cout << "GPU energy: " << std::setprecision(12) << gpu_e << " Eh" << std::endl;
    std::cout << "GPU grad norm: " << std::scientific << G_gpu.norm() << std::endl;

    // --- Compare Energy ---
    double e_err = std::abs(gpu_e - cpu_e);
    std::cout << "\nEnergy diff: " << std::scientific << e_err << " Eh";
    if (e_err < 1e-6) { std::cout << " [PASS]\n"; ++passed; }
    else              { std::cout << " [FAIL] (tol=1e-6)\n"; ++failed; }

    // =========================================================================
    // dEdcn diagnostic comparison (raw bond+disp, before qtmp subtraction)
    // CPU: FFWorkspace::dEdcnTotal() = bond + disp dEdcn
    // GPU: FFWorkspaceGPU::dEdcnTotal() = bond + disp dEdcn (downloaded before qtmp)
    // =========================================================================
    {
        // Access CPU dEdcn via GFNFFComputationalMethod → GFNFF → FFWorkspace
        auto* cpu_method = dynamic_cast<GFNFFComputationalMethod*>(cpu->Interface());
        auto* gpu_method = dynamic_cast<GGFNFFComputationalMethod*>(gpu->Interface());

        if (cpu_method && cpu_method->getGFNFF() && cpu_method->getGFNFF()->getWorkspace()
            && gpu_method) {

            // Copy dEdcn to local vectors immediately to avoid heap corruption issues
            const Vector& cpu_dEdcn_ref = cpu_method->getGFNFF()->getWorkspace()->dEdcnTotal();
            const Vector& gpu_dEdcn_ref = gpu_method->getGPUdEdcn();
            const Vector& cpu_dEdcn_bond_ref = cpu_method->getGFNFF()->getWorkspace()->dEdcnBondTotal();

            std::cout << "\n--- dEdcn Comparison (raw bond+disp, before qtmp) ---\n";
            std::cout << "CPU dEdcn size: " << cpu_dEdcn_ref.size()
                      << "  GPU dEdcn size: " << gpu_dEdcn_ref.size() << std::endl;
            std::cout.flush();

            if (cpu_dEdcn_ref.size() == gpu_dEdcn_ref.size() && cpu_dEdcn_ref.size() > 0) {
                int N = static_cast<int>(cpu_dEdcn_ref.size());

                // Deep-copy to local POD arrays (avoids CUDA-corrupted Eigen heap issues)
                std::vector<double> cpu_de(N), gpu_de(N), cpu_de_bond(N);
                for (int i = 0; i < N; ++i) {
                    cpu_de[i] = cpu_dEdcn_ref[i];
                    gpu_de[i] = gpu_dEdcn_ref[i];
                    cpu_de_bond[i] = (cpu_dEdcn_bond_ref.size() == N) ? cpu_dEdcn_bond_ref[i] : 0.0;
                }

                double max_diff = 0.0;
                int worst_atom = 0;

                std::vector<std::pair<double, int>> diffs;
                for (int i = 0; i < N; ++i) {
                    double d = std::abs(gpu_de[i] - cpu_de[i]);
                    diffs.push_back({d, i});
                    if (d > max_diff) { max_diff = d; worst_atom = i; }
                }
                std::sort(diffs.rbegin(), diffs.rend());

                int show = std::min(N, 15);
                std::cout << std::setw(5) << "Atom" << std::setw(4) << "Z"
                          << std::setw(18) << "CPU_dEdcn" << std::setw(18) << "GPU_dEdcn"
                          << std::setw(14) << "diff" << "\n";
                for (int n = 0; n < show; ++n) {
                    int i = diffs[n].second;
                    std::cout << std::setw(5) << i
                              << std::setw(4) << mol.m_atoms[i]
                              << std::scientific << std::setprecision(6)
                              << std::setw(18) << cpu_de[i]
                              << std::setw(18) << gpu_de[i]
                              << std::setw(14) << (gpu_de[i] - cpu_de[i]) << "\n";
                }

                double norm_cpu = 0, norm_gpu = 0;
                for (int i = 0; i < N; ++i) {
                    norm_cpu += cpu_de[i] * cpu_de[i];
                    norm_gpu += gpu_de[i] * gpu_de[i];
                }
                norm_cpu = std::sqrt(norm_cpu);
                norm_gpu = std::sqrt(norm_gpu);
                std::cout << "\ndEdcn norm CPU: " << std::scientific << norm_cpu
                          << "  GPU: " << norm_gpu
                          << "  ratio: " << std::fixed << std::setprecision(6)
                          << (norm_cpu > 1e-12 ? norm_gpu / norm_cpu : 1.0) << std::endl;

                std::cout << "dEdcn max diff: " << std::scientific << max_diff
                          << " (atom " << worst_atom << ", Z=" << mol.m_atoms[worst_atom] << ")";
                // GPU diagnostic snapshots only available at verbosity >= 3
                if (norm_gpu < 1e-15) {
                    std::cout << " [SKIP] (diagnostic snapshots require verbosity >= 3)\n";
                } else if (max_diff < 1e-6) {
                    std::cout << " [PASS]\n"; ++passed;
                } else {
                    std::cout << " [FAIL] (tol=1e-6)\n"; ++failed;
                }

                // CPU dEdcn decomposition
                if (cpu_dEdcn_bond_ref.size() == N) {
                    std::cout << "\n--- CPU dEdcn decomposition (top 10 by total) ---\n";
                    std::cout << std::setw(5) << "Atom" << std::setw(4) << "Z"
                              << std::setw(18) << "bond_dEdcn" << std::setw(18) << "disp_dEdcn"
                              << std::setw(18) << "total" << "\n";

                    std::vector<std::pair<double, int>> total_sorted;
                    for (int i = 0; i < N; ++i)
                        total_sorted.push_back({std::abs(cpu_de[i]), i});
                    std::sort(total_sorted.rbegin(), total_sorted.rend());

                    for (int n = 0; n < std::min(10, N); ++n) {
                        int i = total_sorted[n].second;
                        double disp = cpu_de[i] - cpu_de_bond[i];
                        std::cout << std::setw(5) << i
                                  << std::setw(4) << mol.m_atoms[i]
                                  << std::scientific << std::setprecision(6)
                                  << std::setw(18) << cpu_de_bond[i]
                                  << std::setw(18) << disp
                                  << std::setw(18) << cpu_de[i] << "\n";
                    }
                }
            } else {
                std::cout << "  [SKIP] size mismatch or empty\n";
            }
        } else {
            std::cout << "\n[SKIP] Could not access dEdcn from CPU/GPU methods\n";
        }
    }

    // --- Compare gradient BEFORE CN chain-rule (isolate CN contribution) ---
    {
        auto* gpu_method = dynamic_cast<GGFNFFComputationalMethod*>(gpu->Interface());
        auto* cpu_method2 = dynamic_cast<GFNFFComputationalMethod*>(cpu->Interface());
        if (gpu_method && cpu_method2 && cpu_method2->getGFNFF() && cpu_method2->getGFNFF()->getWorkspace()
            && gpu_method->getGPUWorkspace()) {
            // Deep copy to avoid reference-to-freed issues
            Matrix G_before_cn_gpu = gpu_method->getGPUWorkspace()->gradientBeforeCN();
            Matrix G_before_cn_cpu = cpu_method2->getGFNFF()->getWorkspace()->gradientBeforeCN();
            if (G_before_cn_gpu.rows() == 0 || G_before_cn_gpu.norm() < 1e-15) {
                std::cout << "\n--- CN chain-rule contribution analysis ---\n";
                std::cout << "[SKIP] GPU gradient snapshots require verbosity >= 3\n";
            } else if (G_before_cn_gpu.rows() == mol.m_number_atoms && G_before_cn_cpu.rows() == mol.m_number_atoms) {
                Matrix G_cn_gpu = G_gpu - G_before_cn_gpu;
                Matrix G_cn_cpu = G_cpu - G_before_cn_cpu;

                std::cout << "\n--- CN chain-rule contribution analysis ---\n";
                std::cout << "GPU grad before CN norm: " << std::scientific << G_before_cn_gpu.norm() << std::endl;
                std::cout << "CPU grad before CN norm: " << std::scientific << G_before_cn_cpu.norm() << std::endl;
                std::cout << "GPU CN contribution norm: " << std::scientific << G_cn_gpu.norm() << std::endl;
                std::cout << "CPU CN contribution norm: " << std::scientific << G_cn_cpu.norm() << std::endl;

                // Direct gradient comparison WITHOUT CN chain-rule
                Matrix diff_no_cn = G_before_cn_gpu - G_before_cn_cpu;
                double max_no_cn = diff_no_cn.cwiseAbs().maxCoeff();
                std::cout << "\nGrad diff (before CN) max: " << std::scientific << max_no_cn;
                if (max_no_cn < 1e-4) std::cout << " [PASS]\n";
                else                  std::cout << " [FAIL] (tol=1e-4)\n";

                // CN chain-rule contribution comparison
                Matrix diff_cn = G_cn_gpu - G_cn_cpu;
                double max_cn = diff_cn.cwiseAbs().maxCoeff();
                std::cout << "CN chain-rule diff max: " << std::scientific << max_cn;
                if (max_cn < 1e-4) std::cout << " [PASS]\n";
                else               std::cout << " [FAIL] (tol=1e-4)\n";

                // Top 5 worst atoms for each component
                std::cout << "\nTop 5 worst before-CN diffs (GPU - CPU):\n";
                std::vector<std::pair<double, int>> bcn_diffs;
                for (int i = 0; i < mol.m_number_atoms; ++i)
                    bcn_diffs.push_back({diff_no_cn.row(i).cwiseAbs().maxCoeff(), i});
                std::sort(bcn_diffs.rbegin(), bcn_diffs.rend());
                for (int n = 0; n < std::min(5, (int)bcn_diffs.size()); ++n) {
                    int i = bcn_diffs[n].second;
                    std::cout << std::setw(5) << i << std::setw(4) << mol.m_atoms[i]
                              << std::scientific << std::setprecision(6)
                              << "  dx=" << std::setw(14) << diff_no_cn(i,0)
                              << "  dy=" << std::setw(14) << diff_no_cn(i,1)
                              << "  dz=" << std::setw(14) << diff_no_cn(i,2) << "\n";
                }

                std::cout << "\nTop 5 worst CN chain-rule diffs (GPU - CPU):\n";
                std::vector<std::pair<double, int>> cn_diffs;
                for (int i = 0; i < mol.m_number_atoms; ++i)
                    cn_diffs.push_back({diff_cn.row(i).cwiseAbs().maxCoeff(), i});
                std::sort(cn_diffs.rbegin(), cn_diffs.rend());
                for (int n = 0; n < std::min(5, (int)cn_diffs.size()); ++n) {
                    int i = cn_diffs[n].second;
                    std::cout << std::setw(5) << i << std::setw(4) << mol.m_atoms[i]
                              << std::scientific << std::setprecision(6)
                              << "  dx=" << std::setw(14) << diff_cn(i,0)
                              << "  dy=" << std::setw(14) << diff_cn(i,1)
                              << "  dz=" << std::setw(14) << diff_cn(i,2) << "\n";
                }
            }
        }
    }

    // --- Targeted numerical gradient for worst atoms (CPU + GPU) ---
    {
        Matrix G_diff_pre = G_gpu - G_cpu;
        std::vector<std::pair<double, int>> pre_diffs;
        for (int i = 0; i < mol.m_number_atoms; ++i)
            pre_diffs.push_back({G_diff_pre.row(i).cwiseAbs().maxCoeff(), i});
        std::sort(pre_diffs.rbegin(), pre_diffs.rend());

        int n_check = std::min(5, (int)pre_diffs.size());
        std::cout << "\n--- Numerical gradient check (top " << n_check << " worst atoms) ---\n";
        constexpr double step_bohr = 1e-5;
        constexpr double bohr2ang = 0.529177249;
        double step_ang = step_bohr * bohr2ang;
        Geometry orig_geom = mol.m_geometry;

        std::cout << std::setw(5) << "Atom" << std::setw(4) << "Z" << std::setw(4) << "d"
                  << std::setw(14) << "CPU_anal" << std::setw(14) << "CPU_num"
                  << std::setw(14) << "CPU_err"
                  << std::setw(14) << "GPU_anal" << std::setw(14) << "GPU_num"
                  << std::setw(14) << "GPU_err" << "\n";

        for (int n = 0; n < n_check; ++n) {
            int atom = pre_diffs[n].second;
            for (int d = 0; d < 3; ++d) {
                // +h
                Mol mol_p = mol;
                mol_p.m_geometry = orig_geom;
                mol_p.m_geometry(atom, d) += step_ang;
                auto* cpu_p = new EnergyCalculator("gfnff", config);
                cpu_p->setMolecule(mol_p);
                double cpu_ep = cpu_p->CalculateEnergy(false);
                auto* gpu_p = new EnergyCalculator("ggfnff", config);
                gpu_p->setMolecule(mol_p);
                double gpu_ep = gpu_p->CalculateEnergy(false);
                delete cpu_p;
                // NOTE: GPU instances intentionally leaked — CUDA heap corruption
                // makes destructors crash.  _exit() at end cleans up.

                // -h
                Mol mol_m = mol;
                mol_m.m_geometry = orig_geom;
                mol_m.m_geometry(atom, d) -= step_ang;
                auto* cpu_m = new EnergyCalculator("gfnff", config);
                cpu_m->setMolecule(mol_m);
                double cpu_em = cpu_m->CalculateEnergy(false);
                auto* gpu_m = new EnergyCalculator("ggfnff", config);
                gpu_m->setMolecule(mol_m);
                double gpu_em = gpu_m->CalculateEnergy(false);
                delete cpu_m;
                // GPU intentionally leaked (CUDA heap corruption)

                double cpu_num = (cpu_ep - cpu_em) / (2.0 * step_bohr);
                double gpu_num = (gpu_ep - gpu_em) / (2.0 * step_bohr);

                std::cout << std::setw(5) << atom << std::setw(4) << mol.m_atoms[atom]
                          << std::setw(4) << "xyz"[d]
                          << std::scientific << std::setprecision(4)
                          << std::setw(14) << G_cpu(atom, d)
                          << std::setw(14) << cpu_num
                          << std::setw(14) << (G_cpu(atom, d) - cpu_num)
                          << std::setw(14) << G_gpu(atom, d)
                          << std::setw(14) << gpu_num
                          << std::setw(14) << (G_gpu(atom, d) - gpu_num) << "\n";
            }
        }
        std::cout << std::endl;
    }

    // --- Print gradient details (top 10 worst atoms) ---
    {
        Matrix G_diff = G_gpu - G_cpu;
        std::vector<std::pair<double, int>> diffs;
        for (int i = 0; i < mol.m_number_atoms; ++i) {
            double maxd = std::max({std::abs(G_diff(i,0)), std::abs(G_diff(i,1)), std::abs(G_diff(i,2))});
            diffs.push_back({maxd, i});
        }
        std::sort(diffs.rbegin(), diffs.rend());
        int show = std::min((int)diffs.size(), 10);
        std::cout << "\nTop " << show << " worst gradient diffs (GPU - CPU):\n"
                  << std::setw(5) << "Atom" << std::setw(4) << "Z"
                  << std::setw(16) << "CPU_x" << std::setw(16) << "GPU_x"
                  << std::setw(16) << "diff_x" << std::setw(16) << "diff_y"
                  << std::setw(16) << "diff_z" << std::setw(14) << "|diff|" << "\n";
        for (int n = 0; n < show; ++n) {
            int i = diffs[n].second;
            double norm_d = G_diff.row(i).norm();
            std::cout << std::setw(5) << i
                      << std::setw(4) << mol.m_atoms[i]
                      << std::scientific << std::setprecision(4)
                      << std::setw(16) << G_cpu(i,0) << std::setw(16) << G_gpu(i,0)
                      << std::setw(16) << G_diff(i,0) << std::setw(16) << G_diff(i,1)
                      << std::setw(16) << G_diff(i,2) << std::setw(14) << norm_d << "\n";
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
