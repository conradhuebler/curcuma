/**
 * GPU analytical vs numerical gradient check for ggfnff.
 * Claude Generated (March 2026)
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

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    CurcumaLogger::set_verbosity(0);
    std::string mol_file = (argc > 1) ? argv[1] : "A.xyz";
    double step_bohr = 1e-5;
    // m_geometry is in Angstrom; convert step to Angstrom for displacement
    constexpr double bohr2ang = 0.529177249;
    double step_ang = step_bohr * bohr2ang;

    std::cout << "GPU ggfnff: Analytical vs Numerical Gradient\n"
              << "Molecule: " << mol_file << ", step=" << step_bohr << " Bohr\n"
              << std::string(60, '=') << std::endl;

    curcuma::Molecule molecule(mol_file);
    Mol mol = molecule.getMolInfo();
    int N = mol.m_number_atoms;
    std::cout << "Atoms: " << N << std::endl;

    json config = {{"verbosity", 0}, {"threads", 1}, {"ggfnff", json::object()}};

    // Analytical gradient
    auto* calc = new EnergyCalculator("ggfnff", config);
    calc->setMolecule(mol);
    double E0 = calc->CalculateEnergy(true);
    Matrix G_anal = calc->Gradient();
    std::cout << "Analytical energy: " << std::setprecision(12) << E0 << " Eh\n";
    std::cout << "Analytical grad norm: " << std::scientific << G_anal.norm() << "\n\n";

    // Numerical gradient
    Matrix G_num = Matrix::Zero(N, 3);
    Geometry orig_geom = mol.m_geometry; // save original

    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            // +h
            Mol mol_p = mol;
            mol_p.m_geometry(i, d) = orig_geom(i, d) + step_ang;
            auto* c_p = new EnergyCalculator("ggfnff", config);
            c_p->setMolecule(mol_p);
            double ep = c_p->CalculateEnergy(false);
            delete c_p;

            // -h
            Mol mol_m = mol;
            mol_m.m_geometry(i, d) = orig_geom(i, d) - step_ang;
            auto* c_m = new EnergyCalculator("ggfnff", config);
            c_m->setMolecule(mol_m);
            double em = c_m->CalculateEnergy(false);
            delete c_m;

            G_num(i, d) = (ep - em) / (2.0 * step_bohr); // Eh/Bohr
        }
    }

    std::cout << "Numerical grad norm: " << std::scientific << G_num.norm() << "\n\n";

    // Compare
    Matrix G_diff = G_anal - G_num;
    double max_diff = G_diff.cwiseAbs().maxCoeff();
    double norm_ratio = G_anal.norm() / G_num.norm();

    std::cout << "Max component diff: " << std::scientific << max_diff << " Eh/Bohr\n";
    std::cout << "Norm ratio (anal/num): " << std::fixed << std::setprecision(6) << norm_ratio << "\n\n";

    // Per-atom details
    std::cout << std::setw(5) << "Atom" << std::setw(14) << "Anal_x" << std::setw(14) << "Num_x"
              << std::setw(14) << "Anal_y" << std::setw(14) << "Num_y"
              << std::setw(14) << "Anal_z" << std::setw(14) << "Num_z"
              << std::setw(12) << "|diff|" << "\n";
    for (int i = 0; i < N; ++i) {
        double nd = G_diff.row(i).norm();
        std::cout << std::setw(5) << i
                  << std::scientific << std::setprecision(4)
                  << std::setw(14) << G_anal(i,0) << std::setw(14) << G_num(i,0)
                  << std::setw(14) << G_anal(i,1) << std::setw(14) << G_num(i,1)
                  << std::setw(14) << G_anal(i,2) << std::setw(14) << G_num(i,2)
                  << std::setw(12) << nd << "\n";
    }

    bool pass = max_diff < 1e-4;
    std::cout << "\nResult: " << (pass ? "PASS" : "FAIL") << " (max_diff=" << std::scientific << max_diff << ")\n";

    std::cout.flush();
    _exit(pass ? 0 : 1);
}

#else
#include <iostream>
int main() { std::cout << "No CUDA\n"; return 0; }
#endif
