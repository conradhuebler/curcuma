/*
 * Native QM gradient dumper (WP8 validation).
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reads an XYZ, runs `hf` (or `hf-3c`) and prints JSON with the energy, the
 * analytic gradient and -- with --fd -- the central finite difference of the
 * SAME method's energy (h in Bohr), all in Eh/Bohr. check_qm_gradient.py
 * compares the analytic gradient with the FD and with PySCF (+ simple-dftd3).
 *
 *   dump_qm_gradient <xyz> [--method hf|hf-3c] [--basis NAME] [--fd] [--h 1e-4]
 *
 * Claude Generated (Sep 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/hf3c_method.h"
#include "src/core/energy_calculators/qm_methods/qm_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: dump_qm_gradient <xyz> [--method hf|hf-3c] [--basis NAME] [--fd] [--h 1e-4]\n";
        return 1;
    }
    std::string method = "hf", basis = "def2-SVP";
    bool fd = false;
    double h = 1.0e-4;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--method" && i + 1 < argc) method = argv[++i];
        else if (a == "--basis" && i + 1 < argc) basis = argv[++i];
        else if (a == "--fd") fd = true;
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
    }
    CurcumaLogger::set_verbosity(0);

    std::map<std::string, int> Z{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
        { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 } };
    std::ifstream f(argv[1]);
    int n;
    f >> n;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    Mol mol;
    mol.m_number_atoms = n;
    mol.m_atoms.resize(n);
    mol.m_geometry = Geometry::Zero(n, 3);
    mol.m_charge = 0;
    mol.m_spin = 0;
    for (int i = 0; i < n; ++i) {
        std::string s;
        f >> s >> mol.m_geometry(i, 0) >> mol.m_geometry(i, 1) >> mol.m_geometry(i, 2);
        mol.m_atoms[i] = Z[s];
    }

    const json cfg = { { "qm", { { "basis", basis }, { "scf_threshold", 1e-10 }, { "scf_max_iterations", 200 } } } };
    auto make = [&]() -> std::unique_ptr<ComputationalMethod> {
        if (method == "hf-3c") return std::make_unique<HF3CMethod>(cfg);
        return std::make_unique<QMMethod>(QMFunctional::HF, cfg);
    };
    auto m = make();
    m->setMolecule(mol);
    const double e = m->calculateEnergy(true);
    // getGradient() is Eh/Angstrom by contract; report Eh/Bohr.
    const Matrix g = m->getGradient() * CurcumaUnit::Length::BOHR_TO_ANGSTROM;

    json out;
    out["energy"] = e;
    out["natoms"] = n;
    std::vector<double> gv(g.data(), g.data() + 3 * n);  // row-major: atom-major
    out["gradient"] = gv;

    if (fd) {
        const double hA = h * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
        std::vector<double> fdv(3 * n);
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) {
                Geometry gp = mol.m_geometry, gm = mol.m_geometry;
                gp(i, k) += hA;
                gm(i, k) -= hA;
                auto mp = make(), mm = make();
                Mol a = mol, b = mol;
                a.m_geometry = gp;
                b.m_geometry = gm;
                mp->setMolecule(a);
                mm->setMolecule(b);
                fdv[3 * i + k] = (mp->calculateEnergy(false) - mm->calculateEnergy(false)) / (2.0 * h);
            }
        out["fd_gradient"] = fdv;
        out["fd_step_bohr"] = h;
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}
