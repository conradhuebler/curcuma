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
 *                    [--ri off|j|jk] [--ri_metric cholesky|eigen] [--cartesian] [--threads N]
 *
 * --method mp2|scs-mp2|sos-mp2 (Oct 2026, MP2 gradient): native canonical MP2 (MP2Method);
 * --frozen_core adds -qm.frozen_core true; --qm KEY VALUE sets any further qm-module key (VALUE parsed as
 * JSON when possible, e.g. --qm scf_direct on, --qm mp2_batch_mb 0.001).
 *
 * --ri (Oct 2026, RI gradient): density-fitted SCF and its analytic gradient; the FD then
 * differentiates the RI energy.
 *
 * Claude Generated (Sep 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/hf3c_method.h"
#include "src/core/energy_calculators/qm_methods/qm_method.h"
#include "src/core/energy_calculators/qm_methods/mp2_method.h"
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
    std::string method = "hf", basis = "def2-SVP", ri = "off", ri_metric = "cholesky";
    bool fd = false, cartesian = false, frozen_core = false;
    json qm_extra = json::object();
    int threads = 0;  // 0: engine default
    double h = 1.0e-4;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--method" && i + 1 < argc) method = argv[++i];
        else if (a == "--basis" && i + 1 < argc) basis = argv[++i];
        else if (a == "--fd") fd = true;
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
        else if (a == "--ri" && i + 1 < argc) ri = argv[++i];
        else if (a == "--ri_metric" && i + 1 < argc) ri_metric = argv[++i];
        else if (a == "--cartesian") cartesian = true;
        else if (a == "--frozen_core") frozen_core = true;
        else if (a == "--qm" && i + 2 < argc) {
            const std::string key = argv[++i], val = argv[++i];
            const json parsed = json::parse(val, nullptr, false);
            qm_extra[key] = parsed.is_discarded() ? json(val) : parsed;
        }
        else if (a == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
    }
    CurcumaLogger::set_verbosity(0);

    std::map<std::string, int> Z{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
        { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 } };
    std::ifstream f(argv[1]);
    int n = 0;
    if (!(f >> n) || n <= 0) {
        std::cerr << "dump_qm_gradient: cannot read " << argv[1] << "\n";
        return 1;
    }
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

    json cfg = { { "qm", { { "basis", basis }, { "scf_threshold", 1e-10 }, { "scf_max_iterations", 200 } } } };
    if (ri != "off") {
        cfg["qm"]["ri"] = ri;
        cfg["qm"]["ri_metric"] = ri_metric;
    }
    if (cartesian) cfg["qm"]["cartesian_d"] = true;
    if (threads > 0) cfg["qm"]["threads"] = threads;
    if (frozen_core) cfg["qm"]["frozen_core"] = true;
    for (auto it = qm_extra.begin(); it != qm_extra.end(); ++it) cfg["qm"][it.key()] = it.value();
    auto make = [&]() -> std::unique_ptr<ComputationalMethod> {
        if (method == "hf-3c") return std::make_unique<HF3CMethod>(cfg);
        if (method == "mp2" || method == "scs-mp2" || method == "sos-mp2") return std::make_unique<MP2Method>(method, cfg);
        return std::make_unique<QMMethod>(QMFunctional::HF, cfg);
    };
    auto m = make();
    m->setMolecule(mol);
    const double e = m->calculateEnergy(true);
    if (m->hasError()) {
        std::cerr << "dump_qm_gradient: " << m->getErrorMessage() << "\n";
        return 2;
    }
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
