/*
 * Native DFT 1-electron integral dumper (WP1 validation suite).
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Standalone binary mirroring dump_tblite_reference.cpp / diag_curcuma_atomic_c6:
 * reads an XYZ, builds the native DFT basis + 1e integrals (S, T, V, Hc=T+V),
 * and emits them as JSON on stdout so scripts/diff_dft_1e.py can compare
 * against the ORCA reference (scripts/dft_1e_reference.py via orca_2json) and
 * an independent Python integral witness, plus run the internal-consistency
 * checks (symmetry, H=T+V, Tr(P*S)=1, nbf match).
 *
 *   dump_dft_1e <input.xyz> [--basis NAME] [--cartesian_d] [--charge Q] [--spin S]
 *
 * Output schema:
 *   { "molecule":{name,natoms,atoms:[{z,x,y,z(Bohr)}]},
 *     "basis","cartesian_d","nbf","num_electrons","nuclear_repulsion",
 *     "S":[[...]],"T":[[...]],"V":[[...]],"H":[[...]] }
 *
 * Claude Generated (WP1). GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/dft.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "src/core/units.h"

#include "external/json.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;
using namespace CurcumaUnit;

namespace {

constexpr double AA_TO_BOHR = 1.0 / 0.529177210903;  // Angstrom -> Bohr

int elementZ(const std::string& s)
{
    static const std::map<std::string, int> m = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6},
        {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    };
    auto it = m.find(s);
    return (it == m.end()) ? -1 : it->second;
}

bool readXYZ(const std::string& path, std::vector<int>& atoms,
             std::vector<double>& coord_ang, std::string& name)
{
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return false; }
    int nat = 0;
    f >> nat;
    std::string line;
    std::getline(f, line);  // finish nat line
    std::getline(f, line);  // comment line -> molecule name
    name = line;
    atoms.clear();
    coord_ang.assign(3 * nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        std::string sym;
        double x, y, z;
        if (!(f >> sym >> x >> y >> z)) { std::cerr << "xyz parse error\n"; return false; }
        int Z = elementZ(sym);
        if (Z < 0) { std::cerr << "unsupported element " << sym << " (WP1: H-Ne)\n"; return false; }
        atoms.push_back(Z);
        coord_ang[3 * i + 0] = x;
        coord_ang[3 * i + 1] = y;
        coord_ang[3 * i + 2] = z;
    }
    return true;
}

json matrixToJson(const Matrix& M)
{
    json J = json::array();
    for (int i = 0; i < M.rows(); ++i) {
        json row = json::array();
        for (int j = 0; j < M.cols(); ++j)
            row.push_back(M(i, j));
        J.push_back(std::move(row));
    }
    return J;
}

}  // namespace

int main(int argc, char** argv)
{
    // Silence the global CurcumaLogger so only the JSON document reaches stdout.
    CurcumaLogger::set_verbosity(0);

    if (argc < 2) {
        std::cerr << "usage: dump_dft_1e <input.xyz> [--basis NAME] [--cartesian_d] [--charge Q] [--spin S]\n";
        return 2;
    }

    std::string xyz = argv[1];
    std::string basis = "def2-SVP";
    bool cartesian_d = false;
    double charge = 0.0;
    int spin = 0;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--basis" && i + 1 < argc) basis = argv[++i];
        else if (a == "--cartesian_d") cartesian_d = true;
        else if (a == "--charge" && i + 1 < argc) charge = std::stod(argv[++i]);
        else if (a == "--spin" && i + 1 < argc) spin = std::stoi(argv[++i]);
    }

    std::vector<int> atoms;
    std::vector<double> coord_ang;
    std::string name;
    if (!readXYZ(xyz, atoms, coord_ang, name)) return 1;
    const int nat = (int)atoms.size();

    json cfg = json::object();
    cfg["basis"] = basis;
    cfg["cartesian_d"] = cartesian_d;

    DFT dft(DFTFunctional::HF, cfg);
    // The int* QMInterface overload stores geometry (Angstrom) and then calls
    // the virtual (no-arg) InitialiseMolecule(), which builds the basis + 1e
    // integrals.
    if (!dft.QMInterface::InitialiseMolecule(atoms.data(), coord_ang.data(), nat, charge, spin)) {
        std::cerr << "DFT InitialiseMolecule failed\n";
        return 1;
    }

    json out;
    json jat = json::array();
    for (int i = 0; i < nat; ++i) {
        jat.push_back({
            {"z", atoms[i]},
            {"x", coord_ang[3 * i + 0] * AA_TO_BOHR},
            {"y", coord_ang[3 * i + 1] * AA_TO_BOHR},
            {"z", coord_ang[3 * i + 2] * AA_TO_BOHR},
        });
    }
    out["molecule"] = {{"name", name}, {"natoms", nat}, {"atoms", jat}};
    out["basis"] = basis;
    out["cartesian_d"] = cartesian_d;
    out["nbf"] = dft.nbf();
    out["num_electrons"] = dft.numElectrons();
    // Nuclear repulsion from the geometry (independent of the integrals).
    // Recompute here for the record (Calculation would set m_total_energy).
    dft.Calculation(false);
    out["nuclear_repulsion"] = dft.TotalEnergy();
    out["S"] = matrixToJson(dft.overlapMatrix());
    out["T"] = matrixToJson(dft.kineticMatrix());
    out["V"] = matrixToJson(dft.nuclearAttractionMatrix());
    out["H"] = matrixToJson(dft.coreHamiltonian());

    std::cout << out.dump() << "\n";
    return 0;
}