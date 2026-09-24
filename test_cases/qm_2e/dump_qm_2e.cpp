/*
 * Native DFT 2-electron integral dumper (WP2 validation suite).
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Standalone binary mirroring dump_dft_1e.cpp: reads an XYZ, builds the native
 * DFT basis + 1e integrals, then builds the 4-centre ERI tensor (chemists'
 * (mu nu | lam sig), McMurchie-Davidson, cartesian 6d) via DFT::cartesianERI(),
 * and the Coulomb (J) / exchange (K) matrices from a dummy closed-shell density
 * P. Emits everything as JSON on stdout so scripts/diff_dft_2e.py can compare
 * against the independent Python MD witness (scripts/dft_2e_python_ints.py)
 * element-wise, plus run the internal-consistency checks (8-fold ERI symmetry,
 * J/K symmetric, Tr(P.J)==Tr(P.K) single-orbital identity).
 *
 * The dummy density is built IDENTICALLY to the Python witness:
 *   P = 2 * c * c^T,  c = X[:,k],  X = S^{-1/2},  k = argmax column norm of X,
 * with X from an ascending-eigenvalue symmetric orthogonalization (Eigen
 * SelfAdjointEigenSolver sorts ascending; the Python witness sorts its Jacobi
 * spectrum ascending to match). Column norms are sign-invariant, so k -- and
 * hence P, J, K -- are solver-independent and match element-wise.
 *
 *   dump_dft_2e <input.xyz> [--basis NAME] [--charge Q] [--spin S]
 *
 * Output schema:
 *   { "molecule":{name,natoms,atoms:[{z,x,y,z(Bohr)}]},
 *     "basis","cartesian_d":true,"nbf","num_electrons","nuclear_repulsion",
 *     "eri_order":"mu_nu_lam_sig","S":[[...]],"ERI":[...flat n^4...],
 *     "P":[[...]],"J":[[...]],"K":[[...]],
 *     "dummy_density":"2*X[:,k]*X[:,k]^T, X=S^-1/2, k=argmax col norm" }
 *
 * Claude Generated (WP2). GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/dft.h"
#include "src/core/energy_calculators/qm_methods/dft_integrals.hpp"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "src/core/units.h"

#include "external/json.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
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
        if (Z < 0) { std::cerr << "unsupported element " << sym << " (WP2: H-Ne)\n"; return false; }
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

// Dummy closed-shell single-orbital density P = 2 c c^T, where c is the
// ones-vector S-orthonormalized by the SCALAR c^T S c = 1 (no eigendecomposition).
// This is solver-independent and robust to degenerate S eigensubspaces (which
// made an argmax-column-of-S^{-1/2} pick ambiguous on symmetric molecules such
// as linear BeH2). P is rank-1, so Tr(P J) == Tr(P K) holds for ANY ERI by
// dummy-index relabeling; Tr(P S) = 2 c^T S c = 2 (one spatial orbital, 2
// electrons). The fixed ones-vector exercises every AO index in J/K. Matches
// the Python witness (dummy_density) element-wise.
Matrix dummyDensity(const Matrix& S)
{
    const int n = (int)S.rows();
    Matrix c = Matrix::Ones(n, 1);
    const double q = (c.transpose() * S * c)(0, 0);  // c^T S c
    if (q <= 0.0)
        throw std::runtime_error("c^T S c non-positive for dummy density");
    c /= std::sqrt(q);
    return 2.0 * c * c.transpose();
}

}  // namespace

int main(int argc, char** argv)
{
    // Silence the global CurcumaLogger so only the JSON document reaches stdout.
    CurcumaLogger::set_verbosity(0);

    if (argc < 2) {
        std::cerr << "usage: dump_dft_2e <input.xyz> [--basis NAME] [--charge Q] [--spin S]\n";
        return 2;
    }

    std::string xyz = argv[1];
    std::string basis = "def2-SVP";
    double charge = 0.0;
    int spin = 0;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--basis" && i + 1 < argc) basis = argv[++i];
        else if (a == "--charge" && i + 1 < argc) charge = std::stod(argv[++i]);
        else if (a == "--spin" && i + 1 < argc) spin = std::stoi(argv[++i]);
    }

    std::vector<int> atoms;
    std::vector<double> coord_ang;
    std::string name;
    if (!readXYZ(xyz, atoms, coord_ang, name)) return 1;
    const int nat = (int)atoms.size();

    // The ERI kernel gate runs in the CARTESIAN 6d basis (the Python witness
    // builds cartesian too), so force cartesian_d=true here regardless of any
    // future default. DFT::cartesianERI() always builds in m_gto_basis (the
    // flat cartesian basis), and overlapMatrix() with cartesian_d=true returns
    // the cartesian S -- the two must be in the same AO order for P/J/K.
    json cfg = json::object();
    cfg["basis"] = basis;
    cfg["cartesian_d"] = true;

    DFT dft(DFTFunctional::HF, cfg);
    if (!dft.QMInterface::InitialiseMolecule(atoms.data(), coord_ang.data(), nat, charge, spin)) {
        std::cerr << "DFT InitialiseMolecule failed\n";
        return 1;
    }

    const int n = dft.nbf();
    const Matrix S = dft.overlapMatrix();  // cartesian (cartesian_d=true)

    const dft1e::ERITensor& eri = dft.cartesianERI();
    if (eri.empty() || eri.n() != n) {
        std::cerr << "ERI build failed or size mismatch (n=" << n << ")\n";
        return 1;
    }

    const Matrix P = dummyDensity(S);
    const Matrix J = dft1e::buildCoulomb(eri, P);
    const Matrix K = dft1e::buildExchange(eri, P);

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
    out["cartesian_d"] = true;
    out["nbf"] = n;
    out["num_electrons"] = dft.numElectrons();
    // Nuclear repulsion from the scaffold Calculation() (sets m_total_energy).
    dft.Calculation(false);
    out["nuclear_repulsion"] = dft.TotalEnergy();
    out["eri_order"] = "mu_nu_lam_sig";
    out["S"] = matrixToJson(S);

    // Flat n^4 ERI in row-major order index = ((mu*n+nu)*n+lam)*n+sig,
    // identical to the Python witness and to ERITensor::operator().
    json eri_flat = json::array();
    eri_flat.get_ref<json::array_t&>().reserve((size_t)n * n * n * n);
    for (int mu = 0; mu < n; ++mu)
        for (int nu = 0; nu < n; ++nu)
            for (int lam = 0; lam < n; ++lam)
                for (int sig = 0; sig < n; ++sig)
                    eri_flat.push_back(eri(mu, nu, lam, sig));
    out["ERI"] = std::move(eri_flat);

    out["P"] = matrixToJson(P);
    out["J"] = matrixToJson(J);
    out["K"] = matrixToJson(K);
    out["dummy_density"] = "2*c*c^T, c=ones normalized so c^T S c=1 (rank-1, no eig)";

    std::cout << out.dump() << "\n";
    return 0;
}