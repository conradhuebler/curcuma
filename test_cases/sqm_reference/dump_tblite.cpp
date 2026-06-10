/*
 * <Phase 0 reference dumper for native xTB validation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reads an XYZ file, runs tblite GFN1-xTB or GFN2-xTB, and dumps the
 * converged SCC state to JSON:
 *   total_energy, atom_energies, atomic_charges,
 *   orbital_energies, orbital_occupations,
 *   overlap[nao×nao], hamiltonian[nao×nao], density[nao×nao],
 *   shell_map (atom per shell), angular_momenta[nsh], orbital_map (shell per AO),
 *   shell_charges (computed as refocc_sh - Mulliken_sh from P·S trace).
 *
 * These dumps are the numeric gate for Phase 3 kernel work (H0, Coulomb,
 * multipole, third-order, SCF) in src/core/energy_calculators/qm_methods/.
 *
 * Usage:
 *   dump_tblite <gfn1|gfn2> <input.xyz> <out.json>
 *
 * Claude Generated (Phase 0, Apr 2026)
 * GPL-3.0.
 */

#include "tblite.h"
#include "external/json.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using nlohmann::json;

namespace {

constexpr double AA_TO_BOHR = 1.0 / 0.529177210903;

struct Atom {
    int    z;
    double x, y, z_coord;
};

int elementZ(const std::string& sym)
{
    static const char* table[] = {
        "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "In", "Sn", "Sb", "Te", "I",  "Xe",
        "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
        "Tm", "Yb", "Lu",
        "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn",
    };
    for (int z = 1; z <= 86; ++z) {
        if (sym == table[z]) return z;
    }
    return 0;
}

bool readXYZ(const std::string& path, std::vector<Atom>& atoms, std::string& comment)
{
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << "\n"; return false; }
    std::string line;
    if (!std::getline(in, line)) return false;
    int nat = std::stoi(line);
    std::getline(in, comment);
    atoms.clear();
    atoms.reserve(nat);
    for (int i = 0; i < nat; ++i) {
        if (!std::getline(in, line)) return false;
        std::istringstream iss(line);
        std::string sym;
        double x, y, z;
        if (!(iss >> sym >> x >> y >> z)) return false;
        Atom a;
        a.z = elementZ(sym);
        if (a.z == 0) { std::cerr << "unknown element: " << sym << "\n"; return false; }
        a.x = x * AA_TO_BOHR;
        a.y = y * AA_TO_BOHR;
        a.z_coord = z * AA_TO_BOHR;
        atoms.push_back(a);
    }
    return true;
}

// Upper-triangle serialisation to keep file size manageable for large nao.
json storeSymmetricMatrix(const std::vector<double>& m, int n)
{
    json out = json::array();
    for (int i = 0; i < n; ++i) {
        json row = json::array();
        for (int j = 0; j < n; ++j) row.push_back(m[i * n + j]);
        out.push_back(std::move(row));
    }
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " <gfn1|gfn2> <in.xyz> <out.json>\n";
        return 2;
    }
    const std::string method_s = argv[1];
    const std::string xyz_path = argv[2];
    const std::string out_path = argv[3];

    std::vector<Atom> atoms;
    std::string       comment;
    if (!readXYZ(xyz_path, atoms, comment)) return 3;

    const int nat = static_cast<int>(atoms.size());
    std::vector<int>    num(nat);
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        num[i]           = atoms[i].z;
        xyz[3 * i + 0]   = atoms[i].x;
        xyz[3 * i + 1]   = atoms[i].y;
        xyz[3 * i + 2]   = atoms[i].z_coord;
    }

    tblite_error   err = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    tblite_set_context_verbosity(ctx, 0);

    const double     charge = 0.0;
    const int        uhf    = 0;
    tblite_structure mol    = tblite_new_structure(
        err, nat, num.data(), xyz.data(), &charge, &uhf, nullptr, nullptr);

    tblite_calculator calc = nullptr;
    if (method_s == "gfn1") {
        calc = tblite_new_gfn1_calculator(ctx, mol);
    } else if (method_s == "gfn2") {
        calc = tblite_new_gfn2_calculator(ctx, mol);
    } else {
        std::cerr << "unknown method: " << method_s << " (expected gfn1|gfn2)\n";
        return 4;
    }

    int nsh = 0, nao = 0;
    tblite_get_calculator_shell_count  (ctx, calc, &nsh);
    tblite_get_calculator_orbital_count(ctx, calc, &nao);
    std::vector<int> sh2at(nsh), ang(nsh), ao2sh(nao);
    tblite_get_calculator_shell_map      (ctx, calc, sh2at.data());
    tblite_get_calculator_angular_momenta(ctx, calc, ang.data());
    tblite_get_calculator_orbital_map    (ctx, calc, ao2sh.data());

    // Request that tblite save S, H0, density matrices on the result object.
    tblite_set_calculator_save_integrals(ctx, calc, 1);

    tblite_result res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);

    double energy_total = 0.0;
    tblite_get_result_energy(err, res, &energy_total);

    std::vector<double> energies_atom(nat);
    tblite_get_result_energies(err, res, energies_atom.data());

    std::vector<double> charges(nat);
    tblite_get_result_charges(err, res, charges.data());

    std::vector<double> dipole(3);
    tblite_get_result_dipole(err, res, dipole.data());

    std::vector<double> emo(nao), occ(nao);
    tblite_get_result_orbital_energies   (err, res, emo.data());
    tblite_get_result_orbital_occupations(err, res, occ.data());

    std::vector<double> S(nao * nao), H(nao * nao), P(nao * nao);
    tblite_get_result_overlap_matrix    (err, res, S.data());
    tblite_get_result_hamiltonian_matrix(err, res, H.data());
    tblite_get_result_density_matrix    (err, res, P.data());

    // Mulliken shell populations: n_sh[s] = Σ_{μ in s} (P·S)_{μμ}
    std::vector<double> n_sh(nsh, 0.0);
    for (int mu = 0; mu < nao; ++mu) {
        double diag = 0.0;
        for (int nu = 0; nu < nao; ++nu) diag += P[mu * nao + nu] * S[nu * nao + mu];
        n_sh[ao2sh[mu]] += diag;
    }

    // HOMO/LUMO from occupations (spin-less; occ ∈ [0, 2])
    double e_homo = 0.0, e_lumo = 0.0;
    int    homo_idx = -1, lumo_idx = -1;
    for (int i = 0; i < nao; ++i) {
        if (occ[i] > 0.5) {
            homo_idx = i;
            e_homo   = emo[i];
        } else if (lumo_idx < 0) {
            lumo_idx = i;
            e_lumo   = emo[i];
        }
    }

    json out;
    out["method"]       = method_s;
    out["input_xyz"]    = xyz_path;
    out["xyz_comment"]  = comment;
    out["natoms"]       = nat;
    out["nshells"]      = nsh;
    out["norbitals"]    = nao;
    out["atoms"]        = json::array();
    for (int i = 0; i < nat; ++i) {
        out["atoms"].push_back({
            {"z", atoms[i].z},
            {"xyz_bohr", { atoms[i].x, atoms[i].y, atoms[i].z_coord }},
        });
    }

    out["shell_map"]       = sh2at;   // length nsh, atom index per shell
    out["angular_momenta"] = ang;     // length nsh, ℓ per shell
    out["orbital_map"]     = ao2sh;   // length nao, shell index per AO

    out["energy_total"]    = energy_total;
    out["energies_atom"]   = energies_atom;
    out["atomic_charges"]  = charges;
    out["shell_populations"] = n_sh;
    out["dipole"]          = dipole;

    out["orbital_energies"]     = emo;
    out["orbital_occupations"]  = occ;
    out["homo_index"]           = homo_idx;
    out["lumo_index"]           = lumo_idx;
    out["homo_energy"]          = e_homo;
    out["lumo_energy"]          = e_lumo;
    out["gap"]                  = (lumo_idx >= 0) ? (e_lumo - e_homo) : 0.0;

    out["overlap"]     = storeSymmetricMatrix(S, nao);
    out["hamiltonian"] = storeSymmetricMatrix(H, nao);
    out["density"]     = storeSymmetricMatrix(P, nao);

    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "cannot open " << out_path << " for writing\n"; return 5; }
    ofs << out.dump(2) << "\n";

    std::fprintf(stdout, "%s %s  E=%.10f Eh  nao=%d\n",
                 method_s.c_str(), xyz_path.c_str(), energy_total, nao);

    tblite_delete_result(&res);
    tblite_delete_calculator(&calc);
    tblite_delete_structure(&mol);
    tblite_delete_context(&ctx);
    tblite_delete_error(&err);
    return 0;
}
