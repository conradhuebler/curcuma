/*
 * Compact tblite reference dumper for the native GFN1/GFN2 validation suite.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Unlike dump_tblite_multipole.cpp (a forensic tool that writes full nao x nao
 * matrices, tens of MB for large systems), this emits ONLY the validation-
 * relevant scalars/vectors so the resulting .ref.json is small and can be
 * committed to git. native-only builds (USE_TBLITE=OFF) then validate against
 * these references without a tblite build, exactly how the GFN-FF suite
 * (test_cases/reference_data/) works.
 *
 * Schema:
 *   { "method", "molecule":{name,natoms,comment,atoms:[{z,x,y,z_coord (Bohr)}]},
 *     "total_energy", "charges":[], "dipole":[x,y,z], "homo_energy","lumo_energy",
 *     "energy_components":{repulsion,dispersion,electronic,interactions,halogen},
 *     "gradient":{norm, values:[[gx,gy,gz],...]} }
 *
 *   dump_tblite_reference <gfn1|gfn2> <input.xyz> [out.json] [--solvent NAME]
 *
 * With --solvent NAME the ALPB implicit-solvation container is attached
 * (tblite_new_alpb_solvation_solvent, refstate=gsolv), so total_energy/gradient
 * become the fully self-consistent solvated values. The solvent name is recorded
 * under "solvent". This is the reference for the native GFN1/GFN2 ALPB path
 * (docs/SQM_SOLVATION_WP.md, WP0.4 / WP1.4).
 *
 * The energy_components block requires the curcuma tblite diagnostic patch
 * (patches/tblite/curcuma-tblite.patch); it is omitted gracefully if absent.
 *
 * API usage mirrors dump_tblite_multipole.cpp (proven build pattern):
 *   tblite_new_structure(err, nat, num, xyz, &charge, &uhf, ...),
 *   tblite_get_singlepoint(ctx, mol, calc, res)  [4-arg form],
 *   orbital count via tblite_get_calculator_orbital_count.
 *
 * Claude Generated (SQM validation suite, 2026-05).
 */

#include "tblite.h"

#include "external/json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

constexpr double AA_TO_BOHR = 1.0 / 0.529177210903;

struct Atom {
    int    z;
    double x, y, z_coord;  // Bohr
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
    for (int z = 1; z <= 86; ++z)
        if (sym == table[z]) return z;
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

#ifdef CURCUMA_TBLITE_ECOMP
// Sum a per-atom energy-component vector to a molecular total. Returns false if
// the getter raised an error (patch absent or container missing), clearing it.
// Only available with the curcuma tblite diagnostic patch
// (patches/tblite/curcuma-tblite.patch); guarded so the dumper builds against a
// vanilla tblite (energy_components is then simply omitted).
bool pullComponent(tblite_error err, tblite_result res, int nat,
                   void (*getter)(tblite_error, tblite_result, double*),
                   double& total)
{
    std::vector<double> v(nat, 0.0);
    getter(err, res, v.data());
    if (tblite_check_error(err)) { tblite_clear_error(err); return false; }
    total = 0.0;
    for (double x : v) total += x;
    return true;
}
#endif  // CURCUMA_TBLITE_ECOMP

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <gfn1|gfn2> <in.xyz> [out.json] [--solvent NAME] [--model alpb|gbsa]\n", argv[0]);
        return 1;
    }
    // Collect positionals + the optional --solvent / --model flags (parse-anywhere).
    std::string solvent;          // empty -> gas phase
    std::string solv_model = "alpb";  // alpb (default) | gbsa; selects tblite version
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--solvent" && i + 1 < argc) { solvent = argv[++i]; }
        else if (a == "--model" && i + 1 < argc) { solv_model = argv[++i]; }
        else pos.push_back(a);
    }
    if (pos.size() < 2) {
        std::fprintf(stderr, "usage: %s <gfn1|gfn2> <in.xyz> [out.json] [--solvent NAME] [--model alpb|gbsa]\n", argv[0]);
        return 1;
    }
    if (solv_model != "alpb" && solv_model != "gbsa") {
        std::fprintf(stderr, "error: --model must be 'alpb' or 'gbsa' (got '%s')\n", solv_model.c_str());
        return 1;
    }
    const std::string method_s = pos[0];
    const std::string xyz_path = pos[1];
    const std::string out_path = (pos.size() >= 3) ? pos[2] : "reference.json";

    std::vector<Atom> atoms;
    std::string comment;
    if (!readXYZ(xyz_path, atoms, comment)) {
        std::fprintf(stderr, "ERROR reading %s\n", xyz_path.c_str());
        return 1;
    }
    const int nat = static_cast<int>(atoms.size());
    std::vector<int>    num(nat);
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        num[i]         = atoms[i].z;
        xyz[3 * i + 0] = atoms[i].x;
        xyz[3 * i + 1] = atoms[i].y;
        xyz[3 * i + 2] = atoms[i].z_coord;
    }

    tblite_error   err = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    tblite_set_context_verbosity(ctx, 0);

    const double charge = 0.0;
    const int    uhf    = 0;
    tblite_structure mol = tblite_new_structure(
        err, nat, num.data(), xyz.data(), &charge, &uhf, nullptr, nullptr);

    tblite_calculator calc = nullptr;
    if (method_s == "gfn2")
        calc = tblite_new_gfn2_calculator(ctx, mol);
    else if (method_s == "gfn1")
        calc = tblite_new_gfn1_calculator(ctx, mol);
    else {
        std::fprintf(stderr, "unknown method %s (expected gfn1|gfn2)\n", method_s.c_str());
        return 1;
    }

    // Attach implicit solvation via the same constructor; the version selects the
    // model+method (ALPB(GFN1)=11, ALPB(GFN2)=12, GBSA(GFN1)=21, GBSA(GFN2)=22).
    // refstate: gsolv=1 -> total_shift = gshift only, matching the native target.
    if (!solvent.empty()) {
        const bool is_gbsa = (solv_model == "gbsa");
        int version;
        if (method_s == "gfn1") version = is_gbsa ? 21 : 11;
        else                    version = is_gbsa ? 22 : 12;
        const int refstate = 1;  // gsolv
        std::vector<char> sbuf(solvent.begin(), solvent.end());
        sbuf.push_back('\0');
        tblite_container cont =
            tblite_new_alpb_solvation_solvent(err, mol, sbuf.data(), version, refstate);
        if (tblite_check_error(err)) {
            std::fprintf(stderr, "tblite %s setup failed for solvent '%s'\n",
                         solv_model.c_str(), solvent.c_str());
            return 1;
        }
        tblite_calculator_push_back(ctx, calc, &cont);
    }

    int nao = 0;
    tblite_get_calculator_orbital_count(ctx, calc, &nao);

    tblite_result res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check_error(err)) {
        std::fprintf(stderr, "tblite singlepoint failed for %s\n", xyz_path.c_str());
        return 1;
    }

    double energy_total = 0.0;
    tblite_get_result_energy(err, res, &energy_total);

    std::vector<double> charges(nat);
    tblite_get_result_charges(err, res, charges.data());

    std::vector<double> dipole(3);
    tblite_get_result_dipole(err, res, dipole.data());

    std::vector<double> emo(nao), occ(nao);
    tblite_get_result_orbital_energies(err, res, emo.data());
    tblite_get_result_orbital_occupations(err, res, occ.data());
    double e_homo = 0.0, e_lumo = 0.0;
    bool have_homo = false, have_lumo = false;
    for (int i = 0; i < nao; ++i) {
        if (occ[i] > 0.5) { e_homo = emo[i]; have_homo = true; }
        else if (!have_lumo) { e_lumo = emo[i]; have_lumo = true; }
    }

    // Gradient (shape [nat][3], Eh/Bohr). get_singlepoint populates it.
    std::vector<double> grad(3 * nat, 0.0);
    tblite_get_result_gradient(err, res, grad.data());
    bool have_grad = !tblite_check_error(err);
    if (!have_grad) tblite_clear_error(err);

    // Energy components (patch-only). Molecular totals only (tblite distributes
    // some terms across containers; combined disp+electronic is the meaningful
    // gate, see diag_curcuma_energy_components header). Requires the curcuma
    // tblite diagnostic patch; without it the block is simply omitted (the
    // solvation validation gates on total_energy, not the decomposition).
    json ec;
#ifdef CURCUMA_TBLITE_ECOMP
    {
        // Capture each component INDEPENDENTLY. A getter that raises -- e.g. the
        // halogen container is legitimately absent for H/C/N/O molecules -- must
        // only drop its own entry, not the whole block. pullComponent clears the
        // tblite error on failure, so subsequent getters still work.
        double v = 0.0;
        if (pullComponent(err, res, nat, tblite_get_result_energy_component_repulsion,    v)) ec["repulsion"]    = v;
        if (pullComponent(err, res, nat, tblite_get_result_energy_component_dispersion,   v)) ec["dispersion"]   = v;
        if (pullComponent(err, res, nat, tblite_get_result_energy_component_electronic,   v)) ec["electronic"]   = v;
        if (pullComponent(err, res, nat, tblite_get_result_energy_component_interactions, v)) ec["interactions"] = v;
        if (pullComponent(err, res, nat, tblite_get_result_energy_component_halogen,      v)) ec["halogen"]      = v;
    }
#endif

    // ----- assemble compact reference -----
    json out;
    out["method"] = method_s;
    if (!solvent.empty()) {
        out["solvent"] = solvent;
        out["solvent_model"] = solv_model;
    }

    json m;
    {
        std::string base = xyz_path;
        auto slash = base.find_last_of("/\\");
        if (slash != std::string::npos) base = base.substr(slash + 1);
        auto dot = base.find_last_of('.');
        if (dot != std::string::npos) base = base.substr(0, dot);
        m["name"] = base;
    }
    m["natoms"]  = nat;
    m["comment"] = comment;
    json jatoms = json::array();
    for (int i = 0; i < nat; ++i)
        jatoms.push_back({{"z", atoms[i].z},
                          {"x", atoms[i].x},
                          {"y", atoms[i].y},
                          {"z_coord", atoms[i].z_coord}});
    m["atoms"] = jatoms;
    out["molecule"] = m;

    out["total_energy"] = energy_total;
    out["charges"]      = charges;
    out["dipole"]       = dipole;
    if (have_homo) out["homo_energy"] = e_homo;
    if (have_lumo) out["lumo_energy"] = e_lumo;
    if (!ec.empty()) out["energy_components"] = ec;

    if (have_grad) {
        double gnorm = 0.0;
        json gvals = json::array();
        for (int i = 0; i < nat; ++i) {
            const double gx = grad[3 * i + 0];
            const double gy = grad[3 * i + 1];
            const double gz = grad[3 * i + 2];
            gnorm += gx * gx + gy * gy + gz * gz;
            gvals.push_back({gx, gy, gz});
        }
        json g;
        g["norm"]   = std::sqrt(gnorm);
        g["values"] = gvals;
        out["gradient"] = g;
    }

    out["_generator"] = "dump_tblite_reference (curcuma SQM validation suite)";

    std::ofstream os(out_path);
    if (!os) {
        std::fprintf(stderr, "ERROR writing %s\n", out_path.c_str());
        return 1;
    }
    os << out.dump(2) << "\n";
    std::fprintf(stderr, "wrote %s  (%s, E=%.8f Eh, %d atoms, grad=%s)\n",
                 out_path.c_str(), method_s.c_str(), energy_total, nat,
                 have_grad ? "yes" : "no");
    return 0;
}
