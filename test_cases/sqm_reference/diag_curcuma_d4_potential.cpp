/*
 * Step 3a diagnostic: curcuma's D4 SCF potential at tblite's converged charges.
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026.
 *
 * Reads a release_tblite dump (geometry + tblite's converged atomic_charges +
 * vat_tblite + gf2.vat_extra) and computes curcuma's dE_D4/dq at *tblite's*
 * Mulliken charges, mirroring xtb_native::calcDispersionEnergy / addDispersionPotential
 * exactly (per_reference_charge=true, setD4CovalentCN(true),
 * setUseD4SingleShotEEQ(false), setTopologyCharges(q_tblite)).
 *
 * Compares against (vat_tblite − gf2.vat_extra), which equals tblite's D4 SCF
 * potential at tblite's converged charges (the multipole-machinery diff is
 * ~1e-7 at fixed density — proven by diff_multipole_potential.py).
 *
 * Outcome interpretation:
 *   max|diff| ≲ 1e-7  → D4 SCF potential is exact; residual is NOT in D4.
 *   max|diff| ~5e-5   → D4 SCF potential is the polarity-driven Fock source.
 *
 * Usage:  diag_curcuma_d4_potential <dump.json> <xyz.xyz> [<dump.json> <xyz.xyz> ...]
 */

#include "src/core/config_manager.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/global.h"

#include "external/json.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using json = nlohmann::json;
static constexpr double AA_TO_AU = 1.0 / 0.529177210903;

static int Z_of(const std::string& s) {
    static const std::map<std::string, int> m = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},
        {"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
        {"S",16},{"Cl",17},{"Ar",18}};
    auto it = m.find(s);
    return (it == m.end()) ? -1 : it->second;
}

static bool readXYZ(const std::string& path,
                    std::vector<int>& atoms, Matrix& geom_ang)
{
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return false; }
    int nat = 0;
    f >> nat;
    std::string line;
    std::getline(f, line);  // rest of first line
    std::getline(f, line);  // comment line
    atoms.clear();
    geom_ang = Matrix::Zero(nat, 3);
    for (int i = 0; i < nat; ++i) {
        std::string sym;
        double x, y, z;
        if (!(f >> sym >> x >> y >> z)) { std::cerr << "xyz parse error at atom " << i << "\n"; return false; }
        atoms.push_back(Z_of(sym));
        geom_ang(i, 0) = x; geom_ang(i, 1) = y; geom_ang(i, 2) = z;
    }
    return true;
}

static void analyse(const std::string& dump_path, const std::string& xyz_path)
{
    std::ifstream f(dump_path);
    if (!f) { std::cerr << "cannot open " << dump_path << "\n"; return; }
    json d;
    f >> d;

    std::vector<int> atoms;
    Matrix geom_ang;
    if (!readXYZ(xyz_path, atoms, geom_ang)) return;
    const int nat = static_cast<int>(atoms.size());
    if (static_cast<int>(d["atomic_charges"].size()) != nat) {
        std::cerr << "atom-count mismatch between dump and xyz\n";
        return;
    }
    const Matrix geom_bohr = geom_ang * AA_TO_AU;

    // tblite's converged Mulliken charges
    Vector q_tblite(nat);
    for (int i = 0; i < nat; ++i) q_tblite(i) = d["atomic_charges"][i].get<double>();

    // Mirror xtb_native::calcDispersionEnergy GFN2 setup
    json cfg;
    cfg["d4_s6"] = 1.0; cfg["d4_s8"] = 2.7; cfg["d4_a1"] = 0.52;
    cfg["d4_a2"] = 5.0; cfg["d4_alp"] = 16.0;
    ConfigManager cm("d4param", cfg);
    D4ParameterGenerator gen(cm);
    gen.setUseD4SingleShotEEQ(false);
    gen.setD4CovalentCN(true);
    gen.setTopologyCharges(q_tblite);

    curcuma::dispersion::D4Params p;
    p.s6 = 1.0; p.s8 = 2.7; p.a1 = 0.52; p.a2 = 5.0;
    p.s9 = 5.0; p.alpha = 16.0;
    p.damping = curcuma::dispersion::DampingFormula::StandardBJ_D4;
    p.per_reference_charge = true;
    curcuma::dispersion::D4Evaluator eval(&gen, p);

    gen.GenerateParameters(atoms, geom_bohr);
    Matrix grad;
    Vector dEdCN, dEdq;
    eval.computeEnergyAndGradient(atoms, geom_bohr,
                                  /*with_gradient=*/true,
                                  grad, dEdCN, dEdq, /*with_dEdq=*/true);

    // tblite's D4 SCF potential implied by the dump:
    //   implied = vat_tblite − gf2.vat_extra
    // (multipole-machinery residual at fixed density is ~1e-7, so this is
    //  effectively tblite's D4 contribution to vat at tblite's converged charges)
    Vector vat_tbl(nat), vat_extra(nat);
    for (int i = 0; i < nat; ++i) {
        vat_tbl(i)   = d["vat_tblite"][i].get<double>();
        vat_extra(i) = d["gf2"]["vat_extra"][i].get<double>();
    }
    const Vector implied = vat_tbl - vat_extra;

    std::cout << dump_path << "  (nat=" << nat << ")\n";
    std::cout << "  atom   curcuma_dE_D4/dq    implied_D4_vat       diff         sum\n";
    double max_diff = 0.0, max_sum = 0.0;
    for (int i = 0; i < nat; ++i) {
        const double diff = dEdq(i) - implied(i);
        const double sum  = dEdq(i) + implied(i);
        if (std::abs(diff) > max_diff) max_diff = std::abs(diff);
        if (std::abs(sum)  > max_sum)  max_sum  = std::abs(sum);
        std::printf("  %3d  %18.10e  %18.10e  %10.2e  %10.2e\n",
                    i, dEdq(i), implied(i), diff, sum);
    }
    std::printf("  max|diff|=%.2e   max|sum|=%.2e%s\n",
                max_diff, max_sum,
                (max_sum < max_diff * 0.5 && max_diff > 1e-9) ? "   <-- sign flip" : "");
}

int main(int argc, char** argv)
{
    if (argc < 3 || (argc - 1) % 2 != 0) {
        std::cerr << "usage: diag_curcuma_d4_potential <dump.json> <xyz.xyz> [<dump.json> <xyz.xyz> ...]\n";
        return 2;
    }
    for (int i = 1; i < argc; i += 2) analyse(argv[i], argv[i + 1]);
    return 0;
}
