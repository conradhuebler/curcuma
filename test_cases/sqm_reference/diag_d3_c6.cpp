/*
 * Dump curcuma's interpolated GFN1 D3 per-pair C6 (and CN) for a molecule.
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026 (SQM WP2: GFN1 D3 C6 localisation).
 *
 * Builds D3ParameterGenerator::createForGFN1(), runs GenerateParameters on the
 * molecule, and prints the per-pair c6 = interpolateC6(Zi,Zj, atom_i, atom_j)
 * together with the coordination numbers. Diffed against s-dftd3's get_atomic_c6
 * (dump_dftd3_atomic_c6) to find where curcuma's interpolated D3 C6 diverges --
 * the term WP2 isolated as the entire GFN1 dispersion residual.
 *
 *   diag_d3_c6 <in.xyz> [out.json]
 *
 * JSON: { "nat":N, "cn":[...], "pairs":[ {"i","j","zi","zj","c6"} ] }
 */
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"

#include "external/json.hpp"

#include <Eigen/Dense>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {
int Z_of(const std::string& s) {
    static const std::map<std::string, int> m = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},
        {"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
        {"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},{"Br",35},{"I",53}};
    auto it = m.find(s);
    return (it == m.end()) ? -1 : it->second;
}
}  // namespace

int main(int argc, char** argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: %s <in.xyz> [out.json]\n", argv[0]); return 2; }
    const std::string xyz = argv[1];
    const std::string out = (argc >= 3) ? argv[2] : "";

    std::ifstream f(xyz);
    if (!f) { std::fprintf(stderr, "cannot open %s\n", xyz.c_str()); return 3; }
    int nat = 0; f >> nat; std::string line; std::getline(f, line); std::getline(f, line);
    std::vector<int> atoms; Eigen::MatrixXd geom(nat, 3);
    for (int i = 0; i < nat; ++i) {
        std::string sym; double x, y, z;
        if (!(f >> sym >> x >> y >> z)) { std::fprintf(stderr, "parse error atom %d\n", i); return 3; }
        int zz = Z_of(sym);
        if (zz < 0) { std::fprintf(stderr, "unknown element %s\n", sym.c_str()); return 3; }
        atoms.push_back(zz); geom(i, 0) = x; geom(i, 1) = y; geom(i, 2) = z;  // Angstrom
    }

    D3ParameterGenerator gen = D3ParameterGenerator::createForGFN1();
    gen.GenerateParameters(atoms, geom);
    json p = gen.getParameters();

    json outj; outj["nat"] = nat;
    json pairs = json::array();
    for (const auto& pr : p["d3_dispersion_pairs"]) {
        json o;
        o["i"]  = pr["i"]; o["j"] = pr["j"];
        o["zi"] = pr["element_i"]; o["zj"] = pr["element_j"];
        o["c6"] = pr["c6"];
        pairs.push_back(o);
    }
    outj["pairs"] = pairs;

    // Per-element reference blocks (curcuma's loaded tables) for the unique
    // elements present: number of references, reference CN, and the self C6
    // block c6(Z,Z,ri,rj). Diffed against s-dftd3 model%ref / model%cn / model%c6.
    json refdata;
    std::map<int, bool> seen;
    for (int z : atoms) {
        if (seen.count(z)) continue;
        seen[z] = true;
        int nref = gen.getNumberofReferences(z);
        json rd;
        rd["nref"] = nref;
        json cn = json::array();
        for (int r = 0; r < 7; ++r) cn.push_back(gen.getReferenceCN(z, r));
        rd["refcn"] = cn;
        json c6self = json::array();
        for (int ri = 0; ri < 7; ++ri) {
            json row = json::array();
            for (int rj = 0; rj < 7; ++rj) row.push_back(gen.getC6(z, z, ri, rj));
            c6self.push_back(row);
        }
        rd["c6self"] = c6self;
        refdata[std::to_string(z)] = rd;
    }
    outj["refdata"] = refdata;

    std::string s = outj.dump(2);
    if (!out.empty()) { std::ofstream os(out); os << s << "\n"; }
    else std::cout << s << "\n";
    return 0;
}
