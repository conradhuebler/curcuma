/*
 * GFN2-D4 WEIGHTED per-atom C6 audit: curcuma weightedC6Gfn2 vs tblite get_atomic_c6.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026. GPL-3.0.
 *
 * Reads dump_dftd4_atomic_c6's JSON (tblite's weighted c6(i,j) + the CN and
 * charges it used) and compares against curcuma's weightedC6Gfn2 evaluated at
 * the SAME charges. Also diffs curcuma's D4 covalent CN against tblite's.
 *
 * This is the step beyond diag_curcuma_d4_c6 (reference C6, ruled out): it
 * audits the WEIGHTED per-atom C6 — the quantity that actually enters the D4
 * energy. A divergence here localises the residual to the CN/charge Gaussian
 * weighting (gwvec). Agreement here rules D4 out entirely, pointing at the
 * electronic half of the disp+elec lump.
 *
 *   diag_curcuma_atomic_c6 [--tol REL] [--quiet] <atomic_c6.json> <xyz.xyz> [...]
 */

#include "src/core/config_manager.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/energy_calculators/dispersion/d4_ncoord.h"
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
        {"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},{"Br",35},{"I",53}};
    auto it = m.find(s);
    return (it == m.end()) ? -1 : it->second;
}

static bool readXYZ(const std::string& path, std::vector<int>& atoms, Matrix& geom_ang) {
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return false; }
    int nat = 0; f >> nat;
    std::string line; std::getline(f, line); std::getline(f, line);
    atoms.clear(); geom_ang = Matrix::Zero(nat, 3);
    for (int i = 0; i < nat; ++i) {
        std::string sym; double x, y, z;
        if (!(f >> sym >> x >> y >> z)) { std::cerr << "xyz parse error\n"; return false; }
        atoms.push_back(Z_of(sym));
        geom_ang(i,0)=x; geom_ang(i,1)=y; geom_ang(i,2)=z;
    }
    return true;
}

// Returns max relative |diff| of the weighted C6, or -1.0 on error.
static double analyse(const std::string& c6_path, const std::string& xyz_path, bool verbose) {
    std::ifstream f(c6_path);
    if (!f) { std::cerr << "cannot open " << c6_path << "\n"; return -1.0; }
    json d;
    try { f >> d; } catch (const std::exception& e) {
        std::cerr << c6_path << ": JSON parse error: " << e.what() << "\n"; return -1.0;
    }
    const int nat = d.at("nat").get<int>();
    std::vector<double> cn_tbl  = d.at("cn").get<std::vector<double>>();
    std::vector<double> qat     = d.at("qat").get<std::vector<double>>();
    std::vector<double> c6_tbl  = d.at("c6").get<std::vector<double>>();
    if ((int)c6_tbl.size() != nat*nat) { std::cerr << "c6 size mismatch\n"; return -1.0; }

    std::vector<int> atoms; Matrix geom_ang;
    if (!readXYZ(xyz_path, atoms, geom_ang)) return -1.0;
    if ((int)atoms.size() != nat) { std::cerr << "nat mismatch xyz vs dump\n"; return -1.0; }
    const Matrix geom_bohr = geom_ang * AA_TO_AU;

    // curcuma's D4 covalent CN (same function GenerateParameters uses, default 30 Bohr cutoff)
    std::vector<double> cn_cur = curcuma::dispersion::computeD4CovalentCN(atoms, geom_bohr);

    // Build the GFN2 generator exactly like xtb_native::calcDispersionEnergy,
    // pinned to tblite's converged charges.
    json cfg; cfg["d4_s6"]=1.0; cfg["d4_s8"]=2.7; cfg["d4_a1"]=0.52; cfg["d4_a2"]=5.0; cfg["d4_alp"]=16.0;
    ConfigManager cm("d4param", cfg);
    D4ParameterGenerator gen(cm);
    gen.setUseD4SingleShotEEQ(false);
    gen.setD4CovalentCN(true);
    Vector qv(nat); for (int i=0;i<nat;++i) qv(i)=qat[i];
    gen.setTopologyCharges(qv);
    gen.GenerateParameters(atoms, geom_bohr);   // sets m_cn_values = cn_cur, rebuilds scaled cache

    // CN comparison
    double max_cn = 0.0; int max_cn_at = -1;
    for (int i = 0; i < nat; ++i) {
        double dcn = std::abs(cn_cur[i] - cn_tbl[i]);
        if (dcn > max_cn) { max_cn = dcn; max_cn_at = i; }
    }

    std::cout << c6_path << "  (nat=" << nat << ")\n";
    if (verbose) {
        std::printf("  %3s %3s %3s %3s  %18s %18s %12s %10s\n",
                    "i","j","Zi","Zj","tblite_c6","curcuma_c6","abs_diff","rel");
    }
    double max_rel = 0.0, max_abs = 0.0;
    int wi=-1, wj=-1;
    int shown = 0;
    for (int i = 0; i < nat; ++i) {
        for (int j = i+1; j < nat; ++j) {
            const double t = c6_tbl[(size_t)i*nat + j];
            const double c = gen.weightedC6Gfn2(atoms[i], atoms[j], i, j, qat[i], qat[j]).c6;
            const double ad = std::abs(c - t);
            const double rel = (std::abs(t) > 1e-12) ? ad/std::abs(t) : ad;
            if (rel > max_rel) { max_rel = rel; wi=i; wj=j; }
            if (ad > max_abs) max_abs = ad;
            if (verbose && shown < 40 && rel > 1e-6) {
                std::printf("  %3d %3d %3d %3d  %18.10e %18.10e %12.3e %10.3e\n",
                            i, j, atoms[i], atoms[j], t, c, c-t, rel);
                ++shown;
            }
        }
    }
    std::printf("  CN: max|diff|=%.3e (atom %d: cur=%.6f tbl=%.6f)\n",
                max_cn, max_cn_at,
                max_cn_at>=0?cn_cur[max_cn_at]:0.0, max_cn_at>=0?cn_tbl[max_cn_at]:0.0);
    std::printf("  weighted C6: max|rel|=%.3e  max|abs|=%.3e  worst pair (%d,%d  Z %d-%d)\n",
                max_rel, max_abs, wi, wj,
                wi>=0?atoms[wi]:-1, wj>=0?atoms[wj]:-1);
    return max_rel;
}

int main(int argc, char** argv) {
    double tol = -1.0; bool verbose = true;
    std::vector<std::string> pos;
    for (int i=1;i<argc;++i){ std::string a=argv[i];
        if (a=="--tol"&&i+1<argc) tol=std::stod(argv[++i]);
        else if (a=="--quiet") verbose=false;
        else if (a=="-h"||a=="--help"){ std::cerr<<"usage: diag_curcuma_atomic_c6 [--tol REL] [--quiet] <atomic_c6.json> <xyz> [...]\n"; return 0; }
        else pos.push_back(a);
    }
    if (pos.size()<2 || pos.size()%2!=0){ std::cerr<<"usage: diag_curcuma_atomic_c6 [--tol REL] [--quiet] <atomic_c6.json> <xyz> [...]\n"; return 2; }
    int rc=0;
    for (size_t i=0;i<pos.size();i+=2){
        double mr = analyse(pos[i], pos[i+1], verbose);
        if (mr<0.0){ rc=2; continue; }
        if (tol>0.0 && mr>tol){ std::printf("  FAIL: max|rel|=%.3e > tol=%.3e\n", mr, tol); rc=1; }
    }
    return rc;
}
