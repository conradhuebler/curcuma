/*
 * GFN2-D4 reference-C6 audit: curcuma's m_c6_flat_cache vs dftd4's model%c6.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026. GPL-3.0.
 *
 * Reads a c6 dump produced by dump_dftd4_c6 (dftd4 new_d4_model ref=gfn2,
 * model%c6(iref,jref,isp,jsp) — geometry/charge-INDEPENDENT reference C6),
 * builds curcuma's GFN2 D4 reference cache the exact way
 * xtb_native::calcDispersionEnergy does (setD4CovalentCN(true), per-reference),
 * and diffs per (Z_i, Z_j, iref, jref).
 *
 * This isolates the model-construction step (zeta-scaled aiw -> Casimir-Polder
 * C6) from the CN/charge weighting. The prior audit (docs/GFN2_D4_STATUS.md)
 * proved all *input* tables bit-identical; the remaining ~5% C-path residual
 * was suspected to live in the constructed reference C6 itself. This tool
 * confirms or refutes that directly.
 *
 *   diag_curcuma_d4_c6 [--tol REL] [--quiet] <c6_dump.json> <xyz.xyz> [...]
 *     --tol REL   exit 1 if any pair's max relative |diff| > REL
 *     --quiet     summary only (suppress per-entry table)
 */

#include "src/core/config_manager.h"
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

// Must mirror D4ParameterGenerator's private layout (d4param_generator.h).
static constexpr int CC_MAX_ELEM = 118;
static constexpr int CC_MAX_REF = 7;
static inline size_t curcumaC6Index(int elem_i, int elem_j, int ref_i, int ref_j) {
    return static_cast<size_t>(elem_i) * CC_MAX_ELEM * CC_MAX_REF * CC_MAX_REF
         + static_cast<size_t>(elem_j) * CC_MAX_REF * CC_MAX_REF
         + static_cast<size_t>(ref_i) * CC_MAX_REF
         + static_cast<size_t>(ref_j);
}

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
    int nat = 0;
    f >> nat;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    atoms.clear();
    geom_ang = Matrix::Zero(nat, 3);
    for (int i = 0; i < nat; ++i) {
        std::string sym; double x, y, z;
        if (!(f >> sym >> x >> y >> z)) { std::cerr << "xyz parse error at atom " << i << "\n"; return false; }
        atoms.push_back(Z_of(sym));
        geom_ang(i, 0) = x; geom_ang(i, 1) = y; geom_ang(i, 2) = z;
    }
    return true;
}

// Returns max relative |diff| over all compared (Zi,Zj,iref,jref), or -1.0 on error.
static double analyse(const std::string& c6_path, const std::string& xyz_path, bool verbose) {
    std::ifstream f(c6_path);
    if (!f) { std::cerr << "cannot open " << c6_path << "\n"; return -1.0; }
    json d;
    try { f >> d; } catch (const std::exception& e) {
        std::cerr << c6_path << ": JSON parse error: " << e.what() << "\n";
        return -1.0;
    }

    const int nid  = d.at("nid").get<int>();
    const int mref = d.at("mref").get<int>();
    std::vector<int> species_z = d.at("species_z").get<std::vector<int>>();
    std::vector<int> nref      = d.at("nref").get<std::vector<int>>();
    std::vector<double> c6_tbl = d.at("c6").get<std::vector<double>>();
    if (static_cast<int>(species_z.size()) != nid || static_cast<int>(nref.size()) != nid) {
        std::cerr << c6_path << ": species_z/nref size mismatch\n"; return -1.0;
    }
    if (static_cast<int>(c6_tbl.size()) != nid * nid * mref * mref) {
        std::cerr << c6_path << ": c6 flat size mismatch\n"; return -1.0;
    }
    auto tblIdx = [&](int isp, int jsp, int iref, int jref) {  // all 0-based
        return ((static_cast<size_t>(isp) * nid + jsp) * mref + iref) * mref + jref;
    };

    // Geometry only needed so the generator can run GenerateParameters; the
    // reference C6 cache it fills is geometry/charge independent.
    std::vector<int> atoms;
    Matrix geom_ang;
    if (!readXYZ(xyz_path, atoms, geom_ang)) return -1.0;
    const int nat = static_cast<int>(atoms.size());
    const Matrix geom_bohr = geom_ang * AA_TO_AU;

    // Mirror xtb_native::calcDispersionEnergy GFN2 setup (same as diag_curcuma_d4_potential).
    json cfg;
    cfg["d4_s6"] = 1.0; cfg["d4_s8"] = 2.7; cfg["d4_a1"] = 0.52;
    cfg["d4_a2"] = 5.0; cfg["d4_alp"] = 16.0;
    ConfigManager cm("d4param", cfg);
    D4ParameterGenerator gen(cm);
    gen.setUseD4SingleShotEEQ(false);
    gen.setD4CovalentCN(true);                       // GFN2 alpha-zeta scaling on
    gen.setTopologyCharges(Vector::Zero(nat));       // c6 ref is charge-independent
    gen.GenerateParameters(atoms, geom_bohr);

    const std::vector<double>& c6_cur = gen.getC6FlatCache();
    const std::vector<int>& refn_cur = gen.getRefN();   // 0-based by Z-1

    std::cout << c6_path << "  (nid=" << nid << ", mref=" << mref << ")\n";
    if (verbose) {
        std::printf("  %3s %3s %4s %4s  %20s %20s %12s %10s\n",
                    "Zi", "Zj", "ir", "jr", "tblite_c6", "curcuma_c6", "abs_diff", "rel");
    }

    double max_rel = 0.0, max_abs = 0.0;
    int worst_zi = -1, worst_zj = -1, worst_ir = -1, worst_jr = -1;
    bool nref_mismatch = false;

    for (int isp = 0; isp < nid; ++isp) {
        const int Zi = species_z[isp];
        const int cur_nref_i = (Zi - 1 < static_cast<int>(refn_cur.size())) ? refn_cur[Zi - 1] : 0;
        if (cur_nref_i != nref[isp]) {
            std::printf("  WARNING: nref mismatch for Z=%d: tblite=%d curcuma=%d\n",
                        Zi, nref[isp], cur_nref_i);
            nref_mismatch = true;
        }
        for (int jsp = 0; jsp < nid; ++jsp) {
            const int Zj = species_z[jsp];
            const int cur_nref_j = (Zj - 1 < static_cast<int>(refn_cur.size())) ? refn_cur[Zj - 1] : 0;
            const int n_i = std::min(nref[isp], cur_nref_i);
            const int n_j = std::min(nref[jsp], cur_nref_j);
            for (int iref = 0; iref < n_i && iref < mref; ++iref) {
                for (int jref = 0; jref < n_j && jref < mref; ++jref) {
                    const double t = c6_tbl[tblIdx(isp, jsp, iref, jref)];
                    const double c = c6_cur[curcumaC6Index(Zi - 1, Zj - 1, iref, jref)];
                    const double ad = std::abs(c - t);
                    const double rel = (std::abs(t) > 1e-12) ? ad / std::abs(t) : ad;
                    if (rel > max_rel) {
                        max_rel = rel; worst_zi = Zi; worst_zj = Zj; worst_ir = iref; worst_jr = jref;
                    }
                    if (ad > max_abs) max_abs = ad;
                    if (verbose) {
                        std::printf("  %3d %3d %4d %4d  %20.12e %20.12e %12.3e %10.3e\n",
                                    Zi, Zj, iref, jref, t, c, c - t, rel);
                    }
                }
            }
        }
    }

    std::printf("  max|rel|=%.3e  max|abs|=%.3e  worst (Zi=%d Zj=%d ir=%d jr=%d)%s\n",
                max_rel, max_abs, worst_zi, worst_zj, worst_ir, worst_jr,
                nref_mismatch ? "   [NREF MISMATCH]" : "");
    return max_rel;
}

int main(int argc, char** argv) {
    double tol = -1.0;
    bool verbose = true;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--tol" && i + 1 < argc) { tol = std::stod(argv[++i]); }
        else if (a == "--quiet") { verbose = false; }
        else if (a == "-h" || a == "--help") {
            std::cerr << "usage: diag_curcuma_d4_c6 [--tol REL] [--quiet]"
                         " <c6_dump.json> <xyz.xyz> [<c6_dump.json> <xyz.xyz> ...]\n";
            return 0;
        }
        else positional.push_back(a);
    }
    if (positional.size() < 2 || positional.size() % 2 != 0) {
        std::cerr << "usage: diag_curcuma_d4_c6 [--tol REL] [--quiet]"
                     " <c6_dump.json> <xyz.xyz> [<c6_dump.json> <xyz.xyz> ...]\n";
        return 2;
    }
    int rc = 0;
    for (size_t i = 0; i < positional.size(); i += 2) {
        const double max_rel = analyse(positional[i], positional[i + 1], verbose);
        if (max_rel < 0.0) { rc = 2; continue; }
        if (tol > 0.0 && max_rel > tol) {
            std::printf("  FAIL: max|rel|=%.3e > tol=%.3e\n", max_rel, tol);
            rc = 1;
        }
    }
    return rc;
}
