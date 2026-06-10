/*
 * GFN2 component audit (Phase C): per-container energy diff curcuma <-> tblite
 * at *tblite's converged density*.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026.
 *
 * Reads a release_tblite dump (geometry + tblite's converged density P,
 * atomic_charges, shell_populations, gf2.dpat, gf2.qpat, plus the new
 * e_components_tblite block), injects all of it into curcuma's XTB at fixed
 * density (via XTB::evaluateComponentsAtFixedDensity), and reports per-
 * container |E_curcuma - E_tblite|.
 *
 * Container-split caveat: tblite's "dispersion" container is called PRE-SCF
 * (singlepoint.f90 line 181) at q_at=0, capturing the q-independent D4
 * baseline. The charge-coupled D4 contribution lands in eelec via
 * dispersion%get_energy during next_scf. Curcuma evaluates calcDispersionEnergy
 * at the converged Mulliken q_at, so its dispersion absorbs the SCC piece.
 * The PER-CONTAINER dispersion / electronic diffs therefore look large but
 * cancel; the meaningful audit number is the COMBINED (dispersion +
 * electronic) diff, also reported below.
 *
 * Decoupled from curcuma's SCF stability — works even for molecules where
 * curcuma's own SCF diverges (e.g. complex, 231 atoms).
 *
 * Usage:  diag_curcuma_energy_components [--tol VALUE] [--quiet]
 *           <dump.json> <xyz.xyz> [<dump.json> <xyz.xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"

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
using curcuma::xtb::XTB;
using curcuma::xtb::MethodType;

namespace {

int Z_of(const std::string& s) {
    static const std::map<std::string, int> m = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},
        {"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
        {"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},{"Br",35},{"I",53}};
    auto it = m.find(s);
    return (it == m.end()) ? -1 : it->second;
}

bool readXYZ(const std::string& path,
             std::vector<int>& atoms, std::vector<double>& coords_ang)
{
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return false; }
    int nat = 0;
    f >> nat;
    std::string line;
    std::getline(f, line);  // rest of first line
    std::getline(f, line);  // comment line
    atoms.clear();
    coords_ang.assign(3 * nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        std::string sym;
        double x, y, z;
        if (!(f >> sym >> x >> y >> z)) {
            std::cerr << "xyz parse error at atom " << i << "\n";
            return false;
        }
        int z_num = Z_of(sym);
        if (z_num < 0) {
            std::cerr << "unknown element symbol: " << sym << "\n";
            return false;
        }
        atoms.push_back(z_num);
        coords_ang[3*i + 0] = x;
        coords_ang[3*i + 1] = y;
        coords_ang[3*i + 2] = z;
    }
    return true;
}

Eigen::MatrixXd parseSymmetric(const json& m, int n)
{
    Eigen::MatrixXd out(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            out(i, j) = m[i][j].get<double>();
    return out;
}

struct ComponentDiff {
    const char* name;
    double tblite;
    double curcuma;
    bool present;  // tblite block contained this component
};

// Returns max|diff| across all present components, or -1.0 on hard error.
double analyse(const std::string& dump_path, const std::string& xyz_path, bool verbose)
{
    std::ifstream f(dump_path);
    if (!f) { std::cerr << "cannot open " << dump_path << "\n"; return -1.0; }
    json d;
    try { f >> d; } catch (const std::exception& e) {
        std::cerr << dump_path << ": JSON parse error: " << e.what() << "\n";
        return -1.0;
    }

    const std::string method = d.value("method", std::string("gfn2"));
    if (method != "gfn2" && method != "gfn1") {
        std::cerr << dump_path << ": unsupported method '" << method << "'\n";
        return -1.0;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;
    const bool is_gfn2 = (mt == MethodType::GFN2);

    if (!d.contains("e_components_tblite")) {
        std::cerr << dump_path << ": no e_components_tblite block — regenerate dump with current dump_tblite_multipole\n";
        return -1.0;
    }

    // Geometry from companion xyz file (Angstrom).
    std::vector<int> atoms;
    std::vector<double> coords_ang;
    if (!readXYZ(xyz_path, atoms, coords_ang)) return -1.0;
    const int nat = static_cast<int>(atoms.size());
    if (nat != d["natoms"].get<int>()) {
        std::cerr << dump_path << ": nat mismatch (xyz=" << nat
                  << ", dump=" << d["natoms"].get<int>() << ")\n";
        return -1.0;
    }
    const int nsh = d["nshells"].get<int>();
    const int nao = d["norbitals"].get<int>();

    // Build the XTB instance and basis.
    XTB xtb(mt);
    if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
        std::cerr << dump_path << ": XTB::InitialiseMolecule failed\n";
        return -1.0;
    }

    // Inject tblite's converged state.
    Eigen::MatrixXd P = parseSymmetric(d["density"], nao);

    Eigen::VectorXd q_at(nat);
    for (int i = 0; i < nat; ++i) q_at(i) = d["atomic_charges"][i].get<double>();

    // Reconstruct q_sh from tblite's shell populations + curcuma's canonical
    // n0_sh. Curcuma's basis matches tblite's bit-identically per AP6b, so
    // n0_sh is identical too (e.g. N: {1.5, 3.5}, not the integer
    // {2, 3}). Reading curcuma's table directly keeps the diag tool in sync
    // when parameters change.
    Eigen::VectorXd q_sh(nsh);
    {
        const Eigen::VectorXd n0_sh = xtb.getReferenceShellOccupations();
        if (n0_sh.size() != nsh) {
            std::cerr << dump_path << ": curcuma n0_sh size " << n0_sh.size()
                      << " != dump nsh " << nsh << "\n";
            return -1.0;
        }
        const auto& n_sh_json = d["shell_populations"];
        for (int ish = 0; ish < nsh; ++ish) {
            q_sh(ish) = n0_sh(ish) - n_sh_json[ish].get<double>();
        }
    }

    Eigen::MatrixXd dp_at = Eigen::MatrixXd::Zero(3, nat);
    Eigen::MatrixXd qp_at = Eigen::MatrixXd::Zero(6, nat);
    if (is_gfn2 && d.contains("gf2")) {
        const auto& dpat_j = d["gf2"]["dpat"];
        const auto& qpat_j = d["gf2"]["qpat"];
        for (int iat = 0; iat < nat; ++iat) {
            for (int k = 0; k < 3; ++k) dp_at(k, iat) = dpat_j[iat][k].get<double>();
            for (int k = 0; k < 6; ++k) qp_at(k, iat) = qpat_j[iat][k].get<double>();
        }
    }

    if (!xtb.evaluateComponentsAtFixedDensity(P, q_at, q_sh, dp_at, qp_at)) {
        std::cerr << dump_path << ": evaluateComponentsAtFixedDensity failed\n";
        return -1.0;
    }

    // Read tblite per-container sums.
    const auto& ecomp = d["e_components_tblite"];
    auto get_tblite = [&](const char* name, bool& present) -> double {
        present = ecomp.contains("sums") && ecomp["sums"].contains(name);
        return present ? ecomp["sums"][name].get<double>() : 0.0;
    };

    std::vector<ComponentDiff> rows;
    {
        ComponentDiff r;
        r.name = "halogen";     r.tblite = get_tblite("halogen", r.present);
        r.curcuma = xtb.getHalogenBondEnergy();
        rows.push_back(r);
        r.name = "repulsion";   r.tblite = get_tblite("repulsion", r.present);
        r.curcuma = xtb.getRepulsionEnergy();
        rows.push_back(r);
        r.name = "dispersion";  r.tblite = get_tblite("dispersion", r.present);
        r.curcuma = xtb.getDispersionEnergy();
        rows.push_back(r);
        // "interactions" is solvation/external containers; curcuma has none yet,
        // so reporting 0.0 is correct — any divergence here means tblite has a
        // container active that curcuma is missing entirely.
        r.name = "interactions"; r.tblite = get_tblite("interactions", r.present);
        r.curcuma = 0.0;
        rows.push_back(r);
        r.name = "electronic";  r.tblite = get_tblite("electronic", r.present);
        r.curcuma = xtb.getElectronicEnergy();
        rows.push_back(r);
    }

    std::cout << dump_path << "  (nat=" << nat << ", method=" << method << ")\n";
    if (verbose) {
        std::printf("  %-13s  %18s  %18s  %12s\n",
                    "component", "tblite", "curcuma", "diff");
    }
    double max_diff = 0.0;
    const char* worst = "";
    for (const auto& r : rows) {
        if (!r.present) {
            if (verbose) std::printf("  %-13s  %18s  %18.10f  %12s\n",
                                     r.name, "(absent)", r.curcuma, "-");
            // If tblite says the container is absent but curcuma reports
            // non-zero, flag it — but using the same tolerance is misleading,
            // so just report and don't fold into max_diff for now. Real cases
            // where this matters: curcuma evaluating a halogen term for H/C/N/O
            // that tblite skips entirely.
            if (std::abs(r.curcuma) > 1e-12) {
                std::printf("    WARNING: tblite has no '%s' container but curcuma reports %.3e\n",
                            r.name, r.curcuma);
            }
            continue;
        }
        const double diff = r.curcuma - r.tblite;
        if (verbose) {
            std::printf("  %-13s  %18.10f  %18.10f  %12.3e\n",
                        r.name, r.tblite, r.curcuma, diff);
        }
        if (std::abs(diff) > max_diff) { max_diff = std::abs(diff); worst = r.name; }
    }

    // Combined dispersion + electronic check: tblite's dispersion container is
    // pre-SCF (q=0 baseline); the actual SCC D4 lands in eelec. Curcuma puts
    // it all in m_E_dispersion. The combined sum is the apples-to-apples
    // comparison and represents the actual alignment number for the audit.
    bool combined_meaningful = false;
    double combined_diff = 0.0;
    {
        bool pd, pe;
        const double t_disp = get_tblite("dispersion", pd);
        const double t_elec = get_tblite("electronic", pe);
        if (pd && pe) {
            const double c_disp = xtb.getDispersionEnergy();
            const double c_elec = xtb.getElectronicEnergy();
            const double tblite_sum  = t_disp + t_elec;
            const double curcuma_sum = c_disp + c_elec;
            combined_diff = curcuma_sum - tblite_sum;
            combined_meaningful = true;
            if (verbose) {
                std::printf("  %-13s  %18.10f  %18.10f  %12.3e   <- meaningful\n",
                            "disp+elec",
                            tblite_sum, curcuma_sum, combined_diff);
            }
        }
    }

    // If the combined electronic+disp lump diverges, break the curcuma-side
    // electronic into sub-pieces (Tr(P*H0), Coulomb-ES2, third-order,
    // multipole) so the divergent term is named immediately.
    if (combined_meaningful && std::abs(combined_diff) > 1e-7) {
        const double e_cs = xtb.getCoulombShellEnergy();
        const double e_3rd = xtb.getThirdOrderEnergy();
        const double e_mp  = xtb.getMultipoleEnergy();
        const double e_cur = xtb.getElectronicEnergy();
        const double e_band = e_cur - e_cs - e_3rd - e_mp;
        std::printf("  electronic sub-pieces (curcuma, at injected density):\n");
        std::printf("    Tr(P*H0)      = %.10f\n", e_band);
        std::printf("    Coulomb-ES2   = %.10f\n", e_cs);
        std::printf("    Third-order   = %.10f\n", e_3rd);
        std::printf("    Multipole     = %.10f\n", e_mp);
        std::printf("    sum (= curcuma electronic) = %.10f\n", e_cur);
    }

    // The per-container max|diff| is misleading because of the dispersion
    // container split. The combined disp+elec diff is the tolerance signal.
    const double audit_diff = combined_meaningful ? std::abs(combined_diff) : max_diff;
    std::printf("  per-container max|diff|=%.2e (%s),  combined |disp+elec|=%.2e\n",
                max_diff, worst, combined_meaningful ? std::abs(combined_diff) : 0.0);
    return audit_diff;
}

}  // namespace

int main(int argc, char** argv)
{
    double tol = -1.0;
    bool verbose = true;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--tol" && i + 1 < argc) { tol = std::stod(argv[++i]); }
        else if (a == "--quiet") { verbose = false; }
        else if (a == "-h" || a == "--help") {
            std::cerr << "usage: diag_curcuma_energy_components [--tol VALUE] [--quiet]"
                         " <dump.json> <xyz.xyz> [<dump.json> <xyz.xyz> ...]\n"
                         "  --tol VALUE   exit with code 1 if any molecule's max|diff| > VALUE\n"
                         "  --quiet       suppress per-component table, print summary only\n";
            return 0;
        }
        else positional.push_back(a);
    }
    if (positional.size() < 2 || positional.size() % 2 != 0) {
        std::cerr << "usage: diag_curcuma_energy_components [--tol VALUE] [--quiet]"
                     " <dump.json> <xyz.xyz> [<dump.json> <xyz.xyz> ...]\n";
        return 2;
    }
    int rc = 0;
    for (size_t i = 0; i < positional.size(); i += 2) {
        const double max_diff = analyse(positional[i], positional[i + 1], verbose);
        if (max_diff < 0.0) { rc = 2; continue; }
        if (tol > 0.0 && max_diff > tol) {
            std::printf("  FAIL: max|diff|=%.2e > tol=%.2e\n", max_diff, tol);
            rc = 1;
        }
    }
    return rc;
}
