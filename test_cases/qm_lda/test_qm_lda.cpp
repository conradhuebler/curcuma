/*
 * Native QM LDA test driver (WP5).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout; check_qm_lda.py compares with reference.json from scripts/qm_lda_pyscf_reference.py):
 *
 *   test_qm_lda kernel <reference.json>
 *       Slater exchange and VWN5 / VWN-RPA correlation (qmxc point kernels) at the densities of reference.json
 *       ("kernel"."rho"): eps, v, and v from a central 4-point finite difference of rho * eps (h = 1e-3 rho).
 *
 *   test_qm_lda energy <xyz> [--grid sg1] [--vwn vwn5] [--threads N] [--ri off|j|jk] [--direct auto|on|off]
 *                     [--cartesian] [--gradient] [--threshold 1e-12]
 *       Closed-shell LDA/def2-SVP SCF (scf_threshold 1e-10, energy 1e-12, orbital gradient 1e-9): E, E_xc, E_x, E_c,
 *       electrons on the grid, iterations; with --gradient also whether the gradient request failed (since WP8 DFT it must not).
 *
 *   test_qm_lda update <xyzA> <xyzB> [options as energy]
 *       Engine on geometry A, then UpdateMolecule to B (same atoms), against a fresh engine on B: the grid and the
 *       V_xc cache must follow the geometry.
 *
 * Claude Generated (Oct 2026, WP5). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/energy_calculators/qm_methods/qm_xc.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string mode, xyz, xyz2, grid = "sg1", vwn = "vwn5", ri = "off", direct = "auto";
    int threads = 4;
    bool cartesian = false, gradient = false;
    double threshold = 1e-12;
};

bool readXYZ(const std::string& path, std::vector<int>& Z, std::vector<double>& ang)
{
    static const std::map<std::string, int> el{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
        { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 } };
    std::ifstream f(path);
    int n = 0;
    if (!(f >> n) || n <= 0) return false;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    for (int i = 0; i < n; ++i) {
        std::string s;
        double x, y, z;
        if (!(f >> s >> x >> y >> z)) return false;
        auto it = el.find(s);
        if (it == el.end()) return false;
        Z.push_back(it->second);
        ang.insert(ang.end(), { x, y, z });
    }
    return true;
}

json engineConfig(const Args& a)
{
    return { { "basis", "def2-SVP" }, { "scf_threshold", 1e-10 }, { "scf_energy_threshold", 1e-12 },
        { "scf_gradient_threshold", 1e-9 }, { "scf_max_iterations", 200 }, { "threads", a.threads },
        { "grid", a.grid }, { "vwn_variant", a.vwn }, { "ri", a.ri }, { "scf_direct", a.direct },
        { "cartesian_d", a.cartesian }, { "xc_density_threshold", a.threshold } };
}

json resultJson(QMEngine& e, double E)
{
    return { { "energy", E }, { "converged", e.scfConverged() }, { "iterations", e.scfIterations() },
        { "exc", e.xcEnergy() }, { "ex", e.xcExchangeEnergy() }, { "ec", e.xcCorrelationEnergy() },
        { "nelec_grid", e.xcGridElectrons() }, { "num_electrons", e.numElectrons() }, { "label", e.xcFunctional().label() } };
}

int modeKernel(const std::string& ref_path)
{
    std::ifstream f(ref_path);
    json ref;
    try {
        f >> ref;
    } catch (const std::exception& ex) {
        std::cerr << "cannot read " << ref_path << ": " << ex.what() << "\n";
        return 1;
    }
    const std::vector<double> rho = ref["kernel"]["rho"].get<std::vector<double>>();
    json out = json::object();
    auto run = [&](const std::string& key, auto fn) {
        std::vector<double> eps, v, vfd;
        for (double r : rho) {
            double e, vv;
            fn(r, e, vv);
            eps.push_back(e);
            v.push_back(vv);
            // d(rho eps)/d rho by a central 4-point difference, h = 1e-3 rho
            const double h = 1e-3 * r;
            auto g = [&](double x) {
                double ee, vx;
                fn(x, ee, vx);
                return x * ee;
            };
            vfd.push_back((-g(r + 2 * h) + 8 * g(r + h) - 8 * g(r - h) + g(r - 2 * h)) / (12 * h));
        }
        out[key] = { { "eps", eps }, { "v", v }, { "v_fd", vfd } };
    };
    run("slater", [](double r, double& e, double& v) { qmxc::slaterExchange(r, e, v); });
    run("vwn5", [](double r, double& e, double& v) { qmxc::vwnCorrelation(r, qmxc::VWNVariant::VWN5, e, v); });
    run("vwn_rpa", [](double r, double& e, double& v) { qmxc::vwnCorrelation(r, qmxc::VWNVariant::VWN_RPA, e, v); });
    out["rho"] = rho;
    std::cout << out.dump(1) << std::endl;
    return 0;
}

std::unique_ptr<QMEngine> makeEngine(const Args& a, const std::string& xyz)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(xyz, Z, ang)) {
        std::cerr << "cannot read " << xyz << "\n";
        return nullptr;
    }
    auto e = std::make_unique<QMEngine>(QMFunctional::LDA, engineConfig(a));
    if (!e->QMInterface::InitialiseMolecule(Z.data(), ang.data(), static_cast<int>(Z.size()), 0.0, 0)) return nullptr;
    return e;
}

int modeEnergy(const Args& a)
{
    auto e = makeEngine(a, a.xyz);
    if (!e) return 1;
    const double E = e->Calculation(a.gradient);
    json out = resultJson(*e, E);
    out["has_gradient"] = e->hasGradient();
    if (a.gradient) {
        out["gradient_failed"] = e->gradientFailed();
        out["gradient_error"] = e->gradientError();
    }
    std::cout << out.dump(1) << std::endl;
    return e->scfConverged() ? 0 : 2;
}

int modeUpdate(const Args& a)
{
    auto e = makeEngine(a, a.xyz);
    auto fresh = makeEngine(a, a.xyz2);
    if (!e || !fresh) return 1;
    std::vector<int> Z;
    std::vector<double> ang;
    readXYZ(a.xyz2, Z, ang);
    const double EA = e->Calculation(false);
    Matrix geo(static_cast<int>(Z.size()), 3);
    for (int i = 0; i < static_cast<int>(Z.size()); ++i)
        for (int k = 0; k < 3; ++k) geo(i, k) = ang[3 * i + k];
    e->QMInterface::UpdateMolecule(geo);
    const double EB_updated = e->Calculation(false);
    const double EB_fresh = fresh->Calculation(false);
    json out = { { "energy_A", EA }, { "energy_B_updated", EB_updated }, { "energy_B_fresh", EB_fresh },
        { "exc_B_updated", e->xcEnergy() }, { "exc_B_fresh", fresh->xcEnergy() },
        { "warm_start", e->lastScfWarmStarted() }, { "converged", e->scfConverged() && fresh->scfConverged() } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 3) {
        std::cerr << "usage: test_qm_lda kernel <reference.json> | energy <xyz> | update <xyzA> <xyzB> [options]\n";
        return 1;
    }
    Args a;
    a.mode = argv[1];
    a.xyz = argv[2];
    int i = 3;
    if (a.mode == "update") {
        if (argc < 4) return 1;
        a.xyz2 = argv[3];
        i = 4;
    }
    for (; i < argc; ++i) {
        const std::string s = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if (s == "--grid") a.grid = next();
        else if (s == "--vwn") a.vwn = next();
        else if (s == "--ri") a.ri = next();
        else if (s == "--direct") a.direct = next();
        else if (s == "--threads") a.threads = std::stoi(next());
        else if (s == "--threshold") a.threshold = std::stod(next());
        else if (s == "--cartesian") a.cartesian = true;
        else if (s == "--gradient") a.gradient = true;
        else {
            std::cerr << "unknown option " << s << "\n";
            return 1;
        }
    }
    try {
        if (a.mode == "kernel") return modeKernel(a.xyz);
        if (a.mode == "energy") return modeEnergy(a);
        if (a.mode == "update") return modeUpdate(a);
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 3;
    }
    std::cerr << "unknown mode " << a.mode << "\n";
    return 1;
}
