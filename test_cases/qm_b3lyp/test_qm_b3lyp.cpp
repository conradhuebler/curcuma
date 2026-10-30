/*
 * Native QM B3LYP test driver (WP7), derived from test_qm_pbe.cpp (WP6).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout; check_qm_b3lyp.py compares with reference.json from scripts/qm_b3lyp_pyscf_reference.py):
 *
 *   test_qm_b3lyp kernel <reference.json>
 *       B88 exchange and LYP correlation (qmxc point kernels) at the (rho, sigma) pairs of reference.json:
 *       eps, v_rho, v_sigma at max(sigma, libxc floor) (as libxc evaluates), eps at the raw sigma, and central
 *       4-point finite differences of rho * eps in rho (h = 1e-3 rho) and in sigma (h = 1e-3 sigma; null at sigma = 0).
 *
 *   test_qm_b3lyp energy <xyz> [--functional b3lyp|pbe|lda|hf] [--vwn vwn5|vwn_rpa] [--grid sg1] [--threads N]
 *                       [--ri off|j|jk] [--direct auto|on|off] [--gradient] [--threshold 1e-12]
 *       Closed-shell KS/def2-SVP SCF (scf_threshold 1e-10, energy 1e-12, orbital gradient 1e-9): E, E_xc, E_x, E_c,
 *       electrons on the grid, iterations, exact-exchange energy -(a_x/4) Tr(PK); with --gradient also whether the
 *       gradient request failed (since WP8 DFT it must not).
 *
 *   test_qm_b3lyp update <xyzA> <xyzB> [options as energy]
 *       Engine on geometry A, then UpdateMolecule to B, against a fresh engine on B.
 *
 *   test_qm_b3lyp vxcfd <xyz> [--grid sg1]
 *       V_xc against a finite difference of E_xc: with the converged density P and a fixed symmetric perturbation D,
 *       (E_xc(P + hD) - E_xc(P - hD)) / 2h (4-point, h = 1e-4) against Tr(V_xc D); also V_xc with 1 and 4 threads
 *       (max |difference|, must be 0: fixed chunking).
 *
 * Claude Generated (Oct 2026, WP7). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/energy_calculators/qm_methods/qm_grid.h"
#include "src/core/energy_calculators/qm_methods/qm_xc.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string mode, xyz, xyz2, functional = "b3lyp", vwn = "vwn5", grid = "sg1", ri = "off", direct = "auto";
    int threads = 4;
    bool gradient = false;
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
        { "grid", a.grid }, { "ri", a.ri }, { "scf_direct", a.direct }, { "xc_density_threshold", a.threshold },
        { "vwn_variant", a.vwn } };
}

json resultJson(QMEngine& e, double E)
{
    return { { "energy", E }, { "converged", e.scfConverged() }, { "iterations", e.scfIterations() },
        { "exc", e.xcEnergy() }, { "ex", e.xcExchangeEnergy() }, { "ec", e.xcCorrelationEnergy() },
        { "nelec_grid", e.xcGridElectrons() }, { "num_electrons", e.numElectrons() },
        { "label", e.xcFunctional().label() }, { "exact_exchange", e.energyComponents()[3] },
        { "a_x", e.xcFunctional().exact_exchange } };
}

/// central 4-point difference of f at x with step h
template <class F>
double fd4(F f, double x, double h)
{
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
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
    const json& k = ref["kernel"];
    const std::vector<double> rho = k["rho"].get<std::vector<double>>();
    const std::vector<double> sigma = k["sigma"].get<std::vector<double>>();
    json out = json::object();
    using KFn = void (*)(double, double, double&, double&, double&);
    auto run = [&](const std::string& key, KFn fn, double floor) {
        std::vector<double> eps, vr, vs, eps_raw, vr_fd;
        json vs_fd = json::array();
        for (size_t i = 0; i < rho.size(); ++i) {
            const double r = rho[i], s = std::max(sigma[i], floor);
            double e, a, b;
            fn(r, s, e, a, b);
            eps.push_back(e);
            vr.push_back(a);
            vs.push_back(b);
            fn(r, sigma[i], e, a, b);
            eps_raw.push_back(e);
            vr_fd.push_back(fd4([&](double x) { double ee, aa, bb; fn(x, s, ee, aa, bb); return x * ee; }, r, 1e-3 * r));
            if (sigma[i] > 0.0)
                vs_fd.push_back(fd4([&](double x) { double ee, aa, bb; fn(r, x, ee, aa, bb); return r * ee; }, s,
                                    1e-3 * s));
            else
                vs_fd.push_back(nullptr);
        }
        out[key] = { { "eps", eps }, { "vrho", vr }, { "vsigma", vs }, { "eps_raw", eps_raw }, { "vrho_fd", vr_fd },
            { "vsigma_fd", vs_fd } };
    };
    run("b88", qmxc::b88Exchange, k["sigma_floor"]["b88"].get<double>());
    run("lyp", qmxc::lypCorrelation, k["sigma_floor"]["lyp"].get<double>());
    // NaN/inf scan: both kernels on a dense grid down to the GGA floor (1e-20) and up to sigma 1e20
    long bad = 0, n_scan = 0;
    for (int i = -80; i <= 8; ++i)
        for (int j = -40; j <= 80; j += 2) {
            const double r = std::pow(10.0, 0.25 * i), s = j == -40 ? 0.0 : std::pow(10.0, 0.25 * j);
            for (KFn fn : { qmxc::b88Exchange, qmxc::lypCorrelation }) {
                double e, a, b;
                fn(r, s, e, a, b);
                ++n_scan;
                if (!std::isfinite(e) || !std::isfinite(a) || !std::isfinite(b)) ++bad;
            }
        }
    out["nonfinite"] = bad;
    out["nonfinite_scanned"] = n_scan;
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
    const QMFunctional fn = a.functional == "lda" ? QMFunctional::LDA
        : a.functional == "pbe"                   ? QMFunctional::PBE
        : a.functional == "hf"                    ? QMFunctional::HF
                                                  : QMFunctional::B3LYP;
    auto e = std::make_unique<QMEngine>(fn, engineConfig(a));
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

int modeVxcFD(const Args& a)
{
    auto e = makeEngine(a, a.xyz);
    if (!e) return 1;
    e->Calculation(false);
    if (!e->scfConverged()) return 2;
    qmgrid::MolecularGrid grid;
    std::string why;
    if (!e->buildGrid(grid, why)) {
        std::cerr << why << "\n";
        return 1;
    }
    const qmgrid::AOEvaluator ao(e->gtoBasis(), e->sphericalTransform());
    const Matrix P = e->density();
    const int n = static_cast<int>(P.rows());
    // fixed symmetric perturbation, entries of size ~0.1 times a smooth pattern (no RNG: reproducible everywhere)
    Matrix D(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) D(i, j) = 0.05 * (std::sin(1.3 * (i + 1) * (j + 1)) + std::cos(0.7 * (i + j + 2)));
    D = Matrix(0.5 * (D + D.transpose()));
    const qmxc::XCFunctional f = e->xcFunctional();
    auto exc = [&](double h) {
        qmxc::XCResult r;
        qmxc::buildVxc(grid, ao, Matrix(P + h * D), f, a.threshold, 4, r);
        return r.exc;
    };
    qmxc::XCResult r4, r1;
    qmxc::buildVxc(grid, ao, P, f, a.threshold, 4, r4);
    qmxc::buildVxc(grid, ao, P, f, a.threshold, 1, r1);
    const double analytic = r4.V.cwiseProduct(D).sum();
    const double h = 1e-4;
    const double fd = fd4(exc, 0.0, h);
    json out = { { "trVD", analytic }, { "fd", fd }, { "rel", std::abs(analytic - fd) / std::abs(analytic) },
        { "threads_maxdiff", (r4.V - r1.V).cwiseAbs().maxCoeff() }, { "exc_threads_diff", r4.exc - r1.exc },
        { "asym", (r4.V - r4.V.transpose()).cwiseAbs().maxCoeff() }, { "npts", grid.size() },
        { "label", f.label() } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 3) {
        std::cerr << "usage: test_qm_b3lyp kernel <reference.json> | energy <xyz> | update <xyzA> <xyzB> | vxcfd <xyz> "
                     "[options]\n";
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
        else if (s == "--functional") a.functional = next();
        else if (s == "--vwn") a.vwn = next();
        else if (s == "--ri") a.ri = next();
        else if (s == "--direct") a.direct = next();
        else if (s == "--threads") a.threads = std::stoi(next());
        else if (s == "--threshold") a.threshold = std::stod(next());
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
        if (a.mode == "vxcfd") return modeVxcFD(a);
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 3;
    }
    std::cerr << "unknown mode " << a.mode << "\n";
    return 1;
}
