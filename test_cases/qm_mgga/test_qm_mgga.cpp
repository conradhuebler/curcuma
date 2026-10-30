/*
 * Native QM meta-GGA test driver (WP11 stage 1: TPSS, SCAN; energies on the CPU).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout; check_qm_mgga.py compares with reference.json from scripts/qm_mgga_pyscf_reference.py):
 *
 *   test_qm_mgga kernel <reference.json>
 *       TPSS-X, TPSS-C, SCAN-X, SCAN-C (qmxc::evaluateKernelMGGA) at the (rho, sigma, tau) points of reference.json,
 *       sigma raised to libxc's floor (as libxc evaluates): eps, v_rho, v_sigma, v_tau, and central 4-point finite
 *       differences of rho * eps in rho, sigma and tau (step 1e-4 times the variable or smaller, so that alpha moves
 *       by at most 1e-5 max(alpha, 0.01), and 10 and 100 times that step, the closest to the analytic value kept;
 *       null where not taken: sigma = 0, points outside the admissible set, with
 *       alpha < 0.01 (a step crosses tau = tau_W), tau < 1e-18 (libxc's tau floor 1e-20) or rho < 1e-10).
 *
 *   test_qm_mgga energy <xyz> [--functional tpss|scan] [--grid sg1] [--threads N] [--ri off|j|jk]
 *                      [--ao-threshold 1e-12] [--gradient]
 *       Closed-shell KS/def2-SVP SCF (scf_threshold 1e-10, energy 1e-12, orbital gradient 1e-9): E, E_xc, E_x, E_c,
 *       electrons on the grid, iterations; with --gradient whether the request failed (WP11 stage 2: it must not).
 *
 *   test_qm_mgga update <xyzA> <xyzB> [options as energy]
 *       Engine on geometry A, then UpdateMolecule to B, against a fresh engine on B.
 *
 *   test_qm_mgga vxcfd <xyz> [--functional tpss|scan] [--grid sg1] [--ao-threshold t]
 *       V_xc against a finite difference of E_xc in the density matrix: converged P, fixed symmetric perturbation D,
 *       (E_xc(P + hD) - E_xc(P - hD)) / 2h (4-point, h = 1e-4) against Tr(V_xc D); V_xc with 1 and 4 threads
 *       (max |difference|, must be 0). With --ao-threshold > 0 the blocks and kept shells of the AO screen.
 *
 * Claude Generated (Oct 2026, WP11 stage 1). GPL-3.0.
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
    std::string mode, xyz, xyz2, functional = "tpss", grid = "sg1", ri = "off";
    int threads = 4;
    bool gradient = false;
    double ao_threshold = 1e-12;
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
        { "scf_gradient_threshold", 1e-9 }, { "scf_max_iterations", 300 }, { "threads", a.threads },
        { "grid", a.grid }, { "ri", a.ri }, { "xc_ao_threshold", a.ao_threshold } };
}

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
    const std::vector<double> tau = k["tau"].get<std::vector<double>>();
    const std::vector<double> alpha = k["alpha"].get<std::vector<double>>();
    const std::vector<std::string> kind = k["kind"].get<std::vector<std::string>>();
    const double floor = k["sigma_floor"].get<double>();
    const std::pair<const char*, qmxc::XCKernel> kernels[] = { { "tpss_x", qmxc::XCKernel::TPSS_X },
        { "tpss_c", qmxc::XCKernel::TPSS_C }, { "scan_x", qmxc::XCKernel::SCAN_X }, { "scan_c", qmxc::XCKernel::SCAN_C } };
    json out = json::object();
    for (const auto& kk : kernels) {
        const qmxc::XCKernel kern = kk.second;
        auto erho = [&](double r, double s, double t) {
            double e, a, b, c;
            qmxc::evaluateKernelMGGA(kern, r, s, t, e, a, b, c);
            return r * e;
        };
        std::vector<double> eps, vr, vs, vt;
        json fr = json::array(), fs = json::array(), ft = json::array();
        for (size_t i = 0; i < rho.size(); ++i) {
            const double r = rho[i], s = std::max(sigma[i], floor), t = tau[i];
            double e, a, b, c;
            qmxc::evaluateKernelMGGA(kern, r, s, t, e, a, b, c);
            eps.push_back(e);
            vr.push_back(a);
            vs.push_back(b);
            vt.push_back(c);
            // FD only where it is meaningful: admissible points with alpha >= 0.01, tau well above libxc's tau floor
            // 1e-20 (below it the kernel is constant in tau), rho >= 1e-10 (at rho = 1e-12 the cancellation in
            // eps_c amplifies rounding so that no step gives 1e-6: measured 4e-5 at the best step, SCAN-C, s = 30)
            const bool fd_ok = kind[i] == "admissible" && alpha[i] >= 0.01 && tau[i] >= 1e-18 && r >= 1e-10;
            if (fd_ok) {
                // steps: relative 1e-4, but small enough that alpha moves by at most 1e-5 max(alpha, 0.01) over one
                // step (alpha = (tau - sigma/(8 rho)) / tau_unif is very sensitive to tau and sigma at large s, and
                // the stencil must not cross tau = tau_W or a SCAN cut-off)
                const double tu = 0.3 * std::cbrt(9.0 * M_PI * M_PI * M_PI * M_PI) * r * std::cbrt(r * r);
                const double tw = s / (8.0 * r);
                const double da = 1e-5 * std::max(alpha[i], 0.01);
                const double hr = std::min(1e-4 * r, da * r / (tw / tu + 5.0 / 3.0 * alpha[i]));
                const double ht = std::min(1e-4 * t, da * tu);
                const double hs = std::min(1e-4 * s, da * 8.0 * r * tu);
                // three steps h, 10 h, 100 h (no larger than 1e-3 x): the one closest to the analytic value is kept
                // (small steps are noise-limited where eps_c cancels at large t, large ones may cross a branch)
                auto best = [](auto fn, double x, double h, double ref) {
                    double out = fd4(fn, x, h);
                    for (double m : { 10.0, 100.0 }) {
                        if (m * h > 1e-3 * x) break;
                        const double v = fd4(fn, x, m * h);
                        if (std::abs(v - ref) < std::abs(out - ref)) out = v;
                    }
                    return out;
                };
                fr.push_back(best([&](double x) { return erho(x, s, t); }, r, hr, a));
                ft.push_back(best([&](double x) { return erho(r, s, x); }, t, ht, c));
                if (sigma[i] > 0.0)
                    fs.push_back(best([&](double x) { return erho(r, x, t); }, s, hs, b));
                else
                    fs.push_back(nullptr);
            } else {
                fr.push_back(nullptr);
                fs.push_back(nullptr);
                ft.push_back(nullptr);
            }
        }
        out[kk.first] = { { "eps", eps }, { "vrho", vr }, { "vsigma", vs }, { "vtau", vt }, { "vrho_fd", fr },
            { "vsigma_fd", fs }, { "vtau_fd", ft } };
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

QMFunctional functionalOf(const std::string& s)
{
    return s == "scan" ? QMFunctional::SCAN : s == "pbe" ? QMFunctional::PBE : QMFunctional::TPSS;
}

std::unique_ptr<QMEngine> makeEngine(const Args& a, const std::string& xyz)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(xyz, Z, ang)) {
        std::cerr << "cannot read " << xyz << "\n";
        return nullptr;
    }
    auto e = std::make_unique<QMEngine>(functionalOf(a.functional), engineConfig(a));
    if (!e->QMInterface::InitialiseMolecule(Z.data(), ang.data(), static_cast<int>(Z.size()), 0.0, 0)) return nullptr;
    return e;
}

int modeEnergy(const Args& a)
{
    auto e = makeEngine(a, a.xyz);
    if (!e) return 1;
    const double E = e->Calculation(a.gradient);
    json out = { { "energy", E }, { "converged", e->scfConverged() }, { "iterations", e->scfIterations() },
        { "exc", e->xcEnergy() }, { "ex", e->xcExchangeEnergy() }, { "ec", e->xcCorrelationEnergy() },
        { "nelec_grid", e->xcGridElectrons() }, { "num_electrons", e->numElectrons() },
        { "label", e->xcFunctional().label() }, { "has_gradient", e->hasGradient() } };
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
    qmgrid::BlockScreen screen;
    if (a.ao_threshold > 0.0) {
        std::vector<int> begin, count;
        qmgrid::sortGridBlocks(grid, 128, begin, count);
        if (!qmgrid::buildBlockScreen(grid, ao, a.ao_threshold, begin, count, screen, why)) {
            std::cerr << why << "\n";
            return 1;
        }
    }
    const Matrix P = e->density();
    const int n = static_cast<int>(P.rows());
    Matrix D(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) D(i, j) = 0.05 * (std::sin(1.3 * (i + 1) * (j + 1)) + std::cos(0.7 * (i + j + 2)));
    D = Matrix(0.5 * (D + D.transpose()));
    const qmxc::XCFunctional f = e->xcFunctional();
    const double thr = 1e-12;
    auto exc = [&](double h) {
        qmxc::XCResult r;
        qmxc::buildVxc(grid, ao, Matrix(P + h * D), f, thr, 4, r, 128, &screen);
        return r.exc;
    };
    qmxc::XCResult r4, r1;
    qmxc::buildVxc(grid, ao, P, f, thr, 4, r4, 128, &screen);
    qmxc::buildVxc(grid, ao, P, f, thr, 1, r1, 128, &screen);
    const double analytic = r4.V.cwiseProduct(D).sum();
    const double fd = fd4(exc, 0.0, 1e-4);
    json out = { { "trVD", analytic }, { "fd", fd }, { "rel", std::abs(analytic - fd) / std::abs(analytic) },
        { "threads_maxdiff", (r4.V - r1.V).cwiseAbs().maxCoeff() }, { "exc_threads_diff", r4.exc - r1.exc },
        { "asym", (r4.V - r4.V.transpose()).cwiseAbs().maxCoeff() }, { "npts", grid.size() }, { "exc", r4.exc },
        { "screened", screen.active() }, { "label", f.label() } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 3) {
        std::cerr << "usage: test_qm_mgga kernel <reference.json> | energy <xyz> | update <xyzA> <xyzB> | vxcfd <xyz> "
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
        else if (s == "--ri") a.ri = next();
        else if (s == "--threads") a.threads = std::stoi(next());
        else if (s == "--ao-threshold") a.ao_threshold = std::stod(next());
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
