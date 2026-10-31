/*
 * Native double-hybrid test driver (WP12 stage 1: B2PLYP, B2GP-PLYP, DSD-*, revDSD-PBEP86; energies, CPU).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout; check_qm_dh.py compares with reference.json from scripts/qm_dh_pyscf_reference.py):
 *
 *   test_qm_dh kernel <reference.json>
 *       P86 (PZ81 local), P86 (VWN5 local) at the (rho, sigma) points, B95 at the (rho, sigma, tau) points of
 *       reference.json, sigma raised to libxc's floor: eps, v_rho, v_sigma (v_tau), and 4-point finite differences of
 *       rho * eps (step 1e-4 x, also 10 and 100 times smaller-than-1e-3 x steps; the closest to the analytic value is
 *       kept; only rho >= 1e-10, sigma > 0, B95 admissible points with tau >= 1e-18).
 *
 *   test_qm_dh energy <xyz> --method <name> [--grid sg1] [--threads N] [--fc] [--pt2 ri|canonical]
 *                     [--p86 pz81|vwn5] [--dispersion auto|none|d3bj|d4] [--gradient] [--ao-threshold t]
 *       DoubleHybridMethod (def2-SVP, SCF 1e-10 / dE 1e-12 / orbital gradient 1e-9): getEnergyDecomposition() plus
 *       the error flag and message (a gradient request must be an error, stage 1).
 *
 *   test_qm_dh update <xyzA> <xyzB> --method <name> [options as energy]
 *       Method on geometry A, then updateGeometry(B), against a fresh method on B.
 *
 * Claude Generated (Oct 2026, WP12 stage 1). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/double_hybrid_method.h"
#include "src/core/energy_calculators/qm_methods/qm_xc.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string mode, xyz, xyz2, method = "b2plyp", grid = "sg1", pt2 = "ri", p86 = "pz81", dispersion = "auto";
    int threads = 4;
    bool fc = false, gradient = false;
    double ao_threshold = 1e-12;
};

bool readMol(const std::string& path, Mol& m)
{
    static const std::map<std::string, int> el{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
        { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 } };
    std::ifstream f(path);
    int n = 0;
    if (!(f >> n) || n <= 0) return false;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    m = Mol();
    m.m_number_atoms = n;
    m.m_geometry = Geometry::Zero(n, 3);
    for (int i = 0; i < n; ++i) {
        std::string s;
        double x, y, z;
        if (!(f >> s >> x >> y >> z)) return false;
        auto it = el.find(s);
        if (it == el.end()) return false;
        m.m_atoms.push_back(it->second);
        m.m_geometry(i, 0) = x;
        m.m_geometry(i, 1) = y;
        m.m_geometry(i, 2) = z;
    }
    m.m_charge = 0;
    m.m_spin = 0;
    return true;
}

json methodConfig(const Args& a)
{
    json qm = { { "basis", "def2-SVP" }, { "scf_threshold", 1e-10 }, { "scf_energy_threshold", 1e-12 },
        { "scf_gradient_threshold", 1e-9 }, { "scf_gradient_threshold_tight", 1e-9 }, { "scf_max_iterations", 300 },
        { "threads", a.threads }, { "grid", a.grid }, { "frozen_core", a.fc }, { "dh_pt2", a.pt2 },
        { "p86_lda", a.p86 }, { "dispersion", a.dispersion }, { "xc_ao_threshold", a.ao_threshold } };
    return { { "qm", qm } };
}

template <class F>
double fd4(F f, double x, double h)
{
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
}

template <class F>
double bestFD(F fn, double x, double h, double ref)
{
    double out = fd4(fn, x, h);
    for (double m : { 10.0, 0.1, 0.01 }) {  // 10 h: noise-limited points (cancellation at low density)
        const double v = fd4(fn, x, m * h);
        if (std::abs(v - ref) < std::abs(out - ref)) out = v;
    }
    return out;
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
    json out = json::object();
    {  // GGA: P86 on PZ81 and on VWN5
        const json& g = ref["kernel"]["gga"];
        const std::vector<double> rho = g["rho"].get<std::vector<double>>(), sigma = g["sigma"].get<std::vector<double>>();
        const std::pair<const char*, qmxc::XCKernel> kernels[] = { { "p86", qmxc::XCKernel::P86_C },
            { "p86vwn", qmxc::XCKernel::P86VWN_C } };
        for (const auto& kk : kernels) {
            const double floor = g[kk.first]["sigma_floor"].get<double>();
            auto erho = [&](double r, double s) {
                double e, a, b;
                qmxc::evaluateKernel(kk.second, r, s, e, a, b);
                return r * e;
            };
            std::vector<double> eps, vr, vs;
            json fr = json::array(), fs = json::array();
            for (size_t i = 0; i < rho.size(); ++i) {
                const double r = rho[i], s = std::max(sigma[i], floor);
                double e, a, b;
                qmxc::evaluateKernel(kk.second, r, s, e, a, b);
                eps.push_back(e);
                vr.push_back(a);
                vs.push_back(b);
                if (r >= 1e-10) {
                    fr.push_back(bestFD([&](double x) { return erho(x, s); }, r, 1e-4 * r, a));
                    if (sigma[i] > 0.0)
                        fs.push_back(bestFD([&](double x) { return erho(r, x); }, s, 1e-4 * s, b));
                    else
                        fs.push_back(nullptr);
                } else {
                    fr.push_back(nullptr);
                    fs.push_back(nullptr);
                }
            }
            out[kk.first] = { { "eps", eps }, { "vrho", vr }, { "vsigma", vs }, { "vrho_fd", fr }, { "vsigma_fd", fs } };
        }
    }
    {  // meta-GGA: B95
        const json& g = ref["kernel"]["mgga"];
        const std::vector<double> rho = g["rho"].get<std::vector<double>>(), sigma = g["sigma"].get<std::vector<double>>(),
                                  tau = g["tau"].get<std::vector<double>>();
        const std::vector<std::string> kind = g["kind"].get<std::vector<std::string>>();
        const double floor = g["b95"]["sigma_floor"].get<double>();
        auto erho = [&](double r, double s, double t) {
            double e, a, b, c;
            qmxc::evaluateKernelMGGA(qmxc::XCKernel::B95_C, r, s, t, e, a, b, c);
            return r * e;
        };
        std::vector<double> eps, vr, vs, vt;
        json fr = json::array(), fs = json::array(), ft = json::array();
        for (size_t i = 0; i < rho.size(); ++i) {
            const double r = rho[i], s = std::max(sigma[i], floor), t = tau[i];
            double e, a, b, c;
            qmxc::evaluateKernelMGGA(qmxc::XCKernel::B95_C, r, s, t, e, a, b, c);
            eps.push_back(e);
            vr.push_back(a);
            vs.push_back(b);
            vt.push_back(c);
            // B95 is smooth in tau (no switch), but the stencil must not cross tau = tau_W (clamp): steps kept small
            const double tw = s / (8.0 * r);
            const bool ok = kind[i] == "admissible" && r >= 1e-10 && t >= 1e-18 && t - tw > 1e-6 * t;
            if (ok) {
                const double gap = (t - tw) / t;  // relative distance to tau_W
                const double hr = 1e-4 * r * std::min(1.0, gap), hs = 1e-4 * s * std::min(1.0, gap),
                             ht = 1e-4 * t * std::min(1.0, gap);
                fr.push_back(bestFD([&](double x) { return erho(x, s, t); }, r, hr, a));
                ft.push_back(bestFD([&](double x) { return erho(r, s, x); }, t, ht, c));
                if (sigma[i] > 0.0)
                    fs.push_back(bestFD([&](double x) { return erho(r, x, t); }, s, hs, b));
                else
                    fs.push_back(nullptr);
            } else {
                fr.push_back(nullptr);
                fs.push_back(nullptr);
                ft.push_back(nullptr);
            }
        }
        out["b95"] = { { "eps", eps }, { "vrho", vr }, { "vsigma", vs }, { "vtau", vt }, { "vrho_fd", fr },
            { "vsigma_fd", fs }, { "vtau_fd", ft } };
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

json runEnergy(const Args& a, const Mol& m)
{
    DoubleHybridMethod dh(a.method, methodConfig(a));
    if (!dh.setMolecule(m)) return { { "error", true }, { "message", dh.getErrorMessage() } };
    const double E = dh.calculateEnergy(a.gradient);
    json out = dh.getEnergyDecomposition();
    out["energy"] = E;
    out["error"] = dh.hasError();
    out["message"] = dh.getErrorMessage();
    out["has_gradient"] = dh.hasGradient();
    return out;
}

int parse(int argc, char** argv, Args& a)
{
    if (argc < 2) return 1;
    a.mode = argv[1];
    std::vector<std::string> pos;
    for (int i = 2; i < argc; ++i) {
        const std::string s = argv[i];
        auto next = [&]() -> std::string { return i + 1 < argc ? argv[++i] : ""; };
        if (s == "--method") a.method = next();
        else if (s == "--grid") a.grid = next();
        else if (s == "--threads") a.threads = std::stoi(next());
        else if (s == "--pt2") a.pt2 = next();
        else if (s == "--p86") a.p86 = next();
        else if (s == "--dispersion") a.dispersion = next();
        else if (s == "--ao-threshold") a.ao_threshold = std::stod(next());
        else if (s == "--fc") a.fc = true;
        else if (s == "--gradient") a.gradient = true;
        else pos.push_back(s);
    }
    if (!pos.empty()) a.xyz = pos[0];
    if (pos.size() > 1) a.xyz2 = pos[1];
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    Args a;
    if (parse(argc, argv, a)) {
        std::cerr << "usage: test_qm_dh kernel <ref.json> | energy <xyz> --method <name> [...] | update <xyzA> <xyzB>\n";
        return 2;
    }
    CurcumaLogger::set_verbosity(0);
    if (a.mode == "kernel") return modeKernel(a.xyz);
    try {
        if (a.mode == "energy") {
            Mol m;
            if (!readMol(a.xyz, m)) return 2;
            std::cout << runEnergy(a, m).dump(1) << std::endl;
            return 0;
        }
        if (a.mode == "update") {
            Mol A, B;
            if (!readMol(a.xyz, A) || !readMol(a.xyz2, B)) return 2;
            DoubleHybridMethod moving(a.method, methodConfig(a));
            moving.setMolecule(A);
            moving.calculateEnergy(false);
            moving.updateGeometry(B.m_geometry);
            const double e_moving = moving.calculateEnergy(false);
            const json fresh = runEnergy(a, B);
            std::cout << json({ { "moving", e_moving }, { "fresh", fresh["energy"] } }).dump(1) << std::endl;
            return 0;
        }
    } catch (const std::exception& ex) {
        std::cout << json({ { "exception", ex.what() } }).dump(1) << std::endl;
        return 0;
    }
    return 2;
}
