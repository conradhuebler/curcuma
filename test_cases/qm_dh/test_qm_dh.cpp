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
 *   test_qm_dh d4 <xyz> --method <functional> [--charge q] [--fd]
 *       WP12 stage 4: the native DFT-D4 term alone (curcuma::dispersion::D4DftDispersion, parameters of
 *       d4DftParameters() for <functional>): energy, two-body / ATM parts, EEQ charges, CN, analytic gradient
 *       (Eh/Bohr) and with --fd the 4-point central difference of the D4 energy (h = 5e-4 Bohr, every coordinate).
 *
 *   test_qm_dh dhd4 <xyz> --method <dh-name> [--p86 ..] [--fc] [--fd]
 *       WP12 stage 4: energy + gradient of the double hybrid (its dispersion per --dispersion), of the same method
 *       with --dispersion none, and of the D4 term alone (the D4 parameters of the method name without "-d4"):
 *       the differences E - E(none) and g - g(none) must be the D4 term; with --fd also the 4-point central
 *       difference of the whole energy (h = 1e-3 Bohr in Angstrom, every coordinate; the grid moves with the atoms).
 *
 *   test_qm_dh ksd4 <xyz> --method <pbe|b3lyp|tpss|scan> [--grid g]
 *       WP12 stage 4: QMMethod <name>-d4 against QMMethod <name> and the D4 term alone (energy, gradient).
 *
 * Claude Generated (Oct 2026, WP12 stage 1). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/double_hybrid_method.h"
#include "src/core/energy_calculators/qm_methods/qm_method.h"
#include "src/core/energy_calculators/qm_methods/qm_xc.h"
#include "src/core/energy_calculators/dispersion/d4_dft.h"
#include "src/core/units.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// Driver defaults = the PySCF/libxc reference convention (PZ81-P86, no frozen core), always passed explicitly to
// DoubleHybridMethod; curcuma's own DH defaults are vwn5 + frozen core since Oct 31 2026.
struct Args {
    std::string mode, xyz, xyz2, method = "b2plyp", grid = "sg1", pt2 = "ri", p86 = "pz81", dispersion = "auto";
    int threads = 4;
    int charge = 0;
    bool fc = false, gradient = false, fd = false;
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
        else if (s == "--fd") a.fd = true;
        else if (s == "--charge") a.charge = std::stoi(next());
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
        if (a.mode == "d4") {
            // WP12 stage 4 (Claude Generated, Oct 2026): the D4 term alone, with analytic and FD gradient
            Mol m;
            if (!readMol(a.xyz, m)) return 2;
            const auto* p = curcuma::dispersion::findD4DftParameters(a.method);
            if (!p) {
                std::cout << json({ { "error", "no D4 parameters for " + a.method } }).dump(1) << std::endl;
                return 0;
            }
            curcuma::dispersion::D4DftDispersion d4(*p);
            const Matrix xb = m.m_geometry * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
            double e = 0.0;
            Matrix g;
            json out;
            if (!d4.compute(m.m_atoms, xb, a.charge, true, e, g)) {
                std::cout << json({ { "error", d4.error() } }).dump(1) << std::endl;
                return 0;
            }
            out["energy"] = e;
            out["two_body"] = d4.twoBodyEnergy();
            out["three_body"] = d4.threeBodyEnergy();
            out["charges"] = std::vector<double>(d4.charges().data(), d4.charges().data() + d4.charges().size());
            out["cn"] = d4.coordinationNumbers();
            std::vector<double> gv, fdv;
            for (int i = 0; i < g.rows(); ++i)
                for (int k = 0; k < 3; ++k) gv.push_back(g(i, k));
            out["gradient"] = gv;
            if (a.fd) {
                curcuma::dispersion::D4DftDispersion d4fd(*p);
                const double h = 5e-4;
                for (int i = 0; i < xb.rows(); ++i)
                    for (int k = 0; k < 3; ++k) {
                        auto en = [&](double x) {
                            Matrix y = xb;
                            y(i, k) = x;
                            double ee = 0.0;
                            Matrix gg;
                            d4fd.compute(m.m_atoms, y, a.charge, false, ee, gg);
                            return ee;
                        };
                        fdv.push_back(fd4(en, xb(i, k), h));
                    }
                out["fd_gradient"] = fdv;
            }
            std::cout << out.dump(1) << std::endl;
            return 0;
        }
        if (a.mode == "dhd4" || a.mode == "ksd4") {
            // WP12 stage 4 (Claude Generated, Oct 2026): the method with D4 = the method without + the D4 term
            Mol m;
            if (!readMol(a.xyz, m)) return 2;
            auto flat = [](const Matrix& g) {
                std::vector<double> v;
                for (int i = 0; i < g.rows(); ++i)
                    for (int k = 0; k < 3; ++k) v.push_back(g(i, k));
                return v;
            };
            std::string key = a.method;
            if (key.size() > 3 && key.compare(key.size() - 3, 3, "-d4") == 0) key.resize(key.size() - 3);
            json out;
            auto run = [&](const std::string& disp, json& o, std::unique_ptr<ComputationalMethod>& keep) {
                Args b = a;
                b.dispersion = disp;
                json cfg = methodConfig(b);
                std::unique_ptr<ComputationalMethod> meth;
                if (a.mode == "dhd4") {
                    meth = std::make_unique<DoubleHybridMethod>(a.method, cfg);
                } else {
                    const QMFunctional f = key == "pbe" ? QMFunctional::PBE : key == "tpss" ? QMFunctional::TPSS
                        : key == "scan" ? QMFunctional::SCAN : QMFunctional::B3LYP;
                    meth = std::make_unique<QMMethod>(f, cfg, false);
                }
                meth->setMolecule(m);
                o["energy"] = meth->calculateEnergy(true);
                o["error"] = meth->hasError();
                o["message"] = meth->getErrorMessage();
                o["gradient"] = flat(meth->getGradient() * CurcumaUnit::Length::BOHR_TO_ANGSTROM);  // Eh/Bohr
                o["name"] = meth->getMethodName();
                keep = std::move(meth);
            };
            std::unique_ptr<ComputationalMethod> with, without;
            json ow, on;
            run(a.mode == "dhd4" ? a.dispersion : "d4", ow, with);
            run("none", on, without);
            out["with"] = ow;
            out["none"] = on;
            const auto* p = curcuma::dispersion::findD4DftParameters(key);
            if (p) {
                curcuma::dispersion::D4DftDispersion d4(*p);
                double e = 0.0;
                Matrix g;
                d4.compute(m.m_atoms, m.m_geometry * CurcumaUnit::Length::ANGSTROM_TO_BOHR, 0, true, e, g);
                out["d4"] = { { "energy", e }, { "gradient", flat(g) } };
            }
            if (a.fd) {
                const double h = 1e-3 * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
                std::vector<double> fdv;
                for (int i = 0; i < m.m_geometry.rows(); ++i)
                    for (int k = 0; k < 3; ++k) {
                        auto en = [&](double x) {
                            Matrix y = m.m_geometry;
                            y(i, k) = x;
                            with->updateGeometry(y);
                            return with->calculateEnergy(false);
                        };
                        fdv.push_back(fd4(en, m.m_geometry(i, k), h) * CurcumaUnit::Length::BOHR_TO_ANGSTROM);
                    }
                out["fd_gradient"] = fdv;
            }
            std::cout << out.dump(1) << std::endl;
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
