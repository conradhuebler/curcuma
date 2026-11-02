/*
 * Spin-polarised XC kernels of the open-shell Kohn-Sham path (stage 5b, qm_xc_spin.cpp) at fixed points.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 *   test_qm_uks_kernel <reference.json>
 *
 * Reads the point set (rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb) of reference.json["kernel"], evaluates every kernel
 * listed there through qmxc::evaluatePointSpin (one-term functional, weight 1, density threshold 0) and prints one JSON
 * object: per kernel eps = e / rho, vrho[2], vsigma[3] per point, plus two self-checks per kernel:
 *   fd_max   largest |v - fd| / (|v| + |fd| + 1e-13 |e| / h) against a central finite difference of e (h = 1e-4 x)
 *   cs_max   closed-shell limit: for the total rho and sigma of every point, the spin form at rho_a = rho_b = rho/2,
 *            sigma_aa = sigma_ab = sigma_bb = sigma/4 against the closed-shell kernel of qm_xc_kernels.h: largest
 *            |d| / (1 + |ref|) of eps, v_rho (v_rho_a) and 2 v_sigma (v_sigma_aa + v_sigma_ab / 2)
 * check_qm_uhf.py --check kernel compares eps/vrho/vsigma with libxc.
 *
 * Claude Generated (Nov 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/qm_xc_spin.h"

#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

using json = nlohmann::json;

namespace {
double relerr(double a, double b) { return std::abs(a - b) / std::max(1e-12, std::max(std::abs(a), std::abs(b))); }
/// |a - b| / (1 + |b|): the scaled metric of the libxc comparison (a.u. values; no relative test of ~0 by cancellation)
double scalederr(double a, double b) { return std::abs(a - b) / (1.0 + std::abs(b)); }
}  // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: test_qm_uks_kernel <reference.json>\n";
        return 1;
    }
    std::ifstream f(argv[1]);
    const json ref = json::parse(f);
    const json& pts = ref["kernel"]["points"];
    const std::map<std::string, qmxc::XCKernel> names = {
        { "LDA_X", qmxc::XCKernel::SlaterX },     { "LDA_C_VWN", qmxc::XCKernel::VWN5_C },
        { "LDA_C_VWN_RPA", qmxc::XCKernel::VWNRPA_C }, { "LDA_C_PW_MOD", qmxc::XCKernel::PW92_C },
        { "GGA_C_PBE", qmxc::XCKernel::PBE_C },   { "GGA_C_LYP", qmxc::XCKernel::LYP_C },
        { "GGA_X_PBE", qmxc::XCKernel::PBE_X },   { "GGA_X_B88", qmxc::XCKernel::B88_X } };
    json out;
    for (const auto& [name, k] : names) {
        qmxc::XCFunctional fn;
        fn.name = name;
        fn.exact_exchange = 0.0;
        fn.terms = { { k, 1.0 } };
        const bool gga = qmxc::isGGAKernel(k);
        const double rho_min = qmxc::kern::effectiveThreshold(gga, 0.0);
        json eps = json::array(), vrho = json::array(), vsig = json::array();
        double fd_max = 0.0, cs_max = 0.0;
        json fd_worst = json::array(), cs_worst = json::array();
        for (const auto& p : pts) {
            double x[5];
            for (int i = 0; i < 5; ++i) x[i] = p[i].get<double>();
            double ex = 0, ec = 0, v[5];
            qmxc::evaluatePointSpin(fn, gga, rho_min, x[0], x[1], x[2], x[3], x[4], ex, ec, v);
            const double rho = x[0] + x[1];
            eps.push_back((ex + ec) / rho);
            vrho.push_back({ v[0], v[1] });
            if (gga) vsig.push_back({ v[2], v[3], v[4] });
            // finite difference of e in each variable (only where the variable is clearly inside its domain)
            for (int i = 0; i < (gga ? 5 : 2); ++i) {
                if (x[i] <= 1e-8 || (i < 2 && x[i] < 1e-6 * rho)) continue;
                const double h = 1e-4 * std::abs(x[i]);
                double y[5];
                for (int j = 0; j < 5; ++j) y[j] = x[j];
                double e[2];
                for (int s = 0; s < 2; ++s) {
                    y[i] = x[i] + (s == 0 ? h : -h);
                    double a = 0, b = 0, w[5];
                    qmxc::evaluatePointSpin(fn, gga, rho_min, y[0], y[1], y[2], y[3], y[4], a, b, w);
                    e[s] = a + b;
                }
                const double fd = (e[0] - e[1]) / (2 * h);
                // relative error with the rounding noise of the difference quotient, 1e-13 |e| / h, in the denominator
                const double noise = 1e-13 * std::max(std::abs(ex + ec), 1e-300) / h;
                // a variable whose whole contribution |v x| is below 1e-6 |e| is not resolvable by a difference quotient
                if ((std::abs(v[i]) > 1e-10 || std::abs(fd) > 1e-10) && std::abs(v[i] * x[i]) > 1e-6 * std::abs(ex + ec)) {
                    const double r = std::abs(v[i] - fd) / (std::abs(v[i]) + std::abs(fd) + noise);
                    if (r > fd_max) fd_worst = { static_cast<int>(eps.size()) - 1, i, v[i], fd };
                    fd_max = std::max(fd_max, r);
                }
            }
            // closed-shell limit at the same total density and sigma: rho_a = rho_b = rho/2,
            // sigma_aa = sigma_ab = sigma_bb = sigma/4; then v_rho_a = v_rho and v_sigma_aa + v_sigma_ab/2 = 2 v_sigma
            {
                const double sigma = gga ? std::max(0.0, x[2] + 2 * x[3] + x[4]) : 0.0;
                double e_cs = 0, vr_cs = 0, vs_cs = 0;
                if (qmxc::kern::evaluateKernel(k, rho, sigma, e_cs, vr_cs, vs_cs) && rho >= rho_min) {
                    double a = 0, b = 0, w[5];
                    qmxc::evaluatePointSpin(fn, gga, rho_min, 0.5 * rho, 0.5 * rho, 0.25 * sigma, 0.25 * sigma, 0.25 * sigma, a, b, w);
                    const double c3[3] = { scalederr((a + b) / rho, e_cs), scalederr(w[0], vr_cs),
                                           gga ? scalederr(w[2] + 0.5 * w[3], 2.0 * vs_cs) : 0.0 };
                    for (int q = 0; q < 3; ++q)
                        if (c3[q] > cs_max) {
                            cs_max = c3[q];
                            cs_worst = { static_cast<int>(eps.size()) - 1, q };
                        }
                }
            }
        }
        out[name] = { { "eps", eps }, { "vrho", vrho }, { "fd_max", fd_max }, { "fd_worst", fd_worst }, { "cs_max", cs_max }, { "cs_worst", cs_worst } };
        if (gga) out[name]["vsigma"] = vsig;
    }
    std::cout << out.dump() << std::endl;
    return 0;
}
