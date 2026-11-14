/*
 * <Unit test for the automatic z threshold of the path-CV reference histogram.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * Synthetic excess distributions z - z_path(s):
 *  - on-path frames: chi^2 with 3 degrees of freedom scaled to median 0.01 A^2, plus a 1 %
 *    outlier tail log-uniform in [0.01, 0.2] A^2;
 *  - far frames (other basin, permuted labels): log-normal around 0.7 A^2 (sigma_ln 0.2).
 * Bimodal case (70 % on-path, 30 % far): >= 99 % of the far frames discarded and >= 99 % of
 * the on-path frames kept. Unimodal case (on-path only): <= 1 % discarded. 8 seeds each.
 */

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "src/capabilities/z_filter.h"

static int failures = 0;
static void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
    if (!ok)
        ++failures;
}

static double onPath(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) < 0.01) // outlier tail, log-uniform in [0.01, 0.2]
        return 0.01 * std::pow(20.0, u(rng));
    std::chi_squared_distribution<double> chi(3.0);
    return chi(rng) * 0.01 / 2.365973884; // median of chi^2_3 = 2.36597
}

static double farFrame(std::mt19937_64& rng)
{
    std::normal_distribution<double> g(0.0, 0.2);
    return 0.7 * std::exp(g(rng));
}

int main()
{
    const int N = 20000;
    const int seeds = 8;
    std::cout << "Bimodal: 70 % on-path + 30 % far, N = " << N << ", " << seeds << " seeds" << std::endl;
    for (int sd = 0; sd < seeds; ++sd) {
        std::mt19937_64 rng(1000 + sd);
        std::vector<double> e;
        std::vector<int> far;
        for (int k = 0; k < N; ++k) {
            const bool f = (k % 10) < 3;
            e.push_back(f ? farFrame(rng) : onPath(rng));
            far.push_back(f);
        }
        const auto r = curcuma::z_filter::chooseThreshold(e);
        long nf = 0, nf_rej = 0, no = 0, no_keep = 0;
        for (int k = 0; k < N; ++k) {
            if (far[k]) {
                ++nf;
                nf_rej += e[k] > r.threshold;
            } else {
                ++no;
                no_keep += e[k] <= r.threshold;
            }
        }
        const double fr = double(nf_rej) / nf, ok = double(no_keep) / no;
        check(r.bimodal && fr >= 0.99 && ok >= 0.99,
            "seed " + std::to_string(sd) + ": threshold " + std::to_string(r.threshold) + " A^2 (" + r.rule
                + "), far discarded " + std::to_string(fr) + ", on-path kept " + std::to_string(ok));
    }
    std::cout << "Unimodal: on-path only, N = " << N << ", " << seeds << " seeds" << std::endl;
    for (int sd = 0; sd < seeds; ++sd) {
        std::mt19937_64 rng(2000 + sd);
        std::vector<double> e;
        for (int k = 0; k < N; ++k)
            e.push_back(onPath(rng));
        const auto r = curcuma::z_filter::chooseThreshold(e);
        long rej = 0;
        for (double v : e)
            rej += v > r.threshold;
        const double fr = double(rej) / N;
        check(!r.bimodal && fr <= 0.01,
            "seed " + std::to_string(sd) + ": threshold " + std::to_string(r.threshold) + " A^2 (" + r.rule
                + "), discarded " + std::to_string(fr));
    }
    // degenerate inputs: no filter
    check(std::isinf(curcuma::z_filter::chooseThreshold(std::vector<double>(10, 0.01)).threshold),
        "fewer than 20 frames -> no filter");
    check(std::isinf(curcuma::z_filter::chooseThreshold(std::vector<double>(100, -0.001)).threshold),
        "median excess <= 0 -> no filter");
    std::cout << (failures ? "FAILED" : "ALL PASSED") << " (" << failures << " failures)" << std::endl;
    return failures ? 1 : 0;
}
