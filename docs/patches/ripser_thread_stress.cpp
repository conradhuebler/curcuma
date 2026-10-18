// Claude Generated (Oct 2026): ripser thread-safety stress test on real molecular patches.
// Build against a ripser.h (header-only): g++ -O2 -std=c++17 -pthread -I<dir with ripser.h> ripser_thread_stress.cpp -o stress
// Run: ./stress conf.xyz <nstruct> <threads> > barcodes.txt   (multi-XYZ input; one patch per atom, radius 5 A, H0+H1)
// Compare the barcodes of a 1-thread run with a many-thread run (cmp); the unpatched header aborts with threads > 1.
// Written for the Ripser fork commit 5512434 (conradhuebler/ripser, function-local statics removed).
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include "ripser.h"
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <map>
#include <cstdio>
using Pt = std::array<double, 3>;
static std::string run_patch(const std::vector<Pt>& P, const std::vector<int>& idx)
{
    const int m = idx.size();
    std::vector<float> d;
    for (int a = 1; a < m; ++a)
        for (int b = 0; b < a; ++b) {
            double s = 0;
            for (int k = 0; k < 3; ++k) { double x = P[idx[a]][k] - P[idx[b]][k]; s += x * x; }
            d.push_back((float)std::sqrt(s));
        }
    compressed_lower_distance_matrix dist(std::move(d));
    value_t enc = std::numeric_limits<value_t>::infinity();
    for (size_t i = 0; i < dist.size(); ++i) {
        value_t r = -std::numeric_limits<value_t>::infinity();
        for (size_t j = 0; j < dist.size(); ++j) r = std::max(r, dist(i, j));
        enc = std::min(enc, r);
    }
    auto res = ripser<compressed_lower_distance_matrix>(std::move(dist), 1, enc, 1.0f, 2).compute_barcodes();
    std::ostringstream o;
    std::map<int, std::vector<std::pair<double, double>>> sorted;
    for (auto& a : res) { auto v = a.second; std::sort(v.begin(), v.end()); sorted[a.first] = v; }
    for (auto& a : sorted) { o << a.first << ":"; for (auto& p : a.second) o << " " << p.first << "," << p.second; o << ";"; }
    return o.str();
}
int main(int argc, char** argv)
{
    // usage: stress conf.xyz nstruct threads
    std::ifstream in(argv[1]);
    int nstruct = std::atoi(argv[2]), threads = std::atoi(argv[3]);
    std::vector<std::vector<Pt>> S;
    std::string line;
    while (std::getline(in, line) && (int)S.size() < nstruct) {
        int n = std::atoi(line.c_str());
        std::getline(in, line);
        std::vector<Pt> P;
        for (int i = 0; i < n; ++i) { std::getline(in, line); std::istringstream l(line); std::string e; Pt p; l >> e >> p[0] >> p[1] >> p[2]; P.push_back(p); }
        S.push_back(P);
    }
    struct Job { int s, c; };
    std::vector<Job> jobs;
    for (int s = 0; s < (int)S.size(); ++s) for (int c = 0; c < (int)S[s].size(); ++c) jobs.push_back({ s, c });
    std::vector<std::string> out(jobs.size());
    auto t0 = std::chrono::steady_clock::now();
    std::atomic<int> next(0);
    auto work = [&]() {
        for (int i = next++; i < (int)jobs.size(); i = next++) {
            const auto& P = S[jobs[i].s];
            std::vector<int> idx;
            for (int j = 0; j < (int)P.size(); ++j) {
                double s = 0; for (int k = 0; k < 3; ++k) { double x = P[j][k] - P[jobs[i].c][k]; s += x * x; }
                if (std::sqrt(s) <= 5.0) idx.push_back(j);
            }
            out[i] = run_patch(P, idx);
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) pool.emplace_back(work);
    for (auto& t : pool) t.join();
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::hash<std::string> h; size_t acc = 0;
    for (auto& s : out) acc = acc * 1000003u + h(s);
    std::fprintf(stderr, "patches %zu threads %d time %.2f s\n", jobs.size(), threads, sec);
    for (auto& s : out) std::cout << s << "\n";
    return 0;
}
