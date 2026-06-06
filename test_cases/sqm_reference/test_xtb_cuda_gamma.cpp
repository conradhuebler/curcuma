/*
 * Stage 3c component test: device Coulomb γ matrix vs the CPU reference.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 3c).
 *
 * Builds the native xTB γ matrix on the CPU (XTB::getGammaMatrix), computes it on
 * the device (k_gamma, Klopman–Ohno with host-precomputed per-shell hardness),
 * and compares elementwise @1e-9. Pairwise (no reductions) → ~ULP agreement.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_gamma [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/curcuma_logger.h"
#include "test_cuda_basis.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using curcuma::xtb::MethodType;
using curcuma::xtb::XTB;
using curcuma::xtb::GpuBasisFlat;
using curcuma::xtb::GpuH0Flat;
using curcuma::xtb::gpu::XtbGpuContext;

namespace {

int Z_of(const std::string& s)
{
    static const std::map<std::string, int> m = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},
        {"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
        {"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},{"Br",35},{"I",53}};
    auto it = m.find(s);
    return (it == m.end()) ? -1 : it->second;
}

bool readXYZ(const std::string& path, std::vector<int>& atoms, std::vector<double>& coords_ang)
{
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return false; }
    int nat = 0;
    f >> nat;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    atoms.clear();
    coords_ang.assign(3 * nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        std::string sym; double x, y, z;
        if (!(f >> sym >> x >> y >> z)) { std::cerr << "xyz parse error\n"; return false; }
        const int zn = Z_of(sym);
        if (zn < 0) { std::cerr << "unknown element: " << sym << "\n"; return false; }
        atoms.push_back(zn);
        coords_ang[3*i + 0] = x; coords_ang[3*i + 1] = y; coords_ang[3*i + 2] = z;
    }
    return true;
}

} // namespace

int main(int argc, char* argv[])
{
    CurcumaLogger::set_verbosity(0);

    double tol = 1.0e-9;
    bool quiet = false;
    std::string method;
    std::vector<std::string> xyz_paths;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--quiet") quiet = true;
        else if (a == "--tol" && i + 1 < argc) tol = std::stod(argv[++i]);
        else if (method.empty() && (a == "gfn1" || a == "gfn2")) method = a;
        else xyz_paths.push_back(a);
    }
    if (method.empty() || xyz_paths.empty()) {
        std::cerr << "usage: test_xtb_cuda_gamma [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_gamma: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        XTB xtb(mt);
        if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);

        const Eigen::MatrixXd& G_cpu = xtb.getGammaMatrix();
        const int nsh = static_cast<int>(G_cpu.rows());

        GpuBasisFlat bf; GpuH0Flat hf;
        xtb.exportGpuBasis(bf, hf);
        if (bf.nsh != nsh) { std::cerr << xyz << ": nsh mismatch\n"; ++failures; continue; }

        auto bd = makeGpuBasisData(bf, hf);
        if (!ctx.beginBasis(bd) || !ctx.computeIntegrals(bf.xyz_bohr.data())) {
            std::cerr << xyz << ": device integral build failed\n"; ++failures; continue;
        }
        std::vector<double> G_gpu(static_cast<size_t>(nsh) * nsh);
        if (!ctx.downloadGamma(G_gpu.data())) {
            std::cerr << xyz << ": device gamma download failed\n"; ++failures; continue;
        }

        double dG = 0.0;
        for (int i = 0; i < nsh; ++i)
            for (int j = 0; j < nsh; ++j) {
                const size_t cm = static_cast<size_t>(i) + static_cast<size_t>(j) * nsh;
                dG = std::max(dG, std::fabs(G_gpu[cm] - G_cpu(i, j)));
            }

        const bool ok = (dG <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s %s: max|dGamma|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(), dG, tol);
    }
    return failures == 0 ? 0 : 1;
}
