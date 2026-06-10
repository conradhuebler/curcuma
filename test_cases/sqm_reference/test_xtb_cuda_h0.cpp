/*
 * Stage 3b component test: device overlap S, bare Hamiltonian H0 and Cholesky L
 * vs the CPU reference.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 3b).
 *
 * Builds the native xTB basis on the CPU (XTB::getOverlap / getBareHamiltonian
 * are the references), computes S/H0 on the device (k_overlap_h0) and L = chol(S)
 * (cusolverDnDpotrf), and compares elementwise. S/H0 are symmetric so the
 * column-major device buffer equals the row-major CPU matrix value at (i,j). The
 * Cholesky factor is validated by reconstructing S = L·Lᵀ from the device L and
 * comparing to the device S (no CPU L getter needed). Pairwise integrals (no
 * reductions) → bit-level agreement expected; gate at 1e-9 for pow/exp ULP.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_h0 [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
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
        std::cerr << "usage: test_xtb_cuda_h0 [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_h0: no usable CUDA device\n";
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
        xtb.Calculation(false);  // build basis + S + H0 (CN/H0 are pre-SCF)

        const Matrix& S_cpu  = xtb.getOverlap();
        const Matrix& H0_cpu = xtb.getBareHamiltonian();
        const int nao = static_cast<int>(S_cpu.rows());

        GpuBasisFlat bf; GpuH0Flat hf;
        xtb.exportGpuBasis(bf, hf);
        if (bf.nao != nao) { std::cerr << xyz << ": nao mismatch\n"; ++failures; continue; }

        auto bd = makeGpuBasisData(bf, hf);
        if (!ctx.beginBasis(bd) || !ctx.computeIntegrals(bf.xyz_bohr.data())) {
            std::cerr << xyz << ": device integral build failed\n"; ++failures; continue;
        }
        const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
        std::vector<double> S_gpu(nn), H0_gpu(nn), L_gpu(nn);
        if (!ctx.downloadOverlap(S_gpu.data()) || !ctx.downloadH0(H0_gpu.data())
            || !ctx.downloadCholesky(L_gpu.data())) {
            std::cerr << xyz << ": device download failed\n"; ++failures; continue;
        }

        // S / H0 elementwise (column-major device vs row-major CPU; both symmetric).
        double dS = 0.0, dH0 = 0.0;
        for (int i = 0; i < nao; ++i)
            for (int j = 0; j < nao; ++j) {
                const size_t cm = static_cast<size_t>(i) + static_cast<size_t>(j) * nao;
                dS  = std::max(dS,  std::fabs(S_gpu[cm]  - S_cpu(i, j)));
                dH0 = std::max(dH0, std::fabs(H0_gpu[cm] - H0_cpu(i, j)));
            }

        // Cholesky consistency: reconstruct S = L·Lᵀ from the device lower L
        // (column-major: L(i,k) at i + k*nao, lower triangle k<=i) and compare to
        // the device S.
        double dChol = 0.0;
        for (int i = 0; i < nao; ++i)
            for (int j = 0; j < nao; ++j) {
                double acc = 0.0;
                const int kmax = std::min(i, j);
                for (int k = 0; k <= kmax; ++k)
                    acc += L_gpu[static_cast<size_t>(i) + static_cast<size_t>(k) * nao]
                         * L_gpu[static_cast<size_t>(j) + static_cast<size_t>(k) * nao];
                const size_t cm = static_cast<size_t>(i) + static_cast<size_t>(j) * nao;
                dChol = std::max(dChol, std::fabs(acc - S_gpu[cm]));
            }

        const bool ok = (dS <= tol) && (dH0 <= tol) && (dChol <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok) {
            std::printf("%s %-28s %s: max|dS|=%.3e max|dH0|=%.3e max|LLt-S|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(),
                        dS, dH0, dChol, tol);
        }
    }
    return failures == 0 ? 0 : 1;
}
