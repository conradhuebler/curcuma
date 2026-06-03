/*
 * Stage 3d component test: device GFN2 multipole integrals (dp_int/qp_int) vs the
 * CPU reference.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 3d).
 *
 * Builds the native GFN2 multipole integrals on the CPU (XTB::getDipoleIntegrals /
 * getQuadrupoleIntegrals via setupMultipole), computes them on the device
 * (k_multipole_ints: global integral → per-column origin shift → traceless
 * transform, using the resident overlap S), and compares elementwise @1e-9.
 * dp_int/qp_int are NOT symmetric; the device buffers are column-major contiguous
 * (3·nn / 6·nn), element (μ,ν) of block k at k·nn + μ + ν·nao.
 *
 * GFN2 only. Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_multipole [--quiet] [--tol V] <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/curcuma_logger.h"
#include "test_cuda_basis.hpp"

#include <array>
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
    std::vector<std::string> xyz_paths;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--quiet") quiet = true;
        else if (a == "--tol" && i + 1 < argc) tol = std::stod(argv[++i]);
        else xyz_paths.push_back(a);
    }
    if (xyz_paths.empty()) {
        std::cerr << "usage: test_xtb_cuda_multipole [--quiet] [--tol V] <xyz> [<xyz> ...]\n";
        return 2;
    }

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_multipole: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        XTB xtb(MethodType::GFN2);
        if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);

        const auto& dp_cpu = xtb.getDipoleIntegrals();      // array<MatrixXd,3>
        const auto& qp_cpu = xtb.getQuadrupoleIntegrals();  // array<MatrixXd,6>
        const int nao = static_cast<int>(dp_cpu[0].rows());

        GpuBasisFlat bf; GpuH0Flat hf;
        xtb.exportGpuBasis(bf, hf);
        if (bf.nao != nao) { std::cerr << xyz << ": nao mismatch\n"; ++failures; continue; }

        auto bd = makeGpuBasisData(bf, hf);
        if (!ctx.beginBasis(bd) || !ctx.computeIntegrals(bf.xyz_bohr.data())) {
            std::cerr << xyz << ": device integral build failed\n"; ++failures; continue;
        }
        const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
        std::vector<double> dp_gpu(3 * nn), qp_gpu(6 * nn);
        if (!ctx.downloadMultipoleInts(dp_gpu.data(), qp_gpu.data())) {
            std::cerr << xyz << ": device multipole download failed\n"; ++failures; continue;
        }

        double dDp = 0.0, dQp = 0.0;
        for (int mu = 0; mu < nao; ++mu)
            for (int nu = 0; nu < nao; ++nu) {
                const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
                for (int k = 0; k < 3; ++k)
                    dDp = std::max(dDp, std::fabs(dp_gpu[k * nn + mn] - dp_cpu[k](mu, nu)));
                for (int k = 0; k < 6; ++k)
                    dQp = std::max(dQp, std::fabs(qp_gpu[k * nn + mn] - qp_cpu[k](mu, nu)));
            }

        const bool ok = (dDp <= tol) && (dQp <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s gfn2: max|dDp|=%.3e max|dQp|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), dDp, dQp, tol);
    }
    return failures == 0 ? 0 : 1;
}
