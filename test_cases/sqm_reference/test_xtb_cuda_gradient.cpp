/*
 * Stage 4a system test: full GFN1 nuclear gradient on the device vs CPU.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 4a).
 *
 * Runs the same molecule twice: once on the bare CPU XTB (Calculation(true) →
 * host calculateGradient) and once with a device-resident SCF + device gradient
 * backend installed (XTB::calculateGradientGpu: repulsion/H0-Pulay/Coulomb on the
 * GPU, dispersion + CN chain-rule on the host). Compares getGradient() (Eh/Å)
 * elementwise. GFN1 exercises the device gradient; GFN2 currently falls back to
 * the CPU gradient (Stage 4b pending) so it is a pure SCF-path regression check.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_gradient [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/xtb_gpu_method.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/curcuma_logger.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using curcuma::xtb::MethodType;
using curcuma::xtb::XTB;
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

    double tol = 1.0e-7;
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
        std::cerr << "usage: test_xtb_cuda_gradient [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_gradient: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        // Tight SCF: the analytic gradient is only as accurate as the SCF fixed
        // point. The GPU device-resident loop (Stage 6) mixes charges on the device
        // (≈ the host Broyden to ~1e-11/step), so it lands at a slightly different
        // point in the 1e-5 default convergence ball than the CPU host Broyden —
        // tightening both to 1e-8 shrinks the ball below the gradient tolerance.
        const double tight = 1.0e-8;

        // CPU reference.
        XTB xtb_cpu(mt);
        if (!xtb_cpu.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": CPU InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb_cpu.setScfThreshold(tight);
        xtb_cpu.Calculation(true);
        const Matrix G_cpu = xtb_cpu.Gradient();

        // GPU path: device-resident SCF + device gradient backend.
        auto backend = createXtbGpuScfBackend(&ctx);
        XTB xtb_gpu(mt);
        xtb_gpu.setGpuScfBackend(backend.get());
        if (!xtb_gpu.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": GPU InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb_gpu.setScfThreshold(tight);
        xtb_gpu.Calculation(true);
        const Matrix G_gpu = xtb_gpu.Gradient();

        if (G_cpu.rows() != nat || G_gpu.rows() != nat) {
            std::cerr << xyz << ": gradient shape mismatch\n"; ++failures; continue;
        }
        double dmax = 0.0;
        for (int i = 0; i < nat; ++i)
            for (int k = 0; k < 3; ++k)
                dmax = std::max(dmax, std::fabs(G_cpu(i, k) - G_gpu(i, k)));

        const bool ok = (dmax <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s %s: max|dGrad|=%.3e Eh/A  tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(), dmax, tol);
    }
    return failures == 0 ? 0 : 1;
}
