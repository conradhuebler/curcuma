/*
 * Stage 3a component test: device CN + self-energy vs the CPU reference.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 3a).
 *
 * Builds the native xTB basis on the CPU, exports the flattened basis, computes
 * the coordination numbers (cn_exp for GFN1 / cn_gfn for GFN2) and CN-shifted
 * shell self-energies on the device (XtbGpuContext), and compares them
 * elementwise against the CPU reference (the same cn_exp/cn_gfn free functions
 * and se = selfenergy − kcn·CN[sh2at]). The integrals are pairwise (no
 * reductions) so the device should match to ~1e-12; the gate is 1e-9 to absorb
 * device exp() ULP differences.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0, so it never
 * false-fails off-GPU.
 *
 * Usage:  test_xtb_cuda_cn [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/parameters/xtb_params_extra.hpp"
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
    std::getline(f, line);  // rest of first line
    std::getline(f, line);  // comment line
    atoms.clear();
    coords_ang.assign(3 * nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        std::string sym; double x, y, z;
        if (!(f >> sym >> x >> y >> z)) {
            std::cerr << "xyz parse error at atom " << i << " in " << path << "\n";
            return false;
        }
        const int zn = Z_of(sym);
        if (zn < 0) { std::cerr << "unknown element: " << sym << "\n"; return false; }
        atoms.push_back(zn);
        coords_ang[3*i + 0] = x;
        coords_ang[3*i + 1] = y;
        coords_ang[3*i + 2] = z;
    }
    return true;
}

double maxAbsDiff(const std::vector<double>& a, const std::vector<double>& b)
{
    double m = 0.0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
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
        std::cerr << "usage: test_xtb_cuda_cn [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;
    const bool is_gfn2 = (mt == MethodType::GFN2);

    // Device handshake: skip cleanly (pass) when no usable CUDA device.
    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_cn: no usable CUDA device\n";
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
            std::cerr << xyz << ": XTB::InitialiseMolecule failed\n"; ++failures; continue;
        }
        // Full single point to guarantee the basis + H0 parameters are built
        // (CN itself is built pre-SCF, so non-convergence does not matter here).
        xtb.Calculation(false);

        GpuBasisFlat bf;
        GpuH0Flat hf;
        xtb.exportGpuBasis(bf, hf);
        const int nsh = bf.nsh;

        // CPU reference: the exact cn_exp/cn_gfn free functions + the self-energy
        // formula, evaluated on the SAME flattened z / xyz the device consumes.
        std::vector<double> cn_cpu = is_gfn2
            ? curcuma::xtb::cn_gfn(bf.z, bf.xyz_bohr)
            : curcuma::xtb::cn_exp(bf.z, bf.xyz_bohr);
        std::vector<double> se_cpu(nsh);
        for (int s = 0; s < nsh; ++s)
            se_cpu[s] = hf.selfenergy[s] - hf.kcn[s] * cn_cpu[bf.sh2at[s]];

        // Device path.
        auto bd = makeGpuBasisData(bf, hf);
        if (!ctx.beginBasis(bd) || !ctx.computeCnSelfEnergy(bf.xyz_bohr.data())) {
            std::cerr << xyz << ": device CN compute failed\n"; ++failures; continue;
        }
        std::vector<double> cn_gpu(nat, 0.0), se_gpu(nsh, 0.0);
        if (!ctx.downloadCn(cn_gpu.data()) || !ctx.downloadSelfEnergy(se_gpu.data())) {
            std::cerr << xyz << ": device download failed\n"; ++failures; continue;
        }

        const double dcn = maxAbsDiff(cn_cpu, cn_gpu);
        const double dse = maxAbsDiff(se_cpu, se_gpu);
        const bool ok = (dcn <= tol) && (dse <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok) {
            std::printf("%s %-28s %s: max|dCN|=%.3e max|dSE|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(), dcn, dse, tol);
        }
    }
    return failures == 0 ? 0 : 1;
}
