/*
 * Stage 5 (Part B2) component test: in-SCF GFN2 D4 atom-potential dE_D4/dq on the
 * device vs the host D4Evaluator, at a frozen set of charges.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 5).
 *
 * Runs a CPU GFN2 SCF to get a converged Mulliken charge set q_ref, then
 * evaluates XTB::computeD4PotentialDedq(q_ref) on (a) the bare CPU XTB (host
 * D4Evaluator per-reference path) and (b) an XTB with the device-resident
 * backend installed (device k_d4_dedq: host-built W/dWq + device 7×7 contraction
 * × BJ disp_sum). Both use the SAME q_ref, so the comparison isolates the kernel
 * from any SCF difference. The contraction is a 49-term sum reassociated per atom,
 * so the device matches the host to ~1e-12 (gate 1e-9).
 *
 * GFN2 only (D4 charge coupling is GFN2). No CUDA device -> SKIP, exit 0.
 *
 * Usage:  test_xtb_cuda_d4pot [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/xtb_gpu_method.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"

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
    if (method != "gfn2" || xyz_paths.empty()) {
        std::cerr << "usage: test_xtb_cuda_d4pot [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]\n";
        return 2;
    }

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_d4pot: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        // CPU reference: converged Mulliken charges + host D4 potential dE/dq.
        XTB xtb_cpu(MethodType::GFN2);
        if (!xtb_cpu.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": CPU InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb_cpu.Calculation(false);
        const Vector q_ref = xtb_cpu.getCharges();
        Vector dEdq_cpu;
        if (!xtb_cpu.computeD4PotentialDedq(q_ref, dEdq_cpu) || dEdq_cpu.size() != nat) {
            std::cerr << xyz << ": CPU computeD4PotentialDedq failed\n"; ++failures; continue;
        }

        // Device path: same q_ref, dE/dq via the resident D4 backend.
        auto backend = createXtbGpuScfBackend(&ctx);
        XTB xtb_gpu(MethodType::GFN2);
        xtb_gpu.setGpuScfBackend(backend.get());
        if (!xtb_gpu.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": GPU InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb_gpu.Calculation(false);  // primes the device D4 reference data (beginDispersion)
        Vector dEdq_gpu;
        if (!xtb_gpu.computeD4PotentialDedq(q_ref, dEdq_gpu) || dEdq_gpu.size() != nat) {
            std::cerr << xyz << ": GPU computeD4PotentialDedq failed\n"; ++failures; continue;
        }

        double dmax = 0.0;
        for (int i = 0; i < nat; ++i) dmax = std::max(dmax, std::fabs(dEdq_cpu(i) - dEdq_gpu(i)));
        const bool ok = (dmax <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s gfn2: max|d(dE/dq)|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), dmax, tol);
    }
    return failures == 0 ? 0 : 1;
}
