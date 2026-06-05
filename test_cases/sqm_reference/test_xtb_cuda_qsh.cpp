/*
 * Stage 6 (S6.2) component test: device shell Mulliken charges from the resident
 * density vs the host updatePopulations.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 6).
 *
 * Shell analogue of test_xtb_cuda_qat. Drives a full device-resident GFN2 SCF
 * (createXtbGpuScfBackend installed on a bare XTB), then reduces the resident AO
 * populations into per-shell charges on the device (XtbGpuContext::
 * residentShellCharges via GpuScfBackend::shellCharges, AO->shell map) and
 * compares them to the host m_wfn.q_sh (= getShellCharges()), which the resident
 * path computed from the SAME downloaded pop_ao. Both derive from one pop_ao, so
 * they agree to ~1e-13 (gate 1e-9 absorbs atomicAdd order).
 *
 * GFN2 only (the resident ao2sh map is set up on the GFN2 device-potential path).
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_qsh [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]
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
        std::cerr << "usage: test_xtb_cuda_qsh [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]\n";
        return 2;
    }

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_qsh: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        // Device-resident GFN2 SCF; the converged AO populations stay resident.
        auto backend = createXtbGpuScfBackend(&ctx);
        XTB xtb(MethodType::GFN2);
        xtb.setGpuScfBackend(backend.get());
        if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);

        // Host reference: q_sh from the resident path (downloaded converged charges).
        const Vector q_host = xtb.getShellCharges();
        const Vector n0_sh  = xtb.getReferenceShellOccupations();
        const int nsh = static_cast<int>(n0_sh.size());

        // Device reduction of the SAME resident populations.
        Vector q_gpu;
        if (!backend->shellCharges(n0_sh, q_gpu) || q_gpu.size() != nsh) {
            std::cerr << xyz << ": device shellCharges failed (resident path inactive?)\n";
            ++failures; continue;
        }
        if (q_host.size() != nsh) {
            std::cerr << xyz << ": host q_sh size mismatch\n"; ++failures; continue;
        }

        double dq = 0.0;
        for (int i = 0; i < nsh; ++i) dq = std::max(dq, std::fabs(q_host(i) - q_gpu(i)));
        const bool ok = (dq <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s gfn2: max|dq_sh|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), dq, tol);
    }
    return failures == 0 ? 0 : 1;
}
