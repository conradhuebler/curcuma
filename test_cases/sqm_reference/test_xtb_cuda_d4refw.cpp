/*
 * Stage 6 (S6.2b) component test: device D4 per-atom reference weights W = gwk(CN)
 * ·zeta(q) and dWq = dW/dq, rebuilt on the device from the SCF charges, vs the
 * host D4ParameterGenerator::buildRefWFlat.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 6).
 *
 * Runs a CPU GFN2 SCF to prime the D4 generator and obtain a converged Mulliken
 * charge set q_ref, then builds the reference weights two ways at the SAME q_ref:
 * (a) the host XTB::buildD4RefWFlat (= the per-iteration host path that fed the
 * device in Stage 5), and (b) the device kernel k_d4_build_refw, primed with the
 * q-independent reference tables via exportD4RefWDeviceData -> ctx.
 * beginDispersionWeights and launched by ctx.dispersionBuildRefW(q_ref). The
 * device exists so the fused resident loop (S6.5) rebuilds W/dWq from the resident
 * charges with no host buildRefWFlat per iteration. Only entries j < nref[i] per
 * atom are physical; the device matches the host to ~1e-12 (gate 1e-9).
 *
 * GFN2 only (D4 charge coupling is GFN2). No CUDA device -> SKIP, exit 0.
 *
 * Usage:  test_xtb_cuda_d4refw [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]
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

// Mirrors D4ParameterGenerator::MAX_REF / XtbGpuContext Impl::D4_MAX_REF.
constexpr int D4_MAX_REF = 7;

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
        std::cerr << "usage: test_xtb_cuda_d4refw [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]\n";
        return 2;
    }

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_d4refw: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        // CPU GFN2 SCF primes m_d4_generator and gives the frozen charge set.
        XTB xtb(MethodType::GFN2);
        if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);
        const Vector q_ref = xtb.getCharges();
        if (q_ref.size() != nat) { std::cerr << xyz << ": charge size mismatch\n"; ++failures; continue; }

        // Host reference: per-atom W/dWq at q_ref (flat, nat*D4_MAX_REF).
        std::vector<double> W_host, dWq_host;
        xtb.buildD4RefWFlat(q_ref, W_host, dWq_host);
        if (static_cast<int>(W_host.size()) != nat * D4_MAX_REF
            || static_cast<int>(dWq_host.size()) != nat * D4_MAX_REF) {
            std::cerr << xyz << ": host buildD4RefWFlat empty/short (D4 generator not prepared?)\n";
            ++failures; continue;
        }

        // q-independent reference tables -> device.
        std::vector<double> cn, gi, zeff, refcn, refcovcn, refq;
        std::vector<int> nref;
        xtb.exportD4RefWDeviceData(cn, gi, zeff, refcn, refcovcn, refq, nref);
        if (static_cast<int>(cn.size()) != nat || static_cast<int>(nref.size()) != nat) {
            std::cerr << xyz << ": exportD4RefWDeviceData failed\n"; ++failures; continue;
        }
        if (!ctx.beginDispersionWeights(nat, cn.data(), gi.data(), zeff.data(), refcn.data(),
                                        refcovcn.data(), refq.data(), nref.data())) {
            std::cerr << xyz << ": device beginDispersionWeights failed\n"; ++failures; continue;
        }

        // Device rebuild at the SAME q_ref.
        std::vector<double> W_gpu(nat * D4_MAX_REF, 0.0), dWq_gpu(nat * D4_MAX_REF, 0.0);
        if (!ctx.dispersionBuildRefW(nat, q_ref.data(), W_gpu.data(), dWq_gpu.data())) {
            std::cerr << xyz << ": device dispersionBuildRefW failed\n"; ++failures; continue;
        }

        // Compare only the physical references (j < nref[i]).
        double dW = 0.0, dWq = 0.0;
        for (int i = 0; i < nat; ++i) {
            const int nr = (nref[i] < D4_MAX_REF) ? nref[i] : D4_MAX_REF;
            for (int j = 0; j < nr; ++j) {
                const int k = i * D4_MAX_REF + j;
                dW  = std::max(dW,  std::fabs(W_host[k]   - W_gpu[k]));
                dWq = std::max(dWq, std::fabs(dWq_host[k] - dWq_gpu[k]));
            }
        }
        const bool ok = (dW <= tol) && (dWq <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s gfn2: max|dW|=%.3e max|d(dWq)|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), dW, dWq, tol);
    }
    return failures == 0 ? 0 : 1;
}
