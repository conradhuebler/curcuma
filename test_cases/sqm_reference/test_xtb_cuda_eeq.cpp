/*
 * Stage 5 (Part A) component test: device single-shot D4 EEQ charge model vs the
 * CPU reference (curcuma::dispersion::D4ChargeModel).
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 5).
 *
 * Builds the native xTB basis on the CPU, takes the flattened atomic numbers +
 * Bohr geometry, resolves the per-atom EEQ parameters via the shared
 * D4ChargeModel::resolveParams, and:
 *   (1) solves the (N+1) augmented EEQ system on the device (XtbGpuContext::
 *       eeqCharges) and compares the atomic charges to D4ChargeModel CPU;
 *   (2) feeds a deterministic synthetic dE/dq into the device adjoint pair-loop
 *       (eeqChargeResponseGradient) and compares the Cartesian charge-response
 *       gradient to D4ChargeModel::addChargeResponseGradient.
 *
 * The LU solve is FP64; the charges match to ~1e-10. The response has an O(N)
 * reduction (CN sum) + an O(N) per-atom pair loop summed in a different order
 * than the CPU's i<j loop, so the gate is 1e-9 to absorb reassociation.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_eeq [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.h"
#include "src/core/energy_calculators/dispersion/d4_charge_model.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"

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
using curcuma::dispersion::D4ChargeModel;

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
        std::cerr << "usage: test_xtb_cuda_eeq [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_eeq: no usable CUDA device\n";
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
        xtb.Calculation(false);

        GpuBasisFlat bf;
        GpuH0Flat hf;
        xtb.exportGpuBasis(bf, hf);

        // Geometry in Bohr (the flattened basis already carries it) as an N×3 Matrix
        // for the CPU reference and a contiguous [3i+k] array for the device.
        Matrix geom_bohr(nat, 3);
        std::vector<double> xyz_flat(3 * nat);
        for (int i = 0; i < nat; ++i)
            for (int k = 0; k < 3; ++k) {
                geom_bohr(i, k) = bf.xyz_bohr[3 * i + k];
                xyz_flat[3 * i + k] = bf.xyz_bohr[3 * i + k];
            }

        // CPU reference: charges + a deterministic synthetic dE/dq response.
        D4ChargeModel eeq;
        Vector q_cpu = eeq.computeCharges(atoms, geom_bohr, 0.0);
        Vector dEdq(nat);
        for (int i = 0; i < nat; ++i) dEdq(i) = 0.01 * std::sin(0.7 * i + 1.0);
        Matrix grad_cpu = Matrix::Zero(nat, 3);
        eeq.addChargeResponseGradient(dEdq, grad_cpu);

        // Device path: identical params via the shared resolver.
        std::vector<double> chi, gam, alp, cnf, rcov;
        D4ChargeModel::resolveParams(atoms, chi, gam, alp, cnf, rcov);
        std::vector<double> q_gpu(nat, 0.0), dedq_h(nat), grad_gpu(3 * nat, 0.0);
        for (int i = 0; i < nat; ++i) dedq_h[i] = dEdq(i);

        if (!ctx.eeqCharges(nat, xyz_flat.data(), chi.data(), gam.data(), alp.data(),
                            cnf.data(), rcov.data(), 0.0, q_gpu.data())) {
            std::cerr << xyz << ": device eeqCharges failed\n"; ++failures; continue;
        }
        if (!ctx.eeqChargeResponseGradient(nat, dedq_h.data(), grad_gpu.data())) {
            std::cerr << xyz << ": device eeqChargeResponseGradient failed\n"; ++failures; continue;
        }

        double dq = 0.0, dgrad = 0.0;
        for (int i = 0; i < nat; ++i) dq = std::max(dq, std::fabs(q_cpu(i) - q_gpu[i]));
        for (int a = 0; a < nat; ++a)
            for (int k = 0; k < 3; ++k)
                dgrad = std::max(dgrad, std::fabs(grad_cpu(a, k) - grad_gpu[3 * a + k]));

        const bool ok = (dq <= tol) && (dgrad <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok) {
            std::printf("%s %-28s %s: max|dq|=%.3e max|dgrad|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(), dq, dgrad, tol);
        }
    }
    return failures == 0 ? 0 : 1;
}
