/*
 * Stage 6 (S6.1) component test: device orbital occupations vs the host
 * occupationsFromEps.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 6).
 *
 * Runs a bare CPU XTB SCF to obtain a realistic ascending eigenvalue spectrum,
 * then fills the occupations two ways at the SAME (eps, electronic temperature,
 * electron count): (a) the authoritative host logic XTB::computeOccupations (=
 * occupationsFromEps — integer closed-shell fill at Tele=0, else a Fermi-level
 * bisection) and (b) the device single-block kernel XtbGpuContext::occupations
 * (k_occupations, the port of occupationsFromEps). Both consume the same eps, so
 * the comparison isolates the kernel from any SCF difference. Pure arithmetic →
 * the occupations agree to ~1e-12 (gate 1e-9); ncol (the last column with
 * occ>1e-12) must match exactly.
 *
 * Method-independent kernel (gfn1 + gfn2 only differ in the spectrum). Without a
 * usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_occ [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
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
    if (method.empty() || xyz_paths.empty()) {
        std::cerr << "usage: test_xtb_cuda_occ [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_occ: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        // CPU SCF to get a realistic ascending eps spectrum + the temperature /
        // electron count the host occupation logic uses.
        XTB xtb(mt);
        if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);
        const Vector eps   = xtb.getOrbitalEnergies();
        const double Tele  = xtb.electronicTemperature();
        const double nelec = static_cast<double>(xtb.getNumElectrons());
        const int n = static_cast<int>(eps.size());
        if (n <= 0) { std::cerr << xyz << ": empty spectrum\n"; ++failures; continue; }

        // Host reference: the SCF's own occupation logic.
        Eigen::VectorXd occ_host;
        int ncol_host = 0;
        xtb.computeOccupations(eps, occ_host, ncol_host);

        // Device kernel at the same eps / Tele / nelec.
        std::vector<double> occ_gpu(n, 0.0);
        int    ncol_gpu = 0;
        double mu_gpu   = 0.0;
        if (!ctx.occupations(eps.data(), n, Tele, nelec, occ_gpu.data(), &ncol_gpu, &mu_gpu)) {
            std::cerr << xyz << ": device occupations failed\n"; ++failures; continue;
        }
        if (occ_host.size() != n) {
            std::cerr << xyz << ": host occ size mismatch\n"; ++failures; continue;
        }

        double docc = 0.0;
        for (int i = 0; i < n; ++i) docc = std::max(docc, std::fabs(occ_host(i) - occ_gpu[i]));
        const bool ok = (docc <= tol) && (ncol_host == ncol_gpu);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s %s: max|d_occ|=%.3e ncol(h/g)=%d/%d tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(),
                        docc, ncol_host, ncol_gpu, tol);
    }
    return failures == 0 ? 0 : 1;
}
