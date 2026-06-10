/*
 * Stage 6 (S6.3) component test: device SCC energy components vs the host energy
 * helpers, at the identical converged resident state.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 6).
 *
 * Drives a full device-resident GFN2 SCF. The fused loop leaves the converged
 * shell/atom charges + multipole moments resident, so the device SCC energy
 * (E_coulomb = 1/2 q_sh^T gamma q_sh, the GFN2 shell third-order, the GFN2
 * multipole energy) is read straight off the resident state via GpuScfBackend::
 * sccEnergy. The host reference is the tblite-validated XTB::
 * evaluateComponentsAtFixedDensity recomputed at the SAME state (the downloaded
 * density / charges + the device-downloaded moments) — a true frozen-state
 * device-vs-host comparison at identical charges, so the device matches the host
 * to ~1e-10 (gate 1e-9), NOT a CPU-vs-GPU SCF fixpoint difference.
 *
 * Note: after a fully-resident GFN2 run the public getE*() members still hold the
 * device per-iteration values, so the host reference must be re-derived with
 * evaluateComponentsAtFixedDensity (which sets m_E_coulomb_shell / m_E_third_order
 * / m_E_multipole on the host) rather than read directly off the GPU run.
 *
 * GFN2 only (the resident loop + device SCC energy are the GFN2 device path).
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_scc_energy [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]
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
        std::cerr << "usage: test_xtb_cuda_scc_energy [--quiet] [--tol V] gfn2 <xyz> [<xyz> ...]\n";
        return 2;
    }

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_scc_energy: no usable CUDA device\n";
        return 0;
    }

    int failures = 0;
    for (const std::string& xyz : xyz_paths) {
        std::vector<int> atoms;
        std::vector<double> coords_ang;
        if (!readXYZ(xyz, atoms, coords_ang)) { ++failures; continue; }
        const int nat = static_cast<int>(atoms.size());

        // Device-resident GFN2 SCF; converged charges/moments stay resident.
        auto backend = createXtbGpuScfBackend(&ctx);
        XTB xtb(MethodType::GFN2);
        xtb.setGpuScfBackend(backend.get());
        if (!xtb.InitialiseMolecule(atoms.data(), coords_ang.data(), nat, 0.0, 0)) {
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);

        // Converged state from the host wavefunction + the device-resident moments.
        const Matrix P    = xtb.getDensity();
        const Vector q_at = xtb.getCharges();
        const Vector q_sh = xtb.getShellCharges();
        const int nsh = static_cast<int>(q_sh.size());
        Eigen::MatrixXd dp_at, qp_at;
        if (!backend->multipoleMoments(dp_at, qp_at)) {
            std::cerr << xyz << ": device multipoleMoments failed (resident path inactive?)\n";
            ++failures; continue;
        }

        // Device SCC energy from the resident state (capture before the host rebuild).
        double ec_gpu = 0.0, e3_gpu = 0.0, emp_gpu = 0.0;
        if (!backend->sccEnergy(nat, nsh, ec_gpu, e3_gpu, emp_gpu)) {
            std::cerr << xyz << ": device sccEnergy failed\n"; ++failures; continue;
        }

        // Host reference: recompute the components at the SAME state. This sets the
        // host m_E_coulomb_shell / m_E_third_order / m_E_multipole, which the public
        // getters then report (the resident run left device values in them).
        if (!xtb.evaluateComponentsAtFixedDensity(P, q_at, q_sh, dp_at, qp_at)) {
            std::cerr << xyz << ": host evaluateComponentsAtFixedDensity failed\n"; ++failures; continue;
        }
        const double ec_host  = xtb.getECoulombShell();
        const double e3_host   = xtb.getEThirdOrder();
        const double emp_host  = xtb.getEMultipole();

        const double de_c = std::fabs(ec_gpu  - ec_host);
        const double de_3 = std::fabs(e3_gpu  - e3_host);
        const double de_m = std::fabs(emp_gpu - emp_host);
        const double dmax = std::max(de_c, std::max(de_3, de_m));
        const bool ok = (dmax <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s gfn2: max|dE_scc|=%.3e (coul=%.2e third=%.2e mp=%.2e) tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), dmax, de_c, de_3, de_m, tol);
    }
    return failures == 0 ? 0 : 1;
}
