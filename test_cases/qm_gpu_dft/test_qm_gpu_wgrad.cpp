/*
 * <Becke weight derivative (grid response) on the device vs the CPU for any number of atoms>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026). The device kernel of cuda/qm_gpu_xc.cu (QMGpuXC::weightGradient) had a limit of 700
 * atoms (shared memory); this test compares it with qmgrid::gridWeightGradient on the structure of --xyz (ctest:
 * test_cases/structures/bulk/peo201-chain.xyz, 1410 atoms, small grid), with a synthetic function f on the grid (smooth, both signs, every 7th point zero to exercise the skip). Only the Becke partition
 * and the grid enter; the AO tables handed to the device are one s function per atom (not used by this check).
 *
 *   test_qm_gpu_wgrad --xyz <file> --grid sg0 [--partition becke|stratmann] [--threads 8] [--tol 1e-12] [--stride 7]
 *
 * Prints one line "WGRAD natoms=... npts=... max|dg|=... max|g|=... rel=... cpu_ms=... dev_ms=..." and the device
 * statistics; exit 0 when rel = max|g_dev - g_cpu| / max|g_cpu| <= tol, 1 otherwise, 77 without a CUDA device/plugin.
 *
 * This program is free software under GPL-3.0
 */

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/energy_calculators/qm_methods/qm_gpu_xc_backend.h"
#include "src/core/energy_calculators/qm_methods/qm_grid.h"
#include "src/core/energy_calculators/qm_methods/qm_xc_kernels.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace {

/// XYZ file (Angstrom) -> atomic numbers and Bohr positions; elements H-Ne. Claude Generated.
bool readXYZ(const std::string& path, std::vector<int>& Z, Matrix& pos)
{
    std::ifstream in(path);
    int n = 0;
    std::string line;
    if (!(in >> n) || n <= 0) return false;
    std::getline(in, line);
    std::getline(in, line);
    static const char* sym[] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne" };
    pos.resize(n, 3);
    for (int a = 0; a < n; ++a) {
        std::string e;
        double x, y, z;
        if (!(in >> e >> x >> y >> z)) return false;
        int zz = 0;
        for (int k = 0; k < 10; ++k)
            if (e == sym[k]) zz = k + 1;
        if (zz == 0) return false;
        Z.push_back(zz);
        const double b = 1.0 / 0.529177210903;  // CODATA 2018 Bohr radius in Angstrom, as src/core/units.h
        pos(a, 0) = x * b;
        pos(a, 1) = y * b;
        pos(a, 2) = z * b;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv)
{
    int threads = 8, stride = 7;
    std::string xyz, level = "sg0", partition = "becke";
    double tol = 1e-12;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string k = argv[i], v = argv[i + 1];
        if (k == "--xyz") xyz = v;
        else if (k == "--grid") level = v;
        else if (k == "--partition") partition = v;
        else if (k == "--threads") threads = std::atoi(v.c_str());
        else if (k == "--tol") tol = std::atof(v.c_str());
        else if (k == "--stride") stride = std::max(2, std::atoi(v.c_str()));
    }
    CurcumaLogger::set_verbosity(0);
    std::vector<int> Z;
    Matrix pos;
    if (!readXYZ(xyz, Z, pos)) {
        std::printf("cannot read %s (xyz, elements H-Ne)\n", xyz.c_str());
        return 1;
    }
    const int na = static_cast<int>(Z.size());
    qmgrid::GridOptions opt;
    opt.level = level;
    opt.partition = partition;
    qmgrid::MolecularGrid grid;
    std::string why;
    const auto t0 = std::chrono::steady_clock::now();
    if (!qmgrid::buildMolecularGrid(Z, pos, opt, grid, why, threads)) {
        std::printf("grid: %s\n", why.c_str());
        return 1;
    }
    const double grid_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const int npts = grid.size();

    // one s function per atom (the device object needs AO tables; the weight derivative does not use them)
    std::vector<GTO::Orbital> basis(na);
    for (int a = 0; a < na; ++a) {
        basis[a].x = pos(a, 0);
        basis[a].y = pos(a, 1);
        basis[a].z = pos(a, 2);
        basis[a].type = GTO::S;
        basis[a].exponents = { 1.0 };
        basis[a].coefficients = { 1.0 };
        basis[a].atom = a;
    }
    const qmgrid::AOEvaluator ao(basis, Matrix());
    qmgrid::FlatShells fs;
    Matrix inv_dist, a_adj;
    bool stratmann = false;
    if (!ao.flatten(fs, why) || !qmgrid::partitionTables(Z, pos, opt, inv_dist, a_adj, stratmann, why)) {
        std::printf("setup: %s\n", why.c_str());
        return 1;
    }
    std::vector<int> sh_atom(fs.nshell);
    for (int s = 0; s < fs.nshell; ++s) sh_atom[s] = basis[fs.cart_first[s]].atom;
    const int kernel = static_cast<int>(qmxc::XCKernel::SlaterX);
    const double weight = 1.0;
    qmint::QMGpuXCSetup su;
    su.nbf = ao.nbf();
    su.npts = npts;
    su.natoms = na;
    su.points = grid.points.data();
    su.weights = grid.weights.data();
    su.quad_weights = grid.quadrature_weights.data();
    su.point_atom = grid.atom.data();
    su.atom_offset = grid.atom_offset.data();
    su.nshell = fs.nshell;
    su.max_ncomp = fs.max_ncomp;
    su.nexps = (int)fs.exps.size();
    su.ncoef = (int)fs.coef.size();
    su.nlmn = (int)fs.lmn.size();
    su.nT = (int)fs.T.size();
    su.ncomp_total = (int)fs.coef_off.size();
    su.sh_xyz = fs.xyz.data();
    su.sh_amin = fs.amin.data();
    su.sh_nprim = fs.nprim.data();
    su.sh_prim_off = fs.prim_off.data();
    su.sh_ncomp = fs.ncomp.data();
    su.sh_comp_off = fs.comp_off.data();
    su.sh_out_first = fs.out_first.data();
    su.sh_out_width = fs.out_width.data();
    su.sh_T_off = fs.T_off.data();
    su.sh_atom = sh_atom.data();
    su.exps = fs.exps.data();
    su.coef = fs.coef.data();
    su.coef_off = fs.coef_off.data();
    su.lmn = fs.lmn.data();
    su.T = nullptr;
    su.nterms = 1;
    su.term_kernel = &kernel;
    su.term_weight = &weight;
    su.rho_min = 1e-12;
    su.atom_xyz = pos.data();
    su.inv_dist = inv_dist.data();
    su.a_adjust = a_adj.data();
    su.stratmann = stratmann ? 1 : 0;
    su.block_points = 1024;  // small AO blocks: the test needs little device memory besides the weight scratch
    if (!gpu_plugin::available("cuda")) {
        std::printf("no CUDA plugin: skipped\n");
        return 77;
    }
    auto dev = gpu_plugin::createQmXC("cuda", su, 0, why);
    if (!dev) {
        std::printf("no device XC object (%s): skipped\n", why.c_str());
        return 77;
    }

    // synthetic f: smooth in space, both signs; zero on every point with p % stride == 0 (stride 7: every 7th point
    // skipped; a large stride keeps only every stride-th point and so shortens the O(natoms^2) per point CPU reference)
    Vector f(npts);
    for (int p = 0; p < npts; ++p) {
        const int A = grid.atom[p];
        const double d = (grid.points.row(p) - pos.row(A)).norm();
        const bool keep = stride <= 7 ? (p % stride != 0) : (p % stride == 1);
        f(p) = !keep ? 0.0 : std::exp(-0.5 * d) * (0.3 + std::sin(1.7 * grid.points(p, 0) + 0.4 * grid.points(p, 2)));
    }
    Matrix gc = Matrix::Zero(na, 3), gd = Matrix::Zero(na, 3);
    const auto c0 = std::chrono::steady_clock::now();
    if (!qmgrid::gridWeightGradient(Z, pos, opt, grid, f, gc, why, threads)) {
        std::printf("CPU weight derivative: %s\n", why.c_str());
        return 1;
    }
    const auto c1 = std::chrono::steady_clock::now();
    if (!dev->weightGradient(f.data(), gd.data())) {
        std::printf("device weight derivative failed\n");
        return 1;
    }
    const auto c2 = std::chrono::steady_clock::now();
    const double dmax = (gd - gc).cwiseAbs().maxCoeff(), gmax = gc.cwiseAbs().maxCoeff();
    const double rel = dmax / std::max(gmax, 1e-300);
    const bool ok = rel <= tol;
    long nonzero = 0;
    for (int p = 0; p < npts; ++p) nonzero += f(p) != 0.0;
    std::printf("WGRAD nonzero_f=%ld ", nonzero);
    std::printf("natoms=%d npts=%d grid=%s partition=%s max|dg|=%.3e max|g|=%.3e rel=%.3e (tol %.0e) grid_ms=%.1f "
                "cpu_ms=%.1f (%d threads) dev_ms=%.1f %s\n",
                na, npts, level.c_str(), partition.c_str(), dmax, gmax, rel, tol, grid_ms,
                std::chrono::duration<double, std::milli>(c1 - c0).count(), threads,
                std::chrono::duration<double, std::milli>(c2 - c1).count(), ok ? "ok" : "FAIL");
    std::printf("device: %s\n", dev->stats().c_str());
    return ok ? 0 : 1;
}
