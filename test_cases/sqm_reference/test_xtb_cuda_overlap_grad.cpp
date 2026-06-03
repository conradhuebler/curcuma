/*
 * Stage 4 component test: device overlap derivative dS/dR vs the CPU reference.
 *
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-06, GPU port Stage 4).
 *
 * The overlap derivative dS_μν/dR_{atom(μ)} (Obara-Saika) is the crux primitive of
 * the nuclear-gradient H0/Pulay term. This test validates the device port
 * (k_overlap_grad → d_cgto_overlap_grad) against the CPU CGTO::cgto_overlap_grad
 * elementwise over every AO pair, before the full gradient is assembled. Pairwise
 * (no reductions) → bit-level agreement expected; gate 1e-9 for pow/exp ULP.
 *
 * Without a usable CUDA device the test prints SKIP and exits 0.
 *
 * Usage:  test_xtb_cuda_overlap_grad [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"
#include "src/core/energy_calculators/qm_methods/STO_CGTO.hpp"
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

// tblite AO ordering, mirrors ao_to_type in xtb_h0.cpp.
int ao_to_type(int ang, int local_ao)
{
    if (ang == 0) return 0;
    if (ang == 1) { static const int p_map[3] = {2, 3, 1}; return p_map[local_ao]; }
    return -1;
}

// Reconstruct a CGTO::Shell for shell `s` from the flattened basis.
CGTO::Shell shellOf(const GpuBasisFlat& bf, int s)
{
    CGTO::Shell sh;
    sh.ang   = bf.ang_sh[s];
    sh.nprim = bf.sh_nprim[s];
    const int off = bf.sh_prim_off[s];
    sh.alpha.assign(bf.prim_alpha.begin() + off, bf.prim_alpha.begin() + off + sh.nprim);
    sh.coeff.assign(bf.prim_coeff.begin() + off, bf.prim_coeff.begin() + off + sh.nprim);
    return sh;
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
        std::cerr << "usage: test_xtb_cuda_overlap_grad [--quiet] [--tol V] <gfn1|gfn2> <xyz> [<xyz> ...]\n";
        return 2;
    }
    const MethodType mt = (method == "gfn1") ? MethodType::GFN1 : MethodType::GFN2;

    XtbGpuContext ctx;
    if (!ctx.ok()) {
        if (!quiet) std::cout << "SKIP test_xtb_cuda_overlap_grad: no usable CUDA device\n";
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
            std::cerr << xyz << ": InitialiseMolecule failed\n"; ++failures; continue;
        }
        xtb.Calculation(false);

        GpuBasisFlat bf; GpuH0Flat hf;
        xtb.exportGpuBasis(bf, hf);
        const int nao = bf.nao;

        auto bd = makeGpuBasisData(bf, hf);
        const size_t nn = static_cast<size_t>(nao) * static_cast<size_t>(nao);
        std::vector<double> dSdR(3 * nn);
        if (!ctx.beginBasis(bd) || !ctx.computeOverlapGrad(bf.xyz_bohr.data(), dSdR.data())) {
            std::cerr << xyz << ": device overlap-grad failed\n"; ++failures; continue;
        }

        // CPU reference: CGTO::cgto_overlap_grad per AO pair (dS/dR of the row atom).
        double dmax = 0.0;
        for (int mu = 0; mu < nao; ++mu) {
            const int isha = bf.ao2sh[mu], iat = bf.ao2at[mu];
            const int ta = ao_to_type(bf.ang_sh[isha], mu - bf.iao_sh[isha]);
            const CGTO::Shell shA = shellOf(bf, isha);
            for (int nu = 0; nu < nao; ++nu) {
                const int ishb = bf.ao2sh[nu], jat = bf.ao2at[nu];
                const int tb = ao_to_type(bf.ang_sh[ishb], nu - bf.iao_sh[ishb]);
                double gcpu[3] = {0, 0, 0};
                if (ta >= 0 && tb >= 0) {
                    const CGTO::Shell shB = shellOf(bf, ishb);
                    CGTO::cgto_overlap_grad(shA, shB,
                        bf.xyz_bohr[3*iat+0], bf.xyz_bohr[3*iat+1], bf.xyz_bohr[3*iat+2],
                        bf.xyz_bohr[3*jat+0], bf.xyz_bohr[3*jat+1], bf.xyz_bohr[3*jat+2],
                        ta, tb, gcpu);
                }
                const size_t mn = static_cast<size_t>(mu) + static_cast<size_t>(nu) * nao;
                for (int k = 0; k < 3; ++k)
                    dmax = std::max(dmax, std::fabs(dSdR[k * nn + mn] - gcpu[k]));
            }
        }

        const bool ok = (dmax <= tol);
        if (!ok) ++failures;
        if (!quiet || !ok)
            std::printf("%s %-28s %s: max|d(dS/dR)|=%.3e tol=%.1e\n",
                        ok ? "PASS" : "FAIL", xyz.c_str(), method.c_str(), dmax, tol);
    }
    return failures == 0 ? 0 : 1;
}
