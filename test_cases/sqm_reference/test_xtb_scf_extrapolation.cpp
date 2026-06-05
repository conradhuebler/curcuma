/*
 * Multi-step SCC extrapolation validation — native xTB (GFN1/GFN2).
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Drives the native XTB SCF along a smooth synthetic trajectory (linear
 * interpolation between a start geometry and its optimised geometry, N+1 frames)
 * and compares three SCC-guess strategies across geometry steps:
 *
 *   none  : the existing 1-step warm-start (reuse only the previous step)
 *   aspc  : Kolafa ASPC predictor over the last order+2 converged steps
 *   gauss : least-squares polynomial extrapolation of the SCC-vector history
 *
 * For each method (GFN1, GFN2) it asserts:
 *   (A) Correctness — the final-frame energy of aspc/gauss equals none within a
 *       tight gate (same SCF fixpoint; the predictor only changes the guess).
 *   (B) Speed-up — the cumulative SCF iteration count of aspc and gauss is
 *       <= none for this smooth path (the whole point of the feature).
 *
 * Both molecule geometries come from companion xyz files (no hardcoded
 * geometry, per the test-suite rule):
 *
 *   test_xtb_scf_extrapolation <start.xyz> <opt.xyz> [nsteps]
 *
 * Claude Generated (SQM SCC extrapolation WP, June 2026). GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/xtb_native.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using curcuma::xtb::MethodType;
using curcuma::xtb::XTB;

namespace {

int Z_of(const std::string& s)
{
    static const std::map<std::string, int> m = {
        { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 }, { "C", 6 },
        { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 }, { "Na", 11 }, { "Mg", 12 },
        { "Al", 13 }, { "Si", 14 }, { "P", 15 }, { "S", 16 }, { "Cl", 17 }, { "Ar", 18 },
        { "K", 19 }, { "Ca", 20 }, { "Br", 35 }, { "I", 53 }
    };
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
    std::getline(f, line); // rest of first line
    std::getline(f, line); // comment line
    atoms.clear();
    coords_ang.assign(3 * nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        std::string sym;
        double x, y, z;
        if (!(f >> sym >> x >> y >> z)) {
            std::cerr << "xyz parse error at atom " << i << "\n";
            return false;
        }
        int zn = Z_of(sym);
        if (zn < 0) { std::cerr << "unknown element symbol: " << sym << "\n"; return false; }
        atoms.push_back(zn);
        coords_ang[3 * i + 0] = x;
        coords_ang[3 * i + 1] = y;
        coords_ang[3 * i + 2] = z;
    }
    return true;
}

struct Result {
    double energy = 0.0;
    long   iters  = 0;
    bool   ok     = false;
};

// Drive the SCF along the interpolated trajectory with the given extrapolation
// strategy and return the final-frame energy + cumulative SCF iterations.
Result run(MethodType mt, std::vector<int> atoms,
           const std::vector<double>& a, const std::vector<double>& b,
           int nsteps, const std::string& mode, int order,
           const std::string& apply = "guess", int correctors = 1)
{
    Result r;
    const int nat = static_cast<int>(atoms.size());
    XTB xtb(mt);
    xtb.setScfThreshold(1.0e-8);  // tight: the converged energies are comparable to ~1e-8
    xtb.setWarmStart(true);
    xtb.setScfExtrapolation(mode);
    xtb.setScfExtrapolationOrder(order);
    xtb.setScfExtrapolationApply(apply);
    xtb.setScfXlbomdCorrectors(correctors);

    if (!xtb.InitialiseMolecule(atoms.data(), a.data(), nat, 0.0, 0)) {
        std::cerr << "InitialiseMolecule failed\n";
        return r;
    }

    for (int s = 0; s <= nsteps; ++s) {
        const double t = static_cast<double>(s) / static_cast<double>(nsteps);
        if (s > 0) {
            Matrix geom = Matrix::Zero(nat, 3);
            for (int i = 0; i < nat; ++i)
                for (int k = 0; k < 3; ++k)
                    geom(i, k) = (1.0 - t) * a[3 * i + k] + t * b[3 * i + k];
            if (!xtb.UpdateMolecule(geom)) {
                std::cerr << "UpdateMolecule failed at step " << s << "\n";
                return r;
            }
        }
        r.energy = xtb.Calculation(false);
        r.iters += xtb.scfIterations();
    }
    r.ok = std::isfinite(r.energy);
    return r;
}

bool checkMethod(MethodType mt, const char* name,
                 const std::vector<int>& atoms,
                 const std::vector<double>& a, const std::vector<double>& b,
                 int nsteps, int order)
{
    const double E_GATE = 1.0e-6; // same SCF fixpoint (loose-ball safety on a tight 1e-8 SCF)
    Result none  = run(mt, atoms, a, b, nsteps, "none",  order);
    Result aspc  = run(mt, atoms, a, b, nsteps, "aspc",  order);
    Result gauss = run(mt, atoms, a, b, nsteps, "gauss", order);

    if (!none.ok || !aspc.ok || !gauss.ok) {
        std::printf("[%s] FAIL: a calculation did not finish\n", name);
        return false;
    }

    const double dE_aspc  = std::fabs(aspc.energy  - none.energy);
    const double dE_gauss = std::fabs(gauss.energy - none.energy);

    std::printf("[%s] final E: none=%.10f  aspc=%.10f (dE=%.2e)  gauss=%.10f (dE=%.2e)\n",
                name, none.energy, aspc.energy, dE_aspc, gauss.energy, dE_gauss);
    std::printf("[%s] SCF iters over %d steps: none=%ld  aspc=%ld (%+ld)  gauss=%ld (%+ld)\n",
                name, nsteps + 1, none.iters,
                aspc.iters, aspc.iters - none.iters,
                gauss.iters, gauss.iters - none.iters);

    bool ok = true;
    if (dE_aspc > E_GATE) {
        std::printf("[%s] FAIL: aspc energy differs by %.2e > %.0e\n", name, dE_aspc, E_GATE);
        ok = false;
    }
    if (dE_gauss > E_GATE) {
        std::printf("[%s] FAIL: gauss energy differs by %.2e > %.0e\n", name, dE_gauss, E_GATE);
        ok = false;
    }
    // Speed-up: on a smooth trajectory the multi-step predictor must not cost
    // more iterations than the 1-step warm-start.
    if (aspc.iters > none.iters) {
        std::printf("[%s] FAIL: aspc used MORE iterations (%ld > %ld)\n", name, aspc.iters, none.iters);
        ok = false;
    }
    if (gauss.iters > none.iters) {
        std::printf("[%s] FAIL: gauss used MORE iterations (%ld > %ld)\n", name, gauss.iters, none.iters);
        ok = false;
    }
    if (ok)
        std::printf("[%s] PASS\n", name);
    return ok;
}

// XL-BOMD (apply=xlbomd) check: the time-reversibly propagated auxiliary density seeds
// the SCF, which still converges, so the energy must (a) stay finite, (b) equal the
// fully-converged (none) energy, and (c) cost fewer SCF iterations (good guess).
bool checkXlbomd(MethodType mt, const char* name,
                 const std::vector<int>& atoms,
                 const std::vector<double>& a, const std::vector<double>& b,
                 int nsteps, int order)
{
    const double E_TRACK = 1.0e-6; // corrector converges -> same fixpoint as none
    Result none = run(mt, atoms, a, b, nsteps, "none",  order);
    Result xl   = run(mt, atoms, a, b, nsteps, "aspc",  order, "xlbomd", 1);

    if (!none.ok || !xl.ok) {
        std::printf("[%s/xlbomd] FAIL: a calculation did not finish (non-finite energy)\n", name);
        return false;
    }
    const double dE = std::fabs(xl.energy - none.energy);
    std::printf("[%s/xlbomd] final E: none=%.10f  xlbomd=%.10f (dE=%.2e)\n",
                name, none.energy, xl.energy, dE);
    std::printf("[%s/xlbomd] SCF iters over %d steps: none=%ld  xlbomd=%ld (%+ld)\n",
                name, nsteps + 1, none.iters, xl.iters, xl.iters - none.iters);

    bool ok = true;
    if (dE > E_TRACK) {
        std::printf("[%s/xlbomd] FAIL: energy tracking %.2e > %.0e\n", name, dE, E_TRACK);
        ok = false;
    }
    if (xl.iters >= none.iters) {
        std::printf("[%s/xlbomd] FAIL: xlbomd did not reduce iterations (%ld >= %ld)\n",
                    name, xl.iters, none.iters);
        ok = false;
    }
    if (ok)
        std::printf("[%s/xlbomd] PASS\n", name);
    return ok;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <start.xyz> <opt.xyz> [nsteps]\n";
        return 2;
    }
    const std::string start_path = argv[1];
    const std::string opt_path   = argv[2];
    const int nsteps = (argc > 3) ? std::atoi(argv[3]) : 14;
    const int order  = 3;

    std::vector<int> a_atoms, b_atoms;
    std::vector<double> a_xyz, b_xyz;
    if (!readXYZ(start_path, a_atoms, a_xyz)) return 2;
    if (!readXYZ(opt_path, b_atoms, b_xyz)) return 2;
    if (a_atoms != b_atoms || a_xyz.size() != b_xyz.size()) {
        std::cerr << "start/opt xyz mismatch (atom order or count)\n";
        return 2;
    }

    bool ok = true;
    ok &= checkMethod(MethodType::GFN1, "gfn1", a_atoms, a_xyz, b_xyz, nsteps, order);
    ok &= checkMethod(MethodType::GFN2, "gfn2", a_atoms, a_xyz, b_xyz, nsteps, order);
    // XL-BOMD (Phase 2, experimental) corrector-only coupling.
    ok &= checkXlbomd(MethodType::GFN1, "gfn1", a_atoms, a_xyz, b_xyz, nsteps, order);
    ok &= checkXlbomd(MethodType::GFN2, "gfn2", a_atoms, a_xyz, b_xyz, nsteps, order);

    std::printf("%s\n", ok ? "ALL PASS" : "FAILURES");
    return ok ? 0 : 1;
}
