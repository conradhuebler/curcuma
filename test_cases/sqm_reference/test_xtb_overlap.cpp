/*
 * <Phase 3.1 overlap validation — native xTB vs tblite dump>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reads one Phase-0 dump, rebuilds the xTB basis locally using
 * curcuma::xtb::gfn1_params / gfn2_params + CGTO::slater_to_gauss
 * (+ orthogonalize for same-ℓ shells), computes the overlap S, and
 * diffs against the dump's S matrix element-wise.
 *
 *   test_xtb_overlap <dump.json>
 *
 * Exit 0 if max |ΔS| < tol (1e-8); non-zero otherwise with a report.
 * The numeric gate for the next kernel (H0) — there is no point fitting
 * H0 if the overlap under it is wrong.
 *
 * Claude Generated (Phase 3.1, Apr 2026)
 * GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/parameters/gfn1_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn2_params.hpp"
#include "src/core/energy_calculators/qm_methods/STO_CGTO.hpp"
#include "src/core/energy_calculators/qm_methods/xtb_ao_utils.hpp"  // X-I1: d block

#include "external/json.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

struct ShellRef {
    int atom;      // atom index this shell lives on
    int ang;       // 0=s, 1=p
    int ao_start;  // first AO index in global AO ordering
    int nao;       // AOs in this shell (1 for s, 3 for p)
    CGTO::Shell cg;
};

// Build the basis for a molecule using curcuma's extracted xTB parameters.
// Matches tblite's basis_type layout: flat list of shells in atom-major order,
// with same-ℓ shells on the same atom orthogonalized Gram-Schmidt.
std::vector<ShellRef> buildBasis(const std::vector<int>& z_atoms,
                                  const std::string&     method)
{
    std::vector<ShellRef> shells;
    int ao_cursor = 0;

    auto nshell = [&](int z) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::nshell[z - 1]
                 : curcuma::xtb::gfn2_params::nshell[z - 1];
    };
    auto angShell = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::ang_shell[z - 1][ish]
                 : curcuma::xtb::gfn2_params::ang_shell[z - 1][ish];
    };
    auto principalN = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::principal_quantum_number[z - 1][ish]
                 : curcuma::xtb::gfn2_params::principal_quantum_number[z - 1][ish];
    };
    auto zeta = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::slater_exponent[z - 1][ish]
                 : curcuma::xtb::gfn2_params::slater_exponent[z - 1][ish];
    };
    auto nprim = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::number_of_primitives[z - 1][ish]
                 : curcuma::xtb::gfn2_params::number_of_primitives[z - 1][ish];
    };

    for (size_t iat = 0; iat < z_atoms.size(); ++iat) {
        const int z  = z_atoms[iat];
        const int ns = nshell(z);

        // Build CGTOs first per-shell, then orthogonalize pairs sharing ℓ.
        std::vector<CGTO::Shell> cgs(ns);
        for (int ish = 0; ish < ns; ++ish) {
            cgs[ish] = CGTO::slater_to_gauss(
                principalN(z, ish),
                angShell(z, ish),
                zeta(z, ish),
                nprim(z, ish));
        }
        for (int ish = 1; ish < ns; ++ish) {
            for (int jsh = 0; jsh < ish; ++jsh) {
                if (cgs[ish].ang == cgs[jsh].ang) {
                    CGTO::orthogonalize(cgs[jsh], cgs[ish]);
                }
            }
        }
        for (int ish = 0; ish < ns; ++ish) {
            ShellRef sr{};
            sr.atom     = static_cast<int>(iat);
            sr.ang      = cgs[ish].ang;
            sr.ao_start = ao_cursor;
            sr.nao      = 2 * cgs[ish].ang + 1;
            sr.cg       = cgs[ish];
            shells.push_back(std::move(sr));
            ao_cursor += sr.nao;
        }
    }
    return shells;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <dump.json>\n", argv[0]);
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 3; }
    json d;
    in >> d;

    const std::string method = d.at("method");
    const int nat = d.at("natoms");
    const int nao = d.at("norbitals");

    std::vector<int> z_atoms;
    z_atoms.reserve(nat);
    std::vector<std::array<double,3>> xyz(nat);
    for (int i = 0; i < nat; ++i) {
        z_atoms.push_back(d["atoms"][i]["z"]);
        xyz[i][0] = d["atoms"][i]["xyz_bohr"][0];
        xyz[i][1] = d["atoms"][i]["xyz_bohr"][1];
        xyz[i][2] = d["atoms"][i]["xyz_bohr"][2];
    }

    auto shells = buildBasis(z_atoms, method);
    int my_nao = 0;
    for (const auto& sr : shells) my_nao += sr.nao;
    if (my_nao != nao) {
        std::fprintf(stderr,
                     "AO count mismatch: curcuma=%d  dump=%d — basis skeleton diverges\n",
                     my_nao, nao);
        return 4;
    }

    // Build the reference overlap matrix from the dump.
    std::vector<double> Sref(static_cast<size_t>(nao) * nao);
    for (int i = 0; i < nao; ++i)
        for (int j = 0; j < nao; ++j)
            Sref[i * nao + j] = d["overlap"][i][j];

    // Compute S via shell-pair blocks, mirroring xtb_h0.cpp: pure s/p pairs use
    // the scalar CGTO::cgto_overlap path, d-touching pairs use the production
    // curcuma::xtb::sphericalOverlapBlock (cartesian->spherical dtrafo, X-I1).
    std::vector<double> S(static_cast<size_t>(nao) * nao, 0.0);
    for (size_t ai = 0; ai < shells.size(); ++ai) {
        const auto& sa = shells[ai];
        const int A = sa.atom;
        for (size_t bi = 0; bi < shells.size(); ++bi) {
            const auto& sb = shells[bi];
            const int B = sb.atom;
            if (sa.ang < 2 && sb.ang < 2) {
                for (int ia = 0; ia < sa.nao; ++ia) {
                    const int t_a = curcuma::xtb::ao_to_type(sa.ang, ia);
                    for (int jb = 0; jb < sb.nao; ++jb) {
                        const int t_b = curcuma::xtb::ao_to_type(sb.ang, jb);
                        const double s = (A == B && ai == bi && t_a == t_b)
                            ? 1.0
                            : CGTO::cgto_overlap(sa.cg, sb.cg,
                                                 xyz[A][0], xyz[A][1], xyz[A][2],
                                                 xyz[B][0], xyz[B][1], xyz[B][2],
                                                 t_a, t_b);
                        S[(sa.ao_start + ia) * nao + (sb.ao_start + jb)] = s;
                    }
                }
            } else {
                double blk[6 * 6];
                curcuma::xtb::sphericalOverlapBlock(
                    sa.cg, sa.ang, sb.cg, sb.ang,
                    xyz[A][0], xyz[A][1], xyz[A][2],
                    xyz[B][0], xyz[B][1], xyz[B][2], blk, 6);
                for (int ia = 0; ia < sa.nao; ++ia)
                    for (int jb = 0; jb < sb.nao; ++jb) {
                        double s = blk[ia * 6 + jb];
                        if (A == B && ai == bi && ia == jb) s = 1.0;
                        S[(sa.ao_start + ia) * nao + (sb.ao_start + jb)] = s;
                    }
            }
        }
    }

    // Diff.
    double max_abs = 0.0, rms = 0.0;
    int    imax = 0, jmax = 0;
    for (int i = 0; i < nao; ++i) {
        for (int j = 0; j < nao; ++j) {
            const double diff = S[i * nao + j] - Sref[i * nao + j];
            if (std::fabs(diff) > max_abs) {
                max_abs = std::fabs(diff);
                imax = i; jmax = j;
            }
            rms += diff * diff;
        }
    }
    rms = std::sqrt(rms / (nao * nao));

    std::printf("%s  nao=%d  max|ΔS|=%.3e at (%d,%d)  S_me=%.10f  S_ref=%.10f  rms=%.3e\n",
                argv[1], nao, max_abs, imax, jmax,
                S[imax * nao + jmax], Sref[imax * nao + jmax], rms);

    const double tol = 1e-8;
    return (max_abs < tol) ? 0 : 1;
}
