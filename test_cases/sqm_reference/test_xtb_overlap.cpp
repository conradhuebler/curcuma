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

// Map a global AO index to (shell_index, type_in_STO_CGTO).
// type encoding: 0=s, 1=px, 2=py, 3=pz.  Only s and p handled so far.
void aoToType(const std::vector<ShellRef>& shells, int ao, int& ish_out, int& type_out)
{
    for (size_t i = 0; i < shells.size(); ++i) {
        const auto& sr = shells[i];
        if (ao >= sr.ao_start && ao < sr.ao_start + sr.nao) {
            ish_out = static_cast<int>(i);
            const int local = ao - sr.ao_start;
            if (sr.ang == 0) {
                type_out = 0;   // s
            } else if (sr.ang == 1) {
                // tblite p-shell AO ordering is (py, pz, px) — spherical m=(-1,0,+1),
                // see external/tblite/src/tblite/integral/multipole.f90 lx table (comment
                // "x (+1), y (-1), z (0) in [-1, 0, 1] sorting").
                static const int p_map[3] = { 2, 3, 1 };   // py, pz, px in STO_CGTO types
                type_out = p_map[local];
            } else {
                type_out = -1;  // d not handled
            }
            return;
        }
    }
    ish_out  = -1;
    type_out = -1;
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

    // Compute S via CGTO; one s/p pair at a time.
    std::vector<double> S(static_cast<size_t>(nao) * nao, 0.0);
    for (int mu = 0; mu < nao; ++mu) {
        int ish_a = -1, type_a = -1;
        aoToType(shells, mu, ish_a, type_a);
        const auto& sa = shells[ish_a];
        for (int nu = 0; nu < nao; ++nu) {
            int ish_b = -1, type_b = -1;
            aoToType(shells, nu, ish_b, type_b);
            const auto& sb = shells[ish_b];
            if (type_a < 0 || type_b < 0) continue; // d unsupported here
            const int ia = sa.atom, ib = sb.atom;
            const double s = CGTO::cgto_overlap(
                sa.cg, sb.cg,
                xyz[ia][0], xyz[ia][1], xyz[ia][2],
                xyz[ib][0], xyz[ib][1], xyz[ib][2],
                type_a, type_b);
            S[mu * nao + nu] = s;
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
