/*
 * <Phase 3.2 H0 validation — native xTB bare Hamiltonian vs tblite dump>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reads one Phase-0 dump, rebuilds the xTB basis and coordination number,
 * assembles H0 per tblite's xtb/h0.f90 formula:
 *
 *   rr_ij = sqrt( sqrt(r²) / (rad_i + rad_j) )
 *   π_ij  = (1 + shpoly_i · rr_ij) · (1 + shpoly_j · rr_ij)
 *   ε_i   = p_selfenergy(i) − p_kcn(i) · CN_atom(i)
 *   H_μν  = S_μν · 0.5·(ε_i + ε_j) · hscale_{ij} · π_ij    (i ≠ j)
 *   H_μν  = S_μν · 0.5·(ε_i + ε_j)                         (i = j, on-atom)
 *
 * hscale_{ij} differs between GFN1 (kpair, valence-aware, s<->p=2.08) and
 * GFN2 (kpair=1, no valence, Slater-ratio prefactor zij = (2√(ζ_iζ_j)/(ζ_i+ζ_j))^0.5).
 *
 *   test_xtb_h0 <dump.json>
 *
 * Exit 0 if max |ΔH| < tol (1e-7 Eh); non-zero otherwise.
 *
 * Claude Generated (Phase 3.2, Apr 2026). GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/parameters/gfn1_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn2_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/xtb_params_extra.hpp"
#include "src/core/energy_calculators/qm_methods/STO_CGTO.hpp"

#include "external/json.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

struct ShellRef {
    int atom;
    int ang;         // angular momentum ℓ
    int ish_local;   // index within atom's shell list
    int ao_start;
    int nao;
    bool valence;    // GFN1 only; GFN2 ignores this
    CGTO::Shell cg;
};

// Look up parameters from the method-specific generated table.
struct MethodTables {
    const int*       nshell;         // [86]
    const int      (*ang_shell)[3];  // [86][3]
    const int      (*principalN)[3]; // [86][3]
    const double   (*zeta)[3];       // [86][3]  slater_exponent
    const int      (*nprim)[3];      // [86][3]
    const double   (*selfenergy)[3]; // [86][3]  (already Hartree)
    const double   (*kcn)[3];        // [86][3]  (already Hartree/CN)
    const double   (*shpoly)[3];     // [86][3]  indexed by angular momentum
};

MethodTables tablesFor(const std::string& method)
{
    if (method == "gfn1") {
        return {
            curcuma::xtb::gfn1_params::nshell,
            curcuma::xtb::gfn1_params::ang_shell,
            curcuma::xtb::gfn1_params::principal_quantum_number,
            curcuma::xtb::gfn1_params::slater_exponent,
            curcuma::xtb::gfn1_params::number_of_primitives,
            curcuma::xtb::gfn1_params::p_selfenergy,
            curcuma::xtb::gfn1_params::p_kcn,
            curcuma::xtb::gfn1_params::p_shpoly,
        };
    }
    return {
        curcuma::xtb::gfn2_params::nshell,
        curcuma::xtb::gfn2_params::ang_shell,
        curcuma::xtb::gfn2_params::principal_quantum_number,
        curcuma::xtb::gfn2_params::slater_exponent,
        curcuma::xtb::gfn2_params::number_of_primitives,
        curcuma::xtb::gfn2_params::p_selfenergy,
        curcuma::xtb::gfn2_params::p_kcn,
        curcuma::xtb::gfn2_params::p_shpoly,
    };
}

std::vector<ShellRef> buildBasis(const std::vector<int>& z_atoms,
                                  const MethodTables&    t)
{
    std::vector<ShellRef> shells;
    int ao_cursor = 0;

    for (size_t iat = 0; iat < z_atoms.size(); ++iat) {
        const int z  = z_atoms[iat];
        const int ns = t.nshell[z - 1];

        std::vector<CGTO::Shell> cgs(ns);
        for (int ish = 0; ish < ns; ++ish) {
            cgs[ish] = CGTO::slater_to_gauss(
                t.principalN[z - 1][ish],
                t.ang_shell[z - 1][ish],
                t.zeta[z - 1][ish],
                t.nprim[z - 1][ish]);
        }
        for (int ish = 1; ish < ns; ++ish) {
            for (int jsh = 0; jsh < ish; ++jsh) {
                if (cgs[ish].ang == cgs[jsh].ang) {
                    CGTO::orthogonalize(cgs[jsh], cgs[ish]);
                }
            }
        }

        // Valence flag: first shell of each ℓ is valence (GFN1 rule).
        int ang_seen[3] = {0, 0, 0};
        for (int ish = 0; ish < ns; ++ish) {
            ShellRef sr{};
            sr.atom       = static_cast<int>(iat);
            sr.ang        = cgs[ish].ang;
            sr.ish_local  = ish;
            sr.ao_start   = ao_cursor;
            sr.nao        = 2 * cgs[ish].ang + 1;
            const int l   = cgs[ish].ang;
            sr.valence    = (ang_seen[l] == 0);
            if (sr.valence) ang_seen[l] = 1;
            sr.cg         = cgs[ish];
            shells.push_back(std::move(sr));
            ao_cursor += sr.nao;
        }
    }
    return shells;
}

void aoToType(const std::vector<ShellRef>& shells, int ao,
              int& ish_out, int& type_out)
{
    for (size_t i = 0; i < shells.size(); ++i) {
        const auto& sr = shells[i];
        if (ao >= sr.ao_start && ao < sr.ao_start + sr.nao) {
            ish_out = static_cast<int>(i);
            const int local = ao - sr.ao_start;
            if (sr.ang == 0)       type_out = 0;
            else if (sr.ang == 1)  { static const int p_map[3] = {2, 3, 1}; type_out = p_map[local]; }
            else                   type_out = -1;
            return;
        }
    }
    ish_out = -1; type_out = -1;
}

// GFN1 hscale: valence-aware, enscale negative, kpair with hardcoded specials.
double hscale_gfn1(int zi, int zj, int il, int jl, bool vi, bool vj)
{
    using namespace curcuma::xtb;
    if (vi && vj) {
        const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
        const double enp = 1.0 + gfn1::enscale * den;
        return gfn1::kpair(zi, zj) * gfn1::kshell(il, jl) * enp;
    } else if (vi && !vj) {
        return 0.5 * (gfn1::kshell(il, il) + gfn1::kdiff);
    } else if (!vi && vj) {
        return 0.5 * (gfn1::kshell(jl, jl) + gfn1::kdiff);
    }
    return gfn1::kdiff;
}

// GFN2 hscale: no valence, kpair=1, extra Slater-ratio prefactor zij.
double hscale_gfn2(int zi, int zj, int il, int jl, double zeta_i, double zeta_j)
{
    using namespace curcuma::xtb;
    const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
    const double enp = 1.0 + gfn2::enscale * den;
    const double km  = gfn2::kpair(zi, zj) * gfn2::kshell(il, jl) * enp;
    const double zij = std::pow(2.0 * std::sqrt(zeta_i * zeta_j) / (zeta_i + zeta_j),
                                gfn2::wexp);
    return zij * km;
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
    json d; in >> d;

    const std::string method = d.at("method");
    const int nat = d.at("natoms");
    const int nao = d.at("norbitals");

    std::vector<int>    z_atoms;
    std::vector<double> xyz_flat(3 * nat);
    z_atoms.reserve(nat);
    for (int i = 0; i < nat; ++i) {
        z_atoms.push_back(d["atoms"][i]["z"]);
        xyz_flat[3*i + 0] = d["atoms"][i]["xyz_bohr"][0];
        xyz_flat[3*i + 1] = d["atoms"][i]["xyz_bohr"][1];
        xyz_flat[3*i + 2] = d["atoms"][i]["xyz_bohr"][2];
    }

    const MethodTables t = tablesFor(method);
    auto shells = buildBasis(z_atoms, t);
    int my_nao = 0;
    for (const auto& sr : shells) my_nao += sr.nao;
    if (my_nao != nao) {
        std::fprintf(stderr, "AO count mismatch: curcuma=%d  dump=%d\n", my_nao, nao);
        return 4;
    }

    // Overlap and reference Hamiltonian from dump.
    std::vector<double> Sref(static_cast<size_t>(nao) * nao);
    std::vector<double> Href(static_cast<size_t>(nao) * nao);
    for (int i = 0; i < nao; ++i) {
        for (int j = 0; j < nao; ++j) {
            Sref[i * nao + j] = d["overlap"][i][j];
            Href[i * nao + j] = d["hamiltonian"][i][j];
        }
    }

    // Coordination number.
    std::vector<double> cn = (method == "gfn1")
        ? curcuma::xtb::cn_exp(z_atoms, xyz_flat)
        : curcuma::xtb::cn_gfn(z_atoms, xyz_flat);

    // Per-shell CN-adjusted self energy, indexed globally by shell number.
    std::vector<double> eps_shell(shells.size());
    for (size_t ish = 0; ish < shells.size(); ++ish) {
        const int iat = shells[ish].atom;
        const int z   = z_atoms[iat];
        const int li  = shells[ish].ish_local;
        eps_shell[ish] = t.selfenergy[z - 1][li] - t.kcn[z - 1][li] * cn[iat];
    }

    // Assemble H0 element-wise.
    std::vector<double> H(static_cast<size_t>(nao) * nao, 0.0);
    for (int mu = 0; mu < nao; ++mu) {
        int ish_a = -1, type_a = -1;
        aoToType(shells, mu, ish_a, type_a);
        const auto& sa = shells[ish_a];
        for (int nu = 0; nu < nao; ++nu) {
            int ish_b = -1, type_b = -1;
            aoToType(shells, nu, ish_b, type_b);
            const auto& sb = shells[ish_b];
            if (type_a < 0 || type_b < 0) continue;
            const int ia = sa.atom, ib = sb.atom;
            const int zi = z_atoms[ia], zj = z_atoms[ib];

            const double avg_eps = 0.5 * (eps_shell[ish_a] + eps_shell[ish_b]);
            double hij;
            if (ia == ib) {
                hij = avg_eps;  // on-atom: shpoly=1, no hscale factor
            } else {
                const double dx = xyz_flat[3*ia+0] - xyz_flat[3*ib+0];
                const double dy = xyz_flat[3*ia+1] - xyz_flat[3*ib+1];
                const double dz = xyz_flat[3*ia+2] - xyz_flat[3*ib+2];
                const double r2 = dx*dx + dy*dy + dz*dz;
                const double rr = std::sqrt(std::sqrt(r2)
                                 / (curcuma::xtb::atomic_rad_au(zi)
                                    + curcuma::xtb::atomic_rad_au(zj)));
                const double shp_i = t.shpoly[zi - 1][sa.ang];
                const double shp_j = t.shpoly[zj - 1][sb.ang];
                const double pi_ij = (1.0 + shp_i * rr) * (1.0 + shp_j * rr);

                double hs;
                if (method == "gfn1") {
                    hs = hscale_gfn1(zi, zj, sa.ang, sb.ang,
                                     sa.valence, sb.valence);
                } else {
                    const double zeta_i = t.zeta[zi - 1][sa.ish_local];
                    const double zeta_j = t.zeta[zj - 1][sb.ish_local];
                    hs = hscale_gfn2(zi, zj, sa.ang, sb.ang, zeta_i, zeta_j);
                }
                hij = avg_eps * hs * pi_ij;
            }
            H[mu * nao + nu] = Sref[mu * nao + nu] * hij;
        }
    }

    // Diff.
    double max_abs = 0.0, rms = 0.0;
    int imax = 0, jmax = 0;
    for (int i = 0; i < nao; ++i) {
        for (int j = 0; j < nao; ++j) {
            const double diff = H[i * nao + j] - Href[i * nao + j];
            if (std::fabs(diff) > max_abs) {
                max_abs = std::fabs(diff);
                imax = i; jmax = j;
            }
            rms += diff * diff;
        }
    }
    rms = std::sqrt(rms / (nao * nao));

    std::printf("%s  nao=%d  max|ΔH|=%.3e at (%d,%d)  H_me=%.10f  H_ref=%.10f  rms=%.3e\n",
                argv[1], nao, max_abs, imax, jmax,
                H[imax * nao + jmax], Href[imax * nao + jmax], rms);

    const double tol = 1e-7;
    return (max_abs < tol) ? 0 : 1;
}
