/*
 * gCP and SRB of HF-3c for every element H-Rn against otool_gcp (ORCA 6.1's gCP program).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * reference_gcp.json (scripts/hf3c_gcp_reference.py): per structure E_gCP and E_SRB (10 decimals, as
 * otool_gcp prints them), the gradient of gCP alone (-level hf/minix -grad) and of gCP + SRB (-level hf3c
 * -grad), Eh/Bohr, and whether otool's gradient is the derivative of its own energy. Checks gcp::energy /
 * gcp::baseEnergy with the HF-3c parameters:
 *   |dE_gCP|, |dE_SRB| <= 2e-8 Eh. otool reads its tables from default-REAL (single precision) literals
 *     (emiss, D3 radii) and converts with autoang 0.5291772083; curcuma keeps the double tables of
 *     simple-dftd3 (the H-Ne oracle, unchanged), measured gap <= 1.2e-8 Eh (Ar cluster).
 *   gradient = 4-point central difference of curcuma's own energy (h = 1e-3 Bohr) to 1e-8 Eh/Bohr, every structure;
 *   gradient vs otool to 1e-8 Eh/Bohr only where otool_grad_consistent (17 of 90 structures are not: pairs with
 *     He/Ne, which have no virtual MINIX functions, and pairs of near-equal Slater exponents, see the script).
 * Claude Generated (Nov 2026). AI-generated, machine-tested. GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/gcp.h"
#include "src/core/units.h"

#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>

using json = nlohmann::json;

static int elementNumber(const std::string& s)
{
    static const char* sym[] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
        "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
        "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
        "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
        "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn" };
    for (int i = 0; i < 86; ++i)
        if (s == sym[i]) return i + 1;
    return 0;
}

static double maxDiff(const Matrix& g, const json& ref)
{
    double d = 0.0;
    for (int i = 0; i < g.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            d = std::max(d, std::abs(g(i, k) - ref[i][k].get<double>()));
    return d;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::printf("usage: test_gcp_elements reference_gcp.json\n");
        return 2;
    }
    std::ifstream f(argv[1]);
    const json ref = json::parse(f);
    const gcp::Parameters p3c = gcp::hf3c();
    gcp::Parameters pgcp = p3c;
    pgcp.base = false;

    int fails = 0, n = 0;
    double worst_e = 0.0, worst_g = 0.0, worst_fd = 0.0;
    for (const auto& [name, s] : ref["structures"].items()) {
        std::vector<int> atoms;
        Matrix xyz(s["atoms"].size(), 3);
        int i = 0;
        for (const auto& a : s["atoms"]) {
            atoms.push_back(elementNumber(a[0].get<std::string>()));
            for (int k = 0; k < 3; ++k)
                xyz(i, k) = a[k + 1].get<double>() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
            ++i;
        }
        Matrix g_gcp, g_tot;
        const double e_gcp = gcp::energy(atoms, xyz, pgcp, &g_gcp);
        const double e_tot = gcp::energy(atoms, xyz, p3c, &g_tot);
        const double e_srb = gcp::baseEnergy(atoms, xyz, p3c);
        const double de_gcp = std::abs(e_gcp - s["gcp"].get<double>());
        const double de_srb = std::abs(e_srb - s["srb"].get<double>());
        const double de_sum = std::abs(e_tot - e_gcp - e_srb);
        const double dg_gcp = maxDiff(g_gcp, s["grad_gcp"]);
        const double dg_tot = maxDiff(g_tot, s["grad_total"]);
        // 4-point central difference of the gCP+SRB energy, Bohr
        double dg_fd = 0.0;
        const double h = 1.0e-3;  // Bohr; smaller steps see the cancellation noise of the closed-form B_k (mix_tm: 1e-7 at 2e-5)
        for (int a = 0; a < xyz.rows(); ++a)
            for (int k = 0; k < 3; ++k) {
                auto e = [&](double d) { Matrix q = xyz; q(a, k) += d; return gcp::energy(atoms, q, p3c); };
                const double fd = (8.0 * (e(h) - e(-h)) - (e(2.0 * h) - e(-2.0 * h))) / (12.0 * h);  // 4-point
                dg_fd = std::max(dg_fd, std::abs(fd - g_tot(a, k)));
            }
        const bool consistent = s["otool_grad_consistent"].get<bool>();
        const bool ok = de_gcp <= 2e-8 && de_srb <= 2e-8 && de_sum <= 1e-14 && dg_fd <= 1e-8
            && (!consistent || (dg_gcp <= 1e-8 && dg_tot <= 1e-8));
        worst_e = std::max({ worst_e, de_gcp, de_srb });
        worst_fd = std::max(worst_fd, dg_fd);
        if (consistent) worst_g = std::max({ worst_g, dg_gcp, dg_tot });
        if (!ok) ++fails;
        ++n;
        std::printf("%-16s gCP %.10f (d %.1e)  SRB %.10f (d %.1e)  |dG| vs FD %.1e, vs otool gCP %.1e total %.1e%s  %s\n",
                    name.c_str(), e_gcp, de_gcp, e_srb, de_srb, dg_fd, dg_gcp, dg_tot,
                    consistent ? "" : " (otool gradient not its energy's derivative, not compared)", ok ? "ok" : "FAIL");
    }
    std::printf("%d structures, %d failed; max |dE| %.1e Eh, max |dG - FD| %.1e, max |dG - otool| %.1e Eh/Bohr "
                "(consistent structures)\n", n, fails, worst_e, worst_fd, worst_g);
    return fails == 0 ? 0 : 1;
}
