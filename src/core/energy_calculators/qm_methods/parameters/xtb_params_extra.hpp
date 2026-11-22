/* -------------------------------------------------------------------------- *
 *  xTB auxiliary parameters and generators needed by the bare Hamiltonian.   *
 *                                                                            *
 *  Tables ported verbatim from tblite:                                       *
 *   - external/tblite/src/tblite/data/covrad.f90      (covalent_rad_2009)    *
 *   - external/tblite/src/tblite/data/paulingen.f90   (pauling_en)           *
 *   - external/tblite/src/tblite/data/atomicrad.f90   (atomic_rad)           *
 *   - external/tblite/src/tblite/xtb/gfn1.f90         (kshell, kpair)        *
 *   - external/tblite/src/tblite/xtb/gfn2.f90         (kshell, wexp)         *
 *   - external/tblite/src/tblite/ncoord/exp.f90       (GFN1 CN, kcn=16)      *
 *   - external/tblite/src/tblite/ncoord/gfn.f90       (GFN2 CN, ka=10 kb=20) *
 *                                                                            *
 *  Values are kept in atomic units where tblite stores them in AU, and in    *
 *  Angstrom where tblite does the aatoau multiplication at table level —     *
 *  the `aatoau()` helper below applies the conversion explicitly.            *
 *                                                                            *
 *  Claude Generated (Phase 3.2, Apr 2026). GPL-3.0.                          *
 * -------------------------------------------------------------------------- */
#pragma once
#ifndef CURCUMA_XTB_PARAMS_EXTRA_H_
#define CURCUMA_XTB_PARAMS_EXTRA_H_

#include "src/core/element_tables.h"
#include "src/core/units.h"
#include <cmath>
#include <vector>

namespace curcuma::xtb {

inline constexpr double AA_TO_AU = CurcumaUnit::Length::ANGSTROM_TO_BOHR;

// ---- Pauling electronegativities (dimensionless, elements 1..86) ----------
// paulingen.f90 lines 40-57. Elements 87+ are dummy 1.50.
inline constexpr const double (&pauling_en)[86] = CurcumaElements::PAULING_EN;  // table lives in src/core/element_tables.h

// ---- Atomic radii (CRC Handbook / Mantina), elements 1..86, in Angstrom ---
// atomicrad.f90 lines 46-61 (AU there; we store in Å and convert via AA_TO_AU).
inline constexpr const double (&atomic_rad_angstrom)[86] = CurcumaElements::ATOMIC_RAD_ANGSTROM;  // table lives in src/core/element_tables.h

// ---- Pyykkö-Atsumi covalent radii (metals -10%), elements 1..86, in Å -----
// covrad.f90 lines 43-60. D3 form is (4/3) * this (also applied below).
inline constexpr const double (&covalent_rad_2009_angstrom)[86] = CurcumaElements::COVALENT_RAD_2009_ANGSTROM;  // table lives in src/core/element_tables.h

inline double atomic_rad_au(int z) {
    return atomic_rad_angstrom[z - 1] * AA_TO_AU;
}
inline double covalent_rad_d3_au(int z) {
    return (4.0 / 3.0) * covalent_rad_2009_angstrom[z - 1] * AA_TO_AU;
}

// ---- GFN1 hscale constants (gfn1.f90 lines 56-58) -------------------------
namespace gfn1 {
    inline constexpr double kdiag[5] = {1.85, 2.25, 2.00, 2.00, 2.00};
    inline constexpr double enscale  = -7.0e-3;
    inline constexpr double kdiff    = 2.85;

    // kshell(jl, il) — gfn1.f90 line 1040: 2.08 for s<->p cross, else average.
    inline constexpr double kshell(int il, int jl) {
        const bool sp_cross = (il == 0 && jl == 1) || (il == 1 && jl == 0);
        return sp_cross ? 2.08 : 0.5 * (kdiag[il] + kdiag[jl]);
    }

    // d-block row (gfn1.f90 lines 756-770).
    inline int dblock_row(int z) {
        if (z > 20 && z < 30) return 1;
        if (z > 38 && z < 48) return 2;
        if (z > 56 && z < 80) return 3;
        return 0;
    }

    // get_pair_param (gfn1.f90 lines 723-754): hardcoded special cases, then
    // d-block average, else 1.0.
    inline double kpair(int izp, int jzp) {
        if (izp == 1 && jzp == 1) return 0.96;
        auto is_pair = [&](int a, int b) {
            return (izp == a && jzp == b) || (izp == b && jzp == a);
        };
        if (is_pair(5, 1))  return 0.95;
        if (is_pair(7, 1))  return 1.04;
        if (is_pair(28, 1)) return 0.90;
        if (is_pair(75, 1)) return 0.80;
        if (is_pair(78, 1)) return 0.80;
        if (is_pair(15, 5)) return 0.97;
        if (is_pair(14, 7)) return 1.01;
        const int itr = dblock_row(izp);
        const int jtr = dblock_row(jzp);
        if (itr > 0 && jtr > 0) {
            static constexpr double kp[3] = {1.1, 1.2, 1.2};
            return 0.5 * (kp[itr - 1] + kp[jtr - 1]);
        }
        return 1.0;
    }
} // namespace gfn1

// ---- GFN2 hscale constants (gfn2.f90 lines 60-61) -------------------------
namespace gfn2 {
    inline constexpr double kdiag[5] = {1.85, 2.23, 2.23, 2.23, 2.23};
    inline constexpr double enscale  = 2.0e-2;
    inline constexpr double wexp     = 0.5;

    // kshell — gfn2.f90 line 1005: 2.0 for p<->d, s<->d cross, else average.
    inline constexpr double kshell(int il, int jl) {
        const bool sd_pd_cross = (jl == 2 && (il == 0 || il == 1))
                              || (il == 2 && (jl == 0 || jl == 1));
        return sd_pd_cross ? 2.0 : 0.5 * (kdiag[il] + kdiag[jl]);
    }

    // GFN2 has no special kpair pairs — kpair(:, :) = 1.0 (gfn2.f90 line 727).
    inline constexpr double kpair(int /*izp*/, int /*jzp*/) { return 1.0; }
} // namespace gfn2

// ---- Coordination number (exp counting for GFN1) --------------------------
// exp.f90 line 103: count = 1/(1 + exp(-kcn*(rc/r - 1))), default kcn=16,
// rc = rcov_d3(i) + rcov_d3(j).
inline std::vector<double>
cn_exp(const std::vector<int>& z_atoms,
       const std::vector<double>& xyz_bohr,
       double kcn = 16.0, double cutoff_au = 25.0)
{
    const int nat = static_cast<int>(z_atoms.size());
    std::vector<double> cn(nat, 0.0);
    std::vector<double> rcov(nat);
    for (int i = 0; i < nat; ++i) rcov[i] = covalent_rad_d3_au(z_atoms[i]);
    const double cutoff2 = cutoff_au * cutoff_au;
    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j < i; ++j) {
            const double dx = xyz_bohr[3*i + 0] - xyz_bohr[3*j + 0];
            const double dy = xyz_bohr[3*i + 1] - xyz_bohr[3*j + 1];
            const double dz = xyz_bohr[3*i + 2] - xyz_bohr[3*j + 2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > cutoff2 || r2 < 1.0e-12) continue;
            const double r  = std::sqrt(r2);
            const double rc = rcov[i] + rcov[j];
            const double c  = 1.0 / (1.0 + std::exp(-kcn * (rc / r - 1.0)));
            cn[i] += c;
            cn[j] += c;   // directed_factor = 1.0
        }
    }
    return cn;
}

// ---- Coordination number (double-exp "gfn" counting for GFN2) -------------
// gfn.f90 line 99: count = gfn(ka=10, r, rc) * gfn(kb=20, r, rc + r_shift=2.0)
inline std::vector<double>
cn_gfn(const std::vector<int>& z_atoms,
       const std::vector<double>& xyz_bohr,
       double cutoff_au = 25.0)
{
    constexpr double ka = 10.0, kb = 20.0, r_shift = 2.0;
    auto gcount = [](double k, double r, double r0) {
        return 1.0 / (1.0 + std::exp(-k * (r0 / r - 1.0)));
    };
    const int nat = static_cast<int>(z_atoms.size());
    std::vector<double> cn(nat, 0.0);
    std::vector<double> rcov(nat);
    for (int i = 0; i < nat; ++i) rcov[i] = covalent_rad_d3_au(z_atoms[i]);
    const double cutoff2 = cutoff_au * cutoff_au;
    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j < i; ++j) {
            const double dx = xyz_bohr[3*i + 0] - xyz_bohr[3*j + 0];
            const double dy = xyz_bohr[3*i + 1] - xyz_bohr[3*j + 1];
            const double dz = xyz_bohr[3*i + 2] - xyz_bohr[3*j + 2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > cutoff2 || r2 < 1.0e-12) continue;
            const double r  = std::sqrt(r2);
            const double rc = rcov[i] + rcov[j];
            const double c  = gcount(ka, r, rc) * gcount(kb, r, rc + r_shift);
            cn[i] += c;
            cn[j] += c;
        }
    }
    return cn;
}

} // namespace curcuma::xtb

#endif
