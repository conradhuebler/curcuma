/*
 * <rev-gfnff P3 "harris": non-self-consistent excess-electron energy correction g(r)>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (Sep 24, 2026) - test_cases/revgfnff/_log/P2P3_HARRIS_STATUS.md.
 * AI-FITTED, machine-tested only; human production testing pending.
 *
 * HAND-MAINTAINED. Not written by scripts/revgfnff_wellfit.py (it lives in its own header so no
 * generator can overwrite it).
 *
 * rev_excess_mode harris leaves the split charges of a perceived 2c-3e pair (Cl2-, F2-) free, so
 * the Coulomb term keeps the EEQ delocalisation energy of that pair - which has the wrong size and
 * the wrong r-trend (P2P3_STATUS.md section 1). This table takes it out again with an additive
 * energy that is a function of the pair distance and of the topological excess-electron count x
 * only:
 *
 *     E_harris = x_ij * g(r_ij),        dE_harris/dr = x_ij * g'(r_ij)
 *
 * x_ij comes from GFNFF::revExcessElectrons() (a topology constant, never the self-consistent
 * charges), so E_harris is never fed back into the charge solve ("Harris-like") and the
 * gradient has no charge or x cross term.
 *
 * g(r) = A - B exp(-c r)   (kcal/mol, r in Angstrom; bounded: g -> A at long r), fitted to
 *     g(r_k) = E_DLPNO-CCSD(T)(r_k) - E_model(r_k; rev_excess_mode flat, rev_excess_kappa 0)
 * i.e. the correction that makes the harris total equal the reference at the fit points. Energies
 * relative to X + X-, DLPNO-CCSD(T)/aug-cc-pVTZ grids ref/E/{cl2m_Cl-Cl-,f2m_F-F-}_dlpno_ccsdt,
 * model flags -gfnff.rev_charge_model sqe -gfnff.rev_sqe_phase1 true -gfnff.rev_excess_electron
 * true (ALL kappa_Z = 0 - g is only valid there: with kappa_Z > 0 the charges localise by
 * themselves and g over-corrects, measured +128 / +225 kcal/mol in react mode). Two point sets:
 *   (1) every static grid point whose topology contains the X-X bond (Cl2- 1.52-2.73 A, n = 11;
 *       F2- 1.44-2.02 A, n = 7), fresh single points, weight 3;
 *   (2) the react-mode breaking scan (-gfnff.topology_mode react, 0.05 A steps) up to the point
 *       where the pair's bond order falls below rev_sqe_bmin (Cl2- to 3.80 A, F2- to 2.50 A),
 *       weight 1 - without them the curve is extrapolated 40 % past the static range and a
 *       linear g over-corrects the react tail by up to +17 kcal/mol.
 * A per c scan with linear least squares in A, B. Measured (P2P3_HARRIS_STATUS.md): static
 * bonded rms 2.57 / 2.93 kcal/mol (leave-one-out 3.02 / 3.52), react breaking rms 1.79 / 2.59.
 * Linear g on set (1) alone: static 2.23 / 2.07, react breaking 5.35 / 4.88.
 */

#pragma once

#include <cmath>
#include <cstddef>

namespace RevHarrisTable {

struct HarrisEntry {
    int z1, z2;     ///< element pair, z1 <= z2
    double A;       ///< kcal/mol (long-range limit of g)
    double B;       ///< kcal/mol
    double c;       ///< 1/Angstrom
};

// ---- BEGIN hand-maintained block ------------------------------------------------------------
inline constexpr HarrisEntry kHarrisEntries[] = {
    { 17, 17, 115.8267243370, 120.0728512838, 0.7429824561 },   // Cl-Cl (Cl2-, DLPNO-CCSD(T))
    { 35, 35, 115.4623886836, 99.6617053206, 0.6835839599 },   // Br-Br (Br2-, DLPNO-CCSD(T), static bonded rms 1.69 LOO 2.14, react break 1.82) X2BR:harris
    { 53, 53, 89.3328189456, 52.6698342651, 0.4756892231 },   // I-I (I2-, DLPNO-CCSD(T), static bonded rms 2.12 LOO 2.61, react break rms 2.02; refit on the final half row, I2_CLF_STATUS 12) X2I:harris
    { 9, 9, 212.2064031706, 238.6608379001, 1.3765664160 },     // F-F   (F2-,  DLPNO-CCSD(T))
    { 9, 17, -11.8546099202, -198.5314773496, 0.0500000000 },    // Cl-F (ClF-, DLPNO-CCSD(T), static bonded rms 8.40 LOO 9.52, react break 6.94; refit on the final half row, c at grid bound = linear limit, I2_CLF_STATUS 12) X2CLF:harris
};
// ---- END hand-maintained block --------------------------------------------------------------
inline constexpr std::size_t kHarrisCount = sizeof(kHarrisEntries) / sizeof(kHarrisEntries[0]);

/// g(r) in Eh for r in Bohr, and dg/dr in Eh/Bohr. Returns false (g = dg = 0) for a pair
/// without a row - such a pair then carries no correction.
inline bool harrisG(int za, int zb, double r_bohr, double& g_eh, double& dgdr_eh_per_bohr)
{
    const int z1 = za < zb ? za : zb;
    const int z2 = za < zb ? zb : za;
    constexpr double kcal = 627.5094740631;       // kcal/mol per Eh (CODATA-2018, src/core/units.h)
    constexpr double bohr_to_ang = 0.529177210903; // CODATA-2018, src/core/units.h
    for (std::size_t k = 0; k < kHarrisCount; ++k) {
        if (kHarrisEntries[k].z1 != z1 || kHarrisEntries[k].z2 != z2)
            continue;
        const HarrisEntry& e = kHarrisEntries[k];
        const double ex = std::exp(-e.c * r_bohr * bohr_to_ang);
        g_eh = (e.A - e.B * ex) / kcal;
        dgdr_eh_per_bohr = e.B * e.c * ex * bohr_to_ang / kcal;   // d/dr of -B exp(-c r)
        return true;
    }
    g_eh = 0.0;
    dgdr_eh_per_bohr = 0.0;
    return false;
}

} // namespace RevHarrisTable
