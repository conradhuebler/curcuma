/*
 * <QMDFF element reference data for Curcuma>
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
 * =============================================================================
 * ACKNOWLEDGMENT OF ORIGINAL WORK
 * =============================================================================
 *
 * The tables in this file are transcribed from the QMDFF force field of
 *   S. Grimme, "A General Quantum Mechanically Derived Force Field for
 *   Structures, Vibrational Frequencies and Non-Covalent Interaction Energies
 *   of Large Molecules", J. Chem. Theory Comput. 2014, 10, 4497-4514.
 *   DOI: 10.1021/ct500573f
 *
 * Fortran reference implementation (removed from xtb in commit b7dbd36,
 * "Refactoring of external drivers (#568)"):
 *   xtb/src/qmdff.f90        -- module xtb_qmdff (ff_ini, valelff, ...)
 *   xtb/src/disp/dftd3.f     -- setr0ab (D3 cutoff radii)
 *   xtb/src/param/atomicrad.f90 -- atomicRad (Mantina/Truhlar radii)
 *
 * Claude Generated (August 2026): C++ transcription of the element tables that
 * the QMDFF energy expression needs. Values are byte-for-byte the Fortran ones.
 */

#pragma once

namespace QMDFFData {

/// Highest element supported by the QMDFF tables (Fortran: max_elem = 94).
inline constexpr int kMaxElem = 94;

/**
 * @brief D3 cutoff radius R0(A,B) in Bohr.
 *
 * Reference: xtb/src/disp/dftd3.f, subroutine setr0ab (packed lower triangle,
 * 4465 = 94*95/2 entries, tabulated in Angstrom and divided by autoang=0.52917726).
 *
 * Used by QMDFF only to set the repulsion range parameter
 *   alpha(A,B) = 16.5 / R0(A,B)^1.5        (qmdff.f90:144)
 *
 * @param zi Atomic number of A (1-based, 1..94)
 * @param zj Atomic number of B (1-based, 1..94)
 * @return R0 in Bohr; 0.0 for out-of-range elements
 */
double r0ab(int zi, int zj);

/**
 * @brief QMDFF repulsion range parameter alpha(A,B) = 16.5 / R0(A,B)^1.5.
 *
 * Reference: qmdff.f90:143-144
 * @return alpha in 1/Bohr; 0.0 for out-of-range elements
 */
double repulsionAlpha(int zi, int zj);

/**
 * @brief Mantina/Truhlar atomic radius in Angstrom.
 *
 * Reference: xtb/src/param/atomicrad.f90 (CRC Handbook, 91st edition).
 * The QMDFF bend/torsion damping (abdamp) uses these radii, see qmdff_terms.h.
 *
 * @param z Atomic number (1..118)
 * @return radius in Angstrom; 1.0 for out-of-range elements
 */
double atomicRadAngstrom(int z);

/**
 * @brief Effective valence-electron number Z_A entering the QMDFF repulsion.
 *
 * Reference: qmdff.f90:1084-1158 (subroutine valelff). Formal valence electron
 * counts, row-wise scaled by an empirical factor (2.35 / 0.95 / 0.75 / 0.65 /
 * 0.60), with hand-fitted overrides for groups 1 and 2.
 *
 * The repulsion prefactor is zab(A,B) = Z_A * Z_B (qmdff.f90:140).
 *
 * @param z Atomic number (1..94)
 * @return Z_A; 0.0 for out-of-range elements
 */
double valenceElectrons(int z);

/**
 * @brief Element-specific hydrogen-bond strength scaling.
 *
 * Reference: qmdff.f90:103-111 (scalehb). Non-zero only for
 * N (0.8), O (0.3), F (0.1), P/S/Cl/Se/Br (2.0).
 * A zero value means "this element is not a hydrogen-bond donor/acceptor".
 */
double scaleHB(int z);

/**
 * @brief Element-specific halogen-bond strength scaling.
 *
 * Reference: qmdff.f90:112-116 (scalexb). Non-zero only for
 * Cl (0.30), Br (0.60), I (0.80), At (1.00).
 */
double scaleXB(int z);

/**
 * @brief Coulomb scaling factor for non-covalent interaction class nk.
 *
 * Reference: qmdff.f90:46/117-128 (eps1). Index nk is the value stored in
 * nci(3,·) of a QMDFF parameter file:
 *   nk = 1 : 1,2 pair (bonded)          -> 0.00
 *   nk = 2 : 1,3 pair (geminal)         -> 0.00
 *   nk = 3 : 1,4 pair                   -> 0.85
 *   nk = 4 : 1,5 pair                   -> 1.00
 *   nk = 5 : 1,6 and beyond / true NCI  -> 1.00
 *   nk = 6 : electrostatics switched off-> 0.00
 *
 * @param nk 1-based class index (1..6); returns 0.0 outside that range
 */
double eps1(int nk);

/**
 * @brief Dispersion/repulsion scaling factor for non-covalent class nk.
 *
 * Reference: qmdff.f90:47/122-127 (eps2): 0.00, 0.00, 0.50, 0.50, 1.00, 1.00.
 * @param nk 1-based class index (1..6); returns 0.0 outside that range
 */
double eps2(int nk);

} // namespace QMDFFData
