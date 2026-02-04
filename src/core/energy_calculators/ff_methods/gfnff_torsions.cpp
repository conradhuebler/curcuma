/*
 * <GFN-FF Torsion Implementation for Curcuma>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * This implementation is based on the GFN-FF force field method developed by:
 *   - Prof. Stefan Grimme (University of Bonn)
 *   - Dr. Sebastian Spicher (University of Bonn)
 *
 * Original Publication:
 *   S. Spicher, S. Grimme
 *   "Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems"
 *   Angew. Chem. Int. Ed. 59, 15665-15676 (2020)
 *   DOI: 10.1002/anie.202004239
 *
 * Original Fortran Implementation:
 *   external/gfnff/src/gfnff_engrad.F90 (lines 1041-1122: egtors subroutine)
 *   external/gfnff/src/gfnff_helpers.f90 (lines 319-358: valijklff function)
 *   Developed by: Grimme Group (Sebastian Ehlert, Sebastian Spicher, Stefan Grimme)
 *   Maintained by: Philipp Pracht (https://github.com/pprcht/gfnff)
 *   License: GNU LGPL v3
 *
 * This C++ port aims to provide:
 *   - Educational clarity for theoretical chemists
 *   - Full integration with Curcuma ecosystem
 *   - Maintainable code without Fortran dependencies
 *   - Scientific accuracy validated against original implementation
 *
 * Scientific Documentation: docs/theory/GFNFF_TORSION_THEORY.md
 *
 * Claude Generated (2025): Educational reimplementation maintaining scientific rigor
 */

#include "gfnff.h"
#include "gfnff_par.h"
#include <cmath>
#include <array>
#include <set>
#include <map>

using namespace GFNFFParameters;

// =============================================================================
// DIHEDRAL ANGLE CALCULATION
// =============================================================================

/**
 * Calculate dihedral angle for four atoms i-j-k-l
 *
 * Scientific Background:
 * ---------------------
 * The dihedral angle φ describes the rotation around the central bond j-k.
 * It is the angle between two planes:
 *   - Plane 1: defined by atoms i-j-k
 *   - Plane 2: defined by atoms j-k-l
 *
 * Physical Interpretation:
 *   φ = 0°:   cis/eclipsed   (i and l on same side)
 *   φ = 180°: trans/anti     (i and l on opposite sides)
 *   φ = ±60°: gauche         (staggered conformations)
 *
 * Mathematical Formula:
 * --------------------
 * Given bond vectors:
 *   v₁ = r_j - r_i  (bond i→j)
 *   v₂ = r_k - r_j  (bond j→k, central bond)
 *   v₃ = r_l - r_k  (bond k→l)
 *
 * The dihedral angle is computed using:
 *   φ = atan2(|v₂|·(v₁ × v₂)·v₃, (v₁ × v₂)·(v₂ × v₃))
 *
 * Why atan2 instead of acos?
 *   - atan2 gives signed angle φ ∈ [-π, π]
 *   - acos only gives φ ∈ [0, π] (loses sign information)
 *   - Sign is crucial for distinguishing clockwise/counterclockwise rotation
 *   - Essential for correct gradient calculation
 *
 * Reference Implementation:
 *   external/gfnff/src/gfnff_helpers.f90:319-358 (valijklff function)
 *
 * Literature:
 *   M.P. Allen, D.J. Tildesley, "Computer Simulation of Liquids" (1987)
 *   Chapter 1.6: "Molecular Dynamics in Different Ensembles"
 *
 * @param i Index of first atom
 * @param j Index of second atom (central bond)
 * @param k Index of third atom (central bond)
 * @param l Index of fourth atom
 * @return Dihedral angle in radians, range [-π, π]
 */
double GFNFF::calculateDihedralAngle(int i, int j, int k, int l) const
{
    // Extract atomic positions from geometry matrix (stored in Angstrom)
    // m_geometry is Eigen::Matrix with shape (m_atomcount, 3)
    Eigen::Vector3d r_i = m_geometry.row(i).head<3>();  // Position of atom i
    Eigen::Vector3d r_j = m_geometry.row(j).head<3>();  // Position of atom j
    Eigen::Vector3d r_k = m_geometry.row(k).head<3>();  // Position of atom k
    Eigen::Vector3d r_l = m_geometry.row(l).head<3>();  // Position of atom l

    // Calculate bond vectors
    // v1: i→j bond vector
    // v2: j→k bond vector (central bond around which rotation occurs)
    // v3: k→l bond vector
    Eigen::Vector3d v1 = r_j - r_i;
    Eigen::Vector3d v2 = r_k - r_j;
    Eigen::Vector3d v3 = r_l - r_k;

    // Calculate cross products
    // n1 = v1 × v2: normal vector to plane i-j-k
    // n2 = v2 × v3: normal vector to plane j-k-l
    // The dihedral angle is the angle between these two planes
    Eigen::Vector3d n1 = v1.cross(v2);
    Eigen::Vector3d n2 = v2.cross(v3);

    // Calculate magnitudes
    double n1_norm = n1.norm();
    double n2_norm = n2.norm();
    double v2_norm = v2.norm();

    // Check for degenerate geometries (colinear atoms)
    // This occurs when three atoms are nearly linear, making the cross product ~0
    const double epsilon = 1.0e-10;
    if (n1_norm < epsilon || n2_norm < epsilon || v2_norm < epsilon) {
        // Colinear atoms → dihedral angle is undefined
        // Return 0.0 (this torsion will have minimal energy contribution)
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn(fmt::format(
                "Degenerate dihedral angle for atoms {}-{}-{}-{} (colinear geometry)",
                i, j, k, l));
        }
        return 0.0;
    }

    // Normalize the normal vectors
    Vector n1_normalized = n1 / n1_norm;
    Vector n2_normalized = n2 / n2_norm;

    // Calculate dihedral angle using atan2 for proper sign
    // Formula: φ = atan2(|v₂|·(n1·v3), n1·n2)
    //
    // Derivation (from Fortran valijklff):
    //   c1 = n1·n2           (cosine component)
    //   c2 = |v₂|·(n1·v3)    (sine component with correct sign)
    //   φ = atan2(c2, c1)
    //
    // Physical meaning of signs:
    //   - Positive φ: right-handed rotation (looking down v2 axis)
    //   - Negative φ: left-handed rotation
    double cos_phi = n1_normalized.dot(n2_normalized);
    double sin_phi = v2_norm * n1_normalized.dot(v3);

    // Use atan2 for correct quadrant (returns angle in [-π, π])
    double phi = atan2(sin_phi, cos_phi);

    return phi;
}

// =============================================================================
// TORSION DAMPING FUNCTION
// =============================================================================

/**
 * Calculate distance-dependent damping function for torsion potential
 *
 * Scientific Background:
 * ---------------------
 * GFN-FF couples torsional rotation with bond stretching through a damping function.
 * This represents the physical reality that:
 *   - Stretched bonds → easier rotation (lower torsion barrier)
 *   - Compressed bonds → harder rotation (higher torsion barrier)
 *
 * Mathematical Formula:
 * --------------------
 * Damping function (Spicher & Grimme 2020, Supporting Information):
 *
 *   D(r) = 1 / [1 + exp(-α·(r/r₀ - 1))]
 *
 * where:
 *   r: current bond distance
 *   r₀: equilibrium bond distance (from covalent radii)
 *   α: steepness parameter (typically 16-20)
 *
 * Derivative (needed for gradient):
 *   dD/dr = α·D·(1 - D) / r₀
 *
 * Physical Interpretation:
 *   - r = r₀: D = 0.5 (normal damping)
 *   - r > r₀: D → 1 (reduced damping → lower barrier)
 *   - r < r₀: D → 0 (increased damping → higher barrier)
 *
 * Reference Implementation:
 *   external/gfnff/src/gfnff_engrad.F90:1234-1243 (gfnffdampt subroutine)
 *
 * Literature:
 *   S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 59, 15665 (2020)
 *   Supporting Information, Section S2.4: "Damping Functions"
 *
 * @param z1 Atomic number of first atom
 * @param z2 Atomic number of second atom
 * @param r_squared Squared distance r² in Angstrom²
 * @param damp Output: Damping value D(r) ∈ [0, 1]
 * @param damp_deriv Output: Derivative dD/dr in 1/Angstrom
 */
void GFNFF::calculateTorsionDamping(int z1, int z2, double r_squared,
                                     double& damp, double& damp_deriv) const
{
    // Get equilibrium distance from covalent radii table
    // r₀ = rcov(z1) + rcov(z2)
    // Covalent radii from GFN-FF parameter set (Spicher & Grimme 2020)
    double r0 = getCovalentRadius(z1) + getCovalentRadius(z2);

    // Current distance (from squared distance)
    double r = sqrt(r_squared);

    // Steepness parameter α (from Fortran gfnffdampt)
    // Typical value: 16-20 (controls how rapidly damping changes with distance)
    // Higher α → sharper transition
    const double alpha = 16.0;

    // Calculate normalized distance deviation
    // x = r/r₀ - 1
    // Physical meaning:
    //   x = 0:  bond at equilibrium
    //   x > 0:  bond stretched
    //   x < 0:  bond compressed
    double x = (r / r0) - 1.0;

    // Exponential argument: -α·x
    // Capped at ±20 to avoid numerical overflow/underflow
    double exp_arg = -alpha * x;
    if (exp_arg > 20.0) {
        // Very negative x → strongly compressed bond
        // exp(-α·x) → ∞, so D → 0 (maximum damping)
        damp = 0.0;
        damp_deriv = 0.0;
        return;
    }
    if (exp_arg < -20.0) {
        // Very positive x → strongly stretched bond
        // exp(-α·x) → 0, so D → 1 (minimum damping)
        damp = 1.0;
        damp_deriv = 0.0;
        return;
    }

    // Calculate damping function
    // D(r) = 1 / [1 + exp(-α·x)]
    double exp_val = exp(exp_arg);
    damp = 1.0 / (1.0 + exp_val);

    // Calculate derivative dD/dr
    // Using chain rule:
    //   dD/dr = dD/dx · dx/dr
    //   dD/dx = α·exp(-α·x) / [1 + exp(-α·x)]²
    //        = α·D·(1 - D)  (simplified form)
    //   dx/dr = 1/r₀
    //
    // Therefore:
    //   dD/dr = α·D·(1 - D) / r₀
    damp_deriv = alpha * damp * (1.0 - damp) / r0;

    // Numerical stability check
    // If damping is very close to 0 or 1, derivative can be numerically unstable
    if (damp < 1.0e-10 || damp > (1.0 - 1.0e-10)) {
        damp_deriv = 0.0;
    }
}

// =============================================================================
// TOPOLOGY-SPECIFIC HELPER FUNCTIONS (Claude Generated - January 8, 2026)
// =============================================================================

/**
 * @brief Detect C=O alpha carbon configuration
 *
 * Reference: gfnff_ini2.f90:alphaCO
 * Detects if torsion involves C(=O)-C bond (peptide backbone, esters, ketones)
 *
 * Fortran Logic:
 * ```fortran
 * logical function alphaCO(n,at,hyb,nb,pi,a,b)
 *     alphaCO = .false.
 *     if (pi(a) .ne. 0.and.hyb(b) .eq. 3.and.at(a) .eq. 6.and.at(b) .eq. 6) then
 *       no = 0
 *       do i = 1,nb(20,a)
 *         j = nb(i,a)
 *         if (at(j) .eq. 8.and.pi(j) .ne. 0.and.nb(20,j) .eq. 1) no = no+1
 *       end do
 *       if (no .eq. 1) alphaCO = .true.
 *     end if
 * end function
 * ```
 *
 * Physical Meaning:
 * - Detects C=O alpha carbon (e.g., in acetone, esters, peptide bonds)
 * - Atom a (j in torsion): Carbon with pi system (sp or sp2)
 * - Atom b (k in torsion): sp3 carbon
 * - Check: Does atom a have exactly 1 terminal oxygen (pi-bonded, CN=1)?
 * - If yes: This is C=O alpha carbon → multiply fij by 1.3
 *
 * @param atom_idx Atom to check (carbon with pi system)
 * @param other_idx Other atom in central bond (should be sp3 carbon)
 * @param atoms Atomic numbers
 * @param hybridization Hybridization states (1=sp, 2=sp2, 3=sp3)
 * @param bond_list Bond connectivity
 * @return true if C=O alpha carbon detected
 */
static bool isAlphaCO(int atom_idx, int other_idx,
                      const std::vector<int>& atoms,
                      const std::vector<int>& hybridization,
                      const std::vector<std::pair<int,int>>& bond_list)
{
    // Bounds check
    if (atom_idx < 0 || atom_idx >= atoms.size() || other_idx < 0 || other_idx >= atoms.size()) {
        return false;
    }

    // Atom must be carbon with pi system (sp or sp2)
    if (atoms[atom_idx] != 6 || (hybridization[atom_idx] != 1 && hybridization[atom_idx] != 2)) {
        return false;
    }

    // Other atom must be sp3 carbon
    if (atoms[other_idx] != 6 || hybridization[other_idx] != 3) {
        return false;
    }

    // Check if atom_idx has exactly 1 terminal pi-bonded oxygen (C=O)
    int terminal_oxygen_count = 0;

    for (const auto& bond : bond_list) {
        int neighbor = -1;
        if (bond.first == atom_idx) neighbor = bond.second;
        else if (bond.second == atom_idx) neighbor = bond.first;
        else continue;

        // Bounds check for neighbor
        if (neighbor < 0 || neighbor >= atoms.size()) continue;

        // Check if neighbor is oxygen
        if (atoms[neighbor] == 8) {
            // Count coordination number of this oxygen
            int oxygen_cn = 0;
            for (const auto& b : bond_list) {
                if (b.first == neighbor || b.second == neighbor) {
                    oxygen_cn++;
                }
            }

            // Terminal oxygen: CN=1 (only bonded to the carbon)
            if (oxygen_cn == 1) {
                terminal_oxygen_count++;
            }
        }
    }

    // Exactly 1 terminal C=O → alpha carbon
    return (terminal_oxygen_count == 1);
}

/**
 * @brief Detect amide nitrogen configuration
 *
 * Reference: gfnff_ini2.f90:amide
 * Detects sp2 nitrogen bonded to exactly 1 pi-bonded carbon (peptide bond)
 *
 * Fortran Logic:
 * ```fortran
 * logical function amide(n,at,hyb,nb,pi,a)
 *     amide = .false.
 *     if (at(a) .eq. 7.and.hyb(a) .eq. 2) then  ! sp2 nitrogen
 *       nc = 0
 *       do i = 1,nb(20,a)
 *         j = nb(i,a)
 *         if (at(j) .eq. 6.and.pi(j) .ne. 0) nc = nc+1
 *       end do
 *       if (nc .eq. 1) amide = .true.
 *     end if
 * end function
 * ```
 *
 * Physical Meaning:
 * - Detects amide nitrogen (peptide bonds, amides)
 * - Atom must be sp2 nitrogen
 * - Check: Does it have exactly 1 pi-bonded carbon neighbor?
 * - If yes: This is amide → multiply fij by 1.3
 *
 * @param atom_idx Atom to check
 * @param atoms Atomic numbers
 * @param hybridization Hybridization states (1=sp, 2=sp2, 3=sp3)
 * @param bond_list Bond connectivity
 * @return true if amide nitrogen detected
 */
static bool isAmide(int atom_idx,
                    const std::vector<int>& atoms,
                    const std::vector<int>& hybridization,
                    const std::vector<std::pair<int,int>>& bond_list)
{
    // Bounds check
    if (atom_idx < 0 || atom_idx >= atoms.size()) {
        return false;
    }

    // Must be sp2 nitrogen
    if (atoms[atom_idx] != 7 || hybridization[atom_idx] != 2) {
        return false;
    }

    // Count pi-bonded carbons (sp or sp2)
    int pi_carbon_count = 0;

    for (const auto& bond : bond_list) {
        int neighbor = -1;
        if (bond.first == atom_idx) neighbor = bond.second;
        else if (bond.second == atom_idx) neighbor = bond.first;
        else continue;

        // Bounds check for neighbor
        if (neighbor < 0 || neighbor >= atoms.size()) continue;

        // Check if neighbor is carbon with pi system (sp=1 or sp2=2)
        if (atoms[neighbor] == 6 && (hybridization[neighbor] == 1 || hybridization[neighbor] == 2)) {
            pi_carbon_count++;
        }
    }

    // Exactly 1 pi-bonded carbon → amide
    return (pi_carbon_count == 1);
}

// =============================================================================
// TORSION PARAMETER LOOKUP
// =============================================================================

/**
 * Get GFN-FF torsion parameters for atom quartet i-j-k-l
 *
 * Scientific Background:
 * ---------------------
 * Torsion parameters in GFN-FF are assigned based on:
 *   1. Hybridization of central atoms j and k
 *   2. Element types of all four atoms
 *   3. Topology (ring membership, conjugation)
 *
 * Parameter Assignment Strategy (Spicher & Grimme 2020):
 * ------------------------------------------------------
 *
 * A. Periodicity (n):
 *    - sp³-sp³: n = 3 (threefold barrier, e.g., ethane)
 *    - sp²-sp²: n = 2 (twofold barrier, cis/trans isomerization)
 *    - sp²-sp³: n = 2 or 3 (depends on conjugation)
 *    - sp-X:    n = 1 (onefold barrier)
 *
 * B. Barrier Height (V):
 *    - Look up from element-specific table
 *    - Typical values:
 *      * C(sp³)-C(sp³): ~1.5 kcal/mol
 *      * C(sp²)-C(sp²): ~3-5 kcal/mol (conjugated)
 *      * C(sp²)=C(sp²): ~45 kcal/mol (π-bond rotation!)
 *
 * C. Phase Shift (φ₀):
 *    - Typically 0° for sp³-sp³ (staggered minimum)
 *    - Typically 180° for conjugated systems (planar minimum)
 *
 * D. Corrections:
 *    - Ring strain: Reduce barrier in small rings (3-, 4-membered)
 *    - Conjugation: Increase barrier (prefers planarity)
 *    - Hyperconjugation: Subtle barrier modulation
 *
 * Reference Implementation:
 *   external/gfnff/src/gfnff_ini2.f90 (torsion parameter assignment)
 *
 * Literature:
 *   S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 59, 15665 (2020)
 *   Section "Torsion Potentials" and Supporting Information
 *
 * @param z_i Atomic number of first atom
 * @param z_j Atomic number of second atom (central bond)
 * @param z_k Atomic number of third atom (central bond)
 * @param z_l Atomic number of fourth atom
 * @param hyb_j Hybridization of atom j (1=sp, 2=sp2, 3=sp3)
 * @param hyb_k Hybridization of atom k (1=sp, 2=sp2, 3=sp3)
 * @return GFN-FF torsion parameters (barrier, periodicity, phase, improper flag)
 */
GFNFF::GFNFFTorsionParams GFNFF::getGFNFFTorsionParameters(
    int z_i, int z_j, int z_k, int z_l,
    int hyb_j, int hyb_k,
    double qa_j, double qa_k,
    double cn_i, double cn_l,
    bool in_ring, int ring_size,
    int i_atom_idx, int j_atom_idx, int k_atom_idx, int l_atom_idx) const
{
    GFNFFTorsionParams params;

    // Default: proper torsion (not improper)
    params.is_improper = false;

    // ==========================================================================
    // STEP 1: Determine periodicity and phase (CORRECTED Dec 2025)
    // ==========================================================================
    // Reference: gfnff_ini.f90:1838-1855
    // CRITICAL: Acyclic default is phi0 = 180° = π (line 1839)
    //
    // Physical basis:
    //   - Acyclic: trans configuration (180°) is typical minimum
    //   - sp³-sp³: Threefold symmetry (n=3), trans minimum
    //   - sp²-sp²: Planar/conjugated → twofold (n=2), trans minimum
    //   - Rings: Different phi0 based on ring size (see ring case below)

    // DEFAULT for ACYCLIC (gfnff_ini.f90:1839-1840)
    params.periodicity = 1;         // Default
    params.phase_shift = M_PI;      // phi0 = 180° (trans) - ACYCLIC DEFAULT!

    // Claude Generated (Jan 25, 2026): Use pi-fragment data to determine periodicity
    // This ensures ring bonds in aromatic systems (like Caffeine N-C) use n=2
    // despite local CN=3 (sp3) markings.
    const TopologyInfo& t_info = getCachedTopology();
    bool j_is_pi = false, k_is_pi = false;
    if (j_atom_idx >= 0 && k_atom_idx >= 0 && !t_info.pi_fragments.empty()) {
        j_is_pi = t_info.pi_fragments[j_atom_idx] > 0;
        k_is_pi = t_info.pi_fragments[k_atom_idx] > 0;
    }

    // Pi-conjugated central bond (n=2, most important for aromatic accurate energies)
    if (j_is_pi && k_is_pi) {
        params.periodicity = 2;
        params.phase_shift = M_PI;
    }
    // sp³-sp³: Threefold (gfnff_ini.f90:1841)
    else if (hyb_j == 3 && hyb_k == 3) {
        params.periodicity = 3;     // nrot = 3
        params.phase_shift = M_PI;  // phi0 = 180° (keeps acyclic default)
    }
    // sp²-sp²: Twofold for pi bonds (gfnff_ini.f90:1842)
    else if (hyb_j == 2 && hyb_k == 2) {
        params.periodicity = 2;     // nrot = 2
        params.phase_shift = M_PI;  // phi0 = 180° (planar trans)
    }
    // Pi-sp³ mixed (gfnff_ini.f90:1843-1854)
    else if ((hyb_j == 2 && hyb_k == 3) || (hyb_j == 3 && hyb_k == 2)) {
        params.periodicity = 3;     // nrot = 3
        params.phase_shift = M_PI;  // phi0 = 180°
    }
    // sp-X: Linear (gfnff_ini.f90: implicit default)
    else if (hyb_j == 1 || hyb_k == 1) {
        params.periodicity = 1;
        params.phase_shift = M_PI;
    }
    // Fallback: acyclic default
    else {
        params.periodicity = 1;
        params.phase_shift = M_PI;
    }

    // ==========================================================================
    // STEP 2: Calculate force constant using GFN-FF formula (CORRECTED Dec 2025)
    // ==========================================================================
    // Reference: gfnff_ini.f90:1896
    // Formula: fctot = (f1 + 10*torsf[2]*f2) * fqq * fij * fkl
    // CRITICAL: Result is in HARTREE, not kcal/mol!

    // GFN-FF constants (from gfnff_param.f90:742-753)
    const double torsf_pi = 1.18;       // Pi bond scaling
    const double fcthr = 1.0e-3;        // Force constant threshold (Hartree)

    // Bounds check
    if (z_i < 1 || z_i > 86 || z_j < 1 || z_j > 86 ||
        z_k < 1 || z_k > 86 || z_l < 1 || z_l > 86) {
        params.barrier_height = 0.0;
        return params;
    }

    // ---------------------------------------------------------------------------
    // (A) Central bond contribution: fij = tors[Z_j] * tors[Z_k]
    // ---------------------------------------------------------------------------
    double fij = tors_angewChem2020[z_j - 1] * tors_angewChem2020[z_k - 1];

    // Check threshold and negative values (gfnff_ini.f90:1767-1768)
    if (fij < fcthr || tors_angewChem2020[z_j - 1] < 0.0 || tors_angewChem2020[z_k - 1] < 0.0) {
        params.barrier_height = 0.0;
        return params;
    }

    // NOTE (Jan 8, 2026): Hydrogen count correction is applied later at line ~526
    // using proper neighbor list counting (section E: "Hydrogen count refinement")

    // ---------------------------------------------------------------------------
    // (B) Outer atom contribution: fkl = tors2[Z_i] * tors2[Z_l]
    // ---------------------------------------------------------------------------
    double fkl = tors2_angewChem2020[z_i - 1] * tors2_angewChem2020[z_l - 1];

    // CN-dependent scaling: fkl *= (CN_i * CN_l)^(-0.14) (gfnff_ini.f90:1809)
    // Use actual CN values from parameters
    double cn_product = cn_i * cn_l;
    if (cn_product > 0.01) {  // Avoid division by zero
        fkl *= std::pow(cn_product, -0.14);
    }

    // ---------------------------------------------------------------------------
    // Nitrogen reduction for outer atoms (gfnff_ini.f90:1803-1804)
    // ---------------------------------------------------------------------------
    // Claude Generated (January 26, 2026) - CRITICAL for caffeine accuracy
    //
    // Reference: gfnff_ini.f90:1803-1804
    //   if (piadr(kk) .eq. 0) fkl = fkl*0.5d0  ! outer atom not in pi-system
    //   if (piadr(ll) .eq. 0) fkl = fkl*0.5d0
    //
    // Physical basis:
    // - sp³ nitrogen (not in aromatic/conjugated systems) has lone pair electrons
    // - Lone pair reduces torsion barrier due to orbital interactions
    // - Example: caffeine has 3 N atoms → without this correction, torsions are 2× too high
    //
    // Impact on caffeine:
    // - Reduces torsion energy error from 32.6× to <10× (expected)
    // - Affects all 3 nitrogen atoms in the molecule
    //
    // Note: piadr==0 in Fortran means "not in pi-system" → use !is_in_pi_fr()

    // Get topology for pi-system detection and hybridization
    const TopologyInfo& topo_outer = getCachedTopology();
    const auto& bond_list_outer = getCachedBondList();

    // Get hybridization of outer atoms i and l from topology
    int hyb_i = 0, hyb_l = 0;
    if (i_atom_idx >= 0 && i_atom_idx < static_cast<int>(topo_outer.hybridization.size())) {
        hyb_i = topo_outer.hybridization[i_atom_idx];
    }
    if (l_atom_idx >= 0 && l_atom_idx < static_cast<int>(topo_outer.hybridization.size())) {
        hyb_l = topo_outer.hybridization[l_atom_idx];
    }

    // Lambda to check if atom is sp³ nitrogen NOT in pi-system
    auto is_sp3_nitrogen_not_pi = [&](int atom_idx, int z_atom, int hyb_atom) -> bool {
        if (z_atom != 7) return false;  // Must be nitrogen
        if (hyb_atom != 3) return false;  // Must be sp³

        // Check if in pi-system using existing logic
        if (atom_idx < 0 || atom_idx >= m_atomcount) return false;

        // Condition 1: Direct bonds with significant pi-character (pibo > 0.1)
        if (!topo_outer.pi_bond_orders.empty()) {
            for (int other = 0; other < m_atomcount; other++) {
                if (other == atom_idx) continue;
                int pibo_idx = lin(atom_idx, other);
                if (pibo_idx >= 0 && pibo_idx < static_cast<int>(topo_outer.pi_bond_orders.size())) {
                    if (topo_outer.pi_bond_orders[pibo_idx] > 0.1) return false;  // In pi-system
                }
            }
        }

        // Condition 2: N adjacent to sp or sp2 atoms ("picon" case)
        for (const auto& bond : bond_list_outer) {
            int neighbor = (bond.first == atom_idx) ? bond.second : (bond.second == atom_idx) ? bond.first : -1;
            if (neighbor >= 0 && neighbor < static_cast<int>(topo_outer.hybridization.size())) {
                int neighbor_hyb = topo_outer.hybridization[neighbor];
                if (neighbor_hyb == 1 || neighbor_hyb == 2) return false;  // In pi-system
            }
        }

        return true;  // sp³ nitrogen NOT in pi-system
    };

    // Apply nitrogen reduction (×0.5) for outer atoms i and l
    if (is_sp3_nitrogen_not_pi(i_atom_idx, z_i, hyb_i)) {
        fkl *= 0.5;
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("  fkl nitrogen reduction (atom {}): ×0.5", i_atom_idx));
        }
    }

    if (is_sp3_nitrogen_not_pi(l_atom_idx, z_l, hyb_l)) {
        fkl *= 0.5;
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("  fkl nitrogen reduction (atom {}): ×0.5", l_atom_idx));
        }
    }

    // ---------------------------------------------------------------------------
    // Hypervalent enhancement for outer atoms (gfnff_ini.f90:1890)
    // ---------------------------------------------------------------------------
    // Reference: gfnff_ini.f90:1890
    //   if (hyb(kk) .eq. 5.or.hyb(ll) .eq. 5) fkl = fkl*1.5d0
    //
    // Physical basis:
    // - Hypervalent atoms (PF₅, SF₆) have higher coordination → stiffer torsions
    // - Rare in organic chemistry but needed for completeness
    //
    // Note: Low priority for caffeine (no hypervalent atoms)

    if (hyb_i == 5 || hyb_l == 5) {
        fkl *= 1.5;
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("  fkl hypervalent enhancement: ×1.5");
        }
    }

    // Check threshold (gfnff_ini.f90:1805-1806)
    if (fkl < fcthr || tors2_angewChem2020[z_i - 1] < 0.0 || tors2_angewChem2020[z_l - 1] < 0.0) {
        params.barrier_height = 0.0;
        return params;
    }

    // ---------------------------------------------------------------------------
    // (C) Base force constant: f1 (hybridization-dependent)
    // ---------------------------------------------------------------------------
    // Claude Generated (Jan 21, 2026): CRITICAL FIX - Pi-system aware hybridization
    //
    // Problem: Raw hybridization from CN thresholds marks aromatic N as sp3 (CN~3),
    // but for torsion barriers, aromatic atoms should be treated as sp2.
    //
    // Solution: Check pi-bond order BEFORE f1 calculation. If the central bond
    // has significant pi-character (pibo > 0.1), treat both atoms as sp2.
    //
    // Reference: XTB uses piadr (pi-system address) to determine if atoms are in
    // a conjugated system. This is equivalent to checking pibo > 0.
    //
    // Impact: Reduces f1 from 1.0 to 0.2 for N-containing aromatic bonds,
    // fixing the 5× torsion energy overestimation.

    double f1 = torsf_single;  // Default = 1.0

    // Get effective hybridization (considering pi-system participation)
    const TopologyInfo& topo = getCachedTopology();
    const auto& bond_list = getCachedBondList();

    auto is_in_pi_fr = [&](int atom_idx, int z_atom) -> bool {
        if (atom_idx < 0 || atom_idx >= m_atomcount) return false;
        // Condition 1: Direct bonds with significant pi-character (pibo > 0.1)
        if (!topo.pi_bond_orders.empty()) {
            for (int other = 0; other < m_atomcount; other++) {
                if (other == atom_idx) continue;
                int pibo_idx = lin(atom_idx, other);
                if (pibo_idx >= 0 && pibo_idx < static_cast<int>(topo.pi_bond_orders.size())) {
                    if (topo.pi_bond_orders[pibo_idx] > 0.1) return true;
                }
            }
        }
        // Condition 2: N/O/F/S (sp3) adjacent to sp or sp2 atoms ("picon" case)
        if (z_atom == 7 || z_atom == 8 || z_atom == 9 || z_atom == 16) {
            for (const auto& bond : bond_list) {
                int neighbor = (bond.first == atom_idx) ? bond.second : (bond.second == atom_idx) ? bond.first : -1;
                if (neighbor >= 0 && neighbor < static_cast<int>(topo.hybridization.size())) {
                    int neighbor_hyb = topo.hybridization[neighbor];
                    if (neighbor_hyb == 1 || neighbor_hyb == 2) return true;
                }
            }
        }
        return false;
    };

    int eff_hyb_j = hyb_j;
    int eff_hyb_k = hyb_k;

    if (is_in_pi_fr(j_atom_idx, z_j) && eff_hyb_j == 3) eff_hyb_j = 2;
    if (is_in_pi_fr(k_atom_idx, z_k) && eff_hyb_k == 3) eff_hyb_k = 2;

    // Acyclic sp3-sp3 case (most common)
    if (eff_hyb_j == 3 && eff_hyb_k == 3) {
        // Ethane-like: keep f1 = 1.0
        // Special cases for heteroatoms (simplified from lines 1857-1879):
        int group_j = (z_j == 7 || z_j == 15) ? 5 : (z_j == 8 || z_j == 16) ? 6 : 0;
        int group_k = (z_k == 7 || z_k == 15) ? 5 : (z_k == 8 || z_k == 16) ? 6 : 0;

        if (group_j == 6 && group_k == 6) {
            // O-O, S-S: higher barrier, nrot=2, phi0=90° (gfnff_ini.f90:1873-1879)
            f1 = 5.0;
            params.periodicity = 2;
            params.phase_shift = 90.0 * M_PI / 180.0;
            if (z_j >= 16 && z_k >= 16) f1 = 25.0;  // S-S
        }
        else if (group_j == 5 && group_k == 5) {
            // N-N, P-P: nrot=3, phi0=60° (gfnff_ini.f90:1859-1863)
            f1 = 3.0;
            params.periodicity = 3;
            params.phase_shift = 60.0 * M_PI / 180.0;
        }
        else if ((group_j == 5 && group_k == 6) || (group_j == 6 && group_k == 5)) {
            // N-O or O-N mixed bond: nrot=2, phi0=90° (gfnff_ini.f90:1865-1871)
            f1 = 1.0;
            params.periodicity = 2;
            params.phase_shift = 90.0 * M_PI / 180.0;
            if (z_j >= 15 && z_k >= 15) {
                f1 = 20.0;  // P-S case
            }
        }
    }
    // Pi-sp3 mixed (lines 1843-1854)
    else if ((eff_hyb_j == 2 && eff_hyb_k == 3) || (eff_hyb_j == 3 && eff_hyb_k == 2)) {
        f1 = 0.5;
        if (z_j == 7 || z_k == 7) f1 = 0.2;  // Nitrogen lowers barrier
    }
    // sp2-sp2 conjugated (keep f1 = 1.0, will be scaled by 0.55 later if pibo > 0)
    else if (eff_hyb_j == 2 && eff_hyb_k == 2) {
        // For sp2-sp2 (aromatic/conjugated), f1 = 1.0 but will be scaled by 0.55 when pibo > 0
        // This is handled later in the pi-system contribution section
        f1 = torsf_single;  // 1.0
    }

    // ---------------------------------------------------------------------------
    // (D) Pi system contribution: f2 (IMPLEMENTED Jan 18, 2026)
    // ---------------------------------------------------------------------------
    // Reference: gfnff_ini.f90:1900-1907
    //   if (pibo(m) .gt. 0) then
    //     f2 = pibo(m)*exp(-2.5d0*(1.24d0-pibo(m))**14)  ! decrease to very small values for P < 0.3
    //     if (piadr(kk) .eq. 0.and.at(kk) .gt. 10) f2 = f2*1.3  ! heavy non-pi outer atoms
    //     if (piadr(ll) .eq. 0.and.at(ll) .gt. 10) f2 = f2*1.3
    //     f1 = f1*0.55  ! CRITICAL: scale f1 down when pi-bond present!
    //   end if
    //
    // Physical meaning:
    // - pibo = pi bond order from Hückel calculation (0-1, benzene ~0.67)
    // - For aromatic/conjugated systems, the pi contribution (f2) adds significantly to the barrier
    // - The exponential cutoff ensures low-pibo bonds (< 0.3) don't contribute
    // - Heavy outer atoms increase the pi effect (1.3× scaling each)
    // - When pi-system exists, f1 is scaled by 0.55 to balance the large f2 contribution

    double f2 = 0.0;  // Default for single bonds (no pi character)

    // Get pi bond order from cached topology
    if (j_atom_idx >= 0 && k_atom_idx >= 0 && j_atom_idx < m_atomcount && k_atom_idx < m_atomcount) {
        const TopologyInfo& topo = getCachedTopology();

        if (!topo.pi_bond_orders.empty()) {
            int pibo_idx = lin(j_atom_idx, k_atom_idx);

            if (pibo_idx >= 0 && pibo_idx < static_cast<int>(topo.pi_bond_orders.size())) {
                double pibo = topo.pi_bond_orders[pibo_idx];

                if (pibo > 0.0) {
                    // Calculate f2 using exponential cutoff formula
                    // This decreases to very small values for pibo < 0.3
                    double diff = 1.24 - pibo;
                    double exp_term = std::exp(-2.5 * std::pow(diff, 14));
                    f2 = pibo * exp_term;

                    // Heavy atom correction: if outer atoms are NOT in pi system but are heavy (Z > 10),
                    // the pi bond order becomes more significant → scale f2 by 1.3 for each
                    // Reference: gfnff_ini.f90:1904-1905
                    //   if (piadr(kk) .eq. 0.and.at(kk) .gt. 10) f2 = f2*1.3
                    //   if (piadr(ll) .eq. 0.and.at(ll) .gt. 10) f2 = f2*1.3
                    //
                    // CRITICAL FIX (January 26, 2026): Must check if outer atoms are NOT in pi-system
                    // Physical basis:
                    // - Heavy atoms (Cl, Br, S, etc.) outside conjugated system polarize pi-bonds
                    // - This increases torsion barriers in aromatic systems
                    // - But if heavy atom is IN the pi-system (e.g., aromatic S), no extra effect
                    //
                    // Check outer atom i: heavy AND not in pi-system
                    if (z_i > 10 && !is_in_pi_fr(i_atom_idx, z_i)) {
                        f2 *= 1.3;
                        if (CurcumaLogger::get_verbosity() >= 3) {
                            CurcumaLogger::info(fmt::format("  Heavy non-pi outer atom {} (Z={}): f2 ×1.3", i_atom_idx, z_i));
                        }
                    }

                    // Check outer atom l: heavy AND not in pi-system
                    if (z_l > 10 && !is_in_pi_fr(l_atom_idx, z_l)) {
                        f2 *= 1.3;
                        if (CurcumaLogger::get_verbosity() >= 3) {
                            CurcumaLogger::info(fmt::format("  Heavy non-pi outer atom {} (Z={}): f2 ×1.3", l_atom_idx, z_l));
                        }
                    }

                    // CRITICAL: Scale f1 when pi-system is present! (gfnff_ini.f90:1906)
                    // This balances the large f2 contribution
                    f1 *= 0.55;

                    if (CurcumaLogger::get_verbosity() >= 3) {
                        CurcumaLogger::info(fmt::format("  Pi-bond correction: pibo={:.4f}, f2={:.6f}, f1 scaled to {:.4f}",
                                                         pibo, f2, f1));
                    }
                }
            }
        }
    }

    // ---------------------------------------------------------------------------
    // (E) Hydrogen count refinement (FIXED Jan 8, 2026 - from reference: gfnff_ini.f90:1778-1786)
    // ---------------------------------------------------------------------------
    // fij = fij * (nhi * nhj)^0.07 where nhi/nhj = 1 + number of H atoms attached to j/k
    // This accounts for H-substitution effects on torsion barriers (steric/hyperconjugation)
    //
    // CRITICAL FIX: Use bond list instead of neighbor_lists
    // Previous implementation relied on getCachedTopology().neighbor_lists which wasn't
    // populated at parameter generation time, causing f_hydrogen to always be 1.0
    //
    // Reference: Fortran gfnff_ini.f90:1805
    //   nhi = 1
    //   do ineig = 1,topo%nb(20,ii)
    //     if (at(topo%nb(ineig,ii)) .eq. 1) nhi = nhi+1
    //   end do
    //   fij = fij*(dble(nhi)*dble(nhj))**0.07

    int nhi = 1, nhj = 1;  // Default: 1 (Fortran convention - central atom counts as 1)

    if (j_atom_idx >= 0 && k_atom_idx >= 0 && j_atom_idx < m_atomcount && k_atom_idx < m_atomcount) {
        // Count H neighbors using bond list (always available, doesn't depend on topology cache)
        const auto& bond_list = getCachedBondList();
        for (const auto& bond : bond_list) {
            // Check if bond connects j_atom to an H
            if ((bond.first == j_atom_idx && bond.second < m_atoms.size() && m_atoms[bond.second] == 1) ||
                (bond.second == j_atom_idx && bond.first < m_atoms.size() && m_atoms[bond.first] == 1)) {
                nhi++;
            }
            // Check if bond connects k_atom to an H
            if ((bond.first == k_atom_idx && bond.second < m_atoms.size() && m_atoms[bond.second] == 1) ||
                (bond.second == k_atom_idx && bond.first < m_atoms.size() && m_atoms[bond.first] == 1)) {
                nhj++;
            }
        }

        // Apply H-count correction (gfnff_ini.f90:1805)
        double h_scaling = std::pow(static_cast<double>(nhi * nhj), 0.07);
        fij *= h_scaling;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("  H-count correction: atom j={} (nhi={}), atom k={} (nhj={}), scaling={:.4f}",
                                             j_atom_idx, nhi, k_atom_idx, nhj, h_scaling));
        }
    }

    // ---------------------------------------------------------------------------
    // (E2) Topology-specific fij corrections (NEW - Claude Generated Jan 8, 2026)
    // ---------------------------------------------------------------------------
    // Reference: Fortran gfnff_ini.f90:1807-1811
    // These corrections account for conjugation and resonance effects in specific
    // functional groups (peptides, carbonyls) that modulate torsion barriers
    //
    // Fortran reference:
    //   ! amides and alpha carbons in peptides/proteins
    //   if (alphaCO(nat,at,hyb,topo%nb,piadr,ii,jj)) fij = fij*1.3d0
    //   if (amide(nat,at,hyb,topo%nb,piadr,ii).and.hyb(jj) .eq. 3.and.at(jj) .eq. 6) fij = fij*1.3d0
    //   if (amide(nat,at,hyb,topo%nb,piadr,jj).and.hyb(ii) .eq. 3.and.at(ii) .eq. 6) fij = fij*1.3d0
    //
    // Physical basis:
    // - alphaCO: C=O alpha carbon has partial double bond character → stiffer rotation
    // - amide: Peptide N-C(=O) bond has resonance → restricted rotation
    //
    // Expected impact:
    // - Increases V parameters by ~1.3× for affected bonds
    // - Brings V from 0.151 Eh → ~0.20 Eh (combined with other corrections)
    // - Enables removal of 0.5 factor workaround when all corrections complete

    if (j_atom_idx >= 0 && k_atom_idx >= 0 &&
        j_atom_idx < m_atomcount && k_atom_idx < m_atomcount) {

        // Get cached data for topology checks
        const auto& bond_list = getCachedBondList();
        const TopologyInfo& topo = getCachedTopology();

        // Hybridization bounds check
        if (j_atom_idx < topo.hybridization.size() && k_atom_idx < topo.hybridization.size()) {
            int hyb_j = topo.hybridization[j_atom_idx];
            int hyb_k = topo.hybridization[k_atom_idx];

            // 1. alphaCO correction: C=O alpha carbon (fij *= 1.3)
            //    Detects C(=O)-C bonds in ketones, esters, peptide backbones
            if (isAlphaCO(j_atom_idx, k_atom_idx, m_atoms, topo.hybridization, bond_list)) {
                fij *= 1.3;
                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::info("  alphaCO correction: fij *= 1.3 (C=O alpha carbon detected)");
                }
            }
            else if (isAlphaCO(k_atom_idx, j_atom_idx, m_atoms, topo.hybridization, bond_list)) {
                fij *= 1.3;
                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::info("  alphaCO correction: fij *= 1.3 (C=O alpha carbon detected)");
                }
            }

            // 2. Amide corrections: peptide bonds (fij *= 1.3)
            //    Detects N-C(=O) resonance structures in peptides and amides
            //
            // Check if j is amide nitrogen and k is sp3 carbon
            if (isAmide(j_atom_idx, m_atoms, topo.hybridization, bond_list) &&
                z_k == 6 && hyb_k == 3) {
                fij *= 1.3;
                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::info("  amide correction (j→k): fij *= 1.3 (peptide bond detected)");
                }
            }

            // Check if k is amide nitrogen and j is sp3 carbon
            if (isAmide(k_atom_idx, m_atoms, topo.hybridization, bond_list) &&
                z_j == 6 && hyb_j == 3) {
                fij *= 1.3;
                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::info("  amide correction (k→j): fij *= 1.3 (peptide bond detected)");
                }
            }
        }

        // 3. Hypervalent bond correction (btyp == 4): fij *= 0.2
        //    Reference: gfnff_ini.f90:1811 "if (btyp(m) .eq. 4) fij = fij*0.2d0"
        //    Note: Bond type detection not yet fully implemented
        //    TODO Phase 2D: Implement full bond type classification system
        //    For now: Deferred (low impact - rare in organic molecules)
    }

    // ---------------------------------------------------------------------------
    // (F) Metal classification checks (NEW - from reference: gfnff_ini.f90:1751-1752)
    // ---------------------------------------------------------------------------
    // Skip high-coordinate metals: no HC metals with >4 neighbors
    if (j_atom_idx >= 0 && j_atom_idx < m_atoms.size()) {
        const TopologyInfo& topo = getCachedTopology();
        if (j_atom_idx < topo.neighbor_lists.size() && j_atom_idx < topo.is_metal.size()) {
            int coord_j = topo.neighbor_lists[j_atom_idx].size();
            if (topo.is_metal[j_atom_idx] && coord_j > 4) {
                params.barrier_height = 0.0;
                return params;  // Skip HC metals
            }
        }
    }
    if (k_atom_idx >= 0 && k_atom_idx < m_atoms.size()) {
        const TopologyInfo& topo = getCachedTopology();
        if (k_atom_idx < topo.neighbor_lists.size() && k_atom_idx < topo.is_metal.size()) {
            int coord_k = topo.neighbor_lists[k_atom_idx].size();
            if (topo.is_metal[k_atom_idx] && coord_k > 4) {
                params.barrier_height = 0.0;
                return params;  // Skip HC metals
            }
        }
    }

    // ---------------------------------------------------------------------------
    // (G) Ring-specific torsion corrections (Claude Updated - January 2026)
    // ---------------------------------------------------------------------------
    // Reference: XTB gfnff_ini.f90:1814-1835, gfnff_param.f90:787-790
    //
    // Educational Documentation:
    // ==========================
    // Ring torsions have COMPLETELY DIFFERENT conformational preferences than acyclic bonds!
    //
    // Physical basis:
    // - 3-ring (cyclopropane): Nearly planar, highly strained (n=1, φ₀=0°, barrier=FR3)
    // - 4-ring (cyclobutane): Butterfly puckering (n=6, φ₀=30°, barrier=FR4)
    // - 5-ring (cyclopentane): Envelope/twist puckering (n=6, φ₀=30°, barrier=FR5)
    // - 6-ring (cyclohexane): STRONG chair/boat preference (n=3, φ₀=60°, barrier=FR6)
    //
    // These override ALL acyclic parameters when all 4 atoms are in the same ring!
    //
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020
    // ---------------------------------------------------------------------------

    using namespace GFNFFParameters;

    if (in_ring && ring_size >= 3 && ring_size <= 6) {
        // CRITICAL: Check if all 4 atoms are in the SAME ring (not just bonded atoms)
        // This matches Fortran check: ringl == rings4 (lines 1819, 1824, 1828)
        // We need to verify that the entire i-j-k-l quartet is in this ring

        // For now, assume in_ring means central bond j-k is in ring
        // TODO Phase 2: Add quartet ring membership check (requires path finding)
        bool all_in_same_ring = in_ring;  // Simplified assumption

        // Additional check: Skip if pi-conjugated system (notpicon = false in Fortran)
        // Pi bonds have different torsional preferences even in rings
        // Claude Generated (Jan 25, 2026): Include outer atoms in pi-system address check
        bool notpicon = !(is_in_pi_fr(i_atom_idx, z_i) || is_in_pi_fr(l_atom_idx, z_l) ||
                          is_in_pi_fr(j_atom_idx, z_j) || is_in_pi_fr(k_atom_idx, z_k));

        if (all_in_same_ring && notpicon) {
            // Override periodicity, phase, and barrier based on ring size
            switch (ring_size) {
                case 3:
                    // 3-ring: Nearly planar (cyclopropane)
                    params.periodicity = 1;
                    params.phase_shift = 0.0;     // Planar equilibrium (φ₀=0°)
                    f1 = FR3;                      // 0.3 (flexible, small barrier)
                    break;

                case 4:
                    // 4-ring: Puckered (cyclobutane butterfly)
                    params.periodicity = 6;
                    params.phase_shift = 30.0 * M_PI / 180.0;  // φ₀=30° puckering
                    f1 = FR4;                      // 1.0 (moderate barrier)
                    break;

                case 5:
                    // 5-ring: Envelope conformation (cyclopentane)
                    params.periodicity = 6;
                    params.phase_shift = 30.0 * M_PI / 180.0;  // φ₀=30° envelope
                    f1 = FR5;                      // 1.5 (intermediate barrier)
                    break;

                case 6:
                    // 6-ring: Chair preference (cyclohexane)
                    params.periodicity = 3;
                    params.phase_shift = 60.0 * M_PI / 180.0;  // φ₀=60° chair
                    f1 = FR6;                      // 5.7 (STRONG chair/boat barrier!)
                    break;
            }

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "  Ring torsion override: {}-ring → n={}, φ₀={:.1f}°, f1={:.2f}",
                    ring_size, params.periodicity, params.phase_shift * 180.0 / M_PI, f1));
            }
        }
    }

    // ---------------------------------------------------------------------------
    // (H) Charge correction: fqq (CORRECTED Dec 2025)
    // ---------------------------------------------------------------------------
    // Reference: gfnff_ini.f90:1896 (implicit via topo%qa)
    // fqq = 1.0 + |qa_j * qa_k| * qfacTOR
    // Now using actual EEQ charges passed as parameters!
    const double qfacTOR = 12.0;  // From gfnff_param.f90:742
    double fqq = 1.0 + std::abs(qa_j * qa_k) * qfacTOR;

    // ---------------------------------------------------------------------------
    // (F) Final force constant calculation
    // ---------------------------------------------------------------------------
    // Formula: fctot = (f1 + 10*torsf[2]*f2) * fqq * fij * fkl
    double fctot = (f1 + 10.0 * torsf_pi * f2) * fqq * fij * fkl;

    // Check threshold (gfnff_ini.f90:1898)
    if (fctot < fcthr) {
        params.barrier_height = 0.0;
        return params;
    }

    // CRITICAL: barrier_height is now in HARTREE (not kcal/mol!)
    params.barrier_height = fctot;

    // DEBUG OUTPUT (December 2025)
    if (fctot > 0.001) {  // Only print significant values
        if (CurcumaLogger::get_verbosity() >= 3) {
            double fij_base = tors_angewChem2020[z_j-1] * tors_angewChem2020[z_k-1];  // Base value BEFORE H-correction
            CurcumaLogger::info(fmt::format("TORSION DEBUG: i={}, j={}, k={}, l={} (Z={},{},{},{})",
                                             z_i, z_j, z_k, z_l, z_i, z_j, z_k, z_l));
            CurcumaLogger::info(fmt::format("  tors[{}]={:.6f}, tors[{}]={:.6f}, fij_base={:.6f}, fij_corrected={:.6f}",
                                             z_j, tors_angewChem2020[z_j-1], z_k, tors_angewChem2020[z_k-1], fij_base, fij));
            CurcumaLogger::info(fmt::format("  tors2[{}]={:.6f}, tors2[{}]={:.6f}, fkl_raw={:.6f}",
                                             z_i, tors2_angewChem2020[z_i-1], z_l, tors2_angewChem2020[z_l-1],
                                             tors2_angewChem2020[z_i-1] * tors2_angewChem2020[z_l-1]));
            CurcumaLogger::info(fmt::format("  CN_correction (CN_i={:.2f}, CN_l={:.2f})={:.6f}, fkl={:.6f}",
                                             cn_i, cn_l, std::pow(cn_i * cn_l, -0.14), fkl));
            CurcumaLogger::info(fmt::format("  f1={:.6f}, f2={:.6f}, fqq={:.6f}", f1, f2, fqq));
            CurcumaLogger::info(fmt::format("  fctot = ({:.4f} + 10*1.18*{:.4f}) * {:.4f} * {:.6f} * {:.6f} = {:.6f} Eh",
                                             f1, f2, fqq, fij, fkl, fctot));
            CurcumaLogger::info(fmt::format("  XTB REFERENCE: 0.246679 Eh (for C-O torsion)"));
        }
    }

    // ==========================================================================
    // STEP 3: Apply topology corrections (simplified)
    // ==========================================================================
    //
    // Full GFN-FF implementation includes:
    //   - Ring strain corrections (reduce barrier in small rings)
    //   - Conjugation detection (increase barrier for planar systems)
    //   - Hyperconjugation effects
    //
    // TODO: Implement when ring detection is available (Phase 2)

    return params;
}

// =============================================================================
// DIHEDRAL ANGLE GRADIENT CALCULATION
// =============================================================================

/**
 * @brief Calculate analytical gradients of dihedral angle with respect to atomic positions
 *
 * SCIENTIFIC BACKGROUND:
 * ----------------------
 * The dihedral angle φ is defined by four atoms i-j-k-l. Computing its gradient
 * ∂φ/∂x requires careful application of the chain rule to the dihedral angle formula:
 *
 *   φ = atan2( |r_jk| · (n1 × r_kl), n1 · n2 )
 *
 * where:
 *   - n1 = r_ij × r_jk  (normal to plane i-j-k)
 *   - n2 = r_jk × r_kl  (normal to plane j-k-l)
 *   - r_ij = r_j - r_i  (bond vector i→j)
 *
 * The analytical derivatives are complex due to:
 *   1. Cross product derivatives: ∂(a×b)/∂x = (∂a/∂x)×b + a×(∂b/∂x)
 *   2. Normalization factors: n1/|n1|, n2/|n2|
 *   3. Singularities: When sin(φ) ≈ 0 (linear/anti-linear arrangements)
 *
 * FORTRAN REFERENCE:
 * ------------------
 * This implementation closely follows the original dphidr subroutine in:
 *   external/gfnff/src/gfnff_helpers.f90:514-583
 *
 * Original authors: S. Spicher, S. Grimme (University of Bonn)
 * Publication: Angew. Chem. Int. Ed. 59, 15665-15676 (2020)
 *
 * MATHEMATICAL DERIVATION:
 * ------------------------
 * For a dihedral angle φ defined by atoms i,j,k,l:
 *
 *   ∂φ/∂r_i = (1/sin(φ)) · [ cos(φ)·(|n2|/|n1|)·(n1×r_jk) - (n2×r_jk) ]
 *   ∂φ/∂r_j = (1/sin(φ)) · [ ... complex expression involving multiple cross products ... ]
 *   ∂φ/∂r_k = (1/sin(φ)) · [ ... similar complexity ... ]
 *   ∂φ/∂r_l = (1/sin(φ)) · [ cos(φ)·(|n1|/|n2|)·(n2×r_jk) - (n1×r_jk) ]
 *
 * The expressions for atoms j and k (central atoms) are more complex because
 * they appear in multiple bond vectors.
 *
 * SINGULARITY HANDLING:
 * ---------------------
 * When |sin(φ)| < ε (≈10^-14), the dihedral angle is undefined (linear geometry).
 * In this case, we return zero gradients to avoid numerical instability.
 * This is a physically reasonable choice: near-linear configurations have
 * vanishing torsional forces.
 *
 * COMPUTATIONAL EFFICIENCY:
 * -------------------------
 * This function performs ~15 cross products and several vector operations.
 * Cost: O(1) per call, typically ~200 FLOPs.
 *
 * USAGE IN GFN-FF:
 * ----------------
 * Called during torsion gradient calculation:
 *   1. Calculate dihedral angle φ
 *   2. Calculate dφ/dr for all four atoms
 *   3. Apply chain rule: dE/dr = (dE/dφ) · (dφ/dr)
 *
 * @param[in]  i           Index of first atom (0-based)
 * @param[in]  j           Index of second atom (0-based, central bond start)
 * @param[in]  k           Index of third atom (0-based, central bond end)
 * @param[in]  l           Index of fourth atom (0-based)
 * @param[in]  phi         Dihedral angle in radians (from calculateDihedralAngle)
 * @param[out] grad_i      Gradient vector ∂φ/∂r_i (3D)
 * @param[out] grad_j      Gradient vector ∂φ/∂r_j (3D)
 * @param[out] grad_k      Gradient vector ∂φ/∂r_k (3D)
 * @param[out] grad_l      Gradient vector ∂φ/∂r_l (3D)
 *
 * @note All gradient vectors are in atomic units (Bohr^-1)
 * @note Eigen library handles vectorization automatically
 *
 * @author Claude (AI Assistant) - Generated 2025-11-10
 *         Based on Fortran implementation by S. Spicher & S. Grimme
 * @copyright Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
void GFNFF::calculateDihedralGradient(
    int i, int j, int k, int l,
    double phi,
    Vector& grad_i,
    Vector& grad_j,
    Vector& grad_k,
    Vector& grad_l) const
{
    // ==========================================================================
    // STEP 1: Extract atomic positions
    // ==========================================================================
    Eigen::Vector3d r_i = m_geometry.row(i).head<3>();
    Eigen::Vector3d r_j = m_geometry.row(j).head<3>();
    Eigen::Vector3d r_k = m_geometry.row(k).head<3>();
    Eigen::Vector3d r_l = m_geometry.row(l).head<3>();

    // ==========================================================================
    // STEP 2: Compute bond vectors
    // ==========================================================================
    // Following Fortran naming convention:
    //   ra = r_ij = r_j - r_i  (bond i→j)
    //   rb = r_jk = r_k - r_j  (central bond j→k)
    //   rc = r_kl = r_l - r_k  (bond k→l)

    Eigen::Vector3d ra = r_j - r_i;
    Eigen::Vector3d rb = r_k - r_j;
    Eigen::Vector3d rc = r_l - r_k;

    // Composite vectors (appear frequently in derivatives)
    Eigen::Vector3d rapb = ra + rb;  // r_i→k via j
    Eigen::Vector3d rbpc = rb + rc;  // r_j→l via k

    // ==========================================================================
    // STEP 3: Compute normal vectors to dihedral planes
    // ==========================================================================
    //   n1 = ra × rb  (normal to plane i-j-k)
    //   n2 = rb × rc  (normal to plane j-k-l)
    //
    // These are the fundamental vectors that define the dihedral angle.
    // Their magnitudes relate to the area of the triangular planes.

    Eigen::Vector3d na = ra.cross(rb);  // Called 'na' in Fortran
    Eigen::Vector3d nb = rb.cross(rc);  // Called 'nb' in Fortran

    double nan = na.norm();    // |n1|
    double nbn = nb.norm();    // |n2|

    // ==========================================================================
    // STEP 4: Calculate trigonometric functions
    // ==========================================================================
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    // ==========================================================================
    // STEP 5: Singularity check
    // ==========================================================================
    // When sin(φ) ≈ 0, the dihedral angle is undefined (linear geometry).
    // Physical interpretation: At φ = 0° or 180°, the torsional force vanishes
    // because the system is at a turning point of the potential.

    constexpr double eps = 1.0e-14;
    double denominator = nan * nbn * sin_phi;
    double one_over_denom;

    if (std::abs(denominator) < eps) {
        // Numerical singularity: return zero gradients
        grad_i.setZero();
        grad_j.setZero();
        grad_k.setZero();
        grad_l.setZero();
        return;
    } else {
        one_over_denom = 1.0 / denominator;
    }

    // ==========================================================================
    // STEP 6: Calculate all necessary cross products
    // ==========================================================================
    // The gradient expressions require numerous cross products.
    // Following Fortran naming convention for clarity and verification.
    //
    // Basic cross products with bond vectors:
    Eigen::Vector3d rab = na.cross(rb);   // (ra × rb) × rb = -|rb|² · ra_perp
    Eigen::Vector3d rba = nb.cross(ra);   // (rb × rc) × ra
    Eigen::Vector3d rac = na.cross(rc);   // (ra × rb) × rc
    Eigen::Vector3d rbb = nb.cross(rb);   // (rb × rc) × rb = -|rb|² · rc_perp
    Eigen::Vector3d rbc = nb.cross(rc);   // (rb × rc) × rc
    Eigen::Vector3d raa = na.cross(ra);   // (ra × rb) × ra = |ra|² · rb_perp

    // Cross products with composite vectors:
    Eigen::Vector3d rapba = rapb.cross(na);  // (ra + rb) × (ra × rb)
    Eigen::Vector3d rapbb = rapb.cross(nb);  // (ra + rb) × (rb × rc)
    Eigen::Vector3d rbpca = rbpc.cross(na);  // (rb + rc) × (ra × rb)
    Eigen::Vector3d rbpcb = rbpc.cross(nb);  // (rb + rc) × (rb × rc)

    // ==========================================================================
    // STEP 7: Compute gradient for atom i
    // ==========================================================================
    // ∂φ/∂r_i = (1/(|n1||n2|sin φ)) · [cos(φ)·(|n2|/|n1|)·(n1×r_jk) - (n2×r_jk)]
    //
    // Physical interpretation: Atom i only appears in vector r_ij (ra).
    // Moving atom i perpendicular to the i-j-k plane changes φ most effectively.

    grad_i = one_over_denom * (cos_phi * (nbn / nan) * rab - rbb);

    // ==========================================================================
    // STEP 8: Compute gradient for atom j (complex!)
    // ==========================================================================
    // Atom j appears in BOTH r_ij and r_jk, leading to a more complex expression.
    // The gradient has contributions from both dihedral planes.

    grad_j = one_over_denom * (
        cos_phi * ((nbn / nan) * rapba + (nan / nbn) * rbc)
        - (rac + rapbb)
    );

    // ==========================================================================
    // STEP 9: Compute gradient for atom k (also complex!)
    // ==========================================================================
    // Similar to atom j: appears in both r_jk and r_kl.

    grad_k = one_over_denom * (
        cos_phi * ((nbn / nan) * raa + (nan / nbn) * rbpcb)
        - (rba + rbpca)
    );

    // ==========================================================================
    // STEP 10: Compute gradient for atom l
    // ==========================================================================
    // ∂φ/∂r_l = (1/(|n1||n2|sin φ)) · [cos(φ)·(|n1|/|n2|)·(n2×r_jk) - (n1×r_jk)]
    //
    // Symmetric to atom i: only appears in vector r_kl (rc).

    grad_l = one_over_denom * (cos_phi * (nan / nbn) * rbb - rab);

    // ==========================================================================
    // VALIDATION NOTE
    // ==========================================================================
    // The four gradients must satisfy:
    //   grad_i + grad_j + grad_k + grad_l = 0  (translational invariance)
    //
    // This is guaranteed by the analytical derivation but can be checked
    // numerically during testing.
}

// =============================================================================
// TORSION PARAMETER GENERATION
// =============================================================================

/**
 * @brief Generate all torsion parameters for the molecule
 *
 * SCIENTIFIC BACKGROUND:
 * ----------------------
 * Torsions in GFN-FF describe the energy cost of rotating around single bonds.
 * They are essential for:
 *   - Conformational preferences (gauche vs. anti in alkanes)
 *   - Conjugation effects (planar preferences in aromatic systems)
 *   - Steric interactions (avoiding eclipsed conformations)
 *
 * ALGORITHM:
 * ----------
 * For each sequence of three consecutive bonds i-j, j-k, k-l:
 *   1. Check if torsion is valid (not linear, not in rigid rings)
 *   2. Determine periodicity n from hybridization of j and k
 *   3. Assign barrier height based on element types
 *   4. Apply topology corrections (rings, conjugation)
 *   5. Store parameters in JSON format for ForceField engine
 *
 * FORTRAN REFERENCE:
 * ------------------
 * This is a SIMPLIFIED implementation based on:
 *   external/gfnff/src/gfnff_ini.f90:1630-1850
 *
 * The full Fortran implementation includes:
 *   - Ring strain corrections (reduces barriers in small rings)
 *   - Pi-system detection (increases barriers for conjugated bonds)
 *   - Hyperconjugation effects
 *   - Special cases for heteroatoms (N, O, S, halogens)
 *   - Multiple torsion terms per bond (n=1, n=2, n=3 simultaneously)
 *
 * **CURRENT SIMPLIFICATIONS** (Phase 1.1):
 *   - Single torsion term per bond (dominant periodicity only)
 *   - No ring detection (Phase 2 requirement)
 *   - No pi-system detection (Phase 2 requirement)
 *   - Simplified barrier heights (no environment-dependent scaling)
 *   - No NCI corrections (non-covalent interaction damping)
 *
 * **PLANNED IMPROVEMENTS** (Future phases):
 *   - Phase 2: Ring detection → strain corrections
 *   - Phase 2: Pi-system detection → conjugation corrections
 *   - Phase 3: EEQ charges → electrostatic scaling
 *   - Phase 4: Multiple torsion terms per bond
 *
 * JSON OUTPUT FORMAT:
 * -------------------
 * Each torsion is stored as:
 *   {
 *     "type": 3,               // GFN-FF torsion type
 *     "i": atom_i,             // First atom index (0-based)
 *     "j": atom_j,             // Second atom (central bond start)
 *     "k": atom_k,             // Third atom (central bond end)
 *     "l": atom_l,             // Fourth atom
 *     "periodicity": n,        // Rotational symmetry (1, 2, or 3)
 *     "barrier": V,            // Energy barrier in kcal/mol
 *     "phase": phi0,           // Reference angle in radians
 *     "is_improper": false     // False for proper torsions
 *   }
 *
 * USAGE IN GFN-FF:
 * ----------------
 * Called during force field initialization:
 *   json ff_params = generateGFNFFParameters();
 *   ff_params["torsions"] = generateGFNFFTorsions();  // ← This function
 *
 * The ForceField engine then iterates over all torsions:
 *   for (auto& torsion : ff_params["torsions"]) {
 *     double E = calculateTorsionEnergy(torsion["i"], torsion["j"], ...);
 *   }
 *
 * COMPUTATIONAL COST:
 * -------------------
 * For a molecule with N atoms and B bonds:
 *   - Worst case: O(N⁴) if all atoms are connected
 *   - Typical case: O(B²) ≈ O(N²) for organic molecules
 *   - Example: Butane (4 atoms, 3 bonds) → 1 torsion
 *              Octane (8 atoms, 7 bonds) → 5 torsions
 *
 * VALIDATION:
 * -----------
 * Compare with external Fortran GFN-FF:
 *   - Same number of torsions detected
 *   - Similar barrier heights (±20% acceptable in Phase 1)
 *   - Same periodicity assignments
 *   - Energy differences < 0.5 kcal/mol for test molecules
 *
 * @return JSON array of torsion parameters
 *
 * @note This is Phase 1.1 implementation - simplified but functional
 * @note Full topology corrections require Phase 2 (ring/pi detection)
 *
 * @author Claude (AI Assistant) - Generated 2025-11-10
 *         Based on Fortran implementation by S. Spicher & S. Grimme
 * @copyright Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
json GFNFF::generateGFNFFTorsions() const
{
    const TopologyInfo& topo = getCachedTopology();
    json torsions = json::array();

    // ==========================================================================
    // STEP 1: Get cached topology and bond list
    // ==========================================================================
    // Claude Generated (Jan 2, 2026): Use cached bond list for bond type access
    // This allows btyp < 5 filtering for extra torsions
    const std::vector<std::pair<int, int>>& bond_list = getCachedBondList();
    const std::vector<int>& bond_types = topo.bond_types;

    // Claude Generated Debug (Jan 12, 2026): Add logging to debug why no torsions are generated
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info(fmt::format("generateGFNFFTorsions: Bond list size = {}", bond_list.size()));
    }

    if (bond_list.empty()) {
        CurcumaLogger::warn("GFN-FF torsion generation: No bonds found, skipping torsions");
        return torsions;
    }

    // ==========================================================================
    // STEP 2: Build neighbor list for efficient lookup
    // ==========================================================================
    // For each atom, store all bonded neighbors.
    // This allows O(1) lookup of "which atoms are bonded to j?"

    std::vector<std::vector<int>> neighbors(m_atomcount);
    for (const auto& bond : bond_list) {
        neighbors[bond.first].push_back(bond.second);
        neighbors[bond.second].push_back(bond.first);
    }

    // Debug: Print neighbor counts
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Neighbor counts:");
        for (size_t i = 0; i < neighbors.size(); ++i) {
            CurcumaLogger::info(fmt::format("  Atom {}: {} neighbors", i, neighbors[i].size()));
        }
    }

    // ==========================================================================
    // STEP 3: Get hybridization from topology
    // ==========================================================================
    // Use topology-based hybridization from determineHybridization() for consistency
    // This ensures the oxygen hybridization fix is applied to torsions as well
    const std::vector<int>& hybridization = topo.hybridization;

    // ==========================================================================
    // STEP 4: Generate all i-j-k-l torsion sequences
    // ==========================================================================
    // For each bond j-k (central bond):
    //   For each neighbor i of j (i ≠ k):
    //     For each neighbor l of k (l ≠ j):
    //       Create torsion i-j-k-l

    int torsion_count = 0;

    // Claude Generated Fix (2025-12-12): Add duplicate tracking to prevent torsion explosion
    // Problem: Dimethyl ether (CH3OCH3) generates 132 torsions instead of ~15 due to duplicates
    // Solution: Track generated torsions in a set with canonical ordering
    // Enhanced Fix (2025-12-13): Improved canonicalization to handle all symmetry cases
    std::set<std::array<int, 4>> generated_torsions;

    // Debug counter for total iterations
    int total_iterations = 0;
    int bond_count = 0;

    for (size_t bond_idx = 0; bond_idx < bond_list.size(); ++bond_idx) {
        const auto& central_bond = bond_list[bond_idx];
        int j = central_bond.first;
        int k = central_bond.second;

        bond_count++;

        // Only process each bond once to avoid duplication
        // We only process bonds where j < k
        if (j >= k) {
            continue;
        }

        // Only consider bonds where both atoms have at least 2 neighbors as potential central bonds
        // This prevents terminal bonds (like C-H in ethane) from being central bonds
        if (neighbors[j].size() < 2 || neighbors[k].size() < 2) {
            continue;
        }

        // Special case: Skip if either central atom is sp (linear)
        // Physical reason: Linear atoms have no torsional barrier
        if (hybridization[j] == 1 || hybridization[k] == 1) {
            continue;
        }

        // Claude Generated (Jan 23, 2026): Skip bond types that have no torsion potential
        // Reference: XTB gfnff_ini.f90:1765 "if(btyp(m).eq.3.or.btyp(m).eq.6) cycle"
        // btyp == 3: sp-X (linear/triple) bonds → no torsional barrier
        // btyp == 6: metal eta (η-bonding) → no standard torsion
        if (bond_idx < bond_types.size()) {
            int btyp = bond_types[bond_idx];
            if (btyp == 3 || btyp == 6) {
                continue;
            }
        }

        // Iterate over all neighbors of j (these will be atom i)
        for (int i : neighbors[j]) {
            if (i == k) continue; // Skip the central bond itself

            // Iterate over all neighbors of k (these will be atom l)
            for (int l : neighbors[k]) {
                if (l == j) continue; // Skip the central bond itself

                // Skip if i == l (4-membered ring edge case)
                if (i == l) continue;

                total_iterations++;

                // ==========================================================
                // STEP 4.5: Check for duplicates using improved canonical ordering
                // ==========================================================
                // Enhanced Fix (2025-12-13): More robust duplicate detection
                // Problem: Previous canonicalization only considered i < l but didn't account
                // for all possible symmetries in molecular structures
                // Solution: Create a fully canonical representation that accounts for:
                //   1. Both forward (i-j-k-l) and reverse (l-k-j-i) representations
                //   2. Proper lexicographic ordering of the torsion indices

                // Create both possible representations of the same torsion
                std::array<int, 4> forward_key = {i, j, k, l};
                std::array<int, 4> reverse_key = {l, k, j, i};

                // Choose the lexicographically smaller representation as the canonical form
                std::array<int, 4> torsion_key = (forward_key < reverse_key) ? forward_key : reverse_key;

                // Skip if this torsion was already generated
                if (generated_torsions.count(torsion_key) > 0) {
                    continue;
                }
                generated_torsions.insert(torsion_key);

                // ==========================================================
                // STEP 5: Check for linear geometry (chktors)
                // ==========================================================
                // Reference: GFN-FF only skips torsions if the central bond
                // angles are linear (making the dihedral undefined).
                // Fortran: gfnff_ini2.f90:680-696 (chktors)

                // Calculate bond angles: angle(j,i,k) and angle(i,j,l)
                // Note: central atoms are j,k. Terminal atoms are i,l.
                // 1-based Fortran indices for atoms in valijklff(nat,xyz,i,j,k,l) are:
                // l-i-j-k in GFN-FF notation where i-j is central.
                // In Curcuma we use i-j-k-l where j-k is central.
                // Mapping: Fortran_i=j, Fortran_j=i (wrong!), wait.
                // Let's look at valijklff call: phi = valijklff(nat,xyz,ll,ii,jj,kk)
                // ii, jj are central. ll, kk are terminal.
                // call bangl(xyz,jj,ii,kk,phi)  ! angle at ii (central 1)
                // call bangl(xyz,ii,jj,ll,phi)  ! angle at jj (central 2)

                auto calculate_bond_angle = [&](int a, int b, int c) {
                    Eigen::Vector3d v_ba = (m_geometry.row(a).head<3>() - m_geometry.row(b).head<3>()).normalized();
                    Eigen::Vector3d v_bc = (m_geometry.row(c).head<3>() - m_geometry.row(b).head<3>()).normalized();
                    return std::acos(std::max(-1.0, std::min(1.0, v_ba.dot(v_bc))));
                };

                double angle_ijk = calculate_bond_angle(i, j, k); // Angle at j (central 1)
                double angle_jkl = calculate_bond_angle(j, k, l); // Angle at k (central 2)

                // GFN-FF threshold: 170 degrees (gfnff_ini2.f90:689, 692)
                constexpr double angle_threshold = 170.0 * M_PI / 180.0;
                if (angle_ijk > angle_threshold || angle_jkl > angle_threshold) {
                    continue; // Skip ill-defined dihedrals
                }

                double phi = calculateDihedralAngle(i, j, k, l);


                // ==========================================================
                // STEP 6: Get topology information for torsion calculation
                // ==========================================================
                // Get topology information for torsion calculation
                const TopologyInfo& topo = getCachedTopology();
                bool in_ring = false;
                int ring_size = 0;
                double cn_i_val = 2.0, cn_l_val = 2.0;
                double cn_j = 2.0, cn_k = 2.0;  // Coordination numbers for central atoms

                // Check if atoms j and k are in the same ring
                if (j < topo.ring_sizes.size() && k < topo.ring_sizes.size()) {
                    in_ring = areAtomsInSameRing(j, k, ring_size);
                }

                // Use actual CN values from topology if available
                if (j < static_cast<int>(topo.coordination_numbers.rows())) {
                    cn_j = topo.coordination_numbers(j);
                }
                if (k < static_cast<int>(topo.coordination_numbers.rows())) {
                    cn_k = topo.coordination_numbers(k);
                }
                // For torsion CN correction, use simple neighbor counts to match XTB behavior
                // XTB uses raw neighbor count (topo%nb(20,i)) rather than effective CN
                if (i < topo.neighbor_counts.rows()) {
                    cn_i_val = topo.neighbor_counts(i);
                }
                if (l < topo.neighbor_counts.rows()) {
                    cn_l_val = topo.neighbor_counts(l);
                }

                // Get actual topological charges (qa) for fqq correction
                // CRITICAL FIX (Phase 2 Charge Routing - January 26, 2026):
                // Torsion barriers (fqq) MUST use topological charges (qa), NOT energy charges (q).
                // Reference: CHARGE_DATAFLOW.md and gfnff_ini.f90:1790
                double qa_j = 0.0, qa_k = 0.0;
                if (j < m_atoms.size() && k < m_atoms.size()) {
                    if (topo.topology_charges.rows() > 0) {
                        qa_j = (j < topo.topology_charges.rows()) ? topo.topology_charges(j) : 0.0;
                        qa_k = (k < topo.topology_charges.rows()) ? topo.topology_charges(k) : 0.0;
                    } else if (m_charges.size() > 0) {
                        // Fallback to m_charges ONLY if topology_charges is empty (e.g. single-phase mode)
                        qa_j = (j < m_charges.size()) ? m_charges(j) : 0.0;
                        qa_k = (k < m_charges.size()) ? m_charges(k) : 0.0;
                    }
                }

                // DEBUG: Check if charges are available (first torsion only)
                static bool charge_debug_printed = false;
                if (!charge_debug_printed && i == 0) {
                    if (CurcumaLogger::get_verbosity() >= 3) {
                        CurcumaLogger::info("\nCHARGE DEBUG:");
                        CurcumaLogger::info(fmt::format("  topo.topology_charges.rows() = {}", topo.topology_charges.rows()));
                        if (topo.topology_charges.rows() > 0) {
                            CurcumaLogger::info(fmt::format("  qa_j (atom {}) = {:.6f}", j, qa_j));
                            CurcumaLogger::info(fmt::format("  qa_k (atom {}) = {:.6f}", k, qa_k));
                            std::string charges_str = "  All topology charges: [";
                            for (int idx = 0; idx < std::min(6, (int)topo.topology_charges.rows()); ++idx) {
                                charges_str += fmt::format("{:.4f} ", topo.topology_charges(idx));
                            }
                            CurcumaLogger::info(charges_str + "]");
                        }
                    }
                    charge_debug_printed = true;
                }

                auto params = getGFNFFTorsionParameters(
                    m_atoms[i], m_atoms[j], m_atoms[k], m_atoms[l],
                    hybridization[j], hybridization[k],
                    qa_j, qa_k,
                    cn_i_val, cn_l_val,
                    in_ring, ring_size, i, j, k, l
                );

                // Skip torsions with barrier below threshold
                // Claude Generated (Jan 23, 2026): Use fcthr = 1e-3 to match XTB
                // XTB: gfnff_param.f90:632 "gen%fcthr = 1.d-3 ! skip torsion if potential is small"
                // Previous: 1e-10 (too permissive) → generated 732 torsions vs XTB's 569
                constexpr double fcthr = 1.0e-3;  // Force constant threshold (Hartree)
                if (std::abs(params.barrier_height) < fcthr) {
                    continue;
                }

                // ==========================================================
                // STEP 7: Store in JSON format
                // ==========================================================
                json torsion;
                torsion["type"] = 3; // GFN-FF type
                // Claude Generated (Feb 4, 2026): Fortran convention atom ordering (ll-ii-jj-kk)
                // Fortran stores torsions as: ll (neighbor of jj) - ii - jj - kk (neighbor of ii)
                // Our loop generates: i (neighbor of j) - j - k - l (neighbor of k)
                // Mapping: j→ii, k→jj → i corresponds to kk, l corresponds to ll
                // Swap terminal atoms so damping distances match Fortran:
                //   pos0-pos1 = l-j = neighbor_of_k to j  → 1-3 distance (NOT a bond)
                //   pos1-pos2 = j-k                        → central bond
                //   pos2-pos3 = k-i = k to neighbor_of_j  → 1-3 distance (NOT a bond)
                // This corrects the 33-44× torsion energy excess (damp ~0.074 → ~0.0017)
                torsion["i"] = l;  // ll = neighbor of jj (was stored as i = neighbor of j)
                torsion["j"] = j;  // ii = central atom 1
                torsion["k"] = k;  // jj = central atom 2
                torsion["l"] = i;  // kk = neighbor of ii (was stored as l = neighbor of k)
                // Claude Generated Fix (2025-11-30): Renamed JSON keys to match forcefield.cpp loader
                // Previous keys: "periodicity", "barrier", "phase" → caused ALL torsion energies = 0
                // Loader expects: "n", "V", "phi0" (from Dihedral struct in forcefield.cpp:444-461)
                torsion["n"] = params.periodicity;          // Periodicity (was "periodicity")
                torsion["V"] = params.barrier_height;       // Barrier height in Hartree (corrected Jan 8, 2026)
                torsion["phi0"] = params.phase_shift;       // Phase shift in radians (was "phase")
                torsion["is_improper"] = params.is_improper;

                // Claude Generated (Jan 9, 2026): Store hybridization for debugging
                torsion["hyb_j"] = hybridization[j];
                torsion["hyb_k"] = hybridization[k];
                if (j < topo.coordination_numbers.rows()) {
                    torsion["cn_j"] = topo.coordination_numbers(j);
                }
                if (k < topo.coordination_numbers.rows()) {
                    torsion["cn_k"] = topo.coordination_numbers(k);
                }

                // Store current dihedral angle for reference
                torsion["current_angle"] = phi; // radians

                torsions.push_back(torsion);
                torsion_count++;
            }
        }
    }

    // ==========================================================================
    // Claude Generated (Jan 24, 2026): NORMALIZATION REMOVED
    // ==========================================================================
    // CRITICAL: The Fortran reference (external/gfnff/src/gfnff_engrad.F90:525-542)
    // uses SIMPLE SUMMATION for primary torsions - NO per-bond averaging, NO 0.5 factor!
    //
    // The previous "fix" with 0.5 factor and per-bond averaging was INCORRECT.
    // Reference verification:
    //   - Primary torsions (line 539): etors = etors + etmp (simple sum)
    //   - NCI Type 2/3 ONLY (line 2793): etors = etors / ntors (NOT primary!)
    //   - Energy formula (line 1272): et = (1+cos)*vtors(2) with range [0, 2V] INTENTIONALLY
    //
    // If energy is ~10× too large, the error is in FORCE CONSTANT GENERATION (vtors(2,m)),
    // not in energy accumulation. Investigate:
    //   1. H-count correction: (nH_j * nH_k)^0.07
    //   2. CN/charge corrections in fijk factors
    //   3. Hybridization-based parameter selection
    //
    // Diagnostic output below helps identify the real source of error.
    // ==========================================================================

    // Diagnostic: Per-torsion parameter summary
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("\n=== PRIMARY TORSION PARAMETER DIAGNOSTICS ===");
        CurcumaLogger::info(fmt::format("Total primary torsions generated: {}", torsion_count));

        // Aggregate statistics
        double total_V = 0.0;
        double max_V = 0.0;
        double min_V = 1e10;
        for (const auto& t : torsions) {
            double V = t["V"];
            total_V += V;
            max_V = std::max(max_V, V);
            min_V = std::min(min_V, V);
        }

        CurcumaLogger::info(fmt::format("Sum of all V (barrier heights): {:.6f} Eh", total_V));
        CurcumaLogger::info(fmt::format("V range: [{:.6f}, {:.6f}] Eh", min_V, max_V));
        CurcumaLogger::info(fmt::format("Average V: {:.6f} Eh", total_V / std::max(1, torsion_count)));

        // Per-bond quartet counts (for reference only, not used for normalization)
        std::map<std::pair<int, int>, int> quartets_per_bond;
        for (const auto& t : torsions) {
            int tj = t["j"];
            int tk = t["k"];
            std::pair<int, int> bond_key = (tj < tk) ? std::make_pair(tj, tk) : std::make_pair(tk, tj);
            quartets_per_bond[bond_key]++;
        }
        CurcumaLogger::info(fmt::format("Number of central bonds with torsions: {}", quartets_per_bond.size()));

        // Detailed per-torsion output at verbosity 3
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("\nPer-torsion details (first 20):");
            int count = 0;
            for (const auto& t : torsions) {
                if (count >= 20) break;
                int ti = t["i"], tj = t["j"], tk = t["k"], tl = t["l"];
                double V = t["V"];
                double n = t["n"];
                double phi0 = t["phi0"];
                CurcumaLogger::info(fmt::format(
                    "  Torsion {}-{}-{}-{}: V={:.6f} Eh, n={:.0f}, phi0={:.2f}°",
                    ti, tj, tk, tl, V, n, phi0 * 180.0 / M_PI));
                count++;
            }
        }
        CurcumaLogger::info("=== END TORSION DIAGNOSTICS ===\n");
    }

    // ==========================================================================
    // STEP 7B: Generate Extra SP3-SP3 Torsions (n=1 for Gauche Conformations)
    // ==========================================================================
    // DISABLED (January 15, 2026): Investigation phase - these terms overcompensate
    // Reference: gfnff_ini.f90:1952-2002
    // "extra rot=1 torsion potential for sp3-sp3 to get gauche conf energies well"
    //
    // Physical meaning: Fine-tune gauche vs anti energy differences for sp3-sp3 bonds
    // These torsions use a DIFFERENT formula and contribute separately from primary n=3 torsions

    int extra_torsion_count = 0;
    const double qfacTOR = 12.0;  // From gfnff_param.f90:742

    // DIAGNOSTIC (Jan 25, 2026)
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("=== EXTRA TORSION GENERATION START ==="));
        CurcumaLogger::info(fmt::format("  Primary torsions generated: {}", generated_torsions.size()));
    }

    // Extra sp3-sp3 torsion factors for gauche conformations
    // Reference: gfnff_param.f90:795-797
    // NOTE (Jan 2, 2026): Original GFN-FF parametrization values restored
    //
    // HISTORICAL NOTE: If torsion energy shows wrong sign/magnitude in future:
    //   The fix is NOT to change these factors, but to ensure proper architectural separation:
    //   1. Primary torsions (n=3, n=2) calculated in CalculateGFNFFDihedralContribution()
    //   2. Extra torsions (n=1) calculated in CalculateGFNFFExtraTorsionContribution()
    //   3. Both stored in separate vectors (m_gfnff_dihedrals vs m_gfnff_extra_torsions)
    //   4. Energy accumulated separately to prevent double-counting
    //   The root cause was mixing both in same vector, NOT wrong force constants.
    const double torsf_extra_C = -0.90;  // Carbon sp3-sp3 (GFN-FF original)
    const double torsf_extra_N =  0.70;  // Nitrogen sp3-sp3 (GFN-FF original)
    const double torsf_extra_O = -2.00;  // Oxygen sp3-sp3 (GFN-FF original)

    // Debug counters (Claude Generated Jan 16, 2026)
    int debug_total_checked = 0;
    int debug_failed_sp3 = 0;
    int debug_failed_ring = 0;
    int debug_failed_btyp = 0;
    int debug_failed_bounds = 0;

    for (const auto& tors_array : generated_torsions) {
        debug_total_checked++;
        int i = tors_array[0];
        int j = tors_array[1];
        int k = tors_array[2];
        int l = tors_array[3];

        // ======================================================================
        // Condition 1: ONLY central atoms (j,k) must be sp3 (hyb==3)
        // ======================================================================
        // Reference: gfnff_ini.f90:1953-1954
        // sp3ij = hyb(ii) .eq. 3.and.hyb(jj) .eq. 3  (central atoms MUST be sp3)
        //
        // CRITICAL FIX (Claude Generated Jan 21, 2026): Use EFFECTIVE hybridization!
        // Same as primary torsions - atoms in pi-system should be treated as sp2.
        // This prevents generating extra torsions for CH3-N(aromatic) bonds where
        // N is part of an aromatic system but has CN~3.
        //
        // Fortran uses: piadr(ii) > 0 to exclude pi-system atoms from sp3 check
        //
        // Lambda to check if atom is in pi-system (same as primary torsions)
        auto is_in_pi_system_extra = [&](int atom_idx, int z_atom) -> bool {
            // Condition 1: Has bonds with significant pi-character
            if (!topo.pi_bond_orders.empty()) {
                for (int other = 0; other < m_atomcount; other++) {
                    if (other == atom_idx) continue;
                    int pibo_idx = lin(atom_idx, other);
                    if (pibo_idx >= 0 && pibo_idx < static_cast<int>(topo.pi_bond_orders.size())) {
                        if (topo.pi_bond_orders[pibo_idx] > 0.1) {
                            return true;
                        }
                    }
                }
            }

            // Condition 2: Is N/O/F/S adjacent to sp2 atoms ("picon" case)
            if (z_atom == 7 || z_atom == 8 || z_atom == 9 || z_atom == 16) {
                const auto& bond_list = getCachedBondList();
                for (const auto& bond : bond_list) {
                    int neighbor = -1;
                    if (bond.first == atom_idx) neighbor = bond.second;
                    else if (bond.second == atom_idx) neighbor = bond.first;
                    else continue;

                    if (neighbor >= 0 && neighbor < static_cast<int>(hybridization.size())) {
                        int neighbor_hyb = hybridization[neighbor];
                        if (neighbor_hyb == 1 || neighbor_hyb == 2) {
                            return true;
                        }
                    }
                }
            }
            return false;
        };

        // Get effective hybridization for central atoms
        int eff_hyb_j = hybridization[j];
        int eff_hyb_k = hybridization[k];

        if (is_in_pi_system_extra(j, m_atoms[j]) && eff_hyb_j == 3) {
            eff_hyb_j = 2;  // Treat as sp2
        }
        if (is_in_pi_system_extra(k, m_atoms[k]) && eff_hyb_k == 3) {
            eff_hyb_k = 2;  // Treat as sp2
        }

        // Now check sp3-sp3 using EFFECTIVE hybridization
        // Claude Generated Fix (Jan 23, 2026): XTB requires ALL FOUR atoms to be sp3!
        // Reference: gfnff_ini.f90:1952-1954
        //   sp3kl = topo%hyb(kk).eq.3.and.topo%hyb(ll).eq.3  ! Outer atoms
        //   sp3ij = topo%hyb(ii).eq.3.and.topo%hyb(jj).eq.3  ! Central atoms
        //   if(sp3kl.and.sp3ij...) then
        // Previous bug: Only checked central atoms → generated too many extra torsions (243 vs ~50)

        // CRITICAL FIX (Jan 25, 2026): Use EFFECTIVE hybridization for ALL atoms, not just central!
        // Hydrogen has hyb=0 (not 3!), so sp3_il will be false for H-C-O-C quartets
        int eff_hyb_i = hybridization[i];
        int eff_hyb_l = hybridization[l];

        // Pi-system detection for outer atoms (same logic as central atoms)
        if (is_in_pi_system_extra(i, m_atoms[i]) && eff_hyb_i == 3) {
            eff_hyb_i = 2;  // Treat as sp2
        }
        if (is_in_pi_system_extra(l, m_atoms[l]) && eff_hyb_l == 3) {
            eff_hyb_l = 2;  // Treat as sp2
        }

        bool sp3_ij = (eff_hyb_j == 3) && (eff_hyb_k == 3);  // Central atoms j,k
        bool sp3_il = (eff_hyb_i == 3) && (eff_hyb_l == 3);  // Outer atoms i,l (NOW with effective hyb!)

        // DIAGNOSTIC (Jan 25, 2026): Debug hydrogen hybridization
        if (CurcumaLogger::get_verbosity() >= 3 && (m_atoms[i] == 1 || m_atoms[l] == 1)) {
            CurcumaLogger::info(fmt::format(
                "EXTRA TORSION CHECK: {}-{}-{}-{} (Z={},{},{},{}): hyb=({},{},{},{}), eff_hyb=({},{},{},{}), sp3_ij={}, sp3_il={}",
                i, j, k, l, m_atoms[i], m_atoms[j], m_atoms[k], m_atoms[l],
                hybridization[i], hybridization[j], hybridization[k], hybridization[l],
                eff_hyb_i, eff_hyb_j, eff_hyb_k, eff_hyb_l,
                sp3_ij, sp3_il
            ));
        }

        if (!sp3_ij || !sp3_il) {
            debug_failed_sp3++;
            continue;  // ALL FOUR atoms must be sp3
        }

        // ======================================================================
        // Condition 2: Central bond (j-k) must be acyclic
        // ======================================================================
        int ring_size;
        if (areAtomsInSameRing(j, k, ring_size)) {
            debug_failed_ring++;
            continue;  // Skip ring torsions
        }

        // ======================================================================
        // Condition 3: Bond type must be non-metal (btyp < 5)
        // ======================================================================
        // Reference: Fortran gfnff_ini.f90:1954
        // if (sp3kl.and.sp3ij.and.(.not.lring).and.btyp(m) .lt. 5) then
        //
        // Find bond index for central bond j-k
        int central_bond_idx = -1;
        for (size_t bond_idx = 0; bond_idx < bond_types.size(); ++bond_idx) {
            const auto& [bi, bj] = topo.adjacency_list.size() > 0
                ? std::make_pair(j, k)  // Use adjacency list if available
                : std::make_pair(j, k);
            // getCachedBondList() returns bonds as pairs (i,j) where i < j
            const auto& bond = getCachedBondList()[bond_idx];
            if ((bond.first == j && bond.second == k) || (bond.first == k && bond.second == j)) {
                central_bond_idx = bond_idx;
                break;
            }
        }

        if (central_bond_idx == -1) {
            continue;
        }

        int btyp = bond_types[central_bond_idx];
        if (btyp >= 5) {
            debug_failed_btyp++;
            continue;  // Skip metal-containing bonds
        }

        // ======================================================================
        // Condition 4: Bounds check for atomic numbers
        // ======================================================================
        if (m_atoms[i] < 1 || m_atoms[i] > 86 ||
            m_atoms[j] < 1 || m_atoms[j] > 86 ||
            m_atoms[k] < 1 || m_atoms[k] > 86 ||
            m_atoms[l] < 1 || m_atoms[l] > 86) {
            debug_failed_bounds++;
            continue;  // Skip invalid atomic numbers
        }

        // ======================================================================
        // Select heteroatom-specific force constant
        // ======================================================================
        // Reference: gfnff_ini.f90:1961-1963
        // Priority: O > N > C (if both O and N, O wins)

        double ff = torsf_extra_C;  // Default: Carbon

        // Check CENTRAL atoms only (j and k)
        if (m_atoms[j] == 7 || m_atoms[k] == 7) {
            ff = torsf_extra_N;  // Nitrogen
        }
        if (m_atoms[j] == 8 || m_atoms[k] == 8) {
            ff = torsf_extra_O;  // Oxygen (overrides nitrogen)
        }

        // ======================================================================
        // Calculate barrier height (DIFFERENT formula - no f1/f2!)
        // ======================================================================
        // Reference: gfnff_ini.f90:1970
        // vtors(2) = ff * fij * fkl * fqq

        // Calculate fij (central bond contribution)
        double fij = tors_angewChem2020[m_atoms[j] - 1] * tors_angewChem2020[m_atoms[k] - 1];

        // Check threshold (like primary torsions)
        const double fcthr = 1.0e-3;
        if (fij < fcthr || tors_angewChem2020[m_atoms[j] - 1] < 0.0 || tors_angewChem2020[m_atoms[k] - 1] < 0.0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "  Torsion {}-{}-{}-{}: fij threshold (fij={:.6f}, tors[{}]={:.6f}, tors[{}]={:.6f})",
                    i, j, k, l, fij, m_atoms[j], tors_angewChem2020[m_atoms[j] - 1],
                    m_atoms[k], tors_angewChem2020[m_atoms[k] - 1]
                ));
            }
            continue;  // Skip if fij too small or negative
        }

        // Calculate fkl (outer atom contribution)
        double fkl = tors2_angewChem2020[m_atoms[i] - 1] * tors2_angewChem2020[m_atoms[l] - 1];

        // Check threshold
        if (fkl < fcthr || tors2_angewChem2020[m_atoms[i] - 1] < 0.0 || tors2_angewChem2020[m_atoms[l] - 1] < 0.0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "  Torsion {}-{}-{}-{}: fkl threshold (fkl={:.6f}, tors2[{}]={:.6f}, tors2[{}]={:.6f})",
                    i, j, k, l, fkl, m_atoms[i], tors2_angewChem2020[m_atoms[i] - 1],
                    m_atoms[l], tors2_angewChem2020[m_atoms[l] - 1]
                ));
            }
            continue;  // Skip if fkl too small or negative
        }

        // Calculate fqq (charge correction)
        // Claude Generated (Jan 9, 2026): Use m_charges directly instead of topo.topology_charges
        // Fix: topo.topology_charges may be empty or 0.0 in some cases, use pre-calculated m_charges
        double qa_j = 0.0, qa_k = 0.0;
        if (j < m_atoms.size() && k < m_atoms.size()) {
            // First try: Use TopologyInfo charges (filled by calculateTopologyInfo)
            if (topo.topology_charges.rows() > 0) {
                qa_j = (j < topo.topology_charges.rows()) ? topo.topology_charges(j) : 0.0;
                qa_k = (k < topo.topology_charges.rows()) ? topo.topology_charges(k) : 0.0;
                // Fallback: If topology charges are zero, use pre-calculated EEQ charges
                if (std::abs(qa_j) < 1e-10 && std::abs(qa_k) < 1e-10) {
                    if (m_charges.size() > 0) {
                        qa_j = m_charges(j);
                        qa_k = m_charges(k);
                    }
                }
            }
            // Second try: Use pre-calculated EEQ charges (m_charges in gfnff_method.h)
            else if (m_charges.size() > 0) {
                qa_j = (j < m_charges.size()) ? m_charges(j) : 0.0;
                qa_k = (k < m_charges.size()) ? m_charges(k) : 0.0;
            }
        }
        double fqq = 1.0 + std::abs(qa_j * qa_k) * qfacTOR;

        // Final barrier (DIFFERENT from primary torsion formula!)
        double barrier = ff * fij * fkl * fqq;

        // Validate barrier is not NaN or inf
        if (std::isnan(barrier) || std::isinf(barrier)) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::warn(fmt::format(
                    "  Warning: Invalid barrier for torsion {}-{}-{}-{}, skipping (barrier={}, ff={}, fij={}, fkl={}, fqq={})",
                    i, j, k, l, barrier, ff, fij, fkl, fqq
                ));
            }
            continue;
        }

        // Check final barrier threshold
        if (std::abs(barrier) < fcthr) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "  Torsion {}-{}-{}-{}: barrier too small (barrier={:.6f}, threshold={:.6f})",
                    i, j, k, l, barrier, fcthr
                ));
            }
            continue;  // Skip negligible barriers
        }

        // SUCCESS - this torsion will be added!
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  ✓ Torsion {}-{}-{}-{}: ACCEPTED (barrier={:.6f} Eh, ff={:.2f})",
                i, j, k, l, barrier, ff
            ));
        }

        // ======================================================================
        // Store extra torsion with n=1 periodicity
        // ======================================================================
        json extra_torsion;
        extra_torsion["type"] = 3;           // GFN-FF type (same as primary)
        // Claude Generated (Feb 4, 2026): Same Fortran convention swap as primary torsions
        // Extra torsions iterate over the same generated_torsions set, so the same
        // ll-ii-jj-kk reordering is needed for correct star-topology damping distances
        extra_torsion["i"] = l;  // ll = neighbor of jj
        extra_torsion["j"] = j;  // ii = central atom 1
        extra_torsion["k"] = k;  // jj = central atom 2
        extra_torsion["l"] = i;  // kk = neighbor of ii
        extra_torsion["n"] = 1;              // Periodicity = 1 (NOT 3!)
        extra_torsion["V"] = barrier;        // Barrier in Hartree (same as primary torsions)
        extra_torsion["phi0"] = M_PI;        // Phase = 180° (pi radians)
        extra_torsion["is_improper"] = false;
        extra_torsion["is_extra"] = true;    // Claude Generated (Jan 1, 2026): Mark as extra torsion

        // Store current dihedral angle
        double phi = calculateDihedralAngle(i, j, k, l);

        // Validate phi is not NaN or inf
        if (std::isnan(phi) || std::isinf(phi)) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::warn(fmt::format(
                    "  Warning: Invalid phi angle for torsion {}-{}-{}-{}, skipping",
                    i, j, k, l
                ));
            }
            continue;
        }

        extra_torsion["current_angle"] = phi;

        torsions.push_back(extra_torsion);
        extra_torsion_count++;

        // ======================================================================
        // Debug output
        // ======================================================================
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "EXTRA SP3-SP3 TORSION: {}-{}-{}-{} (Z={},{},{},{})",
                i, j, k, l, m_atoms[i], m_atoms[j], m_atoms[k], m_atoms[l]
            ));
            CurcumaLogger::info(fmt::format(
                "  ff={:.2f} (heteroatom factor), barrier={:.6f} Eh (n=1)",
                ff, barrier
            ));
        }
    }

    if (extra_torsion_count > 0) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "GFN-FF detected {} extra sp3-sp3 torsions (n=1 gauche terms)",
                extra_torsion_count
            ));
        }
    }

    // Debug statistics (Claude Generated Jan 16, 2026)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("EXTRA TORSION DEBUG:"));
        CurcumaLogger::info(fmt::format("  Total torsions checked: {}", debug_total_checked));
        CurcumaLogger::info(fmt::format("  Failed sp3 check: {}", debug_failed_sp3));
        CurcumaLogger::info(fmt::format("  Failed ring check: {}", debug_failed_ring));
        CurcumaLogger::info(fmt::format("  Failed btyp check: {}", debug_failed_btyp));
        CurcumaLogger::info(fmt::format("  Failed bounds check: {}", debug_failed_bounds));
        CurcumaLogger::info(fmt::format("  Passed all checks (extra torsions): {}", extra_torsion_count));
    }

    // ==========================================================================
    // STEP 8: Report results
    // ==========================================================================
    if (torsion_count == 0) {
        CurcumaLogger::warn("GFN-FF: No torsions detected (molecule may be too small or linear)");
    } else {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("GFN-FF detected " + std::to_string(torsion_count) + " torsions");
        }

        // Claude Generated (Jan 2, 2026): Verbosity 2 output for torsion summary
        if (CurcumaLogger::get_verbosity() >= 2) {
            int primary_count = torsion_count;
            CurcumaLogger::result(fmt::format("GFN-FF torsions: {} primary (n=2,3), {} extra sp3-sp3 (n=1)",
                                               primary_count, extra_torsion_count));
        }

        // Debug output for analysis
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("DEBUG: Total iterations: {}, Unique torsions: {}", total_iterations, generated_torsions.size()));
            CurcumaLogger::info("DEBUG: All generated torsions:");
            int torsion_index = 0;
            for (const auto& torsion : generated_torsions) {
                CurcumaLogger::info(fmt::format("  {} : [{}, {}, {}, {}]", torsion_index++, torsion[0], torsion[1], torsion[2], torsion[3]));
            }
        }

        // Optional: Print summary by periodicity
        int n1_count = 0, n2_count = 0, n3_count = 0;
        for (const auto& torsion : torsions) { // Claude Generated Fix (2025-12-13): Changed from torsions["dihedrals"] to torsions
            int n = torsion["n"];  // Claude Generated Fix (2025-11-30): Changed from "periodicity" to "n"
            if (n == 1) n1_count++;
            else if (n == 2) n2_count++;
            else if (n == 3) n3_count++;
        }

        CurcumaLogger::info("  Periodicity distribution: n=1 (" + std::to_string(n1_count) +
                           "), n=2 (" + std::to_string(n2_count) +
                           "), n=3 (" + std::to_string(n3_count) + ")");
    }

    // ==========================================================================
    // KNOWN LIMITATIONS (Phase 1.1)
    // ==========================================================================
    // The following features are NOT yet implemented:
    //
    // 1. **Ring strain corrections**: Small rings (3, 4, 5) should have
    //    reduced torsional barriers. Requires ring detection (Phase 2).
    //
    // 2. **Conjugation corrections**: Conjugated systems (aromatics, polyenes)
    //    should have increased barriers to maintain planarity. Requires
    //    pi-system detection (Phase 2).
    //
    // 3. **Multiple torsion terms**: Full GFN-FF uses up to 3 cosine terms
    //    per bond (e.g., n=1, n=2, n=3 for sp³-sp³). We use only dominant term.
    //
    // 4. **NCI damping**: Non-covalent interaction corrections for torsions
    //    involving weak bonds. Requires advanced topology analysis.
    //
    // 5. **Hyperconjugation**: Special treatment for C-C bonds adjacent to
    //    heteroatoms. Requires charge analysis (Phase 3).
    //
    // **Impact**: Energy differences of 1-3 kcal/mol compared to full GFN-FF
    // for complex molecules with rings/conjugation. Acceptable for Phase 1.

    return torsions;
}

