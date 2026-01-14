/*
 * <Hückel Solver for GFN-FF π-Bond Order Calculation>
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
 * Claude Generated (January 14, 2026) - Phase 1: Full Hückel Implementation
 *
 * Self-consistent iterative Hückel method for π-bond order calculation.
 * Port from Fortran gfnff_ini.f90:928-1062 and gfnff_qm.f90:39-154.
 *
 * The algorithm performs:
 * 1. Identify π-systems and count π-electrons per atom
 * 2. Build Hückel Hamiltonian with P-dependent off-diagonal coupling
 * 3. Iteratively diagonalize until convergence (max 5 iterations)
 * 4. Extract π-bond orders from density matrix
 *
 * Reference:
 * - Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020, 59, 15665
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <cmath>

/**
 * @brief Triangular indexing for symmetric matrix storage
 *
 * Claude Generated (January 14, 2026)
 * Port from Fortran gfnff_helpers.f90:416-423
 *
 * @param i First index (0-based)
 * @param j Second index (0-based)
 * @return Linear index in triangular array
 */
inline int huckel_lin(int i, int j) {
    int imax = std::max(i, j);
    int imin = std::min(i, j);
    return imin + imax * (imax + 1) / 2;
}

/**
 * @brief Self-consistent iterative Hückel solver for π-bond orders
 *
 * Claude Generated (January 14, 2026)
 *
 * Implements the GFN-FF Hückel method which differs from classical HMO:
 * - P-dependent off-diagonal coupling (prevents over-delocalization)
 * - Charge-dependent diagonal elements
 * - Fermi smearing at 4000K for biradical handling
 * - Element-specific parameters (B, C, N, O, F, S, Cl)
 */
class HuckelSolver {
public:
    /**
     * @brief Default constructor
     */
    HuckelSolver() = default;

    /**
     * @brief Calculate π-bond orders for all bonds using iterative Hückel
     *
     * This is the main entry point. It:
     * 1. Identifies all π-systems from hybridization data
     * 2. For each π-system, performs iterative Hückel calculation
     * 3. Returns π-bond orders in triangular format using huckel_lin(i,j)
     *
     * @param atoms Atomic numbers (Z) for each atom
     * @param hybridization Hybridization states (1=sp, 2=sp², 3=sp³)
     * @param pi_fragments π-system fragment IDs (0 = not in π-system)
     * @param charges EEQ atomic charges
     * @param bonds List of bonded atom pairs (0-based indices)
     * @param distances N×N distance matrix in Bohr
     * @param itag Special atom tags (1=carbene for C, 1=NO₂ for N)
     * @return π-bond orders in triangular format [huckel_lin(i,j)]
     */
    std::vector<double> calculatePiBondOrders(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::vector<int>& pi_fragments,
        const std::vector<double>& charges,
        const std::vector<std::pair<int,int>>& bonds,
        const Eigen::MatrixXd& distances,
        const std::vector<int>& itag = {}
    );

    /**
     * @brief Set verbosity level for debug output
     * @param level 0=silent, 1=minimal, 2=detailed, 3=full debug
     */
    void setVerbosity(int level) { m_verbosity = level; }

private:
    // ========================================
    // Hückel Parameters (from gfnff_param.f90:818-840)
    // ========================================

    // Diagonal elements relative to Carbon (Coulomb integrals)
    // Index by atomic number: hdiag[Z]
    static constexpr double hdiag[18] = {
        0.0,    // 0: placeholder
        0.0,    // 1: H (not in π-system)
        0.0,    // 2: He
        0.0,    // 3: Li
        0.0,    // 4: Be
        -0.5,   // 5: B
        0.0,    // 6: C (reference)
        0.14,   // 7: N
        -0.38,  // 8: O
        -0.29,  // 9: F
        0.0,    // 10: Ne
        0.0,    // 11: Na
        0.0,    // 12: Mg
        0.0,    // 13: Al
        0.0,    // 14: Si
        0.0,    // 15: P
        -0.30,  // 16: S
        -0.30   // 17: Cl
    };

    // Off-diagonal constants (exchange integrals, β parameters)
    // Index by atomic number: hoffdiag[Z]
    static constexpr double hoffdiag[18] = {
        0.0,    // 0: placeholder
        0.0,    // 1: H
        0.0,    // 2: He
        0.0,    // 3: Li
        0.0,    // 4: Be
        0.5,    // 5: B
        1.0,    // 6: C (reference)
        0.66,   // 7: N
        1.1,    // 8: O
        0.23,   // 9: F
        0.0,    // 10: Ne
        0.0,    // 11: Na
        0.0,    // 12: Mg
        0.0,    // 13: Al
        0.0,    // 14: Si
        0.0,    // 15: P
        0.6,    // 16: S
        1.0     // 17: Cl
    };

    // Algorithm control parameters
    static constexpr double hiter = 0.700;       // Iteration mixing (β scaling reduction)
    static constexpr double htriple = 1.45;      // Triple bond β reduction factor
    static constexpr double hueckelp3 = -0.24;   // Diagonal charge dependence
    static constexpr double pilpf = 0.530;       // Lone pair diagonal shift (2e⁻ systems)
    static constexpr int maxhiter = 5;           // Maximum iterations (prevents divergence)
    static constexpr double fermi_temp = 4000.0; // Electronic temperature for Fermi smearing (K)
    static constexpr double conv_threshold = 1e-4; // Energy convergence criterion

    // Physical constants
    static constexpr double boltz_ev = 8.617333262e-5; // Boltzmann constant in eV/K
    static constexpr double hartree_to_ev = 27.2113957; // Hartree to eV conversion

    // ========================================
    // Helper Methods
    // ========================================

    /**
     * @brief Count π-electrons for an atom based on type and hybridization
     *
     * Port from gfnff_ini.f90:959-978
     *
     * @param atom_type Atomic number
     * @param hyb Hybridization (1=sp, 2=sp², 3=sp³)
     * @param tag Special tag (1=carbene for C, 1=NO₂ for N)
     * @return Number of π-electrons contributed by this atom
     */
    int countPiElectrons(int atom_type, int hyb, int tag = 0) const;

    /**
     * @brief Build Hückel Hamiltonian matrix for a π-system
     *
     * Port from gfnff_ini.f90:989-1009
     *
     * H_ii = hdiag[Z] + q*hueckelp3 - (nel-1)*pilpf
     * H_ij = -β * (1 - hiter*(2/3 - P_old_ij))
     *
     * @param pi_atoms Indices of atoms in this π-system (in original numbering)
     * @param pi_atom_map Reverse map: original atom → π-system index (-1 if not in system)
     * @param pi_electrons Number of π-electrons per atom
     * @param atom_types Atomic numbers
     * @param hybridization Hybridization states
     * @param charges Atomic charges
     * @param bonds Bond list
     * @param distances Distance matrix
     * @param P_old Previous density matrix (for P-dependent coupling)
     * @return Hamiltonian matrix
     */
    Eigen::MatrixXd buildHamiltonian(
        const std::vector<int>& pi_atoms,
        const std::vector<int>& pi_atom_map,
        const std::vector<int>& pi_electrons,
        const std::vector<int>& atom_types,
        const std::vector<int>& hybridization,
        const std::vector<double>& charges,
        const std::vector<std::pair<int,int>>& bonds,
        const Eigen::MatrixXd& distances,
        const Eigen::MatrixXd& P_old
    ) const;

    /**
     * @brief Solve eigenvalue problem and compute density matrix
     *
     * Port from gfnff_qm.f90:39-154
     *
     * Steps:
     * 1. Diagonalize H → eigenvalues ε, eigenvectors C
     * 2. Scale energies: ε *= 0.1 * 27.2113957 (→ eV)
     * 3. Compute occupations via Fermi smearing at 4000K
     * 4. Build density matrix: P = C * diag(occ) * C^T
     * 5. Return electronic energy: E = Σ occ_i * ε_i
     *
     * @param H Hamiltonian matrix (input), density matrix (output)
     * @param nel Number of electrons
     * @return Electronic energy (for convergence check)
     */
    double solveAndBuildDensity(Eigen::MatrixXd& H, int nel) const;

    /**
     * @brief Compute Fermi-Dirac occupations at given temperature
     *
     * Port from gfnff_qm.f90:157-219
     *
     * @param eigenvalues Orbital energies in eV
     * @param nel Number of electrons
     * @param temp Electronic temperature in K
     * @return Occupation numbers
     */
    std::vector<double> fermiSmear(
        const Eigen::VectorXd& eigenvalues,
        int nel,
        double temp
    ) const;

    /**
     * @brief Compute density matrix from MO coefficients and occupations
     *
     * Port from gfnff_qm.f90:281-304
     *
     * P_μν = Σ_i n_i * C_μi * C_νi
     *
     * @param C MO coefficient matrix (columns are MOs)
     * @param occ Occupation numbers
     * @return Density matrix
     */
    Eigen::MatrixXd computeDensityMatrix(
        const Eigen::MatrixXd& C,
        const std::vector<double>& occ
    ) const;

    // Verbosity level
    int m_verbosity = 0;
};
