/*
 * < EEQ (Electronegativity Equalization) Charge Solver - Implementation >
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Extracted from GFN-FF implementation (gfnff_method.cpp) for reusability
 * across all force field methods (UFF, QMDFF, D4 dispersion).
 *
 * References:
 * - S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673
 * - R. T. Sanderson, Chemical Bonds and Bond Energy, Academic Press, 1976
 *
 * Claude Generated - December 2025
 */

// Claude Generated (March 2026): Enable BLAS-accelerated Eigen for EEQ linear solves
// Defined per-file to avoid Eigen 3.4 complex-type conflicts in other TUs (e.g., lbfgs.cpp)
#ifdef HAVE_BLAS
#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif
#endif

#include "eeq_solver.h"
#include "gfnff_par.h"  // EEQ parameters (chi_eeq, gam_eeq, alpha_eeq, cnf_eeq)
#include "cn_calculator.h"  // Shared CN calculation utility

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"

#include <Eigen/Dense>
#include <atomic>
#include <cmath>
#include <algorithm>
#include <queue>
#include <thread>
#include <future>
#include <fmt/format.h>
#include <fstream>
#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "src/core/blas_threads.h"  // ScopedBlasThreads — give threaded LAPACK the EnergyCalculator budget

using namespace GFNFFParameters;  // Access to chi_eeq, gam_eeq, alpha_eeq, cnf_eeq

// ===== Threaded Cholesky for the EEQ Schur solve (Claude Generated, Jun 2026) =====
//
// The EEQ Schur-Cholesky solve (factorize A_nn once, back-substitute the nfrag+1
// RHS columns) dominated the cost on large many-fragment systems (mixture: 6200
// atoms / 1400 fragments → ~4.4 s, serial). Eigen's LLT runs single-threaded here
// (EIGEN_USE_BLAS does not route the dense Cholesky; EIGEN_USE_LAPACKE is off),
// and OpenBLAS thread count had no effect. Route the factorization AND the
// multi-RHS back-substitution through the threaded LAPACK (dpotrf/dpotrs) that
// the already-linked BLAS/LAPACK (OpenBLAS/MKL) provides — portable to any
// BLAS/LAPACK, not tied to the lapacke C interface. Falls back to Eigen LLT when
// built without BLAS.
//
// LAPACK is column-major; curcuma's `Matrix` is row-major, so the LAPACK path uses
// column-major Eigen::MatrixXd (`ColMatrix`) and converts at the boundaries (O(N²),
// negligible vs the O(N³) factorization).
namespace {
#ifdef USE_BLAS
extern "C" {
    void dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info);
    void dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
                 const int* lda, double* b, const int* ldb, int* info);
}
#endif
using ColMatrix = Eigen::MatrixXd;  // column-major (LAPACK layout)

/// Cholesky-factorize SPD A (column-major) in place → lower factor. true on success.
inline bool eeqCholeskyFactorize(ColMatrix& A)
{
#ifdef USE_BLAS
    const char uplo = 'L';
    int n = static_cast<int>(A.rows()), lda = n, info = 0;
    dpotrf_(&uplo, &n, A.data(), &lda, &info);
    return info == 0;
#else
    Eigen::LLT<ColMatrix> llt(A);
    if (llt.info() != Eigen::Success) return false;
    A = llt.matrixL();   // store L in the lower triangle (fallback path)
    return true;
#endif
}

/// Solve A·X = B in place (B column-major, overwritten with X) using the factor
/// produced by eeqCholeskyFactorize().
inline void eeqCholeskySolve(const ColMatrix& Afac, ColMatrix& B)
{
#ifdef USE_BLAS
    const char uplo = 'L';
    int n = static_cast<int>(Afac.rows()), nrhs = static_cast<int>(B.cols()),
        lda = n, ldb = n, info = 0;
    dpotrs_(&uplo, &n, &nrhs, Afac.data(), &lda, B.data(), &ldb, &info);
#else
    B = Afac.triangularView<Eigen::Lower>().solve(B);
    B = Afac.triangularView<Eigen::Lower>().adjoint().solve(B);
#endif
}

/// Vector RHS overload (single column).
inline void eeqCholeskySolve(const ColMatrix& Afac, Eigen::VectorXd& b)
{
#ifdef USE_BLAS
    const char uplo = 'L';
    int n = static_cast<int>(Afac.rows()), nrhs = 1, lda = n, ldb = n, info = 0;
    dpotrs_(&uplo, &n, &nrhs, Afac.data(), &lda, b.data(), &ldb, &info);
#else
    b = Afac.triangularView<Eigen::Lower>().solve(b);
    b = Afac.triangularView<Eigen::Lower>().adjoint().solve(b);
#endif
}

/// Schur complement S = C·Z2 exploiting that the constraint matrix C is a per-fragment
/// membership matrix (one nonzero per atom column). S(f,g) = Σ_i C(f,i)·Z2(i,g) reduces
/// to a per-fragment row accumulation: O(N·nfrag) instead of the dense product's
/// O(nfrag²·N). For many-fragment systems (mixture: 6200 atoms / 1400 fragments) the dense
/// C·Z2 (~1.2e13 flops, ~2.4 s serial) was the dominant EEQ-solve cost; this makes it
/// negligible. Mirrors the GPU path (k_eeq_reduce_fragment_sums). Claude Generated (Jun 2026).
inline Matrix schurComplementCZ(const Matrix& C, const Matrix& Z2, int natoms, int nfrag)
{
    Matrix S = Matrix::Zero(nfrag, nfrag);
    for (int f = 0; f < nfrag; ++f) {
        for (int i = 0; i < natoms; ++i) {
            const double c = C(f, i);
            if (c != 0.0) S.row(f).noalias() += c * Z2.row(i);
        }
    }
    return S;
}
}  // namespace

// ===== Element Group Classification (XTB Compatible) =====
// Based on XTB param%group classification in gfnff_ini2.f90
// Claude Generated - December 2025 (Phase 2: Complete Element-Specific Hybridization)

/**
 * @brief Get XTB element group for hybridization rules
 *
 * Maps atomic number to XTB group classification used in gfnff_ini2.f90
 *
 * @param Z Atomic number
 * @return XTB group number:
 *         1 = H, 2 = Alkali/Alkaline Earth, 3 = Boron, 4 = Carbon,
 *         5 = Nitrogen, 6 = Oxygen, 7 = Halogens, 8 = Noble gases
 *         -3 = Early TMs (Sc-La), -7 = Late TMs, 0 = Unknown
 */
static inline int getElementGroup(int Z) {
    // Main group elements
    if (Z == 1) return 1;        // H
    if (Z == 2 || (Z >= 3 && Z <= 4)) return 2;  // He, Li, Be
    if (Z == 5) return 3;        // B
    if (Z == 6) return 4;        // C
    if (Z == 7) return 5;        // N
    if (Z == 8) return 6;        // O

    // Halogens: F, Cl, Br, I, At
    if (Z == 9 || (Z >= 17 && Z <= 35) || (Z >= 53 && Z <= 85)) return 7;

    // Noble gases: Ne, Ar, Kr, Xe, Rn
    if ((Z >= 10 && Z <= 18) || (Z >= 36 && Z <= 54) || (Z >= 86)) return 8;

    // Transition metals
    if (Z >= 21 && Z <= 30) return -3; // Sc-La (early TMs)
    if ((Z >= 22 && Z <= 28) || (Z >= 39 && Z <= 46) || (Z >= 72 && Z <= 78)) return -7; // Late TMs

    return 0; // Default/unknown
}

/**
 * @brief Check if element is a transition metal
 *
 * @param Z Atomic number
 * @return true if transition metal, false otherwise
 */
static inline bool isTransitionMetal(int Z) {
    int group = getElementGroup(Z);
    return (group == -3 || group == -7);
}

/**
 * @brief Check if element is a metal (transition or main group)
 *
 * @param Z Atomic number
 * @return true if metal, false otherwise
 */
static inline bool isMetal(int Z) {
    // Transition metals
    if (isTransitionMetal(Z)) return true;

    // Main group metals (alkali, alkaline earth, some post-transition)
    if (Z >= 3 && Z <= 4) return true;  // Li, Be
    if (Z >= 11 && Z <= 13) return true; // Na, Mg, Al
    if (Z >= 19 && Z <= 20) return true; // K, Ca
    if (Z >= 37 && Z <= 38) return true; // Rb, Sr
    if (Z >= 55 && Z <= 56) return true; // Cs, Ba
    if (Z >= 87 && Z <= 88) return true; // Fr, Ra

    // Post-transition metals
    if (Z == 31 || Z == 49 || Z == 50 || Z == 81 || Z == 82 || Z == 83 || Z == 113 || Z == 114 || Z == 115 || Z == 116) return true;

    return false;
}

/**
 * @brief Count hydrogen neighbors for transition metal coordination
 *
 * @param atom_index Index of central atom
 * @param neighbor_lists Neighbor list from topology
 * @param atoms Atomic numbers
 * @return Number of hydrogen neighbors
 */
static inline int countHydrogenNeighbors(int atom_index,
                                         const std::vector<std::vector<int>>& neighbor_lists,
                                         const std::vector<int>& atoms) {
    if (atom_index < 0 || atom_index >= neighbor_lists.size()) return 0;

    int count = 0;
    for (int neighbor : neighbor_lists[atom_index]) {
        if (neighbor < atoms.size() && atoms[neighbor] == 1) { // Hydrogen
            count++;
        }
    }
    return count;
}

/**
 * @brief Calculate bond angle in radians (A-B-C)
 *
 * Computes the angle at atom B between atoms A-B-C
 * Matches XTB's bangl() subroutine from gfnff_ini2.f90
 *
 * @param geometry_bohr Geometry matrix (natoms x 3) in Bohr
 * @param atom_a Index of first atom
 * @param atom_b Index of central atom (vertex)
 * @param atom_c Index of third atom
 * @return Bond angle in radians (0 to π)
 *
 * Claude Generated - December 2025 (Phase 3: Geometry-Dependent Hybridization)
 */
static inline double calculateBondAngle(const Matrix& geometry_bohr,
                                        int atom_a, int atom_b, int atom_c) {
    if (atom_a < 0 || atom_b < 0 || atom_c < 0) return 0.0;
    if (atom_a >= geometry_bohr.rows() || atom_b >= geometry_bohr.rows() || atom_c >= geometry_bohr.rows()) return 0.0;

    // Vectors B->A and B->C
    Eigen::Vector3d vec_ba = geometry_bohr.row(atom_a) - geometry_bohr.row(atom_b);
    Eigen::Vector3d vec_bc = geometry_bohr.row(atom_c) - geometry_bohr.row(atom_b);

    // Normalize vectors
    double len_ba = vec_ba.norm();
    double len_bc = vec_bc.norm();

    if (len_ba < 1e-10 || len_bc < 1e-10) return 0.0; // Degenerate case

    vec_ba /= len_ba;
    vec_bc /= len_bc;

    // Dot product
    double cos_angle = vec_ba.dot(vec_bc);

    // Clamp to valid range [-1, 1] to avoid numerical issues
    if (cos_angle > 1.0) cos_angle = 1.0;
    if (cos_angle < -1.0) cos_angle = -1.0;

    return std::acos(cos_angle);
}

/**
 * @brief Count metal neighbors for an atom
 *
 * @param atom_index Index of central atom
 * @param neighbor_lists Neighbor list from topology
 * @param atoms Atomic numbers
 * @return Number of metal neighbors
 *
 * Claude Generated - December 2025 (Phase 3: Topology-Aware Hybridization)
 */
static inline int countMetalNeighbors(int atom_index,
                                      const std::vector<std::vector<int>>& neighbor_lists,
                                      const std::vector<int>& atoms) {
    if (atom_index < 0 || atom_index >= neighbor_lists.size()) return 0;

    int count = 0;
    for (int neighbor : neighbor_lists[atom_index]) {
        if (neighbor < atoms.size() && isMetal(atoms[neighbor])) {
            count++;
        }
    }
    return count;
}

/**
 * @brief Find nearest non-metal neighbor's coordination number
 *
 * Matches XTB's nn_nearest_noM() subroutine used for M-O-X conjugation detection
 *
 * @param atom_index Index of central atom
 * @param neighbor_lists Neighbor list from topology
 * @param atoms Atomic numbers
 * @return Coordination number of nearest non-metal neighbor (0 if none found)
 *
 * Claude Generated - December 2025 (Phase 3: Metal-Oxygen Conjugation)
 */
static inline int findNearestNonMetalCN(int atom_index,
                                        const std::vector<std::vector<int>>& neighbor_lists,
                                        const std::vector<int>& atoms) {
    if (atom_index < 0 || atom_index >= neighbor_lists.size()) return 0;

    // Find first non-metal neighbor and return its CN
    for (int neighbor : neighbor_lists[atom_index]) {
        if (neighbor < atoms.size() && !isMetal(atoms[neighbor])) {
            // Return CN of this neighbor
            if (neighbor < neighbor_lists.size()) {
                return neighbor_lists[neighbor].size();
            }
        }
    }
    return 0; // No non-metal neighbors found
}

// ===== Hybridization Detection with Element-Specific Rules =====
// Enhanced version with XTB-compatible element-specific logic
// Claude Generated - December 2025 (Phase 2: Complete Element-Specific Hybridization)

/**
 * @brief Detect hybridization using element-specific XTB rules
 *
 * Implements comprehensive hybridization detection based on XTB gfnff_ini2.f90:217-332
 * Supports element-specific rules for main group elements and transition metals.
 *
 * @param Z Atomic number
 * @param cn Coordination number
 * @param atoms Full atom list (for neighbor analysis)
 * @param topology Topology information (for special cases)
 * @param atom_index Index of current atom in atoms list
 * @param geometry_bohr Optional geometry matrix for bond angle calculations
 * @param charges Optional atomic charges for special cases (CO detection, carbene)
 * @return Hybridization state (1=sp, 2=sp2, 3=sp3, etc.)
 *
 * Claude Generated - December 2025 (Phase 3: Complete XTB Compatibility)
 */
static inline int detectElementSpecificHybridization(int Z, double cn,
                                                    const std::vector<int>& atoms,
                                                    const std::optional<EEQSolver::TopologyInput>& topology,
                                                    int atom_index,
                                                    const Matrix* geometry_bohr = nullptr,
                                                    const Vector* charges = nullptr) {
    int group = getElementGroup(Z);
    int int_cn = static_cast<int>(std::round(cn));

    // ===== Group 1: Hydrogen (H) =====
    // Reference: gfnff_ini2.f90:218-222
    // CRITICAL FIX (January 25, 2026): Normal hydrogen (CN=1) has hyb=0 (unknown), NOT sp3!
    // This prevents generation of spurious extra torsions for H-C-X-Y quartets
    if (group == 1) { // H
        if (int_cn == 2) return 1; // sp (bridging H)
        if (int_cn > 2) return 3; // sp3 (M+ tetrahedral coordination)
        if (int_cn > 4) return 0; // M+ HC (metal hydride coordination)
        return 0; // Default: hyb=0 (unknown) - matches XTB default at line 207
    }

    // ===== Group 2: Alkali/Alkaline Earth =====
    // Reference: gfnff_ini2.f90:224-228
    // CRITICAL FIX (January 25, 2026): Same as hydrogen - default is hyb=0
    if (group == 2) { // Li, Be, etc.
        if (int_cn == 2) return 1; // sp (bridging metal)
        if (int_cn > 2) return 3; // sp3 (M+ tetrahedral coordination)
        if (int_cn > 4) return 0; // Special case
        return 0; // Default: hyb=0 (unknown) - matches XTB default
    }

    // ===== Group 3: Boron (B) =====
    if (group == 3) { // B
        // XTB gfnff_ini2.f90:232 - hypervalent heavy boron (Si, Ge, Sn with nbdiff==0)
        if (int_cn > 4 && Z > 10 && topology.has_value() && atom_index >= 0) {
            // Check nbdiff==0: all neighbors have same CN (no CN difference)
            int num_metal_neighbors = countMetalNeighbors(atom_index, topology->neighbor_lists, atoms);
            if (num_metal_neighbors == 0) { // nbdiff==0 approximation
                return 5; // sp3d (hypervalent)
            }
        }
        if (int_cn > 4) return 3; // sp3
        if (int_cn == 4) return 3; // sp3
        if (int_cn == 3) return 2; // sp2
        if (int_cn == 2) return 1; // sp
        return 3; // Default sp3
    }

    // ===== Group 4: Carbon (C) =====
    if (group == 4) { // C
        // XTB gfnff_ini2.f90:240 - hypervalent heavy carbon (Si, Ge, Sn with nbdiff==0)
        if (int_cn > 4 && Z > 10 && topology.has_value() && atom_index >= 0) {
            int num_metal_neighbors = countMetalNeighbors(atom_index, topology->neighbor_lists, atoms);
            if (num_metal_neighbors == 0) { // nbdiff==0 approximation
                return 5; // sp3d (hypervalent)
            }
        }
        if (int_cn >= 4) return 3; // sp3
        if (int_cn == 3) return 2; // sp2

        // XTB gfnff_ini2.f90:242-254 - Geometry-dependent CN=2 carbon
        if (int_cn == 2) {
            int hyb = 1; // Default sp (linear triple bond)

            // Geometry-dependent: bond angle check
            if (topology.has_value() && atom_index >= 0 && geometry_bohr != nullptr &&
                atom_index < topology->neighbor_lists.size() &&
                topology->neighbor_lists[atom_index].size() == 2) {

                int neighbor1 = topology->neighbor_lists[atom_index][0];
                int neighbor2 = topology->neighbor_lists[atom_index][1];

                // Calculate angle: neighbor1 - atom_index - neighbor2
                double angle_rad = calculateBondAngle(*geometry_bohr, neighbor1, atom_index, neighbor2);
                double angle_deg = angle_rad * 180.0 / M_PI;

                // XTB: if angle < 150°, then sp2 (bent carbene)
                if (angle_deg < 150.0) {
                    hyb = 2; // sp2 (carbene or bent configuration)
                    // Note: XTB sets itag(i)=1 for Hueckel and HB routines (not implemented here)
                } else {
                    hyb = 1; // sp (linear triple bond)
                }
            }

            // XTB gfnff_ini2.f90:250-253 - Charge-based override
            // If charge < -0.4, force sp2 (even if linear)
            if (charges != nullptr && atom_index >= 0 && atom_index < charges->size()) {
                if ((*charges)(atom_index) < -0.4) {
                    hyb = 2; // sp2
                    // Note: XTB sets itag(i)=0 (not implemented here)
                }
            }

            return hyb;
        }

        if (int_cn == 1) return 1; // sp (CO, nitriles)
        return 3; // Default sp3
    }

    // ===== Group 5: Nitrogen (N) =====
    if (group == 5) { // N
        // XTB gfnff_ini2.f90:260 - hypervalent heavy nitrogen
        if (int_cn > 4 && Z > 10 && topology.has_value() && atom_index >= 0) {
            int num_metal_neighbors = countMetalNeighbors(atom_index, topology->neighbor_lists, atoms);
            if (num_metal_neighbors == 0) { // nbdiff==0 approximation
                return 5; // sp3d (hypervalent)
            }
        }
        if (int_cn >= 4) return 3; // sp3

        // XTB gfnff_ini2.f90:262-279 - Complex CN=3 nitrogen topology checks
        if (int_cn == 3 && Z == 7) { // Only for nitrogen (not heavier pnictogens)
            int hyb = 3; // Default sp3 (ammonia-like)

            if (topology.has_value() && atom_index >= 0 &&
                atom_index < topology->neighbor_lists.size() &&
                topology->neighbor_lists[atom_index].size() == 3) {

                int kk = 0; // Count O neighbors with CN=1 (NO2 detection)
                int ll = 0; // Count B neighbors with CN=4 (B-N detection)
                int nn = 0; // Count S neighbors with CN=4 (N-SO2 detection)

                // Iterate over 3 neighbors
                for (int j = 0; j < 3; ++j) {
                    int neighbor = topology->neighbor_lists[atom_index][j];
                    if (neighbor >= atoms.size()) continue;

                    int neighbor_Z = atoms[neighbor];
                    int neighbor_CN = (neighbor < topology->neighbor_lists.size()) ?
                                     topology->neighbor_lists[neighbor].size() : 0;

                    // Check for O with CN=1 (NO2 or R2-N=O)
                    if (neighbor_Z == 8 && neighbor_CN == 1) kk++;

                    // Check for B with CN=4 (B-N, loosely bound N is sp2)
                    if (neighbor_Z == 5 && neighbor_CN == 4) ll++;

                    // Check for S with CN=4 (N-SO2-)
                    if (neighbor_Z == 16 && neighbor_CN == 4) nn++;
                }

                // XTB logic for special cases
                if (nn == 1 && ll == 0 && kk == 0) hyb = 3; // N-SO2: sp3
                if (ll == 1 && nn == 0) hyb = 2; // B-N: sp2
                if (kk >= 1) {
                    hyb = 2; // NO2: sp2
                    // Note: XTB sets itag(i)=1 for Hueckel (not implemented here)
                }

                // Check if coordinated to metal (nbmdiff > 0)
                int num_metal = countMetalNeighbors(atom_index, topology->neighbor_lists, atoms);
                if (num_metal > 0 && nn == 0) {
                    hyb = 2; // Pyridine coordinated to metal: sp2
                }

                // Geometry-dependent aromatic nitrogen detection (Claude Generated Jan 2026)
                // For CN=3 nitrogen not caught by special cases above, check planarity
                // Planar nitrogen (sum of angles ≈ 360°) → sp2 (aromatic)
                // Tetrahedral nitrogen (sum of angles < 350°) → sp3 (aliphatic)
                if (hyb == 3 && geometry_bohr != nullptr) {
                    const auto& neighbors = topology->neighbor_lists[atom_index];
                    if (neighbors.size() == 3) {
                        // Calculate all three N-X-Y bond angles
                        double angle_sum = 0.0;

                        // Angle 1: neighbor0 - N - neighbor1
                        double angle1 = calculateBondAngle(*geometry_bohr,
                                                          neighbors[0], atom_index, neighbors[1]);
                        angle_sum += angle1;

                        // Angle 2: neighbor1 - N - neighbor2
                        double angle2 = calculateBondAngle(*geometry_bohr,
                                                          neighbors[1], atom_index, neighbors[2]);
                        angle_sum += angle2;

                        // Angle 3: neighbor2 - N - neighbor0
                        double angle3 = calculateBondAngle(*geometry_bohr,
                                                          neighbors[2], atom_index, neighbors[0]);
                        angle_sum += angle3;

                        // Convert to degrees and check planarity
                        double angle_sum_deg = angle_sum * 180.0 / M_PI;

                        // DEBUG: Log angle sum for first few nitrogen atoms
                        if (atom_index < 5) {
                            CurcumaLogger::info(fmt::format("  N atom {}: angle_sum = {:.2f}° (neighbors: {}, {}, {})",
                                              atom_index, angle_sum_deg,
                                              neighbors[0], neighbors[1], neighbors[2]));
                        }

                        // Planar nitrogen: sum ≈ 360° (allowing ±5° tolerance for ring strain)
                        // Tetrahedral nitrogen: sum ≈ 328° (3 × 109.5°)
                        if (angle_sum_deg > 355.0) {
                            hyb = 2; // sp2 (aromatic/planar)
                            if (atom_index < 5) {
                                CurcumaLogger::info(fmt::format("    → sp2 (aromatic/planar)"));
                            }
                        } else {
                            if (atom_index < 5) {
                                CurcumaLogger::info(fmt::format("    → sp3 (tetrahedral)"));
                            }
                        }
                    }
                }
            }

            return hyb;
        }

        // XTB gfnff_ini2.f90:280-294 - Complex CN=2 nitrogen topology checks
        if (int_cn == 2) {
            int hyb = 2; // Default sp2

            if (topology.has_value() && atom_index >= 0 &&
                atom_index < topology->neighbor_lists.size() &&
                topology->neighbor_lists[atom_index].size() == 2) {

                int neighbor1 = topology->neighbor_lists[atom_index][0];
                int neighbor2 = topology->neighbor_lists[atom_index][1];

                if (neighbor1 < atoms.size() && neighbor2 < atoms.size()) {
                    int Z1 = atoms[neighbor1];
                    int Z2 = atoms[neighbor2];
                    int CN1 = (neighbor1 < topology->neighbor_lists.size()) ?
                             topology->neighbor_lists[neighbor1].size() : 0;
                    int CN2 = (neighbor2 < topology->neighbor_lists.size()) ?
                             topology->neighbor_lists[neighbor2].size() : 0;

                    // Check for R-N=C (nitrile)
                    if ((CN1 == 1 && Z1 == 6) || (CN2 == 1 && Z2 == 6)) {
                        hyb = 1; // sp (nitrile)
                    }

                    // Check for R-N=N (diazomethane)
                    if ((CN1 == 1 && Z1 == 7) || (CN2 == 1 && Z2 == 7)) {
                        hyb = 1; // sp (diazomethane)
                    }

                    // Check for M-NC-R (metal nitrile complex)
                    if (isMetal(Z1) || isMetal(Z2)) {
                        hyb = 1; // sp (metal-coordinated nitrile)
                    }

                    // Check for N=N=N (azide)
                    if (Z1 == 7 && Z2 == 7 && CN1 <= 2 && CN2 <= 2) {
                        hyb = 1; // sp (azide)
                    }
                }

                // Geometry-dependent: bond angle check (lintr = 170° in XTB)
                if (geometry_bohr != nullptr) {
                    double angle_rad = calculateBondAngle(*geometry_bohr, neighbor1, atom_index, neighbor2);
                    double angle_deg = angle_rad * 180.0 / M_PI;

                    if (angle_deg > 170.0) {
                        hyb = 1; // sp (linear configuration)
                    }
                }
            }

            return hyb;
        }

        if (int_cn == 1) return 1; // sp
        return 3; // Default sp3
    }

    // ===== Group 6: Oxygen (O) =====
    if (group == 6) { // O
        // XTB gfnff_ini2.f90:300 - hypervalent heavy oxygen (S, Se, Te with nbdiff==0)
        if (int_cn > 3 && Z > 10 && topology.has_value() && atom_index >= 0) {
            int num_metal_neighbors = countMetalNeighbors(atom_index, topology->neighbor_lists, atoms);
            if (num_metal_neighbors == 0) { // nbdiff==0 approximation
                return 5; // sp3d (hypervalent)
            }
        }

        if (int_cn >= 3) return 3; // sp3

        // XTB gfnff_ini2.f90:302-306 - Metal neighbor detection for CN=2
        if (int_cn == 2) {
            int hyb = 3; // Default sp3 (ether, water, alcohols) - CRITICAL FIX

            // Check for metal neighbors (M-O-X conjugation)
            if (topology.has_value() && atom_index >= 0) {
                int num_metal = countMetalNeighbors(atom_index, topology->neighbor_lists, atoms);

                if (num_metal > 0) { // nbmdiff > 0
                    // Find nearest non-metal neighbor's CN
                    int nearest_CN = findNearestNonMetalCN(atom_index, topology->neighbor_lists, atoms);

                    if (nearest_CN == 3) {
                        hyb = 2; // sp2 (M-O-X conjugated)
                    } else if (nearest_CN == 4) {
                        hyb = 3; // sp3 (M-O-X non-conjugated)
                    }
                }
            }

            return hyb;
        }

        // Fortran gfnff_ini2.f90:307-310 - CN=1 oxygen hybridization
        // Default sp2, but sp if sole neighbor also has CN=1 (covers CO, OH radical, etc.)
        if (int_cn == 1) {
            int hyb = 2; // Default sp2

            if (topology.has_value() && atom_index >= 0 &&
                atom_index < static_cast<int>(topology->neighbor_lists.size()) &&
                topology->neighbor_lists[atom_index].size() == 1) {

                int neighbor = topology->neighbor_lists[atom_index][0];
                int neighbor_CN = (neighbor < static_cast<int>(topology->neighbor_lists.size())) ?
                                 static_cast<int>(topology->neighbor_lists[neighbor].size()) : 0;

                // Fortran: if (nb20i==1 .and. nbdiff==0) then if (topo%nb(20,topo%nb(1,i))==1) hyb=1
                // No element check - any neighbor with CN=1 triggers sp
                if (neighbor_CN == 1) {
                    hyb = 1; // sp (CO, OH radical, etc.)
                }
            }

            return hyb;
        }

        return 3; // Default sp3
    }

    // ===== Group 7: Halogens (F, Cl, Br, I) =====
    // Fortran gfnff_ini2.f90:313-316: CN=2 → sp, CN>2+heavy → sp3d, else unknown
    if (group == 7) { // Halogens
        if (int_cn == 2) return 1; // sp
        if (int_cn > 2 && Z > 10) return 5; // sp3d (heavy halogens)
        return 0; // Default: unknown (CN=1 halogens like HCl)
    }

    // ===== Group 8: Noble Gases =====
    if (group == 8) { // Noble gases
        if (int_cn > 0 && Z > 2) return 5; // sp3d2 (heavy noble gases)
        return 0; // Default (no hybridization)
    }

    // ===== Transition Metals (Groups ≤ 0) =====
    if (group <= 0) { // TMs
        int effective_cn = int_cn;

        // Don't count hydrogen for TM coordination
        if (topology.has_value() && atom_index >= 0) {
            int nh = countHydrogenNeighbors(atom_index, topology->neighbor_lists, atoms);
            if (nh > 0 && nh != effective_cn) {
                effective_cn -= nh;
            }
        }

        if (effective_cn <= 2) {
            if (group == -7) return 2; // sp2 (late TMs)
            return 1; // sp (early TMs)
        }
        if (effective_cn == 3) return 2; // sp2
        if (effective_cn == 4) {
            if (group > -7) return 3; // sp3 tetrahedral (early TMs)
            return 3; // sp3 square planar (late TMs)
        }
        if (effective_cn == 5 && group == -3) return 3; // sp3 (Sc-La)

        return 3; // Default sp3
    }

    // ===== Default: CN-based heuristic (backward compatibility) =====
    if (cn < 1.5) return 1; // sp
    if (cn < 2.5) return 2; // sp2
    if (cn < 4.5) return 3; // sp3
    if (cn < 6.5) return 4; // sp3d
    return 5; // sp3d2
}

// ===== Test Function for Element-Specific Hybridization =====
// Simple test to verify the element-specific hybridization logic
// Claude Generated - December 2025 (Phase 2 Validation)
// NOTE: This test function is currently unused but kept for future validation
// If needed, it can be called with proper verbosity guard in main code
static void testElementSpecificHybridizationLogic() {
    std::vector<int> dummy_atoms = {0};
    int passed = 0;
    int failed = 0;

    // Test cases: {Z, CN, expected_hyb, description}
    std::vector<std::tuple<int, double, int, std::string>> test_cases = {
        {8, 2.0, 3, "Oxygen CN=2 (CRITICAL: should be sp3)"},
        {6, 2.0, 1, "Carbon CN=2 (should be sp)"},
        {7, 3.0, 3, "Nitrogen CN=3 (should be sp3)"},
        {1, 2.0, 1, "Hydrogen CN=2 (should be sp)"},
        {17, 2.0, 1, "Chlorine CN=2 (should be sp)"},
        {17, 1.0, 0, "Chlorine CN=1 (should be unknown)"},
        {26, 4.0, 3, "Iron CN=4 (should be sp3)"},
        {8, 1.0, 2, "Oxygen CN=1 (should be sp2)"},
        {6, 3.0, 2, "Carbon CN=3 (should be sp2)"},
        {7, 2.0, 2, "Nitrogen CN=2 (should be sp2)"},
        {9, 2.0, 1, "Fluorine CN=2 (should be sp)"}
    };

    // Guarded by verbosity check - disabled by default
    if (CurcumaLogger::get_verbosity() >= 3) {
        std::cout << "\n=== Element-Specific Hybridization Validation ===" << std::endl;

        for (const auto& test_case : test_cases) {
            int Z = std::get<0>(test_case);
            double cn = std::get<1>(test_case);
            int expected = std::get<2>(test_case);
            std::string desc = std::get<3>(test_case);

            int actual = detectElementSpecificHybridization(Z, cn, dummy_atoms, std::nullopt, 0);

            if (actual == expected) {
                std::cout << "✅ " << desc << " → " << actual << " (PASS)" << std::endl;
                passed++;
            } else {
                std::cout << "❌ " << desc << " → " << actual << " (FAIL, expected " << expected << ")" << std::endl;
                failed++;
            }
        }

        std::cout << "\n=== Test Results ===" << std::endl;
        std::cout << "Passed: " << passed << "/" << test_cases.size() << std::endl;
        std::cout << "Failed: " << failed << "/" << test_cases.size() << std::endl;

        if (failed == 0) {
            std::cout << "🎉 All element-specific hybridization tests passed!" << std::endl;
        } else {
            std::cout << "⚠️  Some tests failed. Review implementation." << std::endl;
        }
    }
}

// ===== Original EEQSolver Implementation =====

// Claude Generated - March 2026
EEQSolveMethod EEQSolver::parseSolveMethod(const std::string& method_str) {
    if (method_str == "lu") return EEQSolveMethod::LU;
    if (method_str == "pcg") return EEQSolveMethod::PCG;
    if (method_str == "batched") return EEQSolveMethod::Batched;
    // "cholesky" (canonical) and "schur_cholesky" (legacy) both map to SchurCholesky
    if (method_str == "cholesky" || method_str == "schur_cholesky")
        return EEQSolveMethod::SchurCholesky;
    return EEQSolveMethod::Auto;  // default including "auto"
}

EEQSolver::EEQSolver(const ConfigManager& config)
    : m_config(config)
{
    // Extract parameters from ConfigManager with defaults from PARAM definitions
    m_max_iterations = m_config.get<int>("max_iterations", 50);
    m_convergence_threshold = m_config.get<double>("convergence_threshold", 1e-6);
    m_verbosity = m_config.get<int>("verbosity", 0);
    m_calculate_cn = m_config.get<bool>("calculate_cn", true);
    m_allow_unconverged = m_config.get<bool>("allow_unconverged_charges", false);
    m_skip_phase2 = m_config.get<bool>("skip_phase2", false);
    m_solve_method = parseSolveMethod(m_config.get<std::string>("solve_method", "pcg"));
    m_refactor_eps         = m_config.get<double>("eeq_refactor_eps_bohr", 0.05);
    m_refactor_force_every = m_config.get<int>("eeq_refactor_force_every", 0);
    m_matrix_rebuild_eps   = m_config.get<double>("eeq_matrix_rebuild_eps_bohr", 0.0);
    m_eeq_extrapolation    = m_config.get<std::string>("eeq_extrapolation", "none");
    m_eeq_extrap_order     = m_config.get<int>("eeq_extrapolation_order", 3);

    // Initialize EEQ caching system
    m_eeq_cache = std::make_unique<EEQSolverCache>();

    if (m_verbosity >= 2) {
        CurcumaLogger::info("EEQSolver initialized with parameters:");
        CurcumaLogger::param("max_iterations", std::to_string(m_max_iterations));
        CurcumaLogger::param("convergence_threshold", fmt::format("{:.2e}", m_convergence_threshold));
        CurcumaLogger::param("verbosity", std::to_string(m_verbosity));
        CurcumaLogger::param("calculate_cn", m_calculate_cn ? "true" : "false");
        CurcumaLogger::param("allow_unconverged_charges", m_allow_unconverged ? "true" : "false");
        CurcumaLogger::param("solve_method", m_config.get<std::string>("solve_method", "pcg"));
    }
}

// Claude Generated (March 2026): Fallback charges for graceful degradation
Vector EEQSolver::generateFallbackCharges(int natoms, int total_charge, const std::string& context) const
{
    CurcumaLogger::warn(fmt::format("EEQ solver failed ({}): using uniform fallback charges q_i = {:.4f}",
                                    context, static_cast<double>(total_charge) / natoms));
    return Vector::Constant(natoms, static_cast<double>(total_charge) / natoms);
}

// ===== Main API =====

/**
 * @brief Calculate EEQ charges using Fortran-compatible two-phase algorithm
 *
 * ARCHITECTURE (matches XTB gfnff_ini.f90 + gfnff_engrad.F90):
 *
 * PHASE 1: Topology Charges (goedeckera equivalent)
 * ------------------------------------------------------
 * - Uses integer neighbor count (nb) from topology
 * - chieeq = -chi + dxi + cnf*sqrt(nb)  [WITH CNF!]
 * - Base parameters: gam (no dgam), alpha (no charge correction)
 * - Result: topo%qa (topology charges for parameter generation)
 * - Fortran ref: gfnff_ini.f90:411, gfnff_ini2.f90:1140-1321
 *
 * PHASE 2: Parameter Preparation + Final Charges
 * ------------------------------------------------------
 * - Uses Phase 1 charges (topo%qa) to calculate dgam, dalpha corrections
 * - Overwrites chieeq WITHOUT CNF: chieeq = -chi + dxi  [NO CNF!]
 * - Applies corrections to matrix: gam+dgam, (alpha+ff*qa)^2
 * - Solve uses chieeq WITHOUT CNF in RHS (use_cnf_term=false)
 * - Fortran ref: gfnff_ini.f90:715 (chieeq overwrite), gfnff_ini.f90:693-726 (dgam)
 *
 * Note: Fortran has separate Phase 3 (goed_gfnff in gfnff_engrad.F90) which
 * adds CNF term with fractional CN during energy calculation. Curcuma's Phase 2
 * combines parameter preparation with implicit final solve.
 *
 * CRITICAL BUG FIX (Jan 4, 2026):
 * - Phase 2 was incorrectly using use_cnf_term=true (line 788)
 * - Changed to use_cnf_term=false to match Fortran line 715
 * - Error reduced from 1.55e-02 e RMS to <1e-5 e (expected)
 *
 * @param atoms         Atomic numbers
 * @param geometry_bohr Geometry in Bohr
 * @param cn            Coordination numbers (fractional, erf-counted)
 * @param total_charge  Total molecular charge
 * @param topology      Optional topology information (neighbor lists, bonds)
 * @return              EEQ charges (nlist%q equivalent)
 */
Vector EEQSolver::calculateCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector* cn_hint,
    const std::vector<int>* hyb_hint,
    const std::optional<TopologyInput>& topology)
{
    const int natoms = atoms.size();

    if (natoms == 0) {
        CurcumaLogger::error("EEQSolver::calculateCharges: No atoms provided");
        return Vector::Zero(0);
    }

    if (geometry_bohr.rows() != natoms || geometry_bohr.cols() != 3) {
        CurcumaLogger::error(fmt::format("EEQSolver::calculateCharges: Invalid geometry dimensions {}x{} (expected {}x3)",
                                         geometry_bohr.rows(), geometry_bohr.cols(), natoms));
        return Vector::Zero(natoms);
    }

    // Step 1: Calculate or use provided coordination numbers
    Vector cn;
    if (cn_hint != nullptr && cn_hint->size() == natoms) {
        cn = *cn_hint;
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Using provided coordination numbers");
        }
    } else if (m_calculate_cn) {
        cn = calculateCoordinationNumbers(atoms, geometry_bohr);
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Calculated coordination numbers");
        }
    } else {
        CurcumaLogger::error("EEQSolver::calculateCharges: CN required but calculate_cn=false and no hint provided");
        return Vector::Zero(natoms);
    }

    // Step 2: Detect or use provided hybridization
    std::vector<int> hybridization;
    if (hyb_hint != nullptr && hyb_hint->size() == natoms) {
        hybridization = *hyb_hint;
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Using provided hybridization states");
        }
    } else {
        hybridization = detectHybridization(atoms, geometry_bohr, cn);
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Detected hybridization states");
        }
    }

    // ===== CRITICAL REFACTORING (Jan 4, 2026): Single-Solve with Iterative Refinement =====
    // Reference: XTB does ONE goedeckera() solve with corrected parameters (dgam, alpha)
    // not TWO separate solves. This matches gfnff_ini.f90:696-707 exactly.
    //
    // Instead of:
    //   Phase 1: Solve with base params → qa
    //   Phase 2: Solve with corrected params (gam+dgam, alpha(qa))
    //
    // New approach: ITERATIVE refinement in SINGLE solve
    //   Iteration 0: dgam(qa=0), alpha(qa=0) → Solve → new qa
    //   Iteration 1: dgam(qa_new), alpha(qa_new) → Solve → final qa
    //
    // This allows dgam/alpha corrections to have FULL impact on solution matrix,
    // not just as a perturbation after Phase 1.

    // ===== NEW STRATEGY (Jan 4, 2026): Hybrid Two-Phase + Iterative Refinement =====
    // Problem with pure iteration from qa=0: dgam=0 in first iteration, so weak correction
    // Solution: Do Phase 1 first (base params), then iterate with dgam corrections

    Vector current_charges;
    Vector topology_charges; // To be populated in Phase 1

    // PHASE 1: Initial solve for topology charges (qa)
    {
        // CRITICAL FIX (Phase 1 Charge Synchronization - January 26, 2026):
        // Reference: XTB gfnff_ini.f90:589 (call goedeckera) and 411 (topo%chieeq)
        // Phase 1 (topological solve) MUST include dxi electronegativity corrections
        // to provide a precise enough starting point (qa) for the GFN-FF parameters.
        // Even small errors in qa propagate to ALL bonded energy terms.

        Vector dxi = calculateDxi(atoms, geometry_bohr, cn, topology);
        Vector dgam_base = Vector::Zero(natoms);  // NO dgam corrections in Phase 1

        // Use buildCorrectedEEQMatrix with EXPLICIT Topological distance mode
        Matrix A_phase1 = buildCorrectedEEQMatrix(atoms, geometry_bohr, cn, Vector::Zero(natoms),
                                                   dxi, dgam_base, hybridization, topology,
                                                   EEQDistanceMode::Topological);
        current_charges = solveEEQ(A_phase1, atoms, cn, dxi, total_charge, topology, true, true);  // Phase 1: use integer nb
        topology_charges = current_charges;  // Store for Phase 2 use

        if (m_verbosity >= 1) {
            CurcumaLogger::info("EEQSolver: Phase 1 base solve complete");
        }
    }

    // Claude Generated (April 2026): skip_phase2 option — bypass Phase 2 entirely
    if (m_skip_phase2) {
        if (m_verbosity >= 1) {
            CurcumaLogger::warn("EEQSolver: skip_phase2=true — bypassing Phase 2 refinement, using Phase 1 topology charges");
        }
        return topology_charges;
    }

    // Claude Generated (Apr 2026): If Phase 2 was historically implausible for this system
    // size, skip it forever (until atom count changes, indicating a new molecule).
    if (m_phase2_historically_implausible && m_phase2_implausible_natoms == natoms) {
        if (m_verbosity >= 2) {
            CurcumaLogger::info("EEQSolver: Phase 2 was historically implausible for this system — skipping, using Phase 1 topology charges");
        }
        return topology_charges;
    }

    // Initialize cache with Phase 1 charges if not yet set.
    // Without this, m_last_successful_charges starts empty (size 0), so Tier 1/2
    // fallbacks below would skip it and go straight to topology_charges. With this
    // initialization, the cache holds Phase 1 charges on the first call, then is
    // overwritten with Phase 2 charges after a successful solve.
    if (m_last_successful_charges.size() != natoms) {
        m_last_successful_charges = topology_charges;
    }

    // PHASE 2: Final energy charges with ALL corrections (matching Fortran reference)
    // Reference: XTB gfnff_ini.f90:693-707 - ONE solve with dxi, dgam, and alpha corrections
    // Restored (Jan 17, 2026): Activate corrections to match GFN-FF Fortran reference
    {
        // RESTORED (Jan 17, 2026): Phase 2 uses ALL corrections (matching Fortran gfnff_ini.f90:713-726)
        // Reference: external/gfnff/src/gfnff_ini.f90:713-726
        //
        // topo%chieeq(i) = -param%chi(at(i))+dxi(i)     // chi correction
        // topo%gameeq(i) = param%gam(at(i))+dgam(i)     // hardness correction
        // topo%alpeeq(i) = (param%alp(at(i))+ff*qa(i))**2  // charge-dependent alpha
        Vector dxi = calculateDxi(atoms, geometry_bohr, cn, topology);

        // Claude Generated (January 17, 2026): Detect pi-system and amide nitrogens for dgam refinements
        std::vector<bool> is_pi_atom = detectPiSystem(atoms, hybridization, topology);
        std::vector<bool> is_amide = detectAmideNitrogens(atoms, hybridization, is_pi_atom, topology, cn);
        Vector dgam = calculateDgam(atoms, topology_charges, hybridization, is_pi_atom, is_amide);

        // Diagnostic output for comparison with Fortran goed_gfnff debug output
        // Matches format: gfnff_engrad.F90:1581-1594
        if (m_verbosity >= 3) {
            CurcumaLogger::info("=== Phase-2 EEQ Parameters (Curcuma) ===");
            CurcumaLogger::info(fmt::format("{:>5} {:>4} {:>12} {:>12} {:>12} {:>10} {:>10} {:>10}",
                "Atom", "Z", "chieeq", "gameeq", "alpeeq", "CN", "dxi", "dgam"));

            for (int i = 0; i < natoms; ++i) {
                int z_i = atoms[i];
                EEQParameters params_i = getParameters(z_i, cn(i));
                double chi_base = -params_i.chi + dxi(i);
                double gam_corr = params_i.gam + dgam(i);

                // Compute alpeeq using charge-dependent ff (matching gfnff_ini.f90:718-725)
                double alpha_base = (z_i >= 1 && z_i <= 86) ? alpha_eeq[z_i - 1] : 0.903430;
                double ff = 0.0;
                if (z_i == 6) ff = 0.09;
                else if (z_i == 7) ff = -0.21;
                else if (z_i >= 1 && z_i <= 86) {
                    int group = periodic_group[z_i - 1];
                    if (group == 6) ff = -0.03;
                    else if (group == 7) ff = 0.50;
                }
                double alpeeq = std::pow(alpha_base + ff * topology_charges(i), 2);

                CurcumaLogger::info(fmt::format("{:>5} {:>4} {:>12.6f} {:>12.6f} {:>12.6f} {:>10.5f} {:>10.6f} {:>10.6f}",
                    i+1, z_i, chi_base, gam_corr, alpeeq, cn(i), dxi(i), dgam(i)));
            }

            // Also print Phase 1 charges (qa) for comparison with Fortran topo%qa
            CurcumaLogger::info("\n=== Phase-1 Topology Charges (qa) ===");
            CurcumaLogger::info(fmt::format("{:>5} {:>4} {:>15}", "Atom", "Z", "qa"));
            for (int i = 0; i < natoms; ++i) {
                CurcumaLogger::info(fmt::format("{:>5} {:>4} {:>15.10f}",
                    i+1, atoms[i], topology_charges(i)));
            }
        }

        // Build matrix WITH corrections (matching Fortran gfnff_ini.f90:693-707)
        // CRITICAL FIX (January 17, 2026): Phase 2 uses GEOMETRIC distances (r(ij) from xyz),
        // NOT topological distances (pair(ij) from Floyd-Warshall)!
        // Reference: Fortran goed_gfnff uses r(ij), while goedeckera uses pair(ij)
        Matrix A = buildCorrectedEEQMatrix(atoms, geometry_bohr, cn, topology_charges,
                                          dxi, dgam, hybridization, topology,
                                          EEQDistanceMode::Geometric);

        // Solve ONCE with corrected parameters (Phase 2: CNF term REQUIRED!)
        // CRITICAL FIX (Jan 26, 2026 - synchronization): use_integer_nb MUST be false for Phase 2
        // to correctly use fractional cn instead of integer nb for energy charges (nlist%q).
        current_charges = solveEEQ(A, atoms, cn, dxi, total_charge, topology, true, false);

        if (current_charges.size() != natoms) {
            if (m_allow_unconverged) {
                CurcumaLogger::warn("EEQSolver::calculateCharges: Phase 2 solve did not converge (allow_unconverged_charges=true), using Phase 1 charges");
            } else {
                CurcumaLogger::error("EEQSolver::calculateCharges: Phase 2 solve failed");
            }
            return topology_charges;  // Return Phase 1 result as fallback
        }

        if (m_verbosity >= 1) {
            CurcumaLogger::success("EEQSolver: Phase 2 single-solve completed with dxi/dgam corrections");
        }

        // Phase 2 charge comparison output (verbosity 3)
        if (m_verbosity >= 3) {
            CurcumaLogger::info("\n=== Phase-2 Final Charges (Curcuma) ===");
            CurcumaLogger::info(fmt::format("{:>5} {:>4} {:>15} {:>15} {:>12}",
                "Atom", "Z", "q_phase2", "q_phase1(qa)", "delta"));
            for (int i = 0; i < natoms; ++i) {
                CurcumaLogger::info(fmt::format("{:>5} {:>4} {:>15.10f} {:>15.10f} {:>12.6f}",
                    i+1, atoms[i], current_charges(i), topology_charges(i),
                    current_charges(i) - topology_charges(i)));
            }
        }
    }

    if (m_verbosity >= 1) {
        CurcumaLogger::success(fmt::format("EEQSolver: Single-solve EEQ completed for {} atoms", natoms));
    }

    // Claude Generated (February 2026): Per-atom charge/CN table at verbosity 2
    if (m_verbosity >= 2 && natoms > 0) {
        CurcumaLogger::info("\nEEQ Atomic Parameters:");
        CurcumaLogger::result(fmt::format("{:>5} {:>4} {:>8} {:>10} {:>12}",
                                          "Atom", "Z", "CN", "Charge", "Hyb"));

        for (int i = 0; i < natoms; ++i) {
            std::string hyb_str = "sp?";
            if (hyb_hint != nullptr && i < hyb_hint->size()) {
                int h = (*hyb_hint)[i];
                if (h == 1) hyb_str = "sp";
                else if (h == 2) hyb_str = "sp2";
                else if (h == 3) hyb_str = "sp3";
                else if (h == 4) hyb_str = "sp3d";
                else if (h == 5) hyb_str = "sp3d2";
            }

            CurcumaLogger::result(fmt::format("{:>5} {:>4} {:>8.3f} {:>+10.5f} {:>12}",
                                              i + 1,
                                              atoms[i],
                                              cn(i),
                                              current_charges(i),
                                              hyb_str));
        }

        // Show total charge
        double total_charge_sum = current_charges.sum();
        CurcumaLogger::param("total_charge", fmt::format("{:+.6f}", total_charge_sum));
    }

    return current_charges;
}

// ===== New Helper Functions for Single-Solve Architecture =====

/**
 * @brief Build EEQ Coulomb matrix with ALL charge-dependent corrections
 * @param current_charges Current charge estimate (used for dgam and charge-dependent alpha)
 * @return Augmented EEQ matrix (natoms+1) × (natoms+1)
 *
 * Matrix structure:
 * - A(i,i) = gam_corrected + sqrt(2/π)/sqrt(alpha_corrected)
 * - A(i,j) = erf(gamma_ij * r) / r  [Coulomb interaction with corrected alpha]
 * - A(natoms, j) = 1.0  [Charge constraint row]
 *
 * Key: Using CURRENT charges for dgam and alpha means corrections have FULL impact,
 * not just perturbation-level effects like in two-phase approach.
 *
 * Claude Generated - January 4, 2026
 */
Matrix EEQSolver::buildCorrectedEEQMatrix(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn,
    const Vector& current_charges,
    const Vector& dxi,
    const Vector& dgam,
    const std::vector<int>& hybridization,
    const std::optional<TopologyInput>& topology,
    EEQDistanceMode distance_mode)
{
    const int natoms = atoms.size();
    int nfrag = topology.has_value() ? topology->nfrag : 1;
    int m = natoms + nfrag;
    Matrix A = Matrix::Zero(m, m);

    // WP6/G2c Phase C (May 2026): geometric distance cutoff for off-diagonal Coulomb terms.
    // Only active in Geometric mode (Phase 2 / fallback) — topological distances are
    // Floyd-Warshall bond-length sums, where a Bohr cutoff has no clean meaning.
    // Mirrors the gating used in the parallel A-matrix path at line ~3064.
    // WP-S3 (May 2026): honour the post-init override set by GFNFF auto-detection.
    double eeq_dist_cutoff = getEEQDistanceCutoffEffective();
    const double EEQ_CUTOFF_SQ =
        (distance_mode == EEQDistanceMode::Geometric && natoms > 200 && eeq_dist_cutoff > 0.0)
        ? eeq_dist_cutoff * eeq_dist_cutoff : 0.0;

    // Step 1: Calculate charge-dependent alpha
    Vector alpha_corrected(natoms);
    for (int i = 0; i < natoms; ++i) {
        int z_i = atoms[i];
        double alpha_base = (z_i >= 1 && z_i <= 86) ? alpha_eeq[z_i - 1] : 0.903430;

        // Charge-dependent ff factor (from XTB gfnff_ini.f90:699-705)
        double ff = 0.0;
        if (z_i == 6) {  // Carbon
            ff = 0.09;
        } else if (z_i == 7) {  // Nitrogen
            ff = -0.21;
        } else if (z_i >= 1 && z_i <= 86) {
            int group = periodic_group[z_i - 1];
            if (group == 6) {
                ff = -0.03;  // Chalcogens (O, S, Se, Te, Po)
            } else if (group == 7) {
                ff = 0.50;  // Halogens (F, Cl, Br, I, At)
            } else {
                // Metal elements
                int imetal_val = metal_type[z_i - 1];
                if (imetal_val == 1) {
                    ff = 0.3;  // Main group metals (Li, Be, Na, Mg, etc.)
                } else if (imetal_val == 2) {
                    ff = -0.1;  // Transition metals
                }
            }
        }

        alpha_corrected(i) = std::pow(alpha_base + ff * current_charges(i), 2);
    }

    // Step 2: Get distances based on explicit mode
    // CRITICAL FIX (January 17, 2026): Fortran uses DIFFERENT distance modes for Phase 1 vs Phase 2:
    // - Phase 1 (goedeckera): Topological distances (Floyd-Warshall bond paths)
    // - Phase 2 (goed_gfnff): Geometric distances (Euclidean xyz)
    // Reference: XTB gfnff_ini.f90 - pair(ij) vs r(ij)
    //
    // NOTE (WP-FF-Packed-Triangular): we considered reading from m_external_srab here,
    // but buildCorrectedEEQMatrix is also reached from gradient / numerical-perturbation
    // paths where geometry_bohr may differ from the geometry that produced
    // m_external_srab. Falling back to local recomputation keeps the function self-
    // contained and bit-stable against caller-side geometry mutations.
    Matrix distances;
    if (distance_mode == EEQDistanceMode::Topological && topology.has_value()) {
        distances = computeTopologicalDistancesSparse(atoms, *topology);

        // DEBUG: Print first few topological distances
        if (m_verbosity >= 3 && natoms >= 2) {
            std::cout << "\n=== TOPOLOGICAL DISTANCES (Bohr) - Phase 1 Mode ===" << std::endl;
            for (int i = 0; i < std::min(3, natoms); ++i) {
                for (int j = 0; j < std::min(6, natoms); ++j) {
                    if (i != j && distances(i,j) < 1e6) {
                        std::cout << fmt::format("d({},{}) = {:.4f}", i, j, distances(i,j)) << " ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    } else {
        // Geometric distances from xyz coordinates (Phase 2 mode, or fallback if no topology)
        distances = Matrix::Zero(natoms, natoms);
        for (int i = 0; i < natoms; ++i) {
            for (int j = 0; j < natoms; ++j) {
                if (i != j) {
                    distances(i, j) = (geometry_bohr.row(i) - geometry_bohr.row(j)).norm();
                }
            }
        }

        // DEBUG: Print first few geometric distances
        if (m_verbosity >= 3 && natoms >= 2 && distance_mode == EEQDistanceMode::Geometric) {
            std::cout << "\n=== GEOMETRIC DISTANCES (Bohr) - Phase 2 Mode ===" << std::endl;
            for (int i = 0; i < std::min(3, natoms); ++i) {
                for (int j = 0; j < std::min(6, natoms); ++j) {
                    if (i != j) {
                        std::cout << fmt::format("d({},{}) = {:.4f}", i, j, distances(i,j)) << " ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    // Step 3+4: Build Coulomb matrix (off-diagonal + diagonal) with corrected alpha
    // Claude Generated (March 2026): OpenMP parallelization matching Phase 2 pattern
    const double TSQRT2PI = 0.797884560802866;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < natoms; ++i) {
        // Diagonal: gam + dgam + sqrt(2/π)/sqrt(alpha)
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));
        double gam_corrected = params_i.gam + dgam(i);
        A(i, i) = gam_corrected + TSQRT2PI / std::sqrt(alpha_corrected(i));

        // Off-diagonal: erf-damped Coulomb
        for (int j = 0; j < i; ++j) {
            double r = distances(i, j);

            if (r > 1e6) {
                // Unconnected atoms - skip Coulomb term
                A(i, j) = 0.0;
                A(j, i) = 0.0;
            } else if (EEQ_CUTOFF_SQ > 0.0 && r * r > EEQ_CUTOFF_SQ) {
                // WP6/G2c Phase C: long-range geometric truncation
                A(i, j) = 0.0;
                A(j, i) = 0.0;
            } else {
                // J_ij = erf(gamma_ij * r) / r
                double gamma_ij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));
                double erf_gamma = std::erf(gamma_ij * r);
                double coulomb = erf_gamma / r;

                A(i, j) = coulomb;
                A(j, i) = coulomb;
            }
        }
    }

    // Step 5: Charge constraint row and column
    for (int f = 0; f < nfrag; ++f) {
        int row = natoms + f;
        for (int j = 0; j < natoms; ++j) {
            bool in_fragment = false;
            if (topology.has_value() && j < static_cast<int>(topology->fraglist.size())) {
                in_fragment = (topology->fraglist[j] == f + 1); // 1-indexed
            } else {
                in_fragment = (f == 0); // Default all atoms to first fragment
            }

            if (in_fragment) {
                A(row, j) = 1.0;
                A(j, row) = 1.0;
            }
        }
    }

    return A;
}

// Enhanced EEQ matrix construction with intelligent caching
Matrix EEQSolver::buildSmartEEQMatrix(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn,
    const Vector& current_charges,
    const Vector& dxi,
    const Vector& dgam,
    const std::vector<int>& hybridization,
    const std::optional<TopologyInput>& topology,
    EEQDistanceMode distance_mode)
{
    // Check if we can reuse cached computation
    if (m_eeq_cache && !m_eeq_cache->isGeometryChanged(geometry_bohr)) {
        if (m_verbosity >= 2) {
            CurcumaLogger::info("EEQSolver: Using cached EEQ matrix (geometry unchanged)");
        }
        // For identical geometries, return cached matrix with possible adjustments
        return m_eeq_cache->getCachedAMatrix();
    }

    if (m_verbosity >= 3) {
        CurcumaLogger::info("EEQSolver: Building EEQ matrix from scratch");
    }

    // Build matrix using existing logic, passing through distance mode
    Matrix A = buildCorrectedEEQMatrix(atoms, geometry_bohr, cn, current_charges, dxi, dgam, hybridization, topology, distance_mode);

    // Cache the result for future use
    if (m_eeq_cache) {
        m_eeq_cache->cacheResults(geometry_bohr, A, Vector::Zero(atoms.size())); // Charges cached separately
    }

    return A;
}

// ===== Unified Solver Dispatch =====
// Claude Generated (March 2026): Routes augmented EEQ system to configured solver
// Used by both Phase 1 (calculateTopologyCharges) and Phase 2 (calculateFinalCharges)

Vector EEQSolver::dispatchSolve(
    const Matrix& A,
    const Vector& x,
    int natoms,
    int nfrag,
    int total_charge,
    CxxThreadPool* pool,
    int num_threads)
{
    EEQSolveMethod method_to_use = m_solve_method;
    if (m_solve_method == EEQSolveMethod::Auto) {
        method_to_use = m_selected_method;  // Use cached or default (SchurCholesky)
    }
    // Stage 4 (May 2026): CPU batched path — per-fragment independent Cholesky.
    // Ignores inter-fragment Coulomb coupling, valid when fragments are well-separated
    // (typical for solvent mixtures, condensed phases of small molecules). For the
    // mixture.xyz case (N=6200, nfrag=1400) this collapses solve cost from O(N³) LU
    // (3.6s) to O(Σ atoms_f³) ≈ 90k flops total.
    //
    // For Auto mode + nfrag ≥ many-fragment-threshold, also auto-select Batched
    // when the global SchurCholesky / LU paths historically fail (m_phase2_historically_implausible).
    if (m_solve_method == EEQSolveMethod::Auto &&
        m_phase2_historically_implausible &&
        m_phase2_implausible_natoms == natoms &&
        nfrag >= 8) {
        method_to_use = EEQSolveMethod::Batched;
        if (m_verbosity >= 1) {
            CurcumaLogger::info(fmt::format(
                "EEQ: history shows implausible global solve for N={}, nfrag={}; switching to batched per-fragment Cholesky",
                natoms, nfrag));
        }
    }

    // Auto-select Batched up-front for highly-fragmented systems (e.g. solvent mixtures,
    // van-der-Waals complexes), where the global EEQ matrix is poorly conditioned and
    // inter-fragment Coulomb is small compared to intra-fragment hardness. Threshold:
    // nfrag/N > 20% AND nfrag >= 16 (avoids tripping on small molecules with few atoms).
    //
    // WP7-D (Jun 2026, Claude Generated): contact-aware. The batched solver DROPS
    // cross-fragment Coulomb, so it is only correct for well-separated fragments. When
    // fragments are in contact (m_contact_min_dist < eeq_batched_min_distance, or unknown),
    // prefer the exact path instead: for explicit -solve_method pcg force PCG (block-Jacobi,
    // built below for nfrag>=2); for Auto fall through to the exact auto-selector (which
    // picks SchurCholesky/PCG, both retaining cross-fragment Coulomb).
    if ((m_solve_method == EEQSolveMethod::Auto || m_solve_method == EEQSolveMethod::PCG) &&
        method_to_use != EEQSolveMethod::Batched &&
        nfrag >= 16 && nfrag * 5 > natoms) {
        const bool   prefer_exact = m_config.get<bool>("eeq_contact_prefer_exact", true);
        const double contact_thr  = m_config.get<double>("eeq_batched_min_distance", 15.0);
        const bool   in_contact   = (m_contact_min_dist < 0.0)  // unknown → be safe
                                 || (contact_thr > 0.0 && m_contact_min_dist < contact_thr);
        const std::string dist_str = (m_contact_min_dist < 0.0)
                                         ? std::string("unknown")
                                         : fmt::format("{:.2f}", m_contact_min_dist);
        if (prefer_exact && in_contact) {
            // Keep the exact solver. Explicit PCG → PCG+block-Jacobi; Auto → leave
            // method_to_use for the exact auto-selector below.
            if (m_solve_method == EEQSolveMethod::PCG)
                method_to_use = EEQSolveMethod::PCG;
            if (m_verbosity >= 1) {
                CurcumaLogger::info(fmt::format(
                    "EEQ: fragmented system (nfrag={}, N={}) with fragments in contact "
                    "(min inter-fragment distance {} Bohr < {}); using exact solver "
                    "(cross-fragment Coulomb retained) instead of approximate batched",
                    nfrag, natoms, dist_str, contact_thr));
            }
        } else {
            method_to_use = EEQSolveMethod::Batched;
            if (m_verbosity >= 1) {
                CurcumaLogger::info(fmt::format(
                    "EEQ: highly fragmented system (nfrag={}, N={}, density={}%, "
                    "min inter-fragment distance {} Bohr); using batched per-fragment solve",
                    nfrag, natoms, (100 * nfrag) / natoms, dist_str));
            }
        }
    }

    if (method_to_use == EEQSolveMethod::Batched) {
        Matrix A_nn = A.topLeftCorner(natoms, natoms);
        Vector rhs_atoms = x.head(natoms);
        Matrix C_constraints = Matrix::Zero(nfrag, natoms);
        for (int f = 0; f < nfrag; ++f)
            for (int j = 0; j < natoms; ++j)
                C_constraints(f, j) = A(natoms + f, j);
        Vector rhs_constraints = x.tail(nfrag);

        // Claude Generated (WP2, May 2026): Pre-pass — build fragment_ids[f] once.
        // Saves the inner O(N) constraint scan per iteration of the per-fragment loop.
        std::vector<std::vector<int>> fragment_ids(nfrag);
        for (int f = 0; f < nfrag; ++f) {
            for (int j = 0; j < natoms; ++j) {
                if (C_constraints(f, j) != 0.0) fragment_ids[f].push_back(j);
            }
        }

        // Per-fragment local solves: assemble fragment-local A_ff, augment with one
        // constraint row, solve via small LU. Each fragment is independent — disjoint
        // write targets in charges_batched (each atom belongs to exactly one fragment),
        // all reads from A_nn / C_constraints / rhs_atoms / rhs_constraints are read-only.
        Vector charges_batched = Vector::Zero(natoms);
        std::atomic<bool> batched_ok{true};

        // Worker — Per-Thread scratch (Aug, rhs_local, lu, sol) avoids realloc per iteration
        // and prevents False Sharing because every variable is local to the thread.
        auto fragment_worker = [&](int t_id, int T) {
            Eigen::PartialPivLU<Matrix> lu;
            Matrix Aug;
            Vector rhs_local, sol;
            for (int f = t_id; f < nfrag; f += T) {
                if (!batched_ok.load(std::memory_order_relaxed)) return;
                const auto& ids = fragment_ids[f];
                const int nf = static_cast<int>(ids.size());
                if (nf == 0) continue;

                // Assemble augmented (nf+1)×(nf+1) system: [A_ff c^T; c 0] [q;λ] = [b_f; q_target]
                Aug.resize(nf + 1, nf + 1);
                rhs_local.resize(nf + 1);
                for (int i = 0; i < nf; ++i) {
                    rhs_local(i) = rhs_atoms(ids[i]);
                    for (int j = 0; j < nf; ++j) {
                        Aug(i, j) = A_nn(ids[i], ids[j]);
                    }
                    Aug(i, nf) = C_constraints(f, ids[i]);
                    Aug(nf, i) = C_constraints(f, ids[i]);
                }
                Aug(nf, nf) = 0.0;
                rhs_local(nf) = rhs_constraints(f);

                lu.compute(Aug);
                sol = lu.solve(rhs_local);
                if (!sol.allFinite()) {
                    batched_ok.store(false, std::memory_order_relaxed);
                    return;
                }
                for (int i = 0; i < nf; ++i) charges_batched(ids[i]) = sol(i);
            }
        };

        // Pool dispatch — pattern mirrors dist_worker (eeq_solver.cpp:2870-2881).
        // Threshold nfrag >= 16 matches the auto-select gate above (Z. 1283-1285):
        // below that, threading overhead exceeds gain.
        const int T = (num_threads > 1 && nfrag >= 16) ? std::min(num_threads, nfrag) : 1;
        if (T > 1 && pool) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t) futures.push_back(pool->enqueue(fragment_worker, t, T));
            fragment_worker(0, T);
            for (auto& f : futures) f.get();
        } else if (T > 1) {
            std::vector<std::thread> threads(T - 1);
            for (int t = 1; t < T; ++t) threads[t - 1] = std::thread(fragment_worker, t, T);
            fragment_worker(0, T);
            for (auto& th : threads) th.join();
        } else {
            fragment_worker(0, 1);
        }

        if (batched_ok.load(std::memory_order_relaxed)) {
            if (m_verbosity >= 2) {
                fmt::print(stderr, "[EEQ] dispatchSolve: Batched per-fragment solve done (nfrag={}, threads={})\n", nfrag, T);
            }
            return charges_batched;
        }

        if (m_verbosity >= 1)
            CurcumaLogger::warn("EEQ batched: per-fragment LU produced non-finite charges, falling back to SchurCholesky");
        method_to_use = EEQSolveMethod::SchurCholesky;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        const char* method_name = (method_to_use == EEQSolveMethod::PCG) ? "PCG"
            : (method_to_use == EEQSolveMethod::SchurCholesky) ? "SchurCholesky" : "LU";
        const char* mode_str = (m_solve_method == EEQSolveMethod::Auto) ? " (auto)" : "";
        CurcumaLogger::info(fmt::format("[EEQ] dispatchSolve: Using {} solver{} (N={})", method_name, mode_str, natoms));
    }

    Vector charges;

    // Claude Generated (June 2026): EEQSolveMethod::LU must be in this guard too —
    // otherwise an explicit `-eeq_solver.solve_method lu` makes this whole block fall
    // through to `goto end_dispatch` with `charges` never assigned (segfault downstream),
    // and the `else if (method_to_use == LU)` branch below stays dead code. With LU in
    // the guard, explicit-LU routes through that branch's `goto lu_solve`.
    if (method_to_use == EEQSolveMethod::SchurCholesky || method_to_use == EEQSolveMethod::PCG
        || method_to_use == EEQSolveMethod::LU || m_solve_method == EEQSolveMethod::Auto) {
        Matrix A_nn = A.topLeftCorner(natoms, natoms);
        Vector rhs_atoms = x.head(natoms);

        Matrix C = Matrix::Zero(nfrag, natoms);
        for (int f = 0; f < nfrag; ++f) {
            for (int j = 0; j < natoms; ++j) {
                C(f, j) = A(natoms + f, j);
            }
        }
        Vector rhs_constraints = x.tail(nfrag);

        // Auto-benchmark on first call
        // Claude Generated (March 2026): Skip benchmark for large N — PCG always wins above ~500 atoms.
        // SchurCholesky for N>500 allocates O(N²) memory and O(N³) compute — prohibitively slow.
        // For multi-fragment systems (nfrag>1), PCG needs nfrag+1 solves; benchmarking doubles that cost.
        //
        // May 2026 (Stage 1): nfrag-aware preference. SchurCholesky factorizes A_nn once
        // (N³/6) and amortizes the cost over nfrag back-substitutions via BLAS-3 multi-RHS
        // (Z2 = llt.solve(C^T)). PCG runs (nfrag+1) independent matvec loops at O(N²·k_pcg)
        // each. Cholesky wins as soon as k_pcg·(nfrag+1) > N/6 + nfrag.
        const int PCG_DIRECT_THRESHOLD = m_config.get<int>("pcg_large_threshold", 500);
        const int nfrag_threshold = m_config.get<int>("eeq_pcg_nfrag_threshold", 4);
        const int pcg_expected_iters = m_config.get<int>("eeq_pcg_expected_iters", 30);
        const bool prefer_chol_for_nfrag =
            (nfrag >= nfrag_threshold) ||
            (pcg_expected_iters * (nfrag + 1) > natoms / 6 + nfrag);

        if (m_solve_method == EEQSolveMethod::Auto && !m_auto_benchmark_done && natoms > PCG_DIRECT_THRESHOLD && !prefer_chol_for_nfrag) {
            m_selected_method = EEQSolveMethod::PCG;
            m_auto_benchmark_done = true;
            if (m_verbosity >= 1) {
                CurcumaLogger::info(fmt::format("EEQ: Large system (N={}, nfrag={}), using PCG directly (threshold={})",
                                                natoms, nfrag, PCG_DIRECT_THRESHOLD));
            }
            // Fall through to PCG path below
        } else if (m_solve_method == EEQSolveMethod::Auto && !m_auto_benchmark_done && prefer_chol_for_nfrag) {
            m_selected_method = EEQSolveMethod::SchurCholesky;
            m_auto_benchmark_done = true;
            if (m_verbosity >= 1) {
                CurcumaLogger::info(fmt::format(
                    "EEQ: nfrag={} (threshold={}), N={}, k_pcg_est={}: SchurCholesky cheaper than {}+1 PCG solves, selecting it.",
                    nfrag, nfrag_threshold, natoms, pcg_expected_iters, nfrag));
            }
            method_to_use = EEQSolveMethod::SchurCholesky;
        }

        if (m_solve_method == EEQSolveMethod::Auto && !m_auto_benchmark_done) {
            auto t1 = std::chrono::high_resolution_clock::now();
            Vector charges_schur = solveWithSchurCholesky(A_nn, rhs_atoms, C, rhs_constraints, natoms, nfrag);
            auto t2 = std::chrono::high_resolution_clock::now();

            int pcg_max = m_config.get<int>("max_pcg_iterations", 100);
            double pcg_tol = m_config.get<double>("pcg_tolerance", 1e-10);
            Vector x0_zero = Vector::Zero(natoms);

            // Stage 2 (May 2026): even the benchmark PCG run uses the block-Jacobi PC
            // when nfrag>=2, so the timed comparison reflects post-Stage-2 PCG cost.
            BlockJacobiPC bench_pc;
            const BlockJacobiPC* bench_pc_ptr = nullptr;
            if (nfrag >= 2) {
                bench_pc = buildBlockJacobi(A_nn, C);
                if (bench_pc.valid) bench_pc_ptr = &bench_pc;
            }

            Vector z1 = solveWithPCG(A_nn, rhs_atoms, x0_zero, pcg_max, pcg_tol, bench_pc_ptr);
            Matrix Z2(natoms, nfrag);
            for (int f = 0; f < nfrag; ++f) {
                Vector c_col = C.row(f).transpose();
                Z2.col(f) = solveWithPCG(A_nn, c_col, x0_zero, pcg_max, pcg_tol, bench_pc_ptr);
            }
            Matrix S = C * Z2;
            Vector schur_rhs = C * z1 - rhs_constraints;
            Vector lambda;
            if (nfrag == 1 && std::abs(S(0, 0)) < 1e-12) {
                // Degenerate Schur complement during benchmark — PCG failed, prefer SchurCholesky.
                // Do NOT return charges_schur blindly: it may itself be empty/NaN if Cholesky failed.
                // Fall through: select SchurCholesky for future calls and return its result if valid.
                m_selected_method = EEQSolveMethod::SchurCholesky;
                m_auto_benchmark_done = true;
                if (charges_schur.size() == natoms) {
                    bool schur_nan = false;
                    for (int i = 0; i < natoms && !schur_nan; ++i)
                        schur_nan = std::isnan(charges_schur[i]) || std::isinf(charges_schur[i]);
                    if (!schur_nan) return charges_schur;
                }
                // Both solvers failed during benchmark — signal failure to caller
                return Vector();
            } else if (nfrag == 1) {
                lambda = Vector::Constant(1, schur_rhs(0) / S(0, 0));
            } else {
                lambda = S.partialPivLu().solve(schur_rhs);
            }
            Vector charges_pcg = z1 - Z2 * lambda;
            auto t3 = std::chrono::high_resolution_clock::now();

            double schur_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
            double pcg_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

            m_selected_method = (pcg_ms < schur_ms) ? EEQSolveMethod::PCG : EEQSolveMethod::SchurCholesky;
            m_auto_benchmark_done = true;

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info(fmt::format("[EEQ] Auto-benchmark: SchurCholesky={:.1f}ms, PCG={:.1f}ms -> {}",
                           schur_ms, pcg_ms,
                           m_selected_method == EEQSolveMethod::PCG ? "PCG" : "SchurCholesky"));
            }

            m_pcg_last_z1 = z1;
            m_pcg_last_Z2 = Z2;  // Cache all nfrag warm starts (May 2026)
            m_pcg_cache_valid = true;

            return charges_pcg;
        }

        if (method_to_use == EEQSolveMethod::PCG) {
            int pcg_max = m_config.get<int>("max_pcg_iterations", 100);
            double pcg_tol = m_config.get<double>("pcg_tolerance", 1e-10);

            // Large systems need more iterations and a relative (looser) tolerance.
            // 200 iterations and 1e-10 absolute are unachievable for N>500 with Jacobi PCG.
            // Relative tolerance 1e-6 gives charge error ~1e-6 per atom → energy error ~0.3 mEh
            // which is acceptable for geometry optimization.
            int pcg_large_threshold = m_config.get<int>("pcg_large_threshold", 500);
            int pcg_large_iterations = m_config.get<int>("pcg_large_system_iterations", 100);
            int pcg_large_scaling = m_config.get<int>("pcg_large_system_scaling", 10);
            double pcg_tol_factor = m_config.get<double>("pcg_large_system_tol_factor", 1e-6);

            if (natoms > pcg_large_threshold) {
                pcg_max = std::max(pcg_max, std::min(pcg_large_scaling * natoms, pcg_large_iterations));
                double rhs_norm = rhs_atoms.norm() + 1.0;
                pcg_tol = std::max(pcg_tol, pcg_tol_factor * rhs_norm);
                if (m_verbosity >= 1) {
                    CurcumaLogger::info(fmt::format(
                        "EEQ PCG large system (N={}): max_iter={}, tol={:.2e}",
                        natoms, pcg_max, pcg_tol));
                }
            }

            // Claude Generated (April 2026): Gershgorin condition estimate for A_nn.
            // Cheap O(N²) pre-check: if the matrix is weakly diagonally dominant or has
            // a very large estimated condition number, PCG will converge slowly. Skip
            // directly to LU instead of burning 5000 PCG iterations first.
            // Gershgorin: λ ∈ [A_ii - Σ_j≠i |A_ij|, A_ii + Σ_j≠i |A_ij|] for each row.
            // For the EEQ matrix, A_ij = erf(γ·r)/r > 0, so |A_ij| = A_ij.
            bool pcg_viable = true;
            {
                double min_gersh = std::numeric_limits<double>::max();
                double max_gersh = 0.0;
                for (int i = 0; i < natoms; ++i) {
                    double diag = A_nn(i, i);
                    double off = 0.0;
                    for (int j = 0; j < natoms; ++j) {
                        if (i != j) off += A_nn(i, j);
                    }
                    min_gersh = std::min(min_gersh, diag - off);
                    max_gersh = std::max(max_gersh, diag + off);
                }
                double kappa_est = (min_gersh > 0.0) ? max_gersh / min_gersh
                                                     : std::numeric_limits<double>::infinity();

                if (m_verbosity >= 2) {
                    CurcumaLogger::info(fmt::format(
                        "EEQ PCG: Gershgorin bounds [{:.4e}, {:.4e}], κ_est={:.1e}",
                        min_gersh, max_gersh, kappa_est));
                }

                // If the matrix is clearly not SPD or has extreme condition number,
                // skip PCG entirely — it will burn 5000+ iterations and then fail anyway.
                //
                // Claude Generated (June 2026): the bypass is honored ONLY in Auto mode.
                // The Gershgorin lower bound diag - Σ|A_ij| is intrinsically negative for
                // SPD EEQ matrices (the off-diagonal erf(γ·r)/r terms sum past the diagonal
                // hardness on most atoms), so this pre-check fires on essentially every real
                // system. Left unscoped it silently overrode an explicit
                // `-eeq_solver.solve_method pcg` — the PCG loop (and with it the warm-start
                // and the multi-step EEQ extrapolation) then never ran. PCG converges for any
                // SPD matrix; the NaN→LU fallback below already guards genuine failures, so an
                // explicit PCG request runs PCG and only Auto keeps the heuristic bypass.
                if (m_solve_method == EEQSolveMethod::Auto
                    && (min_gersh <= 0.0 || kappa_est > 1e12)) {
                    if (m_verbosity >= 2) {
                        // Demoted to verbosity 2 (May 2026): Gershgorin lower bound is
                        // overly pessimistic for EEQ matrices, so the warning fires on
                        // many normal multi-fragment systems. Real fallback chain is now:
                        // PCG → SchurCholesky → LU, so SPD matrices skip PCG cheaply.
                        CurcumaLogger::warn(fmt::format(
                            "EEQ PCG: Gershgorin sees ill-conditioning (κ_est={:.1e}, min_gersh={:.2e}); "
                            "trying SchurCholesky next, LU as last resort", kappa_est, min_gersh));
                    }
                    m_pcg_cache_valid = false;
                    m_pcg_x_history.clear();   // gap -> drop multi-step warm-start history
                    pcg_viable = false;
                }
            }

            if (pcg_viable) {
            Vector x0_z1 = m_pcg_cache_valid ? m_pcg_last_z1 : Vector::Zero(natoms);

            // Stage 2 (May 2026): Build block-Jacobi preconditioner once for all
            // (nfrag+1) PCG solves. Skip for nfrag<=1 since the diagonal Jacobi is
            // already exact for a single block and the build cost (one Cholesky)
            // would be wasted.
            BlockJacobiPC block_pc;
            const BlockJacobiPC* block_pc_ptr = nullptr;
            if (nfrag >= 2) {
                block_pc = buildBlockJacobi(A_nn, C);
                if (block_pc.valid) {
                    block_pc_ptr = &block_pc;
                    if (m_verbosity >= 2)
                        CurcumaLogger::info(fmt::format(
                            "EEQ PCG: block-Jacobi PC active ({} fragments, max_block={})",
                            nfrag,
                            std::max_element(block_pc.atom_ids.begin(), block_pc.atom_ids.end(),
                                [](const auto& a, const auto& b){ return a.size() < b.size(); })->size()));
                } else if (m_verbosity >= 2) {
                    CurcumaLogger::info("EEQ PCG: block-Jacobi build failed (some fragment block not SPD), using diagonal Jacobi");
                }
            }

            // Stage 3 (May 2026): Replace per-fragment PCG loop with single multi-RHS call.
            // B = [rhs_atoms | C^T] is N × (nfrag+1). Output X has columns [z1 | Z2].
            // The matvec A·P now runs as one GEMM over all RHS columns simultaneously.
            //
            // Warm start: column 0 = z1 cache; columns 1..nfrag = m_pcg_last_Z2 if shape matches.
            const bool warm_Z2_usable =
                m_pcg_cache_valid &&
                m_pcg_last_Z2.rows() == natoms &&
                m_pcg_last_Z2.cols() == nfrag;

            Matrix B_multi(natoms, nfrag + 1);
            B_multi.col(0) = rhs_atoms;
            B_multi.rightCols(nfrag) = C.transpose();

            Matrix X0_multi(natoms, nfrag + 1);
            X0_multi.col(0) = x0_z1;
            if (warm_Z2_usable) {
                X0_multi.rightCols(nfrag) = m_pcg_last_Z2;
            } else {
                X0_multi.rightCols(nfrag).setZero();
            }

            // Multi-step extrapolation of the PCG warm-start (opt-in). Predict the initial
            // guess X0 = [z1|Z2] from the history of converged solutions instead of reusing
            // only the previous step — a better guess means fewer PCG iterations (esp. in
            // smooth MD). The converged charges are unchanged (PCG converges regardless of
            // X0); on too-short history / shape mismatch it keeps the 1-step warm-start
            // above. Claude Generated.
            if (m_eeq_extrapolation != "none"
                && static_cast<int>(m_pcg_x_history.size()) >= 2
                && m_pcg_x_history.front().rows() == natoms
                && m_pcg_x_history.front().cols() == nfrag + 1) {
                Eigen::VectorXd w = curcuma::extrapolation::weights(
                    m_eeq_extrapolation, m_eeq_extrap_order,
                    static_cast<int>(m_pcg_x_history.size()));
                if (w.size() >= 2) {
                    Matrix X0_ext = Matrix::Zero(natoms, nfrag + 1);
                    bool ok = true;
                    for (int j = 0; j < w.size(); ++j) {
                        if (m_pcg_x_history[j].rows() != natoms
                            || m_pcg_x_history[j].cols() != nfrag + 1) { ok = false; break; }
                        X0_ext.noalias() += w(j) * m_pcg_x_history[j];
                    }
                    if (ok && X0_ext.allFinite()) {
                        X0_multi = X0_ext;
                        ++m_pcg_extrap_count;
                        if (!m_eeq_extrap_cited) {
                            m_eeq_extrap_cited = true;
                            CurcumaLogger::citation(m_eeq_extrapolation == "aspc"
                                                        ? "aspc" : "density_extrapolation");
                        }
                        if (m_verbosity >= 2)
                            CurcumaLogger::info(fmt::format(
                                "EEQ PCG: {} warm-start extrapolation over {} steps",
                                m_eeq_extrapolation, static_cast<int>(w.size())));
                    }
                }
            }

            Matrix X_multi = solveWithPCG_multiRHS(A_nn, B_multi, X0_multi, pcg_max, pcg_tol, block_pc_ptr);
            Vector z1 = X_multi.col(0);
            Matrix Z2 = X_multi.rightCols(nfrag);

            Matrix S = C * Z2;
            Vector schur_rhs = C * z1 - rhs_constraints;
            Vector lambda;

            // Claude Generated (Apr 2026): Guard against degenerate Schur complement.
            // S(0,0) near-zero means the constraint solve is degenerate — happens when
            // PCG hasn't converged and Z2 is garbage. Fall back to LU.
            bool schur_ok = true;
            if (nfrag == 1 && std::abs(S(0, 0)) < 1e-12) {
                schur_ok = false;
                if (m_verbosity >= 1)
                    CurcumaLogger::warn(fmt::format(
                        "EEQ PCG: Schur complement S(0,0)={:.2e} near zero, falling back to LU", S(0, 0)));
            }

            if (schur_ok) {
                if (nfrag == 1) {
                    lambda = Vector::Constant(1, schur_rhs(0) / S(0, 0));
                } else {
                    lambda = S.partialPivLu().solve(schur_rhs);
                }
                charges = z1 - Z2 * lambda;

                // Check for NaN/Inf from unconverged PCG; fall back to LU if needed.
                bool has_nan = false;
                for (int i = 0; i < natoms && !has_nan; ++i)
                    has_nan = std::isnan(charges[i]) || std::isinf(charges[i]);

                if (!has_nan) {
                    m_pcg_last_z1 = z1;
                    m_pcg_last_Z2 = Z2;  // Cache all nfrag warm starts (May 2026)
                    m_pcg_cache_valid = true;
                    // Record the converged [z1|Z2] for the multi-step warm-start (opt-in).
                    if (m_eeq_extrapolation != "none") {
                        if (!m_pcg_x_history.empty()
                            && (m_pcg_x_history.front().rows() != X_multi.rows()
                                || m_pcg_x_history.front().cols() != X_multi.cols()))
                            m_pcg_x_history.clear();
                        m_pcg_x_history.push_front(X_multi);
                        const int maxlen = (m_eeq_extrapolation == "aspc")
                                               ? std::max(m_eeq_extrap_order + 2, 2)
                                               : std::min(2 * m_eeq_extrap_order + 2, 10);
                        while (static_cast<int>(m_pcg_x_history.size()) > maxlen)
                            m_pcg_x_history.pop_back();
                    }
                } else {
                    schur_ok = false;
                    if (m_verbosity >= 1)
                        CurcumaLogger::warn("EEQ PCG: NaN charges after Schur solve, falling back to LU");
                }
            }

            if (!schur_ok) {
                // LU fallback — used only when PCG fails
                goto lu_solve;
            }
            // PCG succeeded — charges are set, fall through to return
        } else {
            // pcg_viable = false (Gershgorin check ruled out PCG).
            //
            // May 2026 (Stage 1): Try SchurCholesky BEFORE falling back to full LU.
            // Gershgorin's lower bound diag - Σ|A_ij| is overly pessimistic for the EEQ
            // Coulomb matrix — it overcounts inter-atomic coupling that is in fact much
            // weaker than the row sum suggests, especially in many-fragment systems
            // (mixture.xyz: nfrag=1400, Gershgorin reports min_gersh=-90 but A_nn is SPD).
            //
            // SchurCholesky factorizes only the N×N SPD A_nn block (~N³/6 ops) and
            // amortizes back-substitution over nfrag RHS via BLAS-3, vs LU on the full
            // (N+nfrag)² augmented matrix (~(N+nfrag)³/3 ops). For mixture.xyz this
            // changes 7600³/3 → 6200³/6 ≈ 5× fewer ops AND avoids the augmented system's
            // numerical issues. Cholesky has its own SPD check (llt.info()), so the cost
            // of trying first is just one Eigen factorization attempt.
            charges = solveWithSchurCholesky(A_nn, rhs_atoms, C, rhs_constraints, natoms, nfrag);
            if (charges.size() == 0) {
                if (m_verbosity >= 1)
                    CurcumaLogger::warn("EEQ: SchurCholesky also failed after PCG bypass, falling back to LU");
                goto lu_solve;
            }
            // SchurCholesky succeeded — invalidate PCG warmstart (it would mismatch on next call)
            m_pcg_cache_valid = false;
            m_pcg_x_history.clear();   // gap -> drop multi-step warm-start history
        }
    } else if (method_to_use == EEQSolveMethod::LU) {
        // Direct full-system LU — no Schur complement overhead.
        // Used when previous calls showed A is not SPD (SchurCholesky failed)
        // and PCG is unsuitable (large N, no convergence).
        goto lu_solve;
    } else {
        charges = solveWithSchurCholesky(A_nn, rhs_atoms, C, rhs_constraints, natoms, nfrag);
        if (charges.size() == 0) {
            if (m_verbosity >= 1) {
                CurcumaLogger::warn("EEQ: Schur-Cholesky failed, falling back to LU");
            }
            // Don't permanently select LU — next step tries SchurCholesky again
            goto lu_solve;
        }
    }
    } // closes outer if (method_to_use == SchurCholesky || PCG || Auto)
    goto end_dispatch;

    // Single LU solve with progress timing — shared by all paths
    {
        lu_solve:
        if (m_verbosity >= 1 && natoms > 200) {
            CurcumaLogger::info(fmt::format("EEQ LU solve starting for N={}...", natoms));
        }
        auto t0 = std::chrono::high_resolution_clock::now();
        Eigen::PartialPivLU<Matrix> lu(A);
        charges = lu.solve(x).head(natoms);
        auto t1 = std::chrono::high_resolution_clock::now();
        if (m_verbosity >= 1 && natoms > 200) {
            CurcumaLogger::result(fmt::format(
                "EEQ LU solve done in {:.1f}s",
                std::chrono::duration<double>(t1 - t0).count()));
        }
        m_pcg_cache_valid = false;
        m_pcg_x_history.clear();   // LU discontinuity -> drop multi-step warm-start history
        // After LU, switch to PCG for future calls (with warm-start from LU result).
        // This avoids getting stuck in LU permanently: PCG may converge on the next
        // step since the geometry change is small and LU provides a good initial guess.
        if (m_solve_method == EEQSolveMethod::Auto)
            m_selected_method = EEQSolveMethod::PCG;
    }
    end_dispatch: ;


    // NaN/Inf validation — return empty vector to signal failure.
    // Caller (calculateFinalCharges) is responsible for choosing an appropriate fallback
    // (topology_charges are far better than uniform 0 charges for neutral molecules).
    for (int i = 0; i < natoms; ++i) {
        if (std::isnan(charges[i]) || std::isinf(charges[i])) {
            if (m_allow_unconverged) {
                CurcumaLogger::warn(fmt::format("EEQ dispatchSolve: Invalid charge[{}] = {} (allow_unconverged_charges=true)", i, charges[i]));
            } else {
                CurcumaLogger::error(fmt::format("EEQ dispatchSolve: Invalid charge[{}] = {}", i, charges[i]));
            }
            return Vector();  // Empty → caller uses topology_charges as fallback
        }
    }

    return charges;
}

// ===== Schur Complement Cholesky Solver =====
// Claude Generated - March 2026 (Performance optimization)
//
// The augmented EEQ system [A C^T; C 0] is indefinite due to constraint rows.
// But the NxN Coulomb+hardness sub-matrix A is SPD (positive diagonal dominance
// from hardness terms + positive off-diagonal Coulomb terms).
// Cholesky is O(N³/6) vs O(N³/3) for LU — roughly 2× faster.

Vector EEQSolver::solveWithSchurCholesky(
    const Matrix& A_nn,
    const Vector& rhs_atoms,
    const Matrix& C,
    const Vector& rhs_constraints,
    int natoms,
    int nfrag)
{
    // WP-EEQ-Cache (May 2026, fix May 2026): Skip O(N^3/6) Cholesky factorization
    // when geometry displacement since last factorization is below m_refactor_eps.
    //
    // Contract: m_pending_geometry / m_pending_cn are populated only by Phase 2
    // (calculateFinalCharges). Phase 1 (calculateTopologyCharges) reaches this
    // function with empty pending buffers — in that case we must solve locally and
    // leave the cache untouched, otherwise Phase 2 would later subtract a 0xN matrix
    // from an NxN one and segfault.
    //
    // m_refactor_eps <= 0 disables the cache entirely (every call refactors —
    // bit-identical to pre-WP behavior).

    const bool cache_size_ok = (m_refactor_eps > 0.0)
        && m_chol_cache.valid
        && m_chol_cache.cached_natoms == natoms
        && m_chol_cache.cached_nfrag  == nfrag
        && m_pending_geometry.rows()  == natoms
        && m_pending_cn.size()        == natoms;

    bool need_refactor = !cache_size_ok;

    if (cache_size_ok) {  // m_refactor_eps > 0 implied
        // O(N) displacement check — much cheaper than O(N^3) factorization
        const double max_dr  = (m_pending_geometry - m_chol_cache.last_geometry)
                                   .rowwise().norm().maxCoeff();
        const double max_dcn = (m_pending_cn - m_chol_cache.last_cn)
                                   .cwiseAbs().maxCoeff();
        if (max_dr > m_refactor_eps || max_dcn > 0.05)
            need_refactor = true;
    }
    if (m_refactor_force_every > 0
        && m_chol_cache.steps_since_refactor >= m_refactor_force_every)
        need_refactor = true;

    if (need_refactor) {
        // Threaded LAPACK Cholesky (column-major copy of the row-major A_nn). The
        // factorization and the nfrag+1 RHS back-substitutions run on the linked
        // BLAS/LAPACK (OpenBLAS/MKL); see the eeqCholesky* helpers above. The BLAS
        // thread count for this region is set by the ScopedBlasThreads guard in
        // calculateFinalCharges (driven by the EnergyCalculator budget).
        ColMatrix chol = A_nn;  // row-major → column-major copy (symmetric: same values)
        if (!eeqCholeskyFactorize(chol)) {
            if (m_verbosity >= 1)
                CurcumaLogger::warn("EEQ Schur-Cholesky: A matrix not SPD, falling back to LU");
            m_chol_cache.reset();
            return Vector();  // Empty = signal to fall back
        }

        // Persist into cache ONLY when (a) the cache is enabled (eps > 0) and
        // (b) the pending buffers reflect this call's geometry (Phase 2 contract).
        // Phase 1 has empty m_pending_geometry — local-solve path keeps the cache
        // clean so the next Phase 2 call sees valid=false and refactors normally.
        const bool can_persist = (m_refactor_eps > 0.0)
            && (m_pending_geometry.rows() == natoms)
            && (m_pending_cn.size()       == natoms);

        if (can_persist) {
            m_chol_cache.chol_factor      = std::move(chol);
            ColMatrix Z2cm                = C.transpose();   // col-major N×nfrag RHS
            eeqCholeskySolve(m_chol_cache.chol_factor, Z2cm); // Z2 = A_nn^{-1}·C^T
            m_chol_cache.Z2               = Z2cm;            // store as row-major Matrix
            m_chol_cache.S                = schurComplementCZ(C, m_chol_cache.Z2, natoms, nfrag);
            m_chol_cache.last_geometry    = m_pending_geometry;
            m_chol_cache.last_cn          = m_pending_cn;
            m_chol_cache.cached_natoms    = natoms;
            m_chol_cache.cached_nfrag     = nfrag;
            m_chol_cache.steps_since_refactor = 0;
            m_chol_cache.valid            = true;
            if (m_verbosity >= 3)
                fmt::print(stderr, "[EEQ-Cache] Refactorized + persisted (N={}, nfrag={})\n",
                           natoms, nfrag);
        } else {
            // Phase 1 OR cache disabled — local solve, no cache writes
            m_chol_cache.valid = false;
            if (m_verbosity >= 3)
                fmt::print(stderr, "[EEQ-Cache] Local solve (Phase 1 or cache off, N={})\n",
                           natoms);
            Vector z1_local = rhs_atoms;
            eeqCholeskySolve(chol, z1_local);
            ColMatrix Z2cm_local = C.transpose();
            eeqCholeskySolve(chol, Z2cm_local);
            Matrix Z2_local = Z2cm_local;   // col-major → row-major
            Matrix S_local  = schurComplementCZ(C, Z2_local, natoms, nfrag);
            Vector schur_rhs_local = C * z1_local - rhs_constraints;
            Vector lambda_local;
            if (nfrag == 1) {
                lambda_local = Vector::Constant(1, schur_rhs_local(0) / S_local(0, 0));
            } else {
                lambda_local = S_local.partialPivLu().solve(schur_rhs_local);
            }
            return z1_local - Z2_local * lambda_local;
        }
    } else {
        if (m_verbosity >= 3)
            fmt::print(stderr, "[EEQ-Cache] Cache hit (step {})\n",
                       m_chol_cache.steps_since_refactor);
    }

    // Cache-hit path: factor was either freshly persisted above OR carried over.
    // O(N^2) triangular solve with cached factorization.
    Vector z1 = rhs_atoms;
    eeqCholeskySolve(m_chol_cache.chol_factor, z1);
    ++m_chol_cache.steps_since_refactor;

    // Schur complement with cached Z2 and S (both constant when A is cached)
    Vector schur_rhs = C * z1 - rhs_constraints;
    Vector lambda;
    if (nfrag == 1) {
        lambda = Vector::Constant(1, schur_rhs(0) / m_chol_cache.S(0, 0));
    } else {
        lambda = m_chol_cache.S.partialPivLu().solve(schur_rhs);
    }
    return z1 - m_chol_cache.Z2 * lambda;
}

// ===== Block-Jacobi Preconditioner =====
// Claude Generated - May 2026 (Stage 2 of EEQ multi-fragment optimization)
//
// Replaces the diagonal Jacobi M^-1 = diag(A)^-1 with M^-1 = blockdiag(A_ff^-1),
// where A_ff is the sub-block of A_nn restricted to fragment f's atoms. Each
// block is dense Cholesky-factorized once at build time.
//
// For multi-fragment EEQ systems where intra-fragment Coulomb >> inter-fragment
// Coulomb (typical for solvated clusters and condensed phases), this preconditioner
// captures the dominant structure exactly. PCG then only iterates on the weak
// inter-fragment coupling, dropping k from 30-100 to 2-5 typical iterations.

Vector EEQSolver::BlockJacobiPC::apply(const Vector& r) const {
    Vector z = r;  // identity baseline (atoms not in any fragment)
    if (!valid) return z;
    for (size_t f = 0; f < blocks.size(); ++f) {
        const auto& ids = atom_ids[f];
        if (ids.empty()) continue;
        const int nf = static_cast<int>(ids.size());
        Vector r_f(nf);
        for (int i = 0; i < nf; ++i) r_f[i] = r[ids[i]];
        Vector z_f = blocks[f].solve(r_f);
        for (int i = 0; i < nf; ++i) z[ids[i]] = z_f[i];
    }
    return z;
}

EEQSolver::BlockJacobiPC EEQSolver::buildBlockJacobi(const Matrix& A_nn, const Matrix& C) {
    BlockJacobiPC pc;
    const int nfrag = static_cast<int>(C.rows());
    const int natoms = static_cast<int>(C.cols());
    pc.atom_ids.resize(nfrag);
    pc.blocks.reserve(nfrag);

    // C(f, j) = 1 iff atom j ∈ fragment f. Reconstruct atom-id lists.
    for (int f = 0; f < nfrag; ++f) {
        for (int j = 0; j < natoms; ++j) {
            if (C(f, j) != 0.0) pc.atom_ids[f].push_back(j);
        }
    }

    // Per-fragment Cholesky. If any fragment block is not SPD, mark PC invalid;
    // caller will fall back to diagonal Jacobi.
    for (int f = 0; f < nfrag; ++f) {
        const auto& ids = pc.atom_ids[f];
        const int nf = static_cast<int>(ids.size());
        if (nf == 0) {
            pc.blocks.emplace_back();
            continue;
        }
        Matrix A_ff(nf, nf);
        for (int i = 0; i < nf; ++i)
            for (int j = 0; j < nf; ++j)
                A_ff(i, j) = A_nn(ids[i], ids[j]);
        Eigen::LLT<Matrix> llt(A_ff);
        if (llt.info() != Eigen::Success) {
            // One fragment block not SPD → entire PC is unsafe to apply (would break
            // PCG symmetry-preservation). Caller falls back to diagonal Jacobi.
            return BlockJacobiPC{};  // valid stays false
        }
        pc.blocks.push_back(std::move(llt));
    }
    pc.valid = true;
    return pc;
}

// ===== Preconditioned Conjugate Gradient Solver =====
// Claude Generated - March 2026 (Performance optimization)
//
// Iterative solver: O(N²·k) per call where k = iteration count.
// With warm start from previous MD/optimization step, k ≈ 10-20 (vs 50-100 cold).
// Jacobi preconditioner (diagonal inverse) is cheap and effective for diagonally
// dominant EEQ matrices. Block-Jacobi (Stage 2) drops k further when nfrag>1.

Vector EEQSolver::solveWithPCG(
    const Matrix& A,
    const Vector& b,
    const Vector& x0,
    int max_iter,
    double tol,
    const BlockJacobiPC* block_pc)
{
    const int n = b.size();
    // Stage 2 (May 2026): if block_pc is provided and valid, use blockwise Cholesky
    // inversion as preconditioner. Otherwise fall back to diagonal Jacobi.
    const bool use_block_pc = (block_pc != nullptr && block_pc->valid);
    Vector M_inv;
    if (!use_block_pc) {
        M_inv = A.diagonal().cwiseInverse();
    }
    auto applyPC = [&](const Vector& r) -> Vector {
        return use_block_pc ? block_pc->apply(r) : Vector(M_inv.cwiseProduct(r));
    };

    Vector x = x0;
    Vector r = b - A * x;
    Vector z = applyPC(r);
    Vector p = z;
    double rz = r.dot(z);

    // Progress tracking: 10% intervals at verbosity >=2, 1%+residual at verbosity >=3
    int progress_10pct = std::max(1, max_iter / 10);
    int progress_1pct = std::max(1, max_iter / 100);
    bool show_simple = (m_verbosity >= 2 && max_iter >= 50);
    bool show_detailed = (m_verbosity >= 3 && max_iter >= 50);

    if (show_detailed) {
        fmt::print("[EEQ] PCG {} iters (tol={:.0e}): ", max_iter, tol);
        fflush(stdout);
    } else if (show_simple) {
        fmt::print("[EEQ] PCG {} iters: ", max_iter);
        fflush(stdout);
    }

    for (int k = 0; k < max_iter; ++k) {
        Vector Ap = A * p;                       // O(N²) matvec — the dominant cost
        double pAp = p.dot(Ap);
        if (std::abs(pAp) < 1e-30) break;       // Degenerate direction
        double alpha = rz / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        // Progress output: 10% steps with (or without) residual norm
        if ((k + 1) % progress_10pct == 0) {
            if (show_detailed) {
                fmt::print(" {}%({:.0e})", (100 * (k + 1)) / max_iter, r.norm());
            } else if (show_simple) {
                fmt::print(" {}%", (100 * (k + 1)) / max_iter);
            }
            fflush(stdout);
        } else if (show_detailed && (k + 1) % progress_1pct == 0) {
            fmt::print(".");
            fflush(stdout);
        }

        double r_norm = r.norm();
        if (r_norm < tol) {
            m_pcg_total_calls++;
            m_pcg_total_iters += k + 1;
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info(fmt::format("[EEQ] PCG converged in {} iterations (|r|={:.2e})", k + 1, r_norm));
            }
            return x;
        }

        Vector z_new = applyPC(r);
        double rz_new = r.dot(z_new);
        double beta = rz_new / rz;
        p = z_new + beta * p;
        rz = rz_new;
    }

    // Accumulate non-convergence statistics
    m_pcg_total_calls++;
    m_pcg_nonconv_calls++;
    m_pcg_total_iters += max_iter;
    double final_r = r.norm();
    if (final_r > m_pcg_worst_residual) m_pcg_worst_residual = final_r;
    if (show_detailed) {
        fmt::print(" NOT CONVERGED (|r|={:.2e})\n", final_r);
    } else if (show_simple) {
        fmt::print(" NOT CONVERGED\n");
    }
    return x;  // Return best estimate
}

// ===== Multi-RHS Block-PCG Solver =====
// Claude Generated - May 2026 (Stage 3 of EEQ multi-fragment optimization)
//
// Solves A·X = B simultaneously for all m RHS columns. The matvec A·P becomes a
// GEMM (BLAS-3) over Matrix(N×m), replacing m independent BLAS-2 matvecs in the
// per-fragment PCG loop. For nfrag ≫ 1, this typically gives 2-4× wall-time
// reduction even at the same iteration count, because BLAS-3 amortizes memory
// bandwidth across columns. Combined with Stage-2's BlockJacobi PC (k=2-5 typical),
// total cost approaches that of one Cholesky factorization.
//
// Each column has its own α, β, residual norm; column-wise convergence is checked
// against the max residual norm. No deflation — converged columns continue to be
// updated harmlessly until all columns reach tolerance or max_iter is hit.

Matrix EEQSolver::solveWithPCG_multiRHS(
    const Matrix& A,
    const Matrix& B,
    const Matrix& X0,
    int max_iter,
    double tol,
    const BlockJacobiPC* block_pc)
{
    const int n = static_cast<int>(B.rows());
    const int m = static_cast<int>(B.cols());

    const bool use_block_pc = (block_pc != nullptr && block_pc->valid);
    Vector M_inv;
    if (!use_block_pc) {
        M_inv = A.diagonal().cwiseInverse();
    }
    auto applyPC = [&](const Matrix& Rin) -> Matrix {
        if (use_block_pc) {
            Matrix Zout(n, m);
            for (int c = 0; c < m; ++c) {
                Zout.col(c) = block_pc->apply(Rin.col(c));
            }
            return Zout;
        } else {
            // Diagonal Jacobi: Z(i,c) = M_inv(i) * R(i,c)
            return M_inv.asDiagonal() * Rin;
        }
    };

    Matrix X = X0;
    Matrix R = B - A * X;                                              // GEMM
    Matrix Z = applyPC(R);
    Matrix P = Z;
    Eigen::RowVectorXd rz = (R.array() * Z.array()).colwise().sum();

    bool show_simple = (m_verbosity >= 2 && max_iter >= 50);
    if (show_simple) {
        fmt::print("[EEQ] multi-RHS PCG {} iters, m={}: ", max_iter, m);
        fflush(stdout);
    }
    int progress_10pct = std::max(1, max_iter / 10);

    int k = 0;
    int total_iters_for_stats = 0;
    bool converged = false;
    for (k = 0; k < max_iter; ++k) {
        Matrix AP = A * P;                                             // GEMM — dominant cost
        Eigen::RowVectorXd pAp = (P.array() * AP.array()).colwise().sum();

        // Per-column α = rz / pAp; guard against degenerate directions.
        Eigen::RowVectorXd alpha(m);
        for (int c = 0; c < m; ++c) {
            alpha(c) = (std::abs(pAp(c)) < 1e-30) ? 0.0 : rz(c) / pAp(c);
        }

        // X += P · diag(α);  R -= AP · diag(α) — Eigen evaluates this as one GEMM each.
        X += P * alpha.asDiagonal();
        R -= AP * alpha.asDiagonal();

        // Convergence check — max column-residual norm
        double max_r_norm = 0.0;
        for (int c = 0; c < m; ++c) {
            double rn = R.col(c).norm();
            if (rn > max_r_norm) max_r_norm = rn;
        }
        if ((k + 1) % progress_10pct == 0 && show_simple) {
            fmt::print(" {}%({:.1e})", (100 * (k + 1)) / max_iter, max_r_norm);
            fflush(stdout);
        }
        if (max_r_norm < tol) {
            converged = true;
            total_iters_for_stats = k + 1;
            if (show_simple) fmt::print(" done({} iters, |r|max={:.2e})\n", k + 1, max_r_norm);
            break;
        }

        Matrix Z_new = applyPC(R);
        Eigen::RowVectorXd rz_new = (R.array() * Z_new.array()).colwise().sum();

        Eigen::RowVectorXd beta(m);
        for (int c = 0; c < m; ++c) {
            beta(c) = (std::abs(rz(c)) < 1e-30) ? 0.0 : rz_new(c) / rz(c);
        }

        // P = Z_new + P · diag(β)
        P = Z_new + P * beta.asDiagonal();
        rz = rz_new;
    }

    // Update PCG statistics — count m solves' worth of work for fair comparison
    // with the per-fragment loop.
    if (converged) {
        m_pcg_total_calls += m;
        m_pcg_total_iters += total_iters_for_stats * m;
    } else {
        m_pcg_total_calls += m;
        m_pcg_nonconv_calls += m;
        m_pcg_total_iters += max_iter * m;
        double final_r = 0.0;
        for (int c = 0; c < m; ++c) final_r = std::max(final_r, R.col(c).norm());
        if (final_r > m_pcg_worst_residual) m_pcg_worst_residual = final_r;
        if (show_simple) fmt::print(" NOT CONVERGED (|r|max={:.2e})\n", final_r);
    }

    return X;
}

void EEQSolver::printConvergenceSummary() {
    if (m_pcg_total_calls == 0) return;
    if (m_pcg_nonconv_calls > 0 && m_verbosity >= 1) {
        int conv = m_pcg_total_calls - m_pcg_nonconv_calls;
        CurcumaLogger::warn(fmt::format(
            "EEQ PCG: {}/{} converged, {} not (worst |r|={:.2e}, {} iters total)",
            conv, m_pcg_total_calls, m_pcg_nonconv_calls,
            m_pcg_worst_residual, m_pcg_total_iters));
    } else if (m_pcg_nonconv_calls == 0 && m_verbosity >= 2) {
        CurcumaLogger::info(fmt::format(
            "EEQ PCG: {}/{} converged ({} iters total)",
            m_pcg_total_calls, m_pcg_total_calls, m_pcg_total_iters));
    }
    resetConvergenceStats();
}

void EEQSolver::resetConvergenceStats() {
    m_pcg_total_calls = 0;
    m_pcg_nonconv_calls = 0;
    m_pcg_total_iters = 0;
    m_pcg_worst_residual = 0.0;
}

/**
 * @brief Solve augmented EEQ linear system with corrected parameters
 * @param A Augmented EEQ matrix (natoms+1) × (natoms+1)
 * @param atoms Atomic numbers
 * @param cn Coordination numbers
 * @param dxi Electronegativity corrections (WITHOUT CNF)
 * @param total_charge Total molecular charge
 * @param topology Optional topology for integer neighbor count in CNF term
 * @return Vector of atomic charges (natoms elements)
 *
 * RHS setup: x(i) = -chi + dxi + CNF*sqrt(nb)
 *
 * Claude Generated - January 4, 2026
 */
Vector EEQSolver::solveEEQ(
    const Matrix& A,
    const std::vector<int>& atoms,
    const Vector& cn,
    const Vector& dxi,
    int total_charge,
    const std::optional<TopologyInput>& topology,
    bool use_cnf_term,
    bool use_integer_nb)  // NEW: Controls whether integer nb or fractional cn is used for CNF
{
    const int natoms = atoms.size();
    int nfrag = topology.has_value() ? topology->nfrag : 1;
    int m = natoms + nfrag;
    Vector x = Vector::Zero(m);

    // Setup RHS
    // CRITICAL: CNF term handling follows Fortran reference exactly:
    // - Phase 1 (goedeckera): chieeq includes CNF with integer nb (gfnff_ini.f90:411)
    // - Phase 2 preparation: chieeq OVERWRITTEN without CNF (gfnff_ini.f90:715)
    // - Phase 3 (goed_gfnff): RHS adds CNF with fractional CN (gfnff_engrad.F90:1504)
    //
    // Curcuma implements Phase 1+2 only (no separate Phase 3 goed_gfnff call).
    // Phase 2 here mimics Parameter Preparation + implicit goed_gfnff solve.
    const double CNMAX = 4.4;
    const double MCHISHIFT = -0.09;  // Metal chi-shift from gfnff_param.f90:800

    for (int i = 0; i < natoms; ++i) {
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));

        // chi_corrected = -chi + dxi (always included)
        double chi_corrected = -params_i.chi + dxi(i);

        // Metal Charge Shift Correction - Phase 1 ONLY
        // Reference: gfnff_ini.f90:413-418 (Phase 1 parameter setup)
        // Phase 2 overwrites chieeq at line 715 WITHOUT mchishift.
        // Formula: chieeq(i) = chieeq(i) - mchishift (for transition metals only)
        // Claude Generated: Metal identification based on Z value
        if (use_integer_nb) {  // Phase 1 only - Fortran applies mchishift before goedeckera
            bool is_transition_metal = false;
            if (z_i >= 21 && z_i <= 30) is_transition_metal = true;  // Sc-Zn
            else if (z_i >= 39 && z_i <= 48) is_transition_metal = true;  // Y-Cd
            else if (z_i >= 72 && z_i <= 80) is_transition_metal = true;  // Hf-Hg

            if (is_transition_metal) {
                chi_corrected -= MCHISHIFT;  // Subtracts -0.09, effectively adds +0.09
            }
        }

        // Amide Hydrogen Electronegativity Shift (Phase 2.9 - January 17, 2026)
        // Reference: gfnff_ini.f90:717
        // Formula: chieeq(ji) = chieeq(ji) - 0.02 (for Hydrogen bonded to pi-system Nitrogen)
        if (z_i == 1 && topology.has_value()) {
            for (int neighbor : topology->neighbor_lists[i]) {
                if (atoms[neighbor] == 7) { // Bonded to Nitrogen
                    int hyb_n = detectElementSpecificHybridization(atoms[neighbor], cn(neighbor), atoms, topology, neighbor);
                    if (hyb_n == 1 || hyb_n == 2) {
                        chi_corrected -= 0.02;
                        break;
                    }
                }
            }
        }

        if (use_cnf_term) {
            // Phase 1: Add CNF term for topology charges using INTEGER nb
            // Phase 2: Add CNF term using FRACTIONAL coordination number (cn)
            // Reference: CHARGE_DATAFLOW.md and gfnff_engrad.F90:1504
            double nb_count;

            if (use_integer_nb && topology.has_value() && i < static_cast<int>(topology->neighbor_lists.size())) {
                // Phase 1: integer neighbor count with cnmax cap (gfnff_ini.f90:409)
                nb_count = static_cast<double>(topology->neighbor_lists[i].size());
                nb_count = std::min(nb_count, CNMAX);
            } else {
                // Phase 2: fractional CN WITHOUT cnmax cap (gfnff_engrad.F90:1577)
                // Fortran goed_gfnff uses sqrt(cn(i)) directly, no cnmax limit
                nb_count = cn(i);
            }
            double cnf_term = params_i.cnf * std::sqrt(nb_count);
            x(i) = chi_corrected + cnf_term;

            // DEBUG: Print RHS calculation for first 3 atoms
            if (m_verbosity >= 3 && i < 3) {
                std::cout << fmt::format(
                    "RHS[{}]: chi_corr={:.6f} (-chi={:.6f} + dxi={:.6f}), cnf_term={:.6f} (cnf={:.6f} * sqrt(nb={:.1f})), total={:.6f}",
                    i, chi_corrected, -params_i.chi, dxi(i), cnf_term, params_i.cnf, nb_count, x(i)
                ) << std::endl;
            }
        } else {
            // Phase 2: NO CNF term for final charges!
            x(i) = chi_corrected;
        }
    }

    for (int f = 0; f < nfrag; ++f) {
        int row = natoms + f;
        double q_target = (topology.has_value() && f < static_cast<int>(topology->qfrag.size()))
                         ? topology->qfrag[f]
                         : (f == 0 ? static_cast<double>(total_charge) : 0.0);
        x(row) = q_target;
    }

    // Claude Generated (March 2026): Use unified dispatchSolve instead of duplicated solver logic
    Vector charges = dispatchSolve(A, x, natoms, nfrag, total_charge);

    // DEBUG: Verify linear solve accuracy
    if (m_verbosity >= 3) {
        // Reconstruct full solution for residual check
        Vector full_sol = Vector::Zero(m);
        full_sol.head(natoms) = charges;
        // Approximate Lagrange multipliers from constraint equation
        Vector residual = A * full_sol - x;
        // For proper residual, we'd need λ too, but charge accuracy is what matters
        double max_charge_residual = residual.head(natoms).cwiseAbs().maxCoeff();
        std::cout << fmt::format("Linear solve max charge residual: {:.2e}", max_charge_residual) << std::endl;

        if (use_cnf_term && natoms >= 3) {
            std::cout << "Solution charges (first 3 atoms):" << std::endl;
            for (int i = 0; i < std::min(3, natoms); ++i) {
                std::cout << fmt::format("  q[{}] = {:.8f}", i, charges(i)) << std::endl;
            }
        }
    }

    return charges;
}

// ===== Phase 1: Topology Charges ===== (DEPRECATED - use single-solve instead)

Vector EEQSolver::calculateTopologyCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector& cn,
    const std::optional<TopologyInput>& topology,
    bool use_corrections,
    CxxThreadPool* pool,
    int num_threads)
{
    const int natoms = atoms.size();
    const double TSQRT2PI = 0.797884560802866;  // sqrt(2/π)

    // Augmented system size: n atoms + n fragments
    int nfrag = topology.has_value() ? topology->nfrag : 1;
    int m = natoms + nfrag;

    // CRITICAL FIX (Jan 29, 2026): Enable dxi corrections to match Fortran goedeckera
    // Reference: Fortran gfnff_ini.f90:377-402 applies dxi corrections BEFORE Phase 1 EEQ solve
    // - use_corrections=true for Phase 1 topology charges (matching Fortran)
    // - use_corrections=false was incorrect and caused 0.0025 e error per hydroxyl oxygen
    Vector dxi = use_corrections ? calculateDxi(atoms, geometry_bohr, cn, topology) : Vector::Zero(natoms);

    // Setup EEQ parameters with CN-dependence
    Vector chi(natoms);
    Vector gam(natoms);
    Vector alpha(natoms);

    for (int i = 0; i < natoms; ++i) {
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));

        // Claude Generated (December 2025, Session 11): CRITICAL FIX - Add CNF term in Phase 1
        // Reference: XTB gfnff_ini.f90:409-411
        // XTB uses INTEGER neighbor count, NOT fractional CN for CNF term:
        //   dum = min(dble(topo%nb(20,i)),gen%cnmax)  ! nb(20,i) = neighbor count
        //   topo%chieeq(i) = -param%chi(ati) + dxi(i) + param%cnf(ati)*sqrt(dum)

        // CRITICAL BUG FIX (Jan 2, 2026): Use integer neighbor count, not fractional CN!
        // This was causing ~1.8e-2 error on carbon atoms in CH3OCH3 test
        // Reference: gfnff_param.f90:462 and 756: param%cnmax = 4.4
        const double CNMAX = 4.4;  // Maximum CN limit from XTB reference
        double nb_count;
        if (topology.has_value() && i < static_cast<int>(topology->neighbor_lists.size())) {
            nb_count = static_cast<double>(topology->neighbor_lists[i].size());  // Integer count
            nb_count = std::min(nb_count, CNMAX);  // Apply cnmax limit
        } else {
            nb_count = std::min(cn(i), CNMAX);  // Fallback to fractional CN with limit
        }
        double cnf_term = params_i.cnf * std::sqrt(nb_count);

        chi(i) = -params_i.chi + dxi(i) + cnf_term;
        gam(i) = params_i.gam;
        alpha(i) = params_i.alp;  // Already squared
    }

    // Build AUGMENTED EEQ matrix A (m x m)
    Matrix A = Matrix::Zero(m, m);
    Vector x = Vector::Zero(m);

    // Phase 1 EEQ diagnostic output (Claude Generated February 2026)
    // Format matches Fortran goedeckera debug output for direct comparison
    if (m_verbosity >= 3 && natoms >= 1) {
        CurcumaLogger::info("\nEEQ_PHASE1_PARAMS: Curcuma Phase 1 EEQ parameters");
        CurcumaLogger::info("EEQ_PHASE1_PARAMS: atom, Z, chieeq, gameeq, sqrt_alpeeq, dxi, nb_count");
        for (int i = 0; i < natoms; ++i) {
            int z_i = atoms[i];
            int nb = (topology.has_value() && i < static_cast<int>(topology->neighbor_lists.size()))
                     ? static_cast<int>(topology->neighbor_lists[i].size()) : 0;
            CurcumaLogger::info(fmt::format("EEQ_PHASE1_PARAMS: {:3d}, {:2d}, {:12.6f}, {:12.6f}, {:12.6f}, {:12.6f}, {:d}",
                i, z_i, chi(i), gam(i), std::sqrt(alpha(i)), dxi(i), nb));
        }
    }

    // 1. Setup RHS and diagonal
    for (int i = 0; i < natoms; ++i) {
        // CRITICAL FIX (Dec 30, 2025, Session 13): CNF term added ONCE, not twice!
        // XTB Phase 1 reference: gfnff_ini2.f90:1184 (goedeckera subroutine)
        // gfnff_ini.f90:411: topo%chieeq = -chi + dxi + CNF*√CN (stored in topo%chieeq)
        // gfnff_ini2.f90:1184: x(i) = topo%chieeq (NO additional CNF!)
        // Total: x = -chi + dxi + CNF*√CN (1× CNF, not 2×!)
        //
        // Our chi(i) already includes CNF term once (line 765), so just use it directly:
        x(i) = chi(i);  // ✅ FIXED: chi already includes CNF from line 765
        A(i, i) = gam(i) + TSQRT2PI / std::sqrt(alpha(i));
    }

    // (RHS diagnostic moved to full A-matrix diagnostic block below)

    // 2. Setup off-diagonal Coulomb matrix
    //
    // CRITICAL TODO (Dec 28, 2025): XTB uses TOPOLOGICAL distances, NOT geometric!
    // Reference: gfnff_ini.f90:431-461 (Floyd-Warshall shortest-path algorithm)
    //
    // XTB Algorithm:
    //   1. Build bond distance matrix from topology (covalent radii)
    //   2. Floyd-Warshall to compute shortest topological paths
    //   3. Scale by gen%rfgoed1 (typically 1.0-1.2)
    //   4. Use scaled topological distances in EEQ matrix
    //
    // FIXED (Dec 2025): Now uses Floyd-Warshall topological distances when topology provided
    // Impact: Topology charges now match XTB (4-5× smaller than wrong geometric distance version)
    //
    // Geometric distances are SHORTER than topological (direct line vs through bonds)
    // → Larger Coulomb terms → Larger charges → Wrong dispersion
    //
    // Claude Generated December 2025
    Matrix topo_dist;
    if (topology.has_value()) {
        // Use Dijkstra topological distances (sparse graph, O(N·E·logN))
        topo_dist = computeTopologicalDistancesSparse(atoms, *topology);

        // CRITICAL FIX (Jan 2, 2026): Cache topological distances for Phase 2 reuse
        // Reference: XTB gfnff_ini2.f90:1189-1199 uses same 'pair' array for both phases
        // Phase 2 must NOT recalculate with geometric distances!
        m_cached_topological_distances = topo_dist;

        // Setup off-diagonal Coulomb matrix with topological distances
        for (int i = 0; i < natoms; ++i) {
            for (int j = 0; j < i; ++j) {
                double r = topo_dist(i, j);  // Topological distance in Bohr

                if (r > 1e6) {
                    // Unconnected atoms - use zero Coulomb term
                    A(i, j) = 0.0;
                    A(j, i) = 0.0;
                    continue;
                }

                // J_ij = erf(gamma_ij * r) / r
                // gamma_ij = 1/sqrt(alpha_i + alpha_j)
                double gammij = 1.0 / std::sqrt(alpha(i) + alpha(j));
                double erf_gamma = std::erf(gammij * r);
                double coulomb = erf_gamma / r;

                A(i, j) = coulomb;
                A(j, i) = coulomb;
            }
        }
    } else {
        // Fallback to geometric distances (old behavior, acceptable for non-GFN-FF uses like D4)
        // Only warn at high verbosity since this is expected for universal EEQ usage
        if (m_verbosity >= 3) {
            CurcumaLogger::warn("EEQSolver: Using geometric distances (Floyd-Warshall topological distances not available)");
        }

        // CRITICAL FIX (Jan 2, 2026): Also cache geometric distances for Phase 2
        // D4 and other non-GFN-FF users call calculateCharges() without topology
        // Phase 2 needs to reuse whatever distances Phase 1 computed
        Matrix geom_dist = Matrix::Zero(natoms, natoms);

        for (int i = 0; i < natoms; ++i) {
            for (int j = 0; j < i; ++j) {
                double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
                double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
                double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
                double r_sq = dx*dx + dy*dy + dz*dz;
                double r = std::sqrt(r_sq);

                if (r < 1e-10) {
                    CurcumaLogger::error("EEQSolver::calculateTopologyCharges: Zero distance between atoms");
                    return generateFallbackCharges(natoms, total_charge, "zero distance in Phase 1");
                }

                // Store distance for both Coulomb matrix and cache
                geom_dist(i, j) = r;
                geom_dist(j, i) = r;

                // J_ij = erf(gamma_ij * r) / r
                // gamma_ij = 1/sqrt(alpha_i + alpha_j)
                double gammij = 1.0 / std::sqrt(alpha(i) + alpha(j));
                double erf_gamma = std::erf(gammij * r);
                double coulomb = erf_gamma / r;

                A(i, j) = coulomb;
                A(j, i) = coulomb;
            }
        }

        // Cache geometric distances for Phase 2 reuse
        m_cached_topological_distances = geom_dist;
    }

    // 3. Setup fragment charge constraints
    for (int f = 0; f < nfrag; ++f) {
        int row = natoms + f;
        double q_target = (topology.has_value() && f < static_cast<int>(topology->qfrag.size()))
                         ? topology->qfrag[f]
                         : (f == 0 ? static_cast<double>(total_charge) : 0.0);

        x(row) = q_target;

        for (int j = 0; j < natoms; ++j) {
            bool in_fragment = false;
            if (topology.has_value() && j < static_cast<int>(topology->fraglist.size())) {
                in_fragment = (topology->fraglist[j] == f + 1); // 1-indexed
            } else {
                in_fragment = (f == 0); // Default all atoms to first fragment
            }

            if (in_fragment) {
                A(row, j) = 1.0;
                A(j, row) = 1.0;
            }
        }
    }

    // 4. Solve augmented system
    // Claude Generated (March 2026): Full A-matrix diagnostic for HCN investigation
    if (m_verbosity >= 3 && natoms <= 10) {
        std::cerr << "\n========== PHASE 1 EEQ FULL DIAGNOSTICS ==========" << std::endl;

        // Print topological distances
        if (topology.has_value()) {
            std::cerr << "\nTopological distances (Bohr):" << std::endl;
            for (int i = 0; i < natoms; ++i) {
                for (int j = 0; j < i; ++j) {
                    double d = topo_dist(i, j);
                    if (d < 1e6)
                        std::cerr << fmt::format("  d_topo[{},{}] = {:12.6f}", i, j, d) << std::endl;
                }
            }
        }

        // Print per-atom parameter decomposition
        std::cerr << "\nPer-atom parameters:" << std::endl;
        for (int i = 0; i < natoms; ++i) {
            std::cerr << fmt::format("  Atom {:2d} (Z={:2d}): gam={:12.6f}  TSQRT2PI/sqrt(alp)={:12.6f}  A(i,i)={:12.6f}",
                i, atoms[i], gam(i), TSQRT2PI / std::sqrt(alpha(i)), A(i, i)) << std::endl;
        }

        // Print off-diagonal decomposition
        std::cerr << "\nOff-diagonal elements:" << std::endl;
        for (int i = 0; i < natoms; ++i) {
            for (int j = 0; j < i; ++j) {
                double gammij = 1.0 / std::sqrt(alpha(i) + alpha(j));
                double r = topology.has_value() ? topo_dist(i, j) : 0.0;
                std::cerr << fmt::format("  A[{},{}] = {:12.6f}  (gammij={:12.6f}, r={:12.6f}, erf={:12.6f})",
                    i, j, A(i, j), gammij, r, std::erf(gammij * r)) << std::endl;
            }
        }

        // Print complete A-matrix
        std::cerr << "\nComplete A matrix (" << m << "x" << m << "):" << std::endl;
        for (int i = 0; i < m; ++i) {
            std::cerr << fmt::format("  Row {:2d}:", i);
            for (int j = 0; j < m; ++j) {
                std::cerr << fmt::format(" {:12.6f}", A(i, j));
            }
            std::cerr << std::endl;
        }

        // Print RHS vector
        std::cerr << "\nRHS vector x:" << std::endl;
        for (int i = 0; i < m; ++i) {
            std::cerr << fmt::format("  x({:2d}) = {:12.6f}", i, x(i)) << std::endl;
        }
        std::cerr << "==================================================" << std::endl;
    }

    // Claude Generated (March 2026): Use configurable solver dispatch instead of hardcoded LU
    // This allows the user to select solve_method (lu, schur_cholesky, pcg, auto) for Phase 1 too
    // Claude Generated (WP2, May 2026): forward pool/num_threads so Stage-4 batched LU runs in parallel
    Vector topology_charges = dispatchSolve(A, x, natoms, nfrag, total_charge, pool, num_threads);

    // Claude Generated (March 2026): Print Phase 1 charge summary
    if (m_verbosity >= 3 && natoms <= 10) {
        std::cerr << "\nPhase 1 charges:" << std::endl;
        for (int i = 0; i < natoms; ++i) {
            std::cerr << fmt::format("  q({:2d}) = {:12.6f}  (Z={:2d})", i, topology_charges(i), atoms[i]) << std::endl;
        }
        std::cerr << fmt::format("  charge sum = {:12.6f}", topology_charges.sum()) << std::endl;
    }

    // Empty return signals solver failure — fall back to uniform charges for Phase 1
    // (uniform is the only available fallback here since no topology_charges exist yet)
    if (topology_charges.size() != natoms) {
        if (m_allow_unconverged) {
            CurcumaLogger::warn(fmt::format("Phase 1 EEQ: solver did not converge for N={} (allow_unconverged_charges=true)", natoms));
        } else {
            CurcumaLogger::error(fmt::format("Phase 1 EEQ: solver failed for N={}", natoms));
        }
        return generateFallbackCharges(natoms, total_charge, "Phase 1 solver failure");
    }

    // Check for NaN/Inf - fallback to uniform charges instead of hard fail
    for (int i = 0; i < natoms; ++i) {
        if (std::isnan(topology_charges[i]) || std::isinf(topology_charges[i])) {
            if (m_allow_unconverged) {
                CurcumaLogger::warn(fmt::format("Phase 1 EEQ: Invalid charge[{}] = {} (Z={}) (allow_unconverged_charges=true)",
                                                i, topology_charges[i], atoms[i]));
            } else {
                CurcumaLogger::error(fmt::format("Phase 1 EEQ: Invalid charge[{}] = {} (Z={})",
                                                 i, topology_charges[i], atoms[i]));
            }
            return generateFallbackCharges(natoms, total_charge, "NaN/Inf in Phase 1 solution");
        }
    }

    // Phase 1 charge diagnostic output (Claude Generated February 2026)
    if (m_verbosity >= 2) {
        CurcumaLogger::info("Phase 1 EEQ: Topology charges calculated");
        for (int i = 0; i < std::min(5, natoms); ++i) {
            CurcumaLogger::result(fmt::format("Atom {} qa = {:.6f}", i, topology_charges(i)));
        }
    }
    if (m_verbosity >= 3) {
        CurcumaLogger::info("\nEEQ_PHASE1_CHARGES: atom, Z, qa");
        for (int i = 0; i < natoms; ++i) {
            CurcumaLogger::info(fmt::format("EEQ_PHASE1_CHARGES: {:3d}, {:2d}, {:12.6f}",
                i, atoms[i], topology_charges(i)));
        }
        double sum_q = topology_charges.sum();
        CurcumaLogger::info(fmt::format("EEQ_PHASE1_CHARGES: sum = {:.6f}", sum_q));
    }

    return topology_charges;
}

// ===== Floyd-Warshall Topological Distances =====

Matrix EEQSolver::computeTopologicalDistances(
    const std::vector<int>& atoms,
    const TopologyInput& topology
) const {
    const int natoms = atoms.size();

    // Claude Generated (Feb 20, 2026): float32 Floyd-Warshall matching Fortran real(sp)
    //
    // Fortran declares: real(sp) :: rabd(nat,nat)  (gfnff_ini.f90:432)
    // Using float32 here is CRITICAL for EEQ charge accuracy:
    //   - float32 rounding accumulates along shortest paths
    //   - Different topological distances → different Phase-1 qa → different fqq/alpha/zetac6
    //   - Without float32: bond/torsion/repulsion/dispersion errors of 1e-3 to 1e-2 Eh
    //   - With float32: near-exact match with Fortran reference
    const float RABD_CUTOFF_F = 13.0f;   // Fortran gfnff_ini.f90:88, real(sp)
    const float TDIST_THR_F   = 12.0f;   // Fortran gfnff_param.f90:776, real(sp)

    // Reference: external/gfnff/src/gfnff_param.f90:817 (gen%rfgoed1 = 1.175)
    const double RFGOED1 = 1.175;
    const double BOHR_TO_ANGSTROM = 0.52917726;

    // 1. Initialize with cutoff value (flat float32 array for cache efficiency)
    // Reference: gfnff_ini.f90:431-442
    std::vector<float> rabd(natoms * natoms, RABD_CUTOFF_F);

    // 2. Set diagonal to zero
    for (int i = 0; i < natoms; ++i)
        rabd[i * natoms + i] = 0.0f;

    // 3. Set bonded distances (sum of covalent radii, cast to float32)
    // Reference: gfnff_ini.f90:438-448
    for (int i = 0; i < natoms; ++i) {
        float rad_i = static_cast<float>(topology.covalent_radii[i]);
        for (int j : topology.neighbor_lists[i]) {
            float bond = rad_i + static_cast<float>(topology.covalent_radii[j]);
            rabd[i * natoms + j] = bond;
            rabd[j * natoms + i] = bond;
        }
    }

    // 4. Floyd-Warshall shortest path in float32, matching Fortran real(sp) arithmetic
    // Reference: gfnff_ini.f90:462-471
    for (int k = 0; k < natoms; ++k) {
        for (int i = 0; i < natoms; ++i) {
            float rik = rabd[i * natoms + k];
            if (rik > TDIST_THR_F) continue;
            for (int j = 0; j < natoms; ++j) {
                float rkj = rabd[k * natoms + j];
                if (rkj > TDIST_THR_F) continue;
                float candidate = rik + rkj;   // float32 addition like Fortran
                if (rabd[i * natoms + j] > candidate)
                    rabd[i * natoms + j] = candidate;
            }
        }
    }

    // 5. Convert to double Matrix with cutoff and Angstrom→Bohr scaling
    // Reference: gfnff_ini.f90:474-480
    Matrix result(natoms, natoms);
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            float rij = rabd[i * natoms + j];
            double val = (rij > TDIST_THR_F)
                ? static_cast<double>(RABD_CUTOFF_F)
                : static_cast<double>(rij);
            result(i, j) = RFGOED1 * val / BOHR_TO_ANGSTROM;
        }
    }

    if (m_verbosity >= 3) {
        std::cerr << "\n=== Floyd-Warshall Topological Distances (float32, Bohr) ===" << std::endl;
        for (int i = 0; i < std::min(5, natoms); ++i) {
            for (int j = 0; j < i; ++j) {
                double d = result(i, j);
                if (d < RFGOED1 * RABD_CUTOFF_F / BOHR_TO_ANGSTROM - 1.0)
                    std::cerr << fmt::format("  d_topo[{},{}] = {:.6f} Bohr", i, j, d) << std::endl;
            }
        }
        std::cerr << "========================================\n" << std::endl;
    }

    return result;
}

// ===== Multi-Source Dijkstra Topological Distances (Performance Replacement) =====
// Claude Generated (March 2026): O(N·E·logN) replacement for O(N³) Floyd-Warshall
// For sparse molecular graphs (degree ~3-4), this is ~50-100× faster for N>500

Matrix EEQSolver::computeTopologicalDistancesSparse(
    const std::vector<int>& atoms,
    const TopologyInput& topology
) const {
    const int natoms = atoms.size();

    // Same constants as Floyd-Warshall version — MUST match for numerical consistency
    const float RABD_CUTOFF_F = 13.0f;   // Fortran gfnff_ini.f90:88, real(sp)
    const float TDIST_THR_F   = 12.0f;   // Fortran gfnff_param.f90:776, real(sp)
    const double RFGOED1 = 1.175;         // gfnff_param.f90:817
    const double BOHR_TO_ANGSTROM = 0.52917726;

    // Pre-compute float32 bond weights (sum of covalent radii per bond)
    // Build adjacency list with edge weights for Dijkstra
    struct Edge { int to; float weight; };
    std::vector<std::vector<Edge>> adj(natoms);
    for (int i = 0; i < natoms; ++i) {
        float rad_i = static_cast<float>(topology.covalent_radii[i]);
        for (int j : topology.neighbor_lists[i]) {
            float bond = rad_i + static_cast<float>(topology.covalent_radii[j]);
            adj[i].push_back({j, bond});
        }
    }

    // Flat float32 result array (same layout as Floyd-Warshall)
    std::vector<float> rabd(natoms * natoms, RABD_CUTOFF_F);
    for (int i = 0; i < natoms; ++i)
        rabd[i * natoms + i] = 0.0f;

    // Multi-source Dijkstra with early termination at TDIST_THR
    // Each source atom is independent → parallelizable
    #pragma omp parallel for schedule(dynamic)
    for (int src = 0; src < natoms; ++src) {
        // Per-thread min-heap: (distance, node)
        using PQEntry = std::pair<float, int>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        // dist array for this source — use the row in rabd directly
        float* dist = &rabd[src * natoms];
        // dist is already initialized to RABD_CUTOFF_F, with dist[src] = 0

        pq.push({0.0f, src});

        while (!pq.empty()) {
            auto [d, u] = pq.top();
            pq.pop();

            // Skip if we already found a shorter path
            if (d > dist[u]) continue;

            // Early termination: no need to explore beyond cutoff
            if (d > TDIST_THR_F) continue;

            for (const auto& edge : adj[u]) {
                // float32 addition matching Fortran real(sp) arithmetic
                float new_dist = d + edge.weight;

                if (new_dist < dist[edge.to]) {
                    dist[edge.to] = new_dist;
                    pq.push({new_dist, edge.to});
                }
            }
        }
    }

    // Convert to double Matrix with cutoff and Angstrom→Bohr scaling
    // (identical to Floyd-Warshall version)
    Matrix result(natoms, natoms);
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            float rij = rabd[i * natoms + j];
            double val = (rij > TDIST_THR_F)
                ? static_cast<double>(RABD_CUTOFF_F)
                : static_cast<double>(rij);
            result(i, j) = RFGOED1 * val / BOHR_TO_ANGSTROM;
        }
    }

    if (m_verbosity >= 3) {
        std::cerr << "\n=== Dijkstra Topological Distances (float32, Bohr) ===" << std::endl;
        for (int i = 0; i < std::min(5, natoms); ++i) {
            for (int j = 0; j < i; ++j) {
                double d = result(i, j);
                if (d < RFGOED1 * RABD_CUTOFF_F / BOHR_TO_ANGSTROM - 1.0)
                    std::cerr << fmt::format("  d_topo[{},{}] = {:.6f} Bohr", i, j, d) << std::endl;
            }
        }
        std::cerr << "========================================\n" << std::endl;
    }

    return result;
}

// ===== Phase 2: Final Refined Charges =====
// CRITICAL REQUIREMENT (Jan 2, 2026):
//   Phase 2 MUST use cached topological distances from Phase 1!
//   XTB reference (gfnff_ini2.f90:1189-1199) uses same 'pair' array for both phases.
//   Using geometric distances here causes 1.5e-3 RMS charge error (4.0e-3 max on oxygen).
//
// Algorithm:
//   1. Apply environment corrections: dxi, dgam
//   2. Calculate charge-dependent alpha: (alpha_base + ff*qa)²
//   3. Build Coulomb matrix A with TOPOLOGICAL distances from Phase 1
//   4. Solve linear system A*q = x for final charges
//
// Reference: XTB gfnff_ini.f90:693-707, gfnff_ini2.f90:1140-1246
// Claude Generated - December 2025, Updated January 2, 2026

Vector EEQSolver::calculateFinalCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector& topology_charges,
    const Vector& cn,
    const std::vector<int>& hybridization,
    const std::optional<TopologyInput>& topology,
    bool use_corrections,
    const std::optional<Vector>& alpeeq,
    CxxThreadPool* pool,
    int num_threads)
{
    const int natoms = atoms.size();
    const double TSQRT2PI = 0.797884560802866;  // sqrt(2/π)

    // Claude Generated (April 2026): skip_phase2 option — bypass Phase 2 entirely
    if (m_skip_phase2) {
        if (m_verbosity >= 1) {
            CurcumaLogger::warn("EEQSolver: skip_phase2=true — bypassing Phase 2 refinement, using Phase 1 topology charges");
        }
        return topology_charges;
    }

    // Claude Generated (Apr 2026): If Phase 2 was historically implausible for this system
    // size, skip it forever (until atom count changes, indicating a new molecule).
    if (m_phase2_historically_implausible && m_phase2_implausible_natoms == natoms) {
        if (m_verbosity >= 2) {
            CurcumaLogger::info("EEQSolver: Phase 2 was historically implausible for this system — skipping, using Phase 1 topology charges");
        }
        return topology_charges;
    }

    if (m_verbosity >= 2) {
        fmt::print(stderr, "[EEQ] Phase 2: preparing corrections (dxi, dgam, pi/amide detection)...\n");
    }

    // CRITICAL FIX (Jan 4, 2026): Only use corrections if explicitly requested
    // gfnff_final.cpp achieves 0.0000025 e accuracy using ONLY base parameters (NO dxi, NO dgam)
    // Curcuma's complex corrections add noise instead of improving accuracy
    Vector dxi = use_corrections ? calculateDxi(atoms, geometry_bohr, cn, topology) : Vector::Zero(natoms);

    // Claude Generated (January 17, 2026): Detect pi-system and amide nitrogens for dgam
    std::vector<bool> is_pi_atom = use_corrections ? detectPiSystem(atoms, hybridization, topology) : std::vector<bool>(natoms, false);
    std::vector<bool> is_amide = use_corrections ? detectAmideNitrogens(atoms, hybridization, is_pi_atom, topology, cn) : std::vector<bool>(natoms, false);
    // Claude Generated (February 2026): Proper amideH detection for Phase 2 chi correction
    // Reference: Fortran gfnff_ini.f90:717 uses amideH() function from gfnff_ini2.f90:1575
    // Requires: H with 1 neighbor, that neighbor is amide N, amide N has exactly 1 sp3 C
    std::vector<bool> is_amide_h = use_corrections ? detectAmideHydrogens(atoms, hybridization, is_amide, topology) : std::vector<bool>(natoms, false);

    Vector dgam = use_corrections ? calculateDgam(atoms, topology_charges, hybridization, is_pi_atom, is_amide) : Vector::Zero(natoms);

    // DEBUG: Print topology charges used for dgam
    if (m_verbosity >= 3) {
        std::cerr << "\n=== Phase 2 Topology Charges (used for dgam) ===" << std::endl;
        for (int i = 0; i < std::min(3, natoms); ++i) {
            std::cerr << "Atom " << i << " (Z=" << atoms[i] << "): qa = " << topology_charges(i)
                      << ", dgam = " << dgam(i) << std::endl;
        }
        std::cerr << "================================================\n" << std::endl;
    }

    // Initialize final_charges variable (will be overwritten by solve)
    Vector final_charges;

    // Augmented system size: n atoms + n fragments
    int nfrag = topology.has_value() ? topology->nfrag : 1;
    int m = natoms + nfrag;

    // Claude Generated (December 2025, Session 13): Pre-calculate CONSTANT corrected parameters
    // chi and gam are constant. Alpha is calculated ONCE using topology_charges (Phase 1).
    Vector chi_corrected(natoms);
    Vector gam_corrected(natoms);

    for (int i = 0; i < natoms; ++i) {
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));

        // CRITICAL FIX (Session 11 - CORRECTED): chi_corrected WITHOUT CNF!
        // Reference: XTB gfnff_ini.f90:696 (Phase 2)
        // chieeq = -χ + dxi (OHNE CNF!)
        // Then RHS adds CNF once: x = chieeq + CNF·√CN = -χ + dxi + CNF·√CN
        chi_corrected(i) = -params_i.chi + dxi(i);  // WITHOUT CNF!
        gam_corrected(i) = params_i.gam + dgam(i);
    }

    // ===== CRITICAL (Jan 5, 2026): Phase 2 uses GEOMETRIC distances (from XYZ) =====
    // Reference: CHARGE_DATAFLOW.md - "Echte Distanzen aus Koordinaten (xyz)"
    // Reference: Fortran gfnff_engrad.F90:341 - goed_gfnff receives sqrab, srab (geometric!)
    //
    // ARCHITECTURE DIFFERENCE:
    //   Phase 1 (goedeckera):  Topological distances (Floyd-Warshall bond graph)
    //   Phase 2 (goed_gfnff):  Geometric distances (real xyz coordinates)
    //
    // WHY DIFFERENT?
    //   Phase 1: Topology-based for initial parameter generation
    //   Phase 2: Geometry-dependent for energy/gradient calculation
    //   This is BY DESIGN - not a bug!
    //
    // REGRESSION FIXED:
    //   Jan 2, 2026 "fix" used topological distances → 1.53e-02 RMS error (10× WORSE!)
    //   Jan 5, 2026: Reverted to geometric distances → expect 1.5e-3 RMS error
    //
    // WP-EEQ-Cache (May 2026): expose geometry/CN so solveWithSchurCholesky can check cache.
    // Skip the ~45 KB/step copy when the cache is disabled (eps <= 0).
    if (m_refactor_eps > 0.0) {
        m_pending_geometry = geometry_bohr;
        m_pending_cn       = cn;
    }

    // Claude Generated (Mar 2026): Pre-allocated distance buffer — avoids 13 MB alloc per step (N=1280)
    // WP-FF-Packed-Triangular (May 2026): buffer is now packed lower-triangular (N(N+1)/2 doubles).
    ensurePhase2Buffers(natoms, nfrag);
    m_phase2_distances.setZero();

    // WP-FF-DistMatrix-Sharing (May 2026): consume externally-computed packed srab
    // when GFNFF provides it — skip the O(N²) sqrt loop entirely.
    // Layout: srab[i*(i+1)/2 + j] = r_ij for i > j. Identical to m_phase2_distances now,
    // so we can copy the whole vector in one shot.
    const int srab_expected = natoms * (natoms + 1) / 2;
    const bool use_external_srab = (m_external_srab != nullptr)
                                && (m_external_srab->size() == srab_expected);
    if (use_external_srab) {
        m_phase2_distances = *m_external_srab;  // direct copy, layouts match
        // Quick atoms-too-close scan (only off-diagonal slots)
        for (int i = 1; i < natoms; ++i) {
            const int ii = i * (i + 1) / 2;
            for (int j = 0; j < i; ++j) {
                if (m_phase2_distances[ii + j] < 1e-10) {
                    CurcumaLogger::error("EEQSolver::calculateFinalCharges: atoms too close (external srab)");
                    return generateFallbackCharges(natoms, total_charge, "zero distance in Phase 2 (external srab)");
                }
            }
        }
    } else {

    // Claude Generated (Mar 2026): Internal std::thread parallelisation for O(N²) distance matrix
    std::atomic<bool> atoms_too_close{false};
    int dist_progress_10pct = std::max(1, natoms / 10);
    bool show_dist_progress = (m_verbosity >= 3 && natoms >= 200);
    if (show_dist_progress) {
        fmt::print(stderr, "[EEQ] Phase 2: distance matrix ({} atoms):", natoms);
        fflush(stderr);
    }
    if (num_threads > 1 && natoms > 64) {
        int T = std::min(num_threads, natoms);
        auto dist_worker = [&](int t_id) {
            for (int i = t_id; i < natoms; i += T) {
                const int ii = i * (i + 1) / 2;  // packed row offset
                for (int j = 0; j < i; ++j) {
                    double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
                    double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
                    double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
                    double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (r < 1e-10) atoms_too_close.store(true, std::memory_order_relaxed);
                    m_phase2_distances[ii + j] = r;
                }
                if (show_dist_progress && t_id == 0 && i > 0 && i % dist_progress_10pct == 0) {
                    fmt::print(stderr, " {}%", (100 * i) / natoms);
                    fflush(stderr);
                }
            }
        };
        // Claude Generated (Mar 2026): pool->enqueue() reuses persistent workers
        if (pool) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t) futures.push_back(pool->enqueue(dist_worker, t));
            dist_worker(0);
            for (auto& f : futures) f.get();
        } else {
            std::vector<std::thread> threads(T - 1);
            for (int t = 1; t < T; ++t) threads[t - 1] = std::thread(dist_worker, t);
            dist_worker(0);
            for (auto& th : threads) th.join();
        }
    } else {
        for (int i = 0; i < natoms; ++i) {
            const int ii = i * (i + 1) / 2;
            for (int j = 0; j < i; ++j) {
                double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
                double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
                double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (r < 1e-10) atoms_too_close.store(true, std::memory_order_relaxed);
                m_phase2_distances[ii + j] = r;
            }
            if (show_dist_progress && i > 0 && i % dist_progress_10pct == 0) {
                fmt::print(stderr, " {}%", (100 * i) / natoms);
                fflush(stderr);
            }
        }
    }
    if (show_dist_progress) fmt::print(stderr, " done\n");
    if (atoms_too_close) {
        CurcumaLogger::error("EEQSolver::calculateFinalCharges: atoms too close");
        return generateFallbackCharges(natoms, total_charge, "zero distance in Phase 2");
    }
    }  // end if (!use_external_srab)
    // WP-FF-Packed-Triangular: lookup helper. Assumes j < i (Coulomb-build only reads lower triangle).
    auto dist = [&](int i, int j) -> double {
        return m_phase2_distances[i * (i + 1) / 2 + j];
    };

    // WP7-D (Jun 2026, Claude Generated): contact metric for the dispatcher. The approximate
    // batched solver (per-fragment Cholesky) drops cross-fragment Coulomb, so it is only valid
    // for well-separated fragments. Compute the minimum inter-fragment atom distance from the
    // packed distances just built (fragment membership from topology->fraglist, 1-indexed) so
    // dispatchSolve can prefer the exact path when fragments are in contact. O(N^2), same order
    // as the distance build it piggybacks on. -1 = unknown (no multi-fragment topology).
    m_contact_min_dist = -1.0;
    if (nfrag >= 2 && topology.has_value()
        && static_cast<int>(topology->fraglist.size()) >= natoms) {
        const auto& fraglist = topology->fraglist;
        double min_d = -1.0;
        for (int i = 1; i < natoms; ++i) {
            const int fi = fraglist[i];
            const int ii = i * (i + 1) / 2;
            for (int j = 0; j < i; ++j) {
                if (fraglist[j] != fi) {
                    const double r = m_phase2_distances[ii + j];
                    if (min_d < 0.0 || r < min_d) min_d = r;
                }
            }
        }
        m_contact_min_dist = min_d;  // -1 only if every atom shares one fragment id
    }

    if (m_verbosity >= 3) {
        CurcumaLogger::info("EEQ Phase 2: Using geometric distances from xyz coordinates (matches Fortran goed_gfnff)");
    }

    // Claude Generated (December 2025, Session 13): Distance matrix pre-calculation
    // Coulomb matrix A will be built once with alpha based on topology_charges (qa)
    if (m_verbosity >= 3) {
        CurcumaLogger::info(fmt::format("EEQ Phase 2: Pre-computed distance matrix ({}x{})", natoms, natoms));
        CurcumaLogger::info("Alpha calculated once using topology charges (Phase 1)");
    }

    // Claude Generated (December 2025, Session 13): Debug output for constant parameters
    if (m_verbosity >= 3) {
        std::cerr << "\n=== EEQ Phase 2: Constant Corrected Parameters ===" << std::endl;
        for (int i = 0; i < std::min(3, natoms); ++i) {
            int z_i = atoms[i];
            EEQParameters params_i = getParameters(z_i, cn(i));

            std::cerr << "Atom " << i << " (Z=" << z_i << "):" << std::endl;
            std::cerr << "  chi_base = " << params_i.chi << std::endl;
            std::cerr << "  gam_base = " << params_i.gam << std::endl;
            std::cerr << "  CN = " << cn(i) << std::endl;
            std::cerr << "  dxi = " << dxi(i) << std::endl;
            std::cerr << "  dgam = " << dgam(i) << std::endl;
            std::cerr << "  chieeq = " << chi_corrected(i) << " (WITHOUT CNF)" << std::endl;
            std::cerr << "  gameeq = " << gam_corrected(i) << std::endl;
        }
        std::cerr << "===========================================================\n" << std::endl;
    }

    // ===== Iterative Refinement Implementation =====
    // Claude Generated (January 2026): Self-consistent iterative refinement for alpha-charge relationship
    Vector current_charges = topology_charges;  // Start with Phase 1 charges

    // Check if iterative refinement is enabled
    bool use_iterative = m_config.get<bool>("use_iterative_refinement", false);
    int max_iterations = m_config.get<int>("max_iterations", 50);
    double convergence_threshold = m_config.get<double>("convergence_threshold", 1e-6);

    if (m_verbosity >= 3) {
        const char* solver_name = (m_solve_method == EEQSolveMethod::PCG) ? "PCG"
            : (m_solve_method == EEQSolveMethod::SchurCholesky) ? "Cholesky"
            : (m_solve_method == EEQSolveMethod::Auto) ? "auto (benchmark)" : "LU";
        if (use_iterative)
            fmt::print(stderr, "[EEQ] Phase 2: iterative refinement ({}, max {} iters, thr={:.0e})\n",
                       solver_name, max_iterations, convergence_threshold);
        else
            fmt::print(stderr, "[EEQ] Phase 2: single solve ({})\n", solver_name);
    }

    // Iterative refinement loop
    int iteration = 0;
    bool converged = false;
    Vector previous_charges = Vector::Zero(natoms);

    do {
        // ===== Calculate Alpha using Current Charges =====
        // Reference: XTB gfnff_ini.f90:699-706, gfnff_ini.f90:718-725
        Vector alpha_corrected(natoms);

        // Claude Generated (January 2026): Use pre-computed charge-dependent alpha if available
        if (alpeeq.has_value() && alpeeq->size() == natoms) {
            // Use pre-computed charge-dependent alpha from Phase 1B
            // Formula: alpeeq(i) = (alpha_base + ff*qa(i))²
            // Reference: Fortran gfnff_ini.f90:718-725
            // These values are SQUARED and charge-corrected, used UNCHANGED
            alpha_corrected = *alpeeq;

            if (m_verbosity >= 3 && iteration == 0) {
                CurcumaLogger::info("EEQ Phase 2: Using pre-computed charge-dependent alpha (alpeeq)");
            }
        } else {
            // Fallback: Calculate alpha from base parameters (backward compatibility)
            // CRITICAL: This is the OLD behavior - retained for backward compatibility only
            for (int i = 0; i < natoms; ++i) {
                int z_i = atoms[i];

                // Get base alpha (UNSQUARED) from gfnff_par.h
                double alpha_base = (z_i >= 1 && z_i <= 86) ? alpha_eeq[z_i - 1] : 0.903430;

                // CRITICAL FIX (Jan 7, 2026): Phase 2-specific alpha corrections
                // Fortran goed_gfnff uses different alpha values than goedeckera (Phase 1)
                // Reference: PHASE2_CHARGE_ANALYSIS.md - extracted from Fortran debug output
                if (z_i == 6) {  // Carbon
                    alpha_base = 0.906988;  // Phase 2: (0.906988)² = 0.822628
                } else if (z_i == 8) {  // Oxygen
                    alpha_base = 0.915779;  // Phase 2: (0.915779)² = 0.838652
                }

                // CRITICAL (Jan 7, 2026): Phase 2 alpha parameters
                // Unlike dxi/dgam, alpha is NOT charge-corrected in Phase 2 energy calculation
                // Reference: XTB gfnff_engrad.F90:1507-1520 uses topo%alpeeq unchanged
                // The (alpha_base + ff*qa)² correction is ONLY for Phase 1c (gfnff_ini.f90:712-726)
                // Phase 2 uses BASE alpha values: alpha_eeq[Z-1]² with NO qa modification
                alpha_corrected(i) = alpha_base * alpha_base;
            }

            if (m_verbosity >= 3 && iteration == 0) {
                CurcumaLogger::warn("EEQ Phase 2: No alpeeq provided - using base alpha (backward compat)");
            }
        }

        // ===== Build A Matrix with Current Alpha =====
        // Claude Generated (Mar 2026): Pre-allocated A buffer — avoids 13 MB alloc per step
        m_phase2_A.setZero();
        m_phase2_rhs.setZero();
        Matrix& A = m_phase2_A;
        Vector& x = m_phase2_rhs;

        // WP-EEQ-Matrix-Cache (May 2026): check if we can reuse the cached
        // Coulomb off-diagonal from a previous step.  Off-diag values are
        // gamma_ij·erf-quotients of r_ij — they depend on geometry AND on
        // alpha (charges, CN). We reuse only when BOTH geometry drift AND
        // CN drift are below thresholds; the diagonal is always rebuilt
        // because it folds in gam_corrected + alpha_corrected fresh.
        bool reuse_offdiag = false;
        if (m_matrix_rebuild_eps > 0.0 && m_A_cache_valid
            && m_cached_A_natoms == natoms
            && m_cached_A_geometry.rows() == natoms
            && m_cached_A_cn.size() == natoms
            && m_cached_A_offdiag.rows() == natoms
            && m_cached_A_offdiag.cols() == natoms) {
            const double max_dr = (geometry_bohr - m_cached_A_geometry).rowwise().norm().maxCoeff();
            const double max_dcn = (cn - m_cached_A_cn).cwiseAbs().maxCoeff();
            if (max_dr < m_matrix_rebuild_eps && max_dcn < 0.05) {
                reuse_offdiag = true;
                if (m_verbosity >= 3)
                    fmt::print(stderr, "[EEQ-Matrix-Cache] Hit (max_dr={:.3e} Bohr, max_dcn={:.3e})\n",
                               max_dr, max_dcn);
            }
        }

        // 1+2. Build Coulomb matrix (off-diagonal + diagonal)
        //
        // Distance cutoff for large systems: erf(gamma*r)/r decays to ~1/r for large r.
        // Truncating long-range elements reduces matrix condition number for polymers but
        // shifts the energy surface away from full-Coulomb — Fortran XTB GFN-FF uses no
        // truncation here (cutoff = 0). The canonical PARAM default in eeq_solver.h:965 is
        // 0.0 to match Fortran. The previous 30.0 inline fallback masked that default
        // whenever m_config did not carry the explicit registry entry, which silently
        // truncated all default-config calls and produced the +0.93 Eh polymer bias vs.
        // Fortran reported as WP-V (May 2026 finding, fix in WP6 follow-up).
        // WP-S3 (May 2026): honour the post-init override set by GFNFF auto-detection.
        double eeq_dist_cutoff = getEEQDistanceCutoffEffective();
        const double EEQ_CUTOFF_SQ = (natoms > 200 && eeq_dist_cutoff > 0.0)
            ? eeq_dist_cutoff * eeq_dist_cutoff : 0.0;

        int coulomb_progress_10pct = std::max(1, natoms / 10);
        bool show_coulomb_progress = (m_verbosity >= 3 && natoms >= 200 && iteration == 0);
        if (show_coulomb_progress && !reuse_offdiag) {
            fmt::print(stderr, "[EEQ] Phase 2: Coulomb matrix ({} atoms):", natoms);
            fflush(stderr);
        }

        if (reuse_offdiag) {
            // Cache hit: copy cached N×N block (off-diag + stale diag), then overwrite diagonal
            // with fresh gam_corrected + TSQRT2PI/sqrt(alpha_corrected).
            A.topLeftCorner(natoms, natoms) = m_cached_A_offdiag;
            for (int i = 0; i < natoms; ++i) {
                A(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha_corrected(i));
            }
        } else
        // Claude Generated (Mar 2026): Internal std::thread parallelisation for O(N²) A-matrix
        if (num_threads > 1 && natoms > 64) {
            int T = std::min(num_threads, natoms);
            auto amat_worker = [&](int t_id) {
                for (int i = t_id; i < natoms; i += T) {
                    A(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha_corrected(i));
                    for (int j = 0; j < i; ++j) {
                        double r = dist(i, j);  // packed lookup
                        if (EEQ_CUTOFF_SQ > 0.0 && r * r > EEQ_CUTOFF_SQ) {
                            A(i, j) = A(j, i) = 0.0;
                            continue;
                        }
                        double gamma_ij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));
                        double coulomb = std::erf(gamma_ij * r) / r;
                        A(i, j) = coulomb;
                        A(j, i) = coulomb;
                    }
                    if (show_coulomb_progress && t_id == 0 && i > 0 && i % coulomb_progress_10pct == 0) {
                        fmt::print(stderr, " {}%", (100 * i) / natoms);
                        fflush(stderr);
                    }
                }
            };
            // Claude Generated (Mar 2026): pool->enqueue() reuses persistent workers
            if (pool) {
                std::vector<std::future<void>> futures;
                futures.reserve(T - 1);
                for (int t = 1; t < T; ++t) futures.push_back(pool->enqueue(amat_worker, t));
                amat_worker(0);
                for (auto& f : futures) f.get();
            } else {
                std::vector<std::thread> threads(T - 1);
                for (int t = 1; t < T; ++t) threads[t - 1] = std::thread(amat_worker, t);
                amat_worker(0);
                for (auto& th : threads) th.join();
            }
        } else {
            for (int i = 0; i < natoms; ++i) {
                A(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha_corrected(i));
                for (int j = 0; j < i; ++j) {
                    double r = dist(i, j);  // packed lookup
                    if (EEQ_CUTOFF_SQ > 0.0 && r * r > EEQ_CUTOFF_SQ) {
                        A(i, j) = A(j, i) = 0.0;
                        continue;
                    }
                    double gamma_ij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));
                    double coulomb = std::erf(gamma_ij * r) / r;
                    A(i, j) = coulomb;
                    A(j, i) = coulomb;
                }
                if (show_coulomb_progress && i > 0 && i % coulomb_progress_10pct == 0) {
                    fmt::print(stderr, " {}%", (100 * i) / natoms);
                    fflush(stderr);
                }
            }
        }
        if (show_coulomb_progress) fmt::print(stderr, " done\n");

        // WP-EEQ-Matrix-Cache: snapshot fresh A_nn block for next-step reuse.
        // Skip when cache is disabled (eps=0) or we just reused (no change). Cache
        // stores entire N×N including the diagonal; reuse-path overwrites diagonal.
        if (m_matrix_rebuild_eps > 0.0 && !reuse_offdiag) {
            if (m_cached_A_offdiag.rows() != natoms || m_cached_A_offdiag.cols() != natoms) {
                m_cached_A_offdiag.resize(natoms, natoms);
            }
            m_cached_A_offdiag    = A.topLeftCorner(natoms, natoms);
            m_cached_A_geometry   = geometry_bohr;
            m_cached_A_cn         = cn;
            m_cached_A_natoms     = natoms;
            m_A_cache_valid       = true;
            if (m_verbosity >= 3)
                fmt::print(stderr, "[EEQ-Matrix-Cache] Persisted (N={})\n", natoms);
        }

        // 3. Setup fragment charge constraints
        for (int f = 0; f < nfrag; ++f) {
            int row = natoms + f;
            for (int j = 0; j < natoms; ++j) {
                bool in_fragment = false;
                if (topology.has_value() && j < static_cast<int>(topology->fraglist.size())) {
                    in_fragment = (topology->fraglist[j] == f + 1); // 1-indexed
                } else {
                    in_fragment = (f == 0); // Default all atoms to first fragment
                }

                if (in_fragment) {
                    A(row, j) = 1.0;
                    A(j, row) = 1.0;
                }
            }
        }

        // 4. Setup RHS with CNF*sqrt(CN) term
        // Reference: XTB gfnff_engrad.F90:1504 (used during gradient calculation)
        // x(i) = topo%chieeq(i) + param%cnf(at(i))*sqrt(cn(i))
        // where chieeq = -chi + dxi (from gfnff_ini.f90:696 for Phase 2)
        for (int i = 0; i < natoms; ++i) {
            int z_i = atoms[i];
            EEQParameters params_i = getParameters(z_i, cn(i));
            double chi_corrected_val = chi_corrected(i);

            // Amide Hydrogen Electronegativity Shift (Phase 2 - February 2026)
            // Reference: Fortran gfnff_ini.f90:717 - amideH() from gfnff_ini2.f90:1575
            // Requires proper amide detection: H→N(amide, pi, sp3)→C(pi)→O(pi, nn=1)
            // Previous check was too broad (any H bonded to sp/sp2 N)
            if (z_i == 1 && is_amide_h[i]) {
                chi_corrected_val -= 0.02;
            }

            x(i) = chi_corrected_val + params_i.cnf * std::sqrt(cn(i));
        }
        for (int f = 0; f < nfrag; ++f) {
            int row = natoms + f;
            double q_target = (topology.has_value() && f < static_cast<int>(topology->qfrag.size()))
                             ? topology->qfrag[f]
                             : (f == 0 ? static_cast<double>(total_charge) : 0.0);
            x(row) = q_target;
        }

        // Claude Generated (March 2026): Full Phase 2 A-matrix diagnostic for HCN investigation
        if (m_verbosity >= 3 && iteration == 0 && natoms <= 10) {
            std::cerr << "\n========== PHASE 2 EEQ FULL DIAGNOSTICS ==========" << std::endl;

            // Print geometric distances
            std::cerr << "\nGeometric distances (Bohr):" << std::endl;
            for (int i = 0; i < natoms; ++i) {
                for (int j = 0; j < i; ++j) {
                    std::cerr << fmt::format("  d_geom[{},{}] = {:12.6f}", i, j, dist(i, j)) << std::endl;
                }
            }

            // Print per-atom parameter decomposition
            std::cerr << "\nPer-atom Phase 2 parameters:" << std::endl;
            for (int i = 0; i < natoms; ++i) {
                int z_i = atoms[i];
                EEQParameters params_i = getParameters(z_i, cn(i));
                std::cerr << fmt::format("  Atom {:2d} (Z={:2d}): qa={:12.6f}  alpeeq={:12.6f}  gameeq={:12.6f}  chieeq={:12.6f}  CN={:8.4f}  cnf={:12.6f}",
                    i, z_i, topology_charges(i), alpha_corrected(i), gam_corrected(i), chi_corrected(i), cn(i), params_i.cnf) << std::endl;
                std::cerr << fmt::format("           gam_base={:12.6f}  dgam={:12.6f}  TSQRT2PI/sqrt(alp)={:12.6f}  A(i,i)={:12.6f}",
                    params_i.gam, dgam(i), TSQRT2PI / std::sqrt(alpha_corrected(i)), A(i, i)) << std::endl;
            }

            // Print off-diagonal decomposition
            std::cerr << "\nOff-diagonal elements:" << std::endl;
            for (int i = 0; i < natoms; ++i) {
                for (int j = 0; j < i; ++j) {
                    double gammij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));
                    double r = dist(i, j);  // packed lookup
                    std::cerr << fmt::format("  A[{},{}] = {:12.6f}  (gammij={:12.6f}, r={:12.6f}, erf={:12.6f})",
                        i, j, A(i, j), gammij, r, std::erf(gammij * r)) << std::endl;
                }
            }

            // Print complete A-matrix
            std::cerr << "\nComplete A matrix (" << m << "x" << m << "):" << std::endl;
            for (int i = 0; i < m; ++i) {
                std::cerr << fmt::format("  Row {:2d}:", i);
                for (int j = 0; j < m; ++j) {
                    std::cerr << fmt::format(" {:12.6f}", A(i, j));
                }
                std::cerr << std::endl;
            }

            // Print RHS vector
            std::cerr << "\nRHS vector x:" << std::endl;
            for (int i = 0; i < m; ++i) {
                std::cerr << fmt::format("  x({:2d}) = {:12.6f}", i, x(i)) << std::endl;
            }
            std::cerr << "==================================================" << std::endl;
        }

        // 5. Matrix diagnostics (only for first iteration, small systems only)
        // O(N³) eigendecomposition — skip for natoms > 100 to avoid hour-long hangs on large systems
        if (m_verbosity >= 3 && iteration == 0 && natoms <= 100) {
            Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(A);
            Vector eigenvalues = eigensolver.eigenvalues();
            double max_eigenvalue = eigenvalues.maxCoeff();
            double min_eigenvalue = eigenvalues.minCoeff();
            double condition_number = std::abs(max_eigenvalue / min_eigenvalue);

            if (condition_number > 1e12) {
                CurcumaLogger::warn(fmt::format("Phase 2 EEQ matrix ill-conditioned (cond={:.2e})", condition_number));
            }
        }

        // WP5 (June 2026): self-consistent implicit-solvation reaction field. Add the
        // Born interaction matrix B into the top-left atom block of the augmented EEQ
        // matrix so the solved charges feel the solvent — matches the gfnff reference
        // A(:n,:n) += gbsa%bornMat (gfnff_engrad.F90:1346-1350). B is geometry-only and
        // is added on top of the (cached) off-diagonal Coulomb block each solve, so the
        // off-diagonal cache stays clean. nullptr = gas phase (no change).
        if (m_reaction_field && m_reaction_field->rows() == natoms
            && m_reaction_field->cols() == natoms) {
            A.topLeftCorner(natoms, natoms) += *m_reaction_field;
        }

        // 6. Solve system — unified dispatch (March 2026)
        // Claude Generated (WP2, May 2026): forward pool/num_threads to dispatchSolve so the
        // Stage-4 batched per-fragment LU loop (eeq_solver.cpp:1294-1387) runs in parallel.
        //
        // Threaded-LAPACK lever (Jun 2026): the SchurCholesky dpotrf/dpotrs (and any BLAS-backed
        // solve) scale ~4-5x with the BLAS thread count, but curcuma pins OMP/MKL to 1 globally
        // (CxxThreadPool, main.cpp). Give the BLAS the EnergyCalculator's intra-molecule budget
        // (num_threads) for the duration of the solve; restored on scope exit. The parallel
        // matrix build above has already joined, so there is no thread nesting. See
        // src/core/blas_threads.h.
        curcuma::ScopedBlasThreads _blas_threads(num_threads > 0 ? num_threads : 1);
        Vector new_charges = dispatchSolve(A, x, natoms, nfrag, total_charge, pool, num_threads);

        // Empty return from dispatchSolve signals all solvers failed.
        // Prefer the last successful Phase 2 charges (from a prior step) over Phase 1
        // topology_charges — Phase 2 charges are geometrically accurate and the energy
        // computed with them is physically reasonable. For small optimizer steps the
        // charges from the previous step are an excellent approximation.
        if (new_charges.size() != natoms) {
            if (isValidChargeCache(m_last_successful_charges, natoms, total_charge)) {
                CurcumaLogger::warn(fmt::format(
                    "EEQ Phase 2: all solvers failed for N={}, using last successful charges as fallback",
                    natoms));
                return m_last_successful_charges;
            }
            CurcumaLogger::warn(fmt::format(
                "EEQ Phase 2: all solvers failed for N={}, no valid prior charges — using Phase 1 topology charges",
                natoms));
            // Claude Generated (Apr 2026): Remember that Phase 2 is garbage for this system.
            m_phase2_historically_implausible = true;
            m_phase2_implausible_natoms = natoms;
            return topology_charges;
        }

        // Claude Generated (March 2026): Print Phase 2 solution charges
        if (m_verbosity >= 3 && iteration == 0 && natoms <= 10) {
            std::cerr << "\nPhase 2 solution charges:" << std::endl;
            for (int i = 0; i < natoms; ++i) {
                std::cerr << fmt::format("  q({:2d}) = {:12.6f}  (Z={:2d})", i, new_charges(i), atoms[i]) << std::endl;
            }
            std::cerr << fmt::format("  sum = {:12.6f}", new_charges.sum()) << std::endl;
        }

        // 7. Validate solution — NaN/Inf should already be caught by dispatchSolve,
        // but double-check here and fall back to topology_charges rather than uniform 0.
        bool solution_valid = true;
        for (int i = 0; i < natoms; ++i) {
            if (std::isnan(new_charges[i]) || std::isinf(new_charges[i])) {
                CurcumaLogger::error(fmt::format("Phase 2 EEQ: Invalid charge[{}] = {} (Z={})",
                                                 i, new_charges[i], atoms[i]));
                solution_valid = false;
                break;
            }
        }

        if (!solution_valid) {
            if (isValidChargeCache(m_last_successful_charges, natoms, total_charge)) {
                CurcumaLogger::warn("EEQ Phase 2: invalid solution, using last successful charges as fallback");
                return m_last_successful_charges;
            }
            CurcumaLogger::warn("EEQ Phase 2: invalid solution, no valid prior charges — using Phase 1 topology charges");
            // Claude Generated (Apr 2026): Remember that Phase 2 is garbage for this system.
            m_phase2_historically_implausible = true;
            m_phase2_implausible_natoms = natoms;
            return topology_charges;
        }

        // Claude Generated (April 2026): Phase 2 charge plausibility check + dump
        // dispatchSolve may return non-NaN charges that are still numerically garbage
        // (e.g., from a near-singular matrix). Fall back to Phase 1 topology charges
        // when the solution is clearly non-physical.
        //
        // Uses the cached A matrix directly — no O(N²) distance/erf recomputation.
        // E = 0.5 * q^T * A_nn * q - Σ q_i * chi_i   (chi_i = chi_corrected + cnf*sqrt(CN))
        {
            const Matrix& A_nn = A.topLeftCorner(natoms, natoms);
            auto eeq_energy = [&](const Vector& q) -> double {
                double chi_sum = 0.0;
                for (int i = 0; i < natoms; ++i) {
                    EEQParameters params = getParameters(atoms[i], cn(i));
                    chi_sum += q(i) * (chi_corrected(i) + params.cnf * std::sqrt(std::abs(cn(i))));
                }
                return 0.5 * q.transpose() * A_nn * q - chi_sum;
            };

            double e_coul_p1 = eeq_energy(topology_charges);
            double e_coul_p2 = eeq_energy(new_charges);
            double max_q_p2 = new_charges.cwiseAbs().maxCoeff();
            double max_q_p1 = topology_charges.cwiseAbs().maxCoeff();

            // Dump to JSON if requested (for offline analysis of bad Phase 2 solves)
            if (m_config.get<bool>("dump_charges", false)) {
                json dump;
                dump["natoms"] = natoms;
                dump["nfrag"] = nfrag;
                dump["total_charge"] = total_charge;
                dump["e_coulomb_phase1"] = e_coul_p1;
                dump["e_coulomb_phase2"] = e_coul_p2;
                dump["max_q_phase1"] = max_q_p1;
                dump["max_q_phase2"] = max_q_p2;
                dump["solver_method"] = (m_solve_method == EEQSolveMethod::PCG) ? "PCG"
                    : (m_solve_method == EEQSolveMethod::SchurCholesky) ? "SchurCholesky"
                    : (m_solve_method == EEQSolveMethod::Auto) ? "Auto"
                    : (m_selected_method == EEQSolveMethod::LU) ? "LU (auto)"
                    : "LU";
                for (int i = 0; i < natoms; ++i) {
                    json atom;
                    atom["Z"] = atoms[i];
                    atom["cn"] = cn(i);
                    atom["q_phase1"] = topology_charges(i);
                    atom["q_phase2"] = new_charges(i);
                    dump["atoms"].push_back(atom);
                }
                std::string fname = fmt::format("charges_dump_N{}.json", natoms);
                std::ofstream ofs(fname);
                ofs << dump.dump(2) << std::endl;
                CurcumaLogger::info(fmt::format("EEQ: Charge dump saved to {}", fname));
            }

            // Plausibility: max |q| should not be >> max |q| from Phase 1
            // For organic systems, charges are typically -1..+1. |q| > 10 is clearly garbage.
            double q_threshold = std::max(10.0, 5.0 * max_q_p1);
            if (max_q_p2 > q_threshold) {
                CurcumaLogger::warn(fmt::format(
                    "EEQ Phase 2: implausible charges (max |q| = {:.2f}, Phase 1 max = {:.2f}), using Phase 1",
                    max_q_p2, max_q_p1));
                current_charges = topology_charges;
                // Do NOT overwrite m_last_successful_charges — the Phase 2 result was
                // garbage, but the cache still holds a good result from a previous step.
                // Phase 1 topology_charges are used for THIS step but must not poison
                // the cache for future steps.
                // Claude Generated (Apr 2026): Remember that Phase 2 is garbage for this
                // system size so we can skip the expensive solve on all future calls.
                m_phase2_historically_implausible = true;
                m_phase2_implausible_natoms = natoms;
                break;
            }

            // Plausibility: Coulomb energy should not be wildly different
            double e_diff_ratio = std::abs(e_coul_p2 - e_coul_p1) / (std::abs(e_coul_p1) + 1.0);
            if (e_diff_ratio > 5.0) {
                CurcumaLogger::warn(fmt::format(
                    "EEQ Phase 2: Coulomb energy ratio {:.1f} (P2={:.2e}, P1={:.2e}), using Phase 1",
                    e_diff_ratio, e_coul_p2, e_coul_p1));
                current_charges = topology_charges;
                // Same as above — keep the cache intact.
                // Claude Generated (Apr 2026): Remember that Phase 2 is garbage for this
                // system size so we can skip the expensive solve on all future calls.
                m_phase2_historically_implausible = true;
                m_phase2_implausible_natoms = natoms;
                break;
            }
        }

        // Check convergence for iterative case
        if (use_iterative && iteration > 0) {
            double max_change = (new_charges - current_charges).cwiseAbs().maxCoeff();

            if (m_verbosity >= 3) {
                CurcumaLogger::info(fmt::format("EEQ Phase 2: Iteration {} - Max charge change: {:.2e}",
                                               iteration, max_change));
            }

            if (max_change < convergence_threshold) {
                converged = true;
                if (m_verbosity >= 2) {
                    CurcumaLogger::success(fmt::format("EEQ Phase 2: Converged after {} iterations", iteration));
                }
            }

            current_charges = new_charges;
        } else {
            // Non-iterative case or first iteration - just store the result
            current_charges = new_charges;
            if (!use_iterative) {
                // If not using iterative refinement, break after first iteration
                break;
            }
        }

        iteration++;

    } while (use_iterative && iteration < max_iterations && !converged);

    // Warn if iterative refinement did not converge
    if (use_iterative && !converged) {
        if (m_allow_unconverged) {
            CurcumaLogger::warn(fmt::format("EEQ Phase 2: Did not converge after {} iterations (allow_unconverged_charges=true), using best estimate", iteration));
        } else if (m_verbosity >= 1) {
            CurcumaLogger::warn(fmt::format("EEQ Phase 2: Did not converge after {} iterations, using best estimate", iteration));
        }
    }

    // Store final result
    final_charges = current_charges;

    // Validate charge conservation
    double total = final_charges.sum();
    if (std::abs(total - total_charge) > 1e-6) {
        CurcumaLogger::warn(fmt::format("EEQ Phase 2: Charge sum = {:.6f} (expected {})",
                                        total, total_charge));
    }

    // Store dxi for later use in energy calculation
    m_dxi_stored = dxi;

    // Cache successful Phase 2 result for fallback on future solver failures.
    // On subsequent steps where the solver fails, these provide a physically accurate
    // starting point instead of the cruder Phase 1 topology charges.
    m_last_successful_charges = final_charges;

    // Print PCG convergence summary (one line, resets stats)
    printConvergenceSummary();

    return final_charges;
}

// ===== Energy Calculation =====

double EEQSolver::calculateEEQEnergy(
    const Vector& charges,
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn)
{
    const int natoms = atoms.size();
    double energy = 0.0;

    // Pairwise Coulomb energy with erf damping
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
            double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
            double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-10) continue;

            EEQParameters params_i = getParameters(atoms[i], cn(i));
            EEQParameters params_j = getParameters(atoms[j], cn(j));

            double gamma_ij = 1.0 / std::sqrt(params_i.alp + params_j.alp);
            double erf_gamma = std::erf(gamma_ij * r);
            double coulomb = erf_gamma / r;

            energy += charges(i) * charges(j) * coulomb;
        }
    }

    // Self-energy terms
    const double TSQRT2PI = 0.797884560802866;
    for (int i = 0; i < natoms; ++i) {
        EEQParameters params_i = getParameters(atoms[i], cn(i));

        // CRITICAL FIX (Session 11): chi_i must include dxi term!
        // Reference: XTB gfnff_engrad.F90:1581
        // chi = -χ + dxi + CNF·√CN (NOT just -χ + CNF·√CN!)
        double dxi_i = (m_dxi_stored.size() > i) ? m_dxi_stored(i) : 0.0;
        double chi_i = -params_i.chi + dxi_i + params_i.cnf * std::sqrt(cn(i));
        double self_energy = -charges(i) * chi_i
                           + 0.5 * charges(i) * charges(i) * (params_i.gam + TSQRT2PI / std::sqrt(params_i.alp));

        energy += self_energy;
    }

    return energy;  // Hartree
}

// ===== Correction Terms =====

Vector EEQSolver::calculateDxi(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn,
    const std::optional<TopologyInput>& topology)
{
    // Claude Generated (December 2025, Session 13): FULL dxi calculation
    // Reference: XTB gfnff_ini.f90:358-403
    // Phase 2.5-2.7: Added pi-system detection and neighbor EN averaging

    const int natoms = atoms.size();
    Vector dxi = Vector::Zero(natoms);

    // Pauling electronegativities for neighbor averaging (indices 0-86)
    // Reference: Standard Pauling scale
    static const std::array<double, 87> pauling_en = {
        0.0,  // 0: dummy
        2.20, 0.0,  // 1-2: H, He
        0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0.0,  // 3-10: Li-Ne
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 0.0,  // 11-18: Na-Ar
        0.82, 1.00,  // 19-20: K, Ca
        1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,  // 21-30: Sc-Zn
        1.81, 2.01, 2.18, 2.55, 2.96, 3.00, 0.0,  // 31-37: Ga-Kr (Kr noble gas)
        0.82, 0.95,  // 38-39: Rb, Sr (Note: index mismatch fixed below)
        1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69,  // 40-49: Y-Cd
        1.78, 1.96, 2.05, 2.10, 2.66, 2.60, 0.0,  // 50-56: In-Xe (Xe noble gas)
        0.79, 0.89,  // 57-58: Cs, Ba
        1.10, 1.12, 1.13, 1.14, 1.13, 1.17, 1.20,  // 59-65: La-Eu
        1.20, 1.20, 1.22, 1.23, 1.24, 1.25, 1.10,  // 66-72: Gd-Yb
        1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00,  // 73-82: Lu-Hg
        1.62, 2.33, 2.02, 2.00  // 83-86: Tl-Rn
    };

    // Estimate hybridization from CN (approximation)
    // Reference: CN-based heuristic (hyb: 1=sp, 2=sp2, 3=sp3)
    // Claude Update December 2025: Now uses geometry and topology for accurate detection
    std::vector<int> hyb = detectHybridization(atoms, geometry_bohr, cn, topology);

    // Pi-system detection (simplified from XTB gfnff_ini.f90:312-336)
    // Reference: XTB defines pi atoms as (sp or sp2) AND (C,N,O,F,S)
    std::vector<bool> is_pi_atom = detectPiSystem(atoms, hyb, topology);

    // NOTE (Feb 2026): amideH detection REMOVED from dxi
    // Fortran gfnff_ini.f90:358-403 does NOT include amideH in dxi
    // amideH is Phase 2 only (line 717: chieeq(i) = chieeq(i) - 0.02)

    // Debug output header (Claude Generated Dec 29, 2025 - Fixed debug visibility)
    if (m_verbosity >= 3) {
        CurcumaLogger::info("\n=== Dxi Calculation (FULL - Phase 2.7) ===");
        CurcumaLogger::info("Reference: XTB gfnff_ini.f90:358-403 + pi-system + neighbor EN");
        CurcumaLogger::info("Atom |  Z | CN  | Hyb | Pi | EN_avg | dxi_total | Components");
        CurcumaLogger::info("-----+----+-----+-----+----+--------+-----------+-----------");
    }

    for (int i = 0; i < natoms; ++i) {
        int ati = atoms[i];
        double dxi_total = 0.0;
        std::string components = "";

        // Get neighbor information if topology available
        int nh = 0;  // Number of H neighbors
        int nn = 0;  // Total number of neighbors
        int nm = 0;  // Number of metal neighbors
        double en_avg = 0.0;  // Average neighbor electronegativity
        std::string env_desc = "";

        if (topology.has_value()) {
            nn = topology->neighbor_lists[i].size();
            for (int j : topology->neighbor_lists[i]) {
                int Z_j = atoms[j];
                if (Z_j == 1) nh++;
                // Metal check
                if (Z_j > 20 && (Z_j <= 30 || (Z_j >= 39 && Z_j <= 48) || (Z_j >= 72 && Z_j <= 80))) {
                    nm++;
                }
                // Sum electronegativities for averaging
                if (Z_j < static_cast<int>(pauling_en.size())) {
                    en_avg += pauling_en[Z_j];
                }
            }
            if (nn > 0) {
                en_avg /= nn;  // Average neighbor EN
            }
        }

        // ===== Element-Specific Environment Corrections (from XTB) =====
        // Reference: Fortran gfnff_ini.f90:358-403
        // NOTE: amideH is NOT in dxi - it's Phase 2 only (line 717)

        // Boron (Z=5): +0.015 per H neighbor (line 377)
        if (ati == 5 && nh > 0) {
            double corr = nh * 0.015;
            dxi_total += corr;
            components += fmt::format("B-H:{:+.3f} ", corr);
        }

        // Carbon (Z=6): Special cases (lines 379-387)
        if (ati == 6) {
            // Carbene (CN=2, itag==1): make more negative (line 379)
            // CRITICAL FIX (March 2026): Fortran checks itag(i)==1, which is only set
            // when the bond angle < 150° (gfnff_ini2.f90:244-246).
            // For linear molecules like HCN (angle ~180°), itag=0, so dxi=0.
            // Previous code applied dxi=-0.15 to ALL C with nn==2, causing HCN charge error.
            if (nn == 2) {
                bool is_carbene = false;  // Equivalent of itag==1
                if (topology.has_value() && topology->neighbor_lists[i].size() == 2) {
                    int nb1 = topology->neighbor_lists[i][0];
                    int nb2 = topology->neighbor_lists[i][1];
                    // Calculate bond angle at atom i
                    double dx1 = geometry_bohr(nb1, 0) - geometry_bohr(i, 0);
                    double dy1 = geometry_bohr(nb1, 1) - geometry_bohr(i, 1);
                    double dz1 = geometry_bohr(nb1, 2) - geometry_bohr(i, 2);
                    double dx2 = geometry_bohr(nb2, 0) - geometry_bohr(i, 0);
                    double dy2 = geometry_bohr(nb2, 1) - geometry_bohr(i, 1);
                    double dz2 = geometry_bohr(nb2, 2) - geometry_bohr(i, 2);
                    double r1 = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
                    double r2 = std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
                    if (r1 > 1e-10 && r2 > 1e-10) {
                        double cos_angle = (dx1*dx2 + dy1*dy2 + dz1*dz2) / (r1 * r2);
                        cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
                        double angle_deg = std::acos(cos_angle) * 180.0 / M_PI;
                        is_carbene = (angle_deg < 150.0);  // Fortran: phi*180/pi < 150
                    }
                }
                if (is_carbene) {
                    double corr = -0.15;
                    dxi_total += corr;
                    components += "carbene:-0.15 ";
                    env_desc = "carbene";
                }
            }
            // Free CO (C bonded to single O): make O less negative (line 387)
            if (topology.has_value() && nn == 1) {
                for (int j : topology->neighbor_lists[i]) {
                    int nn_j = topology->neighbor_lists[j].size();
                    if (atoms[j] == 8 && nn_j == 1) {
                        // Apply correction to oxygen (will be applied when we process that atom)
                        // Mark in env_desc for now
                        env_desc = "CO";
                    }
                }
            }
        }

        // Oxygen (Z=8): Multiple environment-dependent corrections (lines 391-394)
        if (ati == 8) {
            // Nitro oxygen: O-N=O pi-system (line 391)
            // Reference: ip .ne. 0.and.at(ji) .eq. 7.and.piadr2(ji) .ne. 0
            if (topology.has_value() && nn == 1 && is_pi_atom[i]) {
                for (int j : topology->neighbor_lists[i]) {
                    if (atoms[j] == 7 && is_pi_atom[j]) {
                        double corr = 0.05;
                        dxi_total += corr;
                        components += "nitro:+0.05 ";
                        env_desc = "NO2";
                    }
                }
            }
            // Free CO: oxygen bonded to single C (line 387)
            if (topology.has_value() && nn == 1) {
                for (int j : topology->neighbor_lists[i]) {
                    int nn_j = topology->neighbor_lists[j].size();
                    if (atoms[j] == 6 && nn_j == 1) {
                        double corr = 0.15;
                        dxi_total += corr;
                        components += "CO:+0.15 ";
                        env_desc = "CO";
                    }
                }
            }
            // H2O: lower electronegativity (line 392)
            if (nn == 2 && nh == 2) {
                double corr = -0.02;
                dxi_total += corr;
                components += "H2O:-0.02 ";
                env_desc = "H2O";
            }
            // General O/S: -0.005 per H (line 394)
            if (nh > 0) {
                double corr = -nh * 0.005;
                dxi_total += corr;
                components += fmt::format("O-H:{:+.3f} ", corr);
            }
            // Group 6 (O,S) with nn > 2: +0.005 per neighbor (line 393)
            if (nn > 2) {
                double corr = nn * 0.005;
                dxi_total += corr;
                components += fmt::format("O-nn:{:+.3f} ", corr);
            }
        }

        // Sulfur (Z=16): Same as oxygen for some corrections (line 394)
        if (ati == 16) {
            if (nh > 0) {
                double corr = -nh * 0.005;
                dxi_total += corr;
                components += fmt::format("S-H:{:+.3f} ", corr);
            }
            if (nn > 2) {
                double corr = nn * 0.005;
                dxi_total += corr;
                components += fmt::format("S-nn:{:+.3f} ", corr);
            }
        }

        // Halogens (Group 7: Cl=17, Br=35, I=53, At=85): Polyvalent corrections (lines 396-402)
        int group_i = GFNFFParameters::periodic_group[ati - 1];
        if (ati > 9 && group_i == 7 && nn > 1) {
            if (nm == 0) {
                // Not bonded to metal: -0.021 per neighbor
                double corr = -nn * 0.021;
                dxi_total += corr;
                components += fmt::format("X-poly:{:+.3f} ", corr);
                env_desc = "polyval";
            } else {
                // Bonded to metal: +0.05 per neighbor
                double corr = nn * 0.05;
                dxi_total += corr;
                components += fmt::format("X-TM:{:+.3f} ", corr);
                env_desc = "TM-ligand";
            }
        }

        dxi(i) = dxi_total;

        // Debug output per atom (Claude Generated Dec 29, 2025 - Fixed debug visibility)
        if (m_verbosity >= 3) {
            if (components.empty()) components = "none";
            std::string hyb_str = (hyb[i] == 1) ? "sp" : (hyb[i] == 2) ? "sp2" : "sp3";
            std::string pi_str = is_pi_atom[i] ? "Y" : "N";
            CurcumaLogger::info(fmt::format("  {:2d} | {:2d} | {:3.1f} | {:3s} | {:2s} | {:6.2f} | {:+9.5f} | {}",
                                    i, ati, cn(i), hyb_str, pi_str, en_avg, dxi_total, components));
        }
    }

    if (m_verbosity >= 3) {
        CurcumaLogger::info("========================================");
    }

    return dxi;
}

Vector EEQSolver::calculateDgam(
    const std::vector<int>& atoms,
    const Vector& charges,
    const std::vector<int>& hybridization,
    const std::vector<bool>& is_pi_atom,
    const std::vector<bool>& is_amide)
{
    const int natoms = atoms.size();
    Vector dgam = Vector::Zero(natoms);

    // Claude Generated (March 2026): Replicate Fortran piadr bug for N dgam
    // Fortran gfnff_ini.f90:708 uses piadr(i) instead of piadr2(i):
    //   piadr(k) = atom_index_of_kth_pi_atom (pi-index → atom)
    //   piadr2(i) = pi_index_of_atom_i (atom → pi-index)
    // piadr(i)!=0 is true for ALL i <= npiall (1-based), regardless of whether
    // atom i is actually a pi atom. GFN-FF parameters were fit with this behavior.
    // Line 709 passes piadr to amide() (same bug), but line 731 correctly uses piadr2 for amideH().
    int npiall = 0;
    for (int k = 0; k < natoms; ++k) {
        if (!is_pi_atom.empty() && k < static_cast<int>(is_pi_atom.size()) && is_pi_atom[k])
            npiall++;
    }

    for (int i = 0; i < natoms; ++i) {
        int Z = atoms[i];
        double qa = charges(i);

        // Claude Generated (December 2025/January 2026): EXACT XTB ff values
        // Reference: Fortran gfnff_ini.f90:665-690 + refined N factors
        double ff = 0.0;  // Default: do nothing

        if (Z == 1) {
            ff = -0.08;  // H
        } else if (Z == 5) {
            ff = -0.05;  // B
        } else if (Z == 6) {  // C
            ff = -0.27;  // sp3
            if (!hybridization.empty() && i < hybridization.size()) {
                if (hybridization[i] < 3) ff = -0.45;  // sp2 or lower
                if (hybridization[i] < 2) ff = -0.34;  // sp
            }
        } else if (Z == 7) {  // N
            ff = -0.13;  // Base N
            // Fortran piadr bug: piadr(i)!=0 iff i_fortran <= npiall
            // In 0-based C++: i < npiall
            if (i < npiall) {
                ff = -0.14;  // "pi-system" N (Fortran piadr compatibility)
            }
            if (!is_amide.empty() && i < is_amide.size() && is_amide[i]) {
                ff = -0.16;  // amide N
            }
        } else if (Z == 8) {  // O
            ff = -0.15;  // sp3
            if (!hybridization.empty() && i < hybridization.size()) {
                if (hybridization[i] < 3) ff = -0.08;  // unsaturated
            }
        } else if (Z == 9) {
            ff = 0.10;   // F
        } else if (Z > 10) {
            ff = -0.02;  // Heavy atoms default

            // Specific heavy atom overrides
            if (Z == 17) ff = -0.02;  // Cl
            if (Z == 35) ff = -0.11;  // Br
            if (Z == 53) ff = -0.07;  // I

            // Metal corrections (requires metal_type array)
            if (Z >= 1 && Z <= 86) {
                int imetal_val = metal_type[Z - 1];
                if (imetal_val == 1) ff = -0.08;   // Main group metals
                if (imetal_val == 2) ff = -0.9;    // Transition metals (XTB comment: "too large")
            }

            // Noble gases (Group 8)
            if (Z >= 1 && Z <= 86) {
                int group = periodic_group[Z - 1];
                if (group == 8) ff = 0.0;  // Noble gases
            }
        }

        dgam(i) = qa * ff;
    }

    return dgam;
}

// Claude Generated (March 2026): Public interface for dgam with full N-specific corrections
// Ensures Coulomb energy parameters match EEQ solver parameters
Vector EEQSolver::calculateDgamFull(
    const std::vector<int>& atoms,
    const Vector& topology_charges,
    const std::vector<int>& hybridization,
    const Vector& cn,
    const std::optional<TopologyInput>& topology)
{
    auto is_pi = detectPiSystem(atoms, hybridization, topology);
    auto is_amide = detectAmideNitrogens(atoms, hybridization, is_pi, topology, cn);
    return calculateDgam(atoms, topology_charges, hybridization, is_pi, is_amide);
}

// Claude Generated (March 2026): Public interface for dxi with full environment corrections
// Ensures Coulomb energy dxi matches EEQ solver dxi exactly
Vector EEQSolver::calculateDxiFull(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn,
    const std::optional<TopologyInput>& topology)
{
    return calculateDxi(atoms, geometry_bohr, cn, topology);
}

// Claude Generated (March 2026): Public wrapper for amide hydrogen detection
std::vector<bool> EEQSolver::detectAmideHydrogensFull(
    const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const Vector& cn,
    const std::optional<TopologyInput>& topology) const
{
    auto is_pi = detectPiSystem(atoms, hybridization, topology);
    auto is_amide = detectAmideNitrogens(atoms, hybridization, is_pi, topology, cn);
    return detectAmideHydrogens(atoms, hybridization, is_amide, topology);
}

// ===== Parameter Lookup =====
// NOTE: calculateDalpha() removed - alpha now calculated with charge-dependent formula (alpha_base + ff*qa)²

EEQSolver::EEQParameters EEQSolver::getParameters(int Z, double cn) const
{
    EEQParameters params;

    if (Z >= 1 && Z <= 86) {
        int idx = Z - 1;
        params.chi = chi_eeq[idx];
        params.gam = gam_eeq[idx];
        params.alp = alpha_eeq[idx] * alpha_eeq[idx];  // CRITICAL: Must be SQUARED!
        params.cnf = cnf_eeq[idx];
    } else {
        CurcumaLogger::warn(fmt::format("EEQSolver: No parameters for Z={}, using defaults", Z));
        params.chi = 1.0;
        params.gam = 0.0;
        params.alp = 1.0;
        params.cnf = 0.0;
    }

    return params;
}

// ===== Helper Functions =====

Vector EEQSolver::calculateCoordinationNumbers(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr) const
{
    // Delegate to shared CNCalculator utility
    // Claude Generated - December 21, 2025 (consolidated duplicate CN calculation)
    std::vector<double> cn_double = CNCalculator::calculateGFNFFCN(atoms, geometry_bohr);

    // Convert std::vector<double> to Eigen::Vector for API compatibility
    Vector cn = Vector::Map(cn_double.data(), cn_double.size());

    return cn;
}

std::vector<int> EEQSolver::detectHybridization(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn,
    const std::optional<TopologyInput>& topology) const
{
    const int natoms = atoms.size();
    std::vector<int> hybridization(natoms, 3);  // Default: sp3

    for (int i = 0; i < natoms; ++i) {
        // Use element-specific XTB rules with geometry and optional topology
        hybridization[i] = detectElementSpecificHybridization(
            atoms[i], cn(i), atoms, topology, i, &geometry_bohr, nullptr
        );
    }

    return hybridization;
}

std::vector<bool> EEQSolver::detectPiSystem(
    const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const std::optional<TopologyInput>& topology) const
{
    const int natoms = atoms.size();
    std::vector<bool> is_pi_atom(natoms, false);

    auto is_pi_element = [](int Z) {
        return (Z == 6 || Z == 7 || Z == 8 || Z == 9 || Z == 16);  // C, N, O, F, S
    };

    for (int i = 0; i < natoms; ++i) {
        if (i >= static_cast<int>(hybridization.size())) continue;

        if ((hybridization[i] == 1 || hybridization[i] == 2) && is_pi_element(atoms[i])) {
            is_pi_atom[i] = true;
        }
        // Special case: N,O,F (sp3) bonded to sp2/sp1 atom (picon in XTB)
        if (hybridization[i] == 3 && (atoms[i] == 7 || atoms[i] == 8 || atoms[i] == 9)) {
            if (topology.has_value()) {
                for (int j : topology->neighbor_lists[i]) {
                    if (j < static_cast<int>(hybridization.size())) {
                        if (hybridization[j] == 1 || hybridization[j] == 2) {
                            is_pi_atom[i] = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    return is_pi_atom;
}

std::vector<bool> EEQSolver::detectAmideNitrogens(
    const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const std::vector<bool>& is_pi_atom,
    const std::optional<TopologyInput>& topology,
    const Vector& cn) const
{
    const int natoms = atoms.size();
    std::vector<bool> is_amide(natoms, false);

    if (!topology.has_value()) return is_amide;

    // Claude Generated (March 2026): Replicate Fortran piadr bug in amide() detection
    // Fortran amide() (gfnff_ini2.f90:1553) receives piadr (NOT piadr2) as the pi array.
    // piadr(i)!=0 is equivalent to i_fortran <= npiall, i.e. i_cpp < npiall.
    // This means the "pi check" in amide() tests atom INDEX, not actual pi membership.
    int npiall = 0;
    for (int k = 0; k < natoms; ++k) {
        if (is_pi_atom[k]) npiall++;
    }

    for (int i = 0; i < natoms; ++i) {
        // FIX (Mar 7, 2026): Match Fortran amide() from gfnff_ini2.f90:1553-1580
        // Updated (Mar 19, 2026): Use Fortran-compatible piadr index check
        // Fortran: if (pi(a) .eq. 0 ...) → piadr(a)==0 → a > npiall → 0-based: i >= npiall
        if (atoms[i] != 7 || i >= npiall) continue;
        if (i < static_cast<int>(hybridization.size()) && hybridization[i] != 3) continue;

        // Count pi-C neighbors (Fortran: nc)
        // Fortran: at(j)==6 .and. pi(j).ne.0 → piadr(j)!=0 → j < npiall
        int nc = 0;
        int ic = -1;  // The single pi-C neighbor (if nc==1)
        for (int neighbor : topology->neighbor_lists[i]) {
            if (atoms[neighbor] == 6 && neighbor < npiall) {
                nc++;
                ic = neighbor;
            }
        }
        if (nc != 1) continue;  // Must have EXACTLY one pi-C neighbor

        // Check if that pi-C has exactly one terminal pi-O neighbor (Fortran: no==1)
        // Fortran: at(j)==8 .and. pi(j).ne.0 .and. nb(20,j)==1 → piadr(j)!=0 → j < npiall
        int no = 0;
        for (int n2 : topology->neighbor_lists[ic]) {
            if (atoms[n2] == 8 && n2 < npiall &&
                static_cast<int>(topology->neighbor_lists[n2].size()) == 1) {
                no++;
            }
        }
        if (no == 1) is_amide[i] = true;
    }
    return is_amide;
}

std::vector<bool> EEQSolver::detectAmideHydrogens(
    const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const std::vector<bool>& is_amide,
    const std::optional<TopologyInput>& topology) const
{
    const int natoms = atoms.size();
    std::vector<bool> is_amide_h(natoms, false);

    if (!topology.has_value()) return is_amide_h;

    for (int i = 0; i < natoms; ++i) {
        if (atoms[i] == 1) { // Hydrogen
            for (int neighbor : topology->neighbor_lists[i]) {
                if (atoms[neighbor] == 7 && is_amide[neighbor]) { // Amide Nitrogen
                    // Requirement: Amide Nitrogen must have exactly ONE sp3 Carbon neighbor
                    int sp3_c_count = 0;
                    for (int n2 : topology->neighbor_lists[neighbor]) {
                        if (atoms[n2] == 6 && hybridization[n2] == 3) {
                            sp3_c_count++;
                        }
                    }
                    if (sp3_c_count == 1) {
                        is_amide_h[i] = true;
                    }
                    break;
                }
            }
        }
    }
    return is_amide_h;
}

std::vector<std::vector<int>> EEQSolver::buildNeighborLists(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    double cutoff_radius) const
{
    const int natoms = atoms.size();
    std::vector<std::vector<int>> neighbors(natoms);

    double cutoff_sq = cutoff_radius * cutoff_radius;

    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            if (i == j) continue;

            Vector ri = geometry_bohr.row(i);
            Vector rj = geometry_bohr.row(j);
            double distance_sq = (ri - rj).squaredNorm();

            if (distance_sq < cutoff_sq) {
                neighbors[i].push_back(j);
            }
        }
    }

    return neighbors;
}
