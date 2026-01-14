/*
 * <Hückel Solver for GFN-FF π-Bond Order Calculation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (January 14, 2026) - Phase 1: Full Hückel Implementation
 *
 * Port from Fortran:
 * - gfnff_ini.f90:928-1062 (main Hückel algorithm)
 * - gfnff_qm.f90:39-154 (diagonalization and density matrix)
 * - gfnff_param.f90:818-840 (parameters)
 */

#include "huckel_solver.h"
#include "src/core/curcuma_logger.h"

#include <Eigen/Eigenvalues>
#include <fmt/format.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <iostream>

// ============================================================================
// Main Entry Point
// ============================================================================

std::vector<double> HuckelSolver::calculatePiBondOrders(
    const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const std::vector<int>& pi_fragments,
    const std::vector<double>& charges,
    const std::vector<std::pair<int,int>>& bonds,
    const Eigen::MatrixXd& distances,
    const std::vector<int>& itag)
{
    const int natoms = static_cast<int>(atoms.size());
    const int nbonds = static_cast<int>(bonds.size());

    // Initialize π-bond orders to zero
    // Triangular storage: pbo[huckel_lin(i,j)] for atom pair (i,j)
    int ntriangular = natoms * (natoms + 1) / 2;
    std::vector<double> pbo(ntriangular, 0.0);

    // Create itag vector if not provided (all zeros)
    std::vector<int> tags = itag;
    if (tags.empty()) {
        tags.resize(natoms, 0);
    }

    // Find unique π-system IDs (excluding 0 which means "not in π-system")
    std::set<int> pi_system_ids;
    for (int frag : pi_fragments) {
        if (frag > 0) {
            pi_system_ids.insert(frag);
        }
    }

    int picount = static_cast<int>(pi_system_ids.size());
    if (picount == 0) {
        if (m_verbosity >= 2) {
            CurcumaLogger::info("HuckelSolver: No π-systems found, returning zero bond orders");
        }
        return pbo;
    }

    if (m_verbosity >= 1) {
        CurcumaLogger::info(fmt::format("HuckelSolver: Processing {} π-system(s)", picount));
    }

    // Process each π-system separately
    for (int pis : pi_system_ids) {
        // ====================================================================
        // Step 1: Identify atoms in this π-system and count electrons
        // Port from gfnff_ini.f90:946-980
        // ====================================================================

        std::vector<int> pi_atoms;       // Original atom indices in this system
        std::vector<int> pi_atom_map(natoms, -1);  // Original → π-system index
        std::vector<int> pi_electrons(natoms, 0);  // π-electrons per atom

        int nelpi = 0;  // Total π-electrons

        for (int k = 0; k < natoms; k++) {
            if (pi_fragments[k] == pis) {
                int pi_idx = static_cast<int>(pi_atoms.size());
                pi_atoms.push_back(k);
                pi_atom_map[k] = pi_idx;

                // Count π-electrons for this atom
                int nel_atom = countPiElectrons(atoms[k], hybridization[k], tags[k]);
                pi_electrons[k] = nel_atom;
                nelpi += nel_atom;

                // Cap at 2 electrons per atom (as in Fortran)
                if (pi_electrons[k] > 2) {
                    pi_electrons[k] = 2;
                }
            }
        }

        int npi = static_cast<int>(pi_atoms.size());

        // Skip if too small (need at least 2 atoms and 1 electron)
        if (npi < 2 || nelpi < 1) {
            if (m_verbosity >= 2) {
                CurcumaLogger::info(fmt::format("HuckelSolver: Skipping π-system {} (npi={}, nel={})",
                                  pis, npi, nelpi));
            }
            continue;
        }

        if (m_verbosity >= 2) {
            CurcumaLogger::info(fmt::format("HuckelSolver: π-system {} has {} atoms, {} electrons",
                              pis, npi, nelpi));
        }

        // ====================================================================
        // Step 2: Iterative Hückel loop
        // Port from gfnff_ini.f90:985-1025
        // ====================================================================

        double E_old = 0.0;
        Eigen::MatrixXd P_old = Eigen::MatrixXd::Constant(npi, npi, 2.0/3.0);  // Benzene reference

        for (int iter = 0; iter < maxhiter; iter++) {
            // Build Hamiltonian with P-dependent off-diagonal coupling
            Eigen::MatrixXd H = buildHamiltonian(
                pi_atoms, pi_atom_map, pi_electrons,
                atoms, hybridization, charges, bonds, distances,
                P_old
            );

            // Solve eigenvalue problem and get density matrix
            // H is modified in place to contain the density matrix
            double E_new = solveAndBuildDensity(H, nelpi);

            if (m_verbosity >= 3) {
                CurcumaLogger::info(fmt::format("  Iter {}: E = {:.6f}", iter + 1, E_new));
            }

            // Check convergence
            if (std::abs(E_new - E_old) < conv_threshold) {
                if (m_verbosity >= 2) {
                    CurcumaLogger::info(fmt::format("  Converged after {} iterations", iter + 1));
                }
                P_old = H;  // H now contains density matrix
                break;
            }

            P_old = H;  // H now contains density matrix
            E_old = E_new;
        }

        // ====================================================================
        // Step 3: Extract π-bond orders from density matrix
        // Port from gfnff_ini.f90:1043-1055
        // ====================================================================

        for (const auto& [ii, jj] : bonds) {
            int ia = pi_atom_map[ii];
            int ja = pi_atom_map[jj];

            if (ia >= 0 && ja >= 0) {
                // Bond order is the density matrix element
                double bond_order = P_old(ja, ia);
                pbo[huckel_lin(ii, jj)] = bond_order;

                if (m_verbosity >= 3) {
                    CurcumaLogger::info(fmt::format("  Bond {}-{}: pbo = {:.4f}", ii, jj, bond_order));
                }
            }
        }
    }

    return pbo;
}

// ============================================================================
// π-Electron Counting
// ============================================================================

int HuckelSolver::countPiElectrons(int atom_type, int hyb, int tag) const
{
    // Port from gfnff_ini.f90:959-978
    // Returns number of π-electrons contributed by this atom

    int nel = 0;

    switch (atom_type) {
        case 5:  // Boron
            if (hyb == 1) nel = 1;  // B in borane (sp hybridized)
            break;

        case 6:  // Carbon
            if (tag != 1) nel = 1;  // Skip carbenes (itag=1)
            break;

        case 7:  // Nitrogen
            if (hyb == 2 && tag == 1) {
                // NO₂ group: itag=1 avoids odd electron number
                nel = 1;
            } else if (hyb <= 2) {
                nel = 1;  // sp or sp²
            } else if (hyb == 3) {
                nel = 2;  // sp³ (lone pair)
            }
            break;

        case 8:  // Oxygen
            if (hyb == 1) {
                nel = 1;  // sp (carbonyl)
            } else if (hyb == 2) {
                nel = 1;  // sp² (ether)
            } else if (hyb == 3) {
                nel = 2;  // sp³ (lone pairs)
            }
            break;

        case 9:  // Fluorine
            if (hyb != 1) {
                nel = 2;  // Standard F
            } else {
                nel = 3;  // sp hybridized (fluor-furan+ case)
            }
            break;

        case 16:  // Sulfur
            if (hyb == 1) {
                nel = 1;
            } else if (hyb == 2) {
                nel = 1;
            } else if (hyb == 3) {
                nel = 2;
            }
            break;

        case 17:  // Chlorine
            if (hyb == 0) {
                nel = 2;
            } else if (hyb == 1) {
                nel = 3;
            }
            break;

        default:
            // Other elements: no π-electrons
            nel = 0;
            break;
    }

    return nel;
}

// ============================================================================
// Hamiltonian Construction
// ============================================================================

Eigen::MatrixXd HuckelSolver::buildHamiltonian(
    const std::vector<int>& pi_atoms,
    const std::vector<int>& pi_atom_map,
    const std::vector<int>& pi_electrons,
    const std::vector<int>& atom_types,
    const std::vector<int>& hybridization,
    const std::vector<double>& charges,
    const std::vector<std::pair<int,int>>& bonds,
    const Eigen::MatrixXd& distances,
    const Eigen::MatrixXd& P_old) const
{
    int npi = static_cast<int>(pi_atoms.size());
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(npi, npi);

    // ========================================
    // Diagonal elements (Coulomb integrals)
    // H_ii = hdiag[Z] + q*hueckelp3 - (nel-1)*pilpf
    // Port from gfnff_ini.f90:990-993
    // ========================================

    for (int i = 0; i < npi; i++) {
        int atom_idx = pi_atoms[i];
        int Z = atom_types[atom_idx];
        double q = charges[atom_idx];
        int nel = pi_electrons[atom_idx];

        // Get element-specific diagonal value (default 0 for unlisted elements)
        double h_diag = (Z < 18) ? hdiag[Z] : 0.0;

        // Build diagonal element
        H(i, i) = h_diag + q * hueckelp3 - static_cast<double>(nel - 1) * pilpf;
    }

    // ========================================
    // Off-diagonal elements (exchange integrals)
    // H_ij = -β * (1 - hiter*(2/3 - P_old_ij))
    // Port from gfnff_ini.f90:995-1008
    // ========================================

    for (const auto& [ii, jj] : bonds) {
        int ia = pi_atom_map[ii];
        int ja = pi_atom_map[jj];

        if (ia >= 0 && ja >= 0) {
            int Z_i = atom_types[ii];
            int Z_j = atom_types[jj];
            int hyb_i = hybridization[ii];
            int hyb_j = hybridization[jj];

            // Get element-specific β values
            double h_off_i = (Z_i < 18) ? hoffdiag[Z_i] : 1.0;
            double h_off_j = (Z_j < 18) ? hoffdiag[Z_j] : 1.0;

            // Geometric mean of β values
            double beta = std::sqrt(h_off_i * h_off_j);

            // Small distance distortion to break symmetry in degenerate systems (COT)
            // Uses distance in Bohr with 1e-9 scaling factor
            double dist = distances(ii, jj);
            beta -= 1.0e-9 * dist;

            // Iteration damping factor
            double damping = hiter;

            // Triple bond correction (sp hybridized atoms are less conjugated)
            if (hyb_i == 1) damping *= htriple;
            if (hyb_j == 1) damping *= htriple;

            // P-dependent scaling with benzene (P=2/3) as reference
            // When P = 2/3, full coupling; deviations reduce coupling
            double P_old_ij = P_old(ja, ia);
            double H_ij = -beta * (1.0 - damping * (2.0/3.0 - P_old_ij));

            H(ja, ia) = H_ij;
            H(ia, ja) = H_ij;
        }
    }

    return H;
}

// ============================================================================
// Eigenvalue Solver and Density Matrix
// ============================================================================

double HuckelSolver::solveAndBuildDensity(Eigen::MatrixXd& H, int nel) const
{
    // Port from gfnff_qm.f90:39-154

    int ndim = static_cast<int>(H.rows());

    // ========================================
    // Step 1: Diagonalize Hamiltonian
    // ========================================

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
    if (solver.info() != Eigen::Success) {
        CurcumaLogger::error("HuckelSolver: Eigenvalue decomposition failed");
        return 0.0;
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    // ========================================
    // Step 2: Scale energies to eV
    // Port from gfnff_qm.f90:98
    // ========================================

    eigenvalues *= 0.1 * hartree_to_ev;

    // ========================================
    // Step 3: Compute occupations via Fermi smearing
    // ========================================

    std::vector<double> occ = fermiSmear(eigenvalues, nel, fermi_temp);

    // ========================================
    // Step 4: Check for perfect biradical (anti-aromatic)
    // Port from gfnff_qm.f90:119-129
    // ========================================

    int ihomo = nel / 2;
    if (ihomo > 0 && ihomo < ndim) {
        if (std::abs(occ[ihomo - 1] - occ[ihomo]) < 1e-4) {
            // Perfect biradical detected - break symmetry
            if (m_verbosity >= 2) {
                CurcumaLogger::info("  Perfect biradical detected, breaking symmetry");
            }
            std::fill(occ.begin(), occ.end(), 0.0);
            for (int i = 0; i < nel / 2; i++) {
                occ[i] = 2.0;
            }
        }
    }

    // ========================================
    // Step 5: Compute electronic energy
    // E = Σ occ_i * ε_i
    // ========================================

    double E_el = 0.0;
    for (int i = 0; i < ndim; i++) {
        E_el += occ[i] * eigenvalues(i);
    }

    // ========================================
    // Step 6: Build density matrix
    // P = C * diag(occ) * C^T
    // ========================================

    H = computeDensityMatrix(eigenvectors, occ);

    return E_el;
}

// ============================================================================
// Fermi Smearing
// ============================================================================

std::vector<double> HuckelSolver::fermiSmear(
    const Eigen::VectorXd& eigenvalues,
    int nel,
    double temp) const
{
    // Port from gfnff_qm.f90:157-219

    int norbs = static_cast<int>(eigenvalues.size());
    std::vector<double> occ(norbs, 0.0);

    // If temperature is negligible, use integer occupation
    if (temp < 1.0) {
        for (int i = 0; i < nel / 2 && i < norbs; i++) {
            occ[i] = 2.0;
        }
        if (nel % 2 != 0 && nel / 2 < norbs) {
            occ[nel / 2] = 1.0;
        }
        return occ;
    }

    // Boltzmann factor in eV
    double bkt = boltz_ev * temp;

    // Initial guess for Fermi energy
    int ihomo = std::min(nel, norbs) - 1;
    int ilumo = std::min(nel, norbs - 1);
    double e_fermi = 0.5 * (eigenvalues(ihomo) + eigenvalues(ilumo));

    double target_occ = static_cast<double>(nel);

    // Iterative search for Fermi level
    for (int cycle = 0; cycle < 200; cycle++) {
        double total_occ = 0.0;
        double total_deriv = 0.0;

        for (int i = 0; i < norbs; i++) {
            double x = (eigenvalues(i) - e_fermi) / bkt;

            double fermi_func = 0.0;
            double fermi_deriv = 0.0;

            if (x < 50.0) {
                double exp_x = std::exp(x);
                fermi_func = 1.0 / (exp_x + 1.0);
                fermi_deriv = exp_x / (bkt * (exp_x + 1.0) * (exp_x + 1.0));
            }

            occ[i] = fermi_func;
            total_occ += fermi_func;
            total_deriv += fermi_deriv;
        }

        double delta_fermi = (target_occ - total_occ) / total_deriv;
        e_fermi += delta_fermi;

        if (std::abs(target_occ - total_occ) < 1e-9) {
            break;
        }
    }

    return occ;
}

// ============================================================================
// Density Matrix Construction
// ============================================================================

Eigen::MatrixXd HuckelSolver::computeDensityMatrix(
    const Eigen::MatrixXd& C,
    const std::vector<double>& occ) const
{
    // Port from gfnff_qm.f90:281-304
    // P_μν = Σ_i n_i * C_μi * C_νi

    int ndim = static_cast<int>(C.rows());

    // Create scaled coefficient matrix: Ptmp_μi = C_μi * n_i
    Eigen::MatrixXd Ptmp(ndim, ndim);
    for (int i = 0; i < ndim; i++) {
        for (int mu = 0; mu < ndim; mu++) {
            Ptmp(mu, i) = C(mu, i) * occ[i];
        }
    }

    // P = C * Ptmp^T = C * (C * diag(occ))^T
    Eigen::MatrixXd P = C * Ptmp.transpose();

    return P;
}
