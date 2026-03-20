/*
 * <Native GFN1-xTB Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the GFN1-xTB method developed by:
 *   Stefan Grimme, Christoph Bannwarth, Philip Shushkov
 *   Mulliken Center for Theoretical Chemistry, University of Bonn
 *
 * Reference implementation: TBLite (https://github.com/tblite/tblite)
 *   See: src/tblite/xtb/gfn1.f90
 *
 * Original publication:
 *   S. Grimme, C. Bannwarth, P. Shushkov
 *   J. Chem. Theory Comput. 2017, 13, 1989-2009
 *
 * This program is free software under GPL-3.0
 */

#include "gfn1.h"
#include "diis_accelerator.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "ParallelEigenSolver.hpp"

#include <fmt/format.h>
#include <cmath>
#include <algorithm>

using namespace CurcumaUnit;

// Reuse covalent radii from GFN2
extern double getCovalentRadius(int Z);

namespace {
    // GFN1 specific parameters (from original paper)
    const double GFN1_K1 = 16.0;      // CN steepness
    const double GFN1_K2 = 4.0 / 3.0; // CN range decay

    // Halogen bond elements
    inline bool isHalogen(int Z) {
        return (Z == 9 || Z == 17 || Z == 35 || Z == 53 || Z == 85);  // F, Cl, Br, I, At
    }
}

// =================================================================================
// Constructor / Destructor
// =================================================================================

GFN1::GFN1()
    : m_params()
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_repulsion(0.0)
    , m_energy_coulomb(0.0)
    , m_energy_dispersion(0.0)
    , m_energy_halogen_bond(0.0)
    , m_scf_max_iterations(250)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    // Load complete GFN1 parameter database (15 elements, 25 pairs)
    // Claude Generated: Parameter database integration (November 2025)
    if (!m_param_db.loadCompleteGFN1()) {
        CurcumaLogger::error("Failed to load GFN1 parameter database - using legacy parameters");
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing native GFN1-xTB method");
        CurcumaLogger::param("method", "GFN1-xTB");
        CurcumaLogger::param("reference", "Grimme et al. JCTC 2017, 13, 1989");
        CurcumaLogger::param("parameter_database", fmt::format("{} elements, {} pairs (with halogen bond correction)",
                           m_param_db.getNumElements(), m_param_db.getNumPairs()));
    }
}

GFN1::GFN1(const ArrayParameters& params)
    : m_params(params)
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_repulsion(0.0)
    , m_energy_coulomb(0.0)
    , m_energy_dispersion(0.0)
    , m_energy_halogen_bond(0.0)
    , m_scf_max_iterations(250)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    // Load complete GFN1 parameter database
    // Claude Generated: Parameter database integration (November 2025)
    if (!m_param_db.loadCompleteGFN1()) {
        CurcumaLogger::error("Failed to load GFN1 parameter database - using legacy parameters");
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing GFN1-xTB with custom parameters");
        CurcumaLogger::param("parameter_database", fmt::format("{} elements, {} pairs",
                           m_param_db.getNumElements(), m_param_db.getNumPairs()));
    }
}

// =================================================================================
// QMDriver Interface Implementation
// =================================================================================

bool GFN1::InitialiseMolecule()
{
    if (m_atoms.size() == 0) {
        CurcumaLogger::error("No atoms in molecule for GFN1 initialization");
        return false;
    }

    // Validate all atoms are supported
    for (size_t i = 0; i < m_atoms.size(); ++i) {
        if (!m_params.isValidAtom(m_atoms[i])) {
            CurcumaLogger::error(fmt::format("GFN1 parameters not available for element Z={} at atom {}",
                                            m_atoms[i], i+1));
            return false;
        }
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing GFN1 calculation");
        CurcumaLogger::param("atoms", static_cast<int>(m_atoms.size()));
        CurcumaLogger::param("charge", m_charge);
    }

    // Build basis set
    m_nbasis = buildBasisSet();

    if (m_nbasis == 0) {
        CurcumaLogger::error("Failed to build basis set for GFN1");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("basis_functions", m_nbasis);
        CurcumaLogger::param("electrons", m_num_electrons);
    }

    // Calculate coordination numbers
    m_coordination_numbers = calculateCoordinationNumbers();

    // Build overlap matrix
    m_overlap = MakeOverlap(m_basis);

    // Build core Hamiltonian
    m_hamiltonian = MakeH(m_overlap, m_basis);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("GFN1 initialization complete");
    }

    return true;
}

double GFN1::Calculation(bool gradient)
{
    try {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::info("Starting GFN1 energy calculation");
        }

        // Run SCF convergence
        bool scf_success = runSCF();

        if (!scf_success) {
            CurcumaLogger::error("GFN1 SCF did not converge");
            return 0.0;
        }

        // Calculate energy components
        m_energy_electronic = calculateElectronicEnergy();
        m_energy_repulsion = calculateRepulsionEnergy();
        m_energy_coulomb = calculateCoulombEnergy();
        m_energy_dispersion = calculateDispersionEnergy();  // D3 stub
        m_energy_halogen_bond = calculateHalogenBondCorrection();

        // Total energy
        m_total_energy = m_energy_electronic + m_energy_repulsion +
                        m_energy_coulomb + m_energy_dispersion +
                        m_energy_halogen_bond;

        // Level 1+: Results
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success("GFN1 calculation complete");
            CurcumaLogger::energy_abs(m_total_energy, "GFN1_total_energy");
        }

        // Level 2+: Energy decomposition
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Energy decomposition:");
            CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", m_energy_electronic));
            CurcumaLogger::param("repulsion", fmt::format("{:.6f} Eh", m_energy_repulsion));
            CurcumaLogger::param("coulomb", fmt::format("{:.6f} Eh", m_energy_coulomb));
            CurcumaLogger::param("dispersion", fmt::format("{:.6f} Eh (D3 stub)", m_energy_dispersion));
            CurcumaLogger::param("halogen_bond", fmt::format("{:.6f} Eh", m_energy_halogen_bond));

            double homo = getHOMOEnergy();
            double lumo = getLUMOEnergy();
            double gap = getHOMOLUMOGap();

            CurcumaLogger::param("HOMO", fmt::format("{:.4f} eV", homo));
            CurcumaLogger::param("LUMO", fmt::format("{:.4f} eV", lumo));
            CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} eV", gap));
        }

        // Gradient calculation (if requested)
        if (gradient) {
            m_gradient = calculateGradient();
        }

        return m_total_energy;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("GFN1 calculation failed: {}", e.what()));
        return 0.0;
    }
}

// =================================================================================
// Basis Set Construction (same as GFN2)
// =================================================================================

int GFN1::buildBasisSet()
{
    m_basis.clear();
    m_num_electrons = 0;

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];

        // GFN1 uses minimal valence basis (similar to GFN2)
        int n_shells = 1;
        if (Z > 2) n_shells = 2;
        if (Z > 18) n_shells = 3;

        for (int shell = 0; shell < n_shells; ++shell) {
            // Get parameters for this shell from parameter database
            // Use m_param_db (TBLite parameters) if available, else fallback to m_params
            double zeta = 0.0;
            int principal_qn = shell + 1;  // Default: n = l + 1

            if (m_param_db.hasElement(Z) && m_param_db.getElement(Z).shells.count(shell)) {
                const auto& shell_params = m_param_db.getElement(Z).shells.at(shell);
                zeta = shell_params.gexp;  // Use gexp (Gaussian exponent / Yeff)
            } else {
                zeta = m_params.getYeff(Z, shell);  // Legacy fallback
                principal_qn = static_cast<int>(m_params.getPrincipalQN(Z, shell));
            }

            if (zeta < 1.0e-6) continue;  // Skip unparameterized shells

            int l = shell;

            for (int m = -l; m <= l; ++m) {
                STO::Orbital orbital;

                // Set orbital type based on angular momentum
                if (l == 0) {
                    orbital.type = STO::S;
                } else if (l == 1) {
                    // p-orbitals: m = -1 (py), 0 (pz), +1 (px)
                    if (m == -1) orbital.type = STO::PY;
                    else if (m == 0) orbital.type = STO::PZ;
                    else orbital.type = STO::PX;
                } else if (l == 2) {
                    // d-orbitals: would need proper m → type mapping
                    orbital.type = STO::S;  // Placeholder
                } else {
                    orbital.type = STO::S;  // Fallback
                }

                orbital.x = m_geometry(i, 0) / au;
                orbital.y = m_geometry(i, 1) / au;
                orbital.z = m_geometry(i, 2) / au;
                orbital.zeta = zeta;
                orbital.VSIP = 0.0;
                orbital.atom = i;
                orbital.shell = shell;

                m_basis.push_back(orbital);
            }

            // Fix 2 (March 2026): Valence-only electron count, mirrors GFN2 buildBasisSet
            // Bug was: min(2,Z) and min(6,Z-2) counted core electrons → C got 6, O got 8
            // Correct: only valence electrons fill the minimal basis
            if (Z <= 2) { if (shell == 0) m_num_electrons += Z; }
            else if (Z <= 10) {
                if (shell == 0) m_num_electrons += std::min(2, Z - 2);
                else if (shell == 1) m_num_electrons += std::max(0, Z - 4);
            } else if (Z <= 18) {
                if (shell == 0) m_num_electrons += std::min(2, Z - 10);
                else if (shell == 1) m_num_electrons += std::max(0, Z - 12);
            } else {
                if (shell == 0) m_num_electrons += 2;
                else if (shell == 1) m_num_electrons += 6;
            }
        }
    }

    return static_cast<int>(m_basis.size());
}

Matrix GFN1::MakeOverlap(std::vector<STO::Orbital>& basisset)
{
    Matrix S = Matrix::Zero(basisset.size(), basisset.size());

    for (size_t i = 0; i < basisset.size(); ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double overlap = STO::calculateOverlap(basisset[i], basisset[j]);
            S(i, j) = S(j, i) = overlap;
        }
    }

    return S;
}

// =================================================================================
// Hamiltonian Construction (simplified compared to GFN2)
// =================================================================================

Matrix GFN1::MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset)
{
    // Claude Generated: GFN1 Hamiltonian matrix construction
    // Fix 1 (March 2026): Two-pass loop — mirrors GFN2 MakeH fix.
    // Bug was: single nested loop computed off-diagonal H(i,j) before H(i,i) was set,
    // so (H(i,i)+H(j,j))/2 used 0 for the unset diagonal → wrong off-diagonal elements.
    // Solution: pass 1 sets all diagonals, pass 2 uses them for off-diagonals.

    Matrix H = Matrix::Zero(basisset.size(), basisset.size());

    // Pass 1: set all diagonal (self-energy) elements first
    for (size_t i = 0; i < basisset.size(); ++i) {
        int atom_i = basisset[i].atom;
        int Z_i = m_atoms[atom_i];
        double CN_i = m_coordination_numbers(atom_i);
        int shell_i = 0;
        if (basisset[i].type == STO::PX || basisset[i].type == STO::PY || basisset[i].type == STO::PZ)
            shell_i = 1;
        H(i, i) = getSelfEnergy(Z_i, shell_i, CN_i);
    }

    // Pass 2: off-diagonal elements using already-set diagonals
    for (size_t i = 0; i < basisset.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            double dx = basisset[i].x - basisset[j].x;
            double dy = basisset[i].y - basisset[j].y;
            double dz = basisset[i].z - basisset[j].z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz) * au;

            double scale = getHamiltonianScale(basisset[i], basisset[j], distance);
            double h_ij = scale * S(i, j) * (H(i, i) + H(j, j)) / 2.0;
            H(i, j) = H(j, i) = h_ij;
        }
    }

    return H;
}

double GFN1::getSelfEnergy(int element, int shell, double CN) const
{
    // Claude Generated: GFN1 self-energy with shell-resolved parameters
    // Formula: E_ii = E_base + k_CN * CN (simpler than GFN2)
    // Reference: S. Grimme et al., JCTC 2017, 13, 1989
    // Updated November 2025: Uses real TBLite-derived parameters for 15 elements

    double E = 0.0;
    bool use_real_params = false;

    // Try to use real shell-resolved parameters from database
    if (m_param_db.hasElement(element)) {
        const auto& elem_params = m_param_db.getElement(element);

        // Check if this shell exists for this element
        auto shell_it = elem_params.shells.find(shell);
        if (shell_it != elem_params.shells.end()) {
            const auto& shell_params = shell_it->second;

            // GFN1 self-energy with real parameters
            // Fix 4 (March 2026): sign was inverted — kcn is negative, so + gives more
            // negative self-energy for higher CN (physically correct, mirrors GFN2)
            E = shell_params.selfenergy + shell_params.kcn * CN;
            use_real_params = true;

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param(fmt::format("SelfEnergy[Z={},shell={}]", element, shell),
                                   fmt::format("E={:.6f} Eh (CN={:.2f}) [TBLite params]", E, CN));
            }
        }
    }

    // Fallback to legacy parameters if element/shell not in database
    if (!use_real_params) {
        double E_base = -m_params.getHardness(element);
        double k_CN = m_params.getShellHardness(element) * 0.01;

        // Fix 4 (March 2026): same sign correction for legacy fallback
        E = E_base + k_CN * CN;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("SelfEnergy[Z={},shell={}]", element, shell),
                               fmt::format("E={:.6f} Eh (CN={:.2f}) [legacy fallback]", E, CN));
        }
    }

    return E;
}

double GFN1::getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const
{
    // Claude Generated: GFN1 Hamiltonian scaling with pair-specific parameters
    // GFN1 uses simpler pair interactions than GFN2
    // Reference: S. Grimme et al., JCTC 2017, 13, 1989
    // Updated November 2025: Uses real TBLite pair parameters for 25 element pairs

    int atom_i = fi.atom;
    int atom_j = fj.atom;
    int Z_i = m_atoms[atom_i];
    int Z_j = m_atoms[atom_j];

    double EN_i = m_params.getElectronegativity(Z_i);
    double EN_j = m_params.getElectronegativity(Z_j);
    double delta_EN = EN_i - EN_j;

    double zeta_i = fi.zeta;
    double zeta_j = fj.zeta;
    double z_ij = std::sqrt(2.0 * std::sqrt(zeta_i * zeta_j) / (zeta_i + zeta_j));

    // Try to use real pair-specific parameters
    double k_pair = 1.0;
    double k_shell = 1.0;
    bool use_real_params = false;

    if (m_param_db.hasPair(Z_i, Z_j)) {
        const auto& pair_params = m_param_db.getPair(Z_i, Z_j);

        // Use element-pair specific coupling
        k_pair = pair_params.kpair;

        // Determine shell types for k_shell coupling
        int shell_i = 0;  // s-orbital
        int shell_j = 0;
        if (fi.type == STO::PX || fi.type == STO::PY || fi.type == STO::PZ) shell_i = 1;
        if (fj.type == STO::PX || fj.type == STO::PY || fj.type == STO::PZ) shell_j = 1;

        // Select shell-specific coupling
        if (shell_i == 0 && shell_j == 0) {
            k_shell = pair_params.kshell_ss;  // s-s coupling
        } else if (shell_i == 1 && shell_j == 1) {
            k_shell = pair_params.kshell_pp;  // p-p coupling
        } else {
            k_shell = pair_params.kshell_sp;  // s-p coupling
        }

        use_real_params = true;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("HScale[{}-{}]", Z_i, Z_j),
                               fmt::format("k_pair={:.3f}, k_shell={:.3f} [TBLite]",
                                         k_pair, k_shell));
        }
    } else {
        // Fallback: use Alpha parameter as approximation
        k_pair = std::sqrt(m_params.getAlpha(Z_i) * m_params.getAlpha(Z_j)) * 0.1;
        k_shell = 1.0;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("HScale[{}-{}]", Z_i, Z_j),
                               "using legacy approximation");
        }
    }

    double en_factor = 1.0 + 0.02 * delta_EN * delta_EN;

    double r_bohr = distance / au;
    double poly_r = std::exp(-0.5 * (r_bohr - 3.0));

    return z_ij * k_pair * k_shell * en_factor * poly_r;
}

// =================================================================================
// Coordination Numbers (same as GFN2)
// =================================================================================

Vector GFN1::calculateCoordinationNumbers()
{
    // Claude Generated: GFN1 CN calculation (identical to GFN2)

    Vector CN = Vector::Zero(m_atomcount);

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        double R_cov_A = getCovalentRadius(Z_A);

        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;

            int Z_B = m_atoms[B];
            double R_cov_B = getCovalentRadius(Z_B);
            double R_cov = R_cov_A + R_cov_B;

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (R_AB < 1.0e-6) continue;

            double count = 1.0 / (1.0 + std::exp(-GFN1_K1 * (R_cov / R_AB - 1.0)));
            CN(A) += std::pow(count, GFN1_K2);
        }
    }

    return CN;
}

// =================================================================================
// SCF Procedure (same as GFN2)
// =================================================================================

bool GFN1::runSCF()
{
    // Claude Generated (March 2026): GFN1 SCF with DIIS acceleration
    // DIIS (Pulay 1980) is essential for convergence on molecules with >5 atoms.
    // Without it, simple damping oscillates and never converges for C6H6, caffeine, etc.

    m_density = Matrix::Zero(m_nbasis, m_nbasis);
    double prev_energy = 0.0;
    double energy_diff = 0.0;

    // Compute S^(-1/2) once before SCF loop
    ParallelEigenSolver solver(500, 128, 1.0e-10, false);
    solver.setThreadCount(m_threads);
    if (!solver.computeS_1_2(m_overlap, m_S_inv_sqrt, m_threads)) {
        CurcumaLogger::error("Failed to compute S^(-1/2) for GFN1 SCF");
        return false;
    }

    DIISAccelerator diis(8);  // Keep last 8 Fock matrices
    const int diis_start = 2;  // Start DIIS after this many iterations

    for (int iter = 0; iter < m_scf_max_iterations; ++iter) {
        m_fock = buildFockMatrix(m_density);

        // DIIS: store Fock matrix and extrapolate when enough vectors available
        Matrix F_scf = m_fock;
        if (iter >= diis_start) {
            diis.push(m_fock, m_density, m_overlap);
            if (diis.size() >= 2)
                F_scf = diis.extrapolate();
        }

        // Claude Generated (March 2026): Fixed transformMOs=true (was false)
        // Without back-transform, MOs are in orthonormal basis → P = 2*C*C^T gives Tr(PS) ≠ N_elec
        bool success = solver.solveWithPrecalculatedS_1_2(m_S_inv_sqrt, F_scf, m_energies, m_mo, m_threads, true);

        if (!success) {
            CurcumaLogger::error("Eigenvalue solution failed in GFN1 SCF");
            return false;
        }

        Matrix density_new = buildDensityMatrix(m_mo, m_energies);

        // Energy convergence tracking
        double current_energy = 0.5 * (density_new.cwiseProduct(m_hamiltonian + m_fock)).sum();
        if (iter > 0) energy_diff = std::abs(current_energy - prev_energy);
        prev_energy = current_energy;

        double density_change = (density_new - m_density).norm();
        bool converged = (density_change < m_scf_threshold) && (energy_diff < m_scf_threshold * 0.1);

        if (CurcumaLogger::get_verbosity() >= 2 && (iter % 10 == 0 || converged || iter == 0)) {
            CurcumaLogger::info(fmt::format("GFN1 SCF iter {:3d}: E = {:15.8f} Eh, ΔE = {:8.2e}, ΔP = {:8.2e}",
                                          iter + 1, current_energy, energy_diff, density_change));
        }

        if (converged) {
            m_density = density_new;
            m_scf_converged = true;
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::success(fmt::format("GFN1 SCF converged in {} iterations (ΔE = {:.2e})", iter + 1, energy_diff));
            return true;
        }

        // Iteration 0: use full new density (zero → damping creates artificial charges)
        // All other iterations: damping for stability (DIIS handles Fock extrapolation)
        if (iter == 0) {
            m_density = density_new;
        } else {
            m_density = m_scf_damping * density_new + (1.0 - m_scf_damping) * m_density;
        }
    }

    CurcumaLogger::warn(fmt::format("GFN1 SCF did not converge in {} iterations (ΔE = {:.2e})",
                                   m_scf_max_iterations, energy_diff));
    m_scf_converged = false;
    return false;
}

Matrix GFN1::buildFockMatrix(const Matrix& density)
{
    // Claude Generated (March 2026): Self-consistent Fock matrix with Coulomb potential
    // F = H₀ - 0.5 * S * (V[μ] + V[ν])
    // This is essential for SCF self-consistency — without it, the density never responds
    // to charge redistribution and the SCF just returns the H₀ eigenstates.
    // Reference: GFN2::buildFockMatrix pattern (gfn2.cpp:504-536)
    //
    // GFN1 uses atom-resolved Hubbard (gamma_ss for all shells), not shell-resolved.

    Matrix F = m_hamiltonian;
    Matrix PS = density * m_overlap;

    // Compute shell populations per atom
    std::vector<std::map<int, double>> shell_pop(m_atomcount);
    for (int mu = 0; mu < m_nbasis; ++mu) {
        int A = m_basis[mu].atom;
        int l = (m_basis[mu].type == STO::PX || m_basis[mu].type == STO::PY || m_basis[mu].type == STO::PZ) ? 1 : 0;
        shell_pop[A][l] += PS(mu, mu);
    }

    // Delta-charges per shell and per atom using refocc
    std::vector<std::map<int, double>> dq(m_atomcount);
    Vector atomic_charges = Vector::Zero(m_atomcount);
    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        for (auto& [s, pop] : shell_pop[A]) {
            double ref = 0.0;
            if (m_param_db.hasElement(Z_A) && m_param_db.getElement(Z_A).shells.count(s))
                ref = m_param_db.getElement(Z_A).shells.at(s).refocc;
            dq[A][s] = ref - pop;
            atomic_charges(A) += dq[A][s];
        }
    }
    m_charges = atomic_charges;

    // Skip Coulomb contribution for zero density (first SCF iteration)
    if (density.norm() < 1e-10) return F;

    // Compute Coulomb potential V[μ] for each basis function
    // GFN1: atom-resolved gamma_ss for all shells (simpler than GFN2's shell-resolved)
    std::vector<double> V(m_nbasis, 0.0);
    for (int mu = 0; mu < m_nbasis; ++mu) {
        int A = m_basis[mu].atom;
        double gamma_AA = m_param_db.hasElement(m_atoms[A])
                          ? m_param_db.getElement(m_atoms[A]).gamma_ss
                          : m_params.getHardness(m_atoms[A]);

        for (int B = 0; B < m_atomcount; ++B) {
            double R_AB = (A == B) ? 0.0 : (m_geometry.row(A) - m_geometry.row(B)).norm() / au;
            double gamma_BB = m_param_db.hasElement(m_atoms[B])
                              ? m_param_db.getElement(m_atoms[B]).gamma_ss
                              : m_params.getHardness(m_atoms[B]);
            // Klopman-Ohno Coulomb kernel (same formula as corrected calculateCoulombEnergy)
            double gamma_AB = 1.0 / std::sqrt(R_AB * R_AB + std::pow(2.0 / (gamma_AA + gamma_BB), 2.0));
            V[mu] += atomic_charges(B) * gamma_AB;
        }
    }

    // Modify Fock matrix: F(μ,ν) -= 0.5 * S(μ,ν) * (V[μ] + V[ν])
    for (int mu = 0; mu < m_nbasis; ++mu)
        for (int nu = 0; nu < m_nbasis; ++nu)
            F(mu, nu) -= 0.5 * m_overlap(mu, nu) * (V[mu] + V[nu]);

    return F;
}

Matrix GFN1::buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies)
{
    // Claude Generated (March 2026): Two occupation schemes controlled by m_electronic_temperature
    //   etemp == 0  → integer closed-shell (2/0 occupation)
    //   etemp > 0   → Fermi-Dirac smearing (fractional occupation)
    // Default: 300K (consistent with TBLite). Set via -etemp CLI parameter.
    int nbas = mo_coefficients.rows();

    if (m_electronic_temperature <= 0.0) {
        // Standard closed-shell: P = 2 * C_occ * C_occ^T
        int n_occ = m_num_electrons / 2;
        auto C_occ = mo_coefficients.leftCols(n_occ);
        return 2.0 * C_occ * C_occ.transpose();
    }

    // Fermi-Dirac smearing: f_i = 2 / (1 + exp((ε_i - μ) / kT))
    const double kT = m_electronic_temperature * 3.166808e-6;  // K → Hartree
    int n_elec = m_num_electrons;

    // Find Fermi level by bisection
    double mu_lo = mo_energies.minCoeff() - 1.0;
    double mu_hi = mo_energies.maxCoeff() + 1.0;
    for (int bisect = 0; bisect < 100; ++bisect) {
        double mu = 0.5 * (mu_lo + mu_hi);
        double n_sum = 0.0;
        for (int i = 0; i < nbas; ++i) {
            double x = (mo_energies(i) - mu) / kT;
            n_sum += 2.0 / (1.0 + std::exp(std::min(x, 500.0)));
        }
        if (n_sum > n_elec) mu_hi = mu;
        else mu_lo = mu;
        if (mu_hi - mu_lo < 1e-14) break;
    }
    double mu = 0.5 * (mu_lo + mu_hi);

    Matrix P = Matrix::Zero(nbas, nbas);
    for (int i = 0; i < nbas; ++i) {
        double x = (mo_energies(i) - mu) / kT;
        double f_i = 2.0 / (1.0 + std::exp(std::min(x, 500.0)));
        if (f_i < 1e-12) continue;
        for (int mu_idx = 0; mu_idx < nbas; ++mu_idx)
            for (int nu_idx = 0; nu_idx < nbas; ++nu_idx)
                P(mu_idx, nu_idx) += f_i * mo_coefficients(mu_idx, i) * mo_coefficients(nu_idx, i);
    }
    return P;
}

// =================================================================================
// Energy Components
// =================================================================================

double GFN1::calculateElectronicEnergy() const
{
    // GFN1/xTB electronic energy: E_elec = Tr(P * H_0)
    // Claude Generated (March 2026): Vectorized with Eigen cwiseProduct
    return (m_density.cwiseProduct(m_hamiltonian)).sum();
}

double GFN1::calculateRepulsionEnergy() const
{
    // Claude Generated: GFN1 repulsion energy with TBLite parameters
    // Reference: S. Grimme et al., JCTC 2017, 13, 1989
    // Updated November 2025: Uses real rep_alpha and rep_zeff parameters

    double E_rep = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];

        // Try to get real repulsion parameters
        double alpha_A, Z_eff_A;
        if (m_param_db.hasElement(Z_A)) {
            const auto& elem_A = m_param_db.getElement(Z_A);
            alpha_A = elem_A.rep_alpha;
            Z_eff_A = elem_A.rep_zeff;
        } else {
            // Fallback to legacy parameters
            alpha_A = m_params.getMultipoleRadius(Z_A);
            Z_eff_A = m_params.getElectronegativity(Z_A) * 0.5;
        }

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];

            // Try to get real repulsion parameters
            double alpha_B, Z_eff_B;
            if (m_param_db.hasElement(Z_B)) {
                const auto& elem_B = m_param_db.getElement(Z_B);
                alpha_B = elem_B.rep_alpha;
                Z_eff_B = elem_B.rep_zeff;
            } else {
                // Fallback to legacy parameters
                alpha_B = m_params.getMultipoleRadius(Z_B);
                Z_eff_B = m_params.getElectronegativity(Z_B) * 0.5;
            }

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz) / au;

            // Fix 5 (March 2026): correct GFN1/GFN2 repulsion formula — mirrors GFN2
            // Bugs: arithmetic mean of alpha (should be geometric), sum of Zeff (should be product),
            // missing kexp polynomial (R^1.5 for heavy, R^1.0 for H-H)
            double alpha_AB = std::sqrt(alpha_A * alpha_B);   // geometric mean
            double zeff_AB  = Z_eff_A * Z_eff_B;              // product
            double kexp     = (m_atoms[A] > 2 || m_atoms[B] > 2) ? 1.5 : 1.0;
            double V_rep    = zeff_AB * std::exp(-alpha_AB * std::pow(R_AB, kexp)) / R_AB;

            E_rep += V_rep;
        }
    }

    return E_rep;
}

double GFN1::calculateCoulombEnergy() const
{
    // Claude Generated: GFN1 Coulomb energy with TBLite parameters
    // GFN1 uses simpler ES2 model (no ES3 third-order correction)
    // Reference: S. Grimme et al., JCTC 2017, 13, 1989
    //
    // Fix 3 (March 2026): replace Z_A-based charge with refocc-based delta-charge.
    // Bug: charges(A) = Z_A - mulliken_pop used nuclear charge Z_A as reference.
    // For C (Z=6) with ~4 valence electrons this gave q=+2, for O (Z=8) q=+2,
    // causing E_ES2 ~ 0.5*4*gamma per atom — many Hartree of false positive energy.
    // Correct: delta_q[A][s] = refocc[s] - shell_pop[s]  (mirrors GFN2 line 592)

    Matrix PS = m_density * m_overlap;

    // Accumulate shell populations per atom
    std::vector<std::map<int, double>> shell_pop(m_atomcount);
    for (int mu = 0; mu < m_nbasis; ++mu) {
        int A   = m_basis[mu].atom;
        int l   = (m_basis[mu].type == STO::PX || m_basis[mu].type == STO::PY || m_basis[mu].type == STO::PZ) ? 1 : 0;
        shell_pop[A][l] += PS(mu, mu);
    }

    // Delta-charges per shell and per atom using refocc
    std::vector<std::map<int, double>> dq(m_atomcount);
    Vector charges = Vector::Zero(m_atomcount);

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        for (auto& [s, pop] : shell_pop[A]) {
            double ref = 0.0;
            if (m_param_db.hasElement(Z_A) && m_param_db.getElement(Z_A).shells.count(s))
                ref = m_param_db.getElement(Z_A).shells.at(s).refocc;
            dq[A][s] = ref - pop;
            charges(A) += dq[A][s];
        }
    }

    const_cast<GFN1*>(this)->m_charges = charges;

    // ES2: Effective Coulomb using shell-resolved delta-charges
    // GFN1 uses atom-resolved (not shell-resolved) Hubbard kernel
    double E_ES2 = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        double gamma_AA = m_param_db.hasElement(Z_A)
                          ? m_param_db.getElement(Z_A).gamma_ss
                          : m_params.getHardness(Z_A);

        E_ES2 += 0.5 * charges(A) * charges(A) * gamma_AA;

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];
            double gamma_BB = m_param_db.hasElement(Z_B)
                              ? m_param_db.getElement(Z_B).gamma_ss
                              : m_params.getHardness(Z_B);

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz) / au;

            // Claude Generated (March 2026): Fixed Coulomb kernel — must match GFN2::calculateCoulombKernel
            // Old formula: 1/sqrt(R² + 0.5*(γA+γB)) — wrong: on-site gives 1/sqrt(γ), not γ
            // Correct Klopman-Ohno: 1/sqrt(R² + (2/(γA+γB))²) — on-site gives (γA+γB)/2 = γ ✓
            double gamma_AB = 1.0 / std::sqrt(R_AB*R_AB + std::pow(2.0 / (gamma_AA + gamma_BB), 2.0));

            E_ES2 += charges(A) * charges(B) * gamma_AB;
        }
    }

    return E_ES2;
}

double GFN1::calculateDispersionEnergy() const
{
    // Claude Generated: Native D3(BJ) dispersion for GFN1-xTB
    // Reference: Grimme et al. JCTC 2017, 13, 1989
    // Parameters: s6=1.0, s8=2.4, a1=0.63, a2=5.0 (from TBLite gfn1-xtb.toml)
    //
    // Uses the native D3ParameterGenerator from GFN-FF infrastructure.
    // D3 takes geometry in Angstrom, computes BJ-damped dispersion energy in Eh.
    // GFN1 uses numerical gradients, so D3 gradient is automatically included
    // via finite differences — no analytical D3 gradient needed.

    // Lazy initialization of D3 calculator with GFN1 parameters
    if (!m_d3) {
        m_d3 = std::make_unique<D3ParameterGenerator>(D3ParameterGenerator::createForGFN1());
    }

    // Generate D3 parameters (geometry is in Angstrom, as D3 expects)
    m_d3->GenerateParameters(m_atoms, m_geometry);

    double disp_energy = m_d3->getTotalEnergy();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dispersion_d3", fmt::format("{:.6f} Eh", disp_energy));
    }

    return disp_energy;
}

double GFN1::calculateHalogenBondCorrection() const
{
    // Claude Generated: GFN1 halogen bond correction with TBLite parameters
    // Reference: Grimme et al. JCTC 2017, 13, 1989 (Section 2.5)
    // Updated November 2025: Uses real xb_radius and xb_strength parameters
    //
    // Halogen bonding: X···A interaction (X = F, Cl, Br, I; A = N, O, S, etc.)
    // Empirical correction based on geometry and element-specific parameters

    double E_XB = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];

        if (!isHalogen(Z_A)) continue;

        // Get halogen bond parameters for this halogen
        double xb_radius = 0.0;
        double xb_strength = 0.0;
        bool has_xb_params = false;

        if (m_param_db.hasElement(Z_A)) {
            const auto& elem_A = m_param_db.getElement(Z_A);
            if (elem_A.xb_strength > 1.0e-6) {  // Has XB parameters
                xb_radius = elem_A.xb_radius;
                xb_strength = elem_A.xb_strength;
                has_xb_params = true;
            }
        }

        if (!has_xb_params) continue;  // No XB parameters for this element

        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;

            int Z_B = m_atoms[B];

            // Acceptor atoms (N, O, S, P, etc.)
            bool is_acceptor = (Z_B == 7 || Z_B == 8 || Z_B == 15 || Z_B == 16);

            if (!is_acceptor) continue;

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Halogen bond correction with real parameters
            // E_XB = -strength * f(R) where f(R) is distance-dependent damping
            double R_cutoff = xb_radius + getCovalentRadius(Z_B);

            if (R_AB < R_cutoff * 1.3) {
                // Convert strength from kcal/mol to Hartree
                double k_XB = xb_strength / 627.509;  // kcal/mol → Hartree

                // Distance damping function
                double damp = std::exp(-2.0 * (R_AB / R_cutoff - 1.0));

                E_XB -= k_XB * damp;

                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::param(fmt::format("XB[{}-{}]", Z_A, Z_B),
                                       fmt::format("{:.4f} Å, {:.6f} Eh", R_AB, -k_XB * damp));
                }
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2 && std::abs(E_XB) > 1.0e-6) {
        CurcumaLogger::param("E_halogen_bond_total", fmt::format("{:.6f} Eh", E_XB));
    }

    return E_XB;
}

// =================================================================================
// Gradient Calculation (numerical, same as GFN2)
// =================================================================================

// Claude Generated (March 2026): TBLite atomic radii for shpoly distance scaling (Bohr)
// Source: tblite/data/atomicrad.f90 (Mantina et al., CRC Handbook 2010), converted aatoau
namespace {
    const double aatoau = 1.0 / 0.529177249;
    const double GFN1_ATOMIC_RAD[87] = {
        0.0,  // dummy
        0.32*aatoau, 0.37*aatoau, 1.30*aatoau, 0.99*aatoau, 0.84*aatoau, 0.75*aatoau,  // H-C
        0.71*aatoau, 0.64*aatoau, 0.60*aatoau, 0.62*aatoau,  // N-Ne
        1.60*aatoau, 1.40*aatoau, 1.24*aatoau, 1.14*aatoau, 1.09*aatoau, 1.04*aatoau,  // Na-S
        1.00*aatoau, 1.01*aatoau,  // Cl, Ar
        2.00*aatoau, 1.74*aatoau, 1.59*aatoau, 1.48*aatoau, 1.44*aatoau, 1.30*aatoau,  // K-Cr
        1.29*aatoau, 1.24*aatoau, 1.18*aatoau, 1.17*aatoau, 1.22*aatoau, 1.20*aatoau,  // Mn-Zn
        1.23*aatoau, 1.20*aatoau, 1.20*aatoau, 1.18*aatoau, 1.17*aatoau, 1.16*aatoau,  // Ga-Kr
        2.15*aatoau, 1.90*aatoau, 1.76*aatoau, 1.64*aatoau, 1.56*aatoau, 1.46*aatoau,  // Rb-Mo
        1.38*aatoau, 1.36*aatoau, 1.34*aatoau, 1.30*aatoau, 1.36*aatoau, 1.40*aatoau,  // Tc-Cd
        1.42*aatoau, 1.40*aatoau, 1.40*aatoau, 1.37*aatoau, 1.36*aatoau, 1.36*aatoau,  // In-Xe
        2.38*aatoau, 2.06*aatoau, 1.94*aatoau, 1.84*aatoau, 1.90*aatoau, 1.88*aatoau,  // Cs-Nd
        1.86*aatoau, 1.85*aatoau, 1.83*aatoau, 1.82*aatoau, 1.81*aatoau, 1.80*aatoau,  // Pm-Dy
        1.79*aatoau, 1.77*aatoau, 1.77*aatoau, 1.78*aatoau, 1.74*aatoau, 1.64*aatoau,  // Ho-Hf
        1.58*aatoau, 1.50*aatoau, 1.41*aatoau, 1.36*aatoau, 1.32*aatoau, 1.30*aatoau,  // Ta-Ir
        1.30*aatoau, 1.32*aatoau, 1.44*aatoau, 1.45*aatoau, 1.50*aatoau, 1.42*aatoau,  // Pt-Po
        1.48*aatoau, 1.46*aatoau  // At, Rn
    };

    // GFN1 shell polynomials (p_shpoly * 0.01) from TBLite gfn1.f90 lines 248-292
    // Indexed as [Z][shell], Z=1..86, shell=0(s),1(p),2(d)
    // Only first 18 elements shown explicitly; rest can be added as needed
    const double GFN1_SHPOLY[87][3] = {
        {0.0, 0.0, 0.0},  // dummy Z=0
        { 0.000000*0.01,  0.000000*0.01,  0.000000*0.01},  // H
        { 8.084149*0.01,  0.000000*0.01,  0.000000*0.01},  // He
        {-4.102845*0.01,  9.259276*0.01,  0.000000*0.01},  // Li
        {-12.991482*0.01, -1.308797*0.01,  0.000000*0.01},  // Be
        {-7.088823*0.01,  0.655877*0.01,  0.000000*0.01},  // B
        {-7.082170*0.01,  0.812216*0.01,  0.000000*0.01},  // C
        {-12.745585*0.01, -1.428367*0.01,  0.000000*0.01},  // N
        {-13.729047*0.01, -4.453341*0.01,  0.000000*0.01},  // O
        {-3.921613*0.01, -11.422491*0.01,  0.000000*0.01},  // F
        {-2.115896*0.01, -15.124326*0.01,  0.000000*0.01},  // Ne
        {13.188489*0.01,  10.969376*0.01,  0.000000*0.01},  // Na
        {-19.219408*0.01, 18.272922*0.01,  0.000000*0.01},  // Mg
        {-21.085827*0.01, 24.805127*0.01, 26.405814*0.01},  // Al
        {-14.201582*0.01, -3.893343*0.01, 25.499221*0.01},  // Si
        {-16.118985*0.01, -2.241189*0.01, 30.984577*0.01},  // P
        {-16.989922*0.01, -6.067779*0.01, 16.248395*0.01},  // S
        {-9.341919*0.01,  -8.499805*0.01, 13.088867*0.01},  // Cl
        {-0.082808*0.01,  -9.217948*0.01, 12.204172*0.01},  // Ar
        {12.482844*0.01,  22.323655*0.01,  0.000000*0.01},  // K
        {-11.421376*0.01, 14.628284*0.01, 10.129602*0.01},  // Ca
        { 9.522966*0.01,  44.183320*0.01,-36.027863*0.01},  // Sc
        {24.879987*0.01,  18.910954*0.01,-24.908650*0.01},  // Ti
        {-5.301066*0.01,  22.945047*0.01,-29.197847*0.01},  // V
        {-2.432193*0.01,  11.274054*0.01,-22.608167*0.01},  // Cr
        { 1.025345*0.01,   1.834626*0.01,-25.016650*0.01},  // Mn
        {-2.182723*0.01,  11.769535*0.01,-22.920815*0.01},  // Fe
        { 0.815250*0.01,  15.765732*0.01,-21.678930*0.01},  // Co
        {15.160508*0.01,  15.782685*0.01,-26.348820*0.01},  // Ni
        {-3.590501*0.01,   7.413473*0.01,-21.142399*0.01},  // Cu
        {-15.535695*0.01,  4.061664*0.01,  0.000000*0.01},  // Zn
        {-14.584657*0.01,  9.375082*0.01, 19.671655*0.01},  // Ga
        {-12.195371*0.01,-11.374296*0.01,  9.364108*0.01},  // Ge
        {-17.489686*0.01, -6.747956*0.01, 17.858510*0.01},  // As
        {-14.852299*0.01, -9.863477*0.01,  9.556181*0.01},  // Se
        {-17.815502*0.01,-14.058044*0.01,  5.468245*0.01},  // Br
        {-25.437273*0.01,-12.813227*0.01, 10.440712*0.01},  // Kr
        {-7.450752*0.01,  16.670533*0.01,  0.000000*0.01},  // Rb
        {-6.087125*0.01,   2.115262*0.01, 17.076466*0.01},  // Sr
        {10.950764*0.01,  45.679760*0.01,-28.061976*0.01},  // Y
        {44.110231*0.01,  25.863572*0.01,-22.240873*0.01},  // Zr
        {15.379439*0.01,  30.159730*0.01,-25.998052*0.01},  // Nb
        { 5.815301*0.01,  14.527159*0.01,-22.556077*0.01},  // Mo
        {24.977603*0.01,   1.953838*0.01,-23.231470*0.01},  // Tc
        {15.281981*0.01,   1.340798*0.01,-23.099524*0.01},  // Ru
        {24.928002*0.01,  -4.330556*0.01,-19.564083*0.01},  // Rh-correction: moved
        {25.774929*0.01,  -0.704597*0.01,-21.172493*0.01},  // Pd
        {38.415536*0.01,  -0.665483*0.01,-22.169385*0.01},  // Ag
        {-11.443658*0.01, -5.119735*0.01,-11.067532*0.01},  // Cd
        {-6.581368*0.01,   3.995243*0.01,  0.000000*0.01},  // In
        {-2.193199*0.01,   0.060451*0.01,  0.000000*0.01},  // Sn
        {-10.874138*0.01, -6.034796*0.01,  0.000000*0.01},  // Sb
        {-20.410234*0.01, -9.424568*0.01,  0.000000*0.01},  // Te
        {-18.477865*0.01,-14.037423*0.01, 13.809093*0.01},  // I
        {-21.965390*0.01,-12.804436*0.01, 16.836546*0.01},  // Xe
        {-22.139701*0.01,-20.539955*0.01, 17.249637*0.01},  // Cs
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Ba-Pm (56-61)
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Sm-Tb (62-65)
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Dy-Tm (66-69)
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Yb-Hf (70-73)
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Ta-Os (74-77)
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Ir-Pb (78-82)
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},  // Bi-Rn (83-86)
    };

    // GFN1 TBLite constants (gfn1.f90 lines 54-58)
    constexpr double GFN1_REP_KEXP = 1.5;
    constexpr double GFN1_REP_REXP = 1.0;
    constexpr double GFN1_ENSCALE = -7.0e-3;
    constexpr double GFN1_KDIAG[] = {1.85, 2.25, 2.0, 2.0, 2.0};
    constexpr double GFN1_KDIFF = 2.85;

    inline double getGFN1AtomicRad(int Z) {
        if (Z >= 1 && Z <= 86) return GFN1_ATOMIC_RAD[Z];
        return 2.0;  // fallback
    }

    inline double getGFN1ShPoly(int Z, int shell) {
        if (Z >= 1 && Z <= 86 && shell >= 0 && shell <= 2) return GFN1_SHPOLY[Z][shell];
        return 0.0;
    }
}

Matrix GFN1::calculateGradient() const
{
    // Claude Generated (March 2026): Analytical gradients for GFN1-xTB
    // Adapted from GFN2 analytical gradient pattern
    // Components: Pulay forces (H0/S derivatives), repulsion, Coulomb (ES2), CN, D3, halogen bond
    // Reference: TBLite xtb/h0.f90 get_hamiltonian_gradient

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Calculating GFN1 analytical gradients");
    }

    Matrix gradient = Matrix::Zero(m_atomcount, 3);

    // Energy-weighted density matrix W = 2 * sum_{i occ} eps_i * C_i * C_i^T
    Matrix W = Matrix::Zero(m_nbasis, m_nbasis);
    int n_occ = m_num_electrons / 2;
    for (int i = 0; i < n_occ; ++i) {
        double eps = m_energies(i);
        for (int mu = 0; mu < m_nbasis; ++mu)
            for (int nu = 0; nu < m_nbasis; ++nu)
                W(mu, nu) += 2.0 * eps * m_mo(mu, i) * m_mo(nu, i);
    }

    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = A + 1; B < m_atomcount; ++B) {
            double dx = m_geometry(B, 0) - m_geometry(A, 0);
            double dy = m_geometry(B, 1) - m_geometry(A, 1);
            double dz = m_geometry(B, 2) - m_geometry(A, 2);
            double R_AB_ang = std::sqrt(dx*dx + dy*dy + dz*dz);
            double R_AB_bohr = R_AB_ang / au;
            if (R_AB_ang < 1e-10) continue;
            double nx = dx / R_AB_ang, ny = dy / R_AB_ang, nz = dz / R_AB_ang;

            int ZA = m_atoms[A], ZB = m_atoms[B];

            // ===== 1. Repulsion gradient =====
            // V_rep = zeff_AB * exp(-alpha_AB * R^kexp) / R^rexp
            // dV/dR = -(alpha*R^kexp*kexp + rexp) * V / R
            double alpha_A = 0, alpha_B = 0, zeff_A = 0, zeff_B = 0;
            if (m_param_db.hasElement(ZA)) { alpha_A = m_param_db.getElement(ZA).rep_alpha; zeff_A = m_param_db.getElement(ZA).rep_zeff; }
            if (m_param_db.hasElement(ZB)) { alpha_B = m_param_db.getElement(ZB).rep_alpha; zeff_B = m_param_db.getElement(ZB).rep_zeff; }
            double alpha_AB = std::sqrt(alpha_A * alpha_B);
            double zeff_AB = zeff_A * zeff_B;
            double kexp = (ZA > 2 || ZB > 2) ? GFN1_REP_KEXP : 1.0;
            double rexp = GFN1_REP_REXP;
            double r1k = std::pow(R_AB_bohr, kexp);
            double r1r = std::pow(R_AB_bohr, rexp);
            double Vrep = zeff_AB * std::exp(-alpha_AB * r1k) / r1r;
            double dVdR = -(alpha_AB * r1k * kexp + rexp) * Vrep / (R_AB_bohr * R_AB_bohr);
            // dVdR is in Hartree/Bohr, project onto Cartesian via unit vector
            // But R_AB_ang = R_AB_bohr * au, so dV/d(R_ang) = dV/d(R_bohr) / au
            gradient(A, 0) += dVdR * nx / au; gradient(A, 1) += dVdR * ny / au; gradient(A, 2) += dVdR * nz / au;
            gradient(B, 0) -= dVdR * nx / au; gradient(B, 1) -= dVdR * ny / au; gradient(B, 2) -= dVdR * nz / au;

            // ===== 2. Hamiltonian (Pulay) forces =====
            // dE_el/dR = sum_{mu,nu} [2*P*dH0/dR - 2*W*dS/dR]
            double rad_A = getGFN1AtomicRad(ZA), rad_B = getGFN1AtomicRad(ZB);

            for (size_t mu = 0; mu < m_basis.size(); ++mu) {
                if (m_basis[mu].atom != A) continue;
                for (size_t nu = 0; nu < m_basis.size(); ++nu) {
                    if (m_basis[nu].atom != B) continue;

                    double dS_dR = STO::calculateOverlapDerivative(m_basis[mu], m_basis[nu], R_AB_bohr);
                    double h_avg = (m_hamiltonian(mu, mu) + m_hamiltonian(nu, nu)) / 2.0;

                    // Shell types for shpoly lookup
                    int shell_mu = (m_basis[mu].type == STO::PX || m_basis[mu].type == STO::PY || m_basis[mu].type == STO::PZ) ? 1 : 0;
                    int shell_nu = (m_basis[nu].type == STO::PX || m_basis[nu].type == STO::PY || m_basis[nu].type == STO::PZ) ? 1 : 0;

                    // Distance-dependent polynomial (TBLite shpoly formulation)
                    double sh_A = getGFN1ShPoly(ZA, shell_mu);
                    double sh_B = getGFN1ShPoly(ZB, shell_nu);
                    double rr = std::sqrt(R_AB_bohr / (rad_A + rad_B));
                    double poly = (1.0 + sh_A * rr) * (1.0 + sh_B * rr);

                    // Derivative: d(poly)/dR = (sh_A * drr * (1+sh_B*rr) + (1+sh_A*rr) * sh_B * drr)
                    // drr = d/dR sqrt(R/(radA+radB)) = 1 / (2*sqrt(R*(radA+radB)))
                    double drr = 1.0 / (2.0 * std::sqrt(R_AB_bohr * (rad_A + rad_B)));
                    double dpoly = sh_A * drr * (1.0 + sh_B * rr) + (1.0 + sh_A * rr) * sh_B * drr;

                    // Orbital overlap scaling
                    double z_ij = std::sqrt(2.0 * std::sqrt(m_basis[mu].zeta * m_basis[nu].zeta) / (m_basis[mu].zeta + m_basis[nu].zeta));

                    // EN scaling (GFN1: enscale = -7.0e-3)
                    double EN_A = m_params.getElectronegativity(ZA);
                    double EN_B = m_params.getElectronegativity(ZB);
                    double dEN = EN_A - EN_B;

                    // Shell coupling from kdiag
                    double k_s = (GFN1_KDIAG[shell_mu] + GFN1_KDIAG[shell_nu]) / 2.0;

                    // Pair coupling
                    double k_p = 1.0;
                    if (m_param_db.hasPair(ZA, ZB)) {
                        const auto& pair = m_param_db.getPair(ZA, ZB);
                        k_p = pair.kpair;
                        if (shell_mu == 0 && shell_nu == 0) k_s = pair.kshell_ss;
                        else if (shell_mu == 1 && shell_nu == 1) k_s = pair.kshell_pp;
                        else k_s = pair.kshell_sp;
                    }

                    double C = z_ij * k_p * k_s * (1.0 + GFN1_ENSCALE * dEN * dEN);
                    double dH0 = (C * dpoly * m_overlap(mu, nu) + C * poly * dS_dR) * h_avg;
                    double contrib = 2.0 * (m_density(mu, nu) * dH0 - W(mu, nu) * dS_dR);

                    gradient(A, 0) += (contrib / au) * nx; gradient(A, 1) += (contrib / au) * ny; gradient(A, 2) += (contrib / au) * nz;
                    gradient(B, 0) -= (contrib / au) * nx; gradient(B, 1) -= (contrib / au) * ny; gradient(B, 2) -= (contrib / au) * nz;
                }
            }

            // ===== 3. Coulomb (ES2) gradient =====
            // gamma_AB = 1/sqrt(R^2 + 0.5*(gamma_AA + gamma_BB))
            // dgamma/dR = -R / (R^2 + eta)^(3/2) where eta = 0.5*(gAA+gBB)
            double gamma_AA = m_param_db.hasElement(ZA) ? m_param_db.getElement(ZA).gamma_ss : m_params.getHardness(ZA);
            double gamma_BB = m_param_db.hasElement(ZB) ? m_param_db.getElement(ZB).gamma_ss : m_params.getHardness(ZB);
            double eta = 0.5 * (gamma_AA + gamma_BB);
            double dgamma = -R_AB_bohr / std::pow(R_AB_bohr * R_AB_bohr + eta, 1.5);
            double dES2 = (m_charges(A) * m_charges(B) * dgamma) / au;
            gradient(A, 0) += dES2 * nx; gradient(A, 1) += dES2 * ny; gradient(A, 2) += dES2 * nz;
            gradient(B, 0) -= dES2 * nx; gradient(B, 1) -= dES2 * ny; gradient(B, 2) -= dES2 * nz;

            // ===== 4. CN gradient =====
            // CN_A depends on R_AB, self-energy depends on CN: dE/dR = dE/dCN * dCN/dR
            double R_cov = getCovalentRadius(ZA) + getCovalentRadius(ZB);
            double exp_val = std::exp(-GFN1_K1 * (R_cov / R_AB_ang - 1.0));
            double count = 1.0 / (1.0 + exp_val);
            // d(count)/dR = K1 * R_cov * exp / (R^2 * (1+exp)^2)
            // d(count^K2)/dR = K2 * count^(K2-1) * d(count)/dR
            double dcount_dR = (-GFN1_K1 * R_cov * exp_val) / (R_AB_ang * R_AB_ang * (1.0 + exp_val) * (1.0 + exp_val));
            double dCN_dR = GFN1_K2 * std::pow(count, GFN1_K2 - 1.0) * dcount_dR;

            // dE/dCN_A = sum_{mu on A} kcn_A_shell * P(mu,mu)
            // dE/dCN_B = sum_{mu on B} kcn_B_shell * P(mu,mu)
            double dE_dCN_A = 0, dE_dCN_B = 0;
            for (size_t mu = 0; mu < m_basis.size(); ++mu) {
                if (m_basis[mu].atom == A) {
                    int shell = (m_basis[mu].type == STO::PX || m_basis[mu].type == STO::PY || m_basis[mu].type == STO::PZ) ? 1 : 0;
                    double kcn = 0;
                    if (m_param_db.hasElement(ZA) && m_param_db.getElement(ZA).shells.count(shell))
                        kcn = m_param_db.getElement(ZA).shells.at(shell).kcn;
                    dE_dCN_A += kcn * m_density(mu, mu);
                }
                if (m_basis[mu].atom == B) {
                    int shell = (m_basis[mu].type == STO::PX || m_basis[mu].type == STO::PY || m_basis[mu].type == STO::PZ) ? 1 : 0;
                    double kcn = 0;
                    if (m_param_db.hasElement(ZB) && m_param_db.getElement(ZB).shells.count(shell))
                        kcn = m_param_db.getElement(ZB).shells.at(shell).kcn;
                    dE_dCN_B += kcn * m_density(mu, mu);
                }
            }
            double dECN = (dE_dCN_A + dE_dCN_B) * dCN_dR;
            gradient(A, 0) += dECN * nx; gradient(A, 1) += dECN * ny; gradient(A, 2) += dECN * nz;
            gradient(B, 0) -= dECN * nx; gradient(B, 1) -= dECN * ny; gradient(B, 2) -= dECN * nz;

            // ===== 5. Halogen bond gradient =====
            if (isHalogen(ZA) || isHalogen(ZB)) {
                int Z_hal = isHalogen(ZA) ? ZA : ZB;
                int Z_acc = isHalogen(ZA) ? ZB : ZA;
                int idx_hal = isHalogen(ZA) ? A : B;
                int idx_acc = isHalogen(ZA) ? B : A;
                bool is_acceptor = (Z_acc == 7 || Z_acc == 8 || Z_acc == 15 || Z_acc == 16);

                if (is_acceptor && m_param_db.hasElement(Z_hal)) {
                    const auto& elem = m_param_db.getElement(Z_hal);
                    if (elem.xb_strength > 1.0e-6) {
                        double R_cutoff = elem.xb_radius + getCovalentRadius(Z_acc);
                        if (R_AB_ang < R_cutoff * 1.3) {
                            double k_XB = elem.xb_strength / 627.509;
                            double damp = std::exp(-2.0 * (R_AB_ang / R_cutoff - 1.0));
                            // dE_XB/dR = k_XB * (2/R_cutoff) * damp
                            double dXB = k_XB * (2.0 / R_cutoff) * damp;
                            gradient(idx_hal, 0) += dXB * nx; gradient(idx_hal, 1) += dXB * ny; gradient(idx_hal, 2) += dXB * nz;
                            gradient(idx_acc, 0) -= dXB * nx; gradient(idx_acc, 1) -= dXB * ny; gradient(idx_acc, 2) -= dXB * nz;
                        }
                    }
                }
            }
        }
    }

    // ===== 6. D3 dispersion gradient (via D3ParameterGenerator) =====
    // D3 gradient is computed by finite differences on D3 energy
    // (D3ParameterGenerator doesn't have analytical gradient API yet)
    if (m_d3) {
        const double delta = 1.0e-5;  // Angstrom
        for (int atom = 0; atom < m_atomcount; ++atom) {
            for (int coord = 0; coord < 3; ++coord) {
                Matrix geom_p = m_geometry, geom_m = m_geometry;
                geom_p(atom, coord) += delta;
                geom_m(atom, coord) -= delta;
                D3ParameterGenerator d3p(D3ParameterGenerator::createForGFN1());
                D3ParameterGenerator d3m(D3ParameterGenerator::createForGFN1());
                d3p.GenerateParameters(m_atoms, geom_p);
                d3m.GenerateParameters(m_atoms, geom_m);
                gradient(atom, coord) += (d3p.getTotalEnergy() - d3m.getTotalEnergy()) / (2.0 * delta);
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("analytical_gradient_norm", fmt::format("{:.6e} Eh/Å", gradient.norm()));
        for (int i = 0; i < m_atomcount; ++i)
            CurcumaLogger::info(fmt::format("Atom {}: {:10.6f} {:10.6f} {:10.6f}", i+1, gradient(i,0), gradient(i,1), gradient(i,2)));
    }

    return gradient;
}

Matrix GFN1::calculateNumericalGradient(double delta) const
{
    // Claude Generated (March 2026): Numerical gradient for validation
    Matrix gradient = Matrix::Zero(m_atomcount, 3);
    Matrix geom_orig = m_geometry;
    double E0 = m_total_energy;

    for (int atom = 0; atom < m_atomcount; ++atom) {
        for (int coord = 0; coord < 3; ++coord) {
            const_cast<GFN1*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) + delta;
            const_cast<GFN1*>(this)->InitialiseMolecule();
            double E_plus = const_cast<GFN1*>(this)->Calculation(false);

            const_cast<GFN1*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) - delta;
            const_cast<GFN1*>(this)->InitialiseMolecule();
            double E_minus = const_cast<GFN1*>(this)->Calculation(false);

            gradient(atom, coord) = (E_plus - E_minus) / (2.0 * delta);
            const_cast<GFN1*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord);
        }
    }

    const_cast<GFN1*>(this)->m_geometry = geom_orig;
    const_cast<GFN1*>(this)->InitialiseMolecule();
    const_cast<GFN1*>(this)->m_total_energy = E0;

    return gradient / au;
}

// =================================================================================
// Property Access
// =================================================================================

json GFN1::getEnergyDecomposition() const
{
    json decomp;
    decomp["electronic"] = m_energy_electronic;
    decomp["repulsion"] = m_energy_repulsion;
    decomp["coulomb"] = m_energy_coulomb;
    decomp["dispersion"] = m_energy_dispersion;
    decomp["halogen_bond"] = m_energy_halogen_bond;
    decomp["total"] = m_total_energy;

    return decomp;
}

double GFN1::getHOMOLUMOGap() const
{
    return (getLUMOEnergy() - getHOMOEnergy());
}

double GFN1::getHOMOEnergy() const
{
    int n_occ = m_num_electrons / 2;
    if (n_occ == 0 || n_occ > m_nbasis) return 0.0;

    return m_energies(n_occ - 1) * eV2Eh;
}

double GFN1::getLUMOEnergy() const
{
    int n_occ = m_num_electrons / 2;
    if (n_occ >= m_nbasis) return 0.0;

    return m_energies(n_occ) * eV2Eh;
}
