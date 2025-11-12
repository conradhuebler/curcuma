/*
 * <Native GFN2-xTB Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the GFN2-xTB method developed by:
 *   Stefan Grimme, Christoph Bannwarth, Sebastian Ehlert
 *   Mulliken Center for Theoretical Chemistry, University of Bonn
 *
 * Reference implementation: TBLite (https://github.com/tblite/tblite)
 *   See: src/tblite/xtb/gfn2.f90
 *
 * Original publication:
 *   C. Bannwarth, S. Ehlert, S. Grimme
 *   J. Chem. Theory Comput. 2019, 15, 1652-1671
 *
 * This program is free software under GPL-3.0
 */

#include "gfn2.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "ParallelEigenSolver.hpp"

#include <fmt/format.h>
#include <cmath>
#include <algorithm>

using namespace CurcumaUnit;

// =================================================================================
// Constructor / Destructor
// =================================================================================

GFN2::GFN2()
    : m_params()
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_repulsion(0.0)
    , m_energy_coulomb(0.0)
    , m_energy_dispersion(0.0)
    , m_scf_max_iterations(100)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing native GFN2-xTB method");
        CurcumaLogger::param("method", "GFN2-xTB");
        CurcumaLogger::param("reference", "Bannwarth et al. JCTC 2019, 15, 1652");
    }
}

GFN2::GFN2(const ArrayParameters& params)
    : m_params(params)
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_repulsion(0.0)
    , m_energy_coulomb(0.0)
    , m_energy_dispersion(0.0)
    , m_scf_max_iterations(100)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing GFN2-xTB with custom parameters");
    }
}

// =================================================================================
// QMDriver Interface Implementation
// =================================================================================

bool GFN2::InitialiseMolecule()
{
    if (m_atoms.size() == 0) {
        CurcumaLogger::error("No atoms in molecule for GFN2 initialization");
        return false;
    }

    // Validate all atoms are supported
    for (size_t i = 0; i < m_atoms.size(); ++i) {
        if (!m_params.isValidAtom(m_atoms[i])) {
            CurcumaLogger::error(fmt::format("GFN2 parameters not available for element Z={} at atom {}",
                                            m_atoms[i], i+1));
            return false;
        }
    }

    // Level 1+: Initialization info
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing GFN2 calculation");
        CurcumaLogger::param("atoms", static_cast<int>(m_atoms.size()));
        CurcumaLogger::param("charge", m_charge);
        CurcumaLogger::param("spin", m_spin);
    }

    // Build basis set
    m_nbasis = buildBasisSet();

    if (m_nbasis == 0) {
        CurcumaLogger::error("Failed to build basis set");
        return false;
    }

    // Level 2+: Basis set details
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("basis_functions", m_nbasis);
        CurcumaLogger::param("electrons", m_num_electrons);
    }

    // Calculate coordination numbers
    m_coordination_numbers = calculateCoordinationNumbers();

    // Level 3+: CN details
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Coordination numbers calculated");
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            CurcumaLogger::param(fmt::format("CN_atom_{}", i+1),
                               fmt::format("{:.3f}", m_coordination_numbers(i)));
        }
    }

    // Initialize matrices
    m_overlap = Matrix::Zero(m_nbasis, m_nbasis);
    m_hamiltonian = Matrix::Zero(m_nbasis, m_nbasis);
    m_fock = Matrix::Zero(m_nbasis, m_nbasis);
    m_density = Matrix::Zero(m_nbasis, m_nbasis);
    m_gradient = Matrix::Zero(m_atomcount, 3);
    m_mo = Matrix::Zero(m_nbasis, m_nbasis);
    m_energies = Vector::Zero(m_nbasis);
    m_charges = Vector::Zero(m_atomcount);

    return true;
}

double GFN2::Calculation(bool gradient)
{
    if (!m_initialised && m_atoms.size() == 0) {
        CurcumaLogger::error("Molecule not initialized for GFN2 calculation");
        return 0.0;
    }

    // Level 1+: Starting calculation
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Starting GFN2-xTB calculation");
    }

    // Gradient warning
    if (gradient) {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::warn("GFN2 gradients not yet implemented - results will be zero");
        }
    }

    try {
        // Step 1: Build overlap matrix
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 1: Computing overlap matrix");
        }
        m_overlap = MakeOverlap(m_basis);

        // Step 2: Build core Hamiltonian
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 2: Building core Hamiltonian");
        }
        m_hamiltonian = MakeH(m_overlap, m_basis);

        // Step 3: SCF convergence
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 3: Running SCF convergence");
        }
        m_scf_converged = runSCF();

        if (!m_scf_converged) {
            CurcumaLogger::warn("SCF did not converge within maximum iterations");
        }

        // Step 4: Calculate energy components
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 4: Calculating energy components");
        }

        m_energy_electronic = calculateElectronicEnergy();
        m_energy_repulsion = calculateRepulsionEnergy();
        m_energy_coulomb = calculateCoulombEnergy();
        m_energy_dispersion = calculateDispersionEnergy();  // Stub: returns 0

        m_total_energy = m_energy_electronic + m_energy_repulsion +
                        m_energy_coulomb + m_energy_dispersion;

        // Level 1+: Final results
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::energy_abs(m_total_energy, "GFN2 Total Energy");
        }

        // Level 2+: Energy decomposition
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", m_energy_electronic));
            CurcumaLogger::param("repulsion", fmt::format("{:.6f} Eh", m_energy_repulsion));
            CurcumaLogger::param("coulomb", fmt::format("{:.6f} Eh", m_energy_coulomb));
            CurcumaLogger::param("dispersion", fmt::format("{:.6f} Eh (stub)", m_energy_dispersion));

            // Orbital energies
            double homo = getHOMOEnergy();
            double lumo = getLUMOEnergy();
            double gap = getHOMOLUMOGap();

            CurcumaLogger::param("HOMO", fmt::format("{:.4f} eV", homo));
            CurcumaLogger::param("LUMO", fmt::format("{:.4f} eV", lumo));
            CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} eV", gap));
        }

        // Level 3+: Full details
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("SCF_converged", m_scf_converged ? "yes" : "no");
            CurcumaLogger::param("basis_size", m_nbasis);
        }

        // Gradient calculation (if requested and implemented)
        if (gradient) {
            m_gradient = calculateGradient();
        }

        return m_total_energy;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("GFN2 calculation failed: {}", e.what()));
        return 0.0;
    }
}

// =================================================================================
// Basis Set Construction
// =================================================================================

int GFN2::buildBasisSet()
{
    m_basis.clear();
    m_num_electrons = 0;

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];

        // Get number of shells for this element
        // GFN2 uses minimal valence basis:
        // - Period 1 (H, He): s-shell only
        // - Period 2 (Li-Ne): s, p shells
        // - Period 3+: s, p, d shells (for transition metals)

        int n_shells = 1;  // At least s-shell
        if (Z > 2) n_shells = 2;   // s, p shells
        if (Z > 18) n_shells = 3;  // s, p, d shells (simplified)

        for (int shell = 0; shell < n_shells; ++shell) {
            // Get parameters for this shell
            double zeta = m_params.getYeff(Z, shell);
            int principal_qn = static_cast<int>(m_params.getPrincipalQN(Z, shell));

            if (zeta < 1.0e-6) continue;  // Skip unparameterized shells

            // Angular momentum: shell 0 = s, shell 1 = p, shell 2 = d
            int l = shell;

            // Create basis functions for this shell
            for (int m = -l; m <= l; ++m) {
                STO::Orbital orbital;
                orbital.type = (l == 0) ? STO::S : ((l == 1) ? STO::PX : STO::S);  // Simplified
                orbital.x = m_geometry(i, 0) / au;  // Convert Å → Bohr
                orbital.y = m_geometry(i, 1) / au;
                orbital.z = m_geometry(i, 2) / au;
                orbital.zeta = zeta;
                orbital.VSIP = 0.0;  // Not used in GFN2 (uses selfenergy table)
                orbital.atom = i;

                m_basis.push_back(orbital);
            }

            // Add electrons for this shell (simplified)
            if (shell == 0) m_num_electrons += std::min(2, Z);              // s: max 2e
            else if (shell == 1) m_num_electrons += std::min(6, Z - 2);     // p: max 6e
            else if (shell == 2) m_num_electrons += std::min(10, Z - 10);   // d: max 10e
        }
    }

    return static_cast<int>(m_basis.size());
}

// =================================================================================
// Overlap Matrix
// =================================================================================

Matrix GFN2::MakeOverlap(const std::vector<STO::Orbital>& basisset)
{
    Matrix S = Matrix::Zero(basisset.size(), basisset.size());

    // Update basis function positions
    for (size_t i = 0; i < basisset.size(); ++i) {
        // Positions already set in buildBasisSet
    }

    // Calculate overlap integrals
    for (size_t i = 0; i < basisset.size(); ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double overlap = STO::calculateOverlap(basisset[i], basisset[j]);
            S(i, j) = S(j, i) = overlap;
        }
    }

    return S;
}

// =================================================================================
// Hamiltonian Matrix
// =================================================================================

Matrix GFN2::MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset)
{
    Matrix H = Matrix::Zero(basisset.size(), basisset.size());

    // Stub implementation - simplified for now
    // TODO: Implement full GFN2 Hamiltonian construction

    for (size_t i = 0; i < basisset.size(); ++i) {
        for (size_t j = 0; j <= i; ++j) {
            if (i == j) {
                // Diagonal: self-energy
                int atom_i = basisset[i].atom;
                int Z_i = m_atoms[atom_i];
                int shell_i = 0;  // Simplified: determine from orbital type
                double CN_i = m_coordination_numbers(atom_i);

                H(i, j) = getSelfEnergy(Z_i, shell_i, CN_i);
            } else {
                // Off-diagonal: hopping integral (stub)
                H(i, j) = H(j, i) = 0.0;  // TODO: implement getHamiltonianScale
            }
        }
    }

    return H;
}

double GFN2::getSelfEnergy(int element, int shell, double CN) const
{
    // GFN2 Eq. 12a: E_ii = E_base + k_CN * CN
    // Note: ArrayParameters stores energies in eV, convert to Hartree

    // Stub: simplified implementation
    // TODO: Access correct shell-resolved parameters from m_params

    double E_base = -10.0;  // eV (placeholder)
    double k_CN = 0.1;      // eV (placeholder)

    double E_hartree = (E_base + k_CN * CN) / eV2Eh;

    return E_hartree;
}

double GFN2::getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const
{
    // TODO: Implement GFN2 Eq. 12b
    // z_ij * k_pair * k_shell * (1 + 0.02 * ΔEN²) * poly(r) * S_ij

    return 0.0;  // Stub
}

// =================================================================================
// Coordination Numbers
// =================================================================================

Vector GFN2::calculateCoordinationNumbers()
{
    Vector CN = Vector::Zero(m_atomcount);

    // GFN2 parameters (from paper)
    const double k1 = 16.0;      // Steepness
    const double k2 = 4.0 / 3.0; // Range decay exponent

    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;

            // Distance in Ångström
            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Covalent radii (placeholder - should use m_params)
            double R_cov_A = 1.0;  // Å (TODO: get from parameters)
            double R_cov_B = 1.0;  // Å
            double R_cov = R_cov_A + R_cov_B;

            // Counting function: 1 / (1 + exp(-k1 * (R_cov/R - 1)))
            double count = 1.0 / (1.0 + std::exp(-k1 * (R_cov / R_AB - 1.0)));

            // Power k2
            CN(A) += std::pow(count, k2);
        }
    }

    return CN;
}

// =================================================================================
// SCF Procedure
// =================================================================================

bool GFN2::runSCF()
{
    // Initial guess: zero density (simplified)
    // TODO: Implement SAD (Superposition of Atomic Densities) guess
    m_density = Matrix::Zero(m_nbasis, m_nbasis);

    for (int iter = 0; iter < m_scf_max_iterations; ++iter) {
        // 1. Build Fock matrix
        m_fock = buildFockMatrix(m_density);

        // 2. Solve generalized eigenvalue problem: FC = SCE
        ParallelEigenSolver solver(500, 128, 1.0e-10, false);
        solver.setThreadCount(m_threads);

        bool success = solver.solve(m_overlap, m_fock, m_energies, m_mo, m_threads, false);

        if (!success) {
            CurcumaLogger::error("Eigenvalue solution failed in SCF");
            return false;
        }

        // 3. Build new density matrix
        Matrix density_new = buildDensityMatrix(m_mo, m_energies);

        // 4. Check convergence
        double delta_P = (density_new - m_density).norm();

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("SCF_iter_{}", iter+1),
                               fmt::format("ΔP = {:.6e}", delta_P));
        }

        if (delta_P < m_scf_threshold) {
            m_density = density_new;
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::success(fmt::format("SCF converged in {} iterations", iter+1));
            }
            return true;
        }

        // 5. Apply damping
        m_density = m_scf_damping * density_new + (1.0 - m_scf_damping) * m_density;
    }

    return false;  // Did not converge
}

Matrix GFN2::buildFockMatrix(const Matrix& density)
{
    // Simplified: F = H₀ (no Coulomb terms yet)
    // TODO: Add G(P) Coulomb contributions

    return m_hamiltonian;  // Stub
}

Matrix GFN2::buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies)
{
    Matrix P = Matrix::Zero(m_nbasis, m_nbasis);

    int n_occ = m_num_electrons / 2;  // Doubly occupied orbitals

    for (int mu = 0; mu < m_nbasis; ++mu) {
        for (int nu = 0; nu < m_nbasis; ++nu) {
            for (int i = 0; i < n_occ; ++i) {
                P(mu, nu) += 2.0 * mo_coefficients(mu, i) * mo_coefficients(nu, i);
            }
        }
    }

    return P;
}

// =================================================================================
// Energy Components
// =================================================================================

double GFN2::calculateElectronicEnergy() const
{
    // E_elec = Tr(P * (H + F)) / 2
    // Simplified for now

    double E = 0.0;

    for (int mu = 0; mu < m_nbasis; ++mu) {
        for (int nu = 0; nu < m_nbasis; ++nu) {
            E += m_density(mu, nu) * (m_hamiltonian(mu, nu) + m_fock(mu, nu));
        }
    }

    return E / 2.0;
}

double GFN2::calculateRepulsionEnergy() const
{
    // TODO: Implement pairwise exponential repulsion
    // E_rep = ∑_{A<B} Z_eff,A * exp(-α_A * R_AB)

    return 0.0;  // Stub
}

double GFN2::calculateCoulombEnergy() const
{
    // TODO: Implement ES2 + ES3 + AES2
    // See GFN2 paper Eq. 7-11

    return 0.0;  // Stub
}

double GFN2::calculateDispersionEnergy() const
{
    // D4 dispersion stub (separate TODO)
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::warn("D4 dispersion not yet integrated - energy incomplete");
    }

    return 0.0;
}

// =================================================================================
// Gradient Calculation
// =================================================================================

Matrix GFN2::calculateGradient() const
{
    // TODO: Implement analytical gradients
    // Uses Hellmann-Feynman theorem

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::warn("GFN2 gradients not yet implemented");
    }

    return Matrix::Zero(m_atomcount, 3);
}

// =================================================================================
// Property Access
// =================================================================================

json GFN2::getEnergyDecomposition() const
{
    json decomp;
    decomp["electronic"] = m_energy_electronic;
    decomp["repulsion"] = m_energy_repulsion;
    decomp["coulomb"] = m_energy_coulomb;
    decomp["dispersion"] = m_energy_dispersion;
    decomp["total"] = m_total_energy;
    decomp["scf_converged"] = m_scf_converged;

    return decomp;
}

double GFN2::getHOMOLUMOGap() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) return 0.0;

    int homo_index = m_num_electrons / 2 - 1;
    int lumo_index = homo_index + 1;

    if (homo_index < 0 || lumo_index >= m_energies.size()) return 0.0;

    return (m_energies(lumo_index) - m_energies(homo_index)) * eV2Eh;  // Convert to eV
}

double GFN2::getHOMOEnergy() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) return 0.0;

    int homo_index = m_num_electrons / 2 - 1;
    if (homo_index < 0 || homo_index >= m_energies.size()) return 0.0;

    return m_energies(homo_index) * eV2Eh;
}

double GFN2::getLUMOEnergy() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) return 0.0;

    int lumo_index = m_num_electrons / 2;
    if (lumo_index >= m_energies.size()) return 0.0;

    return m_energies(lumo_index) * eV2Eh;
}
