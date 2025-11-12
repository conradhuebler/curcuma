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
// Covalent Radii (Pyykkö 2015, triple bond radii in Ångström)
// Shared by GFN2, GFN1, and other xTB methods
// =================================================================================
// Claude Generated: Pyykkö triple-bond covalent radii for xTB coordination numbers
// Reference: P. Pyykkö, J. Phys. Chem. A 2015, 119, 2326-2337
namespace {
    const double COVALENT_RADII[87] = {
        0.0,    // Dummy for index 0
        0.32,   // H
        0.46,   // He
        1.33,   // Li
        1.02,   // Be
        0.85,   // B
        0.75,   // C
        0.71,   // N
        0.63,   // O
        0.64,   // F
        0.67,   // Ne
        1.55,   // Na
        1.39,   // Mg
        1.26,   // Al
        1.16,   // Si
        1.11,   // P
        1.03,   // S
        0.99,   // Cl
        0.96,   // Ar
        1.96,   // K
        1.71,   // Ca
        1.48,   // Sc
        1.36,   // Ti
        1.34,   // V
        1.22,   // Cr
        1.19,   // Mn
        1.16,   // Fe
        1.11,   // Co
        1.10,   // Ni
        1.12,   // Cu
        1.18,   // Zn
        1.24,   // Ga
        1.21,   // Ge
        1.21,   // As
        1.16,   // Se
        1.14,   // Br
        1.17,   // Kr
        2.10,   // Rb
        1.85,   // Sr
        1.63,   // Y
        1.54,   // Zr
        1.47,   // Nb
        1.38,   // Mo
        1.28,   // Tc
        1.25,   // Ru
        1.25,   // Rh
        1.20,   // Pd
        1.28,   // Ag
        1.36,   // Cd
        1.42,   // In
        1.40,   // Sn
        1.40,   // Sb
        1.36,   // Te
        1.33,   // I
        1.31,   // Xe
        2.32,   // Cs
        1.96,   // Ba
        1.80,   // La
        1.63,   // Ce
        1.76,   // Pr
        1.74,   // Nd
        1.73,   // Pm
        1.72,   // Sm
        1.68,   // Eu
        1.69,   // Gd
        1.68,   // Tb
        1.67,   // Dy
        1.66,   // Ho
        1.65,   // Er
        1.64,   // Tm
        1.70,   // Yb
        1.62,   // Lu
        1.52,   // Hf
        1.46,   // Ta
        1.37,   // W
        1.31,   // Re
        1.29,   // Os
        1.22,   // Ir
        1.23,   // Pt
        1.24,   // Au
        1.33,   // Hg
        1.44,   // Tl
        1.44,   // Pb
        1.51,   // Bi
        1.45    // Po
    };

}

// Make available for GFN1 and other xTB methods
double getCovalentRadius(int Z) {
    if (Z < 1 || Z > 86) return 1.5;  // Default for unsupported elements
    return COVALENT_RADII[Z];
}

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
    // Load complete GFN2 parameter database (26 elements, 48 pairs)
    // Claude Generated: Parameter database integration (November 2025)
    if (!m_param_db.loadCompleteGFN2()) {
        CurcumaLogger::error("Failed to load GFN2 parameter database - using legacy parameters");
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing native GFN2-xTB method");
        CurcumaLogger::param("method", "GFN2-xTB");
        CurcumaLogger::param("reference", "Bannwarth et al. JCTC 2019, 15, 1652");
        CurcumaLogger::param("parameter_database", fmt::format("{} elements, {} pairs",
                           m_param_db.getNumElements(), m_param_db.getNumPairs()));
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
    // Load complete GFN2 parameter database
    // Claude Generated: Parameter database integration (November 2025)
    if (!m_param_db.loadCompleteGFN2()) {
        CurcumaLogger::error("Failed to load GFN2 parameter database - using legacy parameters");
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing GFN2-xTB with custom parameters");
        CurcumaLogger::param("parameter_database", fmt::format("{} elements, {} pairs",
                           m_param_db.getNumElements(), m_param_db.getNumPairs()));
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
    // Claude Generated: GFN2 Hamiltonian matrix construction
    // Implements GFN2 Eq. 12: H_ij = E_i δ_ij + scale_ij * S_ij
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652

    Matrix H = Matrix::Zero(basisset.size(), basisset.size());

    for (size_t i = 0; i < basisset.size(); ++i) {
        int atom_i = basisset[i].atom;
        int Z_i = m_atoms[atom_i];
        double CN_i = m_coordination_numbers(atom_i);

        // Determine shell from orbital type (simplified)
        int shell_i = 0;  // s-orbital
        if (basisset[i].type == STO::PX || basisset[i].type == STO::PY || basisset[i].type == STO::PZ) {
            shell_i = 1;  // p-orbital
        }

        for (size_t j = 0; j <= i; ++j) {
            if (i == j) {
                // Diagonal: self-energy with CN dependence
                H(i, i) = getSelfEnergy(Z_i, shell_i, CN_i);
            } else {
                // Off-diagonal: hopping integral (GFN2 Eq. 12b)
                int atom_j = basisset[j].atom;

                // Calculate interatomic distance
                double dx = basisset[i].x - basisset[j].x;
                double dy = basisset[i].y - basisset[j].y;
                double dz = basisset[i].z - basisset[j].z;
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz) * au;  // Bohr to Å

                // Get scaling factor
                double scale = getHamiltonianScale(basisset[i], basisset[j], distance);

                // H_ij = scale * S_ij
                H(i, j) = H(j, i) = scale * S(i, j);
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("Hamiltonian matrix constructed ({} × {})",
                                       basisset.size(), basisset.size()));
        CurcumaLogger::param("H_diagonal_range",
                           fmt::format("[{:.4f}, {:.4f}] Eh",
                                     H.diagonal().minCoeff(),
                                     H.diagonal().maxCoeff()));
    }

    return H;
}

double GFN2::getSelfEnergy(int element, int shell, double CN) const
{
    // Claude Generated: GFN2 self-energy calculation with shell-resolved parameters
    // Formula (GFN2 Eq. 12a): E_ii = E_base + k_CN * CN
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652
    // Updated November 2025: Uses real TBLite-derived parameters for 26 elements

    double E = 0.0;
    bool use_real_params = false;

    // Try to use real shell-resolved parameters from database
    if (m_param_db.hasElement(element)) {
        const auto& elem_params = m_param_db.getElement(element);

        // Check if this shell exists for this element
        auto shell_it = elem_params.shells.find(shell);
        if (shell_it != elem_params.shells.end()) {
            const auto& shell_params = shell_it->second;

            // GFN2 Eq. 12a with real parameters
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
        // Use chemical hardness as base energy (already in Hartree)
        double E_base = -m_params.getHardness(element);

        // CN shift: use shell_hardness (kpoly) parameter scaled
        double k_CN = m_params.getShellHardness(element) * 0.01;  // Scaled to reasonable magnitude

        // Apply CN shift
        E = E_base + k_CN * CN;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("SelfEnergy[Z={},shell={}]", element, shell),
                               fmt::format("E={:.6f} Eh (CN={:.2f}) [legacy fallback]", E, CN));
        }
    }

    return E;
}

double GFN2::getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const
{
    // Claude Generated: GFN2 Hamiltonian scaling factor with pair-specific parameters
    // Formula (GFN2 Eq. 12b): H_ij = scale * S_ij
    // where scale = z_ij * k_pair * k_shell * (1 + 0.02 * ΔEN²) * poly(r)
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652
    // Updated November 2025: Uses real TBLite pair parameters for 48 element pairs

    int atom_i = fi.atom;
    int atom_j = fj.atom;
    int Z_i = m_atoms[atom_i];
    int Z_j = m_atoms[atom_j];

    // Electronegativity difference (Pauling scale)
    double EN_i = m_params.getElectronegativity(Z_i);
    double EN_j = m_params.getElectronegativity(Z_j);
    double delta_EN = EN_i - EN_j;

    // z_ij: orbital overlap scaling factor
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

    // Electronegativity correction (GFN2 Eq. 12b)
    double en_factor = 1.0 + 0.02 * delta_EN * delta_EN;

    // Distance polynomial correction (simplified exponential decay)
    // TODO: Full polynomial from TBLite (requires radial basis parameters)
    double r_bohr = distance / au;  // Convert Å to Bohr
    double poly_r = std::exp(-0.5 * (r_bohr - 3.0));  // Simplified decay

    // Combined scaling factor
    double scale = z_ij * k_pair * k_shell * en_factor * poly_r;

    return scale;
}

// =================================================================================
// Coordination Numbers
// =================================================================================

Vector GFN2::calculateCoordinationNumbers()
{
    // Claude Generated: GFN2 coordination number calculation
    // Formula (GFN2 Eq. 4): CN_i = ∑_{j≠i} [1 / (1 + exp(-k₁(R_cov,ij/R_ij - 1)))]^(k₂)
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652

    Vector CN = Vector::Zero(m_atomcount);

    // GFN2 parameters (from paper)
    const double k1 = 16.0;      // Steepness
    const double k2 = 4.0 / 3.0; // Range decay exponent

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        double R_cov_A = getCovalentRadius(Z_A);

        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;

            int Z_B = m_atoms[B];
            double R_cov_B = getCovalentRadius(Z_B);
            double R_cov = R_cov_A + R_cov_B;

            // Distance in Ångström
            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Avoid division by zero for overlapping atoms
            if (R_AB < 1.0e-6) continue;

            // Counting function: 1 / (1 + exp(-k1 * (R_cov/R - 1)))
            double count = 1.0 / (1.0 + std::exp(-k1 * (R_cov / R_AB - 1.0)));

            // Power k2
            CN(A) += std::pow(count, k2);
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Coordination numbers calculated:");
        for (int i = 0; i < m_atomcount; ++i) {
            CurcumaLogger::param(fmt::format("CN[atom_{}]", i+1),
                               fmt::format("{:.3f}", CN(i)));
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
    // Claude Generated: GFN2 repulsion energy calculation with TBLite parameters
    // Formula: E_rep = ∑_{A<B} V_rep(R_AB)
    // where V_rep is an exponential repulsion potential
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652
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

            // Interatomic distance in Bohr
            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz) / au;  // Å to Bohr

            // Average repulsion exponent
            double alpha_avg = (alpha_A + alpha_B) / 2.0;

            // Exponential repulsion (GFN2 formula)
            double V_rep = (Z_eff_A + Z_eff_B) * std::exp(-alpha_avg * R_AB) / R_AB;

            E_rep += V_rep;
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("E_repulsion", fmt::format("{:.6f} Eh", E_rep));
    }

    return E_rep;
}

double GFN2::calculateCoulombEnergy() const
{
    // Claude Generated: GFN2 Coulomb energy calculation
    // Three components (GFN2 Eq. 7-11):
    // 1. ES2: Effective Coulomb interaction (Eq. 7)
    // 2. ES3: Third-order onsite correction (Eq. 9)
    // 3. AES2: Anisotropic multipole (Eq. 10-11) - TODO for now
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652
    //
    // TODO: Extract real GFN2 gamma parameters and implement full AES2

    // First, calculate Mulliken charges from density matrix
    // q_A = Z_A - ∑_{μ∈A} (P*S)_μμ
    Vector charges = Vector::Zero(m_atomcount);
    Matrix PS = m_density * m_overlap;

    // Map basis functions to atoms and accumulate populations
    for (int mu = 0; mu < m_nbasis; ++mu) {
        int atom = m_basis[mu].atom;
        charges(atom) += PS(mu, mu);
    }

    // Convert populations to charges: q = Z - electrons
    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        charges(A) = Z_A - charges(A);
    }

    // Store charges for property access
    const_cast<GFN2*>(this)->m_charges = charges;

    // ES2: Effective Coulomb energy (GFN2 Eq. 7)
    // E_ES2 = (1/2) * ∑_{A,B} q_A * q_B * γ_AB
    // Updated November 2025: Uses real gamma_ss parameters from TBLite
    double E_ES2 = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];

        // Try to get real gamma parameter (use gamma_ss for onsite)
        double gamma_AA;
        if (m_param_db.hasElement(Z_A)) {
            const auto& elem_A = m_param_db.getElement(Z_A);
            gamma_AA = elem_A.gamma_ss;  // Onsite Coulomb integral
        } else {
            gamma_AA = m_params.getHardness(Z_A);  // Fallback
        }

        // Onsite (A=B)
        E_ES2 += 0.5 * charges(A) * charges(A) * gamma_AA;

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];

            // Get gamma for atom B
            double gamma_BB;
            if (m_param_db.hasElement(Z_B)) {
                const auto& elem_B = m_param_db.getElement(Z_B);
                gamma_BB = elem_B.gamma_ss;
            } else {
                gamma_BB = m_params.getHardness(Z_B);
            }

            // Interatomic distance in Bohr
            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz) / au;  // Å to Bohr

            // Off-site Coulomb kernel (GFN2 Eq. 7)
            // γ_AB = 1 / √(R² + (η_A + η_B)²)
            // where η is related to gamma (chemical hardness)
            double gamma_AB = 1.0 / std::sqrt(R_AB*R_AB + 0.5 * (gamma_AA + gamma_BB));

            E_ES2 += charges(A) * charges(B) * gamma_AB;
        }
    }

    // ES3: Third-order onsite correction (GFN2 Eq. 9)
    // E_ES3 = (1/6) * ∑_A dU_A * q_A³
    double E_ES3 = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        double dU_A = m_params.getHubbardDerivative(Z_A);  // Hubbard derivative

        E_ES3 += (1.0/6.0) * dU_A * std::pow(charges(A), 3.0);
    }

    // AES2: Anisotropic multipole contributions (GFN2 Eq. 10-11)
    // Simplified dipole-dipole interactions with Tang-Toennies damping
    // TODO: Add full quadrupole terms from TBLite
    double E_AES2 = 0.0;

    // Approximate atomic dipoles from charge distribution
    // μ_A ≈ q_A * R_A (simplified, should use density matrix)
    for (int A = 0; A < m_atomcount - 1; ++A) {
        int Z_A = m_atoms[A];
        double R_A = getCovalentRadius(Z_A);

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];
            double R_B = getCovalentRadius(Z_B);

            // Interatomic distance
            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB_ang = std::sqrt(dx*dx + dy*dy + dz*dz);
            double R_AB = R_AB_ang / au;  // Bohr

            // Approximate dipole magnitude
            double mu_A = std::abs(charges(A)) * R_A / au;
            double mu_B = std::abs(charges(B)) * R_B / au;

            // Tang-Toennies damping function
            double alpha_avg = (m_params.getMultipoleRadius(Z_A) +
                              m_params.getMultipoleRadius(Z_B)) / 2.0;
            double x = alpha_avg * R_AB;
            double damp = 1.0 - std::exp(-x) * (1.0 + x + x*x/2.0 + x*x*x/6.0);

            // Dipole-dipole interaction (simplified isotropic approximation)
            double E_dd = -mu_A * mu_B * damp / (R_AB * R_AB * R_AB);
            E_AES2 += E_dd * 0.1;  // Scale factor (approximate)
        }
    }

    double E_coulomb = E_ES2 + E_ES3 + E_AES2;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("E_ES2", fmt::format("{:.6f} Eh", E_ES2));
        CurcumaLogger::param("E_ES3", fmt::format("{:.6f} Eh", E_ES3));
        CurcumaLogger::param("E_coulomb_total", fmt::format("{:.6f} Eh", E_coulomb));
    }

    return E_coulomb;
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
    // Claude Generated: GFN2 analytical gradient calculation
    // Uses Hellmann-Feynman theorem: dE/dR = Tr(P * dH/dR) + dE_rep/dR + dE_coul/dR
    // Reference: C. Bannwarth et al., JCTC 2019, 15, 1652
    //
    // Simplified implementation - numerical gradients used for now
    // TODO: Full analytical derivatives of overlap, Hamiltonian, and energy components

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Calculating GFN2 gradients (numerical)");
    }

    Matrix gradient = Matrix::Zero(m_atomcount, 3);
    const double delta = 1.0e-5;  // Numerical displacement in Ångström

    // Store original geometry
    Matrix geom_orig = m_geometry;
    double E0 = m_total_energy;

    // Numerical gradients: dE/dR ≈ (E(R+δ) - E(R-δ)) / (2δ)
    for (int atom = 0; atom < m_atomcount; ++atom) {
        for (int coord = 0; coord < 3; ++coord) {
            // Forward displacement
            const_cast<GFN2*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) + delta;
            const_cast<GFN2*>(this)->InitialiseMolecule();
            double E_plus = const_cast<GFN2*>(this)->Calculation(false);

            // Backward displacement
            const_cast<GFN2*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) - delta;
            const_cast<GFN2*>(this)->InitialiseMolecule();
            double E_minus = const_cast<GFN2*>(this)->Calculation(false);

            // Central difference
            gradient(atom, coord) = (E_plus - E_minus) / (2.0 * delta);

            // Restore original coordinate
            const_cast<GFN2*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord);
        }
    }

    // Restore original state
    const_cast<GFN2*>(this)->m_geometry = geom_orig;
    const_cast<GFN2*>(this)->InitialiseMolecule();
    const_cast<GFN2*>(this)->m_total_energy = E0;

    if (CurcumaLogger::get_verbosity() >= 3) {
        double grad_norm = gradient.norm();
        CurcumaLogger::param("gradient_norm", fmt::format("{:.6e} Eh/Å", grad_norm));
        CurcumaLogger::param("max_gradient_component",
                           fmt::format("{:.6e} Eh/Å", gradient.cwiseAbs().maxCoeff()));
    }

    // Convert Eh/Å to Eh/Bohr for consistency
    return gradient / au;
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
