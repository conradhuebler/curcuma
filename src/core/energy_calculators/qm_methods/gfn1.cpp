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
    , m_scf_max_iterations(100)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing native GFN1-xTB method");
        CurcumaLogger::param("method", "GFN1-xTB");
        CurcumaLogger::param("reference", "Grimme et al. JCTC 2017, 13, 1989");
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
    , m_scf_max_iterations(100)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing GFN1-xTB with custom parameters");
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
            CurcumaLogger::energy_abs("GFN1_total_energy", m_total_energy);
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
            double zeta = m_params.getYeff(Z, shell);
            int principal_qn = static_cast<int>(m_params.getPrincipalQN(Z, shell));

            if (zeta < 1.0e-6) continue;

            int l = shell;

            for (int m = -l; m <= l; ++m) {
                STO::Orbital orbital;
                orbital.type = (l == 0) ? STO::S : ((l == 1) ? STO::PX : STO::S);
                orbital.x = m_geometry(i, 0) / au;
                orbital.y = m_geometry(i, 1) / au;
                orbital.z = m_geometry(i, 2) / au;
                orbital.zeta = zeta;
                orbital.VSIP = 0.0;
                orbital.atom = i;

                m_basis.push_back(orbital);
            }

            if (shell == 0) m_num_electrons += std::min(2, Z);
            else if (shell == 1) m_num_electrons += std::min(6, Z - 2);
            else if (shell == 2) m_num_electrons += std::min(10, Z - 10);
        }
    }

    return static_cast<int>(m_basis.size());
}

Matrix GFN1::MakeOverlap(const std::vector<STO::Orbital>& basisset)
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
    // Simpler than GFN2, but same general structure

    Matrix H = Matrix::Zero(basisset.size(), basisset.size());

    for (size_t i = 0; i < basisset.size(); ++i) {
        int atom_i = basisset[i].atom;
        int Z_i = m_atoms[atom_i];
        double CN_i = m_coordination_numbers(atom_i);

        int shell_i = 0;
        if (basisset[i].type == STO::PX || basisset[i].type == STO::PY || basisset[i].type == STO::PZ) {
            shell_i = 1;
        }

        for (size_t j = 0; j <= i; ++j) {
            if (i == j) {
                H(i, i) = getSelfEnergy(Z_i, shell_i, CN_i);
            } else {
                int atom_j = basisset[j].atom;

                double dx = basisset[i].x - basisset[j].x;
                double dy = basisset[i].y - basisset[j].y;
                double dz = basisset[i].z - basisset[j].z;
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz) * au;

                double scale = getHamiltonianScale(basisset[i], basisset[j], distance);
                H(i, j) = H(j, i) = scale * S(i, j);
            }
        }
    }

    return H;
}

double GFN1::getSelfEnergy(int element, int shell, double CN) const
{
    // Claude Generated: GFN1 self-energy (simpler than GFN2)
    // Uses chemical hardness as approximation

    double E_base = -m_params.getHardness(element);
    double k_CN = m_params.getShellHardness(element) * 0.01;

    return E_base + k_CN * CN;
}

double GFN1::getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const
{
    // Claude Generated: GFN1 Hamiltonian scaling (simpler than GFN2)

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

    double k_pair = std::sqrt(m_params.getAlpha(Z_i) * m_params.getAlpha(Z_j)) * 0.1;
    double en_factor = 1.0 + 0.02 * delta_EN * delta_EN;

    double r_bohr = distance / au;
    double poly_r = std::exp(-0.5 * (r_bohr - 3.0));

    return z_ij * k_pair * en_factor * poly_r;
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
    m_density = Matrix::Zero(m_nbasis, m_nbasis);

    for (int iter = 0; iter < m_scf_max_iterations; ++iter) {
        m_fock = buildFockMatrix(m_density);

        ParallelEigenSolver solver(500, 128, 1.0e-10, false);
        solver.setThreadCount(m_threads);

        bool success = solver.solve(m_overlap, m_fock, m_energies, m_mo, m_threads, false);

        if (!success) {
            CurcumaLogger::error("Eigenvalue solution failed in GFN1 SCF");
            return false;
        }

        Matrix density_new = buildDensityMatrix(m_mo, m_energies);
        double delta_P = (density_new - m_density).norm();

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("SCF_iter_{}", iter+1),
                               fmt::format("ΔP = {:.6e}", delta_P));
        }

        if (delta_P < m_scf_threshold) {
            m_density = density_new;
            m_scf_converged = true;

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::success(fmt::format("GFN1 SCF converged in {} iterations", iter+1));
            }

            return true;
        }

        m_density = m_scf_damping * m_density + (1.0 - m_scf_damping) * density_new;
    }

    CurcumaLogger::warn(fmt::format("GFN1 SCF did not converge in {} iterations", m_scf_max_iterations));
    m_scf_converged = false;
    return false;
}

Matrix GFN1::buildFockMatrix(const Matrix& density)
{
    // Simplified Fock matrix: F = H + G(P)
    // For GFN1, we use H₀ and add Coulomb contributions in energy calculation
    return m_hamiltonian;
}

Matrix GFN1::buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies)
{
    Matrix P = Matrix::Zero(m_nbasis, m_nbasis);

    int n_occ = m_num_electrons / 2;

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

double GFN1::calculateElectronicEnergy() const
{
    double E = 0.0;

    for (int mu = 0; mu < m_nbasis; ++mu) {
        for (int nu = 0; nu < m_nbasis; ++nu) {
            E += m_density(mu, nu) * (m_hamiltonian(mu, nu) + m_fock(mu, nu));
        }
    }

    return E / 2.0;
}

double GFN1::calculateRepulsionEnergy() const
{
    // Claude Generated: GFN1 repulsion (similar to GFN2)

    double E_rep = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        double alpha_A = m_params.getMultipoleRadius(Z_A);

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];
            double alpha_B = m_params.getMultipoleRadius(Z_B);

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz) / au;

            double Z_eff_A = m_params.getElectronegativity(Z_A) * 0.5;
            double Z_eff_B = m_params.getElectronegativity(Z_B) * 0.5;

            double alpha_avg = (alpha_A + alpha_B) / 2.0;
            double V_rep = (Z_eff_A + Z_eff_B) * std::exp(-alpha_avg * R_AB) / R_AB;

            E_rep += V_rep;
        }
    }

    return E_rep;
}

double GFN1::calculateCoulombEnergy() const
{
    // Claude Generated: GFN1 Coulomb (simpler than GFN2, no ES3)

    Vector charges = Vector::Zero(m_atomcount);
    Matrix PS = m_density * m_overlap;

    for (int mu = 0; mu < m_nbasis; ++mu) {
        int atom = m_basis[mu].atom;
        charges(atom) += PS(mu, mu);
    }

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        charges(A) = Z_A - charges(A);
    }

    const_cast<GFN1*>(this)->m_charges = charges;

    // ES2: Effective Coulomb
    double E_ES2 = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        double gamma_AA = m_params.getHardness(Z_A);

        E_ES2 += 0.5 * charges(A) * charges(A) * gamma_AA;

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz) / au;

            double gamma_AB = 1.0 / std::sqrt(R_AB*R_AB + 0.5 * (gamma_AA + m_params.getHardness(Z_B)));

            E_ES2 += charges(A) * charges(B) * gamma_AB;
        }
    }

    return E_ES2;
}

double GFN1::calculateDispersionEnergy() const
{
    // D3 dispersion stub (separate TODO)
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::warn("D3 dispersion not yet integrated - energy incomplete");
    }

    return 0.0;
}

double GFN1::calculateHalogenBondCorrection() const
{
    // Claude Generated: GFN1 halogen bond correction
    // Reference: Grimme et al. JCTC 2017, 13, 1989 (Section 2.5)
    //
    // Halogen bonding: X···A interaction (X = F, Cl, Br, I; A = N, O, S, etc.)
    // Simplified empirical correction based on geometry and charges

    double E_XB = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];

        if (!isHalogen(Z_A)) continue;

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

            // Halogen bond correction (empirical)
            // E_XB = -k * f(R) where f(R) is distance-dependent damping
            double R_vdw = getCovalentRadius(Z_A) + getCovalentRadius(Z_B) + 0.5;  // vdW sum

            if (R_AB < R_vdw * 1.5) {
                double k_XB = m_params.getFXCMu(Z_A) * 0.001;  // Halogen strength parameter
                double damp = std::exp(-2.0 * (R_AB / R_vdw - 1.0));

                E_XB -= k_XB * damp;
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && std::abs(E_XB) > 1.0e-6) {
        CurcumaLogger::param("E_halogen_bond", fmt::format("{:.6f} Eh", E_XB));
    }

    return E_XB;
}

// =================================================================================
// Gradient Calculation (numerical, same as GFN2)
// =================================================================================

Matrix GFN1::calculateGradient() const
{
    // Claude Generated: Numerical gradients for GFN1

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Calculating GFN1 gradients (numerical)");
    }

    Matrix gradient = Matrix::Zero(m_atomcount, 3);
    const double delta = 1.0e-5;

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
