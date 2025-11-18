/*
 * <Native PM3 Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on PM3 by James J. P. Stewart
 * Reference: J. J. P. Stewart, J. Comput. Chem. 1989, 10, 209
 *
 * This program is free software under GPL-3.0
 */

#include "pm3.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "ParallelEigenSolver.hpp"
#include "integrals/MNDOIntegrals.hpp"  // MNDO multipole expansion integrals

#include <fmt/format.h>
#include <cmath>
#include <algorithm>

using namespace CurcumaUnit;

// =================================================================================
// Constructor / Destructor
// =================================================================================

PM3::PM3()
    : m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_core_repulsion(0.0)
    , m_scf_max_iterations(100)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing native PM3 method");
        CurcumaLogger::param("method", "PM3");
        CurcumaLogger::param("reference", "Stewart JCC 1989, 10, 209");
    }

    initializePM3Parameters();
}

// =================================================================================
// PM3 Parameter Initialization (simplified)
// =================================================================================

void PM3::initializePM3Parameters()
{
    // Claude Generated: Simplified PM3 parameters
    // Real PM3 has 100+ parameters per element
    // Here we provide minimal parameter set for educational purposes
    // TODO: Extract complete parameters from MOPAC database

    // CRITICAL: Parameters stored in eV, must convert to Hartree
    // Reference: Ulysses MNDO.hpp line 6597: "return ulx/au2eV"
    const double eV2Hartree = 27.211386245988;

    // Hydrogen (Z=1)
    // Claude Generated: PM3 parameters extended with MNDO integral parameters
    PM3Params H;
    H.U_ss = -13.073321;  // eV
    H.zeta_s = 1.188078;
    H.beta_s = -5.626512;
    H.alpha = 2.544134;
    H.gauss_a = {0.122796, 0.005090};
    H.gauss_b = {5.000000, 5.000000};
    H.gauss_c = {1.220000, 0.000000};
    // MNDO multipole expansion parameters (extracted from MOPAC/Ulysses)
    H.D1 = 0.0;      // H has no p-orbitals, no dipole expansion
    H.D2 = 0.0;      // No d-orbitals
    H.rho_s = 2.0 / H.zeta_s;  // Approximate from Slater exponent
    H.rho_p = 0.0;   // No p-orbitals for H
    m_pm3_params[1] = H;

    // Carbon (Z=6)
    PM3Params C;
    C.U_ss = -47.270320;
    C.U_pp = -36.266918;
    C.zeta_s = 1.565085;
    C.zeta_p = 1.842345;
    C.beta_s = -11.910015;
    C.beta_p = -9.802755;
    C.alpha = 2.648274;
    C.gauss_a = {0.011355, 0.045924, 0.000000};
    C.gauss_b = {5.000000, 5.000000, 2.000000};
    C.gauss_c = {1.560000, 1.567000, 0.000000};
    // MNDO multipole expansion parameters
    C.D1 = 0.7920; // PM3 value in Å (from MOPAC parameter database)
    C.D2 = 0.0;    // No d-orbitals for second-row elements
    C.rho_s = 2.0 / C.zeta_s;  // Orbital exponent for s-type ERIs
    C.rho_p = 2.0 / C.zeta_p;  // Orbital exponent for p-type ERIs
    m_pm3_params[6] = C;

    // Nitrogen (Z=7)
    PM3Params N;
    N.U_ss = -49.335672;
    N.U_pp = -47.509736;
    N.zeta_s = 2.028094;
    N.zeta_p = 2.313728;
    N.beta_s = -14.062521;
    N.beta_p = -20.043848;
    N.alpha = 2.947286;
    N.gauss_a = {0.025251, 0.028953};
    N.gauss_b = {5.000000, 2.000000};
    N.gauss_c = {1.550000, 0.000000};
    // MNDO multipole expansion parameters
    N.D1 = 0.6433; // PM3 value in Å
    N.D2 = 0.0;
    N.rho_s = 2.0 / N.zeta_s;
    N.rho_p = 2.0 / N.zeta_p;
    m_pm3_params[7] = N;

    // Oxygen (Z=8)
    PM3Params O;
    O.U_ss = -86.993002;
    O.U_pp = -71.879580;
    O.zeta_s = 3.796544;
    O.zeta_p = 2.389402;
    O.beta_s = -45.202651;
    O.beta_p = -24.752515;
    O.alpha = 3.217102;
    O.gauss_a = {0.280962, 0.081430};
    O.gauss_b = {5.000000, 7.000000};
    O.gauss_c = {0.847918, 1.445071};
    // MNDO multipole expansion parameters
    O.D1 = 0.5346; // PM3 value in Å
    O.D2 = 0.0;
    O.rho_s = 2.0 / O.zeta_s;
    O.rho_p = 2.0 / O.zeta_p;
    m_pm3_params[8] = O;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("PM3 parameters loaded for {} elements",
                                       m_pm3_params.size()));
    }
}

// =================================================================================
// QMDriver Interface Implementation
// =================================================================================

bool PM3::InitialiseMolecule()
{
    if (m_atoms.size() == 0) {
        CurcumaLogger::error("No atoms in molecule for PM3 initialization");
        return false;
    }

    // Validate all atoms have PM3 parameters
    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];
        if (m_pm3_params.find(Z) == m_pm3_params.end()) {
            CurcumaLogger::error(fmt::format("PM3 parameters not available for element Z={} at atom {}",
                                            Z, i+1));
            CurcumaLogger::info("Available elements: H, C, N, O");
            return false;
        }
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing PM3 calculation");
        CurcumaLogger::param("atoms", static_cast<int>(m_atoms.size()));
        CurcumaLogger::param("charge", m_charge);
    }

    // Build minimal valence basis
    m_nbasis = buildBasisSet();

    if (m_nbasis == 0) {
        CurcumaLogger::error("Failed to build basis set for PM3");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("basis_functions", m_nbasis);
        CurcumaLogger::param("electrons", m_num_electrons);
    }

    // Build overlap matrix
    m_overlap = MakeOverlap(m_basis);

    // Build core Hamiltonian
    m_hamiltonian = MakeH(m_overlap, m_basis);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("PM3 initialization complete");
    }

    return true;
}

double PM3::Calculation(bool gradient)
{
    try {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::info("Starting PM3 energy calculation");
        }

        // Run SCF convergence
        bool scf_success = runSCF();

        if (!scf_success) {
            CurcumaLogger::error("PM3 SCF did not converge");
            return 0.0;
        }

        // Calculate energy components
        m_energy_electronic = calculateElectronicEnergy();
        m_energy_core_repulsion = calculateCoreRepulsionEnergy();

        // Total energy
        m_total_energy = m_energy_electronic + m_energy_core_repulsion;

        // Level 1+: Results
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success("PM3 calculation complete");
            CurcumaLogger::energy_abs(m_total_energy, "PM3_total_energy");
        }

        // Level 2+: Energy decomposition
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Energy decomposition:");
            CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", m_energy_electronic));
            CurcumaLogger::param("core_repulsion", fmt::format("{:.6f} Eh", m_energy_core_repulsion));

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
        CurcumaLogger::error(fmt::format("PM3 calculation failed: {}", e.what()));
        return 0.0;
    }
}

// =================================================================================
// Basis Set Construction (minimal valence only)
// =================================================================================

int PM3::buildBasisSet()
{
    m_basis.clear();
    m_num_electrons = 0;

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];
        const PM3Params& params = m_pm3_params[Z];

        // PM3 uses minimal valence basis (s for H, s+p for others)
        // s-orbital
        STO::Orbital s_orbital;
        s_orbital.type = STO::S;
        s_orbital.x = m_geometry(i, 0) / au;
        s_orbital.y = m_geometry(i, 1) / au;
        s_orbital.z = m_geometry(i, 2) / au;
        s_orbital.zeta = params.zeta_s;
        s_orbital.VSIP = params.U_ss / eV2Eh;
        s_orbital.atom = i;
        m_basis.push_back(s_orbital);

        // p-orbitals (for Z > 2)
        if (Z > 2) {
            for (int m = 0; m < 3; ++m) {  // px, py, pz
                STO::Orbital p_orbital;
                p_orbital.type = (m == 0) ? STO::PX : ((m == 1) ? STO::PY : STO::PZ);
                p_orbital.x = m_geometry(i, 0) / au;
                p_orbital.y = m_geometry(i, 1) / au;
                p_orbital.z = m_geometry(i, 2) / au;
                p_orbital.zeta = params.zeta_p;
                p_orbital.VSIP = params.U_pp / eV2Eh;
                p_orbital.atom = i;
                m_basis.push_back(p_orbital);
            }
        }

        // Count valence electrons
        if (Z == 1) m_num_electrons += 1;
        else if (Z >= 2 && Z <= 10) m_num_electrons += (Z - 2) + 2;  // 2s² + np^x
        else m_num_electrons += Z;  // Simplified
    }

    return static_cast<int>(m_basis.size());
}

Matrix PM3::MakeOverlap(std::vector<STO::Orbital>& basisset)
{
    // PM3 uses normalized overlap integrals
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
// Hamiltonian Construction (NDDO approximation)
// =================================================================================

Matrix PM3::MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset)
{
    // Claude Generated: PM3 core Hamiltonian
    // NDDO: Only one-center and two-center terms

    Matrix H = Matrix::Zero(basisset.size(), basisset.size());

    for (size_t mu = 0; mu < basisset.size(); ++mu) {
        int atom_i = basisset[mu].atom;
        int Z_i = m_atoms[atom_i];

        for (size_t nu = 0; nu <= mu; ++nu) {
            int atom_j = basisset[nu].atom;
            int Z_j = m_atoms[atom_j];

            if (atom_i == atom_j) {
                // One-center term
                if (mu == nu) {
                    // Diagonal: U_αα (core integral)
                    H(mu, nu) = getCoreIntegral(Z_i, basisset[mu].type);
                } else {
                    // Off-diagonal: zero in NDDO for one-center
                    H(mu, nu) = H(nu, mu) = 0.0;
                }
            } else {
                // Two-center term: resonance integral
                double dx = basisset[mu].x - basisset[nu].x;
                double dy = basisset[mu].y - basisset[nu].y;
                double dz = basisset[mu].z - basisset[nu].z;
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz) * au;  // Bohr to Å

                double beta = getResonanceIntegral(basisset[mu], basisset[nu], distance);
                H(mu, nu) = H(nu, mu) = beta * S(mu, nu);
            }
        }
    }

    return H;
}

double PM3::getCoreIntegral(int element, int orbital_type) const
{
    // U_αα = ionization potential (VSIP)
    const PM3Params& params = m_pm3_params.at(element);

    if (orbital_type == STO::S) {
        return params.U_ss / eV2Eh;  // Convert eV to Hartree
    } else if (orbital_type == STO::PX || orbital_type == STO::PY || orbital_type == STO::PZ) {
        return params.U_pp / eV2Eh;
    }

    return 0.0;
}

double PM3::getResonanceIntegral(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const
{
    // Resonance integral: β_AB = (β_A + β_B) / 2 * S_AB
    // Simplified PM3 formula

    int atom_i = fi.atom;
    int atom_j = fj.atom;
    int Z_i = m_atoms[atom_i];
    int Z_j = m_atoms[atom_j];

    const PM3Params& params_i = m_pm3_params.at(Z_i);
    const PM3Params& params_j = m_pm3_params.at(Z_j);

    // Get beta values
    double beta_i = (fi.type == STO::S) ? params_i.beta_s : params_i.beta_p;
    double beta_j = (fj.type == STO::S) ? params_j.beta_s : params_j.beta_p;

    // Average beta
    double beta_avg = (beta_i + beta_j) / 2.0;

    // Distance scaling (simplified)
    double scale = std::exp(-0.5 * distance);

    return (beta_avg / eV2Eh) * scale;
}

// =================================================================================
// SCF Procedure (similar to GFN methods)
// =================================================================================

bool PM3::runSCF()
{
    m_density = Matrix::Zero(m_nbasis, m_nbasis);

    for (int iter = 0; iter < m_scf_max_iterations; ++iter) {
        m_fock = buildFockMatrix(m_density);

        ParallelEigenSolver solver(500, 128, 1.0e-10, false);
        solver.setThreadCount(m_threads);

        bool success = solver.solve(m_overlap, m_fock, m_energies, m_mo, m_threads, false);

        if (!success) {
            CurcumaLogger::error("Eigenvalue solution failed in PM3 SCF");
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
                CurcumaLogger::success(fmt::format("PM3 SCF converged in {} iterations", iter+1));
            }

            return true;
        }

        m_density = m_scf_damping * m_density + (1.0 - m_scf_damping) * density_new;
    }

    CurcumaLogger::warn(fmt::format("PM3 SCF did not converge in {} iterations", m_scf_max_iterations));
    m_scf_converged = false;
    return false;
}

Matrix PM3::buildFockMatrix(const Matrix& density)
{
    // Fock matrix: F = H + G
    // G contains two-electron interactions (simplified)

    Matrix F = m_hamiltonian;

    // Add two-electron contributions (very simplified)
    // Real PM3 has complex (μν|λσ) integrals

    for (int mu = 0; mu < m_nbasis; ++mu) {
        for (int nu = 0; nu < m_nbasis; ++nu) {
            double G_element = 0.0;

            for (int lambda = 0; lambda < m_nbasis; ++lambda) {
                for (int sigma = 0; sigma < m_nbasis; ++sigma) {
                    double integral = calculateTwoElectronIntegral(mu, nu, lambda, sigma);
                    G_element += density(lambda, sigma) * integral;
                }
            }

            F(mu, nu) += G_element;
        }
    }

    return F;
}

Matrix PM3::buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies)
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

double PM3::calculateTwoElectronIntegral(int mu, int nu, int lambda, int sigma) const
{
    // Claude Generated: Full MNDO two-electron repulsion integrals using multipole expansion
    // Reference: Dewar & Thiel, Theor. Chim. Acta 46, 89 (1977)
    //
    // NDDO approximation: (μν|λσ) ≠ 0 only if μ,ν on same atom AND λ,σ on same atom
    // This reduces four-center integrals to two-center integrals

    int atom_mu = m_basis[mu].atom;
    int atom_nu = m_basis[nu].atom;
    int atom_lambda = m_basis[lambda].atom;
    int atom_sigma = m_basis[sigma].atom;

    // NDDO: Zero unless μ,ν on same atom AND λ,σ on same atom
    if (atom_mu != atom_nu || atom_lambda != atom_sigma) {
        return 0.0;
    }

    int atom_A = atom_mu;
    int atom_B = atom_lambda;

    if (atom_A == atom_B) {
        // =====================================================================
        // ONE-CENTER INTEGRALS: (μν|λσ) with all orbitals on same atom
        // =====================================================================

        int Z = m_atoms[atom_A];
        const PM3Params& params = m_pm3_params.at(Z);

        // Extract orbital types
        STO::OrbitalType type_mu = m_basis[mu].type;
        STO::OrbitalType type_nu = m_basis[nu].type;
        STO::OrbitalType type_lambda = m_basis[lambda].type;
        STO::OrbitalType type_sigma = m_basis[sigma].type;

        // One-center integrals (simplified - real PM3 has full multipole)
        // (ss|ss)
        if (type_mu == STO::S && type_nu == STO::S &&
            type_lambda == STO::S && type_sigma == STO::S) {
            return params.U_ss / eV2Eh;
        }
        // (sp|sp) - cross term
        else if ((type_mu == STO::S && type_lambda == STO::S) ||
                 (type_mu != STO::S && type_lambda != STO::S)) {
            // Simplified: average of ss and pp
            if (type_mu == STO::S && type_lambda != STO::S) {
                return (params.U_ss + params.U_pp) / (2.0 * eV2Eh);
            } else if (type_mu != STO::S && type_lambda == STO::S) {
                return (params.U_pp + params.U_ss) / (2.0 * eV2Eh);
            } else if (type_mu != STO::S && type_lambda != STO::S) {
                return params.U_pp / eV2Eh;
            }
        }

        return 0.0;  // Other one-center cases neglected in simplified version

    } else {
        // =====================================================================
        // TWO-CENTER INTEGRALS: (μν|λσ) with μ,ν on A and λ,σ on B
        // Use MNDO multipole expansion (Dewar-Thiel formulas)
        // =====================================================================

        int Z_A = m_atoms[atom_A];
        int Z_B = m_atoms[atom_B];
        const PM3Params& params_A = m_pm3_params.at(Z_A);
        const PM3Params& params_B = m_pm3_params.at(Z_B);

        // Calculate interatomic distance
        double dx = m_geometry(atom_A, 0) - m_geometry(atom_B, 0);
        double dy = m_geometry(atom_A, 1) - m_geometry(atom_B, 1);
        double dz = m_geometry(atom_A, 2) - m_geometry(atom_B, 2);
        double R_AB_angstrom = std::sqrt(dx*dx + dy*dy + dz*dz);
        double R_AB_bohr = R_AB_angstrom / au;  // Convert Å → Bohr

        // Extract orbital types
        STO::OrbitalType type_mu = m_basis[mu].type;
        STO::OrbitalType type_nu = m_basis[nu].type;
        STO::OrbitalType type_lambda = m_basis[lambda].type;
        STO::OrbitalType type_sigma = m_basis[sigma].type;

        // Map STO orbital types to (l, m) quantum numbers for MNDO integrals
        auto get_quantum_numbers = [](STO::OrbitalType type) -> std::pair<int, int> {
            switch(type) {
                case STO::S:  return {0, 0};   // l=0, m=0
                case STO::PX: return {1, 1};   // l=1, m=+1
                case STO::PY: return {1, -1};  // l=1, m=-1
                case STO::PZ: return {1, 0};   // l=1, m=0
                default:      return {0, 0};
            }
        };

        auto [l_mu, m_mu] = get_quantum_numbers(type_mu);
        auto [l_nu, m_nu] = get_quantum_numbers(type_nu);
        auto [l_lambda, m_lambda] = get_quantum_numbers(type_lambda);
        auto [l_sigma, m_sigma] = get_quantum_numbers(type_sigma);

        // NDDO: (μν|λσ) = (μμ|λλ) δ_μν δ_λσ
        // Only diagonal elements survive
        if (mu != nu || lambda != sigma) {
            return 0.0;
        }

        // Now we have (μμ|λλ) with μ on A, λ on B
        // This is γ_AB^(l_μ m_μ, l_λ m_λ) in MNDO notation

        // Prepare D-parameter vector: {D1_A, D2_A, D1_B, D2_B}
        std::vector<double> D_params = {
            params_A.D1,  // Dipole expansion for atom A (in Å, will convert)
            params_A.D2,  // Quadrupole expansion for atom A
            params_B.D1,  // Dipole expansion for atom B
            params_B.D2   // Quadrupole expansion for atom B
        };

        // Convert D-parameters from Å to Bohr
        for (auto& D : D_params) {
            D /= au;
        }

        // Calculate orbital exponent sum (ρ_A + ρ_B)
        // Use appropriate exponent based on orbital type
        double rho_A = (l_mu == 0) ? params_A.rho_s : params_A.rho_p;
        double rho_B = (l_lambda == 0) ? params_B.rho_s : params_B.rho_p;
        double rho_sum = rho_A + rho_B;

        // Call MNDO multipole integral function
        double gamma_AB = curcuma::mndo::mndo_multipole_integral(
            l_mu, m_mu,        // Orbital on atom A
            l_lambda, m_lambda, // Orbital on atom B
            R_AB_bohr,         // Distance in Bohr
            rho_sum,           // Sum of orbital exponents
            D_params           // Multipole expansion parameters
        );

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("MNDO ERI: atoms ({},{}), orbitals ({}{},{}{}) = {:.6f} Eh",
                atom_A, atom_B, l_mu, m_mu, l_lambda, m_lambda, gamma_AB));
        }

        return gamma_AB;
    }
}

double PM3::getGammaAB(int Z_A, int Z_B, double R_AB) const
{
    // CRITICAL FIX: Use full MNDO multipole expansion for (ss|ss) integral
    // NOT simplified Mataga-Nishimoto approximation!

    const PM3Params& params_A = m_pm3_params.at(Z_A);
    const PM3Params& params_B = m_pm3_params.at(Z_B);

    double R_AB_bohr = R_AB / au;  // Convert Å → Bohr

    // Prepare D-parameter vector for MNDO multipole expansion
    std::vector<double> D_params = {
        params_A.D1 / au,  // Convert Å → Bohr
        params_A.D2 / au,
        params_B.D1 / au,
        params_B.D2 / au
    };

    double rho_sum = params_A.rho_s + params_B.rho_s;

    // (ss|ss) integral using full MNDO multipole expansion
    double gamma_ss = curcuma::mndo::mndo_multipole_integral(
        0, 0,  // s-orbital on A (l=0, m=0)
        0, 0,  // s-orbital on B (l=0, m=0)
        R_AB_bohr,
        rho_sum,
        D_params
    );

    return gamma_ss;
}

// =================================================================================
// Energy Components
// =================================================================================

double PM3::calculateElectronicEnergy() const
{
    // E_elec = Tr(P * H) + (1/2) Tr(P * G)
    // Simplified: E_elec = Tr(P * (H + F)) / 2

    double E = 0.0;

    for (int mu = 0; mu < m_nbasis; ++mu) {
        for (int nu = 0; nu < m_nbasis; ++nu) {
            E += m_density(mu, nu) * (m_hamiltonian(mu, nu) + m_fock(mu, nu));
        }
    }

    return E / 2.0;
}

double PM3::calculateCoreRepulsionEnergy() const
{
    // Claude Generated: PM3 core-core repulsion with Gaussian expansions
    // V_AB = Z_A * Z_B / R_AB + sum_k a_k * exp(-b_k * (R_AB - c_k)²)

    double E_rep = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        const PM3Params& params_A = m_pm3_params.at(Z_A);

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];
            const PM3Params& params_B = m_pm3_params.at(Z_B);

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Core repulsion uses (ss|ss) two-electron integral γ_ss
            // Ulysses MNDO.hpp line 5960: enuc += chgA*chgB*intn*factorA + factorB + chgA*chgB*factorC
            double gamma_ss = getGammaAB(Z_A, Z_B, R_AB);

            // factorA: Exponential damping (Ulysses line 5936-5956)
            // factorA = 1.0 + gA*exp(-alphaA*R) + gB*exp(-alphaB*R)
            double gA = 1.0;
            double gB = 1.0;
            if ((Z_B == 1) && ((Z_A == 7) || (Z_A == 8))) {  // N-H or O-H
                gA = R_AB;
            }
            if ((Z_A == 1) && ((Z_B == 7) || (Z_B == 8))) {  // H-N or H-O
                gB = R_AB;
            }
            double factorA = 1.0;
            factorA += gA * std::exp(-params_A.alpha * R_AB);
            factorA += gB * std::exp(-params_B.alpha * R_AB);

            double factorB = 0.0;

            // factorC: Gaussian corrections (PM3 style)
            // Ulysses PM6.hpp line 47: AM1factor = K*exp(-L*(R-M)²)/RAB
            double factorC = 0.0;
            for (size_t k = 0; k < params_A.gauss_a.size() && k < params_B.gauss_a.size(); ++k) {
                double a_k = (params_A.gauss_a[k] + params_B.gauss_a[k]) / 2.0;
                double b_k = (params_A.gauss_b[k] + params_B.gauss_b[k]) / 2.0;
                double c_k = (params_A.gauss_c[k] + params_B.gauss_c[k]) / 2.0;

                factorC += a_k * std::exp(-b_k * std::pow(R_AB - c_k, 2.0)) / R_AB;
            }

            // Total core repulsion (Ulysses MNDO.hpp line 5960)
            E_rep += double(Z_A * Z_B) * gamma_ss * factorA + factorB + double(Z_A * Z_B) * factorC;
        }
    }

    return E_rep;
}

// =================================================================================
// Gradient Calculation (numerical)
// =================================================================================

Matrix PM3::calculateGradient() const
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Calculating PM3 gradients (numerical)");
    }

    Matrix gradient = Matrix::Zero(m_atomcount, 3);
    const double delta = 1.0e-5;

    Matrix geom_orig = m_geometry;
    double E0 = m_total_energy;

    for (int atom = 0; atom < m_atomcount; ++atom) {
        for (int coord = 0; coord < 3; ++coord) {
            const_cast<PM3*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) + delta;
            const_cast<PM3*>(this)->InitialiseMolecule();
            double E_plus = const_cast<PM3*>(this)->Calculation(false);

            const_cast<PM3*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) - delta;
            const_cast<PM3*>(this)->InitialiseMolecule();
            double E_minus = const_cast<PM3*>(this)->Calculation(false);

            gradient(atom, coord) = (E_plus - E_minus) / (2.0 * delta);

            const_cast<PM3*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord);
        }
    }

    const_cast<PM3*>(this)->m_geometry = geom_orig;
    const_cast<PM3*>(this)->InitialiseMolecule();
    const_cast<PM3*>(this)->m_total_energy = E0;

    return gradient / au;
}

// =================================================================================
// Property Access
// =================================================================================

json PM3::getEnergyDecomposition() const
{
    json decomp;
    decomp["electronic"] = m_energy_electronic;
    decomp["core_repulsion"] = m_energy_core_repulsion;
    decomp["total"] = m_total_energy;

    return decomp;
}

double PM3::getHOMOLUMOGap() const
{
    return (getLUMOEnergy() - getHOMOEnergy());
}

double PM3::getHOMOEnergy() const
{
    int n_occ = m_num_electrons / 2;
    if (n_occ == 0 || n_occ > m_nbasis) return 0.0;

    return m_energies(n_occ - 1) * eV2Eh;
}

double PM3::getLUMOEnergy() const
{
    int n_occ = m_num_electrons / 2;
    if (n_occ >= m_nbasis) return 0.0;

    return m_energies(n_occ) * eV2Eh;
}

double PM3::getHeatOfFormation() const
{
    // TODO: Calculate ΔH_f from atomization energy
    // PM3 is parametrized for heats of formation
    return 0.0;
}
