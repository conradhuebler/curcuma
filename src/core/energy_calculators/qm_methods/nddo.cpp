/*
 * <Unified NDDO Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Generic NDDO solver for MNDO, AM1, PM3, PM6.
 *
 * Key physics:
 *   Two-center integrals are computed in the diatomic frame (z along A->B bond)
 *   using Dewar-Thiel multipole expansion, then rotated to laboratory frame
 *   via SP (3x3) and PP (6x6) rotation matrices.
 *
 * References:
 *   MNDO: M. J. S. Dewar, W. Thiel, JACS 1977, 99, 4899
 *   AM1:  M. J. S. Dewar et al., JACS 1985, 107, 3902
 *   PM3:  J. J. P. Stewart, J. Comput. Chem. 1989, 10, 209
 *   PM6:  J. J. P. Stewart, J. Mol. Model. 2007, 13, 1173
 *   Rotation: W. Thiel, A. A. Voityuk, Theor. Chim. Acta 81, 391 (1992)
 *
 * Claude Generated: NDDO implementation with diatomic-frame integral rotation
 *
 * This program is free software under GPL-3.0
 */

#include "nddo.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "ParallelEigenSolver.hpp"
#include "integrals/MNDOIntegrals.hpp"

#include <fmt/format.h>
#include <cmath>
#include <algorithm>

using namespace CurcumaUnit;

// =================================================================================
// Constructor
// =================================================================================

NDDO::NDDO(NDDOMethodType method_type)
    : m_method_type(method_type)
    , m_params(&NDDOParams::getParams(method_type))
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_core_repulsion(0.0)
    , m_scf_max_iterations(100)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.4)
    , m_scf_converged(false)
{
    const char* name = NDDOParams::methodName(m_method_type);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Initializing native {} method", name));
        CurcumaLogger::param("method", name);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("{} parameters loaded for {} elements",
                                        name, m_params->size()));
    }
}

// =================================================================================
// QMDriver Interface Implementation
// =================================================================================

bool NDDO::InitialiseMolecule()
{
    const char* name = NDDOParams::methodName(m_method_type);

    if (m_atoms.size() == 0) {
        CurcumaLogger::error(fmt::format("No atoms in molecule for {} initialization", name));
        return false;
    }

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];
        if (m_params->find(Z) == m_params->end()) {
            CurcumaLogger::error(fmt::format("{} parameters not available for element Z={} at atom {}",
                                            name, Z, i+1));
            return false;
        }
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info(fmt::format("Initializing {} calculation", name));
        CurcumaLogger::param("atoms", static_cast<int>(m_atoms.size()));
        CurcumaLogger::param("charge", m_charge);
    }

    m_nbasis = buildBasisSet();

    if (m_nbasis == 0) {
        CurcumaLogger::error(fmt::format("Failed to build basis set for {}", name));
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("basis_functions", m_nbasis);
        CurcumaLogger::param("electrons", m_num_electrons);
    }

    // Precompute two-center integral blocks (rotated to lab frame)
    // and nuclear attraction integrals VAB
    m_gamma_blocks.resize(m_atomcount * m_atomcount);
    m_VAB.resize(m_atomcount);

    for (int A = 0; A < m_atomcount; ++A) {
        int naoA = m_atom_nao[A];
        int npairA = nPairs(naoA);
        m_VAB[A] = Eigen::VectorXd::Zero(npairA);

        for (int B = 0; B < m_atomcount; ++B) {
            int naoB = m_atom_nao[B];
            int npairB = nPairs(naoB);

            if (A == B) {
                // One-center: fill from parametrized Gss/Gpp/Gsp/Gp2/Hsp
                Eigen::MatrixXd block = Eigen::MatrixXd::Zero(npairA, npairA);
                int Z = m_atoms[A];
                for (int mu_l = 0; mu_l < naoA; ++mu_l) {
                    for (int nu_l = mu_l; nu_l < naoA; ++nu_l) {
                        int p1 = pairIndex(nu_l, mu_l);
                        for (int lam_l = 0; lam_l < naoA; ++lam_l) {
                            for (int sig_l = lam_l; sig_l < naoA; ++sig_l) {
                                int p2 = pairIndex(sig_l, lam_l);
                                block(p1, p2) = oneCenterERI(Z, mu_l, nu_l, lam_l, sig_l);
                            }
                        }
                    }
                }
                m_gamma_blocks[A * m_atomcount + A] = block;
            } else if (B > A) {
                // Two-center electron-electron integrals (type=0)
                Eigen::MatrixXd block;
                computeIntegralBlock2C(A, B, block, 0);
                m_gamma_blocks[A * m_atomcount + B] = block;
                // Transpose for B,A access
                m_gamma_blocks[B * m_atomcount + A] = block.transpose();
            }

            if (B != A) {
                // Core (nuclear attraction) integrals for Hcore: type=1
                // VAB for atom A due to atom B: -Z_B * block(mu_nu_on_A, ss_on_B)
                Eigen::MatrixXd coreBlock;
                computeIntegralBlock2C(A, B, coreBlock, 1);
                int Z_B = m_atoms[B];
                int chgB = m_params->at(Z_B).n_valence;
                for (int p = 0; p < npairA; ++p) {
                    // ss on B is pair index 0
                    m_VAB[A](p) -= double(chgB) * coreBlock(p, 0);
                }
            }
        }
    }

    m_overlap = MakeOverlap(m_basis);
    m_hamiltonian = MakeH(m_overlap, m_basis);

    // DEBUG: Print key matrices for small molecules
    if (CurcumaLogger::get_verbosity() >= 3 && m_nbasis <= 10) {
        fmt::print("\n=== DEBUG: Overlap Matrix S ({0}x{0}) ===\n", m_nbasis);
        for (int i = 0; i < m_nbasis; ++i) {
            for (int j = 0; j < m_nbasis; ++j)
                fmt::print(" {:12.6f}", m_overlap(i,j));
            fmt::print("\n");
        }
        fmt::print("\n=== DEBUG: Core Hamiltonian H ({0}x{0}) ===\n", m_nbasis);
        for (int i = 0; i < m_nbasis; ++i) {
            for (int j = 0; j < m_nbasis; ++j)
                fmt::print(" {:12.6f}", m_hamiltonian(i,j));
            fmt::print("\n");
        }
        fmt::print("\n=== DEBUG: One-center gamma block for atom 0 ===\n");
        const auto& g0 = m_gamma_blocks[0];
        for (int i = 0; i < g0.rows(); ++i) {
            for (int j = 0; j < g0.cols(); ++j)
                fmt::print(" {:12.6f}", g0(i,j));
            fmt::print("\n");
        }
        fmt::print("\n=== DEBUG: VAB vectors ===\n");
        for (int A = 0; A < m_atomcount; ++A) {
            fmt::print("VAB[{}]: ", A);
            for (int p = 0; p < m_VAB[A].size(); ++p)
                fmt::print(" {:12.6f}", m_VAB[A](p));
            fmt::print("\n");
        }
        // Print two-center gamma blocks
        for (int A = 0; A < m_atomcount; ++A) {
            for (int B = A+1; B < m_atomcount; ++B) {
                const auto& gAB = m_gamma_blocks[A * m_atomcount + B];
                fmt::print("\n=== DEBUG: gamma_block[{},{}] ({}x{}) ===\n", A, B, gAB.rows(), gAB.cols());
                for (int i = 0; i < gAB.rows(); ++i) {
                    for (int j = 0; j < gAB.cols(); ++j)
                        fmt::print(" {:12.6f}", gAB(i,j));
                    fmt::print("\n");
                }
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("{} initialization complete", name));
    }

    return true;
}

double NDDO::Calculation(bool gradient)
{
    const char* name = NDDOParams::methodName(m_method_type);

    try {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::info(fmt::format("Starting {} energy calculation", name));
        }

        bool scf_success = runSCF();

        if (!scf_success) {
            CurcumaLogger::error(fmt::format("{} SCF did not converge", name));
            return 0.0;
        }

        m_energy_electronic = calculateElectronicEnergy();
        m_energy_core_repulsion = calculateCoreRepulsionEnergy();
        m_total_energy = m_energy_electronic + m_energy_core_repulsion;

        // DEBUG: Print converged density and Fock for small molecules
        if (CurcumaLogger::get_verbosity() >= 3 && m_nbasis <= 10) {
            fmt::print("\n=== DEBUG: Converged Density Matrix P ({0}x{0}) ===\n", m_nbasis);
            for (int i = 0; i < m_nbasis; ++i) {
                for (int j = 0; j < m_nbasis; ++j)
                    fmt::print(" {:12.6f}", m_density(i,j));
                fmt::print("\n");
            }
            fmt::print("\n=== DEBUG: Converged Fock Matrix F ({0}x{0}) ===\n", m_nbasis);
            for (int i = 0; i < m_nbasis; ++i) {
                for (int j = 0; j < m_nbasis; ++j)
                    fmt::print(" {:12.6f}", m_fock(i,j));
                fmt::print("\n");
            }
        }

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format("{} calculation complete", name));
            CurcumaLogger::energy_abs(m_total_energy, fmt::format("{}_total_energy", name));
        }

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

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("=== {} Algorithm Details (Level 3) ===", name));
            CurcumaLogger::param("basis_functions", m_nbasis);
            CurcumaLogger::param("electrons", m_num_electrons);
            CurcumaLogger::param("occupied_orbitals", m_num_electrons / 2);

            CurcumaLogger::info("All orbital energies:");
            for (int i = 0; i < std::min(m_nbasis, 10); ++i) {
                std::string occ_str = (i < m_num_electrons / 2) ? "occ" : "virt";
                CurcumaLogger::param(fmt::format("  MO[{}] ({})", i, occ_str),
                                    fmt::format("{:.4f} eV", m_energies(i) * eV2Eh));
            }

            double trace_P = m_density.trace();
            CurcumaLogger::param("density_matrix_trace", fmt::format("{:.6f}", trace_P));
            CurcumaLogger::info(fmt::format("=== End {} Details ===", name));
        }

        if (gradient) {
            m_gradient = calculateGradient();
        }

        return m_total_energy;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("{} calculation failed: {}", name, e.what()));
        return 0.0;
    }
}

// =================================================================================
// Basis Set Construction
// =================================================================================

int NDDO::buildBasisSet()
{
    m_basis.clear();
    m_num_electrons = 0;
    m_atom_ao_offset.resize(m_atoms.size());
    m_atom_nao.resize(m_atoms.size());

    int offset = 0;
    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];
        const NDDOParams::ElementParams& p = m_params->at(Z);

        m_atom_ao_offset[i] = offset;

        // s-orbital
        STO::Orbital s_orbital;
        s_orbital.type = STO::S;
        s_orbital.x = m_geometry(i, 0);
        s_orbital.y = m_geometry(i, 1);
        s_orbital.z = m_geometry(i, 2);
        s_orbital.zeta = p.zeta_s;
        s_orbital.VSIP = p.U_ss / eV2Eh;
        s_orbital.atom = i;
        s_orbital.principal_n = p.n_principal;
        m_basis.push_back(s_orbital);

        int nao = 1;
        // p-orbitals (Z > 2)
        if (Z > 2) {
            for (int m = 0; m < 3; ++m) {
                STO::Orbital p_orbital;
                p_orbital.type = (m == 0) ? STO::PX : ((m == 1) ? STO::PY : STO::PZ);
                p_orbital.x = m_geometry(i, 0);
                p_orbital.y = m_geometry(i, 1);
                p_orbital.z = m_geometry(i, 2);
                p_orbital.zeta = p.zeta_p;
                p_orbital.VSIP = p.U_pp / eV2Eh;
                p_orbital.atom = i;
                p_orbital.principal_n = p.n_principal;
                m_basis.push_back(p_orbital);
            }
            nao = 4;
        }

        m_atom_nao[i] = nao;
        offset += nao;
        m_num_electrons += p.n_valence;
    }

    return static_cast<int>(m_basis.size());
}

Matrix NDDO::MakeOverlap(std::vector<STO::Orbital>& basisset)
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
// Two-Center Integral Computation with Diatomic Frame Rotation
// Claude Generated: Following Ulysses OvRot.hpp (SPtransf, PPtransf) and MNDO.hpp (IntegralBlock2C)
// =================================================================================

double NDDO::getRho(int atomIdx, int l) const
{
    // Returns rho value in Angstrom for atom at index atomIdx and multipole level l
    // l=0: monopole (rho_s), l=1: dipole (rho_p), l=2: quadrupole (rho_q)
    int Z = m_atoms[atomIdx];
    const NDDOParams::ElementParams& p = m_params->at(Z);
    switch (l) {
        case 0: return p.rho_s;
        case 1: return p.rho_p;
        case 2: return p.rho_q;
        default: return p.rho_s;
    }
}

double NDDO::eri2Center(int L1, int M1, int L2, int M2,
                        int L3, int M3, int L4, int M4,
                        double R_AB, const std::vector<double>& D,
                        int atomA, int atomB, int core) const
{
    // Claude Generated: eri2Center following Ulysses MNDO.hpp lines 5917-5953
    // Computes two-center ERI in diatomic frame using NonZeroMultipoleMom decomposition
    //
    // L1,M1,L2,M2: orbital quantum numbers on atom A (bra pair)
    // L3,M3,L4,M4: orbital quantum numbers on atom B (ket pair)
    // core=0: electron-electron, core=1: use rho_core for B, core=2: use rho_core for A too

    double eri = 0.0;

    // Decompose bra pair (L1,M1)x(L2,M2) into multipole moments
    // Following NonZeroMultipoleMom_sp logic from Ulysses 2ElectronDewar.hpp
    auto getMultipoles = [](int LA, int MA, int LB, int MB) -> std::vector<std::pair<int,int>> {
        std::vector<std::pair<int,int>> result;
        int L1 = LA, M1 = MA, L2 = LB, M2 = MB;
        // Ensure L1 <= L2; if equal, ensure M1 <= M2
        if (L1 > L2) { std::swap(L1, L2); std::swap(M1, M2); }
        if (L1 == L2 && M1 > M2) { std::swap(M1, M2); }

        if (L1 == 0 && L2 == 0) {
            result.push_back({0, 0});  // monopole
        }
        else if (L1 == 0 && L2 == 1) {
            if (M2 == 0) {
                result.push_back({1, 0});  // dipole sigma
            } else {
                result.push_back({1, M2});  // dipole pi (sign from M2)
            }
        }
        else if (L1 == 1 && L2 == 1) {
            if (M1 == 0 && M2 == 0) {
                // sigma-sigma -> monopole + quadrupole_sigma
                result.push_back({0, 0});
                result.push_back({2, 0});
            }
            else if ((M1 == 0 && M2 == 1) || (M1 == 1 && M2 == 0)) {
                result.push_back({2, 1});  // quadrupole pi
            }
            else if ((M1 == 0 && M2 == -1) || (M1 == -1 && M2 == 0)) {
                result.push_back({2, -1});
            }
            else if ((M1 == 1 && M2 == 1) || (M1 == -1 && M2 == -1)) {
                // pi-pi or pib-pib -> monopole + quadrupole_delta
                result.push_back({0, 0});
                result.push_back({2, 2 * M2});  // sign from M2
            }
            else if ((M1 == 1 && M2 == -1) || (M1 == -1 && M2 == 1)) {
                // pi-pib -> quadrupole_xy
                result.push_back({2, 3});
            }
        }
        return result;
    };

    auto chg1 = getMultipoles(L1, M1, L2, M2);
    auto chg2 = getMultipoles(L3, M3, L4, M4);

    for (const auto& [lA, mA] : chg1) {
        for (const auto& [lB, mB] : chg2) {
            // Only even-even or odd-odd mA,mB contribute
            if ((std::abs(mA) % 2 == 0 && std::abs(mB) % 2 == 0) ||
                (std::abs(mA) % 2 != 0 && std::abs(mB) % 2 != 0)) {

                int mA_use = mA, mB_use = mB;
                if (std::abs(mA) != std::abs(mB)) {
                    mA_use = std::abs(mA);
                    mB_use = std::abs(mB);
                }

                // Get rho values based on multipole level
                double rhoA = getRho(atomA, lA) / au;  // Convert Angstrom -> Bohr
                double rhoB = getRho(atomB, lB) / au;
                if (core > 0) { rhoB = getRho(atomB, 0) / au; }      // Core: use rho_s for B
                if (core > 1 || core < 0) { rhoA = getRho(atomA, 0) / au; }

                // Handle m=3 (xy quadrupole): decompose into (2,2)-(2,-2) difference
                if (mA_use != 3) {
                    eri += curcuma::mndo::mndo_multipole_integral(lA, mA_use, lB, mB_use,
                                                                  R_AB, rhoA + rhoB, D);
                } else {
                    eri += 0.5 * (curcuma::mndo::mndo_multipole_integral(2, 2, 2, 2,
                                                                          R_AB, rhoA + rhoB, D)
                                - curcuma::mndo::mndo_multipole_integral(2, 2, 2, -2,
                                                                          R_AB, rhoA + rhoB, D));
                }
            }
        }
    }

    return eri;
}

void NDDO::computeIntegralBlock2C(int atomA, int atomB,
                                   Eigen::MatrixXd& block, int type) const
{
    // Claude Generated: Two-center integral block with diatomic frame rotation
    // Following Ulysses MNDO.hpp IntegralBlock2C (lines 1095-1217) and
    // OvRot.hpp SPtransf (lines 35-55), PPtransf (lines 884-910)
    //
    // Physics:
    // 1. Compute direction cosines for A->B bond direction
    // 2. Build SP rotation matrix (3x3): diatomic sigma/pi/pi_bar -> lab px/py/pz
    // 3. Build PP rotation matrix (6x6): from SP matrix via outer products
    // 4. Compute integrals in diatomic frame using eri2Center
    // 5. Rotate to lab frame

    int naoA = m_atom_nao[atomA];
    int naoB = m_atom_nao[atomB];
    int npairA = nPairs(naoA);
    int npairB = nPairs(naoB);

    block = Eigen::MatrixXd::Zero(npairA, npairB);

    // Step 1: Direction cosines
    double dx = m_geometry(atomB, 0) - m_geometry(atomA, 0);
    double dy = m_geometry(atomB, 1) - m_geometry(atomA, 1);
    double dz = m_geometry(atomB, 2) - m_geometry(atomA, 2);
    double R_ang = std::sqrt(dx*dx + dy*dy + dz*dz);
    double R_bohr = R_ang / au;  // Convert to Bohr for integral evaluation

    // Normalize direction vector
    double cost = dz / R_ang;
    double sint = std::sqrt(dx*dx + dy*dy) / R_ang;
    double cosp = 1.0, sinp = 0.0;
    const double tolerance = 1.0e-8;
    if (std::fabs(sint) > tolerance) {
        cosp = dx / (R_ang * sint);
        sinp = dy / (R_ang * sint);
    }

    // Step 2: SP rotation matrix (3x3)
    // Maps diatomic (sigma, pi_c, pi_s) -> lab (px, py, pz)
    // Row index: 0=px, 1=py, 2=pz; Col index: 0=sigma, 1=pi_c, 2=pi_s
    Eigen::Matrix3d SProt;
    SProt(0, 0) = sint * cosp;   SProt(0, 1) = cost * cosp;   SProt(0, 2) = -sinp;
    SProt(1, 0) = sint * sinp;   SProt(1, 1) = cost * sinp;   SProt(1, 2) = cosp;
    SProt(2, 0) = cost;          SProt(2, 1) = -sint;          SProt(2, 2) = 0.0;

    // Step 3: PP rotation matrix (6x6)
    // Row: orbital pair in lab (pxpx, pxpy, pxpz, pypy, pypz, pzpz)
    // Col: orbital pair in diatomic (sigma-sigma, sigma-pi_c, sigma-pi_s, pi_c-pi_c, pi_c-pi_s, pi_s-pi_s)
    // PProt(k, col) built from SP matrix via outer products
    Eigen::MatrixXd PProt = Eigen::MatrixXd::Zero(6, 6);
    {
        int counter = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j, ++counter) {
                int cnt = 0;
                for (int a = 0; a < 3; ++a) {
                    // Diagonal: SProt(i,a) * SProt(j,a)
                    PProt(counter, a) = SProt(i, a) * SProt(j, a);
                    for (int b = a + 1; b < 3; ++b, ++cnt) {
                        // Cross: SProt(i,a)*SProt(j,b) + SProt(i,b)*SProt(j,a)
                        PProt(counter, 3 + cnt) = SProt(i, a) * SProt(j, b) + SProt(i, b) * SProt(j, a);
                    }
                }
            }
        }
    }

    // Step 4: D-parameters for this atom pair (in Bohr)
    int Z_A = m_atoms[atomA];
    int Z_B = m_atoms[atomB];
    const NDDOParams::ElementParams& pA = m_params->at(Z_A);
    const NDDOParams::ElementParams& pB = m_params->at(Z_B);
    std::vector<double> D = { pA.D1 / au, pA.D2 / au, pB.D1 / au, pB.D2 / au };

    // Step 5: Compute diatomic-frame integrals and rotate to lab frame

    // (ss|ss) - always present
    block(0, 0) = eri2Center(0, 0, 0, 0, 0, 0, 0, 0, R_bohr, D, atomA, atomB, type);

    if (naoA > 1) {
        // === (sp|ss) integrals on atom A ===
        // Diatomic: (sp_sigma|ss)
        double sp_sigma_ss = eri2Center(0, 0, 1, 0, 0, 0, 0, 0, R_bohr, D, atomA, atomB, type);
        // Diatomic: (pp_sigma|ss) and (pp_pi|ss)
        double pp_sigma_ss = eri2Center(1, 0, 1, 0, 0, 0, 0, 0, R_bohr, D, atomA, atomB, type);
        double pp_pi_ss   = eri2Center(1, 1, 1, 1, 0, 0, 0, 0, R_bohr, D, atomA, atomB, type);

        for (int idx = 0; idx < 3; ++idx) {
            // (sp_i|ss): pair index for (s, p_i) = pairIndex(i+1, 0) where i+1 is local index of p_i
            int sp_pair = pairIndex(idx + 1, 0);
            block(sp_pair, 0) = -sp_sigma_ss * SProt(idx, 0);
        }

        for (int k = 0; k < 6; ++k) {
            // (p_i p_j|ss): pair index determined by pp pair ordering
            // k=0: pxpx(pair=pairIndex(1,1)), k=1: pxpy(pair=pairIndex(2,1)),
            // k=2: pxpz(pair=pairIndex(3,1)), k=3: pypy(pair=pairIndex(2,2)),
            // k=4: pypz(pair=pairIndex(3,2)), k=5: pzpz(pair=pairIndex(3,3))
            // But local indices are 1-based for p: local px=1, py=2, pz=3
            static const int pp_local_pairs[6][2] = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}};
            int pp_pair = pairIndex(pp_local_pairs[k][1], pp_local_pairs[k][0]);
            block(pp_pair, 0) = pp_sigma_ss * PProt(k, 0)
                              + pp_pi_ss * (PProt(k, 1) + PProt(k, 2));
        }
    }

    if (naoB > 1 && type == 0) {
        // === (ss|sp) and (ss|pp) integrals on atom B ===
        double ss_sp_sigma = eri2Center(0, 0, 0, 0, 0, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double ss_pp_sigma = eri2Center(0, 0, 0, 0, 1, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double ss_pp_pi    = eri2Center(0, 0, 0, 0, 1, 1, 1, 1, R_bohr, D, atomA, atomB, type);

        for (int idx = 0; idx < 3; ++idx) {
            int sp_pair = pairIndex(idx + 1, 0);
            block(0, sp_pair) = -ss_sp_sigma * SProt(idx, 0);
        }

        static const int pp_local_pairs[6][2] = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}};
        for (int k = 0; k < 6; ++k) {
            int pp_pair = pairIndex(pp_local_pairs[k][1], pp_local_pairs[k][0]);
            block(0, pp_pair) = ss_pp_sigma * PProt(k, 0)
                              + ss_pp_pi * (PProt(k, 1) + PProt(k, 2));
        }
    }

    if (naoA > 1 && naoB > 1 && type == 0) {
        // === (sp|sp) integrals ===
        double sp_sp_ss = eri2Center(0, 0, 1, 0, 0, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double sp_sp_pp = eri2Center(0, 0, 1, 1, 0, 0, 1, 1, R_bohr, D, atomA, atomB, type);

        for (int ib = 0; ib < 3; ++ib) {
            for (int ik = 0; ik < 3; ++ik) {
                int pairA = pairIndex(ib + 1, 0);
                int pairB = pairIndex(ik + 1, 0);
                block(pairA, pairB) = sp_sp_ss * SProt(ib, 0) * SProt(ik, 0)
                                    + sp_sp_pp * (SProt(ib, 1) * SProt(ik, 1) + SProt(ib, 2) * SProt(ik, 2));
            }
        }

        // === (sp|pp) integrals ===
        double sp_pp_ss    = eri2Center(0, 0, 1, 0, 1, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double sp_pp_pp    = eri2Center(0, 0, 1, 0, 1, 1, 1, 1, R_bohr, D, atomA, atomB, type);
        double sp_pp_cross = eri2Center(0, 0, 1, 1, 1, 0, 1, 1, R_bohr, D, atomA, atomB, type);

        static const int pp_local_pairs[6][2] = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}};
        for (int ib = 0; ib < 3; ++ib) {
            for (int ic = 0; ic < 6; ++ic) {
                int pairA = pairIndex(ib + 1, 0);
                int pairB = pairIndex(pp_local_pairs[ic][1], pp_local_pairs[ic][0]);
                block(pairA, pairB) = -SProt(ib, 0) * (sp_pp_ss * PProt(ic, 0) + sp_pp_pp * (PProt(ic, 1) + PProt(ic, 2)))
                                    - sp_pp_cross * (SProt(ib, 1) * PProt(ic, 3) + SProt(ib, 2) * PProt(ic, 4));
            }
        }

        // === (pp|sp) integrals ===
        double pp_sp_ss    = eri2Center(1, 0, 1, 0, 0, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double pp_sp_pp    = eri2Center(1, 1, 1, 1, 0, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double pp_sp_cross = eri2Center(1, 0, 1, 1, 0, 0, 1, 1, R_bohr, D, atomA, atomB, type);

        for (int ir = 0; ir < 6; ++ir) {
            for (int ic = 0; ic < 3; ++ic) {
                int pairA = pairIndex(pp_local_pairs[ir][1], pp_local_pairs[ir][0]);
                int pairB = pairIndex(ic + 1, 0);
                block(pairA, pairB) = -SProt(ic, 0) * (pp_sp_ss * PProt(ir, 0) + pp_sp_pp * (PProt(ir, 1) + PProt(ir, 2)))
                                    - pp_sp_cross * (SProt(ic, 1) * PProt(ir, 3) + SProt(ic, 2) * PProt(ir, 4));
            }
        }

        // === (pp|pp) integrals ===
        double pp_pp_ssss = eri2Center(1, 0, 1, 0, 1, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double pp_pp_ppss = eri2Center(1, 1, 1, 1, 1, 0, 1, 0, R_bohr, D, atomA, atomB, type);
        double pp_pp_sspp = eri2Center(1, 0, 1, 0, 1, 1, 1, 1, R_bohr, D, atomA, atomB, type);
        double pp_pp_pppp = eri2Center(1, 1, 1, 1, 1, 1, 1, 1, R_bohr, D, atomA, atomB, type);
        double pp_pp_psps = eri2Center(1, 1, 1, 0, 1, 1, 1, 0, R_bohr, D, atomA, atomB, type);
        double pp_pp_ppnn = eri2Center(1, 1, 1, 1, 1, -1, 1, -1, R_bohr, D, atomA, atomB, type);
        double pp_pp_npnp = eri2Center(1, -1, 1, 1, 1, -1, 1, 1, R_bohr, D, atomA, atomB, type);

        for (int ir = 0; ir < 6; ++ir) {
            for (int ic = 0; ic < 6; ++ic) {
                int pairA = pairIndex(pp_local_pairs[ir][1], pp_local_pairs[ir][0]);
                int pairB = pairIndex(pp_local_pairs[ic][1], pp_local_pairs[ic][0]);
                block(pairA, pairB) = pp_pp_ssss * PProt(ir, 0) * PProt(ic, 0)
                    + pp_pp_ppss * (PProt(ir, 1) + PProt(ir, 2)) * PProt(ic, 0)
                    + pp_pp_sspp * PProt(ir, 0) * (PProt(ic, 1) + PProt(ic, 2))
                    + pp_pp_pppp * (PProt(ir, 1) * PProt(ic, 1) + PProt(ir, 2) * PProt(ic, 2))
                    + pp_pp_psps * (PProt(ir, 3) * PProt(ic, 3) + PProt(ir, 4) * PProt(ic, 4))
                    + pp_pp_ppnn * (PProt(ir, 1) * PProt(ic, 2) + PProt(ir, 2) * PProt(ic, 1))
                    + pp_pp_npnp * PProt(ir, 5) * PProt(ic, 5);
            }
        }
    }
}

// =================================================================================
// One-Center Two-Electron Integrals from parametrized Gss/Gpp/Gsp/Gp2/Hsp
// =================================================================================

double NDDO::oneCenterERI(int Z, int mu_l, int nu_l, int lam_l, int sig_l) const
{
    // Claude Generated: One-center ERI using parametrized integrals
    // mu_l, nu_l, lam_l, sig_l are LOCAL orbital indices on one atom:
    //   0=s, 1=px, 2=py, 3=pz
    // Returns integral value in Hartree

    const NDDOParams::ElementParams& p = m_params->at(Z);

    auto is_s = [](int idx) { return idx == 0; };

    // (ss|ss) = Gss
    if (is_s(mu_l) && is_s(nu_l) && is_s(lam_l) && is_s(sig_l))
        return p.Gss / eV2Eh;

    // (ss|pp) = Gsp
    if (is_s(mu_l) && is_s(nu_l) && !is_s(lam_l) && lam_l == sig_l)
        return p.Gsp / eV2Eh;

    // (pp|ss) = Gsp
    if (!is_s(mu_l) && mu_l == nu_l && is_s(lam_l) && is_s(sig_l))
        return p.Gsp / eV2Eh;

    // (pp|pp) same direction = Gpp, different = Gp2
    if (!is_s(mu_l) && mu_l == nu_l && !is_s(lam_l) && lam_l == sig_l)
        return (mu_l == lam_l) ? p.Gpp / eV2Eh : p.Gp2 / eV2Eh;

    // (pp'|pp') = Hpp = 0.5*(Gpp - Gp2) — rotational invariance
    if (!is_s(mu_l) && !is_s(nu_l) && mu_l != nu_l &&
        !is_s(lam_l) && !is_s(sig_l) && lam_l != sig_l) {
        if ((mu_l == lam_l && nu_l == sig_l) || (mu_l == sig_l && nu_l == lam_l))
            return NDDOParams::getHpp(p) / eV2Eh;
    }

    // (sp|sp) = Hsp
    bool bra_sp = (is_s(mu_l) != is_s(nu_l));
    bool ket_sp = (is_s(lam_l) != is_s(sig_l));
    if (bra_sp && ket_sp) {
        int p_bra = is_s(mu_l) ? nu_l : mu_l;
        int p_ket = is_s(lam_l) ? sig_l : lam_l;
        if (p_bra == p_ket)
            return p.Hsp / eV2Eh;
    }

    return 0.0;
}

// =================================================================================
// Core Hamiltonian (NDDO approximation with proper off-diagonal one-center terms)
// =================================================================================

Matrix NDDO::MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset)
{
    // Claude Generated: NDDO core Hamiltonian with off-diagonal one-center nuclear attraction
    //
    // Diagonal one-center: H(mu,mu) = U_mu + VAB[pairIndex(mu_local, mu_local)]
    // Off-diagonal one-center (mu,nu on same atom, mu!=nu):
    //   H(mu,nu) = VAB[pairIndex(mu_local, nu_local)]
    //   This includes nuclear attraction from sp pairs — NOT zero!
    // Two-center: H(mu,nu) = beta_avg * S(mu,nu)
    //
    // Reference: Dewar & Thiel, Theor. Chim. Acta 46, 89 (1977), Eq. 11-13

    Matrix H = Matrix::Zero(basisset.size(), basisset.size());

    for (size_t mu = 0; mu < basisset.size(); ++mu) {
        int atom_i = basisset[mu].atom;
        int Z_i = m_atoms[atom_i];
        int offset_i = m_atom_ao_offset[atom_i];
        int mu_local = static_cast<int>(mu) - offset_i;

        for (size_t nu = 0; nu <= mu; ++nu) {
            int atom_j = basisset[nu].atom;

            if (atom_i == atom_j) {
                int nu_local = static_cast<int>(nu) - offset_i;
                int pair_idx = pairIndex(mu_local, nu_local);

                if (mu == nu) {
                    // Diagonal: U_mu + nuclear attraction from all other atoms
                    double H_diag = getCoreIntegral(Z_i, basisset[mu].type);
                    H_diag += m_VAB[atom_i](pair_idx);
                    H(mu, nu) = H_diag;
                } else {
                    // Off-diagonal one-center: nuclear attraction contribution
                    // This is CRUCIAL for sp pairs — provides the sp mixing
                    H(mu, nu) = H(nu, mu) = m_VAB[atom_i](pair_idx);
                }
            } else {
                // Two-center: resonance integral H_mu_nu = beta_avg * S_mu_nu
                double beta = getResonanceIntegral(basisset[mu], basisset[nu], 0.0);
                H(mu, nu) = H(nu, mu) = beta * S(mu, nu);
            }
        }
    }

    return H;
}

double NDDO::getCoreIntegral(int element, int orbital_type) const
{
    const NDDOParams::ElementParams& p = m_params->at(element);

    if (orbital_type == STO::S) {
        return p.U_ss / eV2Eh;
    } else if (orbital_type == STO::PX || orbital_type == STO::PY || orbital_type == STO::PZ) {
        return p.U_pp / eV2Eh;
    }

    return 0.0;
}

double NDDO::getResonanceIntegral(const STO::Orbital& fi, const STO::Orbital& fj, double /*distance*/) const
{
    // NDDO resonance integral: beta_mu_nu = (beta_mu + beta_nu) / 2
    // Reference: Dewar & Thiel, Theor. Chim. Acta 46, 89 (1977), Eq. 13
    int Z_i = m_atoms[fi.atom];
    int Z_j = m_atoms[fj.atom];

    const NDDOParams::ElementParams& pi = m_params->at(Z_i);
    const NDDOParams::ElementParams& pj = m_params->at(Z_j);

    double beta_i = (fi.type == STO::S) ? pi.beta_s : pi.beta_p;
    double beta_j = (fj.type == STO::S) ? pj.beta_s : pj.beta_p;

    double beta_avg = (beta_i + beta_j) / 2.0;

    return beta_avg / eV2Eh;
}

// =================================================================================
// SCF Procedure
// =================================================================================

bool NDDO::runSCF()
{
    // Claude Generated: NDDO SCF with Pulay DIIS acceleration
    // Reference: P. Pulay, Chem. Phys. Lett. 73, 393 (1980)
    //
    // DIIS error vector (ZDO/S=I basis): e = F*P - P*F
    // Stores last diis_max Fock matrices and error vectors
    // Starts DIIS after diis_start simple iterations
    const char* name = NDDOParams::methodName(m_method_type);
    m_density = Matrix::Zero(m_nbasis, m_nbasis);

    constexpr int diis_max = 8;    // DIIS history depth
    constexpr int diis_start = 3;  // iterations before DIIS kicks in

    std::vector<Matrix> diis_fock;
    std::vector<Matrix> diis_error;

    Matrix S_identity = Matrix::Identity(m_nbasis, m_nbasis);

    for (int iter = 0; iter < m_scf_max_iterations; ++iter) {
        m_fock = buildFockMatrix(m_density);

        // DIIS: build error vector e = F*P - P*F, extrapolate Fock matrix
        Matrix fock_use = m_fock;
        if (iter >= diis_start && !diis_fock.empty()) {
            // Error vector: commutator [F, P] in the ZDO (S=I) basis
            Matrix err = m_fock * m_density - m_density * m_fock;

            diis_fock.push_back(m_fock);
            diis_error.push_back(err);

            // Trim history
            if ((int)diis_fock.size() > diis_max) {
                diis_fock.erase(diis_fock.begin());
                diis_error.erase(diis_error.begin());
            }

            int n = static_cast<int>(diis_fock.size());
            // Build DIIS B matrix: B_ij = tr(e_i * e_j)
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n + 1, n + 1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double bij = (diis_error[i].array() * diis_error[j].array()).sum();
                    B(i, j) = B(j, i) = bij;
                }
                B(i, n) = B(n, i) = -1.0;
            }
            B(n, n) = 0.0;

            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
            rhs(n) = -1.0;

            // Solve via Cholesky decomposition (B is symmetric)
            Eigen::VectorXd c = B.ldlt().solve(rhs);

            // Check solution validity
            bool valid = c.head(n).allFinite() && c.head(n).norm() < 1e6;
            if (valid) {
                fock_use = Matrix::Zero(m_nbasis, m_nbasis);
                for (int i = 0; i < n; ++i) {
                    fock_use += c(i) * diis_fock[i];
                }
            }
        } else {
            // Pre-DIIS: store for history
            if (iter == diis_start - 1) {
                Matrix err = m_fock * m_density - m_density * m_fock;
                diis_fock.push_back(m_fock);
                diis_error.push_back(err);
            }
        }

        // NDDO uses standard eigenvalue problem F*C = ε*C (ZDO: S = I)
        ParallelEigenSolver solver(500, 128, 1.0e-10, false);
        solver.setThreadCount(m_threads);

        bool success = solver.solve(S_identity, fock_use, m_energies, m_mo, m_threads, false);

        if (!success) {
            CurcumaLogger::error(fmt::format("Eigenvalue solution failed in {} SCF", name));
            return false;
        }

        Matrix density_new = buildDensityMatrix(m_mo, m_energies);
        double delta_P = (density_new - m_density).norm();

        if (CurcumaLogger::get_verbosity() >= 3) {
            double E_elec_iter = 0.0;
            for (int mu = 0; mu < m_nbasis; ++mu) {
                for (int nu = 0; nu < m_nbasis; ++nu) {
                    E_elec_iter += density_new(mu, nu) * (m_hamiltonian(mu, nu) + m_fock(mu, nu));
                }
            }
            E_elec_iter /= 2.0;

            CurcumaLogger::param(fmt::format("SCF_iter_{}", iter+1),
                               fmt::format("dP = {:.6e}, E_elec = {:.6f} Eh", delta_P, E_elec_iter));
        }

        if (delta_P < m_scf_threshold) {
            m_density = density_new;
            m_scf_converged = true;

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::success(fmt::format("{} SCF converged in {} iterations", name, iter+1));
            }

            return true;
        }

        // Simple damping for density update (pre-DIIS and as fallback)
        m_density = m_scf_damping * m_density + (1.0 - m_scf_damping) * density_new;
    }

    CurcumaLogger::warn(fmt::format("{} SCF did not converge in {} iterations", name, m_scf_max_iterations));
    m_scf_converged = false;
    return false;
}

// =================================================================================
// Fock Matrix using precomputed gamma blocks
// Claude Generated: Following Ulysses MNDO.hpp calcFock (lines 588-675)
// =================================================================================

Matrix NDDO::buildFockMatrix(const Matrix& density)
{
    // NDDO Fock matrix: F = H + G
    //
    // One-center contributions (atoms A, orbitals mu,nu on A):
    //   Diagonal (mu==nu):
    //     F(mu,mu) += sum_{nu,lam on A} P(nu,lam) * [gamma(mu_mu, nu_lam) - 0.5*gamma(mu_nu, mu_lam)]
    //   Off-diagonal (mu!=nu):
    //     F(mu,nu) += sum_{lam,sig on A} P(lam,sig) * [gamma(mu_nu, lam_sig) - 0.5*gamma(mu_lam, nu_sig)]
    //
    // Two-center contributions (A != B):
    //   Coulomb on A: F(mu,nu) += sum_{lam,sig on B} P(lam,sig) * gamma_AB(mu_nu, lam_sig)
    //   Exchange: F(mu,lam) -= 0.5 * sum_{nu on A, sig on B} P(nu,sig) * gamma_AB(mu_nu, lam_sig)
    //
    // gamma_AB indexed by pair indices: gamma_AB(pairIndex(mu_l,nu_l), pairIndex(lam_l,sig_l))

    Matrix F = m_hamiltonian;

    for (int A = 0; A < m_atomcount; ++A) {
        int offA = m_atom_ao_offset[A];
        int naoA = m_atom_nao[A];

        // One-center Fock contributions
        const Eigen::MatrixXd& gammaAA = m_gamma_blocks[A * m_atomcount + A];

        // One-center Fock: following Ulysses MNDO.hpp calcFock lines 604-637
        // F(mu,mu) gets contributions from ALL (nu, lam) pairs on same atom
        // F(mu,nu) with mu!=nu gets additional contributions from (lam, sig) pairs
        for (int imu = 0; imu < naoA; ++imu) {
            int mu = offA + imu;
            int pos_mu_mu = pairIndex(imu, imu);

            for (int inu = 0; inu < naoA; ++inu) {
                int nu = offA + inu;
                int pos_mu_nu = pairIndex(imu, inu);

                for (int ilam = 0; ilam < naoA; ++ilam) {
                    int lam = offA + ilam;
                    int pos_nu_lam = pairIndex(inu, ilam);
                    int pos_mu_lam = pairIndex(imu, ilam);

                    // Always accumulate diagonal: F(mu,mu) += P(nu,lam) * [gamma(mu_mu,nu_lam) - 0.5*gamma(mu_nu,mu_lam)]
                    F(mu, mu) += density(nu, lam) * (gammaAA(pos_mu_mu, pos_nu_lam)
                                                   - 0.5 * gammaAA(pos_mu_nu, pos_mu_lam));

                    // Off-diagonal: F(mu,nu) += P(lam,sig) * [gamma(mu_nu,lam_sig) - 0.5*gamma(mu_lam,nu_sig)]
                    if (imu != inu) {
                        for (int isig = 0; isig < naoA; ++isig) {
                            int sig = offA + isig;
                            int pos_lam_sig = pairIndex(ilam, isig);
                            int pos_nu_sig = pairIndex(inu, isig);
                            F(mu, nu) += density(lam, sig) * (gammaAA(pos_mu_nu, pos_lam_sig)
                                                            - 0.5 * gammaAA(pos_mu_lam, pos_nu_sig));
                        }
                    }
                }
            }
        }

        // Two-center Fock contributions
        for (int B = 0; B < A; ++B) {
            int offB = m_atom_ao_offset[B];
            int naoB = m_atom_nao[B];

            const Eigen::MatrixXd& gammaAB = m_gamma_blocks[A * m_atomcount + B];
            const Eigen::MatrixXd& gammaBA = m_gamma_blocks[B * m_atomcount + A];

            for (int imu = 0; imu < naoA; ++imu) {
                int mu = offA + imu;

                for (int inu = 0; inu < naoA; ++inu) {
                    int nu = offA + inu;
                    int pos_mu_nu = pairIndex(imu, inu);

                    for (int ilam = 0; ilam < naoB; ++ilam) {
                        int lam = offB + ilam;

                        for (int isig = 0; isig < naoB; ++isig) {
                            int sig = offB + isig;
                            int pos_lam_sig = pairIndex(ilam, isig);

                            // Coulomb on atom A from atom B density
                            // F(mu,nu) on A += P(lam,sig) on B * gamma_AB(mu_nu, lam_sig)
                            F(mu, nu) += density(lam, sig) * gammaAB(pos_mu_nu, pos_lam_sig);

                            // Coulomb on atom B from atom A density
                            // F(lam,sig) on B += P(mu,nu) on A * gamma_BA(lam_sig, mu_nu)
                            F(lam, sig) += density(mu, nu) * gammaBA(pos_lam_sig, pos_mu_nu);

                            // Exchange: F(mu,lam) -= 0.5 * P(nu,sig) * gamma_AB(mu_nu, lam_sig)
                            F(mu, lam) -= 0.5 * density(nu, sig) * gammaAB(pos_mu_nu, pos_lam_sig);

                            // Exchange: F(lam,mu) -= 0.5 * P(sig,nu) * gamma_BA(lam_sig, mu_nu)
                            F(lam, mu) -= 0.5 * density(sig, nu) * gammaBA(pos_lam_sig, pos_mu_nu);
                        }
                    }
                }
            }
        }
    }

    return F;
}

Matrix NDDO::buildDensityMatrix(const Matrix& mo_coefficients, const Vector& /*mo_energies*/)
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

double NDDO::calculateElectronicEnergy() const
{
    // E_elec = Tr(P * (H + F)) / 2
    double E = 0.0;

    for (int mu = 0; mu < m_nbasis; ++mu) {
        for (int nu = 0; nu < m_nbasis; ++nu) {
            E += m_density(mu, nu) * (m_hamiltonian(mu, nu) + m_fock(mu, nu));
        }
    }

    return E / 2.0;
}

double NDDO::calculateCoreRepulsionEnergy() const
{
    // Claude Generated: NDDO core-core repulsion
    // MNDO: E_core = Z_A * Z_B * gamma_ss_core * factorA  (no Gaussians)
    // AM1/PM3/PM6: E_core += Z_A * Z_B * factorC  (Gaussian corrections)
    // gamma_ss_core uses rho_core (= rho_s) for both atoms
    // Reference: Ulysses MNDO.hpp line 5880

    const bool use_gaussians = NDDOParams::hasGaussianCorrections(m_method_type);

    double E_rep = 0.0;

    for (int A = 0; A < m_atomcount; ++A) {
        int Z_A = m_atoms[A];
        const NDDOParams::ElementParams& pA = m_params->at(Z_A);
        int chgA = pA.n_valence;

        for (int B = A + 1; B < m_atomcount; ++B) {
            int Z_B = m_atoms[B];
            const NDDOParams::ElementParams& pB = m_params->at(Z_B);
            int chgB = pB.n_valence;

            double dx = m_geometry(A, 0) - m_geometry(B, 0);
            double dy = m_geometry(A, 1) - m_geometry(B, 1);
            double dz = m_geometry(A, 2) - m_geometry(B, 2);
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);
            double R_bohr = R_AB / au;

            // gamma_ss with core rho (both atoms use rho_s) — type=2 in Ulysses
            std::vector<double> D = { pA.D1 / au, pA.D2 / au, pB.D1 / au, pB.D2 / au };
            double gamma_ss = eri2Center(0, 0, 0, 0, 0, 0, 0, 0, R_bohr, D, A, B, 2);

            // factorA: Exponential damping with method-specific g-factor
            // MNDO only: special N-H/O-H treatment (gfactor returns R_Ang for N/O paired with H)
            // AM1/PM3/PM6: gfactor always returns 1.0 — see Ulysses AM1::gfactor (line 7022),
            //              PM3::gfactor (line 8308): both override to return 1.0
            // Reference: Ulysses MNDO.hpp line 5393 (MNDO::gfactor), 7022 (AM1), 8308 (PM3)
            double gA = 1.0;
            double gB = 1.0;
            const bool use_special_gfactor = (m_method_type == NDDOMethodType::MNDO);
            if (use_special_gfactor) {
                if ((Z_B == 1) && ((Z_A == 7) || (Z_A == 8))) {
                    gA = R_AB;
                }
                if ((Z_A == 1) && ((Z_B == 7) || (Z_B == 8))) {
                    gB = R_AB;
                }
            }
            double factorA = 1.0;
            factorA += gA * std::exp(-pA.alpha * R_AB);
            factorA += gB * std::exp(-pB.alpha * R_AB);

            // Gaussian corrections (AM1/PM3/PM6 only)
            // Following Ulysses AM1factor(): RAB is in Bohr internally, AM1factor converts
            // to Angstrom via *dist_Angstrom2au (=0.5292). K values are in eV·Å, L in Å⁻²,
            // M in Å. Formula: factorC += K/(au2eV*R_ang) * exp(-L*(R_ang-M)^2)
            // R_AB here is already in Angstrom (from molecular geometry).
            double factorC = 0.0;
            if (use_gaussians) {
                constexpr double au2eV = 27.211383473452294;
                for (size_t k = 0; k < pA.gauss_K.size(); ++k) {
                    factorC += pA.gauss_K[k]
                               * std::exp(-pA.gauss_L[k] * std::pow(R_AB - pA.gauss_M[k], 2.0))
                               / (au2eV * R_AB);
                }
                for (size_t k = 0; k < pB.gauss_K.size(); ++k) {
                    factorC += pB.gauss_K[k]
                               * std::exp(-pB.gauss_L[k] * std::pow(R_AB - pB.gauss_M[k], 2.0))
                               / (au2eV * R_AB);
                }
            }

            double ZZ = double(chgA * chgB);
            E_rep += ZZ * gamma_ss * factorA + ZZ * factorC;
        }
    }

    return E_rep;
}

// =================================================================================
// Gradient Calculation (numerical)
// =================================================================================

Matrix NDDO::calculateGradient() const
{
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Calculating {} gradients (numerical)",
                                        NDDOParams::methodName(m_method_type)));
    }

    Matrix gradient = Matrix::Zero(m_atomcount, 3);
    const double delta = 1.0e-5;

    Matrix geom_orig = m_geometry;
    double E0 = m_total_energy;

    for (int atom = 0; atom < m_atomcount; ++atom) {
        for (int coord = 0; coord < 3; ++coord) {
            const_cast<NDDO*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) + delta;
            const_cast<NDDO*>(this)->InitialiseMolecule();
            double E_plus = const_cast<NDDO*>(this)->Calculation(false);

            const_cast<NDDO*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord) - delta;
            const_cast<NDDO*>(this)->InitialiseMolecule();
            double E_minus = const_cast<NDDO*>(this)->Calculation(false);

            gradient(atom, coord) = (E_plus - E_minus) / (2.0 * delta);

            const_cast<NDDO*>(this)->m_geometry(atom, coord) = geom_orig(atom, coord);
        }
    }

    const_cast<NDDO*>(this)->m_geometry = geom_orig;
    const_cast<NDDO*>(this)->InitialiseMolecule();
    const_cast<NDDO*>(this)->m_total_energy = E0;

    return gradient / au;
}

// =================================================================================
// Property Access
// =================================================================================

json NDDO::getEnergyDecomposition() const
{
    json decomp;
    decomp["electronic"] = m_energy_electronic;
    decomp["core_repulsion"] = m_energy_core_repulsion;
    decomp["total"] = m_total_energy;
    return decomp;
}

double NDDO::getHOMOLUMOGap() const
{
    return getLUMOEnergy() - getHOMOEnergy();
}

double NDDO::getHOMOEnergy() const
{
    int n_occ = m_num_electrons / 2;
    if (n_occ == 0 || n_occ > m_nbasis) return 0.0;
    return m_energies(n_occ - 1) * eV2Eh;
}

double NDDO::getLUMOEnergy() const
{
    int n_occ = m_num_electrons / 2;
    if (n_occ >= m_nbasis) return 0.0;
    return m_energies(n_occ) * eV2Eh;
}

double NDDO::getHeatOfFormation() const
{
    // TODO: Calculate heat of formation from atomization energy
    return 0.0;
}
