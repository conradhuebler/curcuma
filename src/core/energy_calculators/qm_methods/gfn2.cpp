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
#include "diis_accelerator.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "ParallelEigenSolver.hpp"

#include <fmt/format.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <vector>
#include <map>

using namespace CurcumaUnit;

// =================================================================================
// Covalent Radii (Pyykkö 2015, triple bond radii in Ångström)
// Shared by GFN2, GFN1, and other xTB methods
// =================================================================================
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

double getCovalentRadius(int Z) {
    if (Z < 1 || Z > 86) return 1.5;
    return COVALENT_RADII[Z];
}

GFN2::GFN2()
    : m_params()
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_repulsion(0.0)
    , m_energy_coulomb(0.0)
    , m_energy_dispersion(0.0)
    , m_energy_solvation(0.0)
    , m_solvation(nullptr)
    , m_solvent("none")
    , m_scf_max_iterations(150)      // GFN2 needs more iterations (TBLite default)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.6)             // GFN2 needs stronger damping (TBLite default)
    , m_scf_converged(false)
{
    if (!m_param_db.loadCompleteGFN2()) {
        CurcumaLogger::error("Failed to load GFN2 parameter database");
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Initializing native GFN2-xTB method");
        CurcumaLogger::param("method", "GFN2-xTB");
        CurcumaLogger::param("reference", "Bannwarth et al. JCTC 2019, 15, 1652");
        CurcumaLogger::param("parameter_database", fmt::format("{} elements, {} pairs",
                           m_param_db.getNumElements(), m_param_db.getNumPairs()));
        CurcumaLogger::param("scf_damping", fmt::format("{:.2f}", m_scf_damping));
        CurcumaLogger::param("max_iterations", m_scf_max_iterations);
    }
}

GFN2::GFN2(const ArrayParameters& params)
    : m_params(params)
    , m_nbasis(0)
    , m_energy_electronic(0.0)
    , m_energy_repulsion(0.0)
    , m_energy_coulomb(0.0)
    , m_energy_dispersion(0.0)
    , m_energy_solvation(0.0)
    , m_solvation(nullptr)
    , m_solvent("none")
    , m_scf_max_iterations(150)      // GFN2 needs more iterations (TBLite default)
    , m_scf_threshold(1.0e-6)
    , m_scf_damping(0.6)             // GFN2 needs stronger damping (TBLite default)
    , m_scf_converged(false)
{
    if (!m_param_db.loadCompleteGFN2()) {
        CurcumaLogger::error("Failed to load GFN2 parameter database");
    }
}

bool GFN2::InitialiseMolecule()
{
    if (m_atoms.size() == 0) return false;

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        if (!m_params.isValidAtom(m_atoms[i])) return false;
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing GFN2 calculation");
        CurcumaLogger::param("atoms", static_cast<int>(m_atoms.size()));
        CurcumaLogger::param("charge", m_charge);
        CurcumaLogger::param("spin", m_spin);
    }

    m_nbasis = buildBasisSet();
    if (m_nbasis == 0) return false;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("basis_functions", m_nbasis);
        CurcumaLogger::param("electrons", m_num_electrons);
    }

    m_coordination_numbers = calculateCoordinationNumbers();

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
    if (!m_initialised && m_atoms.size() == 0) return 0.0;

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Starting GFN2-xTB calculation");
    }

    try {
        m_overlap = MakeOverlap(m_basis);
        m_hamiltonian = MakeH(m_overlap, m_basis);
        m_scf_converged = runSCF();

        if (!m_scf_converged) CurcumaLogger::warn("SCF did not converge");

        m_energy_electronic = calculateElectronicEnergy();
        m_energy_repulsion = calculateRepulsionEnergy();
        m_energy_coulomb = calculateCoulombEnergy();
        m_energy_dispersion = calculateDispersionEnergy();

        m_energy_solvation = 0.0;
        if (m_solvation && m_solvent != "none") {
            std::vector<std::array<double, 3>> positions(m_atomcount);
            for (int i = 0; i < m_atomcount; ++i) positions[i] = {m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2)};
            std::vector<double> charges(m_atomcount);
            for (int i = 0; i < m_atomcount; ++i) charges[i] = m_charges(i);
            m_energy_solvation = m_solvation->calculateEnergy(m_atoms, positions, charges);
        }

        m_total_energy = m_energy_electronic + m_energy_repulsion + m_energy_coulomb + m_energy_dispersion + m_energy_solvation;

        // Debug: dump density + orbital info at verbosity 3 for tblite comparison
        if (CurcumaLogger::get_verbosity() >= 3) {
            json st;
            st["method"] = "gfn2";
            st["natoms"] = m_atomcount;
            st["nbasis"] = m_nbasis;
            st["atoms"] = json::array();
            for (int i = 0; i < m_atomcount; ++i)
                st["atoms"].push_back({{"z", m_atoms[i]},
                                       {"xyz_bohr", {m_geometry(i,0), m_geometry(i,1), m_geometry(i,2)}}});
            st["energy_total"] = m_total_energy;
            st["atomic_charges"] = json::array();
            for (int i = 0; i < m_atomcount; ++i) st["atomic_charges"].push_back(m_charges(i));
            st["orbital_energies"] = json::array();
            for (int i = 0; i < m_nbasis; ++i) st["orbital_energies"].push_back(m_energies(i));
            st["density"] = json::array();
            for (int i = 0; i < m_nbasis; ++i) {
                json row = json::array();
                for (int j = 0; j < m_nbasis; ++j) row.push_back(m_density(i, j));
                st["density"].push_back(std::move(row));
            }
            auto dip = calculateAtomicDipoles();
            auto qua = calculateAtomicQuadrupoles();
            st["dpat_simple"] = json::array();
            st["qpat_simple"] = json::array();
            for (int i = 0; i < m_atomcount; ++i) {
                st["dpat_simple"].push_back({dip[i][0], dip[i][1], dip[i][2]});
                st["qpat_simple"].push_back({qua[i][0], qua[i][1], qua[i][2], qua[i][3], qua[i][4], qua[i][5]});
            }
            std::ofstream ofs("gfn2_state.json");
            if (ofs) ofs << st.dump(2);
            CurcumaLogger::info("GFN2 state dumped to gfn2_state.json");
        }

        if (CurcumaLogger::get_verbosity() >= 1) CurcumaLogger::energy_abs(m_total_energy, "GFN2 Total Energy");

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", m_energy_electronic));
            CurcumaLogger::param("repulsion", fmt::format("{:.6f} Eh", m_energy_repulsion));
            CurcumaLogger::param("coulomb", fmt::format("{:.6f} Eh", m_energy_coulomb - m_energy_aes2));
            CurcumaLogger::param("aes2", fmt::format("{:.6f} Eh", m_energy_aes2));
#ifdef USE_D4
            CurcumaLogger::param("dispersion_d4", fmt::format("{:.6f} Eh", m_energy_dispersion));
#else
            CurcumaLogger::param("dispersion", fmt::format("{:.6f} Eh (D4 not compiled)", m_energy_dispersion));
#endif
            CurcumaLogger::param("HOMO", fmt::format("{:.4f} eV", getHOMOEnergy()));
            CurcumaLogger::param("LUMO", fmt::format("{:.4f} eV", getLUMOEnergy()));
            CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} eV", getHOMOLUMOGap()));
        }

        if (gradient) {
            m_gradient = calculateGradient();
            if (m_solvation && m_solvent != "none") {
                std::vector<std::array<double, 3>> positions(m_atomcount), g_solv(m_atomcount);
                std::vector<double> charges(m_atomcount);
                for (int i = 0; i < m_atomcount; ++i) {
                    positions[i] = {m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2)};
                    charges[i] = m_charges(i);
                    g_solv[i] = {0.0, 0.0, 0.0};
                }
                m_solvation->calculateEnergyAndGradients(m_atoms, positions, charges, g_solv);
                for (int i = 0; i < m_atomcount; ++i) for (int k = 0; k < 3; ++k) m_gradient(i, k) += g_solv[i][k];
            }
        }

        return m_total_energy;
    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("GFN2 calculation failed: {}", e.what()));
        return 0.0;
    }
}

int GFN2::buildBasisSet()
{
    m_basis.clear();
    m_num_electrons = 0;

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        int Z = m_atoms[i];
        int n_shells = (Z > 12) ? 3 : (Z > 2) ? 2 : 1;

        for (int shell = 0; shell < n_shells; ++shell) {
            double zeta = 0.0;
            if (m_param_db.hasElement(Z) && m_param_db.getElement(Z).shells.count(shell)) {
                zeta = m_param_db.getElement(Z).shells.at(shell).zeta;
            } else zeta = m_params.getYeff(Z, shell);

            if (zeta < 1.0e-6) continue;

            int l = shell;
            for (int m = -l; m <= l; ++m) {
                STO::Orbital orbital;
                if (l == 0) orbital.type = STO::S;
                else if (l == 1) {
                    if (m == -1) orbital.type = STO::PY;
                    else if (m == 0) orbital.type = STO::PZ;
                    else orbital.type = STO::PX;
                } else orbital.type = STO::S; // Placeholder for d

                orbital.x = m_geometry(i, 0) / au;
                orbital.y = m_geometry(i, 1) / au;
                orbital.z = m_geometry(i, 2) / au;
                orbital.zeta = zeta;
                orbital.atom = i;
                orbital.shell = shell;
                m_basis.push_back(orbital);
            }

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

Matrix GFN2::MakeOverlap(std::vector<STO::Orbital>& basisset)
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

Matrix GFN2::MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset)
{
    Matrix H = Matrix::Zero(basisset.size(), basisset.size());
    for (size_t i = 0; i < basisset.size(); ++i) {
        int atom_i = basisset[i].atom;
        double CN_i = m_coordination_numbers(atom_i);
        H(i, i) = getSelfEnergy(m_atoms[atom_i], basisset[i].shell, CN_i);
    }
    for (size_t i = 0; i < basisset.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            double dx = basisset[i].x - basisset[j].x, dy = basisset[i].y - basisset[j].y, dz = basisset[i].z - basisset[j].z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz) * au;
            double scale = getHamiltonianScale(basisset[i], basisset[j], distance);
            double h_ij = scale * S(i, j) * (H(i, i) + H(j, j)) / 2.0;
            H(i, j) = H(j, i) = h_ij;
        }
    }
    return H;
}

double GFN2::getSelfEnergy(int element, int shell, double CN) const
{
    if (m_param_db.hasElement(element)) {
        const auto& ep = m_param_db.getElement(element);
        if (ep.shells.count(shell)) return ep.shells.at(shell).selfenergy + ep.shells.at(shell).kcn * CN;
    }
    return -m_params.getHardness(element) + m_params.getShellHardness(element) * 0.01 * CN;
}

double GFN2::getShellHubbard(int Z, int shell) const
{
    double g = m_param_db.getShellHubbard(Z, shell);
    return (g > 0.0) ? g : m_params.getHardness(Z);
}

// Claude Generated (March 2026): Fixed to match TBLite gfn2.f90 lines 1002-1006
// Uses global kshell from kdiag, kpair=1.0 for all pairs
double GFN2::getHamiltonianScale(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const
{
    int Zi = m_atoms[fi.atom], Zj = m_atoms[fj.atom];
    double z_ij = std::sqrt(2.0 * std::sqrt(fi.zeta * fj.zeta) / (fi.zeta + fj.zeta));

    // kpair = 1.0 for all pairs in GFN2 (TBLite gfn2.f90 line 727)
    double k_pair = 1.0;

    // kshell from kdiag array (TBLite gfn2.f90 line 61, kshell function line 1002)
    // Special case: s-d and p-d interactions use 2.0
    int li = fi.shell, lj = fj.shell;
    double k_shell;
    if ((li == 2 && (lj == 0 || lj == 1)) || (lj == 2 && (li == 0 || li == 1))) {
        k_shell = 2.0;  // s-d and p-d special case
    } else {
        constexpr double kdiag[] = {GFN2Params::GFN2_KDIAG_S, GFN2Params::GFN2_KDIAG_P, GFN2Params::GFN2_KDIAG_D};
        k_shell = (kdiag[li] + kdiag[lj]) / 2.0;
    }

    double delta_EN = m_params.getElectronegativity(Zi) - m_params.getElectronegativity(Zj);
    double en_factor = 1.0 + 0.02 * delta_EN * delta_EN;

    double r_bohr = distance / au;
    double rad_i = m_param_db.hasElement(Zi) ? m_param_db.getElement(Zi).rad : 2.0;
    double rad_j = m_param_db.hasElement(Zj) ? m_param_db.getElement(Zj).rad : 2.0;
    double shpoly_i = m_param_db.hasElement(Zi) ? m_param_db.getElement(Zi).shells.at(fi.shell).shpoly : 0.0;
    double shpoly_j = m_param_db.hasElement(Zj) ? m_param_db.getElement(Zj).shells.at(fj.shell).shpoly : 0.0;
    double rr = std::sqrt(r_bohr / (rad_i + rad_j));
    double poly_r = (1.0 + shpoly_i * rr) * (1.0 + shpoly_j * rr);

    return z_ij * k_pair * k_shell * en_factor * poly_r;
}

// Claude Generated (March 2026): Fixed to match TBLite ncoord/gfn.f90 lines 82-101
// Uses double-exponential CN (ka*kb product), NOT D3 single-exponential^(4/3)
Vector GFN2::calculateCoordinationNumbers()
{
    Vector CN = Vector::Zero(m_atomcount);
    // r_shift is in Bohr in TBLite; convert to Angstrom for our distance units
    double r_shift_ang = GFN2Params::GFN2_CN_RSHIFT * au;
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            double R_AB = (m_geometry.row(A) - m_geometry.row(B)).norm();  // Angstrom
            if (R_AB < 1.0e-6) continue;
            double R_cov = getCovalentRadius(m_atoms[A]) + getCovalentRadius(m_atoms[B]);  // Angstrom

            // GFN2 double-exponential: count1 * count2
            double count1 = 1.0 / (1.0 + std::exp(-GFN2Params::GFN2_CN_KA * (R_cov / R_AB - 1.0)));
            double count2 = 1.0 / (1.0 + std::exp(-GFN2Params::GFN2_CN_KB * ((R_cov + r_shift_ang) / R_AB - 1.0)));
            CN(A) += count1 * count2;
        }
    }
    return CN;
}

bool GFN2::runSCF()
{
    // Claude Generated (March 2026): GFN2 SCF with DIIS acceleration
    // DIIS (Pulay 1980) is essential for convergence on molecules with >5 atoms.
    // Without it, simple damping cannot converge C6H6, caffeine, etc.

    m_density = Matrix::Zero(m_nbasis, m_nbasis);
    double prev_energy = 0.0;
    double energy_diff = 0.0;

    // Compute S^(-1/2) once before SCF loop
    ParallelEigenSolver solver(500, 128, 1.0e-10, false);
    solver.setThreadCount(m_threads);
    if (!solver.computeS_1_2(m_overlap, m_S_inv_sqrt, m_threads)) {
        CurcumaLogger::error("Failed to compute S^(-1/2) for GFN2 SCF");
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

        bool success = solver.solveWithPrecalculatedS_1_2(m_S_inv_sqrt, F_scf, m_energies, m_mo, m_threads, true);
        if (!success) {
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::warn("SCF: Eigenvalue solver failed at iteration " + std::to_string(iter + 1));
            return false;
        }

        Matrix density_new = buildDensityMatrix(m_mo, m_energies);

        double current_energy = 0.5 * (density_new.cwiseProduct(m_hamiltonian + m_fock)).sum();

        if (iter > 0) energy_diff = std::abs(current_energy - prev_energy);
        prev_energy = current_energy;

        // Check both density and energy convergence
        double density_change = (density_new - m_density).norm();
        bool converged = (density_change < m_scf_threshold) && (energy_diff < m_scf_threshold * 0.1);

        if (CurcumaLogger::get_verbosity() >= 2 && (iter % 10 == 0 || converged || iter == 0)) {
            CurcumaLogger::info(fmt::format("SCF iter {:3d}: E = {:15.8f} Eh, ΔE = {:8.2e}, ΔP = {:8.2e}",
                                          iter + 1, current_energy, energy_diff, density_change));
        }

        if (converged) {
            m_density = density_new;
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::success(fmt::format("SCF converged in {} iterations (ΔE = {:.2e})", iter + 1, energy_diff));
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

    if (CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("SCF did not converge after {} iterations (ΔE = {:.2e})",
                                       m_scf_max_iterations, energy_diff));
    return false;
}

Matrix GFN2::buildFockMatrix(const Matrix& density)
{
    Matrix F = m_hamiltonian;
    Matrix PS = density * m_overlap;
    std::vector<std::map<int, double>> shell_pop(m_atomcount);
    for (int mu = 0; mu < m_nbasis; ++mu) shell_pop[m_basis[mu].atom][m_basis[mu].shell] += PS(mu, mu);

    Vector atomic_charges = Vector::Zero(m_atomcount);
    std::vector<std::map<int, double>> dq(m_atomcount);
    for (int A = 0; A < m_atomcount; ++A) {
        const auto& ep = m_param_db.getElement(m_atoms[A]);
        for (auto const& [s, p] : shell_pop[A]) {
            dq[A][s] = ep.shells.at(s).refocc - p;
            atomic_charges(A) += dq[A][s];
        }
    }
    const_cast<GFN2*>(this)->m_charges = atomic_charges;
    if (density.norm() < 1e-10) return F;

    std::vector<double> V(m_nbasis, 0.0);
    for (int mu = 0; mu < m_nbasis; ++mu) {
        int A = m_basis[mu].atom;
        double g_AA = getShellHubbard(m_atoms[A], m_basis[mu].shell);
        for (int B = 0; B < m_atomcount; ++B) {
            double R_AB = (A == B) ? 0.0 : (m_geometry.row(A) - m_geometry.row(B)).norm() / au;
            for (auto const& [sB, dqB] : dq[B]) {
                double g_AB = calculateCoulombKernel(g_AA, getShellHubbard(m_atoms[B], sB), R_AB);
                V[mu] += dqB * g_AB;
            }
        }
    }
    for (int mu = 0; mu < m_nbasis; ++mu) for (int nu = 0; nu < m_nbasis; ++nu) F(mu, nu) -= 0.5 * m_overlap(mu, nu) * (V[mu] + V[nu]);
    return F;
}

Matrix GFN2::buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies)
{
    // Claude Generated (March 2026): Two occupation schemes controlled by m_electronic_temperature
    //   etemp == 0  → integer closed-shell (2/0 occupation)
    //   etemp > 0   → Fermi-Dirac smearing (fractional occupation)
    int nbas = mo_coefficients.rows();

    if (m_electronic_temperature <= 0.0) {
        int n_occ = m_num_electrons / 2;
        auto C_occ = mo_coefficients.leftCols(n_occ);
        return 2.0 * C_occ * C_occ.transpose();
    }

    const double kT = m_electronic_temperature * 3.166808e-6;  // K → Hartree
    int n_elec = m_num_electrons;

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

double GFN2::calculateElectronicEnergy() const {
    // GFN2 electronic energy: E_el = Tr(P * H_0)
    // Claude Generated (March 2026): Vectorized with Eigen cwiseProduct
    return (m_density.cwiseProduct(m_hamiltonian)).sum();
}

// Claude Generated (March 2026): Fixed to match TBLite effective.f90 lines 162-183
// Bugs fixed: Z_A*Z_B (not sum), sqrt(alpha) (not arithmetic mean), R^kexp
double GFN2::calculateRepulsionEnergy() const {
    double E_rep = 0.0;
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = A + 1; B < m_atomcount; ++B) {
            double R_AB = (m_geometry.row(A) - m_geometry.row(B)).norm() / au;
            const auto &ea = m_param_db.getElement(m_atoms[A]), &eb = m_param_db.getElement(m_atoms[B]);

            double alpha_AB = std::sqrt(ea.rep_alpha * eb.rep_alpha);  // geometric mean
            double zeff_AB = ea.rep_zeff * eb.rep_zeff;                // product, not sum

            // kexp = 1.5 for heavy pairs, 1.0 when both atoms Z <= 2
            double kexp = (m_atoms[A] > 2 || m_atoms[B] > 2)
                          ? GFN2Params::GFN2_REP_KEXP : GFN2Params::GFN2_REP_KEXP_LIGHT;

            double r1k = std::pow(R_AB, kexp);
            double r1r = std::pow(R_AB, GFN2Params::GFN2_REP_REXP);
            E_rep += zeff_AB * std::exp(-alpha_AB * r1k) / r1r;
        }
    }
    return E_rep;
}

double GFN2::calculateCoulombEnergy() const {
    double E_ES2 = 0.0, E_ES3 = 0.0;
    Matrix PS = m_density * m_overlap;
    std::vector<std::map<int, double>> dq(m_atomcount);
    Vector q(m_atomcount);
    for (int mu = 0; mu < m_nbasis; ++mu) {
        int A = m_basis[mu].atom;
        dq[A][m_basis[mu].shell] += PS(mu, mu);
    }
    for (int A = 0; A < m_atomcount; ++A) {
        const auto& ep = m_param_db.getElement(m_atoms[A]);
        for (auto& [s, p] : dq[A]) { p = ep.shells.at(s).refocc - p; q(A) += p; }
        E_ES3 += (1.0/6.0) * ep.hubbard_deriv * std::pow(q(A), 3.0);
    }
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            double R_AB = (A == B) ? 0.0 : (m_geometry.row(A) - m_geometry.row(B)).norm() / au;
            for (auto const& [sA, dqA] : dq[A]) for (auto const& [sB, dqB] : dq[B]) {
                double g_AB = calculateCoulombKernel(getShellHubbard(m_atoms[A], sA), getShellHubbard(m_atoms[B], sB), R_AB);
                E_ES2 += 0.5 * dqA * dqB * g_AB;
            }
        }
    }

    auto dipoles = calculateAtomicDipoles();
    auto quadrupoles = calculateAtomicQuadrupoles();
    const_cast<GFN2*>(this)->m_energy_aes2 = calculateAES2Energy(dipoles, quadrupoles);

    return E_ES2 + E_ES3 + m_energy_aes2;
}

double GFN2::calculateDispersionEnergy() const {
#ifdef USE_D4
    // Claude Generated (March 2026): DFT-D4 dispersion for GFN2
    // Reference: E. Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122
    // GFN2 parameters from TBLite gfn2.f90 line 54:
    //   s6 = 1.0, s8 = 2.7, a1 = 0.52, a2 = 5.0
    //
    // Note: This is a const method, so we need to cast away constness for the D4 interface
    // which modifies internal state during calculation.

    // Get non-const access to D4 interface (lazy initialization)
    GFN2* nonconst_this = const_cast<GFN2*>(this);

    if (!nonconst_this->m_d4) {
        nonconst_this->m_d4 = std::make_unique<DFTD4Interface>();

        // Set GFN2 dispersion parameters
        // These are the standard GFN2-xTB dispersion parameters
        json d4_params;
        d4_params["s6"] = GFN2_D4_S6;        // 1.0
        d4_params["s8"] = GFN2_D4_S8;        // 2.7
        d4_params["a1"] = GFN2_D4_A1;        // 0.52
        d4_params["a2"] = GFN2_D4_A2;        // 5.0
        d4_params["alpha"] = 16.0;           // D4 standard
        d4_params["functional"] = "gfn2";    // GFN2 specific

        nonconst_this->m_d4->UpdateParameters(d4_params);
    }

    // Initialize molecule for D4 calculation
    Mol d4_mol;
    d4_mol.m_number_atoms = m_atomcount;
    d4_mol.m_atoms = m_atoms;
    d4_mol.setGeometry(m_geometry);

    if (!nonconst_this->m_d4->InitialiseMolecule(d4_mol)) {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::warn("D4: Failed to initialize molecule for dispersion calculation");
        }
        return 0.0;
    }

    // Set charge for D4 calculation
    nonconst_this->m_d4->setCharge(m_charge);

    // Calculate D4 dispersion energy
    double dispersion_energy = nonconst_this->m_d4->Calculation(false);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dispersion_d4", fmt::format("{:.6f} Eh", dispersion_energy));
    }

    return dispersion_energy;
#else
    // Claude Generated 2026: Native D4 via the standalone D4Evaluator
    // (src/core/energy_calculators/dispersion/d4_evaluator.h). Uses the
    // same Casimir-Polder C6 reference data that GFN-FF carries, but
    // with GFN2's own scaling parameters (s8=2.7 instead of GFN-FF's 2.0)
    // — the BJ damping math is identical, just the coefficients differ.
    //
    // GFN2-xTB D4 parameters (TBLite/Ulysses-confirmed): s6=1.0, s8=2.7,
    // a1=0.52, a2=5.0, s9=5.0, alpha=16.0.
    //
    // Caveat: zeta-charge scaling (zetac6) is currently treated as a static
    // prefactor as in GFN-FF. A full GFN2-style q-response chain rule
    // (dE/dq · dq/dx via SCF response) is deferred to a follow-up AP — the
    // expected residual is sub-mEh vs TBLite. See dispersion/CLAUDE.md.

    GFN2* nonconst_this = const_cast<GFN2*>(this);

    if (!nonconst_this->m_d4_native) {
        json d4_config;
        d4_config["d4_s6"] = GFN2Params::GFN2_D4_S6;   // 1.0
        d4_config["d4_s8"] = GFN2Params::GFN2_D4_S8;   // 2.7
        d4_config["d4_a1"] = GFN2Params::GFN2_D4_A1;   // 0.52
        d4_config["d4_a2"] = GFN2Params::GFN2_D4_A2;   // 5.0
        d4_config["d4_alp"] = 16.0;                    // GFN2 (Caldeweyher 2019)
        ConfigManager config("d4param", d4_config);
        nonconst_this->m_d4_native = std::make_unique<D4ParameterGenerator>(config);
    }

    if (!nonconst_this->m_d4_evaluator) {
        curcuma::dispersion::D4Params p;
        p.s6 = GFN2Params::GFN2_D4_S6;
        p.s8 = GFN2Params::GFN2_D4_S8;
        p.a1 = GFN2Params::GFN2_D4_A1;
        p.a2 = GFN2Params::GFN2_D4_A2;
        p.s9 = 5.0;                                    // Ulysses GFN2 ATM scaling
        p.alpha = 16.0;
        p.damping = curcuma::dispersion::DampingFormula::StandardBJ_D4;
        nonconst_this->m_d4_evaluator = std::make_unique<curcuma::dispersion::D4Evaluator>(
            nonconst_this->m_d4_native.get(), p);
    }

    // D4 expects geometry in Bohr
    Matrix geometry_bohr = m_geometry / au;
    nonconst_this->m_d4_native->GenerateParameters(m_atoms, geometry_bohr);

    // Compute dispersion gradient + dE/dCN regardless of whether the caller
    // requested a gradient: GFN2's calculateGradient() folds these in if
    // present, and the cost is dominated by the energy loop anyway. The
    // cached Matrix/Vector are reused across SCF cycles for the same geometry.
    // Legacy GFN2 class: q-response not wired here (uses static prefactor).
    // dEdq_unused only binds the reference; with_dEdq defaults to false.
    Vector dEdq_unused;
    double disp_energy = nonconst_this->m_d4_evaluator->computeEnergyAndGradient(
        m_atoms, geometry_bohr,
        /*with_gradient=*/true,
        nonconst_this->m_disp_gradient,
        nonconst_this->m_disp_dEdcn,
        dEdq_unused);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dispersion_d4_native", fmt::format("{:.6f} Eh", disp_energy));
    }

    return disp_energy;
#endif
}

double GFN2::calculateCoulombKernel(double gA, double gB, double R) const {
    return 1.0 / std::sqrt(std::pow(R, 2.0) + std::pow(2.0 / (gA + gB), 2.0));
}

Matrix GFN2::calculateGradient() const
{
    std::vector<Eigen::Vector3d> dipoles = calculateAtomicDipoles();
    std::vector<Eigen::VectorXd> quadrupoles = calculateAtomicQuadrupoles();

    Matrix gradient = Matrix::Zero(m_atomcount, 3);
    Matrix W = Matrix::Zero(m_nbasis, m_nbasis);
    int n_occ = m_num_electrons / 2;
    for (int i = 0; i < n_occ; ++i) {
        double eps = m_energies(i);
        for (int mu = 0; mu < m_nbasis; ++mu) for (int nu = 0; nu < m_nbasis; ++nu)
            W(mu, nu) += 2.0 * eps * m_mo(mu, i) * m_mo(nu, i);
    }

    // Prepare AES2 parameters
    Vector rad(m_atomcount);
    for (int i = 0; i < m_atomcount; ++i) {
        if (m_param_db.hasElement(m_atoms[i]))
            rad(i) = m_param_db.getElement(m_atoms[i]).multipole_rad;
        else
            rad(i) = 2.0;
    }

    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = A + 1; B < m_atomcount; ++B) {
            double dx = m_geometry(B, 0) - m_geometry(A, 0), dy = m_geometry(B, 1) - m_geometry(A, 1), dz = m_geometry(B, 2) - m_geometry(A, 2);
            double R_AB_ang = std::sqrt(dx*dx + dy*dy + dz*dz), R_AB_bohr = R_AB_ang / au;
            if (R_AB_ang < 1e-10) continue;
            double nx = dx / R_AB_ang, ny = dy / R_AB_ang, nz = dz / R_AB_ang;

            // Claude Generated (March 2026): Fixed repulsion gradient to match TBLite effective.f90:232
            // dG = -(alpha*r1k*kexp + rexp) * dE * rij/r2
            const auto &ea = m_param_db.getElement(m_atoms[A]), &eb = m_param_db.getElement(m_atoms[B]);
            double alpha_AB = std::sqrt(ea.rep_alpha * eb.rep_alpha);
            double zeff_AB = ea.rep_zeff * eb.rep_zeff;
            double kexp = (m_atoms[A] > 2 || m_atoms[B] > 2)
                          ? GFN2Params::GFN2_REP_KEXP : GFN2Params::GFN2_REP_KEXP_LIGHT;
            double rexp = GFN2Params::GFN2_REP_REXP;
            double r1k = std::pow(R_AB_bohr, kexp);
            double r1r = std::pow(R_AB_bohr, rexp);
            double dE = zeff_AB * std::exp(-alpha_AB * r1k) / r1r;
            double dG_scalar = -(alpha_AB * r1k * kexp + rexp) * dE / (R_AB_bohr * R_AB_bohr) * R_AB_bohr;
            // dG = dG_scalar * (rij / r2) but we have unit vector nx,ny,nz = rij/R_AB
            // and R_AB_bohr = R_AB_ang/au, so dG_scalar already has rij/r2 factor via the formula
            // Actually: dG = -(alpha*r1k*kexp + rexp) * dE * rij/r2
            // rij/r2 = (rij/R) * (1/R) = n_hat / R_AB_bohr
            double dVdR = -(alpha_AB * r1k * kexp + rexp) * dE / (R_AB_bohr * R_AB_bohr);
            gradient(A, 0) += dVdR * nx / au; gradient(A, 1) += dVdR * ny / au; gradient(A, 2) += dVdR * nz / au;
            gradient(B, 0) -= dVdR * nx / au; gradient(B, 1) -= dVdR * ny / au; gradient(B, 2) -= dVdR * nz / au;

            for (size_t mu = 0; mu < m_basis.size(); ++mu) {
                if (m_basis[mu].atom != A) continue;
                for (size_t nu = 0; nu < m_basis.size(); ++nu) {
                    if (m_basis[nu].atom != B) continue;
                    double dS_dR = STO::calculateOverlapDerivative(m_basis[mu], m_basis[nu], R_AB_bohr);
                    double h_avg = (m_hamiltonian(mu, mu) + m_hamiltonian(nu, nu)) / 2.0;
                    int Zi = m_atoms[A], Zj = m_atoms[B];
                    double rad_i = ea.rad, rad_j = eb.rad;
                    double sh_i = ea.shells.at(m_basis[mu].shell).shpoly, sh_j = eb.shells.at(m_basis[nu].shell).shpoly;
                    double rr = std::sqrt(R_AB_bohr / (rad_i + rad_j));
                    double drr = 1.0 / (2.0 * std::sqrt(R_AB_bohr * (rad_i + rad_j)));
                    double dpoly = sh_i * drr * (1.0 + sh_j * rr) + (1.0 + sh_i * rr) * sh_j * drr;
                    double z_ij = std::sqrt(2.0 * std::sqrt(m_basis[mu].zeta * m_basis[nu].zeta) / (m_basis[mu].zeta + m_basis[nu].zeta));
                    double dEN = m_params.getElectronegativity(Zi) - m_params.getElectronegativity(Zj);
                    // Claude Generated (March 2026): Use global kshell, kpair=1.0
                    double k_p = 1.0;
                    int li_g = m_basis[mu].shell, lj_g = m_basis[nu].shell;
                    double k_s;
                    if ((li_g == 2 && (lj_g == 0 || lj_g == 1)) || (lj_g == 2 && (li_g == 0 || li_g == 1))) {
                        k_s = 2.0;
                    } else {
                        constexpr double kdiag[] = {GFN2Params::GFN2_KDIAG_S, GFN2Params::GFN2_KDIAG_P, GFN2Params::GFN2_KDIAG_D};
                        k_s = (kdiag[li_g] + kdiag[lj_g]) / 2.0;
                    }
                    double C = z_ij * k_p * k_s * (1.0 + 0.02 * dEN * dEN);
                    double dH0 = (C * dpoly * m_overlap(mu, nu) + C * (1.0 + sh_i * rr) * (1.0 + sh_j * rr) * dS_dR) * h_avg;
                    double contrib = 2.0 * (m_density(mu, nu) * dH0 - W(mu, nu) * dS_dR);
                    gradient(A, 0) += (contrib / au) * nx; gradient(A, 1) += (contrib / au) * ny; gradient(A, 2) += (contrib / au) * nz;
                    gradient(B, 0) -= (contrib / au) * nx; gradient(B, 1) -= (contrib / au) * ny; gradient(B, 2) -= (contrib / au) * nz;
                }
            }
            double qA = m_charges(A), qB = m_charges(B);
            double gA = getShellHubbard(m_atoms[A], 0), gB = getShellHubbard(m_atoms[B], 0);
            double dg = -std::pow(std::pow(R_AB_bohr, 2.0) + std::pow(2.0/(gA+gB), 2.0), -1.5) * R_AB_bohr;
            double dES2 = (qA * qB * dg) / au;
            gradient(A, 0) += dES2 * nx; gradient(A, 1) += dES2 * ny; gradient(A, 2) += dES2 * nz;
            gradient(B, 0) -= dES2 * nx; gradient(B, 1) -= dES2 * ny; gradient(B, 2) -= dES2 * nz;

            // Claude Generated (March 2026): Fixed CN gradient for double-exponential
            // d(count1*count2)/dR = dcount1*count2 + count1*dcount2
            double R_cov = getCovalentRadius(m_atoms[A]) + getCovalentRadius(m_atoms[B]);
            double r_shift_ang = GFN2Params::GFN2_CN_RSHIFT * au;
            double rc2 = R_cov + r_shift_ang;
            double exp1 = std::exp(-GFN2Params::GFN2_CN_KA * (R_cov / R_AB_ang - 1.0));
            double exp2 = std::exp(-GFN2Params::GFN2_CN_KB * (rc2 / R_AB_ang - 1.0));
            double count1 = 1.0 / (1.0 + exp1);
            double count2 = 1.0 / (1.0 + exp2);
            double dcount1 = (-GFN2Params::GFN2_CN_KA * R_cov * exp1) / (R_AB_ang * R_AB_ang * (exp1 + 1.0) * (exp1 + 1.0));
            double dcount2 = (-GFN2Params::GFN2_CN_KB * rc2 * exp2) / (R_AB_ang * R_AB_ang * (exp2 + 1.0) * (exp2 + 1.0));
            double dCN = dcount1 * count2 + count1 * dcount2;
            double kA = ea.shells.at(0).kcn, kB = eb.shells.at(0).kcn;
            double sPA = 0, sPB = 0;
            for (size_t mu = 0; mu < m_basis.size(); ++mu) { if (m_basis[mu].atom == A) sPA += m_density(mu, mu); if (m_basis[mu].atom == B) sPB += m_density(mu, mu); }
            double dECN = (sPA * kA + sPB * kB) * dCN;
            gradient(A, 0) += dECN * nx; gradient(A, 1) += dECN * ny; gradient(A, 2) += dECN * nz;
            gradient(B, 0) -= dECN * nx; gradient(B, 1) -= dECN * ny; gradient(B, 2) -= dECN * nz;

            // =====================================================================
            // AES2 Gradient Contributions (Claude Generated March 2026)
            // Reference: TBLite multipole.f90 get_multipole_gradient
            // =====================================================================

            // 1. On-site dipole kernel gradient (dE/dR from dkernel * |mu|^2)
            // This term is zero because dkernel is atom-specific and R-independent

            // 2. Charge-Dipole gradient: E = sum_{A!=B} q_A * mu_B . grad(gamma_AB)
            // dE/dR_AB = q_A * sum_i mu_B^i * d(dgamma^i/dR)
            double qA_val = qA;
            const auto& dipA = dipoles[A], dipB = dipoles[B];
            const auto& quadB = quadrupoles[B];
            (void)quadB; // Suppress unused warning (used in charge-quadrupole section)

            // Charge-Dipole: d/dR [q_A * mu_B . (R_AB / R_AB^3) * f_dmp3]
            double rr3 = R_AB_bohr / (rad(A) + rad(B));
            double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr3, 8.0));
            double dFdR_cd = qA_val / (R_AB_bohr * R_AB_bohr * R_AB_bohr); // d/dR (1/R^2)
            // Simplified: only radial derivative
            double dE_dR_cd = qA_val * dipB.dot(Eigen::Vector3d(nx, ny, nz)) * dFdR_cd * fdmp3;
            gradient(A, 0) += dE_dR_cd * nx; gradient(A, 1) += dE_dR_cd * ny; gradient(A, 2) += dE_dR_cd * nz;
            gradient(B, 0) -= dE_dR_cd * nx; gradient(B, 1) -= dE_dR_cd * ny; gradient(B, 2) -= dE_dR_cd * nz;

            // 3. Dipole-Dipole gradient: E = 0.5 * sum_{A!=B} mu_A . T_AB . mu_B
            // T_AB = dipole-dipole interaction tensor
            double dd_prefactor = 0.5 * fdmp3; // Use same damping for simplicity
            double T_trace = (dipA.dot(dipB) - 3.0 * dipA.dot(Eigen::Vector3d(nx, ny, nz)) * dipB.dot(Eigen::Vector3d(nx, ny, nz))) / (R_AB_bohr * R_AB_bohr * R_AB_bohr);
            double dE_dR_dd = dd_prefactor * T_trace * (-3.0 / R_AB_bohr);
            gradient(A, 0) += dE_dR_dd * nx; gradient(A, 1) += dE_dR_dd * ny; gradient(A, 2) += dE_dR_dd * nz;
            gradient(B, 0) -= dE_dR_dd * nx; gradient(B, 1) -= dE_dR_dd * ny; gradient(B, 2) -= dE_dR_dd * nz;

            // 4. Charge-Quadrupole gradient: E = sum_{A!=B} q_A * Theta_B : grad(grad(gamma_AB))
            // Simplified: only isotropic contribution
            if (quadB.size() > 0) {
                double rr5 = R_AB_bohr / (rad(A) + rad(B));
                double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr5, 10.0));
                double dE_dR_cq = qA_val * (quadB[0] + quadB[3] + quadB[5]) * fdmp5 / (R_AB_bohr * R_AB_bohr * R_AB_bohr * R_AB_bohr);
                gradient(A, 0) += dE_dR_cq * nx; gradient(A, 1) += dE_dR_cq * ny; gradient(A, 2) += dE_dR_cq * nz;
                gradient(B, 0) -= dE_dR_cq * nx; gradient(B, 1) -= dE_dR_cq * ny; gradient(B, 2) -= dE_dR_cq * nz;
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

json GFN2::getEnergyDecomposition() const {
    json d; d["electronic"] = m_energy_electronic; d["repulsion"] = m_energy_repulsion; d["coulomb"] = m_energy_coulomb - m_energy_aes2; d["aes2"] = m_energy_aes2; d["total"] = m_total_energy; d["scf_converged"] = m_scf_converged; return d;
}

double GFN2::getHOMOLUMOGap() const {
    if (m_energies.size() < 2 || m_num_electrons < 2) return 0;
    int h = m_num_electrons/2 - 1; return (m_energies(h+1) - m_energies(h)) * eV2Eh;
}

double GFN2::getHOMOEnergy() const { if (m_energies.size() == 0) return 0; return m_energies(m_num_electrons/2 - 1) * eV2Eh; }
double GFN2::getLUMOEnergy() const { if (m_energies.size() == 0) return 0; return m_energies(m_num_electrons/2) * eV2Eh; }

Matrix GFN2::calculateNumericalGradient(double h) const
{
    // Claude Generated (March 2026): Numerical gradient for testing analytical gradients
    // Reference: Finite difference method, central difference formula

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Calculating numerical gradient with h = {:.4f} Ang", h));
    }

    Matrix num_gradient = Matrix::Zero(m_atomcount, 3);

    // Cast away constness to temporarily modify geometry (const method for API convenience)
    GFN2* nonconst_this = const_cast<GFN2*>(this);
    Geometry original_geometry = nonconst_this->m_geometry;

    for (int A = 0; A < m_atomcount; ++A) {
        for (int k = 0; k < 3; ++k) {
            // Displace atom in +direction
            nonconst_this->m_geometry(A, k) = original_geometry(A, k) + h;
            double E_plus = nonconst_this->Calculation(false);

            // Displace atom in -direction
            nonconst_this->m_geometry(A, k) = original_geometry(A, k) - h;
            double E_minus = nonconst_this->Calculation(false);

            // Central difference: dE/dx ≈ (E(x+h) - E(x-h)) / (2h)
            // Convert from Hartree/Ang to Hartree/Bohr
            num_gradient(A, k) = (E_plus - E_minus) / (2.0 * h) / au;

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Atom {} dim {}: E+ = {:.8f}, E- = {:.8f}, dE/dx = {:.6f} Eh/Bohr",
                                              A+1, k, E_plus, E_minus, num_gradient(A, k)));
            }
        }
    }

    // Restore original geometry
    nonconst_this->m_geometry = original_geometry;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("numerical_gradient_norm", fmt::format("{:.6e} Eh/Ang", num_gradient.norm()));
    }

    return num_gradient;
}

std::vector<Eigen::Vector3d> GFN2::calculateAtomicDipoles() const
{
    std::vector<Eigen::Vector3d> dipoles(m_atomcount, Eigen::Vector3d::Zero());
    for (int A = 0; A < m_atomcount; ++A) {
        if (!m_param_db.hasElement(m_atoms[A])) continue;
        const auto& ep = m_param_db.getElement(m_atoms[A]);
        int s_idx = -1, p_idx[3] = {-1, -1, -1};
        for (int mu = 0; mu < m_nbasis; ++mu) {
            if (m_basis[mu].atom != A) continue;
            if (m_basis[mu].type == STO::S) s_idx = mu;
            else if (m_basis[mu].type == STO::PX) p_idx[0] = mu;
            else if (m_basis[mu].type == STO::PY) p_idx[1] = mu;
            else if (m_basis[mu].type == STO::PZ) p_idx[2] = mu;
        }
        if (s_idx != -1 && ep.dkernel > 1e-6) {
            for (int i = 0; i < 3; ++i) {
                if (p_idx[i] != -1) {
                    dipoles[A][i] = 2.0 * ep.dkernel * m_density(s_idx, p_idx[i]);
                }
            }
        }
    }
    return dipoles;
}

std::vector<Eigen::VectorXd> GFN2::calculateAtomicQuadrupoles() const
{
    std::vector<Eigen::VectorXd> quadrupoles(m_atomcount, Eigen::VectorXd::Zero(6));
    for (int A = 0; A < m_atomcount; ++A) {
        if (!m_param_db.hasElement(m_atoms[A])) continue;
        const auto& ep = m_param_db.getElement(m_atoms[A]);
        int p_idx[3] = {-1, -1, -1};
        for (int mu = 0; mu < m_nbasis; ++mu) {
            if (m_basis[mu].atom != A) continue;
            if (m_basis[mu].type == STO::PX) p_idx[0] = mu;
            else if (m_basis[mu].type == STO::PY) p_idx[1] = mu;
            else if (m_basis[mu].type == STO::PZ) p_idx[2] = mu;
        }

        if (ep.qkernel > 1e-6 && p_idx[0] != -1 && p_idx[1] != -1 && p_idx[2] != -1) {
            double pxx = m_density(p_idx[0], p_idx[0]);
            double pyy = m_density(p_idx[1], p_idx[1]);
            double pzz = m_density(p_idx[2], p_idx[2]);
            double pxy = m_density(p_idx[0], p_idx[1]);
            double pxz = m_density(p_idx[0], p_idx[2]);
            double pyz = m_density(p_idx[1], p_idx[2]);

            // GFN2 Quadrupole moments (Bannwarth 2019 Eq. 11 & TBLite)
            // Order: xx, xy, xz, yy, yz, zz
            quadrupoles[A][0] = ep.qkernel * (2.0 * pxx - pyy - pzz); // xx
            quadrupoles[A][1] = ep.qkernel * 3.0 * pxy;               // xy
            quadrupoles[A][2] = ep.qkernel * 3.0 * pxz;               // xz
            quadrupoles[A][3] = ep.qkernel * (2.0 * pyy - pxx - pzz); // yy
            quadrupoles[A][4] = ep.qkernel * 3.0 * pyz;               // yz
            quadrupoles[A][5] = ep.qkernel * (2.0 * pzz - pxx - pyy); // zz
        }
    }
    return quadrupoles;
}

std::vector<Matrix> GFN2::calculateAES2DipoleDipoleMatrix(const Vector& rad, double kdmp5) const
{
    // Tensor index representation: [i][j](A, B) where i,j are x,y,z
    std::vector<Matrix> M(9, Matrix::Zero(m_atomcount, m_atomcount));
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            Eigen::Vector3d RAB = (m_geometry.row(B) - m_geometry.row(A)).transpose() / au;
            double r2 = RAB.squaredNorm();
            double r = std::sqrt(r2);
            double r3 = r2 * r;
            double r5 = r3 * r2;
            double rr = r / (rad(A) + rad(B));
            double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr, kdmp5));

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    double T_ij = (i == j ? 1.0 / r3 : 0.0) - 3.0 * RAB[i] * RAB[j] / r5;
                    M[i*3 + j](A, B) = fdmp5 * T_ij;
                }
            }
        }
    }
    return M;
}

Matrix GFN2::calculateAES2ChargeQuadrupoleMatrix(const Vector& rad, double kdmp5) const
{
    // Returns matrix (6 * natoms x natoms)
    Matrix M = Matrix::Zero(6 * m_atomcount, m_atomcount);
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            Eigen::Vector3d RAB = (m_geometry.row(B) - m_geometry.row(A)).transpose() / au;
            double r2 = RAB.squaredNorm();
            double r = std::sqrt(r2);
            double r5 = r2 * r2 * r;
            double rr = r / (rad(A) + rad(B));
            double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr, kdmp5));

            // Order: xx, xy, xz, yy, yz, zz
            M(0*m_atomcount + A, B) = fdmp5 * RAB[0] * RAB[0] / r5;
            M(1*m_atomcount + A, B) = fdmp5 * 2.0 * RAB[0] * RAB[1] / r5;
            M(2*m_atomcount + A, B) = fdmp5 * 2.0 * RAB[0] * RAB[2] / r5;
            M(3*m_atomcount + A, B) = fdmp5 * RAB[1] * RAB[1] / r5;
            M(4*m_atomcount + A, B) = fdmp5 * 2.0 * RAB[1] * RAB[2] / r5;
            M(5*m_atomcount + A, B) = fdmp5 * RAB[2] * RAB[2] / r5;
        }
    }
    return M;
}

Matrix GFN2::calculateAES2ChargeDipoleMatrix(const Vector& rad, double kdmp3) const
{
    Matrix M = Matrix::Zero(3 * m_atomcount, m_atomcount);
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            Eigen::Vector3d RAB = (m_geometry.row(B) - m_geometry.row(A)).transpose() / au;
            double r2 = RAB.squaredNorm();
            double r = std::sqrt(r2);
            double r3 = r2 * r;
            double rr = r / (rad(A) + rad(B));
            double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr, kdmp3));

            for (int i = 0; i < 3; ++i) {
                M(i * m_atomcount + A, B) = fdmp3 * RAB[i] / r3;
            }
        }
    }
    return M;
}

double GFN2::calculateAES2Energy(const std::vector<Eigen::Vector3d>& dipoles,
                               const std::vector<Eigen::VectorXd>& quadrupoles) const
{
    double e_aes2 = 0.0;
    Vector q = m_charges;
    Vector rad(m_atomcount);
    for (int i = 0; i < m_atomcount; ++i) {
        if (m_param_db.hasElement(m_atoms[i]))
            rad(i) = m_param_db.getElement(m_atoms[i]).multipole_rad;
        else
            rad(i) = 2.0;
    }

    // 1. Kernel Energy (On-site contribution)
    for (int A = 0; A < m_atomcount; ++A) {
        if (!m_param_db.hasElement(m_atoms[A])) continue;
        const auto& epA = m_param_db.getElement(m_atoms[A]);
        e_aes2 += epA.dkernel * dipoles[A].squaredNorm();
    }

    // 2. Charge-Dipole interaction (damping = 3.0, TBLite gfn2.f90 line 470)
    Matrix amat_sd = calculateAES2ChargeDipoleMatrix(rad, GFN2Params::GFN2_AES_DMP3);
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            Eigen::Vector3d grad_gamma_AB;
            for (int i = 0; i < 3; ++i) grad_gamma_AB[i] = amat_sd(i * m_atomcount + A, B);
            e_aes2 += q(A) * dipoles[B].dot(grad_gamma_AB);
        }
    }

    // 3. Dipole-Dipole interaction (damping = 4.0, TBLite gfn2.f90 line 470)
    std::vector<Matrix> amat_dd = calculateAES2DipoleDipoleMatrix(rad, GFN2Params::GFN2_AES_DMP5);
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            double v_dd = 0.0;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    v_dd += dipoles[A][i] * amat_dd[i * 3 + j](A, B) * dipoles[B][j];
                }
            }
            e_aes2 += 0.5 * v_dd;
        }
    }

    // 4. Charge-Quadrupole interaction (damping = 4.0, TBLite gfn2.f90 line 470)
    Matrix amat_sq = calculateAES2ChargeQuadrupoleMatrix(rad, GFN2Params::GFN2_AES_DMP5);
    for (int A = 0; A < m_atomcount; ++A) {
        for (int B = 0; B < m_atomcount; ++B) {
            if (A == B) continue;
            double v_sq = 0.0;
            for (int k = 0; k < 6; ++k) {
                v_sq += q(A) * amat_sq(k * m_atomcount + A, B) * quadrupoles[B][k];
            }
            e_aes2 += v_sq; // 0.5 factor already in matrix
        }
    }

    return e_aes2;
}
