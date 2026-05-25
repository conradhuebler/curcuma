/*
 * <Unified Native xTB Implementation — GFN1 / GFN2>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Skeleton implementation of curcuma::xtb::XTB.  The class is structured so
 * each SCC-relevant kernel (H0, isotropic Coulomb, third-order, multipole,
 * SCF iterator) lives in its own translation unit mirroring tblite's module
 * layout.  This file holds the public API, state initialisation, and glue;
 * the kernel files (xtb_h0.cpp, xtb_coulomb.cpp, …) are added as Phase 3
 * progresses.
 *
 * Phase 2 status (Apr 2026): skeleton compiles cleanly; Calculation()
 * returns 0 and warns that the native xTB SCC is not yet implemented.  This
 * is intentional — the class is wired into the build so subsequent kernel
 * work links against a real type, but not yet into MethodFactory.
 *
 * Claude Generated: Phase 2 scaffolding
 *
 * This program is free software under GPL-3.0.
 */

#include "xtb_native.h"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include "STO_CGTO.hpp"
#include "src/core/curcuma_logger.h"
#include "src/core/config_manager.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "diis_accelerator.h"

#include <chrono>
#include <stdexcept>

namespace curcuma::xtb {

/* ------------------------------------------------------------------------- *
 *  Lifecycle
 * ------------------------------------------------------------------------- */
XTB::XTB(MethodType method)
    : QMDriver()
    , m_method(method)
{
}

XTB::~XTB() = default;

bool XTB::InitialiseMolecule()
{
    if (m_atomcount <= 0) {
        CurcumaLogger::error("XTB::InitialiseMolecule: no atoms set");
        return false;
    }
    buildBasis();
    buildH0Data();
    buildReferenceOccupations();
    return true;
}

// Non-member helper for convergence check (must be defined before Calculation)
namespace {
    bool checkConvergence_impl(const Vector& q_old, const Vector& q_new,
                                double e_old, double e_new, double thresh)
    {
        const double dq = (q_new - q_old).cwiseAbs().maxCoeff();
        const double de = std::fabs(e_new - e_old);
        return (dq < thresh && de < thresh * 100.0);
    }
}

double XTB::Calculation(bool gradient)
{
    using clock = std::chrono::steady_clock;
    auto ms = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    const int verb = CurcumaLogger::get_verbosity();
    const auto t0 = clock::now();

    if (verb >= 2) {
        CurcumaLogger::header(fmt::format("{} SCF", methodName()));
        CurcumaLogger::param("atoms",    m_atomcount);
        CurcumaLogger::param("shells",   m_basis.nsh);
        CurcumaLogger::param("basis_fns", m_basis.nao);
        CurcumaLogger::param("electrons", static_cast<int>(m_wfn.nocc));
    }

    // 1. Coordination numbers
    Vector cn = computeCoordinationNumbers();
    const auto t_cn = clock::now();

    // 2. Self-energies
    Vector se;
    getSelfEnergies(cn, se);

    // 3. Build overlap and H0
    Matrix S, H0;
    getHamiltonianH0(se, S, H0);
    m_S  = S;
    m_H0 = H0;
    const auto t_h0 = clock::now();

    // 4. Build gamma matrix (once per geometry)
    buildGammaMatrix();
    const auto t_gamma = clock::now();

    // 5. Multipole setup (GFN2 only) — fills m_dp_int, m_qp_int, interaction matrices
    if (m_method == MethodType::GFN2) {
        setupMultipole();
    }
    const auto t_setup = clock::now();

    if (verb >= 3) {
        CurcumaLogger::info("Setup timing:");
        CurcumaLogger::info_fmt("  coordination numbers : {:8.2f} ms", ms(t0, t_cn));
        CurcumaLogger::info_fmt("  overlap + H0         : {:8.2f} ms", ms(t_cn, t_h0));
        CurcumaLogger::info_fmt("  Coulomb gamma matrix : {:8.2f} ms", ms(t_h0, t_gamma));
        if (m_method == MethodType::GFN2)
            CurcumaLogger::info_fmt("  multipole setup      : {:8.2f} ms", ms(t_gamma, t_setup));
    }

    // 6. SCF loop with DIIS acceleration (Pulay 1980/1982)
    m_pot.reset();

    if (verb >= 2) {
        CurcumaLogger::info("SCF iterations (DIIS):");
        CurcumaLogger::info("  iter             E(SCC)/Eh            dE        max|dq|     t/ms");
    }

    const int max_iter   = m_scf_max_iter;
    const double thresh  = m_scf_threshold;
    const int nsh = m_basis.nsh;

    // Initial guess: zero potential → F = H0
    Vector q_sh_old = Vector::Zero(nsh);
    Vector q_sh_new;
    double e_total_old = 0.0;

    // DIIS accelerator: keep last 6 Fock matrices
    // Push from iter >= 1 (P is valid only after first solveEigen call)
    DIISAccelerator diis(6);

    const auto t_scf_start = clock::now();
    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        const auto t_iter0 = clock::now();
        // Reset and build potentials
        m_pot.reset();

        // Isotropic Coulomb: v_sh += γ * q_sh
        addCoulombShellPotential(m_pot);

        // Third-order: add to v_at/v_sh
        addThirdOrderPotential(m_pot);

        // GFN2 multipole
        if (m_method == MethodType::GFN2) {
            addMultipolePotential(m_pot);
        }

        // Expand to AO potential and build Fock
        Matrix F = buildFock(m_H0, m_S, m_pot);

        // DIIS: push F+P into history and extrapolate when ≥2 vectors available.
        // P is only valid from iter >= 1 (after the first solveEigen call).
        if (iter >= 1) {
            diis.push(F, m_wfn.P, m_S);
            if (diis.size() >= 2) {
                F = diis.extrapolate();
            }
        }

        // Diagonalize
        if (!solveEigen(F, m_S)) {
            CurcumaLogger::warn("XTB::Calculation: eigen solve failed at iteration " + std::to_string(iter));
            m_scf_converged = false;
            m_scf_iterations = iter;
            return m_E_total;
        }

        // Update populations
        updatePopulations(m_S);
        q_sh_new = m_wfn.q_sh;

        // Energies for this iteration
        m_E_coulomb_shell = energyCoulombShell();
        m_E_third_order   = energyThirdOrder();
        m_E_multipole     = energyMultipole();

        // Band energy: Tr(P · H0)
        m_E_electronic = (m_wfn.P.cwiseProduct(m_H0)).sum();

        // Total SCC = band + coulomb + third-order + multipole
        const double e_scc = m_E_electronic + m_E_coulomb_shell
                           + m_E_third_order + m_E_multipole;

        // Per-iteration diagnostics
        const double dq = (q_sh_new - q_sh_old).cwiseAbs().maxCoeff();
        const double de = (iter > 0) ? std::fabs(e_scc - e_total_old) : 0.0;
        const double t_iter_ms = ms(t_iter0, clock::now());
        if (verb >= 2) {
            std::string line = fmt::format("  {:4d}   {:18.10f}   {:10.2e}   {:10.2e}   {:6.1f}",
                                           iter, e_scc, de, dq, t_iter_ms);
            if (verb >= 3)
                line += fmt::format("   [DIIS n={} err={:.2e}]", diis.size(), diis.lastErrorNorm());
            CurcumaLogger::info(line);
        }

        // Check convergence (after first iteration)
        if (iter > 0) {
            if (checkConvergence_impl(q_sh_old, q_sh_new,
                                       e_total_old, e_scc, thresh)) {
                m_scf_converged = true;
                break;
            }
        }

        q_sh_old = m_wfn.q_sh;
        e_total_old = e_scc;
    }

    m_scf_iterations = iter + 1;
    const auto t_scf_end = clock::now();

    if (verb >= 1) {
        if (m_scf_converged)
            CurcumaLogger::success_fmt("SCF converged in {} iterations ({:.1f} ms)",
                                       m_scf_iterations, ms(t_scf_start, t_scf_end));
        else
            CurcumaLogger::warn_fmt("SCF NOT converged after {} iterations ({:.1f} ms)",
                                    m_scf_iterations, ms(t_scf_start, t_scf_end));
    }

    // Persist converged Fock matrix for gradient / debug
    m_F = buildFock(m_H0, m_S, m_pot);

    // Final energies
    m_E_electronic    = energyCoulombShell() + energyThirdOrder() + energyMultipole();
    // Band energy: Tr(P · H0) for electronic part
    m_E_electronic   += (m_wfn.P.cwiseProduct(m_H0)).sum();
    m_E_repulsion     = calcRepulsionEnergy();
    m_E_halogen_bond  = calcHalogenBondEnergy();
    m_E_dispersion    = calcDispersionEnergy();

    m_E_total = m_E_electronic + m_E_repulsion
              + m_E_halogen_bond + m_E_dispersion;
    const auto t_energies = clock::now();

    if (verb >= 2) {
        CurcumaLogger::info("Energy decomposition:");
        CurcumaLogger::info(fmt::format("  Electronic   = {: .8f} Eh", m_E_electronic));
        CurcumaLogger::info(fmt::format("  Coulomb ES2  = {: .8f} Eh", m_E_coulomb_shell));
        CurcumaLogger::info(fmt::format("  Third-order  = {: .8f} Eh", m_E_third_order));
        if (m_method == MethodType::GFN2)
            CurcumaLogger::info(fmt::format("  Multipole    = {: .8f} Eh", m_E_multipole));
        CurcumaLogger::info(fmt::format("  Repulsion    = {: .8f} Eh", m_E_repulsion));
        if (m_method == MethodType::GFN1)
            CurcumaLogger::info(fmt::format("  Halogen bond = {: .8f} Eh", m_E_halogen_bond));
        CurcumaLogger::info(fmt::format("  Dispersion   = {: .8f} Eh", m_E_dispersion));
        CurcumaLogger::result(fmt::format("  Total        = {: .8f} Eh", m_E_total));

        // Frontier orbitals
        const double homo = getHOMOEnergy();
        const double lumo = getLUMOEnergy();
        CurcumaLogger::info(fmt::format("HOMO = {: .6f} Eh   LUMO = {: .6f} Eh   gap = {:.4f} eV",
                                        homo, lumo, (lumo - homo) * 27.211386245988));
    }

    // Update QMDriver state for wrapper compatibility
    m_mo = m_wfn.C;
    m_energies = m_wfn.eps;
    m_num_electrons = static_cast<int>(m_wfn.nocc);
    m_coordination_numbers = cn;

    const auto t_grad0 = clock::now();
    if (gradient) {
        // Rebuild potential with converged shell charges before computing gradient.
        // (m_pot was overwritten during SCF; final values must be regenerated.)
        m_pot.reset();
        addCoulombShellPotential(m_pot);
        addThirdOrderPotential(m_pot);
        if (m_method == MethodType::GFN2) addMultipolePotential(m_pot);

        calculateGradient();       // fills m_gradient in Eh/Bohr
        m_gradient /= au;          // Eh/Bohr → Eh/Å: 1 Eh/Bohr = (1/au) Eh/Å
    }
    const auto t_end = clock::now();

    if (verb >= 2) {
        CurcumaLogger::info("Timing breakdown:");
        CurcumaLogger::info_fmt("  setup        : {:8.2f} ms", ms(t0, t_setup));
        CurcumaLogger::info_fmt("  SCF ({:3d} it) : {:8.2f} ms", m_scf_iterations, ms(t_scf_start, t_scf_end));
        CurcumaLogger::info_fmt("  post-SCF E   : {:8.2f} ms", ms(t_scf_end, t_energies));
        if (gradient)
            CurcumaLogger::info_fmt("  gradient     : {:8.2f} ms", ms(t_grad0, t_end));
        CurcumaLogger::info_fmt("  TOTAL        : {:8.2f} ms", ms(t0, t_end));
    }

    return m_E_total;
}

/* ------------------------------------------------------------------------- *
 *  Build-once state (stub implementations)
 * ------------------------------------------------------------------------- */
void XTB::buildBasis()
{
    m_basis = BasisMap{};
    m_basis.nat = m_atomcount;
    m_basis.z.assign(m_atoms.begin(), m_atoms.end());

    int shell_cursor = 0;
    int ao_cursor    = 0;
    m_basis.ish_at.resize(m_atomcount, 0);
    m_basis.nsh_at.resize(m_atomcount, 0);

    auto nshellFor = [this](int z) -> int {
        const int idx = z - 1;
        if (idx < 0 || idx >= gfn1_params::MAX_ELEM) return 0;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::nshell[idx]
                 : gfn2_params::nshell[idx];
    };
    auto angFor = [this](int z, int ish) -> int {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::ang_shell[idx][ish]
                 : gfn2_params::ang_shell[idx][ish];
    };
    auto principalFor = [this](int z, int ish) -> int {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::principal_quantum_number[idx][ish]
                 : gfn2_params::principal_quantum_number[idx][ish];
    };
    auto zetaFor = [this](int z, int ish) -> double {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::slater_exponent[idx][ish]
                 : gfn2_params::slater_exponent[idx][ish];
    };
    auto nprimFor = [this](int z, int ish) -> int {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::number_of_primitives[idx][ish]
                 : gfn2_params::number_of_primitives[idx][ish];
    };

    for (int iat = 0; iat < m_atomcount; ++iat) {
        const int z     = m_atoms[iat];
        const int nshell = nshellFor(z);
        m_basis.ish_at[iat] = shell_cursor;
        m_basis.nsh_at[iat] = nshell;

        // Build CGTO::Shell for each shell on this atom
        std::vector<CGTO::Shell> cgs(nshell);
        for (int ish = 0; ish < nshell; ++ish) {
            cgs[ish] = CGTO::slater_to_gauss(
                principalFor(z, ish),
                angFor(z, ish),
                zetaFor(z, ish),
                nprimFor(z, ish));
        }
        // Gram-Schmidt orthogonalize same-l shells
        for (int ish = 1; ish < nshell; ++ish) {
            for (int jsh = 0; jsh < ish; ++jsh) {
                if (cgs[ish].ang == cgs[jsh].ang) {
                    CGTO::orthogonalize(cgs[jsh], cgs[ish]);
                }
            }
        }

        for (int ish = 0; ish < nshell; ++ish) {
            const int ang = cgs[ish].ang;
            const int nao = 2 * ang + 1;

            CGTOShell cg{};
            cg.ang         = ang;
            cg.principal_n = principalFor(z, ish);
            cg.slater_exp  = zetaFor(z, ish);
            cg.n_prim      = static_cast<int>(cgs[ish].alpha.size());
            cg.alpha       = std::move(cgs[ish].alpha);
            cg.coeff       = std::move(cgs[ish].coeff);
            m_basis.cgto.push_back(std::move(cg));

            m_basis.ang_sh.push_back(ang);
            m_basis.nao_sh.push_back(nao);
            m_basis.iao_sh.push_back(ao_cursor);
            m_basis.sh2at.push_back(iat);
            for (int iao = 0; iao < nao; ++iao) {
                m_basis.ao2at.push_back(iat);
                m_basis.ao2sh.push_back(shell_cursor);
            }
            ao_cursor    += nao;
            shell_cursor += 1;
        }
    }
    m_basis.nsh = shell_cursor;
    m_basis.nao = ao_cursor;
}

void XTB::buildH0Data()
{
    m_h0 = H0Data{};
    m_h0.selfenergy.assign(m_basis.nsh, 0.0);
    m_h0.kcn       .assign(m_basis.nsh, 0.0);
    m_h0.shpoly    .assign(m_basis.nsh, 0.0);
    m_h0.rad       .assign(m_basis.nat, 0.0);
    for (int iat = 0; iat < m_basis.nat; ++iat) {
        m_h0.rad[iat] = atomic_rad_au(m_basis.z[iat]);
    }

    for (int iat = 0; iat < m_basis.nat; ++iat) {
        const int z = m_basis.z[iat];
        for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
            const int sh_idx = m_basis.ish_at[iat] + ish;
            if (m_method == MethodType::GFN1) {
                m_h0.selfenergy[sh_idx] = gfn1_params::p_selfenergy[z - 1][ish];
                m_h0.kcn       [sh_idx] = gfn1_params::p_kcn       [z - 1][ish];
                m_h0.shpoly    [sh_idx] = gfn1_params::p_shpoly    [z - 1][ish];
            } else {
                m_h0.selfenergy[sh_idx] = gfn2_params::p_selfenergy[z - 1][ish];
                m_h0.kcn       [sh_idx] = gfn2_params::p_kcn       [z - 1][ish];
                m_h0.shpoly    [sh_idx] = gfn2_params::p_shpoly    [z - 1][ish];
            }
        }
    }
    // hscale is computed on-the-fly in getHamiltonianH0() (xtb_h0.cpp).
    // m_h0.hscale is kept as zero and not used; the per-shell-pair logic
    // (Slater-ratio prefactor, kpair, kshell, EN scaling) lives in the H0 builder.
    m_h0.hscale = Matrix::Zero(m_basis.nsh, m_basis.nsh);
}

void XTB::buildReferenceOccupations()
{
    m_wfn = Wavefunction{};
    m_wfn.n0_at = Vector::Zero(m_basis.nat);
    m_wfn.n0_sh = Vector::Zero(m_basis.nsh);

    double total = -static_cast<double>(m_charge);
    for (int iat = 0; iat < m_basis.nat; ++iat) {
        const int z = m_basis.z[iat];
        for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
            const int sh_idx = m_basis.ish_at[iat] + ish;
            const double refocc = (m_method == MethodType::GFN1)
                                    ? gfn1_params::reference_occ[z - 1][ish]
                                    : gfn2_params::reference_occ[z - 1][ish];
            m_wfn.n0_sh(sh_idx) += refocc;
            m_wfn.n0_at(iat)    += refocc;
            total               += refocc;
        }
    }
    m_wfn.nocc = total;
    m_wfn.q_at.setZero(m_basis.nat);
    m_wfn.q_sh.setZero(m_basis.nsh);
    if (m_method == MethodType::GFN2) {
        m_wfn.dp_at = Eigen::MatrixXd::Zero(3, m_basis.nat);
        m_wfn.qp_at = Eigen::MatrixXd::Zero(6, m_basis.nat);
    }

    m_pot.v_at = Eigen::VectorXd::Zero(m_basis.nat);
    m_pot.v_sh = Eigen::VectorXd::Zero(m_basis.nsh);
    m_pot.v_ao = Eigen::VectorXd::Zero(m_basis.nao);
    m_pot.v_dp = Eigen::MatrixXd::Zero(3, m_basis.nat);
    m_pot.v_qp = Eigen::MatrixXd::Zero(6, m_basis.nat);
}

/* ------------------------------------------------------------------------- *
 *  Orbital-energy accessors
 * ------------------------------------------------------------------------- */
double XTB::getHOMOEnergy() const
{
    const int homo_idx = static_cast<int>(std::floor(m_wfn.nocc / 2.0)) - 1;
    if (homo_idx < 0 || homo_idx >= m_wfn.eps.size()) return 0.0;
    return m_wfn.eps(homo_idx);
}
double XTB::getLUMOEnergy() const
{
    const int lumo_idx = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
    if (lumo_idx < 0 || lumo_idx >= m_wfn.eps.size()) return 0.0;
    return m_wfn.eps(lumo_idx);
}
double XTB::getHOMOLUMOGap() const { return getLUMOEnergy() - getHOMOEnergy(); }

/* ------------------------------------------------------------------------- *
 *  Geometry update with full cache invalidation
 * ------------------------------------------------------------------------- */
bool XTB::UpdateMolecule(const Matrix& geometry)
{
    // 1. Update base class geometry
    m_geometry = geometry;

    // 2. Invalidate all geometry-dependent cached state
    //    Overlap and H0 depend on atom positions via STO-CGTO overlap integrals
    m_S.resize(0, 0);
    m_H0.resize(0, 0);

    // 3. Fock matrix depends on geometry via H0 + potential
    m_F.resize(0, 0);

    // 4. Gamma matrix depends on interatomic distances (Klopman-Ohno R_AB)
    m_gamma.resize(0, 0);

    // 5. Multipole interaction matrices depend on distances and damping radii
    m_mp_initialized = false;
    for (auto& m : m_mp_amat_sd) m.resize(0, 0);
    for (auto& row : m_mp_amat_dd)
        for (auto& m : row) m.resize(0, 0);
    for (auto& m : m_mp_amat_sq) m.resize(0, 0);
    m_mp_mrad.clear();
    m_mp_dkernel.clear();
    m_mp_qkernel.clear();

    // 5. Reset wavefunction populations (Mulliken charges depend on geometry)
    //    Keep m_wfn.n0_at / n0_sh (reference occupations, geometry-independent)
    m_wfn.q_at.setZero(m_basis.nat);
    m_wfn.q_sh.setZero(m_basis.nsh);
    m_wfn.dp_at.setZero(3, m_basis.nat);
    m_wfn.qp_at.setZero(6, m_basis.nat);

    // 6. Energy components are now stale
    m_E_electronic = m_E_repulsion = m_E_coulomb_shell = 0.0;
    m_E_third_order = m_E_multipole = m_E_halogen_bond = m_E_dispersion = 0.0;
    m_E_total = 0.0;
    m_scf_converged = false;
    m_scf_iterations = 0;

    return true;
}

/* ------------------------------------------------------------------------- *
 *  Energy decomposition accessor
 * ------------------------------------------------------------------------- */
nlohmann::json XTB::getEnergyDecomposition() const
{
    nlohmann::json j;
    j["electronic"]     = m_E_electronic;
    j["coulomb_shell"]  = m_E_coulomb_shell;
    j["third_order"]    = m_E_third_order;
    j["multipole"]      = m_E_multipole;
    j["repulsion"]      = m_E_repulsion;
    j["halogen_bond"]   = m_E_halogen_bond;
    j["dispersion"]     = m_E_dispersion;
    j["total"]          = m_E_total;
    j["scf_converged"]  = m_scf_converged;
    j["scf_iterations"] = m_scf_iterations;

    // AP6a: dump converged Fock matrix (was H0 before; F = H0 + potential)
    if (m_F.rows() > 0 && m_F.cols() > 0) {
        nlohmann::json fmat = nlohmann::json::array();
        for (int i = 0; i < m_F.rows(); ++i) {
            nlohmann::json row = nlohmann::json::array();
            for (int j = 0; j < m_F.cols(); ++j)
                row.push_back(m_F(i, j));
            fmat.push_back(row);
        }
        j["debug_fock"] = fmat;
    }

    // AP6 debug: dump GFN2 multipole quantities
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        nlohmann::json vat = nlohmann::json::array();
        for (int i = 0; i < m_atomcount; ++i) vat.push_back(m_pot.v_at(i));
        j["debug_vat_extra"] = vat;

        nlohmann::json vdp = nlohmann::json::array();
        for (int i = 0; i < m_atomcount; ++i)
            vdp.push_back({m_pot.v_dp(0,i), m_pot.v_dp(1,i), m_pot.v_dp(2,i)});
        j["debug_vdp"] = vdp;

        nlohmann::json vqp = nlohmann::json::array();
        for (int i = 0; i < m_atomcount; ++i)
            vqp.push_back({m_pot.v_qp(0,i), m_pot.v_qp(1,i), m_pot.v_qp(2,i),
                           m_pot.v_qp(3,i), m_pot.v_qp(4,i), m_pot.v_qp(5,i)});
        j["debug_vqp"] = vqp;

        nlohmann::json dpat = nlohmann::json::array();
        for (int i = 0; i < m_atomcount; ++i)
            dpat.push_back({m_wfn.dp_at(0,i), m_wfn.dp_at(1,i), m_wfn.dp_at(2,i)});
        j["debug_dpat"] = dpat;

        nlohmann::json qpat = nlohmann::json::array();
        for (int i = 0; i < m_atomcount; ++i)
            qpat.push_back({m_wfn.qp_at(0,i), m_wfn.qp_at(1,i), m_wfn.qp_at(2,i),
                            m_wfn.qp_at(3,i), m_wfn.qp_at(4,i), m_wfn.qp_at(5,i)});
        j["debug_qpat"] = qpat;
    }
    return j;
}

/* ------------------------------------------------------------------------- *
 *  GFN2 D4 dispersion — native, via curcuma::dispersion::D4Evaluator
 *
 *  Standard BJ damping (Caldeweyher 2019), GFN2 parameters from TBLite
 *  (gfn2.f90 / Ulysses D4par.hpp): s6=1.0, s8=2.7, a1=0.52, a2=5.0,
 *  s9=5.0, alpha=16.0.
 *
 *  The same D4Evaluator class also powers GFN-FF's D4 (with modified-BJ
 *  parameters s6=1, s8=2). See dispersion/CLAUDE.md for the unified
 *  formula derivation.
 *
 *  Caveat: zeta-charge scaling (zetac6) is currently treated as a static
 *  prefactor as in GFN-FF. The full GFN2 q-response chain rule
 *  (∂E_D4/∂q · ∂q/∂x via SCF response) is deferred — expected residual
 *  is sub-mEh on the AP test set (H2O, CH4, CH3OCH3, caffeine).
 *
 *  GFN1 still uses D3 (separate AP) — returns 0 here for now, preserving
 *  the existing GFN1 behaviour at the cost of a small dispersion deficit.
 * ------------------------------------------------------------------------- */
double XTB::calcDispersionEnergy() const
{
    m_disp_gradient_valid = false;

    if (m_method != MethodType::GFN2) {
        // GFN1 D3 integration is a follow-up AP. Returning 0 preserves the
        // current (D3-missing) energy for GFN1 — see qm_methods/CLAUDE.md.
        return 0.0;
    }

    // Lazy initialization (cached across SCF cycles for the same geometry)
    if (!m_d4_generator) {
        nlohmann::json d4_config;
        d4_config["d4_s6"] = 1.0;
        d4_config["d4_s8"] = 2.7;
        d4_config["d4_a1"] = 0.52;
        d4_config["d4_a2"] = 5.0;
        d4_config["d4_alp"] = 16.0;
        ConfigManager cfg("d4param", d4_config);
        m_d4_generator = std::make_unique<::D4ParameterGenerator>(cfg);
    }
    if (!m_d4_evaluator) {
        curcuma::dispersion::D4Params p;
        p.s6 = 1.0;
        p.s8 = 2.7;
        p.a1 = 0.52;
        p.a2 = 5.0;
        p.s9 = 5.0;
        p.alpha = 16.0;
        p.damping = curcuma::dispersion::DampingFormula::StandardBJ_D4;
        m_d4_evaluator = std::make_unique<curcuma::dispersion::D4Evaluator>(
            m_d4_generator.get(), p);
    }

    // Charge source for the q-response:
    //  - "eeq": canonical single-shot dftd4 EEQ (analytical dq/dx, D4ChargeModel).
    //  - "mulliken": the converged GFN2 SCF atomic charges drive zetac6; the
    //    dq/dx response comes from CPSCF (computeMullikenChargeResponse).
    const bool use_mulliken = (m_d4_charge_source == "mulliken");
    m_d4_generator->setUseD4SingleShotEEQ(m_d4_charge_source == "eeq");
    if (use_mulliken) {
        // Feed the SCF Mulliken charges as the zeta charges. setTopologyCharges
        // makes getZetaCharges() return them and GenerateDispersionPairsNative
        // build zetac6 from them; dEdq (Phase 1) is then taken w.r.t. these.
        m_d4_generator->setTopologyCharges(m_wfn.q_at);
    }

    // m_geometry is in Angstrom; D4 expects Bohr.
    Matrix geom_bohr = m_geometry * AA_TO_AU;
    m_d4_generator->GenerateParameters(m_atoms, geom_bohr);

    // Always compute the gradient block — D4 is cheap and the cached
    // m_disp_gradient / m_disp_dEdcn are folded into calculateGradient()
    // when a gradient was requested at the top of Calculation().
    //
    // m_disp_dEdq holds dE_D4/dq (charge-response chain rule, first half). It
    // is populated whenever a charge source is selected; the dq/dx contraction
    // that turns it into a Cartesian contribution is applied below (EEQ) or in
    // the GFN2 CPSCF path (mulliken). With the default "eeq" source this is
    // the dftd4-conform behaviour.
    const bool want_dEdq = !m_d4_charge_source.empty();
    const double E = m_d4_evaluator->computeEnergyAndGradient(
        m_atoms, geom_bohr,
        /*with_gradient=*/true,
        m_disp_gradient, m_disp_dEdcn, m_disp_dEdq, want_dEdq);

    // q-response chain rule: fold Σ_A dE_D4/dq_A · ∂q_A/∂R into the cached
    // geometry gradient. For the "eeq" source this is the analytical single-shot
    // EEQ response (D4ChargeModel); xtb_gradient.cpp then adds m_disp_gradient
    // to the total gradient unchanged.
    if (want_dEdq && m_d4_charge_source == "eeq" && m_disp_dEdq.size() == m_atomcount) {
        m_d4_generator->addChargeResponseGradient(m_disp_dEdq, m_disp_gradient);
    } else if (want_dEdq && m_d4_charge_source == "mulliken"
               && m_disp_dEdq.size() == m_atomcount) {
        // GFN2 CPSCF/Z-vector response ∂q_Mulliken/∂x (Phase 3b, isotropic).
        computeMullikenChargeResponse(m_disp_dEdq, m_disp_gradient);
    }
    m_disp_gradient_valid = true;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dispersion_d4_native", fmt::format("{:.6f} Eh", E));
    }
    return E;
}

/* ------------------------------------------------------------------------- *
 *  Legacy QMDriver hooks — not used yet; we go through buildH0Data() instead.
 * ------------------------------------------------------------------------- */
Matrix XTB::MakeOverlap(std::vector<STO::Orbital>& /*basisset*/)
{
    return Matrix::Identity(m_basis.nao, m_basis.nao);
}
Matrix XTB::MakeH(const Matrix& /*S*/, const std::vector<STO::Orbital>& /*basisset*/)
{
    return Matrix::Zero(m_basis.nao, m_basis.nao);
}

} // namespace curcuma::xtb
