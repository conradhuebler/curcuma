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
#include "src/core/energy_calculators/dispersion/d4_ncoord.h"  // D4 CN gradient (GFN2)
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "diis_accelerator.h"
#include "broyden_mixer.h"

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

/* ------------------------------------------------------------------------- *
 *  EEQ initial guess (Claude Generated).
 *
 *  Seeds shell charges from a single-shot dftd4 EEQ solve (one smooth linear
 *  system, curcuma::dispersion::D4ChargeModel). The per-atom EEQ charge is
 *  partitioned across that atom's shells in proportion to the reference
 *  occupations n0_sh, so Σ_{s∈A} q_sh(s) = q_at(A). Building the iter-0 Fock
 *  from these charges starts the SCF with a physically correct Coulomb shift
 *  and keeps polar systems out of the wrong charge-transfer basin.
 * ------------------------------------------------------------------------- */
bool XTB::seedEEQGuess(Vector& q_sh_out)
{
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;
    if (nat <= 0 || nsh <= 0)
        return false;

    // m_geometry is in Angstrom; the EEQ model expects Bohr.
    const Matrix geom_bohr = m_geometry * AA_TO_AU;
    curcuma::dispersion::D4ChargeModel eeq;
    Vector q_at = eeq.computeCharges(m_atoms, geom_bohr,
                                     static_cast<double>(m_charge));
    if (q_at.size() != nat || !q_at.allFinite())
        return false;

    q_sh_out = Vector::Zero(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int iat = m_basis.sh2at[s];
        const double n0a = m_wfn.n0_at(iat);
        const double w = (n0a > 1.0e-12)
                             ? (m_wfn.n0_sh(s) / n0a)
                             : (1.0 / static_cast<double>(
                                   m_basis.nsh_at[iat] > 0 ? m_basis.nsh_at[iat] : 1));
        q_sh_out(s) = q_at(iat) * w;
    }
    return true;
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
        CurcumaLogger::info(fmt::format("SCF iterations (mode={}, guess={}):",
                                        scfModeName(m_scf_mode), m_scf_guess));
        CurcumaLogger::info("  iter             E(SCC)/Eh            dE        max|dq|     t/ms");
    }

    const int max_iter   = m_scf_max_iter;
    const double thresh  = m_scf_threshold;
    const int nsh = m_basis.nsh;

    // Initial guess. Default ("h0"): zero shell charges → F = H0. The optional
    // "eeq" guess seeds q_sh from a single-shot dftd4 EEQ solve so iter 0 builds
    // the Fock with a physically correct Coulomb shift, avoiding the wrong
    // charge-transfer basin instead of having to recover from it.
    Vector q_sh_old = Vector::Zero(nsh);
    if (m_scf_guess == "eeq") {
        Vector q_sh_guess;
        if (seedEEQGuess(q_sh_guess)) {
            q_sh_old   = q_sh_guess;
            m_wfn.q_sh = q_sh_guess;
            m_wfn.q_at.setZero(m_atomcount);
            for (int s = 0; s < nsh; ++s)
                m_wfn.q_at(m_basis.sh2at[s]) += q_sh_guess(s);
            if (verb >= 2)
                CurcumaLogger::info("SCF initial guess: single-shot EEQ shell charges");
        } else if (verb >= 1) {
            CurcumaLogger::warn("EEQ initial guess unavailable; falling back to bare-H0 guess");
        }
    }
    Vector q_sh_new;
    double e_total_old = 0.0;

    // SCF strategy (selectable). The default DIIS path with diis_start=5 and a
    // 6-Fock subspace reproduces the historic charge-sloshing control: a few
    // strongly damped warmup iterations keep polar / multiple-bond systems in
    // the ground-state basin before Pulay DIIS takes over.
    const double  damp       = m_scf_damping;   // density mixing factor
    const int     diis_start = m_diis_start;    // damped warmup iterations before DIIS
    const ScfMode mode       = m_scf_mode;
    Matrix P_old;
    double dq_prev = 1.0e30;   // last max|dq| — controls the level-shift fade-out

    // DIIS accelerator (history depth configurable). Unused in Plain mode.
    DIISAccelerator diis(m_diis_subspace);
    // Never extrapolate from a pathologically large commutator error — that is
    // exactly what threw `complex` out of the basin at iter 5. The cutoff is
    // generous; a well-behaved SCF never approaches it, so the default DIIS path
    // is bit-for-bit unchanged on the existing test set.
    constexpr double kDiisErrorCutoff = 1.0e3;

    // Broyden mixer on the SCC charge/multipole vector (tblite-style). Unused in
    // the other modes. alpha = m_scf_damping so -scf_damping tunes the seed step.
    BroydenMixer broyden(damp, /*max_history=*/m_diis_subspace > 2 ? 20 : m_diis_subspace);
    const int nat = m_atomcount;

    // Pack the self-consistent SCC vector from m_wfn: shell charges for GFN1;
    // shell charges + atomic dipoles + quadrupoles for GFN2 (mirrors what tblite
    // mixes). Claude Generated.
    auto packSCC = [&]() -> Vector {
        if (m_method == MethodType::GFN2) {
            Vector v = Vector::Zero(nsh + 9 * nat);
            v.head(nsh) = m_wfn.q_sh;
            if (m_wfn.dp_at.cols() == nat)
                for (int k = 0; k < 3; ++k)
                    v.segment(nsh + k * nat, nat) = m_wfn.dp_at.row(k).transpose();
            if (m_wfn.qp_at.cols() == nat)
                for (int k = 0; k < 6; ++k)
                    v.segment(nsh + 3 * nat + k * nat, nat) = m_wfn.qp_at.row(k).transpose();
            return v;
        }
        return m_wfn.q_sh;
    };
    // Unpack an SCC vector back into m_wfn (q_at is rederived from q_sh).
    auto unpackSCC = [&](const Vector& v) {
        m_wfn.q_sh = v.head(nsh);
        m_wfn.q_at = Vector::Zero(nat);
        for (int s = 0; s < nsh; ++s)
            m_wfn.q_at(m_basis.sh2at[s]) += m_wfn.q_sh(s);
        if (m_method == MethodType::GFN2) {
            m_wfn.dp_at = Eigen::MatrixXd::Zero(3, nat);
            m_wfn.qp_at = Eigen::MatrixXd::Zero(6, nat);
            for (int k = 0; k < 3; ++k)
                m_wfn.dp_at.row(k) = v.segment(nsh + k * nat, nat).transpose();
            for (int k = 0; k < 6; ++k)
                m_wfn.qp_at.row(k) = v.segment(nsh + 3 * nat + k * nat, nat).transpose();
        }
    };
    // Ensure the GFN2 multipole moments are sized before the first pack.
    if (mode == ScfMode::Broyden && m_method == MethodType::GFN2) {
        if (m_wfn.dp_at.cols() != nat) m_wfn.dp_at = Eigen::MatrixXd::Zero(3, nat);
        if (m_wfn.qp_at.cols() != nat) m_wfn.qp_at = Eigen::MatrixXd::Zero(6, nat);
    }

    const auto t_scf_start = clock::now();
    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        const auto t_iter0 = clock::now();

        // Broyden mixes the SCC charge vector: capture the input x_in (current
        // m_wfn charges) before the potential is built and the populations are
        // overwritten by the diagonalisation.
        Vector x_in;
        if (mode == ScfMode::Broyden)
            x_in = packSCC();

        // Reset and build potentials
        m_pot.reset();

        // Isotropic Coulomb: v_sh += γ * q_sh
        addCoulombShellPotential(m_pot);

        // Third-order: add to v_at/v_sh
        addThirdOrderPotential(m_pot);

        // GFN2 multipole
        if (m_method == MethodType::GFN2) {
            addMultipolePotential(m_pot);
            // Self-consistent D4: exact per-reference dE_D4/dq into the atom potential.
            addDispersionPotential(m_pot);
        }

        // Expand to AO potential and build Fock
        Matrix F = buildFock(m_H0, m_S, m_pot);

        // --- SCF acceleration strategy (Claude Generated) -----------------
        switch (mode) {
        case ScfMode::Plain:
            // No DIIS, no level shift. Convergence comes entirely from the
            // density mixing applied after the solve below.
            break;

        case ScfMode::Broyden:
            // The Fock matrix is left untouched; Broyden mixes the SCC charge
            // vector after the populations are computed (end of the loop body).
            break;

        case ScfMode::Diis:
            // DIIS only after the damped warmup, and never from a diverging
            // Fock — so the early iterations cannot extrapolate the density out
            // of the right basin.
            if (iter >= diis_start) {
                diis.push(F, m_wfn.P, m_S);
                if (diis.size() >= 2 && diis.lastErrorNorm() < kDiisErrorCutoff)
                    F = diis.extrapolate();
            }
            break;

        case ScfMode::LevelShift:
            // Virtual-orbital level shift layered on top of the density mixing
            // applied after the solve. The shift raises the empty orbitals so
            // the density responds less violently to the Fock update; it is
            // faded out once max|dq| has settled, so the converged eps / density
            // are unshifted. No DIIS — the shift + mixing are the stabiliser.
            if (iter > 0 && P_old.size() > 0 && dq_prev > 1.0e-3)
                F = applyLevelShift(F, m_S, P_old, m_level_shift);
            break;
        }

        // Diagonalize
        if (!solveEigen(F, m_S)) {
            CurcumaLogger::warn("XTB::Calculation: eigen solve failed at iteration " + std::to_string(iter));
            m_scf_converged = false;
            m_scf_iterations = iter;
            return m_E_total;
        }

        // Density damping: P = damp·P_new + (1-damp)·P_old. Plain and LevelShift
        // mix every iteration (Plain's only convergence mechanism; LevelShift
        // layers the Fock shift on top); DIIS mixes only during the warmup, then
        // lets the extrapolation take over. Broyden never mixes the density — it
        // mixes the SCC charge vector at the end of the loop body instead.
        const bool mix = (mode != ScfMode::Broyden)
                         && (mode == ScfMode::Plain || mode == ScfMode::LevelShift
                             || iter < diis_start);
        if (iter > 0 && mix) {
            m_wfn.P = damp * m_wfn.P + (1.0 - damp) * P_old;
        }
        P_old = m_wfn.P;

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
        dq_prev = dq;   // drives the LevelShift fade-out on the next iteration
        const double de = (iter > 0) ? std::fabs(e_scc - e_total_old) : 0.0;
        const double t_iter_ms = ms(t_iter0, clock::now());
        if (verb >= 2) {
            std::string line = fmt::format("  {:4d}   {:18.10f}   {:10.2e}   {:10.2e}   {:6.1f}",
                                           iter, e_scc, de, dq, t_iter_ms);
            if (verb >= 3)
                line += fmt::format("   [DIIS n={} err={:.2e}]", diis.size(), diis.lastErrorNorm());
            CurcumaLogger::info(line);
        }

        // Check convergence (after first iteration). On convergence we break
        // here, so m_wfn keeps the consistent output charges x_out (the Broyden
        // mix below is skipped — no re-mixing of the converged state).
        if (iter > 0) {
            if (checkConvergence_impl(q_sh_old, q_sh_new,
                                       e_total_old, e_scc, thresh)) {
                m_scf_converged = true;
                break;
            }
        }

        // Broyden: mix the SCC charge vector to get the next input. The raw
        // output x_out (in m_wfn after updatePopulations) is combined with the
        // input x_in via the modified-Broyden quasi-Newton update; unpacking
        // writes the next input back into m_wfn so q_sh_old below tracks it.
        if (mode == ScfMode::Broyden) {
            const Vector x_out  = packSCC();
            const Vector x_next = broyden.update(x_in, x_out);
            unpackSCC(x_next);
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
    m_E_dispersion    = calcDispersionEnergy(gradient);

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
 *  GFN2 component audit (Claude Generated): SCF-free per-container energy
 *  evaluation at an externally supplied wavefunction state. Mirrors the
 *  pre-SCF setup steps from Calculation() and the post-SCF energy
 *  accumulators, but skips diagonalisation. See xtb_native.h for the
 *  contract; intended to be driven by diag_curcuma_energy_components.
 * ------------------------------------------------------------------------- */
bool XTB::evaluateComponentsAtFixedDensity(
    const Matrix& P,
    const Vector& q_at,
    const Vector& q_sh,
    const Eigen::MatrixXd& dp_at,
    const Eigen::MatrixXd& qp_at)
{
    const int nat = m_basis.nat;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;
    if (nat <= 0 || nsh <= 0 || nao <= 0) {
        CurcumaLogger::error("XTB::evaluateComponentsAtFixedDensity: basis not built; call InitialiseMolecule() first");
        return false;
    }
    if (P.rows() != nao || P.cols() != nao) {
        CurcumaLogger::error_fmt("XTB::evaluateComponentsAtFixedDensity: P shape ({}x{}) != nao ({})",
                                 P.rows(), P.cols(), nao);
        return false;
    }
    if (q_at.size() != nat) {
        CurcumaLogger::error_fmt("XTB::evaluateComponentsAtFixedDensity: q_at size {} != nat {}",
                                 q_at.size(), nat);
        return false;
    }
    if (q_sh.size() != nsh) {
        CurcumaLogger::error_fmt("XTB::evaluateComponentsAtFixedDensity: q_sh size {} != nsh {}",
                                 q_sh.size(), nsh);
        return false;
    }

    // Pre-SCF setup (mirrors Calculation() steps 1-5 minus the SCF loop).
    Vector cn = computeCoordinationNumbers();
    Vector se;
    getSelfEnergies(cn, se);
    Matrix S, H0;
    getHamiltonianH0(se, S, H0);
    m_S  = S;
    m_H0 = H0;
    buildGammaMatrix();
    if (m_method == MethodType::GFN2) {
        setupMultipole();
    }
    m_coordination_numbers = cn;

    // Inject wavefunction state. q_at/q_sh and the dual n_at/n_sh stay
    // consistent (n = n0 - q) so any helper that reads either form sees the
    // same physical state.
    m_wfn.P    = P;
    m_wfn.q_at = q_at;
    m_wfn.q_sh = q_sh;
    m_wfn.n_at = m_wfn.n0_at - q_at;
    m_wfn.n_sh = m_wfn.n0_sh - q_sh;
    if (m_method == MethodType::GFN2) {
        if (dp_at.rows() != 3 || dp_at.cols() != nat) {
            CurcumaLogger::error_fmt("XTB::evaluateComponentsAtFixedDensity: dp_at shape ({}x{}) != (3x{})",
                                     dp_at.rows(), dp_at.cols(), nat);
            return false;
        }
        if (qp_at.rows() != 6 || qp_at.cols() != nat) {
            CurcumaLogger::error_fmt("XTB::evaluateComponentsAtFixedDensity: qp_at shape ({}x{}) != (6x{})",
                                     qp_at.rows(), qp_at.cols(), nat);
            return false;
        }
        m_wfn.dp_at = dp_at;
        m_wfn.qp_at = qp_at;
    } else {
        m_wfn.dp_at = Eigen::MatrixXd::Zero(3, nat);
        m_wfn.qp_at = Eigen::MatrixXd::Zero(6, nat);
    }

    // Evaluate every per-container energy at the injected state. Order matches
    // Calculation()'s final-energy block (lines 260-268) so the m_E_* fields
    // and the public getters report the same quantities, just at the
    // *injected* density rather than at curcuma's own SCF fixpoint.
    m_E_coulomb_shell = energyCoulombShell();
    m_E_third_order   = energyThirdOrder();
    m_E_multipole     = energyMultipole();
    const double e_band = (m_wfn.P.cwiseProduct(m_H0)).sum();
    m_E_electronic = e_band + m_E_coulomb_shell + m_E_third_order + m_E_multipole;
    m_E_repulsion    = calcRepulsionEnergy();
    m_E_halogen_bond = calcHalogenBondEnergy();
    // Skip the dispersion CPSCF charge-response fold: it needs MO coefficients
    // (m_wfn.C / m_wfn.eps) which are absent because we did not diagonalise.
    // The dispersion energy itself is unaffected.
    m_disp_audit_mode = true;
    m_E_dispersion   = calcDispersionEnergy();
    m_disp_audit_mode = false;
    m_E_total = m_E_electronic + m_E_repulsion
              + m_E_halogen_bond + m_E_dispersion;

    return true;
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
double XTB::calcDispersionEnergy(bool need_gradient) const
{
    m_disp_gradient_valid = false;

    if (m_method != MethodType::GFN2) {
        // GFN1: native D3(BJ) dispersion (s6=1.0, s8=2.4, a1=0.63, a2=5.0).
        // D3ParameterGenerator exposes only the energy, so the geometry
        // gradient is obtained by central finite differences (D3 is cheap;
        // this mirrors the legacy gfn1.cpp path). m_geometry is in Angstrom
        // (the generator converts to Bohr internally). The FD gradient is
        // returned in Eh/Bohr to match the rest of calculateGradient(), which
        // divides the assembled gradient by `au` (Eh/Bohr -> Eh/Angstrom) at
        // the end. The FD is a *full* gradient: regenerating parameters at each
        // displaced geometry also captures the CN-dependence of the C6 values.
        if (!m_d3_generator) {
            m_d3_generator = std::make_unique<::D3ParameterGenerator>(
                ::D3ParameterGenerator::createForGFN1());
        }
        ::D3ParameterGenerator& d3 = *m_d3_generator;
        d3.GenerateParameters(m_atoms, m_geometry);
        const double e_disp = d3.getTotalEnergy();

        // The central-difference geometry gradient regenerates the D3 parameters
        // at 6·N displaced geometries — O(N) full D3 evaluations, prohibitively
        // expensive on large systems (≈ 1400 evaluations × 26k pairs for the
        // 231-atom complex). Compute it only when a gradient is actually
        // requested; a single-point energy never reads m_disp_gradient.
        if (!need_gradient) {
            m_disp_gradient = Matrix::Zero(m_atomcount, 3);
            m_disp_dEdcn = Vector();
            m_disp_gradient_valid = false;
            return e_disp;
        }

        // Central-difference geometry gradient (Eh/Bohr).
        const double h = 1.0e-4;  // Angstrom displacement
        m_disp_gradient = Matrix::Zero(m_atomcount, 3);
        Matrix geom = m_geometry;
        for (int a = 0; a < m_atomcount; ++a) {
            for (int c = 0; c < 3; ++c) {
                const double orig = geom(a, c);
                geom(a, c) = orig + h;
                d3.GenerateParameters(m_atoms, geom);
                const double ep = d3.getTotalEnergy();
                geom(a, c) = orig - h;
                d3.GenerateParameters(m_atoms, geom);
                const double em = d3.getTotalEnergy();
                geom(a, c) = orig;
                // dE/dAngstrom -> dE/dBohr : multiply by au (1 Bohr = au Angstrom)
                m_disp_gradient(a, c) = (ep - em) / (2.0 * h) * au;
            }
        }
        // Restore unperturbed parameters for any later query.
        d3.GenerateParameters(m_atoms, m_geometry);
        // The FD gradient is complete (it already differentiates the CN-dependent
        // C6 values), so there is no separate dE/dCN chain-rule term: keep
        // m_disp_dEdcn empty so the GFN2 CN-fold in xtb_gradient.cpp is skipped.
        m_disp_dEdcn = Vector();
        m_disp_gradient_valid = true;
        return e_disp;
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
        // AP6b: native GFN2 uses the exact dftd4 per-reference charge weighting.
        p.per_reference_charge = true;
        m_d4_evaluator = std::make_unique<curcuma::dispersion::D4Evaluator>(
            m_d4_generator.get(), p);
    }

    // AP6b: the exact per-reference GFN2 D4 (weightedC6Gfn2) drives the C6 charge
    // weighting with the converged SCF Mulliken charges (tblite uses wfn%qat). The
    // dftd4 single-shot EEQ is disabled; getZetaCharges() returns the Mulliken set.
    // The dq/dx gradient response then comes from CPSCF (computeMullikenChargeResponse).
    m_d4_generator->setUseD4SingleShotEEQ(false);
    m_d4_generator->setTopologyCharges(m_wfn.q_at);
    // GFN2: use the dftd4 EN-weighted covalent CN (matches tblite's
    // get_coordination_number(rcov, en)) instead of the GFN-FF log-capped erf-CN.
    m_d4_generator->setD4CovalentCN(true);

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
    const bool want_dEdq = true;
    double E = m_d4_evaluator->computeEnergyAndGradient(
        m_atoms, geom_bohr,
        /*with_gradient=*/true,
        m_disp_gradient, m_disp_dEdcn, m_disp_dEdq, want_dEdq);

    // ATM three-body term (GFN2 s9=5.0). Charge-independent (dftd4 evaluates it at
    // the q=0 reference C6 — matches tblite get_dispersion_nonsc), so no q-response.
    // ACCUMULATES into m_disp_gradient and m_disp_dEdcn (the latter folded together
    // with the 2-body dEdcn via the D4 covalent CN below). Closes the triose C-path
    // remainder (see docs/GFN2_D4_STATUS.md). Claude Generated 2026-05-29.
    E += m_d4_evaluator->computeATM(
        m_atoms, geom_bohr, /*with_gradient=*/true, m_disp_gradient, m_disp_dEdcn);

    // CN chain rule: fold Σ_A dE_D4/dCN_A · ∂CN_A/∂R into the cached gradient using
    // the D4 CN's OWN derivative (EN-weighted erf), keeping the dispersion gradient
    // self-consistent with the D4 CN. m_disp_dEdcn is then cleared so the GFN2
    // Hamiltonian CN loop in calculateGradient() does not double-distribute it
    // through the (different) logistic Hamiltonian-CN derivative.
    if (m_disp_dEdcn.size() == m_atomcount) {
        curcuma::dispersion::addD4CovalentCNGradient(
            m_atoms, geom_bohr, m_disp_dEdcn, m_disp_gradient);
        m_disp_dEdcn = Vector();
    }

    // q-response chain rule: fold Σ_A dE_D4/dq_A · ∂q_A/∂R into the cached geometry
    // gradient. The per-reference path is Mulliken-self-consistent, so ∂q/∂x comes
    // from the GFN2 CPSCF/Z-vector response. Skipped in audit mode (no MO state).
    if (m_disp_dEdq.size() == m_atomcount && !m_disp_audit_mode) {
        computeMullikenChargeResponse(m_disp_dEdq, m_disp_gradient);
    }
    m_disp_gradient_valid = true;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dispersion_d4_native", fmt::format("{:.6f} Eh", E));
    }
    return E;
}

/* ------------------------------------------------------------------------- *
 *  Self-consistent D4 dispersion potential (GFN2, AP6b)
 *
 *  GFN2-xTB couples D4 to the SCF: the C6 are charge-scaled (per-reference
 *  zeta of the SCF Mulliken charges), so E_D4 = E_D4(q) and dE_D4/dq_A is an
 *  extra atom potential in the Fock. tblite adds exactly this (`dispersion%
 *  get_potential`: pot%vat += dE_disp/dq) every SCF iteration. We evaluate the
 *  exact per-reference dE_D4/dq at m_wfn.q_at and add it to pot.v_at (the same
 *  atom potential the multipole back-reaction writes; verified bit-for-bit vs
 *  tblite's pot%vat). GFN2 only. Claude Generated.
 * ------------------------------------------------------------------------- */
void XTB::addDispersionPotential(Potential& pot) const
{
    if (m_method != MethodType::GFN2) return;

    // Lazy init shared with calcDispersionEnergy (idempotent).
    if (!m_d4_generator) {
        nlohmann::json d4_config;
        d4_config["d4_s6"] = 1.0; d4_config["d4_s8"] = 2.7;
        d4_config["d4_a1"] = 0.52; d4_config["d4_a2"] = 5.0; d4_config["d4_alp"] = 16.0;
        ConfigManager cfg("d4param", d4_config);
        m_d4_generator = std::make_unique<::D4ParameterGenerator>(cfg);
    }
    if (!m_d4_evaluator) {
        curcuma::dispersion::D4Params p;
        p.s6 = 1.0; p.s8 = 2.7; p.a1 = 0.52; p.a2 = 5.0; p.s9 = 5.0; p.alpha = 16.0;
        p.damping = curcuma::dispersion::DampingFormula::StandardBJ_D4;
        p.per_reference_charge = true;
        m_d4_evaluator = std::make_unique<curcuma::dispersion::D4Evaluator>(
            m_d4_generator.get(), p);
    }

    // Exact per-reference dE_D4/dq at the current SCF Mulliken charges.
    m_d4_generator->setUseD4SingleShotEEQ(false);
    m_d4_generator->setTopologyCharges(m_wfn.q_at);
    m_d4_generator->setD4CovalentCN(true);  // GFN2 dftd4 EN-weighted covalent CN
    const Matrix geom_bohr = m_geometry * AA_TO_AU;
    m_d4_generator->GenerateParameters(m_atoms, geom_bohr);

    Matrix scratch_grad; Vector scratch_dEdcn, dEdq;
    m_d4_evaluator->computeEnergyAndGradient(
        m_atoms, geom_bohr, /*with_gradient=*/true,
        scratch_grad, scratch_dEdcn, dEdq, /*with_dEdq=*/true);

    if (dEdq.size() != m_atomcount) return;
    for (int A = 0; A < m_atomcount; ++A)
        pot.v_at(A) += dEdq(A);
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
