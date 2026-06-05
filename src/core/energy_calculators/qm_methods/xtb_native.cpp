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
#include "src/core/intra_parallel_context.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/energy_calculators/dispersion/d4_ncoord.h"  // D4 CN gradient (GFN2)
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"  // D3 CN gradient (May 2026)
#include "diis_accelerator.h"
#include "broyden_mixer.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <chrono>
#include <future>
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

/* ------------------------------------------------------------------------- *
 *  Intra-molecule parallelism helpers (Claude Generated).
 *
 *  Curcuma parallelises at the molecule level (one CxxThreadPool worker per
 *  molecule in MD/conformer/batch runs), so these honour the intra-SCF thread
 *  budget ONLY for a single large molecule and otherwise stay serial. See the
 *  auto-gate rationale in xtb_native.h (CxxThreadPoolWorkerFlag, MklThreadScope).
 * ------------------------------------------------------------------------- */
int XTB::effectiveIntraThreads(int work_units) const
{
    if (m_intra_threads <= 1) return 1;
    // Already running under molecule-level parallelism → stay serial (no N^2).
    if (curcuma::intraParallelSuppressed()) return 1;
    // Size guard: keep at least kMinWorkPerThread items per thread; tiny systems
    // are dominated by dispatch overhead (tuned in the benchmark step).
    const int by_work = work_units / kMinWorkPerThread;
    const int t = std::min(m_intra_threads, std::max(1, by_work));
    return t;
}

void XTB::ensurePool(int n_threads) const
{
    if (n_threads <= 1) return;
    if (!m_pool) {
        m_pool = std::make_unique<CxxThreadPool>();
        m_pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    }
    // setActiveThreadCount is idempotent: it only rebuilds workers when the count
    // actually changes, so calling it every region is cheap.
    m_pool->setActiveThreadCount(n_threads);
}

void XTB::parallelStripes(int n_threads,
                          const std::function<void(int, int)>& worker) const
{
    if (n_threads <= 1) { worker(0, 1); return; }
    ensurePool(n_threads);
    std::vector<std::future<void>> futures;
    futures.reserve(n_threads - 1);
    for (int t = 1; t < n_threads; ++t)
        futures.push_back(m_pool->enqueue(worker, t, n_threads));
    worker(0, n_threads);                 // main thread does stripe 0
    for (auto& f : futures) f.get();       // join
}

bool XTB::InitialiseMolecule()
{
    if (m_atomcount <= 0) {
        CurcumaLogger::error("XTB::InitialiseMolecule: no atoms set");
        return false;
    }
    buildBasis();
    buildH0Data();
    buildReferenceOccupations();
    // A fresh molecule invalidates any SCC-extrapolation history from a previous
    // system (different nsh/nat). Claude Generated.
    m_scf_history.clear();
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

    // Stage 5 (Part A): on the device-resident path solve the single-shot EEQ on
    // the GPU (self-contained — no resident SCF state needed yet), so the per-step
    // guess does not pull back to the host. resolveParams is the shared parameter
    // source; the host D4ChargeModel is the fallback on any device failure.
    Vector q_at;
    bool got_device = false;
    if (m_gpu_scf && m_gpu_scf->supportsDeviceEeq()) {
        std::vector<double> chi, gam, alp, cnf, rcov;
        curcuma::dispersion::D4ChargeModel::resolveParams(m_atoms, chi, gam, alp, cnf, rcov);
        std::vector<double> xyz(3 * nat), q(nat, 0.0);
        for (int i = 0; i < nat; ++i) {
            xyz[3 * i + 0] = geom_bohr(i, 0);
            xyz[3 * i + 1] = geom_bohr(i, 1);
            xyz[3 * i + 2] = geom_bohr(i, 2);
        }
        if (m_gpu_scf->eeqCharges(nat, xyz.data(), chi.data(), gam.data(), alp.data(),
                                  cnf.data(), rcov.data(), static_cast<double>(m_charge),
                                  q.data())) {
            q_at = Vector::Zero(nat);
            for (int i = 0; i < nat; ++i) q_at(i) = q[i];
            got_device = q_at.allFinite();
        }
    }
    if (!got_device) {
        curcuma::dispersion::D4ChargeModel eeq;
        q_at = eeq.computeCharges(m_atoms, geom_bohr, static_cast<double>(m_charge));
    }
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

/* ------------------------------------------------------------------------- *
 *  Multi-step SCC extrapolation helpers (Claude Generated).
 *
 *  packSccState / unpackSccState (de)serialise the SCC mixing vector — the
 *  quantity the SCF iterates on. GFN1: shell charges q_sh only. GFN2: the full
 *  vector [q_sh; vec(dp_at); vec(qp_at)] (atomic dipoles 3×nat and quadrupoles
 *  6×nat, column-major flatten), matching exactly what the 1-step warm-start
 *  saves/restores (see xtb_native.cpp warm-start block + UpdateMolecule).
 * ------------------------------------------------------------------------- */
Eigen::VectorXd XTB::packSccState() const
{
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;
    // GFN1 (SCC vector IS q_sh), or a GFN2 state whose multipoles are not yet sized
    // (degenerate): return just q_sh. A shorter-than-packed_len vector is rejected by
    // the size guard in Calculation(), which then falls back to the 1-step warm-start.
    if (m_method != MethodType::GFN2
        || m_wfn.dp_at.cols() != nat || m_wfn.qp_at.cols() != nat) {
        return m_wfn.q_sh;
    }
    Eigen::VectorXd v(nsh + 3 * nat + 6 * nat);
    v.head(nsh) = m_wfn.q_sh;
    // dp_at is 3×nat, qp_at is 6×nat; Eigen is column-major, so .data() is the
    // column-major flatten that Eigen::Map reads back identically.
    v.segment(nsh, 3 * nat)            = Eigen::Map<const Eigen::VectorXd>(m_wfn.dp_at.data(), 3 * nat);
    v.segment(nsh + 3 * nat, 6 * nat)  = Eigen::Map<const Eigen::VectorXd>(m_wfn.qp_at.data(), 6 * nat);
    return v;
}

void XTB::unpackSccState(const Eigen::VectorXd& v)
{
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;
    m_wfn.q_sh = v.head(nsh);
    // Rebuild atomic charges from shell charges (as the warm-start block does).
    m_wfn.q_at.setZero(nat);
    for (int s = 0; s < nsh; ++s)
        m_wfn.q_at(m_basis.sh2at[s]) += m_wfn.q_sh(s);
    if (m_method == MethodType::GFN2 && v.size() == nsh + 9 * nat) {
        m_wfn.dp_at = Eigen::Map<const Eigen::MatrixXd>(v.data() + nsh, 3, nat);
        m_wfn.qp_at = Eigen::Map<const Eigen::MatrixXd>(v.data() + nsh + 3 * nat, 6, nat);
    }
}

/* ------------------------------------------------------------------------- *
 *  Charge-conserving extrapolation weights w_j for P_pred = Σ_j w_j P_hist[j],
 *  history[0] = newest. Σ w_j = 1 (so the predicted SCC vector keeps the total
 *  charge / trace of the history). Returns an empty vector when there is too
 *  little history (caller then keeps the 1-step warm-start). Claude Generated.
 *
 *  aspc : Always Stable Predictor-Corrector (Kolafa 2004). m = min(k+2, n_hist)
 *         points, w_j = (-1)^j C(m, j+1). Time-reversible for fixed-step MD.
 *  gauss: least-squares polynomial of degree d = min(order, m-1) fit to the
 *         history at abscissae x_j = -(j+1) and extrapolated to x = 0. With
 *         w = V·a where (VᵀV)a = e_0 the value at x=0 is Σ w_j y_j; robust for
 *         the irregular steps of a geometry optimisation.
 * ------------------------------------------------------------------------- */
Eigen::VectorXd XTB::extrapolationWeights(const std::string& mode, int order, int n_hist)
{
    if (n_hist < 2)
        return {};

    if (mode == "aspc") {
        const int m = std::min(order + 2, n_hist);
        if (m < 2) return {};
        // Binomial C(m, i) via the multiplicative recurrence (exact for small m).
        Eigen::VectorXd w(m);
        long long c = 1; // C(m,0)
        for (int j = 0; j < m; ++j) {
            // c currently holds C(m, j); advance to C(m, j+1)
            c = c * (m - j) / (j + 1);
            w(j) = (j % 2 == 0 ? 1.0 : -1.0) * static_cast<double>(c);
        }
        return w;
    }

    if (mode == "gauss") {
        const int m = n_hist;
        const int d = std::min(order, m - 1); // polynomial degree actually fit
        if (d < 1) return {};
        Eigen::MatrixXd V(m, d + 1);
        for (int j = 0; j < m; ++j) {
            const double x = -static_cast<double>(j + 1);
            double xp = 1.0;
            for (int p = 0; p <= d; ++p) { V(j, p) = xp; xp *= x; }
        }
        Eigen::VectorXd e0 = Eigen::VectorXd::Zero(d + 1);
        e0(0) = 1.0;
        Eigen::MatrixXd VtV = V.transpose() * V;
        Eigen::VectorXd a = VtV.ldlt().solve(e0);
        Eigen::VectorXd w = V * a;
        if (!w.allFinite())
            return {};
        return w;
    }

    return {}; // unknown mode -> fall back to warm-start
}

double XTB::Calculation(bool gradient)
{
    // Pin MKL to one thread for the whole native SCF (see MklSerialScope in
    // xtb_native.h): the serial iteration over small/medium matrices is faster
    // without MKL's per-call thread team up to at least 231 atoms.
    MklSerialScope mkl_serial;

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

    // Claude Generated (2026-05): new geometry -> the geometry-dependent D4 reference
    // set (CN, Gaussian weights, C6 cache) must be regenerated once this Calculation.
    // Per-SCF addDispersionPotential() then only refreshes the SCF charges. Resetting
    // here is REQUIRED for opt/MD correctness: without it later geometries would reuse
    // the first geometry's stale D4 CN/weights.
    m_d4_prepared = false;
    m_d4_genparams_calls = 0;

    // 1. Coordination numbers
    Vector cn = computeCoordinationNumbers();
    const auto t_cn = clock::now();

    // 2. Self-energies
    Vector se;
    getSelfEnergies(cn, se);

    // 3+4. Overlap S, bare Hamiltonian H0, Cholesky L (=m_X), Coulomb gamma.
    // AP4 (Claude Generated): on the GPU device-resident path the device builds these
    // integrals itself (beginComputed), so build them on the device FIRST and download
    // S/H0/L/gamma into the host matrices, SKIPPING the redundant host
    // getHamiltonianH0/buildOrthonormalizer/buildGammaMatrix (~27 ms/geometry — matters
    // on -opt/-md). The device matrices match the CPU to ~1e-15 (component tests
    // sqm_cuda_*). exportGpuBasis/beginBasis/beginComputed are hoisted here (out of the
    // SCF-setup block below) and their result (gpu_computed) is reused there, so the
    // device build runs exactly once. Any failure / non-resident mode → full host build.
    Matrix S, H0;
    bool gpu_computed = false;        // device integral build succeeded (reused below)
    bool integrals_from_device = false;
    if (m_gpu_scf && m_scf_mode == ScfMode::Broyden) {
        GpuBasisFlat gbf;
        GpuH0Flat    gh0;
        exportGpuBasis(gbf, gh0);
        bool basis_ok = true;
        if (m_gpu_basis_dirty) {
            basis_ok = m_gpu_scf->beginBasis(gbf, gh0);
            if (basis_ok) m_gpu_basis_dirty = false;
        }
        gpu_computed = basis_ok && m_gpu_scf->beginComputed(gbf.xyz_bohr);
        if (gpu_computed) {
            const int nao = m_basis.nao, nsh = m_basis.nsh;
            Eigen::MatrixXd Scm(nao, nao), H0cm(nao, nao), Lcm(nao, nao), Gcm(nsh, nsh);
            if (m_gpu_scf->downloadOverlap(Scm) && m_gpu_scf->downloadH0(H0cm)
                && m_gpu_scf->downloadCholesky(Lcm) && m_gpu_scf->downloadGamma(Gcm)) {
                m_S     = Scm;   // col-major device → row-major Matrix (S symmetric)
                m_H0    = H0cm;  // H0 symmetric
                m_X     = Lcm;   // m_X is column-major Eigen::MatrixXd → direct
                m_gamma = Gcm;
                integrals_from_device = true;
            }
        }
    }
    if (!integrals_from_device) {
        // Full host integral build (CPU path, or device build unavailable).
        getHamiltonianH0(se, S, H0);
        m_S  = S;
        m_H0 = H0;
        // Orthonormalizer X = S^{-1/2}, built once here so every SCF iteration solves
        // the cheap standard eigenproblem instead of re-factorizing the constant S.
        buildOrthonormalizer();
        buildGammaMatrix();   // Coulomb gamma (once per geometry)
    }
    const auto t_h0 = clock::now();
    const auto t_gamma = clock::now();   // S/H0/L/gamma fused above (AP4 device path)

    // 5. Multipole setup (GFN2 only) — fills m_dp_int, m_qp_int, interaction matrices.
    // Always on the host: the AP4 device path uploads dp_int/qp_int, it does not build
    // the AO multipole integrals.
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

    // Context-aware verbosity: in iterative mode (MD/opt) raise the display
    // threshold by one so -verbosity 2 is needed to see per-iteration output.
    // SP mode uses -verbosity 1 as the normal threshold. Claude Generated.
    const int scf_min = m_is_iterative ? 2 : 1;

    if (verb >= scf_min) {
        CurcumaLogger::result(fmt::format("SCF iterations (mode={}, guess={}):",
                                          scfModeName(m_scf_mode), m_scf_guess));
        CurcumaLogger::result("  iter             E(SCC)/Eh            dE        max|dq|     t/ms");
    }

    const int max_iter   = m_scf_max_iter;
    const double thresh  = m_scf_threshold;
    const int nsh = m_basis.nsh;

    // Initial guess. Warm-start: reuse converged shell charges from the previous
    // geometry (saved by UpdateMolecule). Otherwise: "h0" (zero) or "eeq".
    // Claude Generated warm-start path.
    Vector q_sh_old = Vector::Zero(nsh);
    bool guess_set = false;

    // (0) Multi-step SCC extrapolation (opt-in, guess mode). Predict the new-step
    //     SCC vector from the history of converged steps so the SCF starts closer
    //     to the new fixpoint than the previous step alone (fewer iterations, esp.
    //     in MD). Only the initial guess changes; the SCF below still converges to
    //     the same fixpoint within scf_threshold. A poor prediction or too-short
    //     history simply falls through to the 1-step warm-start. Claude Generated.
    const int packed_len = (m_method == MethodType::GFN2) ? (nsh + 9 * m_atomcount) : nsh;
    // XL-BOMD corrector-only coupling is Phase 2 / not yet implemented; warn once and
    // fall back to the safe full-SCF guess mode so the result stays correct. Claude Generated.
    if (m_scf_extrapolation != "none" && m_scf_extrap_apply == "xlbomd" && !m_xlbomd_warned) {
        m_xlbomd_warned = true;
        CurcumaLogger::warn("scf_extrapolation_apply=xlbomd is not yet implemented; "
                            "using the safe full-SCF 'guess' coupling instead");
    }
    if (m_warmstart && m_scf_extrapolation != "none"
        && static_cast<int>(m_scf_history.size()) >= 2) {
        Eigen::VectorXd w = extrapolationWeights(m_scf_extrapolation, m_scf_extrap_order,
                                                 static_cast<int>(m_scf_history.size()));
        if (w.size() >= 2) {
            Eigen::VectorXd pred = Eigen::VectorXd::Zero(packed_len);
            bool ok = true;
            for (int j = 0; j < w.size(); ++j) {
                if (m_scf_history[j].size() != packed_len) { ok = false; break; }
                pred.noalias() += w(j) * m_scf_history[j];
            }
            if (ok && pred.allFinite()) {
                unpackSccState(pred);
                q_sh_old  = m_wfn.q_sh;
                guess_set = true;
                // Cite the extrapolation scheme actually used (deduplicated, so
                // calling it every step is a no-op after the first). Claude Generated.
                CurcumaLogger::citation(m_scf_extrapolation == "aspc"
                                            ? "aspc" : "density_extrapolation");
                if (verb >= scf_min)
                    CurcumaLogger::result(fmt::format(
                        "SCF initial guess: {} extrapolation over {} steps",
                        m_scf_extrapolation, static_cast<int>(w.size())));
            }
        }
    }

    if (guess_set) {
        // extrapolation already populated m_wfn / q_sh_old
    } else if (m_warmstart && m_warmstart_q_sh.size() == nsh) {
        q_sh_old   = m_warmstart_q_sh;
        m_wfn.q_sh = m_warmstart_q_sh;
        m_wfn.q_at.setZero(m_atomcount);
        for (int s = 0; s < nsh; ++s)
            m_wfn.q_at(m_basis.sh2at[s]) += m_warmstart_q_sh(s);
        // GFN2: restore the converged atomic multipoles too (the rest of the SCC
        // vector), so the first Fock is built from the full previous-step state rather
        // than q_sh + zero multipoles. This is what lets the GFN2 opt/md SCF drop toward
        // the GFN1-like few-iteration tail near convergence. Claude Generated.
        if (m_method == MethodType::GFN2
            && m_warmstart_dp_at.cols() == m_atomcount
            && m_warmstart_qp_at.cols() == m_atomcount) {
            m_wfn.dp_at = m_warmstart_dp_at;
            m_wfn.qp_at = m_warmstart_qp_at;
        }
        if (verb >= scf_min)
            CurcumaLogger::result("SCF initial guess: warm-start from previous step");
    } else if (m_scf_guess == "eeq") {
        Vector q_sh_guess;
        if (seedEEQGuess(q_sh_guess)) {
            q_sh_old   = q_sh_guess;
            m_wfn.q_sh = q_sh_guess;
            m_wfn.q_at.setZero(m_atomcount);
            for (int s = 0; s < nsh; ++s)
                m_wfn.q_at(m_basis.sh2at[s]) += q_sh_guess(s);
            if (verb >= scf_min)
                CurcumaLogger::result("SCF initial guess: single-shot EEQ shell charges");
        } else if (verb >= scf_min) {
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

    // DIIS and Broyden are member variables so history can optionally survive
    // across geometry steps (m_keep_diis=true, set via -keep_diis true). Default
    // (m_keep_diis=false): always reset — old error vectors from a different
    // geometry typically slow convergence. Charge warm-start (m_warmstart_q_sh)
    // is independent of this flag and always reused when available. Claude Generated.
    constexpr double kDiisErrorCutoff = 1.0e3;
    if (!m_keep_diis) {
        m_diis    = DIISAccelerator(m_diis_subspace);
        m_broyden = BroydenMixer(damp, /*max_history=*/m_diis_subspace > 2 ? 20 : m_diis_subspace);
    }
    // Local aliases for brevity in the SCF loop body.
    DIISAccelerator& diis    = m_diis;
    BroydenMixer&    broyden = m_broyden;
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

    // Intra-iteration timing buckets (printed at verbosity >= 3). A handful of
    // steady_clock reads per phase — negligible vs the ms-scale work — so they
    // run unconditionally, matching the existing per-iter timer. Claude Generated.
    double acc_pot = 0.0, acc_fock = 0.0, acc_solve = 0.0, acc_mull = 0.0, acc_energy = 0.0;
    double acc_disp = 0.0;   // D4 in-SCF potential subset of acc_pot (GFN2), verbosity 3
    m_t_xfx = m_t_diag = m_t_back = m_t_dens = 0.0;

    // Intra-molecule thread count for the per-iteration eigensolve. The eigensolve
    // (dsygst/dsyevd/dtrsm) is the one region handed to MKL rather than the
    // CxxThreadPool; effectiveIntraThreads() gates it (serial under molecule-level
    // parallelism or for small bases). nao is constant over the SCF. Claude Generated.
    const int eig_threads = effectiveIntraThreads(m_basis.nao);
    if (verb >= 3 && eig_threads > 1)
        CurcumaLogger::info_fmt("Eigensolve MKL threads: {}", eig_threads);

    // Device-resident SCF (Claude Generated, GPU port Stage 2). Enabled with the
    // default Broyden charge mixing and an available lower Cholesky factor L
    // (m_X). begin() uploads the geometry-constant H0/S/L once; the per-iteration
    // heavy linear algebra (Fock build, eigensolve, density, Mulliken-AO, band)
    // then runs on the device with those matrices resident, so only length-nao
    // vectors cross the bus each iteration. GFN2 additionally keeps the multipole
    // integrals resident and adds the anisotropic Fock + atomic moments on the
    // device (Stage 2b); the isotropic potential (third-order, multipole scalar
    // shift, in-SCF D4) is still folded into v_ao on the host. On any failure we
    // silently keep the full CPU SCF. Non-Broyden modes keep the per-iteration
    // GPU eigensolver hook (Stage 1) instead.
    bool use_gpu_resident = false;
    bool gpu_multipole    = false;
    if (m_gpu_scf && mode == ScfMode::Broyden
        && m_X.rows() == m_basis.nao && m_X.cols() == m_basis.nao) {
        // Stage 3: the device-computed integral build (beginBasis once + beginComputed
        // per geometry: CN/S/H0/L on the device, no nao² upload) was hoisted into the
        // integral-build step above (AP4) so S/H0/L/gamma could be downloaded in place
        // of the redundant host build; reuse its result here (the device build runs
        // exactly once). Fall back to the Stage-2 begin() upload path if the device
        // build was unavailable. The device S/H0/L match the CPU to ~1e-15 (component
        // tests sqm_cuda_*), so the converged energy is unchanged. Claude Generated.
        const bool computed = gpu_computed;
        if (m_method == MethodType::GFN1) {
            use_gpu_resident = computed || m_gpu_scf->begin(m_H0, m_S, m_X);
            if (verb >= 2)
                CurcumaLogger::info(!use_gpu_resident
                    ? "SCF: GPU resident begin() failed; running CPU SCF"
                    : (computed
                        ? "SCF: device-resident GFN1 path (CN/S/H0/L built on GPU; no nao^2 upload)"
                        : "SCF: device-resident GFN1 path (H0/S/L uploaded; per-iter vectors only)"));
        } else if (m_method == MethodType::GFN2 && m_mp_initialized
                   && m_gpu_scf->supportsMultipole()) {
            const bool base = computed || m_gpu_scf->begin(m_H0, m_S, m_X);
            // Stage 3d: when the device built the integrals, it also built the
            // multipole dp_int/qp_int — enter the resident multipole loop without
            // uploading the 9 nao² matrices. Otherwise upload the CPU-built ones.
            const bool mp = computed
                ? m_gpu_scf->beginMultipoleComputed()
                : m_gpu_scf->beginMultipole(m_dp_int, m_qp_int, m_basis.ao2at);
            use_gpu_resident = base && mp;
            gpu_multipole = use_gpu_resident;
            if (verb >= 2)
                CurcumaLogger::info(!use_gpu_resident
                    ? "SCF: GPU resident GFN2 begin() failed; running CPU SCF"
                    : (computed
                        ? "SCF: device-resident GFN2 path (CN/S/H0/L on GPU; multipole integrals uploaded)"
                        : "SCF: device-resident GFN2 path (H0/S/L + multipole integrals uploaded)"));
        }
        // Stage 3c: when the device built the integrals, use its Coulomb γ for the
        // host potential (v_sh += γ·q_sh). γ is nsh² (cheap to download once per
        // geometry) and matches the CPU build to ~1e-13; the gemv stays on the
        // host (the rest of v_sh — third-order, D4 — is assembled there too).
        if (use_gpu_resident && computed) {
            Eigen::MatrixXd dev_gamma;
            if (m_gpu_scf->downloadGamma(dev_gamma)
                && dev_gamma.rows() == m_gamma.rows()
                && dev_gamma.cols() == m_gamma.cols())
                m_gamma = dev_gamma;
        }
    }

    // Stage 5 (Part B3/B4, Claude Generated): full device GFN2 potential build.
    // When the backend supports it, the WHOLE per-iteration potential (γ·q_sh +
    // shell third-order + multipole v_dp/v_qp/v_at scalar shift + in-SCF D4) is
    // built on the device, so the SCF loop uploads only the mixed q_sh/dp_at/qp_at
    // (+ the host-built D4 reference weights) instead of v_ao. The geometry-fixed
    // reference data (D4 cache + multipole interaction matrices + third-order
    // hardness) is uploaded once here. Falls back to the host potential build on
    // any failure. GFN2 only (the D4/multipole coupling).
    bool use_device_potential = false;
    if (use_gpu_resident && gpu_multipole && m_method == MethodType::GFN2 && m_mp_initialized
        && m_gpu_scf->supportsDeviceDispersion() && m_gpu_scf->supportsDevicePotential()) {
        // Prime the D4 reference set + device beginDispersion (the m_d4_prepared
        // first-call path); the returned potential is discarded here.
        Vector warm_dedq;
        computeD4PotentialDedq(m_wfn.q_at, warm_dedq);
        // Flatten + upload the multipole interaction matrices (column-major nat²
        // blocks), the on-site XC kernels, and the per-shell third-order hardness
        // Γ_s (= ½·thirdOrderKernelDiag at q_sh=1, GFN2 shell-resolved).
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        std::vector<double> sd(3 * nn), dd(9 * nn), sq(6 * nn), dk(nat), qk(nat), g3(nsh);
        for (int k = 0; k < 3; ++k)
            std::memcpy(sd.data() + static_cast<size_t>(k) * nn, m_mp_amat_sd[k].data(), nn * sizeof(double));
        for (int a = 0; a < 3; ++a)
            for (int bb = 0; bb < 3; ++bb)
                std::memcpy(dd.data() + (static_cast<size_t>(a) * 3 + bb) * nn,
                            m_mp_amat_dd[a][bb].data(), nn * sizeof(double));
        for (int k = 0; k < 6; ++k)
            std::memcpy(sq.data() + static_cast<size_t>(k) * nn, m_mp_amat_sq[k].data(), nn * sizeof(double));
        for (int i = 0; i < nat; ++i) { dk[i] = m_mp_dkernel[i]; qk[i] = m_mp_qkernel[i]; }
        const Vector g3diag = thirdOrderKernelDiag(Vector::Ones(nsh), m_wfn.q_at);  // = 2·Γ_s
        for (int s = 0; s < nsh; ++s) g3[s] = 0.5 * g3diag(s);
        use_device_potential = m_gpu_scf->beginPotential(nat, nsh, sd.data(), dd.data(),
                                                         sq.data(), dk.data(), qk.data(), g3.data());
        if (verb >= 2)
            CurcumaLogger::info(use_device_potential
                ? "SCF: GFN2 potential (Coulomb+third-order+multipole+D4) built on the device"
                : "SCF: device potential build unavailable; host potential build");
    }

    // Stage 6 (S6.5, Claude Generated): fully device-resident SCF loop. When the
    // backend supports it, the WHOLE per-iteration loop body (potential → Fock →
    // eigensolve → occupation → density → charges/moments → SCC energy → Broyden)
    // runs on the device; the host polls only the convergence dq + the 4 energy
    // scalars per step (eps/occ/pop_ao/moments/q_sh stay resident). Upload the
    // q-independent D4 reference tables once (so the device rebuilds W/dWq from the
    // resident charges) and the initial SCC guess + EEQ-guess q_at + reference
    // occupations, then drive m_gpu_scf->residentScfStep below. Falls back to the
    // host-driven device-potential loop on any failure. GFN2 only.
    bool use_resident_loop = false;
    if (use_device_potential && mode == ScfMode::Broyden && m_gpu_scf->supportsResidentLoop()) {
        std::vector<double> rcn, rgi, rzeff, rrefcn, rrefcovcn, rrefq;
        std::vector<int> rnref;
        exportD4RefWDeviceData(rcn, rgi, rzeff, rrefcn, rrefcovcn, rrefq, rnref);
        const int nocc_pairs = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
        const int b_maxhist = (m_diis_subspace > 2) ? 20 : m_diis_subspace;
        use_resident_loop = !rcn.empty()
            && m_gpu_scf->beginDispersionWeights(rcn, rgi, rzeff, rrefcn, rrefcovcn, rrefq, rnref)
            && m_gpu_scf->beginResidentLoop(m_wfn.q_sh, m_wfn.dp_at, m_wfn.qp_at, m_wfn.q_at,
                                            m_wfn.n0_sh, m_wfn.n0_at, m_electronic_temp,
                                            m_wfn.nocc, nocc_pairs, damp, b_maxhist, 0.01);
        if (verb >= 2)
            CurcumaLogger::info(use_resident_loop
                ? "SCF: fully device-resident loop (host polls only dq + energy scalars)"
                : "SCF: resident-loop setup unavailable; host-driven device-potential loop");
    }

    // AP1 (Claude Generated): GPU partial diagonalisation window. The device density
    // only needs the occupied(+buffer) eigenvectors, so the resident eigensolve can
    // compute just the lowest gpu_n_eig pairs (cusolverDnDsyevdx) instead of the full
    // nao spectrum. OPT-IN (m_gpu_partial_diag, -gpu_partial_diag): measured
    // net-neutral on an RTX 5080 (the tridiagonalization, not the eigenvector count,
    // dominates the dense eigensolve), so the default GPU path stays the full solve.
    // Disabled when virtuals are genuinely required: the Mulliken-CPSCF D4 charge
    // response (d4_charge_source=mulliken) and the verbosity>=3 orbital listing both
    // read virtual orbitals. The buffer covers the LUMO plus a few-kT Fermi tail; if
    // the occupied tail ever reaches the window edge (tiny-gap/metallic), the loop
    // widens to the full solve once so the 1e-8 energy is never at risk. Disabled on
    // the device-potential path (its re-solve would need the host D4 weights again).
    int gpu_n_eig = 0;  // 0 = full spectrum
    bool gpu_partial_diag = m_gpu_partial_diag && use_gpu_resident && !use_device_potential
        && (m_d4_charge_source != "mulliken") && (verb < 3);
    if (gpu_partial_diag) {
        const int nocc_orbs = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
        const int buffer = std::max(8, (m_basis.nao + 19) / 20);  // ~5% of nao, >= 8
        gpu_n_eig = std::min(m_basis.nao, nocc_orbs + buffer);
        if (gpu_n_eig >= m_basis.nao) { gpu_n_eig = 0; gpu_partial_diag = false; }
    }
    if (gpu_partial_diag && verb >= 2)
        CurcumaLogger::info("GPU partial diagonalisation: lowest " + std::to_string(gpu_n_eig)
            + " of " + std::to_string(m_basis.nao) + " eigenpairs (occupied " +
            std::to_string(static_cast<int>(std::floor(m_wfn.nocc / 2.0))) + " + buffer)");

    const auto t_scf_start = clock::now();
    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        const auto t_iter0 = clock::now();

        // Stage 6 (S6.5): fully device-resident step. One backend call runs the
        // whole iteration on the GPU; only the convergence dq + the 4 SCC energy
        // scalars come back. The host keeps the iteration counter, the mixed-
        // precision (FP32→FP64) decision, the convergence test, and the verbosity
        // line — the rest (incl. the Broyden mix) is device-resident. Claude Generated.
        if (use_resident_loop) {
            m_eig_fp32 = m_scf_mixed_precision && (dq_prev > m_scf_fp32_threshold);
            double dq = 0.0, eb = 0.0, ecoul = 0.0, ethird = 0.0, emp = 0.0;
            if (!m_gpu_scf->residentScfStep(m_eig_fp32, dq, eb, ecoul, ethird, emp)) {
                CurcumaLogger::warn("XTB::Calculation: GPU resident SCF step failed at iteration "
                                    + std::to_string(iter));
                m_scf_converged = false; m_scf_iterations = iter; return m_E_total;
            }
            m_E_electronic = eb; m_E_coulomb_shell = ecoul;
            m_E_third_order = ethird; m_E_multipole = emp;
            const double e_scc = eb + ecoul + ethird + emp;
            const double de = (iter > 0) ? std::fabs(e_scc - e_total_old) : 0.0;
            dq_prev = dq;
            if (verb >= scf_min) {
                const double t_iter_ms = ms(t_iter0, clock::now());
                CurcumaLogger::result(fmt::format("  {:4d}   {:18.10f}   {:10.2e}   {:10.2e}   {:6.1f}",
                                                  iter, e_scc, de, dq, t_iter_ms));
            }
            if (iter > 0 && !m_eig_fp32 && dq < thresh && de < thresh * 100.0) {
                m_scf_converged = true;
                e_total_old = e_scc;
                break;
            }
            e_total_old = e_scc;
            continue;
        }

        // Broyden mixes the SCC charge vector: capture the input x_in (current
        // m_wfn charges) before the potential is built and the populations are
        // overwritten by the diagonalisation.
        Vector x_in;
        if (mode == ScfMode::Broyden)
            x_in = packSCC();

        // Reset and build potentials. On the device-potential path (Stage 5 B3/B4)
        // the whole potential is built inside the GPU solve below, so the host
        // build is skipped entirely (the per-iteration energies are recomputed
        // from the charges, not from m_pot — see energyCoulombShell/ThirdOrder/
        // Multipole below). m_pot is rebuilt once post-SCF for the gradient.
        if (!use_device_potential) {
            m_pot.reset();

            // Isotropic Coulomb: v_sh += γ * q_sh
            addCoulombShellPotential(m_pot);

            // Third-order: add to v_at/v_sh
            addThirdOrderPotential(m_pot);

            // GFN2 multipole
            if (m_method == MethodType::GFN2) {
                addMultipolePotential(m_pot);
                // Self-consistent D4: exact per-reference dE_D4/dq into the atom potential.
                const auto t_disp0 = clock::now();
                addDispersionPotential(m_pot);
                acc_disp += ms(t_disp0, clock::now());
            }
        }
        const auto t_pot = clock::now();

        // Build Fock, diagonalise, build the density and update populations. The
        // device-resident GFN1 path (Claude Generated, GPU port Stage 2) fuses
        // the Fock build + eigensolve + density + Mulliken-AO onto the GPU with
        // H0/S/L resident — only length-nao vectors cross the bus — while the CPU
        // path below is unchanged. The shared scaffolding (mixing, energies,
        // convergence) follows both branches.
        double gpu_band = 0.0;
        std::chrono::steady_clock::time_point t_fock, t_solve, t_mull;
        if (use_gpu_resident) {
            // Expand shell+atom potential to AO resolution (host; mirrors
            // expand_potential in xtb_scf.cpp). v_ao(μ)=v_sh(sh)+v_at(at) — this
            // isotropic part carries Coulomb + third-order (v_sh) and, for GFN2,
            // the multipole scalar shift + in-SCF D4 (v_at), all built on the host
            // above. The GFN2 anisotropic dipole/quadrupole potential (v_dp/v_qp)
            // is applied in the device Fock build via solveMultipole.
            // Mixed precision (GPU): solve in FP32 while far from convergence
            // (max|dq| above the threshold; iter 0 starts FP32 via dq_prev=1e30),
            // reverting to FP64 near convergence so the converged energy is FP64.
            m_eig_fp32 = m_scf_mixed_precision && (dq_prev > m_scf_fp32_threshold);
            bool solve_ok;
            Eigen::VectorXd v_ao;  // host-expanded potential (non-device-potential path)
            if (use_device_potential) {
                // Stage 5 (B3/B4): the device builds v_sh/v_at/v_dp/v_qp from the
                // mixed q_sh/dp_at/qp_at + the host-built D4 reference weights, then
                // the Fock + eigensolve — no v_ao upload. buildRefWFlat reuses the
                // validated host per-atom CN-Gaussian × zeta weights at the current
                // (mixed) charges.
                const auto t_disp0 = clock::now();
                std::vector<double> W, dWq;
                m_d4_generator->buildRefWFlat(m_atoms, m_wfn.q_at, W, dWq);
                acc_disp += ms(t_disp0, clock::now());
                t_fock = clock::now();
                solve_ok = m_gpu_scf->solvePotential(m_wfn.q_sh, m_wfn.dp_at, m_wfn.qp_at,
                                                     W, dWq, m_wfn.eps, m_eig_fp32, gpu_n_eig);
            } else {
                // Expand shell+atom potential to AO resolution (host; mirrors
                // expand_potential in xtb_scf.cpp). v_ao(μ)=v_sh(sh)+v_at(at) carries
                // Coulomb + third-order (v_sh) and, for GFN2, the multipole scalar
                // shift + in-SCF D4 (v_at), all built on the host above. The GFN2
                // anisotropic v_dp/v_qp is applied in the device Fock via solveMultipole.
                v_ao.resize(m_basis.nao);
                for (int ao = 0; ao < m_basis.nao; ++ao)
                    v_ao(ao) = m_pot.v_sh(m_basis.ao2sh[ao]) + m_pot.v_at(m_basis.ao2at[ao]);
                t_fock = clock::now();   // Fock build is fused into solve() on the device
                // Device: F = H0 − ½·S·(v_ao⊕v_ao) [+ GFN2 multipole]; solve F C = S C ε.
                // AP1: gpu_n_eig>0 solves only the lowest occupied(+buffer) eigenpairs.
                solve_ok = gpu_multipole
                    ? m_gpu_scf->solveMultipole(v_ao, m_pot.v_dp, m_pot.v_qp, m_wfn.eps, m_eig_fp32, gpu_n_eig)
                    : m_gpu_scf->solve(v_ao, m_wfn.eps, m_eig_fp32, gpu_n_eig);
            }
            if (!solve_ok) {
                CurcumaLogger::warn("XTB::Calculation: GPU resident solve failed at iteration "
                                    + std::to_string(iter));
                m_scf_converged = false; m_scf_iterations = iter; return m_E_total;
            }
            // Occupations on the host (shared Fermi/integer logic with the CPU path).
            Eigen::VectorXd occ; int ncol = 0;
            occupationsFromEps(m_wfn.eps, occ, ncol);
            // AP1 validation: ncol reaching the partial window edge means the occupied
            // tail spilled past gpu_n_eig (tiny-gap / metallic). Widen to the full
            // spectrum once and re-occupy so the density/energy stay exact; the rest of
            // the SCF then runs the full solve.
            if (gpu_partial_diag && gpu_n_eig > 0 && ncol >= gpu_n_eig) {
                if (verb >= 1)
                    CurcumaLogger::warn("GPU partial diagonalisation window too small (occupied "
                        "tail reached eigenpair " + std::to_string(gpu_n_eig)
                        + "); switching to the full spectrum for the rest of the SCF.");
                gpu_partial_diag = false; gpu_n_eig = 0;
                const bool re_ok = gpu_multipole
                    ? m_gpu_scf->solveMultipole(v_ao, m_pot.v_dp, m_pot.v_qp, m_wfn.eps, m_eig_fp32, gpu_n_eig)
                    : m_gpu_scf->solve(v_ao, m_wfn.eps, m_eig_fp32, gpu_n_eig);
                if (!re_ok) {
                    CurcumaLogger::warn("XTB::Calculation: GPU resident full re-solve failed at "
                                        "iteration " + std::to_string(iter));
                    m_scf_converged = false; m_scf_iterations = iter; return m_E_total;
                }
                occupationsFromEps(m_wfn.eps, occ, ncol);
            }
            // Device: P = C·diag(occ)·Cᵀ; pop_ao(μ)=Σ_ν P_μν·S_μν; band=Σ P⊙H0.
            Eigen::VectorXd pop_ao;
            if (!m_gpu_scf->density(occ, ncol, pop_ao, gpu_band)) {
                CurcumaLogger::warn("XTB::Calculation: GPU resident density failed at iteration "
                                    + std::to_string(iter));
                m_scf_converged = false; m_scf_iterations = iter; return m_E_total;
            }
            t_solve = clock::now();
            // Broyden never mixes the density, so there is no P damping here.
            updatePopulationsFromPopAo(pop_ao);   // isotropic q_sh / q_at
            // GFN2: atomic dipole/quadrupole moments from the resident density
            // (the multipole block of updatePopulations), into m_wfn.dp_at/qp_at
            // so the shared Broyden pack/mix/unpack and energyMultipole work.
            if (gpu_multipole && !m_gpu_scf->multipoleMoments(m_wfn.dp_at, m_wfn.qp_at)) {
                CurcumaLogger::warn("XTB::Calculation: GPU resident multipole moments failed at iteration "
                                    + std::to_string(iter));
                m_scf_converged = false; m_scf_iterations = iter; return m_E_total;
            }
            q_sh_new = m_wfn.q_sh;
            t_mull = clock::now();
        } else {
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

            t_fock = clock::now();

            // Mixed precision (opt-in, MKL path): solve in FP32 while far from convergence
            // (previous max|dq| above the threshold; iter 0 starts in FP32 via dq_prev=1e30),
            // reverting to FP64 near convergence so the converged energy is FP64. Claude Generated.
            m_eig_fp32 = m_scf_mixed_precision && (dq_prev > m_scf_fp32_threshold);

            // Diagonalize. Let MKL thread the eigensolve (dsygst/dsyevd/dtrsm) for a
            // single large molecule; the surrounding MklSerialScope keeps MKL serial
            // everywhere else so the hand-threaded regions are not oversubscribed.
            bool eig_ok;
            {
                MklThreadScope eig_scope(eig_threads);
                eig_ok = solveEigen(F, m_S);
            }
            if (!eig_ok) {
                CurcumaLogger::warn("XTB::Calculation: eigen solve failed at iteration " + std::to_string(iter));
                m_scf_converged = false;
                m_scf_iterations = iter;
                return m_E_total;
            }
            t_solve = clock::now();

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
            t_mull = clock::now();
        }

        // Energies for this iteration
        m_E_coulomb_shell = energyCoulombShell();
        m_E_third_order   = energyThirdOrder();
        m_E_multipole     = energyMultipole();

        // Band energy: Tr(P · H0). Device-resident path returns it from the GPU
        // (P is not on the host during the loop); CPU path sums P⊙H0 directly.
        m_E_electronic = use_gpu_resident ? gpu_band
                                          : (m_wfn.P.cwiseProduct(m_H0)).sum();

        // Total SCC = band + coulomb + third-order + multipole
        const double e_scc = m_E_electronic + m_E_coulomb_shell
                           + m_E_third_order + m_E_multipole;

        // Accumulate per-phase timings (see breakdown at verbosity >= 3).
        acc_pot    += ms(t_iter0, t_pot);
        acc_fock   += ms(t_pot,   t_fock);
        acc_solve  += ms(t_fock,  t_solve);
        acc_mull   += ms(t_solve, t_mull);
        acc_energy += ms(t_mull,  clock::now());

        // Per-iteration diagnostics
        const double dq = (q_sh_new - q_sh_old).cwiseAbs().maxCoeff();
        dq_prev = dq;   // drives the LevelShift fade-out on the next iteration
        const double de = (iter > 0) ? std::fabs(e_scc - e_total_old) : 0.0;
        const double t_iter_ms = ms(t_iter0, clock::now());
        if (verb >= scf_min) {
            std::string line = fmt::format("  {:4d}   {:18.10f}   {:10.2e}   {:10.2e}   {:6.1f}",
                                           iter, e_scc, de, dq, t_iter_ms);
            if (verb >= 3)
                line += fmt::format("   [DIIS n={} err={:.2e}]", diis.size(), diis.lastErrorNorm());
            CurcumaLogger::result(line);
        }

        // Check convergence (after first iteration). On convergence we break
        // here, so m_wfn keeps the consistent output charges x_out (the Broyden
        // mix below is skipped — no re-mixing of the converged state).
        // Never accept convergence on a mixed-precision (FP32) step — otherwise a
        // fast-converging molecule (e.g. LiH) can declare convergence on an
        // FP32-accurate density (~1e-7 error). Requiring !m_eig_fp32 forces a final
        // FP64 polish: once max|dq| drops below the FP32 threshold the next step is
        // FP64, and that step is the one allowed to converge. No-op when mixed
        // precision is off (m_eig_fp32 always false). Claude Generated.
        if (iter > 0 && !m_eig_fp32) {
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

    if (verb >= scf_min) {
        if (m_scf_converged)
            CurcumaLogger::success_fmt("SCF converged in {} iterations ({:.1f} ms)",
                                       m_scf_iterations, ms(t_scf_start, t_scf_end));
        else
            CurcumaLogger::warn_fmt("SCF NOT converged after {} iterations ({:.1f} ms)",
                                    m_scf_iterations, ms(t_scf_start, t_scf_end));
    }

    // Stage 6 (S6.5): the fused device loop kept the converged charges/moments
    // resident — download them once into the host wavefunction so the post-SCF
    // potential rebuild + energies + gradient below run unchanged. Claude Generated.
    if (use_resident_loop && m_gpu_scf) {
        m_wfn.q_sh.resize(m_basis.nsh);
        m_wfn.q_at.resize(m_atomcount);
        m_wfn.dp_at.resize(3, m_atomcount);
        m_wfn.qp_at.resize(6, m_atomcount);
        m_wfn.eps.resize(m_basis.nao);
        m_gpu_scf->residentLoopCharges(m_wfn.q_sh, m_wfn.q_at, m_wfn.dp_at, m_wfn.qp_at, m_wfn.eps);
    }

    // Device-resident path: download the converged density and MO coefficients
    // once, so the post-SCF energies (band Tr(P·H0)) and the (still CPU) gradient
    // read them from the host. Claude Generated, GPU port Stage 2.
    if (use_gpu_resident && m_gpu_scf)
        m_gpu_scf->finalize(m_wfn.P, m_wfn.C);

    // Persist converged Fock matrix for gradient / debug. On the device-potential
    // path (Stage 5 B3/B4) the host m_pot was skipped during the loop, so rebuild
    // it once at the converged charges before forming m_F (the gradient block below
    // rebuilds it again, but m_F must be consistent here too).
    if (use_device_potential) {
        m_pot.reset();
        addCoulombShellPotential(m_pot);
        addThirdOrderPotential(m_pot);
        addMultipolePotential(m_pot);
        addDispersionPotential(m_pot);
    }
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

        // Stage 4 (Claude Generated): native nuclear gradient on the device when the
        // SCF ran device-resident. GFN1: repulsion/H0-Pulay/Coulomb (1/2/3) on the
        // GPU; GFN2 additionally the multipole-integral Pulay on the GPU and the
        // direct multipole interaction (section 5) on the host. Dispersion (3b) +
        // CN chain-rule (4) on the host. Falls back to the full CPU gradient if the
        // backend has no gradient.
        bool gpu_grad = false;
        if (use_gpu_resident && m_gpu_scf && m_gpu_scf->supportsGradient())
            gpu_grad = calculateGradientGpu();
        if (!gpu_grad)
            calculateGradient();   // fills m_gradient in Eh/Bohr
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

    if (verb >= 3) {
        const int it = m_scf_iterations > 0 ? m_scf_iterations : 1;
        const double solve_sum = m_t_xfx + m_t_diag + m_t_back + m_t_dens;
        CurcumaLogger::info("SCF per-phase breakdown (summed over iterations):");
        CurcumaLogger::info_fmt("  potential build : {:8.2f} ms ({:5.2f}/it)", acc_pot,    acc_pot / it);
        if (m_method == MethodType::GFN2)
            CurcumaLogger::info_fmt("    - of which D4 : {:8.2f} ms ({:5.2f}/it)", acc_disp, acc_disp / it);
        CurcumaLogger::info_fmt("  build Fock      : {:8.2f} ms ({:5.2f}/it)", acc_fock,   acc_fock / it);
        CurcumaLogger::info_fmt("  solve eigen     : {:8.2f} ms ({:5.2f}/it)", acc_solve,  acc_solve / it);
        CurcumaLogger::info_fmt("    - reduce      : {:8.2f} ms ({:5.2f}/it)", m_t_xfx,    m_t_xfx / it);
        CurcumaLogger::info_fmt("    - dsyevd      : {:8.2f} ms ({:5.2f}/it)", m_t_diag,   m_t_diag / it);
        CurcumaLogger::info_fmt("    - back-transf : {:8.2f} ms ({:5.2f}/it)", m_t_back,   m_t_back / it);
        CurcumaLogger::info_fmt("    - density P   : {:8.2f} ms ({:5.2f}/it)", m_t_dens,   m_t_dens / it);
        CurcumaLogger::info_fmt("    (solve detail sum = {:.2f} ms, untimed = {:.2f} ms)",
                                solve_sum, acc_solve - solve_sum);
        CurcumaLogger::info_fmt("  populations     : {:8.2f} ms ({:5.2f}/it)", acc_mull,   acc_mull / it);
        CurcumaLogger::info_fmt("  energy/mix      : {:8.2f} ms ({:5.2f}/it)", acc_energy, acc_energy / it);
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
    m_gpu_basis_dirty = true;   // basis changed → re-upload to the device once
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

    // 5. Reset wavefunction populations (Mulliken charges depend on geometry).
    //    Save the converged shell charges before clearing so Calculation() can
    //    use them as a warm-start guess for the next SCF (Claude Generated).
    // Save converged shell charges for warm-start only when the previous SCF
    // actually converged; zero-initialized charges (from InitialiseMolecule)
    // are a worse starting guess than EEQ and must not be saved. Claude Generated.
    const bool had_converged = m_scf_converged;
    if (m_warmstart && had_converged && m_wfn.q_sh.size() == m_basis.nsh) {
        m_warmstart_q_sh = m_wfn.q_sh;
        // GFN2: carry the converged atomic multipoles too, so Calculation() can restore
        // the FULL SCC vector. Without this the multipoles re-converge from zero each
        // step and the opt/md SCF plateaus (see m_warmstart_dp_at doc). Claude Generated.
        if (m_method == MethodType::GFN2
            && m_wfn.dp_at.cols() == m_basis.nat && m_wfn.qp_at.cols() == m_basis.nat) {
            m_warmstart_dp_at = m_wfn.dp_at;
            m_warmstart_qp_at = m_wfn.qp_at;
        }
        // Multi-step SCC extrapolation (opt-in): record the converged SCC vector so
        // Calculation() can predict the next step from several past steps, not just
        // this one. Pushed before m_wfn is zeroed below. Trimmed to the order the
        // predictor needs. Claude Generated.
        if (m_scf_extrapolation != "none") {
            Eigen::VectorXd state = packSccState();
            if (!m_scf_history.empty() && m_scf_history.front().size() != state.size())
                m_scf_history.clear(); // geometry/basis changed -> stale history
            m_scf_history.push_front(std::move(state));
            const int maxlen = (m_scf_extrapolation == "aspc")
                                   ? std::max(m_scf_extrap_order + 2, 2)
                                   : std::min(2 * m_scf_extrap_order + 2, 10);
            while (static_cast<int>(m_scf_history.size()) > maxlen)
                m_scf_history.pop_back();
        }
    }
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
        // Claude Generated (May 2026): analytical gradient replaces central
        // finite differences. The direct geometry gradient + CN chain rule
        // are computed in one pass through the D3 pair list.
        if (!m_d3_generator) {
            m_d3_generator = std::make_unique<::D3ParameterGenerator>(
                ::D3ParameterGenerator::createForGFN1());
        }
        ::D3ParameterGenerator& d3 = *m_d3_generator;
        // Lightweight prepare (CN + Gaussian weights only) — getEnergyAndGradient
        // computes C6/C8/energy/gradient directly from the references, so the
        // O(N²) JSON pair list that GenerateParameters builds is pure overhead
        // here (≈0.5 s on complex/231). Claude Generated.
        d3.prepareForEnergyGradient(m_atoms, m_geometry);

        if (!need_gradient) {
            m_disp_gradient = Matrix::Zero(m_atomcount, 3);
            m_disp_dEdcn = Vector();
            m_disp_gradient_valid = false;
            Matrix g_unused; Vector cn_unused;
            return d3.getEnergyAndGradient(/*need_gradient=*/false, g_unused, cn_unused);
        }

        // Analytical energy + gradient + dEdcn
        m_disp_gradient = Matrix::Zero(m_atomcount, 3);
        m_disp_dEdcn = Vector::Zero(m_atomcount);
        const double e_disp = d3.getEnergyAndGradient(
            /*need_gradient=*/true, m_disp_gradient, m_disp_dEdcn);

        // CN chain rule: fold Σ_A dE_D3/dCN_A · ∂CN_A/∂R into the cached
        // gradient using the D3 exponential counting function's derivative.
        // This keeps the dispersion gradient self-consistent with the D3 CN.
        // m_disp_dEdcn is then cleared so the GFN1 Hamiltonian CN loop in
        // calculateGradient() does not double-distribute it.
        if (m_disp_dEdcn.size() == m_atomcount) {
            // addD3CNGradient produces Eh/Å (geometry is in Å); scale by au
            // to convert to Eh/Bohr so m_disp_gradient matches the direct term.
            CNCalculator::addD3CNGradient(
                m_atoms, m_geometry, m_disp_dEdcn, m_disp_gradient,
                /*k1=*/16.0, /*k2=*/4.0/3.0, /*distance_unit_to_bohr=*/au);
            m_disp_dEdcn = Vector();
        }

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
        // Native GFN2 consumes the C6 reference cache through D4Evaluator, not the
        // JSON pair list — skip building it (json-per-pair + OpenMP every geometry).
        m_d4_generator->setBuildPairLists(false);
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
    // m_geometry is in Angstrom; D4 expects Bohr.
    Matrix geom_bohr = m_geometry * AA_TO_AU;
    // Claude Generated (2026-05): reuse the geometry-prepared D4 reference set from the
    // SCF loop (addDispersionPotential). Only the converged Mulliken charges are pushed
    // here. If the SCF never prepared it (energy-only path with no potential build),
    // run the full GenerateParameters once. updateCNValuesForGradient (dc6/dCN) is then
    // refreshed inside computeEnergyAndGradient(with_gradient=true) below.
    m_d4_generator->setTopologyCharges(m_wfn.q_at);
    if (!m_d4_prepared) {
        m_d4_generator->setUseD4SingleShotEEQ(false);
        // GFN2: use the dftd4 EN-weighted covalent CN (matches tblite's
        // get_coordination_number(rcov, en)) instead of the GFN-FF log-capped erf-CN.
        m_d4_generator->setD4CovalentCN(true);
        m_d4_generator->GenerateParameters(m_atoms, geom_bohr);
        ++m_d4_genparams_calls;
        m_d4_prepared = true;
    }

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
    using d4clk = std::chrono::steady_clock;
    auto d4ms = [](d4clk::time_point a, d4clk::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    const auto td0 = d4clk::now();
    // WP2-ext: let the D4 pair loop fan out over the same gated intra-thread budget.
    m_d4_evaluator->setThreads(effectiveIntraThreads(m_atomcount));
    double E = m_d4_evaluator->computeEnergyAndGradient(
        m_atoms, geom_bohr,
        /*with_gradient=*/true,
        m_disp_gradient, m_disp_dEdcn, m_disp_dEdq, want_dEdq);
    const auto td1 = d4clk::now();

    // ATM three-body term (GFN2 s9=5.0). Charge-independent (dftd4 evaluates it at
    // the q=0 reference C6 — matches tblite get_dispersion_nonsc), so no q-response.
    // ACCUMULATES into m_disp_gradient and m_disp_dEdcn (the latter folded together
    // with the 2-body dEdcn via the D4 covalent CN below). Closes the triose C-path
    // remainder (see docs/GFN2_D4_STATUS.md). Claude Generated 2026-05-29.
    E += m_d4_evaluator->computeATM(
        m_atoms, geom_bohr, /*with_gradient=*/true, m_disp_gradient, m_disp_dEdcn);
    const auto td2 = d4clk::now();

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
    const auto td3 = d4clk::now();
    if (m_disp_dEdq.size() == m_atomcount && !m_disp_audit_mode) {
        // q-response source (m_d4_charge_source, default "eeq"):
        //   "eeq"      — analytic ∂q/∂x from the single-shot dftd4 EEQ model
        //                (one LU solve + adjoint contraction, ~ms). dftd4-conform
        //                and the validated default. ~100x cheaper than the CPSCF.
        //   "mulliken" — exact ∂q_Mulliken/∂x via the GFN2 CPSCF/Z-vector solve
        //                (computeMullikenChargeResponse). Self-consistent with the
        //                Mulliken-charge D4 energy but expensive (dominates the
        //                post-SCF time on large systems: ~88% of D4 on complex).
        // The energy weighting stays on the SCF Mulliken charges either way, so
        // this choice only affects the gradient, never the energy.
        if (m_d4_charge_source == "mulliken") {
            computeMullikenChargeResponse(m_disp_dEdq, m_disp_gradient);
        } else if (m_gpu_scf && m_gpu_scf->supportsDeviceEeq()) {
            // Stage 5 (Part A): device EEQ charges + adjoint dq/dx response. The
            // GPU writes the per-atom response into a scratch N×3 buffer; we add it
            // to m_disp_gradient. Host D4ChargeModel is the fallback on failure.
            const int nat = m_atomcount;
            std::vector<double> chi, gam, alp, cnf, rcov;
            curcuma::dispersion::D4ChargeModel::resolveParams(m_atoms, chi, gam, alp, cnf, rcov);
            std::vector<double> xyz(3 * nat), qscratch(nat, 0.0), dedq(nat, 0.0),
                grad(3 * nat, 0.0);
            for (int i = 0; i < nat; ++i) {
                xyz[3 * i + 0] = geom_bohr(i, 0);
                xyz[3 * i + 1] = geom_bohr(i, 1);
                xyz[3 * i + 2] = geom_bohr(i, 2);
                dedq[i] = m_disp_dEdq(i);
            }
            const bool ok =
                m_gpu_scf->eeqCharges(nat, xyz.data(), chi.data(), gam.data(), alp.data(),
                                      cnf.data(), rcov.data(), static_cast<double>(m_charge),
                                      qscratch.data())
                && m_gpu_scf->eeqChargeResponse(nat, dedq.data(), grad.data());
            if (ok) {
                for (int a = 0; a < nat; ++a) {
                    m_disp_gradient(a, 0) += grad[3 * a + 0];
                    m_disp_gradient(a, 1) += grad[3 * a + 1];
                    m_disp_gradient(a, 2) += grad[3 * a + 2];
                }
            } else {
                curcuma::dispersion::D4ChargeModel eeq;
                eeq.computeCharges(m_atoms, geom_bohr, static_cast<double>(m_charge));
                if (eeq.valid())
                    eeq.addChargeResponseGradient(m_disp_dEdq, m_disp_gradient);
            }
        } else {
            curcuma::dispersion::D4ChargeModel eeq;
            eeq.computeCharges(m_atoms, geom_bohr, static_cast<double>(m_charge));
            if (eeq.valid())
                eeq.addChargeResponseGradient(m_disp_dEdq, m_disp_gradient);
        }
    }
    const auto td4 = d4clk::now();
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info_fmt("  D4 energy+grad   : {:8.2f} ms", d4ms(td0, td1));
        CurcumaLogger::info_fmt("  D4 ATM 3-body    : {:8.2f} ms", d4ms(td1, td2));
        CurcumaLogger::info_fmt("  D4 CN-grad fold  : {:8.2f} ms", d4ms(td2, td3));
        CurcumaLogger::info_fmt("  D4 q-resp CPSCF  : {:8.2f} ms", d4ms(td3, td4));
    }
    m_disp_gradient_valid = true;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("dispersion_d4_native", fmt::format("{:.6f} Eh", E));
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("d4_genparams_calls_this_geometry",
                             std::to_string(m_d4_genparams_calls));
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
bool XTB::computeD4PotentialDedq(const Vector& q_at, Vector& dEdq) const
{
    if (m_method != MethodType::GFN2) return false;

    // Lazy init shared with calcDispersionEnergy (idempotent).
    if (!m_d4_generator) {
        nlohmann::json d4_config;
        d4_config["d4_s6"] = 1.0; d4_config["d4_s8"] = 2.7;
        d4_config["d4_a1"] = 0.52; d4_config["d4_a2"] = 5.0; d4_config["d4_alp"] = 16.0;
        ConfigManager cfg("d4param", d4_config);
        m_d4_generator = std::make_unique<::D4ParameterGenerator>(cfg);
        // Native GFN2 consumes the C6 reference cache through D4Evaluator, not the
        // JSON pair list — skip building it (json-per-pair + OpenMP every geometry).
        m_d4_generator->setBuildPairLists(false);
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
    //
    // Claude Generated (2026-05): the heavy D4 reference set (CN + Gaussian weights +
    // C6 reference cache) is geometry-dependent and identical across SCF iterations;
    // only the SCF charges change. Run the full GenerateParameters once per geometry
    // (m_d4_prepared); on subsequent iterations only push the new SCF charges and
    // re-evaluate dE_D4/dq against the cached reference set. This removes the per-SCF
    // CN + Gaussian-weight regeneration while keeping the result bit-identical (the
    // weights/CN are geometry-fixed; only the charges enter the per-reference C6
    // scaling). with_gradient stays true because D4Evaluator only fills dEdq when a
    // gradient is requested (do_dEdq = with_gradient && with_dEdq); the scratch
    // Cartesian gradient is discarded here and reassembled once post-SCF.
    const Matrix geom_bohr = m_geometry * AA_TO_AU;
    m_d4_generator->setTopologyCharges(q_at);
    // Stage 5 (Part B2): on the device-resident path run the per-iteration O(N²)
    // dE_D4/dq pair loop on the GPU (the host still builds the per-atom reference
    // weights once per iteration; the device does the contraction + BJ damping).
    const bool use_device_disp = (m_gpu_scf && m_gpu_scf->supportsDeviceDispersion());
    if (!m_d4_prepared) {
        m_d4_generator->setUseD4SingleShotEEQ(false);
        m_d4_generator->setD4CovalentCN(true);  // GFN2 dftd4 EN-weighted covalent CN
        m_d4_generator->GenerateParameters(m_atoms, geom_bohr);
        ++m_d4_genparams_calls;
        m_d4_prepared = true;
        // WP2: new geometry → the evaluator's cached pair list is stale (rebuilt on the
        // next computeEnergyAndGradient). Per-iter calls within this geometry reuse it.
        m_d4_evaluator->invalidatePairCache();
        if (use_device_disp) {
            // Upload the geometry-fixed D4 reference data once per geometry. The
            // BJ params (1, 2.7, 0.52, 5.0) match the evaluator init above; the
            // effective pair cutoff is 50 Bohr (GFNFFDispersion::r_cut).
            const int nat = m_atomcount;
            const std::vector<int>& refn = m_d4_generator->getRefN();
            std::vector<int> Zv(nat), nrefv(nat);
            std::vector<double> sqrtv(nat), xyzv(3 * nat);
            for (int a = 0; a < nat; ++a) {
                const int Z = m_atoms[a];
                Zv[a]    = Z;
                nrefv[a] = (Z >= 1 && (Z - 1) < static_cast<int>(refn.size())) ? refn[Z - 1] : 0;
                sqrtv[a] = m_d4_generator->getSqrtZr4r2(Z);
                xyzv[3 * a + 0] = geom_bohr(a, 0);
                xyzv[3 * a + 1] = geom_bohr(a, 1);
                xyzv[3 * a + 2] = geom_bohr(a, 2);
            }
            const std::vector<double>& c6 = m_d4_generator->getC6FlatCache();
            m_gpu_scf->beginDispersion(nat, Zv.data(), sqrtv.data(), nrefv.data(),
                                       xyzv.data(), c6.data(),
                                       static_cast<int>(c6.size()),
                                       1.0, 2.7, 0.52, 5.0, 50.0);
        }
    }

    bool device_ok = false;
    if (use_device_disp) {
        // Host builds the per-atom reference weights (CN-Gaussian × zeta) once per
        // iteration; the device runs the O(N²) contraction + BJ disp_sum → dEdq.
        std::vector<double> W, dWq;
        m_d4_generator->buildRefWFlat(m_atoms, q_at, W, dWq);
        std::vector<double> dedq(m_atomcount, 0.0);
        if (m_gpu_scf->dispersionDedq(m_atomcount, W.data(), dWq.data(), dedq.data())) {
            dEdq = Vector::Zero(m_atomcount);
            for (int A = 0; A < m_atomcount; ++A) dEdq(A) = dedq[A];
            device_ok = true;
        }
    }
    if (!device_ok) {
        Matrix scratch_grad; Vector scratch_dEdcn;
        // WP2-ext: thread the per-SCF D4 pair loop on the same gated budget.
        m_d4_evaluator->setThreads(effectiveIntraThreads(m_atomcount));
        m_d4_evaluator->computeEnergyAndGradient(
            m_atoms, geom_bohr, /*with_gradient=*/true,
            scratch_grad, scratch_dEdcn, dEdq, /*with_dEdq=*/true);
    }

    return dEdq.size() == m_atomcount;
}

/* ------------------------------------------------------------------------- *
 *  Stage 6 (S6.2b) validation passthroughs (Claude Generated). The host
 *  reference W/dWq and the q-independent per-atom reference tables the device
 *  kernel k_d4_build_refw rebuilds them from. Require a prior GFN2 Calculation
 *  (m_d4_generator prepared by computeD4PotentialDedq / calcDispersionEnergy).
 * ------------------------------------------------------------------------- */
void XTB::buildD4RefWFlat(const Vector& q, std::vector<double>& W,
                          std::vector<double>& dWq) const
{
    W.clear(); dWq.clear();
    if (m_d4_generator) m_d4_generator->buildRefWFlat(m_atoms, q, W, dWq);
}

void XTB::exportD4RefWDeviceData(std::vector<double>& cn, std::vector<double>& gi,
                                 std::vector<double>& zeff, std::vector<double>& refcn,
                                 std::vector<double>& refcovcn, std::vector<double>& refq,
                                 std::vector<int>& nref) const
{
    cn.clear(); gi.clear(); zeff.clear(); refcn.clear(); refcovcn.clear(); refq.clear(); nref.clear();
    if (m_d4_generator)
        m_d4_generator->exportRefWDeviceData(m_atoms, cn, gi, zeff, refcn, refcovcn, refq, nref);
}

/* ------------------------------------------------------------------------- *
 *  Self-consistent D4 dispersion potential (GFN2, AP6b) — thin wrapper that
 *  folds the per-atom dE_D4/dq (computeD4PotentialDedq, device or host) into the
 *  Fock atom-potential. Claude Generated.
 * ------------------------------------------------------------------------- */
void XTB::addDispersionPotential(Potential& pot) const
{
    if (m_method != MethodType::GFN2) return;
    Vector dEdq;
    if (!computeD4PotentialDedq(m_wfn.q_at, dEdq) || dEdq.size() != m_atomcount)
        return;
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
