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
#include "xtb_coulomb.hpp"   // atom_hubbard_deriv_gfn1 (GFN1 resident loop)

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include "STO_CGTO.hpp"
#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"  // Molecule::GetFragments for the fragment SCF guess (Sep 2026)
#include "src/core/citation_registry.h"
#include "src/core/config_manager.h"
#include "src/core/charge_extrapolation.h"
#include "src/core/intra_parallel_context.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/energy_calculators/dispersion/d4_ncoord.h"  // D4 CN gradient (GFN2)
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"  // D3 CN gradient (May 2026)
#include "src/core/energy_calculators/ff_methods/alpb_solvation.h"  // implicit solvation (June 2026)
#include "diis_accelerator.h"
#include "broyden_mixer.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <future>
#include <stdexcept>
#include <thread>
#include <vector>

namespace curcuma::xtb {

namespace {
/// Transpose a square n x n matrix in place, in 64 x 64 tiles over up to 8 host threads
/// (Claude Generated, Sep 29, 2026). Pure element swaps, so the result is exact. Used to
/// turn a column-major device download into curcuma's row-major Matrix without a temporary.
void transposeSquareInPlace(double* a, int n)
{
    constexpr int T = 64;
    const int nb = (n + T - 1) / T;
    const unsigned nt = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
    auto work = [&](unsigned t) {
        for (int bi = static_cast<int>(t); bi < nb; bi += static_cast<int>(nt)) {
            const int i0 = bi * T, i1 = std::min(n, i0 + T);
            for (int bj = bi; bj < nb; ++bj) {
                const int j0 = bj * T, j1 = std::min(n, j0 + T);
                for (int i = i0; i < i1; ++i)
                    for (int j = (bi == bj ? i + 1 : j0); j < j1; ++j)
                        std::swap(a[static_cast<size_t>(i) * n + j], a[static_cast<size_t>(j) * n + i]);
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(work, t);
    for (auto& th : pool) th.join();
}
} // namespace

/* ------------------------------------------------------------------------- *
 *  Lifecycle
 * ------------------------------------------------------------------------- */
XTB::XTB(MethodType method)
    : m_method(method)
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
    // Already running under molecule-level parallelism → stay within the batch budget
    // (1 = serial for CPU batches, no N^2; cores/GPU-slots for GPU batch workers).
    const int cap = curcuma::intraParallelSuppressed() ? curcuma::intraThreadBudget() : m_intra_threads;
    if (cap <= 1) return 1;
    // Size guard: keep at least kMinWorkPerThread items per thread; tiny systems
    // are dominated by dispatch overhead (tuned in the benchmark step).
    const int by_work = work_units / kMinWorkPerThread;
    const int t = std::min(std::min(m_intra_threads, cap), std::max(1, by_work));
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
    // Reset hard-error state for the new molecule (Claude Generated): the B0
    // d-shell guard below sets it, so a stale flag must not leak from a previous
    // (rejected) system into a later valid one.
    m_has_error = false;
    m_error_message.clear();
    m_cn_cache_valid = false;  // X-I3: invalidate CN cache for the new system

    if (m_atomcount <= 0) {
        CurcumaLogger::error("XTB::InitialiseMolecule: no atoms set");
        return false;
    }

    // ── Citation registry (Claude Generated, June 2026) ──
    // Register the method and all sub-component citations, matching xtb output.
    const std::string meth = (m_method == MethodType::GFN1) ? "gfn1" : "gfn2";
    CitationRegistry::cite("xtb");
    CitationRegistry::cite(meth, "xtb");
    if (m_method == MethodType::GFN1)
        CitationRegistry::cite("d3", meth);        // D3(BJ) dispersion (GFN1 only)
    else
        CitationRegistry::cite("d4", meth);        // D4 dispersion (GFN2 only)
    if (m_solvent != "none" && !m_solvent.empty()) {
        const int model = (m_solvent_model > 0) ? m_solvent_model : 3;
        if (model == 3) CitationRegistry::cite("alpb", meth);
        else if (model == 2) { CitationRegistry::cite("alpb", meth); CitationRegistry::cite("still_gb", meth); }
    }
    CitationRegistry::cite("diis", meth);         // SCF convergence (DIIS/Broyden)

    buildBasis();

    // X-I1 (Claude Generated): d-shell basis functions are now supported on the CPU
    // path (cartesian->spherical dtrafo overlap/H0/multipole + gradients, validated to
    // <=1e-8 Eh vs tblite on H2S/PH3/HCl/SiH4). The former B0 loud-error guard was
    // removed here once d was implemented. Main-group d (S/P/Cl/Si/...) is validated;
    // transition metals are enabled but NOT validated (no reference set) -- see
    // docs/SQM_DSHELL_WP.md. GPU device d kernels are pending; d systems fall back to
    // the CPU integral/SCF/gradient path (m_has_dshell gate in buildBasis).

    buildH0Data();
    buildReferenceOccupations();
    // A fresh molecule invalidates any SCC-extrapolation history from a previous
    // system (different nsh/nat). Claude Generated.
    m_scf_history.clear();
    m_xlbomd_aux.clear();
    // A new molecule needs a new SCF preconditioner (-scf_precond frag). Claude Generated (Oct 2026).
    m_scf_precond_frag.reset();

    // ── Implicit solvation (Claude Generated, June 2026) ──
    // Build the self-consistent reaction-field model when a solvent is requested.
    // Created here (atoms known); update() refreshes the geometry-dependent state
    // each Calculation(). solvent_model: 3=ALPB, 2=GBSA (both wired); 1=CPCM is a
    // later WP. GBSA = ALPB with the shape term off (see ALPBSolvation::setUseAlpb).
    m_solvation.reset();
    m_E_solvation = 0.0;
    if (m_solvent != "none" && !m_solvent.empty()) {
        const int model = (m_solvent_model > 0) ? m_solvent_model : 3;  // default ALPB
        const std::string meth = (m_method == MethodType::GFN1) ? "gfn1" : "gfn2";
        if (model == 3 || model == 2) {
            const bool use_alpb = (model == 3);
            const char* tag = use_alpb ? "ALPB" : "GBSA";
            auto solv = std::make_unique<ALPBSolvation>();
            solv->setUseAlpb(use_alpb);
            if (solv->init(m_atoms, m_solvent, meth)) {
                m_solvation = std::move(solv);
            } else {
                CurcumaLogger::warn(std::string("Native solvation: no ") + tag
                                    + " parameters for solvent '" + m_solvent + "' ("
                                    + meth + "); running gas phase");
            }
        } else {
            CurcumaLogger::warn("Native solvation: only solvent_model alpb and gbsa "
                                "are implemented for native GFN (cpcm pending); "
                                "running gas phase");
        }
    }
    return true;
}

// ═══════════════════════════════════════════════════════════════════
// Implicit solvation hooks (Claude Generated, June 2026).
// addSolvationPotential adds the in-SCF Fock contribution v_at += dE_solv/dq;
// energySolvation returns the solvation free energy at the current charges.
// Both are no-ops when no solvent is active. See docs/SQM_SOLVATION_WP.md.
// ═══════════════════════════════════════════════════════════════════
void XTB::addSolvationPotential(Potential& pot) const
{
    if (m_solvation)
        m_solvation->addPotential(m_wfn.q_at, pot.v_at);
}

double XTB::energySolvation() const
{
    return m_solvation ? m_solvation->energy(m_wfn.q_at) : 0.0;
}

// Non-member helper for convergence check (must be defined before Calculation)
namespace {
    // dq is the max-abs residual of the SCC vector the mixer works on: shell charges
    // for GFN1, shell charges + atomic dipole/quadrupole moments for GFN2. tblite
    // converges on the same two conditions (xtb/singlepoint.f90:254-256:
    // `econverged = |sum(eelec) - elast| < econv`, `pconverged = mixer%get_error() < pconv`).
    bool checkConvergence_impl(double dq, double e_old, double e_new, double thresh)
    {
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
 *  Fragment initial guess (Claude Generated, Sep 2026; operator's proposal).
 *
 *  Large clusters of many small molecules (mixture2: 1000 water + 400 urea,
 *  6200 atoms) have a tiny GFN2 HOMO-LUMO gap (~0.1 eV per xtb 6.7.1), and the
 *  first Fock matrix built from EEQ or bare-H0 charges moves up to ~2.5 e
 *  between distant molecules; the SCF then runs away (measured: energy
 *  +2e5 Eh by iteration 3, on CPU and GPU, FP32 and FP64). Starting from the
 *  converged state of every isolated molecule puts iteration 0 next to the
 *  physical solution. Each fragment is converged in its actual geometry, so the
 *  GFN2 atomic dipoles/quadrupoles are already in the lab frame - no rotation
 *  or molecule-type matching needed. Fragments = Molecule::GetFragments(), the
 *  same partition large_system_mode=fragments uses.
 * ------------------------------------------------------------------------- */
bool XTB::seedFragmentGuess(Vector& q_sh_out, Matrix& dp_out, Matrix& qp_out)
{
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;
    if (nat <= 0 || nsh <= 0)
        return false;
    // A net charge cannot be assigned to fragments unambiguously - same limitation as
    // large_system_mode=fragments. The caller falls back to the EEQ guess.
    if (std::abs(static_cast<double>(m_charge)) > 1.0e-12)
        return false;

    const auto t0 = std::chrono::high_resolution_clock::now();

    Mol full;
    full.m_number_atoms = nat;
    full.m_charge = 0;
    full.m_spin = 0.0;
    full.m_atoms = m_atoms;
    full.m_geometry = m_geometry;
    const std::vector<std::vector<int>> fragments = Molecule(full).GetFragments();
    if (fragments.size() <= 1)
        return false;

    q_sh_out = Vector::Zero(nsh);
    dp_out = Matrix::Zero(3, nat);
    qp_out = Matrix::Zero(6, nat);

    // The fragment SCFs are an implementation detail of the guess: keep their output
    // (one SCF table per fragment) out of the log, restore the level afterwards.
    const int saved_verbosity = CurcumaLogger::get_verbosity();
    CurcumaLogger::set_verbosity(0);
    bool ok = true;
    int max_iter = 0, not_converged = 0;
    for (const auto& atoms : fragments) {
        Mol sub;
        sub.m_number_atoms = static_cast<int>(atoms.size());
        sub.m_charge = 0;
        sub.m_spin = 0.0;
        sub.m_geometry = Matrix(atoms.size(), 3);
        for (std::size_t k = 0; k < atoms.size(); ++k) {
            sub.m_atoms.push_back(m_atoms[atoms[k]]);
            sub.m_geometry.row(k) = m_geometry.row(atoms[k]);
        }
        XTB frag(m_method);
        frag.setScfGuess("eeq");   // never recurse into the fragment guess
        frag.setIntraThreads(1);
        if (!frag.InitialiseMolecule(sub)) { ok = false; break; }
        frag.Calculation(false);
        if (!frag.m_wfn.q_sh.allFinite()) { ok = false; break; }
        max_iter = std::max(max_iter, frag.m_scf_iterations);
        if (!frag.m_scf_converged) ++not_converged;
        for (std::size_t k = 0; k < atoms.size(); ++k) {
            const int a = atoms[k];
            const int n = m_basis.nsh_at[a];
            if (frag.m_basis.nsh_at[k] != n) { ok = false; break; }
            q_sh_out.segment(m_basis.ish_at[a], n) = frag.m_wfn.q_sh.segment(frag.m_basis.ish_at[k], n);
            if (m_method == MethodType::GFN2) {
                dp_out.col(a) = frag.m_wfn.dp_at.col(k);
                qp_out.col(a) = frag.m_wfn.qp_at.col(k);
            }
        }
        if (!ok) break;
    }
    CurcumaLogger::set_verbosity(saved_verbosity);
    if (!ok)
        return false;

    const double ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();
    if (saved_verbosity >= 1)
        CurcumaLogger::result(fmt::format(
            "SCF initial guess: {} fragments converged separately ({:.0f} ms, <= {} iterations{})",
            fragments.size(), ms, max_iter,
            not_converged ? fmt::format(", {} not converged", not_converged) : std::string()));
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
    // Shared with the GFN-FF EEQ PCG warm-start (src/core/charge_extrapolation.h).
    return curcuma::extrapolation::weights(mode, order, n_hist);
}

/* ------------------------------------------------------------------------- *
 *  Niklasson dissipative XL-BOMD coefficients (J. Chem. Phys. 130, 214109
 *  (2009), Table I). The integrator is
 *    P_{n+1} = 2 P_n - P_{n-1} + kappa (D[P_n] - P_n) + alpha Σ_{k=0}^{K} c_k P_{n-k},
 *  where D[P_n] is the SCF map applied to P_n. kappa is the spring constant,
 *  alpha scales the weak dissipation that drains accumulated numerical noise, and
 *  the c_k (which sum to zero, so the dissipation does not shift the fixpoint) are
 *  tabulated by the dissipation order K. Claude Generated (Phase 2).
 * ------------------------------------------------------------------------- */
bool XTB::xlbomdCoefficients(int K, double& kappa, double& alpha, std::vector<double>& c)
{
    switch (K) {
    case 3: kappa = 1.69; alpha = 0.150;  c = { -2,  3,   0,  -1 };                 return true;
    case 4: kappa = 1.75; alpha = 0.057;  c = { -3,  6,  -2,  -2,  1 };             return true;
    case 5: kappa = 1.82; alpha = 0.018;  c = { -6, 14,  -8,  -3,  4,  -1 };        return true;
    case 6: kappa = 1.84; alpha = 0.0055; c = { -14, 36, -27, -2, 12,  -6,  1 };    return true;
    case 7: kappa = 1.86; alpha = 0.0016; c = { -36, 99, -88, 11, 32, -25,  8, -1 };return true;
    default: return false;
    }
}

void XTB::setHardError(const std::string& msg)
{
    m_has_error = true;
    if (m_error_message.empty())
        m_error_message = msg;
    CurcumaLogger::error("XTB engine: " + msg);
}

double XTB::Calculation(bool gradient)
{
    // Reset the hard-error state for this evaluation (Claude Generated): genuine
    // solve breakdowns below set it so the wrapper refuses to return E=0 silently.
    m_has_error = false;
    m_error_message.clear();

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

    // 1. Coordination numbers + 2. self-energies. Claude Generated (Oct 2026): on the device-integral path the
    // device computes both itself (beginComputed) and the CN is downloaded from there, so the host O(nat^2) loop
    // is only run when the host builds the integrals (ensure_host_cn) or a host consumer needs it later.
    Vector cn;
    Vector se;
    bool cn_ready = false;
    auto ensure_host_cn = [&]() {
        if (!cn_ready) {
            cn = computeCoordinationNumbers();
            getSelfEnergies(cn, se);
            cn_ready = true;
        }
    };
    const bool cn_from_device_ok = m_gpu_scf && m_scf_mode == ScfMode::Broyden && !std::getenv("CURCUMA_CN_HOST");
    if (!cn_from_device_ok) ensure_host_cn();
    const auto t_cn = clock::now();

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
    // Claude Generated (Sep 29, 2026): device path only - split of the "overlap + H0" bucket into
    // the device integral build and the download of S/H0/L/gamma into the host matrices.
    auto t_dev_built = t_cn;
    bool gpu_computed = false;        // device integral build succeeded (reused below)
    bool integrals_from_device = false;
    // X-I1: device d kernels are backend-gated. CUDA supports them; ROCm/Vulkan
    // do not yet, so for those d systems use the validated CPU integral/SCF/gradient
    // path (atom-based D4/EEQ stay on device).
    const bool gpu_d_ok = (m_gpu_scf && m_gpu_scf->supportsDshell());
    if (m_gpu_scf && m_has_dshell && !gpu_d_ok) {
        static bool warned = false;
        if (!warned && verb >= 1) {
            warned = true;
            CurcumaLogger::warn("native xTB: d-shell system - integrals/SCF/gradient "
                                "run on CPU (device d kernels pending for this backend); "
                                "D4/EEQ stay on device");
        }
    }
    if (m_gpu_scf && m_scf_mode == ScfMode::Broyden && (!m_has_dshell || gpu_d_ok)) {
        GpuBasisFlat gbf;
        GpuH0Flat    gh0;
        exportGpuBasis(gbf, gh0);
        bool basis_ok = true;
        if (m_gpu_basis_dirty) {
            basis_ok = m_gpu_scf->beginBasis(gbf, gh0);
            if (basis_ok) m_gpu_basis_dirty = false;
        }
        gpu_computed = basis_ok && m_gpu_scf->beginComputed(gbf.xyz_bohr);
        t_dev_built = clock::now();
        if (gpu_computed) {
            const int nao = m_basis.nao, nsh = m_basis.nsh;
            // Claude Generated (Sep 29, 2026): download straight into the members. S and H0 come
            // column-major into the row-major m_S/m_H0 and are then transposed in place, which
            // gives exactly the element-wise copy the old temporaries made (m_S(i,j) = S(i,j))
            // without two 1.9 GB temporaries and four serial copies (polymer_2x: 2.4 s -> see
            // docs/GFN2_GPU_COST_PLAN.md stage 1a). Backends without the raw download keep the
            // old path.
            const auto td0 = clock::now();
            bool got = false;
            // Claude Generated (Oct 2026): lazy host copies (ensureHostIntegrals); CURCUMA_HOST_INTS_EAGER=1 restores
            // the immediate download. The host-driven loops, the host gradient and the responses ask for them.
            const bool defer_ints = m_scf_mode == ScfMode::Broyden && !std::getenv("CURCUMA_HOST_INTS_EAGER");
            if (defer_ints) {
                got = true;
                m_host_ints_deferred = true;
                m_S.resize(0, 0); m_H0.resize(0, 0); m_X.resize(0, 0); m_gamma.resize(0, 0);
            } else {
            m_S.resize(nao, nao);
            m_H0.resize(nao, nao);
            if (m_gpu_scf->downloadOverlapInto(m_S.data()) && m_gpu_scf->downloadH0Into(m_H0.data())) {
                transposeSquareInPlace(m_S.data(), nao);
                transposeSquareInPlace(m_H0.data(), nao);
                got = m_gpu_scf->downloadCholesky(m_X) && m_gpu_scf->downloadGamma(m_gamma);
            } else {
                Eigen::MatrixXd Scm(nao, nao), H0cm(nao, nao), Lcm(nao, nao), Gcm(nsh, nsh);
                if (m_gpu_scf->downloadOverlap(Scm) && m_gpu_scf->downloadH0(H0cm)
                    && m_gpu_scf->downloadCholesky(Lcm) && m_gpu_scf->downloadGamma(Gcm)) {
                    m_S     = Scm;   // col-major device → row-major Matrix
                    m_H0    = H0cm;
                    m_X     = Lcm;   // m_X is column-major Eigen::MatrixXd → direct
                    m_gamma = Gcm;
                    got = true;
                }
            }
            }
            if (got) {
                if (verb >= 3)
                    CurcumaLogger::info_fmt("    S/H0/L/gamma download into the host copies: {:.0f} ms",
                                            ms(td0, clock::now()));
                integrals_from_device = true;
                if (!cn_ready) {
                    Vector cn_dev(m_atomcount);
                    if (m_gpu_scf->downloadCn(cn_dev.data(), m_atomcount)) {
                        cn = cn_dev;
                        m_cn_cache = cn;            // the memo of computeCoordinationNumbers() for this geometry
                        m_cn_cache_valid = true;
                        cn_ready = true;
                    }
                }
                if (verb >= 2)
                    CurcumaLogger::info("SCF: integrals built on GPU device "
                                        "(S/H0/Lowdin/gamma downloaded; host build skipped)");
            }
        }
    }
    // B0 (Jul 2026): time the three host build steps separately. They used to share
    // one stamp (t_h0 == t_gamma), so "overlap + H0" silently also contained
    // buildOrthonormalizer() + buildGammaMatrix() and "Coulomb gamma" always read
    // 0.00 ms — which hid how much of the bucket the integral kernel actually owns.
    auto t_h0 = t_cn;
    auto t_ortho = t_cn;
    auto t_gamma = t_cn;
    if (!integrals_from_device) {
        // Full host integral build (CPU path, or device build unavailable).
        ensure_host_cn();
        getHamiltonianH0(se, S, H0);
        m_S  = S;
        m_H0 = H0;
        t_h0 = clock::now();
        // Orthonormalizer X = S^{-1/2}, built once here so every SCF iteration solves
        // the cheap standard eigenproblem instead of re-factorizing the constant S.
        buildOrthonormalizer();
        t_ortho = clock::now();
        buildGammaMatrix();   // Coulomb gamma (once per geometry)
        t_gamma = clock::now();
    } else {
        // Device path: S/H0/Lowdin/gamma were fused in the device build above.
        t_h0 = t_ortho = t_gamma = clock::now();
    }

    // 5. Multipole setup (GFN2 only) — fills m_dp_int, m_qp_int, interaction matrices.
    // Stage 3m (Claude Generated): on the non-resident device GFN2 path (Vulkan/ROCm),
    // build dp_int/qp_int ON the device and download them, skipping the host O(nao²)
    // integral loop. CUDA's resident multipole loop (supportsMultipole) builds AND
    // consumes them on-device, so it is excluded here and keeps the full host build.
    if (m_method == MethodType::GFN2) {
        bool mp_on_device = false;
        // The device built dp_int/qp_int during the integral build (computeIntegrals);
        // download them for the host GFN2 gradient (multipole Pulay term) and skip the
        // host O(nao²) integral loop. Independent of whether the SCF runs the resident
        // multipole loop (Stage 2b) — that is set up in the SCF block below. CUDA's
        // adapter does not override downloadMultipoleInts (its gradient is on-device), so
        // this is a no-op there and CUDA keeps its full host build (unchanged).
        // Claude Generated (Oct 2026, ROCm WP8b): a backend with deferHostCopies() (CUDA: false) skips the download
        // where the deferral below applies, so the host never holds the 9 nao^2 integrals unless a fallback needs them.
        const bool backend_defers_mp = integrals_from_device && m_gpu_scf && m_gpu_scf->deferHostCopies()
            && m_gpu_scf->supportsResidentLoop() && m_scf_mode == ScfMode::Broyden && m_d4_charge_source != "cpscf"
            && (9.0 * static_cast<double>(m_basis.nao) * m_basis.nao * sizeof(double) > 2.0e9
                || !std::getenv("CURCUMA_MP_INTS_EAGER"));
        if (integrals_from_device && m_gpu_scf && !backend_defers_mp
            && m_gpu_scf->downloadMultipoleInts(m_dp_int, m_qp_int)) {
            mp_on_device = true;
            if (verb >= 2)
                CurcumaLogger::info("SCF: GFN2 multipole integrals built on GPU device "
                                    "(dp_int/qp_int downloaded; host build skipped)");
        }
        // Claude Generated (Sep 2026): a large GFN2 system on the CUDA device-resident path
        // never reads the host dense multipole integrals (9 nao^2 doubles, 17 GB and a long
        // O(nao^2) host loop at nao = 15444). Defer them; the host fallbacks (host SCF loop,
        // host gradient, CPSCF D4 response) build them on demand. Small systems (< ~2 GB of
        // integrals) keep the eager build, so their results and debug dumps are unchanged.
        const double mp_bytes = 9.0 * static_cast<double>(m_basis.nao) * m_basis.nao * sizeof(double);
        const bool defer = !mp_on_device && integrals_from_device && m_gpu_scf
            && m_gpu_scf->supportsResidentLoop() && m_scf_mode == ScfMode::Broyden
            && m_d4_charge_source != "cpscf" && (mp_bytes > 2.0e9 || !std::getenv("CURCUMA_MP_INTS_EAGER"));
        // Claude Generated (Oct 2026): the threshold is gone - the resident loop and the device gradient build the
        // AO dipole/quadrupole integrals themselves, so the eager host O(nao^2) loop (712 ms on polymer, nao 3222)
        // was pure host work; the host fallbacks make them on demand (ensureHostMultipoleIntegrals).
        // CURCUMA_MP_INTS_EAGER=1 restores the eager build below 2 GB.
        // Claude Generated (Oct 2026): the nat x nat interaction matrices (host O(nat^2) build, 0.23 s at 2910
        // atoms, plus an 18 nat^2 upload) are built by the device on the CUDA resident path; the host copies are
        // made on demand by the host fallbacks (ensureHostMultipoleMatrices). CURCUMA_MP_AMAT_DEVICE=0 keeps
        // the host build.
        bool defer_amat = integrals_from_device && m_gpu_scf && m_gpu_scf->supportsResidentLoop()
            && m_gpu_scf->supportsStoredMultipoleOnDevice() && m_scf_mode == ScfMode::Broyden;
        if (const char* e = std::getenv("CURCUMA_MP_AMAT_DEVICE")) defer_amat = defer_amat && (e[0] != '0');
        setupMultipole(mp_on_device || defer, defer_amat);
        m_mp_ints_deferred = defer;
        if (defer && verb >= 2)
            CurcumaLogger::info(fmt::format("SCF: host GFN2 multipole integrals deferred ({:.1f} GB; "
                                            "the device-resident loop builds its own)", mp_bytes / 1.0e9));
    }

    // Implicit solvation: refresh the geometry-dependent state (Born radii, SASA,
    // Born matrix, CM5) once per geometry, before the SCF loop. m_geometry is in
    // Angstrom; the model works in Bohr (Claude Generated, June 2026).
    if (m_solvation)
        m_solvation->update(m_atoms, m_geometry * AA_TO_AU);
    const auto t_setup = clock::now();
    // Claude Generated (Oct 2026): host timing of the stretch between the integral setup and the first SCF iteration
    // (CURCUMA_GPU_PROFILE / verbosity 3), to see what is still done on the host before the device loop starts.
    std::vector<std::pair<std::string, double>> pre_t;
    auto pre_tp = clock::now();
    auto pre_mark = [&](const char* what) {
        const auto now = clock::now();
        pre_t.emplace_back(what, ms(pre_tp, now));
        pre_tp = now;
    };

    if (verb >= 3) {
        CurcumaLogger::info("Setup timing:");
        CurcumaLogger::info_fmt("  coordination numbers : {:8.2f} ms", ms(t0, t_cn));
        CurcumaLogger::info_fmt("  overlap + H0         : {:8.2f} ms", ms(t_cn, t_h0));
        if (integrals_from_device) {
            CurcumaLogger::info_fmt("    device build (S/H0/L/gamma) : {:8.2f} ms", ms(t_cn, t_dev_built));
            CurcumaLogger::info_fmt("    host download + copy        : {:8.2f} ms", ms(t_dev_built, t_h0));
        }
        CurcumaLogger::info_fmt("  orthonormalizer      : {:8.2f} ms", ms(t_h0, t_ortho));
        CurcumaLogger::info_fmt("  Coulomb gamma matrix : {:8.2f} ms", ms(t_ortho, t_gamma));
        if (m_method == MethodType::GFN2) {
            CurcumaLogger::info_fmt("  multipole setup      : {:8.2f} ms", ms(t_gamma, t_setup));
            // B0: attribute the multipole bucket to its five sub-phases.
            CurcumaLogger::info_fmt("    AO dp/qp integrals : {:8.2f} ms", m_mp_t_ao_ints);
            CurcumaLogger::info_fmt("    d-shell block      : {:8.2f} ms", m_mp_t_d_block);
            CurcumaLogger::info_fmt("    origin shift       : {:8.2f} ms", m_mp_t_shift);
            CurcumaLogger::info_fmt("    CN + mrad          : {:8.2f} ms", m_mp_t_cn_mrad);
            CurcumaLogger::info_fmt("    interaction mats   : {:8.2f} ms", m_mp_t_amat);
        }
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

    const int nsh = m_basis.nsh;
    const int packed_len = (m_method == MethodType::GFN2) ? (nsh + 9 * m_atomcount) : nsh;

    // XL-BOMD (extended-Lagrangian Born-Oppenheimer MD) — opt-in, EXPERIMENTAL.
    // Each Calculation() is treated as one MD/geometry step. When a time-reversibly
    // propagated auxiliary SCC trajectory (m_xlbomd_aux, length K+1) is available, this
    // step SEEDS the SCF from the auxiliary front() and runs the normal SCF to
    // convergence (the auxiliary is a near-perfect guess, so this is few iterations);
    // the auxiliary is then advanced post-loop by the time-reversible Niklasson
    // dissipative integrator from the CONVERGED density D[P_n]. (A naive "few bare maps,
    // no convergence" corrector is NOT contractive for tight-binding SCF and diverges,
    // so the corrector converges; XL-BOMD's value here is the time-reversible auxiliary
    // for low MD energy drift, plus the guess-driven iteration savings.) The first K+1
    // steps bootstrap by running the normal SCF and seeding the trajectory with
    // converged densities. Claude Generated (Phase 2).
    const bool xlbomd_active = (m_scf_extrapolation != "none" && m_scf_extrap_apply == "xlbomd");
    const int  xlbomd_K = std::min(std::max(m_scf_extrap_order, 3), 7); // dissipation order (table 3..7)
    bool   xlbomd_corrector_step = false;
    Vector xlbomd_Pn;
    if (xlbomd_active && static_cast<int>(m_xlbomd_aux.size()) >= xlbomd_K + 1
        && m_xlbomd_aux.front().size() == packed_len) {
        xlbomd_corrector_step = true;
        xlbomd_Pn = m_xlbomd_aux.front();
    }
    if (xlbomd_active && !m_xlbomd_warned) {
        m_xlbomd_warned = true;
        if (verb >= 1)
            CurcumaLogger::warn("scf_extrapolation_apply=xlbomd is EXPERIMENTAL "
                                "(extended-Lagrangian; MD energy-drift unvalidated). "
                                "Use with scf_extrapolation=aspc and a fixed MD timestep.");
    }

    const int    max_iter = m_scf_max_iter;
    const double thresh   = m_scf_threshold;

    // Initial guess. Warm-start: reuse converged shell charges from the previous
    // geometry (saved by UpdateMolecule). Otherwise: "h0" (zero) or "eeq".
    // Claude Generated warm-start path.
    Vector q_sh_old = Vector::Zero(nsh);
    bool guess_set = false;

    // XL-BOMD corrector step: seed the SCF from the time-reversibly propagated
    // auxiliary density, then converge normally below. Takes precedence over the
    // guess-mode extrapolation. Claude Generated (Phase 2).
    if (xlbomd_corrector_step) {
        unpackSccState(xlbomd_Pn);
        q_sh_old  = m_wfn.q_sh;
        guess_set = true;
        CurcumaLogger::citation("xlbomd");
        if (verb >= scf_min)
            CurcumaLogger::result(fmt::format(
                "SCF guess: XL-BOMD auxiliary density (K={})", xlbomd_K));
    }

    // (0) Multi-step SCC extrapolation (opt-in, guess mode). Predict the new-step
    //     SCC vector from the history of converged steps so the SCF starts closer
    //     to the new fixpoint than the previous step alone (fewer iterations, esp.
    //     in MD). Only the initial guess changes; the SCF below still converges to
    //     the same fixpoint within scf_threshold. A poor prediction or too-short
    //     history simply falls through to the 1-step warm-start. Also runs during the
    //     XL-BOMD bootstrap steps (guess_set still false there). Claude Generated.
    if (!guess_set && m_warmstart && m_scf_extrapolation != "none"
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

    if (m_force_h0_guess) {
        // Runaway retry: start from the bare H0, i.e. do NOT take the extrapolated or
        // warm-started charges — those ARE the runaway solution we are trying to escape.
        // Claude Generated (Sep 2026).
        guess_set = false;
        q_sh_old.setZero(nsh);
        m_wfn.q_sh.setZero(nsh);
        m_wfn.q_at.setZero(m_atomcount);
        if (m_method == MethodType::GFN2) {
            m_wfn.dp_at.setZero(3, m_atomcount);
            m_wfn.qp_at.setZero(6, m_atomcount);
        }
    } else if (guess_set) {
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
    } else if ((m_scf_guess == "eeq" || m_scf_guess == "fragments") && !m_force_h0_guess) {
        bool seeded = false;
        if (m_scf_guess == "fragments") {
            Vector q_sh_guess;
            Matrix dp_guess, qp_guess;
            if (seedFragmentGuess(q_sh_guess, dp_guess, qp_guess)) {
                q_sh_old   = q_sh_guess;
                m_wfn.q_sh = q_sh_guess;
                m_wfn.q_at.setZero(m_atomcount);
                for (int s = 0; s < nsh; ++s)
                    m_wfn.q_at(m_basis.sh2at[s]) += q_sh_guess(s);
                if (m_method == MethodType::GFN2) {
                    m_wfn.dp_at = dp_guess;
                    m_wfn.qp_at = qp_guess;
                }
                seeded = true;
            } else if (verb >= scf_min) {
                CurcumaLogger::warn("Fragment initial guess unavailable (one fragment, charged "
                                    "system or a failed fragment SCF); using the EEQ guess");
            }
        }
        Vector q_sh_guess;
        if (seeded) {
            // fragment guess set above
        } else if (seedEEQGuess(q_sh_guess)) {
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
    // Claude Generated (Sep 2026): FP32 stagnation guard. The mixed-precision phase reverts to
    // FP64 once max|dq| < scf_fp32_threshold, but FP32 eigenvectors carry ~1e-7 relative noise,
    // which shows up as a floor of ~1e-5 in dq for a large system. With the GPU defaults
    // (scf_fp32_threshold = scf_threshold = 1e-5) the SCF then hovers AT that floor and only
    // converges when an iteration happens to dip below it - observed on 2x H200, polymer_2x:
    // 15 vs 21 iterations for runs that differ only by rounding. So: if dq stops improving while
    // still in FP32, switch to FP64 for the rest, where the iteration is exact and convergence
    // is monotone again. Costs nothing when the FP32 phase is converging normally.
    double dq_best_fp32 = 1.0e30;
    int    fp32_stall = 0;
    bool   fp32_exhausted = false;
    // Claude Generated (Oct 2026): CURCUMA_SCF_JAC=<dir> dumps the numerical Jacobian of the SCC
    // fixed-point map after convergence (writeScfJacobianDump). m_scf_jac_capture (set only on the
    // fragment objects of -scf_precond frag, xtb_precond.cpp) computes the same Jacobian into
    // m_scf_jac_result instead; such objects ignore CURCUMA_SCF_JAC. In both modes the whole SCF
    // runs in FP64, so x* is an FP64 fixed point. Neither set: no effect on any numerical path.
    const char* scf_jac_dir = m_scf_jac_capture ? nullptr : std::getenv("CURCUMA_SCF_JAC");
    const bool  scf_jac_dump = (scf_jac_dir != nullptr && scf_jac_dir[0] != '\0');
    const bool  scf_jac_active = scf_jac_dump || m_scf_jac_capture;
    Vector scf_jac_x_star, scf_jac_g_star;   // converged input x* and its output G(x*) (CPU loop)
    double scf_jac_e_scc = 0.0;
    auto fp32_wanted = [&](double dq_last) {
        return !scf_jac_active && m_scf_mixed_precision && !fp32_exhausted
               && (dq_last > m_scf_fp32_threshold);
    };
    // Claude Generated (Sep 2026, measured on H200): the FP32 phase can converge to a fixed point
    // that is not the FP64 one. polymer_2x, nao 15444: FP32 reported max|dq| = 7.1e-6 at an energy
    // 1.1 kcal/mol off, and the next FP64 iteration - taken only because convergence is never
    // accepted on an FP32 step - showed the true residual to be 9.4e-3. Left alone the SCF then
    // bounces between the two precisions. So: if an FP64 iteration finds a residual far above what
    // FP32 last claimed, FP32 has been lying and the rest of the SCF runs in FP64.
    double dq_fp32_last = -1.0;
    auto note_fp64_reality_check = [&](double dq_now) {
        if (m_eig_fp32 || fp32_exhausted || dq_fp32_last < 0.0) return;
        if (m_fp32_false_fixpoint <= 0.0) return;   // check disabled
        if (dq_now > m_fp32_false_fixpoint * std::max(dq_fp32_last, 1.0e-12)) {
            fp32_exhausted = true;
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::warn_fmt("SCF: the FP32 phase converged to a false fixed point "
                                        "(FP32 max|dq| {:.2e}, FP64 says {:.2e}); continuing in FP64",
                                        dq_fp32_last, dq_now);
        }
    };
    // Call after every iteration that ran in FP32, with that iteration's max|dq|.
    auto note_fp32_progress = [&](double dq_now) {
        if (!m_eig_fp32) { note_fp64_reality_check(dq_now); return; }
        dq_fp32_last = dq_now;
        if (dq_now < 0.7 * dq_best_fp32) {     // still making real progress
            dq_best_fp32 = dq_now;
            fp32_stall = 0;
            return;
        }
        if (m_fp32_stall_patience > 0 && ++fp32_stall >= m_fp32_stall_patience) {
            fp32_exhausted = true;
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info_fmt("SCF: FP32 phase stalled at max|dq| = {:.2e} "
                                        "(>= scf_fp32_threshold {:.1e}); continuing in FP64",
                                        dq_now, m_scf_fp32_threshold);
        }
    };

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

    pre_mark("initial guess (EEQ / warm start / fragments)");
    // Intra-iteration timing buckets (printed at verbosity >= 3). A handful of
    // steady_clock reads per phase — negligible vs the ms-scale work — so they
    // run unconditionally, matching the existing per-iter timer. Claude Generated.
    double acc_pot = 0.0, acc_fock = 0.0, acc_solve = 0.0, acc_mull = 0.0, acc_energy = 0.0;
    double acc_e_coul = 0.0, acc_e_third = 0.0, acc_e_mp = 0.0, acc_e_band = 0.0;
    double acc_disp = 0.0;   // D4 in-SCF potential subset of acc_pot (GFN2), verbosity 3
    m_t_xfx = m_t_diag = m_t_back = m_t_dens = m_t_xfx_copy = 0.0;
    m_eig_calls_native = m_eig_calls_fp32 = m_eig_calls_lapack = 0;

    // Intra-molecule thread count for the per-iteration eigensolve. The eigensolve
    // (dsygst/dsyevd/dtrsm) is the one region handed to MKL rather than the
    // CxxThreadPool; effectiveIntraThreads() gates it (serial under molecule-level
    // parallelism or for small bases). nao is constant over the SCF. Claude Generated.
    //
    // `-threads N` decides (operator, Sep 2026): the eigensolve gets the same intra-molecule
    // budget as everything else. It used to be capped at 8 because the reduction ran in FP32
    // and dominated (see docs/SQM_PERFORMANCE.md); with the FP64 reduction the optimum moved
    // and is machine-dependent - polymer/nao 3222 on this 36-core box: 8 threads 27.1 s,
    // 16 threads 23.8 s, 24 threads 26.5 s, 36 threads 28.6 s. The D&C eigensolve is
    // memory-bandwidth-bound, so more threads than memory channels can still lose; that is now
    // a `-threads` choice. `-eigensolver_max_threads N` caps the eigensolve independently of
    // -threads (CURCUMA_EIG_MAX_THREADS is the older spelling and still wins over the flag).
    int eig_cap = m_eig_max_threads;   // 0 = no cap, follow -threads (-eigensolver_max_threads)
    if (const char* env = std::getenv("CURCUMA_EIG_MAX_THREADS")) {
        const int v = std::atoi(env);
        if (v > 0) eig_cap = v;        // the environment variable still wins
    }
    const int eig_intra = effectiveIntraThreads(m_basis.nao);
    const int eig_threads = eig_cap > 0 ? std::min(eig_intra, eig_cap) : eig_intra;
    if (verb >= 3 && eig_threads > 1)
        CurcumaLogger::info_fmt("Eigensolve BLAS/LAPACK threads: {}{}", eig_threads,
                                eig_cap > 0 ? " (eigensolver_max_threads cap)" : " (from -threads)");

    pre_mark("thread budgets, eigensolver settings");
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
    if (m_gpu_scf && mode == ScfMode::Broyden && (!m_has_dshell || m_gpu_scf->supportsDshell())
        && ((m_X.rows() == m_basis.nao && m_X.cols() == m_basis.nao) || m_host_ints_deferred)) {
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
        // (Oct 2026: already in m_gamma when the integrals came from the device - the setup downloaded it; the second
        // nsh^2 download cost 55 ms (GFN2) / 115 ms (GFN1) on peo201-water500 for the same numbers.)
        if (use_gpu_resident && computed && !integrals_from_device) {
            Eigen::MatrixXd dev_gamma;
            if (m_gpu_scf->downloadGamma(dev_gamma)
                && dev_gamma.rows() == m_gamma.rows()
                && dev_gamma.cols() == m_gamma.cols())
                m_gamma = dev_gamma;
        }
    }

    pre_mark("device-resident begin (H0/S/L, gamma download, multipole begin)");
    // Stage 5 (Part B3/B4, Claude Generated): full device GFN2 potential build.
    // When the backend supports it, the WHOLE per-iteration potential (γ·q_sh +
    // shell third-order + multipole v_dp/v_qp/v_at scalar shift + in-SCF D4) is
    // built on the device, so the SCF loop uploads only the mixed q_sh/dp_at/qp_at
    // (+ the host-built D4 reference weights) instead of v_ao. The geometry-fixed
    // reference data (D4 cache + multipole interaction matrices + third-order
    // hardness) is uploaded once here. Falls back to the host potential build on
    // any failure. GFN2 only (the D4/multipole coupling).
    //
    // Solvation on the device-potential path (WP4b, Claude Generated June 2026): the
    // device build adds the in-SCF reaction field v_at += B·q_at when the Born matrix
    // is uploaded (beginSolvation, below). This needs q_solute = q_at, which holds for
    // GFN2 (Mulliken) — deviceBornMatrix() returns the matrix only then, and nullptr
    // for the GFN1 CM5 path. If solvation is active but the device cannot take it
    // (no born matrix / backend lacks supportsDeviceSolvation), block the device path
    // so GFN2+solvent falls back to the host-driven device loop (WP4a), which folds
    // addSolvationPotential into the uploaded v_ao — still correct, still GPU eigensolve.
    const double* solv_born = m_solvation ? m_solvation->deviceBornMatrix() : nullptr;
    const bool solv_blocks_device =
        m_solvation && !(solv_born && m_gpu_scf && m_gpu_scf->supportsDeviceSolvation());
    bool use_device_potential = false;
    // Claude Generated (Sep 2026): whether the 18 nat^2 multipole interaction matrices sit on
    // the device (stored path) and how much that is - used to explain a device-memory failure
    // of the resident loop, which is where the shortage actually surfaces.
    bool   mp_stored_on_device = false;
    double mp_stored_bytes     = 0.0;
    if (use_gpu_resident && gpu_multipole && m_method == MethodType::GFN2 && m_mp_initialized
        && !solv_blocks_device
        && m_gpu_scf->supportsDeviceDispersion() && m_gpu_scf->supportsDevicePotential()) {
        // Prime the D4 reference set + device beginDispersion (the m_d4_prepared
        // first-call path); the returned potential is discarded here.
        Vector warm_dedq;
        computeD4PotentialDedq(m_wfn.q_at, warm_dedq);
        // Flatten + upload the multipole interaction matrices (column-major nat²
        // blocks), the on-site XC kernels, and the per-shell third-order hardness
        // Γ_s (= ½·thirdOrderKernelDiag at q_sh=1, GFN2 shell-resolved).
        const size_t nn = static_cast<size_t>(nat) * static_cast<size_t>(nat);
        std::vector<double> dk(nat), qk(nat), g3(nsh);
        for (int i = 0; i < nat; ++i) { dk[i] = m_mp_dkernel[i]; qk[i] = m_mp_qkernel[i]; }
        const Vector g3diag = thirdOrderKernelDiag(Vector::Ones(nsh), m_wfn.q_at);  // = 2·Γ_s
        for (int s = 0; s < nsh; ++s) g3[s] = 0.5 * g3diag(s);
        // Claude Generated (Sep 2026): above ~1 GB of interaction matrices (nat > ~2700) the
        // device rebuilds the matrix elements per iteration instead of storing 18 nat^2 doubles
        // (7.7 GB at 7320 atoms, plus the same again as host upload copies).
        // `-gpu_multipole_otf on|off` forces the choice (CURCUMA_GPU_MP_OTF=1/0 still wins).
        bool otf = 18.0 * static_cast<double>(nn) * sizeof(double) > 1.0e9;
        mp_stored_bytes = 18.0 * static_cast<double>(nn) * sizeof(double);
        bool otf_forced = true;   // Claude Generated (Oct 2026, ROCm WP11)
        if (m_gpu_mp_otf == "on" || m_gpu_mp_otf == "true")  otf = true;
        else if (m_gpu_mp_otf == "off" || m_gpu_mp_otf == "false") otf = false;
        else otf_forced = false;
        if (const char* e = std::getenv("CURCUMA_GPU_MP_OTF")) { otf = (e[0] == '1'); otf_forced = true; }
        // ROCm WP11: the automatic choice keeps the stored matrices while they fit the backend's own limit.
        if (otf && !otf_forced && m_gpu_scf->preferStoredMultipole(nat)) otf = false;
        if (otf && m_gpu_scf->supportsOnTheFlyMultipole()) {
            std::vector<double> xyzb(3 * static_cast<size_t>(nat));
            for (int i = 0; i < nat; ++i)
                for (int k = 0; k < 3; ++k) xyzb[3 * i + k] = m_geometry(i, k) * AA_TO_AU;
            use_device_potential = m_gpu_scf->beginPotentialOnTheFly(
                nat, nsh, xyzb.data(), m_mp_mrad.data(), gfn2_params::mp_dmp3, gfn2_params::mp_dmp5,
                dk.data(), qk.data(), g3.data());
        } else {
            // Device build of the stored matrices first (no host matrices, no 18 nat^2 upload); the host path
            // below is the fallback and needs the host copies.
            if (m_gpu_scf->supportsStoredMultipoleOnDevice()) {
                std::vector<double> xyzb(3 * static_cast<size_t>(nat));
                for (int i = 0; i < nat; ++i)
                    for (int k = 0; k < 3; ++k) xyzb[3 * i + k] = m_geometry(i, k) * AA_TO_AU;
                use_device_potential = m_gpu_scf->beginPotentialStored(
                    nat, nsh, xyzb.data(), m_mp_mrad.data(), gfn2_params::mp_dmp3, gfn2_params::mp_dmp5,
                    dk.data(), qk.data(), g3.data());
                if (use_device_potential) {
                    mp_stored_on_device = true;
                    if (verb >= 2) CurcumaLogger::info("SCF: GFN2 multipole interaction matrices built on the device");
                }
            }
            if (!use_device_potential) {
            ensureHostMultipoleMatrices();
            std::vector<double> sd(3 * nn), dd(9 * nn), sq(6 * nn);
            for (int k = 0; k < 3; ++k)
                std::memcpy(sd.data() + static_cast<size_t>(k) * nn, m_mp_amat_sd[k].data(), nn * sizeof(double));
            for (int a = 0; a < 3; ++a)
                for (int bb = 0; bb < 3; ++bb)
                    std::memcpy(dd.data() + (static_cast<size_t>(a) * 3 + bb) * nn,
                                m_mp_amat_dd[a][bb].data(), nn * sizeof(double));
            for (int k = 0; k < 6; ++k)
                std::memcpy(sq.data() + static_cast<size_t>(k) * nn, m_mp_amat_sq[k].data(), nn * sizeof(double));
            mp_stored_on_device = true;
            use_device_potential = m_gpu_scf->beginPotential(nat, nsh, sd.data(), dd.data(),
                                                             sq.data(), dk.data(), qk.data(), g3.data());
            // Claude Generated (Sep 2026): the stored path needs 18 nat^2 doubles on the device
            // (7.7 GB at 7320 atoms) and the same again here on the host. If that does not fit,
            // say so and rebuild the matrices per iteration instead of walking into an SCF that
            // cannot run - `-gpu_multipole_otf off` on polymer_2x used to fail at iteration 0
            // with no indication of the cause.
            if (!use_device_potential && m_gpu_scf->supportsOnTheFlyMultipole()) {
                const std::string why = m_gpu_scf->lastError();
                CurcumaLogger::warn("SCF: storing the GFN2 multipole matrices failed"
                                    + (why.empty() ? std::string() : " (" + why + ")")
                                    + "; rebuilding them per iteration instead");
                std::vector<double> xyzb(3 * static_cast<size_t>(nat));
                for (int i = 0; i < nat; ++i)
                    for (int k = 0; k < 3; ++k) xyzb[3 * i + k] = m_geometry(i, k) * AA_TO_AU;
                use_device_potential = m_gpu_scf->beginPotentialOnTheFly(
                    nat, nsh, xyzb.data(), m_mp_mrad.data(), gfn2_params::mp_dmp3,
                    gfn2_params::mp_dmp5, dk.data(), qk.data(), g3.data());
                mp_stored_on_device = !use_device_potential;
            }
            }   // !use_device_potential after the device build
        }
        // WP4b: upload the Born matrix so the device build adds v_at += B·q_at in-SCF.
        // On failure, drop to the host-driven loop (which still applies solvation).
        if (use_device_potential && solv_born) {
            if (!m_gpu_scf->beginSolvation(nat, solv_born)) {
                use_device_potential = false;
                if (verb >= 2)
                    CurcumaLogger::warn("SCF: device solvation upload failed; host potential build");
            }
        }
        if (verb >= 2)
            CurcumaLogger::info(use_device_potential
                ? "SCF: GFN2 potential (Coulomb+third-order+multipole+D4"
                  + std::string(solv_born ? "+solvation" : "") + ") built on the device"
                : "SCF: device potential build unavailable; host potential build");
    }

    // Claude Generated (Oct 2026), ROCm -xtb.rocm_device_multipole (opt-in): in the host-driven
    // resident GFN2 loop (no device potential build) only the multipole potential and energy
    // move to the device. The host interaction matrices are uploaded once per setupMultipole()
    // (m_mp_version); per iteration q_at/dp_at/qp_at go up and v_at/v_dp/v_qp come back, while
    // v_dp/v_qp also stay on the device for the Fock build. Any failure switches this object to
    // the host path for good (m_mp_device_failed); the backend prints one notice.
    bool use_dev_mp = false;
    if (use_gpu_resident && gpu_multipole && !use_device_potential && m_mp_initialized
        && !m_mp_device_failed && m_gpu_scf->supportsDeviceMultipolePotential()) {
        use_dev_mp = m_gpu_scf->beginMultipolePotential(m_mp_amat_sd, m_mp_amat_dd, m_mp_amat_sq,
                                                        m_mp_dkernel, m_mp_qkernel, m_mp_version);
        if (!use_dev_mp) {
            m_mp_device_failed = true;
            m_gpu_scf->noteDeviceMultipoleFallback("upload of the interaction matrices");
        }
        if (verb >= 2)
            CurcumaLogger::info(use_dev_mp
                ? "SCF: GFN2 multipole potential and energy on the device (host interaction matrices)"
                : "SCF: device multipole potential unavailable; host multipole potential");
    }
    // Device step failed: host path from here on (the matrices are still on the host).
    auto dev_mp_failed = [&](const char* step) {
        use_dev_mp = false;
        m_mp_device_failed = true;
        m_gpu_scf->noteDeviceMultipoleFallback(step);
    };

    pre_mark("device potential begin (D4 warm-up, multipole matrices, third order)");
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
        // Oct 16, 2026 (Known Issue #39): the loop's default occupation is the closed-shell filling from
        // the total electron count; an open shell needs the two spin fillings (n_alpha, n_beta) on the
        // device. A backend without them keeps the host-driven loop, whose occupations
        // (occupationsFromEps) are open-shell aware.
        if (use_resident_loop && openShell())
            use_resident_loop = m_gpu_scf->setResidentOpenShell(m_nalpha, m_nbeta);
        if (verb >= 2)
            CurcumaLogger::info(use_resident_loop
                ? "SCF: fully device-resident loop (host polls only dq + energy scalars)"
                : "SCF: resident-loop setup unavailable; host-driven device-potential loop");
    }

    pre_mark("resident loop begin (D4 weights tables, Broyden)");
    // Claude Generated (Oct 2026): GFN1 on the fully device-resident loop (shell charges only; v_sh = gamma q_sh,
    // v_at = Gamma_A q_A^2). Needs the device-built integrals (resident gamma, L). Solvation keeps the host-driven
    // GFN1 device loop unless the backend takes the CM5 reaction field (supportsDeviceSolvationGfn1, ROCm WP6:
    // v_at += B q_at + B cm5 on the device). `-gpu_gfn1_resident false` / CURCUMA_GFN1_RESIDENT=0 switch it off.
    bool gfn1_resident = false;
    const double* gfn1_born = nullptr;   // Claude Generated (Oct 2026, ROCm WP6)
    Vector gfn1_solv_shift;
    if (m_solvation && m_method == MethodType::GFN1 && m_gpu_scf && m_gpu_scf->supportsDeviceSolvationGfn1()) {
        if (const auto* alpb = dynamic_cast<const ALPBSolvation*>(m_solvation.get())) {
            if (alpb->bornMatrix().rows() == nat && alpb->bornMatrix().cols() == nat) {
                gfn1_born = alpb->bornMatrix().data();
                gfn1_solv_shift = Vector::Zero(nat);
                m_solvation->addPotential(Vector::Zero(nat), gfn1_solv_shift);   // B (0 + cm5) = B cm5
            }
        }
    }
    if (!use_resident_loop && m_method == MethodType::GFN1 && use_gpu_resident && gpu_computed
        && mode == ScfMode::Broyden && (!m_solvation || gfn1_born) && m_gpu_scf->supportsResidentLoopGfn1()) {
        bool want = m_gpu_gfn1_resident;
        if (const char* e = std::getenv("CURCUMA_GFN1_RESIDENT")) want = (e[0] != '0');
        if (want) {
            const int nocc_pairs = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
            const int b_maxhist = (m_diis_subspace > 2) ? 20 : m_diis_subspace;
            std::vector<double> g3at(nat);
            for (int i = 0; i < nat; ++i) g3at[i] = coulomb::atom_hubbard_deriv_gfn1(m_basis.z[i]);
            use_resident_loop = m_gpu_scf->beginResidentLoopGfn1(m_wfn.q_sh, m_wfn.q_at, m_wfn.n0_sh, m_wfn.n0_at,
                                                                 g3at, m_electronic_temp, m_wfn.nocc, nocc_pairs,
                                                                 damp, b_maxhist, 0.01);
            if (use_resident_loop && openShell())
                use_resident_loop = m_gpu_scf->setResidentOpenShell(m_nalpha, m_nbeta);
            if (use_resident_loop && gfn1_born) {
                use_resident_loop = m_gpu_scf->beginSolvationGfn1(nat, gfn1_born, gfn1_solv_shift.data());
                if (verb >= 2)
                    CurcumaLogger::info(use_resident_loop
                        ? "SCF: GFN1 solvation reaction field (Born matrix + CM5 shift) on the device"
                        : "SCF: device solvation upload failed; host-driven GFN1 loop with host solvation");
            }
            gfn1_resident = use_resident_loop;
            if (verb >= 2)
                CurcumaLogger::info(use_resident_loop
                    ? "SCF: fully device-resident GFN1 loop (host polls only dq + energy scalars)"
                    : "SCF: GFN1 resident-loop setup unavailable; host-driven device loop");
        }
    }

    pre_mark("GFN1 resident loop begin");
    // AP1 (Claude Generated): GPU partial diagonalisation window. The device density
    // only needs the occupied(+buffer) eigenvectors, so the resident eigensolve can
    // compute just the lowest gpu_n_eig pairs (cusolverDnDsyevdx) instead of the full
    // nao spectrum. OPT-IN (m_gpu_partial_diag, -gpu_partial_diag): measured
    // net-neutral on an RTX 5080 (the tridiagonalization, not the eigenvector count,
    // dominates the dense eigensolve), so the default GPU path stays the full solve.
    // Disabled when virtuals are genuinely required: the explicit Mulliken-CPSCF D4
    // charge response (d4_charge_source=cpscf) and the verbosity>=3 orbital listing both
    // read virtual orbitals. (The default variational D4 q-response needs no virtuals.)
    // The buffer covers the LUMO plus a few-kT Fermi tail; if the occupied tail ever
    // reaches the window edge (tiny-gap/metallic), the loop widens to the full solve once
    // so the 1e-8 energy is never at risk. Disabled on the device-potential path (its
    // re-solve would need the host D4 weights again).
    int gpu_n_eig = 0;  // 0 = full spectrum
    bool gpu_partial_diag = m_gpu_partial_diag && use_gpu_resident && !use_device_potential
        && (m_d4_charge_source != "cpscf") && (verb < 3);
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

    // Deferred host multipole integrals are needed as soon as any iteration runs on the host.
    if (!use_resident_loop)
        ensureHostMultipoleIntegrals();

    // Claude Generated (Oct 2026): Broyden seed preconditioner (xtb_precond.cpp). The scalar seed
    // alpha*F is replaced by M*F with M ~ -J^-1: -scf_precond frag builds M block-wise from the
    // isolated fragments, CURCUMA_SCF_PRECOND_FILE (developer override, wins) reads a dense M.
    // Set after every (re)construction of m_broyden, CPU loop only. The guard (ScfPrecondGuard,
    // -scf_precond_guard / CURCUMA_SCF_PRECOND_GUARD) falls back to the scalar step when the
    // residual does not drop. scf_precond none and no file: nothing is set, the path is unchanged.
    //
    // GPU (Claude Generated, Oct 2026): the fully device-resident CUDA loop mixes on the device,
    // so there M is uploaded once per build (flat blocks) and applied by the device Broyden; the
    // guard runs there on ||F||_2 with the host-side rule (ScfPrecondGuard::firesNorm) as a
    // callback. The resident ROCm loop does the same (Oct 24, 2026, WP5; HipScfBackend maps the indices
    // back to its packSCC order). Every other GPU path (host-driven device loops: ROCm, Vulkan, CUDA
    // without the resident loop) and a multi-GPU SCF split ignore the switch with one warning. Not counted
    // as a GPU fallback (nothing moves to the CPU), so -gpu_strict does not stop on it.
    pre_mark("partial diagonalisation, deferred host integrals");
    ScfPrecondGuard precond_guard;
    struct DeviceGuardState { ScfPrecondGuard g; int fired_at = -1; double fnorm = 0.0; };
    std::shared_ptr<DeviceGuardState> device_guard;   // set only when M runs on the device
    {
        const bool cpu_loop = !use_gpu_resident && !use_resident_loop && !use_device_potential;
        bool device_loop = false;
        std::string no_device_reason;
        if (use_resident_loop) {
            if (!m_gpu_scf->supportsResidentPreconditioner())
                no_device_reason = "this GPU backend has no device preconditioner";
            else if (m_gpu_scf->multiGpuScfSplit(m_basis.nao))
                no_device_reason = "a multi-GPU SCF split (eigensolver or density over several devices) is active";
            else
                device_loop = true;
        } else if (!cpu_loop) {
            no_device_reason = "host-driven GPU SCF path (ROCm, Vulkan, or CUDA without the fully "
                               "device-resident loop); implemented for the CPU loop and the resident CUDA and ROCm loops";
        }
        std::function<Vector(const Vector&)> pc_hook;
        std::string pc_err;
        if (!selectScfPreconditioner(mode == ScfMode::Broyden, cpu_loop, device_loop, no_device_reason,
                                     packed_len, pc_hook, precond_guard, pc_err)) {
            m_scf_converged = false; m_scf_iterations = 0;
            setHardError("SCF preconditioner file unusable: " + pc_err);
            return m_E_total;
        }
        if (pc_hook && device_loop) {
            pc_hook = nullptr;   // the device applies M itself
            const auto t_up0 = clock::now();
            std::string why;
            if (m_scf_precond_dev_src != m_scf_precond_frag) {   // new M: new process-unique tag
                static std::atomic<unsigned long long> next_tag{ 0 };
                m_scf_precond_dev_src = m_scf_precond_frag;
                m_scf_precond_dev_tag = ++next_tag;
            }
            GpuPrecondFlat flat;
            const bool resident = (m_gpu_scf->residentPreconditionerTag() == m_scf_precond_dev_tag);
            if (!resident) {
                flat.N = m_scf_precond_frag->dim();
                if (!m_scf_precond_frag->flatten(flat.idx_off, flat.idx, flat.b_off, flat.B,
                                                 flat.fallback, flat.max_dim)) {
                    why = "M has weighted blocks or more than 2^31 entries";
                } else {
                    // M is built on the host packing (packSCC: dp at nsh + k*nat + atom, qp at
                    // nsh + 3*nat + k*nat + atom); the device packs the moments atom-major
                    // (Eigen 3 x nat / 6 x nat storage: nsh + 3*atom + k, nsh + 3*nat + 6*atom + k).
                    auto to_device = [nsh, nat](int j) {
                        if (j < nsh) return j;
                        if (j < nsh + 3 * nat) { const int r = j - nsh; return nsh + 3 * (r % nat) + r / nat; }
                        const int r = j - nsh - 3 * nat;
                        return nsh + 3 * nat + 6 * (r % nat) + r / nat;
                    };
                    for (int& j : flat.idx) j = to_device(j);
                    std::vector<double> fb(flat.fallback.size());
                    for (int j = 0; j < static_cast<int>(fb.size()); ++j) fb[to_device(j)] = flat.fallback[j];
                    flat.fallback.swap(fb);
                }
            }
            auto st = std::make_shared<DeviceGuardState>();
            st->g = precond_guard;
            GpuBroydenGuardFn guard_fn;
            if (precond_guard.rule != ScfPrecondGuard::Rule::Off)
                guard_fn = [st](int n, double f, bool* store_best) {
                    bool new_best = false;
                    const bool hit = st->g.firesNorm(n, f, &new_best);
                    if (store_best) *store_best = new_best;
                    if (hit) { st->fired_at = n; st->fnorm = f; }
                    return hit;
                };
            if (why.empty()
                && m_gpu_scf->residentSetPreconditioner(resident ? nullptr : &flat, m_scf_precond_dev_tag,
                                                        guard_fn, precond_guard.keep, precond_guard.best)) {
                device_guard = st;
                const char* guard_name = precond_guard.rule == ScfPrecondGuard::Rule::Off
                                             ? "off" : precond_guard.name.c_str();
                if (!resident && verb >= 2)
                    CurcumaLogger::info_fmt("SCF precond: device M uploaded to the resident GPU loop ({} blocks, "
                                            "largest {}, {:.3f} MB, {:.1f} ms), guard '{}'",
                                            m_scf_precond_frag->blockCount(), flat.max_dim,
                                            8.0 * static_cast<double>(flat.B.size() + flat.fallback.size()) / 1048576.0,
                                            ms(t_up0, clock::now()), guard_name);
                else if (resident && verb >= 3)
                    CurcumaLogger::info_fmt("SCF precond: device M reused (resident), guard '{}'", guard_name);
            } else {
                if (why.empty()) why = m_gpu_scf->lastError();
                if (!m_scf_precond_warned)
                    CurcumaLogger::warn("scf_precond 'frag' ignored: the device preconditioner could not be set up"
                                        + (why.empty() ? std::string() : " (" + why + ")")
                                        + "; scalar Broyden seed (scf_damping) is used");
                m_scf_precond_warned = true;
            }
        } else if (pc_hook) {
            broyden.setStepHook(std::move(pc_hook));
        } else if (broyden.hasStepHook()) {   // only with -keep_diis: a hook left from an earlier SCF
            broyden.dropStepHook(true);
        }
    }

    double resident_band = 0.0;   // Tr(P H0) of the last device-resident step
    double resident_ecoul = 0.0, resident_ethird = 0.0, resident_emp = 0.0;  // its SCC energies
    pre_mark("preconditioner, SCF state");
    const auto t_scf_start = clock::now();
    // Claude Generated (Oct 2026): only UpdateMolecule cleared the flag, so a second Calculation on
    // the same geometry (the bare-H0 retry after implausible charges) inherited "converged" from the
    // first SCF and an SCF that ran out of iterations was accepted with exit code 0.
    m_scf_converged = false;
    m_pseudo_c_valid = false; m_pseudo_consecutive = 0; m_pseudo_used = false; m_eig_pseudo = false;   // Oct 2026
    m_pseudo_taken = m_pseudo_rejected = 0;
    int iter;
    if (!use_resident_loop) ensureHostIntegrals();   // the host-driven loops read m_S / m_H0 / m_X
    for (iter = 0; iter < max_iter; ++iter) {
        const auto t_iter0 = clock::now();

        // Stage 6 (S6.5): fully device-resident step. One backend call runs the
        // whole iteration on the GPU; only the convergence dq + the 4 SCC energy
        // scalars come back. The host keeps the iteration counter, the mixed-
        // precision (FP32→FP64) decision, the convergence test, and the verbosity
        // line — the rest (incl. the Broyden mix) is device-resident. Claude Generated.
        if (use_resident_loop) {
            m_eig_fp32 = fp32_wanted(dq_prev);
            double dq = 0.0, eb = 0.0, ecoul = 0.0, ethird = 0.0, emp = 0.0;
            if (!m_gpu_scf->residentScfStep(m_eig_fp32, dq, eb, ecoul, ethird, emp)) {
                std::string why = m_gpu_scf->lastError();
                // Claude Generated (Oct 2026, ROCm WP2): a backend that allows it restarts this SCF from the
                // unchanged initial guess in m_wfn on the host-driven loop (the resident steps never wrote
                // m_wfn or the host mixer); default backends keep the hard error below.
                if (m_gpu_scf->residentLoopFallback("SCF step " + std::to_string(iter)
                                                    + (why.empty() ? std::string() : " (" + why + ")"))) {
                    use_resident_loop = false;
                    ensureHostMultipoleIntegrals();
                    ensureHostIntegrals();
                    dq_prev = 1.0e30; e_total_old = 0.0; q_sh_old = m_wfn.q_sh;
                    dq_best_fp32 = 1.0e30; fp32_stall = 0; fp32_exhausted = false; dq_fp32_last = -1.0;
                    if (verb >= 2)
                        CurcumaLogger::warn("SCF: device-resident step failed at iteration " + std::to_string(iter)
                                            + "; restarting this SCF on the host-driven loop");
                    iter = -1;
                    continue;
                }
                // The classic cause on a large system: the stored multipole matrices fit, and
                // then the loop's own buffers no longer do. Name the flag that avoids it -
                // measured on polymer_2x (nat 7320, 7.7 GB stored) on a 20 GB card, where the
                // failure used to arrive with no indication of the cause. Claude Generated.
                if (iter == 0 && mp_stored_on_device && mp_stored_bytes > 1.0e9) {
                    if (!why.empty()) why += "; ";
                    why += "the stored GFN2 multipole interaction matrices hold "
                        + fmt::format("{:.1f}", mp_stored_bytes / 1073741824.0)
                        + " GB of device memory - rebuild them per iteration instead "
                          "(-gpu_multipole_otf auto, the default above ~2700 atoms)";
                }
                CurcumaLogger::warn("XTB::Calculation: GPU resident SCF step failed at iteration "
                                    + std::to_string(iter)
                                    + (why.empty() ? std::string() : ": " + why));
                m_scf_converged = false; m_scf_iterations = iter;
                setHardError("native GPU-resident SCF solve failed (see warning above) at iteration " + std::to_string(iter));
                return m_E_total;
            }
            m_E_electronic = eb; m_E_coulomb_shell = ecoul;
            resident_band = eb; resident_ecoul = ecoul; resident_ethird = ethird; resident_emp = emp;
            m_E_third_order = ethird; m_E_multipole = emp;
            const double e_scc = eb + ecoul + ethird + emp;
            const double de = (iter > 0) ? std::fabs(e_scc - e_total_old) : 0.0;
            // A non-finite SCC state is a failure, never convergence (Sep 2026, see the CPU loop).
            if (!std::isfinite(dq) || !std::isfinite(e_scc)) {
                CurcumaLogger::warn_fmt("XTB::Calculation: SCF produced non-finite values at iteration {} "
                                        "(E = {}, max|dq| = {}) - diverged", iter, e_scc, dq);
                m_scf_converged = false; m_scf_iterations = iter + 1;
                setHardError(fmt::format("native SCF diverged: non-finite energy/charges at iteration {}", iter));
                return m_E_total;
            }
            note_fp32_progress(dq);
            dq_prev = dq;
            if (verb >= scf_min) {
                const double t_iter_ms = ms(t_iter0, clock::now());
                CurcumaLogger::result(fmt::format("  {:4d}   {:18.10f}   {:10.2e}   {:10.2e}   {:6.1f}",
                                                  iter, e_scc, de, dq, t_iter_ms));
            }
            if (iter > 0 && !m_eig_fp32 && dq < thresh && de < thresh * 100.0
                && (!xlbomd_corrector_step || iter + 1 >= m_scf_xlbomd_correctors)) {
                m_scf_converged = true;
                e_total_old = e_scc;
                break;
            }
            // Claude Generated (Oct 2026): the device guard decided inside this step's Broyden
            // update; report it like the CPU loop (which never checks a converged iteration).
            if (device_guard && device_guard->fired_at == iter) {
                const ScfPrecondGuard& g = device_guard->g;
                CurcumaLogger::warn_fmt("SCF preconditioner guard '{}' fired at iteration {} "
                                        "(||F||_2 = {:.3e}, ||F_0||_2 = {:.3e}): scalar Broyden step "
                                        "for the rest of this SCF ({} history, from {} x; device)",
                                        g.name, iter, device_guard->fnorm, g.f0,
                                        g.keep ? "kept" : "reset", g.best ? "best" : "current");
                if (std::getenv("CURCUMA_SCF_TRACE"))
                    fmt::print(stderr, "SCFTRACE_GUARD iter {} rule {} device\n", iter, g.name);
            }
            e_total_old = e_scc;
            continue;
        }

        // Capture the SCC input vector x_in (current m_wfn charges, and for GFN2 the
        // atomic multipole moments) before the potential is built and the populations
        // are overwritten by the diagonalisation. Broyden mixes it; the convergence
        // test below measures the residual packSCC() - x_in on the SAME vector, which
        // is what tblite does (`pconverged = mixer%get_error() < pconv`,
        // xtb/singlepoint.f90:255). Captured for every mode, not just Broyden.
        Vector x_in = packSCC();

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
                // -xtb.rocm_device_multipole: same assignments as addMultipolePotential
                // (v_at replaced, v_dp/v_qp set), computed on the device. Claude Generated (Oct 2026).
                bool mp_on_dev = false;
                if (use_dev_mp) {
                    mp_on_dev = m_gpu_scf->multipolePotential(m_wfn.q_at, m_wfn.dp_at, m_wfn.qp_at,
                                                              m_pot.v_at, m_pot.v_dp, m_pot.v_qp);
                    if (!mp_on_dev) dev_mp_failed("multipole potential");
                }
                if (!mp_on_dev)
                    addMultipolePotential(m_pot);
                // Self-consistent D4: exact per-reference dE_D4/dq into the atom potential.
                const auto t_disp0 = clock::now();
                addDispersionPotential(m_pot);
                acc_disp += ms(t_disp0, clock::now());
            }

            // Implicit solvation: v_at += (B * q_at). Folded into v_ao by
            // expand_potential below, so it enters the Fock (and the host->GPU
            // upload) automatically. Charge-independent SASA/shift do not appear.
            addSolvationPotential(m_pot);
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
            m_eig_fp32 = fp32_wanted(dq_prev);
            bool solve_ok;
            Eigen::VectorXd v_ao;  // host-expanded potential (non-device-potential path)
            // GFN2 multipole solve. With -xtb.rocm_device_multipole the device already holds
            // v_dp/v_qp (no upload); if that solve fails, the host copies in m_pot are uploaded
            // as before. Claude Generated (Oct 2026).
            auto solve_mp = [&](const Eigen::VectorXd& vao) -> bool {
                if (use_dev_mp) {
                    if (m_gpu_scf->solveMultipoleResident(vao, m_wfn.eps, m_eig_fp32)) return true;
                    dev_mp_failed("resident multipole solve");
                }
                return m_gpu_scf->solveMultipole(vao, m_pot.v_dp, m_pot.v_qp, m_wfn.eps, m_eig_fp32, gpu_n_eig);
            };
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
                // Claude Generated (Oct 2026, ROCm WP1): a backend that allows it continues this SCF with the
                // host potential build (same input charges/moments); default backends keep the hard error.
                if (!solve_ok && m_gpu_scf->devicePotentialFallback("SCF step " + std::to_string(iter))) {
                    use_device_potential = false;
                    m_pot.reset();
                    addCoulombShellPotential(m_pot);
                    addThirdOrderPotential(m_pot);
                    addMultipolePotential(m_pot);
                    addDispersionPotential(m_pot);
                    addSolvationPotential(m_pot);
                    v_ao.resize(m_basis.nao);
                    for (int ao = 0; ao < m_basis.nao; ++ao)
                        v_ao(ao) = m_pot.v_sh(m_basis.ao2sh[ao]) + m_pot.v_at(m_basis.ao2at[ao]);
                    solve_ok = solve_mp(v_ao);
                    if (verb >= 2)
                        CurcumaLogger::warn("SCF: device potential build failed at iteration " + std::to_string(iter)
                                            + "; host potential build for the rest of this SCF");
                }
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
                    ? solve_mp(v_ao)
                    : m_gpu_scf->solve(v_ao, m_wfn.eps, m_eig_fp32, gpu_n_eig);
            }
            if (!solve_ok) {
                CurcumaLogger::warn("XTB::Calculation: GPU resident solve failed at iteration "
                                    + std::to_string(iter));
                m_scf_converged = false; m_scf_iterations = iter;
                setHardError("native GPU-resident SCF solve failed (see warning above) at iteration " + std::to_string(iter));
                return m_E_total;
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
                    ? solve_mp(v_ao)
                    : m_gpu_scf->solve(v_ao, m_wfn.eps, m_eig_fp32, gpu_n_eig);
                if (!re_ok) {
                    CurcumaLogger::warn("XTB::Calculation: GPU resident full re-solve failed at "
                                        "iteration " + std::to_string(iter));
                    m_scf_converged = false; m_scf_iterations = iter;
                setHardError("native GPU-resident SCF solve failed (see warning above) at iteration " + std::to_string(iter));
                return m_E_total;
                }
                occupationsFromEps(m_wfn.eps, occ, ncol);
            }
            // Device: P = C·diag(occ)·Cᵀ; pop_ao(μ)=Σ_ν P_μν·S_μν; band=Σ P⊙H0.
            Eigen::VectorXd pop_ao;
            if (!m_gpu_scf->density(occ, ncol, pop_ao, gpu_band)) {
                CurcumaLogger::warn("XTB::Calculation: GPU resident density failed at iteration "
                                    + std::to_string(iter));
                m_scf_converged = false; m_scf_iterations = iter;
                setHardError("native GPU-resident SCF solve failed (see warning above) at iteration " + std::to_string(iter));
                return m_E_total;
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
                m_scf_converged = false; m_scf_iterations = iter;
                setHardError("native GPU-resident SCF solve failed (see warning above) at iteration " + std::to_string(iter));
                return m_E_total;
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
            m_eig_fp32 = fp32_wanted(dq_prev);
            // Oct 2026 (opt-in -scf_pseudo_diag_cpu): a pseudo-diagonalisation step may replace the
            // eigensolve while the residual is below scf_pseudo_diag_max_dq and still above 10x the
            // convergence threshold, so the last iterations are full solves (exact eigenpairs for the gradient).
            m_eig_pseudo = m_pseudo_cpu && !scf_jac_active && iter > 0
                           && dq_prev < m_pseudo_max_dq && dq_prev > 10.0 * thresh;

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
                setHardError("native SCF eigensolver failed at iteration " + std::to_string(iter)
                             + " (often a singular overlap, e.g. an unsupported basis shell)");
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

        // Energies for this iteration (Claude Generated, Sep 2026: timed individually, the
        // "energy/mix" bucket used to be one opaque number).
        const auto t_e0 = clock::now();
        m_E_coulomb_shell = energyCoulombShell();
        const auto t_e1 = clock::now();
        m_E_third_order   = energyThirdOrder();
        const auto t_e2 = clock::now();
        // -xtb.rocm_device_multipole: E_mp from the device copy of the interaction
        // matrices, host energyMultipole() on failure. Claude Generated (Oct 2026).
        {
            double e_mp_dev = 0.0;
            if (use_dev_mp && m_gpu_scf->multipoleEnergy(m_wfn.q_at, m_wfn.dp_at, m_wfn.qp_at, e_mp_dev)) {
                m_E_multipole = e_mp_dev;
            } else {
                if (use_dev_mp) dev_mp_failed("multipole energy");
                m_E_multipole = energyMultipole();
            }
        }
        const auto t_e3 = clock::now();

        // Band energy: Tr(P · H0). Device-resident path returns it from the GPU
        // (P is not on the host during the loop); CPU path sums P⊙H0 directly.
        m_E_electronic = use_gpu_resident ? gpu_band
                                          : (m_wfn.P.cwiseProduct(m_H0)).sum();

        const auto t_e4 = clock::now();
        acc_e_coul  += ms(t_e0, t_e1);
        acc_e_third += ms(t_e1, t_e2);
        acc_e_mp    += ms(t_e2, t_e3);
        acc_e_band  += ms(t_e3, t_e4);

        // Total SCC = band + coulomb + third-order + multipole
        const double e_scc = m_E_electronic + m_E_coulomb_shell
                           + m_E_third_order + m_E_multipole;

        // Accumulate per-phase timings (see breakdown at verbosity >= 3).
        acc_pot    += ms(t_iter0, t_pot);
        acc_fock   += ms(t_pot,   t_fock);
        acc_solve  += ms(t_fock,  t_solve);
        acc_mull   += ms(t_solve, t_mull);
        acc_energy += ms(t_mull,  clock::now());

        // Per-iteration diagnostics. dq is the residual of the FULL SCC vector the
        // mixer works on: shell charges for GFN1, shell charges + atomic dipole and
        // quadrupole moments for GFN2. Measuring only the shell charges (as this did
        // until Sep 2026) lets GFN2 stop while the moments are still moving — the
        // energy is then still right to O(dP^2) but the GRADIENT only to O(dP).
        // H2 is the extreme case: symmetry pins its shell charges to 0 from the first
        // iteration, so the SCF exited after 2 cycles with the moments unconverged and
        // the analytic gradient came out 0.0330 instead of 0.0211 Eh/Bohr (60 % off,
        // while the energy still matched xtb to 1e-8). Claude Generated.
        const Vector x_out_check = packSCC();
        const double dq = (x_out_check - x_in).cwiseAbs().maxCoeff();
        // Claude Generated (Sep 2026): a non-finite SCC state is a failure, not convergence. With
        // -scf_mode diis on a 1846-atom water/urea cluster the charges turned NaN at iteration 6;
        // maxCoeff() over the NaN vector then came out 0, dE too, and the SCF "converged" to
        // +540710 Eh, which was returned as a valid single point. Test the vectors themselves,
        // not only dq.
        if (!std::isfinite(dq) || !std::isfinite(e_scc) || !x_out_check.allFinite() || !x_in.allFinite()) {
            CurcumaLogger::warn_fmt("XTB::Calculation: SCF produced non-finite values at iteration {} "
                                    "(E = {}) - diverged", iter, e_scc);
            m_scf_converged = false; m_scf_iterations = iter + 1;
            setHardError(fmt::format("native SCF diverged: non-finite energy/charges at iteration {}", iter));
            return m_E_total;
        }
        // Claude Generated (Oct 2026): CURCUMA_SCF_TRACE=1 prints, per CPU SCF iteration, the Broyden
        // residual F = x_out - x_in split into the shell-charge block (q_sh) and the GFN2 multipole
        // block (dp, qp) as max|F| and ||F||_2, plus the phase times in ms, to stderr. Read-only;
        // no numerical path depends on it.
        static const bool scf_trace = (std::getenv("CURCUMA_SCF_TRACE") != nullptr);
        if (scf_trace) {
            const Vector Fres = x_out_check - x_in;
            const Eigen::Index nq = nsh;
            const Eigen::Index nm = Fres.size() - nq;
            const double q_max = Fres.head(nq).cwiseAbs().maxCoeff();
            const double q_nrm = Fres.head(nq).norm();
            const double m_max = (nm > 0) ? Fres.tail(nm).cwiseAbs().maxCoeff() : 0.0;
            const double m_nrm = (nm > 0) ? Fres.tail(nm).norm() : 0.0;
            fmt::print(stderr, "SCFTRACE iter {} nat {} qmax {:.6e} qnorm {:.6e} mpmax {:.6e} mpnorm {:.6e} "
                               "t_pot {:.3f} t_fock {:.3f} t_solve {:.3f} t_mull {:.3f} t_energy {:.3f} t_iter {:.3f}\n",
                       iter, nat, q_max, q_nrm, m_max, m_nrm,
                       ms(t_iter0, t_pot), ms(t_pot, t_fock), ms(t_fock, t_solve),
                       ms(t_solve, t_mull), ms(t_mull, clock::now()), ms(t_iter0, clock::now()));
        }
        note_fp32_progress(dq);
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
        // XL-BOMD corrector: optionally require a minimum number of corrector SCF
        // cycles (scf_xlbomd_correctors) before accepting convergence, polishing the
        // density beyond the loose threshold. Default 1 = no effect. Claude Generated.
        if (iter > 0 && !m_eig_fp32 && !m_pseudo_used
            && (!xlbomd_corrector_step || iter + 1 >= m_scf_xlbomd_correctors)) {
            // dq is the full-SCC-vector residual computed above (shell charges for
            // GFN1, + atomic multipole moments for GFN2), matching what the mixer
            // acts on and what tblite converges on.
            if (checkConvergence_impl(dq, e_total_old, e_scc, thresh)) {
                m_scf_converged = true;
                if (scf_jac_active) {   // keep x* and G(x*) for the dump / capture below
                    scf_jac_x_star = x_in;
                    scf_jac_g_star = x_out_check;
                    scf_jac_e_scc  = e_scc;
                }
                break;
            }
        }

        // Broyden: mix the SCC charge vector to get the next input. The raw
        // output x_out (in m_wfn after updatePopulations) is combined with the
        // input x_in via the modified-Broyden quasi-Newton update; unpacking
        // writes the next input back into m_wfn so q_sh_old below tracks it.
        if (mode == ScfMode::Broyden) {
            const auto t_mix0 = clock::now();
            const Vector x_out  = packSCC();
            // Claude Generated (Oct 2026): preconditioner guard. Uses only the residual norm of
            // this iteration; on a hit M is dropped for the rest of this SCF (scalar step from the
            // current x, history reset or kept). Off (default) or no hook: nothing is computed.
            bool restart_from_best = false;
            if (precond_guard.rule != ScfPrecondGuard::Rule::Off && broyden.hasStepHook()) {
                const Vector Fres = x_out - x_in;
                const double fnorm = Fres.norm();
                if (precond_guard.fires(iter, x_in, Fres)) {
                    broyden.dropStepHook(precond_guard.keep);
                    restart_from_best = precond_guard.best;
                    CurcumaLogger::warn_fmt("SCF preconditioner guard '{}' fired at iteration {} "
                                            "(||F||_2 = {:.3e}, ||F_0||_2 = {:.3e}): scalar Broyden step "
                                            "for the rest of this SCF ({} history, from {} x)",
                                            precond_guard.name, iter, fnorm, precond_guard.f0,
                                            precond_guard.keep ? "kept" : "reset",
                                            precond_guard.best ? "best" : "current");
                    if (scf_trace)
                        fmt::print(stderr, "SCFTRACE_GUARD iter {} rule {}\n", iter, precond_guard.name);
                }
            }
            // With "best", the scalar step is taken from the stored best input x_k and its output
            // x_k + F_k (both from an earlier evaluation of the map).
            const Vector x_next = restart_from_best
                ? broyden.update(precond_guard.x_best, precond_guard.x_best + precond_guard.F_best)
                : broyden.update(x_in, x_out);
            restart_from_best = false;
            unpackSCC(x_next);
            if (scf_trace)
                fmt::print(stderr, "SCFTRACE_MIX iter {} t_mix {:.3f}\n", iter, ms(t_mix0, clock::now()));
        }

        q_sh_old = m_wfn.q_sh;
        e_total_old = e_scc;
    }

    m_scf_iterations = iter + 1;
    const auto t_scf_end = clock::now();

    // XL-BOMD (Phase 2): advance / seed the auxiliary SCC trajectory from the CONVERGED
    // density D[P_n] (= current m_wfn). Claude Generated.
    //  - Corrector step: propagate the auxiliary with the time-reversible Niklasson
    //    dissipative integrator P_{n+1}=2P_n-P_{n-1}+kappa(D-P_n)+alpha Σ c_k P_{n-k}.
    //  - Bootstrap step (aux not yet full): seed the trajectory with the converged density.
    if (xlbomd_active && m_scf_converged) {
        const Vector D_n = packSccState();
        if (D_n.size() == packed_len) {
            if (xlbomd_corrector_step) {
                double kappa = 0.0, alpha = 0.0;
                std::vector<double> cc;
                if (xlbomdCoefficients(xlbomd_K, kappa, alpha, cc)
                    && static_cast<int>(m_xlbomd_aux.size()) >= xlbomd_K + 1) {
                    Vector Pnext = 2.0 * m_xlbomd_aux[0] - m_xlbomd_aux[1]
                                 + kappa * (D_n - m_xlbomd_aux[0]);
                    for (int k = 0; k <= xlbomd_K; ++k)
                        Pnext.noalias() += (alpha * cc[k]) * m_xlbomd_aux[k];
                    if (Pnext.allFinite()) {
                        m_xlbomd_aux.push_front(std::move(Pnext));
                        while (static_cast<int>(m_xlbomd_aux.size()) > xlbomd_K + 1)
                            m_xlbomd_aux.pop_back();
                    }
                }
            } else {
                if (!m_xlbomd_aux.empty() && m_xlbomd_aux.front().size() != packed_len)
                    m_xlbomd_aux.clear();
                m_xlbomd_aux.push_front(D_n);
                while (static_cast<int>(m_xlbomd_aux.size()) > xlbomd_K + 1)
                    m_xlbomd_aux.pop_back();
            }
        }
    }

    if (verb >= scf_min) {
        if (m_scf_converged)
            CurcumaLogger::success_fmt("SCF converged in {} iterations ({:.1f} ms)",
                                       m_scf_iterations, ms(t_scf_start, t_scf_end));
        else
            CurcumaLogger::warn_fmt("SCF NOT converged after {} iterations ({:.1f} ms)",
                                    m_scf_iterations, ms(t_scf_start, t_scf_end));
    }
    if (m_pseudo_cpu && verb >= 2)
        CurcumaLogger::info_fmt("SCF: CPU pseudo-diagonalisation steps: {} taken, {} rejected or gated", m_pseudo_taken, m_pseudo_rejected);

    // Claude Generated (Sep 2026): an unconverged SCF is not a result. It used to be returned
    // like one - level shift on a 1846-atom water/urea cluster ended after 150 iterations at
    // +3080 Eh with exit code 0. xtb and tblite stop with an error here too. The explicit opt-out
    // keeps the old behaviour for anyone who knowingly accepts a nearly converged energy.
    if (!m_scf_converged && !m_scf_allow_unconverged) {
        setHardError(fmt::format("native SCF did not converge within {} iterations (last max|dq| = {:.2e}); "
                                 "-scf_allow_unconverged true accepts the unconverged energy",
                                 m_scf_iterations, dq_prev));
        return m_E_total;
    }

    // Claude Generated (Oct 2026): SCF Jacobian debug dump (CURCUMA_SCF_JAC=<dir>), before any
    // post-SCF step. evalG is one pure evaluation of the SCC fixed-point map x -> G(x) of the CPU
    // loop: unpack x -> potential (Coulomb, third order, GFN2 multipole + in-SCF D4, solvation)
    // -> Fock -> FP64 eigensolve with Fermi occupations -> Mulliken charges and moments -> pack.
    // No mixer, level shift, DIIS or FP32 state enters it. m_wfn / m_pot and the eigensolver
    // timing counters are restored afterwards, so the post-SCF path sees the converged state.
    if (scf_jac_active && m_scf_converged) {
        const bool cpu_loop = !use_gpu_resident && !use_resident_loop && !use_device_potential;
        const bool pure_eig = !m_external_eigensolver && m_eigensolver != "lobpcg"
                              && m_eigensolver != "purify";
        if (!cpu_loop || !pure_eig || mode != ScfMode::Broyden || scf_jac_x_star.size() == 0) {
            if (scf_jac_dump)
                CurcumaLogger::warn("CURCUMA_SCF_JAC: Jacobian dump skipped (needs the CPU SCF loop, "
                                    "-scf_mode broyden and a stateless eigensolver mkl/native)");
        } else {
            const Wavefunction wfn_save = m_wfn;
            const Potential    pot_save = m_pot;
            const bool   fp32_save = m_eig_fp32;
            const double t_save[5] = { m_t_xfx, m_t_diag, m_t_back, m_t_dens, m_t_xfx_copy };
            const int    c_save[3] = { m_eig_calls_native, m_eig_calls_fp32, m_eig_calls_lapack };
            bool eval_failed = false;
            auto evalG = [&](const Vector& x) -> Vector {
                unpackSCC(x);
                m_pot.reset();
                addCoulombShellPotential(m_pot);
                addThirdOrderPotential(m_pot);
                if (m_method == MethodType::GFN2) {
                    addMultipolePotential(m_pot);
                    addDispersionPotential(m_pot);
                }
                addSolvationPotential(m_pot);
                const Matrix F = buildFock(m_H0, m_S, m_pot);
                m_eig_fp32 = false;   // FP64 only in this mode
                bool ok;
                {
                    MklThreadScope eig_scope(eig_threads);
                    ok = solveEigen(F, m_S);
                }
                if (!ok) { eval_failed = true; return Vector::Constant(x.size(), std::nan("")); }
                updatePopulations(m_S);
                return packSCC();
            };
            if (scf_jac_dump) {
                writeScfJacobianDump(scf_jac_dir, scf_jac_x_star, scf_jac_g_star, scf_jac_e_scc, evalG);
                if (eval_failed)
                    CurcumaLogger::warn("CURCUMA_SCF_JAC: an eigensolve failed during the Jacobian dump");
            } else if (m_scf_jac_capture) {
                // Fragment of -scf_precond frag: J at x* by central differences, h = 1e-4 (the
                // default of the dump, with which the prototype measurements were made), or by linear
                // response (-scf_precond_jacobian response, Oct 2026). CURCUMA_SCF_JAC_COMPARE=1
                // (developer switch) computes both and prints their deviation on stderr.
                const bool use_resp = (m_scf_precond_jacobian == "response");
                const bool compare = std::getenv("CURCUMA_SCF_JAC_COMPARE") != nullptr;
                const int dim = static_cast<int>(scf_jac_x_star.size());
                PrecondMatrix J_fd, J_resp;
                double t_fd = 0.0, t_resp = 0.0;
                if (!use_resp || compare) {
                    const auto t0 = clock::now();
                    std::vector<double> J;
                    numericalScfJacobian(scf_jac_x_star, scf_jac_g_star, evalG, 1.0e-4, J, nullptr);
                    t_fd = ms(t0, clock::now());
                    if (!eval_failed) J_fd = Eigen::Map<const PrecondMatrix>(J.data(), dim, dim);
                }
                if (use_resp || compare) {
                    const auto t0 = clock::now();
                    std::string rerr;
                    bool ok;
                    {
                        MklThreadScope eig_scope(eig_threads);
                        ok = computeScfJacobianResponse(scf_jac_x_star, J_resp, rerr);
                    }
                    t_resp = ms(t0, clock::now());
                    if (!ok) J_resp = PrecondMatrix();
                }
                m_scf_jac_result = use_resp ? J_resp : J_fd;
                if (compare && J_fd.rows() == dim && J_resp.rows() == dim) {
                    const PrecondMatrix Bf = Eigen::PartialPivLU<PrecondMatrix>(J_fd).inverse();
                    const PrecondMatrix Br = Eigen::PartialPivLU<PrecondMatrix>(J_resp).inverse();
                    fmt::print(stderr, "SCFJACCMP nat {} dim {} maxabs {:.3e} relfro {:.3e} inv_relfro {:.3e} t_fd_ms {:.1f} t_resp_ms {:.1f}\n",
                               m_atomcount, dim, (J_resp - J_fd).cwiseAbs().maxCoeff(),
                               (J_resp - J_fd).norm() / J_fd.norm(), (Br - Bf).norm() / Bf.norm(), t_fd, t_resp);
                }
            }
            m_wfn = wfn_save;
            m_pot = pot_save;
            m_eig_fp32 = fp32_save;
            m_t_xfx = t_save[0]; m_t_diag = t_save[1]; m_t_back = t_save[2];
            m_t_dens = t_save[3]; m_t_xfx_copy = t_save[4];
            m_eig_calls_native = c_save[0]; m_eig_calls_fp32 = c_save[1]; m_eig_calls_lapack = c_save[2];
        }
    }

    // Claude Generated (Sep 2026): post-SCF host phase timings (printed with CURCUMA_GPU_PROFILE
    // or verbosity 3) - on a 7k-atom GPU run this block took 40 s.
    std::vector<std::pair<const char*, double>> post_t;
    auto post_mark = [&, tp = clock::now()](const char* name) mutable {
        const auto now = clock::now();
        post_t.emplace_back(name, ms(tp, now));
        tp = now;
    };

    // Stage 6 (S6.5): the fused device loop kept the converged charges/moments
    // resident — download them once into the host wavefunction so the post-SCF
    // potential rebuild + energies + gradient below run unchanged. Claude Generated.
    if (use_resident_loop && m_gpu_scf) {
        m_wfn.q_sh.resize(m_basis.nsh);
        m_wfn.q_at.resize(m_atomcount);
        m_wfn.dp_at = Eigen::MatrixXd::Zero(3, m_atomcount);
        m_wfn.qp_at = Eigen::MatrixXd::Zero(6, m_atomcount);
        m_wfn.eps.resize(m_basis.nao);
        m_gpu_scf->residentLoopCharges(m_wfn.q_sh, m_wfn.q_at, m_wfn.dp_at, m_wfn.qp_at, m_wfn.eps);
    }

    // Device-resident path: download the converged density and MO coefficients
    // once, so the post-SCF energies (band Tr(P·H0)) and the (still CPU) gradient
    // read them from the host. Claude Generated, GPU port Stage 2.
    post_mark("download charges/moments/eps");
    // Claude Generated (Sep 2026): a large single point on the device-resident path does not
    // need the dense P and C on the host (band energy comes from the device, the gradient is
    // on the device). Rebuilding and downloading them took 16 s at 7320 atoms; they are fetched
    // on demand (host gradient fallback, property accessors after ensureHostWavefunction).
    //
    // Point 1 (Claude Generated, Sep 2026): the above previously forced finalize() whenever a
    // gradient was requested, even though calculateGradientGpu() (called below with the SAME
    // gate as device_grad_will_run) never reads m_wfn.P/C - it passes pc_resident=true and
    // reads dP/dC straight off the device (xtb_gradient.cpp:859-868). Skip finalize() there too;
    // ensureHostWavefunction() downloads on demand for the few paths that still need the host
    // matrices (host gradient fallback, property accessors) - see xtb_multipole.cpp:44.
    const bool device_grad_will_run = gradient && use_gpu_resident && m_gpu_scf
        && m_gpu_scf->supportsGradient();
    m_wfn_on_device = false;
    if (use_gpu_resident && m_gpu_scf) {
        // The converged P/C stay on the device only for big systems (as before for GFN2: dense multipole integrals
        // above 2 GB, nao > ~5270); the small ones keep the public accessors (getDensity, MolecularOrbitals)
        // filled. The multipole integrals themselves are deferred for every size (see above), which does not
        // need P/C on the host.
        const bool big_wfn = 9.0 * static_cast<double>(m_basis.nao) * m_basis.nao * sizeof(double) > 2.0e9;
        if (use_resident_loop && big_wfn && (m_mp_ints_deferred || gfn1_resident) && (!gradient || device_grad_will_run))
            m_wfn_on_device = true;
        else
            m_gpu_scf->finalize(m_wfn.P, m_wfn.C);
    }
    post_mark("finalize: download P and C");

    // GPU paths (device eigensolve and the fully-resident loop) keep the
    // occupations on-device: solveEigen is the only writer of m_wfn.focc and it
    // is skipped there. Rebuild focc on the host from the downloaded eps so the
    // electronic free-energy term below AND the gradient's energy-weighted
    // density W use the fractional Fermi occupations instead of the integer
    // fallback (xtb_gradient.cpp). occupationsFromEps is bit-faithful to
    // solveEigen, so this is identical to the CPU path. Claude Generated.
    if ((use_gpu_resident || use_resident_loop) && m_wfn.eps.size() == m_basis.nao) {
        Eigen::VectorXd occ; int ncol = 0;
        occupationsFromEps(m_wfn.eps, occ, ncol);
        m_wfn.focc = occ;
    }
    post_mark("host occupations");

    // Persist converged Fock matrix for gradient / debug. On the device-potential
    // path (Stage 5 B3/B4) the host m_pot was skipped during the loop, so rebuild
    // it once at the converged charges before forming m_F (the gradient block below
    // rebuilds it again, but m_F must be consistent here too).
    // (not needed when the wavefunction stays on the device: no host Fock, no host gradient)
    // Claude Generated (Oct 2026, ROCm WP8b): with deferHostCopies() (CUDA: false) and S/H0 still deferred, m_F is
    // skipped below, so this rebuild had no reader (the gradient rebuilds m_pot or takes the device potentials); it
    // would only fetch the host S/H0/L/gamma and build the 18 nat^2 multipole interaction matrices.
    const bool skip_post_pot = use_resident_loop && m_gpu_scf && m_gpu_scf->deferHostCopies()
        && m_host_ints_deferred && m_d4_charge_source != "cpscf";
    if ((use_device_potential || gfn1_resident) && !m_wfn_on_device && !skip_post_pot) {
        m_pot.reset();
        addCoulombShellPotential(m_pot);
        addThirdOrderPotential(m_pot);
        addMultipolePotential(m_pot);
        addDispersionPotential(m_pot);
        addSolvationPotential(m_pot);
    }
    post_mark("host potential rebuild");
    // Claude Generated (Sep 2026): m_F is only read by debug dumps / consistency audits. With
    // deferred host multipole integrals (large system on the GPU) skip the O(nao^2) host
    // rebuild instead of forcing the 17 GB integral build just for it.
    if (m_mp_ints_deferred || m_wfn_on_device || m_host_ints_deferred)
        m_F.resize(0, 0);
    else
        m_F = buildFock(m_H0, m_S, m_pot);
    post_mark("host Fock matrix");

    // Final energies
    // Claude Generated (Sep 2026): the device-resident step already evaluated these SCC energies
    // from the same converged output charges/moments the host would use (7320 atoms: 5 s of
    // O(nat^2) host work). Reused when the wavefunction stays on the device.
    if (m_wfn_on_device || use_resident_loop) {
        m_E_coulomb_shell = resident_ecoul;
        m_E_third_order   = resident_ethird;
        m_E_multipole     = resident_emp;
        m_E_electronic    = resident_ecoul + resident_ethird + resident_emp;
    } else {
        m_E_electronic    = energyCoulombShell() + energyThirdOrder() + energyMultipole();
    }
    post_mark("Coulomb/third-order/multipole energies");
    // Band energy: Tr(P · H0) for electronic part (from the device when P stayed there)
    m_E_electronic   += (m_wfn_on_device || use_resident_loop) ? resident_band : (m_wfn.P.cwiseProduct(m_H0)).sum();
    post_mark("band energy Tr(P H0)");
    // Electronic free-energy (Mermin/Fermi entropy) term g = -T*S. Zero for
    // gapped systems (occ ∈ {0,2}); only bites on small-gap systems with
    // fractional occupations. Folded into the electronic container so it matches
    // tblite/xtb, which report the free energy A = E - T*S. Claude Generated.
    m_E_entropy       = electronicFreeEnergy();
    m_E_electronic   += m_E_entropy;
    m_E_repulsion     = calcRepulsionEnergy();
    m_E_halogen_bond  = calcHalogenBondEnergy();
    post_mark("entropy/repulsion/halogen bond");
    // Exact GFN2 self-consistent-D4 q-response (F5): d4_charge_source="mulliken" (default)
    // handles dE_D4/dq·dq/dR variationally through the gradient charge-Pulay + W — dE_D4/dq
    // is folded into the gradient v_at (as tblite does in dispersion%get_potential), so the
    // −P·(v_a+v_b)·dS overlap-Pulay plus W = Σf·ε·ccᵀ (ε already carry the in-SCF D4 potential)
    // reproduce the full Mulliken charge response with no separate solve. "eeq"/"cpscf" fall
    // back to the explicit single-shot-EEQ / Z-vector response. Set BEFORE calcDispersionEnergy
    // (so its q-response branch is skipped) and BEFORE the gradient m_pot rebuild.
    m_d4_variational_qresp = (m_method == MethodType::GFN2)
        && (m_d4_charge_source == "mulliken")
        && !std::getenv("CURCUMA_D4_NONVARIATIONAL");
    m_E_dispersion    = calcDispersionEnergy(gradient);
    post_mark("dispersion (D4, incl. gradient prep)");
    // Implicit solvation free energy (Born + CDS + shift). Added separately to the
    // total — NOT into Tr(P·H0) — mirroring the Coulomb ES2 handling, so there is no
    // double counting (the in-SCF v_at only polarized the density). Claude Generated.
    m_E_solvation     = energySolvation();

    m_E_total = m_E_electronic + m_E_repulsion
              + m_E_halogen_bond + m_E_dispersion + m_E_solvation;
    post_mark("solvation + total");
    const auto t_energies = clock::now();
    if (verb >= 3 || std::getenv("CURCUMA_GPU_PROFILE")) {
        std::string rep = "Post-SCF host phases:\n";
        for (const auto& [name, t] : post_t)
            rep += fmt::format("  {:<40s} {:10.1f} ms\n", name, t);
        CurcumaLogger::result(rep);
    }

    if (verb >= 2) {
        CurcumaLogger::info("Energy decomposition:");
        CurcumaLogger::info(fmt::format("  Electronic   = {: .8f} Eh", m_E_electronic));
        CurcumaLogger::info(fmt::format("  Coulomb ES2  = {: .8f} Eh", m_E_coulomb_shell));
        CurcumaLogger::info(fmt::format("  Third-order  = {: .8f} Eh", m_E_third_order));
        if (m_method == MethodType::GFN2)
            CurcumaLogger::info(fmt::format("  Multipole    = {: .8f} Eh", m_E_multipole));
        CurcumaLogger::info(fmt::format("  Free E (-TS) = {: .8f} Eh (folded into Electronic)", m_E_entropy));
        CurcumaLogger::info(fmt::format("  Repulsion    = {: .8f} Eh", m_E_repulsion));
        if (m_method == MethodType::GFN1)
            CurcumaLogger::info(fmt::format("  Halogen bond = {: .8f} Eh", m_E_halogen_bond));
        CurcumaLogger::info(fmt::format("  Dispersion   = {: .8f} Eh", m_E_dispersion));
        if (m_solvation)
            CurcumaLogger::info(fmt::format("  Solvation    = {: .8f} Eh", m_E_solvation));
        CurcumaLogger::result(fmt::format("  Total        = {: .8f} Eh", m_E_total));

        // Frontier orbitals
        const double homo = getHOMOEnergy();
        const double lumo = getLUMOEnergy();
        CurcumaLogger::info(fmt::format("HOMO = {: .6f} Eh   LUMO = {: .6f} Eh   gap = {:.4f} eV",
                                        homo, lumo, (lumo - homo) * 27.211386245988));
    }

    // Mirror the converged wavefunction for the wrapper accessors
    if (!m_wfn_on_device) m_mo = m_wfn.C;
    else m_mo.resize(0, 0);
    m_energies = m_wfn.eps;
    m_num_electrons = static_cast<int>(m_wfn.nocc);
    ensure_host_cn();
    m_coordination_numbers = cn;

    const auto t_grad0 = clock::now();
    if (gradient) {
        // Rebuild potential with converged shell charges before computing gradient.
        // (m_pot was overwritten during SCF; final values must be regenerated.)
        const auto t_gp0 = clock::now();
        // Claude Generated (Oct 2026): on the resident loop the device builds the converged potentials itself
        // (convergedPotentials: gamma q_sh, third order, multipole, D4 dE/dq, solvation not supported) and the
        // gradient reads them from there; the host m_pot is only rebuilt for the host fallbacks below.
        // CURCUMA_GRAD_DEVPOT=0 keeps the host rebuild.
        m_dev_potentials = false;
        if (use_resident_loop && m_gpu_scf && !m_solvation
            && (m_method == MethodType::GFN1 || (m_method == MethodType::GFN2 && m_d4_variational_qresp))) {
            const char* env_off = std::getenv("CURCUMA_GRAD_DEVPOT");
            if (!(env_off && env_off[0] == '0'))
                m_dev_potentials = m_gpu_scf->prepareConvergedPotentials(m_method == MethodType::GFN2);
        }
        auto rebuild_host_potentials = [&]() {
            m_pot.reset();
            addCoulombShellPotential(m_pot);
            addThirdOrderPotential(m_pot);
            if (m_method == MethodType::GFN2) addMultipolePotential(m_pot);

        // Exact GFN2 self-consistent-D4 q-response (F5): tblite folds dE_D4/dq into
        // pot%vat (dispersion%get_potential) so the charge-Pulay -P·(v_a+v_b)·dS + the
        // W = Σf·ε·ccᵀ term (ε already carry the in-SCF D4 potential) reproduce the full
        // dE_D4/dq·dq/dR variationally — no separate CPSCF/EEQ response. Mirror it here
        // (m_d4_variational_qresp set above, before calcDispersionEnergy).
        if (m_method == MethodType::GFN2 && m_d4_variational_qresp)
            addDispersionPotential(m_pot);
        // Solvation v_at must be in m_pot here too: calculateGradient's overlap-
        // derivative (Pulay) term uses -P·(v_ao(μ)+v_ao(ν))·dS/dR with the FULL
        // converged potential, and W is built from the solvation-polarized orbital
        // energies — omitting v_solv would make the two inconsistent. The explicit
        // ∂E_solv/∂R is added separately below via m_solvation->addGradient.
        addSolvationPotential(m_pot);
        };
        const auto t_gp1 = clock::now();
        if (!m_dev_potentials) rebuild_host_potentials();
        const auto t_gp2 = clock::now();

        // Stage 4 (Claude Generated): native nuclear gradient on the device when the
        // SCF ran device-resident. GFN1: repulsion/H0-Pulay/Coulomb (1/2/3) on the
        // GPU; GFN2 additionally the multipole-integral Pulay on the GPU and the
        // direct multipole interaction (section 5) on the host. Dispersion (3b) +
        // CN chain-rule (4) on the host. Falls back to the full CPU gradient if the
        // backend has no gradient.
        bool gpu_grad = false;
        const auto t_gp3 = clock::now();
        if (use_gpu_resident && m_gpu_scf && m_gpu_scf->supportsGradient())
            gpu_grad = calculateGradientGpu();
        if (verb >= 3 || std::getenv("CURCUMA_GPU_PROFILE"))
            CurcumaLogger::result(fmt::format(
                "  Gradient host phases: device potentials {}, host potential rebuild {:.1f} ms, "
                "calculateGradientGpu {:.1f} ms",
                m_dev_potentials ? "yes" : "no", ms(t_gp1, t_gp2), ms(t_gp3, clock::now())));
        if (!gpu_grad && m_dev_potentials) {   // host fallback needs the host potentials
            m_dev_potentials = false;
            rebuild_host_potentials();
        }
        if (!gpu_grad) {
            ensureHostWavefunction();
            ensureHostMultipoleIntegrals();
            calculateGradient();   // fills m_gradient in Eh/Bohr
        }

        // Explicit solvation nuclear gradient ∂E_solv/∂R at the converged charges
        // (Born radii / SASA / HB-surface / CM5 derivatives). Added in Eh/Bohr,
        // i.e. in the SAME space as calculateGradient, BEFORE the unit conversion.
        if (m_solvation)
            m_solvation->addGradient(m_atoms, m_geometry * AA_TO_AU, m_wfn.q_at, m_gradient);

        // Debug (Jul 2026, F5): env-gated GFN2 multipole-gradient audit. Runs in
        // Eh/Bohr (before the unit conversion), restores m_gradient on exit.
        if (m_method == MethodType::GFN2 && std::getenv("CURCUMA_MP_GRAD_AUDIT"))
            auditGfn2GradientTerms();

        m_gradient /= au;          // Eh/Bohr → Eh/Å: 1 Eh/Bohr = (1/au) Eh/Å
    }
    const auto t_end = clock::now();

    if (verb >= 2) {
        CurcumaLogger::info("Timing breakdown:");
        CurcumaLogger::info_fmt("  setup        : {:8.2f} ms", ms(t0, t_setup));
        // Claude Generated (Sep 2026): the device uploads/EEQ guess/D4 warm-up between the
        // integral setup and the first SCF iteration were counted nowhere but TOTAL.
        CurcumaLogger::info_fmt("  pre-SCF      : {:8.2f} ms", ms(t_setup, t_scf_start));
        if (verb >= 3 || std::getenv("CURCUMA_GPU_PROFILE")) {
            for (const auto& [name, t] : pre_t)
                CurcumaLogger::info_fmt("    pre-SCF: {:<62s} {:8.2f} ms", name, t);
        }
        CurcumaLogger::info_fmt("  SCF ({:3d} it) : {:8.2f} ms", m_scf_iterations, ms(t_scf_start, t_scf_end));
        CurcumaLogger::info_fmt("  post-SCF E   : {:8.2f} ms", ms(t_scf_end, t_energies));
        if (gradient)
            CurcumaLogger::info_fmt("  gradient     : {:8.2f} ms", ms(t_grad0, t_end));
        CurcumaLogger::info_fmt("  TOTAL        : {:8.2f} ms", ms(t0, t_end));
    }
    if (m_gpu_scf && verb >= 1) {
        const std::string prof = m_gpu_scf->profileReport();
        if (!prof.empty()) CurcumaLogger::result(prof);
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
        if (m_t_xfx_copy > 0.0)
            CurcumaLogger::info_fmt("      - of which F copy (row->col major) : {:8.2f} ms ({:5.2f}/it)",
                                    m_t_xfx_copy, m_t_xfx_copy / it);
        CurcumaLogger::info_fmt("      - eigensolve branch calls: lapack {}, native {}, fp32 {}",
                                m_eig_calls_lapack, m_eig_calls_native, m_eig_calls_fp32);
        if (m_blas_threads > 0)
            CurcumaLogger::info_fmt("      - BLAS/LAPACK threads in the solve : {}", m_blas_threads);
        CurcumaLogger::info_fmt("    - dsyevd      : {:8.2f} ms ({:5.2f}/it)", m_t_diag,   m_t_diag / it);
        CurcumaLogger::info_fmt("    - back-transf : {:8.2f} ms ({:5.2f}/it)", m_t_back,   m_t_back / it);
        CurcumaLogger::info_fmt("    - density P   : {:8.2f} ms ({:5.2f}/it)", m_t_dens,   m_t_dens / it);
        CurcumaLogger::info_fmt("    (solve detail sum = {:.2f} ms, untimed = {:.2f} ms)",
                                solve_sum, acc_solve - solve_sum);
        CurcumaLogger::info_fmt("  populations     : {:8.2f} ms ({:5.2f}/it)", acc_mull,   acc_mull / it);
        CurcumaLogger::info_fmt("  energy/mix      : {:8.2f} ms ({:5.2f}/it)", acc_energy, acc_energy / it);
        CurcumaLogger::info_fmt("    - coulomb {:.2f}/it, third order {:.2f}/it, multipole {:.2f}/it, "
                                "band Tr(P H0) {:.2f}/it, rest (mixing) {:.2f}/it",
                                acc_e_coul / it, acc_e_third / it, acc_e_mp / it, acc_e_band / it,
                                (acc_energy - acc_e_coul - acc_e_third - acc_e_mp - acc_e_band) / it);
    }

    // Reject a converged solution whose charges are physically impossible and redo the
    // whole calculation from the bare-H0 guess. See the note on m_force_h0_guess in the
    // header. The bound is generous: GFN Mulliken charges of a real system stay well
    // inside +-2 e, and a bare monoatomic ion cannot exceed its own formal charge, so
    // 4 e on top of the molecular charge can only be reached by a runaway.
    if (!m_in_scf_retry && m_wfn.q_at.size() == m_atomcount) {
        const double q_max = m_wfn.q_at.cwiseAbs().maxCoeff();
        const double q_bound = 4.0 + std::abs(static_cast<double>(m_charge));
        if (q_max > q_bound) {
            CurcumaLogger::warn_fmt(
                "SCF converged to an implausible charge distribution (max |q| = {:.2f} e > "
                "{:.2f}); repeating from the bare-H0 guess", q_max, q_bound);
            m_in_scf_retry = true;
            m_force_h0_guess = true;
            m_warmstart_q_sh.resize(0);
            m_scf_history.clear();
            const double e_retry = Calculation(gradient);
            m_in_scf_retry = false;
            m_force_h0_guess = false;
            // Claude Generated (Sep 2026): if the bare-H0 retry lands on an impossible charge
            // distribution too, the result is not usable - it used to be returned anyway.
            if (!hasError() && m_wfn.q_at.size() == m_atomcount) {
                const double q_retry = m_wfn.q_at.cwiseAbs().maxCoeff();
                if (!(q_retry <= q_bound))   // also catches NaN
                    setHardError(fmt::format("native SCF converged to an implausible charge distribution "
                                             "(max |q| = {:.2f} e > {:.2f}) also from the bare-H0 guess",
                                             q_retry, q_bound));
            }
            return e_retry;
        }
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
    ensureHostIntegrals(__func__);   // Claude Generated (Oct 2026): host copies of S/H0/L/gamma are lazy
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
 *  auditGfn2GradientTerms() — Claude Generated (Jul 2026, F5 diagnosis)
 *
 *  Isolate each gfn2-specific analytic gradient term (multipole/AES, ES2 Coulomb,
 *  D4 dispersion) via the per-term gates and FD-check it against a frozen-DENSITY
 *  numerical derivative. The FD freezes the density matrix P and holds the atomic
 *  CHARGES q FIXED at q0 — the charge response ∂E/∂q·∂q/∂R belongs to the global
 *  energy-weighted-density term, not the per-term gradient — while recomputing the
 *  multipole MOMENTS from P0 (so the multipole FD carries §5 + AP5b + mrad/CN).
 *  A term whose ratio ≠ 1.000 carries the ~1.3% gfn2 metal gradient residual.
 *  Debug only; leaves m_gradient and the geometry restored on exit.
 * ------------------------------------------------------------------------- */
void XTB::auditGfn2GradientTerms()
{
    ensureHostIntegrals(__func__);   // Claude Generated (Oct 2026): host copies of S/H0/L/gamma are lazy
    if (m_method != MethodType::GFN2 || !m_mp_initialized) {
        CurcumaLogger::warn("auditGfn2GradientTerms: needs a converged GFN2 + multipole state");
        return;
    }
    const int nat = m_atomcount;

    // Consistency probe (Jul 2026, F5): the −Tr(W·dS) gradient term uses W = Cᵒᶜᶜ·diag(f·ε)·Cᵒᶜᶜᵀ,
    // which cancels the density response ONLY if the stored density P equals Cᵒᶜᶜ·diag(f)·Cᵒᶜᶜᵀ.
    // A mismatch (e.g. a Broyden/device density not rebuilt from the final C) breaks that
    // cancellation and would show up exactly as the ~1.3% electronic-response residual.
    if (m_wfn.C.rows() == m_wfn.P.rows() && m_wfn.focc.size() == m_wfn.C.cols()) {
        const int nao_ = m_wfn.P.rows();
        int ncol = 0; for (int i = 0; i < nao_; ++i) if (m_wfn.focc(i) > 1e-12) ncol = i + 1;
        Matrix Prec = m_wfn.C.leftCols(ncol)
                    * m_wfn.focc.head(ncol).asDiagonal()
                    * m_wfn.C.leftCols(ncol).transpose();
        CurcumaLogger::result(fmt::format(
            "  [consistency] ||P - C·f·Cᵀ|| = {:.3e}  (||P||={:.3e})  ncol={}",
            (m_wfn.P - Prec).norm(), m_wfn.P.norm(), ncol));
    }
    // Generalized-eigenproblem residual: do the stored eps/C satisfy F·C = S·C·diag(eps)
    // for the FINAL full Fock m_F? W = C·(f·eps)·Cᵀ is only the correct energy-weighted
    // density (so −Tr(W·dS) cancels the response) if they do. A nonzero residual means the
    // eps used for W are inconsistent with the converged Fock — the −Tr(W·dS) bug.
    if (m_F.rows() == m_wfn.C.rows() && m_S.rows() == m_wfn.C.rows()
        && m_wfn.eps.size() == m_wfn.C.cols()) {
        Matrix lhs = m_F * m_wfn.C;
        Matrix rhs = m_S * m_wfn.C * m_wfn.eps.asDiagonal();
        CurcumaLogger::result(fmt::format(
            "  [consistency] ||F·C - S·C·ε|| = {:.3e}  (||F·C||={:.3e})",
            (lhs - rhs).norm(), lhs.norm()));
    }

    const Matrix geom0   = m_geometry;   // Angstrom
    const Matrix P0      = m_wfn.P;
    const Vector q0      = m_wfn.q_at;    // converged charges — FIXED in the FD
    const Vector qsh0    = m_wfn.q_sh;
    const Eigen::MatrixXd dp0 = m_wfn.dp_at;
    const Eigen::MatrixXd qp0 = m_wfn.qp_at;
    const Matrix g_saved = m_gradient;   // caller's gradient — restored on exit

    // --- analytic per-term gradients (Eh/Bohr) via the per-term gates ---
    auto grad_with = [&](bool mp, bool coul, bool d4, bool rep, bool h0only,
                         bool cnoff = false) -> Matrix {
        m_mp_grad_off = mp; m_coulomb_grad_off = coul; m_d4_grad_off = d4;
        m_repulsion_grad_off = rep; m_sval_h0_only = h0only; m_cnchain_off = cnoff;
        calculateGradient();
        m_mp_grad_off = m_coulomb_grad_off = m_d4_grad_off = false;
        m_repulsion_grad_off = m_sval_h0_only = m_cnchain_off = false;
        return m_gradient;
    };
    const Matrix g_full = grad_with(false, false, false, false, false);
    const Matrix g_mp   = g_full - grad_with(true , false, false, false, false);
    const Matrix g_coul = g_full - grad_with(false, true , false, false, false);
    const Matrix g_d4   = g_full - grad_with(false, false, true , false, false);
    // Band term Tr(P·dH0/dR): all potential terms off + sval reduced to 2·P·h_av
    // (drops −2W and the charge-Pulay), leaving 2·P·h_av·dS + G_shpoly + H0/CN chain.
    const Matrix g_band = grad_with(true, true, true, true, true);
    // Band decomposition (Jul 2026, F5): the CN chain (section 4) OFF isolates the
    // direct band part 2·P·h_av·dS + G_shpoly (frozen-CN); the remainder is the
    // H0/CN chain that flows through dEdcn → section 4.
    const Matrix g_band_nocn = grad_with(true, true, true, true, true, /*cnoff=*/true);
    const Matrix g_band_cn   = g_band - g_band_nocn;

    // --- frozen-P, fixed-q component energies at a geometry ---
    auto e_mp = [&](const Matrix& geom) -> double {          // §5 + AP5b + mrad/CN
        m_geometry = geom;
        Vector cn = computeCoordinationNumbers(); Vector se; getSelfEnergies(cn, se);
        Matrix S, H0; getHamiltonianH0(se, S, H0); m_S = S; m_H0 = H0; m_coordination_numbers = cn;
        buildGammaMatrix(); setupMultipole();
        m_wfn.P = P0; updatePopulations(m_S);                // recompute moments from P0
        m_wfn.q_at = q0; m_wfn.q_sh = qsh0;                  // freeze charges
        return energyMultipole();
    };
    // evaluateComponentsAtFixedDensity injects (P0, q0, moments0) and rebuilds at the
    // current geometry; E_coulomb_shell = explicit ES2 (fixed q), E_dispersion = D4
    // direct+CN (audit mode, no q-response). Moments do not enter either energy.
    auto e_coul = [&](const Matrix& geom) -> double {
        m_geometry = geom; evaluateComponentsAtFixedDensity(P0, q0, qsh0, dp0, qp0);
        return m_E_coulomb_shell;
    };
    auto e_d4 = [&](const Matrix& geom) -> double {
        m_geometry = geom; evaluateComponentsAtFixedDensity(P0, q0, qsh0, dp0, qp0);
        return m_E_dispersion;
    };
    // Band energy Tr(P0·H0(R)) at the frozen density: evaluateComponentsAtFixedDensity
    // rebuilds m_H0 at R and injects m_wfn.P = P0, so this is exactly the band term.
    auto e_band = [&](const Matrix& geom) -> double {
        m_geometry = geom; evaluateComponentsAtFixedDensity(P0, q0, qsh0, dp0, qp0);
        return (m_wfn.P.cwiseProduct(m_H0)).sum();
    };
    // Band energy with the self-energies FROZEN at the R0 coordination number
    // (se0 = getSelfEnergies(cn0)). H0(R) still carries the geometry-dependence of
    // S and the shpoly factor pi_ij, but NOT the CN-dependence of avg_eps. Its FD
    // is therefore the direct band part (2·P·h_av·dS + G_shpoly), matching the
    // analytic g_band_nocn. FD(e_band) − FD(e_band_fixedcn) = the H0/CN chain part.
    m_geometry = geom0;                              // (already geom0 here; explicit)
    const Vector cn0 = computeCoordinationNumbers();
    Vector se0; getSelfEnergies(cn0, se0);
    auto e_band_fixedcn = [&](const Matrix& geom) -> double {
        m_geometry = geom;
        Matrix S, H0; getHamiltonianH0(se0, S, H0);
        return (P0.cwiseProduct(H0)).sum();
    };

    const double h = 1.0e-4;   // Angstrom
    // CRITICAL (Jul 2026, F5): computeCoordinationNumbers() memoises on m_cn_cache,
    // invalidated only by UpdateMolecule/InitialiseMolecule. The FD sets m_geometry
    // DIRECTLY, so the cache must be force-invalidated before each evaluation or the
    // FD silently freezes CN at geom0 (dropping the H0/CN chain from every component).
    auto fd = [&](auto efn) -> Matrix {
        Matrix g = Matrix::Zero(nat, 3);
        for (int i = 0; i < nat; ++i)
            for (int j = 0; j < 3; ++j) {
                Matrix gp = geom0, gm = geom0; gp(i, j) += h; gm(i, j) -= h;
                m_cn_cache_valid = false; const double ep = efn(gp);
                m_cn_cache_valid = false; const double em = efn(gm);
                g(i, j) = (ep - em) / (2.0 * h) / AA_TO_AU;   // Eh/Ang -> Eh/Bohr
            }
        return g;
    };
    const Matrix fd_mp   = fd(e_mp);
    const Matrix fd_coul = fd(e_coul);
    const Matrix fd_d4   = fd(e_d4);
    const Matrix fd_band = fd(e_band);
    const Matrix fd_band_nocn = fd(e_band_fixedcn);     // direct (overlap+shpoly) only
    const Matrix fd_band_cn   = fd_band - fd_band_nocn; // H0/CN chain only

    // NOTE: a frozen-P TOTAL-energy FD is deliberately NOT reported. ∂E_total/∂R|_P omits
    // the −Tr(W·dS) orthonormality term (from the C†SC=I constraint, not from ∂E/∂R|_P), so
    // it does not equal the true gradient and cannot audit the band/charge-response terms.
    // The per-component energies above are free of that constraint, so their FDs are valid.

    // --- restore consistent geom0 state + caller's gradient ---
    m_geometry = geom0;
    m_cn_cache_valid = false;   // FD left a perturbed-geom CN cached; force recompute
    { Vector cn = computeCoordinationNumbers(); Vector se; getSelfEnergies(cn, se);
      Matrix S, H0; getHamiltonianH0(se, S, H0); m_S = S; m_H0 = H0; m_coordination_numbers = cn;
      buildGammaMatrix(); setupMultipole(); m_wfn.P = P0; updatePopulations(m_S);
      m_wfn.q_at = q0; m_wfn.q_sh = qsh0; m_wfn.dp_at = dp0; m_wfn.qp_at = qp0; }
    m_gradient = g_saved;

    auto report = [&](const char* name, const Matrix& a, const Matrix& f) {
        double num = 0.0, maxa = 0.0;
        for (int i = 0; i < nat; ++i) for (int j = 0; j < 3; ++j) {
            const double d = a(i,j) - f(i,j); num += d*d; maxa = std::max(maxa, std::fabs(d)); }
        CurcumaLogger::result(fmt::format(
            "  {:<11} RMS(ana-FD)={:.3e}  max={:.3e}  |ana|={:.3e}  |FD|={:.3e}  ratio={:.5f}",
            name, std::sqrt(num / (3.0 * nat)), maxa, a.norm(), f.norm(),
            std::sqrt(a.squaredNorm() / std::max(1e-30, f.squaredNorm()))));
    };
    CurcumaLogger::result("GFN2 per-term gradient audit [Eh/Bohr]  analytic vs frozen-density(fixed-q) FD   ratio 1.000 = exact");
    report("multipole",   g_mp,   fd_mp);
    report("coulomb-ES2", g_coul, fd_coul);
    report("d4(dir+CN)",  g_d4,   fd_d4);
    report("band Tr(PdH0)", g_band, fd_band);
    report("band:ovlp+shp", g_band_nocn, fd_band_nocn);
    report("band:H0-CN",    g_band_cn,   fd_band_cn);
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
                nprimFor(z, ish),
                m_sto6g_legacy_4sp);
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

    // X-I1: detect d shells. The device integral/SCF/gradient kernels are s/p-only,
    // so d systems use the CPU integral path (gated below). Claude Generated.
    m_has_dshell = false;
    for (int s = 0; s < m_basis.nsh; ++s)
        if (m_basis.ang_sh[s] > 1) { m_has_dshell = true; break; }
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
            // GFN2 stores p_kcn and p_shpoly indexed by ANGULAR MOMENTUM (tblite
            // gfn2.f90: p_kcn(0:2), p_shpoly(0:2)); its shells never repeat an
            // angular momentum, so mapping shell->l reorders the [d,s,p] transition-
            // metal shells correctly (main-group is unaffected: ang_sh==ish).
            // GFN1 must NOT do this: its basis has valence+polarisation shells that
            // share an angular momentum (e.g. H = two s shells with distinct p_kcn),
            // and tblite gfn1.f90 keeps p_kcn shell-indexed (max_shell); indexing by
            // l would collapse those shells. GFN1 keeps its historic [ish] indexing.
            const int l = m_basis.ang_sh[sh_idx];
            if (m_method == MethodType::GFN1) {
                // GFN1: p_selfenergy and p_kcn are shell-indexed (tblite gfn1.f90
                // max_shell); p_shpoly is angular-momentum-indexed (0:2), applied to
                // every shell incl. polarisation (tblite get_shpoly). ang-index shpoly.
                m_h0.selfenergy[sh_idx] = gfn1_params::p_selfenergy[z - 1][ish];
                m_h0.kcn       [sh_idx] = gfn1_params::p_kcn       [z - 1][ish];
                m_h0.shpoly    [sh_idx] = gfn1_params::p_shpoly    [z - 1][l];
            } else {
                m_h0.selfenergy[sh_idx] = gfn2_params::p_selfenergy[z - 1][ish];
                m_h0.kcn       [sh_idx] = gfn2_params::p_kcn       [z - 1][l];
                m_h0.shpoly    [sh_idx] = gfn2_params::p_shpoly    [z - 1][l];
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
        // reference_occ is indexed by ANGULAR MOMENTUM (tblite gfn{1,2}.f90 (0:2,elem),
        // refocc = reference_occ(shell%ang)); indexing by shell position scrambles the
        // [d,s,p] transition-metal shells (main-group [s,p(,d)] is unaffected: ang==pos).
        // GFN1 additionally has valence+polarisation shells that share an angular
        // momentum (only H: two s-shells): the occupation goes to the VALENCE shell
        // (first of each l), polarisation shells get 0 (tblite gfn1 get_reference_occ).
        bool seen_l[3] = {false, false, false};
        for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
            const int sh_idx = m_basis.ish_at[iat] + ish;
            const int l = m_basis.ang_sh[sh_idx];
            double refocc;
            if (m_method == MethodType::GFN1) {
                const bool is_valence = (l < 0 || l > 2) ? true : !seen_l[l];
                refocc = is_valence ? gfn1_params::reference_occ[z - 1][l] : 0.0;
            } else {
                refocc = gfn2_params::reference_occ[z - 1][l];
            }
            if (l >= 0 && l <= 2) seen_l[l] = true;
            m_wfn.n0_sh(sh_idx) += refocc;
            m_wfn.n0_at(iat)    += refocc;
            total               += refocc;
        }
    }
    m_wfn.nocc = total;
    // Split into alpha/beta, verbatim from tblite get_alpha_beta_occupation
    // (src/tblite/wavefunction/type.f90:162-176). m_spin carries multiplicity - 1,
    // i.e. the number of unpaired electrons (tblite nuhf). Closed shell (nuhf = 0)
    // gives nalpha = nbeta = nocc/2, and the open-shell code paths then never run.
    // Claude Generated (Sep 2026).
    {
        const double nuhf = nUhf();
        const double diff = std::min(nuhf, m_wfn.nocc);
        const double ntmp = m_wfn.nocc - diff;
        m_nalpha = ntmp / 2.0 + diff;
        m_nbeta = ntmp / 2.0;
        if (nuhf > 0.0 && CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::param("open_shell_nuhf", static_cast<int>(std::lround(nuhf)));
            CurcumaLogger::param("n_alpha", m_nalpha);
            CurcumaLogger::param("n_beta", m_nbeta);
        }
        // Claude Generated (Oct 2026): an electron count whose parity contradicts -spin
        // (odd electrons without -spin, or -spin 2 on a doublet) gives fractional alpha/beta
        // occupations - tblite's behaviour, but xtb 6.7.1 defaults an odd count to a doublet.
        // CH3 without -spin: -3.56397536 Eh vs xtb -3.56265832. Warn, do not change numbers.
        const long nel = std::lround(m_wfn.nocc);
        if (std::abs(m_wfn.nocc - nel) < 1e-6 && ((nel - std::lround(nuhf)) % 2) != 0
            && CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::warn(fmt::format("native xTB: {} electrons with {} unpaired (-spin) "
                "have inconsistent parity; occupations are fractional. For a radical pass "
                "-spin <multiplicity-1>.", nel, std::lround(nuhf)));
        }
    }
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
    m_host_ints_deferred = false;   // Oct 2026: the lazy host copies belong to the old geometry
    m_mp_amat_deferred = false;

    // 3. Fock matrix depends on geometry via H0 + potential
    m_F.resize(0, 0);

    // 4. Gamma matrix depends on interatomic distances (Klopman-Ohno R_AB)
    m_gamma.resize(0, 0);

    // X-I3: coordination numbers depend on the geometry — drop the memoised CN.
    m_cn_cache_valid = false;

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
    j["entropy"]        = m_E_entropy;   // -T*S, already included in "electronic"
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
void XTB::ensureHostIntegrals(const char* who) const
{
    if (m_host_ints_deferred) const_cast<XTB*>(this)->ensureHostIntegralsImpl(who);
}

void XTB::ensureHostIntegralsImpl(const char* who)
{
    m_host_ints_deferred = false;    // first: the rebuild below calls functions that would ensure again
    if (CurcumaLogger::get_verbosity() >= 3 || std::getenv("CURCUMA_GPU_PROFILE"))
        CurcumaLogger::result(fmt::format("  host copies of S/H0/L/gamma fetched (first consumer: {})", who ? who : "getter"));
    const int nao = m_basis.nao, nsh = m_basis.nsh;
    bool got = false;
    if (m_gpu_scf) {
        m_S.resize(nao, nao);
        m_H0.resize(nao, nao);
        if (m_gpu_scf->downloadOverlapInto(m_S.data()) && m_gpu_scf->downloadH0Into(m_H0.data())) {
            transposeSquareInPlace(m_S.data(), nao);
            transposeSquareInPlace(m_H0.data(), nao);
            got = m_gpu_scf->downloadCholesky(m_X) && m_gpu_scf->downloadGamma(m_gamma)
                  && m_X.rows() == nao && m_gamma.rows() == nsh;
        }
    }
    if (!got) {
        // The device no longer holds this geometry (or a backend without raw downloads): the host build.
        const Vector cn = computeCoordinationNumbers();
        Vector se;
        getSelfEnergies(cn, se);
        Matrix S, H0;
        getHamiltonianH0(se, S, H0);
        m_S = S;
        m_H0 = H0;
        buildOrthonormalizer();
        buildGammaMatrix();
    }
}

// Claude Generated (Oct 2026): the dftd4 covalent CN of the next D4 GenerateParameters() from the device
// (CURCUMA_D4_CN_DEVICE=0 keeps the host loop).
void XTB::primeD4CnFromDevice(const Matrix& geom_bohr) const
{
    if (!m_gpu_scf || !m_d4_generator) return;
    const char* env_off = std::getenv("CURCUMA_D4_CN_DEVICE");
    if (env_off && env_off[0] == '0') return;
    const int nat = m_atomcount;
    std::vector<double> xyz(3 * static_cast<size_t>(nat)), cn(nat, 0.0);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) xyz[3 * i + k] = geom_bohr(i, k);
    const auto tab = curcuma::dispersion::d4CovalentCNTables();
    if (m_gpu_scf->d4CovalentCnDevice(nat, m_atoms.data(), xyz.data(), tab.rad, tab.en, tab.kn, tab.k4, tab.k5,
                                      tab.k6, /*cutoff=*/30.0, cn.data()))
        m_d4_generator->setPrecomputedCovalentCN(std::move(cn));
}

double XTB::calcDispersionEnergy(bool need_gradient) const
{
    m_disp_gradient_valid = false;

    if (m_method != MethodType::GFN2) {
        // GFN1: native D3(BJ) dispersion (s6=1.0, s8=2.4, a1=0.63, a2=5.0).
        // Claude Generated (May 2026): analytical gradient replaces central
        // finite differences. The direct geometry gradient + CN chain rule
        // are computed in one pass through the D3 pair list.
        const auto t_d3init0 = std::chrono::steady_clock::now();
        const bool d3_first = !m_d3_generator;
        if (!m_d3_generator) {
            m_d3_generator = std::make_unique<::D3ParameterGenerator>(
                ::D3ParameterGenerator::createForGFN1());
        }
        if (d3_first && (CurcumaLogger::get_verbosity() >= 3 || std::getenv("CURCUMA_GPU_PROFILE")))
            CurcumaLogger::result(fmt::format("  GFN1 D3 (host): generator construction {:.1f} ms",
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_d3init0).count()));
        ::D3ParameterGenerator& d3 = *m_d3_generator;
        // Lightweight prepare (CN + Gaussian weights only) — getEnergyAndGradient
        // computes C6/C8/energy/gradient directly from the references, so the
        // O(N²) JSON pair list that GenerateParameters builds is pure overhead
        // here (≈0.5 s on complex/231). Claude Generated.
        // Claude Generated (Oct 2026): CN, weights, pair energy / gradient / dE/dCN and the CN chain rule on the
        // device (XtbGpuContext::d3DispersionDevice); the host code below is the fallback and the reference.
        // CURCUMA_D3_DEVICE=0 forces it.
        if (m_gpu_scf) {
            const char* env_off = std::getenv("CURCUMA_D3_DEVICE");
            if (!(env_off && env_off[0] == '0')) {
                const int nat = m_atomcount;
                const auto tb = d3.buildDeviceTables(m_atoms);
                const auto& rcov = CNCalculator::covalentRadiiD3();
                std::vector<double> xyz(3 * static_cast<size_t>(nat)), gflat(need_gradient ? 3 * nat : 0, 0.0);
                for (int i = 0; i < nat; ++i)
                    for (int k = 0; k < 3; ++k) xyz[3 * i + k] = m_geometry(i, k);
                double e_dev = 0.0;
                if (tb.nt > 0 && m_gpu_scf->d3DispersionDevice(
                        nat, m_atoms.data(), xyz.data(), tb.nt, tb.atom_type.data(), tb.c6tab.data(), tb.r6tab.data(),
                        tb.nref.data(), tb.refcn.data(), rcov.data(), static_cast<int>(rcov.size()),
                        tb.s6, tb.s8, tb.a1, tb.a2, /*k1=*/16.0, /*k2=*/4.0 / 3.0,
                        /*cn_cut_ang=*/25.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM,
                        CurcumaUnit::Length::angstrom_to_bohr_or_legacy(1.88972612546), au,
                        need_gradient, &e_dev, gflat.data())) {
                    m_disp_gradient = Matrix::Zero(nat, 3);
                    m_disp_dEdcn = Vector();
                    if (need_gradient) {
                        for (int i = 0; i < nat; ++i)
                            for (int k = 0; k < 3; ++k) m_disp_gradient(i, k) = gflat[3 * i + k];
                    }
                    m_disp_gradient_valid = need_gradient;
                    return e_dev;
                }
            }
        }
        const auto t_d3a = std::chrono::steady_clock::now();
        d3.prepareForEnergyGradient(m_atoms, m_geometry);
        const auto t_d3b = std::chrono::steady_clock::now();

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
        const auto t_d3c = std::chrono::steady_clock::now();

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
        if (CurcumaLogger::get_verbosity() >= 3 || std::getenv("CURCUMA_GPU_PROFILE")) {
            const auto t_d3d = std::chrono::steady_clock::now();
            auto msd = [](std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count(); };
            CurcumaLogger::result(fmt::format("  GFN1 D3 (host): CN + weights {:.1f} ms, energy/gradient/dC6dCN {:.1f} ms, CN chain rule {:.1f} ms",
                                              msd(t_d3a, t_d3b), msd(t_d3b, t_d3c), msd(t_d3c, t_d3d)));
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
        primeD4CnFromDevice(geom_bohr);
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
    double E = 0.0;
    bool disp_grad_device = false;
    // Claude Generated (Oct 2026): everything below on the device in one call (weights, 2-body, ATM with its C6
    // matrices and neighbour list, dE/dCN fold); the host flow after it stays as the fallback and for the
    // backends without the call. CURCUMA_D4_POSTSCF_DEVICE=0 forces the host flow.
    bool d4_all_device = false;
    if (m_gpu_scf && m_gpu_scf->supportsDeviceDispersion() && m_d4_prepared) {
        const char* env_off = std::getenv("CURCUMA_D4_POSTSCF_DEVICE");
        if (!(env_off && env_off[0] == '0')) {
            const int nat = m_atomcount;
            const auto tab = curcuma::dispersion::d4CovalentCNTables();
            std::vector<double> qv(nat), gflat(3 * nat, 0.0), dq(nat, 0.0);
            for (int a = 0; a < nat; ++a) qv[a] = (a < m_wfn.q_at.size()) ? m_wfn.q_at(a) : 0.0;
            double e_dev = 0.0;
            if (m_gpu_scf->d4PostScfDevice(nat, qv.data(), /*with_atm=*/true, 5.0, 0.52, 5.0, 16.0,
                                           m_d4_atm_cutoff, /*cn_cutoff=*/30.0, tab.rad, tab.en,
                                           tab.kn, tab.k4, tab.k5, tab.k6, &e_dev, gflat.data(), dq.data())) {
                m_disp_gradient = Matrix::Zero(nat, 3);
                m_disp_dEdq = Vector::Zero(nat);
                for (int a = 0; a < nat; ++a) {
                    m_disp_gradient(a, 0) = gflat[3 * a + 0];
                    m_disp_gradient(a, 1) = gflat[3 * a + 1];
                    m_disp_gradient(a, 2) = gflat[3 * a + 2];
                    m_disp_dEdq(a) = dq[a];
                }
                m_disp_dEdcn = Vector();   // folded on the device (as the host fold clears it)
                E = e_dev;
                d4_all_device = true;
            }
        }
    }
    // Stage 5 (Part B2): the whole 2-body D4 (energy + gradient + dE/dCN + dE/dq) on the
    // device in one gather, reusing the geometry-fixed reference data uploaded during the SCF
    // (beginDispersion). The host then adds ATM + the CN-distribution + the q-response on top,
    // exactly as for the host evaluator. buildRefWFlat uses the same converged charges
    // (m_wfn.q_at == getZetaCharges, set above) and CN (m_cn_values) the host path uses, so the
    // result matches the host D4Evaluator to FP-reassociation (~1e-12 energy, gradient FD-grade).
    // Returns false → host fallback (CUDA, or if the device reference was not uploaded).
    if (!d4_all_device && m_gpu_scf && m_gpu_scf->supportsDeviceDispersion()) {
        const int nat = m_atomcount;
        std::vector<double> W, dWq, dWc;
        m_d4_generator->buildRefWFlat(m_atoms, m_wfn.q_at, W, dWq, dWc);
        std::vector<double> e_atom(nat, 0.0), gflat(3 * nat, 0.0), dcn(nat, 0.0), dq(nat, 0.0);
        if (m_gpu_scf->dispersionGradient(nat, W.data(), dWq.data(), dWc.data(),
                                          e_atom.data(), gflat.data(), dcn.data(), dq.data())) {
            m_disp_gradient = Matrix::Zero(nat, 3);
            m_disp_dEdcn = Vector::Zero(nat);
            m_disp_dEdq = Vector::Zero(nat);
            double Esum = 0.0;
            for (int a = 0; a < nat; ++a) {
                Esum += e_atom[a];                       // total 2-body E = ½·Σ (gather double-counts)
                m_disp_gradient(a, 0) = gflat[3 * a + 0];
                m_disp_gradient(a, 1) = gflat[3 * a + 1];
                m_disp_gradient(a, 2) = gflat[3 * a + 2];
                m_disp_dEdcn(a) = dcn[a];
                m_disp_dEdq(a)  = dq[a];
            }
            E = 0.5 * Esum;
            disp_grad_device = true;
        }
    }
    if (!disp_grad_device && !d4_all_device) {
        // WP2-ext: let the D4 pair loop fan out over the same gated intra-thread budget.
        m_d4_evaluator->setThreads(effectiveIntraThreads(m_atomcount));
        E = m_d4_evaluator->computeEnergyAndGradient(
            m_atoms, geom_bohr,
            /*with_gradient=*/true,
            m_disp_gradient, m_disp_dEdcn, m_disp_dEdq, want_dEdq);
    }
    const auto td1 = d4clk::now();

    // ATM three-body term (GFN2 s9=5.0). Charge-independent (dftd4 evaluates it at
    // the q=0 reference C6 — matches tblite get_dispersion_nonsc), so no q-response.
    // ACCUMULATES into m_disp_gradient and m_disp_dEdcn (the latter folded together
    // with the 2-body dEdcn via the D4 covalent CN below). Closes the triose C-path
    // remainder (see docs/GFN2_D4_STATUS.md). Claude Generated 2026-05-29.
    // ATM 3-body: device when the 2-body went to the device (so beginDispersion's geometry +
    // √r4r2 are resident) and the backend implements it; else the host evaluator. The device
    // gather sums each unique triple once per member → total ATM E = ⅓·Σ e_atom; grad/dEdcn
    // accumulate on top of the 2-body. Params mirror the D4Params block above (GFN2: s9=5.0,
    // a1=0.52, a2=5.0, alpha=16.0; computeATM cutoff default 25.0 Bohr). Returns false → host.
    bool atm_device = d4_all_device;
    if (!d4_all_device && disp_grad_device && m_gpu_scf) {
        const int nat = m_atomcount;
        std::vector<double> c6f, dc6f;
        const auto ta0 = d4clk::now();
        m_d4_generator->buildAtmC6Flat(m_atoms, c6f, dc6f);
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info_fmt("  D4 ATM C6 build  : {:8.2f} ms (host)", d4ms(ta0, d4clk::now()));
        std::vector<double> e_at(nat, 0.0), g_at(3 * nat, 0.0), dcn_at(nat, 0.0);
        if (m_gpu_scf->dispersionATM(nat, c6f.data(), dc6f.data(),
                                     /*s9=*/5.0, /*a1=*/0.52, /*a2=*/5.0, /*alp=*/16.0,
                                     /*cutoff=*/m_d4_atm_cutoff, e_at.data(), g_at.data(), dcn_at.data())) {
            double esum = 0.0;
            for (int a = 0; a < nat; ++a) {
                esum += e_at[a];
                m_disp_gradient(a, 0) += g_at[3 * a + 0];
                m_disp_gradient(a, 1) += g_at[3 * a + 1];
                m_disp_gradient(a, 2) += g_at[3 * a + 2];
                m_disp_dEdcn(a) += dcn_at[a];
            }
            E += esum / 3.0;   // each triple counted by its 3 members
            atm_device = true;
            if (CurcumaLogger::get_verbosity() >= 3)
                CurcumaLogger::info("  D4 ATM 3-body on the device");
        }
    }
    if (!atm_device) {
        E += m_d4_evaluator->computeATM(
            m_atoms, geom_bohr, /*with_gradient=*/true, m_disp_gradient, m_disp_dEdcn,
            m_d4_atm_cutoff);
    }
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
    // When m_d4_variational_qresp is set (exact GFN2, F5), the dE_D4/dq·dq/dR term is
    // handled by the gradient charge-Pulay + W (dE_D4/dq folded into m_pot.v_at), so the
    // separate eeq/CPSCF response here is skipped to avoid double-counting.
    if (m_disp_dEdq.size() == m_atomcount && !m_disp_audit_mode && !m_d4_variational_qresp) {
        // q-response source (m_d4_charge_source) — reached only when NOT variational:
        //   "eeq"   — analytic ∂q/∂x from the single-shot dftd4 EEQ model (one LU solve
        //             + adjoint contraction, ~ms). dftd4-conform, approximate (uses EEQ
        //             not Mulliken charges): ~1% gradient residual on TM complexes.
        //   "cpscf" — explicit ∂q_Mulliken/∂x via the GFN2 CPSCF/Z-vector solve
        //             (computeMullikenChargeResponse). Exact but expensive (~88% of D4
        //             on complex). Superseded by the default variational path, which
        //             gets the same Mulliken response for free through the charge-Pulay.
        // The energy weighting stays on the SCF Mulliken charges in all cases, so this
        // choice only affects the gradient, never the energy.
        if (m_d4_charge_source == "cpscf") {
            // (never deferred: the deferral is disabled for d4_charge_source=cpscf)
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
        primeD4CnFromDevice(geom_bohr);
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

} // namespace curcuma::xtb
