/*
 * <Fragment-wise SCF preconditioner for the native xTB Broyden mixer - implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (Oct 2026). Mathematics and storage form: xtb_precond.h. This file holds
 *   - ScfPrecondGuard (residual-based fallback to the scalar Broyden seed),
 *   - SchwarzPreconditioner (block storage, M*v blockwise),
 *   - numericalScfJacobian / detectCovalentFragments (free helpers),
 *   - the XTB members selectScfPreconditioner, buildFragmentPreconditioner,
 *     writeScfJacobianDump (CURCUMA_SCF_JAC) and loadScfPrecondMatrix
 *     (CURCUMA_SCF_PRECOND_FILE), moved here from xtb_native.cpp without a change of behaviour.
 *
 * AI-generated, machine-tested only; human production testing pending (docs/SCF_PRECOND.md).
 */

#include "xtb_native.h"
#include "xtb_precond.h"

#include "parameters/xtb_params_extra.hpp"
#include "xtb_coulomb.hpp"   // shell hardness / third-order Gamma for the SCF Jacobian dump

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/units.h"
#include "src/core/intra_parallel_context.h"
#include "src/core/energy_calculators/dispersion/d4_charge_model.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <mutex>
#include <thread>

namespace curcuma::xtb {

/* ------------------------------------------------------------------------- *
 *  ScfPrecondGuard
 * ------------------------------------------------------------------------- */
ScfPrecondGuard::ScfPrecondGuard() = default;
ScfPrecondGuard::~ScfPrecondGuard() = default;
ScfPrecondGuard::ScfPrecondGuard(const ScfPrecondGuard&) = default;
ScfPrecondGuard::ScfPrecondGuard(ScfPrecondGuard&&) noexcept = default;
ScfPrecondGuard& ScfPrecondGuard::operator=(const ScfPrecondGuard&) = default;
ScfPrecondGuard& ScfPrecondGuard::operator=(ScfPrecondGuard&&) noexcept = default;

ScfPrecondGuard ScfPrecondGuard::parse(const std::string& spec, const std::string& origin, bool warn)
{
    ScfPrecondGuard g;
    if (spec.empty() || spec == "off") return g;
    g.name = spec;
    std::vector<std::string> tok;
    std::string cur;
    for (const char* p = spec.c_str(); ; ++p) {
        if (*p == ':' || *p == '\0') { tok.push_back(cur); cur.clear(); if (*p == '\0') break; }
        else cur += *p;
    }
    try {
        std::vector<std::string> num;
        for (size_t i = 1; i < tok.size(); ++i) {
            if (tok[i] == "keep") g.keep = true;
            else if (tok[i] == "best") g.best = true;
            else num.push_back(tok[i]);
        }
        if (tok[0] == "g1" && num.size() == 1) {
            g.rule = Rule::Growth; g.tau = std::stod(num[0]);
        } else if (tok[0] == "g2" && (num.size() == 1 || num.size() == 2)) {
            g.rule = Rule::Stagnation; g.n0 = std::stoi(num[0]);
            if (num.size() == 2) g.ratio = std::stod(num[1]);
            else if (g.n0 == 4) g.ratio = 0.15;
            else if (g.n0 == 6) g.ratio = 0.065;
            else throw std::invalid_argument("g2 needs an explicit ratio for this n0");
        } else {
            throw std::invalid_argument("unknown rule");
        }
    } catch (const std::exception&) {
        if (warn)
            CurcumaLogger::warn(origin + "='" + spec
                                + "' not understood (g1:<tau>[:keep][:best] | g2:<n0>[:<ratio>][:keep][:best]); guard off");
        g.rule = Rule::Off;
    }
    return g;
}

ScfPrecondGuard ScfPrecondGuard::fromEnv()
{
    const char* e = std::getenv("CURCUMA_SCF_PRECOND_GUARD");
    if (!e || e[0] == '\0') return ScfPrecondGuard{};
    return parse(e, "CURCUMA_SCF_PRECOND_GUARD");
}

bool ScfPrecondGuard::fires(int n, const Vector& x, const Vector& F)
{
    bool new_best = false;
    const bool hit = firesNorm(n, F.norm(), &new_best);
    if (new_best) { x_best = x; F_best = F; }
    return hit;
}

bool ScfPrecondGuard::firesNorm(int n, double f, bool* new_best)
{
    // Order as before the split (Oct 2026): the best pair is decided against fmin of the earlier
    // iterations, then fmin is updated.
    if (new_best) *new_best = best && (n == 0 || f < fmin);
    if (n == 0) { f0 = f; fmin = f; return false; }
    bool hit = false;
    if (rule == Rule::Growth)
        hit = (n >= 2 && f > tau * fmin);
    else if (rule == Rule::Stagnation)
        hit = (n == n0 && f >= ratio * f0);
    fmin = std::min(fmin, f);
    return hit;
}

/* ------------------------------------------------------------------------- *
 *  SchwarzPreconditioner: M = sum_i R_i^T D_i B_i R_i + C
 * ------------------------------------------------------------------------- */
SchwarzPreconditioner::Block::Block() = default;
SchwarzPreconditioner::Block::~Block() = default;
SchwarzPreconditioner::Block::Block(Block&&) noexcept = default;
SchwarzPreconditioner::Block& SchwarzPreconditioner::Block::operator=(Block&&) noexcept = default;

SchwarzPreconditioner::SchwarzPreconditioner(int dim, double alpha)
    : m_dim(dim)
    , m_alpha(alpha)
    , m_fallback(Vector::Constant(dim, alpha))
{
}

bool SchwarzPreconditioner::addBlock(Block block)
{
    const int n = static_cast<int>(block.index.size());
    if (n == 0 || block.B.rows() != n || block.B.cols() != n) return false;
    if (block.weight.size() != 0 && block.weight.size() != n) return false;
    for (int i : block.index)
        if (i < 0 || i >= m_dim) return false;
    for (int i : block.index) m_fallback(i) = 0.0;   // covered: C has no entry here
    m_blocks.push_back(std::move(block));
    return true;
}

Vector SchwarzPreconditioner::apply(const Vector& v) const
{
    // C v: alpha on indices outside every block.
    Vector out = m_fallback.cwiseProduct(v);
    Vector vl, y;
    for (const Block& b : m_blocks) {
        const int n = static_cast<int>(b.index.size());
        vl.resize(n);
        for (int k = 0; k < n; ++k) vl(k) = v(b.index[k]);   // R_i v
        y.noalias() = b.B * vl;                                // B_i R_i v
        if (b.weight.size() == 0)
            for (int k = 0; k < n; ++k) out(b.index[k]) += y(k);              // R_i^T B_i R_i v
        else
            for (int k = 0; k < n; ++k) out(b.index[k]) += b.weight(k) * y(k); // R_i^T D_i B_i R_i v
    }
    return out;
}

std::size_t SchwarzPreconditioner::storedDoubles() const
{
    std::size_t n = static_cast<std::size_t>(m_fallback.size());
    for (const Block& b : m_blocks)
        n += static_cast<std::size_t>(b.B.size()) + static_cast<std::size_t>(b.weight.size())
             + b.index.size() / 2;   // int indices, counted as half a double
    return n;
}

GpuPrecondFlat::GpuPrecondFlat() = default;
GpuPrecondFlat::~GpuPrecondFlat() = default;

bool SchwarzPreconditioner::flatten(std::vector<int>& idx_off, std::vector<int>& idx,
                                    std::vector<int>& b_off, std::vector<double>& B,
                                    std::vector<double>& fallback, int& max_dim) const
{
    idx_off.assign(1, 0);
    b_off.assign(1, 0);
    idx.clear();
    B.clear();
    max_dim = 0;
    std::size_t nB = 0, nI = 0;
    for (const Block& b : m_blocks) {
        if (b.weight.size() != 0) return false;   // D_i != 1: not supported by the device kernel
        nB += static_cast<std::size_t>(b.B.size());
        nI += b.index.size();
    }
    if (nB > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || nI > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        return false;
    idx.reserve(nI);
    B.reserve(nB);
    for (const Block& b : m_blocks) {
        const int n = static_cast<int>(b.index.size());
        max_dim = std::max(max_dim, n);
        idx.insert(idx.end(), b.index.begin(), b.index.end());
        B.insert(B.end(), b.B.data(), b.B.data() + b.B.size());   // PrecondMatrix is row-major
        idx_off.push_back(static_cast<int>(idx.size()));
        b_off.push_back(static_cast<int>(B.size()));
    }
    fallback.assign(m_fallback.data(), m_fallback.data() + m_fallback.size());
    return true;
}

/* ------------------------------------------------------------------------- *
 *  Free helpers
 * ------------------------------------------------------------------------- */
void numericalScfJacobian(const Vector& x_star, const Vector& g_star,
                          const std::function<Vector(const Vector&)>& evalG, double h,
                          std::vector<double>& J, std::vector<double>* J_fwd)
{
    const int dim = static_cast<int>(x_star.size());
    J.assign(static_cast<size_t>(dim) * dim, 0.0);
    if (J_fwd) J_fwd->assign(static_cast<size_t>(dim) * dim, 0.0);
    Vector xp = x_star;
    for (int k = 0; k < dim; ++k) {
        xp(k) = x_star(k) + h;
        const Vector gp = evalG(xp);
        xp(k) = x_star(k) - h;
        const Vector gm = evalG(xp);
        xp(k) = x_star(k);
        for (int i = 0; i < dim; ++i) {
            const double d = (i == k) ? 1.0 : 0.0;
            J[static_cast<size_t>(i) * dim + k] = (gp(i) - gm(i)) / (2.0 * h) - d;
            if (J_fwd) (*J_fwd)[static_cast<size_t>(i) * dim + k] = (gp(i) - g_star(i)) / h - d;
        }
    }
}

std::vector<std::vector<int>> detectCovalentFragments(const std::vector<int>& atoms,
                                                      const Matrix& geometry)
{
    const int nat = static_cast<int>(atoms.size());
    std::vector<double> rcov(nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        const int z = atoms[i];
        rcov[i] = (z > 0 && z < static_cast<int>(Elements::CovalentRadius.size()))
                      ? Elements::CovalentRadius[z] : 1.5;
    }
    // Union-find over the bond graph.
    std::vector<int> parent(nat);
    for (int i = 0; i < nat; ++i) parent[i] = i;
    auto root = [&](int i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    };
    for (int i = 0; i < nat; ++i)
        for (int j = i + 1; j < nat; ++j) {
            const double cut = 1.25 * (rcov[i] + rcov[j]);
            if ((geometry.row(i) - geometry.row(j)).squaredNorm() < cut * cut) {
                const int a = root(i), b = root(j);
                if (a != b) parent[std::max(a, b)] = std::min(a, b);
            }
        }
    std::vector<int> frag_of(nat, -1);
    std::vector<std::vector<int>> frags;
    for (int i = 0; i < nat; ++i) {
        const int r = root(i);
        if (frag_of[r] < 0) { frag_of[r] = static_cast<int>(frags.size()); frags.emplace_back(); }
        frags[frag_of[r]].push_back(i);   // i ascending: each fragment is sorted
    }
    return frags;
}

/* ------------------------------------------------------------------------- *
 *  XTB::selectScfPreconditioner
 * ------------------------------------------------------------------------- */
/**
 * @brief Decide which seed replaces alpha*F in this SCF and set up the guard
 *        (Claude Generated, Oct 2026).
 *
 * Order: CURCUMA_SCF_PRECOND_FILE (developer override, dense M from a file; prototype
 * behaviour) > scf_precond (none | frag). With scf_precond none and no file the
 * function returns immediately and touches nothing, so the default SCF path is unchanged.
 * Conditions under which 'frag' is ignored (one warning per XTB object): not GFN2,
 * not the Broyden mode, neither the CPU loop nor a device loop that applies M itself
 * (@p device_loop, the fully device-resident CUDA loop; @p no_device_reason names why a GPU
 * path cannot), total charge != 0, open shell or odd electron count. With @p device_loop the
 * returned hook only signals "M is ready" (m_scf_precond_frag); the caller uploads M to the
 * device instead of installing the hook (Claude Generated, Oct 2026). The file matrix stays
 * CPU-only.
 *
 * @return false only for a hard error (unusable CURCUMA_SCF_PRECOND_FILE), with @p err set.
 */
bool XTB::selectScfPreconditioner(bool broyden_mode, bool cpu_loop, bool device_loop,
                                  const std::string& no_device_reason, int dim,
                                  std::function<Vector(const Vector&)>& hook,
                                  ScfPrecondGuard& guard, std::string& err)
{
    hook = nullptr;
    // Developer override: dense M from a file (behaviour of the Oct 2026 prototype).
    if (const char* pc_file = std::getenv("CURCUMA_SCF_PRECOND_FILE"); pc_file && pc_file[0] != '\0') {
        if (!broyden_mode) {
            CurcumaLogger::warn("CURCUMA_SCF_PRECOND_FILE ignored: needs the Broyden SCF mode");
        } else if (!cpu_loop) {
            CurcumaLogger::warn("CURCUMA_SCF_PRECOND_FILE ignored: a GPU SCF path is active "
                                "(the dense file matrix is implemented for the CPU SCF loop only)");
        } else {
            if (!loadScfPrecondMatrix(pc_file, dim, err)) {
                CurcumaLogger::error("CURCUMA_SCF_PRECOND_FILE: " + err);
                return false;
            }
            std::shared_ptr<const PrecondMatrix> M = m_scf_precond_M;
            hook = [M](const Vector& v) -> Vector { return (*M) * v; };
            guard = ScfPrecondGuard::fromEnv();
        }
        return true;
    }

    if (m_scf_precond.empty() || m_scf_precond == "none")
        return true;

    auto warn_once = [&](const std::string& why) {
        if (!m_scf_precond_warned)
            CurcumaLogger::warn("scf_precond '" + m_scf_precond + "' ignored: " + why
                                + "; scalar Broyden seed (scf_damping) is used");
        m_scf_precond_warned = true;
    };
    if (m_scf_precond != "frag") {
        warn_once("unknown value (none | frag)");
        return true;
    }
    if (m_method != MethodType::GFN2) { warn_once("implemented for GFN2 only"); return true; }
    if (!broyden_mode) { warn_once("needs -scf_mode broyden"); return true; }
    if (!cpu_loop && !device_loop) { warn_once(no_device_reason); return true; }
    if (m_charge != 0) { warn_once(fmt::format("total charge {} (neutral systems only)", m_charge)); return true; }
    const long long nel = std::llround(m_wfn.nocc);
    if (openShell() || (nel % 2) != 0) {
        warn_once("open shell or odd electron count (closed shells only)");
        return true;
    }

    // Reuse M unless the molecule changed or the atoms moved more than scf_precond_rebuild_rms.
    bool rebuild = !m_scf_precond_frag || m_scf_precond_frag->dim() != dim
                   || m_scf_precond_atoms != m_atoms
                   || m_scf_precond_geom.rows() != m_geometry.rows();
    double rms = 0.0;
    if (!rebuild) {
        rms = std::sqrt((m_geometry - m_scf_precond_geom).rowwise().squaredNorm().mean());
        if (m_scf_precond_rebuild_rms > 0.0 && rms > m_scf_precond_rebuild_rms) {
            rebuild = true;
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info_fmt("SCF precond: RMS displacement {:.3f} A > scf_precond_rebuild_rms "
                                        "{:.3f} A, rebuilding M", rms, m_scf_precond_rebuild_rms);
        }
    }
    if (rebuild) {
        std::string build_err;
        if (!buildFragmentPreconditioner(build_err)) {
            m_scf_precond_frag.reset();
            CurcumaLogger::warn("scf_precond frag: building M failed (" + build_err
                                + "); scalar Broyden seed for this SCF");
            return true;
        }
    }
    if (m_scf_precond_frag->blockCount() == 0) {
        // Every fragment is a single atom, an ion, odd-electron or failed: M would be alpha*I.
        warn_once("no fragment qualifies for a block (single atoms, ions, odd-electron, charged or too large fragments)");
        return true;
    }
    std::shared_ptr<const SchwarzPreconditioner> P = m_scf_precond_frag;
    hook = [P](const Vector& v) -> Vector { return P->apply(v); };
    if (const char* e = std::getenv("CURCUMA_SCF_PRECOND_GUARD"); e && e[0] != '\0')
        guard = ScfPrecondGuard::fromEnv();
    else {
        // Parsed every SCF (cheap); an unreadable spec warns once per object and spec.
        guard = ScfPrecondGuard::parse(m_scf_precond_guard, "scf_precond_guard",
                                       m_scf_precond_guard != m_scf_precond_guard_checked);
        m_scf_precond_guard_checked = m_scf_precond_guard;
    }
    return true;
}

/* ------------------------------------------------------------------------- *
 *  XTB::buildFragmentPreconditioner
 * ------------------------------------------------------------------------- */
namespace {

/// Claude Generated (Oct 2026): silences the logger while fragment SCFs run. CurcumaLogger keeps ONE process-wide
/// verbosity, so concurrent builds (Hessian or batch workers that each build an M) must not save and restore it
/// independently, else the verbosity can stay at 0: the first scope saves, the last one restores.
struct QuietLogScope {
    static std::mutex& mtx() { static std::mutex m; return m; }
    static int& depth() { static int d = 0; return d; }
    static int& saved() { static int v = 0; return v; }
    QuietLogScope()
    {
        std::lock_guard<std::mutex> lk(mtx());
        if (depth()++ == 0) { saved() = CurcumaLogger::get_verbosity(); CurcumaLogger::set_verbosity(0); }
    }
    ~QuietLogScope()
    {
        std::lock_guard<std::mutex> lk(mtx());
        if (--depth() == 0) CurcumaLogger::set_verbosity(saved());
    }
    QuietLogScope(const QuietLogScope&) = delete;
    QuietLogScope& operator=(const QuietLogScope&) = delete;
};
} // namespace

/**
 * @brief Build M = sum_f R_f^T (-J_f^-1) R_f + alpha on the rest (Claude Generated, Oct 2026).
 *
 * Fragments: detectCovalentFragments. A fragment gets its own block unless it is a single atom,
 * has an odd number of valence electrons as a neutral species, or carries a net charge in the
 * single-shot EEQ model (|sum q_EEQ| > 0.5 e) - then its indices keep the scalar alpha
 * (scf_damping). Every included fragment runs its own neutral closed-shell GFN2 SCF (FP64,
 * default settings, the parent's electronic temperature, gas phase) and then the central-
 * difference Jacobian of its SCC map at its converged state (m_scf_jac_capture, h = 1e-4).
 * The fragments run in parallel over the -threads budget (std::thread, one XTB per fragment,
 * each serial); every result goes to its own slot, so the outcome does not depend on the thread
 * count or the scheduling.
 *
 * @param err reason on failure
 * @return true if m_scf_precond_frag holds the new M
 */
bool XTB::buildFragmentPreconditioner(std::string& err)
{
    using clk = std::chrono::steady_clock;
    const auto t0 = clk::now();
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const int dim = nsh + 9 * nat;
    const double alpha = m_scf_damping;

    const std::vector<std::vector<int>> frags = detectCovalentFragments(m_atoms, m_geometry);

    // Net charge per fragment from the single-shot EEQ model (same model as the default guess).
    Vector q_eeq;
    {
        curcuma::dispersion::D4ChargeModel eeq;
        q_eeq = eeq.computeCharges(m_atoms, m_geometry * AA_TO_AU, static_cast<double>(m_charge));
    }
    const bool have_eeq = (q_eeq.size() == nat && q_eeq.allFinite());

    enum class Skip { None, SingleAtom, OddElectrons, Charged, TooLarge };
    std::vector<Skip> skip(frags.size(), Skip::None);
    std::vector<int> jobs;
    for (size_t f = 0; f < frags.size(); ++f) {
        const auto& a = frags[f];
        double nel = 0.0, q = 0.0;
        for (int i : a) { nel += m_wfn.n0_at(i); if (have_eeq) q += q_eeq(i); }
        if (a.size() == 1) skip[f] = Skip::SingleAtom;
        else if (m_scf_precond_max_frag_atoms > 0 && a.size() > static_cast<size_t>(m_scf_precond_max_frag_atoms)) skip[f] = Skip::TooLarge;
        else if ((std::llround(nel) % 2) != 0) skip[f] = Skip::OddElectrons;
        else if (!have_eeq || std::fabs(q) > 0.5) skip[f] = Skip::Charged;
        else jobs.push_back(static_cast<int>(f));
    }

    // Elements outside the measured set (H, C, N, O): allowed, one warning per process.
    {
        std::set<int> other;
        for (int f : jobs)
            for (int i : frags[f])
                if (m_atoms[i] != 1 && m_atoms[i] != 6 && m_atoms[i] != 7 && m_atoms[i] != 8)
                    other.insert(m_atoms[i]);
        static std::atomic<bool> warned_elements{ false };
        if (!other.empty() && !warned_elements.exchange(true)) {
            std::string list;
            for (int z : other)
                list += (list.empty() ? "" : " ")
                        + (z < static_cast<int>(Elements::ElementAbbr.size()) ? Elements::ElementAbbr[z]
                                                                              : std::to_string(z));
            CurcumaLogger::warn("scf_precond frag: fragments with elements outside the measured set "
                                "(H C N O): " + list + " - allowed, but not measured (machine-tested only)");
        }
    }

    // Per-fragment Jacobian -> block B_f = -J_f^-1.
    struct Result { bool ok = false; std::string why; PrecondMatrix B; };
    std::vector<Result> res(jobs.size());
    const double kT_parent = m_electronic_temp;
    auto run_job = [&](size_t j) {
        Result& r = res[j];
        const auto& atoms = frags[jobs[j]];
        try {
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
            frag.setScfGuess("eeq");
            frag.setIntraThreads(1);
            frag.setElectronicTemperature(kT_parent);
            frag.setD4AtmCutoff(m_d4_atm_cutoff);
            frag.setSto6gLegacy4sp(m_sto6g_legacy_4sp);
            frag.m_scf_jac_capture = true;
            if (!frag.InitialiseMolecule(sub)) { r.why = "fragment initialisation failed"; return; }
            frag.Calculation(false);
            if (frag.hasError() || !frag.m_scf_converged) { r.why = "fragment SCF failed"; return; }
            const int fdim = frag.m_basis.nsh + 9 * static_cast<int>(atoms.size());
            if (frag.m_scf_jac_result.rows() != fdim) { r.why = "fragment Jacobian missing"; return; }
            for (std::size_t k = 0; k < atoms.size(); ++k)
                if (frag.m_basis.nsh_at[k] != m_basis.nsh_at[atoms[k]]) { r.why = "shell count mismatch"; return; }
            Eigen::PartialPivLU<PrecondMatrix> lu(frag.m_scf_jac_result);
            r.B = -lu.inverse();
            if (!r.B.allFinite()) { r.why = "singular fragment Jacobian"; return; }
            r.ok = true;
        } catch (const std::exception& e) {
            r.why = e.what();
        }
    };

    {
        QuietLogScope quiet;   // fragment SCF tables are an implementation detail
        const int nthr = curcuma::intraParallelSuppressed()
                             ? 1 : std::max(1, std::min<int>(m_intra_threads, static_cast<int>(jobs.size())));
        std::atomic<size_t> next{ 0 };
        auto worker = [&]() {
            curcuma::SuppressIntraParallel serial;   // each fragment SCF stays serial
            for (size_t j; (j = next.fetch_add(1)) < jobs.size();) run_job(j);
        };
        if (nthr <= 1) {
            worker();
        } else {
            std::vector<std::thread> pool;
            for (int t = 0; t < nthr; ++t) pool.emplace_back(worker);
            for (auto& th : pool) th.join();
        }
    }
    const int saved_verbosity = CurcumaLogger::get_verbosity();   // restored by the scope above

    // Assemble the Schwarz blocks (non-overlapping, D_f = 1).
    auto P = std::make_shared<SchwarzPreconditioner>(dim, alpha);
    int n_block = 0, n_failed = 0, n_scalar_atoms = 0;
    std::string first_fail;
    for (size_t j = 0; j < jobs.size(); ++j) {
        const auto& atoms = frags[jobs[j]];
        if (!res[j].ok) {
            ++n_failed; n_scalar_atoms += static_cast<int>(atoms.size());
            if (first_fail.empty()) first_fail = res[j].why;
            continue;
        }
        SchwarzPreconditioner::Block b;
        // Local order of the fragment SCC vector: shells atom by atom, then dp (k = x,y,z) and
        // qp (k = xx,xy,yy,xz,yz,zz) component-wise over the fragment's atoms.
        for (int a : atoms)
            for (int s = 0; s < m_basis.nsh_at[a]; ++s) b.index.push_back(m_basis.ish_at[a] + s);
        for (int k = 0; k < 3; ++k)
            for (int a : atoms) b.index.push_back(nsh + k * nat + a);
        for (int k = 0; k < 6; ++k)
            for (int a : atoms) b.index.push_back(nsh + 3 * nat + k * nat + a);
        b.B = std::move(res[j].B);
        if (!P->addBlock(std::move(b))) {
            err = "inconsistent block for fragment " + std::to_string(jobs[j]);
            return false;
        }
        ++n_block;
    }
    int n_skip[5] = { 0, 0, 0, 0, 0 };
    for (size_t f = 0; f < frags.size(); ++f)
        if (skip[f] != Skip::None) {
            ++n_skip[static_cast<int>(skip[f])];
            n_scalar_atoms += static_cast<int>(frags[f].size());
        }
    if (n_skip[4] > 0) {
        static std::atomic<bool> warned_large{ false };
        if (!warned_large.exchange(true))
            CurcumaLogger::warn(fmt::format("scf_precond frag: {} fragment(s) with more than {} atoms keep the scalar seed "
                                            "(a fragment Jacobian costs 2*dim SCC-map evaluations, ~2 min for 212 atoms on one thread; raise scf_precond_max_fragment_atoms, 0 = no limit)",
                                            n_skip[4], m_scf_precond_max_frag_atoms));
    }
    if (n_failed > 0)
        CurcumaLogger::warn(fmt::format("scf_precond frag: {} of {} fragment Jacobians failed (first: {}); "
                                        "their atoms keep the scalar seed", n_failed, jobs.size(), first_fail));

    m_scf_precond_frag = P;
    m_scf_precond_geom = m_geometry;
    m_scf_precond_atoms = m_atoms;
    const double t_ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();
    if (saved_verbosity >= 2) {
        CurcumaLogger::info(fmt::format("SCF precond: frag, {} fragments, dim {}, built in {:.1f} ms",
                                        frags.size(), dim, t_ms));
        CurcumaLogger::info(fmt::format("SCF precond: {} blocks, scalar seed on {} atoms (single atoms {}, "
                                        "odd electrons {}, charged {}, too large {}, failed {}), memory {:.3f} MB",
                                        n_block, n_scalar_atoms, n_skip[1], n_skip[2], n_skip[3], n_skip[4], n_failed,
                                        8.0 * static_cast<double>(P->storedDoubles()) / 1048576.0));
    }
    return true;
}

/* ------------------------------------------------------------------------- *
 *  CURCUMA_SCF_JAC: Jacobian dump (moved from xtb_native.cpp, unchanged output)
 * ------------------------------------------------------------------------- */
/**
 * @brief Write the numerical Jacobian of the SCC fixed-point map plus preconditioner model data
 *        (Claude Generated, Oct 2026; debug path, active only with CURCUMA_SCF_JAC=<dir>).
 *
 * Scientific context: the native SCF is the fixed-point iteration x = G(x) on the SCC vector
 * x = [q_sh, dp_at, qp_at] (GFN2; GFN1: q_sh only) that the Broyden mixer mixes. With the
 * residual F(x) = G(x) - x the Jacobian is J = dF/dx = G'(x) - I. Its spectrum near x* decides
 * how fast a (preconditioned) quasi-Newton mixer converges; preconditioner candidates such as
 * M = eta (eta + gamma)^-1 (shell hardness eta, shell Coulomb matrix gamma) or an LDOS model
 * built from the MO coefficients approximate -J^-1. This routine writes J by central finite
 * differences of real evaluations of G, plus eta, gamma, the third-order Gamma, the GFN2
 * multipole kernels and the converged orbitals, overlap and occupations.
 *
 * Files (all raw float64 / int32 little endian, matrices row-major, shapes in jac_meta.json):
 *   x_star.bin, g_star.bin (dim), J.bin and J_fwd.bin (dim x dim, J[i,k] = dG_i/dx_k - delta_ik),
 *   lin_delta.bin, lin_g.bin (nlin x dim: random perturbation delta and G(x*+delta)),
 *   eta_sh.bin (nsh), gamma_sh.bin (nsh x nsh), third_sh.bin (nsh), [GFN2] mp_*.bin,
 *   eps.bin, focc.bin (nao), C.bin (nao x nao, C[mu,i]), S.bin (nao x nao), ao2sh/ao2at (int32).
 * Environment: CURCUMA_SCF_JAC_H (step h, default 1e-4), CURCUMA_SCF_JAC_LINTEST (number of
 * linearity-test perturbations, default 5), CURCUMA_SCF_JAC_LINSCALE (|delta|/|x*|, default 1e-3).
 *
 * @param dir     output directory (created if missing)
 * @param x_star  SCC input vector of the converged iteration
 * @param g_loop  SCC output vector of that iteration, as the SCF loop computed it
 * @param e_scc   SCC energy of that iteration (Eh)
 * @param evalG   one pure evaluation of the fixed-point map (modifies m_wfn / m_pot; the caller restores)
 */
void XTB::writeScfJacobianDump(const std::string& dir, const Vector& x_star, const Vector& g_loop,
                               double e_scc, const std::function<Vector(const Vector&)>& evalG)
{
    namespace fs = std::filesystem;
    using clk = std::chrono::steady_clock;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        CurcumaLogger::warn("CURCUMA_SCF_JAC: cannot create directory " + dir + ": " + ec.message());
        return;
    }
    const int  nat  = m_atomcount;
    const int  nsh  = m_basis.nsh;
    const int  nao  = m_basis.nao;
    const int  dim  = static_cast<int>(x_star.size());
    const bool gfn2 = (m_method == MethodType::GFN2);

    double h = 1.0e-4;
    if (const char* e = std::getenv("CURCUMA_SCF_JAC_H")) {
        const double v = std::strtod(e, nullptr);
        if (v > 0.0 && std::isfinite(v)) h = v;
    }
    int nlin = 5;
    if (const char* e = std::getenv("CURCUMA_SCF_JAC_LINTEST")) {
        const int v = std::atoi(e);
        if (v >= 0) nlin = v;
    }
    double lin_scale = 1.0e-3;   // |delta| / |x*| of the linearity-test perturbations
    if (const char* e = std::getenv("CURCUMA_SCF_JAC_LINSCALE")) {
        const double v = std::strtod(e, nullptr);
        if (v > 0.0 && std::isfinite(v)) lin_scale = v;
    }

    bool io_ok = true;
    auto writeD = [&](const std::string& name, const double* p, size_t n) {
        std::ofstream f(fs::path(dir) / name, std::ios::binary);
        if (!f) { io_ok = false; return; }
        f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n * sizeof(double)));
        if (!f) io_ok = false;
    };
    auto writeI = [&](const std::string& name, const std::vector<int>& v) {
        std::vector<int32_t> b(v.begin(), v.end());
        std::ofstream f(fs::path(dir) / name, std::ios::binary);
        if (!f) { io_ok = false; return; }
        f.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size() * sizeof(int32_t)));
        if (!f) io_ok = false;
    };
    auto writeVec = [&](const std::string& name, const Vector& v) {
        writeD(name, v.data(), static_cast<size_t>(v.size()));
    };
    // Any Eigen matrix (row- or column-major storage) as row-major M[i,j].
    auto writeMat = [&](const std::string& name, const auto& M) {
        const Eigen::Index r = M.rows(), c = M.cols();
        std::vector<double> b(static_cast<size_t>(r * c));
        for (Eigen::Index i = 0; i < r; ++i)
            for (Eigen::Index j = 0; j < c; ++j) b[static_cast<size_t>(i * c + j)] = M(i, j);
        writeD(name, b.data(), b.size());
    };

    // (4) LDOS data of the converged state (m_wfn on entry: orbitals of the last iteration, whose
    // input was x*). Fermi level recomputed with the bisection of solveEigen (closed shell, T > 0).
    const double kT = m_electronic_temp * 3.166808e-6;   // K -> Eh, same constant as solveEigen
    double mu_fermi = std::nan("");
    if (!openShell() && m_electronic_temp > 0.0 && m_wfn.eps.size() == nao) {
        double mu_lo = m_wfn.eps.minCoeff() - 1.0, mu_hi = m_wfn.eps.maxCoeff() + 1.0;
        for (int bisect = 0; bisect < 100; ++bisect) {
            const double mu = 0.5 * (mu_lo + mu_hi);
            double n_sum = 0.0;
            for (int i = 0; i < nao; ++i)
                n_sum += 2.0 / (1.0 + std::exp(std::min((m_wfn.eps(i) - mu) / kT, 500.0)));
            if (n_sum > m_wfn.nocc) mu_hi = mu; else mu_lo = mu;
            if (mu_hi - mu_lo < 1e-14) break;
        }
        mu_fermi = 0.5 * (mu_lo + mu_hi);
    }
    const bool have_orb = (m_wfn.eps.size() == nao && m_wfn.C.rows() == nao && m_wfn.C.cols() == nao);
    if (have_orb) {
        writeVec("eps.bin", m_wfn.eps);
        writeMat("C.bin", m_wfn.C);
        if (m_wfn.focc.size() == nao) writeVec("focc.bin", m_wfn.focc);
    }
    writeMat("S.bin", m_S);
    writeI("ao2sh.bin", m_basis.ao2sh);
    writeI("ao2at.bin", m_basis.ao2at);

    // (3) Model quantities. eta_sh: shell hardness as it enters the second-order term (diagonal of
    // gamma, hubbard(Z) * shell_hubbard(l, Z), xtb_coulomb.hpp); gamma_sh: m_gamma as used in
    // v_sh += gamma q_sh; third_sh: Gamma per shell with E3 = sum_s q_s^3 Gamma_s / 3 (GFN2 shell-
    // resolved; GFN1 atom-resolved, the atom value repeated on its shells, E3 over q_at).
    const auto cmeth = gfn2 ? coulomb::Method::GFN2 : coulomb::Method::GFN1;
    Vector eta_sh(nsh), third_sh(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int z = m_basis.z[m_basis.sh2at[s]];
        const int l = m_basis.ang_sh[s];
        eta_sh(s)   = coulomb::shell_hardness(cmeth, z, l);
        third_sh(s) = gfn2 ? coulomb::shell_hubbard_deriv_gfn2(z, l) : coulomb::atom_hubbard_deriv_gfn1(z);
    }
    writeVec("eta_sh.bin", eta_sh);
    writeMat("gamma_sh.bin", m_gamma);
    writeVec("third_sh.bin", third_sh);
    const bool have_mp = gfn2 && m_mp_initialized && static_cast<int>(m_mp_dkernel.size()) == nat
                         && m_mp_amat_sd[0].rows() == nat;
    if (have_mp) {
        writeD("mp_dkernel.bin", m_mp_dkernel.data(), m_mp_dkernel.size());
        writeD("mp_qkernel.bin", m_mp_qkernel.data(), m_mp_qkernel.size());
        writeD("mp_mrad.bin", m_mp_mrad.data(), m_mp_mrad.size());
        Eigen::MatrixXd sd(3 * nat, nat), dd(9 * nat, nat), sq(6 * nat, nat);
        for (int k = 0; k < 3; ++k) sd.middleRows(k * nat, nat) = m_mp_amat_sd[k];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) dd.middleRows((3 * a + b) * nat, nat) = m_mp_amat_dd[a][b];
        for (int k = 0; k < 6; ++k) sq.middleRows(k * nat, nat) = m_mp_amat_sq[k];
        writeMat("mp_amat_sd.bin", sd);
        writeMat("mp_amat_dd.bin", dd);
        writeMat("mp_amat_sq.bin", sq);
    }

    // (1) x*, G(x*) and the purity check: two evaluations in a row, and against the loop's output.
    const Vector g1 = evalG(x_star);
    const Vector g2 = evalG(x_star);
    const double noise_repeat = (g1 - g2).cwiseAbs().maxCoeff();
    const double noise_loop   = (g_loop.size() == dim) ? (g1 - g_loop).cwiseAbs().maxCoeff() : std::nan("");
    const double resid_max    = (g1 - x_star).cwiseAbs().maxCoeff();
    writeVec("x_star.bin", x_star);
    writeVec("g_star.bin", g1);

    // (2) J[i,k] = (G_i(x*+h e_k) - G_i(x*-h e_k)) / 2h - delta_ik, row-major (i = output row).
    // J_fwd from the same evaluations: (G_i(x*+h e_k) - G_i(x*)) / h - delta_ik.
    const auto t_jac0 = clk::now();
    std::vector<double> J, Jf;
    numericalScfJacobian(x_star, g1, evalG, h, J, &Jf);
    const double t_jac_ms = std::chrono::duration<double, std::milli>(clk::now() - t_jac0).count();
    writeD("J.bin", J.data(), J.size());
    writeD("J_fwd.bin", Jf.data(), Jf.size());

    // Linearity test data: random perturbations of norm lin_scale * |x*| (default 1e-3) and real evaluations G(x*+delta).
    std::mt19937_64 rng(20261010ULL);
    std::normal_distribution<double> nd(0.0, 1.0);
    const double dnorm = lin_scale * (x_star.norm() > 0.0 ? x_star.norm() : 1.0);
    std::vector<double> lin_d(static_cast<size_t>(nlin) * dim), lin_g(static_cast<size_t>(nlin) * dim);
    for (int t = 0; t < nlin; ++t) {
        Vector d(dim);
        for (int i = 0; i < dim; ++i) d(i) = nd(rng);
        d *= dnorm / d.norm();
        const Vector g = evalG(x_star + d);
        for (int i = 0; i < dim; ++i) {
            lin_d[static_cast<size_t>(t) * dim + i] = d(i);
            lin_g[static_cast<size_t>(t) * dim + i] = g(i);
        }
    }
    if (nlin > 0) {
        writeD("lin_delta.bin", lin_d.data(), lin_d.size());
        writeD("lin_g.bin", lin_g.data(), lin_g.size());
    }

    // Meta data.
    nlohmann::json meta;
    meta["method"] = gfn2 ? "gfn2" : "gfn1";
    meta["natoms"] = nat; meta["nsh"] = nsh; meta["nao"] = nao; meta["dim"] = dim;
    meta["atomic_numbers"] = m_basis.z;
    meta["shell_atom"] = m_basis.sh2at;
    meta["shell_l"] = m_basis.ang_sh;
    nlohmann::json blocks;
    blocks["q_sh"] = { { "offset", 0 }, { "size", nsh },
                       { "layout", "shell charge q_s = n0_s - n_s (e), shell order as shell_atom" } };
    if (gfn2) {
        blocks["dp"] = { { "offset", nsh }, { "size", 3 * nat },
                         { "layout", "atomic dipole, index offset + k*natoms + atom, k = x,y,z (a.u.)" } };
        blocks["qp"] = { { "offset", nsh + 3 * nat }, { "size", 6 * nat },
                         { "layout", "atomic traceless quadrupole, index offset + k*natoms + atom, "
                                     "k = xx,xy,yy,xz,yz,zz (a.u.)" } };
    }
    meta["blocks"] = blocks;
    meta["h"] = h;
    meta["e_scc_last_iteration_Eh"] = e_scc;
    meta["electronic_temperature_K"] = m_electronic_temp;
    meta["kT_Eh"] = kT;
    meta["fermi_level_Eh"] = std::isfinite(mu_fermi) ? nlohmann::json(mu_fermi) : nlohmann::json(nullptr);
    meta["nelectrons"] = m_wfn.nocc;
    meta["open_shell"] = openShell();
    meta["noise_G_repeat_maxabs"] = noise_repeat;
    meta["noise_G_vs_loop_maxabs"] = noise_loop;
    meta["residual_G_minus_x_maxabs"] = resid_max;
    meta["jacobian_time_ms"] = t_jac_ms;
    meta["jacobian_evaluations"] = 2 * dim;
    meta["lin_test_count"] = nlin;
    meta["lin_test_delta_norm"] = dnorm;
    meta["J_convention"] = "J.bin row-major float64 dim x dim; J[i,k] = (G_i(x*+h e_k) - G_i(x*-h e_k))/(2h) "
                           "- delta_ik = dF_i/dx_k with F = G(x) - x; row i = output component, column k = "
                           "perturbed input component. J_fwd.bin: same with the forward difference "
                           "(G_i(x*+h e_k) - G_i(x*))/h - delta_ik.";
    meta["files"] = {
        { "x_star.bin", "float64 [dim], converged SCC input of the last iteration" },
        { "g_star.bin", "float64 [dim], G(x*) evaluated by the dump (first of the two purity evaluations)" },
        { "J.bin", "float64 [dim, dim]" }, { "J_fwd.bin", "float64 [dim, dim]" },
        { "lin_delta.bin", "float64 [lin_test_count, dim], random delta" },
        { "lin_g.bin", "float64 [lin_test_count, dim], G(x*+delta)" },
        { "eta_sh.bin", "float64 [nsh], shell hardness = diagonal of gamma (Eh)" },
        { "gamma_sh.bin", "float64 [nsh, nsh], shell Coulomb matrix, v_sh = gamma q_sh (Eh)" },
        { "third_sh.bin", gfn2 ? "float64 [nsh], Gamma_s, E3 = sum_s q_s^3 Gamma_s / 3, v_s = Gamma_s q_s^2"
                               : "float64 [nsh], atom Gamma_A on each shell of A, E3 = sum_A q_A^3 Gamma_A / 3" },
        { "eps.bin", "float64 [nao], orbital energies (Eh)" },
        { "focc.bin", "float64 [nao], occupations (0..2)" },
        { "C.bin", "float64 [nao, nao], C[mu, i], AO mu, MO i; C^T S C = I" },
        { "S.bin", "float64 [nao, nao], AO overlap" },
        { "ao2sh.bin", "int32 [nao], AO -> shell" }, { "ao2at.bin", "int32 [nao], AO -> atom" },
    };
    if (have_mp) {
        meta["files"]["mp_dkernel.bin"] = "float64 [natoms], on-site dipole XC kernel";
        meta["files"]["mp_qkernel.bin"] = "float64 [natoms], on-site quadrupole XC kernel";
        meta["files"]["mp_mrad.bin"]    = "float64 [natoms], multipole damping radius";
        meta["files"]["mp_amat_sd.bin"] = "float64 [3, natoms, natoms], m_mp_amat_sd[k](i, j)";
        meta["files"]["mp_amat_dd.bin"] = "float64 [3, 3, natoms, natoms], m_mp_amat_dd[a][b](i, j)";
        meta["files"]["mp_amat_sq.bin"] = "float64 [6, natoms, natoms], m_mp_amat_sq[k](i, j)";
    }
    meta["have_orbitals"] = have_orb;
    {
        std::ofstream f(fs::path(dir) / "jac_meta.json");
        if (f) f << meta.dump(2) << "\n"; else io_ok = false;
    }
    if (!io_ok)
        CurcumaLogger::warn("CURCUMA_SCF_JAC: writing a file in " + dir + " failed");
    fmt::print(stderr, "SCFJAC dim {} h {:.1e} noise_repeat {:.3e} noise_vs_loop {:.3e} resid {:.3e} "
                       "t_jac_ms {:.1f} evals {} dir {}\n",
               dim, h, noise_repeat, noise_loop, resid_max, t_jac_ms, 2 * dim, dir);
}

/* ------------------------------------------------------------------------- *
 *  CURCUMA_SCF_PRECOND_FILE: dense M from a file (moved from xtb_native.cpp)
 * ------------------------------------------------------------------------- */
/**
 * @brief Load the Broyden preconditioner matrix M from a raw float64 file
 *        (Claude Generated, Oct 2026; debug path, active only with CURCUMA_SCF_PRECOND_FILE).
 *
 * The modified Broyden mixer seeds its inverse-Jacobian estimate with the scalar step alpha*F.
 * With a step hook, M*F replaces alpha*F (BroydenMixer::setStepHook); a good choice is
 * M ~ -J^-1 with J = G' - I the Jacobian of the SCC residual F(x) = G(x) - x
 * (G' as written to J.bin by CURCUMA_SCF_JAC). M = alpha*I reproduces the scalar path.
 *
 * File format: dim*dim little-endian float64, row-major (entry (i,k) = row i, column k), where
 * dim is the length of the SCC vector x = [q_sh | dp | qp] (packSCC in Calculation()). The matrix
 * is cached per XTB object and read again only when the path or the file mtime change.
 *
 * @param path file name
 * @param dim  expected dimension (length of the SCC vector)
 * @param err  error text on failure
 * @return true if m_scf_precond_M holds a dim x dim matrix
 */
bool XTB::loadScfPrecondMatrix(const std::string& path, int dim, std::string& err)
{
    std::error_code ec;
    const auto mtime = std::filesystem::last_write_time(path, ec);
    if (ec) {
        err = "cannot stat " + path + ": " + ec.message();
        return false;
    }
    const long long mtime_ticks = static_cast<long long>(mtime.time_since_epoch().count());
    if (!(m_scf_precond_M && m_scf_precond_path == path && m_scf_precond_mtime == mtime_ticks)) {
        const auto bytes = std::filesystem::file_size(path, ec);
        if (ec) {
            err = "cannot read size of " + path + ": " + ec.message();
            return false;
        }
        if (bytes == 0 || bytes % sizeof(double) != 0) {
            err = fmt::format("{}: size {} bytes is not a multiple of 8 (dim*dim float64 expected)",
                              path, bytes);
            return false;
        }
        const std::uintmax_t n = bytes / sizeof(double);
        const auto n_dim = static_cast<std::uintmax_t>(std::llround(std::sqrt(static_cast<double>(n))));
        if (n_dim * n_dim != n) {
            err = fmt::format("{}: {} float64 entries is not a square dim*dim matrix", path, n);
            return false;
        }
        auto M = std::make_shared<PrecondMatrix>(static_cast<Eigen::Index>(n_dim),
                                                 static_cast<Eigen::Index>(n_dim));
        std::ifstream in(path, std::ios::binary);
        if (!in || !in.read(reinterpret_cast<char*>(M->data()), static_cast<std::streamsize>(bytes))) {
            err = "reading " + path + " failed";
            return false;
        }
        // The file is little-endian; this build reads it natively (x86-64/aarch64 are little-endian).
        static_assert(sizeof(double) == 8, "float64 expected");
        m_scf_precond_M = std::move(M);
        m_scf_precond_path = path;
        m_scf_precond_mtime = mtime_ticks;
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format("SCF precond: M loaded (dim {}) from {}", n_dim, path));
    }
    if (m_scf_precond_M->rows() != dim) {
        err = fmt::format("{}: matrix dim {} does not match the SCC vector length {}",
                          path, m_scf_precond_M->rows(), dim);
        return false;
    }
    return true;
}

} // namespace curcuma::xtb
