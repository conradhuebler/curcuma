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

bool SchwarzPreconditioner::addBlock(Block block, const std::vector<char>& covers)
{
    const int n = static_cast<int>(block.index.size());
    if (n == 0 || block.B.rows() != n || block.B.cols() != n) return false;
    if (block.weight.size() != 0 && block.weight.size() != n) return false;
    if (!covers.empty() && static_cast<int>(covers.size()) != n) return false;
    for (int i : block.index)
        if (i < 0 || i >= m_dim) return false;
    for (int k = 0; k < n; ++k)   // covered: C has no entry here
        if (covers.empty() || covers[k]) m_fallback(block.index[k]) = 0.0;
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

namespace {
/// Pyykko-Atsumi single-bond covalent radius (Angstrom), 1.5 for elements outside the table.
double covalentRadius(int z)
{
    return (z > 0 && z < static_cast<int>(Elements::CovalentRadius.size())) ? Elements::CovalentRadius[z] : 1.5;
}
} // namespace

std::vector<std::vector<int>> detectCovalentFragments(const std::vector<int>& atoms,
                                                      const Matrix& geometry)
{
    const int nat = static_cast<int>(atoms.size());
    std::vector<double> rcov(nat, 0.0);
    for (int i = 0; i < nat; ++i) rcov[i] = covalentRadius(atoms[i]);
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

ScfPrecondSegment::ScfPrecondSegment() = default;
ScfPrecondSegment::~ScfPrecondSegment() = default;
ScfPrecondSegment::ScfPrecondSegment(const ScfPrecondSegment&) = default;
ScfPrecondSegment::ScfPrecondSegment(ScfPrecondSegment&&) noexcept = default;
ScfPrecondSegment& ScfPrecondSegment::operator=(const ScfPrecondSegment&) = default;
ScfPrecondSegment& ScfPrecondSegment::operator=(ScfPrecondSegment&&) noexcept = default;

std::vector<ScfPrecondSegment> segmentFragment(const std::vector<int>& atoms, const Matrix& geometry,
                                               const std::vector<int>& frag, int segment_heavy_atoms,
                                               double radius)
{
    std::vector<ScfPrecondSegment> segs;
    const int nf = static_cast<int>(frag.size());
    if (nf == 0 || segment_heavy_atoms <= 0) return segs;

    // Local positions 0..nf-1 refer to frag[]. Bond graph inside the fragment (same criterion as
    // detectCovalentFragments), neighbours in ascending order.
    std::vector<std::vector<int>> nbr(nf);
    for (int a = 0; a < nf; ++a)
        for (int b = a + 1; b < nf; ++b) {
            const double cut = 1.25 * (covalentRadius(atoms[frag[a]]) + covalentRadius(atoms[frag[b]]));
            if ((geometry.row(frag[a]) - geometry.row(frag[b])).squaredNorm() < cut * cut) {
                nbr[a].push_back(b);
                nbr[b].push_back(a);
            }
        }
    std::vector<char> heavy(nf, 0);
    bool any_heavy = false;
    for (int a = 0; a < nf; ++a) { heavy[a] = atoms[frag[a]] > 1; any_heavy = any_heavy || heavy[a]; }
    if (!any_heavy) std::fill(heavy.begin(), heavy.end(), 1);

    // Ring bonds = non-bridges of the heavy-atom bond graph (Tarjan low-link, iterative DFS).
    // bridge_to[a] holds the neighbours b with a-b a bridge.
    std::vector<std::set<int>> bridge_to(nf);
    {
        std::vector<int> tin(nf, -1), low(nf, 0), parent(nf, -1);
        std::vector<std::size_t> it(nf, 0);
        int timer = 0;
        for (int s = 0; s < nf; ++s) {
            if (!heavy[s] || tin[s] >= 0) continue;
            std::vector<int> stack{ s };
            tin[s] = low[s] = timer++;
            while (!stack.empty()) {
                const int a = stack.back();
                if (it[a] < nbr[a].size()) {
                    const int b = nbr[a][it[a]++];
                    if (!heavy[b] || b == parent[a]) continue;
                    if (tin[b] >= 0) { low[a] = std::min(low[a], tin[b]); continue; }
                    parent[b] = a;
                    tin[b] = low[b] = timer++;
                    stack.push_back(b);
                } else {
                    stack.pop_back();
                    const int p = parent[a];
                    if (p >= 0) {
                        low[p] = std::min(low[p], low[a]);
                        if (low[a] > tin[p]) { bridge_to[p].insert(a); bridge_to[a].insert(p); }
                    }
                }
            }
        }
    }
    // A heavy-heavy bond may be cut only if it is a bridge (not in a ring) and a single bond:
    // r_ij / (r_cov,i + r_cov,j) >= 0.95 (Pyykko-Atsumi single-bond radii; measured ratios: single
    // ~1.0, amide C-N ~0.94, aromatic ~0.92, C=O ~0.88). Atoms joined by uncuttable bonds form one
    // unit (union-find) that is never split: a hydrogen cap there would not give a closed shell.
    auto cuttable = [&](int a, int b) {
        if (!bridge_to[a].count(b)) return false;
        const double r = (geometry.row(frag[a]) - geometry.row(frag[b])).norm();
        return r / (covalentRadius(atoms[frag[a]]) + covalentRadius(atoms[frag[b]])) >= 0.95;
    };
    std::vector<int> unit(nf);
    for (int a = 0; a < nf; ++a) unit[a] = a;
    auto root = [&](int i) {
        while (unit[i] != i) { unit[i] = unit[unit[i]]; i = unit[i]; }
        return i;
    };
    for (int a = 0; a < nf; ++a)
        if (heavy[a])
            for (int b : nbr[a])
                if (b > a && heavy[b] && !cuttable(a, b)) {
                    const int ra = root(a), rb = root(b);
                    if (ra != rb) unit[std::max(ra, rb)] = std::min(ra, rb);   // root = smallest index
                }
    std::vector<std::vector<int>> members(nf);   // heavy atoms per unit root, ascending
    for (int a = 0; a < nf; ++a)
        if (heavy[a]) members[root(a)].push_back(a);

    // Units in breadth-first order over the bonds between units (start at the smallest index,
    // neighbours in ascending order of their member atoms); an unreachable unit starts a new search.
    std::vector<int> order;   // unit roots
    std::vector<char> seen(nf, 0);
    for (int s = 0; s < nf; ++s) {
        if (!heavy[s] || seen[root(s)]) continue;
        std::vector<int> queue{ root(s) };
        seen[root(s)] = 1;
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const int u = queue[q];
            order.push_back(u);
            std::vector<int> next;
            for (int a : members[u])
                for (int b : nbr[a])
                    if (heavy[b] && !seen[root(b)]) { seen[root(b)] = 1; next.push_back(root(b)); }
            std::sort(next.begin(), next.end());   // root = smallest member index
            queue.insert(queue.end(), next.begin(), next.end());
        }
    }

    // Every light atom (H) joins its nearest heavy atom.
    std::vector<int> host(nf, -1);
    for (int a = 0; a < nf; ++a) {
        if (heavy[a]) { host[a] = a; continue; }
        double best = std::numeric_limits<double>::max();
        for (int b = 0; b < nf; ++b)
            if (heavy[b]) {
                const double d = (geometry.row(frag[a]) - geometry.row(frag[b])).squaredNorm();
                if (d < best) { best = d; host[a] = b; }
            }
    }
    std::vector<std::vector<int>> light_of(nf);
    for (int a = 0; a < nf; ++a)
        if (!heavy[a]) light_of[host[a]].push_back(a);

    const double r2 = radius * radius;
    const double rcov_h = covalentRadius(1);
    for (std::size_t g0 = 0; g0 < order.size();) {
        // Fill the segment with whole units until it holds >= segment_heavy_atoms heavy atoms.
        std::vector<int> own_heavy;
        std::size_t g1 = g0;
        while (g1 < order.size() && static_cast<int>(own_heavy.size()) < segment_heavy_atoms) {
            own_heavy.insert(own_heavy.end(), members[order[g1]].begin(), members[order[g1]].end());
            ++g1;
        }
        g0 = g1;
        std::vector<char> own(nf, 0), dom(nf, 0);
        for (int a : own_heavy) {
            own[a] = 1;
            for (int l : light_of[a]) own[l] = 1;
        }
        // Domain: every unit with a heavy atom within radius of an own heavy atom, with the light atoms.
        for (int b = 0; b < nf; ++b) {
            if (!heavy[b] || dom[b]) continue;
            for (int o : own_heavy)
                if ((geometry.row(frag[b]) - geometry.row(frag[o])).squaredNorm() <= r2) {
                    for (int m : members[root(b)]) {
                        dom[m] = 1;
                        for (int l : light_of[m]) dom[l] = 1;
                    }
                    break;
                }
        }
        ScfPrecondSegment seg;
        std::vector<Eigen::RowVector3d> caps;
        for (int a = 0; a < nf; ++a) {
            if (own[a]) seg.own.push_back(frag[a]);
            if (!dom[a]) continue;
            seg.domain.push_back(frag[a]);
            if (!heavy[a]) continue;
            for (int k : nbr[a]) {   // cut bonds: heavy-heavy, between units, hence cuttable single bonds
                if (dom[k] || !heavy[k]) continue;
                const Eigen::RowVector3d ra = geometry.row(frag[a]);
                const Eigen::RowVector3d u = (geometry.row(frag[k]) - geometry.row(frag[a])).normalized();
                caps.push_back(ra + u * (covalentRadius(atoms[frag[a]]) + rcov_h));
            }
        }
        seg.caps = Matrix(static_cast<Eigen::Index>(caps.size()), 3);
        for (std::size_t c = 0; c < caps.size(); ++c) seg.caps.row(static_cast<Eigen::Index>(c)) = caps[c];
        segs.push_back(std::move(seg));
    }
    return segs;
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
 * Conditions under which 'frag' is ignored (one warning per XTB object): neither GFN1 nor GFN2
 * (GFN1: CPU loop and the backends with supportsResidentPreconditionerGfn1, i.e. ROCm; its SCC vector is q_sh alone),
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
    if (m_method != MethodType::GFN2 && m_method != MethodType::GFN1) { warn_once("implemented for GFN1 and GFN2 only"); return true; }
    if (!broyden_mode) { warn_once("needs -scf_mode broyden"); return true; }
    if (!cpu_loop && !device_loop) { warn_once(no_device_reason); return true; }
    // GFN1 device loop only where the backend says so (ROCm WP5, Claude Generated Oct 2026).
    if (!cpu_loop && m_method != MethodType::GFN2
        && !(m_gpu_scf && m_gpu_scf->supportsResidentPreconditionerGfn1())) { warn_once("the device preconditioner is implemented for GFN2 only (GFN1: CPU SCF loop)"); return true; }
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
    if (rebuild && adoptSharedPreconditioner()) {
        rebuild = false;   // another XTB object of this process built M for the same molecule (Hessian workers)
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

/// Setters with input checks (Claude Generated, Oct 2026): one warning per object for invalid values.
void XTB::setScfPrecondMaxFragmentAtoms(int n)
{
    if (n < 0) {
        if (!m_scf_precond_input_warned)
            CurcumaLogger::warn(fmt::format("scf_precond_max_fragment_atoms {} is negative; the default 100 is used (0 = no limit)", n));
        m_scf_precond_input_warned = true;
        m_scf_precond_max_frag_atoms = 100;
        return;
    }
    m_scf_precond_max_frag_atoms = n;
}

void XTB::setScfPrecondSegmentRadius(double r)
{
    if (!(r > 0.0)) {
        if (!m_scf_precond_input_warned)
            CurcumaLogger::warn(fmt::format("scf_precond_segment_radius {} must be positive; {:.2f} A is kept", r, m_scf_precond_segment_radius));
        m_scf_precond_input_warned = true;
        return;
    }
    m_scf_precond_segment_radius = r;
}

/* ------------------------------------------------------------------------- *
 *  Process-wide sharing of M between XTB objects (Claude Generated, Oct 2026)
 * ------------------------------------------------------------------------- */
namespace {
/// One shared M: the key holds every input M depends on besides the geometry; the entry lives only
/// as long as some XTB object holds M (weak_ptr), so the cache never grows beyond the live objects.
struct SharedPrecondEntry {
    std::string key;
    Matrix geom;
    std::weak_ptr<const SchwarzPreconditioner> M;
};
std::mutex& sharedPrecondMutex() { static std::mutex m; return m; }
std::vector<SharedPrecondEntry>& sharedPrecondCache() { static std::vector<SharedPrecondEntry> c; return c; }
} // namespace

/// Key of M without the geometry: method, charge, atoms and every parameter that enters the build.
std::string XTB::scfPrecondShareKey() const
{
    std::string k = fmt::format("{}|{}|{:.17g}|{:.17g}|{}|{:.17g}|{}|{}|{:.17g}|{}|", static_cast<int>(m_method), m_charge,
                                m_electronic_temp, m_d4_atm_cutoff, m_sto6g_legacy_4sp ? 1 : 0, m_scf_damping,
                                m_scf_precond_max_frag_atoms, m_scf_precond_segment_atoms, m_scf_precond_segment_radius,
                                m_scf_precond_jacobian);
    for (int z : m_atoms) k += std::to_string(z) + ",";
    return k;
}

/**
 * @brief Take over an M that another live XTB object of this process built for the same molecule
 *        (Claude Generated, Oct 2026).
 *
 * The semi-numerical Hessian gives every worker thread its own XTB object; without sharing each of
 * them built its own M at its first displaced geometry. Only objects marked with setScfPrecondShare(true)
 * (the Hessian's workers and its reference calculator) take part, so other parallel users of XTB
 * (ConfSearch, batch optimisations) keep their own M and their results do not depend on thread timing. An entry is used if its key matches
 * (scfPrecondShareKey) and the RMS displacement from its build geometry is at most
 * scf_precond_rebuild_rms, or 0.01 A when that is 0 (Hessian steps are 0.005 Bohr on one atom).
 * M is immutable (shared_ptr<const>), so sharing is thread-safe; m_scf_precond_geom becomes the
 * build geometry of the shared M, so the rebuild rule keeps measuring from there.
 */
bool XTB::adoptSharedPreconditioner()
{
    if (!m_scf_precond_share) return false;
    const std::string key = scfPrecondShareKey();
    const double tol = (m_scf_precond_rebuild_rms > 0.0) ? m_scf_precond_rebuild_rms : 0.01;
    std::lock_guard<std::mutex> lk(sharedPrecondMutex());
    auto& cache = sharedPrecondCache();
    cache.erase(std::remove_if(cache.begin(), cache.end(),
                               [](const SharedPrecondEntry& e) { return e.M.expired(); }), cache.end());
    for (const SharedPrecondEntry& e : cache) {
        if (e.key != key || e.geom.rows() != m_geometry.rows()) continue;
        const double rms = std::sqrt((m_geometry - e.geom).rowwise().squaredNorm().mean());
        if (rms > tol) continue;
        std::shared_ptr<const SchwarzPreconditioner> M = e.M.lock();
        if (!M) continue;
        m_scf_precond_frag = M;
        m_scf_precond_geom = e.geom;
        m_scf_precond_atoms = m_atoms;
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info_fmt("SCF precond: M shared with another calculator of this process (RMS {:.4f} A from its build geometry)", rms);
        return true;
    }
    return false;
}

/// Register a freshly built M for adoptSharedPreconditioner (Claude Generated, Oct 2026).
void XTB::publishSharedPreconditioner()
{
    if (!m_scf_precond_share) return;
    std::lock_guard<std::mutex> lk(sharedPrecondMutex());
    sharedPrecondCache().push_back({ scfPrecondShareKey(), m_scf_precond_geom, m_scf_precond_frag });
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

/// Claude Generated (Oct 2026): one local cluster of the fragment preconditioner build - real atoms (global
/// indices, block order) plus hydrogen cap positions; own[k] marks the real atoms whose rows enter M (all for a
/// whole fragment). Special members out of line so that -Winline stays quiet.
struct PrecondJob {
    int frag = -1;
    bool segment = false;
    std::vector<int> real;
    Matrix caps;
    std::vector<char> own;
    PrecondJob();
    ~PrecondJob();
    PrecondJob(PrecondJob&&) noexcept;
    PrecondJob& operator=(PrecondJob&&) = delete;
};
PrecondJob::PrecondJob() = default;
PrecondJob::~PrecondJob() = default;
PrecondJob::PrecondJob(PrecondJob&&) noexcept = default;
} // namespace

/**
 * @brief Build M = sum_i R_i^T D_i (-J_i^-1) R_i + alpha on the rest (Claude Generated, Oct 2026).
 *
 * Fragments: detectCovalentFragments. A fragment gets its own block unless it is a single atom,
 * has an odd number of valence electrons as a neutral species, or carries a net charge in the
 * single-shot EEQ model (|sum q_EEQ| > 0.5 e) - then its indices keep the scalar alpha
 * (scf_damping). Every included fragment runs its own neutral closed-shell GFN1/GFN2 SCF (FP64,
 * default settings, the parent's electronic temperature, gas phase) and then the central-
 * difference Jacobian of its SCC map at its converged state (m_scf_jac_capture, h = 1e-4).
 *
 * A fragment with more than scf_precond_max_fragment_atoms atoms keeps the scalar seed, unless
 * scf_precond_segment_atoms > 0: then it is cut into segments (segmentFragment) and each segment
 * contributes the rows of its own atoms of -J_loc^-1, J_loc the Jacobian of its hydrogen-capped,
 * neutral local cluster (restricted additive Schwarz, D_i = 1 on own components, 0 elsewhere).
 * A cluster with an odd electron count (valence electrons + one per cap) is skipped; its own atoms
 * keep alpha. If the domain is the whole fragment the block is exact (no caps).
 *
 * All jobs (whole fragments and segments) run in one pool over the -threads budget (std::thread,
 * one XTB per job, each serial); every result goes to its own slot, so the outcome does not depend
 * on the thread count or the scheduling.
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
    const int nmp = (m_method == MethodType::GFN2) ? 9 : 0;   // GFN2 mixes dp (3) and qp (6) per atom, GFN1 only q_sh
    const int dim = nsh + nmp * nat;
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
    const bool segmenting = (m_scf_precond_segment_atoms > 0);
    std::vector<Skip> skip(frags.size(), Skip::None);
    std::vector<PrecondJob> jobs;
    int n_seg_frags = 0, n_segments = 0, n_seg_odd = 0, n_seg_odd_atoms = 0;
    for (size_t f = 0; f < frags.size(); ++f) {
        const auto& a = frags[f];
        double nel = 0.0, q = 0.0;
        for (int i : a) { nel += m_wfn.n0_at(i); if (have_eeq) q += q_eeq(i); }
        const bool large = m_scf_precond_max_frag_atoms > 0 && a.size() > static_cast<size_t>(m_scf_precond_max_frag_atoms);
        if (a.size() == 1) skip[f] = Skip::SingleAtom;
        else if (large && !segmenting) skip[f] = Skip::TooLarge;
        else if ((std::llround(nel) % 2) != 0) skip[f] = Skip::OddElectrons;
        else if (!have_eeq || std::fabs(q) > 0.5) skip[f] = Skip::Charged;
        else if (!large) {
            PrecondJob j;
            j.frag = static_cast<int>(f);
            j.real = a;
            j.own.assign(a.size(), 1);
            jobs.push_back(std::move(j));
        } else {
            ++n_seg_frags;
            for (ScfPrecondSegment& s : segmentFragment(m_atoms, m_geometry, a, m_scf_precond_segment_atoms,
                                                        m_scf_precond_segment_radius)) {
                ++n_segments;
                double nel_c = static_cast<double>(s.caps.rows());   // one electron per H cap
                for (int i : s.domain) nel_c += m_wfn.n0_at(i);
                if ((std::llround(nel_c) % 2) != 0) {
                    ++n_seg_odd;
                    n_seg_odd_atoms += static_cast<int>(s.own.size());
                    continue;
                }
                PrecondJob j;
                j.frag = static_cast<int>(f);
                j.segment = true;
                j.own.assign(s.domain.size(), 0);
                for (size_t k = 0, o = 0; k < s.domain.size(); ++k) {   // both lists ascending
                    while (o < s.own.size() && s.own[o] < s.domain[k]) ++o;
                    if (o < s.own.size() && s.own[o] == s.domain[k]) j.own[k] = 1;
                }
                j.real = std::move(s.domain);
                j.caps = std::move(s.caps);
                jobs.push_back(std::move(j));
            }
        }
    }

    // Elements outside the measured set (H, C, N, O): allowed, one warning per process.
    {
        std::set<int> other;
        for (const PrecondJob& j : jobs)
            for (int i : j.real)
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

    // Per-job Jacobian -> block B = -J^-1, restricted to the real atoms, rows of non-own atoms zero.
    struct Result { bool ok = false; std::string why; PrecondMatrix B; };
    std::vector<Result> res(jobs.size());
    const double kT_parent = m_electronic_temp;
    // Cluster Jacobian: central differences (fd) or linear response (Oct 2026); unknown value -> fd, one warning.
    std::string jac_mode = m_scf_precond_jacobian;
    if (jac_mode != "fd" && jac_mode != "response") {
        static std::atomic<bool> warned_jac{ false };
        if (!warned_jac.exchange(true))
            CurcumaLogger::warn("scf_precond_jacobian '" + jac_mode + "' unknown (fd | response); fd is used");
        jac_mode = "fd";
    }
    auto run_job = [&](size_t jj) {
        Result& r = res[jj];
        const PrecondJob& job = jobs[jj];
        const int nreal = static_cast<int>(job.real.size());
        const int ncap = static_cast<int>(job.caps.rows());
        const int nsub = nreal + ncap;
        try {
            Mol sub;
            sub.m_number_atoms = nsub;
            sub.m_charge = 0;
            sub.m_spin = 0.0;
            sub.m_geometry = Matrix(nsub, 3);
            for (int k = 0; k < nreal; ++k) {
                sub.m_atoms.push_back(m_atoms[job.real[k]]);
                sub.m_geometry.row(k) = m_geometry.row(job.real[k]);
            }
            for (int c = 0; c < ncap; ++c) {   // caps last: their rows/columns are dropped below
                sub.m_atoms.push_back(1);
                sub.m_geometry.row(nreal + c) = job.caps.row(c);
            }
            XTB frag(m_method);
            frag.setScfGuess("eeq");
            frag.setIntraThreads(1);
            frag.setElectronicTemperature(kT_parent);
            frag.setD4AtmCutoff(m_d4_atm_cutoff);
            frag.setSto6gLegacy4sp(m_sto6g_legacy_4sp);
            frag.setScfPrecondJacobian(jac_mode);
            frag.m_scf_jac_capture = true;
            if (!frag.InitialiseMolecule(sub)) { r.why = "fragment initialisation failed"; return; }
            frag.Calculation(false);
            if (frag.hasError() || !frag.m_scf_converged) { r.why = "fragment SCF failed"; return; }
            const int fnsh = frag.m_basis.nsh;
            const int fdim = fnsh + nmp * nsub;
            if (frag.m_scf_jac_result.rows() != fdim) { r.why = "fragment Jacobian missing"; return; }
            for (int k = 0; k < nreal; ++k)
                if (frag.m_basis.nsh_at[k] != m_basis.nsh_at[job.real[k]]) { r.why = "shell count mismatch"; return; }
            Eigen::PartialPivLU<PrecondMatrix> lu(frag.m_scf_jac_result);
            const PrecondMatrix Bloc = -lu.inverse();
            if (!Bloc.allFinite()) { r.why = "singular fragment Jacobian"; return; }
            if (ncap == 0 && !job.segment) { r.B = Bloc; r.ok = true; return; }   // whole fragment: local = block order
            // Local SCC positions of the real atoms in block order (shells atom by atom, then dp
            // and qp component-wise over the real atoms) and their own flags.
            std::vector<int> sel;
            std::vector<char> own_row;
            for (int k = 0; k < nreal; ++k)
                for (int s = 0; s < frag.m_basis.nsh_at[k]; ++s) { sel.push_back(frag.m_basis.ish_at[k] + s); own_row.push_back(job.own[k]); }
            for (int c = 0; c < nmp; ++c)   // GFN2: dp k = 0..2, qp k = 0..5 at offset nsh + c*nsub
                for (int k = 0; k < nreal; ++k) { sel.push_back(fnsh + c * nsub + k); own_row.push_back(job.own[k]); }
            const int n = static_cast<int>(sel.size());
            r.B = PrecondMatrix::Zero(n, n);
            for (int i = 0; i < n; ++i)
                if (own_row[i])
                    for (int k = 0; k < n; ++k) r.B(i, k) = Bloc(sel[i], sel[k]);
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

    // Assemble the Schwarz blocks (own rows disjoint, D_i = 1 on them).
    auto P = std::make_shared<SchwarzPreconditioner>(dim, alpha);
    int n_block = 0, n_failed = 0, n_scalar_atoms = n_seg_odd_atoms;
    long long seg_cluster_atoms = 0;
    int n_seg_jobs = 0;
    std::string first_fail;
    for (size_t j = 0; j < jobs.size(); ++j) {
        const PrecondJob& job = jobs[j];
        if (job.segment) {
            ++n_seg_jobs;
            seg_cluster_atoms += static_cast<long long>(job.real.size()) + job.caps.rows();
        }
        int n_own = 0;
        for (char o : job.own) n_own += o ? 1 : 0;
        if (!res[j].ok) {
            ++n_failed; n_scalar_atoms += n_own;
            if (first_fail.empty()) first_fail = res[j].why;
            continue;
        }
        SchwarzPreconditioner::Block b;
        std::vector<char> covers;
        // Order of the local SCC vector: shells atom by atom, then (GFN2) dp (k = x,y,z) and
        // qp (k = xx,xy,yy,xz,yz,zz) component-wise over the real atoms.
        for (size_t k = 0; k < job.real.size(); ++k) {
            const int a = job.real[k];
            for (int s = 0; s < m_basis.nsh_at[a]; ++s) { b.index.push_back(m_basis.ish_at[a] + s); covers.push_back(job.own[k]); }
        }
        for (int k = 0; k < nmp; ++k)   // GFN2 only: dp (k < 3) at nsh + k*nat, qp at nsh + 3*nat + (k-3)*nat
            for (size_t m = 0; m < job.real.size(); ++m) { b.index.push_back(nsh + k * nat + job.real[m]); covers.push_back(job.own[m]); }
        b.B = std::move(res[j].B);
        if (!P->addBlock(std::move(b), covers)) {
            err = "inconsistent block for fragment " + std::to_string(job.frag);
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
                                            "(a fragment Jacobian costs 2*dim SCC-map evaluations, ~2 min for 212 atoms on one thread; raise scf_precond_max_fragment_atoms, 0 = no limit, "
                                            "or segment them with scf_precond_segment_atoms)",
                                            n_skip[4], m_scf_precond_max_frag_atoms));
    }
    if (n_failed > 0)
        CurcumaLogger::warn(fmt::format("scf_precond frag: {} of {} block Jacobians failed (first: {}); "
                                        "their atoms keep the scalar seed", n_failed, jobs.size(), first_fail));

    m_scf_precond_frag = P;
    m_scf_precond_geom = m_geometry;
    m_scf_precond_atoms = m_atoms;
    publishSharedPreconditioner();
    const double t_ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();
    if (saved_verbosity >= 2) {
        CurcumaLogger::info(fmt::format("SCF precond: frag, {} fragments, dim {}, Jacobian {}, built in {:.1f} ms",
                                        frags.size(), dim, jac_mode, t_ms));
        CurcumaLogger::info(fmt::format("SCF precond: {} blocks, scalar seed on {} atoms (single atoms {}, "
                                        "odd electrons {}, charged {}, too large {}, failed {}), memory {:.3f} MB",
                                        n_block, n_scalar_atoms, n_skip[1], n_skip[2], n_skip[3], n_skip[4], n_failed,
                                        8.0 * static_cast<double>(P->storedDoubles()) / 1048576.0));
        if (n_seg_frags > 0)
            CurcumaLogger::info(fmt::format("SCF precond: {} fragment(s) segmented (RAS, {} heavy atoms per segment, R = {:.2f} A): "
                                            "{} segments, mean cluster {:.1f} atoms incl. caps, {} skipped (odd electrons, {} atoms on the scalar seed)",
                                            n_seg_frags, m_scf_precond_segment_atoms, m_scf_precond_segment_radius, n_segments,
                                            n_seg_jobs > 0 ? static_cast<double>(seg_cluster_atoms) / n_seg_jobs : 0.0, n_seg_odd, n_seg_odd_atoms));
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
    ensureHostIntegrals(__func__);   // Claude Generated (Oct 2026): host copies of S/H0/L/gamma are lazy
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
    if (gfn2 && m_mp_initialized) ensureHostMultipoleMatrices();   // developer dump: needs the host copies
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
    // Linear-response Jacobian (Oct 2026) for comparison with J.bin; closed shell only.
    double t_resp_ms = -1.0;
    bool have_resp = false;
    if (!openShell()) {
        const auto t_r0 = clk::now();
        PrecondMatrix Jr;
        std::string rerr;
        have_resp = computeScfJacobianResponse(x_star, Jr, rerr);
        t_resp_ms = std::chrono::duration<double, std::milli>(clk::now() - t_r0).count();
        if (have_resp) writeD("J_response.bin", Jr.data(), static_cast<size_t>(Jr.size()));   // row-major like J.bin
        else CurcumaLogger::warn("CURCUMA_SCF_JAC: " + rerr);
    }

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
    meta["jacobian_response_time_ms"] = t_resp_ms;
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
        { "J_response.bin", "float64 [dim, dim], linear-response Jacobian (computeScfJacobianResponse), same convention as J.bin; closed shell only" },
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
