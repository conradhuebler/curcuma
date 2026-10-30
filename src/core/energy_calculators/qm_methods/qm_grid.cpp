/*
 * <Native QM: molecular DFT quadrature grid and AO values on grid points (WP4)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas and literature: see qm_grid.h. Claude Generated (Oct 2026, WP4).
 *
 * This program is free software under GPL-3.0
 */

#include "qm_grid.h"
#include "qm_grid_tables.hpp"
#include "qm_xc_kernels.h"  // partition polynomials and AO components shared with cuda/qm_gpu_xc.cu

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace qmgrid {

GridOptions::GridOptions() = default;
GridOptions::GridOptions(const GridOptions&) = default;
GridOptions& GridOptions::operator=(const GridOptions&) = default;
GridOptions::~GridOptions() = default;
MolecularGrid::MolecularGrid() = default;
MolecularGrid::MolecularGrid(const MolecularGrid&) = default;
MolecularGrid& MolecularGrid::operator=(const MolecularGrid&) = default;
MolecularGrid::~MolecularGrid() = default;
AOBlock::AOBlock() = default;
AOBlock::~AOBlock() = default;
AOEvaluator::Shell::Shell() = default;
AOEvaluator::Shell::Shell(const Shell&) = default;
AOEvaluator::Shell::Shell(Shell&&) noexcept = default;
AOEvaluator::Shell& AOEvaluator::Shell::operator=(const Shell&) = default;
AOEvaluator::Shell& AOEvaluator::Shell::operator=(Shell&&) noexcept = default;
AOEvaluator::Shell::~Shell() = default;
FlatShells::FlatShells() = default;
FlatShells::~FlatShells() = default;
AOEvaluator::AOEvaluator() = default;
AOEvaluator::AOEvaluator(const AOEvaluator&) = default;
AOEvaluator& AOEvaluator::operator=(const AOEvaluator&) = default;
AOEvaluator::~AOEvaluator() = default;

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Octahedral orbit expansion (Laikov's gen_oh) in the point order of pyscf.dft.LebedevGrid.SphGenOh.
void sphGenOh(int code, double a, double b, double v, std::vector<double>& out)
{
    auto add = [&](double x, double y, double z) {
        out.push_back(x);
        out.push_back(y);
        out.push_back(z);
        out.push_back(v);
    };
    switch (code) {
    case 0:  // (1,0,0) and permutations: 6 points
        add(1, 0, 0); add(-1, 0, 0); add(0, 1, 0); add(0, -1, 0); add(0, 0, 1); add(0, 0, -1);
        break;
    case 1: {  // (0,a,a), a = 1/sqrt(2): 12 points
        a = std::sqrt(0.5);
        add(0, a, a); add(0, -a, a); add(0, a, -a); add(0, -a, -a);
        add(a, 0, a); add(-a, 0, a); add(a, 0, -a); add(-a, 0, -a);
        add(a, a, 0); add(-a, a, 0); add(a, -a, 0); add(-a, -a, 0);
    } break;
    case 2: {  // (a,a,a), a = 1/sqrt(3): 8 points
        a = std::sqrt(1.0 / 3.0);
        add(a, a, a); add(-a, a, a); add(a, -a, a); add(-a, -a, a);
        add(a, a, -a); add(-a, a, -a); add(a, -a, -a); add(-a, -a, -a);
    } break;
    case 3: {  // (a,a,b), b = sqrt(1 - 2 a^2): 24 points
        b = std::sqrt(1.0 - 2.0 * a * a);
        add(a, a, b); add(-a, a, b); add(a, -a, b); add(-a, -a, b);
        add(a, a, -b); add(-a, a, -b); add(a, -a, -b); add(-a, -a, -b);
        add(a, b, a); add(-a, b, a); add(a, -b, a); add(-a, -b, a);
        add(a, b, -a); add(-a, b, -a); add(a, -b, -a); add(-a, -b, -a);
        add(b, a, a); add(-b, a, a); add(b, -a, a); add(-b, -a, a);
        add(b, a, -a); add(-b, a, -a); add(b, -a, -a); add(-b, -a, -a);
    } break;
    case 4: {  // (a,b,0), b = sqrt(1 - a^2): 24 points
        b = std::sqrt(1.0 - a * a);
        add(a, b, 0); add(-a, b, 0); add(a, -b, 0); add(-a, -b, 0);
        add(b, a, 0); add(-b, a, 0); add(b, -a, 0); add(-b, -a, 0);
        add(a, 0, b); add(-a, 0, b); add(a, 0, -b); add(-a, 0, -b);
        add(b, 0, a); add(-b, 0, a); add(b, 0, -a); add(-b, 0, -a);
        add(0, a, b); add(0, -a, b); add(0, a, -b); add(0, -a, -b);
        add(0, b, a); add(0, -b, a); add(0, b, -a); add(0, -b, -a);
    } break;
    case 5: {  // (a,b,c), c = sqrt(1 - a^2 - b^2): 48 points
        const double c = std::sqrt(1.0 - a * a - b * b);
        const double t[6][3] = { { a, b, c }, { a, c, b }, { b, a, c }, { b, c, a }, { c, a, b }, { c, b, a } };
        for (const auto& p : t) {
            add(p[0], p[1], p[2]); add(-p[0], p[1], p[2]); add(p[0], -p[1], p[2]); add(-p[0], -p[1], p[2]);
            add(p[0], p[1], -p[2]); add(-p[0], p[1], -p[2]); add(p[0], -p[1], -p[2]); add(-p[0], -p[1], -p[2]);
        }
    } break;
    default:
        break;
    }
}

const tables::LebedevGridDef* findLebedev(int npts)
{
    for (int k = 0; k < tables::kNumLebedevGrids; ++k)
        if (tables::kLebedevGrids[k].npts == npts) return &tables::kLebedevGrids[k];
    return nullptr;
}

std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

// Partition polynomials (Becke, Stratmann) and their derivatives: qm_xc_kernels.h (shared with the GPU).
using kern::beckeIterate;
using kern::stratmannIterate;
using kern::beckeIterateDerivative;
using kern::stratmannIterateDerivative;

/// 1/R_ij and the size-adjustment parameters a_ij of the partition (shared by buildMolecularGrid and
/// gridWeightGradient). Claude Generated (WP4; split out for WP8).
bool partitionPairs(const std::vector<int>& Z, const Matrix& posBohr, const std::string& adjust, Matrix& inv_dist,
                    Matrix& a, std::string& why)
{
    const int natoms = static_cast<int>(Z.size());
    inv_dist = Matrix::Zero(natoms, natoms);
    a = Matrix::Zero(natoms, natoms);
    for (int i = 0; i < natoms; ++i)
        for (int j = 0; j < natoms; ++j) {
            if (i == j) continue;
            const double d = (posBohr.row(i) - posBohr.row(j)).norm();
            if (d < 1.0e-8) {
                why = "molecular grid: atoms " + std::to_string(i + 1) + " and " + std::to_string(j + 1) + " coincide";
                return false;
            }
            inv_dist(i, j) = 1.0 / d;
            if (adjust != "none") {
                double ri = tables::kBraggRadiiBohr[Z[i]], rj = tables::kBraggRadiiBohr[Z[j]];
                if (adjust == "treutler") {
                    ri = std::sqrt(ri);
                    rj = std::sqrt(rj);
                }
                // pyscf: rr[i,j] = rad_i / rad_j, a = 1/4 (rr^T - rr) -> a[i,j] = 1/4 (rad_j/rad_i - rad_i/rad_j)
                double aij = 0.25 * (rj / ri - ri / rj);
                aij = std::min(0.5, std::max(-0.5, aij));
                a(i, j) = aij;
            }
        }
    return true;
}

}  // namespace

// =====================================================================================================
// Level names, radial and angular grids
// =====================================================================================================

bool parseGridLevel(const std::string& level_in, int& n_radial, int& n_angular, std::string& why)
{
    const std::string level = lower(level_in);
    // Point counts of xcDFT read_grid.f90 (SG-0 ... SG-3 sizes), used here WITHOUT pruning.
    if (level == "sg0") { n_radial = 23; n_angular = 170; }
    else if (level == "sg1") { n_radial = 50; n_angular = 194; }
    else if (level == "sg2") { n_radial = 75; n_angular = 302; }
    else if (level == "sg3") { n_radial = 99; n_angular = 590; }
    else {
        const auto x = level.find('x');
        if (x == std::string::npos || x == 0 || x + 1 >= level.size()) {
            why = "unknown grid '" + level_in + "' (use sg0|sg1|sg2|sg3 or <n_radial>x<n_angular>, e.g. 75x302)";
            return false;
        }
        try {
            size_t p1 = 0, p2 = 0;
            n_radial = std::stoi(level.substr(0, x), &p1);
            n_angular = std::stoi(level.substr(x + 1), &p2);
            if (p1 != x || p2 != level.size() - x - 1) throw 0;
        } catch (...) {
            why = "cannot read grid '" + level_in + "' as <n_radial>x<n_angular>";
            return false;
        }
    }
    if (n_radial < 2 || n_radial > 1000) {
        why = "grid '" + level_in + "': number of radial points must be 2..1000";
        return false;
    }
    if (!findLebedev(n_angular)) {
        why = "grid '" + level_in + "': " + std::to_string(n_angular) + " is not a Lebedev grid size (6, 14, 26, 38, "
              "50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, ... 5810)";
        return false;
    }
    return true;
}

bool radialGrid(const std::string& scheme_in, int n, int Z, std::vector<double>& r, std::vector<double>& w,
                std::string& why)
{
    if (Z < 1 || Z > tables::kMaxTableZ) {
        why = "radial grid: element Z = " + std::to_string(Z) + " outside the tables (1.." +
              std::to_string(tables::kMaxTableZ) + ")";
        return false;
    }
    const std::string scheme = lower(scheme_in);
    r.assign(n, 0.0);
    w.assign(n, 0.0);
    if (scheme == "treutler") {
        // Treutler-Ahlrichs M4 (JCP 102, 346 (1995), eq. 19 with alpha = 0.6), Chebyshev 2nd kind nodes
        // x_i = cos(i pi/(n+1)); r = xi/ln2 (1+x)^0.6 ln(2/(1-x)). dr carries the Chebyshev weight
        // pi/(n+1) sin(theta) (the sin^2 of the 2nd-kind rule divided by sqrt(1-x^2)). Same arithmetic as
        // pyscf.dft.radi.treutler_ahlrichs; stored with ascending r (pyscf reverses the same way).
        const double xi = tables::kTreutlerXi[Z];
        const double step = kPi / (n + 1);
        const double ln2 = xi / std::log(2.0);
        for (int i = 0; i < n; ++i) {
            const double x = std::cos((i + 1) * step);
            const double ri = -ln2 * std::pow(1.0 + x, 0.6) * std::log((1.0 - x) / 2.0);
            const double dr = step * std::sin((i + 1) * step) * ln2 * std::pow(1.0 + x, 0.6) *
                              (-0.6 / (1.0 + x) * std::log((1.0 - x) / 2.0) + 1.0 / (1.0 - x));
            r[n - 1 - i] = ri;
            w[n - 1 - i] = 4.0 * kPi * ri * ri * dr;
        }
        return true;
    }
    if (scheme == "euler_maclaurin") {
        // Murray, Handy, Laming, Mol. Phys. 78, 997 (1993), m_R = 2 (xcDFT dft_grid.f EulMac):
        // r_i = R (i/(n+1-i))^2, r^2 dr = 2 R^3 (n+1) i^5 / (n+1-i)^7, R = Bragg radius (Bohr).
        const double R = tables::kBraggRadiiBohr[Z];
        for (int i = 1; i <= n; ++i) {
            const double q = static_cast<double>(n + 1 - i);
            const double ri = R * (i / q) * (i / q);
            const double r2dr = 2.0 * R * R * R * (n + 1) * std::pow(static_cast<double>(i), 5) / std::pow(q, 7);
            r[i - 1] = ri;
            w[i - 1] = 4.0 * kPi * r2dr;
        }
        return true;
    }
    why = "unknown radial grid '" + scheme_in + "' (use treutler|euler_maclaurin)";
    return false;
}

bool lebedevGrid(int npts, Matrix& xyz, Vector& w)
{
    const tables::LebedevGridDef* def = findLebedev(npts);
    if (!def) return false;
    std::vector<double> buf;
    buf.reserve(4 * npts);
    for (int k = 0; k < def->count; ++k) {
        const tables::LebedevOrbit& o = tables::kLebedevOrbits[def->first + k];
        sphGenOh(o.code, o.a, o.b, o.v, buf);
    }
    if (static_cast<int>(buf.size()) != 4 * npts) return false;  // table inconsistent
    xyz.resize(npts, 3);
    w.resize(npts);
    for (int p = 0; p < npts; ++p) {
        xyz(p, 0) = buf[4 * p];
        xyz(p, 1) = buf[4 * p + 1];
        xyz(p, 2) = buf[4 * p + 2];
        w(p) = buf[4 * p + 3];
    }
    return true;
}

int lebedevDegree(int npts)
{
    const tables::LebedevGridDef* def = findLebedev(npts);
    return def ? def->degree : 0;
}

// =====================================================================================================
// Molecular grid with Becke partition
// =====================================================================================================

bool buildMolecularGrid(const std::vector<int>& Z, const Matrix& posBohr, const GridOptions& opt,
                        MolecularGrid& grid, std::string& why, int threads)
{
    const int natoms = static_cast<int>(Z.size());
    if (natoms == 0 || posBohr.rows() != natoms || posBohr.cols() != 3) {
        why = "molecular grid: no atoms or position matrix does not match";
        return false;
    }
    int nr = 0, na = 0;
    if (!parseGridLevel(opt.level, nr, na, why)) return false;
    const std::string partition = lower(opt.partition);
    const std::string adjust = lower(opt.radii_adjust);
    if (partition != "becke" && partition != "stratmann") {
        why = "unknown grid partition '" + opt.partition + "' (use becke|stratmann)";
        return false;
    }
    if (adjust != "treutler" && adjust != "becke" && adjust != "none") {
        why = "unknown grid radii adjustment '" + opt.radii_adjust + "' (use treutler|becke|none)";
        return false;
    }
    for (int z : Z)
        if (z < 1 || z > tables::kMaxTableZ) {
            why = "molecular grid: element Z = " + std::to_string(z) + " outside the grid tables";
            return false;
        }

    Matrix ang;
    Vector angw;
    lebedevGrid(na, ang, angw);

    // Atomic grids (one per atom; identical for equal elements).
    const int per_atom = nr * na;
    const int npts = natoms * per_atom;
    grid.n_radial = nr;
    grid.n_angular = na;
    grid.points.resize(npts, 3);
    grid.weights.resize(npts);
    grid.quadrature_weights.resize(npts);
    grid.atom.assign(npts, 0);
    grid.atom_offset.assign(natoms + 1, 0);
    for (int A = 0; A < natoms; ++A) {
        std::vector<double> r, wr;
        if (!radialGrid(opt.radial, nr, Z[A], r, wr, why)) return false;
        const int off = A * per_atom;
        grid.atom_offset[A] = off;
        for (int i = 0; i < nr; ++i)
            for (int j = 0; j < na; ++j) {
                const int p = off + i * na + j;
                for (int k = 0; k < 3; ++k) grid.points(p, k) = ang(j, k) * r[i] + posBohr(A, k);
                grid.quadrature_weights(p) = wr[i] * angw(j);
                grid.atom[p] = A;
            }
    }
    grid.atom_offset[natoms] = npts;

    // Interatomic distances and size-adjustment parameters a_ij (pyscf radi.*_atomic_radii_adjust).
    Matrix inv_dist, a;
    if (!partitionPairs(Z, posBohr, adjust, inv_dist, a, why)) return false;
    const bool stratmann = (partition == "stratmann");

    // Becke weights: P_A(r) = prod_{B != A} s(nu_AB), w = W_quad * P_owner / sum_B P_B. Loop order (i, j < i) as
    // pyscf get_partition, so the rounding is the same.
#pragma omp parallel for schedule(static) num_threads(threads > 0 ? threads : 1)
    for (int p = 0; p < npts; ++p) {
        std::vector<double> dist(natoms), P(natoms, 1.0);
        for (int A = 0; A < natoms; ++A) {
            const double dx = grid.points(p, 0) - posBohr(A, 0);
            const double dy = grid.points(p, 1) - posBohr(A, 1);
            const double dz = grid.points(p, 2) - posBohr(A, 2);
            dist[A] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        for (int i = 0; i < natoms; ++i)
            for (int j = 0; j < i; ++j) {
                double g = inv_dist(i, j) * (dist[i] - dist[j]);  // mu_ij
                g = g + a(i, j) * (1.0 - g * g);                  // nu_ij (size adjustment)
                g = stratmann ? stratmannIterate(g) : beckeIterate(g);
                P[i] *= 0.5 * (1.0 - g);
                P[j] *= 0.5 * (1.0 + g);
            }
        double sum = 0.0;
        for (int A = 0; A < natoms; ++A) sum += P[A];
        grid.weights(p) = grid.quadrature_weights(p) * P[grid.atom[p]] * (1.0 / sum);
    }
    return true;
}

// =====================================================================================================
// AO values on grid points
// =====================================================================================================

namespace {
void powersOf(GTO::OrbitalType t, int& l, int& m, int& n)
{
    l = m = n = 0;
    switch (t) {
    case GTO::PX: l = 1; break;
    case GTO::PY: m = 1; break;
    case GTO::PZ: n = 1; break;
    case GTO::DXX: l = 2; break;
    case GTO::DYY: m = 2; break;
    case GTO::DZZ: n = 2; break;
    case GTO::DXY: l = 1; m = 1; break;
    case GTO::DXZ: l = 1; n = 1; break;
    case GTO::DYZ: m = 1; n = 1; break;
    default: break;
    }
}

}  // namespace

AOEvaluator::AOEvaluator(const std::vector<GTO::Orbital>& basis, const Matrix& Q)
{
    m_ncart = static_cast<int>(basis.size());
    for (int mu = 0; mu < m_ncart; ++mu) {
        const GTO::Orbital& o = basis[mu];
        int l, m, n;
        powersOf(o.type, l, m, n);
        bool join = false;
        if (!m_shells.empty()) {
            const Shell& s = m_shells.back();
            join = (s.x == o.x && s.y == o.y && s.z == o.z && s.exps == o.exponents &&
                    s.l[0] + s.m[0] + s.n[0] == l + m + n);
        }
        if (!join) {
            Shell s;
            s.x = o.x;
            s.y = o.y;
            s.z = o.z;
            s.exps = o.exponents;
            s.amin = s.exps.empty() ? 0.0 : *std::min_element(s.exps.begin(), s.exps.end());
            m_shells.push_back(std::move(s));
        }
        Shell& s = m_shells.back();
        s.cart.push_back(mu);
        s.l.push_back(l);
        s.m.push_back(m);
        s.n.push_back(n);
        s.coef.insert(s.coef.end(), o.coefficients.begin(), o.coefficients.end());
    }
    m_transform = Q.size() > 0;
    m_nbf = m_transform ? static_cast<int>(Q.cols()) : m_ncart;
    if (m_transform) {
        m_q_cols.assign(m_nbf, {});
        for (int c = 0; c < Q.rows(); ++c)
            for (int a = 0; a < Q.cols(); ++a)
                if (Q(c, a) != 0.0) m_q_cols[a].emplace_back(c, Q(c, a));
    }
    // AO screening (Claude Generated, Oct 2026): active range per shell, as flatten() derives it
    const int nsh = static_cast<int>(m_shells.size());
    m_out_first.assign(nsh, 0);
    m_out_width.assign(nsh, 0);
    m_shell_local = true;
    if (!m_transform) {
        for (int s = 0; s < nsh; ++s) {
            m_out_first[s] = m_shells[s].cart.empty() ? 0 : m_shells[s].cart[0];
            m_out_width[s] = static_cast<int>(m_shells[s].cart.size());
        }
        return;
    }
    std::vector<int> cart_shell(m_ncart, -1);
    for (int s = 0; s < nsh; ++s)
        for (int c : m_shells[s].cart) cart_shell[c] = s;
    std::vector<int> col_shell(m_nbf, -1);
    for (int a = 0; a < m_nbf; ++a)
        for (const auto& [c, qv] : m_q_cols[a]) {
            if (col_shell[a] >= 0 && col_shell[a] != cart_shell[c]) m_shell_local = false;
            col_shell[a] = cart_shell[c];
        }
    std::vector<int> lo(nsh, m_nbf), hi(nsh, -1);
    for (int a = 0; a < m_nbf; ++a) {
        if (col_shell[a] < 0) {
            m_shell_local = false;
            continue;
        }
        lo[col_shell[a]] = std::min(lo[col_shell[a]], a);
        hi[col_shell[a]] = std::max(hi[col_shell[a]], a);
    }
    for (int s = 0; s < nsh; ++s) {
        if (hi[s] < lo[s]) {
            m_shell_local = false;
            continue;
        }
        m_out_first[s] = lo[s];
        m_out_width[s] = hi[s] - lo[s] + 1;
        for (int a = lo[s]; a <= hi[s]; ++a)
            if (col_shell[a] != s) m_shell_local = false;  // the columns of a shell must be contiguous
    }
}

void AOEvaluator::shellCentre(int s, double& x, double& y, double& z) const
{
    x = m_shells[s].x;
    y = m_shells[s].y;
    z = m_shells[s].z;
}

// Claude Generated (Oct 2026, AO screening): reach of every shell, formula in qm_grid.h.
std::vector<double> AOEvaluator::shellRadii(double threshold) const
{
    const int nsh = static_cast<int>(m_shells.size());
    std::vector<double> R(nsh, std::numeric_limits<double>::infinity());
    if (!(threshold > 0.0)) return R;
    for (int s = 0; s < nsh; ++s) {
        const Shell& sh = m_shells[s];
        const int np = static_cast<int>(sh.exps.size());
        if (np == 0 || sh.cart.empty() || !(sh.amin > 0.0)) continue;  // no exponent: never screened
        const int L = sh.l[0] + sh.m[0] + sh.n[0];
        std::vector<double> c(np, 0.0);
        for (size_t comp = 0; comp < sh.cart.size(); ++comp)
            for (int k = 0; k < np; ++k) c[k] = std::max(c[k], std::abs(sh.coef[comp * np + k]));
        double t = 1.0;
        if (m_transform) {
            t = 0.0;
            for (int w = 0; w < m_out_width[s]; ++w) {
                double sum = 0.0;
                for (const auto& [ci, qv] : m_q_cols[m_out_first[s] + w]) sum += std::abs(qv);
                t = std::max(t, sum);
            }
        }
        auto rpow = [](double r, int k) { return k < 0 ? 0.0 : std::pow(r, k); };
        auto bound = [&](double r) {
            double b0 = 0.0, b1 = 0.0, b2 = 0.0;
            for (int k = 0; k < np; ++k) {
                const double a = sh.exps[k], e = c[k] * std::exp(-a * r * r);
                b0 += rpow(r, L) * e;
                b1 += (L * rpow(r, L - 1) + 2.0 * a * rpow(r, L + 1)) * e;
                b2 += (L * (L - 1) * rpow(r, L - 2) + 2.0 * a * (2 * L + 1) * rpow(r, L) + 4.0 * a * a * rpow(r, L + 2)) * e;
            }
            return t * std::max(b0, std::max(b1, b2));
        };
        const double rtop = std::sqrt(80.0 / sh.amin);  // evaluate() skips the shell beyond this radius
        // inwards in coarse steps of rtop/200 to the first point with bound >= threshold, then in fine steps of
        // rtop/4000 inside the last coarse interval
        const int ncoarse = 200, nfine = 20;
        const double H = rtop / ncoarse, h = H / nfine;
        int i = ncoarse;
        while (i > 0 && bound(H * i) < threshold) --i;
        if (i == ncoarse) {  // the bound reaches the threshold at the evaluate() cutoff: the cutoff is the reach
            R[s] = rtop;
            continue;
        }
        int j = nfine;
        while (j > 0 && bound(H * i + h * j) < threshold) --j;
        R[s] = std::min(rtop, H * i + h * j + 2.0 * h);
    }
    return R;
}

// Claude Generated (Oct 2026, AO screening): evaluate() restricted to the listed shells, compact columns.
void AOEvaluator::evaluateShells(const Matrix& points, int begin, int count, int deriv, const int* shells, int nsh,
                                 AOBlock& out, std::vector<int>& cols) const
{
    const bool grad = deriv >= 1;
    const bool hess = deriv >= 2;
    std::vector<int> ccol(nsh), ocol(nsh);
    int nc = 0, no = 0;
    cols.clear();
    for (int j = 0; j < nsh; ++j) {
        const int s = shells[j];
        ccol[j] = nc;
        nc += static_cast<int>(m_shells[s].cart.size());
        ocol[j] = no;
        no += m_out_width[s];
        for (int w = 0; w < m_out_width[s]; ++w) cols.push_back(m_out_first[s] + w);
    }
    Matrix cbuf[10];
    Matrix* dst[10] = { &out.phi, &out.dx, &out.dy, &out.dz, &out.dxx, &out.dxy, &out.dxz, &out.dyy, &out.dyz, &out.dzz };
    Matrix* v[10];
    for (int k = 0; k < 10; ++k) v[k] = m_transform ? &cbuf[k] : dst[k];
    const int nmat = hess ? 10 : (grad ? 4 : 1);
    for (int k = 0; k < nmat; ++k) v[k]->setZero(count, nc);
    std::vector<double> e, ae, aae;
    for (int j = 0; j < nsh; ++j) {
        const Shell& s = m_shells[shells[j]];
        const int nprim = static_cast<int>(s.exps.size());
        e.resize(nprim);
        ae.resize(nprim);
        aae.resize(nprim);
        for (int q = 0; q < count; ++q) {
            const int p = begin + q;
            const double x = points(p, 0) - s.x, y = points(p, 1) - s.y, z = points(p, 2) - s.z;
            const double r2 = x * x + y * y + z * z;
            if (s.amin * r2 > 80.0) continue;
            for (int k = 0; k < nprim; ++k) {
                e[k] = std::exp(-s.exps[k] * r2);
                ae[k] = s.exps[k] * e[k];
                aae[k] = s.exps[k] * ae[k];
            }
            for (size_t c = 0; c < s.cart.size(); ++c) {
                const double* cc = &s.coef[c * nprim];
                double R0 = 0.0, R1 = 0.0, R2 = 0.0;
                for (int k = 0; k < nprim; ++k) {
                    R0 += cc[k] * e[k];
                    R1 += cc[k] * ae[k];
                }
                if (hess)
                    for (int k = 0; k < nprim; ++k) R2 += cc[k] * aae[k];
                double val[10];
                kern::aoComponent(s.l[c], s.m[c], s.n[c], x, y, z, R0, R1, R2, deriv >= 2 ? 2 : (grad ? 1 : 0), val);
                for (int k = 0; k < nmat; ++k) (*v[k])(q, ccol[j] + static_cast<int>(c)) = val[k];
            }
        }
    }
    for (int k = nmat; k < 10; ++k) dst[k]->resize(0, 0);
    if (!m_transform) return;
    for (int k = 0; k < nmat; ++k) {
        Matrix& d = *dst[k];
        const Matrix& src = cbuf[k];
        d.setZero(count, no);
        for (int j = 0; j < nsh; ++j) {
            const int s = shells[j], c0 = m_shells[s].cart[0];
            for (int w = 0; w < m_out_width[s]; ++w)
                for (const auto& [c, qv] : m_q_cols[m_out_first[s] + w]) {
                    const int lc = ccol[j] + (c - c0);
                    for (int q = 0; q < count; ++q) d(q, ocol[j] + w) += qv * src(q, lc);
                }
        }
    }
}

// =====================================================================================================
// AO screening (Claude Generated, Oct 2026): compact blocks and shell lists per block
// =====================================================================================================

BlockScreen::BlockScreen() = default;
BlockScreen::BlockScreen(const BlockScreen&) = default;
BlockScreen& BlockScreen::operator=(const BlockScreen&) = default;
BlockScreen::~BlockScreen() = default;

namespace {
/// Recursive median split of idx[lo, hi) along the longest bounding-box axis into ceil(n / block) leaves.
void splitPoints(std::vector<int>& idx, const Matrix& pts, int lo, int hi, int block, std::vector<int>& begin,
                 std::vector<int>& count)
{
    const int n = hi - lo;
    if (n <= block) {
        if (n > 0) {
            begin.push_back(lo);
            count.push_back(n);
        }
        return;
    }
    double mn[3] = { 1e300, 1e300, 1e300 }, mx[3] = { -1e300, -1e300, -1e300 };
    for (int i = lo; i < hi; ++i)
        for (int k = 0; k < 3; ++k) {
            mn[k] = std::min(mn[k], pts(idx[i], k));
            mx[k] = std::max(mx[k], pts(idx[i], k));
        }
    int axis = 0;
    for (int k = 1; k < 3; ++k)
        if (mx[k] - mn[k] > mx[axis] - mn[axis]) axis = k;
    const int k = (n + block - 1) / block, kl = k / 2;
    const int nl = static_cast<int>((long long)n * kl / k);
    std::nth_element(idx.begin() + lo, idx.begin() + lo + nl, idx.begin() + hi, [&](int a, int b) {
        const double va = pts(a, axis), vb = pts(b, axis);
        return va < vb || (va == vb && a < b);
    });
    splitPoints(idx, pts, lo, lo + nl, block, begin, count);
    splitPoints(idx, pts, lo + nl, hi, block, begin, count);
}
}  // namespace

void sortGridBlocks(MolecularGrid& grid, int block, std::vector<int>& begin, std::vector<int>& count)
{
    begin.clear();
    count.clear();
    if (block < 1) block = 128;
    const int npts = grid.size();
    std::vector<int> idx(npts);
    for (int p = 0; p < npts; ++p) idx[p] = p;
    const int natoms = static_cast<int>(grid.atom_offset.size()) - 1;
    for (int A = 0; A < natoms; ++A) splitPoints(idx, grid.points, grid.atom_offset[A], grid.atom_offset[A + 1], block, begin, count);
    const bool has_quad = grid.quadrature_weights.size() == grid.weights.size();
    Matrix pts(npts, 3);
    Vector w(npts), qw(has_quad ? npts : 0);
    std::vector<int> at(npts);
    for (int i = 0; i < npts; ++i) {
        const int p = idx[i];
        pts.row(i) = grid.points.row(p);
        w(i) = grid.weights(p);
        if (has_quad) qw(i) = grid.quadrature_weights(p);
        at[i] = grid.atom[p];
    }
    grid.points = std::move(pts);
    grid.weights = std::move(w);
    if (has_quad) grid.quadrature_weights = std::move(qw);
    grid.atom = std::move(at);
}

bool buildBlockScreen(const MolecularGrid& grid, const AOEvaluator& ao, double threshold, const std::vector<int>& begin,
                      const std::vector<int>& count, BlockScreen& out, std::string& why)
{
    out = BlockScreen();
    if (!(threshold > 0.0)) return true;
    if (!ao.shellLocal()) {
        why = "AO screening: the AO transform mixes shells";
        return false;
    }
    const std::vector<double> R = ao.shellRadii(threshold);
    const int nsh = ao.nshell();
    std::vector<double> cx(nsh), cy(nsh), cz(nsh);
    std::vector<int> width(nsh);
    for (int s = 0; s < nsh; ++s) {
        ao.shellCentre(s, cx[s], cy[s], cz[s]);
        width[s] = ao.shellFunctions(s);
    }
    out.threshold = threshold;
    out.begin = begin;
    out.count = count;
    out.shell_off.assign(begin.size() + 1, 0);
    for (size_t b = 0; b < begin.size(); ++b) {
        double mn[3] = { 1e300, 1e300, 1e300 }, mx[3] = { -1e300, -1e300, -1e300 };
        for (int p = begin[b]; p < begin[b] + count[b]; ++p)
            for (int k = 0; k < 3; ++k) {
                mn[k] = std::min(mn[k], grid.points(p, k));
                mx[k] = std::max(mx[k], grid.points(p, k));
            }
        for (int s = 0; s < nsh; ++s) {
            const double c[3] = { cx[s], cy[s], cz[s] };
            double d2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double d = std::max(0.0, std::max(mn[k] - c[k], c[k] - mx[k]));
                d2 += d * d;
            }
            if (d2 <= R[s] * R[s]) {
                out.shells.push_back(s);
                out.fn_kept += width[s];
            }
        }
        out.fn_total += ao.nbf();
        out.shell_off[b + 1] = static_cast<int>(out.shells.size());
    }
    return true;
}

void AOEvaluator::evaluate(const Matrix& points, int begin, int count, int deriv, AOBlock& out) const
{
    const bool grad = deriv >= 1;
    const bool hess = deriv >= 2;  // WP8: second derivatives for the GGA nuclear gradient
    // Cartesian values go straight into `out` without a transform, else into local buffers.
    Matrix cbuf[10];
    Matrix* dst[10] = { &out.phi, &out.dx, &out.dy, &out.dz, &out.dxx, &out.dxy, &out.dxz, &out.dyy, &out.dyz, &out.dzz };
    Matrix* v[10];
    for (int k = 0; k < 10; ++k) v[k] = m_transform ? &cbuf[k] : dst[k];
    const int nmat = hess ? 10 : (grad ? 4 : 1);
    for (int k = 0; k < nmat; ++k) v[k]->setZero(count, m_ncart);
    std::vector<double> e, ae, aae;
    for (const Shell& s : m_shells) {
        const int nprim = static_cast<int>(s.exps.size());
        e.resize(nprim);
        ae.resize(nprim);
        aae.resize(nprim);
        for (int q = 0; q < count; ++q) {
            const int p = begin + q;
            const double x = points(p, 0) - s.x, y = points(p, 1) - s.y, z = points(p, 2) - s.z;
            const double r2 = x * x + y * y + z * z;
            if (s.amin * r2 > 80.0) continue;
            for (int k = 0; k < nprim; ++k) {
                e[k] = std::exp(-s.exps[k] * r2);
                ae[k] = s.exps[k] * e[k];
                aae[k] = s.exps[k] * ae[k];
            }
            for (size_t c = 0; c < s.cart.size(); ++c) {
                const double* cc = &s.coef[c * nprim];
                double R0 = 0.0, R1 = 0.0, R2 = 0.0;  // sum c e,  sum c a e,  sum c a^2 e
                for (int k = 0; k < nprim; ++k) {
                    R0 += cc[k] * e[k];
                    R1 += cc[k] * ae[k];
                }
                if (hess)
                    for (int k = 0; k < nprim; ++k) R2 += cc[k] * aae[k];
                const int mu = s.cart[c];
                double val[10];
                kern::aoComponent(s.l[c], s.m[c], s.n[c], x, y, z, R0, R1, R2, deriv >= 2 ? 2 : (grad ? 1 : 0), val);
                for (int k = 0; k < nmat; ++k) (*v[k])(q, mu) = val[k];
            }
        }
    }
    for (int k = nmat; k < 10; ++k) dst[k]->resize(0, 0);
    if (!m_transform) return;
    // Active basis: phi_a = sum_c Q(c, a) phi_c (Q is block diagonal per shell; only its nonzeros are used).
    auto transform = [&](const Matrix& src, Matrix& d) {
        d.setZero(count, m_nbf);
        for (int a = 0; a < m_nbf; ++a)
            for (const auto& [c, qv] : m_q_cols[a])
                for (int q = 0; q < count; ++q) d(q, a) += qv * src(q, c);
    };
    for (int k = 0; k < nmat; ++k) transform(cbuf[k], *dst[k]);
}

// Claude Generated (Oct 2026, GPU KS-DFT): flat tables of the shells for cuda/qm_gpu_xc.cu.
bool AOEvaluator::flatten(FlatShells& f, std::string& why) const
{
    f = FlatShells();
    f.nbf = m_nbf;
    f.ncart = m_ncart;
    f.nshell = static_cast<int>(m_shells.size());
    // active column -> shell that owns it (all its nonzeros must lie in that shell's cartesian rows)
    std::vector<int> cart_shell(m_ncart, -1);
    for (int si = 0; si < f.nshell; ++si)
        for (int c : m_shells[si].cart) cart_shell[c] = si;
    for (int si = 0; si < f.nshell; ++si) {
        const Shell& s = m_shells[si];
        const int nc = static_cast<int>(s.cart.size());
        const int np = static_cast<int>(s.exps.size());
        f.xyz.insert(f.xyz.end(), { s.x, s.y, s.z });
        f.amin.push_back(s.amin);
        f.nprim.push_back(np);
        f.prim_off.push_back(static_cast<int>(f.exps.size()));
        f.exps.insert(f.exps.end(), s.exps.begin(), s.exps.end());
        f.ncomp.push_back(nc);
        f.comp_off.push_back(static_cast<int>(f.coef_off.size()));
        f.cart_first.push_back(s.cart.empty() ? 0 : s.cart[0]);
        f.max_ncomp = std::max(f.max_ncomp, nc);
        for (int c = 0; c < nc; ++c) {
            if (c > 0 && s.cart[c] != s.cart[c - 1] + 1) {
                why = "AO tables: cartesian components of a shell are not consecutive";
                return false;
            }
            f.lmn.insert(f.lmn.end(), { s.l[c], s.m[c], s.n[c] });
            f.coef_off.push_back(static_cast<int>(f.coef.size()));
            f.coef.insert(f.coef.end(), s.coef.begin() + (size_t)c * np, s.coef.begin() + (size_t)(c + 1) * np);
        }
        if (!m_transform) {
            f.out_first.push_back(s.cart.empty() ? 0 : s.cart[0]);
            f.out_width.push_back(nc);
            f.T_off.push_back(-1);
            continue;
        }
        int lo = m_nbf, hi = -1;
        for (int a = 0; a < m_nbf; ++a)
            for (const auto& [c, qv] : m_q_cols[a])
                if (cart_shell[c] == si) {
                    lo = std::min(lo, a);
                    hi = std::max(hi, a);
                }
        if (hi < lo) {
            why = "AO tables: shell without an active function";
            return false;
        }
        f.out_first.push_back(lo);
        f.out_width.push_back(hi - lo + 1);
        f.T_off.push_back(static_cast<int>(f.T.size()));
        std::vector<double> T((size_t)nc * (hi - lo + 1), 0.0);
        for (int a = lo; a <= hi; ++a)
            for (const auto& [c, qv] : m_q_cols[a]) {
                if (cart_shell[c] != si) {
                    why = "AO tables: the spherical transform mixes shells";
                    return false;
                }
                T[(size_t)(c - s.cart[0]) * (hi - lo + 1) + (a - lo)] = qv;
            }
        f.T.insert(f.T.end(), T.begin(), T.end());
    }
    return true;
}

bool partitionTables(const std::vector<int>& Z, const Matrix& posBohr, const GridOptions& opt, Matrix& inv_dist,
                     Matrix& a, bool& stratmann, std::string& why)
{
    const std::string partition = lower(opt.partition);
    stratmann = (partition == "stratmann");
    if (partition != "becke" && !stratmann) {
        why = "unknown grid partition '" + opt.partition + "'";
        return false;
    }
    return partitionPairs(Z, posBohr, lower(opt.radii_adjust), inv_dist, a, why);
}

void densityOnBlock(const AOBlock& ao, const Matrix& P, Vector& rho, Matrix* grad)
{
    const Matrix Pphi = ao.phi * P;  // (P phi)_pm, P symmetric
    const int nb = static_cast<int>(ao.phi.rows());
    rho.resize(nb);
    for (int p = 0; p < nb; ++p) rho(p) = ao.phi.row(p).dot(Pphi.row(p));
    if (grad) {
        grad->resize(nb, 3);
        for (int p = 0; p < nb; ++p) {
            (*grad)(p, 0) = 2.0 * ao.dx.row(p).dot(Pphi.row(p));
            (*grad)(p, 1) = 2.0 * ao.dy.row(p).dot(Pphi.row(p));
            (*grad)(p, 2) = 2.0 * ao.dz.row(p).dot(Pphi.row(p));
        }
    }
}

double integrateDensity(const MolecularGrid& grid, const AOEvaluator& ao, const Matrix& P, int threads, int block)
{
    const int npts = grid.size();
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    std::vector<double> partial(nblocks, 0.0);
#pragma omp parallel for schedule(dynamic) num_threads(threads > 0 ? threads : 1)
    for (int b = 0; b < nblocks; ++b) {
        const int begin = b * block;
        const int count = std::min(block, npts - begin);
        AOBlock blk;
        ao.evaluate(grid.points, begin, count, 0, blk);
        Vector rho;
        densityOnBlock(blk, P, rho);
        double s = 0.0;
        for (int q = 0; q < count; ++q) s += grid.weights(begin + q) * rho(q);
        partial[b] = s;
    }
    double sum = 0.0;
    for (double v : partial) sum += v;
    return sum;
}


// =====================================================================================================
// WP8: nuclear derivative of the partition weights (grid response)
// =====================================================================================================
//
// w_p = W_p P_B(r_p) / Z(r_p),  Z = sum_C P_C,  P_C = prod_{D != C} s_CD,  s_CD = 1/2 (1 - h(nu_CD)),
// nu_CD = mu_CD + a_CD (1 - mu_CD^2),  mu_CD = (|r - R_C| - |r - R_D|) / R_CD.
// For A != B (point fixed in space):  d mu_CD / d R_C = -u_C / R_CD - mu_CD (R_C - R_D) / R_CD^2,
//                                    d mu_CD / d R_D = +u_D / R_CD + mu_CD (R_C - R_D) / R_CD^2,
// u_C = (r - R_C) / |r - R_C|;  d s_CD = -1/2 h'(nu) (1 - 2 a_CD mu_CD) d mu_CD;
// d w / d R_A = W (d P_B / Z - P_B / Z^2 sum_C d P_C),  d w / d R_B = - sum_{A != B} d w / d R_A.
// Claude Generated (Oct 2026, WP8).
bool gridWeightGradient(const std::vector<int>& Z, const Matrix& posBohr, const GridOptions& opt,
                        const MolecularGrid& grid, const Vector& f, Matrix& grad, std::string& why, int threads)
{
    const int natoms = static_cast<int>(Z.size());
    const int npts = grid.size();
    if (f.size() != npts || posBohr.rows() != natoms) {
        why = "grid weight gradient: sizes do not match";
        return false;
    }
    const std::string partition = lower(opt.partition);
    const std::string adjust = lower(opt.radii_adjust);
    Matrix inv_dist, a;
    if (!partitionPairs(Z, posBohr, adjust, inv_dist, a, why)) return false;
    const bool stratmann = (partition == "stratmann");
    if (grad.rows() != natoms || grad.cols() != 3) grad = Matrix::Zero(natoms, 3);
    if (natoms < 2) return true;

    // fixed chunking (independent of the thread count), chunks added in order
    const int nchunks = std::max(1, std::min(64, npts / 256));
    std::vector<Matrix> gc(nchunks);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int p0 = static_cast<int>((long long)npts * c / nchunks);
        const int p1 = static_cast<int>((long long)npts * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        std::vector<double> dist(natoms), ux(natoms), uy(natoms), uz(natoms), P(natoms);
        Matrix S(natoms, natoms), dS(natoms, natoms), excl(natoms, natoms), MU(natoms, natoms);
        std::vector<double> pre(natoms + 1), suf(natoms + 1);
        std::vector<double> dP(3 * (size_t)natoms);
        for (int p = p0; p < p1; ++p) {
            if (f(p) == 0.0) continue;
            const int B = grid.atom[p];
            for (int A = 0; A < natoms; ++A) {
                const double dx = grid.points(p, 0) - posBohr(A, 0);
                const double dy = grid.points(p, 1) - posBohr(A, 1);
                const double dz = grid.points(p, 2) - posBohr(A, 2);
                dist[A] = std::sqrt(dx * dx + dy * dy + dz * dz);
                const double inv = dist[A] > 0.0 ? 1.0 / dist[A] : 0.0;
                ux[A] = dx * inv;
                uy[A] = dy * inv;
                uz[A] = dz * inv;
            }
            // cell functions s_ij and d s_ij / d mu_ij for every ordered pair (s_ji = 1 - s_ij, mu_ji = -mu_ij)
            for (int i = 0; i < natoms; ++i) {
                S(i, i) = 1.0;
                dS(i, i) = 0.0;
                for (int j = 0; j < i; ++j) {
                    const double mu = inv_dist(i, j) * (dist[i] - dist[j]);
                    const double nu = mu + a(i, j) * (1.0 - mu * mu);
                    const double h = stratmann ? stratmannIterate(nu) : beckeIterate(nu);
                    const double dh = (stratmann ? stratmannIterateDerivative(nu) : beckeIterateDerivative(nu))
                                    * (1.0 - 2.0 * a(i, j) * mu);
                    S(i, j) = 0.5 * (1.0 - h);
                    S(j, i) = 0.5 * (1.0 + h);
                    dS(i, j) = -0.5 * dh;  // d s_ij / d mu_ij
                    dS(j, i) = -0.5 * dh;  // d s_ji / d mu_ji (mu_ji = -mu_ij, s_ji = 1 - s_ij)
                    MU(i, j) = mu;
                    MU(j, i) = -mu;
                }
            }
            // P_i and P_i without the factor s_ij (prefix/suffix products, no division)
            double Zs = 0.0;
            for (int i = 0; i < natoms; ++i) {
                pre[0] = 1.0;
                for (int j = 0; j < natoms; ++j) pre[j + 1] = pre[j] * (j == i ? 1.0 : S(i, j));
                suf[natoms] = 1.0;
                for (int j = natoms - 1; j >= 0; --j) suf[j] = suf[j + 1] * (j == i ? 1.0 : S(i, j));
                P[i] = pre[natoms];
                for (int j = 0; j < natoms; ++j) excl(i, j) = (j == i) ? 0.0 : pre[j] * suf[j + 1];
                Zs += P[i];
            }
            if (!(Zs > 0.0)) continue;
            const double fw = f(p) * grid.quadrature_weights(p);
            const double invZ = 1.0 / Zs, PB = P[B];
            for (int A = 0; A < natoms; ++A) {
                if (A == B) continue;
                // d P_C / d R_A: C = A (all pairs (A, D)) and C != A (the pair (C, A))
                std::fill(dP.begin(), dP.end(), 0.0);
                for (int D = 0; D < natoms; ++D) {
                    if (D == A) continue;
                    // d mu_AD / d R_A = -u_A / R_AD - mu_AD (R_A - R_D) / R_AD^2
                    const double iR = inv_dist(A, D), mu = MU(A, D);
                    double dmu[3];
                    const double uA[3] = { ux[A], uy[A], uz[A] };
                    for (int k = 0; k < 3; ++k)
                        dmu[k] = -uA[k] * iR - mu * (posBohr(A, k) - posBohr(D, k)) * iR * iR;
                    // factor s_AD in P_A, factor s_DA in P_D (d mu_DA = -d mu_AD)
                    const double cA = excl(A, D) * dS(A, D), cD = -excl(D, A) * dS(D, A);
                    for (int k = 0; k < 3; ++k) {
                        dP[3 * A + k] += cA * dmu[k];
                        dP[3 * D + k] += cD * dmu[k];
                    }
                }
                for (int k = 0; k < 3; ++k) {
                    double sum = 0.0;
                    for (int C = 0; C < natoms; ++C) sum += dP[3 * C + k];
                    const double dw = fw * (dP[3 * B + k] * invZ - PB * invZ * invZ * sum);
                    g(A, k) += dw;
                    g(B, k) -= dw;
                }
            }
        }
        gc[c] = std::move(g);
    }
    for (int c = 0; c < nchunks; ++c) grad += gc[c];
    return true;
}

}  // namespace qmgrid

