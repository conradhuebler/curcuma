/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Homogeneous long-range tail of the periodic GFN-FF D4 dispersion (WP-PERIODIC-NONBONDED A3c) -
 * Claude Generated (Nov 2026).
 *
 * E_tail = -(2 pi / V) [ S6/(3 r_c^3) + 2 S8/(5 r_c^5) ], S6 = sum_ij zeta_i zeta_j C6_ij,
 * S8 = sum_ij zeta_i zeta_j r4r2ij C6_ij over all ordered pairs incl. i = j
 * (D4ParameterGenerator::homogeneousDispersionTail, derivation there).
 *
 * (a) Formula: (a1) prefactor against a quadrature of (N^2/2)(4 pi/V) int r^2 e(r) dr for one species
 *     with constant C6 (large-r form, 1e-12 relative); the damped form is reported against the error
 *     bound (1/3)(R0/r_c)^6, (5/13)(R0/r_c)^8. (a2) factorised sums on a water box against the explicit
 *     O(N^2) double sum over getChargeWeightedC6 (1e-12 relative).
 * (b) Central finite difference (Angstrom, h = 2e-4 / 1e-4) of the total dispersion (pair + tail) and of
 *     the tail alone for the atoms nearest a face, against G(disp on) - G(disp off) and
 *     G(tail on) - G(tail off); error quadratic in h.
 * (c) No cell: tail 0, energy bitwise equal with dispersion_tail false. Cell with dispersion_tail false:
 *     tail 0, pair term bitwise unchanged, E_tot(on) - E_tot(off) == E_tail.
 * (d) r_c dependence of E_pair + E_tail (r_c = 8, 10, 12 Bohr) on a dense random water packing
 *     (L = 28 Bohr, NOT an equilibrated liquid: random sequential placement, no MD).
 * (e) 2x1x1 supercell: E_tail and E_tail per atom equal to the cell (1e-12 relative).
 *
 * Measurement mode (`test_pbc_dispersion_tail --measure`, not run by ctest): 2x2x2 replica against the
 * small cell for the random packing and for a simple cubic water lattice, and the tail of a 12^3-water
 * lattice at r_c = L/2 (comparison with the A3a extrapolation). Water geometry: water-dimer.roo2.8 of
 * the structure library (registry "H2O_dimer"), moved rigidly only. Exit 0 when every check passes.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}
static std::string sci(double v)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", v);
    return b;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
static const double EH2KCAL = 627.509474;
static const double PI = 3.14159265358979323846;

static Mol dimerAngstrom()
{
    return TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
}

/// one water (O at the origin), Angstrom
static Eigen::MatrixXd waterTemplate(std::vector<int>& atoms)
{
    const Mol d = dimerAngstrom();
    Eigen::MatrixXd w(3, 3);
    atoms.clear();
    for (int i = 0; i < 3; ++i) {
        w.row(i) = d.m_geometry.row(i) - d.m_geometry.row(0);
        atoms.push_back(d.m_atoms[i]);
    }
    return w;
}

/// Random sequential placement of rigid waters in a cube (minimum image), random orientation.
/// min_oo: O-O distance, min_any: any intermolecular atom pair (Angstrom). Not a liquid.
static Mol randomWaterBox(int nwat, double L_ang, unsigned seed, double min_oo, double min_any)
{
    std::vector<int> wa;
    const Eigen::MatrixXd w = waterTemplate(wa);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Mol m = dimerAngstrom();
    m.m_atoms.clear();
    m.m_geometry = Geometry(3 * nwat, 3);
    int placed = 0;
    long tries = 0;
    auto mic = [&](Eigen::Vector3d v) {
        for (int k = 0; k < 3; ++k)
            v[k] -= L_ang * std::rint(v[k] / L_ang);
        return v;
    };
    while (placed < nwat && tries < 2000000) {
        ++tries;
        Eigen::Vector3d p(u(rng) * L_ang, u(rng) * L_ang, u(rng) * L_ang);
        Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
        qr.normalize();
        const Eigen::Matrix3d R = qr.toRotationMatrix();
        Eigen::Vector3d x[3];
        for (int i = 0; i < 3; ++i)
            x[i] = p + R * w.row(i).transpose();
        bool ok = true;
        for (int q = 0; q < placed && ok; ++q) {
            const Eigen::Vector3d oq = m.m_geometry.row(3 * q).transpose();
            if (mic(p - oq).norm() < min_oo) {
                ok = false;
                break;
            }
            for (int a = 0; a < 3 && ok; ++a)
                for (int b = 0; b < 3; ++b)
                    if (mic(x[a] - m.m_geometry.row(3 * q + b).transpose()).norm() < min_any) {
                        ok = false;
                        break;
                    }
        }
        if (!ok)
            continue;
        for (int i = 0; i < 3; ++i) {
            m.m_geometry.row(3 * placed + i) = x[i].transpose();
            m.m_atoms.push_back(wa[i]);
        }
        ++placed;
    }
    m.m_number_atoms = 3 * placed;
    m.m_geometry.conservativeResize(3 * placed, 3);
    check(placed == nwat, "random box placed all " + std::to_string(nwat) + " molecules (" + std::to_string(placed) + ")");
    return m;
}

/// Simple cubic lattice of waters (n^3 sites, spacing a Angstrom), random orientation. A lattice, not a liquid.
static Mol latticeWaterBox(int n, double a_ang, unsigned seed)
{
    std::vector<int> wa;
    const Eigen::MatrixXd w = waterTemplate(wa);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Mol m = dimerAngstrom();
    m.m_atoms.clear();
    m.m_geometry = Geometry(3 * n * n * n, 3);
    int k = 0;
    for (int ix = 0; ix < n; ++ix)
        for (int iy = 0; iy < n; ++iy)
            for (int iz = 0; iz < n; ++iz) {
                const Eigen::Vector3d p(a_ang * (ix + 0.5), a_ang * (iy + 0.5), a_ang * (iz + 0.5));
                Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
                qr.normalize();
                const Eigen::Matrix3d R = qr.toRotationMatrix();
                for (int i = 0; i < 3; ++i) {
                    m.m_geometry.row(k) = (p + R * w.row(i).transpose()).transpose();
                    m.m_atoms.push_back(wa[i]);
                    ++k;
                }
            }
    m.m_number_atoms = k;
    return m;
}

/// replica nx x ny x nz of a cubic cell (L Angstrom)
static Mol replicate(const Mol& box, double L, int nx, int ny, int nz)
{
    Mol sc = box;
    const int n = box.m_number_atoms;
    const int nr = nx * ny * nz;
    sc.m_geometry = Geometry(nr * n, 3);
    sc.m_atoms.clear();
    int r = 0;
    for (int a = 0; a < nx; ++a)
        for (int b = 0; b < ny; ++b)
            for (int c = 0; c < nz; ++c, ++r)
                for (int i = 0; i < n; ++i) {
                    sc.m_geometry.row(r * n + i) = box.m_geometry.row(i);
                    sc.m_geometry(r * n + i, 0) += a * L;
                    sc.m_geometry(r * n + i, 1) += b * L;
                    sc.m_geometry(r * n + i, 2) += c * L;
                    sc.m_atoms.push_back(box.m_atoms[i]);
                }
    sc.m_number_atoms = nr * n;
    sc.m_has_pbc = true;
    Eigen::Matrix3d cell = Eigen::Matrix3d::Zero();
    cell(0, 0) = nx * L;
    cell(1, 1) = ny * L;
    cell(2, 2) = nz * L;
    sc.m_unit_cell = cell;
    return sc;
}

// ---------------------------------------------------------------------------------------------
struct Run {
    bool ok = false;
    double e_disp = 0.0, e_tail = 0.0, e_tot = 0.0;
    double s6 = 0.0, s8 = 0.0, rc = 0.0, volume = 0.0;
    double s6_explicit = 0.0, s8_explicit = 0.0;  // only with explicit = true
    double e_tail_damped = 0.0;                    // explicit, with the R0 damping kept (diagnostic)
    Matrix grad;  // Eh/Bohr
    int natoms = 0;
};

static json boxFlags(double rc)
{
    json e = json::object();
    e["hbond"] = false;      // HB/XB search radii 15.8/21.2 Bohr exceed L/2 of the small cells
    e["repulsion"] = false;  // repulsion list 22 Bohr likewise
    if (rc > 0.0)
        e["dispersion_cutoff_bohr"] = rc;
    return e;
}

static double radialQuad(double rc, double c6, double zeta2, double r4r2, double r0, bool damped);

/// mol: cell taken from mol.m_unit_cell if mol.m_has_pbc, else none
static Run runGFNFF(const Mol& mol, const json& extra, bool gradient, int threads = 1, bool explicit_sums = false)
{
    Run r;
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    try {
        GFNFF ff(c);
        if (!ff.InitialiseMolecule(mol))
            return r;
        r.e_tot = ff.Calculation(gradient);
        r.e_disp = ff.DispersionEnergy();
        r.e_tail = ff.DispersionTailEnergy();
        if (gradient)
            r.grad = ff.Gradient();
        r.natoms = mol.m_number_atoms;
        r.rc = ff.dispersionRadiiBohr().eval;
        if (mol.m_has_pbc) {
            r.volume = mol.m_unit_cell(0, 0) * mol.m_unit_cell(1, 1) * mol.m_unit_cell(2, 2) * A2B * A2B * A2B;
            D4ParameterGenerator* gen = ff.getD4Generator();
            if (gen)
                gen->homogeneousDispersionTail(r.rc, r.volume, nullptr, &r.s6, &r.s8);
            if (gen && explicit_sums) {
                const auto& zeta = gen->atomZeta();
                const int n = mol.m_number_atoms;
                double s6 = 0.0, s8 = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) {
                        const int Zi = mol.m_atoms[i], Zj = mol.m_atoms[j];
                        const double c6 = gen->getChargeWeightedC6(Zi, Zj, i, j);
                        const double r4r2 = 3.0 * gen->getSqrtZr4r2(Zi) * gen->getSqrtZr4r2(Zj);
                        s6 += zeta[i] * zeta[j] * c6;
                        s8 += zeta[i] * zeta[j] * r4r2 * c6;
                    }
                r.s6_explicit = s6;
                r.s8_explicit = s8;
                // diagnostic: the same g(r) = 1 tail with the damping kept, per element pair by quadrature
                // (R0 = 0.58 sqrt(r4r2) + 4.80, the GFN-FF defaults of the pair list)
                std::map<std::pair<int, int>, double> jtab;
                double sd = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) {
                        const int Zi = mol.m_atoms[i], Zj = mol.m_atoms[j];
                        const double r4r2 = 3.0 * gen->getSqrtZr4r2(Zi) * gen->getSqrtZr4r2(Zj);
                        auto key = std::make_pair(std::min(Zi, Zj), std::max(Zi, Zj));
                        auto it = jtab.find(key);
                        if (it == jtab.end()) {
                            const double r0 = 0.58 * std::sqrt(r4r2) + 4.80;
                            // radialQuad(..., zeta2 = 1, c6 = 1) = -(J6 + 2 r4r2 J8)
                            it = jtab.emplace(key, -radialQuad(r.rc, 1.0, 1.0, r4r2, r0, true)).first;
                        }
                        sd += zeta[i] * zeta[j] * gen->getChargeWeightedC6(Zi, Zj, i, j) * it->second;
                    }
                r.e_tail_damped = -(2.0 * PI / r.volume) * sd;
            }
        }
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("GFNFF threw: %s\n", e.what());
    }
    return r;
}

static Mol withCell(Mol m, double L_ang)
{
    m.m_has_pbc = true;
    m.m_unit_cell = Eigen::Matrix3d::Identity() * L_ang;
    return m;
}

static const double kBoxBohr = 24.0;  // the A3b box: 21 waters, min O-O 3 A

static Mol a3bBox()
{
    return randomWaterBox(21, kBoxBohr / A2B, 20261103, 3.0, 0.0);
}

// ---------------------------------------------------------------------------------------------
// (a1) prefactor against a radial quadrature, one species with constant C6
static double radialQuad(double rc, double c6, double zeta2, double r4r2, double r0, bool damped)
{
    // int_{rc}^inf r^2 e(r) dr with r = rc/u, dr = rc/u^2 du, u in (0, 1]; composite Simpson, 20000 intervals
    const int M = 20000;
    double sum = 0.0;
    for (int k = 0; k <= M; ++k) {
        const double u = static_cast<double>(k) / M;
        double f = 0.0;
        if (u > 0.0) {
            const double r = rc / u;
            const double r2 = r * r, r6 = r2 * r2 * r2, r8 = r6 * r2;
            const double r06 = std::pow(r0, 6), r08 = std::pow(r0, 8);
            const double e = damped ? -zeta2 * c6 * (1.0 / (r6 + r06) + 2.0 * r4r2 / (r8 + r08))
                                    : -zeta2 * c6 * (1.0 / r6 + 2.0 * r4r2 / r8);
            f = r2 * e * rc / (u * u);
        }
        const double wgt = (k == 0 || k == M) ? 1.0 : (k % 2 ? 4.0 : 2.0);
        sum += wgt * f;
    }
    return sum / (3.0 * M);
}

static void testPrefactor()
{
    const int N = 50;
    const double c6 = 10.0, zeta = 1.1, r4r2 = 20.2, V = kBoxBohr * kBoxBohr * kBoxBohr;
    const double r0 = 0.58 * std::sqrt(r4r2) + 4.80;
    for (double rc : { 11.5, 35.0 }) {
        const double s6 = double(N) * N * zeta * zeta * c6;
        const double s8 = double(N) * N * zeta * zeta * r4r2 * c6;
        const double e_formula = D4ParameterGenerator::homogeneousTailFromSums(s6, s8, rc, V);
        // E_tail = (1/2) sum_i sum_j (1/V) 4 pi int r^2 e dr = (N^2/2)(4 pi / V) int
        const double e_quad = 0.5 * N * N * (4.0 * PI / V) * radialQuad(rc, c6, zeta * zeta, r4r2, r0, false);
        const double rel = std::abs(e_formula - e_quad) / std::abs(e_quad);
        check(rel < 1e-12, "prefactor at r_c " + std::to_string(rc) + ": formula == quadrature (" + sci(rel) + ")");
        const double e_damped = 0.5 * N * N * (4.0 * PI / V) * radialQuad(rc, c6, zeta * zeta, r4r2, r0, true);
        const double rel_d = (e_formula - e_damped) / e_damped;  // > 0: large-r form over-binds
        const double bound = std::max(std::pow(r0 / rc, 6) / 3.0, 5.0 / 13.0 * std::pow(r0 / rc, 8));
        check(rel_d > 0.0 && rel_d <= bound, "damping error at r_c " + std::to_string(rc) + " within bound (" + sci(rel_d) + " <= " + sci(bound) + ")");
        std::printf("(a1) r_c %.1f: formula %.15e, quadrature rel %.1e; damped form differs by %.3e (bound %.3e, R0 %.3f Bohr)\n",
            rc, e_formula, rel, rel_d, bound, r0);
    }
}

// (a2) factorisation against the explicit double sum
static void testFactorisation()
{
    const Mol box = withCell(a3bBox(), kBoxBohr / A2B);
    const Run r = runGFNFF(box, boxFlags(11.5), false, 1, true);
    check(r.ok, "(a2) box runs");
    if (!r.ok)
        return;
    const double d6 = std::abs(r.s6 - r.s6_explicit) / std::abs(r.s6_explicit);
    const double d8 = std::abs(r.s8 - r.s8_explicit) / std::abs(r.s8_explicit);
    const double e_ex = D4ParameterGenerator::homogeneousTailFromSums(r.s6_explicit, r.s8_explicit, r.rc, r.volume);
    const double de = std::abs(r.e_tail - e_ex) / std::abs(e_ex);
    check(d6 < 1e-12 && d8 < 1e-12, "(a2) factorised S6/S8 == explicit double sum (" + sci(d6) + ", " + sci(d8) + ")");
    check(de < 1e-12, "(a2) E_tail == explicit (" + sci(de) + ")");
    check(r.e_tail < 0.0, "(a2) E_tail negative");
    const double n = r.natoms;
    std::printf("(a2) %d atoms, r_c %.2f: S6 %.10e (explicit rel %.1e), S8 %.10e (rel %.1e), <zeta^2 C6> %.4f au, "
                "E_tail %.12e Eh (rel %.1e), E_pair %.12e, C8 share of tail %.4f\n",
        r.natoms, r.rc, r.s6, d6, r.s8, d8, r.s6 / (n * n), r.e_tail, de, r.e_disp,
        D4ParameterGenerator::homogeneousTailFromSums(0.0, r.s8, r.rc, r.volume) / r.e_tail);
}

// (b) finite difference
static void testFiniteDifference()
{
    const double L = kBoxBohr / A2B;
    const Mol box = withCell(a3bBox(), L);
    const json on_flags = boxFlags(0.0);  // default cut-off reduced to L/2 (12 Bohr)
    json tail_off = on_flags;
    tail_off["dispersion_tail"] = false;
    json disp_off = on_flags;
    disp_off["dispersion"] = false;
    const Run on = runGFNFF(box, on_flags, true);
    const Run toff = runGFNFF(box, tail_off, true);
    const Run doff = runGFNFF(box, disp_off, true);
    check(on.ok && toff.ok && doff.ok, "(b) reference runs");
    if (!on.ok || !toff.ok || !doff.ok)
        return;
    check(doff.e_tail == 0.0 && toff.e_tail == 0.0, "(b) no tail with dispersion or dispersion_tail off");
    const Matrix g_disp = on.grad - doff.grad;  // pair + tail, Eh/Bohr
    const Matrix g_tail = on.grad - toff.grad;  // tail only
    std::vector<std::pair<double, int>> dist;
    for (int i = 0; i < box.m_number_atoms; ++i) {
        double dmin = 1e300;
        for (int k = 0; k < 3; ++k) {
            const double x = box.m_geometry(i, k) - L * std::floor(box.m_geometry(i, k) / L);
            dmin = std::min({ dmin, x, L - x });
        }
        dist.push_back({ dmin, i });
    }
    std::sort(dist.begin(), dist.end());
    const int natoms_fd = 6;
    double max_g = 0.0, max_t = 0.0;
    auto fdError = [&](double h, double& err_tail) {
        double err = 0.0;
        err_tail = 0.0;
        for (int s = 0; s < natoms_fd; ++s) {
            const int a = dist[s].second;
            for (int k = 0; k < 3; ++k) {
                Mol p = box, m = box;
                p.m_geometry(a, k) += h;
                m.m_geometry(a, k) -= h;
                const Run rp = runGFNFF(p, on_flags, false), rm = runGFNFF(m, on_flags, false);
                const double fd = ((rp.e_disp + rp.e_tail) - (rm.e_disp + rm.e_tail)) / (2.0 * h);  // Eh/A
                const double fd_t = (rp.e_tail - rm.e_tail) / (2.0 * h);
                err = std::max(err, std::abs(fd - g_disp(a, k) * A2B));
                err_tail = std::max(err_tail, std::abs(fd_t - g_tail(a, k) * A2B));
                max_g = std::max(max_g, std::abs(g_disp(a, k) * A2B));
                max_t = std::max(max_t, std::abs(g_tail(a, k) * A2B));
            }
        }
        return err;
    };
    double t2 = 0.0, t1 = 0.0;
    const double e2 = fdError(2e-4, t2), e1 = fdError(1e-4, t1);
    check(e1 < 1e-8, "(b) pair + tail gradient vs central FD, max error " + sci(e1) + " Eh/A");
    check(e1 < e2 / 3.0, "(b) FD error quadratic in h (" + sci(e2) + " -> " + sci(e1) + ")");
    check(t1 < 1e-3 * max_t + 1e-13, "(b) tail gradient vs FD of E_tail (" + sci(t1) + " of max " + sci(max_t) + ")");
    check(max_t > 1e-9, "(b) tail force non-zero (test not vacuous): " + sci(max_t));
    double tail_all = g_tail.cwiseAbs().maxCoeff() * A2B, disp_all = g_disp.cwiseAbs().maxCoeff() * A2B;
    std::printf("(b) E_pair %.12e, E_tail %.12e Eh; 6 atoms %.2f..%.2f A from a face: max |g_disp| %.3e, max |g_tail| %.3e Eh/A; "
                "|an - FD| total %.2e (h 2e-4) %.2e (h 1e-4), tail %.2e / %.2e; all atoms max |g_tail| %.3e of max |g_disp| %.3e\n",
        on.e_disp, on.e_tail, dist[0].first, dist[natoms_fd - 1].first, max_g, max_t, e2, e1, t2, t1, tail_all, disp_all);
}

// (c) tail vanishes without a cell and when switched off
static void testSwitch()
{
    const Mol plain = a3bBox();
    json off = boxFlags(0.0);
    off["dispersion_tail"] = false;
    const Run n_on = runGFNFF(plain, boxFlags(0.0), true);
    const Run n_off = runGFNFF(plain, off, true);
    check(n_on.ok && n_off.ok, "(c) no-cell runs");
    check(n_on.e_tail == 0.0, "(c) no cell: E_tail == 0");
    check(n_on.e_tot == n_off.e_tot && n_on.e_disp == n_off.e_disp, "(c) no cell: energy bitwise equal with dispersion_tail false");
    check(n_on.grad.rows() == n_off.grad.rows() && (n_on.grad - n_off.grad).cwiseAbs().maxCoeff() == 0.0,
        "(c) no cell: gradient bitwise equal with dispersion_tail false");
    const Mol box = withCell(plain, kBoxBohr / A2B);
    const Run c_on = runGFNFF(box, boxFlags(0.0), true);
    const Run c_off = runGFNFF(box, off, true);
    check(c_on.ok && c_off.ok, "(c) cell runs");
    check(c_off.e_tail == 0.0 && c_on.e_tail < 0.0, "(c) cell: tail only when switched on");
    check(c_on.e_disp == c_off.e_disp, "(c) cell: pair term bitwise unchanged by the tail");
    const double d = std::abs((c_on.e_tot - c_off.e_tot) - c_on.e_tail);
    check(d < 1e-13, "(c) cell: E_tot(on) - E_tot(off) == E_tail (" + sci(d) + ")");
    std::printf("(c) no cell E_tail %.1e, |E(on)-E(off)| %.1e; cell: E_tail %.12e, |dE_tot - E_tail| %.1e\n",
        n_on.e_tail, std::abs(n_on.e_tot - n_off.e_tot), c_on.e_tail, d);
}

// (d) r_c dependence on a dense random packing
static const double kDenseBohr = 28.0;
static Mol denseBox()
{
    // liquid water 0.0334 / A^3 -> 108.4 molecules in (28 Bohr)^3; 100 placed (0.92 of liquid density)
    return withCell(randomWaterBox(100, kDenseBohr / A2B, 4711, 2.6, 1.7), kDenseBohr / A2B);
}

static void testCutoffDependence()
{
    const Mol box = denseBox();
    double e_pair[3], e_sum[3], e_tl[3], e_dmp[3];
    const double rcs[3] = { 8.0, 10.0, 12.0 };
    for (int k = 0; k < 3; ++k) {
        const Run r = runGFNFF(box, boxFlags(rcs[k]), false, 1, true);
        check(r.ok, "(d) runs at r_c " + std::to_string(rcs[k]));
        e_pair[k] = r.e_disp;
        e_tl[k] = r.e_tail;
        e_sum[k] = r.e_disp + r.e_tail;
        e_dmp[k] = r.e_disp + r.e_tail_damped;
    }
    const double dp = e_pair[2] - e_pair[0], ds = e_sum[2] - e_sum[0];
    std::printf("(d) dense random packing %d atoms, L %.1f Bohr (not a liquid):\n", box.m_number_atoms, kDenseBohr);
    for (int k = 0; k < 3; ++k)
        std::printf("    r_c %4.1f: E_pair %.10f  E_tail %.10f  sum %.10f Eh (sum - sum(12) %+.4f kcal/mol; "
                    "with damped tail (diagnostic) %+.4f kcal/mol)\n",
            rcs[k], e_pair[k], e_tl[k], e_sum[k], (e_sum[k] - e_sum[2]) * EH2KCAL, (e_dmp[k] - e_dmp[2]) * EH2KCAL);
    std::printf("    8 -> 12 Bohr: E_pair moves %.4f kcal/mol, E_pair + E_tail %.4f kcal/mol\n", dp * EH2KCAL, ds * EH2KCAL);
    check(std::abs(ds) < 0.25 * std::abs(dp), "(d) E_pair + E_tail depends less on r_c than E_pair (" + sci(ds) + " vs " + sci(dp) + ")");
}

// (e) 2x1x1 supercell
static void testSupercell()
{
    const double L = kBoxBohr / A2B;
    const Mol cell = withCell(a3bBox(), L);
    const Mol sc = replicate(cell, L, 2, 1, 1);
    const Run r1 = runGFNFF(cell, boxFlags(11.5), false);
    const Run r2 = runGFNFF(sc, boxFlags(11.5), false);
    check(r1.ok && r2.ok, "(e) supercell runs");
    const double rel = std::abs(r2.e_tail - 2.0 * r1.e_tail) / std::abs(r1.e_tail);
    const double pa1 = r1.e_tail / r1.natoms, pa2 = r2.e_tail / r2.natoms;
    check(rel < 1e-12, "(e) 2x1x1: E_tail == 2 x cell (" + sci(rel) + ")");
    std::printf("(e) 2x1x1: E_tail %.15e vs 2 x cell %.15e (rel %.1e); per atom %.6e / %.6e Eh; pair term %.15e vs %.15e\n",
        r2.e_tail, 2.0 * r1.e_tail, rel, pa1, pa2, r2.e_disp, 2.0 * r1.e_disp);
}

// ---------------------------------------------------------------------------------------------
// measurement mode
static void replicaCompare(const char* name, const Mol& small, double L_ang)
{
    const auto t0 = std::chrono::steady_clock::now();
    const Run rs = runGFNFF(small, boxFlags(0.0), false, 8);
    const Mol rep = replicate(small, L_ang, 2, 2, 2);
    const Run rr = runGFNFF(rep, boxFlags(0.0), false, 8);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!rs.ok || !rr.ok) {
        std::printf("%s: run failed\n", name);
        return;
    }
    const double es = rs.e_disp + rs.e_tail, er = (rr.e_disp + rr.e_tail) / 8.0;
    std::printf("%s: %d atoms, L %.3f Bohr | small r_c %.3f: pair %.10f tail %.10f | replica %d atoms r_c %.3f (per cell): pair %.10f tail %.10f\n",
        name, rs.natoms, L_ang * A2B, rs.rc, rs.e_disp, rs.e_tail, rr.natoms, rr.rc, rr.e_disp / 8.0, rr.e_tail / 8.0);
    std::printf("    shell L/2..L: explicit (replica pair - small pair) %+.6f kcal/mol, small-cell tail - replica tail %+.6f kcal/mol\n",
        (rr.e_disp / 8.0 - rs.e_disp) * EH2KCAL, (rs.e_tail - rr.e_tail / 8.0) * EH2KCAL);
    std::printf("    (pair+tail) small - replica %+.6f kcal/mol (%+.3e kcal/mol per atom), small-cell tail %.6f kcal/mol; %.0f s\n",
        (es - er) * EH2KCAL, (es - er) * EH2KCAL / rs.natoms, rs.e_tail * EH2KCAL, secs);
}

static void measure()
{
    const Mol dense = denseBox();
    replicaCompare("random packing", dense, kDenseBohr / A2B);
    const double a = 3.1;  // Angstrom: 0.0336 molecules / A^3, liquid density
    replicaCompare("cubic lattice 5^3", withCell(latticeWaterBox(5, a, 99), 5 * a), 5 * a);
    // 12^3 = 1728 waters, L = 37.2 A, r_c = L/2: tail only needs C6 and the cell
    {
        const Mol big = withCell(latticeWaterBox(12, a, 7), 12 * a);
        const Run r = runGFNFF(big, boxFlags(0.0), false, 8);
        if (r.ok) {
            const double n = r.natoms;
            std::printf("cubic lattice 12^3 (%d atoms, L %.2f Bohr, r_c %.3f): E_tail %.6f kcal/mol, per atom %.3e kcal/mol, "
                        "<zeta^2 C6> = S6/N^2 %.4f au, C8 share %.4f, E_pair %.6f Eh, tail/pair %.4e\n",
                r.natoms, 12 * a * A2B, r.rc, r.e_tail * EH2KCAL, r.e_tail * EH2KCAL / n, r.s6 / (n * n),
                D4ParameterGenerator::homogeneousTailFromSums(0.0, r.s8, r.rc, r.volume) / r.e_tail, r.e_disp, r.e_tail / r.e_disp);
        } else
            std::printf("cubic lattice 12^3: run failed\n");
    }
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc > 1 && std::strcmp(argv[1], "--measure") == 0) {
        measure();
        return g_failures == 0 ? 0 : 1;
    }
    testPrefactor();
    testFactorisation();
    testFiniteDifference();
    testSwitch();
    testCutoffDependence();
    testSupercell();
    if (g_failures == 0)
        std::printf("pbc_dispersion_tail: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
