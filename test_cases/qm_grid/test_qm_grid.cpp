/*
 * Native QM DFT grid test driver (WP4).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout; check_qm_grid.py compares with reference.json from scripts/qm_grid_pyscf_reference.py):
 *
 *   test_qm_grid lebedev
 *       For every tabulated Lebedev grid: max |sum_j w_j Y_lm(u_j) - sqrt(4 pi) delta_l0| over all real spherical
 *       harmonics with l <= degree, the same for l = degree + 1 (must NOT be exact: shows the degree is not higher),
 *       and for grids up to 590 points the orthonormality error max |<Y_lm|Y_l'm'> - delta| for l, l' <= degree/2.
 *
 *   test_qm_grid grid <xyz> [--grid sg1] [--radial treutler] [--partition becke] [--adjust treutler]
 *                [--stride K] [--dump-all FILE]
 *       Grid totals (points, sum of weights, integral of sum_A exp(-|r - R_A|^2) and of a p-type test function),
 *       per-atom weight sums and every K-th point (index, atom, i_rad, j_ang, x, y, z, w, w_quad) for a pointwise
 *       comparison.
 *
 *   test_qm_grid density <xyz> [--grid sg1] [--threads N] [--stride K]
 *       HF/def2-SVP SCF (scf_threshold 1e-10), integral of rho on the grid (electron count) and rho, grad rho at
 *       every K-th grid point.
 *
 *   test_qm_grid nelec <xyz> [--levels sg0,sg1,sg2,sg3] [--threads N]
 *       Electron count deviation per grid level from one HF density (the WP4 table).
 *
 *   test_qm_grid aofd <xyz> [--cartesian]
 *       AO gradients vs 5-point finite differences of the AO values (h = 1e-4 Bohr) at every 97th grid point,
 *       error relative to max(1, |grad phi|).
 *
 * Claude Generated (Oct 2026, WP4). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/energy_calculators/qm_methods/qm_grid.h"
#include "src/core/units.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr double kPi = 3.14159265358979323846;

struct Args {
    std::string mode, xyz, grid = "sg1", radial = "treutler", partition = "becke", adjust = "treutler";
    std::string levels = "sg0,sg1,sg2,sg3", dump_all;
    int stride = 0, threads = 4;
    bool cartesian = false;
};

bool readXYZ(const std::string& path, std::vector<int>& Z, std::vector<double>& ang)
{
    static const std::map<std::string, int> el{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
        { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 } };
    std::ifstream f(path);
    int n = 0;
    if (!(f >> n) || n <= 0) return false;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    for (int i = 0; i < n; ++i) {
        std::string s;
        double x, y, z;
        if (!(f >> s >> x >> y >> z)) return false;
        auto it = el.find(s);
        if (it == el.end()) return false;
        Z.push_back(it->second);
        ang.insert(ang.end(), { x, y, z });
    }
    return true;
}

/// HF engine on the molecule (tight SCF); nullptr on failure.
std::unique_ptr<QMEngine> runHF(const std::vector<int>& Z, std::vector<double>& ang, const Args& a,
                                const std::string& grid)
{
    json cfg = { { "basis", "def2-SVP" }, { "scf_threshold", 1e-10 }, { "scf_energy_threshold", 1e-12 },
        { "scf_gradient_threshold", 1e-9 }, { "scf_max_iterations", 200 }, { "threads", a.threads }, { "grid", grid },
        { "grid_radial", a.radial }, { "grid_partition", a.partition }, { "grid_radii_adjust", a.adjust },
        { "cartesian_d", a.cartesian } };
    auto e = std::make_unique<QMEngine>(QMFunctional::HF, cfg);
    std::vector<int> z = Z;
    if (!e->QMInterface::InitialiseMolecule(z.data(), ang.data(), static_cast<int>(Z.size()), 0.0, 0)) return nullptr;
    e->Calculation(false);
    if (!e->scfConverged()) return nullptr;
    return e;
}

/// Real spherical harmonics up to lmax at unit vector u, index l*l + l + m (orthonormal on the sphere).
void realYlm(int lmax, double x, double y, double z, std::vector<double>& Y)
{
    Y.assign((lmax + 1) * (lmax + 1), 0.0);
    const double ct = std::max(-1.0, std::min(1.0, z));
    const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
    const double phi = std::atan2(y, x);
    // fully normalised associated Legendre functions Pbar_l^m (Y = Pbar e^{i m phi})
    std::vector<std::vector<double>> P(lmax + 1, std::vector<double>(lmax + 1, 0.0));
    P[0][0] = std::sqrt(1.0 / (4.0 * kPi));
    for (int m = 1; m <= lmax; ++m) P[m][m] = -std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * st * P[m - 1][m - 1];
    for (int m = 0; m < lmax; ++m) P[m + 1][m] = std::sqrt(2.0 * m + 3.0) * ct * P[m][m];
    for (int m = 0; m <= lmax; ++m)
        for (int l = m + 2; l <= lmax; ++l) {
            const double alm = std::sqrt((4.0 * l * l - 1.0) / (double(l) * l - double(m) * m));
            const double al1 = std::sqrt((4.0 * (l - 1) * (l - 1) - 1.0) / (double(l - 1) * (l - 1) - double(m) * m));
            P[l][m] = alm * (ct * P[l - 1][m] - P[l - 2][m] / al1);
        }
    for (int l = 0; l <= lmax; ++l) {
        Y[l * l + l] = P[l][0];
        for (int m = 1; m <= l; ++m) {
            Y[l * l + l + m] = std::sqrt(2.0) * P[l][m] * std::cos(m * phi);
            Y[l * l + l - m] = std::sqrt(2.0) * P[l][m] * std::sin(m * phi);
        }
    }
}

int modeLebedev()
{
    json out = json::array();
    const int sizes[] = { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202,
        1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 };
    for (int n : sizes) {
        Matrix u;
        Vector w;
        if (!qmgrid::lebedevGrid(n, u, w)) {
            std::cerr << "lebedev grid " << n << " missing\n";
            return 1;
        }
        const int L = qmgrid::lebedevDegree(n);
        const int lmax = L + 1;
        std::vector<double> integ((lmax + 1) * (lmax + 1), 0.0), Y;
        const bool ortho = n <= 590;
        const int lh = L / 2;
        const int nh = (lh + 1) * (lh + 1);
        Matrix S = Matrix::Zero(nh, nh);
        double wsum = 0.0;
        for (int p = 0; p < n; ++p) {
            realYlm(lmax, u(p, 0), u(p, 1), u(p, 2), Y);
            const double wp = 4.0 * kPi * w(p);
            wsum += w(p);
            for (size_t k = 0; k < Y.size(); ++k) integ[k] += wp * Y[k];
            if (ortho)
                for (int a = 0; a < nh; ++a)
                    for (int b = 0; b <= a; ++b) S(a, b) += wp * Y[a] * Y[b];
        }
        double err_le = 0.0, err_next = 0.0, err_ortho = 0.0;
        for (int l = 0; l <= lmax; ++l)
            for (int m = -l; m <= l; ++m) {
                const double ref = (l == 0) ? std::sqrt(4.0 * kPi) : 0.0;
                const double d = std::abs(integ[l * l + l + m] - ref);
                if (l <= L) err_le = std::max(err_le, d);
                else err_next = std::max(err_next, d);
            }
        if (ortho)
            for (int a = 0; a < nh; ++a)
                for (int b = 0; b <= a; ++b) err_ortho = std::max(err_ortho, std::abs(S(a, b) - (a == b ? 1.0 : 0.0)));
        json r = { { "npts", n }, { "degree", L }, { "weight_sum_minus_1", wsum - 1.0 }, { "max_err_l_le_degree", err_le },
            { "max_err_l_degree_plus_1", err_next } };
        if (ortho) {
            r["ortho_lmax"] = lh;
            r["max_err_orthonormality"] = err_ortho;
        }
        out.push_back(r);
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

json pointJson(const qmgrid::MolecularGrid& g, int p)
{
    const int A = g.atom[p];
    const int local = p - g.atom_offset[A];
    return { p, A, local / g.n_angular, local % g.n_angular, g.points(p, 0), g.points(p, 1), g.points(p, 2),
        g.weights(p), g.quadrature_weights(p) };
}

int modeGrid(const Args& a)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    const int nat = static_cast<int>(Z.size());
    Matrix pos(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) pos(i, k) = ang[3 * i + k] * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    qmgrid::GridOptions opt;
    opt.level = a.grid;
    opt.radial = a.radial;
    opt.partition = a.partition;
    opt.radii_adjust = a.adjust;
    qmgrid::MolecularGrid g;
    std::string why;
    const auto t0 = std::chrono::steady_clock::now();
    if (!qmgrid::buildMolecularGrid(Z, pos, opt, g, why, a.threads)) {
        std::cerr << "grid: " << why << "\n";
        return 2;
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    double wsum = 0.0, gauss = 0.0, pfun = 0.0;
    std::vector<double> atom_w(nat, 0.0);
    for (int p = 0; p < g.size(); ++p) {
        const double w = g.weights(p);
        wsum += w;
        atom_w[g.atom[p]] += w;
        for (int A = 0; A < nat; ++A) {
            const double dx = g.points(p, 0) - pos(A, 0), dy = g.points(p, 1) - pos(A, 1), dz = g.points(p, 2) - pos(A, 2);
            const double r2 = dx * dx + dy * dy + dz * dz;
            gauss += w * std::exp(-r2);
            pfun += w * dz * dz * std::exp(-2.0 * r2);  // integral = pi^(3/2) / 2^(7/2) per atom
        }
    }
    json out = { { "npts", g.size() }, { "n_radial", g.n_radial }, { "n_angular", g.n_angular }, { "weight_sum", wsum },
        { "gauss_integral", gauss }, { "gauss_exact", nat * std::pow(kPi, 1.5) }, { "pz2_integral", pfun },
        { "pz2_exact", nat * std::pow(kPi, 1.5) / std::pow(2.0, 3.5) }, { "atom_weight_sums", atom_w },
        { "build_ms", ms } };
    if (a.stride > 0) {
        json pts = json::array();
        for (int p = 0; p < g.size(); p += a.stride) pts.push_back(pointJson(g, p));
        out["points"] = pts;
        out["stride"] = a.stride;
    }
    if (!a.dump_all.empty()) {
        FILE* f = std::fopen(a.dump_all.c_str(), "w");
        if (!f) return 3;
        for (int p = 0; p < g.size(); ++p)
            std::fprintf(f, "%d %d %.17g %.17g %.17g %.17g\n", g.atom[p], p - g.atom_offset[g.atom[p]], g.points(p, 0),
                         g.points(p, 1), g.points(p, 2), g.weights(p));
        std::fclose(f);
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeDensity(const Args& a)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    auto e = runHF(Z, ang, a, a.grid);
    if (!e) {
        std::cerr << "density: HF failed\n";
        return 2;
    }
    qmgrid::MolecularGrid g;
    std::string why;
    if (!e->buildGrid(g, why)) {
        std::cerr << "density: " << why << "\n";
        return 2;
    }
    qmgrid::AOEvaluator ao(e->gtoBasis(), e->sphericalTransform());
    const double n_el = qmgrid::integrateDensity(g, ao, e->density(), a.threads);
    json out = { { "energy", e->TotalEnergy() }, { "num_electrons", e->numElectrons() }, { "npts", g.size() },
        { "nelec_grid", n_el } };
    if (a.stride > 0) {
        json pts = json::array();
        for (int p = 0; p < g.size(); p += a.stride) {
            qmgrid::AOBlock blk;
            ao.evaluate(g.points, p, 1, 1, blk);
            Vector rho;
            Matrix grad;
            qmgrid::densityOnBlock(blk, e->density(), rho, &grad);
            pts.push_back({ p, rho(0), grad(0, 0), grad(0, 1), grad(0, 2) });
        }
        out["points"] = pts;
        out["stride"] = a.stride;
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeNelec(const Args& a)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    auto e = runHF(Z, ang, a, "sg1");
    if (!e) {
        std::cerr << "nelec: HF failed\n";
        return 2;
    }
    qmgrid::AOEvaluator ao(e->gtoBasis(), e->sphericalTransform());
    json out = { { "energy", e->TotalEnergy() }, { "num_electrons", e->numElectrons() }, { "nbf", e->nbf() } };
    json rows = json::array();
    std::stringstream ss(a.levels);
    std::string lv;
    while (std::getline(ss, lv, ',')) {
        qmgrid::GridOptions opt;
        opt.level = lv;
        opt.radial = a.radial;
        opt.partition = a.partition;
        opt.radii_adjust = a.adjust;
        qmgrid::MolecularGrid g;
        std::string why;
        const auto t0 = std::chrono::steady_clock::now();
        if (!qmgrid::buildMolecularGrid(e->atomicNumbers(), e->atomPositionsBohr(), opt, g, why, a.threads)) {
            std::cerr << "nelec: " << why << "\n";
            return 2;
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double n = qmgrid::integrateDensity(g, ao, e->density(), a.threads);
        const auto t2 = std::chrono::steady_clock::now();
        rows.push_back({ { "level", lv }, { "npts", g.size() }, { "nelec_grid", n },
            { "deviation", n - e->numElectrons() },
            { "grid_ms", std::chrono::duration<double, std::milli>(t1 - t0).count() },
            { "rho_ms", std::chrono::duration<double, std::milli>(t2 - t1).count() } });
    }
    out["levels"] = rows;
    out["threads"] = a.threads;
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeAOFD(const Args& a)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    json cfg = { { "basis", "def2-SVP" }, { "cartesian_d", a.cartesian } };
    QMEngine e(QMFunctional::HF, cfg);
    std::vector<int> z = Z;
    if (!e.QMInterface::InitialiseMolecule(z.data(), ang.data(), static_cast<int>(Z.size()), 0.0, 0)) return 2;
    qmgrid::MolecularGrid g;
    std::string why;
    if (!e.buildGrid(g, why)) return 2;
    qmgrid::AOEvaluator ao(e.gtoBasis(), e.sphericalTransform());
    // 5-point central stencil, h = 1e-4 Bohr: error h^4 f^(5)/30. The error is reported relative to
    // max(1, |grad phi|) because tight core functions (a ~ 5e3) have large derivatives near their nucleus.
    const double h = 1.0e-4;
    double max_err = 0.0, max_val = 0.0;
    int npts = 0;
    for (int p = 0; p < g.size(); p += 97) {
        Matrix pt(13, 3);
        for (int r = 0; r < 13; ++r) pt.row(r) = g.points.row(p);
        for (int k = 0; k < 3; ++k) {
            pt(1 + 4 * k, k) += h;
            pt(2 + 4 * k, k) -= h;
            pt(3 + 4 * k, k) += 2.0 * h;
            pt(4 + 4 * k, k) -= 2.0 * h;
        }
        qmgrid::AOBlock b0, bs;
        ao.evaluate(pt, 0, 1, 1, b0);
        ao.evaluate(pt, 0, 13, 0, bs);
        for (int mu = 0; mu < ao.nbf(); ++mu) {
            const double an[3] = { b0.dx(0, mu), b0.dy(0, mu), b0.dz(0, mu) };
            for (int k = 0; k < 3; ++k) {
                const double fd = (8.0 * (bs.phi(1 + 4 * k, mu) - bs.phi(2 + 4 * k, mu)) -
                                      (bs.phi(3 + 4 * k, mu) - bs.phi(4 + 4 * k, mu))) / (12.0 * h);
                max_err = std::max(max_err, std::abs(fd - an[k]) / std::max(1.0, std::abs(an[k])));
                max_val = std::max(max_val, std::abs(an[k]));
            }
        }
        ++npts;
    }
    json out = { { "points_checked", npts }, { "nbf", ao.nbf() }, { "max_rel_err", max_err }, { "max_abs_grad", max_val } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}
}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 2) {
        std::cerr << "usage: test_qm_grid lebedev | grid|density|nelec|aofd <xyz> [options]\n";
        return 1;
    }
    Args a;
    a.mode = argv[1];
    int i = 2;
    if (a.mode != "lebedev") {
        if (argc < 3) return 1;
        a.xyz = argv[2];
        i = 3;
    }
    for (; i < argc; ++i) {
        const std::string s = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if (s == "--grid") a.grid = next();
        else if (s == "--radial") a.radial = next();
        else if (s == "--partition") a.partition = next();
        else if (s == "--adjust") a.adjust = next();
        else if (s == "--levels") a.levels = next();
        else if (s == "--stride") a.stride = std::stoi(next());
        else if (s == "--threads") a.threads = std::stoi(next());
        else if (s == "--dump-all") a.dump_all = next();
        else if (s == "--cartesian") a.cartesian = true;
        else {
            std::cerr << "unknown option " << s << "\n";
            return 1;
        }
    }
    if (a.mode == "lebedev") return modeLebedev();
    if (a.mode == "grid") return modeGrid(a);
    if (a.mode == "density") return modeDensity(a);
    if (a.mode == "nelec") return modeNelec(a);
    if (a.mode == "aofd") return modeAOFD(a);
    std::cerr << "unknown mode " << a.mode << "\n";
    return 1;
}
