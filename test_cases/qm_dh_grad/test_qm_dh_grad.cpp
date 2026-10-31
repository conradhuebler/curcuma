/*
 * Native double-hybrid gradient test driver (WP12 stage 2: analytic gradient, CPU).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout, gradients in Eh/Bohr atom-major; check_qm_dh_grad.py applies the gates):
 *
 *   test_qm_dh_grad kernel <reference.json>
 *       At the (rho, sigma[, tau]) points of reference.json["kernel"]: for B88_X, LYP_C, PBE_X, P86_C, P86VWN_C, B95_C
 *       the jet first derivatives (qm_xc_fxc_kernels.h) next to the production potentials (qm_xc_kernels.h), the jet second
 *       derivatives and a central 4-point finite difference of the production potentials (test tool only).
 *
 *   test_qm_dh_grad fxc <xyz> --method M [options]
 *       After the SCF: F_xc[D] (qmxc::buildFxc) against (V_xc[P + h D] - V_xc[P - h D]) / 2h (4-point) for a fixed
 *       symmetric D (finite difference = test tool only, never the production path); also the directional XC gradient
 *       (qmxc::buildXCResponseGradient: basis, points, weight energy density) against the finite difference of
 *       qmxc::buildXCGradient in direction D.
 *
 *   test_qm_dh_grad grad <xyz> --method M [--pt2 ri|canonical] [--fc] [--grid sg1] [--response 1|0] [--threads N]
 *                        [--p86 pz81|vwn5] [--fd moving|frozen] [--h 1e-3] [--fd_atoms N]
 *       DoubleHybridMethod energy and analytic gradient, parts, sum over atoms, torque; with --fd the central 4-point
 *       finite difference of the total energy (moving grid = with grid response, frozen grid = without).
 *
 *   test_qm_dh_grad rotate|threads <xyz> [options]   rotated/translated geometry; 1 vs 4 vs 8 threads (bitwise)
 *   test_qm_dh_grad update <xyzA> <xyzB> [options]   updateGeometry A -> B against a fresh method on B
 *
 * Claude Generated (Oct 2026, WP12 stage 2). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/double_hybrid_method.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/energy_calculators/qm_methods/qm_grid.h"
#include "src/core/energy_calculators/qm_methods/qm_xc_fxc.h"
#include "src/core/units.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// Driver defaults = the PySCF/libxc reference convention (PZ81-P86, no frozen core), always passed explicitly to
// DoubleHybridMethod; curcuma's own DH defaults are vwn5 + frozen core since Oct 31 2026.
struct Args {
    std::string mode, xyz, xyz2, method = "b2plyp", grid = "sg1", pt2 = "ri", p86 = "pz81", fd;
    int threads = 4, fd_atoms = -1;
    bool fc = false, response = true;
    double h = 1.0e-3;  // Bohr
    Args();
    Args(const Args&);
    ~Args();  // out of line: the implicit inline one trips -Winline
};
Args::Args() = default;
Args::Args(const Args&) = default;
Args::~Args() = default;

bool readMol(const std::string& path, Mol& m)
{
    static const std::map<std::string, int> el{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
        { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 } };
    std::ifstream f(path);
    int n = 0;
    if (!(f >> n) || n <= 0) return false;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    m = Mol();
    m.m_number_atoms = n;
    m.m_geometry = Geometry::Zero(n, 3);
    for (int i = 0; i < n; ++i) {
        std::string s;
        double x, y, z;
        if (!(f >> s >> x >> y >> z)) return false;
        auto it = el.find(s);
        if (it == el.end()) return false;
        m.m_atoms.push_back(it->second);
        m.m_geometry(i, 0) = x;
        m.m_geometry(i, 1) = y;
        m.m_geometry(i, 2) = z;
    }
    m.m_charge = 0;
    m.m_spin = 0;
    return true;
}

json methodConfig(const Args& a, int threads)
{
    json qm = { { "basis", "def2-SVP" }, { "scf_threshold", 1e-10 }, { "scf_energy_threshold", 1e-12 },
        { "scf_gradient_threshold", 1e-9 }, { "scf_gradient_threshold_tight", 1e-9 }, { "scf_max_iterations", 300 },
        { "threads", threads }, { "grid", a.grid }, { "frozen_core", a.fc }, { "dh_pt2", a.pt2 }, { "p86_lda", a.p86 },
        { "grid_response", a.response } };
    return { { "qm", qm } };
}

std::vector<double> flat(const Matrix& g)
{
    std::vector<double> v;
    for (int i = 0; i < g.rows(); ++i)
        for (int k = 0; k < 3; ++k) v.push_back(g(i, k));
    return v;
}

/// Energy and gradient (Eh/Bohr) of the method at geometry m; false on error.
bool gradientAt(const Args& a, const Mol& m, int threads, double& E, Matrix& g, json* parts = nullptr, json* info = nullptr)
{
    DoubleHybridMethod dh(a.method, methodConfig(a, threads));
    if (!dh.setMolecule(m)) return false;
    E = dh.calculateEnergy(true);
    if (dh.hasError()) {
        std::cerr << "gradient failed: " << dh.getErrorMessage() << "\n";
        return false;
    }
    g = dh.gradientBohr();
    if (parts) {
        const char* names[6] = { "one_electron", "two_electron", "nuclear", "xc_basis", "grid_response", "dispersion" };
        for (size_t k = 0; k < dh.gradientParts().size() && k < 6; ++k) (*parts)[names[k]] = flat(dh.gradientParts()[k]);
    }
    if (info) {
        const qmint::MP2Result& r = dh.pt2GradientResult();
        *info = { { "zvector_iterations", r.zvector_iterations }, { "zvector_residual", r.zvector_residual },
            { "e_os_gradient_path", r.e_os }, { "e_ss_gradient_path", r.e_ss }, { "e_os", dh.pt2Result().e_os },
            { "e_ss", dh.pt2Result().e_ss }, { "decomposition", dh.getEnergyDecomposition() },
            { "t_zvector_ms", r.t_grad_zvector_ms }, { "t_deriv2e_ms", r.t_grad_deriv2e_ms }, { "info", r.info } };
    }
    return true;
}

json sumAndTorque(const Matrix& g, const Mol& m)
{
    const int n = static_cast<int>(g.rows());
    double c[3] = { 0, 0, 0 }, s[3] = { 0, 0, 0 }, t[3] = { 0, 0, 0 };
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) c[k] += m.m_geometry(i, k) / n;
    for (int i = 0; i < n; ++i) {
        double r[3];
        for (int k = 0; k < 3; ++k) {
            r[k] = CurcumaUnit::Length::ANGSTROM_TO_BOHR * (m.m_geometry(i, k) - c[k]);
            s[k] += g(i, k);
        }
        t[0] += r[1] * g(i, 2) - r[2] * g(i, 1);
        t[1] += r[2] * g(i, 0) - r[0] * g(i, 2);
        t[2] += r[0] * g(i, 1) - r[1] * g(i, 0);
    }
    return { { "force_sum", { s[0], s[1], s[2] } }, { "torque", { t[0], t[1], t[2] } } };
}

int modeGrad(const Args& a)
{
    Mol m;
    if (!readMol(a.xyz, m)) return 1;
    const int n = m.AtomCount();
    double E = 0.0;
    Matrix g;
    json parts = json::object(), info;
    if (!gradientAt(a, m, a.threads, E, g, &parts, &info)) return 2;
    json out = { { "energy", E }, { "natoms", n }, { "gradient", flat(g) }, { "parts", parts }, { "info", info },
        { "grid_response", a.response } };
    out.update(sumAndTorque(g, m));
    if (!a.fd.empty()) {
        qmgrid::MolecularGrid grid0;
        if (a.fd == "frozen") {
            DoubleHybridMethod d0(a.method, methodConfig(a, a.threads));
            d0.setMolecule(m);
            std::string why;
            if (!d0.engine().buildGrid(grid0, why)) return 1;
        }
        Args af = a;
        af.response = false;
        auto energy = [&](int atom, int k, double dB) {
            Mol x = m;
            x.m_geometry(atom, k) += dB * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
            DoubleHybridMethod dh(a.method, methodConfig(af, a.threads));
            dh.setMolecule(x);
            if (a.fd == "frozen") dh.engine().freezeXCGrid(grid0);
            const double v = dh.calculateEnergy(false);
            if (dh.hasError()) std::cerr << "FD energy failed: " << dh.getErrorMessage() << "\n";
            return v;
        };
        const int nd = a.fd_atoms > 0 ? std::min(a.fd_atoms, n) : n;
        std::vector<double> fd;
        for (int i = 0; i < nd; ++i)
            for (int k = 0; k < 3; ++k) {
                const double h = a.h;
                fd.push_back((-energy(i, k, 2 * h) + 8 * energy(i, k, h) - 8 * energy(i, k, -h) + energy(i, k, -2 * h))
                             / (12 * h));
            }
        out["fd_gradient"] = fd;
        out["fd_mode"] = a.fd;
        out["fd_step_bohr"] = a.h;
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeRotate(const Args& a)
{
    Mol m;
    if (!readMol(a.xyz, m)) return 1;
    const int n = m.AtomCount();
    auto rot = [](double ax, double ay, double az) {
        Eigen::Matrix3d Rx, Ry, Rz;
        Rx << 1, 0, 0, 0, std::cos(ax), -std::sin(ax), 0, std::sin(ax), std::cos(ax);
        Ry << std::cos(ay), 0, std::sin(ay), 0, 1, 0, -std::sin(ay), 0, std::cos(ay);
        Rz << std::cos(az), -std::sin(az), 0, std::sin(az), std::cos(az), 0, 0, 0, 1;
        return Eigen::Matrix3d(Rz * Ry * Rx);
    };
    const Eigen::Matrix3d R = rot(1.1, 0.4, 0.7);
    Mol mr = m, mt = m;
    for (int i = 0; i < n; ++i) {
        const Eigen::Vector3d v(m.m_geometry(i, 0), m.m_geometry(i, 1), m.m_geometry(i, 2));
        const Eigen::Vector3d w = R * v;
        for (int k = 0; k < 3; ++k) {
            mr.m_geometry(i, k) = w(k);
            mt.m_geometry(i, k) = m.m_geometry(i, k) + (k == 0 ? 0.37 : k == 1 ? -0.21 : 0.53);
        }
    }
    double E0, Er, Et;
    Matrix g0, gr, gt;
    if (!gradientAt(a, m, a.threads, E0, g0) || !gradientAt(a, mr, a.threads, Er, gr) || !gradientAt(a, mt, a.threads, Et, gt))
        return 2;
    Matrix gback(n, 3);
    for (int i = 0; i < n; ++i) {
        const Eigen::Vector3d v(gr(i, 0), gr(i, 1), gr(i, 2));
        const Eigen::Vector3d w = R.transpose() * v;
        for (int k = 0; k < 3; ++k) gback(i, k) = w(k);
    }
    json out = { { "energy", E0 }, { "energy_rotated", Er }, { "energy_translated", Et },
        { "rotation_max_dg", (gback - g0).cwiseAbs().maxCoeff() }, { "translation_max_dg", (gt - g0).cwiseAbs().maxCoeff() },
        { "max_g", g0.cwiseAbs().maxCoeff() } };
    out["rotated"] = sumAndTorque(gr, mr);
    out["original"] = sumAndTorque(g0, m);
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeThreads(const Args& a)
{
    Mol m;
    if (!readMol(a.xyz, m)) return 1;
    // 1 vs 4 vs 8 threads (8 added Oct 31 2026 with the fixed chunks of mp2GradientHost)
    double E1, E4, E8;
    Matrix g1, g4, g8;
    if (!gradientAt(a, m, 1, E1, g1) || !gradientAt(a, m, 4, E4, g4) || !gradientAt(a, m, 8, E8, g8)) return 2;
    json out = { { "energy_1", E1 }, { "energy_4", E4 }, { "energy_8", E8 },
        { "max_dg", std::max((g1 - g4).cwiseAbs().maxCoeff(), (g1 - g8).cwiseAbs().maxCoeff()) },
        { "bitwise", (g1.array() == g4.array()).all() && (g1.array() == g8.array()).all() && E1 == E4 && E1 == E8 } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeUpdate(const Args& a)
{
    Mol A, B;
    if (!readMol(a.xyz, A) || !readMol(a.xyz2, B)) return 1;
    DoubleHybridMethod dh(a.method, methodConfig(a, a.threads));
    dh.setMolecule(A);
    dh.calculateEnergy(true);
    dh.updateGeometry(B.m_geometry);
    const double EBu = dh.calculateEnergy(true);
    const Matrix gBu = dh.gradientBohr();
    double EBf;
    Matrix gBf;
    if (dh.hasError() || !gradientAt(a, B, a.threads, EBf, gBf)) return 2;
    json out = { { "energy_B_updated", EBu }, { "energy_B_fresh", EBf }, { "max_dg", (gBu - gBf).cwiseAbs().maxCoeff() } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

// ---- kernel jets vs production kernels and FD of the production potentials ----
void productionV(qmxc::XCKernel k, double r, double s, double t, double v[3])
{
    double e = 0.0;
    v[0] = v[1] = v[2] = 0.0;
    if (qmxc::kern::isMGGA(k))
        qmxc::kern::evaluateKernelMGGA(k, r, s, t, e, v[0], v[1], v[2]);
    else
        qmxc::kern::evaluateKernel(k, r, s, e, v[0], v[1]);
}

int modeKernel(const std::string& ref_path)
{
    std::ifstream f(ref_path);
    json ref;
    try {
        f >> ref;
    } catch (const std::exception& ex) {
        std::cerr << "cannot read " << ref_path << ": " << ex.what() << "\n";
        return 1;
    }
    const json& K = ref["kernel"];
    const std::vector<double> rho = K["rho"].get<std::vector<double>>(), sigma = K["sigma"].get<std::vector<double>>(),
                              tau = K["tau"].get<std::vector<double>>();
    const std::pair<const char*, qmxc::XCKernel> kernels[] = { { "b88", qmxc::XCKernel::B88_X },
        { "lyp", qmxc::XCKernel::LYP_C }, { "pbex", qmxc::XCKernel::PBE_X }, { "p86", qmxc::XCKernel::P86_C },
        { "p86vwn", qmxc::XCKernel::P86VWN_C }, { "b95", qmxc::XCKernel::B95_C } };
    json out = json::object();
    for (const auto& kk : kernels) {
        const bool mg = qmxc::kern::isMGGA(kk.second);
        json jv = json::array(), pv = json::array(), jh = json::array(), fh = json::array(), fn = json::array();
        for (size_t i = 0; i < rho.size(); ++i) {
            const double r = rho[i], s = sigma[i], t = mg ? tau[i] : 0.0;
            qmxc::fxc::Jet e;
            qmxc::fxc::kernelJet(kk.second, r, s, t, e);
            double v[3];
            productionV(kk.second, r, s, t, v);
            jv.push_back({ e.g[0], e.g[1], e.g[2] });
            pv.push_back({ v[0], v[1], v[2] });
            jh.push_back({ e.h[0], e.h[1], e.h[2], e.h[3], e.h[4], e.h[5] });
            // FD of the production potentials: column y of the Hessian = d v / d y (4-point). Relative steps 1e-2 ..
            // 1e-5 of the variable; per entry the one closest to the jet value is kept (at small sigma the short
            // steps are rounding-limited, the long ones truncation-limited); a test tool, not the production path.
            const double x0[3] = { r, s, t };
            double H[3][3] = { { 0 } }, Nz[3][3] = { { 0 } };  // Nz: rounding-noise estimate of the chosen FD entry
            const int nv = mg ? 3 : 2;
            bool ok = true;
            const int I3[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
            for (int y = 0; y < nv; ++y) {
                bool first = true;
                for (double rel : { 1e-2, 1e-3, 1e-4, 1e-5 }) {
                    const double hstep = rel * std::max(std::abs(x0[y]), 1e-30);
                    double vp1[3], vm1[3], vp2[3], vm2[3];
                    double x[3] = { r, s, t };
                    x[y] = x0[y] + hstep; productionV(kk.second, x[0], x[1], x[2], vp1);
                    x[y] = x0[y] - hstep; productionV(kk.second, x[0], x[1], x[2], vm1);
                    x[y] = x0[y] + 2 * hstep; productionV(kk.second, x[0], x[1], x[2], vp2);
                    x[y] = x0[y] - 2 * hstep; productionV(kk.second, x[0], x[1], x[2], vm2);
                    if (x0[y] - 2 * hstep <= 0.0) ok = false;
                    for (int q = 0; q < 3; ++q) {
                        const double d = (-vp2[q] + 8 * vp1[q] - 8 * vm1[q] + vm2[q]) / (12 * hstep);
                        if (first || std::abs(d - e.h[I3[q][y]]) < std::abs(H[q][y] - e.h[I3[q][y]])) {
                            H[q][y] = d;
                            const double vm = std::max({ std::abs(vp1[q]), std::abs(vm1[q]), std::abs(vp2[q]), std::abs(vm2[q]) });
                            Nz[q][y] = 18.0 * 2.2e-16 * vm / (12.0 * hstep);
                        }
                    }
                    first = false;
                }
            }
            if (ok) {
                fh.push_back({ H[0][0], H[0][1], H[0][2], H[1][1], H[1][2], H[2][2] });
                fn.push_back({ Nz[0][0], Nz[0][1], Nz[0][2], Nz[1][1], Nz[1][2], Nz[2][2] });
            } else {
                fh.push_back(nullptr);
                fn.push_back(nullptr);
            }
        }
        out[kk.first] = { { "jet_v", jv }, { "prod_v", pv }, { "jet_h", jh }, { "fd_h", fh }, { "fd_noise", fn } };
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}

// ---- F_xc[D] and the directional XC gradient against finite differences of the production builders ----
int modeFxc(const Args& a)
{
    Mol m;
    if (!readMol(a.xyz, m)) return 1;
    Args ar = a;
    DoubleHybridMethod dh(a.method, methodConfig(ar, a.threads));
    if (!dh.setMolecule(m)) return 2;
    dh.calculateEnergy(false);
    if (dh.hasError()) return 2;
    QMEngine& eng = dh.engine();
    qmgrid::MolecularGrid grid;
    std::string why;
    if (!eng.buildGrid(grid, why)) return 2;
    const qmgrid::AOEvaluator ao(eng.gtoBasis(), eng.sphericalTransform());
    const Matrix P = eng.density();
    const int n = (int)P.rows();
    // fixed symmetric direction: D_mn = 0.05 cos(m + 2n) cos(n + 2m) / (1 + |m - n|) scaled with sqrt(|P_mm P_nn|) + 1e-3
    Matrix D(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            D(i, j) = 0.05 * std::cos(i + 2.0 * j) * std::cos(j + 2.0 * i) / (1.0 + std::abs(i - j))
                      * (std::sqrt(std::abs(P(i, i) * P(j, j))) + 1e-3);
    const qmxc::XCFunctional& f = eng.xcFunctional();
    const double thr = 1e-12;
    Matrix V1;
    if (!qmxc::buildFxc(grid, ao, P, D, f, thr, a.threads, V1, nullptr)) return 3;
    const double h = 1e-4;
    auto vxc = [&](double t) {
        qmxc::XCResult r;
        qmxc::buildVxc(grid, ao, Matrix(P + t * D), f, thr, a.threads, r, 128, nullptr);
        return r.V;
    };
    const Matrix Vfd = (-vxc(2 * h) + 8 * vxc(h) - 8 * vxc(-h) + vxc(-2 * h)) / (12 * h);
    // directional XC gradient
    std::vector<int> bf_atom(n, 0);
    const Matrix& Q = eng.sphericalTransform();
    for (int c = 0; c < n; ++c) {
        if (Q.size() == 0) {
            bf_atom[c] = eng.gtoBasis()[c].atom;
            continue;
        }
        for (int r = 0; r < Q.rows(); ++r)
            if (Q(r, c) != 0.0) {
                bf_atom[c] = eng.gtoBasis()[r].atom;
                break;
            }
    }
    const int nat = m.AtomCount();
    qmxc::XCGradient gl;
    if (!qmxc::buildXCResponseGradient(grid, ao, P, D, f, thr, a.threads, bf_atom, nat, gl, nullptr)) return 3;
    auto xcg = [&](double t) {
        qmxc::XCGradient g;
        qmxc::buildXCGradient(grid, ao, Matrix(P + t * D), f, thr, a.threads, bf_atom, nat, g, 128, nullptr);
        return g;
    };
    const qmxc::XCGradient gp2 = xcg(2 * h), gp1 = xcg(h), gm1 = xcg(-h), gm2 = xcg(-2 * h);
    const Matrix bfd = (-gp2.basis + 8 * gp1.basis - 8 * gm1.basis + gm2.basis) / (12 * h);
    const Matrix pfd = (-gp2.points + 8 * gp1.points - 8 * gm1.points + gm2.points) / (12 * h);
    const Vector efd = (-gp2.edens + 8 * gp1.edens - 8 * gm1.edens + gm2.edens) / (12 * h);
    json out = { { "fxc_max_diff", (V1 - Vfd).cwiseAbs().maxCoeff() }, { "fxc_max", V1.cwiseAbs().maxCoeff() },
        { "basis_max_diff", (gl.basis - bfd).cwiseAbs().maxCoeff() }, { "basis_max", gl.basis.cwiseAbs().maxCoeff() },
        { "points_max_diff", (gl.points - pfd).cwiseAbs().maxCoeff() }, { "points_max", gl.points.cwiseAbs().maxCoeff() },
        { "edens_max_diff", (gl.edens - efd).cwiseAbs().maxCoeff() }, { "edens_max", gl.edens.cwiseAbs().maxCoeff() },
        { "fd_step", h } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

/// Energy scan along one cartesian coordinate with the grid frozen in space (diagnostic of energy discontinuities,
/// e.g. the r_s = 1 seam of PZ81): E(x0 + k h), k = -n..n, atom --fd_atoms (index), component 0..2 from --threads.
int modeScan(const Args& a)
{
    Mol m;
    if (!readMol(a.xyz, m)) return 1;
    const int atom = std::max(0, a.fd_atoms), comp = 2, nsteps = 20;
    qmgrid::MolecularGrid grid0;
    {
        DoubleHybridMethod d0(a.method, methodConfig(a, 4));
        d0.setMolecule(m);
        std::string why;
        if (!d0.engine().buildGrid(grid0, why)) return 1;
    }
    Args af = a;
    af.response = false;
    json E = json::array();
    for (int k = -nsteps; k <= nsteps; ++k) {
        Mol x = m;
        x.m_geometry(atom, comp) += k * a.h * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
        DoubleHybridMethod dh(a.method, methodConfig(af, 4));
        dh.setMolecule(x);
        dh.engine().freezeXCGrid(grid0);
        E.push_back(dh.calculateEnergy(false));
    }
    std::cout << json({ { "step_bohr", a.h }, { "atom", atom }, { "energies", E } }).dump(1) << std::endl;
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 3) {
        std::cerr << "usage: test_qm_dh_grad kernel <ref.json> | fxc|grad|rotate|threads <xyz> | update <xyzA> <xyzB> [options]\n";
        return 1;
    }
    Args a;
    a.mode = argv[1];
    a.xyz = argv[2];
    int i = 3;
    if (a.mode == "update") {
        if (argc < 4) return 1;
        a.xyz2 = argv[3];
        i = 4;
    }
    for (; i < argc; ++i) {
        const std::string s = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if (s == "--method") a.method = next();
        else if (s == "--grid") a.grid = next();
        else if (s == "--pt2") a.pt2 = next();
        else if (s == "--p86") a.p86 = next();
        else if (s == "--fc") a.fc = true;
        else if (s == "--threads") a.threads = std::stoi(next());
        else if (s == "--response") a.response = std::stoi(next()) != 0;
        else if (s == "--fd") a.fd = next();
        else if (s == "--h") a.h = std::stod(next());
        else if (s == "--fd_atoms") a.fd_atoms = std::stoi(next());
        else {
            std::cerr << "unknown option " << s << "\n";
            return 1;
        }
    }
    try {
        if (a.mode == "kernel") return modeKernel(a.xyz);
        if (a.mode == "fxc") return modeFxc(a);
        if (a.mode == "grad") return modeGrad(a);
        if (a.mode == "rotate") return modeRotate(a);
        if (a.mode == "threads") return modeThreads(a);
        if (a.mode == "update") return modeUpdate(a);
        if (a.mode == "scan") return modeScan(a);
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 3;
    }
    std::cerr << "unknown mode " << a.mode << "\n";
    return 1;
}
