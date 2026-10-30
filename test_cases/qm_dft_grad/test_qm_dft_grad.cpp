/*
 * Native QM Kohn-Sham gradient test driver (WP8 DFT), derived from test_qm_b3lyp.cpp (WP7).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Modes (JSON on stdout, gradients in Eh/Bohr, atom-major; check_qm_dft_grad.py compares with reference.json from
 * scripts/qm_dft_grad_pyscf_reference.py):
 *
 *   test_qm_dft_grad grad <xyz> [--functional lda|pbe|b3lyp|tpss|scan|hf] [--vwn vwn5|vwn_rpa] [--grid sg1] [--ri off|j|jk]
 *                         [--response 1|0] [--threads N] [--fd moving|frozen] [--h 2e-3] [--fd_atoms N]
 *       KS/def2-SVP SCF (scf_threshold 1e-10, energy 1e-12, orbital gradient 1e-9) and its analytic gradient, parts
 *       {1e, 2e, nuc, xc basis, grid response}, sum over atoms and total torque sum_A (R_A - R_c) x g_A.
 *       --fd moving: central 4-point finite difference of the energy (grid rebuilt at every displaced geometry,
 *                    i.e. moving with the atoms) -> must equal the gradient WITH grid response;
 *       --fd frozen: the same with the grid of the reference geometry kept fixed in space (QMEngine::freezeXCGrid)
 *                    -> must equal the gradient WITHOUT grid response (use --response 0).
 *       --fd_atoms N: only the first N atoms are displaced (cost).
 *
 *   test_qm_dft_grad rotate <xyz> [options as grad]
 *       Gradient at the geometry rotated by fixed Euler angles (the grid does not rotate with the molecule), rotated
 *       back, against the gradient at the original geometry; also the energy difference.
 *
 *   test_qm_dft_grad threads <xyz> [options]   gradient with 1 and 4 threads (max |difference|, bitwise target 0)
 *   test_qm_dft_grad update <xyzA> <xyzB> [options]   UpdateMolecule A -> B against a fresh engine on B
 *
 * Claude Generated (Oct 2026, WP8). GPL-3.0.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/energy_calculators/qm_methods/qm_grid.h"
#include "src/core/units.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string mode, xyz, xyz2, functional = "b3lyp", vwn = "vwn5", grid = "sg1", ri = "off", fd;
    int threads = 4, fd_atoms = -1;
    bool response = true;
    double h = 2.0e-3;  // Bohr
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

json engineConfig(const Args& a, int threads)
{
    return { { "basis", "def2-SVP" }, { "scf_threshold", 1e-10 }, { "scf_energy_threshold", 1e-12 },
        { "scf_gradient_threshold", 1e-9 }, { "scf_gradient_threshold_tight", 1e-9 }, { "scf_max_iterations", 200 },
        { "threads", threads }, { "grid", a.grid }, { "ri", a.ri }, { "vwn_variant", a.vwn },
        { "grid_response", a.response } };
}

QMFunctional functionalOf(const std::string& s)
{
    // tpss/scan: WP11 stage 2 (meta-GGA gradient), Claude Generated (Oct 2026)
    return s == "lda" ? QMFunctional::LDA : s == "pbe" ? QMFunctional::PBE : s == "hf" ? QMFunctional::HF
         : s == "tpss" ? QMFunctional::TPSS : s == "scan" ? QMFunctional::SCAN : QMFunctional::B3LYP;
}

std::unique_ptr<QMEngine> makeEngine(const Args& a, const std::vector<int>& Z, const std::vector<double>& ang,
                                     int threads)
{
    auto e = std::make_unique<QMEngine>(functionalOf(a.functional), engineConfig(a, threads));
    std::vector<double> xyz = ang;
    std::vector<int> z = Z;
    if (!e->QMInterface::InitialiseMolecule(z.data(), xyz.data(), static_cast<int>(Z.size()), 0.0, 0)) return nullptr;
    return e;
}

std::vector<double> flat(const Matrix& g)
{
    std::vector<double> v;
    for (int i = 0; i < g.rows(); ++i)
        for (int k = 0; k < 3; ++k) v.push_back(g(i, k));
    return v;
}

/// Gradient (Eh/Bohr) and energy at the geometry `ang`; empty matrix on failure.
bool gradientAt(const Args& a, const std::vector<int>& Z, const std::vector<double>& ang, int threads, double& E,
                Matrix& g, json* parts = nullptr)
{
    auto e = makeEngine(a, Z, ang, threads);
    if (!e) return false;
    E = e->Calculation(true);
    if (!e->scfConverged() || e->gradientFailed()) {
        std::cerr << "gradient failed: " << e->gradientError() << "\n";
        return false;
    }
    g = e->gradientBohr();
    if (parts) {
        const char* names[5] = { "one_electron", "two_electron", "nuclear", "xc_basis", "grid_response" };
        for (size_t k = 0; k < e->gradientParts().size() && k < 5; ++k) (*parts)[names[k]] = flat(e->gradientParts()[k]);
    }
    return true;
}

json sumAndTorque(const Matrix& g, const std::vector<double>& ang)
{
    const int n = static_cast<int>(g.rows());
    double c[3] = { 0, 0, 0 }, s[3] = { 0, 0, 0 }, t[3] = { 0, 0, 0 };
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) c[k] += ang[3 * i + k] / n;
    for (int i = 0; i < n; ++i) {
        double r[3];
        for (int k = 0; k < 3; ++k) {
            r[k] = CurcumaUnit::Length::ANGSTROM_TO_BOHR * (ang[3 * i + k] - c[k]);
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
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    const int n = static_cast<int>(Z.size());
    double E = 0.0;
    Matrix g;
    json parts = json::object();
    if (!gradientAt(a, Z, ang, a.threads, E, g, &parts)) return 2;
    json out = { { "energy", E }, { "natoms", n }, { "gradient", flat(g) }, { "parts", parts },
        { "grid_response", a.response } };
    out.update(sumAndTorque(g, ang));
    if (!a.fd.empty()) {
        // grid of the reference geometry for --fd frozen
        qmgrid::MolecularGrid grid0;
        if (a.fd == "frozen") {
            auto e0 = makeEngine(a, Z, ang, a.threads);
            std::string why;
            if (!e0 || !e0->buildGrid(grid0, why)) return 1;
        }
        Args af = a;
        af.response = false;  // energies only; the frozen grid forbids the response
        auto energy = [&](int atom, int k, double dB) {
            std::vector<double> x = ang;
            x[3 * atom + k] += dB * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
            auto e = makeEngine(af, Z, x, a.threads);
            if (a.fd == "frozen") e->freezeXCGrid(grid0);
            const double v = e->Calculation(false);
            if (!e->scfConverged()) std::cerr << "FD SCF not converged\n";
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
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    const int n = static_cast<int>(Z.size());
    // R = Rz(0.7) Ry(0.4) Rx(1.1), fixed angles (no RNG)
    auto rot = [](double ax, double ay, double az) {
        Eigen::Matrix3d Rx, Ry, Rz;
        Rx << 1, 0, 0, 0, std::cos(ax), -std::sin(ax), 0, std::sin(ax), std::cos(ax);
        Ry << std::cos(ay), 0, std::sin(ay), 0, 1, 0, -std::sin(ay), 0, std::cos(ay);
        Rz << std::cos(az), -std::sin(az), 0, std::sin(az), std::cos(az), 0, 0, 0, 1;
        return Eigen::Matrix3d(Rz * Ry * Rx);
    };
    const Eigen::Matrix3d R = rot(1.1, 0.4, 0.7);
    std::vector<double> rang(3 * n), tang(3 * n);
    for (int i = 0; i < n; ++i) {
        const Eigen::Vector3d v(ang[3 * i], ang[3 * i + 1], ang[3 * i + 2]);
        const Eigen::Vector3d w = R * v;
        for (int k = 0; k < 3; ++k) {
            rang[3 * i + k] = w(k);
            tang[3 * i + k] = ang[3 * i + k] + (k == 0 ? 0.37 : k == 1 ? -0.21 : 0.53);  // translation (Angstrom)
        }
    }
    double E0, Er, Et;
    Matrix g0, gr, gt;
    if (!gradientAt(a, Z, ang, a.threads, E0, g0) || !gradientAt(a, Z, rang, a.threads, Er, gr)
        || !gradientAt(a, Z, tang, a.threads, Et, gt))
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
    out["rotated"] = sumAndTorque(gr, rang);
    out["original"] = sumAndTorque(g0, ang);
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeThreads(const Args& a)
{
    std::vector<int> Z;
    std::vector<double> ang;
    if (!readXYZ(a.xyz, Z, ang)) return 1;
    double E1, E4;
    Matrix g1, g4;
    if (!gradientAt(a, Z, ang, 1, E1, g1) || !gradientAt(a, Z, ang, 4, E4, g4)) return 2;
    json out = { { "energy_1", E1 }, { "energy_4", E4 }, { "max_dg", (g1 - g4).cwiseAbs().maxCoeff() },
        { "bitwise", (g1.array() == g4.array()).all() && E1 == E4 } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

int modeUpdate(const Args& a)
{
    std::vector<int> ZA, ZB;
    std::vector<double> angA, angB;
    if (!readXYZ(a.xyz, ZA, angA) || !readXYZ(a.xyz2, ZB, angB) || ZA != ZB) return 1;
    auto e = makeEngine(a, ZA, angA, a.threads);
    if (!e) return 1;
    e->Calculation(true);
    const Matrix gA = e->gradientBohr();
    const int n = static_cast<int>(ZB.size());
    Matrix geo(n, 3);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) geo(i, k) = angB[3 * i + k];
    e->QMInterface::UpdateMolecule(geo);
    const double EBu = e->Calculation(true);
    const Matrix gBu = e->gradientBohr();
    double EBf;
    Matrix gBf;
    if (!gradientAt(a, ZB, angB, a.threads, EBf, gBf)) return 2;
    json out = { { "energy_B_updated", EBu }, { "energy_B_fresh", EBf }, { "max_dg", (gBu - gBf).cwiseAbs().maxCoeff() },
        { "max_g_A_minus_B", (gA - gBf).cwiseAbs().maxCoeff() }, { "converged", e->scfConverged() } };
    std::cout << out.dump(1) << std::endl;
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    if (argc < 3) {
        std::cerr << "usage: test_qm_dft_grad grad|rotate|threads <xyz> | update <xyzA> <xyzB> [options]\n";
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
        if (s == "--grid") a.grid = next();
        else if (s == "--functional") a.functional = next();
        else if (s == "--vwn") a.vwn = next();
        else if (s == "--ri") a.ri = next();
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
        if (a.mode == "grad") return modeGrad(a);
        if (a.mode == "rotate") return modeRotate(a);
        if (a.mode == "threads") return modeThreads(a);
        if (a.mode == "update") return modeUpdate(a);
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 3;
    }
    std::cerr << "unknown mode " << a.mode << "\n";
    return 1;
}
