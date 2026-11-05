/*
 * Native QM analytic RHF Hessian dumper (validation of qm_hessian.cpp).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reads an XYZ, runs `hf`, computes the analytic Hessian (QMEngine::computeHessian) and prints JSON (Eh/Bohr^2):
 *
 *   dump_qm_hessian <xyz> [--basis NAME] [--fd] [--h 5e-4] [--stencil 2|4] [--threads N] [--pieces]
 *                   [--newton] [--qm KEY VALUE] [--cphf_threshold 1e-10]
 *
 *  --fd       central finite difference of the analytic gradient (QMMethod::getGradient, the fully assembled
 *             surface, Eh/Angstrom by contract, converted) with the step h in Bohr; 2- or 4-point stencil
 *  --pieces   piecewise checks of the integral routines with a FIXED density: the J - K/2 digest against the
 *             stored-tensor J/K, S^x and (T+V)^x against finite differences of the 1e matrices, the derivative
 *             Fock matrices against finite differences of J - K/2, and the skeleton second derivatives against
 *             finite differences of qmint::gradientOneElectron / gradientTwoElectron
 *  --newton   Newton-Raphson minimisation with the analytic gradient and Hessian (translations and rotations
 *             projected out) before the Hessian is reported; the run prints |g| per step
 *
 * Always reported: Hessian symmetry (before symmetrisation), translational invariance (sum over atoms of every
 * row block) and rotational invariance H (L r) = L g (L the generator of an infinitesimal rotation), which holds
 * at every geometry, the six smallest |eigenvalues|, timings.
 *
 * Claude Generated (Nov 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/qm_hessian.h"
#include "src/core/energy_calculators/qm_methods/qm_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace {

std::vector<GTO::Orbital> shiftedBasis(const std::vector<GTO::Orbital>& b, int atom, int k, double h)
{
    std::vector<GTO::Orbital> out = b;
    for (auto& o : out)
        if (o.atom == atom) {
            if (k == 0) o.x += h;
            if (k == 1) o.y += h;
            if (k == 2) o.z += h;
        }
    return out;
}

double maxAbs(const Matrix& m) { return m.size() ? m.cwiseAbs().maxCoeff() : 0.0; }

std::vector<double> flat(const Matrix& m)
{
    std::vector<double> v;
    for (int i = 0; i < m.rows(); ++i)
        for (int j = 0; j < m.cols(); ++j) v.push_back(m(i, j));
    return v;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: dump_qm_hessian <xyz> [--basis NAME] [--fd] [--h 5e-4] [--stencil 2|4] [--pieces] [--newton]\n";
        return 1;
    }
    std::string basis = "def2-SVP";
    bool fd = false, pieces = false, newton = false;
    int threads = 0, stencil = 2;
    double h = 5.0e-4;
    qmhess::HessianOptions hopt;
    json qm_extra = json::object();
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--basis" && i + 1 < argc) basis = argv[++i];
        else if (a == "--fd") fd = true;
        else if (a == "--pieces") pieces = true;
        else if (a == "--newton") newton = true;
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
        else if (a == "--stencil" && i + 1 < argc) stencil = std::stoi(argv[++i]);
        else if (a == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
        else if (a == "--cphf_threshold" && i + 1 < argc) hopt.cphf_threshold = std::stod(argv[++i]);
        else if (a == "--qm" && i + 2 < argc) {
            const std::string key = argv[++i], val = argv[++i];
            const json parsed = json::parse(val, nullptr, false);
            qm_extra[key] = parsed.is_discarded() ? json(val) : parsed;
        }
    }
    CurcumaLogger::set_verbosity(0);

    std::map<std::string, int> Z{ { "H", 1 }, { "He", 2 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 }, { "C", 6 },
        { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Ne", 10 }, { "Na", 11 }, { "Mg", 12 }, { "Al", 13 }, { "Si", 14 },
        { "P", 15 }, { "S", 16 }, { "Cl", 17 }, { "Ar", 18 } };
    std::ifstream f(argv[1]);
    int n = 0;
    if (!(f >> n) || n <= 0) {
        std::cerr << "dump_qm_hessian: cannot read " << argv[1] << "\n";
        return 1;
    }
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    Mol mol;
    mol.m_number_atoms = n;
    mol.m_atoms.resize(n);
    mol.m_geometry = Geometry::Zero(n, 3);
    mol.m_charge = 0;
    mol.m_spin = 0;
    for (int i = 0; i < n; ++i) {
        std::string s;
        f >> s >> mol.m_geometry(i, 0) >> mol.m_geometry(i, 1) >> mol.m_geometry(i, 2);
        if (!Z.count(s)) {
            std::cerr << "dump_qm_hessian: element " << s << " not in H-Ar\n";
            return 1;
        }
        mol.m_atoms[i] = Z[s];
    }
    // tight SCF: the Hessian and the FD of the gradient are first order in the orbital error
    json cfg = { { "qm", { { "basis", basis }, { "scf_threshold", 1e-10 }, { "scf_energy_threshold", 1e-12 },
                            { "scf_gradient_threshold_tight", 1e-10 }, { "scf_max_iterations", 300 } } } };
    if (threads > 0) cfg["qm"]["threads"] = threads;
    for (auto it = qm_extra.begin(); it != qm_extra.end(); ++it) cfg["qm"][it.key()] = it.value();
    const double B2A = CurcumaUnit::Length::BOHR_TO_ANGSTROM;
    auto gradientAt = [&](const Geometry& geom, double& e) -> Matrix {
        QMMethod m(QMFunctional::HF, cfg);
        Mol x = mol;
        x.m_geometry = geom;
        m.setMolecule(x);
        e = m.calculateEnergy(true);
        if (m.hasError()) throw std::runtime_error(m.getErrorMessage());
        return m.getGradient() * B2A;  // Eh/Angstrom -> Eh/Bohr
    };
    using clk = std::chrono::steady_clock;
    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };

    json out;
    out["natoms"] = n;
    out["basis"] = basis;
    try {
        QMMethod m(QMFunctional::HF, cfg);
        m.setMolecule(mol);
        double e = m.calculateEnergy(true);
        if (m.hasError()) throw std::runtime_error(m.getErrorMessage());
        Matrix g = m.getGradient() * B2A;
        Matrix Hb;
        qmhess::HessianInfo info;

        auto translRot = [&](const Geometry& geomA) {
            // 3N x 6: translations, rotations about the origin (Bohr coordinates)
            Matrix TR = Matrix::Zero(3 * n, 6);
            for (int i = 0; i < n; ++i) {
                const Eigen::Vector3d r = geomA.row(i).transpose() / B2A;
                for (int k = 0; k < 3; ++k) {
                    TR(3 * i + k, k) = 1.0;
                    Eigen::Vector3d ek = Eigen::Vector3d::Zero();
                    ek(k) = 1.0;
                    const Eigen::Vector3d lr = ek.cross(r);
                    for (int c = 0; c < 3; ++c) TR(3 * i + c, 3 + k) = lr(c);
                }
            }
            return TR;
        };

        if (newton) {
            std::vector<double> gnorm;
            Geometry geom = mol.m_geometry;
            for (int it = 0; it < 12; ++it) {
                QMMethod mm(QMFunctional::HF, cfg);
                Mol x = mol;
                x.m_geometry = geom;
                mm.setMolecule(x);
                mm.calculateEnergy(true);
                if (mm.hasError()) throw std::runtime_error(mm.getErrorMessage());
                const Matrix gg = mm.getGradient() * B2A;
                gnorm.push_back(gg.norm());
                if (gg.norm() < 1e-9) break;
                Matrix HH;
                qmhess::HessianInfo ii;
                if (!mm.engine()->computeHessian(hopt, HH, ii)) throw std::runtime_error(ii.error);
                // Newton step in the complement of translations and rotations (rank 5 for a linear molecule)
                Matrix TR = translRot(geom);
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(TR, Eigen::ComputeFullU);
                int rank = 0;
                for (int c = 0; c < 6; ++c) rank += svd.singularValues()(c) > 1e-8 * svd.singularValues()(0) ? 1 : 0;
                const Eigen::MatrixXd Vc = svd.matrixU().rightCols(3 * n - rank);
                Eigen::VectorXd gv(3 * n);
                for (int i = 0; i < n; ++i)
                    for (int k = 0; k < 3; ++k) gv(3 * i + k) = gg(i, k);
                const Eigen::MatrixXd Hr = Vc.transpose() * HH * Vc;
                const Eigen::VectorXd step = -Vc * Hr.ldlt().solve(Vc.transpose() * gv);
                for (int i = 0; i < n; ++i)
                    for (int k = 0; k < 3; ++k) geom(i, k) += step(3 * i + k) * B2A;
            }
            out["newton_gnorm"] = gnorm;
            mol.m_geometry = geom;
            std::vector<double> gv;
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; ++k) gv.push_back(geom(i, k));
            out["newton_geometry_angstrom"] = gv;
            m.setMolecule(mol);
            e = m.calculateEnergy(true);
            g = m.getGradient() * B2A;
        }

        const auto th0 = clk::now();
        if (!m.engine()->computeHessian(hopt, Hb, info)) throw std::runtime_error(info.error);
        const double t_hess = ms(th0, clk::now());
        out["energy"] = e;
        out["gradient"] = flat(g);
        out["hessian"] = flat(Hb);
        out["hessian_ms"] = t_hess;
        out["info"] = { { "cphf_iterations", info.cphf_iterations }, { "cphf_residual", info.cphf_residual },
            { "max_asymmetry", info.max_asymmetry }, { "jk_builds", info.jk_builds }, { "t_1e_ms", info.t_1e_ms },
            { "t_2e_ms", info.t_2e_ms }, { "t_cphf_ms", info.t_cphf_ms }, { "nbf", m.engine()->nbf() } };

        // invariances
        const Matrix TR = translRot(mol.m_geometry);
        Eigen::VectorXd gv(3 * n);
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) gv(3 * i + k) = g(i, k);
        double trans = 0.0, rot = 0.0;
        for (int c = 0; c < 3; ++c) trans = std::max(trans, (Hb * TR.col(c)).cwiseAbs().maxCoeff());
        for (int c = 0; c < 3; ++c) {
            // E(R r) = E(r) => H (L r) = -L^T g = L g per atom (L the cross-product generator)
            Eigen::VectorXd Lg(3 * n);
            Eigen::Vector3d ek = Eigen::Vector3d::Zero();
            ek(c) = 1.0;
            for (int i = 0; i < n; ++i) {
                const Eigen::Vector3d gi = gv.segment<3>(3 * i);
                Lg.segment<3>(3 * i) = ek.cross(gi);
            }
            rot = std::max(rot, (Hb * TR.col(3 + c) - Lg).cwiseAbs().maxCoeff());
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hb);
        std::vector<double> ev(es.eigenvalues().data(), es.eigenvalues().data() + 3 * n);
        std::vector<double> absev;
        for (double x : ev) absev.push_back(std::abs(x));
        std::sort(absev.begin(), absev.end());
        out["translation_residual"] = trans;
        out["rotation_residual"] = rot;
        out["eigenvalues"] = ev;
        out["smallest_abs_eigenvalues"] = std::vector<double>(absev.begin(), absev.begin() + std::min<size_t>(6, absev.size()));
        out["gradient_norm"] = g.norm();
        {
            std::vector<double> gv2;
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; ++k) gv2.push_back(mol.m_geometry(i, k));
            out["geometry_angstrom"] = gv2;
        }

        if (pieces) {
            QMEngine* E = m.engine();
            const std::vector<GTO::Orbital>& bas = E->gtoBasis();
            const Matrix& Q = E->sphericalTransform();
            const Matrix P = E->density();
            const Matrix C = E->MolecularOrbitals();
            const Vector eps = E->Energies();
            const int nocc = E->numElectrons() / 2;
            Matrix W = Matrix::Zero(P.rows(), P.cols());
            for (int i = 0; i < nocc; ++i) W += 2.0 * eps(i) * (C.col(i) * C.col(i).transpose());
            const Matrix Pc = Q.size() ? Matrix(Q * P * Q.transpose()) : P;
            const Matrix Wc = Q.size() ? Matrix(Q * W * Q.transpose()) : W;
            const Matrix pos = E->atomPositionsBohr();
            const std::vector<int>& Zc = E->coreCharges();
            json pj;
            // (1) digest of the derivative driver with undifferentiated integrals vs stored-tensor J - K/2
            {
                const qmint::ERITensor eri = qmint::buildERI(bas, 1, 0.0);
                const Matrix G0 = qmint::buildCoulomb(eri, Pc) - 0.5 * qmint::buildExchange(eri, Pc);
                const Matrix G1 = qmhess::fockTwoElectronReference(bas, Pc, threads > 0 ? threads : 1);
                pj["digest_vs_tensor_JK"] = maxAbs(G1 - G0);
            }
            std::vector<Matrix> S1, H1, G1;
            Matrix h1e, h2e;
            qmhess::oneElectronDerivatives(bas, Zc, pos, Pc, Wc, S1, H1, h1e, threads > 0 ? threads : 1);
            qmhess::twoElectronDerivatives(bas, Pc, n, G1, h2e, threads > 0 ? threads : 1);
            const double hp = 1e-4;
            double dS = 0, dH = 0, dG = 0, d1 = 0, d2 = 0, mS = 0, mH = 0, mG = 0, m1 = 0, m2 = 0;
            for (int A = 0; A < n; ++A)
                for (int k = 0; k < 3; ++k) {
                    const int x = 3 * A + k;
                    Matrix posp = pos, posm = pos;
                    posp(A, k) += hp;
                    posm(A, k) -= hp;
                    const auto bp = shiftedBasis(bas, A, k, hp), bm = shiftedBasis(bas, A, k, -hp);
                    const Matrix Sfd = (qmint::buildOverlap(bp) - qmint::buildOverlap(bm)) / (2 * hp);
                    const Matrix Hfd = (qmint::buildCoreHamiltonian(bp, Zc, posp) - qmint::buildCoreHamiltonian(bm, Zc, posm)) / (2 * hp);
                    const qmint::ERITensor ep = qmint::buildERI(bp, 1, 0.0), em = qmint::buildERI(bm, 1, 0.0);
                    const Matrix Gfd = ((qmint::buildCoulomb(ep, Pc) - 0.5 * qmint::buildExchange(ep, Pc))
                                        - (qmint::buildCoulomb(em, Pc) - 0.5 * qmint::buildExchange(em, Pc))) / (2 * hp);
                    dS = std::max(dS, maxAbs(S1[x] - Sfd));
                    dH = std::max(dH, maxAbs(H1[x] - Hfd));
                    dG = std::max(dG, maxAbs(G1[x] - Gfd));
                    mS = std::max(mS, maxAbs(Sfd));
                    mH = std::max(mH, maxAbs(Hfd));
                    mG = std::max(mG, maxAbs(Gfd));
                    // skeleton second derivatives: d/dx of the fixed-density gradient contractions
                    const Matrix g1p = qmint::gradientOneElectron(bp, Zc, posp, Pc, Wc), g1m = qmint::gradientOneElectron(bm, Zc, posm, Pc, Wc);
                    const Matrix g2p = qmint::gradientTwoElectron(bp, Pc, n), g2m = qmint::gradientTwoElectron(bm, Pc, n);
                    for (int B = 0; B < n; ++B)
                        for (int l = 0; l < 3; ++l) {
                            const double f1 = (g1p(B, l) - g1m(B, l)) / (2 * hp), f2 = (g2p(B, l) - g2m(B, l)) / (2 * hp);
                            d1 = std::max(d1, std::abs(h1e(x, 3 * B + l) - f1));
                            d2 = std::max(d2, std::abs(h2e(x, 3 * B + l) - f2));
                            m1 = std::max(m1, std::abs(f1));
                            m2 = std::max(m2, std::abs(f2));
                        }
                }
            pj["S1_vs_fd"] = dS;
            pj["S1_max"] = mS;
            pj["H1_vs_fd"] = dH;
            pj["H1_max"] = mH;
            pj["G1_vs_fd"] = dG;
            pj["G1_max"] = mG;
            pj["hess1e_vs_fd"] = d1;
            pj["hess1e_max"] = m1;
            pj["hess2e_vs_fd"] = d2;
            pj["hess2e_max"] = m2;
            pj["fd_step_bohr"] = hp;
            out["pieces"] = pj;
        }

        if (fd) {
            const auto tf0 = clk::now();
            const double hA = h * B2A;
            Matrix Hfd = Matrix::Zero(3 * n, 3 * n);
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; ++k) {
                    auto gAt = [&](double s) {
                        Geometry gg = mol.m_geometry;
                        gg(i, k) += s * hA;
                        double ee;
                        return gradientAt(gg, ee);
                    };
                    Matrix col;
                    if (stencil == 4)
                        col = (-gAt(2) + 8.0 * gAt(1) - 8.0 * gAt(-1) + gAt(-2)) / (12.0 * h);
                    else
                        col = (gAt(1) - gAt(-1)) / (2.0 * h);
                    for (int j = 0; j < n; ++j)
                        for (int l = 0; l < 3; ++l) Hfd(3 * i + k, 3 * j + l) = col(j, l);
                }
            out["fd_hessian"] = flat(Hfd);
            out["fd_ms"] = ms(tf0, clk::now());
            out["fd_step_bohr"] = h;
            out["fd_stencil"] = stencil;
            out["max_abs_diff_fd"] = maxAbs(Hb - Hfd);
            out["max_abs_diff_fd_sym"] = maxAbs(Hb - 0.5 * (Hfd + Hfd.transpose()));
            out["max_abs_hessian"] = maxAbs(Hb);
        }
    } catch (const std::exception& ex) {
        std::cerr << "dump_qm_hessian: " << ex.what() << "\n";
        return 2;
    }
    std::cout << out.dump(1) << std::endl;
    return 0;
}
