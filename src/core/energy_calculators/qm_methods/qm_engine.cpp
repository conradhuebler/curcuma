/*
 * <Native ab-initio QM Engine Implementation (HF / KS-DFT)>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * WP1: real 1-electron GTO integrals (overlap S, kinetic T, nuclear attraction V)
 * over a contracted cartesian Gaussian basis (def2-SVP, H-Ne). InitialiseMolecule
 * resolves and parses the basis file, builds the flat contracted basis (Bohr
 * centers, pre-normalized), and assembles S/T/V/Hc = T+V (with a cartesian 6d ->
 * spherical 5d transform when cartesian_d is false, matching ORCA def2-SVP).
 * Calculation() still returns the nuclear repulsion energy only -- the SCF loop
 * and XC functional arrive in WP3-WP7. MakeOverlap/MakeH return the cached 1e
 * matrices for the (future) SCF driver.
 *
 * Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90
 * Second reference: ORCA 6.1 (def2-SVP).
 *
 * Claude Generated: WP1 1-electron GTO integrals
 *
 * This program is free software under GPL-3.0
 */

#include "qm_engine.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fmt/format.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace CurcumaUnit;

namespace {
// Lowercase display name for a functional (used for logging / method name).
std::string functionalName(QMFunctional f)
{
    switch (f) {
        case QMFunctional::HF:    return "hf";
        case QMFunctional::LDA:   return "lda";
        case QMFunctional::PBE:   return "pbe";
        case QMFunctional::B3LYP: return "b3lyp";
    }
    return "qm";
}

// Element number -> standard symbol. WP1 supports H (1) through Ne (10) only,
// the elements present in the shipped def2-SVP basis file.
const char* elementSymbol(int Z)
{
    static const char* kSymbols[] = {
        nullptr, "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"
    };
    if (Z < 1 || Z > 10) return nullptr;
    return kSymbols[Z];
}

bool fileExists(const std::string& path)
{
    std::ifstream f(path);
    return f.good();
}

// Resolve a basis-set name (e.g. "def2-SVP") to a .dat file path.
// Order: explicit path/.dat as-is; $CURCUMA_QM_BASIS/<name>.dat;
// $CURCUMA_DATA/<name>.dat; source-relative
// <CURCUMA_SOURCE_DIR>/src/core/energy_calculators/qm_methods/<name>.dat;
// cwd-relative <name>.dat. Returns the first existing candidate, else the
// source-relative path (parseBasisSetFile will throw a clear error if absent).
std::string resolveBasisFile(const std::string& name)
{
    // Explicit path or already a .dat file -> use as-is if it exists.
    if (name.find('/') != std::string::npos) {
        if (fileExists(name)) return name;
    }
    std::string stem = name;
    if (stem.size() >= 4 && stem.compare(stem.size() - 4, 4, ".dat") == 0)
        stem = stem.substr(0, stem.size() - 4);

    std::vector<std::string> candidates;
    // CURCUMA_QM_BASIS (CURCUMA_DFT_BASIS is the pre-Sep-2026 name, still read)
    for (const char* var : { "CURCUMA_QM_BASIS", "CURCUMA_DFT_BASIS" }) {
        if (const char* env = std::getenv(var))
            if (env[0] != '\0') candidates.push_back(std::string(env) + "/" + stem + ".dat");
    }
    if (const char* env = std::getenv("CURCUMA_DATA")) {
        if (env[0] != '\0') candidates.push_back(std::string(env) + "/" + stem + ".dat");
    }
#ifdef CURCUMA_SOURCE_DIR
    candidates.push_back(std::string(CURCUMA_SOURCE_DIR) +
                         "/src/core/energy_calculators/qm_methods/" + stem + ".dat");
#endif
    candidates.push_back(stem + ".dat");  // cwd-relative fallback

    for (const auto& p : candidates)
        if (fileExists(p)) return p;

#ifdef CURCUMA_SOURCE_DIR
    return std::string(CURCUMA_SOURCE_DIR) +
           "/src/core/energy_calculators/qm_methods/" + stem + ".dat";
#else
    return stem + ".dat";
#endif
}
}  // namespace

// =================================================================================
// Constructor
// =================================================================================

QMEngine::QMEngine(QMFunctional functional, const json& config)
    : m_functional(functional)
{
    if (config.contains("basis") && config["basis"].is_string())
        m_basis_name = config["basis"].get<std::string>();
    if (config.contains("cartesian_d") && config["cartesian_d"].is_boolean())
        m_cartesian_d = config["cartesian_d"].get<bool>();
    // is_number(), not is_number_integer(): the CLI and the registry defaults can
    // deliver integers as floating point (1.0), and the old test then silently kept
    // QMDriver's default of 4 threads whatever -qm.threads said (Sep 2026).
    if (config.contains("threads") && config["threads"].is_number())
        m_threads = std::max(1, static_cast<int>(std::lround(config["threads"].get<double>())));
    if (config.contains("eri_screening") && config["eri_screening"].is_number())
        m_eri_screening = config["eri_screening"].get<double>();
    if (config.contains("scf_max_iterations") && config["scf_max_iterations"].is_number())
        m_scf_max_iter = static_cast<int>(std::lround(config["scf_max_iterations"].get<double>()));
    if (config.contains("scf_threshold") && config["scf_threshold"].is_number())
        m_scf_threshold = config["scf_threshold"].get<double>();
    if (config.contains("scf_mode") && config["scf_mode"].is_string())
        m_scf_mode = config["scf_mode"].get<std::string>();
    if (config.contains("scf_guess") && config["scf_guess"].is_string())
        m_scf_guess = config["scf_guess"].get<std::string>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Initializing native QM engine (functional={}, basis={})",
                                        functionalName(m_functional), m_basis_name));
        CurcumaLogger::param("functional", functionalName(m_functional));
        CurcumaLogger::param("basis", m_basis_name);
        CurcumaLogger::param("cartesian_d", m_cartesian_d ? "true" : "false");
    }
}

// =================================================================================
// QMDriver Interface
// =================================================================================

bool QMEngine::InitialiseMolecule()
{
    // Geometry/atoms/charge are already loaded by QMInterface::InitialiseMolecule(Mol).
    if (m_atoms.empty()) {
        CurcumaLogger::error("QM: no atoms in molecule for initialization");
        return false;
    }
    if (m_integrals_ready) return true;  // already built (WP1: static geometry)

    // Resolve + parse the basis file (cached).
    try {
        m_basis_file = resolveBasisFile(m_basis_name);
        if (m_basis_map_cache.empty())
            m_basis_map_cache = BasisSetParser::parseBasisSetFile(m_basis_file);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("QM: basis setup failed: ") + e.what());
        return false;
    }

    // Build the flat contracted cartesian basis (Bohr centers, pre-normalized).
    m_gto_basis.clear();
    for (int i = 0; i < m_atomcount; ++i) {
        const int Z = m_atoms[i];
        const char* sym = elementSymbol(Z);
        if (sym == nullptr) {
            CurcumaLogger::error(fmt::format(
                "QM: element Z={} is outside the WP1 H-Ne scope", Z));
            return false;
        }
        auto it = m_basis_map_cache.find(sym);
        if (it == m_basis_map_cache.end()) {
            CurcumaLogger::error(fmt::format(
                "QM: no basis data for element {} in {}", sym, m_basis_file));
            return false;
        }
        const double bx = Length::angstrom_to_bohr(m_geometry(i, 0));
        const double by = Length::angstrom_to_bohr(m_geometry(i, 1));
        const double bz = Length::angstrom_to_bohr(m_geometry(i, 2));
        auto orbs = BasisSetParser::createGTOFromBasis(it->second, bx, by, bz, i);
        m_gto_basis.insert(m_gto_basis.end(), orbs.begin(), orbs.end());
    }

    // Number of electrons (nuclear charge sum minus total molecular charge).
    int nelec = 0;
    for (int Z : m_atoms) nelec += Z;
    nelec -= static_cast<int>(std::round(m_charge));
    m_num_electrons = nelec;

    buildOneElectronIntegrals();
    m_integrals_ready = true;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "QM: basis {} -> {} basis functions ({}), {} electrons",
            m_basis_file, m_nbf, m_cartesian_d ? "cartesian 6d" : "spherical 5d",
            m_num_electrons));
    }
    return true;
}

// Claude Generated (Sep 2026): the base UpdateMolecule() only returned true, so a
// second geometry (opt step, scan point, trajectory frame) reused the first
// geometry's integrals and SCF -- m_integrals_ready was never cleared.
bool QMEngine::UpdateMolecule()
{
    m_integrals_ready = false;
    m_eri_ready = false;
    m_eri_active_ready = false;
    m_scf_ready = false;
    m_scf_converged = false;
    return InitialiseMolecule();
}

double QMEngine::Calculation(bool gradient)
{
    // Build the 1e integrals if not done yet (the SCF uses them).
    if (!m_integrals_ready) InitialiseMolecule();

    const double e_nn = calculateCoreRepulsionEnergy();

    // Only the HF functional has a real SCF (rung 666, exact exchange only). The
    // DFT functionals (LDA/PBE/B3LYP) need V_xc (WP4 grid -> WP5-WP7); they stay on
    // the scaffold (E_nn + 0) with a clear notice rather than a silently wrong HF
    // energy. Honest per CLAUDE.md.
    if (m_functional != QMFunctional::HF) {
        m_total_energy = e_nn;
        m_scf_converged = false;
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::result(fmt::format(
                "native QM ({}) -- V_xc not yet implemented (WP5-WP7); "
                "returning nuclear repulsion only", functionalName(m_functional)));
            CurcumaLogger::energy_abs(m_total_energy, "nuclear repulsion energy");
        }
        return m_total_energy;
    }

    // --- HF: real closed-shell SCF (WP3) ---
    if (!m_scf_ready) {
        const bool ok = runSCF();
        m_scf_ready = true;
        if (!ok) {
            CurcumaLogger::error(fmt::format(
                "QM HF SCF did not converge in {} iterations", m_scf_max_iter));
        }
    }
    m_total_energy = m_e_elec + e_nn;

    // WP8: analytic gradient from the converged density (energy-only callers skip it).
    if (gradient) {
        if (m_scf_converged)
            computeGradient();
        else {
            m_gradient = Matrix::Zero(m_atomcount, 3);
            CurcumaLogger::warn("QM: SCF not converged -- no gradient computed (zero returned)");
        }
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::result(fmt::format("native QM -- HF/{} SCF {} ({} iterations)",
            m_basis_name, m_scf_converged ? "converged" : "NOT converged",
            m_scf_iterations));
        if (m_scf_converged && CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("E_elec", fmt::format("{:.10f} Eh", m_e_elec));
            CurcumaLogger::param("E_nn  ", fmt::format("{:.10f} Eh", e_nn));
        }
        CurcumaLogger::energy_abs(m_total_energy, "total HF energy");
    }
    return m_total_energy;
}

std::string QMEngine::getMethodNameStr() const
{
    return functionalName(m_functional);
}

// =================================================================================
// QMDriver pure-virtual hooks
// =================================================================================

Matrix QMEngine::MakeOverlap(Basisset& basisset)
{
    (void)basisset;  // DFT uses its own m_gto_basis, not the STO Basisset typedef
    return m_S;
}

Matrix QMEngine::MakeH(const Matrix& S, const Basisset& basisset)
{
    (void)S;
    (void)basisset;
    return m_H;  // Hc = T + V (nuclear repulsion is separate)
}

// =================================================================================
// 1-electron integral assembly (WP1)
// =================================================================================

void QMEngine::buildOneElectronIntegrals()
{
    // Atom positions in Bohr (the integrals use atomic units).
    Matrix atomPosBohr = Matrix::Zero(m_atomcount, 3);
    for (int i = 0; i < m_atomcount; ++i)
        for (int k = 0; k < 3; ++k)
            atomPosBohr(i, k) = Length::angstrom_to_bohr(m_geometry(i, k));

    const Matrix Scart = qmint::buildOverlap(m_gto_basis);
    const Matrix Tcart = qmint::buildKinetic(m_gto_basis);
    const Matrix Vcart = qmint::buildNuclearAttraction(m_gto_basis, m_atoms, atomPosBohr);
    const Matrix Hcart = Tcart + Vcart;

    // New geometry -> invalidate the WP3 active-basis ERI and SCF caches.
    m_eri_active_ready = false;
    m_scf_ready = false;
    m_Q = Matrix();  // rebuilt in the spherical branch; stays empty for cartesian_d

    if (m_cartesian_d) {
        m_S = Scart;
        m_T = Tcart;
        m_V = Vcart;
        m_H = Hcart;
        m_nbf = static_cast<int>(Scart.rows());
        return;
    }

    // Spherical 5d (ORCA def2-SVP default): M_sph = Q^T M_cart Q. s/p pass
    // through; each 6d shell collapses to 5d. If the basis has no d, Q is empty
    // and the cartesian matrices are kept unchanged. Q is cached in m_Q so the
    // WP3 SCF can build the active-basis ERI (applySphericalTransformERI) without
    // recomputing it.
    const Matrix Q = qmint::buildSphericalTransform(m_gto_basis, Scart);
    m_Q = Q;
    m_eri_active_ready = false;  // geometry change invalidates the active ERI
    if (Q.size() == 0) {
        m_S = Scart;
        m_T = Tcart;
        m_V = Vcart;
        m_H = Hcart;
        m_nbf = static_cast<int>(Scart.rows());
        return;
    }
    m_S = qmint::applySphericalTransform(Scart, Q);
    m_T = qmint::applySphericalTransform(Tcart, Q);
    m_V = qmint::applySphericalTransform(Vcart, Q);
    m_H = qmint::applySphericalTransform(Hcart, Q);
    m_nbf = static_cast<int>(Q.cols());
}

// =================================================================================
// 2-electron integral assembly (WP2) -- lazy
// =================================================================================

const qmint::ERITensor& QMEngine::cartesianERI() const
{
    // Built on first request and cached. NOT called by Calculation() (the
    // scaffold -sp path returns E_nn only), so the ERI cost is paid only by the
    // WP2 dumper and (later) the WP3 SCF. Built in the cartesian basis
    // (m_gto_basis); the spherical 5d transform is the caller's job via
    // qmint::applySphericalTransformERI.
    if (!m_eri_ready) {
        if (!m_integrals_ready) {
            // InitialiseMolecule is non-const; the lazy build requires the 1e
            // basis to exist. Guard: if not ready, return the empty tensor
            // (caller must have called InitialiseMolecule first).
            return m_eri_cart;
        }
        const auto t0 = std::chrono::steady_clock::now();
        m_eri_cart = qmint::buildERI(m_gto_basis, m_threads, m_eri_screening);
        m_eri_ready = true;
        if (CurcumaLogger::get_verbosity() >= 2) {
            const int nc = (int)m_gto_basis.size();
            CurcumaLogger::info(fmt::format(
                "QM: built 4-centre ERI (cartesian, {} basis functions, {:.1f} MB) in {:.1f} ms ({} threads)",
                nc, (double)nc * nc * nc * nc * 8.0 / 1.0e6,
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(),
                m_threads));
        }
    }
    return m_eri_cart;
}

// =================================================================================
// WP8: analytic RHF nuclear gradient (Claude Generated, Sep 2026)
// =================================================================================

// dE/dR_A of E_nn = sum_{i<j} Z_i Z_j / R_ij:  -Z_A Z_B (R_A - R_B) / R_AB^3  [Eh/Bohr]
Matrix QMEngine::calculateCoreRepulsionGradient() const
{
    Matrix g = Matrix::Zero(m_atomcount, 3);
    for (int i = 0; i < m_atomcount; ++i)
        for (int j = i + 1; j < m_atomcount; ++j) {
            double d[3], r2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                d[k] = Length::angstrom_to_bohr(m_geometry(i, k) - m_geometry(j, k));
                r2 += d[k] * d[k];
            }
            if (r2 < 1.0e-24) continue;
            const double f = -(double)m_atoms[i] * m_atoms[j] / (r2 * std::sqrt(r2));
            for (int k = 0; k < 3; ++k) {
                g(i, k) += f * d[k];
                g(j, k) -= f * d[k];
            }
        }
    return g;
}

// Gradient = one-electron (T, V incl. Hellmann-Feynman, -W dS) + two-electron
// + nuclear repulsion; see the WP8 block in qm_integrals.cpp for the formulas.
// The integrals are differentiated in the CARTESIAN basis, so the active-basis
// density P and energy-weighted density W are taken back with the 6d->5d map Q:
// E = Tr(P_sph H_sph) = Tr(P_sph Q^T H_cart Q) = Tr((Q P_sph Q^T) H_cart). Q only
// mixes components of one shell on one atom, so it does not depend on geometry.
bool QMEngine::computeGradient()
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n_occ = m_num_electrons / 2;
    Matrix W = Matrix::Zero(m_nbf, m_nbf);
    for (int i = 0; i < n_occ; ++i)
        W += 2.0 * m_energies(i) * (m_mo.col(i) * m_mo.col(i).transpose());

    Matrix Pc = m_density, Wc = W;
    if (m_Q.size() != 0) {
        Pc = m_Q * m_density * m_Q.transpose();
        Wc = m_Q * W * m_Q.transpose();
    }

    Matrix posBohr = Matrix::Zero(m_atomcount, 3);
    for (int i = 0; i < m_atomcount; ++i)
        for (int k = 0; k < 3; ++k)
            posBohr(i, k) = Length::angstrom_to_bohr(m_geometry(i, k));

    const Matrix g1 = qmint::gradientOneElectron(m_gto_basis, m_atoms, posBohr, Pc, Wc, m_threads);
    const auto t1 = std::chrono::steady_clock::now();
    const Matrix g2 = qmint::gradientTwoElectron(m_gto_basis, Pc, m_atomcount, m_threads);
    const auto t2 = std::chrono::steady_clock::now();
    const Matrix gn = calculateCoreRepulsionGradient();
    m_gradient = g1 + g2 + gn;
    m_gradient_parts = { g1, g2, gn };

    if (CurcumaLogger::get_verbosity() >= 2) {
        auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        CurcumaLogger::info(fmt::format(
            "QM: analytic gradient |g| = {:.6e} Eh/Bohr (1e {:.1f} ms, 2e {:.1f} ms, {} threads)",
            m_gradient.norm(), ms(t0, t1), ms(t1, t2), m_threads));
    }
    return true;
}

// =================================================================================
// Energy components
// =================================================================================

double QMEngine::calculateCoreRepulsionEnergy() const
{
    // E_nn = sum_{i<j} Z_i * Z_j / R_ij  in atomic units (Hartree).
    // m_atoms holds nuclear charges (element numbers); m_geometry is in Angstrom.
    double e_nn = 0.0;
    const int n = m_atomcount;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double dx = m_geometry(i, 0) - m_geometry(j, 0);
            const double dy = m_geometry(i, 1) - m_geometry(j, 1);
            const double dz = m_geometry(i, 2) - m_geometry(j, 2);
            const double r_ang = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r_ang < 1.0e-12) continue;  // skip coincident atoms
            const double r_bohr = Length::angstrom_to_bohr(r_ang);
            e_nn += (m_atoms[i] * m_atoms[j]) / r_bohr;
        }
    }
    return e_nn;
}