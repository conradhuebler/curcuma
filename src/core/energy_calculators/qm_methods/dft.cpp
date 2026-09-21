/*
 * <Native KS-DFT Engine Implementation>
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

#include "dft.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fmt/format.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace CurcumaUnit;

namespace {
// Lowercase display name for a functional (used for logging / method name).
std::string functionalName(DFTFunctional f)
{
    switch (f) {
        case DFTFunctional::HF:    return "hf";
        case DFTFunctional::LDA:   return "lda";
        case DFTFunctional::PBE:   return "pbe";
        case DFTFunctional::B3LYP: return "b3lyp";
    }
    return "dft";
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
// Order: explicit path/.dat as-is; $CURCUMA_DFT_BASIS/<name>.dat;
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
    if (const char* env = std::getenv("CURCUMA_DFT_BASIS")) {
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

DFT::DFT(DFTFunctional functional, const json& config)
    : m_functional(functional)
{
    if (config.contains("basis") && config["basis"].is_string())
        m_basis_name = config["basis"].get<std::string>();
    if (config.contains("cartesian_d") && config["cartesian_d"].is_boolean())
        m_cartesian_d = config["cartesian_d"].get<bool>();
    if (config.contains("threads") && config["threads"].is_number_integer())
        m_threads = config["threads"].get<int>();
    if (config.contains("scf_max_iterations") && config["scf_max_iterations"].is_number_integer())
        m_scf_max_iter = config["scf_max_iterations"].get<int>();
    if (config.contains("scf_threshold") && config["scf_threshold"].is_number())
        m_scf_threshold = config["scf_threshold"].get<double>();
    if (config.contains("scf_mode") && config["scf_mode"].is_string())
        m_scf_mode = config["scf_mode"].get<std::string>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Initializing native DFT engine (functional={}, basis={})",
                                        functionalName(m_functional), m_basis_name));
        CurcumaLogger::param("functional", functionalName(m_functional));
        CurcumaLogger::param("basis", m_basis_name);
        CurcumaLogger::param("cartesian_d", m_cartesian_d ? "true" : "false");
    }
}

// =================================================================================
// QMDriver Interface
// =================================================================================

bool DFT::InitialiseMolecule()
{
    // Geometry/atoms/charge are already loaded by QMInterface::InitialiseMolecule(Mol).
    if (m_atoms.empty()) {
        CurcumaLogger::error("DFT: no atoms in molecule for initialization");
        return false;
    }
    if (m_integrals_ready) return true;  // already built (WP1: static geometry)

    // Resolve + parse the basis file (cached).
    try {
        m_basis_file = resolveBasisFile(m_basis_name);
        if (m_basis_map_cache.empty())
            m_basis_map_cache = BasisSetParser::parseBasisSetFile(m_basis_file);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("DFT: basis setup failed: ") + e.what());
        return false;
    }

    // Build the flat contracted cartesian basis (Bohr centers, pre-normalized).
    m_gto_basis.clear();
    for (int i = 0; i < m_atomcount; ++i) {
        const int Z = m_atoms[i];
        const char* sym = elementSymbol(Z);
        if (sym == nullptr) {
            CurcumaLogger::error(fmt::format(
                "DFT: element Z={} is outside the WP1 H-Ne scope", Z));
            return false;
        }
        auto it = m_basis_map_cache.find(sym);
        if (it == m_basis_map_cache.end()) {
            CurcumaLogger::error(fmt::format(
                "DFT: no basis data for element {} in {}", sym, m_basis_file));
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
            "DFT: basis {} -> {} basis functions ({}), {} electrons",
            m_basis_file, m_nbf, m_cartesian_d ? "cartesian 6d" : "spherical 5d",
            m_num_electrons));
    }
    return true;
}

double DFT::Calculation(bool gradient)
{
    (void)gradient;  // WP8 brings the analytic gradient

    // Build the 1e integrals if not done yet (the SCF uses them).
    if (!m_integrals_ready) InitialiseMolecule();

    const double e_nn = calculateCoreRepulsionEnergy();

    // Only the HF functional has a real SCF (rung 666, exact exchange only). The
    // DFT functionals (LDA/PBE/B3LYP) need V_xc (WP4 grid -> WP5-WP7); they stay on
    // the scaffold (E_nn + 0) with a clear notice rather than a silently wrong HF
    // energy. Honest per CLAUDE.md.
    if (m_functional != DFTFunctional::HF) {
        m_total_energy = e_nn;
        m_scf_converged = false;
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::result(fmt::format(
                "native DFT ({}) -- V_xc not yet implemented (WP5-WP7); "
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
                "DFT HF SCF did not converge in {} iterations", m_scf_max_iter));
        }
    }
    m_total_energy = m_e_elec + e_nn;

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::result(fmt::format("native DFT -- HF/{} SCF {} ({} iterations)",
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

std::string DFT::getMethodNameStr() const
{
    return functionalName(m_functional);
}

// =================================================================================
// QMDriver pure-virtual hooks
// =================================================================================

Matrix DFT::MakeOverlap(Basisset& basisset)
{
    (void)basisset;  // DFT uses its own m_gto_basis, not the STO Basisset typedef
    return m_S;
}

Matrix DFT::MakeH(const Matrix& S, const Basisset& basisset)
{
    (void)S;
    (void)basisset;
    return m_H;  // Hc = T + V (nuclear repulsion is separate)
}

// =================================================================================
// 1-electron integral assembly (WP1)
// =================================================================================

void DFT::buildOneElectronIntegrals()
{
    // Atom positions in Bohr (the integrals use atomic units).
    Matrix atomPosBohr = Matrix::Zero(m_atomcount, 3);
    for (int i = 0; i < m_atomcount; ++i)
        for (int k = 0; k < 3; ++k)
            atomPosBohr(i, k) = Length::angstrom_to_bohr(m_geometry(i, k));

    const Matrix Scart = dft1e::buildOverlap(m_gto_basis);
    const Matrix Tcart = dft1e::buildKinetic(m_gto_basis);
    const Matrix Vcart = dft1e::buildNuclearAttraction(m_gto_basis, m_atoms, atomPosBohr);
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
    const Matrix Q = dft1e::buildSphericalTransform(m_gto_basis, Scart);
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
    m_S = dft1e::applySphericalTransform(Scart, Q);
    m_T = dft1e::applySphericalTransform(Tcart, Q);
    m_V = dft1e::applySphericalTransform(Vcart, Q);
    m_H = dft1e::applySphericalTransform(Hcart, Q);
    m_nbf = static_cast<int>(Q.cols());
}

// =================================================================================
// 2-electron integral assembly (WP2) -- lazy
// =================================================================================

const dft1e::ERITensor& DFT::cartesianERI() const
{
    // Built on first request and cached. NOT called by Calculation() (the
    // scaffold -sp path returns E_nn only), so the ERI cost is paid only by the
    // WP2 dumper and (later) the WP3 SCF. Built in the cartesian basis
    // (m_gto_basis); the spherical 5d transform is the caller's job via
    // dft1e::applySphericalTransformERI.
    if (!m_eri_ready) {
        if (!m_integrals_ready) {
            // InitialiseMolecule is non-const; the lazy build requires the 1e
            // basis to exist. Guard: if not ready, return the empty tensor
            // (caller must have called InitialiseMolecule first).
            return m_eri_cart;
        }
        m_eri_cart = dft1e::buildERI(m_gto_basis);
        m_eri_ready = true;
        if (CurcumaLogger::get_verbosity() >= 2) {
            const int nc = (int)m_gto_basis.size();
            CurcumaLogger::info(fmt::format(
                "DFT: built 4-centre ERI (cartesian, {} basis functions, {} entries)",
                nc, (size_t)nc * nc * nc * nc));
        }
    }
    return m_eri_cart;
}

// =================================================================================
// Energy components
// =================================================================================

double DFT::calculateCoreRepulsionEnergy() const
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