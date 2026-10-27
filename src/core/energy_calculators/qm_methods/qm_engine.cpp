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
#include "src/core/blas_threads.h"
#include "src/core/curcuma_logger.h"
#include "qm_citations.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"
#include "src/core/units.h"

#include <fmt/format.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
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
    if (config.contains("gpu") && config["gpu"].is_string()) {
        m_gpu = config["gpu"].get<std::string>();
        std::transform(m_gpu.begin(), m_gpu.end(), m_gpu.begin(), ::tolower);
        // -gpu auto was removed on master (Oct 26, 2026, operator decision): the backend is always named.
        if (m_gpu == "auto") {
            CurcumaLogger::warn("QM: -gpu auto was removed; name the backend (-gpu cuda|rocm|vulkan) or omit -gpu "
                                "for the CPU. Using CPU.");
            curcuma::reportGpuFallback("QM: -gpu auto was removed, calculation on the CPU", m_gpu);
            m_gpu = "none";
        }
    }
    if (config.contains("gpu_eri") && config["gpu_eri"].is_boolean())
        m_gpu_eri = config["gpu_eri"].get<bool>();
    if (config.contains("gpu_direct")) {
        // Accept a CLI boolean too (-qm.gpu_direct true/false).
        if (config["gpu_direct"].is_boolean())
            m_gpu_direct = config["gpu_direct"].get<bool>() ? "on" : "off";
        else if (config["gpu_direct"].is_string())
            m_gpu_direct = config["gpu_direct"].get<std::string>();
    }
    if (config.contains("gpu_scf")) {
        // Accept a CLI boolean too (-qm.gpu_scf true/false). Claude Generated (Oct 2026, Q-G4).
        if (config["gpu_scf"].is_boolean())
            m_gpu_scf = config["gpu_scf"].get<bool>() ? "on" : "off";
        else if (config["gpu_scf"].is_string())
            m_gpu_scf = config["gpu_scf"].get<std::string>();
    }
    if (config.contains("gpu_grad")) {
        // Accept a CLI boolean too (-qm.gpu_grad true/false). Claude Generated (Oct 2026, Q-G7).
        if (config["gpu_grad"].is_boolean())
            m_gpu_grad = config["gpu_grad"].get<bool>() ? "on" : "off";
        else if (config["gpu_grad"].is_string())
            m_gpu_grad = config["gpu_grad"].get<std::string>();
    }
    if (config.contains("gpu_grad1e")) {
        // Accept a CLI boolean too (-qm.gpu_grad1e true/false). Claude Generated (Oct 2026, Q-G8).
        if (config["gpu_grad1e"].is_boolean())
            m_gpu_grad1e = config["gpu_grad1e"].get<bool>() ? "on" : "off";
        else if (config["gpu_grad1e"].is_string())
            m_gpu_grad1e = config["gpu_grad1e"].get<std::string>();
    }
    if (config.contains("gpu_scf_min_nbf") && config["gpu_scf_min_nbf"].is_number())
        m_gpu_scf_min_nbf = static_cast<int>(std::lround(config["gpu_scf_min_nbf"].get<double>()));
    if (config.contains("gpu_max_memory_mb") && config["gpu_max_memory_mb"].is_number())
        m_gpu_max_memory_mb = config["gpu_max_memory_mb"].get<double>();
    if (config.contains("gpu_device") && config["gpu_device"].is_number())
        m_gpu_device = std::max(0, static_cast<int>(std::lround(config["gpu_device"].get<double>())));
    if (config.contains("eri_screening") && config["eri_screening"].is_number())
        m_eri_screening = config["eri_screening"].get<double>();
    if (config.contains("scf_max_iterations") && config["scf_max_iterations"].is_number())
        m_scf_max_iter = static_cast<int>(std::lround(config["scf_max_iterations"].get<double>()));
    if (config.contains("scf_threshold") && config["scf_threshold"].is_number())
        m_scf_threshold = config["scf_threshold"].get<double>();
    if (config.contains("scf_convergence") && config["scf_convergence"].is_string())
        m_scf_convergence = config["scf_convergence"].get<std::string>();
    if (config.contains("scf_energy_threshold") && config["scf_energy_threshold"].is_number())
        m_scf_energy_threshold = config["scf_energy_threshold"].get<double>();
    if (config.contains("scf_gradient_threshold") && config["scf_gradient_threshold"].is_number())
        m_scf_gradient_threshold = config["scf_gradient_threshold"].get<double>();
    if (config.contains("scf_gradient_threshold_tight") && config["scf_gradient_threshold_tight"].is_number())
        m_scf_gradient_threshold_tight = config["scf_gradient_threshold_tight"].get<double>();
    if (m_scf_convergence != "auto" && m_scf_convergence != "density") {
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.scf_convergence '{}' (use auto|density); using auto", m_scf_convergence));
        m_scf_convergence = "auto";
    }
    if (config.contains("scf_warm_start") && config["scf_warm_start"].is_boolean())
        m_scf_warm_start = config["scf_warm_start"].get<bool>();
    // Oct 27 2026: scf_diis_start / scf_diis_subspace (QMMethod::engineConfig maps the old spellings -qm.diis_start /
    // -qm.diis_subspace onto them). A bare top-level "diis_start" is NOT read: it is the geometry optimiser's.
    if (config.contains("scf_diis_start") && config["scf_diis_start"].is_number())
        m_diis_start = std::max(0, static_cast<int>(std::lround(config["scf_diis_start"].get<double>())));
    if (config.contains("scf_diis_subspace") && config["scf_diis_subspace"].is_number())
        m_diis_subspace = std::max(2, static_cast<int>(std::lround(config["scf_diis_subspace"].get<double>())));
    if (config.contains("scf_mode") && config["scf_mode"].is_string()) {
        m_scf_mode = config["scf_mode"].get<std::string>();
        std::transform(m_scf_mode.begin(), m_scf_mode.end(), m_scf_mode.begin(), ::tolower);
        if (m_scf_mode == "auto")
            m_scf_mode = "diis";  // the shared scf module's "method default" marker
        // Oct 27 2026: an SCF mode of another solver (broyden, level-shift: gfn1/gfn2) is an error, as for NDDO,
        // instead of a warning at run time and silent DIIS.
        if (m_scf_mode != "diis" && m_scf_mode != "adiis" && m_scf_mode != "plain")
            throw std::invalid_argument(fmt::format(
                "scf_mode '{}' is not available for native HF ({}); use diis, adiis or plain",
                config["scf_mode"].get<std::string>(), functionalName(m_functional)));
    }
    if (config.contains("scf_damping") && config["scf_damping"].is_number()) {
        const double d = config["scf_damping"].get<double>();
        if (!(d > 0.0 && d <= 1.0))
            throw std::invalid_argument(fmt::format("scf_damping {} is outside (0, 1] for native HF", d));
        m_scf_damping = d;
    }
    if (config.contains("scf_verified_stop") && config["scf_verified_stop"].is_boolean())
        m_scf_verified_stop = config["scf_verified_stop"].get<bool>();
    if (config.contains("scf_extrapolation") && config["scf_extrapolation"].is_string()) {
        m_scf_extrapolation = config["scf_extrapolation"].get<std::string>();
        std::transform(m_scf_extrapolation.begin(), m_scf_extrapolation.end(), m_scf_extrapolation.begin(), ::tolower);
        if (m_scf_extrapolation != "none" && m_scf_extrapolation != "aspc" && m_scf_extrapolation != "gauss")
            throw std::invalid_argument(fmt::format(
                "scf_extrapolation '{}' is not available for native HF; use none, aspc or gauss", m_scf_extrapolation));
    }
    if (config.contains("scf_extrapolation_order") && config["scf_extrapolation_order"].is_number())
        m_scf_extrap_order = std::max(0, static_cast<int>(std::lround(config["scf_extrapolation_order"].get<double>())));
    if (config.contains("scf_allow_unconverged")) {
        const json& v = config["scf_allow_unconverged"];
        if (v.is_boolean())
            m_scf_allow_unconverged = v.get<bool>();
        else if (v.is_number())
            m_scf_allow_unconverged = v.get<double>() != 0.0;
        else if (v.is_string())
            m_scf_allow_unconverged = (v.get<std::string>() == "true" || v.get<std::string>() == "1");
    }
    if (config.contains("scf_guess") && config["scf_guess"].is_string())
        m_scf_guess = config["scf_guess"].get<std::string>();
    if (config.contains("scf_early_screening") && config["scf_early_screening"].is_number())
        m_scf_early_screening = config["scf_early_screening"].get<double>();
    if (config.contains("scf_early_switch") && config["scf_early_switch"].is_number())
        m_scf_early_switch = config["scf_early_switch"].get<double>();
    if (config.contains("gpu_eri_precision") && config["gpu_eri_precision"].is_string())
        m_gpu_eri_precision = config["gpu_eri_precision"].get<std::string>();
    if (config.contains("gpu_fp32_switch") && config["gpu_fp32_switch"].is_number())
        m_gpu_fp32_switch = config["gpu_fp32_switch"].get<double>();
    if (config.contains("eri_prim_density") && config["eri_prim_density"].is_boolean())
        m_eri_prim_density = config["eri_prim_density"].get<bool>();
    if (config.contains("scf_direct")) {
        // Accept a CLI boolean too (-qm.scf_direct true/false).
        if (config["scf_direct"].is_boolean())
            m_scf_direct = config["scf_direct"].get<bool>() ? "on" : "off";
        else if (config["scf_direct"].is_string())
            m_scf_direct = config["scf_direct"].get<std::string>();
    }
    if (config.contains("eri_max_memory_mb") && config["eri_max_memory_mb"].is_number())
        m_eri_max_memory_mb = config["eri_max_memory_mb"].get<double>();
    readRIConfig(config);  // WP10 -qm.ri* (qm_ri.cpp)

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("Initializing native QM engine (functional={}, basis={})",
                                        functionalName(m_functional), m_basis_name));
        CurcumaLogger::param("functional", functionalName(m_functional));
        CurcumaLogger::param("basis", m_basis_name);
        CurcumaLogger::param("cartesian_d", m_cartesian_d ? "true" : "false");
    }
}

QMEngine::~QMEngine() = default;

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
        qmcite::orbitalBasis(m_basis_name);   // literature of the basis set actually loaded
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
    m_direct = qmint::DirectJK();
    m_gpu_jk.reset();
    m_gpu_djk.reset();
    m_gpu_jk_checked = false;
    m_scf_ready = false;
    m_scf_converged = false;
    ++m_geom_serial;  // a new geometry of the same molecule (density history of -qm.scf_extrapolation)
    return InitialiseMolecule();
}

// Claude Generated (Sep 2026): a reused engine object must not keep the previous
// molecule's integrals (InitialiseMolecule() returns early once they exist) nor
// warm-start from its orbitals.
void QMEngine::resetForNewMolecule()
{
    m_integrals_ready = false;
    m_eri_ready = false;
    m_eri_active_ready = false;
    m_direct = qmint::DirectJK();
    m_gpu_jk.reset();
    m_gpu_djk.reset();
    m_gpu_jk_checked = false;
    m_scf_ready = false;
    m_scf_converged = false;
    m_gpu_scf_state.reset();
    m_warm_C = Matrix();
    m_warm_atoms.clear();
    m_dens_hist.clear();
    m_dens_hist_serial = -1;
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
    // Claude Generated (Oct 25 2026): an analytic gradient (or a correlated method, setTightOrbitals) needs the
    // tight orbital-gradient threshold; an energy-only SCF of this geometry is then continued from its density.
    const bool want_tight = gradient || m_scf_tight_orbitals;
    const bool continue_tight = m_scf_ready && m_scf_converged && want_tight && !m_scf_tight_done
        && m_scf_convergence == "auto";
    if (!m_scf_ready || continue_tight) {
        m_scf_tight_now = want_tight;
        m_scf_continue = continue_tight;
        m_scf_tight_done = false;
        const bool ok = runSCF();
        m_scf_continue = false;
        m_scf_ready = true;
        if (!ok && !riSetupFailed()) {  // WP10: a failed RI setup has reported its own error
            // Oct 27 2026: master's rule (SCF_MODES.md "Failure handling"): the wrappers turn this into an error
            // unless -scf_allow_unconverged true, which keeps the energy of the last density with a warning.
            if (m_scf_allow_unconverged)
                CurcumaLogger::warn(fmt::format(
                    "QM HF SCF did not converge in {} iterations; energy of the last density accepted "
                    "(-scf_allow_unconverged true)", m_scf_max_iter));
            else
                CurcumaLogger::error(fmt::format(
                    "QM HF SCF did not converge in {} iterations", m_scf_max_iter));
        }
        if (riActive()) riReport();
    }
    m_total_energy = m_e_elec + e_nn;

    // WP8: analytic gradient from the converged density (energy-only callers skip it).
    // With -qm.ri j|jk the two-electron part is the derivative of the RI energy (qm_ri.cpp section 6);
    // if it cannot be formed, m_gradient_failed is set and QMMethod reports an error (no number returned).
    m_gradient_failed = false;
    m_gradient_error.clear();
    if (gradient) {
        if (m_scf_converged)
            computeGradient();
        else {
            // Oct 27 2026: a failed gradient, not a zero one (a zero gradient let -opt report convergence)
            m_gradient = Matrix::Zero(m_atomcount, 3);
            m_gradient_failed = true;
            m_gradient_error = "QM: SCF not converged -- no analytic gradient";
        }
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::result(fmt::format("native QM -- HF/{}{} SCF {} ({} iterations)",
            m_basis_name, riActive() ? (m_ri_mode == "j" ? " RI-J" : " RI-JK") : "",
            m_scf_converged ? "converged" : "NOT converged",
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

    // OpenMP team for the row-parallel 1e builders (bit-identical for any count; Claude Generated, Oct 2026)
    curcuma::ScopedBlasThreads omp_guard(std::max(1, m_threads));
    const Matrix Scart = qmint::buildOverlap(m_gto_basis);
    const Matrix Tcart = qmint::buildKinetic(m_gto_basis);
    const Matrix Vcart = qmint::buildNuclearAttraction(m_gto_basis, m_atoms, atomPosBohr);
    const Matrix Hcart = Tcart + Vcart;

    // New geometry -> invalidate the WP3 active-basis ERI and SCF caches.
    resetRI();  // WP10: fitted tensor belongs to the old geometry
    m_eri_active_ready = false;
    m_direct = qmint::DirectJK();
    m_gpu_jk.reset();
    m_gpu_djk.reset();
    m_gpu_jk_checked = false;
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
    m_direct = qmint::DirectJK();
    m_gpu_jk.reset();
    m_gpu_djk.reset();
    m_gpu_jk_checked = false;
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

    // RI (Oct 2026): the two-electron part comes from riGradientTwoElectron, so the exact-path device
    // 2e gradient (Q-G7) is not asked for; the 1e part (Q-G8) is the same as without RI.
    const bool ri = riActive();
    // Q-G7/Q-G8: the device gradient parts read one set of tables (cartesian, Hermite tables raised by
    // one on both shells), built once per gradient on the host when at least one part runs there.
    const std::string be1 = gpuGradBackend(m_gpu_grad1e, "gpu_grad1e", "qm_grad1e", m_gpu_grad1e_failed);
    const std::string be2 = ri ? std::string() : gpuGradBackend(m_gpu_grad, "gpu_grad", "qm_grad", m_gpu_grad_failed);
    qmint::EriGpuTables tables;
    std::string tab_why;
    bool tables_ok = false;
    const auto t_tab0 = std::chrono::steady_clock::now();
    if (!be1.empty() || !be2.empty())
        tables_ok = qmint::buildEriGradGpuTables(m_gto_basis, m_eri_screening, tables, tab_why);
    const auto t_tab = std::chrono::steady_clock::now();

    // One-electron part (T, V incl. Hellmann-Feynman, -W dS): device (Q-G8) or CPU.
    Matrix g1;
    std::string dev1_info;
    const bool on_device1 = !be1.empty()
        && gpuGradientOneElectron(tables_ok ? &tables : nullptr, tab_why, be1, posBohr, Pc, Wc, g1, dev1_info);
    if (!on_device1) g1 = qmint::gradientOneElectron(m_gto_basis, m_atoms, posBohr, Pc, Wc, m_threads);
    const auto t1 = std::chrono::steady_clock::now();
    // Two-electron part: RI (-qm.ri j|jk), device (Q-G7) or CPU.
    Matrix g2;
    std::string dev_info;
    bool on_device = false;
    if (ri) {
        if (!riGradientTwoElectron(Pc, g2, dev_info)) {
            m_gradient = Matrix::Zero(m_atomcount, 3);
            m_gradient_failed = true;
            m_gradient_error = dev_info.empty() ? std::string("RI gradient failed") : dev_info;
            CurcumaLogger::error("QM: " + m_gradient_error);
            return false;
        }
        on_device = dev_info.find("device") != std::string::npos;
    } else {
        on_device = !be2.empty() && gpuGradientTwoElectron(tables_ok ? &tables : nullptr, tab_why, be2, Pc, g2, dev_info);
        if (!on_device) g2 = qmint::gradientTwoElectron(m_gto_basis, Pc, m_atomcount, m_threads);
    }
    const auto t2 = std::chrono::steady_clock::now();
    const Matrix gn = calculateCoreRepulsionGradient();
    m_gradient = g1 + g2 + gn;
    m_gradient_parts = { g1, g2, gn };

    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    // Diagnostic (Q-G7/Q-G8): -opt silences the method; CURCUMA_QM_SCF_TRACE=1 (the Q-G4 switch) also
    // prints one ASCII line per gradient with the place and time of each part.
    if (const char* tr = std::getenv("CURCUMA_QM_SCF_TRACE"); tr && tr[0] == '1')
        fmt::print(stderr, "QM gradient trace: 1e={} 2e={} tables_ms={:.1f} 2e_ms={:.1f} 1e_ms={:.1f}\n",
                   on_device1 ? "device" : "cpu", on_device ? "device" : "cpu", ms(t_tab0, t_tab), ms(t1, t2), ms(t_tab, t1));
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "QM: analytic gradient |g| = {:.6e} Eh/Bohr (device tables {:.1f} ms; 1e {:.1f} ms {}; 2e {:.1f} ms {}; {} threads)",
            m_gradient.norm(), ms(t_tab0, t_tab), ms(t_tab, t1), on_device1 ? "on the device: " + dev1_info : "CPU",
            ms(t1, t2), ri ? dev_info : (on_device ? "on the device: " + dev_info : "CPU"), m_threads));
    }
    return true;
}

// Claude Generated (Oct 2026, Q-G7/Q-G8): which backend a device gradient part runs on, "" = CPU.
std::string QMEngine::gpuGradBackend(const std::string& mode, const char* option, const char* tag, bool& failed) const
{
    if (!gpuRequested() || failed) return std::string();
    if (mode == "off" || mode == "false") return std::string();
    if (mode != "auto" && mode != "on" && mode != "true") {
        static bool warned[2] = { false, false };
        bool& w = warned[std::string(tag) == "qm_grad1e" ? 1 : 0];
        if (!w && CurcumaLogger::get_verbosity() >= 1) {
            w = true;
            CurcumaLogger::warn(fmt::format("QM: unknown -qm.{} '{}' (use auto|on|off); using auto", option, mode));
        }
    }
    const std::string backend = m_gpu;
    if (backend == "none" || !gpu_plugin::available(backend)) {
        failed = true;
        curcuma::reportGpuFallback(tag, "no GPU plugin for -gpu " + m_gpu);
        return std::string();
    }
    return backend;
}

// Claude Generated (Oct 2026, Q-G7): two-electron gradient on the device. Tables are rebuilt per
// gradient (they depend on the geometry; the gradient is formed once per geometry) and shared with
// the 1e part (Q-G8). The device sees the same cartesian density as qmint::gradientTwoElectron and the
// same canonical quartets; the only difference is the density-weighted quartet screen (threshold
// -qm.eri_screening / 100, 0 = all quartets, as the CPU). CURCUMA_QM_GPU_GRAD_CHECK=1 also runs the CPU
// reference and prints one comparison line (test_cases/qm_gpu/check_qm_gpu_grad.py parses it).
bool QMEngine::gpuGradientTwoElectron(const qmint::EriGpuTables* tables, const std::string& why_tables,
                                      const std::string& backend, const Matrix& Pc, Matrix& g2, std::string& info) const
{
    const auto t0 = std::chrono::steady_clock::now();
    std::string why = why_tables;
    std::unique_ptr<qmint::QMGpuGrad2e> dev;
    if (tables) {
        qmint::QMGpuGradExtra extra;
        extra.natoms = m_atomcount;
        extra.sh_atom = tables->sh_atom.data();
        extra.nprimpair = (int)tables->pp_ea.size();
        extra.pp_ea = tables->pp_ea.data();
        extra.pp_eb = tables->pp_eb.data();
        dev = gpu_plugin::createQmGrad2e(backend, tables->view(), extra, m_gpu_device, why);
    }
    const auto t1 = std::chrono::steady_clock::now();
    Matrix g = Matrix::Zero(m_atomcount, 3);
    const Matrix P = Pc;  // row-major, contiguous
    // The derivative integrals are not bounded by the Schwarz factors of the undifferentiated ones
    // (the centre derivative brings factors 2 alpha), so the quartet screen runs two orders tighter
    // than the ERI screen: benzene/def2-SVP 2.6e-10 (1e-12) vs 6.2e-12 Eh/Bohr (1e-14) from the
    // unscreened CPU gradient, for 5 % more quartets. 0 keeps every quartet (as the CPU).
    const double grad_thr = m_eri_screening * 1.0e-2;
    if (!dev || !dev->compute(P.data(), grad_thr, g.data())) {
        m_gpu_grad_failed = true;
        curcuma::reportGpuFallback("qm_grad", "device 2e gradient " + std::string(dev ? "kernel failed" : "not possible: " + why));
        return false;
    }
    g2 = g;
    const auto t2 = std::chrono::steady_clock::now();
    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    info = fmt::format("{} ({}; setup {:.1f} ms)", dev->name(), dev->stats(), ms(t0, t1));
    if (const char* chk = std::getenv("CURCUMA_QM_GPU_GRAD_CHECK"); chk && chk[0] == '1') {
        const auto c0 = std::chrono::steady_clock::now();
        const Matrix gc = qmint::gradientTwoElectron(m_gto_basis, Pc, m_atomcount, m_threads);
        const double cpu_ms = ms(c0, std::chrono::steady_clock::now());
        const Eigen::RowVector3d sd = g2.colwise().sum(), sc = gc.colwise().sum();
        CurcumaLogger::result(fmt::format(
            "QM GPU grad check: natoms={} max|dev-cpu|={:.3e} max|cpu|={:.3e} sum_dev=({:.3e},{:.3e},{:.3e}) "
            "sum_cpu=({:.3e},{:.3e},{:.3e}) quartets={} dev_ms={:.1f} cpu_ms={:.1f}",
            m_atomcount, (g2 - gc).cwiseAbs().maxCoeff(), gc.cwiseAbs().maxCoeff(), sd(0), sd(1), sd(2), sc(0), sc(1),
            sc(2), dev->lastQuartets(), ms(t0, t2), cpu_ms));
    }
    return true;
}

// Claude Generated (Oct 2026, Q-G8): one-electron gradient on the device from the Q-G7 tables plus
// shell centres, nuclei and nuclear charges. Same quantity as qmint::gradientOneElectron (no
// screening on either side). CURCUMA_QM_GPU_GRAD_CHECK=1 also runs the CPU 1e gradient from the same
// P, W and prints "QM GPU grad1e check: ..." (parsed by test_cases/qm_gpu/check_qm_gpu_grad.py).
bool QMEngine::gpuGradientOneElectron(const qmint::EriGpuTables* tables, const std::string& why_tables,
                                      const std::string& backend, const Matrix& posBohr, const Matrix& Pc,
                                      const Matrix& Wc, Matrix& g1, std::string& info) const
{
    const auto t0 = std::chrono::steady_clock::now();
    std::string why = why_tables;
    std::unique_ptr<qmint::QMGpuGrad1e> dev;
    if (tables) {
        const int ns = (int)tables->sh_atom.size();
        std::vector<double> sh_xyz(3 * (size_t)ns), atom_xyz(3 * (size_t)m_atomcount), atom_Z(m_atomcount);
        for (int s = 0; s < ns; ++s) {
            // cartesian tables: sh_out_first is the index of the shell's first function in m_gto_basis
            const GTO::Orbital& o = m_gto_basis[tables->sh_out_first[s]];
            sh_xyz[3 * s] = o.x;
            sh_xyz[3 * s + 1] = o.y;
            sh_xyz[3 * s + 2] = o.z;
        }
        for (int i = 0; i < m_atomcount; ++i) {
            for (int k = 0; k < 3; ++k) atom_xyz[3 * i + k] = posBohr(i, k);
            atom_Z[i] = (double)m_atoms[i];
        }
        qmint::QMGpuGrad1eExtra extra;
        extra.base.natoms = m_atomcount;
        extra.base.sh_atom = tables->sh_atom.data();
        extra.base.nprimpair = (int)tables->pp_ea.size();
        extra.base.pp_ea = tables->pp_ea.data();
        extra.base.pp_eb = tables->pp_eb.data();
        extra.sh_xyz = sh_xyz.data();
        extra.atom_xyz = atom_xyz.data();
        extra.atom_Z = atom_Z.data();
        dev = gpu_plugin::createQmGrad1e(backend, tables->view(), extra, m_gpu_device, why);
    }
    const auto t1 = std::chrono::steady_clock::now();
    Matrix g = Matrix::Zero(m_atomcount, 3);
    const Matrix P = Pc, W = Wc;  // row-major, contiguous
    if (!dev || !dev->compute(P.data(), W.data(), g.data())) {
        m_gpu_grad1e_failed = true;
        curcuma::reportGpuFallback("qm_grad1e", "device 1e gradient " + std::string(dev ? "kernel failed" : "not possible: " + why));
        return false;
    }
    g1 = g;
    const auto t2 = std::chrono::steady_clock::now();
    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    info = fmt::format("{} ({}; setup {:.1f} ms)", dev->name(), dev->stats(), ms(t0, t1));
    if (const char* chk = std::getenv("CURCUMA_QM_GPU_GRAD_CHECK"); chk && chk[0] == '1') {
        const auto c0 = std::chrono::steady_clock::now();
        const Matrix gc = qmint::gradientOneElectron(m_gto_basis, m_atoms, posBohr, Pc, Wc, m_threads);
        const double cpu_ms = ms(c0, std::chrono::steady_clock::now());
        const Eigen::RowVector3d sd = g1.colwise().sum(), sc = gc.colwise().sum();
        CurcumaLogger::result(fmt::format(
            "QM GPU grad1e check: natoms={} max|dev-cpu|={:.3e} max|cpu|={:.3e} sum_dev=({:.3e},{:.3e},{:.3e}) "
            "sum_cpu=({:.3e},{:.3e},{:.3e}) dev_ms={:.1f} cpu_ms={:.1f}",
            m_atomcount, (g1 - gc).cwiseAbs().maxCoeff(), gc.cwiseAbs().maxCoeff(), sd(0), sd(1), sd(2), sc(0), sc(1),
            sc(2), ms(t0, t2), cpu_ms));
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