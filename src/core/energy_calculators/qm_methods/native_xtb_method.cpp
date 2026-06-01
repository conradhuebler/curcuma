/*
 * <Native xTB Method Wrapper Implementation — unified GFN1 / GFN2>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated: replaces gfn1_method.cpp / gfn2_method.cpp with one
 * method-parametrized adapter. The body is the former GFN2Method logic with the
 * method type, name, and default config supplied at construction.
 */

#include "native_xtb_method.h"

#include "src/tools/general.h"          // MergeJson
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>
#include <algorithm>
#include <fstream>
#include <iomanip>

using curcuma::xtb::MethodType;

// ---------------------------------------------------------------------------
// Construction / configuration
// ---------------------------------------------------------------------------
NativeXtbMethod::NativeXtbMethod(MethodType method, const json& config)
    : m_method(method)
{
    m_parameters = MergeJson(getDefaultConfig(method), config);
    m_xtb = std::make_unique<curcuma::xtb::XTB>(method);
    applyConfig();

    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format("{} initialized (native xTB)", getMethodName()));
}

std::string NativeXtbMethod::getMethodName() const
{
    return m_method == MethodType::GFN1 ? "gfn1" : "gfn2";
}

json NativeXtbMethod::getDefaultConfig(MethodType method)
{
    // SCF-convergence defaults are identical for GFN1 and GFN2; only the
    // dispersion correction (and the GFN2-only D4 charge source) differ.
    json cfg = {
        { "scf_max_iterations", 100 },   // Maximum SCF iterations
        { "scf_threshold", 1.0e-5 },     // Convergence threshold on max|dq_shell|.
                                         // Energy bit-identical to 1e-6 (<1e-8 Eh) but
                                         // converges ~10-20% fewer iterations. The default
                                         // eeq D4 gradient is insensitive; MD/opt or the
                                         // opt-in mulliken response may tighten via
                                         // -scf_threshold.
        { "scf_damping", 0.4 },          // Density damping factor
        { "scf_mode", "broyden" },       // broyden(default) | diis | plain | level-shift
        { "scf_guess", "eeq" },          // Initial charge guess: eeq(default) | h0
        { "eigensolver", "mkl" },        // Eigensolve backend: mkl(default, dsyevd) | native/dnc
                                         // (self-contained Householder+QL, no LAPACK; WP4)
                                         // eeq seeds shell charges from a single-shot dftd4 EEQ
                                         // solve; on complex(231) it cuts SCF iters gfn1 35->16,
                                         // gfn2 34->22 (energy bit-identical). Falls back to h0
                                         // if the EEQ solve fails. Claude Generated 2026-06.
        { "diis_start", 5 },             // Damped warmup iterations before DIIS
        { "diis_subspace", 6 },          // DIIS history depth
        { "level_shift", 0.2 },          // Virtual-orbital shift (Eh), level-shift mode
        { "threads", 1 },                // Single-threaded by default
        { "print_orbitals", false },     // Print energy decomposition at verbosity >= 2
        { "warm_start", true },          // Reuse converged charges as SCF guess (harmless for SP)
        { "keep_diis", false }           // Preserve DIIS/Broyden history across geometry steps
    };
    if (method == MethodType::GFN1) {
        cfg["dispersion"] = "d3";        // GFN1: D3(BJ)
    } else {
        cfg["dispersion"] = "d4";        // GFN2: D4
        cfg["d4_charge_source"] = "eeq"; // D4 zeta charges: "eeq" | "mulliken"
        cfg["save_orbitals"] = false;
    }
    return cfg;
}

void NativeXtbMethod::applyConfig()
{
    if (!m_xtb) return;

    // D4 charge-response source ("eeq" default, or "mulliken" via CPSCF). The CLI
    // flag -d4_charge_source auto-routes to the "xtb" scope; fall back to a
    // top-level key, then the default. No-op for GFN1 (D3 ignores it).
    std::string d4src = "eeq";
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object()
        && m_parameters["xtb"].contains("d4_charge_source"))
        d4src = m_parameters["xtb"]["d4_charge_source"].get<std::string>();
    else
        d4src = m_parameters.value("d4_charge_source", std::string("eeq"));
    m_xtb->setD4ChargeSource(d4src);

    // SCF-convergence settings (mode, guess, damping, DIIS, level shift) from the
    // controller, reading the "xtb" scope with a top-level fallback.
    curcuma::xtb::applyXtbScfConfig(*m_xtb, m_parameters);
}

void NativeXtbMethod::handleError(const std::string& operation)
{
    m_has_error = true;
    m_error_message = fmt::format("{} error during {}", getMethodName(), operation);
}

bool NativeXtbMethod::isMoleculeSupported(const Mol& mol)
{
    // Native GFN1/GFN2 parameters cover Z = 1..86 (H to Rn).
    for (int z : mol.m_atoms)
        if (z < 1 || z > 86)
            return false;
    return true;
}

// ---------------------------------------------------------------------------
// Core molecular setup and calculation
// ---------------------------------------------------------------------------
bool NativeXtbMethod::setMolecule(const Mol& mol)
{
    try {
        if (!isMoleculeSupported(mol)) {
            m_has_error = true;
            m_error_message = fmt::format("Molecule contains elements unsupported by {} (Z must be 1..86)",
                                          getMethodName());
            return false;
        }

        m_molecule = mol;
        m_calculation_done = false;
        m_initialized = true;
        clearError();

        if (!m_xtb->QMInterface::InitialiseMolecule(mol)) {
            handleError("molecule initialization");
            return false;
        }
        return true;

    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("{} setMolecule failed: {}", getMethodName(), e.what());
        return false;
    }
}

bool NativeXtbMethod::updateGeometry(const Matrix& geometry)
{
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return false;
    }

    try {
        m_molecule.m_geometry = geometry;
        m_calculation_done = false;
        clearError();

        if (!m_xtb->UpdateMolecule(geometry)) {
            handleError("geometry update");
            return false;
        }
        return true;

    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("{} updateGeometry failed: {}", getMethodName(), e.what());
        return false;
    }
}

double NativeXtbMethod::calculateEnergy(bool gradient)
{
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return 0.0;
    }

    try {
        clearError();

        // Apply controller settings each call so geometry steps in -opt / MD pick
        // up the configuration (D4 charge source + SCF convergence settings).
        applyConfig();

        m_last_energy = m_xtb->Calculation(gradient);
        m_calculation_done = true;

        if (m_parameters.value("print_orbitals", false)
            && CurcumaLogger::get_verbosity() >= 2) {
            json decomp = m_xtb->getEnergyDecomposition();
            CurcumaLogger::info(fmt::format("{} Energy Decomposition:", getMethodName()));
            CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", decomp["electronic"].get<double>()));
            CurcumaLogger::param("repulsion",  fmt::format("{:.6f} Eh", decomp["repulsion"].get<double>()));
            CurcumaLogger::param("coulomb",    fmt::format("{:.6f} Eh", decomp["coulomb"].get<double>()));
            CurcumaLogger::param("dispersion", fmt::format("{:.6f} Eh", decomp["dispersion"].get<double>()));
        }

        return m_last_energy;

    } catch (const std::exception& e) {
        handleError(fmt::format("energy calculation: {}", e.what()));
        return 0.0;
    }
}

// ---------------------------------------------------------------------------
// Warm-start / iterative-mode controls (Claude Generated)
// ---------------------------------------------------------------------------
void NativeXtbMethod::setWarmStart(bool on)
{
    if (m_xtb) m_xtb->setWarmStart(on);
}

void NativeXtbMethod::setIterativeMode(bool on)
{
    if (m_xtb) m_xtb->setIterativeMode(on);
}

// ---------------------------------------------------------------------------
// Property access
// ---------------------------------------------------------------------------
Matrix NativeXtbMethod::getGradient() const
{
    if (!m_calculation_done || !m_xtb)
        return Matrix::Zero(m_molecule.AtomCount(), 3);

    Matrix grad = m_xtb->Gradient();
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::param("native_xtb_gradient_norm", fmt::format("{:.6e} Eh/Å", grad.norm()));
    return grad;
}

Vector NativeXtbMethod::getCharges() const
{
    if (!m_calculation_done || !m_xtb)
        return Vector::Zero(m_molecule.AtomCount());
    return m_xtb->getPartialCharges();
}

void NativeXtbMethod::setParameters(const json& params)
{
    m_parameters = MergeJson(m_parameters, params);
    applyConfig();
}

// ---------------------------------------------------------------------------
// Advanced (QM) features — all delegate to the native solver
// ---------------------------------------------------------------------------
Vector NativeXtbMethod::getOrbitalEnergies() const
{
    if (!m_calculation_done || !m_xtb) return Vector{};
    try { return m_xtb->Energies(); } catch (const std::exception&) { return Vector{}; }
}

Matrix NativeXtbMethod::getMolecularOrbitals() const
{
    if (!m_calculation_done || !m_xtb) return Matrix{};
    try { return m_xtb->MolecularOrbitals(); } catch (const std::exception&) { return Matrix{}; }
}

int NativeXtbMethod::getNumElectrons() const
{
    if (!m_calculation_done || !m_xtb) return 0;
    try { return m_xtb->NumElectrons(); } catch (const std::exception&) { return 0; }
}

double NativeXtbMethod::getHOMOLUMOGap() const
{
    if (!m_calculation_done || !m_xtb) return 0.0;
    try { return m_xtb->getHOMOLUMOGap(); } catch (const std::exception&) { return 0.0; }
}

double NativeXtbMethod::getHOMOEnergy() const
{
    if (!m_calculation_done || !m_xtb) return 0.0;
    try { return m_xtb->getHOMOEnergy(); } catch (const std::exception&) { return 0.0; }
}

double NativeXtbMethod::getLUMOEnergy() const
{
    if (!m_calculation_done || !m_xtb) return 0.0;
    try { return m_xtb->getLUMOEnergy(); } catch (const std::exception&) { return 0.0; }
}

Vector NativeXtbMethod::getCoordinationNumbers() const
{
    if (!m_calculation_done || !m_xtb) return Vector{};
    try { return m_xtb->getCoordinationNumbers(); } catch (const std::exception&) { return Vector{}; }
}

json NativeXtbMethod::getEnergyDecomposition() const
{
    if (!m_calculation_done || !m_xtb) return json{};
    try { return m_xtb->getEnergyDecomposition(); } catch (const std::exception&) { return json{}; }
}

bool NativeXtbMethod::saveToFile(const std::string& filename) const
{
    if (!m_calculation_done) return false;

    try {
        json output_data;
        output_data["method"] = getMethodName();
        output_data["energy"] = m_last_energy;
        output_data["energy_decomposition"] = getEnergyDecomposition();
        output_data["orbital_energies"] = getOrbitalEnergies();
        output_data["coordination_numbers"] = getCoordinationNumbers();
        output_data["num_electrons"] = getNumElectrons();
        output_data["homo_energy"] = getHOMOEnergy();
        output_data["lumo_energy"] = getLUMOEnergy();
        output_data["homo_lumo_gap"] = getHOMOLUMOGap();
        output_data["parameters"] = m_parameters;

        std::ofstream file(filename);
        if (file.is_open()) {
            file << std::setw(2) << output_data << std::endl;
            file.close();
            return true;
        }
    } catch (const std::exception&) {
        return false;
    }
    return false;
}
