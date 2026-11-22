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
#include "src/core/solvation/ddcosmo.h"
#include "src/core/citation_registry.h"
#include "src/core/energy_calculators/ff_methods/alpb_parameters.h"

#include <cstdio>
#include <cstdlib>

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
    // CPCM-X (Claude Generated, Nov 2026): a second solver for the conductor COSMO SCF; m_xtb stays in the gas phase.
    m_cpcmx = (curcuma::xtb::solventModelCode(scopedString("solvent_model", "none")) == 5);
    if (m_cpcmx) {
        m_xtb_conductor = std::make_unique<curcuma::xtb::XTB>(method);
        m_cpcmx_solvent = curcuma::solvation::cpcmx::solventName(scopedString("solvent", "none"));
        if (method == MethodType::GFN1)   // same warning as xtb prog/main.F90
            CurcumaLogger::warn("CPCM-X was parametrized for GFN2-xTB. The results are probably inaccurate with other methods.");
        CitationRegistry::cite("cpcm_x", getMethodName());
    }
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
    // The SCF iteration limit is not set here: XTB keeps its own default (150), changed
    // only by scf_max_iterations (module scf) or -xtb.max_iterations (applyXtbScfConfig).
    json cfg = {
        { "scf_threshold", 2.0e-5 },    // Convergence threshold on max|dq_shell| (2e-5 since
                                         // Nov 2026, was 1e-5). Energy error < 0.01 uEh, rms
                                         // force error 5e-4 mEh/A on complex(231); ~6% fewer
                                         // iterations than 1e-5. The default eeq D4 gradient
                                         // is insensitive; MD/opt or the opt-in mulliken
                                         // response may tighten via -scf_threshold.
                                         // See docs/SCF_THRESHOLD_2026-11.md.
        { "scf_criterion", "default" },  // default | xtb (xtb 6.7.1 convergence rule, CPU loop)
        { "scf_accuracy", 1.0 },         // factor on the xtb limits (scf_criterion xtb)
        { "scf_damping", 0.4 },          // Density damping factor
        { "scf_mode", "broyden" },       // broyden(default) | diis | plain | level-shift
        { "scf_guess", "eeq" },          // Initial charge guess: eeq(default) | h0
        { "sto6g_legacy_4sp", false },   // false(default) = tblite's STO-6G 4s/4p tables;
                                         // true = xtb's older ones (K, Ca, Ge-Kr only).
                                         // The two references disagree here and nowhere
                                         // else; tblite's fit the exact Slater function
                                         // better. See STO_CGTO.hpp.
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
        { "keep_diis", false },          // Preserve DIIS/Broyden history across geometry steps
        { "scf_extrapolation", "none" }, // Multi-step SCC extrapolation: none|aspc|gauss (opt-in)
        { "scf_extrapolation_order", 3 },// ASPC order k (history k+2) / Gauss polynomial degree
        { "scf_extrapolation_apply", "guess" }, // guess (full SCF) | xlbomd (corrector-only, experimental)
        { "scf_xlbomd_correctors", 1 }   // Corrector SCF maps per step (xlbomd only)
    };
    if (method == MethodType::GFN1) {
        cfg["dispersion"] = "d3";        // GFN1: D3(BJ)
    } else {
        cfg["dispersion"] = "d4";        // GFN2: D4
        cfg["d4_charge_source"] = "mulliken"; // D4 q-response: "mulliken"(variational, exact) | "eeq" | "cpscf"
        cfg["save_orbitals"] = false;
    }
    return cfg;
}

void NativeXtbMethod::applyConfig()
{
    if (!m_xtb) return;

    // D4 charge-response source (default "mulliken" = variational charge-Pulay, exact;
    // "eeq" = single-shot EEQ, approximate; "cpscf" = explicit Z-vector solve). The CLI
    // flag -d4_charge_source auto-routes to the "xtb" scope; fall back to a top-level
    // key, then the default. No-op for GFN1 (D3 ignores it).
    std::string d4src = "mulliken";
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object()
        && m_parameters["xtb"].contains("d4_charge_source"))
        d4src = m_parameters["xtb"]["d4_charge_source"].get<std::string>();
    else
        d4src = m_parameters.value("d4_charge_source", std::string("mulliken"));
    m_xtb->setD4ChargeSource(d4src);

    // SCF-convergence settings (mode, guess, damping, DIIS, level shift) from the
    // controller, reading the "xtb" scope with a top-level fallback.
    curcuma::xtb::applyXtbScfConfig(*m_xtb, m_parameters);

    // CPCM-X: the conductor solver gets the same settings with solvent_model cpcmx (conductor COSMO set), the main
    // solver runs in the gas phase (xtb: cpxcalc = calc without solvation). Claude Generated (Nov 2026).
    if (m_cpcmx && m_xtb_conductor) {
        m_xtb_conductor->setD4ChargeSource(d4src);
        curcuma::xtb::applyXtbScfConfig(*m_xtb_conductor, m_parameters);
        m_xtb_conductor->setSolventModel(5);
        m_xtb->setSolvent("none");
        m_xtb->setSolventModel(0);
        m_xtb->setSolventEpsilon(-1.0);
    }
}

std::string NativeXtbMethod::scopedString(const char* key, const std::string& fallback) const
{
    const json* v = nullptr;
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object() && m_parameters["xtb"].contains(key))
        v = &m_parameters["xtb"][key];
    else if (m_parameters.contains(key))
        v = &m_parameters[key];
    if (!v) return fallback;
    if (v->is_string()) return v->get<std::string>();
    if (v->is_number()) return v->dump();
    return fallback;
}

bool NativeXtbMethod::scopedBool(const char* key, bool fallback) const
{
    const json* v = nullptr;
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object() && m_parameters["xtb"].contains(key))
        v = &m_parameters["xtb"][key];
    else if (m_parameters.contains(key))
        v = &m_parameters[key];
    if (!v) return fallback;
    if (v->is_boolean()) return v->get<bool>();
    if (v->is_number()) return v->get<double>() != 0.0;
    if (v->is_string()) return *v == "true" || *v == "1";
    return fallback;
}

void NativeXtbMethod::handleError(const std::string& operation)
{
    m_has_error = true;
    m_error_message = fmt::format("{} error during {}", getMethodName(), operation);
}

std::string NativeXtbMethod::largeSystemModeString() const
{
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object()
        && m_parameters["xtb"].contains("large_system_mode"))
        return m_parameters["xtb"]["large_system_mode"].get<std::string>();
    return m_parameters.value("large_system_mode", std::string("none"));
}

std::string NativeXtbMethod::eigensolverString() const
{
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object()
        && m_parameters["xtb"].contains("eigensolver"))
        return m_parameters["xtb"]["eigensolver"].get<std::string>();
    return m_parameters.value("eigensolver", std::string("mkl"));
}

std::string NativeXtbMethod::scfPrecondString() const
{
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object()
        && m_parameters["xtb"].contains("scf_precond") && m_parameters["xtb"]["scf_precond"].is_string())
        return m_parameters["xtb"]["scf_precond"].get<std::string>();
    if (m_parameters.contains("scf_precond") && m_parameters["scf_precond"].is_string())
        return m_parameters["scf_precond"].get<std::string>();
    return "none";
}

double NativeXtbMethod::electronicTemperature() const
{
    if (m_parameters.contains("xtb") && m_parameters["xtb"].is_object()
        && m_parameters["xtb"].contains("electronic_temperature"))
        return m_parameters["xtb"]["electronic_temperature"].get<double>();
    return m_parameters.value("electronic_temperature", 300.0);
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

        // large_system_mode path: when mode != none, delegate to the
        // fragmentation driver (it owns its own per-fragment XTB instances).
        const std::string mode_str = largeSystemModeString();
        const auto mode = curcuma::xtb::parseLargeSystemMode(mode_str);
        if (m_cpcmx && mode != curcuma::xtb::LargeSystemMode::None) {
            m_has_error = true;
            m_error_message = "solvent_model cpcmx is not available with -large_system_mode " + mode_str;
            return false;
        }
        if (mode != curcuma::xtb::LargeSystemMode::None) {
            // Hard-error: large_system_mode=fragments|dc + -eigensolver=purify
            // requires -electronic_temperature 0 (purification is 0-K by
            // construction). large_system_mode=sparse is unaffected (sparse IS
            // the eigensolver replacement, no m_eigensolver involved).
            if (mode != curcuma::xtb::LargeSystemMode::Sparse
                && eigensolverString() == "purify"
                && electronicTemperature() > 0.0) {
                m_has_error = true;
                m_error_message = fmt::format(
                    "{}: -eigensolver=purify is 0-K only; set -electronic_temperature 0 "
                    "(currently {:.1f} K) or use -eigensolver={}",
                    getMethodName(), electronicTemperature(),
                    mode == curcuma::xtb::LargeSystemMode::Fragments
                        ? "mkl|native|lobpcg" : "mkl|native|lobpcg");
                return false;
            }
            // -scf_precond acts on the dense CPU Broyden loop only; the large-system drivers own
            // their XTB objects, so strip the switch from their config (Claude Generated, Oct 2026).
            json driver_cfg = m_parameters;
            const std::string precond = scfPrecondString();
            if (!precond.empty() && precond != "none") {
                CurcumaLogger::warn(fmt::format("scf_precond '{}' ignored: -large_system_mode '{}' is active "
                                                "(dense CPU SCF only); scalar Broyden seed is used",
                                                precond, mode_str));
                if (driver_cfg.contains("xtb") && driver_cfg["xtb"].is_object())
                    driver_cfg["xtb"]["scf_precond"] = "none";
                driver_cfg["scf_precond"] = "none";
            }
            m_c1_driver = std::make_unique<curcuma::xtb::FragmentScfDriver>(m_method, driver_cfg);
            m_c1_driver->setIntraThreads(m_thread_count);
            if (!m_c1_driver->setMolecule(mol)) {
                handleError("large_system driver molecule initialization");
                return false;
            }
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format("{}: -large_system_mode '{}' active",
                                                getMethodName(), mode_str));
            return true;
        }
        m_c1_driver.reset();

        if (!m_xtb->QMInterface::InitialiseMolecule(mol)) {
            // Surface the engine's specific reason (e.g. unsupported d-shell element)
            // rather than a generic message. Claude Generated.
            handleError(m_xtb->hasError() ? m_xtb->errorMessage()
                                          : std::string("molecule initialization"));
            return false;
        }
        if (m_cpcmx) {
            if (!m_xtb_conductor->QMInterface::InitialiseMolecule(mol)) {
                handleError(m_xtb_conductor->hasError() ? m_xtb_conductor->errorMessage()
                                                        : std::string("CPCM-X conductor COSMO initialization"));
                return false;
            }
            if (!m_cpcmx_solvent_data) {   // solvent sigma segments from the database, once per method object
                const std::string dir = curcuma::solvation::cpcmx::findDatabase(scopedString("cpcmx_db", ""));
                auto data = std::make_unique<curcuma::solvation::cpcmx::Molecule>();
                std::string err;
                if (dir.empty() || !curcuma::solvation::cpcmx::loadSolvent(dir, m_cpcmx_solvent, *data, err)) {
                    m_has_error = true;
                    m_error_message = dir.empty() ? std::string("CPCM-X solvent database not found (set -xtb.cpcmx_db)") : err;
                    return false;
                }
                m_cpcmx_solvent_data = std::move(data);
            }
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

        if (m_c1_driver) {
            if (!m_c1_driver->updateGeometry(geometry)) {
                handleError("large_system driver geometry update");
                return false;
            }
            return true;
        }

        if (!m_xtb->UpdateMolecule(geometry)) {
            handleError("geometry update");
            return false;
        }
        if (m_xtb_conductor && !m_xtb_conductor->UpdateMolecule(geometry)) {
            handleError("CPCM-X conductor geometry update");
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

        // large_system_mode path: delegate the whole calculation to the driver.
        if (m_c1_driver) {
            m_last_energy = m_c1_driver->calculate(gradient);
            m_calculation_done = true;
            if (m_c1_driver->hasError()) {
                handleError(m_c1_driver->errorMessage());
                return 0.0;
            }
            return m_last_energy;
        }

        // Apply controller settings each call so geometry steps in -opt / MD pick
        // up the configuration (D4 charge source + SCF convergence settings).
        applyConfig();

        if (m_cpcmx) {
            // xtb has no CPCM-X gradient either (cpx.F90 computes energies only); refuse instead of returning the
            // gas-phase gradient with a solvated energy.
            if (gradient) {
                m_has_error = true;
                m_error_message = "solvent_model cpcmx has no analytic gradient (as xtb --cpcmx): single points only "
                                  "(-sp without -dump_gradient); for -opt, -md, -hessian use cosmo, cpcm or alpb";
                return 0.0;
            }
            m_last_energy = calculateCpcmx();
            m_calculation_done = !m_has_error;
            return m_has_error ? 0.0 : m_last_energy;
        }

        m_last_energy = m_xtb->Calculation(gradient);
        m_calculation_done = true;

        // Fail-loud (Claude Generated): a genuine solve breakdown (singular overlap,
        // eigensolver/density failure, NaN) must not be reported as a valid E=0. The
        // engine sets a hard-error flag distinct from benign max-iter non-convergence;
        // surface it so EnergyCalculator refuses the result instead of handing E=0 to
        // an optimiser/MD.
        if (m_xtb->hasError()) {
            handleError(m_xtb->errorMessage());
            return 0.0;
        }

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

/**
 * CPCM-X single point as xtb 6.7.1 --cpcmx (prog/main.F90 + solv/cpx.F90):
 *   1. SCF in the conductor COSMO set (xtb --cosmo inf): E_COSMO and the surface segments (xtb.cosmo)
 *   2. gas-phase SCF: E_gas (xtb: singlepoint of the calculator without solvation)
 *   3. post-SCF CPCM-X: average_charge, init_bonding, solv('crs', 298.15 K, 500 cycles, 1e-4 kcal/mol),
 *      state_correction, cds(probe 0.4); E = E_gas + dG_solv.
 * Claude Generated (Nov 2026).
 */
double NativeXtbMethod::calculateCpcmx()
{
    namespace cx = curcuma::solvation::cpcmx;
    const double e_cond = m_xtb_conductor->Calculation(false);
    if (m_xtb_conductor->hasError()) {
        handleError("CPCM-X conductor SCF: " + m_xtb_conductor->errorMessage());
        return 0.0;
    }
    const auto* cosmo = dynamic_cast<const DdCosmoSolvation*>(m_xtb_conductor->solvationModel());
    DdCosmoSolvation::Segments seg;
    if (!cosmo || !cosmo->segments(m_xtb_conductor->getCharges(), seg)) {
        handleError("CPCM-X: no conductor COSMO segments");
        return 0.0;
    }
    const double e_gas = m_xtb->Calculation(false);
    if (m_xtb->hasError()) {
        handleError("CPCM-X gas-phase SCF: " + m_xtb->errorMessage());
        return 0.0;
    }

    // Diagnostic dump (env CURCUMA_CPCMX_DUMP=<file>): the segment lines of xtb's xtb.cosmo (cosmo.f90 writeCosmoFile
    // format '(2i5, 7(1x, f14.9))') for a line-by-line comparison with xtb --cpcmx. No effect on the result.
    if (const char* dump = std::getenv("CURCUMA_CPCMX_DUMP"); dump && *dump) {
        if (std::FILE* f = std::fopen(dump, "w")) {
            const double autoaa = ALPBParameters::autoaa;
            for (size_t i = 0; i < seg.atom.size(); ++i)
                std::fprintf(f, "%5zu%5d %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n", i + 1, seg.atom[i] + 1,
                             seg.position(0, i), seg.position(1, i), seg.position(2, i), seg.charge(i), seg.area(i),
                             seg.charge(i) / seg.area(i), seg.potential(i) / autoaa);
            std::fclose(f);
        }
    }

    cx::ConductorData data;
    data.segment_atom = seg.atom;
    data.segment_xyz = seg.position;
    data.charge = seg.charge;
    data.area = seg.area;
    data.potential = seg.potential;
    data.z = m_molecule.m_atoms;
    data.atom_xyz = seg.atom_position.transpose();
    data.energy = e_cond;
    cx::Molecule solute = cx::soluteFromConductor(data, scopedBool("cpcmx_reference_rounding", true));
    solute.energy_gas = e_gas;
    cx::Molecule solvent = *m_cpcmx_solvent_data;
    std::string err;
    cx::OxygenRadiusRule rule;
    if (!cx::parseOxygenRadiusRule(scopedString("cpcmx_smd_oxygen_radius", "atom8"), rule)) {
        handleError("CPCM-X: unknown cpcmx_smd_oxygen_radius '" + scopedString("cpcmx_smd_oxygen_radius", "")
                    + "' (atom8, oxygen or none)");
        return 0.0;
    }
    if (!cx::calculate(solute, solvent, m_cpcmx_solvent, 298.15, 500, 1.0e-4, 0.4, rule, m_cpcmx_result, err)) {
        handleError(err);
        return 0.0;
    }
    m_cpcmx_conductor_energy = e_cond;
    m_cpcmx_gas_energy = e_gas;
    const auto& r = m_cpcmx_result;
    const double total = e_gas + r.dG();
    constexpr double kcal = 627.509;   // CPCM-X autokcal
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("CPCM-X post-SCF solvation ({}, xtb parameter set, port of CPCM-X 1.1.0)", m_cpcmx_solvent));
        CurcumaLogger::param("conductor COSMO energy", fmt::format("{:.10f} Eh ({} segments)", e_cond, seg.atom.size()));
        CurcumaLogger::param("Ideal state (dG_IS)", fmt::format("{:.10f} Eh {:.5f} kcal/mol", r.dG_is, r.dG_is * kcal));
        CurcumaLogger::param("Averaging correction (dG_av)", fmt::format("{:.10f} Eh {:.5f} kcal/mol", r.dG_cc, r.dG_cc * kcal));
        CurcumaLogger::param("Restoring free energy (dG_res)", fmt::format("{:.10f} Eh {:.5f} kcal/mol", r.dG_res, r.dG_res * kcal));
        CurcumaLogger::param("SMD contribution (dG_CDS)", fmt::format("{:.10f} Eh {:.5f} kcal/mol", r.dG_smd, r.dG_smd * kcal));
        CurcumaLogger::param("Standard state correction (dG_corr)", fmt::format("{:.10f} Eh {:.5f} kcal/mol", r.dG_ss, r.dG_ss * kcal));
        CurcumaLogger::param("Empirical shift (dG_shift)", fmt::format("{:.10f} Eh", r.dG_shift));
        CurcumaLogger::param("solvent profile cycles / N_ring", fmt::format("{} / {}", r.solvent_iterations, r.near));
    }
    if (CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::result(fmt::format("CPCM-X ({}): dG_solv = {:.10f} Eh ({:.5f} kcal/mol), E_gas = {:.10f} Eh, "
                                          "total free energy = {:.10f} Eh",
                                          m_cpcmx_solvent, r.dG(), r.dG() * kcal, e_gas, total));
    return total;
}

// ---------------------------------------------------------------------------
// Warm-start / iterative-mode controls (Claude Generated)
// ---------------------------------------------------------------------------
void NativeXtbMethod::setWarmStart(bool on)
{
    if (m_xtb) m_xtb->setWarmStart(on);
}


void NativeXtbMethod::resetWarmStart()
{
    if (m_xtb) m_xtb->resetWarmStart();
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
    if (m_c1_driver)
        return m_calculation_done ? m_c1_driver->getGradient()
                                  : Matrix::Zero(m_molecule.AtomCount(), 3);

    if (!m_calculation_done || !m_xtb)
        return Matrix::Zero(m_molecule.AtomCount(), 3);

    Matrix grad = m_xtb->Gradient();
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::param("native_xtb_gradient_norm", fmt::format("{:.6e} Eh/Angstrom", grad.norm()));
    return grad;
}

Vector NativeXtbMethod::getCharges() const
{
    if (m_c1_driver)
        return m_calculation_done ? m_c1_driver->getCharges()
                                  : Vector::Zero(m_molecule.AtomCount());
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
    if (m_c1_driver)
        return m_calculation_done ? m_c1_driver->getEnergyDecomposition() : json{};
    if (!m_calculation_done || !m_xtb) return json{};
    try {
        json d = m_xtb->getEnergyDecomposition();
        if (m_cpcmx && m_calculation_done) {   // CPCM-X contributions (Eh); the other terms are those of the gas phase
            const auto& r = m_cpcmx_result;
            d["cpcmx_dG_is"] = r.dG_is;
            d["cpcmx_dG_av"] = r.dG_cc;
            d["cpcmx_dG_res"] = r.dG_res;
            d["cpcmx_dG_cds"] = r.dG_smd;
            d["cpcmx_dG_ss"] = r.dG_ss;
            d["cpcmx_dG_shift"] = r.dG_shift;
            d["cpcmx_dG_solv"] = r.dG();
            d["cpcmx_E_gas"] = m_cpcmx_gas_energy;
            d["cpcmx_E_conductor"] = m_cpcmx_conductor_energy;
        }
        return d;
    } catch (const std::exception&) { return json{}; }
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
