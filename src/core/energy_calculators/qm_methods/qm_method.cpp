/*
 * <Native ab-initio QM Method Wrapper Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated: Native KS-DFT method wrapper for MethodFactory integration,
 *                   renamed DFTMethod -> QMMethod (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#include "qm_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/method_factory.h"
#include "src/core/parameter_registry.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <algorithm>

json QMMethod::engineConfig(const json& config)
{
    json full_config = ParameterRegistry::getInstance().getDefaultJson("qm");
    // The CLI auto-routes every -qm.<param> (and its flat form) into the "qm"
    // module scope, i.e. controller["qm"], which is the ONLY place they land --
    // so the scope has to be merged explicitly. Reading the controller's top level
    // alone (as this wrapper once did) silently dropped every scoped flag.
    // Order: top level (flat/legacy fallback) < "dft" (the pre-Sep-2026 module
    // name, kept so old command lines and -import_config files still work) < "qm".
    //
    // Oct 27 2026 (master's standard for the shared SCF knobs, docs/SCF_MODES.md): the shared module "scf"
    // (scf_parameters.h; -scf_threshold, -scf.scf_threshold, -scf_max_iterations, ...) sits between "dft" and
    // "qm", as for gfn1/gfn2 (applyXtbScfConfig: method scope > scf > default). Keys of the qm and dft scopes
    // are mapped onto their canonical qm names, so the old spellings -qm.diis_start / -qm.diis_subspace keep
    // working; the top level is not alias-resolved (its diis_start / max_iterations are the optimiser's).
    auto& reg = ParameterRegistry::getInstance();
    auto merge_scope = [&](const json& scope) {
        for (const auto& [key, value] : scope.items()) {
            const std::string canon = reg.resolveAlias("qm", key);
            full_config[canon.empty() ? key : canon] = value;
        }
    };
    if (!config.empty())
        full_config.merge_patch(config);
    if (config.contains("dft") && config["dft"].is_object())
        merge_scope(config["dft"]);
    if (config.contains("scf") && config["scf"].is_object()) {
        json shared = config;  // aliases resolved and "method default" markers dropped (idempotent)
        MethodFactory::applySharedScfScope(shared, "hf");
        if (shared.contains("scf") && shared["scf"].is_object())
            for (const auto& [key, value] : shared["scf"].items())
                full_config[key] = value;  // the scf names are the qm names (scf_threshold, scf_mode, ...)
    }
    if (config.contains("qm") && config["qm"].is_object())
        merge_scope(config["qm"]);
    return full_config;
}

QMMethod::QMMethod(QMFunctional functional, const json& config)
    : m_engine(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_engine = std::make_unique<QMEngine>(functional, engineConfig(config));
    m_method_name = m_engine->getMethodNameStr();

    // Convert to lowercase for MethodFactory consistency
    std::transform(m_method_name.begin(), m_method_name.end(),
                   m_method_name.begin(), ::tolower);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("QMMethod initialized: " + m_engine->getMethodNameStr());
    }
}

bool QMMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    m_engine->resetForNewMolecule();  // a reused object must not keep the old molecule
    return m_engine->QMInterface::InitialiseMolecule(mol);
}

bool QMMethod::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;
    // QMInterface::UpdateMolecule(Matrix) stores the geometry and calls the
    // engine's UpdateMolecule(), which invalidates the integral/SCF caches.
    return m_engine->UpdateMolecule(geometry);
}

double QMMethod::calculateEnergy(bool gradient)
{
    m_has_error = false;
    m_error_message.clear();
    m_last_energy = m_engine->Calculation(gradient);
    m_calculation_done = true;
    if (m_engine->riSetupFailed()) {
        // Claude Generated (Oct 2026, WP10): no energy without the fit the user asked for.
        m_has_error = true;
        m_error_message = "density fitting (-qm.ri " + m_engine->riMode() + ") setup failed, see the error above";
    } else if (m_engine->getFunctional() == QMFunctional::HF && !m_engine->scfConverged()
               && !m_engine->scfAllowUnconverged()) {
        // Claude Generated (Oct 27 2026): master's rule for every SCF method (docs/SCF_MODES.md "Failure handling"):
        // an unconverged SCF is an error, not a result. Before, HF printed an error and returned the energy (exit 0).
        m_has_error = true;
        m_error_message = fmt::format("native HF SCF did not converge within {} iterations; "
                                      "-scf_allow_unconverged true accepts the unconverged energy",
                                      m_engine->scfMaxIter());
    } else if (gradient && m_engine->gradientFailed()) {
        // Claude Generated (Oct 2026): fail loud instead of handing back a zero gradient
        // (RI gradient not possible, e.g. memory limit for the 3-index densities).
        m_has_error = true;
        m_error_message = m_engine->gradientError();
    }
    return m_last_energy;
}

Matrix QMMethod::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("QMMethod: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }
    if (!m_engine->hasGradient())
        return Matrix::Zero(m_molecule.AtomCount(), 3);  // lda/pbe/b3lyp: no V_xc, no gradient
    // The engine works in Eh/Bohr; ComputationalMethod::getGradient() is Eh/Angstrom
    // (dE/dx[A] = dE/dx[Bohr] * Bohr-per-Angstrom), see Known Issue #28.
    return m_engine->gradientBohr() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
}

Vector QMMethod::getCharges() const
{
    // No population analysis yet.
    return Vector::Zero(m_molecule.AtomCount());
}

json QMMethod::getEnergyDecomposition() const
{
    json decomp;
    decomp["total"] = m_last_energy;
    if (m_engine->getFunctional() == QMFunctional::HF && m_engine->scfConverged()) {
        // {ET, EV, EJ, Ex, E_elec}; E_nn is the remainder of the total.
        const std::vector<double> c = m_engine->energyComponents();
        decomp["kinetic"] = c[0];
        decomp["nuclear_attraction"] = c[1];
        decomp["coulomb"] = c[2];
        decomp["exchange"] = c[3];
        decomp["electronic"] = c[4];
        decomp["nuclear_repulsion"] = m_last_energy - c[4];
    } else {
        // DFT functionals without V_xc yet: the total IS the nuclear repulsion.
        decomp["nuclear_repulsion"] = m_last_energy;
        decomp["electronic"] = 0.0;
    }
    return decomp;
}
