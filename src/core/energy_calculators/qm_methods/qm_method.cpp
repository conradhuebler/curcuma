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
#include "src/core/parameter_registry.h"
#include "src/core/units.h"

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
    if (!config.empty())
        full_config.merge_patch(config);
    if (config.contains("dft") && config["dft"].is_object())
        full_config.merge_patch(config["dft"]);
    if (config.contains("qm") && config["qm"].is_object())
        full_config.merge_patch(config["qm"]);
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
    m_last_energy = m_engine->Calculation(gradient);
    m_calculation_done = true;
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
