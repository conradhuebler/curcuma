/*
 * <Unified NDDO Method Wrapper Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated: Unified NDDO method wrapper for MethodFactory integration
 *
 * This program is free software under GPL-3.0
 */

#include "nddo_method.h"
#include "src/core/curcuma_logger.h"

#include <algorithm>

NDDOMethod::NDDOMethod(NDDOMethodType type, const json& config)
    : m_nddo(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_nddo = std::make_unique<NDDO>(type);
    m_method_name = m_nddo->getMethodNameStr();

    // Convert to lowercase for MethodFactory consistency
    std::transform(m_method_name.begin(), m_method_name.end(),
                   m_method_name.begin(), ::tolower);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("NDDOMethod initialized: " + m_nddo->getMethodNameStr());
    }
}

bool NDDOMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    return m_nddo->QMInterface::InitialiseMolecule(mol);
}

bool NDDOMethod::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;
    return m_nddo->UpdateMolecule(geometry);
}

double NDDOMethod::calculateEnergy(bool gradient)
{
    m_last_energy = m_nddo->Calculation(gradient);
    m_calculation_done = true;
    return m_last_energy;
}

Matrix NDDOMethod::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("NDDOMethod: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }
    return m_nddo->Gradient();
}

Vector NDDOMethod::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }
    return m_nddo->getPartialCharges();
}

json NDDOMethod::getEnergyDecomposition() const
{
    return m_nddo->getEnergyDecomposition();
}

json NDDOMethod::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_damping", 0.4 },
        { "threads", 1 }
    };
}
