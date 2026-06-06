/*
 * <MNDO Method Wrapper Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * This program is free software under GPL-3.0
 *
 * Claude Generated: MNDO method wrapper for MethodFactory integration
 */

#include "mndo_method.h"
#include "src/core/curcuma_logger.h"

MNDOMethod::MNDOMethod(const json& config)
    : m_mndo(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_mndo = std::make_unique<MNDO>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("MNDOMethod initialized with native implementation");
    }
}

bool MNDOMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    // Initialize with molecule data (explicit base class call to avoid name hiding)
    return m_mndo->QMInterface::InitialiseMolecule(mol);
}

bool MNDOMethod::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;

    // Update geometry only
    return m_mndo->UpdateMolecule(geometry);
}

double MNDOMethod::calculateEnergy(bool gradient)
{
    m_last_energy = m_mndo->Calculation(gradient);
    m_calculation_done = true;

    return m_last_energy;
}

Matrix MNDOMethod::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("MNDOMethod: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }

    return m_mndo->Gradient();
}

Vector MNDOMethod::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }

    return m_mndo->getPartialCharges();
}

json MNDOMethod::getEnergyDecomposition() const
{
    return m_mndo->getEnergyDecomposition();
}

json MNDOMethod::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_damping", 0.4 },
        { "threads", 1 }
    };
}
