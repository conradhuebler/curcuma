/*
 * <AM1 Method Wrapper Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * This program is free software under GPL-3.0
 *
 * Claude Generated: AM1 method wrapper for MethodFactory integration
 */

#include "am1_method.h"
#include "src/core/curcuma_logger.h"

AM1Method::AM1Method(const json& config)
    : m_am1(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_am1 = std::make_unique<AM1>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("AM1Method initialized with native implementation");
    }
}

bool AM1Method::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    // Initialize with molecule data (explicit base class call to avoid name hiding)
    return m_am1->QMInterface::InitialiseMolecule(mol);
}

bool AM1Method::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;

    // Update geometry only
    return m_am1->UpdateMolecule(geometry);
}

double AM1Method::calculateEnergy(bool gradient)
{
    m_last_energy = m_am1->Calculation(gradient);
    m_calculation_done = true;

    return m_last_energy;
}

Matrix AM1Method::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("AM1Method: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }

    return m_am1->Gradient();
}

Vector AM1Method::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }

    return m_am1->getPartialCharges();
}

json AM1Method::getEnergyDecomposition() const
{
    return m_am1->getEnergyDecomposition();
}

json AM1Method::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_damping", 0.4 },
        { "threads", 1 }
    };
}
