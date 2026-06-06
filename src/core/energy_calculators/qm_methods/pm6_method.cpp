/*
 * <PM6 Method Wrapper Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * This program is free software under GPL-3.0
 *
 * Claude Generated: PM6 method wrapper for MethodFactory integration
 */

#include "pm6_method.h"
#include "src/core/curcuma_logger.h"

PM6Method::PM6Method(const json& config)
    : m_pm6(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_pm6 = std::make_unique<PM6>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("PM6Method initialized with native implementation");
    }
}

bool PM6Method::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    // Initialize with molecule data (explicit base class call to avoid name hiding)
    return m_pm6->QMInterface::InitialiseMolecule(mol);
}

bool PM6Method::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;

    // Update geometry only
    return m_pm6->UpdateMolecule(geometry);
}

double PM6Method::calculateEnergy(bool gradient)
{
    m_last_energy = m_pm6->Calculation(gradient);
    m_calculation_done = true;

    return m_last_energy;
}

Matrix PM6Method::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("PM6Method: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }

    return m_pm6->Gradient();
}

Vector PM6Method::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }

    return m_pm6->getPartialCharges();
}

json PM6Method::getEnergyDecomposition() const
{
    return m_pm6->getEnergyDecomposition();
}

json PM6Method::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_damping", 0.4 },
        { "threads", 1 }
    };
}
