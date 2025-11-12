/*
 * <PM3 Method Wrapper Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * This program is free software under GPL-3.0
 */

#include "pm3_method.h"
#include "src/core/curcuma_logger.h"

PM3Method::PM3Method(const json& config)
    : m_pm3(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_pm3 = std::make_unique<PM3>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("PM3Method initialized with native implementation");
    }
}

bool PM3Method::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    m_pm3->setMolecule(mol);
    return m_pm3->InitialiseMolecule();
}

bool PM3Method::updateGeometry(const Matrix& geometry)
{
    m_pm3->updateGeometry(geometry);
    m_calculation_done = false;

    return m_pm3->InitialiseMolecule();
}

double PM3Method::calculateEnergy(bool gradient)
{
    m_last_energy = m_pm3->Calculation(gradient);
    m_calculation_done = true;

    return m_last_energy;
}

Matrix PM3Method::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("PM3Method: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }

    return m_pm3->Gradient();
}

Vector PM3Method::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }

    return m_pm3->getPartialCharges();
}

json PM3Method::getEnergyDecomposition() const
{
    return m_pm3->getEnergyDecomposition();
}

json PM3Method::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_damping", 0.4 },
        { "threads", 1 }
    };
}
