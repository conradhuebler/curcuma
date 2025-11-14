/*
 * <GFN1 Method Wrapper Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#include "gfn1_method.h"
#include "src/core/curcuma_logger.h"

GFN1Method::GFN1Method(const json& config)
    : m_gfn1(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    // Merge with defaults
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_gfn1 = std::make_unique<GFN1>();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("GFN1Method initialized with native implementation");
    }
}

bool GFN1Method::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    // Initialize with molecule data (explicit base class call to avoid name hiding)
    return m_gfn1->QMInterface::InitialiseMolecule(mol);
}

bool GFN1Method::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;

    // Update geometry only
    return m_gfn1->UpdateMolecule(geometry);
}

double GFN1Method::calculateEnergy(bool gradient)
{
    m_last_energy = m_gfn1->Calculation(gradient);
    m_calculation_done = true;

    return m_last_energy;
}

Matrix GFN1Method::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("GFN1Method: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }

    return m_gfn1->Gradient();
}

Vector GFN1Method::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }

    return m_gfn1->getPartialCharges();
}

Vector GFN1Method::getCoordinationNumbers() const
{
    return m_gfn1->getCoordinationNumbers();
}

json GFN1Method::getEnergyDecomposition() const
{
    return m_gfn1->getEnergyDecomposition();
}

bool GFN1Method::saveToFile(const std::string& filename) const
{
    // TODO: Implement if needed
    return false;
}

json GFN1Method::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_damping", 0.4 },
        { "threads", 1 },
        { "dispersion", "d3" }
    };
}
