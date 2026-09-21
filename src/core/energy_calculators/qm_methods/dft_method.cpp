/*
 * <Native KS-DFT Method Wrapper Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated: Native KS-DFT method wrapper for MethodFactory integration
 *
 * This program is free software under GPL-3.0
 */

#include "dft_method.h"
#include "src/core/curcuma_logger.h"

#include <algorithm>

DFTMethod::DFTMethod(DFTFunctional functional, const json& config)
    : m_dft(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) {
        full_config.merge_patch(config);
    }

    m_dft = std::make_unique<DFT>(functional, full_config);
    m_method_name = m_dft->getMethodNameStr();

    // Convert to lowercase for MethodFactory consistency
    std::transform(m_method_name.begin(), m_method_name.end(),
                   m_method_name.begin(), ::tolower);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("DFTMethod initialized: " + m_dft->getMethodNameStr());
    }
}

bool DFTMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;

    return m_dft->QMInterface::InitialiseMolecule(mol);
}

bool DFTMethod::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;
    return m_dft->UpdateMolecule(geometry);
}

double DFTMethod::calculateEnergy(bool gradient)
{
    m_last_energy = m_dft->Calculation(gradient);
    m_calculation_done = true;
    return m_last_energy;
}

Matrix DFTMethod::getGradient() const
{
    if (!m_calculation_done) {
        CurcumaLogger::warn("DFTMethod: No calculation done yet");
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    }
    // WP0 scaffold: analytic gradient arrives in WP8.
    return Matrix::Zero(m_molecule.AtomCount(), 3);
}

Vector DFTMethod::getCharges() const
{
    if (!m_calculation_done) {
        return Vector::Zero(m_molecule.AtomCount());
    }
    // WP0 scaffold: no population analysis yet.
    return Vector::Zero(m_molecule.AtomCount());
}

json DFTMethod::getEnergyDecomposition() const
{
    // WP0 scaffold: only the nuclear repulsion component is meaningful.
    json decomp;
    decomp["nuclear_repulsion"] = m_last_energy;
    decomp["electronic"] = 0.0;
    decomp["total"] = m_last_energy;
    return decomp;
}

json DFTMethod::getDefaultConfig()
{
    return json{
        { "basis", "def2-SVP" },
        { "grid", "sg1" },
        { "scf_max_iterations", 100 },
        { "scf_threshold", 1.0e-6 },
        { "scf_mode", "diis" },
        { "threads", 1 },
        { "cartesian_d", false }
    };
}