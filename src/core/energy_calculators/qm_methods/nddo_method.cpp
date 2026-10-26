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
#include <cmath>
#include <stdexcept>
#include <fmt/format.h>

NDDOMethod::NDDOMethod(NDDOMethodType type, const json& config)
    : m_nddo(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_nddo = std::make_unique<NDDO>(type);
    m_method_name = m_nddo->getMethodNameStr();
    applyScfConfig(config);

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
    m_error_message.clear();
    m_last_energy = m_nddo->Calculation(gradient);
    m_calculation_done = true;
    // Claude Generated (Oct 2026): an unconverged or non-finite SCF is reported as an error
    // (EnergyCalculator then refuses the result), never as E = 0.
    if (!m_nddo->scfError().empty())
        m_error_message = m_nddo->scfError();
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

/**
 * @brief Hand the shared SCF settings (registry module "scf", scf_parameters.h) to NDDO.
 * Claude Generated (Oct 2026). Reads config["scf"] (aliases resolved and "method default"
 * markers dropped by MethodFactory::applySharedScfScope), then the same names at top level.
 * Absent keys keep the NDDO defaults: 100 iterations, ||dP||_F < 1e-6, mixing 0.6, DIIS from
 * iteration 3 over 8 Fock matrices. scf_mode broyden/level-shift exist only for gfn1/gfn2.
 */
void NDDOMethod::applyScfConfig(const json& config)
{
    // scf_diis_start / scf_diis_subspace carry the scf_ prefix because optimiser configs carry
    // their own top-level diis_start, which must not reach the SCF.
    auto find = [&](const char* key) -> const json* {
        if (config.contains("scf") && config["scf"].is_object() && config["scf"].contains(key))
            return &config["scf"][key];
        const std::string k(key);
        if (config.is_object() && config.contains(key))
            return &config[key];
        return nullptr;
    };
    auto as_int = [](const json& v) { return static_cast<int>(std::lround(v.get<double>())); };   // CLI numbers arrive as double

    if (const json* v = find("scf_max_iterations"); v && v->is_number() && as_int(*v) > 0)
        m_nddo->setScfMaxIterations(as_int(*v));
    if (const json* v = find("scf_threshold"); v && v->is_number() && v->get<double>() > 0.0)
        m_nddo->setScfThreshold(v->get<double>());
    if (const json* v = find("scf_damping"); v && v->is_number() && v->get<double>() > 0.0) {
        const double d = v->get<double>();
        if (d > 1.0)
            throw std::invalid_argument(fmt::format("scf_damping {} is outside (0, 1] for {}", d, m_method_name));
        m_nddo->setScfDamping(d);
    }
    if (const json* v = find("scf_diis_start"); v && v->is_number() && as_int(*v) >= 0)
        m_nddo->setDiisStart(as_int(*v));
    if (const json* v = find("scf_diis_subspace"); v && v->is_number() && as_int(*v) > 0)
        m_nddo->setDiisSubspace(std::max(2, as_int(*v)));
    if (const json* v = find("scf_allow_unconverged"); v) {
        if (v->is_boolean())
            m_nddo->setScfAllowUnconverged(v->get<bool>());
        else if (v->is_number())
            m_nddo->setScfAllowUnconverged(v->get<double>() != 0.0);
    }
    if (const json* v = find("scf_mode"); v && v->is_string()) {
        std::string mode = v->get<std::string>();
        std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
        if (mode == "diis" || mode == "pulay" || mode == "auto")
            m_nddo->setScfUseDiis(true);
        else if (mode == "plain" || mode == "damped")
            m_nddo->setScfUseDiis(false);
        else
            throw std::invalid_argument(fmt::format("scf_mode '{}' is not available for {}; use diis or plain",
                                                    v->get<std::string>(), m_method_name));
    }
}
