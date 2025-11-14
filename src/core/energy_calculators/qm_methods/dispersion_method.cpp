/*
 * < Dispersion Correction Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "src/core/energy_calculators/qm_methods/interface/abstract_interface.h"

#include "dispersion_method.h"
#include "src/core/config_manager.h"

DispersionMethod::DispersionMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name), m_calculation_done(false), m_last_energy(0.0) {

    // Wrap JSON in ConfigManager for interface constructors
    if (method_name == "d3") {
#ifdef USE_D3
        ConfigManager d3_config("dftd3", config);
        m_dispersion = std::make_unique<DFTD3Interface>(d3_config);
#endif
    } else if (method_name == "d4") {
#ifdef USE_D4
        ConfigManager d4_config("dftd4", config);
        m_dispersion = std::make_unique<DFTD4Interface>(d4_config);
#endif
    }
    m_parameters = config;
}

bool DispersionMethod::setMolecule(const Mol& mol) {
    m_molecule = mol;
    m_initialized = true;
    if (m_dispersion) {
        return m_dispersion->InitialiseMolecule(mol);
    }
    return false;
}

bool DispersionMethod::updateGeometry(const Matrix& geometry) {
    if (m_dispersion) {
        return m_dispersion->UpdateMolecule(geometry);
    }
    return false;
}

double DispersionMethod::calculateEnergy(bool gradient)
{
    if (m_dispersion) {
        m_last_energy = m_dispersion->Calculation(gradient);
        m_calculation_done = true;
        return m_last_energy;
    }
    return 0.0;
}

Matrix DispersionMethod::getGradient() const {
    if (m_dispersion) {
        return m_dispersion->Gradient();
    }
    return Matrix::Zero(1, 3);
}

Vector DispersionMethod::getCharges() const {
    if (m_dispersion) {
        return m_dispersion->Charges();
    }
    return Vector::Zero(1);
}

Vector DispersionMethod::getBondOrders() const {
    return Vector{};
}

Position DispersionMethod::getDipole() const {
    if (m_dispersion) {
        Vector dipole = m_dispersion->Dipole();
        if (dipole.size() >= 3) {
            return Position{dipole(0), dipole(1), dipole(2)};
        }
    }
    return Position{0.0, 0.0, 0.0};
}

void DispersionMethod::setThreadCount(int threads) { m_thread_count = threads; }
void DispersionMethod::setParameters(const json& params) { m_parameters = params; }
json DispersionMethod::getParameters() const { return m_parameters; }
bool DispersionMethod::hasError() const { return m_has_error; }
void DispersionMethod::clearError() { m_has_error = false; m_error_message.clear(); }
std::string DispersionMethod::getErrorMessage() const { return m_error_message; }

bool DispersionMethod::isD3Available() {
#ifdef USE_D3
    return true;
#else
    return false;
#endif
}

bool DispersionMethod::isD4Available() {
#ifdef USE_D4
    return true;
#else
    return false;
#endif
}

std::vector<std::string> DispersionMethod::getSupportedMethods() {
    std::vector<std::string> methods;
    if (isD3Available()) methods.push_back("d3");
    if (isD4Available()) methods.push_back("d4");
    return methods;
}