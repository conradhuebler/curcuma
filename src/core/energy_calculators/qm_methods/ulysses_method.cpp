/*
 * < Ulysses Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "ulysses_method.h"

UlyssesMethod::UlyssesMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name), m_calculation_done(false), m_last_energy(0.0) {
    std::cout << "[DEBUG] UlyssesMethod constructor called with method: " << method_name << std::endl;
#ifdef USE_ULYSSES
    try {
        std::cout << "[DEBUG] Creating UlyssesInterface..." << std::endl;
        std::cout << "[DEBUG] Config passed to UlyssesInterface: " << config.dump(2) << std::endl;
        m_ulysses = std::make_unique<UlyssesInterface>(config);
        std::cout << "[DEBUG] UlyssesInterface created successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "[ERROR] UlyssesInterface creation failed: " << e.what() << std::endl;
        throw;
    }
#else
    std::cout << "[ERROR] USE_ULYSSES not defined - Ulysses not available" << std::endl;
#endif
    m_parameters = config;
    std::cout << "[DEBUG] UlyssesMethod constructor completed" << std::endl;
}

bool UlyssesMethod::setMolecule(const Mol& mol) {
    m_molecule = mol;
    m_initialized = true;
#ifdef USE_ULYSSES
    return m_ulysses->QMInterface::InitialiseMolecule(mol);
#else
    (void)mol;
    return false;
#endif
}

bool UlyssesMethod::updateGeometry(const Matrix& geometry) {
#ifdef USE_ULYSSES
    return m_ulysses->QMInterface::UpdateMolecule(geometry);
#else
    (void)geometry;
    return false;
#endif
}

double UlyssesMethod::calculateEnergy(bool gradient)
{
#ifdef USE_ULYSSES
    m_last_energy = m_ulysses->Calculation(gradient);
    m_calculation_done = true;
    return m_last_energy;
#else
    return 0.0;
#endif
}

Matrix UlyssesMethod::getGradient() const {
#ifdef USE_ULYSSES
    return m_ulysses->Gradient();
#else
    return Matrix::Zero(1, 3);
#endif
}

Vector UlyssesMethod::getCharges() const {
#ifdef USE_ULYSSES
    return m_ulysses->Charges();
#else
    return Vector::Zero(1);
#endif
}

Vector UlyssesMethod::getBondOrders() const {
#ifdef USE_ULYSSES
    return m_ulysses->BondOrders();
#else
    return Vector{};
#endif
}

Position UlyssesMethod::getDipole() const {
#ifdef USE_ULYSSES
    Vector dipole = m_ulysses->Dipole();
    if (dipole.size() >= 3) {
        return Position{dipole(0), dipole(1), dipole(2)};
    }
#endif
    return Position{0.0, 0.0, 0.0};
}

void UlyssesMethod::setThreadCount(int threads) { m_thread_count = threads; }
void UlyssesMethod::setParameters(const json& params) { m_parameters = params; }
json UlyssesMethod::getParameters() const { return m_parameters; }
bool UlyssesMethod::hasError() const { return m_has_error; }
void UlyssesMethod::clearError() { m_has_error = false; m_error_message.clear(); }
std::string UlyssesMethod::getErrorMessage() const { return m_error_message; }

Vector UlyssesMethod::getOrbitalEnergies() const {
#ifdef USE_ULYSSES
    return m_ulysses->OrbitalEnergies();
#else
    return Vector{};
#endif
}

Vector UlyssesMethod::getOrbitalOccupations() const {
#ifdef USE_ULYSSES
    return m_ulysses->OrbitalOccupations();
#else
    return Vector{};
#endif
}

bool UlyssesMethod::isAvailable() {
#ifdef USE_ULYSSES
    return true;
#else
    return false;
#endif
}

std::vector<std::string> UlyssesMethod::getSupportedMethods() {
    return {"ugfn2", "GFN2L", "pm3", "PM3PDDG", "MNDOPDDG", "PM3BP", "RM1", "AM1", "MNDO", "MNDOd", "pm6"};
}