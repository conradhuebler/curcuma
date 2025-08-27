/*
 * < XTB Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "xtb_method.h"
#include "src/tools/general.h"

XTBMethod::XTBMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name), m_calculation_done(false), m_last_energy(0.0) {
#ifdef USE_XTB
    m_xtb = std::make_unique<XTBInterface>(config);
#endif
    m_parameters = config;
}

bool XTBMethod::setMolecule(const Mol& mol) {
    m_molecule = mol;
    m_initialized = true;
#ifdef USE_XTB
    return m_xtb->InitialiseMolecule(mol);
#else
    return false;
#endif
}

bool XTBMethod::updateGeometry(const Matrix& geometry) {
#ifdef USE_XTB
    return m_xtb->UpdateMolecule(geometry);
#else
    (void)geometry; // Suppress unused parameter warning
    return false;
#endif
}

double XTBMethod::calculateEnergy(bool gradient)
{
#ifdef USE_XTB
    m_last_energy = m_xtb->Calculation(gradient);
    m_calculation_done = true;
    return m_last_energy;
#else
    return 0.0;
#endif
}

Matrix XTBMethod::getGradient() const {
#ifdef USE_XTB
    return m_xtb->Gradient();
#else
    return Matrix::Zero(1, 3);
#endif
}

Vector XTBMethod::getCharges() const {
#ifdef USE_XTB
    return m_xtb->Charges();
#else
    return Vector::Zero(1);
#endif
}

Vector XTBMethod::getBondOrders() const {
#ifdef USE_XTB
    return m_xtb->BondOrders();
#else
    return Vector{};
#endif
}

Position XTBMethod::getDipole() const {
#ifdef USE_XTB
    Vector dipole = m_xtb->Dipole();
    if (dipole.size() >= 3) {
        return Position{dipole(0), dipole(1), dipole(2)};
    }
#endif
    return Position{0.0, 0.0, 0.0};
}

bool XTBMethod::isThreadSafe() const { return true; }
void XTBMethod::setThreadCount(int threads) { m_thread_count = threads; }
void XTBMethod::setParameters(const json& params) { m_parameters = params; }
json XTBMethod::getParameters() const { return m_parameters; }
bool XTBMethod::hasError() const { return m_has_error; }
void XTBMethod::clearError() { m_has_error = false; m_error_message.clear(); }
std::string XTBMethod::getErrorMessage() const { return m_error_message; }

Vector XTBMethod::getOrbitalEnergies() const {
#ifdef USE_XTB
    return m_xtb->OrbitalEnergies();
#else
    return Vector{};
#endif
}

Vector XTBMethod::getOrbitalOccupations() const {
#ifdef USE_XTB
    return m_xtb->OrbitalOccupations();
#else
    return Vector{};
#endif
}

bool XTBMethod::isAvailable() {
#ifdef USE_XTB
    return true;
#else
    return false;
#endif
}

std::vector<std::string> XTBMethod::getSupportedMethods() {
    return {"gfnff", "xtb-gfn1", "xtb-gfn2"};
}

bool XTBMethod::saveToFile(const std::string& filename) const {
    return false; // Stub
}