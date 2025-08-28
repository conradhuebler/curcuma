/*
 * < TBLite Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "src/core/curcuma_logger.h"
#include "src/tools/general.h"

#include "tblite_method.h"

TBLiteMethod::TBLiteMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    CurcumaLogger::info("TBLiteMethod constructor called");
    CurcumaLogger::param("method_name", method_name);

#ifdef USE_TBLITE
    try {
        CurcumaLogger::info("Creating TBLiteInterface with method: " + method_name);
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("config", config.dump(2));
        }
        m_tblite = std::make_unique<TBLiteInterface>(config);

        // CRITICAL: Set the specific method name in TBLiteInterface
        CurcumaLogger::info("Setting TBLite method to: " + method_name);
        m_tblite->setMethod(method_name);
        CurcumaLogger::success("TBLiteInterface created and configured successfully for method: " + method_name);
    } catch (const std::exception& e) {
        CurcumaLogger::error("TBLiteInterface creation failed for method '" + method_name + "': " + std::string(e.what()));
        throw;
    }
#else
    CurcumaLogger::error("USE_TBLITE not defined - TBLite not available for method: " + method_name);
    throw std::runtime_error("TBLite not available in this build");
#endif
    m_parameters = config;
}

bool TBLiteMethod::setMolecule(const Mol& mol) {
    m_molecule = mol;
    m_initialized = true;
#ifdef USE_TBLITE
    return m_tblite->InitialiseMolecule(mol);
#else
    return false;
#endif
}

bool TBLiteMethod::updateGeometry(const Matrix& geometry) {
#ifdef USE_TBLITE
    return m_tblite->UpdateMolecule(geometry);
#else
    return false;
#endif
}

double TBLiteMethod::calculateEnergy(bool gradient)
{
#ifdef USE_TBLITE
    m_last_energy = m_tblite->Calculation(gradient);
    m_calculation_done = true;
    return m_last_energy;
#else
    return 0.0;
#endif
}

Matrix TBLiteMethod::getGradient() const {
#ifdef USE_TBLITE
    return m_tblite->Gradient();
#else
    return Matrix::Zero(1, 3);
#endif
}

Vector TBLiteMethod::getCharges() const {
#ifdef USE_TBLITE
    return m_tblite->Charges();
#else
    return Vector::Zero(1);
#endif
}

Vector TBLiteMethod::getBondOrders() const {
#ifdef USE_TBLITE
    return m_tblite->BondOrders();
#else
    return Vector{};
#endif
}

Position TBLiteMethod::getDipole() const {
#ifdef USE_TBLITE
    Vector dipole = m_tblite->Dipole();
    if (dipole.size() >= 3) {
        return Position{dipole(0), dipole(1), dipole(2)};
    }
#endif
    return Position{0.0, 0.0, 0.0};
}

void TBLiteMethod::setThreadCount(int threads) { m_thread_count = threads; }
void TBLiteMethod::setParameters(const json& params) { m_parameters = params; }
json TBLiteMethod::getParameters() const { return m_parameters; }
bool TBLiteMethod::hasError() const { return m_has_error; }
void TBLiteMethod::clearError() { m_has_error = false; m_error_message.clear(); }
std::string TBLiteMethod::getErrorMessage() const { return m_error_message; }

Vector TBLiteMethod::getOrbitalEnergies() const {
#ifdef USE_TBLITE
    return m_tblite->OrbitalEnergies();
#else
    return Vector{};
#endif
}

Vector TBLiteMethod::getOrbitalOccupations() const {
#ifdef USE_TBLITE
    return m_tblite->OrbitalOccupations();
#else
    return Vector{};
#endif
}

bool TBLiteMethod::isAvailable() {
#ifdef USE_TBLITE
    return true;
#else
    return false;
#endif
}

std::vector<std::string> TBLiteMethod::getSupportedMethods() {
    return {"ipea1", "gfn1", "gfn2"};
}