/**
 * @file orca_method.cpp
 * @brief ORCA method wrapper implementation
 *
 * Claude Generated - June 2026
 */

#include "orca_method.h"
#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/tools/general.h"

#include <fmt/format.h>

OrcaMethod::OrcaMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name)
    , m_orca_keyword(methodToOrcaKeyword(method_name))
{
    m_parameters = config;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("OrcaMethod constructor called");
        CurcumaLogger::param("method_name", method_name);
        CurcumaLogger::param("orca_keyword", m_orca_keyword);
    }

    // Create ConfigManager and OrcaInterface
    ConfigManager orca_config("orca", config);
    m_orca = std::make_unique<OrcaInterface>(orca_config);
}

std::string OrcaMethod::methodToOrcaKeyword(const std::string& method)
{
    std::string m = method;
    std::transform(m.begin(), m.end(), m.begin(), ::tolower);

    if (m == "hf-3c")     return "HF-3c";
    if (m == "b97-3c")    return "B97-3c";
    if (m == "r2scan-3c") return "r2SCAN-3c";
    if (m == "pbeh-3c")   return "PBEh-3c";
    if (m == "orca")      return "";  // custom input

    // Unknown: pass through as-is (user might know what they're doing)
    return method;
}

bool OrcaMethod::initializeOrca()
{
    if (!m_orca) {
        m_has_error = true;
        m_error_message = "OrcaInterface not initialized";
        return false;
    }
    return true;
}

void OrcaMethod::updateMoleculeGeometry(const Matrix& geometry)
{
    if (geometry.rows() == 0) return;

    m_molecule.m_geometry.resize(geometry.rows(), 3);
    for (int i = 0; i < geometry.rows(); ++i) {
        m_molecule.m_geometry(i, 0) = geometry(i, 0);
        m_molecule.m_geometry(i, 1) = geometry(i, 1);
        m_molecule.m_geometry(i, 2) = geometry(i, 2);
    }
}

bool OrcaMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_initialized = true;

    if (!initializeOrca()) return false;

    // Apply charge/multiplicity overrides if set
    if (m_charge_override != 0)
        m_molecule.m_charge = m_charge_override;
    if (m_multiplicity_override > 1)
        m_molecule.m_spin = (m_multiplicity_override - 1) / 2.0;

    bool ok = m_orca->generateInput(m_molecule, m_orca_keyword);
    if (!ok) {
        m_has_error = m_orca->hasError();
        m_error_message = m_orca->getErrorMessage();
    }
    return ok;
}

bool OrcaMethod::updateGeometry(const Matrix& geometry)
{
    updateMoleculeGeometry(geometry);

    if (!m_orca) return false;

    bool ok = m_orca->updateGeometryInInput(m_molecule);
    if (!ok) {
        m_has_error = m_orca->hasError();
        m_error_message = m_orca->getErrorMessage();
    }
    return ok;
}

double OrcaMethod::calculateEnergy(bool gradient)
{
    if (!m_orca) {
        m_has_error = true;
        m_error_message = "OrcaInterface not initialized";
        return 0.0;
    }

    CitationRegistry::cite("orca");

    m_last_energy = m_orca->calculate(gradient, CurcumaLogger::get_verbosity());
    m_calculation_done = true;

    if (m_orca->hasError()) {
        m_has_error = true;
        m_error_message = m_orca->getErrorMessage();
    }

    return m_last_energy;
}

Matrix OrcaMethod::getGradient() const
{
    if (!m_orca) return Matrix::Zero(1, 3);
    return m_orca->getGradient();
}

Vector OrcaMethod::getCharges() const
{
    if (!m_orca) return Vector::Zero(1);
    return m_orca->getCharges();
}

Vector OrcaMethod::getBondOrders() const
{
    // ORCA property.json does not expose bond orders in a simple format.
    return Vector{};
}

Position OrcaMethod::getDipole() const
{
    if (!m_orca) return Position{0.0, 0.0, 0.0};
    return m_orca->getDipole();
}

void OrcaMethod::setThreadCount(int threads)
{
    m_thread_count = threads;
    // ORCA parallelism is via %pal nprocs in the input file;
    // update the parameter so the next input generation uses it.
    m_parameters["orca_nprocs"] = threads;
    if (m_orca) {
        // Re-create interface with updated parameters
        ConfigManager orca_config("orca", m_parameters);
        m_orca = std::make_unique<OrcaInterface>(orca_config);
    }
}

void OrcaMethod::setParameters(const json& params)
{
    m_parameters = params;
    if (m_orca) {
        ConfigManager orca_config("orca", m_parameters);
        m_orca = std::make_unique<OrcaInterface>(orca_config);
    }
}

json OrcaMethod::getParameters() const
{
    return m_parameters;
}

bool OrcaMethod::hasError() const
{
    return m_has_error;
}

void OrcaMethod::clearError()
{
    m_has_error = false;
    m_error_message.clear();
    if (m_orca) m_orca->clearError();
}

std::string OrcaMethod::getErrorMessage() const
{
    return m_error_message;
}

json OrcaMethod::getEnergyDecomposition() const
{
    // ORCA composite methods do not expose energy decomposition in a
    // machine-readable format via the text/JSON output we consume.
    json energy_json = {
        {"Bond", 0.0},
        {"Angle", 0.0},
        {"Torsion", 0.0},
        {"Inversion", 0.0},
        {"Dispersion", 0.0},
        {"Coulomb", 0.0},
        {"HBond", 0.0},
        {"XBond", 0.0},
        {"ATM", 0.0},
        {"BATM", 0.0}
    };
    return energy_json;
}

void OrcaMethod::setCustomKeywords(const std::string& keywords)
{
    m_parameters["orca_keywords"] = keywords;
    if (m_orca) {
        ConfigManager orca_config("orca", m_parameters);
        m_orca = std::make_unique<OrcaInterface>(orca_config);
    }
}

std::string OrcaMethod::getInputPath() const
{
    if (!m_orca) return "";
    return m_orca->getInputPath();
}

std::string OrcaMethod::getOutputPath() const
{
    if (!m_orca) return "";
    return m_orca->getOutputPath();
}

bool OrcaMethod::isAvailable()
{
    return OrcaInterface::checkOrcaExecutable();
}

std::vector<std::string> OrcaMethod::getSupportedMethods()
{
    return {"hf-3c", "b97-3c", "r2scan-3c", "pbeh-3c", "orca"};
}
