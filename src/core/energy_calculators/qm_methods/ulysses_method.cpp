/*
 * < Ulysses Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "ulysses_method.h"
#include "src/core/curcuma_logger.h"
#include "src/core/config_manager.h"
#include <algorithm>
#include <cctype>

// Claude Generated: Helper function to parse Ulysses method names with correction suffixes
std::pair<std::string, std::string> parseMethodName(const std::string& method_name)
{
    std::string method_lower = method_name;
    std::transform(method_lower.begin(), method_lower.end(), method_lower.begin(), ::tolower);

    // C++17 compatible string suffix checking
    if (method_lower.length() >= 6 && method_lower.substr(method_lower.length() - 6) == "-d3h4x") {
        std::string base = method_lower.substr(0, method_lower.length() - 6);
        return { base, "D3H4X" };
    } else if (method_lower.length() >= 5 && method_lower.substr(method_lower.length() - 5) == "-d3h+") {
        std::string base = method_lower.substr(0, method_lower.length() - 5);
        return { base, "D3H+" };
    }

    // No correction suffix found
    return { method_lower, "0" };
}

UlyssesMethod::UlyssesMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name), m_calculation_done(false), m_last_energy(0.0) {
    CurcumaLogger::info("UlyssesMethod constructor called");
    CurcumaLogger::param("method", method_name);

    // Parse method name for correction suffixes
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Parsing Ulysses method name");
        CurcumaLogger::param("input_method_name", method_name);
    }

    auto [base_method, correction] = parseMethodName(method_name);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Ulysses method parsed");
        CurcumaLogger::param("base_method", base_method);
        CurcumaLogger::param("correction", correction);
    }

#ifdef USE_ULYSSES
    try {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Creating UlyssesInterface");
            CurcumaLogger::param("config", config.dump(2));
        }

        // Create enhanced config with parsed correction (wrap JSON in ConfigManager)
        json enhanced_config = config;
        enhanced_config["base_method"] = base_method;
        enhanced_config["corecorrection"] = correction;

        ConfigManager config_mgr("ulysses", enhanced_config);
        m_ulysses = std::make_unique<UlyssesInterface>(config_mgr);
        CurcumaLogger::success("UlyssesInterface created successfully");
    } catch (const std::exception& e) {
        CurcumaLogger::error("UlyssesInterface creation failed: " + std::string(e.what()));
        throw;
    }
#else
    CurcumaLogger::error("USE_ULYSSES not defined - Ulysses not available");
#endif
    m_parameters = config;
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
    // Claude Generated: Complete list of Ulysses methods with correction mode support
    std::vector<std::string> base_methods = {
        "ugfn2", "pm6", "am1", "pm3", "mndo", "mndod",
        "rm1", "pm3pddg", "mndopddg", "pm3bp"
    };

    std::vector<std::string> all_methods;

    // Add base methods (without corrections)
    for (const auto& base : base_methods) {
        all_methods.push_back(base);
    }

    // Add methods with D3H4X corrections (skip ugfn2 as it's special)
    for (const auto& base : base_methods) {
        if (base != "ugfn2") {
            all_methods.push_back(base + "-d3h4x");
        }
    }

    // Add methods with D3H+ corrections (skip ugfn2 as it's special)
    for (const auto& base : base_methods) {
        if (base != "ugfn2") {
            all_methods.push_back(base + "-d3h+");
        }
    }

    return all_methods;
}