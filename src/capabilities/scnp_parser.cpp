/*
 * <SCNP Input Parser - Molsim format to Curcuma JSON converter>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "scnp_parser.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <cmath>
#include <iomanip>

ScnpInputParser::ScnpInputParser(bool silent)
    : m_silent(silent)
{
    initializeParameterMappings();

    if (!m_silent) {
        CurcumaLogger::info("SCNP/Molsim input parser initialized");
        CurcumaLogger::info("Supporting Fortran namelist format with automatic parameter documentation");
    }
}

bool ScnpInputParser::isScnpFormat(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    bool found_namelist = false;
    int lines_checked = 0;

    // Check first 10 lines for namelist indicators
    while (getline(file, line) && lines_checked < 10) {
        // Remove whitespace and convert to lowercase
        std::string trimmed = line;
        trimmed.erase(std::remove_if(trimmed.begin(), trimmed.end(), ::isspace), trimmed.end());
        std::transform(trimmed.begin(), trimmed.end(), trimmed.begin(), ::tolower);

        // Look for Fortran namelist indicators
        if (trimmed.find("&nml") == 0 ||
            trimmed.find("txtitle") != std::string::npos ||
            trimmed.find("txmethod") != std::string::npos ||
            trimmed.find("txensemb") != std::string::npos ||
            trimmed.find("nstep") != std::string::npos) {
            found_namelist = true;
            break;
        }
        lines_checked++;
    }

    file.close();
    return found_namelist;
}

json ScnpInputParser::parseScnpFile(const std::string& filename)
{
    if (!m_silent) {
        CurcumaLogger::info_fmt("Parsing SCNP/molsim input file: {}", filename);
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        CurcumaLogger::error_fmt("Could not open SCNP input file: {}", filename);
        return json{};
    }

    // Read entire file content
    std::string content;
    std::string line;
    while (getline(file, line)) {
        content += line + "\n";
    }
    file.close();

    // Parse the content
    json raw_params = parseNamelist(content);

    if (!m_silent) {
        CurcumaLogger::info_fmt("Parsed {} SCNP parameters", raw_params.size());
    }

    // Document all found parameters
    for (const auto& [key, value] : raw_params.items()) {
        if (m_parameter_db.find(key) == m_parameter_db.end()) {
            documentUnknownParameter(key, value);
        }
    }

    // Convert to Curcuma format
    json curcuma_config = convertToStandardConfig(raw_params);

    // Generate warnings
    generateWarnings(raw_params);

    // Document parsing results
    if (!m_silent) {
        documentParameters();
    }

    return curcuma_config;
}

json ScnpInputParser::parseNamelist(const std::string& content)
{
    json params;
    std::string line;
    std::istringstream stream(content);
    bool in_namelist = false;

    while (getline(stream, line)) {
        // Remove comments (starting with ! or #)
        size_t comment_pos = line.find_first_of("!#");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line.empty()) continue;

        // Check for namelist start
        if (line.find("&") == 0) {
            in_namelist = true;
            continue;
        }

        // Check for namelist end
        if (line.find("/") == 0 || line.find("&end") == 0) {
            in_namelist = false;
            continue;
        }

        // Parse parameter lines
        if (in_namelist || line.find("=") != std::string::npos) {
            auto [param_name, param_value] = parseParameterLine(line);
            if (!param_name.empty()) {
                params[param_name] = param_value;
            }
        }
    }

    return params;
}

std::pair<std::string, json> ScnpInputParser::parseParameterLine(const std::string& line)
{
    size_t equals_pos = line.find("=");
    if (equals_pos == std::string::npos) {
        return {"", json{}};
    }

    // Extract parameter name
    std::string name = line.substr(0, equals_pos);
    name.erase(0, name.find_first_not_of(" \t"));
    name.erase(name.find_last_not_of(" \t") + 1);

    // Extract value
    std::string value_str = line.substr(equals_pos + 1);
    value_str.erase(0, value_str.find_first_not_of(" \t"));
    value_str.erase(value_str.find_last_not_of(" \t,") + 1);

    // Remove trailing comma if present
    if (!value_str.empty() && value_str.back() == ',') {
        value_str.pop_back();
    }

    json value = inferType(value_str);
    return {name, value};
}

json ScnpInputParser::inferType(const std::string& value)
{
    if (value.empty()) return json{};

    // Remove quotes for string detection
    std::string trimmed = value;
    if ((trimmed.front() == '\'' && trimmed.back() == '\'') ||
        (trimmed.front() == '"' && trimmed.back() == '"')) {
        return trimmed.substr(1, trimmed.length() - 2); // String
    }

    // Boolean detection
    std::string lower_value = value;
    std::transform(lower_value.begin(), lower_value.end(), lower_value.begin(), ::tolower);
    if (lower_value == "true" || lower_value == ".true." || lower_value == "t") return true;
    if (lower_value == "false" || lower_value == ".false." || lower_value == "f") return false;

    // Number detection (including expressions)
    try {
        if (value.find("*") != std::string::npos ||
            value.find("/") != std::string::npos ||
            value.find("+") != std::string::npos ||
            value.find("-") != std::string::npos) {
            // Mathematical expression
            return evaluateExpression(value);
        } else if (value.find(".") != std::string::npos) {
            // Floating point
            return std::stod(value);
        } else {
            // Integer
            return std::stoi(value);
        }
    } catch (const std::exception&) {
        // Fallback to string if number parsing fails
        return value;
    }
}

double ScnpInputParser::evaluateExpression(const std::string& expr)
{
    // Simple expression evaluator for basic math operations
    // Handle cases like "3*1200.0"

    std::string clean_expr = expr;
    clean_expr.erase(std::remove_if(clean_expr.begin(), clean_expr.end(), ::isspace), clean_expr.end());

    // Find operator
    size_t op_pos = clean_expr.find_first_of("*/+-");
    if (op_pos == std::string::npos) {
        return std::stod(clean_expr);
    }

    char op = clean_expr[op_pos];
    double left = std::stod(clean_expr.substr(0, op_pos));
    double right = std::stod(clean_expr.substr(op_pos + 1));

    switch (op) {
        case '*': return left * right;
        case '/': return left / right;
        case '+': return left + right;
        case '-': return left - right;
        default: return left;
    }
}

void ScnpInputParser::initializeParameterMappings()
{
    // Core simulation parameters - Claude Generated
    m_known_mappings["txtitle"] = "simulation_title";
    m_known_mappings["txmethod"] = "method_type";
    m_known_mappings["txmode"] = "simulation_mode";
    m_known_mappings["txensemb"] = "ensemble";
    m_known_mappings["txbc"] = "boundary_conditions";
    m_known_mappings["txstart"] = "starting_configuration";

    // Physical parameters
    m_known_mappings["temp"] = "temperature";
    m_known_mappings["nstep1"] = "steps";
    m_known_mappings["nstep"] = "steps";
    m_known_mappings["boxlen"] = "box_length";
    m_known_mappings["press"] = "pressure";

    // Output and control
    m_known_mappings["nfreq"] = "output_frequency";
    m_known_mappings["nsave"] = "save_frequency";
    m_known_mappings["nprint"] = "print_frequency";
    m_known_mappings["seed"] = "seed";

    // Initialize parameter database with known parameters
    ParameterInfo info;

    // txtitle
    info = {"txtitle", "Simulation title/description", "simulation_title", "Casino MC Simulation", true, "scnp", {}};
    m_parameter_db["txtitle"] = info;

    // txmethod
    info = {"txmethod", "Simulation method (MC, MD, etc.)", "method_type", "mc", true, "scnp",
            {"Curcuma only supports Monte Carlo (mc) method"}};
    m_parameter_db["txmethod"] = info;

    // txensemb
    info = {"txensemb", "Statistical ensemble (NVT, NPT, etc.)", "ensemble", "nvt", true, "scnp",
            {"Curcuma primarily supports NVT ensemble"}};
    m_parameter_db["txensemb"] = info;

    // temp
    info = {"temp", "Temperature in Kelvin", "temperature", 298.15, true, "scnp", {}};
    m_parameter_db["temp"] = info;

    // nstep1/nstep
    info = {"nstep1", "Number of simulation steps", "steps", 10000, true, "scnp", {}};
    m_parameter_db["nstep1"] = info;
    info.name = "nstep";
    m_parameter_db["nstep"] = info;

    // boxlen
    info = {"boxlen", "Box length (cubic or array)", "box_length", 1200.0, true, "scnp",
            {"Curcuma uses cell parameters in Molecule class"}};
    m_parameter_db["boxlen"] = info;
}

void ScnpInputParser::documentUnknownParameter(const std::string& name, const json& value)
{
    if (m_parameter_db.find(name) != m_parameter_db.end()) return;

    ParameterInfo info;
    info.name = name;
    info.is_known = false;
    info.source_format = "scnp";
    info.default_value = value;

    // Generate educated guesses based on parameter names
    std::string lower_name = name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);

    if (lower_name.find("temp") != std::string::npos) {
        info.probable_meaning = "Temperature-related parameter";
        info.curcuma_equivalent = "temperature";
    } else if (lower_name.find("step") != std::string::npos) {
        info.probable_meaning = "Step count or step size parameter";
        info.curcuma_equivalent = "steps";
    } else if (lower_name.find("freq") != std::string::npos) {
        info.probable_meaning = "Output or update frequency";
        info.curcuma_equivalent = "output_frequency";
    } else if (lower_name.find("box") != std::string::npos || lower_name.find("cell") != std::string::npos) {
        info.probable_meaning = "Box/cell dimension parameter";
        info.curcuma_equivalent = "box_length";
    } else if (lower_name.find("press") != std::string::npos) {
        info.probable_meaning = "Pressure-related parameter";
        info.curcuma_equivalent = "pressure";
        info.notes.push_back("Pressure control not fully implemented in Curcuma");
    } else if (lower_name.find("cut") != std::string::npos) {
        info.probable_meaning = "Cutoff distance parameter";
        info.curcuma_equivalent = "cutoff_distance";
    } else if (lower_name.find("seed") != std::string::npos || lower_name.find("rand") != std::string::npos) {
        info.probable_meaning = "Random number seed";
        info.curcuma_equivalent = "seed";
    } else {
        info.probable_meaning = "Unknown parameter - please consult molsim documentation";
        info.curcuma_equivalent = "unknown";
        info.notes.push_back("This parameter was not recognized - using as-is with default value");
    }

    m_parameter_db[name] = info;
    m_unknown_parameters.push_back(name);
}

std::string ScnpInputParser::mapToStandardParameter(const std::string& scnp_param)
{
    auto it = m_known_mappings.find(scnp_param);
    if (it != m_known_mappings.end()) {
        return it->second;
    }
    return scnp_param; // Use original name if no mapping found
}

json ScnpInputParser::convertParameter(const std::string& scnp_param, const json& value)
{
    // Special conversion cases
    if (scnp_param == "txmethod") {
        if (value.is_string() && value.get<std::string>() != "mc") {
            m_warnings.push_back("Only Monte Carlo (mc) method supported - forcing MC simulation");
            return "mc";
        }
        return "mc"; // Force MC method
    }

    if (scnp_param == "txensemb") {
        if (value.is_string() && value.get<std::string>() != "nvt") {
            m_warnings.push_back("Only NVT ensemble fully supported - some features may not work");
        }
        return value; // Keep original but warn
    }

    if (scnp_param == "boxlen") {
        // Convert box length to appropriate format
        if (value.is_number()) {
            return value; // Single cubic box dimension
        }
        // For arrays, would need more complex parsing
        return value;
    }

    return value; // Default: no conversion
}

json ScnpInputParser::convertToStandardConfig(const json& scnp_params)
{
    json curcuma_config;

    // Start with default Casino configuration
    curcuma_config = {
        {"steps", 10000},
        {"temperature", 300.0},
        {"step_size", 0.1},
        {"method", "uff"},
        {"move_strategy", "single_atom"},
        {"output_frequency", 100},
        {"verbose", true}
    };

    // Convert known parameters
    for (const auto& [scnp_name, scnp_value] : scnp_params.items()) {
        std::string curcuma_name = mapToStandardParameter(scnp_name);
        json converted_value = convertParameter(scnp_name, scnp_value);

        // Map to standard Curcuma parameters
        if (curcuma_name == "temperature") {
            curcuma_config["temperature"] = converted_value;
        } else if (curcuma_name == "steps") {
            curcuma_config["steps"] = converted_value;
        } else if (curcuma_name == "output_frequency") {
            curcuma_config["output_frequency"] = converted_value;
        } else if (curcuma_name == "seed") {
            curcuma_config["seed"] = converted_value;
        } else if (curcuma_name == "box_length") {
            // Store box information for later use in molecule setup
            curcuma_config["box_length"] = converted_value;
        } else {
            // Store unknown parameters in a special section for documentation
            if (!curcuma_config.contains("scnp_original")) {
                curcuma_config["scnp_original"] = json{};
            }
            curcuma_config["scnp_original"][scnp_name] = converted_value;
        }
    }

    return curcuma_config;
}

void ScnpInputParser::generateWarnings(const json& scnp_params)
{
    // Check for unsupported features
    for (const auto& [key, value] : scnp_params.items()) {
        std::string lower_key = key;
        std::transform(lower_key.begin(), lower_key.end(), lower_key.begin(), ::tolower);

        if (lower_key.find("md") != std::string::npos && lower_key.find("method") != std::string::npos) {
            m_warnings.push_back("Molecular dynamics not supported - converting to Monte Carlo");
        }

        if (lower_key.find("npt") != std::string::npos || lower_key.find("press") != std::string::npos) {
            m_warnings.push_back("NPT ensemble and pressure control not fully implemented");
        }

        if (lower_key.find("ewald") != std::string::npos || lower_key.find("pme") != std::string::npos) {
            m_warnings.push_back("Long-range electrostatics (Ewald/PME) not implemented");
        }
    }
}

void ScnpInputParser::documentParameters() const
{
    if (m_parameter_db.empty()) return;

    CurcumaLogger::info("=== SCNP Parameter Documentation ===");

    // Known parameters
    CurcumaLogger::info("Known Parameters:");
    for (const auto& [name, info] : m_parameter_db) {
        if (info.is_known) {
            CurcumaLogger::param(name, fmt::format("{} → {}", info.probable_meaning, info.curcuma_equivalent));
        }
    }

    // Unknown parameters
    if (!m_unknown_parameters.empty()) {
        CurcumaLogger::warn("Unknown Parameters (educated guesses):");
        for (const std::string& name : m_unknown_parameters) {
            const auto& info = m_parameter_db.at(name);
            CurcumaLogger::param(name, fmt::format("Probably: {} → {}", info.probable_meaning, info.curcuma_equivalent));
        }
    }

    // Warnings
    if (!m_warnings.empty()) {
        CurcumaLogger::warn("Conversion Warnings:");
        for (const std::string& warning : m_warnings) {
            CurcumaLogger::warn(warning);
        }
    }
}

ParameterInfo ScnpInputParser::getParameterInfo(const std::string& param_name) const
{
    auto it = m_parameter_db.find(param_name);
    if (it != m_parameter_db.end()) {
        return it->second;
    }

    // Return empty info for unknown parameters
    ParameterInfo empty_info;
    empty_info.name = param_name;
    empty_info.is_known = false;
    empty_info.probable_meaning = "Unknown parameter";
    return empty_info;
}

void ScnpInputParser::printParameterMappingGuide() const
{
    std::cout << "SCNP/Molsim to Curcuma Parameter Mapping Guide" << std::endl;
    std::cout << "=" << std::string(50, '=') << std::endl;
    std::cout << std::endl;
    std::cout << std::left << std::setw(15) << "SCNP Parameter"
              << std::setw(25) << "Curcuma Parameter"
              << "Description" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (const auto& [name, info] : m_parameter_db) {
        if (info.is_known) {
            std::cout << std::left << std::setw(15) << name
                      << std::setw(25) << info.curcuma_equivalent
                      << info.probable_meaning << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "Usage Example:" << std::endl;
    std::cout << "  curcuma -casino structure.xyz -input scnp.in" << std::endl;
    std::cout << "  # Automatically converts SCNP parameters to Curcuma format" << std::endl;
    std::cout << std::endl;
}

// Factory function implementation
json parseInputFile(const std::string& filename, bool silent)
{
    if (ScnpInputParser::isScnpFormat(filename)) {
        if (!silent) {
            CurcumaLogger::info("Detected SCNP/molsim format - parsing with SCNP parser");
        }
        ScnpInputParser parser(silent);
        return parser.parseScnpFile(filename);
    } else {
        if (!silent) {
            CurcumaLogger::info("Standard JSON format detected");
        }
        // Fall back to standard JSON parsing
        std::ifstream file(filename);
        if (file.is_open()) {
            json config;
            file >> config;
            return config;
        }
    }

    return json{};
}