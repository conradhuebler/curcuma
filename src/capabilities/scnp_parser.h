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

#pragma once

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <regex>

#include "src/core/curcuma_logger.h"
#include "json.hpp"

using json = nlohmann::json;

// Parameter documentation structure - Claude Generated
struct ParameterInfo {
    std::string name;                    // Original parameter name
    std::string probable_meaning;        // What the parameter probably does
    std::string curcuma_equivalent;      // Corresponding Curcuma parameter
    json default_value;                  // Safe default value
    bool is_known;                       // Whether we understand this parameter
    std::string source_format;           // "scnp" or "molsim"
    std::vector<std::string> notes;      // Additional notes or warnings
};

/*! \brief SCNP/Molsim Input Parser - Claude Generated
 *
 * Parses molsim scnp.in format files and converts them to Curcuma JSON configuration.
 * Provides automatic documentation for unknown parameters and intelligent mapping
 * of molsim concepts to Curcuma Monte Carlo simulation parameters.
 *
 * Supported formats:
 * - Fortran namelist format (&nmlSystem ... /)
 * - Key-value pairs with automatic type detection
 * - Mathematical expressions (e.g., "3*1200.0")
 * - String values with quotes
 *
 * Features:
 * - Automatic parameter documentation
 * - Unknown parameter handling with safe defaults
 * - Type inference and conversion
 * - Warning system for unsupported features
 * - Comprehensive mapping documentation
 */
class ScnpInputParser
{
public:
    ScnpInputParser(bool silent = false);
    ~ScnpInputParser() = default;

    /*! \brief Parse SCNP input file and convert to Curcuma JSON */
    json parseScnpFile(const std::string& filename);

    /*! \brief Convert raw SCNP parameters to Curcuma configuration */
    json convertToStandardConfig(const json& scnp_params);

    /*! \brief Auto-detect if file is in SCNP/molsim format */
    static bool isScnpFormat(const std::string& filename);

    /*! \brief Document all parameters found during parsing */
    void documentParameters() const;

    /*! \brief Get parameter documentation for specific parameter */
    ParameterInfo getParameterInfo(const std::string& param_name) const;

    /*! \brief Print comprehensive parameter mapping guide */
    void printParameterMappingGuide() const;

    /*! \brief Get all documented parameters */
    const std::map<std::string, ParameterInfo>& getParameterDatabase() const { return m_parameter_db; }

    /*! \brief Clear documented parameters */
    void clearParameterDatabase() { m_parameter_db.clear(); m_unknown_parameters.clear(); }

private:
    /*! \brief Parse Fortran namelist format */
    json parseNamelist(const std::string& content);

    /*! \brief Parse single parameter line */
    std::pair<std::string, json> parseParameterLine(const std::string& line);

    /*! \brief Evaluate mathematical expressions */
    double evaluateExpression(const std::string& expr);

    /*! \brief Infer JSON type from string value */
    json inferType(const std::string& value);

    /*! \brief Initialize known parameter mappings */
    void initializeParameterMappings();

    /*! \brief Document unknown parameter with educated guess */
    void documentUnknownParameter(const std::string& name, const json& value);

    /*! \brief Map SCNP parameter to Curcuma equivalent */
    std::string mapToStandardParameter(const std::string& scnp_param);

    /*! \brief Validate and apply parameter conversion */
    json convertParameter(const std::string& scnp_param, const json& value);

    /*! \brief Generate warnings for unsupported features */
    void generateWarnings(const json& scnp_params);

    // Core data
    bool m_silent;
    std::map<std::string, ParameterInfo> m_parameter_db;
    std::vector<std::string> m_unknown_parameters;
    std::vector<std::string> m_warnings;

    // Known parameter mappings: SCNP → Curcuma
    std::map<std::string, std::string> m_known_mappings;
    std::map<std::string, json> m_default_values;
};

/*! \brief Factory function to detect and parse various input formats - Claude Generated */
json parseInputFile(const std::string& filename, bool silent = false);