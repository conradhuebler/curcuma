/*
 * <Centralized parameter validation utilities>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/curcuma_logger.h"
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>

namespace ParamValidation {

/*
 * Phase 3: Centralized Parameter Validation - Claude Generated 2025
 *
 * Reusable validation functions with educational error messages.
 * All functions show warnings but allow auto-correction to maintain
 * backward compatibility while alerting users to potential issues.
 */

/**
 * @brief Validate that a numeric value is positive
 * @param param_name Parameter name for error messages
 * @param value The value to validate
 * @param default_value Default to use if validation fails
 * @return Either the value or default_value if validation failed
 */
template<typename T>
inline T validate_positive(const std::string& param_name, T value, T default_value)
{
    if (value <= 0) {
        CurcumaLogger::warn_fmt(
            "Parameter '{}' must be positive, got {}. Using default: {}",
            param_name, value, default_value);
        return default_value;
    }
    return value;
}

/**
 * @brief Validate that a numeric value is within range
 * @param param_name Parameter name for error messages
 * @param value The value to validate
 * @param min Minimum allowed value (inclusive)
 * @param max Maximum allowed value (inclusive)
 * @param default_value Default to use if validation fails
 * @return Either the value or default_value if validation failed
 */
template<typename T>
inline T validate_range(
    const std::string& param_name,
    T value,
    T min,
    T max,
    T default_value)
{
    if (value < min || value > max) {
        CurcumaLogger::warn_fmt(
            "Parameter '{}' must be in range [{}, {}], got {}. Using default: {}",
            param_name, min, max, value, default_value);
        return default_value;
    }
    return value;
}

/**
 * @brief Validate that a string choice is valid
 * @param param_name Parameter name for error messages
 * @param value The value to validate
 * @param valid_choices List of allowed values
 * @param default_value Default to use if validation fails
 * @return Either the value or default_value if validation failed
 */
inline std::string validate_choice(
    const std::string& param_name,
    const std::string& value,
    const std::vector<std::string>& valid_choices,
    const std::string& default_value)
{
    auto it = std::find(valid_choices.begin(), valid_choices.end(), value);
    if (it == valid_choices.end()) {
        std::string choices_str;
        for (size_t i = 0; i < valid_choices.size(); ++i) {
            if (i > 0) choices_str += ", ";
            choices_str += "'" + valid_choices[i] + "'";
        }
        CurcumaLogger::warn_fmt(
            "Invalid value for '{}': '{}'. Valid choices: [{}]. Using default: '{}'",
            param_name, value, choices_str, default_value);
        return default_value;
    }
    return value;
}

/**
 * @brief Validate that a file exists (for file parameters)
 * @param param_name Parameter name for error messages
 * @param filepath Path to check
 * @return true if file exists, false otherwise
 */
inline bool validate_file_exists(const std::string& param_name, const std::string& filepath)
{
    if (!std::filesystem::exists(filepath)) {
        CurcumaLogger::error_fmt(
            "File specified for '{}' not found: {}", param_name, filepath);
        return false;
    }
    return true;
}

/**
 * @brief Validate interdependent parameters
 * @param param1_name First parameter name
 * @param param1_enabled Whether first parameter is enabled
 * @param param2_name Second parameter name
 * @param param2_enabled Whether second parameter is enabled
 * @return true if validation passed, false otherwise
 *
 * Example: scattering_per_frame requires scattering_enabled=true
 */
inline bool validate_dependency(
    const std::string& param1_name,
    bool param1_enabled,
    const std::string& param2_name,
    bool param2_enabled)
{
    if (param1_enabled && !param2_enabled) {
        CurcumaLogger::warn_fmt(
            "Parameter '{}' requires '{}=true'. Disabling '{}'.",
            param1_name, param2_name, param1_name);
        return false;
    }
    return true;
}

} // namespace ParamValidation
