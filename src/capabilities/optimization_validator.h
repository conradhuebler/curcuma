/*
 * <Optimization Parameter Validator - JSON Schema Validation>
 * Copyright (C) 2025 Claude AI - Generated Code
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

#include "optimizer_interface.h"
#include <string>
#include <vector>

namespace Optimization {

/**
 * @brief Validation result container - Claude Generated
 */
struct ValidationResult {
    bool is_valid = true;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    void addError(const std::string& error)
    {
        errors.push_back(error);
        is_valid = false;
    }

    void addWarning(const std::string& warning)
    {
        warnings.push_back(warning);
    }

    std::string getErrorSummary() const;
    std::string getWarningSummary() const;
    std::string getFullReport() const;
};

/**
 * @brief Parameter bounds and validation rules - Claude Generated
 */
struct ParameterBounds {
    double min_value = -std::numeric_limits<double>::infinity();
    double max_value = std::numeric_limits<double>::infinity();
    bool required = false;
    std::string description;
    std::string unit;

    ParameterBounds() = default;
    ParameterBounds(double min, double max, bool req = false,
        const std::string& desc = "", const std::string& u = "")
        : min_value(min)
        , max_value(max)
        , required(req)
        , description(desc)
        , unit(u)
    {
    }
};

/**
 * @brief Optimization parameter validator - Claude Generated
 * Provides comprehensive validation with bounds checking and type safety
 */
class OptimizationValidator {
public:
    /**
     * @brief Validate configuration for specific optimizer type
     */
    static ValidationResult validate(const json& config, OptimizerType type);

    /**
     * @brief Validate general optimization parameters
     */
    static ValidationResult validateGeneral(const json& config);

    /**
     * @brief Validate LBFGS-specific parameters
     */
    static ValidationResult validateLBFGS(const json& config);

    /**
     * @brief Validate convergence criteria
     */
    static ValidationResult validateConvergence(const json& config);

    /**
     * @brief Validate molecular properties
     */
    static ValidationResult validateMolecularProperties(const json& config);

    /**
     * @brief Get parameter bounds for specific parameter
     */
    static ParameterBounds getParameterBounds(const std::string& parameter_name, OptimizerType type);

    /**
     * @brief Get all parameter bounds for optimizer type
     */
    static std::map<std::string, ParameterBounds> getAllParameterBounds(OptimizerType type);

    /**
     * @brief Sanitize configuration (fix out-of-bounds values)
     */
    static json sanitizeConfiguration(const json& config, OptimizerType type,
        ValidationResult* warnings = nullptr);

    /**
     * @brief Check parameter compatibility between optimizer types
     */
    static ValidationResult checkCompatibility(const json& config,
        OptimizerType from_type,
        OptimizerType to_type);

private:
    // Parameter bounds definitions
    static const std::map<std::string, ParameterBounds> s_general_bounds;
    static const std::map<std::string, ParameterBounds> s_lbfgs_bounds;
    static const std::map<std::string, ParameterBounds> s_convergence_bounds;
    static const std::map<std::string, ParameterBounds> s_molecular_bounds;

    // Helper methods
    static void validateParameter(const json& config, const std::string& key,
        const ParameterBounds& bounds, ValidationResult& result);
    static void validateIntegerParameter(const json& config, const std::string& key,
        int min_val, int max_val, bool required,
        const std::string& description, ValidationResult& result);
    static void validateBooleanParameter(const json& config, const std::string& key,
        bool required, const std::string& description,
        ValidationResult& result);
    static void validateStringParameter(const json& config, const std::string& key,
        const std::vector<std::string>& valid_values,
        bool required, const std::string& description,
        ValidationResult& result);
};

} // namespace Optimization