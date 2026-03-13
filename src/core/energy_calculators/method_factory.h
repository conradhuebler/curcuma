/*
 * < Method Factory for Computational Methods >
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "computational_method.h"
#include "json.hpp"

#include <memory>
#include <string>
#include <vector>

using json = nlohmann::json;

/**
 * @brief Factory for creating computational methods with priority-based resolution
 *
 * Maps method names to ComputationalMethod instances. Priority methods like "gfn2"
 * try multiple providers (TBLite > Ulysses > XTB) based on compilation flags.
 *
 * Claude Generated: Simplified from registry-based to direct if/else dispatch
 */
class MethodFactory {
public:
    /**
     * @brief Create computational method based on method name and configuration
     * @param method_name Method identifier (e.g., "gfn2", "eht", "uff", "gfnff")
     * @param config JSON configuration for the method
     * @return Unique pointer to computational method
     * @throws MethodCreationException if method cannot be created
     */
    static std::unique_ptr<ComputationalMethod> create(
        const std::string& method_name,
        const json& config = json{}
    );

    /**
     * @brief Get list of available methods based on compilation flags
     */
    static std::vector<std::string> getAvailableMethods();

    /**
     * @brief Check if specific method is available
     */
    static bool isMethodAvailable(const std::string& method_name);

    /**
     * @brief Get information about method resolution
     */
    static json getMethodInfo(const std::string& method_name);

    /**
     * @brief Print available methods and their providers
     */
    static void printAvailableMethods();

private:
    // Priority-based method creation (multiple providers with fallback)
    static std::unique_ptr<ComputationalMethod> createGFN2(const json& config);
    static std::unique_ptr<ComputationalMethod> createGFN1(const json& config);
    static std::unique_ptr<ComputationalMethod> createIPEA1(const json& config);

    // =================================================================================
    // Explicit Method Creation (single provider)
    // =================================================================================

    static std::unique_ptr<ComputationalMethod> createEHT(const json& config);
    static std::unique_ptr<ComputationalMethod> createGFNFF(const json& config);
    static std::unique_ptr<ComputationalMethod> createForceField(const json& config);
    static std::unique_ptr<ComputationalMethod> createUlyssesExplicit(const std::string& method, const json& config);

    // Explicit XTB method creation (xtb-gfn1, xtb-gfn2, gfnff)
    static std::unique_ptr<ComputationalMethod> createXTBExplicit(const std::string& method, const json& config);
    static std::unique_ptr<ComputationalMethod> createUlyssesExplicit(const std::string& method, const json& config);
    static std::unique_ptr<ComputationalMethod> createDFTD3(const json& config);
    static std::unique_ptr<ComputationalMethod> createDFTD4(const json& config);

    // Compilation flag checks
    static bool hasTBLite();
    static bool hasXTB();
    static bool hasUlysses();
    static bool hasGFNFF();
    static bool hasD3();
    static bool hasD4();
    static bool checkCompilationFlag(const std::string& flag);

    // Method classification helpers
    static bool isUlyssesMethod(const std::string& method);

    // Method lists
    static const std::vector<std::string> m_ff_methods;
    static const std::vector<std::string> m_tblite_methods;
    static const std::vector<std::string> m_xtb_methods;
    static const std::vector<std::string> m_ulysses_methods;
};

/**
 * @brief Exception thrown when method creation fails
 */
class MethodCreationException : public std::runtime_error {
public:
    explicit MethodCreationException(const std::string& message)
        : std::runtime_error(message) {}
};
