/*
 * < Method Factory for Computational Methods >
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

#include "computational_method.h"
#include "json.hpp"

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <functional>

using json = nlohmann::json;

/**
 * @brief Factory for creating computational methods with priority-based resolution
 * 
 * This factory implements the method resolution logic from the original 
 * EnergyCalculator::SwitchMethod(), maintaining the hierarchical priority
 * system where methods like "gfn2" can be provided by multiple libraries
 * (TBLite > Ulysses > XTB) based on compilation flags and availability.
 * 
 * Claude Generated: Big-Bang refactoring of EnergyCalculator method dispatch
 */
class MethodFactory {
public:
    /**
     * @brief Create computational method based on method name and configuration
     * @param method_name Method identifier (e.g., "gfn2", "eht", "uff", "cgfnff")
     * @param config JSON configuration for the method
     * @return Unique pointer to computational method or nullptr on failure
     */
    static std::unique_ptr<ComputationalMethod> create(
        const std::string& method_name, 
        const json& config = json{}
    );
    
    /**
     * @brief Get list of available methods based on compilation flags
     * @return Vector of available method names
     */
    static std::vector<std::string> getAvailableMethods();
    
    /**
     * @brief Check if specific method is available
     * @param method_name Method to check
     * @return true if method can be created
     */
    static bool isMethodAvailable(const std::string& method_name);
    
    /**
     * @brief Get information about method resolution
     * @param method_name Method to query
     * @return JSON with resolution details (provider, priority, etc.)
     */
    static json getMethodInfo(const std::string& method_name);
    
    /**
     * @brief Print available methods and their providers (for debugging)
     */
    static void printAvailableMethods();
    
private:
    // =================================================================================
    // Priority-based Method Resolution (from original SwitchMethod logic)
    // =================================================================================
    
    /**
     * @brief Create GFN2 method with priority: TBLite > Ulysses > XTB > Native
     */
    static std::unique_ptr<ComputationalMethod> createGFN2(const json& config);
    
    /**
     * @brief Create GFN1 method with priority: TBLite > XTB > Native
     */
    static std::unique_ptr<ComputationalMethod> createGFN1(const json& config);
    
    /**
     * @brief Create IPEA1 method (TBLite only)
     */
    static std::unique_ptr<ComputationalMethod> createIPEA1(const json& config);
    
    // =================================================================================
    // Explicit Method Creation (single provider)
    // =================================================================================
    
    /**
     * @brief Create Extended Hückel Theory method
     */
    static std::unique_ptr<ComputationalMethod> createEHT(const json& config);
    
    /**
     * @brief Create native GFN-FF method (cgfnff)
     */
    static std::unique_ptr<ComputationalMethod> createCGFNFF(const json& config);
    
    /**
     * @brief Create UFF/QMDFF force field methods
     */
    static std::unique_ptr<ComputationalMethod> createForceField(const json& config);
    
    /**
     * @brief Create explicit Ulysses method (ugfn2, pm3, am1, etc.)
     */
    static std::unique_ptr<ComputationalMethod> createUlyssesExplicit(const std::string& method, const json& config);
    
    /**
     * @brief Create explicit XTB method (xtb-gfn1, xtb-gfn2, gfnff)
     */
    static std::unique_ptr<ComputationalMethod> createXTBExplicit(const std::string& method, const json& config);
    
    /**
     * @brief Create DFT-D3 dispersion correction
     */
    static std::unique_ptr<ComputationalMethod> createDFTD3(const json& config);
    
    /**
     * @brief Create DFT-D4 dispersion correction
     */
    static std::unique_ptr<ComputationalMethod> createDFTD4(const json& config);
    
    // =================================================================================
    // Compilation Flag Checks (from original isCompiled logic)
    // =================================================================================
    
    static bool hasTBLite();
    static bool hasXTB();
    static bool hasUlysses();
    static bool hasD3();
    static bool hasD4();
    
    // =================================================================================
    // Method Registration System
    // =================================================================================
    
    /**
     * @brief Method priority entry for shared methods
     */
    struct MethodPriority {
        std::string method_name;
        std::vector<std::pair<std::string, std::function<std::unique_ptr<ComputationalMethod>(const json&)>>> priorities;
        // Format: {{"provider_name", creation_function}, ...} in priority order
    };
    
    /**
     * @brief Explicit method entry for single-provider methods
     */
    struct ExplicitMethod {
        std::string method_name;
        std::function<std::unique_ptr<ComputationalMethod>(const json&)> creator;
        std::string provider_name;
        bool requires_compilation_flag;
        std::string compilation_flag;
    };
    
    /**
     * @brief Get priority method registry
     */
    static const std::vector<MethodPriority>& getPriorityMethods();
    
    /**
     * @brief Get explicit method registry  
     */
    static const std::vector<ExplicitMethod>& getExplicitMethods();
    
    /**
     * @brief Check if method name matches any pattern in list
     */
    static bool matchesMethodList(const std::string& method, const std::vector<std::string>& method_list);
    
    // =================================================================================
    // Method Lists (from original EnergyCalculator)
    // =================================================================================
    
    static const std::vector<std::string> m_ff_methods;      // {"uff", "uff-d3", "qmdff", "cgfnff"}
    static const std::vector<std::string> m_tblite_methods;  // {"ipea1", "gfn1", "gfn2"}
    static const std::vector<std::string> m_xtb_methods;     // {"gfnff", "xtb-gfn1", "xtb-gfn2"}
    static const std::vector<std::string> m_ulysses_methods; // {"ugfn2", "GFN2L", "pm3", "PM3PDDG", ...}
    static const std::vector<std::string> m_d3_methods;      // {"d3"}
    static const std::vector<std::string> m_d4_methods;      // {"d4"}
};

/**
 * @brief Exception thrown when method creation fails
 */
class MethodCreationException : public std::runtime_error {
public:
    explicit MethodCreationException(const std::string& message) 
        : std::runtime_error(message) {}
};

/**
 * @brief Utility macros for method availability checking
 */
#define CURCUMA_HAS_TBLITE() MethodFactory::hasTBLite()
#define CURCUMA_HAS_XTB() MethodFactory::hasXTB()  
#define CURCUMA_HAS_ULYSSES() MethodFactory::hasUlysses()
#define CURCUMA_HAS_D3() MethodFactory::hasD3()
#define CURCUMA_HAS_D4() MethodFactory::hasD4()