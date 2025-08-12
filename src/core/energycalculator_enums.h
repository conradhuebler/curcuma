/*
 * EnergyCalculator Enums and Type Definitions
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Refactoring to replace magic numbers with descriptive enums
 *
 * This header defines the method types and enums used by EnergyCalculator
 * to replace hardcoded magic numbers with descriptive names.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

/**
 * @brief Enum for computational methods supported by EnergyCalculator
 *
 * Replaces the magic numbers (0-9) used in switch statements with
 * descriptive names for better code readability and maintainability.
 */
enum class MethodType : int {
    FORCEFIELD = 0, ///< Force field methods (UFF, QMDFF, cgfnff)
    TBLITE = 1, ///< TBLite tight-binding DFT methods
    XTB = 2, ///< Extended tight-binding methods
    ULYSSES = 3, ///< Semi-empirical methods via Ulysses
    DFTD3 = 4, ///< DFT-D3 dispersion correction
    DFTD4 = 5, ///< DFT-D4 dispersion correction
    EHT = 6, ///< Extended Hückel Theory
    QMDFF_LEGACY = 7, ///< Legacy QMDFF (deprecated)
    UFF_LEGACY = 8, ///< Legacy UFF (deprecated)
    // NOTE: Case 9 (cgfnff) should be moved to FORCEFIELD (case 0)
    CGFNFF_LEGACY = 9, ///< TODO: Move cgfnff to ForceField system
    UNKNOWN = -1 ///< Unknown/unsupported method
};

/**
 * @brief Method detection configuration
 *
 * Maps method names to their corresponding MethodType for automatic
 * method detection from string identifiers.
 */
struct MethodConfig {
    MethodType type;
    std::vector<std::string> method_names;
    std::string description;
    bool requires_compilation_flag;
    std::string compilation_flag;

    MethodConfig(MethodType t, std::vector<std::string> names,
        const std::string& desc, bool requires_flag = false,
        const std::string& flag = "")
        : type(t)
        , method_names(names)
        , description(desc)
        , requires_compilation_flag(requires_flag)
        , compilation_flag(flag)
    {
    }
};

/**
 * @brief Method registry for centralized method configuration
 *
 * Contains all supported methods with their configuration details
 * for automatic detection and validation.
 */
class MethodRegistry {
public:
    static const std::vector<MethodConfig>& getAllMethods()
    {
        static std::vector<MethodConfig> methods = {
            // TODO: Move cgfnff from case 9 to FORCEFIELD (case 0)
            { MethodType::CGFNFF_LEGACY, { "cgfnff" }, "Native Curcuma GFN-FF (should be in ForceField)" },

            { MethodType::TBLITE, { "ipea1", "gfn1", "gfn2" }, "TBLite tight-binding DFT", true, "USE_TBLITE" },

            { MethodType::XTB, { "gfnff", "xtb-gfn1", "xtb-gfn2" }, "Extended tight-binding methods", true, "USE_XTB" },

            // TODO: Verify these Ulysses method names - taken from current code but need validation
            { MethodType::ULYSSES, { "ugfn2", "GFN2L", "pm3", "PM3PDDG", "MNDOPDDG", "PM3BP", "RM1", "AM1", "MNDO", "MNDOd", "pm6" },
                "Semi-empirical methods via Ulysses", true, "USE_ULYSSES" },

            { MethodType::DFTD3, { "d3" }, "DFT-D3 dispersion correction", true, "USE_D3" },
            { MethodType::DFTD4, { "d4" }, "DFT-D4 dispersion correction", true, "USE_D4" },
            { MethodType::EHT, { "eht" }, "Extended Hückel Theory" },
            { MethodType::QMDFF_LEGACY, { "fqmdff" }, "Legacy QMDFF (deprecated)" },

            // TODO: Add cgfnff here when moved from case 9
            { MethodType::FORCEFIELD, { "uff", "uff-d3", "qmdff" /* TODO: add "cgfnff" */ }, "Force field methods" },

            { MethodType::UFF_LEGACY, { "fuff" }, "Legacy UFF (deprecated)" }
        };
        return methods;
    }

    /**
     * @brief Detect method type from string identifier
     * @param method_name Method name to detect
     * @return MethodType corresponding to the method name
     */
    static MethodType detectMethodType(const std::string& method_name)
    {
        for (const auto& config : getAllMethods()) {
            for (const auto& name : config.method_names) {
                if (name.find(method_name) != std::string::npos || method_name.find(name) != std::string::npos) {
                    return config.type;
                }
            }
        }
        return MethodType::UNKNOWN;
    }

    /**
     * @brief Get method configuration by type
     * @param type MethodType to look up
     * @return MethodConfig for the specified type
     */
    static const MethodConfig* getMethodConfig(MethodType type)
    {
        for (const auto& config : getAllMethods()) {
            if (config.type == type) {
                return &config;
            }
        }
        return nullptr;
    }

    // TODO: Verify and update method names, especially:
    // 1. Check Ulysses documentation for correct method identifiers
    // 2. Validate TBLite method names against actual implementation
    // 3. Confirm XTB method naming convention
    // 4. Move cgfnff from case 9 to ForceField system (case 0)
};